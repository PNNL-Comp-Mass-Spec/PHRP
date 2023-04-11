// This class reads a DMS-based parameter file for MS-GF+ or MSPathFinder
// to extract the dynamic and static modification information
// Example param file contents:
//
// #Modifications (see below for examples)
// StaticMod=C2H3NO,   C,   fix, any,         Carbamidomethyl           # Fixed Carbamidomethyl C
//
// DynamicMod=O1,      M,   opt, any,         Oxidation                 # Oxidation M
// DynamicMod=HO3P,    STY, opt, any,         Phospho                   # Phosphorylation STY
// DynamicMod=H-1,     C,   opt, any,         Dehydro                   # Dehydro C
// DynamicMod=C2H2O,   *,   opt, Prot-N-term, Acetyl                    # Acetylation Protein N-term (C2H2O can be replaced with "H(2) C(2) O")
//
//
// Note that DMS uses this information to create a Mods.txt file that is provided to MS-GF+ or MSPathFinder
// When doing this, the static mods defs have ",fix," while the dynamic mod defs have ',opt'
// Example contents of an auto-generated Mods.txt file:
//
// # Static mods
// C2H3NO,C,fix,any,Carbamidomethyl     # Fixed Carbamidomethyl C
//
// # Dynamic mods
// O1,M,opt,any,Oxidation             # Oxidation M
// HO3P,STY,opt,any,Phospho           # Phosphorylation STY
// H-1,C,opt,any,Dehydro              # Dehydro C
// C2H2O,*,opt,Prot-N-term,Acetyl     # Acetylation Protein N-term (C2H2O can be replaced with "H(2) C(2) O")
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 7/16/2015
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.IO;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PHRPReader
{
    /// <summary>
    /// This class reads a DMS-based parameter file for MS-GF+ or MSPathFinder to extract the dynamic and static modification information
    /// </summary>
    /// <remarks>See above for an example parameter file</remarks>
    public class MSGFPlusParamFileModExtractor : PRISM.EventNotifier
    {
        // Ignore Spelling: Acetyl, Acetylation, Acetylhexosamine, Carbamidomethyl, Carbamidomethylation, cterm, defs, Dehydro
        // Ignore Spelling: Hexosam, Hydroxyproline, itrac, nterm, num, phospho, Phosphorylated, phosphorylation, Prot, protcterm, protnterm, UniMod

        /// <summary>
        /// Unknown MS-GF+ mod symbols
        /// </summary>
        public const char UNKNOWN_MSGFPlus_MOD_SYMBOL = '?';

        /// <summary>
        /// Static mod parameter file keyword
        /// </summary>
        public const string PARAM_TAG_MOD_STATIC = "StaticMod";

        /// <summary>
        /// Dynamic mod parameter file keyword
        /// </summary>
        public const string PARAM_TAG_MOD_DYNAMIC = "DynamicMod";

        /// <summary>
        /// Custom amino acid definition parameter file keyword
        /// </summary>
        public const string PARAM_TAG_CUSTOM_AA = "CustomAA";

        // ReSharper disable CommentTypo

        /// <summary>
        /// DIA-NN parameter indicating that static Cys carbamidomethyl (+57.021) was used (aka iodoacetamide alkylation)
        /// </summary>
        /// <remarks>Only used by DIA-NN</remarks>
        public const string PARAM_TAG_STATIC_CYS_CARBAMIDOMETHYL = "StaticCysCarbamidomethyl";

        // ReSharper restore CommentTypo

        private const string MSGFPLUS_COMMENT_CHAR = "#";

        /// <summary>
        /// MS-GF+ modification type
        /// </summary>
        public enum MSGFPlusModType
        {
            /// <summary>
            /// Unknown
            /// </summary>
            Unknown = 0,

            /// <summary>
            /// Dynamic
            /// </summary>
            DynamicMod = 1,

            /// <summary>
            /// Static
            /// </summary>
            StaticMod = 2,

            /// <summary>
            /// N-terminal peptide dynamic
            /// </summary>
            DynNTermPeptide = 3,

            /// <summary>
            /// C-terminal peptide dynamic
            /// </summary>
            DynCTermPeptide = 4,

            /// <summary>
            /// N-terminal protein dynamic
            /// </summary>
            DynNTermProtein = 5,

            /// <summary>
            /// C-terminal protein dynamic
            /// </summary>
            DynCTermProtein = 6,

            /// <summary>
            /// Custom amino acid definition
            /// </summary>
            CustomAA = 7
        }

        /// <summary>
        /// Modification specification formats
        /// </summary>
        public enum ModSpecFormats
        {
            /// <summary>
            /// MS-GF+, MSPathFinder, and DIA-NN
            /// </summary>
            /// <remarks>
            /// Examples:
            ///            # Mass or Formula, Residues, ModType, Position, Name
            ///   DynamicMod=O1,         M, opt, any, Oxidation         # Oxidized methionine
            ///   DynamicMod=C2H3N1O1,   C, opt, any, Carbamidomethyl   # Carbamidomethyl C
            ///   DynamicMod=199.955,    C, opt, any, Mercury           # Mercury Hg(II) adduct (+199.955)
            ///   StaticMod=229.1629,    *,  fix, N-term,    TMT6plex
            /// </remarks>
            MSGFPlusAndMSPathFinder = 0,

            /// <summary>
            /// TopPIC
            /// </summary>
            /// <remarks>
            /// Examples:
            ///           # Name,            Mass,  Residues, Position, UniModID
            ///   StaticMod=Carbamidomethyl, 57.021464,  C,   any,      4
            ///   DynamicMod=Acetyl,         42.010565,  K,   any,      1      # Acetyl
            /// </remarks>
            TopPIC = 1,

            /// <summary>
            /// DIA-NN
            /// </summary>
            /// <remarks>
            /// Examples:
            ///            # Name,        Mass,  Residues
            ///   DynamicMod=UniMod:35,   15.994915,  M       # Oxidized methionine
            ///   DynamicMod=UniMod:1,    42.010565,  *n      # Acetylation protein N-term
            ///   DynamicMod=UniMod:21,   79.966331,  STY     # Phosphorylated STY
            /// </remarks>
            DiaNN = 2,
        }

        /// <summary>
        /// Tracks dynamic and static modification details
        /// Also tracks Custom amino acids
        /// </summary>
        /// <remarks>
        /// Notes when tracking information for custom amino acids
        ///   ModName:    Name associated with the custom amino acid
        ///   ModMass:    Composition string, for example C5H7N1O2S0 for Hydroxyproline
        ///   ModMassVal: Computed mass of the composition string
        ///   Residues:   Single letter abbreviation for the custom amino acid, for example J or X
        ///   ModType:    MSGFPlusModType.CustomAA
        ///   ModSymbol:  ?   (a question mark; not used)
        /// </remarks>
        public struct ModInfo
        {
            /// <summary>
            /// True if this is an isobaric mod (e.g. TMT or iTRAQ)
            /// </summary>
            public bool IsobaricMod;

            /// <summary>
            /// Mod name (read from the parameter file)
            /// </summary>
            /// <remarks>
            /// The mod name isn't used by MS-GF+, but it is used by MSPathFinder and MaxQuant
            /// </remarks>
            public string ModName;

            /// <summary>
            /// Mod mass, stored as a string since reading from a text file and writing to a text file.  Also, can be a mass or an empirical formula
            /// </summary>
            public string ModMass;

            /// <summary>
            /// Modification mass
            /// </summary>
            public double ModMassVal;

            /// <summary>
            /// Affected residues
            /// </summary>
            /// <remarks>
            /// One letter amino acid symbols, or one of the symbols used by DMS
            /// to track N-terminal and C-Terminal peptide/protein residues
            /// (e.g., constant N_TERMINAL_PEPTIDE_SYMBOL_DMS in the <see cref="AminoAcidModInfo"/> class)
            /// </remarks>
            public string Residues;

            /// <summary>
            /// Modification type, as tracked internally by PHRP
            /// </summary>
            /// <remarks>
            /// PHRP does not support static mods at the N or C terminus that only apply to specific residues
            /// When these are found, the ModType tracked in-memory is changed to DynamicMod, DynNTermPeptide, DynCTermPeptide, DynNTermProtein, or DynCTermProtein
            /// </remarks>
            public MSGFPlusModType ModType;

            /// <summary>
            /// Modification type, as defined in the parameter file
            /// </summary>
            /// <remarks>
            /// As explained in <see cref="ModType"/>, this modification type could differ from ModType
            /// </remarks>
            public MSGFPlusModType ModTypeInParameterFile;

            /// <summary>
            /// Modification symbol: *, #, @, ... ; dash if a static mod
            /// </summary>
            public char ModSymbol;

            /// <summary>
            /// Modification name, without residues
            /// </summary>
            /// <remarks>
            /// Only used for MaxQuant results; empty string for other tools
            /// </remarks>
            public string ShortName;

            /// <summary>
            /// Mod type, name, mass, residues
            /// </summary>
            public override string ToString()
            {
                return string.Format("{0} {1}, {2}; {3}", ModType, ModName, ModMass, Residues);
            }
        }

        private string mErrorMessage;
        private readonly string mToolName;

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage => mErrorMessage;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="toolName">
        /// Search engine name, typically MS-GF+
        /// This name is only used in log messages
        /// </param>
        public MSGFPlusParamFileModExtractor(string toolName)
        {
            mErrorMessage = string.Empty;
            if (string.IsNullOrWhiteSpace(toolName))
            {
                mToolName = "Search engine";
            }
            else
            {
                mToolName = toolName;
            }
        }

        private double ComputeMass(string empiricalFormula)
        {
            // Originally only C, H, N, O, S, and P were allowed
            // We now support any element symbol
            //
            // Format is: C[Num]H[Num]N[Num]O[Num]S[Num]P[Num]
            // 	- C (Carbon), H (Hydrogen), N (Nitrogen), O (Oxygen), S (Sulfur) and P (Phosphorus) are allowed.
            // 	- Atom can be omitted.
            // 	- Negative numbers are allowed.
            // 	- Examples:
            //     C2H2O1
            //     C+2H+3N+1O+1
            //     H2C1O1
            //     H-1N-1O
            //     C3H6N2O0S1

            if (string.Equals(empiricalFormula, "HexNAc", StringComparison.OrdinalIgnoreCase))
            {
                // This is a special-case modification that MS-GF+ and MSPathFinder recognize
                // It is listed in DMS as Hexosam, which means N-Acetylhexosamine
                // It is tracked by UniMod as HexNAc
                return 203.079376;
            }

            EmpiricalFormula empiricalFormulaInstance;

            try
            {
                empiricalFormulaInstance = PeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormula);
            }
            catch (Exception ex)
            {
                ReportError(ex.Message);
                return 0;
            }

            var monoisotopicMass = PeptideMassCalculator.ComputeMonoisotopicMass(empiricalFormulaInstance.ElementCounts, out var unknownSymbols);

            if (unknownSymbols?.Count > 0)
            {
                var errMsg = "Error parsing empirical formula '" + empiricalFormula + "', ";

                if (unknownSymbols.Count == 1)
                {
                    ReportError(errMsg + "unknown element " + unknownSymbols[0]);
                }
                else
                {
                    ReportError(errMsg + "unknown elements " + string.Join(", ", unknownSymbols));
                }

                return 0;
            }

            return monoisotopicMass;
        }

        /// <summary>
        /// Extracts mod info from a MS-GF+, MSPathFinder, MSFragger, or DIA-NN param file, or from a MSGFPlus_Mods.txt file (previously MSGFDB_Mods.txt)
        /// </summary>
        /// <param name="paramFilePath"></param>
        /// <param name="modSpecFormat"></param>
        /// <param name="modList"></param>
        /// <returns>True if success; false if a problem</returns>
        public bool ExtractModInfoFromParamFile(string paramFilePath, ModSpecFormats modSpecFormat, out List<ModInfo> modList)
        {
            var tagNamesToFind = new List<string> {
                PARAM_TAG_MOD_STATIC,
                PARAM_TAG_MOD_DYNAMIC,
                PARAM_TAG_CUSTOM_AA,                 // "CustomAA" is only used by MS-GF+
                PARAM_TAG_STATIC_CYS_CARBAMIDOMETHYL // "StaticCysCarbamidomethyl" is only used by DIA-NN
            };

            // Initialization
            modList = new List<ModInfo>();

            var unnamedModID = 0;
            mErrorMessage = string.Empty;

            try
            {
                if (string.IsNullOrEmpty(paramFilePath))
                {
                    ReportError(mToolName + " Parameter File name not defined; unable to extract mod info");
                    return false;
                }

                var paramFile = new FileInfo(paramFilePath);
                if (!paramFile.Exists)
                {
                    ReportError(mToolName + " param file not found: " + paramFilePath);
                    return false;
                }

                // Read the contents of the parameter (or mods) file
                using var reader = new StreamReader(new FileStream(paramFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var dataLine = reader.ReadLine();

                    if (string.IsNullOrWhiteSpace(dataLine))
                        continue;

                    var trimmedLine = dataLine.Trim();

                    var modSpec = string.Empty;

                    if (trimmedLine.StartsWith(MSGFPLUS_COMMENT_CHAR))
                    {
                        // Comment line (starts with #)
                        // Skip it
                        continue;
                    }

                    var modType = MSGFPlusModType.Unknown;

                    foreach (var tagName in tagNamesToFind)
                    {
                        // Check whether the line starts with StaticMod, DynamicMod, or CustomAA
                        modSpec = ValidateIsValidModSpec(trimmedLine, tagName);
                        if (string.IsNullOrEmpty(modSpec))
                            continue;

                        // Known tag found; lineIn was something like this:
                        //   StaticMod=C2H3N1O1,C,fix,any,Carbamidomethylation
                        //   DynamicMod=C2H3NO, *,  opt, N-term,   Carbamidomethylation
                        //   CustomAA=C5H7N1O2S0,J,custom,P,Hydroxylation     # Hydroxyproline
                        //   StaticCysCarbamidomethyl=True

                        // And modSpec will now be something like this:
                        //   C2H3N1O1,C,fix,any,Carbamidomethylation
                        //   C2H3NO, *,  opt, N-term,   Carbamidomethylation
                        //   C5H7N1O2S0,J,custom,P,Hydroxylation     # Hydroxyproline
                        //   True


                        if (tagName.Equals(PARAM_TAG_MOD_STATIC, StringComparison.OrdinalIgnoreCase))
                        {
                            modType = MSGFPlusModType.StaticMod;
                        }
                        else if (tagName.Equals(PARAM_TAG_MOD_DYNAMIC, StringComparison.OrdinalIgnoreCase))
                        {
                            modType = MSGFPlusModType.DynamicMod;
                        }
                        else if (tagName.Equals(PARAM_TAG_CUSTOM_AA, StringComparison.OrdinalIgnoreCase))
                        {
                            modType = MSGFPlusModType.CustomAA;
                        }
                        else if (tagName.Equals(PARAM_TAG_STATIC_CYS_CARBAMIDOMETHYL, StringComparison.OrdinalIgnoreCase))
                        {
                            if (bool.TryParse(modSpec, out var cysAlkEnabled))
                            {
                                if (cysAlkEnabled)
                                {
                                    // The parameter file has "StaticCysCarbamidomethyl=True"
                                    modType = MSGFPlusModType.StaticMod;

                                    // Change modSpec from "True" to the DIA-NN style modification info for static sys Carbamidomethyl
                                    modSpec = "UniMod:4, 57.021465, C";
                                }
                                else
                                {
                                    // The parameter file has "StaticCysCarbamidomethyl=False"
                                    modSpec = string.Empty;
                                }
                            }
                            else
                            {
                                ReportWarning(string.Format(
                                    "{0} is not followed by True or False in the DIA-NN parameter file: {1}",
                                    PARAM_TAG_STATIC_CYS_CARBAMIDOMETHYL, dataLine));

                                modSpec = string.Empty;
                            }
                        }
                        else
                        {
                            modType = MSGFPlusModType.Unknown;
                        }

                        break;
                    }

                    if (string.IsNullOrEmpty(modSpec))
                    {
                        // The line does not start with StaticMod, DynamicMod, CustomAA, or StaticCysCarbamidomethyl
                        // This method also supports MSGFPlus_Mods.txt files, which specify mods with ,opt, or ,fix,
                        var lineInNoSpaces = TrimComment(trimmedLine).Replace(" ", string.Empty);

                        if (lineInNoSpaces.Contains(",opt,"))
                        {
                            modType = MSGFPlusModType.DynamicMod;
                            modSpec = lineInNoSpaces;
                        }
                        else if (lineInNoSpaces.Contains(",fix,"))
                        {
                            modType = MSGFPlusModType.StaticMod;
                            modSpec = lineInNoSpaces;
                        }
                        else if (lineInNoSpaces.Contains(",custom,"))
                        {
                            modType = MSGFPlusModType.CustomAA;
                            modSpec = lineInNoSpaces;
                        }
                    }

                    if (string.IsNullOrEmpty(modSpec))
                    {
                        continue;
                    }

                    if (modSpec.Contains("="))
                    {
                        // Unrecognized tag name
                        ReportError(string.Format(
                            "Mod spec '{0}' contains an unknown keyword before the equals sign; see parameter file {1}",
                            modSpec, Path.GetFileName(paramFilePath)));
                        return false;
                    }

                    // Modification definition line found
                    // Split the line on commas
                    var splitLine = modSpec.Split(',');

                    // Define the minimum number of parts in the mod spec
                    int minimumModDefParts;

                    switch (modSpecFormat)
                    {
                        case ModSpecFormats.MSGFPlusAndMSPathFinder:
                        case ModSpecFormats.TopPIC:
                            minimumModDefParts = 5;
                            break;

                        case ModSpecFormats.DiaNN:
                            minimumModDefParts = 3;
                            break;

                        default:
                            // Unrecognized format
                            ReportError(string.Format("Mod spec format {0} not recognized; unable to extract mods from {1}",
                                modSpecFormat.ToString(), Path.GetFileName(paramFilePath)));
                            return false;
                    }

                    if (splitLine.Length < minimumModDefParts)
                    {
                        continue;
                    }

                    switch (modSpecFormat)
                    {
                        case ModSpecFormats.MSGFPlusAndMSPathFinder:
                            if (ParseModSpecMSGFPlus(paramFilePath, splitLine, ref unnamedModID, out var udtMSGFPlusModInfo))
                            {
                                modList.Add(udtMSGFPlusModInfo);
                            }
                            break;

                        case ModSpecFormats.TopPIC:
                            if (ParseModSpecTopPIC(paramFilePath, splitLine, modType, ref unnamedModID, out var udtTopPICModInfo))
                            {
                                modList.Add(udtTopPICModInfo);
                            }
                            break;

                        case ModSpecFormats.DiaNN:
                            if (ParseModSpecDiaNN(paramFilePath, splitLine, modType, ref unnamedModID, out var udtDiannModInfo))
                            {
                                modList.Add(udtDiannModInfo);
                            }
                            break;
                    }
                }

                Console.WriteLine();
            }
            catch (Exception ex)
            {
                ReportError(string.Format(
                                "Error reading mod info from the {0} parameter file ({1}): {2}",
                                mToolName, Path.GetFileName(paramFilePath), ex.Message));
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse the mod spec definition from a DIA-NN parameter file
        /// </summary>
        /// <param name="paramFilePath"></param>
        /// <param name="splitLine"></param>
        /// <param name="modType">MSGFPlusModType.DynamicMod or MSGFPlusModType.StaticMod</param>
        /// <param name="unnamedModID"></param>
        /// <param name="udtModInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseModSpecDiaNN(
            string paramFilePath,
            IReadOnlyList<string> splitLine,
            MSGFPlusModType modType,
            ref int unnamedModID,
            out ModInfo udtModInfo)
        {
            // Modification definition format:
            //   ModName:   UniMod ID (e.g. UniMod:35) or a custom name
            //   Mass:      Mod mass
            //   Residues:  affected residues *

            udtModInfo = new ModInfo();

            try
            {
                udtModInfo.ModMass = splitLine[1].Trim();

                // udtModInfo.ModMass should be a number
                if (!double.TryParse(udtModInfo.ModMass, out udtModInfo.ModMassVal))
                {
                    ReportWarning("Non-numeric mod mass in the " + mToolName + " parameter file: " + string.Join(",", splitLine));
                    return false;
                }

                udtModInfo.Residues = splitLine[2].Trim();
                udtModInfo.ModSymbol = UNKNOWN_MSGFPlus_MOD_SYMBOL;

                udtModInfo.ModType = modType;
                udtModInfo.ModTypeInParameterFile = udtModInfo.ModType;

                if (udtModInfo.Residues == "*n")
                {
                    // N-terminal protein mod
                    if (udtModInfo.ModType == MSGFPlusModType.StaticMod)
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                    }

                    udtModInfo.Residues = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                    if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                        udtModInfo.ModType = MSGFPlusModType.DynNTermPeptide;
                }

                udtModInfo.ModName = ParseModSpecGetName(splitLine[0], ref unnamedModID);
                udtModInfo.ShortName = string.Empty;

                return true;
            }
            catch (Exception ex)
            {
                ReportError(string.Format(
                                "Error reading mod info from the {0} parameter file ({1}): {2}",
                                mToolName, Path.GetFileName(paramFilePath), ex.Message));

                return false;
            }
        }

        /// <summary>
        /// Parse the mod spec definition from a MS-GF+ or MSPathFinder parameter file
        /// </summary>
        /// <param name="paramFilePath"></param>
        /// <param name="splitLine"></param>
        /// <param name="unnamedModID"></param>
        /// <param name="udtModInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseModSpecMSGFPlus(
            string paramFilePath,
            IReadOnlyList<string> splitLine,
            ref int unnamedModID,
            out ModInfo udtModInfo)
        {
            // Modification definition format:
            //   Mass or EmpiricalFormula:   Mod mass, or empirical formula (like C2H3N1O1 or H-2O-1)
            //   Residues:  affected residues or N/C-terminus or *
            //   ModType:   fix or opt
            //   Position:  any, N-term, Prot-N-term, etc.; see MSGFPlusModType.CustomAA
            //   ModName:   UniMod name or custom name

            // Custom amino acids using a different format:
            //   EmpiricalFormula
            //   Symbol:     Single letter abbreviation for the custom amino acid: B, J, O, U, X, Z
            //   Type:       always "custom"
            //   Unused:     unused field, but a letter must be present
            //   Name:       Name associated with the custom amino acid

            udtModInfo = new ModInfo();

            try
            {
                udtModInfo.ModMass = splitLine[0].Trim();

                // udtModInfo.ModMass could be a number, or could be an empirical formula
                // First try to parse out a number
                if (!double.TryParse(udtModInfo.ModMass, out udtModInfo.ModMassVal))
                {
                    // Not a number
                    // Mod (or custom AA) is specified as an empirical formula
                    // Compute the mass
                    udtModInfo.ModMassVal = ComputeMass(udtModInfo.ModMass);
                }

                udtModInfo.Residues = splitLine[1].Trim();
                udtModInfo.ModSymbol = UNKNOWN_MSGFPlus_MOD_SYMBOL;

                switch (splitLine[2].Trim().ToLower())
                {
                    case "opt":
                        udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                        break;

                    case "fix":
                        udtModInfo.ModType = MSGFPlusModType.StaticMod;
                        break;

                    case "custom":
                        udtModInfo.ModType = MSGFPlusModType.CustomAA;
                        break;

                    default:
                        ReportWarning(string.Format(
                                          "Unrecognized Mod Type {0} in the {1} parameter file; should be 'opt', 'fix', or 'custom'; will assume 'opt'",
                                          splitLine[2], mToolName));

                        udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                        break;
                }

                udtModInfo.ModTypeInParameterFile = udtModInfo.ModType;

                if (udtModInfo.ModType != MSGFPlusModType.CustomAA)
                {
                    switch (splitLine[3].Trim().ToLower().Replace("-", string.Empty))
                    {
                        case "any":
                            // Leave .ModType unchanged; this is a static or dynamic mod (fix or opt)
                            break;

                        case "nterm":
                            if (udtModInfo.ModType == MSGFPlusModType.StaticMod && udtModInfo.Residues != "*")
                            {
                                // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                            }

                            udtModInfo.Residues = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                            if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                                udtModInfo.ModType = MSGFPlusModType.DynNTermPeptide;

                            break;

                        case "cterm":
                            if (udtModInfo.ModType == MSGFPlusModType.StaticMod && udtModInfo.Residues != "*")
                            {
                                // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                            }
                            udtModInfo.Residues = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                            if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                                udtModInfo.ModType = MSGFPlusModType.DynCTermPeptide;

                            break;

                        case "protnterm":
                            // Includes Prot-N-Term, Prot-n-Term, ProtNTerm, etc.
                            if (udtModInfo.ModType == MSGFPlusModType.StaticMod && udtModInfo.Residues != "*")
                            {
                                // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                            }
                            udtModInfo.Residues = AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                            if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                                udtModInfo.ModType = MSGFPlusModType.DynNTermProtein;

                            break;

                        case "protcterm":
                            // Includes Prot-C-Term, Prot-c-Term, ProtCterm, etc.
                            if (udtModInfo.ModType == MSGFPlusModType.StaticMod && udtModInfo.Residues != "*")
                            {
                                // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                            }
                            udtModInfo.Residues = AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                            if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                                udtModInfo.ModType = MSGFPlusModType.DynCTermProtein;

                            break;

                        default:
                            ReportWarning(string.Format(
                                              "Unrecognized Mod Position {0} in the {1} parameter file; " +
                                              "should be 'any', 'N-term', 'C-term', 'Prot-N-term', or 'Prot-C-term'",
                                              splitLine[3], mToolName));
                            break;
                    }
                }

                udtModInfo.ModName = ParseModSpecGetName(splitLine[4], ref unnamedModID);
                udtModInfo.ShortName = string.Empty;

                // ReSharper disable once ConvertSwitchStatementToSwitchExpression
                switch (udtModInfo.ModName)
                {
                    // ReSharper disable StringLiteralTypo

                    // UniMod names
                    case "iTRAQ4plex" or "iTRAQ8plex" or
                         "TMT" or "TMT2plex" or "TMT6plex" or "TMTpro" or "TMTpro_zero" or
                         "cysTMT" or "cysTMT6plex" or "iodoTMT" or "iodoTMT6plex":
                        udtModInfo.IsobaricMod = true;
                        break;

                    // ReSharper restore StringLiteralTypo

                    // DMS Names
                    case "itrac" or "iTRAQ8" or "TMT0Tag" or "TMT2Tag" or "TMT6Tag" or "TMT16Tag" or "iodoTMT6" or "TMT6Gly2":
                        udtModInfo.IsobaricMod = true;
                        break;

                    default:
                        udtModInfo.IsobaricMod = false;
                        break;
                }

                return true;
            }
            catch (Exception ex)
            {
                ReportError(string.Format(
                                "Error reading mod info from the {0} parameter file ({1}): {2}",
                                mToolName, Path.GetFileName(paramFilePath), ex.Message));
                return false;
            }
        }

        /// <summary>
        /// Parse the mod spec definition from a TopPIC parameter file
        /// </summary>
        /// <param name="paramFilePath"></param>
        /// <param name="splitLine"></param>
        /// <param name="modType">MSGFPlusModType.DynamicMod or MSGFPlusModType.StaticMod</param>
        /// <param name="unnamedModID"></param>
        /// <param name="udtModInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseModSpecTopPIC(
            string paramFilePath,
            IReadOnlyList<string> splitLine,
            MSGFPlusModType modType,
            ref int unnamedModID,
            out ModInfo udtModInfo)
        {
            // Modification definition format:
            //   ModName:   UniMod name or custom name
            //   Mass:      Mod mass
            //   Residues:  affected residues *
            //   Position:  any, N-term, or C-term
            //   UniModID:  UniMod ID, or -1 if not in UniMod

            udtModInfo = new ModInfo();

            try
            {
                udtModInfo.ModMass = splitLine[1].Trim();

                // udtModInfo.ModMass should be a number
                if (!double.TryParse(udtModInfo.ModMass, out udtModInfo.ModMassVal))
                {
                    ReportWarning("Non-numeric mod mass in the " + mToolName + " parameter file: " + string.Join(",", splitLine));
                    return false;
                }

                udtModInfo.Residues = splitLine[2].Trim();
                udtModInfo.ModSymbol = UNKNOWN_MSGFPlus_MOD_SYMBOL;

                udtModInfo.ModType = modType;
                udtModInfo.ModTypeInParameterFile = udtModInfo.ModType;

                switch (splitLine[3].Trim().ToLower().Replace("-", string.Empty))
                {
                    case "any":
                        // Leave .ModType unchanged; this is a static or dynamic mod
                        break;

                    case "nterm":
                        if (udtModInfo.ModType == MSGFPlusModType.StaticMod && udtModInfo.Residues != "*")
                        {
                            // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                            udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                        }
                        udtModInfo.Residues = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                            udtModInfo.ModType = MSGFPlusModType.DynNTermPeptide;

                        break;

                    case "cterm":
                        if (udtModInfo.ModType == MSGFPlusModType.StaticMod && udtModInfo.Residues != "*")
                        {
                            // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                            udtModInfo.ModType = MSGFPlusModType.DynamicMod;
                        }
                        udtModInfo.Residues = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        if (udtModInfo.ModType == MSGFPlusModType.DynamicMod)
                            udtModInfo.ModType = MSGFPlusModType.DynCTermPeptide;

                        break;

                    default:
                        ReportWarning(string.Format(
                                          "Unrecognized Mod Position {0} in the {1} parameter file; " +
                                          "should be 'any', 'N-term', or 'C-term'",
                                          splitLine[3], mToolName));
                        break;
                }

                udtModInfo.ModName = ParseModSpecGetName(splitLine[0], ref unnamedModID);
                udtModInfo.ShortName = string.Empty;

                // splitLine[4] has UniModID

                return true;
            }
            catch (Exception ex)
            {
                ReportError(string.Format(
                                "Error reading mod info from the {0} parameter file ({1}): {2}",
                                mToolName, Path.GetFileName(paramFilePath), ex.Message));

                return false;
            }
        }

        private string ParseModSpecGetName(string modName, ref int unnamedModID)
        {
            if (string.IsNullOrWhiteSpace(modName))
            {
                unnamedModID++;
                return "UnnamedMod" + unnamedModID;
            }

            return modName.Trim();
        }

        private void ReportError(string message)
        {
            mErrorMessage = message;
            OnErrorEvent(message);
        }

        private void ReportWarning(string message)
        {
            OnWarningEvent(message);
        }

        /// <summary>
        /// Resolve MS-GF+, MSPathFinder, TopPIC, or DIA-NN mods with modification definitions
        /// </summary>
        /// <param name="modList"></param>
        /// <param name="peptideMods"></param>
        public void ResolveMSGFPlusModsWithModDefinitions(List<ModInfo> modList, PeptideModificationContainer peptideMods)
        {
            if (modList == null)
                return;

            for (var index = 0; index <= modList.Count - 1; index++)
            {
                var udtModInfo = modList[index];
                int resIndexStart;
                int resIndexEnd;

                if (udtModInfo.Residues.Length > 0)
                {
                    resIndexStart = 0;
                    resIndexEnd = udtModInfo.Residues.Length - 1;
                }
                else
                {
                    resIndexStart = -1;
                    resIndexEnd = -1;
                }

                for (var residueIndex = resIndexStart; residueIndex <= resIndexEnd; residueIndex++)
                {
                    char chTargetResidue;
                    if (residueIndex >= 0)
                    {
                        chTargetResidue = udtModInfo.Residues[residueIndex];
                        if (chTargetResidue == '*')
                        {
                            // This is a terminal mod, and MSGFDB lists the target residue as * for terminal mods
                            // This program requires that chTargetResidue be Nothing
                            chTargetResidue = default;
                        }
                    }
                    else
                    {
                        chTargetResidue = default;
                    }

                    var modType = ModificationDefinition.ResidueModificationType.DynamicMod;
                    AminoAcidModInfo.ResidueTerminusState residueTerminusState;

                    if (udtModInfo.ModType == MSGFPlusModType.DynNTermPeptide)
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                    }
                    else if (udtModInfo.ModType == MSGFPlusModType.DynCTermPeptide)
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                    }
                    else if (udtModInfo.ModType == MSGFPlusModType.DynNTermProtein)
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus;
                    }
                    else if (udtModInfo.ModType == MSGFPlusModType.DynCTermProtein)
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus;
                    }
                    else
                    {
                        switch (chTargetResidue)
                        {
                            case AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS:
                                residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                                if (udtModInfo.ModType == MSGFPlusModType.StaticMod)
                                    modType = ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod;
                                break;

                            case AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS:
                                residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                                if (udtModInfo.ModType == MSGFPlusModType.StaticMod)
                                    modType = ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod;
                                break;

                            case AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS:
                                residueTerminusState = AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus;
                                if (udtModInfo.ModType == MSGFPlusModType.StaticMod)
                                    modType = ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod;
                                break;

                            case AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS:
                                residueTerminusState = AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus;
                                if (udtModInfo.ModType == MSGFPlusModType.StaticMod)
                                    modType = ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod;
                                break;

                            default:
                                residueTerminusState = AminoAcidModInfo.ResidueTerminusState.None;
                                if (udtModInfo.ModType == MSGFPlusModType.StaticMod)
                                    modType = ModificationDefinition.ResidueModificationType.StaticMod;
                                break;
                        }
                    }

                    var modificationDefinition = peptideMods.LookupModificationDefinitionByMassAndModType(
                        udtModInfo.ModMassVal, modType, chTargetResidue, residueTerminusState, out _, true);

                    if (residueIndex == resIndexStart)
                    {
                        // Update the Mod Symbol
                        udtModInfo.ModSymbol = modificationDefinition.ModificationSymbol;
                    }
                }

                modList[index] = udtModInfo;
            }
        }

        private static string TrimComment(string value)
        {
            // Look for the MS-GF+ comment character
            var commentCharIndex = value.IndexOf(MSGFPLUS_COMMENT_CHAR, StringComparison.Ordinal);

            if (commentCharIndex > 0)
            {
                // Trim off the comment
                return value.Substring(0, commentCharIndex).Trim();
            }

            return value.Trim();
        }

        /// <param name="modTag"></param>
        /// <returns>Modification definition if found, otherwise an empty string</returns>
        private string ValidateIsValidModSpec(string lineIn, string modTag)
        {
            if (!lineIn.StartsWith(modTag, StringComparison.OrdinalIgnoreCase))
                return string.Empty;

            var kvSetting = SynFileReaderBaseClass.ParseKeyValueSetting(lineIn, '=', "#");

            if (string.IsNullOrEmpty(kvSetting.Value) || string.Equals(kvSetting.Value, "none", StringComparison.OrdinalIgnoreCase))
            {
                // Not a valid mod spec
                return string.Empty;
            }

            // Mod spec found
            // Note that there may be spaces before or after the commas separating the mod spec fields
            return kvSetting.Value;
        }
    }
}
