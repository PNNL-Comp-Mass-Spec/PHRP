// This class reads a DMS-based parameter file for MSGF+ or MSPathFinder
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
// Note that DMS uses this information to create a Mods.txt file that is provided to MSGF+ or MSPathFinder
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
using System.Linq;

namespace PHRPReader
{
    /// <summary>
    /// This class reads a DMS-based parameter file for MSGF+ or MSPathFinder to extract the dynamic and static modification information
    /// </summary>
    /// <remarks>See above for an example parameter file</remarks>
    public class clsMSGFPlusParamFileModExtractor : PRISM.EventNotifier
    {
        #region "Constants and Enums"

        /// <summary>
        /// Unknown MSGF+ mod symbols
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

        private const string MSGFPLUS_COMMENT_CHAR = "#";

        /// <summary>
        /// MSGF+ modification type
        /// </summary>
        public enum eMSGFDBModType
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

        #endregion

        #region "Structures"

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
        ///   ModType:    eMSGFDBModType.CustomAA
        ///   ModSymbol:  ?   (a question mark; not used)
        /// </remarks>
        public struct udtModInfoType
        {
            /// <summary>
            /// Mod name (read from the parameter file) isn't used by MSGF+, but it is used by MSPathFinder
            /// </summary>
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
            public string Residues;

            /// <summary>
            /// Modification type
            /// </summary>
            public eMSGFDBModType ModType;

            /// <summary>
            /// Modification symbol: *, #, @, ... ; dash if a static mod
            /// </summary>
            public char ModSymbol;

            /// <summary>
            /// Mod type, name, mass, residues
            /// </summary>
            /// <returns></returns>
            public override string ToString()
            {
                return ModType + " " + ModName + ", " + ModMass + "; " + Residues;
            }
        }
        #endregion

        #region "Classwide Variables"
        private string mErrorMessage;
        private readonly string mToolName;
        #endregion

        #region "Properties"

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage => mErrorMessage;

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="toolName">
        /// Search engine name, typically MSGF+
        /// This name is only used in log messages
        /// </param>
        /// <remarks></remarks>
        public clsMSGFPlusParamFileModExtractor(string toolName)
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

            if (string.Equals(empiricalFormula, "HexNAc", StringComparison.InvariantCultureIgnoreCase))
            {
                // This is a special-case modification that MSGF+ and MSPathFinder recognize
                // It is listed in DMS as Hexosam, which means N-Acetylhexosamine
                // It is tracked by UniMod as HexNAc
                return 203.079376;
            }

            clsEmpiricalFormula empiricalFormulaInstance;

            try
            {
                empiricalFormulaInstance = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormula);
            }
            catch (Exception ex)
            {
                ReportError(ex.Message);
                return 0;
            }

            var monoisotopicMass = clsPeptideMassCalculator.ComputeMonoisotopicMass(empiricalFormulaInstance.ElementCounts, out var unknownSymbols);

            if (unknownSymbols != null && unknownSymbols.Count > 0)
            {
                var errMsg = "Error parsing empirical formula '" + empiricalFormula + "', ";
                if (unknownSymbols.Count == 1)
                {
                    ReportError(errMsg + "unknown element " + unknownSymbols.First());
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
        /// Extracts mod info from either a MSGF+ or MSPathFinder param file or from a MSGFPlus_Mods.txt file (previously MSGFDB_Mods.txt)
        /// </summary>
        /// <param name="paramFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True if success; false if a problem</returns>
        /// <remarks></remarks>
        public bool ExtractModInfoFromParamFile(string paramFilePath, out List<udtModInfoType> modInfo)
        {
            var tagNamesToFind = new List<string> {
                PARAM_TAG_MOD_STATIC,
                PARAM_TAG_MOD_DYNAMIC,
                PARAM_TAG_CUSTOM_AA };

            // Initialization
            modInfo = new List<udtModInfoType>();

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
                using (var reader = new StreamReader(new FileStream(paramFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
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

                        foreach (var tagName in tagNamesToFind)
                        {
                            modSpec = ValidateIsValidModSpec(trimmedLine, tagName);
                            if (!string.IsNullOrEmpty(modSpec))
                            {
                                // Known tag found; lineIn was something like this:
                                //   StaticMod=C2H3N1O1,C,fix,any,Carbamidomethylation
                                //   DynamicMod=C2H3NO, *,  opt, N-term,   Carbamidomethylation
                                //   CustomAA=C5H7N1O2S0,J,custom,P,Hydroxylation     # Hydroxyproline
                                //
                                // And modSpec will now be something like this:
                                //   C2H3N1O1,C,fix,any,Carbamidomethylation
                                //   C2H3NO, *,  opt, N-term,   Carbamidomethylation
                                //   C5H7N1O2S0,J,custom,P,Hydroxylation     # Hydroxyproline

                                break;
                            }
                        }

                        if (string.IsNullOrEmpty(modSpec))
                        {
                            var lineInNoSpaces = TrimComment(trimmedLine).Replace(" ", "");
                            if (lineInNoSpaces.Contains(",opt,") || lineInNoSpaces.Contains(",fix,") || lineInNoSpaces.Contains(",custom,"))
                            {
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
                            ReportError("Mod spec '" + modSpec + "' contains an unknown keyword before the equals sign; see parameter file " + Path.GetFileName(paramFilePath));
                            return false;
                        }
                        // Modification definition line found

                        // Split the line on commas
                        var splitLine = modSpec.Split(',');

                        if (splitLine.Length < 5)
                        {
                            continue;
                        }


                        // Notes when tracking information for custom amino acids
                        //   ModName:    Name associated with the custom amino acid
                        //   ModMass:    Composition string, for example C5H7N1O2S0 for Hydroxyproline
                        //   ModMassVal: Computed mass of the composition string
                        //   Residues:   Single letter abbreviation for the custom amino acid, for example J or X
                        //   ModType:    eMSGFDBModType.CustomAA
                        //   ModSymbol:  ?   (a question mark; not used)

                        var udtModInfo = new udtModInfoType {
                            ModMass = splitLine[0].Trim()
                        };

                        // .ModMass could be a number, or could be an empirical formula
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
                                udtModInfo.ModType = eMSGFDBModType.DynamicMod;
                                break;
                            case "fix":
                                udtModInfo.ModType = eMSGFDBModType.StaticMod;
                                break;
                            case "custom":
                                udtModInfo.ModType = eMSGFDBModType.CustomAA;
                                break;
                            default:
                                ReportWarning("Unrecognized Mod Type in the " + mToolName + " parameter file; should be 'opt' or 'fix'");
                                udtModInfo.ModType = eMSGFDBModType.DynamicMod;
                                break;
                        }

                        if (udtModInfo.ModType != eMSGFDBModType.CustomAA)
                        {
                            switch (splitLine[3].Trim().ToLower().Replace("-", string.Empty))
                            {
                                case "any":
                                    break;
                                // Leave .ModType unchanged; this is a static or dynamic mod (fix or opt)
                                case "nterm":
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod && udtModInfo.Residues != "*")
                                    {
                                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        udtModInfo.ModType = eMSGFDBModType.DynamicMod;
                                    }
                                    udtModInfo.Residues = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                                    if (udtModInfo.ModType == eMSGFDBModType.DynamicMod)
                                        udtModInfo.ModType = eMSGFDBModType.DynNTermPeptide;

                                    break;
                                case "cterm":
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod && udtModInfo.Residues != "*")
                                    {
                                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        udtModInfo.ModType = eMSGFDBModType.DynamicMod;
                                    }
                                    udtModInfo.Residues = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                                    if (udtModInfo.ModType == eMSGFDBModType.DynamicMod)
                                        udtModInfo.ModType = eMSGFDBModType.DynCTermPeptide;

                                    break;
                                case "protnterm":
                                    // Includes Prot-N-Term, Prot-n-Term, ProtNTerm, etc.
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod && udtModInfo.Residues != "*")
                                    {
                                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        udtModInfo.ModType = eMSGFDBModType.DynamicMod;
                                    }
                                    udtModInfo.Residues = clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                                    if (udtModInfo.ModType == eMSGFDBModType.DynamicMod)
                                        udtModInfo.ModType = eMSGFDBModType.DynNTermProtein;

                                    break;
                                case "protcterm":
                                    // Includes Prot-C-Term, Prot-c-Term, ProtCterm, etc.
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod && udtModInfo.Residues != "*")
                                    {
                                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        udtModInfo.ModType = eMSGFDBModType.DynamicMod;
                                    }
                                    udtModInfo.Residues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                                    if (udtModInfo.ModType == eMSGFDBModType.DynamicMod)
                                        udtModInfo.ModType = eMSGFDBModType.DynCTermProtein;

                                    break;
                                default:
                                    ReportWarning("Unrecognized Mod Type in the " + mToolName + " parameter file; should be 'any', 'N-term', 'C-term', 'Prot-N-term', or 'Prot-C-term'");
                                    break;
                            }
                        }

                        udtModInfo.ModName = splitLine[4].Trim();
                        if (string.IsNullOrEmpty(udtModInfo.ModName))
                        {
                            unnamedModID += 1;
                            udtModInfo.ModName = "UnnamedMod" + unnamedModID.ToString();
                        }

                        modInfo.Add(udtModInfo);
                    }
                }

                Console.WriteLine();
            }
            catch (Exception ex)
            {
                ReportError("Error reading mod info the " + mToolName + " parameter file (" + Path.GetFileName(paramFilePath) + "): " + ex.Message);
                return false;
            }

            return true;
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
        /// Resolve MSGF+ mods with mod definitions
        /// </summary>
        /// <param name="modInfo"></param>
        /// <param name="oPeptideMods"></param>
        public void ResolveMSGFDBModsWithModDefinitions(List<udtModInfoType> modInfo, clsPeptideModificationContainer oPeptideMods)
        {
            if (modInfo != null)
            {
                // Call .LookupModificationDefinitionByMass for each entry in msgfdbModInfo

                for (var index = 0; index <= modInfo.Count - 1; index++)
                {
                    var udtModInfo = modInfo[index];
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
                                chTargetResidue = default(char);
                            }
                        }
                        else
                        {
                            chTargetResidue = default(char);
                        }

                        var eModType = clsModificationDefinition.eModificationTypeConstants.DynamicMod;
                        clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState;

                        if (udtModInfo.ModType == eMSGFDBModType.DynNTermPeptide)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                        }
                        else if (udtModInfo.ModType == eMSGFDBModType.DynCTermPeptide)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                        }
                        else if (udtModInfo.ModType == eMSGFDBModType.DynNTermProtein)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                        }
                        else if (udtModInfo.ModType == eMSGFDBModType.DynCTermProtein)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                        }
                        else
                        {
                            switch (chTargetResidue)
                            {
                                case clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS:
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod)
                                        eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                                    break;
                                case clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS:
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod)
                                        eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                                    break;
                                case clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS:
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod)
                                        eModType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod;
                                    break;
                                case clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS:
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod)
                                        eModType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod;
                                    break;
                                default:
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                                    if (udtModInfo.ModType == eMSGFDBModType.StaticMod)
                                        eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod;
                                    break;
                            }
                        }

                        var modificationDefinition = oPeptideMods.LookupModificationDefinitionByMassAndModType(
                            udtModInfo.ModMassVal, eModType, chTargetResidue, eResidueTerminusState, out _, true);

                        if (residueIndex == resIndexStart)
                        {
                            // Update the Mod Symbol
                            udtModInfo.ModSymbol = modificationDefinition.ModificationSymbol;
                        }
                    }

                    modInfo[index] = udtModInfo;
                }
            }
        }

        private static string TrimComment(string value)
        {
            // Look for the MSGF+ comment character
            var commentCharIndex = value.IndexOf(MSGFPLUS_COMMENT_CHAR, StringComparison.Ordinal);

            if (commentCharIndex > 0)
            {
                // Trim off the comment
                return value.Substring(0, commentCharIndex).Trim();
            }

            return value.Trim();
        }

        private string ValidateIsValidModSpec(string lineIn, string modTag)
        {
            var modSpec = string.Empty;

            if (lineIn.StartsWith(modTag, StringComparison.InvariantCultureIgnoreCase))
            {
                var kvSetting = clsPHRPParser.ParseKeyValueSetting(lineIn, '=', "#");

                if (string.IsNullOrEmpty(kvSetting.Value) || kvSetting.Value.ToLower() == "none")
                {
                    // Not a valid mod spec
                    modSpec = string.Empty;
                }
                else
                {
                    // Mod spec found
                    // Note that there may be spaces before or after the commas separating the mod spec fields
                    modSpec = kvSetting.Value;
                }
            }

            return modSpec;
        }
    }
}
