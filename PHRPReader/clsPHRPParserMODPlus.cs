//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 05/18/2015
//
// This class parses data lines from modp_syn.txt files
//
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for MOD+
    /// </summary>
    public class clsPHRPParserMODPlus : clsPHRPParser
    {
        #region "Constants"

#pragma warning disable 1591
        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Spectrum_Index = "Spectrum_Index";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_NTT = "NTT";
        public const string DATA_COLUMN_Modification_Annotation = "Modification_Annotation";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Peptide_Position = "Peptide_Position";
        public const string DATA_COLUMN_Score = "Score";
        public const string DATA_COLUMN_Probability = "Probability";
        public const string DATA_COLUMN_Rank_Score = "Rank_Score";
        public const string DATA_COLUMN_QValue = "QValue";

        public const string FILENAME_SUFFIX_SYN = "_modp_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_modp_fht.txt";

        private const string MODPlus_SEARCH_ENGINE_NAME = "MODPlus";


        /// <summary>
        /// These columns correspond to the Synopsis file created by clsMODPlusResultsProcessor
        /// </summary>
        public enum MODPlusSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Spectrum_Index = 2,
            Charge = 3,
            PrecursorMZ = 4,
            DelM = 5,                            // Precursor error, in Da
            DelM_PPM = 6,                        // Precursor error, in ppm
            MH = 7,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
            Peptide = 8,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
            NTT = 9,
            ModificationAnnotation = 10,
            Protein = 11,
            Peptide_Position = 12,
            Score = 13,
            Probability = 14,
            Rank_Score = 15,
            QValue = 16
        }

#pragma warning restore 1591


        #endregion

        #region "Properties"

        /// <summary>
        /// First hits file
        /// </summary>
        public override string PHRPFirstHitsFileName => GetPHRPFirstHitsFileName(mDatasetName);

        /// <summary>
        /// Mod summary file
        /// </summary>
        public override string PHRPModSummaryFileName => GetPHRPModSummaryFileName(mDatasetName);

        /// <summary>
        /// Peptide to protein map file
        /// </summary>
        public override string PHRPPepToProteinMapFileName => GetPHRPPepToProteinMapFileName(mDatasetName);

        /// <summary>
        /// Protein mods file
        /// </summary>
        public override string PHRPProteinModsFileName => GetPHRPProteinModsFileName(mDatasetName);

        /// <summary>
        /// Synopsis file
        /// </summary>
        public override string PHRPSynopsisFileName => GetPHRPSynopsisFileName(mDatasetName);

        /// <summary>
        /// Result to sequence map file
        /// </summary>
        public override string PHRPResultToSeqMapFileName => GetPHRPResultToSeqMapFileName(mDatasetName);

        /// <summary>
        /// Sequence info file
        /// </summary>
        public override string PHRPSeqInfoFileName => GetPHRPSeqInfoFileName(mDatasetName);

        /// <summary>
        /// Sequence to protein map file
        /// </summary>
        public override string PHRPSeqToProteinMapFileName => GetPHRPSeqToProteinMapFileName(mDatasetName);

        /// <summary>
        /// Search engine name
        /// </summary>
        public override string SearchEngineName => GetSearchEngineName();

        #endregion

        /// <summary>
        /// Constructor; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <remarks></remarks>
        public clsPHRPParserMODPlus(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        /// <remarks></remarks>
        public clsPHRPParserMODPlus(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MODPlus, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMODPlus(string datasetName, string inputFilePath, clsPHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MODPlus, startupOptions)
        {
        }

        /// <summary>
        /// Get the header names in the PHRP synopsis or first hits file for this tool
        /// </summary>
        /// <returns></returns>
        protected override List<string> GetColumnHeaderNames()
        {
            var headerNames = new List<string>();
            headerNames.AddRange(GetColumnHeaderNamesAndIDs().Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <returns></returns>
        public static SortedDictionary<string, MODPlusSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            var headerColumns = new SortedDictionary<string, MODPlusSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {DATA_COLUMN_ResultID, MODPlusSynFileColumns.ResultID},
                {DATA_COLUMN_Scan, MODPlusSynFileColumns.Scan},
                {DATA_COLUMN_Spectrum_Index, MODPlusSynFileColumns.Spectrum_Index},
                {DATA_COLUMN_Charge, MODPlusSynFileColumns.Charge},
                {DATA_COLUMN_PrecursorMZ, MODPlusSynFileColumns.PrecursorMZ},
                {DATA_COLUMN_DelM, MODPlusSynFileColumns.DelM},
                {DATA_COLUMN_DelM_PPM, MODPlusSynFileColumns.DelM_PPM},
                {DATA_COLUMN_MH, MODPlusSynFileColumns.MH},
                {DATA_COLUMN_Peptide, MODPlusSynFileColumns.Peptide},
                {DATA_COLUMN_NTT, MODPlusSynFileColumns.NTT},
                {DATA_COLUMN_Modification_Annotation, MODPlusSynFileColumns.ModificationAnnotation},
                {DATA_COLUMN_Protein, MODPlusSynFileColumns.Protein},
                {DATA_COLUMN_Peptide_Position, MODPlusSynFileColumns.Peptide_Position},
                {DATA_COLUMN_Score, MODPlusSynFileColumns.Score},
                {DATA_COLUMN_Probability, MODPlusSynFileColumns.Probability},
                {DATA_COLUMN_Rank_Score, MODPlusSynFileColumns.Rank_Score},
                {DATA_COLUMN_QValue, MODPlusSynFileColumns.QValue}
            };

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MODPlusSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MODPlusSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // MODPlus does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_modp_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_modp_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_modp_syn_ProteinMods.txt";
        }

        /// <summary>
        /// Default Synopsis file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_SYN;
        }

        /// <summary>
        /// Default ResultToSeq map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(string datasetName)
        {
            return datasetName + "_modp_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_modp_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_modp_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MODPlus_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MODp parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out clsSearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new clsSearchEngineParameters(MODPlus_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, clsSearchEngineParameters searchEngineParams)
        {
            try
            {
                // For MS-GF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                // But MODPlus does not have a _ModDefs.txt file because it performs a blind search
                // The user can define static mods on any of the residues, plus the peptide terminii; check for these now

                var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

                // Read the contents of the parameter file
                var doc = new XmlDocument();
                doc.Load(paramFilePath);

                var dbNodes = doc.SelectNodes("/search/database");
                if (dbNodes != null && dbNodes.Count > 0)
                {
                    searchEngineParams.FastaFilePath = GetAttribute(dbNodes[0], "local_path");
                }

                var enzymeNodes = doc.SelectNodes("/search/enzyme_rule");
                if (enzymeNodes != null && enzymeNodes.Count > 0)
                {
                    searchEngineParams.Enzyme = GetAttribute(enzymeNodes[0], "name");
                }

                var instrumentResolutionNodes = doc.SelectNodes("/search/instrument_resolution");
                if (instrumentResolutionNodes != null && instrumentResolutionNodes.Count > 0)
                {
                    searchEngineParams.PrecursorMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes[0], "ms"));

                    // ReSharper disable once StringLiteralTypo
                    searchEngineParams.FragmentMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes[0], "msms"));
                }

                var enzymeConstraintNodes = doc.SelectNodes("/search/parameters/enzyme_constraint");
                if (enzymeConstraintNodes != null && enzymeConstraintNodes.Count > 0)
                {
                    if (int.TryParse(GetAttribute(enzymeConstraintNodes[0], "max_miss_cleavages"), out var maxNumberInternalCleavages))
                    {
                        searchEngineParams.MaxNumberInternalCleavages = maxNumberInternalCleavages;
                    }
                    if (int.TryParse(GetAttribute(enzymeConstraintNodes[0], "min_number_termini"), out var minNumberTermini))
                    {
                        searchEngineParams.MinNumberTermini = minNumberTermini;
                    }
                }

                var massTolParamNode = doc.SelectNodes("/search/parameters/peptide_mass_tol");

                if (massTolParamNode != null && massTolParamNode.Count > 0)
                {
                    var tolerance = GetAttribute(massTolParamNode[0], "value");
                    var massUnits = GetAttribute(massTolParamNode[0], "unit");

                    if (massUnits.ToLower() == "ppm")
                    {
                        // Parent mass tolerance, in ppm

                        if (double.TryParse(tolerance, out var precursorMassTolerancePpm))
                        {
                            searchEngineParams.PrecursorMassTolerancePpm = precursorMassTolerancePpm;
                            searchEngineParams.PrecursorMassToleranceDa = clsPeptideMassCalculator.PPMToMass(searchEngineParams.PrecursorMassTolerancePpm, 2000);
                        }
                    }
                    else if (massUnits.ToLower() == "da")
                    {
                        if (double.TryParse(tolerance, out var precursorMassToleranceDa))
                        {
                            searchEngineParams.PrecursorMassToleranceDa = precursorMassToleranceDa;
                        }

                        // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                        searchEngineParams.PrecursorMassToleranceDa = clsPeptideMassCalculator.MassToPPM(searchEngineParams.PrecursorMassToleranceDa, 2000);
                    }
                }

                var modNodes = doc.SelectNodes("/search/modifications/fixed/mod");
                if (modNodes == null || modNodes.Count <= 0)
                    return true;

                // Store the fixed mods

                foreach (XmlNode node in modNodes)
                {
                    // var modName = GetAttribute(node, "name");
                    var residue = GetAttribute(node, "site").Trim();
                    // var modPosition = GetAttribute(node, "position");

                    // ReSharper disable once StringLiteralTypo
                    var modMass = GetAttribute(node, "massdiff");

                    // Replace N-Term or C-Term with < or >
                    if (residue.ToLower() == "n-term")
                        residue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                    if (residue.ToLower() == "c-term")
                        residue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                    if (!double.TryParse(modMass, out var modMassDa))
                        continue;

                    if (Math.Abs(modMassDa - 0) < float.Epsilon)
                        continue;

                    var eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod;
                    if (residue == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || residue == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                    }

                    var modDef = new clsModificationDefinition(
                        clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                        modMassDa,
                        residue,
                        eModType,
                        "Mod" + modMassDa.ToString("0"));

                    searchEngineParams.AddModification(modDef);
                }

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }

        private string ConvertResolutionModeToMassType(string resolutionType)
        {
            if (resolutionType == "high")
            {
                return clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
            }

            if (resolutionType == "low")
            {
                return clsSearchEngineParameters.MASS_TYPE_AVERAGE;
            }

            return "unknown";
        }

        private string GetAttribute(XmlNode node, string attributeName)
        {
            if (node.Attributes != null && node.Attributes.Count > 0)
            {
                try
                {
                    var attribute = node.Attributes[attributeName];

                    if (attribute != null)
                    {
                        return attribute.Value;
                    }
                }
                catch (Exception ex)
                {
                    // Ignore errors
                    Console.WriteLine("Attribute lookup error: " + ex.Message);
                }
            }

            return string.Empty;
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out clsPSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;

            var columns = line.Split('\t');

            var success = false;

            psm = new clsPSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_ResultID, mColumnHeaders, 0);
                    psm.ScoreRank = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Rank_Score, mColumnHeaders, 1);

                    var peptide = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Peptide, mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                    var protein = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Protein, mColumnHeaders);
                    psm.AddProtein(protein);

                    var precursorMZ = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                    psm.MassErrorDa = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM, mColumnHeaders);
                    psm.MassErrorPPM = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, DATA_COLUMN_Spectrum_Index);

                    AddScore(psm, columns, DATA_COLUMN_MH);

                    AddScore(psm, columns, DATA_COLUMN_Modification_Annotation);
                    AddScore(psm, columns, DATA_COLUMN_Peptide_Position);

                    AddScore(psm, columns, DATA_COLUMN_Score);
                    AddScore(psm, columns, DATA_COLUMN_Probability);
                    AddScore(psm, columns, DATA_COLUMN_QValue);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MODPlus data file: " + ex.Message);
            }

            return success;
        }
    }
}
