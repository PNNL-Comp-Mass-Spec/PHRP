//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 03/19/2019
//
// This class parses data lines from toppic_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP parser for TopPIC
    /// </summary>
    public class TopPICSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: toppic, prsm, Frag, Da, Prot

        #region "Constants and Enums"

#pragma warning disable 1591

        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Prsm_ID = "Prsm_ID";
        public const string DATA_COLUMN_Spectrum_ID = "Spectrum_ID";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_Proteoform_ID = "Proteoform_ID";
        public const string DATA_COLUMN_Feature_Intensity = "Feature_Intensity";
        public const string DATA_COLUMN_Feature_Score = "Feature_Score";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_ResidueStart = "ResidueStart";
        public const string DATA_COLUMN_ResidueEnd = "ResidueEnd";
        public const string DATA_COLUMN_Unexpected_Mod_Count = "Unexpected_Mod_Count";
        public const string DATA_COLUMN_Peak_Count = "Peak_Count";
        public const string DATA_COLUMN_Matched_Peak_Count = "Matched_Peak_Count";
        public const string DATA_COLUMN_Matched_Fragment_Ion_Count = "Matched_Fragment_Ion_Count";
        public const string DATA_COLUMN_PValue = "PValue";
        public const string DATA_COLUMN_Rank_PValue = "Rank_PValue";
        public const string DATA_COLUMN_EValue = "EValue";
        public const string DATA_COLUMN_QValue = "QValue";
        public const string DATA_COLUMN_Proteoform_FDR = "Proteoform_FDR";
        public const string DATA_COLUMN_Proteoform_QValue = "Proteoform_QValue";
        public const string DATA_COLUMN_FragMethod = "FragMethod";
        public const string DATA_COLUMN_Variable_PTMs = "Variable_PTMs";

        public const string FILENAME_SUFFIX_SYN = "_toppic_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_toppic_fht.txt";

        private const string TopPIC_SEARCH_ENGINE_NAME = "TopPIC";

        /// <summary>
        /// These columns correspond to the Synopsis file created by TopPICResultsProcessor
        /// </summary>
        public enum TopPICSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Prsm_ID = 2,
            Spectrum_ID = 3,
            FragMethod = 4,
            Charge = 5,
            PrecursorMZ = 6,
            DelM = 7,                            // Precursor error, in Da
            DelMPPM = 8,                         // Precursor error, in ppm
            MH = 9,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
            Peptide = 10,                        // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
            Proteoform_ID = 11,
            Feature_Intensity = 12,
            Feature_Score = 13,
            Protein = 14,                        // Protein Name
            ResidueStart = 15,
            ResidueEnd = 16,
            Unexpected_Mod_Count = 17,
            Peak_Count = 18,
            Matched_Peak_Count = 19,
            Matched_Fragment_Ion_Count = 20,
            PValue = 21,
            Rank_PValue = 22,
            EValue = 23,
            QValue = 24,
            Proteoform_QValue = 25,
            VariablePTMs = 26
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
        public TopPICSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public TopPICSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.TopPIC, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public TopPICSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.TopPIC, startupOptions)
        {
        }

        /// <summary>
        /// Get the header names in the PHRP synopsis or first hits file for this tool
        /// </summary>
        /// <returns>List of header names</returns>
        protected override List<string> GetColumnHeaderNames()
        {
            var headerNames = new List<string>();
            headerNames.AddRange(GetColumnHeaderNamesAndIDs(true).Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, TopPICSynFileColumns> GetColumnHeaderNamesAndIDs(bool includeLegacyNames)
        {
            var headerColumns = new SortedDictionary<string, TopPICSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {DATA_COLUMN_ResultID, TopPICSynFileColumns.ResultID},
                {DATA_COLUMN_Scan, TopPICSynFileColumns.Scan},
                {DATA_COLUMN_Prsm_ID, TopPICSynFileColumns.Prsm_ID},
                {DATA_COLUMN_Spectrum_ID, TopPICSynFileColumns.Spectrum_ID},
                {DATA_COLUMN_FragMethod, TopPICSynFileColumns.FragMethod},
                {DATA_COLUMN_Charge, TopPICSynFileColumns.Charge},
                {DATA_COLUMN_PrecursorMZ, TopPICSynFileColumns.PrecursorMZ},
                {DATA_COLUMN_DelM, TopPICSynFileColumns.DelM},
                {DATA_COLUMN_DelM_PPM, TopPICSynFileColumns.DelMPPM},
                {DATA_COLUMN_MH, TopPICSynFileColumns.MH},
                {DATA_COLUMN_Peptide, TopPICSynFileColumns.Peptide},
                {DATA_COLUMN_Proteoform_ID, TopPICSynFileColumns.Proteoform_ID},
                {DATA_COLUMN_Feature_Intensity, TopPICSynFileColumns.Feature_Intensity},
                {DATA_COLUMN_Feature_Score, TopPICSynFileColumns.Feature_Score},
                {DATA_COLUMN_Protein, TopPICSynFileColumns.Protein},
                {DATA_COLUMN_ResidueStart, TopPICSynFileColumns.ResidueStart},
                {DATA_COLUMN_ResidueEnd, TopPICSynFileColumns.ResidueEnd},
                {DATA_COLUMN_Unexpected_Mod_Count, TopPICSynFileColumns.Unexpected_Mod_Count},
                {DATA_COLUMN_Peak_Count, TopPICSynFileColumns.Peak_Count},
                {DATA_COLUMN_Matched_Peak_Count, TopPICSynFileColumns.Matched_Peak_Count},
                {DATA_COLUMN_Matched_Fragment_Ion_Count, TopPICSynFileColumns.Matched_Fragment_Ion_Count},
                {DATA_COLUMN_PValue, TopPICSynFileColumns.PValue},
                {DATA_COLUMN_Rank_PValue, TopPICSynFileColumns.Rank_PValue},
                {DATA_COLUMN_EValue, TopPICSynFileColumns.EValue},
                {DATA_COLUMN_QValue, TopPICSynFileColumns.QValue}
            };

            if (!includeLegacyNames)
            {
                headerColumns.Add(DATA_COLUMN_Proteoform_QValue, TopPICSynFileColumns.Proteoform_QValue);
                headerColumns.Add(DATA_COLUMN_Variable_PTMs, TopPICSynFileColumns.VariablePTMs);

                return headerColumns;
            }

            headerColumns.Add(DATA_COLUMN_Proteoform_FDR, TopPICSynFileColumns.Proteoform_QValue);
            headerColumns.Add(DATA_COLUMN_Proteoform_QValue, TopPICSynFileColumns.Proteoform_QValue);
            headerColumns.Add(DATA_COLUMN_Variable_PTMs, TopPICSynFileColumns.VariablePTMs);

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum TopPICSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<TopPICSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs(true);
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_FHT;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            // ReSharper disable once StringLiteralTypo
            return datasetName + "_toppic_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            // ReSharper disable once StringLiteralTypo
            return datasetName + "_toppic_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_toppic_syn_ProteinMods.txt";
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
            return datasetName + "_toppic_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_toppic_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_toppic_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return TopPIC_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified TopPIC parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(TopPIC_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                const PeptideHitResultTypes resultType = PeptideHitResultTypes.TopPIC;
                var success = ReadKeyValuePairSearchEngineParamFile(TopPIC_SEARCH_ENGINE_NAME, searchEngineParamFileName, resultType, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                searchEngineParams.Enzyme = "no_enzyme";
                searchEngineParams.MinNumberTermini = 0;

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                searchEngineParams.PrecursorMassToleranceDa = MSGFPlusSynFileReader.DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM, resultType);
                searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM details</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;
            const string NOT_FOUND = "==SCORE_NOT_FOUND==";

            var columns = line.Split('\t');

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_ResultID, mColumnHeaders, 0);
                psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_Rank_PValue, mColumnHeaders, 1);

                var peptide = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_Peptide, mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = Convert.ToInt16(ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                var protein = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_Protein, mColumnHeaders);
                psm.AddProtein(protein);

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0);
                psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_DelM, mColumnHeaders);
                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                psm.MSGFSpecEValue = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_EValue, mColumnHeaders);

                // Store the remaining scores
                AddScore(psm, columns, DATA_COLUMN_Prsm_ID);
                AddScore(psm, columns, DATA_COLUMN_Spectrum_ID);
                AddScore(psm, columns, DATA_COLUMN_FragMethod);

                AddScore(psm, columns, DATA_COLUMN_MH);

                AddScore(psm, columns, DATA_COLUMN_Proteoform_ID);
                AddScore(psm, columns, DATA_COLUMN_Feature_Intensity);
                AddScore(psm, columns, DATA_COLUMN_Feature_Score);
                AddScore(psm, columns, DATA_COLUMN_ResidueStart);
                AddScore(psm, columns, DATA_COLUMN_ResidueEnd);

                AddScore(psm, columns, DATA_COLUMN_Unexpected_Mod_Count);
                AddScore(psm, columns, DATA_COLUMN_Peak_Count);
                AddScore(psm, columns, DATA_COLUMN_Matched_Peak_Count);
                AddScore(psm, columns, DATA_COLUMN_Matched_Fragment_Ion_Count);

                AddScore(psm, columns, DATA_COLUMN_PValue);
                AddScore(psm, columns, DATA_COLUMN_EValue);
                AddScore(psm, columns, DATA_COLUMN_QValue);

                var proteoformQValue = ReaderFactory.LookupColumnValue(columns, DATA_COLUMN_Proteoform_QValue, mColumnHeaders, NOT_FOUND);
                if (proteoformQValue != NOT_FOUND)
                {
                    psm.SetScore(DATA_COLUMN_Proteoform_QValue, proteoformQValue);
                }
                else
                {
                    AddScore(psm, columns, DATA_COLUMN_Proteoform_FDR);
                }

                AddScore(psm, columns, DATA_COLUMN_Variable_PTMs);

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the TopPIC data file: " + ex.Message);
                return false;
            }
        }
    }
}
