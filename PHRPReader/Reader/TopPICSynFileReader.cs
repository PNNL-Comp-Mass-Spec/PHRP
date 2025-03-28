﻿//*********************************************************************************************************
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
    /// PHRP SynFileReader for TopPIC
    /// </summary>
    public class TopPICSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: Da, FHT, Frag, Prot, prsm, psm, toppic

        /// <summary>
        /// TopPIC synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_toppic_syn.txt";

        /// <summary>
        /// TopPIC first hits file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_FHT = "_toppic_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string TopPIC_SEARCH_ENGINE_NAME = "TopPIC";

        /// <summary>
        /// Mapping from enum to synopsis file column name for TopPIC
        /// </summary>
        private static readonly Dictionary<TopPICSynFileColumns, string> mSynopsisFileColumn = new();

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
                { "ResultID", TopPICSynFileColumns.ResultID },
                { "Scan", TopPICSynFileColumns.Scan },
                { "Prsm_ID", TopPICSynFileColumns.Prsm_ID },
                { "Spectrum_ID", TopPICSynFileColumns.Spectrum_ID },
                { "FragMethod", TopPICSynFileColumns.FragMethod },
                { "Charge", TopPICSynFileColumns.Charge },
                { "PrecursorMZ", TopPICSynFileColumns.PrecursorMZ },
                { "DelM", TopPICSynFileColumns.DelM },
                { "DelM_PPM", TopPICSynFileColumns.DelMPPM },
                { "MH", TopPICSynFileColumns.MH },
                { "Peptide", TopPICSynFileColumns.Peptide },
                { "Proteoform_ID", TopPICSynFileColumns.Proteoform_ID },
                { "Feature_Intensity", TopPICSynFileColumns.Feature_Intensity },
                { "Feature_Score", TopPICSynFileColumns.Feature_Score },
                { "Feature_Apex_Time", TopPICSynFileColumns.Feature_Apex_Time },
                { "Protein_Count", TopPICSynFileColumns.Protein_Count },
                { "Protein", TopPICSynFileColumns.Protein },
                { "Protein_N-terminal_Form", TopPICSynFileColumns.Protein_Nterminal_Form },
                { "ResidueStart", TopPICSynFileColumns.ResidueStart },
                { "ResidueEnd", TopPICSynFileColumns.ResidueEnd },
                { "Unexpected_Mod_Count", TopPICSynFileColumns.Unexpected_Mod_Count },
                { "Peak_Count", TopPICSynFileColumns.Peak_Count },
                { "Matched_Peak_Count", TopPICSynFileColumns.Matched_Peak_Count },
                { "Matched_Fragment_Ion_Count", TopPICSynFileColumns.Matched_Fragment_Ion_Count },
                { "PValue", TopPICSynFileColumns.PValue },
                { "Rank_PValue", TopPICSynFileColumns.Rank_PValue },
                { "EValue", TopPICSynFileColumns.EValue },
                { "QValue", TopPICSynFileColumns.QValue },
            };

            if (!includeLegacyNames)
            {
                headerColumns.Add("Proteoform_QValue", TopPICSynFileColumns.Proteoform_QValue);
                headerColumns.Add("Variable_PTMs", TopPICSynFileColumns.VariablePTMs);
                return headerColumns;
            }

            headerColumns.Add("Proteoform_FDR", TopPICSynFileColumns.Proteoform_QValue);
            headerColumns.Add("Proteoform_QValue", TopPICSynFileColumns.Proteoform_QValue);
            headerColumns.Add("Variable_PTMs", TopPICSynFileColumns.VariablePTMs);

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum TopPICSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames">List of header names</param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<TopPICSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs(true);
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column">Column enum</param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(TopPICSynFileColumns column)
        {
            if (mSynopsisFileColumn.Count > 0)
            {
                return mSynopsisFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs(false))
            {
                mSynopsisFileColumn.Add(item.Value, item.Key);
            }

            return mSynopsisFileColumn[column];
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_FHT;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            // ReSharper disable once StringLiteralTypo
            return datasetName + "_toppic_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            // ReSharper disable once StringLiteralTypo
            return datasetName + "_toppic_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_toppic_syn_ProteinMods.txt";
        }

        /// <summary>
        /// Default Synopsis file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_SYN;
        }

        /// <summary>
        /// Default ResultToSeq map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(string datasetName)
        {
            return datasetName + "_toppic_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_toppic_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
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
        /// <param name="searchEngineParamFileName">TopPIC parameter file name</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
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
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM details</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;
            const string NOT_FOUND = "==SCORE_NOT_FOUND==";

            var columns = line.Split('\t');

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.ResultID), mColumnHeaders, 0);
                psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.Rank_PValue), mColumnHeaders, 1);

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.Charge), mColumnHeaders, 0);

                var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.Protein), mColumnHeaders);
                psm.AddProtein(protein);

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);
                psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.DelM), mColumnHeaders);
                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.DelMPPM), mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                psm.MSGFSpecEValue = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(TopPICSynFileColumns.EValue), mColumnHeaders);

                // Store the remaining scores
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Prsm_ID));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Spectrum_ID));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.FragMethod));

                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.MH));

                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Proteoform_ID));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Feature_Intensity));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Feature_Score));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.ResidueStart));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.ResidueEnd));

                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Unexpected_Mod_Count));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Peak_Count));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Matched_Peak_Count));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.Matched_Fragment_Ion_Count));

                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.PValue));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.EValue));
                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.QValue));

                var qValueColumn = GetColumnNameByID(TopPICSynFileColumns.Proteoform_QValue);

                var proteoformQValue = ReaderFactory.LookupColumnValue(columns, qValueColumn, mColumnHeaders, NOT_FOUND);

                if (proteoformQValue != NOT_FOUND)
                {
                    psm.SetScore(qValueColumn, proteoformQValue);
                }
                else
                {
                    AddScore(psm, columns, "Proteoform_FDR");
                }

                AddScore(psm, columns, GetColumnNameByID(TopPICSynFileColumns.VariablePTMs));

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
