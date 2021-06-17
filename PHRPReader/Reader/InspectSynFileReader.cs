//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from InSpecT _inspect_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for InSpecT
    /// </summary>
    public class InspectSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: chymotrypsin, Da, FilePos, PepToProt

        /// <summary>
        /// InSpecT synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_inspect_syn.txt";

        /// <summary>
        /// InSpecT first hits file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_FHT = "_inspect_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string INS_SEARCH_ENGINE_NAME = "InSpecT";

        /// <summary>
        /// Mapping from enum to synopsis file column name for InSpecT
        /// </summary>
        private static readonly Dictionary<InspectSynFileColumns, string> mSynopsisFileColumn = new();

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
        public InspectSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public InspectSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.Inspect, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public InspectSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.Inspect, startupOptions)
        {
        }

        /// <summary>
        /// Determines the precursor mass tolerance
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <returns>Precursor tolerance, in Da</returns>
        private double DeterminePrecursorMassTolerance(SearchEngineParameters searchEngineParams, out double tolerancePPM)
        {
            double toleranceDa = 0;
            double tolerance = 0;

            tolerancePPM = 0;

            if (searchEngineParams.Parameters.TryGetValue("ParentPPM", out var toleranceText))
            {
                // Parent mass tolerance, in ppm
                double.TryParse(toleranceText, out tolerancePPM);
                // Convert from PPM to Dalton (assuming a mass of 2000 m/z)
                tolerance = PeptideMassCalculator.PPMToMass(tolerancePPM, 2000);
            }

            if (searchEngineParams.Parameters.TryGetValue("PMTolerance", out toleranceText))
            {
                // Parent mass tolerance, in Dalton
                double.TryParse(toleranceText, out toleranceDa);

                // Convert from Dalton to PPM (assuming a mass of 2000 m/z)
                tolerancePPM = PeptideMassCalculator.MassToPPM(toleranceDa, 2000);
            }

            return Math.Max(tolerance, toleranceDa);
        }

        /// <summary>
        /// Get the header names in the PHRP synopsis or first hits file for this tool
        /// </summary>
        /// <returns>List of header names</returns>
        protected override List<string> GetColumnHeaderNames()
        {
            var headerNames = new List<string>();
            headerNames.AddRange(GetColumnHeaderNamesAndIDs().Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, InspectSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", InspectSynFileColumns.ResultID },
                { "Scan", InspectSynFileColumns.Scan },
                { "Peptide", InspectSynFileColumns.Peptide },
                { "Protein", InspectSynFileColumns.Protein },
                { "Charge", InspectSynFileColumns.Charge },
                { "MQScore", InspectSynFileColumns.MQScore },
                { "Length", InspectSynFileColumns.Length },
                { "TotalPRMScore", InspectSynFileColumns.TotalPRMScore },
                { "MedianPRMScore", InspectSynFileColumns.MedianPRMScore },
                { "FractionY", InspectSynFileColumns.FractionY },
                { "FractionB", InspectSynFileColumns.FractionB },
                { "Intensity", InspectSynFileColumns.Intensity },
                { "NTT", InspectSynFileColumns.NTT },
                { "PValue", InspectSynFileColumns.PValue },
                { "FScore", InspectSynFileColumns.FScore },
                { "DeltaScore", InspectSynFileColumns.DeltaScore },
                { "DeltaScoreOther", InspectSynFileColumns.DeltaScoreOther },
                { "DeltaNormMQScore", InspectSynFileColumns.DeltaNormMQScore },
                { "DeltaNormTotalPRMScore", InspectSynFileColumns.DeltaNormTotalPRMScore },
                { "RankTotalPRMScore", InspectSynFileColumns.RankTotalPRMScore },
                { "RankFScore", InspectSynFileColumns.RankFScore },
                { "MH", InspectSynFileColumns.MH },
                { "RecordNumber", InspectSynFileColumns.RecordNumber },
                { "DBFilePos", InspectSynFileColumns.DBFilePos },
                { "SpecFilePos", InspectSynFileColumns.SpecFilePos },
                { "PrecursorMZ", InspectSynFileColumns.PrecursorMZ },
                { "PrecursorError", InspectSynFileColumns.PrecursorError },
                { "DelM_PPM", InspectSynFileColumns.DelM_PPM }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum InspectSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<InspectSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(InspectSynFileColumns column)
        {
            if (mSynopsisFileColumn.Count > 0)
            {
                return mSynopsisFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mSynopsisFileColumn.Add(item.Value, item.Key);
            }

            return mSynopsisFileColumn[column];
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
            return datasetName + "_inspect_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_inspect_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_inspect_syn_ProteinMods.txt";
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
            return datasetName + "_inspect_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_inspect_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_inspect_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return INS_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified InSpecT parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(INS_SEARCH_ENGINE_NAME, mModInfo);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            var success = false;

            try
            {
                success = ReadKeyValuePairSearchEngineParamFile(INS_SEARCH_ENGINE_NAME, searchEngineParamFileName, PeptideHitResultTypes.Inspect, searchEngineParams);

                if (success)
                {
                    // Determine the enzyme name
                    if (searchEngineParams.Parameters.TryGetValue("protease", out var settingValue))
                    {
                        switch (settingValue.ToLower())
                        {
                            case "trypsin":
                                searchEngineParams.Enzyme = "trypsin";
                                break;
                            case "none":
                                searchEngineParams.Enzyme = "no_enzyme";
                                break;
                            case "chymotrypsin":
                                searchEngineParams.Enzyme = "chymotrypsin";
                                break;
                            default:
                                if (!string.IsNullOrEmpty(settingValue))
                                {
                                    searchEngineParams.Enzyme = settingValue;
                                }
                                break;
                        }
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM);
                    searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;
                }
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
            }

            return success;
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

            var columns = line.Split('\t');

            var success = false;

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.ResultID), mColumnHeaders, 0);
                    psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.RankTotalPRMScore), mColumnHeaders, 0);

                    var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.Peptide), mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.Charge), mColumnHeaders, 0);

                    var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.Protein), mColumnHeaders);
                    psm.AddProtein(protein);

                    var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                    psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.PrecursorError), mColumnHeaders);
                    psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(InspectSynFileColumns.DelM_PPM), mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.MQScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.TotalPRMScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.MedianPRMScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.PValue));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.FScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.DeltaScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.DeltaScoreOther));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.DeltaNormMQScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.DeltaNormTotalPRMScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.RankTotalPRMScore));
                    AddScore(psm, columns, GetColumnNameByID(InspectSynFileColumns.RankFScore));
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the InSpecT data file: " + ex.Message);
            }

            return success;
        }
    }
}
