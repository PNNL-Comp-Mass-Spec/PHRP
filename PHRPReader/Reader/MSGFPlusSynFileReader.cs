//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from MS-GF+ msgfplus_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for MS-GF+
    /// </summary>
    public class MSGFPlusSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: enzymeid, Da, Daltons, Frag, novo, msgfdb, peptidomics, prot, nnet, ntt, tda, tryptic, Za, Validator
        // Ignore Spelling: Arg, Chymotrypsin, Glu, Lys

        /// <summary>
        /// Legacy SpecProb column name for tool MSGFDB
        /// </summary>
        public const string MSGFDB_SpecProb = "MSGFDB_SpecProb";

        /// <summary>
        /// Legacy Rank SpecProb column name for tool MSGFDB
        /// </summary>
        private const string MSGFDB_RankSpecProb = "Rank_MSGFDB_SpecProb";

        /// <summary>
        /// Legacy PValue column name for tool MSGFDB
        /// </summary>
        private const string MSGFDB_PValue = "PValue";

        /// <summary>
        /// Legacy FDR column name for tool MSGFDB
        /// </summary>
        /// <remarks>
        /// Only present if a Target/Decoy (TDA) search was used
        /// </remarks>
        private const string MSGFDB_FDR = "FDR";

        /// <summary>
        /// Legacy PepFDR column name for tool MSGFDB
        /// </summary>
        /// <remarks>
        /// Only valid if a Target/Decoy (TDA) search was used; if EFDR is present, will contain 1 for every row
        /// </remarks>
        private const string MSGFDB_PepFDR = "PepFDR";

        /// <summary>
        /// MS-GF+ synopsis file suffix
        /// </summary>
        /// <remarks>
        /// Changed from _msgfdb_syn.txt to _msgfplus_syn.txt in November 2016
        /// </remarks>
        public const string FILENAME_SUFFIX_SYN = "_msgfplus_syn.txt";

        /// <summary>
        /// MS-GF+ first hits file suffix
        /// </summary>
        /// <remarks>
        /// Changed from _msgfdb_fht.txt to _msgfplus_fht.txt in November 2016
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_msgfplus_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string MSGFPLUS_SEARCH_ENGINE_NAME = "MS-GF+";

        /// <summary>
        /// Precursor mass tolerance parameter name
        /// </summary>
        public const string PRECURSOR_TOLERANCE_PARAM_NAME = "PrecursorMassTolerance";

        /// <summary>
        /// Precursor mass tolerance parameter name synonym
        /// </summary>
        public const string PRECURSOR_TOLERANCE_PARAM_NAME_SYNONYM = "PMTolerance";

        /// <summary>
        /// Charge carrier mass parameter name
        /// </summary>
        public const string CHARGE_CARRIER_MASS_PARAM_NAME = "ChargeCarrierMass";

        /// <summary>
        /// Mapping from enum to synopsis file column name for legacy tool MSGFDB
        /// </summary>
        private static readonly Dictionary<MSGFDBSynFileColumns, string> mMSGFDBColumns = new();

        /// <summary>
        /// Mapping from enum to synopsis file column name for MS-GF+
        /// </summary>
        private static readonly Dictionary<MSGFPlusSynFileColumns, string> mSynopsisFileColumn = new();

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
        public MSGFPlusSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MSGFPlusSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MSGFPlus, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MSGFPlusSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MSGFPlus, startupOptions)
        {
        }

        /// <summary>
        /// Determines the precursor mass tolerance for either MS-GF+, MSPathFinder, or TopPIC
        /// </summary>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <param name="resultType">PHRP result type</param>
        /// <returns>Precursor tolerance, in Da</returns>
        public static double DeterminePrecursorMassTolerance(
            SearchEngineParameters searchEngineParams,
            out double tolerancePPM,
            PeptideHitResultTypes resultType)
        {
            var extraToleranceMatcherWithUnits = new Regex("([0-9.]+)([A-Za-z]+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);
            var extraToleranceMatcherNoUnits = new Regex("([0-9.]+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            string tolerance;
            double toleranceDa = 0;

            tolerancePPM = 0;

            if (resultType == PeptideHitResultTypes.TopPIC)
            {
                // TopPIC
                if (!searchEngineParams.Parameters.TryGetValue("ErrorTolerance", out tolerance))
                {
                    return toleranceDa;
                }
            }
            else
            {
                // MS-GF+ or MSPathFinder
                if (!searchEngineParams.Parameters.TryGetValue(PRECURSOR_TOLERANCE_PARAM_NAME, out tolerance))
                {
                    if (!searchEngineParams.Parameters.TryGetValue(PRECURSOR_TOLERANCE_PARAM_NAME_SYNONYM, out tolerance))
                    {
                        return toleranceDa;
                    }
                }
            }

            // Parent mass tolerance
            // Might contain two values, separated by a comma
            foreach (var item in tolerance.Split(','))
            {
                if (item.Trim().StartsWith("#"))
                    continue;

                Match match;
                if (resultType is PeptideHitResultTypes.MSPathFinder or PeptideHitResultTypes.TopPIC)
                {
                    match = extraToleranceMatcherNoUnits.Match(item);
                }
                else
                {
                    match = extraToleranceMatcherWithUnits.Match(item);
                }

                if (!match.Success)
                    continue;

                if (!double.TryParse(match.Groups[1].Value, out var toleranceCurrent))
                    continue;

                if (resultType is PeptideHitResultTypes.MSPathFinder or PeptideHitResultTypes.TopPIC)
                {
                    // Units are always ppm
                    tolerancePPM = toleranceCurrent;
                    toleranceCurrent = PeptideMassCalculator.PPMToMass(toleranceCurrent, 2000);
                }
                else if (match.Groups.Count > 1 && match.Groups[2].Value.IndexOf("ppm", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    // Ppm
                    // Convert from PPM to Daltons (assuming a mass of 2000 m/z)
                    tolerancePPM = toleranceCurrent;
                    toleranceCurrent = PeptideMassCalculator.PPMToMass(toleranceCurrent, 2000);
                }

                toleranceDa = Math.Max(toleranceDa, toleranceCurrent);
            }

            if (Math.Abs(tolerancePPM) < float.Epsilon && Math.Abs(toleranceDa) > float.Epsilon)
            {
                tolerancePPM = PeptideMassCalculator.MassToPPM(toleranceDa, 2000);
            }

            return toleranceDa;
        }

        /// <summary>
        /// Look for MS-GF+ parameter ChargeCarrierMass
        /// If defined, update chargeCarrierMass with the associated mass value and return True
        /// Otherwise return false
        /// </summary>
        /// <remarks>This method is used by PHRPMassErrorValidator in the Analysis Manager</remarks>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <param name="chargeCarrierMass">Charge carrier mass</param>
        /// <returns>True if successful, false if an error</returns>
        public static bool GetCustomChargeCarrierMass(SearchEngineParameters searchEngineParams, out double chargeCarrierMass)
        {
            if (searchEngineParams.Parameters.TryGetValue(CHARGE_CARRIER_MASS_PARAM_NAME, out var value))
            {
                if (double.TryParse(value, out chargeCarrierMass))
                {
                    return true;
                }
            }

            chargeCarrierMass = PeptideMassCalculator.MASS_PROTON;
            return false;
        }

        /// <summary>
        /// Get the headerAddHeaderColumn names in the PHRP synopsis or first hits file for this tool
        /// </summary>
        /// <returns>List of header names, including legacy columns</returns>
        protected override List<string> GetColumnHeaderNames()
        {
            var headerNames = new List<string>();
            headerNames.AddRange(GetColumnHeaderNamesAndIDs(true, false).Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <remarks>This includes headers for synopsis files from both MSGFDB and MS-GF+</remarks>
        /// <param name="includeLegacyNames">When true, include legacy column names</param>
        /// <param name="includeExtraColumns">When true, include IMS_Scan and IMS_Drift_Time columns</param>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, MSGFPlusSynFileColumns> GetColumnHeaderNamesAndIDs(bool includeLegacyNames, bool includeExtraColumns)
        {
            var headerColumns = new SortedDictionary<string, MSGFPlusSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", MSGFPlusSynFileColumns.ResultID },
                { "Scan", MSGFPlusSynFileColumns.Scan },
                { "FragMethod", MSGFPlusSynFileColumns.FragMethod },
                { "SpecIndex", MSGFPlusSynFileColumns.SpecIndex },
                { "Charge", MSGFPlusSynFileColumns.Charge },
                { "PrecursorMZ", MSGFPlusSynFileColumns.PrecursorMZ },
                { "DelM", MSGFPlusSynFileColumns.DelM },
                { "DelM_PPM", MSGFPlusSynFileColumns.DelMPPM },
                { "MH", MSGFPlusSynFileColumns.MH },
                { "Peptide", MSGFPlusSynFileColumns.Peptide },
                { "Protein", MSGFPlusSynFileColumns.Protein },
                { "NTT", MSGFPlusSynFileColumns.NTT },
                { "DeNovoScore", MSGFPlusSynFileColumns.DeNovoScore },
                { "MSGFScore", MSGFPlusSynFileColumns.MSGFScore },
                { "MSGFDB_SpecEValue", MSGFPlusSynFileColumns.SpecEValue },
                { "Rank_MSGFDB_SpecEValue", MSGFPlusSynFileColumns.RankSpecEValue },
                { "EValue", MSGFPlusSynFileColumns.EValue },
                { "QValue", MSGFPlusSynFileColumns.QValue },
                { "PepQValue", MSGFPlusSynFileColumns.PepQValue },
                { "EFDR", MSGFPlusSynFileColumns.EFDR },
                { "IsotopeError", MSGFPlusSynFileColumns.IsotopeError }
            };

            if (includeExtraColumns)
            {
                headerColumns.Add("IMS_Scan", MSGFPlusSynFileColumns.IMSScan);
                headerColumns.Add("IMS_Drift_Time", MSGFPlusSynFileColumns.IMSDriftTime);
            }

            if (!includeLegacyNames)
                return headerColumns;

            headerColumns.Add(MSGFDB_SpecProb, MSGFPlusSynFileColumns.SpecEValue);
            headerColumns.Add(MSGFDB_RankSpecProb, MSGFPlusSynFileColumns.RankSpecEValue);
            headerColumns.Add(MSGFDB_PValue, MSGFPlusSynFileColumns.EValue);
            headerColumns.Add(MSGFDB_FDR, MSGFPlusSynFileColumns.QValue);
            headerColumns.Add(MSGFDB_PepFDR, MSGFPlusSynFileColumns.PepQValue);

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSGFPlusSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames">List of header names</param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MSGFPlusSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs(true, true);
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column">Column enum</param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MSGFPlusSynFileColumns column)
        {
            if (mSynopsisFileColumn.Count > 0)
            {
                return mSynopsisFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs(false, true))
            {
                mSynopsisFileColumn.Add(item.Value, item.Key);
            }

            return mSynopsisFileColumn[column];
        }

        /// <summary>
        /// Get the MSGFDB column name associated with the given column (for legacy tool MSGFDB, not MS-GF+)
        /// </summary>
        /// <param name="column">Column enum</param>
        /// <returns>Column name</returns>
        public static string GetMSGFDBColumnNameByID(MSGFDBSynFileColumns column)
        {
            if (mMSGFDBColumns.Count > 0)
            {
                return mMSGFDBColumns[column];
            }

            mMSGFDBColumns.Add(MSGFDBSynFileColumns.SpecProb, MSGFDB_SpecProb);
            mMSGFDBColumns.Add(MSGFDBSynFileColumns.RankSpecProb, MSGFDB_RankSpecProb);
            mMSGFDBColumns.Add(MSGFDBSynFileColumns.PValue, MSGFDB_PValue);
            mMSGFDBColumns.Add(MSGFDBSynFileColumns.FDR, MSGFDB_FDR);
            mMSGFDBColumns.Add(MSGFDBSynFileColumns.PepFDR, MSGFDB_PepFDR);

            return mMSGFDBColumns[column];
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
            return datasetName + "_msgfplus_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msgfplus_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_msgfplus_syn_ProteinMods.txt";
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
            return datasetName + "_msgfplus_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_msgfplus_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msgfplus_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MSGFPLUS_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MS-GF+ (previously MSGFDB) parameter file
        /// </summary>
        /// <param name="searchEngineParamFilePath">MS-GF+ parameter file name</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFilePath, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MSGFPLUS_SEARCH_ENGINE_NAME, mModInfo);

            var success = ReadSearchEngineParamFile(searchEngineParamFilePath, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFilePath, SearchEngineParameters searchEngineParams)
        {
            try
            {
                mPeptideMassCalculator.ResetAminoAcidMasses();

                var success = ReadKeyValuePairSearchEngineParamFile(MSGFPLUS_SEARCH_ENGINE_NAME, searchEngineParamFilePath, PeptideHitResultTypes.MSGFPlus, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                // Determine the enzyme name
                if (searchEngineParams.Parameters.TryGetValue("enzymeid", out var settingValue) &&
                    int.TryParse(settingValue, out var value))
                {
                    searchEngineParams.Enzyme = value switch
                    {
                        0 => "no_enzyme",
                        1 => "trypsin",
                        2 => "Chymotrypsin",
                        3 => "Lys-C",
                        4 => "Lys-N",
                        5 => "Glu-C",
                        6 => "Arg-C",
                        7 => "Asp-N",
                        8 => "alphaLP",
                        9 => "no_enzyme_peptidomics",
                        _ => "unknown_enzyme",
                    };
                }

                // Determine the cleavage specificity
                if (searchEngineParams.Parameters.TryGetValue("nnet", out settingValue))
                {
                    // NNET means number of non-enzymatic termini

                    if (int.TryParse(settingValue, out value))
                    {
                        searchEngineParams.MinNumberTermini = value switch
                        {
                            0 => 2, // Fully-tryptic
                            1 => 1, // Partially-tryptic
                            _ => 0, // No-enzyme search
                        };
                    }
                }
                else
                {
                    // MS-GF+ uses ntt instead of nnet; thus look for ntt

                    if (searchEngineParams.Parameters.TryGetValue("ntt", out settingValue))
                    {
                        // NTT means number of tolerable termini

                        if (int.TryParse(settingValue, out value))
                        {
                            searchEngineParams.MinNumberTermini = value switch
                            {
                                0 => 0, // No-enzyme search
                                1 => 1, // Partially-tryptic
                                _ => 2, // Fully-tryptic
                            };
                        }
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM, PeptideHitResultTypes.MSGFPlus);
                searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;

                // Look for Custom Amino Acid definitions
                if (searchEngineParams.Parameters.All(paramEntry => paramEntry.Key != MSGFPlusParamFileModExtractor.PARAM_TAG_CUSTOM_AA))
                {
                    // No custom amino acid entries
                    return true;
                }

                // Store the Custom Amino Acid info
                // Need to use a different parsing method to extract it
                success = UpdateMassCalculatorMasses(searchEngineParamFilePath);

                // Look for a custom charge carrier mass
                if (GetCustomChargeCarrierMass(searchEngineParams, out var customChargeCarrierMass))
                {
                    ShowMessage(string.Format("Using a charge carrier mass of {0:F3} Da", customChargeCarrierMass));
                    mPeptideMassCalculator.ChargeCarrierMass = customChargeCarrierMass;
                }

                return success;
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

            psm = new PSM();

            try
            {
                var columns = line.Split('\t');

                var msgfPlusResults = ReaderFactory.LookupColumnIndex(GetColumnNameByID(MSGFPlusSynFileColumns.SpecEValue), mColumnHeaders) >= 0;

                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.ResultID), mColumnHeaders, 0);

                if (msgfPlusResults)
                {
                    psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.RankSpecEValue), mColumnHeaders, 1);
                }
                else
                {
                    psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.RankSpecProb), mColumnHeaders, 1);
                }

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.Charge), mColumnHeaders, 0);

                var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.Protein), mColumnHeaders);
                psm.AddProtein(protein);

                psm.CollisionMode = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.FragMethod), mColumnHeaders, "n/a");

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);
                psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.DelM), mColumnHeaders);
                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.DelMPPM), mColumnHeaders);

                if (msgfPlusResults)
                {
                    psm.MSGFSpecEValue = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSGFPlusSynFileColumns.SpecEValue), mColumnHeaders);
                }
                else
                {
                    psm.MSGFSpecEValue = ReaderFactory.LookupColumnValue(columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.SpecProb), mColumnHeaders);
                }

                if (psm.MSGFSpecEValue.Length > 13)
                {
                    // Attempt to shorten the SpecEValue
                    if (double.TryParse(psm.MSGFSpecEValue, out var specEValue))
                    {
                        psm.MSGFSpecEValue = specEValue.ToString("0.0000000E-00");
                    }
                }

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                // Store the remaining scores
                AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.DeNovoScore));

                AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.MSGFScore));

                if (msgfPlusResults)
                {
                    AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.SpecEValue));
                    AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.RankSpecEValue));
                    AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.EValue));
                    AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.QValue));
                    AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.PepQValue));
                    AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.IsotopeError));

                    // Duplicate the score values to provide backwards compatibility
                    if (psm.TryGetScore(GetColumnNameByID(MSGFPlusSynFileColumns.SpecEValue), out var value))
                        psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.SpecProb), value);

                    if (psm.TryGetScore(GetColumnNameByID(MSGFPlusSynFileColumns.RankSpecEValue), out value))
                        psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.RankSpecProb), value);

                    if (psm.TryGetScore(GetColumnNameByID(MSGFPlusSynFileColumns.QValue), out value))
                        psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.FDR), value);

                    if (psm.TryGetScore(GetColumnNameByID(MSGFPlusSynFileColumns.PepQValue), out value))
                        psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PepFDR), value);

                    var pValueStored = false;

                    if (psm.TryGetScore(GetColumnNameByID(MSGFPlusSynFileColumns.EValue), out var eValueText))
                    {
                        if (psm.TryGetScore(GetColumnNameByID(MSGFPlusSynFileColumns.SpecEValue), out var specEValueText))
                        {
                            // Compute PValue using EValue and SpecEValue
                            if (double.TryParse(eValueText, out var eValue))
                            {
                                if (double.TryParse(specEValueText, out var specEValue))
                                {
                                    if (specEValue > 0)
                                    {
                                        var n = eValue / specEValue;
                                        var pValue = 1 - Math.Pow(1 - specEValue, n);

                                        if (Math.Abs(pValue) <= double.Epsilon)
                                        {
                                            psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PValue), "0");
                                        }
                                        else
                                        {
                                            psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PValue), pValue.ToString("0.00000E-00"));
                                        }

                                        pValueStored = true;
                                    }
                                }
                            }
                        }

                        if (!pValueStored)
                        {
                            // Store E-value as P-value (these values are not identical, and will only be close for high-confidence results, i.e. results with FDR < 2%)
                            psm.SetScore(GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PValue), eValueText);
                        }
                    }
                }
                else
                {
                    AddScore(psm, columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.SpecProb));
                    AddScore(psm, columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.RankSpecProb));
                    AddScore(psm, columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PValue));
                    AddScore(psm, columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.FDR));
                    AddScore(psm, columns, GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PepFDR));
                }

                AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.EFDR)); // This column will not be present if a Target/Decoy (TDA) search was performed

                AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.IMSScan));

                AddScore(psm, columns, GetColumnNameByID(MSGFPlusSynFileColumns.IMSDriftTime));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MS-GF+ data file: " + ex.Message);
                return false;
            }
        }

        private bool UpdateMassCalculatorMasses(string searchEngineParamFilePath)
        {
            var modFileProcessor = new MSGFPlusParamFileModExtractor("MS-GF+");
            RegisterEvents(modFileProcessor);

            var success = UpdateMassCalculatorMasses(searchEngineParamFilePath, modFileProcessor, mPeptideMassCalculator, out var localErrorMsg);

            if (!string.IsNullOrWhiteSpace(localErrorMsg) && string.IsNullOrWhiteSpace(mErrorMessage))
            {
                ReportError(localErrorMsg);
            }

            return success;
        }

        /// <summary>
        /// Look for custom amino acid definitions in the MS-GF+ parameter file
        /// If any are found, update the amino acid mass values in the PeptideMassCalculator instance
        /// </summary>
        /// <param name="searchEngineParamFilePath">MS-GF+ parameter file path</param>
        /// <param name="modFileProcessor">Modification parser</param>
        /// <param name="peptideMassCalculator">Peptide mass calculator</param>
        /// <param name="errorMessage">Output: error message</param>
        /// <returns>True if successful, false if an error</returns>
        public static bool UpdateMassCalculatorMasses(
            string searchEngineParamFilePath,
            MSGFPlusParamFileModExtractor modFileProcessor,
            PeptideMassCalculator peptideMassCalculator,
            out string errorMessage)
        {
            if (modFileProcessor == null)
            {
                throw new ObjectDisposedException("modFileProcessor is not initialized");
            }

            errorMessage = string.Empty;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                searchEngineParamFilePath,
                MSGFPlusParamFileModExtractor.ModSpecFormats.MSGFPlusAndMSPathFinder,
                out var modDef);

            if (!success)
            {
                errorMessage = modFileProcessor.ErrorMessage;
                return false;
            }

            var customAminoAcidDefs = (from item in modDef where item.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.CustomAA select item).ToList();

            if (customAminoAcidDefs.Count == 0)
            {
                // There are no custom amino acids
                return true;
            }

            foreach (var customAADef in customAminoAcidDefs)
            {
                var aminoAcidSymbol = customAADef.Residues[0];
                var empiricalFormula = customAADef.ModMass;
                var aminoAcidMass = customAADef.ModMassVal;

                try
                {
                    var elementalComposition = PeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormula);

                    peptideMassCalculator.SetAminoAcidMass(aminoAcidSymbol, aminoAcidMass);
                    peptideMassCalculator.SetAminoAcidAtomCounts(aminoAcidSymbol, elementalComposition);
                }
                catch (Exception ex)
                {
                    errorMessage = ex.Message;
                    return false;
                }
            }

            return true;
        }
    }
}
