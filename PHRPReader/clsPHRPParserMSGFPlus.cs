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

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for MS-GF+
    /// </summary>
    public class clsPHRPParserMSGFPlus : clsPHRPParser
    {
        #region "Constants"

#pragma warning disable 1591
        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_FragMethod = "FragMethod";
        public const string DATA_COLUMN_SpecIndex = "SpecIndex";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_NTT = "NTT";
        public const string DATA_COLUMN_DeNovoScore = "DeNovoScore";
        public const string DATA_COLUMN_MSGFScore = "MSGFScore";

        public const string DATA_COLUMN_MSGFDB_SpecProb = "MSGFDB_SpecProb";           // MSGFDB
        public const string DATA_COLUMN_Rank_MSGFDB_SpecProb = "Rank_MSGFDB_SpecProb"; // MSGFDB

        public const string DATA_COLUMN_MSGFPlus_SpecEValue = "MSGFDB_SpecEValue";           // MS-GF+
        public const string DATA_COLUMN_Rank_MSGFPlus_SpecEValue = "Rank_MSGFDB_SpecEValue"; // MS-GF+

        public const string DATA_COLUMN_PValue = "PValue"; // MSGFDB
        public const string DATA_COLUMN_EValue = "EValue"; // MS-GF+

        public const string DATA_COLUMN_FDR = "FDR";        // MSGFDB; Only present if a Target/Decoy (TDA) search was used
        public const string DATA_COLUMN_PepFDR = "PepFDR";  // MSGFDB; Only valid if a Target/Decoy (TDA) search was used; if EFDR is present, will contain 1 for every row

        public const string DATA_COLUMN_QValue = "QValue";       // MS-GF+ reports QValue instead of FDR
        public const string DATA_COLUMN_PepQValue = "PepQValue"; // MS-GF+ reports pepQValue instead of PepFDR

        public const string DATA_COLUMN_EFDR = "EFDR";  // Only present if a Target/Decoy (TDA) search was not used

        public const string DATA_COLUMN_IMS_Scan = "IMS_Scan";
        public const string DATA_COLUMN_IMS_Drift_Time = "IMS_Drift_Time";

        public const string DATA_COLUMN_Isotope_Error = "IsotopeError"; // Only reported by MS-GF+

        // These suffixes were changed from_msgfdb to _msgfplus in November 2016
        public const string FILENAME_SUFFIX_SYN = "_msgfplus_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_msgfplus_fht.txt";

        // ReSharper disable once CommentTypo
        // Renamed from "MS-GFDB" to "MS-GF+" in November 2016
        private const string MSGFPLUS_SEARCH_ENGINE_NAME = "MS-GF+";

        public const string CHARGE_CARRIER_MASS_PARAM_NAME = "ChargeCarrierMass";

        /// <summary>
        /// These columns correspond to the Synopsis file created by clsMSGFPlusResultsProcessor
        /// </summary>
        public enum MSGFPlusSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            FragMethod = 2,
            SpecIndex = 3,
            Charge = 4,
            PrecursorMZ = 5,
            DelM = 6,                            // Precursor error, in Da; if the search used a tolerance less than 0.5 Da or less than 500 ppm, this value is computed from the DelMPPM value
            DelMPPM = 7,                         // Precursor error, in ppm; corrected for isotope selection errors
            MH = 8,                              // Theoretical monoisotopic peptide mass (computed by PHRP)
            Peptide = 9,                         // This is the sequence with prefix and suffix residues and also with modification symbols
            Protein = 10,                        // Protein Name (remove description)
            NTT = 11,                            // Number of tryptic terminii
            DeNovoScore = 12,
            MSGFScore = 13,
            SpecProb_EValue = 14,
            RankSpecProb = 15,                   // Rank 1 means lowest SpecEValue, 2 means next higher score, etc. (ties get the same rank)
            PValue_EValue = 16,
            FDR_QValue = 17,                     // Only present if searched using -tda 1
            PepFDR_PepQValue = 18,               // Only present if searched using -tda 1
            EFDR = 19,                           // Only present if did not search using -tda 1
            // ReSharper disable UnusedMember.Global
            IMSScan = 20,                        // Only present for MSGFDB_IMS results
            IMSDriftTime = 21,                   // Only present for MSGFDB_IMS results
            // ReSharper restore UnusedMember.Global
            IsotopeError = 22
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
        public clsPHRPParserMSGFPlus(string datasetName, string inputFilePath)
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
        public clsPHRPParserMSGFPlus(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MSGFPlus, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMSGFPlus(string datasetName, string inputFilePath, clsPHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MSGFPlus, startupOptions)
        {
        }

        /// <summary>
        /// Determines the precursor mass tolerance for either MS-GF+, MSPathFinder, or TopPIC
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <param name="resultType"></param>
        /// <returns>Precursor tolerance, in Da</returns>
        /// <remarks></remarks>
        public static double DeterminePrecursorMassTolerance(
            clsSearchEngineParameters searchEngineParams,
            out double tolerancePPM,
            clsPHRPReader.ePeptideHitResultType resultType)
        {
            var reExtraToleranceWithUnits = new Regex(@"([0-9.]+)([A-Za-z]+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);
            var reExtraToleranceNoUnits = new Regex(@"([0-9.]+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            string tolerance;
            double toleranceDa = 0;

            tolerancePPM = 0;

            if (resultType == clsPHRPReader.ePeptideHitResultType.TopPIC)
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
                if (!searchEngineParams.Parameters.TryGetValue("PrecursorMassTolerance", out tolerance))
                {
                    if (!searchEngineParams.Parameters.TryGetValue("PMTolerance", out tolerance))
                    {
                        return toleranceDa;
                    }
                }
            }

            // Parent mass tolerance
            // Might contain two values, separated by a comma
            var toleranceSplit = tolerance.Split(',');

            foreach (var item in toleranceSplit)
            {
                if (item.Trim().StartsWith("#"))
                    continue;

                Match reMatch;
                if (resultType == clsPHRPReader.ePeptideHitResultType.MSPathFinder ||
                    resultType == clsPHRPReader.ePeptideHitResultType.TopPIC)
                {
                    reMatch = reExtraToleranceNoUnits.Match(item);
                }
                else
                {
                    reMatch = reExtraToleranceWithUnits.Match(item);
                }

                if (!reMatch.Success)
                    continue;

                if (!double.TryParse(reMatch.Groups[1].Value, out var toleranceCurrent))
                    continue;

                if (resultType == clsPHRPReader.ePeptideHitResultType.MSPathFinder ||
                    resultType == clsPHRPReader.ePeptideHitResultType.TopPIC)
                {
                    // Units are always ppm
                    tolerancePPM = toleranceCurrent;
                    toleranceCurrent = clsPeptideMassCalculator.PPMToMass(toleranceCurrent, 2000);
                }
                else if (reMatch.Groups.Count > 1 && reMatch.Groups[2].Value.ToLower().Contains("ppm"))
                {
                    // Ppm
                    // Convert from PPM to dalton (assuming a mass of 2000 m/z)
                    tolerancePPM = toleranceCurrent;
                    toleranceCurrent = clsPeptideMassCalculator.PPMToMass(toleranceCurrent, 2000);
                }

                toleranceDa = Math.Max(toleranceDa, toleranceCurrent);
            }

            if (Math.Abs(tolerancePPM) < float.Epsilon & Math.Abs(toleranceDa) > float.Epsilon)
            {
                tolerancePPM = clsPeptideMassCalculator.MassToPPM(toleranceDa, 2000);
            }

            return toleranceDa;
        }

        /// <summary>
        /// Look for MS-GF+ parameter ChargeCarrierMass
        /// If defined, update chargeCarrierMass with the associated mass value and return True
        /// Otherwise return false
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <param name="chargeCarrierMass"></param>
        /// <returns></returns>
        /// <remarks>This function is used by clsPHRPMassErrorValidator in the Analysis Manager</remarks>
        public static bool GetCustomChargeCarrierMass(clsSearchEngineParameters searchEngineParams, out double chargeCarrierMass)
        {
            if (searchEngineParams.Parameters.TryGetValue(CHARGE_CARRIER_MASS_PARAM_NAME, out var value))
            {
                if (double.TryParse(value, out chargeCarrierMass))
                {
                    return true;
                }
            }

            chargeCarrierMass = clsPeptideMassCalculator.MASS_PROTON;
            return false;
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
        /// <remarks>This includes headers for synopsis files from both MSGFDB and MS-GF+</remarks>
        public static SortedDictionary<string, MSGFPlusSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            var headerColumns = new SortedDictionary<string, MSGFPlusSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {DATA_COLUMN_ResultID, MSGFPlusSynFileColumns.ResultID},
                {DATA_COLUMN_Scan, MSGFPlusSynFileColumns.Scan},
                {DATA_COLUMN_FragMethod, MSGFPlusSynFileColumns.FragMethod},
                {DATA_COLUMN_SpecIndex, MSGFPlusSynFileColumns.SpecIndex},
                {DATA_COLUMN_Charge, MSGFPlusSynFileColumns.Charge},
                {DATA_COLUMN_PrecursorMZ, MSGFPlusSynFileColumns.PrecursorMZ},
                {DATA_COLUMN_DelM, MSGFPlusSynFileColumns.DelM},
                {DATA_COLUMN_DelM_PPM, MSGFPlusSynFileColumns.DelMPPM},
                {DATA_COLUMN_MH, MSGFPlusSynFileColumns.MH},
                {DATA_COLUMN_Peptide, MSGFPlusSynFileColumns.Peptide},
                {DATA_COLUMN_Protein, MSGFPlusSynFileColumns.Protein},
                {DATA_COLUMN_NTT, MSGFPlusSynFileColumns.NTT},
                {DATA_COLUMN_DeNovoScore, MSGFPlusSynFileColumns.DeNovoScore},
                {DATA_COLUMN_MSGFScore, MSGFPlusSynFileColumns.MSGFScore},
                {DATA_COLUMN_MSGFDB_SpecProb, MSGFPlusSynFileColumns.SpecProb_EValue},
                {DATA_COLUMN_MSGFPlus_SpecEValue, MSGFPlusSynFileColumns.SpecProb_EValue},
                {DATA_COLUMN_Rank_MSGFDB_SpecProb, MSGFPlusSynFileColumns.RankSpecProb},
                {DATA_COLUMN_Rank_MSGFPlus_SpecEValue, MSGFPlusSynFileColumns.RankSpecProb},
                {DATA_COLUMN_PValue, MSGFPlusSynFileColumns.PValue_EValue},
                {DATA_COLUMN_EValue, MSGFPlusSynFileColumns.PValue_EValue},
                {DATA_COLUMN_FDR, MSGFPlusSynFileColumns.FDR_QValue},
                {DATA_COLUMN_QValue, MSGFPlusSynFileColumns.FDR_QValue},
                {DATA_COLUMN_PepFDR, MSGFPlusSynFileColumns.PepFDR_PepQValue},
                {DATA_COLUMN_PepQValue, MSGFPlusSynFileColumns.PepFDR_PepQValue},
                {DATA_COLUMN_EFDR, MSGFPlusSynFileColumns.EFDR},
                {DATA_COLUMN_Isotope_Error, MSGFPlusSynFileColumns.IsotopeError}
            };

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSGFPlusSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MSGFPlusSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
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
            return datasetName + FILENAME_SUFFIX_FHT;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_msgfplus_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msgfplus_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_msgfplus_syn_ProteinMods.txt";
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
            return datasetName + "_msgfplus_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_msgfplus_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
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
        /// Parses the specified MSGFDB (aka MS-GF+) parameter file
        /// </summary>
        /// <param name="searchEngineParamFilePath"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string searchEngineParamFilePath, out clsSearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new clsSearchEngineParameters(MSGFPLUS_SEARCH_ENGINE_NAME, mModInfo);

            var success = ReadSearchEngineParamFile(searchEngineParamFilePath, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFilePath, clsSearchEngineParameters searchEngineParams)
        {
            try
            {
                mPeptideMassCalculator.ResetAminoAcidMasses();

                var success = ReadKeyValuePairSearchEngineParamFile(MSGFPLUS_SEARCH_ENGINE_NAME, searchEngineParamFilePath, clsPHRPReader.ePeptideHitResultType.MSGFPlus, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                int value;

                // Determine the enzyme name
                if (searchEngineParams.Parameters.TryGetValue("enzymeid", out var settingValue))
                {
                    if (int.TryParse(settingValue, out value))
                    {
                        switch (value)
                        {
                            case 0:
                                searchEngineParams.Enzyme = "no_enzyme";
                                break;
                            case 1:
                                searchEngineParams.Enzyme = "trypsin";
                                break;
                            case 2:
                                searchEngineParams.Enzyme = "Chymotrypsin";
                                break;
                            case 3:
                                searchEngineParams.Enzyme = "Lys-C";
                                break;
                            case 4:
                                searchEngineParams.Enzyme = "Lys-N";
                                break;
                            case 5:
                                searchEngineParams.Enzyme = "Glu-C";
                                break;
                            case 6:
                                searchEngineParams.Enzyme = "Arg-C";
                                break;
                            case 7:
                                searchEngineParams.Enzyme = "Asp-N";
                                break;
                            case 8:
                                searchEngineParams.Enzyme = "alphaLP";
                                break;
                            case 9:
                                searchEngineParams.Enzyme = "no_enzyme_peptidomics";
                                break;
                            default:
                                searchEngineParams.Enzyme = "unknown_enzyme";
                                break;
                        }
                    }
                }

                // Determine the cleavage specificity
                if (searchEngineParams.Parameters.TryGetValue("nnet", out settingValue))
                {
                    // NNET means number of non-enzymatic terminii

                    if (int.TryParse(settingValue, out value))
                    {
                        switch (value)
                        {
                            case 0:
                                // Fully-tryptic
                                searchEngineParams.MinNumberTermini = 2;
                                break;
                            case 1:
                                // Partially-tryptic
                                searchEngineParams.MinNumberTermini = 1;
                                break;
                            default:
                                // No-enzyme search
                                searchEngineParams.MinNumberTermini = 0;
                                break;
                        }
                    }
                }
                else
                {
                    // MS-GF+ uses ntt instead of nnet; thus look for ntt

                    if (searchEngineParams.Parameters.TryGetValue("ntt", out settingValue))
                    {
                        // NTT means number of tolerable terminii

                        if (int.TryParse(settingValue, out value))
                        {
                            switch (value)
                            {
                                case 0:
                                    // No-enzyme search
                                    searchEngineParams.MinNumberTermini = 0;
                                    break;
                                case 1:
                                    // Partially-tryptic
                                    searchEngineParams.MinNumberTermini = 1;
                                    break;
                                default:
                                    // Fully-tryptic
                                    searchEngineParams.MinNumberTermini = 2;
                                    break;
                            }
                        }
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM, clsPHRPReader.ePeptideHitResultType.MSGFPlus);
                searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;

                // Look for Custom Amino Acid definitions
                if (searchEngineParams.Parameters.All(paramEntry => paramEntry.Key != clsMSGFPlusParamFileModExtractor.PARAM_TAG_CUSTOM_AA))
                {
                    // No custom amino acid entries
                    return true;
                }

                // Store the Custom Amino Acid info
                // Need to use a different parsing function to extract it
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
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out clsPSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;

            psm = new clsPSM();

            try
            {
                var columns = line.Split('\t');

                bool msgfPlusResults;
                if (clsPHRPReader.LookupColumnIndex(DATA_COLUMN_MSGFPlus_SpecEValue, mColumnHeaders) >= 0)
                {
                    msgfPlusResults = true;
                }
                else
                {
                    msgfPlusResults = false;
                }

                psm.DataLineText = line;
                psm.ScanNumber = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_ResultID, mColumnHeaders, 0);

                if (msgfPlusResults)
                {
                    psm.ScoreRank = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Rank_MSGFPlus_SpecEValue, mColumnHeaders, 1);
                }
                else
                {
                    psm.ScoreRank = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Rank_MSGFDB_SpecProb, mColumnHeaders, 1);
                }

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

                psm.CollisionMode = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_FragMethod, mColumnHeaders, "n/a");

                var precursorMZ = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0);
                psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                psm.MassErrorDa = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM, mColumnHeaders);
                psm.MassErrorPPM = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                if (msgfPlusResults)
                {
                    psm.MSGFSpecEValue = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_MSGFPlus_SpecEValue, mColumnHeaders);
                }
                else
                {
                    psm.MSGFSpecEValue = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_MSGFDB_SpecProb, mColumnHeaders);
                }

                if (psm.MSGFSpecEValue.Length > 13)
                {
                    // Attempt to shorten the SpecEValue value
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
                AddScore(psm, columns, DATA_COLUMN_DeNovoScore);

                AddScore(psm, columns, DATA_COLUMN_MSGFScore);

                if (msgfPlusResults)
                {
                    AddScore(psm, columns, DATA_COLUMN_MSGFPlus_SpecEValue);
                    AddScore(psm, columns, DATA_COLUMN_Rank_MSGFPlus_SpecEValue);
                    AddScore(psm, columns, DATA_COLUMN_EValue);
                    AddScore(psm, columns, DATA_COLUMN_QValue);
                    AddScore(psm, columns, DATA_COLUMN_PepQValue);
                    AddScore(psm, columns, DATA_COLUMN_Isotope_Error);

                    // Duplicate the score values to provide backwards compatibility
                    if (psm.TryGetScore(DATA_COLUMN_MSGFPlus_SpecEValue, out var value))
                        psm.SetScore(DATA_COLUMN_MSGFDB_SpecProb, value);
                    if (psm.TryGetScore(DATA_COLUMN_Rank_MSGFPlus_SpecEValue, out value))
                        psm.SetScore(DATA_COLUMN_Rank_MSGFDB_SpecProb, value);
                    if (psm.TryGetScore(DATA_COLUMN_QValue, out value))
                        psm.SetScore(DATA_COLUMN_FDR, value);
                    if (psm.TryGetScore(DATA_COLUMN_PepQValue, out value))
                        psm.SetScore(DATA_COLUMN_PepFDR, value);

                    var pValueStored = false;

                    if (psm.TryGetScore(DATA_COLUMN_EValue, out var eValueText))
                    {
                        if (psm.TryGetScore(DATA_COLUMN_MSGFPlus_SpecEValue, out var specEValueText))
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
                                            psm.SetScore(DATA_COLUMN_PValue, "0");
                                        }
                                        else
                                        {
                                            psm.SetScore(DATA_COLUMN_PValue, pValue.ToString("0.00000E-00"));
                                        }

                                        pValueStored = true;
                                    }
                                }
                            }
                        }

                        if (!pValueStored)
                        {
                            // Store E-value as P-value (these values are not identical, and will only be close for high-confidence results, i.e. results with FDR < 2%)
                            psm.SetScore(DATA_COLUMN_PValue, eValueText);
                        }
                    }
                }
                else
                {
                    AddScore(psm, columns, DATA_COLUMN_MSGFDB_SpecProb);
                    AddScore(psm, columns, DATA_COLUMN_Rank_MSGFDB_SpecProb);
                    AddScore(psm, columns, DATA_COLUMN_PValue);
                    AddScore(psm, columns, DATA_COLUMN_FDR);
                    AddScore(psm, columns, DATA_COLUMN_PepFDR);
                }

                AddScore(psm, columns, DATA_COLUMN_EFDR); // This column will not be present if a Target/Decoy (TDA) search was performed

                AddScore(psm, columns, DATA_COLUMN_IMS_Scan);

                AddScore(psm, columns, DATA_COLUMN_IMS_Drift_Time);

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MSGFDB data file: " + ex.Message);
                return false;
            }

        }

        private bool UpdateMassCalculatorMasses(string searchEngineParamFilePath)
        {
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MS-GF+");
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
        /// <param name="searchEngineParamFilePath"></param>
        /// <param name="modFileProcessor"></param>
        /// <param name="peptideMassCalculator"></param>
        /// <param name="errorMessage"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static bool UpdateMassCalculatorMasses(
            string searchEngineParamFilePath,
            clsMSGFPlusParamFileModExtractor modFileProcessor,
            clsPeptideMassCalculator peptideMassCalculator,
            out string errorMessage)
        {
            if (modFileProcessor == null)
            {
                throw new ObjectDisposedException("modFileProcessor is not initialized");
            }

            errorMessage = string.Empty;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                searchEngineParamFilePath,
                clsMSGFPlusParamFileModExtractor.ModSpecFormats.MSGFPlusAndMSPathFinder,
                out var modInfo);

            if (!success)
            {
                errorMessage = modFileProcessor.ErrorMessage;
                return false;
            }

            var customAminoAcidDefs = (from item in modInfo where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.CustomAA select item).ToList();
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
                    var elementalComposition = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormula);

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
