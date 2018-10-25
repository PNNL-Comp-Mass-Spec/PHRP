//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from MSGFDB msgfdb_syn.txt files
//
//*********************************************************************************************************
using System;
using System.Linq;
using System.Text.RegularExpressions;

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for MSGF+
    /// </summary>
    public class clsPHRPParserMSGFDB : clsPHRPParser
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

        public const string DATA_COLUMN_MSGFPlus_SpecEValue = "MSGFDB_SpecEValue";           // MSGF+
        public const string DATA_COLUMN_Rank_MSGFPlus_SpecEValue = "Rank_MSGFDB_SpecEValue"; // MSGF+

        public const string DATA_COLUMN_PValue = "PValue"; // MSGFDB
        public const string DATA_COLUMN_EValue = "EValue"; // MSGF+

        public const string DATA_COLUMN_FDR = "FDR";        // MSGFDB; Only present if a Target/Decoy (TDA) search was used
        public const string DATA_COLUMN_PepFDR = "PepFDR";  // MSGFDB; Only valid if a Target/Decoy (TDA) search was used; if EFDR is present, will contain 1 for every row

        public const string DATA_COLUMN_QValue = "QValue";       // MSGF+ reports QValue instead of FDR
        public const string DATA_COLUMN_PepQValue = "PepQValue"; // MSGF+ reports pepQValue instead of PepFDR

        public const string DATA_COLUMN_EFDR = "EFDR";  // Only present if a Target/Decoy (TDA) search was not used

        public const string DATA_COLUMN_IMS_Scan = "IMS_Scan";
        public const string DATA_COLUMN_IMS_Drift_Time = "IMS_Drift_Time";

        public const string DATA_COLUMN_Isotope_Error = "IsotopeError"; // Only reported by MSGF+

        // These suffixes were changed from_msgfdb to _msgfplus in November 2016
        public const string FILENAME_SUFFIX_SYN = "_msgfplus_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_msgfplus_fht.txt";

        // ReSharper disable once CommentTypo
        // Renamed from "MS-GFDB" to "MS-GF+" in November 2016
        private const string MSGFDB_SEARCH_ENGINE_NAME = "MS-GF+";

        public const string CHARGE_CARRIER_MASS_PARAM_NAME = "ChargeCarrierMass";
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
        public clsPHRPParserMSGFDB(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, then load the ModSummary file and SeqInfo files</param>
        /// <remarks></remarks>
        public clsPHRPParserMSGFDB(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMSGFDB(string datasetName, string inputFilePath, clsPHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB, startupOptions)
        {
        }

        /// <summary>
        /// Define column header names
        /// </summary>
        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_ResultID);
            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_FragMethod);
            AddHeaderColumn(DATA_COLUMN_SpecIndex);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_PrecursorMZ);
            AddHeaderColumn(DATA_COLUMN_DelM);
            AddHeaderColumn(DATA_COLUMN_DelM_PPM);
            AddHeaderColumn(DATA_COLUMN_MH);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_NTT);
            AddHeaderColumn(DATA_COLUMN_DeNovoScore);
            AddHeaderColumn(DATA_COLUMN_MSGFScore);
            AddHeaderColumn(DATA_COLUMN_MSGFDB_SpecProb);
            AddHeaderColumn(DATA_COLUMN_Rank_MSGFDB_SpecProb);
            AddHeaderColumn(DATA_COLUMN_PValue);
            AddHeaderColumn(DATA_COLUMN_FDR);
            AddHeaderColumn(DATA_COLUMN_EFDR);
            AddHeaderColumn(DATA_COLUMN_PepFDR);

            // Add the MSGF+ columns
            AddHeaderColumn(DATA_COLUMN_MSGFPlus_SpecEValue);
            AddHeaderColumn(DATA_COLUMN_Rank_MSGFPlus_SpecEValue);
            AddHeaderColumn(DATA_COLUMN_EValue);

            AddHeaderColumn(DATA_COLUMN_QValue);
            AddHeaderColumn(DATA_COLUMN_PepQValue);

            AddHeaderColumn(DATA_COLUMN_IMS_Scan);
            AddHeaderColumn(DATA_COLUMN_IMS_Drift_Time);
            AddHeaderColumn(DATA_COLUMN_Isotope_Error);
        }

        /// <summary>
        /// Determines the precursor mass tolerance for either MSGF+ or MSPathFinder
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

            double toleranceDa = 0;

            tolerancePPM = 0;

            if (!searchEngineParams.Parameters.TryGetValue("PMTolerance", out var tolerance))
            {
                return toleranceDa;
            }

            // Parent mass tolerance
            // Might contain two values, separated by a comma
            var toleranceSplit = tolerance.Split(',');

            foreach (var item in toleranceSplit)
            {
                if (item.Trim().StartsWith("#"))
                    continue;

                Match reMatch;
                if (resultType == clsPHRPReader.ePeptideHitResultType.MSPathFinder)
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

                if (resultType == clsPHRPReader.ePeptideHitResultType.MSPathFinder)
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
        /// Look for MSGF+ parameter ChargeCarrierMass
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
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_FHT;
        }

        /// <summary>
        /// Defeault ModSummary file name for the given dataset
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
            return MSGFDB_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MSGFDB (aka MS-GF+) parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out clsSearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new clsSearchEngineParameters(MSGFDB_SEARCH_ENGINE_NAME, mModInfo);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, clsSearchEngineParameters searchEngineParams)
        {
            try
            {
                mPeptideMassCalculator.ResetAminoAcidMasses();

                var success = ReadKeyValuePairSearchEngineParamFile(MSGFDB_SEARCH_ENGINE_NAME, searchEngineParamFileName, clsPHRPReader.ePeptideHitResultType.MSGFDB, searchEngineParams);

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
                    // MSGF+ uses ntt instead of nnet; thus look for ntt

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
                searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM, clsPHRPReader.ePeptideHitResultType.MSGFDB);
                searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;

                // Look for Custom Amino Acid definitions
                if (searchEngineParams.Parameters.All(paramEntry => paramEntry.Key != clsMSGFPlusParamFileModExtractor.PARAM_TAG_CUSTOM_AA))
                {
                    // No custom amino acid entries
                    return true;
                }

                // Store the Custom Amino Acid info
                // Need to use a different parsing function to extract it
                success = UpdateMassCalculatorMasses(searchEngineParamFileName);

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
        /// <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out clsPSM psm, bool fastReadMode)
        {
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
                psm.ScanNumber = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, -100);
                if (psm.ScanNumber == -100)
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

        private bool UpdateMassCalculatorMasses(string searchEngineParamFileName)
        {
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MSGF+");
            RegisterEvents(modFileProcessor);

            var success = UpdateMassCalculatorMasses(searchEngineParamFileName, modFileProcessor, mPeptideMassCalculator, out var localErrorMsg);

            if (!string.IsNullOrWhiteSpace(localErrorMsg) && string.IsNullOrWhiteSpace(mErrorMessage))
            {
                ReportError(localErrorMsg);
            }

            return success;
        }

        /// <summary>
        /// Look for custom amino acid definitions in the MSGF+ parameter file
        /// If any are found, update the amino acid mass values in the PeptideMassCalculator instance
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="modFileProcessor"></param>
        /// <param name="peptideMassCalculator"></param>
        /// <param name="errorMessage"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static bool UpdateMassCalculatorMasses(
            string searchEngineParamFileName,
            clsMSGFPlusParamFileModExtractor modFileProcessor,
            clsPeptideMassCalculator peptideMassCalculator,
            out string errorMessage)
        {
            if (modFileProcessor == null)
            {
                throw new ObjectDisposedException("modFileProcessor is not initialized");
            }

            errorMessage = string.Empty;

            // Note that this call will initialize modInfo
            var success = modFileProcessor.ExtractModInfoFromParamFile(searchEngineParamFileName, out var modInfo);
            if (!success)
            {
                errorMessage = modFileProcessor.ErrorMessage;
                return false;
            }

            var customAminoAcidDefs = (from item in modInfo where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA select item).ToList();
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
