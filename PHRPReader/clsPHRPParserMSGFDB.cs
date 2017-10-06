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
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;

namespace PHRPReader
{
    public class clsPHRPParserMSGFDB : clsPHRPParser
    {
        #region "Constants"
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

        // Renamed from "MS-GFDB" to "MS-GF+" in November 2016
        private const string MSGFDB_SEARCH_ENGINE_NAME = "MS-GF+";

        public const string CHARGE_CARRIER_MASS_PARAM_NAME = "ChargeCarrierMass";

        #endregion

        #region "Properties"

        public override string PHRPFirstHitsFileName
        {
            get { return GetPHRPFirstHitsFileName(mDatasetName); }
        }

        public override string PHRPModSummaryFileName
        {
            get { return GetPHRPModSummaryFileName(mDatasetName); }
        }

        public override string PHRPPepToProteinMapFileName
        {
            get { return GetPHRPPepToProteinMapFileName(mDatasetName); }
        }

        public override string PHRPProteinModsFileName
        {
            get { return GetPHRPProteinModsFileName(mDatasetName); }
        }

        public override string PHRPSynopsisFileName
        {
            get { return GetPHRPSynopsisFileName(mDatasetName); }
        }

        public override string PHRPResultToSeqMapFileName
        {
            get { return GetPHRPResultToSeqMapFileName(mDatasetName); }
        }

        public override string PHRPSeqInfoFileName
        {
            get { return GetPHRPSeqInfoFileName(mDatasetName); }
        }

        public override string PHRPSeqToProteinMapFileName
        {
            get { return GetPHRPSeqToProteinMapFileName(mDatasetName); }
        }

        public override string SearchEngineName
        {
            get { return GetSearchEngineName(); }
        }

        #endregion

        /// <summary>
        /// Constructor; assumes blnLoadModsAndSeqInfo=True
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <remarks></remarks>
        public clsPHRPParserMSGFDB(string strDatasetName, string strInputFilePath)
            : this(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="blnLoadModsAndSeqInfo">If True, then load the ModSummary file and SeqInfo files</param>
        /// <remarks></remarks>
        public clsPHRPParserMSGFDB(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMSGFDB(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB, startupOptions)
        {
        }

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
        /// <param name="objSearchEngineParams"></param>
        /// <param name="dblTolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <param name="resultType"></param>
        /// <returns>Precursor tolerance, in Da</returns>
        /// <remarks></remarks>
        public static double DeterminePrecursorMassTolerance(
            clsSearchEngineParameters objSearchEngineParams,
            out double dblTolerancePPM,
            clsPHRPReader.ePeptideHitResultType resultType)
        {
            var strTolerance = string.Empty;
            string[] strToleranceSplit = null;

            Match reMatch = default(Match);
            var reExtraToleranceWithUnits = new Regex(@"([0-9.]+)([A-Za-z]+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);
            var reExtraToleranceNoUnits = new Regex(@"([0-9.]+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            double dblToleranceDa = 0;

            dblTolerancePPM = 0;

            if (!objSearchEngineParams.Parameters.TryGetValue("PMTolerance", out strTolerance))
            {
                return dblToleranceDa;
            }

            // Parent mass tolerance
            // Might contain two values, separated by a comma
            strToleranceSplit = strTolerance.Split(',');

            if (strToleranceSplit == null)
            {
                return dblToleranceDa;
            }

            foreach (var strItem in strToleranceSplit)
            {
                if (strItem.Trim().StartsWith("#"))
                    continue;

                if (resultType == clsPHRPReader.ePeptideHitResultType.MSPathFinder)
                {
                    reMatch = reExtraToleranceNoUnits.Match(strItem);
                }
                else
                {
                    reMatch = reExtraToleranceWithUnits.Match(strItem);
                }

                if (!reMatch.Success)
                    continue;

                double dblToleranceCurrent = 0;

                if (!double.TryParse(reMatch.Groups[1].Value, out dblToleranceCurrent))
                    continue;

                if (resultType == clsPHRPReader.ePeptideHitResultType.MSPathFinder)
                {
                    // Units are always ppm
                    dblTolerancePPM = dblToleranceCurrent;
                    dblToleranceCurrent = clsPeptideMassCalculator.PPMToMass(dblToleranceCurrent, 2000);
                }
                else if (reMatch.Groups.Count > 1 && reMatch.Groups[2].Value.ToLower().Contains("ppm"))
                {
                    // Ppm
                    // Convert from PPM to dalton (assuming a mass of 2000 m/z)
                    dblTolerancePPM = dblToleranceCurrent;
                    dblToleranceCurrent = clsPeptideMassCalculator.PPMToMass(dblToleranceCurrent, 2000);
                }

                dblToleranceDa = Math.Max(dblToleranceDa, dblToleranceCurrent);
            }

            if (Math.Abs(dblTolerancePPM) < float.Epsilon & Math.Abs(dblToleranceDa) > float.Epsilon)
            {
                dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000);
            }

            return dblToleranceDa;
        }

        /// <summary>
        /// Look for MSGF+ parameter ChargeCarrierMass
        /// If defined, update chargeCarrierMass with the associated mass value and return True
        /// Otherwise return false
        /// </summary>
        /// <param name="objSearchEngineParams"></param>
        /// <param name="chargeCarrierMass"></param>
        /// <returns></returns>
        /// <remarks>This function is used by clsPHRPMassErrorValidator in the Analysis Manager</remarks>
        public static bool GetCustomChargeCarrierMass(clsSearchEngineParameters objSearchEngineParams, out double chargeCarrierMass)
        {
            string strValue = null;
            if (objSearchEngineParams.Parameters.TryGetValue(CHARGE_CARRIER_MASS_PARAM_NAME, out strValue))
            {
                if (double.TryParse(strValue, out chargeCarrierMass))
                {
                    return true;
                }
            }

            chargeCarrierMass = clsPeptideMassCalculator.MASS_PROTON;
            return false;
        }

        public static string GetPHRPFirstHitsFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_FHT;
        }

        public static string GetPHRPModSummaryFileName(string strDatasetName)
        {
            return strDatasetName + "_msgfplus_syn_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_msgfplus_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_msgfplus_syn_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_msgfplus_syn_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_msgfplus_syn_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_msgfplus_syn_SeqToProteinMap.txt";
        }

        public static string GetSearchEngineName()
        {
            return MSGFDB_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MSGFDB (aka MS-GF+) parameter file
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            bool blnSuccess = false;

            objSearchEngineParams = new clsSearchEngineParameters(MSGFDB_SEARCH_ENGINE_NAME, mModInfo);

            blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private bool ReadSearchEngineParamFile(string strSearchEngineParamFileName, clsSearchEngineParameters objSearchEngineParams)
        {
            try
            {
                mPeptideMassCalculator.ResetAminoAcidMasses();

                var success = ReadKeyValuePairSearchEngineParamFile(MSGFDB_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, clsPHRPReader.ePeptideHitResultType.MSGFDB, objSearchEngineParams);

                if (!success)
                {
                    return false;
                }

                var strSettingValue = string.Empty;
                int intValue = 0;

                // Determine the enzyme name
                if (objSearchEngineParams.Parameters.TryGetValue("enzymeid", out strSettingValue))
                {
                    if (int.TryParse(strSettingValue, out intValue))
                    {
                        switch (intValue)
                        {
                            case 0:
                                objSearchEngineParams.Enzyme = "no_enzyme";
                                break;
                            case 1:
                                objSearchEngineParams.Enzyme = "trypsin";
                                break;
                            case 2:
                                objSearchEngineParams.Enzyme = "Chymotrypsin";
                                break;
                            case 3:
                                objSearchEngineParams.Enzyme = "Lys-C";
                                break;
                            case 4:
                                objSearchEngineParams.Enzyme = "Lys-N";
                                break;
                            case 5:
                                objSearchEngineParams.Enzyme = "Glu-C";
                                break;
                            case 6:
                                objSearchEngineParams.Enzyme = "Arg-C";
                                break;
                            case 7:
                                objSearchEngineParams.Enzyme = "Asp-N";
                                break;
                            case 8:
                                objSearchEngineParams.Enzyme = "alphaLP";
                                break;
                            case 9:
                                objSearchEngineParams.Enzyme = "no_enzyme_peptidomics";
                                break;
                            default:
                                objSearchEngineParams.Enzyme = "unknown_enzyme";
                                break;
                        }
                    }
                }

                // Determine the cleavage specificity
                if (objSearchEngineParams.Parameters.TryGetValue("nnet", out strSettingValue))
                {
                    // NNET means number of non-enzymatic terminii

                    if (int.TryParse(strSettingValue, out intValue))
                    {
                        switch (intValue)
                        {
                            case 0:
                                // Fully-tryptic
                                objSearchEngineParams.MinNumberTermini = 2;
                                break;
                            case 1:
                                // Partially-tryptic
                                objSearchEngineParams.MinNumberTermini = 1;
                                break;
                            default:
                                // No-enzyme search
                                objSearchEngineParams.MinNumberTermini = 0;
                                break;
                        }
                    }
                }
                else
                {
                    // MSGF+ uses ntt instead of nnet; thus look for ntt

                    if (objSearchEngineParams.Parameters.TryGetValue("ntt", out strSettingValue))
                    {
                        // NTT means number of tolerable terminii

                        if (int.TryParse(strSettingValue, out intValue))
                        {
                            switch (intValue)
                            {
                                case 0:
                                    // No-enzyme search
                                    objSearchEngineParams.MinNumberTermini = 0;
                                    break;
                                case 1:
                                    // Partially-tryptic
                                    objSearchEngineParams.MinNumberTermini = 1;
                                    break;
                                default:
                                    // Fully-tryptic
                                    objSearchEngineParams.MinNumberTermini = 2;
                                    break;
                            }
                        }
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                double dblTolerancePPM = 0;
                objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, out dblTolerancePPM, clsPHRPReader.ePeptideHitResultType.MSGFDB);
                objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM;

                // Look for Custom Amino Acid definitions
                if (!objSearchEngineParams.Parameters.Any(paramEntry => paramEntry.Key == clsMSGFPlusParamFileModExtractor.PARAM_TAG_CUSTOMAA))
                {
                    // No custom amino acid entries
                    return true;
                }

                // Store the Custom Amino Acid info
                // Need to use a different parsing function to extract it
                success = UpdateMassCalculatorMasses(strSearchEngineParamFileName);

                // Look for a custom charge carrier mass
                double customChargeCarrierMass = 0;
                if (GetCustomChargeCarrierMass(objSearchEngineParams, out customChargeCarrierMass))
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
        /// <param name="strLine">Data line</param>
        /// <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="objPSM">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string strLine, int intLinesRead, out clsPSM objPSM, bool fastReadMode)
        {
            var strColumns = strLine.Split('\t');
            string strPeptide = null;
            string strProtein = null;

            double dblPrecursorMZ = 0;

            bool blnMSGFPlusResults = false;
            bool blnSuccess = false;

            objPSM = new clsPSM();

            try
            {
                if (clsPHRPReader.LookupColumnIndex(DATA_COLUMN_MSGFPlus_SpecEValue, mColumnHeaders) >= 0)
                {
                    blnMSGFPlusResults = true;
                }
                else
                {
                    blnMSGFPlusResults = false;
                }

                objPSM.DataLineText = strLine;
                objPSM.ScanNumber = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100);
                if (objPSM.ScanNumber == -100)
                {
                    // Data line is not valid
                }
                else
                {
                    objPSM.ResultID = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0);

                    if (blnMSGFPlusResults)
                    {
                        objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Rank_MSGFPlus_SpecEValue, mColumnHeaders, 1);
                    }
                    else
                    {
                        objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb, mColumnHeaders, 1);
                    }

                    strPeptide = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders);

                    if (fastReadMode)
                    {
                        objPSM.SetPeptide(strPeptide, blnUpdateCleanSequence: false);
                    }
                    else
                    {
                        objPSM.SetPeptide(strPeptide, mCleavageStateCalculator);
                    }

                    objPSM.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                    strProtein = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders);
                    objPSM.AddProtein(strProtein);

                    objPSM.CollisionMode = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_FragMethod, mColumnHeaders, "n/a");

                    dblPrecursorMZ = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0);
                    objPSM.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, objPSM.Charge, 0);

                    objPSM.MassErrorDa = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders);
                    objPSM.MassErrorPPM = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                    if (blnMSGFPlusResults)
                    {
                        objPSM.MSGFSpecEValue = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_MSGFPlus_SpecEValue, mColumnHeaders);
                    }
                    else
                    {
                        objPSM.MSGFSpecEValue = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_MSGFDB_SpecProb, mColumnHeaders);
                    }

                    if (objPSM.MSGFSpecEValue.Length > 13)
                    {
                        // Attempt to shorten the SpecEValue value
                        double dblSpecEValue;
                        if (double.TryParse(objPSM.MSGFSpecEValue, out dblSpecEValue))
                        {
                            objPSM.MSGFSpecEValue = dblSpecEValue.ToString("0.0000000E-00");
                        }
                    }
                    blnSuccess = true;
                }

                if (blnSuccess)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(objPSM);
                    }

                    // Store the remaining scores
                    AddScore(objPSM, strColumns, DATA_COLUMN_DeNovoScore);

                    AddScore(objPSM, strColumns, DATA_COLUMN_MSGFScore);

                    if (blnMSGFPlusResults)
                    {
                        AddScore(objPSM, strColumns, DATA_COLUMN_MSGFPlus_SpecEValue);
                        AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFPlus_SpecEValue);
                        AddScore(objPSM, strColumns, DATA_COLUMN_EValue);
                        AddScore(objPSM, strColumns, DATA_COLUMN_QValue);
                        AddScore(objPSM, strColumns, DATA_COLUMN_PepQValue);
                        AddScore(objPSM, strColumns, DATA_COLUMN_Isotope_Error);

                        string strValue = string.Empty;

                        // Duplicate the score values to provide backwards compatibility
                        if (objPSM.TryGetScore(DATA_COLUMN_MSGFPlus_SpecEValue, out strValue))
                            objPSM.SetScore(DATA_COLUMN_MSGFDB_SpecProb, strValue);
                        if (objPSM.TryGetScore(DATA_COLUMN_Rank_MSGFPlus_SpecEValue, out strValue))
                            objPSM.SetScore(DATA_COLUMN_Rank_MSGFDB_SpecProb, strValue);
                        if (objPSM.TryGetScore(DATA_COLUMN_QValue, out strValue))
                            objPSM.SetScore(DATA_COLUMN_FDR, strValue);
                        if (objPSM.TryGetScore(DATA_COLUMN_PepQValue, out strValue))
                            objPSM.SetScore(DATA_COLUMN_PepFDR, strValue);

                        var strEValue = string.Empty;

                        bool blnPValueStored = false;

                        if (objPSM.TryGetScore(DATA_COLUMN_EValue, out strEValue))
                        {
                            string strSpecEValue = string.Empty;
                            if (objPSM.TryGetScore(DATA_COLUMN_MSGFPlus_SpecEValue, out strSpecEValue))
                            {
                                // Compute PValue using EValue and SpecEValue
                                double dblEValue = 0;
                                if (double.TryParse(strEValue, out dblEValue))
                                {
                                    double dblSpecEValue = 0;
                                    if (double.TryParse(strSpecEValue, out dblSpecEValue))
                                    {
                                        if (dblSpecEValue > 0)
                                        {
                                            double dblN = dblEValue / dblSpecEValue;
                                            double dblPValue = 1 - Math.Pow((1 - dblSpecEValue), dblN);

                                            if (Math.Abs(dblPValue) <= double.Epsilon)
                                            {
                                                objPSM.SetScore(DATA_COLUMN_PValue, "0");
                                            }
                                            else
                                            {
                                                objPSM.SetScore(DATA_COLUMN_PValue, dblPValue.ToString("0.00000E-00"));
                                            }

                                            blnPValueStored = true;
                                        }
                                    }
                                }
                            }

                            if (!blnPValueStored)
                            {
                                // Store E-value as P-value (these values are not identical, and will only be close for high-confidence results, i.e. results with FDR < 2%)
                                objPSM.SetScore(DATA_COLUMN_PValue, strEValue);
                            }
                        }
                    }
                    else
                    {
                        AddScore(objPSM, strColumns, DATA_COLUMN_MSGFDB_SpecProb);
                        AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb);
                        AddScore(objPSM, strColumns, DATA_COLUMN_PValue);
                        AddScore(objPSM, strColumns, DATA_COLUMN_FDR);
                        AddScore(objPSM, strColumns, DATA_COLUMN_PepFDR);
                    }

                    AddScore(objPSM, strColumns, DATA_COLUMN_EFDR);      // This column will not be present if a Target/Decoy (TDA) search was performed

                    AddScore(objPSM, strColumns, DATA_COLUMN_IMS_Scan);

                    AddScore(objPSM, strColumns, DATA_COLUMN_IMS_Drift_Time);
                }
            }
            catch (Exception ex)
            {
                base.ReportError("Error parsing line " + intLinesRead + " in the MSGFDB data file: " + ex.Message);
            }

            return blnSuccess;
        }

        private bool UpdateMassCalculatorMasses(string strSearchEngineParamFileName)
        {
            string localErrorMsg = string.Empty;
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MSGF+");

            modFileProcessor.ErrorOccurred += ModExtractorErrorHandler;
            modFileProcessor.WarningMessageEvent += ModExtractorWarningHandler;

            var success = UpdateMassCalculatorMasses(strSearchEngineParamFileName, modFileProcessor, mPeptideMassCalculator, out localErrorMsg);

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
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="modFileProcessor"></param>
        /// <param name="peptideMassCalculator"></param>
        /// <param name="errorMessage"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static bool UpdateMassCalculatorMasses(
            string strSearchEngineParamFileName,
            clsMSGFPlusParamFileModExtractor modFileProcessor,
            clsPeptideMassCalculator peptideMassCalculator,
            out string errorMessage)
        {
            if (modFileProcessor == null)
            {
                throw new ObjectDisposedException("modFileProcessor is not initialized");
            }

            errorMessage = string.Empty;

            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo = null;

            // Note that this call will initialize lstModInfo
            var success = modFileProcessor.ExtractModInfoFromParamFile(strSearchEngineParamFileName, out lstModInfo);
            if (!success)
            {
                errorMessage = modFileProcessor.ErrorMessage;
                return false;
            }

            var customAminoAcidDefs = (from item in lstModInfo where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA select item).ToList();
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

        #region "Event Handlers"
        private void ModExtractorErrorHandler(string errMsg)
        {
            ReportError(errMsg);
        }

        private void ModExtractorWarningHandler(string warningMsg)
        {
            ReportError(warningMsg);
        }
        #endregion
    }
}
