//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from Inspect _inspect_syn.txt files
//
//*********************************************************************************************************
using System;
using System.Runtime.InteropServices;

namespace PHRPReader
{
    public class clsPHRPParserInspect : clsPHRPParser
    {
        #region "Constants"
        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_MQScore = "MQScore";
        public const string DATA_COLUMN_Length = "Length";
        public const string DATA_COLUMN_TotalPRMScore = "TotalPRMScore";
        public const string DATA_COLUMN_MedianPRMScore = "MedianPRMScore";
        public const string DATA_COLUMN_FractionY = "FractionY";
        public const string DATA_COLUMN_FractionB = "FractionB";
        public const string DATA_COLUMN_Intensity = "Intensity";
        public const string DATA_COLUMN_NTT = "NTT";
        public const string DATA_COLUMN_PValue = "PValue";
        public const string DATA_COLUMN_FScore = "FScore";
        public const string DATA_COLUMN_DeltaScore = "DeltaScore";
        public const string DATA_COLUMN_DeltaScoreOther = "DeltaScoreOther";
        public const string DATA_COLUMN_DeltaNormMQScore = "DeltaNormMQScore";
        public const string DATA_COLUMN_DeltaNormTotalPRMScore = "DeltaNormTotalPRMScore";
        public const string DATA_COLUMN_RankTotalPRMScore = "RankTotalPRMScore";
        public const string DATA_COLUMN_RankFScore = "RankFScore";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_RecordNumber = "RecordNumber";
        public const string DATA_COLUMN_DBFilePos = "DBFilePos";
        public const string DATA_COLUMN_SpecFilePos = "SpecFilePos";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_PrecursorError = "PrecursorError";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";

        public const string FILENAME_SUFFIX_SYN = "_inspect_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_inspect_fht.txt";

        private const string INS_SEARCH_ENGINE_NAME = "Inspect";
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
        public clsPHRPParserInspect(string strDatasetName, string strInputFilePath)
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
        public clsPHRPParserInspect(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.Inspect, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserInspect(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.Inspect, startupOptions)
        {
        }

        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_ResultID);
            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_MQScore);
            AddHeaderColumn(DATA_COLUMN_Length);
            AddHeaderColumn(DATA_COLUMN_TotalPRMScore);
            AddHeaderColumn(DATA_COLUMN_MedianPRMScore);
            AddHeaderColumn(DATA_COLUMN_FractionY);
            AddHeaderColumn(DATA_COLUMN_FractionB);
            AddHeaderColumn(DATA_COLUMN_Intensity);
            AddHeaderColumn(DATA_COLUMN_NTT);
            AddHeaderColumn(DATA_COLUMN_PValue);
            AddHeaderColumn(DATA_COLUMN_FScore);
            AddHeaderColumn(DATA_COLUMN_DeltaScore);
            AddHeaderColumn(DATA_COLUMN_DeltaScoreOther);
            AddHeaderColumn(DATA_COLUMN_DeltaNormMQScore);
            AddHeaderColumn(DATA_COLUMN_DeltaNormTotalPRMScore);
            AddHeaderColumn(DATA_COLUMN_RankTotalPRMScore);
            AddHeaderColumn(DATA_COLUMN_RankFScore);
            AddHeaderColumn(DATA_COLUMN_MH);
            AddHeaderColumn(DATA_COLUMN_RecordNumber);
            AddHeaderColumn(DATA_COLUMN_DBFilePos);
            AddHeaderColumn(DATA_COLUMN_SpecFilePos);
            AddHeaderColumn(DATA_COLUMN_PrecursorMZ);
            AddHeaderColumn(DATA_COLUMN_PrecursorError);
            AddHeaderColumn(DATA_COLUMN_DelM_PPM);
        }

        /// <summary>
        /// Determines the precursor mass tolerance
        /// </summary>
        /// <param name="objSearchEngineParams"></param>
        /// <param name="dblTolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <returns>Precursor tolerance, in Da</returns>
        /// <remarks></remarks>
        private double DeterminePrecursorMassTolerance(clsSearchEngineParameters objSearchEngineParams, out double dblTolerancePPM)
        {
            string strTolerance = string.Empty;

            double dblToleranceDa = 0;
            double dblTolerance = 0;

            dblTolerancePPM = 0;

            if (objSearchEngineParams.Parameters.TryGetValue("ParentPPM", out strTolerance))
            {
                // Parent mass tolerance, in ppm
                double.TryParse(strTolerance, out dblTolerancePPM);
                // Convert from PPM to dalton (assuming a mass of 2000 m/z)
                dblTolerance = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000);
            }

            if (objSearchEngineParams.Parameters.TryGetValue("PMTolerance", out strTolerance))
            {
                // Parent mass tolerance, in Da
                double.TryParse(strTolerance, out dblToleranceDa);

                // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000);
            }

            dblTolerance = Math.Max(dblTolerance, dblToleranceDa);

            return dblTolerance;
        }

        public static string GetPHRPFirstHitsFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_FHT;
        }

        public static string GetPHRPModSummaryFileName(string strDatasetName)
        {
            return strDatasetName + "_inspect_syn_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_inspect_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_inspect_syn_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_inspect_syn_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_inspect_syn_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_inspect_syn_SeqToProteinMap.txt";
        }

        public static string GetSearchEngineName()
        {
            return INS_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified Inspect parameter file
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            bool blnSuccess = false;

            objSearchEngineParams = new clsSearchEngineParameters(INS_SEARCH_ENGINE_NAME, mModInfo);

            blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private bool ReadSearchEngineParamFile(string strSearchEngineParamFileName, clsSearchEngineParameters objSearchEngineParams)
        {
            string strSettingValue = string.Empty;
            bool blnSuccess = false;

            try
            {
                blnSuccess = ReadKeyValuePairSearchEngineParamFile(INS_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, clsPHRPReader.ePeptideHitResultType.Inspect, objSearchEngineParams);

                if (blnSuccess)
                {
                    // Determine the enzyme name
                    if (objSearchEngineParams.Parameters.TryGetValue("protease", out strSettingValue))
                    {
                        switch (strSettingValue.ToLower())
                        {
                            case "trypsin":
                                objSearchEngineParams.Enzyme = "trypsin";
                                break;
                            case "none":
                                objSearchEngineParams.Enzyme = "no_enzyme";
                                break;
                            case "chymotrypsin":
                                objSearchEngineParams.Enzyme = "chymotrypsin";
                                break;
                            default:
                                if (!string.IsNullOrEmpty(strSettingValue))
                                {
                                    objSearchEngineParams.Enzyme = strSettingValue;
                                }
                                break;
                        }
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    double dblTolerancePPM = 0;
                    objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, out dblTolerancePPM);
                    objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM;
                }
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
            }

            return blnSuccess;
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
            string[] strColumns = strLine.Split('\t');
            string strPeptide = null;
            string strProtein = null;
            double dblPrecursorMZ = 0;

            bool blnSuccess = false;

            objPSM = new clsPSM();

            try
            {
                objPSM.DataLineText = strLine;
                objPSM.ScanNumber = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100);
                if (objPSM.ScanNumber == -100)
                {
                    // Data line is not valid
                }
                else
                {
                    objPSM.ResultID = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0);
                    objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_RankTotalPRMScore, mColumnHeaders, 0);

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

                    dblPrecursorMZ = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0);
                    objPSM.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, objPSM.Charge, 0);

                    objPSM.MassErrorDa = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_PrecursorError, mColumnHeaders);
                    objPSM.MassErrorPPM = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                    blnSuccess = true;
                }

                if (blnSuccess)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(objPSM);
                    }

                    // Store the remaining scores
                    AddScore(objPSM, strColumns, DATA_COLUMN_MQScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_TotalPRMScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_MedianPRMScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_PValue);
                    AddScore(objPSM, strColumns, DATA_COLUMN_FScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DeltaScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DeltaScoreOther);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DeltaNormMQScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DeltaNormTotalPRMScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_RankTotalPRMScore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_RankFScore);
                }
            }
            catch (Exception ex)
            {
                base.ReportError("Error parsing line " + intLinesRead + " in the Inspect data file: " + ex.Message);
            }

            return blnSuccess;
        }
    }
}
