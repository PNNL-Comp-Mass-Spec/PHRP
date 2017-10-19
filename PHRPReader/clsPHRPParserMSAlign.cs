//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 11/28/2012
//
// This class parses data lines from MSAlign msalign_syn.txt files
//
//*********************************************************************************************************
using System;

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for MSAlign
    /// </summary>
    public class clsPHRPParserMSAlign : clsPHRPParser
    {
        #region "Constants"

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
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Protein_Mass = "Protein_Mass";
        public const string DATA_COLUMN_Unexpected_Mod_Count = "Unexpected_Mod_Count";
        public const string DATA_COLUMN_Peak_Count = "Peak_Count";
        public const string DATA_COLUMN_Matched_Peak_Count = "Matched_Peak_Count";
        public const string DATA_COLUMN_Matched_Fragment_Ion_Count = "Matched_Fragment_Ion_Count";
        public const string DATA_COLUMN_PValue = "PValue";
        public const string DATA_COLUMN_Rank_PValue = "Rank_PValue";
        public const string DATA_COLUMN_EValue = "EValue";
        public const string DATA_COLUMN_FDR = "FDR";
        public const string DATA_COLUMN_Species_ID = "Species_ID";
        public const string DATA_COLUMN_FragMethod = "FragMethod";

        public const string FILENAME_SUFFIX_SYN = "_msalign_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_msalign_fht.txt";

        private const string MSAlign_SEARCH_ENGINE_NAME = "MSAlign";
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
        /// Constructor; assumes blnLoadModsAndSeqInfo=True
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <remarks></remarks>
        public clsPHRPParserMSAlign(string strDatasetName, string strInputFilePath)
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
        public clsPHRPParserMSAlign(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MSAlign, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMSAlign(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MSAlign, startupOptions)
        {
        }

        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_ResultID);
            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_Prsm_ID);
            AddHeaderColumn(DATA_COLUMN_Spectrum_ID);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_PrecursorMZ);
            AddHeaderColumn(DATA_COLUMN_DelM);
            AddHeaderColumn(DATA_COLUMN_DelM_PPM);
            AddHeaderColumn(DATA_COLUMN_MH);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_Protein_Mass);
            AddHeaderColumn(DATA_COLUMN_Unexpected_Mod_Count);
            AddHeaderColumn(DATA_COLUMN_Peak_Count);
            AddHeaderColumn(DATA_COLUMN_Matched_Peak_Count);
            AddHeaderColumn(DATA_COLUMN_Matched_Fragment_Ion_Count);
            AddHeaderColumn(DATA_COLUMN_PValue);
            AddHeaderColumn(DATA_COLUMN_Rank_PValue);
            AddHeaderColumn(DATA_COLUMN_EValue);
            AddHeaderColumn(DATA_COLUMN_FDR);
            AddHeaderColumn(DATA_COLUMN_Species_ID);
            AddHeaderColumn(DATA_COLUMN_FragMethod);
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
            var strTolerance = string.Empty;

            double dblToleranceDa = 0;

            dblTolerancePPM = 0;

            if (objSearchEngineParams.Parameters.TryGetValue("errorTolerance", out strTolerance))
            {
                // Parent mass tolerance, in ppm
                if (double.TryParse(strTolerance, out dblTolerancePPM))
                {
                    dblToleranceDa = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000);
                }
            }

            return dblToleranceDa;
        }

        public static string GetPHRPFirstHitsFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_FHT;
        }

        public static string GetPHRPModSummaryFileName(string strDatasetName)
        {
            return strDatasetName + "_msalign_syn_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_msalign_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_msalign_syn_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_msalign_syn_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_msalign_syn_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_msalign_syn_SeqToProteinMap.txt";
        }

        public static string GetSearchEngineName()
        {
            return MSAlign_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MSAlign parameter file
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            var blnSuccess = false;

            objSearchEngineParams = new clsSearchEngineParameters(MSAlign_SEARCH_ENGINE_NAME);

            blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private bool ReadSearchEngineParamFile(string strSearchEngineParamFileName, clsSearchEngineParameters objSearchEngineParams)
        {
            var strSettingValue = string.Empty;
            var objModDef = default(clsModificationDefinition);
            var blnSuccess = false;

            try
            {
                blnSuccess = ReadKeyValuePairSearchEngineParamFile(MSAlign_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, clsPHRPReader.ePeptideHitResultType.MSAlign, objSearchEngineParams);

                if (blnSuccess)
                {
                    // For MSGF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                    // But MSAlign does not have a _ModDefs.txt file because it performs a blind search
                    // The user can define static mods on cysteine; check for these now

                    if (objSearchEngineParams.Parameters.TryGetValue("cysteineProtection", out strSettingValue))
                    {
                        switch (strSettingValue)
                        {
                            case "C57":
                                objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 57.0215, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcet");
                                objSearchEngineParams.AddModification(objModDef);
                                break;
                            case "C58":
                                objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 58.0055, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcid");
                                objSearchEngineParams.AddModification(objModDef);
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
            var strColumns = strLine.Split('\t');
            string strPeptide = null;
            string strProtein = null;

            double dblPrecursorMZ = 0;

            var blnSuccess = false;

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
                    objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Rank_PValue, mColumnHeaders, 1);

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

                    objPSM.MassErrorDa = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders);
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
                    AddScore(objPSM, strColumns, DATA_COLUMN_Prsm_ID);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Spectrum_ID);

                    AddScore(objPSM, strColumns, DATA_COLUMN_MH);

                    AddScore(objPSM, strColumns, DATA_COLUMN_Protein_Mass);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Unexpected_Mod_Count);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Peak_Count);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Matched_Peak_Count);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Matched_Fragment_Ion_Count);

                    AddScore(objPSM, strColumns, DATA_COLUMN_PValue);
                    AddScore(objPSM, strColumns, DATA_COLUMN_EValue);
                    AddScore(objPSM, strColumns, DATA_COLUMN_FDR);

                    AddScore(objPSM, strColumns, DATA_COLUMN_Species_ID);
                    AddScore(objPSM, strColumns, DATA_COLUMN_FragMethod);
                }
            }
            catch (Exception ex)
            {
                base.ReportError("Error parsing line " + intLinesRead + " in the MSAlign data file: " + ex.Message);
            }

            return blnSuccess;
        }
    }
}
