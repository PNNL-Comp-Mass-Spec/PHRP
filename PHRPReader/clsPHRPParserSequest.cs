//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from Sequest _syn.txt and _fht.txt ifles
//
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace PHRPReader
{
    public class clsPHRPParserSequest : clsPHRPParser
    {
        #region "Constants"
        public const string DATA_COLUMN_HitNum = "HitNum";
        public const string DATA_COLUMN_ScanNum = "ScanNum";
        public const string DATA_COLUMN_ScanCount = "ScanCount";
        public const string DATA_COLUMN_ChargeState = "ChargeState";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_XCorr = "XCorr";
        public const string DATA_COLUMN_DelCn = "DelCn";
        public const string DATA_COLUMN_Sp = "Sp";
        public const string DATA_COLUMN_Reference = "Reference";
        public const string DATA_COLUMN_MultiProtein = "MultiProtein";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_DelCn2 = "DelCn2";
        public const string DATA_COLUMN_RankSp = "RankSp";
        public const string DATA_COLUMN_RankXc = "RankXc";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_XcRatio = "XcRatio";
        public const string DATA_COLUMN_PassFilt = "PassFilt";
        public const string DATA_COLUMN_MScore = "MScore";
        public const string DATA_COLUMN_Ions_Observed = "Ions_Observed";
        public const string DATA_COLUMN_Ions_Expected = "Ions_Expected";
        public const string DATA_COLUMN_NumTrypticEnds = "NumTrypticEnds";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";

        public const string FILENAME_SUFFIX_SYN = "_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_fht.txt";

        private const string SEQ_SEARCH_ENGINE_NAME = "SEQUEST";
        #endregion

        #region "Properties"

        public override string PHRPFirstHitsFileName => GetPHRPFirstHitsFileName(mDatasetName);

        public override string PHRPModSummaryFileName => GetPHRPModSummaryFileName(mDatasetName);

        public override string PHRPPepToProteinMapFileName => GetPHRPPepToProteinMapFileName(mDatasetName);

        public override string PHRPProteinModsFileName => GetPHRPProteinModsFileName(mDatasetName);

        public override string PHRPSynopsisFileName => GetPHRPSynopsisFileName(mDatasetName);

        public override string PHRPResultToSeqMapFileName => GetPHRPResultToSeqMapFileName(mDatasetName);

        public override string PHRPSeqInfoFileName => GetPHRPSeqInfoFileName(mDatasetName);

        public override string PHRPSeqToProteinMapFileName => GetPHRPSeqToProteinMapFileName(mDatasetName);

        public override string SearchEngineName => GetSearchEngineName();

        #endregion

        /// <summary>
        /// Constructor; assumes blnLoadModsAndSeqInfo=True
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <remarks></remarks>
        public clsPHRPParserSequest(string strDatasetName, string strInputFilePath)
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
        public clsPHRPParserSequest(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.Sequest, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserSequest(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.Sequest, startupOptions)
        {
        }

        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_HitNum);
            AddHeaderColumn(DATA_COLUMN_ScanNum);
            AddHeaderColumn(DATA_COLUMN_ScanCount);
            AddHeaderColumn(DATA_COLUMN_ChargeState);
            AddHeaderColumn(DATA_COLUMN_MH);
            AddHeaderColumn(DATA_COLUMN_XCorr);
            AddHeaderColumn(DATA_COLUMN_DelCn);
            AddHeaderColumn(DATA_COLUMN_Sp);
            AddHeaderColumn(DATA_COLUMN_Reference);
            AddHeaderColumn(DATA_COLUMN_MultiProtein);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_DelCn2);
            AddHeaderColumn(DATA_COLUMN_RankSp);
            AddHeaderColumn(DATA_COLUMN_RankXc);
            AddHeaderColumn(DATA_COLUMN_DelM);
            AddHeaderColumn(DATA_COLUMN_XcRatio);

            mColumnHeaders.Add(DATA_COLUMN_PassFilt, -1);
            mColumnHeaders.Add(DATA_COLUMN_MScore, -1);

            AddHeaderColumn(DATA_COLUMN_Ions_Observed);
            AddHeaderColumn(DATA_COLUMN_Ions_Expected);

            AddHeaderColumn(DATA_COLUMN_NumTrypticEnds);
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
            var strPeptideMassTolerance = string.Empty;
            var strPeptideMassUnits = string.Empty;

            double dblToleranceDa = 0;
            double dblValue = 0;

            dblTolerancePPM = 0;

            if (objSearchEngineParams.Parameters.TryGetValue("peptide_mass_tolerance", out strPeptideMassTolerance))
            {
                if (double.TryParse(strPeptideMassTolerance, out dblValue))
                {
                    // Determine the mass units
                    // 0 means Da, 1 means mmu, 2 means ppm
                    var intUnits = 0;

                    if (objSearchEngineParams.Parameters.TryGetValue("peptide_mass_units", out strPeptideMassUnits))
                    {
                        if (!string.IsNullOrEmpty(strPeptideMassUnits))
                        {
                            int.TryParse(strPeptideMassUnits, out intUnits);
                        }
                    }

                    if (intUnits == 2)
                    {
                        // Tolerance is in ppm; convert to Da at 2000 m/z
                        dblTolerancePPM = dblValue;

                        dblToleranceDa = clsPeptideMassCalculator.PPMToMass(dblValue, 2000);
                    }
                    else if (intUnits == 1)
                    {
                        // Tolerance is in milli mass units
                        dblToleranceDa = dblValue / 1000.0;

                        // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                        dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000);
                    }
                    else
                    {
                        // Tolerance is in daltons
                        dblToleranceDa = dblValue;

                        // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                        dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000);
                    }
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
            return strDatasetName + "_syn_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_syn_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_syn_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_syn_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_syn_SeqToProteinMap.txt";
        }

        public static string GetSearchEngineName()
        {
            return SEQ_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified Sequest parameter file
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            var blnSuccess = false;

            objSearchEngineParams = new clsSearchEngineParameters(SEQ_SEARCH_ENGINE_NAME, mModInfo);
            objSearchEngineParams.Enzyme = "trypsin";

            blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private bool ReadSearchEngineParamFile(string strSearchEngineParamFileName, clsSearchEngineParameters objSearchEngineParams)
        {
            string strParamFilePath = null;
            string strLineIn = null;

            string strFastaFilePath = null;
            string strSettingValue = null;

            var intCharIndex = 0;
            var intValue = 0;

            var kvSetting = default(KeyValuePair<string, string>);

            var reEnzymeSpecificity = new Regex(@"^\S+\s(\d)\s\d\s.+", RegexOptions.Compiled | RegexOptions.IgnoreCase);
            var reMatch = default(Match);

            var blnSuccess = false;

            try
            {
                strParamFilePath = Path.Combine(mInputFolderPath, strSearchEngineParamFileName);

                if (!File.Exists(strParamFilePath))
                {
                    ReportError("Sequest param file not found: " + strParamFilePath);
                }
                else
                {
                    using (var srInFile = new StreamReader(new FileStream(strParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!srInFile.EndOfStream)
                        {
                            strLineIn = srInFile.ReadLine().TrimStart();

                            if (!string.IsNullOrWhiteSpace(strLineIn) && !strLineIn.StartsWith(";") && !strLineIn.StartsWith("[") && strLineIn.Contains("="))
                            {
                                // Split the line on the equals sign
                                kvSetting = ParseKeyValueSetting(strLineIn, '=');

                                if (!string.IsNullOrEmpty(kvSetting.Key))
                                {
                                    // Trim off any text that occurs after a semicolon in kvSetting.Value
                                    strSettingValue = kvSetting.Value;
                                    intCharIndex = strSettingValue.IndexOf(';');
                                    if (intCharIndex > 0)
                                    {
                                        strSettingValue = strSettingValue.Substring(intCharIndex).Trim();
                                    }

                                    objSearchEngineParams.AddUpdateParameter(kvSetting.Key, strSettingValue);

                                    switch (kvSetting.Key.ToLower())
                                    {
                                        case "first_database_name":
                                        case "database_name":
                                            try
                                            {
                                                strFastaFilePath = Path.Combine("C:\\Database", Path.GetFileName(strSettingValue));
                                            }
                                            catch (Exception)
                                            {
                                                strFastaFilePath = strSettingValue;
                                            }
                                            objSearchEngineParams.FastaFilePath = strFastaFilePath;

                                            break;
                                        case "mass_type_parent":
                                            if (strSettingValue == "0")
                                            {
                                                // Average mass
                                                objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_AVERAGE;
                                            }
                                            else
                                            {
                                                // Monoisotopic mass
                                                objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
                                            }

                                            break;
                                        case "mass_type_fragment":
                                            if (strSettingValue == "0")
                                            {
                                                // Average mass
                                                objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_AVERAGE;
                                            }
                                            else
                                            {
                                                // Monoisotopic mass
                                                objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
                                            }

                                            break;
                                        case "max_num_internal_cleavage_sites":
                                            if (int.TryParse(strSettingValue, out intValue))
                                            {
                                                objSearchEngineParams.MaxNumberInternalCleavages = intValue;
                                            }

                                            break;
                                        case "enzyme_info":
                                            // Used in new-style sequest parameter files

                                            // Examples include:
                                            // Fully-tryptic:     Trypsin(KR) 1 1 KR -
                                            // Partially-tryptic: Trypsin(KR) 2 1 KR -
                                            // No-enzyme:         No_Enzyme(-) 0 0 - -
                                            //
                                            objSearchEngineParams.Enzyme = "trypsin";

                                            if (strSettingValue.StartsWith("no_enzyme", StringComparison.InvariantCultureIgnoreCase))
                                            {
                                                objSearchEngineParams.MinNumberTermini = 0;
                                            }
                                            else
                                            {
                                                // Parse out the cleavage specificity number
                                                // This is the first number after the closing parenthesis in the above examples
                                                reMatch = reEnzymeSpecificity.Match(strSettingValue);
                                                if (reMatch != null && reMatch.Success)
                                                {
                                                    if (int.TryParse(reMatch.Groups[1].Value, out intValue))
                                                    {
                                                        objSearchEngineParams.MinNumberTermini = intValue;
                                                    }
                                                }
                                            }

                                            break;
                                        case "enzyme_number":
                                            // Used in old-style sequest parameter files
                                            if (int.TryParse(strSettingValue, out intValue))
                                            {
                                                if (intValue == 0)
                                                {
                                                    // No-enzyme
                                                    objSearchEngineParams.Enzyme = "trypsin";
                                                    objSearchEngineParams.MinNumberTermini = 0;
                                                }
                                                else
                                                {
                                                    switch (intValue)
                                                    {
                                                        case 1:
                                                            objSearchEngineParams.Enzyme = "trypsin";
                                                            break;
                                                        case 2:
                                                            objSearchEngineParams.Enzyme = "trypsin_modified";
                                                            break;
                                                        case 3:
                                                            objSearchEngineParams.Enzyme = "Chymotrypsin";
                                                            break;
                                                        case 4:
                                                            objSearchEngineParams.Enzyme = "Chymotrypsin_modified";
                                                            break;
                                                        case 5:
                                                            objSearchEngineParams.Enzyme = "Clostripain";
                                                            break;
                                                        case 6:
                                                            objSearchEngineParams.Enzyme = "Cyanogen_Bromide";
                                                            break;
                                                        case 7:
                                                            objSearchEngineParams.Enzyme = "IodosoBenzoate";
                                                            break;
                                                        case 8:
                                                            objSearchEngineParams.Enzyme = "Proline_Endopept";
                                                            break;
                                                        case 9:
                                                            objSearchEngineParams.Enzyme = "Staph_Protease";
                                                            break;
                                                        case 10:
                                                            objSearchEngineParams.Enzyme = "Trypsin_K";
                                                            break;
                                                        case 11:
                                                            objSearchEngineParams.Enzyme = "Trypsin_R";
                                                            break;
                                                        case 12:
                                                            objSearchEngineParams.Enzyme = "GluC";
                                                            break;
                                                        case 13:
                                                            objSearchEngineParams.Enzyme = "LysC";
                                                            break;
                                                        case 14:
                                                            objSearchEngineParams.Enzyme = "AspN";
                                                            break;
                                                        case 15:
                                                            objSearchEngineParams.Enzyme = "Elastase";
                                                            break;
                                                        case 16:
                                                            objSearchEngineParams.Enzyme = "Elastase/Tryp/Chymo";
                                                            break;
                                                        default:
                                                            objSearchEngineParams.Enzyme = "Unknown";
                                                            break;
                                                    }
                                                    objSearchEngineParams.MinNumberTermini = 2;
                                                }
                                            }
                                            break;
                                    }
                                }
                            }
                        }
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    double dblTolerancePPM = 0;
                    objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, out dblTolerancePPM);
                    objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM;

                    blnSuccess = true;
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

            double dblPrecursorMH = 0;
            double dblMassErrorDa = 0;

            var blnSuccess = false;

            objPSM = new clsPSM();

            try
            {
                objPSM.DataLineText = strLine;
                objPSM.ScanNumber = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_ScanNum, mColumnHeaders, -100);
                if (objPSM.ScanNumber == -100)
                {
                    // Data line is not valid
                }
                else
                {
                    objPSM.ResultID = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_HitNum, mColumnHeaders, 0);
                    objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_RankXc, mColumnHeaders, 1);

                    strPeptide = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders);

                    if (fastReadMode)
                    {
                        objPSM.SetPeptide(strPeptide, blnUpdateCleanSequence: false);
                    }
                    else
                    {
                        objPSM.SetPeptide(strPeptide, mCleavageStateCalculator);
                    }

                    objPSM.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_ChargeState, mColumnHeaders, 0));

                    strProtein = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Reference, mColumnHeaders);
                    objPSM.AddProtein(strProtein);

                    // Note that the MH value listed in Sequest files is not the precursor MH but is instead the theoretical (computed) MH of the peptide
                    // We'll update this value below using dblMassErrorDa
                    // We'll further update this value using the ScanStatsEx data
                    dblPrecursorMH = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_MH, mColumnHeaders, 0.0);
                    objPSM.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPrecursorMH, 1, 0);

                    objPSM.MassErrorDa = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders);
                    if (double.TryParse(objPSM.MassErrorDa, out dblMassErrorDa))
                    {
                        // Adjust the precursor mass
                        objPSM.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPrecursorMH - dblMassErrorDa, 1, 0);
                    }

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
                    AddScore(objPSM, strColumns, DATA_COLUMN_XCorr);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DelCn);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Sp);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DelCn2);
                    AddScore(objPSM, strColumns, DATA_COLUMN_RankSp);
                    AddScore(objPSM, strColumns, DATA_COLUMN_RankXc);
                    AddScore(objPSM, strColumns, DATA_COLUMN_XcRatio);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Ions_Observed);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Ions_Expected);
                    AddScore(objPSM, strColumns, DATA_COLUMN_NumTrypticEnds);
                }
            }
            catch (Exception ex)
            {
                base.ReportError("Error parsing line " + intLinesRead + " in the Sequest data file: " + ex.Message);
            }

            return blnSuccess;
        }
    }
}
