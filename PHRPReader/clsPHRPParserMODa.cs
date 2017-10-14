//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/02/2014
//
// This class parses data lines from moda_syn.txt files
//
//*********************************************************************************************************
using System;
using System.Collections.Generic;

namespace PHRPReader
{
    public class clsPHRPParserMODa : clsPHRPParser
    {
        #region "Constants"
        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Spectrum_Index = "Spectrum_Index";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Score = "Score";
        public const string DATA_COLUMN_Probability = "Probability";
        public const string DATA_COLUMN_Rank_Probability = "Rank_Probability";
        public const string DATA_COLUMN_Peptide_Position = "Peptide_Position";
        public const string DATA_COLUMN_QValue = "QValue";

        public const string FILENAME_SUFFIX_SYN = "_moda_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_moda_fht.txt";

        private const string MODa_SEARCH_ENGINE_NAME = "MODa";
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
        public clsPHRPParserMODa(string strDatasetName, string strInputFilePath)
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
        public clsPHRPParserMODa(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MODa, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMODa(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MODa, startupOptions)
        {
        }

        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_ResultID);
            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_Spectrum_Index);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_PrecursorMZ);
            AddHeaderColumn(DATA_COLUMN_DelM);
            AddHeaderColumn(DATA_COLUMN_DelM_PPM);
            AddHeaderColumn(DATA_COLUMN_MH);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_Score);
            AddHeaderColumn(DATA_COLUMN_Probability);
            AddHeaderColumn(DATA_COLUMN_Rank_Probability);
            AddHeaderColumn(DATA_COLUMN_Peptide_Position);
            AddHeaderColumn(DATA_COLUMN_QValue);
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

            if (objSearchEngineParams.Parameters.TryGetValue("PPMTolerance", out strTolerance))
            {
                // Parent mass tolerance, in ppm
                if (double.TryParse(strTolerance, out dblTolerancePPM))
                {
                    dblToleranceDa = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000);
                }
            }
            else if (objSearchEngineParams.Parameters.TryGetValue("PeptTolerance", out strTolerance))
            {
                // Parent mass tolerance, in Da
                double.TryParse(strTolerance, out dblToleranceDa);

                // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000);
            }

            return dblToleranceDa;
        }

        public static string GetPHRPFirstHitsFileName(string strDatasetName)
        {
            // MODa does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        public static string GetPHRPModSummaryFileName(string strDatasetName)
        {
            return strDatasetName + "_moda_syn_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_moda_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_moda_syn_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_moda_syn_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_moda_syn_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_moda_syn_SeqToProteinMap.txt";
        }

        public static string GetSearchEngineName()
        {
            return MODa_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MODa parameter file
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            var blnSuccess = false;

            objSearchEngineParams = new clsSearchEngineParameters(MODa_SEARCH_ENGINE_NAME);

            blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private bool ReadSearchEngineParamFile(string strSearchEngineParamFileName, clsSearchEngineParameters objSearchEngineParams)
        {
            var strSettingValue = string.Empty;
            var objModDef = default(clsModificationDefinition);

            try
            {
                var blnSuccess = ReadKeyValuePairSearchEngineParamFile(MODa_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, clsPHRPReader.ePeptideHitResultType.MODa, objSearchEngineParams);

                if (!blnSuccess)
                {
                    return false;
                }

                // For MSGF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                // But MODa does not have a _ModDefs.txt file because it performs a blind search
                // The user can define static mods on any of the residues, plus the peptide terminii; check for these now

                var lstResiduesToFind = new List<string> {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

                // This dictionary tracks the static mod names we will look for
                // It is populated using the amino acid letters in lstResiduesToFind, plus also the N and T terminus tags
                var dctResiduesAndSymbols = new Dictionary<string, string>();

                foreach (var residueSymbol in lstResiduesToFind)
                {
                    dctResiduesAndSymbols.Add(residueSymbol, residueSymbol);
                }

                dctResiduesAndSymbols.Add("ADD_NTerm", clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString());
                dctResiduesAndSymbols.Add("ADD_CTerm", clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString());

                foreach (var residueSpec in dctResiduesAndSymbols)
                {
                    var strKey = "ADD_" + residueSpec.Key;

                    if (objSearchEngineParams.Parameters.TryGetValue(strKey, out strSettingValue))
                    {
                        double modMassDa = 0;

                        if (double.TryParse(strSettingValue, out modMassDa))
                        {
                            if (Math.Abs(modMassDa) > float.Epsilon)
                            {
                                var eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod;
                                if (residueSpec.Value == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || residueSpec.Value == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                                {
                                    eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                                }

                                objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, residueSpec.Value, eModType, "Mod" + modMassDa.ToString("0"));
                                objSearchEngineParams.AddModification(objModDef);
                            }
                        }
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                double dblTolerancePPM = 0;
                objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, out dblTolerancePPM);
                objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM;

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
                    objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Rank_Probability, mColumnHeaders, 1);

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
                    AddScore(objPSM, strColumns, DATA_COLUMN_Spectrum_Index);

                    AddScore(objPSM, strColumns, DATA_COLUMN_MH);

                    AddScore(objPSM, strColumns, DATA_COLUMN_Score);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Probability);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Position);
                    AddScore(objPSM, strColumns, DATA_COLUMN_QValue);
                }
            }
            catch (Exception ex)
            {
                base.ReportError("Error parsing line " + intLinesRead + " in the MODa data file: " + ex.Message);
            }

            return blnSuccess;
        }
    }
}
