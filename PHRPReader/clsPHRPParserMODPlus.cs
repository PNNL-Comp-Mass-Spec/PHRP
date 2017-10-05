//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 05/18/2015
//
// This class parses data lines from modp_syn.txt files
//
//*********************************************************************************************************
using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Xml;

namespace PHRPReader
{
    public class clsPHRPParserMODPlus : clsPHRPParser
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
        public const string DATA_COLUMN_NTT = "NTT";
        public const string DATA_COLUMN_Modification_Annotation = "Modification_Annotation";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Peptide_Position = "Peptide_Position";
        public const string DATA_COLUMN_Score = "Score";
        public const string DATA_COLUMN_Probability = "Probability";
        public const string DATA_COLUMN_Rank_Score = "Rank_Score";
        public const string DATA_COLUMN_QValue = "QValue";

        public const string FILENAME_SUFFIX_SYN = "_modp_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_modp_fht.txt";

        private const string MODplus_SEARCH_ENGINE_NAME = "MODplus";

        private const string PEPTIDE_MASS_TOL_PPM = "PeptideMassTolPPM";
        private const string PEPTIDE_MASS_TOL_DA = "PeptideMassTolDa";
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
        public clsPHRPParserMODPlus(string strDatasetName, string strInputFilePath)
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
        public clsPHRPParserMODPlus(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MODPlus, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMODPlus(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.MODPlus, startupOptions)
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
            AddHeaderColumn(DATA_COLUMN_NTT);
            AddHeaderColumn(DATA_COLUMN_Modification_Annotation);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_Peptide_Position);
            AddHeaderColumn(DATA_COLUMN_Score);
            AddHeaderColumn(DATA_COLUMN_Probability);
            AddHeaderColumn(DATA_COLUMN_Rank_Score);
            AddHeaderColumn(DATA_COLUMN_QValue);
        }

        public static string GetPHRPFirstHitsFileName(string strDatasetName)
        {
            // MODplus does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        public static string GetPHRPModSummaryFileName(string strDatasetName)
        {
            return strDatasetName + "_modp_syn_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_modp_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_modp_syn_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_modp_syn_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_modp_syn_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_modp_syn_SeqToProteinMap.txt";
        }

        public static string GetSearchEngineName()
        {
            return MODplus_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MODp parameter file
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            bool blnSuccess = false;

            objSearchEngineParams = new clsSearchEngineParameters(MODplus_SEARCH_ENGINE_NAME);

            blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private bool ReadSearchEngineParamFile(string strSearchEngineParamFileName, clsSearchEngineParameters objSearchEngineParams)
        {
            try
            {
                // For MSGF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                // But MODplus does not have a _ModDefs.txt file because it performs a blind search
                // The user can define static mods on any of the residues, plus the peptide terminii; check for these now

                var strParamFilePath = Path.Combine(mInputFolderPath, strSearchEngineParamFileName);

                // Read the contents of the parameter file
                var doc = new XmlDocument();
                doc.Load(strParamFilePath);

                var dbNodes = doc.SelectNodes("/search/database");
                if (dbNodes.Count > 0)
                {
                    objSearchEngineParams.FastaFilePath = GetAttribute(dbNodes[0], "local_path");
                }

                var enzymeNodes = doc.SelectNodes("/search/enzyme_rule");
                if (enzymeNodes.Count > 0)
                {
                    objSearchEngineParams.Enzyme = GetAttribute(enzymeNodes[0], "name");
                }

                var instrumentResolutionNodes = doc.SelectNodes("/search/instrument_resolution");
                if (instrumentResolutionNodes.Count > 0)
                {
                    objSearchEngineParams.PrecursorMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes[0], "ms"));
                    objSearchEngineParams.FragmentMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes[0], "msms"));
                }

                var enzymeConstraintNodes = doc.SelectNodes("/search/parameters/enzyme_constraint");
                if (enzymeConstraintNodes.Count > 0)
                {
                    if (int.TryParse(GetAttribute(enzymeConstraintNodes[0], "max_miss_cleavages"), out var maxNumberInternalCleavages))
                    {
                        objSearchEngineParams.MaxNumberInternalCleavages = maxNumberInternalCleavages;
                    }
                    if (int.TryParse(GetAttribute(enzymeConstraintNodes[0], "min_number_termini"), out var minNumberTermini))
                    {
                        objSearchEngineParams.MinNumberTermini = minNumberTermini;
                    }
                }

                var massTolParamNode = doc.SelectNodes("/search/parameters/peptide_mass_tol");

                if (massTolParamNode.Count > 0)
                {
                    var strTolerance = GetAttribute(massTolParamNode[0], "value");
                    var massUnits = GetAttribute(massTolParamNode[0], "unit");

                    if (massUnits.ToLower() == "ppm")
                    {
                        // Parent mass tolerance, in ppm

                        if (double.TryParse(strTolerance, out var precursorMassTolerancePpm))
                        {
                            objSearchEngineParams.PrecursorMassTolerancePpm = precursorMassTolerancePpm;
                            objSearchEngineParams.PrecursorMassToleranceDa = clsPeptideMassCalculator.PPMToMass(objSearchEngineParams.PrecursorMassTolerancePpm, 2000);
                        }
                    }
                    else if (massUnits.ToLower() == "da")
                    {
                        if (double.TryParse(strTolerance, out var precursorMassToleranceDa))
                        {
                            objSearchEngineParams.PrecursorMassToleranceDa = precursorMassToleranceDa;
                        }

                        // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                        objSearchEngineParams.PrecursorMassToleranceDa = clsPeptideMassCalculator.MassToPPM(objSearchEngineParams.PrecursorMassToleranceDa, 2000);
                    }
                }

                var modNodes = doc.SelectNodes("/search/modifications/fixed/mod");
                if (modNodes.Count > 0)
                {
                    // Store the fixed mods

                    foreach (XmlNode node in modNodes)
                    {
                        var modName = GetAttribute(node, "name");
                        var strResidue = GetAttribute(node, "site").Trim();
                        var modPosition = GetAttribute(node, "position");
                        var modMass = GetAttribute(node, "massdiff");

                        // Replace N-Term or C-Term with < or >
                        if (strResidue.ToLower() == "n-term")
                            strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        if (strResidue.ToLower() == "c-term")
                            strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                        double modMassDa = 0;

                        if (double.TryParse(modMass, out modMassDa))
                        {
                            if (Math.Abs(modMassDa - 0) > float.Epsilon)
                            {
                                var eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod;
                                if (strResidue == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || strResidue == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                                {
                                    eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                                }

                                var objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, strResidue, eModType, "Mod" + modMassDa.ToString("0"));
                                objSearchEngineParams.AddModification(objModDef);
                            }
                        }
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }

        private string ConvertResolutionModeToMassType(string resolutionType)
        {
            if (resolutionType == "high")
            {
                return clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
            }

            if (resolutionType == "low")
            {
                return clsSearchEngineParameters.MASS_TYPE_AVERAGE;
            }

            return "unknown";
        }

        private string GetAttribute(XmlNode node, string attributeName)
        {
            if (node.Attributes.Count > 0)
            {
                try
                {
                    var attribute = node.Attributes[attributeName];

                    if ((attribute != null))
                    {
                        return attribute.Value;
                    }
                }
                catch (Exception ex)
                {
                    // Ignore errors
                    Console.WriteLine("Attribute lookup error: " + ex.Message);
                }
            }

            return string.Empty;
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
                    objPSM.ScoreRank = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Rank_Score, mColumnHeaders, 1);

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

                    AddScore(objPSM, strColumns, DATA_COLUMN_Modification_Annotation);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Position);

                    AddScore(objPSM, strColumns, DATA_COLUMN_Score);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Probability);
                    AddScore(objPSM, strColumns, DATA_COLUMN_QValue);
                }
            }
            catch (Exception ex)
            {
                base.ReportError("Error parsing line " + intLinesRead + " in the MODPlus data file: " + ex.Message);
            }

            return blnSuccess;
        }
    }
}
