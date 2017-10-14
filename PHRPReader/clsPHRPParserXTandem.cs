//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from Xtandem _xt.txt files
//
// Note: in order to fully extract the search parameters, you need to have these files in the same folder as the _xt.txt file
//	The file passed to LoadSearchEngineParameters()  (typically input.xml)
//	The taxonomy file that it references             (typically taxonomy.xml)
//   The default input file defined in input.xml      (typically default_input.xml)
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;

namespace PHRPReader
{
    public class clsPHRPParserXTandem : clsPHRPParser
    {
        #region "Constants"
        public const string DATA_COLUMN_Result_ID = "Result_ID";
        public const string DATA_COLUMN_Group_ID = "Group_ID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_Peptide_MH = "Peptide_MH";
        public const string DATA_COLUMN_Peptide_Hyperscore = "Peptide_Hyperscore";
        public const string DATA_COLUMN_Peptide_Expectation_Value_LogE = "Peptide_Expectation_Value_Log(e)";
        public const string DATA_COLUMN_Multiple_Protein_Count = "Multiple_Protein_Count";
        public const string DATA_COLUMN_Peptide_Sequence = "Peptide_Sequence";
        public const string DATA_COLUMN_DeltaCn2 = "DeltaCn2";
        public const string DATA_COLUMN_y_score = "y_score";
        public const string DATA_COLUMN_y_ions = "y_ions";
        public const string DATA_COLUMN_b_score = "b_score";
        public const string DATA_COLUMN_b_ions = "b_ions";
        public const string DATA_COLUMN_Delta_Mass = "Delta_Mass";
        public const string DATA_COLUMN_Peptide_Intensity_LogI = "Peptide_Intensity_Log(I)";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";

        public const string FILENAME_SUFFIX_SYN = "_xt.txt";
        public const string FILENAME_SUFFIX_FHT = "_xt.txt";

        private const string XT_SEARCH_ENGINE_NAME = "X! Tandem";

        private const string TAXONOMY_INFO_KEY_NAME = "list path, taxonomy information";
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
        public clsPHRPParserXTandem(string strDatasetName, string strInputFilePath)
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
        public clsPHRPParserXTandem(string strDatasetName, string strInputFilePath, bool blnLoadModsAndSeqInfo)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.XTandem, blnLoadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserXTandem(string strDatasetName, string strInputFilePath, clsPHRPStartupOptions startupOptions)
            : base(strDatasetName, strInputFilePath, clsPHRPReader.ePeptideHitResultType.XTandem, startupOptions)
        {
        }

        private static string AppendToString(string strText, string strAppend)
        {
            if (string.IsNullOrEmpty(strText))
            {
                return strAppend;
            }

            return strText + "; " + strAppend;
        }

        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_Result_ID);
            AddHeaderColumn(DATA_COLUMN_Group_ID);
            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_Peptide_MH);
            AddHeaderColumn(DATA_COLUMN_Peptide_Hyperscore);
            AddHeaderColumn(DATA_COLUMN_Peptide_Expectation_Value_LogE);
            AddHeaderColumn(DATA_COLUMN_Multiple_Protein_Count);
            AddHeaderColumn(DATA_COLUMN_Peptide_Sequence);
            AddHeaderColumn(DATA_COLUMN_DeltaCn2);
            AddHeaderColumn(DATA_COLUMN_y_score);
            AddHeaderColumn(DATA_COLUMN_y_ions);
            AddHeaderColumn(DATA_COLUMN_b_score);
            AddHeaderColumn(DATA_COLUMN_b_ions);
            AddHeaderColumn(DATA_COLUMN_Delta_Mass);
            AddHeaderColumn(DATA_COLUMN_Peptide_Intensity_LogI);
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
            var blnPPM = false;

            double dblTolerancePlus = 0;
            double dblToleranceMinus = 0;

            if (objSearchEngineParams.Parameters.TryGetValue("spectrum, out parent monoisotopic mass error units", out var strUnits))
            {
                if (strUnits.ToLower().Trim() == "ppm")
                    blnPPM = true;
            }

            if (objSearchEngineParams.Parameters.TryGetValue("spectrum, out parent monoisotopic mass error minus", out var strTolerance))
            {
                double.TryParse(strTolerance, out dblToleranceMinus);
            }

            if (objSearchEngineParams.Parameters.TryGetValue("spectrum, out parent monoisotopic mass error plus", out strTolerance))
            {
                double.TryParse(strTolerance, out dblTolerancePlus);
            }

            var dblTolerance = Math.Max(dblToleranceMinus, dblTolerancePlus);
            if (blnPPM)
            {
                dblTolerancePPM = dblTolerance;

                // Convert from PPM to dalton (assuming a mass of 2000 m/z)
                dblTolerance = clsPeptideMassCalculator.PPMToMass(dblTolerance, 2000);
            }
            else
            {
                // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblTolerance, 2000);
            }

            return dblTolerance;
        }

        public static string GetPHRPFirstHitsFileName(string strDatasetName)
        {
            // X!Tandem does not have a first-hits file; just the _xt.txt file
            return string.Empty;
        }

        public static string GetPHRPModSummaryFileName(string strDatasetName)
        {
            return strDatasetName + "_xt_ModSummary.txt";
        }

        public static string GetPHRPPepToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_xt_PepToProtMapMTS.txt";
        }

        public static string GetPHRPProteinModsFileName(string strDatasetName)
        {
            return strDatasetName + "_xt_ProteinMods.txt";
        }

        public static string GetPHRPSynopsisFileName(string strDatasetName)
        {
            return strDatasetName + FILENAME_SUFFIX_SYN;
        }

        public static string GetPHRPResultToSeqMapFileName(string strDatasetName)
        {
            return strDatasetName + "_xt_ResultToSeqMap.txt";
        }

        public static string GetPHRPSeqInfoFileName(string strDatasetName)
        {
            return strDatasetName + "_xt_SeqInfo.txt";
        }

        public static string GetPHRPSeqToProteinMapFileName(string strDatasetName)
        {
            return strDatasetName + "_xt_SeqToProteinMap.txt";
        }

        public static List<string> GetAdditionalSearchEngineParamFileNames(string strSearchEngineParamFilePath)
        {
            var lstFileNames = new List<string>();

            var strErrorMessage = string.Empty;

            try
            {
                if (!File.Exists(strSearchEngineParamFilePath))
                {
                    lstFileNames.Add("default_input.xml  (Not Confirmed)");
                    lstFileNames.Add("taxonomy.xml  (Not Confirmed)");
                }
                else
                {
                    var paramFile = new FileInfo(strSearchEngineParamFilePath);
                    var objSearchEngineParams = new clsSearchEngineParameters(XT_SEARCH_ENGINE_NAME);

                    try
                    {
                        var strDefaultParamsFilename = GetXTandemDefaultParamsFilename(strSearchEngineParamFilePath);
                        if (!string.IsNullOrEmpty(strDefaultParamsFilename))
                        {
                            lstFileNames.Add(strDefaultParamsFilename);
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error in GetXTandemDefaultParamsFilename: " + ex.Message);
                    }

                    try
                    {
                        ParseXTandemParamFileWork(paramFile.DirectoryName, paramFile.Name, objSearchEngineParams, false, true, ref strErrorMessage);

                        if (objSearchEngineParams.Parameters.TryGetValue(TAXONOMY_INFO_KEY_NAME, out var strTaxonomyFilename))
                        {
                            lstFileNames.Add(Path.GetFileName(strTaxonomyFilename));
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error in ParseXTandemParamFileWork: " + ex.Message);
                    }

                    if (!string.IsNullOrEmpty(strErrorMessage))
                    {
                        Console.WriteLine(strErrorMessage);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception in GetAdditionalSearchEngineParamFileNames: " + ex.Message);
            }

            return lstFileNames;
        }

        private static string GetFastaFileFromTaxonomyFile(string strInputFolderPath, string strTaxononomyFilename, out string strErrorMessage)
        {
            var strFastaFile = string.Empty;

            strErrorMessage = string.Empty;

            try
            {
                var strTaxonomyFilePath = Path.Combine(strInputFolderPath, strTaxononomyFilename);
                if (!File.Exists(strTaxonomyFilePath))
                {
                    strErrorMessage = AppendToString(strErrorMessage, "Warning, taxonomy file not found: " + strTaxonomyFilePath);
                }
                else
                {
                    // Open the XML file and look for the "file" element with attribute "peptide"
                    using (var objXMLReader = new XmlTextReader(strTaxonomyFilePath))
                    {
                        while (objXMLReader.Read())
                        {
                            XMLTextReaderSkipWhitespace(objXMLReader);
                            if (objXMLReader.ReadState != ReadState.Interactive)
                                break;

                            if (objXMLReader.NodeType == XmlNodeType.Element)
                            {
                                if (objXMLReader.Name.ToLower() == "file")
                                {
                                    var strFileFormat = XMLTextReaderGetAttributeValue(objXMLReader, "format", string.Empty);

                                    if (strFileFormat == "peptide")
                                    {
                                        strFastaFile = XMLTextReaderGetAttributeValue(objXMLReader, "URL", string.Empty);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                strErrorMessage = AppendToString(strErrorMessage, "Error in GetFastaFileFromTaxonomyFile: " + ex.Message);
            }

            return strFastaFile;
        }

        public static string GetSearchEngineName()
        {
            return XT_SEARCH_ENGINE_NAME;
        }

        private static string GetXTandemDefaultParamsFilename(string strParamFilePath)
        {
            var strDefaultParamsFilename = string.Empty;

            // Open the XML file and look for the "list path, default parameters" entry
            using (var objXMLReader = new XmlTextReader(strParamFilePath))
            {
                while (MoveToNextInputParam(objXMLReader, out var kvSetting))
                {
                    if (kvSetting.Key == "list path, default parameters")
                    {
                        strDefaultParamsFilename = string.Copy(kvSetting.Value);
                        break;
                    }
                }
            }

            return strDefaultParamsFilename;
        }

        /// <summary>
        /// Parses the specified X!Tandem parameter file
        /// Note that the file specified by parameter "list path, default parameters" will also be auto-parsed (if found in folder mInputFolderPath)
        /// </summary>
        /// <param name="strSearchEngineParamFileName"></param>
        /// <param name="objSearchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams)
        {
            objSearchEngineParams = new clsSearchEngineParameters(XT_SEARCH_ENGINE_NAME, mModInfo);

            return ParseXTandemParamFile(strSearchEngineParamFileName, objSearchEngineParams, blnLookForDefaultParamsFileName: true);
        }

        public bool ParseXTandemParamFile(string strParamFileName, clsSearchEngineParameters objSearchEngineParams, bool blnLookForDefaultParamsFileName)
        {
            return ParseXTandemParamFile(strParamFileName, objSearchEngineParams, blnLookForDefaultParamsFileName, blnDetermineFastaFileNameUsingTaxonomyFile: true);
        }

        public bool ParseXTandemParamFile(string strParamFileName, clsSearchEngineParameters objSearchEngineParams, bool blnLookForDefaultParamsFileName, bool blnDetermineFastaFileNameUsingTaxonomyFile)
        {
            var strErrorMessage = string.Empty;

            var blnSuccess = false;

            try
            {
                var strParamFilePath = Path.Combine(mInputFolderPath, strParamFileName);

                if (!File.Exists(strParamFilePath))
                {
                    ReportError("X!Tandem param file not found: " + strParamFilePath);
                }
                else
                {
                    try
                    {
                        blnSuccess = ParseXTandemParamFileWork(mInputFolderPath, strParamFileName, objSearchEngineParams, blnDetermineFastaFileNameUsingTaxonomyFile, blnLookForDefaultParamsFileName, ref strErrorMessage);
                    }
                    catch (Exception ex)
                    {
                        ReportError("Error in ParseXTandemParamFileWork: " + ex.Message);
                    }

                    if (!string.IsNullOrEmpty(strErrorMessage))
                    {
                        ReportError(strErrorMessage);
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, out var dblTolerancePPM);
                    objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM;
                }
            }
            catch (Exception ex)
            {
                ReportError("Error in ParseXTandemParamFile: " + ex.Message);
            }

            ReadSearchEngineVersion(mPeptideHitResultType, objSearchEngineParams);

            return blnSuccess;
        }

        private static bool ParseXTandemParamFileWork(
            string strInputFolderPath,
            string strParamFileName,
            clsSearchEngineParameters objSearchEngineParams,
            bool blnDetermineFastaFileNameUsingTaxonomyFile,
            bool blnLookForDefaultParamsFileName,
            ref string strErrorMessage)
        {
            // Note: Do not put a Try/Catch block in this function

            var strParamFilePath = Path.Combine(strInputFolderPath, strParamFileName);

            if (blnLookForDefaultParamsFileName)
            {
                string strDefaultParamsFilename;
                try
                {
                    strDefaultParamsFilename = GetXTandemDefaultParamsFilename(strParamFilePath);
                }
                catch (Exception ex)
                {
                    strErrorMessage = AppendToString(strErrorMessage, "Error in GetXTandemDefaultParamsFilename: " + ex.Message);
                    strDefaultParamsFilename = string.Empty;
                }

                if (!string.IsNullOrEmpty(strDefaultParamsFilename))
                {
                    // Read the parameters from the default parameters file and store them in objSearchEngineParams
                    // Do this by recursively calling this function

                    // First confirm that the file exists
                    if (File.Exists(Path.Combine(strInputFolderPath, strDefaultParamsFilename)))
                    {
                        ParseXTandemParamFileWork(strInputFolderPath, strDefaultParamsFilename, objSearchEngineParams, blnDetermineFastaFileNameUsingTaxonomyFile, false, ref strErrorMessage);
                    }
                    else
                    {
                        strErrorMessage = AppendToString(strErrorMessage, "Warning, file not found: " + strDefaultParamsFilename);
                    }
                }
            }

            // Now read the parameters in strParamFilePath
            using (var objXMLReader = new XmlTextReader(strParamFilePath))
            {
                while (MoveToNextInputParam(objXMLReader, out var kvSetting))
                {
                    if (!string.IsNullOrEmpty(kvSetting.Key))
                    {
                        objSearchEngineParams.AddUpdateParameter(kvSetting.Key, kvSetting.Value);

                        switch (kvSetting.Key)
                        {
                            case TAXONOMY_INFO_KEY_NAME:

                                if (blnDetermineFastaFileNameUsingTaxonomyFile)
                                {
                                    // Open the taxonomy file to determine the fasta file used
                                    var strSetting = GetFastaFileFromTaxonomyFile(strInputFolderPath, Path.GetFileName(kvSetting.Value), out strErrorMessage);

                                    if (!string.IsNullOrEmpty(strSetting))
                                    {
                                        objSearchEngineParams.FastaFilePath = strSetting;
                                    }
                                }

                                break;
                            case "spectrum, fragment mass type":
                                objSearchEngineParams.FragmentMassType = kvSetting.Value;

                                break;
                            case "scoring, maximum missed cleavage sites":
                                if (int.TryParse(kvSetting.Value, out var intValue))
                                {
                                    objSearchEngineParams.MaxNumberInternalCleavages = intValue;
                                }

                                break;
                        }
                    }
                }
            }

            return true;
        }

        private static bool MoveToNextInputParam(XmlReader objXMLReader, out KeyValuePair<string, string> kvParameter)
        {

            while (objXMLReader.Read())
            {
                XMLTextReaderSkipWhitespace(objXMLReader);
                if (objXMLReader.ReadState != ReadState.Interactive)
                    break;

                if (objXMLReader.NodeType != XmlNodeType.Element)
                    continue;

                if (objXMLReader.Name.ToLower() != "note")
                    continue;

                var strNoteType = XMLTextReaderGetAttributeValue(objXMLReader, "type", string.Empty);

                if (strNoteType != "input") continue;

                var strParamName = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);

                if (objXMLReader.Read())
                {
                    // Read the note's inner text
                    var strValue = XMLTextReaderGetInnerText(objXMLReader);

                    kvParameter = new KeyValuePair<string, string>(strParamName, strValue);
                    return true;
                }
            }

            kvParameter = default(KeyValuePair<string, string>);
            return false;
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
                    objPSM.ResultID = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Result_ID, mColumnHeaders, 0);

                    // X!Tandem only tracks the top-ranked peptide for each spectrum
                    objPSM.ScoreRank = 1;

                    var strPeptide = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Peptide_Sequence, mColumnHeaders);

                    if (fastReadMode)
                    {
                        objPSM.SetPeptide(strPeptide, blnUpdateCleanSequence: false);
                    }
                    else
                    {
                        objPSM.SetPeptide(strPeptide, mCleavageStateCalculator);
                    }

                    objPSM.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                    // Lookup the protein name(s) using mResultIDToProteins
                    if (mResultIDToProteins.TryGetValue(objPSM.ResultID, out var lstProteinsForResultID))
                    {
                        foreach (var strProtein in lstProteinsForResultID)
                        {
                            objPSM.AddProtein(strProtein);
                        }
                    }

                    // The Peptide_MH value listed in X!Tandem files is the theoretical (computed) MH of the peptide
                    // We'll update this value below using dblMassErrorDa
                    // We'll further update this value using the ScanStatsEx data
                    var dblPeptideMH = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Peptide_MH, mColumnHeaders, 0.0);
                    objPSM.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPeptideMH, 1, 0);

                    objPSM.MassErrorDa = clsPHRPReader.LookupColumnValue(strColumns, DATA_COLUMN_Delta_Mass, mColumnHeaders);
                    if (double.TryParse(objPSM.MassErrorDa, out var dblMassErrorDa))
                    {
                        // Adjust the precursor mass
                        objPSM.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPeptideMH - dblMassErrorDa, 1, 0);
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
                    AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Hyperscore);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Expectation_Value_LogE);
                    AddScore(objPSM, strColumns, DATA_COLUMN_DeltaCn2);
                    AddScore(objPSM, strColumns, DATA_COLUMN_y_score);
                    AddScore(objPSM, strColumns, DATA_COLUMN_y_ions);
                    AddScore(objPSM, strColumns, DATA_COLUMN_b_score);
                    AddScore(objPSM, strColumns, DATA_COLUMN_b_ions);
                    AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Intensity_LogI);

                    // This is the base-10 log of the expectation value
                    if (double.TryParse(objPSM.GetScore(DATA_COLUMN_Peptide_Expectation_Value_LogE), out var dblLogEValue))
                    {
                        // Record the original E-value
                        objPSM.SetScore("Peptide_Expectation_Value", Math.Pow(10, dblLogEValue).ToString("0.00e+000"));
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error parsing line " + intLinesRead + " in the X!Tandem data file", ex);
            }

            return blnSuccess;
        }

        private static string XMLTextReaderGetAttributeValue(XmlReader objXMLReader, string strAttributeName, string strValueIfMissing)
        {
            objXMLReader.MoveToAttribute(strAttributeName);
            if (objXMLReader.ReadAttributeValue())
            {
                return objXMLReader.Value;
            }

            return string.Copy(strValueIfMissing);
        }

        private static string XMLTextReaderGetInnerText(XmlReader objXMLReader)
        {
            var strValue = string.Empty;
            bool blnSuccess;

            if (objXMLReader.NodeType == XmlNodeType.Element)
            {
                // Advance the reader so that we can read the value
                blnSuccess = objXMLReader.Read();
            }
            else
            {
                blnSuccess = true;
            }

            if (blnSuccess && objXMLReader.NodeType != XmlNodeType.Whitespace & objXMLReader.HasValue)
            {
                strValue = objXMLReader.Value;
            }

            return strValue;
        }

        private static void XMLTextReaderSkipWhitespace(XmlReader objXMLReader)
        {
            if (objXMLReader.NodeType == XmlNodeType.Whitespace)
            {
                // Whitespace; read the next node
                objXMLReader.Read();
            }
        }
    }
}
