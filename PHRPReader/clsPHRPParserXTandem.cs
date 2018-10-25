//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from Xtandem _xt.txt files
//
// Note: in order to fully extract the search parameters, you need to have these files in the same directory as the _xt.txt file
//	The file passed to LoadSearchEngineParameters()  (typically input.xml)
//	The taxonomy file that it references             (typically taxonomy.xml)
//   The default input file defined in input.xml      (typically default_input.xml)
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using PRISM;

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for X!Tandem
    /// </summary>
    public class clsPHRPParserXTandem : clsPHRPParser
    {
        #region "Constants"

#pragma warning disable 1591
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
        public clsPHRPParserXTandem(string datasetName, string inputFilePath)
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
        public clsPHRPParserXTandem(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.XTandem, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserXTandem(string datasetName, string inputFilePath, clsPHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.XTandem, startupOptions)
        {
        }

        private static string AppendToString(string text, string append)
        {
            if (string.IsNullOrEmpty(text))
            {
                return append;
            }

            return text + "; " + append;
        }

        /// <summary>
        /// Define column header names
        /// </summary>
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
        /// <param name="searchEngineParams"></param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <returns>Precursor tolerance, in Da</returns>
        /// <remarks></remarks>
        private double DeterminePrecursorMassTolerance(clsSearchEngineParameters searchEngineParams, out double tolerancePPM)
        {
            var isPPM = false;

            double tolerancePlus = 0;
            double toleranceMinus = 0;

            if (searchEngineParams.Parameters.TryGetValue("spectrum, out parent monoisotopic mass error units", out var units))
            {
                if (units.ToLower().Trim() == "ppm")
                    isPPM = true;
            }

            if (searchEngineParams.Parameters.TryGetValue("spectrum, out parent monoisotopic mass error minus", out var toleranceText))
            {
                double.TryParse(toleranceText, out toleranceMinus);
            }

            if (searchEngineParams.Parameters.TryGetValue("spectrum, out parent monoisotopic mass error plus", out toleranceText))
            {
                double.TryParse(toleranceText, out tolerancePlus);
            }

            var tolerance = Math.Max(toleranceMinus, tolerancePlus);
            if (isPPM)
            {
                tolerancePPM = tolerance;

                // Convert from PPM to dalton (assuming a mass of 2000 m/z)
                tolerance = clsPeptideMassCalculator.PPMToMass(tolerance, 2000);
            }
            else
            {
                // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                tolerancePPM = clsPeptideMassCalculator.MassToPPM(tolerance, 2000);
            }

            return tolerance;
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // X!Tandem does not have a first-hits file; just the _xt.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_xt_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_xt_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_xt_ProteinMods.txt";
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
            return datasetName + "_xt_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_xt_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_xt_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Additional search engine parameter file names
        /// </summary>
        /// <param name="searchEngineParamFilePath"></param>
        /// <returns></returns>
        public static List<string> GetAdditionalSearchEngineParamFileNames(string searchEngineParamFilePath)
        {
            var fileNames = new List<string>();

            var errorMessage = string.Empty;

            try
            {
                if (!File.Exists(searchEngineParamFilePath))
                {
                    fileNames.Add("default_input.xml  (Not Confirmed)");
                    fileNames.Add("taxonomy.xml  (Not Confirmed)");
                }
                else
                {
                    var paramFile = new FileInfo(searchEngineParamFilePath);
                    var searchEngineParams = new clsSearchEngineParameters(XT_SEARCH_ENGINE_NAME);

                    try
                    {
                        var defaultParamsFilename = GetXTandemDefaultParamsFilename(searchEngineParamFilePath);
                        if (!string.IsNullOrEmpty(defaultParamsFilename))
                        {
                            fileNames.Add(defaultParamsFilename);
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error in GetXTandemDefaultParamsFilename: " + ex.Message);
                    }

                    try
                    {
                        if (paramFile.Directory == null)
                        {
                            ConsoleMsgUtils.ShowWarning("Unable to determine the parent directory of " + paramFile.FullName);
                        }
                        else
                        {

                            ParseXTandemParamFileWork(paramFile.Directory.FullName, paramFile.Name, searchEngineParams, false, true,
                                                      ref errorMessage);

                            if (searchEngineParams.Parameters.TryGetValue(TAXONOMY_INFO_KEY_NAME, out var taxonomyFilename))
                            {
                                fileNames.Add(Path.GetFileName(taxonomyFilename));
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("Error in ParseXTandemParamFileWork: " + ex.Message);
                    }

                    if (!string.IsNullOrEmpty(errorMessage))
                    {
                        Console.WriteLine(errorMessage);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Exception in GetAdditionalSearchEngineParamFileNames: " + ex.Message);
            }

            return fileNames;
        }

        private static string GetFastaFileFromTaxonomyFile(string inputDirectoryPath, string taxonomyFilename, out string errorMessage)
        {
            var fastaFile = string.Empty;

            errorMessage = string.Empty;

            try
            {
                var taxonomyFilePath = Path.Combine(inputDirectoryPath, taxonomyFilename);
                if (!File.Exists(taxonomyFilePath))
                {
                    errorMessage = AppendToString(errorMessage, "Warning, taxonomy file not found: " + taxonomyFilePath);
                }
                else
                {
                    // Open the XML file and look for the "file" element with attribute "peptide"
                    using (var xmlReader = new XmlTextReader(taxonomyFilePath))
                    {
                        while (xmlReader.Read())
                        {
                            XMLTextReaderSkipWhitespace(xmlReader);
                            if (xmlReader.ReadState != ReadState.Interactive)
                                break;

                            if (xmlReader.NodeType == XmlNodeType.Element)
                            {
                                if (xmlReader.Name.ToLower() == "file")
                                {
                                    var fileFormat = XMLTextReaderGetAttributeValue(xmlReader, "format", string.Empty);

                                    if (fileFormat == "peptide")
                                    {
                                        fastaFile = XMLTextReaderGetAttributeValue(xmlReader, "URL", string.Empty);
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
                errorMessage = AppendToString(errorMessage, "Error in GetFastaFileFromTaxonomyFile: " + ex.Message);
            }

            return fastaFile;
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return XT_SEARCH_ENGINE_NAME;
        }

        private static string GetXTandemDefaultParamsFilename(string paramFilePath)
        {
            var defaultParamsFilename = string.Empty;

            // Open the XML file and look for the "list path, default parameters" entry
            using (var xmlReader = new XmlTextReader(paramFilePath))
            {
                while (MoveToNextInputParam(xmlReader, out var kvSetting))
                {
                    if (kvSetting.Key == "list path, default parameters")
                    {
                        defaultParamsFilename = string.Copy(kvSetting.Value);
                        break;
                    }
                }
            }

            return defaultParamsFilename;
        }

        /// <summary>
        /// Parses the specified X!Tandem parameter file
        /// Note that the file specified by parameter "list path, default parameters" will also be auto-parsed (if found in directory mInputDirectoryPath)
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out clsSearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new clsSearchEngineParameters(XT_SEARCH_ENGINE_NAME, mModInfo);

            return ParseXTandemParamFile(searchEngineParamFileName, searchEngineParams, lookForDefaultParamsFileName: true);
        }

        /// <summary>
        /// Parse an X!Tandem parameter file
        /// </summary>
        /// <param name="paramFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <param name="lookForDefaultParamsFileName"></param>
        /// <param name="determineFastaFileNameUsingTaxonomyFile"></param>
        /// <returns></returns>
        public bool ParseXTandemParamFile(
            string paramFileName,
            clsSearchEngineParameters searchEngineParams,
            bool lookForDefaultParamsFileName,
            bool determineFastaFileNameUsingTaxonomyFile = true)
        {
            var errorMessage = string.Empty;

            var success = false;

            try
            {
                var paramFilePath = Path.Combine(mInputDirectoryPath, paramFileName);

                if (!File.Exists(paramFilePath))
                {
                    ReportError("X!Tandem param file not found: " + paramFilePath);
                }
                else
                {
                    try
                    {
                        success = ParseXTandemParamFileWork(mInputDirectoryPath, paramFileName, searchEngineParams, determineFastaFileNameUsingTaxonomyFile, lookForDefaultParamsFileName, ref errorMessage);
                    }
                    catch (Exception ex)
                    {
                        ReportError("Error in ParseXTandemParamFileWork: " + ex.Message);
                    }

                    if (!string.IsNullOrEmpty(errorMessage))
                    {
                        ReportError(errorMessage);
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM);
                    searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;
                }
            }
            catch (Exception ex)
            {
                ReportError("Error in ParseXTandemParamFile: " + ex.Message);
            }

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private static bool ParseXTandemParamFileWork(
            string inputDirectoryPath,
            string paramFileName,
            clsSearchEngineParameters searchEngineParams,
            bool determineFastaFileNameUsingTaxonomyFile,
            bool lookForDefaultParamsFileName,
            ref string errorMessage)
        {
            // Note: Do not put a Try/Catch block in this function

            var paramFilePath = Path.Combine(inputDirectoryPath, paramFileName);

            if (lookForDefaultParamsFileName)
            {
                string defaultParamsFilename;
                try
                {
                    defaultParamsFilename = GetXTandemDefaultParamsFilename(paramFilePath);
                }
                catch (Exception ex)
                {
                    errorMessage = AppendToString(errorMessage, "Error in GetXTandemDefaultParamsFilename: " + ex.Message);
                    defaultParamsFilename = string.Empty;
                }

                if (!string.IsNullOrEmpty(defaultParamsFilename))
                {
                    // Read the parameters from the default parameters file and store them in searchEngineParams
                    // Do this by recursively calling this function

                    // First confirm that the file exists
                    if (File.Exists(Path.Combine(inputDirectoryPath, defaultParamsFilename)))
                    {
                        ParseXTandemParamFileWork(inputDirectoryPath, defaultParamsFilename, searchEngineParams, determineFastaFileNameUsingTaxonomyFile, false, ref errorMessage);
                    }
                    else
                    {
                        errorMessage = AppendToString(errorMessage, "Warning, file not found: " + defaultParamsFilename);
                    }
                }
            }

            // Now read the parameters in paramFilePath
            using (var xmlReader = new XmlTextReader(paramFilePath))
            {
                while (MoveToNextInputParam(xmlReader, out var kvSetting))
                {
                    if (!string.IsNullOrEmpty(kvSetting.Key))
                    {
                        searchEngineParams.AddUpdateParameter(kvSetting.Key, kvSetting.Value);

                        switch (kvSetting.Key)
                        {
                            case TAXONOMY_INFO_KEY_NAME:

                                if (determineFastaFileNameUsingTaxonomyFile)
                                {
                                    // Open the taxonomy file to determine the fasta file used
                                    var setting = GetFastaFileFromTaxonomyFile(inputDirectoryPath, Path.GetFileName(kvSetting.Value), out errorMessage);

                                    if (!string.IsNullOrEmpty(setting))
                                    {
                                        searchEngineParams.FastaFilePath = setting;
                                    }
                                }

                                break;
                            case "spectrum, fragment mass type":
                                searchEngineParams.FragmentMassType = kvSetting.Value;

                                break;
                            case "scoring, maximum missed cleavage sites":
                                if (int.TryParse(kvSetting.Value, out var value))
                                {
                                    searchEngineParams.MaxNumberInternalCleavages = value;
                                }

                                break;
                        }
                    }
                }
            }

            return true;
        }

        private static bool MoveToNextInputParam(XmlReader xmlReader, out KeyValuePair<string, string> kvParameter)
        {

            while (xmlReader.Read())
            {
                XMLTextReaderSkipWhitespace(xmlReader);
                if (xmlReader.ReadState != ReadState.Interactive)
                    break;

                if (xmlReader.NodeType != XmlNodeType.Element)
                    continue;

                if (xmlReader.Name.ToLower() != "note")
                    continue;

                var noteType = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty);

                if (noteType != "input") continue;

                var paramName = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);

                if (xmlReader.Read())
                {
                    // Read the note's inner text
                    var value = XMLTextReaderGetInnerText(xmlReader);

                    kvParameter = new KeyValuePair<string, string>(paramName, value);
                    return true;
                }
            }

            kvParameter = default(KeyValuePair<string, string>);
            return false;
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
            var columns = line.Split('\t');

            var success = false;

            psm = new clsPSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, -100);
                if (psm.ScanNumber == -100)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Result_ID, mColumnHeaders, 0);

                    // X!Tandem only tracks the top-ranked peptide for each spectrum
                    psm.ScoreRank = 1;

                    var peptide = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Peptide_Sequence, mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                    // Lookup the protein name(s) using mResultIDToProteins
                    if (mResultIDToProteins.TryGetValue(psm.ResultID, out var proteinsForResultID))
                    {
                        foreach (var protein in proteinsForResultID)
                        {
                            psm.AddProtein(protein);
                        }
                    }

                    // The Peptide_MH value listed in X!Tandem files is the theoretical (computed) MH of the peptide
                    // We'll update this value below using massErrorDa
                    // We'll further update this value using the ScanStatsEx data
                    var peptideMH = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Peptide_MH, mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(peptideMH, 1, 0);

                    psm.MassErrorDa = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Delta_Mass, mColumnHeaders);
                    if (double.TryParse(psm.MassErrorDa, out var massErrorDa))
                    {
                        // Adjust the precursor mass
                        psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(peptideMH - massErrorDa, 1, 0);
                    }

                    psm.MassErrorPPM = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, DATA_COLUMN_Peptide_Hyperscore);
                    AddScore(psm, columns, DATA_COLUMN_Peptide_Expectation_Value_LogE);
                    AddScore(psm, columns, DATA_COLUMN_DeltaCn2);
                    AddScore(psm, columns, DATA_COLUMN_y_score);
                    AddScore(psm, columns, DATA_COLUMN_y_ions);
                    AddScore(psm, columns, DATA_COLUMN_b_score);
                    AddScore(psm, columns, DATA_COLUMN_b_ions);
                    AddScore(psm, columns, DATA_COLUMN_Peptide_Intensity_LogI);

                    // This is the base-10 log of the expectation value
                    if (double.TryParse(psm.GetScore(DATA_COLUMN_Peptide_Expectation_Value_LogE), out var logEValue))
                    {
                        // Record the original E-value
                        psm.SetScore("Peptide_Expectation_Value", Math.Pow(10, logEValue).ToString("0.00e+000"));
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error parsing line " + linesRead + " in the X!Tandem data file", ex);
            }

            return success;
        }

        private static string XMLTextReaderGetAttributeValue(XmlReader xmlReader, string attributeName, string valueIfMissing)
        {
            xmlReader.MoveToAttribute(attributeName);
            if (xmlReader.ReadAttributeValue())
            {
                return xmlReader.Value;
            }

            return string.Copy(valueIfMissing);
        }

        private static string XMLTextReaderGetInnerText(XmlReader xmlReader)
        {
            var value = string.Empty;
            bool success;

            if (xmlReader.NodeType == XmlNodeType.Element)
            {
                // Advance the reader so that we can read the value
                success = xmlReader.Read();
            }
            else
            {
                success = true;
            }

            if (success && xmlReader.NodeType != XmlNodeType.Whitespace & xmlReader.HasValue)
            {
                value = xmlReader.Value;
            }

            return value;
        }

        private static void XMLTextReaderSkipWhitespace(XmlReader xmlReader)
        {
            if (xmlReader.NodeType == XmlNodeType.Whitespace)
            {
                // Whitespace; read the next node
                xmlReader.Read();
            }
        }
    }
}
