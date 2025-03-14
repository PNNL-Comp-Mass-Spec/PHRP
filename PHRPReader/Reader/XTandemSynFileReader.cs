﻿//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from X!Tandem _xt.txt files
//
// Note: in order to fully extract the search parameters, you need to have these files in the same directory as the _xt.txt file
//	The file passed to LoadSearchEngineParameters()  (typically input.xml)
//	The taxonomy file that it references             (typically taxonomy.xml)
//  The default input file defined in input.xml      (typically default_input.xml)
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using PHRPReader.Data;
using PRISM;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for X!Tandem
    /// </summary>
    public class XTandemSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: Da, Fasta, FHT, Hyperscore, ProtMap, psm, xt

        /// <summary>
        /// X!Tandem synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_xt.txt";

        /// <summary>
        /// X!Tandem first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for MaxQuant results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_xt.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string XT_SEARCH_ENGINE_NAME = "X! Tandem";

        private const string TAXONOMY_INFO_KEY_NAME = "list path, taxonomy information";

        /// <summary>
        /// Mapping from enum to synopsis file column name for X!Tandem
        /// </summary>
        private static readonly Dictionary<XTandemSynFileColumns, string> mSynopsisFileColumn = new();

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

        /// <summary>
        /// Constructor; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        public XTandemSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public XTandemSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.XTandem, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public XTandemSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.XTandem, startupOptions)
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
        /// Determines the precursor mass tolerance
        /// </summary>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <returns>Precursor tolerance, in Da</returns>
        private double DeterminePrecursorMassTolerance(SearchEngineParameters searchEngineParams, out double tolerancePPM)
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

                // Convert from PPM to Dalton (assuming a mass of 2000 m/z)
                tolerance = PeptideMassCalculator.PPMToMass(tolerance, 2000);
            }
            else
            {
                // Convert from Dalton to PPM (assuming a mass of 2000 m/z)
                tolerancePPM = PeptideMassCalculator.MassToPPM(tolerance, 2000);
            }

            return tolerance;
        }

        /// <summary>
        /// Get the header names in the PHRP synopsis or first hits file for this tool
        /// </summary>
        /// <returns>List of header names</returns>
        protected override List<string> GetColumnHeaderNames()
        {
            var headerNames = new List<string>();
            headerNames.AddRange(GetColumnHeaderNamesAndIDs().Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, XTandemSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, XTandemSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "Result_ID", XTandemSynFileColumns.ResultID },
                { "Group_ID", XTandemSynFileColumns.GroupID },
                { "Scan", XTandemSynFileColumns.Scan },
                { "Charge", XTandemSynFileColumns.Charge },
                { "Peptide_MH", XTandemSynFileColumns.MH },
                { "Peptide_Hyperscore", XTandemSynFileColumns.Hyperscore },
                { "Peptide_Expectation_Value_Log(e)", XTandemSynFileColumns.EValue },
                { "Multiple_Protein_Count", XTandemSynFileColumns.ProteinCount },
                { "Peptide_Sequence", XTandemSynFileColumns.Peptide },
                { "DeltaCn2", XTandemSynFileColumns.DeltaCn2 },
                { "y_score", XTandemSynFileColumns.YScore },
                { "y_ions", XTandemSynFileColumns.YIons },
                { "b_score", XTandemSynFileColumns.BScore },
                { "b_ions", XTandemSynFileColumns.BIons },
                { "Delta_Mass", XTandemSynFileColumns.DelM },
                { "Peptide_Intensity_Log(I)", XTandemSynFileColumns.Intensity },
                { "DelM_PPM", XTandemSynFileColumns.DelMPPM },
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum XTandemSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames">List of header names</param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<XTandemSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column">Column enum</param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(XTandemSynFileColumns column)
        {
            if (mSynopsisFileColumn.Count > 0)
            {
                return mSynopsisFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mSynopsisFileColumn.Add(item.Value, item.Key);
            }

            return mSynopsisFileColumn[column];
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Empty string, since X!Tandem does not have a first-hits file; just the _xt.txt file</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_xt_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_xt_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_xt_ProteinMods.txt";
        }

        /// <summary>
        /// Default Synopsis file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_SYN;
        }

        /// <summary>
        /// Default ResultToSeq map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(string datasetName)
        {
            return datasetName + "_xt_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_xt_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_xt_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Additional search engine parameter file names
        /// </summary>
        /// <param name="searchEngineParamFilePath">Search engine parameter file path</param>
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
                    var searchEngineParams = new SearchEngineParameters(XT_SEARCH_ENGINE_NAME);

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
            errorMessage = string.Empty;

            try
            {
                var taxonomyFilePath = Path.Combine(inputDirectoryPath, taxonomyFilename);

                if (!File.Exists(taxonomyFilePath))
                {
                    errorMessage = AppendToString(errorMessage, "Warning, taxonomy file not found: " + taxonomyFilePath);
                    return string.Empty;
                }

                // Open the XML file and look for the "file" element with attribute "peptide"
                using var xmlReader = new XmlTextReader(taxonomyFilePath);

                while (xmlReader.Read())
                {
                    XMLTextReaderSkipWhitespace(xmlReader);
                    if (xmlReader.ReadState != ReadState.Interactive)
                        break;

                    if (xmlReader.NodeType == XmlNodeType.Element)
                    {
                        if (string.Equals(xmlReader.Name, "file", StringComparison.OrdinalIgnoreCase))
                        {
                            var fileFormat = XMLTextReaderGetAttributeValue(xmlReader, "format", string.Empty);

                            if (fileFormat == "peptide")
                            {
                                return XMLTextReaderGetAttributeValue(xmlReader, "URL", string.Empty);
                            }
                        }
                    }
                }

                return string.Empty;
            }
            catch (Exception ex)
            {
                errorMessage = AppendToString(errorMessage, "Error in GetFastaFileFromTaxonomyFile: " + ex.Message);
                return string.Empty;
            }
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
            using var xmlReader = new XmlTextReader(paramFilePath);

            while (MoveToNextInputParam(xmlReader, out var kvSetting))
            {
                if (kvSetting.Key == "list path, default parameters")
                {
                    defaultParamsFilename = kvSetting.Value;
                    break;
                }
            }

            return defaultParamsFilename;
        }

        /// <summary>
        /// Parses the specified X!Tandem parameter file
        /// Note that the file specified by parameter "list path, default parameters" will also be auto-parsed (if found in directory mInputDirectoryPath)
        /// </summary>
        /// <param name="searchEngineParamFileName">X!Tandem parameter file name</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(XT_SEARCH_ENGINE_NAME, mModInfo);

            return ParseXTandemParamFile(searchEngineParamFileName, searchEngineParams, lookForDefaultParamsFileName: true);
        }

        /// <summary>
        /// Parse an X!Tandem parameter file
        /// </summary>
        /// <param name="paramFileName">X!Tandem parameter file name</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <param name="lookForDefaultParamsFileName">When true, look for the default X!Tandem parameter file</param>
        /// <param name="determineFastaFileNameUsingTaxonomyFile">When true, determine the FASTA file name using the taxonomy file</param>
        /// <returns>True if successful, false if an error</returns>
        public bool ParseXTandemParamFile(
            string paramFileName,
            SearchEngineParameters searchEngineParams,
            bool lookForDefaultParamsFileName,
            bool determineFastaFileNameUsingTaxonomyFile = true)
        {
            var errorMessage = string.Empty;

            var success = false;

            try
            {
                var paramFilePath = Path.Combine(InputDirectoryPath, paramFileName);

                if (!File.Exists(paramFilePath))
                {
                    ReportError("X!Tandem param file not found: " + paramFilePath);
                }
                else
                {
                    try
                    {
                        success = ParseXTandemParamFileWork(InputDirectoryPath, paramFileName, searchEngineParams, determineFastaFileNameUsingTaxonomyFile, lookForDefaultParamsFileName, ref errorMessage);
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
            SearchEngineParameters searchEngineParams,
            bool determineFastaFileNameUsingTaxonomyFile,
            bool lookForDefaultParamsFileName,
            ref string errorMessage)
        {
            // Note: Do not put a Try/Catch block in this method

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
                    // Do this by recursively calling this method

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
            using var xmlReader = new XmlTextReader(paramFilePath);

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
                                // Open the taxonomy file to determine the FASTA file used
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

                if (!string.Equals(xmlReader.Name, "note", StringComparison.OrdinalIgnoreCase))
                    continue;

                var noteType = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty);

                if (noteType != "input")
                    continue;

                var paramName = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);

                if (xmlReader.Read())
                {
                    // Read the note's inner text
                    var value = XMLTextReaderGetInnerText(xmlReader);

                    kvParameter = new KeyValuePair<string, string>(paramName, value);
                    return true;
                }
            }

            kvParameter = default;
            return false;
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM details</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;

            var columns = line.Split('\t');

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.ResultID), mColumnHeaders, 0);

                // X!Tandem only tracks the top-ranked peptide for each spectrum
                psm.ScoreRank = 1;

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.Charge), mColumnHeaders, 0);

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
                var peptideMH = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.MH), mColumnHeaders, 0.0);
                psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(peptideMH, 1, 0);

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.DelM), mColumnHeaders);
                if (double.TryParse(psm.MassErrorDa, out var massErrorDa))
                {
                    // Adjust the precursor mass
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(peptideMH - massErrorDa, 1, 0);
                }

                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(XTandemSynFileColumns.DelMPPM), mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                // Store the remaining scores
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.Hyperscore));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.EValue));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.DeltaCn2));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.YScore));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.YIons));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.BScore));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.BIons));
                AddScore(psm, columns, GetColumnNameByID(XTandemSynFileColumns.Intensity));

                // This is the base-10 log of the expectation value
                if (double.TryParse(psm.GetScore(GetColumnNameByID(XTandemSynFileColumns.EValue)), out var logEValue))
                {
                    // Record the original E-value
                    psm.SetScore("Peptide_Expectation_Value", Math.Pow(10, logEValue).ToString("0.00e+000"));
                }

                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error parsing line " + linesRead + " in the X!Tandem data file", ex);
                return false;
            }
        }

        private static string XMLTextReaderGetAttributeValue(XmlReader xmlReader, string attributeName, string valueIfMissing)
        {
            xmlReader.MoveToAttribute(attributeName);
            if (xmlReader.ReadAttributeValue())
            {
                return xmlReader.Value;
            }

            return valueIfMissing;
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

            if (success && xmlReader.NodeType != XmlNodeType.Whitespace && xmlReader.HasValue)
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
