using System;
using System.Collections.Generic;

namespace PHRPReader
{
    // ReSharper disable UnusedMember.Global
#pragma warning disable RCS1163 // Unused parameter.
#pragma warning disable IDE0060 // Remove unused parameter

    /// <summary>
    /// Legacy PHRP Reader
    /// </summary>
    [Obsolete("Legacy PHRP Reader; superseded by the ReaderFactory class")]
    public class clsPHRPReader : PRISM.EventNotifier, IDisposable
    {
        // Ignore Spelling: xt, msx, fht, Ss, Za, msgfdb, MODa, moda, modp, kv, toppic, mspath, msa, modplus, msp, prot, tpc

        #region "Constants"

        /// <summary>
        /// Symbol used by PHRP to indicate a protein terminus
        /// </summary>
        public const char PROTEIN_TERMINUS_SYMBOL_PHRP = '-';

#pragma warning disable 1591

        public const string MSGF_RESULT_COLUMN_SpectrumFile = "#SpectrumFile";
        public const string MSGF_RESULT_COLUMN_Title = "Title";
        public const string MSGF_RESULT_COLUMN_Annotation = "Annotation";

        public const string XT_RESULT_TO_SEQ_MAP_SUFFIX = "_xt_ResultToSeqMap.txt";
        public const string XT_SEQ_TO_PROTEIN_MAP_SUFFIX = "_xt_SeqToProteinMap.txt";

        public const string DOT_RAW_EXTENSION = ".raw";
        public const string DOT_MZXML_EXTENSION = ".mzXML";

        public const string MSGF_RESULT_FILENAME_SUFFIX = "_MSGF.txt";
        public const string SCAN_STATS_FILENAME_SUFFIX = "_ScanStats.txt";
        public const string EXTENDED_SCAN_STATS_FILENAME_SUFFIX = "_ScanStatsEx.txt";

#pragma warning restore 1591

        #endregion

        #region "Properties"

        /// <summary>
        /// Returns True if the input file was successfully opened and data remains to be read
        /// </summary>
        /// <returns>True if the file is readable</returns>
        public bool CanRead { get; } = false;

        // public clsPSM CurrentPSM => mPSMCurrent;

        // public clsSeqInfo CurrentPSMSeqInfo

        /// <summary>
        /// Dataset name (auto-determined based on the input filename)
        /// </summary>
        public string DatasetName { get; } = string.Empty;

        /// <summary>
        /// If True, will display messages at the console
        /// </summary>
        public bool EchoMessagesToConsole { get; set; }

        /// <summary>
        /// Cached error messages
        /// </summary>
        public List<string> ErrorMessages { get; }

        /// <summary>
        /// Current error message
        /// </summary>
        public string ErrorMessage { get; } = string.Empty;

        /// <summary>
        /// Used to enable fast read mode when calling MoveNext
        /// When FastReadMode is True, you should call FinalizeCurrentPSM after calling MoveNext to populate the remaining fields if the peptide is a peptide of interest
        /// </summary>
        /// <remarks>Once FastReadMode is enabled it cannot be turned off (this is a safety measure due to how data is cached)</remarks>
        public bool FastReadMode { get; set; }

        /// <summary>
        /// If True, looks for and loads the modification definitions from the _ModSummary.txt file associated with the input file
        /// Also reads the SeqInfo and related files
        /// </summary>
        public bool LoadModsAndSeqInfo => false;

        /// <summary>
        /// If true, loads the MSGF SpecProb values from the _MSGF.txt file associated with the input file
        /// </summary>
        public bool LoadMSGFResults => false;

        /// <summary>
        /// If True, loads the MASIC _ScanStats.txt file
        /// </summary>
        public bool LoadScanStatsData => false;

        /// <summary>
        /// The maximum number of proteins that will be tracked for each PSM
        /// </summary>
        public int MaxProteinsPerPSM => 0;

        /// <summary>
        /// Returns True if the ModSummary file was successfully loaded
        /// </summary>
        public bool ModSummaryFileLoaded { get; } = false;

        /// <summary>
        /// Peptide hit result type; SEQUEST, XTandem, Inspect, MSGFPlus, etc.
        /// </summary>
        public PeptideHitResultTypes PeptideHitResultType => PeptideHitResultTypes.Unknown;

        /// <summary>
        /// Returns a number between 0 and 100 indicating the percentage of the source file that has been read
        /// </summary>
        public float PercentComplete => 0;

        // public clsPHRPParser PHRPParser { get; private set; }

        /// <summary>
        /// Returns the cached mapping between ResultID and SeqID
        /// </summary>
        public SortedList<int, int> ResultToSeqMap => null;

        // public SortedList<int, clsSeqInfo> SeqInfo

        // public SortedList<int, List<clsProteinInfo>> SeqToProteinMap

        /// <summary>
        /// When True, skips near-duplicate lines in the PHRP data file (lines with the same peptide in the same scan, but different protein names)
        /// </summary>
        public bool SkipDuplicatePSMs { get; set; }

        /// <summary>
        /// Cached warning messages
        /// </summary>
        public List<string> WarningMessages { get; }

        #endregion

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
        public clsPHRPReader(string inputFilePath)
            : this(inputFilePath, PeptideHitResultTypes.Unknown, loadModsAndSeqInfo: true, loadMSGFResults: true, loadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="resultType">Source file PeptideHit result type</param>
        /// <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
        public clsPHRPReader(string inputFilePath, PeptideHitResultTypes resultType)
            : this(inputFilePath, resultType, loadModsAndSeqInfo: true, loadMSGFResults: true, loadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="loadModsAndSeqInfo">If True, looks for and auto-loads the modification definitions from the _ModSummary.txt file</param>
        /// <param name="loadMSGFResults">If True, looks for and auto-loads the MSGF results from the _msg.txt file</param>
        public clsPHRPReader(string inputFilePath, bool loadModsAndSeqInfo, bool loadMSGFResults)
            : this(inputFilePath, PeptideHitResultTypes.Unknown, loadModsAndSeqInfo, loadMSGFResults, loadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="loadModsAndSeqInfo">If True, looks for and auto-loads the modification definitions from the _ModSummary.txt file</param>
        /// <param name="loadMSGFResults">If True, looks for and auto-loads the MSGF results from the _msg.txt file</param>
        /// <param name="loadScanStats">If True, looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
        public clsPHRPReader(string inputFilePath, bool loadModsAndSeqInfo, bool loadMSGFResults, bool loadScanStats)
            : this(inputFilePath, PeptideHitResultTypes.Unknown, loadModsAndSeqInfo, loadMSGFResults, loadScanStats)
        {
        }

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="startupOptions">Startup options</param>
        public clsPHRPReader(string inputFilePath, clsPHRPStartupOptions startupOptions)
            : this(inputFilePath, PeptideHitResultTypes.Unknown, startupOptions)
        {
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="resultType">Source file PeptideHit result type</param>
        /// <param name="loadModsAndSeqInfo">If True, looks for and auto-loads the modification definitions from the _ModSummary.txt file</param>
        /// <param name="loadMSGFResults">If True, looks for and auto-loads the MSGF results from the _msg.txt file</param>
        public clsPHRPReader(string inputFilePath, PeptideHitResultTypes resultType, bool loadModsAndSeqInfo, bool loadMSGFResults)
            : this(inputFilePath, resultType, loadModsAndSeqInfo, loadMSGFResults, loadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="resultType">Source file PeptideHit result type</param>
        /// <param name="loadModsAndSeqInfo">If True, looks for and auto-loads the modification definitions from the _ModSummary.txt file</param>
        /// <param name="loadMSGFResults">If True, looks for and auto-loads the MSGF results from the _msg.txt file</param>
        /// <param name="loadScanStats">If True, looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
        public clsPHRPReader(string inputFilePath, PeptideHitResultTypes resultType, bool loadModsAndSeqInfo, bool loadMSGFResults,
                             bool loadScanStats)
        {
            ErrorMessages = new List<string>();
            WarningMessages = new List<string>();

            throw new Exception("Class clsPHRPReader is obsolete; replace with the ReaderFactory class");
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="resultType">Source file PeptideHit result type</param>
        /// <param name="startupOptions">Startup options</param>
        public clsPHRPReader(string inputFilePath, PeptideHitResultTypes resultType, clsPHRPStartupOptions startupOptions)
        {
            ErrorMessages = new List<string>();
            WarningMessages = new List<string>();

            throw new Exception("Class clsPHRPReader is obsolete; replace with the ReaderFactory class");
        }

        /// <summary>
        /// Updates filePath to have _msgfdb instead of _msgfplus if basePHRPFileName contains _msgfdb
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="basePHRPFileName"></param>
        [Obsolete("Superseded by ReaderFactory.AutoSwitchToLegacyMSGFDBIfRequired")]
        public static string AutoSwitchToLegacyMSGFDBIfRequired(string filePath, string basePHRPFileName)
        {
            return ReaderFactory.AutoSwitchToLegacyMSGFDBIfRequired(filePath, basePHRPFileName);
        }

        /// <summary>
        /// Clear any cached error messages
        /// </summary>
        public void ClearErrors()
        {
            ErrorMessages.Clear();
        }

        /// <summary>
        /// Clear any cached warning messages
        /// </summary>
        public void ClearWarnings()
        {
            WarningMessages.Clear();
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineBestInputFile")]
        public static string AutoDetermineBestInputFile(string inputDirectoryPath)
        {
            return AutoDetermineBestInputFile(inputDirectoryPath, out _);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="matchedResultType">Output parameter: the result type of the best result file found</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineBestInputFile")]
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, out PeptideHitResultTypes matchedResultType)
        {
            return ReaderFactory.AutoDetermineBestInputFile(inputDirectoryPath, out matchedResultType);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineBestInputFile")]
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, string datasetName)
        {
            return AutoDetermineBestInputFile(inputDirectoryPath, datasetName, out _);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="matchedResultType">Output parameter: the result type of the best result file found</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineBestInputFile")]
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, string datasetName,
                                                        out PeptideHitResultTypes matchedResultType)
        {
            var datasetNames = new List<string>
            {
                datasetName
            };

            return AutoDetermineBestInputFile(inputDirectoryPath, datasetNames, out matchedResultType);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the given list of datasets in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="datasetNames">List of dataset names to search for</param>
        /// <param name="matchedResultType">Output parameter: the result type of the best result file found</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineBestInputFile")]
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, List<string> datasetNames,
                                                        out PeptideHitResultTypes matchedResultType)
        {
            return ReaderFactory.AutoDetermineBestInputFile(inputDirectoryPath, datasetNames, out matchedResultType);
        }

        /// <summary>
        /// Auto-determine the dataset name using the input file path
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns>Dataset name</returns>
        /// <remarks>Returns an empty string if unable to determine the dataset name</remarks>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineDatasetName")]
        public static string AutoDetermineDatasetName(string filePath)
        {
            var resultType = AutoDetermineResultType(filePath);
            return AutoDetermineDatasetName(filePath, resultType);
        }

        /// <summary>
        /// Auto-determine the dataset name using the input file path and specified PeptideHit result type
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="resultType"></param>
        /// <returns>Dataset name</returns>
        /// <remarks>Returns an empty string if unable to determine the dataset name</remarks>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineDatasetName")]
        public static string AutoDetermineDatasetName(string filePath, PeptideHitResultTypes resultType)
        {
            return ReaderFactory.AutoDetermineDatasetName(filePath, resultType);
        }

        /// <summary>
        /// Determine the PeptideHit result type given the input file path
        /// </summary>
        /// <param name="filePath"></param>
        [Obsolete("Superseded by ReaderFactory.AutoDetermineResultType")]
        public static PeptideHitResultTypes AutoDetermineResultType(string filePath)
        {
            return ReaderFactory.AutoDetermineResultType(filePath);
        }

        /// <summary>
        /// Find the ModSummary file for the given input file
        /// </summary>
        /// <param name="peptideHitResultType">PHRP Result Type of the input file</param>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputDirectoryPath">Input directory (can be an empty string if inputFileName is a full path)</param>
        /// <param name="inputFileName">Name or path of the input file, e.g. Dataset_msgfplus_syn.txt or Dataset_syn.txt</param>
        /// <param name="modSummaryFileNamePreferred">Output: preferred mod summary filename (based on whether a _syn.txt or _fht.txt file is present)</param>
        [Obsolete("Superseded by ReaderFactory.FindModSummaryFile")]
        public static string FindModSummaryFile(
            PeptideHitResultTypes peptideHitResultType,
            string datasetName,
            string inputDirectoryPath,
            string inputFileName,
            out string modSummaryFileNamePreferred)
        {
            modSummaryFileNamePreferred = string.Empty;
            return string.Empty;
        }

        /// <summary>
        /// Find the ModSummary file for the given input file
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory (can be an empty string if inputFileName is a full path)</param>
        /// <param name="inputFileName">Name or path of the input file, e.g. Dataset_msgfplus_syn.txt or Dataset_syn.txt</param>
        /// <param name="modSummaryFileName">Expected mod summary filename</param>
        /// <param name="modSummaryFileNamePreferred">Output: preferred mod summary filename (based on whether a _syn.txt or _fht.txt file is present)</param>
        /// <returns>Mod summary file path if found; otherwise, an empty string</returns>
        [Obsolete("Superseded by ReaderFactory.FindModSummaryFile")]
        public static string FindModSummaryFile(
            string inputDirectoryPath,
            string inputFileName,
            string modSummaryFileName,
            out string modSummaryFileNamePreferred)
        {
            return FindPHRPFile(inputDirectoryPath, inputFileName, modSummaryFileName, out modSummaryFileNamePreferred);
        }

        /// <summary>
        /// Find the ResultToSeqMap file for the given input file
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory (can be an empty string if inputFileName is a full path)</param>
        /// <param name="inputFileName">Name or path of the input file, e.g. Dataset_msgfplus_syn.txt or Dataset_syn.txt</param>
        /// <param name="resultToSeqMapFileName">Expected ResultToSeqMap filename</param>
        /// <param name="resultToSeqMapFileNamePreferred">Output: preferred ResultToSeqMap filename (based on whether a _syn.txt or _fht.txt file is present)</param>
        /// <returns>Mod summary file path if found; otherwise, an empty string</returns>
        [Obsolete("Superseded by ReaderFactory.FindResultToSeqMapFile")]
        public static string FindResultToSeqMapFile(
            string inputDirectoryPath,
            string inputFileName,
            string resultToSeqMapFileName,
            out string resultToSeqMapFileNamePreferred)
        {
            return FindPHRPFile(inputDirectoryPath, inputFileName, resultToSeqMapFileName, out resultToSeqMapFileNamePreferred);
        }

        /// <summary>
        /// Find the given PHRP result file for the given input file
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory (can be an empty string if inputFileName is a full path)</param>
        /// <param name="inputFileName">Name or path of the input file, e.g. Dataset_msgfplus_syn.txt or Dataset_syn.txt</param>
        /// <param name="fileNameToFind">Expected PHRP result filename</param>
        /// <param name="preferredName">Output: preferred PHRP result filename (based on whether a _syn.txt or _fht.txt file is present)</param>
        /// <returns>Mod summary file path if found; otherwise, an empty string</returns>
        [Obsolete("Superseded by ReaderFactory.FindPHRPFile")]
        public static string FindPHRPFile(
            string inputDirectoryPath,
            string inputFileName,
            string fileNameToFind,
            out string preferredName)
        {
            return ReaderFactory.FindPHRPFile(inputDirectoryPath, inputFileName, fileNameToFind, out preferredName);
        }

        /// <summary>
        /// Returns the filename of the MSGF file that corresponds to synopsisOrFirstHitsFileName
        /// </summary>
        /// <param name="synopsisOrFirstHitsFileName">Filename (or full path) to the synopsis or first-hits file</param>
        [Obsolete("Superseded by ReaderFactory.GetMSGFFileName")]
        public static string GetMSGFFileName(string synopsisOrFirstHitsFileName)
        {
            return ReaderFactory.GetMSGFFileName(synopsisOrFirstHitsFileName);
        }

        /// <summary>
        /// Get the peptide hit result type for the given result type name
        /// </summary>
        /// <param name="resultTypeName"></param>
        [Obsolete("Superseded by ReaderFactory.GetPeptideHitResultType")]
        public static PeptideHitResultTypes GetPeptideHitResultType(string resultTypeName)
        {
            return ReaderFactory.GetPeptideHitResultType(resultTypeName);
        }

        /// <summary>
        /// Get the list of auxiliary file suffixes
        /// </summary>
        [Obsolete("Superseded by ReaderFactory.GetPHRPAuxiliaryFileSuffixes")]
        public static List<string> GetPHRPAuxiliaryFileSuffixes()
        {
            return ReaderFactory.GetPHRPAuxiliaryFileSuffixes();
        }

        /// <summary>
        /// Returns the default first-hits file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPFirstHitsFileName")]
        public static string GetPHRPFirstHitsFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPFirstHitsFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default ModSummary file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPModSummaryFileName")]
        public static string GetPHRPModSummaryFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPModSummaryFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default PepToProtMap file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPPepToProteinMapFileName")]
        public static string GetPHRPPepToProteinMapFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPPepToProteinMapFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default ProteinMods file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPProteinModsFileName")]
        public static string GetPHRPProteinModsFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPProteinModsFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default Synopsis file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPSynopsisFileName")]
        public static string GetPHRPSynopsisFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPSynopsisFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default ResultToSeq Map file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPResultToSeqMapFileName")]
        public static string GetPHRPResultToSeqMapFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPResultToSeqMapFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default SeqInfo file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPSeqInfoFileName")]
        public static string GetPHRPSeqInfoFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPSeqInfoFileName(resultType, datasetName);
        }

        /// <summary>
        /// Returns the default SeqToProtein Map file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetPHRPSeqToProteinMapFileName")]
        public static string GetPHRPSeqToProteinMapFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            return ReaderFactory.GetPHRPSeqToProteinMapFileName(resultType, datasetName);
        }

        /// <summary>
        /// Get the ScanStats filename for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetScanStatsFilename")]
        public static string GetScanStatsFilename(string datasetName)
        {
            return ReaderFactory.GetScanStatsFilename(datasetName);
        }

        /// <summary>
        /// Get the extended ScanStats filename for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetExtendedScanStatsFilename")]
        public static string GetExtendedScanStatsFilename(string datasetName)
        {
            return ReaderFactory.GetExtendedScanStatsFilename(datasetName);
        }

        /// <summary>
        /// Get the tool version info filename for the given analysis tool
        /// </summary>
        /// <param name="resultType"></param>
        /// <returns>Filename</returns>
        [Obsolete("Superseded by ReaderFactory.GetToolVersionInfoFilename")]
        public static string GetToolVersionInfoFilename(PeptideHitResultTypes resultType)
        {
            return ReaderFactory.GetToolVersionInfoFilename(resultType);
        }

        /// <summary>
        /// Returns true if the character is a letter between A and Z or a and z
        /// </summary>
        /// <param name="chChar">Character to examine</param>
        /// <remarks>The Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters; use this function to only return True for normal letters between A and Z</remarks>
        [Obsolete("Superseded by ReaderFactory.IsLetterAtoZ")]
        public static bool IsLetterAtoZ(char chChar)
        {
            return ReaderFactory.IsLetterAtoZ(chChar);
        }

        /// <summary>
        /// Examines the string to determine if it is numeric
        /// </summary>
        /// <param name="data"></param>
        /// <returns>True if a number, otherwise false</returns>
        [Obsolete("Superseded by ReaderFactory.IsNumber")]
        public static bool IsNumber(string data)
        {
            return ReaderFactory.IsNumber(data);
        }

        /// <summary>
        /// Returns the index of the indicated column, as tracked by columnHeaders
        /// </summary>
        /// <param name="columnName"></param>
        /// <param name="columnHeaders"></param>
        /// <returns>Column index, or -1 if not found</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnIndex")]
        public static int LookupColumnIndex(string columnName, SortedDictionary<string, int> columnHeaders)
        {
            return ReaderFactory.LookupColumnIndex(columnName, columnHeaders);
        }

        /// <summary>
        /// Returns the index of the indicated column, as tracked by columnHeaders
        /// </summary>
        /// <param name="columnEnum"></param>
        /// <param name="columnHeaders"></param>
        /// <returns>Column index, or -1 if not found</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnIndex")]
        public static int LookupColumnIndex(Enum columnEnum, SortedDictionary<Enum, int> columnHeaders)
        {
            return ReaderFactory.LookupColumnIndex(columnEnum, columnHeaders);
        }

        /// <summary>
        /// Returns the string stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; an empty string if the specific column name is not recognized</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnValue")]
        public static string LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders)
        {
            return ReaderFactory.LookupColumnValue(dataColumns, columnName, columnHeaders);
        }

        /// <summary>
        /// Returns the string stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; valueIfMissing if the specific column name is not recognized</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnValue")]
        public static string LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders,
                                               string valueIfMissing)
        {
            return ReaderFactory.LookupColumnValue(dataColumns, columnName, columnHeaders, valueIfMissing);
        }

        /// <summary>
        /// Returns the string stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; valueIfMissing if the specific column name is not recognized</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnValue")]
        public static string LookupColumnValue(string[] dataColumns, Enum columnEnum, SortedDictionary<Enum, int> columnHeaders,
                                               string valueIfMissing)
        {
            return ReaderFactory.LookupColumnValue(dataColumns, columnEnum, columnHeaders, valueIfMissing);
        }

        /// <summary>
        /// Returns the value stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnValue")]
        public static int LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders, int valueIfMissing)
        {
            return ReaderFactory.LookupColumnValue(dataColumns, columnName, columnHeaders, valueIfMissing);
        }

        /// <summary>
        /// Returns the value stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
        [Obsolete("Superseded by ReaderFactory.LookupColumnValue")]
        public static double LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders,
                                               double valueIfMissing)
        {
            return ReaderFactory.LookupColumnValue(dataColumns, columnName, columnHeaders, valueIfMissing);
        }

        /// <summary>
        /// Updates the column name to column index mapping in columnHeaders
        /// </summary>
        /// <param name="dataColumns">Column names read from the input file</param>
        /// <param name="columnHeaders">Column mapping dictionary object to update</param>
        /// <remarks>The SortedDictionary object should be instantiated using a case-insensitive comparer, i.e. (StringComparer.OrdinalIgnoreCase)</remarks>
        [Obsolete("Superseded by ReaderFactory.ParseColumnHeaders")]
        public static void ParseColumnHeaders(string[] dataColumns, SortedDictionary<string, int> columnHeaders)
        {
            ReaderFactory.ParseColumnHeaders(dataColumns, columnHeaders);
            return;
        }

        /// <summary>
        /// Reads the next line from a synopsis file or first hits file
        /// </summary>
        /// <returns>True if a line was read, false if not more data is available</returns>
        /// <remarks>When FastReadMode is True, you should call FinalizeCurrentPSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
        [Obsolete("Superseded by ReaderFactory.MoveNext")]
        public bool MoveNext()
        {
            return false;
        }

        /// <summary>
        /// When FastReadMode is True, first call MoveNext to read the peptide scores.
        /// Then, if the peptide is a peptide of interest, call this function to finalize any processing steps that were skipped.
        /// </summary>
        [Obsolete("Superseded by ReaderFactory.FinalizeCurrentPSM")]
        public void FinalizeCurrentPSM()
        {
            throw new Exception("Method FinalizeCurrentPSM is obsolete");
        }

        #region "IDisposable Support"

        /// <summary>
        /// Used to detect redundant calls
        /// </summary>
#pragma warning disable 414
        private bool mDisposedValue;
#pragma warning restore 414

        /// <summary>
        /// Dispose of this class
        /// </summary>
        /// <param name="disposing"></param>
        protected virtual void Dispose(bool disposing)
        {
            mDisposedValue = true;
        }

        /// <summary>
        /// This code added by Visual Studio to correctly implement the disposable pattern.
        /// </summary>
        public void Dispose()
        {
            // Do not change this code.  Put cleanup code in Dispose(disposing As Boolean) above.
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        #endregion
    }

    /// <summary>
    /// Legacy PHRP startup options
    /// </summary>
    [Obsolete("Legacy Startup Options; superseded by StartupOptions")]
    public class clsPHRPStartupOptions
    {
        /// <summary>
        /// If true, load the modification and SeqInfo data
        /// </summary>
        public bool LoadModsAndSeqInfo { get; set; }

        /// <summary>
        /// If true, load MSGF results (not MS-GF+)
        /// </summary>
        public bool LoadMSGFResults { get; set; }

        /// <summary>
        /// If true, load ScanStats data
        /// </summary>
        public bool LoadScanStatsData { get; set; }

        /// <summary>
        /// Maximum number of proteins to associate with each PSM
        /// </summary>
        /// <remarks>Set to 0 to load all proteins</remarks>
        public int MaxProteinsPerPSM { get; set; }

        /// <summary>
        /// Use this to override the default peptide mass calculator class;
        /// this is useful if custom amino acids are in use
        /// </summary>
        public PeptideMassCalculator PeptideMassCalculator { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public clsPHRPStartupOptions()
        {
            LoadModsAndSeqInfo = true;
            LoadMSGFResults = true;
            LoadScanStatsData = false;
            MaxProteinsPerPSM = 0;      // 0 means to load all proteins
        }
    }

#pragma warning restore IDE0060 // Remove unused parameter
#pragma warning restore RCS1163 // Unused parameter.
}
