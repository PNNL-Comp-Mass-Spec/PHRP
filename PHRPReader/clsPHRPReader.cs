//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

// ReSharper disable UnusedMember.Global

namespace PHRPReader
{
    /// <summary>
    ///  This class reads a tab-delimited text file (created by the Peptide File Extractor or by PHRP)
    ///  and returns the data for each peptide hit search result
    ///
    ///  It also integrates MSGF results with the peptide hit search results
    ///  And, it integrates scan stats values (to determine elution time)
    /// </summary>
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

        // This RegEx is used to extract parent ion m/z from a filter string that does not contain msx
        // ${ParentMZ} will hold the last parent ion m/z found
        // For example, 756.71 in FTMS + p NSI d Full ms3 850.70@cid35.00 756.71@cid35.00 [195.00-2000.00]
        private const string PARENT_ION_ONLY_NON_MSX_REGEX = @"[Mm][Ss]\d*[^\[\r\n]* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*(\[[^\]\r\n]\])?";

        // This RegEx is used to extract parent ion m/z from a filter string that does contain msx
        // ${ParentMZ} will hold the first parent ion m/z found (the first parent ion m/z corresponds to the highest peak)
        // For example, 636.04 in FTMS + p NSI Full msx ms2 636.04@hcd28.00 641.04@hcd28.00 654.05@hcd28.00 [88.00-1355.00]
        private const string PARENT_ION_ONLY_MSX_REGEX = @"[Mm][Ss]\d* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*[^\[\r\n]*(\[[^\]\r\n]+\])?";

#pragma warning disable 1591

        /// <summary>
        /// Peptide hit results type
        /// </summary>
        [Obsolete("Superseded by enum PeptideHitResultTypes")]
        public enum ePeptideHitResultType
        {
            Unknown = 0,
            Sequest = 1,
            XTandem = 2,
            Inspect = 3,
            [Obsolete("Use MSGFPlus")]
            MSGFDB = 4,
            MSGFPlus = 4,      // Aka MS-GF+
            MSAlign = 5,
            MODa = 6,
            MODPlus = 7,
            MSPathFinder = 8,
            TopPIC = 9
        }

        /// <summary>
        /// Peptide hit results type
        /// </summary>
        public enum PeptideHitResultTypes
        {
            Unknown = 0,
            Sequest = 1,
            XTandem = 2,
            Inspect = 3,
            [Obsolete("Use MSGFPlus")]
            MSGFDB = 4,
            MSGFPlus = 4,      // Aka MS-GF+
            MSAlign = 5,
            MODa = 6,
            MODPlus = 7,
            MSPathFinder = 8,
            TopPIC = 9
        }

        /// <summary>
        /// PHRP Reader error codes
        /// </summary>
        public enum PHRPReaderErrorCodes
        {
            NoError = 0,
            InvalidInputFilePath = 1,
            InputFileFormatNotRecognized = 2,
            RequiredInputFileNotFound = 3,
            MissingRawOrMzXmlFile = 4,
            MSGFProgramNotFound = 5,
            UnspecifiedError = -1
        }

#pragma warning restore 1591

        #endregion

        #region "Module variables"

        private string mInputFilePath;

        private string mInputDirectoryPath;

        private readonly clsPHRPStartupOptions mStartupOptions;

        private bool mInitialized;

        /// <summary>
        /// When set to true, calls to MoveNext will read the next data line, but will skip several additional processing steps for performance reasons
        /// </summary>
        /// <remarks>If the peptide is a peptide of interest, you must call FinalizeCurrentPSM after calling .MoveNext()</remarks>
        private bool mFastReadMode;

        private StreamReader mSourceFile;

        private int mSourceFileLineCount;

        private int mSourceFileLinesRead;

        private readonly clsPeptideMassCalculator mPeptideMassCalculator;

        // This dictionary contains mod symbols as the key and modification definition as the values
        private readonly SortedDictionary<char, clsModificationDefinition> mDynamicMods;

        // This dictionary contains amino acid names as the key and the corresponding mod modification (or mod modifications)
        private readonly SortedDictionary<string, List<clsModificationDefinition>> mStaticMods;

        // This dictionary tracks the MSGFSpecEvalue values for each entry in the source file
        // The keys are Result_ID and the string is MSGFSpecEValue (stored as string to preserve formatting)
        private Dictionary<int, string> mMSGFCachedResults;

        // This dictionary tracks scan stats values, in particular elution time
        //The keys are ScanNumber and values are clsScanStatsInfo objects
        private Dictionary<int, clsScanStatsInfo> mScanStats;

        // This dictionary tracks extended scan stats values, including parent ion mz (via MonoisotopicMZ)and collision mode
        //The keys are ScanNumber and values are clsScanStatsExInfo objects
        private Dictionary<int, clsScanStatsExInfo> mScanStatsEx;

        private clsPSM mPSMCurrent;
        private bool mPSMCurrentFinalized;

        private bool mExtendedScanStatsValid;
        private clsScanStatsExInfo mExtendedScanStatsInfo;

        private bool mHeaderLineParsed;
        private bool mCachedLineAvailable;
        private string mCachedLine;
        private clsPSM mCachedPSM;

        private PHRPReaderErrorCodes mLocalErrorCode;

        /// <summary>
        /// RegEx to extract parent ions from filter strings that do not have Full msx
        /// </summary>
        /// <remarks>Shared (aka static) only to speed up unit tests</remarks>
        private static readonly Regex mFindParentIonOnlyNonMsx = new Regex(PARENT_ION_ONLY_NON_MSX_REGEX, RegexOptions.Compiled | RegexOptions.IgnoreCase);

        /// <summary>
        /// RegEx to extract parent ions from filter strings that have Full msx
        /// </summary>
        /// <remarks>Shared (aka static) only to speed up unit tests</remarks>
        private static readonly Regex mFindParentIonOnlyMsx = new Regex(PARENT_ION_ONLY_MSX_REGEX, RegexOptions.Compiled | RegexOptions.IgnoreCase);

        #endregion

        #region "Properties"

        /// <summary>
        /// Returns True if the input file was successfully opened and data remains to be read
        /// </summary>
        /// <returns>True if the file is readable</returns>
        public bool CanRead { get; private set; }

        /// <summary>
        /// Returns the most recently loaded PSM
        /// </summary>
        public clsPSM CurrentPSM => mPSMCurrent;

        /// <summary>
        /// Returns the most recently loaded PSM's sequence info (if available)
        /// </summary>
        public clsSeqInfo CurrentPSMSeqInfo
        {
            get
            {
                if (mPSMCurrent == null || PHRPParser.SeqInfo == null)
                {
                    return null;
                }

                PHRPParser.SeqInfo.TryGetValue(mPSMCurrent.SeqID, out var seqInfo);
                return seqInfo;
            }
        }

        /// <summary>
        /// Dataset name (auto-determined based on the input filename)
        /// </summary>
        public string DatasetName { get; private set; }

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
        public string ErrorMessage { get; private set; } = string.Empty;

        /// <summary>
        /// Used to enable fast read mode when calling MoveNext
        /// When FastReadMode is True, you should call FinalizeCurrentPSM after calling MoveNext to populate the remaining fields if the peptide is a peptide of interest
        /// </summary>
        /// <remarks>Once FastReadMode is enabled it cannot be turned off (this is a safety measure due to how data is cached)</remarks>
        public bool FastReadMode
        {
            get => mFastReadMode;
            set
            {
                if (value)
                {
                    mFastReadMode = true;
                }
            }
        }

        /// <summary>
        /// If True, looks for and loads the modification definitions from the _ModSummary.txt file associated with the input file
        /// Also reads the SeqInfo and related files
        /// </summary>
        public bool LoadModsAndSeqInfo => mStartupOptions.LoadModsAndSeqInfo;

        /// <summary>
        /// If true, loads the MSGF SpecProb values from the _MSGF.txt file associated with the input file
        /// </summary>
        public bool LoadMSGFResults => mStartupOptions.LoadMSGFResults;

        /// <summary>
        /// If True, loads the MASIC _ScanStats.txt file
        /// </summary>
        public bool LoadScanStatsData => mStartupOptions.LoadScanStatsData;

        /// <summary>
        /// The maximum number of proteins that will be tracked for each PSM
        /// </summary>
        public int MaxProteinsPerPSM
        {
            get => mStartupOptions.MaxProteinsPerPSM;
            set => mStartupOptions.MaxProteinsPerPSM = value;
        }

        /// <summary>
        /// Returns True if the ModSummary file was successfully loaded
        /// </summary>
        public bool ModSummaryFileLoaded { get; private set; }

        /// <summary>
        /// Peptide hit result type; Sequest, XTandem, Inspect, or MSGFPlus (aka MS-GF+)
        /// </summary>
        public PeptideHitResultTypes PeptideHitResultType
        {
            get
            {
                return (PHRPParser?.PeptideHitResultType) ?? PeptideHitResultTypes.Unknown;
            }
        }

        /// <summary>
        /// Returns a number between 0 and 100 indicating the percentage of the source file that has been read
        /// </summary>
        public float PercentComplete
        {
            get
            {
                return mSourceFileLineCount > 0 ? mSourceFileLinesRead / Convert.ToSingle(mSourceFileLineCount) * 100f : 0;
            }
        }

        /// <summary>
        /// Returns the PHRP Parser object
        /// </summary>
        public clsPHRPParser PHRPParser { get; private set; }

        /// <summary>
        /// Returns the cached mapping between ResultID and SeqID
        /// </summary>
        public SortedList<int, int> ResultToSeqMap
        {
            get
            {
                return PHRPParser == null ? new SortedList<int, int>() : PHRPParser.ResultToSeqMap;
            }
        }

        /// <summary>
        /// Returns the cached sequence info, where key is SeqID
        /// </summary>
        public SortedList<int, clsSeqInfo> SeqInfo
        {
            get
            {
                return PHRPParser == null ? new SortedList<int, clsSeqInfo>() : PHRPParser.SeqInfo;
            }
        }

        /// <summary>
        /// Returns the cached sequence to protein map information
        /// </summary>
        public SortedList<int, List<clsProteinInfo>> SeqToProteinMap
        {
            get
            {
                return PHRPParser == null ? new SortedList<int, List<clsProteinInfo>>() : PHRPParser.SeqToProteinMap;
            }
        }

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
        public clsPHRPReader(string inputFilePath, PeptideHitResultTypes resultType, bool loadModsAndSeqInfo, bool loadMSGFResults, bool loadScanStats)
        {
            var startupOptions = new clsPHRPStartupOptions
            {
                LoadModsAndSeqInfo = loadModsAndSeqInfo,
                LoadMSGFResults = loadMSGFResults,
                LoadScanStatsData = loadScanStats
            };

            mStartupOptions = startupOptions;

            mMSGFCachedResults = new Dictionary<int, string>();

            mDynamicMods = new SortedDictionary<char, clsModificationDefinition>();
            mStaticMods = new SortedDictionary<string, List<clsModificationDefinition>>();

            mPeptideMassCalculator = new clsPeptideMassCalculator();

            ErrorMessages = new List<string>();
            WarningMessages = new List<string>();

            InitializeClass(inputFilePath, resultType);
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="resultType">Source file PeptideHit result type</param>
        /// <param name="startupOptions">Startup options</param>
        public clsPHRPReader(string inputFilePath, PeptideHitResultTypes resultType, clsPHRPStartupOptions startupOptions)
        {
            mStartupOptions = startupOptions ?? throw new ArgumentNullException(nameof(startupOptions));

            mMSGFCachedResults = new Dictionary<int, string>();

            mDynamicMods = new SortedDictionary<char, clsModificationDefinition>();
            mStaticMods = new SortedDictionary<string, List<clsModificationDefinition>>();

            mPeptideMassCalculator = startupOptions.PeptideMassCalculator ?? new clsPeptideMassCalculator();

            ErrorMessages = new List<string>();
            WarningMessages = new List<string>();

            InitializeClass(inputFilePath, resultType);
        }

        /// <summary>
        /// Updates filePath to have _fht instead of _syn if filePath contains_syn yet basePHRPFileName contains _fht
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="basePHRPFileName"></param>
        private static string AutoSwitchToFHTIfRequired(string filePath, string basePHRPFileName)
        {
            if (string.IsNullOrEmpty(basePHRPFileName))
            {
                return filePath;
            }

            var basePHRPFile = new FileInfo(basePHRPFileName);
            if (basePHRPFile.Name.IndexOf("_fht", StringComparison.OrdinalIgnoreCase) >= 0)
            {
                // basePHRPFileName is first-hits-file based

                var firstHitsFile = new FileInfo(filePath);
                var synIndex = firstHitsFile.Name.LastIndexOf("_syn", StringComparison.OrdinalIgnoreCase);
                if (synIndex > 0)
                {
                    // filePath is synopsis-file based
                    // Change filePath to contain _fht instead of _syn

                    var filePathFHT = firstHitsFile.Name.Substring(0, synIndex) + "_fht" + firstHitsFile.Name.Substring(synIndex + "_syn".Length);

                    if (Path.IsPathRooted(filePath) && firstHitsFile.Directory != null)
                    {
                        return Path.Combine(firstHitsFile.Directory.FullName, filePathFHT);
                    }

                    return filePathFHT;
                }
            }

            return filePath;
        }

        /// <summary>
        /// Updates filePath to have _msgfdb instead of _msgfplus if basePHRPFileName contains _msgfdb
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="basePHRPFileName"></param>
        public static string AutoSwitchToLegacyMSGFDBIfRequired(string filePath, string basePHRPFileName)
        {
            var basePHRPFile = new FileInfo(basePHRPFileName);
            if (basePHRPFile.Name.IndexOf("_msgfdb", StringComparison.OrdinalIgnoreCase) >= 0)
            {
                var dataFile = new FileInfo(filePath);
                var charIndex = dataFile.Name.LastIndexOf("_msgfplus", StringComparison.OrdinalIgnoreCase);
                if (charIndex > 0)
                {
                    // filePath has _msgfplus but should have _msgfdb

                    var filePathNew = dataFile.Name.Substring(0, charIndex) + "_msgfdb" + dataFile.Name.Substring(charIndex + "_msgfplus".Length);

                    if (Path.IsPathRooted(filePath) && dataFile.Directory != null)
                    {
                        return Path.Combine(dataFile.Directory.FullName, filePathNew);
                    }

                    return filePathNew;
                }
            }

            return filePath;
        }

        /// <summary>
        /// Clear any cached error messages
        /// </summary>
        public void ClearErrors()
        {
            ErrorMessages.Clear();
            PHRPParser?.ClearErrors();
        }

        private int CountLines(string textFilePath)
        {
            int totalLines;

            try
            {
                totalLines = 0;
                using (var reader = new StreamReader(new FileStream(textFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        reader.ReadLine();
                        totalLines++;
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error counting the lines in " + Path.GetFileName(textFilePath) + ": " + ex.Message, ex);
            }

            return totalLines;
        }

        /// <summary>
        /// Clear any cached warning messages
        /// </summary>
        public void ClearWarnings()
        {
            WarningMessages.Clear();
            PHRPParser?.ClearWarnings();
        }

        /// <summary>
        /// Initialize the class
        /// </summary>
        /// <param name="inputFilePath">Input file to read</param>
        /// <param name="resultType">Source file PeptideHit result type</param>
        private void InitializeClass(string inputFilePath, PeptideHitResultTypes resultType)
        {
            mInitialized = false;

            InitializeMemberVariables();

            InitializeReader(inputFilePath, resultType);

            mInitialized = true;
        }

        private void InitializeMemberVariables()
        {
            DatasetName = string.Empty;
            mInputFilePath = string.Empty;
            mInputDirectoryPath = string.Empty;

            CanRead = false;
            ModSummaryFileLoaded = false;

            SkipDuplicatePSMs = true;

            EchoMessagesToConsole = false;

            ErrorMessage = string.Empty;
            mLocalErrorCode = PHRPReaderErrorCodes.NoError;

            mSourceFileLineCount = 0;
        }

        private void InitializeReader(string inputFilePath, PeptideHitResultTypes resultType)
        {
            var modSummaryFilePath = string.Empty;

            try
            {
                if (string.IsNullOrEmpty(inputFilePath))
                {
                    ReportError("Input file name is empty");
                    SetLocalErrorCode(PHRPReaderErrorCodes.InvalidInputFilePath);
                    if (!mInitialized)
                        throw new FileNotFoundException(ErrorMessage);
                    return;
                }

                // Confirm that the source file exists
                // Make sure inputFilePath points to a valid file
                var inputFile = new FileInfo(inputFilePath);

                if (inputFile.Directory == null)
                {
                    ReportError("Unable to determine the parent directory of " + inputFile.FullName);
                    SetLocalErrorCode(PHRPReaderErrorCodes.InvalidInputFilePath);
                    if (!mInitialized)
                        throw new FileNotFoundException(ErrorMessage);
                    return;
                }

                mInputDirectoryPath = inputFile.Directory.FullName;
                mInputFilePath = inputFile.FullName;

                if (!inputFile.Exists)
                {
                    ReportError("Input file not found: " + inputFilePath);
                    if (inputFilePath.Contains(".."))
                    {
                        ReportWarning("Absolute path: " + inputFile.DirectoryName);
                    }
                    SetLocalErrorCode(PHRPReaderErrorCodes.InvalidInputFilePath);
                    if (!mInitialized)
                        throw new FileNotFoundException(ErrorMessage);
                    return;
                }

                // Note that the following populates mDatasetName
                var success = ValidateInputFiles(inputFilePath, ref resultType, ref modSummaryFilePath);
                if (!success)
                {
                    SetLocalErrorCode(PHRPReaderErrorCodes.RequiredInputFileNotFound, true);
                    if (!mInitialized)
                        throw new FileNotFoundException(ErrorMessage);
                    return;
                }

                // Open the input file for reading
                // Note that this will also load the MSGFSpecEValue info and ScanStats info
                success = InitializeParser(resultType);

                if (success && mStartupOptions.LoadModsAndSeqInfo)
                {
                    // Read the PHRP Mod Summary File to populate mDynamicMods and mStaticMods
                    // Note that the PHRPParser also loads the ModSummary file, and that mDynamicMods and mStaticMods are only used if the _SeqInfo.txt file is not found
                    success = ReadModSummaryFile(modSummaryFilePath, mDynamicMods, mStaticMods);
                    if (!success)
                    {
                        ModSummaryFileLoaded = false;
                        success = true;
                    }
                    else
                    {
                        ModSummaryFileLoaded = true;
                    }
                }

                if (success && mStartupOptions.LoadMSGFResults)
                {
                    // Cache the MSGF values (if present)
                    ReadAndCacheMSGFData();
                }

                if (success && mStartupOptions.LoadScanStatsData)
                {
                    // Cache the Scan Stats values (if present)
                    ReadScanStatsData();
                    ReadExtendedScanStatsData();
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in InitializeReader", ex);
                if (!mInitialized)
                    throw new Exception(ErrorMessage, ex);
            }
        }

        private bool InitializeParser(PeptideHitResultTypes resultType)
        {
            var success = true;
            var datasetName = string.Copy(DatasetName);

            try
            {
                if (string.IsNullOrEmpty(datasetName))
                {
                    if (mStartupOptions.LoadModsAndSeqInfo)
                    {
                        ReportError("Dataset name is undefined; unable to continue since loading ModsAndSeqInfo");
                        return false;
                    }

                    if (mStartupOptions.LoadMSGFResults)
                    {
                        ReportError("Dataset name is undefined; unable to continue since loading MSGF results");
                        return false;
                    }

                    if (mStartupOptions.LoadScanStatsData)
                    {
                        ReportError("Dataset name is undefined; unable to continue since loading ScanStatsData");
                        return false;
                    }

                    datasetName = "Unknown_Dataset";
                }

                // Initialize some tracking variables
                mMSGFCachedResults.Clear();

                mPSMCurrent = new clsPSM();

                mSourceFileLinesRead = 0;
                mHeaderLineParsed = false;
                mCachedLineAvailable = false;
                mCachedLine = string.Empty;

                // Open the peptide-hit result file (from PHRP) for reading
                // Instantiate the appropriate PHRP Parser
                switch (resultType)
                {
                    case PeptideHitResultTypes.Sequest:
                        PHRPParser = new clsPHRPParserSequest(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.XTandem:
                        // Note that Result to Protein mapping will be auto-loaded during instantiation of mPHRPParser
                        PHRPParser = new clsPHRPParserXTandem(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.Inspect:
                        PHRPParser = new clsPHRPParserInspect(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.MSGFPlus:
                        // MS-GF+
                        PHRPParser = new clsPHRPParserMSGFPlus(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.MSAlign:
                        PHRPParser = new clsPHRPParserMSAlign(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.MODa:
                        PHRPParser = new clsPHRPParserMODa(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.MODPlus:
                        PHRPParser = new clsPHRPParserMODPlus(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.MSPathFinder:
                        PHRPParser = new clsPHRPParserMSPathFinder(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    case PeptideHitResultTypes.TopPIC:
                        PHRPParser = new clsPHRPParserTopPIC(datasetName, mInputFilePath, mStartupOptions);
                        break;

                    default:
                        //Should never get here; invalid result type specified
                        ReportError("Invalid PeptideHit ResultType specified: " + resultType);
                        success = false;
                        break;
                }

                if (!success)
                {
                    return false;
                }

                // Attach the event handlers
                RegisterEvents(PHRPParser);

                // Report any errors cached during instantiation of mPHRPParser
                foreach (var message in PHRPParser.ErrorMessages)
                {
                    ReportError(message);
                }

                // Report any warnings cached during instantiation of mPHRPParser
                foreach (var message in PHRPParser.WarningMessages)
                {
                    ReportWarning(message);
                }

                PHRPParser.ClearErrors();
                PHRPParser.ClearWarnings();

                // Open the data file and count the number of lines so that we can compute progress
                mSourceFileLineCount = CountLines(mInputFilePath);

                // Open the data file for reading
                mSourceFile = new StreamReader(new FileStream(mInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                CanRead = true;

                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error in clsPHRPReader.InitializeParser", ex);
                if (!mInitialized)
                    throw new Exception(ErrorMessage, ex);

                return false;
            }
        }

        /// <summary>
        /// Add a PHRP file to find when auto-determining the input file type
        /// </summary>
        /// <param name="filesToFind"></param>
        /// <param name="resultType"></param>
        /// <param name="getPhrpFilenameMethod"></param>
        /// <param name="datasetName"></param>
        private static void AddFileToFind(
            ICollection<KeyValuePair<string, PeptideHitResultTypes>> filesToFind,
            PeptideHitResultTypes resultType,
            Func<string, string> getPhrpFilenameMethod,
            string datasetName)
        {
            var phrpFilename = getPhrpFilenameMethod(datasetName);
            filesToFind.Add(new KeyValuePair<string, PeptideHitResultTypes>(phrpFilename, resultType));
        }

        /// <summary>
        /// Add a PHRP file to find when auto-determining the input file type
        /// </summary>
        /// <param name="filesToFind"></param>
        /// <param name="resultType"></param>
        /// <param name="phrpFilename"></param>
        private static void AddFileToFind(
            ICollection<KeyValuePair<string, PeptideHitResultTypes>> filesToFind,
            PeptideHitResultTypes resultType,
            string phrpFilename)
        {
            filesToFind.Add(new KeyValuePair<string, PeptideHitResultTypes>(phrpFilename, resultType));
        }

        /// <summary>
        /// Add a filename suffix to check when auto-determining the input file type
        /// </summary>
        /// <param name="suffixesToCheck"></param>
        /// <param name="resultType"></param>
        /// <param name="suffix"></param>
        private static void AddSuffixToCheck(
            ICollection<KeyValuePair<string, PeptideHitResultTypes>> suffixesToCheck,
            PeptideHitResultTypes resultType,
            string suffix)
        {
            suffixesToCheck.Add(new KeyValuePair<string, PeptideHitResultTypes>(suffix, resultType));
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
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
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, out PeptideHitResultTypes matchedResultType)
        {
            // Find candidate dataset names in inputDirectoryPath

            var datasetNames = new SortedSet<string>(StringComparer.OrdinalIgnoreCase);
            var filesToFind = new List<string>();

            if (string.IsNullOrWhiteSpace(inputDirectoryPath))
            {
                throw new DirectoryNotFoundException("Input directory path is empty");
            }

            var inputDirectory = new DirectoryInfo(inputDirectoryPath);
            if (!inputDirectory.Exists)
            {
                throw new DirectoryNotFoundException("Input directory not found: " + inputDirectoryPath);
            }

            // MS-GF+
            filesToFind.Add(clsPHRPParserMSGFPlus.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserMSGFPlus.FILENAME_SUFFIX_FHT);

            // MS-GF+ prior to November 2016
            filesToFind.Add(GetLegacyMSGFPlusName(clsPHRPParserMSGFPlus.FILENAME_SUFFIX_SYN));
            filesToFind.Add(GetLegacyMSGFPlusName(clsPHRPParserMSGFPlus.FILENAME_SUFFIX_FHT));

            // X!Tandem (only has _xt.txt files)
            filesToFind.Add(clsPHRPParserXTandem.FILENAME_SUFFIX_SYN);

            // MSAlign
            filesToFind.Add(clsPHRPParserMSAlign.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserMSAlign.FILENAME_SUFFIX_FHT);

            // Inspect
            filesToFind.Add(clsPHRPParserInspect.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserInspect.FILENAME_SUFFIX_FHT);

            // MODa
            filesToFind.Add(clsPHRPParserMODa.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserMODa.FILENAME_SUFFIX_FHT);

            // MODPlus
            filesToFind.Add(clsPHRPParserMODPlus.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserMODPlus.FILENAME_SUFFIX_FHT);

            // MSPathFinder
            filesToFind.Add(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_FHT);

            // TopPIC
            filesToFind.Add(clsPHRPParserTopPIC.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserTopPIC.FILENAME_SUFFIX_FHT);

            // *****************
            // ** Important: Sequest needs to be added last since files simply end in _syn.txt or _fht.txt)
            // *****************
            // Sequest
            filesToFind.Add(clsPHRPParserSequest.FILENAME_SUFFIX_SYN);
            filesToFind.Add(clsPHRPParserSequest.FILENAME_SUFFIX_FHT);

            foreach (var fileSpec in filesToFind)
            {
                foreach (var dataFile in inputDirectory.GetFiles("*" + fileSpec))
                {
                    var dataset = dataFile.Name;

                    var charIndex = dataset.ToLower().IndexOf(fileSpec, StringComparison.Ordinal);
                    if (charIndex > 0)
                    {
                        dataset = dataset.Substring(0, charIndex);

                        if (!datasetNames.Contains(dataset))
                        {
                            datasetNames.Add(dataset);
                        }
                    }
                }
            }

            if (datasetNames.Count == 0)
            {
                Console.WriteLine("Did not find any files matching the expected filename suffixes");
                Console.WriteLine("Looked for the following in " + inputDirectoryPath);
                foreach (var fileSpec in filesToFind)
                {
                    Console.WriteLine("  " + fileSpec);
                }
                matchedResultType = PeptideHitResultTypes.Unknown;
                return string.Empty;
            }

            return AutoDetermineBestInputFile(inputDirectoryPath, datasetNames.ToList(), out matchedResultType);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified directory
        /// If both the _syn.txt and _fht.txt files are present, chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
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
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, string datasetName,
            out PeptideHitResultTypes matchedResultType)
        {
            var datasetNames = new List<string> {
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
        public static string AutoDetermineBestInputFile(string inputDirectoryPath, List<string> datasetNames,
            out PeptideHitResultTypes matchedResultType)
        {
            // This list contains the standard PHRP file suffixes
            var auxiliaryFileSuffixes = GetPHRPAuxiliaryFileSuffixes();

            // The key in this variable is the full path to the best Synopsis or First hits file and the value is the number of PHRP-related auxiliary files
            var kvBestSynOrFHTFile = new KeyValuePair<string, int>(string.Empty, 0);

            // Set the matched result type to Unknown for now
            matchedResultType = PeptideHitResultTypes.Unknown;

            if (string.IsNullOrWhiteSpace(inputDirectoryPath))
            {
                throw new DirectoryNotFoundException("Input directory path is empty");
            }

            var inputDirectory = new DirectoryInfo(inputDirectoryPath);
            if (!inputDirectory.Exists)
            {
                throw new DirectoryNotFoundException("Input directory not found: " + inputDirectoryPath);
            }

            if (datasetNames == null || datasetNames.Count == 0)
            {
                throw new ArgumentException("List datasetNames cannot be empty; cannot determine the best input file");
            }

            // Construct a list of the files to search for
            // Items in this list are KeyValuePairs where the key is a filename to look for and the value is a PeptideHitResultTypes
            var filesToFind = new List<KeyValuePair<string, PeptideHitResultTypes>>();

            foreach (var dataset in datasetNames)
            {
                // MS-GF+
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSGFPlus, clsPHRPParserMSGFPlus.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSGFPlus, clsPHRPParserMSGFPlus.GetPHRPFirstHitsFileName, dataset);

                // MS-GF+ prior to November 2016
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSGFPlus, GetLegacyMSGFPlusName(clsPHRPParserMSGFPlus.GetPHRPSynopsisFileName(dataset)));
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSGFPlus, GetLegacyMSGFPlusName(clsPHRPParserMSGFPlus.GetPHRPFirstHitsFileName(dataset)));

                AddFileToFind(filesToFind, PeptideHitResultTypes.XTandem, clsPHRPParserXTandem.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.XTandem, clsPHRPParserXTandem.GetPHRPFirstHitsFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSAlign, clsPHRPParserMSAlign.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSAlign, clsPHRPParserMSAlign.GetPHRPFirstHitsFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MODa, clsPHRPParserMODa.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MODa, clsPHRPParserMODa.GetPHRPFirstHitsFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MODPlus, clsPHRPParserMODPlus.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MODPlus, clsPHRPParserMODPlus.GetPHRPFirstHitsFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSPathFinder, clsPHRPParserMSPathFinder.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.MSPathFinder, clsPHRPParserMSPathFinder.GetPHRPFirstHitsFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.TopPIC, clsPHRPParserTopPIC.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.TopPIC, clsPHRPParserTopPIC.GetPHRPFirstHitsFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.Inspect, clsPHRPParserInspect.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.Inspect, clsPHRPParserInspect.GetPHRPFirstHitsFileName, dataset);

                // Sequest (needs to be added last since files simply end in _syn.txt or _fht.txt)
                AddFileToFind(filesToFind, PeptideHitResultTypes.Sequest, clsPHRPParserSequest.GetPHRPSynopsisFileName, dataset);
                AddFileToFind(filesToFind, PeptideHitResultTypes.Sequest, clsPHRPParserSequest.GetPHRPFirstHitsFileName, dataset);
            }

            foreach (var kvFileToFind in filesToFind)
            {
                if (string.IsNullOrEmpty(kvFileToFind.Key))
                    continue;

                var synOrFHTFile = new FileInfo(Path.Combine(inputDirectory.FullName, kvFileToFind.Key));

                if (!synOrFHTFile.Exists || synOrFHTFile.Directory == null)
                    continue;

                // Match found
                // Look for PHRP-related auxiliary files
                var auxFileCount = 0;
                var baseName = Path.Combine(synOrFHTFile.Directory.FullName, Path.GetFileNameWithoutExtension(synOrFHTFile.Name));

                foreach (var suffix in auxiliaryFileSuffixes)
                {
                    if (File.Exists(baseName + suffix))
                    {
                        auxFileCount++;
                    }
                }

                if (string.IsNullOrEmpty(kvBestSynOrFHTFile.Key) || auxFileCount > kvBestSynOrFHTFile.Value)
                {
                    kvBestSynOrFHTFile = new KeyValuePair<string, int>(synOrFHTFile.FullName, auxFileCount);
                    matchedResultType = kvFileToFind.Value;
                }
            }

            if (!string.IsNullOrWhiteSpace(kvBestSynOrFHTFile.Key))
                return kvBestSynOrFHTFile.Key;

            if (datasetNames.Count == 1)
            {
                Console.WriteLine("Could not find a Synopsis or First Hits file for dataset " + datasetNames.First());
            }
            else
            {
                Console.WriteLine("Could not find a Synopsis or First Hits file for any of the candidate datasets");
            }

            Console.WriteLine("Looked for the following files:");
            foreach (var fileName in filesToFind)
            {
                if (!string.IsNullOrWhiteSpace(fileName.Key))
                {
                    Console.WriteLine("  " + fileName.Key);
                }
            }

            // kvBestSynOrFHTFile should now contain the PHRP result file with the most auxiliary files
            return kvBestSynOrFHTFile.Key;
        }

        /// <summary>
        /// Auto-determine the dataset name using the input file path
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns>Dataset name</returns>
        /// <remarks>Returns an empty string if unable to determine the dataset name</remarks>
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
        public static string AutoDetermineDatasetName(string filePath, PeptideHitResultTypes resultType)
        {
            var datasetName = string.Empty;

            var inputFileName = Path.GetFileNameWithoutExtension(filePath);

            if (string.IsNullOrWhiteSpace(inputFileName))
                return string.Empty;

            switch (resultType)
            {
                case PeptideHitResultTypes.Sequest:
                case PeptideHitResultTypes.Inspect:
                case PeptideHitResultTypes.MSGFPlus:
                case PeptideHitResultTypes.MSAlign:
                case PeptideHitResultTypes.MODa:
                case PeptideHitResultTypes.MODPlus:
                case PeptideHitResultTypes.MSPathFinder:
                case PeptideHitResultTypes.TopPIC:

                    if (inputFileName.EndsWith("_fht", StringComparison.OrdinalIgnoreCase) ||
                        inputFileName.EndsWith("_syn", StringComparison.OrdinalIgnoreCase))
                    {
                        datasetName = inputFileName.Substring(0, inputFileName.Length - 4);

                        if (resultType == PeptideHitResultTypes.Inspect)
                        {
                            if (datasetName.EndsWith("_inspect", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_inspect".Length);
                            }
                        }
                        else if (resultType == PeptideHitResultTypes.MSGFPlus)
                        {
                            if (datasetName.EndsWith("_msgfplus", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_msgfplus".Length);
                            }
                            else if (datasetName.EndsWith("_msgfdb", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_msgfdb".Length);
                            }
                        }
                        else if (resultType == PeptideHitResultTypes.MSAlign)
                        {
                            if (datasetName.EndsWith("_msalign", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_msalign".Length);
                            }
                        }
                        else if (resultType == PeptideHitResultTypes.MODa)
                        {
                            if (datasetName.EndsWith("_moda", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_moda".Length);
                            }
                        }
                        else if (resultType == PeptideHitResultTypes.MODPlus)
                        {
                            if (datasetName.EndsWith("_modp", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_modp".Length);
                            }
                        }
                        else if (resultType == PeptideHitResultTypes.MSPathFinder)
                        {
                            if (datasetName.EndsWith("_mspath", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_mspath".Length);
                            }
                        }
                        else if (resultType == PeptideHitResultTypes.TopPIC)
                        {
                            if (datasetName.EndsWith("_toppic", StringComparison.OrdinalIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_toppic".Length);
                            }
                        }
                    }
                    break;

                case PeptideHitResultTypes.XTandem:
                    if (inputFileName.EndsWith("_xt", StringComparison.OrdinalIgnoreCase))
                    {
                        datasetName = inputFileName.Substring(0, inputFileName.Length - 3);
                    }
                    break;
            }

            if (string.IsNullOrEmpty(datasetName))
            {
                if (AutoTrimExtraSuffix(filePath, out var filePathTrimmed))
                {
                    datasetName = AutoDetermineDatasetName(filePathTrimmed, resultType);
                }
            }

            return datasetName;
        }

        /// <summary>
        /// Determine the PeptideHit result type given the input file path
        /// </summary>
        /// <param name="filePath"></param>
        public static PeptideHitResultTypes AutoDetermineResultType(string filePath)
        {
            const string LEGACY_MSGFPLUS_SUFFIX_SYN = "_msgfdb_syn.txt";
            const string LEGACY_MSGFPLUS_SUFFIX_FHT = "_msgfdb_fht.txt";

            var resultType = PeptideHitResultTypes.Unknown;

            var filePathLCase = filePath.ToLower();

            var suffixesToCheck = new List<KeyValuePair<string, PeptideHitResultTypes>>();
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.XTandem, clsPHRPParserXTandem.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSGFPlus, clsPHRPParserMSGFPlus.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSGFPlus, clsPHRPParserMSGFPlus.FILENAME_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSGFPlus, LEGACY_MSGFPLUS_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSGFPlus, LEGACY_MSGFPLUS_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSAlign, clsPHRPParserMSAlign.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSAlign, clsPHRPParserMSAlign.FILENAME_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MODa, clsPHRPParserMODa.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MODa, clsPHRPParserMODa.FILENAME_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MODPlus, clsPHRPParserMODPlus.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MODPlus, clsPHRPParserMODPlus.FILENAME_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSPathFinder, clsPHRPParserMSPathFinder.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.MSPathFinder, clsPHRPParserMSPathFinder.FILENAME_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.TopPIC, clsPHRPParserTopPIC.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.TopPIC, clsPHRPParserTopPIC.FILENAME_SUFFIX_FHT);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.Inspect, clsPHRPParserInspect.FILENAME_SUFFIX_SYN);
            AddSuffixToCheck(suffixesToCheck, PeptideHitResultTypes.Inspect, clsPHRPParserInspect.FILENAME_SUFFIX_FHT);

            foreach (var item in suffixesToCheck)
            {
                if (filePathLCase.EndsWith(item.Key, StringComparison.OrdinalIgnoreCase))
                {
                    resultType = item.Value;
                    break;
                }
            }

            if (resultType != PeptideHitResultTypes.Unknown)
            {
                return resultType;
            }

            // Open the file and read the header line to determine if this is a Sequest file, Inspect file, MSGFDB, or something else
            if (!File.Exists(filePath))
            {
                // File doesn't exist; assume MSGFPlus
                return PeptideHitResultTypes.MSGFPlus;
            }

            using (var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                if (!reader.EndOfStream)
                {
                    var headerLine = reader.ReadLine();

                    if (LineContainsValues(headerLine, clsPHRPParserInspect.DATA_COLUMN_MQScore,
                                           clsPHRPParserInspect.DATA_COLUMN_TotalPRMScore))
                    {
                        resultType = PeptideHitResultTypes.Inspect;
                    }
                    else if (LineContainsValues(headerLine, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFScore,
                                                clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFDB_SpecProb) ||
                             LineContainsValues(headerLine, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFScore,
                                                clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFPlus_SpecEValue) ||
                             LineContainsValues(headerLine, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFScore,
                                                clsPHRPParserMSGFPlus.DATA_COLUMN_DeNovoScore))
                    {
                        resultType = PeptideHitResultTypes.MSGFPlus;
                    }
                    else if (LineContainsValues(headerLine, clsPHRPParserSequest.DATA_COLUMN_XCorr, clsPHRPParserSequest.DATA_COLUMN_DelCn))
                    {
                        resultType = PeptideHitResultTypes.Sequest;
                    }
                }
            }

            if (resultType != PeptideHitResultTypes.Unknown)
                return resultType;

            if (AutoTrimExtraSuffix(filePath, out var filePathTrimmed))
            {
                resultType = AutoDetermineResultType(filePathTrimmed);
            }

            return resultType;
        }

        private static bool AutoTrimExtraSuffix(string filePath, out string filePathTrimmed)
        {
            // Check whether filePath ends in other known PHRP extensions
            var extraSuffixes = GetPHRPAuxiliaryFileSuffixes();

            foreach (var suffix in extraSuffixes)
            {
                if (filePath.EndsWith(suffix, StringComparison.OrdinalIgnoreCase))
                {
                    filePathTrimmed = filePath.Substring(0, filePath.Length - suffix.Length) + Path.GetExtension(filePath);
                    return true;
                }
            }

            filePathTrimmed = string.Empty;

            return false;
        }

        private readonly StringBuilder mNewPeptide = new StringBuilder();

        /// <summary>
        /// Look for dynamic mod symbols in the peptide sequence; replace with the corresponding mod masses
        /// Note that if the _SeqInfo.txt file is available, this function will not be used
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="peptideWithNumericMods">Peptide with numeric mods (output)</param>
        /// <param name="peptideMods">List of modified amino acids (output)</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>peptideWithNumericMods will look like R.TDM+15.9949ESALPVTVLSAEDIAK.T</remarks>
        private bool ConvertModsToNumericMods(string peptide,
            out string peptideWithNumericMods,
            out List<clsAminoAcidModInfo> peptideMods)
        {
            var residueLocInPeptide = 0;

            peptideMods = new List<clsAminoAcidModInfo>();
            peptideWithNumericMods = string.Empty;

            try
            {
                if (mDynamicMods.Count == 0 && mStaticMods.Count == 0)
                {
                    // No mods are defined; simply update peptideWithNumericMods to be peptide
                    peptideWithNumericMods = peptide;
                    return true;
                }

                mNewPeptide.Length = 0;
                var peptideLength = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(peptide, true).Length;

                var indexStart = 0;
                var indexEnd = peptide.Length - 1;

                if (peptide.Length >= 4)
                {
                    if (peptide[1] == '.')
                    {
                        // Peptide is of the form R.HRDTGILDSIGR.F
                        // Skip the first two characters
                        indexStart = 2;
                    }

                    if (peptide[peptide.Length - 2] == '.')
                    {
                        // Peptide is of the form R.HRDTGILDSIGR.F
                        // Skip the last two characters
                        indexEnd = peptide.Length - 3;
                    }
                }

                var index = 0;
                var mostRecentResidue = '.';
                while (index < peptide.Length)
                {
                    if (index < indexStart || index > indexEnd)
                    {
                        // We're before or after the primary peptide sequence; simply append the character
                        mNewPeptide.Append(peptide[index]);
                    }
                    else
                    {
                        var residueTerminusState = clsAminoAcidModInfo.ResidueTerminusStateConstants.None;

                        if (IsLetterAtoZ(peptide[index]))
                        {
                            mostRecentResidue = peptide[index];
                            residueLocInPeptide++;

                            if (residueLocInPeptide == 1)
                            {
                                residueTerminusState = clsAminoAcidModInfo.ResidueTerminusStateConstants.PeptideNTerminus;
                            }
                            else if (residueLocInPeptide == peptideLength)
                            {
                                residueTerminusState = clsAminoAcidModInfo.ResidueTerminusStateConstants.PeptideCTerminus;
                            }
                            else
                            {
                                residueTerminusState = clsAminoAcidModInfo.ResidueTerminusStateConstants.None;
                            }

                            // Character is a letter; append it
                            mNewPeptide.Append(mostRecentResidue);

                            if (mStaticMods.Count > 0)
                            {
                                // See if it is present in mStaticMods (this is a case-sensitive search)
                                AddStaticModIfPresent(mStaticMods, mostRecentResidue, residueLocInPeptide, residueTerminusState, mNewPeptide, peptideMods);

                                if (index == indexStart && mStaticMods.Count > 0)
                                {
                                    // We're at the N-terminus of the peptide
                                    // Possibly add a static N-terminal peptide mod (for example, iTRAQ8, which is 304.2022 Da)
                                    AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS, residueLocInPeptide, clsAminoAcidModInfo.ResidueTerminusStateConstants.PeptideNTerminus, mNewPeptide, peptideMods);

                                    if (peptide.StartsWith(PROTEIN_TERMINUS_SYMBOL_PHRP.ToString()))
                                    {
                                        // We're at the N-terminus of the protein
                                        // Possibly add a static N-terminal protein mod
                                        AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS, residueLocInPeptide, clsAminoAcidModInfo.ResidueTerminusStateConstants.ProteinNTerminus, mNewPeptide, peptideMods);
                                    }
                                }
                            }
                        }
                        else
                        {
                            // Not a letter; see if it is present in mDynamicMods
                            AddDynamicModIfPresent(mDynamicMods, mostRecentResidue, peptide[index], residueLocInPeptide, residueTerminusState, mNewPeptide, peptideMods);
                        }

                        if (index == indexEnd && mStaticMods.Count > 0)
                        {
                            // Possibly add a static C-terminal peptide mod
                            AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS, residueLocInPeptide, clsAminoAcidModInfo.ResidueTerminusStateConstants.PeptideCTerminus, mNewPeptide, peptideMods);

                            if (peptide.EndsWith(PROTEIN_TERMINUS_SYMBOL_PHRP.ToString()))
                            {
                                // We're at the C-terminus of the protein
                                // Possibly add a static C-terminal protein mod
                                AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS, residueLocInPeptide, clsAminoAcidModInfo.ResidueTerminusStateConstants.ProteinCTerminus, mNewPeptide, peptideMods);
                            }
                        }
                    }
                    index++;
                }

                peptideWithNumericMods = mNewPeptide.ToString();
            }
            catch (Exception ex)
            {
                HandleException("Error adding dynamic and static mod masses to peptide " + peptide, ex);
                return false;
            }

            return true;
        }

        private void AddDynamicModIfPresent(
            IReadOnlyDictionary<char, clsModificationDefinition> mods,
            char residue,
            char modSymbol,
            int ResidueLocInPeptide,
            clsAminoAcidModInfo.ResidueTerminusStateConstants ResidueTerminusState,
            StringBuilder sbNewPeptide,
            ICollection<clsAminoAcidModInfo> peptideMods)
        {
            if (mods.TryGetValue(modSymbol, out var modDef))
            {
                // Mod mass found for dynamic mod symbol; append the mod
                sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(modDef.ModificationMass, 4));
                peptideMods.Add(new clsAminoAcidModInfo(residue, ResidueLocInPeptide, ResidueTerminusState, modDef));
            }
        }

        private void AddStaticModIfPresent(
            IReadOnlyDictionary<string, List<clsModificationDefinition>> mods,
            char residue,
            int ResidueLocInPeptide,
            clsAminoAcidModInfo.ResidueTerminusStateConstants ResidueTerminusState,
            StringBuilder sbNewPeptide,
            ICollection<clsAminoAcidModInfo> peptideMods)
        {
            if (mods.TryGetValue(residue.ToString(), out var modDefs))
            {
                // Static mod applies to this residue; append the mod (add a plus sign if it doesn't start with a minus sign)

                foreach (var modDef in modDefs)
                {
                    sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(modDef.ModificationMass, 4));
                    peptideMods.Add(new clsAminoAcidModInfo(residue, ResidueLocInPeptide, ResidueTerminusState, modDef));
                }
            }
        }

        /// <summary>
        /// Find the ModSummary file for the given input file
        /// </summary>
        /// <param name="peptideHitResultType">PHRP Result Type of the input file</param>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputDirectoryPath">Input directory (can be an empty string if inputFileName is a full path)</param>
        /// <param name="inputFileName">Name or path of the input file, e.g. Dataset_msgfplus_syn.txt or Dataset_syn.txt</param>
        /// <param name="modSummaryFileNamePreferred">Output: preferred mod summary filename (based on whether a _syn.txt or _fht.txt file is present)</param>
        public static string FindModSummaryFile(
            PeptideHitResultTypes peptideHitResultType,
            string datasetName,
            string inputDirectoryPath,
            string inputFileName,
            out string modSummaryFileNamePreferred)
        {
            var modSummaryFileName = GetPHRPModSummaryFileName(peptideHitResultType, datasetName);
            if (string.IsNullOrEmpty(modSummaryFileName))
            {
                modSummaryFileNamePreferred = string.Empty;
                return string.Empty;
            }

            return FindModSummaryFile(inputDirectoryPath, inputFileName, modSummaryFileName, out modSummaryFileNamePreferred);
        }

        /// <summary>
        /// Find the ModSummary file for the given input file
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory (can be an empty string if inputFileName is a full path)</param>
        /// <param name="inputFileName">Name or path of the input file, e.g. Dataset_msgfplus_syn.txt or Dataset_syn.txt</param>
        /// <param name="modSummaryFileName">Expected mod summary filename</param>
        /// <param name="modSummaryFileNamePreferred">Output: preferred mod summary filename (based on whether a _syn.txt or _fht.txt file is present)</param>
        /// <returns>Mod summary file path if found; otherwise, an empty string</returns>
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
        public static string FindPHRPFile(
            string inputDirectoryPath,
            string inputFileName,
            string fileNameToFind,
            out string preferredName)
        {
            preferredName = string.Empty;

            try
            {
                if (string.IsNullOrEmpty(fileNameToFind))
                {
                    return string.Empty;
                }

                if (string.IsNullOrWhiteSpace(inputDirectoryPath))
                {
                    var inputFile = new FileInfo(inputFileName);
                    inputDirectoryPath = inputFile.DirectoryName;
                }

                var namesToFind = new List<string>();

                fileNameToFind = AutoSwitchToLegacyMSGFDBIfRequired(fileNameToFind, inputFileName);
                preferredName = AutoSwitchToFHTIfRequired(fileNameToFind, inputFileName);

                namesToFind.Add(preferredName);
                namesToFind.Add(fileNameToFind);

                if (preferredName.Contains("_msgfplus_"))
                    namesToFind.Add(preferredName.Replace("_msgfplus_", "_"));

                if (fileNameToFind.Contains("_msgfplus_"))
                    namesToFind.Add(fileNameToFind.Replace("_msgfplus_", "_"));

                foreach (var nameToFind in namesToFind)
                {
                    if (string.IsNullOrWhiteSpace(inputDirectoryPath))
                    {
                        if (File.Exists(nameToFind))
                            return nameToFind;

                        continue;
                    }

                    var fileToFind = new FileInfo(Path.Combine(inputDirectoryPath, nameToFind));
                    if (fileToFind.Exists)
                    {
                        return fileToFind.FullName;
                    }
                }
            }
            catch (Exception ex)
            {
                PRISM.ConsoleMsgUtils.ShowWarning("Exception in FindModSummaryFile: {0}", ex.Message);
            }

            return string.Empty;
        }

        /// <summary>
        /// Determines the collision mode using the Scan Type name
        /// </summary>
        /// <param name="scanTypeName"></param>
        private string GetCollisionMode(string scanTypeName)
        {
            // Typical scan types, along with DMS usage count as of 4/26/2012

            // Scan Type     Usage Count
            // CID-HMSn        2,800 spectra
            // CID-MSn        68,335
            // CID-SRM       258,479
            // ETD-HMSn          974
            // ETD-MSn           736
            // GC-MS           2,646
            // HCD-HMSn        6,767
            // HMS            73,431
            // HMSn               82
            // MRM_Full_NL         1
            // MS             23,201
            // MSn             5,229
            // PQD-HMSn            1
            // PQD-MSn            51
            // Q1MS              216
            // Q3MS              363
            // SA_ETD-HMSn       368
            // SA_ETD-MSn      2,925
            // SIM ms            306
            // SRM            95,207
            // Zoom-MS            29

            var collisionMode = scanTypeName.ToUpper();
            if (collisionMode.StartsWith("SA_"))
            {
                collisionMode = collisionMode.Substring(3);
            }

            if (collisionMode.StartsWith("CID"))
            {
                return "CID";
            }

            if (collisionMode.StartsWith("ETD"))
            {
                return "ETD";
            }

            if (collisionMode.StartsWith("HCD"))
            {
                return "HCD";
            }

            if (collisionMode.StartsWith("PQD"))
            {
                return "PQD";
            }

            return collisionMode.StartsWith("SID") ? "SID" : string.Empty;
        }

        private static string GetLegacyMSGFPlusName(string msgfPlusName)
        {
            return msgfPlusName.Replace("_msgfplus_", "_msgfdb_");
        }

        /// <summary>
        /// Returns the filename of the MSGF file that corresponds to synopsisOrFirstHitsFileName
        /// </summary>
        /// <param name="synopsisOrFirstHitsFileName">Filename (or full path) to the synopsis or first-hits file</param>
        public static string GetMSGFFileName(string synopsisOrFirstHitsFileName)
        {
            return Path.GetFileNameWithoutExtension(synopsisOrFirstHitsFileName) + MSGF_RESULT_FILENAME_SUFFIX;
        }

        /// <summary>
        /// Get the peptide hit result type for the given result type name
        /// </summary>
        /// <param name="ResultTypeName"></param>
        public static PeptideHitResultTypes GetPeptideHitResultType(string ResultTypeName)
        {
            switch (ResultTypeName.ToLower())
            {
                case "peptide_hit":
                    return PeptideHitResultTypes.Sequest;

                case "xt_peptide_hit":
                    return PeptideHitResultTypes.XTandem;

                case "in_peptide_hit":
                    return PeptideHitResultTypes.Inspect;

                case "msg_peptide_hit":
                    return PeptideHitResultTypes.MSGFPlus;

                case "msa_peptide_hit":
                    return PeptideHitResultTypes.MSAlign;

                case "moda_peptide_hit":
                    return PeptideHitResultTypes.MODa;

                case "modplus_peptide_hit":
                    return PeptideHitResultTypes.MODPlus;

                case "msp_peptide_hit":
                    return PeptideHitResultTypes.MSPathFinder;

                case "tpc_peptide_hit":
                    return PeptideHitResultTypes.TopPIC;

                default:
                    return PeptideHitResultTypes.Unknown;
            }
        }

        /// <summary>
        /// Get the list of auxiliary file suffixes
        /// </summary>
        public static List<string> GetPHRPAuxiliaryFileSuffixes()
        {
            var auxSuffixes = new List<string>
            {
                "_ResultToSeqMap.txt",
                "_SeqToProteinMap.txt",
                "_SeqInfo.txt",
                "_MSGF.txt",
                "_peptides.txt",
                "_ProteinMods.txt",
                "_ModDetails.txt",
                "_ModSummary.txt"
            };

            return auxSuffixes;
        }

        private static clsPHRPParser mCachedParser;
        private static PeptideHitResultTypes mCachedParserType = PeptideHitResultTypes.Unknown;
        private static string mCachedDataset = string.Empty;

        private static clsPHRPParser GetPHRPFileFreeParser(PeptideHitResultTypes resultType, string datasetName)
        {
            if (mCachedParserType != PeptideHitResultTypes.Unknown && mCachedParserType == resultType && mCachedDataset == datasetName)
            {
                return mCachedParser;
            }

            switch (resultType)
            {
                case PeptideHitResultTypes.Sequest:
                    mCachedParser = new clsPHRPParserSequest(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.XTandem:
                    mCachedParser = new clsPHRPParserXTandem(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.Inspect:
                    mCachedParser = new clsPHRPParserInspect(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.MSGFPlus:
                    mCachedParser = new clsPHRPParserMSGFPlus(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.MSAlign:
                    mCachedParser = new clsPHRPParserMSAlign(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.MODa:
                    mCachedParser = new clsPHRPParserMODa(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.MODPlus:
                    mCachedParser = new clsPHRPParserMODPlus(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.MSPathFinder:
                    mCachedParser = new clsPHRPParserMSPathFinder(datasetName, string.Empty);
                    break;

                case PeptideHitResultTypes.TopPIC:
                    mCachedParser = new clsPHRPParserTopPIC(datasetName, string.Empty);
                    break;

                default:
                    throw new Exception("Unsupported PeptideHitResultTypes value: " + resultType);
            }

            mCachedParserType = resultType;
            mCachedDataset = string.Copy(datasetName);

            return mCachedParser;
        }

        /// <summary>
        /// Returns the default first-hits file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var phrpResultsFileName = parser.PHRPFirstHitsFileName;

            return phrpResultsFileName;
        }

        /// <summary>
        /// Returns the default ModSummary file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var phrpModSummaryFileName = parser.PHRPModSummaryFileName;

            return phrpModSummaryFileName;
        }

        /// <summary>
        /// Returns the default PepToProtMap file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var pepToProteinMapFileName = parser.PHRPPepToProteinMapFileName;

            return pepToProteinMapFileName;
        }

        /// <summary>
        /// Returns the default ProteinMods file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var proteinModsFileName = parser.PHRPProteinModsFileName;

            return proteinModsFileName;
        }

        /// <summary>
        /// Returns the default Synopsis file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var phrpResultsFileName = parser.PHRPSynopsisFileName;

            return phrpResultsFileName;
        }

        /// <summary>
        /// Returns the default ResultToSeq Map file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var resultToSeqMapFilename = parser.PHRPResultToSeqMapFileName;

            return resultToSeqMapFilename;
        }

        /// <summary>
        /// Returns the default SeqInfo file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var seqInfoFilename = parser.PHRPSeqInfoFileName;

            return seqInfoFilename;
        }

        /// <summary>
        /// Returns the default SeqToProtein Map file name for the given PeptideHit result type
        /// </summary>
        /// <param name="resultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(PeptideHitResultTypes resultType, string datasetName)
        {
            var parser = GetPHRPFileFreeParser(resultType, datasetName);
            var seqToProteinMapFileName = parser.PHRPSeqToProteinMapFileName;

            return seqToProteinMapFileName;
        }

        /// <summary>
        /// Get the ScanStats filename for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetScanStatsFilename(string datasetName)
        {
            return datasetName + SCAN_STATS_FILENAME_SUFFIX;
        }

        /// <summary>
        /// Get the extended ScanStats filename for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetExtendedScanStatsFilename(string datasetName)
        {
            return datasetName + EXTENDED_SCAN_STATS_FILENAME_SUFFIX;
        }

        /// <summary>
        /// Get the tool version info filename for the given analysis tool
        /// </summary>
        /// <param name="resultType"></param>
        /// <returns>Filename</returns>
        public static string GetToolVersionInfoFilename(PeptideHitResultTypes resultType)
        {
            switch (resultType)
            {
                case PeptideHitResultTypes.Sequest:
                    return "Tool_Version_Info_Sequest.txt";

                case PeptideHitResultTypes.XTandem:
                    return "Tool_Version_Info_XTandem.txt";

                case PeptideHitResultTypes.Inspect:
                    return "Tool_Version_Info_Inspect.txt";

                case PeptideHitResultTypes.MSGFPlus:
                    // Changed from "Tool_Version_Info_MSGFDB.txt" to "Tool_Version_Info_MSGFPlus.txt" in November 2016
                    return "Tool_Version_Info_MSGFPlus.txt";

                case PeptideHitResultTypes.MSAlign:
                    return "Tool_Version_Info_MSAlign.txt";

                case PeptideHitResultTypes.MODa:
                    return "Tool_Version_Info_MODa.txt";

                case PeptideHitResultTypes.MODPlus:
                    return "Tool_Version_Info_MODPlus.txt";

                case PeptideHitResultTypes.MSPathFinder:
                    return "Tool_Version_Info_MSPathFinder.txt";

                case PeptideHitResultTypes.TopPIC:
                    return "Tool_Version_Info_TopPIC.txt";

                default:
                    return string.Empty;
            }
        }

        private void HandleException(string baseMessage, Exception ex)
        {
            if (string.IsNullOrEmpty(baseMessage))
            {
                baseMessage = "Error";
            }

            ReportError(baseMessage + ": " + ex.Message);
        }

        private static readonly Regex RegexIsLetter = new Regex("[A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Returns true if the character is a letter between A and Z or a and z
        /// </summary>
        /// <param name="chChar">Character to examine</param>
        /// <remarks>The Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters; use this function to only return True for normal letters between A and Z</remarks>
        public static bool IsLetterAtoZ(char chChar)
        {
            return RegexIsLetter.IsMatch(chChar.ToString());
        }

        /// <summary>
        /// Examines the string to determine if it is numeric
        /// </summary>
        /// <param name="data"></param>
        /// <returns>True if a number, otherwise false</returns>
        public static bool IsNumber(string data)
        {
            try
            {
                if (double.TryParse(data, out _))
                {
                    return true;
                }

                if (int.TryParse(data, out _))
                {
                    return true;
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }

            return false;
        }

        private static bool LineContainsValues(string dataLine, params string[] valuesToFind)
        {
            var matchCount = 0;

            foreach (var item in valuesToFind)
            {
                if (dataLine.IndexOf(item, StringComparison.OrdinalIgnoreCase) > -1)
                {
                    matchCount++;
                }
            }

            return matchCount == valuesToFind.Length;
        }

        /// <summary>
        /// Returns the index of the indicated column, as tracked by columnHeaders
        /// </summary>
        /// <param name="columnName"></param>
        /// <param name="columnHeaders"></param>
        /// <returns>Column index, or -1 if not found</returns>
        public static int LookupColumnIndex(string columnName, SortedDictionary<string, int> columnHeaders)
        {
            if (columnHeaders.TryGetValue(columnName, out var colIndex))
            {
                if (colIndex >= 0)
                {
                    return colIndex;
                }
            }

            return -1;
        }

        /// <summary>
        /// Returns the index of the indicated column, as tracked by columnHeaders
        /// </summary>
        /// <param name="columnEnum"></param>
        /// <param name="columnHeaders"></param>
        /// <returns>Column index, or -1 if not found</returns>
        public static int LookupColumnIndex(Enum columnEnum, SortedDictionary<Enum, int> columnHeaders)
        {
            if (columnHeaders.TryGetValue(columnEnum, out var colIndex))
            {
                if (colIndex >= 0)
                {
                    return colIndex;
                }
            }

            return -1;
        }

        /// <summary>
        /// Returns the string stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; an empty string if the specific column name is not recognized</returns>
        public static string LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders)
        {
            return LookupColumnValue(dataColumns, columnName, columnHeaders, string.Empty);
        }

        /// <summary>
        /// Returns the string stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; valueIfMissing if the specific column name is not recognized</returns>
        public static string LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders, string valueIfMissing)
        {
            if (dataColumns != null)
            {
                var colIndex = LookupColumnIndex(columnName, columnHeaders);
                if (colIndex >= 0 && colIndex < dataColumns.Length)
                {
                    return string.IsNullOrWhiteSpace(dataColumns[colIndex]) ? string.Empty : dataColumns[colIndex];
                }
            }

            // If we get here, return valueIfMissing
            return valueIfMissing;
        }

        /// <summary>
        /// Returns the string stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; valueIfMissing if the specific column name is not recognized</returns>
        public static string LookupColumnValue(string[] dataColumns, Enum columnEnum, SortedDictionary<Enum, int> columnHeaders, string valueIfMissing)
        {
            if (dataColumns != null)
            {
                var colIndex = LookupColumnIndex(columnEnum, columnHeaders);
                if (colIndex >= 0 && colIndex < dataColumns.Length)
                {
                    return string.IsNullOrWhiteSpace(dataColumns[colIndex]) ? string.Empty : dataColumns[colIndex];
                }
            }

            // If we get here, return valueIfMissing
            return valueIfMissing;
        }

        /// <summary>
        /// Returns the value stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
        public static int LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders, int valueIfMissing)
        {
            var valueText = LookupColumnValue(dataColumns, columnName, columnHeaders, valueIfMissing.ToString());

            int.TryParse(valueText, out var value);

            return value;
        }

        /// <summary>
        /// Returns the value stored in the given named column (using columnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
        public static double LookupColumnValue(string[] dataColumns, string columnName, SortedDictionary<string, int> columnHeaders, double valueIfMissing)
        {
            var valueText = LookupColumnValue(dataColumns, columnName, columnHeaders, valueIfMissing.ToString(CultureInfo.InvariantCulture));

            double.TryParse(valueText, out var value);

            return value;
        }

        /// <summary>
        /// Updates the column name to column index mapping in columnHeaders
        /// </summary>
        /// <param name="dataColumns">Column names read from the input file</param>
        /// <param name="columnHeaders">Column mapping dictionary object to update</param>
        /// <remarks>The SortedDictionary object should be instantiated using a case-insensitive comparer, i.e. (StringComparer.OrdinalIgnoreCase)</remarks>
        public static void ParseColumnHeaders(string[] dataColumns, SortedDictionary<string, int> columnHeaders)
        {
            // Reset the column indices in columnHeaders
            if (columnHeaders.Count > 0)
            {
                var keys = new string[columnHeaders.Count];
                columnHeaders.Keys.CopyTo(keys, 0);

                foreach (var key in keys)
                {
                    columnHeaders[key] = -1;
                }
            }

            for (var index = 0; index <= dataColumns.Length - 1; index++)
            {
                if (columnHeaders.ContainsKey(dataColumns[index]))
                {
                    // Update the index associated with this column name
                    columnHeaders[dataColumns[index]] = index;
                }
            }
        }

        /// <summary>
        /// Reads the next line from a synopsis file or first hits file
        /// </summary>
        /// <returns>True if a line was read, false if not more data is available</returns>
        /// <remarks>When FastReadMode is True, you should call FinalizeCurrentPSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
        public bool MoveNext()
        {
            var lineIn = string.Empty;

            var success = false;
            var usingCachedPSM = false;

            if (mCachedLineAvailable)
            {
                lineIn = mCachedLine;
                mCachedLineAvailable = false;
                mPSMCurrent = mCachedPSM;
                usingCachedPSM = true;
                success = true;
                mPSMCurrentFinalized = false;
                mExtendedScanStatsValid = false;
            }
            else if (!mSourceFile.EndOfStream)
            {
                lineIn = mSourceFile.ReadLine();
                mSourceFileLinesRead++;
                success = true;
            }
            else
            {
                CanRead = false;
            }

            if (!success || string.IsNullOrEmpty(lineIn))
            {
                return false;
            }

            if (!mHeaderLineParsed)
            {
                var splitLine = lineIn.Split('\t');
                if (!IsNumber(splitLine[0]))
                {
                    // Parse the header line to confirm the column ordering
                    PHRPParser.ParseColumnHeaders(splitLine);

                    mHeaderLineParsed = true;
                    return MoveNext();
                }

                mHeaderLineParsed = true;
            }

            if (!usingCachedPSM)
            {
                mPSMCurrent = null;

                success = PHRPParser.ParsePHRPDataLine(lineIn, mSourceFileLinesRead, out mPSMCurrent, mFastReadMode);

                if (mPSMCurrent == null)
                    mPSMCurrent = new clsPSM();

                mPSMCurrentFinalized = false;
                mExtendedScanStatsValid = false;
            }

            if (!success)
            {
                return false;
            }

            // The PHRPParser will update .PeptideWithNumericMods if the _SeqInfo.txt file is loaded
            // If it wasn't loaded, this class can update .PeptideWithNumericMods and .PeptideMods
            // by inferring the mods using mDynamicMods and mStaticMods (which were populated using the PHRP ModSummary file)
            if (!mFastReadMode && mStartupOptions.LoadModsAndSeqInfo && string.IsNullOrEmpty(mPSMCurrent.PeptideWithNumericMods))
            {
                MarkupPeptideWithMods();
            }

            var scanStatsValid = TryGetScanStats(mPSMCurrent.ScanNumber, out var scanStatsInfo);
            mExtendedScanStatsValid = TryGetExtendedScanStats(mPSMCurrent.ScanNumber, out mExtendedScanStatsInfo);

            if (scanStatsValid)
            {
                // Update the elution time
                mPSMCurrent.ElutionTimeMinutes = scanStatsInfo.ScanTimeMinutes;
            }

            if (string.IsNullOrEmpty(mPSMCurrent.CollisionMode) || mPSMCurrent.CollisionMode == clsPSM.UNKNOWN_COLLISION_MODE)
            {
                // Determine the ScanTypeName using the ScanStats or ExtendedScanStats info
                if (scanStatsValid && !string.IsNullOrEmpty(scanStatsInfo.ScanTypeName))
                {
                    mPSMCurrent.CollisionMode = GetCollisionMode(scanStatsInfo.ScanTypeName);
                }

                if (string.IsNullOrEmpty(mPSMCurrent.CollisionMode) && mExtendedScanStatsValid && mExtendedScanStatsInfo != null)
                {
                    // Scan type still not determined, but Extended Scan Stats data is available
                    if (!string.IsNullOrEmpty(mExtendedScanStatsInfo.CollisionMode))
                    {
                        // Check for Collision mode being "0"
                        // This is often the case for the first scan in a Thermo .Raw file
                        if (mExtendedScanStatsInfo.CollisionMode != "0")
                        {
                            mPSMCurrent.CollisionMode = mExtendedScanStatsInfo.CollisionMode.ToUpper();
                        }
                    }
                }
            }

            if (!mFastReadMode)
            {
                if (PHRPParser.PeptideHitResultType == PeptideHitResultTypes.Sequest || PHRPParser.PeptideHitResultType == PeptideHitResultTypes.XTandem)
                {
                    ComputePrecursorNeutralMass();
                }
            }

            if (mMSGFCachedResults?.Count > 0)
            {
                if (mMSGFCachedResults.TryGetValue(mPSMCurrent.ResultID, out var specEValueText))
                {
                    mPSMCurrent.MSGFSpecEValue = specEValueText;
                    if (specEValueText.Length > 12)
                    {
                        // Attempt to shorten the SpecEValue
                        if (double.TryParse(specEValueText, out var specEValue))
                        {
                            mPSMCurrent.MSGFSpecEValue = specEValue.ToString("0.00000E-00");
                        }
                    }
                }
            }

            if (SkipDuplicatePSMs)
            {
                // Read the next line and check whether it's the same hit, but a different protein
                var readNext = true;
                while (readNext && !mSourceFile.EndOfStream)
                {
                    lineIn = mSourceFile.ReadLine();
                    mSourceFileLinesRead++;

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    PHRPParser.ParsePHRPDataLine(lineIn, mSourceFileLinesRead, out var newPSM, mFastReadMode);

                    // Check for duplicate lines
                    // If this line is a duplicate of the previous line, skip it
                    // This happens in Sequest _syn.txt files where the line is repeated for all protein matches
                    // It can also happen in MS-GF+ results, though the prefix and suffix residues could differ for the same peptide, depending on the protein context

                    var isDuplicate = false;

                    if (mPSMCurrent.ScanNumber == newPSM.ScanNumber && mPSMCurrent.Charge == newPSM.Charge)
                    {
                        if (string.Equals(mPSMCurrent.Peptide, newPSM.Peptide))
                        {
                            isDuplicate = true;
                        }
                        else
                        {
                            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(mPSMCurrent.Peptide, out var peptide1, out _, out _) &&
                                clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(newPSM.Peptide, out var peptide2, out _, out _))
                            {
                                if (string.Equals(peptide1, peptide2))
                                {
                                    isDuplicate = true;
                                }
                            }
                        }
                    }

                    if (isDuplicate)
                    {
                        // Update the protein list
                        var additionalProteins = newPSM.Proteins.Except(mPSMCurrent.Proteins, StringComparer.OrdinalIgnoreCase).ToList();
                        if (additionalProteins.Count > 0)
                        {
                            foreach (var item in additionalProteins)
                                mPSMCurrent.AddProtein(item);
                        }
                    }
                    else
                    {
                        readNext = false;
                        mCachedLine = string.Copy(lineIn);
                        mCachedLineAvailable = true;
                        mCachedPSM = newPSM;
                    }
                }
            }

            return true;
        }

        private void ComputePrecursorNeutralMass()
        {
            if (!mExtendedScanStatsValid || mExtendedScanStatsInfo == null)
                return;

            double monoisotopicPrecursorMass = 0;

            // Try to extract out the precursor m/z value from the "Scan Filter Text" field

            if (ExtractParentIonMzFromFilterText(mExtendedScanStatsInfo.ScanFilterText, out var parentIonMZ))
            {
                if (parentIonMZ > 0)
                {
                    monoisotopicPrecursorMass = mPeptideMassCalculator.ConvoluteMass(parentIonMZ, mPSMCurrent.Charge, 0);
                }
            }

            if (Math.Abs(monoisotopicPrecursorMass) < double.Epsilon)
            {
                if (mExtendedScanStatsInfo.MonoisotopicMZ > 0)
                {
                    // Determine the precursor m/z value using the Monoisotopic m/z value reported by the instrument
                    monoisotopicPrecursorMass = mPeptideMassCalculator.ConvoluteMass(mExtendedScanStatsInfo.MonoisotopicMZ, mPSMCurrent.Charge, 0);
                }
            }

            if (monoisotopicPrecursorMass > 0)
            {
                mPSMCurrent.PrecursorNeutralMass = monoisotopicPrecursorMass;
            }
        }

        /// <summary>
        /// This function extracts the Parent Ion m/z from the filter string
        /// </summary>
        /// <param name="filterText"></param>
        /// <param name="parentIonMz"></param>
        /// <returns>True if parsing successful</returns>
        /// <remarks>The original version of this code is in ThermoRawFileReader.XRawFileIO.ExtractParentIonMZFromFilterText(string, out double)</remarks>
        private bool ExtractParentIonMzFromFilterText(string filterText, out double parentIonMz)
        {
            Regex matcher;

            if (filterText.IndexOf("msx", StringComparison.OrdinalIgnoreCase) >= 0)
            {
                matcher = mFindParentIonOnlyMsx;
            }
            else
            {
                matcher = mFindParentIonOnlyNonMsx;
            }

            var match = matcher.Match(filterText);

            if (match.Success)
            {
                var parentIonMzText = match.Groups["ParentMZ"].Value;

                var success = double.TryParse(parentIonMzText, out parentIonMz);
                return success;
            }

            parentIonMz = 0;
            return false;
        }

        private void MarkupPeptideWithMods()
        {
            // Markup the peptide with the dynamic and static mods

            var success = ConvertModsToNumericMods(mPSMCurrent.Peptide.Trim(), out var peptideWithMods, out var peptideMods);
            if (success)
            {
                double totalModMass = 0;

                mPSMCurrent.PeptideWithNumericMods = peptideWithMods;
                mPSMCurrent.ClearModifiedResidues();
                foreach (var modEntry in peptideMods)
                {
                    mPSMCurrent.AddModifiedResidue(modEntry);
                    totalModMass += modEntry.ModDefinition.ModificationMass;
                }

                if (Math.Abs(mPSMCurrent.PeptideMonoisotopicMass) < double.Epsilon)
                {
                    mPSMCurrent.PeptideMonoisotopicMass = mPeptideMassCalculator.ComputeSequenceMass(mPSMCurrent.PeptideCleanSequence) + totalModMass;
                }
            }
        }

        /// <summary>
        /// When FastReadMode is True, first call MoveNext to read the peptide scores.
        /// Then, if the peptide is a peptide of interest, call this function to finalize any processing steps that were skipped.
        /// </summary>
        public void FinalizeCurrentPSM()
        {
            if (mPSMCurrentFinalized)
                return;

            // Determine the clean sequence and cleavage state, and update the Seq_ID fields
            PHRPParser.FinalizePSM(mPSMCurrent);

            MarkupPeptideWithMods();

            ComputePrecursorNeutralMass();

            mPSMCurrentFinalized = true;
        }

        private void ReadAndCacheMSGFData()
        {
            try
            {
                var msgfFilePath = GetMSGFFileName(mInputFilePath);
                msgfFilePath = Path.Combine(mInputDirectoryPath, msgfFilePath);

                ReadAndCacheMSGFData(msgfFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Exception determining MSGF file path", ex);
            }
        }

        private void ReadAndCacheMSGFData(string msgfFilePath)
        {
            try
            {
                msgfFilePath = AutoSwitchToLegacyMSGFDBIfRequired(msgfFilePath, mInputFilePath);
                msgfFilePath = AutoSwitchToFHTIfRequired(msgfFilePath, mInputFilePath);

                if (File.Exists(msgfFilePath))
                {
                    var msgfReader = new clsMSGFResultsReader();
                    mMSGFCachedResults = msgfReader.ReadMSGFData(msgfFilePath);

                    if (msgfReader.ErrorMessage.Length > 0)
                    {
                        ReportError("Error reading MSGF data: " + msgfReader.ErrorMessage);
                    }
                }
                else
                {
                    ReportWarning("MSGF file not found: " + msgfFilePath);
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading MSGF file", ex);
            }
        }

        /// <summary>
        /// Reads the data in modSummaryFilePath.  Populates dynamicMods and staticMods with the modification definitions
        /// </summary>
        /// <param name="modSummaryFilePath">Path to the PHRP Mod Summary file to read</param>
        /// <param name="dynamicMods">List with mod symbols as the key and the corresponding mod mass</param>
        /// <param name="staticMods">List with amino acid names as the key and the corresponding mod mass</param>
        /// <returns>True if success; false if an error</returns>
        private bool ReadModSummaryFile(
            string modSummaryFilePath,
            IDictionary<char, clsModificationDefinition> dynamicMods,
            IDictionary<string, List<clsModificationDefinition>> staticMods)
        {
            try
            {
                // Clear dynamicMods and staticMods (should have been instantiated by the calling function)
                dynamicMods.Clear();
                staticMods.Clear();

                if (string.IsNullOrEmpty(modSummaryFilePath))
                {
                    ReportError("ModSummaryFile path is empty; unable to continue");
                    return false;
                }

                if (!File.Exists(modSummaryFilePath))
                {
                    ReportError("ModSummary file not found: " + modSummaryFilePath);
                    return false;
                }

                ShowMessage("Reading the PHRP ModSummary file");

                var modSummaryReader = new clsPHRPModSummaryReader(modSummaryFilePath);
                var success = modSummaryReader.Success;
                if (!success)
                    return false;

                if (modSummaryReader.ModificationDefs.Count == 0)
                    return true;

                var duplicateModSymbolCounts = new Dictionary<char, int>();

                foreach (var modDef in modSummaryReader.ModificationDefs)
                {
                    var modMass = modSummaryReader.GetModificationMassAsText(modDef.MassCorrectionTag);

                    switch (modDef.ModificationType)
                    {
                        case clsModificationDefinition.ModificationTypeConstants.StaticMod:
                        case clsModificationDefinition.ModificationTypeConstants.TerminalPeptideStaticMod:
                        case clsModificationDefinition.ModificationTypeConstants.ProteinTerminusStaticMod:

                            // "S", "T", or "P"
                            // Static residue mod, peptide terminus static mod, or protein terminus static mod
                            // Note that < and > mean peptide N and C terminus (clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS and clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                            // Note that [ and ] mean protein N and C terminus (clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS and clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)

                            // This mod could apply to multiple residues, so need to process each character in targetResidues
                            foreach (var chChar in modDef.TargetResidues)
                            {
                                try
                                {
                                    if (staticMods.TryGetValue(chChar.ToString(), out var modDefs))
                                    {
                                        if (!modDefs.Contains(modDef))
                                        {
                                            // Residue is already present in staticMods; this is unusual, but we'll allow it
                                            // We'll log a warning, but continue
                                            ShowMessage("Warning: Residue '" + chChar + "' has more than one static mod defined; " +
                                                        "this is not typically used, but will be allowed");
                                            modDefs.Add(modDef);
                                        }
                                    }
                                    else
                                    {
                                        modDefs = new List<clsModificationDefinition>
                                        {
                                            modDef
                                        };
                                        staticMods.Add(chChar.ToString(), modDefs);
                                    }
                                }
                                catch (Exception ex)
                                {
                                    HandleException("Exception adding static mod for " + chChar + " with ModMass=" + modMass, ex);
                                }
                            }

                            break;
                        case clsModificationDefinition.ModificationTypeConstants.DynamicMod:
                            // Dynamic residue mod (Includes mod type "D")
                            // Note that < and > mean peptide N and C terminus (clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS and clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)

                            try
                            {
                                if (dynamicMods.ContainsKey(modDef.ModificationSymbol))
                                {
                                    // Mod symbol already present in dynamicMods; this is unexpected, but can happen with MODa, which finds lots of mods
                                    if (duplicateModSymbolCounts.TryGetValue(modDef.ModificationSymbol, out var duplicateCount))
                                    {
                                        duplicateModSymbolCounts[modDef.ModificationSymbol] = duplicateCount + 1;
                                    }
                                    else
                                    {
                                        duplicateModSymbolCounts.Add(modDef.ModificationSymbol, 1);

                                        // We'll log a warning, but continue
                                        ShowMessage("Warning: Dynamic mod symbol '" + modDef.ModificationSymbol + "' is already defined; " +
                                                    "it cannot have more than one associated mod mass (duplicate has ModMass=" + modMass + ")");
                                    }
                                }
                                else
                                {
                                    dynamicMods.Add(modDef.ModificationSymbol, modDef);
                                }
                            }
                            catch (Exception ex)
                            {
                                HandleException("Exception adding dynamic mod for " + modDef.ModificationSymbol + " with ModMass=" + modMass, ex);
                            }

                            break;

                        case clsModificationDefinition.ModificationTypeConstants.IsotopicMod:
                            // Isotopic mods are not supported by this class
                            // However, do not log a warning since these are rarely used
                            break;

                        case clsModificationDefinition.ModificationTypeConstants.UnknownType:
                            // Unknown type; just ignore it
                            break;
                    }
                }

                foreach (var item in duplicateModSymbolCounts)
                {
                    if (item.Value > 1)
                    {
                        ShowMessage(string.Format(" Duplicate count for symbol '{0}': {1}", item.Key, item.Value));
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading PHRP Mod Summary file", ex);
                return false;
            }

            return true;
        }

        private void ReadScanStatsData()
        {
            try
            {
                var scanStatsFilePath = GetScanStatsFilename(DatasetName);
                scanStatsFilePath = Path.Combine(mInputDirectoryPath, scanStatsFilePath);

                ReadScanStatsData(scanStatsFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Exception determining Scan Stats file path", ex);
            }
        }

        private void ReadScanStatsData(string scanStatsFilePath)
        {
            try
            {
                if (File.Exists(scanStatsFilePath))
                {
                    var scanStatsReader = new clsScanStatsReader();
                    mScanStats = scanStatsReader.ReadScanStatsData(scanStatsFilePath);

                    if (scanStatsReader.ErrorMessage.Length > 0)
                    {
                        ReportError("Error reading ScanStats data: " + scanStatsReader.ErrorMessage);
                    }
                }
                else
                {
                    ReportWarning("ScanStats file not found: " + scanStatsFilePath);
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading Scan Stats file", ex);
            }
        }

        private void ReadExtendedScanStatsData()
        {
            try
            {
                var extendedScanStatsFilePath = GetExtendedScanStatsFilename(DatasetName);
                extendedScanStatsFilePath = Path.Combine(mInputDirectoryPath, extendedScanStatsFilePath);

                ReadExtendedScanStatsData(extendedScanStatsFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Exception determining Scan Stats file path", ex);
            }
        }

        private void ReadExtendedScanStatsData(string extendedScanStatsFilePath)
        {
            try
            {
                if (File.Exists(extendedScanStatsFilePath))
                {
                    var extendedScanStatsReader = new clsExtendedScanStatsReader();
                    mScanStatsEx = extendedScanStatsReader.ReadExtendedScanStatsData(extendedScanStatsFilePath);

                    if (extendedScanStatsReader.ErrorMessage.Length > 0)
                    {
                        ReportError("Error reading Extended ScanStats data: " + extendedScanStatsReader.ErrorMessage);
                    }
                }
                else
                {
                    // Note: we do not need to raise a warning for MSGFDB results since the extended scan stats file isn't needed
                    if (PHRPParser.PeptideHitResultType != PeptideHitResultTypes.MSGFPlus)
                    {
                        ReportWarning("Extended ScanStats file not found: " + extendedScanStatsFilePath);
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading Extended Scan Stats file", ex);
            }
        }

        private void ReportError(string message)
        {
            ErrorMessage = message;
            if (EchoMessagesToConsole)
                Console.WriteLine(message);
            ErrorMessages.Add(message);

            OnErrorEvent(message);
        }

        private void ReportWarning(string message)
        {
            if (EchoMessagesToConsole)
                Console.WriteLine(message);
            WarningMessages.Add(message);

            OnWarningEvent(message);
        }

        private void SetLocalErrorCode(PHRPReaderErrorCodes newErrorCode, bool leaveExistingErrorCodeUnchanged = false)
        {
            if (leaveExistingErrorCodeUnchanged && mLocalErrorCode != PHRPReaderErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mLocalErrorCode = newErrorCode;
            }
        }

        private void ShowMessage(string message)
        {
            if (EchoMessagesToConsole)
                Console.WriteLine(message);

            OnStatusEvent(message);
        }

        private bool TryGetScanStats(int scanNumber,
            out clsScanStatsInfo scanStatsInfo)
        {
            if (mScanStats?.Count > 0)
            {
                if (mScanStats.TryGetValue(scanNumber, out scanStatsInfo))
                {
                    return true;
                }
            }
            scanStatsInfo = null;
            return false;
        }

        private bool TryGetExtendedScanStats(int scanNumber,
            out clsScanStatsExInfo extendedScanStatsInfo)
        {
            if (mScanStatsEx != null && mScanStats.Count > 0)
            {
                if (mScanStatsEx.TryGetValue(scanNumber, out extendedScanStatsInfo))
                {
                    return true;
                }
            }
            extendedScanStatsInfo = null;
            return false;
        }

        private bool ValidateInputFiles(string inputFilePath, ref PeptideHitResultTypes resultType, ref string modSummaryFilePath)
        {
            var inputFile = new FileInfo(inputFilePath);
            if (!inputFile.Exists)
            {
                SetLocalErrorCode(PHRPReaderErrorCodes.InvalidInputFilePath);
                ReportError("Input file not found: " + inputFilePath);
                if (inputFilePath.Contains(".."))
                {
                    ReportWarning("Absolute path: " + inputFile.DirectoryName);
                }

                return false;
            }

            // Try to auto-determine the result type if it is not specified
            if (resultType == PeptideHitResultTypes.Unknown)
            {
                resultType = AutoDetermineResultType(inputFilePath);
            }

            if (resultType == PeptideHitResultTypes.Unknown)
            {
                SetLocalErrorCode(PHRPReaderErrorCodes.InputFileFormatNotRecognized);
                ReportError("Error: Unable to auto-determine file format for " + inputFilePath);
                return false;
            }

            // Extract the dataset name from the input file path
            DatasetName = AutoDetermineDatasetName(inputFilePath, resultType);
            if (string.IsNullOrEmpty(DatasetName))
            {
                if (mStartupOptions.LoadModsAndSeqInfo || mStartupOptions.LoadMSGFResults || mStartupOptions.LoadScanStatsData)
                {
                    ReportError("Error: Unable to auto-determine the dataset name from the input file name: " + inputFilePath);
                    SetLocalErrorCode(PHRPReaderErrorCodes.InputFileFormatNotRecognized);
                    return false;
                }

                ReportWarning("Unable to auto-determine the dataset name from the input file name; this is not a critical error since not reading related files: " + inputFilePath);
            }

            if (mStartupOptions.LoadModsAndSeqInfo)
            {
                var modSummaryFileName = GetPHRPModSummaryFileName(resultType, DatasetName);
                if (string.IsNullOrEmpty(modSummaryFileName))
                {
                    ReportError(string.Format("Could not determine the ModSummaryFile name for dataset {0} " +
                                              "given the PHRP result type {1}",
                                              DatasetName, resultType.ToString()));
                    SetLocalErrorCode(PHRPReaderErrorCodes.RequiredInputFileNotFound);
                    return false;
                }

                modSummaryFilePath = FindModSummaryFile(inputFile.DirectoryName,
                                                        inputFile.Name,
                                                        modSummaryFileName,
                                                        out _);

                if (!ValidateRequiredFileExists("ModSummary file", modSummaryFilePath) && inputFile.Name.IndexOf("_fht", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    SetLocalErrorCode(PHRPReaderErrorCodes.RequiredInputFileNotFound);
                    return false;
                }
            }
            else
            {
                modSummaryFilePath = string.Empty;
            }

            return true;
        }

        private bool ValidateRequiredFileExists(string fileDescription, string filePath, bool reportErrors = true)
        {
            if (string.IsNullOrEmpty(filePath))
            {
                if (reportErrors)
                {
                    ReportError(fileDescription + " is not defined");
                }
                return false;
            }

            if (!File.Exists(filePath))
            {
                if (reportErrors)
                {
                    ReportError(fileDescription + " not found: " + filePath);
                }
                return false;
            }

            return true;
        }

        #region "IDisposable Support"

        /// <summary>
        /// Used to detect redundant calls
        /// </summary>
        private bool disposedValue;

        /// <summary>
        /// Dispose of this class
        /// </summary>
        /// <param name="disposing"></param>
        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    mSourceFile?.Close();
                }
            }
            disposedValue = true;
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
}
