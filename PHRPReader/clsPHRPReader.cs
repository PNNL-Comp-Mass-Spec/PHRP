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

namespace PHRPReader
{
    /// <summary>
    ///  This class reads a tab-delimited text file (created by the Peptide File Extractor or by PHRP)
    ///  and returns the data for each peptide hit search result
    ///
    ///  It also integrates MSGF results with the peptide hit search results
    ///  And, it integrates scan stats values (to determine elution time)
    /// </summary>
    public class clsPHRPReader : PRISM.clsEventNotifier, IDisposable
    {
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
        //  For example, 756.71 in FTMS + p NSI d Full ms3 850.70@cid35.00 756.71@cid35.00 [195.00-2000.00]
        private const string PARENTION_ONLY_NONMSX_REGEX = @"[Mm][Ss]\d*[^\[\r\n]* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*(\[[^\]\r\n]\])?";

        //  This RegEx is used to extract parent ion m/z from a filter string that does contain msx
        //  ${ParentMZ} will hold the first parent ion m/z found (the first parent ion m/z corresponds to the highest peak)
        //  For example, 636.04 in FTMS + p NSI Full msx ms2 636.04@hcd28.00 641.04@hcd28.00 654.05@hcd28.00 [88.00-1355.00]
        private const string PARENTION_ONLY_MSX_REGEX = @"[Mm][Ss]\d* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*[^\[\r\n]*(\[[^\]\r\n]+\])?";

#pragma warning disable 1591

        /// <summary>
        /// Peptide hit results type
        /// </summary>
        public enum ePeptideHitResultType
        {
            Unknown = 0,
            Sequest = 1,
            XTandem = 2,
            Inspect = 3,
            MSGFDB = 4,      // Aka MSGF+
            MSAlign = 5,
            MODa = 6,
            MODPlus = 7,
            MSPathFinder = 8
        }

        /// <summary>
        /// PHRP Reader error codes
        /// </summary>
        public enum ePHRPReaderErrorCodes
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
        private string mDatasetName;
        private string mInputFilePath;
        private string mInputFolderPath;

        private bool mSkipDuplicatePSMs;

        private readonly clsPHRPStartupOptions mStartupOptions;

        private bool mEchoMessagesToConsole;

        private bool mCanRead;
        private bool mInitialized;
        private bool mModSummaryFileLoaded;

        /// <summary>
        /// When set to true, calls to MoveNext will read the next data line, but will skip several additional processing steps for performance reasons
        /// </summary>
        /// <remarks>If the peptide is a peptide of interest, you must call FinalizeCurrentPSM after calling .MoveNext()</remarks>
        private bool mFastReadMode;

        private StreamReader mSourceFile;
        private int mSourceFileLineCount;
        private int mSourceFileLinesRead;

        private clsPHRPParser mPHRPParser;
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

        private readonly List<string> mErrorMessages;
        private readonly List<string> mWarningMessages;

        private string mErrorMessage = string.Empty;
        private ePHRPReaderErrorCodes mLocalErrorCode;

        /// <summary>
        /// RegEx to extract parent ions from filter strings that do not have Full msx
        /// </summary>
        /// <remarks>Shared (aka static) only to speed up unit tests</remarks>
        private static readonly Regex mParentIonMzMatchNonMsx = new Regex(PARENTION_ONLY_NONMSX_REGEX, RegexOptions.Compiled | RegexOptions.IgnoreCase);

        /// <summary>
        /// RegEx to extract parent ions from filter strings that have Full msx
        /// </summary>
        /// <remarks>Shared (aka static) only to speed up unit tests</remarks>
        private static readonly Regex mParentIonMzMatchMsx = new Regex(PARENTION_ONLY_MSX_REGEX, RegexOptions.Compiled | RegexOptions.IgnoreCase);

        #endregion

        #region "Properties"

        /// <summary>
        /// Returns True if the input file was successfully opened and data remains to be read
        /// </summary>
        /// <value></value>
        /// <returns>True if the file is readable</returns>
        /// <remarks></remarks>
        public bool CanRead => mCanRead;

        /// <summary>
        /// Returns the most recently loaded PSM
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPSM CurrentPSM => mPSMCurrent;

        /// <summary>
        /// Returns the most recently loaded PSM's sequence info (if available)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsSeqInfo CurrentPSMSeqInfo
        {
            get
            {
                if (mPSMCurrent == null || mPHRPParser.SeqInfo == null)
                {
                    return null;
                }

                mPHRPParser.SeqInfo.TryGetValue(mPSMCurrent.SeqID, out var oSeqInfo);
                return oSeqInfo;
            }
        }

        /// <summary>
        /// Dataset name (auto-determined based on the input filename)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string DatasetName => mDatasetName;

        /// <summary>
        /// If True, then will display messages at the console
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool EchoMessagesToConsole
        {
            get => mEchoMessagesToConsole;
            set => mEchoMessagesToConsole = value;
        }

        /// <summary>
        /// Cached error messages
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public List<string> ErrorMessages => mErrorMessages;

        /// <summary>
        /// Current error message
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ErrorMessage => mErrorMessage;

        /// <summary>
        /// Used to enable fast read mode when calling MoveNext
        /// When FastReadMode is True, you should call FinalizeCurrentPSM after calling MoveNext to populate the remaining fields if the peptide is a peptide of interest
        /// </summary>
        /// <value></value>
        /// <returns></returns>
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
        /// If True, then looks for and loads the modification definitions from the _ModSummary.txt file associated with the input file
        /// Also reads the SeqInfo and related files
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool LoadModsAndSeqInfo => mStartupOptions.LoadModsAndSeqInfo;

        /// <summary>
        /// If true, then loads the MSGF SpecProb values from the _MSGF.txt file associated with the input file
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool LoadMSGFResults => mStartupOptions.LoadMSGFResults;

        /// <summary>
        /// If True, then loads the MASIC _ScanStats.txt file
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool LoadScanStatsData => mStartupOptions.LoadScanStatsData;

        /// <summary>
        /// The maximum number of proteins that will be tracked for each PSM
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int MaxProteinsPerPSM
        {
            get => mStartupOptions.MaxProteinsPerPSM;
            set => mStartupOptions.MaxProteinsPerPSM = value;
        }

        /// <summary>
        /// Returns True if the ModSummary file was successfully loaded
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool ModSummaryFileLoaded => mModSummaryFileLoaded;

        /// <summary>
        /// Peptide hit result type; Sequest, XTandem, Inspect, or MSGFDB (aka MSGF+)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public ePeptideHitResultType PeptideHitResultType
        {
            get
            {
                if (mPHRPParser == null)
                {
                    return ePeptideHitResultType.Unknown;
                }

                return mPHRPParser.PeptideHitResultType;
            }
        }

        /// <summary>
        /// Returns a number between 0 and 100 indicating the percentage of the source file that has been read
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public float PercentComplete
        {
            get
            {
                if (mSourceFileLineCount > 0)
                {
                    return mSourceFileLinesRead / Convert.ToSingle(mSourceFileLineCount) * 100f;
                }

                return 0;
            }
        }

        /// <summary>
        /// Returns the PHRP Parser object
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPHRPParser PHRPParser => mPHRPParser;

        /// <summary>
        /// Returns the cached mapping between ResultID and SeqID
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedList<int, int> ResultToSeqMap
        {
            get
            {
                if (mPHRPParser == null)
                {
                    return new SortedList<int, int>();
                }

                return mPHRPParser.ResultToSeqMap;
            }
        }

        /// <summary>
        /// Returns the cached sequence info, where key is SeqID
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedList<int, clsSeqInfo> SeqInfo
        {
            get
            {
                if (mPHRPParser == null)
                {
                    return new SortedList<int, clsSeqInfo>();
                }

                return mPHRPParser.SeqInfo;
            }
        }

        /// <summary>
        /// Returns the cached sequence to protein map information
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedList<int, List<clsProteinInfo>> SeqToProteinMap
        {
            get
            {
                if (mPHRPParser == null)
                {
                    return new SortedList<int, List<clsProteinInfo>>();
                }

                return mPHRPParser.SeqToProteinMap;
            }
        }

        /// <summary>
        /// When True, then skips near-duplicate lines in the PHRP data file (lines with the same peptide in the same scan, but different protein names)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SkipDuplicatePSMs
        {
            get => mSkipDuplicatePSMs;
            set => mSkipDuplicatePSMs = value;
        }

        /// <summary>
        /// Cached warning messages
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public List<string> WarningMessages => mWarningMessages;

        #endregion

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
        public clsPHRPReader(string strInputFilePath)
            : this(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo: true, blnLoadMSGFResults: true, blnLoadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="eResultType">Source file PeptideHit result type</param>
        /// <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
        public clsPHRPReader(string strInputFilePath, ePeptideHitResultType eResultType)
            : this(strInputFilePath, eResultType, blnLoadModsAndSeqInfo: true, blnLoadMSGFResults: true, blnLoadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
        /// <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
        /// <remarks></remarks>
        public clsPHRPReader(string strInputFilePath, bool blnLoadModsAndSeqInfo, bool blnLoadMSGFResults)
            : this(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
        /// <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
        /// <param name="blnLoadScanStats">If True, then looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
        /// <remarks></remarks>
        public clsPHRPReader(string strInputFilePath, bool blnLoadModsAndSeqInfo, bool blnLoadMSGFResults, bool blnLoadScanStats)
            : this(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats)
        {
        }

        /// <summary>
        /// Constructor that auto-determines the PeptideHit result type based on the filename
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="oStartupOptions">Startup options</param>
        /// <remarks></remarks>
        public clsPHRPReader(string strInputFilePath, clsPHRPStartupOptions oStartupOptions)
            : this(strInputFilePath, ePeptideHitResultType.Unknown, oStartupOptions)
        {
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// ''' <param name="eResultType">Source file PeptideHit result type</param>
        /// <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
        /// <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
        /// <remarks></remarks>
        public clsPHRPReader(string strInputFilePath, ePeptideHitResultType eResultType, bool blnLoadModsAndSeqInfo, bool blnLoadMSGFResults)
            : this(strInputFilePath, eResultType, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats: false)
        {
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="eResultType">Source file PeptideHit result type</param>
        /// <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
        /// <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
        /// <param name="blnLoadScanStats">If True, then looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
        /// <remarks></remarks>
        public clsPHRPReader(string strInputFilePath, ePeptideHitResultType eResultType, bool blnLoadModsAndSeqInfo, bool blnLoadMSGFResults, bool blnLoadScanStats)
        {
            var startupOptions = new clsPHRPStartupOptions
            {
                LoadModsAndSeqInfo = blnLoadModsAndSeqInfo,
                LoadMSGFResults = blnLoadMSGFResults,
                LoadScanStatsData = blnLoadScanStats
            };

            mStartupOptions = startupOptions;

            mMSGFCachedResults = new Dictionary<int, string>();

            mDynamicMods = new SortedDictionary<char, clsModificationDefinition>();
            mStaticMods = new SortedDictionary<string, List<clsModificationDefinition>>();

            mPeptideMassCalculator = new clsPeptideMassCalculator();

            mErrorMessages = new List<string>();
            mWarningMessages = new List<string>();

            InitializeClass(strInputFilePath, eResultType);
        }

        /// <summary>
        /// Constructor where the PeptideHit result type is explicitly set
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="eResultType">Source file PeptideHit result type</param>
        /// <param name="startupOptions">Startup options</param>
        /// <remarks></remarks>
        public clsPHRPReader(string strInputFilePath, ePeptideHitResultType eResultType, clsPHRPStartupOptions startupOptions)
        {
            mStartupOptions = startupOptions ?? throw new ArgumentNullException(nameof(startupOptions));

            mMSGFCachedResults = new Dictionary<int, string>();

            mDynamicMods = new SortedDictionary<char, clsModificationDefinition>();
            mStaticMods = new SortedDictionary<string, List<clsModificationDefinition>>();

            mPeptideMassCalculator = startupOptions.PeptideMassCalculator ?? new clsPeptideMassCalculator();

            mErrorMessages = new List<string>();
            mWarningMessages = new List<string>();

            InitializeClass(strInputFilePath, eResultType);
        }

        /// <summary>
        /// Updates strFilePath to have _fht instead of _syn if strFilePath contains_syn yet strBasePHRPFileName contains _fht
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="strBasePHRPFileName"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static string AutoSwitchToFHTIfRequired(string strFilePath, string strBasePHRPFileName)
        {
            if (string.IsNullOrEmpty(strBasePHRPFileName))
            {
                return strFilePath;
            }

            var basePHRPFile = new FileInfo(strBasePHRPFileName);
            if (basePHRPFile.Name.ToLower().Contains("_fht"))
            {
                // strBasePHRPFileName is first-hits-file based

                var firstHitsFile= new FileInfo(strFilePath);
                var synIndex = firstHitsFile.Name.LastIndexOf("_syn", StringComparison.InvariantCultureIgnoreCase);
                if (synIndex > 0)
                {
                    // strFilePath is synopsis-file based
                    // Change strFilePath to contain _fht instead of _syn

                    var strFilePathFHT = firstHitsFile.Name.Substring(0, synIndex) + "_fht" + firstHitsFile.Name.Substring(synIndex + "_syn".Length);

                    if (Path.IsPathRooted(strFilePath))
                    {
                        return Path.Combine(firstHitsFile.DirectoryName, strFilePathFHT);
                    }

                    return strFilePathFHT;
                }
            }

            return strFilePath;
        }

        /// <summary>
        /// Updates strFilePath to have _msgfdb instead of _msgfplus if strBasePHRPFileName contains _msgfdb
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="strBasePHRPFileName"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static string AutoSwitchToLegacyMSGFDBIfRequired(string strFilePath, string strBasePHRPFileName)
        {
            var basePHRPFile = new FileInfo(strBasePHRPFileName);
            if (basePHRPFile.Name.ToLower().Contains("_msgfdb"))
            {
                var dataFile = new FileInfo(strFilePath);
                var charIndex = dataFile.Name.LastIndexOf("_msgfplus", StringComparison.InvariantCultureIgnoreCase);
                if (charIndex > 0)
                {
                    // strFilePath has _msgfplus but should have _msgfdb

                    var strFilePathNew = dataFile.Name.Substring(0, charIndex) + "_msgfdb" + dataFile.Name.Substring(charIndex + "_msgfplus".Length);

                    if (Path.IsPathRooted(strFilePath))
                    {
                        return Path.Combine(dataFile.DirectoryName, strFilePathNew);
                    }

                    return strFilePathNew;
                }
            }

            return strFilePath;
        }

        /// <summary>
        /// Clear any cached error messages
        /// </summary>
        /// <remarks></remarks>
        public void ClearErrors()
        {
            mErrorMessages.Clear();
            mPHRPParser?.ClearErrors();
        }

        private int CountLines(string strTextFilePath)
        {
            int intTotalLines;

            try
            {
                intTotalLines = 0;
                using (var srReader = new StreamReader(new FileStream(strTextFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srReader.EndOfStream)
                    {
                        srReader.ReadLine();
                        intTotalLines += 1;
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error counting the lines in " + Path.GetFileName(strTextFilePath) + ": " + ex.Message, ex);
            }

            return intTotalLines;
        }

        /// <summary>
        /// Clear any cached warning messages
        /// </summary>
        /// <remarks></remarks>
        public void ClearWarnings()
        {
            mWarningMessages.Clear();
            mPHRPParser?.ClearWarnings();
        }

        /// <summary>
        /// Initialize the class
        /// </summary>
        /// <param name="strInputFilePath">Input file to read</param>
        /// <param name="eResultType">Source file PeptideHit result type</param>
        /// <remarks></remarks>
        private void InitializeClass(string strInputFilePath, ePeptideHitResultType eResultType)
        {
            mInitialized = false;

            InitializeMemberVariables();

            InitializeReader(strInputFilePath, eResultType);

            mInitialized = true;
        }

        private void InitializeMemberVariables()
        {
            mDatasetName = string.Empty;
            mInputFilePath = string.Empty;
            mInputFolderPath = string.Empty;

            mCanRead = false;
            mModSummaryFileLoaded = false;

            mSkipDuplicatePSMs = true;

            mEchoMessagesToConsole = false;

            mErrorMessage = string.Empty;
            mLocalErrorCode = ePHRPReaderErrorCodes.NoError;

            mSourceFileLineCount = 0;
        }

        private void InitializeReader(string strInputFilePath, ePeptideHitResultType eResultType)
        {
            var strModSummaryFilePath = string.Empty;

            try
            {
                if (string.IsNullOrEmpty(strInputFilePath))
                {
                    ReportError("Input file name is empty");
                    SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath);
                    if (!mInitialized)
                        throw new FileNotFoundException(mErrorMessage);
                    return;
                }

                // Confirm that the source file exists
                // Make sure strInputFilePath points to a valid file
                var inputFile = new FileInfo(strInputFilePath);

                mInputFolderPath = inputFile.DirectoryName;
                mInputFilePath = inputFile.FullName;

                if (!inputFile.Exists)
                {
                    ReportError("Input file not found: " + strInputFilePath);
                    SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath);
                    if (!mInitialized)
                        throw new FileNotFoundException(mErrorMessage);
                    return;
                }

                // Note that the following populates mDatasetName
                var blnSuccess = ValidateInputFiles(strInputFilePath, ref eResultType, ref strModSummaryFilePath);
                if (!blnSuccess)
                {
                    SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound, true);
                    if (!mInitialized)
                        throw new FileNotFoundException(mErrorMessage);
                    return;
                }

                // Open the input file for reading
                // Note that this will also load the MSGFSpecEValue info and ScanStats info
                blnSuccess = InitializeParser(eResultType);

                if (blnSuccess && mStartupOptions.LoadModsAndSeqInfo)
                {
                    // Read the PHRP Mod Summary File to populate mDynamicMods and mStaticMods
                    // Note that the PHRPParser also loads the ModSummary file, and that mDynamicMods and mStaticMods are only used if the _SeqInfo.txt file is not found
                    blnSuccess = ReadModSummaryFile(strModSummaryFilePath, mDynamicMods, mStaticMods);
                    if (!blnSuccess)
                    {
                        mModSummaryFileLoaded = false;
                        blnSuccess = true;
                    }
                    else
                    {
                        mModSummaryFileLoaded = true;
                    }
                }

                if (blnSuccess && mStartupOptions.LoadMSGFResults)
                {
                    // Cache the MSGF values (if present)
                    ReadAndCacheMSGFData();
                }

                if (blnSuccess && mStartupOptions.LoadScanStatsData)
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
                    throw new Exception(mErrorMessage, ex);
            }
        }

        private bool InitializeParser(ePeptideHitResultType eResultType)
        {
            var blnSuccess = true;
            var datasetName = string.Copy(mDatasetName);

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
                // Instantiate the appropriare PHRP Parser
                switch (eResultType)
                {
                    case ePeptideHitResultType.Sequest:
                        mPHRPParser = new clsPHRPParserSequest(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.XTandem:
                        // Note that Result to Protein mapping will be auto-loaded during instantiation of mPHRPParser
                        mPHRPParser = new clsPHRPParserXTandem(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.Inspect:
                        mPHRPParser = new clsPHRPParserInspect(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.MSGFDB:
                        // MSGF+
                        mPHRPParser = new clsPHRPParserMSGFDB(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.MSAlign:
                        mPHRPParser = new clsPHRPParserMSAlign(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.MODa:
                        mPHRPParser = new clsPHRPParserMODa(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.MODPlus:
                        mPHRPParser = new clsPHRPParserMODPlus(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    case ePeptideHitResultType.MSPathFinder:
                        mPHRPParser = new clsPHRPParserMSPathFinder(datasetName, mInputFilePath, mStartupOptions);

                        break;
                    default:
                        //Should never get here; invalid result type specified
                        ReportError("Invalid PeptideHit ResultType specified: " + eResultType);
                        blnSuccess = false;
                        break;
                }

                if (!blnSuccess)
                {
                    return false;
                }

                // Attach the event handlers
                RegisterEvents(mPHRPParser);

                // Report any errors cached during instantiation of mPHRPParser
                foreach (var strMessage in mPHRPParser.ErrorMessages)
                {
                    ReportError(strMessage);
                }

                // Report any warnings cached during instantiation of mPHRPParser
                foreach (var strMessage in mPHRPParser.WarningMessages)
                {
                    ReportWarning(strMessage);
                }

                mPHRPParser.ClearErrors();
                mPHRPParser.ClearWarnings();

                // Open the data file and count the number of lines so that we can compute progress
                mSourceFileLineCount = CountLines(mInputFilePath);

                // Open the data file for reading
                mSourceFile = new StreamReader(new FileStream(mInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                mCanRead = true;

                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error in InitializeParser", ex);
                if (!mInitialized)
                    throw new Exception(mErrorMessage, ex);

                return false;
            }

        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified folder
        /// If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputFolderPath">Input folder path</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        /// <remarks></remarks>
        public static string AutoDetermineBestInputFile(string inputFolderPath)
        {
            return AutoDetermineBestInputFile(inputFolderPath, out _);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified folder
        /// If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="inputFolderPath">Input folder path</param>
        /// <param name="eMatchedResultType">Output parameter: the result type of the best result file found</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        /// <remarks></remarks>
        public static string AutoDetermineBestInputFile(string inputFolderPath,
            out ePeptideHitResultType eMatchedResultType)
        {
            // Find candidate dataset names in strInputFolderPath

            var lstDatasetNames = new SortedSet<string>(StringComparer.CurrentCultureIgnoreCase);
            var lstFileSpec = new List<string>();

            if (string.IsNullOrWhiteSpace(inputFolderPath))
            {
                throw new DirectoryNotFoundException("Input folder path is empty");
            }

            var inputFolder = new DirectoryInfo(inputFolderPath);
            if (!inputFolder.Exists)
            {
                throw new DirectoryNotFoundException("Input folder not found: " + inputFolderPath);
            }

            // MSGF+
            lstFileSpec.Add(clsPHRPParserMSGFDB.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserMSGFDB.FILENAME_SUFFIX_FHT);

            // MSGF+ prior to November 2016
            lstFileSpec.Add(GetLegacyMSGFPlusName(clsPHRPParserMSGFDB.FILENAME_SUFFIX_SYN));
            lstFileSpec.Add(GetLegacyMSGFPlusName(clsPHRPParserMSGFDB.FILENAME_SUFFIX_FHT));

            // X!Tandem (only has _xt.txt files)
            lstFileSpec.Add(clsPHRPParserXTandem.FILENAME_SUFFIX_SYN);

            // MSAlign
            lstFileSpec.Add(clsPHRPParserMSAlign.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserMSAlign.FILENAME_SUFFIX_FHT);

            // Inspect
            lstFileSpec.Add(clsPHRPParserInspect.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserInspect.FILENAME_SUFFIX_FHT);

            // MODa
            lstFileSpec.Add(clsPHRPParserMODa.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserMODa.FILENAME_SUFFIX_FHT);

            // MODPlus
            lstFileSpec.Add(clsPHRPParserMODPlus.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserMODPlus.FILENAME_SUFFIX_FHT);

            // MSPathFinder
            lstFileSpec.Add(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_FHT);

            // *****************
            // ** Important: Sequest needs to be added last since files simply end in _syn.txt or _fht.txt)
            // *****************
            // Sequest
            lstFileSpec.Add(clsPHRPParserSequest.FILENAME_SUFFIX_SYN);
            lstFileSpec.Add(clsPHRPParserSequest.FILENAME_SUFFIX_FHT);

            foreach (var strFileSpec in lstFileSpec)
            {
                foreach (var dataFile in inputFolder.GetFiles("*" + strFileSpec))
                {
                    var strDataset = dataFile.Name;

                    var intCharIndex = strDataset.ToLower().IndexOf(strFileSpec, StringComparison.Ordinal);
                    if (intCharIndex > 0)
                    {
                        strDataset = strDataset.Substring(0, intCharIndex);

                        if (!lstDatasetNames.Contains(strDataset))
                        {
                            lstDatasetNames.Add(strDataset);
                        }
                    }
                }
            }

            if (lstDatasetNames.Count == 0)
            {
                Console.WriteLine("Did not find any files matching the expected filename suffixes");
                Console.WriteLine("Looked for the following in " + inputFolderPath);
                foreach (var strFileSpec in lstFileSpec)
                {
                    Console.WriteLine("  " + strFileSpec);
                }
                eMatchedResultType = ePeptideHitResultType.Unknown;
                return string.Empty;
            }

            return AutoDetermineBestInputFile(inputFolderPath, lstDatasetNames.ToList(), out eMatchedResultType);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified folder
        /// If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="strInputFolderPath">Input folder path</param>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        /// <remarks></remarks>
        public static string AutoDetermineBestInputFile(string strInputFolderPath, string datasetName)
        {
            return AutoDetermineBestInputFile(strInputFolderPath, datasetName, out _);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified folder
        /// If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="strInputFolderPath">Input folder path</param>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="eMatchedResultType">Output parameter: the result type of the best result file found</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        /// <remarks></remarks>
        public static string AutoDetermineBestInputFile(string strInputFolderPath, string datasetName,
            out ePeptideHitResultType eMatchedResultType)
        {
            var lstDatasetNames = new List<string> {
                datasetName
            };

            return AutoDetermineBestInputFile(strInputFolderPath, lstDatasetNames, out eMatchedResultType);
        }

        /// <summary>
        /// Looks for a valid _syn.txt or _fht.txt file for the given list of datasets in the specified folder
        /// If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <param name="strInputFolderPath">Input folder path</param>
        /// <param name="lstDatasetNames">List of dataset names to search for</param>
        /// <param name="eMatchedResultType">Output parameter: the result type of the best result file found</param>
        /// <returns>The full path to the most appropriate Synopsis or First hits file</returns>
        /// <remarks></remarks>
        public static string AutoDetermineBestInputFile(string strInputFolderPath, List<string> lstDatasetNames,
            out ePeptideHitResultType eMatchedResultType)
        {
            // This list contains the standard PHRP file suffixes
            var lstAuxiliaryFileSuffixes = GetPHRPAuxiliaryFileSuffixes();

            // The key in this variable is the full path to the best Synopsis or First hits file and the value is the number of PHRP-related auxiliary files
            var kvBestSynOrFHTFile = new KeyValuePair<string, int>(string.Empty, 0);

            // Set the matched result type to Unknown for now
            eMatchedResultType = ePeptideHitResultType.Unknown;

            if (string.IsNullOrWhiteSpace(strInputFolderPath))
            {
                throw new DirectoryNotFoundException("Input folder path is empty");
            }

            var inputFolder = new DirectoryInfo(strInputFolderPath);
            if (!inputFolder.Exists)
            {
                throw new DirectoryNotFoundException("Input folder not found: " + strInputFolderPath);
            }

            if (lstDatasetNames == null || lstDatasetNames.Count == 0)
            {
                throw new ArgumentException("List lstDatasetNames cannot be empty; cannot determine the best input file");
            }

            // Construct a list of the files to search for
            // Items in this list are KeyValuePairs where the key is a filename to look for and the value is a PeptideHitResultType
            var lstFilesToFind = new List<KeyValuePair<string, ePeptideHitResultType>>();

            foreach (var strDataset in lstDatasetNames)
            {
                // MSGF+
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMSGFDB.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MSGFDB));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMSGFDB.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MSGFDB));

                // MSGF+ prior to November 2016
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(GetLegacyMSGFPlusName(clsPHRPParserMSGFDB.GetPHRPSynopsisFileName(strDataset)), ePeptideHitResultType.MSGFDB));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(GetLegacyMSGFPlusName(clsPHRPParserMSGFDB.GetPHRPFirstHitsFileName(strDataset)), ePeptideHitResultType.MSGFDB));

                // X!Tandem
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserXTandem.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.XTandem));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserXTandem.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.XTandem));

                // MSAlign
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMSAlign.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MSAlign));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMSAlign.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MSAlign));

                // MODa
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMODa.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MODa));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMODa.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MODa));

                // MODPlus
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMODPlus.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MODPlus));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMODPlus.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MODPlus));

                // MSPathFinder
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMSPathFinder.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MSPathFinder));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserMSPathFinder.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MSPathFinder));

                // Inspect
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserInspect.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.Inspect));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserInspect.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.Inspect));

                // Sequest (needs to be added last since files simply end in _syn.txt or _fht.txt)
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserSequest.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.Sequest));
                lstFilesToFind.Add(new KeyValuePair<string, ePeptideHitResultType>(clsPHRPParserSequest.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.Sequest));
            }

            foreach (var kvFileToFind in lstFilesToFind)
            {
                if (!string.IsNullOrEmpty(kvFileToFind.Key))
                {
                    var synOrFHTFile = new FileInfo(Path.Combine(inputFolder.FullName, kvFileToFind.Key));

                    if (synOrFHTFile.Exists)
                    {
                        // Match found
                        // Look for PHRP-related auxiliary files

                        var intAuxFileCount = 0;
                        var strBaseName = Path.Combine(synOrFHTFile.Directory.FullName, Path.GetFileNameWithoutExtension(synOrFHTFile.Name));

                        foreach (var strSuffix in lstAuxiliaryFileSuffixes)
                        {
                            if (File.Exists(strBaseName + strSuffix))
                            {
                                intAuxFileCount += 1;
                            }
                        }

                        if (string.IsNullOrEmpty(kvBestSynOrFHTFile.Key) || intAuxFileCount > kvBestSynOrFHTFile.Value)
                        {
                            kvBestSynOrFHTFile = new KeyValuePair<string, int>(synOrFHTFile.FullName, intAuxFileCount);
                            eMatchedResultType = kvFileToFind.Value;
                        }
                    }
                }
            }

            if (string.IsNullOrWhiteSpace(kvBestSynOrFHTFile.Key))
            {
                if (lstDatasetNames.Count == 1)
                {
                    Console.WriteLine("Could not find a Synopsis or First Hits file for dataset " + lstDatasetNames.First());
                }
                else
                {
                    Console.WriteLine("Could not find a Synopsis or First Hits file for any of the candidate datasets");
                }

                Console.WriteLine("Looked for the following files:");
                foreach (var fileName in lstFilesToFind)
                {
                    if (!string.IsNullOrWhiteSpace(fileName.Key))
                    {
                        Console.WriteLine("  " + fileName.Key);
                    }
                }
            }

            // kvBestSynOrFHTFile should now contain the PHRP result file with the most auxiliary files
            return kvBestSynOrFHTFile.Key;
        }

        /// <summary>
        /// Auto-determine the dataset name using the input file path
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <returns>Dataset name</returns>
        /// <remarks>Returns an empty string if unable to determine the dataset name</remarks>
        public static string AutoDetermineDatasetName(string strFilePath)
        {
            var eResultType = AutoDetermineResultType(strFilePath);
            return AutoDetermineDatasetName(strFilePath, eResultType);
        }

        /// <summary>
        /// Auto-determine the dataset name using the input file path and specified PeptideHit result type
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="eResultType"></param>
        /// <returns>Dataset name</returns>
        /// <remarks>Returns an empty string if unable to determine the dataset name</remarks>
        public static string AutoDetermineDatasetName(string strFilePath, ePeptideHitResultType eResultType)
        {
            var datasetName = string.Empty;

            var strInputFileName = Path.GetFileNameWithoutExtension(strFilePath);

            if (string.IsNullOrWhiteSpace(strInputFileName))
                return string.Empty;

            switch (eResultType)
            {
                case ePeptideHitResultType.Sequest:
                case ePeptideHitResultType.Inspect:
                case ePeptideHitResultType.MSGFDB:
                case ePeptideHitResultType.MSAlign:
                case ePeptideHitResultType.MODa:
                case ePeptideHitResultType.MODPlus:
                case ePeptideHitResultType.MSPathFinder:

                    if (strInputFileName.EndsWith("_fht", StringComparison.InvariantCultureIgnoreCase) ||
                        strInputFileName.EndsWith("_syn", StringComparison.InvariantCultureIgnoreCase))
                    {
                        datasetName = strInputFileName.Substring(0, strInputFileName.Length - 4);

                        if (eResultType == ePeptideHitResultType.Inspect)
                        {
                            if (datasetName.EndsWith("_inspect", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_inspect".Length);
                            }
                        }
                        else if (eResultType == ePeptideHitResultType.MSGFDB)
                        {
                            if (datasetName.EndsWith("_msgfplus", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_msgfplus".Length);
                            }
                            else if (datasetName.EndsWith("_msgfdb", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_msgfdb".Length);
                            }
                        }
                        else if (eResultType == ePeptideHitResultType.MSAlign)
                        {
                            if (datasetName.EndsWith("_msalign", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_msalign".Length);
                            }
                        }
                        else if (eResultType == ePeptideHitResultType.MODa)
                        {
                            if (datasetName.EndsWith("_moda", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_moda".Length);
                            }
                        }
                        else if (eResultType == ePeptideHitResultType.MODPlus)
                        {
                            if (datasetName.EndsWith("_modp", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_modp".Length);
                            }
                        }
                        else if (eResultType == ePeptideHitResultType.MSPathFinder)
                        {
                            if (datasetName.EndsWith("_mspath", StringComparison.InvariantCultureIgnoreCase))
                            {
                                datasetName = datasetName.Substring(0, datasetName.Length - "_mspath".Length);
                            }
                        }
                    }

                    break;
                case ePeptideHitResultType.XTandem:
                    if (strInputFileName.EndsWith("_xt", StringComparison.InvariantCultureIgnoreCase))
                    {
                        datasetName = strInputFileName.Substring(0, strInputFileName.Length - 3);
                    }

                    break;
            }

            if (string.IsNullOrEmpty(datasetName))
            {
                if (AutoTrimExtraSuffix(strFilePath, out var strFilePathTrimmed))
                {
                    datasetName = AutoDetermineDatasetName(strFilePathTrimmed, eResultType);
                }
            }

            return datasetName;
        }

        /// <summary>
        /// Determine the PeptideHit result type given the input file path
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static ePeptideHitResultType AutoDetermineResultType(string strFilePath)
        {
            const string LEGACY_MSGFPLUS_SUFFIX_SYN = "_msgfdb_syn.txt";
            const string LEGACY_MSGFPLUS_SUFFIX_FHT = "_msgfdb_fht.txt";

            var eResultType = ePeptideHitResultType.Unknown;

            var strFilePathLCase = strFilePath.ToLower();

            if (strFilePathLCase.EndsWith(clsPHRPParserXTandem.FILENAME_SUFFIX_SYN))
            {
                eResultType = ePeptideHitResultType.XTandem;
            }
            else
            {
                if (strFilePathLCase.EndsWith(clsPHRPParserMSGFDB.FILENAME_SUFFIX_SYN) || strFilePathLCase.EndsWith(clsPHRPParserMSGFDB.FILENAME_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.MSGFDB;
                }
                else if (strFilePathLCase.EndsWith(LEGACY_MSGFPLUS_SUFFIX_SYN) || strFilePathLCase.EndsWith(LEGACY_MSGFPLUS_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.MSGFDB;
                }
                else if (strFilePathLCase.EndsWith(clsPHRPParserMSAlign.FILENAME_SUFFIX_SYN) || strFilePathLCase.EndsWith(clsPHRPParserMSAlign.FILENAME_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.MSAlign;
                }
                else if (strFilePathLCase.EndsWith(clsPHRPParserMODa.FILENAME_SUFFIX_SYN) || strFilePathLCase.EndsWith(clsPHRPParserMODa.FILENAME_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.MODa;
                }
                else if (strFilePathLCase.EndsWith(clsPHRPParserMODPlus.FILENAME_SUFFIX_SYN) || strFilePathLCase.EndsWith(clsPHRPParserMODPlus.FILENAME_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.MODPlus;
                }
                else if (strFilePathLCase.EndsWith(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_SYN) || strFilePathLCase.EndsWith(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.MSPathFinder;
                }
                else if (strFilePathLCase.EndsWith(clsPHRPParserInspect.FILENAME_SUFFIX_SYN) || strFilePathLCase.EndsWith(clsPHRPParserInspect.FILENAME_SUFFIX_FHT))
                {
                    eResultType = ePeptideHitResultType.Inspect;
                }
                else
                {
                    // Open the file and read the header line to determine if this is a Sequest file, Inspect file, MSGFDB, or something else

                    if (!File.Exists(strFilePath))
                    {
                        // File doesn't exist; assume Sequest
                        eResultType = ePeptideHitResultType.Sequest;
                    }
                    else
                    {
                        using (var srInFile = new StreamReader(new FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                        {
                            if (!srInFile.EndOfStream)
                            {
                                var strHeaderLine = srInFile.ReadLine();

                                if (LineContainsValues(strHeaderLine, clsPHRPParserInspect.DATA_COLUMN_MQScore, clsPHRPParserInspect.DATA_COLUMN_TotalPRMScore))
                                {
                                    eResultType = ePeptideHitResultType.Inspect;
                                }
                                else if (LineContainsValues(strHeaderLine, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb) ||
                                         LineContainsValues(strHeaderLine, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue) ||
                                         LineContainsValues(strHeaderLine, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore))
                                {
                                    eResultType = ePeptideHitResultType.MSGFDB;
                                }
                                else if (LineContainsValues(strHeaderLine, clsPHRPParserSequest.DATA_COLUMN_XCorr, clsPHRPParserSequest.DATA_COLUMN_DelCn))
                                {
                                    eResultType = ePeptideHitResultType.Sequest;
                                }
                            }
                        }
                    }
                }
            }

            if (eResultType == ePeptideHitResultType.Unknown)
            {
                if (AutoTrimExtraSuffix(strFilePath, out var strFilePathTrimmed))
                {
                    eResultType = AutoDetermineResultType(strFilePathTrimmed);
                }
            }

            return eResultType;
        }

        private static bool AutoTrimExtraSuffix(string strFilePath, out string strFilePathTrimmed)
        {
            // Check whether strfilePathLCase ends in other known PHRP extensions
            var lstExtraSuffixes = GetPHRPAuxiliaryFileSuffixes();

            foreach (var strSuffix in lstExtraSuffixes)
            {
                if (strFilePath.EndsWith(strSuffix, StringComparison.InvariantCultureIgnoreCase))
                {
                    strFilePathTrimmed = strFilePath.Substring(0, strFilePath.Length - strSuffix.Length) + ".txt";
                    return true;
                }
            }

            strFilePathTrimmed = string.Empty;

            return false;
        }

        private readonly StringBuilder mNewPeptide = new StringBuilder();

        /// <summary>
        /// Look for dynamic mod symbols in the peptide sequence; replace with the corresponding mod masses
        /// Note that if the _SeqInfo.txt file is available, then this function will not be used
        /// </summary>
        /// <param name="strPeptide"></param>
        /// <param name="strPeptideWithNumericMods">Peptide with numeric mods (output)</param>
        /// <param name="lstPeptideMods">List of modified amino acids (output)</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>strPeptideWithNumericMods will look like R.TDM+15.9949ESALPVTVLSAEDIAK.T</remarks>
        private bool ConvertModsToNumericMods(string strPeptide,
            out string strPeptideWithNumericMods,
            out List<clsAminoAcidModInfo> lstPeptideMods)
        {
            var intResidueLocInPeptide = 0;

            lstPeptideMods = new List<clsAminoAcidModInfo>();
            strPeptideWithNumericMods = string.Empty;

            try
            {
                if (mDynamicMods.Count == 0 && mStaticMods.Count == 0)
                {
                    // No mods are defined; simply update strPeptideWithNumericMods to be strPeptide
                    strPeptideWithNumericMods = strPeptide;
                    return true;
                }

                mNewPeptide.Length = 0;
                var intPeptideLength = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strPeptide, true).Length;

                var intIndexStart = 0;
                var intIndexEnd = strPeptide.Length - 1;

                if (strPeptide.Length >= 4)
                {
                    if (strPeptide[1] == '.')
                    {
                        // Peptide is of the form R.HRDTGILDSIGR.F
                        // Skip the first two characters
                        intIndexStart = 2;
                    }

                    if (strPeptide[strPeptide.Length - 2] == '.')
                    {
                        // Peptide is of the form R.HRDTGILDSIGR.F
                        // Skip the last two characters
                        intIndexEnd = strPeptide.Length - 3;
                    }
                }

                var intIndex = 0;
                var chMostRecentResidue = '.';
                while (intIndex < strPeptide.Length)
                {
                    if (intIndex < intIndexStart || intIndex > intIndexEnd)
                    {
                        // We're before or after the primary peptide sequence; simply append the character
                        mNewPeptide.Append(strPeptide[intIndex]);
                    }
                    else
                    {
                        var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;

                        if (IsLetterAtoZ(strPeptide[intIndex]))
                        {
                            chMostRecentResidue = strPeptide[intIndex];
                            intResidueLocInPeptide += 1;

                            if (intResidueLocInPeptide == 1)
                            {
                                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                            }
                            else if (intResidueLocInPeptide == intPeptideLength)
                            {
                                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                            }
                            else
                            {
                                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                            }

                            // Character is a letter; append it
                            mNewPeptide.Append(chMostRecentResidue);

                            if (mStaticMods.Count > 0)
                            {
                                // See if it is present in mStaticMods (this is a case-sensitive search)
                                AddStaticModIfPresent(mStaticMods, chMostRecentResidue, intResidueLocInPeptide, eResidueTerminusState, mNewPeptide, lstPeptideMods);

                                if (intIndex == intIndexStart && mStaticMods.Count > 0)
                                {
                                    // We're at the N-terminus of the peptide
                                    // Possibly add a static N-terminal peptide mod (for example, iTRAQ8, which is 304.2022 Da)
                                    AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus, mNewPeptide, lstPeptideMods);

                                    if (strPeptide.StartsWith(PROTEIN_TERMINUS_SYMBOL_PHRP.ToString()))
                                    {
                                        // We're at the N-terminus of the protein
                                        // Possibly add a static N-terminal protein mod
                                        AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus, mNewPeptide, lstPeptideMods);
                                    }
                                }
                            }
                        }
                        else
                        {
                            // Not a letter; see if it is present in mDynamicMods
                            AddDynamicModIfPresent(mDynamicMods, chMostRecentResidue, strPeptide[intIndex], intResidueLocInPeptide, eResidueTerminusState, mNewPeptide, lstPeptideMods);
                        }

                        if (intIndex == intIndexEnd && mStaticMods.Count > 0)
                        {
                            // Possibly add a static C-terminal peptide mod
                            AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus, mNewPeptide, lstPeptideMods);

                            if (strPeptide.EndsWith(PROTEIN_TERMINUS_SYMBOL_PHRP.ToString()))
                            {
                                // We're at the C-terminus of the protein
                                // Possibly add a static C-terminal protein mod
                                AddStaticModIfPresent(mStaticMods, clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus, mNewPeptide, lstPeptideMods);
                            }
                        }
                    }
                    intIndex += 1;
                }

                strPeptideWithNumericMods = mNewPeptide.ToString();
            }
            catch (Exception ex)
            {
                HandleException("Error adding dynamic and static mod masses to peptide " + strPeptide, ex);
                return false;
            }

            return true;
        }

        private void AddDynamicModIfPresent(
            IReadOnlyDictionary<char, clsModificationDefinition> objMods,
            char chResidue,
            char chModSymbol,
            int ResidueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants ResidueTerminusState,
            StringBuilder sbNewPeptide,
            ICollection<clsAminoAcidModInfo> lstPeptideMods)
        {
            if (objMods.TryGetValue(chModSymbol, out var objModDef))
            {
                // Mod mass found for dynamic mod symbol; append the mod
                sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(objModDef.ModificationMass, 4));
                lstPeptideMods.Add(new clsAminoAcidModInfo(chResidue, ResidueLocInPeptide, ResidueTerminusState, objModDef));
            }
        }

        private void AddStaticModIfPresent(
            IReadOnlyDictionary<string, List<clsModificationDefinition>> objMods,
            char chResidue,
            int ResidueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants ResidueTerminusState,
            StringBuilder sbNewPeptide,
            ICollection<clsAminoAcidModInfo> lstPeptideMods)
        {
            if (objMods.TryGetValue(chResidue.ToString(), out var lstModDefs))
            {
                // Static mod applies to this residue; append the mod (add a plus sign if it doesn't start with a minus sign)

                foreach (var objModDef in lstModDefs)
                {
                    sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(objModDef.ModificationMass, 4));
                    lstPeptideMods.Add(new clsAminoAcidModInfo(chResidue, ResidueLocInPeptide, ResidueTerminusState, objModDef));
                }
            }
        }

        /// <summary>
        /// Determines the collision mode using the Scan Type name
        /// </summary>
        /// <param name="strScanTypeName"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private string GetCollisionMode(string strScanTypeName)
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

            var strCollisionMode = strScanTypeName.ToUpper();
            if (strCollisionMode.StartsWith("SA_"))
            {
                strCollisionMode = strCollisionMode.Substring(3);
            }

            if (strCollisionMode.StartsWith("CID"))
            {
                return "CID";
            }

            if (strCollisionMode.StartsWith("ETD"))
            {
                return "ETD";
            }

            if (strCollisionMode.StartsWith("HCD"))
            {
                return "HCD";
            }

            if (strCollisionMode.StartsWith("PQD"))
            {
                return "PQD";
            }

            if (strCollisionMode.StartsWith("SID"))
            {
                return "SID";
            }

            return "";
        }

        private static string GetLegacyMSGFPlusName(string msgfPlusName)
        {
            return msgfPlusName.Replace("_msgfplus_", "_msgfdb_");
        }

        /// <summary>
        /// Returns the filename of the MSGF file that corresponds to strSynopsisOrFirstHitsFileName
        /// </summary>
        /// <param name="strSynopsisOrFirstHitsFileName">Filename (or full path) to the synopsis or first-hits file</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static string GetMSGFFileName(string strSynopsisOrFirstHitsFileName)
        {
            return Path.GetFileNameWithoutExtension(strSynopsisOrFirstHitsFileName) + MSGF_RESULT_FILENAME_SUFFIX;
        }

        /// <summary>
        /// Get the peptide hit result type for the given result type name
        /// </summary>
        /// <param name="ResultTypeName"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static ePeptideHitResultType GetPeptideHitResultType(string ResultTypeName)
        {
            switch (ResultTypeName.ToLower())
            {
                //case "Peptide_Hit".ToLower():
                case "peptide_hit":

                    return ePeptideHitResultType.Sequest;
                //case "XT_Peptide_Hit".ToLower():
                case "xt_peptide_hit":

                    return ePeptideHitResultType.XTandem;
                //case "IN_Peptide_Hit".ToLower():
                case "in_peptide_hit":

                    return ePeptideHitResultType.Inspect;
                //case "MSG_Peptide_Hit".ToLower():
                case "msg_peptide_hit":

                    return ePeptideHitResultType.MSGFDB;
                //case "MSA_Peptide_Hit".ToLower():
                case "msa_peptide_hit":

                    return ePeptideHitResultType.MSAlign;
                //case "MODa_Peptide_Hit".ToLower():
                case "moda_peptide_hit":

                    return ePeptideHitResultType.MODa;
                //case "MODPlus_Peptide_Hit".ToLower():
                case "modplus_peptide_hit":

                    return ePeptideHitResultType.MODPlus;
                //case "MSP_Peptide_Hit".ToLower():
                case "msp_peptide_hit":

                    return ePeptideHitResultType.MSPathFinder;
                default:
                    return ePeptideHitResultType.Unknown;
            }
        }

        /// <summary>
        /// Get the list of auxiliary file suffixes
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public static List<string> GetPHRPAuxiliaryFileSuffixes()
        {
            var lstAuxSuffixes = new List<string>
            {
                "_ResultToSeqMap.txt",
                "_SeqToProteinMap.txt",
                "_SeqInfo.txt",
                "_MSGF.txt",
                "_ProteinMods.txt",
                "_ModDetails.txt",
                "_ModSummary.txt"
            };

            return lstAuxSuffixes;
        }

        private static clsPHRPParser oCachedParser;
        private static ePeptideHitResultType eCachedResultType = ePeptideHitResultType.Unknown;
        private static string strCachedDataset = string.Empty;

        private static clsPHRPParser GetPHRPFileFreeParser(ePeptideHitResultType eResultType, string datasetName)
        {
            if (eCachedResultType != ePeptideHitResultType.Unknown && eCachedResultType == eResultType && strCachedDataset == datasetName)
            {
                return oCachedParser;
            }

            switch (eResultType)
            {
                case ePeptideHitResultType.Sequest:
                    oCachedParser = new clsPHRPParserSequest(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.XTandem:
                    oCachedParser = new clsPHRPParserXTandem(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.Inspect:
                    oCachedParser = new clsPHRPParserInspect(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.MSGFDB:
                    oCachedParser = new clsPHRPParserMSGFDB(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.MSAlign:
                    oCachedParser = new clsPHRPParserMSAlign(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.MODa:
                    oCachedParser = new clsPHRPParserMODa(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.MODPlus:
                    oCachedParser = new clsPHRPParserMODPlus(datasetName, string.Empty);

                    break;
                case ePeptideHitResultType.MSPathFinder:
                    oCachedParser = new clsPHRPParserMSPathFinder(datasetName, string.Empty);

                    break;
                default:
                    throw new Exception("Unsupported ePeptideHitResultType value: " + eResultType);
            }

            eCachedResultType = eResultType;
            strCachedDataset = string.Copy(datasetName);

            return oCachedParser;
        }

        /// <summary>
        /// Returns the default first-hits file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strPHRPResultsFileName = oParser.PHRPFirstHitsFileName;

            return strPHRPResultsFileName;
        }

        /// <summary>
        /// Returns the default ModSummary file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strPHRPModSummaryFileName = oParser.PHRPModSummaryFileName;

            return strPHRPModSummaryFileName;
        }

        /// <summary>
        /// Returns the default PepToProtMap file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strPepToProteinMapFileName = oParser.PHRPPepToProteinMapFileName;

            return strPepToProteinMapFileName;
        }

        /// <summary>
        /// Returns the default ProteinMods file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strProteinModsFileName = oParser.PHRPProteinModsFileName;

            return strProteinModsFileName;
        }

        /// <summary>
        /// Returns the default Synopsis file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strPHRPResultsFileName = oParser.PHRPSynopsisFileName;

            return strPHRPResultsFileName;
        }

        /// <summary>
        /// Returns the default ResultToSeq Map file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strResultToSeqMapFilename = oParser.PHRPResultToSeqMapFileName;

            return strResultToSeqMapFilename;
        }

        /// <summary>
        /// Returns the default SeqInfo file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strSeqInfoFilename = oParser.PHRPSeqInfoFileName;

            return strSeqInfoFilename;
        }

        /// <summary>
        /// Returns the default SeqToProtein Map file name for the given PeptideHit result type
        /// </summary>
        /// <param name="eResultType"></param>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(ePeptideHitResultType eResultType, string datasetName)
        {
            var oParser = GetPHRPFileFreeParser(eResultType, datasetName);
            var strSeqToProteinMapFileName = oParser.PHRPSeqToProteinMapFileName;

            return strSeqToProteinMapFileName;
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
        /// <param name="eResultType"></param>
        /// <returns>Filename</returns>
        public static string GetToolVersionInfoFilename(ePeptideHitResultType eResultType)
        {
            var strToolVersionInfoFilename = string.Empty;

            switch (eResultType)
            {
                case ePeptideHitResultType.Sequest:
                    strToolVersionInfoFilename = "Tool_Version_Info_Sequest.txt";

                    break;
                case ePeptideHitResultType.XTandem:
                    strToolVersionInfoFilename = "Tool_Version_Info_XTandem.txt";

                    break;
                case ePeptideHitResultType.Inspect:
                    strToolVersionInfoFilename = "Tool_Version_Info_Inspect.txt";

                    break;
                case ePeptideHitResultType.MSGFDB:
                    // Changed from "Tool_Version_Info_MSGFDB.txt" to "Tool_Version_Info_MSGFPlus.txt" in November 2016
                    strToolVersionInfoFilename = "Tool_Version_Info_MSGFPlus.txt";

                    break;
                case ePeptideHitResultType.MSAlign:
                    strToolVersionInfoFilename = "Tool_Version_Info_MSAlign.txt";

                    break;
                case ePeptideHitResultType.MODa:
                    strToolVersionInfoFilename = "Tool_Version_Info_MODa.txt";

                    break;
                case ePeptideHitResultType.MODPlus:
                    strToolVersionInfoFilename = "Tool_Version_Info_MODPlus.txt";

                    break;
                case ePeptideHitResultType.MSPathFinder:
                    strToolVersionInfoFilename = "Tool_Version_Info_MSPathFinder.txt";

                    break;
            }

            return strToolVersionInfoFilename;
        }

        private void HandleException(string strBaseMessage, Exception ex)
        {
            if (string.IsNullOrEmpty(strBaseMessage))
            {
                strBaseMessage = "Error";
            }

            ReportError(strBaseMessage + ": " + ex.Message);
        }

        private static readonly Regex RegexIsLetter = new Regex(@"[A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Returns true if the character is a letter between A and Z or a and z
        /// </summary>
        /// <param name="chChar">Character to examine</param>
        /// <returns></returns>
        /// <remarks>The Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters; use this function to only return True for normal letters between A and Z</remarks>
        public static bool IsLetterAtoZ(char chChar)
        {
            if (RegexIsLetter.IsMatch(chChar.ToString()))
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Examines the string to determine if it is numeric
        /// </summary>
        /// <param name="strData"></param>
        /// <returns>True if a number, otherwise false</returns>
        public static bool IsNumber(string strData)
        {
            try
            {
                if (double.TryParse(strData, out _))
                {
                    return true;
                }

                if (int.TryParse(strData, out _))
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

        private static bool LineContainsValues(string strDataLine, params string[] lstValuesToFind)
        {
            var intMatchCount = 0;

            foreach (var item in lstValuesToFind)
            {
                if (strDataLine.IndexOf(item, StringComparison.CurrentCultureIgnoreCase) > -1)
                {
                    intMatchCount += 1;
                }
            }

            if (intMatchCount == lstValuesToFind.Length)
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Returns the index of the indicated column, as tracked by objColumnHeaders
        /// </summary>
        /// <param name="strColumnName"></param>
        /// <param name="objColumnHeaders"></param>
        /// <returns>Column index, or -1 if not found</returns>
        /// <remarks></remarks>
        public static int LookupColumnIndex(string strColumnName, SortedDictionary<string, int> objColumnHeaders)
        {
            if (objColumnHeaders.TryGetValue(strColumnName, out var intColIndex))
            {
                if (intColIndex >= 0)
                {
                    return intColIndex;
                }
            }

            return -1;
        }

        /// <summary>
        /// Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; an empty string if the specific column name is not recognized</returns>
        /// <remarks></remarks>
        public static string LookupColumnValue(string[] strColumns, string strColumnName, SortedDictionary<string, int> objColumnHeaders)
        {
            return LookupColumnValue(strColumns, strColumnName, objColumnHeaders, string.Empty);
        }

        /// <summary>
        /// Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The text in the specified column; strValueIfMissing if the specific column name is not recognized</returns>
        /// <remarks></remarks>
        public static string LookupColumnValue(string[] strColumns, string strColumnName, SortedDictionary<string, int> objColumnHeaders, string strValueIfMissing)
        {
            if (strColumns != null)
            {
                var intColIndex = LookupColumnIndex(strColumnName, objColumnHeaders);
                if (intColIndex >= 0 && intColIndex < strColumns.Length)
                {
                    if (string.IsNullOrWhiteSpace(strColumns[intColIndex]))
                    {
                        return string.Empty;
                    }

                    return strColumns[intColIndex];
                }
            }

            // If we get here, return strValueIfMissing
            return strValueIfMissing;
        }

        /// <summary>
        /// Returns the value stored in the given named column (using objColumnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
        /// <remarks></remarks>
        public static int LookupColumnValue(string[] strColumns, string strColumnName, SortedDictionary<string, int> objColumnHeaders, int valueIfMissing)
        {
            var strValue = LookupColumnValue(strColumns, strColumnName, objColumnHeaders, valueIfMissing.ToString());

            int.TryParse(strValue, out var intValue);

            return intValue;
        }

        /// <summary>
        /// Returns the value stored in the given named column (using objColumnHeaders to dereference column name with column index)
        /// </summary>
        /// <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
        /// <remarks></remarks>
        public static double LookupColumnValue(string[] strColumns, string strColumnName, SortedDictionary<string, int> objColumnHeaders, double ValueIfMissing)
        {
            var strValue = LookupColumnValue(strColumns, strColumnName, objColumnHeaders, ValueIfMissing.ToString(CultureInfo.InvariantCulture));

            double.TryParse(strValue, out var dblValue);

            return dblValue;
        }

        /// <summary>
        /// Updates the column name to column index mapping in objColumnHeaders
        /// </summary>
        /// <param name="strColumns">Column names read from the input file</param>
        /// <param name="objColumnHeaders">Column mapping dictionary object to update</param>
        /// <remarks>The SortedDictionary object should be instantiated using a case-insensitive comparer, i.e. (StringComparer.CurrentCultureIgnoreCase)</remarks>
        public static void ParseColumnHeaders(string[] strColumns, SortedDictionary<string, int> objColumnHeaders)
        {
            // Reset the column indices in objColumnHeaders
            if (objColumnHeaders.Count > 0)
            {
                var strKeys = new string[objColumnHeaders.Count];
                objColumnHeaders.Keys.CopyTo(strKeys, 0);

                foreach (var strKey in strKeys)
                {
                    objColumnHeaders[strKey] = -1;
                }
            }

            for (var intIndex = 0; intIndex <= strColumns.Length - 1; intIndex++)
            {
                if (objColumnHeaders.ContainsKey(strColumns[intIndex]))
                {
                    // Update the index associated with this column name
                    objColumnHeaders[strColumns[intIndex]] = intIndex;
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
            var strLineIn = string.Empty;

            var blnSuccess = false;
            var blnUsingCachedPSM = false;

            if (mCachedLineAvailable)
            {
                strLineIn = mCachedLine;
                mCachedLineAvailable = false;
                mPSMCurrent = mCachedPSM;
                blnUsingCachedPSM = true;
                blnSuccess = true;
                mPSMCurrentFinalized = false;
                mExtendedScanStatsValid = false;
            }
            else if (!mSourceFile.EndOfStream)
            {
                strLineIn = mSourceFile.ReadLine();
                mSourceFileLinesRead += 1;
                blnSuccess = true;
            }
            else
            {
                mCanRead = false;
            }

            if (!blnSuccess || string.IsNullOrEmpty(strLineIn))
            {
                return false;
            }

            if (!mHeaderLineParsed)
            {
                var strSplitLine = strLineIn.Split('\t');
                if (!IsNumber(strSplitLine[0])) {
                    // Parse the header line to confirm the column ordering
                    mPHRPParser.ParseColumnHeaders(strSplitLine);

                    mHeaderLineParsed = true;
                    return MoveNext();
                }

                mHeaderLineParsed = true;
            }

            if (!blnUsingCachedPSM) {
                mPSMCurrent = null;

                blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mSourceFileLinesRead, out mPSMCurrent, mFastReadMode);

                if (mPSMCurrent == null)
                    mPSMCurrent = new clsPSM();

                mPSMCurrentFinalized = false;
                mExtendedScanStatsValid = false;
            }

            if (!blnSuccess) {
                return false;
            }

            var blnMatchFound = true;

            // The PHRPParser will update .PeptideWithNumericMods if the _SeqInfo.txt file is loaded
            // If it wasn't loaded, then this class can update .PeptideWithNumericMods and .PeptideMods
            // by inferring the mods using mDynamicMods and mStaticMods (which were populated using the PHRP ModSummary file)
            if (!mFastReadMode && mStartupOptions.LoadModsAndSeqInfo && string.IsNullOrEmpty(mPSMCurrent.PeptideWithNumericMods)) {
                MarkupPeptideWithMods();
            }

            var blnScanStatsValid = TryGetScanStats(mPSMCurrent.ScanNumber, out var objScanStatsInfo);
            mExtendedScanStatsValid = TryGetExtendedScanStats(mPSMCurrent.ScanNumber, out mExtendedScanStatsInfo);

            if (blnScanStatsValid) {
                // Update the elution time
                mPSMCurrent.ElutionTimeMinutes = objScanStatsInfo.ScanTimeMinutes;
            }

            if (string.IsNullOrEmpty(mPSMCurrent.CollisionMode) || mPSMCurrent.CollisionMode == clsPSM.UNKNOWN_COLLISION_MODE) {
                // Determine the ScanTypeName using the the ScanStats or ExtendedScanStats info
                if (blnScanStatsValid && !string.IsNullOrEmpty(objScanStatsInfo.ScanTypeName)) {
                    mPSMCurrent.CollisionMode = GetCollisionMode(objScanStatsInfo.ScanTypeName);
                }

                if (string.IsNullOrEmpty(mPSMCurrent.CollisionMode) && mExtendedScanStatsValid && mExtendedScanStatsInfo != null) {
                    // Scan type still not determined, but Extended Scan Stats data is available
                    if (!string.IsNullOrEmpty(mExtendedScanStatsInfo.CollisionMode)) {
                        // Check for Collision mode being "0"
                        // This is often the case for the first scan in a Thermo .Raw file
                        if (mExtendedScanStatsInfo.CollisionMode != "0") {
                            mPSMCurrent.CollisionMode = mExtendedScanStatsInfo.CollisionMode.ToUpper();
                        }
                    }
                }
            }

            if (!mFastReadMode) {
                if (mPHRPParser.PeptideHitResultType == ePeptideHitResultType.Sequest || mPHRPParser.PeptideHitResultType == ePeptideHitResultType.XTandem) {
                    ComputePrecursorNeutralMass();
                }
            }

            if (mMSGFCachedResults != null && mMSGFCachedResults.Count > 0) {
                if (mMSGFCachedResults.TryGetValue(mPSMCurrent.ResultID, out var specEValueText)) {
                    mPSMCurrent.MSGFSpecEValue = specEValueText;
                    if (specEValueText.Length > 12)
                    {
                        // Attempt to shorten the SpecEValue value
                        if (double.TryParse(specEValueText, out var specEValue)) {
                            mPSMCurrent.MSGFSpecEValue = specEValue.ToString("0.00000E-00");
                        }
                    }
                }
            }

            if (mSkipDuplicatePSMs) {
                // Read the next line and check whether it's the same hit, but a different protein
                var blnReadNext = true;
                while (blnReadNext && !mSourceFile.EndOfStream) {
                    strLineIn = mSourceFile.ReadLine();
                    mSourceFileLinesRead += 1;

                    if (string.IsNullOrEmpty(strLineIn))
                        continue;

                    mPHRPParser.ParsePHRPDataLine(strLineIn, mSourceFileLinesRead, out var objNewPSM, mFastReadMode);

                    // Check for duplicate lines
                    // If this line is a duplicate of the previous line, then skip it
                    // This happens in Sequest _syn.txt files where the line is repeated for all protein matches
                    // It can also happen in MSGF+ results, though the prefix and suffix residues could differ for the same peptide, depending on the protein context

                    var isDuplicate = false;

                    if (mPSMCurrent.ScanNumber == objNewPSM.ScanNumber && mPSMCurrent.Charge == objNewPSM.Charge) {
                        if (string.Equals(mPSMCurrent.Peptide, objNewPSM.Peptide)) {
                            isDuplicate = true;
                        } else {
                            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(mPSMCurrent.Peptide, out var strPeptide1, out _, out _) &&
                                clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(objNewPSM.Peptide, out var strPeptide2, out _, out _)) {
                                if (string.Equals(strPeptide1, strPeptide2)) {
                                    isDuplicate = true;
                                }
                            }
                        }
                    }

                    if (isDuplicate) {
                        // Update the protein list
                        var addnlProteins = objNewPSM.Proteins.Except(mPSMCurrent.Proteins, StringComparer.CurrentCultureIgnoreCase).ToList();
                        if (addnlProteins.Any())
                        {
                            foreach (var item in addnlProteins)
                                mPSMCurrent.AddProtein(item);
                        }
                    } else {
                        blnReadNext = false;
                        mCachedLine = string.Copy(strLineIn);
                        mCachedLineAvailable = true;
                        mCachedPSM = objNewPSM;
                    }
                }
            }

            return blnMatchFound;
        }

        private void ComputePrecursorNeutralMass()
        {
            if (!mExtendedScanStatsValid || mExtendedScanStatsInfo == null)
                return;

            double dblMonoisotopicPrecursorMass = 0;

            // Try to extract out the precursor m/z value from the "Scan Filter Text" field

            if (ExtractParentIonMzFromFilterText(mExtendedScanStatsInfo.ScanFilterText, out var dblParentIonMZ))
            {
                if (dblParentIonMZ > 0)
                {
                    dblMonoisotopicPrecursorMass = mPeptideMassCalculator.ConvoluteMass(dblParentIonMZ, mPSMCurrent.Charge, 0);
                }
            }

            if (Math.Abs(dblMonoisotopicPrecursorMass) < double.Epsilon)
            {
                if (mExtendedScanStatsInfo.MonoisotopicMZ > 0)
                {
                    // Determine the precursor m/z value using the Monoisotopic m/z value reported by the instrument
                    dblMonoisotopicPrecursorMass = mPeptideMassCalculator.ConvoluteMass(mExtendedScanStatsInfo.MonoisotopicMZ, mPSMCurrent.Charge, 0);
                }
            }

            if (dblMonoisotopicPrecursorMass > 0)
            {
                mPSMCurrent.PrecursorNeutralMass = dblMonoisotopicPrecursorMass;
            }
        }

        /// <summary>
        /// This function extracts the Parent Ion m/z from the filter string
        /// </summary>
        /// <param name="filterText"></param>
        /// <param name="parentIonMz"></param>
        /// <returns>True if parsing successful</returns>
        /// <remarks>The original version of this code (C#) is in ThermoRawFileReader.XRawFileIO.ExtractParentIonMZFromFilterText(string, out double)</remarks>
        private bool ExtractParentIonMzFromFilterText(string filterText, out double parentIonMz)
        {
            Regex matcher;

            if (filterText.ToLower().Contains("msx"))
            {
                matcher = mParentIonMzMatchMsx;
            }
            else
            {
                matcher = mParentIonMzMatchNonMsx;
            }

            var reMatch = matcher.Match(filterText);

            if (reMatch.Success)
            {
                var parentIonMzText = reMatch.Groups["ParentMZ"].Value;
                var success = double.TryParse(parentIonMzText, out parentIonMz);
                return success;
            }

            parentIonMz = 0;
            return false;
        }

        private void MarkupPeptideWithMods()
        {
            // Markup the peptide with the dynamic and static mods

            var blnSuccess = ConvertModsToNumericMods(mPSMCurrent.Peptide.Trim(), out var strPeptideWithMods, out var lstPeptideMods);
            if (blnSuccess)
            {
                double dblTotalModMass = 0;

                mPSMCurrent.PeptideWithNumericMods = strPeptideWithMods;
                mPSMCurrent.ClearModifiedResidues();
                foreach (var objModEntry in lstPeptideMods)
                {
                    mPSMCurrent.AddModifiedResidue(objModEntry);
                    dblTotalModMass += objModEntry.ModDefinition.ModificationMass;
                }

                if (Math.Abs(mPSMCurrent.PeptideMonoisotopicMass) < double.Epsilon)
                {
                    mPSMCurrent.PeptideMonoisotopicMass = mPeptideMassCalculator.ComputeSequenceMass(mPSMCurrent.PeptideCleanSequence) + dblTotalModMass;
                }
            }
        }

        /// <summary>
        /// When FastReadMode is True, you first call MoveNext to read the peptide scores, then if the peptide
        /// is a peptide of interest, you call this function to finalize any processing steps that were skipped
        /// </summary>
        /// <remarks></remarks>
        public void FinalizeCurrentPSM()
        {
            if (mPSMCurrentFinalized)
                return;

            // Determine the clean sequence and cleavage state, and update the Seq_ID fields
            mPHRPParser.FinalizePSM(mPSMCurrent);

            MarkupPeptideWithMods();

            ComputePrecursorNeutralMass();

            mPSMCurrentFinalized = true;
        }

        private void ReadAndCacheMSGFData()
        {

            try
            {
                var strMSGFFilePath = GetMSGFFileName(mInputFilePath);
                strMSGFFilePath = Path.Combine(mInputFolderPath, strMSGFFilePath);

                ReadAndCacheMSGFData(strMSGFFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Exception determining MSGF file path", ex);
            }

        }

        private void ReadAndCacheMSGFData(string strMSGFFilePath)
        {

            try
            {
                strMSGFFilePath = AutoSwitchToLegacyMSGFDBIfRequired(strMSGFFilePath, mInputFilePath);
                strMSGFFilePath = AutoSwitchToFHTIfRequired(strMSGFFilePath, mInputFilePath);

                if (File.Exists(strMSGFFilePath))
                {
                    var oMSGFReader = new clsMSGFResultsReader();
                    mMSGFCachedResults = oMSGFReader.ReadMSGFData(strMSGFFilePath);

                    if (oMSGFReader.ErrorMessage.Length > 0)
                    {
                        ReportError("Error reading MSGF data: " + oMSGFReader.ErrorMessage);
                    }
                }
                else
                {
                    ReportWarning("MSGF file not found: " + strMSGFFilePath);
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading MSGF file", ex);
            }

        }

        /// <summary>
        /// Reads the data in strModSummaryFilePath.  Populates objDynamicMods and objStaticMods with the modification definitions
        /// </summary>
        /// <param name="strModSummaryFilePath">Path to the PHRP Mod Summary file to read</param>
        /// <param name="objDynamicMods">List with mod symbols as the key and the corresponding mod mass</param>
        /// <param name="objStaticMods">List with amino acid names as the key and the corresponding mod mass</param>
        /// <returns>True if success; false if an error</returns>
        private bool ReadModSummaryFile(
            string strModSummaryFilePath,
            IDictionary<char, clsModificationDefinition> objDynamicMods,
            IDictionary<string, List<clsModificationDefinition>> objStaticMods)
        {
            try
            {
                // Clear objDynamicMods and objStaticMods (should have been instantiated by the calling function)
                objDynamicMods.Clear();
                objStaticMods.Clear();

                if (string.IsNullOrEmpty(strModSummaryFilePath))
                {
                    ReportError("ModSummaryFile path is empty; unable to continue");
                    return false;
                }

                if (!File.Exists(strModSummaryFilePath))
                {
                    ReportError("ModSummary file not found: " + strModSummaryFilePath);
                    return false;
                }

                ShowMessage("Reading the PHRP ModSummary file");

                var objModSummaryReader = new clsPHRPModSummaryReader(strModSummaryFilePath);
                var blnSuccess = objModSummaryReader.Success;

                if (blnSuccess && objModSummaryReader.ModificationDefs.Count > 0)
                {
                    foreach (var objModDef in objModSummaryReader.ModificationDefs)
                    {
                        var strModMass = objModSummaryReader.GetModificationMassAsText(objModDef.MassCorrectionTag);

                        switch (objModDef.ModificationType)
                        {
                            case clsModificationDefinition.eModificationTypeConstants.StaticMod:
                            case clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod:
                            case clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod:

                                // "S", "T", or "P"
                                // Static residue mod, peptide terminus static mod, or protein terminus static mod
                                // Note that < and > mean peptide N and C terminus (clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS and clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                                // Note that [ and ] mean protein N and C terminus (clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS and clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)

                                // This mod could apply to multiple residues, so need to process each character in strTargetResidues
                                foreach (var chChar in objModDef.TargetResidues)
                                {
                                    try
                                    {
                                        if (objStaticMods.TryGetValue(chChar.ToString(), out var lstModDefs))
                                        {
                                            if (!lstModDefs.Contains(objModDef))
                                            {
                                                // Residue is already present in objStaticMods; this is unusual, but we'll allow it
                                                // We'll log a warning, but continue
                                                ShowMessage("Warning: Residue '" + chChar + "' has more than one static mod defined; this is not typically used, but will be allowed");
                                                lstModDefs.Add(objModDef);
                                            }
                                        }
                                        else
                                        {
                                            lstModDefs = new List<clsModificationDefinition> {
                                                objModDef
                                            };
                                            objStaticMods.Add(chChar.ToString(), lstModDefs);
                                        }
                                    }
                                    catch (Exception ex)
                                    {
                                        HandleException("Exception adding static mod for " + chChar + " with ModMass=" + strModMass, ex);
                                    }
                                }

                                break;
                            case clsModificationDefinition.eModificationTypeConstants.DynamicMod:
                                // Dynamic residue mod (Includes mod type "D")
                                // Note that < and > mean peptide N and C terminus (clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS and clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)

                                try
                                {
                                    if (objDynamicMods.ContainsKey(objModDef.ModificationSymbol))
                                    {
                                        // Mod symbol already present in objDynamicMods; this is unexpected
                                        // We'll log a warning, but continue
                                        ShowMessage("Warning: Dynamic mod symbol '" + objModDef.ModificationSymbol + "' is already defined; it cannot have more than one associated mod mass (duplicate has ModMass=" + strModMass + ")");
                                    }
                                    else
                                    {
                                        objDynamicMods.Add(objModDef.ModificationSymbol, objModDef);
                                    }
                                }
                                catch (Exception ex)
                                {
                                    HandleException("Exception adding dynamic mod for " + objModDef.ModificationSymbol + " with ModMass=" + strModMass, ex);
                                }

                                break;

                            case clsModificationDefinition.eModificationTypeConstants.IsotopicMod:
                                break;
                            // Isotopic mods are not supported by this class
                            // However, do not log a warning since these are rarely used

                            case clsModificationDefinition.eModificationTypeConstants.UnknownType:
                                // Unknown type; just ignore it
                                break;
                        }
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
                var strScanStatsFilePath = GetScanStatsFilename(mDatasetName);
                strScanStatsFilePath = Path.Combine(mInputFolderPath, strScanStatsFilePath);

                ReadScanStatsData(strScanStatsFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Exception determining Scan Stats file path", ex);
            }

        }

        private void ReadScanStatsData(string strScanStatsFilePath)
        {
            try
            {
                if (File.Exists(strScanStatsFilePath))
                {
                    var oScanStatsReader = new clsScanStatsReader();
                    mScanStats = oScanStatsReader.ReadScanStatsData(strScanStatsFilePath);

                    if (oScanStatsReader.ErrorMessage.Length > 0)
                    {
                        ReportError("Error reading ScanStats data: " + oScanStatsReader.ErrorMessage);
                    }
                }
                else
                {
                    ReportWarning("ScanStats file not found: " + strScanStatsFilePath);
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
                var strExtendedScanStatsFilePath = GetExtendedScanStatsFilename(mDatasetName);
                strExtendedScanStatsFilePath = Path.Combine(mInputFolderPath, strExtendedScanStatsFilePath);

                ReadExtendedScanStatsData(strExtendedScanStatsFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Exception determining Scan Stats file path", ex);
            }

        }

        private void ReadExtendedScanStatsData(string strExtendedScanStatsFilePath)
        {
            try
            {
                if (File.Exists(strExtendedScanStatsFilePath))
                {
                    var oExtendedScanStatsReader = new clsExtendedScanStatsReader();
                    mScanStatsEx = oExtendedScanStatsReader.ReadExtendedScanStatsData(strExtendedScanStatsFilePath);

                    if (oExtendedScanStatsReader.ErrorMessage.Length > 0)
                    {
                        ReportError("Error reading Extended ScanStats data: " + oExtendedScanStatsReader.ErrorMessage);
                    }
                }
                else
                {
                    // Note: we do not need to raise a warning for MSGFDB results since the extended scan stats file isn't needed
                    if (mPHRPParser.PeptideHitResultType != ePeptideHitResultType.MSGFDB)
                    {
                        ReportWarning("Extended ScanStats file not found: " + strExtendedScanStatsFilePath);
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
            mErrorMessage = message;
            if (mEchoMessagesToConsole)
                Console.WriteLine(message);
            mErrorMessages.Add(message);

            OnErrorEvent(message);
        }

        private void ReportWarning(string message)
        {
            if (mEchoMessagesToConsole)
                Console.WriteLine(message);
            mWarningMessages.Add(message);

            OnWarningEvent(message);
        }

        private void SetLocalErrorCode(ePHRPReaderErrorCodes eNewErrorCode, bool blnLeaveExistingErrorCodeUnchanged = false)
        {
            if (blnLeaveExistingErrorCodeUnchanged && mLocalErrorCode != ePHRPReaderErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mLocalErrorCode = eNewErrorCode;
            }
        }

        private void ShowMessage(string message)
        {
            if (mEchoMessagesToConsole)
                Console.WriteLine(message);

            OnStatusEvent(message);
        }

        private bool TryGetScanStats(int intScanNumber,
            out clsScanStatsInfo objScanStatsInfo)
        {
            if (mScanStats != null && mScanStats.Count > 0)
            {
                if (mScanStats.TryGetValue(intScanNumber, out objScanStatsInfo))
                {
                    return true;
                }
            }
            objScanStatsInfo = null;
            return false;
        }

        private bool TryGetExtendedScanStats(int intScanNumber,
            out clsScanStatsExInfo objExtendedScanStatsInfo)
        {
            if (mScanStatsEx != null && mScanStats.Count > 0)
            {
                if (mScanStatsEx.TryGetValue(intScanNumber, out objExtendedScanStatsInfo))
                {
                    return true;
                }
            }
            objExtendedScanStatsInfo = null;
            return false;
        }

        private bool ValidateInputFiles(string strInputFilePath, ref ePeptideHitResultType eResultType, ref string strModSummaryFilePath)
        {

            var inputFile = new FileInfo(strInputFilePath);
            if (!inputFile.Exists)
            {
                SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath);
                ReportError("Input file not found: " + strInputFilePath);
                return false;
            }

            // Try to auto-determine the result type if it is not specified
            if (eResultType == ePeptideHitResultType.Unknown)
            {
                eResultType = AutoDetermineResultType(strInputFilePath);
            }

            if (eResultType == ePeptideHitResultType.Unknown)
            {
                SetLocalErrorCode(ePHRPReaderErrorCodes.InputFileFormatNotRecognized);
                ReportError("Error: Unable to auto-determine file format for " + strInputFilePath);
                return false;
            }

            // Extract the dataset name from the input file path
            mDatasetName = AutoDetermineDatasetName(strInputFilePath, eResultType);
            if (string.IsNullOrEmpty(mDatasetName))
            {
                if (mStartupOptions.LoadModsAndSeqInfo || mStartupOptions.LoadMSGFResults || mStartupOptions.LoadScanStatsData)
                {
                    ReportError("Error: Unable to auto-determine the dataset name from the input file name: " + strInputFilePath);
                    SetLocalErrorCode(ePHRPReaderErrorCodes.InputFileFormatNotRecognized);
                    return false;
                }

                ReportWarning("Unable to auto-determine the dataset name from the input file name; this is not a critical error since not reading related files: " + strInputFilePath);
            }

            if (mStartupOptions.LoadModsAndSeqInfo)
            {
                strModSummaryFilePath = GetPHRPModSummaryFileName(eResultType, mDatasetName);
                strModSummaryFilePath = Path.Combine(inputFile.DirectoryName, strModSummaryFilePath);

                strModSummaryFilePath = AutoSwitchToLegacyMSGFDBIfRequired(strModSummaryFilePath, inputFile.Name);
                var strModSummaryFilePathPreferred = AutoSwitchToFHTIfRequired(strModSummaryFilePath, inputFile.Name);
                if (strModSummaryFilePath != strModSummaryFilePathPreferred && File.Exists(strModSummaryFilePathPreferred))
                {
                    strModSummaryFilePath = strModSummaryFilePathPreferred;
                }

                if (!ValidateRequiredFileExists("ModSummary file", strModSummaryFilePath) && inputFile.Name.ToLower().Contains("_fht"))
                {
                    SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound);
                    return false;
                }
            }
            else
            {
                strModSummaryFilePath = string.Empty;
            }

            return true;
        }

        private bool ValidateRequiredFileExists(string strFileDescription, string strFilePath, bool blnReportErrors = true)
        {
            if (string.IsNullOrEmpty(strFilePath))
            {
                if (blnReportErrors)
                {
                    ReportError(strFileDescription + " is not defined");
                }
                return false;
            }

            if (!File.Exists(strFilePath))
            {
                if (blnReportErrors)
                {
                    ReportError(strFileDescription + " not found: " + strFilePath);
                }
                return false;
            }

            return true;
        }

        #region "IDisposable Support"
        private bool disposedValue; // To detect redundant calls

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
