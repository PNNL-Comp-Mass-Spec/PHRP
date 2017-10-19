//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class is the base class for classes used to parse PHRP data lines
// It must be derived by a sub-class customized for the specific analysis tool (Sequest, X!Tandem, Inspect, etc.)
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser base class
    /// </summary>
    public abstract class clsPHRPParser : PRISM.clsEventNotifier
    {
        #region "Structures"

        /// <summary>
        /// Tracks ambigious modifications
        /// </summary>
        protected struct udtAmbiguousModInfo
        {
            /// <summary>
            /// First residue the mod could apply to
            /// </summary>
            public int ResidueStart;

            /// <summary>
            /// Last residue the mod could apply to
            /// </summary>
            public int ResidueEnd;

            /// <summary>
            /// Modification mass (as a string)
            /// </summary>
            public string ModMassString;
        }
        #endregion

        #region "Module variables"

        /// <summary>
        /// Dataset name
        /// </summary>
        protected string mDatasetName;

        /// <summary>
        /// Input file path
        /// </summary>
        private string mInputFilePath;

        /// <summary>
        /// Input folder path
        /// </summary>
        protected string mInputFolderPath;

        /// <summary>
        /// True if initialized
        /// </summary>
        private bool mInitialized;

        /// <summary>
        /// Column headers in the synopsis file and first hits file
        /// </summary>
        protected SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Error message
        /// </summary>
        protected string mErrorMessage = string.Empty;

        /// <summary>
        /// Cleavage state calculator
        /// </summary>
        protected readonly clsPeptideCleavageStateCalculator mCleavageStateCalculator;

        /// <summary>
        /// Peptide mass calculator
        /// </summary>
        protected readonly clsPeptideMassCalculator mPeptideMassCalculator;

        /// <summary>
        /// PHRP result type
        /// </summary>
        protected clsPHRPReader.ePeptideHitResultType mPeptideHitResultType;

        /// <summary>
        /// Modification info
        /// </summary>
        protected List<clsModificationDefinition> mModInfo;

        private SortedList<int, int> mResultToSeqMap;
        private SortedList<int, clsSeqInfo> mSeqInfo;
        private SortedList<int, List<clsProteinInfo>> mSeqToProteinMap;
        private Dictionary<string, clsPepToProteinMapInfo> mPepToProteinMap;

        /// <summary>
        /// Protein Names for each ResultID
        /// </summary>
        protected readonly SortedList<int, List<string>> mResultIDToProteins;

        private readonly List<string> mErrorMessages;
        private readonly List<string> mWarningMessages;

        #endregion

        #region "Properties"

        /// <summary>
        /// Cached error messages
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public List<string> ErrorMessages => mErrorMessages;

        /// <summary>
        /// Input file path
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string InputFilePath => mInputFilePath;

        /// <summary>
        /// Input folder path
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string InputFolderPath => mInputFolderPath;

        /// <summary>
        /// Maximum number of proteins to associate with each PSM
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>0 means to load all proteins</remarks>
        public int MaxProteinsPerPSM { get; set; }

        /// <summary>
        /// Peptide hit result type; Sequest, XTandem, Inspect, or MSGFDB
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPHRPReader.ePeptideHitResultType PeptideHitResultType => mPeptideHitResultType;

        /// <summary>
        /// Peptide to protein map file name
        /// </summary>
        /// <returns></returns>
        public Dictionary<string, clsPepToProteinMapInfo> PepToProteinMap => mPepToProteinMap;

        /// <summary>
        /// Returns the cached mapping between ResultID and SeqID
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedList<int, int> ResultToSeqMap => mResultToSeqMap;

        /// <summary>
        /// Returns the cached sequence info, where key is SeqID
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedList<int, clsSeqInfo> SeqInfo => mSeqInfo;

        /// <summary>
        /// Returns the cached sequence to protein map information
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedList<int, List<clsProteinInfo>> SeqToProteinMap => mSeqToProteinMap;

        /// <summary>
        /// Cached warning messages
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public List<string> WarningMessages => mWarningMessages;

        #endregion

        #region "Properties overridden by derived classes"

        /// <summary>
        /// First hits file
        /// </summary>
        public abstract string PHRPFirstHitsFileName { get; }

        /// <summary>
        /// Mod summary file
        /// </summary>
        public abstract string PHRPModSummaryFileName { get; }

        /// <summary>
        /// Peptide to protein map file
        /// </summary>
        public abstract string PHRPPepToProteinMapFileName { get; }

        /// <summary>
        /// Protein mods file
        /// </summary>
        public abstract string PHRPProteinModsFileName { get; }

        /// <summary>
        /// Synopsis file
        /// </summary>
        public abstract string PHRPSynopsisFileName { get; }

        /// <summary>
        /// Result to sequence map file
        /// </summary>
        public abstract string PHRPResultToSeqMapFileName { get; }

        /// <summary>
        /// Sequence info file
        /// </summary>
        public abstract string PHRPSeqInfoFileName { get; }

        /// <summary>
        /// Sequence to protein map file
        /// </summary>
        public abstract string PHRPSeqToProteinMapFileName { get; }

        /// <summary>
        /// Search engine name
        /// </summary>
        public abstract string SearchEngineName { get; }

        #endregion

        /// <summary>
        /// Initialize the parser for the given dataset, input file, and result type
        /// </summary>
        /// <param name="strDatasetName">Dataset Name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>If strInputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
        protected clsPHRPParser(string strDatasetName, string strInputFilePath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            bool blnLoadModsAndSeqInfo)
        {
            mErrorMessages = new List<string>();
            mWarningMessages = new List<string>();

            mResultIDToProteins = new SortedList<int, List<string>>();

            mCleavageStateCalculator = new clsPeptideCleavageStateCalculator();

            mPeptideMassCalculator = new clsPeptideMassCalculator();

            var startupOptions = new clsPHRPStartupOptions {LoadModsAndSeqInfo = blnLoadModsAndSeqInfo};

            InitializeParser(strDatasetName, strInputFilePath, ePeptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the parser for the given dataset, input file, and result type
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="ePeptideHitResultType"></param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
        /// <remarks>If strInputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
        protected clsPHRPParser(string strDatasetName, string strInputFilePath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType, clsPHRPStartupOptions startupOptions)
        {
            mErrorMessages = new List<string>();
            mWarningMessages = new List<string>();

            mResultIDToProteins = new SortedList<int, List<string>>();

            mCleavageStateCalculator = new clsPeptideCleavageStateCalculator();

            mPeptideMassCalculator = startupOptions.PeptideMassCalculator ?? new clsPeptideMassCalculator();

            InitializeParser(strDatasetName, strInputFilePath, ePeptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the parser for the given dataset and input file
        /// </summary>
        /// <param name="strDatasetName">Dataset Name</param>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="startupOptions">Startup options</param>
        /// <remarks>If strInputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable
        /// startupOptions.LoadModsAndSeqInfo controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read
        /// Setting startupOptions.MaxProteinsPerPSM to a non-zero value will limit the number of proteins that are tracked
        /// </remarks>
        private void InitializeParser(string strDatasetName, string strInputFilePath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType, clsPHRPStartupOptions startupOptions)
        {
            if (string.IsNullOrWhiteSpace(strDatasetName))
                strDatasetName = "Undefined";
            mDatasetName = strDatasetName;

            mPeptideHitResultType = ePeptideHitResultType;

            mMaxProteinsPerPSM = startupOptions.MaxProteinsPerPSM;

            var blnIsSynopsisFile = false;

            if (string.IsNullOrEmpty(strInputFilePath))
            {
                // User instantiated the class without a filename
                // Functions that solely require a dataset name will be callable, but cannot call functions that read a data line
                mInputFilePath = string.Empty;
                mInputFolderPath = string.Empty;

                startupOptions.LoadModsAndSeqInfo = false;
            }
            else
            {
                var inputFile = new FileInfo(strInputFilePath);
                mInputFilePath = inputFile.FullName;
                mInputFolderPath = inputFile.DirectoryName;

                var expectedSynopsisName = clsPHRPReader.GetPHRPSynopsisFileName(mPeptideHitResultType, mDatasetName);
                expectedSynopsisName = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(expectedSynopsisName, inputFile.Name);

                if (string.Equals(inputFile.Name, expectedSynopsisName, StringComparison.CurrentCultureIgnoreCase))
                {
                    blnIsSynopsisFile = true;
                }
            }

            mErrorMessage = string.Empty;

            // Initialize the column mapping object
            // Using a case-insensitive comparer
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase);

            // Initialize the tracking lists
            // These will get updated via the call to objReader.GetProteinMapping
            mResultToSeqMap = new SortedList<int, int>();
            mSeqInfo = new SortedList<int, clsSeqInfo>();
            mSeqToProteinMap = new SortedList<int, List<clsProteinInfo>>();
            mPepToProteinMap = new Dictionary<string, clsPepToProteinMapInfo>();

            if (startupOptions.LoadModsAndSeqInfo)
            {
                // Read the ModSummary file (if it exists)
                LoadModSummary();
            }

            if (startupOptions.LoadModsAndSeqInfo)
            {
                // Read the ResultToSeqMapInfo (if the files exist)
                if (blnIsSynopsisFile)
                {
                    // Assume the files exist
                    LoadSeqInfo();
                }
                else
                {
                    // Only continue if the fht versions exists

                    var strResultToSeqMapFilePath = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName);
                    var blnSeqInfoLoaded = false;

                    if (!string.IsNullOrEmpty(strResultToSeqMapFilePath))
                    {
                        strResultToSeqMapFilePath = Path.Combine(mInputFolderPath, strResultToSeqMapFilePath);
                        strResultToSeqMapFilePath = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strResultToSeqMapFilePath, mInputFilePath);
                        strResultToSeqMapFilePath = clsPHRPReader.AutoSwitchToFHTIfRequired(strResultToSeqMapFilePath, mInputFilePath);

                        if (File.Exists(strResultToSeqMapFilePath))
                        {
                            blnSeqInfoLoaded = LoadSeqInfo();
                        }
                    }

                    if (!blnSeqInfoLoaded)
                    {
                        if (string.IsNullOrEmpty(strResultToSeqMapFilePath))
                        {
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file and unable to determine the ResultToSeqMapFilename using clsPHRPReader.GetPHRPResultToSeqMapFileName()");
                        }
                        else
                        {
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file but the ResultToSeqMap file does not exist: " + strResultToSeqMapFilePath);
                        }
                    }
                }
            }

            // The following will be overridden by a derived form of this class
            DefineColumnHeaders();

            mInitialized = true;
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name; assumes blnLoadModsAndSeqInfo=True
        /// </summary>
        /// <param name="strInputFilePath">Input file path</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from strInputFilePath</remarks>
        public static clsPHRPParser GetParser(string strInputFilePath)
        {
            return GetParser(strInputFilePath, true);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name
        /// </summary>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from strInputFilePath</remarks>
        public static clsPHRPParser GetParser(string strInputFilePath, bool blnLoadModsAndSeqInfo)
        {
            var ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath);

            if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + strInputFilePath);
            }

            var strDatasetName = clsPHRPReader.AutoDetermineDatasetName(strInputFilePath);
            if (string.IsNullOrEmpty(strDatasetName))
            {
                throw new Exception("Unable to auto-determine the Dataset Name for " + strInputFilePath);
            }

            return GetParser(strInputFilePath, strDatasetName, ePeptideHitResultType, blnLoadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name
        /// </summary>
        /// <param name="strInputFilePath">Input file path</param>
        /// ''' <param name="strDatasetName">Dataset Name</param>
        /// <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type from strInputFilePath</remarks>
        public static clsPHRPParser GetParser(string strInputFilePath, string strDatasetName, bool blnLoadModsAndSeqInfo)
        {
            var ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath);

            if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + strInputFilePath);
            }

            return GetParser(strInputFilePath, strDatasetName, ePeptideHitResultType, blnLoadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on ePeptideHitResultType
        /// </summary>
        /// <param name="strInputFilePath">Input file path</param>
        /// <param name="strDatasetName">Dataset Name</param>
        /// <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks></remarks>
        public static clsPHRPParser GetParser(string strInputFilePath, string strDatasetName, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            bool blnLoadModsAndSeqInfo)
        {
            switch (ePeptideHitResultType)
            {
                case clsPHRPReader.ePeptideHitResultType.Inspect:
                    return new clsPHRPParserInspect(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MSAlign:
                    return new clsPHRPParserMSAlign(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MSGFDB:
                    return new clsPHRPParserMSGFDB(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.Sequest:
                    return new clsPHRPParserSequest(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.XTandem:
                    return new clsPHRPParserXTandem(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MODa:
                    return new clsPHRPParserMODa(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MODPlus:
                    return new clsPHRPParserMODPlus(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo);

                default:
                    throw new Exception("Unrecognized value for PeptideHitResultType: " + ePeptideHitResultType.ToString());
            }
        }

        #region "Functions overridden by derived classes"

        protected abstract void DefineColumnHeaders();

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="strLine">Data line</param>
        /// <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="objPSM">clsPSM object (output)</param>
        /// <returns>True if success, false if an error</returns>
        public bool ParsePHRPDataLine(string strLine, int intLinesRead, out clsPSM objPSM)
        {
            return ParsePHRPDataLine(strLine, intLinesRead, out objPSM, fastReadMode: false);
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="strLine">Data line</param>
        /// <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="objPSM">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
        public abstract bool ParsePHRPDataLine(string strLine, int intLinesRead, out clsPSM objPSM, bool fastReadMode);

        /// <summary>
        /// Parses the specified parameter file
        /// Also reads the Tool_Version_Info file in the same folder (if present)
        /// </summary>
        /// <param name="strSearchEngineParamFileName">Name of the parameter file to parse (must reside in InputFolderPath)</param>
        /// <param name="objSearchEngineParams">Search engine parameters class (output)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public abstract bool LoadSearchEngineParameters(string strSearchEngineParamFileName, out clsSearchEngineParameters objSearchEngineParams);

        #endregion

        protected void AddHeaderColumn(string strColumnName)
        {
            mColumnHeaders.Add(strColumnName, mColumnHeaders.Count);
        }

        protected void AddScore(clsPSM objPSM, string[] strColumns, string strScoreColumnName)
        {
            const string NOT_FOUND = "==SCORE_NOT_FOUND==";

            var strValue = clsPHRPReader.LookupColumnValue(strColumns, strScoreColumnName, mColumnHeaders, NOT_FOUND);

            if (strValue != NOT_FOUND)
            {
                objPSM.SetScore(strScoreColumnName, strValue);
            }
        }

        /// <summary>
        /// Clear any cached error messages
        /// </summary>
        /// <remarks></remarks>
        public void ClearErrors()
        {
            mErrorMessages.Clear();
        }

        /// <summary>
        /// Clear any cached warning messages
        /// </summary>
        /// <remarks></remarks>
        public void ClearWarnings()
        {
            mWarningMessages.Clear();
        }

        private readonly StringBuilder mNewPeptide = new StringBuilder();

        private string ConvertModsToNumericMods(string strCleanSequence, IReadOnlyCollection<clsAminoAcidModInfo> lstModifiedResidues)
        {
            mNewPeptide.Length = 0;

            if (lstModifiedResidues == null || lstModifiedResidues.Count == 0)
            {
                return strCleanSequence;
            }

            for (var intIndex = 0; intIndex <= strCleanSequence.Length - 1; intIndex++)
            {
                mNewPeptide.Append(strCleanSequence[intIndex]);

                foreach (var objModInfo in lstModifiedResidues)
                {
                    if (objModInfo.ResidueLocInPeptide == intIndex + 1)
                    {
                        mNewPeptide.Append(NumToStringPlusMinus(objModInfo.ModDefinition.ModificationMass, 4));
                    }
                }
            }

            return mNewPeptide.ToString();
        }

        /// <summary>
        /// Look for ambiguous mods in strSequenceWithMods
        /// For example, -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q
        /// </summary>
        /// <param name="strSequenceWithMods"></param>
        /// <returns></returns>
        /// <remarks>List of ambiguous mods, where the keys are the start residues and the values are the ambiguous mod info</remarks>
        private SortedList<int, udtAmbiguousModInfo> ExtractAmbiguousMods(string strSequenceWithMods)
        {
            if (!clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, out var strPrimarySequence, out _, out _))
            {
                strPrimarySequence = string.Copy(strSequenceWithMods);
            }

            var lstAmbiguousMods = new SortedList<int, udtAmbiguousModInfo>();

            var intResidueNumber = 0;
            var blnParsingAmbiguousMod = false;
            var udtCurrentMod = default(udtAmbiguousModInfo);

            for (var intCharIndex = 0; intCharIndex <= strPrimarySequence.Length - 1; intCharIndex++)
            {
                if (clsPHRPReader.IsLetterAtoZ(strPrimarySequence[intCharIndex]))
                {
                    // Found a letter
                    intResidueNumber += 1;

                    if (intCharIndex > 0 && strPrimarySequence[intCharIndex - 1] == '(')
                    {
                        // Found an ambiguous mod
                        if (!blnParsingAmbiguousMod)
                        {
                            blnParsingAmbiguousMod = true;
                            udtCurrentMod.ResidueStart = intResidueNumber;
                            udtCurrentMod.ResidueEnd = intResidueNumber;
                            udtCurrentMod.ModMassString = string.Empty;
                        }
                    }
                }
                else if (blnParsingAmbiguousMod)
                {
                    // Found a non-letter, and we are parsing an ambiguous mod

                    udtCurrentMod.ResidueEnd = intResidueNumber;
                    blnParsingAmbiguousMod = false;

                    // The mod mass should be next, in the form [-30.09]
                    // Parse out the mod mass
                    if (intCharIndex < strPrimarySequence.Length - 2)
                    {
                        if (strPrimarySequence[intCharIndex + 1] == '[')
                        {
                            var strModMassString = strPrimarySequence.Substring(intCharIndex + 2);
                            var bracketIndex = strModMassString.IndexOf(']');
                            if (bracketIndex > 0)
                            {
                                // Valid ambiguous mod found; store it
                                strModMassString = strModMassString.Substring(0, bracketIndex);
                                udtCurrentMod.ModMassString = string.Copy(strModMassString);

                                lstAmbiguousMods.Add(udtCurrentMod.ResidueStart, udtCurrentMod);
                            }
                        }
                    }
                }
            }

            return lstAmbiguousMods;
        }

        public void FinalizePSM(clsPSM objPSM)
        {
            objPSM.UpdateCleanSequence();

            objPSM.UpdateCleavageInfo(mCleavageStateCalculator);

            UpdatePSMUsingSeqInfo(objPSM);
        }

        private static KeyValuePair<string, string> GetMODaStaticModSetting(KeyValuePair<string, string> kvSetting, out string warningMessage)
        {
            var strKey = kvSetting.Key;
            var strValue = kvSetting.Value;

            if (!string.Equals(strKey, "ADD", StringComparison.CurrentCultureIgnoreCase))
            {
                throw new Exception("Key name is not ADD; this is not a MODa Static Mod Setting");
            }

            var commaIndex = strValue.IndexOf(',');

            if (commaIndex > 0)
            {
                var strResidue = strValue.Substring(0, commaIndex).Trim();
                strValue = strValue.Substring(commaIndex + 1).Trim();

                // Update Key to look like ADD_A or ADD_NTerm
                strKey = strKey + "_" + strResidue;

                kvSetting = new KeyValuePair<string, string>(strKey, strValue);
                warningMessage = string.Empty;
            }
            else
            {
                warningMessage = "Value for MODa keyword ADD does not contain a comma";
            }

            return kvSetting;
        }

        protected void HandleException(string strBaseMessage, Exception ex)
        {
            if (string.IsNullOrEmpty(strBaseMessage))
            {
                strBaseMessage = "Error";
            }

            ReportError(strBaseMessage + ": " + ex.Message);
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

        /// <summary>
        /// Reads the data in strModSummaryFilePath.  Populates mModInfo with the modification names, masses, and affected residues
        /// </summary>
        /// <returns>True if success; false if an error</returns>
        private void LoadModSummary()
        {
            try
            {
                var strModSummaryFilePath = clsPHRPReader.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName);
                if (string.IsNullOrEmpty(strModSummaryFilePath))
                {
                    ReportWarning("ModSummaryFile path is empty; unable to continue");
                    return;
                }

                strModSummaryFilePath = Path.Combine(mInputFolderPath, strModSummaryFilePath);

                strModSummaryFilePath = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strModSummaryFilePath, mInputFilePath);
                var strModSummaryFilePathPreferred = clsPHRPReader.AutoSwitchToFHTIfRequired(strModSummaryFilePath, mInputFilePath);
                if (strModSummaryFilePath != strModSummaryFilePathPreferred && File.Exists(strModSummaryFilePathPreferred))
                {
                    strModSummaryFilePath = strModSummaryFilePathPreferred;
                }

                if (!File.Exists(strModSummaryFilePath))
                {
                    ReportWarning("ModSummary file not found: " + strModSummaryFilePath);
                    return;
                }

                ShowMessage("Reading the PHRP ModSummary file");

                var objModSummaryReader = new clsPHRPModSummaryReader(strModSummaryFilePath);
                var blnSuccess = objModSummaryReader.Success;

                if (blnSuccess)
                {
                    mModInfo = objModSummaryReader.ModificationDefs;
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading PHRP Mod Summary file", ex);
            }

        }

        private bool LoadSeqInfo()
        {
            bool blnSuccess;

            var entriesParsed = 0;
            var dtLastProgress = DateTime.UtcNow;
            var blnNotifyComplete = false;

            try
            {
                ShowMessage("Reading the PHRP SeqInfo file");

                // Instantiate the reader
                var objReader =
                    new clsPHRPSeqMapReader(mDatasetName, mInputFolderPath, mPeptideHitResultType, mInputFilePath)
                    {
                        MaxProteinsPerSeqID = mMaxProteinsPerPSM
                    };


                // Read the files
                blnSuccess = objReader.GetProteinMapping(out mResultToSeqMap, out mSeqToProteinMap, out mSeqInfo, out mPepToProteinMap);

                if (!blnSuccess)
                {
                    ReportWarning(objReader.ErrorMessage);
                }

                mResultIDToProteins.Clear();

                if (blnSuccess)
                {
                    // Populate mResultIDToProteins

                    foreach (var objItem in mResultToSeqMap)
                    {
                        List<string> lstProteinsForResultID;

                        if (mSeqToProteinMap.TryGetValue(objItem.Value, out var lstProteinsForSeqID))
                        {
                            lstProteinsForResultID = (from objProtein in lstProteinsForSeqID select objProtein.ProteinName).ToList();
                            // lstProteinsForResultID = new List<string>(lstProteinsForSeqID.Count);
                            // foreach (clsProteinInfo objProtein in lstProteinsForSeqID)
                            // {
                            //     if (!lstProteinsForResultID.Contains(objProtein.ProteinName))
                            //     {
                            //         lstProteinsForResultID.Add(objProtein.ProteinName);
                            //     }
                            // }
                        }
                        else
                        {
                            lstProteinsForResultID = new List<string>();
                        }

                        if (mMaxProteinsPerPSM > 0 && lstProteinsForResultID.Count > mMaxProteinsPerPSM)
                        {
                            // Only add a subset of the proteins in lstProteinsForResultID
                            var lstProteinSubset = lstProteinsForResultID.Take(mMaxProteinsPerPSM).OrderBy(item => item).ToList();
                            mResultIDToProteins.Add(objItem.Key, lstProteinSubset);
                        }
                        else
                        {
                            mResultIDToProteins.Add(objItem.Key, lstProteinsForResultID);
                        }

                        entriesParsed += 1;
                        if (DateTime.UtcNow.Subtract(dtLastProgress).TotalSeconds >= 5)
                        {
                            var pctComplete = entriesParsed / Convert.ToDouble(mResultToSeqMap.Count) * 100;
                            Console.WriteLine(" ... associating proteins with sequences: " + pctComplete.ToString("0.0") + "% complete");
                            dtLastProgress = DateTime.UtcNow;
                            blnNotifyComplete = true;
                        }
                    }

                    if (blnNotifyComplete)
                    {
                        Console.WriteLine(" ... associating proteins with sequences: 100% complete");
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error loading PHRP Seq Info", ex);
                blnSuccess = false;
                if (!mInitialized)
                    throw new Exception(mErrorMessage, ex);
            }

            return blnSuccess;
        }

        /// <summary>
        /// Formats a number so that it begins with a + sign if positive or a - sign if negative
        /// Rounds the number to the specified number of digits, trimming off trailing zeros
        /// Example output: +79.9663 or -17.016
        /// </summary>
        /// <param name="value"></param>
        /// <param name="digitsOfPrecision"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static string NumToStringPlusMinus(double value, int digitsOfPrecision)
        {
            var strFormatString = "+0;-0";
            if (digitsOfPrecision > 0)
            {
                strFormatString = "+0." + new string('0', digitsOfPrecision) + ";-0." + new string('0', digitsOfPrecision);
            }

            var strValue = value.ToString(strFormatString).TrimEnd('0');

            if (strValue.EndsWith("."))
            {
                // Ends in a decimal point; remove the decimal point
                strValue = strValue.TrimEnd('.');
            }

            return strValue;
        }

        /// <summary>
        /// Parse the column names in strSplitLine and update the local column header mapping
        /// </summary>
        /// <param name="strSplitLine"></param>
        /// <remarks></remarks>
        public void ParseColumnHeaders(string[] strSplitLine)
        {
            clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders);
        }

        /// <summary>
        /// Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
        /// </summary>
        /// <param name="strText"></param>
        /// <param name="chDelimiter"></param>
        /// <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
        /// <remarks>Automatically trims whitespace</remarks>
        public static KeyValuePair<string, string> ParseKeyValueSetting(string strText, char chDelimiter)
        {
            return ParseKeyValueSetting(strText, chDelimiter, string.Empty);
        }

        /// <summary>
        /// Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
        /// </summary>
        /// <param name="strText"></param>
        /// <param name="chDelimiter"></param>
        /// <param name="strCommentChar">If defined, then looks for this character in the value portion of the setting and removes that character plus any text after it</param>
        /// <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
        /// <remarks>Automatically trims whitespace</remarks>
        public static KeyValuePair<string, string> ParseKeyValueSetting(string strText, char chDelimiter, string strCommentChar)
        {
            if (string.IsNullOrEmpty(strText))
                return new KeyValuePair<string, string>(string.Empty, string.Empty);

            var intCharIndex = strText.IndexOf(chDelimiter);
            if (intCharIndex <= 0)
                return new KeyValuePair<string, string>(string.Empty, string.Empty);

            var strKey = strText.Substring(0, intCharIndex).Trim();
            string strValue;
            if (intCharIndex < strText.Length - 1)
            {
                strValue = strText.Substring(intCharIndex + 1).Trim();

                if (!string.IsNullOrEmpty(strCommentChar))
                {
                    // Look for the comment character
                    var commentCharIndex = strValue.IndexOf(strCommentChar, StringComparison.Ordinal);
                    if (commentCharIndex > 0)
                    {
                        // Trim off the comment
                        strValue = strValue.Substring(0, commentCharIndex).Trim();
                    }
                }
            }
            else
            {
                strValue = string.Empty;
            }

            var kvSetting = new KeyValuePair<string, string>(strKey, strValue);
            return kvSetting;
        }

        /// <summary>
        /// Read a Search Engine parameter file where settings are stored as key/value pairs
        /// </summary>
        /// <param name="strSearchEngineName">Search engine name (e.g. MSGF+)</param>
        /// <param name="strSearchEngineParamFileName">Search engine parameter file name (must exist in mInputFolderPath)</param>
        /// <param name="ePeptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
        /// <param name="objSearchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
        /// <returns>True if success, false if an error</returns>
        protected bool ReadKeyValuePairSearchEngineParamFile(
            string strSearchEngineName,
            string strSearchEngineParamFileName,
            clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            clsSearchEngineParameters objSearchEngineParams)
        {
            var paramFilePath = Path.Combine(mInputFolderPath, strSearchEngineParamFileName);

            var success = ReadKeyValuePairSearchEngineParamFile(strSearchEngineName, paramFilePath, ePeptideHitResultType, objSearchEngineParams,
                out var errorMessage, out var warningMessage);

            if (!string.IsNullOrWhiteSpace(errorMessage))
            {
                ReportError(errorMessage);
            }

            if (!string.IsNullOrWhiteSpace(warningMessage))
            {
                ReportWarning(warningMessage);
            }

            return success;
        }

        /// <summary>
        /// Read a Search Engine parameter file where settings are stored as key/value pairs
        /// </summary>
        /// <param name="searchEngineName">Search engine name (e.g. MS-GF+)</param>
        /// <param name="paramFilePath">Search engine parameter file path</param>
        /// <param name="ePeptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
        /// <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
        /// <param name="errorMessage">Output: error message</param>
        /// <param name="warningMessage">Output: warning message</param>
        /// <returns>True if success, false if an error</returns>
        public static bool ReadKeyValuePairSearchEngineParamFile(
            string searchEngineName,
            string paramFilePath,
            clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            clsSearchEngineParameters searchEngineParams,
            out string errorMessage,
            out string warningMessage)
        {
            errorMessage = string.Empty;
            warningMessage = string.Empty;

            try
            {
                if (string.IsNullOrWhiteSpace(searchEngineName))
                    searchEngineName = "?? Unknown tool ??";

                if (!File.Exists(paramFilePath))
                {
                    errorMessage = searchEngineName + " param file not found: " + paramFilePath;
                    return false;
                }

                searchEngineParams.UpdateSearchEngineParamFilePath(paramFilePath);

                using (var srInFile = new StreamReader(new FileStream(paramFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        var lineIn = srInFile.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        var dataLine = lineIn.TrimStart();

                        if (string.IsNullOrWhiteSpace(dataLine) || dataLine.StartsWith("#") || !dataLine.Contains('='))
                        {
                            continue;
                        }

                        // Split the line on the equals sign
                        var kvSetting = ParseKeyValueSetting(dataLine, '=', "#");

                        if (string.IsNullOrEmpty(kvSetting.Key))
                        {
                            continue;
                        }

                        if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.MODa &&
                            string.Equals(kvSetting.Key, "add", StringComparison.CurrentCultureIgnoreCase))
                        {
                            // ModA defines all of its static modifications with the ADD keyword
                            // Split the value at the comma and create a new setting entry with the residue name
                            kvSetting = GetMODaStaticModSetting(kvSetting, out var warningMessageAddon);
                            if (!string.IsNullOrWhiteSpace(warningMessageAddon))
                            {
                                if (string.IsNullOrWhiteSpace(warningMessage))
                                {
                                    warningMessage = string.Copy(warningMessageAddon);
                                }
                                else if (!warningMessage.Contains(warningMessageAddon))
                                {
                                    warningMessage += "; " + string.Copy(warningMessageAddon);
                                }
                            }
                        }

                        searchEngineParams.AddUpdateParameter(kvSetting);
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                errorMessage = string.Format("Error in ReadKeyValuePairSearchEngineParamFile for {0}, param file {1}: {2}", searchEngineName,
                    Path.GetFileName(paramFilePath), ex.Message);
                return false;
            }
        }

        protected bool ReadSearchEngineVersion(
            clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            clsSearchEngineParameters objSearchEngineParams)
        {
            var blnSuccess = false;

            try
            {
                // Read the Tool_Version_Info file to determine the analysis time and the tool version
                var toolVersionInfoFilePath = Path.Combine(mInputFolderPath, clsPHRPReader.GetToolVersionInfoFilename(ePeptideHitResultType));

                if (!File.Exists(toolVersionInfoFilePath) && ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.MSGFDB)
                {
                    // This could be an older MSGF+ job; check for a _MSGFDB.txt tool version file
                    var alternativeVersionInfoFilePath = Path.Combine(mInputFolderPath, "Tool_Version_Info_MSGFDB.txt");
                    if (File.Exists(alternativeVersionInfoFilePath))
                    {
                        toolVersionInfoFilePath = alternativeVersionInfoFilePath;
                    }
                }

                if (!File.Exists(toolVersionInfoFilePath))
                {
                    ReportWarning("Tool version info file not found: " + toolVersionInfoFilePath);
                    return false;
                }

                var strSearchEngineVersion = "Unknown";
                var dtSearchDate = new DateTime(1980, 1, 1);
                var blnValidDate = false;
                var blnValidVersion = false;

                using (var srInFile = new StreamReader(new FileStream(toolVersionInfoFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        var lineIn = srInFile.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        var dataLine = lineIn.TrimStart();

                        // Split the line on a colon
                        var kvSetting = ParseKeyValueSetting(dataLine, ':');

                        switch (kvSetting.Key.ToLower())
                        {
                            case "date":
                                blnValidDate = DateTime.TryParse(kvSetting.Value, out dtSearchDate);

                                break;
                            case "toolversioninfo":
                                if (!string.IsNullOrEmpty(kvSetting.Value))
                                {
                                    strSearchEngineVersion = string.Copy(kvSetting.Value);
                                    blnValidVersion = true;
                                }
                                else
                                {
                                    // The next line contains the search engine version
                                    if (!srInFile.EndOfStream)
                                    {
                                        var searchEngineLine = srInFile.ReadLine();
                                        if (!string.IsNullOrEmpty(searchEngineLine))
                                        {
                                            strSearchEngineVersion = string.Copy(searchEngineLine.Trim());
                                            blnValidVersion = true;
                                        }
                                    }
                                }
                                break;
                            default:
                                // Ignore the line
                                break;
                        }
                    }
                }

                if (!blnValidDate)
                {
                    ReportError("Date line not found in the ToolVersionInfo file");
                }
                else if (!blnValidVersion)
                {
                    ReportError("ToolVersionInfo line not found in the ToolVersionInfo file");
                }
                else
                {
                    blnSuccess = true;
                }

                objSearchEngineParams.UpdateSearchEngineVersion(strSearchEngineVersion);
                objSearchEngineParams.UpdateSearchDate(dtSearchDate);

                return blnSuccess;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineVersion: " + ex.Message);
                return false;
            }

        }

        protected void ReportError(string message)
        {
            mErrorMessage = message;
            mErrorMessages.Add(message);
            OnErrorEvent(message);
        }

        protected void ReportWarning(string message)
        {
            mWarningMessages.Add(message);
            OnWarningEvent(message);
        }

        protected void ShowMessage(string message)
        {
            OnStatusEvent(message);
        }

        private void StoreModInfo(clsPSM objPSM, clsSeqInfo objSeqInfo)
        {

            var lstNTerminalModsAdded = new List<string>();
            var lstCTerminalModsAdded = new List<string>();
            var intPeptideResidueCount = objPSM.PeptideCleanSequence.Length;

            objPSM.PeptideMonoisotopicMass = objSeqInfo.MonoisotopicMass;

            objPSM.ClearModifiedResidues();

            if (objSeqInfo.ModCount <= 0)
                return;

            // Split objSeqInfo.ModDescription on the comma character
            var strMods = objSeqInfo.ModDescription.Split(',');

            if (strMods.Length <= 0)
                return;

            // Parse objPSM.Peptide to look for ambiguous mods, for example -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q

            var lstAmbiguousMods = ExtractAmbiguousMods(objPSM.Peptide);

            // ReSharper disable once UseImplicitlyTypedVariableEvident

            for (var intModIndex = 0; intModIndex <= strMods.Length - 1; intModIndex++)
            {
                // Split strMods on the colon characters
                var kvModDetails = ParseKeyValueSetting(strMods[intModIndex], ':');

                if (string.IsNullOrEmpty(kvModDetails.Key) || string.IsNullOrEmpty(kvModDetails.Value))
                    continue;

                var strMassCorrectionTag = kvModDetails.Key;
                if (!int.TryParse(kvModDetails.Value, out var intResidueLoc))
                    continue;

                // Find the modification definition in mModInfo
                // Note that a given mass correction tag might be present multiple times in mModInfo, since it could be used as both a static peptide mod and a static peptide terminal mod
                // Thus, if intResidueLoc = 1 or intResidueLoc = objPSM.PeptideCleanSequence.Length then we'll first look for a peptide or protein terminal static mod

                bool blnFavorTerminalMods;
                clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState;

                if (intResidueLoc == 1)
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                    if (lstNTerminalModsAdded.Contains(strMassCorrectionTag))
                    {
                        // We have likely already added this modification as an N-terminal mod, thus, don't favor terminal mods this time
                        // An example is an iTraq peptide where there is a K at the N-terminus
                        // It gets modified with iTraq twice: once because of the N-terminus and once because of Lysine
                        // For example, R.K+144.102063+144.102063TGSY+79.9663GALAEITASK+144.102063.E
                        blnFavorTerminalMods = false;
                    }
                    else
                    {
                        blnFavorTerminalMods = true;
                    }
                }
                else if (intResidueLoc == intPeptideResidueCount)
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                    if (lstCTerminalModsAdded.Contains(strMassCorrectionTag))
                    {
                        blnFavorTerminalMods = false;
                    }
                    else
                    {
                        blnFavorTerminalMods = true;
                    }
                }
                else
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                    blnFavorTerminalMods = false;
                }

                clsModificationDefinition objMatchedModDef;
                bool blnMatchFound;

                if (mModInfo == null)
                {
                    objMatchedModDef = new clsModificationDefinition {MassCorrectionTag = strMassCorrectionTag};
                    blnMatchFound = true;
                }
                else
                {
                    blnMatchFound = UpdatePSMFindMatchingModInfo(strMassCorrectionTag, blnFavorTerminalMods, eResidueTerminusState,
                                                                 out objMatchedModDef);
                }

                if (blnMatchFound)
                {
                    var lstMatches = (from item in lstAmbiguousMods where item.Key == intResidueLoc select item.Value).ToList();

                    if (lstMatches.Count > 0)
                    {
                        // Ambiguous modification
                        objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence[intResidueLoc - 1], intResidueLoc,
                                                  eResidueTerminusState, objMatchedModDef, lstMatches.First().ResidueEnd);
                    }
                    else
                    {
                        // Normal, non-ambiguous modified residue
                        objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence[intResidueLoc - 1], intResidueLoc,
                                                  eResidueTerminusState, objMatchedModDef);
                    }

                    if (intResidueLoc == 1)
                    {
                        lstNTerminalModsAdded.Add(objMatchedModDef.MassCorrectionTag);
                    }
                    else if (intResidueLoc == intPeptideResidueCount)
                    {
                        lstCTerminalModsAdded.Add(objMatchedModDef.MassCorrectionTag);
                    }
                }
                else
                {
                    // Could not find a valid entry in mModInfo
                    ReportError("Unrecognized mass correction tag found in the SeqInfo file: " + strMassCorrectionTag);
                }
            }
        }

        /// <summary>
        /// Updates the theoretical (computed) monoisotopic mass of objPSM using mResultToSeqMap and mSeqInfo
        /// Also updates the modification info
        /// Also updates SeqID
        /// </summary>
        /// <param name="objPSM"></param>
        /// <returns>True if success, False if objPSM.ResultID is not found in mResultToSeqMap</returns>
        /// <remarks></remarks>
        protected bool UpdatePSMUsingSeqInfo(clsPSM objPSM)
        {
            var blnSuccess = false;

            // First determine the modified residues present in this peptide
            if (mResultToSeqMap != null && mResultToSeqMap.Count > 0)
            {
                if (mResultToSeqMap.TryGetValue(objPSM.ResultID, out var intSeqID))
                {
                    objPSM.SeqID = intSeqID;

                    if (mSeqInfo.TryGetValue(intSeqID, out var objSeqInfo))
                    {
                        StoreModInfo(objPSM, objSeqInfo);
                        blnSuccess = true;
                    }

                    // Lookup the protein details using mSeqToProteinMap
                    if (mSeqToProteinMap.TryGetValue(intSeqID, out var lstProteinDetails))
                    {
                        foreach (var oProtein in lstProteinDetails)
                        {
                            objPSM.AddProteinDetail(oProtein);
                        }
                    }

                    // Make sure all of the proteins in objPSM.Proteins are defined in objPSM.ProteinDetails
                    var addnlProteins1 = objPSM.Proteins.Except(objPSM.ProteinDetails.Keys, StringComparer.CurrentCultureIgnoreCase).ToList();

                    foreach (var proteinName in addnlProteins1)
                    {
                        if (mMaxProteinsPerPSM > 0 && objPSM.ProteinDetails.Count > mMaxProteinsPerPSM)
                        {
                            // Maximum number of proteins reached (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)
                            break;
                        }

                        var oProtein = new clsProteinInfo(proteinName, 0, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None);
                        objPSM.ProteinDetails.Add(proteinName, oProtein);
                    }

                    // Make sure all of the proteins in objPSM.ProteinDetails are defined in objPSM.Proteins
                    var addnlProteins2 = (from item in objPSM.ProteinDetails select item.Key).Except(objPSM.Proteins, StringComparer.CurrentCultureIgnoreCase).ToList();

                    var additionThresholdCheck = mMaxProteinsPerPSM;
                    if (additionThresholdCheck < int.MaxValue)
                    {
                        additionThresholdCheck += 1;
                    }

                    if (mMaxProteinsPerPSM > 0 && objPSM.Proteins.Count + addnlProteins2.Count > additionThresholdCheck)
                    {
                        // Maximum number of proteins will be reached; only add a subset of the proteins in addnlProteins2
                        // (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)

                        foreach (var oProtein in addnlProteins2)
                        {
                            if (mMaxProteinsPerPSM > 0 && objPSM.Proteins.Count >= mMaxProteinsPerPSM)
                            {
                                // Maximum number of proteins reached
                                break;
                            }
                            objPSM.Proteins.Add(oProtein);
                        }
                    }
                    else
                    {
                        objPSM.Proteins.AddRange(addnlProteins2);
                    }

                    if (mPepToProteinMap.Count > 0)
                    {
                        // Make sure the residue start/end locations are up-to-date in objPSM.ProteinDetails

                        if (mPepToProteinMap.TryGetValue(objPSM.PeptideCleanSequence, out var oPepToProteinMapInfo))
                        {
                            foreach (var oProtein in objPSM.ProteinDetails)
                            {
                                // Find the matching protein in oPepToProteinMapInfo

                                if (oPepToProteinMapInfo.ProteinMapInfo.TryGetValue(oProtein.Key, out var lstLocations))
                                {
                                    var udtFirstLocation = lstLocations.First();
                                    oProtein.Value.UpdateLocationInProtein(udtFirstLocation.ResidueStart, udtFirstLocation.ResidueEnd);
                                }
                            }
                        }
                    }
                }
            }

            if (blnSuccess)
            {
                if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(objPSM.Peptide, out _, out var strPrefix, out var strSuffix))
                {
                    objPSM.PeptideWithNumericMods = strPrefix + "." + ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues) + "." + strSuffix;
                }
                else
                {
                    objPSM.PeptideWithNumericMods = ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues);
                }
            }

            return blnSuccess;
        }

        private bool UpdatePSMFindMatchingModInfo(
            string strMassCorrectionTag,
            bool blnFavorTerminalMods,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            out clsModificationDefinition objMatchedModDef)
        {
            objMatchedModDef = new clsModificationDefinition();

            if (mModInfo == null)
            {
                return false;
            }


            var lstMatchedDefs = new List<clsModificationDefinition>();

            foreach (var objMod in mModInfo)
            {
                if (string.Equals(strMassCorrectionTag, objMod.MassCorrectionTag, StringComparison.InvariantCultureIgnoreCase))
                {
                    lstMatchedDefs.Add(objMod);
                }
            }

            var blnMatchFound = false;

            if (lstMatchedDefs.Count > 0)
            {
                while (!blnMatchFound)
                {
                    if (blnFavorTerminalMods)
                    {
                        // Look for an entry in lstMatchedDefs that is a terminal mod
                        foreach (var objMod in lstMatchedDefs)
                        {
                            if (objMod.ModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod ||
                                objMod.ModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod)
                            {
                                if (eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus &&
                                    (objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) || objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)))
                                {
                                    blnMatchFound = true;
                                    objMatchedModDef = objMod;
                                    break;
                                }

                                if (eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus &&
                                    (objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS) || objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)))
                                {
                                    blnMatchFound = true;
                                    objMatchedModDef = objMod;
                                    break;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Look for an entry in lstMatchedDefs that is not a terminal mod
                        foreach (var objMod in lstMatchedDefs)
                        {
                            if (!(objMod.ModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod ||
                                  objMod.ModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod))
                            {
                                blnMatchFound = true;
                                objMatchedModDef = objMod;
                                break;
                            }
                        }
                    }

                    if (!blnMatchFound)
                    {
                        if (blnFavorTerminalMods)
                        {
                            blnFavorTerminalMods = false;
                        }
                        else
                        {
                            // Still no match found (this shouldn't happen); use the first entry in lstMatchedDefs
                            objMatchedModDef = lstMatchedDefs[0];
                            blnMatchFound = true;
                        }
                    }
                }
            }

            return blnMatchFound;
        }
    }
}