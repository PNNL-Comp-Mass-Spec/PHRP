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
    public abstract class clsPHRPParser : PRISM.EventNotifier
    {
        #region "Structures"

        /// <summary>
        /// Tracks ambiguous modifications
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
        /// Input directory path
        /// </summary>
        protected string mInputDirectoryPath;

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
        /// Input directory path
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string InputDirectoryPath => mInputDirectoryPath;

        /// <summary>
        /// Input directory path
        /// </summary>
        [Obsolete("Use InputDirectoryPath")]
        public string InputFolderPath => mInputDirectoryPath;

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
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>If inputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
        protected clsPHRPParser(string datasetName, string inputFilePath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            bool loadModsAndSeqInfo)
        {
            mErrorMessages = new List<string>();
            mWarningMessages = new List<string>();

            mResultIDToProteins = new SortedList<int, List<string>>();

            mCleavageStateCalculator = new clsPeptideCleavageStateCalculator();

            mPeptideMassCalculator = new clsPeptideMassCalculator();

            var startupOptions = new clsPHRPStartupOptions { LoadModsAndSeqInfo = loadModsAndSeqInfo };

            InitializeParser(datasetName, inputFilePath, ePeptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the parser for the given dataset, input file, and result type
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="ePeptideHitResultType"></param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks>If inputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
        protected clsPHRPParser(string datasetName, string inputFilePath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType, clsPHRPStartupOptions startupOptions)
        {
            mErrorMessages = new List<string>();
            mWarningMessages = new List<string>();

            mResultIDToProteins = new SortedList<int, List<string>>();

            mCleavageStateCalculator = new clsPeptideCleavageStateCalculator();

            mPeptideMassCalculator = startupOptions.PeptideMassCalculator ?? new clsPeptideMassCalculator();

            InitializeParser(datasetName, inputFilePath, ePeptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the parser for the given dataset and input file
        /// </summary>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="startupOptions">Startup options</param>
        /// <remarks>If inputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable
        /// startupOptions.LoadModsAndSeqInfo controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read
        /// Setting startupOptions.MaxProteinsPerPSM to a non-zero value will limit the number of proteins that are tracked
        /// </remarks>
        private void InitializeParser(string datasetName, string inputFilePath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType, clsPHRPStartupOptions startupOptions)
        {
            if (string.IsNullOrWhiteSpace(datasetName))
                datasetName = "Undefined";

            mDatasetName = datasetName;

            mPeptideHitResultType = ePeptideHitResultType;

            MaxProteinsPerPSM = startupOptions.MaxProteinsPerPSM;

            var isSynopsisFile = false;

            if (string.IsNullOrEmpty(inputFilePath))
            {
                // User instantiated the class without a filename
                // Functions that solely require a dataset name will be callable, but cannot call functions that read a data line
                mInputFilePath = string.Empty;
                mInputDirectoryPath = string.Empty;

                startupOptions.LoadModsAndSeqInfo = false;
            }
            else
            {
                var inputFile = new FileInfo(inputFilePath);
                mInputFilePath = inputFile.FullName;
                if (inputFile.Directory != null)
                {
                    mInputDirectoryPath = inputFile.Directory.FullName;
                }

                var expectedSynopsisName = clsPHRPReader.GetPHRPSynopsisFileName(mPeptideHitResultType, mDatasetName);
                expectedSynopsisName = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(expectedSynopsisName, inputFile.Name);

                if (string.Equals(inputFile.Name, expectedSynopsisName, StringComparison.CurrentCultureIgnoreCase))
                {
                    isSynopsisFile = true;
                }
            }

            mErrorMessage = string.Empty;

            // Initialize the column mapping object
            // Using a case-insensitive comparer
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase);

            // Initialize the tracking lists
            // These will get updated via the call to reader.GetProteinMapping
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
                if (isSynopsisFile)
                {
                    // Assume the files exist
                    LoadSeqInfo();
                }
                else
                {
                    // Only continue if the fht versions exists

                    var resultToSeqMapFilePath = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName);
                    var seqInfoLoaded = false;

                    if (!string.IsNullOrEmpty(resultToSeqMapFilePath))
                    {
                        if (!string.IsNullOrWhiteSpace(mInputDirectoryPath))
                            resultToSeqMapFilePath = Path.Combine(mInputDirectoryPath, resultToSeqMapFilePath);

                        resultToSeqMapFilePath = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(resultToSeqMapFilePath, mInputFilePath);
                        resultToSeqMapFilePath = clsPHRPReader.AutoSwitchToFHTIfRequired(resultToSeqMapFilePath, mInputFilePath);

                        if (File.Exists(resultToSeqMapFilePath))
                        {
                            seqInfoLoaded = LoadSeqInfo();
                        }
                    }

                    if (!seqInfoLoaded)
                    {
                        if (string.IsNullOrEmpty(resultToSeqMapFilePath))
                        {
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file and unable to determine the ResultToSeqMapFilename using clsPHRPReader.GetPHRPResultToSeqMapFileName()");
                        }
                        else
                        {
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file but the ResultToSeqMap file does not exist: " + resultToSeqMapFilePath);
                        }
                    }
                }
            }

            // The following will be overridden by a derived form of this class
            DefineColumnHeaders();

            mInitialized = true;
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from inputFilePath</remarks>
        public static clsPHRPParser GetParser(string inputFilePath)
        {
            return GetParser(inputFilePath, true);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from inputFilePath</remarks>
        public static clsPHRPParser GetParser(string inputFilePath, bool loadModsAndSeqInfo)
        {
            var ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(inputFilePath);

            if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + inputFilePath);
            }

            var datasetName = clsPHRPReader.AutoDetermineDatasetName(inputFilePath);
            if (string.IsNullOrEmpty(datasetName))
            {
                throw new Exception("Unable to auto-determine the Dataset Name for " + inputFilePath);
            }

            return GetParser(inputFilePath, datasetName, ePeptideHitResultType, loadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// ''' <param name="datasetName">Dataset Name</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type from inputFilePath</remarks>
        public static clsPHRPParser GetParser(string inputFilePath, string datasetName, bool loadModsAndSeqInfo)
        {
            var ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(inputFilePath);

            if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + inputFilePath);
            }

            return GetParser(inputFilePath, datasetName, ePeptideHitResultType, loadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on ePeptideHitResultType
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks></remarks>
        public static clsPHRPParser GetParser(string inputFilePath, string datasetName, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            bool loadModsAndSeqInfo)
        {
            switch (ePeptideHitResultType)
            {
                case clsPHRPReader.ePeptideHitResultType.Inspect:
                    return new clsPHRPParserInspect(datasetName, inputFilePath, loadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MSAlign:
                    return new clsPHRPParserMSAlign(datasetName, inputFilePath, loadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MSGFDB:
                    return new clsPHRPParserMSGFDB(datasetName, inputFilePath, loadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.Sequest:
                    return new clsPHRPParserSequest(datasetName, inputFilePath, loadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.XTandem:
                    return new clsPHRPParserXTandem(datasetName, inputFilePath, loadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MODa:
                    return new clsPHRPParserMODa(datasetName, inputFilePath, loadModsAndSeqInfo);

                case clsPHRPReader.ePeptideHitResultType.MODPlus:
                    return new clsPHRPParserMODPlus(datasetName, inputFilePath, loadModsAndSeqInfo);

                default:
                    throw new Exception("Unrecognized value for PeptideHitResultType: " + ePeptideHitResultType.ToString());
            }
        }

        #region "Functions overridden by derived classes"

        /// <summary>
        /// Define column header names
        /// </summary>
        protected abstract void DefineColumnHeaders();

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">clsPSM object (output)</param>
        /// <returns>True if success, false if an error</returns>
        public bool ParsePHRPDataLine(string line, int linesRead, out clsPSM psm)
        {
            return ParsePHRPDataLine(line, linesRead, out psm, fastReadMode: false);
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
        public abstract bool ParsePHRPDataLine(string line, int linesRead, out clsPSM psm, bool fastReadMode);

        /// <summary>
        /// Parses the specified parameter file
        /// Also reads the Tool_Version_Info file in the same directory (if present)
        /// </summary>
        /// <param name="searchEngineParamFileName">Name of the parameter file to parse (must reside in InputDirectoryPath)</param>
        /// <param name="searchEngineParams">Search engine parameters class (output)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public abstract bool LoadSearchEngineParameters(string searchEngineParamFileName, out clsSearchEngineParameters searchEngineParams);

        #endregion

        /// <summary>
        /// Add a header column
        /// </summary>
        /// <param name="columnName"></param>
        protected void AddHeaderColumn(string columnName)
        {
            mColumnHeaders.Add(columnName, mColumnHeaders.Count);
        }

        /// <summary>
        /// Add a score to a PSM
        /// </summary>
        /// <param name="psm"></param>
        /// <param name="columns"></param>
        /// <param name="scoreColumnName"></param>
        protected void AddScore(clsPSM psm, string[] columns, string scoreColumnName)
        {
            const string NOT_FOUND = "==SCORE_NOT_FOUND==";

            var value = clsPHRPReader.LookupColumnValue(columns, scoreColumnName, mColumnHeaders, NOT_FOUND);

            if (value != NOT_FOUND)
            {
                psm.SetScore(scoreColumnName, value);
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

        private string ConvertModsToNumericMods(string cleanSequence, IReadOnlyCollection<clsAminoAcidModInfo> modifiedResidues)
        {
            mNewPeptide.Length = 0;

            if (modifiedResidues == null || modifiedResidues.Count == 0)
            {
                return cleanSequence;
            }

            for (var index = 0; index <= cleanSequence.Length - 1; index++)
            {
                mNewPeptide.Append(cleanSequence[index]);

                foreach (var modInfo in modifiedResidues)
                {
                    if (modInfo.ResidueLocInPeptide == index + 1)
                    {
                        mNewPeptide.Append(NumToStringPlusMinus(modInfo.ModDefinition.ModificationMass, 4));
                    }
                }
            }

            return mNewPeptide.ToString();
        }

        /// <summary>
        /// Look for ambiguous mods in sequenceWithMods
        /// For example, -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q
        /// </summary>
        /// <param name="sequenceWithMods"></param>
        /// <returns></returns>
        /// <remarks>List of ambiguous mods, where the keys are the start residues and the values are the ambiguous mod info</remarks>
        private SortedList<int, udtAmbiguousModInfo> ExtractAmbiguousMods(string sequenceWithMods)
        {
            if (!clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out _, out _))
            {
                primarySequence = string.Copy(sequenceWithMods);
            }

            var ambiguousMods = new SortedList<int, udtAmbiguousModInfo>();

            var residueNumber = 0;
            var parsingAmbiguousMod = false;
            var udtCurrentMod = default(udtAmbiguousModInfo);

            for (var charIndex = 0; charIndex <= primarySequence.Length - 1; charIndex++)
            {
                if (clsPHRPReader.IsLetterAtoZ(primarySequence[charIndex]))
                {
                    // Found a letter
                    residueNumber += 1;

                    if (charIndex > 0 && primarySequence[charIndex - 1] == '(')
                    {
                        // Found an ambiguous mod
                        if (!parsingAmbiguousMod)
                        {
                            parsingAmbiguousMod = true;
                            udtCurrentMod.ResidueStart = residueNumber;
                            udtCurrentMod.ResidueEnd = residueNumber;
                            udtCurrentMod.ModMassString = string.Empty;
                        }
                    }
                }
                else if (parsingAmbiguousMod)
                {
                    // Found a non-letter, and we are parsing an ambiguous mod

                    udtCurrentMod.ResidueEnd = residueNumber;
                    parsingAmbiguousMod = false;

                    // The mod mass should be next, in the form [-30.09]
                    // Parse out the mod mass
                    if (charIndex < primarySequence.Length - 2)
                    {
                        if (primarySequence[charIndex + 1] == '[')
                        {
                            var modMassString = primarySequence.Substring(charIndex + 2);
                            var bracketIndex = modMassString.IndexOf(']');
                            if (bracketIndex > 0)
                            {
                                // Valid ambiguous mod found; store it
                                modMassString = modMassString.Substring(0, bracketIndex);
                                udtCurrentMod.ModMassString = string.Copy(modMassString);

                                ambiguousMods.Add(udtCurrentMod.ResidueStart, udtCurrentMod);
                            }
                        }
                    }
                }
            }

            return ambiguousMods;
        }

        /// <summary>
        /// Finalize the PSM by updating the clean sequence, updating mod info, and updating the sequence info
        /// </summary>
        /// <param name="psm"></param>
        public void FinalizePSM(clsPSM psm)
        {
            psm.UpdateCleanSequence();

            psm.UpdateCleavageInfo(mCleavageStateCalculator);

            UpdatePSMUsingSeqInfo(psm);
        }

        private static KeyValuePair<string, string> GetMODaStaticModSetting(KeyValuePair<string, string> kvSetting, out string warningMessage)
        {
            var key = kvSetting.Key;
            var value = kvSetting.Value;

            if (!string.Equals(key, "ADD", StringComparison.CurrentCultureIgnoreCase))
            {
                throw new Exception("Key name is not ADD; this is not a MODa Static Mod Setting");
            }

            var commaIndex = value.IndexOf(',');

            if (commaIndex > 0)
            {
                var residue = value.Substring(0, commaIndex).Trim();
                value = value.Substring(commaIndex + 1).Trim();

                // Update Key to look like ADD_A or ADD_NTerm
                key = key + "_" + residue;

                kvSetting = new KeyValuePair<string, string>(key, value);
                warningMessage = string.Empty;
            }
            else
            {
                warningMessage = "Value for MODa keyword ADD does not contain a comma";
            }

            return kvSetting;
        }

        /// <summary>
        /// Report an exception as an error
        /// </summary>
        /// <param name="baseMessage"></param>
        /// <param name="ex"></param>
        protected void HandleException(string baseMessage, Exception ex)
        {
            if (string.IsNullOrEmpty(baseMessage))
            {
                baseMessage = "Error";
            }

            ReportError(baseMessage + ": " + ex.Message);
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

        /// <summary>
        /// Reads the data in modSummaryFilePath.  Populates mModInfo with the modification names, masses, and affected residues
        /// </summary>
        /// <returns>True if success; false if an error</returns>
        private void LoadModSummary()
        {
            try
            {
                var modSummaryFilePath = clsPHRPReader.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName);
                if (string.IsNullOrEmpty(modSummaryFilePath))
                {
                    ReportWarning("ModSummaryFile path is empty; unable to continue");
                    return;
                }

                modSummaryFilePath = Path.Combine(mInputDirectoryPath, modSummaryFilePath);

                modSummaryFilePath = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(modSummaryFilePath, mInputFilePath);
                var modSummaryFilePathPreferred = clsPHRPReader.AutoSwitchToFHTIfRequired(modSummaryFilePath, mInputFilePath);
                if (modSummaryFilePath != modSummaryFilePathPreferred && File.Exists(modSummaryFilePathPreferred))
                {
                    modSummaryFilePath = modSummaryFilePathPreferred;
                }

                if (!File.Exists(modSummaryFilePath))
                {
                    ReportWarning("ModSummary file not found: " + modSummaryFilePath);
                    return;
                }

                ShowMessage("Reading the PHRP ModSummary file");

                var modSummaryReader = new clsPHRPModSummaryReader(modSummaryFilePath);
                var success = modSummaryReader.Success;

                if (success)
                {
                    mModInfo = modSummaryReader.ModificationDefs;
                }
            }
            catch (Exception ex)
            {
                HandleException("Exception reading PHRP Mod Summary file", ex);
            }

        }

        private bool LoadSeqInfo()
        {
            bool success;

            var entriesParsed = 0;
            var lastProgress = DateTime.UtcNow;
            var notifyComplete = false;

            try
            {
                ShowMessage("Reading the PHRP SeqInfo file");

                // Instantiate the reader
                var reader =
                    new clsPHRPSeqMapReader(mDatasetName, mInputDirectoryPath, mPeptideHitResultType, mInputFilePath)
                    {
                        MaxProteinsPerSeqID = MaxProteinsPerPSM
                    };


                // Read the files
                success = reader.GetProteinMapping(out mResultToSeqMap, out mSeqToProteinMap, out mSeqInfo, out mPepToProteinMap);

                if (!success)
                {
                    ReportWarning(reader.ErrorMessage);
                }

                mResultIDToProteins.Clear();

                if (success)
                {
                    // Populate mResultIDToProteins

                    foreach (var mapItem in mResultToSeqMap)
                    {
                        List<string> proteinsForResultID;

                        if (mSeqToProteinMap.TryGetValue(mapItem.Value, out var proteinsForSeqID))
                        {
                            proteinsForResultID = (from protein in proteinsForSeqID select protein.ProteinName).ToList();
                            // proteinsForResultID = new List<string>(proteinsForSeqID.Count);
                            // foreach (clsProteinInfo protein in proteinsForSeqID)
                            // {
                            //     if (!proteinsForResultID.Contains(protein.ProteinName))
                            //     {
                            //         proteinsForResultID.Add(protein.ProteinName);
                            //     }
                            // }
                        }
                        else
                        {
                            proteinsForResultID = new List<string>();
                        }

                        if (MaxProteinsPerPSM > 0 && proteinsForResultID.Count > MaxProteinsPerPSM)
                        {
                            // Only add a subset of the proteins in proteinsForResultID
                            var proteinSubset = proteinsForResultID.Take(MaxProteinsPerPSM).OrderBy(item => item).ToList();
                            mResultIDToProteins.Add(mapItem.Key, proteinSubset);
                        }
                        else
                        {
                            mResultIDToProteins.Add(mapItem.Key, proteinsForResultID);
                        }

                        entriesParsed += 1;
                        if (DateTime.UtcNow.Subtract(lastProgress).TotalSeconds >= 5)
                        {
                            var pctComplete = entriesParsed / Convert.ToDouble(mResultToSeqMap.Count) * 100;
                            Console.WriteLine(" ... associating proteins with sequences: " + pctComplete.ToString("0.0") + "% complete");
                            lastProgress = DateTime.UtcNow;
                            notifyComplete = true;
                        }
                    }

                    if (notifyComplete)
                    {
                        Console.WriteLine(" ... associating proteins with sequences: 100% complete");
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error loading PHRP Seq Info", ex);
                success = false;
                if (!mInitialized)
                    throw new Exception(mErrorMessage, ex);
            }

            return success;
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
            var formatString = "+0;-0";
            if (digitsOfPrecision > 0)
            {
                formatString = "+0." + new string('0', digitsOfPrecision) + ";-0." + new string('0', digitsOfPrecision);
            }

            var valueText = value.ToString(formatString).TrimEnd('0');

            if (valueText.EndsWith("."))
            {
                // Ends in a decimal point; remove the decimal point
                valueText = valueText.TrimEnd('.');
            }

            return valueText;
        }

        /// <summary>
        /// Parse the column names in splitLine and update the local column header mapping
        /// </summary>
        /// <param name="splitLine"></param>
        /// <remarks></remarks>
        public void ParseColumnHeaders(string[] splitLine)
        {
            clsPHRPReader.ParseColumnHeaders(splitLine, mColumnHeaders);
        }

        /// <summary>
        /// Splits text on text, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
        /// </summary>
        /// <param name="text"></param>
        /// <param name="chDelimiter"></param>
        /// <returns>KeyValuePair with key and value from text; key and value will be empty if chDelimiter was not found</returns>
        /// <remarks>Automatically trims whitespace</remarks>
        public static KeyValuePair<string, string> ParseKeyValueSetting(string text, char chDelimiter)
        {
            return ParseKeyValueSetting(text, chDelimiter, string.Empty);
        }

        /// <summary>
        /// Splits text on text, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
        /// </summary>
        /// <param name="text"></param>
        /// <param name="chDelimiter"></param>
        /// <param name="commentChar">If defined, then looks for this character in the value portion of the setting and removes that character plus any text after it</param>
        /// <returns>KeyValuePair with key and value from text; key and value will be empty if chDelimiter was not found</returns>
        /// <remarks>Automatically trims whitespace</remarks>
        public static KeyValuePair<string, string> ParseKeyValueSetting(string text, char chDelimiter, string commentChar)
        {
            if (string.IsNullOrEmpty(text))
                return new KeyValuePair<string, string>(string.Empty, string.Empty);

            var charIndex = text.IndexOf(chDelimiter);
            if (charIndex <= 0)
                return new KeyValuePair<string, string>(string.Empty, string.Empty);

            var key = text.Substring(0, charIndex).Trim();
            string value;
            if (charIndex < text.Length - 1)
            {
                value = text.Substring(charIndex + 1).Trim();

                if (!string.IsNullOrEmpty(commentChar))
                {
                    // Look for the comment character
                    var commentCharIndex = value.IndexOf(commentChar, StringComparison.Ordinal);
                    if (commentCharIndex > 0)
                    {
                        // Trim off the comment
                        value = value.Substring(0, commentCharIndex).Trim();
                    }
                }
            }
            else
            {
                value = string.Empty;
            }

            var kvSetting = new KeyValuePair<string, string>(key, value);
            return kvSetting;
        }

        /// <summary>
        /// Read a Search Engine parameter file where settings are stored as key/value pairs
        /// </summary>
        /// <param name="searchEngineName">Search engine name (e.g. MSGF+)</param>
        /// <param name="searchEngineParamFileName">Search engine parameter file name (must exist in mInputDirectoryPath)</param>
        /// <param name="ePeptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
        /// <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
        /// <returns>True if success, false if an error</returns>
        protected bool ReadKeyValuePairSearchEngineParamFile(
            string searchEngineName,
            string searchEngineParamFileName,
            clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            clsSearchEngineParameters searchEngineParams)
        {
            var paramFilePath = Path.Combine(mInputDirectoryPath, searchEngineParamFileName);

            var success = ReadKeyValuePairSearchEngineParamFile(searchEngineName, paramFilePath, ePeptideHitResultType, searchEngineParams,
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

                using (var reader = new StreamReader(new FileStream(paramFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
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

        /// <summary>
        /// Determine the search engine version using a Tool_Version_Info file
        /// </summary>
        /// <param name="ePeptideHitResultType"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns></returns>
        protected bool ReadSearchEngineVersion(
            clsPHRPReader.ePeptideHitResultType ePeptideHitResultType,
            clsSearchEngineParameters searchEngineParams)
        {
            var success = false;

            try
            {
                // Read the Tool_Version_Info file to determine the analysis time and the tool version
                var toolVersionInfoFilePath = Path.Combine(mInputDirectoryPath, clsPHRPReader.GetToolVersionInfoFilename(ePeptideHitResultType));

                if (!File.Exists(toolVersionInfoFilePath) && ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.MSGFDB)
                {
                    // This could be an older MSGF+ job; check for a _MSGFDB.txt tool version file
                    var alternativeVersionInfoFilePath = Path.Combine(mInputDirectoryPath, "Tool_Version_Info_MSGFDB.txt");
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

                var searchEngineVersion = "Unknown";
                var searchDate = new DateTime(1980, 1, 1);
                var validDate = false;
                var validVersion = false;

                using (var reader = new StreamReader(new FileStream(toolVersionInfoFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        var dataLine = lineIn.TrimStart();

                        // Split the line on a colon
                        var kvSetting = ParseKeyValueSetting(dataLine, ':');

                        switch (kvSetting.Key.ToLower())
                        {
                            case "date":
                                validDate = DateTime.TryParse(kvSetting.Value, out searchDate);

                                break;
                            case "toolversioninfo":
                                if (!string.IsNullOrEmpty(kvSetting.Value))
                                {
                                    searchEngineVersion = string.Copy(kvSetting.Value);
                                    validVersion = true;
                                }
                                else
                                {
                                    // The next line contains the search engine version
                                    if (!srInFile.EndOfStream)
                                    {
                                        var searchEngineLine = srInFile.ReadLine();
                                        if (!string.IsNullOrEmpty(searchEngineLine))
                                        {
                                            searchEngineVersion = string.Copy(searchEngineLine.Trim());
                                            validVersion = true;
                                        }
                                    }
                                }
                                break;
                        }
                    }
                }

                if (!validDate)
                {
                    ReportError("Date line not found in the ToolVersionInfo file");
                }
                else if (!validVersion)
                {
                    ReportError("ToolVersionInfo line not found in the ToolVersionInfo file");
                }
                else
                {
                    success = true;
                }

                searchEngineParams.UpdateSearchEngineVersion(searchEngineVersion);
                searchEngineParams.UpdateSearchDate(searchDate);

                return success;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineVersion: " + ex.Message);
                return false;
            }

        }

        /// <summary>
        /// Report an error
        /// </summary>
        /// <param name="message"></param>
        protected void ReportError(string message)
        {
            mErrorMessage = message;
            mErrorMessages.Add(message);
            OnErrorEvent(message);
        }

        /// <summary>
        /// Report a warning
        /// </summary>
        /// <param name="message"></param>
        protected void ReportWarning(string message)
        {
            mWarningMessages.Add(message);
            OnWarningEvent(message);
        }

        /// <summary>
        /// Report a status message
        /// </summary>
        /// <param name="message"></param>
        protected void ShowMessage(string message)
        {
            OnStatusEvent(message);
        }

        private void StoreModInfo(clsPSM currentPSM, clsSeqInfo seqInfo)
        {

            var nTerminalModsAdded = new List<string>();
            var cTerminalModsAdded = new List<string>();
            var peptideResidueCount = currentPSM.PeptideCleanSequence.Length;

            currentPSM.PeptideMonoisotopicMass = seqInfo.MonoisotopicMass;

            currentPSM.ClearModifiedResidues();

            if (seqInfo.ModCount <= 0)
                return;

            // Split seqInfo.ModDescription on the comma character
            var mods = seqInfo.ModDescription.Split(',');

            if (mods.Length <= 0)
                return;

            // Parse currentPSM.Peptide to look for ambiguous mods, for example -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q

            var ambiguousMods = ExtractAmbiguousMods(currentPSM.Peptide);

            // ReSharper disable once UseImplicitlyTypedVariableEvident

            for (var modIndex = 0; modIndex <= mods.Length - 1; modIndex++)
            {
                // Split mods on the colon characters
                var kvModDetails = ParseKeyValueSetting(mods[modIndex], ':');

                if (string.IsNullOrEmpty(kvModDetails.Key) || string.IsNullOrEmpty(kvModDetails.Value))
                    continue;

                var massCorrectionTag = kvModDetails.Key;
                if (!int.TryParse(kvModDetails.Value, out var residueLoc))
                    continue;

                // Find the modification definition in mModInfo
                // Note that a given mass correction tag might be present multiple times in mModInfo, since it could be used as both a static peptide mod and a static peptide terminal mod
                // Thus, if residueLoc = 1 or residueLoc = currentPSM.PeptideCleanSequence.Length then we'll first look for a peptide or protein terminal static mod

                bool favorTerminalMods;
                clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState;

                if (residueLoc == 1)
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                    if (nTerminalModsAdded.Contains(massCorrectionTag))
                    {
                        // We have likely already added this modification as an N-terminal mod, thus, don't favor terminal mods this time
                        // An example is an iTraq peptide where there is a K at the N-terminus
                        // It gets modified with iTraq twice: once because of the N-terminus and once because of Lysine
                        // For example, R.K+144.102063+144.102063TGSY+79.9663GALAEITASK+144.102063.E
                        favorTerminalMods = false;
                    }
                    else
                    {
                        favorTerminalMods = true;
                    }
                }
                else if (residueLoc == peptideResidueCount)
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                    if (cTerminalModsAdded.Contains(massCorrectionTag))
                    {
                        favorTerminalMods = false;
                    }
                    else
                    {
                        favorTerminalMods = true;
                    }
                }
                else
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                    favorTerminalMods = false;
                }

                clsModificationDefinition matchedModDef;
                bool matchFound;

                if (mModInfo == null)
                {
                    matchedModDef = new clsModificationDefinition { MassCorrectionTag = massCorrectionTag };
                    matchFound = true;
                }
                else
                {
                    matchFound = UpdatePSMFindMatchingModInfo(massCorrectionTag, favorTerminalMods, eResidueTerminusState,
                                                                 out matchedModDef);
                }

                if (matchFound)
                {
                    var matches = (from item in ambiguousMods where item.Key == residueLoc select item.Value).ToList();

                    if (matches.Count > 0)
                    {
                        // Ambiguous modification
                        currentPSM.AddModifiedResidue(currentPSM.PeptideCleanSequence[residueLoc - 1], residueLoc,
                                                  eResidueTerminusState, matchedModDef, matches.First().ResidueEnd);
                    }
                    else
                    {
                        // Normal, non-ambiguous modified residue
                        currentPSM.AddModifiedResidue(currentPSM.PeptideCleanSequence[residueLoc - 1], residueLoc,
                                                  eResidueTerminusState, matchedModDef);
                    }

                    if (residueLoc == 1)
                    {
                        nTerminalModsAdded.Add(matchedModDef.MassCorrectionTag);
                    }
                    else if (residueLoc == peptideResidueCount)
                    {
                        cTerminalModsAdded.Add(matchedModDef.MassCorrectionTag);
                    }
                }
                else
                {
                    // Could not find a valid entry in mModInfo
                    ReportError("Unrecognized mass correction tag found in the SeqInfo file: " + massCorrectionTag);
                }
            }
        }

        /// <summary>
        /// Updates the theoretical (computed) monoisotopic mass of currentPSM using mResultToSeqMap and mSeqInfo
        /// Also updates the modification info
        /// Also updates SeqID
        /// </summary>
        /// <param name="currentPSM"></param>
        /// <returns>True if success, False if currentPSM.ResultID is not found in mResultToSeqMap</returns>
        /// <remarks></remarks>
        protected bool UpdatePSMUsingSeqInfo(clsPSM currentPSM)
        {
            var success = false;

            // First determine the modified residues present in this peptide
            if (mResultToSeqMap != null && mResultToSeqMap.Count > 0)
            {
                if (mResultToSeqMap.TryGetValue(currentPSM.ResultID, out var seqID))
                {
                    currentPSM.SeqID = seqID;

                    if (mSeqInfo.TryGetValue(seqID, out var seqInfo))
                    {
                        StoreModInfo(currentPSM, seqInfo);
                        success = true;
                    }

                    // Lookup the protein details using mSeqToProteinMap
                    if (mSeqToProteinMap.TryGetValue(seqID, out var proteinDetails))
                    {
                        foreach (var protein in proteinDetails)
                        {
                            currentPSM.AddProteinDetail(protein);
                        }
                    }

                    // Make sure all of the proteins in currentPSM.Proteins are defined in currentPSM.ProteinDetails
                    var additionalProteins1 = currentPSM.Proteins.Except(currentPSM.ProteinDetails.Keys, StringComparer.CurrentCultureIgnoreCase).ToList();

                    foreach (var proteinName in additionalProteins1)
                    {
                        if (MaxProteinsPerPSM > 0 && currentPSM.ProteinDetails.Count > MaxProteinsPerPSM)
                        {
                            // Maximum number of proteins reached (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)
                            break;
                        }

                        var proteinInfo = new clsProteinInfo(proteinName, 0, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None);
                        currentPSM.AddProtein(proteinInfo);
                    }

                    if (mPepToProteinMap.Count > 0)
                    {
                        // Make sure the residue start/end locations are up-to-date in currentPSM.ProteinDetails

                        if (mPepToProteinMap.TryGetValue(currentPSM.PeptideCleanSequence, out var pepToProteinMapInfo))
                        {
                            foreach (var protein in currentPSM.ProteinDetails)
                            {
                                // Find the matching protein in oPepToProteinMapInfo

                                if (pepToProteinMapInfo.ProteinMapInfo.TryGetValue(protein.Key, out var locations))
                                {
                                    var udtFirstLocation = locations.First();
                                    protein.Value.UpdateLocationInProtein(udtFirstLocation.ResidueStart, udtFirstLocation.ResidueEnd);
                                }
                            }
                        }
                    }
                }
            }

            if (success)
            {
                if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(currentPSM.Peptide, out _, out var prefix, out var suffix))
                {
                    currentPSM.PeptideWithNumericMods = prefix + "." + ConvertModsToNumericMods(currentPSM.PeptideCleanSequence, currentPSM.ModifiedResidues) + "." + suffix;
                }
                else
                {
                    currentPSM.PeptideWithNumericMods = ConvertModsToNumericMods(currentPSM.PeptideCleanSequence, currentPSM.ModifiedResidues);
                }
            }

            return success;
        }

        private bool UpdatePSMFindMatchingModInfo(
            string massCorrectionTag,
            bool favorTerminalMods,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            out clsModificationDefinition matchedModDef)
        {
            matchedModDef = new clsModificationDefinition();

            if (mModInfo == null)
            {
                return false;
            }


            var matchedDefs = new List<clsModificationDefinition>();

            foreach (var mod in mModInfo)
            {
                if (string.Equals(massCorrectionTag, mod.MassCorrectionTag, StringComparison.InvariantCultureIgnoreCase))
                {
                    matchedDefs.Add(mod);
                }
            }

            var matchFound = false;

            if (matchedDefs.Count > 0)
            {
                while (!matchFound)
                {
                    if (favorTerminalMods)
                    {
                        // Look for an entry in matchedDefs that is a terminal mod
                        foreach (var mod in matchedDefs)
                        {
                            if (mod.ModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod ||
                                mod.ModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod)
                            {
                                if (eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus &&
                                    (mod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) || mod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)))
                                {
                                    matchFound = true;
                                    matchedModDef = mod;
                                    break;
                                }

                                if (eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus &&
                                    (mod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS) || mod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)))
                                {
                                    matchFound = true;
                                    matchedModDef = mod;
                                    break;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Look for an entry in matchedDefs that is not a terminal mod
                        foreach (var mod in matchedDefs)
                        {
                            if (!(mod.ModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod ||
                                  mod.ModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod))
                            {
                                matchFound = true;
                                matchedModDef = mod;
                                break;
                            }
                        }
                    }

                    if (!matchFound)
                    {
                        if (favorTerminalMods)
                        {
                            favorTerminalMods = false;
                        }
                        else
                        {
                            // Still no match found (this shouldn't happen); use the first entry in matchedDefs
                            matchedModDef = matchedDefs[0];
                            matchFound = true;
                        }
                    }
                }
            }

            return matchFound;
        }
    }
}