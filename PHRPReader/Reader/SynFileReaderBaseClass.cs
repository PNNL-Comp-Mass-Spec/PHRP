//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class is the base class for classes used to parse PHRP data lines
// It must be inherited by a class customized for the specific analysis tool (MS-GF+, X!Tandem, MSPathFinder, etc.)
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using PHRPReader.Data;

// ReSharper disable UnusedMember.Global

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP synopsis file reader base class
    /// </summary>
    public abstract class SynFileReaderBaseClass : PRISM.EventNotifier
    {
        // Ignore Spelling: Defs, iTraq, MODa, Sequest

        /// <summary>
        /// Tracks ambiguous modifications
        /// </summary>
        protected struct AmbiguousModInfo
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
            /// Modification mass (as a string) or modification name
            /// </summary>
            public string ModMassString;

            /// <summary>
            /// Show the modification mass
            /// </summary>
            public override string ToString()
            {
                return ModMassString;
            }
        }

        /// <summary>
        /// Dataset name
        /// </summary>
        protected string mDatasetName;

        /// <summary>
        /// True if initialized
        /// </summary>
        private bool mInitialized;

        /// <summary>
        /// Column headers in the synopsis file and first hits file
        /// </summary>
        protected readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Error message
        /// </summary>
        protected string mErrorMessage = string.Empty;

        /// <summary>
        /// Cleavage state calculator
        /// </summary>
        protected readonly PeptideCleavageStateCalculator mCleavageStateCalculator;

        /// <summary>
        /// Peptide mass calculator
        /// </summary>
        protected readonly PeptideMassCalculator mPeptideMassCalculator;

        /// <summary>
        /// PHRP result type
        /// </summary>
        protected PeptideHitResultTypes mPeptideHitResultType;

        /// <summary>
        /// Modification info
        /// </summary>
        protected List<ModificationDefinition> mModInfo;

        /// <summary>
        /// Protein Names for each ResultID
        /// </summary>
        protected readonly SortedList<int, List<string>> mResultIDToProteins;

        /// <summary>
        /// Cached error messages
        /// </summary>
        public List<string> ErrorMessages { get; }

        /// <summary>
        /// Input file path
        /// </summary>
        public string InputFilePath { get; private set; } = string.Empty;

        /// <summary>
        /// Input directory path
        /// </summary>
        public string InputDirectoryPath { get; private set; } = string.Empty;

        /// <summary>
        /// Input directory path
        /// </summary>
        [Obsolete("Use InputDirectoryPath")]
        public string InputFolderPath => InputDirectoryPath;

        /// <summary>
        /// Maximum number of proteins to associate with each PSM
        /// </summary>
        /// <remarks>0 means to load all proteins</remarks>
        public int MaxProteinsPerPSM { get; set; }

        /// <summary>
        /// Peptide hit result type; Sequest, XTandem, InSpecT, MSGFPlus, etc.
        /// </summary>
        public PeptideHitResultTypes PeptideHitResultType => mPeptideHitResultType;

        /// <summary>
        /// Peptide to protein map file name
        /// </summary>
        public Dictionary<string, PepToProteinMapInfo> PepToProteinMap { get; }

        /// <summary>
        /// Returns the cached mapping between ResultID and SeqID
        /// </summary>
        public SortedList<int, int> ResultToSeqMap { get; }

        /// <summary>
        /// Returns the cached sequence info, where key is SeqID
        /// </summary>
        public SortedList<int, SequenceInfo> SeqInfo { get; }

        /// <summary>
        /// Returns the cached sequence to protein map information
        /// </summary>
        public SortedList<int, List<ProteinInfo>> SeqToProteinMap { get; }

        /// <summary>
        /// Cached warning messages
        /// </summary>
        public List<string> WarningMessages { get; }

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

        /// <summary>
        /// Initialize the SynFileReader for the given dataset, input file, and result type
        /// </summary>
        /// <remarks>
        /// If inputFilePath is an empty string, the methods that solely depend on dataset name will be callable,
        /// but data related methods will not be callable
        /// </remarks>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="peptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        protected SynFileReaderBaseClass(
            string datasetName,
            string inputFilePath,
            PeptideHitResultTypes peptideHitResultType,
            bool loadModsAndSeqInfo)
        {
            ErrorMessages = new List<string>();
            WarningMessages = new List<string>();

            mResultIDToProteins = new SortedList<int, List<string>>();

            mCleavageStateCalculator = new PeptideCleavageStateCalculator();

            mPeptideMassCalculator = new PeptideMassCalculator();

            // Initialize the column mapping object
            // Using a case-insensitive comparer
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);

            // Initialize the tracking lists
            // These will get populated via the call to reader.GetProteinMapping
            ResultToSeqMap = new SortedList<int, int>();
            SeqInfo = new SortedList<int, SequenceInfo>();
            SeqToProteinMap = new SortedList<int, List<ProteinInfo>>();
            PepToProteinMap = new Dictionary<string, PepToProteinMapInfo>();

            var startupOptions = new StartupOptions { LoadModsAndSeqInfo = loadModsAndSeqInfo };

            InitializeReader(datasetName, inputFilePath, peptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the SynFileReader for the given dataset, input file, and result type
        /// </summary>
        /// <remarks>
        /// If inputFilePath is an empty string, the methods that solely depend on dataset name will be callable,
        /// but data related methods will not be callable
        /// </remarks>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="peptideHitResultType">PHRP results type</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        protected SynFileReaderBaseClass(string datasetName, string inputFilePath, PeptideHitResultTypes peptideHitResultType, StartupOptions startupOptions)
        {
            ErrorMessages = new List<string>();
            WarningMessages = new List<string>();

            mResultIDToProteins = new SortedList<int, List<string>>();

            mCleavageStateCalculator = new PeptideCleavageStateCalculator();

            mPeptideMassCalculator = startupOptions.PeptideMassCalculator ?? new PeptideMassCalculator();

            // Initialize the column mapping object
            // Using a case-insensitive comparer
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);

            // Initialize the tracking lists
            // These will get populated via the call to reader.GetProteinMapping
            ResultToSeqMap = new SortedList<int, int>();
            SeqInfo = new SortedList<int, SequenceInfo>();
            SeqToProteinMap = new SortedList<int, List<ProteinInfo>>();
            PepToProteinMap = new Dictionary<string, PepToProteinMapInfo>();

            InitializeReader(datasetName, inputFilePath, peptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the SynFileReader for the given dataset and input file
        /// </summary>
        /// <remarks>
        /// If inputFilePath is an empty string, the methods that solely depend on dataset name will be callable, but data related methods will not be callable
        /// startupOptions.LoadModsAndSeqInfo controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read
        /// Setting startupOptions.MaxProteinsPerPSM to a non-zero value will limit the number of proteins that are tracked
        /// </remarks>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="peptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="startupOptions">Startup options</param>
        private void InitializeReader(string datasetName, string inputFilePath, PeptideHitResultTypes peptideHitResultType, StartupOptions startupOptions)
        {
            if (string.IsNullOrWhiteSpace(datasetName))
                datasetName = "Undefined";

            mDatasetName = datasetName;

            mPeptideHitResultType = peptideHitResultType;

            MaxProteinsPerPSM = startupOptions.MaxProteinsPerPSM;

            bool isSynopsisFile;

            if (string.IsNullOrEmpty(inputFilePath))
            {
                // User instantiated the class without a filename
                // Methods that solely require a dataset name will be callable, but cannot call methods that read a data line

                InputFilePath = string.Empty;
                InputDirectoryPath = string.Empty;

                startupOptions.LoadModsAndSeqInfo = false;
                isSynopsisFile = false;
            }
            else
            {
                var inputFile = new FileInfo(inputFilePath);
                InputFilePath = inputFile.FullName;
                if (inputFile.Directory != null)
                {
                    InputDirectoryPath = inputFile.Directory.FullName;
                }

                if (startupOptions.DisableOpeningInputFiles)
                {
                    // File loading is disabled
                    // Methods that solely require a dataset name will be callable, but cannot call methods that read a data line

                    isSynopsisFile = false;
                }
                else
                {
                    var phrpSynopsisName = ReaderFactory.GetPHRPSynopsisFileName(mPeptideHitResultType, mDatasetName);
                    var expectedSynopsisName = ReaderFactory.AutoSwitchToLegacyMSGFDBIfRequired(phrpSynopsisName, inputFile.Name);

                    isSynopsisFile = string.Equals(inputFile.Name, expectedSynopsisName, StringComparison.OrdinalIgnoreCase) ||
                                     inputFile.Name.EndsWith("_syn.txt", StringComparison.OrdinalIgnoreCase);
                }
            }

            mErrorMessage = string.Empty;

            // Initialize the column mapping object
            // Using a case-insensitive comparer
            mColumnHeaders.Clear();

            // Initialize the tracking lists
            // These will get updated via the call to reader.GetProteinMapping
            ResultToSeqMap.Clear();
            SeqInfo.Clear();
            SeqToProteinMap.Clear();
            PepToProteinMap.Clear();

            if (startupOptions.LoadModsAndSeqInfo)
            {
                // Read the ModSummary file (if it exists)
                LoadModSummary();
            }

            if (startupOptions.LoadModsAndSeqInfo)
            {
                // Read the ResultToSeqMapInfo (if the files exist)

                var resultToSeqMapFilename = ReaderFactory.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName);
                var resultToSeqMapFilePathPreferred = string.Empty;
                var resultToSeqMapFilePath = string.Empty;

                if (!string.IsNullOrEmpty(resultToSeqMapFilename))
                {
                    resultToSeqMapFilePath = ReaderFactory.FindResultToSeqMapFile(InputDirectoryPath,
                                                                                  InputFilePath,
                                                                                  resultToSeqMapFilename,
                                                                                  out resultToSeqMapFilePathPreferred);
                }

                if (isSynopsisFile && !string.IsNullOrWhiteSpace(resultToSeqMapFilePath))
                {
                    LoadSeqInfo();
                }
                else if (!isSynopsisFile && !string.IsNullOrWhiteSpace(resultToSeqMapFilePath))
                {
                    // Processing an FHT file, and the corresponding ResultToSeqMap file does exist

                    var seqInfoLoaded = LoadSeqInfo();

                    if (!seqInfoLoaded)
                    {
                        if (string.IsNullOrEmpty(resultToSeqMapFilePathPreferred))
                        {
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file and unable to determine the ResultToSeqMapFilename using ReaderFactory.GetPHRPResultToSeqMapFileName()");
                        }
                        else
                        {
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file but the ResultToSeqMap file does not exist: " + resultToSeqMapFilePathPreferred);
                        }
                    }
                }
            }

            DefineColumnHeaders();

            mInitialized = true;
        }

        /// <summary>
        /// Returns the appropriate SynFileReader class based on the input file name; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from inputFilePath</remarks>
        /// <param name="inputFilePath">Input file path</param>
        public static SynFileReaderBaseClass GetReader(string inputFilePath)
        {
            return GetReader(inputFilePath, true);
        }

        /// <summary>
        /// Returns the appropriate SynFileReader class based on the input file name
        /// </summary>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from inputFilePath</remarks>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        public static SynFileReaderBaseClass GetReader(string inputFilePath, bool loadModsAndSeqInfo)
        {
            var peptideHitResultType = ReaderFactory.AutoDetermineResultType(inputFilePath);

            if (peptideHitResultType == PeptideHitResultTypes.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + inputFilePath);
            }

            var datasetName = ReaderFactory.AutoDetermineDatasetName(inputFilePath);

            if (string.IsNullOrEmpty(datasetName))
            {
                throw new Exception("Unable to auto-determine the Dataset Name for " + inputFilePath);
            }

            return GetReader(inputFilePath, datasetName, peptideHitResultType, loadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate SynFileReader class based on the input file name
        /// </summary>
        /// <remarks>Throws an exception if unable to auto-determine the input file type from inputFilePath</remarks>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        public static SynFileReaderBaseClass GetReader(string inputFilePath, string datasetName, bool loadModsAndSeqInfo)
        {
            var peptideHitResultType = ReaderFactory.AutoDetermineResultType(inputFilePath);

            if (peptideHitResultType == PeptideHitResultTypes.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + inputFilePath);
            }

            return GetReader(inputFilePath, datasetName, peptideHitResultType, loadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate SynFileReader class based on PeptideHitResultType
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="peptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        public static SynFileReaderBaseClass GetReader(string inputFilePath, string datasetName, PeptideHitResultTypes peptideHitResultType,
            bool loadModsAndSeqInfo)
        {
            return peptideHitResultType switch
            {
                PeptideHitResultTypes.Inspect => new InspectSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                PeptideHitResultTypes.MSAlign => new MSAlignSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                PeptideHitResultTypes.MSGFPlus => new MSGFPlusSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                PeptideHitResultTypes.Sequest => new SequestSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                PeptideHitResultTypes.XTandem => new XTandemSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                PeptideHitResultTypes.MODa => new MODaSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                PeptideHitResultTypes.MODPlus => new MODPlusSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo),
                _ => throw new Exception("Unrecognized value for PeptideHitResultType: " + peptideHitResultType)
            };
        }

        /// <summary>
        /// Define header names for the PHRP synopsis or first hits file for the given tool
        /// </summary>
        protected virtual void DefineColumnHeaders()
        {
            // Define the default column mapping
            var headerNames = GetColumnHeaderNames();

            DefineColumnHeaders(mColumnHeaders, headerNames);
        }

        /// <summary>
        /// List of header names for the PHRP synopsis or first hits file for the given tool
        /// </summary>
        protected abstract List<string> GetColumnHeaderNames();

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM info</param>
        /// <returns>True if successful, false if an error</returns>
        public bool ParsePHRPDataLine(string line, int linesRead, out PSM psm)
        {
            return ParsePHRPDataLine(line, linesRead, out psm, fastReadMode: false);
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM info</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        public abstract bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode);

        /// <summary>
        /// Parses the specified parameter file
        /// Also reads the Tool_Version_Info file in the same directory (if present)
        /// </summary>
        /// <param name="searchEngineParamFileName">Name of the parameter file to parse (must reside in InputDirectoryPath)</param>
        /// <param name="searchEngineParams">Output: Search engine parameters class</param>
        /// <returns>True if successful, false if an error</returns>
        public abstract bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams);

        /// <summary>
        /// Add a score to a PSM
        /// </summary>
        /// <param name="psm">PSM</param>
        /// <param name="columns">column data</param>
        /// <param name="scoreColumnName">Score column name</param>
        protected void AddScore(PSM psm, string[] columns, string scoreColumnName)
        {
            const string NOT_FOUND = "==SCORE_NOT_FOUND==";

            var value = ReaderFactory.LookupColumnValue(columns, scoreColumnName, mColumnHeaders, NOT_FOUND);

            if (value != NOT_FOUND)
            {
                psm.SetScore(scoreColumnName, value);
            }
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

        private readonly StringBuilder mNewPeptide = new();

        private string ConvertModsToNumericMods(string cleanSequence, IReadOnlyCollection<AminoAcidModInfo> modifiedResidues)
        {
            mNewPeptide.Length = 0;

            if (modifiedResidues == null || modifiedResidues.Count == 0)
            {
                return cleanSequence;
            }

            for (var index = 0; index <= cleanSequence.Length - 1; index++)
            {
                mNewPeptide.Append(cleanSequence[index]);

                foreach (var modDef in modifiedResidues)
                {
                    if (modDef.ResidueLocInPeptide == index + 1)
                    {
                        mNewPeptide.Append(NumToStringPlusMinus(modDef.ModDefinition.ModificationMass, 4));
                    }
                }
            }

            return mNewPeptide.ToString();
        }

        /// <summary>
        /// Define the default column mapping
        /// </summary>
        public static void DefineColumnHeaders(SortedDictionary<string, int> columnHeaders, List<string> headerNames)
        {
            columnHeaders.Clear();

            foreach (var headerName in headerNames)
            {
                columnHeaders.Add(headerName, columnHeaders.Count);
            }
        }

        /// <summary>
        /// Look for ambiguous mods in sequenceWithMods
        /// For example, -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q
        /// </summary>
        /// <param name="sequenceWithMods">Peptide sequence, with modification masses</param>
        /// <returns>List of ambiguous mods, where the keys are the start residues and the values are the ambiguous mod info</returns>
        private SortedList<int, AmbiguousModInfo> ExtractAmbiguousMods(string sequenceWithMods)
        {
            if (!PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out _, out _))
            {
                primarySequence = sequenceWithMods;
            }

            var ambiguousMods = new SortedList<int, AmbiguousModInfo>();

            var residueNumber = 0;
            var parsingAmbiguousMod = false;
            var udtCurrentMod = default(AmbiguousModInfo);

            for (var charIndex = 0; charIndex <= primarySequence.Length - 1; charIndex++)
            {
                if (ReaderFactory.IsLetterAtoZ(primarySequence[charIndex]))
                {
                    // Found a letter
                    residueNumber++;

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

                    // The mod mass should be next, in the form [-30.09] or [Acetyl]
                    // Parse out the mod mass or mod name
                    if (charIndex < primarySequence.Length - 2)
                    {
                        if (primarySequence[charIndex + 1] == '[')
                        {
                            var bracketIndex = primarySequence.IndexOf(']', charIndex + 2);

                            if (bracketIndex > 0)
                            {
                                // Valid ambiguous mod found; store it
                                udtCurrentMod.ModMassString = primarySequence.Substring(charIndex + 2, bracketIndex - charIndex - 2);

                                ambiguousMods.Add(udtCurrentMod.ResidueStart, udtCurrentMod);

                                charIndex = bracketIndex;
                            }
                            else
                            {
                                OnWarningEvent("Opening bracket at index {0} does not have a matching closing bracket: {1}",
                                    charIndex + 1, primarySequence);
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
        /// <param name="psm">PSM</param>
        public void FinalizePSM(PSM psm)
        {
            psm.UpdateCleanSequence();

            psm.UpdateCleavageInfo(mCleavageStateCalculator);

            UpdatePSMUsingSeqInfo(psm);
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by headerColumnInfo
        /// Populates a dictionary mapping a PHRPReader enum to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames">List of column names from the header line of a data file</param>
        /// <param name="headerColumnInfo">Dictionary mapping standard header column names to a PHRPReader enum (e.g. MSGFPlusSynFileColumns)</param>
        /// <returns>Dictionary mapping the PHRPReader enum value to the column index in headerNames (0-based column index)</returns>
        internal static Dictionary<T, int> GetColumnMapFromHeaderLine<T>(List<string> headerNames, SortedDictionary<string, T> headerColumnInfo)
        {
            var columnNameToIndexMap = new Dictionary<T, int>();

            for (var index = 0; index < headerNames.Count; index++)
            {
                var headerName = headerNames[index];

                foreach (var headerColumn in headerColumnInfo)
                {
                    if (!headerName.Equals(headerColumn.Key, StringComparison.OrdinalIgnoreCase))
                        continue;

                    if (columnNameToIndexMap.ContainsKey(headerColumn.Value))
                    {
                        // Header name is present more than once; ignore the duplicate entry
                        break;
                    }

                    columnNameToIndexMap.Add(headerColumn.Value, index);
                    break;
                }
            }

            return columnNameToIndexMap;
        }

        private static KeyValuePair<string, string> GetMODaStaticModSetting(KeyValuePair<string, string> kvSetting, out string warningMessage)
        {
            var key = kvSetting.Key;
            var value = kvSetting.Value;

            if (!string.Equals(key, "ADD", StringComparison.OrdinalIgnoreCase))
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
        /// <param name="baseMessage">Base message</param>
        /// <param name="ex">Exception</param>
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
        /// <param name="data">Value to parse</param>
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
        /// Reads the data in modSummaryFilePath. Populates mModInfo with the modification names, masses, and affected residues
        /// </summary>
        /// <returns>True if success; false if an error</returns>
        private void LoadModSummary()
        {
            try
            {
                var modSummaryFileName = ReaderFactory.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName);

                if (string.IsNullOrEmpty(modSummaryFileName))
                {
                    ReportWarning("ModSummaryFile name is empty; unable to continue");
                    return;
                }

                var modSummaryFilePath = ReaderFactory.FindModSummaryFile(InputDirectoryPath,
                                                                          InputFilePath,
                                                                          modSummaryFileName,
                                                                          out var modSummaryFileNamePreferred);

                if (string.IsNullOrWhiteSpace(modSummaryFilePath) || !File.Exists(modSummaryFilePath))
                {
                    // ModSummary file not found, expected name: Dataset_msgfplus_syn_ModSummary.txt
                    ReportWarning("ModSummary file not found: " + modSummaryFileNamePreferred);
                    return;
                }

                ShowMessage("Reading the PHRP ModSummary file");

                var modSummaryReader = new PHRPModSummaryReader(modSummaryFilePath);
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
                    new PHRPSeqMapReader(mDatasetName, InputDirectoryPath, mPeptideHitResultType, InputFilePath)
                    {
                        MaxProteinsPerSeqID = MaxProteinsPerPSM
                    };

                // Read the files
                success = reader.GetProteinMapping(ResultToSeqMap, SeqToProteinMap, SeqInfo, PepToProteinMap);

                if (!success)
                {
                    ReportWarning(reader.ErrorMessage);
                }

                mResultIDToProteins.Clear();

                if (success)
                {
                    // Populate mResultIDToProteins

                    foreach (var mapItem in ResultToSeqMap)
                    {
                        List<string> proteinsForResultID;

                        if (SeqToProteinMap.TryGetValue(mapItem.Value, out var proteinsForSeqID))
                        {
                            proteinsForResultID = (from protein in proteinsForSeqID select protein.ProteinName).ToList();
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

                        entriesParsed++;
                        if (DateTime.UtcNow.Subtract(lastProgress).TotalSeconds >= 5)
                        {
                            var pctComplete = entriesParsed / (double)ResultToSeqMap.Count * 100;
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
        /// <param name="value">Number to format</param>
        /// <param name="digitsOfPrecision">Digits of precision</param>
        /// <returns>Formatted number</returns>
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
        /// <param name="splitLine">Array of column names</param>
        public void ParseColumnHeaders(string[] splitLine)
        {
            ReaderFactory.ParseColumnHeaders(splitLine, mColumnHeaders);
        }

        /// <summary>
        /// Convert a list of delimited items to an enumerable list, trimming whitespace from each item
        /// </summary>
        /// <param name="delimitedList">Delimited list of items</param>
        /// <param name="delimiter">Delimiter</param>
        /// <param name="includeEmptyItems">If true, include empty items in the list</param>
        /// <remarks>
        /// If delimitedList is "Value1;;Value2" and includeEmptyItems is false, a 2 item list will be returned
        /// If delimitedList is "Value1;;Value2" and includeEmptyItems is true,  a 3 item list will be returned
        /// </remarks>
        /// <returns>List of items</returns>
        public static IEnumerable<string> ParseDelimitedList(string delimitedList, char delimiter, bool includeEmptyItems = false)
        {
            return includeEmptyItems
                ? delimitedList.Split(delimiter).Select(item => item.Trim()).ToList()
                : (from item in delimitedList.Split(delimiter) where !string.IsNullOrWhiteSpace(item) select item.Trim()).ToList();
        }

        /// <summary>
        /// Splits text on text, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
        /// </summary>
        /// <remarks>Automatically trims whitespace</remarks>
        /// <param name="text">Text to parse</param>
        /// <param name="chDelimiter">Delimiter</param>
        /// <returns>KeyValuePair with key and value from text; key and value will be empty if chDelimiter was not found</returns>
        public static KeyValuePair<string, string> ParseKeyValueSetting(string text, char chDelimiter)
        {
            return ParseKeyValueSetting(text, chDelimiter, string.Empty);
        }

        /// <summary>
        /// Splits text on text, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
        /// </summary>
        /// <remarks>Automatically trims whitespace</remarks>
        /// <param name="text">Text to parse</param>
        /// <param name="chDelimiter">Delimiter</param>
        /// <param name="commentChar">If defined, looks for this character in the value portion of the setting and removes that character plus any text after it</param>
        /// <returns>KeyValuePair with key and value from text; key and value will be empty if chDelimiter was not found</returns>
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

            return new KeyValuePair<string, string>(key, value);
        }

        /// <summary>
        /// Read a Search Engine parameter file where settings are stored as key/value pairs
        /// </summary>
        /// <param name="searchEngineName">Search engine name (e.g. MS-GF+)</param>
        /// <param name="searchEngineParamFileName">Search engine parameter file name (must exist in InputDirectoryPath)</param>
        /// <param name="peptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
        /// <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling method)</param>
        /// <returns>True if successful, false if an error</returns>
        protected bool ReadKeyValuePairSearchEngineParamFile(
            string searchEngineName,
            string searchEngineParamFileName,
            PeptideHitResultTypes peptideHitResultType,
            SearchEngineParameters searchEngineParams)
        {
            var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

            var success = ReadKeyValuePairSearchEngineParamFile(
                searchEngineName, paramFilePath, peptideHitResultType, searchEngineParams,
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
        /// <param name="peptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
        /// <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling method)</param>
        /// <param name="errorMessage">Output: error message</param>
        /// <param name="warningMessage">Output: warning message</param>
        /// <returns>True if successful, false if an error</returns>
        public static bool ReadKeyValuePairSearchEngineParamFile(
            string searchEngineName,
            string paramFilePath,
            PeptideHitResultTypes peptideHitResultType,
            SearchEngineParameters searchEngineParams,
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

                using var reader = new StreamReader(new FileStream(paramFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                // The delimiter is typically an equals sign, but for InSpecT it is a comma
                var keyValueDelimiter = peptideHitResultType == PeptideHitResultTypes.Inspect ? ',' : '=';

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    var dataLine = lineIn.TrimStart();

                    if (string.IsNullOrWhiteSpace(dataLine) || dataLine.StartsWith("#") || !dataLine.Contains(keyValueDelimiter))
                    {
                        continue;
                    }

                    // Split the line on the equals sign (for InSpecT, split on the first comma)
                    var kvSetting = ParseKeyValueSetting(dataLine, keyValueDelimiter, "#");

                    if (string.IsNullOrEmpty(kvSetting.Key))
                    {
                        continue;
                    }

                    if (peptideHitResultType == PeptideHitResultTypes.MODa &&
                        string.Equals(kvSetting.Key, "add", StringComparison.OrdinalIgnoreCase))
                    {
                        // ModA defines all of its static modifications with the ADD keyword
                        // Split the value at the comma and create a new setting entry with the residue name
                        kvSetting = GetMODaStaticModSetting(kvSetting, out var warningMessageAddon);
                        if (!string.IsNullOrWhiteSpace(warningMessageAddon))
                        {
                            if (string.IsNullOrWhiteSpace(warningMessage))
                            {
                                warningMessage = warningMessageAddon;
                            }
                            else if (!warningMessage.Contains(warningMessageAddon))
                            {
                                warningMessage += "; " + warningMessageAddon;
                            }
                        }
                    }

                    searchEngineParams.AddUpdateParameter(kvSetting);
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
        /// <param name="peptideHitResultType">PHRP result type</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <returns>True if successful, false if an error</returns>
        protected bool ReadSearchEngineVersion(
            PeptideHitResultTypes peptideHitResultType,
            SearchEngineParameters searchEngineParams)
        {
            var success = false;

            try
            {
                // Read the Tool_Version_Info file to determine the analysis time and the tool version
                var toolVersionInfoFilePath = string.Empty;

                var toolVersionInfoFilenames = ReaderFactory.GetToolVersionInfoFilenames(peptideHitResultType);

                if (toolVersionInfoFilenames.Count == 0)
                {
                    ReportWarning("GetToolVersionInfoFilenames returned an empty list for result type " + peptideHitResultType);
                    return false;
                }

                foreach (var toolVersionInfoFile in toolVersionInfoFilenames)
                {
                    toolVersionInfoFilePath = Path.Combine(InputDirectoryPath, toolVersionInfoFile);

                    if (File.Exists(toolVersionInfoFilePath) || peptideHitResultType != PeptideHitResultTypes.MSGFPlus)
                    {
                        break;
                    }

                    // This could be an older MS-GF+ job; check for a _MSGFDB.txt tool version file
                    var alternativeVersionInfoFilePath = Path.Combine(InputDirectoryPath, "Tool_Version_Info_MSGFDB.txt");

                    // ReSharper disable once InvertIf
                    if (File.Exists(alternativeVersionInfoFilePath))
                    {
                        toolVersionInfoFilePath = alternativeVersionInfoFilePath;
                        break;
                    }
                }

                if (string.IsNullOrWhiteSpace(toolVersionInfoFilePath) || !File.Exists(toolVersionInfoFilePath))
                {
                    ReportWarning("Tool version info file not found: " + toolVersionInfoFilePath);
                    return false;
                }

                var searchEngineVersion = "Unknown";
                var searchDate = new DateTime(1980, 1, 1);
                var validDate = false;
                var validVersion = false;

                using var reader = new StreamReader(new FileStream(toolVersionInfoFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    var dataLine = lineIn.TrimStart();

                    // Split the line on a colon
                    var kvSetting = ParseKeyValueSetting(dataLine, ':');

                    if (kvSetting.Key.Equals("date", StringComparison.OrdinalIgnoreCase))
                    {
                        validDate = DateTime.TryParse(kvSetting.Value, out searchDate);
                    }
                    else if (kvSetting.Key.Equals("ToolVersionInfo", StringComparison.OrdinalIgnoreCase))
                    {
                        if (!string.IsNullOrEmpty(kvSetting.Value))
                        {
                            searchEngineVersion = kvSetting.Value;
                            validVersion = true;
                        }
                        else
                        {
                            // The next line contains the search engine version
                            if (!reader.EndOfStream)
                            {
                                var searchEngineLine = reader.ReadLine();

                                if (!string.IsNullOrEmpty(searchEngineLine))
                                {
                                    searchEngineVersion = searchEngineLine.Trim();
                                    validVersion = true;
                                }
                            }
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
        /// <param name="message">Error message</param>
        protected void ReportError(string message)
        {
            mErrorMessage = message;
            ErrorMessages.Add(message);
            OnErrorEvent(message);
        }

        /// <summary>
        /// Report a warning
        /// </summary>
        /// <param name="message">Warning message</param>
        protected void ReportWarning(string message)
        {
            WarningMessages.Add(message);
            OnWarningEvent(message);
        }

        /// <summary>
        /// Report a status message
        /// </summary>
        /// <param name="message">Status message</param>
        protected void ShowMessage(string message)
        {
            OnStatusEvent(message);
        }

        private void StoreModInfo(PSM currentPSM, SequenceInfo seqInfo)
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

            if (mods.Length == 0)
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
                // Thus, if residueLoc = 1 or residueLoc = currentPSM.PeptideCleanSequence.Length, we'll first look for a peptide or protein terminal static mod

                bool favorTerminalMods;
                AminoAcidModInfo.ResidueTerminusState residueTerminusState;

                if (residueLoc == 1)
                {
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
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
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                    favorTerminalMods = !cTerminalModsAdded.Contains(massCorrectionTag);
                }
                else
                {
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusState.None;
                    favorTerminalMods = false;
                }

                ModificationDefinition matchedModDef;
                bool matchFound;

                if (mModInfo == null)
                {
                    matchedModDef = new ModificationDefinition { MassCorrectionTag = massCorrectionTag };
                    matchFound = true;
                }
                else
                {
                    matchFound = UpdatePSMFindMatchingModInfo(
                        massCorrectionTag, favorTerminalMods,
                        residueTerminusState, out matchedModDef);
                }

                if (!matchFound)
                {
                    // Could not find a valid entry in mModInfo
                    matchedModDef = new ModificationDefinition
                    {
                        MassCorrectionTag = massCorrectionTag
                    };

                    switch (matchedModDef.MassCorrectionTag)
                    {
                        case "MinusH2O":
                            matchedModDef.ModificationMass = -18.010565;
                            break;

                        default:
                            ReportError("Unrecognized mass correction tag found in the SeqInfo file: " + massCorrectionTag);
                            break;
                    }
                }

                var matches = (from item in ambiguousMods where item.Key == residueLoc select item.Value).ToList();

                if (matches.Count > 0)
                {
                    // Ambiguous modification
                    currentPSM.AddModifiedResidue(
                        currentPSM.PeptideCleanSequence[residueLoc - 1], residueLoc,
                        residueTerminusState, matchedModDef, matches.First().ResidueEnd);
                }
                else
                {
                    // Normal, non-ambiguous modified residue
                    currentPSM.AddModifiedResidue(
                        currentPSM.PeptideCleanSequence[residueLoc - 1], residueLoc,
                        residueTerminusState, matchedModDef);
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
        }

        /// <summary>
        /// Updates the theoretical (computed) monoisotopic mass of currentPSM using mResultToSeqMap and mSeqInfo
        /// Also updates the modification info
        /// Also updates SeqID
        /// </summary>
        /// <param name="currentPSM">Current PSM</param>
        /// <returns>True if successful, False if currentPSM.ResultID is not found in mResultToSeqMap</returns>
        protected bool UpdatePSMUsingSeqInfo(PSM currentPSM)
        {
            if (ResultToSeqMap == null || ResultToSeqMap.Count == 0)
            {
                return false;
            }

            // First determine the modified residues present in this peptide
            if (!ResultToSeqMap.TryGetValue(currentPSM.ResultID, out var seqID))
            {
                return false;
            }

            currentPSM.SeqID = seqID;

            if (SeqInfo.TryGetValue(seqID, out var seqInfo))
            {
                StoreModInfo(currentPSM, seqInfo);
            }

            // Lookup the protein details using mSeqToProteinMap
            if (SeqToProteinMap.TryGetValue(seqID, out var proteinDetails))
            {
                foreach (var protein in proteinDetails)
                {
                    currentPSM.AddProteinDetail(protein);
                }
            }

            // Make sure all the proteins in currentPSM.Proteins are defined in currentPSM.ProteinDetails
            foreach (var proteinName in currentPSM.Proteins.Except(currentPSM.ProteinDetails.Keys, StringComparer.OrdinalIgnoreCase).ToList())
            {
                if (MaxProteinsPerPSM > 0 && currentPSM.ProteinDetails.Count > MaxProteinsPerPSM)
                {
                    // Maximum number of proteins reached (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)
                    break;
                }

                var proteinInfo = new ProteinInfo(proteinName, 0, PeptideCleavageStateCalculator.PeptideCleavageState.NonSpecific,
                                                     PeptideCleavageStateCalculator.PeptideTerminusState.None);
                currentPSM.AddProtein(proteinInfo);
            }

            if (PepToProteinMap.Count > 0)
            {
                // Make sure the residue start/end locations are up-to-date in currentPSM.ProteinDetails

                if (PepToProteinMap.TryGetValue(currentPSM.PeptideCleanSequence, out var pepToProteinMapInfo))
                {
                    foreach (var protein in currentPSM.ProteinDetails)
                    {
                        // Find the matching protein in pepToProteinMapInfo
                        if (pepToProteinMapInfo.ProteinMapInfo.TryGetValue(protein.Key, out var locations))
                        {
                            var udtFirstLocation = locations.First();
                            protein.Value.UpdateLocationInProtein(udtFirstLocation.ResidueStart, udtFirstLocation.ResidueEnd);
                        }
                    }
                }
            }

            if (PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(currentPSM.Peptide, out _, out var prefix, out var suffix))
            {
                currentPSM.PeptideWithNumericMods = prefix + "." + ConvertModsToNumericMods(currentPSM.PeptideCleanSequence, currentPSM.ModifiedResidues) + "." + suffix;
            }
            else
            {
                currentPSM.PeptideWithNumericMods = ConvertModsToNumericMods(currentPSM.PeptideCleanSequence, currentPSM.ModifiedResidues);
            }

            return true;
        }

        private bool UpdatePSMFindMatchingModInfo(
            string massCorrectionTag,
            bool favorTerminalMods,
            AminoAcidModInfo.ResidueTerminusState residueTerminusState,
            out ModificationDefinition matchedModDef)
        {
            matchedModDef = new ModificationDefinition();

            if (mModInfo == null)
            {
                return false;
            }

            var matchedDefs = new List<ModificationDefinition>();

            foreach (var mod in mModInfo)
            {
                if (string.Equals(massCorrectionTag, mod.MassCorrectionTag, StringComparison.OrdinalIgnoreCase))
                {
                    matchedDefs.Add(mod);
                }
            }

            if (matchedDefs.Count == 0 && massCorrectionTag.Equals("MinusH2O"))
            {
                // ReSharper disable once ForeachCanBePartlyConvertedToQueryUsingAnotherGetEnumerator
                foreach (var mod in mModInfo)
                {
                    if (Math.Abs(mod.ModificationMass + 18.010565) > 0.001)
                        continue;

                    OnWarningEvent("Mod {0} not found in mModInfo by name, but was found by modification mass");
                    matchedDefs.Add(mod);
                    break;
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
                            if (mod.ModificationType is ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod or ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod)
                            {
                                if (residueTerminusState == AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus &&
                                    (mod.TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) || mod.TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)))
                                {
                                    matchFound = true;
                                    matchedModDef = mod;
                                    break;
                                }

                                if (residueTerminusState == AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus &&
                                    (mod.TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS) || mod.TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)))
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
                            if (!(mod.ModificationType is ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod or ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod))
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
