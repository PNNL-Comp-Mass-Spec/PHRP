﻿//*********************************************************************************************************
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

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP parser base class
    /// </summary>
    public abstract class SynFileReaderBaseClass : PRISM.EventNotifier
    {
        // Ignore Spelling: MODa, iTraq, Defs

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
        protected PHRPReader.PeptideHitResultTypes mPeptideHitResultType;

        /// <summary>
        /// Modification info
        /// </summary>
        protected List<ModificationDefinition> mModInfo;

        /// <summary>
        /// Protein Names for each ResultID
        /// </summary>
        protected readonly SortedList<int, List<string>> mResultIDToProteins;

        #endregion

        #region "Properties"

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
        /// Peptide hit result type; Sequest, XTandem, Inspect, MSGFPlus, etc.
        /// </summary>
        public PHRPReader.PeptideHitResultTypes PeptideHitResultType => mPeptideHitResultType;

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
        /// <param name="peptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>If inputFilePath is an empty string, the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
        protected SynFileReaderBaseClass(string datasetName, string inputFilePath, PHRPReader.PeptideHitResultTypes peptideHitResultType,
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

            var startupOptions = new PHRPStartupOptions { LoadModsAndSeqInfo = loadModsAndSeqInfo };

            InitializeParser(datasetName, inputFilePath, peptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the parser for the given dataset, input file, and result type
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="peptideHitResultType"></param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks>If inputFilePath is an empty string, the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
        protected SynFileReaderBaseClass(string datasetName, string inputFilePath, PHRPReader.PeptideHitResultTypes peptideHitResultType, PHRPStartupOptions startupOptions)
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

            InitializeParser(datasetName, inputFilePath, peptideHitResultType, startupOptions);
        }

        /// <summary>
        /// Initialize the parser for the given dataset and input file
        /// </summary>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="peptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="startupOptions">Startup options</param>
        /// <remarks>
        /// If inputFilePath is an empty string,  the functions that solely depend on dataset name will be callable, but data related functions will not be callable
        /// startupOptions.LoadModsAndSeqInfo controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read
        /// Setting startupOptions.MaxProteinsPerPSM to a non-zero value will limit the number of proteins that are tracked
        /// </remarks>
        private void InitializeParser(string datasetName, string inputFilePath, PHRPReader.PeptideHitResultTypes peptideHitResultType, PHRPStartupOptions startupOptions)
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
                // Functions that solely require a dataset name will be callable, but cannot call functions that read a data line
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

                var phrpSynopsisName = PHRPReader.GetPHRPSynopsisFileName(mPeptideHitResultType, mDatasetName);
                var expectedSynopsisName = PHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(phrpSynopsisName, inputFile.Name);

                isSynopsisFile = string.Equals(inputFile.Name, expectedSynopsisName, StringComparison.OrdinalIgnoreCase) ||
                                 inputFile.Name.EndsWith("_syn.txt", StringComparison.OrdinalIgnoreCase);
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

                var resultToSeqMapFilename = PHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName);
                var resultToSeqMapFilePathPreferred = string.Empty;
                var resultToSeqMapFilePath = string.Empty;

                if (!string.IsNullOrEmpty(resultToSeqMapFilename))
                {
                    resultToSeqMapFilePath = PHRPReader.FindResultToSeqMapFile(InputDirectoryPath,
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
                            ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file and unable to determine the ResultToSeqMapFilename using clsPHRPReader.GetPHRPResultToSeqMapFileName()");
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
        /// Returns the appropriate PHRPParser class based on the input file name; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from inputFilePath</remarks>
        public static SynFileReaderBaseClass GetParser(string inputFilePath)
        {
            return GetParser(inputFilePath, true);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from inputFilePath</remarks>
        public static SynFileReaderBaseClass GetParser(string inputFilePath, bool loadModsAndSeqInfo)
        {
            var peptideHitResultType = PHRPReader.AutoDetermineResultType(inputFilePath);

            if (peptideHitResultType == PHRPReader.PeptideHitResultTypes.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + inputFilePath);
            }

            var datasetName = PHRPReader.AutoDetermineDatasetName(inputFilePath);
            if (string.IsNullOrEmpty(datasetName))
            {
                throw new Exception("Unable to auto-determine the Dataset Name for " + inputFilePath);
            }

            return GetParser(inputFilePath, datasetName, peptideHitResultType, loadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on the input file name
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        /// <remarks>Throws an exception if unable to auto-determine the input file type from inputFilePath</remarks>
        public static SynFileReaderBaseClass GetParser(string inputFilePath, string datasetName, bool loadModsAndSeqInfo)
        {
            var peptideHitResultType = PHRPReader.AutoDetermineResultType(inputFilePath);

            if (peptideHitResultType == PHRPReader.PeptideHitResultTypes.Unknown)
            {
                throw new Exception("Unable to auto-determine the PeptideHitResultType for " + inputFilePath);
            }

            return GetParser(inputFilePath, datasetName, peptideHitResultType, loadModsAndSeqInfo);
        }

        /// <summary>
        /// Returns the appropriate PHRPParser class based on PeptideHitResultType
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="datasetName">Dataset Name</param>
        /// <param name="peptideHitResultType">Peptide Hit Results file type</param>
        /// <param name="loadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
        public static SynFileReaderBaseClass GetParser(string inputFilePath, string datasetName, PHRPReader.PeptideHitResultTypes peptideHitResultType,
            bool loadModsAndSeqInfo)
        {
            switch (peptideHitResultType)
            {
                case PHRPReader.PeptideHitResultTypes.Inspect:
                    return new InspectSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                case PHRPReader.PeptideHitResultTypes.MSAlign:
                    return new MSAlignSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                case PHRPReader.PeptideHitResultTypes.MSGFPlus:
                    return new MSGFPlusSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                case PHRPReader.PeptideHitResultTypes.Sequest:
                    return new SequestSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                case PHRPReader.PeptideHitResultTypes.XTandem:
                    return new XTandemSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                case PHRPReader.PeptideHitResultTypes.MODa:
                    return new MODaSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                case PHRPReader.PeptideHitResultTypes.MODPlus:
                    return new MODPlusSynFileReader(datasetName, inputFilePath, loadModsAndSeqInfo);

                default:
                    throw new Exception("Unrecognized value for PeptideHitResultType: " + peptideHitResultType.ToString());
            }
        }

        #region "Functions overridden by derived classes"

        /// <summary>
        /// Define header names for the PHRP synopsis or first hits file for the given tool
        /// </summary>
        protected virtual void DefineColumnHeaders()
        {
            // Define the default column mapping
            var headerNames = GetColumnHeaderNames();

            mColumnHeaders.Clear();

            foreach (var headerName in headerNames)
            {
                AddHeaderColumn(headerName);
            }
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
        /// <param name="psm">clsPSM object (output)</param>
        /// <returns>True if successful, false if an error</returns>
        public bool ParsePHRPDataLine(string line, int linesRead, out PSM psm)
        {
            return ParsePHRPDataLine(line, linesRead, out psm, fastReadMode: false);
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
        public abstract bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode);

        /// <summary>
        /// Parses the specified parameter file
        /// Also reads the Tool_Version_Info file in the same directory (if present)
        /// </summary>
        /// <param name="searchEngineParamFileName">Name of the parameter file to parse (must reside in InputDirectoryPath)</param>
        /// <param name="searchEngineParams">Search engine parameters class (output)</param>
        /// <returns>True if successful, false if an error</returns>
        public abstract bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams);

        #endregion

        /// <summary>
        /// Add a PHRP synopsis or first hits header column to mColumnHeaders
        /// </summary>
        /// <param name="columnName"></param>
        /// <remarks>
        /// The column index will be set to mColumnHeaders.Count
        /// That value will be updated by ParseColumnHeaders
        /// </remarks>
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
        protected void AddScore(PSM psm, string[] columns, string scoreColumnName)
        {
            const string NOT_FOUND = "==SCORE_NOT_FOUND==";

            var value = PHRPReader.LookupColumnValue(columns, scoreColumnName, mColumnHeaders, NOT_FOUND);

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

        private readonly StringBuilder mNewPeptide = new StringBuilder();

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
        /// <returns>List of ambiguous mods, where the keys are the start residues and the values are the ambiguous mod info</returns>
        private SortedList<int, udtAmbiguousModInfo> ExtractAmbiguousMods(string sequenceWithMods)
        {
            if (!PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out _, out _))
            {
                primarySequence = string.Copy(sequenceWithMods);
            }

            var ambiguousMods = new SortedList<int, udtAmbiguousModInfo>();

            var residueNumber = 0;
            var parsingAmbiguousMod = false;
            var udtCurrentMod = default(udtAmbiguousModInfo);

            for (var charIndex = 0; charIndex <= primarySequence.Length - 1; charIndex++)
            {
                if (PHRPReader.IsLetterAtoZ(primarySequence[charIndex]))
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
        protected static Dictionary<T, int> GetColumnMapFromHeaderLine<T>(List<string> headerNames, SortedDictionary<string, T> headerColumnInfo)
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
                var modSummaryFileName = PHRPReader.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName);
                if (string.IsNullOrEmpty(modSummaryFileName))
                {
                    ReportWarning("ModSummaryFile name is empty; unable to continue");
                    return;
                }

                var modSummaryFilePath = PHRPReader.FindModSummaryFile(InputDirectoryPath,
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
                            var pctComplete = entriesParsed / Convert.ToDouble(ResultToSeqMap.Count) * 100;
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
        /// <param name="splitLine"></param>
        public void ParseColumnHeaders(string[] splitLine)
        {
            PHRPReader.ParseColumnHeaders(splitLine, mColumnHeaders);
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
        /// <param name="commentChar">If defined, looks for this character in the value portion of the setting and removes that character plus any text after it</param>
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
        /// <param name="searchEngineName">Search engine name (e.g. MS-GF+)</param>
        /// <param name="searchEngineParamFileName">Search engine parameter file name (must exist in InputDirectoryPath)</param>
        /// <param name="peptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
        /// <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
        /// <returns>True if successful, false if an error</returns>
        protected bool ReadKeyValuePairSearchEngineParamFile(
            string searchEngineName,
            string searchEngineParamFileName,
            PHRPReader.PeptideHitResultTypes peptideHitResultType,
            SearchEngineParameters searchEngineParams)
        {
            var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

            var success = ReadKeyValuePairSearchEngineParamFile(searchEngineName, paramFilePath, peptideHitResultType, searchEngineParams,
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
        /// <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
        /// <param name="errorMessage">Output: error message</param>
        /// <param name="warningMessage">Output: warning message</param>
        /// <returns>True if successful, false if an error</returns>
        public static bool ReadKeyValuePairSearchEngineParamFile(
            string searchEngineName,
            string paramFilePath,
            PHRPReader.PeptideHitResultTypes peptideHitResultType,
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

                        if (peptideHitResultType == PHRPReader.PeptideHitResultTypes.MODa &&
                            string.Equals(kvSetting.Key, "add", StringComparison.OrdinalIgnoreCase))
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
        /// <param name="peptideHitResultType"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        protected bool ReadSearchEngineVersion(
            PHRPReader.PeptideHitResultTypes peptideHitResultType,
            SearchEngineParameters searchEngineParams)
        {
            var success = false;

            try
            {
                // Read the Tool_Version_Info file to determine the analysis time and the tool version
                var toolVersionInfoFilePath = Path.Combine(InputDirectoryPath, PHRPReader.GetToolVersionInfoFilename(peptideHitResultType));

                if (!File.Exists(toolVersionInfoFilePath) && peptideHitResultType == PHRPReader.PeptideHitResultTypes.MSGFPlus)
                {
                    // This could be an older MS-GF+ job; check for a _MSGFDB.txt tool version file
                    var alternativeVersionInfoFilePath = Path.Combine(InputDirectoryPath, "Tool_Version_Info_MSGFDB.txt");
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

                        if (kvSetting.Key.Equals("date", StringComparison.OrdinalIgnoreCase))
                        {
                            validDate = DateTime.TryParse(kvSetting.Value, out searchDate);
                        }
                        else if (kvSetting.Key.Equals("ToolVersionInfo", StringComparison.OrdinalIgnoreCase))
                        {
                            if (!string.IsNullOrEmpty(kvSetting.Value))
                            {
                                searchEngineVersion = string.Copy(kvSetting.Value);
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
                                        searchEngineVersion = string.Copy(searchEngineLine.Trim());
                                        validVersion = true;
                                    }
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
        /// <param name="message"></param>
        protected void ReportError(string message)
        {
            mErrorMessage = message;
            ErrorMessages.Add(message);
            OnErrorEvent(message);
        }

        /// <summary>
        /// Report a warning
        /// </summary>
        /// <param name="message"></param>
        protected void ReportWarning(string message)
        {
            WarningMessages.Add(message);
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
                AminoAcidModInfo.ResidueTerminusStateConstants residueTerminusState;

                if (residueLoc == 1)
                {
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusStateConstants.PeptideNTerminus;
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
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusStateConstants.PeptideCTerminus;
                    favorTerminalMods = !cTerminalModsAdded.Contains(massCorrectionTag);
                }
                else
                {
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusStateConstants.None;
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
                    matchFound = UpdatePSMFindMatchingModInfo(massCorrectionTag, favorTerminalMods, residueTerminusState,
                                                                 out matchedModDef);
                }

                if (matchFound)
                {
                    var matches = (from item in ambiguousMods where item.Key == residueLoc select item.Value).ToList();

                    if (matches.Count > 0)
                    {
                        // Ambiguous modification
                        currentPSM.AddModifiedResidue(currentPSM.PeptideCleanSequence[residueLoc - 1], residueLoc,
                                                  residueTerminusState, matchedModDef, matches.First().ResidueEnd);
                    }
                    else
                    {
                        // Normal, non-ambiguous modified residue
                        currentPSM.AddModifiedResidue(currentPSM.PeptideCleanSequence[residueLoc - 1], residueLoc,
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

            // Make sure all of the proteins in currentPSM.Proteins are defined in currentPSM.ProteinDetails
            var additionalProteins1 = currentPSM.Proteins.Except(currentPSM.ProteinDetails.Keys, StringComparer.OrdinalIgnoreCase).ToList();

            foreach (var proteinName in additionalProteins1)
            {
                if (MaxProteinsPerPSM > 0 && currentPSM.ProteinDetails.Count > MaxProteinsPerPSM)
                {
                    // Maximum number of proteins reached (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)
                    break;
                }

                var proteinInfo = new ProteinInfo(proteinName, 0, PeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific,
                                                     PeptideCleavageStateCalculator.PeptideTerminusStateConstants.None);
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
            AminoAcidModInfo.ResidueTerminusStateConstants residueTerminusState,
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
                            if (mod.ModificationType == ModificationDefinition.ModificationTypeConstants.TerminalPeptideStaticMod ||
                                mod.ModificationType == ModificationDefinition.ModificationTypeConstants.ProteinTerminusStaticMod)
                            {
                                if (residueTerminusState == AminoAcidModInfo.ResidueTerminusStateConstants.PeptideNTerminus &&
                                    (mod.TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) || mod.TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)))
                                {
                                    matchFound = true;
                                    matchedModDef = mod;
                                    break;
                                }

                                if (residueTerminusState == AminoAcidModInfo.ResidueTerminusStateConstants.PeptideCTerminus &&
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
                            if (!(mod.ModificationType == ModificationDefinition.ModificationTypeConstants.TerminalPeptideStaticMod ||
                                  mod.ModificationType == ModificationDefinition.ModificationTypeConstants.ProteinTerminusStaticMod))
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