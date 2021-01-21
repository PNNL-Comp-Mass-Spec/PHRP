// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or https://www.pnnl.gov/sysbio/ or https://panomics.pnnl.gov/
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Threading;
using PeptideToProteinMapEngine;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class can be used as a base class for peptide hit results processor classes
    /// </summary>
    public abstract class clsPHRPBaseClass : PRISM.EventNotifier
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        protected clsPHRPBaseClass()
        {
            FileDate = "October 25, 2020";

            mPeptideSeqMassCalculator = new clsPeptideMassCalculator { ChargeCarrierMass = clsPeptideMassCalculator.MASS_PROTON };

            // Initialize mPeptideMods
            mPeptideMods = new clsPeptideModificationContainer();

            // Initialize mUniqueSequences
            mUniqueSequences = new clsUniqueSequencesContainer();

            // Initialize mSeqToProteinMap
            mSeqToProteinMap = new SortedSet<string>();

            // Define a RegEx to replace all of the non-letter characters
            mReplaceSymbols = new Regex("[^A-Za-z]", RegexOptions.Compiled);

            mProteinNameOrder = new Dictionary<string, int>();

            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        private const string UNIQUE_SEQ_TO_PROTEIN_MAP_SEP = "_";

        private const string COLUMN_NAME_UNIQUE_SEQ_ID = "Unique_Seq_ID";
        private const string COLUMN_NAME_PROTEIN_NAME = "Protein_Name";
        protected const string COLUMN_NAME_RESULTID = "ResultID";
        protected const string COLUMN_NAME_PEPTIDE = "Peptide";
        private const string COLUMN_NAME_RESIDUE = "Residue";
        private const string COLUMN_NAME_PROTEIN_RESIDUE_NUMBER = "Protein_Residue_Num";
        private const string COLUMN_NAME_RESIDUE_MOD_NAME = "Mod_Name";
        private const string COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER = "Peptide_Residue_Num";
        private const string COLUMN_NAME_MSGF_SPECPROB = "MSGF_SpecProb";

        public const string XTANDEM_RESULTS_FILE_SUFFIX = "_xt.xml";
        public const string SEQUEST_SYNOPSIS_FILE_SUFFIX = "_syn.txt";
        public const string SEQUEST_FIRST_HITS_FILE_SUFFIX = "_fht.txt";

        public const string INSPECT_RESULTS_FILE_SUFFIX = "_inspect.txt";
        public const string INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX = "_fht.txt";
        public const string INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX = "_Fscore_fht.txt";

        public const string MSGFDB_RESULTS_FILE_SUFFIX = "_msgfdb.txt";

        public const string MSALIGN_RESULTS_FILE_SUFFIX = "_MSAlign_ResultTable.txt";

        public const string MODa_RESULTS_FILE_SUFFIX = "_moda.id.txt";
        public const string MODPlus_RESULTS_FILE_SUFFIX = "_modp.id.txt";
        public const string MSPathFinder_RESULTS_FILE_SUFFIX = "_IcTda.tsv";
        public const string TopPIC_RESULTS_FILE_SUFFIX = "_TopPIC_PrSMs.txt";

        public const string FILENAME_SUFFIX_RESULT_TO_SEQ_MAP = "_ResultToSeqMap.txt";
        public const string FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP = "_SeqToProteinMap.txt";

        public const string FILENAME_SUFFIX_SEQ_INFO = "_SeqInfo.txt";
        public const string FILENAME_SUFFIX_MOD_DETAILS = "_ModDetails.txt";
        public const string FILENAME_SUFFIX_MOD_SUMMARY = "_ModSummary.txt";

        public const string FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING = "_PepToProtMap";
        public const string FILENAME_SUFFIX_PROTEIN_MODS = "_ProteinMods.txt";
        public const string FILENAME_SUFFIX_MSGF = "_MSGF.txt";

        protected const float PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE = 90;
        private const float PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE = 95;

        private const string PROTEIN_NAME_NO_MATCH = "__NoMatch__";

        public enum ePeptideHitResultsFileFormatConstants
        {
            AutoDetermine = 0,

            /// <summary>
            /// SEQUEST synopsis hits file
            /// </summary>
            SequestSynopsisFile = 1,

            /// <summary>
            /// SEQUEST first hits file
            /// </summary>
            SequestFirstHitsFile = 2,

            /// <summary>
            /// X!Tandem
            /// </summary>
            XTandemXMLFile = 3,

            /// <summary>
            /// Inspect
            /// </summary>
            InspectTXTFile = 4,

            /// <summary>
            /// MSGFDB, MS-GF+ (MSGF+)
            /// </summary>
            MSGFPlusTXTFile = 5,

            /// <summary>
            /// MSAlign
            /// </summary>
            MSAlignTXTFile = 6,

            /// <summary>
            /// MODa
            /// </summary>
            MODaTXTFile = 7,

            /// <summary>
            /// MODPlus
            /// </summary>
            MODPlusTXTFile = 8,

            /// <summary>
            /// MSPathFinder
            /// </summary>
            MSPathFinderTSVFile = 9,

            /// <summary>
            /// TopPIC
            /// </summary>
            // ReSharper disable once IdentifierTypo
            TopPICTXTFile = 10
        }

        public enum ePHRPErrorCodes
        {
            NoError = 0,
            InvalidInputFilePath = 1,
            InvalidOutputDirectoryPath = 2,
            ParameterFileNotFound = 3,
            MassCorrectionTagsFileNotFound = 4,
            ModificationDefinitionFileNotFound = 5,

            ErrorReadingInputFile = 6,
            ErrorCreatingOutputFiles = 7,
            ErrorReadingParameterFile = 8,
            ErrorReadingMassCorrectionTagsFile = 9,
            ErrorReadingModificationDefinitionsFile = 10,

            FilePathError = 11,
            UnspecifiedError = -1
        }

        #endregion

        #region "Structures"
        protected struct udtSearchOptionModificationInfoType
        {
            public int SortOrder;
            public double ModificationMass;
            public string TargetResidues;
            public clsModificationDefinition.eModificationTypeConstants ModificationType;

            public override string ToString()
            {
                return ModificationType + ": " + ModificationMass + " @ " + TargetResidues;
            }
        }

        internal struct udtModNameAndResidueLocType
        {
            public string ModName;
            public int ResidueLocInPeptide;

            public override string ToString()
            {
                return ResidueLocInPeptide + ": " + ModName;
            }
        }

        protected struct udtPepToProteinMappingType
        {
            public string Peptide;
            public string Protein;
            public int ResidueStart;
            public int ResidueEnd;

            public override string ToString()
            {
                return Peptide + ", " + Protein;
            }
        }

        #endregion

        #region "Class wide Variables"

        protected ePHRPErrorCodes mErrorCode = ePHRPErrorCodes.NoError;
        protected string mErrorMessage = string.Empty;

        protected readonly clsPeptideMassCalculator mPeptideSeqMassCalculator;

        protected readonly clsPeptideModificationContainer mPeptideMods;
        private readonly clsUniqueSequencesContainer mUniqueSequences;
        private readonly SortedSet<string> mSeqToProteinMap;

        private StreamWriter mResultToSeqMapFile;
        private StreamWriter mSeqInfoFile;
        private StreamWriter mModDetailsFile;
        private StreamWriter mSeqToProteinMapFile;

        private int mNextPeptideToProteinMapperLevel;

        /// <summary>
        /// Tracks the protein names in the order that they are listed in the FASTA file
        /// Keys are protein Names, values are a sequentially assigned integer
        /// </summary>
        protected readonly Dictionary<string, int> mProteinNameOrder;

        private readonly Regex mReplaceSymbols;

        #endregion

        #region "Progress Events and Variables"

        public event ProgressResetEventHandler ProgressReset;
        public delegate void ProgressResetEventHandler();

        public event ProgressCompleteEventHandler ProgressComplete;
        public delegate void ProgressCompleteEventHandler();

        protected string mProgressStepDescription = string.Empty;

        /// <summary>
        /// Ranges from 0 to 100, but can contain decimal percentage values
        /// </summary>
        protected float mProgressPercentComplete;
        #endregion

        #region "Properties"

        public bool AbortProcessing { get; set; }

        public bool CreateModificationSummaryFile { get; set; }

        public bool CreateInspectFirstHitsFile { get; set; }

        public bool CreateInspectSynopsisFile { get; set; }

        /// <summary>
        /// Create protein mods file
        /// </summary>
        /// <returns></returns>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the Fasta file specified by mFastaFilePath</remarks>
        public bool CreateProteinModsFile { get; set; }

        public clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType EnzymeMatchSpec { get; set; }

        public ePHRPErrorCodes ErrorCode => mErrorCode;

        public string ErrorMessage => GetErrorMessage();

        public string FastaFilePath { get; set; }

        public string FileVersion => GetVersionForExecutingAssembly();

        public string FileDate { get; protected set; }

        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        /// <summary>
        /// Inspect synopsis file p-value threshold
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float InspectSynopsisFilePValueThreshold { get; set; }

        public string MassCorrectionTagsFilePath { get; set; }

        public string ModificationDefinitionsFilePath { get; set; }

        /// <summary>
        /// Used by clsMODaResultsProcessor and clsMODPlusResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Higher probability are higher confidence results</remarks>
        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        /// <summary>
        /// Used by MSAlign and TopPIC
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float MSAlignAndTopPICSynopsisFilePValueThreshold { get; set; }

        /// <summary>
        /// MSGFDB synopsis file p-value threshold
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower p-values are higher confidence results; renamed to EValue in MS-GF+</remarks>
        [Obsolete("Use MSGFDBSynopsisFileEValueThreshold")]
        public float MSGFDBSynopsisFilePValueThreshold
        {
            get => MSGFDBSynopsisFileEValueThreshold;
            set => MSGFDBSynopsisFileEValueThreshold = value;
        }

        /// <summary>
        /// MSGFDB synopsis file specEValue threshold
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower SpecProb values are higher confidence results; renamed to SpecEValue in MS-GF+</remarks>
        [Obsolete("Use MSGFDBSynopsisFileSpecEValueThreshold")]
        public float MSGFDBSynopsisFileSpecProbThreshold
        {
            get => MSGFDBSynopsisFileSpecEValueThreshold;
            set => MSGFDBSynopsisFileSpecEValueThreshold = value;
        }

        /// <summary>
        /// Used by clsMSGFPlusResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower E-values are higher confidence results</remarks>
        [Obsolete("Use MSGFPlusSynopsisFileEValueThreshold")]
        public float MSGFDBSynopsisFileEValueThreshold
        {
            get => MSGFPlusSynopsisFileEValueThreshold;
            set => MSGFPlusSynopsisFileEValueThreshold = value;
        }

        /// <summary>
        /// clsMSGFPlusResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower SpecEValue values are higher confidence results</remarks>
        [Obsolete("Use MSGFPlusSynopsisFileSpecEValueThreshold")]
        public float MSGFDBSynopsisFileSpecEValueThreshold
        {
            get => MSGFPlusSynopsisFileSpecEValueThreshold;
            set => MSGFPlusSynopsisFileSpecEValueThreshold = value;
        }

        /// <summary>
        /// Used by clsMSGFPlusResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower E-values are higher confidence results</remarks>
        public float MSGFPlusSynopsisFileEValueThreshold { get; set; }

        /// <summary>
        /// clsMSGFPlusResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower SpecEValue values are higher confidence results</remarks>
        public float MSGFPlusSynopsisFileSpecEValueThreshold { get; set; }

        /// <summary>
        /// Typical non-zero value is 17.0027387
        /// </summary>
        /// <returns></returns>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideCTerminusMassChange { get; set; }

        /// <summary>
        ///  typical non-zero value is 1.0078246
        /// </summary>
        /// <returns></returns>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideNTerminusMassChange { get; set; }

        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        public string ProgressStepDescription => mProgressStepDescription;

        // ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public float ProgressPercentComplete => Convert.ToSingle(Math.Round(mProgressPercentComplete, 2));

        /// <summary>
        /// Search tool parameter file path
        /// </summary>
        /// <returns></returns>
        /// <remarks>Used by clsInSpecTResultsProcessor and clsMSGFPlusResultsProcessor (aka SearchEngineParamFileName)</remarks>
        public string SearchToolParameterFilePath { get; set; }

        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        #endregion

        public void AbortProcessingNow()
        {
            AbortProcessing = true;
        }

        public static string AutoDefinePeptideHitResultsFilePath(
            ePeptideHitResultsFileFormatConstants ePeptideHitResultFileFormat,
            string sourceDirectoryPath,
            string baseName)
        {
            if (string.IsNullOrEmpty(baseName))
                return AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath);

            switch (ePeptideHitResultFileFormat)
            {
                case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile:
                    return Path.Combine(sourceDirectoryPath, baseName + SEQUEST_FIRST_HITS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile:
                    return Path.Combine(sourceDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.XTandemXMLFile:
                    return Path.Combine(sourceDirectoryPath, baseName + XTANDEM_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.InspectTXTFile:
                    return Path.Combine(sourceDirectoryPath, baseName + INSPECT_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.MSGFPlusTXTFile:
                    return Path.Combine(sourceDirectoryPath, baseName + MSGFDB_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.MSAlignTXTFile:
                    return Path.Combine(sourceDirectoryPath, baseName + MSALIGN_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.MODaTXTFile:
                    return Path.Combine(sourceDirectoryPath, baseName + MODa_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.MODPlusTXTFile:
                    return Path.Combine(sourceDirectoryPath, baseName + MODPlus_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile:
                    return Path.Combine(sourceDirectoryPath, baseName + MSPathFinder_RESULTS_FILE_SUFFIX);

                case ePeptideHitResultsFileFormatConstants.TopPICTXTFile:
                    return Path.Combine(sourceDirectoryPath, baseName + TopPIC_RESULTS_FILE_SUFFIX);

            }

            return AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath);
        }

        public static string AutoDefinePeptideHitResultsFilePath(string sourceDirectoryPath)
        {
            // Looks for a file ending in _syn.txt, _fht.txt, _xt.xml, or _inspect.txt in directory sourceDirectoryPath
            // Returns the first matching file found

            var matchSpec = string.Empty;

            try
            {
                for (var index = 0; index <= 3; index++)
                {
                    switch (index)
                    {
                        case 0:
                            matchSpec = "*" + SEQUEST_SYNOPSIS_FILE_SUFFIX;
                            break;
                        case 1:
                            matchSpec = "*" + SEQUEST_FIRST_HITS_FILE_SUFFIX;
                            break;
                        case 2:
                            matchSpec = "*" + XTANDEM_RESULTS_FILE_SUFFIX;
                            break;
                        case 3:
                            matchSpec = "*" + INSPECT_RESULTS_FILE_SUFFIX;
                            break;
                    }

                    var sourceDirectory = new DirectoryInfo(sourceDirectoryPath);
                    foreach (var resultsFile in sourceDirectory.GetFiles(matchSpec))
                    {
                        // If we get here, a match was found; return its path
                        return resultsFile.FullName;
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }

            // No match; return empty
            return string.Empty;
        }

        protected bool CacheProteinNamesFromFasta()
        {
            if (string.IsNullOrWhiteSpace(FastaFilePath))
            {
                // Nothing to do
                return true;
            }

            mProteinNameOrder.Clear();
            var reExtractProteinName = new Regex(@"^>([^ ]+)", RegexOptions.Compiled);

            ReportMessage("Caching protein names from the FASTA file");

            try
            {
                var proteinNumber = 0;

                using (var reader = new StreamReader(new FileStream(FastaFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(lineIn))
                            continue;

                        var reMatch = reExtractProteinName.Match(lineIn);

                        if (!reMatch.Success)
                            continue;
                        var proteinName = reMatch.Groups[1].Value;

                        if (mProteinNameOrder.ContainsKey(proteinName))
                            continue;

                        proteinNumber += 1;

                        mProteinNameOrder.Add(proteinName, proteinNumber);
                    }
                }

                ReportMessage(string.Format("Cached {0:N0} proteins", mProteinNameOrder.Count));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error caching protein names from fasta file " + Path.GetFileName(FastaFilePath) + ": " + ex.Message);
                return false;
            }
        }

        protected bool CheckSeqToProteinMapDefined(int uniqueSeqID, string proteinName)
        {
            // Returns True if the sequence to protein map was already defined
            // Returns False if the mapping was not defined (will also update mSeqToProteinMap)

            bool existingMapFound;

            try
            {
                if (proteinName == null)
                    proteinName = string.Empty;

                var key = uniqueSeqID + UNIQUE_SEQ_TO_PROTEIN_MAP_SEP + proteinName;

                if (mSeqToProteinMap.Contains(key))
                {
                    existingMapFound = true;
                }
                else
                {
                    mSeqToProteinMap.Add(key);
                    existingMapFound = false;
                }
            }
            catch (Exception)
            {
                existingMapFound = false;
            }

            return existingMapFound;
        }

        protected int CIntSafe(string value, int defaultValue)
        {
            try
            {
                // Note: Integer.Parse() fails if value contains a decimal point, even if it is "8.000"
                // Thus, we're converting to a double first, and then rounding
                return (int)Math.Round(Convert.ToDouble(value));
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        protected double CDblSafe(string value, double defaultValue)
        {
            try
            {
                return double.Parse(value);
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        protected float CSngSafe(string value, float defaultValue)
        {
            try
            {
                return float.Parse(value);
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        /// <summary>
        /// Validate the input file and output directory paths
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <returns>True if success, False if failure</returns>
        protected bool CleanupFilePaths(ref string inputFilePath, ref string outputDirectoryPath)
        {

            try
            {
                // Make sure inputFilePath points to a valid file
                var inputFile = new FileInfo(inputFilePath);

                if (!inputFile.Exists)
                {
                    SetErrorMessage("Input file not found: " + inputFilePath);
                    if (inputFilePath.Contains(".."))
                    {
                        ReportWarning("Absolute path: " + inputFile.DirectoryName);
                    }
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                if (string.IsNullOrWhiteSpace(outputDirectoryPath))
                {
                    // Define outputDirectoryPath based on inputFilePath
                    outputDirectoryPath = inputFile.DirectoryName;
                }

                // Make sure outputDirectoryPath points to a directory
                DirectoryInfo outputDirectory;

                if (string.IsNullOrWhiteSpace(outputDirectoryPath))
                {
                    outputDirectory = inputFile.Directory;
                }
                else
                {
                    outputDirectory = new DirectoryInfo(outputDirectoryPath);
                }

                if (outputDirectory == null)
                {
                    outputDirectoryPath = ".";
                    return true;
                }

                if (outputDirectory.Exists)
                    return true;

                // outputDirectoryPath points to a non-existent directory; attempt to create it
                try
                {
                    outputDirectory.Create();
                }
                catch (Exception ex2)
                {
                    SetErrorMessage("Invalid output directory: " + outputDirectoryPath, ex2);
                    SetErrorCode(ePHRPErrorCodes.InvalidOutputDirectoryPath);
                    return false;
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error cleaning up the file paths: " + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.FilePathError);
                return false;
            }
        }

        /// <summary>
        /// Collapses a list of strings to a tab-delimited line of text
        /// </summary>
        /// <param name="fields"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected string CollapseList(List<string> fields)
        {
            return string.Join("\t", fields);
        }

        /// <summary>
        /// Examine the extension on filePath to determine the file format
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static ePeptideHitResultsFileFormatConstants DetermineResultsFileFormat(string filePath)
        {

            if (string.IsNullOrWhiteSpace(filePath))
                return ePeptideHitResultsFileFormatConstants.AutoDetermine;

            var extensionLCase = Path.GetExtension(filePath).ToLower();
            var baseFileName = Path.GetFileNameWithoutExtension(filePath).ToLower();

            if (extensionLCase == ".xml")
            {
                return ePeptideHitResultsFileFormatConstants.XTandemXMLFile;
            }

            if (baseFileName.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile;
            }

            if (baseFileName.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.SequestSynopsisFile;
            }

            if (baseFileName.EndsWith(clsInSpecTResultsProcessor.FILENAME_SUFFIX_INSPECT_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.InspectTXTFile;
            }

            if (baseFileName.EndsWith(clsMSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSGFPlusTXTFile;
            }

            if (baseFileName.EndsWith(clsMSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSGFPlusTXTFile;
            }
            if (baseFileName.EndsWith(clsMSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSAlignTXTFile;
            }

            if (baseFileName.EndsWith(clsMODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MODaTXTFile;
            }

            if (baseFileName.EndsWith(clsMODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MODPlusTXTFile;
            }

            if (baseFileName.EndsWith(clsMSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile;
            }

            if (baseFileName.EndsWith(clsTopPICResultsProcessor.FILENAME_SUFFIX_TopPIC_PRSMs_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.TopPICTXTFile;
            }

            if (extensionLCase == ".tsv")
            {
                // Assume this is an MS-GF+ TSV file
                return ePeptideHitResultsFileFormatConstants.MSGFPlusTXTFile;
            }

            // Unknown extension
            return ePeptideHitResultsFileFormatConstants.AutoDetermine;
        }

        protected void CloseSequenceOutputFiles()
        {
            try
            {
                if (mResultToSeqMapFile != null)
                {
                    mResultToSeqMapFile.Close();
                    mResultToSeqMapFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }

            try
            {
                if (mSeqInfoFile != null)
                {
                    mSeqInfoFile.Close();
                    mSeqInfoFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }

            try
            {
                if (mModDetailsFile != null)
                {
                    mModDetailsFile.Close();
                    mModDetailsFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }

            try
            {
                if (mSeqToProteinMapFile != null)
                {
                    mSeqToProteinMapFile.Close();
                    mSeqToProteinMapFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }
        }

        protected void ComputePseudoPeptideLocInProtein(clsSearchResultsBaseClass searchResult)
        {
            // Set these to 1 and 10000 since MSGFDB, Sequest, and Inspect results files do not contain protein sequence information
            // If we find later that the peptide sequence spans the length of the protein, we'll revise .ProteinSeqResidueNumberEnd as needed
            searchResult.ProteinSeqResidueNumberStart = 1;
            searchResult.ProteinSeqResidueNumberEnd = 10000;

            if (searchResult.PeptidePreResidues.Trim().EndsWith(clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST.ToString()))
            {
                // The peptide is at the N-Terminus of the protein
                searchResult.PeptideLocInProteinStart = searchResult.ProteinSeqResidueNumberStart;
                searchResult.PeptideLocInProteinEnd = searchResult.PeptideLocInProteinStart + searchResult.PeptideCleanSequence.Length - 1;

                if (searchResult.PeptidePostResidues.Trim()[0] == clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST)
                {
                    // The peptide spans the entire length of the protein
                    searchResult.ProteinSeqResidueNumberEnd = searchResult.PeptideLocInProteinEnd;
                }
                else
                {
                    if (searchResult.PeptideLocInProteinEnd > searchResult.ProteinSeqResidueNumberEnd)
                    {
                        // The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                        searchResult.ProteinSeqResidueNumberEnd = searchResult.PeptideLocInProteinEnd + 1;
                    }
                }
            }
            else if (searchResult.PeptidePostResidues.Trim().StartsWith(clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST.ToString()))
            {
                // The peptide is at the C-Terminus of the protein
                searchResult.PeptideLocInProteinEnd = searchResult.ProteinSeqResidueNumberEnd;
                searchResult.PeptideLocInProteinStart = searchResult.PeptideLocInProteinEnd - searchResult.PeptideCleanSequence.Length + 1;

                if (searchResult.PeptideLocInProteinStart < searchResult.ProteinSeqResidueNumberStart)
                {
                    // The peptide is more than 10000 characters long; this is highly unlikely
                    searchResult.ProteinSeqResidueNumberEnd = searchResult.ProteinSeqResidueNumberStart + 1 + searchResult.PeptideCleanSequence.Length;
                    searchResult.PeptideLocInProteinEnd = searchResult.ProteinSeqResidueNumberEnd;
                    searchResult.PeptideLocInProteinStart = searchResult.PeptideLocInProteinEnd - searchResult.PeptideCleanSequence.Length + 1;
                }
            }
            else
            {
                searchResult.PeptideLocInProteinStart = searchResult.ProteinSeqResidueNumberStart + 1;
                searchResult.PeptideLocInProteinEnd = searchResult.PeptideLocInProteinStart + searchResult.PeptideCleanSequence.Length - 1;

                if (searchResult.PeptideLocInProteinEnd > searchResult.ProteinSeqResidueNumberEnd)
                {
                    // The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                    searchResult.ProteinSeqResidueNumberEnd = searchResult.PeptideLocInProteinEnd + 1;
                }
            }
        }

        protected virtual string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            return ConstructPepToProteinMapFilePathWork(inputFilePath, outputDirectoryPath, mts);
        }

        private string ConstructPepToProteinMapFilePathWork(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var pepToProteinMapBaseName = Path.GetFileNameWithoutExtension(inputFilePath);

            string pepToProteinMapFileName;
            if (mts)
            {
                pepToProteinMapFileName = pepToProteinMapBaseName + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + "MTS.txt";
            }
            else
            {
                pepToProteinMapFileName = pepToProteinMapBaseName + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + ".txt";
            }

            string pepToProteinMapFilePath;

            if (string.IsNullOrWhiteSpace(outputDirectoryPath))
            {
                var inputFile = new FileInfo(inputFilePath);
                if (string.IsNullOrWhiteSpace(inputFile.DirectoryName))
                {
                    pepToProteinMapFilePath = pepToProteinMapFileName;
                }
                else
                {
                    pepToProteinMapFilePath = Path.Combine(inputFile.DirectoryName, pepToProteinMapFileName);
                }
            }
            else
            {
                pepToProteinMapFilePath = Path.Combine(outputDirectoryPath, pepToProteinMapFileName);
            }

            return pepToProteinMapFilePath;
        }

        protected string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts, List<string> suffixesToFind, int charsToRemove)
        {
            var baseName = Path.GetFileNameWithoutExtension(inputFilePath);
            if (string.IsNullOrEmpty(baseName))
                return string.Empty;

            foreach (var item in suffixesToFind)
            {
                if (!baseName.EndsWith(item, StringComparison.OrdinalIgnoreCase)) continue;

                // baseName matches something like Dataset_msgfplus_syn
                var baseNameTrimmed = baseName.Substring(0, baseName.Length - charsToRemove);
                return ConstructPepToProteinMapFilePathWork(baseNameTrimmed, outputDirectoryPath, mts);
            }

            return ConstructPepToProteinMapFilePathWork(baseName, outputDirectoryPath, mts);
        }

        /// <summary>
        /// Use the PeptideToProteinMapEngine to create the Peptide to Protein map file for the file or files in sourcePHRPDataFiles
        /// </summary>
        /// <param name="sourcePHRPDataFiles"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreatePepToProteinMapFile(List<string> sourcePHRPDataFiles, string mtsPepToProteinMapFilePath)
        {
            var success = false;

            try
            {
                if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because mtsPepToProteinMapFilePath is empty; likely a programming bug");
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                if (string.IsNullOrEmpty(FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File Path is not defined");
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                if (!File.Exists(FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File was not found: " + FastaFilePath);
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                // Verify that the fasta file is not a DNA-sequence based fasta file
                success = ValidateProteinFastaFile(FastaFilePath);
                if (!success)
                {
                    return false;
                }

                Console.WriteLine();
                Console.WriteLine();
                UpdateProgress("Creating Peptide to Protein Map file", PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE);

                // Initialize items
                var mtsPepToProteinMapFile = new FileInfo(mtsPepToProteinMapFilePath);
                string outputDirectoryPath;

                if (string.IsNullOrWhiteSpace(mtsPepToProteinMapFile.DirectoryName))
                    outputDirectoryPath = string.Empty;
                else
                    outputDirectoryPath = mtsPepToProteinMapFile.DirectoryName;

                var peptideToProteinMapResults = new SortedSet<string>();

                var peptideToProteinMapper = new clsPeptideToProteinMapEngine
                {
                    DeleteTempFiles = true,
                    IgnoreILDifferences = false,
                    InspectParameterFilePath = string.Empty,
                    LogMessagesToFile = false,
                    MatchPeptidePrefixAndSuffixToProtein = false,
                    OutputProteinSequence = false,
                    PeptideInputFileFormat = clsPeptideToProteinMapEngine.ePeptideInputFileFormatConstants.PHRPFile,
                    PeptideFileSkipFirstLine = false,
                    ProteinDataRemoveSymbolCharacters = true,
                    ProteinInputFilePath = FastaFilePath,
                    SaveProteinToPeptideMappingFile = true,
                    SearchAllProteinsForPeptideSequence = true,
                    SearchAllProteinsSkipCoverageComputationSteps = true
                };

                RegisterEvents(peptideToProteinMapper);

                // Handle progress updates using PeptideToProteinMapper_ProgressChanged instead of OnProgressUpdate
                peptideToProteinMapper.ProgressUpdate -= OnProgressUpdate;
                peptideToProteinMapper.ProgressUpdate += PeptideToProteinMapper_ProgressChanged;
                peptideToProteinMapper.SkipConsoleWriteIfNoProgressListener = true;

                using (var writer = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    foreach (var inputFilePath in sourcePHRPDataFiles)
                    {
                        var resultsFilePath = Path.GetFileNameWithoutExtension(inputFilePath) +
                            clsPeptideToProteinMapEngine.FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING;

                        resultsFilePath = Path.Combine(outputDirectoryPath, resultsFilePath);

                        // Make sure the results file doesn't already exist
                        DeleteFileIgnoreErrors(resultsFilePath);

                        peptideToProteinMapper.ProgressUpdate += PeptideToProteinMapper_ProgressChanged;
                        mNextPeptideToProteinMapperLevel = 25;

                        success = peptideToProteinMapper.ProcessFile(inputFilePath, outputDirectoryPath, string.Empty, true);

                        peptideToProteinMapper.ProgressUpdate -= PeptideToProteinMapper_ProgressChanged;

                        if (success)
                        {
                            if (!File.Exists(resultsFilePath))
                            {
                                SetErrorMessage("Peptide to protein mapping file was not created for " + inputFilePath);
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                success = false;
                                break;
                            }
                            success = ValidatePeptideToProteinMapResults(resultsFilePath, IgnorePeptideToProteinMapperErrors);
                        }
                        else
                        {
                            if (string.IsNullOrWhiteSpace(peptideToProteinMapper.GetErrorMessage()) && peptideToProteinMapper.StatusMessage.IndexOf("error", StringComparison.OrdinalIgnoreCase) >= 0)
                            {
                                SetErrorMessage("Error running clsPeptideToProteinMapEngine: " + peptideToProteinMapper.StatusMessage);
                            }
                            else
                            {
                                if (peptideToProteinMapper.StatusMessage.Length > 0)
                                {
                                    SetErrorMessage("clsPeptideToProteinMapEngine status: " + peptideToProteinMapper.StatusMessage);
                                }
                                SetErrorMessage("Error running clsPeptideToProteinMapEngine: " + peptideToProteinMapper.GetErrorMessage());
                            }

                            if (IgnorePeptideToProteinMapperErrors)
                            {
                                ReportWarning("Ignoring protein mapping error since 'IgnorePeptideToProteinMapperErrors' = True");

                                if (File.Exists(resultsFilePath))
                                {
                                    success = ValidatePeptideToProteinMapResults(resultsFilePath, IgnorePeptideToProteinMapperErrors);
                                }
                                else
                                {
                                    mErrorMessage = string.Empty;
                                    mErrorCode = ePHRPErrorCodes.NoError;
                                    success = true;
                                }
                            }
                            else
                            {
                                SetErrorMessage("Error in CreatePepToProteinMapFile");
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                            }
                        }

                        if (!File.Exists(resultsFilePath))
                        {
                            continue;
                        }

                        // Read the newly created file and append new entries to mtsPepToProteinMapFilePath
                        using (var reader = new StreamReader(new FileStream(resultsFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                        {
                            while (!reader.EndOfStream)
                            {
                                var lineIn = reader.ReadLine();

                                if (string.IsNullOrWhiteSpace(lineIn)) continue;

                                var splitLine = lineIn.Split(new[] { '\t' }, 2);
                                if (splitLine.Length < 2) continue;

                                var peptideAndProteinKey = splitLine[0] + "_" + splitLine[1];

                                if (!peptideToProteinMapResults.Contains(peptideAndProteinKey))
                                {
                                    peptideToProteinMapResults.Add(peptideAndProteinKey);
                                    writer.WriteLine(lineIn);
                                }
                            }
                        }

                        // Delete the interim results file
                        DeleteFileIgnoreErrors(resultsFilePath);
                    }
                }

                peptideToProteinMapper.CloseLogFileNow();
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreatePepToProteinMapFile:" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
            }

            return success;
        }

        /// <summary>
        /// Create the protein mod details file for the specified PHRP data file
        /// </summary>
        /// <param name="phrpDataFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool CreateProteinModDetailsFile(string phrpDataFilePath, string outputDirectoryPath)
        {
            var success = false;

            try
            {
                var inputFile = new FileInfo(phrpDataFilePath);

                var sourcePHRPDataFiles = new List<string> {
                    inputFile.FullName
                };

                var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, outputDirectoryPath, mts: true);

                success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);

                if (success)
                {
                    success = CreateProteinModDetailsFile(phrpDataFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Unknown);
                }
            }
            catch (Exception ex)
            {
                ReportWarning("Error in CreateProteinModDetailsFile: " + ex.Message);
            }

            return success;
        }

        public bool CreateProteinModDetailsFile(
            string phrpDataFilePath,
            string outputDirectoryPath,
            string mtsPepToProteinMapFilePath,
            clsPHRPReader.ePeptideHitResultType ePHRPResultType)
        {
            try
            {
                Console.WriteLine();

                var progressAtStart = mProgressPercentComplete;
                UpdateProgress("Creating the Protein Mod Details file", PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE);

                // Confirm that the PHRP data file exists
                var phrpDataFile = new FileInfo(phrpDataFilePath);
                if (!phrpDataFile.Exists)
                {
                    SetErrorMessage("PHRP data file not found in CreateProteinModDetailsFile: " + phrpDataFilePath);
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                // Confirm that the _PepToProtMapMTS.txt file exists
                if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath))
                {
                    SetErrorMessage("Cannot create the ProteinMods file because mtsPepToProteinMapFilePath is empty");
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                // Initialize pepToProteinMapping
                var pepToProteinMapping = new List<udtPepToProteinMappingType>();

                // Read the _PepToProtMapMTS file
                var success = LoadPeptideToProteinMapInfo(mtsPepToProteinMapFilePath, pepToProteinMapping, out _);
                if (!success)
                {
                    return false;
                }

                // Assure that pepToProteinMapping is sorted on peptide
                if (pepToProteinMapping.Count > 1)
                {
                    pepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                var proteinModsFilePath = ReplaceFilenameSuffix(phrpDataFile, FILENAME_SUFFIX_PROTEIN_MODS);
                if (!string.IsNullOrEmpty(outputDirectoryPath))
                {
                    var proteinModsFile = Path.GetFileName(proteinModsFilePath);
                    if (string.IsNullOrEmpty(proteinModsFile))
                    {
                        proteinModsFile = Path.GetFileNameWithoutExtension(phrpDataFile.Name) + FILENAME_SUFFIX_PROTEIN_MODS;
                    }
                    proteinModsFilePath = Path.Combine(outputDirectoryPath, proteinModsFile);
                }

                var psmCount = 0;
                var psmCountSkippedSinceReversedOrScrambledProtein = 0;

                // Create a ProteinMods file parallel to the PHRP file
                using (var writer = new StreamWriter(new FileStream(proteinModsFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    // Write the header line
                    writer.WriteLine(COLUMN_NAME_RESULTID + "\t" +
                                     COLUMN_NAME_PEPTIDE + "\t" +
                                     COLUMN_NAME_UNIQUE_SEQ_ID + "\t" +
                                     COLUMN_NAME_PROTEIN_NAME + "\t" +
                                     COLUMN_NAME_RESIDUE + "\t" +
                                     COLUMN_NAME_PROTEIN_RESIDUE_NUMBER + "\t" +
                                     COLUMN_NAME_RESIDUE_MOD_NAME + "\t" +
                                     COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER + "\t" +
                                     COLUMN_NAME_MSGF_SPECPROB);

                    var loadMSGFResults = ePHRPResultType != clsPHRPReader.ePeptideHitResultType.MSGFPlus;

                    // Update the Mass Calculator to use the one tracked by this class
                    // (since this class's calculator knows about custom amino acids and custom charge carriers)
                    var startupOptions = new clsPHRPStartupOptions
                    {
                        LoadModsAndSeqInfo = true,
                        LoadMSGFResults = loadMSGFResults,
                        LoadScanStatsData = false,
                        PeptideMassCalculator = mPeptideSeqMassCalculator
                    };

                    using (var reader = new clsPHRPReader(phrpDataFilePath, ePHRPResultType, startupOptions))
                    {
                        reader.EchoMessagesToConsole = false;
                        reader.SkipDuplicatePSMs = true;

                        foreach (var errorMessage in reader.ErrorMessages)
                        {
                            OnErrorEvent(errorMessage);
                        }
                        RegisterEvents(reader);

                        foreach (var warningMessage in reader.WarningMessages)
                        {
                            var msg = warningMessage;
                            if (warningMessage.StartsWith("MSGF file not found", StringComparison.OrdinalIgnoreCase))
                            {
                                msg = "MSGF file not found; column " + COLUMN_NAME_MSGF_SPECPROB + " will not have any data";
                            }
                            ReportWarning(msg);
                        }

                        reader.ClearErrors();
                        reader.ClearWarnings();

                        var peptidesNotFoundInPepToProtMapping = 0;
                        while (reader.MoveNext())
                        {
                            // Use binary search to find this peptide in pepToProteinMapping
                            var pepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(pepToProteinMapping, reader.CurrentPSM.Peptide);

                            if (pepToProteinMapIndex >= 0)
                            {
                                do
                                {
                                    psmCount += 1;

                                    var skipProtein = false;
                                    if (!ProteinModsFileIncludesReversedProteins)
                                    {
                                        skipProtein = IsReversedProtein(pepToProteinMapping[pepToProteinMapIndex].Protein);
                                        if (skipProtein)
                                        {
                                            psmCountSkippedSinceReversedOrScrambledProtein += 1;
                                        }
                                    }

                                    if (!skipProtein)
                                    {
                                        WriteModDetailsEntry(reader,
                                                             writer,
                                                             pepToProteinMapping,
                                                             pepToProteinMapIndex,
                                                             ref psmCountSkippedSinceReversedOrScrambledProtein);
                                    }

                                    pepToProteinMapIndex += 1;
                                } while (pepToProteinMapIndex < pepToProteinMapping.Count &&
                                         reader.CurrentPSM.Peptide == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                            }
                            else
                            {
                                peptidesNotFoundInPepToProtMapping++;
                                ShowPeriodicWarning(peptidesNotFoundInPepToProtMapping, 10,
                                                    "Peptide not found in pepToProteinMapping: " + reader.CurrentPSM.Peptide);
                            }

                            UpdateProgress(progressAtStart + reader.PercentComplete * (100 - progressAtStart) / 100);
                        }

                        if (psmCount > 0)
                        {
                            if (psmCountSkippedSinceReversedOrScrambledProtein == psmCount)
                            {
                                Console.WriteLine();
                                ReportWarning("All PSMs map to reversed or scrambled proteins; the _ProteinMods.txt file is empty");
                            }
                            else if (psmCountSkippedSinceReversedOrScrambledProtein > 0)
                            {
                                Console.WriteLine();
                                Console.WriteLine("Note: skipped {0:N0} / {1:N0} PSMs that map to reversed or scrambled proteins " +
                                                  "while creating the _ProteinMods.txt file",
                                                  psmCountSkippedSinceReversedOrScrambledProtein, psmCount);
                            }
                        }

                        if (peptidesNotFoundInPepToProtMapping > 10)
                        {
                            Console.WriteLine("Note: {0:N0} peptides were not found in pepToProteinMapping", peptidesNotFoundInPepToProtMapping);
                        }
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateProteinModDetailsFile:" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }

        }

        protected void DeleteFileIgnoreErrors(string filePath)
        {
            try
            {
                if (File.Exists(filePath))
                {
                    Thread.Sleep(200);
                    File.Delete(filePath);
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        protected void ExpandListIfRequired<T>(List<T> items, int countToAdd, int largeListThreshold = 1000000)
        {
            if (items.Count > largeListThreshold && items.Count + countToAdd > items.Capacity)
            {
                // .NET by default will double the size of the list to accomodate these new items
                // Instead, expand the list by 20% of the current size
                items.Capacity = items.Capacity + Convert.ToInt32(items.Count / 5);
            }
        }

        private readonly IComparer<udtPepToProteinMappingType> peptideSearchComparer = new PepToProteinMappingPeptideSearchComparer();

        protected int FindFirstMatchInPepToProteinMapping(List<udtPepToProteinMappingType> pepToProteinMapping, string peptideToFind)
        {
            // Use binary search to find this peptide in pepToProteinMapping
            var udtItemToFind = new udtPepToProteinMappingType
            {
                Peptide = peptideToFind
            };

            var pepToProteinMapIndex = pepToProteinMapping.BinarySearch(udtItemToFind, peptideSearchComparer);

            if (pepToProteinMapIndex > 0)
            {
                // Step Backward until the first match is found
                while (pepToProteinMapIndex > 0 && pepToProteinMapping[pepToProteinMapIndex - 1].Peptide == peptideToFind)
                {
                    pepToProteinMapIndex -= 1;
                }
            }

            return pepToProteinMapIndex;
        }

        protected string GetCleanSequence(string sequenceWithMods)
        {
            return GetCleanSequence(sequenceWithMods, out _, out _);
        }

        protected string GetCleanSequence(string sequenceWithMods, out string prefix, out string suffix)
        {
            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out prefix, out suffix))
            {
                // Remove any non-letter characters
                primarySequence = mReplaceSymbols.Replace(primarySequence, string.Empty);
            }
            else
            {
                // Sequence does not have prefix or suffix letters; use sequenceWithMods
                primarySequence = mReplaceSymbols.Replace(sequenceWithMods, string.Empty);
            }

            return primarySequence;
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to String.Empty
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out string value)
        {
            return GetColumnValue(splitLine, columnIndex, out value, string.Empty);
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to 0
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out int value)
        {
            return GetColumnValue(splitLine, columnIndex, out value, 0);
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out string value, string valueIfMissing)
        {
            if (columnIndex >= 0 && columnIndex < splitLine.Length)
            {
                value = string.Copy(splitLine[columnIndex]);
                return true;
            }

            value = string.Copy(valueIfMissing);
            return false;
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out int value, int valueIfMissing)
        {
            if (GetColumnValue(splitLine, columnIndex, out var valueText, valueIfMissing.ToString()))
            {
                if (int.TryParse(valueText, out value))
                {
                    return true;
                }

                value = valueIfMissing;
                return false;
            }

            value = valueIfMissing;
            return false;
        }

        /// <summary>
        /// Get the error message, or an empty string if no error
        /// </summary>
        /// <returns></returns>
        protected string GetErrorMessage()
        {
            string message;

            switch (ErrorCode)
            {
                case ePHRPErrorCodes.NoError:
                    message = string.Empty;
                    break;
                case ePHRPErrorCodes.InvalidInputFilePath:
                    message = "Invalid input file path";
                    break;
                case ePHRPErrorCodes.InvalidOutputDirectoryPath:
                    message = "Invalid output directory path";
                    break;
                case ePHRPErrorCodes.ParameterFileNotFound:
                    message = "Parameter file not found";
                    break;
                case ePHRPErrorCodes.MassCorrectionTagsFileNotFound:
                    message = "Mass correction tags file not found";
                    break;
                case ePHRPErrorCodes.ModificationDefinitionFileNotFound:
                    message = "Modification definition file not found";

                    break;
                case ePHRPErrorCodes.ErrorReadingInputFile:
                    message = "Error reading input file";
                    break;
                case ePHRPErrorCodes.ErrorCreatingOutputFiles:
                    message = "Error creating output files";
                    break;
                case ePHRPErrorCodes.ErrorReadingParameterFile:
                    message = "Invalid parameter file";
                    break;
                case ePHRPErrorCodes.ErrorReadingMassCorrectionTagsFile:
                    message = "Error reading mass correction tags file";
                    break;
                case ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile:
                    message = "Error reading modification definitions file";

                    break;
                case ePHRPErrorCodes.FilePathError:
                    message = "General file path error";
                    break;
                case ePHRPErrorCodes.UnspecifiedError:
                    message = "Unspecified error";
                    break;
                default:
                    // This shouldn't happen
                    message = "Unknown error state";
                    break;
            }

            if (mErrorMessage.Length > 0)
            {
                if (message.Length > 0)
                {
                    message += "; ";
                }
                message += mErrorMessage;
            }

            return message;
        }

        private string GetVersionForExecutingAssembly()
        {
            string version;

            try
            {
                version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            }
            catch (Exception)
            {
                version = "??.??.??.??";
            }

            return version;
        }

        private void InitializeLocalVariables()
        {
            mErrorCode = ePHRPErrorCodes.NoError;
            mErrorMessage = string.Empty;
            WarnMissingParameterFileSection = false;

            CreateModificationSummaryFile = true;
            CreateProteinModsFile = false;
            FastaFilePath = string.Empty;
            IgnorePeptideToProteinMapperErrors = false;
            ProteinModsFileIncludesReversedProteins = false;
            UseExistingMTSPepToProteinMapFile = false;

            CreateInspectFirstHitsFile = true;
            CreateInspectSynopsisFile = true;

            MassCorrectionTagsFilePath = string.Empty;
            ModificationDefinitionsFilePath = string.Empty;
            SearchToolParameterFilePath = string.Empty;

            InspectSynopsisFilePValueThreshold = clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MODaMODPlusSynopsisFileProbabilityThreshold = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            MSAlignAndTopPICSynopsisFilePValueThreshold = clsMSAlignResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MSGFPlusSynopsisFileEValueThreshold = clsMSGFPlusResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MSGFPlusSynopsisFileSpecEValueThreshold = clsMSGFPlusResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            EnzymeMatchSpec = clsPeptideCleavageStateCalculator.GetDefaultEnzymeMatchSpec();

            PeptideNTerminusMassChange = clsPeptideMassCalculator.DEFAULT_N_TERMINUS_MASS_CHANGE;
            PeptideCTerminusMassChange = clsPeptideMassCalculator.DEFAULT_C_TERMINUS_MASS_CHANGE;
        }

        /// <summary>
        /// Initializes the StreamWriter objects using baseOutputFilePath as a base name and replacing the suffix with the default suffix names
        /// </summary>
        /// <param name="baseOutputFilePath"></param>
        /// <returns>True if success; does not catch errors; they will be thrown to the calling function if they occur</returns>
        protected bool InitializeSequenceOutputFiles(string baseOutputFilePath)
        {
            var outputFileInfo = new FileInfo(baseOutputFilePath);

            // Initialize the file paths based on baseOutputFilePath
            var resultToSeqMapFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_RESULT_TO_SEQ_MAP);
            var seqInfoFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_SEQ_INFO);
            var modDetailsFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_MOD_DETAILS);
            var seqToProteinMapFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP);

            // Clear the unique sequences container
            mUniqueSequences.Clear();

            // Clear the sequence to protein map
            mSeqToProteinMap.Clear();

            var resultToSeqMapHeaders = new List<string>
            {
                "Result_ID",
                COLUMN_NAME_UNIQUE_SEQ_ID
            };

            var seqInfoHeaders = new List<string>
            {
                COLUMN_NAME_UNIQUE_SEQ_ID,
                "Mod_Count",
                "Mod_Description",
                "Monoisotopic_Mass"
            };

            var modDetailsHeaders = new List<string>
            {
                COLUMN_NAME_UNIQUE_SEQ_ID,
                "Mass_Correction_Tag",
                "Position"
            };

            var seqToProteinMapHeaders = new List<string>
            {
                COLUMN_NAME_UNIQUE_SEQ_ID,
                "Cleavage_State",
                "Terminus_State",
                COLUMN_NAME_PROTEIN_NAME,
                "Protein_Expectation_Value_Log(e)",
                "Protein_Intensity_Log(I)"
            };


            // Initialize the ResultToSeqMap file
            mResultToSeqMapFile = new StreamWriter(resultToSeqMapFilePath);
            mResultToSeqMapFile.WriteLine(CollapseList(resultToSeqMapHeaders));

            // Initialize the SeqInfo file
            mSeqInfoFile = new StreamWriter(seqInfoFilePath, false);
            mSeqInfoFile.WriteLine(CollapseList(seqInfoHeaders));

            // Initialize the ModDetails file
            mModDetailsFile = new StreamWriter(modDetailsFilePath);
            mModDetailsFile.WriteLine(CollapseList(modDetailsHeaders));

            // Initialize the SeqToProtein map file
            mSeqToProteinMapFile = new StreamWriter(seqToProteinMapFilePath, false);
            mSeqToProteinMapFile.WriteLine(CollapseList(seqToProteinMapHeaders));

            return true;
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

        protected bool IsReversedProtein(string proteinName)
        {
            if (proteinName.StartsWith("reversed_", StringComparison.OrdinalIgnoreCase))
            {
                // Used in DMS-generated protein collections
                return true;
            }

            if (proteinName.StartsWith("REV_", StringComparison.OrdinalIgnoreCase))
            {
                // Used by MSGFDB
                return true;
            }

            if (proteinName.StartsWith("scrambled_", StringComparison.OrdinalIgnoreCase))
            {
                // Used in DMS-generated protein collections
                return true;
            }
            if (proteinName.StartsWith("xxx_", StringComparison.OrdinalIgnoreCase))
            {
                // Used by MS-GF+ and MSFragger
                return true;
            }

            if (proteinName.StartsWith("xxx.", StringComparison.OrdinalIgnoreCase))
            {
                // Used by Inspect
                return true;
            }

            if (proteinName.EndsWith(":reversed", StringComparison.OrdinalIgnoreCase))
            {
                // Used by X!Tandem
                return true;
            }

            return false;
        }

        protected virtual bool LoadParameterFileSettings(string parameterFilePath)
        {
            const string OPTIONS_SECTION = "PeptideHitResultsProcessorOptions";

            var settingsFile = new PRISM.XmlSettingsFileAccessor();

            try
            {
                if (string.IsNullOrWhiteSpace(parameterFilePath))
                {
                    // No parameter file specified; nothing to load
                    return true;
                }

                if (!File.Exists(parameterFilePath))
                {
                    // See if parameterFilePath points to a file in the same directory as the application
                    var appDirPath = PRISM.FileProcessor.ProcessFilesOrDirectoriesBase.GetAppDirectoryPath();
                    if (string.IsNullOrWhiteSpace(appDirPath))
                    {
                        SetErrorCode(ePHRPErrorCodes.ParameterFileNotFound);
                        return false;
                    }

                    parameterFilePath = Path.Combine(appDirPath, Path.GetFileName(parameterFilePath));
                    if (!File.Exists(parameterFilePath))
                    {
                        SetErrorCode(ePHRPErrorCodes.ParameterFileNotFound);
                        return false;
                    }
                }

                if (settingsFile.LoadSettings(parameterFilePath))
                {
                    if (!settingsFile.SectionPresent(OPTIONS_SECTION))
                    {
                        // Section OPTIONS_SECTION was not found in the parameter file; warn the user if mWarnMissingParameterFileSection = True
                        if (WarnMissingParameterFileSection)
                        {
                            SetErrorMessage("The node '<section name=\"" + OPTIONS_SECTION + "\"> was not found in the parameter file: " + parameterFilePath);
                            SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile);
                            return false;
                        }
                        return true;
                    }

                    MassCorrectionTagsFilePath = settingsFile.GetParam(OPTIONS_SECTION, "MassCorrectionTagsFilePath", MassCorrectionTagsFilePath);
                    ModificationDefinitionsFilePath = settingsFile.GetParam(OPTIONS_SECTION, "ModificationDefinitionsFilePath", ModificationDefinitionsFilePath);
                    SearchToolParameterFilePath = settingsFile.GetParam(OPTIONS_SECTION, "SearchToolParameterFilePath", SearchToolParameterFilePath);

                    CreateModificationSummaryFile = settingsFile.GetParam(OPTIONS_SECTION, "CreateModificationSummaryFile", CreateModificationSummaryFile);

                    CreateProteinModsFile = settingsFile.GetParam(OPTIONS_SECTION, "CreateProteinModsFile", CreateProteinModsFile);
                    FastaFilePath = settingsFile.GetParam(OPTIONS_SECTION, "FastaFilePath", FastaFilePath);
                    ProteinModsFileIncludesReversedProteins = settingsFile.GetParam(OPTIONS_SECTION, "ProteinModsFileIncludesReversedProteins", ProteinModsFileIncludesReversedProteins);
                    UseExistingMTSPepToProteinMapFile = settingsFile.GetParam(OPTIONS_SECTION, "UseExistingMTSPepToProteinMapFile", UseExistingMTSPepToProteinMapFile);

                    var leftResidueRegEx = string.Copy(EnzymeMatchSpec.LeftResidueRegEx);
                    var rightResidueRegEx = string.Copy(EnzymeMatchSpec.RightResidueRegEx);

                    leftResidueRegEx = settingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecLeftResidue", leftResidueRegEx, out var valueNotPresent);
                    if (!valueNotPresent)
                    {
                        rightResidueRegEx = settingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecRightResidue", rightResidueRegEx, out valueNotPresent);

                        if (!valueNotPresent)
                        {
                            EnzymeMatchSpec = new clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType(leftResidueRegEx, rightResidueRegEx);
                        }
                    }

                    PeptideNTerminusMassChange = settingsFile.GetParam(OPTIONS_SECTION, "PeptideNTerminusMassChange", PeptideNTerminusMassChange);
                    PeptideCTerminusMassChange = settingsFile.GetParam(OPTIONS_SECTION, "PeptideCTerminusMassChange", PeptideCTerminusMassChange);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadParameterFileSettings:" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Load the PeptideToProteinMap information
        /// </summary>
        /// <param name="pepToProteinMapFilePath">File to read</param>
        /// <param name="pepToProteinMapping">Output parameter: peptide to protein mapping (calling function must pre-initialize the list)</param>
        /// <param name="headerLine">Output parameter: Header line text</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool LoadPeptideToProteinMapInfo(
            string pepToProteinMapFilePath,
            List<udtPepToProteinMappingType> pepToProteinMapping,
            out string headerLine)
        {
            headerLine = string.Empty;

            bool success;

            try
            {
                // Initialize the output parameters
                pepToProteinMapping.Clear();
                headerLine = string.Empty;

                if (string.IsNullOrWhiteSpace(pepToProteinMapFilePath))
                {
                    SetErrorMessage("Warning: PepToProteinMap file is not defined");
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    SetErrorMessage("Warning: PepToProteinMap file does not exist: " + pepToProteinMapFilePath);
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                // Open proteinToPeptideMappingFilePath for reading
                using (var reader = new StreamReader(new FileStream(pepToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    var linesRead = 0;
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(lineIn))
                            continue;

                        var dataLine = lineIn.Trim();
                        if (dataLine.Length == 0)
                            continue;

                        linesRead++;

                        // Split the line on tabs
                        var splitLine = dataLine.TrimEnd().Split('\t');

                        if (splitLine.Length >= 4)
                        {
                            if (linesRead == 1 && !int.TryParse(splitLine[2], out _))
                            {
                                // Header line; cache it
                                headerLine = string.Copy(dataLine);
                            }
                            else
                            {
                                var pepToProteinMappingEntry = new udtPepToProteinMappingType
                                {
                                    Peptide = string.Copy(splitLine[0]),
                                    Protein = string.Copy(splitLine[1])
                                };
                                int.TryParse(splitLine[2], out pepToProteinMappingEntry.ResidueStart);
                                int.TryParse(splitLine[3], out pepToProteinMappingEntry.ResidueEnd);

                                ExpandListIfRequired(pepToProteinMapping, 1);

                                pepToProteinMapping.Add(pepToProteinMappingEntry);
                            }
                        }
                    }
                }

                success = true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the Peptide to Protein Map File (" + Path.GetFileName(pepToProteinMapFilePath) + "): " + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                success = false;
            }

            return success;
        }

        protected string MassErrorToString(double massErrorDa)
        {
            if (Math.Abs(massErrorDa) < 0.000001)
                return "0";

            return Math.Abs(massErrorDa) < 0.0001 ?
                PRISM.StringUtilities.DblToString(massErrorDa, 6, 0.0000001) :
                PRISM.StringUtilities.DblToString(massErrorDa, 5, 0.000001);
        }

        protected void OperationComplete()
        {
            ProgressComplete?.Invoke();
        }

        public bool ProcessFile(string inputFilePath, string outputDirectoryPath)
        {
            return ProcessFile(inputFilePath, outputDirectoryPath, string.Empty);
        }

        public abstract bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath);

        /// <summary>
        /// Appends newSuffix to the base name of the original file, then returns a full path to the file using the directory associated with originalFilePath
        /// Note that newSuffix may contain a file extension though it does not have to
        /// If newSuffix does not contain an extension, the path returned will end in the same extension as originalFilePath
        /// </summary>
        /// <param name="originalFile"></param>
        /// <param name="newSuffix"></param>
        /// <returns></returns>
        protected string ReplaceFilenameSuffix(FileInfo originalFile, string newSuffix)
        {
            // Keep track of the original extension on originalFilePath
            var originalExtension = originalFile.Extension;

            // Make sure newSuffix is not nothing
            if (newSuffix == null)
                newSuffix = string.Empty;

            // Obtain the filename, without its extension
            var newFileBaseName = Path.GetFileNameWithoutExtension(originalFile.Name);

            // Append newSuffix to newFileBaseName
            string newFileName;
            if (Path.HasExtension(newSuffix))
            {
                newFileName = newFileBaseName + newSuffix;
            }
            else
            {
                newFileName = newFileBaseName + newSuffix + originalExtension;
            }

            if (string.IsNullOrWhiteSpace(originalFile.DirectoryName))
                return newFileName;

            var newFilePath = Path.Combine(originalFile.DirectoryName, newFileName);

            return newFilePath;
        }

        protected void ReportError(string errMsg, bool throwException = false, Exception ex = null)
        {
            SetErrorMessage(errMsg);

            if (throwException)
            {
                if (ex == null)
                {
                    throw new Exception(errMsg);
                }

                throw new Exception(errMsg, ex);
            }
        }

        protected void ReportMessage(string message)
        {
            OnStatusEvent(message);
        }

        protected void ReportWarning(string message)
        {
            OnWarningEvent(message);
        }

        public bool ResetMassCorrectionTagsAndModificationDefinitions()
        {
            var fileNotFound = false;

            // Note: If mMassCorrectionTagsFilePath is blank, the mass correction tags will be reset to the defaults and success will be True
            var success = mPeptideMods.ReadMassCorrectionTagsFile(MassCorrectionTagsFilePath, ref fileNotFound);
            if (!success)
            {
                if (fileNotFound)
                {
                    SetErrorCode(ePHRPErrorCodes.MassCorrectionTagsFileNotFound);
                }
                else
                {
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingMassCorrectionTagsFile);
                }
            }

            // Note: If mModificationDefinitionsFilePath is blank, the modifications will be cleared and success will be True
            success = mPeptideMods.ReadModificationDefinitionsFile(ModificationDefinitionsFilePath, ref fileNotFound);
            if (!success)
            {
                if (fileNotFound)
                {
                    SetErrorCode(ePHRPErrorCodes.ModificationDefinitionFileNotFound);
                    ReportWarning("File not found: " + ModificationDefinitionsFilePath);
                }
                else
                {
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                }
            }

            return success;
        }

        protected void ResetProgress()
        {
            ResetProgress(string.Empty);
        }

        protected void ResetProgress(string progressStepDescription, bool echoToConsole = false)
        {
            mProgressStepDescription = string.Copy(progressStepDescription);
            mProgressPercentComplete = 0;
            ProgressReset?.Invoke();

            if (echoToConsole)
            {
                Console.WriteLine();
                Console.WriteLine();
                Console.WriteLine(ProgressStepDescription);
            }
        }

        protected void SaveModificationSummaryFile(string modificationSummaryFilePath)
        {
            using (var writer = new StreamWriter(modificationSummaryFilePath, false))
            {
                // Write the header line
                var headerNames = new List<string>
                {
                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol,
                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass,
                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Target_Residues,
                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Type,
                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Mass_Correction_Tag,
                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Occurrence_Count
                };

                writer.WriteLine(CollapseList(headerNames));

                for (var index = 0; index <= mPeptideMods.ModificationCount - 1; index++)
                {
                    var modInfo = mPeptideMods.GetModificationByIndex(index);
                    var searchResult = modInfo;
                    if (searchResult.OccurrenceCount <= 0 && searchResult.UnknownModAutoDefined)
                        continue;

                    var data = new List<string>
                    {
                        searchResult.ModificationSymbol.ToString(),
                        searchResult.ModificationMass.ToString(CultureInfo.InvariantCulture),
                        searchResult.TargetResidues,
                        clsModificationDefinition.ModificationTypeToModificationSymbol(searchResult.ModificationType).ToString(),
                        searchResult.MassCorrectionTag,
                        searchResult.OccurrenceCount.ToString()
                    };
                    writer.WriteLine(CollapseList(data));
                }
            }
        }

        protected void SaveResultsFileEntrySeqInfo(clsSearchResultsBaseClass searchResult, bool updateResultToSeqMapFile)
        {
            // Note: Be sure to call Me.InitializeOutputFiles before calling this function
            // updateResultToSeqMapFile should be set to True only for the first protein of each peptide in each group

            // This ID is assigned using a SortedSet containing mPeptideCleanSequence and mPeptideModDescription
            var uniqueSeqID = mUniqueSequences.GetNextUniqueSequenceID(
                searchResult.PeptideCleanSequence,
                searchResult.PeptideModDescription,
                out var existingSequenceFound);

            if (updateResultToSeqMapFile)
            {
                // Write a new entry to the ResultToSeqMap file
                var seqMapData = new List<string>
                {
                    searchResult.ResultID.ToString(),
                    uniqueSeqID.ToString()
                };
                mResultToSeqMapFile.WriteLine(CollapseList(seqMapData));

                // Only write this entry to the SeqInfo and ModDetails files if existingSequenceFound is False

                if (!existingSequenceFound)
                {
                    // Write a new entry to the SeqInfo file
                    var seqInfoData = new List<string>
                    {
                        uniqueSeqID.ToString(),
                        searchResult.SearchResultModificationCount.ToString(),
                        searchResult.PeptideModDescription,
                        PRISM.StringUtilities.DblToString(searchResult.PeptideMonoisotopicMass, 5, 0.000001)
                    };
                    mSeqInfoFile.WriteLine(CollapseList(seqInfoData));

                    if (searchResult.SearchResultModificationCount > 0)
                    {
                        var udtModNameAndResidueLoc = new udtModNameAndResidueLocType[searchResult.SearchResultModificationCount];
                        var pointerArray = new int[searchResult.SearchResultModificationCount];

                        if (searchResult.SearchResultModificationCount == 1)
                        {
                            pointerArray[0] = 0;
                        }
                        else
                        {
                            // Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                            for (var index = 0; index <= searchResult.SearchResultModificationCount - 1; index++)
                            {
                                var resultModDetails = searchResult.GetSearchResultModDetailsByIndex(index);
                                udtModNameAndResidueLoc[index].ResidueLocInPeptide = resultModDetails.ResidueLocInPeptide;
                                udtModNameAndResidueLoc[index].ModName = resultModDetails.ModDefinition.MassCorrectionTag;
                                pointerArray[index] = index;
                            }

                            Array.Sort(udtModNameAndResidueLoc, pointerArray, new IModNameAndResidueLocComparer());
                        }

                        // Write out the modifications to the ModDetails file
                        // Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
                        for (var index = 0; index <= searchResult.SearchResultModificationCount - 1; index++)
                        {
                            var resultModDetails = searchResult.GetSearchResultModDetailsByIndex(pointerArray[index]);

                            var modDetailsData = new List<string>
                            {
                                uniqueSeqID.ToString(),
                                resultModDetails.ModDefinition.MassCorrectionTag,
                                resultModDetails.ResidueLocInPeptide.ToString()
                            };

                            mModDetailsFile.WriteLine(CollapseList(modDetailsData));
                        }
                    }
                }
            }

            // Write a new entry to the SeqToProteinMap file if not yet defined
            if (!CheckSeqToProteinMapDefined(uniqueSeqID, searchResult.ProteinName))
            {
                var seqToProteinData = new List<string>
                {
                    uniqueSeqID.ToString(),
                    Convert.ToInt32(searchResult.PeptideCleavageState).ToString(),
                    Convert.ToInt32(searchResult.PeptideTerminusState).ToString(),
                    searchResult.ProteinName,
                    searchResult.ProteinExpectationValue,
                    searchResult.ProteinIntensity
                };

                mSeqToProteinMapFile.WriteLine(CollapseList(seqToProteinData));
            }
        }

        protected void SetErrorCode(ePHRPErrorCodes eNewErrorCode)
        {
            SetErrorCode(eNewErrorCode, false);
        }

        protected void SetErrorCode(ePHRPErrorCodes eNewErrorCode, bool leaveExistingErrorCodeUnchanged)
        {
            if (leaveExistingErrorCodeUnchanged && mErrorCode != ePHRPErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mErrorCode = eNewErrorCode;
            }
        }

        protected void SetErrorMessage(string message, Exception ex = null)
        {
            if (message == null)
                message = string.Empty;

            mErrorMessage = message;
            if (message.Length > 0)
            {
                OnErrorEvent(message, ex);
            }
        }

        protected void ShowPeriodicWarning(int warningCount, int thresholdCountAlwaysShow, string warningMessage)
        {
            if (warningCount <= thresholdCountAlwaysShow ||
                warningCount < 1000 && warningCount % 100 == 0 ||
                warningCount < 10000 && warningCount % 1000 == 0 ||
                warningCount < 100000 && warningCount % 10000 == 0 ||
                warningCount < 1000000 && warningCount % 100000 == 0)
            {
                ReportWarning(warningMessage);
            }
        }

        /// <summary>
        /// If valueText is 0.0, returns 0
        /// If otherwise returns valueText
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="valueText"></param>
        /// <returns></returns>
        protected string TrimZeroIfNotFirstID(int resultID, string valueText)
        {
            return resultID > 1 ? TrimZero(valueText) : valueText;
        }

        /// <summary>
        /// If valueText is 0.0, returns 0
        /// If otherwise returns valueText
        /// </summary>
        /// <param name="valueText"></param>
        /// <returns></returns>
        private string TrimZero(string valueText)
        {
            return valueText.Equals("0.0") ? "0" : valueText;
        }

        /// <summary>
        /// Return the text up to (but not including) the first space in proteinNameAndDescription
        /// </summary>
        /// <param name="proteinNameAndDescription"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected virtual string TruncateProteinName(string proteinNameAndDescription)
        {
            var index = proteinNameAndDescription.IndexOf(' ');
            if (index > 0)
            {
                return proteinNameAndDescription.Substring(0, index);
            }

            return proteinNameAndDescription;
        }

        protected void UpdatePepToProteinMapPeptide(List<udtPepToProteinMappingType> pepToProteinMapping, int index, string peptide)
        {
            var udtItem = pepToProteinMapping[index];
            udtItem.Peptide = peptide;
            pepToProteinMapping[index] = udtItem;
        }

        protected void UpdateProgress(string progressStepDescription)
        {
            UpdateProgress(progressStepDescription, mProgressPercentComplete);
        }

        protected void UpdateProgress(float percentComplete)
        {
            UpdateProgress(ProgressStepDescription, percentComplete);
        }

        protected void UpdateProgress(string progressStepDescription, float percentComplete)
        {
            mProgressStepDescription = string.Copy(progressStepDescription);
            if (percentComplete < 0)
            {
                percentComplete = 0;
            }
            else if (percentComplete > 100)
            {
                percentComplete = 100;
            }
            mProgressPercentComplete = percentComplete;

            OnProgressUpdate(ProgressStepDescription, ProgressPercentComplete);
        }

        /// <summary>
        /// Validate that the specified file exists and has at least one tab-delimited row with a numeric value in the first column
        /// </summary>
        /// <param name="filePath">Path to the file</param>
        /// <param name="fileDescription">File description, e.g. Synopsis</param>
        /// <param name="errorMessage"></param>
        /// <returns>True if the file has data; otherwise false</returns>
        /// <remarks></remarks>
        public static bool ValidateFileHasData(string filePath, string fileDescription, out string errorMessage)
        {
            const int numericDataColIndex = 0;
            return ValidateFileHasData(filePath, fileDescription, out errorMessage, numericDataColIndex);
        }

        /// <summary>
        /// Validate that the specified file exists and has at least one tab-delimited row with a numeric value
        /// </summary>
        /// <param name="filePath">Path to the file</param>
        /// <param name="fileDescription">File description, e.g. Synopsis</param>
        /// <param name="errorMessage"></param>
        /// <param name="numericDataColIndex">Index of the numeric data column; use -1 to simply look for any text in the file</param>
        /// <returns>True if the file has data; otherwise false</returns>
        /// <remarks></remarks>
        public static bool ValidateFileHasData(string filePath, string fileDescription, out string errorMessage, int numericDataColIndex)
        {

            var dataFound = false;

            errorMessage = string.Empty;

            try
            {
                var fileInfo = new FileInfo(filePath);

                if (!fileInfo.Exists)
                {
                    errorMessage = fileDescription + " file not found: " + fileInfo.Name;
                    return false;
                }

                if (fileInfo.Length == 0)
                {
                    errorMessage = fileDescription + " file is empty (zero-bytes)";
                    return false;
                }

                // Open the file and confirm it has data rows
                using (var reader = new StreamReader(new FileStream(fileInfo.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream && !dataFound)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        if (numericDataColIndex < 0)
                        {
                            dataFound = true;
                        }
                        else
                        {
                            // Split on the tab character and check if the first column is numeric
                            var splitLine = lineIn.Split('\t');

                            if (splitLine.Length <= numericDataColIndex)
                                continue;

                            if (double.TryParse(splitLine[numericDataColIndex], out _))
                            {
                                dataFound = true;
                            }
                        }
                    }
                }

                if (!dataFound)
                {
                    errorMessage = fileDescription + " is empty (no data)";
                }

            }
            catch (Exception)
            {
                errorMessage = "Exception validating " + fileDescription + " file";
                return false;
            }

            return dataFound;

        }

        /// <summary>
        /// Compare the two mass values; warn the user if more than 0.1 Da apart (slightly larger threshold if over 5000 Da)
        /// </summary>
        /// <param name="toolName">Tool name</param>
        /// <param name="peptide">Peptide sequence</param>
        /// <param name="peptideMonoMassFromPHRP">Peptide monoisotopic mass, as computed by PHRP</param>
        /// <param name="peptideMonoMassFromTool">Peptide monoisotopic mass, as computed by the search engine</param>
        /// <param name="deltaMassWarningCount">Keeps track of the number of times the monoisotopic mass values differ more than the threshold</param>
        protected void ValidateMatchingMonoisotopicMass(
            string toolName,
            string peptide,
            double peptideMonoMassFromPHRP,
            double peptideMonoMassFromTool,
            ref int deltaMassWarningCount)
        {
            var massDiffThreshold = peptideMonoMassFromTool / 5000 / 10;
            if (massDiffThreshold < 0.1)
                massDiffThreshold = 0.1;

            if (Math.Abs(peptideMonoMassFromPHRP - peptideMonoMassFromTool) <= massDiffThreshold)
                return;

            // Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da
            // (or by a slightly larger value if over 5000 Da)

            // This is unexpected

            string first30Residues;
            if (peptide.Length < 27)
            {
                first30Residues = peptide;
            }
            else
            {
                first30Residues = peptide.Substring(0, 27) + "...";
            }

            deltaMassWarningCount++;
            ShowPeriodicWarning(deltaMassWarningCount,
                                10,
                                string.Format(
                                    "The monoisotopic mass computed by PHRP is more than {0:F2} Da away from " +
                                    "the mass computed by {1}: {2:F4} vs. {3:F4}; peptide {4}",
                                    massDiffThreshold, toolName, peptideMonoMassFromPHRP, peptideMonoMassFromTool, first30Residues));

        }

        private bool ValidatePeptideToProteinMapResults(string peptideToProteinMapFilePath, bool ignorePeptideToProteinMapperErrors)
        {
            bool success;

            var peptideCount = 0;
            var peptideCountNoMatch = 0;
            var linesRead = 0;
            var chSplitChars = new[] { '\t' };

            try
            {
                // Validate that none of the results in peptideToProteinMapFilePath has protein name PROTEIN_NAME_NO_MATCH ( __NoMatch__ )

                var lastPeptide = string.Empty;

                using (var reader = new StreamReader(new FileStream(peptideToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        linesRead += 1;

                        if (linesRead <= 1 || string.IsNullOrEmpty(lineIn))
                            continue;

                        var splitLine = lineIn.Split(chSplitChars, 2);
                        if (splitLine.Length <= 0)
                            continue;

                        if (splitLine[0] != lastPeptide)
                        {
                            peptideCount += 1;
                            lastPeptide = string.Copy(splitLine[0]);
                        }

                        if (lineIn.Contains(PROTEIN_NAME_NO_MATCH))
                        {
                            peptideCountNoMatch += 1;
                        }
                    }
                }

                if (peptideCount == 0)
                {
                    SetErrorMessage("Peptide to protein mapping file is empty: " + peptideToProteinMapFilePath);
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    success = false;
                }
                else if (peptideCountNoMatch == 0)
                {
                    success = true;
                }
                else
                {
                    // Value between 0 and 100
                    var errorPercent = peptideCountNoMatch / (double)peptideCount * 100.0;

                    var message = string.Format("{0:0.00}% of the entries ({1:N0} / {2:N0}) in the peptide to protein map file ({3}) " +
                                                "did not match to a protein in the FASTA file ({4})",
                                                errorPercent, peptideCountNoMatch, peptideCount,
                                                Path.GetFileName(peptideToProteinMapFilePath),
                                                Path.GetFileName(FastaFilePath));

                    if (ignorePeptideToProteinMapperErrors || errorPercent < 0.1)
                    {
                        ReportWarning(message);
                        success = true;
                    }
                    else
                    {
                        SetErrorMessage(message);
                        success = false;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ValidatePeptideToProteinMapResults:" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        protected void ValidatePHRPReaderSupportFiles(string phrpDataFilePath, string outputDirectoryPath)
        {

            try
            {
                if (string.IsNullOrWhiteSpace(outputDirectoryPath))
                    return;

                var phrpDataFile = new FileInfo(phrpDataFilePath);
                var outputDirectory = new DirectoryInfo(outputDirectoryPath);

                if (string.Equals(phrpDataFile.DirectoryName, outputDirectory.FullName, StringComparison.OrdinalIgnoreCase))
                    return;

                var msgfFileName = Path.GetFileName(ReplaceFilenameSuffix(phrpDataFile, FILENAME_SUFFIX_MSGF));

                if (string.IsNullOrWhiteSpace(msgfFileName))
                    return;

                string sourcePath;

                if (string.IsNullOrWhiteSpace(phrpDataFile.DirectoryName))
                    sourcePath = msgfFileName;
                else
                    sourcePath = Path.Combine(phrpDataFile.DirectoryName, msgfFileName);

                var targetPath = Path.Combine(outputDirectory.FullName, msgfFileName);

                if (File.Exists(sourcePath) && !File.Exists(targetPath))
                {
                    File.Copy(sourcePath, targetPath);
                }
            }
            catch (Exception ex)
            {
                ReportWarning("Error in ValidatePHRPReaderSupportFiles: " + ex.Message);
            }
        }

        protected bool ValidateProteinFastaFile(string fastaFilePath)
        {
            var success = ValidateProteinFastaFile(fastaFilePath, out var warningMessage);

            if (!success)
            {
                ReportWarning(warningMessage);
            }

            return success;
        }

        public static bool ValidateProteinFastaFile(string fastaFilePath, out string warningMessage)
        {
            // This RegEx looks for standard amino acids, skipping A, T, C, and G
            var reDefiniteAminoAcid = new Regex(@"[DEFHIKLMNPQRSVWY]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            // This RegEx looks for A, T, C, and G
            var rePotentialNucleicAcid = new Regex(@"[ATCG]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            // This matches any letter
            var reLetter = new Regex(@"[A-Z]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            var validProteinCount = 0;
            var invalidProteinCount = 0;

            try
            {
                warningMessage = string.Empty;

                if (string.IsNullOrEmpty(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "fastaFilePath is not defined in ValidateProteinFastaFile";
                    return false;
                }

                if (!File.Exists(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "Fasta file not found: " + fastaFilePath;
                    return false;
                }

                var fastaFile = new ProteinFileReader.FastaFileReader();
                if (!fastaFile.OpenFile(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "Error opening the fasta file: " + fastaFilePath;
                    return false;
                }

                // Read the first 500 proteins and confirm that each contains amino acid residues
                while (fastaFile.ReadNextProteinEntry())
                {
                    var definiteAminoAcidCount = reDefiniteAminoAcid.Matches(fastaFile.ProteinSequence).Count;
                    var potentialNucleicAcidCount = rePotentialNucleicAcid.Matches(fastaFile.ProteinSequence).Count;
                    var letterCount = reLetter.Matches(fastaFile.ProteinSequence).Count;

                    if (definiteAminoAcidCount > 0.1 * letterCount)
                    {
                        validProteinCount += 1;
                    }
                    else if (potentialNucleicAcidCount > 0.95 * letterCount)
                    {
                        invalidProteinCount += 1;
                    }

                    if (validProteinCount + invalidProteinCount >= 500)
                    {
                        break;
                    }
                }

                if (validProteinCount < invalidProteinCount)
                {
                    Console.WriteLine();
                    warningMessage = "Fasta file contains Nucleic Acids, not Amino Acids: " + Path.GetFileName(fastaFilePath);
                    return false;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine();
                warningMessage = "Exception in ValidateProteinFastaFile: " + ex.Message;
                return false;
            }

            return true;
        }

        private void WriteModDetailsEntry(
            clsPHRPReader reader,
            TextWriter writer,
            IReadOnlyList<udtPepToProteinMappingType> pepToProteinMapping,
            int pepToProteinMapIndex,
            ref int psmCountSkippedSinceReversedOrScrambledProtein)
        {
            foreach (var mod in reader.CurrentPSM.ModifiedResidues)
            {
                var residueLocInProtein = pepToProteinMapping[pepToProteinMapIndex].ResidueStart + mod.ResidueLocInPeptide - 1;
                string residue;

                if (IsLetterAtoZ(mod.Residue))
                {
                    residue = mod.Residue.ToString();
                }
                else
                {
                    var cleanSequence = reader.CurrentPSM.PeptideCleanSequence;

                    if (mod.ResidueLocInPeptide < 1)
                    {
                        // This shouldn't be the case, but we'll check for it anyway
                        residue = cleanSequence.Substring(0, 1);
                    }
                    else if (mod.ResidueLocInPeptide > cleanSequence.Length)
                    {
                        // This shouldn't be the case, but we'll check for it anyway
                        residue = cleanSequence.Substring(cleanSequence.Length - 1, 1);
                    }
                    else
                    {
                        residue = cleanSequence.Substring(mod.ResidueLocInPeptide - 1, 1);
                    }
                }

                if (pepToProteinMapping[pepToProteinMapIndex].Protein == PROTEIN_NAME_NO_MATCH && IsReversedProtein(reader.CurrentPSM.ProteinFirst))
                {
                    // Skip this result
                    psmCountSkippedSinceReversedOrScrambledProtein += 1;
                }
                else
                {
                    writer.WriteLine(reader.CurrentPSM.ResultID + "\t" +
                                     reader.CurrentPSM.Peptide + "\t" +
                                     reader.CurrentPSM.SeqID + "\t" +
                                     pepToProteinMapping[pepToProteinMapIndex].Protein + "\t" +
                                     residue + "\t" +
                                     residueLocInProtein + "\t" +
                                     mod.ModDefinition.MassCorrectionTag + "\t" +
                                     mod.ResidueLocInPeptide + "\t" +
                                     reader.CurrentPSM.MSGFSpecEValue);
                }
            }

        }

        #region "PeptideToProteinMapper Event Handlers"

        private void PeptideToProteinMapper_ProgressChanged(string taskDescription, float percentComplete)
        {
            if (percentComplete >= mNextPeptideToProteinMapperLevel)
            {
                mNextPeptideToProteinMapperLevel += 25;
                UpdateProgress(PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE + percentComplete * (PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE - PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE) / 100);

                if (!HasEventListenerProgressUpdate)
                {
                    Console.WriteLine(" PeptideToProteinMapper is " + percentComplete.ToString("0") + "% complete");
                }
            }
        }

        #endregion

        #region "IComparer classes"

        protected class ISearchOptionModificationInfoComparer : IComparer<udtSearchOptionModificationInfoType>
        {
            public int Compare(udtSearchOptionModificationInfoType x, udtSearchOptionModificationInfoType y)
            {
                if (x.SortOrder > y.SortOrder)
                {
                    return 1;
                }

                if (x.SortOrder < y.SortOrder)
                {
                    return -1;
                }

                if (x.ModificationMass > y.ModificationMass)
                {
                    return 1;
                }

                if (x.ModificationMass < y.ModificationMass)
                {
                    return -1;
                }

                return 0;
            }
        }

        internal class IModNameAndResidueLocComparer : IComparer<udtModNameAndResidueLocType>
        {
            public int Compare(udtModNameAndResidueLocType x, udtModNameAndResidueLocType y)
            {
                if (x.ResidueLocInPeptide > y.ResidueLocInPeptide)
                {
                    return 1;
                }

                if (x.ResidueLocInPeptide < y.ResidueLocInPeptide)
                {
                    return -1;
                }

                if (x.ModName == null)
                    x.ModName = string.Empty;

                if (y.ModName == null)
                    y.ModName = string.Empty;

                return string.Compare(x.ModName, y.ModName, StringComparison.Ordinal);
            }
        }

        protected class PepToProteinMappingComparer : IComparer<udtPepToProteinMappingType>
        {
            public int Compare(udtPepToProteinMappingType x, udtPepToProteinMappingType y)
            {
                var result = string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
                if (result == 0)
                {
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }
                return result;
            }
        }

        private class PepToProteinMappingPeptideSearchComparer : IComparer<udtPepToProteinMappingType>
        {
            public int Compare(udtPepToProteinMappingType x, udtPepToProteinMappingType y)
            {
                return string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
            }
        }

        #endregion

    }
}
