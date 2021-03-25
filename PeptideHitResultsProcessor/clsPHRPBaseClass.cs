using System;
using System.IO;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// Old base class for peptide hit results processor classes
    /// </summary>
    [Obsolete("Superseded by Processor.PHRPBaseClass")]
    public abstract class clsPHRPBaseClass : PRISM.EventNotifier
    {
        // Ignore Spelling: Da, A-Za-z, Fscore, prot, mts, MSFragger, xxx

        /// <summary>
        /// Constructor
        /// </summary>
        protected clsPHRPBaseClass()
        {
            //FileDate = "February 3, 2021";

            //mPeptideSeqMassCalculator = new clsPeptideMassCalculator { ChargeCarrierMass = clsPeptideMassCalculator.MASS_PROTON };

            //// Initialize mPeptideMods
            //mPeptideMods = new clsPeptideModificationContainer();

            //// Initialize mUniqueSequences
            //mUniqueSequences = new clsUniqueSequencesContainer();

            //// Initialize mSeqToProteinMap
            //mSeqToProteinMap = new SortedSet<string>();

            //// Define a RegEx to replace all of the non-letter characters
            //mReplaceSymbols = new Regex("[^A-Za-z]", RegexOptions.Compiled);

            //mProteinNameOrder = new Dictionary<string, int>();

            //InitializeLocalVariables();
        }

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

        [Obsolete("Superseded by ResultsFileFormat")]
        public enum PeptideHitResultsFileFormatConstants
        {
            AutoDetermine = 0,
            SequestSynopsisFile = 1,
            SequestFirstHitsFile = 2,
            XTandemXMLFile = 3,
            InspectTXTFile = 4,
            MSGFPlusTXTFile = 5,
            MSAlignTXTFile = 6,
            MODaTXTFile = 7,
            MODPlusTXTFile = 8,
            MSPathFinderTSVFile = 9,
            TopPICTXTFile = 10
        }

        [Obsolete("Superseded by PHRPErrorCode")]
        public enum PHRPErrorCodes
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

        public bool AbortProcessing { get; set; }

        public bool CreateModificationSummaryFile { get; set; }

        public bool CreateFirstHitsFile { get; set; }

        public bool CreateSynopsisFile { get; set; }

        /// <summary>
        /// Create protein mods file
        /// </summary>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the Fasta file specified by mFastaFilePath</remarks>
        public bool CreateProteinModsFile { get; set; }

        // public clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType EnzymeMatchSpec { get; set; }

        public PHRPErrorCodes ErrorCode => 0;

        public string ErrorMessage => string.Empty;

        public string FastaFilePath { get; set; }

        public string FileVersion => string.Empty;

        public string FileDate { get; protected set; }

        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        /// <summary>
        /// Inspect synopsis file p-value threshold
        /// </summary>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float InspectSynopsisFilePValueThreshold { get; set; }

        public string MassCorrectionTagsFilePath { get; set; }

        public string ModificationDefinitionsFilePath { get; set; }

        /// <summary>
        /// Used by clsMODaResultsProcessor and clsMODPlusResultsProcessor
        /// </summary>
        /// <remarks>Higher probability are higher confidence results</remarks>
        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        /// <summary>
        /// Used by MSAlign and TopPIC
        /// </summary>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float MSAlignAndTopPICSynopsisFilePValueThreshold { get; set; }

        /// <summary>
        /// Used by clsMSGFPlusResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <remarks>Lower E-values are higher confidence results</remarks>
        public float MSGFPlusSynopsisFileEValueThreshold { get; set; }

        /// <summary>
        /// clsMSGFPlusResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <remarks>Lower SpecEValue values are higher confidence results</remarks>
        public float MSGFPlusSynopsisFileSpecEValueThreshold { get; set; }

        /// <summary>
        /// Typical non-zero value is 17.0027387
        /// </summary>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideCTerminusMassChange { get; set; }

        /// <summary>
        /// Typical non-zero value is 1.0078246
        /// </summary>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideNTerminusMassChange { get; set; }

        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        public string ProgressStepDescription => string.Empty;

        // ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public float ProgressPercentComplete => 0;

        /// <summary>
        /// Search tool parameter file path
        /// </summary>
        /// <remarks>Used by clsInSpecTResultsProcessor and clsMSGFPlusResultsProcessor (aka SearchEngineParamFileName)</remarks>
        public string SearchToolParameterFilePath { get; set; }

        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        public void AbortProcessingNow()
        {
            AbortProcessing = true;
        }

        public static string AutoDefinePeptideHitResultsFilePath(
            PeptideHitResultsFileFormatConstants peptideHitResultFileFormat,
            string sourceDirectoryPath,
            string baseName)
        {
            if (string.IsNullOrEmpty(baseName))
                return AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath);

            return peptideHitResultFileFormat switch
            {
                PeptideHitResultsFileFormatConstants.SequestFirstHitsFile => Path.Combine(sourceDirectoryPath, baseName + SEQUEST_FIRST_HITS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.SequestSynopsisFile => Path.Combine(sourceDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.XTandemXMLFile => Path.Combine(sourceDirectoryPath, baseName + XTANDEM_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.InspectTXTFile => Path.Combine(sourceDirectoryPath, baseName + INSPECT_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.MSGFPlusTXTFile => Path.Combine(sourceDirectoryPath, baseName + MSGFDB_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.MSAlignTXTFile => Path.Combine(sourceDirectoryPath, baseName + MSALIGN_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.MODaTXTFile => Path.Combine(sourceDirectoryPath, baseName + MODa_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.MODPlusTXTFile => Path.Combine(sourceDirectoryPath, baseName + MODPlus_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.MSPathFinderTSVFile => Path.Combine(sourceDirectoryPath, baseName + MSPathFinder_RESULTS_FILE_SUFFIX),
                PeptideHitResultsFileFormatConstants.TopPICTXTFile => Path.Combine(sourceDirectoryPath, baseName + TopPIC_RESULTS_FILE_SUFFIX),
                _ => AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath)
            };
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
                    matchSpec = index switch
                    {
                        0 => "*" + SEQUEST_SYNOPSIS_FILE_SUFFIX,
                        1 => "*" + SEQUEST_FIRST_HITS_FILE_SUFFIX,
                        2 => "*" + XTANDEM_RESULTS_FILE_SUFFIX,
                        3 => "*" + INSPECT_RESULTS_FILE_SUFFIX,
                        _ => matchSpec
                    };

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
    }
}