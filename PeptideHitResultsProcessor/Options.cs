﻿using PeptideHitResultsProcessor.Processor;
using PHRPReader;
using PRISM;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// PHRP Options
    /// </summary>
    public class PHRPOptions
    {
        // Ignore Spelling: conf, Dia, DMS, Fasta, Fragger, Hyperscore, MaxQuant, MODa, prepended, Quant, txt

        /// <summary>
        /// Default database connection string
        /// </summary>
        public const string DMS_CONNECTION_STRING = "host=prismdb2.emsl.pnl.gov;Port=5432;Database=dms;UserId=dmsreader;Password=dms4fun";

        /// <summary>
        /// Input file path
        /// </summary>
        /// <remarks>This path can contain wildcard characters, e.g. C:\*.raw</remarks>
        [Option("InputFilePath", "InputDataFilePath", "I",
            ArgPosition = 1, Required = true, HelpShowsDefault = false, IsInputFilePath = true,
            HelpText = "The name of a file or directory to process; the path can contain the wildcard character *\n" +
                       "Either define this at the command line using /I or in a parameter file")]
        public string InputFilePath { get; set; }

        /// <summary>
        /// Output directory path
        /// </summary>
        [Option("OutputDirectoryPath", "OutputDirectory", "O", HelpShowsDefault = false,
            HelpText = "Output directory name (or full path)\n" +
                       "If omitted, the output files will be created in the program directory")]
        public string OutputDirectoryPath { get; set; }

        /// <summary>
        /// Base name for output files
        /// </summary>
        [Option("OutputFileBaseName", "OutputBaseName", "BaseName", "B", HelpShowsDefault = false,
            HelpText = "Output file base name\n" +
                       "If an empty string, the output file base name will be auto-defined based on the input file name.\n" +
                       "For MaxQuant and MSFragger results, if the input file has multiple datasets, " +
                       "the base name will be either the base name defined here, " +
                       "or the longest string in common for the dataset names.\n" +
                       "This parameter should not be used when the input file spec is a wildcard.")]
        public string OutputFileBaseName { get; set; }

        /// <summary>
        /// Alternate base path for input files
        /// </summary>
        [Option("AlternateBasePath", "AltBasePath", "BasePath", HelpShowsDefault = false,
            HelpText = "Alternative base path for input files\n" +
                       "The program first looks for input files using the path provided for a given parameter. " +
                       "If the file is not found, the base path specified here is prepended to the provided path " +
                       "to look for the file in an alternative location " +
                       @"(not valid for rooted paths like C:\Work\Dataset.tsv or \\Server\Share\Dataset.tsv)")]
        public string AlternateBasePath { get; set; }

        /// <summary>
        /// Typically a Key=Value parameter file path, but can also be a legacy XML-based parameter file
        /// </summary>
        /// <remarks>
        /// This property is intended for use when using PeptideToProteinMapEngine.dll along with a parameter file
        /// For PeptideHitResultsProcRunner.exe, specify the Key=Value parameter file using /Conf or /P
        /// </remarks>
        public string ParameterFilePath { get; set; }

        /// <summary>
        /// When true, recurse subdirectories
        /// </summary>
        /// <remarks>
        /// This will be auto-set to true if MaxLevelsToRecurse is defined in the parameter file or if /S is used at the command line
        /// This functionality is enabled by the ArgExistsProperty option
        /// </remarks>
        public bool RecurseDirectories { get; set; }

        /// <summary>
        /// Process files in subdirectories
        /// </summary>
        [Option("MaxLevelsToRecurse", "S", ArgExistsProperty = nameof(RecurseDirectories),
            HelpShowsDefault = false, SecondaryArg = true,
            HelpText = "If supplied, process all valid files in the input directory and subdirectories\n" +
                       "Include a number after /S (like /S:2) to limit the level of subdirectories to examine (0 means to recurse infinitely)\n" +
                       "The equivalent notation in a parameter file is MaxLevelsToRecurse=2")]
        public int MaxLevelsToRecurse { get; set; }

        /// <summary>
        /// When true, create a log file
        /// </summary>
        /// <remarks>
        /// This is auto-set to true if /L or /LogDir is used at the command line, or if
        /// LogfilePath or LogDirectory have a name defined in a parameter file
        /// </remarks>
        public bool LogMessagesToFile { get; set; }

        /// <summary>
        /// Log file path
        /// </summary>
        [Option("LogFilePath", "LogFile", "Log", "L",
            ArgExistsProperty = nameof(LogMessagesToFile), HelpShowsDefault = false,
            HelpText = "Log file path.\n" +
                       "Use /L at the command line to log messages to a file whose name is auto-defined using the current date, " +
                       "or use /L:LogFileName.txt to specify the name.\n" +
                       "In a Key=Value parameter file, define a file name or path to enable logging to a file.")]
        public string LogFilePath { get; set; }

        /// <summary>
        /// Log directory path
        /// </summary>
        [Option("LogDirectoryPath", "LogDirectory", "LogDir",
            ArgExistsProperty = nameof(LogMessagesToFile), HelpShowsDefault = false,
            HelpText = "The directory where the log file should be written")]
        public string LogDirectoryPath { get; set; }

        /// <summary>
        /// Mass correction tags file path
        /// </summary>
        [Option("MassCorrectionTagsFile", "MassCorrectionTags", "T",
            HelpShowsDefault = false, IsInputFilePath = true,
            HelpText = "Tab-delimited text file that lists modification names and masses.\n" +
                       "The first column has the modification name (aka mass correction tag name), " +
                       "the second column has the monoisotopic mass, and the optional third column shows the affected atom (for isotopic mods).\n" +
                       "An example file is at https://github.com/PNNL-Comp-Mass-Spec/PHRP/blob/master/Data/Mass_Correction_Tags.txt")]
        public string MassCorrectionTagsFilePath { get; set; }

        /// <summary>
        /// Modification definitions file path (DMS mod names)
        /// </summary>
        [Option("ModificationDefinitionsFile", "ModDefsFile", "M",
            HelpShowsDefault = false, IsInputFilePath = true,
            HelpText = "Tab-delimited text file that defines modification symbols, mod masses, and mod names.\n" +
                       "The first column has the modification symbol, the second column shows the modification mass, " +
                       "and the optional third column lists the residues that can be modified with the given mass " +
                       "(1 letter residue symbols, no need to separate with commas or spaces).\n" +
                       "If the file has a header line, it can also include columns listing the modification type, " +
                       "mass correction tag name, UniMod name, and MaxQuant mod name.\n" +
                       "An example file is at https://github.com/PNNL-Comp-Mass-Spec/PHRP/blob/master/Data/Example_ModDefs.txt")]
        public string ModificationDefinitionsFilePath { get; set; }

        /// <summary>
        /// Search tool parameter file path (aka SearchEngineParamFileName)
        /// </summary>
        /// <remarks>Used by MSGFPlusResultsProcessor, MaxQuantResultsProcessor, and others</remarks>
        [Option("SearchToolParameterFile", "ToolParamFile", "N",
            HelpShowsDefault = false, IsInputFilePath = true,
            HelpText = "The parameter file provided to the search tool.\n" +
                       "This is used when processing results from MS-GF+, MSPathFinder, MaxQuant, MODa, MODPlus, MSAlign, MSFragger, TopPIC, and InSpecT.\n" +
                       "For MaxQuant, provide either an XML-based parameter file (root element is <MaxQuantParams>) " +
                       "or provide the parameters.txt file created in the txt results directory.\n" +
                       "The XML-based parameter file is preferred, since it is required to allow PHRP " +
                       "to accurately compute monoisotopic masses of peptides identified by MaxQuant.")]
        public string SearchToolParameterFilePath { get; set; }

        /// <summary>
        /// FASTA file path
        /// </summary>
        [Option("FastaFile", "FASTA", "F",
            HelpShowsDefault = false, IsInputFilePath = true,
            HelpText = "FASTA file path. The order of the proteins in the FASTA file " +
                       "dictates which protein is listed for each peptide in the First Hits file.")]
        public string FastaFilePath { get; set; }

        /// <summary>
        /// Create modification summary file
        /// </summary>
        [Option("CreateModificationSummaryFile", "CreateModSummaryFile",
            HelpText = "When true, create the _ModSummary.txt file")]
        public bool CreateModificationSummaryFile { get; set; }

        /// <summary>
        /// Create protein mods file
        /// </summary>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the FASTA file specified by FastaFilePath</remarks>
        [Option("CreateProteinModsFile", "ProteinMods",
            HelpText = "When true, create the _ProteinMods.txt file.\n" +
                       "This requires that either an existing _PepToProtMapMTS.txt file exists, " +
                       "or that the FASTA file be defined using /F or FastaFilePath.")]
        public bool CreateProteinModsFile { get; set; }

        /// <summary>
        /// If true, assume that the _ProteinMods.txt file has reversed protein sequences
        /// </summary>
        [Option("ProteinModsFileIncludesReversedProteins", "ProteinModsIncludeReversed",
            HelpText = "Set this to true if an existing _ProteinMods.txt file has reversed protein sequences, " +
                       "or if the FASTA file has reversed proteins.\n" +
                       "If false, will skip reversed proteins when creating the _ProteinMods.txt file.")]
        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        /// <summary>
        /// Use existing MTS peptide to protein map file
        /// </summary>
        [Option("UseExistingPepToProteinMapFile",
            HelpText = "When true, look for an existing _PepToProtMap.txt file; if not found, it will be created using the FASTA file.")]
        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        /// <summary>
        /// Create protein mods using existing PHRP data file
        /// </summary>
        [Option("CreateProteinModsUsingPHRPDataFile", "CreateProteinModsViaPHRP",
            HelpText = "When true, create the _ProteinMods.txt file using existing PHRP data files.\n" +
                       "This requires that either an existing _PepToProtMapMTS.txt file exist, " +
                       "or that the FASTA file be defined using /F or FastaFilePath.")]
        public bool CreateProteinModsUsingPHRPDataFile { get; set; }

        /// <summary>
        /// When true, ignore peptide to protein mapping errors
        /// </summary>
        /// <remarks>
        /// Example error that may be reported if this is false:
        /// 0.43% of the entries (96 / 22,127) in the peptide to protein map file (Dataset_maxq_syn_PepToProtMap.txt)
        /// did not match to a protein in the FASTA file (Proteins_2020-10-21.fasta)
        /// </remarks>
        [Option("IgnorePepToProtMapErrors",
            HelpText = "When true, ignore peptide to protein mapping errors")]
        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        /// <summary>
        /// Create first hits file
        /// </summary>
        [Option("CreateFirstHitsFile", "FHT",
            HelpText = "When true, create the first hits file (_fht.txt)")]
        public bool CreateFirstHitsFile { get; set; }

        /// <summary>
        /// Create synopsis file
        /// </summary>
        [Option("CreateSynopsisFile", "Syn",
            HelpText = "When true, create the synopsis file (_syn.txt)")]
        public bool CreateSynopsisFile { get; set; }

        // Default filter thresholds, by tool

        // | Tool         | Filter Rules                 |
        // |--------------|------------------------------|
        // | DIA-NN       | Keep if QValue < 0.10        |
        // |              |      or CScore < 0.25        |
        // |--------------|------------------------------|
        // | MaxQuant     | Keep if Andromeda Score > 50 |
        // |              |      or PepValue < 0.01      |
        // |--------------|------------------------------|
        // | MSFragger    | Keep if E-Value < 0.75       |
        // |              |      or Hyperscore > 20      |
        // |--------------|------------------------------|
        // | MODa and     | Keep if probability > 0.05   |
        // | MODPlus      |                              |
        // |--------------|------------------------------|
        // | MSGFPlus     | Keep if Spec_E-value < 5E-7  |
        // |              |      or E-Value < 0.75       |
        // |              |      or Q-Value < 0.01  (1%) |
        // |--------------|------------------------------|
        // | MSPathFinder | Keep if SpecEValue < 5E-07   |
        // |              |      or QValue < 0.10        |
        // |--------------|------------------------------|
        // | TopPIC       | Keep if p-value < 0.95       |
        // |              |      or E-value < 0.95       |
        // |--------------|------------------------------|
        // |  X!Tandem    | Keep all results             |
        // |              |                              |
        // |--------------|------------------------------|

        /// <summary>
        /// DIA-NN Q-Value threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Q-Value is below the DiaNNQValue threshold, or if its confidence score (CScore) is above the DiaNNCScore threshold
        /// </remarks>
        [Option("DiaNNQValueThreshold", "DiaNNQValue",
            HelpText = "When processing DIA-NN results, the Q-Value threshold used to determine which peptides are written to the synopsis file.\n" +
                       "A PSM is stored if its Q-Value is below the DiaNNQValue threshold, or if its confidence score (CScore) is above the DiaNNCScore threshold.")]
        public double DiaNNQValueThreshold { get; set; }

        /// <summary>
        /// DIA-NN confidence score (CScore) threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Q-Value is below the DiaNNQValue threshold, or if its confidence score (CScore) is above the DiaNNCScore threshold
        /// </remarks>
        [Option("DiaNNCScoreThreshold", "DiaNNCScore",
            HelpText = "When processing DIA-NN results, the confidence score threshold used to determine which peptides are written to the synopsis file.\n" +
                       "A PSM is stored if its Q-Value is below the DiaNNQValue threshold, or if its confidence score (CScore) is above the DiaNNCScore threshold.")]
        public double DiaNNConfidenceScoreThreshold { get; set; }

        /// <summary>
        /// MaxQuant Andromeda score threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
        /// </remarks>
        [Option("MaxQuantAndromedaScoreThreshold", "MaxQScore",
            HelpText = "When processing MaxQuant results, the Andromeda score threshold used to determine which peptides are written to the synopsis file.\n" +
                       "A PSM is stored if its Andromeda score is above the MaxQScore threshold, or if its PEP score is below the MaxQPEP threshold.")]
        public int MaxQuantAndromedaScoreThreshold { get; set; }

        /// <summary>
        /// MaxQuant Posterior Error Probability (PEP) score threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
        /// </remarks>
        [Option("MaxQuantPosteriorErrorProbabilityThreshold", "MaxQPEP",
            HelpText = "When processing MaxQuant results, the Posterior Error Probability (PEP) score threshold used to determine which peptides are written to the synopsis file.\n" +
                       "A PSM is stored if its Andromeda score is above the MaxQScore threshold, or if its PEP score is below the MaxQPEP threshold.")]
        public float MaxQuantPosteriorErrorProbabilityThreshold { get; set; }

        /// <summary>
        /// MSFragger Hyperscore threshold to use when creating the synopsis file
        /// </summary>
        [Option("MSFraggerHyperscoreThreshold", "MaxHyperscore",
            HelpText = "When processing MSFragger results, the Expectation score (E-Value) used to determine which peptides are written to the synopsis file\n" +
                       "Smaller E-values are higher confidence results.\n" +
                       "Filter passing peptides have E-value less than 0.75 or Hyperscore greater than 20")]
        public int MSFraggerHyperscoreThreshold { get; set; }

        /// <summary>
        /// MSGFPlusResultsProcessor and MSPathFinderResultsProcessor
        /// </summary>
        /// <remarks>Smaller SpecEValue values are higher confidence results</remarks>
        [Option("MSGFPlusSynopsisFileSpecEValueThreshold", "MSGFPlusSpecEValue",
            HelpText = "When processing MS-GF+ results, the spec E-value threshold used to determine which peptides are written to the synopsis file.\n" +
                       "Smaller spec E-values are higher confidence results.")]
        public float MSGFPlusSynopsisFileSpecEValueThreshold { get; set; }

        /// <summary>
        /// Used by MSGFPlusResultsProcessor and MSFragger
        /// </summary>
        /// <remarks>Smaller E-values are higher confidence results</remarks>
        [Option("MSGFPlusSynopsisFileEValueThreshold", "MSGFPlusEValue",
            HelpText = "When processing MS-GF+ results, the E-value threshold used to determine which peptides are written to the synopsis file.\n" +
                       "Smaller E-values are higher confidence results.\n" +
                       "Filter passing peptides have Spec E-value less than 5E-7 or E-Value (EValue) less than 0.75 or Q-Value (QValue) less than 1%.")]
        public float MSGFPlusSynopsisFileEValueThreshold { get; set; }

        /// <summary>
        /// Used by MODaResultsProcessor and MODPlusResultsProcessor
        /// </summary>
        /// <remarks>Larger probability values are higher confidence results</remarks>
        [Option("MODaMODPlusSynopsisFileProbabilityThreshold", "SynProb",
            HelpText = "When processing a MODPlus or MODa results, the probability threshold used to determine which peptides are written to the synopsis file.\n" +
                       "Larger probability values are higher confidence results, thus the default of 0.05 is a very loose filter.")]
        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        /// <summary>
        /// Used by MSAlign and TopPIC
        /// </summary>
        /// <remarks>Smaller p-values are higher confidence results</remarks>
        [Option("MSAlignAndTopPICSynopsisFilePValueThreshold", "SynPValue",
            HelpText = "When processing a MODPlus or MODa results, the p-value threshold used to determine which peptides are written to the synopsis file.\n" +
                       "Smaller p-values are higher confidence results, thus the default of 0.95 is a very loose filter.")]
        public float MSAlignAndTopPICSynopsisFilePValueThreshold { get; set; }

        /// <summary>
        /// DMS database connection string (only used if the computer is on the pnl.gov domain)
        /// </summary>
        /// <remarks>Set this to an empty string or to "false" to disable contacting DMS</remarks>
        [Option("DMSConnectionString", "DB",
            HelpText = "DMS database connection string. Set this to an empty string or to 'false' to disable contacting DMS.")]
        public string DMSConnectionString { get; set; }

        /// <summary>
        /// Enzyme match specification
        /// </summary>
        public PeptideCleavageStateCalculator.EnzymeMatchSpecInfo EnzymeMatchSpec { get; internal set; }

        private string mEnzymeMatchSpecLeftResidue;
        private string mEnzymeMatchSpecRightResidue;

        /// <summary>
        /// RegEx for the residue to the left of a cleavage point
        /// </summary>
        [Option("EnzymeMatchSpecLeftResidue",
            HelpText = "Regular expression for the residue to the left of cleavage points for the given enzyme",
            SecondaryArg = true)]
        public string EnzymeMatchSpecLeftResidue
        {
            get => mEnzymeMatchSpecLeftResidue;
            set
            {
                mEnzymeMatchSpecLeftResidue = value;
                UpdateEnzymeMatchSpecIfDefined();
            }
        }

        /// <summary>
        /// RegEx for the residue to the right of a cleavage point
        /// </summary>
        [Option("EnzymeMatchSpecRightResidue",
            HelpText = "Regular expression for the residue to the right of cleavage points for the given enzyme",
            SecondaryArg = true)]
        public string EnzymeMatchSpecRightResidue
        {
            get => mEnzymeMatchSpecRightResidue;
            set
            {
                mEnzymeMatchSpecRightResidue = value;
                UpdateEnzymeMatchSpecIfDefined();
            }
        }

        /// <summary>
        /// Mass to add to the N-terminus of peptides
        /// </summary>
        /// <remarks>Typical non-zero value is 1.0078246</remarks>
        [Option("PeptideNTerminusMassChange",
            HelpText = "Peptide N-terminus mass to add to peptides; ignored if 0",
            SecondaryArg = true)]
        public double PeptideNTerminusMassChange { get; set; }

        /// <summary>
        /// Mass to add to the C-terminus of peptides
        /// </summary>
        /// <remarks>Typical non-zero value is 17.0027387</remarks>
        [Option("PeptideCTerminusMassChange",
            HelpText = "Peptide C-terminus mass to add to peptides; ignored if 0",
            SecondaryArg = true)]
        public double PeptideCTerminusMassChange { get; set; }

        /// <summary>
        /// InSpecT synopsis file p-value threshold
        /// </summary>
        /// <remarks>
        /// Smaller p-values are higher confidence results
        /// Peptides with a TotalPRMScore >= 50 or an FScore >= 0 will also be included in the synopsis file
        /// </remarks>
        [Option("InspectSynopsisFilePValueThreshold", "InsSynPValue",
            HelpText = "When processing InSpecT results, the PValue threshold used to determine which peptides are written to the synopsis file",
            Hidden = true)]
        public float InspectSynopsisFilePValueThreshold { get; set; }

        /// <summary>
        /// Legacy XML-based parameter file
        /// </summary>
        [Option("XmlParameterFile", "XmlParamFile",
            HelpText = "Legacy XML-based parameter file; deprecated with Peptide Hit Results Processor 3.1",
            Hidden = true)]
        public string XmlParameterFile { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public PHRPOptions()
        {
            MassCorrectionTagsFilePath = string.Empty;
            ModificationDefinitionsFilePath = string.Empty;
            SearchToolParameterFilePath = string.Empty;
            FastaFilePath = string.Empty;

            CreateModificationSummaryFile = true;

            CreateProteinModsFile = false;
            ProteinModsFileIncludesReversedProteins = false;
            UseExistingMTSPepToProteinMapFile = false;
            CreateProteinModsUsingPHRPDataFile = false;

            IgnorePeptideToProteinMapperErrors = false;

            CreateFirstHitsFile = true;
            CreateSynopsisFile = true;

            DiaNNQValueThreshold = DiaNNResultsProcessor.DEFAULT_QVALUE_THRESHOLD;
            DiaNNConfidenceScoreThreshold = DiaNNResultsProcessor.DEFAULT_CONFIDENCE_SCORE_THRESHOLD;

            MaxQuantAndromedaScoreThreshold = MaxQuantResultsProcessor.DEFAULT_ANDROMEDA_SCORE_THRESHOLD;
            MaxQuantPosteriorErrorProbabilityThreshold = MaxQuantResultsProcessor.DEFAULT_PEP_THRESHOLD;

            MSFraggerHyperscoreThreshold = MSFraggerResultsProcessor.DEFAULT_HYPERSCORE_THRESHOLD;

            MSGFPlusSynopsisFileEValueThreshold = MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MSGFPlusSynopsisFileSpecEValueThreshold = MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            MODaMODPlusSynopsisFileProbabilityThreshold = MODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            MSAlignAndTopPICSynopsisFilePValueThreshold = MSAlignResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            DMSConnectionString = DMS_CONNECTION_STRING;

            EnzymeMatchSpecLeftResidue = PeptideCleavageStateCalculator.TRYPSIN_LEFT_RESIDUE_REGEX;
            EnzymeMatchSpecRightResidue = PeptideCleavageStateCalculator.TRYPSIN_RIGHT_RESIDUE_REGEX;

            // Note that updating EnzymeMatchSpecLeftResidue and EnzymeMatchSpecRightResidue will result in EnzymeMatchSpec being auto-defined
            // Thus, this should always evaluate to false
            if (string.IsNullOrWhiteSpace(EnzymeMatchSpec.LeftResidueRegEx) || string.IsNullOrWhiteSpace(EnzymeMatchSpec.RightResidueRegEx))
            {
                EnzymeMatchSpec = PeptideCleavageStateCalculator.GetDefaultEnzymeMatchSpec();
            }

            PeptideNTerminusMassChange = PeptideMassCalculator.DEFAULT_N_TERMINUS_MASS_CHANGE;
            PeptideCTerminusMassChange = PeptideMassCalculator.DEFAULT_C_TERMINUS_MASS_CHANGE;

            InspectSynopsisFilePValueThreshold = InSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            XmlParameterFile = string.Empty;
        }

        /// <summary>
        /// Copy options from sourceOptions to the current instance of PHRPOptions
        /// </summary>
        /// <param name="sourceOptions">Source options</param>
        public void UpdateAll(PHRPOptions sourceOptions)
        {
            InputFilePath = sourceOptions.InputFilePath;

            OutputDirectoryPath = sourceOptions.OutputDirectoryPath;
            OutputFileBaseName = sourceOptions.OutputFileBaseName;

            AlternateBasePath = sourceOptions.AlternateBasePath;
            ParameterFilePath = sourceOptions.ParameterFilePath;

            RecurseDirectories = sourceOptions.RecurseDirectories;
            MaxLevelsToRecurse = sourceOptions.MaxLevelsToRecurse;

            LogFilePath = sourceOptions.LogFilePath;
            LogDirectoryPath = sourceOptions.LogDirectoryPath;
            LogMessagesToFile = sourceOptions.LogMessagesToFile;

            MassCorrectionTagsFilePath = sourceOptions.MassCorrectionTagsFilePath;
            ModificationDefinitionsFilePath = sourceOptions.ModificationDefinitionsFilePath;
            SearchToolParameterFilePath = sourceOptions.SearchToolParameterFilePath;
            FastaFilePath = sourceOptions.FastaFilePath;

            CreateModificationSummaryFile = sourceOptions.CreateModificationSummaryFile;

            CreateProteinModsFile = sourceOptions.CreateProteinModsFile;
            ProteinModsFileIncludesReversedProteins = sourceOptions.ProteinModsFileIncludesReversedProteins;
            UseExistingMTSPepToProteinMapFile = sourceOptions.UseExistingMTSPepToProteinMapFile;

            CreateFirstHitsFile = sourceOptions.CreateFirstHitsFile;
            CreateSynopsisFile = sourceOptions.CreateSynopsisFile;

            DMSConnectionString = sourceOptions.DMSConnectionString;

            EnzymeMatchSpecLeftResidue = sourceOptions.EnzymeMatchSpecLeftResidue;
            EnzymeMatchSpecRightResidue = sourceOptions.EnzymeMatchSpecLeftResidue;

            PeptideNTerminusMassChange = sourceOptions.PeptideNTerminusMassChange;
            PeptideCTerminusMassChange = sourceOptions.PeptideCTerminusMassChange;

            IgnorePeptideToProteinMapperErrors = sourceOptions.IgnorePeptideToProteinMapperErrors;

            InspectSynopsisFilePValueThreshold = sourceOptions.InspectSynopsisFilePValueThreshold;
            XmlParameterFile = sourceOptions.XmlParameterFile;

            MaxQuantAndromedaScoreThreshold = sourceOptions.MaxQuantAndromedaScoreThreshold;
            MaxQuantPosteriorErrorProbabilityThreshold = sourceOptions.MaxQuantPosteriorErrorProbabilityThreshold;

            MSFraggerHyperscoreThreshold = sourceOptions.MSFraggerHyperscoreThreshold;

            MODaMODPlusSynopsisFileProbabilityThreshold = sourceOptions.MODaMODPlusSynopsisFileProbabilityThreshold;

            MSAlignAndTopPICSynopsisFilePValueThreshold = sourceOptions.MSAlignAndTopPICSynopsisFilePValueThreshold;

            MSGFPlusSynopsisFileEValueThreshold = sourceOptions.MSGFPlusSynopsisFileEValueThreshold;

            MSGFPlusSynopsisFileSpecEValueThreshold = sourceOptions.MSGFPlusSynopsisFileSpecEValueThreshold;
        }

        private void UpdateEnzymeMatchSpecIfDefined()
        {
            if (string.IsNullOrWhiteSpace(mEnzymeMatchSpecLeftResidue) || string.IsNullOrWhiteSpace(mEnzymeMatchSpecRightResidue))
            {
                return;
            }

            EnzymeMatchSpec = new PeptideCleavageStateCalculator.EnzymeMatchSpecInfo(mEnzymeMatchSpecLeftResidue, mEnzymeMatchSpecRightResidue);
        }

        /// <summary>
        /// Verify that an input file is defined
        /// </summary>
        /// <returns>True if options are valid, otherwise false</returns>
        public bool Validate()
        {
            // ReSharper disable once InvertIf
            if (string.IsNullOrWhiteSpace(InputFilePath))
            {
                ConsoleMsgUtils.ShowWarning("Input file not defined; cannot continue");
                return false;
            }

            return true;
        }
    }
}
