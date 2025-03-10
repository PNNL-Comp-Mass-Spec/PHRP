using System;

namespace PeptideHitResultsProcessor
{
    // Ignore Spelling: MaxQuant, MODa

    internal enum MaxQuantModPosition
    {
        Anywhere = 0,
        AnyNterm = 1,
        AnyCterm = 2,
        NotCterm = 3,
        ProteinNterm = 4,
        ProteinCterm = 5
    }

    internal enum MaxQuantModType
    {
        Standard = 0,
        IsobaricLabel = 1,
        Label = 2,
        NeuCodeLabel = 3,
        Glycan = 4,
        AaSubstitution = 5,
        CleavedCrosslink = 6,
        SequenceBasedModifier = 7
    }

    [Obsolete("Deprecated with MaxQuant v2.4.0")]
    internal enum MaxQuantTerminusType
    {
        None = 0,
        NTerminus = 1
    }

    /// <summary>
    /// PHRP Error code
    /// </summary>
    public enum PHRPErrorCode
    {
        /// <summary>
        /// No error
        /// </summary>
        NoError = 0,

        /// <summary>
        /// Invalid input file path
        /// </summary>
        InvalidInputFilePath = 1,

        /// <summary>
        /// Invalid output directory path
        /// </summary>
        InvalidOutputDirectoryPath = 2,

        /// <summary>
        /// Parameter file not found
        /// </summary>
        ParameterFileNotFound = 3,

        /// <summary>
        /// Mass correction tags file not found
        /// </summary>
        MassCorrectionTagsFileNotFound = 4,

        /// <summary>
        /// Modification definition file not found
        /// </summary>
        ModificationDefinitionFileNotFound = 5,

        /// <summary>
        /// Error reading the input file
        /// </summary>
        ErrorReadingInputFile = 6,

        /// <summary>
        /// Error creating output files
        /// </summary>
        ErrorCreatingOutputFiles = 7,

        /// <summary>
        /// Error reading the parameter file
        /// </summary>
        ErrorReadingParameterFile = 8,

        /// <summary>
        /// Error reading the mass correction tags file
        /// </summary>
        ErrorReadingMassCorrectionTagsFile = 9,

        /// <summary>
        /// Error reading the modification definitions file
        /// </summary>
        ErrorReadingModificationDefinitionsFile = 10,

        /// <summary>
        /// File path error
        /// </summary>
        FilePathError = 11,

        /// <summary>
        /// Unspecified error
        /// </summary>
        UnspecifiedError = -1
    }

    /// <summary>
    /// MS/MS search tool results file format
    /// </summary>
    public enum ResultsFileFormat
    {
        /// <summary>
        /// Auto-determine the result type
        /// </summary>
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
        /// InSpecT
        /// </summary>
        InspectTXTFile = 4,

        /// <summary>
        /// MSGFDB or MS-GF+ (MSGF+) .tsv file
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
        TopPICTXTFile = 10,

        /// <summary>
        /// MaxQuant
        /// </summary>
        MaxQuantTXTFile = 11,

        /// <summary>
        /// MSFragger
        /// </summary>
        MSFraggerTSVFile = 12,

        /// <summary>
        /// DIA-NN 1.x
        /// </summary>
        DiannTSVFile = 13,

        /// <summary>
        /// DIA-NN 2.x
        /// </summary>
        DiannParquetFile = 14
    }
}
