namespace PeptideHitResultsProcessor
{
    // Ignore Spelling: MODa

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
}
