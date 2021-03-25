using System;

namespace PHRPReader
{
    /// <summary>
    /// PHRP Reader error codes
    /// </summary>
    public enum PHRPReaderErrorCodes
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
        /// Input file format not recognized
        /// </summary>
        InputFileFormatNotRecognized = 2,

        /// <summary>
        /// Required input file not found
        /// </summary>
        RequiredInputFileNotFound = 3,

        /// <summary>
        /// Missing .raw or .mzXML file
        /// </summary>
        MissingRawOrMzXmlFile = 4,

        /// <summary>
        /// MSGF program not found
        /// </summary>
        MSGFProgramNotFound = 5,

        /// <summary>
        /// Unspecified error
        /// </summary>
        UnspecifiedError = -1
    }

    /// <summary>
    /// Peptide hit results type
    /// </summary>
    public enum PeptideHitResultTypes
    {
        /// <summary>
        /// Unknown file type
        /// </summary>
        Unknown = 0,

        /// <summary>
        /// SEQUEST
        /// </summary>
        Sequest = 1,

        /// <summary>
        /// XTandem
        /// </summary>
        XTandem = 2,

        /// <summary>
        /// Inspect
        /// </summary>
        Inspect = 3,

        /// <summary>
        /// MSGFDB
        /// </summary>
        [Obsolete("Use MSGFPlus")]
        MSGFDB = 4,

        /// <summary>
        /// MSGFPlus
        /// </summary>
        /// <remarks>Aka MS-GF+</remarks>
        MSGFPlus = 4,

        /// <summary>
        /// MSAlign
        /// </summary>
        MSAlign = 5,

        /// <summary>
        /// MODa
        /// </summary>
        MODa = 6,

        /// <summary>
        /// MODPlus
        /// </summary>
        MODPlus = 7,

        /// <summary>
        /// MSPathFinder
        /// </summary>
        MSPathFinder = 8,

        /// <summary>
        /// TopPIC
        /// </summary>
        TopPIC = 9
    }
}
