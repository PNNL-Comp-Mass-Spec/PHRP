using System;

namespace PHRPReader
{
    public class Enums
    {
#pragma warning disable 1591

        /// <summary>
        /// Peptide hit results type
        /// </summary>
        public enum PeptideHitResultTypes
        {
            Unknown = 0,
            Sequest = 1,
            XTandem = 2,
            Inspect = 3,
            [Obsolete("Use MSGFPlus")]
            MSGFDB = 4,
            MSGFPlus = 4,      // Aka MS-GF+
            MSAlign = 5,
            MODa = 6,
            MODPlus = 7,
            MSPathFinder = 8,
            TopPIC = 9
        }

        /// <summary>
        /// PHRP Reader error codes
        /// </summary>
        public enum PHRPReaderErrorCodes
        {
            NoError = 0,
            InvalidInputFilePath = 1,
            InputFileFormatNotRecognized = 2,
            RequiredInputFileNotFound = 3,
            MissingRawOrMzXmlFile = 4,
            MSGFProgramNotFound = 5,
            UnspecifiedError = -1
        }

#pragma warning restore 1591
    }
}
