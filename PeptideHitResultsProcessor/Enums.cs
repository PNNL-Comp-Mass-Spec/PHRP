using System;

namespace PeptideHitResultsProcessor
{
    public static class Enums
    {
        // Ignore Spelling: MODa

        public enum ResultsFileFormat
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

        public enum PHRPErrorCode
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
    }
}
