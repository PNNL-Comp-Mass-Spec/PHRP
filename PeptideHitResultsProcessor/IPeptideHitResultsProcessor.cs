using System;
using System.Collections.Generic;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// Defines minimum required functionality for classes that will process peptide hit results files
    /// </summary>
    /// <remarks>Note that the peptide hit results file format will determined using AnalysisToolName</remarks>
    public abstract class IPeptideHitResultsProcessor : PRISM.EventNotifier
    {
        #region "Enums"

        /// <summary>
        /// /Return values for Start and Abort functions
        /// </summary>
        public enum ProcessResults
        {
            /// <summary>
            /// Spectra file creation aborted
            /// </summary>
            PH_ABORTED = -2,

            /// <summary>
            /// Operation failed
            /// </summary>
            PH_FAILURE = -1,

            /// <summary>
            /// Operation succeeded
            /// </summary>
            PH_SUCCESS = 0,
        }

        /// <summary>
        /// Return value for status property
        /// </summary>
        public enum ProcessStatus
        {
            /// <summary>
            /// Plugin initialization in progress
            /// </summary>
            PH_STARTING,

            /// <summary>
            /// Plugin is attempting to do its job
            /// </summary>
            PH_RUNNING,

            /// <summary>
            /// Plugin successfully completed its job
            /// </summary>
            PH_COMPLETE,

            /// <summary>
            /// There was an error somewhere
            /// </summary>
            PH_ERROR,

            /// <summary>
            /// An ABORT command has been received; plugin shutdown in progress
            /// </summary>
            PH_ABORTING
        }
        #endregion

        #region "Structures"
        public struct InitializationParams
        {
            /// <summary>
            /// Source directory path
            /// </summary>
            public string SourceDirectoryPath;

            /// <summary>
            /// Output directory path
            /// </summary>
            public string OutputDirectoryPath;

            /// <summary>
            /// Source directory path
            /// </summary>
            [Obsolete("Use SourceDirectoryPath")]
            public string SourceFolderPath
            {
                get => SourceDirectoryPath;
                set => SourceDirectoryPath = value;
            }

            /// <summary>
            /// Output directory path
            /// </summary>
            [Obsolete("Use OutputDirectoryPath")]
            public string OutputFolderPath
            {
                get => OutputDirectoryPath;
                set => OutputDirectoryPath = value;
            }

            /// <summary>
            /// Peptide hit results filename
            /// </summary>
            /// <remarks>If this is empty, it will be auto-defined using: DatasetName + TOOL_NAME_RESULTS_FILE_SUFFIX</remarks>
            public string PeptideHitResultsFileName;

            /// <summary>
            /// Mass correction tags filename
            /// </summary>
            /// <remarks>
            /// If this is empty, it will be auto-defined as "Mass_Correction_Tags.txt"
            /// </remarks>
            public string MassCorrectionTagsFileName;

            /// <summary>
            /// Modification definitions filename
            /// </summary>
            /// <remarks>
            /// If this is empty, it will be auto-defined using:
            /// Path.GetFileNameWithoutExtension(ParameterFileName) + MODIFICATION_DEFINITIONS_FILE_SUFFIX
            /// where MODIFICATION_DEFINITIONS_FILE_SUFFIX = "_ModDefs.txt"
            /// </remarks>
            public string ModificationDefinitionsFileName;

            /// <summary>
            /// Miscellaneous parameters
            /// </summary>
            public Dictionary<string, string> MiscParams;

            /// <summary>
            /// Debug level
            /// </summary>
            /// <remarks>0=minimum, 5=maximum verbosity</remarks>
            public int DebugLevel;

            /// <summary>
            /// Tool name
            /// </summary>
            public string AnalysisToolName;

            /// <summary>
            /// Dataset name
            /// </summary>
            public string DatasetName;

            /// <summary>
            /// Parameter file for the job
            /// </summary>
            public string ParameterFileName;

            /// <summary>
            /// Contains XML settings for PeptideHitResultsProcessor
            /// </summary>
            public string SettingsFileName;

            /// <summary>
            /// When true, create the Inspect synopsis file
            /// </summary>
            public bool CreateSynopsisFile;

            /// <summary>
            /// When true, create the Inspect first hits file
            /// </summary>
            public bool CreateFirstHitsFile;
        }
        #endregion

        #region "Properties"

        /// <summary>
        /// Input: path to directory containing the peptide hit results file
        /// </summary>
        public virtual string SourceDirectoryPath { get; set; }

        /// <summary>
        /// Input: path to directory where processed file is to be placed
        /// </summary>
        public virtual string OutputDirectoryPath { get; set; }

        /// <summary>
        /// Input: path to directory containing the peptide hit results file
        /// </summary>
        [Obsolete("Use SourceDirectoryPath")]
        public virtual string SourceFolderPath
        {
            get => SourceDirectoryPath;
            set => SourceDirectoryPath = value;
        }

        /// <summary>
        /// Input: path to directory where processed file is to be placed
        /// </summary>
        [Obsolete("Use OutputDirectoryPath")]
        public virtual string OutputFolderPath
        {
            get => OutputDirectoryPath;
            set => OutputDirectoryPath = value;
        }

        /// <summary>
        /// Input: Filename containing the peptide hit results;
        /// See comment for PeptideHitResultsFileName in Structure InitializationParams above
        /// </summary>
        public virtual string PeptideHitResultsFileName { get; set; }

        /// <summary>
        /// Input: Filename containing the global mass correction tags list;
        /// See comment for PeptideHitResultsFileName in Structure InitializationParams above
        /// </summary>
        public virtual string MassCorrectionTagsFileName { get; set; }

        /// <summary>
        /// Input: Filename containing modification definitions for this peptide hit results file;
        /// See comment for PeptideHitResultsFileName in Structure InitializationParams above
        /// </summary>
        public virtual string ModificationDefinitionsFileName { get; set; }

        public virtual string AnalysisToolName { get; set; }

        public virtual string DatasetName { get; set; }

        /// <summary>
        /// Peptide search tool parameter file name
        /// </summary>
        public virtual string ParameterFileName { get; set; }

        /// <summary>
        /// XML settings file with section PeptideHitResultsProcessorOptions
        /// </summary>
        public virtual string SettingsFileName { get; set; }

        public virtual bool CreateSynopsisFile { get; set; }

        public virtual bool CreateFirstHitsFile { get; set; }

        /// <summary>
        /// For passing miscellaneous parameters (not presently used)
        /// </summary>
        [Obsolete("Unused")]
        public virtual Dictionary<string, string> MiscParams { get; set; }

        /// <summary>
        /// Allows calling program to get current status
        /// </summary>
        public abstract ProcessStatus Status { get; }

        /// <summary>
        /// Allows calling program to determine if processing succeeded
        /// </summary>
        public abstract ProcessResults Results { get; }

        /// <summary>
        /// Error message describing any errors encountered
        /// </summary>
        public abstract string ErrMsg { get; }

        /// <summary>
        /// Progress indicator, value between 0 and 100
        /// </summary>
        public abstract float PercentComplete { get; }

        /// <summary>
        /// Allows control of debug information verbosity
        /// </summary>
        /// <remarks>0=minimum, 5=maximum verbosity</remarks>
        public virtual int DebugLevel { get; set; }

        #endregion

        #region "Methods"

        /// <summary>
        /// Initializes parameters. Must be called before executing Start()
        /// </summary>
        /// <param name="options"></param>
        public abstract void Setup(InitializationParams options);

        /// <summary>
        /// Starts the spectra file creation process
        /// </summary>
        public abstract ProcessStatus Start();

        /// <summary>
        /// Aborts spectra file creation
        /// </summary>
        public abstract ProcessStatus Abort();

        #endregion
    }
}
