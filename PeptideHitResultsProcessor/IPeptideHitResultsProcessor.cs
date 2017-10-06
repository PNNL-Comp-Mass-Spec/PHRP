using System.Collections.Generic;

namespace PeptideHitResultsProcessor
{
    public abstract class IPeptideHitResultsProcessor
    {
        //Defines minimum required functionality for classes that will process peptide hit results files
        //Note that the peptide hit results file format will determined using AnalysisToolName

        #region "Enums"
        public enum ProcessResults
        {
            //Return values for Start and Abort functions
            PH_SUCCESS = 0,   //Operation succeeded
            PH_FAILURE = -1,  //Operation failed
            PH_ABORTED = -2   //Spectra file creation aborted
        }

        public enum ProcessStatus
        {
            //Return value for status property
            PH_STARTING,  //Plugin initialization in progress
            PH_RUNNING,   //Plugin is attempting to do its job
            PH_COMPLETE,  //Plugin successfully completed its job
            PH_ERROR,     //There was an error somewhere
            PH_ABORTING   //An ABORT command has been received; plugin shutdown in progress
        }
        #endregion

        #region "Structures"
        public struct InitializationParams
        {
            public string SourceFolderPath;
            public string OutputFolderPath;
            public string PeptideHitResultsFileName;             // If this is empty then it will be auto-defined using: DatasetName & XTANDEM_RESULTS_FILE_SUFFIX
            public string MassCorrectionTagsFileName;            // If this is empty then it will be auto-defined using: Const DEFAULT_MASS_CORRECTION_TAGS_FILENAME = "Mass_Correction_Tags.txt"
            public string ModificationDefinitionsFileName;       // If this is empty then it will be auto-defined using: Path.GetFileNameWithoutExtension(ParameterFileName) & MODIFICATION_DEFINITIONS_FILE_SUFFIX where MODIFICATION_DEFINITIONS_FILE_SUFFIX = "_ModDefs.txt"
            public Dictionary<string, string> MiscParams;
            public int DebugLevel;

            public string AnalysisToolName;
            public string DatasetName;
            public string ParameterFileName;                     // Parameter file for the job
            public string SettingsFileName;                      // Contains XML settings for PeptideHitResultsProcessor

            public bool CreateInspectSynopsisFile;
            public bool CreateInspectFirstHitsFile;
        }
        #endregion

        #region "Events"
        public abstract event ErrorOccurredEventHandler ErrorOccurred;
        public delegate void ErrorOccurredEventHandler(string strMessage);
        public abstract event DebugEventEventHandler DebugEvent;
        public delegate void DebugEventEventHandler(string strMessage);

        // PercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public abstract event ProgressChangedEventHandler ProgressChanged;
        public delegate void ProgressChangedEventHandler(string taskDescription, float percentComplete);
        #endregion

        #region "Properties"
        public abstract string SourceFolderPath { get; set; } // (in) - path to folder containing the peptide hit results file
        public abstract string OutputFolderPath { get; set; } // (in) - path to folder where processed file is to be placed
        public abstract string PeptideHitResultsFileName { get; set; } // (in) - Filename containing the peptide hit results; See comment for PeptideHitResultsFileName in Structure InitializationParams above
        public abstract string MassCorrectionTagsFileName { get; set; } // (in) - Filename containing the global mass correction tags list; See comment for PeptideHitResultsFileName in Structure InitializationParams above
        public abstract string ModificationDefinitionsFileName { get; set; } // (in) - Filename containing modification definitions for this peptide hit results file; See comment for PeptideHitResultsFileName in Structure InitializationParams above

        public abstract string AnalysisToolName { get; set; }
        public abstract string DatasetName { get; set; }
        public abstract string ParameterFileName { get; set; }
        public abstract string SettingsFileName { get; set; }

        public abstract bool CreateInspectSynopsisFile { get; set; }
        public abstract bool CreateInspectFirstHitsFile { get; set; }

        public abstract Dictionary<string, string> MiscParams { get; set; }    //For passing miscelleneous parameters (not presently used)
        public abstract ProcessStatus Status { get; } //Allows calling program to get current status
        public abstract ProcessResults Results { get; }  //Allows calling program to determine if processing succeeded
        public abstract string ErrMsg { get; }  //Error message describing any errors encountered
        public abstract float PercentComplete { get; }   // Progress indicator, value between 0 and 100
        public abstract int DebugLevel { get; set; } //Allows control of debug information verbosity; 0=minimum, 5=maximum verbosity

        #endregion

        #region "Methods"
        public abstract void Setup(InitializationParams InitParams);  //Initializes parameters. Must be called before executing Start()

        public abstract ProcessStatus Start();  //Starts the spectra file creation process

        public abstract ProcessStatus Abort();  //Aborts spectra file creation
        #endregion
    }
}
