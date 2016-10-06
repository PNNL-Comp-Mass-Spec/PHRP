
Public Interface IPeptideHitResultsProcessor
    'Defines minimum required functionality for classes that will process peptide hit results files
    'Note that the peptide hit results file format will determined using AnalysisToolName

#Region "Enums"
    Enum ProcessResults
        'Return values for Start and Abort functions
        PH_SUCCESS = 0  'Operation succeeded
        PH_FAILURE = -1  'Operation failed
        PH_ABORTED = -2  'Spectra file creation aborted
    End Enum

    Enum ProcessStatus
        'Return value for status property
        PH_STARTING  'Plugin initialization in progress
        PH_RUNNING  'Plugin is attempting to do its job
        PH_COMPLETE  'Plugin successfully completed its job
        PH_ERROR  'There was an error somewhere
        PH_ABORTING  'An ABORT command has been received; plugin shutdown in progress
    End Enum
#End Region

#Region "Structures"
    Structure InitializationParams
        Dim SourceFolderPath As String
        Dim OutputFolderPath As String
        Dim PeptideHitResultsFileName As String             ' If this is empty then it will be auto-defined using: DatasetName & XTANDEM_RESULTS_FILE_SUFFIX
        Dim MassCorrectionTagsFileName As String            ' If this is empty then it will be auto-defined using: Const DEFAULT_MASS_CORRECTION_TAGS_FILENAME = "Mass_Correction_Tags.txt"
        Dim ModificationDefinitionsFileName As String       ' If this is empty then it will be auto-defined using: Path.GetFileNameWithoutExtension(ParameterFileName) & MODIFICATION_DEFINITIONS_FILE_SUFFIX where MODIFICATION_DEFINITIONS_FILE_SUFFIX = "_ModDefs.txt"
        Dim MiscParams As Dictionary(Of String, String)
        Dim DebugLevel As Integer

        Dim AnalysisToolName As String
        Dim DatasetName As String
        Dim ParameterFileName As String                     ' Parameter file for the job
        Dim SettingsFileName As String                      ' Contains XML settings for PeptideHitResultsProcessor

        Dim CreateInspectSynopsisFile As Boolean
        Dim CreateInspectFirstHitsFile As Boolean
    End Structure
#End Region

#Region "Events"
    Event ErrorOccurred(strMessage As String)
    Event DebugEvent(strMessage As String)

    ' PercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Event ProgressChanged(taskDescription As String, percentComplete As Single)
#End Region

#Region "Properties"
    Property SourceFolderPath() As String ' (in) – path to folder containing the peptide hit results file
    Property OutputFolderPath() As String ' (in) – path to folder where processed file is to be placed
    Property PeptideHitResultsFileName() As String ' (in) - Filename containing the peptide hit results; See comment for PeptideHitResultsFileName in Structure InitializationParams above
    Property MassCorrectionTagsFileName() As String ' (in) - Filename containing the global mass correction tags list; See comment for PeptideHitResultsFileName in Structure InitializationParams above
    Property ModificationDefinitionsFileName() As String ' (in) - Filename containing modification definitions for this peptide hit results file; See comment for PeptideHitResultsFileName in Structure InitializationParams above

    Property AnalysisToolName() As String
    Property DatasetName() As String
    Property ParameterFileName() As String
    Property SettingsFileName() As String

    Property CreateInspectSynopsisFile() As Boolean
    Property CreateInspectFirstHitsFile() As Boolean

    WriteOnly Property MiscParams() As Dictionary(Of String, String)    'For passing miscelleneous parameters (not presently used)
    ReadOnly Property Status() As ProcessStatus 'Allows calling program to get current status
    ReadOnly Property Results() As ProcessResults  'Allows calling program to determine if processing succeeded
    ReadOnly Property ErrMsg() As String  'Error message describing any errors encountered
    ReadOnly Property PercentComplete() As Single   ' Progress indicator, value between 0 and 100
    Property DebugLevel() As Integer 'Allows control of debug information verbosity; 0=minimum, 5=maximum verbosity

#End Region

#Region "Methods"
    Sub Setup(InitParams As InitializationParams)  'Initializes parameters. Must be called before executing Start()

    Function Start() As ProcessStatus  'Starts the spectra file creation process

    Function Abort() As ProcessStatus  'Aborts spectra file creation
#End Region

End Interface
