Option Strict On

' This class implements the IPeptideHitResultsProcessor interface 
' It uses class clsXTandemResultsProcessor to read an XTandem search results XML file 
'  or class clsSequestResultsProcessor to read a Sequest Synopsis or First Hits file
'
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Copyright 2006, Battelle Memorial Institute
' Started January 6, 2006

Imports System.IO

Public Class clsAnalysisManagerPeptideHitResultsProcessor
    Implements IPeptideHitResultsProcessor


#Region "Constants and enums"
    Protected Const DEFAULT_MASS_CORRECTION_TAGS_FILENAME As String = "Mass_Correction_Tags.txt"
    Protected Const MODIFICATION_DEFINITIONS_FILE_SUFFIX As String = "_ModDefs.txt"
#End Region

#Region "Classwide variables"

    Protected m_SourceFolderPath As String = String.Empty
    Protected m_OutFolderPath As String = String.Empty
    Protected m_PeptideHitResultsFileName As String = String.Empty
    Protected m_PeptideHitResultsFileFormat As clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine

    Protected m_MassCorrectionTagsFileName As String = String.Empty
    Protected m_ModificationDefinitionsFileName As String = String.Empty

	Protected m_MiscParams As Dictionary(Of String, String)
    Protected m_DebugLevel As Integer = 0

    Protected m_AnalysisToolName As String = String.Empty
    Protected m_DSName As String = String.Empty
    Protected m_ParameterFileName As String = String.Empty      ' Peptide search tool parameter file name
    Protected m_SettingsFileName As String = String.Empty       ' XML settings file with section PeptideHitResultsProcessorOptions

    Protected m_ParameterFilePath As String = String.Empty      ' Peptide search tool parameter file name
    Protected m_SettingsFilePath As String = String.Empty       ' XML settings file with section PeptideHitResultsProcessorOptions

    Protected m_CreateInspectFirstHitsFile As Boolean
    Protected m_CreateInspectSynopsisFile As Boolean

    Protected m_PeptideHitResultsFilePath As String = String.Empty
    Protected m_MassCorrectionTagsFilePath As String = String.Empty
    Protected m_ModificationDefinitionsFilePath As String = String.Empty

    Protected m_ErrMsg As String = String.Empty

    Protected m_Status As IPeptideHitResultsProcessor.ProcessStatus
    Protected m_Results As IPeptideHitResultsProcessor.ProcessResults

    Protected WithEvents m_PeptideHitResultsProcessor As clsPHRPBaseClass

	Protected m_thThread As Threading.Thread

#End Region

#Region "Events"
    Public Event ErrorOccurred(strMessage As String) Implements IPeptideHitResultsProcessor.ErrorOccurred
    Public Event DebugEvent(strMessage As String) Implements IPeptideHitResultsProcessor.DebugEvent

    ' PercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Public Event ProgressChanged(taskDescription As String, percentComplete As Single) Implements IPeptideHitResultsProcessor.ProgressChanged

#End Region

#Region "Properties"
    Public Property AnalysisToolName() As String Implements IPeptideHitResultsProcessor.AnalysisToolName
        Get
            Return m_AnalysisToolName
        End Get
        Set(Value As String)
            m_AnalysisToolName = Value
        End Set
    End Property

    Public Property CreateInspectFirstHitsFile() As Boolean Implements IPeptideHitResultsProcessor.CreateInspectFirstHitsFile
        Get
            Return m_CreateInspectFirstHitsFile
        End Get
        Set(value As Boolean)
            m_CreateInspectFirstHitsFile = value
        End Set
    End Property

    Public Property CreateInspectSynopsisFile() As Boolean Implements IPeptideHitResultsProcessor.CreateInspectSynopsisFile
        Get
            Return m_CreateInspectSynopsisFile
        End Get
        Set(value As Boolean)
            m_CreateInspectSynopsisFile = value
        End Set
    End Property

    Public Property DatasetName() As String Implements IPeptideHitResultsProcessor.DatasetName
        Get
            Return m_DSName
        End Get
        Set(Value As String)
            m_DSName = Value
        End Set
    End Property

    Public Property DebugLevel() As Integer Implements IPeptideHitResultsProcessor.DebugLevel
        Get
            Return m_DebugLevel
        End Get
        Set(Value As Integer)
            m_DebugLevel = Value
        End Set
    End Property

    Public ReadOnly Property ErrMsg() As String Implements IPeptideHitResultsProcessor.ErrMsg
        Get
            Return m_ErrMsg
        End Get
    End Property

    Public Property MassCorrectionTagsFileName() As String Implements IPeptideHitResultsProcessor.MassCorrectionTagsFileName
        Get
            Return m_MassCorrectionTagsFileName
        End Get
        Set(Value As String)
            m_MassCorrectionTagsFileName = Value
        End Set
    End Property

    Public WriteOnly Property MiscParams() As Dictionary(Of String, String) Implements IPeptideHitResultsProcessor.MiscParams
        Set(Value As Dictionary(Of String, String))
            m_MiscParams = Value
        End Set
    End Property

    Public Property ModificationDefinitionsFileName() As String Implements IPeptideHitResultsProcessor.ModificationDefinitionsFileName
        Get
            Return m_ModificationDefinitionsFileName
        End Get
        Set(Value As String)
            m_ModificationDefinitionsFileName = Value
        End Set
    End Property

    Public Property OutputFolderPath() As String Implements IPeptideHitResultsProcessor.OutputFolderPath
        Get
            Return m_OutFolderPath
        End Get
        Set(Value As String)
            m_OutFolderPath = Value
        End Set
    End Property

    Public Property ParameterFileName() As String Implements IPeptideHitResultsProcessor.ParameterFileName
        Get
            Return m_ParameterFileName
        End Get
        Set(Value As String)
            m_ParameterFileName = Value
        End Set
    End Property

    Public Property PeptideHitResultsFileName() As String Implements IPeptideHitResultsProcessor.PeptideHitResultsFileName
        Get
            Return m_PeptideHitResultsFileName
        End Get
        Set(Value As String)
            m_PeptideHitResultsFileName = Value
        End Set
    End Property

    Public ReadOnly Property PercentComplete() As Single Implements IPeptideHitResultsProcessor.PercentComplete
        Get
            If Not m_PeptideHitResultsProcessor Is Nothing Then
                Return m_PeptideHitResultsProcessor.ProgressPercentComplete
            End If
        End Get
    End Property

    Public Property SettingsFileName() As String Implements IPeptideHitResultsProcessor.SettingsFileName
        Get
            Return m_SettingsFileName
        End Get
        Set(Value As String)
            m_SettingsFileName = Value
        End Set
    End Property

    Public Property SourceFolderPath() As String Implements IPeptideHitResultsProcessor.SourceFolderPath
        Get
            Return m_SourceFolderPath
        End Get
        Set(Value As String)
            m_SourceFolderPath = Value
        End Set
    End Property

    Public ReadOnly Property Status() As IPeptideHitResultsProcessor.ProcessStatus Implements IPeptideHitResultsProcessor.Status
        Get
            Return m_Status
        End Get
    End Property

    Public ReadOnly Property Results() As IPeptideHitResultsProcessor.ProcessResults Implements IPeptideHitResultsProcessor.Results
        Get
            Return m_Results
        End Get
    End Property
#End Region

    Public Function Abort() As IPeptideHitResultsProcessor.ProcessStatus Implements IPeptideHitResultsProcessor.Abort
        If Not m_PeptideHitResultsProcessor Is Nothing Then
            m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ABORTING
            m_PeptideHitResultsProcessor.AbortProcessingNow()
        End If
    End Function

    Public Sub Setup(InitParams As IPeptideHitResultsProcessor.InitializationParams) Implements IPeptideHitResultsProcessor.Setup

        'Copies all input data required for plugin operation to appropriate memory variables
        With InitParams
            m_SourceFolderPath = .SourceFolderPath
            m_OutFolderPath = .OutputFolderPath

            m_PeptideHitResultsFileName = .PeptideHitResultsFileName
            m_MassCorrectionTagsFileName = .MassCorrectionTagsFileName
            m_ModificationDefinitionsFileName = .ModificationDefinitionsFileName

            m_MiscParams = .MiscParams
            m_DebugLevel = .DebugLevel

            m_AnalysisToolName = .AnalysisToolName
            m_DSName = .DatasetName
            m_ParameterFileName = .ParameterFileName
            m_SettingsFileName = .SettingsFileName

            m_CreateInspectFirstHitsFile = .CreateInspectFirstHitsFile
            m_CreateInspectSynopsisFile = .CreateInspectSynopsisFile

        End With

    End Sub

    Public Function Start() As IPeptideHitResultsProcessor.ProcessStatus Implements IPeptideHitResultsProcessor.Start

        m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_STARTING

        'Verify necessary files are in specified locations
        If Not InitSetup() Then
            m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_FAILURE
            m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR
            Return m_Status
        End If

        ' Process the results file (the process runs in a separate thread)
        m_Status = ProcessPeptideHitResultsFile()

        If m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR Then
            m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_FAILURE
        End If

        Return m_Status

    End Function

    ''' <summary>
    ''' Initializes m_PeptideHitResultsProcessor then starts a separate thread to process the file
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Overridable Function ProcessPeptideHitResultsFile() As IPeptideHitResultsProcessor.ProcessStatus

        Try
            'Initialize m_PeptideHitResultsProcessor
            Select Case m_PeptideHitResultsFileFormat
                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.XTandemXMLFile
                    m_PeptideHitResultsProcessor = New clsXTandemResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile, clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestSynopsisFile
                    m_PeptideHitResultsProcessor = New clsSequestResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InSpectTXTFile
                    m_PeptideHitResultsProcessor = New clsInSpecTResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile
                    m_PeptideHitResultsProcessor = New clsMSGFDBResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSAlignTXTFile
                    m_PeptideHitResultsProcessor = New clsMSAlignResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODaTXTFile
                    m_PeptideHitResultsProcessor = New clsMODaResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODPlusTXTFile
                    m_PeptideHitResultsProcessor = New clsMODaResultsProcessor

                Case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile
                    m_PeptideHitResultsProcessor = New clsMSPathFinderResultsProcessor

                Case Else
                    ' Unknown format; cannot continue
                    LogErrors("ProcessPeptideHitResultsFile", "Unknown peptide hit results file format: " & m_PeptideHitResultsFileFormat.ToString, Nothing)
                    m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR
                    Exit Try
            End Select

            ' Define the auxiliary file paths
            With m_PeptideHitResultsProcessor
                .MassCorrectionTagsFilePath = m_MassCorrectionTagsFilePath
                .ModificationDefinitionsFilePath = m_ModificationDefinitionsFilePath
                .SearchToolParameterFilePath = m_ParameterFilePath
                .InspectSynopsisFilePValueThreshold = PeptideHitResultsProcessor.clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD

                .CreateInspectFirstHitsFile = m_CreateInspectFirstHitsFile
                .CreateInspectSynopsisFile = m_CreateInspectSynopsisFile
            End With

            m_thThread = New System.Threading.Thread(AddressOf ProcessPeptideHitResultsFileWork)
            m_thThread.Start()

            m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_RUNNING

        Catch ex As Exception
            LogErrors("ProcessPeptideHitResultsFile", "Error initializing and running m_PeptideHitResultsProcessor", ex)
            m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR
        End Try

        Return m_Status
    End Function

    Protected Overridable Sub ProcessPeptideHitResultsFileWork()

        Dim blnSuccess As Boolean

        Try
            blnSuccess = m_PeptideHitResultsProcessor.ProcessFile(m_PeptideHitResultsFilePath, m_OutFolderPath, m_SettingsFilePath)
        Catch ex As Exception
            blnSuccess = False
        End Try

        If blnSuccess Then
            m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_COMPLETE
        Else
            If m_PeptideHitResultsProcessor.AbortProcessing Then
                LogErrors("ProcessPeptideHitResultsFileWork", "Processing aborted", Nothing)
                m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_ABORTED
                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ABORTING
            Else
                LogErrors("ProcessPeptideHitResultsFileWork", m_PeptideHitResultsProcessor.ErrorMessage(), Nothing)
                m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_FAILURE
                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR
            End If
        End If

    End Sub

    Protected Overridable Function InitSetup() As Boolean

        'Initializes module variables and verifies mandatory parameters have been propery specified

        'Output folder name
        If m_OutFolderPath = "" Then
            m_ErrMsg = "Output folder path not specified"
            Return False
        End If

        'Source folder name
        If m_SourceFolderPath = "" Then
            m_ErrMsg = "Source folder not specified"
            Return False
        End If

        If Me.DebugLevel >= 3 Then
            RaiseEvent DebugEvent("Setup params: OutFolderPath = " & m_OutFolderPath)
            RaiseEvent DebugEvent("Setup params: SourceFolderPath = " & m_SourceFolderPath)
        End If

        'Source directory exists?
        If Not VerifyDirExists(m_SourceFolderPath) Then Return False 'Error msg handled by VerifyDirExists

        'Output directory exists?
        If Not VerifyDirExists(m_OutFolderPath) Then Return False 'Error msg handled by VerifyDirExists

        'Analysis tool name defined?
        If String.IsNullOrWhiteSpace(m_AnalysisToolName) Then
            m_ErrMsg = "Analysis tool name not specified"
            Return False
        End If

        'Dataset name defined?
        If String.IsNullOrWhiteSpace(m_DSName) Then
            m_ErrMsg = "Dataset name not specified"
            Return False
        End If

        'Settings file name defined?
        If String.IsNullOrWhiteSpace(m_SettingsFileName) Then
            m_ErrMsg = "Settings file name not specified"
            Return False
        End If

        'Parameter file name defined?
        If String.IsNullOrWhiteSpace(m_ParameterFileName) Then
            m_ErrMsg = "Parameter file name not specified"
            Return False
        End If

        'Define the parameter file path; this is passed as the search tool parameter file
        m_ParameterFilePath = Path.Combine(m_SourceFolderPath, m_ParameterFileName)

        'Define the settings file path; this is passed as the parameter file name to m_PeptideHitResultsProcessor
        m_SettingsFilePath = Path.Combine(m_SourceFolderPath, m_SettingsFileName)

        'Define the peptide hit results format based on the analysis tool name
        If m_AnalysisToolName.IndexOf("xtandem", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.XTandemXMLFile

        ElseIf m_AnalysisToolName.IndexOf("sequest", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestSynopsisFile

        ElseIf m_AnalysisToolName.IndexOf("inspect", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InSpectTXTFile

        ElseIf m_AnalysisToolName.IndexOf("msgfdb", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile

        ElseIf m_AnalysisToolName.IndexOf("msalign", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSAlignTXTFile

        ElseIf m_AnalysisToolName.IndexOf("moda", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODaTXTFile

        ElseIf m_AnalysisToolName.IndexOf("modplus", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODPlusTXTFile

        ElseIf m_AnalysisToolName.IndexOf("mspathfinder", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile

        ElseIf m_AnalysisToolName.IndexOf("dataextractor", StringComparison.InvariantCultureIgnoreCase) >= 0 Then
            ' Data Extractor step-tool; we'll need to auto-determine the results format
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine

        Else
            ' Unrecognized analysis tool name
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine
        End If

        If Me.DebugLevel >= 3 Then
            RaiseEvent DebugEvent("Setup params: AnalysisToolName = " & m_AnalysisToolName)
            RaiseEvent DebugEvent("Setup params: PeptideHitResultsFileFormat = " & m_PeptideHitResultsFileFormat.ToString)

            RaiseEvent DebugEvent("Setup params: DSName = " & m_DSName)
            RaiseEvent DebugEvent("Setup params: SettingsFilePath = " & m_SettingsFilePath)
            RaiseEvent DebugEvent("Setup params: ParameterFilePath = " & m_ParameterFilePath)
        End If

        'Define the peptide hit results file name
        If String.IsNullOrWhiteSpace(m_PeptideHitResultsFileName) Then
            m_PeptideHitResultsFilePath = clsPHRPBaseClass.AutoDefinePeptideHitResultsFilePath(m_PeptideHitResultsFileFormat, m_SourceFolderPath, m_DSName)
        Else
            m_PeptideHitResultsFilePath = Path.Combine(m_SourceFolderPath, m_PeptideHitResultsFileName)
        End If

        If Me.DebugLevel >= 3 Then
            RaiseEvent DebugEvent("Setup params: PeptideHitResultsFilePath = " & m_PeptideHitResultsFilePath)
        End If

        'Now that m_PeptideHitResultsFilePath has been determined, if m_PeptideHitResultsFileFormat is .AutoDetermine then try to determine the correct format
        If m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine Then
            m_PeptideHitResultsFileFormat = clsPHRPBaseClass.DetermineResultsFileFormat(m_PeptideHitResultsFilePath)
        End If

        'Define the mass correction tags file path
        If String.IsNullOrWhiteSpace(m_MassCorrectionTagsFileName) Then
            m_MassCorrectionTagsFilePath = Path.Combine(m_SourceFolderPath, DEFAULT_MASS_CORRECTION_TAGS_FILENAME)
        Else
            m_MassCorrectionTagsFilePath = Path.Combine(m_SourceFolderPath, m_MassCorrectionTagsFileName)
        End If

        'Define the modification definitions file path
        If String.IsNullOrWhiteSpace(m_ModificationDefinitionsFileName) Then
            m_ModificationDefinitionsFilePath = Path.Combine(m_SourceFolderPath, Path.GetFileNameWithoutExtension(m_ParameterFileName) & MODIFICATION_DEFINITIONS_FILE_SUFFIX)
        Else
            m_ModificationDefinitionsFilePath = Path.Combine(m_SourceFolderPath, m_ModificationDefinitionsFileName)
        End If

        If Me.DebugLevel >= 3 Then
            RaiseEvent DebugEvent("Setup params: PeptideHitResultsFileFormat = " & m_PeptideHitResultsFileFormat.ToString)
            RaiseEvent DebugEvent("Setup params: MassCorrectionTagsFilePath = " & m_MassCorrectionTagsFilePath)
            RaiseEvent DebugEvent("Setup params: ModificationDefinitionsFilePath = " & m_ModificationDefinitionsFilePath)
        End If

        'Parameter file exists?
        If Not VerifyFileExists(m_ParameterFilePath) Then Return False 'Error msg handled by VerifyFileExists

        'Settings file exists?
        If Not VerifyFileExists(m_SettingsFilePath) Then Return False 'Error msg handled by VerifyFileExists

        'Peptide hit results file exists?
        If Not VerifyFileExists(m_PeptideHitResultsFilePath) Then Return False 'Error msg handled by VerifyFileExists

        'Modification definitions file exists?
        If Not VerifyFileExists(m_ModificationDefinitionsFilePath) Then Return False 'Error msg handled by VerifyFileExists

        'Mass correction tags file exists?
        If Not VerifyFileExists(m_MassCorrectionTagsFilePath) Then Return False 'Error msg handled by VerifyFileExists

        'If we got here, everything's OK
        Return True

    End Function

    Private Sub LogErrors(strSource As String, strMessage As String, ex As Exception, Optional blnLogLocalOnly As Boolean = True)

        m_ErrMsg = String.Copy(strMessage).Replace(ControlChars.NewLine, "; ")

        If ex Is Nothing Then
            ex = New System.Exception("Error")
        Else
            If Not ex.Message Is Nothing AndAlso ex.Message.Length > 0 Then
                m_ErrMsg &= "; " & ex.Message
            End If
        End If

        Trace.WriteLine(System.DateTime.Now().ToLongTimeString & "; " & m_ErrMsg, strSource)
        Console.WriteLine(System.DateTime.Now().ToLongTimeString & "; " & m_ErrMsg, strSource)

        RaiseEvent ErrorOccurred(m_ErrMsg)

    End Sub

    Protected Sub UpdateProgress(strProgressStepDescription As String, sngPercentComplete As Single)
        Static strProgressStepDescriptionSaved As String = String.Empty
        Static sngProgressPercentComplete As Single = 0

        Dim blnDescriptionChanged As Boolean = False

        If strProgressStepDescription <> strProgressStepDescriptionSaved Then
            blnDescriptionChanged = True
        End If

        strProgressStepDescriptionSaved = String.Copy(strProgressStepDescription)
        If sngPercentComplete < 0 Then
            sngPercentComplete = 0
        ElseIf sngPercentComplete > 100 Then
            sngPercentComplete = 100
        End If
        sngProgressPercentComplete = sngPercentComplete

        If blnDescriptionChanged And Me.DebugLevel >= 2 Then
            If Math.Abs(sngProgressPercentComplete) < Single.Epsilon Then
                RaiseEvent DebugEvent(strProgressStepDescriptionSaved)
            Else
                RaiseEvent DebugEvent(strProgressStepDescriptionSaved & " (" & sngProgressPercentComplete.ToString("0.0") & "% complete)")
            End If
        End If

        RaiseEvent ProgressChanged(strProgressStepDescription, sngPercentComplete)
    End Sub

    Protected Overridable Function VerifyDirExists(TestDir As String) As Boolean

        'Verifies that the specified directory exists
        If Directory.Exists(TestDir) Then
            m_ErrMsg = ""
            Return True
        Else
            m_ErrMsg = "Directory " & TestDir & " not found"
            Return False
        End If

    End Function

    Protected Overridable Function VerifyFileExists(TestFile As String) As Boolean
        'Verifies specified file exists
        If File.Exists(TestFile) Then
            m_ErrMsg = ""
            Return True
        Else
            m_ErrMsg = "File " & TestFile & " not found"
            Return False
        End If

    End Function

    Private Sub m_PeptideHitResultsProcessor_ErrorOccurred(ErrorMessage As String) Handles m_PeptideHitResultsProcessor.ErrorOccurred
        LogErrors("PeptideHitResultsProcessor", ErrorMessage, Nothing, True)
    End Sub

    Private Sub mPeptideHitResultsProcessor_ProgressChanged(taskDescription As String, percentComplete As Single) Handles m_PeptideHitResultsProcessor.ProgressChanged
        UpdateProgress(taskDescription, percentComplete)
    End Sub

    Private Sub mPeptideHitResultsProcessor_ProgressComplete() Handles m_PeptideHitResultsProcessor.ProgressComplete
        ' OperationComplete()
    End Sub

    Private Sub mPeptideHitResultsProcessor_ProgressReset() Handles m_PeptideHitResultsProcessor.ProgressReset
        ' ResetProgress()
    End Sub

End Class
