'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class reads a tab-delimited text file (created by the Peptide File Extractor or by PHRP)
' and returns the data for each peptide hit search result
'
' It also integrates MSGF results with the peptide hit search results
' And, it integrates scan stats values (to determine elution time)
' 
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsModificationDefinition
Imports System.IO
Imports System.Text.RegularExpressions

Public Class clsPHRPReader
    Implements IDisposable

#Region "Constants"

    Public Const N_TERMINAL_PEPTIDE_SYMBOL_DMS As Char = "<"c
    Public Const C_TERMINAL_PEPTIDE_SYMBOL_DMS As Char = ">"c
    Public Const N_TERMINAL_PROTEIN_SYMBOL_DMS As Char = "["c
    Public Const C_TERMINAL_PROTEIN_SYMBOL_DMS As Char = "]"c
    Public Const PROTEIN_TERMINUS_SYMBOL_PHRP As Char = "-"c

    Public Const MSGF_RESULT_COLUMN_SpectrumFile As String = "#SpectrumFile"
    Public Const MSGF_RESULT_COLUMN_Title As String = "Title"
    Public Const MSGF_RESULT_COLUMN_Annotation As String = "Annotation"

    Public Const XT_RESULT_TO_SEQ_MAP_SUFFIX As String = "_xt_ResultToSeqMap.txt"
    Public Const XT_SEQ_TO_PROTEIN_MAP_SUFFIX As String = "_xt_SeqToProteinMap.txt"

    Public Const DOT_RAW_EXTENSION As String = ".raw"
    Public Const DOT_MZXML_EXTENSION As String = ".mzXML"

    Public Const MSGF_RESULT_FILENAME_SUFFIX As String = "_MSGF.txt"
    Public Const SCAN_STATS_FILENAME_SUFFIX As String = "_ScanStats.txt"
    Public Const EXTENDED_SCAN_STATS_FILENAME_SUFFIX As String = "_ScanStatsEx.txt"

    ' This RegEx is used to extract parent ion m/z from a filter string that does not contain msx
    ' ${ParentMZ} will hold the last parent ion m/z found
    '  For example, 756.71 in FTMS + p NSI d Full ms3 850.70@cid35.00 756.71@cid35.00 [195.00-2000.00]
    Private Const PARENTION_ONLY_NONMSX_REGEX As String = "[Mm][Ss]\d*[^\[\r\n]* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*(\[[^\]\r\n]\])?"

    '  This RegEx is used to extract parent ion m/z from a filter string that does contain msx
    '  ${ParentMZ} will hold the first parent ion m/z found (the first parent ion m/z corresponds to the highest peak)
    '  For example, 636.04 in FTMS + p NSI Full msx ms2 636.04@hcd28.00 641.04@hcd28.00 654.05@hcd28.00 [88.00-1355.00]
    Private Const PARENTION_ONLY_MSX_REGEX As String = "[Mm][Ss]\d* (?<ParentMZ>[0-9.]+)@?[A-Za-z]*\d*\.?\d*[^\[\r\n]*(\[[^\]\r\n]+\])?"

    Public Enum ePeptideHitResultType
        Unknown = 0
        Sequest = 1
        XTandem = 2
        Inspect = 3
        MSGFDB = 4      ' Aka MSGF+
        MSAlign = 5
        MODa = 6
        MODPlus = 7
        MSPathFinder = 8
    End Enum

    Public Enum ePHRPReaderErrorCodes As Integer
        NoError = 0
        InvalidInputFilePath = 1
        InputFileFormatNotRecognized = 2
        RequiredInputFileNotFound = 3
        MissingRawOrMzXmlFile = 4
        MSGFProgramNotFound = 5
        UnspecifiedError = -1
    End Enum

#End Region

#Region "Module variables"
    Private mDatasetName As String
    Private mInputFilePath As String
    Private mInputFolderPath As String

    Private mSkipDuplicatePSMs As Boolean

    Private ReadOnly mStartupOptions As clsPHRPStartupOptions

    Private mEchoMessagesToConsole As Boolean

    Private mCanRead As Boolean
    Private mInitialized As Boolean
    Private mModSummaryFileLoaded As Boolean

    ''' <summary>
    ''' When set to true, calls to MoveNext will read the next data line, but will skip several additional processing steps for performance reasons
    ''' </summary>
    ''' <remarks>If the peptide is a peptide of interest, you must call FinalizeCurrentPSM after calling .MoveNext()</remarks>
    Private mFastReadMode As Boolean

    Private mSourceFile As StreamReader
    Private mSourceFileLineCount As Integer
    Private mSourceFileLinesRead As Integer

    Private mPHRPParser As clsPHRPParser
    Private ReadOnly mPeptideMassCalculator As clsPeptideMassCalculator

    ' This dictionary contains mod symbols as the key and modification definition as the values
    Private ReadOnly mDynamicMods As SortedDictionary(Of Char, clsModificationDefinition)

    ' This dictionary contains amino acid names as the key and the corresponding mod modification (or mod modifications) 
    Private ReadOnly mStaticMods As SortedDictionary(Of String, List(Of clsModificationDefinition))

    ' This dictionary tracks the MSGFSpecProb values for each entry in the source file
    ' The keys are Result_ID and the string is MSGFSpecProb (stored as string to preserve formatting)
    Private mMSGFCachedResults As Dictionary(Of Integer, String)

    ' This dictionary tracks scan stats values, in particular elution time
    'The keys are ScanNumber and values are clsScanStatsInfo objects
    Private mScanStats As Dictionary(Of Integer, clsScanStatsInfo)

    ' This dictionary tracks extended scan stats values, including parent ion mz (via MonoisotopicMZ)and collision mode
    'The keys are ScanNumber and values are clsScanStatsExInfo objects
    Private mScanStatsEx As Dictionary(Of Integer, clsScanStatsExInfo)

    Private mPSMCurrent As clsPSM
    Private mPSMCurrentFinalized As Boolean

    Private mExtendedScanStatsValid As Boolean
    Private ReadOnly mExtendedScanStatsInfo As clsScanStatsExInfo

    Private mHeaderLineParsed As Boolean
    Private mCachedLineAvailable As Boolean
    Private mCachedLine As String
    Private mCachedPSM As clsPSM

    Private ReadOnly mErrorMessages As List(Of String)
    Private ReadOnly mWarningMessages As List(Of String)

    Private mErrorMessage As String = String.Empty
    Private mLocalErrorCode As ePHRPReaderErrorCodes

    ''' <summary>
    ''' RegEx to extract parent ions from filter strings that do not have Full msx
    ''' </summary>
    ''' <remarks>Shared (aka static) only to speed up unit tests</remarks>
    Private Shared ReadOnly mParentIonMzMatchNonMsx As Regex = New Regex(PARENTION_ONLY_NONMSX_REGEX, RegexOptions.Compiled Or RegexOptions.IgnoreCase)

    ''' <summary>
    ''' RegEx to extract parent ions from filter strings that have Full msx
    ''' </summary>
    ''' <remarks>Shared (aka static) only to speed up unit tests</remarks>
    Private Shared ReadOnly mParentIonMzMatchMsx As Regex = New Regex(PARENTION_ONLY_MSX_REGEX, RegexOptions.Compiled Or RegexOptions.IgnoreCase)

#End Region

#Region "Events"
    Public Event MessageEvent(strMessage As String)
    Public Event ErrorEvent(strErrorMessage As String)
    Public Event WarningEvent(strWarningMessage As String)
#End Region

#Region "Properties"

    ''' <summary>
    ''' Returns True if the input file was successfully opened and data remains to be read
    ''' </summary>
    ''' <value></value>
    ''' <returns>True if the file is readable</returns>
    ''' <remarks></remarks>
    Public ReadOnly Property CanRead() As Boolean
        Get
            Return mCanRead
        End Get
    End Property

    ''' <summary>
    ''' Returns the most recently loaded PSM
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property CurrentPSM As clsPSM
        Get
            Return mPSMCurrent
        End Get
    End Property

    ''' <summary>
    ''' Returns the most recently loaded PSM's sequence info (if available)
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property CurrentPSMSeqInfo As clsSeqInfo
        Get
            If mPSMCurrent Is Nothing OrElse mPHRPParser.SeqInfo Is Nothing Then
                Return Nothing
            Else
                Dim oSeqInfo As clsSeqInfo = Nothing
                mPHRPParser.SeqInfo.TryGetValue(mPSMCurrent.SeqID, oSeqInfo)
                Return oSeqInfo
            End If
        End Get
    End Property

    ''' <summary>
    ''' Dataset name (auto-determined based on the input filename)
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property DatasetName As String
        Get
            Return mDatasetName
        End Get
    End Property

    ''' <summary>
    ''' If True, then will display messages at the console
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property EchoMessagesToConsole As Boolean
        Get
            Return mEchoMessagesToConsole
        End Get
        Set(value As Boolean)
            mEchoMessagesToConsole = value
        End Set
    End Property

    ''' <summary>
    ''' Cached error messages
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ErrorMessages() As List(Of String)
        Get
            Return mErrorMessages
        End Get
    End Property

    ''' <summary>
    ''' Current error message
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ErrorMessage() As String
        Get
            Return mErrorMessage
        End Get
    End Property

    ''' <summary>
    ''' Used to enable fast read mode when calling MoveNext
    ''' When FastReadMode is True, you should call FinalizeCurrentPSM after calling MoveNext to populate the remaining fields if the peptide is a peptide of interest
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Once FastReadMode is enabled it cannot be turned off (this is a safety measure due to how data is cached)</remarks>
    Public Property FastReadMode() As Boolean
        Get
            Return mFastReadMode
        End Get
        Set(value As Boolean)
            If value Then
                mFastReadMode = True
            End If
        End Set
    End Property

    ''' <summary>
    ''' If True, then looks for and loads the modification definitions from the _ModSummary.txt file associated with the input file
    ''' Also reads the SeqInfo and related files
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property LoadModsAndSeqInfo() As Boolean
        Get
            Return mStartupOptions.LoadModsAndSeqInfo
        End Get
    End Property

    ''' <summary>
    ''' If true, then loads the MSGF SpecProb values from the _MSGF.txt file associated with the input file
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property LoadMSGFResults() As Boolean
        Get
            Return mStartupOptions.LoadMSGFResults
        End Get
    End Property

    ''' <summary>
    ''' If True, then loads the MASIC _ScanStats.txt file
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property LoadScanStatsData() As Boolean
        Get
            Return mStartupOptions.LoadScanStatsData
        End Get
    End Property

    ''' <summary>
    ''' The maximum number of proteins that will be tracked for each PSM
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property MaxProteinsPerPSM() As Integer
        Get
            Return mStartupOptions.MaxProteinsPerPSM
        End Get
        Set(value As Integer)
            mStartupOptions.MaxProteinsPerPSM = value
        End Set
    End Property

    ''' <summary>
    ''' Returns True if the ModSummary file was successfully loaded
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ModSummaryFileLoaded As Boolean
        Get
            Return mModSummaryFileLoaded
        End Get
    End Property

    ''' <summary>
    ''' Peptide hit result type; Sequest, XTandem, Inspect, or MSGFDB
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property PeptideHitResultType As ePeptideHitResultType
        Get
            If mPHRPParser Is Nothing Then
                Return ePeptideHitResultType.Unknown
            Else
                Return mPHRPParser.PeptideHitResultType
            End If
        End Get
    End Property

    ''' <summary>
    ''' Returns a number between 0 and 100 indicating the percentage of the source file that has been read
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property PercentComplete() As Single
        Get
            If mSourceFileLineCount > 0 Then
                Return mSourceFileLinesRead / CSng(mSourceFileLineCount) * 100.0!
            Else
                Return 0
            End If

        End Get
    End Property

    ''' <summary>
    ''' Returns the PHRP Parser object
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property PHRPParser() As clsPHRPParser
        Get
            Return mPHRPParser
        End Get
    End Property

    ''' <summary>
    ''' Returns the cached mapping between ResultID and SeqID
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ResultToSeqMap() As SortedList(Of Integer, Integer)
        Get
            If mPHRPParser Is Nothing Then
                Return New SortedList(Of Integer, Integer)
            Else
                Return mPHRPParser.ResultToSeqMap
            End If

        End Get
    End Property

    ''' <summary>
    ''' Returns the cached sequence info, where key is SeqID
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property SeqInfo() As SortedList(Of Integer, clsSeqInfo)
        Get
            If mPHRPParser Is Nothing Then
                Return New SortedList(Of Integer, clsSeqInfo)
            Else
                Return mPHRPParser.SeqInfo
            End If
        End Get
    End Property

    ''' <summary>
    ''' Returns the cached sequence to protein map information
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property SeqToProteinMap() As SortedList(Of Integer, List(Of clsProteinInfo))
        Get
            If mPHRPParser Is Nothing Then
                Return New SortedList(Of Integer, List(Of clsProteinInfo))
            Else
                Return mPHRPParser.SeqToProteinMap
            End If
        End Get
    End Property

    ''' <summary>
    ''' When True, then skips near-duplicate lines in the PHRP data file (lines with the same peptide in the same scan, but different protein names)
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property SkipDuplicatePSMs() As Boolean
        Get
            Return mSkipDuplicatePSMs
        End Get
        Set(value As Boolean)
            mSkipDuplicatePSMs = value
        End Set
    End Property

    ''' <summary>
    ''' Cached warning messages
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property WarningMessages() As List(Of String)
        Get
            Return mWarningMessages
        End Get
    End Property

#End Region

    ''' <summary>
    ''' Constructor that auto-determines the PeptideHit result type based on the filename
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
    Public Sub New(strInputFilePath As String)
        Me.New(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo:=True, blnLoadMSGFResults:=True, blnLoadScanStats:=False)
    End Sub

    ''' <summary>
    ''' Constructor where the PeptideHit result type is explicitly set
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="eResultType">Source file PeptideHit result type</param>
    ''' <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
    Public Sub New(strInputFilePath As String, eResultType As ePeptideHitResultType)
        Me.New(strInputFilePath, eResultType, blnLoadModsAndSeqInfo:=True, blnLoadMSGFResults:=True, blnLoadScanStats:=False)
    End Sub

    ''' <summary>
    ''' Constructor that auto-determines the PeptideHit result type based on the filename
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
    ''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
    ''' <remarks></remarks>
    Public Sub New(strInputFilePath As String, blnLoadModsAndSeqInfo As Boolean, blnLoadMSGFResults As Boolean)
        Me.New(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats:=False)
    End Sub

    ''' <summary>
    ''' Constructor that auto-determines the PeptideHit result type based on the filename
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
    ''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
    ''' <param name="blnLoadScanStats">If True, then looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
    ''' <remarks></remarks>
    Public Sub New(strInputFilePath As String, blnLoadModsAndSeqInfo As Boolean, blnLoadMSGFResults As Boolean, blnLoadScanStats As Boolean)
        Me.New(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats)
    End Sub

    ''' <summary>
    ''' Constructor that auto-determines the PeptideHit result type based on the filename
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="oStartupOptions">Startup options</param>
    ''' <remarks></remarks>
    Public Sub New(strInputFilePath As String, oStartupOptions As clsPHRPStartupOptions)
        Me.New(strInputFilePath, ePeptideHitResultType.Unknown, oStartupOptions)
    End Sub

    ''' <summary>
    ''' Constructor where the PeptideHit result type is explicitly set
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' ''' <param name="eResultType">Source file PeptideHit result type</param>
    ''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
    ''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
    ''' <remarks></remarks>
    Public Sub New(strInputFilePath As String, eResultType As ePeptideHitResultType, blnLoadModsAndSeqInfo As Boolean, blnLoadMSGFResults As Boolean)
        Me.New(strInputFilePath, eResultType, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats:=False)
    End Sub

    ''' <summary>
    ''' Constructor where the PeptideHit result type is explicitly set
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="eResultType">Source file PeptideHit result type</param>
    ''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
    ''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
    ''' <param name="blnLoadScanStats">If True, then looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
    ''' <remarks></remarks>
    Public Sub New(strInputFilePath As String, eResultType As ePeptideHitResultType, blnLoadModsAndSeqInfo As Boolean, blnLoadMSGFResults As Boolean, blnLoadScanStats As Boolean)

        Dim oStartupOptions = New clsPHRPStartupOptions()
        With oStartupOptions
            .LoadModsAndSeqInfo = blnLoadModsAndSeqInfo
            .LoadMSGFResults = blnLoadMSGFResults
            .LoadScanStatsData = blnLoadScanStats
        End With
        mStartupOptions = oStartupOptions

        mMSGFCachedResults = New Dictionary(Of Integer, String)

        mDynamicMods = New SortedDictionary(Of Char, clsModificationDefinition)
        mStaticMods = New SortedDictionary(Of String, List(Of clsModificationDefinition))

        mPeptideMassCalculator = New clsPeptideMassCalculator()

        mErrorMessages = New List(Of String)
        mWarningMessages = New List(Of String)

        InitializeClass(strInputFilePath, eResultType, mStartupOptions)

    End Sub

    ''' <summary>
    ''' Constructor where the PeptideHit result type is explicitly set
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="eResultType">Source file PeptideHit result type</param>
    ''' <param name="oStartupOptions">Startup options</param>
    ''' <remarks></remarks>
    Public Sub New(strInputFilePath As String, eResultType As ePeptideHitResultType, oStartupOptions As clsPHRPStartupOptions)

        If oStartupOptions Is Nothing Then
            Throw New ArgumentNullException(NameOf(oStartupOptions))
        End If

        mStartupOptions = oStartupOptions

        mMSGFCachedResults = New Dictionary(Of Integer, String)

        mDynamicMods = New SortedDictionary(Of Char, clsModificationDefinition)
        mStaticMods = New SortedDictionary(Of String, List(Of clsModificationDefinition))

        If oStartupOptions.PeptideMassCalculator Is Nothing Then
            mPeptideMassCalculator = New clsPeptideMassCalculator()
        Else
            mPeptideMassCalculator = oStartupOptions.PeptideMassCalculator
        End If

        mErrorMessages = New List(Of String)
        mWarningMessages = New List(Of String)

        InitializeClass(strInputFilePath, eResultType, mStartupOptions)

    End Sub

    ''' <summary>
    ''' Updates strFilePath to have _fht instead of _syn if strFilePath contains_syn yet strBasePHRPFileName contains _fht
    ''' </summary>
    ''' <param name="strFilePath"></param>
    ''' <param name="strBasePHRPFileName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function AutoSwitchToFHTIfRequired(strFilePath As String, strBasePHRPFileName As String) As String

        Dim fiFileInfo As FileInfo
        Dim intSynIndex As Integer
        Dim strFilePathFHT As String

        If String.IsNullOrEmpty(strBasePHRPFileName) Then
            Return strFilePath
        End If

        fiFileInfo = New FileInfo(strBasePHRPFileName)
        If fiFileInfo.Name.ToLower().IndexOf("_fht", StringComparison.Ordinal) > 0 Then
            ' strBasePHRPFileName is first-hits-file based

            fiFileInfo = New FileInfo(strFilePath)
            intSynIndex = fiFileInfo.Name.ToLower().IndexOf("_syn", StringComparison.Ordinal)
            If intSynIndex > 0 Then
                ' strFilePath is synopsis-file based
                ' Change strFilePath to contain _fht instead of _syn

                strFilePathFHT = fiFileInfo.Name.Substring(0, intSynIndex) & "_fht" & fiFileInfo.Name.Substring(intSynIndex + 4)

                If Path.IsPathRooted(strFilePath) Then
                    Return Path.Combine(fiFileInfo.DirectoryName, strFilePathFHT)
                Else
                    Return strFilePathFHT
                End If

            End If

        End If

        Return strFilePath
    End Function

    ''' <summary>
    ''' Clear any cached error messages
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearErrors()
        mErrorMessages.Clear()
        If Not mPHRPParser Is Nothing Then
            mPHRPParser.ClearErrors()
        End If
    End Sub

    Private Function CountLines(strTextFilePath As String) As Integer

        Dim intTotalLines As Integer

        Try
            intTotalLines = 0
            Using srReader = New StreamReader(New FileStream(strTextFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                While Not srReader.EndOfStream
                    srReader.ReadLine()
                    intTotalLines += 1
                End While
            End Using

        Catch ex As Exception
            Throw New Exception("Error counting the lines in " & Path.GetFileName(strTextFilePath) & ": " & ex.Message, ex)
        End Try

        Return intTotalLines

    End Function

    ''' <summary>
    ''' Clear any cached warning messages
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearWarnings()
        mWarningMessages.Clear()
        If Not mPHRPParser Is Nothing Then
            mPHRPParser.ClearWarnings()
        End If
    End Sub

    ''' <summary>
    ''' Initialize the class
    ''' </summary>
    ''' <param name="strInputFilePath">Input file to read</param>
    ''' <param name="eResultType">Source file PeptideHit result type</param>
    ''' <param name="oStartupOptions">Startup options</param>
    ''' <remarks></remarks>
    Private Sub InitializeClass(strInputFilePath As String, eResultType As ePeptideHitResultType, oStartupOptions As clsPHRPStartupOptions)

        mInitialized = False

        InitializeMemberVariables(oStartupOptions)

        InitializeReader(strInputFilePath, eResultType)

        mInitialized = True
    End Sub

    Private Sub InitializeMemberVariables(oStartupOptions As clsPHRPStartupOptions)
        mDatasetName = String.Empty
        mInputFilePath = String.Empty
        mInputFolderPath = String.Empty

        mCanRead = False
        mModSummaryFileLoaded = False

        mSkipDuplicatePSMs = True

        mEchoMessagesToConsole = False

        mErrorMessage = String.Empty
        mLocalErrorCode = ePHRPReaderErrorCodes.NoError

        mSourceFileLineCount = 0

    End Sub

    Private Function InitializeReader(strInputFilePath As String, eResultType As ePeptideHitResultType) As Boolean

        Dim strModSummaryFilePath As String = String.Empty

        Dim blnSuccess As Boolean

        Try
            If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
                ReportError("Input file name is empty")
                SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
                If Not mInitialized Then Throw New FileNotFoundException(mErrorMessage)
                Return False
            Else
                ' Confirm that the source file exists
                ' Make sure strInputFilePath points to a valid file
                Dim fiFileInfo = New FileInfo(strInputFilePath)

                mInputFolderPath = fiFileInfo.DirectoryName
                mInputFilePath = fiFileInfo.FullName

                If Not fiFileInfo.Exists Then
                    ReportError("Input file not found: " & strInputFilePath)
                    SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
                    If Not mInitialized Then Throw New FileNotFoundException(mErrorMessage)
                    Return False
                Else

                    ' Note that the following populates mDatasetName
                    blnSuccess = ValidateInputFiles(strInputFilePath, eResultType, strModSummaryFilePath)
                    If Not blnSuccess Then
                        SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound, True)
                        If Not mInitialized Then Throw New FileNotFoundException(mErrorMessage)
                        Return False
                    End If

                    ' Open the input file for reading
                    ' Note that this will also load the MSGFSpecProb info and ScanStats info
                    blnSuccess = InitializeParser(eResultType)

                    If blnSuccess AndAlso mStartupOptions.LoadModsAndSeqInfo Then
                        ' Read the PHRP Mod Summary File to populate mDynamicMods and mStaticMods
                        ' Note that the PHRPParser also loads the ModSummary file, and that mDynamicMods and mStaticMods are only used if the _SeqInfo.txt file is not found
                        blnSuccess = ReadModSummaryFile(strModSummaryFilePath, mDynamicMods, mStaticMods)
                        If Not blnSuccess Then
                            mModSummaryFileLoaded = False
                            blnSuccess = True
                        Else
                            mModSummaryFileLoaded = True
                        End If
                    End If

                    If blnSuccess AndAlso mStartupOptions.LoadMSGFResults Then
                        ' Cache the MSGF values (if present)
                        ReadAndCacheMSGFData()
                    End If

                    If blnSuccess AndAlso mStartupOptions.LoadScanStatsData Then
                        ' Cache the Scan Stats values (if present)
                        ReadScanStatsData()
                        ReadExtendedScanStatsData()
                    End If
                End If

            End If

        Catch ex As Exception
            HandleException("Error in InitializeReader", ex)
            If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
        End Try

        Return blnSuccess

    End Function

    Private Function InitializeParser(eResultType As ePeptideHitResultType) As Boolean

        Dim blnSuccess = True
        Dim strDatasetName As String = String.Copy(mDatasetName)

        Try
            If String.IsNullOrEmpty(strDatasetName) Then
                If mStartupOptions.LoadModsAndSeqInfo Then
                    ReportError("Dataset name is undefined; unable to continue since loading ModsAndSeqInfo")
                    Return False
                ElseIf mStartupOptions.LoadMSGFResults Then
                    ReportError("Dataset name is undefined; unable to continue since loading MSGF results")
                    Return False
                ElseIf mStartupOptions.LoadScanStatsData Then
                    ReportError("Dataset name is undefined; unable to continue since loading ScanStatsData")
                    Return False
                Else
                    strDatasetName = "Unknown_Dataset"
                End If
            End If

            ' Initialize some tracking variables
            mMSGFCachedResults.Clear()

            mPSMCurrent = New clsPSM()

            mSourceFileLinesRead = 0
            mHeaderLineParsed = False
            mCachedLineAvailable = False
            mCachedLine = String.Empty

            ' Open the peptide-hit result file (from PHRP) for reading
            ' Instantiate the appropriare PHRP Parser
            Select Case eResultType
                Case ePeptideHitResultType.Sequest
                    mPHRPParser = New clsPHRPParserSequest(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.XTandem
                    ' Note that Result to Protein mapping will be auto-loaded during instantiation of mPHRPParser
                    mPHRPParser = New clsPHRPParserXTandem(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.Inspect
                    mPHRPParser = New clsPHRPParserInspect(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.MSGFDB
                    mPHRPParser = New clsPHRPParserMSGFDB(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.MSAlign
                    mPHRPParser = New clsPHRPParserMSAlign(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.MODa
                    mPHRPParser = New clsPHRPParserMODa(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.MODPlus
                    mPHRPParser = New clsPHRPParserMODPlus(strDatasetName, mInputFilePath, mStartupOptions)

                Case ePeptideHitResultType.MSPathFinder
                    mPHRPParser = New clsPHRPParserMSPathFinder(strDatasetName, mInputFilePath, mStartupOptions)

                Case Else
                    'Should never get here; invalid result type specified
                    ReportError("Invalid PeptideHit ResultType specified: " & eResultType)
                    blnSuccess = False
            End Select

            If blnSuccess Then

                ' Attach the event handlers
                AddHandler mPHRPParser.ErrorEvent, AddressOf mPHRPParser_ErrorEvent
                AddHandler mPHRPParser.MessageEvent, AddressOf mPHRPParser_MessageEvent
                AddHandler mPHRPParser.WarningEvent, AddressOf mPHRPParser_WarningEvent

                ' Report any errors cached during instantiation of mPHRPParser
                For Each strMessage In mPHRPParser.ErrorMessages
                    ReportError(strMessage)
                Next

                ' Report any warnings cached during instantiation of mPHRPParser
                For Each strMessage In mPHRPParser.WarningMessages
                    ReportWarning(strMessage)
                Next

                mPHRPParser.ClearErrors()
                mPHRPParser.ClearWarnings()

                ' Open the data file and count the number of lines so that we can compute progress
                mSourceFileLineCount = CountLines(mInputFilePath)

                ' Open the data file for reading
                mSourceFile = New StreamReader(New FileStream(mInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                mCanRead = True

            End If

        Catch ex As Exception
            HandleException("Error in InitializeParser", ex)
            If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified folder
    ''' If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
    ''' </summary>
    ''' <param name="strInputFolderPath">Input folder path</param>
    ''' <returns>The full path to the most appropriate Synopsis or First hits file</returns>
    ''' <remarks></remarks>
    Public Shared Function AutoDetermineBestInputFile(strInputFolderPath As String) As String
        Dim eMatchedResultType = ePeptideHitResultType.Unknown
        Return AutoDetermineBestInputFile(strInputFolderPath, eMatchedResultType)
    End Function

    ''' <summary>
    ''' Looks for a valid _syn.txt or _fht.txt file for any dataset in the specified folder
    ''' If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
    ''' </summary>
    ''' <param name="strInputFolderPath">Input folder path</param>
    ''' <param name="eMatchedResultType">Output parameter: the result type of the best result file found</param>
    ''' <returns>The full path to the most appropriate Synopsis or First hits file</returns>
    ''' <remarks></remarks>
    Public Shared Function AutoDetermineBestInputFile(
       strInputFolderPath As String,
       <Out()> ByRef eMatchedResultType As ePeptideHitResultType) As String

        ' Find candidate dataset names in strInputFolderPath

        Dim fiInputFolder As DirectoryInfo
        Dim lstDatasetNames = New SortedSet(Of String)(StringComparer.CurrentCultureIgnoreCase)
        Dim lstFileSpec = New List(Of String)

        Dim strDataset As String
        Dim intCharIndex As Integer

        If String.IsNullOrWhiteSpace(strInputFolderPath) Then
            Throw New DirectoryNotFoundException("Input folder path is empty")
        End If

        fiInputFolder = New DirectoryInfo(strInputFolderPath)
        If Not fiInputFolder.Exists Then
            Throw New DirectoryNotFoundException("Input folder not found: " & strInputFolderPath)
        End If

        ' MSGF+
        lstFileSpec.Add(clsPHRPParserMSGFDB.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserMSGFDB.FILENAME_SUFFIX_FHT)

        ' X!Tandem (only has _xt.txt files)
        lstFileSpec.Add(clsPHRPParserXTandem.FILENAME_SUFFIX_SYN)

        ' MSAlign
        lstFileSpec.Add(clsPHRPParserMSAlign.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserMSAlign.FILENAME_SUFFIX_FHT)

        ' Inspect
        lstFileSpec.Add(clsPHRPParserInspect.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserInspect.FILENAME_SUFFIX_FHT)

        ' MODa
        lstFileSpec.Add(clsPHRPParserMODa.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserMODa.FILENAME_SUFFIX_FHT)

        ' MODPlus
        lstFileSpec.Add(clsPHRPParserMODPlus.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserMODPlus.FILENAME_SUFFIX_FHT)

        ' MSPathFinder
        lstFileSpec.Add(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_FHT)

        ' *****************
        ' ** Important: Sequest needs to be added last since files simply end in _syn.txt or _fht.txt)
        ' *****************
        ' Sequest
        lstFileSpec.Add(clsPHRPParserSequest.FILENAME_SUFFIX_SYN)
        lstFileSpec.Add(clsPHRPParserSequest.FILENAME_SUFFIX_FHT)


        For Each strFileSpec In lstFileSpec

            For Each fiFile In fiInputFolder.GetFiles("*" & strFileSpec)
                strDataset = fiFile.Name

                intCharIndex = strDataset.ToLower().IndexOf(strFileSpec, StringComparison.Ordinal)
                If intCharIndex > 0 Then
                    strDataset = strDataset.Substring(0, intCharIndex)

                    If Not lstDatasetNames.Contains(strDataset) Then
                        lstDatasetNames.Add(strDataset)
                    End If
                End If
            Next

        Next

        Return AutoDetermineBestInputFile(strInputFolderPath, lstDatasetNames.ToList(), eMatchedResultType)

    End Function

    ''' <summary>
    ''' Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified folder
    ''' If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
    ''' </summary>
    ''' <param name="strInputFolderPath">Input folder path</param>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <returns>The full path to the most appropriate Synopsis or First hits file</returns>
    ''' <remarks></remarks>
    Public Shared Function AutoDetermineBestInputFile(strInputFolderPath As String, strDatasetName As String) As String
        Dim eMatchedResultType = ePeptideHitResultType.Unknown
        Return AutoDetermineBestInputFile(strInputFolderPath, strDatasetName, eMatchedResultType)
    End Function

    ''' <summary>
    ''' Looks for a valid _syn.txt or _fht.txt file for the specified dataset in the specified folder
    ''' If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
    ''' </summary>
    ''' <param name="strInputFolderPath">Input folder path</param>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="eMatchedResultType">Output parameter: the result type of the best result file found</param>
    ''' <returns>The full path to the most appropriate Synopsis or First hits file</returns>
    ''' <remarks></remarks>
    Public Shared Function AutoDetermineBestInputFile(strInputFolderPath As String, strDatasetName As String, <Out()> ByRef eMatchedResultType As ePeptideHitResultType) As String
        Dim lstDatasetNames = New List(Of String)
        lstDatasetNames.Add(strDatasetName)

        Return AutoDetermineBestInputFile(strInputFolderPath, lstDatasetNames, eMatchedResultType)
    End Function

    ''' <summary>
    ''' Looks for a valid _syn.txt or _fht.txt file for the given list of datasets in the specified folder
    ''' If both the _syn.txt and _fht.txt files are present, then chooses the file with _ResultToSeqMap.txt and _SeqInfo.txt files
    ''' </summary>
    ''' <param name="strInputFolderPath">Input folder path</param>
    ''' <param name="lstDatasetNames">List of dataset names to search for</param>
    ''' <param name="eMatchedResultType">Output parameter: the result type of the best result file found</param>
    ''' <returns>The full path to the most appropriate Synopsis or First hits file</returns>
    ''' <remarks></remarks>
    Public Shared Function AutoDetermineBestInputFile(
       strInputFolderPath As String,
       lstDatasetNames As List(Of String),
       <Out()> ByRef eMatchedResultType As ePeptideHitResultType) As String

        Dim fiInputFolder As DirectoryInfo

        ' Items in this list are KeyValuePairs where the key is a filename to look for and the value is a PeptideHitResultType
        Dim lstFilesToFind As List(Of KeyValuePair(Of String, ePeptideHitResultType))

        ' This list contains the standard PHRP file suffixes
        Dim lstAuxiliaryFileSuffixes As List(Of String) = GetPHRPAuxiliaryFileSuffixes()

        ' The key in this variable is the full path to the best Synopsis or First hits file and the value is the number of PHRP-related auxiliary files
        Dim kvBestSynOrFHTFile As KeyValuePair(Of String, Integer)
        kvBestSynOrFHTFile = New KeyValuePair(Of String, Integer)(String.Empty, 0)

        ' Set the matched result type to Unknown for now
        eMatchedResultType = ePeptideHitResultType.Unknown

        If String.IsNullOrWhiteSpace(strInputFolderPath) Then
            Throw New DirectoryNotFoundException("Input folder path is empty")
        End If

        fiInputFolder = New DirectoryInfo(strInputFolderPath)
        If Not fiInputFolder.Exists Then
            Throw New DirectoryNotFoundException("Input folder not found: " & strInputFolderPath)
        End If

        If lstDatasetNames Is Nothing OrElse lstDatasetNames.Count = 0 Then
            Throw New ArgumentException("List lstDatasetNames cannot be empty; cannot determine the best input file")
        End If

        ' Construct a list of the files to search for
        lstFilesToFind = New List(Of KeyValuePair(Of String, ePeptideHitResultType))

        For Each strDataset As String In lstDatasetNames
            ' MSGF+
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMSGFDB.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MSGFDB))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMSGFDB.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MSGFDB))

            ' X!Tandem
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserXTandem.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.XTandem))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserXTandem.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.XTandem))

            ' MSAlign
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMSAlign.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MSAlign))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMSAlign.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MSAlign))

            ' MODa
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMODa.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MODa))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMODa.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MODa))

            ' MODPlus
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMODPlus.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MODPlus))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMODPlus.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MODPlus))

            ' MSPathFinder
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMSPathFinder.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.MSPathFinder))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserMSPathFinder.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.MSPathFinder))

            ' Inspect
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserInspect.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.Inspect))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserInspect.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.Inspect))

            ' Sequest (needs to be added last since files simply end in _syn.txt or _fht.txt)
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserSequest.GetPHRPSynopsisFileName(strDataset), ePeptideHitResultType.Sequest))
            lstFilesToFind.Add(New KeyValuePair(Of String, ePeptideHitResultType)(clsPHRPParserSequest.GetPHRPFirstHitsFileName(strDataset), ePeptideHitResultType.Sequest))
        Next

        For Each kvFileToFind As KeyValuePair(Of String, ePeptideHitResultType) In lstFilesToFind
            If Not String.IsNullOrEmpty(kvFileToFind.Key) Then
                Dim fiSynOrFHTFile = New FileInfo(Path.Combine(fiInputFolder.FullName, kvFileToFind.Key))

                If fiSynOrFHTFile.Exists Then
                    ' Match found
                    ' Look for PHRP-related auxiliary files

                    Dim intAuxFileCount = 0
                    Dim strBaseName As String
                    strBaseName = Path.Combine(fiSynOrFHTFile.Directory.FullName, Path.GetFileNameWithoutExtension(fiSynOrFHTFile.Name))

                    For Each strSuffix As String In lstAuxiliaryFileSuffixes
                        If File.Exists(strBaseName & strSuffix) Then
                            intAuxFileCount += 1
                        End If
                    Next

                    If String.IsNullOrEmpty(kvBestSynOrFHTFile.Key) OrElse intAuxFileCount > kvBestSynOrFHTFile.Value Then
                        kvBestSynOrFHTFile = New KeyValuePair(Of String, Integer)(fiSynOrFHTFile.FullName, intAuxFileCount)
                        eMatchedResultType = kvFileToFind.Value
                    End If
                End If
            End If
        Next

        If String.IsNullOrWhiteSpace(kvBestSynOrFHTFile.Key) Then
            If lstDatasetNames.Count = 1 Then
                Console.WriteLine("Could not find a Synopsis or First Hits file for dataset " & lstDatasetNames.First())
            Else
                Console.WriteLine("Could not find a Synopsis or First Hits file for any of the candidate datasets")
            End If

            Console.WriteLine("Looked for the following files:")
            For Each fileName In lstFilesToFind
                If Not String.IsNullOrWhiteSpace(fileName.Key) Then
                    Console.WriteLine("  " & fileName.Key)
                End If
            Next
        End If

        ' kvBestSynOrFHTFile should now contain the PHRP result file with the most auxiliary files
        Return kvBestSynOrFHTFile.Key

    End Function

    ''' <summary>
    ''' Auto-determine the dataset name using the input file path
    ''' </summary>
    ''' <param name="strFilePath"></param>
    ''' <returns>Dataset name</returns>
    ''' <remarks>Returns an empty string if unable to determine the dataset name</remarks>
    Public Shared Function AutoDetermineDatasetName(strFilePath As String) As String
        Dim eResultType As ePeptideHitResultType

        eResultType = AutoDetermineResultType(strFilePath)
        Return AutoDetermineDatasetName(strFilePath, eResultType)

    End Function

    ''' <summary>
    ''' Auto-determine the dataset name using the input file path and specified PeptideHit result type
    ''' </summary>
    ''' <param name="strFilePath"></param>
    ''' <param name="eResultType"></param>
    ''' <returns>Dataset name</returns>
    ''' <remarks>Returns an empty string if unable to determine the dataset name</remarks>
    Public Shared Function AutoDetermineDatasetName(strFilePath As String, eResultType As ePeptideHitResultType) As String

        Dim strDatasetName As String = String.Empty
        Dim strInputFileName As String

        strInputFileName = Path.GetFileNameWithoutExtension(strFilePath)

        Select Case eResultType
            Case ePeptideHitResultType.Sequest,
                 ePeptideHitResultType.Inspect,
                 ePeptideHitResultType.MSGFDB,
                 ePeptideHitResultType.MSAlign,
                 ePeptideHitResultType.MODa,
                 ePeptideHitResultType.MODPlus,
                 ePeptideHitResultType.MSPathFinder

                If strInputFileName.ToLower.EndsWith("_fht") OrElse strInputFileName.ToLower.EndsWith("_syn") Then
                    strDatasetName = strInputFileName.Substring(0, strInputFileName.Length - 4)

                    If eResultType = ePeptideHitResultType.Inspect Then
                        If strDatasetName.ToLower.EndsWith("_inspect") Then
                            strDatasetName = strDatasetName.Substring(0, strDatasetName.Length - "_inspect".Length)
                        End If

                    ElseIf eResultType = ePeptideHitResultType.MSGFDB Then
                        If strDatasetName.ToLower.EndsWith("_msgfdb") Then
                            strDatasetName = strDatasetName.Substring(0, strDatasetName.Length - "_msgfdb".Length)
                        End If

                    ElseIf eResultType = ePeptideHitResultType.MSAlign Then
                        If strDatasetName.ToLower.EndsWith("_msalign") Then
                            strDatasetName = strDatasetName.Substring(0, strDatasetName.Length - "_msalign".Length)
                        End If

                    ElseIf eResultType = ePeptideHitResultType.MODa Then
                        If strDatasetName.ToLower.EndsWith("_moda") Then
                            strDatasetName = strDatasetName.Substring(0, strDatasetName.Length - "_moda".Length)
                        End If

                    ElseIf eResultType = ePeptideHitResultType.MODPlus Then
                        If strDatasetName.ToLower.EndsWith("_modp") Then
                            strDatasetName = strDatasetName.Substring(0, strDatasetName.Length - "_modp".Length)
                        End If

                    ElseIf eResultType = ePeptideHitResultType.MSPathFinder Then
                        If strDatasetName.ToLower.EndsWith("_mspath") Then
                            strDatasetName = strDatasetName.Substring(0, strDatasetName.Length - "_mspath".Length)
                        End If

                    End If
                End If

            Case ePeptideHitResultType.XTandem
                If strInputFileName.ToLower.EndsWith("_xt") Then
                    strDatasetName = strInputFileName.Substring(0, strInputFileName.Length - 3)
                End If

            Case Else
                ' Unknown format

        End Select

        If String.IsNullOrEmpty(strDatasetName) Then
            Dim strFilePathTrimmed As String = String.Empty
            If AutoTrimExtraSuffix(strFilePath, strFilePathTrimmed) Then
                strDatasetName = AutoDetermineDatasetName(strFilePathTrimmed, eResultType)
            End If
        End If

        Return strDatasetName

    End Function

    ''' <summary>
    ''' Determine the PeptideHit result type given the input file path
    ''' </summary>
    ''' <param name="strFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function AutoDetermineResultType(strFilePath As String) As ePeptideHitResultType

        Dim strFilePathLCase As String

        Dim strHeaderLine As String

        Dim eResultType As ePeptideHitResultType
        eResultType = ePeptideHitResultType.Unknown

        strFilePathLCase = strFilePath.ToLower()

        Try
            If strFilePathLCase.EndsWith(clsPHRPParserXTandem.FILENAME_SUFFIX_SYN) Then
                eResultType = ePeptideHitResultType.XTandem
            Else
                If strFilePathLCase.EndsWith(clsPHRPParserMSGFDB.FILENAME_SUFFIX_SYN) OrElse strFilePathLCase.EndsWith(clsPHRPParserMSGFDB.FILENAME_SUFFIX_FHT) Then
                    eResultType = ePeptideHitResultType.MSGFDB

                ElseIf strFilePathLCase.EndsWith(clsPHRPParserMSAlign.FILENAME_SUFFIX_SYN) OrElse strFilePathLCase.EndsWith(clsPHRPParserMSAlign.FILENAME_SUFFIX_FHT) Then
                    eResultType = ePeptideHitResultType.MSAlign

                ElseIf strFilePathLCase.EndsWith(clsPHRPParserMODa.FILENAME_SUFFIX_SYN) OrElse strFilePathLCase.EndsWith(clsPHRPParserMODa.FILENAME_SUFFIX_FHT) Then
                    eResultType = ePeptideHitResultType.MODa

                ElseIf strFilePathLCase.EndsWith(clsPHRPParserMODPlus.FILENAME_SUFFIX_SYN) OrElse strFilePathLCase.EndsWith(clsPHRPParserMODPlus.FILENAME_SUFFIX_FHT) Then
                    eResultType = ePeptideHitResultType.MODPlus

                ElseIf strFilePathLCase.EndsWith(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_SYN) OrElse strFilePathLCase.EndsWith(clsPHRPParserMSPathFinder.FILENAME_SUFFIX_FHT) Then
                    eResultType = ePeptideHitResultType.MSPathFinder

                ElseIf strFilePathLCase.EndsWith(clsPHRPParserInspect.FILENAME_SUFFIX_SYN) OrElse strFilePathLCase.EndsWith(clsPHRPParserInspect.FILENAME_SUFFIX_FHT) Then
                    eResultType = ePeptideHitResultType.Inspect

                Else
                    ' Open the file and read the header line to determine if this is a Sequest file, Inspect file, MSGFDB, or something else

                    If Not File.Exists(strFilePath) Then
                        ' File doesn't exist; assume Sequest
                        eResultType = ePeptideHitResultType.Sequest
                    Else

                        Using srInFile = New StreamReader(New FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                            If Not srInFile.EndOfStream Then
                                strHeaderLine = srInFile.ReadLine()

                                If LineContainsValues(strHeaderLine, clsPHRPParserInspect.DATA_COLUMN_MQScore, clsPHRPParserInspect.DATA_COLUMN_TotalPRMScore) Then

                                    eResultType = ePeptideHitResultType.Inspect

                                ElseIf _
                                  LineContainsValues(strHeaderLine, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb) OrElse
                                  LineContainsValues(strHeaderLine, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecEValue) OrElse
                                  LineContainsValues(strHeaderLine, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore) Then

                                    eResultType = ePeptideHitResultType.MSGFDB

                                ElseIf LineContainsValues(strHeaderLine, clsPHRPParserSequest.DATA_COLUMN_XCorr, clsPHRPParserSequest.DATA_COLUMN_DelCn) Then

                                    eResultType = ePeptideHitResultType.Sequest

                                End If
                            End If

                        End Using

                    End If

                End If
            End If

        Catch ex As Exception
            Throw
        End Try

        If eResultType = ePeptideHitResultType.Unknown Then
            Dim strFilePathTrimmed As String = String.Empty
            If AutoTrimExtraSuffix(strFilePath, strFilePathTrimmed) Then
                eResultType = AutoDetermineResultType(strFilePathTrimmed)
            End If
        End If

        Return eResultType

    End Function

    Private Shared Function AutoTrimExtraSuffix(strFilePath As String, ByRef strFilePathTrimmed As String) As Boolean

        ' Check whether strfilePathLCase ends in other known PHRP extensions
        Dim lstExtraSuffixes As List(Of String)
        lstExtraSuffixes = GetPHRPAuxiliaryFileSuffixes()

        For Each strSuffix As String In lstExtraSuffixes
            If strFilePath.ToLower().EndsWith(strSuffix.ToLower()) Then

                strFilePathTrimmed = strFilePath.Substring(0, strFilePath.Length - strSuffix.Length) & ".txt"
                Return True

            End If
        Next

        Return False

    End Function

    ''' <summary>
    ''' Look for dynamic mod symbols in the peptide sequence; replace with the corresponding mod masses
    ''' Note that if the _SeqInfo.txt file is available, then this function will not be used
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="strPeptideWithNumericMods">Peptide with numeric mods (output)</param>
    ''' <param name="lstPeptideMods">List of modified amino acids (output)</param>
    ''' <returns>True if success, false if an error</returns>
    ''' <remarks>strPeptideWithNumericMods will look like R.TDM+15.9949ESALPVTVLSAEDIAK.T</remarks>
    Private Function ConvertModsToNumericMods(
      strPeptide As String,
      <Out()> ByRef strPeptideWithNumericMods As String,
      <Out()> ByRef lstPeptideMods As List(Of clsAminoAcidModInfo)) As Boolean

        Static sbNewPeptide As New Text.StringBuilder

        Dim intPeptideLength As Integer
        Dim chMostRecentResidue As Char
        Dim intResidueLocInPeptide = 0
        Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

        Dim intIndex As Integer
        Dim intIndexStart As Integer
        Dim intIndexEnd As Integer

        lstPeptideMods = New List(Of clsAminoAcidModInfo)
        strPeptideWithNumericMods = String.Empty

        Try


            If mDynamicMods.Count = 0 AndAlso mStaticMods.Count = 0 Then
                ' No mods are defined; simply update strPeptideWithNumericMods to be strPeptide
                strPeptideWithNumericMods = strPeptide
                Return True
            End If

            sbNewPeptide.Length = 0
            intPeptideLength = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strPeptide, True).Length

            intIndexStart = 0
            intIndexEnd = strPeptide.Length - 1

            If strPeptide.Length >= 4 Then
                If strPeptide.Chars(1) = "." Then
                    ' Peptide is of the form R.HRDTGILDSIGR.F
                    ' Skip the first two characters
                    intIndexStart = 2
                End If

                If strPeptide.Chars(strPeptide.Length - 2) = "." Then
                    ' Peptide is of the form R.HRDTGILDSIGR.F
                    ' Skip the last two characters
                    intIndexEnd = strPeptide.Length - 3
                End If

            End If

            intIndex = 0
            chMostRecentResidue = "."c
            Do While intIndex < strPeptide.Length
                If intIndex < intIndexStart OrElse intIndex > intIndexEnd Then
                    ' We're before or after the primary peptide sequence; simply append the character
                    sbNewPeptide.Append(strPeptide.Chars(intIndex))
                Else
                    If IsLetterAtoZ(strPeptide.Chars(intIndex)) Then
                        chMostRecentResidue = strPeptide.Chars(intIndex)
                        intResidueLocInPeptide += 1
                        If intResidueLocInPeptide = 1 Then
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
                        ElseIf intResidueLocInPeptide = intPeptideLength Then
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
                        Else
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
                        End If

                        ' Character is a letter; append it
                        sbNewPeptide.Append(chMostRecentResidue)

                        If mStaticMods.Count > 0 Then

                            ' See if it is present in mStaticMods (this is a case-sensitive search)
                            AddStaticModIfPresent(mStaticMods, chMostRecentResidue, intResidueLocInPeptide, eResidueTerminusState, sbNewPeptide, lstPeptideMods)

                            If intIndex = intIndexStart AndAlso mStaticMods.Count > 0 Then
                                ' We're at the N-terminus of the peptide
                                ' Possibly add a static N-terminal peptide mod (for example, iTRAQ8, which is 304.2022 Da)
                                AddStaticModIfPresent(mStaticMods, N_TERMINAL_PEPTIDE_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus, sbNewPeptide, lstPeptideMods)

                                If strPeptide.StartsWith(PROTEIN_TERMINUS_SYMBOL_PHRP) Then
                                    ' We're at the N-terminus of the protein
                                    ' Possibly add a static N-terminal protein mod
                                    AddStaticModIfPresent(mStaticMods, N_TERMINAL_PROTEIN_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus, sbNewPeptide, lstPeptideMods)
                                End If
                            End If
                        End If

                    Else
                        ' Not a letter; see if it is present in mDynamicMods
                        AddDynamicModIfPresent(mDynamicMods, chMostRecentResidue, strPeptide.Chars(intIndex), intResidueLocInPeptide, eResidueTerminusState, sbNewPeptide, lstPeptideMods)
                    End If

                    If intIndex = intIndexEnd AndAlso mStaticMods.Count > 0 Then
                        ' Possibly add a static C-terminal peptide mod
                        AddStaticModIfPresent(mStaticMods, C_TERMINAL_PEPTIDE_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus, sbNewPeptide, lstPeptideMods)

                        If strPeptide.EndsWith(PROTEIN_TERMINUS_SYMBOL_PHRP) Then
                            ' We're at the C-terminus of the protein
                            ' Possibly add a static C-terminal protein mod
                            AddStaticModIfPresent(mStaticMods, C_TERMINAL_PROTEIN_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus, sbNewPeptide, lstPeptideMods)
                        End If

                    End If

                End If
                intIndex += 1
            Loop

            strPeptideWithNumericMods = sbNewPeptide.ToString

        Catch ex As Exception
            HandleException("Error adding dynamic and static mod masses to peptide " & strPeptide, ex)
            Return False
        End Try

        Return True

    End Function

    Private Sub AddDynamicModIfPresent(
      objMods As SortedDictionary(Of Char, clsModificationDefinition),
      chResidue As Char,
      chModSymbol As Char,
      ResidueLocInPeptide As Integer,
      ResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
      sbNewPeptide As Text.StringBuilder,
      lstPeptideMods As List(Of clsAminoAcidModInfo))

        Dim objModDef As clsModificationDefinition = Nothing

        If objMods.TryGetValue(chModSymbol, objModDef) Then
            ' Mod mass found for dynamic mod symbol; append the mod
            sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(objModDef.ModificationMass, 4))
            lstPeptideMods.Add(New clsAminoAcidModInfo(chResidue, ResidueLocInPeptide, ResidueTerminusState, objModDef))
        End If

    End Sub

    Private Sub AddStaticModIfPresent(
       objMods As SortedDictionary(Of String, List(Of clsModificationDefinition)),
       chResidue As Char,
       ResidueLocInPeptide As Integer,
       ResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
       sbNewPeptide As Text.StringBuilder,
       lstPeptideMods As List(Of clsAminoAcidModInfo))

        Dim lstModDefs As List(Of clsModificationDefinition) = Nothing

        If objMods.TryGetValue(chResidue, lstModDefs) Then
            ' Static mod applies to this residue; append the mod (add a plus sign if it doesn't start with a minus sign)

            For Each objModDef In lstModDefs
                sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(objModDef.ModificationMass, 4))
                lstPeptideMods.Add(New clsAminoAcidModInfo(chResidue, ResidueLocInPeptide, ResidueTerminusState, objModDef))
            Next

        End If

    End Sub

    ''' <summary>
    ''' Determines the collision mode using the Scan Type name
    ''' </summary>
    ''' <param name="strScanTypeName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function GetCollisionMode(strScanTypeName As String) As String
        ' Typical scan types, along with DMS usage count as of 4/26/2012

        ' Scan Type     Usage Count
        ' CID-HMSn        2,800 spectra
        ' CID-MSn        68,335
        ' CID-SRM       258,479
        ' ETD-HMSn          974
        ' ETD-MSn           736
        ' GC-MS           2,646
        ' HCD-HMSn        6,767
        ' HMS            73,431
        ' HMSn               82
        ' MRM_Full_NL         1
        ' MS             23,201
        ' MSn             5,229
        ' PQD-HMSn            1
        ' PQD-MSn            51
        ' Q1MS              216
        ' Q3MS              363
        ' SA_ETD-HMSn       368
        ' SA_ETD-MSn      2,925
        ' SIM ms            306
        ' SRM            95,207
        ' Zoom-MS            29


        Dim strCollisionMode As String

        strCollisionMode = strScanTypeName.ToUpper()
        If strCollisionMode.StartsWith("SA_") Then
            strCollisionMode = strCollisionMode.Substring(3)
        End If

        If strCollisionMode.StartsWith("CID") Then
            Return "CID"
        ElseIf strCollisionMode.StartsWith("ETD") Then
            Return "ETD"
        ElseIf strCollisionMode.StartsWith("HCD") Then
            Return "HCD"
        ElseIf strCollisionMode.StartsWith("PQD") Then
            Return "PQD"
        ElseIf strCollisionMode.StartsWith("SID") Then
            Return "SID"
        Else
            Return ""
        End If

    End Function

    ''' <summary>
    ''' Returns the filename of the MSGF file that corresponds to strSynopsisOrFirstHitsFileName
    ''' </summary>
    ''' <param name="strSynopsisOrFirstHitsFileName">Filename (or full path) to the synopsis or first-hits file</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetMSGFFileName(strSynopsisOrFirstHitsFileName As String) As String
        Return Path.GetFileNameWithoutExtension(strSynopsisOrFirstHitsFileName) & MSGF_RESULT_FILENAME_SUFFIX
    End Function

    ''' <summary>
    ''' Get the peptide hit result type for the given result type name
    ''' </summary>
    ''' <param name="ResultTypeName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPeptideHitResultType(ResultTypeName As String) As ePeptideHitResultType
        Select Case ResultTypeName.ToLower()
            Case "Peptide_Hit".ToLower
                Return ePeptideHitResultType.Sequest

            Case "XT_Peptide_Hit".ToLower
                Return ePeptideHitResultType.XTandem

            Case "IN_Peptide_Hit".ToLower
                Return ePeptideHitResultType.Inspect

            Case "MSG_Peptide_Hit".ToLower
                Return ePeptideHitResultType.MSGFDB

            Case "MSA_Peptide_Hit".ToLower
                Return ePeptideHitResultType.MSAlign

            Case "MODa_Peptide_Hit".ToLower
                Return ePeptideHitResultType.MODa

            Case "MODPlus_Peptide_Hit".ToLower
                Return ePeptideHitResultType.MODPlus

            Case "MSP_Peptide_Hit".ToLower
                Return ePeptideHitResultType.MSPathFinder

            Case Else
                Return ePeptideHitResultType.Unknown
        End Select
    End Function

    ''' <summary>
    ''' Get the list of auxiliary file suffixes
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPAuxiliaryFileSuffixes() As List(Of String)
        Dim lstAuxSuffixes As List(Of String)
        lstAuxSuffixes = New List(Of String)

        lstAuxSuffixes.Add("_ResultToSeqMap.txt")
        lstAuxSuffixes.Add("_SeqToProteinMap.txt")
        lstAuxSuffixes.Add("_SeqInfo.txt")
        lstAuxSuffixes.Add("_MSGF.txt")
        lstAuxSuffixes.Add("_ProteinMods.txt")
        lstAuxSuffixes.Add("_ModDetails.txt")
        lstAuxSuffixes.Add("_ModSummary.txt")

        Return lstAuxSuffixes

    End Function

    Private Shared Function GetPHRPFileFreeParser(eResultType As ePeptideHitResultType, strDatasetName As String) As clsPHRPParser

        Static oCachedParser As clsPHRPParser

        ' ReSharper disable once UseImplicitlyTypedVariableEvident
        Static eCachedResultType As ePeptideHitResultType = ePeptideHitResultType.Unknown

        Static strCachedDataset As String = String.Empty

        If eCachedResultType <> ePeptideHitResultType.Unknown AndAlso
           eCachedResultType = eResultType AndAlso
           strCachedDataset = strDatasetName Then

            Return oCachedParser
        End If

        Select Case eResultType
            Case ePeptideHitResultType.Sequest
                oCachedParser = New clsPHRPParserSequest(strDatasetName, String.Empty)

            Case ePeptideHitResultType.XTandem
                oCachedParser = New clsPHRPParserXTandem(strDatasetName, String.Empty)

            Case ePeptideHitResultType.Inspect
                oCachedParser = New clsPHRPParserInspect(strDatasetName, String.Empty)

            Case ePeptideHitResultType.MSGFDB
                oCachedParser = New clsPHRPParserMSGFDB(strDatasetName, String.Empty)

            Case ePeptideHitResultType.MSAlign
                oCachedParser = New clsPHRPParserMSAlign(strDatasetName, String.Empty)

            Case ePeptideHitResultType.MODa
                oCachedParser = New clsPHRPParserMODa(strDatasetName, String.Empty)

            Case ePeptideHitResultType.MODPlus
                oCachedParser = New clsPHRPParserMODPlus(strDatasetName, String.Empty)

            Case ePeptideHitResultType.MSPathFinder
                oCachedParser = New clsPHRPParserMSPathFinder(strDatasetName, String.Empty)

            Case Else
                Throw New Exception("Unsupported ePeptideHitResultType value: " & eResultType)
        End Select

        eCachedResultType = eResultType
        strCachedDataset = String.Copy(strDatasetName)

        Return oCachedParser

    End Function

    ''' <summary>
    ''' Returns the default first-hits file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPFirstHitsFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strPHRPResultsFileName = oParser.PHRPFirstHitsFileName()

        Return strPHRPResultsFileName

    End Function

    ''' <summary>
    ''' Returns the default ModSummary file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPModSummaryFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strPHRPModSummaryFileName = oParser.PHRPModSummaryFileName()

        Return strPHRPModSummaryFileName

    End Function

    ''' <summary>
    ''' Returns the default PepToProtMap file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPPepToProteinMapFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strPepToProteinMapFileName = oParser.PHRPPepToProteinMapFileName()

        Return strPepToProteinMapFileName

    End Function

    ''' <summary>
    ''' Returns the default ProteinMods file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPProteinModsFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strProteinModsFileName = oParser.PHRPProteinModsFileName()

        Return strProteinModsFileName

    End Function

    ''' <summary>
    ''' Returns the default Synopsis file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPSynopsisFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strPHRPResultsFileName = oParser.PHRPSynopsisFileName()

        Return strPHRPResultsFileName

    End Function

    ''' <summary>
    ''' Returns the default ResultToSeq Map file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPResultToSeqMapFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strResultToSeqMapFilename = oParser.PHRPResultToSeqMapFileName()

        Return strResultToSeqMapFilename

    End Function

    ''' <summary>
    ''' Returns the default SeqInfo file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPSeqInfoFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strSeqInfoFilename = oParser.PHRPSeqInfoFileName()

        Return strSeqInfoFilename

    End Function

    ''' <summary>
    ''' Returns the default SeqToProtein Map file name for the given PeptideHit result type
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetPHRPSeqToProteinMapFileName(eResultType As ePeptideHitResultType, strDatasetName As String) As String

        Dim oParser As clsPHRPParser = GetPHRPFileFreeParser(eResultType, strDatasetName)
        Dim strSeqToProteinMapFileName = oParser.PHRPSeqToProteinMapFileName()

        Return strSeqToProteinMapFileName

    End Function

    ''' <summary>
    ''' Get the ScanStats filename for the given dataset
    ''' </summary>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetScanStatsFilename(strDatasetName As String) As String
        Return strDatasetName & SCAN_STATS_FILENAME_SUFFIX
    End Function

    ''' <summary>
    ''' Get the extended ScanStats filename for the given dataset
    ''' </summary>
    ''' <param name="strDatasetName"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetExtendedScanStatsFilename(strDatasetName As String) As String
        Return strDatasetName & EXTENDED_SCAN_STATS_FILENAME_SUFFIX
    End Function

    ''' <summary>
    ''' Get the tool version info filename for the given analysis tool
    ''' </summary>
    ''' <param name="eResultType"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetToolVersionInfoFilename(eResultType As ePeptideHitResultType) As String

        Dim strToolVersionInfoFilename As String = String.Empty

        Select Case eResultType
            Case ePeptideHitResultType.Sequest
                strToolVersionInfoFilename = "Tool_Version_Info_Sequest.txt"

            Case ePeptideHitResultType.XTandem
                strToolVersionInfoFilename = "Tool_Version_Info_XTandem.txt"

            Case ePeptideHitResultType.Inspect
                strToolVersionInfoFilename = "Tool_Version_Info_Inspect.txt"

            Case ePeptideHitResultType.MSGFDB
                strToolVersionInfoFilename = "Tool_Version_Info_MSGFDB.txt"

            Case ePeptideHitResultType.MSAlign
                strToolVersionInfoFilename = "Tool_Version_Info_MSAlign.txt"

            Case ePeptideHitResultType.MODa
                strToolVersionInfoFilename = "Tool_Version_Info_MODa.txt"

            Case ePeptideHitResultType.MODPlus
                strToolVersionInfoFilename = "Tool_Version_Info_MODPlus.txt"

            Case ePeptideHitResultType.MSPathFinder
                strToolVersionInfoFilename = "Tool_Version_Info_MSPathFinder.txt"

        End Select

        Return strToolVersionInfoFilename
    End Function

    Private Sub HandleException(strBaseMessage As String, ex As Exception)
        If String.IsNullOrEmpty(strBaseMessage) Then
            strBaseMessage = "Error"
        End If

        ReportError(strBaseMessage & ": " & ex.Message)

    End Sub

    ''' <summary>
    ''' Returns true if the character is a letter between A and Z or a and z
    ''' </summary>
    ''' <param name="chChar">Character to examine</param>
    ''' <returns></returns>
    ''' <remarks>The Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters; use this function to only return True for normal letters between A and Z</remarks>
    Public Shared Function IsLetterAtoZ(chChar As Char) As Boolean

        Static reIsLetter As Regex = New Regex("[A-Za-z]", RegexOptions.Compiled)

        If reIsLetter.IsMatch(chChar) Then
            Return True
        Else
            Return False
        End If

    End Function

    ''' <summary>
    ''' Examines the string to determine if it is numeric
    ''' </summary>
    ''' <param name="strData"></param>
    ''' <returns>True if a number, otherwise false</returns>
    Public Shared Function IsNumber(strData As String) As Boolean
        Try
            If Double.TryParse(strData, 0) Then
                Return True
            ElseIf Integer.TryParse(strData, 0) Then
                Return True
            End If
        Catch ex As Exception
            ' Ignore errors here
        End Try

        Return False
    End Function

    Private Shared Function LineContainsValues(strDataLine As String, ParamArray lstValuesToFind As String()) As Boolean
        Dim intMatchCount = 0

        For Each item In lstValuesToFind
            If strDataLine.IndexOf(item, StringComparison.CurrentCultureIgnoreCase) > -1 Then
                intMatchCount += 1
            End If
        Next

        If intMatchCount = lstValuesToFind.Count Then
            Return True
        Else
            Return False
        End If

    End Function

    ''' <summary>
    ''' Returns the index of the indicated column, as tracked by objColumnHeaders
    ''' </summary>
    ''' <param name="strColumnName"></param>
    ''' <param name="objColumnHeaders"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function LookupColumnIndex(strColumnName As String, objColumnHeaders As SortedDictionary(Of String, Integer)) As Integer

        Dim intColIndex As Integer

        If objColumnHeaders.TryGetValue(strColumnName, intColIndex) Then
            If intColIndex >= 0 Then
                Return intColIndex
            End If
        End If

        Return -1

    End Function

    ''' <summary>
    ''' Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
    ''' </summary>
    ''' <returns>The text in the specified column; an empty string if the specific column name is not recognized</returns>
    ''' <remarks></remarks>
    Public Shared Function LookupColumnValue(
      strColumns() As String,
      strColumnName As String,
      objColumnHeaders As SortedDictionary(Of String, Integer)) As String

        Return LookupColumnValue(strColumns, strColumnName, objColumnHeaders, String.Empty)
    End Function

    ''' <summary>
    ''' Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
    ''' </summary>
    ''' <returns>The text in the specified column; strValueIfMissing if the specific column name is not recognized</returns>
    ''' <remarks></remarks>
    Public Shared Function LookupColumnValue(
      strColumns() As String,
      strColumnName As String,
      objColumnHeaders As SortedDictionary(Of String, Integer),
      strValueIfMissing As String) As String

        Dim intColIndex As Integer

        If Not strColumns Is Nothing Then
            intColIndex = LookupColumnIndex(strColumnName, objColumnHeaders)
            If intColIndex >= 0 AndAlso intColIndex < strColumns.Length Then
                If String.IsNullOrWhiteSpace(strColumns(intColIndex)) Then
                    Return String.Empty
                Else
                    Return strColumns(intColIndex)
                End If

            End If
        End If

        ' If we get here, return strValueIfMissing
        Return strValueIfMissing

    End Function

    ''' <summary>
    ''' Returns the value stored in the given named column (using objColumnHeaders to dereference column name with column index)
    ''' </summary>
    ''' <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
    ''' <remarks></remarks>
    Public Shared Function LookupColumnValue(
      strColumns() As String,
      strColumnName As String,
      objColumnHeaders As SortedDictionary(Of String, Integer),
      ValueIfMissing As Integer) As Integer

        Dim strValue As String
        Dim intValue As Integer

        strValue = LookupColumnValue(strColumns, strColumnName, objColumnHeaders, ValueIfMissing.ToString)

        Integer.TryParse(strValue, intValue)

        Return intValue

    End Function

    ''' <summary>
    ''' Returns the value stored in the given named column (using objColumnHeaders to dereference column name with column index)
    ''' </summary>
    ''' <returns>The number in the specified column; 0 if the specific column name is not recognized or the column does not contain a number</returns>
    ''' <remarks></remarks>
    Public Shared Function LookupColumnValue(
      strColumns() As String,
      strColumnName As String,
      objColumnHeaders As SortedDictionary(Of String, Integer),
      ValueIfMissing As Double) As Double

        Dim strValue As String
        Dim dblValue As Double

        strValue = LookupColumnValue(strColumns, strColumnName, objColumnHeaders, ValueIfMissing.ToString)

        Double.TryParse(strValue, dblValue)

        Return dblValue

    End Function

    ''' <summary>
    ''' Updates the column name to column index mapping in objColumnHeaders
    ''' </summary>
    ''' <param name="strColumns">Column names read from the input file</param>
    ''' <param name="objColumnHeaders">Column mapping dictionary object to update</param>
    ''' <remarks>The SortedDictionary object should be instantiated using a case-insensitive comparer, i.e. (StringComparer.CurrentCultureIgnoreCase)</remarks>
    Public Shared Sub ParseColumnHeaders(strColumns() As String, objColumnHeaders As SortedDictionary(Of String, Integer))

        Dim intIndex As Integer

        ' Reset the column indices in objColumnHeaders
        If objColumnHeaders.Count > 0 Then
            Dim strKeys() As String
            ReDim strKeys(objColumnHeaders.Count - 1)
            objColumnHeaders.Keys.CopyTo(strKeys, 0)

            For Each strKey As String In strKeys
                objColumnHeaders(strKey) = -1
            Next
        End If

        For intIndex = 0 To strColumns.Length - 1
            If objColumnHeaders.ContainsKey(strColumns(intIndex)) Then
                ' Update the index associated with this column name
                objColumnHeaders(strColumns(intIndex)) = intIndex
            Else
                ' Ignore this column
            End If
        Next intIndex

    End Sub

    ''' <summary>
    ''' Reads the next line from a synopsis file or first hits file
    ''' </summary>
    ''' <returns>True if a line was read, false if not more data is available</returns>
    ''' <remarks>When FastReadMode is True, you should call FinalizeCurrentPSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
    Public Function MoveNext() As Boolean

        Dim strLineIn As String = String.Empty

        Dim blnSuccess = False
        Dim blnMatchFound = False
        Dim blnUsingCachedPSM = False

        If mCachedLineAvailable Then
            strLineIn = mCachedLine
            mCachedLineAvailable = False
            mPSMCurrent = mCachedPSM
            blnUsingCachedPSM = True
            blnSuccess = True
            mPSMCurrentFinalized = False
            mExtendedScanStatsValid = False
        ElseIf Not mSourceFile.EndOfStream Then
            strLineIn = mSourceFile.ReadLine()
            mSourceFileLinesRead += 1
            blnSuccess = True
        Else
            mCanRead = False
        End If

        If Not blnSuccess OrElse String.IsNullOrEmpty(strLineIn) Then
            Return False
        End If

        If Not mHeaderLineParsed Then
            Dim strSplitLine() = strLineIn.Split(ControlChars.Tab)
            If Not IsNumber(strSplitLine(0)) Then
                ' Parse the header line to confirm the column ordering
                mPHRPParser.ParseColumnHeaders(strSplitLine)

                mHeaderLineParsed = True
                Return MoveNext()
            End If

            mHeaderLineParsed = True
        End If

        If Not blnUsingCachedPSM Then
            mPSMCurrent = New clsPSM()
            blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mSourceFileLinesRead, mPSMCurrent, mFastReadMode)
            mPSMCurrentFinalized = False
            mExtendedScanStatsValid = False
        End If

        If Not blnSuccess Then
            Return False
        End If

        blnMatchFound = True

        ' The PHRPParser will update .PeptideWithNumericMods if the _SeqInfo.txt file is loaded
        ' If it wasn't loaded, then this class can update .PeptideWithNumericMods and .PeptideMods 
        ' by inferring the mods using mDynamicMods and mStaticMods (which were populated using the PHRP ModSummary file)
        If Not mFastReadMode AndAlso mStartupOptions.LoadModsAndSeqInfo AndAlso String.IsNullOrEmpty(mPSMCurrent.PeptideWithNumericMods) Then
            MarkupPeptideWithMods()
        End If

        Dim objScanStatsInfo As clsScanStatsInfo = Nothing

        Dim blnScanStatsValid = TryGetScanStats(mPSMCurrent.ScanNumber, objScanStatsInfo)
        mExtendedScanStatsValid = TryGetExtendedScanStats(mPSMCurrent.ScanNumber, mExtendedScanStatsInfo)

        If blnScanStatsValid Then
            ' Update the elution time
            mPSMCurrent.ElutionTimeMinutes = objScanStatsInfo.ScanTimeMinutes
        End If

        If String.IsNullOrEmpty(mPSMCurrent.CollisionMode) OrElse mPSMCurrent.CollisionMode = clsPSM.UNKNOWN_COLLISION_MODE Then

            ' Determine the ScanTypeName using the the ScanStats or ExtendedScanStats info
            If blnScanStatsValid AndAlso Not String.IsNullOrEmpty(objScanStatsInfo.ScanTypeName) Then
                mPSMCurrent.CollisionMode = GetCollisionMode(objScanStatsInfo.ScanTypeName)
            End If

            If String.IsNullOrEmpty(mPSMCurrent.CollisionMode) AndAlso mExtendedScanStatsValid Then
                ' Scan type still not determined, but Extended Scan Stats data is available
                If Not String.IsNullOrEmpty(mExtendedScanStatsInfo.CollisionMode) Then
                    ' Check for Collision mode being "0"
                    ' This is often the case for the first scan in a Thermo .Raw file
                    If mExtendedScanStatsInfo.CollisionMode <> "0" Then
                        mPSMCurrent.CollisionMode = mExtendedScanStatsInfo.CollisionMode.ToUpper()
                    End If
                End If
            End If
        End If

        If Not mFastReadMode Then
            If mPHRPParser.PeptideHitResultType = ePeptideHitResultType.Sequest OrElse mPHRPParser.PeptideHitResultType = ePeptideHitResultType.XTandem Then
                ComputePrecursorNeutralMass()
            End If
        End If

        If Not mMSGFCachedResults Is Nothing AndAlso mMSGFCachedResults.Count > 0 Then
            Dim strMSGFSpecProb As String = String.Empty
            Dim dblSpecProb As Double

            If mMSGFCachedResults.TryGetValue(mPSMCurrent.ResultID, strMSGFSpecProb) Then
                mPSMCurrent.MSGFSpecProb = strMSGFSpecProb
                If strMSGFSpecProb.Length > 12 Then
                    ' Attempt to shorten the SpecProb value
                    If Double.TryParse(strMSGFSpecProb, dblSpecProb) Then
                        mPSMCurrent.MSGFSpecProb = dblSpecProb.ToString("0.00000E-00")
                    End If
                End If
            End If
        End If

        If mSkipDuplicatePSMs Then

            ' Read the next line and check whether it's the same hit, but a different protein
            Dim blnReadNext = True
            Do While blnReadNext AndAlso Not mSourceFile.EndOfStream
                strLineIn = mSourceFile.ReadLine()
                mSourceFileLinesRead += 1

                If Not String.IsNullOrEmpty(strLineIn) Then

                    Dim objNewPSM As New clsPSM()
                    blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mSourceFileLinesRead, objNewPSM, mFastReadMode)

                    ' Check for duplicate lines
                    ' If this line is a duplicate of the previous line, then skip it
                    ' This happens in Sequest _syn.txt files where the line is repeated for all protein matches
                    ' It can also happen in MSGF+ results, though the prefix and suffix residues could differ for the same peptide, depending on the protein context
                    With mPSMCurrent

                        Dim isDuplicate = False

                        If .ScanNumber = objNewPSM.ScanNumber AndAlso
                           .Charge = objNewPSM.Charge Then

                            If String.Equals(.Peptide, objNewPSM.Peptide) Then
                                isDuplicate = True
                            Else
                                Dim strPeptide1 As String = String.Empty
                                Dim strPeptide2 As String = String.Empty

                                If clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(.Peptide, strPeptide1, "", "") AndAlso
                                   clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(objNewPSM.Peptide, strPeptide2, "", "") Then
                                    If String.Equals(strPeptide1, strPeptide2) Then
                                        isDuplicate = True
                                    End If
                                End If
                            End If
                        End If

                        If isDuplicate Then
                            ' Update the protein list
                            Dim addnlProteins = objNewPSM.Proteins.Except(.Proteins, StringComparer.CurrentCultureIgnoreCase)
                            If addnlProteins.Count > 0 Then
                                .Proteins.AddRange(addnlProteins)
                            End If
                        Else
                            blnReadNext = False
                            mCachedLine = String.Copy(strLineIn)
                            mCachedLineAvailable = True
                            mCachedPSM = objNewPSM
                        End If

                    End With
                End If
            Loop

        End If

        Return blnMatchFound

    End Function

    Private Sub ComputePrecursorNeutralMass()
        Dim dblMonoisotopicPrecursorMass As Double

        If Not mExtendedScanStatsValid Then Exit Sub

        dblMonoisotopicPrecursorMass = 0

        ' Try to extract out the precursor m/z value from the "Scan Filter Text" field
        Dim dblParentIonMZ As Double

        If ExtractParentIonMzFromFilterText(mExtendedScanStatsInfo.ScanFilterText, dblParentIonMZ) Then
            If dblParentIonMZ > 0 Then
                dblMonoisotopicPrecursorMass = mPeptideMassCalculator.ConvoluteMass(
                    dblParentIonMZ, mPSMCurrent.Charge, 0)
            End If
        End If

        If Math.Abs(dblMonoisotopicPrecursorMass) < Double.Epsilon Then
            If mExtendedScanStatsInfo.MonoisotopicMZ > 0 Then
                ' Determine the precursor m/z value using the Monoisotopic m/z value reported by the instrument
                dblMonoisotopicPrecursorMass = mPeptideMassCalculator.ConvoluteMass(
                    mExtendedScanStatsInfo.MonoisotopicMZ, mPSMCurrent.Charge, 0)
            End If
        End If

        If dblMonoisotopicPrecursorMass > 0 Then
            mPSMCurrent.PrecursorNeutralMass = dblMonoisotopicPrecursorMass
        End If

    End Sub

    ''' <summary>
    ''' This function extracts the Parent Ion m/z from the filter string
    ''' </summary>
    ''' <param name="filterText"></param>
    ''' <param name="parentIonMz"></param>
    ''' <returns>True if parsing successful</returns>
    ''' <remarks>The original version of this code (C#) is in ThermoRawFileReader.XRawFileIO.ExtractParentIonMZFromFilterText(string, out double)</remarks>
    Private Function ExtractParentIonMzFromFilterText(filterText As String, <Out> ByRef parentIonMz As Double) As Boolean
        Dim matcher As Regex

        If filterText.ToLower().Contains("msx") Then
            matcher = mParentIonMzMatchMsx
        Else
            matcher = mParentIonMzMatchNonMsx
        End If

        Dim reMatch = matcher.Match(filterText)

        If reMatch.Success Then
            Dim parentIonMzText = reMatch.Groups("ParentMZ").Value
            Dim success = Double.TryParse(parentIonMzText, parentIonMz)
            Return success
        End If

        parentIonMz = 0
        Return False

    End Function

    Private Sub MarkupPeptideWithMods()
        Dim blnSuccess As Boolean

        ' Markup the peptide with the dynamic and static mods
        Dim strPeptideWithMods As String = String.Empty
        Dim lstPeptideMods As List(Of clsAminoAcidModInfo) = Nothing

        blnSuccess = ConvertModsToNumericMods(mPSMCurrent.Peptide.Trim, strPeptideWithMods, lstPeptideMods)
        If blnSuccess Then
            Dim dblTotalModMass As Double = 0

            mPSMCurrent.PeptideWithNumericMods = strPeptideWithMods
            mPSMCurrent.ClearModifiedResidues()
            For Each objModEntry In lstPeptideMods
                mPSMCurrent.AddModifiedResidue(objModEntry)
                dblTotalModMass += objModEntry.ModDefinition.ModificationMass
            Next

            If Math.Abs(mPSMCurrent.PeptideMonoisotopicMass) < Double.Epsilon Then
                mPSMCurrent.PeptideMonoisotopicMass = mPeptideMassCalculator.ComputeSequenceMass(mPSMCurrent.PeptideCleanSequence) + dblTotalModMass
            End If

        End If
    End Sub

    ''' <summary>
    ''' When FastReadMode is True, you first call MoveNext to read the peptide scores, then if the peptide
    ''' is a peptide of interest, you call this function to finalize any processing steps that were skipped
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub FinalizeCurrentPSM()

        If mPSMCurrentFinalized Then Exit Sub

        ' Determine the clean sequence and cleavage state, and update the Seq_ID fields
        mPHRPParser.FinalizePSM(mPSMCurrent)

        MarkupPeptideWithMods()

        ComputePrecursorNeutralMass()

        mPSMCurrentFinalized = True

    End Sub

    Private Function ReadAndCacheMSGFData() As Boolean

        Dim strMSGFFilePath As String
        Dim blnSuccess = False

        Try
            strMSGFFilePath = GetMSGFFileName(mInputFilePath)
            strMSGFFilePath = Path.Combine(mInputFolderPath, strMSGFFilePath)

            blnSuccess = ReadAndCacheMSGFData(strMSGFFilePath)

        Catch ex As Exception
            HandleException("Exception determining MSGF file path", ex)
        End Try

        Return blnSuccess

    End Function

    Private Function ReadAndCacheMSGFData(strMSGFFilePath As String) As Boolean
        Dim blnSuccess = False

        Try
            strMSGFFilePath = AutoSwitchToFHTIfRequired(strMSGFFilePath, mInputFilePath)

            If File.Exists(strMSGFFilePath) Then

                Dim oMSGFReader As New clsMSGFResultsReader()
                mMSGFCachedResults = oMSGFReader.ReadMSGFData(strMSGFFilePath)

                If oMSGFReader.ErrorMessage.Length > 0 Then
                    ReportError("Error reading MSGF data: " & oMSGFReader.ErrorMessage)
                Else
                    blnSuccess = True
                End If
            Else
                ReportWarning("MSGF file not found: " & strMSGFFilePath)
            End If
        Catch ex As Exception
            HandleException("Exception reading MSGF file", ex)
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Reads the data in strModSummaryFilePath.  Populates objDynamicMods and objStaticMods with the modification definitions
    ''' </summary>
    ''' <param name="strModSummaryFilePath">Path to the PHRP Mod Summary file to read</param>
    ''' <param name="objDynamicMods">List with mod symbols as the key and the corresponding mod mass</param>
    ''' <param name="objStaticMods">List with amino acid names as the key and the corresponding mod mass</param>
    ''' <returns>True if success; false if an error</returns>
    Private Function ReadModSummaryFile(
      strModSummaryFilePath As String,
      objDynamicMods As SortedDictionary(Of Char, clsModificationDefinition),
      objStaticMods As SortedDictionary(Of String, List(Of clsModificationDefinition))) As Boolean

        Dim objModSummaryReader As clsPHRPModSummaryReader

        Dim lstModDefs As List(Of clsModificationDefinition) = Nothing

        Dim strModMass As String
        Dim blnSuccess As Boolean

        Try
            ' Clear objDynamicMods and objStaticMods (should have been instantiated by the calling function)
            objDynamicMods.Clear()
            objStaticMods.Clear()

            If String.IsNullOrEmpty(strModSummaryFilePath) Then
                ReportError("ModSummaryFile path is empty; unable to continue")
                Return False
            ElseIf Not File.Exists(strModSummaryFilePath) Then
                ReportError("ModSummary file not found: " & strModSummaryFilePath)
                Return False
            End If

            ShowMessage("Reading the PHRP ModSummary file")

            objModSummaryReader = New clsPHRPModSummaryReader(strModSummaryFilePath)
            blnSuccess = objModSummaryReader.Success

            If blnSuccess AndAlso objModSummaryReader.ModificationDefs.Count > 0 Then

                For Each objModDef In objModSummaryReader.ModificationDefs

                    strModMass = objModSummaryReader.GetModificationMassAsText(objModDef.MassCorrectionTag)

                    Select Case objModDef.ModificationType
                        Case eModificationTypeConstants.StaticMod, eModificationTypeConstants.TerminalPeptideStaticMod, eModificationTypeConstants.ProteinTerminusStaticMod

                            ' "S", "T", or "P"
                            ' Static residue mod, peptide terminus static mod, or protein terminus static mod
                            ' Note that < and > mean peptide N and C terminus (N_TERMINAL_PEPTIDE_SYMBOL_DMS and C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                            ' Note that [ and ] mean protein N and C terminus (N_TERMINAL_PROTEIN_SYMBOL_DMS and C_TERMINAL_PROTEIN_SYMBOL_DMS)

                            ' This mod could apply to multiple residues, so need to process each character in strTargetResidues
                            For Each chChar In objModDef.TargetResidues
                                Try

                                    If objStaticMods.TryGetValue(chChar, lstModDefs) Then
                                        If Not lstModDefs.Contains(objModDef) Then
                                            ' Residue is already present in objStaticMods; this is unusual, but we'll allow it
                                            ' We'll log a warning, but continue
                                            ShowMessage("Warning: Residue '" & chChar & "' has more than one static mod defined; this is not typically used, but will be allowed")
                                            lstModDefs.Add(objModDef)
                                        End If
                                    Else
                                        lstModDefs = New List(Of clsModificationDefinition)
                                        lstModDefs.Add(objModDef)
                                        objStaticMods.Add(chChar, lstModDefs)
                                    End If

                                Catch ex As Exception
                                    HandleException("Exception adding static mod for " & chChar & " with ModMass=" & strModMass, ex)
                                End Try
                            Next
                        Case eModificationTypeConstants.DynamicMod
                            ' Dynamic residue mod (Includes mod type "D")
                            ' Note that < and > mean peptide N and C terminus (N_TERMINAL_PEPTIDE_SYMBOL_DMS and C_TERMINAL_PEPTIDE_SYMBOL_DMS)

                            Try
                                If objDynamicMods.ContainsKey(objModDef.ModificationSymbol) Then
                                    ' Mod symbol already present in objDynamicMods; this is unexpected
                                    ' We'll log a warning, but continue
                                    ShowMessage("Warning: Dynamic mod symbol '" & objModDef.ModificationSymbol & "' is already defined; it cannot have more than one associated mod mass (duplicate has ModMass=" & strModMass & ")")
                                Else
                                    objDynamicMods.Add(objModDef.ModificationSymbol, objModDef)
                                End If

                            Catch ex As Exception
                                HandleException("Exception adding dynamic mod for " & objModDef.ModificationSymbol & " with ModMass=" & strModMass, ex)
                            End Try


                        Case eModificationTypeConstants.IsotopicMod
                            ' Isotopic mods are not supported by this class
                            ' However, do not log a warning since these are rarely used

                        Case eModificationTypeConstants.UnknownType
                            ' Unknown type; just ignore it

                        Case Else
                            ' Unknown type; just ignore it

                    End Select
                Next

            End If

        Catch ex As Exception
            HandleException("Exception reading PHRP Mod Summary file", ex)
            Return False
        End Try

        Return True

    End Function

    Private Function ReadScanStatsData() As Boolean

        Dim strScanStatsFilePath As String
        Dim blnSuccess = False

        Try
            strScanStatsFilePath = GetScanStatsFilename(mDatasetName)
            strScanStatsFilePath = Path.Combine(mInputFolderPath, strScanStatsFilePath)

            blnSuccess = ReadScanStatsData(strScanStatsFilePath)

        Catch ex As Exception
            HandleException("Exception determining Scan Stats file path", ex)
        End Try

        Return blnSuccess

    End Function

    Private Function ReadScanStatsData(strScanStatsFilePath As String) As Boolean

        Dim blnSuccess = False

        Try
            If File.Exists(strScanStatsFilePath) Then

                Dim oScanStatsReader As New clsScanStatsReader()
                mScanStats = oScanStatsReader.ReadScanStatsData(strScanStatsFilePath)

                If oScanStatsReader.ErrorMessage.Length > 0 Then
                    ReportError("Error reading ScanStats data: " & oScanStatsReader.ErrorMessage)
                Else
                    blnSuccess = True
                End If
            Else
                ReportWarning("ScanStats file not found: " & strScanStatsFilePath)
            End If

        Catch ex As Exception
            HandleException("Exception reading Scan Stats file", ex)
        End Try

        Return blnSuccess

    End Function

    Private Function ReadExtendedScanStatsData() As Boolean

        Dim strExtendedScanStatsFilePath As String
        Dim blnSuccess = False

        Try
            strExtendedScanStatsFilePath = GetExtendedScanStatsFilename(mDatasetName)
            strExtendedScanStatsFilePath = Path.Combine(mInputFolderPath, strExtendedScanStatsFilePath)

            blnSuccess = ReadExtendedScanStatsData(strExtendedScanStatsFilePath)

        Catch ex As Exception
            HandleException("Exception determining Scan Stats file path", ex)
        End Try

        Return blnSuccess

    End Function

    Private Function ReadExtendedScanStatsData(strExtendedScanStatsFilePath As String) As Boolean

        Dim blnSuccess = False

        Try
            If File.Exists(strExtendedScanStatsFilePath) Then

                Dim oExtendedScanStatsReader As New clsExtendedScanStatsReader()
                mScanStatsEx = oExtendedScanStatsReader.ReadExtendedScanStatsData(strExtendedScanStatsFilePath)

                If oExtendedScanStatsReader.ErrorMessage.Length > 0 Then
                    ReportError("Error reading Extended ScanStats data: " & oExtendedScanStatsReader.ErrorMessage)
                Else
                    blnSuccess = True
                End If
            Else
                ' Note: we do not need to raise a warning for MSGFDB results since the extended scan stats file isn't needed
                If mPHRPParser.PeptideHitResultType <> ePeptideHitResultType.MSGFDB Then
                    ReportWarning("Extended ScanStats file not found: " & strExtendedScanStatsFilePath)
                End If
            End If

        Catch ex As Exception
            HandleException("Exception reading Extended Scan Stats file", ex)
        End Try

        Return blnSuccess

    End Function

    Private Sub ReportError(strErrorMessage As String)
        mErrorMessage = strErrorMessage
        If mEchoMessagesToConsole Then Console.WriteLine(strErrorMessage)
        mErrorMessages.Add(strErrorMessage)

        RaiseEvent ErrorEvent(strErrorMessage)
    End Sub

    Private Sub ReportWarning(strWarningMessage As String)
        If mEchoMessagesToConsole Then Console.WriteLine(strWarningMessage)
        mWarningMessages.Add(strWarningMessage)

        RaiseEvent WarningEvent(strWarningMessage)
    End Sub

    Private Sub SetLocalErrorCode(eNewErrorCode As ePHRPReaderErrorCodes)
        SetLocalErrorCode(eNewErrorCode, False)
    End Sub

    Private Sub SetLocalErrorCode(eNewErrorCode As ePHRPReaderErrorCodes, blnLeaveExistingErrorCodeUnchanged As Boolean)

        If blnLeaveExistingErrorCodeUnchanged AndAlso mLocalErrorCode <> ePHRPReaderErrorCodes.NoError Then
            ' An error code is already defined; do not change it
        Else
            mLocalErrorCode = eNewErrorCode
        End If

    End Sub

    Private Sub ShowMessage(strMessage As String)
        If mEchoMessagesToConsole Then Console.WriteLine(strMessage)
        RaiseEvent MessageEvent(strMessage)
    End Sub

    Private Function TryGetScanStats(intScanNumber As Integer, <Out()> ByRef objScanStatsInfo As clsScanStatsInfo) As Boolean
        If Not mScanStats Is Nothing AndAlso mScanStats.Count > 0 Then
            If mScanStats.TryGetValue(intScanNumber, objScanStatsInfo) Then
                Return True
            End If
        End If
        objScanStatsInfo = Nothing
        Return False
    End Function

    Private Function TryGetExtendedScanStats(intScanNumber As Integer, <Out()> ByRef objExtendedScanStatsInfo As clsScanStatsExInfo) As Boolean
        If Not mScanStatsEx Is Nothing AndAlso mScanStats.Count > 0 Then
            If mScanStatsEx.TryGetValue(intScanNumber, objExtendedScanStatsInfo) Then
                Return True
            End If
        End If
        objExtendedScanStatsInfo = Nothing
        Return False
    End Function

    Private Function ValidateInputFiles(
      strInputFilePath As String,
      ByRef eResultType As ePeptideHitResultType,
      ByRef strModSummaryFilePath As String) As Boolean

        Dim fiFileInfo As FileInfo

        fiFileInfo = New FileInfo(strInputFilePath)
        If Not fiFileInfo.Exists Then
            SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
            ReportError("Input file not found: " & strInputFilePath)
            Return False
        End If

        ' Try to auto-determine the result type if it is not specified
        If eResultType = ePeptideHitResultType.Unknown Then
            eResultType = AutoDetermineResultType(strInputFilePath)
        End If

        If eResultType = ePeptideHitResultType.Unknown Then
            SetLocalErrorCode(ePHRPReaderErrorCodes.InputFileFormatNotRecognized)
            ReportError("Error: Unable to auto-determine file format for " & strInputFilePath)
            Return False
        End If

        ' Extract the dataset name from the input file path
        mDatasetName = AutoDetermineDatasetName(strInputFilePath, eResultType)
        If mDatasetName Is Nothing OrElse mDatasetName.Length = 0 Then
            If mStartupOptions.LoadModsAndSeqInfo OrElse mStartupOptions.LoadMSGFResults OrElse mStartupOptions.LoadScanStatsData Then
                ReportError("Error: Unable to auto-determine the dataset name from the input file name: " & strInputFilePath)
                SetLocalErrorCode(ePHRPReaderErrorCodes.InputFileFormatNotRecognized)
                Return False
            Else
                ReportWarning("Unable to auto-determine the dataset name from the input file name; this is not a critical error since not reading related files: " & strInputFilePath)
            End If
        End If

        If mStartupOptions.LoadModsAndSeqInfo Then
            strModSummaryFilePath = GetPHRPModSummaryFileName(eResultType, mDatasetName)
            strModSummaryFilePath = Path.Combine(fiFileInfo.DirectoryName, strModSummaryFilePath)

            Dim strModSummaryFilePathPreferred As String
            strModSummaryFilePathPreferred = AutoSwitchToFHTIfRequired(strModSummaryFilePath, fiFileInfo.Name)
            If strModSummaryFilePath <> strModSummaryFilePathPreferred AndAlso File.Exists(strModSummaryFilePathPreferred) Then
                strModSummaryFilePath = strModSummaryFilePathPreferred
            End If

            If Not ValidateRequiredFileExists("ModSummary file", strModSummaryFilePath) AndAlso fiFileInfo.Name.ToLower().Contains("_fht") Then
                SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound)
                Return False
            End If
        Else
            strModSummaryFilePath = String.Empty
        End If

        Return True

    End Function

    Private Function ValidateRequiredFileExists(strFileDescription As String, strFilePath As String) As Boolean
        Return ValidateRequiredFileExists(strFileDescription, strFilePath, True)
    End Function

    Private Function ValidateRequiredFileExists(strFileDescription As String, strFilePath As String, blnReportErrors As Boolean) As Boolean

        If strFilePath Is Nothing OrElse strFilePath.Length = 0 Then
            If blnReportErrors Then
                ReportError(strFileDescription & " is not defined")
            End If
            Return False

        ElseIf Not File.Exists(strFilePath) Then
            If blnReportErrors Then
                ReportError(strFileDescription & " not found: " & strFilePath)
            End If
            Return False

        End If

        Return True

    End Function

    Private Sub mPHRPParser_ErrorEvent(strErrorMessage As String)
        ReportError(strErrorMessage)
    End Sub

    Private Sub mPHRPParser_MessageEvent(strMessage As String)
        ShowMessage(strMessage)
    End Sub

    Private Sub mPHRPParser_WarningEvent(strWarningMessage As String)
        ReportWarning(strWarningMessage)
    End Sub

#Region "IDisposable Support"
    Private disposedValue As Boolean ' To detect redundant calls

    ' IDisposable
    Protected Overridable Sub Dispose(disposing As Boolean)
        If Not Me.disposedValue Then
            If disposing Then
                If Not mSourceFile Is Nothing Then
                    mSourceFile.Close()
                End If
            End If

        End If
        Me.disposedValue = True
    End Sub

    ' This code added by Visual Basic to correctly implement the disposable pattern.
    Public Sub Dispose() Implements IDisposable.Dispose
        ' Do not change this code.  Put cleanup code in Dispose(disposing As Boolean) above.
        Dispose(True)
        GC.SuppressFinalize(Me)
    End Sub
#End Region

End Class
