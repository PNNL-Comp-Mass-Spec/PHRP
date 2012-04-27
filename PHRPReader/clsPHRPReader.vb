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

Imports System.Collections.Generic
Imports PHRPReader.clsModificationDefinition

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

	Public Enum ePeptideHitResultType
		Unknown = 0
		Sequest = 1
		XTandem = 2
		Inspect = 3
		MSGFDB = 4
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
	Protected mDatasetName As String
	Protected mInputFilePath As String
	Protected mInputFolderPath As String

	Protected mSkipDuplicatePSMs As Boolean

	Protected mLoadModsAndSeqInfo As Boolean
	Protected mLoadMSGFResults As Boolean
	Protected mLoadScanStatsData As Boolean

	Protected mEchoMessagesToConsole As Boolean

	Protected mCanRead As Boolean
	Protected mInitialized As Boolean
	Protected mModSummaryFileLoaded As Boolean

	Protected mSourceFile As System.IO.StreamReader
	Protected mSourceFileLineCount As Integer
	Protected mSourceFileLinesRead As Integer

	Protected WithEvents mPHRPParser As clsPHRPParser
	Protected mPeptideMassCalculator As clsPeptideMassCalculator

	' This dictionary contains mod symbols as the key and modification definition as the values
	Protected mDynamicMods As SortedDictionary(Of Char, clsModificationDefinition)

	' This dictionary contains amino acid names as the key and the corresponding mod modification (or mod modifications) 
	Protected mStaticMods As SortedDictionary(Of String, List(Of clsModificationDefinition))

	' This dictionary tracks the MSGFSpecProb values for each entry in the source file
	' The keys are Result_ID and the string is MSGFSpecProb (stored as string to preserve formatting)
	Protected mMSGFCachedResults As Dictionary(Of Integer, String)

	' This dictionary tracks scan stats values, in particular elution time
	'The keys are ScanNumber and values are clsScanStatsInfo objects
	Protected mScanStats As Dictionary(Of Integer, clsScanStatsInfo)

	' This dictionary tracks extended scan stats values, in particular parent ion mz (via MonoisotopicMZ)
	'The keys are ScanNumber and values are clsScanStatsExInfo objects
	Protected mScanStatsEx As Dictionary(Of Integer, clsScanStatsExInfo)

	Protected mPSMCurrent As clsPSM

	Protected mHeaderLineParsed As Boolean
	Protected mCachedLineAvailable As Boolean
	Protected mCachedLine As String

	Protected mErrorMessages As System.Collections.Generic.List(Of String)
	Protected mWarningMessages As System.Collections.Generic.List(Of String)

	Protected mErrorMessage As String = String.Empty
	Protected mLocalErrorCode As ePHRPReaderErrorCodes

#End Region

#Region "Events"
	Public Event MessageEvent(ByVal strMessage As String)
	Public Event ErrorEvent(ByVal strErrorMessage As String)
	Public Event WarningEvent(ByVal strWarningMessage As String)
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
	Public ReadOnly Property ErrorMessages() As System.Collections.Generic.List(Of String)
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
	''' If True, then looks for and loads the modification definitions from the _ModSummary.txt file associated with the input file
	''' Also reads the SeqInfo and related files
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property LoadModsAndSeqInfo() As Boolean
		Get
			Return mLoadModsAndSeqInfo
		End Get
		Set(value As Boolean)
			mLoadModsAndSeqInfo = value
		End Set
	End Property

	''' <summary>
	''' If true, then loads the MSGF SpecProb values from the _MSGF.txt file associated with the input file
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property LoadMSGFResults() As Boolean
		Get
			Return mLoadMSGFResults
		End Get
		Set(value As Boolean)
			mLoadMSGFResults = value
		End Set
	End Property

	''' <summary>
	''' If True, then loads the MASIC _ScanStats.txt file
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property LoadScanStatsData() As Boolean
		Get
			Return mLoadScanStatsData
		End Get
		Set(value As Boolean)
			mLoadScanStatsData = value
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
	Public ReadOnly Property WarningMessages() As System.Collections.Generic.List(Of String)
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
	Public Sub New(ByVal strInputFilePath As String)
		Me.New(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo:=True, blnLoadMSGFResults:=True, blnLoadScanStats:=False)
	End Sub

	''' <summary>
	''' Constructor where the PeptideHit result type is explicitly set
	''' </summary>
	''' <param name="strInputFilePath">Input file to read</param>
	''' <param name="eResultType">Source file PeptideHit result type</param>
	''' <remarks>Sets LoadModSummaryFile to True and LoadMSGFResults to true</remarks>
	Public Sub New(ByVal strInputFilePath As String, eResultType As ePeptideHitResultType)
		Me.New(strInputFilePath, eResultType, blnLoadModsAndSeqInfo:=True, blnLoadMSGFResults:=True, blnLoadScanStats:=False)
	End Sub

	''' <summary>
	''' Constructor that auto-determines the PeptideHit result type based on the filename
	''' </summary>
	''' <param name="strInputFilePath">Input file to read</param>
	''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
	''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strInputFilePath As String, ByVal blnLoadModsAndSeqInfo As Boolean, ByVal blnLoadMSGFResults As Boolean)
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
	Public Sub New(ByVal strInputFilePath As String, ByVal blnLoadModsAndSeqInfo As Boolean, ByVal blnLoadMSGFResults As Boolean, ByVal blnLoadScanStats As Boolean)
		Me.New(strInputFilePath, ePeptideHitResultType.Unknown, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats)
	End Sub

	''' <summary>
	''' Constructor where the PeptideHit result type is explicitly set
	''' </summary>
	''' <param name="strInputFilePath">Input file to read</param>
	''' ''' <param name="eResultType">Source file PeptideHit result type</param>
	''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
	''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strInputFilePath As String, eResultType As ePeptideHitResultType, ByVal blnLoadModsAndSeqInfo As Boolean, ByVal blnLoadMSGFResults As Boolean)
		Me.New(strInputFilePath, eResultType, blnLoadModsAndSeqInfo, blnLoadMSGFResults, blnLoadScanStats:=False)
	End Sub

	''' <summary>
	''' Constructor where the PeptideHit result type is explicitly set
	''' </summary>
	''' <param name="strInputFilePath">Input file to read</param>
	''' ''' <param name="eResultType">Source file PeptideHit result type</param>
	''' <param name="blnLoadModsAndSeqInfo">If True, then looks for and auto-loads the modification definitions from the _moddefs.txt file</param>
	''' <param name="blnLoadMSGFResults">If True, then looks for and auto-loads the MSGF results from the _msg.txt file</param>
	''' <param name="blnLoadScanStats">If True, then looks for and auto-loads the MASIC scan stats files (used to determine collision mode and to refine the precursor m/z values)</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strInputFilePath As String, eResultType As ePeptideHitResultType, ByVal blnLoadModsAndSeqInfo As Boolean, ByVal blnLoadMSGFResults As Boolean, ByVal blnLoadScanStats As Boolean)

		mInitialized = False

		Me.InitializeMemberVariables()

		mLoadModsAndSeqInfo = blnLoadModsAndSeqInfo
		mLoadMSGFResults = blnLoadMSGFResults
		mLoadScanStatsData = blnLoadScanStats

		InitializeReader(strInputFilePath, eResultType)

		mInitialized = True
	End Sub

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

	Protected Function CountLines(strTextFilePath As String) As Integer

		Dim intTotalLines As Integer

		Try
			intTotalLines = 0
			Using srReader As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strTextFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.ReadWrite))
				While srReader.Peek > -1
					srReader.ReadLine()
					intTotalLines += 1
				End While
			End Using

		Catch ex As Exception
			Throw New Exception("Error counting the lines in " & System.IO.Path.GetFileName(strTextFilePath) & ": " & ex.Message, ex)
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

	Private Sub InitializeMemberVariables()
		mDatasetName = String.Empty
		mInputFilePath = String.Empty
		mInputFolderPath = String.Empty

		mCanRead = False
		mModSummaryFileLoaded = False

		mSkipDuplicatePSMs = True

		mEchoMessagesToConsole = False

		mLoadModsAndSeqInfo = True
		mLoadMSGFResults = True
		mLoadScanStatsData = False

		mErrorMessage = String.Empty
		mLocalErrorCode = ePHRPReaderErrorCodes.NoError

		mMSGFCachedResults = New Dictionary(Of Integer, String)

		mDynamicMods = New SortedDictionary(Of Char, clsModificationDefinition)
		mStaticMods = New SortedDictionary(Of String, List(Of clsModificationDefinition))

		mPeptideMassCalculator = New clsPeptideMassCalculator()

		mErrorMessages = New System.Collections.Generic.List(Of String)
		mWarningMessages = New System.Collections.Generic.List(Of String)

		mSourceFileLineCount = 0

	End Sub

	Protected Function InitializeReader(ByVal strInputFilePath As String, ByVal eResultType As ePeptideHitResultType) As Boolean

		Dim strModSummaryFilePath As String = String.Empty

		Dim blnSuccess As Boolean

		Try
			If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
				ReportError("Input file name is empty")
				SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
				If Not mInitialized Then Throw New System.IO.FileNotFoundException(mErrorMessage)
				Return False
			Else
				' Confirm that the source file exists
				' Make sure strInputFilePath points to a valid file
				Dim fiFileInfo As System.IO.FileInfo
				fiFileInfo = New System.IO.FileInfo(strInputFilePath)

				mInputFolderPath = fiFileInfo.DirectoryName
				mInputFilePath = fiFileInfo.FullName

				If Not fiFileInfo.Exists Then
					ReportError("Input file not found: " & strInputFilePath)
					SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
					If Not mInitialized Then Throw New System.IO.FileNotFoundException(mErrorMessage)
					Return False
				Else

					' Note that the following populates mDatasetName
					blnSuccess = ValidateInputFiles(strInputFilePath, eResultType, strModSummaryFilePath)
					If Not blnSuccess Then
						SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound, True)
						If Not mInitialized Then Throw New System.IO.FileNotFoundException(mErrorMessage)
						Return False
					End If

					' Open the input file for reading
					' Note that this will also load the MSGFSpecProb info and ScanStats info
					blnSuccess = InitializeParser(eResultType)

					If blnSuccess AndAlso mLoadModsAndSeqInfo Then
						' Read the PHRP Mod Summary File to populate mDynamicMods and mStaticMods
						' Note that the PHRPParser also loads the ModSummary file, and that mDynamicMods and mStaticMods are only used if the _SeqInfo.txt file is not found
						blnSuccess = ReadModSummaryFile(strModSummaryFilePath, mDynamicMods, mStaticMods)
						If Not blnSuccess Then
							SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound, True)
							If Not mInitialized Then Throw New System.IO.FileNotFoundException(mErrorMessage)
							Return False
						Else
							mModSummaryFileLoaded = True
						End If
					End If

					If blnSuccess AndAlso mLoadMSGFResults Then
						' Cache the MSGF values (if present)
						ReadAndCacheMSGFData()
					End If

					If blnSuccess AndAlso mLoadScanStatsData Then
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

	Protected Function InitializeParser(ByVal eResultType As ePeptideHitResultType) As Boolean

		Dim blnSuccess As Boolean = True

		Try
			If String.IsNullOrEmpty(mDatasetName) Then
				ReportError("Dataset name is undefined; unable to continue")
				Return False
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

					' Convert Sequest results to input format required for MSGF
					mPHRPParser = New clsPHRPParserSequest(mDatasetName, mInputFilePath, mLoadModsAndSeqInfo)

				Case ePeptideHitResultType.XTandem

					' Convert X!Tandem results to input format required for MSGF
					' Note that Result to Protein mapping will be auto-loaded during instantiation of mPHRPParser
					mPHRPParser = New clsPHRPParserXTandem(mDatasetName, mInputFilePath, mLoadModsAndSeqInfo)

				Case ePeptideHitResultType.Inspect

					' Convert Inspect results to input format required for MSGF
					mPHRPParser = New clsPHRPParserInspect(mDatasetName, mInputFilePath, mLoadModsAndSeqInfo)

				Case ePeptideHitResultType.MSGFDB

					' Convert MSGFDB results to input format required for MSGF
					mPHRPParser = New clsPHRPParserMSGFDB(mDatasetName, mInputFilePath, mLoadModsAndSeqInfo)

				Case Else
					'Should never get here; invalid result type specified
					ReportError("Invalid PeptideHit ResultType specified: " & eResultType)
					blnSuccess = False
			End Select

			If blnSuccess Then

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
				mSourceFile = New System.IO.StreamReader(New System.IO.FileStream(mInputFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))
				mCanRead = True

			End If

		Catch ex As Exception
			HandleException("Error in InitializeParser", ex)
			If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
		End Try

		Return blnSuccess

	End Function

	''' <summary>
	''' Auto-determine the dataset name using the input file path
	''' </summary>
	''' <param name="strFilePath"></param>
	''' <returns>Dataset name</returns>
	''' <remarks>Returns an empty string if unable to determine the dataset name</remarks>
	Public Shared Function AutoDetermineDatasetName(ByVal strFilePath As String) As String
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
	Public Shared Function AutoDetermineDatasetName(ByVal strFilePath As String, ByVal eResultType As ePeptideHitResultType) As String

		Dim strDatasetName As String = String.Empty
		Dim strInputFileName As String

		strInputFileName = System.IO.Path.GetFileNameWithoutExtension(strFilePath)

		Select Case eResultType
			Case ePeptideHitResultType.Sequest, ePeptideHitResultType.Inspect, ePeptideHitResultType.MSGFDB
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
	Public Shared Function AutoDetermineResultType(ByVal strFilePath As String) As ePeptideHitResultType

		Dim strFilePathLCase As String

		Dim strHeaderLine As String

		Dim eResultType As ePeptideHitResultType
		eResultType = ePeptideHitResultType.Unknown

		strFilePathLCase = strFilePath.ToLower()

		Try
			If strFilePathLCase.EndsWith("_xt.txt") Then
				eResultType = ePeptideHitResultType.XTandem
			Else
				If strFilePathLCase.EndsWith("_msgfdb_syn.txt") OrElse strFilePathLCase.EndsWith("_msgfdb_fht.txt") Then
					eResultType = ePeptideHitResultType.MSGFDB

				ElseIf strFilePathLCase.EndsWith("_inspect_syn.txt") OrElse strFilePathLCase.EndsWith("_inspect_fht.txt") Then
					eResultType = ePeptideHitResultType.Inspect

				ElseIf strFilePathLCase.EndsWith("_syn.txt") OrElse strFilePathLCase.EndsWith("_fht.txt") Then
					' Open the file and read the header line to determine if this is a Sequest file, Inspect file, MSGFDB, or something else

					If Not System.IO.File.Exists(strFilePath) Then
						' File doesn't exist; assume Sequest
						eResultType = ePeptideHitResultType.Sequest
					Else

						Try
							Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.ReadWrite))

								If srInFile.Peek() >= 0 Then
									strHeaderLine = srInFile.ReadLine().ToLower

									If strHeaderLine.Contains(clsPHRPParserInspect.DATA_COLUMN_MQScore.ToLower) AndAlso _
									   strHeaderLine.Contains(clsPHRPParserInspect.DATA_COLUMN_TotalPRMScore.ToLower) Then
										eResultType = ePeptideHitResultType.Inspect

									ElseIf strHeaderLine.Contains(clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore) AndAlso _
									   strHeaderLine.Contains(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore) Then
										eResultType = ePeptideHitResultType.MSGFDB

									ElseIf strHeaderLine.Contains(clsPHRPParserSequest.DATA_COLUMN_XCorr.ToLower) AndAlso _
									   strHeaderLine.Contains(clsPHRPParserSequest.DATA_COLUMN_DelCn.ToLower) Then
										eResultType = ePeptideHitResultType.Sequest

									End If
								End If

							End Using

						Catch ex As Exception
							' Error reading file; assume Sequest
							eResultType = ePeptideHitResultType.Sequest
						End Try

					End If

				End If
			End If

		Catch ex As Exception
			Throw ex
		End Try

		If eResultType = ePeptideHitResultType.Unknown Then
			Dim strFilePathTrimmed As String = String.Empty
			If AutoTrimExtraSuffix(strFilePath, strFilePathTrimmed) Then
				eResultType = AutoDetermineResultType(strFilePathTrimmed)
			End If
		End If

		Return eResultType

	End Function

	Protected Shared Function AutoTrimExtraSuffix(ByVal strFilePath As String, ByRef strFilePathTrimmed As String) As Boolean

		' Check whether strfilePathLCase ends in other known PHRP extensions
		Dim lstExtraSuffixes As List(Of String)
		lstExtraSuffixes = New List(Of String)

		lstExtraSuffixes.Add("_ResultToSeqMap.txt")
		lstExtraSuffixes.Add("_SeqToProteinMap.txt")
		lstExtraSuffixes.Add("_SeqInfo.txt")
		lstExtraSuffixes.Add("_MSGF.txt")

		For Each strSuffix As String In lstExtraSuffixes
			If strFilePath.ToLower().EndsWith(strSuffix.ToLower()) Then

				strFilePathTrimmed = strFilePath.Substring(0, strFilePath.Length - strSuffix.Length) & ".txt"
				Return True

				Exit For
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
	Protected Function ConvertModsToNumericMods(ByVal strPeptide As String, ByRef strPeptideWithNumericMods As String, ByRef lstPeptideMods As List(Of clsAminoAcidModInfo)) As Boolean

		Static sbNewPeptide As New System.Text.StringBuilder

		Dim intPeptideLength As Integer
		Dim chMostRecentResidue As Char
		Dim intResidueLocInPeptide As Integer = 0
		Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

		Dim intIndex As Integer
		Dim intIndexStart As Integer
		Dim intIndexEnd As Integer

		Try
			lstPeptideMods.Clear()

			If mDynamicMods.Count = 0 AndAlso mStaticMods.Count = 0 Then
				' No mods are defined; simply update strPeptideWithNumericMods to be strPeptide
				strPeptideWithNumericMods = strPeptide
				Return True
			End If

			strPeptideWithNumericMods = String.Empty
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
					If Char.IsLetter(strPeptide.Chars(intIndex)) Then
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
								AddStaticModIfPresent(mStaticMods, clsPHRPReader.N_TERMINAL_PEPTIDE_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus, sbNewPeptide, lstPeptideMods)

								If strPeptide.StartsWith(clsPHRPReader.PROTEIN_TERMINUS_SYMBOL_PHRP) Then
									' We're at the N-terminus of the protein
									' Possibly add a static N-terminal protein mod
									AddStaticModIfPresent(mStaticMods, clsPHRPReader.N_TERMINAL_PROTEIN_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus, sbNewPeptide, lstPeptideMods)
								End If
							End If
						End If

					Else
						' Not a letter; see if it is present in mDynamicMods
						AddDynamicModIfPresent(mDynamicMods, chMostRecentResidue, strPeptide.Chars(intIndex), intResidueLocInPeptide, eResidueTerminusState, sbNewPeptide, lstPeptideMods)
					End If

					If intIndex = intIndexEnd AndAlso mStaticMods.Count > 0 Then
						' Possibly add a static C-terminal peptide mod
						AddStaticModIfPresent(mStaticMods, clsPHRPReader.C_TERMINAL_PEPTIDE_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus, sbNewPeptide, lstPeptideMods)

						If strPeptide.EndsWith(clsPHRPReader.PROTEIN_TERMINUS_SYMBOL_PHRP) Then
							' We're at the C-terminus of the protein
							' Possibly add a static C-terminal protein mod
							AddStaticModIfPresent(mStaticMods, clsPHRPReader.C_TERMINAL_PROTEIN_SYMBOL_DMS, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus, sbNewPeptide, lstPeptideMods)
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

	Protected Sub AddDynamicModIfPresent(ByRef objMods As SortedDictionary(Of Char, clsModificationDefinition), _
	  ByVal chResidue As Char, _
	  ByVal chModSymbol As Char, _
	  ByVal ResidueLocInPeptide As Integer, _
	  ByVal ResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
	  ByRef sbNewPeptide As System.Text.StringBuilder, _
	  ByRef lstPeptideMods As List(Of clsAminoAcidModInfo))

		Dim objModDef As clsModificationDefinition = Nothing

		If objMods.TryGetValue(chModSymbol, objModDef) Then
			' Mod mass found for dynamic mod symbol; append the mod
			sbNewPeptide.Append(clsPHRPParser.NumToStringPlusMinus(objModDef.ModificationMass, 4))
			lstPeptideMods.Add(New clsAminoAcidModInfo(chResidue, ResidueLocInPeptide, ResidueTerminusState, objModDef))
		End If

	End Sub

	Protected Sub AddStaticModIfPresent(ByRef objMods As SortedDictionary(Of String, List(Of clsModificationDefinition)), _
	   ByVal chResidue As Char, _
	   ByVal ResidueLocInPeptide As Integer, _
	   ByVal ResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
	   ByRef sbNewPeptide As System.Text.StringBuilder, _
	   ByRef lstPeptideMods As List(Of clsAminoAcidModInfo))

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
	Protected Function GetCollisionMode(ByVal strScanTypeName As String) As String
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
	''' Returns the default first-hits file name for the given PeptideHit result type
	''' </summary>
	''' <param name="eResultType"></param>
	''' <param name="strDatasetName"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetPHRPFirstHitsFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strPHRPResultsFileName As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest: _fht.txt
				strPHRPResultsFileName = clsPHRPParserSequest.GetPHRPFirstHitsFileName(strDatasetName)

			Case ePeptideHitResultType.XTandem
				' X!Tandem does not have a first-hits file; strPHRPResultsFileName will be an empty string
				strPHRPResultsFileName = clsPHRPParserXTandem.GetPHRPFirstHitsFileName(strDatasetName)

			Case ePeptideHitResultType.Inspect
				' Inspect: _inspect_fht.txt
				strPHRPResultsFileName = clsPHRPParserInspect.GetPHRPFirstHitsFileName(strDatasetName)

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB: _msgfdb_fht.txt
				strPHRPResultsFileName = clsPHRPParserMSGFDB.GetPHRPFirstHitsFileName(strDatasetName)

		End Select

		Return strPHRPResultsFileName

	End Function

	''' <summary>
	''' Returns the default ModSummary file name for the given PeptideHit result type
	''' </summary>
	''' <param name="eResultType"></param>
	''' <param name="strDatasetName"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetPHRPModSummaryFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strPHRPModSummaryFileName As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest: _syn.txt
				strPHRPModSummaryFileName = clsPHRPParserSequest.GetPHRPModSummaryFileName(strDatasetName)

			Case ePeptideHitResultType.XTandem
				' X!Tandem: _xt.txt
				strPHRPModSummaryFileName = clsPHRPParserXTandem.GetPHRPModSummaryFileName(strDatasetName)

			Case ePeptideHitResultType.Inspect
				' Inspect: _inspect_syn.txt
				strPHRPModSummaryFileName = clsPHRPParserInspect.GetPHRPModSummaryFileName(strDatasetName)

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB: _msgfdb_syn.txt
				strPHRPModSummaryFileName = clsPHRPParserMSGFDB.GetPHRPModSummaryFileName(strDatasetName)

		End Select

		Return strPHRPModSummaryFileName

	End Function

	''' <summary>
	''' Returns the default Synopsis file name for the given PeptideHit result type
	''' </summary>
	''' <param name="eResultType"></param>
	''' <param name="strDatasetName"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetPHRPSynopsisFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strPHRPResultsFileName As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest: _syn.txt
				strPHRPResultsFileName = clsPHRPParserSequest.GetPHRPSynopsisFileName(strDatasetName)

			Case ePeptideHitResultType.XTandem
				' X!Tandem: _xt.txt
				strPHRPResultsFileName = clsPHRPParserXTandem.GetPHRPSynopsisFileName(strDatasetName)

			Case ePeptideHitResultType.Inspect
				' Inspect: _inspect_syn.txt
				strPHRPResultsFileName = clsPHRPParserInspect.GetPHRPSynopsisFileName(strDatasetName)

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB: _msgfdb_syn.txt
				strPHRPResultsFileName = clsPHRPParserMSGFDB.GetPHRPSynopsisFileName(strDatasetName)

		End Select

		Return strPHRPResultsFileName

	End Function

	''' <summary>
	''' Returns the default ResultToSeq Map file name for the given PeptideHit result type
	''' </summary>
	''' <param name="eResultType"></param>
	''' <param name="strDatasetName"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strPHRPResultsFileName As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest: _syn.txt
				strPHRPResultsFileName = clsPHRPParserSequest.GetPHRPResultToSeqMapFileName(strDatasetName)

			Case ePeptideHitResultType.XTandem
				' X!Tandem: _xt.txt
				strPHRPResultsFileName = clsPHRPParserXTandem.GetPHRPResultToSeqMapFileName(strDatasetName)

			Case ePeptideHitResultType.Inspect
				' Inspect: _inspect_syn.txt
				strPHRPResultsFileName = clsPHRPParserInspect.GetPHRPResultToSeqMapFileName(strDatasetName)

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB: _msgfdb_syn.txt
				strPHRPResultsFileName = clsPHRPParserMSGFDB.GetPHRPResultToSeqMapFileName(strDatasetName)

		End Select

		Return strPHRPResultsFileName

	End Function

	''' <summary>
	''' Returns the default SeqInfo file name for the given PeptideHit result type
	''' </summary>
	''' <param name="eResultType"></param>
	''' <param name="strDatasetName"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetPHRPSeqInfoFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strSeqInfoFilename As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest: _syn.txt
				strSeqInfoFilename = clsPHRPParserSequest.GetPHRPSeqInfoFileName(strDatasetName)

			Case ePeptideHitResultType.XTandem
				' X!Tandem: _xt.txt
				strSeqInfoFilename = clsPHRPParserXTandem.GetPHRPSeqInfoFileName(strDatasetName)

			Case ePeptideHitResultType.Inspect
				' Inspect: _inspect_syn.txt
				strSeqInfoFilename = clsPHRPParserInspect.GetPHRPSeqInfoFileName(strDatasetName)

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB: _msgfdb_syn.txt
				strSeqInfoFilename = clsPHRPParserMSGFDB.GetPHRPSeqInfoFileName(strDatasetName)

		End Select

		Return strSeqInfoFilename

	End Function

	''' <summary>
	''' Returns the default SeqToProtein Map file name for the given PeptideHit result type
	''' </summary>
	''' <param name="eResultType"></param>
	''' <param name="strDatasetName"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strSeqToProteinMapFileName As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest: _syn.txt
				strSeqToProteinMapFileName = clsPHRPParserSequest.GetPHRPSeqToProteinMapFileName(strDatasetName)

			Case ePeptideHitResultType.XTandem
				' X!Tandem: _xt.txt
				strSeqToProteinMapFileName = clsPHRPParserXTandem.GetPHRPSeqToProteinMapFileName(strDatasetName)

			Case ePeptideHitResultType.Inspect
				' Inspect: _inspect_syn.txt
				strSeqToProteinMapFileName = clsPHRPParserInspect.GetPHRPSeqToProteinMapFileName(strDatasetName)

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB: _msgfdb_syn.txt
				strSeqToProteinMapFileName = clsPHRPParserMSGFDB.GetPHRPSeqToProteinMapFileName(strDatasetName)

		End Select

		Return strSeqToProteinMapFileName

	End Function

	Protected Sub HandleException(ByVal strBaseMessage As String, ByVal ex As System.Exception)
		If String.IsNullOrEmpty(strBaseMessage) Then
			strBaseMessage = "Error"
		End If

		ReportError(strBaseMessage & ": " & ex.Message)

	End Sub

	''' <summary>
	''' Examines the string to determine if it is numeric
	''' </summary>
	''' <param name="strData"></param>
	''' <returns>True if a number, otherwise false</returns>
	Public Shared Function IsNumber(ByVal strData As String) As Boolean

		If Double.TryParse(strData, 0) Then
			Return True
		ElseIf Integer.TryParse(strData, 0) Then
			Return True
		End If

		Return False

	End Function

	''' <summary>
	''' Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
	''' </summary>
	''' <returns>The text in the specified column; an empty string if the specific column name is not recognized</returns>
	''' <remarks></remarks>
	Public Shared Function LookupColumnValue(ByRef strColumns() As String, _
	  ByVal strColumnName As String, _
	  ByRef objColumnHeaders As SortedDictionary(Of String, Integer)) As String

		Return LookupColumnValue(strColumns, strColumnName, objColumnHeaders, String.Empty)
	End Function

	''' <summary>
	''' Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
	''' </summary>
	''' <returns>The text in the specified column; strValueIfMissing if the specific column name is not recognized</returns>
	''' <remarks></remarks>
	Public Shared Function LookupColumnValue(ByRef strColumns() As String, _
	  ByVal strColumnName As String, _
	  ByRef objColumnHeaders As SortedDictionary(Of String, Integer), _
	  ByVal strValueIfMissing As String) As String

		Dim intColIndex As Integer

		If Not strColumns Is Nothing Then
			If objColumnHeaders.TryGetValue(strColumnName, intColIndex) Then
				If intColIndex >= 0 AndAlso intColIndex < strColumns.Length Then
					If String.IsNullOrWhiteSpace(strColumns(intColIndex)) Then
						Return String.Empty
					Else
						Return strColumns(intColIndex)
					End If
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
	Public Shared Function LookupColumnValue(ByRef strColumns() As String, _
	  ByVal strColumnName As String, _
	  ByRef objColumnHeaders As SortedDictionary(Of String, Integer), _
	  ByVal ValueIfMissing As Integer) As Integer

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
	Public Shared Function LookupColumnValue(ByRef strColumns() As String, _
	  ByVal strColumnName As String, _
	  ByRef objColumnHeaders As SortedDictionary(Of String, Integer), _
	  ByVal ValueIfMissing As Double) As Double

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
	Public Shared Sub ParseColumnHeaders(ByVal strColumns() As String, ByRef objColumnHeaders As SortedDictionary(Of String, Integer))

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
	''' <remarks></remarks>
	Public Function MoveNext() As Boolean

		Dim strLineIn As String = String.Empty
		Dim strSplitLine() As String

		Dim blnSkipLine As Boolean = False
		Dim blnSuccess As Boolean = False
		Dim blnMatchFound As Boolean = False

		Dim blnScanStatsValid As Boolean
		Dim objScanStatsInfo As clsScanStatsInfo = Nothing

		Dim blnExtendedScanStatsValid As Boolean
		Dim objExtendedScanStatsInfo As clsScanStatsExInfo = Nothing

		Dim dblMonoisotopicPrecursorMass As Double

		If mCachedLineAvailable Then
			strLineIn = mCachedLine
			mCachedLineAvailable = False
			blnSuccess = True
		ElseIf mSourceFile.Peek > -1 Then
			strLineIn = mSourceFile.ReadLine()
			mSourceFileLinesRead += 1
			blnSuccess = True
		Else
			mCanRead = False
		End If

		If blnSuccess Then
			blnSkipLine = False

			If Not String.IsNullOrEmpty(strLineIn) Then
				strSplitLine = strLineIn.Split(ControlChars.Tab)

				If Not mHeaderLineParsed Then
					If Not clsPHRPReader.IsNumber(strSplitLine(0)) Then
						' Parse the header line to confirm the column ordering
						mPHRPParser.ParseColumnHeaders(strSplitLine)

						mHeaderLineParsed = True
						Return Me.MoveNext()
					End If

					mHeaderLineParsed = True
				End If

				If Not blnSkipLine Then

					mPSMCurrent = New clsPSM()
					blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mSourceFileLinesRead, mPSMCurrent)

					If blnSuccess Then
						blnMatchFound = True

						' The PHRPParser will update .PeptideWithNumericMods if the _SeqInfo.txt file is loaded
						' If it wasn't loaded, then this class can update .PeptideWithNumericMods and .PeptideMods 
						' by inferring the mods using mDynamicMods and mStaticMods (which were populated using the PHRP ModSummary file)
						If mLoadModsAndSeqInfo AndAlso String.IsNullOrEmpty(mPSMCurrent.PeptideWithNumericMods) Then

							' Markup the peptide with the dynamic and static mods
							Dim strPeptideWithMods As String = String.Empty
							Dim lstPeptideMods As List(Of clsAminoAcidModInfo)
							lstPeptideMods = New List(Of clsAminoAcidModInfo)

							blnSuccess = ConvertModsToNumericMods(mPSMCurrent.Peptide.Trim, strPeptideWithMods, lstPeptideMods)
							If blnSuccess Then
								Dim dblTotalModMass As Double = 0

								mPSMCurrent.PeptideWithNumericMods = strPeptideWithMods
								mPSMCurrent.ClearModifiedResidues()
								For Each objModEntry In lstPeptideMods
									mPSMCurrent.AddModifiedResidue(objModEntry)
									dblTotalModMass += objModEntry.ModDefinition.ModificationMass
								Next

								If mPSMCurrent.PeptideMonoisotopicMass = 0 Then
									mPSMCurrent.PeptideMonoisotopicMass = mPeptideMassCalculator.ComputeSequenceMass(mPSMCurrent.PeptideCleanSequence) + dblTotalModMass
								End If

							End If
						End If

						blnScanStatsValid = TryGetScanStats(mPSMCurrent.ScanNumber, objScanStatsInfo)
						blnExtendedScanStatsValid = TryGetExtendedScanStats(mPSMCurrent.ScanNumber, objExtendedScanStatsInfo)

						If blnScanStatsValid Then
							' Update the elution time
							mPSMCurrent.ElutionTimeMinutes = objScanStatsInfo.ScanTimeMinutes
						End If

						If String.IsNullOrEmpty(mPSMCurrent.CollisionMode) OrElse mPSMCurrent.CollisionMode = clsPSM.UNKNOWN_COLLISION_MODE Then

							' Determine the ScanTypeName using the the ScanStats or ExtendedScanStats info
							If blnScanStatsValid AndAlso Not String.IsNullOrEmpty(objScanStatsInfo.ScanTypeName) Then
								mPSMCurrent.CollisionMode = GetCollisionMode(objScanStatsInfo.ScanTypeName)
							End If

							If String.IsNullOrEmpty(mPSMCurrent.CollisionMode) AndAlso blnExtendedScanStatsValid Then
								' Scan type still not determined, but Extended Scan Stats data is available
								If Not String.IsNullOrEmpty(objExtendedScanStatsInfo.CollisionMode) Then
									' Check for Collision mode being "0"
									' This is often the case for the first scan in a Thermo .Raw file
									If objExtendedScanStatsInfo.CollisionMode <> "0" Then
										mPSMCurrent.CollisionMode = objExtendedScanStatsInfo.CollisionMode.ToUpper()
									End If
								End If
							End If
						End If

						If mPHRPParser.PeptideHitResultType = ePeptideHitResultType.Sequest OrElse mPHRPParser.PeptideHitResultType = ePeptideHitResultType.XTandem Then

							dblMonoisotopicPrecursorMass = 0
							If blnExtendedScanStatsValid Then
								' Try to extract out the precursor m/z value from the "Scan Filter Text" field
								Dim dblParentIonMZ As Double
								Dim intMSLevel As Integer
								Dim strCollisionMode As String = String.Empty

								If ThermoRawFileReaderDLL.FinniganFileIO.XRawFileIO.ExtractParentIonMZFromFilterText(objExtendedScanStatsInfo.ScanFilterText, dblParentIonMZ, intMSLevel, strCollisionMode) Then
									If dblParentIonMZ > 0 Then
										dblMonoisotopicPrecursorMass = clsPeptideMassCalculator.ConvoluteMass(dblParentIonMZ, mPSMCurrent.Charge, 0)
									End If
								End If

							End If

							If dblMonoisotopicPrecursorMass = 0 AndAlso blnExtendedScanStatsValid Then
								If objExtendedScanStatsInfo.MonoisotopicMZ > 0 Then
									' Determine the precursor m/z value using the Monoisotopic m/z value reported by the instrument
									dblMonoisotopicPrecursorMass = clsPeptideMassCalculator.ConvoluteMass(objExtendedScanStatsInfo.MonoisotopicMZ, mPSMCurrent.Charge, 0)
								End If
							End If

							If dblMonoisotopicPrecursorMass > 0 Then
								mPSMCurrent.PrecursorNeutralMass = dblMonoisotopicPrecursorMass
							End If

						End If

						If Not mMSGFCachedResults Is Nothing AndAlso mMSGFCachedResults.Count > 0 Then
							Dim strMSGFSpecProb As String = String.Empty
							Dim dblSpecProb As Double

							If mMSGFCachedResults.TryGetValue(mPSMCurrent.ResultID, strMSGFSpecProb) Then
								mPSMCurrent.MSGFSpecProb = strMSGFSpecProb
								If strMSGFSpecProb.Length > 13 Then
									' Attempt to shorten the SpecProb value
									If Double.TryParse(strMSGFSpecProb, dblSpecProb) Then
										mPSMCurrent.MSGFSpecProb = dblSpecProb.ToString("0.0000000E-00")
									End If
								End If
							End If
						End If

						If mSkipDuplicatePSMs Then

							' Read the next line and check whether it's the same hit, but a different protein
							Dim blnReadNext As Boolean = True
							Do While blnReadNext AndAlso mSourceFile.Peek > -1
								strLineIn = mSourceFile.ReadLine()
								mSourceFileLinesRead += 1

								If Not String.IsNullOrEmpty(strLineIn) Then

									Dim objNewPSM As New clsPSM
									blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mSourceFileLinesRead, objNewPSM)

									' Check for duplicate lines
									' If this line is a duplicate of the previous line, then skip it
									' This happens in Sequest _syn.txt files where the line is repeated for all protein matches
									With mPSMCurrent
										If .ScanNumber = objNewPSM.ScanNumber AndAlso _
										   .Charge = objNewPSM.Charge AndAlso _
										   .Peptide = objNewPSM.Peptide Then

											' Yes, this is a duplicate
											' Update the protein list
											For Each strProtein As String In objNewPSM.Proteins
												.AddProtein(strProtein)
											Next
										Else
											blnReadNext = False
											mCachedLine = String.Copy(strLineIn)
											mCachedLineAvailable = True
										End If
									End With
								End If
							Loop

						End If

						blnSuccess = True

					End If

				End If
			End If

		End If

		Return blnMatchFound

	End Function

	Protected Function ReadAndCacheMSGFData() As Boolean

		Dim strMSGFFilePath As String
		Dim blnSuccess As Boolean = False

		Try
			strMSGFFilePath = System.IO.Path.GetFileNameWithoutExtension(mInputFilePath) & MSGF_RESULT_FILENAME_SUFFIX
			strMSGFFilePath = System.IO.Path.Combine(mInputFolderPath, strMSGFFilePath)

			blnSuccess = ReadAndCacheMSGFData(strMSGFFilePath)

		Catch ex As Exception
			HandleException("Exception determining MSGF file path", ex)
		End Try

		Return blnSuccess

	End Function

	Protected Function ReadAndCacheMSGFData(ByVal strMSGFFilePath As String) As Boolean
		Dim blnSuccess As Boolean = False

		Try
			If System.IO.File.Exists(strMSGFFilePath) Then

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
	Protected Function ReadModSummaryFile(ByVal strModSummaryFilePath As String, _
	  ByRef objDynamicMods As SortedDictionary(Of Char, clsModificationDefinition), _
	  ByRef objStaticMods As SortedDictionary(Of String, List(Of clsModificationDefinition))) As Boolean

		Dim objModSummaryReader As clsPHRPModSummaryReader

		Dim lstModDefs As List(Of clsModificationDefinition) = Nothing

		Dim strModMass As String
		Dim blnSuccess As Boolean

		Try
			If String.IsNullOrEmpty(strModSummaryFilePath) Then
				ReportError("ModSummaryFile path is empty; unable to continue")
				Return False
			End If

			ShowMessage("Reading the PHRP ModSummary file")

			' Clear objDynamicMods and objStaticMods (should have been instantiated by the calling function)
			objDynamicMods.Clear()
			objStaticMods.Clear()

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

	Protected Function ReadScanStatsData() As Boolean

		Dim strScanStatsFilePath As String
		Dim blnSuccess As Boolean = False

		Try
			strScanStatsFilePath = mDatasetName & SCAN_STATS_FILENAME_SUFFIX
			strScanStatsFilePath = System.IO.Path.Combine(mInputFolderPath, strScanStatsFilePath)

			blnSuccess = ReadScanStatsData(strScanStatsFilePath)

		Catch ex As Exception
			HandleException("Exception determining Scan Stats file path", ex)
		End Try

		Return blnSuccess

	End Function

	Protected Function ReadScanStatsData(ByVal strScanStatsFilePath As String) As Boolean

		Dim blnSuccess As Boolean = False

		Try
			If System.IO.File.Exists(strScanStatsFilePath) Then

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

	Protected Function ReadExtendedScanStatsData() As Boolean

		Dim strExtendedScanStatsFilePath As String
		Dim blnSuccess As Boolean = False

		Try
			strExtendedScanStatsFilePath = mDatasetName & EXTENDED_SCAN_STATS_FILENAME_SUFFIX
			strExtendedScanStatsFilePath = System.IO.Path.Combine(mInputFolderPath, strExtendedScanStatsFilePath)

			blnSuccess = ReadExtendedScanStatsData(strExtendedScanStatsFilePath)

		Catch ex As Exception
			HandleException("Exception determining Scan Stats file path", ex)
		End Try

		Return blnSuccess

	End Function

	Protected Function ReadExtendedScanStatsData(ByVal strExtendedScanStatsFilePath As String) As Boolean

		Dim blnSuccess As Boolean = False

		Try
			If System.IO.File.Exists(strExtendedScanStatsFilePath) Then

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


	Protected Sub ReportError(ByVal strErrorMessage As String)
		mErrorMessage = strErrorMessage
		If mEchoMessagesToConsole Then Console.WriteLine(strErrorMessage)
		mErrorMessages.Add(strErrorMessage)

		RaiseEvent ErrorEvent(strErrorMessage)
	End Sub

	Protected Sub ReportWarning(ByVal strWarningMessage As String)
		If mEchoMessagesToConsole Then Console.WriteLine(strWarningMessage)
		mWarningMessages.Add(strWarningMessage)

		RaiseEvent WarningEvent(strWarningMessage)
	End Sub

	Private Sub SetLocalErrorCode(ByVal eNewErrorCode As ePHRPReaderErrorCodes)
		SetLocalErrorCode(eNewErrorCode, False)
	End Sub

	Private Sub SetLocalErrorCode(ByVal eNewErrorCode As ePHRPReaderErrorCodes, ByVal blnLeaveExistingErrorCodeUnchanged As Boolean)

		If blnLeaveExistingErrorCodeUnchanged AndAlso mLocalErrorCode <> ePHRPReaderErrorCodes.NoError Then
			' An error code is already defined; do not change it
		Else
			mLocalErrorCode = eNewErrorCode
		End If

	End Sub

	Protected Sub ShowMessage(ByVal strMessage As String)
		If mEchoMessagesToConsole Then Console.WriteLine(strMessage)
		RaiseEvent MessageEvent(strMessage)
	End Sub

	Protected Function TryGetScanStats(ByVal intScanNumber As Integer, ByRef objScanStatsInfo As clsScanStatsInfo) As Boolean
		If Not mScanStats Is Nothing AndAlso mScanStats.Count > 0 Then
			If mScanStats.TryGetValue(intScanNumber, objScanStatsInfo) Then
				Return True
			End If
		End If
		Return False
	End Function

	Protected Function TryGetExtendedScanStats(ByVal intScanNumber As Integer, ByRef objExtendedScanStatsInfo As clsScanStatsExInfo) As Boolean
		If Not mScanStatsEx Is Nothing AndAlso mScanStats.Count > 0 Then
			If mScanStatsEx.TryGetValue(intScanNumber, objExtendedScanStatsInfo) Then
				Return True
			End If
		End If
		Return False
	End Function

	Private Function ValidateInputFiles(ByVal strInputFilePath As String, _
	  ByRef eResultType As ePeptideHitResultType, _
	  ByRef strModSummaryFilePath As String) As Boolean

		Dim fiFileInfo As System.IO.FileInfo

		fiFileInfo = New System.IO.FileInfo(strInputFilePath)
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
			ReportError("Error: Unable to auto-determine the dataset name frome the input file name: " & strInputFilePath)
			SetLocalErrorCode(ePHRPReaderErrorCodes.InputFileFormatNotRecognized)
			Return False
		End If

		If mLoadModsAndSeqInfo Then
			strModSummaryFilePath = GetPHRPModSummaryFileName(eResultType, mDatasetName)
			strModSummaryFilePath = System.IO.Path.Combine(fiFileInfo.DirectoryName, strModSummaryFilePath)
			If Not ValidateRequiredFileExists("ModSummary file", strModSummaryFilePath) Then
				SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound)
				Return False
			End If
		Else
			strModSummaryFilePath = String.Empty
		End If

		Return True

	End Function

	Private Function ValidateRequiredFileExists(ByVal strFileDescription As String, ByVal strFilePath As String) As Boolean
		Return ValidateRequiredFileExists(strFileDescription, strFilePath, True)
	End Function

	Private Function ValidateRequiredFileExists(ByVal strFileDescription As String, ByVal strFilePath As String, ByVal blnReportErrors As Boolean) As Boolean

		If strFilePath Is Nothing OrElse strFilePath.Length = 0 Then
			If blnReportErrors Then
				ReportError(strFileDescription & " is not defined")
			End If
			Return False

		ElseIf Not System.IO.File.Exists(strFilePath) Then
			If blnReportErrors Then
				ReportError(strFileDescription & " not found: " & strFilePath)
			End If
			Return False

		End If

		Return True

	End Function

	Private Sub mPHRPParser_ErrorEvent(ByVal strErrorMessage As String) Handles mPHRPParser.ErrorEvent
		ReportError(strErrorMessage)
	End Sub

	Private Sub mPHRPParser_MessageEvent(ByVal strMessage As String) Handles mPHRPParser.MessageEvent
		ShowMessage(strMessage)
	End Sub

	Private Sub mPHRPParser_WarningEvent(ByVal strWarningMessage As String) Handles mPHRPParser.WarningEvent
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
		' Do not change this code.  Put cleanup code in Dispose(ByVal disposing As Boolean) above.
		Dispose(True)
		GC.SuppressFinalize(Me)
	End Sub
#End Region

End Class
