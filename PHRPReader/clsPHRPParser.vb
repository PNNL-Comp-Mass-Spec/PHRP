'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class is the base class for classes used to parse PHRP data lines
' It must be derived by a sub-class customized for the specific analysis tool (Sequest, X!Tandem, Inspect, etc.)
'
'*********************************************************************************************************

Option Strict On

Imports System.IO

Public MustInherit Class clsPHRPParser

#Region "Structures"
	Protected Structure udtAmbiguousModInfo
		Public ResidueStart As Integer
		Public ResidueEnd As Integer
		Public ModMassString As String
	End Structure
#End Region

#Region "Module variables"

	Protected mDatasetName As String
	Protected mInputFilePath As String
	Protected mInputFolderPath As String
	Protected mInitialized As Boolean

	' Column headers in the synopsis file and first hits file
	Protected mColumnHeaders As SortedDictionary(Of String, Integer)

	Protected mErrorMessage As String = String.Empty

	Protected mCleavageStateCalculator As clsPeptideCleavageStateCalculator
	Protected mPeptideMassCalculator As clsPeptideMassCalculator

	Protected mPeptideHitResultType As clsPHRPReader.ePeptideHitResultType

	Protected mModInfo As List(Of clsModificationDefinition)
	Protected mModInfoLoaded As Boolean

	Protected mResultToSeqMap As SortedList(Of Integer, Integer)
	Protected mSeqInfo As SortedList(Of Integer, clsSeqInfo)
	Protected mSeqToProteinMap As SortedList(Of Integer, List(Of clsProteinInfo))
	Protected mPepToProteinMap As Dictionary(Of String, clsPepToProteinMapInfo)

	' This List tracks the Protein Names for each ResultID
	Protected mResultIDToProteins As SortedList(Of Integer, SortedSet(Of String))

	Protected mErrorMessages As List(Of String)
	Protected mWarningMessages As List(Of String)

#End Region

#Region "Events"
	Public Event MessageEvent(ByVal strMessage As String)
	Public Event ErrorEvent(ByVal strErrorMessage As String)
	Public Event WarningEvent(ByVal strWarningMessage As String)
#End Region

#Region "Properties"

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
	''' Input file path
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property InputFilePath As String
		Get
			Return mInputFilePath
		End Get
	End Property

	''' <summary>
	''' Input folder path
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property InputFolderPath As String
		Get
			Return mInputFolderPath
		End Get
	End Property

	''' <summary>
	''' Peptide hit result type; Sequest, XTandem, Inspect, or MSGFDB
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property PeptideHitResultType As clsPHRPReader.ePeptideHitResultType
		Get
			Return mPeptideHitResultType
		End Get
	End Property

	Public ReadOnly Property PepToProteinMap() As Dictionary(Of String, clsPepToProteinMapInfo)
		Get
			Return mPepToProteinMap
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
			Return mResultToSeqMap
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
			Return mSeqInfo
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
			Return mSeqToProteinMap
		End Get
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
	''' Initialize the parser for the given dataset and input file
	''' </summary>
	''' <param name="strDatasetName">Dataset Name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
	''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
	''' <remarks></remarks>
	Protected Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, ByVal blnLoadModsAndSeqInfo As Boolean)

		mErrorMessages = New List(Of String)
		mWarningMessages = New List(Of String)

		mDatasetName = strDatasetName
		mPeptideHitResultType = ePeptideHitResultType

		Dim fiFileInfo As FileInfo = New FileInfo(strInputFilePath)
		mInputFilePath = fiFileInfo.FullName
		mInputFolderPath = fiFileInfo.DirectoryName

		mErrorMessage = String.Empty

		Dim blnIsSynopsisFile As Boolean
		If fiFileInfo.Name.ToLower() = clsPHRPReader.GetPHRPSynopsisFileName(mPeptideHitResultType, mDatasetName).ToLower() Then
			blnIsSynopsisFile = True
		Else
			blnIsSynopsisFile = False
		End If

		' Initialize the column mapping object
		' Using a case-insensitive comparer
		mColumnHeaders = New SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

		mCleavageStateCalculator = New clsPeptideCleavageStateCalculator()
		mPeptideMassCalculator = New clsPeptideMassCalculator()

		' Initialize the tracking lists
		mResultToSeqMap = New SortedList(Of Integer, Integer)
		mSeqInfo = New SortedList(Of Integer, clsSeqInfo)
		mSeqToProteinMap = New SortedList(Of Integer, List(Of clsProteinInfo))
		mPepToProteinMap = New Dictionary(Of String, clsPepToProteinMapInfo)()

		mResultIDToProteins = New SortedList(Of Integer, SortedSet(Of String))

		If blnLoadModsAndSeqInfo Then
			' Read the ModSummary file (if it exists)
			mModInfoLoaded = LoadModSummary()
		Else
			mModInfoLoaded = False
		End If

		If blnLoadModsAndSeqInfo Then
			' Read the ResultToSeqMapInfo (if the files exist)			
			If blnIsSynopsisFile Then
				' Assume the files exist
				LoadSeqInfo()
			Else
				' Only continue if the fht versions exists

				Dim strResultToSeqMapFilePath As String = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName)
				Dim blnSeqInfoLoaded As Boolean = False

				If Not String.IsNullOrEmpty(strResultToSeqMapFilePath) Then
					strResultToSeqMapFilePath = Path.Combine(mInputFolderPath, strResultToSeqMapFilePath)
					strResultToSeqMapFilePath = clsPHRPReader.AutoSwitchToFHTIfRequired(strResultToSeqMapFilePath, mInputFilePath)

					If File.Exists(strResultToSeqMapFilePath) Then
						blnSeqInfoLoaded = LoadSeqInfo()
					End If
				End If

				If Not blnSeqInfoLoaded Then
					If String.IsNullOrEmpty(strResultToSeqMapFilePath) Then
						ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file and unable to determine the ResultToSeqMapFilename using clsPHRPReader.GetPHRPResultToSeqMapFileName()")
					Else
						ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file but the ResultToSeqMap file does not exist: " & strResultToSeqMapFilePath)
					End If
				End If

			End If
		End If

		' The following will be overridden by a derived form of this class
		DefineColumnHeaders()

	End Sub

	''' <summary>
	''' Returns the appropriate PHRPParser class based on the input file name; assumes blnLoadModsAndSeqInfo=True
	''' </summary>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from strInputFilePath</remarks>
	Public Shared Function GetParser(ByVal strInputFilePath As String) As clsPHRPParser
		Return GetParser(strInputFilePath, True)
	End Function

	''' <summary>
	''' Returns the appropriate PHRPParser class based on the input file name
	''' </summary>
	''' <param name="strInputFilePath">Input file path</param>
	''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
	''' <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from strInputFilePath</remarks>
	Public Shared Function GetParser(ByVal strInputFilePath As String, ByVal blnLoadModsAndSeqInfo As Boolean) As clsPHRPParser
		Dim ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath)

		If ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Unknown Then
			Throw New Exception("Unable to auto-determine the PeptideHitResultType for " & strInputFilePath)
		End If

		Dim strDatasetName = clsPHRPReader.AutoDetermineDatasetName(strInputFilePath)
		If String.IsNullOrEmpty(strDatasetName) Then
			Throw New Exception("Unable to auto-determine the Dataset Name for " & strInputFilePath)
		End If

		Return GetParser(strInputFilePath, strDatasetName, ePeptideHitResultType, blnLoadModsAndSeqInfo)
	End Function

	''' <summary>
	''' Returns the appropriate PHRPParser class based on the input file name
	''' </summary>
	''' <param name="strInputFilePath">Input file path</param>
	''' ''' <param name="strDatasetName">Dataset Name</param>
	''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
	''' <remarks>Throws an exception if unable to auto-determine the input file type from strInputFilePath</remarks>
	Public Shared Function GetParser(ByVal strInputFilePath As String, ByVal strDatasetName As String, ByVal blnLoadModsAndSeqInfo As Boolean) As clsPHRPParser
		Dim ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath)

		If ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Unknown Then
			Throw New Exception("Unable to auto-determine the PeptideHitResultType for " & strInputFilePath)
		End If

		Return GetParser(strInputFilePath, strDatasetName, ePeptideHitResultType, blnLoadModsAndSeqInfo)
	End Function

	''' <summary>
	''' Returns the appropriate PHRPParser class based on ePeptideHitResultType
	''' </summary>
	''' <param name="strInputFilePath">Input file path</param>
	''' <param name="strDatasetName">Dataset Name</param>
	''' <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
	''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
	''' <remarks></remarks>
	Public Shared Function GetParser(ByVal strInputFilePath As String, ByVal strDatasetName As String, ByVal ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, ByVal blnLoadModsAndSeqInfo As Boolean) As clsPHRPParser
		Select Case ePeptideHitResultType
			Case clsPHRPReader.ePeptideHitResultType.Inspect
				Return New clsPHRPParserInspect(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

			Case clsPHRPReader.ePeptideHitResultType.MSAlign
				Return New clsPHRPParserMSAlign(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

			Case clsPHRPReader.ePeptideHitResultType.MSGFDB
				Return New clsPHRPParserMSGFDB(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

			Case clsPHRPReader.ePeptideHitResultType.Sequest
				Return New clsPHRPParserSequest(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

			Case clsPHRPReader.ePeptideHitResultType.XTandem
				Return New clsPHRPParserXTandem(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

			Case Else
				Throw New Exception("Unrecognized value for PeptideHitResultType: " & ePeptideHitResultType.ToString())
		End Select

	End Function

#Region "Functions overridden by derived classes"
	Protected MustOverride Sub DefineColumnHeaders()

	''' <summary>
	''' Parse the data line read from a PHRP results file
	''' </summary>
	''' <param name="strLine">Data line</param>
	''' <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
	''' <param name="objPSM">clsPSM object (output)</param>
	''' <returns>True if success, false if an error</returns>
	Public MustOverride Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean

	''' <summary>
	''' Parses the specified parameter file
	''' Also reads the Tool_Version_Info file in the same folder (if present)
	''' </summary>
	''' <param name="strSearchEngineParamFileName">Name of the parameter file to parse (must reside in InputFolderPath)</param>
	''' <param name="objSearchEngineParams">Search engine parameters class (output)</param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public MustOverride Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
#End Region

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub AddScore(ByRef objPSM As clsPSM, ByRef strColumns() As String, ByVal strScoreColumnName As String)
		Const NOT_FOUND As String = "==SCORE_NOT_FOUND=="

		Dim strValue As String
		strValue = clsPHRPReader.LookupColumnValue(strColumns, strScoreColumnName, mColumnHeaders, NOT_FOUND)

		If strValue <> NOT_FOUND Then
			objPSM.SetScore(strScoreColumnName, strValue)
		End If

	End Sub

	''' <summary>
	''' Clear any cached error messages
	''' </summary>
	''' <remarks></remarks>
	Public Sub ClearErrors()
		mErrorMessages.Clear()
	End Sub

	''' <summary>
	''' Clear any cached warning messages
	''' </summary>
	''' <remarks></remarks>
	Public Sub ClearWarnings()
		mWarningMessages.Clear()
	End Sub

	Protected Function ConvertModsToNumericMods(ByVal strCleanSequence As String, ByRef lstModifiedResidues As List(Of clsAminoAcidModInfo)) As String
		Static sbNewPeptide As New Text.StringBuilder

		sbNewPeptide.Length = 0

		If lstModifiedResidues Is Nothing OrElse lstModifiedResidues.Count = 0 Then
			Return strCleanSequence
		End If

		For intIndex = 0 To strCleanSequence.Length - 1
			sbNewPeptide.Append(strCleanSequence.Chars(intIndex))

			For Each objModInfo In lstModifiedResidues
				If objModInfo.ResidueLocInPeptide = intIndex + 1 Then
					sbNewPeptide.Append(NumToStringPlusMinus(objModInfo.ModDefinition.ModificationMass, 4))
				End If
			Next
		Next

		Return sbNewPeptide.ToString()

	End Function

	''' <summary>
	''' Look for ambiguous mods in strSequenceWithMods
	''' For example, -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q
	''' </summary>
	''' <param name="strSequenceWithMods"></param>
	''' <returns></returns>
	''' <remarks>List of ambiguous mods, where the keys are the start residues and the values are the ambiguous mod info</remarks>
	Protected Function ExtractAmbiguousMods(ByVal strSequenceWithMods As String) As SortedList(Of Integer, udtAmbiguousModInfo)

		Dim strPrimarySequence As String = String.Empty
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty

		If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then
			strPrimarySequence = String.Copy(strSequenceWithMods)
		End If

		Dim lstAmbiguousMods = New SortedList(Of Integer, udtAmbiguousModInfo)

		Dim intResidueNumber As Integer = 0
		Dim blnParsingAmbiguousMod As Boolean = False
		Dim udtCurrentMod As udtAmbiguousModInfo

		For intCharIndex = 0 To strPrimarySequence.Length - 1
			If clsPHRPReader.IsLetterAtoZ(strPrimarySequence.Chars(intCharIndex)) Then
				' Found a letter
				intResidueNumber += 1

				If intCharIndex > 0 AndAlso strPrimarySequence.Chars(intCharIndex - 1) = "("c Then
					' Found an ambiguous mod
					If Not blnParsingAmbiguousMod Then
						blnParsingAmbiguousMod = True
						udtCurrentMod.ResidueStart = intResidueNumber
						udtCurrentMod.ResidueEnd = intResidueNumber
						udtCurrentMod.ModMassString = String.Empty
					End If
				End If

			ElseIf blnParsingAmbiguousMod Then
				' Found a non-letter, and we are parsing an ambiguous mod

				udtCurrentMod.ResidueEnd = intResidueNumber
				blnParsingAmbiguousMod = False

				' The mod mass should be next, in the form [-30.09]
				' Parse out the mod mass
				If intCharIndex < strPrimarySequence.Length - 2 Then
					If strPrimarySequence.Chars(intCharIndex + 1) = "[" Then
						Dim strModMassString = strPrimarySequence.Substring(intCharIndex + 2)
						Dim bracketIndex = strModMassString.IndexOf("]"c)
						If bracketIndex > 0 Then
							' Valid ambiguous mod found; store it
							strModMassString = strModMassString.Substring(0, bracketIndex)
							udtCurrentMod.ModMassString = String.Copy(strModMassString)

							lstAmbiguousMods.Add(udtCurrentMod.ResidueStart, udtCurrentMod)
						End If
					End If
				End If
			End If
		Next

		Return lstAmbiguousMods

	End Function

	Protected Sub HandleException(ByVal strBaseMessage As String, ByVal ex As Exception)
		If String.IsNullOrEmpty(strBaseMessage) Then
			strBaseMessage = "Error"
		End If

		ReportError(strBaseMessage & ": " & ex.Message)
	End Sub

	''' <summary>
	''' Reads the data in strModSummaryFilePath.  Populates mModInfo with the modification names, masses, and affected residues
	''' </summary>
	''' <returns>True if success; false if an error</returns>
	Protected Function LoadModSummary() As Boolean

		Dim objModSummaryReader As clsPHRPModSummaryReader

		Dim strModSummaryFilePath As String
		Dim strModSummaryFilePathPreferred As String

		Dim blnSuccess As Boolean

		Try
			strModSummaryFilePath = clsPHRPReader.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName)
			If String.IsNullOrEmpty(strModSummaryFilePath) Then
				ReportWarning("ModSummaryFile path is empty; unable to continue")
				Return False
			End If

			strModSummaryFilePath = Path.Combine(mInputFolderPath, strModSummaryFilePath)
			strModSummaryFilePathPreferred = clsPHRPReader.AutoSwitchToFHTIfRequired(strModSummaryFilePath, mInputFilePath)
			If strModSummaryFilePath <> strModSummaryFilePathPreferred AndAlso File.Exists(strModSummaryFilePathPreferred) Then
				strModSummaryFilePath = strModSummaryFilePathPreferred
			End If

			If Not File.Exists(strModSummaryFilePath) Then
				ReportWarning("ModSummary file not found: " & strModSummaryFilePath)
				Return False
			End If

			ShowMessage("Reading the PHRP ModSummary file")

			objModSummaryReader = New clsPHRPModSummaryReader(strModSummaryFilePath)
			blnSuccess = objModSummaryReader.Success

			If blnSuccess Then
				mModInfo = objModSummaryReader.ModificationDefs
			End If

		Catch ex As Exception
			HandleException("Exception reading PHRP Mod Summary file", ex)
			Return False
		End Try

		Return blnSuccess

	End Function

	Protected Function LoadSeqInfo() As Boolean

		Dim blnSuccess As Boolean
		Dim objReader As clsPHRPSeqMapReader

		Try

			ShowMessage("Reading the PHRP SeqInfo file")

			' Instantiate the reader
			objReader = New clsPHRPSeqMapReader(mDatasetName, mInputFolderPath, mPeptideHitResultType, mInputFilePath)

			' Read the files
			blnSuccess = objReader.GetProteinMapping(mResultToSeqMap, mSeqToProteinMap, mSeqInfo, mPepToProteinMap)

			If Not blnSuccess Then
				ReportWarning(objReader.ErrorMessage)
			End If

			mResultIDToProteins.Clear()

			If blnSuccess Then
				' Populate mResultIDToProteins

				For Each objItem As KeyValuePair(Of Integer, Integer) In mResultToSeqMap

					'intResultID = objItem.Key
					'intSeqID = objItem.Value

					Dim lstProteinsForSeqID As List(Of clsProteinInfo) = Nothing
					Dim lstProteinsForResultID = New SortedSet(Of String)

					If mSeqToProteinMap.TryGetValue(objItem.Value, lstProteinsForSeqID) Then

						For Each objProtein As clsProteinInfo In lstProteinsForSeqID
							If Not lstProteinsForResultID.Contains(objProtein.ProteinName) Then
								lstProteinsForResultID.Add(objProtein.ProteinName)
							End If
						Next

					End If

					mResultIDToProteins.Add(objItem.Key, lstProteinsForResultID)

				Next

			End If

		Catch ex As Exception
			HandleException("Error loading PHRP Seq Info", ex)
			blnSuccess = False
			If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
		End Try

		Return blnSuccess

	End Function

	''' <summary>
	''' Formats a number so that it begins with a + sign if positive or a - sign if negative
	''' Rounds the number to the specified number of digits, trimming off trailing zeros
	''' Example output: +79.9663 or -17.016
	''' </summary>
	''' <param name="Value"></param>
	''' <param name="DigitsOfPrecision"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function NumToStringPlusMinus(ByVal Value As Double, ByVal DigitsOfPrecision As Integer) As String

		Dim strFormatString As String = "+0;-0"
		If DigitsOfPrecision > 0 Then
			strFormatString = "+0." & New String("0"c, DigitsOfPrecision) & ";-0." & New String("0"c, DigitsOfPrecision)
		End If

		Dim strValue As String = Value.ToString(strFormatString).TrimEnd("0"c)

		If strValue.EndsWith("."c) Then
			' Ends in a decimal point; remove the decimal point
			strValue = strValue.TrimEnd("."c)
		End If

		Return strValue

	End Function

	''' <summary>
	''' Parse the column names in strSplitLine and update the local column header mapping
	''' </summary>
	''' <param name="strSplitLine"></param>
	''' <remarks></remarks>
	Public Sub ParseColumnHeaders(ByRef strSplitLine() As String)
		clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders)
	End Sub

	''' <summary>
	''' Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
	''' </summary>
	''' <param name="strText"></param>
	''' <param name="chDelimiter"></param>
	''' <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
	''' <remarks>Automatically trims whitespace</remarks>
	Public Shared Function ParseKeyValueSetting(ByVal strText As String, ByVal chDelimiter As Char) As KeyValuePair(Of String, String)
		Return ParseKeyValueSetting(strText, chDelimiter, String.Empty)
	End Function

	''' <summary>
	''' Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
	''' </summary>
	''' <param name="strText"></param>
	''' <param name="chDelimiter"></param>
	''' <param name="strCommentChar">If defined, then looks for this character in the value portion of the setting and removes that character plus any text after it</param>
	''' <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
	''' <remarks>Automatically trims whitespace</remarks>
	Public Shared Function ParseKeyValueSetting(ByVal strText As String, ByVal chDelimiter As Char, ByVal strCommentChar As String) As KeyValuePair(Of String, String)
		Dim kvSetting As KeyValuePair(Of String, String)
		Dim strKey As String
		Dim strValue As String
		Dim intCharIndex As Integer

		If Not String.IsNullOrEmpty(strText) Then
			intCharIndex = strText.IndexOf(chDelimiter)
			If intCharIndex > 0 Then
				strKey = strText.Substring(0, intCharIndex).Trim()
				If intCharIndex < strText.Length - 1 Then
					strValue = strText.Substring(intCharIndex + 1).Trim()

					If Not String.IsNullOrEmpty(strCommentChar) Then
						' Look for the comment character
						Dim intCommentCharIndex As Integer
						intCommentCharIndex = strValue.IndexOf(strCommentChar, StringComparison.Ordinal)
						If intCommentCharIndex > 0 Then
							' Trim off the comment
							strValue = strValue.Substring(0, intCommentCharIndex).Trim()
						End If
					End If

				Else
					strValue = String.Empty
				End If

				kvSetting = New KeyValuePair(Of String, String)(strKey, strValue)
				Return kvSetting
			End If
		End If

		Return New KeyValuePair(Of String, String)(String.Empty, String.Empty)

	End Function

	Protected Function ReadKeyValuePairSearchEngineParamFile(ByVal strSearchEngineName As String, ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
		Dim strParamFilePath As String
		Dim blnSuccess As Boolean

		Dim strLineIn As String

		Dim kvSetting As KeyValuePair(Of String, String)

		Try
			If String.IsNullOrWhiteSpace(strSearchEngineName) Then strSearchEngineName = "?? Unknown tool ??"

			strParamFilePath = Path.Combine(mInputFolderPath, strSearchEngineParamFileName)

			If Not File.Exists(strParamFilePath) Then
				ReportError(strSearchEngineName & " param file not found: " & strParamFilePath)
			Else
				Using srInFile As StreamReader = New StreamReader(New FileStream(strParamFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))

					While srInFile.Peek > -1
						strLineIn = srInFile.ReadLine().TrimStart()

						If Not String.IsNullOrWhiteSpace(strLineIn) AndAlso Not strLineIn.StartsWith("#") AndAlso strLineIn.Contains("="c) Then

							' Split the line on the equals sign
							kvSetting = ParseKeyValueSetting(strLineIn, "="c, "#")

							If Not String.IsNullOrEmpty(kvSetting.Key) Then
								objSearchEngineParams.AddUpdateParameter(kvSetting)
							End If
						End If

					End While
				End Using

				blnSuccess = True

			End If
		Catch ex As Exception
			ReportError("Error in ReadKeyValuePairSearchEngineParamFile for " & strSearchEngineName & ": " & ex.Message)
		End Try

		Return blnSuccess

	End Function

	Protected Function ReadSearchEngineVersion(ByVal strFolderPath As String, ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

		Dim strToolVersionInfoFilePath As String
		Dim strLineIn As String

		Dim strSearchEngineVersion As String
		Dim dtSearchDate As DateTime

		Dim kvSetting As KeyValuePair(Of String, String)

		Dim blnValidDate As Boolean
		Dim blnValidVersion As Boolean
		Dim blnSuccess As Boolean = False

		Try
			' Read the Tool_Version_Info file to determine the analysis time and the tool version
			strToolVersionInfoFilePath = Path.Combine(mInputFolderPath, clsPHRPReader.GetToolVersionInfoFilename(ePeptideHitResultType))

			If Not File.Exists(strToolVersionInfoFilePath) Then
				ReportWarning("Tool version info file not found: " & strToolVersionInfoFilePath)
			Else
				strSearchEngineVersion = "Unknown"
				dtSearchDate = New DateTime(1980, 1, 1)

				Using srInFile As StreamReader = New StreamReader(New FileStream(strToolVersionInfoFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))

					While srInFile.Peek > -1
						strLineIn = srInFile.ReadLine().TrimStart()

						' Split the line on a colon
						kvSetting = ParseKeyValueSetting(strLineIn, ":"c)

						Select Case kvSetting.Key.ToLower()
							Case "date"
								blnValidDate = DateTime.TryParse(kvSetting.Value, dtSearchDate)

							Case "toolversioninfo"
								If Not String.IsNullOrEmpty(kvSetting.Value) Then
									strSearchEngineVersion = String.Copy(kvSetting.Value)
									blnValidVersion = True
								Else
									' The next line contains the search engine version
									If srInFile.Peek > -1 Then
										strLineIn = srInFile.ReadLine().TrimStart()
										strSearchEngineVersion = String.Copy(strLineIn)
										blnValidVersion = True
									End If
								End If
							Case Else
								' Ignore the line
						End Select
					End While

				End Using

				If Not blnValidDate Then
					ReportError("Date line not found in the ToolVersionInfo file")
					blnSuccess = False
				ElseIf Not blnValidVersion Then
					ReportError("ToolVersionInfo line not found in the ToolVersionInfo file")
					blnSuccess = False
				Else
					blnSuccess = True
				End If

				objSearchEngineParams.UpdateSearchEngineVersion(strSearchEngineVersion)
				objSearchEngineParams.UpdateSearchDate(dtSearchDate)

			End If

		Catch ex As Exception
			ReportError("Error in ReadSearchEngineVersion: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

	Protected Sub ReportError(ByVal strErrorMessage As String)
		mErrorMessage = strErrorMessage
		mErrorMessages.Add(strErrorMessage)
		RaiseEvent ErrorEvent(strErrorMessage)
	End Sub

	Protected Sub ReportWarning(ByVal strWarningMessage As String)
		mWarningMessages.Add(strWarningMessage)
		RaiseEvent WarningEvent(strWarningMessage)
	End Sub

	Protected Sub ShowMessage(ByVal strMessage As String)
		RaiseEvent MessageEvent(strMessage)
	End Sub

	Protected Sub StoreModInfo(ByRef objPSM As clsPSM, ByVal objSeqInfo As clsSeqInfo)
		
		Dim strMods() As String
		Dim kvModDetails As KeyValuePair(Of String, String)

		Dim strMassCorrectionTag As String
		Dim intResidueLoc As Integer

		Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

		Dim blnMatchFound As Boolean

		Dim blnFavorTerminalMods As Boolean

		Dim lstNTerminalModsAdded = New List(Of String)
		Dim lstCTerminalModsAdded = New List(Of String)
		Dim intPeptideResidueCount = objPSM.PeptideCleanSequence.Length

		objPSM.PeptideMonoisotopicMass = objSeqInfo.MonoisotopicMass

		objPSM.ClearModifiedResidues()

		If objSeqInfo.ModCount > 0 Then
			' Split objSeqInfo.ModDescription on the comma character
			strMods = objSeqInfo.ModDescription.Split(","c)

			If Not strMods Is Nothing AndAlso strMods.Count > 0 Then

				' Parse objPSM.Peptide to look for ambiguous mods, for example -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q

				Dim lstAmbiguousMods = ExtractAmbiguousMods(objPSM.Peptide)

				For intModIndex As Integer = 0 To strMods.Count - 1

					' Split strMods on the colon characters
					kvModDetails = ParseKeyValueSetting(strMods(intModIndex), ":"c)

					If Not String.IsNullOrEmpty(kvModDetails.Key) AndAlso Not String.IsNullOrEmpty(kvModDetails.Value) Then
						strMassCorrectionTag = kvModDetails.Key
						If Integer.TryParse(kvModDetails.Value, intResidueLoc) Then
							' Find the modification definition in mModInfo
							' Note that a given mass correction tag might be present multiple times in mModInfo, since it could be used as both a static peptide mod and a static peptide terminal mod
							' Thus, if intResidueLoc = 1 or intResidueLoc = objPSM.PeptideCleanSequence.Length then we'll first look for a peptide or protein terminal static mod

							If intResidueLoc = 1 Then
								eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
								If lstNTerminalModsAdded.Contains(strMassCorrectionTag) Then
									' We have likely already added this modification as an N-terminal mod, thus, don't favor terminal mods this time
									' An example is an iTraq peptide where there is a K at the N-terminus
									' It gets modified with iTraq twice: once because of the N-terminus and once because of Lysine
									' For example, R.K+144.102063+144.102063TGSY+79.9663GALAEITASK+144.102063.E
									blnFavorTerminalMods = False
								Else
									blnFavorTerminalMods = True
								End If

							ElseIf intResidueLoc = intPeptideResidueCount Then
								eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
								If lstCTerminalModsAdded.Contains(strMassCorrectionTag) Then
									blnFavorTerminalMods = False
								Else
									blnFavorTerminalMods = True
								End If

							Else
								eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
								blnFavorTerminalMods = False
							End If

							Dim objMatchedModDef As clsModificationDefinition = Nothing

							If mModInfo Is Nothing Then
								objMatchedModDef = New clsModificationDefinition()
								objMatchedModDef.MassCorrectionTag = strMassCorrectionTag
								blnMatchFound = True
							Else
								blnMatchFound = UpdatePSMFindMatchingModInfo(strMassCorrectionTag, blnFavorTerminalMods, eResidueTerminusState, objMatchedModDef)
							End If

							If blnMatchFound Then
								Dim lstMatches = (From item In lstAmbiguousMods Where item.Key = intResidueLoc Select item.Value).ToList()

								If lstMatches.Count > 0 Then
									' Ambiguous modification
									objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence.Chars(intResidueLoc - 1), intResidueLoc, eResidueTerminusState, objMatchedModDef, lstMatches.First.ResidueEnd)
								Else
									' Normal, non-ambiguous modified residue
									objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence.Chars(intResidueLoc - 1), intResidueLoc, eResidueTerminusState, objMatchedModDef)
								End If

								If intResidueLoc = 1 Then
									lstNTerminalModsAdded.Add(objMatchedModDef.MassCorrectionTag)
								ElseIf intResidueLoc = intPeptideResidueCount Then
									lstCTerminalModsAdded.Add(objMatchedModDef.MassCorrectionTag)
								End If
							Else
								' Could not find a valid entry in mModInfo
								ReportError("Unrecognized mass correction tag found in the SeqInfo file: " & strMassCorrectionTag)
							End If

						End If
					End If
				Next intModIndex
			End If

		End If

	End Sub

	''' <summary>
	''' Updates the theoretical (computed) monoisotopic mass of objPSM using mResultToSeqMap and mSeqInfo
	''' Also updates the modification info
	''' Also updates SeqID
	''' </summary>
	''' <param name="objPSM"></param>
	''' <returns>True if success, False if objPSM.ResultID is not found in mResultToSeqMap</returns>
	''' <remarks></remarks>
	Protected Function UpdatePSMUsingSeqInfo(ByRef objPSM As clsPSM) As Boolean
		Dim intSeqID As Integer
		Dim objSeqInfo As clsSeqInfo = Nothing

		Dim blnSuccess As Boolean

		blnSuccess = False

		' First determine the modified residues present in this peptide
		If Not mResultToSeqMap Is Nothing AndAlso mResultToSeqMap.Count > 0 Then
			If mResultToSeqMap.TryGetValue(objPSM.ResultID, intSeqID) Then

				objPSM.SeqID = intSeqID

				If mSeqInfo.TryGetValue(intSeqID, objSeqInfo) Then
					StoreModInfo(objPSM, objSeqInfo)
					blnSuccess = True
				End If

				' Lookup the protein details using mSeqToProteinMap
				Dim lstProteinDetails As List(Of clsProteinInfo) = Nothing
				If mSeqToProteinMap.TryGetValue(intSeqID, lstProteinDetails) Then
					For Each oProtein In lstProteinDetails
						objPSM.AddProteinDetail(oProtein)
					Next
				End If

				' Make sure all of the proteins in objPSM.Proteins are defined in objPSM.ProteinDetails
				For Each proteinName In objPSM.Proteins
					Dim blnMatchFound As Boolean = False
					For Each oProtein In objPSM.ProteinDetails
						If String.Equals(oProtein.ProteinName, proteinName, StringComparison.CurrentCultureIgnoreCase) Then
							blnMatchFound = True
							Exit For
						End If
					Next

					If Not blnMatchFound Then
						Dim oProtein = New clsProteinInfo(proteinName, 0, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)
						objPSM.ProteinDetails.Add(oProtein)
					End If
				Next

				' Make sure all of the proteins in objPSM.ProteinDetails are defined in objPSM.Proteins
				For Each oProtein In objPSM.ProteinDetails
					If Not objPSM.Proteins.Contains(oProtein.ProteinName, StringComparer.CurrentCultureIgnoreCase) Then
						objPSM.Proteins.Add(oProtein.ProteinName)
					End If
				Next

				If mPepToProteinMap.Count > 0 Then
					Dim oPepToProteinMapInfo As clsPepToProteinMapInfo = Nothing
					If mPepToProteinMap.TryGetValue(objPSM.PeptideCleanSequence, oPepToProteinMapInfo) Then

						For Each oProtein In objPSM.ProteinDetails

							' Find the matching protein in oPepToProteinMapInfo
							For Each udtProteinMapInfo In oPepToProteinMapInfo.ProteinMapInfo
								If String.Equals(udtProteinMapInfo.Protein, oProtein.ProteinName, StringComparison.CurrentCulture) Then
									oProtein.UpdateLocationInProtein(udtProteinMapInfo.ResidueStart, udtProteinMapInfo.ResidueEnd)
								End If
							Next
						Next

					End If
				End If

			End If
		End If

		If blnSuccess Then
			Dim strPrimarySequence As String = String.Empty
			Dim strPrefix As String = String.Empty
			Dim strSuffix As String = String.Empty

			If clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(objPSM.Peptide, strPrimarySequence, strPrefix, strSuffix) Then
				objPSM.PeptideWithNumericMods = strPrefix & "." & ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues) & "." & strSuffix
			Else
				objPSM.PeptideWithNumericMods = ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues)
			End If

		End If

		Return blnSuccess
	End Function

	Protected Function UpdatePSMFindMatchingModInfo( _
	  ByVal strMassCorrectionTag As String, _
	  ByVal blnFavorTerminalMods As Boolean, _
	  ByVal eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
	  ByRef objMatchedModDef As clsModificationDefinition) As Boolean

		If mModInfo Is Nothing Then Return False

		Dim blnMatchFound As Boolean

		Dim lstMatchedDefs As List(Of clsModificationDefinition)
		lstMatchedDefs = New List(Of clsModificationDefinition)

		For Each objMod In mModInfo
			If strMassCorrectionTag.ToLower() = objMod.MassCorrectionTag.ToLower() Then
				lstMatchedDefs.Add(objMod)
			End If
		Next

		blnMatchFound = False
		If lstMatchedDefs.Count > 0 Then

			Do

				If blnFavorTerminalMods Then
					' Look for an entry in lstMatchedDefs that is a terminal mod
					For Each objMod In lstMatchedDefs
						If objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse _
						   objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod Then

							If eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus AndAlso _
							  (objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
								blnMatchFound = True
								objMatchedModDef = objMod
								Exit For
							End If

							If eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus AndAlso _
							  (objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
								blnMatchFound = True
								objMatchedModDef = objMod
								Exit For
							End If
						End If
					Next
				Else
					' Look for an entry in lstMatchedDefs that is not a terminal mod
					For Each objMod In lstMatchedDefs
						If Not (objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse _
						  objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod) Then
							blnMatchFound = True
							objMatchedModDef = objMod
							Exit For
						End If
					Next
				End If

				If Not blnMatchFound Then
					If blnFavorTerminalMods Then
						blnFavorTerminalMods = False
					Else
						' Still no match found (this shouldn't happen); use the first entry in lstMatchedDefs
						objMatchedModDef = lstMatchedDefs(0)
						blnMatchFound = True
					End If
				End If

			Loop While Not blnMatchFound

		End If

		Return blnMatchFound

	End Function

End Class
