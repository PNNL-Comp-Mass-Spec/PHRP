﻿'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class reads a tab-delimited text file (created by the Peptide File Extractor or by PHRP)
' and returns the data for each peptide hit search result
'
' It also integrates MSGF results with the peptide hit search results
' 
' The class must be derived by a sub-class customized for the specific analysis tool (Sequest, X!Tandem, Inspect, etc.)
'
'*********************************************************************************************************

Option Strict On

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

	Protected Const MOD_SUMMARY_COLUMN_Modification_Symbol As String = "Modification_Symbol"
	Protected Const MOD_SUMMARY_COLUMN_Modification_Mass As String = "Modification_Mass"
	Protected Const MOD_SUMMARY_COLUMN_Target_Residues As String = "Target_Residues"
	Protected Const MOD_SUMMARY_COLUMN_Modification_Type As String = "Modification_Type"
	Protected Const MOD_SUMMARY_COLUMN_Mass_Correction_Tag As String = "Mass_Correction_Tag"
	Protected Const MOD_SUMMARY_COLUMN_Occurence_Count As String = "Occurence_Count"

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

	Protected mSourceFile As System.IO.StreamReader
	Protected mPHRPParser As clsPHRPParser

	' This list contains mod symbols as the key and the corresponding mod mass (stored as a string 
	'  to retain the same number of Sig Figs as the _ModSummary.txt file)
	' This dictionary object will use case-sensitive searching
	Protected mDynamicMods As System.Collections.Generic.SortedDictionary(Of String, String)

	' This dictionary contains amino acid names as the key and the corresponding mod mass (stored as a string 
	'  to retain the same number of Sig Figs as the _ModSummary.txt file)
	' This dictionary object will use case-sensitive searching
	Protected mStaticMods As System.Collections.Generic.SortedDictionary(Of String, String)

	' This dictionary tracks the MSGFSpecProb values for each entry in the source file
	' The key is Result_ID and the string is MSGFSpecProb (stored as string to preserve formatting)
	Protected mMSGFCachedResults As System.Collections.Generic.Dictionary(Of Integer, String)

	Protected mPSMCurrent As clsPSM

	Protected mLinesRead As Integer
	Protected mHeaderLineParsed As Boolean
	Protected mCachedLineAvailable As Boolean
	Protected mCachedLine As String

	Protected mErrorMessage As String = String.Empty
	Protected mLocalErrorCode As ePHRPReaderErrorCodes

#End Region

#Region "Events"
	Public Event MessageEvent(ByVal strMessage As String)
	Public Event ErrorEvent(ByVal strErrorMessage As String)
	Public Event WarningEvent(ByVal strWarningMessage As String)
#End Region

#Region "Properties"

	Public ReadOnly Property ErrorMessage() As String
		Get
			Return mErrorMessage
		End Get
	End Property

#End Region

	Public Sub New(ByVal strInputFilePath As String)
		mErrorMessage = String.Empty
		mLocalErrorCode = ePHRPReaderErrorCodes.NoError

		mMSGFCachedResults = New System.Collections.Generic.Dictionary(Of Integer, String)

		Dim blnSuccess As Boolean
		blnSuccess = InitializeReader(strInputFilePath)

	End Sub

	Protected Function InitializeReader(ByVal strInputFilePath As String) As Boolean

		Dim blnSuccess As Boolean = False
		Dim eResultType As ePeptideHitResultType

		Dim strSearchToolParamFilePath As String = String.Empty
		Dim strModSummaryFilePath As String = String.Empty

		Dim fiFileInfo As System.IO.FileInfo

		mErrorMessage = String.Empty
		SetLocalErrorCode(ePHRPReaderErrorCodes.NoError)

		Try
			If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
				ReportError("Input file name is empty")
				SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
			Else
				' Confirm that the source file exists
				' Make sure strInputFilePath points to a valid file
				fiFileInfo = New System.IO.FileInfo(strInputFilePath)

				mInputFolderPath = fiFileInfo.DirectoryName
				mInputFilePath = strInputFilePath

				If Not fiFileInfo.Exists Then
					ReportError("Input file not found: " & strInputFilePath)
					SetLocalErrorCode(ePHRPReaderErrorCodes.InvalidInputFilePath)
				Else

					' Note that the following populates mDatasetName
					blnSuccess = ValidateInputFiles(strInputFilePath, eResultType, strSearchToolParamFilePath, strModSummaryFilePath)

					If Not blnSuccess Then
						SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound, True)
						Return False
					End If

					mInputFolderPath = fiFileInfo.DirectoryName

					' Open the input file for reading
					blnSuccess = InitializeReader(strInputFilePath, eResultType)

				End If
			End If


		Catch ex As Exception
			HandleException("Error in ProcessFile", ex)
		End Try

		Return blnSuccess
	End Function

	Protected Function InitializeReader(ByVal strInputFilePath As String, ByVal eResultType As ePeptideHitResultType) As Boolean

		Dim strModSummaryFilePath As String

		Dim blnSuccess As Boolean

		' This dictionary contains mod symbols as the key and the corresponding mod mass (stored as a string 
		'  to retain the same number of Sig Figs as the _ModSummary.txt file)
		' This dictionary object will use case-sensitive searching
		Dim objDynamicMods As New System.Collections.Generic.SortedDictionary(Of String, String)

		' This dictionary contains amino acid names as the key and the corresponding mod mass (stored as a string 
		'  to retain the same number of Sig Figs as the _ModSummary.txt file)
		' This dictionary object will use case-sensitive searching
		Dim objStaticMods As New System.Collections.Generic.SortedDictionary(Of String, String)


		' Read the PHRP Mod Summary File
		strModSummaryFilePath = System.IO.Path.Combine(strInputFilePath, GetModSummaryFileName(eResultType, mDatasetName))
		blnSuccess = ReadModSummaryFile(strModSummaryFilePath, objDynamicMods, objStaticMods)

		If blnSuccess Then

			' Open the peptide-hit result file (from PHRP) for reading
			Select Case eResultType
				Case ePeptideHitResultType.Sequest

					' Convert Sequest results to input format required for MSGF
					mPHRPParser = New clsPHRPParserSequest(mDatasetName, mInputFilePath)

				Case ePeptideHitResultType.XTandem

					' Convert X!Tandem results to input format required for MSGF
					mPHRPParser = New clsPHRPParserXTandem(mDatasetName, mInputFilePath)


					' Make sure a few more files are present for X!Tandem so that we can extract the protein names and include these in the MSGF files
					Dim strFileToCheck As String
					strFileToCheck = System.IO.Path.Combine(mInputFolderPath, mDatasetName & clsPHRPReader.XT_RESULT_TO_SEQ_MAP_SUFFIX)
					If Not ValidateRequiredFileExists("X!Tandem Result to Seq Map file", strFileToCheck, False) Then
						ShowMessage("Warning: X!Tandem Result to Seq Map file not found (" & System.IO.Path.GetFileName(strFileToCheck) & "); PepXML file will not contain protein names")
					Else
						strFileToCheck = System.IO.Path.Combine(mInputFolderPath, mDatasetName & clsPHRPReader.XT_SEQ_TO_PROTEIN_MAP_SUFFIX)
						If Not ValidateRequiredFileExists("X!Tandem Seq to Protein Map file", strFileToCheck, False) Then
							ShowMessage("Warning: X!Tandem Seq to Protein Map file not found (" & System.IO.Path.GetFileName(strFileToCheck) & "); ; PepXML file will not contain protein names")
						End If
					End If


				Case ePeptideHitResultType.Inspect

					' Convert Inspect results to input format required for MSGF
					mPHRPParser = New clsPHRPParserInspect(mDatasetName, mInputFilePath)

				Case ePeptideHitResultType.MSGFDB

					' Convert MSGFDB results to input format required for MSGF
					mPHRPParser = New clsPHRPParserMSGFDB(mDatasetName, mInputFilePath)

				Case Else
					'Should never get here; invalid result type specified
					ReportError("Invalid PeptideHit ResultType specified: " & eResultType)
					blnSuccess = False
			End Select

			If blnSuccess Then

				blnSuccess = OpenDataFile()

			End If

		End If

		Return blnSuccess
	End Function

	Protected Function AutoDetermineDatasetName(ByVal strFilePath As String, ByVal eResultType As ePeptideHitResultType) As String

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

		Return strDatasetName

	End Function

	Protected Function AutoDetermineResultType(ByVal strFilePath As String) As ePeptideHitResultType

		Dim strFilePathLcase As String
		Dim srInFile As System.IO.StreamReader

		Dim strHeaderLine As String

		Dim eResultType As ePeptideHitResultType
		eResultType = ePeptideHitResultType.Unknown

		strFilePathLcase = strFilePath.ToLower

		Try
			If strFilePathLcase.ToLower.EndsWith("_xt.txt") Then
				eResultType = ePeptideHitResultType.XTandem
			Else
				If strFilePathLcase.EndsWith("_msgfdb_syn.txt") OrElse strFilePathLcase.EndsWith("_msgfdb_fht.txt") Then
					eResultType = ePeptideHitResultType.MSGFDB

				ElseIf strFilePathLcase.EndsWith("_inspect_syn.txt") OrElse strFilePathLcase.EndsWith("_inspect_fht.txt") Then
					eResultType = ePeptideHitResultType.Inspect

				ElseIf strFilePathLcase.EndsWith("_syn.txt") OrElse strFilePathLcase.EndsWith("_fht.txt") Then
					' Open the file and read the header line to determine if this is a Sequest file, Inspect file, MSGFDB, or something else

					srInFile = New System.IO.StreamReader(New System.IO.FileStream(strFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

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

					srInFile.Close()
				End If
			End If

		Catch ex As Exception
			HandleException("Error in AutoDetermineResultType", ex)
		End Try

		Return eResultType
	End Function

	Public Shared Function GetModSummaryFileName(ByVal eResultType As ePeptideHitResultType, ByVal strDatasetName As String) As String

		Dim strModSummaryName As String = String.Empty

		Select Case eResultType
			Case ePeptideHitResultType.Sequest
				' Sequest
				strModSummaryName = strDatasetName & "_syn_ModSummary.txt"

			Case ePeptideHitResultType.XTandem
				' X!Tandem
				strModSummaryName = strDatasetName & "_xt_ModSummary.txt"

			Case ePeptideHitResultType.Inspect
				' Inspect
				strModSummaryName = strDatasetName & "_inspect_syn_ModSummary.txt"

			Case ePeptideHitResultType.MSGFDB
				' MSGFDB
				strModSummaryName = strDatasetName & "_msgfdb_syn_ModSummary.txt"

		End Select

		Return strModSummaryName

	End Function

	Protected Sub HandleException(ByVal strBaseMessage As String, ByVal ex As System.Exception)
		If String.IsNullOrEmpty(strBaseMessage) Then
			strBaseMessage = "Error"
		End If

		ReportError(strBaseMessage & ": " & ex.Message)

	End Sub

	''' <summary>
	''' Look for dynamic mod symbols in the peptide sequence; replace with the corresponding mod masses
	''' </summary>
	''' <returns>True if success, false if an error</returns>
	''' <remarks></remarks>
	Protected Function AddDynamicAndStaticMods(ByVal strPeptide As String, ByRef strPeptideWithMods As String) As Boolean

		Static sbNewPeptide As New System.Text.StringBuilder

		Dim intIndex As Integer
		Dim intIndexStart As Integer
		Dim intIndexEnd As Integer

		Try
			If mDynamicMods.Count = 0 AndAlso mStaticMods.Count = 0 Then
				' No mods are defined; simply update strPeptideWithMods to be strPeptide
				strPeptideWithMods = strPeptide
				Return True
			End If

			strPeptideWithMods = String.Empty
			sbNewPeptide.Length = 0

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
			Do While intIndex < strPeptide.Length
				If intIndex < intIndexStart OrElse intIndex > intIndexEnd Then
					' We're before or after the primary peptide sequence; simply append the character
					sbNewPeptide.Append(strPeptide.Chars(intIndex))
				Else
					If Char.IsLetter(strPeptide.Chars(intIndex)) Then
						' Character is a letter; append it
						sbNewPeptide.Append(strPeptide.Chars(intIndex))

						' See if it is present in mStaticMods (this is a case-sensitive search)
						AddModIfPresent(mStaticMods, strPeptide.Chars(intIndex), sbNewPeptide)

						If intIndex = intIndexStart AndAlso mStaticMods.Count > 0 Then
							' We're at the N-terminus of the peptide
							' Possibly add a static N-terminal peptide mod (for example, iTRAQ8, which is 304.2022 DA)
							AddModIfPresent(mStaticMods, clsPHRPReader.N_TERMINAL_PEPTIDE_SYMBOL_DMS, sbNewPeptide)

							If strPeptide.StartsWith(clsPHRPReader.PROTEIN_TERMINUS_SYMBOL_PHRP) Then
								' We're at the N-terminus of the protein
								' Possibly add a static N-terminal protein mod
								AddModIfPresent(mStaticMods, clsPHRPReader.N_TERMINAL_PROTEIN_SYMBOL_DMS, sbNewPeptide)
							End If
						End If
					Else
						' Not a letter; see if it is present in mDynamicMods
						AddModIfPresent(mDynamicMods, strPeptide.Chars(intIndex), sbNewPeptide)
					End If

					If intIndex = intIndexEnd AndAlso mStaticMods.Count > 0 Then
						' Possibly add a static C-terminal peptide mod
						AddModIfPresent(mStaticMods, clsPHRPReader.C_TERMINAL_PEPTIDE_SYMBOL_DMS, sbNewPeptide)

						If strPeptide.EndsWith(clsPHRPReader.PROTEIN_TERMINUS_SYMBOL_PHRP) Then
							' We're at the C-terminus of the protein
							' Possibly add a static C-terminal protein mod
							AddModIfPresent(mStaticMods, clsPHRPReader.C_TERMINAL_PROTEIN_SYMBOL_DMS, sbNewPeptide)
						End If

					End If

				End If
				intIndex += 1
			Loop

			strPeptideWithMods = sbNewPeptide.ToString

		Catch ex As Exception
			HandleException("Error adding dynamic and static mods to peptide " & strPeptide, ex)
			Return False
		End Try

		Return True

	End Function

	Protected Sub AddModIfPresent(ByRef objMods As System.Collections.Generic.SortedDictionary(Of String, String), _
	   ByVal chResidue As Char, _
	   ByRef sbNewPeptide As System.Text.StringBuilder)

		Dim strModMass As String = String.Empty

		If objMods.TryGetValue(chResidue, strModMass) Then
			' Static mod applies to this residue; append the mod (add a plus sign if it doesn't start with a minus sign)
			If strModMass.StartsWith("-") Then
				sbNewPeptide.Append(strModMass)
			Else
				sbNewPeptide.Append("+" & strModMass)
			End If
		End If

	End Sub

	''' <summary>
	''' Opens the data file for reading
	''' </summary>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Overridable Function OpenDataFile() As Boolean

		Dim blnSuccess As Boolean = False

		Try
			If String.IsNullOrEmpty(mDatasetName) Then
				ReportError("Dataset name is undefined; unable to continue")
				Return False
			End If

			' Initialize some tracking variables
			mMSGFCachedResults.Clear()

			mPSMCurrent = New clsPSM()

			mLinesRead = 0
			mHeaderLineParsed = False
			mCachedLineAvailable = False
			mCachedLine = String.Empty

			If Not String.IsNullOrEmpty(mInputFilePath) AndAlso System.IO.File.Exists(mInputFilePath) Then

				' Cache the MSGF values (if present)
				blnSuccess = ReadAndCacheMSGFData()

				' Open the data file for reading
				mSourceFile = New System.IO.StreamReader(New System.IO.FileStream(mInputFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))
				blnSuccess = True
			End If

		Catch ex As Exception
			HandleException("Error opening the PHRP data file", ex)
		End Try

		Return blnSuccess

	End Function

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
	  ByRef objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)) As String

		Return LookupColumnValue(strColumns, strColumnName, objColumnHeaders, String.Empty)
	End Function

	''' <summary>
	''' Returns the string stored in the given named column (using objColumnHeaders to dereference column name with column index)
	''' </summary>
	''' <returns>The text in the specified column; strValueIfMissing if the specific column name is not recognized</returns>
	''' <remarks></remarks>
	Public Shared Function LookupColumnValue(ByRef strColumns() As String, _
	  ByVal strColumnName As String, _
	  ByRef objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer), _
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
	  ByRef objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer), _
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
	  ByRef objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer), _
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
	Public Shared Sub ParseColumnHeaders(ByVal strColumns() As String, _
	 ByRef objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer))

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

		If mCachedLineAvailable Then
			strLineIn = mCachedLine
			mCachedLineAvailable = False
			blnSuccess = True
		ElseIf mSourceFile.Peek > -1 Then
			strLineIn = mSourceFile.ReadLine()
			mLinesRead += 1
			blnSuccess = True
		End If

		If blnSuccess Then
			blnSkipLine = False

			If Not String.IsNullOrEmpty(strLineIn) Then
				strSplitLine = strLineIn.Split(ControlChars.Tab)

				If Not mHeaderLineParsed Then
					If Not clsPHRPReader.IsNumber(strSplitLine(0)) Then
						' Parse the header line to confirm the column ordering
						mPHRPParser.ParseColumnHeaders(strSplitLine)
						blnSkipLine = True
					End If

					mHeaderLineParsed = True
				End If

				If Not blnSkipLine Then

					blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mLinesRead, mPSMCurrent)

					If blnSuccess Then
						blnMatchFound = True

						' Markup the peptide with the dynamic and static mods
						Dim strPeptideWithMods As String = String.Empty
						blnSuccess = AddDynamicAndStaticMods(mPSMCurrent.Peptide.Trim, strPeptideWithMods)
						If blnSuccess Then
							mPSMCurrent.PeptideWithNumericMods = strPeptideWithMods
						End If

						Dim strMSGFSpecProb As String = String.Empty
						If mMSGFCachedResults.TryGetValue(mPSMCurrent.ResultID, strMSGFSpecProb) Then
							mPSMCurrent.MSGFSpecProb = strMSGFSpecProb
						End If

						' Read the next line and check whether it's the same hit, but a different protein
						Dim blnReadNext As Boolean = True
						Do While blnReadNext AndAlso mSourceFile.Peek > -1
							strLineIn = mSourceFile.ReadLine()
							mLinesRead += 1

							If Not String.IsNullOrEmpty(strLineIn) Then

								Dim objNewPSM As New clsPSM
								blnSuccess = mPHRPParser.ParsePHRPDataLine(strLineIn, mLinesRead, objNewPSM)

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
			mMSGFCachedResults.Clear()

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
	  ByRef objDynamicMods As System.Collections.Generic.SortedDictionary(Of String, String), _
	  ByRef objStaticMods As System.Collections.Generic.SortedDictionary(Of String, String)) As Boolean

		Dim srModSummaryFile As System.IO.StreamReader
		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)

		Dim intLinesRead As Integer
		Dim intIndex As Integer

		Dim strModSymbol As String
		Dim strModMass As String
		Dim strTargetResidues As String
		Dim strModType As String

		Dim blnSkipLine As Boolean
		Dim blnHeaderLineParsed As Boolean

		Try
			If String.IsNullOrEmpty(strModSummaryFilePath) Then
				ReportError("ModSummaryFile path is empty; unable to continue")
				Return False
			End If

			ShowMessage("Reading the PHRP ModSummary file")

			' Initialize the column mapping
			' Using a case-insensitive comparer
			objColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

			' Define the default column mapping
			objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Symbol, 0)
			objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Mass, 1)
			objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Target_Residues, 2)
			objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Type, 3)
			objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Mass_Correction_Tag, 4)
			objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Occurence_Count, 5)

			' Clear objDynamicMods and objStaticMods (should have been instantiated by the calling function)
			objDynamicMods.Clear()
			objStaticMods.Clear()


			' Read the data from the ModSummary.txt file
			' The first line is typically a header line:
			' Modification_Symbol	Modification_Mass	Target_Residues	Modification_Type	Mass_Correction_Tag	Occurence_Count

			srModSummaryFile = New System.IO.StreamReader(New System.IO.FileStream(strModSummaryFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))

			blnHeaderLineParsed = False
			intLinesRead = 0

			Do While srModSummaryFile.Peek >= 0
				strLineIn = srModSummaryFile.ReadLine
				intLinesRead += 1
				blnSkipLine = False

				If Not String.IsNullOrEmpty(strLineIn) Then
					strSplitLine = strLineIn.Split(ControlChars.Tab)

					If Not blnHeaderLineParsed Then
						If strSplitLine(0).ToLower() = MOD_SUMMARY_COLUMN_Modification_Symbol.ToLower Then
							' Parse the header line to confirm the column ordering
							clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders)
							blnSkipLine = True
						End If

						blnHeaderLineParsed = True
					End If

					If Not blnSkipLine AndAlso strSplitLine.Length >= 4 Then
						strModSymbol = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Symbol, objColumnHeaders)
						strModMass = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Mass, objColumnHeaders)
						strTargetResidues = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Target_Residues, objColumnHeaders)
						strModType = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Type, objColumnHeaders)

						Select Case strModType.ToUpper()
							Case "S", "T", "P"
								' Static residue mod, peptide terminus static mod, or protein terminus static mod
								' Note that < and > mean peptide N and C terminus (N_TERMINAL_PEPTIDE_SYMBOL_DMS and C_TERMINAL_PEPTIDE_SYMBOL_DMS)
								' Note that [ and ] mean protein N and C terminus (N_TERMINAL_PROTEIN_SYMBOL_DMS and C_TERMINAL_PROTEIN_SYMBOL_DMS)

								' This mod could apply to multiple residues, so need to process each character in strTargetResidues
								For intIndex = 0 To strTargetResidues.Length - 1
									Try
										If objStaticMods.ContainsKey(strTargetResidues.Chars(intIndex)) Then
											' Residue is already present in objStaticMods; this is unexpected
											' We'll log a warning, but continue
											ShowMessage("Warning: Residue '" & strTargetResidues.Chars(intIndex) & "' has more than one static mod defined; this is not allowed (duplicate has ModMass=" & strModMass & ")")
										Else
											objStaticMods.Add(strTargetResidues.Chars(intIndex), strModMass)
										End If

									Catch ex As Exception
										HandleException("Exception adding static mod for " & strTargetResidues.Chars(intIndex) & " with ModMass=" & strModMass, ex)
									End Try
								Next intIndex

							Case Else
								' Dynamic residue mod (Includes mod type "D")
								' Note that < and > mean peptide N and C terminus (N_TERMINAL_PEPTIDE_SYMBOL_DMS and C_TERMINAL_PEPTIDE_SYMBOL_DMS)

								Try
									If objDynamicMods.ContainsKey(strModSymbol) Then
										' Mod symbol already present in objDynamicMods; this is unexpected
										' We'll log a warning, but continue
										ShowMessage("Warning: Dynamic mod symbol '" & strModSymbol & "' is already defined; it cannot have more than one associated mod mass (duplicate has ModMass=" & strModMass & ")")
									Else
										objDynamicMods.Add(strModSymbol, strModMass)
									End If

								Catch ex As Exception
									HandleException("Exception adding dynamic mod for " & strModSymbol & " with ModMass=" & strModMass, ex)
								End Try

						End Select
					End If
				End If

			Loop

			srModSummaryFile.Close()

		Catch ex As Exception
			HandleException("Exception reading PHRP Mod Summary file", ex)
			Return False
		End Try

		Return True

	End Function


	Protected Sub ReportError(ByVal strErrorMessage As String)
		mErrorMessage = strErrorMessage
		RaiseEvent ErrorEvent(strErrorMessage)
	End Sub

	Protected Sub ReportWarning(ByVal strWarningMessage As String)
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
		RaiseEvent MessageEvent(strMessage)
	End Sub


	Private Function ValidateInputFiles(ByVal strInputFilePath As String, _
	  ByRef eResultType As ePeptideHitResultType, _
	  ByRef strSearchToolParamFilePath As String, _
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

		strModSummaryFilePath = GetModSummaryFileName(eResultType, mDatasetName)
		strModSummaryFilePath = System.IO.Path.Combine(fiFileInfo.DirectoryName, strModSummaryFilePath)
		If Not ValidateRequiredFileExists("ModSummary file", strModSummaryFilePath) Then
			SetLocalErrorCode(ePHRPReaderErrorCodes.RequiredInputFileNotFound)
			Return False
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