'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class parses data lines from Xtandem _xt.txt files
'
'*********************************************************************************************************

Option Strict On

Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserXTandem
	Inherits clsPHRPParser

	Public Const DATA_COLUMN_Result_ID As String = "Result_ID"
	Public Const DATA_COLUMN_Group_ID As String = "Group_ID"
	Public Const DATA_COLUMN_Scan As String = "Scan"
	Public Const DATA_COLUMN_Charge As String = "Charge"
	Public Const DATA_COLUMN_Peptide_MH As String = "Peptide_MH"
	Public Const DATA_COLUMN_Peptide_Hyperscore As String = "Peptide_Hyperscore"
	Public Const DATA_COLUMN_Peptide_Expectation_Value_LogE As String = "Peptide_Expectation_Value_Log(e)"
	Public Const DATA_COLUMN_Multiple_Protein_Count As String = "Multiple_Protein_Count"
	Public Const DATA_COLUMN_Peptide_Sequence As String = "Peptide_Sequence"
	Public Const DATA_COLUMN_DeltaCn2 As String = "DeltaCn2"
	Public Const DATA_COLUMN_y_score As String = "y_score"
	Public Const DATA_COLUMN_y_ions As String = "y_ions"
	Public Const DATA_COLUMN_b_score As String = "b_score"
	Public Const DATA_COLUMN_b_ions As String = "b_ions"
	Public Const DATA_COLUMN_Delta_Mass As String = "Delta_Mass"
	Public Const DATA_COLUMN_Peptide_Intensity_LogI As String = "Peptide_Intensity_Log(I)"
	Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"

	Protected Const XT_SEQ_PROT_MAP_COLUMN_Unique_Seq_ID As String = "Unique_Seq_ID"
	Protected Const XT_SEQ_PROT_MAP_COLUMN_Cleavage_State As String = "Cleavage_State"
	Protected Const XT_SEQ_PROT_MAP_COLUMN_Terminus_State As String = "Terminus_State"
	Protected Const XT_SEQ_PROT_MAP_COLUMN_Protein_Name As String = "Protein_Name"
	Protected Const XT_SEQ_PROT_MAP_COLUMN_Protein_EValue As String = "Protein_Expectation_Value_Log(e)"
	Protected Const XT_SEQ_PROT_MAP_COLUMN_Protein_Intensity As String = "Protein_Intensity_Log(I)"

	Protected mProteinByResultID As System.Collections.Generic.SortedList(Of Integer, String)

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		MyBase.New(strDatasetName, strInputFilePath)

		mProteinByResultID = New System.Collections.Generic.SortedList(Of Integer, String)
		LoadXTandemResultProteins()

		mInitialized = True
	End Sub

	Protected Overrides Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_Result_ID)
		AddHeaderColumn(DATA_COLUMN_Group_ID)
		AddHeaderColumn(DATA_COLUMN_Scan)
		AddHeaderColumn(DATA_COLUMN_Charge)
		AddHeaderColumn(DATA_COLUMN_Peptide_MH)
		AddHeaderColumn(DATA_COLUMN_Peptide_Hyperscore)
		AddHeaderColumn(DATA_COLUMN_Peptide_Expectation_Value_LogE)
		AddHeaderColumn(DATA_COLUMN_Multiple_Protein_Count)
		AddHeaderColumn(DATA_COLUMN_Peptide_Sequence)
		AddHeaderColumn(DATA_COLUMN_DeltaCn2)
		AddHeaderColumn(DATA_COLUMN_y_score)
		AddHeaderColumn(DATA_COLUMN_y_ions)
		AddHeaderColumn(DATA_COLUMN_b_score)
		AddHeaderColumn(DATA_COLUMN_b_ions)
		AddHeaderColumn(DATA_COLUMN_Delta_Mass)
		AddHeaderColumn(DATA_COLUMN_Peptide_Intensity_LogI)
		AddHeaderColumn(DATA_COLUMN_DelM_PPM)

	End Sub

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		' X!Tandem does not have a first-hits file; just the _xt.txt file
		Return String.Empty
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt.txt"
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_SeqToProteinMap.txt"
	End Function

	Protected Function LoadXTandemResultProteins() As Boolean

		' Tracks the ResultIDs that map to each SeqID
		Dim objSeqToResultMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))

		Dim blnSuccess As Boolean

		Try
			' Initialize the tracking lists
			objSeqToResultMap = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))

			mProteinByResultID.Clear()

			blnSuccess = LoadXTandemResultToSeqMapping(objSeqToResultMap)

			If blnSuccess Then
				blnSuccess = LoadXTandemSeqToProteinMapping(objSeqToResultMap, mProteinByResultID)
			End If

		Catch ex As Exception
			HandleException("Error loading X!Tandem protein results", ex)
			blnSuccess = False
			If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
		End Try

		Return blnSuccess

	End Function

	Protected Function LoadXTandemResultToSeqMapping(ByRef objSeqToResultMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))) As Boolean

		Dim objResultIDList As System.Collections.Generic.List(Of Integer) = Nothing

		Dim strFilePath As String
		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim intLinesRead As Integer

		Dim intResultID As Integer
		Dim intSeqID As Integer

		Try

			' Read the data from the result to sequence map file
			strFilePath = System.IO.Path.Combine(mInputFolderPath, mDatasetName & XT_RESULT_TO_SEQ_MAP_SUFFIX)
			If Not System.IO.File.Exists(strFilePath) Then
				Return False
			End If

			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))

				intLinesRead = 0

				Do While srInFile.Peek >= 0
					strLineIn = srInFile.ReadLine
					intLinesRead += 1

					If Not String.IsNullOrEmpty(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If strSplitLine.Length >= 2 Then

							' Parse out the numbers from the first two columns 
							' (the first line of the file is the header line, and it will get skipped)
							If Integer.TryParse(strSplitLine(0), intResultID) Then
								If Integer.TryParse(strSplitLine(1), intSeqID) Then

									If objSeqToResultMap.TryGetValue(intSeqID, objResultIDList) Then
										objResultIDList.Add(intResultID)
									Else
										objResultIDList = New System.Collections.Generic.List(Of Integer)
										objResultIDList.Add(intResultID)
										objSeqToResultMap.Add(intSeqID, objResultIDList)
									End If
								End If
							End If

						End If
					End If
				Loop

			End Using

		Catch ex As Exception
			HandleException("Error reading X!Tandem result to seq map file", ex)
			Return False
		End Try

		Return True

	End Function

	Protected Function LoadXTandemSeqToProteinMapping(ByRef objSeqToResultMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer)), _
	  ByRef objProteinByResultID As System.Collections.Generic.SortedList(Of Integer, String)) As Boolean

		Dim objResultIDList As System.Collections.Generic.List(Of Integer) = Nothing

		Dim objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)

		Dim strFilePath As String
		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim strProtein As String
		Dim intLinesRead As Integer

		Dim intResultID As Integer
		Dim intSeqID As Integer
		Dim intSeqIDPrevious As Integer

		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Try

			' Initialize the column mapping
			' Using a case-insensitive comparer
			objColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

			' Define the default column mapping
			objColumnHeaders.Add(XT_SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, 0)
			objColumnHeaders.Add(XT_SEQ_PROT_MAP_COLUMN_Cleavage_State, 1)
			objColumnHeaders.Add(XT_SEQ_PROT_MAP_COLUMN_Terminus_State, 2)
			objColumnHeaders.Add(XT_SEQ_PROT_MAP_COLUMN_Protein_Name, 3)
			objColumnHeaders.Add(XT_SEQ_PROT_MAP_COLUMN_Protein_EValue, 4)
			objColumnHeaders.Add(XT_SEQ_PROT_MAP_COLUMN_Protein_Intensity, 5)

			' Read the data from the sequence to protein map file
			strFilePath = System.IO.Path.Combine(mInputFolderPath, mDatasetName & XT_SEQ_TO_PROTEIN_MAP_SUFFIX)
			If Not System.IO.File.Exists(strFilePath) Then
				Return False
			End If

			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))


				intLinesRead = 0
				intSeqIDPrevious = 0

				Do While srInFile.Peek >= 0
					strLineIn = srInFile.ReadLine
					intLinesRead += 1
					blnSkipLine = False

					If Not String.IsNullOrEmpty(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If Not blnHeaderLineParsed Then
							If strSplitLine(0).ToLower() = XT_SEQ_PROT_MAP_COLUMN_Unique_Seq_ID.ToLower() Then
								' Parse the header line to confirm the column ordering
								clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders)
								blnSkipLine = True
							End If

							blnHeaderLineParsed = True
						End If

						If Not blnSkipLine AndAlso strSplitLine.Length >= 3 Then

							If Integer.TryParse(strSplitLine(0), intSeqID) Then
								If intSeqID <> intSeqIDPrevious Then
									strProtein = clsPHRPReader.LookupColumnValue(strSplitLine, XT_SEQ_PROT_MAP_COLUMN_Protein_Name, objColumnHeaders, String.Empty)

									If Not String.IsNullOrEmpty(strProtein) Then
										' Find the ResultIDs in objResultToSeqMap() that have sequence ID intSeqID
										If objSeqToResultMap.TryGetValue(intSeqID, objResultIDList) Then

											For Each intResultID In objResultIDList
												If Not objProteinByResultID.ContainsKey(intResultID) Then
													objProteinByResultID.Add(intResultID, strProtein)
												End If
											Next

										End If

									End If

								End If
							End If

						End If

					End If
				Loop

			End Using

		Catch ex As Exception
			HandleException("Error reading X!Tandem seq to protein map file", ex)
			Return False
		End Try

		Return True

	End Function

	''' <summary>
	''' Parse the data line read from a PHRP results file
	''' </summary>
	''' <param name="strLine">Data line</param>
	''' <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
	''' <param name="objPSM">clsPSM object (output)</param>
	''' <returns>True if success, false if an error</returns>
	Public Overrides Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean

		Dim strColumns() As String = strLine.Split(ControlChars.Tab)
		Dim strPeptide As String
		Dim strProtein As String = String.Empty

		Dim dblPrecursorMH As Double

		Dim blnSuccess As Boolean

		Try

			If objPSM Is Nothing Then
				objPSM = New clsPSM
			Else
				objPSM.Clear()
			End If

			With objPSM
				.ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
				If .ScanNumber = -100 Then
					' Data line is not valid
				Else
					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_Result_ID, mColumnHeaders, 0)
					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide_Sequence, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					' Lookup the protein name using mProteinByResultID
					If mProteinByResultID.TryGetValue(.ResultID, strProtein) Then
						.AddProtein(strProtein)
					End If

					dblPrecursorMH = LookupColumnValue(strColumns, DATA_COLUMN_Peptide_MH, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMH, 1, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_Delta_Mass, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				' Store the remaining scores
				AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Hyperscore)
				AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Expectation_Value_LogE)
				AddScore(objPSM, strColumns, DATA_COLUMN_DeltaCn2)
				AddScore(objPSM, strColumns, DATA_COLUMN_y_score)
				AddScore(objPSM, strColumns, DATA_COLUMN_y_ions)
				AddScore(objPSM, strColumns, DATA_COLUMN_b_score)
				AddScore(objPSM, strColumns, DATA_COLUMN_b_ions)
				AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Intensity_LogI)
			End If

		Catch ex As Exception
			HandleException("Error parsing line " & intLinesRead & " in the X!Tandem data file", ex)
		End Try

		Return blnSuccess

	End Function

End Class
