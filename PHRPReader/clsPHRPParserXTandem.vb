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

	' This List tracks the Protein Names for each ResultID
	Protected mResultIDToProteins As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		MyBase.New(strDatasetName, strInputFilePath)

		mResultIDToProteins = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))
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

		Dim objResultToSeqMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))
		Dim objSeqToProteinMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))

		Dim blnSuccess As Boolean
		Dim objReader As clsPHRPSeqMapReader

		Try

			mResultIDToProteins.Clear()

			' Initialize the tracking lists
			objResultToSeqMap = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))
			objSeqToProteinMap = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))

			' Instantiate the reader
			objReader = New clsPHRPSeqMapReader(mDatasetName, mInputFolderPath, ePeptideHitResultType.XTandem)

			' Read the files
			blnSuccess = objReader.GetProteinMapping(objResultToSeqMap, objSeqToProteinMap)

			If blnSuccess Then
				' Populate mResultIDToProteins

				Dim intResultID As Integer
				Dim intSeqID As Integer

				Dim lstProteinsForSeqID As System.Collections.Generic.List(Of String) = Nothing
				Dim lstProteinsForResultID As System.Collections.Generic.List(Of String) = Nothing

				For Each objItem As System.Collections.Generic.KeyValuePair(Of Integer, System.Collections.Generic.List(Of Integer)) In objResultToSeqMap

					intResultID = objItem.Key
					lstProteinsForResultID = New System.Collections.Generic.List(Of String)

					For Each intSeqID In objItem.Value()
						If objSeqToProteinMap.TryGetValue(intSeqID, lstProteinsForSeqID) Then

							For Each strProtein As String In lstProteinsForSeqID
								If Not lstProteinsForResultID.Contains(strProtein) Then
									lstProteinsForResultID.Add(strProtein)
								End If
							Next

						End If
					Next

					mResultIDToProteins.Add(intResultID, lstProteinsForResultID)

				Next

			End If

		Catch ex As Exception
			HandleException("Error loading X!Tandem protein results", ex)
			blnSuccess = False
			If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
		End Try

		Return blnSuccess

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
		Dim lstProteinsForResultID As System.Collections.Generic.List(Of String) = Nothing

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

					' Lookup the protein name(s) using mResultIDToProteins
					If mResultIDToProteins.TryGetValue(.ResultID, lstProteinsForResultID) Then
						For Each strProtein As String In lstProteinsForResultID
							.AddProtein(strProtein)
						Next
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
