Option Strict On

Public Class clsMSGFResultsReader

#Region "Constants"
	Public Const DATA_COLUMN_ResultID As String = "Result_ID"
	Public Const DATA_COLUMN_Scan As String = "Scan"
	Public Const DATA_COLUMN_Charge As String = "Charge"
	Public Const DATA_COLUMN_Protein As String = "Protein"
	Public Const DATA_COLUMN_Peptide As String = "Peptide"
	Public Const DATA_COLUMN_SpecProb As String = "SpecProb"
	Public Const DATA_COLUMN_Notes As String = "Notes"
#End Region

#Region "Class-wide variables"
	' Column headers
	Protected mColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)
	Protected mErrorMessage As String = String.Empty
#End Region

	Public ReadOnly Property ErrorMessage As String
		Get
			If String.IsNullOrEmpty(mErrorMessage) Then
				Return String.Empty
			Else
				Return mErrorMessage
			End If
		End Get
	End Property

	Public Sub New()
		mColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)
	End Sub

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_ResultID)
		AddHeaderColumn(DATA_COLUMN_Scan)
		AddHeaderColumn(DATA_COLUMN_Charge)
		AddHeaderColumn(DATA_COLUMN_Protein)
		AddHeaderColumn(DATA_COLUMN_Peptide)
		AddHeaderColumn(DATA_COLUMN_SpecProb)
		AddHeaderColumn(DATA_COLUMN_Notes)

	End Sub

	Public Function ReadMSGFData(ByVal strInputFilePath As String) As System.Collections.Generic.Dictionary(Of Integer, String)

		Dim lstMSGFData As System.Collections.Generic.Dictionary(Of Integer, String)
		lstMSGFData = New System.Collections.Generic.Dictionary(Of Integer, String)

		Dim strLineIn As String
		Dim strSplitLine() As String
		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Dim intLinesRead As Integer
		Dim intResultID As Integer
		Dim strMSGFSpecProb As String

		Try
			DefineColumnHeaders()
			intLinesRead = 0
			mErrorMessage = String.Empty

			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))
				Do While srInFile.Peek() > -1
					strLineIn = srInFile.ReadLine()
					intLinesRead += 1
					blnSkipLine = False

					If Not String.IsNullOrWhiteSpace(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If Not blnHeaderLineParsed Then
							If Not clsPHRPReader.IsNumber(strSplitLine(0)) Then
								' Parse the header line to confirm the column ordering
								clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders)
								blnSkipLine = True
							End If

							blnHeaderLineParsed = True
						End If

						If Not blnSkipLine AndAlso strSplitLine.Length >= 4 Then

							intResultID = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ResultID, mColumnHeaders, -1)

							If intResultID >= 0 Then
								strMSGFSpecProb = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_SpecProb, mColumnHeaders)

								If Not String.IsNullOrEmpty(strMSGFSpecProb) AndAlso Not lstMSGFData.ContainsKey(intResultID) Then
									lstMSGFData.Add(intResultID, strMSGFSpecProb)
								End If
							End If

						End If
					End If

				Loop
			End Using

		Catch ex As Exception
			mErrorMessage = "Error reading the MSGF data: " & ex.Message
		End Try

		Return lstMSGFData
	End Function

End Class
