'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/03/2012
'
' This class reads MSGF scores from a tab-delimited _msgf.txt file
'
'*********************************************************************************************************

Option Strict On
Imports System.IO

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
    Private ReadOnly mColumnHeaders As SortedDictionary(Of String, Integer)
    Private mErrorMessage As String = String.Empty
#End Region

	''' <summary>
	''' Error message
	''' </summary>
	Public ReadOnly Property ErrorMessage As String
		Get
			If String.IsNullOrEmpty(mErrorMessage) Then
				Return String.Empty
			Else
				Return mErrorMessage
			End If
		End Get
	End Property

	''' <summary>
	''' Constructor
	''' </summary>
	''' <remarks></remarks>
	Public Sub New()
		mColumnHeaders = New SortedDictionary(Of String, Integer)
	End Sub

    Private Sub AddHeaderColumn(strColumnName As String)
        mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
    End Sub

    Private Sub DefineColumnHeaders()

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

    ''' <summary>
    ''' Open a tab-delimited MSGF results file and read the data
    ''' </summary>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <returns>A Dictionary where keys are ResultID and values are MSGF_SpecProb values (stored as strings)</returns>
    Public Function ReadMSGFData(strInputFilePath As String) As Dictionary(Of Integer, String)

        Dim lstMSGFData As Dictionary(Of Integer, String)
        lstMSGFData = New Dictionary(Of Integer, String)

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

            Using srInFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                Do While Not srInFile.EndOfStream
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
