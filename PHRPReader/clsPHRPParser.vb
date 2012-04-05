'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class is the base class for classes used to parse PHRP data lines
'
'*********************************************************************************************************

Option Strict On

Public MustInherit Class clsPHRPParser

#Region "Module variables"

	Protected mDatasetName As String
	Protected mInputFilePath As String
	Protected mInputFolderPath As String

	' Column headers in the synopsis file and first hits file
	Protected mColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)

	Protected mErrorMessage As String = String.Empty

	Protected mCleavageStateCalculator As clsPeptideCleavageStateCalculator
	Protected mPeptideMassCalculator As clsPeptideMassCalculator

#End Region

#Region "Events"
	Public Event MessageEvent(ByVal strMessage As String)
	Public Event ErrorEvent(ByVal strErrorMessage As String)
	Public Event WarningEvent(ByVal strWarningMessage As String)
#End Region

	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)

		mDatasetName = strDatasetName
		mInputFilePath = strInputFilePath

		Dim fiFileInfo As System.IO.FileInfo = New System.IO.FileInfo(strInputFilePath)
		mInputFolderPath = fiFileInfo.DirectoryName

		mErrorMessage = String.Empty

		' Initialize the column mapping object
		' Using a case-insensitive comparer
		mColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

		mCleavageStateCalculator = New clsPeptideCleavageStateCalculator()
		mPeptideMassCalculator = New clsPeptideMassCalculator()

		' The following will be overridden by a derived form of this class
		DefineColumnHeaders()

	End Sub

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub AddScore(ByRef objPSM As clsPSM, ByRef strColumns() As String, ByVal strScoreColumnName As String)

		objPSM.SetScore(strScoreColumnName, PHRPReader.clsPHRPReader.LookupColumnValue(strColumns, strScoreColumnName, mColumnHeaders))

	End Sub

	Protected MustOverride Sub DefineColumnHeaders()

	Protected Sub HandleException(ByVal strBaseMessage As String, ByVal ex As System.Exception)
		If String.IsNullOrEmpty(strBaseMessage) Then
			strBaseMessage = "Error"
		End If

		ReportError(strBaseMessage & ": " & ex.Message)
	End Sub

	Public Sub ParseColumnHeaders(ByRef strSplitLine() As String)
		clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders)
	End Sub

	Public MustOverride Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean

	Protected Sub ReportError(ByVal strErrorMessage As String)
		mErrorMessage = strErrorMessage
		RaiseEvent ErrorEvent(strErrorMessage)
	End Sub

	Protected Sub ReportWarning(ByVal strWarningMessage As String)
		RaiseEvent WarningEvent(strWarningMessage)
	End Sub

	Protected Sub ShowMessage(ByVal strMessage As String)
		RaiseEvent MessageEvent(strMessage)
	End Sub


End Class
