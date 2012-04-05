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

Public MustInherit Class clsPHRPParser

#Region "Module variables"

	Protected mDatasetName As String
	Protected mInputFilePath As String
	Protected mInputFolderPath As String
	Protected mInitialized As Boolean

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

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset Name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
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

#Region "Functions overridden by derived classes"
	Protected MustOverride Sub DefineColumnHeaders()
	Public MustOverride Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean
#End Region

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub AddScore(ByRef objPSM As clsPSM, ByRef strColumns() As String, ByVal strScoreColumnName As String)
		Const NOT_FOUND As String = "==SCORE_NOT_FOUND=="

		Dim strValue As String
		strValue = PHRPReader.clsPHRPReader.LookupColumnValue(strColumns, strScoreColumnName, mColumnHeaders, NOT_FOUND)

		If strValue <> NOT_FOUND Then
			objPSM.SetScore(strScoreColumnName, strValue)
		End If

	End Sub

	Protected Sub HandleException(ByVal strBaseMessage As String, ByVal ex As System.Exception)
		If String.IsNullOrEmpty(strBaseMessage) Then
			strBaseMessage = "Error"
		End If

		ReportError(strBaseMessage & ": " & ex.Message)
	End Sub

	''' <summary>
	''' Parse the column names in strSplitLine and update the local column header mapping
	''' </summary>
	''' <param name="strSplitLine"></param>
	''' <remarks></remarks>
	Public Sub ParseColumnHeaders(ByRef strSplitLine() As String)
		clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders)
	End Sub

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
