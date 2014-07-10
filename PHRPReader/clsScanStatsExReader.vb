'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/20/2012
'
' This class reads MASIC Extended ScanStats data from a tab-delimited _ScanStatsEx.txt file
'
'*********************************************************************************************************

Option Strict On

Imports System.IO

Public Class clsExtendedScanStatsReader

#Region "Constants"
	Public Const DATA_COLUMN_Dataset As String = "Dataset"
	Public Const DATA_COLUMN_ScanNumber As String = "ScanNumber"
	Public Const DATA_COLUMN_IonInjectionTime As String = "Ion Injection Time (ms)"
	Public Const DATA_COLUMN_ScanEvent As String = "Scan Event"
	Public Const DATA_COLUMN_MasterIndex As String = "Master Index"
	Public Const DATA_COLUMN_ElapsedScanTime As String = "Elapsed Scan Time (sec)"
	Public Const DATA_COLUMN_ChargeState As String = "Charge State"
	Public Const DATA_COLUMN_MonoisotopicMZ As String = "Monoisotopic M/Z"
	Public Const DATA_COLUMN_MS2IsolationWidth As String = "MS2 Isolation Width"
	Public Const DATA_COLUMN_FTAnalyzerSettings As String = "FT Analyzer Settings"
	Public Const DATA_COLUMN_FTAnalyzerMessage As String = "FT Analyzer Message"
	Public Const DATA_COLUMN_FTResolution As String = "FT Resolution"
	Public Const DATA_COLUMN_ConversionParameterB As String = "Conversion Parameter B"
	Public Const DATA_COLUMN_ConversionParameterC As String = "Conversion Parameter C"
	Public Const DATA_COLUMN_ConversionParameterD As String = "Conversion Parameter D"
	Public Const DATA_COLUMN_ConversionParameterE As String = "Conversion Parameter E"
	Public Const DATA_COLUMN_CollisionMode As String = "Collision Mode"
	Public Const DATA_COLUMN_ScanFilterText As String = "Scan Filter Text"
	Public Const DATA_COLUMN_SourceVoltage As String = "Source Voltage (kV)"
	Public Const DATA_COLUMN_Source_Current As String = "Source Current (uA)"
#End Region

#Region "Class-wide variables"
	' Column headers
	Protected mColumnHeaders As SortedDictionary(Of String, Integer)
	Protected mErrorMessage As String = String.Empty
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

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_Dataset)
		AddHeaderColumn(DATA_COLUMN_ScanNumber)
		AddHeaderColumn(DATA_COLUMN_IonInjectionTime)
		AddHeaderColumn(DATA_COLUMN_ScanEvent)
		AddHeaderColumn(DATA_COLUMN_MasterIndex)
		AddHeaderColumn(DATA_COLUMN_ElapsedScanTime)
		AddHeaderColumn(DATA_COLUMN_ChargeState)
		AddHeaderColumn(DATA_COLUMN_MonoisotopicMZ)
		AddHeaderColumn(DATA_COLUMN_MS2IsolationWidth)
		AddHeaderColumn(DATA_COLUMN_FTAnalyzerSettings)
		AddHeaderColumn(DATA_COLUMN_FTAnalyzerMessage)
		AddHeaderColumn(DATA_COLUMN_FTResolution)
		AddHeaderColumn(DATA_COLUMN_ConversionParameterB)
		AddHeaderColumn(DATA_COLUMN_ConversionParameterC)
		AddHeaderColumn(DATA_COLUMN_ConversionParameterD)
		AddHeaderColumn(DATA_COLUMN_ConversionParameterE)
		AddHeaderColumn(DATA_COLUMN_CollisionMode)
		AddHeaderColumn(DATA_COLUMN_ScanFilterText)
		AddHeaderColumn(DATA_COLUMN_SourceVoltage)
		AddHeaderColumn(DATA_COLUMN_Source_Current)

	End Sub

	''' <summary>
	''' Open a tab-delimited _ScanStatsEx.txt file and read the data
	''' </summary>
	''' <param name="strInputFilePath">Input file path</param>
	''' <returns>A Dictionary where keys are ScanNumber and values are clsScanStatsInfo objects</returns>
	Public Function ReadExtendedScanStatsData(ByVal strInputFilePath As String) As Dictionary(Of Integer, clsScanStatsExInfo)

		Dim lstScanStats As Dictionary(Of Integer, clsScanStatsExInfo)
		lstScanStats = New Dictionary(Of Integer, clsScanStatsExInfo)

		Dim strLineIn As String
		Dim strSplitLine() As String
		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Dim intLinesRead As Integer
		Dim intScanNumber As Integer

		Try
			DefineColumnHeaders()
			intLinesRead = 0
			mErrorMessage = String.Empty

			Using srInFile As StreamReader = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
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

							intScanNumber = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanNumber, mColumnHeaders, -1)

							If intScanNumber >= 0 AndAlso Not lstScanStats.ContainsKey(intScanNumber) Then

								Dim objScanStatsInfo As clsScanStatsExInfo
								objScanStatsInfo = New clsScanStatsExInfo(intScanNumber)

								objScanStatsInfo.IonInjectionTime = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_IonInjectionTime, mColumnHeaders, 0.0#)
								objScanStatsInfo.ScanEvent = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanEvent, mColumnHeaders, 0)
								objScanStatsInfo.MasterIndex = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_MasterIndex, mColumnHeaders, 0)
								objScanStatsInfo.ElapsedScanTime = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ElapsedScanTime, mColumnHeaders, 0.0#)
								objScanStatsInfo.ChargeState = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ChargeState, mColumnHeaders, 0)
								objScanStatsInfo.MonoisotopicMZ = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_MonoisotopicMZ, mColumnHeaders, 0.0#)
								objScanStatsInfo.MS2IsolationWidth = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_MS2IsolationWidth, mColumnHeaders, 0.0#)
								objScanStatsInfo.FTAnalyzerSettings = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_FTAnalyzerSettings, mColumnHeaders)
								objScanStatsInfo.FTAnalyzerMessage = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_FTAnalyzerMessage, mColumnHeaders)
								objScanStatsInfo.FTResolution = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_FTResolution, mColumnHeaders, 0.0#)
								objScanStatsInfo.ConversionParameterB = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterB, mColumnHeaders, 0.0#)
								objScanStatsInfo.ConversionParameterC = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterC, mColumnHeaders, 0.0#)
								objScanStatsInfo.ConversionParameterD = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterD, mColumnHeaders, 0.0#)
								objScanStatsInfo.ConversionParameterE = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterE, mColumnHeaders, 0.0#)
								objScanStatsInfo.CollisionMode = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_CollisionMode, mColumnHeaders)
								objScanStatsInfo.ScanFilterText = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanFilterText, mColumnHeaders)
								objScanStatsInfo.SourceVoltage = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_SourceVoltage, mColumnHeaders, 0.0#)
								objScanStatsInfo.Source_Current = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_Source_Current, mColumnHeaders, 0.0#)

								lstScanStats.Add(intScanNumber, objScanStatsInfo)
							End If

						End If
					End If

				Loop
			End Using

		Catch ex As Exception
			mErrorMessage = "Error reading the ScanStatsEx data: " & ex.Message
		End Try

		Return lstScanStats

	End Function



End Class
