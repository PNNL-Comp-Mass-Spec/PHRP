
'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/20/2012
'
' This class reads MASIC ScanStats data from a tab-delimited _ScanStats.txt file
'
'*********************************************************************************************************

Option Strict On

Public Class clsScanStatsReader

#Region "Constants"
	Public Const DATA_COLUMN_Dataset As String = "Dataset"
	Public Const DATA_COLUMN_ScanNumber As String = "ScanNumber"
	Public Const DATA_COLUMN_ScanTime As String = "ScanTime"
	Public Const DATA_COLUMN_ScanType As String = "ScanType"
	Public Const DATA_COLUMN_TotalIonIntensity As String = "TotalIonIntensity"
	Public Const DATA_COLUMN_BasePeakIntensity As String = "BasePeakIntensity"
	Public Const DATA_COLUMN_BasePeakMZ As String = "BasePeakMZ"
	Public Const DATA_COLUMN_BasePeakSignalToNoiseRatio As String = "BasePeakSignalToNoiseRatio"
	Public Const DATA_COLUMN_IonCount As String = "IonCount"
	Public Const DATA_COLUMN_IonCountRaw As String = "IonCountRaw"
	Public Const DATA_COLUMN_ScanTypeName As String = "ScanTypeName"
#End Region

#Region "Class-wide variables"
	' Column headers
	Protected mColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)
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
		mColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)
	End Sub

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_Dataset)
		AddHeaderColumn(DATA_COLUMN_ScanNumber)
		AddHeaderColumn(DATA_COLUMN_ScanTime)
		AddHeaderColumn(DATA_COLUMN_ScanType)
		AddHeaderColumn(DATA_COLUMN_TotalIonIntensity)
		AddHeaderColumn(DATA_COLUMN_BasePeakIntensity)
		AddHeaderColumn(DATA_COLUMN_BasePeakMZ)
		AddHeaderColumn(DATA_COLUMN_BasePeakSignalToNoiseRatio)
		AddHeaderColumn(DATA_COLUMN_IonCount)
		AddHeaderColumn(DATA_COLUMN_IonCountRaw)
		AddHeaderColumn(DATA_COLUMN_ScanTypeName)

	End Sub

	''' <summary>
	''' Open a tab-delimited Scan_Stats file and read the data
	''' </summary>
	''' <param name="strInputFilePath">Input file path</param>
	''' <returns>A Dictionary where keys are ScanNumber and values are clsScanStatsInfo objects</returns>
	Public Function ReadScanStatsData(ByVal strInputFilePath As String) As System.Collections.Generic.Dictionary(Of Integer, clsScanStatsInfo)

		Dim lstScanStats As System.Collections.Generic.Dictionary(Of Integer, clsScanStatsInfo)
		lstScanStats = New System.Collections.Generic.Dictionary(Of Integer, clsScanStatsInfo)

		Dim strLineIn As String
		Dim strSplitLine() As String
		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Dim intLinesRead As Integer
		Dim intScanNumber As Integer
		Dim sngScanTimeMinutes As Single
		Dim intScanType As Integer

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

							intScanNumber = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanNumber, mColumnHeaders, -1)
							sngScanTimeMinutes = CSng(clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanTime, mColumnHeaders, 0.0#))
							intScanType = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanType, mColumnHeaders, 0)

							If intScanNumber >= 0 AndAlso Not lstScanStats.ContainsKey(intScanNumber) Then

								Dim objScanStatsInfo As clsScanStatsInfo
								objScanStatsInfo = New clsScanStatsInfo(intScanNumber, sngScanTimeMinutes, intScanType)

								objScanStatsInfo.TotalIonIntensity = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_TotalIonIntensity, mColumnHeaders, 0.0#)
								objScanStatsInfo.BasePeakIntensity = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_BasePeakIntensity, mColumnHeaders, 0.0#)
								objScanStatsInfo.BasePeakMZ = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_BasePeakMZ, mColumnHeaders, 0.0#)
								objScanStatsInfo.BasePeakSignalToNoiseRatio = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_BasePeakSignalToNoiseRatio, mColumnHeaders, 0.0#)
								objScanStatsInfo.IonCount = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_IonCount, mColumnHeaders, 0)
								objScanStatsInfo.IonCountRaw = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_IonCountRaw, mColumnHeaders, 0)
								objScanStatsInfo.ScanTypeName = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanTypeName, mColumnHeaders)

								lstScanStats.Add(intScanNumber, objScanStatsInfo)
							End If

						End If
					End If

				Loop
			End Using

		Catch ex As Exception
			mErrorMessage = "Error reading the ScanStats data: " & ex.Message
		End Try

		Return lstScanStats

	End Function

End Class
