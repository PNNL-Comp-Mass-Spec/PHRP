'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 07/20/2010
'
' This class parses data lines from Inspect _inspect_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserInspect
	Inherits clsPHRPParser

	Public Const DATA_COLUMN_ResultID As String = "ResultID"
	Public Const DATA_COLUMN_Scan As String = "Scan"
	Public Const DATA_COLUMN_Peptide As String = "Peptide"
	Public Const DATA_COLUMN_Protein As String = "Protein"
	Public Const DATA_COLUMN_Charge As String = "Charge"
	Public Const DATA_COLUMN_MQScore As String = "MQScore"
	Public Const DATA_COLUMN_Length As String = "Length"
	Public Const DATA_COLUMN_TotalPRMScore As String = "TotalPRMScore"
	Public Const DATA_COLUMN_MedianPRMScore As String = "MedianPRMScore"
	Public Const DATA_COLUMN_FractionY As String = "FractionY"
	Public Const DATA_COLUMN_FractionB As String = "FractionB"
	Public Const DATA_COLUMN_Intensity As String = "Intensity"
	Public Const DATA_COLUMN_NTT As String = "NTT"
	Public Const DATA_COLUMN_PValue As String = "PValue"
	Public Const DATA_COLUMN_FScore As String = "FScore"
	Public Const DATA_COLUMN_DeltaScore As String = "DeltaScore"
	Public Const DATA_COLUMN_DeltaScoreOther As String = "DeltaScoreOther"
	Public Const DATA_COLUMN_DeltaNormMQScore As String = "DeltaNormMQScore"
	Public Const DATA_COLUMN_DeltaNormTotalPRMScore As String = "DeltaNormTotalPRMScore"
	Public Const DATA_COLUMN_RankTotalPRMScore As String = "RankTotalPRMScore"
	Public Const DATA_COLUMN_RankFScore As String = "RankFScore"
	Public Const DATA_COLUMN_MH As String = "MH"
	Public Const DATA_COLUMN_RecordNumber As String = "RecordNumber"
	Public Const DATA_COLUMN_DBFilePos As String = "DBFilePos"
	Public Const DATA_COLUMN_SpecFilePos As String = "SpecFilePos"
	Public Const DATA_COLUMN_PrecursorMZ As String = "PrecursorMZ"
	Public Const DATA_COLUMN_PrecursorError As String = "PrecursorError"
	Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"

	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		MyBase.New(strDatasetName, strInputFilePath)
	End Sub

	Protected Overrides Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_ResultID)
		AddHeaderColumn(DATA_COLUMN_Scan)
		AddHeaderColumn(DATA_COLUMN_Peptide)
		AddHeaderColumn(DATA_COLUMN_Protein)
		AddHeaderColumn(DATA_COLUMN_Charge)
		AddHeaderColumn(DATA_COLUMN_MQScore)
		AddHeaderColumn(DATA_COLUMN_Length)
		AddHeaderColumn(DATA_COLUMN_TotalPRMScore)
		AddHeaderColumn(DATA_COLUMN_MedianPRMScore)
		AddHeaderColumn(DATA_COLUMN_FractionY)
		AddHeaderColumn(DATA_COLUMN_FractionB)
		AddHeaderColumn(DATA_COLUMN_Intensity)
		AddHeaderColumn(DATA_COLUMN_NTT)
		AddHeaderColumn(DATA_COLUMN_PValue)
		AddHeaderColumn(DATA_COLUMN_FScore)
		AddHeaderColumn(DATA_COLUMN_DeltaScore)
		AddHeaderColumn(DATA_COLUMN_DeltaScoreOther)
		AddHeaderColumn(DATA_COLUMN_DeltaNormMQScore)
		AddHeaderColumn(DATA_COLUMN_DeltaNormTotalPRMScore)
		AddHeaderColumn(DATA_COLUMN_RankTotalPRMScore)
		AddHeaderColumn(DATA_COLUMN_RankFScore)
		AddHeaderColumn(DATA_COLUMN_MH)
		AddHeaderColumn(DATA_COLUMN_RecordNumber)
		AddHeaderColumn(DATA_COLUMN_DBFilePos)
		AddHeaderColumn(DATA_COLUMN_SpecFilePos)
		AddHeaderColumn(DATA_COLUMN_PrecursorMZ)
		AddHeaderColumn(DATA_COLUMN_PrecursorError)
		AddHeaderColumn(DATA_COLUMN_DelM_PPM)

	End Sub

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_fht.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn.txt"
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_SeqToProteinMap.txt"
	End Function

	Public Overrides Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean

		Dim strColumns() As String = strLine.Split(ControlChars.Tab)
		Dim strProtein As String
		Dim dblPrecursorMZ As Double

		Dim blnSuccess As Boolean

		Try

			objPSM.Clear()

			With objPSM
				.ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
				If .ScanNumber = -100 Then
					' Data line is not valid
				Else
					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0)
					.Peptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)
					.UpdateCleavageInfo(mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
					.AddProtein(strProtein)

					dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorError, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					blnSuccess = True
				End If
			End With


			If blnSuccess Then
				' Store the remaining scores

				AddScore(objPSM, strColumns, DATA_COLUMN_MQScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_TotalPRMScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_MedianPRMScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_PValue)
				AddScore(objPSM, strColumns, DATA_COLUMN_FScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_DeltaScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_DeltaScoreOther)
				AddScore(objPSM, strColumns, DATA_COLUMN_DeltaNormMQScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_DeltaNormTotalPRMScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_RankTotalPRMScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_RankFScore)

			End If



		Catch ex As Exception
			MyBase.ReportError("Error parsing line " & intLinesRead & " in the Inspect data file: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

End Class
