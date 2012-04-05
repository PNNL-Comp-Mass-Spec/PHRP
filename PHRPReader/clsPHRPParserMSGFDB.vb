'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 07/20/2010
'
' This class parses data lines from MSGFDB msgfdb_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserMSGFDB
	Inherits clsPHRPParser

	Public Const DATA_COLUMN_ResultID As String = "ResultID"
	Public Const DATA_COLUMN_Scan As String = "Scan"
	Public Const DATA_COLUMN_FragMethod As String = "FragMethod"
	Public Const DATA_COLUMN_SpecIndex As String = "SpecIndex"
	Public Const DATA_COLUMN_Charge As String = "Charge"
	Public Const DATA_COLUMN_PrecursorMZ As String = "PrecursorMZ"
	Public Const DATA_COLUMN_DelM As String = "DelM"
	Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"
	Public Const DATA_COLUMN_MH As String = "MH"
	Public Const DATA_COLUMN_Peptide As String = "Peptide"
	Public Const DATA_COLUMN_Protein As String = "Protein"
	Public Const DATA_COLUMN_NTT As String = "NTT"
	Public Const DATA_COLUMN_DeNovoScore As String = "DeNovoScore"
	Public Const DATA_COLUMN_MSGFScore As String = "MSGFScore"
	Public Const DATA_COLUMN_MSGFDB_SpecProb As String = "MSGFDB_SpecProb"
	Public Const DATA_COLUMN_Rank_MSGFDB_SpecProb As String = "Rank_MSGFDB_SpecProb"
	Public Const DATA_COLUMN_PValue As String = "PValue"
	Public Const DATA_COLUMN_FDR As String = "FDR"
	Public Const DATA_COLUMN_PepFDR As String = "PepFDR"

	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		MyBase.New(strDatasetName, strInputFilePath)
	End Sub


	Protected Overrides Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping	
		AddHeaderColumn(DATA_COLUMN_ResultID)
		AddHeaderColumn(DATA_COLUMN_Scan)
		AddHeaderColumn(DATA_COLUMN_FragMethod)
		AddHeaderColumn(DATA_COLUMN_SpecIndex)
		AddHeaderColumn(DATA_COLUMN_Charge)
		AddHeaderColumn(DATA_COLUMN_PrecursorMZ)
		AddHeaderColumn(DATA_COLUMN_DelM)
		AddHeaderColumn(DATA_COLUMN_DelM_PPM)
		AddHeaderColumn(DATA_COLUMN_MH)
		AddHeaderColumn(DATA_COLUMN_Peptide)
		AddHeaderColumn(DATA_COLUMN_Protein)
		AddHeaderColumn(DATA_COLUMN_NTT)
		AddHeaderColumn(DATA_COLUMN_DeNovoScore)
		AddHeaderColumn(DATA_COLUMN_MSGFScore)
		AddHeaderColumn(DATA_COLUMN_MSGFDB_SpecProb)
		AddHeaderColumn(DATA_COLUMN_Rank_MSGFDB_SpecProb)
		AddHeaderColumn(DATA_COLUMN_PValue)
		AddHeaderColumn(DATA_COLUMN_FDR)
		AddHeaderColumn(DATA_COLUMN_PepFDR)

	End Sub

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_fht.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn.txt"
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_SeqToProteinMap.txt"
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

					.CollisionMode = LookupColumnValue(strColumns, DATA_COLUMN_FragMethod, mColumnHeaders, "n/a")

					dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					.MSGFSpecProb = LookupColumnValue(strColumns, DATA_COLUMN_MSGFScore, mColumnHeaders)

					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				' Store the remaining scores
				AddScore(objPSM, strColumns, DATA_COLUMN_DeNovoScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_MSGFDB_SpecProb)
				AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb)
				AddScore(objPSM, strColumns, DATA_COLUMN_PValue)
				AddScore(objPSM, strColumns, DATA_COLUMN_FDR)
				AddScore(objPSM, strColumns, DATA_COLUMN_PepFDR)
			End If

		Catch ex As Exception
			MyBase.ReportError("Error parsing line " & intLinesRead & " in the MSGFDB data file: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

End Class
