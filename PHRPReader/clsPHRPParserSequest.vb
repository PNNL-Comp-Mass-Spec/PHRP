'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class parses data lines from Sequest _syn.txt and _fht.txt ifles
'
'*********************************************************************************************************

Option Strict On

Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserSequest
	Inherits clsPHRPParser

	Public Const DATA_COLUMN_HitNum As String = "HitNum"
	Public Const DATA_COLUMN_ScanNum As String = "ScanNum"
	Public Const DATA_COLUMN_ScanCount As String = "ScanCount"
	Public Const DATA_COLUMN_ChargeState As String = "ChargeState"
	Public Const DATA_COLUMN_MH As String = "MH"
	Public Const DATA_COLUMN_XCorr As String = "XCorr"
	Public Const DATA_COLUMN_DelCn As String = "DelCn"
	Public Const DATA_COLUMN_Sp As String = "Sp"
	Public Const DATA_COLUMN_Reference As String = "Reference"
	Public Const DATA_COLUMN_MultiProtein As String = "MultiProtein"
	Public Const DATA_COLUMN_Peptide As String = "Peptide"
	Public Const DATA_COLUMN_DelCn2 As String = "DelCn2"
	Public Const DATA_COLUMN_RankSp As String = "RankSp"
	Public Const DATA_COLUMN_RankXc As String = "RankXc"
	Public Const DATA_COLUMN_DelM As String = "DelM"
	Public Const DATA_COLUMN_XcRatio As String = "XcRatio"
	Public Const DATA_COLUMN_PassFilt As String = "PassFilt"
	Public Const DATA_COLUMN_MScore As String = "MScore"
	Public Const DATA_COLUMN_Ions_Observed As String = "Ions_Observed"
	Public Const DATA_COLUMN_Ions_Expected As String = "Ions_Expected"
	Public Const DATA_COLUMN_NumTrypticEnds As String = "NumTrypticEnds"
	Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		MyBase.New(strDatasetName, strInputFilePath)
		mInitialized = True
	End Sub

	Protected Overrides Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_HitNum)
		AddHeaderColumn(DATA_COLUMN_ScanNum)
		AddHeaderColumn(DATA_COLUMN_ScanCount)
		AddHeaderColumn(DATA_COLUMN_ChargeState)
		AddHeaderColumn(DATA_COLUMN_MH)
		AddHeaderColumn(DATA_COLUMN_XCorr)
		AddHeaderColumn(DATA_COLUMN_DelCn)
		AddHeaderColumn(DATA_COLUMN_Sp)
		AddHeaderColumn(DATA_COLUMN_Reference)
		AddHeaderColumn(DATA_COLUMN_MultiProtein)
		AddHeaderColumn(DATA_COLUMN_Peptide)
		AddHeaderColumn(DATA_COLUMN_DelCn2)
		AddHeaderColumn(DATA_COLUMN_RankSp)
		AddHeaderColumn(DATA_COLUMN_RankXc)
		AddHeaderColumn(DATA_COLUMN_DelM)
		AddHeaderColumn(DATA_COLUMN_XcRatio)

		mColumnHeaders.Add(DATA_COLUMN_PassFilt, -1)
		mColumnHeaders.Add(DATA_COLUMN_MScore, -1)

		AddHeaderColumn(DATA_COLUMN_Ions_Observed)
		AddHeaderColumn(DATA_COLUMN_Ions_Expected)

		AddHeaderColumn(DATA_COLUMN_NumTrypticEnds)
		AddHeaderColumn(DATA_COLUMN_DelM_PPM)

	End Sub

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_fht.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn.txt"
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn_SeqToProteinMap.txt"
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
		Dim strProtein As String

		Dim dblPrecursorMH As Double
		Dim dblPrecursorMZ As Double
		Dim dblMassErrorPPM As Double
		Dim dblMassErrorDa As Double

		Dim blnSuccess As Boolean

		Try

			If objPSM Is Nothing Then
				objPSM = New clsPSM
			Else
				objPSM.Clear()
			End If

			With objPSM
				.ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_ScanNum, mColumnHeaders, -100)
				If .ScanNumber = -100 Then
					' Data line is not valid
				Else
					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_HitNum, mColumnHeaders, 0)
					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_ChargeState, mColumnHeaders, 0), Short)

					strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Reference, mColumnHeaders)
					.AddProtein(strProtein)

					' Note that the MH value listed in Sequest files is not the precursor MH but is instead the theoretical MH of the peptide
					dblPrecursorMH = LookupColumnValue(strColumns, DATA_COLUMN_MH, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMH, .Charge, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					If Double.TryParse(.MassErrorPPM, dblMassErrorPPM) Then
						dblPrecursorMZ = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMH, 1, .Charge)
						dblMassErrorDa = clsPeptideMassCalculator.PPMToMass(dblMassErrorPPM, dblPrecursorMZ)
						' Could store this in .MassErrorDa = dblMassErrorDa.ToString("0.00000")
					End If

					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				' Store the remaining scores

				AddScore(objPSM, strColumns, DATA_COLUMN_XCorr)
				AddScore(objPSM, strColumns, DATA_COLUMN_DelCn)
				AddScore(objPSM, strColumns, DATA_COLUMN_Sp)
				AddScore(objPSM, strColumns, DATA_COLUMN_DelCn2)
				AddScore(objPSM, strColumns, DATA_COLUMN_RankSp)
				AddScore(objPSM, strColumns, DATA_COLUMN_RankXc)
				AddScore(objPSM, strColumns, DATA_COLUMN_XcRatio)
				AddScore(objPSM, strColumns, DATA_COLUMN_Ions_Observed)
				AddScore(objPSM, strColumns, DATA_COLUMN_Ions_Expected)

			End If


		Catch ex As Exception
			MyBase.ReportError("Error parsing line " & intLinesRead & " in the Sequest data file: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

End Class
