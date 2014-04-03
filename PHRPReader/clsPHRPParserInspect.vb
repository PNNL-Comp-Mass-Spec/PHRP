'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class parses data lines from Inspect _inspect_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserInspect
	Inherits clsPHRPParser

#Region "Constants"
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

	Public Const FILENAME_SUFFIX_SYN As String = "_inspect_syn.txt"
	Public Const FILENAME_SUFFIX_FHT As String = "_inspect_fht.txt"

	Protected Const INS_SEARCH_ENGINE_NAME As String = "Inspect"
#End Region

	''' <summary>
	''' Constructor; assumes blnLoadModsAndSeqInfo=True
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		Me.New(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo:=True)
	End Sub

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <param name="blnLoadModsAndSeqInfo">If True, then load the ModSummary file and SeqInfo files</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal blnLoadModsAndSeqInfo As Boolean)
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.Inspect, blnLoadModsAndSeqInfo)
		mInitialized = True
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

	''' <summary>
	''' Determines the precursor mass tolerance
	''' </summary>
	''' <param name="objSearchEngineParams"></param>
	''' <param name="dblTolerancePPM">Precursor mass tolerance, in ppm</param>
	''' <returns>Precursor tolerance, in Da</returns>
	''' <remarks></remarks>
	Protected Function DeterminePrecursorMassTolerance(ByRef objSearchEngineParams As clsSearchEngineParameters, <Out()> ByRef dblTolerancePPM As Double) As Double
		Dim strTolerance As String = String.Empty

		Dim dblToleranceDa As Double = 0
		Dim dblTolerance As Double = 0

		dblTolerancePPM = 0

		If objSearchEngineParams.Parameters.TryGetValue("ParentPPM", strTolerance) Then
			' Parent mass tolerance, in ppm
			Double.TryParse(strTolerance, dblTolerancePPM)
			' Convert from PPM to dalton (assuming a mass of 2000 m/z)
			dblTolerance = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000)
		End If

		If objSearchEngineParams.Parameters.TryGetValue("PMTolerance", strTolerance) Then
			' Parent mass tolerance, in Da
			Double.TryParse(strTolerance, dblToleranceDa)

			' Convert from dalton to PPM (assuming a mass of 2000 m/z)
			dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000)
		End If

		dblTolerance = Math.Max(dblTolerance, dblToleranceDa)

		Return dblTolerance

	End Function

	Public Overloads Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_FHT
	End Function

	Public Overloads Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_ModSummary.txt"
	End Function

	Public Overloads Shared Function GetPHRPPepToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_PepToProtMapMTS.txt"
	End Function

	Public Overloads Shared Function GetPHRPProteinModsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_ProteinMods.txt"
	End Function

	Public Overloads Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_SYN
	End Function

	Public Overloads Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_ResultToSeqMap.txt"
	End Function

	Public Overloads Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_SeqInfo.txt"
	End Function

	Public Overloads Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_inspect_syn_SeqToProteinMap.txt"
	End Function

	Public Overloads Shared Function GetSearchEngineName() As String
		Return INS_SEARCH_ENGINE_NAME
	End Function

	''' <summary>
	''' Parses the specified Inspect parameter file
	''' </summary>
	''' <param name="strSearchEngineParamFileName"></param>
	''' <param name="objSearchEngineParams"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

		Dim blnSuccess As Boolean

		objSearchEngineParams = New clsSearchEngineParameters(INS_SEARCH_ENGINE_NAME, mModInfo)

		blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

		ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

		Return blnSuccess

	End Function

	Protected Function ReadSearchEngineParamFile(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

		Dim strSettingValue As String = String.Empty
		Dim blnSuccess As Boolean

		Try
			blnSuccess = ReadKeyValuePairSearchEngineParamFile(INS_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, ePeptideHitResultType.Inspect, objSearchEngineParams)

			If blnSuccess Then

				' Determine the enzyme name
				If objSearchEngineParams.Parameters.TryGetValue("protease", strSettingValue) Then
					Select Case strSettingValue.ToLower()
						Case "trypsin"
							objSearchEngineParams.Enzyme = "trypsin"
						Case "none"
							objSearchEngineParams.Enzyme = "no_enzyme"
						Case "chymotrypsin"
							objSearchEngineParams.Enzyme = "chymotrypsin"
						Case Else
							If Not String.IsNullOrEmpty(strSettingValue) Then
								objSearchEngineParams.Enzyme = strSettingValue
							End If
					End Select
				End If
			
				' Determine the precursor mass tolerance (will store 0 if a problem or not found)
				Dim dblTolerancePPM As Double
				objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, dblTolerancePPM)
				objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM
			End If

		Catch ex As Exception
			ReportError("Error in ReadSearchEngineParamFile: " & ex.Message)
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
		Dim strProtein As String
		Dim dblPrecursorMZ As Double

		Dim blnSuccess As Boolean

		Try

			If objPSM Is Nothing Then
				objPSM = New clsPSM
			Else
				objPSM.Clear()
			End If

			With objPSM
				.DataLineText = strLine
				.ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
				If .ScanNumber = -100 Then
					' Data line is not valid
				Else
					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0)
					.ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_RankTotalPRMScore, mColumnHeaders, 0)

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

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
				UpdatePSMUsingSeqInfo(objPSM)

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
