﻿'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 11/28/2012
'
' This class parses data lines from MSAlign msalign_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserMSAlign
	Inherits clsPHRPParser

#Region "Constants"
	Public Const DATA_COLUMN_ResultID As String = "ResultID"
	Public Const DATA_COLUMN_Scan As String = "Scan"
	Public Const DATA_COLUMN_Prsm_ID As String = "Prsm_ID"
	Public Const DATA_COLUMN_Spectrum_ID As String = "Spectrum_ID"
	Public Const DATA_COLUMN_Charge As String = "Charge"
	Public Const DATA_COLUMN_PrecursorMZ As String = "PrecursorMZ"
	Public Const DATA_COLUMN_DelM As String = "DelM"
	Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"
	Public Const DATA_COLUMN_MH As String = "MH"
	Public Const DATA_COLUMN_Peptide As String = "Peptide"
	Public Const DATA_COLUMN_Protein As String = "Protein"
	Public Const DATA_COLUMN_Protein_Mass As String = "Protein_Mass"
	Public Const DATA_COLUMN_Unexpected_Mod_Count As String = "Unexpected_Mod_Count"
	Public Const DATA_COLUMN_Peak_Count As String = "Peak_Count"
	Public Const DATA_COLUMN_Matched_Peak_Count As String = "Matched_Peak_Count"
	Public Const DATA_COLUMN_Matched_Fragment_Ion_Count As String = "Matched_Fragment_Ion_Count"
	Public Const DATA_COLUMN_PValue As String = "PValue"
	Public Const DATA_COLUMN_Rank_PValue As String = "Rank_PValue"
	Public Const DATA_COLUMN_EValue As String = "EValue"
	Public Const DATA_COLUMN_FDR As String = "FDR"

	Public Const FILENAME_SUFFIX_SYN As String = "_msalign_syn.txt"
	Public Const FILENAME_SUFFIX_FHT As String = "_msalign_fht.txt"

	Protected Const MSAlign_SEARCH_ENGINE_NAME As String = "MSAlign"
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
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSAlign, blnLoadModsAndSeqInfo)
		mInitialized = True
	End Sub

	Protected Overrides Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping	
		AddHeaderColumn(DATA_COLUMN_ResultID)
		AddHeaderColumn(DATA_COLUMN_Scan)
		AddHeaderColumn(DATA_COLUMN_Prsm_ID)
		AddHeaderColumn(DATA_COLUMN_Spectrum_ID)
		AddHeaderColumn(DATA_COLUMN_Charge)
		AddHeaderColumn(DATA_COLUMN_PrecursorMZ)
		AddHeaderColumn(DATA_COLUMN_DelM)
		AddHeaderColumn(DATA_COLUMN_DelM_PPM)
		AddHeaderColumn(DATA_COLUMN_MH)
		AddHeaderColumn(DATA_COLUMN_Peptide)
		AddHeaderColumn(DATA_COLUMN_Protein)
		AddHeaderColumn(DATA_COLUMN_Protein_Mass)
		AddHeaderColumn(DATA_COLUMN_Unexpected_Mod_Count)
		AddHeaderColumn(DATA_COLUMN_Peak_Count)
		AddHeaderColumn(DATA_COLUMN_Matched_Peak_Count)
		AddHeaderColumn(DATA_COLUMN_Matched_Fragment_Ion_Count)
		AddHeaderColumn(DATA_COLUMN_PValue)
		AddHeaderColumn(DATA_COLUMN_Rank_PValue)
		AddHeaderColumn(DATA_COLUMN_EValue)
		AddHeaderColumn(DATA_COLUMN_FDR)

	End Sub

	Protected Function DeterminePrecursorMassTolerance(ByRef objSearchEngineParams As clsSearchEngineParameters) As Double
		Dim strTolerance As String = String.Empty

		Dim dblTolerancePPM As Double
		Dim dblToleranceDa As Double = 0

		If objSearchEngineParams.Parameters.TryGetValue("errorTolerance", strTolerance) Then
			' Parent mass tolerance, in ppm
			If Double.TryParse(strTolerance, dblTolerancePPM) Then
				dblToleranceDa = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000)
			End If
		End If

		Return dblToleranceDa

	End Function

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_FHT
	End Function

	Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msalign_syn_ModSummary.txt"
	End Function

	Public Shared Function GetPHRPProteinModsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msalign_syn_ProteinMods.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_SYN
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msalign_syn_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msalign_syn_SeqInfo.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msalign_syn_SeqToProteinMap.txt"
	End Function

	Public Shared Function GetSearchEngineName() As String
		Return MSAlign_SEARCH_ENGINE_NAME
	End Function

	''' <summary>
	''' Parses the specified MSAlign parameter file
	''' </summary>
	''' <param name="strSearchEngineParamFileName"></param>
	''' <param name="objSearchEngineParams"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

		Dim blnSuccess As Boolean

		objSearchEngineParams = New clsSearchEngineParameters(MSAlign_SEARCH_ENGINE_NAME)

		blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

		ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

		Return blnSuccess

	End Function

	Protected Function ReadSearchEngineParamFile(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
		Dim strSettingValue As String = String.Empty
		Dim objModDef As clsModificationDefinition
		Dim blnSuccess As Boolean

		Try
			blnSuccess = ReadKeyValuePairSearchEngineParamFile(MSAlign_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, objSearchEngineParams)

			If blnSuccess Then
				If objSearchEngineParams.Parameters.TryGetValue("cysteineProtection", strSettingValue) Then

					Select Case strSettingValue
						Case "C57"
							objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 57.0215, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcet")
							objSearchEngineParams.AddModification(objModDef)
						Case "C58"
							objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 58.0055, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcid")
							objSearchEngineParams.AddModification(objModDef)
					End Select

				End If

				' Determine the precursor mass tolerance (will store 0 if a problem or not found)
				objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams)
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
				.ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
				If .ScanNumber = -100 Then
					' Data line is not valid
				Else

					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0)
					.ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_Rank_PValue, mColumnHeaders, 1)

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
					.AddProtein(strProtein)

					dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				UpdatePSMUsingSeqInfo(objPSM)

				' Store the remaining scores
				AddScore(objPSM, strColumns, DATA_COLUMN_Prsm_ID)
				AddScore(objPSM, strColumns, DATA_COLUMN_Spectrum_ID)

				AddScore(objPSM, strColumns, DATA_COLUMN_MH)

				AddScore(objPSM, strColumns, DATA_COLUMN_Protein_Mass)
				AddScore(objPSM, strColumns, DATA_COLUMN_Unexpected_Mod_Count)
				AddScore(objPSM, strColumns, DATA_COLUMN_Peak_Count)
				AddScore(objPSM, strColumns, DATA_COLUMN_Matched_Peak_Count)
				AddScore(objPSM, strColumns, DATA_COLUMN_Matched_Fragment_Ion_Count)

				AddScore(objPSM, strColumns, DATA_COLUMN_PValue)
				AddScore(objPSM, strColumns, DATA_COLUMN_EValue)
				AddScore(objPSM, strColumns, DATA_COLUMN_FDR)

			End If

		Catch ex As Exception
			MyBase.ReportError("Error parsing line " & intLinesRead & " in the MSAlign data file: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

End Class