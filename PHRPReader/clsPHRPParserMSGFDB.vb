'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class parses data lines from MSGFDB msgfdb_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader
Imports System.Text.RegularExpressions

Public Class clsPHRPParserMSGFDB
	Inherits clsPHRPParser

#Region "Constants"
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

	Public Const DATA_COLUMN_MSGFDB_SpecProb As String = "MSGFDB_SpecProb"					' MSGFDB
	Public Const DATA_COLUMN_Rank_MSGFDB_SpecProb As String = "Rank_MSGFDB_SpecProb"		' MSGFDB

	Public Const DATA_COLUMN_MSGFDB_SpecEValue As String = "MSGFDB_SpecEValue"				' MSGF+
	Public Const DATA_COLUMN_Rank_MSGFDB_SpecEValue As String = "Rank_MSGFDB_SpecEValue"	' MSGF+

	Public Const DATA_COLUMN_PValue As String = "PValue"		' MSGFDB
	Public Const DATA_COLUMN_EValue As String = "EValue"		' MSGF+

	Public Const DATA_COLUMN_FDR As String = "FDR"							' MSGFDB; Only present if a Target/Decoy (TDA) search was used
	Public Const DATA_COLUMN_PepFDR As String = "PepFDR"					' MSGFDB; Only valid if a Target/Decoy (TDA) search was used; if EFDR is present, will contain 1 for every row

	Public Const DATA_COLUMN_QValue As String = "QValue"					' MSGF+ reports QValue instead of FDR
	Public Const DATA_COLUMN_PepQValue As String = "PepQValue"				' MSGF+ reports pepQValue instead of PepFDR

	Public Const DATA_COLUMN_EFDR As String = "EFDR"						' Only present if a Target/Decoy (TDA) search was not used

	Public Const DATA_COLUMN_IMS_Scan As String = "IMS_Scan"
	Public Const DATA_COLUMN_IMS_Drift_Time As String = "IMS_Drift_Time"

	Public Const DATA_COLUMN_Isotope_Error As String = "IsotopeError"		' Only reported by MSGF+

	Public Const FILENAME_SUFFIX_SYN As String = "_msgfdb_syn.txt"
	Public Const FILENAME_SUFFIX_FHT As String = "_msgfdb_fht.txt"

	Protected Const MSGFDB_SEARCH_ENGINE_NAME As String = "MS-GFDB"
#End Region

#Region "Properties"

	Public Overrides ReadOnly Property PHRPFirstHitsFileName() As String
		Get
			Return GetPHRPFirstHitsFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPModSummaryFileName() As String
		Get
			Return GetPHRPModSummaryFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPPepToProteinMapFileName() As String
		Get
			Return GetPHRPPepToProteinMapFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPProteinModsFileName() As String
		Get
			Return GetPHRPProteinModsFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPSynopsisFileName() As String
		Get
			Return GetPHRPSynopsisFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPResultToSeqMapFileName() As String
		Get
			Return GetPHRPResultToSeqMapFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPSeqInfoFileName() As String
		Get
			Return GetPHRPSeqInfoFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property PHRPSeqToProteinMapFileName() As String
		Get
			Return GetPHRPSeqToProteinMapFileName(mDatasetName)
		End Get
	End Property

	Public Overrides ReadOnly Property SearchEngineName() As String
		Get
			Return GetSearchEngineName()
		End Get
	End Property

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
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSGFDB, blnLoadModsAndSeqInfo)
	End Sub

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal startupOptions As clsPHRPStartupOptions)
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSGFDB, startupOptions)
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
		AddHeaderColumn(DATA_COLUMN_EFDR)
		AddHeaderColumn(DATA_COLUMN_PepFDR)

		' Add the MSGF+ columns
		AddHeaderColumn(DATA_COLUMN_MSGFDB_SpecEValue)
		AddHeaderColumn(DATA_COLUMN_Rank_MSGFDB_SpecEValue)
		AddHeaderColumn(DATA_COLUMN_EValue)

		AddHeaderColumn(DATA_COLUMN_QValue)
		AddHeaderColumn(DATA_COLUMN_PepQValue)

		AddHeaderColumn(DATA_COLUMN_IMS_Scan)
		AddHeaderColumn(DATA_COLUMN_IMS_Drift_Time)
		AddHeaderColumn(DATA_COLUMN_Isotope_Error)

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
		Dim strToleranceSplit As String()

		Dim reExtraTolerance As Regex
		Dim reMatch As Match
		reExtraTolerance = New Regex("([0-9.]+)([A-Za-z]+)", RegexOptions.Compiled Or RegexOptions.IgnoreCase)

		Dim dblToleranceDa As Double
		Dim dblToleranceCurrent As Double

		dblTolerancePPM = 0

		If objSearchEngineParams.Parameters.TryGetValue("PMTolerance", strTolerance) Then
			' Parent mass tolerance
			' Might contain two values, separated by a comma
			strToleranceSplit = strTolerance.Split(","c)

			If Not strToleranceSplit Is Nothing Then
				For Each strItem As String In strToleranceSplit
					If Not strItem.Trim.StartsWith("#") Then
						reMatch = reExtraTolerance.Match(strItem)

						If reMatch.Success Then
							If Double.TryParse(reMatch.Groups(1).Value, dblToleranceCurrent) Then
								If reMatch.Groups(2).Value.ToLower().Contains("ppm") Then
									' Ppm
									' Convert from PPM to dalton (assuming a mass of 2000 m/z)
									dblToleranceCurrent = clsPeptideMassCalculator.PPMToMass(dblToleranceCurrent, 2000)
								End If

								dblToleranceDa = Math.Max(dblToleranceDa, dblToleranceCurrent)
							End If
						End If
					End If
				Next
			End If

		End If


		Return dblToleranceDa

	End Function

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_FHT
	End Function

	Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_ModSummary.txt"
	End Function

	Public Shared Function GetPHRPPepToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_PepToProtMapMTS.txt"
	End Function

	Public Shared Function GetPHRPProteinModsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_ProteinMods.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_SYN
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_SeqInfo.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_SeqToProteinMap.txt"
	End Function

	Public Shared Function GetSearchEngineName() As String
		Return MSGFDB_SEARCH_ENGINE_NAME
	End Function

	''' <summary>
	''' Parses the specified MSGFDB (aka MS-GF+) parameter file
	''' </summary>
	''' <param name="strSearchEngineParamFileName"></param>
	''' <param name="objSearchEngineParams"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

		Dim blnSuccess As Boolean

		objSearchEngineParams = New clsSearchEngineParameters(MSGFDB_SEARCH_ENGINE_NAME, mModInfo)

		blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

		ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

		Return blnSuccess

	End Function

	Protected Function ReadSearchEngineParamFile(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
		Dim strSettingValue As String = String.Empty
		Dim intValue As Integer
		Dim blnSuccess As Boolean

		Try
			blnSuccess = ReadKeyValuePairSearchEngineParamFile(MSGFDB_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, ePeptideHitResultType.MSGFDB, objSearchEngineParams)

			If blnSuccess Then

				' Determine the enzyme name
				If objSearchEngineParams.Parameters.TryGetValue("enzymeid", strSettingValue) Then
					If Integer.TryParse(strSettingValue, intValue) Then
						Select Case intValue
							Case 0 : objSearchEngineParams.Enzyme = "no_enzyme"
							Case 1 : objSearchEngineParams.Enzyme = "trypsin"
							Case 2 : objSearchEngineParams.Enzyme = "Chymotrypsin"
							Case 3 : objSearchEngineParams.Enzyme = "Lys-C "
							Case 4 : objSearchEngineParams.Enzyme = "Lys-N  "
							Case 5 : objSearchEngineParams.Enzyme = "Glu-C  "
							Case 6 : objSearchEngineParams.Enzyme = "Arg-C "
							Case 7 : objSearchEngineParams.Enzyme = "Asp-N"
						End Select
					End If
				End If

				' Determine the cleavage specificity
				If objSearchEngineParams.Parameters.TryGetValue("nnet", strSettingValue) Then
					' NNET means number of non-enzymatic terminii

					If Integer.TryParse(strSettingValue, intValue) Then
						Select Case intValue
							Case 0
								' Fully-tryptic
								objSearchEngineParams.MinNumberTermini = 2
							Case 1
								' Partially-tryptic
								objSearchEngineParams.MinNumberTermini = 1
							Case Else
								' No-enzyme search
								objSearchEngineParams.MinNumberTermini = 0
						End Select

					End If
				Else
					' MSGF+ uses ntt instead of nnet; thus look for ntt

					If objSearchEngineParams.Parameters.TryGetValue("ntt", strSettingValue) Then
						' NTT means number of tolerable terminii

						If Integer.TryParse(strSettingValue, intValue) Then
							Select Case intValue
								Case 0
									' No-enzyme search
									objSearchEngineParams.MinNumberTermini = 0
								Case 1
									' Partially-tryptic
									objSearchEngineParams.MinNumberTermini = 1
								Case Else
									' Fully-tryptic
									objSearchEngineParams.MinNumberTermini = 2
							End Select

						End If
					End If

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
	''' <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
	''' <returns>True if success, false if an error</returns>
	''' <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
	Public Overrides Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM, ByVal fastReadMode As Boolean) As Boolean

		Dim strColumns() As String = strLine.Split(ControlChars.Tab)
		Dim strPeptide As String
		Dim strProtein As String

		Dim dblPrecursorMZ As Double
		Dim dblSpecProb As Double

		Dim blnMSGFPlusResults As Boolean
		Dim blnSuccess As Boolean

		Try

			If LookupColumnIndex(DATA_COLUMN_MSGFDB_SpecEValue, mColumnHeaders) >= 0 Then
				blnMSGFPlusResults = True
			Else
				blnMSGFPlusResults = False
			End If

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
					.ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb, mColumnHeaders, 1)

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)

					If fastReadMode Then
						.SetPeptide(strPeptide, blnUpdateCleanSequence:=False)
					Else
						.SetPeptide(strPeptide, mCleavageStateCalculator)
					End If

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
					.AddProtein(strProtein)

					.CollisionMode = LookupColumnValue(strColumns, DATA_COLUMN_FragMethod, mColumnHeaders, "n/a")

					dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					If blnMSGFPlusResults Then
						.MSGFSpecProb = LookupColumnValue(strColumns, DATA_COLUMN_MSGFDB_SpecEValue, mColumnHeaders)
					Else
						.MSGFSpecProb = LookupColumnValue(strColumns, DATA_COLUMN_MSGFDB_SpecProb, mColumnHeaders)
					End If

					If .MSGFSpecProb.Length > 13 Then
						' Attempt to shorten the SpecProb value
						If Double.TryParse(.MSGFSpecProb, dblSpecProb) Then
							.MSGFSpecProb = dblSpecProb.ToString("0.0000000E-00")
						End If

					End If
					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				If Not fastReadMode Then
					UpdatePSMUsingSeqInfo(objPSM)
				End If

				' Store the remaining scores
				AddScore(objPSM, strColumns, DATA_COLUMN_DeNovoScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_MSGFScore)

				If blnMSGFPlusResults Then

					AddScore(objPSM, strColumns, DATA_COLUMN_MSGFDB_SpecEValue)
					AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFDB_SpecEValue)
					AddScore(objPSM, strColumns, DATA_COLUMN_EValue)
					AddScore(objPSM, strColumns, DATA_COLUMN_QValue)
					AddScore(objPSM, strColumns, DATA_COLUMN_PepQValue)
					AddScore(objPSM, strColumns, DATA_COLUMN_Isotope_Error)

					' Duplicate the score values to provide backwards compatibility
					Dim strValue As String = String.Empty

					If objPSM.TryGetScore(DATA_COLUMN_MSGFDB_SpecEValue, strValue) Then objPSM.SetScore(DATA_COLUMN_MSGFDB_SpecProb, strValue)
					If objPSM.TryGetScore(DATA_COLUMN_Rank_MSGFDB_SpecEValue, strValue) Then objPSM.SetScore(DATA_COLUMN_Rank_MSGFDB_SpecProb, strValue)
					If objPSM.TryGetScore(DATA_COLUMN_QValue, strValue) Then objPSM.SetScore(DATA_COLUMN_FDR, strValue)
					If objPSM.TryGetScore(DATA_COLUMN_PepQValue, strValue) Then objPSM.SetScore(DATA_COLUMN_PepFDR, strValue)

					Dim strEValue As String = String.Empty
					Dim dblEValue As Double

					Dim strSpecEValue As String = String.Empty
					Dim dblSpecEValue As Double

					Dim blnPValueStored As Boolean = False

					If objPSM.TryGetScore(DATA_COLUMN_EValue, strEValue) Then
						If objPSM.TryGetScore(DATA_COLUMN_MSGFDB_SpecEValue, strSpecEValue) Then
							' Compute PValue using EValue and SpecEValue
							If Double.TryParse(strEValue, dblEValue) Then
								If Double.TryParse(strSpecEValue, dblSpecEValue) Then
									If dblSpecEValue > 0 Then
										Dim dblN As Double = dblEValue / dblSpecEValue
										Dim dblPValue As Double = 1 - (1 - dblSpecEValue) ^ dblN

										If Math.Abs(dblPValue) <= Double.Epsilon Then
											objPSM.SetScore(DATA_COLUMN_PValue, "0")
										Else
											objPSM.SetScore(DATA_COLUMN_PValue, dblPValue.ToString("0.00000E-00"))
										End If

										blnPValueStored = True
									End If
								End If
							End If
						End If

						If Not blnPValueStored Then
							' Store E-value as P-value (these values are not identical, and will only be close for high-confidence results, i.e. results with FDR < 2%)
							objPSM.SetScore(DATA_COLUMN_PValue, strEValue)
						End If

					End If

				Else
					AddScore(objPSM, strColumns, DATA_COLUMN_MSGFDB_SpecProb)
					AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb)
					AddScore(objPSM, strColumns, DATA_COLUMN_PValue)
					AddScore(objPSM, strColumns, DATA_COLUMN_FDR)
					AddScore(objPSM, strColumns, DATA_COLUMN_PepFDR)
				End If

				AddScore(objPSM, strColumns, DATA_COLUMN_EFDR)		' This column will not be present if a Target/Decoy (TDA) search was performed

				AddScore(objPSM, strColumns, DATA_COLUMN_IMS_Scan)
				AddScore(objPSM, strColumns, DATA_COLUMN_IMS_Drift_Time)

			End If

		Catch ex As Exception
			MyBase.ReportError("Error parsing line " & intLinesRead & " in the MSGFDB data file: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

End Class
