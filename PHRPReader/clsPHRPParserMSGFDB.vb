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

Imports PHRPReader.clsPHRPReader

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
	Public Const DATA_COLUMN_MSGFDB_SpecProb As String = "MSGFDB_SpecProb"
	Public Const DATA_COLUMN_Rank_MSGFDB_SpecProb As String = "Rank_MSGFDB_SpecProb"
	Public Const DATA_COLUMN_PValue As String = "PValue"
	Public Const DATA_COLUMN_FDR As String = "FDR"
	Public Const DATA_COLUMN_PepFDR As String = "PepFDR"

	Protected Const MSGFDB_SEARCH_ENGINE_NAME As String = "MS-GF+"
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
		mInitialized = True
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

	Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn_ModSummary.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_msgfdb_syn.txt"
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

	''' <summary>
	''' Parses the specified MSGFDB (aka MS-GF+) parameter file
	''' </summary>
	''' <param name="strSearchEngineParamFileName"></param>
	''' <param name="objSearchEngineParams"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
		Dim strParamFilePath As String
		Dim blnSuccess As Boolean

		Dim strLineIn As String
		Dim strSettingValue As String

		Dim intCharIndex As Integer
		Dim intValue As Integer

		Dim kvSetting As System.Collections.Generic.KeyValuePair(Of String, String)

		Try
			objSearchEngineParams = New clsSearchEngineParameters(MSGFDB_SEARCH_ENGINE_NAME, mModInfo)

			strParamFilePath = System.IO.Path.Combine(mInputFolderPath, strSearchEngineParamFileName)

			If Not System.IO.File.Exists(strParamFilePath) Then
				ReportError("MSGF DB param file not found: " & strParamFilePath)
			Else
				Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strParamFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					While srInFile.Peek > -1
						strLineIn = srInFile.ReadLine().TrimStart()

						If Not String.IsNullOrWhiteSpace(strLineIn) AndAlso Not strLineIn.StartsWith("#") AndAlso strLineIn.Contains("="c) Then

							' Split the line on the equals sign
							kvSetting = ParseKeyValueSetting(strLineIn, "="c)

							If Not String.IsNullOrEmpty(kvSetting.Key) Then

								' Trim off any text that occurs after a # in kvSetting.Value
								strSettingValue = kvSetting.Value
								intCharIndex = strSettingValue.IndexOf("#"c)
								If intCharIndex > 0 Then
									strSettingValue = strSettingValue.Substring(intCharIndex).Trim()
								End If

								objSearchEngineParams.AddUpdateParameter(kvSetting.Key, strSettingValue)

								Select Case kvSetting.Key.ToLower()
									Case "enzymeid"

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

									Case "nnet"
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
								End Select
							End If
						End If

					End While
				End Using

				blnSuccess = True

			End If
		Catch ex As Exception
			ReportError("Error in LoadSearchEngineParameters: " & ex.Message)
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
		Dim dblSpecProb As Double

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
					.ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb, mColumnHeaders, 1)

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
					.AddProtein(strProtein)

					.CollisionMode = LookupColumnValue(strColumns, DATA_COLUMN_FragMethod, mColumnHeaders, "n/a")

					dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					.MSGFSpecProb = LookupColumnValue(strColumns, DATA_COLUMN_MSGFDB_SpecProb, mColumnHeaders)

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
				UpdatePSMUsingSeqInfo(objPSM)

				' Store the remaining scores
				AddScore(objPSM, strColumns, DATA_COLUMN_DeNovoScore)
				AddScore(objPSM, strColumns, DATA_COLUMN_MSGFScore)
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
