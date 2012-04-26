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

#Region "Constants"
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

	Protected Const SEQ_SEARCH_ENGINE_NAME As String = "SEQUEST"
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
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.Sequest, blnLoadModsAndSeqInfo)
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

	Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn_ModSummary.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn.txt"
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn_SeqInfo.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_syn_SeqToProteinMap.txt"
	End Function

	''' <summary>
	''' Parses the specified Sequest parameter file
	''' </summary>
	''' <param name="strSearchEngineParamFileName"></param>
	''' <param name="objSearchEngineParams"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
		Dim strParamFilePath As String
		Dim blnSuccess As Boolean
		Dim strLineIn As String

		Dim strFastaFilePath As String
		Dim strSettingValue As String

		Dim intCharIndex As Integer
		Dim intValue As Integer

		Dim kvSetting As System.Collections.Generic.KeyValuePair(Of String, String)

		Dim reEnzymeSpecificity As System.Text.RegularExpressions.Regex = New System.Text.RegularExpressions.Regex("^\S+\s(\d)\s\d\s.+", Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.IgnoreCase)
		Dim reMatch As System.Text.RegularExpressions.Match

		Try
			objSearchEngineParams = New clsSearchEngineParameters(SEQ_SEARCH_ENGINE_NAME, mModInfo)
			objSearchEngineParams.Enzyme = "trypsin"

			strParamFilePath = System.IO.Path.Combine(mInputFolderPath, strSearchEngineParamFileName)

			If Not System.IO.File.Exists(strParamFilePath) Then
				ReportError("Sequest param file not found: " & strParamFilePath)
			Else
				Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strParamFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					While srInFile.Peek > -1
						strLineIn = srInFile.ReadLine().TrimStart()

						If Not String.IsNullOrWhiteSpace(strLineIn) AndAlso Not strLineIn.StartsWith(";") AndAlso Not strLineIn.StartsWith("[") AndAlso strLineIn.Contains("="c) Then

							' Split the line on the equals sign
							kvSetting = ParseKeyValueSetting(strLineIn, "="c)

							If Not String.IsNullOrEmpty(kvSetting.Key) Then

								' Trim off any text that occurs after a semicolon in kvSetting.Value
								strSettingValue = kvSetting.Value
								intCharIndex = strSettingValue.IndexOf(";"c)
								If intCharIndex > 0 Then
									strSettingValue = strSettingValue.Substring(intCharIndex).Trim()
								End If

								objSearchEngineParams.AddUpdateParameter(kvSetting.Key, strSettingValue)

								Select Case kvSetting.Key.ToLower()
									Case "first_database_name", "database_name"
										Try
											strFastaFilePath = System.IO.Path.Combine("C:\Database", System.IO.Path.GetFileName(strSettingValue))
										Catch ex As Exception
											strFastaFilePath = strSettingValue
										End Try
										objSearchEngineParams.FastaFilePath = strFastaFilePath

									Case "mass_type_parent"
										If strSettingValue = "0" Then
											' Average mass
											objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_AVERAGE
										Else
											' Monoisotopic mass
											objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC
										End If

									Case "mass_type_fragment"
										If strSettingValue = "0" Then
											' Average mass
											objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_AVERAGE
										Else
											' Monoisotopic mass
											objSearchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC
										End If

									Case "max_num_internal_cleavage_sites"
										If Integer.TryParse(strSettingValue, intValue) Then
											objSearchEngineParams.MaxNumberInternalCleavages = intValue
										End If

									Case "enzyme_info"
										' Used in new-style sequest parameter files

										' Examples include:
										' Fully-tryptic:     Trypsin(KR) 1 1 KR -
										' Partially-tryptic: Trypsin(KR) 2 1 KR -
										' No-enzyme:         No_Enzyme(-) 0 0 - -
										' 
										objSearchEngineParams.Enzyme = "trypsin"

										If strSettingValue.ToLower().StartsWith("no_enzyme") Then
											objSearchEngineParams.MinNumberTermini = 0
										Else
											' Parse out the cleavage specificity number
											' This is the first number after the closing parenthesis in the above examples
											reMatch = reEnzymeSpecificity.Match(strSettingValue)
											If Not reMatch Is Nothing AndAlso reMatch.Success Then
												If Integer.TryParse(reMatch.Groups(1).Value, intValue) Then
													objSearchEngineParams.MinNumberTermini = intValue
												End If
											End If
										End If

									Case "enzyme_number"
										' Used in old-style sequest parameter files
										If Integer.TryParse(strSettingValue, intValue) Then
											If intValue = 0 Then
												' No-enzyme
												objSearchEngineParams.Enzyme = "trypsin"
												objSearchEngineParams.MinNumberTermini = 0
											Else
												Select Case intValue
													Case 1 : objSearchEngineParams.Enzyme = "trypsin"
													Case 2 : objSearchEngineParams.Enzyme = "trypsin_modified"
													Case 3 : objSearchEngineParams.Enzyme = "Chymotrypsin"
													Case 4 : objSearchEngineParams.Enzyme = "Chymotrypsin__modified"
													Case 5 : objSearchEngineParams.Enzyme = "Clostripain"
													Case 6 : objSearchEngineParams.Enzyme = "Cyanogen_Bromide"
													Case 7 : objSearchEngineParams.Enzyme = "IodosoBenzoate"
													Case 8 : objSearchEngineParams.Enzyme = "Proline_Endopept"
													Case 9 : objSearchEngineParams.Enzyme = "Staph_Protease"
													Case 10 : objSearchEngineParams.Enzyme = "Trypsin_K"
													Case 11 : objSearchEngineParams.Enzyme = "Trypsin_R"
													Case 12 : objSearchEngineParams.Enzyme = "GluC"
													Case 13 : objSearchEngineParams.Enzyme = "LysC"
													Case 14 : objSearchEngineParams.Enzyme = "AspN"
													Case 15 : objSearchEngineParams.Enzyme = "Elastase"
													Case 16 : objSearchEngineParams.Enzyme = "Elastase/Tryp/Chymo"
													Case Else : objSearchEngineParams.Enzyme = "Unknown"
												End Select
												objSearchEngineParams.MinNumberTermini = 2
											End If

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

		Dim dblPrecursorMH As Double		
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
					.ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_RankXc, mColumnHeaders, 1)

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_ChargeState, mColumnHeaders, 0), Short)

					strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Reference, mColumnHeaders)
					.AddProtein(strProtein)

					' Note that the MH value listed in Sequest files is not the precursor MH but is instead the theoretical (computed) MH of the peptide
					' We'll update this value below using dblMassErrorDa
					' We'll further update this value using the ScanStatsEx data
					dblPrecursorMH = LookupColumnValue(strColumns, DATA_COLUMN_MH, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMH, 1, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
					If Double.TryParse(.MassErrorDa, dblMassErrorDa) Then
						' Adjust the precursor mass
						.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMH - dblMassErrorDa, 1, 0)
					End If

					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				UpdatePSMUsingSeqInfo(objPSM)

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
				AddScore(objPSM, strColumns, DATA_COLUMN_NumTrypticEnds)
			End If


		Catch ex As Exception
			MyBase.ReportError("Error parsing line " & intLinesRead & " in the Sequest data file: " & ex.Message)
		End Try

		Return blnSuccess

	End Function

End Class
