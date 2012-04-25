'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class parses data lines from Xtandem _xt.txt files
'
' Note: in order to fully extract the search parameters, you need to have these files in the same folder as the _xt.txt file
'	The file passed to LoadSearchEngineParameters()  (typically input.xml)
'	The taxonomy file that it references             (typically taxonomy.xml)
'   The default input file defined in input.xml      (typically default_input.xml)
'*********************************************************************************************************

Option Strict On

Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserXTandem
	Inherits clsPHRPParser

	Public Const DATA_COLUMN_Result_ID As String = "Result_ID"
	Public Const DATA_COLUMN_Group_ID As String = "Group_ID"
	Public Const DATA_COLUMN_Scan As String = "Scan"
	Public Const DATA_COLUMN_Charge As String = "Charge"
	Public Const DATA_COLUMN_Peptide_MH As String = "Peptide_MH"
	Public Const DATA_COLUMN_Peptide_Hyperscore As String = "Peptide_Hyperscore"
	Public Const DATA_COLUMN_Peptide_Expectation_Value_LogE As String = "Peptide_Expectation_Value_Log(e)"
	Public Const DATA_COLUMN_Multiple_Protein_Count As String = "Multiple_Protein_Count"
	Public Const DATA_COLUMN_Peptide_Sequence As String = "Peptide_Sequence"
	Public Const DATA_COLUMN_DeltaCn2 As String = "DeltaCn2"
	Public Const DATA_COLUMN_y_score As String = "y_score"
	Public Const DATA_COLUMN_y_ions As String = "y_ions"
	Public Const DATA_COLUMN_b_score As String = "b_score"
	Public Const DATA_COLUMN_b_ions As String = "b_ions"
	Public Const DATA_COLUMN_Delta_Mass As String = "Delta_Mass"
	Public Const DATA_COLUMN_Peptide_Intensity_LogI As String = "Peptide_Intensity_Log(I)"
	Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"

	Protected Const XT_SEARCH_ENGINE_NAME As String = "X! Tandem"

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.XTandem)
		mInitialized = True
	End Sub

	Protected Overrides Sub DefineColumnHeaders()

		mColumnHeaders.Clear()

		' Define the default column mapping
		AddHeaderColumn(DATA_COLUMN_Result_ID)
		AddHeaderColumn(DATA_COLUMN_Group_ID)
		AddHeaderColumn(DATA_COLUMN_Scan)
		AddHeaderColumn(DATA_COLUMN_Charge)
		AddHeaderColumn(DATA_COLUMN_Peptide_MH)
		AddHeaderColumn(DATA_COLUMN_Peptide_Hyperscore)
		AddHeaderColumn(DATA_COLUMN_Peptide_Expectation_Value_LogE)
		AddHeaderColumn(DATA_COLUMN_Multiple_Protein_Count)
		AddHeaderColumn(DATA_COLUMN_Peptide_Sequence)
		AddHeaderColumn(DATA_COLUMN_DeltaCn2)
		AddHeaderColumn(DATA_COLUMN_y_score)
		AddHeaderColumn(DATA_COLUMN_y_ions)
		AddHeaderColumn(DATA_COLUMN_b_score)
		AddHeaderColumn(DATA_COLUMN_b_ions)
		AddHeaderColumn(DATA_COLUMN_Delta_Mass)
		AddHeaderColumn(DATA_COLUMN_Peptide_Intensity_LogI)
		AddHeaderColumn(DATA_COLUMN_DelM_PPM)

	End Sub

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		' X!Tandem does not have a first-hits file; just the _xt.txt file
		Return String.Empty
	End Function

	Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_ModSummary.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt.txt"
	End Function

	Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_ResultToSeqMap.txt"
	End Function

	Public Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_SeqInfo.txt"
	End Function

	Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_SeqToProteinMap.txt"
	End Function

	Protected Function GetFastaFileFromTaxonomyFile(ByVal strTaxononomyFilename As String) As String

		Dim strTaxonomyFilePath As String
		Dim strFastaFile As String = String.Empty

		Dim strFileFormat As String

		Dim kvSetting As System.Collections.Generic.KeyValuePair(Of String, String) = Nothing

		Try
			strTaxonomyFilePath = System.IO.Path.Combine(mInputFolderPath, strTaxononomyFilename)
			If Not System.IO.File.Exists(strTaxonomyFilePath) Then
				ReportError("Taxonomy file not found: " & strTaxonomyFilePath)
			Else

				' Open the XML file and look for the "file" element with attribute "peptide"
				Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(strTaxonomyFilePath)

					Do While objXMLReader.Read()
						XMLTextReaderSkipWhitespace(objXMLReader)
						If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

						If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
							If objXMLReader.Name.ToLower() = "file" Then

								strFileFormat = XMLTextReaderGetAttributeValue(objXMLReader, "format", String.Empty)

								If strFileFormat = "peptide" Then
									strFastaFile = XMLTextReaderGetAttributeValue(objXMLReader, "URL", String.Empty)
									Exit Do
								End If
							End If
						End If
					Loop

				End Using

			End If
		Catch ex As Exception
			ReportError("Error in GetFastaFileFromTaxonomyFile: " & ex.Message)
		End Try

		Return strFastaFile

	End Function

	Protected Function GetXTandemDefaultParamsFilename(ByVal strParamFilePath As String) As String

		Dim strDefaultParamsFilename As String = String.Empty
		Dim kvSetting As System.Collections.Generic.KeyValuePair(Of String, String) = Nothing

		Try

			' Open the XML file and look for the "list path, default parameters" entry
			Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(strParamFilePath)

				Do While MoveToNextInputParam(objXMLReader, kvSetting)
					If kvSetting.Key = "list path, default parameters" Then
						strDefaultParamsFilename = String.Copy(kvSetting.Value)
						Exit Do
					End If
				Loop

			End Using

		Catch ex As Exception
			ReportError("Error in GetXTandemDefaultParamsFilename: " & ex.Message)
		End Try

		Return strDefaultParamsFilename

	End Function

	''' <summary>
	''' Parses the specified X!Tandem parameter file
	''' Note that the file specified by parameter "list path, default parameters" will also be auto-parsed (if found in folder mInputFolderPath)
	''' </summary>
	''' <param name="strSearchEngineParamFileName"></param>
	''' <param name="objSearchEngineParams"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

		objSearchEngineParams = New clsSearchEngineParameters(XT_SEARCH_ENGINE_NAME, mModInfo)

		Return ParseXTandemParamFile(strSearchEngineParamFileName, objSearchEngineParams, blnLookForDefaultParamsFileName:=True)

	End Function

	Public Function ParseXTandemParamFile(ByVal strParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters, blnLookForDefaultParamsFileName As Boolean) As Boolean

		Dim strParamFilePath As String
		Dim strDefaultParamsFilename As String

		Dim kvSetting As System.Collections.Generic.KeyValuePair(Of String, String) = Nothing
		Dim strSetting As String
		Dim intValue As Integer

		Dim blnSuccess As Boolean

		Try

			strParamFilePath = System.IO.Path.Combine(mInputFolderPath, strParamFileName)

			If Not System.IO.File.Exists(strParamFilePath) Then
				ReportError("X!Tandem param file not found: " & strParamFilePath)
			Else

				If blnLookForDefaultParamsFileName Then
					strDefaultParamsFilename = GetXTandemDefaultParamsFilename(strParamFilePath)

					If Not String.IsNullOrEmpty(strDefaultParamsFilename) Then
						' Read the parameters from the default parameters file and store them in objSearchEngineParams
						' Do this by recursively calling this function
						ParseXTandemParamFile(strDefaultParamsFilename, objSearchEngineParams, blnLookForDefaultParamsFileName:=False)
					End If
				End If

				' Now read the parameters in strParamFilePath
				Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(strParamFilePath)

					Do While MoveToNextInputParam(objXMLReader, kvSetting)

						If Not String.IsNullOrEmpty(kvSetting.Key) Then
							objSearchEngineParams.AddUpdateParameter(kvSetting.Key, kvSetting.Value)

							Select Case kvSetting.Key
								Case "list path, taxonomy information"
									' Open the taxonomy file to determine the fasta file used
									strSetting = GetFastaFileFromTaxonomyFile(kvSetting.Value)

									If Not String.IsNullOrEmpty(strSetting) Then
										objSearchEngineParams.FastaFilePath = strSetting
									End If

								Case "spectrum, fragment mass type"
									objSearchEngineParams.FragmentMassType = kvSetting.Value

								Case "scoring, maximum missed cleavage sites"
									If Integer.TryParse(kvSetting.Value, intValue) Then
										objSearchEngineParams.MaxNumberInternalCleavages = intValue
									End If

							End Select
						End If
					Loop

				End Using

				blnSuccess = True

			End If
		Catch ex As Exception
			ReportError("Error in LoadSearchEngineParameters: " & ex.Message)
		End Try

		Return blnSuccess
	End Function

	Protected Function MoveToNextInputParam(ByRef objXMLReader As System.Xml.XmlTextReader, ByRef kvParameter As System.Collections.Generic.KeyValuePair(Of String, String)) As Boolean

		Dim strNoteType As String
		Dim strParamName As String
		Dim strValue As String

		Dim blnMatchFound As Boolean
		blnMatchFound = False

		Do While objXMLReader.Read()
			XMLTextReaderSkipWhitespace(objXMLReader)
			If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

			If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
				If objXMLReader.Name.ToLower() = "note" Then

					strNoteType = XMLTextReaderGetAttributeValue(objXMLReader, "type", String.Empty)

					If strNoteType = "input" Then

						strParamName = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)

						If objXMLReader.Read() Then
							' Read the note's inner text
							strValue = XMLTextReaderGetInnerText(objXMLReader)

							kvParameter = New System.Collections.Generic.KeyValuePair(Of String, String)(strParamName, strValue)

						End If
					End If
				End If
			End If
		Loop

		Return blnMatchFound

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
		Dim lstProteinsForResultID As System.Collections.Generic.List(Of String) = Nothing

		Dim dblPrecursorMH As Double

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
					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_Result_ID, mColumnHeaders, 0)

					' X!Tandem only tracks the top-ranked peptide for each spectrum
					.ScoreRank = 1

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide_Sequence, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					' Lookup the protein name(s) using mResultIDToProteins
					If mResultIDToProteins.TryGetValue(.ResultID, lstProteinsForResultID) Then
						For Each strProtein As String In lstProteinsForResultID
							.AddProtein(strProtein)
						Next
					End If

					dblPrecursorMH = LookupColumnValue(strColumns, DATA_COLUMN_Peptide_MH, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMH, 1, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_Delta_Mass, mColumnHeaders)
					.MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

					blnSuccess = True
				End If
			End With

			If blnSuccess Then
				UpdatePSMUsingSeqInfo(objPSM)

				' Store the remaining scores
				AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Hyperscore)
				AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Expectation_Value_LogE)
				AddScore(objPSM, strColumns, DATA_COLUMN_DeltaCn2)
				AddScore(objPSM, strColumns, DATA_COLUMN_y_score)
				AddScore(objPSM, strColumns, DATA_COLUMN_y_ions)
				AddScore(objPSM, strColumns, DATA_COLUMN_b_score)
				AddScore(objPSM, strColumns, DATA_COLUMN_b_ions)
				AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Intensity_LogI)
			End If

		Catch ex As Exception
			HandleException("Error parsing line " & intLinesRead & " in the X!Tandem data file", ex)
		End Try

		Return blnSuccess

	End Function

	Private Function XMLTextReaderGetAttributeValue(ByRef objXMLReader As System.Xml.XmlTextReader, ByVal strAttributeName As String, ByVal strValueIfMissing As String) As String
		objXMLReader.MoveToAttribute(strAttributeName)
		If objXMLReader.ReadAttributeValue() Then
			Return objXMLReader.Value
		Else
			Return String.Copy(strValueIfMissing)
		End If
	End Function

	Private Function XMLTextReaderGetInnerText(ByRef objXMLReader As System.Xml.XmlTextReader) As String
		Dim strValue As String = String.Empty
		Dim blnSuccess As Boolean

		If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
			' Advance the reader so that we can read the value
			blnSuccess = objXMLReader.Read()
		Else
			blnSuccess = True
		End If

		If blnSuccess AndAlso Not objXMLReader.NodeType = Xml.XmlNodeType.Whitespace And objXMLReader.HasValue Then
			strValue = objXMLReader.Value
		End If

		Return strValue
	End Function

	Private Sub XMLTextReaderSkipWhitespace(ByRef objXMLReader As System.Xml.XmlTextReader)
		If objXMLReader.NodeType = Xml.XmlNodeType.Whitespace Then
			' Whitspace; read the next node
			objXMLReader.Read()
		End If
	End Sub
End Class
