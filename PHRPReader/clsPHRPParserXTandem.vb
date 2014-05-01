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

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader
Imports System.IO

Public Class clsPHRPParserXTandem
	Inherits clsPHRPParser

#Region "Constants"
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

	Public Const FILENAME_SUFFIX_SYN As String = "_xt.txt"
	Public Const FILENAME_SUFFIX_FHT As String = "_xt.txt"

	Protected Const XT_SEARCH_ENGINE_NAME As String = "X! Tandem"

	Protected Const TAXONOMY_INFO_KEY_NAME As String = "list path, taxonomy information"
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
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.XTandem, blnLoadModsAndSeqInfo)
	End Sub

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal startupOptions As clsPHRPStartupOptions)
		MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.XTandem, startupOptions)
	End Sub

	Protected Shared Function AppendToString(ByVal strText As String, ByVal strAppend As String) As String

		If String.IsNullOrEmpty(strText) Then
			Return strAppend
		Else
			Return strText & "; " & strAppend
		End If

	End Function

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

	''' <summary>
	''' Determines the precursor mass tolerance
	''' </summary>
	''' <param name="objSearchEngineParams"></param>
	''' <param name="dblTolerancePPM">Precursor mass tolerance, in ppm</param>
	''' <returns>Precursor tolerance, in Da</returns>
	''' <remarks></remarks>
	Protected Function DeterminePrecursorMassTolerance(ByRef objSearchEngineParams As clsSearchEngineParameters, <Out()> ByRef dblTolerancePPM As Double) As Double
		Dim strTolerance As String = String.Empty
		Dim strUnits As String = String.Empty
		Dim blnPPM As Boolean = False

		Dim dblTolerancePlus As Double = 0
		Dim dblToleranceMinus As Double = 0
		Dim dblTolerance As Double

		If objSearchEngineParams.Parameters.TryGetValue("spectrum, parent monoisotopic mass error units", strUnits) Then
			If strUnits.ToLower().Trim() = "ppm" Then blnPPM = True
		End If

		If objSearchEngineParams.Parameters.TryGetValue("spectrum, parent monoisotopic mass error minus", strTolerance) Then
			Double.TryParse(strTolerance, dblToleranceMinus)
		End If

		If objSearchEngineParams.Parameters.TryGetValue("spectrum, parent monoisotopic mass error plus", strTolerance) Then
			Double.TryParse(strTolerance, dblTolerancePlus)
		End If

		dblTolerance = Math.Max(dblToleranceMinus, dblTolerancePlus)
		If blnPPM Then
			dblTolerancePPM = dblTolerance

			' Convert from PPM to dalton (assuming a mass of 2000 m/z)
			dblTolerance = clsPeptideMassCalculator.PPMToMass(dblTolerance, 2000)
		Else
			' Convert from dalton to PPM (assuming a mass of 2000 m/z)
			dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblTolerance, 2000)
		End If

		Return dblTolerance

	End Function

	Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
		' X!Tandem does not have a first-hits file; just the _xt.txt file
		Return String.Empty
	End Function

	Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_ModSummary.txt"
	End Function

	Public Shared Function GetPHRPPepToProteinMapFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_PepToProtMapMTS.txt"
	End Function

	Public Shared Function GetPHRPProteinModsFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & "_xt_ProteinMods.txt"
	End Function

	Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
		Return strDatasetName & FILENAME_SUFFIX_SYN
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

	Public Shared Function GetAdditionalSearchEngineParamFileNames(ByVal strSearchEngineParamFilePath As String) As List(Of String)
		Dim lstFileNames As List(Of String)
		lstFileNames = New List(Of String)

		Dim strDefaultParamsFilename As String
		Dim strTaxonomyFilename As String = String.Empty
		Dim strErrorMessage As String = String.Empty

		Dim fiFileInfo As FileInfo
		Dim objSearchEngineParams As clsSearchEngineParameters

		Try
			If Not File.Exists(strSearchEngineParamFilePath) Then
				lstFileNames.Add("default_input.xml  (Not Confirmed)")
				lstFileNames.Add("taxonomy.xml  (Not Confirmed)")
			Else
				fiFileInfo = New FileInfo(strSearchEngineParamFilePath)
				objSearchEngineParams = New clsSearchEngineParameters(XT_SEARCH_ENGINE_NAME)

				Try
					strDefaultParamsFilename = GetXTandemDefaultParamsFilename(strSearchEngineParamFilePath)
					If Not String.IsNullOrEmpty(strDefaultParamsFilename) Then
						lstFileNames.Add(strDefaultParamsFilename)
					End If
				Catch ex As Exception
					Console.WriteLine("Error in GetXTandemDefaultParamsFilename: " & ex.Message)
					strDefaultParamsFilename = String.Empty
				End Try

				Try

					ParseXTandemParamFileWork(fiFileInfo.DirectoryName, fiFileInfo.Name, objSearchEngineParams, blnDetermineFastaFileNameUsingTaxonomyFile:=False, blnLookForDefaultParamsFileName:=True, strErrorMessage:=strErrorMessage)

					If objSearchEngineParams.Parameters.TryGetValue(TAXONOMY_INFO_KEY_NAME, strTaxonomyFilename) Then
						lstFileNames.Add(Path.GetFileName(strTaxonomyFilename))
					End If

				Catch ex As Exception
					Console.WriteLine("Error in ParseXTandemParamFileWork: " & ex.Message)
				End Try

				If Not String.IsNullOrEmpty(strErrorMessage) Then
					Console.WriteLine(strErrorMessage)
				End If

			End If

		Catch ex As Exception
			Console.WriteLine("Exception in GetAdditionalSearchEngineParamFileNames: " & ex.Message)
		End Try

		Return lstFileNames

	End Function

	Protected Shared Function GetFastaFileFromTaxonomyFile(ByVal strInputFolderPath As String, ByVal strTaxononomyFilename As String, ByRef strErrorMessage As String) As String

		Dim strTaxonomyFilePath As String
		Dim strFastaFile As String = String.Empty

		Dim strFileFormat As String

		Dim kvSetting As KeyValuePair(Of String, String) = Nothing

		Try
			strTaxonomyFilePath = Path.Combine(strInputFolderPath, strTaxononomyFilename)
			If Not File.Exists(strTaxonomyFilePath) Then
				strErrorMessage = AppendToString(strErrorMessage, "Warning, taxonomy file not found: " & strTaxonomyFilePath)
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
			strErrorMessage = AppendToString(strErrorMessage, "Error in GetFastaFileFromTaxonomyFile: " & ex.Message)
		End Try

		Return strFastaFile

	End Function

	Public Shared Function GetSearchEngineName() As String
		Return XT_SEARCH_ENGINE_NAME
	End Function

	Protected Shared Function GetXTandemDefaultParamsFilename(ByVal strParamFilePath As String) As String

		Dim strDefaultParamsFilename As String = String.Empty
		Dim kvSetting As KeyValuePair(Of String, String) = Nothing

		' Open the XML file and look for the "list path, default parameters" entry
		Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(strParamFilePath)

			Do While MoveToNextInputParam(objXMLReader, kvSetting)
				If kvSetting.Key = "list path, default parameters" Then
					strDefaultParamsFilename = String.Copy(kvSetting.Value)
					Exit Do
				End If
			Loop

		End Using

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

	Public Function ParseXTandemParamFile(ByVal strParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters, ByVal blnLookForDefaultParamsFileName As Boolean) As Boolean
		Return ParseXTandemParamFile(strParamFileName, objSearchEngineParams, blnLookForDefaultParamsFileName, blnDetermineFastaFileNameUsingTaxonomyFile:=True)
	End Function

	Public Function ParseXTandemParamFile(ByVal strParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters, ByVal blnLookForDefaultParamsFileName As Boolean, ByVal blnDetermineFastaFileNameUsingTaxonomyFile As Boolean) As Boolean

		Dim strParamFilePath As String
		Dim strErrorMessage As String = String.Empty

		Dim blnSuccess As Boolean

		Try
			strParamFilePath = Path.Combine(mInputFolderPath, strParamFileName)

			If Not File.Exists(strParamFilePath) Then
				ReportError("X!Tandem param file not found: " & strParamFilePath)
			Else

				Try
					blnSuccess = ParseXTandemParamFileWork(mInputFolderPath, strParamFileName, objSearchEngineParams, blnDetermineFastaFileNameUsingTaxonomyFile, blnLookForDefaultParamsFileName, strErrorMessage)
				Catch ex As Exception
					ReportError("Error in ParseXTandemParamFileWork: " & ex.Message)
				End Try

				If Not String.IsNullOrEmpty(strErrorMessage) Then
					ReportError(strErrorMessage)
				End If

				' Determine the precursor mass tolerance (will store 0 if a problem or not found)
				Dim dblTolerancePPM As Double
				objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, dblTolerancePPM)
				objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM

			End If

		Catch ex As Exception
			ReportError("Error in ParseXTandemParamFile: " & ex.Message)
		End Try

		ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

		Return blnSuccess

	End Function

	Protected Shared Function ParseXTandemParamFileWork(ByVal strInputFolderPath As String, ByVal strParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters, ByVal blnDetermineFastaFileNameUsingTaxonomyFile As Boolean, ByVal blnLookForDefaultParamsFileName As Boolean, ByRef strErrorMessage As String) As Boolean

		' Note: Do not put a Try/Catch block in this function

		Dim strParamFilePath As String
		Dim strDefaultParamsFilename As String

		Dim kvSetting As KeyValuePair(Of String, String) = Nothing
		Dim strSetting As String
		Dim intValue As Integer

		Dim blnSuccess As Boolean

		strParamFilePath = Path.Combine(strInputFolderPath, strParamFileName)

		If blnLookForDefaultParamsFileName Then
			Try
				strDefaultParamsFilename = GetXTandemDefaultParamsFilename(strParamFilePath)
			Catch ex As Exception
				strErrorMessage = AppendToString(strErrorMessage, "Error in GetXTandemDefaultParamsFilename: " & ex.Message)
				strDefaultParamsFilename = String.Empty
			End Try

			If Not String.IsNullOrEmpty(strDefaultParamsFilename) Then
				' Read the parameters from the default parameters file and store them in objSearchEngineParams
				' Do this by recursively calling this function

				' First confirm that the file exists
				If File.Exists(Path.Combine(strInputFolderPath, strDefaultParamsFilename)) Then
					ParseXTandemParamFileWork(strInputFolderPath, strDefaultParamsFilename, objSearchEngineParams, blnDetermineFastaFileNameUsingTaxonomyFile, blnLookForDefaultParamsFileName:=False, strErrorMessage:=strErrorMessage)
				Else
					strErrorMessage = AppendToString(strErrorMessage, "Warning, file not found: " & strDefaultParamsFilename)
				End If
			End If
		End If

		' Now read the parameters in strParamFilePath
		Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(strParamFilePath)

			Do While MoveToNextInputParam(objXMLReader, kvSetting)

				If Not String.IsNullOrEmpty(kvSetting.Key) Then
					objSearchEngineParams.AddUpdateParameter(kvSetting.Key, kvSetting.Value)

					Select Case kvSetting.Key
						Case TAXONOMY_INFO_KEY_NAME

							If blnDetermineFastaFileNameUsingTaxonomyFile Then
								' Open the taxonomy file to determine the fasta file used
								strSetting = GetFastaFileFromTaxonomyFile(strInputFolderPath, Path.GetFileName(kvSetting.Value), strErrorMessage)

								If Not String.IsNullOrEmpty(strSetting) Then
									objSearchEngineParams.FastaFilePath = strSetting
								End If
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

		Return blnSuccess

	End Function

	Protected Shared Function MoveToNextInputParam(ByRef objXMLReader As System.Xml.XmlTextReader, ByRef kvParameter As KeyValuePair(Of String, String)) As Boolean

		Dim strNoteType As String
		Dim strParamName As String
		Dim strValue As String

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

							kvParameter = New KeyValuePair(Of String, String)(strParamName, strValue)
							Return True
						End If
					End If
				End If
			End If
		Loop

		Return False

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

		Dim dblPeptideMH As Double
		Dim dblMassErrorDa As Double

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
					.ResultID = LookupColumnValue(strColumns, DATA_COLUMN_Result_ID, mColumnHeaders, 0)

					' X!Tandem only tracks the top-ranked peptide for each spectrum
					.ScoreRank = 1

					strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide_Sequence, mColumnHeaders)
					.SetPeptide(strPeptide, mCleavageStateCalculator)

					.Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

					' Lookup the protein name(s) using mResultIDToProteins
					Dim lstProteinsForResultID As List(Of String) = Nothing
					If mResultIDToProteins.TryGetValue(.ResultID, lstProteinsForResultID) Then
						For Each strProtein As String In lstProteinsForResultID
							.AddProtein(strProtein)
						Next
					End If

					' The Peptide_MH value listed in X!Tandem files is the theoretical (computed) MH of the peptide
					' We'll update this value below using dblMassErrorDa
					' We'll further update this value using the ScanStatsEx data
					dblPeptideMH = LookupColumnValue(strColumns, DATA_COLUMN_Peptide_MH, mColumnHeaders, 0.0#)
					.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPeptideMH, 1, 0)

					.MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_Delta_Mass, mColumnHeaders)
					If Double.TryParse(.MassErrorDa, dblMassErrorDa) Then
						' Adjust the precursor mass
						.PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPeptideMH - dblMassErrorDa, 1, 0)
					End If

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

				' This is the base-10 log of the expectation value
				Dim dblLogEValue As Double
				If Double.TryParse(objPSM.GetScore(DATA_COLUMN_Peptide_Expectation_Value_LogE), dblLogEValue) Then
					' Record the original E-value
					objPSM.SetScore("Peptide_Expectation_Value", Math.Pow(10, dblLogEValue).ToString("0.00e+000"))
				End If

			End If

		Catch ex As Exception
			HandleException("Error parsing line " & intLinesRead & " in the X!Tandem data file", ex)
		End Try

		Return blnSuccess

	End Function

	Private Shared Function XMLTextReaderGetAttributeValue(ByRef objXMLReader As System.Xml.XmlTextReader, ByVal strAttributeName As String, ByVal strValueIfMissing As String) As String
		objXMLReader.MoveToAttribute(strAttributeName)
		If objXMLReader.ReadAttributeValue() Then
			Return objXMLReader.Value
		Else
			Return String.Copy(strValueIfMissing)
		End If
	End Function

	Private Shared Function XMLTextReaderGetInnerText(ByRef objXMLReader As System.Xml.XmlTextReader) As String
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

	Private Shared Sub XMLTextReaderSkipWhitespace(ByRef objXMLReader As System.Xml.XmlTextReader)
		If objXMLReader.NodeType = Xml.XmlNodeType.Whitespace Then
			' Whitspace; read the next node
			objXMLReader.Read()
		End If
	End Sub
End Class
