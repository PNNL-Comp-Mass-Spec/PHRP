Option Strict On

Imports PHRPReader.clsPeptideCleavageStateCalculator
Imports System.IO

Public Class clsPHRPSeqMapReader

#Region "Constants"
	Public Const SEQ_PROT_MAP_COLUMN_Unique_Seq_ID As String = "Unique_Seq_ID"
	Public Const SEQ_PROT_MAP_COLUMN_Cleavage_State As String = "Cleavage_State"
	Public Const SEQ_PROT_MAP_COLUMN_Terminus_State As String = "Terminus_State"
	Public Const SEQ_PROT_MAP_COLUMN_Protein_Name As String = "Protein_Name"

	Public Const SEQ_PROT_MAP_COLUMN_Protein_EValue As String = "Protein_Expectation_Value_Log(e)"		' Only used by X!Tandem
	Public Const SEQ_PROT_MAP_COLUMN_Protein_Intensity As String = "Protein_Intensity_Log(I)"			' Only used by X!Tandem

	Public Const SEQ_INFO_COLUMN_Unique_Seq_ID As String = "Unique_Seq_ID"
	Public Const SEQ_INFO_COLUMN_Mod_Count As String = "Mod_Count"
	Public Const SEQ_INFO_COLUMN_Mod_Description As String = "Mod_Description"
	Public Const SEQ_INFO_COLUMN_Monoisotopic_Mass As String = "Monoisotopic_Mass"

#End Region

#Region "Module-wide variables"
	Protected mDatasetName As String
	Protected mInputFolderPath As String

	Protected mResultToSeqMapFilename As String
	Protected mSeqToProteinMapFilename As String
	Protected mSeqInfoFilename As String

	Protected mPeptideHitResultType As clsPHRPReader.ePeptideHitResultType

	Protected mErrorMessage As String = String.Empty
#End Region

#Region "Properties"

	Public ReadOnly Property DatasetName As String
		Get
			Return mDatasetName
		End Get
	End Property

	Public ReadOnly Property ErrorMessage As String
		Get
			Return mErrorMessage
		End Get
	End Property
	Public ReadOnly Property InputFolderPath As String
		Get
			Return mInputFolderPath
		End Get
	End Property

	Public ReadOnly Property PeptideHitResultType As clsPHRPReader.ePeptideHitResultType
		Get
			Return mPeptideHitResultType
		End Get
	End Property

	Public ReadOnly Property ResultToSeqMapFilename As String
		Get
			Return mResultToSeqMapFilename
		End Get
	End Property

	Public ReadOnly Property SeqToProteinMapFilename As String
		Get
			Return mSeqToProteinMapFilename
		End Get
	End Property

#End Region

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFolderPath">Input file path</param>
	''' <param name="ePeptideHitResultType">Peptide Hit result type</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFolderPath As String, ByVal ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType)
		Me.New(strDatasetName, strInputFolderPath, ePeptideHitResultType, PHRPReader.clsPHRPReader.GetPHRPSynopsisFileName(ePeptideHitResultType, strDatasetName))
	End Sub

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset name</param>
	''' <param name="strInputFolderPath">Input file path</param>
	''' <param name="ePeptideHitResultType">Peptide Hit result type</param>
	''' <param name="strPHRPDataFileName">The base PHRP data file name; used when calling AutoSwitchToFHTIfRequired</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFolderPath As String, ByVal ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, ByVal strPHRPDataFileName As String)
		mDatasetName = strDatasetName

		If String.IsNullOrEmpty(mDatasetName) Then
			mErrorMessage = "Dataset name cannot be empty"
			Throw New Exception(mErrorMessage)
		End If

		mInputFolderPath = strInputFolderPath
		If String.IsNullOrEmpty(mInputFolderPath) Then
			mInputFolderPath = String.Empty
		End If

		mPeptideHitResultType = ePeptideHitResultType

		mResultToSeqMapFilename = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName)
		If String.IsNullOrEmpty(mResultToSeqMapFilename) Then
			mErrorMessage = "Unable to determine ResultToSeqMap filename for PeptideHitResultType: " & mPeptideHitResultType.ToString()
			Throw New Exception(mErrorMessage)
		Else
			mResultToSeqMapFilename = IO.Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mResultToSeqMapFilename, strPHRPDataFileName))
		End If

		mSeqToProteinMapFilename = clsPHRPReader.GetPHRPSeqToProteinMapFileName(mPeptideHitResultType, mDatasetName)
		If String.IsNullOrEmpty(mSeqToProteinMapFilename) Then
			mErrorMessage = "Unable to determine SeqToProteinMap filename for PeptideHitResultType: " & mPeptideHitResultType.ToString()
			Throw New Exception(mErrorMessage)
		Else
			mSeqToProteinMapFilename = IO.Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mSeqToProteinMapFilename, strPHRPDataFileName))
		End If

		mSeqInfoFilename = clsPHRPReader.GetPHRPSeqInfoFileName(mPeptideHitResultType, mDatasetName)
		If String.IsNullOrEmpty(mSeqInfoFilename) Then
			mErrorMessage = "Unable to determine SeqInfo filename for PeptideHitResultType: " & mPeptideHitResultType.ToString()
			Throw New Exception(mErrorMessage)
		Else
			mSeqInfoFilename = IO.Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mSeqInfoFilename, strPHRPDataFileName))
		End If

	End Sub

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strInputFolderPath">Input folder path</param>
	''' <param name="strResultToSeqMapFilename">ResultToSeqMap filename</param>
	''' <param name="strSeqInfoFilename">SeqInfo filename</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strInputFolderPath As String, ByVal strResultToSeqMapFilename As String, ByVal strSeqToProteinMapFilename As String, ByVal strSeqInfoFilename As String)
		mInputFolderPath = strInputFolderPath

		If String.IsNullOrEmpty(mInputFolderPath) Then
			mInputFolderPath = String.Empty
		End If

		mPeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strResultToSeqMapFilename)
		If mPeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Unknown Then
			mErrorMessage = "Unable to auto-determine the PepthideHit result type based on filename " & strResultToSeqMapFilename
			Throw New Exception(mErrorMessage)
		End If

		mDatasetName = clsPHRPReader.AutoDetermineDatasetName(strResultToSeqMapFilename)
		If String.IsNullOrEmpty(mDatasetName) Then
			mErrorMessage = "Unable to auto-determine the dataset name using filename '" & strResultToSeqMapFilename & "'"
			Throw New Exception(mErrorMessage)
		End If

		mResultToSeqMapFilename = strResultToSeqMapFilename
		mSeqToProteinMapFilename = strSeqToProteinMapFilename
		mSeqInfoFilename = strSeqInfoFilename

	End Sub

	''' <summary>
	''' Load the mapping between ResultID and Protein Name
	''' </summary>
	''' <param name="lstResultToSeqMap">ResultToSeqMap list (output); keys are ResultID, Values as SeqID</param>
	''' <param name="lstSeqToProteinMap">SeqToProteinMap list (output); keys are SeqID, Values are list of clsProteinInfo objects</param>
	''' <param name="lstSeqInfo">SeqInfo list (output); keys are SeqID, Values are seq details stored in clsSeqInfo objects</param>
	''' <returns>True if success, false if an error</returns>
	''' <remarks></remarks>
	Public Function GetProteinMapping( _
	  ByRef lstResultToSeqMap As SortedList(Of Integer, Integer), _
	  ByRef lstSeqToProteinMap As SortedList(Of Integer, List(Of clsProteinInfo)), _
	  ByRef lstSeqInfo As SortedList(Of Integer, clsSeqInfo)) As Boolean

		Dim blnSuccess As Boolean
		Dim strFilePath As String

		' Note: do not put a Try/Catch handler in this function
		'       Instead, allow LoadResultToSeqMapping or LoadSeqToProteinMapping to raise exceptions

		If lstResultToSeqMap Is Nothing Then
			lstResultToSeqMap = New SortedList(Of Integer, Integer)
		Else
			lstResultToSeqMap.Clear()
		End If

		If lstSeqToProteinMap Is Nothing Then
			lstSeqToProteinMap = New SortedList(Of Integer, List(Of clsProteinInfo))
		Else
			lstSeqToProteinMap.Clear()
		End If

		If lstSeqInfo Is Nothing Then
			lstSeqInfo = New SortedList(Of Integer, clsSeqInfo)
		Else
			lstSeqInfo.Clear()
		End If

		If String.IsNullOrEmpty(mResultToSeqMapFilename) Then
			blnSuccess = False
		Else
			strFilePath = Path.Combine(mInputFolderPath, mResultToSeqMapFilename)
			If Not File.Exists(strFilePath) Then
				mErrorMessage = "SeqInfo file not found: " & strFilePath
				blnSuccess = False
			Else
				blnSuccess = LoadResultToSeqMapping(strFilePath, lstResultToSeqMap)
			End If
		End If

		If blnSuccess Then

			If Not String.IsNullOrEmpty(mSeqInfoFilename) Then
				strFilePath = Path.Combine(mInputFolderPath, mSeqInfoFilename)
				If File.Exists(strFilePath) Then
					LoadSeqInfo(strFilePath, lstSeqInfo)
				End If
			End If

			If Not String.IsNullOrEmpty(mSeqToProteinMapFilename) Then
				strFilePath = Path.Combine(mInputFolderPath, mSeqToProteinMapFilename)
				If Not File.Exists(strFilePath) Then
					mErrorMessage = "SeqInfo file not found: " & strFilePath
					blnSuccess = False
				Else
					blnSuccess = LoadSeqToProteinMapping(strFilePath, lstSeqToProteinMap)
				End If
			Else
				blnSuccess = False
			End If

		End If

		Return blnSuccess

	End Function

	''' <summary>
	''' Load the Result to Seq mapping using the specified PHRP result file
	''' </summary>
	''' <param name="strFilePath"></param>
	''' <param name="lstResultToSeqMap"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Function LoadResultToSeqMapping(ByVal strFilePath As String, ByRef lstResultToSeqMap As SortedList(Of Integer, Integer)) As Boolean

		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim intResultID As Integer
		Dim intSeqID As Integer

		Try

			' Read the data from the result to sequence map file
			Using srInFile As StreamReader = New StreamReader(New FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))

				Do While srInFile.Peek > -1
					strLineIn = srInFile.ReadLine

					If Not String.IsNullOrEmpty(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If strSplitLine.Length >= 2 Then

							' Parse out the numbers from the first two columns 
							' (the first line of the file is the header line, and it will get skipped)
							If Integer.TryParse(strSplitLine(0), intResultID) Then
								If Integer.TryParse(strSplitLine(1), intSeqID) Then

									If Not lstResultToSeqMap.ContainsKey(intResultID) Then
										lstResultToSeqMap.Add(intResultID, intSeqID)
									End If
								End If
							End If

						End If
					End If
				Loop

			End Using


		Catch ex As Exception
			Throw New Exception("Exception loading Result to Seq Mapping from " & Path.GetFileName(strFilePath) & ": " & ex.Message)
		End Try

		Return True

	End Function

	Protected Function LoadSeqInfo(ByVal strFilePath As String, ByRef lstSeqInfo As SortedList(Of Integer, clsSeqInfo)) As Boolean

		Dim objColumnHeaders As SortedDictionary(Of String, Integer)

		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim intSeqID As Integer
		Dim intModCount As Integer
		Dim strModDescription As String
		Dim dblMonoisotopicMass As Double

		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Try

			' Initialize the column mapping
			' Using a case-insensitive comparer
			objColumnHeaders = New SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

			' Define the default column mapping
			objColumnHeaders.Add(SEQ_INFO_COLUMN_Unique_Seq_ID, 0)
			objColumnHeaders.Add(SEQ_INFO_COLUMN_Mod_Count, 1)
			objColumnHeaders.Add(SEQ_INFO_COLUMN_Mod_Description, 2)
			objColumnHeaders.Add(SEQ_INFO_COLUMN_Monoisotopic_Mass, 3)

			' Read the data from the sequence info file
			Using srInFile As StreamReader = New StreamReader(New FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))

				Do While srInFile.Peek > -1
					strLineIn = srInFile.ReadLine
					blnSkipLine = False

					If Not String.IsNullOrEmpty(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If Not blnHeaderLineParsed Then
							If strSplitLine(0).ToLower() = SEQ_INFO_COLUMN_Unique_Seq_ID.ToLower() Then
								' Parse the header line to confirm the column ordering
								clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders)
								blnSkipLine = True
							End If

							blnHeaderLineParsed = True
						End If

						If Not blnSkipLine AndAlso strSplitLine.Length >= 3 Then

							If Integer.TryParse(strSplitLine(0), intSeqID) Then

								intModCount = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_INFO_COLUMN_Mod_Count, objColumnHeaders, 0)
								strModDescription = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_INFO_COLUMN_Mod_Description, objColumnHeaders, String.Empty)
								dblMonoisotopicMass = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_INFO_COLUMN_Monoisotopic_Mass, objColumnHeaders, 0.0#)

								If Not lstSeqInfo.ContainsKey(intSeqID) Then
									lstSeqInfo.Add(intSeqID, New clsSeqInfo(intSeqID, dblMonoisotopicMass, intModCount, strModDescription))
								End If

							End If

						End If

					End If
				Loop

			End Using

		Catch ex As Exception
			Throw New Exception("Exception loading Seq Info from " & Path.GetFileName(strFilePath) & ": " & ex.Message)
		End Try

		Return True


	End Function

	''' <summary>
	''' Load the Sequence to Protein mapping using the specified PHRP result file
	''' </summary>
	''' <param name="strFilePath"></param>
	''' <param name="lstSeqToProteinMap"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Function LoadSeqToProteinMapping(ByVal strFilePath As String, _
	  ByRef lstSeqToProteinMap As SortedList(Of Integer, List(Of clsProteinInfo))) As Boolean

		Dim lstProteins As List(Of clsProteinInfo) = Nothing

		Dim objColumnHeaders As SortedDictionary(Of String, Integer)

		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim objProteinInfo As clsProteinInfo

		Dim strProteinName As String
		Dim intSeqID As Integer
		Dim eCleavageState As ePeptideCleavageStateConstants
		Dim eTerminusState As ePeptideTerminusStateConstants

		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Try

			' Initialize the column mapping
			' Using a case-insensitive comparer
			objColumnHeaders = New SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

			' Define the default column mapping
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, 0)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Cleavage_State, 1)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Terminus_State, 2)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_Name, 3)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_EValue, 4)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_Intensity, 5)

			' Read the data from the sequence to protein map file
			Using srInFile As StreamReader = New StreamReader(New FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))

				Do While srInFile.Peek > -1
					strLineIn = srInFile.ReadLine
					blnSkipLine = False

					If Not String.IsNullOrEmpty(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If Not blnHeaderLineParsed Then
							If strSplitLine(0).ToLower() = SEQ_PROT_MAP_COLUMN_Unique_Seq_ID.ToLower() Then
								' Parse the header line to confirm the column ordering
								clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders)
								blnSkipLine = True
							End If

							blnHeaderLineParsed = True
						End If

						If Not blnSkipLine AndAlso strSplitLine.Length >= 3 Then

							If Integer.TryParse(strSplitLine(0), intSeqID) Then

								strProteinName = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Protein_Name, objColumnHeaders, String.Empty)

								If Not String.IsNullOrEmpty(strProteinName) Then

									eCleavageState = CType(clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Cleavage_State, objColumnHeaders, 0), ePeptideCleavageStateConstants)
									eTerminusState = CType(clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Terminus_State, objColumnHeaders, 0), ePeptideTerminusStateConstants)

									objProteinInfo = New clsProteinInfo(strProteinName, intSeqID, eCleavageState, eTerminusState)

									If lstSeqToProteinMap.TryGetValue(intSeqID, lstProteins) Then
										' Sequence already exists in lstSeqToProteinMap; add the new protein info
										lstProteins.Add(objProteinInfo)
									Else
										' New Sequence ID
										lstProteins = New List(Of clsProteinInfo)
										lstProteins.Add(objProteinInfo)
										lstSeqToProteinMap.Add(intSeqID, lstProteins)
									End If

								End If

							End If

						End If

					End If
				Loop

			End Using

		Catch ex As Exception
			Throw New Exception("Exception loading Seq to Protein Mapping from " & Path.GetFileName(strFilePath) & ": " & ex.Message)
		End Try

		Return True

	End Function

End Class
