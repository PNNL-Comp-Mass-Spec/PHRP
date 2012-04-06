Option Strict On

Public Class clsPHRPSeqMapReader

#Region "Constants"
	Public Const SEQ_PROT_MAP_COLUMN_Unique_Seq_ID As String = "Unique_Seq_ID"
	Public Const SEQ_PROT_MAP_COLUMN_Cleavage_State As String = "Cleavage_State"
	Public Const SEQ_PROT_MAP_COLUMN_Terminus_State As String = "Terminus_State"
	Public Const SEQ_PROT_MAP_COLUMN_Protein_Name As String = "Protein_Name"

	Public Const SEQ_PROT_MAP_COLUMN_Protein_EValue As String = "Protein_Expectation_Value_Log(e)"		' Only used by X!Tandem
	Public Const SEQ_PROT_MAP_COLUMN_Protein_Intensity As String = "Protein_Intensity_Log(I)"			' Only used by X!Tandem
#End Region

#Region "Module-wide variables"
	Protected mDatasetName As String
	Protected mInputFolderPath As String

	Protected mResultToSeqMapFilename As String
	Protected mSeqToProteinMapFilename As String

	Protected mPeptideHitResultType As clsPHRPReader.ePeptideHitResultType
#End Region

#Region "Properties"

	Public ReadOnly Property DatasetName As String
		Get
			Return mDatasetName
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
		mDatasetName = strDatasetName
		If String.IsNullOrEmpty(mDatasetName) Then
			Throw New Exception("Dataset name cannot be empty")
		End If

		mInputFolderPath = strInputFolderPath
		If String.IsNullOrEmpty(mInputFolderPath) Then
			mInputFolderPath = String.Empty
		End If

		mPeptideHitResultType = ePeptideHitResultType

		mResultToSeqMapFilename = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName)
		If String.IsNullOrEmpty(mResultToSeqMapFilename) Then
			Throw New Exception("Unable to determine ResultToSeqMap filename for PeptideHitResultType: " & mPeptideHitResultType.ToString())
		End If

		mSeqToProteinMapFilename = clsPHRPReader.GetPHRPSeqToProteinMapFileName(mPeptideHitResultType, mDatasetName)
		If String.IsNullOrEmpty(mSeqToProteinMapFilename) Then
			Throw New Exception("Unable to determine SeqToProteinMap filename for PeptideHitResultType: " & mPeptideHitResultType.ToString())
		End If

	End Sub

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strInputFolderPath">Input folder path</param>
	''' <param name="strResultToSeqMapFilename">ResultToSeqMap filename</param>
	''' <param name="strSeqToProteinMapFilename">SeqToProteinMap filename</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strInputFolderPath As String, ByVal strResultToSeqMapFilename As String, ByVal strSeqToProteinMapFilename As String)
		mInputFolderPath = strInputFolderPath
		If String.IsNullOrEmpty(mInputFolderPath) Then
			mInputFolderPath = String.Empty
		End If

		mPeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strResultToSeqMapFilename)
		If mPeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Unknown Then
			Throw New Exception("Unable to auto-determine the PepthideHit result type based on filename " & strResultToSeqMapFilename)
		End If

		mDatasetName = clsPHRPReader.AutoDetermineDatasetName(strResultToSeqMapFilename)
		If String.IsNullOrEmpty(mDatasetName) Then
			Throw New Exception("Dataset name cannot be empty")
		End If

		mResultToSeqMapFilename = strResultToSeqMapFilename
		mSeqToProteinMapFilename = strSeqToProteinMapFilename

	End Sub

	''' <summary>
	''' Load the mapping between ResultID and Protein Name
	''' </summary>
	''' <param name="objResultToSeqMap">ResultToSeqMap list (output); keys are ResultID, Values as SeqID</param>
	''' ''' <param name="objSeqToProteinMap">SeqToProteinMap list (output); keys are SeqID, Values are Protein Name</param>
	''' <returns>True if success, false if an error</returns>
	''' <remarks></remarks>
	Public Function GetProteinMapping( _
	 ByRef objResultToSeqMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer)), _
	 ByRef objSeqToProteinMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))) As Boolean

		Dim blnSuccess As Boolean

		' Note: do not put a Try/Catch handler in this function
		'       Instead, allow LoadResultToSeqMapping or LoadSeqToProteinMapping to raise exceptions

		If objResultToSeqMap Is Nothing Then
			objResultToSeqMap = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))
		Else
			objResultToSeqMap.Clear()
		End If

		If objSeqToProteinMap Is Nothing Then
			objSeqToProteinMap = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))
		Else
			objSeqToProteinMap.Clear()
		End If

		If String.IsNullOrEmpty(mResultToSeqMapFilename) Then
			blnSuccess = False
		Else
			blnSuccess = LoadResultToSeqMapping(System.IO.Path.Combine(mInputFolderPath, mResultToSeqMapFilename), objResultToSeqMap)
		End If

		If blnSuccess AndAlso Not String.IsNullOrEmpty(mSeqToProteinMapFilename) Then
			blnSuccess = LoadSeqToProteinMapping(System.IO.Path.Combine(mInputFolderPath, mSeqToProteinMapFilename), objSeqToProteinMap)
		End If

		Return blnSuccess

	End Function

	''' <summary>
	''' Load the Result to Seq mapping using the specified PHRP result file
	''' </summary>
	''' <param name="strFilePath"></param>
	''' <param name="objResultToSeqMap"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Function LoadResultToSeqMapping(ByVal strFilePath As String, ByRef objResultToSeqMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of Integer))) As Boolean

		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim intResultID As Integer
		Dim intSeqID As Integer

		Dim lstSeqIDs As System.Collections.Generic.List(Of Integer) = Nothing

		Try

			' Read the data from the result to sequence map file
			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))

				Do While srInFile.Peek >= 0
					strLineIn = srInFile.ReadLine

					If Not String.IsNullOrEmpty(strLineIn) Then
						strSplitLine = strLineIn.Split(ControlChars.Tab)

						If strSplitLine.Length >= 2 Then

							' Parse out the numbers from the first two columns 
							' (the first line of the file is the header line, and it will get skipped)
							If Integer.TryParse(strSplitLine(0), intResultID) Then
								If Integer.TryParse(strSplitLine(1), intSeqID) Then

									If objResultToSeqMap.TryGetValue(intResultID, lstSeqIDs) Then
										If Not lstSeqIDs.Contains(intSeqID) Then
											lstSeqIDs.Add(intSeqID)
											objResultToSeqMap(intResultID) = lstSeqIDs
										End If
									Else
										lstSeqIDs = New System.Collections.Generic.List(Of Integer)
										lstSeqIDs.Add(intSeqID)
										objResultToSeqMap.Add(intResultID, lstSeqIDs)
									End If

								End If
							End If

						End If
					End If
				Loop

			End Using


		Catch ex As Exception
			Throw New Exception("Exception loading Result to Seq Mapping from " & System.IO.Path.GetFileName(strFilePath) & ": " & ex.Message)
		End Try

		Return True

	End Function

	''' <summary>
	''' Load the Sequence to Protein mapping using the specified PHRP result file
	''' </summary>
	''' <param name="strFilePath"></param>
	''' <param name="objSeqToProteinMap"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Function LoadSeqToProteinMapping(ByVal strFilePath As String, _
	  ByRef objSeqToProteinMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))) As Boolean

		Dim lstProteins As System.Collections.Generic.List(Of String) = Nothing

		Dim objColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)

		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim strProtein As String

		Dim intSeqID As Integer

		Dim blnHeaderLineParsed As Boolean
		Dim blnSkipLine As Boolean

		Try

			' Initialize the column mapping
			' Using a case-insensitive comparer
			objColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

			' Define the default column mapping
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, 0)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Cleavage_State, 1)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Terminus_State, 2)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_Name, 3)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_EValue, 4)
			objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_Intensity, 5)

			' Read the data from the sequence to protein map file
			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))

				Do While srInFile.Peek >= 0
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

								strProtein = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Protein_Name, objColumnHeaders, String.Empty)

								If Not String.IsNullOrEmpty(strProtein) Then

									If objSeqToProteinMap.TryGetValue(intSeqID, lstProteins) Then
										If Not lstProteins.Contains(strProtein) Then
											lstProteins.Add(strProtein)
											objSeqToProteinMap(intSeqID) = lstProteins
										End If
									Else
										lstProteins = New System.Collections.Generic.List(Of String)
										lstProteins.Add(strProtein)
										objSeqToProteinMap.Add(intSeqID, lstProteins)
									End If

								End If

							End If

						End If

					End If
				Loop

			End Using

		Catch ex As Exception
			Throw New Exception("Exception loading Seq to Protein Mapping from " & System.IO.Path.GetFileName(strFilePath) & ": " & ex.Message)
		End Try

		Return True

	End Function



End Class
