'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class is the base class for classes used to parse PHRP data lines
' It must be derived by a sub-class customized for the specific analysis tool (Sequest, X!Tandem, Inspect, etc.)
'
'*********************************************************************************************************

Option Strict On

Public MustInherit Class clsPHRPParser

#Region "Module variables"

	Protected mDatasetName As String
	Protected mInputFilePath As String
	Protected mInputFolderPath As String
	Protected mInitialized As Boolean

	' Column headers in the synopsis file and first hits file
	Protected mColumnHeaders As System.Collections.Generic.SortedDictionary(Of String, Integer)

	Protected mErrorMessage As String = String.Empty

	Protected mCleavageStateCalculator As clsPeptideCleavageStateCalculator
	Protected mPeptideMassCalculator As clsPeptideMassCalculator

	Protected mPeptideHitResultType As clsPHRPReader.ePeptideHitResultType

	Protected mModInfo As System.Collections.Generic.List(Of clsModificationDefinition)
	Protected mModInfoLoaded As Boolean

	Protected mResultToSeqMap As System.Collections.Generic.SortedList(Of Integer, Integer)
	Protected mSeqInfo As System.Collections.Generic.SortedList(Of Integer, clsSeqInfo)
	Protected mSeqToProteinMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of clsProteinInfo))

	' This List tracks the Protein Names for each ResultID
	Protected mResultIDToProteins As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))

	Protected mErrorMessages As System.Collections.Generic.List(Of String)
	Protected mWarningMessages As System.Collections.Generic.List(Of String)

#End Region

#Region "Events"
	Public Event MessageEvent(ByVal strMessage As String)
	Public Event ErrorEvent(ByVal strErrorMessage As String)
	Public Event WarningEvent(ByVal strWarningMessage As String)
#End Region

#Region "Properties"

	''' <summary>
	''' Cached error messages
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property ErrorMessages() As System.Collections.Generic.List(Of String)
		Get
			Return mErrorMessages
		End Get
	End Property

	''' <summary>
	''' Input file path
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property InputFilePath As String
		Get
			Return mInputFilePath
		End Get
	End Property

	''' <summary>
	''' Input folder path
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property InputFolderPath As String
		Get
			Return mInputFolderPath
		End Get
	End Property

	''' <summary>
	''' Peptide hit result type; Sequest, XTandem, Inspect, or MSGFDB
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property PeptideHitResultType As clsPHRPReader.ePeptideHitResultType
		Get
			Return mPeptideHitResultType
		End Get
	End Property

	''' <summary>
	''' Returns the cached mapping between ResultID and SeqID
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property ResultToSeqMap() As System.Collections.Generic.SortedList(Of Integer, Integer)
		Get
			Return mResultToSeqMap
		End Get
	End Property

	''' <summary>
	''' Returns the cached sequence info, where key is SeqID
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property SeqInfo() As System.Collections.Generic.SortedList(Of Integer, clsSeqInfo)
		Get
			Return mSeqInfo
		End Get
	End Property

	''' <summary>
	''' Returns the cached sequence to protein map information
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property SeqToProteinMap() As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of clsProteinInfo))
		Get
			Return mSeqToProteinMap
		End Get
	End Property

	''' <summary>
	''' Cached warning messages
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property WarningMessages() As System.Collections.Generic.List(Of String)
		Get
			Return mWarningMessages
		End Get
	End Property

#End Region

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset Name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, ByVal blnLoadModsAndSeqInfo As Boolean)

		mErrorMessages = New System.Collections.Generic.List(Of String)
		mWarningMessages = New System.Collections.Generic.List(Of String)

		mDatasetName = strDatasetName
		mPeptideHitResultType = ePeptideHitResultType

		Dim fiFileInfo As System.IO.FileInfo = New System.IO.FileInfo(strInputFilePath)
		mInputFilePath = fiFileInfo.FullName
		mInputFolderPath = fiFileInfo.DirectoryName

		mErrorMessage = String.Empty

		' Initialize the column mapping object
		' Using a case-insensitive comparer
		mColumnHeaders = New System.Collections.Generic.SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

		mCleavageStateCalculator = New clsPeptideCleavageStateCalculator()
		mPeptideMassCalculator = New clsPeptideMassCalculator()

		' Initialize the tracking lists
		mResultToSeqMap = New System.Collections.Generic.SortedList(Of Integer, Integer)
		mSeqInfo = New System.Collections.Generic.SortedList(Of Integer, clsSeqInfo)
		mSeqToProteinMap = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of clsProteinInfo))

		mResultIDToProteins = New System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))

		If blnLoadModsAndSeqInfo Then
			' Read the ModSummary file (if it exists)
			mModInfoLoaded = LoadModSummary()
		Else
			mModInfoLoaded = False
		End If

		If mModInfoLoaded AndAlso blnLoadModsAndSeqInfo Then
			' Read the ResultToSeqMapInfo (if the files exist)
			LoadSeqInfo()
		End If

		' The following will be overridden by a derived form of this class
		DefineColumnHeaders()

	End Sub

#Region "Functions overridden by derived classes"
	Protected MustOverride Sub DefineColumnHeaders()

	''' <summary>
	''' Parse the data line read from a PHRP results file
	''' </summary>
	''' <param name="strLine">Data line</param>
	''' <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
	''' <param name="objPSM">clsPSM object (output)</param>
	''' <returns>True if success, false if an error</returns>
	Public MustOverride Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean

	''' <summary>
	''' Parses the specified parameter file
	''' </summary>
	''' <param name="strSearchEngineParamFileName">Name of the parameter file to parse (must reside in InputFolderPath)</param>
	''' <param name="objSearchEngineParams">Search engine parameters class (output)</param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public MustOverride Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean
#End Region

	Protected Sub AddHeaderColumn(ByVal strColumnName As String)
		mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
	End Sub

	Protected Sub AddScore(ByRef objPSM As clsPSM, ByRef strColumns() As String, ByVal strScoreColumnName As String)
		Const NOT_FOUND As String = "==SCORE_NOT_FOUND=="

		Dim strValue As String
		strValue = PHRPReader.clsPHRPReader.LookupColumnValue(strColumns, strScoreColumnName, mColumnHeaders, NOT_FOUND)

		If strValue <> NOT_FOUND Then
			objPSM.SetScore(strScoreColumnName, strValue)
		End If

	End Sub

	''' <summary>
	''' Clear any cached error messages
	''' </summary>
	''' <remarks></remarks>
	Public Sub ClearErrors()
		mErrorMessages.Clear()
	End Sub

	''' <summary>
	''' Clear any cached warning messages
	''' </summary>
	''' <remarks></remarks>
	Public Sub ClearWarnings()
		mWarningMessages.Clear()
	End Sub

	Protected Function ConvertModsToNumericMods(ByVal strCleanSequence As String, ByRef lstModifiedResidues As List(Of clsAminoAcidModInfo)) As String
		Static sbNewPeptide As New System.Text.StringBuilder

		sbNewPeptide.Length = 0

		If lstModifiedResidues Is Nothing OrElse lstModifiedResidues.Count = 0 Then
			Return strCleanSequence
		End If

		For intIndex = 0 To strCleanSequence.Length - 1
			sbNewPeptide.Append(strCleanSequence.Chars(intIndex))

			For Each objModInfo In lstModifiedResidues
				If objModInfo.ResidueLocInPeptide = intIndex + 1 Then
					If objModInfo.ModDefinition.ModificationMass < 0 Then
						sbNewPeptide.Append(objModInfo.ModDefinition.ModificationMassAsText)
					Else
						sbNewPeptide.Append("+" & objModInfo.ModDefinition.ModificationMassAsText)
					End If
				End If
			Next
		Next

		Return sbNewPeptide.ToString()

	End Function

	Protected Sub HandleException(ByVal strBaseMessage As String, ByVal ex As System.Exception)
		If String.IsNullOrEmpty(strBaseMessage) Then
			strBaseMessage = "Error"
		End If

		ReportError(strBaseMessage & ": " & ex.Message)
	End Sub

	''' <summary>
	''' Reads the data in strModSummaryFilePath.  Populates mModInfo with the modification names, masses, and affected residues
	''' </summary>
	''' <returns>True if success; false if an error</returns>
	Protected Function LoadModSummary() As Boolean

		Dim objModSummaryReader As clsPHRPModSummaryReader

		Dim strModSummaryFilePath As String

		Dim blnSuccess As Boolean

		Try
			strModSummaryFilePath = clsPHRPReader.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName)
			If String.IsNullOrEmpty(strModSummaryFilePath) Then
				ReportWarning("ModSummaryFile path is empty; unable to continue")
				Return False
			End If

			strModSummaryFilePath = System.IO.Path.Combine(mInputFolderPath, strModSummaryFilePath)
			If Not System.IO.File.Exists(strModSummaryFilePath) Then
				ReportWarning("ModSummary file not found: " & strModSummaryFilePath)
				Return False
			End If

			ShowMessage("Reading the PHRP ModSummary file")

			objModSummaryReader = New clsPHRPModSummaryReader(strModSummaryFilePath)
			blnSuccess = objModSummaryReader.Success

			If blnSuccess Then
				mModInfo = objModSummaryReader.ModificationDefs
			End If

		Catch ex As Exception
			HandleException("Exception reading PHRP Mod Summary file", ex)
			Return False
		End Try

		Return blnSuccess

	End Function

	Protected Function LoadSeqInfo() As Boolean

		Dim blnSuccess As Boolean
		Dim objReader As clsPHRPSeqMapReader

		Try

			ShowMessage("Reading the PHRP SeqInfo file")

			' Instantiate the reader
			objReader = New clsPHRPSeqMapReader(mDatasetName, mInputFolderPath, mPeptideHitResultType)

			' Read the files
			blnSuccess = objReader.GetProteinMapping(mResultToSeqMap, mSeqToProteinMap, mSeqInfo)

			If Not blnSuccess Then
				ReportWarning(objReader.ErrorMessage)
			End If

			mResultIDToProteins.Clear()

			If blnSuccess Then
				' Populate mResultIDToProteins

				Dim intResultID As Integer
				Dim intSeqID As Integer

				For Each objItem As System.Collections.Generic.KeyValuePair(Of Integer, Integer) In mResultToSeqMap

					intResultID = objItem.Key
					intSeqID = objItem.Value

					Dim lstProteinsForSeqID As System.Collections.Generic.List(Of clsProteinInfo) = Nothing
					Dim lstProteinsForResultID As System.Collections.Generic.List(Of String) = Nothing
					lstProteinsForResultID = New System.Collections.Generic.List(Of String)

					If mSeqToProteinMap.TryGetValue(intSeqID, lstProteinsForSeqID) Then

						For Each objProtein As clsProteinInfo In lstProteinsForSeqID
							If Not lstProteinsForResultID.Contains(objProtein.ProteinName) Then
								lstProteinsForResultID.Add(objProtein.ProteinName)
							End If
						Next

					End If

					mResultIDToProteins.Add(intResultID, lstProteinsForResultID)

				Next

			End If

		Catch ex As Exception
			HandleException("Error loading PHRP Seq Info", ex)
			blnSuccess = False
			If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
		End Try

		Return blnSuccess

	End Function

	''' <summary>
	''' Parse the column names in strSplitLine and update the local column header mapping
	''' </summary>
	''' <param name="strSplitLine"></param>
	''' <remarks></remarks>
	Public Sub ParseColumnHeaders(ByRef strSplitLine() As String)
		clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders)
	End Sub

	''' <summary>
	''' Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
	''' </summary>
	''' <param name="strText"></param>
	''' <param name="chDelimiter"></param>
	''' <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
	''' <remarks>Automatically trims whitespace</remarks>
	Protected Function ParseKeyValueSetting(ByVal strText As String, ByVal chDelimiter As Char) As System.Collections.Generic.KeyValuePair(Of String, String)
		Dim strSplitLine() As String
		Dim kvSetting As System.Collections.Generic.KeyValuePair(Of String, String)

		If Not String.IsNullOrEmpty(strText) Then
			strSplitLine = strText.Split(chDelimiter)

			If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
				If strSplitLine.Length = 1 Then
					kvSetting = New System.Collections.Generic.KeyValuePair(Of String, String)(strSplitLine(0).Trim(), String.Empty)
				Else
					kvSetting = New System.Collections.Generic.KeyValuePair(Of String, String)(strSplitLine(0).Trim(), strSplitLine(1).Trim())
				End If
				Return kvSetting
			End If
		End If

		Return New System.Collections.Generic.KeyValuePair(Of String, String)(String.Empty, String.Empty)

	End Function

	Protected Sub ReportError(ByVal strErrorMessage As String)
		mErrorMessage = strErrorMessage
		mErrorMessages.Add(strErrorMessage)
		RaiseEvent ErrorEvent(strErrorMessage)
	End Sub

	Protected Sub ReportWarning(ByVal strWarningMessage As String)
		mWarningMessages.Add(strWarningMessage)
		RaiseEvent WarningEvent(strWarningMessage)
	End Sub

	Protected Sub ShowMessage(ByVal strMessage As String)
		RaiseEvent MessageEvent(strMessage)
	End Sub

	''' <summary>
	''' Updates the theoretical (computed) monoisotopic mass of objPSM using mResultToSeqMap and mSeqInfo
	''' Also updates the modification info
	''' Also updates SeqID
	''' </summary>
	''' <param name="objPSM"></param>
	''' <returns>True if success, False if objPSM.ResultID is not found in mResultToSeqMap</returns>
	''' <remarks></remarks>
	Protected Function UpdatePSMUsingSeqInfo(ByRef objPSM As clsPSM) As Boolean
		Dim intSeqID As Integer
		Dim objSeqInfo As clsSeqInfo = Nothing

		Dim strMods() As String
		Dim kvModDetails As System.Collections.Generic.KeyValuePair(Of String, String)

		Dim strMassCorrectionTag As String
		Dim intResidueLoc As Integer
		Dim intPeptideResidueCount As Integer

		Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

		Dim blnMatchFound As Boolean

		Dim blnFavorTerminalMods As Boolean
		Dim blnSuccess As Boolean

		blnSuccess = False

		' First determine the modified residues present in this peptide
		If Not mResultToSeqMap Is Nothing AndAlso mResultToSeqMap.Count > 0 Then
			If mResultToSeqMap.TryGetValue(objPSM.ResultID, intSeqID) Then

				objPSM.SeqID = intSeqID
				intPeptideResidueCount = objPSM.PeptideCleanSequence.Length

				If mSeqInfo.TryGetValue(intSeqID, objSeqInfo) Then
					objPSM.PeptideMonoisotopicMass = objSeqInfo.MonoisotopicMass

					objPSM.ClearModifiedResidues()

					If objSeqInfo.ModCount > 0 Then
						' Split objSeqInfo.ModDescription on the comma character
						strMods = objSeqInfo.ModDescription.Split(","c)

						If Not strMods Is Nothing AndAlso strMods.Count > 0 Then
							For intModIndex As Integer = 0 To strMods.Count - 1

								' Split strMods on the colon characters
								kvModDetails = ParseKeyValueSetting(strMods(intModIndex), ":"c)

								If Not String.IsNullOrEmpty(kvModDetails.Key) AndAlso Not String.IsNullOrEmpty(kvModDetails.Value) Then
									strMassCorrectionTag = kvModDetails.Key
									If Integer.TryParse(kvModDetails.Value, intResidueLoc) Then
										' Find the modification definition in mModInfo
										' Note that a given mass correction tag might be present multiple times in mModInfo, since it could be used as both a static peptide mod and a static peptide terminal mod
										' Thus, if intResidueLoc = 1 or intResidueLoc = objPSM.PeptideCleanSequence.Length then we'll first look for a peptide or protein terminal static mod

										If intResidueLoc = 1 Then
											eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
											blnFavorTerminalMods = True
										ElseIf intResidueLoc = intPeptideResidueCount Then
											eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
											blnFavorTerminalMods = True
										Else
											eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
											blnFavorTerminalMods = False
										End If

										Dim objMatchedModDef As clsModificationDefinition = Nothing
										blnMatchFound = UpdatePSMFindMatchingModInfo(strMassCorrectionTag, blnFavorTerminalMods, eResidueTerminusState, objMatchedModDef)

										If blnMatchFound Then
											objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence.Chars(intResidueLoc - 1), intResidueLoc, eResidueTerminusState, objMatchedModDef)
										Else
											' Could not find a valid entry in mModInfo
											ReportError("Unrecognized mass correction tag found in the SeqInfo file: " & strMassCorrectionTag)
										End If

									End If
								End If
							Next intModIndex
						End If

					End If

					blnSuccess = True
				End If
			End If
		End If

		If blnSuccess Then
			Dim strPrimarySequence As String = String.Empty
			Dim strPrefix As String = String.Empty
			Dim strSuffix As String = String.Empty

			If clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(objPSM.Peptide, strPrimarySequence, strPrefix, strSuffix) Then
				objPSM.PeptideWithNumericMods = strPrefix & "." & ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues) & "." & strSuffix
			Else
				objPSM.PeptideWithNumericMods = ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues)
			End If

		End If

		Return blnSuccess
	End Function

	Protected Function UpdatePSMFindMatchingModInfo( _
	  ByVal strMassCorrectionTag As String, _
	  ByVal blnFavorTerminalMods As Boolean, _
	  ByVal eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
	  ByRef objMatchedModDef As clsModificationDefinition) As Boolean

		Dim blnMatchFound As Boolean

		blnMatchFound = False
		Do
			For Each objMod In mModInfo
				If strMassCorrectionTag = objMod.MassCorrectionTag Then
					If blnFavorTerminalMods Then
						If objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse _
						   objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod Then

							If eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus AndAlso _
							  (objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
								blnMatchFound = True
							End If

							If eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus AndAlso _
							  (objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
								blnMatchFound = True
							End If

						End If
					Else
						blnMatchFound = True
					End If

					If blnMatchFound Then
						' Match found
						objMatchedModDef = objMod
						Exit For
					End If
				End If
			Next

			If Not blnMatchFound AndAlso blnFavorTerminalMods Then
				' Change FavorTerminalMods to false and try again
				blnFavorTerminalMods = False
			Else
				Exit Do
			End If
		Loop While Not blnMatchFound

		Return blnMatchFound

	End Function
End Class
