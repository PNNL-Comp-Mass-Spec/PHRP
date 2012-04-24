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
Imports PeptideHitResultsProcessor.clsPeptideModificationContainer

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

	Protected mModInfo As System.Collections.Generic.List(Of PeptideHitResultsProcessor.clsModificationDefinition)
	Protected mModInfoLoaded As Boolean

	Protected mResultToSeqMap As System.Collections.Generic.SortedList(Of Integer, Integer)
	Protected mSeqInfo As System.Collections.Generic.SortedList(Of Integer, clsSeqInfo)
	Protected mSeqToProteinMap As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of clsProteinInfo))

	' This List tracks the Protein Names for each ResultID
	Protected mResultIDToProteins As System.Collections.Generic.SortedList(Of Integer, System.Collections.Generic.List(Of String))

#End Region

#Region "Events"
	Public Event MessageEvent(ByVal strMessage As String)
	Public Event ErrorEvent(ByVal strErrorMessage As String)
	Public Event WarningEvent(ByVal strWarningMessage As String)
#End Region

#Region "Properties"
	Public ReadOnly Property PeptideHitResultType As clsPHRPReader.ePeptideHitResultType
		Get
			Return mPeptideHitResultType
		End Get
	End Property
#End Region

	''' <summary>
	''' Constructor
	''' </summary>
	''' <param name="strDatasetName">Dataset Name</param>
	''' <param name="strInputFilePath">Input file path</param>
	''' <remarks></remarks>
	Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType)

		mDatasetName = strDatasetName
		mInputFilePath = strInputFilePath
		mPeptideHitResultType = ePeptideHitResultType

		Dim fiFileInfo As System.IO.FileInfo = New System.IO.FileInfo(strInputFilePath)
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

		' Read the ModSummary file (if it exists)
		mModInfoLoaded = LoadModSummary()

		' Read the ResultToSeqMapInfo (if the files exist)
		LoadSeqInfo()

		' The following will be overridden by a derived form of this class
		DefineColumnHeaders()

	End Sub

#Region "Functions overridden by derived classes"
	Protected MustOverride Sub DefineColumnHeaders()
	Public MustOverride Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, ByRef objPSM As clsPSM) As Boolean
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

	Protected Function ConvertModsToNumericMods(ByVal strCleanSequence As String, ByRef lstModifiedResidues As List(Of PeptideHitResultsProcessor.clsAminoAcidModInfo)) As String
		Static sbNewPeptide As New System.Text.StringBuilder

		sbNewPeptide.Length = 0

		For intIndex = 0 To strCleanSequence.Length - 1
			sbNewPeptide.Append(strCleanSequence.Chars(0))

			For Each objModInfo In lstModifiedResidues
				If objModInfo.ResidueLocInPeptide = intIndex + 1 Then
					sbNewPeptide.Append(objModInfo.ModDefinition.ModificationMass.ToString("+0.00000;-0.00000"))
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
				ReportError("ModSummaryFile path is empty; unable to continue")
				Return False
			End If

			If Not System.IO.File.Exists(strModSummaryFilePath) Then
				ReportError("ModSummaryFile not found: " & strModSummaryFilePath)
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

			mResultIDToProteins.Clear()

			If blnSuccess Then
				' Populate mResultIDToProteins

				Dim intResultID As Integer
				Dim intSeqID As Integer

				Dim lstProteinsForSeqID As System.Collections.Generic.List(Of clsProteinInfo) = Nothing
				Dim lstProteinsForResultID As System.Collections.Generic.List(Of String) = Nothing

				For Each objItem As System.Collections.Generic.KeyValuePair(Of Integer, Integer) In mResultToSeqMap

					intResultID = objItem.Key
					intSeqID = objItem.Value
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

	Protected Sub ReportError(ByVal strErrorMessage As String)
		mErrorMessage = strErrorMessage
		RaiseEvent ErrorEvent(strErrorMessage)
	End Sub

	Protected Sub ReportWarning(ByVal strWarningMessage As String)
		RaiseEvent WarningEvent(strWarningMessage)
	End Sub

	Protected Sub ShowMessage(ByVal strMessage As String)
		RaiseEvent MessageEvent(strMessage)
	End Sub

	''' <summary>
	''' Updates the theoretical (computed) monoisotopic mass of objPSM using mResultToSeqMap and mSeqInfo
	''' Also updates the modification info
	''' </summary>
	''' <param name="objPSM"></param>
	''' <returns>True if success, False if objPSM.ResultID is not found in mResultToSeqMap</returns>
	''' <remarks></remarks>
	Protected Function UpdatePSMUsingSeqInfo(ByRef objPSM As clsPSM) As Boolean
		Dim intSeqID As Integer
		Dim objSeqInfo As clsSeqInfo = Nothing

		Dim strMods() As String
		Dim strModDetails() As String

		Dim strMassCorrectionTag As String
		Dim intResidueLoc As Integer
		Dim intPeptideResidueCount As Integer

		Dim eResidueTerminusState As eResidueTerminusStateConstants

		Dim objMatchedModDef As PeptideHitResultsProcessor.clsModificationDefinition = Nothing
		Dim blnMatchFound As Boolean

		Dim blnFavorTerminalMods As Boolean
		Dim blnSuccess As Boolean

		intPeptideResidueCount = objPSM.PeptideCleanSequence.Length
		blnSuccess = False

		' First determine the modified residues present in this peptide
		If mResultToSeqMap.TryGetValue(objPSM.ResultID, intSeqID) Then
			If mSeqInfo.TryGetValue(intSeqID, objSeqInfo) Then
				objPSM.PeptideMonoisotopicMass = objSeqInfo.MonoisotopicMass

				objPSM.ClearModifiedResidues()

				If objSeqInfo.ModCount > 0 Then
					' Split objSeqInfo.ModDescription on the comma character
					strMods = objSeqInfo.ModDescription.Split(","c)

					If Not strMods Is Nothing AndAlso strMods.Count > 0 Then
						For intModIndex As Integer = 0 To strMods.Count - 1

							' Split strMods on the colon characters
							strModDetails = strMods(intModIndex).Split(":"c)

							If Not strModDetails Is Nothing AndAlso strModDetails.Count = 2 Then
								strMassCorrectionTag = strModDetails(0)
								If Integer.TryParse(strModDetails(1), intResidueLoc) Then
									' Find the modification definition in mModInfo
									' Note that a given mass correction tag might be present multiple times in mModInfo, since it could be used as both a static peptide mod and a static peptide terminal mod
									' Thus, if intResidueLoc = 1 or intResidueLoc = objPSM.PeptideCleanSequence.Length then we'll first look for a peptide or protein terminal static mod

									If intResidueLoc = 1 Then
										eResidueTerminusState = eResidueTerminusStateConstants.PeptideNTerminus
										blnFavorTerminalMods = True
									ElseIf intResidueLoc = intPeptideResidueCount Then
										eResidueTerminusState = eResidueTerminusStateConstants.PeptideCTerminus
										blnFavorTerminalMods = True
									Else
										eResidueTerminusState = eResidueTerminusStateConstants.None
										blnFavorTerminalMods = False
									End If

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

		If blnSuccess Then
			objPSM.PeptideWithNumericMods = ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues)
		End If

		Return blnSuccess
	End Function

	Protected Function UpdatePSMFindMatchingModInfo(ByVal strMassCorrectionTag As String, ByVal blnFavorTerminalMods As Boolean, ByVal eResidueTerminusState As eResidueTerminusStateConstants, ByRef objMatchedModDef As PeptideHitResultsProcessor.clsModificationDefinition) As Boolean

		Dim blnMatchFound As Boolean

		blnMatchFound = False
		Do
			For Each objMod In mModInfo
				If strMassCorrectionTag = objMod.MassCorrectionTag Then
					If blnFavorTerminalMods Then
						If objMod.ModificationType = PeptideHitResultsProcessor.clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse _
						   objMod.ModificationType = PeptideHitResultsProcessor.clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod Then

							If eResidueTerminusState = eResidueTerminusStateConstants.PeptideNTerminus AndAlso _
							  (objMod.TargetResiduesContain(N_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(N_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
								blnMatchFound = True
							End If

							If eResidueTerminusState = eResidueTerminusStateConstants.PeptideCTerminus AndAlso _
							  (objMod.TargetResiduesContain(C_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(C_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
								blnMatchFound = True
							End If

						End If
					Else
						blnMatchFound = True
					End If

					If blnFavorTerminalMods Then
						' Match found
						objMatchedModDef = objMod
						blnMatchFound = True
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
