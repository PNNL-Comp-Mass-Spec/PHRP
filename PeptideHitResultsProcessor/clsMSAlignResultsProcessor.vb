﻿Option Strict On

' This class reads in an MSalign results file (txt format) and creates 
' a tab-delimited text file with the data. 
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 11/28/2012
'
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports PHRPReader

Public Class clsMSAlignResultsProcessor
	Inherits clsPHRPBaseClass

	Public Sub New()
		MyBase.New()
		MyBase.mFileDate = "November 28, 2012"
		InitializeLocalVariables()
	End Sub


#Region "Constants and Enums"

	Public Const FILENAME_SUFFIX_MSALIGN_FILE As String = "_MSAlign_ResultTable"

	Public Const N_TERMINUS_SYMBOL_MSALIGN As String = "."
	Public Const C_TERMINUS_SYMBOL_MSALIGN As String = "."

	Private Const UNKNOWN_MSALIGN_MOD_SYMBOL As Char = "?"c

	Public Const DEFAULT_SYN_FILE_PVALUE_THRESHOLD As Single = 0.95

	Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

	Private Const MSALIGN_MOD_MASS_REGEX As String = "\[([+-]*[0-9\.]+)\]"

	Private Const REGEX_OPTIONS As Text.RegularExpressions.RegexOptions = Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.Singleline Or Text.RegularExpressions.RegexOptions.IgnoreCase

	' These columns correspond to the tab-delimited file created directly by MSAlign
	Protected Const MSAlignResultsFileColCount As Integer = 20
	Public Enum eMSAlignResultsFileColumns As Integer
		Data_file_name = 0
		Prsm_ID = 1
		Spectrum_ID = 2
		Scans = 3
		Peaks = 4
		Charge = 5
		Precursor_mass = 6				' Monoisotopic mass value of the observed precursor_mz
		Adjusted_precursor_mass = 7		' Theoretical monoisotopic mass of the peptide (including mods)
		Protein_ID = 8
		Protein_name = 9				' Protein name and description
		Protein_mass = 10
		First_residue = 11
		Last_residue = 12
		Peptide = 13
		Unexpected_modifications = 14
		Matched_peaks = 15
		Matched_fragment_ions = 16
		Pvalue = 17
		Evalue = 18
		FDR = 19
	End Enum

	' These columns correspond to the Synopsis file created by this class
	Protected Const MSAlignSynFileColCount As Integer = 20
	Public Enum eMSAlignSynFileColumns As Integer
		ResultID = 0
		Scan = 1
		Prsm_ID = 2
		Spectrum_ID = 3
		Charge = 4
		PrecursorMZ = 5
		DelM = 6							' Precursor error, in Da
		DelMPPM = 7							' Precursor error, in ppm
		MH = 8								' Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
		Peptide = 9							' This is the sequence with prefix and suffix residues and also with modification symbols
		Protein = 10						' Protein Name (remove description)
		Protein_Mass = 11
		Unexpected_Mod_Count = 12
		Peak_Count = 13
		Matched_Peak_Count = 14
		Matched_Fragment_Ion_Count = 15
		PValue = 16
		Rank_PValue = 17
		EValue = 18
		FDR = 19
	End Enum

#End Region


#Region "Structures"
	Protected Structure udtMSAlignSearchResultType

		Public Data_file_name As String
		Public Scans As String
		Public ScanNum As Integer
		Public Prsm_ID As String
		Public Spectrum_ID As String
		Public Peaks As String
		Public Charge As String
		Public ChargeNum As Short
		Public Precursor_mass As String				' Monoisotopic mass value of the observed precursor_mz
		Public PrecursorMZ As String				' Computed from Precursor_mass
		Public Adjusted_precursor_mass As String	' Theoretical monoisotopic mass of the peptide (including mods), as computed by MSAlign
		Public MH As String							' Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
		Public DelM As String						' Computed using Precursor_mass - Adjusted_precursor_mass
		Public DelM_PPM As String					' Computed using DelM and Adjusted_precursor_mass
		Public Protein_ID As String
		Public Protein As String
		Public Protein_mass As String
		Public First_residue As String
		Public Last_residue As String
		Public Peptide As String
		Public Unexpected_modifications As String
		Public Matched_peaks As String
		Public Matched_fragment_ions As String
		Public Pvalue As String
		Public PValueNum As Double
		Public Evalue As String
		Public FDR As String
		Public RankPValue As Integer

		Public Sub Clear()
			Data_file_name = String.Empty
			Prsm_ID = String.Empty
			Spectrum_ID = String.Empty
			Scans = String.Empty
			ScanNum = 0
			Peaks = String.Empty
			Charge = String.Empty
			ChargeNum = 0
			Precursor_mass = String.Empty
			Adjusted_precursor_mass = String.Empty
			MH = String.Empty
			DelM = String.Empty
			DelM_PPM = String.Empty
			Protein_ID = String.Empty
			Protein = String.Empty
			Protein_mass = String.Empty
			First_residue = String.Empty
			Last_residue = String.Empty
			Peptide = String.Empty
			Unexpected_modifications = String.Empty
			Matched_peaks = String.Empty
			Matched_fragment_ions = String.Empty
			Pvalue = String.Empty
			PValueNum = 0
			Evalue = String.Empty
			FDR = String.Empty
			RankPValue = 0
		End Sub
	End Structure

#End Region

#Region "Classwide Variables"

#End Region

	Private Sub AddCurrentRecordToSearchResults(ByRef intCurrentScanResultsCount As Integer, _
	  ByRef udtSearchResultsCurrentScan() As udtMSAlignSearchResultType, _
	  ByRef lstSearchResults As System.Collections.Generic.List(Of udtMSAlignSearchResultType), _
	  ByRef strErrorLog As String)

		For Each udtSearchResult As udtMSAlignSearchResultType In lstSearchResults

			If intCurrentScanResultsCount >= udtSearchResultsCurrentScan.Length Then
				ReDim Preserve udtSearchResultsCurrentScan(udtSearchResultsCurrentScan.Length * 2 - 1)
			End If

			udtSearchResultsCurrentScan(intCurrentScanResultsCount) = udtSearchResult
			intCurrentScanResultsCount += 1

		Next

	End Sub

	''' <summary>
	''' Step through .PeptideSequenceWithMods
	''' For each residue, check if a static mod is defined that affects that residue
	''' For each mod mass, determine the modification and add to objSearchResult
	''' </summary>
	''' <param name="objSearchResult"></param>
	''' <param name="blnUpdateModOccurrenceCounts"></param>
	''' <remarks></remarks>
	Private Sub AddDynamicAndStaticResidueMods(ByRef objSearchResult As clsSearchResultsMSAlign, ByVal blnUpdateModOccurrenceCounts As Boolean)
		Const NO_RESIDUE As Char = "-"c

		Dim intIndex As Integer, intModIndex As Integer
		Dim chChar As Char
		Dim objModificationDefinition As clsModificationDefinition

		Dim strSequence As String

		Dim blnParsingModMass As Boolean
		Dim strModMassDigits As String = String.Empty

		Dim chMostRecentResidue As Char
		Dim intResidueLocInPeptide As Integer

		Dim chAmbiguousResidue As Char
		Dim intAmbiguousResidueLocInPeptide As Integer

		Dim blnClearAmbiguousResidue As Boolean
		Dim blnStoreAmbiguousResidue As Boolean

		Dim blnSuccess As Boolean

		blnParsingModMass = False
		strModMassDigits = String.Empty

		chMostRecentResidue = NO_RESIDUE
		intResidueLocInPeptide = 0

		chAmbiguousResidue = NO_RESIDUE
		intAmbiguousResidueLocInPeptide = 0

		blnClearAmbiguousResidue = False
		blnStoreAmbiguousResidue = False


		With objSearchResult
			strSequence = .PeptideSequenceWithMods
			For intIndex = 0 To strSequence.Length - 1
				chChar = strSequence.Chars(intIndex)

				If Char.IsLetter(chChar) Then
					chMostRecentResidue = chChar
					intResidueLocInPeptide += 1

					If blnStoreAmbiguousResidue Then
						chAmbiguousResidue = chChar
						intAmbiguousResidueLocInPeptide = intResidueLocInPeptide
						blnStoreAmbiguousResidue = False

					ElseIf blnClearAmbiguousResidue Then
						chAmbiguousResidue = NO_RESIDUE
						blnClearAmbiguousResidue = False

					End If

					For intModIndex = 0 To mPeptideMods.ModificationCount - 1
						If mPeptideMods.GetModificationTypeByIndex(intModIndex) = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
							objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)

							If objModificationDefinition.TargetResiduesContain(chChar) Then
								' Match found; add this modification
								.SearchResultAddModification(objModificationDefinition, chChar, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
							End If
						End If
					Next intModIndex

				ElseIf chChar = "("c Then
					' Start of a mod group
					blnStoreAmbiguousResidue = True

				ElseIf chChar = ")"c Then
					' End of a mod group
					blnClearAmbiguousResidue = True

				ElseIf chChar = "["c Then
					' Mod Mass Start
					strModMassDigits = String.Empty
					blnParsingModMass = True

				ElseIf chChar = "]"c Then
					' Mod Mass End

					If blnParsingModMass Then
						Dim chResidueForMod As Char
						Dim intResidueLocForMod As Integer
						Dim dblModMass As Double

						If chAmbiguousResidue = NO_RESIDUE Then
							chResidueForMod = chMostRecentResidue
							intResidueLocForMod = intResidueLocInPeptide
						Else
							' Ambigous mod
							' We'll associate it with the first residue of the mod group
							chResidueForMod = chAmbiguousResidue
							intResidueLocForMod = intAmbiguousResidueLocInPeptide
						End If

						If Double.TryParse(strModMassDigits, dblModMass) Then
							If intResidueLocForMod = 0 Then
								' Modification is at the peptide N-terminus
								intResidueLocForMod = 1
							End If

							blnSuccess = .SearchResultAddModification(dblModMass, chResidueForMod, intResidueLocForMod, .DetermineResidueTerminusState(intResidueLocForMod), blnUpdateModOccurrenceCounts)

							If Not blnSuccess Then
								Dim strErrorMessage As String = .ErrorMessage
								If String.IsNullOrEmpty(strErrorMessage) Then
									strErrorMessage = "SearchResultAddDynamicModification returned false for mod mass " & strModMassDigits
								End If
								SetErrorMessage(strErrorMessage & "; ResultID = " & .ResultID)

							End If
						End If

						blnParsingModMass = False
					End If

				ElseIf blnParsingModMass Then
					strModMassDigits &= chChar

				Else
					' Unrecognized symbol; ignore it

				End If

			Next intIndex
		End With
	End Sub

	Private Function AddModificationsAndComputeMass(ByRef objSearchResult As clsSearchResultsMSAlign, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
		Const ALLOW_DUPLICATE_MOD_ON_TERMINUS As Boolean = True

		Dim blnSuccess As Boolean

		Try
			' If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
			objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

			' Parse .PeptideSequenceWithMods to determine the modified residues present
			AddDynamicAndStaticResidueMods(objSearchResult, blnUpdateModOccurrenceCounts)

			' Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
			' Since Inspect allows a terminal peptide residue to be modified twice, we'll allow that to happen,
			'  even though, biologically, that's typically not possible
			' However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus 
			'  (where two COOH groups are present)
			objSearchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts)

			' Compute the monoisotopic mass for this peptide
			objSearchResult.ComputeMonoisotopicMass()

			' Populate .PeptideModDescription
			objSearchResult.UpdateModDescription()

			blnSuccess = True
		Catch ex As Exception
			blnSuccess = False
		End Try

		Return blnSuccess

	End Function

	''' <summary>
	''' Ranks each entry (calling procedure should have already sorted the data by Scan, Charge, and SpecProb)
	''' </summary>
	''' <param name="udtSearchResults"></param>
	''' <param name="intStartIndex"></param>
	''' <param name="intEndIndex"></param>
	''' <remarks></remarks>
	Private Sub AssignRankAndDeltaNormValues(ByRef udtSearchResults() As udtMSAlignSearchResultType, _
	  ByVal intStartIndex As Integer, _
	  ByVal intEndIndex As Integer)

		Dim intIndex As Integer

		Dim intLastCharge As Integer
		Dim dblLastValue As Double

		Dim intCurrentRank As Integer

		For intIndex = intStartIndex To intEndIndex
			If intIndex = intStartIndex OrElse udtSearchResults(intIndex).ChargeNum <> intLastCharge Then
				intLastCharge = udtSearchResults(intIndex).ChargeNum
				dblLastValue = udtSearchResults(intIndex).PValueNum
				intCurrentRank = 1
			Else
				If udtSearchResults(intIndex).PValueNum <> dblLastValue Then
					dblLastValue = udtSearchResults(intIndex).PValueNum
					intCurrentRank += 1
				End If
			End If

			udtSearchResults(intIndex).RankPValue = intCurrentRank
		Next intIndex

	End Sub

	Protected Function ComputePeptideMass(ByVal strPeptide As String, ByVal dblTotalModMass As Double) As Double

		Dim strCleanSequence As String
		Dim dblMass As Double

		strCleanSequence = GetCleanSequence(strPeptide)

		Dim objMassCalculator As clsPeptideMassCalculator = New clsPeptideMassCalculator()

		dblMass = objMassCalculator.ComputeSequenceMass(strCleanSequence)

		If dblTotalModMass <> 0 Then
			dblMass += dblTotalModMass
		End If

		Return dblMass

	End Function

	''' <summary>
	''' Computes the total of all modification masses defined for the peptide
	''' </summary>
	''' <param name="strPeptide">Peptide sequence, with mod masses in the form [23.5432]</param>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Function ComputeTotalModMass(ByVal strPeptide As String) As Double

		Static reModMassRegEx As New System.Text.RegularExpressions.Regex(MSALIGN_MOD_MASS_REGEX, REGEX_OPTIONS)

		Dim dblTotalModMass As Double
		Dim dblModMassFound As Double

		dblTotalModMass = 0

		For Each reMatch As System.Text.RegularExpressions.Match In reModMassRegEx.Matches(strPeptide)
			If Double.TryParse(reMatch.Groups(1).Value, dblModMassFound) Then
				dblTotalModMass += dblModMassFound
			End If
		Next

		Return dblTotalModMass

	End Function

	Protected Overrides Function ConstructPepToProteinMapFilePath(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal MTS As Boolean) As String
		Dim strPepToProteinMapFilePath As String = String.Empty

		strPepToProteinMapFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
		If strPepToProteinMapFilePath.ToLower().EndsWith("_msalign_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_msalign_fht") Then
			' Remove _syn or _fht
			strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4)
		End If

		Return MyBase.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS)

	End Function

	''' <summary>
	''' This routine creates a first hits file or synopsis file from the output from MSAlign
	''' The synopsis file includes every result with a p-value below a set threshold or a SpecProb below a certain threshold
	''' The first-hits file includes the results with the lowest SpecProb (for each scan and charge)
	''' </summary>
	''' <param name="strInputFilePath"></param>
	''' <param name="strOutputFilePath"></param>
	''' <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Protected Function CreateSynResultsFile(ByVal strInputFilePath As String, _
	   ByVal strOutputFilePath As String, _
	   Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean

		Dim strLineIn As String
		Dim protein As String = String.Empty

		Dim lstSearchResults As New System.Collections.Generic.List(Of udtMSAlignSearchResultType)

		Dim intSearchResultsCount As Integer
		Dim udtSearchResults() As udtMSAlignSearchResultType
		Dim sngPercentComplete As Single

		Dim intFilteredSearchResultCount As Integer
		Dim udtFilteredSearchResults() As udtMSAlignSearchResultType

		Dim blnHeaderParsed As Boolean

		Dim intColumnMapping() As Integer = Nothing

		Dim intResultsProcessed As Integer

		Dim blnSuccess As Boolean
		Dim blnValidSearchResult As Boolean

		Dim strErrorLog As String = String.Empty

		Try

			Try
				' Open the input file and parse it
				' Initialize the stream reader and the stream Text writer
				Using srDataFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					Using swResultFile As System.IO.StreamWriter = New System.IO.StreamWriter(New System.IO.FileStream(strOutputFilePath, IO.FileMode.Create, IO.FileAccess.Write, IO.FileShare.Read))

						Dim ioResultFile As System.IO.FileInfo
						ioResultFile = New System.IO.FileInfo(strOutputFilePath)

						strErrorLog = String.Empty
						intResultsProcessed = 0
						blnHeaderParsed = False

						' Initialize array that will hold all of the records in the MSAlign result file
						intSearchResultsCount = 0
						ReDim udtSearchResults(999)

						' Initialize the array that will hold all of the records that will ultimately be written out to disk
						intFilteredSearchResultCount = 0
						ReDim udtFilteredSearchResults(999)

						' Parse the input file
						Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing
							strLineIn = srDataFile.ReadLine()
							If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then

								If Not blnHeaderParsed Then
									blnSuccess = ParseMSAlignResultsFileHeaderLine(strLineIn, intColumnMapping)
									If Not blnSuccess Then
										' Error parsing header
										SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
										Exit Try
									End If
									blnHeaderParsed = True

									' Write the header line
									WriteSynFHTFileHeader(swResultFile, strErrorLog)
								Else

									blnValidSearchResult = ParseMSAlignResultsFileEntry(strLineIn, lstSearchResults, strErrorLog, intResultsProcessed, intColumnMapping)

									If blnValidSearchResult Then
										AddCurrentRecordToSearchResults(intSearchResultsCount, udtSearchResults, lstSearchResults, strErrorLog)
									End If

									' Update the progress
									sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
									If mCreateProteinModsFile Then
										sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
									End If
									UpdateProgress(sngPercentComplete)

									intResultsProcessed += 1
								End If
							End If
						Loop

						' Sort the SearchResults by scan, charge, and ascending SpecProb
						Array.Sort(udtSearchResults, 0, intSearchResultsCount, New MSAlignSearchResultsComparerScanChargePValuePeptide)

						' Now filter the data

						' Initialize variables
						Dim intStartIndex As Integer = 0
						Dim intEndIndex As Integer

						intStartIndex = 0
						intEndIndex = 0
						Do While intStartIndex < intSearchResultsCount
							intEndIndex = intStartIndex
							Do While intEndIndex + 1 < intSearchResultsCount AndAlso udtSearchResults(intEndIndex + 1).ScanNum = udtSearchResults(intStartIndex).ScanNum
								intEndIndex += 1
							Loop

							' Store the results for this scan
							StoreSynMatches(udtSearchResults, intStartIndex, intEndIndex, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)

							intStartIndex = intEndIndex + 1
						Loop

						' Sort the data in udtFilteredSearchResults then write out to disk
						SortAndWriteFilteredSearchResults(swResultFile, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
					End Using
				End Using

				' Inform the user if any errors occurred
				If strErrorLog.Length > 0 Then
					SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
				End If

				blnSuccess = True

			Catch ex As Exception
				SetErrorMessage(ex.Message)
				SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
				blnSuccess = False
			End Try
		Catch ex As Exception
			SetErrorMessage(ex.Message)
			SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
			blnSuccess = False
		End Try

		Return blnSuccess

	End Function

	Protected Function ExtractModInfoFromMSAlignParamFile(ByVal strMSAlignParamFilePath As String, ByRef lstModInfo As Generic.List(Of clsModificationDefinition)) As Boolean

		Dim strLineIn As String
		Dim kvSetting As Generic.KeyValuePair(Of String, String)

		Dim objModDef As clsModificationDefinition

		Dim blnSuccess As Boolean = False

		Try
			' Initialize the modification list
			If lstModInfo Is Nothing Then
				lstModInfo = New Generic.List(Of clsModificationDefinition)
			Else
				lstModInfo.Clear()
			End If

			If String.IsNullOrEmpty(strMSAlignParamFilePath) Then
				SetErrorMessage("MSAlign Parameter File name not defined; unable to extract mod info")
				SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
				Return False
			End If

			If Not System.IO.File.Exists(strMSAlignParamFilePath) Then
				SetErrorMessage("MSAlign param file not found: " & strMSAlignParamFilePath)
			Else
				' Read the contents of the parameter (or mods) file
				Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strMSAlignParamFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					Do While srInFile.Peek <> -1
						strLineIn = srInFile.ReadLine().Trim()

						If strLineIn.Length > 0 Then

							If strLineIn.StartsWith("#"c) Then
								' Comment line; skip it
							Else
								' Split the line on the equals sign
								kvSetting = PHRPReader.clsPHRPParser.ParseKeyValueSetting(strLineIn, "="c, "#")

								If kvSetting.Key.ToLower() = "cysteineProtection".ToLower() Then

									Select Case kvSetting.Value.ToUpper()
										Case "C57"
											objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 57.0215, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcet")
											lstModInfo.Add(objModDef)

										Case "C58"
											objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 58.0055, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcid")
											lstModInfo.Add(objModDef)

									End Select

								End If

							End If

						End If
					Loop
				End Using

				Console.WriteLine()

				blnSuccess = True

			End If
		Catch ex As Exception
			SetErrorMessage("Error reading the MSAlign parameter file (" & System.IO.Path.GetFileName(strMSAlignParamFilePath) & "): " & ex.Message)
			SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
			blnSuccess = False
		End Try

		Return blnSuccess

	End Function

	Private Sub InitializeLocalVariables()
		' Nothing to do at present
	End Sub


	Protected Function ParseMSAlignSynopsisFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByRef lstPepToProteinMapping As Generic.List(Of udtPepToProteinMappingType), ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean) As Boolean
		' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

		Dim strPreviousPValue As String

		' Note that MSAlign synopsis files are normally sorted on PValue value, ascending
		' In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
		'  we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
		' Although this was a possiblity with Inspect, it likely never occurs for MSAlign
		'  But, we'll keep the check in place just in case

		Dim htPeptidesFoundForPValueLevel As Hashtable

		Dim strKey As String

		Dim strLineIn As String
		Dim strModificationSummaryFilePath As String

		Dim objSearchResult As clsSearchResultsMSAlign

		Dim intResultsProcessed As Integer
		Dim intPepToProteinMapIndex As Integer
		Dim sngPercentComplete As Single

		Dim blnHeaderParsed As Boolean
		Dim intColumnMapping() As Integer = Nothing

		Dim strCurrentPeptideWithMods As String = String.Empty
		Dim strCurrentProtein As String

		Dim blnSuccess As Boolean
		Dim blnValidSearchResult As Boolean
		Dim blnFirstMatchForGroup As Boolean

		Dim strErrorLog As String = String.Empty

		Try
			' Possibly reset the mass correction tags and Mod Definitions
			If blnResetMassCorrectionTagsAndModificationDefinitions Then
				ResetMassCorrectionTagsAndModificationDefinitions()
			End If

			' Reset .OccurrenceCount
			mPeptideMods.ResetOccurrenceCountStats()

			' Initialize objSearchResult
			objSearchResult = New clsSearchResultsMSAlign(mPeptideMods)

			' Initialize htPeptidesFoundForPValueLevel
			htPeptidesFoundForPValueLevel = New Hashtable
			strPreviousPValue = String.Empty

			' Assure that lstPepToProteinMapping is sorted on peptide
			If lstPepToProteinMapping.Count > 1 Then
				lstPepToProteinMapping.Sort(New PepToProteinMappingComparer)
			End If

			Try
				objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(mEnzymeMatchSpec, mPeptideNTerminusMassChange, mPeptideCTerminusMassChange)

				' Open the input file and parse it
				' Initialize the stream reader
				Using srDataFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					strErrorLog = String.Empty
					intResultsProcessed = 0
					blnHeaderParsed = False

					' Create the output files
					Dim strBaseOutputFilePath As String = System.IO.Path.Combine(strOutputFolderPath, System.IO.Path.GetFileName(strInputFilePath))
					blnSuccess = MyBase.InitializeSequenceOutputFiles(strBaseOutputFilePath)

					' Parse the input file
					Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing

						strLineIn = srDataFile.ReadLine()
						If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then

							If Not blnHeaderParsed Then
								blnSuccess = ParseMSAlignSynFileHeaderLine(strLineIn, intColumnMapping)
								If Not blnSuccess Then
									' Error parsing header
									SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
									Exit Try
								End If
								blnHeaderParsed = True
							Else

								blnValidSearchResult = ParseMSAlignSynFileEntry(strLineIn, objSearchResult, strErrorLog, _
								  intResultsProcessed, intColumnMapping, _
								  strCurrentPeptideWithMods)

								If blnValidSearchResult Then
									strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge

									If objSearchResult.PValue = strPreviousPValue Then
										' New result has the same PValue as the previous result
										' See if htPeptidesFoundForPValueLevel contains the peptide, scan and charge

										If htPeptidesFoundForPValueLevel.ContainsKey(strKey) Then
											blnFirstMatchForGroup = False
										Else
											htPeptidesFoundForPValueLevel.Add(strKey, 1)
											blnFirstMatchForGroup = True
										End If

									Else
										' New PValue
										' Reset htPeptidesFoundForScan
										htPeptidesFoundForPValueLevel.Clear()

										' Update strPreviousPValue
										strPreviousPValue = objSearchResult.PValue

										' Append a new entry to htPeptidesFoundForScan
										htPeptidesFoundForPValueLevel.Add(strKey, 1)
										blnFirstMatchForGroup = True
									End If

									blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup)
									If Not blnSuccess Then
										If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
											strErrorLog &= "Error adding modifications to sequence at RowIndex '" & objSearchResult.ResultID & "'" & ControlChars.NewLine
										End If
									End If

									MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResult, clsSearchResultsBaseClass), blnFirstMatchForGroup)

									If lstPepToProteinMapping.Count > 0 Then
										' Add the additional proteins for this peptide

										' Use binary search to find this peptide in lstPepToProteinMapping
										intPepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(lstPepToProteinMapping, strCurrentPeptideWithMods)

										If intPepToProteinMapIndex >= 0 Then
											' Call MyBase.SaveResultsFileEntrySeqInfo for each entry in lstPepToProteinMapping() for peptide , skipping objSearchResult.ProteinName
											strCurrentProtein = String.Copy(objSearchResult.ProteinName)
											Do
												If lstPepToProteinMapping(intPepToProteinMapIndex).Protein <> strCurrentProtein Then
													objSearchResult.ProteinName = String.Copy(lstPepToProteinMapping(intPepToProteinMapIndex).Protein)
													MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResult, clsSearchResultsBaseClass), False)
												End If

												intPepToProteinMapIndex += 1
											Loop While intPepToProteinMapIndex < lstPepToProteinMapping.Count AndAlso strCurrentPeptideWithMods = lstPepToProteinMapping(intPepToProteinMapIndex).Peptide
										Else
											' Match not found; this is unexpected
											ReportWarning("no match for '" & strCurrentPeptideWithMods & "' in lstPepToProteinMapping")
										End If
									End If

								End If

								' Update the progress
								sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
								If mCreateProteinModsFile Then
									sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
								End If
								UpdateProgress(sngPercentComplete)

								intResultsProcessed += 1
							End If
						End If

					Loop

				End Using

				If mCreateModificationSummaryFile Then
					' Create the modification summary file
					Dim fiInputFile As System.IO.FileInfo = New System.IO.FileInfo(strInputFilePath)
					strModificationSummaryFilePath = System.IO.Path.GetFileName(MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_MOD_SUMMARY))
					strModificationSummaryFilePath = System.IO.Path.Combine(strOutputFolderPath, strModificationSummaryFilePath)

					SaveModificationSummaryFile(strModificationSummaryFilePath)
				End If

				' Inform the user if any errors occurred
				If strErrorLog.Length > 0 Then
					SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
				End If

				blnSuccess = True

			Catch ex As Exception
				SetErrorMessage(ex.Message)
				SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
				blnSuccess = False
			Finally
				MyBase.CloseSequenceOutputFiles()
			End Try
		Catch ex As Exception
			SetErrorMessage(ex.Message)
			SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
			blnSuccess = False
		End Try

		Return blnSuccess


	End Function

	Private Function ParseMSAlignResultsFileEntry(ByRef strLineIn As String, _
	   ByRef lstSearchResults As System.Collections.Generic.List(Of udtMSAlignSearchResultType), _
	   ByRef strErrorLog As String, _
	   ByVal intResultsProcessed As Integer, _
	   ByRef intColumnMapping() As Integer) As Boolean

		' Parses an entry from the MSAlign results file

		Dim udtSearchResult As udtMSAlignSearchResultType = New udtMSAlignSearchResultType
		Dim strSplitLine() As String = Nothing

		Dim dblPrecursorMonoMass As Double				' Observed m/z, converted to monoisotopic mass
		Dim dblPeptideMonoMassMSAlign As Double		' Theoretical peptide monoisotopic mass, including mods, as computed by MSAlign
		Dim dblPeptideMonoMassPHRP As Double		' Theoretical peptide monoisotopic mass, including mods, as computed by PHRP

		Dim dblPrecursorMZ As Double
		Dim dblDelM As Double

		Dim dblTotalModMass As Double

		Dim blnValidSearchResult As Boolean

		Try
			' Set this to False for now
			blnValidSearchResult = False

			' Reset lstSearchResults
			If lstSearchResults Is Nothing Then
				lstSearchResults = New System.Collections.Generic.List(Of udtMSAlignSearchResultType)
			Else
				lstSearchResults.Clear()
			End If

			udtSearchResult.Clear()
			strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

			If strSplitLine.Length >= 13 Then

				With udtSearchResult
					If Not GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Data_file_name), .Data_file_name) Then
						Throw New EvaluateException("Data_file_name column is missing or invalid")
					End If
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Prsm_ID), .Prsm_ID)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Spectrum_ID), .Spectrum_ID)

					If Not GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Scans), .Scans) Then
						Throw New EvaluateException("Scan(s) column is missing or invalid")
					End If

					If Not Integer.TryParse(.Scans, .ScanNum) Then
						' .Scans likely has a list of scan numbers; extract the first scan number from .scans
						Dim strScanNumberDigits As String = String.Empty
						For Each chChar As Char In .Scans
							If Char.IsDigit(chChar) Then
								strScanNumberDigits &= chChar
							End If
						Next
						If Not Integer.TryParse(strScanNumberDigits, .ScanNum) Then
							ReportWarning("Error parsing out the scan number from the scan list; could not find an integer: " & .Scans)
							.ScanNum = 0
						End If
					End If

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Peaks), .Peaks)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Charge), .Charge)
					.ChargeNum = CShort(CIntSafe(.Charge, 0))

					' Monoisotopic mass value of the observed precursor_mz
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Precursor_mass), .Precursor_mass)

					If Double.TryParse(.Precursor_mass, dblPrecursorMonoMass) Then
						If .ChargeNum > 0 Then
							dblPrecursorMZ = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMonoMass, 0, .ChargeNum)
							.PrecursorMZ = NumToString(dblPrecursorMZ, 6, True)
						End If
					End If

					' Theoretical monoisotopic mass of the peptide (including mods), as computed by MSAlign
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Adjusted_precursor_mass), .Adjusted_precursor_mass)

					' Compute DelM and DelM_PPM
					If Double.TryParse(.Adjusted_precursor_mass, dblPeptideMonoMassMSAlign) Then
						dblDelM = dblPrecursorMonoMass - dblPeptideMonoMassMSAlign
						.DelM = NumToString(dblDelM, 6, True)

						If dblPrecursorMZ > 0 Then
							.DelM_PPM = NumToString(clsPeptideMassCalculator.MassToPPM(dblDelM, dblPrecursorMZ), 5, True)
						Else
							.DelM_PPM = NumToString(clsPeptideMassCalculator.MassToPPM(dblDelM, 1000), 5, True)
						End If

					End If

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Protein_ID), .Protein_ID)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Protein_name), .Protein)
					.Protein = TruncateProteinName(.Protein)

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Protein_mass), .Protein_mass)

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.First_residue), .First_residue)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Last_residue), .Last_residue)

					If Not GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Peptide), .Peptide) Then
						Throw New EvaluateException("Peptide column is missing or invalid")
					End If

					' Add the standard terminus symbols to the peptide sequence
					.Peptide = ReplaceTerminus(.Peptide)

					' Parse the sequence to determine the total mod mass
					' Note that we do not remove any of the mod symbols since MSAlign identifies mods by mass alone, and since mods can ambiguously apply to residues
					dblTotalModMass = ComputeTotalModMass(.Peptide)

					' Compute monoisotopic mass of the peptide
					dblPeptideMonoMassPHRP = ComputePeptideMass(.Peptide, dblTotalModMass)

					If Math.Abs(dblPeptideMonoMassPHRP - dblPeptideMonoMassMSAlign) > 0.1 Then
						' Computed monoisotopic mass values differ by more than 0.1 Da; this is unexpected
						Dim strFirst30Residues As String
						If .Peptide.Length < 27 Then
							strFirst30Residues = .Peptide
						Else
							strFirst30Residues = .Peptide.Substring(0, 27) & "..."
						End If
						ReportWarning("The monoisotopic mass computed by PHRP is more than 0.1 Da away from the mass computed by MSAlign: " & dblPeptideMonoMassPHRP.ToString("0.0000") & " vs. " & dblPeptideMonoMassMSAlign.ToString("0.0000") & "; peptide " & strFirst30Residues)
					End If

					' Store the monoisotopic MH value in .MH; note that this is (M+H)+
					.MH = NumToString(clsPeptideMassCalculator.ConvoluteMass(dblPeptideMonoMassPHRP, 0, 1), 6, True)

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Unexpected_modifications), .Unexpected_modifications)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Matched_peaks), .Matched_peaks)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Matched_fragment_ions), .Matched_fragment_ions)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Pvalue), .Pvalue)
					If Not Double.TryParse(.Pvalue, .PValueNum) Then .PValueNum = 0

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.Evalue), .Evalue)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignResultsFileColumns.FDR), .FDR)

				End With

				' Append udtSearchResult to lstSearchResults
				lstSearchResults.Add(udtSearchResult)

				blnValidSearchResult = True
			End If

		Catch ex As Exception
			' Error parsing this row from the MSAlign results file
			If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
				If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
					strErrorLog &= "Error parsing MSAlign Results in ParseMSAlignResultsFileEntry for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
				Else
					strErrorLog &= "Error parsing MSAlign Results in ParseMSAlignResultsFileEntry" & ControlChars.NewLine
				End If
			End If
			blnValidSearchResult = False
		End Try

		Return blnValidSearchResult

	End Function

	Private Function ParseMSAlignResultsFileHeaderLine(ByVal strLineIn As String, ByRef intColumnMapping() As Integer) As Boolean

		' Parse the header line

		' The expected header from MSAlign is:
		' Data_file_name    Prsm_ID    Spectrum_ID    Scan(s)    #peaks    Charge    Precursor_mass    Adjusted_precursor_mass    Protein_ID    Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions    P-value    E-value    FDR


		Dim strSplitLine() As String
		Dim eResultFileColumn As eMSAlignResultsFileColumns
		Dim lstColumnNames As System.Collections.Generic.SortedDictionary(Of String, eMSAlignResultsFileColumns)
		lstColumnNames = New System.Collections.Generic.SortedDictionary(Of String, eMSAlignResultsFileColumns)(StringComparer.CurrentCultureIgnoreCase)

		ReDim intColumnMapping(MSAlignResultsFileColCount - 1)

		lstColumnNames.Add("Data_file_name", eMSAlignResultsFileColumns.Data_file_name)
		lstColumnNames.Add("Prsm_ID", eMSAlignResultsFileColumns.Prsm_ID)
		lstColumnNames.Add("Spectrum_ID", eMSAlignResultsFileColumns.Spectrum_ID)
		lstColumnNames.Add("Scan(s)", eMSAlignResultsFileColumns.Scans)
		lstColumnNames.Add("#peaks", eMSAlignResultsFileColumns.Peaks)
		lstColumnNames.Add("Charge", eMSAlignResultsFileColumns.Charge)
		lstColumnNames.Add("Precursor_mass", eMSAlignResultsFileColumns.Precursor_mass)
		lstColumnNames.Add("Adjusted_precursor_mass", eMSAlignResultsFileColumns.Adjusted_precursor_mass)
		lstColumnNames.Add("Protein_ID", eMSAlignResultsFileColumns.Protein_ID)
		lstColumnNames.Add("Protein_name", eMSAlignResultsFileColumns.Protein_name)
		lstColumnNames.Add("Protein_mass", eMSAlignResultsFileColumns.Protein_mass)
		lstColumnNames.Add("First_residue", eMSAlignResultsFileColumns.First_residue)
		lstColumnNames.Add("Last_residue", eMSAlignResultsFileColumns.Last_residue)
		lstColumnNames.Add("Peptide", eMSAlignResultsFileColumns.Peptide)
		lstColumnNames.Add("#unexpected_modifications", eMSAlignResultsFileColumns.Unexpected_modifications)
		lstColumnNames.Add("#matched_peaks", eMSAlignResultsFileColumns.Matched_peaks)
		lstColumnNames.Add("#matched_fragment_ions", eMSAlignResultsFileColumns.Matched_fragment_ions)
		lstColumnNames.Add("P-value", eMSAlignResultsFileColumns.Pvalue)
		lstColumnNames.Add("E-value", eMSAlignResultsFileColumns.Evalue)
		lstColumnNames.Add("FDR", eMSAlignResultsFileColumns.FDR)

		Try
			' Initialize each entry in intColumnMapping to -1
			For intIndex As Integer = 0 To intColumnMapping.Length - 1
				intColumnMapping(intIndex) = -1
			Next

			strSplitLine = strLineIn.Split(ControlChars.Tab)
			For intIndex As Integer = 0 To strSplitLine.Length - 1
				If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
					' Recognized column name; update intColumnMapping
					intColumnMapping(eResultFileColumn) = intIndex
				Else
					' Unrecognized column name
					Console.WriteLine("Warning: Unrecognized column header name '" & strSplitLine(intIndex) & "' in ParseMSAlignResultsFileHeaderLine")
				End If
			Next

		Catch ex As Exception
			SetErrorMessage("Error parsing header in MSAlign results file: " & ex.Message)
			Return False
		End Try

		Return True

	End Function

	Private Function ParseMSAlignSynFileHeaderLine(ByVal strLineIn As String, ByRef intColumnMapping() As Integer) As Boolean

		' Parse the header line

		Dim strSplitLine() As String
		Dim eResultFileColumn As eMSAlignSynFileColumns
		Dim lstColumnNames As System.Collections.Generic.SortedDictionary(Of String, eMSAlignSynFileColumns)
		lstColumnNames = New System.Collections.Generic.SortedDictionary(Of String, eMSAlignSynFileColumns)(StringComparer.CurrentCultureIgnoreCase)

		ReDim intColumnMapping(MSAlignSynFileColCount - 1)

		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_ResultID, eMSAlignSynFileColumns.ResultID)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Scan, eMSAlignSynFileColumns.Scan)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Prsm_ID, eMSAlignSynFileColumns.Prsm_ID)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Spectrum_ID, eMSAlignSynFileColumns.Spectrum_ID)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Charge, eMSAlignSynFileColumns.Charge)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_PrecursorMZ, eMSAlignSynFileColumns.PrecursorMZ)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_DelM, eMSAlignSynFileColumns.DelM)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_DelM_PPM, eMSAlignSynFileColumns.DelMPPM)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_MH, eMSAlignSynFileColumns.MH)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Peptide, eMSAlignSynFileColumns.Peptide)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Protein, eMSAlignSynFileColumns.Protein)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Protein_Mass, eMSAlignSynFileColumns.Protein_Mass)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Unexpected_Mod_Count, eMSAlignSynFileColumns.Unexpected_Mod_Count)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Peak_Count, eMSAlignSynFileColumns.Peak_Count)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Matched_Peak_Count, eMSAlignSynFileColumns.Matched_Peak_Count)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Matched_Fragment_Ion_Count, eMSAlignSynFileColumns.Matched_Fragment_Ion_Count)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_PValue, eMSAlignSynFileColumns.PValue)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Rank_PValue, eMSAlignSynFileColumns.Rank_PValue)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_EValue, eMSAlignSynFileColumns.EValue)
		lstColumnNames.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_FDR, eMSAlignSynFileColumns.FDR)

		Try
			' Initialize each entry in intColumnMapping to -1
			For intIndex As Integer = 0 To intColumnMapping.Length - 1
				intColumnMapping(intIndex) = -1
			Next

			strSplitLine = strLineIn.Split(ControlChars.Tab)
			For intIndex As Integer = 0 To strSplitLine.Length - 1
				If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
					' Recognized column name; update intColumnMapping
					intColumnMapping(eResultFileColumn) = intIndex
				End If

			Next

		Catch ex As Exception
			SetErrorMessage("Error parsing header in MSAlign synopsis file: " & ex.Message)
			Return False
		End Try

		Return True

	End Function

	Private Function ParseMSAlignSynFileEntry(ByRef strLineIn As String, _
	  ByRef objSearchResult As clsSearchResultsMSAlign, _
	  ByRef strErrorLog As String, _
	  ByVal intResultsProcessed As Integer, _
	  ByRef intColumnMapping() As Integer, _
	  ByRef strPeptideSequenceWithMods As String) As Boolean

		' Parses an entry from the MSAlign Synopsis file

		Dim strSplitLine() As String = Nothing
		Dim strValue As String = String.Empty

		Dim blnValidSearchResult As Boolean

		Try
			' Set this to False for now
			blnValidSearchResult = False

			' Reset objSearchResult
			objSearchResult.Clear()
			strPeptideSequenceWithMods = String.Empty

			strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

			If strSplitLine.Length >= 15 Then

				With objSearchResult
					If Not GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.ResultID), strValue) Then
						If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
							strErrorLog &= "Error reading ResultID value from MSAlign Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
						End If
						Exit Try
					End If

					.ResultID = Integer.Parse(strValue)

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Scan), .Scan)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Charge), .Charge)

					If Not GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Peptide), strPeptideSequenceWithMods) Then
						If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
							strErrorLog &= "Error reading Peptide sequence value from MSAlign Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
						End If
						Exit Try
					End If

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Protein), .ProteinName)
					.MultipleProteinCount = "0"

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.DelM), .MSAlignComputedDelM)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.DelMPPM), .MSAlignComputedDelMPPM)

					.PeptideDeltaMass = .MSAlignComputedDelM

					' Note: .PeptideDeltaMass is stored in the MSAlign results file as "Observed_Mass - Theoretical_Mass"
					' However, in MTS .peptideDeltaMass is "Theoretical - Observed"
					' Therefore, we will negate .peptideDeltaMass
					Try
						.PeptideDeltaMass = (-Double.Parse(.PeptideDeltaMass)).ToString
					Catch ex As Exception
						' Error; Leave .peptideDeltaMass unchanged
					End Try

					' Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
					.SetPeptideSequenceWithMods(strPeptideSequenceWithMods, True, True)

				End With

				Dim objSearchResultBase As clsSearchResultsBaseClass
				objSearchResultBase = DirectCast(objSearchResult, clsSearchResultsBaseClass)

				MyBase.ComputePseudoPeptideLocInProtein(objSearchResultBase)

				With objSearchResult

					' Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
					' If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide 
					' will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
					.ComputePeptideCleavageStateInProtein()

					' Read the remaining data values
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Prsm_ID), .Prsm_ID)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Spectrum_ID), .Spectrum_ID)

					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.PrecursorMZ), .Precursor_mz)
					
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.MH), .ParentIonMH)
					
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Protein_Mass), .Protein_Mass)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Unexpected_Mod_Count), .Unexpected_Mod_Count)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Peak_Count), .Peak_Count)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Matched_Peak_Count), .Matched_Peak_Count)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Matched_Fragment_Ion_Count), .Matched_Fragment_Ion_Count)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.PValue), .PValue)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.Rank_PValue), .Rank_PValue)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.EValue), .EValue)
					GetColumnValue(strSplitLine, intColumnMapping(eMSAlignSynFileColumns.FDR), .FDR)

				End With

				blnValidSearchResult = True
			End If

		Catch ex As Exception
			' Error parsing this row from the synopsis or first hits file
			If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
				If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
					strErrorLog &= "Error parsing MSAlign Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
				Else
					strErrorLog &= "Error parsing MSAlign Results in ParseMSAlignSynFileEntry" & ControlChars.NewLine
				End If
			End If
			blnValidSearchResult = False
		End Try

		Return blnValidSearchResult

	End Function

	''' <summary>
	''' Main processing function
	''' </summary>
	''' <param name="strInputFilePath">MSAlign results file</param>
	''' <param name="strOutputFolderPath">Output folder</param>
	''' <param name="strParameterFilePath">Parameter file</param>
	''' <returns>True if success, False if failure</returns>
	Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean

		Dim ioInputFile As System.IO.FileInfo

		Dim strBaseName As String = String.Empty
		Dim strSynOutputFilePath As String = String.Empty

		Dim lstMSAlignModInfo As Generic.List(Of clsModificationDefinition)
		Dim lstPepToProteinMapping As Generic.List(Of udtPepToProteinMappingType)
		Dim strMTSPepToProteinMapFilePath As String = String.Empty


		Dim blnSuccess As Boolean

		If Not LoadParameterFileSettings(strParameterFilePath) Then
			SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, True)
			Return False
		End If

		Try
			If String.IsNullOrWhiteSpace(strInputFilePath) Then
				SetErrorMessage("Input file name is empty")
				SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath)
			Else

				blnSuccess = ResetMassCorrectionTagsAndModificationDefinitions()
				If Not blnSuccess Then
					Exit Try
				End If

				MyBase.ResetProgress("Parsing " & System.IO.Path.GetFileName(strInputFilePath))

				If CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
					Try
						' Obtain the full path to the input file
						ioInputFile = New System.IO.FileInfo(strInputFilePath)

						lstMSAlignModInfo = New Generic.List(Of clsModificationDefinition)
						lstPepToProteinMapping = New Generic.List(Of udtPepToProteinMappingType)

						' Load the MSAlign Parameter File so that we can determine whether Cysteine residues are statically modified
						ExtractModInfoFromMSAlignParamFile(mSearchToolParameterFilePath, lstMSAlignModInfo)

						' Resolve the mods in lstMSAlignModInfo with the ModDefs mods
						ResolveMSAlignModsWithModDefinitions(lstMSAlignModInfo)

						' Define the base output filename using strInputFilePath
						strBaseName = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)

						' Auto-replace "_MSAlign_ResultTable" with "_msalign"
						If strBaseName.ToLower().EndsWith("_MSAlign_ResultTable".ToLower()) Then
							strBaseName = strBaseName.Substring(0, strBaseName.Length - "_MSAlign_ResultTable".Length) & "_msalign"
						End If

						' Do not create a first-hits file for MSAlign results

						' Create the synopsis output file
						MyBase.ResetProgress("Creating the SYN file")
						Console.WriteLine()
						Console.WriteLine()
						Console.WriteLine(MyBase.ProgressStepDescription)

						strSynOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

						blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath)

						' Create the other PHRP-specific files
						MyBase.ResetProgress("Creating the PHRP files for " & System.IO.Path.GetFileName(strSynOutputFilePath))
						Console.WriteLine()
						Console.WriteLine()
						Console.WriteLine(MyBase.ProgressStepDescription)

						' Now parse the _syn.txt file that we just created to next create the other PHRP files
						blnSuccess = ParseMSAlignSynopsisFile(strSynOutputFilePath, strOutputFolderPath, lstPepToProteinMapping, False)

						If blnSuccess AndAlso mCreateProteinModsFile Then
							' Create the MTSPepToProteinMap file

							strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS:=True)

							Dim lstSourcePHRPDataFiles As Generic.List(Of String) = New Generic.List(Of String)

							If Not String.IsNullOrEmpty(strSynOutputFilePath) Then
								lstSourcePHRPDataFiles.Add(strSynOutputFilePath)
							End If

							If lstSourcePHRPDataFiles.Count = 0 Then
								SetErrorMessage("Cannot call CreatePepToProteinMapFile since lstSourcePHRPDataFiles is empty")
								SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
								blnSuccess = False
							Else
								If System.IO.File.Exists(strMTSPepToProteinMapFilePath) AndAlso mUseExistingMTSPepToProteinMapFile Then
									blnSuccess = True
								Else
									' Auto-change mIgnorePeptideToProteinMapperErrors to True
									' We only do this since a small number of peptides reported by MSAlign don't perfectly match the fasta file
									mIgnorePeptideToProteinMapperErrors = True
									blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
								End If
							End If

						End If

						If blnSuccess Then
							' If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
							MyBase.ValidatePHRPReaderSupportFiles(IO.Path.Combine(ioInputFile.DirectoryName, IO.Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath)

							' Create the Protein Mods file
							blnSuccess = MyBase.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSAlign)
						End If


						If blnSuccess Then
							MyBase.OperationComplete()
						End If

					Catch ex As Exception
						SetErrorMessage("Error calling CreateSynResultsFile: " & ex.Message)
						SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
					End Try
				End If
			End If
		Catch ex As Exception
			SetErrorMessage("Error in ProcessFile:" & ex.Message)
			SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.UnspecifiedError)
		End Try

		Return blnSuccess

	End Function

	Protected Function ReplaceTerminus(ByVal strPeptide As String) As String

		If strPeptide.StartsWith(N_TERMINUS_SYMBOL_MSALIGN) Then
			strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST & "." & strPeptide.Substring(N_TERMINUS_SYMBOL_MSALIGN.Length)
		End If

		If strPeptide.EndsWith(C_TERMINUS_SYMBOL_MSALIGN) Then
			strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_MSALIGN.Length) & "." & clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
		End If

		Return strPeptide

	End Function

	Protected Sub ResolveMSAlignModsWithModDefinitions(ByRef lstMSAlignModInfo As Generic.List(Of clsModificationDefinition))

		Dim blnExistingModFound As Boolean
		Dim objModDef As clsModificationDefinition

		If Not lstMSAlignModInfo Is Nothing Then

			' Call .LookupModificationDefinitionByMass for each entry in lstMSAlignModInfo
			For Each objModInfo As clsModificationDefinition In lstMSAlignModInfo
				If String.IsNullOrEmpty(objModInfo.TargetResidues) Then
					objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, Nothing, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnExistingModFound, True)
				Else
					For Each chTargetResidue As Char In objModInfo.TargetResidues
						objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, chTargetResidue, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnExistingModFound, True)
					Next
				End If
			Next

		End If

	End Sub

	Protected Sub StoreSearchResult(ByRef udtSearchResult As udtMSAlignSearchResultType, _
	  ByRef intFilteredSearchResultCount As Integer, _
	  ByRef udtFilteredSearchResults() As udtMSAlignSearchResultType, _
	  ByRef strErrorLog As String)

		If intFilteredSearchResultCount = udtFilteredSearchResults.Length Then
			ReDim Preserve udtFilteredSearchResults(udtFilteredSearchResults.Length * 2 - 1)
		End If

		udtFilteredSearchResults(intFilteredSearchResultCount) = udtSearchResult
		intFilteredSearchResultCount += 1

	End Sub

	Private Sub SortAndWriteFilteredSearchResults(ByRef swResultFile As System.IO.StreamWriter, _
	 ByVal intFilteredSearchResultCount As Integer, _
	 ByRef udtFilteredSearchResults() As udtMSAlignSearchResultType, _
	 ByRef strErrorLog As String)

		Dim intIndex As Integer

		' Sort udtFilteredSearchResults by ascending PVAlue, ascending scan, ascending charge, ascending peptide, and ascending protein
		Array.Sort(udtFilteredSearchResults, 0, intFilteredSearchResultCount, New MSAlignSearchResultsComparerPValueScanChargePeptide)

		For intIndex = 0 To intFilteredSearchResultCount - 1
			WriteSearchResultToFile(intIndex + 1, swResultFile, udtFilteredSearchResults(intIndex), strErrorLog)
		Next intIndex

	End Sub

	Private Sub StoreSynMatches(ByRef udtSearchResults() As udtMSAlignSearchResultType, _
	 ByVal intStartIndex As Integer, _
	 ByVal intEndIndex As Integer, _
	 ByRef intFilteredSearchResultCount As Integer, _
	 ByRef udtFilteredSearchResults() As udtMSAlignSearchResultType, _
	 ByRef strErrorLog As String)

		Dim intIndex As Integer
		Dim intCurrentCharge As Short = Short.MinValue

		AssignRankAndDeltaNormValues(udtSearchResults, intStartIndex, intEndIndex)

		' The calling procedure already sorted by scan, charge, and SpecProb; no need to re-sort

		' Now store or write out the matches that pass the filters
		For intIndex = intStartIndex To intEndIndex
			If udtSearchResults(intIndex).PValueNum <= mMSAlignSynopsisFilePValueThreshold Then
				StoreSearchResult(udtSearchResults(intIndex), intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
			End If
		Next intIndex

	End Sub

	Private Sub WriteSynFHTFileHeader(ByRef swResultFile As System.IO.StreamWriter, _
	  ByRef strErrorLog As String)

		' Write out the header line for synopsis / first hits files
		Try
			Dim lstData As New System.Collections.Generic.List(Of String)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_ResultID)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Scan)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Prsm_ID)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Spectrum_ID)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Charge)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_PrecursorMZ)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_DelM)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_DelM_PPM)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_MH)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Peptide)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Protein)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Protein_Mass)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Unexpected_Mod_Count)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Peak_Count)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Matched_Peak_Count)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Matched_Fragment_Ion_Count)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_PValue)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_Rank_PValue)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_EValue)
			lstData.Add(PHRPReader.clsPHRPParserMSAlign.DATA_COLUMN_FDR)

			swResultFile.WriteLine(CollapseList(lstData))

		Catch ex As Exception
			If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
				strErrorLog &= "Error writing synopsis / first hits header" & ControlChars.NewLine
			End If
		End Try

	End Sub

	''' <summary>
	''' Writes an entry to the synopsis file
	''' </summary>
	''' <param name="intResultID"></param>
	''' <param name="swResultFile"></param>
	''' <param name="udtSearchResult"></param>
	''' <param name="strErrorLog"></param>
	''' <remarks></remarks>
	Private Sub WriteSearchResultToFile(ByVal intResultID As Integer, _
	   ByRef swResultFile As System.IO.StreamWriter, _
	   ByRef udtSearchResult As udtMSAlignSearchResultType, _
	   ByRef strErrorLog As String)

		Try

			' Primary Columns
			'
			' MSAlign
			' ResultID  Scan  Prsm_ID  Spectrum_ID  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  Protein_Mass  Unexpected_Mod_Count  Peak_Count  Matched_Peak_Count  Matched_Fragment_Ion_Count  PValue  Rank_PValue  EValue  FDR 

			Dim lstData As New System.Collections.Generic.List(Of String)
			lstData.Add(intResultID.ToString)
			lstData.Add(udtSearchResult.ScanNum.ToString)
			lstData.Add(udtSearchResult.Prsm_ID)
			lstData.Add(udtSearchResult.Spectrum_ID)
			lstData.Add(udtSearchResult.Charge)
			lstData.Add(udtSearchResult.PrecursorMZ)
			lstData.Add(udtSearchResult.DelM)
			lstData.Add(udtSearchResult.DelM_PPM)
			lstData.Add(udtSearchResult.MH)
			lstData.Add(udtSearchResult.Peptide)
			lstData.Add(udtSearchResult.Protein)
			lstData.Add(udtSearchResult.Protein_mass)
			lstData.Add(udtSearchResult.Unexpected_modifications)	' Unexpected_Mod_Count
			lstData.Add(udtSearchResult.Peaks)						' Peak_count
			lstData.Add(udtSearchResult.Matched_peaks)				' Matched_Peak_Count
			lstData.Add(udtSearchResult.Matched_fragment_ions)		' Matched_Fragment_Ion_Count
			lstData.Add(udtSearchResult.Pvalue)
			lstData.Add(udtSearchResult.RankPValue.ToString())
			lstData.Add(udtSearchResult.Evalue)
			lstData.Add(udtSearchResult.FDR)

			swResultFile.WriteLine(CollapseList(lstData))

		Catch ex As Exception
			If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
				strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
			End If
		End Try

	End Sub

#Region "IComparer Classes"

	Protected Class MSAlignSearchResultsComparerScanChargePValuePeptide
		Implements System.Collections.Generic.IComparer(Of udtMSAlignSearchResultType)

		Public Function Compare(x As udtMSAlignSearchResultType, y As udtMSAlignSearchResultType) As Integer Implements System.Collections.Generic.IComparer(Of udtMSAlignSearchResultType).Compare

			If x.ScanNum > y.ScanNum Then
				Return 1
			ElseIf x.ScanNum < y.ScanNum Then
				Return -1
			Else
				' Scan is the same, check charge
				If x.ChargeNum > y.ChargeNum Then
					Return 1
				ElseIf x.ChargeNum < y.ChargeNum Then
					Return -1
				Else
					' Charge is the same; check Pvalue
					If x.Pvalue > y.Pvalue Then
						Return 1
					ElseIf x.Pvalue < y.Pvalue Then
						Return -1
					Else
						' Pvalue is the same; check peptide
						If x.Peptide > y.Peptide Then
							Return 1
						ElseIf x.Peptide < y.Peptide Then
							Return -1
						Else
							' Peptide is the same, check Protein
							If x.Protein > y.Protein Then
								Return 1
							ElseIf x.Protein < y.Protein Then
								Return -1
							Else
								Return 0
							End If
						End If
					End If
				End If
			End If

		End Function

	End Class

	Protected Class MSAlignSearchResultsComparerPValueScanChargePeptide
		Implements System.Collections.Generic.IComparer(Of udtMSAlignSearchResultType)

		Public Function Compare(x As udtMSAlignSearchResultType, y As udtMSAlignSearchResultType) As Integer Implements System.Collections.Generic.IComparer(Of udtMSAlignSearchResultType).Compare

			If x.Pvalue > y.Pvalue Then
				Return 1
			ElseIf x.Pvalue < y.Pvalue Then
				Return -1
			Else
				' Pvalue is the same; check scan number
				If x.ScanNum > y.ScanNum Then
					Return 1
				ElseIf x.ScanNum < y.ScanNum Then
					Return -1
				Else
					' Scan is the same, check charge
					If x.ChargeNum > y.ChargeNum Then
						Return 1
					ElseIf x.ChargeNum < y.ChargeNum Then
						Return -1
					Else
						' Charge is the same; check peptide
						If x.Peptide > y.Peptide Then
							Return 1
						ElseIf x.Peptide < y.Peptide Then
							Return -1
						Else
							' Peptide is the same, check Protein
							If x.Protein > y.Protein Then
								Return 1
							ElseIf x.Protein < y.Protein Then
								Return -1
							Else
								Return 0
							End If
						End If
					End If
				End If
			End If

		End Function

	End Class

#End Region

End Class