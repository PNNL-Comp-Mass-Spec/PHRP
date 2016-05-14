Option Strict On

' This class reads in an MODa results file (txt format) and creates 
' a tab-delimited text file with the data. 
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 04/01/2014
'
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports PHRPReader
Imports System.IO
Imports System.Text.RegularExpressions

Public Class clsMODaResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "May 19, 2015"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MODA_FILE As String = "_moda.id"

    Public Const N_TERMINUS_SYMBOL_MODA As String = "-"
    Public Const C_TERMINUS_SYMBOL_MODA As String = "-"

    Public Const DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD As Single = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    ' Note that as of April 2014, all mod masses reported by MODa are simply integers, meaning matching a trailing period is not necessary
    Private Const MODA_MOD_MASS_REGEX As String = "([+-][0-9.]+)"

    Private Const MODA_MASS_DIGITS_OF_PRECISION As Byte = 0

    Private Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline Or RegexOptions.IgnoreCase

    ' These columns correspond to the tab-delimited file (_moda.id.txt) created by MODa's anal_moda.jar file
    Private Const MODaResultsFileColCount As Integer = 11
    Public Enum eMODaResultsFileColumns As Integer
        SpectrumFileName = 0
        SpectrumIndex = 1
        ObservedMonoMass = 2
        Charge = 3
        CalculatedMonoMass = 4
        DeltaMass = 5
        Score = 6
        Probability = 7
        Peptide = 8
        Protein = 9
        PeptidePosition = 10
    End Enum

    ' These columns correspond to the Synopsis file created by this class
    Private Const MODaSynFileColCount As Integer = 15
    Public Enum eMODaSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        Spectrum_Index = 2
        Charge = 3
        PrecursorMZ = 4
        DelM = 5                            ' Precursor error, in Da
        DelM_PPM = 6                        ' Precursor error, in ppm
        MH = 7                              ' Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
        Peptide = 8                         ' This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
        Protein = 9                         ' Protein Name
        Score = 10
        Probability = 11
        Rank_Probability = 12
        Peptide_Position = 13
        QValue = 14
    End Enum

#End Region

#Region "Structures"
    ' This data structure holds rows read from the tab-delimited file (_moda.id.txt) created by MODa's anal_moda.jar file
    Private Structure udtMODaSearchResultType

        Public SpectrumFileName As String
        Public SpectrumIndex As String
        Public ScanNum As Integer               ' Determined by looking for SpectrumIndex in the _mgf_IndexToScanMap.txt file
        Public Precursor_mass As String         ' Uncharged monoisotopic mass value of the observed precursor_mz, reported as ObservedMonoMass by MODa 
        Public PrecursorMZ As String            ' Computed from ObservedMonoMass
        Public Charge As String
        Public ChargeNum As Short
        Public CalculatedMonoMass As String     ' Theoretical monoisotopic mass of the peptide (including mods), as computed by MODa
        Public DeltaMass As String              ' Computed by MODa
        Public MH As String                     ' Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
        Public DelM As String                   ' Computed using Precursor_mass - CalculatedMonoMass
        Public DelM_PPM As String               ' Computed using DelM and CalculatedMonoMass	
        Public Score As String
        Public Probability As String            ' Higher values are better
        Public ProbabilityNum As Double         ' Higher values are better
        Public RankProbability As Integer
        Public Peptide As String
        Public Protein As String
        Public PeptidePosition As String        ' Protein start/stop residues of the peptide, e.g. 108~115
        Public FDR As Double                    ' Computed by this class
        Public QValue As Double                 ' Computed by this class

        Public Sub Clear()
            SpectrumFileName = String.Empty
            SpectrumIndex = String.Empty
            ScanNum = 0
            Precursor_mass = String.Empty
            PrecursorMZ = String.Empty
            Charge = String.Empty
            ChargeNum = 0
            CalculatedMonoMass = String.Empty
            DeltaMass = String.Empty
            MH = String.Empty
            DelM = String.Empty
            DelM_PPM = String.Empty
            Score = String.Empty
            Probability = String.Empty
            ProbabilityNum = 0
            RankProbability = 0
            Peptide = String.Empty
            Protein = String.Empty
            PeptidePosition = String.Empty
            FDR = 0
            QValue = 0
        End Sub
    End Structure

#End Region

#Region "Classwide Variables"
    Private mSpectrumIndexToScanMap As Dictionary(Of Integer, Integer)
#End Region

    ''' <summary>
    ''' Step through .PeptideSequenceWithMods
    ''' For each residue, check if a static mod is defined that affects that residue
    ''' For each mod mass, determine the modification and add to objSearchResult
    ''' </summary>
    ''' <param name="objSearchResult"></param>
    ''' <param name="blnUpdateModOccurrenceCounts"></param>
    ''' <remarks></remarks>
    Private Sub AddDynamicAndStaticResidueMods(objSearchResult As clsSearchResultsMODa, blnUpdateModOccurrenceCounts As Boolean)
        Const NO_RESIDUE = "-"c

        Dim intIndex As Integer, intModIndex As Integer
        Dim chChar As Char
        Dim objModificationDefinition As clsModificationDefinition

        Dim strSequence As String

        Dim blnParsingModMass As Boolean

        Dim chMostRecentResidue As Char
        Dim intResidueLocInPeptide As Integer

        blnParsingModMass = False
        Dim strModMassDigits = String.Empty

        chMostRecentResidue = NO_RESIDUE
        intResidueLocInPeptide = 0

        strSequence = objSearchResult.PeptideSequenceWithMods
        For intIndex = 0 To strSequence.Length - 1
            chChar = strSequence.Chars(intIndex)

            If IsLetterAtoZ(chChar) Then

                If blnParsingModMass Then
                    ' Associate the mod mass in strModMassDigits with the previous residue
                    AssociateDynamicModWithResidue(objSearchResult, chMostRecentResidue, intResidueLocInPeptide, strModMassDigits, blnUpdateModOccurrenceCounts)
                    blnParsingModMass = False
                End If

                chMostRecentResidue = chChar
                intResidueLocInPeptide += 1

                ' Look for static mods to associate with this residue
                For intModIndex = 0 To mPeptideMods.ModificationCount - 1
                    If mPeptideMods.GetModificationTypeByIndex(intModIndex) = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                        objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)

                        If objModificationDefinition.TargetResiduesContain(chChar) Then
                            ' Match found; add this modification
                            objSearchResult.SearchResultAddModification(objModificationDefinition, chChar, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                        End If
                    End If
                Next intModIndex

            Else
                Dim blnIsNumberChar As Boolean = (chChar = "+"c OrElse chChar = "-"c OrElse Char.IsDigit(chChar))

                If blnParsingModMass Then
                    If blnIsNumberChar OrElse chChar = "."c Then
                        strModMassDigits &= chChar
                    End If
                ElseIf blnIsNumberChar Then
                    ' Mod Mass Start
                    strModMassDigits = chChar
                    blnParsingModMass = True
                Else
                    ' Unrecognized symbol; ignore it
                End If

            End If

        Next intIndex

        If blnParsingModMass Then
            ' Associate the mod mass in strModMassDigits with the previous residue
            AssociateDynamicModWithResidue(objSearchResult, chMostRecentResidue, intResidueLocInPeptide, strModMassDigits, blnUpdateModOccurrenceCounts)
        End If

    End Sub

    Private Function AddModificationsAndComputeMass(objSearchResult As clsSearchResultsMODa, blnUpdateModOccurrenceCounts As Boolean) As Boolean
        Const ALLOW_DUPLICATE_MOD_ON_TERMINUS = True

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

    '' This function was an experiment to compute better DelM_PPM values
    '' by reading the synopsis file with PHRPReader and re-computing the DelM_PPM values based on the monoisotopic mass values computed for the sequences
    '' It turned out to not be required, since the DelM_PPM values reported by MODa are quite accurate (despite the fact that it reports integer mod mass values)
    'Private Function AppendDelMPPMRefinedToSynFile(strSynOutputFilePath As String) As Boolean

    '   Const SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED As String = "DelM_PPM_Refined"

    '	Dim blnSuccess As Boolean

    '	Try

    '		' Keys in this dictionary are ResultID values from the synopsis file
    '		' Values are refined DelM_PPM values
    '		Dim dctRefinedDelMPPMErrors = New Dictionary(Of Integer, Double)

    '		Using objReader As New clsPHRPReader(strSynOutputFilePath, clsPHRPReader.ePeptideHitResultType.MODa, blnLoadModsAndSeqInfo:=True, blnLoadMSGFResults:=False)
    '			objReader.EchoMessagesToConsole = True
    '			objReader.SkipDuplicatePSMs = True
    '			objReader.SkipDuplicatePSMs = False

    '			For Each strErrorMessage As String In objReader.ErrorMessages
    '				SetErrorMessage(strErrorMessage)
    '			Next

    '			For Each strWarningMessage As String In objReader.WarningMessages
    '				ReportWarning(strWarningMessage)
    '			Next

    '			objReader.ClearErrors()
    '			objReader.ClearWarnings()

    '			AddHandler objReader.ErrorEvent, AddressOf PHRPReader_ErrorEvent
    '			AddHandler objReader.WarningEvent, AddressOf PHRPReader_WarningEvent

    '			Do While objReader.MoveNext()

    '				Dim oPSM = objReader.CurrentPSM
    '				Dim oSeqInfo = objReader.CurrentPSMSeqInfo()

    '				If Not oSeqInfo Is Nothing Then
    '					Dim dblDelM = oPSM.PrecursorNeutralMass - oSeqInfo.MonoisotopicMass

    '					Dim dblPeptideDeltaMassRefinedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, oPSM.PrecursorNeutralMass, True, oSeqInfo.MonoisotopicMass)

    '					Dim dblOriginalDelMPPM As Double
    '					If Double.TryParse(oPSM.MassErrorPPM, dblOriginalDelMPPM) Then
    '						If Math.Abs(dblPeptideDeltaMassRefinedPpm - dblOriginalDelMPPM) > 2 Then
    '							Console.WriteLine("Computed a refined DelMPPM value: " & dblPeptideDeltaMassRefinedPpm.ToString("0.0") & " vs. " & dblOriginalDelMPPM.ToString("0.0"))
    '						End If
    '					End If

    '					dctRefinedDelMPPMErrors.Add(oPSM.ResultID, dblPeptideDeltaMassRefinedPpm)
    '				End If

    '			Loop

    '			RemoveHandler objReader.ErrorEvent, AddressOf PHRPReader_ErrorEvent
    '			RemoveHandler objReader.WarningEvent, AddressOf PHRPReader_WarningEvent

    '		End Using

    '		Dim strSynOutputFilePathNew = strSynOutputFilePath & ".refinedDelMPPM"
    '		Dim blnHeadersParsed As Boolean = False
    '		Dim blnSwapFiles As Boolean = True

    '		Using srDataFile = New StreamReader(New FileStream(strSynOutputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
    '			Using swOutfile = New StreamWriter(New FileStream(strSynOutputFilePathNew, FileMode.Create, FileAccess.Write, FileShare.Read))

    '				Do While Not srDataFile.EndOfStream
    '					Dim strLineIn = srDataFile.ReadLine

    '					Dim strSplitLine = strLineIn.Split(ControlChars.Tab)

    '					If Not blnHeadersParsed Then
    '						If strSplitLine.Contains(SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED) Then
    '							' This file already has the refined DelM_PPM column
    '							' Do not update it
    '							blnSwapFiles = False
    '							Exit Do
    '						End If
    '						swOutfile.WriteLine(strLineIn & ControlChars.Tab & SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED)
    '						blnHeadersParsed = True
    '					Else
    '						Dim resultID As Integer
    '						Dim strDelMPPMRefined As String = String.Empty

    '						If Integer.TryParse(strSplitLine(0), resultID) Then

    '							Dim delMPPMRefined As Double
    '							If dctRefinedDelMPPMErrors.TryGetValue(resultID, delMPPMRefined) Then
    '								strDelMPPMRefined = NumToString(delMPPMRefined, 5, True)
    '							End If

    '						End If

    '						swOutfile.WriteLine(strLineIn & ControlChars.Tab & strDelMPPMRefined)
    '					End If
    '				Loop
    '			End Using
    '		End Using

    '		If blnSwapFiles Then
    '			Threading.Thread.Sleep(150)

    '			Try
    '				' Replace the original synopsis file with the updated one

    '				File.Delete(strSynOutputFilePath)
    '				Threading.Thread.Sleep(150)

    '				File.Move(strSynOutputFilePathNew, strSynOutputFilePath)

    '				blnSuccess = True

    '			Catch ex As Exception
    '				SetErrorMessage("Exception adding column " & SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED & " to the synopsis file: " & ex.Message)
    '				blnSuccess = False
    '			End Try

    '		End If

    '	Catch ex As Exception
    '		SetErrorMessage(ex.Message)
    '		SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
    '		blnSuccess = False
    '	End Try

    '	Return blnSuccess

    'End Function

    Private Sub AssociateDynamicModWithResidue(
      objSearchResult As clsSearchResultsMODa,
      chMostRecentResidue As Char,
      intResidueLocInPeptide As Integer,
      strModMassDigits As String,
      blnUpdateModOccurrenceCounts As Boolean)

        Dim blnSuccess As Boolean

        Dim chResidueForMod As Char
        Dim intResidueLocForMod As Integer
        Dim dblModMass As Double

        chResidueForMod = chMostRecentResidue
        intResidueLocForMod = intResidueLocInPeptide

        If Double.TryParse(strModMassDigits, dblModMass) Then
            If intResidueLocForMod = 0 Then
                ' Modification is at the peptide N-terminus
                intResidueLocForMod = 1
            End If

            blnSuccess = objSearchResult.SearchResultAddModification(dblModMass, chResidueForMod, intResidueLocForMod, objSearchResult.DetermineResidueTerminusState(intResidueLocForMod), blnUpdateModOccurrenceCounts, MODA_MASS_DIGITS_OF_PRECISION)

            If Not blnSuccess Then
                Dim strErrorMessage As String = objSearchResult.ErrorMessage
                If String.IsNullOrEmpty(strErrorMessage) Then
                    strErrorMessage = "SearchResultAddDynamicModification returned false for mod mass " & strModMassDigits
                End If
                SetErrorMessage(strErrorMessage & "; ResultID = " & objSearchResult.ResultID)
            End If
        End If

    End Sub

    ''' <summary>
    ''' Ranks each entry  assumes all of the data is from the same scan)
    ''' </summary>
    ''' <param name="lstSearchResults"></param>
    ''' <param name="intStartIndex"></param>
    ''' <param name="intEndIndex"></param>
    ''' <remarks></remarks>
    Private Sub AssignRankAndDeltaNormValues(
      lstSearchResults As List(Of udtMODaSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer)

        ' Prior to September 2014 ranks were assign per charge state per scan; 
        ' Ranks are now assigned per scan (across all charge states)

        ' Duplicate a portion of lstSearchResults so that we can sort by descending Probability

        Dim dctResultsSubset = New Dictionary(Of Integer, udtMODaSearchResultType)
        For intIndex = intStartIndex To intEndIndex
            dctResultsSubset.Add(intIndex, lstSearchResults(intIndex))
        Next

        Dim lstResultsByProbability = (From item In dctResultsSubset Select item Order By item.Value.ProbabilityNum Descending).ToList()

        Dim dblLastValue As Double
        Dim intCurrentRank As Integer = -1

        For Each entry In lstResultsByProbability
            Dim oResult = lstSearchResults(entry.Key)

            If intCurrentRank < 0 Then
                dblLastValue = oResult.ProbabilityNum
                intCurrentRank = 1
            Else
                If Math.Abs(oResult.ProbabilityNum - dblLastValue) > Double.Epsilon Then
                    dblLastValue = oResult.ProbabilityNum
                    intCurrentRank += 1
                End If
            End If

            oResult.RankProbability = intCurrentRank
            lstSearchResults(entry.Key) = oResult
        Next

    End Sub

    Private Function AssureInteger(strInteger As String, intDefaultValue As Integer) As String

        Dim intValue As Integer
        Dim dblValue As Double

        If strInteger.EndsWith(".0") Then strInteger = strInteger.Substring(0, strInteger.Length - 2)

        If Integer.TryParse(strInteger, intValue) Then
            Return intValue.ToString()
        ElseIf Double.TryParse(strInteger, dblValue) Then
            Return dblValue.ToString("0")
        Else
            Return intDefaultValue.ToString()
        End If

    End Function

    Private Function ComputePeptideMass(strPeptide As String, dblTotalModMass As Double) As Double

        Dim strCleanSequence = GetCleanSequence(strPeptide)

        Dim dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence)

        If Math.Abs(dblTotalModMass) > Double.Epsilon Then
            dblMass += dblTotalModMass
        End If

        Return dblMass

    End Function

    ''' <summary>
    ''' Computes the total of all modification masses defined for the peptide
    ''' </summary>
    ''' <param name="strPeptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function ComputeTotalModMass(strPeptide As String) As Double

        Static reModMassRegEx As New Regex(MODA_MOD_MASS_REGEX, REGEX_OPTIONS)

        Dim dblTotalModMass As Double = 0

        Dim strPrimarySequence As String = String.Empty
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, strPrimarySequence, strPrefix, strSuffix)

        ' Parse the dynamic mods reported by MODa
        For Each reMatch As Match In reModMassRegEx.Matches(strPrimarySequence)
            Dim dblModMassFound As Double
            ' We use .TrimEnd() because the matched mod mass will end in a period if this mod applies to the final residue in a peptide
            If Double.TryParse(reMatch.Groups(1).Value.TrimEnd("."c), dblModMassFound) Then
                dblTotalModMass += dblModMassFound
            End If
        Next

        ' Now look for static mods 
        ' First determine the index of the last residue in strPrimarySequence
        Dim intIndexLastChar As Integer = strPrimarySequence.Length

        For intIndex = strPrimarySequence.Length - 1 To 0 Step -1
            If IsLetterAtoZ(strPrimarySequence.Chars(intIndex)) Then
                intIndexLastChar = intIndex
                Exit For
            End If
        Next

        For intIndex = 0 To strPrimarySequence.Length - 1
            Dim chChar = strPrimarySequence.Chars(intIndex)

            If IsLetterAtoZ(chChar) Then

                ' Look for static mods to associate with this residue
                For intModIndex = 0 To mPeptideMods.ModificationCount - 1
                    If mPeptideMods.GetModificationTypeByIndex(intModIndex) = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                        Dim objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)
                        Dim blnMatchFound = objModificationDefinition.TargetResiduesContain(chChar)

                        If Not blnMatchFound AndAlso intIndex = 0 Then
                            blnMatchFound = objModificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS)
                        End If

                        If Not blnMatchFound AndAlso intIndex = intIndexLastChar Then
                            blnMatchFound = objModificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                        End If

                        If blnMatchFound Then
                            dblTotalModMass += objModificationDefinition.ModificationMass
                        End If
                    End If
                Next intModIndex
            End If

        Next

        Return dblTotalModMass

    End Function

    Protected Overrides Function ConstructPepToProteinMapFilePath(strInputFilePath As String, strOutputFolderPath As String, MTS As Boolean) As String

        Dim strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath)
        If strPepToProteinMapFilePath.ToLower().EndsWith("_MODa_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_MODa_fht") Then
            ' Remove _syn or _fht
            strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4)
        End If

        Return MyBase.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS)

    End Function

    ''' <summary>
    ''' This routine creates a first hits file or synopsis file from the output from MODa
    ''' The synopsis file includes every result with a probability above a set threshold
    ''' The first-hits file includes the result with the highest probability (for each scan and charge)
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function CreateSynResultsFile(
      strInputFilePath As String,
      strOutputFilePath As String) As Boolean

        Try
            Dim intColumnMapping() As Integer = Nothing
            Dim strErrorLog = String.Empty

            ' Open the input file and parse it
            ' Initialize the stream reader and the stream Text writer
            Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)),
                  swResultFile = New StreamWriter(New FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                Dim intResultsProcessed = 0

                ' Initialize the list that will hold all of the records in the MODa result file
                Dim lstSearchResultsUnfiltered = New List(Of udtMODaSearchResultType)

                ' Initialize the list that will hold all of the records that will ultimately be written out to disk
                Dim lstFilteredSearchResults = New List(Of udtMODaSearchResultType)

                ' Parse the input file
                Do While Not srDataFile.EndOfStream And Not MyBase.AbortProcessing
                    Dim strLineIn = srDataFile.ReadLine()
                    Dim blnSkipLine = False

                    If String.IsNullOrWhiteSpace(strLineIn) Then
                        Continue Do
                    End If

                    If intResultsProcessed = 0 Then
                        ' The first line might be a header line
                        ' However, as of April 2014, MODa id.txt files do not have a header line

                        blnSkipLine = ParseMODaResultsFileHeaderLine(strLineIn, intColumnMapping)

                        ' Write the header line to the output file
                        WriteSynFHTFileHeader(swResultFile, strErrorLog)
                    End If

                    If Not blnSkipLine Then

                        Dim udtSearchResult = New udtMODaSearchResultType

                        Dim blnValidSearchResult = ParseMODaResultsFileEntry(strLineIn, udtSearchResult, strErrorLog, intColumnMapping)

                        If blnValidSearchResult Then
                            lstSearchResultsUnfiltered.Add(udtSearchResult)
                        End If

                        ' Update the progress
                        Dim sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
                        If mCreateProteinModsFile Then
                            sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
                        End If
                        UpdateProgress(sngPercentComplete)

                    End If

                    intResultsProcessed += 1

                Loop

                ' Sort the SearchResults by scan, charge, and Descending probability
                lstSearchResultsUnfiltered.Sort(New MODaSearchResultsComparerScanChargeProbabilityPeptide)

                ' Now filter the data

                ' Initialize variables
                Dim intStartIndex = 0
                Dim intEndIndex As Integer

                Do While intStartIndex < lstSearchResultsUnfiltered.Count
                    intEndIndex = intStartIndex
                    Do While intEndIndex + 1 < lstSearchResultsUnfiltered.Count AndAlso lstSearchResultsUnfiltered(intEndIndex + 1).ScanNum = lstSearchResultsUnfiltered(intStartIndex).ScanNum
                        intEndIndex += 1
                    Loop

                    ' Store the results for this scan
                    StoreSynMatches(lstSearchResultsUnfiltered, intStartIndex, intEndIndex, lstFilteredSearchResults)

                    intStartIndex = intEndIndex + 1
                Loop

                ' Sort the data in udtFilteredSearchResults then write out to disk
                SortAndWriteFilteredSearchResults(swResultFile, lstFilteredSearchResults, strErrorLog)

            End Using

            ' Inform the user if any errors occurred
            If strErrorLog.Length > 0 Then
                SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
            End If

            Return True

        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            Return False
        End Try

    End Function

    ''' <summary>
    ''' Load the static mods defined in the MODa parameter file
    ''' </summary>
    ''' <param name="strMODaParamFilePath"></param>
    ''' <param name="lstModInfo"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function ExtractModInfoFromMODaParamFile(strMODaParamFilePath As String, ByRef lstModInfo As List(Of clsModificationDefinition)) As Boolean

        Dim strLineIn As String
        Dim kvSetting As KeyValuePair(Of String, String)

        Dim objModDef As clsModificationDefinition

        Dim blnSuccess = False

        Try
            ' Initialize the modification list
            If lstModInfo Is Nothing Then
                lstModInfo = New List(Of clsModificationDefinition)
            Else
                lstModInfo.Clear()
            End If

            If String.IsNullOrEmpty(strMODaParamFilePath) Then
                SetErrorMessage("MODa Parameter File name not defined; unable to extract mod info")
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                Return False
            End If

            If Not File.Exists(strMODaParamFilePath) Then
                SetErrorMessage("MODa param file not found: " & strMODaParamFilePath)
            Else
                ' Read the contents of the parameter (or mods) file
                Using srInFile = New StreamReader(New FileStream(strMODaParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                    Do While Not srInFile.EndOfStream
                        strLineIn = srInFile.ReadLine().Trim()

                        If strLineIn.Length = 0 Then
                            Continue Do
                        End If

                        If strLineIn.StartsWith("#"c) Then
                            ' Comment line; skip it
                            Continue Do
                        End If

                        ' Split the line on the equals sign
                        kvSetting = clsPHRPParser.ParseKeyValueSetting(strLineIn, "="c, "#")

                        If String.IsNullOrEmpty(kvSetting.Key) Then
                            Continue Do
                        End If

                        If String.Equals(kvSetting.Key, "add", StringComparison.InvariantCultureIgnoreCase) Then
                            ' ModA defines all of its static modifications with the ADD keyword
                            ' Split the value at the comma and create a new setting entry with the residue name

                            Dim strValue = kvSetting.Value
                            Dim commaIndex = strValue.IndexOf(","c)

                            Dim strResidue = strValue.Substring(0, commaIndex).Trim()
                            strValue = strValue.Substring(commaIndex + 1).Trim()

                            ' Replace NTerm or CTerm with < or >
                            If strResidue.ToLower() = "nterm" Then strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS
                            If strResidue.ToLower() = "cterm" Then strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS

                            Dim modMass As Double = 0
                            If Double.TryParse(strValue, modMass) Then
                                If Math.Abs(modMass - 0) > Single.Epsilon Then

                                    Dim strMassCorrectionTag As String = mPeptideMods.LookupMassCorrectionTagByMass(modMass)

                                    objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMass, strResidue, clsModificationDefinition.eModificationTypeConstants.StaticMod, strMassCorrectionTag)
                                    lstModInfo.Add(objModDef)

                                End If
                            End If

                        End If

                    Loop
                End Using

                Console.WriteLine()

                blnSuccess = True

            End If
        Catch ex As Exception
            SetErrorMessage("Error reading the MODa parameter file (" & Path.GetFileName(strMODaParamFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Sub InitializeLocalVariables()
        mSpectrumIndexToScanMap = New Dictionary(Of Integer, Integer)
    End Sub

    Private Function LoadMGFIndexToScanMapFile(fiInputFile As FileInfo) As Boolean

        Dim indexToScanMapFilePath As String = String.Empty

        Try
            mSpectrumIndexToScanMap.Clear()

            ' Look for the IndexToScanMap file that corresponds to fiInputFile
            Dim lstScanMapFiles As List(Of FileInfo)
            Dim matchIndex = fiInputFile.Name.LastIndexOf("_moda", StringComparison.Ordinal)
            Dim sourceFileDescription As String

            If matchIndex > 0 Then
                Dim datasetName = fiInputFile.Name.Substring(0, matchIndex)
                lstScanMapFiles = fiInputFile.Directory.GetFiles(datasetName & "*mgf_IndexToScanMap*").ToList()
                sourceFileDescription = " dataset " & datasetName
            Else
                ' Source file does not have "_moda" in the name
                ' Look for any mgf_IndexToScanMap file
                lstScanMapFiles = fiInputFile.Directory.GetFiles("*mgf_IndexToScanMap*").ToList()
                sourceFileDescription = fiInputFile.Name
            End If

            If lstScanMapFiles.Count = 1 Then
                indexToScanMapFilePath = lstScanMapFiles.First.FullName
            ElseIf lstScanMapFiles.Count = 0 Then
                ReportWarning("Did not find a mgf_IndexToScanMap file for " & sourceFileDescription & " in folder " & fiInputFile.Directory.FullName & "; scan numbers will be 0 in the synopsis file")
                Return False
            Else
                ReportWarning("Found more than one potential mgf_IndexToScanMap file for " & sourceFileDescription & " in folder " & fiInputFile.Directory.FullName & " scan numbers will be 0 in the synopsis file")
                Return False
            End If


            Dim fiSourceFile = New FileInfo(indexToScanMapFilePath)

            If Not fiSourceFile.Exists Then
                ReportWarning("MGF Index to Scan Map file not found; scan numbers will be 0 in the synopsis file: " & indexToScanMapFilePath)
                Return False
            End If

            Using srMapFile = New StreamReader(New FileStream(fiSourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                While Not srMapFile.EndOfStream
                    Dim strLineIn = srMapFile.ReadLine()
                    If String.IsNullOrEmpty(strLineIn) Then Continue While

                    Dim strSplitLine = strLineIn.Split(ControlChars.Tab)
                    If strSplitLine.Length >= 3 Then
                        Dim spectrumIndex As Integer
                        Dim scanStart As Integer
                        Dim scanEnd As Integer

                        If Integer.TryParse(strSplitLine(0), spectrumIndex) Then
                            If Integer.TryParse(strSplitLine(1), scanStart) Then
                                Integer.TryParse(strSplitLine(2), scanEnd)

                                mSpectrumIndexToScanMap.Add(spectrumIndex, scanStart)

                            End If
                        End If
                    End If
                End While
            End Using

        Catch ex As Exception
            SetErrorMessage("Error reading the MGF Index to Scan Map file (" & Path.GetFileName(indexToScanMapFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            Return False
        End Try

        Return True

    End Function

    Private Function LookupScanBySpectrumIndex(spectrumIndex As Integer) As Integer

        Dim scanNumber As Integer
        If mSpectrumIndexToScanMap.TryGetValue(spectrumIndex, scanNumber) Then
            Return scanNumber
        End If

        Return 0
    End Function

    Private Function ParseMODaSynopsisFile(
      strInputFilePath As String,
      strOutputFolderPath As String,
      lstPepToProteinMapping As List(Of udtPepToProteinMappingType),
      blnResetMassCorrectionTagsAndModificationDefinitions As Boolean) As Boolean

        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim strPreviousProbability As String

        ' Note that MODa synopsis files are normally sorted on Probability value, ascending
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
        '  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

        Dim intColumnMapping() As Integer = Nothing
        Dim blnSuccess As Boolean
        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            ' Initialize objSearchResult
            Dim objSearchResult = New clsSearchResultsMODa(mPeptideMods)

            ' Initialize htPeptidesFoundForProbabilityLevel
            Dim htPeptidesFoundForProbabilityLevel = New Hashtable
            strPreviousProbability = String.Empty

            ' Assure that lstPepToProteinMapping is sorted on peptide
            If lstPepToProteinMapping.Count > 1 Then
                lstPepToProteinMapping.Sort(New PepToProteinMappingComparer)
            End If

            Try
                objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(mEnzymeMatchSpec, mPeptideNTerminusMassChange, mPeptideCTerminusMassChange)

                Dim strErrorLog = String.Empty

                ' Open the input file and parse it
                ' Initialize the stream reader
                Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                    Dim intResultsProcessed = 0
                    Dim blnHeaderParsed = False

                    ' Create the output files
                    Dim strBaseOutputFilePath As String = Path.Combine(strOutputFolderPath, Path.GetFileName(strInputFilePath))
                    blnSuccess = MyBase.InitializeSequenceOutputFiles(strBaseOutputFilePath)

                    ' Parse the input file
                    Do While Not srDataFile.EndOfStream And Not MyBase.AbortProcessing

                        Dim strLineIn = srDataFile.ReadLine()
                        If String.IsNullOrWhiteSpace(strLineIn) Then
                            Continue Do
                        End If

                        If Not blnHeaderParsed Then
                            blnSuccess = ParseMODaSynFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True
                            Continue Do
                        End If

                        Dim strCurrentPeptideWithMods As String = String.Empty

                        Dim blnValidSearchResult = ParseMODaSynFileEntry(
                          strLineIn, objSearchResult, strErrorLog,
                          intResultsProcessed, intColumnMapping,
                          strCurrentPeptideWithMods)

                        If Not blnValidSearchResult Then
                            Continue Do
                        End If

                        Dim strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge
                        Dim blnFirstMatchForGroup As Boolean

                        If objSearchResult.Probability = strPreviousProbability Then
                            ' New result has the same Probability as the previous result
                            ' See if htPeptidesFoundForProbabilityLevel contains the peptide, scan and charge

                            If htPeptidesFoundForProbabilityLevel.ContainsKey(strKey) Then
                                blnFirstMatchForGroup = False
                            Else
                                htPeptidesFoundForProbabilityLevel.Add(strKey, 1)
                                blnFirstMatchForGroup = True
                            End If

                        Else
                            ' New Probability
                            ' Reset htPeptidesFoundForProbabilityLevel
                            htPeptidesFoundForProbabilityLevel.Clear()

                            ' Update strPreviousProbability
                            strPreviousProbability = objSearchResult.Probability

                            ' Append a new entry to htPeptidesFoundForProbabilityLevel
                            htPeptidesFoundForProbabilityLevel.Add(strKey, 1)
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
                            Dim intPepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(lstPepToProteinMapping, strCurrentPeptideWithMods)

                            If intPepToProteinMapIndex >= 0 Then
                                ' Call MyBase.SaveResultsFileEntrySeqInfo for each entry in lstPepToProteinMapping() for peptide , skipping objSearchResult.ProteinName
                                Dim strCurrentProtein = String.Copy(objSearchResult.ProteinName)
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

                        ' Update the progress
                        Dim sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
                        If mCreateProteinModsFile Then
                            sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
                        End If
                        UpdateProgress(sngPercentComplete)

                        intResultsProcessed += 1

                    Loop

                End Using

                If mCreateModificationSummaryFile Then
                    ' Create the modification summary file
                    Dim fiInputFile = New FileInfo(strInputFilePath)
                    Dim strModificationSummaryFilePath = Path.GetFileName(MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_MOD_SUMMARY))
                    strModificationSummaryFilePath = Path.Combine(strOutputFolderPath, strModificationSummaryFilePath)

                    SaveModificationSummaryFile(strModificationSummaryFilePath)
                End If

                ' Inform the user if any errors occurred
                If strErrorLog.Length > 0 Then
                    SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
                End If

                blnSuccess = True

            Catch ex As Exception
                SetErrorMessage(ex.Message)
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
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

    Private Function ParseMODaResultsFileEntry(
      strLineIn As String,
      ByRef udtSearchResult As udtMODaSearchResultType,
      ByRef strErrorLog As String,
      intColumnMapping() As Integer) As Boolean

        ' Parses an entry from the MODa results file

        Dim rowIndex = "?"
        Dim dblPrecursorMonoMass As Double          ' Observed m/z, converted to monoisotopic mass
        Dim dblPeptideMonoMassMODa As Double        ' Theoretical peptide monoisotopic mass, including mods, as computed by MODa
        Dim dblPeptideMonoMassPHRP As Double        ' Theoretical peptide monoisotopic mass, including mods, as computed by PHRP

        Dim dblPrecursorMZ As Double
        Dim dblDelM As Double

        Dim dblTotalModMass As Double

        Dim blnValidSearchResult As Boolean

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            udtSearchResult.Clear()
            Dim strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length >= 11 Then

                With udtSearchResult
                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.SpectrumFileName), .SpectrumFileName)

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.SpectrumIndex), .SpectrumIndex) Then
                        ReportError("Index column is missing or invalid", True)
                    Else
                        rowIndex = .SpectrumIndex
                    End If

                    Dim spectrumIndex As Integer
                    If Not Integer.TryParse(.SpectrumIndex, spectrumIndex) Then
                        ReportError("Index column is not numeric", True)
                    End If
                    .ScanNum = LookupScanBySpectrumIndex(spectrumIndex)
                    If .ScanNum = 0 Then
                        ReportWarning("Error, could not resolve spectrumIndex to Scan Number: " & spectrumIndex)
                    End If

                    ' Monoisotopic mass value of the observed precursor_mz
                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.ObservedMonoMass), .Precursor_mass)
                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.Charge), .Charge)
                    .ChargeNum = CShort(CIntSafe(.Charge, 0))

                    If Double.TryParse(.Precursor_mass, dblPrecursorMonoMass) Then
                        If .ChargeNum > 0 Then
                            dblPrecursorMZ = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMonoMass, 0, .ChargeNum)
                            .PrecursorMZ = NumToString(dblPrecursorMZ, 6, True)
                        End If
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.CalculatedMonoMass), .CalculatedMonoMass)
                    Double.TryParse(.CalculatedMonoMass, dblPeptideMonoMassMODa)

                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.DeltaMass), .DeltaMass)
                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.Score), .Score)

                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.Probability), .Probability)
                    If Not Double.TryParse(.Probability, .ProbabilityNum) Then .ProbabilityNum = 0

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.Peptide), .Peptide) Then
                        ReportError("Peptide column is missing or invalid", True)
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.Protein), .Protein)
                    .Protein = TruncateProteinName(.Protein)
                    GetColumnValue(strSplitLine, intColumnMapping(eMODaResultsFileColumns.PeptidePosition), .PeptidePosition)

                    ' Parse the sequence to determine the total mod mass
                    ' Note that we do not remove any of the mod symbols since MODa identifies mods by mass alone
                    ' Note that static mods are implied (thus are not explicitly displayed by MODa)
                    dblTotalModMass = ComputeTotalModMass(.Peptide)

                    ' Compute monoisotopic mass of the peptide
                    dblPeptideMonoMassPHRP = ComputePeptideMass(.Peptide, dblTotalModMass)

                    If Math.Abs(dblPeptideMonoMassMODa) < Double.Epsilon Then
                        dblPeptideMonoMassMODa = dblPeptideMonoMassPHRP
                    End If

                    Dim dblMassDiffThreshold As Double = dblPeptideMonoMassMODa / 50000
                    If dblMassDiffThreshold < 0.1 Then dblMassDiffThreshold = 0.1

                    If Math.Abs(dblPeptideMonoMassPHRP - dblPeptideMonoMassMODa) > dblMassDiffThreshold Then
                        ' Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da or by a slightly larger value if over 5000 Da; this is unexpected
                        Dim strFirst30Residues As String
                        If .Peptide.Length < 27 Then
                            strFirst30Residues = .Peptide
                        Else
                            strFirst30Residues = .Peptide.Substring(0, 27) & "..."
                        End If
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " & dblMassDiffThreshold.ToString("0.00") & " Da away from the mass computed by MODa: " & dblPeptideMonoMassPHRP.ToString("0.0000") & " vs. " & dblPeptideMonoMassMODa.ToString("0.0000") & "; peptide " & strFirst30Residues)
                    End If

                    If dblPeptideMonoMassMODa > 0 Then
                        ' Compute DelM and DelM_PPM
                        dblDelM = dblPrecursorMonoMass - dblPeptideMonoMassMODa
                        .DelM = NumToString(dblDelM, 6, True)

                        Dim dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, dblPrecursorMonoMass, True, dblPeptideMonoMassMODa)

                        .DelM_PPM = NumToString(dblPeptideDeltaMassCorrectedPpm, 5, True)
                    End If

                    ' Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    .MH = NumToString(clsPeptideMassCalculator.ConvoluteMass(dblPeptideMonoMassPHRP, 0, 1), 6, True)

                    Dim dblProbability As Double

                    If .Probability.ToLower() = "infinity" Then
                        .Probability = "0"
                    ElseIf Not String.IsNullOrEmpty(.Probability) And Not Double.TryParse(.Probability, dblProbability) Then
                        .Probability = ""
                    End If

                End With

                blnValidSearchResult = True
            End If

        Catch ex As Exception
            ' Error parsing this row from the MODa results file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not String.IsNullOrEmpty(rowIndex) Then
                    strErrorLog &= "Error parsing MODa Results in ParseMODaResultsFileEntry for RowIndex '" & rowIndex & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MODa Results in ParseMODaResultsFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <param name="strLineIn"></param>
    ''' <param name="intColumnMapping"></param>
    ''' <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
    ''' <remarks></remarks>
    Private Function ParseMODaResultsFileHeaderLine(strLineIn As String, ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        ' The expected column order from MODa:
        '   SpectrumFile	Index	ObservedMonoMass	Charge	CalculatedMonoMass	DeltaMass	Score	Probability	Peptide	Protein	PeptidePosition

        Dim lstColumnNames As SortedDictionary(Of String, eMODaResultsFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMODaResultsFileColumns)(StringComparer.InvariantCultureIgnoreCase)

        ReDim intColumnMapping(MODaResultsFileColCount - 1)

        lstColumnNames.Add("SpectrumFile", eMODaResultsFileColumns.SpectrumFileName)
        lstColumnNames.Add("Index", eMODaResultsFileColumns.SpectrumIndex)
        lstColumnNames.Add("ObservedMW", eMODaResultsFileColumns.ObservedMonoMass)
        lstColumnNames.Add("Charge", eMODaResultsFileColumns.Charge)
        lstColumnNames.Add("CalculatedMW", eMODaResultsFileColumns.CalculatedMonoMass)
        lstColumnNames.Add("DeltaMass", eMODaResultsFileColumns.DeltaMass)
        lstColumnNames.Add("Score", eMODaResultsFileColumns.Score)
        lstColumnNames.Add("Probability", eMODaResultsFileColumns.Probability)
        lstColumnNames.Add("Peptide", eMODaResultsFileColumns.Peptide)
        lstColumnNames.Add("Protein", eMODaResultsFileColumns.Protein)
        lstColumnNames.Add("PeptidePosition", eMODaResultsFileColumns.PeptidePosition)

        Try
            ' Initialize each entry in intColumnMapping to -1
            For intIndex = 0 To intColumnMapping.Length - 1
                intColumnMapping(intIndex) = -1
            Next

            Dim strSplitLine = strLineIn.Split(ControlChars.Tab)
            Dim blnUseDefaultHeaders = False

            Dim value As Integer
            If strSplitLine.Length >= 2 Then
                If Integer.TryParse(strSplitLine(1), value) Then
                    ' Second column has a number; this is not a header line					
                    blnUseDefaultHeaders = True
                Else

                    For intIndex = 0 To strSplitLine.Length - 1
                        Dim eResultFileColumn As eMODaResultsFileColumns

                        If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
                            ' Recognized column name; update intColumnMapping
                            intColumnMapping(eResultFileColumn) = intIndex
                            blnUseDefaultHeaders = False
                        Else
                            ' Unrecognized column name
                            Console.WriteLine("Warning: Unrecognized column header name '" & strSplitLine(intIndex) & "' in ParseMODaResultsFileHeaderLine")
                        End If
                    Next

                End If

            End If

            If blnUseDefaultHeaders Then
                ' Use default column mappings
                For intIndex = 0 To intColumnMapping.Length - 1
                    intColumnMapping(intIndex) = intIndex
                Next

                ' This is not a header line; return false
                Return False
            End If

        Catch ex As Exception
            SetErrorMessage("Error parsing header in MODa results file: " & ex.Message)
            Return False
        End Try

        ' Header line found and parsed; return true
        Return True

    End Function

    Private Function ParseMODaSynFileHeaderLine(strLineIn As String, ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMODaSynFileColumns
        Dim lstColumnNames As SortedDictionary(Of String, eMODaSynFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMODaSynFileColumns)(StringComparer.InvariantCultureIgnoreCase)

        ReDim intColumnMapping(MODaSynFileColCount - 1)

        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_ResultID, eMODaSynFileColumns.ResultID)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Scan, eMODaSynFileColumns.Scan)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Spectrum_Index, eMODaSynFileColumns.Spectrum_Index)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Charge, eMODaSynFileColumns.Charge)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_PrecursorMZ, eMODaSynFileColumns.PrecursorMZ)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_DelM, eMODaSynFileColumns.DelM)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_DelM_PPM, eMODaSynFileColumns.DelM_PPM)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_MH, eMODaSynFileColumns.MH)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Peptide, eMODaSynFileColumns.Peptide)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Protein, eMODaSynFileColumns.Protein)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Score, eMODaSynFileColumns.Score)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Probability, eMODaSynFileColumns.Probability)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Rank_Probability, eMODaSynFileColumns.Rank_Probability)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_Peptide_Position, eMODaSynFileColumns.Peptide_Position)
        lstColumnNames.Add(clsPHRPParserMODa.DATA_COLUMN_QValue, eMODaSynFileColumns.QValue)

        Try
            ' Initialize each entry in intColumnMapping to -1
            For intIndex = 0 To intColumnMapping.Length - 1
                intColumnMapping(intIndex) = -1
            Next

            strSplitLine = strLineIn.Split(ControlChars.Tab)
            For intIndex = 0 To strSplitLine.Length - 1
                If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
                    ' Recognized column name; update intColumnMapping
                    intColumnMapping(eResultFileColumn) = intIndex
                End If

            Next

        Catch ex As Exception
            SetErrorMessage("Error parsing header in MODa synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMODaSynFileEntry(
      strLineIn As String,
      objSearchResult As clsSearchResultsMODa,
      ByRef strErrorLog As String,
      intResultsProcessed As Integer,
      ByRef intColumnMapping() As Integer,
      ByRef strPeptideSequenceWithMods As String) As Boolean

        ' Parses an entry from the MODa Synopsis file

        Dim strSplitLine() As String = Nothing

        Try

            ' Reset objSearchResult
            objSearchResult.Clear()
            strPeptideSequenceWithMods = String.Empty

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length < 13 Then
                Return False
            End If

            With objSearchResult
                Dim strValue As String = Nothing
                If Not GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.ResultID), strValue) Then
                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                        strErrorLog &= "Error reading ResultID value from MODa Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                    End If
                    Exit Try
                End If

                .ResultID = Integer.Parse(strValue)

                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Scan), .Scan)
                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Charge), .Charge)

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Peptide), strPeptideSequenceWithMods) Then
                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                        strErrorLog &= "Error reading Peptide sequence value from MODa Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                    End If
                    Exit Try
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Protein), .ProteinName)
                .MultipleProteinCount = "0"

                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.DelM), .MODaComputedDelM)
                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.DelM_PPM), .MODaComputedDelMPPM)

                .PeptideDeltaMass = .MODaComputedDelM

                ' Note: .PeptideDeltaMass is stored in the MODa results file as "Observed_Mass - Theoretical_Mass"
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
                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Spectrum_Index), .Spectrum_Index)

                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.PrecursorMZ), .Precursor_mz)

                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.MH), .ParentIonMH)

                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Score), .MODaScore)
                GetColumnValue(strSplitLine, intColumnMapping(eMODaSynFileColumns.Probability), .Probability)

            End With

            Return True

        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MODa Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MODa Results in ParseMODaSynFileEntry" & ControlChars.NewLine
                End If
            End If
        End Try

        Return False

    End Function

    ''' <summary>
    ''' Main processing function
    ''' </summary>
    ''' <param name="strInputFilePath">MODa results file (Dataset_moda.id.txt)</param>
    ''' <param name="strOutputFolderPath">Output folder</param>
    ''' <param name="strParameterFilePath">Parameter file</param>
    ''' <returns>True if success, False if failure</returns>
    Public Overloads Overrides Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean

        Dim strBaseName As String = String.Empty
        Dim strSynOutputFilePath As String = String.Empty

        Dim lstMODaModInfo As List(Of clsModificationDefinition)
        Dim lstPepToProteinMapping As List(Of udtPepToProteinMappingType)

        Dim blnSuccess As Boolean

        If Not LoadParameterFileSettings(strParameterFilePath) Then
            SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, True)
            Return False
        End If

        Try
            If String.IsNullOrWhiteSpace(strInputFilePath) Then
                SetErrorMessage("Input file name is empty")
                SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath)
                Return False
            End If

            blnSuccess = ResetMassCorrectionTagsAndModificationDefinitions()
            If Not blnSuccess Then
                Return False
            End If

            MyBase.ResetProgress("Parsing " & Path.GetFileName(strInputFilePath))

            If Not CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                Return False
            End If

            Try
                ' Obtain the full path to the input file
                Dim fiInputFile = New FileInfo(strInputFilePath)

                lstMODaModInfo = New List(Of clsModificationDefinition)
                lstPepToProteinMapping = New List(Of udtPepToProteinMappingType)

                ' Load the MODa Parameter File to look for any static mods
                ExtractModInfoFromMODaParamFile(mSearchToolParameterFilePath, lstMODaModInfo)

                ' Resolve the mods in lstMODaModInfo with the ModDefs mods
                ResolveMODaModsWithModDefinitions(lstMODaModInfo)

                ' Define the base output filename using strInputFilePath
                strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath)

                ' Auto-replace "_moda.id" with "_moda"
                If strBaseName.ToLower().EndsWith("_moda.id") Then
                    strBaseName = strBaseName.Substring(0, strBaseName.Length - "_moda.id".Length) & "_moda"
                End If

                ' Load the MSG IndexToScanMap file (if it exists)
                LoadMGFIndexToScanMapFile(fiInputFile)

                ' Do not create a first-hits file for MODa results

                ' Create the synopsis output file
                MyBase.ResetProgress("Creating the SYN file", True)

                ' The synopsis file name will be of the form BasePath_moda_syn.txt
                strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath)

                ' Create the other PHRP-specific files
                MyBase.ResetProgress("Creating the PHRP files for " & Path.GetFileName(strSynOutputFilePath), True)

                ' Now parse the _syn.txt file that we just created to next create the other PHRP files
                blnSuccess = ParseMODaSynopsisFile(strSynOutputFilePath, strOutputFolderPath, lstPepToProteinMapping, False)

                ' This step is not necessary
                'If blnSuccess Then
                '	blnSuccess = AppendDelMPPMRefinedToSynFile(strSynOutputFilePath)
                'End If

                ' Remove all items from lstPepToProteinMapping to reduce memory overhead
                lstPepToProteinMapping.Clear()
                lstPepToProteinMapping.TrimExcess()

                If blnSuccess AndAlso mCreateProteinModsFile Then
                    blnSuccess = CreateProteinModsFileWork(strBaseName, fiInputFile, strSynOutputFilePath, strOutputFolderPath)
                End If

                If blnSuccess Then
                    MyBase.OperationComplete()
                End If

            Catch ex As Exception
                SetErrorMessage("Error in clsMODaResultsProcessor.ProcessFile (2):  " & ex.Message)
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
            End Try

        Catch ex As Exception
            SetErrorMessage("Error in ProcessFile (1):" & ex.Message)
            SetErrorCode(ePHRPErrorCodes.UnspecifiedError)
        End Try

        Return blnSuccess

    End Function

    Private Function CreateProteinModsFileWork(
       strBaseName As String,
       fiInputFile As FileInfo,
       strSynOutputFilePath As String,
       strOutputFolderPath As String) As Boolean

        Dim blnSuccess As Boolean
        Dim strMTSPepToProteinMapFilePath As String

        ' Create the MTSPepToProteinMap file

        strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS:=True)

        Dim lstSourcePHRPDataFiles = New List(Of String)

        If Not String.IsNullOrEmpty(strSynOutputFilePath) Then
            lstSourcePHRPDataFiles.Add(strSynOutputFilePath)
        End If

        If lstSourcePHRPDataFiles.Count = 0 Then
            SetErrorMessage("Cannot call CreatePepToProteinMapFile since lstSourcePHRPDataFiles is empty")
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        Else
            If File.Exists(strMTSPepToProteinMapFilePath) AndAlso mUseExistingMTSPepToProteinMapFile Then
                blnSuccess = True
            Else
                ' Auto-change mIgnorePeptideToProteinMapperErrors to True
                ' We only do this since a small number of peptides reported by MODa don't perfectly match the fasta file
                mIgnorePeptideToProteinMapperErrors = True
                blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
                If Not blnSuccess Then
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False")
                End If
            End If
        End If

        If blnSuccess Then
            ' If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
            MyBase.ValidatePHRPReaderSupportFiles(Path.Combine(fiInputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath)

            ' Create the Protein Mods file
            blnSuccess = MyBase.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MODa)
        End If

        If Not blnSuccess Then
            ' Do not treat this as a fatal error
            blnSuccess = True
        End If
        Return blnSuccess
    End Function

    Private Sub ResolveMODaModsWithModDefinitions(ByRef lstMODaModInfo As List(Of clsModificationDefinition))

        Dim blnExistingModFound As Boolean
        Dim objModDef As clsModificationDefinition

        If Not lstMODaModInfo Is Nothing Then

            ' Call .LookupModificationDefinitionByMass for each entry in lstMODaModInfo
            For Each objModInfo As clsModificationDefinition In lstMODaModInfo
                If String.IsNullOrEmpty(objModInfo.TargetResidues) Then
                    objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, Nothing, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnExistingModFound, True, MODA_MASS_DIGITS_OF_PRECISION)
                Else
                    For Each chTargetResidue As Char In objModInfo.TargetResidues
                        objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, chTargetResidue, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnExistingModFound, True, MODA_MASS_DIGITS_OF_PRECISION)
                    Next
                End If
            Next

        End If

    End Sub

    Private Sub SortAndWriteFilteredSearchResults(
      swResultFile As StreamWriter,
      lstFilteredSearchResults As List(Of udtMODaSearchResultType),
      ByRef strErrorLog As String)

        Dim intIndex As Integer

        ' Sort udtFilteredSearchResults by descending probability, ascending scan, ascending charge, ascending peptide, and ascending protein
        lstFilteredSearchResults.Sort(New MODaSearchResultsComparerProbabilityScanChargePeptide)

        ' Compute FDR values then assign QValues
        ComputeQValues(lstFilteredSearchResults)

        For intIndex = 0 To lstFilteredSearchResults.Count - 1
            WriteSearchResultToFile(intIndex + 1, swResultFile, lstFilteredSearchResults(intIndex), strErrorLog)
        Next intIndex

    End Sub

    ''' <summary>
    ''' Compute FDR values then assign QValues
    ''' </summary>
    ''' <param name="lstSearchResults"></param>
    ''' <remarks>Assumes the data is sorted by descending probability using MODaSearchResultsComparerProbabilityScanChargePeptide</remarks>
    Private Sub ComputeQValues(lstSearchResults As List(Of udtMODaSearchResultType))

        Dim forwardPeptideCount = 0
        Dim reversePeptideCount = 0
        Dim intIndex = 0

        While intIndex < lstSearchResults.Count

            ' Check for entries with multiple proteins listed
            Dim intIndexEnd = intIndex
            Do While intIndexEnd + 1 < lstSearchResults.Count
                If lstSearchResults(intIndex).ScanNum = lstSearchResults(intIndexEnd + 1).ScanNum AndAlso
                   lstSearchResults(intIndex).ChargeNum = lstSearchResults(intIndexEnd + 1).ChargeNum AndAlso
                   lstSearchResults(intIndex).Peptide = lstSearchResults(intIndexEnd + 1).Peptide Then
                    intIndexEnd += 1
                Else
                    Exit Do
                End If
            Loop

            Dim isReverse = True

            ' Look for non-reverse proteins
            For intIndexCheck = intIndex To intIndexEnd
                If Not IsReversedProtein(lstSearchResults(intIndexCheck).Protein) Then
                    isReverse = False
                    Exit For
                End If
            Next

            If isReverse Then
                reversePeptideCount += 1
            Else
                forwardPeptideCount += 1
            End If

            Dim dblFDR As Double = 1

            If forwardPeptideCount > 0 Then
                dblFDR = reversePeptideCount / CDbl(forwardPeptideCount)
            End If

            ' Store the FDR values
            For intIndexStore = intIndex To intIndexEnd
                Dim udtResult = lstSearchResults(intIndexStore)
                udtResult.FDR = dblFDR

                lstSearchResults(intIndexStore) = udtResult
            Next

            intIndex = intIndexEnd + 1
        End While

        ' Now compute QValues
        ' We step through the list, from the worst scoring result to the best result
        ' The first QValue is the FDR of the final entry
        ' The next QValue is the minimum of (QValue, CurrentFDR)

        Dim dblQValue = lstSearchResults.Last().FDR
        If dblQValue > 1 Then dblQValue = 1

        For intIndex = lstSearchResults.Count - 1 To 0 Step -1
            Dim udtResult = lstSearchResults(intIndex)

            dblQValue = Math.Min(dblQValue, udtResult.FDR)
            udtResult.QValue = dblQValue

            lstSearchResults(intIndex) = udtResult
        Next

    End Sub

    Private Sub StoreSynMatches(
      lstSearchResults As List(Of udtMODaSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer,
      lstFilteredSearchResults As List(Of udtMODaSearchResultType))

        Dim intIndex As Integer

        AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure already sorted by scan, charge, and SpecProb; no need to re-sort

        ' Now store or write out the matches that pass the filters
        For intIndex = intStartIndex To intEndIndex
            If lstSearchResults(intIndex).ProbabilityNum >= mMODaMODPlusSynopsisFileProbabilityThreshold Then
                lstFilteredSearchResults.Add(lstSearchResults(intIndex))
            End If
        Next intIndex

    End Sub

    Private Sub WriteSynFHTFileHeader(
      swResultFile As StreamWriter,
      ByRef strErrorLog As String)

        ' Write out the header line for synopsis / first hits files
        Try
            Dim lstData As New List(Of String)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_ResultID)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Scan)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Spectrum_Index)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Charge)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_PrecursorMZ)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_DelM)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_DelM_PPM)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_MH)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Peptide)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Protein)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Score)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Probability)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Rank_Probability)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_Peptide_Position)
            lstData.Add(clsPHRPParserMODa.DATA_COLUMN_QValue)

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
    Private Sub WriteSearchResultToFile(
      intResultID As Integer,
       swResultFile As StreamWriter,
       udtSearchResult As udtMODaSearchResultType,
       ByRef strErrorLog As String)

        Try

            ' Primary Columns
            '
            ' MODa
            ' ResultID   Scan   Spectrum_Index   Charge   PrecursorMZ   DelM   DelM_PPM   MH   Peptide   Protein   Score   Probability   Rank_Probability   PeptidePosition      QValue

            Dim lstData As New List(Of String)
            lstData.Add(intResultID.ToString())
            lstData.Add(udtSearchResult.ScanNum.ToString)
            lstData.Add(udtSearchResult.SpectrumIndex)
            lstData.Add(udtSearchResult.Charge)
            lstData.Add(udtSearchResult.PrecursorMZ)
            lstData.Add(udtSearchResult.DelM)
            lstData.Add(udtSearchResult.DelM_PPM)
            lstData.Add(udtSearchResult.MH)
            lstData.Add(udtSearchResult.Peptide)
            lstData.Add(udtSearchResult.Protein)
            lstData.Add(udtSearchResult.Score)
            lstData.Add(udtSearchResult.Probability)
            lstData.Add(udtSearchResult.RankProbability.ToString())
            lstData.Add(udtSearchResult.PeptidePosition)
            lstData.Add(udtSearchResult.QValue.ToString("0.000000"))

            swResultFile.WriteLine(CollapseList(lstData))

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

#Region "IComparer Classes"

    Private Class MODaSearchResultsComparerScanChargeProbabilityPeptide
        Implements IComparer(Of udtMODaSearchResultType)

        Public Function Compare(x As udtMODaSearchResultType, y As udtMODaSearchResultType) As Integer Implements IComparer(Of udtMODaSearchResultType).Compare

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
                    ' Charge is the same; check ProbabilityNum
                    If x.ProbabilityNum < y.ProbabilityNum Then
                        Return 1
                    ElseIf x.ProbabilityNum > y.ProbabilityNum Then
                        Return -1
                    Else
                        ' Probability is the same; check peptide
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

    Private Class MODaSearchResultsComparerProbabilityScanChargePeptide
        Implements IComparer(Of udtMODaSearchResultType)

        Public Function Compare(x As udtMODaSearchResultType, y As udtMODaSearchResultType) As Integer Implements IComparer(Of udtMODaSearchResultType).Compare

            If x.ProbabilityNum < y.ProbabilityNum Then
                Return 1
            ElseIf x.ProbabilityNum > y.ProbabilityNum Then
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
