Option Strict On

' This class reads in an MSGF_DB results file (txt format) and creates 
' a tab-delimited text file with the data.  It will insert modification symbols
' into the peptide sequences for modified peptides.
'
' The modification definition information is determined from the MSGF-DB parameter file
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 8/12/2011
'
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports System.IO
Imports System.Runtime.InteropServices
Imports PHRPReader
Imports System.Text.RegularExpressions

Public Class clsMSGFDBResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "May 14, 2016"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MSGFDB_FILE As String = "_msgfdb"
    Public Const FILENAME_SUFFIX_MSGFPLUS_FILE As String = "_msgfplus"

    Public Const N_TERMINUS_SYMBOL_MSGFDB As String = "_."
    Public Const C_TERMINUS_SYMBOL_MSGFDB As String = "._"

    ' Filter passing peptides have MSGFDB_SpecEValue <= 0.0001 Or EValue <= DEFAULT_SYN_FILE_PVALUE_THRESHOLD
    ' These filters are also used by MSPathFinder
    Public Const DEFAULT_SYN_FILE_MSGF_SPECPROB_THRESHOLD As Single = 0.0001
    Public Const DEFAULT_SYN_FILE_PVALUE_THRESHOLD As Single = 0.95

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    ' Match mod masses (positive or negative) at start, e.g. 
    ' +57.021HWWTLTTDRINK         matches +57.021
    ' -57.021+42.011HWWTLTTDRINK  matches -57.021+42.011 (two separate mods)
    ' +42.011MDHTPQSQLK           matches +42.011
    Private Const MSGFDB_NTERMINAL_MOD_MASS_REGEX As String = "^([0-9\.\+\-]+)"

    ' Match mod masses (positive or negative) at end, e.g. 
    ' FAACPLTCE+14.0157VS+79.9663+14.0157   matches +79.9663+14.0157
    Private Const MSGFDB_CTERMINAL_MOD_MASS_REGEX As String = "([0-9\.\+\-]+)$"

    Private Const MSGFDB_MOD_MASS_REGEX As String = "([+-][0-9\.]+)"
    Private Const PROTEIN_AND_TERM_SYMBOLS_REGEX As String = "([^;]+)\(pre=(.),post=(.)\)"

    Private Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline Or RegexOptions.IgnoreCase

    ' These columns correspond to the tab-delimited file created directly by MSGF-DB
    Protected Const MSGFDBResultsFileColCount As Integer = 20
    Public Enum eMSGFDBResultsFileColumns As Integer
        SpectrumFile = 0
        SpecIndex = 1               ' SpecID in MSGF+
        Scan = 2
        FragMethod = 3
        PrecursorMZ = 4
        PMErrorDa = 5               ' Corresponds to PMError(Da)
        PMErrorPPM = 6              ' Corresponds to PMError(ppm)
        Charge = 7
        Peptide = 8
        Protein = 9
        DeNovoScore = 10
        MSGFScore = 11
        SpecProb_EValue = 12
        PValue_EValue = 13
        FDR_QValue = 14             ' Only present if searched using -tda 1
        PepFDR_PepQValue = 15       ' Only present if searched using -tda 1
        EFDR = 16                   ' Only present if did not search using -tda 1
        IMSScan = 17                ' Only present for MSGFDB_IMS results
        IMSDriftTime = 18           ' Only present for MSGFDB_IMS results
        IsotopeError = 19           ' Only reported by MSGF+
    End Enum

    ' These columns correspond to the Synopsis and First-Hits files created by this class
    Protected Const MSGFDBSynFileColCount As Integer = 23
    Public Enum eMSFDBSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        FragMethod = 2
        SpecIndex = 3
        Charge = 4
        PrecursorMZ = 5
        DelM = 6                            ' Precursor error, in Da; if the search used a tolerance less than 0.5 Da or less than 500 ppm, then this value is computed from the DelMPPM value
        DelMPPM = 7                         ' Precursor error, in ppm; corrected for isotope selection errors
        MH = 8                              ' Theoretical monoisotopic peptide mass (computed by PHRP)
        Peptide = 9                         ' This is the sequence with prefix and suffix residues and also with modification symbols
        Protein = 10                        ' Protein Name (remove description)
        NTT = 11                            ' Number of tryptic terminii
        DeNovoScore = 12
        MSGFScore = 13
        SpecProb_EValue = 14
        RankSpecProb = 15                   ' Rank 1 means lowest SpecProb, 2 means next higher score, etc. (ties get the same rank)
        PValue_EValue = 16
        FDR_QValue = 17                     ' Only present if searched using -tda 1
        PepFDR_PepQValue = 18               ' Only present if searched using -tda 1
        EFDR = 19                           ' Only present if did not search using -tda 1
        IMSScan = 20                        ' Only present for MSGFDB_IMS results
        IMSDriftTime = 21                   ' Only present for MSGFDB_IMS results
        IsotopeError = 22
    End Enum

    Protected Enum eFilteredOutputFileTypeConstants As Integer
        SynFile = 0
        FHTFile = 1
    End Enum
#End Region

#Region "Structures"
    Protected Structure udtMSGFDBSearchResultType

        Public SpectrumFileName As String
        Public SpecIndex As String
        Public Scan As String
        Public ScanNum As Integer
        Public FragMethod As String
        Public PrecursorMZ As String
        Public PMErrorDa As String              ' Corresponds to PMError(Da); MSGFDB stores this value as Observed - Theoretical
        Public PMErrorPPM As String             ' Corresponds to PMError(ppm); MSGFDB stores this value as Observed - Theoretical
        Public MH As String
        Public Charge As String
        Public ChargeNum As Short
        Public Peptide As String
        Public Protein As String
        Public NTT As String
        Public DeNovoScore As String
        Public MSGFScore As String
        Public SpecProb As String                   ' Smaller values are better scores (e.g. 1E-9 is better than 1E-6); holds SpecEValue for MSGF+
        Public SpecProbNum As Double
        Public PValue As String                     ' Smaller values are better scores (e.g. 1E-7 is better than 1E-3); holds EValue for MSGF+
        Public PValueNum As Double
        Public FDR As String                    ' Holds FDR when a target/decoy search was used; holds EFDR when a non-decoy search was used; holds QValue for MSGF+
        Public PepFDR As String                 ' Only used when target/decoy search was used; holds PepQValue for MSGF+
        Public RankSpecProb As Integer
        Public IMSScan As Integer
        Public IMSDriftTime As String
        Public IsotopeError As Integer          ' Only used by MSGF+

        Public Sub Clear()
            SpectrumFileName = String.Empty
            SpecIndex = String.Empty
            ScanNum = 0
            FragMethod = String.Empty
            PrecursorMZ = String.Empty
            PMErrorDa = String.Empty
            PMErrorPPM = String.Empty
            MH = String.Empty
            Charge = String.Empty
            ChargeNum = 0
            Peptide = String.Empty
            Protein = String.Empty
            NTT = String.Empty
            DeNovoScore = String.Empty
            MSGFScore = String.Empty
            SpecProb = String.Empty
            SpecProbNum = 0
            PValue = String.Empty
            PValueNum = 0
            FDR = String.Empty
            PepFDR = String.Empty
            RankSpecProb = 0
            IMSScan = 0
            IMSDriftTime = String.Empty
            IsotopeError = 0
        End Sub
    End Structure

    Protected Structure udtScanGroupInfoType
        Public ScanGroupID As Integer
        Public Charge As Short
        Public Scan As Integer
    End Structure

    Protected Structure udtTerminusCharsType
        Public NTerm As Char
        Public CTerm As Char
    End Structure

    Protected Structure udtParentMassToleranceType
        ' Given a tolerance of 20ppm, we would have ToleranceLeft=20, ToleranceRight=20, and ToleranceIsPPM=True
        ' Given a tolerance of 0.5Da,2.5Da, we would have ToleranceLeft=0.5, ToleranceRight=2.5, and ToleranceIsPPM=False
        Public ToleranceLeft As Double
        Public ToleranceRight As Double
        Public IsPPM As Boolean

        Public Sub Clear()
            ToleranceLeft = 0
            ToleranceRight = 0
            IsPPM = False
        End Sub
    End Structure

#End Region

#Region "Classwide Variables"
    Protected mPeptideCleavageStateCalculator As clsPeptideCleavageStateCalculator

    Protected mParentMassToleranceInfo As udtParentMassToleranceType

    Protected mPrecursorMassErrorWarningCount As Integer

#End Region

    ''' <summary>
    ''' Step through .PeptideSequenceWithMods
    ''' For each residue, check if a static mod is defined that affects that residue
    ''' For each mod symbol, determine the modification and add to objSearchResult
    ''' </summary>
    ''' <param name="objSearchResult"></param>
    ''' <param name="blnUpdateModOccurrenceCounts"></param>
    ''' <remarks></remarks>
    Private Sub AddDynamicAndStaticResidueMods(objSearchResult As clsSearchResultsMSGFDB, blnUpdateModOccurrenceCounts As Boolean)
        Dim intIndex As Integer, intModIndex As Integer
        Dim chChar As Char
        Dim objModificationDefinition As clsModificationDefinition

        Dim strSequence As String
        Dim chMostRecentLetter As Char
        Dim intResidueLocInPeptide As Integer
        Dim blnSuccess As Boolean

        chMostRecentLetter = "-"c
        intResidueLocInPeptide = 0

        With objSearchResult
            strSequence = .PeptideSequenceWithMods

            For intIndex = 0 To strSequence.Length - 1
                chChar = strSequence.Chars(intIndex)

                If IsLetterAtoZ(chChar) Then
                    chMostRecentLetter = chChar
                    intResidueLocInPeptide += 1

                    For intModIndex = 0 To mPeptideMods.ModificationCount - 1
                        If mPeptideMods.GetModificationTypeByIndex(intModIndex) = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                            objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)

                            If objModificationDefinition.TargetResiduesContain(chChar) Then
                                ' Match found; add this modification
                                .SearchResultAddModification(objModificationDefinition, chChar, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                            End If
                        End If
                    Next intModIndex
                ElseIf IsLetterAtoZ(chMostRecentLetter) Then
                    blnSuccess = .SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                    If Not blnSuccess Then
                        Dim strErrorMessage As String = .ErrorMessage
                        If String.IsNullOrEmpty(strErrorMessage) Then
                            strErrorMessage = "SearchResultAddDynamicModification returned false for symbol " & chChar
                        End If
                        SetErrorMessage(strErrorMessage & "; ResultID = " & .ResultID)
                    End If
                Else
                    ' We found a modification symbol but chMostRecentLetter is not a letter
                    ' Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                End If

            Next intIndex
        End With
    End Sub

    ''' <summary>
    ''' Adds or updates the prefix and suffix residues to the peptide, as defined in kvProteinInfo
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="kvProteinInfo"></param>
    ''' <returns>Peptide sequence with N-terminal and C-Terminal residues</returns>
    ''' <remarks></remarks>
    Protected Function AddUpdatePrefixAndSuffixResidues(strPeptide As String, kvProteinInfo As KeyValuePair(Of String, udtTerminusCharsType)) As String

        If strPeptide.IndexOf("."c) < 0 Then
            Return kvProteinInfo.Value.NTerm & "." & strPeptide & "." & kvProteinInfo.Value.CTerm
        Else
            Dim strPeptideNew As String

            If strPeptide.Length >= 2 Then
                If strPeptide.Chars(1) = "."c Then
                    ' Peptide already has the N-terminal residue
                    ' Replace it using kvProteinInfo
                    strPeptideNew = kvProteinInfo.Value.NTerm & "." & strPeptide.Substring(2)

                ElseIf strPeptide.Chars(0) = "."c Then
                    strPeptideNew = kvProteinInfo.Value.NTerm & strPeptide
                Else
                    strPeptideNew = kvProteinInfo.Value.NTerm & "." & strPeptide
                End If
            Else
                strPeptideNew = String.Copy(strPeptide)
            End If

            If strPeptideNew.Length >= 4 Then
                If strPeptideNew.Chars(strPeptideNew.Length - 2) = "."c Then
                    ' Peptide already has the C-terminal residue
                    ' Replace it using kvProteinInfo
                    strPeptideNew = strPeptideNew.Substring(0, strPeptideNew.Length - 2) & "." & kvProteinInfo.Value.CTerm

                ElseIf strPeptideNew.Chars(strPeptideNew.Length - 1) = "."c Then
                    strPeptideNew = strPeptideNew & kvProteinInfo.Value.CTerm
                Else
                    strPeptideNew = strPeptideNew & "." & kvProteinInfo.Value.CTerm
                End If
            End If

            Return strPeptideNew
        End If

    End Function

    Protected Sub AppendToScanGroupDetails(
      lstScanGroupDetails As List(Of udtScanGroupInfoType),
      htScanGroupCombo As Dictionary(Of String, Boolean),
      udtScanGroupInfo As udtScanGroupInfoType,
      ByRef intCurrentScanGroupID As Integer,
      ByRef intNextScanGroupID As Integer)

        Dim strChargeScanComboText As String

        strChargeScanComboText = udtScanGroupInfo.Charge.ToString & "_" & udtScanGroupInfo.Scan.ToString()

        If Not htScanGroupCombo.ContainsKey(strChargeScanComboText) Then
            If intCurrentScanGroupID < 0 Then
                intCurrentScanGroupID = intNextScanGroupID
                intNextScanGroupID += 1
            End If

            udtScanGroupInfo.ScanGroupID = intCurrentScanGroupID

            lstScanGroupDetails.Add(udtScanGroupInfo)
            htScanGroupCombo.Add(strChargeScanComboText, True)
        End If

    End Sub

    Private Sub AppendToSearchResults(
       lstSearchResults As List(Of udtMSGFDBSearchResultType),
       udtSearchResult As udtMSGFDBSearchResultType,
       lstProteinInfo As Dictionary(Of String, udtTerminusCharsType))

        If lstProteinInfo.Count = 0 Then
            lstSearchResults.Add(udtSearchResult)
        Else
            For Each kvEntry As KeyValuePair(Of String, udtTerminusCharsType) In lstProteinInfo
                udtSearchResult.Protein = kvEntry.Key
                udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, kvEntry)

                lstSearchResults.Add(udtSearchResult)
            Next
        End If

    End Sub

    ''' <summary>
    ''' Ranks each entry (assumes all of the data is from the same scan)
    ''' </summary>
    ''' <param name="lstSearchResults">Search results</param>
    ''' <param name="intStartIndex">Start index for data in this scan</param>
    ''' <param name="intEndIndex">End index for data in this scan</param>
    ''' <remarks></remarks>
    Private Sub AssignRankAndDeltaNormValues(
      lstSearchResults As IList(Of udtMSGFDBSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer)

        ' Prior to September 2014 ranks were assign per charge state per scan; 
        ' Ranks are now assigned per scan (across all charge states)

        If intStartIndex = intEndIndex Then
            ' Only one result
            Dim currentResult = lstSearchResults(intStartIndex)
            currentResult.RankSpecProb = 1
            lstSearchResults(intStartIndex) = currentResult
            Return
        End If

        ' Duplicate a portion of udtSearchResults so that we can sort by ascending Spectral Probability

        Dim dctResultsSubset = New Dictionary(Of Integer, udtMSGFDBSearchResultType)
        For intIndex = intStartIndex To intEndIndex
            dctResultsSubset.Add(intIndex, lstSearchResults(intIndex))
        Next

        Dim lstResultsBySpecProb = (From item In dctResultsSubset Select item Order By item.Value.SpecProbNum).ToList()

        Dim dblLastValue As Double
        Dim intCurrentRank As Integer = -1

        For Each entry In lstResultsBySpecProb
            Dim currentResult = lstSearchResults(entry.Key)

            If intCurrentRank < 0 Then
                dblLastValue = currentResult.SpecProbNum
                intCurrentRank = 1
            Else
                If Math.Abs(currentResult.SpecProbNum - dblLastValue) > Double.Epsilon Then
                    dblLastValue = currentResult.SpecProbNum
                    intCurrentRank += 1
                End If
            End If

            currentResult.RankSpecProb = intCurrentRank

            ' Because this is a list of structs, we have to copy currentResult back into the current position in lstSearchResults
            lstSearchResults(entry.Key) = currentResult
        Next

    End Sub

    Private Function AddModificationsAndComputeMass(objSearchResult As clsSearchResultsMSGFDB, blnUpdateModOccurrenceCounts As Boolean) As Boolean
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

    Protected Function ComputeCleaveageState(strSequenceWithMods As String) As Short

        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Dim eCleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants

        Dim strCleanSequence As String

        ' Remove any non-letter characters when before .ComputeCleavageState()
        strCleanSequence = GetCleanSequence(strSequenceWithMods, strPrefix, strSuffix)

        eCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(strCleanSequence, strPrefix, strSuffix)

        Return CShort(eCleavageState)

    End Function

    ''' <summary>
    ''' This function should only be called when column PMError(Da) is present (and PMError(ppm) is not present)
    ''' </summary>
    ''' <param name="dblPrecursorErrorDa">Mass error (Observed - theoretical)</param>
    ''' <param name="dblPrecursorMZ"></param>
    ''' <param name="intCharge"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ComputeDelMCorrectedPPM(
      dblPrecursorErrorDa As Double,
      dblPrecursorMZ As Double,
      intCharge As Integer,
      dblPeptideMonoisotopicMass As Double,
      blnAdjustPrecursorMassForC13 As Boolean) As Double

        Dim dblPeptideDeltaMassCorrectedPpm As Double

        Dim dblPrecursorMonoMass As Double

        ' Compute the original value for the precursor monoisotopic mass
        dblPrecursorMonoMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, intCharge, 0)

        dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, dblPeptideMonoisotopicMass)

        Return dblPeptideDeltaMassCorrectedPpm

    End Function

    Protected Function ComputePeptideMass(strPeptide As String, dblTotalModMass As Double) As Double

        Dim strCleanSequence As String
        Dim dblMass As Double

        strCleanSequence = GetCleanSequence(strPeptide)

        dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence)
        dblMass += dblTotalModMass

        Return dblMass

    End Function

    Protected Overrides Function ConstructPepToProteinMapFilePath(strInputFilePath As String, strOutputFolderPath As String, MTS As Boolean) As String

        Dim strPepToProteinMapFilePath As String = Path.GetFileNameWithoutExtension(strInputFilePath)

        If strPepToProteinMapFilePath.ToLower().EndsWith("_msgfdb_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_msgfdb_fht") Then
            ' Remove _syn or _fht
            strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4)
        End If

        Return MyBase.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS)

    End Function

    ''' <summary>
    ''' Parses the digits in strModDigits to convert them to one or more modification symbols
    ''' </summary>
    ''' <param name="strModDigits">Example: +57.021 or +79.9663+14.0157 or -18.0106</param>
    ''' <param name="strModSymbols"></param>
    ''' <returns>True if success; false if a problem</returns>
    ''' <remarks></remarks>
    Protected Function ConvertMGSFModMassesToSymbols(
      currentResidue As String,
      strModDigits As String,
      <Out()> ByRef strModSymbols As String,
      lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      blnNterminalMod As Boolean,
      blnPossibleCTerminalMod As Boolean,
      <Out()> ByRef dblModMassFound As Double,
      <Out()> ByRef blnIsStaticMod As Boolean) As Boolean

        Static reModMassRegEx As New Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS)

        Dim reMatches As MatchCollection

        Dim blnMatchFound As Boolean
        Dim blnTestMod As Boolean

        Dim intBestMatchIndex As Integer
        Dim dblBestMassDiff As Double
        Dim intModSymbolsFound = 0
        Dim chSymbolBestMatch As Char = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
        Dim residuesBestBatch = String.Empty

        strModSymbols = String.Empty
        dblModMassFound = 0
        blnIsStaticMod = False

        reMatches = reModMassRegEx.Matches(strModDigits)

        For Each reMatch As Match In reMatches
            Dim strModMass As String
            strModMass = reMatch.Value

            ' Convert strModMass to a mass value
            Dim dblModMass As Double
            dblModMass = Double.Parse(strModMass)
            dblModMassFound += dblModMass

            blnMatchFound = False
            Do
                intBestMatchIndex = -1

                ' Step through the known modifications to find the closest match
                For intIndex = 0 To lstMSGFDBModInfo.Count - 1

                    blnTestMod = True

                    If blnNterminalMod Then
                        ' Only test N-terminal mods
                        If Not (
                           lstMSGFDBModInfo(intIndex).ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynNTermPeptide OrElse
                           lstMSGFDBModInfo(intIndex).ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynNTermProtein) Then
                            blnTestMod = False
                        End If

                    ElseIf Not blnPossibleCTerminalMod Then
                        ' Skip C-terminal mods since we're not at the C-terminus
                        If lstMSGFDBModInfo(intIndex).ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynCTermPeptide OrElse
                           lstMSGFDBModInfo(intIndex).ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynCTermProtein Then
                            blnTestMod = False
                        End If

                    End If

                    If blnTestMod Then
                        Dim dblCandidateMassDiff = Math.Abs(lstMSGFDBModInfo(intIndex).ModMassVal - dblModMass)
                        If dblCandidateMassDiff < 0.25 Then
                            ' Possible match found
                            Dim updateCandidate = False

                            If intBestMatchIndex < 0 OrElse
                               dblCandidateMassDiff < dblBestMassDiff Then
                                updateCandidate = True

                            ElseIf (Math.Abs(dblCandidateMassDiff - dblBestMassDiff) < Single.Epsilon AndAlso
                                    chSymbolBestMatch = "-"c AndAlso
                                    lstMSGFDBModInfo(intIndex).ModSymbol <> "-"c) Then

                                ' Masses are the same, but the existing candidate is a static mod

                                ' If this new one is a dynamic mod, switch to it, but only if residuesBestBatch does not contain the residue while
                                ' the residues for the new candidate does contain the residue

                                If Not residuesBestBatch.Contains(currentResidue) AndAlso
                                   lstMSGFDBModInfo(intIndex).Residues.Contains(currentResidue) Then
                                    updateCandidate = True
                                End If
                            End If

                            If updateCandidate Then
                                intBestMatchIndex = intIndex
                                dblBestMassDiff = dblCandidateMassDiff
                                chSymbolBestMatch = lstMSGFDBModInfo(intIndex).ModSymbol
                                residuesBestBatch = String.Copy(lstMSGFDBModInfo(intIndex).Residues)
                                blnMatchFound = True
                            End If
                        End If

                    End If

                Next

                If Not blnMatchFound AndAlso blnNterminalMod Then
                    ' Set this to false, then search again
                    blnNterminalMod = False
                ElseIf Not blnMatchFound AndAlso Not blnPossibleCTerminalMod Then
                    ' Set this to true, then search again
                    blnPossibleCTerminalMod = True
                Else
                    Exit Do
                End If

            Loop While Not blnMatchFound

            If blnMatchFound Then
                ' Match found; use the mod symbol
                strModSymbols &= lstMSGFDBModInfo(intBestMatchIndex).ModSymbol
                intModSymbolsFound += 1

                If lstMSGFDBModInfo(intBestMatchIndex).ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.StaticMod Then
                    blnIsStaticMod = True
                End If
            Else
                ' Match not found; use the mass value
                strModSymbols &= strModMass
            End If

        Next reMatch

        If intModSymbolsFound > 0 Then
            Return True
        Else
            Return False
        End If


    End Function

    ''' <summary>
    ''' This routine creates a first hits file or synopsis file from the output from MSGF-DB
    ''' The synopsis file includes every result with a p-value below a set threshold or a SpecProb below a certain threshold
    ''' The first-hits file includes the results with the lowest SpecProb (for each scan and charge)
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFilePath"></param>
    ''' <param name="strScanGroupFilePath"></param>
    ''' <param name="lstMSGFDBModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
    ''' <param name="blnMSGFPlus">Output parameter: this function will set this to True if we're processing MSGF+ results</param>
    ''' <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CreateFHTorSYNResultsFile(
      strInputFilePath As String,
      strOutputFilePath As String,
      strScanGroupFilePath As String,
      lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      <Out()> ByRef blnMSGFPlus As Boolean,
      lstSpecIdToIndex As Dictionary(Of String, Integer),
      eFilteredOutputFileType As eFilteredOutputFileTypeConstants) As Boolean

        Dim strLineIn As String

        Dim lstSearchResultsCurrentScan As New List(Of udtMSGFDBSearchResultType)
        Dim lstSearchResultsUnfiltered As New List(Of udtMSGFDBSearchResultType)

        Dim sngPercentComplete As Single

        Dim blnHeaderParsed As Boolean
        Dim blnIncludeFDRandPepFDR = False
        Dim blnIncludeEFDR = False
        Dim blnIncludeIMSFields = False

        Dim intColumnMapping() As Integer = Nothing

        Dim intNextScanGroupID As Integer
        Dim lstScanGroupDetails As List(Of udtScanGroupInfoType)
        Dim htScanGroupCombo As Dictionary(Of String, Boolean)

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean

        Dim strErrorLog As String

        blnMSGFPlus = False

        Try
            lstScanGroupDetails = New List(Of udtScanGroupInfoType)
            htScanGroupCombo = New Dictionary(Of String, Boolean)

            mPrecursorMassErrorWarningCount = 0

            ' Look for custom amino acids
            Dim lstCustomAA = (From item In lstMSGFDBModInfo Where item.ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA Select item).ToList()
            For Each customAADef In lstCustomAA
                Dim aminoAcidSymbol = customAADef.Residues(0)
                Dim empiricalFormula = customAADef.ModMass
                Dim aminoAcidMass = customAADef.ModMassVal

                Try
                    Dim elementalComposition = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormula)
                    Dim atomCounts = clsPeptideMassCalculator.ConvertElementalCompositionToAtomCounts(elementalComposition)

                    mPeptideSeqMassCalculator.SetAminoAcidMass(aminoAcidSymbol, aminoAcidMass)
                    mPeptideSeqMassCalculator.SetAminoAcidAtomCounts(aminoAcidSymbol, atomCounts)

                Catch ex As Exception
                    ReportError(ex.Message)
                End Try

            Next

            Try
                ' Open the input file and parse it
                ' Initialize the stream reader and the stream Text writer
                Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)),
                      swResultFile = New StreamWriter(New FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                    strErrorLog = String.Empty
                    blnHeaderParsed = False

                    intNextScanGroupID = 1
                    lstScanGroupDetails.Clear()
                    htScanGroupCombo.Clear()

                    ' Initialize the array that will hold all of the records that will ultimately be written out to disk
                    Dim lstFilteredSearchResults = New List(Of udtMSGFDBSearchResultType)

                    ' Parse the input file
                    Do While Not srDataFile.EndOfStream And Not MyBase.AbortProcessing
                        strLineIn = srDataFile.ReadLine()
                        If String.IsNullOrWhiteSpace(strLineIn) Then
                            Continue Do
                        End If

                        If Not blnHeaderParsed Then
                            blnSuccess = ParseMSGFDBResultsFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True

                            blnIncludeFDRandPepFDR = False
                            blnIncludeEFDR = False

                            If intColumnMapping(eMSGFDBResultsFileColumns.FDR_QValue) >= 0 OrElse intColumnMapping(eMSGFDBResultsFileColumns.PepFDR_PepQValue) >= 0 Then
                                blnIncludeFDRandPepFDR = True
                            ElseIf intColumnMapping(eMSGFDBResultsFileColumns.EFDR) >= 0 Then
                                blnIncludeEFDR = True
                            End If

                            If intColumnMapping(eMSGFDBResultsFileColumns.IMSDriftTime) >= 0 Then
                                blnIncludeIMSFields = True
                            End If

                            If intColumnMapping(eMSGFDBResultsFileColumns.IsotopeError) >= 0 Then
                                blnMSGFPlus = True
                            End If

                            ' Write the header line
                            WriteSynFHTFileHeader(swResultFile, strErrorLog, blnIncludeFDRandPepFDR, blnIncludeEFDR, blnIncludeIMSFields, blnMSGFPlus)
                            Continue Do
                        End If

                        blnValidSearchResult = ParseMSGFDBResultsFileEntry(
                            strLineIn, blnMSGFPlus, lstMSGFDBModInfo,
                            lstSearchResultsCurrentScan, strErrorLog,
                            intColumnMapping, intNextScanGroupID, lstScanGroupDetails, htScanGroupCombo, lstSpecIdToIndex)

                        If blnValidSearchResult Then
                            ExpandListIfRequired(lstSearchResultsUnfiltered, lstSearchResultsCurrentScan.Count)
                            lstSearchResultsUnfiltered.AddRange(lstSearchResultsCurrentScan)
                        End If

                        ' Update the progress
                        sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
                        If mCreateProteinModsFile Then
                            sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
                        End If
                        UpdateProgress(sngPercentComplete)

                    Loop

                    lstSearchResultsUnfiltered.TrimExcess()

                    ' Sort the SearchResults by scan, charge, and ascending SpecProb
                    lstSearchResultsUnfiltered.Sort(New MSGFDBSearchResultsComparerScanChargeSpecProbPeptide)

                    ' Now filter the data

                    ' Initialize variables
                    Dim intStartIndex = 0
                    Dim intEndIndex As Integer

                    Do While intStartIndex < lstSearchResultsUnfiltered.Count
                        intEndIndex = intStartIndex
                        Do While intEndIndex + 1 < lstSearchResultsUnfiltered.Count AndAlso
                            lstSearchResultsUnfiltered(intEndIndex + 1).ScanNum = lstSearchResultsUnfiltered(intStartIndex).ScanNum
                            intEndIndex += 1
                        Loop

                        ' Store the results for this scan
                        If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
                            StoreSynMatches(lstSearchResultsUnfiltered, intStartIndex, intEndIndex, lstFilteredSearchResults)
                        Else
                            StoreTopFHTMatch(lstSearchResultsUnfiltered, intStartIndex, intEndIndex, lstFilteredSearchResults)
                        End If

                        intStartIndex = intEndIndex + 1
                    Loop

                    ' Sort the data in udtFilteredSearchResults then write out to disk
                    SortAndWriteFilteredSearchResults(swResultFile, lstFilteredSearchResults, strErrorLog, blnIncludeFDRandPepFDR, blnIncludeEFDR, blnIncludeIMSFields, blnMSGFPlus)
                End Using

                ' Write out the scan group info
                If Not String.IsNullOrEmpty(strScanGroupFilePath) Then
                    StoreScanGroupInfo(strScanGroupFilePath, lstScanGroupDetails)
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
            End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Extracts mod info from either a MSGF-DB param file or from a MSGFDB_Mods.txt file
    ''' </summary>
    ''' <param name="strMSGFDBParamFilePath"></param>
    ''' <param name="lstModInfo"></param>
    ''' <returns>True if success; false if a problem</returns>
    ''' <remarks></remarks>
    Private Function ExtractModInfoFromParamFile(
       strMSGFDBParamFilePath As String,
       <Out()> ByRef lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Boolean

        Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor("MSGF+")

        AddHandler modFileProcessor.ErrorOccurred, AddressOf ModExtractorErrorHandler
        AddHandler modFileProcessor.WarningMessageEvent, AddressOf ModExtractorWarningHandler

        ' Note that this call will initialize lstModInfo
        Dim success = modFileProcessor.ExtractModInfoFromParamFile(strMSGFDBParamFilePath, lstModInfo)

        If Not success OrElse mErrorCode <> ePHRPErrorCodes.NoError Then
            If mErrorCode = ePHRPErrorCodes.NoError Then
                SetErrorMessage("Unknown error extracting the modification definitions from the MSGF+ parameter file")
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            End If
            Return False
        End If

        modFileProcessor.ResolveMSGFDBModsWithModDefinitions(lstModInfo, mPeptideMods)

        Return True

    End Function

    ''' <summary>
    ''' Extracts parent mass tolerance from a MSGF-DB parameter file
    ''' </summary>
    ''' <param name="strMSGFDBParamFilePath"></param>	
    ''' <returns>Parent mass tolerance info.  Tolerances will be 0 if an error occurs</returns>
    ''' <remarks></remarks>
    Private Function ExtractParentMassToleranceFromParamFile(strMSGFDBParamFilePath As String) As udtParentMassToleranceType

        Const PM_TOLERANCE_TAG = "PMTolerance"

        Dim strLineIn As String
        Dim strSplitLine As String()

        Dim udtParentMassToleranceInfo As udtParentMassToleranceType

        Try
            udtParentMassToleranceInfo.Clear()

            If String.IsNullOrEmpty(strMSGFDBParamFilePath) Then
                SetErrorMessage("MSGFDB Parameter File name not defined; unable to extract parent mass tolerance info")
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                Return udtParentMassToleranceInfo
            End If

            ' Read the contents of the parameter file
            Using srInFile = New StreamReader(New FileStream(strMSGFDBParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                Do While Not srInFile.EndOfStream
                    strLineIn = srInFile.ReadLine().Trim()

                    If String.IsNullOrWhiteSpace(strLineIn) Then Continue Do

                    If strLineIn.StartsWith("#"c) Then
                        ' Comment line; skip it
                        Continue Do
                    End If

                    If Not strLineIn.ToLower.StartsWith(PM_TOLERANCE_TAG.ToLower()) Then
                        Continue Do
                    End If

                    Dim kvSetting As KeyValuePair(Of String, String)
                    kvSetting = clsPHRPParser.ParseKeyValueSetting(strLineIn, "="c, "#")

                    If String.IsNullOrEmpty(kvSetting.Value) Then
                        Exit Do
                    End If

                    ' Parent ion tolerance line found

                    ' Split the line on commas
                    strSplitLine = kvSetting.Value.Split(","c)

                    Dim dblTolerance As Double
                    Dim blnIsPPM As Boolean

                    If strSplitLine.Count = 1 Then

                        If ParseParentMassTolerance(strSplitLine(0), dblTolerance, blnIsPPM) Then
                            udtParentMassToleranceInfo.ToleranceLeft = dblTolerance
                            udtParentMassToleranceInfo.ToleranceRight = dblTolerance
                            udtParentMassToleranceInfo.IsPPM = blnIsPPM
                        End If

                    ElseIf strSplitLine.Count > 1 Then
                        If ParseParentMassTolerance(strSplitLine(0), dblTolerance, blnIsPPM) Then
                            udtParentMassToleranceInfo.ToleranceLeft = dblTolerance
                            udtParentMassToleranceInfo.IsPPM = blnIsPPM

                            If ParseParentMassTolerance(strSplitLine(1), dblTolerance, blnIsPPM) Then
                                udtParentMassToleranceInfo.ToleranceRight = dblTolerance
                            Else
                                udtParentMassToleranceInfo.ToleranceRight = udtParentMassToleranceInfo.ToleranceLeft
                            End If

                        End If
                    End If

                    Exit Do

                Loop

            End Using

            Console.WriteLine()

        Catch ex As Exception
            SetErrorMessage("Error reading ParentMass tolerance from the MSGF-DB parameter file (" & Path.GetFileName(strMSGFDBParamFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
        End Try

        Return udtParentMassToleranceInfo

    End Function

    Private Sub InitializeLocalVariables()

        ' Initialize mPeptideCleavageStateCalculator
        If mPeptideCleavageStateCalculator Is Nothing Then
            mPeptideCleavageStateCalculator = New clsPeptideCleavageStateCalculator
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin)
        End If

    End Sub

    ''' <summary>
    ''' Load the PeptideToProteinMap information; in addition, creates the _msgfdb_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
    ''' </summary>
    ''' <param name="strPepToProteinMapFilePath"></param>
    ''' <param name="strOutputFolderPath"></param>
    ''' <param name="lstMSGFDBModInfo"></param>
    ''' <param name="blnMSGFPlus">Should be set to True if processing MSGF+ results</param>
    ''' <param name="lstPepToProteinMapping"></param>
    ''' <param name="strMTSPepToProteinMapFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function LoadPeptideToProteinMapInfoMSGFDB(
      strPepToProteinMapFilePath As String,
      strOutputFolderPath As String,
      lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      blnMSGFPlus As Boolean,
      lstPepToProteinMapping As List(Of udtPepToProteinMappingType),
      <Out()> ByRef strMTSPepToProteinMapFilePath As String) As Boolean

        Dim strHeaderLine As String = String.Empty
        Dim strMTSCompatiblePeptide As String

        ' Not used by this function but required for the call to ReplaceMSGFModTextWithSymbol
        Dim dblTotalModMass As Double = 0

        Dim blnSuccess As Boolean

        strMTSPepToProteinMapFilePath = String.Empty

        Try

            If String.IsNullOrWhiteSpace(strPepToProteinMapFilePath) Then
                Console.WriteLine()
                ReportWarning("PepToProteinMap file is not defined")
                Return False
            ElseIf Not File.Exists(strPepToProteinMapFilePath) Then
                Console.WriteLine()
                ReportWarning("PepToProteinMap file does not exist: " & strPepToProteinMapFilePath)
                Return False
            End If

            ' Initialize lstPepToProteinMapping
            lstPepToProteinMapping.Clear()

            ' Read the data in strProteinToPeptideMappingFilePath
            blnSuccess = LoadPeptideToProteinMapInfo(strPepToProteinMapFilePath, lstPepToProteinMapping, strHeaderLine)

            If blnSuccess Then
                strMTSPepToProteinMapFilePath = Path.Combine(strOutputFolderPath, Path.GetFileNameWithoutExtension(strPepToProteinMapFilePath) & "MTS.txt")

                Using swOutFile = New StreamWriter(New FileStream(strMTSPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
                    If Not String.IsNullOrEmpty(strHeaderLine) Then
                        ' Header line
                        swOutFile.WriteLine(strHeaderLine)
                    End If

                    For intIndex = 0 To lstPepToProteinMapping.Count - 1
                        ' Replace any mod text names in the peptide sequence with the appropriate mod symbols
                        ' In addition, replace the * terminus symbols with dashes
                        strMTSCompatiblePeptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(lstPepToProteinMapping(intIndex).Peptide), lstMSGFDBModInfo, blnMSGFPlus, dblTotalModMass)

                        If lstPepToProteinMapping(intIndex).Peptide <> strMTSCompatiblePeptide Then
                            UpdatePepToProteinMapPeptide(lstPepToProteinMapping, intIndex, strMTSCompatiblePeptide)
                        End If

                        swOutFile.WriteLine(
                          lstPepToProteinMapping(intIndex).Peptide & ControlChars.Tab &
                          lstPepToProteinMapping(intIndex).Protein & ControlChars.Tab &
                          lstPepToProteinMapping(intIndex).ResidueStart & ControlChars.Tab &
                          lstPepToProteinMapping(intIndex).ResidueEnd)

                    Next

                End Using

            End If

        Catch ex As Exception
            SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" & Path.GetFileName(strMTSPepToProteinMapFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Sub ModExtractorErrorHandler(errMsg As String)
        SetErrorMessage(errMsg)
        SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
    End Sub

    Private Sub ModExtractorWarningHandler(warningMsg As String)
        ReportWarning(warningMsg)
    End Sub

    Protected Function ParseMSGFDBSynopsisFile(
      strInputFilePath As String,
      strOutputFolderPath As String,
      lstPepToProteinMapping As List(Of udtPepToProteinMappingType),
      blnResetMassCorrectionTagsAndModificationDefinitions As Boolean) As Boolean

        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim strPreviousSpecProb As String

        ' Note that MSGF-DB synopsis files are normally sorted on SpecProb value, ascending
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
        '  we will keep track of the scan, charge, and peptide information parsed for each unique SpecProb encountered
        ' Although this was a possiblity with Inspect, it likely never occurs for MSGF-DB
        '  But, we'll keep the check in place just in case

        Dim htPeptidesFoundForSpecProbLevel As Hashtable

        Dim strKey As String

        Dim strLineIn As String
        Dim strModificationSummaryFilePath As String

        Dim objSearchResult As clsSearchResultsMSGFDB

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

        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            ' Initialize objSearchResult
            objSearchResult = New clsSearchResultsMSGFDB(mPeptideMods, mPeptideSeqMassCalculator)

            ' Initialize htPeptidesFoundForSpecProbLevel
            htPeptidesFoundForSpecProbLevel = New Hashtable
            strPreviousSpecProb = String.Empty

            ' Assure that lstPepToProteinMapping is sorted on peptide
            If lstPepToProteinMapping.Count > 1 Then
                lstPepToProteinMapping.Sort(New PepToProteinMappingComparer)
            End If

            Try
                objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(mEnzymeMatchSpec, mPeptideNTerminusMassChange, mPeptideCTerminusMassChange)

                Dim strErrorLog As String = String.Empty

                ' Open the input file and parse it
                ' Initialize the stream reader
                Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                    intResultsProcessed = 0
                    blnHeaderParsed = False

                    ' Create the output files
                    Dim strBaseOutputFilePath As String = Path.Combine(strOutputFolderPath, Path.GetFileName(strInputFilePath))
                    blnSuccess = MyBase.InitializeSequenceOutputFiles(strBaseOutputFilePath)

                    ' Parse the input file
                    Do While Not srDataFile.EndOfStream And Not MyBase.AbortProcessing

                        strLineIn = srDataFile.ReadLine()
                        If String.IsNullOrWhiteSpace(strLineIn) Then
                            Continue Do
                        End If

                        If Not blnHeaderParsed Then
                            blnSuccess = ParseMSGFDBSynFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True
                            Continue Do
                        End If

                        blnValidSearchResult = ParseMSGFDBSynFileEntry(strLineIn, objSearchResult, strErrorLog,
                          intResultsProcessed, intColumnMapping,
                          strCurrentPeptideWithMods)

                        If blnValidSearchResult Then
                            strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge

                            If objSearchResult.SpecProb = strPreviousSpecProb Then
                                ' New result has the same SpecProb as the previous result
                                ' See if htPeptidesFoundForSpecProbLevel contains the peptide, scan and charge

                                If htPeptidesFoundForSpecProbLevel.ContainsKey(strKey) Then
                                    blnFirstMatchForGroup = False
                                Else
                                    htPeptidesFoundForSpecProbLevel.Add(strKey, 1)
                                    blnFirstMatchForGroup = True
                                End If

                            Else
                                ' New SpecProb
                                ' Reset htPeptidesFoundForSpecProbLevel
                                htPeptidesFoundForSpecProbLevel.Clear()

                                ' Update strPreviousSpecProb
                                strPreviousSpecProb = objSearchResult.SpecProb

                                ' Append a new entry to htPeptidesFoundForSpecProbLevel
                                htPeptidesFoundForSpecProbLevel.Add(strKey, 1)
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

                    Loop

                End Using

                If mCreateModificationSummaryFile Then
                    ' Create the modification summary file
                    Dim fiInputFile = New FileInfo(strInputFilePath)
                    strModificationSummaryFilePath = Path.GetFileName(MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_MOD_SUMMARY))
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

    Private Function ParseMSGFDBResultsFileEntry(
      strLineIn As String,
      blnMSGFPlus As Boolean,
      lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      lstSearchResultsCurrentScan As List(Of udtMSGFDBSearchResultType),
      ByRef strErrorLog As String,
      intColumnMapping() As Integer,
      ByRef intNextScanGroupID As Integer,
      lstScanGroupDetails As List(Of udtScanGroupInfoType),
      htScanGroupCombo As Dictionary(Of String, Boolean),
      lstSpecIdToIndex As Dictionary(Of String, Integer)) As Boolean

        ' Parses an entry from the MSGF-DB results file

        Dim udtSearchResult = New udtMSGFDBSearchResultType
        Dim strSplitLine() As String = Nothing

        Dim intScanCount As Integer
        Dim strSplitResult() As String = Nothing

        Dim udtMergedScanInfo() As udtMSGFDBSearchResultType = Nothing

        Dim lstProteinInfo As Dictionary(Of String, udtTerminusCharsType)

        Dim blnValidSearchResult As Boolean
        Dim intSlashIndex As Integer
        Dim blnTargetDecoyFDRValid As Boolean

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            lstProteinInfo = New Dictionary(Of String, udtTerminusCharsType)

            ' Reset lstSearchResults
            lstSearchResultsCurrentScan.Clear()

            udtSearchResult.Clear()
            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length >= 13 Then

                With udtSearchResult
                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpectrumFile), .SpectrumFileName) Then
                        ReportError("SpectrumFile column is missing or invalid", True)
                    End If
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpecIndex), .SpecIndex)

                    If blnMSGFPlus Then
                        Dim intSpecIndex As Integer
                        Dim blnGenerateSpecIndex = True

                        If Not Integer.TryParse(.SpecIndex, intSpecIndex) Then
                            ' MSGF+ includes text in the SpecID column, for example: "controllerType=0 controllerNumber=1 scan=6390" or "index=4323"
                            ' Need to convert these to an integer

                            If .SpecIndex.StartsWith("index=") Then
                                .SpecIndex = .SpecIndex.Substring("index=".Length)
                                If Integer.TryParse(.SpecIndex, intSpecIndex) Then
                                    blnGenerateSpecIndex = False
                                End If
                            End If

                            If blnGenerateSpecIndex Then
                                If Not lstSpecIdToIndex.TryGetValue(.SpecIndex, intSpecIndex) Then
                                    intSpecIndex = lstSpecIdToIndex.Count + 1
                                    lstSpecIdToIndex.Add(.SpecIndex, intSpecIndex)
                                End If

                                .SpecIndex = intSpecIndex.ToString()
                            End If

                        End If
                    End If

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Scan), .Scan) Then
                        ReportError("Scan column is missing or invalid", True)
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.FragMethod), .FragMethod)

                    intSlashIndex = .Scan.IndexOf("/"c)
                    If intSlashIndex > 0 Then
                        ' This is a merged spectrum and thus scan number looks like: 3010/3011/3012
                        ' Split the Scan list on the slash
                        ' Later in this function, we'll append lstSearchResults with this scan plus the other scans

                        strSplitResult = .Scan.Split("/"c)
                        intScanCount = strSplitResult.Length
                        ReDim udtMergedScanInfo(intScanCount - 1)

                        For intIndex = 0 To intScanCount - 1
                            udtMergedScanInfo(intIndex) = New udtMSGFDBSearchResultType
                            udtMergedScanInfo(intIndex).Clear()
                            udtMergedScanInfo(intIndex).Scan = strSplitResult(intIndex)
                            udtMergedScanInfo(intIndex).ScanNum = CIntSafe(strSplitResult(intIndex), 0)
                        Next

                        ' Now split SpecIndex and store in udtMergedScanInfo
                        strSplitResult = .SpecIndex.Split("/"c)

                        For intIndex = 0 To strSplitResult.Length - 1
                            If intIndex >= udtMergedScanInfo.Length Then
                                ' There are more entries for SpecIndex than there are for Scan#; this is unexpected
                                Exit For
                            End If
                            udtMergedScanInfo(intIndex).SpecIndex = strSplitResult(intIndex)
                        Next

                        ' Now split FragMethod and store in udtMergedScanInfo
                        strSplitResult = .FragMethod.Split("/"c)

                        For intIndex = 0 To strSplitResult.Length - 1
                            If intIndex >= udtMergedScanInfo.Length Then
                                ' There are more entries for FragMethod than there are for Scan#; this is unexpected
                                Exit For
                            End If
                            udtMergedScanInfo(intIndex).FragMethod = strSplitResult(intIndex)
                        Next

                    Else
                        .ScanNum = CIntSafe(.Scan, 0)
                        intScanCount = 1
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PrecursorMZ), .PrecursorMZ)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Charge), .Charge)
                    .ChargeNum = CShort(CIntSafe(.Charge, 0))

                    ' Precursor mass error could be in PPM or Da
                    '   In MSGFDB, the header line will have PMError(ppm)        or PMError(Da)
                    '   In MSGF+,  the header line will have PrecursorError(ppm) or PrecursorError(Da)
                    Dim dblPrecursorErrorDa As Double

                    If intColumnMapping(eMSGFDBResultsFileColumns.PMErrorPPM) >= 0 Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PMErrorPPM), .PMErrorPPM)
                    Else
                        GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PMErrorDa), .PMErrorDa)
                        dblPrecursorErrorDa = CDblSafe(.PMErrorDa, 0)
                        .PMErrorPPM = String.Empty              ' We'll populate this column later in this function
                    End If


                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Peptide), .Peptide) Then
                        ReportError("Peptide column is missing or invalid", True)
                    End If


                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Protein), .Protein)

                    ' MSGF+ .tsv files may have a semicolon separated list of protein names; check for this
                    .Protein = SplitProteinList(.Protein, lstProteinInfo)

                    If lstProteinInfo.Count > 0 Then
                        ' Need to add the prefix and suffix residues
                        .Peptide = AddUpdatePrefixAndSuffixResidues(.Peptide, lstProteinInfo.First)
                    End If

                    ' Replace any mod text values in the peptide sequence with the appropriate mod symbols
                    ' In addition, replace the terminus symbols with dashes
                    Dim dblTotalModMass As Double
                    .Peptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(.Peptide), lstMSGFDBModInfo, blnMSGFPlus, dblTotalModMass)

                    ' Compute monoisotopic mass of the peptide
                    Dim dblPeptideMonoisotopicMass = ComputePeptideMass(.Peptide, dblTotalModMass)

                    ' Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    .MH = NumToString(clsPeptideMassCalculator.ConvoluteMass(dblPeptideMonoisotopicMass, 0, 1), 6, True)

                    If Not String.IsNullOrEmpty(.PMErrorPPM) Then

                        ' Convert the ppm-based PM Error to Da-based

                        Dim dblPMErrorPPM As Double
                        Dim dblPrecursorMZ As Double

                        If Double.TryParse(.PrecursorMZ, dblPrecursorMZ) Then
                            ' Note that since .PMErrorPPM is present, the Precursor m/z is a C13-corrected m/z value
                            ' In other words, it may not be the actual m/z selected for fragmentation.

                            If Double.TryParse(.PMErrorPPM, dblPMErrorPPM) Then

                                If mParentMassToleranceInfo.IsPPM AndAlso
                                  (dblPMErrorPPM < -mParentMassToleranceInfo.ToleranceLeft * 1.5 OrElse
                                   dblPMErrorPPM > mParentMassToleranceInfo.ToleranceRight * 1.5) Then

                                    ' PPM error computed by MSGF+ is more than 1.5-fold larger than the ppm-based parent ion tolerance; don't trust the value computed by MSGF+

                                    mPrecursorMassErrorWarningCount += 1
                                    If mPrecursorMassErrorWarningCount <= 10 Then
                                        ReportWarning("Precursor mass error computed by MSGF+ is 1.5-fold larger than search tolerance: " & .PMErrorPPM & " vs. " & mParentMassToleranceInfo.ToleranceLeft.ToString("0") & "ppm," & mParentMassToleranceInfo.ToleranceRight.ToString("0") & "ppm")
                                        If mPrecursorMassErrorWarningCount = 10 Then
                                            ReportWarning("Additional mass errors will not be reported")
                                        End If
                                    End If


                                    Dim dblPrecursorMonoMass As Double
                                    dblPrecursorMonoMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .ChargeNum, 0)

                                    dblPrecursorErrorDa = dblPrecursorMonoMass - dblPeptideMonoisotopicMass

                                    .PMErrorPPM = String.Empty

                                Else

                                    dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPMErrorPPM, dblPeptideMonoisotopicMass)

                                    ' Note that this will be a C13-corrected precursor error; not the true precursor error
                                    .PMErrorDa = NumToString(dblPrecursorErrorDa, 6, True)
                                End If

                            End If
                        End If
                    End If

                    If String.IsNullOrEmpty(.PMErrorPPM) Then

                        Dim dblPrecursorMZ As Double
                        If Double.TryParse(.PrecursorMZ, dblPrecursorMZ) Then
                            Dim dblPeptideDeltaMassCorrectedPpm As Double

                            dblPeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMZ,
                             .ChargeNum, dblPeptideMonoisotopicMass, True)

                            .PMErrorPPM = NumToString(dblPeptideDeltaMassCorrectedPpm, 5, True)

                            If String.IsNullOrEmpty(.PMErrorDa) Then
                                dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPeptideDeltaMassCorrectedPpm, dblPeptideMonoisotopicMass)

                                ' Note that this will be a C13-corrected precursor error; not the true precursor error
                                .PMErrorDa = NumToString(dblPrecursorErrorDa, 6, True)
                            End If

                        End If

                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.DeNovoScore), .DeNovoScore)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.MSGFScore), .MSGFScore)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpecProb_EValue), .SpecProb)
                    If Not Double.TryParse(.SpecProb, .SpecProbNum) Then .SpecProbNum = 0

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PValue_EValue), .PValue)
                    If Not Double.TryParse(.PValue, .PValueNum) Then .PValueNum = 0

                    blnTargetDecoyFDRValid = GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.FDR_QValue), .FDR)
                    If blnTargetDecoyFDRValid Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PepFDR_PepQValue), .PepFDR)
                    Else
                        GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.EFDR), .FDR)
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.IsotopeError), .IsotopeError)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.IMSScan), .IMSScan)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.IMSDriftTime), .IMSDriftTime)

                    .NTT = ComputeCleaveageState(.Peptide).ToString()
                End With

                Dim udtScanGroupInfo As udtScanGroupInfoType
                Dim intCurrentScanGroupID As Integer = -1

                udtScanGroupInfo.Charge = udtSearchResult.ChargeNum

                If intScanCount > 1 Then
                    ' Append one entry to lstSearchResults for each item in udtMergedScanInfo()

                    For intIndex = 0 To intScanCount - 1
                        udtSearchResult.Scan = udtMergedScanInfo(intIndex).Scan
                        udtSearchResult.ScanNum = udtMergedScanInfo(intIndex).ScanNum

                        udtSearchResult.SpecIndex = udtMergedScanInfo(intIndex).SpecIndex
                        udtSearchResult.FragMethod = udtMergedScanInfo(intIndex).FragMethod

                        AppendToSearchResults(lstSearchResultsCurrentScan, udtSearchResult, lstProteinInfo)

                        ' Append an entry to lstScanGroupDetails
                        udtScanGroupInfo.Scan = udtSearchResult.ScanNum
                        AppendToScanGroupDetails(lstScanGroupDetails, htScanGroupCombo, udtScanGroupInfo, intCurrentScanGroupID, intNextScanGroupID)
                    Next
                Else
                    ' This is not a merged result; simply append udtSearchResult to lstSearchResults
                    AppendToSearchResults(lstSearchResultsCurrentScan, udtSearchResult, lstProteinInfo)

                    ' Also append an entry to lstScanGroupDetails
                    udtScanGroupInfo.Scan = udtSearchResult.ScanNum
                    AppendToScanGroupDetails(lstScanGroupDetails, htScanGroupCombo, udtScanGroupInfo, intCurrentScanGroupID, intNextScanGroupID)

                End If

                blnValidSearchResult = True
            End If

        Catch ex As Exception
            ' Error parsing this row from the MSGF-DB results file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MSGF-DB Results in ParseMSGFDBResultsFileEntry for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MSGF-DB Results in ParseMSGFDBResultsFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    Private Function ParseMSGFDBResultsFileHeaderLine(strLineIn As String, ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        ' The expected header from MSGFDB is:
        ' #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(Da)           Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR
        ' or                                                                                                  
        ' #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(ppm)          Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR

        ' The expected header from MSGF+ is:
        ' #SpecFile    SpecID       ScanNum   FragMethod    Precursor    IsotopeError    PrecursorError(Da)    Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue
        ' or
        ' #SpecFile    SpecID       ScanNum   FragMethod    Precursor    IsotopeError    PrecursorError(ppm)   Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMSGFDBResultsFileColumns
        Dim lstColumnNames As SortedDictionary(Of String, eMSGFDBResultsFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMSGFDBResultsFileColumns)(StringComparer.InvariantCultureIgnoreCase)

        ReDim intColumnMapping(MSGFDBResultsFileColCount - 1)

        lstColumnNames.Add("#SpecFile", eMSGFDBResultsFileColumns.SpectrumFile)

        lstColumnNames.Add("SpecIndex", eMSGFDBResultsFileColumns.SpecIndex)
        lstColumnNames.Add("SpecID", eMSGFDBResultsFileColumns.SpecIndex)

        lstColumnNames.Add("Scan#", eMSGFDBResultsFileColumns.Scan)
        lstColumnNames.Add("ScanNum", eMSGFDBResultsFileColumns.Scan)

        lstColumnNames.Add("FragMethod", eMSGFDBResultsFileColumns.FragMethod)
        lstColumnNames.Add("Precursor", eMSGFDBResultsFileColumns.PrecursorMZ)
        lstColumnNames.Add("IsotopeError", eMSGFDBResultsFileColumns.IsotopeError)

        lstColumnNames.Add("PMError(Da)", eMSGFDBResultsFileColumns.PMErrorDa)
        lstColumnNames.Add("PrecursorError(Da)", eMSGFDBResultsFileColumns.PMErrorDa)

        lstColumnNames.Add("PMError(ppm)", eMSGFDBResultsFileColumns.PMErrorPPM)
        lstColumnNames.Add("PrecursorError(ppm)", eMSGFDBResultsFileColumns.PMErrorPPM)

        lstColumnNames.Add("Charge", eMSGFDBResultsFileColumns.Charge)
        lstColumnNames.Add("Peptide", eMSGFDBResultsFileColumns.Peptide)
        lstColumnNames.Add("Protein", eMSGFDBResultsFileColumns.Protein)
        lstColumnNames.Add("DeNovoScore", eMSGFDBResultsFileColumns.DeNovoScore)
        lstColumnNames.Add("MSGFScore", eMSGFDBResultsFileColumns.MSGFScore)

        lstColumnNames.Add("SpecProb", eMSGFDBResultsFileColumns.SpecProb_EValue)
        lstColumnNames.Add("SpecEValue", eMSGFDBResultsFileColumns.SpecProb_EValue)

        lstColumnNames.Add("P-value", eMSGFDBResultsFileColumns.PValue_EValue)
        lstColumnNames.Add("EValue", eMSGFDBResultsFileColumns.PValue_EValue)

        lstColumnNames.Add("FDR", eMSGFDBResultsFileColumns.FDR_QValue)
        lstColumnNames.Add("QValue", eMSGFDBResultsFileColumns.FDR_QValue)

        lstColumnNames.Add("PepFDR", eMSGFDBResultsFileColumns.PepFDR_PepQValue)
        lstColumnNames.Add("PepQValue", eMSGFDBResultsFileColumns.PepFDR_PepQValue)

        lstColumnNames.Add("EFDR", eMSGFDBResultsFileColumns.EFDR)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Scan, eMSGFDBResultsFileColumns.IMSScan)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Drift_Time, eMSGFDBResultsFileColumns.IMSDriftTime)

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
                Else
                    ' Unrecognized column name
                    Console.WriteLine("Warning: Unrecognized column header name '" & strSplitLine(intIndex) & "' in ParseMSGFDBResultsFileHeaderLine")
                End If
            Next

        Catch ex As Exception
            SetErrorMessage("Error parsing header in MSGFDB results file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMSGFDBSynFileHeaderLine(strLineIn As String, ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMSFDBSynFileColumns
        Dim lstColumnNames As SortedDictionary(Of String, eMSFDBSynFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMSFDBSynFileColumns)(StringComparer.InvariantCultureIgnoreCase)

        ReDim intColumnMapping(MSGFDBSynFileColCount - 1)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_ResultID, eMSFDBSynFileColumns.ResultID)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Scan, eMSFDBSynFileColumns.Scan)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_FragMethod, eMSFDBSynFileColumns.FragMethod)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_SpecIndex, eMSFDBSynFileColumns.SpecIndex)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Charge, eMSFDBSynFileColumns.Charge)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PrecursorMZ, eMSFDBSynFileColumns.PrecursorMZ)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_DelM, eMSFDBSynFileColumns.DelM)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_DelM_PPM, eMSFDBSynFileColumns.DelMPPM)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MH, eMSFDBSynFileColumns.MH)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Peptide, eMSFDBSynFileColumns.Peptide)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Protein, eMSFDBSynFileColumns.Protein)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_NTT, eMSFDBSynFileColumns.NTT)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore, eMSFDBSynFileColumns.DeNovoScore)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, eMSFDBSynFileColumns.MSGFScore)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb, eMSFDBSynFileColumns.SpecProb_EValue)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecEValue, eMSFDBSynFileColumns.SpecProb_EValue)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecProb, eMSFDBSynFileColumns.RankSpecProb)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecEValue, eMSFDBSynFileColumns.RankSpecProb)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PValue, eMSFDBSynFileColumns.PValue_EValue)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EValue, eMSFDBSynFileColumns.PValue_EValue)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_FDR, eMSFDBSynFileColumns.FDR_QValue)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_QValue, eMSFDBSynFileColumns.FDR_QValue)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR, eMSFDBSynFileColumns.PepFDR_PepQValue)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepQValue, eMSFDBSynFileColumns.PepFDR_PepQValue)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EFDR, eMSFDBSynFileColumns.EFDR)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Isotope_Error, eMSFDBSynFileColumns.IsotopeError)

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
            SetErrorMessage("Error parsing header in MSGFDB synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMSGFDBSynFileEntry(
      strLineIn As String,
      objSearchResult As clsSearchResultsMSGFDB,
      ByRef strErrorLog As String,
      intResultsProcessed As Integer,
      ByRef intColumnMapping() As Integer,
      <Out()> ByRef strPeptideSequenceWithMods As String) As Boolean

        ' Parses an entry from the MSGFDB Synopsis file

        Dim strSplitLine() As String = Nothing
        Dim strValue As String = String.Empty

        Dim blnValidSearchResult As Boolean
        Dim blnTargetDecoyFDRValid As Boolean

        strPeptideSequenceWithMods = String.Empty

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset objSearchResult
            objSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length >= 15 Then

                With objSearchResult
                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.ResultID), strValue) Then
                        If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                            strErrorLog &= "Error reading ResultID value from MSGFDB Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                        End If
                        Exit Try
                    End If

                    .ResultID = Integer.Parse(strValue)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.Scan), .Scan)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.Charge), .Charge)

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.Peptide), strPeptideSequenceWithMods) Then
                        If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                            strErrorLog &= "Error reading Peptide sequence value from MSGFDB Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                        End If
                        Exit Try
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.Protein), .ProteinName)
                    .MultipleProteinCount = "0"

                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.DelM), .MSGFDbComputedDelM)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.DelMPPM), .MSGFDbComputedDelMPPM)

                    .PeptideDeltaMass = .MSGFDbComputedDelM

                    ' Note: .PeptideDeltaMass is stored in the MSGF-DB results file as "Observed_Mass - Theoretical_Mass"
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
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.FragMethod), .FragMethod)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PrecursorMZ), .PrecursorMZ)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.MH), .PeptideMH)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.NTT), .NTT)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.DeNovoScore), .DeNovoScore)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.MSGFScore), .MSGFScore)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.SpecProb_EValue), .SpecProb)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.RankSpecProb), .RankSpecProb)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PValue_EValue), .PValue)

                    blnTargetDecoyFDRValid = GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.FDR_QValue), .FDR)
                    If blnTargetDecoyFDRValid Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PepFDR_PepQValue), .PepFDR)
                    Else
                        GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.EFDR), .FDR)
                    End If

                    If intColumnMapping(eMSFDBSynFileColumns.IsotopeError) >= 0 Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.IsotopeError), .IsotopeError)
                        .MSGFPlusResults = True
                    Else
                        .MSGFPlusResults = False
                    End If

                    ' Compute PrecursorMH using PrecursorMZ and charge
                    Dim dblPrecursorMZ As Double
                    Dim intCharge As Integer
                    If Double.TryParse(.PrecursorMZ, dblPrecursorMZ) Then
                        If Integer.TryParse(.Charge, intCharge) Then
                            .ParentIonMH = NumToString(clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, intCharge, 1), 6, True)
                        End If
                    End If

                End With

                blnValidSearchResult = True
            End If

        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MSGFDB Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MSGFDB Results in ParseMSGFDBSynFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    Protected Function ParseParentMassTolerance(strToleranceText As String, <Out()> ByRef dblTolerance As Double, <Out()> ByRef blnIsPPM As Boolean) As Boolean
        dblTolerance = 0
        blnIsPPM = False

        strToleranceText = strToleranceText.ToLower().Trim()

        If strToleranceText.EndsWith("da") Then
            strToleranceText = strToleranceText.Substring(0, strToleranceText.Length - 2)
            blnIsPPM = False

        ElseIf strToleranceText.EndsWith("ppm") Then
            strToleranceText = strToleranceText.Substring(0, strToleranceText.Length - 3)
            blnIsPPM = True

        Else
            Return False
        End If


        If Double.TryParse(strToleranceText, dblTolerance) Then
            Return True
        Else
            Return False
        End If

    End Function

    ''' <summary>
    ''' Main processing function
    ''' </summary>
    ''' <param name="strInputFilePath">MSGFDB results file</param>
    ''' <param name="strOutputFolderPath">Output folder</param>
    ''' <param name="strParameterFilePath">Parameter file for data processing</param>
    ''' <returns>True if success, False if failure</returns>
    ''' <remarks>Use SearchToolParameterFilePath to define the search engine parameter file</remarks>
    Public Overloads Overrides Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean

        Dim strFhtOutputFilePath As String = String.Empty
        Dim strSynOutputFilePath As String = String.Empty
        Dim strPepToProteinMapFilePath As String
        Dim strScanGroupFilePath As String

        Dim lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)
        Dim lstPepToProteinMapping As List(Of udtPepToProteinMappingType)
        Dim strMTSPepToProteinMapFilePath As String = String.Empty

        Dim blnMSGFPlus As Boolean
        Dim lstSpecIdToIndex As Dictionary(Of String, Integer)

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

            mPeptideSeqMassCalculator.ResetAminoAcidMasses()

            lstSpecIdToIndex = New Dictionary(Of String, Integer)

            MyBase.ResetProgress("Parsing " & Path.GetFileName(strInputFilePath))

            If Not CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                Return False
            End If

            Try
                ' Obtain the full path to the input file
                Dim fiInputFile = New FileInfo(strInputFilePath)

                lstMSGFDBModInfo = New List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)
                lstPepToProteinMapping = New List(Of udtPepToProteinMappingType)

                ' Load the MSGF-DB Parameter File so that we can determine the modification names and masses
                ' If the MSGFDB_Mods.txt file was defined, then the mod symbols in that file will be used to define the mod symbols in lstMSGFDBModInfo 
                Dim success = ExtractModInfoFromParamFile(mSearchToolParameterFilePath, lstMSGFDBModInfo)
                If Not success Then
                    Return False
                End If

                mParentMassToleranceInfo = ExtractParentMassToleranceFromParamFile(mSearchToolParameterFilePath)

                Dim query = From item In lstMSGFDBModInfo Where item.ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA
                If query.Any() Then
                    ' Custom amino acids are defined; read their values and update the mass calculator

                    Dim localErrorMsg As String = String.Empty
                    Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor("MSGF+")

                    AddHandler modFileProcessor.ErrorOccurred, AddressOf ModExtractorErrorHandler
                    AddHandler modFileProcessor.WarningMessageEvent, AddressOf ModExtractorWarningHandler

                    clsPHRPParserMSGFDB.UpdateMassCalculatorMasses(mSearchToolParameterFilePath, modFileProcessor, mPeptideSeqMassCalculator, localErrorMsg)

                    If Not String.IsNullOrWhiteSpace(localErrorMsg) AndAlso String.IsNullOrWhiteSpace(mErrorMessage) Then
                        ReportError(localErrorMsg)
                    End If

                End If

                ' Define the base output filename using strInputFilePath
                Dim strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath)

                ' Auto-replace "_msgfplus" with "_msgfdb"
                If strBaseName.ToLower().EndsWith("_msgfplus") Then
                    strBaseName = strBaseName.Substring(0, strBaseName.Length - "_msgfplus".Length) & "_msgfdb"
                End If

                If MyBase.mCreateInspectOrMSGFDbFirstHitsFile Then

                    ' Read the FASTA file to cache the protein names in memory
                    ' These will be used when creating the first hits file
                    If Not CacheProteinNamesFromFasta() Then
                        Return False
                    End If

                    ' Create the first hits output file
                    MyBase.ResetProgress("Creating the FHT file", True)

                    strFhtOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_FIRST_HITS_FILE_SUFFIX)

                    strScanGroupFilePath = String.Empty

                    blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strFhtOutputFilePath, strScanGroupFilePath, lstMSGFDBModInfo, blnMSGFPlus, lstSpecIdToIndex, eFilteredOutputFileTypeConstants.FHTFile)

                End If

                If MyBase.mCreateInspectOrMSGFDbSynopsisFile Then

                    ' Create the synopsis output file
                    MyBase.ResetProgress("Creating the SYN file", True)

                    ' The synopsis file name will be of the form BasePath_msgfdb_syn.txt
                    strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                    strScanGroupFilePath = Path.Combine(strOutputFolderPath, strBaseName & "_ScanGroupInfo.txt")

                    blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, strScanGroupFilePath, lstMSGFDBModInfo, blnMSGFPlus, lstSpecIdToIndex, eFilteredOutputFileTypeConstants.SynFile)

                    ' Load the PeptideToProteinMap information; if the file doesn't exist, then a warning will be displayed, but processing will continue
                    ' LoadPeptideToProteinMapInfoMSGFDB also creates _msgfdb_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols							
                    strPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(Path.Combine(strOutputFolderPath, strBaseName) & ".txt", strOutputFolderPath, MTS:=False)

                    MyBase.ResetProgress("Loading the PepToProtein map file: " & Path.GetFileName(strPepToProteinMapFilePath), True)

                    LoadPeptideToProteinMapInfoMSGFDB(strPepToProteinMapFilePath, strOutputFolderPath, lstMSGFDBModInfo, blnMSGFPlus, lstPepToProteinMapping, strMTSPepToProteinMapFilePath)

                    ' Create the other PHRP-specific files
                    MyBase.ResetProgress("Creating the PHRP files for " & Path.GetFileName(strSynOutputFilePath), True)

                    ' Now parse the _syn.txt file that we just created to next create the other PHRP files
                    blnSuccess = ParseMSGFDBSynopsisFile(strSynOutputFilePath, strOutputFolderPath, lstPepToProteinMapping, False)

                    ' Remove all items from lstPepToProteinMapping to reduce memory overhead
                    lstPepToProteinMapping.Clear()
                    lstPepToProteinMapping.TrimExcess()

                    If blnSuccess AndAlso mCreateProteinModsFile Then
                        blnSuccess = CreateProteinModsFileWork(strBaseName, fiInputFile, strFhtOutputFilePath, strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath)
                    End If

                End If

                If blnSuccess Then
                    MyBase.OperationComplete()
                End If

            Catch ex As Exception
                SetErrorMessage("Error in clsMSGFDBResultsProcessor.ProcessFile (2):  " & ex.Message)
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
            End Try

        Catch ex As Exception
            SetErrorMessage("Error in clsMSGFDBResultsProcessor.ProcessFile (1):" & ex.Message)
            SetErrorCode(ePHRPErrorCodes.UnspecifiedError)
        End Try

        Return blnSuccess

    End Function

    Private Function CreateProteinModsFileWork(
       strBaseName As String,
       fiInputFile As FileInfo,
       strFhtOutputFilePath As String,
       strSynOutputFilePath As String,
       strOutputFolderPath As String,
       strMTSPepToProteinMapFilePath As String
       ) As Boolean

        Dim blnSuccess = True

        If String.IsNullOrEmpty(strMTSPepToProteinMapFilePath) OrElse Not File.Exists(strMTSPepToProteinMapFilePath) Then
            ' MTSPepToProteinMap file not found; auto-create it

            If String.IsNullOrEmpty(strMTSPepToProteinMapFilePath) Then
                strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS:=True)
            End If

            Dim lstSourcePHRPDataFiles = New List(Of String)

            If Not String.IsNullOrEmpty(strFhtOutputFilePath) Then
                lstSourcePHRPDataFiles.Add(strFhtOutputFilePath)
            End If

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
                    ' We only do this for MSGFDB since it often includes reverse protein peptides in the results even though the FASTA file often does not have reverse proteins
                    mIgnorePeptideToProteinMapperErrors = True
                    blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
                    If Not blnSuccess Then
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False")
                    End If
                End If
            End If

        End If

        If blnSuccess Then
            ' If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
            MyBase.ValidatePHRPReaderSupportFiles(Path.Combine(fiInputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath)

            ' Create the Protein Mods file
            blnSuccess = MyBase.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB)
        End If

        If Not blnSuccess Then
            ' Do not treat this as a fatal error
            blnSuccess = True
        End If

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Replaces modification masses in peptide sequences with modification symbols (uses case-sensitive comparisons)
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="lstMSGFDBModInfo">This function assumes that each entry in lstMSGFDBModInfo has both .ModName and .ModSymbol defined</param>
    ''' <param name="blnMSGFPlus">Should be set to True if processing MSGF+ results</param>
    ''' <param name="dblTotalModMass">Output parameter: total mass of all modifications</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ReplaceMSGFModTextWithSymbol(
      strPeptide As String,
      lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      blnMSGFPlus As Boolean,
      <Out()> ByRef dblTotalModMass As Double) As String

        Dim intIndex As Integer
        Dim intIndexFirstResidue As Integer
        Dim intIndexLastResidue As Integer

        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Dim strModSymbols As String = String.Empty

        Dim blnNterminalMod As Boolean
        Dim blnPossibleCTerminalMod As Boolean
        Dim blnIsStaticMod As Boolean

        Static reNTerminalModMassRegEx As New Regex(MSGFDB_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS)
        Static reModMassRegEx As New Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS)

        Dim reMatch As Match

        Dim dblModMassFound As Double

        ' Reset the total mod mass
        dblTotalModMass = 0

        ' Remove the prefix and suffix residues
        If strPeptide.Length >= 4 Then
            If strPeptide.Chars(1) = "."c AndAlso
               strPeptide.Chars(strPeptide.Length - 2) = "."c Then

                strPrefix = strPeptide.Substring(0, 2)
                strSuffix = strPeptide.Substring(strPeptide.Length - 2, 2)

                strPeptide = strPeptide.Substring(2, strPeptide.Length - 4)
            End If
        End If

        ' strPeptide should now be the primary peptide sequence, without the prefix or suffix residues

        ' First look for dynamic N-terminal mods (NTermPeptide or NTermProtein)
        ' This RegEx will match one or more mods, all at the N-terminus
        reMatch = reNTerminalModMassRegEx.Match(strPeptide)
        If reMatch IsNot Nothing AndAlso reMatch.Success Then

            blnNterminalMod = True
            blnPossibleCTerminalMod = False
            blnIsStaticMod = False

            ' Convert the mod mass (or masses) to one or more mod symbols
            If ConvertMGSFModMassesToSymbols("-", reMatch.Groups(1).Value, strModSymbols, lstMSGFDBModInfo, blnNterminalMod, blnPossibleCTerminalMod, dblModMassFound, blnIsStaticMod) Then

                ' Replace the mod digits with the mod symbols

                strPeptide = ReplaceMSGFModTextWithMatchedSymbol(strPeptide, reMatch.Groups(1), strModSymbols, blnMSGFPlus, blnIsStaticMod)
                dblTotalModMass += dblModMassFound

            End If
        End If

        ' Next, step through the peptide and parse each mod mass that follows a residue
        ' Any mod mass at the end must be considered a C-terminal mod 

        ' Need to start at the first letter
        ' If we had N-terminal mods, they're currently notated like this: _.+42.011MDHTPQSQLK.L or _.+42.011+57.021MNDR.Q
        ' We want things to look like this: -.#MDHTPQSQLK.L or -.#*MNDRQLNHR.S

        ' In MSGFDB, static mods do not have a mod mass listed
        ' In MSGF+,  static mods do have a mod mass listed
        ' Regardless, we do not add mod symbols for static mods, but we do increment dblTotalModMass

        ' Find the index of the last residue
        intIndex = strPeptide.Length - 1
        Do While intIndex > 0 AndAlso Not IsLetterAtoZ(strPeptide.Chars(intIndex))
            intIndex -= 1
        Loop
        intIndexLastResidue = intIndex

        ' Find the index of the first residue
        intIndex = 0
        Do While intIndex < strPeptide.Length AndAlso Not IsLetterAtoZ(strPeptide.Chars(intIndex))
            intIndex += 1
        Loop
        intIndexFirstResidue = intIndex

        Dim currentResidue = "-"

        Do While intIndex < strPeptide.Length

            If IsLetterAtoZ(strPeptide.Chars(intIndex)) Then
                currentResidue = strPeptide.Chars(intIndex)

                Dim objModificationDefinition As clsModificationDefinition

                If Not blnMSGFPlus Then

                    ' Look for static mods that should be applied to this residue (only applies to MSGFDB, not MSGF+)
                    For intModIndex = 0 To mPeptideMods.ModificationCount - 1
                        Dim eModificationType As clsModificationDefinition.eModificationTypeConstants
                        eModificationType = mPeptideMods.GetModificationTypeByIndex(intModIndex)

                        If eModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                            objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)

                            If objModificationDefinition.TargetResiduesContain(strPeptide.Chars(intIndex)) Then
                                ' Match found; update dblTotalModMass but do not add a static mod symbol
                                dblTotalModMass += objModificationDefinition.ModificationMass
                            End If

                        ElseIf intIndex = intIndexFirstResidue Then
                            If eModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod AndAlso strPrefix = "_" Then
                                ' N-terminal protein static mod
                                objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)
                                dblTotalModMass += objModificationDefinition.ModificationMass
                            ElseIf eModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod Then
                                ' N-terminal peptide static mod
                                objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)
                                dblTotalModMass += objModificationDefinition.ModificationMass
                            End If
                        End If
                    Next intModIndex

                End If

                intIndex += 1

                If intIndex = intIndexLastResidue Then blnPossibleCTerminalMod = True

            Else
                ' Found a mod; find the extent of the mod digits
                reMatch = reModMassRegEx.Match(strPeptide, intIndex)

                ' Note that blnPossibleCTerminalMod will be set to True once we hit the last residue
                ' Assure blnNterminalMod is false
                blnNterminalMod = False

                ' Convert the mod mass (or masses) to one or more mod symbols
                If ConvertMGSFModMassesToSymbols(currentResidue, reMatch.Groups(1).Value, strModSymbols, lstMSGFDBModInfo,
                  blnNterminalMod, blnPossibleCTerminalMod, dblModMassFound, blnIsStaticMod) Then

                    strPeptide = ReplaceMSGFModTextWithMatchedSymbol(strPeptide, reMatch.Groups(1), strModSymbols, blnMSGFPlus, blnIsStaticMod)
                    dblTotalModMass += dblModMassFound

                    If blnMSGFPlus AndAlso blnIsStaticMod Then
                        ' MSGF+ shows mod masses for static mods
                        ' Thus, we have removed the static mod mass and did not add a mod symbol
                        ' Therefore, leave intIndex unchanged
                    Else
                        intIndex += strModSymbols.Length
                    End If

                Else
                    intIndex += reMatch.Groups(1).Value.Length
                End If

            End If

        Loop

        ' If any N-terminal mods were present, we need to move them to after the first residue
        ' in other words, switch from #MDHTPQSQLK to M#DHTPQSQLK
        '                          or #*MNDRQLNHR to M#*NDRQLNHR

        ' Update intIndexFirstResidue
        intIndexFirstResidue = 0
        Do While intIndexFirstResidue < strPeptide.Length AndAlso Not IsLetterAtoZ(strPeptide.Chars(intIndexFirstResidue))
            intIndexFirstResidue += 1
        Loop

        If intIndexFirstResidue > 0 AndAlso intIndexFirstResidue < strPeptide.Length Then
            Dim strPeptideNew As String
            strPeptideNew = strPeptide.Chars(intIndexFirstResidue) & strPeptide.Substring(0, intIndexFirstResidue)
            If intIndexFirstResidue < strPeptide.Length - 1 Then
                strPeptideNew &= strPeptide.Substring(intIndexFirstResidue + 1)
            End If
            strPeptide = String.Copy(strPeptideNew)
        End If

        Return strPrefix & strPeptide & strSuffix

    End Function
    
    Protected Function ReplaceMSGFModTextWithMatchedSymbol(
      strPeptide As String,
      reGroup As Group,
      strModSymbols As String,
      blnMSGFPlus As Boolean,
      blnIsStaticMod As Boolean) As String

        Dim strPeptideNew As String

        If reGroup.Index > 0 Then
            strPeptideNew = strPeptide.Substring(0, reGroup.Index)
        Else
            strPeptideNew = String.Empty
        End If

        If blnMSGFPlus AndAlso blnIsStaticMod Then
            ' MSGF+ shows mod masses for static mods
            ' However, for consistency with other PHRP results, we do not add a symbol to the peptide for this static mod
        Else
            strPeptideNew &= strModSymbols
        End If

        If reGroup.Index + reGroup.Length < strPeptide.Length Then
            strPeptideNew &= strPeptide.Substring(reGroup.Index + reGroup.Length)
        End If

        Return strPeptideNew

    End Function

    Protected Function ReplaceTerminus(strPeptide As String) As String

        If strPeptide.StartsWith(N_TERMINUS_SYMBOL_MSGFDB) Then
            strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST & "." & strPeptide.Substring(N_TERMINUS_SYMBOL_MSGFDB.Length)
        End If

        If strPeptide.EndsWith(C_TERMINUS_SYMBOL_MSGFDB) Then
            strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_MSGFDB.Length) & "." & clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
        End If

        Return strPeptide

    End Function

    ''' <summary>
    ''' Examines strProteinList to look for a semi-colon separated list of proteins and terminus symbols, for example
    ''' AT1G26570.1(pre=K,post=N);AT3G29360.1(pre=K,post=N);AT3G29360.2(pre=K,post=N)
    ''' </summary>
    ''' <param name="strProteinList">Protein list to examine</param>
    ''' <param name="lstProteinInfo">Protein information, if it is of the form ProteinName(pre=X,post=Y)</param>
    ''' <returns>The name of the first protein</returns>
    ''' <remarks></remarks>
    Protected Function SplitProteinList(strProteinList As String, lstProteinInfo As Dictionary(Of String, udtTerminusCharsType)) As String

        Static reProteinInfo As New Regex(PROTEIN_AND_TERM_SYMBOLS_REGEX, REGEX_OPTIONS)

        Dim reMatches As MatchCollection

        lstProteinInfo.Clear()

        reMatches = reProteinInfo.Matches(strProteinList)

        If reMatches.Count = 0 Then
            ' No match; likely just one protein
            Return TruncateProteinName(strProteinList)
        Else
            For Each reMatch As Match In reMatches
                Dim strProteinName As String
                Dim udtTerminusChars As udtTerminusCharsType

                strProteinName = TruncateProteinName(reMatch.Groups(1).Value)
                udtTerminusChars.NTerm = reMatch.Groups(2).Value.Chars(0)
                udtTerminusChars.CTerm = reMatch.Groups(3).Value.Chars(0)

                If lstProteinInfo.ContainsKey(strProteinName) Then
                    ' Skip this protein since it's already present
                Else
                    lstProteinInfo.Add(strProteinName, udtTerminusChars)
                End If

            Next

            Return lstProteinInfo.First.Key
        End If


    End Function

    Private Sub SortAndWriteFilteredSearchResults(
      swResultFile As StreamWriter,
      lstFilteredSearchResults As List(Of udtMSGFDBSearchResultType),
      ByRef strErrorLog As String,
      blnIncludeFDRandPepFDR As Boolean,
      blnIncludeEFDR As Boolean,
      blnIncludeIMSFields As Boolean,
      blnMSGFPlus As Boolean)

        ' Sort udtFilteredSearchResults by ascending SpecProb, ascending scan, ascending charge, ascending peptide, and ascending protein
        lstFilteredSearchResults.Sort(New MSGFDBSearchResultsComparerSpecProbScanChargePeptide)

        For intIndex = 0 To lstFilteredSearchResults.Count - 1
            WriteSearchResultToFile(intIndex + 1, swResultFile, lstFilteredSearchResults(intIndex), strErrorLog, blnIncludeFDRandPepFDR, blnIncludeEFDR, blnIncludeIMSFields, blnMSGFPlus)
        Next intIndex

    End Sub

    Private Sub StoreScanGroupInfo(strScanGroupFilePath As String, lstScanGroupDetails As List(Of udtScanGroupInfoType))

        Dim intScanGroupIDPrevious As Integer
        Dim blnCreateFile As Boolean

        Try

            ' Only create the ScanGroup file if one or more scan groups exist
            ' Step through lstScanGroupDetails to check for this
            intScanGroupIDPrevious = -1
            blnCreateFile = False
            For Each udtScanGroupInfo As udtScanGroupInfoType In lstScanGroupDetails
                If udtScanGroupInfo.ScanGroupID = intScanGroupIDPrevious Then
                    blnCreateFile = True
                    Exit For
                End If
                intScanGroupIDPrevious = udtScanGroupInfo.ScanGroupID
            Next

            If blnCreateFile Then
                Using swScanGroupFile = New StreamWriter(New FileStream(strScanGroupFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                    swScanGroupFile.WriteLine("Scan_Group_ID" & ControlChars.Tab & "Charge" & ControlChars.Tab & "Scan")

                    For Each udtScanGroupInfo As udtScanGroupInfoType In lstScanGroupDetails
                        With udtScanGroupInfo
                            swScanGroupFile.WriteLine(.ScanGroupID & ControlChars.Tab & .Charge & ControlChars.Tab & .Scan)
                        End With
                    Next

                End Using
            End If

        Catch ex As Exception
            SetErrorMessage("Error creating ScanGroupInfo file: " & ex.Message)
        End Try

    End Sub

    ''' <summary>
    ''' Stores the first hits file matches for a single scan
    ''' </summary>
    ''' <param name="lstSearchResults">Search results</param>
    ''' <param name="intStartIndex">Start index for data in this scan</param>
    ''' <param name="intEndIndex">End index for data in this scan</param>
    ''' <param name="lstFilteredSearchResults">Output parmaeter: the actual filtered search results</param>
    ''' <remarks></remarks>
    Private Sub StoreTopFHTMatch(
      lstSearchResults As IList(Of udtMSGFDBSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer,
      lstFilteredSearchResults As List(Of udtMSGFDBSearchResultType))

        AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure should have already sorted by scan, charge, and SpecProb; no need to re-sort

        ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1)

        ' Now store the first match for each charge for this scan
        ' When storing, we use the protein name that occurred first in the FASTA file

        Dim udtCurrentResult = lstSearchResults(intStartIndex)
        Dim intCurrentCharge As Short = udtCurrentResult.ChargeNum
        Dim currentProteinNumber = Int32.MaxValue

        For intIndex = intStartIndex To intEndIndex
            If intCurrentCharge <> lstSearchResults(intIndex).ChargeNum Then
                ' New charge state
                ' Store udtCurrentResult (from the previous charge state)
                lstFilteredSearchResults.Add(udtCurrentResult)

                udtCurrentResult = lstSearchResults(intIndex)
                intCurrentCharge = udtCurrentResult.ChargeNum
                currentProteinNumber = Int32.MaxValue
            End If

            ' Lookup the protein number (to make sure we use the protein name that occurs first in the FASTA file)
            Dim proteinNumber As Integer
            Dim candidateName = lstSearchResults(intIndex).Protein

            If mProteinNameOrder.TryGetValue(candidateName, proteinNumber) Then
                If proteinNumber < currentProteinNumber Then
                    currentProteinNumber = proteinNumber
                    udtCurrentResult.Protein = candidateName
                End If
            Else
                ' Protein not found in mProteinNameOrder
                ' It's likely a reverse-hit protein
            End If

        Next intIndex

        ' Store udtCurrentResult (from the previous charge state)
        lstFilteredSearchResults.Add(udtCurrentResult)

        ' Alternative method using clsEngineResultsMSGFDB
        ' Dim lstFilteredSearchResultsNew = clsEngineResultUtilities.StoreTopFHTMatch(mProteinNameOrder, lstSearchResults, intStartIndex, intEndIndex)

        ' For Each newResult In lstFilteredSearchResultsNew
        '     lstFilteredSearchResults.Add(CType(lstFilteredSearchResultsNew, clsEngineResultsMSGFDB))
        ' Next

    End Sub


    ''' <summary>
    ''' Stores the synopsis file matches for a single scan
    ''' </summary>
    ''' <param name="lstSearchResults">Search results</param>
    ''' <param name="intStartIndex">Start index for data in this scan</param>
    ''' <param name="intEndIndex">End index for data in this scan</param>
    ''' <param name="lstFilteredSearchResults">Output parmaeter: the actual filtered search results</param>
    ''' <remarks></remarks>
    Private Sub StoreSynMatches(
      lstSearchResults As List(Of udtMSGFDBSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer,
      lstFilteredSearchResults As List(Of udtMSGFDBSearchResultType))

        Dim intIndex As Integer

        AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure already sorted by scan, charge, and SpecProb; no need to re-sort

        ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1)

        ' Now store or write out the matches that pass the filters
        ' By default, filter passing peptides have MSGFDB_SpecEValue <= 0.0001 Or EValue <= DEFAULT_SYN_FILE_PVALUE_THRESHOLD
        For intIndex = intStartIndex To intEndIndex
            If lstSearchResults(intIndex).PValueNum <= mMSGFDBSynopsisFilePValueThreshold OrElse
               lstSearchResults(intIndex).SpecProbNum <= mMSGFDBSynopsisFileSpecProbThreshold Then
                lstFilteredSearchResults.Add(lstSearchResults(intIndex))
            End If
        Next intIndex

    End Sub

    Private Sub WriteSynFHTFileHeader(
      swResultFile As StreamWriter,
      ByRef strErrorLog As String,
      blnIncludeFDRandPepFDR As Boolean,
      blnIncludeEFDR As Boolean,
      blnIncludeIMSFields As Boolean,
      blnMSGFPlus As Boolean)

        ' Write out the header line for synopsis / first hits files
        Try
            Dim lstData As New List(Of String)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_ResultID)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Scan)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_FragMethod)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_SpecIndex)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Charge)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PrecursorMZ)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_DelM)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_DelM_PPM)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MH)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Peptide)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Protein)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_NTT)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore)
            lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore)

            If blnMSGFPlus Then
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecEValue)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecEValue)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EValue)
            Else
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecProb)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PValue)
            End If

            If blnIncludeFDRandPepFDR Then

                If blnMSGFPlus Then
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_QValue)
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepQValue)
                Else
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_FDR)
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR)
                End If

            ElseIf blnIncludeEFDR Then
                ' Note that we'll write out a "1" for "PepFDR" for every result
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EFDR)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR)
            End If

            If blnMSGFPlus Then
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Isotope_Error)
            End If

            If blnIncludeIMSFields Then
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Scan)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Drift_Time)
            End If

            swResultFile.WriteLine(CollapseList(lstData))

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits header" & ControlChars.NewLine
            End If
        End Try

    End Sub

    ''' <summary>
    ''' Writes an entry to a synopsis or first hits file
    ''' </summary>
    ''' <param name="intResultID"></param>
    ''' <param name="swResultFile"></param>
    ''' <param name="udtSearchResult"></param>
    ''' <param name="strErrorLog"></param>
    ''' <param name="blnIncludeFDRandPepFDR"></param>
    ''' <param name="blnIncludeEFDR"></param>
    ''' <param name="blnIncludeIMSFields"></param>
    ''' <param name="blnMSGFPlus"></param>
    ''' <remarks></remarks>
    Private Sub WriteSearchResultToFile(
      intResultID As Integer,
      swResultFile As StreamWriter,
      udtSearchResult As udtMSGFDBSearchResultType,
      ByRef strErrorLog As String,
      blnIncludeFDRandPepFDR As Boolean,
      blnIncludeEFDR As Boolean,
      blnIncludeIMSFields As Boolean,
      blnMSGFPlus As Boolean)

        Try

            ' Primary Columns (other columns are added in certain circumstances):
            '
            ' MSGFDB
            ' ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecProb    Rank_MSGFDB_SpecProb    PValue  FDR     PepFDR

            ' MSGF+
            ' ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecEValue  Rank_MSGFDB_SpecEValue  EValue  QValue  PepQValue  IsotopeError

            Dim lstData As New List(Of String)
            lstData.Add(intResultID.ToString())
            lstData.Add(udtSearchResult.Scan)
            lstData.Add(udtSearchResult.FragMethod)
            lstData.Add(udtSearchResult.SpecIndex)
            lstData.Add(udtSearchResult.Charge)
            lstData.Add(udtSearchResult.PrecursorMZ)
            lstData.Add(udtSearchResult.PMErrorDa)
            lstData.Add(udtSearchResult.PMErrorPPM)
            lstData.Add(udtSearchResult.MH)
            lstData.Add(udtSearchResult.Peptide)
            lstData.Add(udtSearchResult.Protein)
            lstData.Add(udtSearchResult.NTT)
            lstData.Add(udtSearchResult.DeNovoScore)
            lstData.Add(udtSearchResult.MSGFScore)
            lstData.Add(udtSearchResult.SpecProb)
            lstData.Add(udtSearchResult.RankSpecProb.ToString)
            lstData.Add(udtSearchResult.PValue)

            If blnIncludeFDRandPepFDR Then
                lstData.Add(udtSearchResult.FDR)
                lstData.Add(udtSearchResult.PepFDR)
            ElseIf blnIncludeEFDR Then
                lstData.Add(udtSearchResult.FDR)
                lstData.Add("1")
            End If

            If blnMSGFPlus Then
                lstData.Add(udtSearchResult.IsotopeError.ToString)
            End If

            If blnIncludeIMSFields Then
                lstData.Add(udtSearchResult.IMSScan.ToString)
                lstData.Add(udtSearchResult.IMSDriftTime)
            End If

            swResultFile.WriteLine(CollapseList(lstData))

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

#Region "IComparer Classes"

    Protected Class MSGFDBSearchResultsComparerScanChargeSpecProbPeptide
        Implements IComparer(Of udtMSGFDBSearchResultType)

        Public Function Compare(x As udtMSGFDBSearchResultType, y As udtMSGFDBSearchResultType) As Integer Implements IComparer(Of udtMSGFDBSearchResultType).Compare

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
                    ' Charge is the same; check SpecProb
                    If x.SpecProbNum > y.SpecProbNum Then
                        Return 1
                    ElseIf x.SpecProbNum < y.SpecProbNum Then
                        Return -1
                    Else
                        ' SpecProb is the same; check peptide
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

    Protected Class MSGFDBSearchResultsComparerSpecProbScanChargePeptide
        Implements IComparer(Of udtMSGFDBSearchResultType)

        Public Function Compare(x As udtMSGFDBSearchResultType, y As udtMSGFDBSearchResultType) As Integer Implements IComparer(Of udtMSGFDBSearchResultType).Compare

            If x.SpecProbNum > y.SpecProbNum Then
                Return 1
            ElseIf x.SpecProbNum < y.SpecProbNum Then
                Return -1
            Else
                ' SpecProbNum is the same; check scan number
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
