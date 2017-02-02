Option Strict On

' This class reads in an MSGF_DB results file (txt format) and creates 
' a tab-delimited text file with the data.  It will insert modification symbols
' into the peptide sequences for modified peptides.
'
' The modification definition information is determined from the MSGF+ parameter file
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

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "November 22, 2016"
        mModMassRegEx = New Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS)

        mPeptideCleavageStateCalculator = New clsPeptideCleavageStateCalculator()
        mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin)

        mNumericModErrors = 0
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MSGFDB_FILE As String = "_msgfdb"
    Public Const FILENAME_SUFFIX_MSGFPLUS_FILE As String = "_msgfplus"

    Public Const N_TERMINUS_SYMBOL_MSGFDB As String = "_."
    Public Const C_TERMINUS_SYMBOL_MSGFDB As String = "._"

    <Obsolete("Used by MSGF-DB; renamed to SpecEValue in MSGF+")>
    Public Const DEFAULT_SYN_FILE_MSGF_SPECPROB_THRESHOLD As Single = 0.0000005

    <Obsolete("Used by MSGF-DB; renamed to EValue in MSGF+")>
    Public Const DEFAULT_SYN_FILE_PVALUE_THRESHOLD As Single = 0.75

    ''' <summary>
    ''' Filter passing peptides have MSGFDB_SpecEValue less than 5E-7 Or EValue less than 0.75 or QValue less than 10%
    ''' This filter is also used by MSPathFinder
    ''' </summary>
    Public Const DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD As Single = 0.0000005

    ''' <summary>
    ''' Filter passing peptides have MSGFDB_SpecEValue less than 5E-7 Or EValue less than 0.75 or QValue less than 10%
    ''' This filter is also used by MSPathFinder
    ''' </summary>
    Public Const DEFAULT_SYN_FILE_EVALUE_THRESHOLD As Single = 0.75

    Private Const SEARCH_ENGINE_NAME As String = "MSGF+"

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

    ' These columns correspond to the tab-delimited file created directly by MSGF+
    Private Const MSGFDBResultsFileColCount As Integer = 20
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
    Private Const MSGFDBSynFileColCount As Integer = 23
    Public Enum eMSFDBSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        FragMethod = 2
        SpecIndex = 3
        Charge = 4
        PrecursorMZ = 5
        DelM = 6                            ' Precursor error, in Da; if the search used a tolerance less than 0.5 Da or less than 500 ppm, this value is computed from the DelMPPM value
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

    Private Enum eFilteredOutputFileTypeConstants As Integer
        SynFile = 0
        FHTFile = 1
    End Enum
#End Region

#Region "Structures"
    Private Structure udtMSGFDBSearchResultType

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
        Public Peptide As String                ' Peptide sequence, including prefix, suffix, and any mod symbols or mod masses
        Public Protein As String
        Public NTT As String
        Public DeNovoScore As String
        Public MSGFScore As String
        Public SpecEValue As String              ' Smaller values are better scores (e.g. 1E-9 is better than 1E-6); MSGF+ renamed this from SpecProb to SpecEValue
        Public SpecEValueNum As Double
        Public EValue As String                     ' Smaller values are better scores (e.g. 1E-7 is better than 1E-3); MSGF+ renamed this from PValue to EValue
        Public EValueNum As Double
        Public QValue As String                     ' Holds FDR when a target/decoy search was used; holds EFDR when a non-decoy search was used; holds QValue for MSGF+
        Public QValueNum As Double                  ' Numeric equivalent of QValue
        Public PepQValue As String                  ' Only used when target/decoy search was used; holds PepQValue for MSGF+
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
            SpecEValue = String.Empty
            SpecEValueNum = 0
            EValue = String.Empty
            EValueNum = 0
            QValue = String.Empty
            QValueNum = 0
            PepQValue = String.Empty
            RankSpecProb = 0
            IMSScan = 0
            IMSDriftTime = String.Empty
            IsotopeError = 0
        End Sub
    End Structure

    Private Structure udtScanGroupInfoType
        Public ScanGroupID As Integer
        Public Charge As Short
        Public Scan As Integer
    End Structure

    Private Structure udtTerminusCharsType
        Public NTerm As Char
        Public CTerm As Char
    End Structure

    Private Structure udtParentMassToleranceType
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

        Public Overrides Function ToString() As String
            Dim units As String
            Dim equivalenceThreshold As Double

            If IsPPM Then
                units = "ppm"
                equivalenceThreshold = 0.01
            Else
                units = "Da"
                equivalenceThreshold = 0.0001
            End If

            If Math.Abs(ToleranceLeft - ToleranceRight) < equivalenceThreshold Then
                Return "+/-" & ToleranceLeft & " " & units
            Else
                Return "-" & ToleranceRight & ", +" & ToleranceLeft & " " & units
            End If

        End Function
    End Structure

#End Region

#Region "Classwide Variables"
    Private ReadOnly mPeptideCleavageStateCalculator As clsPeptideCleavageStateCalculator

    Private mParentMassToleranceInfo As udtParentMassToleranceType

    Private mPrecursorMassErrorWarningCount As Integer

    ''' <summary>
    ''' Looks for numeric mods in MSGF+ results
    ''' For example, +14.016 in K.LQVPAGK+14.016ANPSPPIGPALGQR.G
    ''' </summary>
    ''' <remarks></remarks>
    Private ReadOnly mModMassRegEx As Regex

    Private mNumericModErrors As Integer

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

        strSequence = objSearchResult.PeptideSequenceWithMods

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
                            objSearchResult.SearchResultAddModification(objModificationDefinition, chChar, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                        End If
                    End If
                Next intModIndex
            ElseIf IsLetterAtoZ(chMostRecentLetter) Then
                blnSuccess = objSearchResult.SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                If Not blnSuccess Then
                    Dim strErrorMessage As String = objSearchResult.ErrorMessage
                    If String.IsNullOrEmpty(strErrorMessage) Then
                        strErrorMessage = "SearchResultAddDynamicModification returned false for symbol " & chChar
                    End If
                    SetErrorMessage(strErrorMessage & "; ResultID = " & objSearchResult.ResultID)
                End If
            Else
                ' We found a modification symbol but chMostRecentLetter is not a letter
                ' Therefore, this modification symbol is at the beginning of the string; ignore the symbol
            End If

        Next intIndex

    End Sub

    ''' <summary>
    ''' Adds or updates the prefix and suffix residues to the peptide, as defined in kvProteinInfo
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="kvProteinInfo"></param>
    ''' <returns>Peptide sequence with N-terminal and C-Terminal residues</returns>
    ''' <remarks></remarks>
    Private Function AddUpdatePrefixAndSuffixResidues(strPeptide As String, kvProteinInfo As KeyValuePair(Of String, udtTerminusCharsType)) As String

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

    Private Sub AppendToScanGroupDetails(
      lstScanGroupDetails As ICollection(Of udtScanGroupInfoType),
      htScanGroupCombo As IDictionary(Of String, Boolean),
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
       lstSearchResults As ICollection(Of udtMSGFDBSearchResultType),
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

        ' Prior to September 2014 ranks were assigned per charge state per scan; 
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

        Dim lstResultsBySpecProb = (From item In dctResultsSubset Select item Order By item.Value.SpecEValueNum).ToList()

        Dim dblLastValue As Double
        Dim intCurrentRank As Integer = -1

        For Each entry In lstResultsBySpecProb
            Dim currentResult = lstSearchResults(entry.Key)

            If intCurrentRank < 0 Then
                dblLastValue = currentResult.SpecEValueNum
                intCurrentRank = 1
            Else
                If Math.Abs(currentResult.SpecEValueNum - dblLastValue) > Double.Epsilon Then
                    dblLastValue = currentResult.SpecEValueNum
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

            ' Make sure .PeptideSequenceWithMods does not have any generic mod masses
            ' It should only have mod symbols
            Dim reMatch = mModMassRegEx.Match(objSearchResult.PeptideSequenceWithMods)
            If reMatch.Success Then
                ' Modification mass did not have a symbol associated with it in the _ModDefs.txt file
                ' We could try to handle this, listing the modification mass in place of the modification symbol in the _ModDetails.txt file, but will
                ' instead abort processing

                mNumericModErrors += 1

                If mNumericModErrors < 250 Then
                    Dim localErrorMessage = "Search result contains a numeric mod mass that could not be associated with a modification symbol; ResultID = " & objSearchResult.ResultID & ", ModMass = " & reMatch.Value.ToString
                    SetErrorMessage(localErrorMessage)
                ElseIf mNumericModErrors = 250 Then
                    SetErrorMessage("Too many numeric mod mass results have been found; suppressing further logging")
                End If

                Return False
            End If

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

    Private Function ComputeCleaveageState(strSequenceWithMods As String) As Short

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
    Private Function ComputeDelMCorrectedPPM(
      dblPrecursorErrorDa As Double,
      dblPrecursorMZ As Double,
      intCharge As Integer,
      dblPeptideMonoisotopicMass As Double,
      blnAdjustPrecursorMassForC13 As Boolean) As Double

        Dim dblPeptideDeltaMassCorrectedPpm As Double

        Dim dblPrecursorMonoMass As Double

        ' Compute the original value for the precursor monoisotopic mass
        dblPrecursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, intCharge, 0)

        dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, dblPeptideMonoisotopicMass)

        Return dblPeptideDeltaMassCorrectedPpm

    End Function

    ''' <summary>
    ''' Compute the monoisotopic mass of the peptide
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="dblTotalModMass"></param>
    ''' <returns></returns>
    Private Function ComputePeptideMass(strPeptide As String, dblTotalModMass As Double) As Double

        Dim strCleanSequence As String
        Dim dblMass As Double

        strCleanSequence = GetCleanSequence(strPeptide)

        dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence)
        dblMass += dblTotalModMass

        Return dblMass

    End Function

    Protected Overrides Function ConstructPepToProteinMapFilePath(strInputFilePath As String, strOutputFolderPath As String, MTS As Boolean) As String

        Dim strPepToProteinMapFilePath As String = Path.GetFileNameWithoutExtension(strInputFilePath)

        If strPepToProteinMapFilePath.ToLower().EndsWith("_msgfplus_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_msgfplus_fht") OrElse
           strPepToProteinMapFilePath.ToLower().EndsWith("_msgfdb_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_msgfdb_fht") Then
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
    Private Function ConvertMGSFModMassesToSymbols(
      currentResidue As String,
      strModDigits As String,
      <Out()> ByRef strModSymbols As String,
      lstMSGFDBModInfo As IReadOnlyList(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      blnNterminalMod As Boolean,
      blnPossibleCTerminalMod As Boolean,
      <Out()> ByRef dblModMassFound As Double,
      <Out()> ByRef blnIsStaticMod As Boolean) As Boolean

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

        reMatches = mModMassRegEx.Matches(strModDigits)

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
    ''' This routine creates a first hits file or synopsis file from the output from MSGF+
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
    Private Function CreateFHTorSYNResultsFile(
      strInputFilePath As String,
      strOutputFilePath As String,
      strScanGroupFilePath As String,
      lstMSGFDBModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      <Out()> ByRef blnMSGFPlus As Boolean,
      lstSpecIdToIndex As IDictionary(Of String, Integer),
      eFilteredOutputFileType As eFilteredOutputFileTypeConstants) As Boolean

        Dim strLineIn As String

        Dim lstSearchResultsCurrentScan As New List(Of udtMSGFDBSearchResultType)
        Dim lstSearchResultsPrefiltered As New List(Of udtMSGFDBSearchResultType)

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
                Dim empiricalFormulaString = customAADef.ModMass
                Dim aminoAcidMass = customAADef.ModMassVal

                Try
                    Dim elementalComposition = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormulaString)

                    mPeptideSeqMassCalculator.SetAminoAcidMass(aminoAcidSymbol, aminoAcidMass)
                    mPeptideSeqMassCalculator.SetAminoAcidAtomCounts(aminoAcidSymbol, elementalComposition)

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

                    ' Initialize a dictionary that tracks the peptide sequence for each combo of scan and charge
                    ' Keys are Scan_Charge, values track the clean sequence, the associated protein name, and the protein number for that name
                    ' Note that we can only track protein numbers if the FASTA file path was provided at the command line
                    Dim scanChargeFirstHit = New Dictionary(Of String, clsFirstHitInfo)

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

                        If blnValidSearchResult AndAlso lstSearchResultsCurrentScan.Count > 0 Then
                            If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
                                ' Synopsis file
                                blnValidSearchResult = MSGFPlusResultPassesSynFilter(lstSearchResultsCurrentScan(0))
                            Else

                                ' First Hits file
                                Dim scanChargeKey = lstSearchResultsCurrentScan(0).Scan & "_" & lstSearchResultsCurrentScan(0).Charge
                                Dim firstHitPeptide As clsFirstHitInfo = Nothing

                                If scanChargeFirstHit.TryGetValue(scanChargeKey, firstHitPeptide) Then
                                    ' A result has already been stored for this scan/charge combo
                                    blnValidSearchResult = False

                                    ' Possibly update the associated protein name
                                    Dim strNewPrefix As String = Nothing
                                    Dim strNewSuffix As String = Nothing
                                    If firstHitPeptide.CleanSequence.Equals(GetCleanSequence(lstSearchResultsCurrentScan(0).Peptide, strNewPrefix, strNewSuffix)) Then
                                        Dim bestProtein = GetBestProteinName(firstHitPeptide.ProteinName, firstHitPeptide.ProteinNumber, lstSearchResultsCurrentScan(0).Protein)
                                        If bestProtein.Value < firstHitPeptide.ProteinNumber Then
                                            firstHitPeptide.ProteinName = bestProtein.Key
                                            firstHitPeptide.ProteinNumber = bestProtein.Value
                                            firstHitPeptide.UpdatePrefixAndSuffix(strNewPrefix, strNewSuffix)
                                        End If
                                    End If

                                Else
                                    firstHitPeptide = New clsFirstHitInfo(lstSearchResultsCurrentScan(0).Peptide, GetCleanSequence(lstSearchResultsCurrentScan(0).Peptide)) With {
                                        .ProteinName = lstSearchResultsCurrentScan(0).Protein,
                                        .ProteinNumber = Int32.MaxValue
                                    }

                                    Dim proteinNumber As Integer
                                    If mProteinNameOrder.TryGetValue(lstSearchResultsCurrentScan(0).Protein, proteinNumber) Then
                                        firstHitPeptide.ProteinNumber = proteinNumber
                                    End If

                                    scanChargeFirstHit.Add(scanChargeKey, firstHitPeptide)
                                End If
                            End If

                            If blnValidSearchResult Then
                                ExpandListIfRequired(lstSearchResultsPrefiltered, lstSearchResultsCurrentScan.Count)
                                lstSearchResultsPrefiltered.AddRange(lstSearchResultsCurrentScan)
                            End If

                        End If

                        ' Update the progress
                        sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
                        If CreateProteinModsFile Then
                            sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
                        End If
                        UpdateProgress(sngPercentComplete)

                    Loop

                    lstSearchResultsPrefiltered.TrimExcess()

                    ' Sort the SearchResults by scan, charge, and ascending SpecProb
                    lstSearchResultsPrefiltered.Sort(New MSGFDBSearchResultsComparerScanChargeSpecProbPeptide)

                    If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.FHTFile Then

                        ' Update the protein names in lstSearchResultsPrefiltered using scanChargeFirstHit
                        ' This step is likely unnecessary, thus the "Unexpected code reached" message below
                        For intIndex = 0 To lstSearchResultsPrefiltered.Count - 1

                            Dim scanChargeKey = lstSearchResultsPrefiltered(intIndex).Scan & "_" & lstSearchResultsPrefiltered(intIndex).Charge
                            Dim firstHitPeptide As clsFirstHitInfo = Nothing

                            If scanChargeFirstHit.TryGetValue(scanChargeKey, firstHitPeptide) Then
                                If Not lstSearchResultsPrefiltered(intIndex).Protein.Equals(firstHitPeptide.ProteinName) Then
                                    Console.WriteLine("Unexpected code reached; possible logic error in clsMSGFDBResultsProcessor.CreateFHTorSYNCResultsFile")

                                    If firstHitPeptide.CleanSequence.Equals(GetCleanSequence(lstSearchResultsPrefiltered(intIndex).Peptide)) Then
                                        Dim updatedSearchResult = lstSearchResultsPrefiltered(intIndex)
                                        updatedSearchResult.Peptide = firstHitPeptide.SequenceWithModsAndContext
                                        updatedSearchResult.Protein = String.Copy(firstHitPeptide.ProteinName)
                                        lstSearchResultsPrefiltered(intIndex) = updatedSearchResult
                                    Else
                                        Console.WriteLine(String.Format("Possible programming bug; " &
                                                                        "mix of peptides tracked for a given scan/charge combo when caching data for First Hits files; " &
                                                                        "see scan_charge {0}", scanChargeKey))
                                    End If

                                End If
                            End If

                        Next

                    End If

                    ' Now filter the data and store in lstFilteredSearchResults
                    ' Due to code updates in October 2016, lstSearchResultsPrefiltered already has filtered data
                    Dim intStartIndex = 0
                    Dim intEndIndex As Integer

                    Do While intStartIndex < lstSearchResultsPrefiltered.Count
                        intEndIndex = intStartIndex
                        ' Find all of the peptides with the same scan number
                        Do While intEndIndex + 1 < lstSearchResultsPrefiltered.Count AndAlso
                            lstSearchResultsPrefiltered(intEndIndex + 1).ScanNum = lstSearchResultsPrefiltered(intStartIndex).ScanNum
                            intEndIndex += 1
                        Loop

                        ' Store the results for this scan
                        If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
                            StoreSynMatches(lstSearchResultsPrefiltered, intStartIndex, intEndIndex, lstFilteredSearchResults)
                        Else
                            StoreTopFHTMatch(lstSearchResultsPrefiltered, intStartIndex, intEndIndex, lstFilteredSearchResults)
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
    ''' Extracts mod info from either a MSGF+ param file or from a MSGFPlus_Mods.txt file (previously MSGFDB_Mods.txt)
    ''' </summary>
    ''' <param name="strMSGFDBParamFilePath"></param>
    ''' <param name="lstModInfo"></param>
    ''' <returns>True if success; false if a problem</returns>
    ''' <remarks></remarks>
    Private Function ExtractModInfoFromParamFile(
       strMSGFDBParamFilePath As String,
       <Out()> ByRef lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Boolean

        Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME)

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
    ''' Extracts parent mass tolerance from the parameters loaded from an MSGF+ parameter file
    ''' </summary>
    ''' <param name="objSearchEngineParams"></param>	
    ''' <returns>Parent mass tolerance info.  Tolerances will be 0 if an error occurs</returns>
    ''' <remarks></remarks>
    Private Function ExtractParentMassToleranceFromParamFile(objSearchEngineParams As clsSearchEngineParameters) As udtParentMassToleranceType

        Const PM_TOLERANCE_TAG = "PMTolerance"

        Dim udtParentMassToleranceInfo = New udtParentMassToleranceType()

        Try
            udtParentMassToleranceInfo.Clear()

            Dim strValue As String = String.Empty
            If objSearchEngineParams.Parameters.TryGetValue(PM_TOLERANCE_TAG, strValue) Then

                ' Parent ion tolerance line found

                ' Split the line on commas
                Dim strSplitLine = strValue.Split(","c)

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
            End If

            Console.WriteLine()

        Catch ex As Exception
            SetErrorMessage(String.Format("Error parsing the ParentMass tolerance from the MSGF+ parameter file ({0}): {1}",
                                          Path.GetFileName(objSearchEngineParams.SearchEngineParamFilePath), ex.Message))
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
        End Try

        Return udtParentMassToleranceInfo

    End Function

    ''' <summary>
    ''' Look for candidateProteinName in mProteinNameOrder
    ''' If found, and if the proteinNumber for that protein is less than currentProteinNumber, 
    ''' return a KeyValuePair with candidateProteinName and the proteinNumber for that protein
    ''' Otherwise, return a KeyValuePair with currentProteinName and currentProteinNumber
    ''' </summary>
    ''' <param name="currentProteinName"></param>
    ''' <param name="currentProteinNumber"></param>
    ''' <param name="candidateProteinName"></param>
    ''' <returns></returns>
    Private Function GetBestProteinName(currentProteinName As String, currentProteinNumber As Integer, candidateProteinName As String) As KeyValuePair(Of String, Integer)

        If mProteinNameOrder.Count > 0 Then

            ' Lookup the protein number (to make sure we use the protein name that occurs first in the FASTA file)
            ' Only possible if the user provided the path to the FASTA file
            Dim proteinNumber As Integer

            If mProteinNameOrder.TryGetValue(candidateProteinName, proteinNumber) Then
                If proteinNumber < currentProteinNumber Then
                    ' A better protein name has been found (or this is the first protein and we just determined the protein number to associate with it)
                    Return New KeyValuePair(Of String, Integer)(candidateProteinName, proteinNumber)
                End If
            Else
                ' Protein not found in mProteinNameOrder
                ' It's likely a reverse-hit protein
            End If
        End If

        ' A better protein name was not found; return the current info
        Return New KeyValuePair(Of String, Integer)(currentProteinName, currentProteinNumber)

    End Function

    ''' <summary>
    ''' Load the PeptideToProteinMap information; in addition, creates the _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
    ''' </summary>
    ''' <param name="strPepToProteinMapFilePath"></param>
    ''' <param name="strOutputFolderPath"></param>
    ''' <param name="lstMSGFDBModInfo"></param>
    ''' <param name="blnMSGFPlus">Should be set to True if processing MSGF+ results</param>
    ''' <param name="lstPepToProteinMapping"></param>
    ''' <param name="strMTSPepToProteinMapFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function LoadPeptideToProteinMapInfoMSGFDB(
      strPepToProteinMapFilePath As String,
      strOutputFolderPath As String,
      lstMSGFDBModInfo As IReadOnlyList(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
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
            End If

            If Not File.Exists(strPepToProteinMapFilePath) Then
                Dim strPepToProteinMapAlternate = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strPepToProteinMapFilePath, "Dataset_msgfdb.txt")
                If File.Exists(strPepToProteinMapAlternate) Then
                    strPepToProteinMapFilePath = strPepToProteinMapAlternate
                End If
            End If

            If Not File.Exists(strPepToProteinMapFilePath) Then
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

    ''' <summary>
    ''' Load the MSGF+ parameter file and updates settings
    ''' </summary>
    ''' <param name="msgfPlusParamFilePath"></param>
    ''' <returns>
    ''' True if success, false if an error.  
    ''' Returns True if msgfPlusParamFilePath is empty
    ''' Returns False if the paramFilePath is defined but the file is not found or cannot be parsed</returns>
    Private Function LoadSearchEngineParamFile(msgfPlusParamFilePath As String) As Boolean

        If String.IsNullOrWhiteSpace(msgfPlusParamFilePath) Then
            ReportWarning("MSGF+ parameter file is not defined. Unable to extract parent mass tolerance info or custom charge carrier masses")
            Return True
        End If

        Dim objSearchEngineParams = New clsSearchEngineParameters(SEARCH_ENGINE_NAME)

        Dim localErrorMessage As String = Nothing
        Dim localWarningMessage As String = Nothing

        Dim success = clsPHRPParser.ReadKeyValuePairSearchEngineParamFile(
            SEARCH_ENGINE_NAME, msgfPlusParamFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB,
            objSearchEngineParams, localErrorMessage, localWarningMessage)

        If Not String.IsNullOrWhiteSpace(localErrorMessage) Then
            ReportError(localErrorMessage)
            Return False
        End If

        If Not String.IsNullOrWhiteSpace(localWarningMessage) Then
            ReportWarning(localWarningMessage)
        End If

        If objSearchEngineParams Is Nothing OrElse objSearchEngineParams.Parameters.Count = 0 Then
            SetErrorMessage("MSGF+ parameter file is empty; unable to extract parent mass tolerance info")
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            Return False
        End If

        ' Parse the PMTolerance setting
        mParentMassToleranceInfo = ExtractParentMassToleranceFromParamFile(objSearchEngineParams)

        ' Parse the ChargeCarrierMass setting
        Dim customChargeCarrierMass As Double
        If clsPHRPParserMSGFDB.GetCustomChargeCarrierMass(objSearchEngineParams, customChargeCarrierMass) Then
            ReportMessage(String.Format("Using a charge carrier mass of {0:F3} Da", customChargeCarrierMass))
            mPeptideSeqMassCalculator.ChargeCarrierMass = customChargeCarrierMass
        End If

        Return success

    End Function

    Private Sub ModExtractorErrorHandler(errMsg As String)
        SetErrorMessage(errMsg)
        SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
    End Sub

    Private Sub ModExtractorWarningHandler(warningMsg As String)
        ReportWarning(warningMsg)
    End Sub

    Private Function MSGFPlusResultPassesSynFilter(udtMsgfdbSearchResultType As udtMSGFDBSearchResultType) As Boolean
        If udtMsgfdbSearchResultType.EValueNum <= MSGFDBSynopsisFileEValueThreshold OrElse
               udtMsgfdbSearchResultType.SpecEValueNum <= MSGFDBSynopsisFileSpecEValueThreshold OrElse
               udtMsgfdbSearchResultType.QValueNum > 0 AndAlso udtMsgfdbSearchResultType.QValueNum < 0.01 Then
            Return True
        End If

        Return False
    End Function

    Private Function ParseMSGFDBSynopsisFile(
      strInputFilePath As String,
      strOutputFolderPath As String,
      lstPepToProteinMapping As List(Of udtPepToProteinMappingType),
      blnResetMassCorrectionTagsAndModificationDefinitions As Boolean) As Boolean

        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim strPreviousSpecProb As String

        ' Note that MSGF+ synopsis files are normally sorted on SpecProb value, ascending
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
        '  we will keep track of the scan, charge, and peptide information parsed for each unique SpecProb encountered
        ' Although this was a possiblity with Inspect, it likely never occurs for MSGF+
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

        Dim successOverall = True

        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            mNumericModErrors = 0

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
                objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange)

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

                            If objSearchResult.SpecEValue = strPreviousSpecProb Then
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
                                strPreviousSpecProb = objSearchResult.SpecEValue

                                ' Append a new entry to htPeptidesFoundForSpecProbLevel
                                htPeptidesFoundForSpecProbLevel.Add(strKey, 1)
                                blnFirstMatchForGroup = True
                            End If

                            blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup)
                            If Not blnSuccess Then
                                successOverall = False
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
                        If CreateProteinModsFile Then
                            sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
                        End If
                        UpdateProgress(sngPercentComplete)

                        intResultsProcessed += 1

                    Loop

                End Using

                If CreateModificationSummaryFile Then
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

                blnSuccess = successOverall

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
      lstMSGFDBModInfo As IReadOnlyList(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      lstSearchResultsCurrentScan As ICollection(Of udtMSGFDBSearchResultType),
      ByRef strErrorLog As String,
      intColumnMapping() As Integer,
      ByRef intNextScanGroupID As Integer,
      lstScanGroupDetails As ICollection(Of udtScanGroupInfoType),
      htScanGroupCombo As IDictionary(Of String, Boolean),
      lstSpecIdToIndex As IDictionary(Of String, Integer)) As Boolean

        ' Parses an entry from the MSGF+ results file

        Dim udtSearchResult = New udtMSGFDBSearchResultType
        Dim strSplitLine() As String = Nothing

        Dim intScanCount As Integer
        Dim strSplitResult() As String

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

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpectrumFile), udtSearchResult.SpectrumFileName) Then
                    ReportError("SpectrumFile column is missing or invalid", True)
                End If
                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpecIndex), udtSearchResult.SpecIndex)

                If blnMSGFPlus Then
                    Dim intSpecIndex As Integer
                    Dim blnGenerateSpecIndex = True

                    If Not Integer.TryParse(udtSearchResult.SpecIndex, intSpecIndex) Then
                        ' MSGF+ includes text in the SpecID column, for example: "controllerType=0 controllerNumber=1 scan=6390" or "index=4323"
                        ' Need to convert these to an integer

                        If udtSearchResult.SpecIndex.StartsWith("index=") Then
                            udtSearchResult.SpecIndex = udtSearchResult.SpecIndex.Substring("index=".Length)
                            If Integer.TryParse(udtSearchResult.SpecIndex, intSpecIndex) Then
                                blnGenerateSpecIndex = False
                            End If
                        End If

                        If blnGenerateSpecIndex Then
                            If Not lstSpecIdToIndex.TryGetValue(udtSearchResult.SpecIndex, intSpecIndex) Then
                                intSpecIndex = lstSpecIdToIndex.Count + 1
                                lstSpecIdToIndex.Add(udtSearchResult.SpecIndex, intSpecIndex)
                            End If

                            udtSearchResult.SpecIndex = intSpecIndex.ToString()
                        End If

                    End If
                End If

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Scan), udtSearchResult.Scan) Then
                    ReportError("Scan column is missing or invalid", True)
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.FragMethod), udtSearchResult.FragMethod)

                intSlashIndex = udtSearchResult.Scan.IndexOf("/"c)
                If intSlashIndex > 0 Then
                    ' This is a merged spectrum and thus scan number looks like: 3010/3011/3012
                    ' Split the Scan list on the slash
                    ' Later in this function, we'll append lstSearchResults with this scan plus the other scans

                    strSplitResult = udtSearchResult.Scan.Split("/"c)
                    intScanCount = strSplitResult.Length
                    ReDim udtMergedScanInfo(intScanCount - 1)

                    For intIndex = 0 To intScanCount - 1
                        udtMergedScanInfo(intIndex) = New udtMSGFDBSearchResultType
                        udtMergedScanInfo(intIndex).Clear()
                        udtMergedScanInfo(intIndex).Scan = strSplitResult(intIndex)
                        udtMergedScanInfo(intIndex).ScanNum = CIntSafe(strSplitResult(intIndex), 0)
                    Next

                    ' Now split SpecIndex and store in udtMergedScanInfo
                    strSplitResult = udtSearchResult.SpecIndex.Split("/"c)

                    For intIndex = 0 To strSplitResult.Length - 1
                        If intIndex >= udtMergedScanInfo.Length Then
                            ' There are more entries for SpecIndex than there are for Scan#; this is unexpected
                            Exit For
                        End If
                        udtMergedScanInfo(intIndex).SpecIndex = strSplitResult(intIndex)
                    Next

                    ' Now split FragMethod and store in udtMergedScanInfo
                    strSplitResult = udtSearchResult.FragMethod.Split("/"c)

                    For intIndex = 0 To strSplitResult.Length - 1
                        If intIndex >= udtMergedScanInfo.Length Then
                            ' There are more entries for FragMethod than there are for Scan#; this is unexpected
                            Exit For
                        End If
                        udtMergedScanInfo(intIndex).FragMethod = strSplitResult(intIndex)
                    Next

                Else
                    udtSearchResult.ScanNum = CIntSafe(udtSearchResult.Scan, 0)
                    intScanCount = 1
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PrecursorMZ), udtSearchResult.PrecursorMZ)

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Charge), udtSearchResult.Charge)
                udtSearchResult.ChargeNum = CShort(CIntSafe(udtSearchResult.Charge, 0))

                ' Precursor mass error could be in PPM or Da
                '   In MSGFDB, the header line will have PMError(ppm)        or PMError(Da)
                '   In MSGF+,  the header line will have PrecursorError(ppm) or PrecursorError(Da)
                Dim dblPrecursorErrorDa As Double

                If intColumnMapping(eMSGFDBResultsFileColumns.PMErrorPPM) >= 0 Then
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PMErrorPPM), udtSearchResult.PMErrorPPM)
                Else
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PMErrorDa), udtSearchResult.PMErrorDa)
                    dblPrecursorErrorDa = CDblSafe(udtSearchResult.PMErrorDa, 0)
                    udtSearchResult.PMErrorPPM = String.Empty              ' We'll populate this column later in this function
                End If


                If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Peptide), udtSearchResult.Peptide) Then
                    ReportError("Peptide column is missing or invalid", True)
                End If


                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Protein), udtSearchResult.Protein)

                ' MSGF+ .tsv files may have a semicolon separated list of protein names; check for this
                udtSearchResult.Protein = SplitProteinList(udtSearchResult.Protein, lstProteinInfo)

                If lstProteinInfo.Count > 0 Then
                    ' Need to add the prefix and suffix residues
                    udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, lstProteinInfo.First)
                End If

                ' Replace any mod text values in the peptide sequence with the appropriate mod symbols
                ' In addition, replace the terminus symbols with dashes
                Dim dblTotalModMass As Double
                udtSearchResult.Peptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(udtSearchResult.Peptide), lstMSGFDBModInfo, blnMSGFPlus, dblTotalModMass)

                ' Compute monoisotopic mass of the peptide
                Dim dblPeptideMonoisotopicMass = ComputePeptideMass(udtSearchResult.Peptide, dblTotalModMass)

                ' Store the monoisotopic MH value in .MH
                ' This is (M+H)+ when the charge carrier is a proton
                udtSearchResult.MH = NumToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPeptideMonoisotopicMass, 0, 1), 6, True)

                If Not String.IsNullOrEmpty(udtSearchResult.PMErrorPPM) Then

                    ' Convert the ppm-based PM Error to Da-based

                    Dim dblPMErrorPPM As Double
                    Dim dblPrecursorMZ As Double

                    If Double.TryParse(udtSearchResult.PrecursorMZ, dblPrecursorMZ) Then
                        ' Note that since .PMErrorPPM is present, the Precursor m/z is a C13-corrected m/z value
                        ' In other words, it may not be the actual m/z selected for fragmentation.

                        If Double.TryParse(udtSearchResult.PMErrorPPM, dblPMErrorPPM) Then

                            If mParentMassToleranceInfo.IsPPM AndAlso
                                (dblPMErrorPPM < -mParentMassToleranceInfo.ToleranceLeft * 1.5 OrElse
                                 dblPMErrorPPM > mParentMassToleranceInfo.ToleranceRight * 1.5) Then

                                ' PPM error computed by MSGF+ is more than 1.5-fold larger than the ppm-based parent ion tolerance; don't trust the value computed by MSGF+

                                mPrecursorMassErrorWarningCount += 1
                                If mPrecursorMassErrorWarningCount <= 10 Then
                                    ReportWarning("Precursor mass error computed by MSGF+ is 1.5-fold larger than search tolerance: " & udtSearchResult.PMErrorPPM & " vs. " & mParentMassToleranceInfo.ToleranceLeft.ToString("0") & "ppm," & mParentMassToleranceInfo.ToleranceRight.ToString("0") & "ppm")
                                    If mPrecursorMassErrorWarningCount = 10 Then
                                        ReportWarning("Additional mass errors will not be reported")
                                    End If
                                End If


                                Dim dblPrecursorMonoMass As Double
                                dblPrecursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, udtSearchResult.ChargeNum, 0)

                                dblPrecursorErrorDa = dblPrecursorMonoMass - dblPeptideMonoisotopicMass

                                udtSearchResult.PMErrorPPM = String.Empty

                            Else

                                dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPMErrorPPM, dblPeptideMonoisotopicMass)

                                ' Note that this will be a C13-corrected precursor error; not the true precursor error
                                udtSearchResult.PMErrorDa = NumToString(dblPrecursorErrorDa, 6, True)
                            End If

                        End If
                    End If
                End If

                If String.IsNullOrEmpty(udtSearchResult.PMErrorPPM) Then

                    Dim dblPrecursorMZ As Double
                    If Double.TryParse(udtSearchResult.PrecursorMZ, dblPrecursorMZ) Then
                        Dim dblPeptideDeltaMassCorrectedPpm As Double

                        dblPeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMZ,
                            udtSearchResult.ChargeNum, dblPeptideMonoisotopicMass, True)

                        udtSearchResult.PMErrorPPM = NumToString(dblPeptideDeltaMassCorrectedPpm, 5, True)

                        If String.IsNullOrEmpty(udtSearchResult.PMErrorDa) Then
                            dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPeptideDeltaMassCorrectedPpm, dblPeptideMonoisotopicMass)

                            ' Note that this will be a C13-corrected precursor error; not the true precursor error
                            udtSearchResult.PMErrorDa = NumToString(dblPrecursorErrorDa, 6, True)
                        End If

                    End If

                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.DeNovoScore), udtSearchResult.DeNovoScore)
                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.MSGFScore), udtSearchResult.MSGFScore)
                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpecProb_EValue), udtSearchResult.SpecEValue)
                If Not Double.TryParse(udtSearchResult.SpecEValue, udtSearchResult.SpecEValueNum) Then udtSearchResult.SpecEValueNum = 0

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PValue_EValue), udtSearchResult.EValue)
                If Not Double.TryParse(udtSearchResult.EValue, udtSearchResult.EValueNum) Then udtSearchResult.EValueNum = 0

                blnTargetDecoyFDRValid = GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.FDR_QValue), udtSearchResult.QValue)
                If Not Double.TryParse(udtSearchResult.QValue, udtSearchResult.QValueNum) Then udtSearchResult.QValueNum = 0

                If blnTargetDecoyFDRValid Then
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PepFDR_PepQValue), udtSearchResult.PepQValue)
                Else
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.EFDR), udtSearchResult.QValue)
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.IsotopeError), udtSearchResult.IsotopeError)

                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.IMSScan), udtSearchResult.IMSScan)
                GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.IMSDriftTime), udtSearchResult.IMSDriftTime)

                udtSearchResult.NTT = ComputeCleaveageState(udtSearchResult.Peptide).ToString()

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
            ' Error parsing this row from the MSGF+ results file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MSGF+ Results in ParseMSGFDBResultsFileEntry for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MSGF+ Results in ParseMSGFDBResultsFileEntry" & ControlChars.NewLine
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
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue, eMSFDBSynFileColumns.SpecProb_EValue)

        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecProb, eMSFDBSynFileColumns.RankSpecProb)
        lstColumnNames.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue, eMSFDBSynFileColumns.RankSpecProb)

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

                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.DelM), .MSGFPlusComputedDelM)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.DelMPPM), .MSGFPlusComputedDelMPPM)

                    .PeptideDeltaMass = .MSGFPlusComputedDelM

                    ' Note: .PeptideDeltaMass is stored in the MSGF+ results file as "Observed_Mass - Theoretical_Mass"
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
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.SpecProb_EValue), .SpecEValue)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.RankSpecProb), .RankSpecEValue)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PValue_EValue), .EValue)

                    blnTargetDecoyFDRValid = GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.FDR_QValue), .QValue)
                    If blnTargetDecoyFDRValid Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PepFDR_PepQValue), .PepQValue)
                    Else
                        GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.EFDR), .QValue)
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
                            .ParentIonMH = NumToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, intCharge, 1), 6, True)
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

    Private Function ParseParentMassTolerance(strToleranceText As String, <Out()> ByRef dblTolerance As Double, <Out()> ByRef blnIsPPM As Boolean) As Boolean
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
        Dim strSynOutputFilePath As String
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

                ' Load the MSGF+ Parameter File so that we can determine the modification names and masses
                ' If the MSGFPlus_Mods.txt or MSGFDB_Mods.txt file was defined, the mod symbols in that file will be used to define the mod symbols in lstMSGFDBModInfo 
                Dim success = ExtractModInfoFromParamFile(SearchToolParameterFilePath, lstMSGFDBModInfo)
                If Not success Then
                    Return False
                End If

                If Not LoadSearchEngineParamFile(SearchToolParameterFilePath) Then
                    Return False
                End If


                Dim query = From item In lstMSGFDBModInfo Where item.ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA
                If query.Any() Then
                    ' Custom amino acids are defined; read their values and update the mass calculator

                    Dim localErrorMsg As String = String.Empty
                    Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME)

                    AddHandler modFileProcessor.ErrorOccurred, AddressOf ModExtractorErrorHandler
                    AddHandler modFileProcessor.WarningMessageEvent, AddressOf ModExtractorWarningHandler

                    clsPHRPParserMSGFDB.UpdateMassCalculatorMasses(SearchToolParameterFilePath, modFileProcessor, mPeptideSeqMassCalculator, localErrorMsg)

                    If Not String.IsNullOrWhiteSpace(localErrorMsg) AndAlso String.IsNullOrWhiteSpace(mErrorMessage) Then
                        ReportError(localErrorMsg)
                    End If

                End If

                ' Define the base output filename using strInputFilePath
                Dim strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath)

                ' Auto-replace "_msgfdb" with "_msgfplus"
                If strBaseName.ToLower().EndsWith("_msgfdb") Then
                    strBaseName = strBaseName.Substring(0, strBaseName.Length - "_msgfdb".Length) & "_msgfplus"
                End If

                If MyBase.CreateInspectFirstHitsFile Then

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

                If MyBase.CreateInspectSynopsisFile Then

                    ' Create the synopsis output file
                    MyBase.ResetProgress("Creating the SYN file", True)

                    ' The synopsis file name will be of the form BasePath_msgfplus_syn.txt
                    strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                    strScanGroupFilePath = Path.Combine(strOutputFolderPath, strBaseName & "_ScanGroupInfo.txt")

                    blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, strScanGroupFilePath, lstMSGFDBModInfo, blnMSGFPlus, lstSpecIdToIndex, eFilteredOutputFileTypeConstants.SynFile)

                    ' Load the PeptideToProteinMap information; if the file doesn't exist, a warning will be displayed, but processing will continue
                    ' LoadPeptideToProteinMapInfoMSGFDB also creates _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols							
                    strPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(Path.Combine(strOutputFolderPath, strBaseName) & ".txt", strOutputFolderPath, MTS:=False)

                    MyBase.ResetProgress("Loading the PepToProtein map file: " & Path.GetFileName(strPepToProteinMapFilePath), True)

                    LoadPeptideToProteinMapInfoMSGFDB(strPepToProteinMapFilePath, strOutputFolderPath, lstMSGFDBModInfo, blnMSGFPlus, lstPepToProteinMapping, strMTSPepToProteinMapFilePath)

                    ' Create the other PHRP-specific files
                    MyBase.ResetProgress("Creating the PHRP files for " & Path.GetFileName(strSynOutputFilePath), True)

                    ' Now parse the _syn.txt file that we just created to create the other PHRP files
                    blnSuccess = ParseMSGFDBSynopsisFile(strSynOutputFilePath, strOutputFolderPath, lstPepToProteinMapping, False)

                    ' Remove all items from lstPepToProteinMapping to reduce memory overhead
                    lstPepToProteinMapping.Clear()
                    lstPepToProteinMapping.TrimExcess()

                    If blnSuccess AndAlso CreateProteinModsFile Then
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
                If File.Exists(strMTSPepToProteinMapFilePath) AndAlso UseExistingMTSPepToProteinMapFile Then
                    blnSuccess = True
                Else
                    ' Auto-change mIgnorePeptideToProteinMapperErrors to True
                    ' We only do this for MSGFDB since it often includes reverse protein peptides in the results even though the FASTA file often does not have reverse proteins
                    IgnorePeptideToProteinMapperErrors = True
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
    Private Function ReplaceMSGFModTextWithSymbol(
      strPeptide As String,
      lstMSGFDBModInfo As IReadOnlyList(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
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

    Private Function ReplaceMSGFModTextWithMatchedSymbol(
      strPeptide As String,
      reGroup As Capture,
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

    Private Function ReplaceTerminus(strPeptide As String) As String

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
    Private Function SplitProteinList(strProteinList As String, lstProteinInfo As IDictionary(Of String, udtTerminusCharsType)) As String

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
      swResultFile As TextWriter,
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
    ''' <param name="lstFilteredSearchResults">Filtered search results</param>
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
        Dim currentPeptide = GetCleanSequence(udtCurrentResult.Peptide)

        For intIndex = intStartIndex To intEndIndex
            If intCurrentCharge <> lstSearchResults(intIndex).ChargeNum Then
                ' New charge state
                ' Store udtCurrentResult (from the previous charge state)
                lstFilteredSearchResults.Add(udtCurrentResult)

                udtCurrentResult = lstSearchResults(intIndex)
                intCurrentCharge = udtCurrentResult.ChargeNum
                currentProteinNumber = Int32.MaxValue
                currentPeptide = GetCleanSequence(udtCurrentResult.Peptide)
            End If

            Dim newPeptide = GetCleanSequence(lstSearchResults(intIndex).Peptide)
            If currentPeptide.Equals(newPeptide) Then
                Dim bestProtein = GetBestProteinName(udtCurrentResult.Protein, currentProteinNumber, lstSearchResults(intIndex).Protein)
                If bestProtein.Value < currentProteinNumber Then
                    currentProteinNumber = bestProtein.Value
                    If Not udtCurrentResult.Protein.Equals(bestProtein.Key) Then
                        udtCurrentResult.Protein = bestProtein.Key
                    End If
                End If
            End If

        Next intIndex

        ' Store udtCurrentResult (from the previous charge state)
        lstFilteredSearchResults.Add(udtCurrentResult)

    End Sub

    ''' <summary>
    ''' Stores the synopsis file matches for a single scan
    ''' </summary>
    ''' <param name="lstSearchResults">Search results</param>
    ''' <param name="intStartIndex">Start index for data in this scan</param>
    ''' <param name="intEndIndex">End index for data in this scan</param>
    ''' <param name="lstFilteredSearchResults">Filtered search results</param>
    ''' <remarks></remarks>
    Private Sub StoreSynMatches(
      lstSearchResults As IList(Of udtMSGFDBSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer,
      lstFilteredSearchResults As List(Of udtMSGFDBSearchResultType))

        Dim intIndex As Integer

        AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure already sorted by scan, charge, and SpecProb; no need to re-sort

        ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1)

        ' Now store or write out the matches that pass the filters
        ' By default, filter passing peptides have MSGFDB_SpecEValue <= 5E-7 Or EValue less than 0.75 or QValue less than 1% (but not 0)
        For intIndex = intStartIndex To intEndIndex
            If MSGFPlusResultPassesSynFilter(lstSearchResults(intIndex)) Then
                lstFilteredSearchResults.Add(lstSearchResults(intIndex))
            End If
        Next intIndex

    End Sub

    Private Sub WriteSynFHTFileHeader(
      swResultFile As TextWriter,
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
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue)
                lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue)
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
      swResultFile As TextWriter,
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
            lstData.Add(udtSearchResult.SpecEValue)
            lstData.Add(udtSearchResult.RankSpecProb.ToString)
            lstData.Add(udtSearchResult.EValue)

            If blnIncludeFDRandPepFDR Then
                lstData.Add(udtSearchResult.QValue)
                lstData.Add(udtSearchResult.PepQValue)
            ElseIf blnIncludeEFDR Then
                lstData.Add(udtSearchResult.QValue)
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

    Private Class MSGFDBSearchResultsComparerScanChargeSpecProbPeptide
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
                    ' Charge is the same; check SpecEValue
                    If x.SpecEValueNum > y.SpecEValueNum Then
                        Return 1
                    ElseIf x.SpecEValueNum < y.SpecEValueNum Then
                        Return -1
                    Else
                        ' SpecEValue is the same; check peptide
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

    Private Class MSGFDBSearchResultsComparerSpecProbScanChargePeptide
        Implements IComparer(Of udtMSGFDBSearchResultType)

        Public Function Compare(x As udtMSGFDBSearchResultType, y As udtMSGFDBSearchResultType) As Integer Implements IComparer(Of udtMSGFDBSearchResultType).Compare

            If x.SpecEValueNum > y.SpecEValueNum Then
                Return 1
            ElseIf x.SpecEValueNum < y.SpecEValueNum Then
                Return -1
            Else
                ' SpecEValueNum is the same; check scan number
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
