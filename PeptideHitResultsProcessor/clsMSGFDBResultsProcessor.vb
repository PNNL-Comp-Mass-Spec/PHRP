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

Public Class clsMSGFDBResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "September 7, 2011"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MSGFDB_FILE As String = "_msgfdb"
    Public Const FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING As String = "_PepToProtMap"

    Public Const N_TERMINUS_SYMBOL_MSGFDB As String = "_."
    Public Const C_TERMINUS_SYMBOL_MSGFDB As String = "._"

    Private Const UNKNOWN_MSGFDB_MOD_SYMBOL As Char = "?"c

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
    Private Const MSGFDB_CTERMINAL_MOD_MASS_REGEX As String =  "([0-9\.\+\-]+)$"

    Private Const MSGFDB_MOD_MASS_REGEX As String = "([+-][0-9\.]+)"

    Private Const REGEX_OPTIONS As Text.RegularExpressions.RegexOptions = Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.Singleline Or Text.RegularExpressions.RegexOptions.IgnoreCase

    ' These columns correspond to the tab-delimited file created directly by MSGF-DB
    Protected Const MSGFDBResultsFileColCount As Integer = 16
    Public Enum eMSGFDBResultsFileColumns As Integer
        SpectrumFile = 0
        SpecIndex = 1
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
        SpecProb = 12
        PValue = 13
        FDR = 14                    ' Only present if searched using -tda 1
        PepFDR = 15                 ' Only present if searched using -tda 1
    End Enum

    ' These columns correspond to the Synopsis and First-Hits files created by this class
    Protected Const MSGFDBSynFileColCount As Integer = 19
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
        SpecProb = 14
        RankSpecProb = 15                   ' Rank 1 means lowest SpecProb, 2 means next higher score, etc. (ties get the same rank)
        PValue = 16
        FDR = 17                            ' Only present if searched using -tda 1
        PepFDR = 18                         ' Only present if searched using -tda 1
    End Enum

    Protected Enum eMSGFDBModType As Integer
        Unknown = 0
        DynamicMod = 1
        StaticMod = 2
        DynNTermPeptide = 3
        DynCTermPeptide = 4
        DynNTermProtein = 5
        DynCTermProtein = 6
    End Enum

    Protected Enum eFilteredOutputFileTypeConstants As Integer
        SynFile = 0
        FHTFile = 1
    End Enum
#End Region

#Region "Structures"
    Protected Structure udtMSGFDBSearchResultType

        Public SpectrumFile As String
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
        Public SpecProb As String                   ' Smaller values are better scores (e.g. 1E-9 is better than 1E-6)
        Public SpecProbNum As Double
        Public PValue As String                     ' Smaller values are better scores (e.g. 1E-7 is better than 1E-3)
        Public PValueNum As Double
        Public FDR As String
        Public PepFDR As String
        Public RankSpecProb As Integer

        Public Sub Clear()
            SpectrumFile = String.Empty
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
        End Sub
    End Structure

    Protected Structure udtModInfoType
        Public ModName As String            ' Mod name isn't used by MSGF-DB, though it is defined in the parameter file
        Public ModMass As String            ' Storing as a string since reading from a text file and writing to a text file
        Public ModMassVal As Double
        Public Residues As String
        Public ModType As eMSGFDBModType
        Public ModSymbol As String
    End Structure

    Protected Structure udtPepToProteinMappingType
        Public Peptide As String
        Public Protein As String
        Public ResidueStart As Integer
        Public ResidueEnd As Integer
    End Structure

#End Region

#Region "Classwide Variables"
    Protected mPeptideCleavageStateCalculator As clsPeptideCleavageStateCalculator
    Protected mReplaceSymbols As System.Text.RegularExpressions.Regex
#End Region

#Region "Properties"
#End Region

    Private Sub AddCurrentRecordToSearchResults(ByRef intCurrentScanResultsCount As Integer, _
                                                ByRef udtSearchResultsCurrentScan() As udtMSGFDBSearchResultType, _
                                                ByRef udtSearchResult As udtMSGFDBSearchResultType, _
                                                ByRef strErrorLog As String)

        If intCurrentScanResultsCount >= udtSearchResultsCurrentScan.Length Then
            ReDim Preserve udtSearchResultsCurrentScan(udtSearchResultsCurrentScan.Length * 2 - 1)
        End If

        udtSearchResultsCurrentScan(intCurrentScanResultsCount) = udtSearchResult
        intCurrentScanResultsCount += 1

    End Sub

    ''' <summary>
    ''' Step through .PeptideSequenceWithMods
    ''' For each residue, check if a static mod is defined that affects that residue
    ''' For each mod symbol, determine the modification and add to objSearchResult
    ''' </summary>
    ''' <param name="objSearchResult"></param>
    ''' <param name="blnUpdateModOccurrenceCounts"></param>
    ''' <remarks></remarks>
    Private Sub AddDynamicAndStaticResidueMods(ByRef objSearchResult As clsSearchResultsMSGFDB, ByVal blnUpdateModOccurrenceCounts As Boolean)
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

                If Char.IsLetter(chChar) Then
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
                ElseIf Char.IsLetter(chMostRecentLetter) Then
                    blnSuccess = .SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                    If Not blnSuccess Then
                        If String.IsNullOrEmpty(mErrorMessage) Then
                            mErrorMessage = "SearchResultAddDynamicModification returned false for symbol " & chChar
                        End If
                        SetErrorMessage(mErrorMessage)
                    End If
                Else
                    ' We found a modification symbol but chMostRecentLetter is not a letter
                    ' Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                End If

            Next intIndex
        End With
    End Sub

    ''' <summary>
    ''' Ranks each entry (calling procedure should have already sorted the data by Scan, Charge, and SpecProb
    ''' </summary>
    ''' <param name="udtSearchResults"></param>
    ''' <param name="intStartIndex"></param>
    ''' <param name="intEndIndex"></param>
    ''' <remarks></remarks>
    Private Sub AssignRankAndDeltaNormValues(ByRef udtSearchResults() As udtMSGFDBSearchResultType, _
                                             ByVal intStartIndex As Integer, _
                                             ByVal intEndIndex As Integer)

        Dim intIndex As Integer

        Dim intLastCharge As Integer
        Dim dblLastValue As Double

        Dim intCurrentRank As Integer

        For intIndex = intStartIndex To intEndIndex
            If intIndex = intStartIndex OrElse udtSearchResults(intIndex).ChargeNum <> intLastCharge Then
                intLastCharge = udtSearchResults(intIndex).ChargeNum
                dblLastValue = udtSearchResults(intIndex).SpecProbNum
                intCurrentRank = 1
            Else
                If udtSearchResults(intIndex).SpecProbNum <> dblLastValue Then
                    dblLastValue = udtSearchResults(intIndex).SpecProbNum
                    intCurrentRank += 1
                End If
            End If

            udtSearchResults(intIndex).RankSpecProb = intCurrentRank
        Next intIndex

    End Sub

    Private Function AddModificationsAndComputeMass(ByRef objSearchResult As clsSearchResultsMSGFDB, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
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

    Protected Function ComputeCleaveageState(ByVal strSequenceWithMods As String) As Short

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
    Protected Function ComputeDelMCorrectedPPM(ByVal dblPrecursorErrorDa As Double, ByVal dblPrecursorMZ As Double, _
                                               ByVal intCharge As Integer, ByVal dblPeptideMonoisotopicMass As Double, _
                                               ByVal blnAdjustPrecursorMassForC13 As Boolean) As Double

        Dim dblPeptideDeltaMassCorrectedPpm As Double


        Dim intCorrectionCount As Integer = 0

        Dim dblPrecursorMonoMass As Double

        ' Compute the original value for the precursor monoisotopic mass
        dblPrecursorMonoMass = dblPrecursorMZ * intCharge - intCharge * clsPeptideMassCalculator.MASS_PROTON

        dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPublic(dblPrecursorErrorDa, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, dblPeptideMonoisotopicMass)

        Return dblPeptideDeltaMassCorrectedPpm

    End Function


    Protected Function ComputePeptideMass(ByVal strPeptide As String, ByVal dblTotalModMass As Double) As Double

        Dim strCleanSequence As String
        Dim dblMass As Double

        strCleanSequence = GetCleanSequence(strPeptide)

        Dim objMassCalculator As clsPeptideMassCalculator = New clsPeptideMassCalculator()

        dblMass = objMassCalculator.ComputeSequenceMass(strCleanSequence)
        dblmass += dblTotalModMass

        Return dblMass

    End Function

    ''' <summary>
    ''' Parses the digits in strModValues to convert them to one or more modification symbols
    ''' </summary>
    ''' <param name="strModDigits">Example: +57.021 or +79.9663+14.0157 or -18.0106</param>
    ''' <param name="strModSymbols"></param>
    ''' <returns>True if success; false if a problem</returns>
    ''' <remarks></remarks>
    Protected Function ConvertMGSFModMassesToSymbols(ByVal strModDigits As String, _
                                                     ByRef strModSymbols As String, _
                                                     ByRef udtMSGFDBModInfo() As udtModInfoType, _
                                                     ByVal blnNterminalMod As Boolean, _
                                                     ByVal blnPossibleCTerminalMod As Boolean, _
                                                     ByRef dblModMassFound As Double) As Boolean

        Static reModMassRegEx As New System.Text.RegularExpressions.Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS)

        Dim reMatches As System.Text.RegularExpressions.MatchCollection
        Dim blnSuccess As Boolean = False

        Dim blnMatchFound As Boolean
        Dim blnTestMod As Boolean

        Dim intBestMatchIndex As Integer
        Dim dblBestMassDiff As Double
        Dim dblCandidateMassDiff As Double
        Dim intModSymbolsFound As Integer = 0

        strModSymbols = String.Empty
        dblModMassFound = 0

        reMatches = reModMassRegEx.Matches(strModDigits)

        For Each reMatch As System.Text.RegularExpressions.Match In reMatches
            Dim strModMass As String
            strModMass = reMatch.Value

            ' Convert strModMass to a mass value
            Dim dblModMass As Double
            dblModMass = Double.Parse(strModMass)
            dblModMassFound += dblModMass

            blnMatchFound = False
            Do
                intBestMatchIndex = -1
                dblCandidateMassDiff = Double.MaxValue

                ' Step through the known modifications to find the closest match
                For intIndex As Integer = 0 To udtMSGFDBModInfo.Length - 1

                    If blnNterminalMod Then
                        ' Only test N-terminal mods
                        If udtMSGFDBModInfo(intIndex).ModType = eMSGFDBModType.DynNTermPeptide OrElse _
                            udtMSGFDBModInfo(intIndex).ModType = eMSGFDBModType.DynNTermProtein Then
                            blnTestMod = True
                        Else
                            blnTestMod = False
                        End If

                    ElseIf blnPossibleCTerminalMod Then
                        ' Only test C-terminal mods
                        If udtMSGFDBModInfo(intIndex).ModType = eMSGFDBModType.DynCTermPeptide OrElse _
                            udtMSGFDBModInfo(intIndex).ModType = eMSGFDBModType.DynCTermProtein Then
                            blnTestMod = True
                        Else
                            blnTestMod = False
                        End If

                    Else
                        blnTestMod = True
                    End If

                    If blnTestMod Then
                        dblCandidateMassDiff = Math.Abs(udtMSGFDBModInfo(intIndex).ModMassVal - dblModMass)
                        If dblCandidateMassDiff < 0.25 Then
                            ' Possible match found
                            If intBestMatchIndex < 0 OrElse dblCandidateMassDiff < dblBestMassDiff Then
                                intBestMatchIndex = intIndex
                                dblBestMassDiff = dblCandidateMassDiff
                                blnMatchFound = True
                            End If
                        End If

                    End If

                Next

                If Not blnMatchFound AndAlso (blnNterminalMod Or blnPossibleCTerminalMod) Then
                    ' Set these to false, then search again
                    blnNterminalMod = False
                    blnPossibleCTerminalMod = False
                Else
                    Exit Do
                End If

            Loop While Not blnMatchFound

            If blnMatchFound Then
                ' Match found; use the mod symbol
                strModSymbols &= udtMSGFDBModInfo(intBestMatchIndex).ModSymbol
                intModSymbolsFound += 1
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
    ''' <param name="udtMSGFDBModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
    ''' <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
    ''' <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CreateFHTorSYNResultsFile(ByVal strInputFilePath As String, _
                                                 ByVal strOutputFilePath As String, _
                                                 ByRef udtMSGFDBModInfo() As udtModInfoType, _
                                                 ByVal eFilteredOutputFileType As eFilteredOutputFileTypeConstants, _
                                                 Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean

        Dim srDataFile As System.IO.StreamReader
        Dim swResultFile As System.IO.StreamWriter

        Dim intPreviousScan As Integer

        Dim strLineIn As String
        Dim protein As String = String.Empty

        Dim udtSearchResult As udtMSGFDBSearchResultType
        udtSearchResult.Clear()

        Dim intSearchResultsCount As Integer
        Dim udtSearchResults() As udtMSGFDBSearchResultType

        Dim intFilteredSearchResultCount As Integer
        Dim udtFilteredSearchResults() As udtMSGFDBSearchResultType

        Dim blnHeaderParsed As Boolean
        Dim blnIncludeFDRandPepFDR As Boolean = False

        Dim intColumnMapping() As Integer = Nothing

        Dim intResultsProcessed As Integer
        Dim intResultID As Integer = 0

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean

        Dim strErrorLog As String = String.Empty

        Try

            Try
                ' Open the input file and parse it
                ' Initialize the stream reader and the stream Text writer
                srDataFile = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

                swResultFile = New System.IO.StreamWriter(New System.IO.FileStream(strOutputFilePath, IO.FileMode.Create, IO.FileAccess.Write, IO.FileShare.Read))

                strErrorLog = String.Empty
                intResultsProcessed = 0
                blnHeaderParsed = False

                ' Initialize array that will hold all of the records in the MSGFDB result file
                intSearchResultsCount = 0
                ReDim udtSearchResults(999)

                ' Initialize the array that will hold all of the records that will ultimately be written out to disk
                intFilteredSearchResultCount = 0
                ReDim udtFilteredSearchResults(999)

                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing
                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then

                        If Not blnHeaderParsed Then
                            blnSuccess = ParseMSGFDBResultsFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True

                            If intColumnMapping(eMSGFDBResultsFileColumns.FDR) >= 0 Or intColumnMapping(eMSGFDBResultsFileColumns.PepFDR) >= 0 Then
                                blnIncludeFDRandPepFDR = True
                            Else
                                blnIncludeFDRandPepFDR = False
                            End If

                            ' Write the header line
                            WriteSynFHTFileHeader(swResultFile, strErrorLog, blnIncludeFDRandPepFDR)
                        Else

                            blnValidSearchResult = ParseMSGFDBResultsFileEntry(strLineIn, udtMSGFDBModInfo, udtSearchResult, _
                                                                               strErrorLog, intResultsProcessed, _
                                                                               intColumnMapping)

                            If blnValidSearchResult Then
                                AddCurrentRecordToSearchResults(intSearchResultsCount, udtSearchResults, udtSearchResult, strErrorLog)
                            End If

                            ' Update the progress
                            UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))
                            intResultsProcessed += 1
                        End If
                    End If
                Loop

                ' Sort the SearchResults by scan, charge, and ascending SpecProb
                Array.Sort(udtSearchResults, 0, intSearchResultsCount, New MSGFDBSearchResultsComparerScanChargeSpecProbPeptide)

                ' Now filter the data

                ' Initialize variables
                Dim intStartIndex As Integer = 0
                Dim intEndIndex As Integer

                If intSearchResultsCount > 0 Then
                    intPreviousScan = udtSearchResults(0).ScanNum
                Else
                    intPreviousScan = Int32.MinValue
                End If

                For intIndex As Integer = 0 To intSearchResultsCount - 1

                    If (intPreviousScan <> udtSearchResults(intIndex).ScanNum) OrElse _
                        intIndex = intSearchResultsCount - 1 Then

                        If intIndex = intSearchResultsCount - 1 Then
                            intEndIndex = intIndex
                        Else
                            intEndIndex = intIndex - 1
                        End If

                        ' New scan encountered or last result; call either StoreTopFHTMatch or StoreSynMatches to store the results
                        If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
                            StoreSynMatches(swResultFile, intResultID, udtSearchResults, intStartIndex, intEndIndex, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                        Else
                            StoreTopFHTMatch(swResultFile, intResultID, udtSearchResults, intStartIndex, intEndIndex, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                        End If

                        intPreviousScan = udtSearchResults(intIndex).ScanNum
                        intStartIndex = intIndex
                    End If
                Next

                ' Sort the data in udtFilteredSearchResults then write out to disk
                SortAndWriteFilteredSearchResults(swResultFile, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog, blnIncludeFDRandPepFDR)

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
                If Not srDataFile Is Nothing Then
                    srDataFile.Close()
                    srDataFile = Nothing
                End If
                If Not swResultFile Is Nothing Then
                    swResultFile.Close()
                    swResultFile = Nothing
                End If
            End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Function ComputeMass(ByVal strEmpiricalformula As String) As Double
        ' CompositionStr (C[Num]H[Num]N[Num]O[Num]S[Num]P[Num])
        ' 	- C (Carbon), H (Hydrogen), N (Nitrogen), O (Oxygen), S (Sulfer) and P (Phosphorus) are allowed.
        ' 	- Atom can be omitted. The sequence of atoms must be followed. 
        ' 	- Negative numbers are allowed.
        ' 	- E.g. C2H2O1 (valid), H2C1O1 (invalid) 

        Static reAtomicFormulaRegEx As New System.Text.RegularExpressions.Regex("[CHNOSP][+-]?\d*", REGEX_OPTIONS)

        Dim reMatches As Text.RegularExpressions.MatchCollection
        Dim strElement As String
        Dim intCount As Integer
        Dim dblMass As Double = 0

        reMatches = reAtomicFormulaRegEx.Matches(strEmpiricalformula)

        If reMatches.Count > 0 Then

            For Each reMatch As System.Text.RegularExpressions.Match In reMatches

                strElement = reMatch.Value.Chars(0)
                If reMatch.Value.Length > 1 Then
                    If Not Integer.TryParse(reMatch.Value.Substring(1), intCount) Then
                        SetErrorMessage("Error parsing empirical formula '" & strEmpiricalformula & "', number not found in " & reMatch.Value)
                        SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                        dblMass = 0
                        Exit For
                    End If
                Else
                    intCount = 1
                End If

                Select Case strElement.ToUpper()
                    Case "C" : dblMass += intCount * 12
                    Case "H" : dblMass += intCount * clsPeptideMassCalculator.MASS_HYDROGEN
                    Case "N" : dblMass += intCount * 14.003074
                    Case "O" : dblMass += intCount * 15.994915
                    Case "S" : dblMass += intCount * 31.972072
                    Case "P" : dblMass += intCount * 30.973763
                    Case Else
                        ' Unknown element
                        SetErrorMessage("Error parsing empirical formula '" & strEmpiricalformula & "', unknown element " & strElement)
                        SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                        dblMass = 0
                        Exit For
                End Select
            Next
        End If

        Return dblMass

    End Function

    Private Function CIntSafe(ByVal strValue As String, ByVal intDefaultValue As Integer) As Integer
        Try
            ' Note: Integer.Parse() fails if strValue contains a decimal point, even if it is "8.000"
            ' Thus, we're using CInt() instead
            Return CInt(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return intDefaultValue
        End Try
    End Function

    Private Function CDblSafe(ByVal strValue As String, ByVal dblgDefaultValue As Double) As Double
        Try
            Return Double.Parse(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return dblgDefaultValue
        End Try
    End Function

    ''' <summary>
    ''' Extracts mod info from either a MSGF-DB param file or from a MSGFDB_Mods.txt file
    ''' </summary>
    ''' <param name="strMSGFDBParamFilePath"></param>
    ''' <param name="udtModList"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function ExtractModInfoFromMSGFDBModsFile(ByVal strMSGFDBParamFilePath As String, ByRef udtModList() As udtModInfoType) As Boolean

        Const STATIC_MOD_TAG As String = "StaticMod"
        Const DYNAMIC_MOD_TAG As String = "DynamicMod"

        Dim srInFile As System.IO.StreamReader = Nothing

        Dim strLineIn As String
        Dim strSplitLine As String()
        Dim strModSpec As String
        Dim intPoundIndex As Integer

        Dim intModCount As Integer
        Dim intUnnamedModID As Integer

        Dim blnSuccess As Boolean = False

        Try
            ' Initialize udtModList and intUnnamedModID
            intModCount = 0
            ReDim udtModList(-1)

            intUnnamedModID = 0

            If String.IsNullOrEmpty(strMSGFDBParamFilePath) Then
                SetErrorMessage("MSGFDB  Parameter File name not defined; unable to extract mod info")
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                Return False
            End If

            ' Read the contents of the parameter (or mods) file
            srInFile = New System.IO.StreamReader(New System.IO.FileStream(strMSGFDBParamFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

            Do While srInFile.Peek <> -1
                strLineIn = srInFile.ReadLine.Trim()

                If strLineIn.Length > 0 Then

                    strModSpec = String.Empty

                    If strLineIn.StartsWith("#"c) Then
                        ' Comment line; skip it
                    Else
                        intPoundIndex = strLineIn.IndexOf("#"c)
                        If intPoundIndex > 0 Then
                            ' Trim out the comment from the end
                            strLineIn = strLineIn.Substring(0, intPoundIndex).Trim()
                        End If

                        strModSpec = ValidateIsValidModSpec(strLineIn, STATIC_MOD_TAG)
                        If Not String.IsNullOrEmpty(strModSpec) Then
                            ' Line resembles StaticMod=C2H3N1O1,C,fix,any,Carbamidomethylation
                        Else
                            strModSpec = ValidateIsValidModSpec(strLineIn, DYNAMIC_MOD_TAG)
                            If Not String.IsNullOrEmpty(strModSpec) Then
                                ' Line resembles DynamicMod=C2H3NO, *,  opt, N-term,   Carbamidomethylation
                            Else
                                If strLineIn.Contains(",opt,") OrElse strLineIn.Contains(",fix,") Then
                                    strModSpec = strLineIn
                                End If
                            End If
                        End If

                    End If

                    If Not String.IsNullOrEmpty(strModSpec) Then

                        ' Modification definition line found

                        ' Split the line on commas
                        strSplitLine = strModSpec.Split(","c)

                        If strSplitLine.Length >= 5 Then

                            If udtModList.Length = 0 Then
                                ReDim udtModList(0)
                            ElseIf intModCount >= udtModList.Length Then
                                ReDim Preserve udtModList(udtModList.Length * 2 - 1)
                            End If

                            With udtModList(intModCount)
                                .ModMass = strSplitLine(0).Trim()

                                If Not Double.TryParse(.ModMass, .ModMassVal) Then
                                    ' Mod is specified as an empirical formula
                                    ' Compute the mass
                                    .ModMassVal = ComputeMass(.ModMass)
                                End If

                                .Residues = strSplitLine(1).Trim()
                                .ModSymbol = UNKNOWN_MSGFDB_MOD_SYMBOL

                                Select Case strSplitLine(2).Trim().ToLower()
                                    Case "opt"
                                        .ModType = eMSGFDBModType.DynamicMod
                                    Case "fix"
                                        .ModType = eMSGFDBModType.StaticMod
                                    Case Else
                                        SetErrorMessage("Warning: Unrecognized Mod Type in the MSGFDB parameter file; should be 'opt' or 'fix'")
                                        .ModType = eMSGFDBModType.DynamicMod
                                End Select

                                Select Case strSplitLine(3).Trim().ToLower().Replace("-", String.Empty)
                                    Case "any"
                                        ' Leave .ModType unchanged; this is a static or dynamic mod (fix or opt)
                                    Case "nterm"
                                        .ModType = eMSGFDBModType.DynNTermPeptide
                                    Case "cterm"
                                        .ModType = eMSGFDBModType.DynCTermPeptide
                                    Case "protnterm"
                                        ' Includes Prot-N-Term, Prot-n-Term, ProtNTerm, etc.
                                        .ModType = eMSGFDBModType.DynNTermProtein
                                    Case "protcterm"
                                        ' Includes Prot-C-Term, Prot-c-Term, ProtCterm, etc.
                                        .ModType = eMSGFDBModType.DynCTermProtein
                                    Case Else
                                        SetErrorMessage("Warning: Unrecognized Mod Type in the MSGFDB parameter file; should be 'any', 'N-term', 'C-term', 'Prot-N-term', or 'Prot-C-term'")
                                End Select

                                .ModName = strSplitLine(4).Trim()
                                If String.IsNullOrEmpty(.ModName) Then
                                    intUnnamedModID += 1
                                    .ModName = "UnnamedMod" & intUnnamedModID.ToString
                                End If

                            End With

                            intModCount += 1
                        End If
                    End If
                End If
            Loop

            ' Shrink udtModList to the appropriate length
            ReDim Preserve udtModList(intModCount - 1)

            Console.WriteLine()

            blnSuccess = True

        Catch ex As Exception
            SetErrorMessage("Error reading the MSGF-DB parameter file (" & System.IO.Path.GetFileName(strMSGFDBParamFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            blnSuccess = False
        Finally
            If Not srInFile Is Nothing Then
                srInFile.Close()
            End If
        End Try

        Return blnSuccess

    End Function

    Protected Function GetCleanSequence(ByVal strSequenceWithMods As String) As String
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Return GetCleanSequence(strSequenceWithMods, strPrefix, strSuffix)
    End Function

    Protected Function GetCleanSequence(ByVal strSequenceWithMods As String, ByVal strPrefix As String, ByVal strSuffix As String) As String

        Dim strPrimarySequence As String = String.Empty
        strPrefix = String.Empty
        strSuffix = String.Empty

        If clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then

            ' Remove any non-letter characters when calling .ComputeCleavageState()

            strPrimarySequence = mReplaceSymbols.Replace(strPrimarySequence, String.Empty)

        Else
            ' Unable to determine cleavage-state
            strPrimarySequence = mReplaceSymbols.Replace(strSequenceWithMods, String.Empty)
        End If

        Return strPrimarySequence

    End Function

    Private Sub InitializeLocalVariables()

        ' Initialize mPeptideCleavageStateCalculator
        If mPeptideCleavageStateCalculator Is Nothing Then
            mPeptideCleavageStateCalculator = New clsPeptideCleavageStateCalculator
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin)
        End If

        ' Define a RegEx to replace all of the non-letter characters
        mReplaceSymbols = New System.Text.RegularExpressions.Regex("[^A-Za-z]", System.Text.RegularExpressions.RegexOptions.Compiled)

    End Sub

    ''' <summary>
    ''' Load the PeptideToProteinMap information; in addition, writes out an updated _inspect_PepToProtMap.txt file with the new mod symbols and corrected terminii symbols
    ''' </summary>
    ''' <param name="strPepToProteinMapFilePath"></param>
    ''' <param name="strOutputFolderPath"></param>
    ''' <param name="udtMSGFDBModInfo"></param>
    ''' <param name="udtPepToProteinMapping"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function LoadPeptideToProteinMapInfo(ByVal strPepToProteinMapFilePath As String, _
                                                   ByVal strOutputFolderPath As String, _
                                                   ByRef udtMSGFDBModInfo() As udtModInfoType, _
                                                   ByRef udtPepToProteinMapping() As udtPepToProteinMappingType) As Boolean


        Dim srInFile As System.IO.StreamReader
        Dim swOutFile As System.IO.StreamWriter

        Dim strMTSPepToProteinMapFilePath As String = String.Empty
        Dim strMTSCompatiblePeptide As String

        Dim strLineIn As String
        Dim strSplitLine As String()

        Dim intLinesRead As Integer
        Dim intValue As Integer

        Dim intProteinCount As Integer
        Dim blnSuccess As Boolean = False

        Dim dblTotalModMass As Double = 0

        Try
            ' Initialize udtPepToProteinMapping
            If strPepToProteinMapFilePath Is Nothing OrElse strPepToProteinMapFilePath.Length = 0 Then
                SetErrorMessage("Warning: PepToProteinMap file is not defined")
                Return False
            ElseIf Not System.IO.File.Exists(strPepToProteinMapFilePath) Then
                SetErrorMessage("Warning: PepToProteinMap file does not exist")
                Return False
            End If

            intProteinCount = 0
            ReDim udtPepToProteinMapping(999)

            ' Open strProteinToPeptideMappingFilePath for reading
            srInFile = New System.IO.StreamReader(New System.IO.FileStream(strPepToProteinMapFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

            strMTSPepToProteinMapFilePath = System.IO.Path.Combine(strOutputFolderPath, System.IO.Path.GetFileNameWithoutExtension(strPepToProteinMapFilePath) & "MTS.txt")
            swOutFile = New System.IO.StreamWriter(New System.IO.FileStream(strMTSPepToProteinMapFilePath, IO.FileMode.Create, IO.FileAccess.Write, IO.FileShare.Read))

            intLinesRead = 0
            Do While srInFile.Peek <> -1
                strLineIn = srInFile.ReadLine

                strLineIn = strLineIn.Trim

                If strLineIn.Length > 0 Then

                    ' Split the line on tabs
                    strSplitLine = strLineIn.Split(ControlChars.Tab)

                    If strSplitLine.Length >= 4 Then

                        If intLinesRead = 0 AndAlso Not Integer.TryParse(strSplitLine(2), intValue) Then
                            ' Header line
                            swOutFile.WriteLine(strLineIn)
                        Else
                            If intProteinCount >= udtPepToProteinMapping.Length Then
                                ReDim Preserve udtPepToProteinMapping(udtPepToProteinMapping.Length * 2 - 1)
                            End If

                            ' Replace any mod text names in the peptide sequence with the appropriate mod symbols
                            ' In addition, replace the * terminus symbols with dashes
                            strMTSCompatiblePeptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(strSplitLine(0)), udtMSGFDBModInfo, dblTotalModMass)

                            With udtPepToProteinMapping(intProteinCount)
                                .Peptide = strMTSCompatiblePeptide
                                .Protein = String.Copy(strSplitLine(1))
                                Integer.TryParse(strSplitLine(2), .ResidueStart)
                                Integer.TryParse(strSplitLine(3), .ResidueEnd)
                            End With


                            swOutFile.WriteLine(strMTSCompatiblePeptide & ControlChars.Tab & _
                                                strSplitLine(1) & ControlChars.Tab & _
                                                strSplitLine(2) & ControlChars.Tab & _
                                                strSplitLine(3))

                            intProteinCount += 1
                        End If
                    End If
                End If

            Loop

            ' Shrink udtPepToProteinMapping to the appropriate length
            ReDim Preserve udtPepToProteinMapping(intProteinCount - 1)

            Console.WriteLine()

            blnSuccess = True

        Catch ex As Exception
            SetErrorMessage("Error reading the Peptide to Protein Map File (" & System.IO.Path.GetFileName(strPepToProteinMapFilePath) & ") " & _
                            "and writing the new map file (" & System.IO.Path.GetFileName(strMTSPepToProteinMapFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
            blnSuccess = False
        Finally
            If Not srInFile Is Nothing Then
                srInFile.Close()
            End If
            If Not swOutFile Is Nothing Then
                swOutFile.Close()
            End If
        End Try

        Return blnSuccess

    End Function

    Protected Function ParseMSGFDBSynopsisFile(ByVal strInputFilePath As String, ByRef udtPepToProteinMapping() As udtPepToProteinMappingType, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim srDataFile As System.IO.StreamReader

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

        Dim blnHeaderParsed As Boolean
        Dim intColumnMapping() As Integer = Nothing

        Dim strCurrentPeptideWithMods As String = String.Empty
        Dim strCurrentProtein As String

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean
        Dim blnFirstMatchForGroup As Boolean

        Dim strErrorLog As String = String.Empty

        Dim objPeptideSearchComparer As PepToProteinMappingPeptideSearchComparer

        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            ' Initialize objSearchResult
            objSearchResult = New clsSearchResultsMSGFDB(mPeptideMods)

            ' Initialize htPeptidesFoundForSpecProbLevel
            htPeptidesFoundForSpecProbLevel = New Hashtable
            strPreviousSpecProb = String.Empty

            ' Initialize objPeptideSearchComparer
            objPeptideSearchComparer = New PepToProteinMappingPeptideSearchComparer

            ' Assure that udtPepToProteinMapping is sorted on peptide
            If udtPepToProteinMapping.Length > 1 Then
                Array.Sort(udtPepToProteinMapping, New PepToProteinMappingComparer)
            End If

            Try
                UpdateSearchResultEnzymeAndTerminusInfo(objSearchResult)

                ' Open the input file and parse it
                ' Initialize the stream reader
                srDataFile = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

                strErrorLog = String.Empty
                intResultsProcessed = 0
                blnHeaderParsed = False

                ' Create the output files
                blnSuccess = MyBase.InitializeSequenceOutputFiles(strInputFilePath)

                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing

                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then

                        If Not blnHeaderParsed Then
                            blnSuccess = ParseMSGFDBSynFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True
                        Else

                            blnValidSearchResult = ParseMSGFDBSynFileEntry(strLineIn, objSearchResult, strErrorLog, _
                                                                           intResultsProcessed, intColumnMapping, _
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
                                    ' Reset htPeptidesFoundForScan
                                    htPeptidesFoundForSpecProbLevel.Clear()

                                    ' Update strPreviousSpecProb
                                    strPreviousSpecProb = objSearchResult.SpecProb

                                    ' Append a new entry to htPeptidesFoundForScan
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


                                If udtPepToProteinMapping.Length > 0 Then
                                    ' Add the additional proteins for this peptide

                                    ' Use binary search to find this peptide in udtPepToProteinMapping
                                    intPepToProteinMapIndex = Array.BinarySearch(udtPepToProteinMapping, strCurrentPeptideWithMods, objPeptideSearchComparer)

                                    If intPepToProteinMapIndex >= 0 Then
                                        ' Step Backward until the first match is found
                                        Do While intPepToProteinMapIndex > 0 AndAlso udtPepToProteinMapping(intPepToProteinMapIndex - 1).Peptide = strCurrentPeptideWithMods
                                            intPepToProteinMapIndex -= 1
                                        Loop

                                        ' Call MyBase.SaveResultsFileEntrySeqInfo for each entry in udtPepToProteinMapping() for peptide , skipping objSearchResult.ProteinName
                                        strCurrentProtein = String.Copy(objSearchResult.ProteinName)
                                        Do
                                            If udtPepToProteinMapping(intPepToProteinMapIndex).Protein <> strCurrentProtein Then
                                                objSearchResult.ProteinName = String.Copy(udtPepToProteinMapping(intPepToProteinMapIndex).Protein)
                                                MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResult, clsSearchResultsBaseClass), False)
                                            End If

                                            intPepToProteinMapIndex += 1
                                        Loop While intPepToProteinMapIndex < udtPepToProteinMapping.Length AndAlso strCurrentPeptideWithMods = udtPepToProteinMapping(intPepToProteinMapIndex).Peptide
                                    Else
                                        ' Match not found; this is unexpected
                                        Console.WriteLine("Warning: no match for '" & strCurrentPeptideWithMods & "' in udtPepToProteinMapping")
                                    End If
                                End If

                            End If

                            ' ToDo: Write out this line to an updated _Syn.txt file since we have now computed .PeptideDeltaMassCorrectedPpm


                            ' Update the progress
                            UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))

                            intResultsProcessed += 1
                        End If
                    End If

                Loop

                If mCreateModificationSummaryFile Then
                    ' Create the modification summary file
                    strModificationSummaryFilePath = MyBase.ReplaceFilenameSuffix(strInputFilePath, "", FILENAME_SUFFIX_MOD_SUMMARY)
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
                If Not srDataFile Is Nothing Then
                    srDataFile.Close()
                    srDataFile = Nothing
                End If
                MyBase.CloseSequenceOutputFiles()
            End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess


    End Function

    Private Function ParseMSGFDBResultsFileEntry(ByRef strLineIn As String, _
                                                 ByRef udtMSGFDBModInfo() As udtModInfoType, _
                                                 ByRef udtSearchResult As udtMSGFDBSearchResultType, _
                                                 ByRef strErrorLog As String, _
                                                 ByVal intResultsProcessed As Integer, _
                                                 ByRef intColumnMapping() As Integer) As Boolean

        ' Parses an entry from the MSGF-DB results file

        Dim strSplitLine() As String = Nothing

        Dim blnValidSearchResult As Boolean

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset udtSearchResult
            udtSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length >= 13 Then

                With udtSearchResult
                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpectrumFile), .SpectrumFile) Then
                        Throw New EvaluateException("SpectrumFile column is missing or invalid")
                    End If
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpecIndex), .SpecIndex)

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Scan), .Scan) Then
                        Throw New EvaluateException("Scan column is missing or invalid")
                    End If
                    .ScanNum = CIntSafe(.Scan, 0)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.FragMethod), .FragMethod)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PrecursorMZ), .PrecursorMZ)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Charge), .Charge)
                    .ChargeNum = CShort(CIntSafe(.Charge, 0))

                    ' PMError could be in PPM or Da; header line will have PMError(ppm) or PMError(Da)
                    Dim dblPrecursorErrorDa As Double

                    If intColumnMapping(eMSGFDBResultsFileColumns.PMErrorPPM) >= 0 Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PMErrorPPM), .PMErrorPPM)
                    Else
                        GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PMErrorDa), .PMErrorDa)
                        dblPrecursorErrorDa = CDblSafe(.PMErrorDa, 0)
                        .PMErrorPPM = String.Empty              ' We'll populate this column later in this function
                    End If


                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Peptide), .Peptide) Then
                        Throw New EvaluateException("Peptide column is missing or invalid")
                    End If

                    ' Replace any mod text values in the peptide sequence with the appropriate mod symbols
                    ' In addition, replace the _ terminus symbols with dashes
                    Dim dblTotalModMass As Double
                    .Peptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(.Peptide), udtMSGFDBModInfo, dblTotalModMass)

                    ' Compute monoisotopic mass of the peptide
                    Dim dblPeptideMonoisotopicMass As Double
                    dblPeptideMonoisotopicMass = ComputePeptideMass(.Peptide, dblTotalModMass)


                    ' Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    .MH = NumToString(clsPeptideMassCalculator.ConvoluteMass(dblPeptideMonoisotopicMass, 0, 1), 5, True)

                    If Not String.IsNullOrEmpty(.PMErrorPPM) Then

                        ' Convert the ppm-based PM Error to Da-based

                        Dim dblPMErrorPPM As Double
                        Dim dblPrecursorMZ As Double

                        If Double.TryParse(.PrecursorMZ, dblPrecursorMZ) Then
                            ' Note that since .PMErrorPPM is present, the Precursor m/z is a C13-corrected m/z value
                            ' In other words, it may not be the actual m/z selected for fragmentation.

                            If Double.TryParse(.PMErrorPPM, dblPMErrorPPM) Then
                                dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPMErrorPPM, dblPeptideMonoisotopicMass)

                                ' Note that this will be a C13-corrected precursor error; not the true precursor error
                                .PMErrorDa = NumToString(dblPrecursorErrorDa, 5, True)
                            End If
                        End If

                    Else

                        Dim dblPrecursorMZ As Double
                        If Double.TryParse(.PrecursorMZ, dblPrecursorMZ) Then
                            Dim dblPeptideDeltaMassCorrectedPpm As Double

                            dblPeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMZ, _
                                                                                      .ChargeNum, dblPeptideMonoisotopicMass, True)

                            .PMErrorPPM = NumToString(dblPeptideDeltaMassCorrectedPpm, 4, True)

                        End If

                    End If


                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.Protein), .Protein)
                    .Protein = TruncateProteinName(.Protein)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.DeNovoScore), .DeNovoScore)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.MSGFScore), .MSGFScore)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.SpecProb), .SpecProb)
                    If Not Double.TryParse(.SpecProb, .SpecProbNum) Then .SpecProbNum = 0

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PValue), .PValue)
                    If Not Double.TryParse(.PValue, .PValueNum) Then .PValueNum = 0

                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.FDR), .FDR)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSGFDBResultsFileColumns.PepFDR), .PepFDR)

                    .NTT = ComputeCleaveageState(.Peptide).ToString()
                End With

                blnValidSearchResult = True
            End If

        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MSGF-DB Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MSGF-DB Results in ParseMSGFDBResultsFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    Private Function ParseMSGFDBResultsFileHeaderLine(ByVal strLineIn As String, _
                                                      ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        ' The expected header line is:
        ' #SpecFile	SpecIndex	Scan#	FragMethod	Precursor	PMError(Da)	    Charge	Peptide	Protein	DeNovoScore	MSGFScore	SpecProb	P-value
        ' or
        ' #SpecFile	SpecIndex	Scan#	FragMethod	Precursor	PMError(ppm)	Charge	Peptide	Protein	DeNovoScore	MSGFScore	SpecProb	P-value

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMSGFDBResultsFileColumns
        Dim lstColumnNames As System.Collections.Generic.SortedDictionary(Of String, eMSGFDBResultsFileColumns)
        lstColumnNames = New System.Collections.Generic.SortedDictionary(Of String, eMSGFDBResultsFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(MSGFDBResultsFileColCount - 1)

        lstColumnNames.Add("#SpecFile", eMSGFDBResultsFileColumns.SpectrumFile)
        lstColumnNames.Add("SpecIndex", eMSGFDBResultsFileColumns.SpecIndex)
        lstColumnNames.Add("Scan#", eMSGFDBResultsFileColumns.Scan)
        lstColumnNames.Add("FragMethod", eMSGFDBResultsFileColumns.FragMethod)
        lstColumnNames.Add("Precursor", eMSGFDBResultsFileColumns.PrecursorMZ)
        lstColumnNames.Add("PMError(Da)", eMSGFDBResultsFileColumns.PMErrorDa)
        lstColumnNames.Add("PMError(ppm)", eMSGFDBResultsFileColumns.PMErrorPPM)
        lstColumnNames.Add("Charge", eMSGFDBResultsFileColumns.Charge)
        lstColumnNames.Add("Peptide", eMSGFDBResultsFileColumns.Peptide)
        lstColumnNames.Add("Protein", eMSGFDBResultsFileColumns.Protein)
        lstColumnNames.Add("DeNovoScore", eMSGFDBResultsFileColumns.DeNovoScore)
        lstColumnNames.Add("MSGFScore", eMSGFDBResultsFileColumns.MSGFScore)
        lstColumnNames.Add("SpecProb", eMSGFDBResultsFileColumns.SpecProb)
        lstColumnNames.Add("P-value", eMSGFDBResultsFileColumns.PValue)
        lstColumnNames.Add("FDR", eMSGFDBResultsFileColumns.FDR)
        lstColumnNames.Add("PepFDR", eMSGFDBResultsFileColumns.PepFDR)

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
            SetErrorMessage("Error parsing header in MSGFDB results file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMSGFDBSynFileHeaderLine(ByVal strLineIn As String, _
                                                  ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMSFDBSynFileColumns
        Dim lstColumnNames As System.Collections.Generic.SortedDictionary(Of String, eMSFDBSynFileColumns)
        lstColumnNames = New System.Collections.Generic.SortedDictionary(Of String, eMSFDBSynFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(MSGFDBSynFileColCount - 1)

        lstColumnNames.Add("ResultID", eMSFDBSynFileColumns.ResultID)
        lstColumnNames.Add("Scan", eMSFDBSynFileColumns.Scan)
        lstColumnNames.Add("FragMethod", eMSFDBSynFileColumns.FragMethod)
        lstColumnNames.Add("SpecIndex", eMSFDBSynFileColumns.SpecIndex)
        lstColumnNames.Add("Charge", eMSFDBSynFileColumns.Charge)
        lstColumnNames.Add("PrecursorMZ", eMSFDBSynFileColumns.PrecursorMZ)
        lstColumnNames.Add("DelM", eMSFDBSynFileColumns.DelM)
        lstColumnNames.Add("DelM_PPM", eMSFDBSynFileColumns.DelMPPM)
        lstColumnNames.Add("MH", eMSFDBSynFileColumns.MH)
        lstColumnNames.Add("Peptide", eMSFDBSynFileColumns.Peptide)
        lstColumnNames.Add("Protein", eMSFDBSynFileColumns.Protein)
        lstColumnNames.Add("NTT", eMSFDBSynFileColumns.NTT)
        lstColumnNames.Add("DeNovoScore", eMSFDBSynFileColumns.DeNovoScore)
        lstColumnNames.Add("MSGFScore", eMSFDBSynFileColumns.MSGFScore)
        lstColumnNames.Add("MSGFDB_SpecProb", eMSFDBSynFileColumns.SpecProb)
        lstColumnNames.Add("Rank_MSGFDB_SpecProb", eMSFDBSynFileColumns.RankSpecProb)
        lstColumnNames.Add("PValue", eMSFDBSynFileColumns.PValue)
        lstColumnNames.Add("FDR", eMSFDBSynFileColumns.FDR)
        lstColumnNames.Add("PepFDR", eMSFDBSynFileColumns.PepFDR)

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
            SetErrorMessage("Error parsing header in MSGFDB synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMSGFDBSynFileEntry(ByRef strLineIn As String, _
                                             ByRef objSearchResult As clsSearchResultsMSGFDB, _
                                             ByRef strErrorLog As String, _
                                             ByVal intResultsProcessed As Integer, _
                                             ByRef intColumnMapping() As Integer, _
                                             ByRef strPeptideSequenceWithMods As String) As Boolean

        ' Parses an entry from the MSGFDB Synopsis file

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
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.SpecProb), .SpecProb)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.RankSpecProb), .RankSpecProb)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PValue), .PValue)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.FDR), .FDR)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSFDBSynFileColumns.PepFDR), .PepFDR)

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

    ' Main processing function
    Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean
        ' Returns True if success, False if failure

        Dim ioInputFile As System.IO.FileInfo

        Dim strInputFilePathFull As String
        Dim strOutputFilePath As String
        Dim strSynOutputFilePath As String
        Dim strPepToProteinMapFilePath As String

        Dim udtMSGFDBModInfo() As udtModInfoType
        Dim udtPepToProteinMapping() As udtPepToProteinMappingType

        Dim blnSuccess As Boolean

        If Not LoadParameterFileSettings(strParameterFilePath) Then
            SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, True)
            Return False
        End If

        Try
            If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
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
                        strInputFilePathFull = ioInputFile.FullName

                        ReDim udtMSGFDBModInfo(-1)
                        ReDim udtPepToProteinMapping(-1)

                        ' Load the MSGF-DB Parameter File so that we can determine the modification names and masses
                        ' The user is allowed to specify either the MSGF-DB parameter file or the MSGFDB_Mods.txt file
                        If Not ExtractModInfoFromMSGFDBModsFile(mSearchToolParameterFilePath, udtMSGFDBModInfo) Then
                            If udtMSGFDBModInfo Is Nothing OrElse udtMSGFDBModInfo.Length = 0 Then
                                ReDim udtMSGFDBModInfo(-1)
                            End If
                        End If

                        ' Resolve the mods in mInspectModInfo with the ModDefs mods
                        ResolveInspectModsWithModDefinitions(udtMSGFDBModInfo)

                        If MyBase.mCreateInspectOrMSGFDbFirstHitsFile Then

                            ' Create the first hits output file
                            MyBase.ResetProgress("Creating the FHT file")
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            strOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                            strOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strOutputFilePath & SEQUEST_FIRST_HITS_FILE_SUFFIX)

                            blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strOutputFilePath, udtMSGFDBModInfo, eFilteredOutputFileTypeConstants.FHTFile)

                        End If

                        If MyBase.mCreateInspectOrMSGFDbSynopsisFile Then

                            ' Create the synopsis output file
                            MyBase.ResetProgress("Creating the SYN file")
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            'Define the synopsis output file name based on strInputFilePath
                            strSynOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                            strSynOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strSynOutputFilePath & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                            blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, udtMSGFDBModInfo, eFilteredOutputFileTypeConstants.SynFile)

                            ' Load the PeptideToProteinMap information; if any modified peptides are present, then write out an updated _inspect_PepToProtMap.txt file with the new mod symbols (file will be named _PepToProtMapMTS.txt)
                            ' If the file doesn't exist, then a warning will be displayed, but processing will continue
                            strPepToProteinMapFilePath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(strInputFilePathFull), System.IO.Path.GetFileNameWithoutExtension(strInputFilePathFull) & FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING & ".txt")

                            MyBase.ResetProgress("Loading the PepToProtein map file: " & System.IO.Path.GetFileName(strPepToProteinMapFilePath))
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            LoadPeptideToProteinMapInfo(strPepToProteinMapFilePath, strOutputFolderPath, udtMSGFDBModInfo, udtPepToProteinMapping)


                            ' Create the other PHRP-specific files
                            MyBase.ResetProgress("Creating the PHRP files for " & System.IO.Path.GetFileName(strSynOutputFilePath))
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            ' Now parse the _syn.txt file that we just created to next create the other PHRP files
                            blnSuccess = ParseMSGFDBSynopsisFile(strSynOutputFilePath, udtPepToProteinMapping, False)

                        End If

                        If blnSuccess Then
                            MyBase.OperationComplete()
                        End If

                    Catch ex As Exception
                        SetErrorMessage("Error calling CreateFHTorSYNResultsFile: " & ex.Message)
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

    ''' <summary>
    ''' Replaces modification name text in peptide sequences with modification symbols (uses case-sensitive comparisons)
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="udtMSGFDBModInfo">This function assumes that each entry in udtMSGFDBModInfo() has both .ModName and .ModSymbol defined</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ReplaceMSGFModTextWithSymbol(ByVal strPeptide As String, _
                                                    ByRef udtMSGFDBModInfo() As udtModInfoType, _
                                                    ByRef dblTotalModMass As Double) As String

        Dim intIndex As Integer
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Dim strModSymbols As String = String.Empty

        Dim blnNterminalMod As Boolean
        Dim blnPossibleCTerminalMod As Boolean

        Static reNTerminalModMassRegEx As New System.Text.RegularExpressions.Regex(MSGFDB_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS)
        Static reCTerminalModMassRegEx As New System.Text.RegularExpressions.Regex(MSGFDB_CTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS)
        Static reModMassRegEx As New System.Text.RegularExpressions.Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS)

        Dim reMatch As System.Text.RegularExpressions.Match

        Dim dblModMassFound As Double

        ' Reset the total mod mass
        dblTotalModMass = 0

        If strPeptide.Length >= 4 Then
            If strPeptide.Chars(1) = "."c AndAlso _
               strPeptide.Chars(strPeptide.Length - 2) = "."c Then
                strPrefix = strPeptide.Substring(0, 2)
                strSuffix = strPeptide.Substring(strPeptide.Length - 2, 2)

                strPeptide = strPeptide.Substring(2, strPeptide.Length - 4)
            End If
        End If

        ' strPeptide should now be the clean peptide, without the prefix or suffix residues

        ' First look for N-terminal mods (NTermPeptide or NTermProtein)
        ' This RegEx will match one or more mods, all at the N-terminus
        reMatch = reNTerminalModMassRegEx.Match(strPeptide)
        If reMatch IsNot Nothing AndAlso reMatch.Success Then

            blnNterminalMod = True
            blnPossibleCTerminalMod = False

            ' Convert the mod mass (or masses) to one or more mod symbols
            If ConvertMGSFModMassesToSymbols(reMatch.Groups(1).Value, strModSymbols, udtMSGFDBModInfo, _
                                             blnNterminalMod, blnPossibleCTerminalMod, dblModMassFound) Then

                ' Replace the mod digits with the mod symbols

                strPeptide = ReplaceMSGFModTextWithMatchedSymbol(strPeptide, reMatch.Groups(1), strModSymbols)
                dblTotalModMass += dblModMassFound

            End If
        End If

        ' Next, step through the peptide and parse each mod mass that follows a residue
        ' Any mod mass at the end must be considered a C-terminal mod 

        ' Need to start at the first letter
        ' If we had N-terminal mods, they're currently notated like this: _.+42.011MDHTPQSQLK.L or _.+42.011+57.021MNDR.Q)
        ' We want things to look like this: -.#MDHTPQSQLK.L or -.#*MNDRQLNHR.S

        ' Static mods will not have a mod mass listed; we do not add mod symbols for static mods, but we do increment dblTotalModMass

        intIndex = 0
        Do While intIndex < strPeptide.Length AndAlso Not Char.IsLetter(strPeptide.Chars(intIndex))
            intIndex += 1
        Loop

        Do While intIndex < strPeptide.Length

            If Char.IsLetter(strPeptide.Chars(intIndex)) Then

                Dim objModificationDefinition As clsModificationDefinition

                For intModIndex As Integer = 0 To mPeptideMods.ModificationCount - 1
                    If mPeptideMods.GetModificationTypeByIndex(intModIndex) = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                        objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)

                        If objModificationDefinition.TargetResiduesContain(strPeptide.Chars(intIndex)) Then
                            ' Match found; update dblTotalModMass but do not add a static mod symbol
                            dblTotalModMass += objModificationDefinition.ModificationMass
                            Exit For
                        End If
                    End If
                Next intModIndex

                intIndex += 1

            Else
                ' Found a mod; find the extent of the mod digits
                reMatch = reModMassRegEx.Match(strPeptide, intIndex)

                blnNterminalMod = False
                If intIndex + reMatch.Groups(1).Value.Length >= strPeptide.Length Then
                    blnPossibleCTerminalMod = True
                Else
                    blnPossibleCTerminalMod = False
                End If

                ' Convert the mod mass (or masses) to one or more mod symbols
                If ConvertMGSFModMassesToSymbols(reMatch.Groups(1).Value, strModSymbols, udtMSGFDBModInfo, _
                                                 blnNterminalMod, blnPossibleCTerminalMod, dblModMassFound) Then

                    ' Replace the mod digits with the mod symbols

                    strPeptide = ReplaceMSGFModTextWithMatchedSymbol(strPeptide, reMatch.Groups(1), strModSymbols)
                    dblTotalModMass += dblModMassFound

                    intIndex += strModSymbols.Length

                Else
                    intIndex += reMatch.Groups(1).Value.Length
                End If

            End If

        Loop

        ' If any N-terminal mods were present, we need to move them to after the first residue
        ' in other words, switch from #MDHTPQSQLK to M#DHTPQSQLK
        '                          or #*MNDRQLNHR to M#*NDRQLNHR

        Dim intIndexFirstLetter As Integer = 0
        Do While intIndexFirstLetter <= strPeptide.Length AndAlso Not Char.IsLetter(strPeptide.Chars(intIndexFirstLetter))
            intIndexFirstLetter += 1
        Loop

        If intIndexFirstLetter > 0 AndAlso intIndexFirstLetter < strPeptide.Length Then
            Dim strPeptideNew As String
            strPeptideNew = strPeptide.Chars(intIndexFirstLetter) & strPeptide.Substring(0, intIndexFirstLetter)
            If intIndexFirstLetter < strPeptide.Length - 1 Then
                strPeptideNew &= strPeptide.Substring(intIndexFirstLetter + 1)
            End If
            strPeptide = String.Copy(strPeptideNew)
        End If

        Return strPrefix & strPeptide & strSuffix

    End Function


    Protected Function ReplaceMSGFModTextWithMatchedSymbol(ByVal strPeptide As String, _
                                                           ByRef reGroup As System.Text.RegularExpressions.Group, _
                                                           ByVal strModSymbols As String) As String
        Dim strPeptideNew As String

        If reGroup.Index > 0 Then
            strPeptideNew = strPeptide.Substring(0, reGroup.Index)
        Else
            strPeptideNew = String.Empty
        End If

        strPeptideNew &= strModSymbols

        If reGroup.Index + reGroup.Length < strPeptide.Length Then
            strPeptideNew &= strPeptide.Substring(reGroup.Index + reGroup.Length)
        End If

        Return strPeptideNew

    End Function

    Protected Function ReplaceTerminus(ByVal strPeptide As String) As String

        If strPeptide.StartsWith(N_TERMINUS_SYMBOL_MSGFDB) Then
            strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST & "." & strPeptide.Substring(N_TERMINUS_SYMBOL_MSGFDB.Length)
        End If

        If strPeptide.EndsWith(C_TERMINUS_SYMBOL_MSGFDB) Then
            strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_MSGFDB.Length) & "." & clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
        End If

        Return strPeptide

    End Function

    Protected Sub ResolveInspectModsWithModDefinitions(ByRef udtMSGFDBModInfo() As udtModInfoType)

        Dim intIndex As Integer
        Dim intResidueIndex As Integer
        Dim intResIndexStart As Integer
        Dim intResIndexEnd As Integer

        Dim chTargetResidue As Char
        Dim eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants
        Dim blnExistingModFound As Boolean

        Dim objModificationDefinition As clsModificationDefinition

        If Not udtMSGFDBModInfo Is Nothing Then
            For intIndex = 0 To udtMSGFDBModInfo.Length - 1
                ' Call .LookupModificationDefinitionByMass for each entry in udtMSGFDBModInfo

                With udtMSGFDBModInfo(intIndex)

                    If .Residues.Length > 0 Then
                        intResIndexStart = 0
                        intResIndexEnd = .Residues.Length - 1
                    Else
                        intResIndexStart = -1
                        intResIndexEnd = -1
                    End If

                    For intResidueIndex = intResIndexStart To intResIndexEnd
                        If intResidueIndex >= 0 Then
                            chTargetResidue = .Residues.Chars(intResidueIndex)
                            If chTargetResidue = "*"c Then
                                ' This is a terminal mod, and MSGFDB lists the target residue as * for terminal mods
                                ' This program requires that chTargetResidue be Nothing
                                chTargetResidue = Nothing
                            End If
                        Else
                            chTargetResidue = Nothing
                        End If

                        If .ModType = eMSGFDBModType.DynNTermPeptide Then
                            eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideNTerminus
                        ElseIf .ModType = eMSGFDBModType.DynCTermPeptide Then
                            eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideCTerminus
                        ElseIf .ModType = eMSGFDBModType.DynNTermProtein Then
                            eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinNTerminus
                        ElseIf .ModType = eMSGFDBModType.DynCTermProtein Then
                            eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinCTerminus
                        Else
                            eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.None
                        End If
                        blnExistingModFound = False

                        objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(.ModMassVal, chTargetResidue, eResidueTerminusState, blnExistingModFound, True)

                        If intResidueIndex = intResIndexStart Then
                            .ModSymbol = objModificationDefinition.ModificationSymbol
                        End If

                    Next intResidueIndex


                End With

            Next intIndex
        End If
    End Sub

    Protected Sub StoreSearchResult(ByRef swResultFile As System.IO.StreamWriter, _
                                    ByRef intResultID As Integer, _
                                    ByRef udtSearchResult As udtMSGFDBSearchResultType, _
                                    ByRef intFilteredSearchResultCount As Integer, _
                                    ByRef udtFilteredSearchResults() As udtMSGFDBSearchResultType, _
                                    ByRef strErrorLog As String)

        If intFilteredSearchResultCount = udtFilteredSearchResults.Length Then
            ReDim Preserve udtFilteredSearchResults(udtFilteredSearchResults.Length * 2 - 1)
        End If

        udtFilteredSearchResults(intFilteredSearchResultCount) = udtSearchResult
        intFilteredSearchResultCount += 1

    End Sub

    Private Sub SortAndWriteFilteredSearchResults(ByRef swResultFile As System.IO.StreamWriter, _
                                                  ByVal intFilteredSearchResultCount As Integer, _
                                                  ByRef udtFilteredSearchResults() As udtMSGFDBSearchResultType, _
                                                  ByRef strErrorLog As String, _
                                                  ByVal blnIncludeFDRandPepFDR As Boolean)

        Dim intIndex As Integer

        ' Sort udtFilteredSearchResults by ascending SpecProb, ascending scan, ascending charge, ascending peptide, and ascending protein
        Array.Sort(udtFilteredSearchResults, 0, intFilteredSearchResultCount, New MSGFDBSearchResultsComparerSpecProbScanChargePeptide)

        For intIndex = 0 To intFilteredSearchResultCount - 1
            WriteSearchResultToFile(intIndex + 1, swResultFile, udtFilteredSearchResults(intIndex), strErrorLog, blnIncludeFDRandPepFDR)
        Next intIndex

    End Sub

    Private Sub StoreTopFHTMatch(ByRef swResultFile As System.IO.StreamWriter, _
                                 ByRef intResultID As Integer, _
                                 ByRef udtSearchResults() As udtMSGFDBSearchResultType, _
                                 ByVal intStartIndex As Integer, _
                                 ByVal intEndIndex As Integer, _
                                 ByRef intFilteredSearchResultCount As Integer, _
                                 ByRef udtFilteredSearchResults() As udtMSGFDBSearchResultType, _
                                 ByRef strErrorLog As String)

        Dim intIndex As Integer
        Dim intCurrentCharge As Short = Short.MinValue

        AssignRankAndDeltaNormValues(udtSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure already sorted by scan, charge, and SpecProb; no need to re-sort

        ' Now store or write out the first match for each charge for this scan
        For intIndex = intStartIndex To intEndIndex
            If intIndex = intStartIndex OrElse intCurrentCharge <> udtSearchResults(intIndex).ChargeNum Then
                StoreSearchResult(swResultFile, intResultID, udtSearchResults(intIndex), intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                intCurrentCharge = udtSearchResults(intIndex).ChargeNum
            End If
        Next intIndex

    End Sub

    Private Sub StoreSynMatches(ByRef swResultFile As System.IO.StreamWriter, _
                                ByRef intResultID As Integer, _
                                ByRef udtSearchResults() As udtMSGFDBSearchResultType, _
                                ByVal intStartIndex As Integer, _
                                ByVal intEndIndex As Integer, _
                                ByRef intFilteredSearchResultCount As Integer, _
                                ByRef udtFilteredSearchResults() As udtMSGFDBSearchResultType, _
                                ByRef strErrorLog As String)

        Dim intIndex As Integer
        Dim intCurrentCharge As Short = Short.MinValue

        AssignRankAndDeltaNormValues(udtSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure already sorted by scan, charge, and SpecProb; no need to re-sort

        ' Now store or write out the matches that pass the filters
        For intIndex = intStartIndex To intEndIndex
            If udtSearchResults(intIndex).PValueNum <= mMSGFDBSynopsisFilePValueThreshold OrElse _
               udtSearchResults(intIndex).SpecProbNum <= mMSGFDBSynopsisFileSpecProbThreshold Then
                StoreSearchResult(swResultFile, intResultID, udtSearchResults(intIndex), intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
            End If
        Next intIndex

    End Sub

    ''' <summary>
    ''' Return the text up to (but not including) the first space in strProteinNameAndDescription
    ''' </summary>
    ''' <param name="strProteinNameAndDescription"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function TruncateProteinName(ByVal strProteinNameAndDescription As String) As String

        Dim intIndex As Integer

        intIndex = strProteinNameAndDescription.IndexOf(" "c)
        If intIndex > 0 Then
            Return strProteinNameAndDescription.Substring(0, intIndex)
        Else
            Return strProteinNameAndDescription
        End If

    End Function

    Private Sub UpdateSearchResultEnzymeAndTerminusInfo(ByRef objSearchResult As clsSearchResultsMSGFDB)
        With objSearchResult
            .SetEnzymeMatchSpec(mEnzymeMatchSpec)

            ' Update the N-Terminus and/or C-Terminus masses if those in the XML file are significantly different than the defaults
            If mPeptideNTerminusMassChange <> 0 Then
                .UpdatePeptideNTerminusMass(mPeptideNTerminusMassChange)
            End If

            If mPeptideCTerminusMassChange <> 0 Then
                .UpdatePeptideCTerminusMass(mPeptideCTerminusMassChange)
            End If
        End With
    End Sub

    Private Function ValidateIsValidModSpec(ByVal strLineIn As String, ByVal strModTag As String) As String

        Dim strModSpec As String = String.Empty

        If strLineIn.ToLower.StartsWith(strModTag.ToLower()) Then
            strModSpec = strLineIn.Substring(strModTag.Length)

            If strModSpec.StartsWith("="c) Then
                ' Mod spec found
                strModSpec = strModSpec.Substring(1)

                If strModSpec.ToLower() = "none" Then
                    ' None means no-mods
                    strModSpec = String.Empty
                End If

            Else
                strModSpec = String.Empty
            End If
        End If

        Return strModSpec

    End Function

    Private Sub WriteSynFHTFileHeader(ByRef swResultFile As System.IO.StreamWriter, _
                                      ByRef strErrorLog As String, _
                                      ByVal blnIncludeFDRandPepFDR As Boolean)

        ' Write out the header line for synopsis / first hits files
        Try
            Dim strSuffix As String = String.Empty
            If blnIncludeFDRandPepFDR Then
                strSuffix = ControlChars.Tab & "FDR" & ControlChars.Tab & "PepFDR"
            End If

            swResultFile.WriteLine("ResultID" & ControlChars.Tab & _
                                   "Scan" & ControlChars.Tab & _
                                   "FragMethod" & ControlChars.Tab & _
                                   "SpecIndex" & ControlChars.Tab & _
                                   "Charge" & ControlChars.Tab & _
                                   "PrecursorMZ" & ControlChars.Tab & _
                                   "DelM" & ControlChars.Tab & _
                                   "DelM_PPM" & ControlChars.Tab & _
                                   "MH" & ControlChars.Tab & _
                                   "Peptide" & ControlChars.Tab & _
                                   "Protein" & ControlChars.Tab & _
                                   "NTT" & ControlChars.Tab & _
                                   "DeNovoScore" & ControlChars.Tab & _
                                   "MSGFScore" & ControlChars.Tab & _
                                   "MSGFDB_SpecProb" & ControlChars.Tab & _
                                   "Rank_MSGFDB_SpecProb" & ControlChars.Tab & _
                                   "PValue" & _
                                   strSuffix)


        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits header" & ControlChars.NewLine
            End If
        End Try

    End Sub

    Private Sub WriteSearchResultToFile(ByVal intResultID As Integer, _
                                        ByRef swResultFile As System.IO.StreamWriter, _
                                        ByRef udtSearchResult As udtMSGFDBSearchResultType, _
                                        ByRef strErrorLog As String, _
                                        ByVal blnIncludeFDRandPepFDR As Boolean)

        ' Writes an entry to a synopsis or first hits file
        Try
            ' Columns:
            ' ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecProb  Rank_MSGFDB_SpecProb  PValue  FDR  PepFDR

            Dim strSuffix As String = String.Empty
            If blnIncludeFDRandPepFDR Then
                strSuffix = ControlChars.Tab & udtSearchResult.FDR & ControlChars.Tab & udtSearchResult.PepFDR
            End If

            swResultFile.WriteLine(intResultID.ToString & ControlChars.Tab & _
                                   udtSearchResult.Scan & ControlChars.Tab & _
                                   udtSearchResult.FragMethod & ControlChars.Tab & _
                                   udtSearchResult.SpecIndex & ControlChars.Tab & _
                                   udtSearchResult.Charge & ControlChars.Tab & _
                                   udtSearchResult.PrecursorMZ & ControlChars.Tab & _
                                   udtSearchResult.PMErrorDa & ControlChars.Tab & _
                                   udtSearchResult.PMErrorPPM & ControlChars.Tab & _
                                   udtSearchResult.MH & ControlChars.Tab & _
                                   udtSearchResult.Peptide & ControlChars.Tab & _
                                   udtSearchResult.Protein & ControlChars.Tab & _
                                   udtSearchResult.NTT & ControlChars.Tab & _
                                   udtSearchResult.DeNovoScore & ControlChars.Tab & _
                                   udtSearchResult.MSGFScore & ControlChars.Tab & _
                                   udtSearchResult.SpecProb & ControlChars.Tab & _
                                   udtSearchResult.RankSpecProb & ControlChars.Tab & _
                                   udtSearchResult.PValue & _
                                   strSuffix)

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

#Region "IComparer Classes"

    Protected Class MSGFDBSearchResultsComparerScanChargeSpecProbPeptide
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtMSGFDBSearchResultType = DirectCast(x, udtMSGFDBSearchResultType)
            Dim yData As udtMSGFDBSearchResultType = DirectCast(y, udtMSGFDBSearchResultType)

            If xData.ScanNum > yData.ScanNum Then
                Return 1
            ElseIf xData.ScanNum < yData.ScanNum Then
                Return -1
            Else
                ' Scan is the same, check charge
                If xData.ChargeNum > yData.ChargeNum Then
                    Return 1
                ElseIf xData.ChargeNum < yData.ChargeNum Then
                    Return -1
                Else
                    ' Charge is the same; check SpecProb
                    If xData.SpecProbNum > yData.SpecProbNum Then
                        Return 1
                    ElseIf xData.SpecProbNum < yData.SpecProbNum Then
                        Return -1
                    Else
                        ' SpecProb is the same; check peptide
                        If xData.Peptide > yData.Peptide Then
                            Return 1
                        ElseIf xData.Peptide < yData.Peptide Then
                            Return -1
                        Else
                            ' Peptide is the same, check Protein
                            If xData.Protein > yData.Protein Then
                                Return 1
                            ElseIf xData.Protein < yData.Protein Then
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
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtMSGFDBSearchResultType = DirectCast(x, udtMSGFDBSearchResultType)
            Dim yData As udtMSGFDBSearchResultType = DirectCast(y, udtMSGFDBSearchResultType)

            If xData.SpecProbNum > yData.SpecProbNum Then
                Return 1
            ElseIf xData.SpecProbNum < yData.SpecProbNum Then
                Return -1
            Else
                ' SpecProbNum is the same; check scan number
                If xData.ScanNum > yData.ScanNum Then
                    Return 1
                ElseIf xData.ScanNum < yData.ScanNum Then
                    Return -1
                Else
                    ' Scan is the same, check charge
                    If xData.ChargeNum > yData.ChargeNum Then
                        Return 1
                    ElseIf xData.ChargeNum < yData.ChargeNum Then
                        Return -1
                    Else
                        ' Charge is the same; check peptide
                        If xData.Peptide > yData.Peptide Then
                            Return 1
                        ElseIf xData.Peptide < yData.Peptide Then
                            Return -1
                        Else
                            ' Peptide is the same, check Protein
                            If xData.Protein > yData.Protein Then
                                Return 1
                            ElseIf xData.Protein < yData.Protein Then
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

    Protected Class PepToProteinMappingComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtPepToProteinMappingType = DirectCast(x, udtPepToProteinMappingType)
            Dim yData As udtPepToProteinMappingType = DirectCast(y, udtPepToProteinMappingType)

            If xData.Peptide > yData.Peptide Then
                Return 1
            ElseIf xData.Peptide < yData.Peptide Then
                Return -1
            Else
                If xData.Protein > yData.Protein Then
                    Return 1
                ElseIf xData.Protein < yData.Protein Then
                    Return -1
                Else
                    Return 0
                End If
            End If

        End Function
    End Class

    Protected Class PepToProteinMappingPeptideSearchComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtPepToProteinMappingType = DirectCast(x, udtPepToProteinMappingType)
            Dim strPeptide As String = DirectCast(y, String)

            If xData.Peptide > strPeptide Then
                Return 1
            ElseIf xData.Peptide < strPeptide Then
                Return -1
            Else

                Return 0
            End If

        End Function
    End Class

#End Region

End Class