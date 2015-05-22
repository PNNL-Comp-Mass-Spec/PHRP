Option Strict On

' This class reads in an MODPlus results file (txt format) and creates 
' a tab-delimited text file with the data. 
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 05/15/2015
'
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports PHRPReader
Imports System.IO
Imports System.Runtime.InteropServices
Imports System.Text.RegularExpressions
Imports System.Xml

Public Class clsMODPlusResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "May 22, 2015"
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MODPlus_FILE As String = "_modp.id"

    Public Const N_TERMINUS_SYMBOL_MODPlus As String = "-"
    Public Const C_TERMINUS_SYMBOL_MODPlus As String = "-"

    ' This is used for filtering both MODa and MODPlus results
    Public Const DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD As Single = 0.05

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Private Const MODPlus_MOD_MASS_REGEX As String = "([+-][0-9.]+)"

    Private Const MODPlus_MASS_DIGITS_OF_PRECISION As Byte = 3

    Private Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline Or RegexOptions.IgnoreCase

    ' These columns correspond to the tab-delimited file (_MODPlus.id.txt) created by MODPlus's tda_plus.jar file
    Protected Const MODPlusResultsFileColCount As Integer = 13
    Public Enum eMODPlusResultsFileColumns As Integer
        SpectrumFileName = 0
        SpectrumIndex = 1
        ScanNumber = 2
        ObservedMonoMass = 3
        Charge = 4
        CalculatedMonoMass = 5
        DeltaMass = 6
        Score = 7
        Probability = 8
        Peptide = 9
        NTT = 10
        ProteinAndPeptidePositionList = 11
        ModificationAnnotation = 12
    End Enum

    ' These columns correspond to the Synopsis file created by this class
    Protected Const MODPlusSynFileColCount As Integer = 17
    Public Enum eMODPlusSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        Spectrum_Index = 2
        Charge = 3
        PrecursorMZ = 4
        DelM = 5                            ' Precursor error, in Da
        DelM_PPM = 6                        ' Precursor error, in ppm
        MH = 7                              ' Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
        Peptide = 8                         ' This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
        NTT = 9
        ModificationAnnotation = 10
        Protein = 11
        Peptide_Position = 12
        Score = 13
        Probability = 14
        Rank_Score = 15
        QValue = 16
    End Enum

#End Region

#Region "Structures"
    ' This data structure holds rows read from the tab-delimited file (_MODPlus.id.txt) created by MODPlus's tda_plus.jar file
    Protected Structure udtMODPlusSearchResultType
        Public SpectrumFileName As String
        Public SpectrumIndex As String
        Public ScanNum As Integer
        Public Precursor_mass As String         ' Uncharged monoisotopic mass value of the observed precursor_mz, reported as ObservedMW by MODPlus 
        Public PrecursorMZ As String            ' Computed by this class from ObservedMonoMass
        Public Charge As String
        Public ChargeNum As Short
        Public CalculatedMonoMass As String     ' Theoretical monoisotopic mass of the peptide (including mods), as computed by MODPlus
        Public DeltaMass As String              ' Computed by MODPlus
        Public MH As String                     ' Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
        Public DelM As String                   ' Computed by this class using Precursor_mass - CalculatedMonoMass
        Public DelM_PPM As String               ' Computed by this class using DelM and CalculatedMonoMass	
        Public Score As String
        Public ScoreNum As Double
        Public Probability As String            ' Higher values are better
        Public ProbabilityNum As Double         ' Higher values are better
        Public RankScore As Integer
        Public Peptide As String
        Public NTT As String                    ' Number of Tryptic Terminii
        Public ProteinList As String            ' One or more protein entries of the form ref|YP_003651515.1[K.196~206.Q(2)] where the text in brackets is the start/stop residues of the peptide; Multiple entries will be separated by semicolons, e.g. ref|YP_003651515.1[K.196~206.Q(2)];ref|YP_003201491.1[K.223~233.Q(2)];ref|YP_003313784.1[K.266~276.Q(2)]
        Public ModificationAnnotation As String
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
            ScoreNum = 0
            Probability = String.Empty
            ProbabilityNum = 0
            RankScore = 0
            Peptide = String.Empty
            NTT = String.Empty
            ProteinList = String.Empty
            ModificationAnnotation = String.Empty
            FDR = 0
            QValue = 0
        End Sub

        Public Overrides Function ToString() As String
            Return Probability & ": " & Peptide
        End Function

    End Structure

#End Region

#Region "Classwide Variables"
    Protected mProteinNamePositionSplit As Regex
#End Region

    ''' <summary>
    ''' Step through .PeptideSequenceWithMods
    ''' For each residue, check if a static mod is defined that affects that residue
    ''' For each mod mass, determine the modification and add to objSearchResult
    ''' </summary>
    ''' <param name="objSearchResult"></param>
    ''' <param name="blnUpdateModOccurrenceCounts"></param>
    ''' <remarks></remarks>
    Private Sub AddDynamicAndStaticResidueMods(ByVal objSearchResult As clsSearchResultsMODPlus, ByVal blnUpdateModOccurrenceCounts As Boolean)
        Const NO_RESIDUE As Char = "-"c

        Dim intIndex As Integer, intModIndex As Integer
        Dim chChar As Char
        Dim objModificationDefinition As clsModificationDefinition

        Dim strSequence As String

        Dim blnParsingModMass As Boolean
        Dim strModMassDigits As String = String.Empty

        Dim chMostRecentResidue As Char
        Dim intResidueLocInPeptide As Integer

        blnParsingModMass = False
        strModMassDigits = String.Empty

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

    Private Function AddModificationsAndComputeMass(ByVal objSearchResult As clsSearchResultsMODPlus, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
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

    Private Sub AssociateDynamicModWithResidue(
      ByVal objSearchResult As clsSearchResultsMODPlus,
      ByVal chMostRecentResidue As Char,
      ByVal intResidueLocInPeptide As Integer,
      ByVal strModMassDigits As String,
      ByVal blnUpdateModOccurrenceCounts As Boolean)

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

            blnSuccess = objSearchResult.SearchResultAddModification(dblModMass, chResidueForMod, intResidueLocForMod, objSearchResult.DetermineResidueTerminusState(intResidueLocForMod), blnUpdateModOccurrenceCounts, MODPlus_MASS_DIGITS_OF_PRECISION)

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
    ''' Ranks each entry assumes all of the data is from the same scan)
    ''' </summary>
    ''' <param name="lstSearchResults"></param>
    ''' <param name="intStartIndex"></param>
    ''' <param name="intEndIndex"></param>
    ''' <remarks></remarks>
    Private Sub AssignRankAndDeltaNormValues(
      ByVal lstSearchResults As List(Of udtMODPlusSearchResultType),
      ByVal intStartIndex As Integer,
      ByVal intEndIndex As Integer)

        ' Prior to September 2014 ranks were assign per charge state per scan; 
        ' Ranks are now assigned per scan (across all charge states)

        ' Duplicate a portion of lstSearchResults so that we can sort by descending Probability

        Dim dctResultsSubset = New Dictionary(Of Integer, udtMODPlusSearchResultType)
        For intIndex = intStartIndex To intEndIndex
            dctResultsSubset.Add(intIndex, lstSearchResults(intIndex))
        Next

        Dim lstResultsByScore = (From item In dctResultsSubset Select item Order By item.Value.ScoreNum Descending).ToList()

        Dim dblLastValue As Double
        Dim intCurrentRank As Integer = -1

        For Each entry In lstResultsByScore
            Dim oResult = lstSearchResults(entry.Key)

            If intCurrentRank < 0 Then
                dblLastValue = oResult.ScoreNum
                intCurrentRank = 1
            Else
                If Math.Abs(oResult.ScoreNum - dblLastValue) > Double.Epsilon Then
                    dblLastValue = oResult.ScoreNum
                    intCurrentRank += 1
                End If
            End If

            oResult.RankScore = intCurrentRank
            lstSearchResults(entry.Key) = oResult
        Next

    End Sub

    Protected Function AssureInteger(ByVal strInteger As String, ByVal intDefaultValue As Integer) As String

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

    Protected Function ComputePeptideMass(ByVal strPeptide As String, ByVal dblTotalModMass As Double) As Double

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
    Protected Function ComputeTotalModMass(ByVal strPeptide As String) As Double

        Static reModMassRegEx As New Regex(MODPlus_MOD_MASS_REGEX, REGEX_OPTIONS)

        Dim dblTotalModMass As Double = 0

        Dim strPrimarySequence As String = String.Empty
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, strPrimarySequence, strPrefix, strSuffix)

        ' Parse the dynamic mods reported by MODPlus
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

    Protected Overrides Function ConstructPepToProteinMapFilePath(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal MTS As Boolean) As String

        Dim strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath)
        If strPepToProteinMapFilePath.ToLower().EndsWith("_MODPlus_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_MODPlus_fht") Then
            ' Remove _syn or _fht
            strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4)
        End If

        Return MyBase.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS)

    End Function

    ''' <summary>
    ''' This routine creates a first hits file or synopsis file from the output from MODPlus
    ''' The synopsis file includes every result with a probability above a set threshold
    ''' The first-hits file includes the result with the highest probability (for each scan and charge)
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFilePath"></param>
    ''' <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CreateSynResultsFile(
       strInputFilePath As String,
       strOutputFilePath As String,
       Optional blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean

        Try
            Dim intColumnMapping() As Integer = Nothing
            Dim strErrorLog = String.Empty

            ' Open the input file and parse it
            ' Initialize the stream reader and the stream Text writer
            Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.Read)),
                  swResultFile = New StreamWriter(New FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                Dim headerParsed = False

                ' Initialize the list that will hold all of the records in the MODPlus result file
                Dim lstSearchResultsUnfiltered = New List(Of udtMODPlusSearchResultType)

                ' Initialize the list that will hold all of the records that will ultimately be written out to disk
                Dim lstFilteredSearchResults = New List(Of udtMODPlusSearchResultType)

                ' Parse the input file
                Do While Not srDataFile.EndOfStream And Not MyBase.AbortProcessing
                    Dim strLineIn = srDataFile.ReadLine()

                    If String.IsNullOrWhiteSpace(strLineIn) Then
                        Continue Do
                    End If

                    If Not headerParsed Then
                        ' Parse the header line

                        Dim blnsuccess = ParseMODPlusResultsFileHeaderLine(strLineIn, intColumnMapping)
                        If Not blnsuccess Then
                            If String.IsNullOrEmpty(mErrorMessage) Then
                                SetErrorMessage("Invalid header line in " & Path.GetFileName(strInputFilePath))
                            End If
                            Return False
                        End If

                        ' Write the header line to the output file
                        WriteSynFHTFileHeader(swResultFile, strErrorLog)

                        headerParsed = True
                        Continue Do
                    End If


                    Dim udtSearchResult = New udtMODPlusSearchResultType

                    Dim blnValidSearchResult = ParseMODPlusResultsFileEntry(strLineIn, udtSearchResult, strErrorLog, intColumnMapping)

                    If blnValidSearchResult Then
                        lstSearchResultsUnfiltered.Add(udtSearchResult)
                    End If

                    ' Update the progress
                    Dim sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
                    If mCreateProteinModsFile Then
                        sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
                    End If
                    UpdateProgress(sngPercentComplete)

                Loop

                ' Sort the SearchResults by scan, charge, and descending score
                lstSearchResultsUnfiltered.Sort(New MODPlusSearchResultsComparerScanChargeScorePeptide)

                ' Now filter the data

                ' Initialize variables
                Dim intStartIndex = 0
                Dim intEndIndex As Integer

                intStartIndex = 0
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
    ''' Load the static mods defined in the MODPlus parameter file
    ''' </summary>
    ''' <param name="strMODPlusParamFilePath"></param>
    ''' <param name="lstModInfo"></param>
    ''' <returns></returns>
    ''' <remarks>We don't care about the dynamic mods because there are so many possible mods.  We'll add each dynamic mod as we encounter it in the results</remarks>
    Protected Function ExtractModInfoFromMODPlusParamFile(ByVal strMODPlusParamFilePath As String, ByRef lstModInfo As List(Of clsModificationDefinition)) As Boolean

        Try
            ' Initialize the modification list
            If lstModInfo Is Nothing Then
                lstModInfo = New List(Of clsModificationDefinition)
            Else
                lstModInfo.Clear()
            End If

            If String.IsNullOrEmpty(strMODPlusParamFilePath) Then
                SetErrorMessage("MODPlus Parameter File name not defined; unable to extract mod info")
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                Return False
            End If

            Dim fiParamFile = New FileInfo(strMODPlusParamFilePath)
            If Not fiParamFile.Exists Then
                SetErrorMessage("MODPlus param file not found: " & strMODPlusParamFilePath)
            End If

            ' Read the contents of the parameter file
            Dim doc = New XmlDocument()
            doc.Load(fiParamFile.FullName)

            Dim nodeList = doc.SelectNodes("/search/modifications/fixed/mod")
            If nodeList.Count > 0 Then
                ' Store the fixed mods

                For Each node As XmlNode In nodeList
                    Dim modName = node.Attributes("name").Value
                    Dim strResidue = node.Attributes("site").Value.Trim()
                    Dim modPosition = node.Attributes("position").Value
                    Dim modMass = node.Attributes("massdiff").Value

                    ' Replace N-Term or C-Term with < or >
                    If strResidue.ToLower() = "n-term" Then strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS
                    If strResidue.ToLower() = "c-term" Then strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS

                    Dim modMassDa As Double
                    If Double.TryParse(modMass, modMassDa) Then

                        If Math.Abs(modMassDa - 0) > Single.Epsilon Then

                            Dim strMassCorrectionTag As String = mPeptideMods.LookupMassCorrectionTagByMass(modMassDa)

                            Dim eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod
                            If strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS OrElse strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                                eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                            End If

                            Dim objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, strResidue, eModType, strMassCorrectionTag)
                            lstModInfo.Add(objModDef)

                        End If
                    End If

                Next

            End If

            Console.WriteLine()

            Return True

        Catch ex As Exception
            SetErrorMessage("Error reading the MODPlus parameter file (" & Path.GetFileName(strMODPlusParamFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            Return False
        End Try

    End Function

    ''' <summary>
    ''' Parse the Synopsis file to create the other PHRP-compatible files
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFolderPath"></param>
    ''' <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ParseMODPlusSynopsisFile(
      ByVal strInputFilePath As String,
      ByVal strOutputFolderPath As String,
      ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean) As Boolean

        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim strPreviousProbability As String

        ' Note that MODPlus synopsis files are normally sorted on Probability value, ascending
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
            Dim objSearchResult = New clsSearchResultsMODPlus(mPeptideMods)

            ' Initialize htPeptidesFoundForProbabilityLevel
            Dim htPeptidesFoundForProbabilityLevel = New Hashtable
            Dim blnFirstMatchForGroup As Boolean

            strPreviousProbability = String.Empty

            Dim strErrorLog = String.Empty

            Try
                objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(mEnzymeMatchSpec, mPeptideNTerminusMassChange, mPeptideCTerminusMassChange)

                ' Open the input file and parse it
                ' Initialize the stream reader
                Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.Read))

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
                            blnSuccess = ParseMODPlusSynFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True
                            Continue Do
                        End If

                        Dim strCurrentPeptideWithMods As String = String.Empty

                        Dim blnValidSearchResult = ParseMODPlusSynFileEntry(
                          strLineIn, objSearchResult, strErrorLog,
                          intResultsProcessed, intColumnMapping,
                          strCurrentPeptideWithMods)

                        If Not blnValidSearchResult Then
                            Continue Do
                        End If

                        Dim strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge

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
                            ' Reset htPeptidesFoundForScan
                            htPeptidesFoundForProbabilityLevel.Clear()

                            ' Update strPreviousProbability
                            strPreviousProbability = objSearchResult.Probability

                            ' Append a new entry to htPeptidesFoundForScan
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

    ''' <summary>
    ''' Parse a MODPlus results line while creating the MODPlus sysnopsis file
    ''' </summary>
    ''' <param name="strLineIn"></param>
    ''' <param name="udtSearchResult"></param>
    ''' <param name="strErrorLog"></param>
    ''' <param name="intColumnMapping"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function ParseMODPlusResultsFileEntry(
      ByVal strLineIn As String,
      ByRef udtSearchResult As udtMODPlusSearchResultType,
      ByRef strErrorLog As String,
      ByVal intColumnMapping() As Integer) As Boolean

        ' Parses an entry from the MODPlus results file

        Dim rowIndex = "?"
        Dim dblPrecursorMonoMass As Double          ' Observed m/z, converted to monoisotopic mass
        Dim dblPeptideMonoMassMODPlus As Double     ' Theoretical peptide monoisotopic mass, including mods, as computed by MODPlus
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
                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.SpectrumFileName), .SpectrumFileName)

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.SpectrumIndex), .SpectrumIndex) Then
                        Throw New EvaluateException("Index column is missing or invalid")
                    Else
                        rowIndex = .SpectrumIndex
                    End If

                    Dim spectrumIndex As Integer
                    If Not Integer.TryParse(.SpectrumIndex, spectrumIndex) Then
                        Throw New EvaluateException("Index column is not numeric")
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.ScanNumber), .ScanNum)

                    ' Monoisotopic mass value of the observed precursor_mz
                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.ObservedMonoMass), .Precursor_mass)
                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.Charge), .Charge)
                    .ChargeNum = CShort(CIntSafe(.Charge, 0))

                    If Double.TryParse(.Precursor_mass, dblPrecursorMonoMass) Then
                        If .ChargeNum > 0 Then
                            dblPrecursorMZ = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMonoMass, 0, .ChargeNum)
                            .PrecursorMZ = NumToString(dblPrecursorMZ, 6, True)
                        End If
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.CalculatedMonoMass), .CalculatedMonoMass)
                    Double.TryParse(.CalculatedMonoMass, dblPeptideMonoMassMODPlus)

                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.DeltaMass), .DeltaMass)
                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.Score), .Score)
                    If Not Double.TryParse(.Score, .ScoreNum) Then .ScoreNum = 0

                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.Probability), .Probability)
                    If Not Double.TryParse(.Probability, .ProbabilityNum) Then .ProbabilityNum = 0

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.Peptide), .Peptide) Then
                        Throw New EvaluateException("Peptide column is missing or invalid")
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.NTT), .NTT)

                    If strSplitLine.Length > eMODPlusResultsFileColumns.ProteinAndPeptidePositionList Then
                        GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.ProteinAndPeptidePositionList), .ProteinList)

                        ' The protein column will have both the protein name and the peptide position
                        ' For example, ref|YP_001038741.1[R.67~78.L(2)]
                        ' It may have multiple proteins listed, separated by semicolons
                        ' We will split the list on semicolons in function ParseMODPlusSynFileEntry

                        If Not .ProteinList.Contains("["c) Then
                            ' This is likely a reverse-hit protein
                            .ModificationAnnotation = String.Copy(.ProteinList)
                            .ProteinList = String.Empty
                        Else
                            If strSplitLine.Length > eMODPlusResultsFileColumns.ModificationAnnotation Then
                                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusResultsFileColumns.ModificationAnnotation), .ModificationAnnotation)
                            End If
                        End If
                    End If

                    ' Parse the sequence to determine the total mod mass
                    ' Note that we do not remove any of the mod symbols since MODPlus identifies mods by mass alone
                    ' Note that static mods are implied (thus are not explicitly displayed by MODPlus)
                    dblTotalModMass = ComputeTotalModMass(.Peptide)

                    ' Compute monoisotopic mass of the peptide
                    dblPeptideMonoMassPHRP = ComputePeptideMass(.Peptide, dblTotalModMass)

                    ' Only override dblPeptideMonoMassMODPlus if it is 0
                    If Math.Abs(dblPeptideMonoMassMODPlus) < Double.Epsilon Then
                        dblPeptideMonoMassMODPlus = dblPeptideMonoMassPHRP
                    End If

                    Dim dblMassDiffThreshold As Double = dblPeptideMonoMassMODPlus / 50000
                    If dblMassDiffThreshold < 0.1 Then dblMassDiffThreshold = 0.1

                    If Math.Abs(dblPeptideMonoMassPHRP - dblPeptideMonoMassMODPlus) > dblMassDiffThreshold Then
                        ' Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da or by a slightly larger value if over 5000 Da; this is unexpected
                        Dim strFirst30Residues As String
                        If .Peptide.Length < 27 Then
                            strFirst30Residues = .Peptide
                        Else
                            strFirst30Residues = .Peptide.Substring(0, 27) & "..."
                        End If
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " & dblMassDiffThreshold.ToString("0.00") & " Da away from the mass computed by MODPlus: " & dblPeptideMonoMassPHRP.ToString("0.0000") & " vs. " & dblPeptideMonoMassMODPlus.ToString("0.0000") & "; peptide " & strFirst30Residues)
                    End If

                    If dblPeptideMonoMassMODPlus > 0 Then
                        ' Compute DelM and DelM_PPM
                        dblDelM = dblPrecursorMonoMass - dblPeptideMonoMassMODPlus
                        .DelM = NumToString(dblDelM, 6, True)

                        Dim dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, dblPrecursorMonoMass, True, dblPeptideMonoMassMODPlus)

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
            ' Error parsing this row from the MODPlus results file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not String.IsNullOrEmpty(rowIndex) Then
                    strErrorLog &= "Error parsing MODPlus Results in ParseMODPlusResultsFileEntry for RowIndex '" & rowIndex & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MODPlus Results in ParseMODPlusResultsFileEntry" & ControlChars.NewLine
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
    Private Function ParseMODPlusResultsFileHeaderLine(ByVal strLineIn As String, <Out()> ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        ' The expected column order from MODPlus:
        '   SpectrumFile   Index   ScanNo   ObservedMW   Charge   CalculatedMW   DeltaMass   Score   Probability   Peptide   NTT    Protein   ModificationAnnotation

        Dim lstColumnNames As SortedDictionary(Of String, eMODPlusResultsFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMODPlusResultsFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(MODPlusResultsFileColCount - 1)

        lstColumnNames.Add("SpectrumFile", eMODPlusResultsFileColumns.SpectrumFileName)
        lstColumnNames.Add("Index", eMODPlusResultsFileColumns.SpectrumIndex)
        lstColumnNames.Add("ScanNo", eMODPlusResultsFileColumns.ScanNumber)
        lstColumnNames.Add("ObservedMW", eMODPlusResultsFileColumns.ObservedMonoMass)
        lstColumnNames.Add("Charge", eMODPlusResultsFileColumns.Charge)
        lstColumnNames.Add("CalculatedMW", eMODPlusResultsFileColumns.CalculatedMonoMass)
        lstColumnNames.Add("DeltaMass", eMODPlusResultsFileColumns.DeltaMass)
        lstColumnNames.Add("Score", eMODPlusResultsFileColumns.Score)
        lstColumnNames.Add("Probability", eMODPlusResultsFileColumns.Probability)
        lstColumnNames.Add("Peptide", eMODPlusResultsFileColumns.Peptide)
        lstColumnNames.Add("NTT", eMODPlusResultsFileColumns.NTT)
        lstColumnNames.Add("Protein", eMODPlusResultsFileColumns.ProteinAndPeptidePositionList)
        lstColumnNames.Add("ModificationAnnotation", eMODPlusResultsFileColumns.ModificationAnnotation)

        Try
            ' Initialize each entry in intColumnMapping to -1
            For intIndex As Integer = 0 To intColumnMapping.Length - 1
                intColumnMapping(intIndex) = -1
            Next

            Dim strSplitLine = strLineIn.Split(ControlChars.Tab)
            Dim blnUseDefaultHeaders As Boolean = False

            Dim value As Integer
            If strSplitLine.Length >= 2 Then
                If Integer.TryParse(strSplitLine(1), value) Then
                    ' Second column has a number; this is not a header line					
                    blnUseDefaultHeaders = True
                Else

                    For intIndex As Integer = 0 To strSplitLine.Length - 1
                        Dim eResultFileColumn As eMODPlusResultsFileColumns

                        If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
                            ' Recognized column name; update intColumnMapping
                            intColumnMapping(eResultFileColumn) = intIndex
                            blnUseDefaultHeaders = False
                        Else
                            ' Unrecognized column name
                            Console.WriteLine("Warning: Unrecognized column header name '" & strSplitLine(intIndex) & "' in ParseMODPlusResultsFileHeaderLine")
                        End If
                    Next

                End If

            End If

            If blnUseDefaultHeaders Then
                ' Use default column mappings
                For intIndex As Integer = 0 To intColumnMapping.Length - 1
                    intColumnMapping(intIndex) = intIndex
                Next

                ' This is not a header line; return false
                Return False
            End If

        Catch ex As Exception
            SetErrorMessage("Error parsing header in MODPlus results file: " & ex.Message)
            Return False
        End Try

        ' Header line found and parsed; return true
        Return True

    End Function

    Private Function ParseMODPlusSynFileHeaderLine(ByVal strLineIn As String, <Out()> ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMODPlusSynFileColumns
        Dim lstColumnNames As SortedDictionary(Of String, eMODPlusSynFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMODPlusSynFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(MODPlusSynFileColCount - 1)

        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_ResultID, eMODPlusSynFileColumns.ResultID)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Scan, eMODPlusSynFileColumns.Scan)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Spectrum_Index, eMODPlusSynFileColumns.Spectrum_Index)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Charge, eMODPlusSynFileColumns.Charge)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_PrecursorMZ, eMODPlusSynFileColumns.PrecursorMZ)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_DelM, eMODPlusSynFileColumns.DelM)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_DelM_PPM, eMODPlusSynFileColumns.DelM_PPM)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_MH, eMODPlusSynFileColumns.MH)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Peptide, eMODPlusSynFileColumns.Peptide)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_NTT, eMODPlusSynFileColumns.NTT)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Modification_Annotation, eMODPlusSynFileColumns.ModificationAnnotation)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Protein, eMODPlusSynFileColumns.Protein)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Peptide_Position, eMODPlusSynFileColumns.Peptide_Position)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Score, eMODPlusSynFileColumns.Score)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Probability, eMODPlusSynFileColumns.Probability)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_Rank_Score, eMODPlusSynFileColumns.Rank_Score)
        lstColumnNames.Add(clsPHRPParserMODPlus.DATA_COLUMN_QValue, eMODPlusSynFileColumns.QValue)

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
            SetErrorMessage("Error parsing header in MODPlus synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMODPlusSynFileEntry(
      ByVal strLineIn As String,
      ByVal objSearchResult As clsSearchResultsMODPlus,
      ByRef strErrorLog As String,
      ByVal intResultsProcessed As Integer,
      ByRef intColumnMapping() As Integer,
      <Out()> ByRef strPeptideSequenceWithMods As String) As Boolean

        ' Parses an entry from the MODPlus Synopsis file

        Dim strSplitLine As String() = Nothing

        strPeptideSequenceWithMods = String.Empty

        Try

            ' Reset objSearchResult
            objSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length < 13 Then
                Return False
            End If

            With objSearchResult
                Dim strValue As String = Nothing

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.ResultID), strValue) Then
                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                        strErrorLog &= "Error reading ResultID value from MODPlus Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                    End If
                    Exit Try
                End If

                .ResultID = Integer.Parse(strValue)

                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Scan), .Scan)
                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Charge), .Charge)

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Peptide), strPeptideSequenceWithMods) Then
                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                        strErrorLog &= "Error reading Peptide sequence value from MODPlus Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                    End If
                    Exit Try
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Protein), .ProteinName)
                .MultipleProteinCount = "0"

                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.DelM), .MODPlusComputedDelM)
                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.DelM_PPM), .MODPlusComputedDelMPPM)

                .PeptideDeltaMass = .MODPlusComputedDelM

                ' Note: .PeptideDeltaMass is stored in the MODPlus results file as "Observed_Mass - Theoretical_Mass"
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
                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Spectrum_Index), .Spectrum_Index)

                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.PrecursorMZ), .Precursor_mz)

                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.MH), .ParentIonMH)

                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Score), .MODPlusScore)
                GetColumnValue(strSplitLine, intColumnMapping(eMODPlusSynFileColumns.Probability), .Probability)

            End With

            Return True

        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MODPlus Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MODPlus Results in ParseMODPlusSynFileEntry" & ControlChars.NewLine
                End If
            End If
        End Try

        Return False

    End Function

    ''' <summary>
    ''' Main processing function
    ''' </summary>
    ''' <param name="strInputFilePath">MODPlus results file (Dataset_MODPlus.id.txt)</param>
    ''' <param name="strOutputFolderPath">Output folder</param>
    ''' <param name="strParameterFilePath">Parameter file</param>
    ''' <returns>True if success, False if failure</returns>
    Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean

        Dim strBaseName As String = String.Empty
        Dim strSynOutputFilePath As String = String.Empty

        Dim lstMODPlusModInfo As List(Of clsModificationDefinition)


        Dim blnSuccess As Boolean

        If Not LoadParameterFileSettings(strParameterFilePath) Then
            SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, True)
            Return False
        End If

        Try
            If String.IsNullOrWhiteSpace(strInputFilePath) Then
                SetErrorMessage("Input file name is empty")
                SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath)
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

                lstMODPlusModInfo = New List(Of clsModificationDefinition)

                ' Load the MODPlus Parameter File to look for any static mods
                ExtractModInfoFromMODPlusParamFile(mSearchToolParameterFilePath, lstMODPlusModInfo)

                ' Resolve the mods in lstMODPlusModInfo with the ModDefs mods
                ResolveMODPlusModsWithModDefinitions(lstMODPlusModInfo)

                ' Define the base output filename using strInputFilePath
                strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath)

                ' Auto-replace "modp.id" with "_modp"
                If strBaseName.ToLower().EndsWith("_modp.id") Then
                    strBaseName = strBaseName.Substring(0, strBaseName.Length - "_modp.id".Length) & "_modp"
                End If

                ' Do not create a first-hits file for MODPlus results

                ' Create the synopsis output file
                MyBase.ResetProgress("Creating the SYN file")
                Console.WriteLine()
                Console.WriteLine()
                Console.WriteLine(MyBase.ProgressStepDescription)

                strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath, True)

                ' Create the other PHRP-specific files
                MyBase.ResetProgress("Creating the PHRP files for " & Path.GetFileName(strSynOutputFilePath))
                Console.WriteLine()
                Console.WriteLine()
                Console.WriteLine(MyBase.ProgressStepDescription)

                ' Now parse the _syn.txt file that we just created to next create the other PHRP files
                blnSuccess = ParseMODPlusSynopsisFile(strSynOutputFilePath, strOutputFolderPath, False)

                If blnSuccess AndAlso mCreateProteinModsFile Then
                    blnSuccess = CreateProteinModsFileWork(strBaseName, fiInputFile, strSynOutputFilePath, strOutputFolderPath)
                End If

                If blnSuccess Then
                    MyBase.OperationComplete()
                End If

            Catch ex As Exception
                SetErrorMessage("Error in clsMODPlusResultsProcessor.ProcessFile (2):  " & ex.Message)
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
            End Try

        Catch ex As Exception
            SetErrorMessage("Error in ProcessFile (1):" & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.UnspecifiedError)
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
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        Else
            If File.Exists(strMTSPepToProteinMapFilePath) AndAlso mUseExistingMTSPepToProteinMapFile Then
                blnSuccess = True
            Else
                ' Auto-change mIgnorePeptideToProteinMapperErrors to True
                ' We only do this because some peptides reported by MODPlus may not match the fasta file (due to amino acid substitutions)
                mIgnorePeptideToProteinMapperErrors = True
                blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
                If Not blnSuccess Then
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False")
                End If
            End If
        End If

        If blnSuccess Then
            ' If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
            MyBase.ValidatePHRPReaderSupportFiles(IO.Path.Combine(fiInputFile.DirectoryName, IO.Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath)

            ' Create the Protein Mods file
            blnSuccess = MyBase.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MODPlus)
        End If

        If Not blnSuccess Then
            ' Do not treat this as a fatal error
            blnSuccess = True
        End If
        Return blnSuccess
    End Function

    Protected Sub ResolveMODPlusModsWithModDefinitions(ByRef lstMODPlusModInfo As List(Of clsModificationDefinition))

        Dim blnExistingModFound As Boolean
        Dim objModDef As clsModificationDefinition

        If Not lstMODPlusModInfo Is Nothing Then

            ' Call .LookupModificationDefinitionByMass for each entry in lstMODPlusModInfo
            For Each objModInfo As clsModificationDefinition In lstMODPlusModInfo
                If String.IsNullOrEmpty(objModInfo.TargetResidues) Then
                    objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, Nothing, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnExistingModFound, True, MODPlus_MASS_DIGITS_OF_PRECISION)
                Else
                    For Each chTargetResidue As Char In objModInfo.TargetResidues
                        objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, chTargetResidue, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnExistingModFound, True, MODPlus_MASS_DIGITS_OF_PRECISION)
                    Next
                End If
            Next

        End If

    End Sub

    Private Sub SortAndWriteFilteredSearchResults(
      ByVal swResultFile As StreamWriter,
      ByVal lstFilteredSearchResults As List(Of udtMODPlusSearchResultType),
      ByRef strErrorLog As String)

        ' Sort udtFilteredSearchResults by descending score, ascending scan, ascending charge, ascending peptide, and ascending protein
        lstFilteredSearchResults.Sort(New MODPlusSearchResultsComparerScoreScanChargePeptide)

        ' Compute FDR values then assign QValues
        ComputeQValues(lstFilteredSearchResults)

        If mProteinNamePositionSplit Is Nothing Then
            mProteinNamePositionSplit = New Regex("(.+)\[([^\]]+)\]", RegexOptions.Compiled)
        End If

        Dim resultID = 1
        For Each result In lstFilteredSearchResults

            Dim proteinList = result.ProteinList.Split(";"c)

            If proteinList.Length = 0 Then
                ' This code should not be reached
                WriteSearchResultToFile(resultID, swResultFile, result, "Unknown_Protein", String.Empty, strErrorLog)
                resultID += 1
            End If

            For Each proteinEntry In proteinList
                Dim proteinName As String
                Dim peptidePosition As String

                Dim reMatch = mProteinNamePositionSplit.Match(proteinEntry)
                If reMatch.Success Then
                    proteinName = reMatch.Groups(1).Value
                    peptidePosition = reMatch.Groups(2).Value
                Else
                    proteinName = String.Copy(proteinEntry)
                    peptidePosition = String.Empty
                End If

                WriteSearchResultToFile(resultID, swResultFile, result, proteinName, peptidePosition, strErrorLog)
                resultID += 1

            Next

        Next

    End Sub

    ''' <summary>
    ''' Compute FDR values then assign QValues
    ''' </summary>
    ''' <param name="lstSearchResults"></param>
    ''' <remarks>Assumes the data is sorted by descending score using MODPlusSearchResultsComparerScoreScanChargePeptide</remarks>
    Private Sub ComputeQValues(ByVal lstSearchResults As List(Of udtMODPlusSearchResultType))

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
                Dim proteinList = lstSearchResults(intIndexCheck).ProteinList.Split(";"c)

                For Each proteinEntry In proteinList
                    If Not IsReversedProtein(proteinEntry) Then
                        isReverse = False
                        Exit For
                    End If
                Next                
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
     ByVal lstSearchResults As List(Of udtMODPlusSearchResultType),
     ByVal intStartIndex As Integer,
     ByVal intEndIndex As Integer,
     ByVal lstFilteredSearchResults As List(Of udtMODPlusSearchResultType))

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
      ByVal swResultFile As StreamWriter,
      ByRef strErrorLog As String)

        ' Write out the header line for synopsis / first hits files
        Try
            Dim lstData As New List(Of String)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_ResultID)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Scan)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Spectrum_Index)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Charge)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_PrecursorMZ)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_DelM)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_DelM_PPM)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_MH)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Peptide)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_NTT)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Modification_Annotation)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Protein)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Peptide_Position)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Score)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Probability)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_Rank_Score)
            lstData.Add(clsPHRPParserMODPlus.DATA_COLUMN_QValue)

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
      ByVal intResultID As Integer,
      ByVal swResultFile As StreamWriter,
      ByRef udtSearchResult As udtMODPlusSearchResultType,
      ByVal proteinName As String,
      ByVal peptidePosition As String,
      ByRef strErrorLog As String)

        Try

            ' Primary Columns
            '
            ' MODPlus
            ' ResultID	Scan	Spectrum_Index	Charge	PrecursorMZ	DelM	DelM_PPM	MH	Peptide	NTT	ModificationAnnotation	Protein	Peptide_Position	Score	Probability	Rank_Probability   QValue

            Dim lstData As New List(Of String)
            lstData.Add(intResultID.ToString)
            lstData.Add(udtSearchResult.ScanNum.ToString)
            lstData.Add(udtSearchResult.SpectrumIndex)
            lstData.Add(udtSearchResult.Charge)
            lstData.Add(udtSearchResult.PrecursorMZ)
            lstData.Add(udtSearchResult.DelM)
            lstData.Add(udtSearchResult.DelM_PPM)
            lstData.Add(udtSearchResult.MH)
            lstData.Add(udtSearchResult.Peptide)
            lstData.Add(udtSearchResult.NTT)
            lstData.Add(udtSearchResult.ModificationAnnotation)
            lstData.Add(proteinName)
            lstData.Add(peptidePosition)
            lstData.Add(udtSearchResult.Score)
            lstData.Add(udtSearchResult.Probability)
            lstData.Add(udtSearchResult.RankScore.ToString())
            lstData.Add(udtSearchResult.QValue.ToString("0.000000"))

            swResultFile.WriteLine(CollapseList(lstData))

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

#Region "IComparer Classes"

    Protected Class MODPlusSearchResultsComparerScanChargeScorePeptide
        Implements IComparer(Of udtMODPlusSearchResultType)

        Public Function Compare(x As udtMODPlusSearchResultType, y As udtMODPlusSearchResultType) As Integer Implements IComparer(Of udtMODPlusSearchResultType).Compare

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
                    ' Charge is the same; check ScoreNum
                    If x.ScoreNum < y.ScoreNum Then
                        Return 1
                    ElseIf x.ScoreNum > y.ScoreNum Then
                        Return -1
                    Else
                        ' Probability is the same; check peptide
                        If x.Peptide > y.Peptide Then
                            Return 1
                        ElseIf x.Peptide < y.Peptide Then
                            Return -1
                        Else
                            ' Peptide is the same, check Protein
                            If x.ProteinList > y.ProteinList Then
                                Return 1
                            ElseIf x.ProteinList < y.ProteinList Then
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

    Protected Class MODPlusSearchResultsComparerScoreScanChargePeptide
        Implements IComparer(Of udtMODPlusSearchResultType)

        Public Function Compare(x As udtMODPlusSearchResultType, y As udtMODPlusSearchResultType) As Integer Implements IComparer(Of udtMODPlusSearchResultType).Compare

            If x.ScoreNum < y.ScoreNum Then
                Return 1
            ElseIf x.ScoreNum > y.ScoreNum Then
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
                            If x.ProteinList > y.ProteinList Then
                                Return 1
                            ElseIf x.ProteinList < y.ProteinList Then
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
