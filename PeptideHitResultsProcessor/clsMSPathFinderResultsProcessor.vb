Option Strict On

' This class reads in a MSPathFinder results file (_IcTda.tsv) and creates 
' a tab-delimited text file with the data. 
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 05/01/2015
'
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports PHRPReader
Imports System.IO
Imports System.Runtime.InteropServices
Imports System.Text.RegularExpressions

Public Class clsMSPathFinderResultsProcessor
    Inherits clsPHRPBaseClass

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "June 22, 2015"

        mGetModName = New Regex("(.+) (\d+)", RegexOptions.Compiled)
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MSPathFinder_FILE As String = "_IcTda"

    Public Const N_TERMINUS_SYMBOL_MSPATHFINDER As String = "-"
    Public Const C_TERMINUS_SYMBOL_MSPATHFINDER As String = "-"

    Public Const DEFAULT_SYN_FILE_QVALUE_THRESHOLD As Single = 0.05

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Private Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline Or RegexOptions.IgnoreCase

    ' These columns correspond to the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
    Protected Const MSPathFinderResultsFileColCount As Integer = 17
    Public Enum eMSPathFinderResultsFileColumns As Integer
        Scan = 0
        PrefixResidue = 1
        Sequence = 2
        SuffixResidue = 3
        Modifications = 4
        Composition = 5
        Protein = 6
        ProteinDesc = 7
        ProteinLength = 8
        ResidueStart = 9
        ResidueEnd = 10
        Charge = 11
        MostAbundantIsotopeMz = 12
        CalculatedMonoMass = 13
        NumMatchedFragments = 14
        QValue = 15
        PepQValue = 16
        ' Future: SpecProb = 17
    End Enum

    ' These columns correspond to the Synopsis file created by this class
    Protected Const MSPathFinderSynFileColCount As Integer = 16
    Public Enum eMSPathFinderSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        Charge = 2
        MostAbundantIsotopeMz = 3
        Mass = 4
        Sequence = 5                ' PrefixLetter.Sequence.SuffixLetter
        Modifications = 6
        Composition = 7
        Protein = 8
        ProteinDesc = 9
        ProteinLength = 10
        ResidueStart = 11
        ResidueEnd = 12
        MatchedFragments = 13
        QValue = 14
        PepQValue = 15
        ' Future: SpecProb = 16
    End Enum

#End Region

#Region "Structures"
    ' This data structure holds rows read from the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
    Protected Structure udtMSPathFinderSearchResultType

        Public Scan As String
        Public ScanNum As Integer
        Public PrefixResidue As String
        Public Sequence As String
        Public SuffixResidue As String
        Public Modifications As String
        Public Composition As String
        Public Protein As String
        Public ProteinDesc As String
        Public ProteinLength As String
        Public ResidueStart As String
        Public ResidueEnd As String
        Public Charge As String
        Public ChargeNum As Short
        Public MostAbundantIsotopeMz As String      ' As reported by MSPathfinder
        Public CalculatedMonoMass As String         ' Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MSPathFinder
        Public CalculatedMonoMassPHRP As Double     ' Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by PHRP
        Public NumMatchedFragments As String
        Public QValue As String                     ' FDR, at the scan level
        Public QValueNum As Double
        Public PepQValue As String                  ' FDR, at the peptide level

        ' Future: Public SpecProb As String

        ' The following are typically defined for other search engines, but are not used for MSPathFinder
        '   Public DelM As String                   ' Computed using Precursor_mass - CalculatedMonoMass
        '   Public DelM_PPM As String               ' Computed using DelM and CalculatedMonoMass	

        Public Sub Clear()
            Scan = String.Empty
            ScanNum = 0
            PrefixResidue = String.Empty
            Sequence = String.Empty
            SuffixResidue = String.Empty
            Modifications = String.Empty
            Composition = String.Empty
            Protein = String.Empty
            ProteinDesc = String.Empty
            ProteinLength = String.Empty
            ResidueStart = String.Empty
            ResidueEnd = String.Empty
            Charge = String.Empty
            ChargeNum = 0
            MostAbundantIsotopeMz = String.Empty
            CalculatedMonoMass = String.Empty
            CalculatedMonoMassPHRP = 0
            NumMatchedFragments = String.Empty
            QValue = String.Empty
            QValueNum = 0
            PepQValue = String.Empty
            ' Future: SpecProb = String.empty

            ' Unused at present: MH = String.Empty
            ' Unused at present: DelM = String.Empty
            ' Unused at present: DelM_PPM = String.Empty

        End Sub
    End Structure

#End Region

#Region "Classwide Variables"

    Private ReadOnly mGetModName As Regex
#End Region


    ''' <summary>
    ''' Step through the Modifications and associate each modification with the residues
    ''' For each residue, check if a static mod is defined that affects that residue
    ''' For each mod mass, determine the modification and add to objSearchResult
    ''' </summary>
    ''' <param name="objSearchResult"></param>
    ''' <param name="blnUpdateModOccurrenceCounts"></param>
    ''' <remarks></remarks>
    Private Sub AddModificationsToResidues(
       objSearchResult As clsSearchResultsMSPathFinder,
       blnUpdateModOccurrenceCounts As Boolean,
       lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType))

        If String.IsNullOrWhiteSpace(objSearchResult.Modifications) Then
            Return
        End If

        Dim lstMods = objSearchResult.Modifications.Split(","c)
        Dim finalResidueLoc = objSearchResult.PeptideCleanSequence.Length

        For Each modEntry In lstMods

            ' Find the mod in the list of known modifications (loaded from the MSPathFinder parameter file)
            Dim matchFound = False

            ' Obtain the mod name, for example "Dehydro" from "Dehydro 52"
            Dim reMatch = mGetModName.Match(modEntry)

            If Not reMatch.Success Then
                ReportError("Mod entry does not have a name separated by a number: " & modEntry, False)
                Continue For
            End If

            Dim modName = reMatch.Groups(1).Value
            Dim residueNumber = reMatch.Groups(2).Value

            For Each modDef As clsMSGFPlusParamFileModExtractor.udtModInfoType In lstModInfo
                If String.Equals(modDef.ModName, modName, StringComparison.CurrentCultureIgnoreCase) Then

                    Dim intResidueLocInPeptide As Integer
                    If Not Integer.TryParse(residueNumber, intResidueLocInPeptide) Then
                        ReportError("Mod entry does not have a number after the name: " & modEntry, False)
                        Continue For
                    End If

                    Dim eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
                    If intResidueLocInPeptide <= 1 Then
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
                    ElseIf intResidueLocInPeptide >= finalResidueLoc Then
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
                    End If

                    ' Now that we know the terminus position, assure that intResidueLocInPeptide is 1 not 0
                    If intResidueLocInPeptide < 1 Then
                        intResidueLocInPeptide = 1
                    ElseIf intResidueLocInPeptide > finalResidueLoc Then
                        intResidueLocInPeptide = finalResidueLoc
                    End If

                    Dim chMostRecentResidue = "-"c

                    If intResidueLocInPeptide >= 1 AndAlso intResidueLocInPeptide <= finalResidueLoc Then
                        chMostRecentResidue = objSearchResult.PeptideCleanSequence.Chars(intResidueLocInPeptide - 1)
                    End If

                    ' Associate the mod with the given residue
                    objSearchResult.SearchResultAddModification(modDef.ModMassVal, chMostRecentResidue, intResidueLocInPeptide, eResidueTerminusState, blnUpdateModOccurrenceCounts)

                    matchFound = True
                    Exit For

                End If
            Next

            If Not matchFound Then
                ReportError("Mod name " & modName & " was not defined in the MSPathFinder parameter file; cannot determine mod mass", False)
            End If
        Next


    End Sub

    Private Function AddModificationsAndComputeMass(
       objSearchResult As clsSearchResultsMSPathFinder,
       blnUpdateModOccurrenceCounts As Boolean,
       lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Boolean

        Dim blnSuccess As Boolean

        Try

            ' For other tools, we would add IsotopicMods here
            ' This is not supported for MSPathFinder
            ' 
            ' objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

            ' Parse .Modifications to determine the modified residues present
            AddModificationsToResidues(objSearchResult, blnUpdateModOccurrenceCounts, lstModInfo)

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

    Protected Function ComputePeptideMass(strCleanSequence As String, dblTotalModMass As Double) As Double

        Dim dblMass As Double

        dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence)
        dblMass += dblTotalModMass

        Return dblMass

    End Function

    ''' <summary>
    ''' Computes the total of all modifications defined for the sequence
    ''' </summary>
    ''' <param name="modificationList">Comma separated list of modifications, e.g. Dehydro 52,Dehydro 63</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ComputeTotalModMass(
      cleanSequence As String,
      modificationList As String,
      lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Double

        If String.IsNullOrWhiteSpace(modificationList) Then
            Return 0
        End If

        Dim dblTotalModMass As Double = 0
        Dim lstMods = modificationList.Split(","c)

        For Each modEntry In lstMods

            ' Convert the mod name to a mass value
            Dim matchFound = False

            ' Obtain the mod name, for example "Dehydro" from "Dehydro 52"
            Dim reMatch = mGetModName.Match(modEntry)

            If Not reMatch.Success Then
                ReportError("Mod entry does not have a name separated by a number: " & modEntry, True)
            End If

            Dim modName = reMatch.Groups(1).Value

            For Each modDef As clsMSGFPlusParamFileModExtractor.udtModInfoType In lstModInfo
                If String.Equals(modDef.ModName, modName, StringComparison.CurrentCultureIgnoreCase) Then
                    dblTotalModMass += modDef.ModMassVal
                    matchFound = True
                    Exit For
                End If
            Next

            If Not matchFound Then
                ReportError("Mod name " & modName & " was not defined in the MSPathFinder parameter file; cannot determine mod mass", True)
            End If
        Next

        Return dblTotalModMass

    End Function

    Private Function CreateProteinModsFileWork(
     strBaseName As String,
     fiInputFile As FileInfo,
     strSynOutputFilePath As String,
     strOutputFolderPath As String) As Boolean

        Dim blnSuccess As Boolean

        ' Create the MTSPepToProteinMap file

        Dim strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS:=True)

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
                blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
                If Not blnSuccess Then
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False")
                End If
            End If
        End If

        If blnSuccess Then
            ' If necessary, copy various PHRPReader support files to the output folder
            MyBase.ValidatePHRPReaderSupportFiles(IO.Path.Combine(fiInputFile.DirectoryName, IO.Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath)

            ' Create the Protein Mods file
            blnSuccess = MyBase.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSPathFinder)
        End If

        If Not blnSuccess Then
            ' Do not treat this as a fatal error
            blnSuccess = True
        End If

        Return blnSuccess

    End Function

    ''' <summary>
    ''' This routine creates a first hits file or synopsis file from the output from MSPathFinder
    ''' The synopsis file includes every result with a probability above a set threshold
    ''' The first-hits file includes the result with the highest probability (for each scan and charge)
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CreateSynResultsFile(
      strInputFilePath As String,
      strOutputFilePath As String,
      lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Boolean

        Try
            Dim intColumnMapping() As Integer = Nothing
            Dim strErrorLog = String.Empty

            ' Open the input file and parse it
            ' Initialize the stream reader and the stream Text writer
            Using srDataFile = New StreamReader(New FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)),
                  swResultFile = New StreamWriter(New FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                Dim headerParsed = False
                Dim rowNumber = 0

                ' Initialize the list that will hold all of the records in the MSPathFinder result file
                Dim lstSearchResultsUnfiltered = New List(Of udtMSPathFinderSearchResultType)

                ' Initialize the list that will hold all of the records that will ultimately be written out to disk
                Dim lstFilteredSearchResults = New List(Of udtMSPathFinderSearchResultType)

                ' Parse the input file
                Do While Not srDataFile.EndOfStream And Not MyBase.AbortProcessing
                    Dim strLineIn = srDataFile.ReadLine()
                    rowNumber += 1

                    If String.IsNullOrWhiteSpace(strLineIn) Then
                        Continue Do
                    End If

                    If Not headerParsed Then
                        ' Parse the header line

                        Dim blnSuccess = ParseMSPathFinderResultsFileHeaderLine(strLineIn, intColumnMapping)
                        If Not blnSuccess Then
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

                    Dim udtSearchResult = New udtMSPathFinderSearchResultType

                    Dim blnValidSearchResult = ParseMSPathFinderResultsFileEntry(strLineIn, udtSearchResult, strErrorLog, intColumnMapping, lstModInfo, rowNumber)

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

                ' Sort the SearchResults by scan, charge, and ascending QValue (future, ascending SpecProb)
                lstSearchResultsUnfiltered.Sort(New MSPathFinderSearchResultsComparerScanChargeScorePeptide)

                ' Now filter the data

                ' Initialize variables
                Dim intStartIndex = 0
                Dim intEndIndex As Integer

                intStartIndex = 0
                Do While intStartIndex < lstSearchResultsUnfiltered.Count
                    intEndIndex = intStartIndex
                    Do While intEndIndex + 1 < lstSearchResultsUnfiltered.Count AndAlso
                             lstSearchResultsUnfiltered(intEndIndex + 1).ScanNum = lstSearchResultsUnfiltered(intStartIndex).ScanNum
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

    Private Function ExtractModInfoFromParamFile(
       strMSGFDBParamFilePath As String,
       <Out()> ByRef lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Boolean

        ' The DMS-based parameter file for MSPathFinder uses the same formatting as MSGF+

        Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor("MSPathFinder")

        AddHandler modFileProcessor.ErrorOccurred, AddressOf ModExtractorErrorHandler
        AddHandler modFileProcessor.WarningMessageEvent, AddressOf ModExtractorWarningHandler

        ' Note that this call will intialize lstModInfo
        Dim success = modFileProcessor.ExtractModInfoFromParamFile(strMSGFDBParamFilePath, lstModInfo)

        If Not success OrElse mErrorCode <> ePHRPErrorCodes.NoError Then
            If mErrorCode = ePHRPErrorCodes.NoError Then
                SetErrorMessage("Unknown error extracting the modification definitions from the MSPathFinder parameter file")
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            End If
            Return False
        End If

        modFileProcessor.ResolveMSGFDBModsWithModDefinitions(lstModInfo, mPeptideMods)

        Return True

    End Function

    ''' <summary>
    ''' Parse the Synopsis file to create the other PHRP-compatible files
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFolderPath"></param>
    ''' <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ParseMSPathfinderSynopsisFile(
      strInputFilePath As String,
      strOutputFolderPath As String,
      blnResetMassCorrectionTagsAndModificationDefinitions As Boolean,
      lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType)) As Boolean

        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim strPreviousQValue As String

        ' Note that ParseMSPathfinderSynopsisFile synopsis files are normally sorted on Probability value, ascending
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
            Dim objSearchResult = New clsSearchResultsMSPathFinder(mPeptideMods)

            ' Initialize htPeptidesFoundForQValue
            Dim htPeptidesFoundForQValue = New Hashtable
            Dim blnFirstMatchForGroup As Boolean

            strPreviousQValue = String.Empty

            Dim strErrorLog = String.Empty

            Try
                objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(mEnzymeMatchSpec, mPeptideNTerminusMassChange, mPeptideCTerminusMassChange)

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
                            blnSuccess = ParseMSPathFinderSynFileHeaderLine(strLineIn, intColumnMapping)
                            If Not blnSuccess Then
                                ' Error parsing header
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                                Exit Try
                            End If
                            blnHeaderParsed = True
                            Continue Do
                        End If

                        Dim strCurrentPeptideWithMods As String = String.Empty

                        Dim blnValidSearchResult = ParseMSPathFinderSynFileEntry(
                          strLineIn, objSearchResult, strErrorLog,
                          intResultsProcessed, intColumnMapping,
                          strCurrentPeptideWithMods)

                        If Not blnValidSearchResult Then
                            Continue Do
                        End If

                        Dim strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge

                        If objSearchResult.QValue = strPreviousQValue Then
                            ' New result has the same QValue as the previous result
                            ' See if htPeptidesFoundForQValue contains the peptide, scan and charge

                            If htPeptidesFoundForQValue.ContainsKey(strKey) Then
                                blnFirstMatchForGroup = False
                            Else
                                htPeptidesFoundForQValue.Add(strKey, 1)
                                blnFirstMatchForGroup = True
                            End If

                        Else
                            ' New QValue
                            ' Reset htPeptidesFoundForQValue
                            htPeptidesFoundForQValue.Clear()

                            ' Update strPreviousQValue
                            strPreviousQValue = objSearchResult.QValue

                            ' Append a new entry to htPeptidesFoundForQValue
                            htPeptidesFoundForQValue.Add(strKey, 1)
                            blnFirstMatchForGroup = True
                        End If

                        blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup, lstModInfo)
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
    ''' Parse a MSPathFinder results line while creating the MSPathFinder synopsis file
    ''' </summary>
    ''' <param name="strLineIn"></param>
    ''' <param name="udtSearchResult"></param>
    ''' <param name="strErrorLog"></param>
    ''' <param name="intColumnMapping"></param>
    ''' <param name="lstModInfo"></param>
    ''' <param name="rowNumber">Row number (used for error reporting)</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function ParseMSPathFinderResultsFileEntry(
      strLineIn As String,
      ByRef udtSearchResult As udtMSPathFinderSearchResultType,
      ByRef strErrorLog As String,
      intColumnMapping() As Integer,
      lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType),
      rowNumber As Integer) As Boolean

        ' Parses an entry from the MSPathFinder results file

        Dim dblSequenceMonoMassMSPathFinder As Double    ' Theoretical peptide monoisotopic mass, including mods, as computed by MSPathFinder

        Dim blnValidSearchResult As Boolean

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            udtSearchResult.Clear()
            Dim strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length >= 11 Then

                With udtSearchResult

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.Scan), .Scan) Then
                        ReportError("Scan column is missing or invalid in row " & rowNumber, True)
                    End If

                    If Not Integer.TryParse(.Scan, .ScanNum) Then
                        ReportError("Scan column is not numeric in row " & rowNumber, True)
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.Charge), .Charge)
                    .ChargeNum = CShort(CIntSafe(.Charge, 0))

                    ' Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MSPathFinder
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.CalculatedMonoMass), .CalculatedMonoMass)
                    Double.TryParse(.CalculatedMonoMass, dblSequenceMonoMassMSPathFinder)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.PrefixResidue), .PrefixResidue)

                    If Not GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.Sequence), .Sequence) Then
                        ReportError("Sequence column is missing or invalid in row " & rowNumber, True)
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.SuffixResidue), .SuffixResidue)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.Modifications), .Modifications)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.Composition), .Composition)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.Protein), .Protein)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.ProteinDesc), .ProteinDesc)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.ProteinLength), .ProteinLength)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.ResidueEnd), .ResidueEnd)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.ResidueStart), .ResidueStart)

                    ' Parse the list of modified residues to determine the total mod mass
                    Dim dblTotalModMass = ComputeTotalModMass(.Sequence, .Modifications, lstModInfo)

                    ' Compute monoisotopic mass of the peptide
                    .CalculatedMonoMassPHRP = ComputePeptideMass(.Sequence, dblTotalModMass)

                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.MostAbundantIsotopeMz), .MostAbundantIsotopeMz)
                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.NumMatchedFragments), .NumMatchedFragments)

                    If GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.QValue), .QValue) Then
                        Double.TryParse(.QValue, .QValueNum)
                    End If

                    GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.PepQValue), .PepQValue)

                    ' Future:
                    'If GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderResultsFileColumns.SpecProb), .SpecProb) Then
                    '    Double.TryParse(.SpecProb, .SpecProbNum)
                    'End If

                End With

                blnValidSearchResult = True
            End If

        Catch ex As Exception
            ' Error parsing this row from the MassMSPathFinder results file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error parsing MassMSPathFinder Results in ParseMSPathFinderResultsFileEntry for Row " & rowNumber & ControlChars.NewLine
            End If

            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    ''' <summary>
    ''' Parse the MSPathFinder results file header line
    ''' </summary>
    ''' <param name="strLineIn"></param>
    ''' <param name="intColumnMapping"></param>
    ''' <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
    ''' <remarks></remarks>
    Private Function ParseMSPathFinderResultsFileHeaderLine(strLineIn As String, <Out()> ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        ' The expected column order from MassMSPathFinder:
        '   Scan	Pre	Sequence	Post	Modifications	Composition	ProteinName	ProteinDesc	ProteinLength	Start	End	Charge	MostAbundantIsotopeMz	Mass	#MatchedFragments	QValue	PepQValue


        Dim lstColumnNames As SortedDictionary(Of String, eMSPathFinderResultsFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMSPathFinderResultsFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(MSPathFinderResultsFileColCount - 1)

        lstColumnNames.Add("Scan", eMSPathFinderResultsFileColumns.Scan)
        lstColumnNames.Add("Pre", eMSPathFinderResultsFileColumns.PrefixResidue)
        lstColumnNames.Add("Sequence", eMSPathFinderResultsFileColumns.Sequence)
        lstColumnNames.Add("Post", eMSPathFinderResultsFileColumns.SuffixResidue)
        lstColumnNames.Add("Modifications", eMSPathFinderResultsFileColumns.Modifications)
        lstColumnNames.Add("Composition", eMSPathFinderResultsFileColumns.Composition)
        lstColumnNames.Add("ProteinName", eMSPathFinderResultsFileColumns.Protein)
        lstColumnNames.Add("ProteinDesc", eMSPathFinderResultsFileColumns.ProteinDesc)
        lstColumnNames.Add("ProteinLength", eMSPathFinderResultsFileColumns.ProteinLength)
        lstColumnNames.Add("Start", eMSPathFinderResultsFileColumns.ResidueStart)
        lstColumnNames.Add("End", eMSPathFinderResultsFileColumns.ResidueEnd)
        lstColumnNames.Add("Charge", eMSPathFinderResultsFileColumns.Charge)
        lstColumnNames.Add("MostAbundantIsotopeMz", eMSPathFinderResultsFileColumns.MostAbundantIsotopeMz)
        lstColumnNames.Add("Mass", eMSPathFinderResultsFileColumns.CalculatedMonoMass)
        lstColumnNames.Add("#MatchedFragments", eMSPathFinderResultsFileColumns.NumMatchedFragments)
        lstColumnNames.Add("QValue", eMSPathFinderResultsFileColumns.QValue)
        lstColumnNames.Add("PepQValue", eMSPathFinderResultsFileColumns.PepQValue)

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
                        Dim eResultFileColumn As eMSPathFinderResultsFileColumns

                        If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
                            ' Recognized column name; update intColumnMapping
                            intColumnMapping(eResultFileColumn) = intIndex
                            blnUseDefaultHeaders = False
                        Else
                            ' Unrecognized column name
                            Console.WriteLine("Warning: Unrecognized column header name '" & strSplitLine(intIndex) & "' in ParseMSPathFinderResultsFileHeaderLine")
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
            SetErrorMessage("Error parsing header in MSPathFinder results file: " & ex.Message)
            Return False
        End Try

        ' Header line found and parsed; return true
        Return True

    End Function

    Private Function ParseMSPathFinderSynFileHeaderLine(strLineIn As String, <Out()> ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eMSPathFinderSynFileColumns
        Dim lstColumnNames As SortedDictionary(Of String, eMSPathFinderSynFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eMSPathFinderSynFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(MSPathFinderSynFileColCount - 1)

        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ResultID, eMSPathFinderSynFileColumns.ResultID)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Scan, eMSPathFinderSynFileColumns.Scan)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Charge, eMSPathFinderSynFileColumns.Charge)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_MostAbundantIsotopeMz, eMSPathFinderSynFileColumns.MostAbundantIsotopeMz)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Mass, eMSPathFinderSynFileColumns.Mass)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Sequence, eMSPathFinderSynFileColumns.Sequence)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Modifications, eMSPathFinderSynFileColumns.Modifications)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Composition, eMSPathFinderSynFileColumns.Composition)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Protein, eMSPathFinderSynFileColumns.Protein)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinDesc, eMSPathFinderSynFileColumns.ProteinDesc)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinLength, eMSPathFinderSynFileColumns.ProteinLength)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueStart, eMSPathFinderSynFileColumns.ResidueStart)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueEnd, eMSPathFinderSynFileColumns.ResidueEnd)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_MatchedFragments, eMSPathFinderSynFileColumns.MatchedFragments)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_QValue, eMSPathFinderSynFileColumns.QValue)
        lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_PepQValue, eMSPathFinderSynFileColumns.PepQValue)
        ' Future: lstColumnNames.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_SpecProb, eMSPathFinderSynFileColumns.SpecProb)

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
            SetErrorMessage("Error parsing header in MSPathFinder synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Function ParseMSPathFinderSynFileEntry(
      strLineIn As String,
      objSearchResult As clsSearchResultsMSPathFinder,
      ByRef strErrorLog As String,
      intResultsProcessed As Integer,
      ByRef intColumnMapping() As Integer,
      <Out()> ByRef strPeptideSequence As String) As Boolean

        ' Parses an entry from the MSPathFinder Synopsis file

        Dim strSplitLine As String() = Nothing

        strPeptideSequence = String.Empty

        Try

            ' Reset objSearchResult
            objSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length < 15 Then
                Return False
            End If

            With objSearchResult
                Dim strValue As String = Nothing

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.ResultID), strValue) Then
                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                        strErrorLog &= "Error reading ResultID value from MSPathFinder Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                    End If
                    Exit Try
                End If

                .ResultID = Integer.Parse(strValue)

                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.Scan), .Scan)
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.Charge), .Charge)

                If Not GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.Sequence), strPeptideSequence) Then
                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                        strErrorLog &= "Error reading Peptide sequence value from MSPathFinder Results line " & (intResultsProcessed + 1).ToString() & ControlChars.NewLine
                    End If
                    Exit Try
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.Protein), .ProteinName)
                .MultipleProteinCount = "0"

                .PeptideDeltaMass = "0"

                ' Note that MSPathFinder sequences don't actually have mod symbols; that informtion is tracked via strModifications

                ' Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                .SetPeptideSequenceWithMods(strPeptideSequence, True, True)

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
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.MostAbundantIsotopeMz), .MostAbundantIsotopeMz)

                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.Modifications), .Modifications)
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.Composition), .Composition)

                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.ProteinDesc), .ProteinDesc)
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.ProteinLength), .ProteinLength)

                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.ResidueStart), .ResidueStart)
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.ResidueEnd), .ResidueEnd)
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.MatchedFragments), .MatchedFragments)

                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.QValue), .QValue)
                GetColumnValue(strSplitLine, intColumnMapping(eMSPathFinderSynFileColumns.PepQValue), .PepQValue)

                ' Update the peptide location in the protein
                If Not String.IsNullOrEmpty(.ResidueStart) Then
                    Integer.TryParse(.ResidueStart, .PeptideLocInProteinStart)
                End If

                If Not String.IsNullOrEmpty(.ResidueEnd) Then
                    Integer.TryParse(.ResidueEnd, .PeptideLocInProteinEnd)
                End If

            End With

            Return True

        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing MSPathFinder Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing MSPathFinder Results in ParseMSPathFinderSynFileEntry" & ControlChars.NewLine
                End If
            End If
        End Try

        Return False

    End Function

    ''' <summary>
    ''' Main processing function
    ''' </summary>
    ''' <param name="strInputFilePath">MSPathFinder results file (Dataset_IcTda.tsv)</param>
    ''' <param name="strOutputFolderPath">Output folder</param>
    ''' <param name="strParameterFilePath">Parameter file</param>
    ''' <returns>True if success, False if failure</returns>
    Public Overrides Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean

        Dim strBaseName As String = String.Empty
        Dim strSynOutputFilePath As String = String.Empty

        Dim lstModInfo As List(Of clsMSGFPlusParamFileModExtractor.udtModInfoType) = Nothing

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

                ' Load the MSPathFinder Parameter File so that we can determine the modification names and masses
                ' Note that this call will intialize lstModInfo
                Dim success = ExtractModInfoFromParamFile(mSearchToolParameterFilePath, lstModInfo)
                If Not success Then
                    Return False
                End If

                ' Define the base output filename using strInputFilePath
                strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath)

                ' Auto-replace "_IcTda.tsv" with "_mspath"
                If strBaseName.ToLower().EndsWith("_IcTda".ToLower()) Then
                    strBaseName = strBaseName.Substring(0, strBaseName.Length - "_IcTda".Length) & "_mspath"
                End If

                ' Do not create a first-hits file for MSPathFinder results

                ' Create the synopsis output file
                MyBase.ResetProgress("Creating the SYN file")
                Console.WriteLine()
                Console.WriteLine()
                Console.WriteLine(MyBase.ProgressStepDescription)

                ' The synopsis file name will be of the form BasePath_mspath_syn.txt
                strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath, lstModInfo)

                ' Create the other PHRP-specific files
                MyBase.ResetProgress("Creating the PHRP files for " & Path.GetFileName(strSynOutputFilePath))
                Console.WriteLine()
                Console.WriteLine()
                Console.WriteLine(MyBase.ProgressStepDescription)

                ' Now parse the _syn.txt file that we just created to next create the other PHRP files
                blnSuccess = ParseMSPathfinderSynopsisFile(strSynOutputFilePath, strOutputFolderPath, False, lstModInfo)

                If blnSuccess AndAlso mCreateProteinModsFile Then
                    blnSuccess = CreateProteinModsFileWork(strBaseName, fiInputFile, strSynOutputFilePath, strOutputFolderPath)
                End If

                If blnSuccess Then
                    MyBase.OperationComplete()
                End If

            Catch ex As Exception
                SetErrorMessage("Error in clsMSPathFinderResultsProcessor.ProcessFile (2):  " & ex.Message)
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
            End Try

        Catch ex As Exception
            SetErrorMessage("Error in ProcessFile (1):" & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.UnspecifiedError)
        End Try

        Return blnSuccess

    End Function

    Private Sub SortAndWriteFilteredSearchResults(
      swResultFile As StreamWriter,
      lstFilteredSearchResults As List(Of udtMSPathFinderSearchResultType),
      ByRef strErrorLog As String)

        ' Sort udtFilteredSearchResults by ascending QValue, scan, peptide, and protein     
        Dim query = From item In lstFilteredSearchResults Order By item.QValueNum, item.ScanNum, item.Sequence, item.Protein Select item

        Dim intIndex = 1
        For Each result In query
            WriteSearchResultToFile(intIndex, swResultFile, result, strErrorLog)
            intIndex += 1
        Next

    End Sub


    ''' <summary>
    ''' Stores the synopsis file matches for a single scan (typically there will only be one result)
    ''' </summary>
    ''' <param name="lstSearchResults">Search results</param>
    ''' <param name="intStartIndex">Start index for data in this scan</param>
    ''' <param name="intEndIndex">End index for data in this scan</param>
    ''' <param name="lstFilteredSearchResults">Output parmaeter: the actual filtered search results</param>
    ''' <remarks></remarks>
    Private Sub StoreSynMatches(
      lstSearchResults As List(Of udtMSPathFinderSearchResultType),
      intStartIndex As Integer,
      intEndIndex As Integer,
      lstFilteredSearchResults As List(Of udtMSPathFinderSearchResultType))

        Dim intIndex As Integer

        ' If there was more than one result, we could rank them by score
        ' AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex)

        ' The calling procedure already sorted by scan, charge, and Score; no need to re-sort

        ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1)

        ' Now store the matches that pass the filters (will filter on SpecProb once it is implemented)        
        For intIndex = intStartIndex To intEndIndex
            'If lstSearchResults(intIndex).PValueNum <= mMSGFDBSynopsisFilePValueThreshold OrElse
            '   lstSearchResults(intIndex).SpecProbNum <= mMSGFDBSynopsisFileSpecProbThreshold Then
            lstFilteredSearchResults.Add(lstSearchResults(intIndex))
            'End If
        Next intIndex

    End Sub
    Private Sub WriteSynFHTFileHeader(
      swResultFile As StreamWriter,
      ByRef strErrorLog As String)

        ' Write out the header line for synopsis / first hits files
        Try
            Dim lstData As New List(Of String)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ResultID)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Scan)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Charge)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_MostAbundantIsotopeMz)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Mass)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Sequence)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Modifications)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Composition)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_Protein)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinDesc)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinLength)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueStart)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueEnd)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_MatchedFragments)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_QValue)
            lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_PepQValue)

            ' Future: lstData.Add(clsPHRPParserMSPathFinder.DATA_COLUMN_SpecProb)

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
    ''' <remarks></remarks>
    Private Sub WriteSearchResultToFile(
      intResultID As Integer,
      swResultFile As StreamWriter,
      ByRef udtSearchResult As udtMSPathFinderSearchResultType,
      ByRef strErrorLog As String)

        Try

            Dim lstData As New List(Of String)
            lstData.Add(intResultID.ToString())
            lstData.Add(udtSearchResult.Scan)
            lstData.Add(udtSearchResult.Charge)
            lstData.Add(udtSearchResult.MostAbundantIsotopeMz)
            lstData.Add(udtSearchResult.CalculatedMonoMass)
            lstData.Add(udtSearchResult.PrefixResidue & "." & udtSearchResult.Sequence & "." & udtSearchResult.SuffixResidue)
            lstData.Add(udtSearchResult.Modifications)
            lstData.Add(udtSearchResult.Composition)
            lstData.Add(udtSearchResult.Protein)
            lstData.Add(udtSearchResult.ProteinDesc)
            lstData.Add(udtSearchResult.ProteinLength)
            lstData.Add(udtSearchResult.ResidueStart)
            lstData.Add(udtSearchResult.ResidueEnd)
            lstData.Add(udtSearchResult.NumMatchedFragments)
            lstData.Add(udtSearchResult.QValue)
            lstData.Add(udtSearchResult.PepQValue)
            ' Future: lstData.Add(udtSearchResult.SpecProb)

            swResultFile.WriteLine(CollapseList(lstData))

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

#Region "Event Handlers"
    Private Sub ModExtractorErrorHandler(errMsg As String)
        SetErrorMessage(errMsg)
        SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
    End Sub

    Private Sub ModExtractorWarningHandler(warningMsg As String)
        ReportWarning(warningMsg)
    End Sub
#End Region


#Region "IComparer Classes"

    Protected Class MSPathFinderSearchResultsComparerScanChargeScorePeptide
        Implements IComparer(Of udtMSPathFinderSearchResultType)

        Public Function Compare(x As udtMSPathFinderSearchResultType, y As udtMSPathFinderSearchResultType) As Integer Implements IComparer(Of udtMSPathFinderSearchResultType).Compare

            ' First sort on Scan
            If x.ScanNum > y.ScanNum Then
                Return 1
            ElseIf x.ScanNum < y.ScanNum Then
                Return -1
            Else
                ' Scan is the same; check QValue (future: SpecProb)
                If x.QValueNum > y.QValueNum Then
                    Return 1
                ElseIf x.QValueNum < y.QValueNum Then
                    Return -1
                Else
                    ' QValue is the same; check sequence
                    If x.Sequence > y.Sequence Then
                        Return 1
                    ElseIf x.Sequence < y.Sequence Then
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

        End Function

    End Class

#End Region

End Class
