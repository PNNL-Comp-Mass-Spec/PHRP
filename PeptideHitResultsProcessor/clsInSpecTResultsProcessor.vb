Option Strict On

' This class reads in an InSpecT results file (txt format) and creates 
' a tab-delimited text file with the data.  It will insert modification symbols
' into the peptide sequences for modified peptides.
'
' The modification definition information is determined from the InSpecT Input 
' Parameters section at the end of the InSpecT results file.
'
' -------------------------------------------------------------------------------
' Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
' Program started August 12, 2008
'
' E-mail: john.sandoval@pnl.gov
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'
' Notice: This computer software was prepared by Battelle Memorial Institute, 
' hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
' Department of Energy (DOE).  All rights in the computer software are reserved 
' by DOE on behalf of the United States Government and the Contractor as 
' provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
' SOFTWARE.  This notice including this sentence must appear on any copies of 
' this computer software.

Public Class clsInSpecTResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "September 26, 2008"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_INSPECT_FILE As String = "_inspect"

    Public Const N_TERMINUS_SYMBOL_INSPECT As String = "*."
    Public Const C_TERMINUS_SYMBOL_INSPECT As String = ".*"

    'Private Const FSCORE_THRESHOLD As Single = 0

    ' When writing the synopsis file, we keep data with a pValue below 0.2 Or a TotalPRMScore over 50; this is an OR, not an AND
    Private Const PVALUE_THRESHOLD As Single = 0.2
    Private Const TOTALPRMSCORE_THRESHOLD As Single = 50

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Private Const DTA_FILENAME_SCAN_NUMBER_REGEX As String = "(\d+)\.\d+\.\d+\.dta"
    Private Const REGEX_OPTIONS As Text.RegularExpressions.RegexOptions = Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.Singleline Or Text.RegularExpressions.RegexOptions.IgnoreCase

    Public Enum eInspectResultsFileColumns As Integer
        SpectrumFile = 0
        Scan = 1
        Annotation = 2
        Protein = 3
        Charge = 4
        MQScore = 5
        Length = 6
        TotalPRMScore = 7
        MedianPRMScore = 8
        FractionY = 9
        FractionB = 10
        Intensity = 11
        NTT = 12
        pvalue = 13
        FScore = 14
        DeltaScore = 15
        DeltaScoreOther = 16
        RecordNumber = 17
        DBFilePos = 18
        SpecFilePos = 19
    End Enum

    Public Enum eInspectSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        Peptide = 2
        Protein = 3
        Charge = 4
        MQScore = 5
        Length = 6
        TotalPRMScore = 7
        MedianPRMScore = 8
        FractionY = 9
        FractionB = 10
        Intensity = 11
        NTT = 12
        PValue = 13
        FScore = 14
        DeltaScore = 15
        DeltaScoreOther = 16
        RecordNumber = 17
        DBFilePos = 18
        SpecFilePos = 19
    End Enum
#End Region

#Region "Structures"
    Protected Structure udtInspectSearchResultType
        Public SpectrumFile As String
        Public Scan As String
        Public ScanNum As Integer
        Public PeptideAnnotation As String
        Public Protein As String
        Public Charge As String
        Public ChargeNum As Short
        Public MQScore As String
        Public Length As String
        Public TotalPRMScore As String          ' Higher values are better scores
        Public TotalPRMScoreNum As Single       ' We store the value of the string for quick reference when sorting
        Public MedianPRMScore As String
        Public FractionY As String
        Public FractionB As String
        Public Intensity As String
        Public NTT As String
        Public pValue As String                 ' Lower values are better scores
        Public PValueNum As Single
        Public FScore As String                 ' Higher values are better scores
        Public FScoreNum As Single
        Public DeltaScore As String
        Public DeltaScoreOther As String
        Public RecordNumber As String
        Public DBFilePos As String
        Public SpecFilePos As String

        Public Sub Clear()
            SpectrumFile = String.Empty
            Scan = String.Empty
            ScanNum = 0
            PeptideAnnotation = String.Empty
            Protein = String.Empty
            Charge = String.Empty
            ChargeNum = 0
            MQScore = String.Empty
            Length = String.Empty
            TotalPRMScore = String.Empty
            TotalPRMScoreNum = 0
            MedianPRMScore = String.Empty
            FractionY = String.Empty
            FractionB = String.Empty
            Intensity = String.Empty
            NTT = String.Empty
            pValue = String.Empty
            PValueNum = 0
            FScore = String.Empty
            FScoreNum = 0
            DeltaScore = String.Empty
            DeltaScoreOther = String.Empty
            RecordNumber = String.Empty
            DBFilePos = String.Empty
            SpecFilePos = String.Empty
        End Sub
    End Structure

#End Region

#Region "Classwide Variables"
    Protected mSortFHTandSynFiles As Boolean
#End Region

#Region "Properties"
    Public Property SortFHTandSynFiles() As Boolean
        Get
            Return mSortFHTandSynFiles
        End Get
        Set(ByVal value As Boolean)
            mSortFHTandSynFiles = value
        End Set
    End Property
#End Region

    Private Sub AddCurrentRecordToSearchResults(ByRef intCurrentScanResultsCount As Integer, _
                                                     ByRef udtSearchResultsCurrentScan() As udtInspectSearchResultType, _
                                                     ByRef udtSearchResult As udtInspectSearchResultType, _
                                                     ByRef strErrorLog As String)

        If intCurrentScanResultsCount >= udtSearchResultsCurrentScan.Length Then
            ReDim Preserve udtSearchResultsCurrentScan(udtSearchResultsCurrentScan.Length * 2 - 1)
        End If

        udtSearchResultsCurrentScan(intCurrentScanResultsCount) = udtSearchResult
        intCurrentScanResultsCount += 1

    End Sub

    Private Sub AddDynamicAndStaticResidueMods(ByRef objSearchResult As clsSearchResultsInSpecT, ByVal blnUpdateModOccurrenceCounts As Boolean)
        ' Step through .PeptideSequenceWithMods
        ' For each residue, check if a static mod is defined that affects that residue
        ' For each mod symbol, determine the modification and add to objSearchResult

        Dim intIndex As Integer, intModIndex As Integer
        Dim chChar As Char
        Dim objModificationDefinition As clsModificationDefinition

        Dim strSequence As String
        Dim chMostRecentLetter As Char
        Dim intResidueLocInPeptide As Integer

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
                    .SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                Else
                    ' We found a modification symbol but chMostRecentLetter is not a letter
                    ' Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                End If

            Next intIndex
        End With
    End Sub

    Private Function AddModificationsAndComputeMass(ByRef objSearchResult As clsSearchResultsInSpecT, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
        Const ALLOW_DUPLICATE_MOD_ON_TERMINUS As Boolean = True

        Dim blnSuccess As Boolean

        Try
            ' If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
            objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

            ' Parse .PeptideSequenceWithMods to determine the modified residues present
            AddDynamicAndStaticResidueMods(objSearchResult, blnUpdateModOccurrenceCounts)

            ' Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
            ' Since Sequest allows a terminal peptide residue to be modified twice, we'll allow that to happen,
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

    Protected Function CreateFHTorSYNResultsFile(ByVal strInputFilePath As String, _
                                                 ByVal strOutputFilePath As String, _
                                                 ByVal blnWriteSynFile As Boolean, _
                                                 Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' This routine creates a first hits file or synopsis file from the output from InSpecT

        ' If blnWriteSynFile = False, then writes the first-hits file, which writes the record with the highest TotalPRMScore
        ' If blnWriteSynFile = True, then writes the synopsis file, which writes every record with p-value below a set threshold

        Dim srDataFile As System.IO.StreamReader
        Dim swResultFile As System.IO.StreamWriter

        Dim intPreviousScan As Integer

        Dim strLineIn As String
        Dim protein As String = String.Empty

        Dim udtSearchResult As udtInspectSearchResultType

        Dim intCurrentScanResultsCount As Integer
        Dim udtSearchResultsCurrentScan() As udtInspectSearchResultType

        Dim intFilteredSearchResultCount As Integer
        Dim udtFilteredSearchResults() As udtInspectSearchResultType

        Dim intResultsProcessed As Integer
        Dim intResultID As Integer = 0

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean

        Dim strErrorLog As String

        Try
            ' Initialize variables
            intPreviousScan = Int32.MinValue

            Try
                ' Open the input file and parse it
                ' Initialize the stream reader and the stream Text writer
                srDataFile = New System.IO.StreamReader(strInputFilePath)

                swResultFile = New System.IO.StreamWriter(strOutputFilePath)

                strErrorLog = String.Empty
                intResultsProcessed = 0

                ' Initialize array that will hold all of the records for a given scan
                intCurrentScanResultsCount = 0
                ReDim udtSearchResultsCurrentScan(9)

                ' Initialize the array that will hold all of the records that will ultimately be written out to disk
                intFilteredSearchResultCount = 0
                ReDim udtFilteredSearchResults(999)

                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing
                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then
                        ' Initialize udtSearchResult
                        udtSearchResult.Clear()

                        blnValidSearchResult = ParseInSpecTResultsFileEntry(strLineIn, udtSearchResult, strErrorLog, intResultsProcessed)

                        If Not blnValidSearchResult Then
                            If intResultsProcessed = 0 Then
                                ' This is the first line; write it as the header
                                WriteFileHeader(strLineIn, swResultFile, strErrorLog)
                            End If
                        Else
                            If blnWriteSynFile Then
                                ' Synopsis file
                                ' If udtSearchResult.FScoreNum >= FSCORE_THRESHOLD OrElse _
                                If udtSearchResult.PValueNum <= PVALUE_THRESHOLD OrElse _
                                   udtSearchResult.TotalPRMScoreNum >= TOTALPRMSCORE_THRESHOLD Then
                                    StoreOrWriteSearchResult(swResultFile, intResultID, udtSearchResult, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                                End If
                            Else
                                ' First-hits file
                                If intPreviousScan <> Int32.MinValue AndAlso intPreviousScan <> udtSearchResult.ScanNum Then
                                    ' New scan encountered; sort and filter the data in udtSearchResultsCurrentScan, then call StoreOrWriteSearchResult
                                    StoreTopFHTMatch(swResultFile, intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                                    intCurrentScanResultsCount = 0
                                End If

                                AddCurrentRecordToSearchResults(intCurrentScanResultsCount, udtSearchResultsCurrentScan, udtSearchResult, strErrorLog)

                                intPreviousScan = udtSearchResult.ScanNum
                            End If

                        End If

                        ' Update the progress
                        UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))
                        intResultsProcessed += 1
                    End If
                Loop

                If Not blnWriteSynFile Then
                    ' Store the last record
                    If intCurrentScanResultsCount > 0 Then
                        StoreTopFHTMatch(swResultFile, intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                        intCurrentScanResultsCount = 0
                    End If
                End If

                If mSortFHTandSynFiles Then
                    ' Sort the data in udtFilteredSearchResults then write out to disk
                    SortAndWriteFilteredSearchResults(swResultFile, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
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

    Private Function CIntSafe(ByVal strValue As String, ByVal intDefaultValue As Integer) As Integer
        Try
            Return Integer.Parse(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return intDefaultValue
        End Try
    End Function

    Private Function CSngSafe(ByVal strValue As String, ByVal sngDefaultValue As Single) As Single
        Try
            Return Single.Parse(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return sngDefaultValue
        End Try
    End Function

    Private Function ExtractScanNumFromDTAName(ByVal spectrumFile As String) As String
        Static reScanNumberRegEx As New System.Text.RegularExpressions.Regex(DTA_FILENAME_SCAN_NUMBER_REGEX, REGEX_OPTIONS)

        Dim scanNum As String = String.Empty

        ' See if strValue resembles a .Dta file name
        ' For example, "MyDataset.300.300.2.dta"

        Try
            With reScanNumberRegEx.Match(spectrumFile)
                If .Success AndAlso .Groups.Count > 1 Then
                    scanNum = .Groups(1).Value
                End If
            End With
        Catch ex As Exception
            ' Ignore errors here
            scanNum = "0"
        End Try

        Return scanNum

    End Function

    Private Sub InitializeLocalVariables()
        mSortFHTandSynFiles = True
    End Sub


    Protected Function ParseInSpecTResultsFile(ByVal strInputFilePath As String, ByVal strOutputFilePath As String, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim srDataFile As System.IO.StreamReader

        Dim strPreviousTotalPRMScore As String

        ' Note that Inspect synopsis files are normally sorted on TotalPRMScore descending
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file,
        '  we will keep track of the scan, charge, and peptide information parsed for each unique TotalPRMScore encountered

        Dim htPeptidesFoundForTotalPRMScoreLevel As Hashtable

        Dim strKey As String

        Dim strLineIn As String
        Dim strModificationSummaryFilePath As String

        Dim objSearchResult As clsSearchResultsInSpecT

        Dim intResultsProcessed As Integer

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean
        Dim blnFirstMatchForGroup As Boolean

        Dim strErrorLog As String

        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            ' Initialize objSearchResult
            objSearchResult = New clsSearchResultsInSpecT(mPeptideMods)

            ' Initialize htPeptidesFoundForTotalPRMScoreLevel
            htPeptidesFoundForTotalPRMScoreLevel = New Hashtable
            strPreviousTotalPRMScore = String.Empty

            Try
                UpdateSearchResultEnzymeAndTerminusInfo(objSearchResult)

                ' Open the input file and parse it
                ' Initialize the stream reader
                srDataFile = New System.IO.StreamReader(strInputFilePath)

                strErrorLog = String.Empty
                intResultsProcessed = 0

                ' Create the output files
                blnSuccess = MyBase.InitializeSequenceOutputFiles(strInputFilePath)

                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing

                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then
                        blnValidSearchResult = ParseInSpectSynFileEntry(strLineIn, objSearchResult, strErrorLog, intResultsProcessed)

                        If blnValidSearchResult Then
                            strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge

                            If objSearchResult.TotalPRMScore = strPreviousTotalPRMScore Then
                                ' New result has the same TotalPRMScore as the previous results
                                ' See if htPeptidesFoundForTotalPRMScoreLevel contains the peptide, scan, charge, and MH

                                If htPeptidesFoundForTotalPRMScoreLevel.ContainsKey(strKey) Then
                                    blnFirstMatchForGroup = False
                                Else
                                    htPeptidesFoundForTotalPRMScoreLevel.Add(strKey, 1)
                                    blnFirstMatchForGroup = True
                                End If

                            Else
                                ' New TotalPRMScore
                                ' Reset htPeptidesFoundForScan
                                htPeptidesFoundForTotalPRMScoreLevel.Clear()

                                ' Update strPreviousTotalPRMScore
                                strPreviousTotalPRMScore = objSearchResult.TotalPRMScore

                                ' Append a new entry to htPeptidesFoundForScan
                                htPeptidesFoundForTotalPRMScoreLevel.Add(strKey, 1)
                                blnFirstMatchForGroup = True
                            End If


                            blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup)
                            If Not blnSuccess Then
                                If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                                    strErrorLog &= "Error adding modifications to sequence at RowIndex '" & objSearchResult.ResultID & "'" & ControlChars.NewLine
                                End If
                            End If
                            MyBase.SaveResultsFileEntrySeqInfo(CType(objSearchResult, clsSearchResultsBaseClass), blnFirstMatchForGroup)
                        End If

                        ' Update the progress
                        UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))

                        intResultsProcessed += 1
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

    Private Function ParseInSpecTResultsFileEntry(ByRef strLineIn As String, _
                                                  ByRef udtSearchResult As udtInspectSearchResultType, _
                                                  ByRef strErrorLog As String, _
                                                  ByVal intResultsProcessed As Integer) As Boolean

        ' Parses the entries in an Inspect results file
        ' The expected header line is:
        ' #SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos	

        Dim strSplitLine() As String

        Dim blnValidSearchResult As Boolean
        blnValidSearchResult = False

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset udtSearchResult
            udtSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If intResultsProcessed = 0 Then
                ' This is the first line of the file; it may be a header row
                ' Determine this by seeing if any of the first three columns contains a number
                If Not (clsPHRPBaseClass.IsNumber(strSplitLine(0)) OrElse _
                        clsPHRPBaseClass.IsNumber(strSplitLine(1)) OrElse _
                        clsPHRPBaseClass.IsNumber(strSplitLine(2))) Then
                    ' This is a header line; ignore it
                    blnValidSearchResult = False
                    Exit Try
                End If
            End If

            With udtSearchResult
                .SpectrumFile = strSplitLine(eInspectResultsFileColumns.SpectrumFile)
                If strSplitLine(eInspectResultsFileColumns.Scan) = "0" Then
                    .Scan = ExtractScanNumFromDTAName(.SpectrumFile)
                Else
                    .Scan = strSplitLine(eInspectResultsFileColumns.Scan)
                End If
                .ScanNum = CIntSafe(.Scan, 0)

                .PeptideAnnotation = ReplaceTerminus(strSplitLine(eInspectResultsFileColumns.Annotation))
                .Protein = TruncateProteinName(strSplitLine(eInspectResultsFileColumns.Protein))

                .Charge = strSplitLine(eInspectResultsFileColumns.Charge)
                .ChargeNum = CShort(CIntSafe(.Charge, 0))

                .MQScore = strSplitLine(eInspectResultsFileColumns.MQScore)
                .Length = strSplitLine(eInspectResultsFileColumns.Length)

                .TotalPRMScore = strSplitLine(eInspectResultsFileColumns.TotalPRMScore)
                .TotalPRMScoreNum = CSngSafe(.TotalPRMScore, 0)

                .MedianPRMScore = strSplitLine(eInspectResultsFileColumns.MedianPRMScore)
                .FractionY = strSplitLine(eInspectResultsFileColumns.FractionY)
                .FractionB = strSplitLine(eInspectResultsFileColumns.FractionB)
                .Intensity = strSplitLine(eInspectResultsFileColumns.Intensity)
                .NTT = strSplitLine(eInspectResultsFileColumns.NTT)

                .pValue = strSplitLine(eInspectResultsFileColumns.pvalue)
                .PValueNum = CSngSafe(.pValue, 0)

                .FScore = strSplitLine(eInspectResultsFileColumns.FScore)
                .FScoreNum = CSngSafe(.FScore, 0)

                .DeltaScore = strSplitLine(eInspectResultsFileColumns.DeltaScore)
                .DeltaScoreOther = strSplitLine(eInspectResultsFileColumns.DeltaScoreOther)
                .RecordNumber = strSplitLine(eInspectResultsFileColumns.RecordNumber)
                .DBFilePos = strSplitLine(eInspectResultsFileColumns.DBFilePos)
                .SpecFilePos = strSplitLine(eInspectResultsFileColumns.SpecFilePos)
            End With

            blnValidSearchResult = True
        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing InSpecT Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing InSpecT Results in ParseSequestResultsFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    Private Function ParseInSpectSynFileEntry(ByRef strLineIn As String, _
                                              ByRef objSearchResult As clsSearchResultsInSpecT, _
                                              ByRef strErrorLog As String, _
                                              ByVal intResultsProcessed As Integer) As Boolean

        ' Parses the entries in an Inspect Synopsis file
        ' The expected header line is:
        ' ResultID	Scan	Peptide	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	PValue	FScore	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos

        Dim strSplitLine() As String
        Dim strPeptideSequenceWithMods As String

        Dim blnValidSearchResult As Boolean
        blnValidSearchResult = False

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset objSearchResult
            objSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If intResultsProcessed = 0 Then
                ' This is the first line of the file; it may be a header row
                ' Determine this by seeing if any of the first three columns contains a number
                If Not (clsPHRPBaseClass.IsNumber(strSplitLine(0)) OrElse _
                        clsPHRPBaseClass.IsNumber(strSplitLine(1)) OrElse _
                        clsPHRPBaseClass.IsNumber(strSplitLine(2))) Then
                    ' This is a header line; ignore it
                    blnValidSearchResult = False
                    Exit Try
                End If
            End If

            With objSearchResult
                .ResultID = Integer.Parse(strSplitLine(eInspectSynFileColumns.ResultID))
                .Scan = strSplitLine(eInspectSynFileColumns.Scan)
                .Charge = strSplitLine(eInspectSynFileColumns.Charge)
                .PeptideMH = "0"
                .ProteinName = strSplitLine(eInspectSynFileColumns.Protein)
                .MultipleProteinCount = "0"

                ' Set these to 1 and 10000 since Inspect results files do not contain protein sequence information
                ' If we find later that the peptide sequence spans the length of the protein, we'll revise .ProteinSeqResidueNumberEnd as needed
                .ProteinSeqResidueNumberStart = 1
                .ProteinSeqResidueNumberEnd = 10000

                strPeptideSequenceWithMods = strSplitLine(eInspectSynFileColumns.Peptide)

                ' Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                .SetPeptideSequenceWithMods(strPeptideSequenceWithMods, True, True)

                If .PeptidePreResidues.Trim.EndsWith(clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST) Then
                    ' The peptide is at the N-Terminus of the protein
                    .PeptideLocInProteinStart = .ProteinSeqResidueNumberStart
                    .PeptideLocInProteinEnd = .PeptideLocInProteinStart + .PeptideCleanSequence.Length - 1

                    If .PeptidePostResidues.Trim.Chars(0) = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST Then
                        ' The peptide spans the entire length of the protein
                        .ProteinSeqResidueNumberEnd = .PeptideLocInProteinEnd
                    Else
                        If .PeptideLocInProteinEnd > .ProteinSeqResidueNumberEnd Then
                            ' The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                            .ProteinSeqResidueNumberEnd = .PeptideLocInProteinEnd + 1
                        End If
                    End If
                ElseIf .PeptidePostResidues.Trim.StartsWith(clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST) Then
                    ' The peptide is at the C-Terminus of the protein
                    .PeptideLocInProteinEnd = .ProteinSeqResidueNumberEnd
                    .PeptideLocInProteinStart = .PeptideLocInProteinEnd - .PeptideCleanSequence.Length + 1

                    If .PeptideLocInProteinStart < .ProteinSeqResidueNumberStart Then
                        ' The peptide is more than 10000 characters long; this is highly unlikely
                        .ProteinSeqResidueNumberEnd = .ProteinSeqResidueNumberStart + 1 + .PeptideCleanSequence.Length
                        .PeptideLocInProteinEnd = .ProteinSeqResidueNumberEnd
                        .PeptideLocInProteinStart = .PeptideLocInProteinEnd - .PeptideCleanSequence.Length + 1
                    End If
                Else
                    .PeptideLocInProteinStart = .ProteinSeqResidueNumberStart + 1
                    .PeptideLocInProteinEnd = .PeptideLocInProteinStart + .PeptideCleanSequence.Length - 1

                    If .PeptideLocInProteinEnd > .ProteinSeqResidueNumberEnd Then
                        ' The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                        .ProteinSeqResidueNumberEnd = .PeptideLocInProteinEnd + 1
                    End If
                End If

                ' Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                ' If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide 
                ' will all be based on the first protein since Sequest only outputs the prefix and suffix letters for the first protein
                .ComputePeptideCleavageStateInProtein()

                .PeptideDeltaMass = "0"

                .MQScore = strSplitLine(eInspectSynFileColumns.MQScore)
                .Length = strSplitLine(eInspectSynFileColumns.Length)
                .TotalPRMScore = strSplitLine(eInspectSynFileColumns.TotalPRMScore)
                .MedianPRMScore = strSplitLine(eInspectSynFileColumns.MedianPRMScore)
                .FractionY = strSplitLine(eInspectSynFileColumns.FractionY)
                .FractionB = strSplitLine(eInspectSynFileColumns.FractionB)
                .Intensity = strSplitLine(eInspectSynFileColumns.Intensity)
                .NTT = strSplitLine(eInspectSynFileColumns.NTT)
                .pValue = strSplitLine(eInspectSynFileColumns.PValue)
                .FScore = strSplitLine(eInspectSynFileColumns.FScore)
                .DeltaScore = strSplitLine(eInspectSynFileColumns.DeltaScore)
                .DeltaScoreOther = strSplitLine(eInspectSynFileColumns.DeltaScoreOther)
                .RecordNumber = strSplitLine(eInspectSynFileColumns.RecordNumber)
                .DBFilePos = strSplitLine(eInspectSynFileColumns.DBFilePos)
                .SpecFilePos = strSplitLine(eInspectSynFileColumns.SpecFilePos)
            End With

            blnValidSearchResult = True
        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing InSpecT Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing InSpecT Results in ParseSequestResultsFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    ' Main processing function
    Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean
        ' Returns True if success, False if failure

        Dim ioFile As System.IO.FileInfo

        Dim strInputFilePathFull As String
        Dim strOutputFilePath As String
        Dim strSynOutputFilePath As String

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

                If CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                    Try
                        ' Obtain the full path to the input file
                        ioFile = New System.IO.FileInfo(strInputFilePath)
                        strInputFilePathFull = ioFile.FullName

                        ' Create the first hits output file
                        MyBase.ResetProgress("Creating the FHT file")
                        Console.WriteLine()
                        Console.WriteLine(MyBase.ProgressStepDescription)

                        strOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                        strOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strOutputFilePath & SEQUEST_FIRST_HITS_FILE_SUFFIX)

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strOutputFilePath, False)

                        ' Create the synopsis output file
                        MyBase.ResetProgress("Creating the SYN file")
                        Console.WriteLine()
                        Console.WriteLine()
                        Console.WriteLine(MyBase.ProgressStepDescription)

                        'Define the synopsis output file name based on strInputFilePath
                        strSynOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                        strSynOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strSynOutputFilePath & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, True)

                        ' Create the other PHRP-specific files
                        MyBase.ResetProgress("Creating the PHRP files for " & System.IO.Path.GetFileName(strSynOutputFilePath))
                        Console.WriteLine()
                        Console.WriteLine()
                        Console.WriteLine(MyBase.ProgressStepDescription)

                        blnSuccess = ParseInSpecTResultsFile(strSynOutputFilePath, strOutputFilePath, False)

                    Catch ex As Exception
                        SetErrorMessage("Error calling CreateFHTorSYNResultsFile" & ex.Message)
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

    Private Function ReplaceTerminus(ByVal inpString As String) As String

        If inpString.Trim().StartsWith(N_TERMINUS_SYMBOL_INSPECT) Then
            inpString = "-." & inpString.Substring(N_TERMINUS_SYMBOL_INSPECT.Length)
        End If

        If inpString.Trim().EndsWith(C_TERMINUS_SYMBOL_INSPECT) Then
            inpString = inpString.Substring(0, inpString.Length - C_TERMINUS_SYMBOL_INSPECT.Length) & ".-"
        End If

        Return inpString

    End Function

    Protected Sub StoreOrWriteSearchResult(ByRef swResultFile As System.IO.StreamWriter, _
                                           ByRef intResultID As Integer, _
                                           ByRef udtSearchResult As udtInspectSearchResultType, _
                                           ByRef intFilteredSearchResultCount As Integer, _
                                           ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                           ByRef strErrorLog As String)
        If mSortFHTandSynFiles Then
            If intFilteredSearchResultCount = udtFilteredSearchResults.Length Then
                ReDim Preserve udtFilteredSearchResults(udtFilteredSearchResults.Length * 2 - 1)
            End If

            udtFilteredSearchResults(intFilteredSearchResultCount) = udtSearchResult
            intFilteredSearchResultCount += 1
        Else
            intResultID += 1
            WriteSearchResultToFile(intResultID, swResultFile, udtSearchResult, strErrorLog)
        End If
    End Sub

    Private Sub SortAndWriteFilteredSearchResults(ByRef swResultFile As System.IO.StreamWriter, _
                                                  ByVal intFilteredSearchResultCount As Integer, _
                                                  ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                                  ByRef strErrorLog As String)

        Dim intIndex As Integer

        ' Sort udtFilteredSearchResults by descending TotalPRMScore, ascending scan, ascending charge, ascending peptide, and ascending protein
        Array.Sort(udtFilteredSearchResults, 0, intFilteredSearchResultCount, New InspectSearchResultsComparerTotalPRMDescScanChargePeptide)

        For intIndex = 0 To intFilteredSearchResultCount - 1
            WriteSearchResultToFile(intIndex + 1, swResultFile, udtFilteredSearchResults(intIndex), strErrorLog)
        Next intIndex

    End Sub

    Private Sub StoreTopFHTMatch(ByRef swResultFile As System.IO.StreamWriter, _
                                      ByRef intResultID As Integer, _
                                      ByVal intCurrentScanResultsCount As Integer, _
                                      ByRef udtSearchResultsCurrentScan() As udtInspectSearchResultType, _
                                      ByRef intFilteredSearchResultCount As Integer, _
                                      ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                      ByRef strErrorLog As String)

        Dim intIndex As Integer
        Dim intCurrentCharge As Short = Short.MinValue

        ' Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending PValue
        ' All of the data in udtSearchResultsCurrentScan should have the same scan number
        Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, New InspectSearchResultsComparerScanChargeTotalPRMDescPValueDesc)

        ' Now store or write out the first match for each charge for this scan
        For intIndex = 0 To intCurrentScanResultsCount - 1
            If intCurrentCharge = Short.MinValue OrElse intCurrentCharge <> udtSearchResultsCurrentScan(intIndex).ChargeNum Then
                StoreOrWriteSearchResult(swResultFile, intResultID, udtSearchResultsCurrentScan(intIndex), intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                intCurrentCharge = udtSearchResultsCurrentScan(intIndex).ChargeNum
            End If
        Next intIndex

    End Sub

    Private Function TruncateProteinName(ByVal strProteinNameAndDescription As String) As String

        Dim intIndex As Integer

        intIndex = strProteinNameAndDescription.IndexOf(" "c)
        If intIndex > 0 Then
            strProteinNameAndDescription = strProteinNameAndDescription.Substring(0, intIndex)
        End If

        Return strProteinNameAndDescription

    End Function

    Private Sub UpdateSearchResultEnzymeAndTerminusInfo(ByRef objSearchResult As clsSearchResultsInSpecT)
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

    Private Sub WriteFileHeader(ByRef strLineIn As String, _
                                     ByRef swResultFile As System.IO.StreamWriter, _
                                     ByRef strErrorLog As String)
        Dim strSplitLine As String()

        'Replace some header names with new name
        Try
            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            swResultFile.WriteLine("ResultID" & ControlChars.Tab & _
                                   "Scan" & ControlChars.Tab & _
                                   "Peptide" & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.Protein) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.Charge) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.MQScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.Length) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.TotalPRMScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.MedianPRMScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.FractionY) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.FractionB) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.Intensity) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.NTT) & ControlChars.Tab & _
                                   "PValue" & ControlChars.Tab & _
                                   "FScore" & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.DeltaScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.DeltaScoreOther) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.RecordNumber) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.DBFilePos) & ControlChars.Tab & _
                                   strSplitLine(eInspectResultsFileColumns.SpecFilePos))

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing first hits header" & ControlChars.NewLine
            End If
        End Try

    End Sub

    Private Sub WriteSearchResultToFile(ByVal intResultID As Integer, _
                                     ByRef swResultFile As System.IO.StreamWriter, _
                                     ByRef udtSearchResult As udtInspectSearchResultType, _
                                     ByRef strErrorLog As String)

        Try
            swResultFile.WriteLine(intResultID.ToString & ControlChars.Tab & _
                                   udtSearchResult.Scan & ControlChars.Tab & _
                                   udtSearchResult.PeptideAnnotation & ControlChars.Tab & _
                                   udtSearchResult.Protein & ControlChars.Tab & _
                                   udtSearchResult.Charge & ControlChars.Tab & _
                                   udtSearchResult.MQScore & ControlChars.Tab & _
                                   udtSearchResult.Length & ControlChars.Tab & _
                                   udtSearchResult.TotalPRMScore & ControlChars.Tab & _
                                   udtSearchResult.MedianPRMScore & ControlChars.Tab & _
                                   udtSearchResult.FractionY & ControlChars.Tab & _
                                   udtSearchResult.FractionB & ControlChars.Tab & _
                                   udtSearchResult.Intensity & ControlChars.Tab & _
                                   udtSearchResult.NTT & ControlChars.Tab & _
                                   udtSearchResult.pValue & ControlChars.Tab & _
                                   udtSearchResult.FScore & ControlChars.Tab & _
                                   udtSearchResult.DeltaScore & ControlChars.Tab & _
                                   udtSearchResult.DeltaScoreOther & ControlChars.Tab & _
                                   udtSearchResult.RecordNumber & ControlChars.Tab & _
                                   udtSearchResult.DBFilePos & ControlChars.Tab & _
                                   udtSearchResult.SpecFilePos & ControlChars.Tab)

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

    Protected Class InspectSearchResultsComparerTotalPRMDescScanChargePeptide
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtInspectSearchResultType = CType(x, udtInspectSearchResultType)
            Dim yData As udtInspectSearchResultType = CType(y, udtInspectSearchResultType)

            If xData.TotalPRMScoreNum > yData.TotalPRMScoreNum Then
                Return -1
            ElseIf xData.TotalPRMScoreNum < yData.TotalPRMScoreNum Then
                Return 1
            Else
                ' TotalPRMScore is the same; check scan number
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
                        If xData.PeptideAnnotation > yData.PeptideAnnotation Then
                            Return 1
                        ElseIf xData.PeptideAnnotation < yData.PeptideAnnotation Then
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

    Protected Class InspectSearchResultsComparerScanChargeTotalPRMDescPValueDesc
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtInspectSearchResultType = CType(x, udtInspectSearchResultType)
            Dim yData As udtInspectSearchResultType = CType(y, udtInspectSearchResultType)

            If xData.ScanNum > yData.ScanNum Then
                Return 1
            ElseIf xData.ScanNum < yData.ScanNum Then
                Return -1
            Else
                If xData.ChargeNum > yData.ChargeNum Then
                    Return 1
                ElseIf xData.ChargeNum < yData.ChargeNum Then
                    Return -1
                Else
                    ' Charge is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                    If xData.TotalPRMScoreNum > yData.TotalPRMScoreNum Then
                        Return -1
                    ElseIf xData.TotalPRMScoreNum < yData.TotalPRMScoreNum Then
                        Return 1
                    Else
                        ' TotalPRMScore is the same; check P-Value (sort on ascending p-value since lower p-values are better)
                        If xData.PValueNum > yData.PValueNum Then
                            Return 1
                        ElseIf xData.PValueNum < yData.PValueNum Then
                            Return -1
                        Else
                            Return 0
                        End If
                    End If
                End If
            End If

        End Function
    End Class
End Class
