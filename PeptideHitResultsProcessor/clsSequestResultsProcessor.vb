Option Strict On

' This class reads in a synopsis or first hits file (tab-delimited representation
' of the out file data, created from STARSuite Extractor) and creates a new file
' containing columns for cleavage and terminus state, modification information, and
' the monoisotopic mass of each peptide.  The data in the new file is linked to the
' original file by the Row ID number in the original file.  The user must provide a
' modification definition file that specifies the dynamic and/or static modifications
' used in the search.
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 2, 2006
'
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

Public Class clsSequestResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "March 6, 2006"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
    Public Const FILENAME_SUFFIX_FIRST_HITS_FILE As String = "_fht"
    Public Const FILENAME_SUFFIX_SYNOPSIS_FILE As String = "_syn"

    Private Const SYNOPSIS_OR_FIRST_HITS_FILE_COLUMN_COUNT_EXPECTED As Integer = 19
    Private Const ADDITIONAL_COLUMN_COUNT_APPENDED As Integer = 5
    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Public Enum eSynopsisFileColumns As Integer
        RowIndex = 0
        Scan = 1
        NumScans = 2
        Charge = 3
        PeptideMH = 4
        XCorr = 5
        DeltaCn = 6
        Sp = 7
        ProteinName = 8                 ' Aka Reference
        MultipleProteinCount = 9        ' Aka MO = MultipleORFCount; this is 0 if the peptide is in just one protein; 1 if in 2 proteins, etc.
        PeptideSequence = 10            ' This is the sequence with prefix and suffix residues and also with modification symbols
        DeltaCn2 = 11
        RankSP = 12
        RankXC = 13
        DelM = 14
        XcRatio = 15
        PassFilt = 16
        MScore = 17
        NTT = 18                        ' Number of tryptic terminii
        Cleavage_State = 19             ' This column and the ones after it are computed by this program and appened to the input file or saved in a new file
        Terminus_State = 20
        Mod_Count = 21
        Mod_Description = 22
        Monoisotopic_Mass = 23
    End Enum
#End Region

#Region "Structures"
#End Region

#Region "Classwide Variables"

#End Region

#Region "Properties"
#End Region

    Private Sub AddDynamicAndStaticResidueMods(ByRef objSearchResult As clsSearchResultsSequest, ByVal blnUpdateModOccurrenceCounts As Boolean)
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

    Private Function AddModificationsAndComputeMass(ByRef objSearchResult As clsSearchResultsSequest, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
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

    Private Sub InitializeLocalVariables()
        ' Reserved for future use
    End Sub

    Protected Function ParseSynopsisOrFirstHitsFile(ByVal strInputFilePath As String, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim srDataFile As System.IO.StreamReader

        Dim strPreviousXCorr As String

        ' Note that synopsis files are normally sorted on XCorr descending, with lines
        '  duplicated when  peptide search results are mapped to multiple proteins
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file,
        '  we will keep track of the scan, charge, and peptide information parsed for each unique XCorr encountered

        Dim htPeptidesFoundForXCorrLevel As Hashtable

        Dim strKey As String

        Dim strLineIn As String
        Dim strModificationSummaryFilePath As String

        Dim objSearchResult As clsSearchResultsSequest

        Dim intIndex As Integer
        Dim intResultsProcessed As Integer

        Dim objModificationDefinition As clsModificationDefinition

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
            objSearchResult = New clsSearchResultsSequest(mPeptideMods)

            ' Initialize htPeptidesFoundForXCorrLevel
            htPeptidesFoundForXCorrLevel = New Hashtable
            strPreviousXCorr = String.Empty

            Try
                UpdateSearchResultEnzymeAndTerminusInfo(objSearchResult)

                ' Open the input file and parse it
                ' Initialize the stream reader and the XML Text Reader
                srDataFile = New System.IO.StreamReader(strInputFilePath)

                strErrorLog = String.Empty
                intResultsProcessed = 0

                ' Create the output files
                blnSuccess = MyBase.InitializeSequenceOutputFiles(strInputFilePath)

                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing

                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then
                        blnValidSearchResult = ParseSequestResultsFileEntry(strLineIn, objSearchResult, strErrorLog, intResultsProcessed)

                        If blnValidSearchResult Then
                            strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.NumScans & "_" & objSearchResult.Charge & "_" & objSearchResult.PeptideMH

                            If objSearchResult.PeptideXCorr = strPreviousXCorr Then
                                ' New result has the same XCorr as the previous results
                                ' See if htPeptidesFoundForXCorrLevel contains the peptide, scan, charge, and MH

                                If htPeptidesFoundForXCorrLevel.ContainsKey(strKey) Then
                                    blnFirstMatchForGroup = False
                                Else
                                    htPeptidesFoundForXCorrLevel.Add(strKey, 1)
                                    blnFirstMatchForGroup = True
                                End If

                            Else
                                ' New XCorr
                                ' Reset htPeptidesFoundForScan
                                htPeptidesFoundForXCorrLevel.Clear()

                                ' Update strPreviousXCorr
                                strPreviousXCorr = objSearchResult.PeptideXCorr

                                ' Append a new entry to htPeptidesFoundForScan
                                htPeptidesFoundForXCorrLevel.Add(strKey, 1)
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

    Private Function ParseSequestResultsFileEntry(ByRef strLineIn As String, ByRef objSearchResult As clsSearchResultsSequest, ByRef strErrorLog As String, ByVal intResultsProcessed As Integer) As Boolean

        Dim strSplitLine() As String
        Dim strPeptideSequenceWithMods As String
        Dim strPeptidePreResidues As String
        Dim strPeptidePostResidues As String

        Dim strPrimarySequence As String

        Dim blnValidSearchResult As Boolean
        blnValidSearchResult = False

        ' The following are the headers in the synopsis file
        ' "RowIndex"
        ' "Scan"
        ' "NumScans"
        ' "Charge"
        ' "MH"
        ' "XCorr"
        ' "DeltaCn"
        ' "Sp"
        ' "Reference"
        ' "MO"                  ' Multiple protein count: 0 if the peptide is in 1 protein, 1 if the peptide is in 2 proteins, etc.
        ' "Peptide"
        ' "DeltaCn2"
        ' "RankSP"
        ' "RankXC"
        ' "DelM"
        ' "XcRatio"
        ' "PassFilt"
        ' "MScore"
        ' "NTT"


        ' "Cleavage_State"
        ' "Terminus_State"
        ' "Mod_Count"
        ' "Mod_Description"
        ' "Monoisotopic_Mass"
        ' "Cleavage_State"
        ' "Terminus_State"
        ' "Mod_Count"
        ' "Mod_Description"
        ' "Monoisotopic_Mass"

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset objSearchResult
            objSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)
            If strSplitLine.Length < SYNOPSIS_OR_FIRST_HITS_FILE_COLUMN_COUNT_EXPECTED Then
                Exit Try
            End If

            If intResultsProcessed = 0 Then
                ' This is the first line of the file; it may be a header row
                ' Determine this by seeing if any of the first three columns contains a number
                If Not (clsPHRPBaseClass.IsNumber(strSplitLine(0)) OrElse _
                   clsPHRPBaseClass.IsNumber(strSplitLine(1)) OrElse _
                   clsPHRPBaseClass.IsNumber(strSplitLine(2))) Then
                    ' This is a header line; skip it
                    blnValidSearchResult = False
                    Exit Try
                End If
            End If

            With objSearchResult
                .ResultID = Integer.Parse(strSplitLine(eSynopsisFileColumns.RowIndex))
                .Scan = strSplitLine(eSynopsisFileColumns.Scan)
                .NumScans = strSplitLine(eSynopsisFileColumns.NumScans)
                .Charge = strSplitLine(eSynopsisFileColumns.Charge)
                .PeptideMH = strSplitLine(eSynopsisFileColumns.PeptideMH)
                .PeptideXCorr = strSplitLine(eSynopsisFileColumns.XCorr)
                .PeptideDeltaCn = strSplitLine(eSynopsisFileColumns.DeltaCn)
                .PeptideSp = strSplitLine(eSynopsisFileColumns.Sp)
                .ProteinName = strSplitLine(eSynopsisFileColumns.ProteinName)
                .MultipleProteinCount = strSplitLine(eSynopsisFileColumns.MultipleProteinCount)

                ' Set these to 1 and 10000 since Sequest results files do not contain protein sequence information
                ' If we find later that the peptide sequence spans the length of the protein, we'll revise .ProteinSeqResidueNumberEnd as needed
                .ProteinSeqResidueNumberStart = 1
                .ProteinSeqResidueNumberEnd = 10000

                strPeptideSequenceWithMods = strSplitLine(eSynopsisFileColumns.PeptideSequence)

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

                .PeptideDeltaCn2 = strSplitLine(eSynopsisFileColumns.DeltaCn2)
                .PeptideRankSP = strSplitLine(eSynopsisFileColumns.RankSP)
                .PeptideRankXC = strSplitLine(eSynopsisFileColumns.RankXC)
                .PeptideDeltaMass = strSplitLine(eSynopsisFileColumns.DelM)
                .PeptideXcRatio = strSplitLine(eSynopsisFileColumns.XcRatio)
                .PeptidePassFilt = strSplitLine(eSynopsisFileColumns.PassFilt)
                .PeptideMScore = strSplitLine(eSynopsisFileColumns.MScore)
                .PeptideNTT = strSplitLine(eSynopsisFileColumns.NTT)
            End With

            blnValidSearchResult = True
        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing Sequest Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing Sequest Results in ParseSequestResultsFileEntry" & ControlChars.NewLine
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
        Dim ioOutputFile As System.IO.FileInfo

        Dim strInputFilePathFull As String

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

                MyBase.mProgressStepDescription = "Parsing " & System.IO.Path.GetFileName(strInputFilePath)
                MyBase.ResetProgress()

                If CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                    Try
                        ' Obtain the full path to the input file
                        ioInputFile = New System.IO.FileInfo(strInputFilePath)
                        strInputFilePathFull = ioInputFile.FullName

                        blnSuccess = ParseSynopsisOrFirstHitsFile(strInputFilePathFull, False)
                    Catch ex As Exception
                        SetErrorMessage("Error calling ParseSynopsisOrFirstHitsFile" & ex.message)
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

    Private Sub UpdateSearchResultEnzymeAndTerminusInfo(ByRef objSearchResult As clsSearchResultsSequest)
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

    ''Private Sub SaveSequestResultsFileEntry(ByRef objSearchResult As clsSearchResultsSequest, ByRef swSynopsisOutputFile As System.IO.StreamWriter)

    ''    ' Write the results to the output file
    ''    With objSearchResult
    ''        swSeqInfoFile.WriteLine(.ResultID & SEP_CHAR & _
    ''                            CInt(.PeptideCleavageState).ToString & SEP_CHAR & _
    ''                            CInt(.PeptideTerminusState).ToString & SEP_CHAR & _
    ''                            .SearchResultModificationCount & SEP_CHAR & _
    ''                            .PeptideModDescription & SEP_CHAR & _
    ''                            .PeptideMonoisotopicMass.ToString("0.000000000"))

    ''        '' Legacy code to re-create the synopsis file
    ''        '' swSynopsisOutputFile.WriteLine(.ResultID & SEP_CHAR & _
    ''        ''                    .Scan & SEP_CHAR & _
    ''        ''                    .NumScans & SEP_CHAR & _
    ''        ''                    .Charge & SEP_CHAR & _
    ''        ''                    .PeptideMH & SEP_CHAR & _
    ''        ''                    .PeptideXCorr & SEP_CHAR & _
    ''        ''                    .PeptideDeltaCn & SEP_CHAR & _
    ''        ''                    .PeptideSp & SEP_CHAR & _
    ''        ''                    .ProteinName & SEP_CHAR & _
    ''        ''                    .MultipleProteinCount & SEP_CHAR & _
    ''        ''                    .SequenceWithPrefixAndSuffix(True) & SEP_CHAR & _
    ''        ''                    .PeptideDeltaCn2 & SEP_CHAR & _
    ''        ''                    .PeptideRankSP & SEP_CHAR & _
    ''        ''                    .PeptideRankXC & SEP_CHAR & _
    ''        ''                    .PeptideDeltaMass & SEP_CHAR & _
    ''        ''                    .PeptideXcRatio & SEP_CHAR & _
    ''        ''                    .PeptidePassFilt & SEP_CHAR & _
    ''        ''                    .PeptideMScore & SEP_CHAR & _
    ''        ''                    .PeptideNTT & SEP_CHAR & _
    ''        ''                    CInt(.PeptideCleavageState).ToString & SEP_CHAR & _
    ''        ''                    CInt(.PeptideTerminusState).ToString & SEP_CHAR & _
    ''        ''                    .SearchResultModificationCount & SEP_CHAR & _
    ''        ''                    .PeptideModDescription & SEP_CHAR & _
    ''        ''                    .PeptideMonoisotopicMass.ToString("0.000000000"))

    ''    End With

    ''End Sub

End Class
