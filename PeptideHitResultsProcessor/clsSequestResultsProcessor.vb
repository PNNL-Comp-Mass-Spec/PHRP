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
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

Imports PHRPReader
Imports System.IO

Public Class clsSequestResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
		MyBase.mFileDate = "June 28, 2013"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
    Public Const FILENAME_SUFFIX_FIRST_HITS_FILE As String = "_fht"
    Public Const FILENAME_SUFFIX_SYNOPSIS_FILE As String = "_syn"

    Private Const SEQUEST_SYN_FILE_MIN_COL_COUNT As Integer = 5
    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    ' These columns correspond to the tab-delimited file created directly by MSGF-DB
    Protected Const SequestSynopsisFileColCount As Integer = 27
    Public Enum eSequestSynopsisFileColumns As Integer
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
        PassFilt = 16                   ' Legacy/unused
        MScore = 17                     ' Legacy/unused
        NTT = 18                        ' Number of tryptic terminii
        IonsObserved = 19               ' Added in August 2011
        IonsExpected = 20               ' Added in August 2011 
        DelMPPM = 21                    ' Added in August 2011
        Cleavage_State = 22             ' This column and the ones after it are computed by this program and appended to the input file or saved in a new file
        Terminus_State = 23
        Mod_Count = 24
        Mod_Description = 25
        Monoisotopic_Mass = 26
    End Enum

#End Region

#Region "Structures"
#End Region

#Region "Classwide Variables"

#End Region

#Region "Properties"
#End Region

    Private Function AddDynamicAndStaticResidueMods(
      objSearchResult As clsSearchResultsSequest,
      blnUpdateModOccurrenceCounts As Boolean) As Boolean

        ' Step through .PeptideSequenceWithMods
        ' For each residue, check if a static mod is defined that affects that residue
        ' For each mod symbol, determine the modification and add to objSearchResult

        Dim intIndex As Integer, intModIndex As Integer
        Dim chChar As Char
        Dim objModificationDefinition As clsModificationDefinition

        Dim strSequence As String
        Dim chMostRecentLetter As Char
        Dim intResidueLocInPeptide As Integer

        Dim blnSuccess As Boolean

        chMostRecentLetter = "-"c
        intResidueLocInPeptide = 0

        ' Assume success for now
        blnSuccess = True

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
                                blnSuccess = .SearchResultAddModification(objModificationDefinition, chChar, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)

                                If Not blnSuccess Then
                                    ' Error adding this static mod
                                    SetErrorCode(ePHRPErrorCodes.UnspecifiedError)
                                    mErrorMessage = "Error calling objSearchResult.SearchResultAddModification for peptide '" & strSequence & "': " & .ErrorMessage
                                    Exit For
                                End If
                            End If
                        End If
                    Next intModIndex
                ElseIf IsLetterAtoZ(chMostRecentLetter) Then
                    blnSuccess = .SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)

                    If Not blnSuccess Then
                        ' Error adding this dynamic mod
                        SetErrorCode(ePHRPErrorCodes.UnspecifiedError)
                        mErrorMessage = "Error calling objSearchResult.SearchResultAddDynamicModification for peptide '" & strSequence & "': " & .ErrorMessage
                        Exit For
                    End If

                Else
                    ' We found a modification symbol but chMostRecentLetter is not a letter
                    ' Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                End If

            Next intIndex
        End With

        Return blnSuccess
    End Function

    Private Function AddModificationsAndComputeMass(objSearchResult As clsSearchResultsSequest, blnUpdateModOccurrenceCounts As Boolean) As Boolean
        Const ALLOW_DUPLICATE_MOD_ON_TERMINUS = True

        Dim blnSuccess As Boolean

        Try
            ' Assume success for now
            blnSuccess = True

            ' If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
            objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

            ' Parse .PeptideSequenceWithMods to determine the modified residues present
            If Not AddDynamicAndStaticResidueMods(objSearchResult, blnUpdateModOccurrenceCounts) Then
                blnSuccess = False
            End If

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

        Catch ex As Exception
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Overrides Function ConstructPepToProteinMapFilePath(strInputFilePath As String, strOutputFolderPath As String, MTS As Boolean) As String
        Dim strPepToProteinMapFilePath As String = String.Empty

        strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath)
        If strPepToProteinMapFilePath.ToLower().EndsWith("_syn") OrElse strPepToProteinMapFilePath.ToLower().EndsWith("_fht") Then
            ' Remove _syn or _fht
            strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4)
        End If

        Return MyBase.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS)

    End Function

    Private Sub InitializeLocalVariables()
        ' Reserved for future use
    End Sub

    Protected Function ParseSynopsisOrFirstHitsFile(strInputFilePath As String, strOutputFolderPath As String, blnResetMassCorrectionTagsAndModificationDefinitions As Boolean) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim strPreviousXCorr As String

        ' Note that synopsis files are normally sorted on XCorr descending, with lines
        '  duplicated when peptide search results are mapped to multiple proteins
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file,
        '  we will keep track of the scan, charge, and peptide information parsed for each unique XCorr encountered

        Dim htPeptidesFoundForXCorrLevel As Hashtable

        Dim strKey As String

        Dim strLineIn As String
        Dim strModificationSummaryFilePath As String

        Dim objSearchResult As clsSearchResultsSequest

        Dim intResultsProcessed As Integer
        Dim sngPercentComplete As Single

        Dim blnSuccess As Boolean
        Dim blnDataLine As Boolean
        Dim blnValidSearchResult As Boolean
        Dim blnFirstMatchForGroup As Boolean

        Dim blnHeaderParsed As Boolean
        Dim intColumnMapping() As Integer = Nothing

        Dim strErrorLog As String = String.Empty

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
                objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(mEnzymeMatchSpec, mPeptideNTerminusMassChange, mPeptideCTerminusMassChange)

                ' Open the input file and parse it
                ' Initialize the stream reader
                Using srDataFile As StreamReader = New StreamReader(strInputFilePath)

                    strErrorLog = String.Empty
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

                        blnDataLine = True

                        If Not blnHeaderParsed Then
                            blnSuccess = ParseSequestSynFileHeaderLine(strLineIn, intColumnMapping)
                            If blnSuccess Then
                                blnDataLine = False
                            Else
                                ' Error parsing header; assume this is a data line
                                blnDataLine = True
                            End If
                            blnHeaderParsed = True
                        End If

                        If blnDataLine Then
                            blnValidSearchResult = ParseSequestResultsFileEntry(strLineIn, intColumnMapping, objSearchResult, strErrorLog)
                        Else
                            blnValidSearchResult = False
                        End If

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
                                ' Reset htPeptidesFoundForXCorrLevel
                                htPeptidesFoundForXCorrLevel.Clear()

                                ' Update strPreviousXCorr
                                strPreviousXCorr = objSearchResult.PeptideXCorr

                                ' Append a new entry to htPeptidesFoundForXCorrLevel
                                htPeptidesFoundForXCorrLevel.Add(strKey, 1)
                                blnFirstMatchForGroup = True
                            End If


                            blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup)
                            If Not blnSuccess Then
                                If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                                    strErrorLog &= "Error adding modifications to sequence at RowIndex '" & objSearchResult.ResultID & "'"
                                    If Not mErrorMessage Is Nothing AndAlso mErrorMessage.Length > 0 Then
                                        strErrorLog &= ": " & mErrorMessage
                                        mErrorMessage = String.Empty
                                    End If
                                    strErrorLog &= ControlChars.NewLine
                                End If
                            End If
                            MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResult, clsSearchResultsBaseClass), blnFirstMatchForGroup)
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
                    Dim fiInputFile As FileInfo = New FileInfo(strInputFilePath)
                    strModificationSummaryFilePath = Path.GetFileName(MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_MOD_SUMMARY))
                    strModificationSummaryFilePath = Path.Combine(strOutputFolderPath, strModificationSummaryFilePath)

                    SaveModificationSummaryFile(strModificationSummaryFilePath)
                End If

                ' Inform the user if any errors occurred
                If strErrorLog.Length > 0 Then
                    SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
                    blnSuccess = False
                Else
                    blnSuccess = True
                End If

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

    Private Function ParseSequestResultsFileEntry(
      ByRef strLineIn As String,
      ByRef intColumnMapping() As Integer,
      objSearchResult As clsSearchResultsSequest,
      ByRef strErrorLog As String) As Boolean

        Dim strSplitLine() As String = Nothing
        Dim strPeptideSequenceWithMods As String = String.Empty

        Dim blnValidSearchResult As Boolean
        blnValidSearchResult = False

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset objSearchResult
            objSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)
            If strSplitLine.Length < SEQUEST_SYN_FILE_MIN_COL_COUNT Then
                Exit Try
            End If

            With objSearchResult
                If Not GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.RowIndex), .ResultID) Then
                    ReportError("RowIndex column is missing or invalid", True)
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.Scan), .Scan)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.NumScans), .NumScans)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.Charge), .Charge)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.PeptideMH), .PeptideMH)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.XCorr), .PeptideXCorr)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.DeltaCn), .PeptideDeltaCn)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.Sp), .PeptideSp)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.ProteinName), .ProteinName)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.MultipleProteinCount), .MultipleProteinCount)

                If Not GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.PeptideSequence), strPeptideSequenceWithMods) Then
                    ReportError("Peptide column is missing or invalid", True)
                End If

                ' Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                .SetPeptideSequenceWithMods(strPeptideSequenceWithMods, True, True)

            End With

            Dim objSearchResultBase As clsSearchResultsBaseClass
            objSearchResultBase = DirectCast(objSearchResult, clsSearchResultsBaseClass)

            MyBase.ComputePseudoPeptideLocInProtein(objSearchResultBase)

            With objSearchResult

                ' Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                ' If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide 
                ' will all be based on the first protein since Sequest only outputs the prefix and suffix letters for the first protein
                .ComputePeptideCleavageStateInProtein()

                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.DeltaCn2), .PeptideDeltaCn2)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.RankSP), .PeptideRankSP)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.RankXC), .PeptideRankXC)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.DelM), .PeptideDeltaMass)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.XcRatio), .PeptideXcRatio)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.PassFilt), .PeptidePassFilt)           ' Legacy/Unused
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.MScore), .PeptideMScore)               ' Legacy/Unused
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.NTT), .PeptideNTT)

                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.IonsObserved), .IonsObserved)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.IonsExpected), .IonsExpected)
                GetColumnValue(strSplitLine, intColumnMapping(eSequestSynopsisFileColumns.DelMPPM), .DelMPPM)

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

    ''' <summary>
    ''' Main processing function
    ''' </summary>
    ''' <param name="strInputFilePath">Sequest Synopsis or First-hits file</param>
    ''' <param name="strOutputFolderPath">Output folder</param>
    ''' <param name="strParameterFilePath">Parameter file</param>
    ''' <returns>True if success, False if failure</returns>
    Public Overloads Overrides Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean

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

            ' Set this to true since Sequest param files can have the same mod mass on different residues, and those residues may use different symbols
            mPeptideMods.ConsiderModSymbolWhenFindingIdenticalMods = True

            blnSuccess = ResetMassCorrectionTagsAndModificationDefinitions()
            If Not blnSuccess Then
                Return False
            End If

            MyBase.ResetProgress("Parsing " & Path.GetFileName(strInputFilePath))

            If Not CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                Return False
            End If

            ' Obtain the full path to the input file
            Dim fiInputFile = New FileInfo(strInputFilePath)

            Try
                blnSuccess = ParseSynopsisOrFirstHitsFile(fiInputFile.FullName, strOutputFolderPath, False)
            Catch ex As Exception
                SetErrorMessage("Error calling ParseSynopsisOrFirstHitsFile" & ex.Message)
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
                blnSuccess = False
            End Try

            If blnSuccess AndAlso mCreateProteinModsFile Then
                blnSuccess = CreateProteinModsFileWork(fiInputFile, strOutputFolderPath)
            End If

            If blnSuccess Then
                MyBase.OperationComplete()
            End If

        Catch ex As Exception
            SetErrorMessage("Error in clsSequestResultsProcessor.ProcessFile:  " & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.UnspecifiedError)
        End Try

        Return blnSuccess

    End Function

    Private Function CreateProteinModsFileWork(fiInputFile As FileInfo, strOutputFolderPath As String) As Boolean
        Dim blnSuccess As Boolean

        ' First create the MTS PepToProteinMap file using fiInputFile
        ' Will also look for the first hits or synopsis file and use that too if it is present

        Dim lstSourcePHRPDataFiles = New List(Of String)
        Dim strMTSPepToProteinMapFilePath As String = String.Empty

        Dim strAdditionalFile As String
        Dim strInputFileBaseName As String = Path.GetFileNameWithoutExtension(fiInputFile.Name)

        lstSourcePHRPDataFiles.Add(fiInputFile.FullName)
        If strInputFileBaseName.ToLower().EndsWith(FILENAME_SUFFIX_SYNOPSIS_FILE) Then
            strAdditionalFile = MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_FIRST_HITS_FILE)
            If File.Exists(strAdditionalFile) Then
                lstSourcePHRPDataFiles.Add(strAdditionalFile)
            End If

        ElseIf strInputFileBaseName.ToLower().EndsWith(FILENAME_SUFFIX_FIRST_HITS_FILE) Then
            strAdditionalFile = MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_SYNOPSIS_FILE)
            If File.Exists(strAdditionalFile) Then
                lstSourcePHRPDataFiles.Add(strAdditionalFile)
            End If
        End If

        strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(fiInputFile.FullName, strOutputFolderPath, MTS:=True)

        If File.Exists(strMTSPepToProteinMapFilePath) AndAlso mUseExistingMTSPepToProteinMapFile Then
            blnSuccess = True
        Else
            ' Mapping file does not exist
            blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
            If Not blnSuccess Then
                ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False")
            End If
        End If

        If blnSuccess Then
            ' If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
            MyBase.ValidatePHRPReaderSupportFiles(fiInputFile.FullName, strOutputFolderPath)

            ' Now create the Protein Mods file
            blnSuccess = MyBase.CreateProteinModDetailsFile(fiInputFile.FullName, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Sequest)
        End If

        If Not blnSuccess Then
            ' Do not treat this as a fatal error
            blnSuccess = True
        End If

        Return blnSuccess

    End Function

    Private Function ParseSequestSynFileHeaderLine(
      strLineIn As String,
      ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eSequestSynopsisFileColumns
        Dim lstColumnNames As SortedDictionary(Of String, eSequestSynopsisFileColumns)
        lstColumnNames = New SortedDictionary(Of String, eSequestSynopsisFileColumns)(StringComparer.InvariantCultureIgnoreCase)

        ReDim intColumnMapping(SequestSynopsisFileColCount - 1)

        lstColumnNames.Add("HitNum", eSequestSynopsisFileColumns.RowIndex)
        lstColumnNames.Add("ScanNum", eSequestSynopsisFileColumns.Scan)
        lstColumnNames.Add("ScanCount", eSequestSynopsisFileColumns.NumScans)
        lstColumnNames.Add("ChargeState", eSequestSynopsisFileColumns.Charge)
        lstColumnNames.Add("MH", eSequestSynopsisFileColumns.PeptideMH)
        lstColumnNames.Add("XCorr", eSequestSynopsisFileColumns.XCorr)
        lstColumnNames.Add("DelCn", eSequestSynopsisFileColumns.DeltaCn)
        lstColumnNames.Add("Sp", eSequestSynopsisFileColumns.Sp)
        lstColumnNames.Add("Reference", eSequestSynopsisFileColumns.ProteinName)
        lstColumnNames.Add("MultiProtein", eSequestSynopsisFileColumns.MultipleProteinCount)                     ' Multiple protein count: 0 if the peptide is in 1 protein, 1 if the peptide is in 2 proteins, etc.
        lstColumnNames.Add("Peptide", eSequestSynopsisFileColumns.PeptideSequence)
        lstColumnNames.Add("DelCn2", eSequestSynopsisFileColumns.DeltaCn2)
        lstColumnNames.Add("RankSP", eSequestSynopsisFileColumns.RankSP)
        lstColumnNames.Add("RankXC", eSequestSynopsisFileColumns.RankXC)
        lstColumnNames.Add("DelM", eSequestSynopsisFileColumns.DelM)
        lstColumnNames.Add("XcRatio", eSequestSynopsisFileColumns.XcRatio)
        lstColumnNames.Add("PassFilt", eSequestSynopsisFileColumns.PassFilt)                ' Legacy/unused
        lstColumnNames.Add("MScore", eSequestSynopsisFileColumns.MScore)                    ' Legacy/unused
        lstColumnNames.Add("NumTrypticEnds", eSequestSynopsisFileColumns.NTT)
        lstColumnNames.Add("Ions_Observed", eSequestSynopsisFileColumns.IonsObserved)
        lstColumnNames.Add("Ions_Expected", eSequestSynopsisFileColumns.IonsExpected)
        lstColumnNames.Add("DelM_PPM", eSequestSynopsisFileColumns.DelMPPM)

        ' The following columns are computed by this program and appended to the input file or saved in a new file
        lstColumnNames.Add("Cleavage_State", eSequestSynopsisFileColumns.Cleavage_State)
        lstColumnNames.Add("Terminus_State", eSequestSynopsisFileColumns.Terminus_State)
        lstColumnNames.Add("Mod_Count", eSequestSynopsisFileColumns.Mod_Count)
        lstColumnNames.Add("Mod_Description", eSequestSynopsisFileColumns.Mod_Description)
        lstColumnNames.Add("Monoisotopic_Mass", eSequestSynopsisFileColumns.Monoisotopic_Mass)

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
            SetErrorMessage("Error parsing header in Sequest synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

End Class
