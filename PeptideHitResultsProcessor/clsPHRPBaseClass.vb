Option Strict On

' This class can be used as a base class for peptide hit results processor classes
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
' Started January 6, 2006
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
Imports PHRPReader.clsPeptideCleavageStateCalculator
Imports System.IO
Imports System.Text.RegularExpressions

Public MustInherit Class clsPHRPBaseClass

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub New()
        mFileDate = "May 19, 2015"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
	Protected Const SEP_CHAR As Char = ControlChars.Tab
	Protected Const UNIQUE_SEQ_TO_PROTEIN_MAP_SEP As String = "_"c

	Protected Const COLUMN_NAME_UNIQUE_SEQ_ID As String = "Unique_Seq_ID"
	Protected Const COLUMN_NAME_PROTEIN_NAME As String = "Protein_Name"
	Protected Const COLUMN_NAME_RESULTID As String = "ResultID"
	Protected Const COLUMN_NAME_PEPTIDE As String = "Peptide"
	Protected Const COLUMN_NAME_RESIDUE As String = "Residue"
	Protected Const COLUMN_NAME_PROTEIN_RESIDUE_NUMBER As String = "Protein_Residue_Num"
	Protected Const COLUMN_NAME_RESIDUE_MOD_NAME As String = "Mod_Name"
	Protected Const COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER As String = "Peptide_Residue_Num"
	Protected Const COLUMN_NAME_MSGF_SPECPROB As String = "MSGF_SpecProb"

	Public Const XTANDEM_RESULTS_FILE_SUFFIX As String = "_xt.xml"
	Public Const SEQUEST_SYNOPSIS_FILE_SUFFIX As String = "_syn.txt"
	Public Const SEQUEST_FIRST_HITS_FILE_SUFFIX As String = "_fht.txt"

	Public Const INSPECT_RESULTS_FILE_SUFFIX As String = "_inspect.txt"
	Public Const INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX As String = "_fht.txt"
	Public Const INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX As String = "_Fscore_fht.txt"

	Public Const MSGFDB_RESULTS_FILE_SUFFIX As String = "_msgfdb.txt"

	Public Const MSALIGN_RESULTS_FILE_SUFFIX As String = "_MSAlign_ResultTable.txt"

    Public Const MODa_RESULTS_FILE_SUFFIX As String = "_moda.id.txt"
    Public Const MODPlus_RESULTS_FILE_SUFFIX As String = "_modp.id.txt"
    Public Const MSPathFinder_RESULTS_FILE_SUFFIX As String = "_IcTda.tsv"

	Public Const FILENAME_SUFFIX_RESULT_TO_SEQ_MAP As String = "_ResultToSeqMap.txt"
	Public Const FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP As String = "_SeqToProteinMap.txt"

	Public Const FILENAME_SUFFIX_SEQ_INFO As String = "_SeqInfo.txt"
	Public Const FILENAME_SUFFIX_MOD_DETAILS As String = "_ModDetails.txt"
	Public Const FILENAME_SUFFIX_MOD_SUMMARY As String = "_ModSummary.txt"

	Public Const FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING As String = "_PepToProtMap"
	Public Const FILENAME_SUFFIX_PROTEIN_MODS As String = "_ProteinMods.txt"
	Public Const FILENAME_SUFFIX_MSGF As String = "_MSGF.txt"

	Protected Const PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE As Single = 90
	Protected Const PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE As Single = 95

	Public Enum ePeptideHitResultsFileFormatConstants As Integer
		AutoDetermine = 0
		SequestSynopsisFile = 1
		SequestFirstHitsFile = 2
		XTandemXMLFile = 3
		InSpectTXTFile = 4
		MSGFDbTXTFile = 5
		MSAlignTXTFile = 6
        MODaTXTFile = 7
        MODPlusTXTFile = 8
        MSPathFinderTSVFile = 9
    End Enum

    Public Enum ePHRPErrorCodes
        NoError = 0
        InvalidInputFilePath = 1
        InvalidOutputFolderPath = 2
        ParameterFileNotFound = 3
        MassCorrectionTagsFileNotFound = 4
        ModificationDefinitionFileNotFound = 5

        ErrorReadingInputFile = 6
        ErrorCreatingOutputFiles = 7
        ErrorReadingParameterFile = 8
        ErrorReadingMassCorrectionTagsFile = 9
        ErrorReadingModificationDefinitionsFile = 10

        FilePathError = 11
        UnspecifiedError = -1
    End Enum
#End Region

#Region "Structures"
    Protected Structure udtSearchOptionModificationInfoType
        Dim SortOrder As Integer
        Dim ModificationMass As Double
        Dim TargetResidues As String
        Public ModificationType As clsModificationDefinition.eModificationTypeConstants
    End Structure

    Friend Structure udtModNameAndResidueLocType
        Dim ModName As String
        Dim ResidueLocInPeptide As Integer
    End Structure

    Protected Structure udtPepToProteinMappingType
        Public Peptide As String
        Public Protein As String
        Public ResidueStart As Integer
        Public ResidueEnd As Integer
    End Structure

#End Region

#Region "Classwide Variables"

    Protected mErrorCode As ePHRPErrorCodes = ePHRPErrorCodes.NoError
    Protected mErrorMessage As String = String.Empty
    Protected mWarnMissingParameterFileSection As Boolean

    Protected mFileDate As String = String.Empty
    Protected mAbortProcessing As Boolean

    Protected mCreateModificationSummaryFile As Boolean

    ' Note: If this is true and the _PepToProtMap.txt file isn't found then it will be created using the the Fasta file specified by mFastaFilePath
    Protected mCreateProteinModsFile As Boolean
    Protected mFastaFilePath As String
    Protected mIgnorePeptideToProteinMapperErrors As Boolean
    Protected mProteinModsFileIncludesReversedProteins As Boolean
    Protected mUseExistingMTSPepToProteinMapFile As Boolean

    ' The following two options are only used by the clsInSpecTResultsProcessor
    Protected mCreateInspectOrMSGFDbSynopsisFile As Boolean
    Protected mCreateInspectOrMSGFDbFirstHitsFile As Boolean

    Protected mMassCorrectionTagsFilePath As String
    Protected mModificationDefinitionsFilePath As String
    Protected mSearchToolParameterFilePath As String            ' Used by clsInSpecTResultsProcessor and clsMSGFDBResultsProcessor

    Protected mInspectSynopsisFilePValueThreshold As Single     ' Only used by clsInSpecTResultsProcessor; note that lower p-values are higher confidence results

    Protected mMODaMODPlusSynopsisFileProbabilityThreshold As Single   ' Used by clsMODaResultsProcessor and clsMODPlusResultsProcessor; note that higher probability are higher confidence results

    Protected mMSAlignSynopsisFilePValueThreshold As Single     ' Only used by clsMSAlignResultsProcessor; note that lower p-values are higher confidence results

    Protected mMSGFDBSynopsisFilePValueThreshold As Single      ' Used by clsMSGFDBResultsProcessor and clsMSPathFinderResultsProcessor; note that lower p-values are higher confidence results
    Protected mMSGFDBSynopsisFileSpecProbThreshold As Single    ' Used by clsMSGFDBResultsProcessor and clsMSPathFinderResultsProcessor; note that lower SpecProb values are higher confidence results

    Protected mEnzymeMatchSpec As udtEnzymeMatchSpecType
    Protected mPeptideNTerminusMassChange As Double             ' This is ignored if equal to 0; typical non-zero value is 1.0078246
    Protected mPeptideCTerminusMassChange As Double             ' This is ignored if equal to 0; typical non-zero value is 17.0027387

    Protected mPeptideSeqMassCalculator As clsPeptideMassCalculator

    Protected mPeptideMods As clsPeptideModificationContainer
    Protected mUniqueSequences As clsUniqueSequencesContainer
    Protected mSeqToProteinMap As Hashtable

    Protected mResultToSeqMapFile As StreamWriter
    Protected mSeqInfoFile As StreamWriter
    Protected mModDetailsFile As StreamWriter
    Protected mSeqToProteinMapFile As StreamWriter

    Protected mNextPeptideToProteinMapperLevel As Integer = 0

    Protected mReplaceSymbols As Regex
#End Region

#Region "Progress Events and Variables"
    Public Event ProgressReset()
    Public Event ProgressChanged(taskDescription As String, percentComplete As Single)     ' PercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Public Event ProgressComplete()

    Public Event ErrorOccurred(ErrMessage As String)
    Public Event WarningMessageEvent(WarningMessage As String)

    Protected mProgressStepDescription As String = String.Empty
    Protected mProgressPercentComplete As Single        ' Ranges from 0 to 100, but can contain decimal percentage values
#End Region

#Region "Properties"
    Public Property AbortProcessing() As Boolean
        Get
            Return mAbortProcessing
        End Get
        Set(Value As Boolean)
            mAbortProcessing = Value
        End Set
    End Property

    Public Property CreateModificationSummaryFile() As Boolean
        Get
            Return mCreateModificationSummaryFile
        End Get
        Set(Value As Boolean)
            mCreateModificationSummaryFile = Value
        End Set
    End Property

    Public Property CreateInspectFirstHitsFile() As Boolean
        Get
            Return mCreateInspectOrMSGFDbFirstHitsFile
        End Get
        Set(value As Boolean)
            mCreateInspectOrMSGFDbFirstHitsFile = value
        End Set
    End Property

    Public Property CreateInspectSynopsisFile() As Boolean
        Get
            Return mCreateInspectOrMSGFDbSynopsisFile
        End Get
        Set(value As Boolean)
            mCreateInspectOrMSGFDbSynopsisFile = value
        End Set
    End Property

    Public Property CreateMSGFDBFirstHitsFile() As Boolean
        Get
            Return mCreateInspectOrMSGFDbFirstHitsFile
        End Get
        Set(value As Boolean)
            mCreateInspectOrMSGFDbFirstHitsFile = value
        End Set
    End Property

    Public Property CreateMSGFDBSynopsisFile() As Boolean
        Get
            Return mCreateInspectOrMSGFDbSynopsisFile
        End Get
        Set(value As Boolean)
            mCreateInspectOrMSGFDbSynopsisFile = value
        End Set
    End Property

    Public Property CreateProteinModsFile() As Boolean
        Get
            Return mCreateProteinModsFile
        End Get
        Set(value As Boolean)
            mCreateProteinModsFile = value
        End Set
    End Property

    Public Property EnzymeMatchSpec() As udtEnzymeMatchSpecType
        Get
            Return mEnzymeMatchSpec
        End Get
        Set(Value As udtEnzymeMatchSpecType)
            mEnzymeMatchSpec = Value
        End Set
    End Property

    Public ReadOnly Property ErrorCode() As ePHRPErrorCodes
        Get
            Return mErrorCode
        End Get
    End Property

    Public ReadOnly Property ErrorMessage() As String
        Get
            Return GetErrorMessage()
        End Get
    End Property

    Public Property FastaFilePath() As String
        Get
            Return mFastaFilePath
        End Get
        Set(value As String)
            mFastaFilePath = value
        End Set
    End Property

    Public ReadOnly Property FileVersion() As String
        Get
            FileVersion = GetVersionForExecutingAssembly()
        End Get
    End Property

    Public ReadOnly Property FileDate() As String
        Get
            FileDate = mFileDate
        End Get
    End Property

    Public Property IgnorePeptideToProteinMapperErrors As Boolean
        Get
            Return mIgnorePeptideToProteinMapperErrors
        End Get
        Set(value As Boolean)
            mIgnorePeptideToProteinMapperErrors = value
        End Set
    End Property

    Public Property InspectSynopsisFilePValueThreshold() As Single
        Get
            Return mInspectSynopsisFilePValueThreshold
        End Get
        Set(value As Single)
            mInspectSynopsisFilePValueThreshold = value
        End Set
    End Property

    Public Property MassCorrectionTagsFilePath() As String
        Get
            Return mMassCorrectionTagsFilePath
        End Get
        Set(Value As String)
            mMassCorrectionTagsFilePath = Value
        End Set
    End Property

    Public Property ModificationDefinitionsFilePath() As String
        Get
            Return mModificationDefinitionsFilePath
        End Get
        Set(Value As String)
            mModificationDefinitionsFilePath = Value
        End Set
    End Property

    Public Property MODaMODPlusSynopsisFileProbabilityThreshold As Single
        Get
            Return mMODaMODPlusSynopsisFileProbabilityThreshold
        End Get
        Set(value As Single)
            mMODaMODPlusSynopsisFileProbabilityThreshold = value
        End Set
    End Property

    Public Property MSAlignSynopsisFilePValueThreshold() As Single
        Get
            Return mMSAlignSynopsisFilePValueThreshold
        End Get
        Set(value As Single)
            mMSAlignSynopsisFilePValueThreshold = value
        End Set
    End Property

    Public Property MSGFDBSynopsisFilePValueThreshold() As Single
        Get
            Return mMSGFDBSynopsisFilePValueThreshold
        End Get
        Set(value As Single)
            mMSGFDBSynopsisFilePValueThreshold = value
        End Set
    End Property

    Public Property MSGFDBSynopsisFileSpecProbThreshold As Single
        Get
            Return mMSGFDBSynopsisFileSpecProbThreshold
        End Get
        Set(value As Single)
            mMSGFDBSynopsisFileSpecProbThreshold = value
        End Set
    End Property

    Public Property PeptideCTerminusMassChange() As Double
        Get
            Return mPeptideCTerminusMassChange
        End Get
        Set(Value As Double)
            mPeptideCTerminusMassChange = Value
        End Set
    End Property

    Public Property PeptideNTerminusMassChange() As Double
        Get
            Return mPeptideNTerminusMassChange
        End Get
        Set(Value As Double)
            mPeptideNTerminusMassChange = Value
        End Set
    End Property

    Public Property ProteinModsFileIncludesReversedProteins() As Boolean
        Get
            Return mProteinModsFileIncludesReversedProteins
        End Get
        Set(value As Boolean)
            mProteinModsFileIncludesReversedProteins = value
        End Set
    End Property

    Public Overridable ReadOnly Property ProgressStepDescription() As String
        Get
            Return mProgressStepDescription
        End Get
    End Property

    ' ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Public ReadOnly Property ProgressPercentComplete() As Single
        Get
            Return CType(Math.Round(mProgressPercentComplete, 2), Single)
        End Get
    End Property

    Public Property SearchToolParameterFilePath() As String
        Get
            Return mSearchToolParameterFilePath
        End Get
        Set(value As String)
            mSearchToolParameterFilePath = value
        End Set
    End Property

    Public Property UseExistingMTSPepToProteinMapFile As Boolean
        Get
            Return mUseExistingMTSPepToProteinMapFile
        End Get
        Set(value As Boolean)
            mUseExistingMTSPepToProteinMapFile = False
        End Set
    End Property

    Public Property WarnMissingParameterFileSection() As Boolean
        Get
            Return mWarnMissingParameterFileSection
        End Get
        Set(Value As Boolean)
            mWarnMissingParameterFileSection = Value
        End Set
    End Property
#End Region

    Public Sub AbortProcessingNow()
        mAbortProcessing = True
    End Sub

    Public Shared Function AutoDefinePeptideHitResultsFilePath(ePeptideHitResultFileFormat As ePeptideHitResultsFileFormatConstants,
      strSourceFolderPath As String,
      strBaseName As String) As String

        If Not strBaseName Is Nothing AndAlso strBaseName.Length > 0 Then
            Select Case ePeptideHitResultFileFormat
                Case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & SEQUEST_FIRST_HITS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.XTandemXMLFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & XTANDEM_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.InSpectTXTFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & INSPECT_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & MSGFDB_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.MSAlignTXTFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & MSALIGN_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.MODaTXTFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & MODa_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.MODPlusTXTFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & MODPlus_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile
                    Return Path.Combine(strSourceFolderPath, strBaseName & MSPathFinder_RESULTS_FILE_SUFFIX)

                Case Else
                    ' Includes ePeptideHitResultsFileFormatConstants.AutoDetermine
                    ' Call AutoDefinePeptideHitResultsFilePath below sending it only strSourceFolderPath
            End Select
        End If

        Return AutoDefinePeptideHitResultsFilePath(strSourceFolderPath)
    End Function

    Public Shared Function AutoDefinePeptideHitResultsFilePath(strSourceFolderPath As String) As String
        ' Looks for a file ending in _syn.txt, _fht.txt, _xt.xml, or _inspect.txt in folder strSourceFolderPath
        ' Returns the first matching file found

        Dim ioFolderInfo As DirectoryInfo
        Dim ioFileInfo As FileInfo
        Dim intIndex As Integer

        Dim strMatchSpec As String = String.Empty

        Try
            For intIndex = 0 To 3
                Select Case intIndex
                    Case 0
                        strMatchSpec = "*" & SEQUEST_SYNOPSIS_FILE_SUFFIX
                    Case 1
                        strMatchSpec = "*" & SEQUEST_FIRST_HITS_FILE_SUFFIX
                    Case 2
                        strMatchSpec = "*" & XTANDEM_RESULTS_FILE_SUFFIX
                    Case 3
                        strMatchSpec = "*" & INSPECT_RESULTS_FILE_SUFFIX
                End Select

                ioFolderInfo = New DirectoryInfo(strSourceFolderPath)
                For Each ioFileInfo In ioFolderInfo.GetFiles(strMatchSpec)
                    ' If we get here, a match was found; return its path
                    Return ioFileInfo.FullName
                Next ioFileInfo
            Next intIndex
        Catch ex As Exception
            ' Ignore errors here
        End Try

        ' No match; return empty
        Return String.Empty
    End Function

    Protected Function CheckSeqToProteinMapDefined(intUniqueSeqID As Integer, strProteinName As String) As Boolean
        ' Returns True if the sequence to protein map was already defined
        ' Returns False if the mapping was not defined (will also update mSeqToProteinMap)

        Dim strKey As String
        Dim blnExistingMapFound As Boolean

        blnExistingMapFound = False

        Try
            If strProteinName Is Nothing Then strProteinName = String.Empty

            strKey = intUniqueSeqID.ToString & UNIQUE_SEQ_TO_PROTEIN_MAP_SEP & strProteinName

            If mSeqToProteinMap.ContainsKey(strKey) Then
                blnExistingMapFound = True
            Else
                mSeqToProteinMap.Add(strKey, 1)
                blnExistingMapFound = False
            End If
        Catch ex As Exception
            blnExistingMapFound = False
        End Try

        Return blnExistingMapFound
    End Function

    Protected Function CIntSafe(strValue As String, intDefaultValue As Integer) As Integer
        Try
            ' Note: Integer.Parse() fails if strValue contains a decimal point, even if it is "8.000"
            ' Thus, we're using CInt() instead
            Return CInt(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return intDefaultValue
        End Try
    End Function

    Protected Function CDblSafe(strValue As String, dblDefaultValue As Double) As Double
        Try
            Return Double.Parse(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return dblDefaultValue
        End Try
    End Function

    Protected Function CSngSafe(strValue As String, sngDefaultValue As Single) As Single
        Try
            Return Single.Parse(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return sngDefaultValue
        End Try
    End Function

    Protected Function CleanupFilePaths(ByRef strInputFilePath As String, ByRef strOutputFolderPath As String) As Boolean
        ' Returns True if success, False if failure

        Dim ioFileInfo As FileInfo
        Dim ioFolder As DirectoryInfo

        Try
            ' Make sure strInputFilePath points to a valid file
            ioFileInfo = New FileInfo(strInputFilePath)

            If Not ioFileInfo.Exists() Then
                SetErrorMessage("Input file not found: " & strInputFilePath)
                SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath)
                Return False
            Else
                If String.IsNullOrWhiteSpace(strOutputFolderPath) Then
                    ' Define strOutputFolderPath based on strInputFilePath
                    strOutputFolderPath = ioFileInfo.DirectoryName
                End If

                ' Make sure strOutputFolderPath points to a folder
                ioFolder = New DirectoryInfo(strOutputFolderPath)

                If Not ioFolder.Exists() Then
                    ' strOutputFolderPath points to a non-existent folder; attempt to create it
                    Try
                        ioFolder.Create()
                    Catch ex As Exception
                        SetErrorMessage("Invalid output folder: " & strOutputFolderPath)
                        SetErrorCode(ePHRPErrorCodes.InvalidOutputFolderPath)
                        Return False
                    End Try
                End If

                Return True
            End If

        Catch ex As Exception
            SetErrorMessage("Error cleaning up the file paths: " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.FilePathError)
            Return False
        End Try

    End Function

    ''' <summary>
    ''' Collapses a list of strings to a tab-delimited line of text
    ''' </summary>
    ''' <param name="lstFields"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CollapseList(lstFields As List(Of String)) As String

        Return String.Join(ControlChars.Tab, lstFields)

    End Function

    Public Shared Function DetermineResultsFileFormat(strFilePath As String) As ePeptideHitResultsFileFormatConstants
        ' Examine the extension on strFilePath to determine the file format

        Dim strExtensionLCase As String
        Dim strBaseFileNameLCase As String

        strExtensionLCase = Path.GetExtension(strFilePath).ToLower()
        strBaseFileNameLCase = Path.GetFileNameWithoutExtension(strFilePath).ToLower()

        If strExtensionLCase = ".xml" Then
            Return ePeptideHitResultsFileFormatConstants.XTandemXMLFile

        ElseIf strBaseFileNameLCase.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile

        ElseIf strBaseFileNameLCase.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.SequestSynopsisFile

        ElseIf strBaseFileNameLCase.EndsWith(clsInSpecTResultsProcessor.FILENAME_SUFFIX_INSPECT_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.InSpectTXTFile

        ElseIf strBaseFileNameLCase.EndsWith(clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile

        ElseIf strBaseFileNameLCase.EndsWith(clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile

        ElseIf strBaseFileNameLCase.EndsWith(clsMSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.MSAlignTXTFile

        ElseIf strBaseFileNameLCase.EndsWith(clsMODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.MODaTXTFile

        ElseIf strBaseFileNameLCase.EndsWith(clsMODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.MODPlusTXTFile

        ElseIf strBaseFileNameLCase.EndsWith(clsMSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_File.ToLower()) Then
            Return ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile

        ElseIf strExtensionLCase = ".tsv" Then
            ' Assume this is an MSGF+ TSV file
            Return ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile

        Else
            ' Unknown extension
            Return ePeptideHitResultsFileFormatConstants.AutoDetermine
        End If
    End Function

    Protected Sub CloseSequenceOutputFiles()
        Try
            If Not mResultToSeqMapFile Is Nothing Then
                mResultToSeqMapFile.Close()
                mResultToSeqMapFile = Nothing
            End If
        Catch ex As Exception
        End Try

        Try
            If Not mSeqInfoFile Is Nothing Then
                mSeqInfoFile.Close()
                mSeqInfoFile = Nothing
            End If
        Catch ex As Exception
        End Try

        Try
            If Not mModDetailsFile Is Nothing Then
                mModDetailsFile.Close()
                mModDetailsFile = Nothing
            End If
        Catch ex As Exception
        End Try

        Try
            If Not mSeqToProteinMapFile Is Nothing Then
                mSeqToProteinMapFile.Close()
                mSeqToProteinMapFile = Nothing
            End If
        Catch ex As Exception
        End Try
    End Sub

    Protected Sub ComputePseudoPeptideLocInProtein(objSearchResult As clsSearchResultsBaseClass)

        With objSearchResult

            ' Set these to 1 and 10000 since MSGFDB, Sequest, and Inspect results files do not contain protein sequence information
            ' If we find later that the peptide sequence spans the length of the protein, we'll revise .ProteinSeqResidueNumberEnd as needed
            .ProteinSeqResidueNumberStart = 1
            .ProteinSeqResidueNumberEnd = 10000

            If .PeptidePreResidues.Trim.EndsWith(TERMINUS_SYMBOL_SEQUEST) Then
                ' The peptide is at the N-Terminus of the protein
                .PeptideLocInProteinStart = .ProteinSeqResidueNumberStart
                .PeptideLocInProteinEnd = .PeptideLocInProteinStart + .PeptideCleanSequence.Length - 1

                If .PeptidePostResidues.Trim.Chars(0) = TERMINUS_SYMBOL_SEQUEST Then
                    ' The peptide spans the entire length of the protein
                    .ProteinSeqResidueNumberEnd = .PeptideLocInProteinEnd
                Else
                    If .PeptideLocInProteinEnd > .ProteinSeqResidueNumberEnd Then
                        ' The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                        .ProteinSeqResidueNumberEnd = .PeptideLocInProteinEnd + 1
                    End If
                End If
            ElseIf .PeptidePostResidues.Trim.StartsWith(TERMINUS_SYMBOL_SEQUEST) Then
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

        End With

    End Sub

    Protected Overridable Function ConstructPepToProteinMapFilePath(strInputFilePath As String, strOutputFolderPath As String, MTS As Boolean) As String

        Dim strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath)

        If MTS Then
            strPepToProteinMapFilePath &= FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING & "MTS.txt"
        Else
            strPepToProteinMapFilePath &= FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING & ".txt"
        End If

        If String.IsNullOrWhiteSpace(strOutputFolderPath) Then
            Dim fiInputFile = New FileInfo(strInputFilePath)
            strOutputFolderPath = fiInputFile.DirectoryName
        End If

        strPepToProteinMapFilePath = Path.Combine(strOutputFolderPath, strPepToProteinMapFilePath)

        Return strPepToProteinMapFilePath

    End Function

    ''' <summary>
    ''' Use the PeptideToProteinMapEngine to create the Peptide to Protein map file for the file or files in lstSourcePHRPDataFiles
    ''' </summary>
    ''' <param name="lstSourcePHRPDataFiles"></param>
    ''' <param name="strMTSPepToProteinMapFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CreatePepToProteinMapFile(lstSourcePHRPDataFiles As List(Of String), strMTSPepToProteinMapFilePath As String) As Boolean

        Dim objPeptideToProteinMapper As PeptideToProteinMapEngine.clsPeptideToProteinMapEngine

        Dim strLineIn As String
        Dim strSplitLine() As String
        Dim strPeptideAndProteinKey As String

        Dim htPeptideToProteinMapResults As Collections.Hashtable

        Dim blnSuccess As Boolean = False

        Try

            If String.IsNullOrEmpty(strMTSPepToProteinMapFilePath) Then
                SetErrorMessage("Cannot create the PepToProtein map file because strMTSPepToProteinMapFilePath is empty; likely a programming bug")
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
                Return False
            End If

            If String.IsNullOrEmpty(mFastaFilePath) Then
                SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File Path is not defined")
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
                Return False
            ElseIf Not IO.File.Exists(mFastaFilePath) Then
                SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File was not found: " & mFastaFilePath)
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
                Return False
            End If

            ' Verify that the fasta file is not a DNA-sequence based fasta file
            blnSuccess = ValidateProteinFastaFile(mFastaFilePath)
            If Not blnSuccess Then
                Return False
            End If

            Console.WriteLine()
            Console.WriteLine()
            UpdateProgress("Creating Peptide to Protein Map file", PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE)

            ' Initialize items
            Dim fiMTSPepToProteinMapFile As FileInfo = New FileInfo(strMTSPepToProteinMapFilePath)
            Dim strOutputFolderPath As String = fiMTSPepToProteinMapFile.DirectoryName

            fiMTSPepToProteinMapFile = New FileInfo(strMTSPepToProteinMapFilePath)
            strOutputFolderPath = fiMTSPepToProteinMapFile.DirectoryName
            htPeptideToProteinMapResults = New Collections.Hashtable

            objPeptideToProteinMapper = New PeptideToProteinMapEngine.clsPeptideToProteinMapEngine

            With objPeptideToProteinMapper
                .DeleteTempFiles = True
                .IgnoreILDifferences = False
                .InspectParameterFilePath = String.Empty
                .LogMessagesToFile = False
                .MatchPeptidePrefixAndSuffixToProtein = False
                .OutputProteinSequence = False
                .PeptideInputFileFormat = PeptideToProteinMapEngine.clsPeptideToProteinMapEngine.ePeptideInputFileFormatConstants.PHRPFile
                .PeptideFileSkipFirstLine = False
                .ProteinDataRemoveSymbolCharacters = True
                .ProteinInputFilePath = mFastaFilePath
                .SaveProteinToPeptideMappingFile = True
                .SearchAllProteinsForPeptideSequence = True
                .SearchAllProteinsSkipCoverageComputationSteps = True
                .ShowMessages = True
            End With

            Using swMTSpepToProteinMapFile = New StreamWriter(New FileStream(strMTSPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                For Each strInputFilePath As String In lstSourcePHRPDataFiles

                    Dim strResultsFilePath As String
                    strResultsFilePath = Path.GetFileNameWithoutExtension(strInputFilePath) & PeptideToProteinMapEngine.clsPeptideToProteinMapEngine.FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING
                    strResultsFilePath = Path.Combine(strOutputFolderPath, strResultsFilePath)

                    ' Make sure the results file doesn't already exist
                    DeleteFileIgnoreErrors(strResultsFilePath)

                    AddHandler objPeptideToProteinMapper.ProgressChanged, AddressOf PeptideToProteinMapper_ProgressChanged
                    mNextPeptideToProteinMapperLevel = 25

                    blnSuccess = objPeptideToProteinMapper.ProcessFile(strInputFilePath, strOutputFolderPath, String.Empty, True)

                    RemoveHandler objPeptideToProteinMapper.ProgressChanged, AddressOf PeptideToProteinMapper_ProgressChanged

                    If blnSuccess Then
                        If Not File.Exists(strResultsFilePath) Then
                            SetErrorMessage("Peptide to protein mapping file was not created for " & strInputFilePath)
                            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                            blnSuccess = False
                            Exit For
                        Else
                            blnSuccess = ValidatePeptideToProteinMapResults(strResultsFilePath, mIgnorePeptideToProteinMapperErrors)
                        End If
                    Else
                        If String.IsNullOrWhiteSpace(objPeptideToProteinMapper.GetErrorMessage) AndAlso objPeptideToProteinMapper.StatusMessage.ToLower().Contains("error") Then
                            SetErrorMessage("Error running clsPeptideToProteinMapEngine: " & objPeptideToProteinMapper.StatusMessage)
                        Else
                            If objPeptideToProteinMapper.StatusMessage.Length > 0 Then
                                SetErrorMessage("clsPeptideToProteinMapEngine status: " & objPeptideToProteinMapper.StatusMessage)
                            End If
                            SetErrorMessage("Error running clsPeptideToProteinMapEngine: " & objPeptideToProteinMapper.GetErrorMessage())
                        End If

                        If mIgnorePeptideToProteinMapperErrors Then

                            ReportWarning("Ignoring protein mapping error since 'IgnorePeptideToProteinMapperErrors' = True")

                            If File.Exists(strResultsFilePath) Then
                                blnSuccess = ValidatePeptideToProteinMapResults(strResultsFilePath, mIgnorePeptideToProteinMapperErrors)
                            Else
                                mErrorMessage = String.Empty
                                mErrorCode = ePHRPErrorCodes.NoError
                                blnSuccess = True
                            End If

                        Else
                            SetErrorMessage("Error in CreatePepToProteinMapFile")
                            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                            blnSuccess = False
                        End If
                    End If


                    If Not File.Exists(strResultsFilePath) Then
                        Continue For
                    End If

                    ' Read the newly created file and append new entries to strMTSPepToProteinMapFilePath
                    Using srResultsFile = New StreamReader(New FileStream(strResultsFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))
                        While Not srResultsFile.EndOfStream
                            strLineIn = srResultsFile.ReadLine()

                            If Not String.IsNullOrWhiteSpace(strLineIn) Then
                                strSplitLine = Split(strLineIn, ControlChars.Tab, 2)
                                If strSplitLine.Length >= 2 Then
                                    strPeptideAndProteinKey = strSplitLine(0) & "_" & strSplitLine(1)

                                    If Not htPeptideToProteinMapResults.ContainsKey(strPeptideAndProteinKey) Then
                                        htPeptideToProteinMapResults.Add(strPeptideAndProteinKey, 0)
                                        swMTSpepToProteinMapFile.WriteLine(strLineIn)
                                    End If
                                End If

                            End If
                        End While

                    End Using

                    ' Delete the interim results file
                    DeleteFileIgnoreErrors(strResultsFilePath)


                Next

            End Using

            objPeptideToProteinMapper.CloseLogFileNow()


        Catch ex As Exception
            SetErrorMessage("Error in CreatePepToProteinMapFile:" & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Create the protein mod details file for the specified PHRP data file
    ''' </summary>
    ''' <param name="strPHRPDataFilePath"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function CreateProteinModDetailsFile(strPHRPDataFilePath As String, strOutputFolderPath As String) As Boolean

        Dim strMTSPepToProteinMapFilePath As String

        Dim blnSuccess As Boolean = False

        Try
            Dim fiInputFile = New FileInfo(strPHRPDataFilePath)

            Dim lstSourcePHRPDataFiles As List(Of String) = New List(Of String)
            lstSourcePHRPDataFiles.Add(fiInputFile.FullName)

            strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(fiInputFile.FullName, strOutputFolderPath, MTS:=True)

            blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)

            If blnSuccess Then
                blnSuccess = CreateProteinModDetailsFile(strPHRPDataFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Unknown)
            End If

        Catch ex As Exception
            ReportWarning("Error in CreateProteinModDetailsFile: " & ex.Message)
        End Try

        Return blnSuccess

    End Function

    Public Function CreateProteinModDetailsFile(
       strPHRPDataFilePath As String,
       strOutputFolderPath As String,
       strMTSPepToProteinMapFilePath As String,
       ePHRPResultType As PHRPReader.clsPHRPReader.ePeptideHitResultType) As Boolean

        Dim intPepToProteinMapIndex As Integer

        Dim strHeaderLine As String = String.Empty

        Dim fiPHRPDataFile As FileInfo
        Dim strProteinModsFilePath As String

        Dim strResidue As String
        Dim strCleanSequence As String
        Dim intResidueLocInProtein As Integer

        Dim intPSMCount As Integer
        Dim intPSMCountSkippedSinceReversedOrScrambledProtein As Integer
        Dim sngProgressAtStart As Single

        Dim blnSkipProtein As Boolean
        Dim blnSuccess As Boolean = False

        Try
            Console.WriteLine()

            sngProgressAtStart = mProgressPercentComplete
            UpdateProgress("Creating the Protein Mod Details file", PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE)

            ' Confirm that the PHRP data file exists
            fiPHRPDataFile = New FileInfo(strPHRPDataFilePath)
            If Not fiPHRPDataFile.Exists() Then
                SetErrorMessage("PHRP data file not found in CreateProteinModDetailsFile: " & strPHRPDataFilePath)
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
                Return False
            End If

            ' Confirm that the _PepToProtMapMTS.txt file exists
            If String.IsNullOrEmpty(strMTSPepToProteinMapFilePath) Then
                SetErrorMessage("Cannot create the ProteinMods file because strMTSPepToProteinMapFilePath is empty")
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
                Return False
            End If

            ' Initialize lstPepToProteinMapping
            Dim lstPepToProteinMapping = New List(Of udtPepToProteinMappingType)

            ' Read the _PepToProtMapMTS file
            blnSuccess = LoadPeptideToProteinMapInfo(strMTSPepToProteinMapFilePath, lstPepToProteinMapping, strHeaderLine)
            If Not blnSuccess Then
                Return False
            End If

            ' Assure that lstPepToProteinMapping is sorted on peptide
            If lstPepToProteinMapping.Count > 1 Then
                lstPepToProteinMapping.Sort(New PepToProteinMappingComparer)
            End If

            strProteinModsFilePath = ReplaceFilenameSuffix(fiPHRPDataFile, FILENAME_SUFFIX_PROTEIN_MODS)
            If Not String.IsNullOrEmpty(strOutputFolderPath) Then
                strProteinModsFilePath = IO.Path.Combine(strOutputFolderPath, IO.Path.GetFileName(strProteinModsFilePath))
            End If

            intPSMCount = 0
            intPSMCountSkippedSinceReversedOrScrambledProtein = 0

            ' Create a ProteinMods file parallel to the PHRP file
            Using swProteinModsFile = New StreamWriter(New FileStream(strProteinModsFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                ' Write the header line
                swProteinModsFile.WriteLine(
                   COLUMN_NAME_RESULTID & ControlChars.Tab &
                   COLUMN_NAME_PEPTIDE & ControlChars.Tab &
                   COLUMN_NAME_UNIQUE_SEQ_ID & ControlChars.Tab &
                   COLUMN_NAME_PROTEIN_NAME & ControlChars.Tab &
                   COLUMN_NAME_RESIDUE & ControlChars.Tab &
                   COLUMN_NAME_PROTEIN_RESIDUE_NUMBER & ControlChars.Tab &
                   COLUMN_NAME_RESIDUE_MOD_NAME & ControlChars.Tab &
                   COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER & ControlChars.Tab &
                   COLUMN_NAME_MSGF_SPECPROB)

                Dim blnLoadMSGFResults = ePHRPResultType <> clsPHRPReader.ePeptideHitResultType.MSGFDB

                Using objReader As New clsPHRPReader(strPHRPDataFilePath, ePHRPResultType, blnLoadModsAndSeqInfo:=True, blnLoadMSGFResults:=blnLoadMSGFResults, blnLoadScanStats:=False)
                    objReader.EchoMessagesToConsole = True
                    objReader.SkipDuplicatePSMs = True

                    For Each strErrorMessage As String In objReader.ErrorMessages
                        RaiseEvent ErrorOccurred(strErrorMessage)
                    Next

                    For Each strWarningMessage As String In objReader.WarningMessages
                        If strWarningMessage.StartsWith("MSGF file not found") Then
                            strWarningMessage = "MSGF file not found; column " & COLUMN_NAME_MSGF_SPECPROB & " will not have any data"
                        End If
                        ReportWarning(strWarningMessage)
                    Next

                    objReader.ClearErrors()
                    objReader.ClearWarnings()

                    AddHandler objReader.ErrorEvent, AddressOf PHRPReader_ErrorEvent
                    AddHandler objReader.WarningEvent, AddressOf PHRPReader_WarningEvent

                    Do While objReader.MoveNext()

                        ' Use binary search to find this peptide in lstPepToProteinMapping
                        intPepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(lstPepToProteinMapping, objReader.CurrentPSM.Peptide)

                        If intPepToProteinMapIndex >= 0 Then

                            Do
                                intPSMCount += 1

                                blnSkipProtein = False
                                If Not mProteinModsFileIncludesReversedProteins Then
                                    blnSkipProtein = IsReversedProtein(lstPepToProteinMapping(intPepToProteinMapIndex).Protein)
                                    If blnSkipProtein Then
                                        intPSMCountSkippedSinceReversedOrScrambledProtein += 1
                                    End If
                                End If

                                If Not blnSkipProtein Then

                                    For Each objMod As PHRPReader.clsAminoAcidModInfo In objReader.CurrentPSM.ModifiedResidues

                                        intResidueLocInProtein = lstPepToProteinMapping(intPepToProteinMapIndex).ResidueStart + objMod.ResidueLocInPeptide - 1

                                        If IsLetterAtoZ(objMod.Residue) Then
                                            strResidue = objMod.Residue
                                        Else
                                            strCleanSequence = objReader.CurrentPSM.PeptideCleanSequence

                                            If objMod.ResidueLocInPeptide < 1 Then
                                                ' This shouldn't be the case, but we'll check for it anyway
                                                strResidue = strCleanSequence.Substring(0, 1)

                                            ElseIf objMod.ResidueLocInPeptide > strCleanSequence.Length Then
                                                ' This shouldn't be the case, but we'll check for it anyway
                                                strResidue = strCleanSequence.Substring(strCleanSequence.Length - 1, 1)

                                            Else
                                                strResidue = strCleanSequence.Substring(objMod.ResidueLocInPeptide - 1, 1)
                                            End If

                                        End If

                                        If lstPepToProteinMapping(intPepToProteinMapIndex).Protein = "__NoMatch__" AndAlso IsReversedProtein(objReader.CurrentPSM.ProteinFirst) Then
                                            ' Skip this result
                                            intPSMCountSkippedSinceReversedOrScrambledProtein += 1
                                        Else
                                            swProteinModsFile.WriteLine(
                                              objReader.CurrentPSM.ResultID & ControlChars.Tab &
                                              objReader.CurrentPSM.Peptide & ControlChars.Tab &
                                              objReader.CurrentPSM.SeqID & ControlChars.Tab &
                                              lstPepToProteinMapping(intPepToProteinMapIndex).Protein & ControlChars.Tab &
                                              strResidue & ControlChars.Tab &
                                              intResidueLocInProtein.ToString() & ControlChars.Tab &
                                              objMod.ModDefinition.MassCorrectionTag & ControlChars.Tab &
                                              objMod.ResidueLocInPeptide.ToString() & ControlChars.Tab &
                                              objReader.CurrentPSM.MSGFSpecProb)
                                        End If

                                    Next
                                End If

                                intPepToProteinMapIndex += 1
                            Loop While intPepToProteinMapIndex < lstPepToProteinMapping.Count AndAlso objReader.CurrentPSM.Peptide = lstPepToProteinMapping(intPepToProteinMapIndex).Peptide

                        Else
                            ReportWarning("Peptide not found in lstPepToProteinMapping: " & objReader.CurrentPSM.Peptide)
                        End If

                        UpdateProgress(sngProgressAtStart + objReader.PercentComplete * (100 - sngProgressAtStart) / 100)
                    Loop

                    RemoveHandler objReader.ErrorEvent, AddressOf PHRPReader_ErrorEvent
                    RemoveHandler objReader.WarningEvent, AddressOf PHRPReader_WarningEvent

                    If intPSMCount > 0 Then
                        If intPSMCountSkippedSinceReversedOrScrambledProtein = intPSMCount Then
                            Console.WriteLine()
                            ReportWarning("All PSMs map to reversed or scrambled proteins; the _ProteinMods.txt file is empty")
                        ElseIf intPSMCountSkippedSinceReversedOrScrambledProtein > 0 Then
                            Console.WriteLine()
                            Console.WriteLine("Note: skipped " & intPSMCountSkippedSinceReversedOrScrambledProtein & " / " & intPSMCount & " PSMs that map to reversed or scrambled proteins while creating the _ProteinMods.txt file")
                        End If
                    End If

                End Using
            End Using

            blnSuccess = True

        Catch ex As Exception
            SetErrorMessage("Error in CreateProteinModDetailsFile:" & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
        End Try

        Return blnSuccess

    End Function

    Protected Sub DeleteFileIgnoreErrors(strFilePath As String)

        Try
            If IO.File.Exists(strFilePath) Then
                Threading.Thread.Sleep(200)
                File.Delete(strFilePath)
            End If

        Catch ex As Exception
            ' Ignore errors here
        End Try
    End Sub

    Protected Sub ExpandListIfRequired(Of T)(lstItems As List(Of T), countToAdd As Integer, Optional largeListThreshold As Integer = 1000000)

        If lstItems.Count > largeListThreshold AndAlso lstItems.Count + countToAdd > lstItems.Capacity Then
            ' .NET by default will double the size of the list to accomodate these new items
            ' Instead, expand the list by 20% of the current size
            lstItems.Capacity = lstItems.Capacity + CInt(lstItems.Count / 5)
        End If
    End Sub

    Protected Function FindFirstMatchInPepToProteinMapping(ByRef lstPepToProteinMapping As List(Of udtPepToProteinMappingType), strPeptideToFind As String) As Integer

        Static objPeptideSearchComparer As PepToProteinMappingPeptideSearchComparer = New PepToProteinMappingPeptideSearchComparer()

        Dim intPepToProteinMapIndex As Integer

        ' Use binary search to find this peptide in lstPepToProteinMapping
        Dim udtItemToFind As udtPepToProteinMappingType = New udtPepToProteinMappingType
        udtItemToFind.Peptide = strPeptideToFind

        intPepToProteinMapIndex = lstPepToProteinMapping.BinarySearch(udtItemToFind, objPeptideSearchComparer)

        If intPepToProteinMapIndex > 0 Then
            ' Step Backward until the first match is found
            Do While intPepToProteinMapIndex > 0 AndAlso lstPepToProteinMapping(intPepToProteinMapIndex - 1).Peptide = strPeptideToFind
                intPepToProteinMapIndex -= 1
            Loop
        End If

        Return intPepToProteinMapIndex

    End Function

    Protected Function GetCleanSequence(strSequenceWithMods As String) As String
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Return GetCleanSequence(strSequenceWithMods, strPrefix, strSuffix)
    End Function

    Protected Function GetCleanSequence(strSequenceWithMods As String, ByRef strPrefix As String, ByRef strSuffix As String) As String

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

    ''' <summary>
    ''' If intColumnIndex is >= 0 then updates strValue with the value at strSplitLine(intColumnIndex)
    ''' Otherwise, updates strValue to String.Empty
    ''' </summary>
    ''' <returns>True if intColumnIndex >= 0</returns>
    ''' <remarks></remarks>
    Protected Function GetColumnValue(ByRef strSplitLine() As String, intColumnIndex As Integer, ByRef strValue As String) As Boolean
        Return GetColumnValue(strSplitLine, intColumnIndex, strValue, String.Empty)
    End Function

    ''' <summary>
    ''' If intColumnIndex is >= 0 then updates intValue with the value at strSplitLine(intColumnIndex)
    ''' Otherwise, updates intValue to 0
    ''' </summary>
    ''' <returns>True if intColumnIndex >= 0</returns>
    ''' <remarks></remarks>
    Protected Function GetColumnValue(ByRef strSplitLine() As String, intColumnIndex As Integer, ByRef intValue As Integer) As Boolean
        Return GetColumnValue(strSplitLine, intColumnIndex, intValue, 0)
    End Function

    ''' <summary>
    ''' If intColumnIndex is >= 0 then updates strValue with the value at strSplitLine(intColumnIndex)
    ''' Otherwise, updates strValue to strValueIfMissing
    ''' </summary>
    ''' <returns>True if intColumnIndex >= 0</returns>
    ''' <remarks></remarks>
    Protected Function GetColumnValue(ByRef strSplitLine() As String, intColumnIndex As Integer, ByRef strValue As String, strValueIfMissing As String) As Boolean
        If intColumnIndex >= 0 AndAlso intColumnIndex < strSplitLine.Length Then
            strValue = String.Copy(strSplitLine(intColumnIndex))
            Return True
        Else
            strValue = String.Copy(strValueIfMissing)
            Return False
        End If
    End Function

    ''' <summary>
    ''' If intColumnIndex is >= 0 then updates intValue with the value at strSplitLine(intColumnIndex)
    ''' Otherwise, updates strValue to intValueIfMissing
    ''' </summary>
    ''' <returns>True if intColumnIndex >= 0</returns>
    ''' <remarks></remarks>
    Protected Function GetColumnValue(ByRef strSplitLine() As String, intColumnIndex As Integer, ByRef intValue As Integer, intValueIfMissing As Integer) As Boolean
        Dim strValue As String = String.Empty

        If GetColumnValue(strSplitLine, intColumnIndex, strValue, intValueIfMissing.ToString) Then
            If Integer.TryParse(strValue, intValue) Then
                Return True
            Else
                intValue = intValueIfMissing
                Return False
            End If
        Else
            intValue = intValueIfMissing
            Return False
        End If
    End Function

    Protected Function GetErrorMessage() As String
        ' Returns String.Empty if no error

        Dim strMessage As String

        Select Case Me.ErrorCode
            Case ePHRPErrorCodes.NoError
                strMessage = String.Empty
            Case ePHRPErrorCodes.InvalidInputFilePath
                strMessage = "Invalid input file path"
            Case ePHRPErrorCodes.InvalidOutputFolderPath
                strMessage = "Invalid output folder path"
            Case ePHRPErrorCodes.ParameterFileNotFound
                strMessage = "Parameter file not found"
            Case ePHRPErrorCodes.MassCorrectionTagsFileNotFound
                strMessage = "Mass correction tags file not found"
            Case ePHRPErrorCodes.ModificationDefinitionFileNotFound
                strMessage = "Modification definition file not found"

            Case ePHRPErrorCodes.ErrorReadingInputFile
                strMessage = "Error reading input file"
            Case ePHRPErrorCodes.ErrorCreatingOutputFiles
                strMessage = "Error creating output files"
            Case ePHRPErrorCodes.ErrorReadingParameterFile
                strMessage = "Invalid parameter file"
            Case ePHRPErrorCodes.ErrorReadingMassCorrectionTagsFile
                strMessage = "Error reading mass correction tags file"
            Case ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile
                strMessage = "Error reading modification definitions file"

            Case ePHRPErrorCodes.FilePathError
                strMessage = "General file path error"
            Case ePHRPErrorCodes.UnspecifiedError
                strMessage = "Unspecified error"
            Case Else
                ' This shouldn't happen
                strMessage = "Unknown error state"
        End Select

        If mErrorMessage.Length > 0 Then
            If strMessage.Length > 0 Then
                strMessage &= "; "
            End If
            strMessage &= mErrorMessage
        End If

        Return strMessage

    End Function

    Private Function GetVersionForExecutingAssembly() As String

        Dim strVersion As String

        Try
            strVersion = System.Reflection.Assembly.GetExecutingAssembly.GetName.Version.ToString()
        Catch ex As Exception
            strVersion = "??.??.??.??"
        End Try

        Return strVersion

    End Function

    Private Sub InitializeLocalVariables()

        mErrorCode = ePHRPErrorCodes.NoError
        mErrorMessage = String.Empty
        mWarnMissingParameterFileSection = False

        mCreateModificationSummaryFile = True
        mCreateProteinModsFile = False
        mFastaFilePath = String.Empty
        mIgnorePeptideToProteinMapperErrors = False
        mProteinModsFileIncludesReversedProteins = False
        mUseExistingMTSPepToProteinMapFile = False

        mCreateInspectOrMSGFDbFirstHitsFile = True
        mCreateInspectOrMSGFDbSynopsisFile = True

        mMassCorrectionTagsFilePath = String.Empty
        mModificationDefinitionsFilePath = String.Empty
        mSearchToolParameterFilePath = String.Empty

        mInspectSynopsisFilePValueThreshold = clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD

        mMODaMODPlusSynopsisFileProbabilityThreshold = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD

        mMSAlignSynopsisFilePValueThreshold = clsMSAlignResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD

        mMSGFDBSynopsisFilePValueThreshold = clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD
        mMSGFDBSynopsisFileSpecProbThreshold = clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPECPROB_THRESHOLD

        mEnzymeMatchSpec = GetDefaultEnzymeMatchSpec()

        mPeptideNTerminusMassChange = clsPeptideMassCalculator.DEFAULT_N_TERMINUS_MASS_CHANGE
        mPeptideCTerminusMassChange = clsPeptideMassCalculator.DEFAULT_C_TERMINUS_MASS_CHANGE

        mPeptideSeqMassCalculator = New clsPeptideMassCalculator()

        ' Initialize mPeptideMods
        mPeptideMods = New clsPeptideModificationContainer

        ' Initialize mUniqueSequences
        mUniqueSequences = New clsUniqueSequencesContainer

        ' Initialize mSeqToProteinMap
        mSeqToProteinMap = New Hashtable

        ' Define a RegEx to replace all of the non-letter characters
        mReplaceSymbols = New Regex("[^A-Za-z]", RegexOptions.Compiled)

    End Sub

    Protected Function InitializeSequenceOutputFiles(strBaseOutputFilePath As String) As Boolean

        ' Initializes the StreamWriter objects using strBaseOutputFilePath as a base name and replacing the suffix with the default suffix names
        ' Returns True if success; does not catch errors; they will be thrown to the calling function if they occur
        Dim fiOutputFileInfo As FileInfo
        Dim strResultToSeqMapFilePath As String
        Dim strSeqInfoFilePath As String
        Dim strModDetailsFilePath As String
        Dim strSeqToProteinMapFilePath As String

        fiOutputFileInfo = New FileInfo(strBaseOutputFilePath)

        ' Initialize the file paths based on strBaseOutputFilePath
        strResultToSeqMapFilePath = ReplaceFilenameSuffix(fiOutputFileInfo, FILENAME_SUFFIX_RESULT_TO_SEQ_MAP)
        strSeqInfoFilePath = ReplaceFilenameSuffix(fiOutputFileInfo, FILENAME_SUFFIX_SEQ_INFO)
        strModDetailsFilePath = ReplaceFilenameSuffix(fiOutputFileInfo, FILENAME_SUFFIX_MOD_DETAILS)
        strSeqToProteinMapFilePath = ReplaceFilenameSuffix(fiOutputFileInfo, FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP)

        ' Clear the unique sequences container
        mUniqueSequences.Clear()

        ' Clear the sequence to protein map
        mSeqToProteinMap.Clear()

        ' Initialize the ResultToSeqMap file
        mResultToSeqMapFile = New StreamWriter(strResultToSeqMapFilePath)
        mResultToSeqMapFile.WriteLine("Result_ID" & SEP_CHAR & COLUMN_NAME_UNIQUE_SEQ_ID)


        ' Initialize the SeqInfo file
        mSeqInfoFile = New StreamWriter(strSeqInfoFilePath, False)
        mSeqInfoFile.WriteLine(COLUMN_NAME_UNIQUE_SEQ_ID & SEP_CHAR &
          "Mod_Count" & SEP_CHAR &
          "Mod_Description" & SEP_CHAR &
          "Monoisotopic_Mass")

        ' Initialize the ModDetails file
        mModDetailsFile = New StreamWriter(strModDetailsFilePath)
        mModDetailsFile.WriteLine(COLUMN_NAME_UNIQUE_SEQ_ID & SEP_CHAR &
          "Mass_Correction_Tag" & SEP_CHAR &
          "Position")


        ' Initialize the SeqToProtein map file
        mSeqToProteinMapFile = New StreamWriter(strSeqToProteinMapFilePath, False)
        mSeqToProteinMapFile.WriteLine(COLUMN_NAME_UNIQUE_SEQ_ID & SEP_CHAR &
         "Cleavage_State" & SEP_CHAR &
         "Terminus_State" & SEP_CHAR &
         COLUMN_NAME_PROTEIN_NAME & SEP_CHAR &
         "Protein_Expectation_Value_Log(e)" & SEP_CHAR &
         "Protein_Intensity_Log(I)")

        Return True

    End Function

    ''' <summary>
    ''' Returns true if the character is a letter between A and Z or a and z
    ''' </summary>
    ''' <param name="chChar">Character to examine</param>
    ''' <returns></returns>
    ''' <remarks>The Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters; use this function to only return True for normal letters between A and Z</remarks>
    Public Shared Function IsLetterAtoZ(chChar As Char) As Boolean

        Static reIsLetter As Regex = New Regex("[A-Za-z]", RegexOptions.Compiled)

        If reIsLetter.IsMatch(chChar) Then
            Return True
        Else
            Return False
        End If

    End Function

    ''' <summary>
    ''' Examines the string to determine if it is numeric
    ''' </summary>
    ''' <param name="strData"></param>
    ''' <returns>True if a number, otherwise false</returns>
    Public Shared Function IsNumber(strData As String) As Boolean
        Try
            If Double.TryParse(strData, 0) Then
                Return True
            ElseIf Integer.TryParse(strData, 0) Then
                Return True
            End If
        Catch ex As Exception
            ' Ignore errors here
        End Try

        Return False
    End Function

    Protected Function IsReversedProtein(strProteinName As String) As Boolean

        Dim strProteinNameLCase As String = strProteinName.ToLower()

        If strProteinNameLCase.StartsWith("reversed_") Then
            ' Used in DMS-generated protein collections
            Return True
        ElseIf strProteinName.StartsWith("REV_") Then
            ' Used by MSGFDB
            Return True
        ElseIf strProteinNameLCase.StartsWith("scrambled_") Then
            ' Used in DMS-generated protein collections
            Return True
        ElseIf strProteinNameLCase.StartsWith("xxx_") Then
            ' Used by MSGF+
            Return True
        ElseIf strProteinNameLCase.StartsWith("xxx.") Then
            ' Used by Inspect
            Return True
        ElseIf strProteinNameLCase.EndsWith(":reversed") Then
            ' Used by X!Tandem
            Return True
        End If

        Return False

    End Function

    Protected Overridable Function LoadParameterFileSettings(strParameterFilePath As String) As Boolean

        Const OPTIONS_SECTION = "PeptideHitResultsProcessorOptions"

        Dim objSettingsFile As New XmlSettingsFileAccessor

        Dim strLeftResidueRegEx As String, strRightResidueRegEx As String
        Dim blnValueNotPresent As Boolean

        Try

            If String.IsNullOrWhiteSpace(strParameterFilePath) Then
                ' No parameter file specified; nothing to load
                Return True
            End If

            If Not File.Exists(strParameterFilePath) Then
                ' See if strParameterFilePath points to a file in the same directory as the application
                strParameterFilePath = Path.Combine(Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location), Path.GetFileName(strParameterFilePath))
                If Not File.Exists(strParameterFilePath) Then
                    SetErrorCode(ePHRPErrorCodes.ParameterFileNotFound)
                    Return False
                End If
            End If

            If objSettingsFile.LoadSettings(strParameterFilePath) Then
                If Not objSettingsFile.SectionPresent(OPTIONS_SECTION) Then
                    ' Section OPTIONS_SECTION was not found in the parameter file; warn the user if mWarnMissingParameterFileSection = True
                    If mWarnMissingParameterFileSection Then
                        SetErrorMessage("The node '<section name=""" & OPTIONS_SECTION & """> was not found in the parameter file: " & strParameterFilePath)
                        SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile)
                        Return False
                    Else
                        Return True
                    End If
                Else
                    mMassCorrectionTagsFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "MassCorrectionTagsFilePath", mMassCorrectionTagsFilePath)
                    mModificationDefinitionsFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "ModificationDefinitionsFilePath", mModificationDefinitionsFilePath)
                    mSearchToolParameterFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "SearchToolParameterFilePath", mSearchToolParameterFilePath)

                    mCreateModificationSummaryFile = objSettingsFile.GetParam(OPTIONS_SECTION, "CreateModificationSummaryFile", mCreateModificationSummaryFile)

                    mCreateProteinModsFile = objSettingsFile.GetParam(OPTIONS_SECTION, "CreateProteinModsFile", mCreateProteinModsFile)
                    mFastaFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "FastaFilePath", mFastaFilePath)
                    mProteinModsFileIncludesReversedProteins = objSettingsFile.GetParam(OPTIONS_SECTION, "ProteinModsFileIncludesReversedProteins", mProteinModsFileIncludesReversedProteins)
                    mUseExistingMTSPepToProteinMapFile = objSettingsFile.GetParam(OPTIONS_SECTION, "UseExistingMTSPepToProteinMapFile", mUseExistingMTSPepToProteinMapFile)

                    With mEnzymeMatchSpec
                        strLeftResidueRegEx = String.Copy(.LeftResidueRegEx)
                        strRightResidueRegEx = String.Copy(.RightResidueRegEx)
                    End With

                    strLeftResidueRegEx = objSettingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecLeftResidue", strLeftResidueRegEx, blnValueNotPresent)
                    If Not blnValueNotPresent Then
                        strRightResidueRegEx = objSettingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecRightResidue", strRightResidueRegEx, blnValueNotPresent)

                        If Not blnValueNotPresent Then
                            With mEnzymeMatchSpec
                                .LeftResidueRegEx = strLeftResidueRegEx
                                .RightResidueRegEx = strRightResidueRegEx
                            End With
                        End If
                    End If

                    mPeptideNTerminusMassChange = objSettingsFile.GetParam(OPTIONS_SECTION, "PeptideNTerminusMassChange", mPeptideNTerminusMassChange)
                    mPeptideCTerminusMassChange = objSettingsFile.GetParam(OPTIONS_SECTION, "PeptideCTerminusMassChange", mPeptideCTerminusMassChange)

                End If
            End If

        Catch ex As Exception
            SetErrorMessage("Error in LoadParameterFileSettings:" & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingParameterFile)
            Return False
        End Try

        Return True

    End Function

    ''' <summary>
    ''' Load the PeptideToProteinMap information
    ''' </summary>
    ''' <param name="strPepToProteinMapFilePath">File to read</param>
    ''' <param name="lstPepToProteinMapping">Output parameter: peptide to protein mapping (calling function must pre-initialize the list)</param>
    ''' <param name="strHeaderLine">Output parameter: Header line text</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function LoadPeptideToProteinMapInfo(
       strPepToProteinMapFilePath As String,
       lstPepToProteinMapping As List(Of udtPepToProteinMappingType),
       ByRef strHeaderLine As String) As Boolean

        Dim strLineIn As String
        Dim strSplitLine As String()

        Dim intLinesRead As Integer
        Dim intValue As Integer

        Dim blnSuccess As Boolean = False

        Try

            ' Initialize the output parameters
            lstPepToProteinMapping.Clear()
            strHeaderLine = String.Empty

            If String.IsNullOrWhiteSpace(strPepToProteinMapFilePath) Then
                SetErrorMessage("Warning: PepToProteinMap file is not defined")
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.InvalidInputFilePath)
                Return False
            ElseIf Not File.Exists(strPepToProteinMapFilePath) Then
                SetErrorMessage("Warning: PepToProteinMap file does not exist: " & strPepToProteinMapFilePath)
                SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.InvalidInputFilePath)
                Return False
            End If

            ' Open strProteinToPeptideMappingFilePath for reading
            Using srInFile = New StreamReader(New FileStream(strPepToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                intLinesRead = 0
                Do While Not srInFile.EndOfStream
                    strLineIn = srInFile.ReadLine().Trim()

                    If strLineIn.Length > 0 Then

                        ' Split the line on tabs
                        strSplitLine = strLineIn.Split(ControlChars.Tab)

                        If strSplitLine.Length >= 4 Then

                            If intLinesRead = 0 AndAlso Not Integer.TryParse(strSplitLine(2), intValue) Then
                                ' Header line; cache it
                                strHeaderLine = String.Copy(strLineIn)
                            Else
                                Dim udtPepToProteinMappingEntry As udtPepToProteinMappingType
                                With udtPepToProteinMappingEntry
                                    .Peptide = String.Copy(strSplitLine(0))
                                    .Protein = String.Copy(strSplitLine(1))
                                    Integer.TryParse(strSplitLine(2), .ResidueStart)
                                    Integer.TryParse(strSplitLine(3), .ResidueEnd)
                                End With

                                ExpandListIfRequired(lstPepToProteinMapping, 1)

                                lstPepToProteinMapping.Add(udtPepToProteinMappingEntry)

                            End If
                        End If
                    End If

                Loop

            End Using

            blnSuccess = True

        Catch ex As Exception
            SetErrorMessage("Error reading the Peptide to Protein Map File (" & Path.GetFileName(strPepToProteinMapFilePath) & "): " & ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Function NumToString(intNumber As Integer, intDigitsOfPrecision As Integer, blnRemoveDecimalsWhenZero As Boolean) As String
        Return intNumber.ToString()
    End Function

    Protected Function NumToString(sngNumber As Single, intDigitsOfPrecision As Integer, blnRemoveDecimalsWhenZero As Boolean) As String
        Return NumToString(CDbl(sngNumber), intDigitsOfPrecision, blnRemoveDecimalsWhenZero)
    End Function

    Protected Function NumToString(dblNumber As Double, intDigitsOfPrecision As Integer, blnRemoveDecimalsWhenZero As Boolean) As String
        Static strFormatString As String = "0"
        Static intFormatStringPrecision As Integer = 0

        If blnRemoveDecimalsWhenZero AndAlso Math.Abs(dblNumber - 0) < Double.Epsilon Then
            Return "0"
        Else
            If intFormatStringPrecision <> intDigitsOfPrecision Then
                ' Update strFormatString
                If intDigitsOfPrecision <= 0 Then
                    strFormatString = "0"
                Else
                    strFormatString = "0." & New String("0"c, intDigitsOfPrecision)
                End If
                intFormatStringPrecision = intDigitsOfPrecision
            End If

            Return dblNumber.ToString(strFormatString)
        End If

    End Function

    Protected Sub OperationComplete()
        RaiseEvent ProgressComplete()
    End Sub

    Public Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String) As Boolean
        Return ProcessFile(strInputFilePath, strOutputFolderPath, String.Empty)
    End Function

    Public MustOverride Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean

    Protected Function ReplaceFilenameSuffix(fiOriginalFile As FileInfo, strNewSuffix As String) As String
        ' Appends strNewSuffix to the base name of the original file, then returns a full path to the file using the folder associated with strOriginalFilePath
        ' Note that strNewSuffix may contain a file extension though it does not have to
        '  If strNewSuffix does not contain an extension, then the path returned will end in the same extension as strOriginalFilePath

        Dim strNewFileName As String
        Dim strOriginalExtension As String

        ' Keep track of the original extension on strOriginalFilePath
        strOriginalExtension = fiOriginalFile.Extension

        ' Make sure strNewSuffix is not nothing
        If strNewSuffix Is Nothing Then strNewSuffix = String.Empty

        ' Obtain the filename, without its extension
        strNewFileName = Path.GetFileNameWithoutExtension(fiOriginalFile.Name)

        ' Append strNewSuffix to strNewFileName
        If Path.HasExtension(strNewSuffix) Then
            strNewFileName &= strNewSuffix
        Else
            strNewFileName &= strNewSuffix & strOriginalExtension
        End If

        strNewFileName = Path.Combine(fiOriginalFile.DirectoryName, strNewFileName)

        Return strNewFileName

    End Function

    Protected Sub ReportError(errMsg As String, Optional throwException As Boolean = False, Optional innerException As Exception = Nothing)

        SetErrorMessage(errMsg)

        If throwException Then
            If innerException Is Nothing Then
                Throw New Exception(errMsg)
            Else
                Throw New Exception(errMsg, innerException)
            End If

        End If

    End Sub

    Protected Sub ReportWarning(strMessage As String)
        RaiseEvent WarningMessageEvent(strMessage)
    End Sub

    Public Function ResetMassCorrectionTagsAndModificationDefinitions() As Boolean
        Dim blnSuccess As Boolean
        Dim blnFileNotFound As Boolean

        ' Note: If mMassCorrectionTagsFilePath is blank then the mass correction tags will be reset to the defaults and blnSuccess will be True
        blnSuccess = mPeptideMods.ReadMassCorrectionTagsFile(mMassCorrectionTagsFilePath, blnFileNotFound)
        If Not blnSuccess Then
            If blnFileNotFound Then
                SetErrorCode(ePHRPErrorCodes.MassCorrectionTagsFileNotFound)
            Else
                SetErrorCode(ePHRPErrorCodes.ErrorReadingMassCorrectionTagsFile)
            End If
        End If

        ' Note: If mModificationDefinitionsFilePath is blank, then the modifications will be cleared and blnSuccess will be True
        blnSuccess = mPeptideMods.ReadModificationDefinitionsFile(mModificationDefinitionsFilePath, blnFileNotFound)
        If Not blnSuccess Then
            If blnFileNotFound Then
                SetErrorCode(ePHRPErrorCodes.ModificationDefinitionFileNotFound)
                ReportWarning("File not found: " & mModificationDefinitionsFilePath)
            Else
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            End If
        End If

        Return blnSuccess

    End Function

    Protected Sub ResetProgress()
        ResetProgress(String.Empty)
    End Sub

    Protected Sub ResetProgress(strProgressStepDescription As String)
        mProgressStepDescription = String.Copy(strProgressStepDescription)
        mProgressPercentComplete = 0
        RaiseEvent ProgressReset()
    End Sub

    Protected Sub SaveModificationSummaryFile(strModificationSummaryFilePath As String)
        Dim intIndex As Integer

        Try
            Using swOutFile As StreamWriter = New StreamWriter(strModificationSummaryFilePath, False)

                ' Write the header line
                swOutFile.WriteLine(
                  clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol & SEP_CHAR &
                  clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass & SEP_CHAR &
                  clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Target_Residues & SEP_CHAR &
                  clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Type & SEP_CHAR &
                  clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Mass_Correction_Tag & SEP_CHAR &
                  clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Occurrence_Count)

                For intIndex = 0 To mPeptideMods.ModificationCount - 1
                    Dim oModInfo = mPeptideMods.GetModificationByIndex(intIndex)
                    With oModInfo
                        If .OccurrenceCount > 0 OrElse Not .UnknownModAutoDefined Then
                            swOutFile.WriteLine(.ModificationSymbol & SEP_CHAR &
                              .ModificationMass.ToString & SEP_CHAR &
                              .TargetResidues & SEP_CHAR &
                              clsModificationDefinition.ModificationTypeToModificationSymbol(.ModificationType) & SEP_CHAR &
                              .MassCorrectionTag & SEP_CHAR &
                              .OccurrenceCount.ToString)
                        End If
                    End With
                Next intIndex
            End Using

        Catch ex As Exception
            Throw ex
        End Try

    End Sub

    Protected Sub SaveResultsFileEntrySeqInfo(objSearchResult As clsSearchResultsBaseClass, blnUpdateResultToSeqMapFile As Boolean)
        ' Note: Be sure to call Me.InitializeOutputFiles before calling this function
        ' blnUpdateResultToSeqMapFile should be set to True only for the first protein of each peptide in each group

        Dim intIndex As Integer
        Dim intUniqueSeqID As Integer   ' This ID is assigned using a hashtable containing mPeptideCleanSequence and mPeptideModDescription
        Dim blnExistingSequenceFound As Boolean

        Dim udtModNameAndResidueLoc() As udtModNameAndResidueLocType
        Dim intPointerArray() As Integer

        With objSearchResult
            intUniqueSeqID = mUniqueSequences.GetNextUniqueSequenceID(.PeptideCleanSequence, .PeptideModDescription, blnExistingSequenceFound)

            If blnUpdateResultToSeqMapFile Then
                ' Write a new entry to the ResultToSeqMap file
                mResultToSeqMapFile.WriteLine(.ResultID.ToString & SEP_CHAR & intUniqueSeqID.ToString)

                ' Only write this entry to the SeqInfo and ModDetails files if blnExistingSequenceFound is False
                If Not blnExistingSequenceFound Then

                    ' Write a new entry to the SeqInfo file
                    mSeqInfoFile.WriteLine(
                     intUniqueSeqID.ToString & SEP_CHAR &
                     .SearchResultModificationCount & SEP_CHAR &
                     .PeptideModDescription & SEP_CHAR &
                     .PeptideMonoisotopicMass.ToString("0.0000000"))


                    If .SearchResultModificationCount > 0 Then
                        ReDim udtModNameAndResidueLoc(.SearchResultModificationCount - 1)
                        ReDim intPointerArray(.SearchResultModificationCount - 1)

                        If .SearchResultModificationCount = 1 Then
                            intPointerArray(0) = 0
                        Else
                            ' Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                            For intIndex = 0 To .SearchResultModificationCount - 1
                                With .GetSearchResultModDetailsByIndex(intIndex)
                                    udtModNameAndResidueLoc(intIndex).ResidueLocInPeptide = .ResidueLocInPeptide
                                    udtModNameAndResidueLoc(intIndex).ModName = .ModDefinition.MassCorrectionTag
                                End With
                                intPointerArray(intIndex) = intIndex
                            Next intIndex

                            Array.Sort(udtModNameAndResidueLoc, intPointerArray, New IModNameAndResidueLocComparer)
                        End If

                        ' Write out the modifications to the ModDetails file
                        ' Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
                        For intIndex = 0 To .SearchResultModificationCount - 1
                            With .GetSearchResultModDetailsByIndex(intPointerArray(intIndex))
                                mModDetailsFile.WriteLine(
                                 intUniqueSeqID.ToString & SEP_CHAR &
                                 .ModDefinition.MassCorrectionTag & SEP_CHAR &
                                 .ResidueLocInPeptide.ToString)
                            End With
                        Next intIndex
                    End If
                End If
            End If

            ' Write a new entry to the SeqToProteinMap file if not yet defined
            If Not CheckSeqToProteinMapDefined(intUniqueSeqID, .ProteinName) Then
                mSeqToProteinMapFile.WriteLine(
                   intUniqueSeqID.ToString & SEP_CHAR &
                   CInt(.PeptideCleavageState).ToString & SEP_CHAR &
                   CInt(.PeptideTerminusState).ToString & SEP_CHAR &
                   .ProteinName & SEP_CHAR &
                   .ProteinExpectationValue & SEP_CHAR &
                   .ProteinIntensity)
            End If
        End With

    End Sub

    Protected Sub SetErrorCode(eNewErrorCode As ePHRPErrorCodes)
        SetErrorCode(eNewErrorCode, False)
    End Sub

    Protected Sub SetErrorCode(eNewErrorCode As ePHRPErrorCodes, blnLeaveExistingErrorCodeUnchanged As Boolean)
        If blnLeaveExistingErrorCodeUnchanged AndAlso mErrorCode <> ePHRPErrorCodes.NoError Then
            ' An error code is already defined; do not change it
        Else
            mErrorCode = eNewErrorCode
        End If
    End Sub

    Protected Sub SetErrorMessage(strMessage As String)
        If strMessage Is Nothing Then strMessage = String.Empty
        mErrorMessage = String.Copy(strMessage)
        If strMessage.Length > 0 Then
            RaiseEvent ErrorOccurred(mErrorMessage)
            Console.WriteLine(mErrorMessage)
        End If
    End Sub

    ''' <summary>
    ''' Return the text up to (but not including) the first space in strProteinNameAndDescription
    ''' </summary>
    ''' <param name="strProteinNameAndDescription"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Overridable Function TruncateProteinName(strProteinNameAndDescription As String) As String

        Dim intIndex As Integer

        intIndex = strProteinNameAndDescription.IndexOf(" "c)
        If intIndex > 0 Then
            Return strProteinNameAndDescription.Substring(0, intIndex)
        Else
            Return strProteinNameAndDescription
        End If

    End Function

    Protected Sub UpdatePepToProteinMapPeptide(ByRef lstPepToProteinMapping As List(Of udtPepToProteinMappingType), intIndex As Integer, strPeptide As String)
        Dim udtItem As udtPepToProteinMappingType
        udtItem = lstPepToProteinMapping(intIndex)
        udtItem.Peptide = strPeptide
        lstPepToProteinMapping(intIndex) = udtItem
    End Sub

    Protected Sub UpdateProgress(strProgressStepDescription As String)
        UpdateProgress(strProgressStepDescription, mProgressPercentComplete)
    End Sub

    Protected Sub UpdateProgress(sngPercentComplete As Single)
        UpdateProgress(Me.ProgressStepDescription, sngPercentComplete)
    End Sub

    Protected Sub UpdateProgress(strProgressStepDescription As String, sngPercentComplete As Single)
        mProgressStepDescription = String.Copy(strProgressStepDescription)
        If sngPercentComplete < 0 Then
            sngPercentComplete = 0
        ElseIf sngPercentComplete > 100 Then
            sngPercentComplete = 100
        End If
        mProgressPercentComplete = sngPercentComplete

        RaiseEvent ProgressChanged(Me.ProgressStepDescription, Me.ProgressPercentComplete)
    End Sub

    Private Function ValidatePeptideToProteinMapResults(strPeptideToProteinMapFilePath As String, blnIgnorePeptideToProteinMapperErrors As Boolean) As Boolean

        Const PROTEIN_NAME_NO_MATCH As String = "__NoMatch__"

        Dim blnSuccess As Boolean = False

        Dim intPeptideCount As Integer = 0
        Dim intPeptideCountNoMatch As Integer = 0
        Dim intLinesRead As Integer = 0
        Dim chSplitChars() As Char = New Char() {ControlChars.Tab}

        Try
            ' Validate that none of the results in strPeptideToProteinMapFilePath has protein name PROTEIN_NAME_NO_MATCH

            Dim strLineIn As String
            Dim strLastPeptide As String = String.Empty
            Dim strSplitLine() As String

            Using srInFile = New StreamReader(New FileStream(strPeptideToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                Do While Not srInFile.EndOfStream
                    strLineIn = srInFile.ReadLine()
                    intLinesRead += 1

                    If intLinesRead > 1 AndAlso Not String.IsNullOrEmpty(strLineIn) Then

                        strSplitLine = strLineIn.Split(chSplitChars, 2)
                        If strSplitLine.Count > 0 Then
                            If strSplitLine(0) <> strLastPeptide Then
                                intPeptideCount += 1
                                strLastPeptide = String.Copy(strSplitLine(0))
                            End If

                            If strLineIn.Contains(PROTEIN_NAME_NO_MATCH) Then
                                intPeptideCountNoMatch += 1
                            End If
                        End If

                    End If
                Loop

            End Using

            If intPeptideCount = 0 Then
                SetErrorMessage("Peptide to protein mapping file is empty: " & strPeptideToProteinMapFilePath)
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
                blnSuccess = False

            ElseIf intPeptideCountNoMatch = 0 Then
                blnSuccess = True

            Else
                Dim dblErrorPercent As Double   ' Value between 0 and 100
                dblErrorPercent = intPeptideCountNoMatch / intPeptideCount * 100.0

                Dim strErrorMessage As String
                strErrorMessage = dblErrorPercent.ToString("0.00") & "% of the entries (" & intPeptideCountNoMatch & " / " & intPeptideCount & ") in the peptide to protein map file (" & IO.Path.GetFileName(strPeptideToProteinMapFilePath) & ") did not match to a protein in the FASTA file (" & IO.Path.GetFileName(mFastaFilePath) & ")"

                If blnIgnorePeptideToProteinMapperErrors OrElse dblErrorPercent < 0.1 Then
                    ReportWarning(strErrorMessage)
                    blnSuccess = True
                Else
                    SetErrorMessage(strErrorMessage)
                    blnSuccess = False
                End If
            End If

        Catch ex As Exception
            SetErrorMessage("Error in ValidatePeptideToProteinMapResults:" & ex.Message)
            SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Sub ValidatePHRPReaderSupportFiles(strPHRPDataFilePath As String, strOutputFolderPath As String)

        Dim ioPHRPFile As FileInfo
        Dim ioOutputFolder As DirectoryInfo
        Dim strMSGFFileName As String
        Dim strSourcePath As String
        Dim strTargetPath As String

        Try
            If Not String.IsNullOrWhiteSpace(strOutputFolderPath) Then

                ioPHRPFile = New FileInfo(strPHRPDataFilePath)
                ioOutputFolder = New DirectoryInfo(strOutputFolderPath)

                If ioPHRPFile.DirectoryName.ToLower() <> ioOutputFolder.FullName.ToLower() Then
                    strMSGFFileName = IO.Path.GetFileName(ReplaceFilenameSuffix(ioPHRPFile, FILENAME_SUFFIX_MSGF))

                    strSourcePath = IO.Path.Combine(ioPHRPFile.DirectoryName, strMSGFFileName)
                    strTargetPath = IO.Path.Combine(ioOutputFolder.FullName, strMSGFFileName)

                    If IO.File.Exists(strSourcePath) And Not IO.File.Exists(strTargetPath) Then
                        IO.File.Copy(strSourcePath, strTargetPath)
                    End If

                End If

            End If

        Catch ex As Exception
            ReportWarning("Error in ValidatePHRPReaderSupportFiles: " & ex.Message)
        End Try

    End Sub

    Protected Function ValidateProteinFastaFile(strFastaFilePath As String) As Boolean

        Dim blnSuccess As Boolean
        Dim strWarningMessage As String = String.Empty
        blnSuccess = ValidateProteinFastaFile(strFastaFilePath, strWarningMessage)

        If Not blnSuccess Then
            ReportWarning(strWarningMessage)
        End If

        Return blnSuccess
    End Function

    Public Shared Function ValidateProteinFastaFile(strFastaFilePath As String, ByRef strWarningMessage As String) As Boolean

        Dim objFastaFile As ProteinFileReader.FastaFileReader

        ' This RegEx looks for standard amino acids, skipping A, T, C, and G
        Dim reDefiniteAminoAcid As Regex = New Regex("[DEFHIKLMNPQRSVWY]", RegexOptions.IgnoreCase Or RegexOptions.Compiled)

        ' This RegEx looks for A, T, C, and G
        Dim rePotentialNucleicAcid As Regex = New Regex("[ATCG]", RegexOptions.IgnoreCase Or RegexOptions.Compiled)

        ' This matches any letter
        Dim reLetter As Regex = New Regex("[A-Z]", RegexOptions.IgnoreCase Or RegexOptions.Compiled)

        Dim intDefiniteAminoAcidCount As Integer
        Dim intPotentialNucleicAcidCount As Integer
        Dim intLetterCount As Integer

        Dim intValidProteinCount As Integer = 0
        Dim intInvalidProteinCount As Integer = 0

        Try
            strWarningMessage = String.Empty

            If String.IsNullOrEmpty(strFastaFilePath) Then
                Console.WriteLine()
                strWarningMessage = "strFastaFilePath is not defined in ValidateProteinFastaFile"
                Return False
            ElseIf Not IO.File.Exists(strFastaFilePath) Then
                Console.WriteLine()
                strWarningMessage = "Fasta file not found: " & strFastaFilePath
                Return False
            End If

            objFastaFile = New ProteinFileReader.FastaFileReader()
            If Not objFastaFile.OpenFile(strFastaFilePath) Then
                Console.WriteLine()
                strWarningMessage = "Error opening the fasta file: " & strFastaFilePath
                Return False
            End If

            ' Read the first 500 proteins and confirm that each contains amino acid residues
            Do While objFastaFile.ReadNextProteinEntry()
                intDefiniteAminoAcidCount = reDefiniteAminoAcid.Matches(objFastaFile.ProteinSequence).Count
                intPotentialNucleicAcidCount = rePotentialNucleicAcid.Matches(objFastaFile.ProteinSequence).Count
                intLetterCount = reLetter.Matches(objFastaFile.ProteinSequence).Count

                If intDefiniteAminoAcidCount > 0.1 * intLetterCount Then
                    intValidProteinCount += 1
                ElseIf intPotentialNucleicAcidCount > 0.95 * intLetterCount Then
                    intInvalidProteinCount += 1
                End If

                If intValidProteinCount + intInvalidProteinCount >= 500 Then
                    Exit Do
                End If
            Loop

            If intValidProteinCount < intInvalidProteinCount Then
                Console.WriteLine()
                strWarningMessage = "Fasta file contains Nucleic Acids, not Amino Acids: " & IO.Path.GetFileName(strFastaFilePath)
                Return False
            End If

        Catch ex As Exception
            Console.WriteLine()
            strWarningMessage = "Exception in ValidateProteinFastaFile: " & ex.Message
            Return False
        End Try

        Return True

    End Function

#Region "PHRPReader Event Handlers"

    Protected Sub PHRPReader_ErrorEvent(strErrorMessage As String)
        RaiseEvent ErrorOccurred(strErrorMessage)
    End Sub

    Protected Sub PHRPReader_WarningEvent(strWarningMessage As String)
        ReportWarning(strWarningMessage)
    End Sub
#End Region


#Region "PeptideToProteinMapper Event Handlers"

    Private Sub PeptideToProteinMapper_ProgressChanged(taskDescription As String, percentComplete As Single)
        If percentComplete >= mNextPeptideToProteinMapperLevel Then
            mNextPeptideToProteinMapperLevel += 25
            UpdateProgress(PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE + percentComplete * (PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE - PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE) / 100)
            Console.WriteLine(" PeptideToProteinMapper is " & percentComplete.ToString("0") & "% complete")
        End If
    End Sub

#End Region

#Region "IComparer classes"

	Protected Class ISearchOptionModificationInfoComparer
		Implements IComparer(Of udtSearchOptionModificationInfoType)

		Public Function Compare(x As udtSearchOptionModificationInfoType, y As udtSearchOptionModificationInfoType) As Integer Implements IComparer(Of udtSearchOptionModificationInfoType).Compare

			If x.SortOrder > y.SortOrder Then
				Return 1
			ElseIf x.SortOrder < y.SortOrder Then
				Return -1
			Else
				If x.ModificationMass > y.ModificationMass Then
					Return 1
				ElseIf x.ModificationMass < y.ModificationMass Then
					Return -1
				Else
					Return 0
				End If
			End If
		End Function

	End Class

	Friend Class IModNameAndResidueLocComparer
		Implements IComparer(Of udtModNameAndResidueLocType)

		Public Function Compare(x As udtModNameAndResidueLocType, y As udtModNameAndResidueLocType) As Integer Implements IComparer(Of udtModNameAndResidueLocType).Compare

			If x.ResidueLocInPeptide > y.ResidueLocInPeptide Then
				Return 1
			ElseIf x.ResidueLocInPeptide < y.ResidueLocInPeptide Then
				Return -1
			Else
				If x.ModName Is Nothing Then x.ModName = String.Empty
				If y.ModName Is Nothing Then y.ModName = String.Empty

				If x.ModName > y.ModName Then
					Return 1
				ElseIf x.ModName < y.ModName Then
					Return -1
				Else
					Return 0
				End If
			End If
		End Function

	End Class

	Protected Class PepToProteinMappingComparer
		Implements IComparer(Of udtPepToProteinMappingType)

		Public Function Compare(x As udtPepToProteinMappingType, y As udtPepToProteinMappingType) As Integer Implements IComparer(Of udtPepToProteinMappingType).Compare

			If x.Peptide > y.Peptide Then
				Return 1
			ElseIf x.Peptide < y.Peptide Then
				Return -1
			Else
				If x.Protein > y.Protein Then
					Return 1
				ElseIf x.Protein < y.Protein Then
					Return -1
				Else
					Return 0
				End If
			End If

		End Function

	End Class

	Protected Class PepToProteinMappingPeptideSearchComparer
		Implements IComparer(Of udtPepToProteinMappingType)

		Public Function Compare(x As clsPHRPBaseClass.udtPepToProteinMappingType, y As clsPHRPBaseClass.udtPepToProteinMappingType) As Integer Implements IComparer(Of clsPHRPBaseClass.udtPepToProteinMappingType).Compare

			If x.Peptide > y.Peptide Then
				Return 1
			ElseIf x.Peptide < y.Peptide Then
				Return -1
			Else
				Return 0
			End If

		End Function

	End Class

#End Region

End Class
