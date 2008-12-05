Option Strict On

' This class can be used as a base class for peptide hit results processor classes
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
' Started January 6, 2006
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

Public MustInherit Class clsPHRPBaseClass

    Public Sub New()
        mFileDate = "December 5, 2008"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
    Protected Const SEP_CHAR As Char = ControlChars.Tab
    Protected Const UNIQUE_SEQ_TO_PROTEIN_MAP_SEP As String = "_"c

    Public Const XTANDEM_RESULTS_FILE_SUFFIX As String = "_xt.xml"
    Public Const SEQUEST_SYNOPSIS_FILE_SUFFIX As String = "_syn.txt"
    Public Const SEQUEST_FIRST_HITS_FILE_SUFFIX As String = "_fht.txt"

    Public Const INSPECT_RESULTS_FILE_SUFFIX As String = "_inspect.txt"
    Public Const INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX As String = "_fht.txt"
    Public Const INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX As String = "_Fscore_fht.txt"

    Public Const FILENAME_SUFFIX_RESULT_TO_SEQ_MAP As String = "_ResultToSeqMap.txt"
    Public Const FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP As String = "_SeqToProteinMap.txt"

    Public Const FILENAME_SUFFIX_SEQ_INFO As String = "_SeqInfo.txt"
    Public Const FILENAME_SUFFIX_MOD_DETAILS As String = "_ModDetails.txt"
    Public Const FILENAME_SUFFIX_MOD_SUMMARY As String = "_ModSummary.txt"


    Public Enum ePeptideHitResultsFileFormatConstants As Integer
        AutoDetermine = 0
        SequestSynopsisFile = 1
        SequestFirstHitsFile = 2
        XTandemXMLFile = 3
        InSpectTXTFile = 4
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
    Friend Structure udtSearchOptionModificationInfoType
        Dim SortOrder As Integer
        Dim ModificationMass As Double
        Dim TargetResidues As String
        Public ModificationType As clsModificationDefinition.eModificationTypeConstants
    End Structure

    Friend Structure udtModNameAndResidueLocType
        Dim ModName As String
        Dim ResidueLocInPeptide As Integer
    End Structure
#End Region

#Region "Classwide Variables"

    Protected mErrorCode As ePHRPErrorCodes = ePHRPErrorCodes.NoError
    Protected mErrorMessage As String = String.Empty
    Protected mWarnMissingParameterFileSection As Boolean

    Protected mFileDate As String = String.Empty
    Protected mAbortProcessing As Boolean

    Protected mCreateModificationSummaryFile As Boolean

    Protected mMassCorrectionTagsFilePath As String
    Protected mModificationDefinitionsFilePath As String
    Protected mSearchToolParameterFilePath As String            ' At present, only used by clsInSpecTResultsProcessor

    Protected mEnzymeMatchSpec As clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType
    Protected mPeptideNTerminusMassChange As Double             ' This is ignored if equal to 0; typical non-zero value is 1.0078246
    Protected mPeptideCTerminusMassChange As Double             ' This is ignored if equal to 0; typical non-zero value is 17.0027387

    Protected mPeptideMods As clsPeptideModificationContainer
    Protected mUniqueSequences As clsUniqueSequencesContainer
    Protected mSeqToProteinMap As Hashtable

    Protected mResultToSeqMapFile As System.IO.StreamWriter
    Protected mSeqInfoFile As System.IO.StreamWriter
    Protected mModDetailsFile As System.IO.StreamWriter
    Protected mSeqToProteinMapFile As System.IO.StreamWriter

#End Region

#Region "Progress Events and Variables"
    Public Event ProgressReset()
    Public Event ProgressChanged(ByVal taskDescription As String, ByVal percentComplete As Single)     ' PercentComplete ranges from 0 to 100, but can contain decimal percentage values
    Public Event ProgressComplete()

    Public Event ErrorOccurred(ByVal ErrorMessage As String)

    Protected mProgressStepDescription As String = String.Empty
    Protected mProgressPercentComplete As Single        ' Ranges from 0 to 100, but can contain decimal percentage values
#End Region

#Region "Properties"
    Public Property AbortProcessing() As Boolean
        Get
            Return mAbortProcessing
        End Get
        Set(ByVal Value As Boolean)
            mAbortProcessing = Value
        End Set
    End Property

    Public Property CreateModificationSummaryFile() As Boolean
        Get
            Return mCreateModificationSummaryFile
        End Get
        Set(ByVal Value As Boolean)
            mCreateModificationSummaryFile = Value
        End Set
    End Property

    Public Property EnzymeMatchSpec() As clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType
        Get
            Return mEnzymeMatchSpec
        End Get
        Set(ByVal Value As clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType)
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

    Public Property MassCorrectionTagsFilePath() As String
        Get
            Return mMassCorrectionTagsFilePath
        End Get
        Set(ByVal Value As String)
            mMassCorrectionTagsFilePath = Value
        End Set
    End Property

    Public Property ModificationDefinitionsFilePath() As String
        Get
            Return mModificationDefinitionsFilePath
        End Get
        Set(ByVal Value As String)
            mModificationDefinitionsFilePath = Value
        End Set
    End Property

    Public Property PeptideCTerminusMassChange() As Double
        Get
            Return mPeptideCTerminusMassChange
        End Get
        Set(ByVal Value As Double)
            mPeptideCTerminusMassChange = Value
        End Set
    End Property

    Public Property PeptideNTerminusMassChange() As Double
        Get
            Return mPeptideNTerminusMassChange
        End Get
        Set(ByVal Value As Double)
            mPeptideNTerminusMassChange = Value
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
        Set(ByVal value As String)
            mSearchToolParameterFilePath = value
        End Set
    End Property
    Public Property WarnMissingParameterFileSection() As Boolean
        Get
            Return mWarnMissingParameterFileSection
        End Get
        Set(ByVal Value As Boolean)
            mWarnMissingParameterFileSection = Value
        End Set
    End Property
#End Region

    Public Sub AbortProcessingNow()
        mAbortProcessing = True
    End Sub

    Protected Function CheckSeqToProteinMapDefined(ByVal intUniqueSeqID As Integer, ByVal strProteinName As String) As Boolean
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

    Public Shared Function AutoDefinePeptideHitResultsFilePath(ByVal ePeptideHitResultFileFormat As ePeptideHitResultsFileFormatConstants, ByVal strSourceFolderPath As String, ByVal strBaseName As String) As String

        If Not strBaseName Is Nothing AndAlso strBaseName.Length > 0 Then
            Select Case ePeptideHitResultFileFormat
                Case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile
                    Return System.IO.Path.Combine(strSourceFolderPath, strBaseName & SEQUEST_FIRST_HITS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile
                    Return System.IO.Path.Combine(strSourceFolderPath, strBaseName & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.XTandemXMLFile
                    Return System.IO.Path.Combine(strSourceFolderPath, strBaseName & XTANDEM_RESULTS_FILE_SUFFIX)

                Case ePeptideHitResultsFileFormatConstants.InSpectTXTFile
                    Return System.IO.Path.Combine(strSourceFolderPath, strBaseName & INSPECT_RESULTS_FILE_SUFFIX)

                Case Else
                    ' Includes ePeptideHitResultsFileFormatConstants.AutoDetermine
                    ' Call AutoDefinePeptideHitResultsFilePath below sending it only strSourceFolderPath
            End Select
        End If

        Return AutoDefinePeptideHitResultsFilePath(strSourceFolderPath)
    End Function

    Public Shared Function AutoDefinePeptideHitResultsFilePath(ByVal strSourceFolderPath As String) As String
        ' Looks for a file ending in _syn.txt, _fht.txt, _xt.xml, or _inspect.txt in folder strSourceFolderPath
        ' Returns the first matching file found

        Dim ioFolderInfo As System.IO.DirectoryInfo
        Dim ioFileInfo As System.IO.FileInfo
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

                ioFolderInfo = New System.IO.DirectoryInfo(strSourceFolderPath)
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

    Protected Function CleanupFilePaths(ByRef strInputFilePath As String, ByRef strOutputFolderPath As String) As Boolean
        ' Returns True if success, False if failure

        Dim ioFileInfo As System.IO.FileInfo
        Dim ioFolder As System.IO.DirectoryInfo

        Try
            ' Make sure strInputFilePath points to a valid file
            ioFileInfo = New System.IO.FileInfo(strInputFilePath)

            If Not ioFileInfo.Exists() Then
                SetErrorMessage("Input file not found: " & strInputFilePath)
                SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath)
                Return False
            Else
                If strOutputFolderPath Is Nothing OrElse strOutputFolderPath.Length = 0 Then
                    ' Define strOutputFolderPath based on strInputFilePath
                    strOutputFolderPath = ioFileInfo.DirectoryName
                End If

                ' Make sure strOutputFolderPath points to a folder
                ioFolder = New System.IO.DirectoryInfo(strOutputFolderPath)

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

    Public Shared Function DetermineResultsFileFormat(ByVal strFilePath As String) As ePeptideHitResultsFileFormatConstants
        ' Examine the extension on strFilePath to determine the file format

        If System.IO.Path.GetExtension(strFilePath).ToLower = ".xml" Then
            Return ePeptideHitResultsFileFormatConstants.XTandemXMLFile

        ElseIf System.IO.Path.GetFileNameWithoutExtension(strFilePath).ToLower.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE.ToLower) Then
            Return ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile

        ElseIf System.IO.Path.GetFileNameWithoutExtension(strFilePath).ToLower.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE.ToLower) Then
            Return ePeptideHitResultsFileFormatConstants.SequestSynopsisFile

        ElseIf System.IO.Path.GetFileNameWithoutExtension(strFilePath).ToLower.EndsWith(clsInSpecTResultsProcessor.FILENAME_SUFFIX_INSPECT_FILE.ToLower) Then
            Return ePeptideHitResultsFileFormatConstants.InSpectTXTFile

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

        mMassCorrectionTagsFilePath = String.Empty
        mModificationDefinitionsFilePath = String.Empty
        mSearchToolParameterFilePath = String.Empty

        mEnzymeMatchSpec = clsPeptideCleavageStateCalculator.GetDefaultEnzymeMatchSpec()

        mPeptideNTerminusMassChange = clsPeptideMassCalculator.DEFAULT_N_TERMINUS_MASS_CHANGE
        mPeptideCTerminusMassChange = clsPeptideMassCalculator.DEFAULT_C_TERMINUS_MASS_CHANGE

        ' Initialize mPeptideMods
        mPeptideMods = New clsPeptideModificationContainer

        ' Initialize mUniqueSequences
        mUniqueSequences = New clsUniqueSequencesContainer

        ' Initialize mSeqToProteinMap
        mSeqToProteinMap = New Hashtable

    End Sub

    Protected Function InitializeSequenceOutputFiles(ByVal strBaseOutputFilePath As String) As Boolean

        ' Initializes the StreamWriter objects using strBaseOutputFilePath as a base name and replacing the suffix with the default suffix names
        ' Returns True if success; does not catch errors; they will be thrown to the calling function if they occur

        Dim strResultToSeqMapFilePath As String
        Dim strSeqInfoFilePath As String
        Dim strModDetailsFilePath As String
        Dim strSeqToProteinMapFilePath As String

        ' Initialize the file paths based on strBaseOutputFilePath
        strResultToSeqMapFilePath = ReplaceFilenameSuffix(strBaseOutputFilePath, "", FILENAME_SUFFIX_RESULT_TO_SEQ_MAP)
        strSeqInfoFilePath = ReplaceFilenameSuffix(strBaseOutputFilePath, "", FILENAME_SUFFIX_SEQ_INFO)
        strModDetailsFilePath = ReplaceFilenameSuffix(strBaseOutputFilePath, "", FILENAME_SUFFIX_MOD_DETAILS)
        strSeqToProteinMapFilePath = ReplaceFilenameSuffix(strBaseOutputFilePath, "", FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP)


        ' Clear the unique sequences container
        mUniqueSequences.Clear()

        ' Clear the sequence to protein map
        mSeqToProteinMap.Clear()

        ' Initialize the ResultToSeqMap file
        mResultToSeqMapFile = New System.IO.StreamWriter(strResultToSeqMapFilePath)
        mResultToSeqMapFile.WriteLine("Result_ID" & SEP_CHAR & "Unique_Seq_ID")


        ' Initialize the SeqInfo file
        mSeqInfoFile = New System.IO.StreamWriter(strSeqInfoFilePath, False)
        mSeqInfoFile.WriteLine("Unique_Seq_ID" & SEP_CHAR & _
                               "Mod_Count" & SEP_CHAR & _
                               "Mod_Description" & SEP_CHAR & _
                               "Monoisotopic_Mass")

        ' Initialize the ModDetails file
        mModDetailsFile = New System.IO.StreamWriter(strModDetailsFilePath)
        mModDetailsFile.WriteLine("Unique_Seq_ID" & SEP_CHAR & _
                                  "Mass_Correction_Tag" & SEP_CHAR & _
                                  "Position")


        ' Initialize the SeqToProtein map file
        mSeqToProteinMapFile = New System.IO.StreamWriter(strSeqToProteinMapFilePath, False)
        mSeqToProteinMapFile.WriteLine("Unique_Seq_ID" & SEP_CHAR & _
                                       "Cleavage_State" & SEP_CHAR & _
                                       "Terminus_State" & SEP_CHAR & _
                                       "Protein_Name" & SEP_CHAR & _
                                       "Protein_Expectation_Value_Log(e)" & SEP_CHAR & _
                                       "Protein_Intensity_Log(I)")

        Return True

    End Function

    Public Shared Function IsNumber(ByVal strValue As String) As Boolean
        Dim objFormatProvider As System.Globalization.NumberFormatInfo
        Try
            Return Double.TryParse(strValue, Globalization.NumberStyles.Any, objFormatProvider, 0)
        Catch ex As Exception
            Return False
        End Try
    End Function

    Protected Overridable Function LoadParameterFileSettings(ByVal strParameterFilePath As String) As Boolean

        Const OPTIONS_SECTION As String = "PeptideHitResultsProcessorOptions"

        Dim objSettingsFile As New XmlSettingsFileAccessor

        Dim strLeftResidueRegEx As String, strRightResidueRegEx As String
        Dim blnValueNotPresent As Boolean

        Try

            If strParameterFilePath Is Nothing OrElse strParameterFilePath.Length = 0 Then
                ' No parameter file specified; nothing to load
                Return True
            End If

            If Not System.IO.File.Exists(strParameterFilePath) Then
                ' See if strParameterFilePath points to a file in the same directory as the application
                strParameterFilePath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location), System.IO.Path.GetFileName(strParameterFilePath))
                If Not System.IO.File.Exists(strParameterFilePath) Then
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

    Protected Sub OperationComplete()
        RaiseEvent ProgressComplete()
    End Sub

    Public Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String) As Boolean
        Return ProcessFile(strInputFilePath, strOutputFolderPath, String.Empty)
    End Function

    Public MustOverride Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean

    Protected Function ReplaceFilenameSuffix(ByVal strOriginalFilePath As String, ByVal strCurrentSuffixExpected As String, ByVal strNewSuffix As String) As String
        ' Examines strOriginalFilePath to see if it ends with strCurrentSuffixExpected (it normally should if you call this function)
        ' If found, it removes strCurrentSuffixExpected
        ' It then appends strNewSuffix and returns a full path to the file using the folder associated with strOriginalFilePath
        ' Note that strCurrentSuffixExpected and strNewSuffix may each contain a file extension though they do not have to
        '  If strNewSuffix does not contain an extension, then the path returned will end in the same extension as strOriginalFilePath

        Dim strNewFileName As String
        Dim strOriginalExtension As String

        ' Keep track of the original extension on strOriginalFilePath
        strOriginalExtension = System.IO.Path.GetExtension(strOriginalFilePath)

        ' Make sure strCurrentSuffixExpected and strNewSuffix are not nothing
        If strCurrentSuffixExpected Is Nothing Then strCurrentSuffixExpected = String.Empty
        If strNewSuffix Is Nothing Then strNewSuffix = String.Empty

        ' See if strCurrentSuffixExpected contains an extension
        If strCurrentSuffixExpected.Length > 0 AndAlso System.IO.Path.HasExtension(strCurrentSuffixExpected) Then
            strNewFileName = System.IO.Path.GetFileName(strOriginalFilePath)
        Else
            strNewFileName = System.IO.Path.GetFileNameWithoutExtension(strOriginalFilePath)
        End If

        ' If strNewFileName ends in strCurrentSuffixExpected, then remove strCurrentSuffixExpected from strNewFileName
        If strCurrentSuffixExpected.Length > 0 AndAlso strNewFileName.ToLower.EndsWith(strCurrentSuffixExpected.ToLower) Then
            strNewFileName = strNewFileName.Substring(0, strNewFileName.Length - strCurrentSuffixExpected.Length)
        Else
            strNewFileName = System.IO.Path.GetFileNameWithoutExtension(strOriginalFilePath)
        End If

        ' Append strNewSuffix to strNewFileName
        If System.IO.Path.HasExtension(strNewSuffix) Then
            strNewFileName &= strNewSuffix
        Else
            strNewFileName &= strNewSuffix & strOriginalExtension
        End If

        If System.IO.Path.IsPathRooted(strOriginalFilePath) Then
            strNewFileName = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(strOriginalFilePath), strNewFileName)
        End If

        Return strNewFileName
    End Function

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
            Else
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
            End If
        End If

        Return blnSuccess

    End Function

    Protected Sub ResetProgress()
        ResetProgress(String.Empty)
    End Sub

    Protected Sub ResetProgress(ByVal strProgressStepDescription As String)
        mProgressStepDescription = String.Copy(strProgressStepDescription)
        mProgressPercentComplete = 0
        RaiseEvent ProgressReset()
    End Sub

    Protected Sub SaveModificationSummaryFile(ByVal strModificationSummaryFilePath As String)
        Dim intIndex As Integer
        Dim swOutFile As System.IO.StreamWriter

        Try
            swOutFile = New System.IO.StreamWriter(strModificationSummaryFilePath, False)

            ' Write the header line
            swOutFile.WriteLine("Modification_Symbol" & SEP_CHAR & _
                                "Modification_Mass" & SEP_CHAR & _
                                "Target_Residues" & SEP_CHAR & _
                                "Modification_Type" & SEP_CHAR & _
                                "Mass_Correction_Tag" & SEP_CHAR & _
                                "Occurence_Count")

            For intIndex = 0 To mPeptideMods.ModificationCount - 1
                With mPeptideMods.GetModificationByIndex(intIndex)
                    If .OccurrenceCount > 0 OrElse Not .UnknownModAutoDefined Then
                        swOutFile.WriteLine(.ModificationSymbol & SEP_CHAR & _
                                            .ModificationMass.ToString & SEP_CHAR & _
                                            .TargetResidues & SEP_CHAR & _
                                            mPeptideMods.ModificationTypeToModificationSymbol(.ModificationType) & SEP_CHAR & _
                                            .MassCorrectionTag & SEP_CHAR & _
                                            .OccurrenceCount.ToString)
                    End If
                End With
            Next intIndex

        Catch ex As Exception
            Throw ex
        Finally
            If Not swOutFile Is Nothing Then
                swOutFile.Close()
            End If
        End Try

    End Sub

    Protected Sub SaveResultsFileEntrySeqInfo(ByRef objSearchResult As clsSearchResultsBaseClass, ByVal blnUpdateResultToSeqMapFile As Boolean)
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
                    mSeqInfoFile.WriteLine( _
                                    intUniqueSeqID.ToString & SEP_CHAR & _
                                    .SearchResultModificationCount & SEP_CHAR & _
                                    .PeptideModDescription & SEP_CHAR & _
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
                                mModDetailsFile.WriteLine( _
                                                intUniqueSeqID.ToString & SEP_CHAR & _
                                                .ModDefinition.MassCorrectionTag & SEP_CHAR & _
                                                .ResidueLocInPeptide.ToString)
                            End With
                        Next intIndex
                    End If
                End If
            End If

            ' Write a new entry to the SeqToProteinMap file if not yet defined
            If Not CheckSeqToProteinMapDefined(intUniqueSeqID, .ProteinName) Then
                mSeqToProteinMapFile.WriteLine( _
                                        intUniqueSeqID.ToString & SEP_CHAR & _
                                        CInt(.PeptideCleavageState).ToString & SEP_CHAR & _
                                        CInt(.PeptideTerminusState).ToString & SEP_CHAR & _
                                        .ProteinName & SEP_CHAR & _
                                        .ProteinExpectationValue & SEP_CHAR & _
                                        .ProteinIntensity)
            End If
        End With

    End Sub

    Protected Sub SetErrorCode(ByVal eNewErrorCode As ePHRPErrorCodes)
        SetErrorCode(eNewErrorCode, False)
    End Sub

    Protected Sub SetErrorCode(ByVal eNewErrorCode As ePHRPErrorCodes, ByVal blnLeaveExistingErrorCodeUnchanged As Boolean)
        If blnLeaveExistingErrorCodeUnchanged AndAlso mErrorCode <> ePHRPErrorCodes.NoError Then
            ' An error code is already defined; do not change it
        Else
            mErrorCode = eNewErrorCode
        End If
    End Sub

    Protected Sub SetErrorMessage(ByVal strMessage As String)
        If strMessage Is Nothing Then strMessage = String.Empty
        mErrorMessage = String.Copy(strMessage)
        If strMessage.Length > 0 Then
            RaiseEvent ErrorOccurred(mErrorMessage)
            Console.WriteLine(mErrorMessage)
        End If
    End Sub

    Protected Sub UpdateProgress(ByVal strProgressStepDescription As String)
        UpdateProgress(strProgressStepDescription, mProgressPercentComplete)
    End Sub

    Protected Sub UpdateProgress(ByVal sngPercentComplete As Single)
        UpdateProgress(Me.ProgressStepDescription, sngPercentComplete)
    End Sub

    Protected Sub UpdateProgress(ByVal strProgressStepDescription As String, ByVal sngPercentComplete As Single)
        mProgressStepDescription = String.Copy(strProgressStepDescription)
        If sngPercentComplete < 0 Then
            sngPercentComplete = 0
        ElseIf sngPercentComplete > 100 Then
            sngPercentComplete = 100
        End If
        mProgressPercentComplete = sngPercentComplete

        RaiseEvent ProgressChanged(Me.ProgressStepDescription, Me.ProgressPercentComplete)
    End Sub

    Protected Class ISearchOptionModificationInfoComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtSearchOptionModificationInfoType = CType(x, udtSearchOptionModificationInfoType)
            Dim yData As udtSearchOptionModificationInfoType = CType(y, udtSearchOptionModificationInfoType)

            If xData.SortOrder > yData.SortOrder Then
                Return 1
            ElseIf xData.SortOrder < yData.SortOrder Then
                Return -1
            Else
                If xData.ModificationMass > yData.ModificationMass Then
                    Return 1
                ElseIf xData.ModificationMass < yData.ModificationMass Then
                    Return -1
                Else
                    Return 0
                End If
            End If
        End Function
    End Class

    Friend Class IModNameAndResidueLocComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtModNameAndResidueLocType = CType(x, udtModNameAndResidueLocType)
            Dim yData As udtModNameAndResidueLocType = CType(y, udtModNameAndResidueLocType)

            If xData.ResidueLocInPeptide > yData.ResidueLocInPeptide Then
                Return 1
            ElseIf xData.ResidueLocInPeptide < yData.ResidueLocInPeptide Then
                Return -1
            Else
                If xData.ModName Is Nothing Then xData.ModName = String.Empty
                If yData.ModName Is Nothing Then yData.ModName = String.Empty

                If xData.ModName > yData.ModName Then
                    Return 1
                ElseIf xData.ModName < yData.ModName Then
                    Return -1
                Else
                    Return 0
                End If
            End If
        End Function
    End Class
End Class
