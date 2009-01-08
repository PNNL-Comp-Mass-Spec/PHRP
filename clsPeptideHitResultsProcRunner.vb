Option Strict On
Imports PeptideHitResultsProcessor.clsPHRPBaseClass

' This class calls clsSequestSynopsisFileProcessor or clsXTandemResultsConverter
' to process the files to determine the modifications present for each peptide,
' along with other information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 6, 2006
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

Public Class clsPeptideHitResultsProcRunner
    Inherits clsProcessFilesBaseClass

    Public Sub New()
        MyBase.mFileDate = "December 9, 2008"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

    Protected Const DMS_CONNECTION_STRING_DEFAULT As String = "Data Source=gigasax;Initial Catalog=DMS5_T3;User=dmsreader;Password=dms4fun"

    ' Error codes specialized for this class
    Public Enum eResultsProcessorErrorCodes As Integer
        NoError = 0
        UnspecifiedError = -1
    End Enum

#End Region

#Region "Structures"

#End Region

#Region "Classwide Variables"
    Protected mPeptideHitResultsFileFormat As ePeptideHitResultsFileFormatConstants

    Protected mObtainModificationDefinitionsFromDMS As Boolean
    Protected mDMSConnectionString As String

    Protected mMassCorrectionTagsFilePath As String
    Protected mModificationDefinitionsFilePath As String
    Protected mSearchToolParameterFilePath As String
    Protected mInspectSynopsisFilePValueThreshold As Single

    Protected WithEvents mPeptideHitResultsProcessor As PeptideHitResultsProcessor.clsPHRPBaseClass

    Protected mWarnMissingParameterFileSection As Boolean
    Protected mLocalErrorCode As eResultsProcessorErrorCodes
#End Region

#Region "Properties"

    Public Property InspectSynopsisFilePValueThreshold() As Single
        Get
            Return mInspectSynopsisFilePValueThreshold
        End Get
        Set(ByVal value As Single)
            mInspectSynopsisFilePValueThreshold = value
        End Set
    End Property
    Public ReadOnly Property LocalErrorCode() As eResultsProcessorErrorCodes
        Get
            Return mLocalErrorCode
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

    Public Property PeptideHitResultsFileFormat() As ePeptideHitResultsFileFormatConstants
        Get
            Return mPeptideHitResultsFileFormat
        End Get
        Set(ByVal Value As ePeptideHitResultsFileFormatConstants)
            mPeptideHitResultsFileFormat = Value
        End Set
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

    Public Overrides Function GetDefaultExtensionsToParse() As String()
        Dim strExtensionsToParse(1) As String

        ' Note: If mPeptideHitResultsFileFormat = .AutoDetermine
        '  then this class will only parse .txt files if they match 
        ' PeptideHitResultsProcessor.clsSequestResultsProcessor.SEQUEST_FIRST_HITS_FILE_SUFFIX or
        ' PeptideHitResultsProcessor.clsSequestResultsProcessor.SEQUEST_SYNOPSIS_FILE_SUFFIX

        strExtensionsToParse(0) = ".txt"
        strExtensionsToParse(1) = ".xml"

        Return strExtensionsToParse

    End Function

    Public Overrides Function GetErrorMessage() As String
        ' Returns "" if no error

        Dim strErrorMessage As String

        If MyBase.ErrorCode = clsProcessFilesBaseClass.eProcessFilesErrorCodes.LocalizedError Or _
           MyBase.ErrorCode = clsProcessFilesBaseClass.eProcessFilesErrorCodes.NoError Then
            Select Case mLocalErrorCode
                Case eResultsProcessorErrorCodes.NoError
                    strErrorMessage = ""
                Case eResultsProcessorErrorCodes.UnspecifiedError
                    strErrorMessage = "Unspecified localized error"
                Case Else
                    ' This shouldn't happen
                    strErrorMessage = "Unknown error state"
            End Select
        Else
            strErrorMessage = MyBase.GetBaseClassErrorMessage()
        End If

        Return strErrorMessage
    End Function

    Private Sub InitializeLocalVariables()
        MyBase.ShowMessages = False

        mPeptideHitResultsFileFormat = ePeptideHitResultsFileFormatConstants.AutoDetermine

        mMassCorrectionTagsFilePath = String.Empty
        mModificationDefinitionsFilePath = String.Empty
        mSearchToolParameterFilePath = String.Empty

        mInspectSynopsisFilePValueThreshold = PeptideHitResultsProcessor.clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD

        mWarnMissingParameterFileSection = True

        mLocalErrorCode = eResultsProcessorErrorCodes.NoError

    End Sub

    Private Sub LoadModificationInfoFromDMS()

        Dim blnSuccess As Boolean

        blnSuccess = False

        ' ToDo: Contact DMS to get the modification information
        ' The results of this query will need to be filtered to get the info for just this analysis job

        'SELECT D.Dataset_Num, 
        '    AJ.AJ_jobID, PFMI.Local_Symbol, 
        '    PFMI.Monoisotopic_Mass_Correction, PFMI.Residue_Symbol, 
        '    PFMI.Mod_Type_Symbol, PFMI.Mass_Correction_Tag
        'FROM dbo.T_Analysis_Job AJ INNER JOIN
        '    dbo.T_Dataset D ON 
        '    AJ.AJ_datasetID = D.Dataset_ID LEFT
        '     OUTER JOIN
        '    dbo.V_Param_File_Mass_Mod_Info PFMI ON 
        '    AJ.AJ_parmFileName = PFMI.Param_File_Name
        'WHERE (D.Dataset_Num = 'QC_05_2_a_24Oct05_Doc_0508-08')
        'ORDER BY AJ.AJ_jobID, PFMI.Local_Symbol

        'SELECT D.Dataset_Num, 
        '    AJ.AJ_jobID, PFMI.Local_Symbol, 
        '    PFMI.Monoisotopic_Mass_Correction, PFMI.Residue_Symbol, 
        '    PFMI.Mod_Type_Symbol, PFMI.Mass_Correction_Tag
        'FROM dbo.T_Analysis_Job AJ INNER JOIN
        '    dbo.T_Dataset D ON 
        '    AJ.AJ_datasetID = D.Dataset_ID LEFT
        '     OUTER JOIN
        '    dbo.V_Param_File_Mass_Mod_Info PFMI ON 
        '    AJ.AJ_parmFileName = PFMI.Param_File_Name
        'WHERE (AJ.AJ_jobID = 47703)
        'ORDER BY AJ.AJ_jobID, PFMI.Local_Symbol

    End Sub

    Private Function LoadParameterFileSettings(ByVal strParameterFilePath As String) As Boolean

        Const OPTIONS_SECTION As String = "PeptideHitResultsProcRunner"

        Dim objSettingsFile As New XmlSettingsFileAccessor

        Dim intValue As Integer

        Try

            If strParameterFilePath Is Nothing OrElse strParameterFilePath.Length = 0 Then
                ' No parameter file specified; nothing to load
                Return True
            End If

            If Not System.IO.File.Exists(strParameterFilePath) Then
                ' See if strParameterFilePath points to a file in the same directory as the application
                strParameterFilePath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location), System.IO.Path.GetFileName(strParameterFilePath))
                If Not System.IO.File.Exists(strParameterFilePath) Then
                    MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.ParameterFileNotFound)
                    Return False
                End If
            End If

            If objSettingsFile.LoadSettings(strParameterFilePath) Then
                If Not objSettingsFile.SectionPresent(OPTIONS_SECTION) Then
                    ShowErrorMessage("The node '<section name=""" & OPTIONS_SECTION & """> was not found in the parameter file: " & strParameterFilePath)
                    MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.InvalidParameterFile)
                    Return False
                Else
                    mDMSConnectionString = objSettingsFile.GetParam(OPTIONS_SECTION, "DMSConnectionString", DMS_CONNECTION_STRING_DEFAULT)
                    mObtainModificationDefinitionsFromDMS = objSettingsFile.GetParam(OPTIONS_SECTION, "ObtainModificationDefinitionsFromDMS", mObtainModificationDefinitionsFromDMS)

                    intValue = objSettingsFile.GetParam(OPTIONS_SECTION, "PeptideHitResultsFileFormat", CInt(mPeptideHitResultsFileFormat))
                    Try
                        mPeptideHitResultsFileFormat = CType(intValue, ePeptideHitResultsFileFormatConstants)
                    Catch ex As Exception
                        mPeptideHitResultsFileFormat = ePeptideHitResultsFileFormatConstants.AutoDetermine
                    End Try
                End If
            End If

        Catch ex As Exception
            HandleException("Error in LoadParameterFileSettings", ex)
            Return False
        End Try

        Return True

    End Function

    ' Main processing function
    Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String, ByVal blnResetErrorCode As Boolean) As Boolean
        ' Returns True if success, False if failure

        Dim strStatusMessage As String
        Dim strMessage As String

        Dim ePeptideHitResultsFormat As ePeptideHitResultsFileFormatConstants

        Dim blnSuccess As Boolean

        If blnResetErrorCode Then
            SetLocalErrorCode(eResultsProcessorErrorCodes.NoError)
        End If

        If Not LoadParameterFileSettings(strParameterFilePath) Then
            strStatusMessage = "Parameter file load error: " & strParameterFilePath
            ShowErrorMessage(strStatusMessage)

            If MyBase.ErrorCode = clsProcessFilesBaseClass.eProcessFilesErrorCodes.NoError Then
                MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.InvalidParameterFile)
            End If
            Return False
        End If

        Try
            If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
                ShowMessage("Input file name is empty")
                MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.InvalidInputFilePath)
            Else
                ' Note that CleanupFilePaths() will update mOutputFolderPath, which is used by LogMessage()
                If Not CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                    MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.FilePathError)
                Else
                    MyBase.mProgressStepDescription = "Parsing " & System.IO.Path.GetFileName(strInputFilePath)
                    LogMessage(MyBase.mProgressStepDescription)
                    MyBase.ResetProgress()

                    Try
                        If mPeptideHitResultsFileFormat = ePeptideHitResultsFileFormatConstants.AutoDetermine Then
                            ePeptideHitResultsFormat = PeptideHitResultsProcessor.clsPHRPBaseClass.DetermineResultsFileFormat(strInputFilePath)
                        Else
                            ePeptideHitResultsFormat = mPeptideHitResultsFileFormat
                        End If

                        If ePeptideHitResultsFormat = ePeptideHitResultsFileFormatConstants.AutoDetermine Then
                            ' If ePeptideHitResultsFormat is still AutoDetermine that means we couldn't figure out the format
                            ' Return True if strInputFilePath ends in .txt, otherwise return false
                            ' Determine the file type and call the appropriate class
                            If System.IO.Path.GetExtension(strInputFilePath).ToLower = ".txt" Then
                                ' General text file; we won't process it but we will set blnSuccess = True
                                blnSuccess = True
                            Else
                                blnSuccess = False
                            End If

                            strMessage = "Warning: Could not determine the format of the input file.  It must end in " & PeptideHitResultsProcessor.clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE & ".txt, " & PeptideHitResultsProcessor.clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE & ".txt, or .xml"
                            ShowMessage(strMessage)
                        Else
                            Select Case ePeptideHitResultsFormat
                                Case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile
                                    mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsSequestResultsProcessor

                                Case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile
                                    mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsSequestResultsProcessor

                                Case ePeptideHitResultsFileFormatConstants.XTandemXMLFile
                                    mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsXTandemResultsProcessor

                                Case ePeptideHitResultsFileFormatConstants.InSpectTXTFile
                                    mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsInSpecTResultsProcessor

                                Case Else
                                    ' Unknown format
                                    blnSuccess = False
                            End Select

                            If Not mPeptideHitResultsProcessor Is Nothing Then
                                If mObtainModificationDefinitionsFromDMS Then
                                    LoadModificationInfoFromDMS()
                                End If
                                With mPeptideHitResultsProcessor
                                    .MassCorrectionTagsFilePath = mMassCorrectionTagsFilePath
                                    .ModificationDefinitionsFilePath = mModificationDefinitionsFilePath
                                    .SearchToolParameterFilePath = mSearchToolParameterFilePath
                                    .InspectSynopsisFilePValueThreshold = mInspectSynopsisFilePValueThreshold

                                    .WarnMissingParameterFileSection = mWarnMissingParameterFileSection
                                End With
                                blnSuccess = mPeptideHitResultsProcessor.ProcessFile(strInputFilePath, strOutputFolderPath, strParameterFilePath)
                                If Not blnSuccess Then
                                    ShowErrorMessage(mPeptideHitResultsProcessor.ErrorMessage)
                                Else
                                    LogMessage("Processing Complete")
                                    OperationComplete()
                                End If
                            End If

                        End If

                    Catch ex As Exception
                        HandleException("Error calling ParseXTandemResultsFile", ex)
                    End Try
                End If
            End If
        Catch ex As Exception
            HandleException("Error in ProcessFile", ex)
        End Try

        Return blnSuccess

    End Function

    Private Sub SetLocalErrorCode(ByVal eNewErrorCode As eResultsProcessorErrorCodes)
        SetLocalErrorCode(eNewErrorCode, False)
    End Sub

    Private Sub SetLocalErrorCode(ByVal eNewErrorCode As eResultsProcessorErrorCodes, ByVal blnLeaveExistingErrorCodeUnchanged As Boolean)

        If blnLeaveExistingErrorCodeUnchanged AndAlso mLocalErrorCode <> eResultsProcessorErrorCodes.NoError Then
            ' An error code is already defined; do not change it
        Else
            mLocalErrorCode = eNewErrorCode

            If eNewErrorCode = eResultsProcessorErrorCodes.NoError Then
                If MyBase.ErrorCode = clsProcessFilesBaseClass.eProcessFilesErrorCodes.LocalizedError Then
                    MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.NoError)
                End If
            Else
                MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.LocalizedError)
            End If
        End If

    End Sub

    Private Sub mPeptideHitResultsProcessor_ErrorOccurred(ByVal ErrorMessage As String) Handles mPeptideHitResultsProcessor.ErrorOccurred
        LogMessage(ErrorMessage, eMessageTypeConstants.ErrorMsg)
    End Sub

    Private Sub mPeptideHitResultsProcessor_ProgressChanged(ByVal taskDescription As String, ByVal percentComplete As Single) Handles mPeptideHitResultsProcessor.ProgressChanged
        UpdateProgress(taskDescription, percentComplete)
    End Sub

End Class
