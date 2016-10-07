Option Strict On

' This class calls clsSequestSynopsisFileProcessor or clsXTandemResultsConverter
' to process the files to determine the modifications present for each peptide,
' along with other information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 6, 2006
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

Imports PeptideHitResultsProcessor.clsPHRPBaseClass
Imports System.IO

Public Class clsPeptideHitResultsProcRunner
    Inherits clsProcessFilesBaseClass

    Public Sub New()
        MyBase.mFileDate = PROGRAM_DATE
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

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

    Protected mUseExistingMTSPepToProteinMapFile As Boolean

    Protected WithEvents mPeptideHitResultsProcessor As PeptideHitResultsProcessor.clsPHRPBaseClass

    Protected mLocalErrorCode As eResultsProcessorErrorCodes
#End Region

#Region "Properties"

    Public Property CreateInspectOrMSGFDbFirstHitsFile As Boolean = True

    Public Property CreateInspectOrMSGFDbSynopsisFile As Boolean = True

    Public Property CreateProteinModsFile As Boolean

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks>
    ''' Setting this to true assumes the input file is a valid PHRP data file
    ''' Consequently, the code will only try to create the _ProteinMods.txt file, it will not re-create the PHRP data files
    ''' When this is True, then mCreateProteinModsFile is assumed to be true
    ''' </remarks>
    Public Property CreateProteinModsUsingPHRPDataFile As Boolean

    Public Property FastaFilePath As String

    Public Property IgnorePeptideToProteinMapperErrors As Boolean

    Public Property InspectSynopsisFilePValueThreshold As Single

    Public ReadOnly Property LocalErrorCode() As eResultsProcessorErrorCodes
        Get
            Return mLocalErrorCode
        End Get
    End Property

    Public Property MassCorrectionTagsFilePath As String

    Public Property ModificationDefinitionsFilePath As String

    Public Property MODaMODPlusSynopsisFileProbabilityThreshold As Single

    Public Property MsgfPlusEValueThreshold As Single

    Public Property MsgfPlusSpecEValueThreshold As Single

    Public Property ProteinModsFileIncludesReversedProteins As Boolean

    Public Property SearchToolParameterFilePath As String

    ''' <summary>
    ''' 
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks>If this is true and the _PepToProtMap.txt file isn't found then it will be created using the the Fasta file specified by mFastaFilePath</remarks>
    Public Property UseExistingMTSPepToProteinMapFile As Boolean

    Public Property WarnMissingParameterFileSection As Boolean

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

        If MyBase.ErrorCode = clsProcessFilesBaseClass.eProcessFilesErrorCodes.LocalizedError Or
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

        MassCorrectionTagsFilePath = String.Empty
        ModificationDefinitionsFilePath = String.Empty
        SearchToolParameterFilePath = String.Empty

        CreateProteinModsFile = False
        FastaFilePath = String.Empty
        IgnorePeptideToProteinMapperErrors = False
        ProteinModsFileIncludesReversedProteins = False
        mUseExistingMTSPepToProteinMapFile = False

        CreateInspectOrMSGFDbFirstHitsFile = False
        CreateInspectOrMSGFDbSynopsisFile = False
        InspectSynopsisFilePValueThreshold = PeptideHitResultsProcessor.clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD

        MODaMODPlusSynopsisFileProbabilityThreshold = PeptideHitResultsProcessor.clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD

        MsgfPlusEValueThreshold = PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD
        MsgfPlusSpecEValueThreshold = PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD

        WarnMissingParameterFileSection = True

        mLocalErrorCode = eResultsProcessorErrorCodes.NoError

    End Sub

    Private Sub InitializePeptideHitResultsProcessor(strInputFilePath As String)

        If mObtainModificationDefinitionsFromDMS Then
            LoadModificationInfoFromDMS()
        End If

        Dim fiSourceFile As FileInfo
        fiSourceFile = New FileInfo(strInputFilePath)

        With mPeptideHitResultsProcessor
            .MassCorrectionTagsFilePath = ResolveFilePath(fiSourceFile.DirectoryName, MassCorrectionTagsFilePath)
            .ModificationDefinitionsFilePath = ResolveFilePath(fiSourceFile.DirectoryName, ModificationDefinitionsFilePath)
            .SearchToolParameterFilePath = ResolveFilePath(fiSourceFile.DirectoryName, SearchToolParameterFilePath)

            .CreateProteinModsFile = CreateProteinModsFile
            .FastaFilePath = FastaFilePath
            .IgnorePeptideToProteinMapperErrors = IgnorePeptideToProteinMapperErrors
            .ProteinModsFileIncludesReversedProteins = ProteinModsFileIncludesReversedProteins
            .UseExistingMTSPepToProteinMapFile = mUseExistingMTSPepToProteinMapFile

            .CreateInspectFirstHitsFile = CreateInspectOrMSGFDbFirstHitsFile
            .CreateInspectSynopsisFile = CreateInspectOrMSGFDbSynopsisFile
            .InspectSynopsisFilePValueThreshold = InspectSynopsisFilePValueThreshold

            .MODaMODPlusSynopsisFileProbabilityThreshold = MODaMODPlusSynopsisFileProbabilityThreshold

            .MSGFDBSynopsisFileEValueThreshold = MsgfPlusEValueThreshold
            .MSGFDBSynopsisFileSpecEValueThreshold = MsgfPlusSpecEValueThreshold

            .WarnMissingParameterFileSection = WarnMissingParameterFileSection
        End With
    End Sub

    Private Sub LoadModificationInfoFromDMS()

        Dim blnSuccess As Boolean = False

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

    Private Function LoadParameterFileSettings(strParameterFilePath As String) As Boolean

        Const OPTIONS_SECTION As String = "PeptideHitResultsProcRunner"

        Dim objSettingsFile As New XmlSettingsFileAccessor

        Dim intValue As Integer

        Try

            If String.IsNullOrWhiteSpace(strParameterFilePath) Then
                ' No parameter file specified; nothing to load
                Return True
            End If

            If Not File.Exists(strParameterFilePath) Then
                ' See if strParameterFilePath points to a file in the same directory as the application
                strParameterFilePath = Path.Combine(Path.GetDirectoryName(Reflection.Assembly.GetExecutingAssembly().Location), Path.GetFileName(strParameterFilePath))
                If Not File.Exists(strParameterFilePath) Then
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
    Public Overloads Overrides Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String, blnResetErrorCode As Boolean) As Boolean
        ' Returns True if success, False if failure

        Dim strStatusMessage As String
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
            If String.IsNullOrWhiteSpace(strInputFilePath) Then
                ShowErrorMessage("Input file name is empty")
                MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.InvalidInputFilePath)
            Else
                ' Note that CleanupFilePaths() will update mOutputFolderPath, which is used by LogMessage()
                If Not CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                    MyBase.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.FilePathError)
                Else
                    MyBase.mProgressStepDescription = "Parsing " & Path.GetFileName(strInputFilePath)
                    LogMessage(MyBase.mProgressStepDescription)
                    MyBase.ResetProgress()

                    If CreateProteinModsUsingPHRPDataFile Then
                        blnSuccess = StartCreateProteinModsViaPHRPData(strInputFilePath, strOutputFolderPath)
                    Else
                        blnSuccess = StartPHRP(strInputFilePath, strOutputFolderPath, strParameterFilePath)
                    End If

                End If
            End If
        Catch ex As Exception
            HandleException("Error in ProcessFile", ex)
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Looks for file strFileNameOrPath in the current working directory
    ''' If not found, then looks in strSourceFolderPath
    ''' </summary>
    ''' <param name="strSourceFolderPath">Path to the folder containing the input file</param>
    ''' <param name="strFileNameOrPath">File to find (either filename or full file path)</param>
    ''' <returns>The path to the file if found, or strFileNameOrPath if not found</returns>
    ''' <remarks></remarks>
    Protected Function ResolveFilePath(strSourceFolderPath As String, strFileNameOrPath As String) As String

        If File.Exists(strFileNameOrPath) Then
            Return strFileNameOrPath
        Else
            Dim strNewPath As String
            strNewPath = Path.Combine(strSourceFolderPath, Path.GetFileName(strFileNameOrPath))
            If File.Exists(strNewPath) Then
                Return strNewPath
            End If
        End If

        Return strFileNameOrPath
    End Function

    Private Sub SetLocalErrorCode(eNewErrorCode As eResultsProcessorErrorCodes)
        SetLocalErrorCode(eNewErrorCode, False)
    End Sub

    Private Sub SetLocalErrorCode(eNewErrorCode As eResultsProcessorErrorCodes, blnLeaveExistingErrorCodeUnchanged As Boolean)

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

    Private Function StartCreateProteinModsViaPHRPData(strInputFilePath As String, strOutputFolderPath As String) As Boolean

        Dim strMessage As String
        Dim ePeptideHitResultType As PHRPReader.clsPHRPReader.ePeptideHitResultType

        Dim blnSuccess As Boolean

        Try
            Select Case mPeptideHitResultsFileFormat
                Case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.Sequest
                    LogMessage("Detected SEQUEST First Hits file")

                Case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.Sequest
                    LogMessage("Detected SEQUEST Synopsis file")

                Case ePeptideHitResultsFileFormatConstants.XTandemXMLFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.XTandem
                    LogMessage("Detected X!Tandem XML file")

                Case ePeptideHitResultsFileFormatConstants.InSpectTXTFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.Inspect
                    LogMessage("Detected Inspect results file")

                Case ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.MSGFDB
                    LogMessage("Detected MSGFDB (or MSGF+) results file")

                Case ePeptideHitResultsFileFormatConstants.MSAlignTXTFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.MSAlign
                    LogMessage("Detected MSAlign results file")

                Case ePeptideHitResultsFileFormatConstants.MODPlusTXTFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.MODPlus
                    LogMessage("Detected MODPlus results file")

                Case ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.MSPathFinder
                    LogMessage("Detected MSPathfinder results file")


                Case Else
                    ' Includes ePeptideHitResultsFileFormatConstants.AutoDetermine
                    ePeptideHitResultType = PHRPReader.clsPHRPReader.AutoDetermineResultType(strInputFilePath)
            End Select

            If ePeptideHitResultType = PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
                blnSuccess = False

                strMessage = "Error: Could not determine the format of the PHRP data file: " & strInputFilePath
                ShowErrorMessage(strMessage)
            Else
                Select Case ePeptideHitResultType
                    Case PHRPReader.clsPHRPReader.ePeptideHitResultType.Sequest
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsSequestResultsProcessor

                    Case PHRPReader.clsPHRPReader.ePeptideHitResultType.XTandem
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsXTandemResultsProcessor

                    Case PHRPReader.clsPHRPReader.ePeptideHitResultType.Inspect
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsInSpecTResultsProcessor

                    Case PHRPReader.clsPHRPReader.ePeptideHitResultType.MSGFDB
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMSGFDBResultsProcessor

                    Case PHRPReader.clsPHRPReader.ePeptideHitResultType.MSAlign
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMSAlignResultsProcessor

                    Case PHRPReader.clsPHRPReader.ePeptideHitResultType.MODa
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMODaResultsProcessor

                    Case Else
                        ' Unknown format
                        strMessage = "Error: Unrecognized value for ePeptideHitResultType: " & ePeptideHitResultType.ToString()
                        ShowErrorMessage(strMessage)
                        blnSuccess = False
                End Select

                If Not mPeptideHitResultsProcessor Is Nothing Then
                    InitializePeptideHitResultsProcessor(strInputFilePath)

                    blnSuccess = mPeptideHitResultsProcessor.CreateProteinModDetailsFile(strInputFilePath, strOutputFolderPath)

                    If Not blnSuccess Then
                        ShowErrorMessage(mPeptideHitResultsProcessor.ErrorMessage)
                    Else
                        LogMessage("Processing Complete")
                        OperationComplete()
                    End If
                End If

            End If

        Catch ex As Exception
            HandleException("Error calling CreateProteinModDetailsFile in CreateProteinModsViaPHRPData", ex)
        End Try

        Return blnSuccess

    End Function

    Private Function StartPHRP(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean

        Dim strMessage As String
        Dim ePeptideHitResultsFormat As ePeptideHitResultsFileFormatConstants

        Dim blnSuccess As Boolean

        Try
            If mPeptideHitResultsFileFormat = ePeptideHitResultsFileFormatConstants.AutoDetermine Then
                ePeptideHitResultsFormat = PeptideHitResultsProcessor.clsPHRPBaseClass.DetermineResultsFileFormat(strInputFilePath)
            Else
                ePeptideHitResultsFormat = mPeptideHitResultsFileFormat
            End If

            If ePeptideHitResultsFormat = ePeptideHitResultsFileFormatConstants.AutoDetermine Then
                ' If ePeptideHitResultsFormat is still AutoDetermine that means we couldn't figure out the format
                blnSuccess = False

                strMessage = "Error: Could not determine the format of the input file.  It must end in " &
                  PeptideHitResultsProcessor.clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE & ".txt, " &
                  PeptideHitResultsProcessor.clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE & ".txt, .xml (for X!Tandem), " &
                  PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE & ".txt, " &
                  PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE & ".tsv, " &
                  PeptideHitResultsProcessor.clsMSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE & ".txt, " &
                  PeptideHitResultsProcessor.clsMODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE & ".txt, " &
                  PeptideHitResultsProcessor.clsMODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE & ".txt, or" &
                  PeptideHitResultsProcessor.clsMSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE & ".tsv"

                ShowErrorMessage(strMessage)
            Else

                Select Case ePeptideHitResultsFormat
                    Case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsSequestResultsProcessor
                        LogMessage("Detected SEQUEST First Hits file")

                    Case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsSequestResultsProcessor
                        LogMessage("Detected SEQUEST Synopsis file")

                    Case ePeptideHitResultsFileFormatConstants.XTandemXMLFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsXTandemResultsProcessor
                        LogMessage("Detected X!Tandem XML file")

                    Case ePeptideHitResultsFileFormatConstants.InSpectTXTFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsInSpecTResultsProcessor
                        LogMessage("Detected Inspect results file")

                    Case ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMSGFDBResultsProcessor
                        LogMessage("Detected MSGFDB (or MSGF+) results file")

                    Case ePeptideHitResultsFileFormatConstants.MSAlignTXTFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMSAlignResultsProcessor
                        LogMessage("Detected MSAlign results file")

                    Case ePeptideHitResultsFileFormatConstants.MODaTXTFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMODaResultsProcessor
                        LogMessage("Detected MODa results file")

                    Case ePeptideHitResultsFileFormatConstants.MODPlusTXTFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMODPlusResultsProcessor
                        LogMessage("Detected MODPlus results file")

                    Case ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile
                        mPeptideHitResultsProcessor = New PeptideHitResultsProcessor.clsMSPathFinderResultsProcessor
                        LogMessage("Detected MSPathFinder results file")
                    Case Else
                        ' Unknown format
                        blnSuccess = False
                End Select

                If Not mPeptideHitResultsProcessor Is Nothing Then
                    InitializePeptideHitResultsProcessor(strInputFilePath)

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
            HandleException("Error calling ProcessFile in StartPHRP", ex)
        End Try

        Return blnSuccess

    End Function

    Private Sub mPeptideHitResultsProcessor_ErrorOccurred(ErrorMessage As String) Handles mPeptideHitResultsProcessor.ErrorOccurred
        LogMessage(ErrorMessage, eMessageTypeConstants.ErrorMsg)
    End Sub

    Private Sub mPeptideHitResultsProcessor_MessageEvent(message As String) Handles mPeptideHitResultsProcessor.MessageEvent
        LogMessage(message)
    End Sub

    Private Sub mPeptideHitResultsProcessor_ProgressChanged(taskDescription As String, percentComplete As Single) Handles mPeptideHitResultsProcessor.ProgressChanged
        UpdateProgress(taskDescription, percentComplete)
    End Sub

    Private Sub mPeptideHitResultsProcessor_WarningMessageEvent(WarningMessage As String) Handles mPeptideHitResultsProcessor.WarningMessageEvent
        LogMessage(WarningMessage, eMessageTypeConstants.Warning)
    End Sub

End Class
