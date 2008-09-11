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
        MyBase.mFileDate = "August 19, 2008"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
    Public Const FILENAME_SUFFIX_INSPECT_FILE As String = "_inspect"

    Public Const TERMINUS_SYMBOL_INSPECT_B As String = "*."
    Public Const TERMINUS_SYMBOL_INSPECT_E As String = ".*"
    Private Const FSCORE_THRESHOLD As Single = 0
    Private Const TOTALPRMSCORE_THRESHOLD As Single = 50

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Public Enum eInspectFileColumns As Integer
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

#End Region

#Region "Structures"
#End Region

#Region "Classwide Variables"
    Protected mNextResultID As Integer
    Protected m_tmpSearchTable As New DataTable
#End Region

#Region "Properties"
#End Region

    Private Sub InitializeLocalVariables()
        ' Note: This function is called from ParseInSpecTResultsFile()
        ' These variables will therefore be reset for each InSpecT TXT file analyzed
        mNextResultID = 1

    End Sub


    Protected Function ParseInSpecTResultsFile(ByVal strInputFilePath As String, ByVal strOutputFilePath As String, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

        Dim blnSuccess As Boolean

        Try

            'Just a placeholder


        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Function CreateFHTResultsFile(ByVal strInputFilePath As String, ByVal strOutputFilePath As String, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' This routine creates a first hits file from the output from InSpecT
        ' The first hits files writes the record with the highest TotalPRMScore

        Dim srDataFile As System.IO.StreamReader
        Dim srResultFile As System.IO.StreamWriter

        Dim strPreviousScan As String
        Dim strPreviousScanState As String

        Dim strLineIn As String
        Dim protein As String = ""

        Dim objSearchResult As clsSearchResultsInSpecT
        Dim tmpObjSearchResult As clsSearchResultsInSpecT

        Dim intResultsProcessed As Integer
        Dim intResultID As Integer = 0

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean

        Dim strErrorLog As String

        Try
            ' Initialize objSearchResult
            objSearchResult = New clsSearchResultsInSpecT(mPeptideMods)
            tmpObjSearchResult = New clsSearchResultsInSpecT(mPeptideMods)

            ' Initialize variables
            strPreviousScan = String.Empty
            strPreviousScanState = String.Empty

            Try
                ' Open the input file and parse it
                ' Initialize the stream reader and the stream Text writer
                srDataFile = New System.IO.StreamReader(strInputFilePath)

                srResultFile = New System.IO.StreamWriter(strOutputFilePath)

                strErrorLog = String.Empty
                intResultsProcessed = 0

                'Create temporary table to hold records
                createTmpTable(strErrorLog)

                Dim firstDataRec As Boolean = True
                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing
                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then
                        blnValidSearchResult = ParseInSpecTResultsFileEntry(strLineIn, objSearchResult, strErrorLog, intResultsProcessed)
                        If Not blnValidSearchResult Then
                            'assume this is the first line, so write it as the header
                            WriteFileHeader(strLineIn, srResultFile, objSearchResult, strErrorLog)
                        Else
                            If firstDataRec Then
                                saveCurrentRecordIntoTempTable(objSearchResult, strErrorLog)
                            ElseIf strPreviousScan <> objSearchResult.Scan And Not firstDataRec Then
                                determineBestRecord(srResultFile, intResultID, strErrorLog)
                                m_tmpSearchTable.Clear()
                                saveCurrentRecordIntoTempTable(objSearchResult, strErrorLog)
                            Else
                                saveCurrentRecordIntoTempTable(objSearchResult, strErrorLog)
                            End If
                            strPreviousScan = objSearchResult.Scan
                            firstDataRec = False
                        End If
                        'check to see if there is any data left
                        'if not, then write the last record
                        If srDataFile.Peek = -1 Then
                            'write the record
                            determineBestRecord(srResultFile, intResultID, strErrorLog)
                        End If
                        ' Update the progress
                        UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))
                        intResultsProcessed += 1
                    End If
                Loop

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
                If Not srResultFile Is Nothing Then
                    srResultFile.Close()
                    srResultFile = Nothing
                End If
            End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Function CreateSYNResultsFile(ByVal strInputFilePath As String, ByVal strOutputFilePath As String, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' This routine creates a synopsis file using the output from InSpecT
        ' The synopsis files writes every record with p-value below a set threshold

        Dim srDataFile As System.IO.StreamReader
        Dim srResultFile As System.IO.StreamWriter

        Dim strPreviousScan As String

        Dim strLineIn As String
        Dim protein As String = ""

        Dim objSearchResult As clsSearchResultsInSpecT

        Dim intResultsProcessed As Integer
        Dim intResultID As Integer = 0

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean

        Dim strErrorLog As String

        Try
            ' Initialize objSearchResult
            objSearchResult = New clsSearchResultsInSpecT(mPeptideMods)

            ' Initialize variables
            strPreviousScan = String.Empty

            Try
                ' Open the input file and parse it
                ' Initialize the stream reader and the stream Text writer
                srDataFile = New System.IO.StreamReader(strInputFilePath)

                srResultFile = New System.IO.StreamWriter(strOutputFilePath)

                strErrorLog = String.Empty
                intResultsProcessed = 0
                Dim firstDataRec As Boolean = True
                ' Parse the input file
                Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing
                    strLineIn = srDataFile.ReadLine
                    If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then
                        blnValidSearchResult = ParseInSpecTResultsFileEntry(strLineIn, objSearchResult, strErrorLog, intResultsProcessed)
                        If Not blnValidSearchResult Then
                            'assume this is the first line, so write it as the header
                            WriteFileHeader(strLineIn, srResultFile, objSearchResult, strErrorLog)
                        End If

                        If blnValidSearchResult Then
                            If CSng(objSearchResult.FScore) >= FSCORE_THRESHOLD Or CSng(objSearchResult.TotalPRMScore) >= TOTALPRMSCORE_THRESHOLD Then
                                intResultID = intResultID + 1
                                'write the record
                                WriteFileRecord(intResultID, srResultFile, objSearchResult, strErrorLog)
                            End If
                        End If
                        ' Update the progress
                        UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))

                        intResultsProcessed += 1
                    End If
                Loop

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
                If Not srResultFile Is Nothing Then
                    srResultFile.Close()
                    srResultFile = Nothing
                End If
            End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Function determineBestRecord(ByRef rsltFile As System.IO.StreamWriter, ByRef intResultID As Integer, ByRef strErrorLog As String) As Boolean
        Dim tmpTblObjSearchResult As clsSearchResultsInSpecT
        tmpTblObjSearchResult = New clsSearchResultsInSpecT(mPeptideMods)

        Dim expression As String = "Scan > 0"
        ' Sort descending by column named CompanyName.
        Dim sortOrder As String = "Scan ASC, Charge ASC"
        Dim foundRows() As DataRow

        Dim strPrevCharge As String
        strPrevCharge = String.Empty

        Dim tmpScan As String = String.Empty
        ' Use the Select method to find all rows matching the filter.
        foundRows = m_tmpSearchTable.Select(expression, sortOrder)

        Dim i As Integer
        ' Print column 0 of each returned row.
        For i = 0 To foundRows.GetUpperBound(0)
            tmpScan = Str(foundRows(i)("Scan"))
            If i = 0 Then
                strPrevCharge = Str(foundRows(i)("Charge"))
                saveCurrentRecordIntoTempRecordTable(tmpTblObjSearchResult, foundRows, i, strErrorLog)
            ElseIf Str(foundRows(i)("Charge")) <> strPrevCharge Then
                intResultID = intResultID + 1
                'write the record
                WriteFileRecord(intResultID, rsltFile, tmpTblObjSearchResult, strErrorLog)
                saveCurrentRecordIntoTempRecordTable(tmpTblObjSearchResult, foundRows, i, strErrorLog)
            Else
                If CSng(foundRows(i)("TotalPRMScore")) > CSng(tmpTblObjSearchResult.TotalPRMScore) Then
                    saveCurrentRecordIntoTempRecordTable(tmpTblObjSearchResult, foundRows, i, strErrorLog)
                End If
            End If
            strPrevCharge = Str(foundRows(i)("Charge"))
        Next i
        intResultID = intResultID + 1
        WriteFileRecord(intResultID, rsltFile, tmpTblObjSearchResult, strErrorLog)

    End Function

    Private Function saveCurrentRecordIntoTempRecordTable(ByRef tmpSearchRecord As clsSearchResultsInSpecT, ByRef curSearchRecord() As DataRow, ByVal rowIdx As Integer, ByRef strErrorLog As String) As Boolean
        Dim tmpstring As String = String.Empty
        Try
            tmpSearchRecord.Scan = CStr(curSearchRecord(rowIdx)("Scan"))
            tmpSearchRecord.Annotation = CStr(curSearchRecord(rowIdx)("Annotation"))
            tmpSearchRecord.Protein = CStr(curSearchRecord(rowIdx)("Protein"))
            tmpSearchRecord.Charge = CStr(curSearchRecord(rowIdx)("Charge"))
            tmpSearchRecord.MQScore = CStr(curSearchRecord(rowIdx)("MQScore"))
            tmpSearchRecord.Length = CStr(curSearchRecord(rowIdx)("Length"))
            tmpSearchRecord.TotalPRMScore = CStr(curSearchRecord(rowIdx)("TotalPRMScore"))
            tmpSearchRecord.MedianPRMScore = CStr(curSearchRecord(rowIdx)("MedianPRMScore"))
            tmpSearchRecord.FractionY = CStr(curSearchRecord(rowIdx)("FractionY"))
            tmpSearchRecord.FractionB = CStr(curSearchRecord(rowIdx)("FractionB"))
            tmpSearchRecord.Intensity = CStr(curSearchRecord(rowIdx)("Intensity"))
            tmpSearchRecord.NTT = CStr(curSearchRecord(rowIdx)("NTT"))
            tmpSearchRecord.pValue = CStr(curSearchRecord(rowIdx)("pValue"))
            tmpSearchRecord.FScore = CStr(curSearchRecord(rowIdx)("FScore"))
            tmpSearchRecord.DeltaScore = CStr(curSearchRecord(rowIdx)("DeltaScore"))
            tmpSearchRecord.DeltaScoreOther = CStr(curSearchRecord(rowIdx)("DeltaScoreOther"))
            tmpSearchRecord.RecordNumber = CStr(curSearchRecord(rowIdx)("RecordNumber"))
            tmpSearchRecord.DBFilePos = CStr(curSearchRecord(rowIdx)("DBFilePos"))
            tmpSearchRecord.SpecFilePos = CStr(curSearchRecord(rowIdx)("SpecFilePos"))
        Catch ex As Exception
            strErrorLog &= "Error saving current record to temp record" & ControlChars.NewLine
            Return False
        End Try

        Return True

    End Function

    Private Function createTmpTable(ByRef strErrorLog As String) As Boolean

        Try
            m_tmpSearchTable.Columns.Add("Scan", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("Annotation", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("Protein", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("Charge", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("MQScore", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("Length", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("TotalPRMScore", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("MedianPRMScore", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("FractionY", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("FractionB", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("Intensity", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("NTT", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("pValue", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("FScore", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("DeltaScore", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("DeltaScoreOther", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("RecordNumber", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("DBFilePos", Type.GetType("System.String"))
            m_tmpSearchTable.Columns.Add("SpecFilePos", Type.GetType("System.String"))
        Catch ex As Exception
            strErrorLog &= "Error saving current record to temp record" & ControlChars.NewLine
            Return False
        End Try

        Return True

    End Function

    Private Function saveCurrentRecordIntoTempTable(ByRef curSearchRecord As clsSearchResultsInSpecT, ByRef strErrorLog As String) As Boolean
        Dim workRow As DataRow = m_tmpSearchTable.NewRow()
        Try
            workRow("Scan") = curSearchRecord.Scan
            workRow("Annotation") = curSearchRecord.Annotation
            workRow("Protein") = curSearchRecord.Protein
            workRow("Charge") = curSearchRecord.Charge
            workRow("MQScore") = curSearchRecord.MQScore
            workRow("Length") = curSearchRecord.Length
            workRow("TotalPRMScore") = curSearchRecord.TotalPRMScore
            workRow("MedianPRMScore") = curSearchRecord.MedianPRMScore
            workRow("FractionY") = curSearchRecord.FractionY
            workRow("FractionB") = curSearchRecord.FractionB
            workRow("Intensity") = curSearchRecord.Intensity
            workRow("NTT") = curSearchRecord.NTT
            workRow("pValue") = curSearchRecord.pValue
            workRow("FScore") = curSearchRecord.FScore
            workRow("DeltaScore") = curSearchRecord.DeltaScore
            workRow("DeltaScoreOther") = curSearchRecord.DeltaScoreOther
            workRow("RecordNumber") = curSearchRecord.RecordNumber
            workRow("DBFilePos") = curSearchRecord.DBFilePos
            workRow("SpecFilePos") = curSearchRecord.SpecFilePos
            m_tmpSearchTable.Rows.Add(workRow)
        Catch ex As Exception
            strErrorLog &= "Error saving current record to temp record" & ControlChars.NewLine
            Return False
        End Try

        Return True

    End Function

    Private Function WriteFileRecord(ByVal intResultID As Integer, ByRef srResultFile As System.IO.StreamWriter, ByRef rsltRecord As clsSearchResultsInSpecT, ByRef strErrorLog As String) As Boolean

        Try
            srResultFile.WriteLine(Str(intResultID) & ControlChars.Tab & _
                                   rsltRecord.Scan & ControlChars.Tab & _
                                   rsltRecord.Annotation & ControlChars.Tab & _
                                   rsltRecord.Protein & ControlChars.Tab & _
                                   rsltRecord.Charge & ControlChars.Tab & _
                                   rsltRecord.MQScore & ControlChars.Tab & _
                                   rsltRecord.Length & ControlChars.Tab & _
                                   rsltRecord.TotalPRMScore & ControlChars.Tab & _
                                   rsltRecord.MedianPRMScore & ControlChars.Tab & _
                                   rsltRecord.FractionY & ControlChars.Tab & _
                                   rsltRecord.FractionB & ControlChars.Tab & _
                                   rsltRecord.Intensity & ControlChars.Tab & _
                                   rsltRecord.NTT & ControlChars.Tab & _
                                   rsltRecord.pValue & ControlChars.Tab & _
                                   rsltRecord.FScore & ControlChars.Tab & _
                                   rsltRecord.DeltaScore & ControlChars.Tab & _
                                   rsltRecord.DeltaScoreOther & ControlChars.Tab & _
                                   rsltRecord.RecordNumber & ControlChars.Tab & _
                                   rsltRecord.DBFilePos & ControlChars.Tab & _
                                   rsltRecord.SpecFilePos & ControlChars.Tab)

        Catch ex As Exception
            strErrorLog &= "Error writing first hits record" & ControlChars.NewLine
            Return False
        End Try

        Return True

    End Function

    Private Function WriteFileHeader(ByRef strLineIn As String, ByRef srResultFile As System.IO.StreamWriter, ByRef clsSearchResultsInSpecT As clsSearchResultsInSpecT, ByRef strErrorLog As String) As Boolean
        Dim strSplitLine As String()
        Dim tmpStr As String = ""

        'Replace some header names with new name
        Try
            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)
            tmpStr = strSplitLine(eInspectFileColumns.Scan)
            srResultFile.WriteLine("ResultID" & ControlChars.Tab & _
                                   "Scan" & ControlChars.Tab & _
                                   "Peptide" & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.Protein) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.Charge) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.MQScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.Length) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.TotalPRMScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.MedianPRMScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.FractionY) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.FractionB) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.Intensity) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.NTT) & ControlChars.Tab & _
                                   "PValue" & ControlChars.Tab & _
                                   "FScore" & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.DeltaScore) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.DeltaScoreOther) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.RecordNumber) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.DBFilePos) & ControlChars.Tab & _
                                   strSplitLine(eInspectFileColumns.SpecFilePos))

        Catch ex As Exception
            strErrorLog &= "Error writing first hits header" & ControlChars.NewLine
            Return False
        End Try

        Return True

    End Function

    Private Function ParseInSpecTResultsFileEntry(ByRef strLineIn As String, ByRef objSearchResult As clsSearchResultsInSpecT, ByRef strErrorLog As String, ByVal intResultsProcessed As Integer) As Boolean

        Dim strSplitLine() As String

        Dim blnValidSearchResult As Boolean
        blnValidSearchResult = False

        ' The following are the headers in the inpsect file
        'SpectrumFile
        'Scan
        'Annotation (Peptide)
        'Protein
        'Charge
        'MQScore
        'Length
        'TotalPRMScore
        'MedianPRMScore
        'FractionY
        'FractionB
        'Intensity
        'NTT
        'p-value
        'F-Score
        'DeltaScore
        'DeltaScoreOther
        'RecordNumber
        'DBFilePos
        'SpecFilePos

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
                .SpectrumFile = strSplitLine(eInspectFileColumns.SpectrumFile)
                .Scan = strSplitLine(eInspectFileColumns.Scan)
                .Annotation = ReplaceTerminus(strSplitLine(eInspectFileColumns.Annotation))
                .Protein = TruncateProteinName(strSplitLine(eInspectFileColumns.Protein))
                .Charge = strSplitLine(eInspectFileColumns.Charge)
                .MQScore = strSplitLine(eInspectFileColumns.MQScore)
                .Length = strSplitLine(eInspectFileColumns.Length)
                .TotalPRMScore = strSplitLine(eInspectFileColumns.TotalPRMScore)
                .MedianPRMScore = strSplitLine(eInspectFileColumns.MedianPRMScore)
                .FractionY = strSplitLine(eInspectFileColumns.FractionY)
                .FractionB = strSplitLine(eInspectFileColumns.FractionB)
                .Intensity = strSplitLine(eInspectFileColumns.Intensity)
                .NTT = strSplitLine(eInspectFileColumns.NTT)
                .pValue = strSplitLine(eInspectFileColumns.pvalue)
                .FScore = strSplitLine(eInspectFileColumns.FScore)
                .DeltaScore = strSplitLine(eInspectFileColumns.DeltaScore)
                .DeltaScoreOther = strSplitLine(eInspectFileColumns.DeltaScoreOther)
                .RecordNumber = strSplitLine(eInspectFileColumns.RecordNumber)
                .DBFilePos = strSplitLine(eInspectFileColumns.DBFilePos)
                .SpecFilePos = strSplitLine(eInspectFileColumns.SpecFilePos)
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

                MyBase.ResetProgress("Parsing " & System.IO.Path.GetFileName(strInputFilePath))

                If CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
                    Try
                        ' Obtain the full path to the input file
                        ioFile = New System.IO.FileInfo(strInputFilePath)
                        strInputFilePathFull = ioFile.FullName

                        ' Define the first hits output file name based on strInputFilePath
                        strOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                        strOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strOutputFilePath & SEQUEST_FIRST_HITS_FILE_SUFFIX)

                        blnSuccess = CreateFHTResultsFile(strInputFilePath, strOutputFilePath)

                        'Define the synopsis output file name based on strInputFilePath
                        strSynOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                        strSynOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strSynOutputFilePath & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                        blnSuccess = CreateSYNResultsFile(strInputFilePath, strSynOutputFilePath)

                        '                        blnSuccess = ParseInSpecTResultsFile(strInputFilePathFull, strOutputFilePath, False)

                    Catch ex As Exception
                        SetErrorMessage("Error calling ParseInSpecTResultsFile" & ex.Message)
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

        If inpString.Trim().StartsWith(TERMINUS_SYMBOL_INSPECT_B) Then
            inpString = "-." + inpString.Substring((Len(inpString) - (Len(inpString) - 2)))
        End If

        If inpString.Trim().StartsWith(TERMINUS_SYMBOL_INSPECT_E) Then
            inpString = ".-" + inpString.Substring((Len(inpString) - (Len(inpString) - 2)))
        End If

        If inpString.Trim().EndsWith(TERMINUS_SYMBOL_INSPECT_B) Then
            inpString = inpString.Substring(0, (Len(inpString) - 2)) & "-."
        End If

        If inpString.Trim().EndsWith(TERMINUS_SYMBOL_INSPECT_E) Then
            inpString = inpString.Substring(0, (Len(inpString) - 2)) & ".-"
        End If

        Return inpString

    End Function

    Private Function TruncateProteinName(ByVal strProteinNameAndDescription As String) As String

        Dim intIndex As Integer

        intIndex = strProteinNameAndDescription.IndexOf(" "c)
        If intIndex > 0 Then
            strProteinNameAndDescription = strProteinNameAndDescription.Substring(0, intIndex)
        End If

        Return strProteinNameAndDescription

    End Function

End Class
