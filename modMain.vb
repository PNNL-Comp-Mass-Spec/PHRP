Option Strict On

' This program processes XTandem search results or Sequest search results to
' determine the modifications present, determine the cleaveage and terminus state
' of each peptide, and compute the monoisotopic mass of each peptide.  See 
' clsSequestSynopsisFileProcessor and clsXTandemResultsConverter for 
' additional information
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

Module modMain
    Public Const PROGRAM_DATE As String = "October 6, 2006"

    Private mInputFilePath As String
    Private mOutputFolderName As String                         ' Optional
    Private mParameterFilePath As String                        ' Optional

    Private mMassCorrectionTagsFilePath As String               ' Optional
    Private mModificationDefinitionsFilePath As String          ' Optional

    Private mOutputFolderAlternatePath As String                ' Optional
    Private mRecreateFolderHierarchyInAlternatePath As Boolean  ' Optional

    Private mRecurseFolders As Boolean
    Private mRecurseFoldersMaxLevels As Integer

    Private mQuietMode As Boolean

    Private WithEvents mPeptideHitResultsProcRunner As clsPeptideHitResultsProcRunner
    Private mLastProgressReportTime As System.DateTime
    Private mLastProgressReportValue As Integer

    Private Sub DisplayProgressPercent(ByVal intPercentComplete As Integer, ByVal blnAddCarriageReturn As Boolean)
        If blnAddCarriageReturn Then
            Console.WriteLine()
        End If
        If intPercentComplete > 100 Then intPercentComplete = 100
        Console.Write("Processing: " & intPercentComplete.ToString & "% ")
        If blnAddCarriageReturn Then
            Console.WriteLine()
        End If
    End Sub

    Public Function Main() As Integer
        ' Returns 0 if no error, error code if an error

        Dim intReturnCode As Integer
        Dim objParseCommandLine As New SharedVBNetRoutines.clsParseCommandLine
        Dim blnProceed As Boolean

        ''Dim objTest As New clsPeptideCleavageStateCalculator
        ''Dim strSeq, strCleanSeq, strPrefix, strSuffix As String
        ''strSeq = "A.BCDE.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "-.BCDE.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "A.BCDE.-" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "A.B.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "A.BCDE" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "BCDE.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "FA.BCDE.FG" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "BCDE." : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = ".BCDE" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = ".F." : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "F..E" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "AF..EF" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
        ''strSeq = "AFF..EF" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)

        ''Return 0

        intReturnCode = 0
        mInputFilePath = String.Empty
        mOutputFolderName = String.Empty
        mParameterFilePath = String.Empty

        mMassCorrectionTagsFilePath = String.Empty
        mModificationDefinitionsFilePath = String.Empty

        mRecurseFolders = False
        mRecurseFoldersMaxLevels = 0

        Try
            blnProceed = False
            If objParseCommandLine.ParseCommandLine Then
                If SetOptionsUsingCommandLineParameters(objParseCommandLine) Then blnProceed = True
            End If

            If Not blnProceed OrElse objParseCommandLine.NeedToShowHelp OrElse objParseCommandLine.ParameterCount = 0 OrElse mInputFilePath.Length = 0 Then
                ShowProgramHelp()
                intReturnCode = -1
            Else
                mPeptideHitResultsProcRunner = New clsPeptideHitResultsProcRunner
                mPeptideHitResultsProcRunner.ShowMessages = Not mQuietMode

                With mPeptideHitResultsProcRunner
                    ' Note: These options will get overridden if defined in the parameter file
                    .MassCorrectionTagsFilePath = mMassCorrectionTagsFilePath
                    .ModificationDefinitionsFilePath = mModificationDefinitionsFilePath
                    .WarnMissingParameterFileSection = True
                End With

                If mRecurseFolders Then
                    If mPeptideHitResultsProcRunner.ProcessFilesAndRecurseFolders(mInputFilePath, mOutputFolderName, mOutputFolderAlternatePath, mRecreateFolderHierarchyInAlternatePath, mParameterFilePath, mRecurseFoldersMaxLevels) Then
                        intReturnCode = 0
                    Else
                        intReturnCode = mPeptideHitResultsProcRunner.ErrorCode
                    End If
                Else
                    If mPeptideHitResultsProcRunner.ProcessFilesWildcard(mInputFilePath, mOutputFolderName, mParameterFilePath) Then
                        intReturnCode = 0
                    Else
                        intReturnCode = mPeptideHitResultsProcRunner.ErrorCode
                        If intReturnCode <> 0 AndAlso Not mQuietMode Then
                            MsgBox("Error while processing: " & mPeptideHitResultsProcRunner.GetErrorMessage(), MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
                        End If
                    End If
                End If

                DisplayProgressPercent(mLastProgressReportValue, True)
            End If

        Catch ex As Exception
            If mQuietMode Then
                Throw ex
            Else
                MsgBox("Error occurred in modMain->Main: " & ControlChars.NewLine & ex.Message, MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
            End If
            intReturnCode = -1
        End Try

        Return intReturnCode

    End Function

    Private Function SetOptionsUsingCommandLineParameters(ByVal objParseCommandLine As SharedVBNetRoutines.clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

        Dim strValue As String
        Dim strValidParameters() As String = New String() {"I", "O", "P", "M", "T", "S", "A", "R", "Q"}

        Try
            ' Make sure no invalid parameters are present
            If objParseCommandLine.InvalidParametersPresent(strValidParameters) Then
                Return False
            Else
                With objParseCommandLine
                    ' Query objParseCommandLine to see if various parameters are present
                    If .RetrieveValueForParameter("I", strValue) Then mInputFilePath = strValue
                    If .RetrieveValueForParameter("O", strValue) Then mOutputFolderName = strValue
                    If .RetrieveValueForParameter("P", strValue) Then mParameterFilePath = strValue
                    If .RetrieveValueForParameter("M", strValue) Then mModificationDefinitionsFilePath = strValue
                    If .RetrieveValueForParameter("T", strValue) Then mMassCorrectionTagsFilePath = strValue

                    If .RetrieveValueForParameter("S", strValue) Then
                        mRecurseFolders = True
                        If IsNumeric(strValue) Then
                            mRecurseFoldersMaxLevels = CInt(strValue)
                        End If
                    End If
                    If .RetrieveValueForParameter("A", strValue) Then mOutputFolderAlternatePath = strValue
                    If .RetrieveValueForParameter("R", strValue) Then mRecreateFolderHierarchyInAlternatePath = True

                    If .RetrieveValueForParameter("Q", strValue) Then mQuietMode = True
                End With

                Return True
            End If

        Catch ex As Exception
            If mQuietMode Then
                Throw New System.Exception("Error parsing the command line parameters", ex)
            Else
                MsgBox("Error parsing the command line parameters: " & ControlChars.NewLine & ex.Message, MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
            End If
        End Try

    End Function

    Private Sub ShowProgramHelp()

        Dim strSyntax As String
        Dim ioPath As System.IO.Path

        Try

            strSyntax = "This program reads in an XTandem results file (XML format) or Sequest Synopsis/First Hits file and creates a tab-delimited text file with the data.  "
            strSyntax &= "It will insert modification symbols into the peptide sequences for modified peptides.  Parallel files will be created containing sequence info and modification details.  "
            strSyntax &= "The user can optionally provide a modification definition file which specifies the symbol to use for each modification mass." & ControlChars.NewLine & ControlChars.NewLine
            strSyntax &= "Program syntax:" & ControlChars.NewLine & ioPath.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location)
            strSyntax &= " /I:InputFilePath_xt.xml [/O:OutputFolderPath] [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath] [/T:MassCorrectionTagsFilePath] [/S:[MaxLevel]] [/A:AlternateOutputFolderPath] [/R] [/W] [/D] [/Q]" & ControlChars.NewLine & ControlChars.NewLine
            strSyntax &= "The input file should be an XTandem Results file (_xt.xml), a Sequest Synopsis File (_syn.txt), or a Sequest First Hits file (_fht.txt)." & ControlChars.NewLine
            strSyntax &= "The output folder switch is optional.  If omitted, the output file will be created in the same folder as the input file." & ControlChars.NewLine
            strSyntax &= "The parameter file path is optional.  If included, it should point to a valid XML parameter file." & ControlChars.NewLine
            strSyntax &= "Use /M to specify the file containing the modification definitions.  This file should be tab delimited, with the first column containing the modification symbol, the second column containing the modification mass, plus optionally a third column listing the residues that can be modified with the given mass (1 letter residue symbols, no need to separated with commas or spaces)." & ControlChars.NewLine
            strSyntax &= "Use /T to specify the file containing the mass correction tag info.  This file should be tab delimited, with the first column containing the mass correction tag name and the second column containing the mass (the name cannot contain commas or colons and can be, at most, 8 characters long)." & ControlChars.NewLine
            strSyntax &= "Use /S to process all valid files in the input folder and subfolders. Include a number after /S (like /S:2) to limit the level of subfolders to examine." & ControlChars.NewLine
            strSyntax &= "When using /S, you can redirect the output of the results using /A." & ControlChars.NewLine
            strSyntax &= "When using /S, you can use /R to re-create the input folder hierarchy in the alternate output folder (if defined)." & ControlChars.NewLine
            strSyntax &= "The optional /Q switch will suppress all error messages." & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "This is version " & System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com" & ControlChars.NewLine
            strSyntax &= "Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Licensed under the Apache License, Version 2.0; you may not use this file except in compliance with the License.  "
            strSyntax &= "You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0" & ControlChars.NewLine & ControlChars.NewLine

            strSyntax &= "Notice: This computer software was prepared by Battelle Memorial Institute, "
            strSyntax &= "hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the "
            strSyntax &= "Department of Energy (DOE).  All rights in the computer software are reserved "
            strSyntax &= "by DOE on behalf of the United States Government and the Contractor as "
            strSyntax &= "provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY "
            strSyntax &= "WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS "
            strSyntax &= "SOFTWARE.  This notice including this sentence must appear on any copies of "
            strSyntax &= "this computer software." & ControlChars.NewLine

            If Not mQuietMode Then
                MsgBox(strSyntax, MsgBoxStyle.Information Or MsgBoxStyle.OKOnly, "Syntax")
            End If

        Catch ex As Exception
            If mQuietMode Then
                Throw New System.Exception("Error displaying the program syntax", ex)
            Else
                MsgBox("Error displaying the program syntax: " & ControlChars.NewLine & ex.Message, MsgBoxStyle.Exclamation Or MsgBoxStyle.OKOnly, "Error")
            End If
        End Try

    End Sub

    Private Sub mXTandemResultsConverter_ProgressChanged(ByVal taskDescription As String, ByVal percentComplete As Single) Handles mPeptideHitResultsProcRunner.ProgressChanged
        Const PERCENT_REPORT_INTERVAL As Integer = 25
        Const PROGRESS_DOT_INTERVAL_MSEC As Integer = 250

        If percentComplete >= mLastProgressReportValue Then
            If mLastProgressReportValue > 0 Then
                Console.WriteLine()
            End If
            DisplayProgressPercent(mLastProgressReportValue, False)
            mLastProgressReportValue += PERCENT_REPORT_INTERVAL
            mLastProgressReportTime = DateTime.Now
        Else
            If DateTime.Now.Subtract(mLastProgressReportTime).TotalMilliseconds > PROGRESS_DOT_INTERVAL_MSEC Then
                mLastProgressReportTime = DateTime.Now
                Console.Write(".")
            End If
        End If
    End Sub

    Private Sub mXTandemResultsConverter_ProgressReset() Handles mPeptideHitResultsProcRunner.ProgressReset
        mLastProgressReportTime = DateTime.Now
        mLastProgressReportValue = 0
    End Sub
End Module
