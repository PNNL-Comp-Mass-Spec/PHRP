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
    Public Const PROGRAM_DATE As String = "October 4, 2008"

    Private mInputFilePath As String
    Private mOutputFolderName As String                         ' Optional
    Private mParameterFilePath As String                        ' Optional

    Private mMassCorrectionTagsFilePath As String               ' Optional
    Private mModificationDefinitionsFilePath As String          ' Optional
    Private mSearchToolParameterFilePath As String              ' Optional

    Private mOutputFolderAlternatePath As String                ' Optional
    Private mRecreateFolderHierarchyInAlternatePath As Boolean  ' Optional

    Private mRecurseFolders As Boolean
    Private mRecurseFoldersMaxLevels As Integer

    Private mLogMessagesToFile As Boolean
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
        Dim objParseCommandLine As New clsParseCommandLine
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
        mSearchToolParameterFilePath = String.Empty

        mRecurseFolders = False
        mRecurseFoldersMaxLevels = 0

        mQuietMode = False
        mLogMessagesToFile = False

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

                With mPeptideHitResultsProcRunner
                    .ShowMessages = Not mQuietMode
                    .LogMessagesToFile = mLogMessagesToFile

                    ' Note: These options will get overridden if defined in the parameter file
                    .MassCorrectionTagsFilePath = mMassCorrectionTagsFilePath
                    .ModificationDefinitionsFilePath = mModificationDefinitionsFilePath
                    .SearchToolParameterFilePath = mSearchToolParameterFilePath

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
                            Console.WriteLine("Error while processing: " & mPeptideHitResultsProcRunner.GetErrorMessage())
                        End If
                    End If
                End If

                DisplayProgressPercent(mLastProgressReportValue, True)
            End If

        Catch ex As Exception
            If mQuietMode Then
                Throw ex
            Else
                Console.WriteLine("Error occurred in modMain->Main: " & ControlChars.NewLine & ex.Message)
            End If
            intReturnCode = -1
        End Try

        Return intReturnCode

    End Function

    Private Function SetOptionsUsingCommandLineParameters(ByVal objParseCommandLine As clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

        Dim strValue As String = String.Empty
        Dim strValidParameters() As String = New String() {"I", "O", "P", "M", "T", "N", "S", "A", "R", "Q"}

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
                    If .RetrieveValueForParameter("N", strValue) Then mSearchToolParameterFilePath = strValue

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
                Console.WriteLine("Error parsing the command line parameters: " & ControlChars.NewLine & ex.Message)
            End If
        End Try

    End Function

    Private Sub ShowProgramHelp()

        Dim strSyntax As String

        Try

            Console.WriteLine("This program reads in an XTandem results file (XML format), Sequest Synopsis/First Hits file, or Inspect search result file, and creates a tab-delimited text file with the data.  ")
            Console.WriteLine("It will insert modification symbols into the peptide sequences for modified peptides.  Parallel files will be created containing sequence info and modification details.  ")
            Console.WriteLine("The user can optionally provide a modification definition file which specifies the symbol to use for each modification mass.")
            Console.WriteLine()
            Console.WriteLine("Program syntax:" & ControlChars.NewLine & System.IO.Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location) & _
                                        " /I:InputFilePath_xt.xml [/O:OutputFolderPath]")
            Console.WriteLine(" [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath]")
            Console.WriteLine(" [/T:MassCorrectionTagsFilePath] [/N:SearchToolParameterFilePath]")
            Console.WriteLine(" [/S:[MaxLevel]] [/A:AlternateOutputFolderPath] [/R] [/W] [/D] [/Q]")
            Console.WriteLine()
            Console.WriteLine("The input file should be an XTandem Results file (_xt.xml), a Sequest Synopsis File (_syn.txt), a Sequest First Hits file (_fht.txt), or an Inspect results file (_inspect.txt).")
            Console.WriteLine("The output folder switch is optional.  If omitted, the output file will be created in the same folder as the input file.")
            Console.WriteLine("The parameter file path is optional.  If included, it should point to a valid XML parameter file.")
            Console.WriteLine()
            Console.WriteLine("Use /M to specify the file containing the modification definitions.  This file should be tab delimited, with the first column containing the modification symbol, the second column containing the modification mass, plus optionally a third column listing the residues that can be modified with the given mass (1 letter residue symbols, no need to separated with commas or spaces).")
            Console.WriteLine("Use /T to specify the file containing the mass correction tag info.  This file should be tab delimited, with the first column containing the mass correction tag name and the second column containing the mass (the name cannot contain commas or colons and can be, at most, 8 characters long).")
            Console.WriteLine("Use /N to specify the parameter file provided to the search tool.  This is only used when processing Inspect files.")
            Console.WriteLine()
            Console.WriteLine("Use /S to process all valid files in the input folder and subfolders. Include a number after /S (like /S:2) to limit the level of subfolders to examine.")
            Console.WriteLine("When using /S, you can redirect the output of the results using /A.")
            Console.WriteLine("When using /S, you can use /R to re-create the input folder hierarchy in the alternate output folder (if defined).")
            Console.WriteLine("The optional /Q switch will suppress all error messages.")
            Console.WriteLine()

            Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006")
            Console.WriteLine()

            Console.WriteLine("This is version " & System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")")
            Console.WriteLine()

            Console.WriteLine("E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com")
            Console.WriteLine("Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/")
            Console.WriteLine()

            Console.WriteLine("Licensed under the Apache License, Version 2.0; you may not use this file except in compliance with the License.  " & _
                              "You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0")
            Console.WriteLine()

            Console.WriteLine("Notice: This computer software was prepared by Battelle Memorial Institute, " & _
                              "hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the " & _
                              "Department of Energy (DOE).  All rights in the computer software are reserved " & _
                              "by DOE on behalf of the United States Government and the Contractor as " & _
                              "provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY " & _
                              "WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS " & _
                              "SOFTWARE.  This notice including this sentence must appear on any copies of " & _
                              "this computer software.")

        Catch ex As Exception
            Console.WriteLine("Error displaying the program syntax: " & ex.Message)
        End Try

    End Sub

    Private Sub mPeptideHitResultsProcRunner_ProgressChanged(ByVal taskDescription As String, ByVal percentComplete As Single) Handles mPeptideHitResultsProcRunner.ProgressChanged
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

    Private Sub mPeptideHitResultsProcRunner_ProgressReset() Handles mPeptideHitResultsProcRunner.ProgressReset
        mLastProgressReportTime = DateTime.Now
        mLastProgressReportValue = 0
    End Sub
End Module
