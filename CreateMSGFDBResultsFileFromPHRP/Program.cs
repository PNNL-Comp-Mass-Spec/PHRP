Option Strict On

Imports System.IO
Imports System.Reflection
Imports System.Text
Imports System.Text.RegularExpressions
Imports System.Threading
Imports PHRPReader
Imports PRISM

' This program reads a PHRP-compatible _msgfdb_fht.txt file and creates the
' equivalent tab-delimited _msgfplus.tsv file that would have been created
' by MSGFPlus when converting the .mzIdentML file to a .tsv file
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started March 14, 2013
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/ or http://panomics.pnnl.gov/
' -------------------------------------------------------------------------------

Module modMain
    Public Const PROGRAM_DATE As String = "March 14, 2013"
    Private Const MASS_C13 As Double = 1.00335483

    Private WithEvents mPHRPReader As clsPHRPReader

    Private mInputFilePath As String
    Private mOutputFilePath As String

    Public Function Main() As Integer

        Dim intReturnCode As Integer
        Dim commandLineParser As New clsParseCommandLine
        Dim blnProceed As Boolean

        intReturnCode = 0
        mInputFilePath = String.Empty
        mOutputFilePath = String.Empty

        Try
            blnProceed = False
            If commandLineParser.ParseCommandLine Then
                If SetOptionsUsingCommandLineParameters(commandLineParser) Then blnProceed = True
            End If

            If Not blnProceed OrElse
               commandLineParser.NeedToShowHelp OrElse
               commandLineParser.ParameterCount + commandLineParser.NonSwitchParameterCount = 0 OrElse
               mInputFilePath.Length = 0 Then
                ShowProgramHelp()
                intReturnCode = -1
            Else

                ConvertFile()

            End If

        Catch ex As Exception
            ShowErrorMessage("Error occurred in modMain->Main", ex)
            intReturnCode = -1
        End Try

        Return intReturnCode

    End Function

    Private Function CleanupPeptide(strPeptide As String) As String

        ' ReSharper disable once UseImplicitlyTypedVariableEvident
        Static reFindItraq As Regex = New Regex("^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", RegexOptions.Compiled Or RegexOptions.IgnoreCase)

        Dim strPrimarySequence As String = String.Empty
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Dim reMatch As Match

        If clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, strPrimarySequence, strPrefix, strSuffix) Then
            ' Look for an N-terminal iTraq mod
            reMatch = reFindItraq.Match(strPrimarySequence)

            If reMatch.Success Then
                strPeptide = strPrefix & "." & reMatch.Groups(2).Value & reMatch.Groups(1).Value & reMatch.Groups(3).Value & "." & strSuffix
            End If
        End If

        Return strPeptide
    End Function

    Private Function ConvertFile() As Boolean

        Try
            Dim fiInputFile = New FileInfo(mInputFilePath)

            If Not fiInputFile.Exists Then
                ShowErrorMessage("Input file not found: " + fiInputFile.FullName)
                Return False
            End If

            If String.IsNullOrEmpty(mOutputFilePath) Then
                ' Auto-define the output file

                mOutputFilePath = Path.GetFileNameWithoutExtension(fiInputFile.Name)
                If mOutputFilePath.ToLower.EndsWith("_msgfplus_fht") OrElse mOutputFilePath.EndsWith("_msgfplus_syn") Then
                    mOutputFilePath = mOutputFilePath.Substring(0, mOutputFilePath.Length - 11)
                End If
                mOutputFilePath = Path.Combine(fiInputFile.Directory.FullName, mOutputFilePath & "_msgfplus.tsv")
            End If

            mPHRPReader = New clsPHRPReader(fiInputFile.FullName, clsPHRPReader.ePeptideHitResultType.Unknown, True, False, False)
            mPHRPReader.EchoMessagesToConsole = False
            mPHRPReader.SkipDuplicatePSMs = False

            If Not mPHRPReader.CanRead Then
                ShowErrorMessage("Aborting since PHRPReader is not ready: " + mPHRPReader.ErrorMessage)
                Return False
            End If

            Dim oMassCalculator = New clsPeptideMassCalculator()

            Using swOutFile = New StreamWriter(New FileStream(mOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))

                Dim lstValues = New List(Of String)
                Dim intPSMsRead = 0

                ' Write the header line
                swOutFile.WriteLine(FlattenList(New List(Of String) From {"#SpecFile", "SpecID", "ScanNum", "FragMethod", "Precursor", "IsotopeError", "PrecursorError(ppm)", "Charge", "Peptide", "Protein", "DeNovoScore", "MSGFScore", "SpecEValue", "EValue", "QValue", "PepQValue"}))

                Dim strMassErrorPPM As String
                Dim intIsotopeErrorComputed As Integer
                Dim strIsotopeError As String

                Do While mPHRPReader.MoveNext()
                    Dim oPsm As clsPSM = mPHRPReader.CurrentPSM
                    intPSMsRead += 1
                    lstValues.Clear()

                    intIsotopeErrorComputed = 0
                    strMassErrorPPM = GetCorrectedMassErrorPPM(oPsm, intIsotopeErrorComputed)

                    lstValues.Add(mPHRPReader.DatasetName & "_dta.txt")                                             ' #SpecFile
                    lstValues.Add("index=" & intPSMsRead)                                                           ' SpecID
                    lstValues.Add(oPsm.ScanNumber.ToString())                                                       ' ScanNum
                    lstValues.Add(oPsm.CollisionMode)                                                               ' FragMethod
                    lstValues.Add(GetPrecursorMZ(oMassCalculator, oPsm))                                                             ' Precursor m/z

                    strIsotopeError = GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_Isotope_Error, "0")
                    If strIsotopeError = "0" And intIsotopeErrorComputed <> 0 Then
                        strIsotopeError = intIsotopeErrorComputed.ToString()
                    End If

                    lstValues.Add(strIsotopeError)                                                                  ' IsotopeError
                    lstValues.Add(strMassErrorPPM)                                                                  ' PrecursorError(ppm)
                    lstValues.Add(oPsm.Charge.ToString())                                                           ' Charge
                    lstValues.Add(CleanupPeptide(oPsm.PeptideWithNumericMods))                                      ' Peptide
                    lstValues.Add(oPsm.ProteinFirst)                                                                ' Protein
                    lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore, "0"))                 ' DeNovoScore
                    lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, "0"))                   ' MSGFScore
                    lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue, "0"))           ' SpecEValue
                    lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_EValue, "0"))                      ' EValue
                    lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_QValue, "0"))                      ' QValue
                    lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_PepQValue, "0"))                   ' PepQValue


                    swOutFile.WriteLine(FlattenList(lstValues))
                Loop

            End Using

            Console.WriteLine("Created file " & mOutputFilePath)

        Catch ex As Exception
            ShowErrorMessage("Error occurred in modMain->ConvertFile", ex)
            Return False
        End Try

        Return True
    End Function

    Private Function FlattenList(lstValues As List(Of String)) As String
        Return FlattenList(lstValues, ControlChars.Tab)
    End Function

    Private Function FlattenList(lstValues As List(Of String), chSepChar As Char) As String
        Dim sbOutline = New StringBuilder()

        For intIndex = 0 To lstValues.Count - 1
            If intIndex > 0 Then
                sbOutline.Append(chSepChar)
            End If
            sbOutline.Append(lstValues(intIndex))
        Next

        Return sbOutline.ToString()
    End Function

    Private Function GetAppVersion() As String
        'Return System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")"

        Return Assembly.GetExecutingAssembly().GetName().Version.ToString() & " (" & PROGRAM_DATE & ")"
    End Function

    Private Function GetCorrectedMassErrorPPM(oPsm As clsPSM, ByRef intIsotopeError As Integer) As String

        Dim dblDelM As Double
        Dim dblMassErrorPPM As Double = 0
        intIsotopeError = 0

        If Double.TryParse(oPsm.MassErrorDa, dblDelM) Then

            ' Examine dblDelM to determine which isotope was chosen
            If dblDelM >= -0.5 Then
                ' This is the typical case
                Do While dblDelM > 0.5
                    dblDelM -= MASS_C13
                    intIsotopeError += 1
                Loop
            Else
                ' This happens less often; but we'll still account for it
                ' In this case, intCorrectionCount will be negative
                Do While dblDelM < -0.5
                    dblDelM += MASS_C13
                    intIsotopeError -= 1
                Loop

            End If

            dblMassErrorPPM = clsPeptideMassCalculator.MassToPPM(dblDelM, oPsm.PrecursorNeutralMass)
        End If

        Return dblMassErrorPPM.ToString("0.0000")
    End Function

    Private Function GetPrecursorMZ(oMassCalculator As clsPeptideMassCalculator, oPsm As clsPSM) As String
        Return oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge).ToString()
    End Function

    Private Function GetScore(oPsm As clsPSM, strScoreName As String, strValueIfMissing As String) As String
        Dim strScoreValue As String = String.Empty

        If Not oPsm.TryGetScore(strScoreName, strScoreValue) Then
            strScoreValue = strValueIfMissing
        End If

        Return strScoreValue

    End Function

    Private Function SetOptionsUsingCommandLineParameters(commandLineParser As clsParseCommandLine) As Boolean
        ' Returns True if no problems; otherwise, returns false

        Dim strValue As String = String.Empty
        Dim lstValidParameters = New List(Of String) From {"I", "O"}

        Try
            ' Make sure no invalid parameters are present
            If commandLineParser.InvalidParametersPresent(lstValidParameters) Then
                ShowErrorMessage("Invalid command line parameters",
                  (From item In commandLineParser.InvalidParameters(lstValidParameters) Select "/" + item).ToList())
                Return False
            Else
                With commandLineParser
                    ' Query commandLineParser to see if various parameters are present
                    If .RetrieveValueForParameter("I", strValue) Then
                        mInputFilePath = String.Copy(strValue)
                    ElseIf .NonSwitchParameterCount > 0 Then
                        mInputFilePath = .RetrieveNonSwitchParameter(0)
                    End If

                    If .RetrieveValueForParameter("O", strValue) Then mOutputFilePath = String.Copy(strValue)

                End With

                Return True
            End If

        Catch ex As Exception
            ShowErrorMessage("Error parsing the command line parameters", ex)
        End Try

        Return False

    End Function

    Private Sub ShowErrorMessage(message As String, Optional ex As Exception = Nothing)
        ConsoleMsgUtils.ShowError(message, ex)
    End Sub

    Private Sub ShowErrorMessage(title As String, errorMessages As List(Of String))
        ConsoleMsgUtils.ShowErrors(title, errorMessages)
    End Sub

    Private Sub ShowProgramHelp()

        Try

            Console.WriteLine("This program reads a PHRP-compatible _msgfplus_fht.txt file and creates the equivalent tab-delimited _msgfplus.tsv file that would have been created by MSGFPlus when converting the .mzIdentML file to a .tsv file")
            Console.WriteLine()
            Console.WriteLine("Program syntax:" & ControlChars.NewLine & Path.GetFileName(Assembly.GetExecutingAssembly().Location) & _
             " InputFilePath [/O:OutputFilePath]")
            Console.WriteLine()

            Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006")
            Console.WriteLine("Version: " & GetAppVersion())

            Console.WriteLine()

            Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com")
            Console.WriteLine("Website: http://omics.pnl.gov/ or http://panomics.pnnl.gov/")

            ' Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
            Thread.Sleep(750)

        Catch ex As Exception
            ShowErrorMessage("Error displaying the program syntax", ex)
        End Try

    End Sub

    Private Sub PHRPReader_ErrorEvent(message As String, ex As Exception) Handles mPHRPReader.ErrorEvent
        ShowErrorMessage(message, ex)
    End Sub

    Private Sub PHRPReader_StatusEvent(message As String) Handles mPHRPReader.StatusEvent
        Console.WriteLine(message)
    End Sub

    Private Sub PHRPReader_WarningEvent(strWarningMessage As String) Handles mPHRPReader.WarningEvent
        Console.WriteLine("Warning: " & strWarningMessage)
    End Sub
End Module
