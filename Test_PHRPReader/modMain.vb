Imports System.IO
Imports System.Text
Imports System.Text.RegularExpressions
Imports PHRPReader

Module modMain

    Private WithEvents mPHRPReader As clsPHRPReader

    Public Sub Main()

        'Const strSequestSynFilePath = "Seq201304121552_Auto934225\Firestone_Soil_07_18_05APR13_Frodo_12-12-04_syn.txt"
        'Const strMSGFPlusSynFilePath = "MSG201304261714_Auto938181\FSFA-299b_25Apr13_Methow_13-02-13_msgfdb_syn.txt"

        Const strSequestFolder = "\\proto-7\VOrbi05\2013_2\Firestone_Soil_07_18_05APR13_Frodo_12-12-04\Seq201304121552_Auto934225"
        ' 		Const strMSGFPlusFolder = "MSG201304261714_Auto938181"

        ' Const strMSGFPlusFolder = "\\proto-7\VOrbiETD03\2015_1\proteogeomics_32_crude_heavy_peptides_200f_25Feb15_Tiger_15-01-26\MSG201503091410_Auto1169297"
        Const strMSGFPlusFolder = "C:\DMS_WorkDir"

        'Const strXTandemFolder = "\\proto-7\VOrbiETD01\2013_3\QC_Shew_13_04_pt1_1_2_27Jun13_Leopard_13-05-20\XTM201307011524_Auto958319"

        'Const strMSAlignFolder = "\\proto-9\VOrbiETD02\2014_1\Synocho_D2_2\MSA201402281500_Auto1030272"

        Dim strSynOrFHTFile As String
        Dim eMatchedResultType As clsPHRPReader.ePeptideHitResultType

        If False Then
            Console.WriteLine()
            strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder, eMatchedResultType)
            If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> clsPHRPReader.ePeptideHitResultType.Unknown Then
                TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=True)
            End If
        End If

        Dim diMSGFPlusFolder = New DirectoryInfo(strMSGFPlusFolder)
        If Not diMSGFPlusFolder.Exists Then
            Console.WriteLine("Warning, Folder not found: " & strMSGFPlusFolder)
        End If

        If True And diMSGFPlusFolder.Exists Then
            Console.WriteLine()
            strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strMSGFPlusFolder, eMatchedResultType)
            If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> clsPHRPReader.ePeptideHitResultType.Unknown Then
                TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=False)
            End If
        End If

        If True And diMSGFPlusFolder.Exists Then
            Console.WriteLine()
            strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strMSGFPlusFolder, eMatchedResultType)
            If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> clsPHRPReader.ePeptideHitResultType.Unknown Then
                TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=True)
            End If
        End If

        If True And diMSGFPlusFolder.Exists Then
            ' Look for an MSGF+ parameter file to parse
            Dim lstFiles = diMSGFPlusFolder.GetFiles("MSGFDB*.txt")
            If lstFiles.Count > 0 Then
                TestMSGFPlusParamFileParsing(lstFiles.First.FullName)
            End If
        End If

        Dim diSequestFolder = New DirectoryInfo(strSequestFolder)
        If Not diSequestFolder.Exists Then
            Console.WriteLine("Warning, Folder not found: " & strSequestFolder)
        End If

        If True Or Not diSequestFolder.Exists Then Return

        Dim dtStartTimeNoSkipDup As DateTime
        Dim dtEndTimeNoSkipDup As DateTime

        Dim dtStartTimeSkipDup As DateTime
        Dim dtEndTimeSkipDup As DateTime

        Console.WriteLine()
        strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder)
        If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> clsPHRPReader.ePeptideHitResultType.Unknown Then
            dtStartTimeNoSkipDup = DateTime.UtcNow()
            TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=False)
            dtEndTimeNoSkipDup = DateTime.UtcNow()
        End If

        Console.WriteLine()
        strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder, eMatchedResultType)
        If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> clsPHRPReader.ePeptideHitResultType.Unknown Then
            dtStartTimeSkipDup = DateTime.UtcNow()
            TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=True)
            dtEndTimeSkipDup = DateTime.UtcNow()
        End If

        Console.WriteLine()

        Console.WriteLine("Elapsed time (Keep Duplicates): " & dtEndTimeNoSkipDup.Subtract(dtStartTimeNoSkipDup).TotalSeconds.ToString("0.0") & " seconds")
        Console.WriteLine("Elapsed time (Skip Duplicates): " & dtEndTimeSkipDup.Subtract(dtStartTimeSkipDup).TotalSeconds.ToString("0.0") & " seconds")
    End Sub

    Private Sub TestMSGFPlusParamFileParsing(msgfPlusParamFilePath As String)

        Dim localErrorMsg As String = String.Empty
        Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor("MSGF+")

        AddHandler modFileProcessor.ErrorOccurred, AddressOf ModExtractorErrorHandler
        AddHandler modFileProcessor.WarningMessageEvent, AddressOf ModExtractorWarningHandler

        Dim peptideMassCalculator = New clsPeptideMassCalculator()
        Dim success = clsPHRPParserMSGFDB.UpdateMassCalculatorMasses(msgfPlusParamFilePath, modFileProcessor, peptideMassCalculator, localErrorMsg)

        Dim udtModInfo = New List(Of clsPeptideMassCalculator.udtPeptideSequenceModInfoType)

        Dim monoMass = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", udtModInfo)

        Console.WriteLine("Mono mass of PEPTIDE: " & monoMass.ToString("0.0000"))

        Dim udtModifiedResidue = New clsPeptideMassCalculator.udtPeptideSequenceModInfoType
        udtModifiedResidue.ResidueLocInPeptide = 4
        udtModifiedResidue.ModificationMass = 79.966
        udtModInfo.Add(udtModifiedResidue)

        Dim monoMassModified = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", udtModInfo)

        Console.WriteLine("Mono mass of PEPT*IDE: " & monoMassModified.ToString("0.0000"))

        Debug.Assert(Math.Abs(monoMass - 799.359926865) < 0.0000001)

        Debug.Assert(Math.Abs(monoMassModified - 879.325926865) < 0.0000001)
    End Sub

    Private Sub TestPHRPReader(strSynOrFHTFile As String, blnSkipDuplicates As Boolean)

        Dim fiInputFile As FileInfo
        fiInputFile = New FileInfo(strSynOrFHTFile)

        Console.WriteLine("Instantiating reader")
        Dim oStartupOptions = New clsPHRPStartupOptions()
        With oStartupOptions
            .LoadModsAndSeqInfo = True
            .LoadMSGFResults = True
            .LoadScanStatsData = False
            .MaxProteinsPerPSM = 100
        End With

        mPHRPReader = New clsPHRPReader(fiInputFile.FullName, clsPHRPReader.ePeptideHitResultType.Unknown, oStartupOptions)
        mPHRPReader.EchoMessagesToConsole = False
        mPHRPReader.SkipDuplicatePSMs = blnSkipDuplicates

        ' Check for any load errors
        If mPHRPReader.ErrorMessages.Count > 0 Then
            Console.WriteLine("Error(s) instantiating the reader:")
            For Each errorMessage In mPHRPReader.ErrorMessages
                Console.WriteLine("  " & errorMessage)
            Next
        End If

        Const fastReadEnabled = True
        mPHRPReader.FastReadMode = fastReadEnabled

        Dim oMassCalculator = New clsPeptideMassCalculator()

        If Not mPHRPReader.CanRead Then
            Console.WriteLine("Aborting since PHRPReader is not ready: " + mPHRPReader.ErrorMessage)
            Return
        End If

        Dim lstValues = New List(Of String)
        Dim intIsotopeErrorComputed As Integer
        Dim strMassErrorPPM As String

        Dim intPSMsRead = 0
        Dim intModifiedPSMsRead = 0

        Dim dctCachedValues = New Dictionary(Of Integer, clsPSM)

        Console.WriteLine("Reading data")

        Do While mPHRPReader.MoveNext()

            Dim oPsm As clsPSM = mPHRPReader.CurrentPSM

            intPSMsRead += 1
            lstValues.Clear()

            mPHRPReader.FinalizeCurrentPSM()

            Dim primarySequence = String.Empty
            Dim prefix = String.Empty
            Dim suffix = String.Empty
            clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(oPsm.Peptide, primarySequence, prefix, suffix)


            intIsotopeErrorComputed = 0
            strMassErrorPPM = GetCorrectedMassErrorPPM(oPsm, intIsotopeErrorComputed)

            lstValues.Add(mPHRPReader.DatasetName & "_dta.txt")                                             ' #SpecFile
            lstValues.Add("index=" & intPSMsRead)                                                           ' SpecID
            lstValues.Add(oPsm.ScanNumber.ToString())                                                       ' ScanNum
            lstValues.Add(oPsm.CollisionMode)                                                               ' FragMethod
            lstValues.Add(oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge).ToString())      ' Precursor m/z

            lstValues.Add(strMassErrorPPM)                                                                  ' PrecursorError(ppm)
            lstValues.Add(oPsm.Charge.ToString())                                                           ' Charge
            lstValues.Add(oPsm.NumTrypticTerminii.ToString())                                               ' Tryptic state (0, 1, or 2)
            lstValues.Add(CleanupPeptide(oPsm.PeptideWithNumericMods))                                      ' Peptide

            If oPsm.SeqID <= 0 Then
                lstValues.Add("**" & oPsm.SeqID & "**")                                                         ' SeqID is undefined
            Else
                lstValues.Add(oPsm.SeqID.ToString())                                                            ' SeqID
            End If

            lstValues.Add(oPsm.ProteinFirst)                                                                ' Protein First

            If oPsm.ProteinDetails.Count > 0 Then
                Dim oFirstProteinDetail = oPsm.ProteinDetails.First                                         ' Protein Details first

                If Not String.Equals(oPsm.ProteinFirst, oFirstProteinDetail.Key) Then
                    lstValues.Add(oFirstProteinDetail.Key)
                Else
                    lstValues.Add("<Match>")
                End If
                lstValues.Add(oFirstProteinDetail.Value.ResidueStart.ToString())
                lstValues.Add(oFirstProteinDetail.Value.ResidueEnd.ToString())
            End If

            Dim strXCorr = GetScore(oPsm, clsPHRPParserSequest.DATA_COLUMN_XCorr, "0")
            lstValues.Add(strXCorr)                                                             ' XCorr

            lstValues.Add(GetScore(oPsm, clsPHRPParserSequest.DATA_COLUMN_Sp, "0"))              ' SP
            lstValues.Add(oPsm.MSGFSpecProb)                                                                ' MSGF SpecProb
            lstValues.Add(GetScore(oPsm, clsPHRPParserSequest.DATA_COLUMN_DelCn2, "0"))          ' DelCn2

            lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_PValue, "0"))           ' PValue
            lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_EValue, "0"))           ' EValue
            lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue, "0"))           ' SpecEValue
            lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_FDR, "1"))              ' FDR

            If oPsm.PeptideCleanSequence = "QQIEESTSDYDKEK" Then
                Console.WriteLine(oPsm.Peptide & " in scan " & oPsm.ScanNumber)

                Dim parentIonMZ = oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge)

                Console.WriteLine("ParentIonMZ   = " & parentIonMZ)
                Console.WriteLine("PeptideWithNumericMods   = " & oPsm.PeptideWithNumericMods)
            End If

            If oPsm.ModifiedResidues.Count > 0 Then
                intModifiedPSMsRead += 1

                If intModifiedPSMsRead Mod 500 = 0 Then
                    Console.WriteLine("PeptideWithNumericMods   = " & oPsm.PeptideWithNumericMods)
                    For Each modifiedResidue In oPsm.ModifiedResidues
                        Console.WriteLine("  " & modifiedResidue.Residue & modifiedResidue.EndResidueLocInPeptide & ": " & modifiedResidue.ModDefinition.ModificationMassAsText)
                    Next
                End If

                Dim dblPeptideMassRecomputed = oMassCalculator.ComputeSequenceMassNumericMods(oPsm.PeptideWithNumericMods)
                If Math.Abs(oPsm.PeptideMonoisotopicMass - dblPeptideMassRecomputed) > 0.1 Then
                    Console.WriteLine("  Peptide mass disagreement: " & (oPsm.PeptideMonoisotopicMass - dblPeptideMassRecomputed).ToString("0.0000000"))
                End If
            End If


            Dim strFlattened = FlattenList(lstValues)

            If intPSMsRead Mod 1000 = 0 Then
                'Console.WriteLine(intPSMsRead.ToString().PadRight(8) & " " & oPsm.Peptide.PadRight(40) & "   " & strXCorr)
                Console.WriteLine(strFlattened)
            End If

            dctCachedValues.Add(intPSMsRead, oPsm)

        Loop

    End Sub

    Private Function CleanupPeptide(strPeptide As String) As String

        Static reFindItraq As Regex = New Regex("^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", RegexOptions.Compiled Or RegexOptions.IgnoreCase)

        Dim strPrimarySequence = String.Empty
        Dim strPrefix = String.Empty
        Dim strSuffix = String.Empty

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

    Private Function GetCorrectedMassErrorPPM(oPsm As clsPSM, ByRef intIsotopeError As Integer) As String

        Const MASS_C13 = 1.00335483

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

    Private Function GetScore(oPsm As clsPSM, strScoreName As String, strValueIfMissing As String) As String
        Dim strScoreValue = String.Empty

        If Not oPsm.TryGetScore(strScoreName, strScoreValue) Then
            strScoreValue = strValueIfMissing
        End If

        Return strScoreValue

    End Function

    Private Sub mPHRPReader_ErrorEvent(strErrorMessage As String) Handles mPHRPReader.ErrorEvent
        Console.WriteLine("Error: " & strErrorMessage)
    End Sub

    Private Sub mPHRPReader_MessageEvent(strMessage As String) Handles mPHRPReader.MessageEvent
        Console.WriteLine(strMessage)
    End Sub

    Private Sub mPHRPReader_WarningEvent(strWarningMessage As String) Handles mPHRPReader.WarningEvent
        Console.WriteLine("Warning: " & strWarningMessage)
    End Sub


#Region "Event Handlers"
    Private Sub ModExtractorErrorHandler(errMsg As String)
        Console.WriteLine("Error: " & errMsg)
    End Sub

    Private Sub ModExtractorWarningHandler(warningMsg As String)
        Console.WriteLine("Warning: " & warningMsg)
    End Sub
#End Region
End Module