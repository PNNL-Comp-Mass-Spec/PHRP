Module modMain

	Private WithEvents mPHRPReader As PHRPReader.clsPHRPReader

	Public Sub Main()


		'Const strSequestSynFilePath As String = "Seq201304121552_Auto934225\Firestone_Soil_07_18_05APR13_Frodo_12-12-04_syn.txt"
		'Const strMSGFPlusSynFilePath As String = "MSG201304261714_Auto938181\FSFA-299b_25Apr13_Methow_13-02-13_msgfdb_syn.txt"

		Const strSequestFolder As String = "\\proto-7\VOrbi05\2013_2\Firestone_Soil_07_18_05APR13_Frodo_12-12-04\Seq201304121552_Auto934225"
		Const strMSGFPlusFolder As String = "MSG201304261714_Auto938181"

		' Const strHugeResultsFolder As String = "E:\DMS_WorkDir"

		'Const strXTandemFolder As String = "\\proto-7\VOrbiETD01\2013_3\QC_Shew_13_04_pt1_1_2_27Jun13_Leopard_13-05-20\XTM201307011524_Auto958319"

		'Const strMSAlignFolder As String = "\\proto-9\VOrbiETD02\2014_1\Synocho_D2_2\MSA201402281500_Auto1030272"

		Dim strSynOrFHTFile As String
		Dim eMatchedResultType As PHRPReader.clsPHRPReader.ePeptideHitResultType

		Console.WriteLine()
		strSynOrFHTFile = PHRPReader.clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder, eMatchedResultType)
		If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
			TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=True)
		End If

		Console.WriteLine()
		strSynOrFHTFile = PHRPReader.clsPHRPReader.AutoDetermineBestInputFile(strMSGFPlusFolder, eMatchedResultType)
		If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
			TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=False)
		End If

		Console.WriteLine()
		strSynOrFHTFile = PHRPReader.clsPHRPReader.AutoDetermineBestInputFile(strMSGFPlusFolder, eMatchedResultType)
		If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
			TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=True)
		End If

		'Console.WriteLine()
		'strSynOrFHTFile = PHRPReader.clsPHRPReader.AutoDetermineBestInputFile(strXTandemFolder)
		'If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
		'	TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=False)
		'End If


		Dim dtStartTimeNoSkipDup As DateTime
		Dim dtEndTimeNoSkipDup As DateTime

		Dim dtStartTimeSkipDup As DateTime
		Dim dtEndTimeSkipDup As DateTime

		Console.WriteLine()
		strSynOrFHTFile = PHRPReader.clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder)
		If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
			dtStartTimeNoSkipDup = System.DateTime.UtcNow()
			TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=False)
			dtEndTimeNoSkipDup = System.DateTime.UtcNow()
		End If

		Console.WriteLine()
		strSynOrFHTFile = PHRPReader.clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder, eMatchedResultType)
		If Not String.IsNullOrEmpty(strSynOrFHTFile) AndAlso eMatchedResultType <> PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown Then
			dtStartTimeSkipDup = System.DateTime.UtcNow()
			TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates:=True)
			dtEndTimeSkipDup = System.DateTime.UtcNow()
		End If
		
		Console.WriteLine()

		Console.WriteLine("Elapsed time (Keep Duplicates): " & dtEndTimeNoSkipDup.Subtract(dtStartTimeNoSkipDup).TotalSeconds.ToString("0.0") & " seconds")
		Console.WriteLine("Elapsed time (Skip Duplicates): " & dtEndTimeSkipDup.Subtract(dtStartTimeSkipDup).TotalSeconds.ToString("0.0") & " seconds")


	End Sub

	Private Sub TestPHRPReader(ByVal strSequestSynFilePath As String, ByVal blnSkipDuplicates As Boolean)

		Dim fiInputFile As IO.FileInfo
		fiInputFile = New IO.FileInfo(strSequestSynFilePath)

		Console.WriteLine("Instantiating reader")
		Dim oStartupOptions = New PHRPReader.clsPHRPStartupOptions()
		With oStartupOptions
			.LoadModsAndSeqInfo = True
			.LoadMSGFResults = True
			.LoadScanStatsData = False
			.MaxProteinsPerPSM = 100
		End With

		mPHRPReader = New PHRPReader.clsPHRPReader(fiInputFile.FullName, PHRPReader.clsPHRPReader.ePeptideHitResultType.Unknown, oStartupOptions)
		mPHRPReader.EchoMessagesToConsole = False
		mPHRPReader.SkipDuplicatePSMs = blnSkipDuplicates

		Dim oMassCalculator = New PHRPReader.clsPeptideMassCalculator()

		If Not mPHRPReader.CanRead Then
			Console.WriteLine("Aborting since PHRPReader is not ready: " + mPHRPReader.ErrorMessage)
			Return
		End If

		Dim lstValues As Generic.List(Of String) = New Generic.List(Of String)
		Dim intIsotopeErrorComputed As Integer
		Dim strMassErrorPPM As String
		Dim intPSMsRead As Integer = 0

		Dim dctCachedValues As Generic.Dictionary(Of Integer, PHRPReader.clsPSM) = New Generic.Dictionary(Of Integer, PHRPReader.clsPSM)

		Console.WriteLine("Reading data")
		Do While mPHRPReader.MoveNext()

			Dim oPsm As PHRPReader.clsPSM = mPHRPReader.CurrentPSM

			intPSMsRead += 1
			lstValues.Clear()

			intIsotopeErrorComputed = 0
			strMassErrorPPM = GetCorrectedMassErrorPPM(oPsm, intIsotopeErrorComputed)

			lstValues.Add(mPHRPReader.DatasetName & "_dta.txt")												' #SpecFile
			lstValues.Add("index=" & intPSMsRead)															' SpecID
			lstValues.Add(oPsm.ScanNumber.ToString())														' ScanNum
			lstValues.Add(oPsm.CollisionMode)																' FragMethod
			lstValues.Add(PHRPReader.clsPeptideMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge).ToString())		' Precursor m/z

			lstValues.Add(strMassErrorPPM)																	' PrecursorError(ppm)
			lstValues.Add(oPsm.Charge.ToString())															' Charge
			lstValues.Add(oPsm.NumTrypticTerminii.ToString())												' Tryptic state (0, 1, or 2)
			lstValues.Add(CleanupPeptide(oPsm.PeptideWithNumericMods))										' Peptide

			If oPsm.SeqID <= 0 Then
				lstValues.Add("**" & oPsm.SeqID & "**")															' SeqID is undefined
			Else
				lstValues.Add(oPsm.SeqID.ToString())															' SeqID
			End If

			lstValues.Add(oPsm.ProteinFirst)																' Protein First

			If oPsm.ProteinDetails.Count > 0 Then
				Dim oFirstProteinDetail = oPsm.ProteinDetails.First											' Protein Details first

				If Not String.Equals(oPsm.ProteinFirst, oFirstProteinDetail.Key) Then
					lstValues.Add(oFirstProteinDetail.Key)
				Else
					lstValues.Add("<Match>")
				End If
				lstValues.Add(oFirstProteinDetail.Value.ResidueStart.ToString())
				lstValues.Add(oFirstProteinDetail.Value.ResidueEnd.ToString())
			End If

			Dim strXCorr As String = GetScore(oPsm, PHRPReader.clsPHRPParserSequest.DATA_COLUMN_XCorr, "0")
			lstValues.Add(strXCorr)																' XCorr

			lstValues.Add(GetScore(oPsm, PHRPReader.clsPHRPParserSequest.DATA_COLUMN_Sp, "0"))				' SP
			lstValues.Add(oPsm.MSGFSpecProb)																' MSGF SpecProb
			lstValues.Add(GetScore(oPsm, PHRPReader.clsPHRPParserSequest.DATA_COLUMN_DelCn2, "0"))			' DelCn2

			lstValues.Add(GetScore(oPsm, PHRPReader.clsPHRPParserMSGFDB.DATA_COLUMN_PValue, "0"))			' PValue
			lstValues.Add(GetScore(oPsm, PHRPReader.clsPHRPParserMSGFDB.DATA_COLUMN_EValue, "0"))			' EValue
			lstValues.Add(GetScore(oPsm, PHRPReader.clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecEValue, "0"))			' SpecEValue
			lstValues.Add(GetScore(oPsm, PHRPReader.clsPHRPParserMSGFDB.DATA_COLUMN_FDR, "1"))				' FDR

			If oPsm.PeptideCleanSequence = "QQIEESTSDYDKEK" Then
				Console.WriteLine(oPsm.Peptide & " in scan " & oPsm.ScanNumber)

				Dim parentIonMZ = PHRPReader.clsPeptideMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge)

				Console.WriteLine("ParentIonMZ   = " & parentIonMZ)
				Console.WriteLine("PeptideWithNumericMods   = " & oPsm.PeptideWithNumericMods)

			End If

			If oPsm.ModifiedResidues.Count > 0 Then
				Dim dblPeptideMassRecomputed = oMassCalculator.ComputeSequenceMassNumericMods(oPsm.PeptideWithNumericMods)
				If Math.Abs(oPsm.PeptideMonoisotopicMass - dblPeptideMassRecomputed) > 0.1 Then
					Console.WriteLine("  Peptide mass disagreement: " & (oPsm.PeptideMonoisotopicMass - dblPeptideMassRecomputed).ToString("0.0000000"))
				End If
			End If


			Dim strFlattened As String = FlattenList(lstValues)

			If intPSMsRead Mod 10000 = 0 Then
				'Console.WriteLine(intPSMsRead.ToString().PadRight(8) & " " & oPsm.Peptide.PadRight(40) & "   " & strXCorr)
				Console.WriteLine(strFlattened)
			End If

			dctCachedValues.Add(intPSMsRead, oPsm)

		Loop

	End Sub

	Private Function CleanupPeptide(ByVal strPeptide As String) As String

		Static reFindItraq As System.Text.RegularExpressions.Regex = New System.Text.RegularExpressions.Regex("^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.IgnoreCase)

		Dim strPrimarySequence As String = String.Empty
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty

		Dim reMatch As System.Text.RegularExpressions.Match

		If PHRPReader.clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, strPrimarySequence, strPrefix, strSuffix) Then
			' Look for an N-terminal iTraq mod
			reMatch = reFindItraq.Match(strPrimarySequence)

			If reMatch.Success Then
				strPeptide = strPrefix & "." & reMatch.Groups(2).Value & reMatch.Groups(1).Value & reMatch.Groups(3).Value & "." & strSuffix
			End If
		End If

		Return strPeptide
	End Function

	Private Function FlattenList(ByVal lstValues As Generic.List(Of String)) As String
		Return FlattenList(lstValues, ControlChars.Tab)
	End Function

	Private Function FlattenList(ByVal lstValues As Generic.List(Of String), ByVal chSepChar As Char) As String
		Dim sbOutline As System.Text.StringBuilder = New System.Text.StringBuilder()

		For intIndex As Integer = 0 To lstValues.Count - 1
			If intIndex > 0 Then
				sbOutline.Append(chSepChar)
			End If
			sbOutline.Append(lstValues(intIndex))
		Next

		Return sbOutline.ToString()
	End Function

	Private Function GetCorrectedMassErrorPPM(ByVal oPsm As PHRPReader.clsPSM, ByRef intIsotopeError As Integer) As String

		Const MASS_C13 As Double = 1.00335483

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

			dblMassErrorPPM = PHRPReader.clsPeptideMassCalculator.MassToPPM(dblDelM, oPsm.PrecursorNeutralMass)
		End If

		Return dblMassErrorPPM.ToString("0.0000")

	End Function

	Private Function GetScore(ByVal oPsm As PHRPReader.clsPSM, ByVal strScoreName As String, ByVal strValueIfMissing As String) As String
		Dim strScoreValue As String = String.Empty

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

End Module