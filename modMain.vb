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
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 

Module modMain
	Public Const PROGRAM_DATE As String = "November 27, 2012"

	Private mInputFilePath As String
	Private mOutputFolderName As String							' Optional
	Private mParameterFilePath As String						' Optional

	Private mMassCorrectionTagsFilePath As String				' Optional
	Private mModificationDefinitionsFilePath As String			' Optional
	Private mSearchToolParameterFilePath As String				' Optional

	' Note: If this is true and the _PepToProtMap.txt file isn't found then it will be created using the the Fasta file specified by mFastaFilePath
	Private mCreateProteinModsFile As Boolean
	Private mFastaFilePath As String
	Private mIgnorePeptideToProteinMapperErrors As Boolean
	Private mProteinModsFileIncludesReversedProteins As Boolean
	Private mUseExistingMTSPepToProteinMapFile As Boolean

	' Setting this to true assumes the input file is a valid PHRP data file
	' Consequently, the code will only try to create the _ProteinMods.txt file, it will not re-create the PHRP data files
	Private mCreateProteinModsUsingPHRPDataFile As Boolean

	Private mCreateInspectOrMSGFDBFirstHitsFile As Boolean
	Private mCreateInspectOrMSGFDBSynopsisFile As Boolean

	Private mInspectSynopsisFilePValueThreshold As Single		' Optional

	Private mMSGFDBSynopsisFilePValueThreshold As Single
	Private mMSGFDBSynopsisFileSpecProbThreshold As Single

	Private mOutputFolderAlternatePath As String				' Optional
	Private mRecreateFolderHierarchyInAlternatePath As Boolean	' Optional

	Private mRecurseFolders As Boolean
	Private mRecurseFoldersMaxLevels As Integer

	Private mLogMessagesToFile As Boolean
	Private mLogFilePath As String = String.Empty
	Private mQuietMode As Boolean

	Private WithEvents mPeptideHitResultsProcRunner As clsPeptideHitResultsProcRunner
	Private mLastProgressReportTime As System.DateTime
	Private mLastProgressReportValue As Integer

	Public Function Main() As Integer
		' Returns 0 if no error, error code if an error

		Dim intReturnCode As Integer
		Dim objParseCommandLine As New clsParseCommandLine
		Dim blnProceed As Boolean

		''Dim objTest As New PeptideHitResultsProcessor.clsPeptideCleavageStateCalculator
		''Dim strSeq, strCleanSeq, strPrefix, strSuffix As String
		''strSeq = "A.BCDE.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "-.BCDE.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "A.BCDE.-" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "A.B.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "A.BCDE" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "BCDE.F" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "FA.BCDE.FG" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "BCDE" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "BCDE." : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = ".BCDE" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = ".F." : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "F..E" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "AF..EF" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "AFF..EF" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "E.TGMLTQKFARSLGMLAVDNQARV.." : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "..TGMLTQKFARSLGMLAVDNQARV.R" : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''strSeq = "..TGMLTQKFARSLGMLAVDNQARV.." : objTest.SplitPrefixAndSuffixFromSequence(strSeq, strCleanSeq, strPrefix, strSuffix) : Console.WriteLine(strSeq & " -> " & strPrefix & "." & strCleanSeq & "." & strSuffix)
		''
		''Return 0

		intReturnCode = 0
		mInputFilePath = String.Empty
		mOutputFolderName = String.Empty
		mParameterFilePath = String.Empty

		mMassCorrectionTagsFilePath = String.Empty
		mModificationDefinitionsFilePath = String.Empty
		mSearchToolParameterFilePath = String.Empty

		mCreateProteinModsFile = False
		mFastaFilePath = String.Empty
		mIgnorePeptideToProteinMapperErrors = False
		mProteinModsFileIncludesReversedProteins = False
		mUseExistingMTSPepToProteinMapFile = False

		mCreateProteinModsUsingPHRPDataFile = False

		' These should default to True
		mCreateInspectOrMSGFDBFirstHitsFile = True
		mCreateInspectOrMSGFDBSynopsisFile = True
		mInspectSynopsisFilePValueThreshold = PeptideHitResultsProcessor.clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD

		mMSGFDBSynopsisFilePValueThreshold = PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD
		mMSGFDBSynopsisFileSpecProbThreshold = PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPECPROB_THRESHOLD

		mRecurseFolders = False
		mRecurseFoldersMaxLevels = 0

		mQuietMode = False
		mLogMessagesToFile = False
		mLogFilePath = String.Empty

		Try
			blnProceed = False
			If objParseCommandLine.ParseCommandLine Then
				If SetOptionsUsingCommandLineParameters(objParseCommandLine) Then blnProceed = True
			End If

			If Not blnProceed OrElse _
			   objParseCommandLine.NeedToShowHelp OrElse _
			   objParseCommandLine.ParameterCount + objParseCommandLine.NonSwitchParameterCount = 0 OrElse _
			   mInputFilePath.Length = 0 Then
				ShowProgramHelp()
				intReturnCode = -1
			Else
				mPeptideHitResultsProcRunner = New clsPeptideHitResultsProcRunner

				With mPeptideHitResultsProcRunner
					.ShowMessages = Not mQuietMode
					.LogMessagesToFile = mLogMessagesToFile
					.LogFilePath = mLogFilePath

					' Note: These options will get overridden if defined in the parameter file
					.MassCorrectionTagsFilePath = mMassCorrectionTagsFilePath
					.ModificationDefinitionsFilePath = mModificationDefinitionsFilePath
					.SearchToolParameterFilePath = mSearchToolParameterFilePath

					.WarnMissingParameterFileSection = True

					.CreateProteinModsFile = mCreateProteinModsFile
					.FastaFilePath = mFastaFilePath
					.IgnorePeptideToProteinMapperErrors = mIgnorePeptideToProteinMapperErrors
					.ProteinModsFileIncludesReversedProteins = mProteinModsFileIncludesReversedProteins
					.UseExistingMTSPepToProteinMapFile = mUseExistingMTSPepToProteinMapFile

					.CreateProteinModsUsingPHRPDataFile = mCreateProteinModsUsingPHRPDataFile

					.CreateInspectOrMSGFDbFirstHitsFile = mCreateInspectOrMSGFDBFirstHitsFile
					.CreateInspectOrMSGFDbSynopsisFile = mCreateInspectOrMSGFDBSynopsisFile
					.InspectSynopsisFilePValueThreshold = mInspectSynopsisFilePValueThreshold
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
							ShowErrorMessage("Error while processing: " & mPeptideHitResultsProcRunner.GetErrorMessage())
						End If
					End If
				End If

				DisplayProgressPercent(mLastProgressReportValue, True)
			End If

		Catch ex As Exception
			ShowErrorMessage("Error occurred in modMain->Main: " & System.Environment.NewLine & ex.Message)
			intReturnCode = -1
		End Try

		Return intReturnCode

	End Function

	Private Sub DisplayProgressPercent(ByVal intPercentComplete As Integer, ByVal blnAddCarriageReturn As Boolean)

		If blnAddCarriageReturn Then
			Console.WriteLine()
		End If
		If intPercentComplete > 100 Then intPercentComplete = 100
		Console.Write("Processing: " & intPercentComplete.ToString() & "% ")

		If blnAddCarriageReturn Then
			Console.WriteLine()
		End If
	End Sub

	Private Function GetAppVersion() As String
		'Return System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")"

		Return System.Reflection.Assembly.GetExecutingAssembly().GetName().Version.ToString() & " (" & PROGRAM_DATE & ")"
	End Function

	''' <summary>
	''' Parse out True/False or Yes/No or T/F or Y/N or 1/0 from strValue
	''' </summary>
	''' <param name="strValue">Text to parse</param>
	''' <param name="blnValue">Output parameter</param>
	''' <returns>True if successfully parsed strValue; the result of the parse is in blnValue</returns>
	''' <remarks></remarks>
	Private Function ParseBoolean(ByVal strValue As String, ByRef blnValue As Boolean) As Boolean

		If String.IsNullOrEmpty(strValue) Then Return False

		If Boolean.TryParse(strValue, blnValue) Then
			Return True
		Else
			Select Case strValue.ToUpper().Chars(0)
				Case "T"c, "Y"c, "1"c
					' True or Yes or 1
					blnValue = True
					Return True
				Case "F"c, "N"c, "0"c
					' False or No or 0
					blnValue = False
					Return True
			End Select
		End If

		Return False

	End Function

	Private Function SetOptionsUsingCommandLineParameters(ByVal objParseCommandLine As clsParseCommandLine) As Boolean
		' Returns True if no problems; otherwise, returns false

		Dim strValue As String = String.Empty
		Dim sngValue As Single
		Dim intValue As Integer
		Dim blnValue As Boolean
		Dim strValidParameters() As String = New String() {"I", "O", "P", "M", "T", "N", "ProteinMods", "F", "IgnorePepToProtMapErrors", "ProteinModsViaPHRP", "ProteinModsIncludeReversed", "SynPvalue", "InsFHT", "InsSyn", "S", "A", "R", "L", "Q"}

		Try
			' Make sure no invalid parameters are present
			If objParseCommandLine.InvalidParametersPresent(strValidParameters) Then
				Return False
			Else
				With objParseCommandLine
					' Query objParseCommandLine to see if various parameters are present
					If .RetrieveValueForParameter("I", strValue) Then
						mInputFilePath = String.Copy(strValue)
					ElseIf .NonSwitchParameterCount > 0 Then
						mInputFilePath = .RetrieveNonSwitchParameter(0)
					End If

					If .RetrieveValueForParameter("O", strValue) Then mOutputFolderName = String.Copy(strValue)
					If .RetrieveValueForParameter("P", strValue) Then mParameterFilePath = String.Copy(strValue)
					If .RetrieveValueForParameter("M", strValue) Then mModificationDefinitionsFilePath = String.Copy(strValue)
					If .RetrieveValueForParameter("T", strValue) Then mMassCorrectionTagsFilePath = String.Copy(strValue)
					If .RetrieveValueForParameter("N", strValue) Then mSearchToolParameterFilePath = String.Copy(strValue)

					If .RetrieveValueForParameter("ProteinMods", strValue) Then
						mCreateProteinModsFile = True
					End If

					If .RetrieveValueForParameter("ProteinModsViaPHRP", strValue) Then
						mCreateProteinModsUsingPHRPDataFile = True
					End If

					If .RetrieveValueForParameter("F", strValue) Then mFastaFilePath = String.Copy(strValue)
					If .RetrieveValueForParameter("IgnorePepToProtMapErrors", strValue) Then mIgnorePeptideToProteinMapperErrors = True
					If .RetrieveValueForParameter("ProteinModsIncludeReversed", strValue) Then mProteinModsFileIncludesReversedProteins = True
					If .RetrieveValueForParameter("UseExistingPepToProteinMapFile", strValue) Then mUseExistingMTSPepToProteinMapFile = True

					If .RetrieveValueForParameter("InsFHT", strValue) Then
						If ParseBoolean(strValue, blnValue) Then
							mCreateInspectOrMSGFDBFirstHitsFile = blnValue
						End If
					End If

					If .RetrieveValueForParameter("InsSyn", strValue) Then
						If ParseBoolean(strValue, blnValue) Then
							mCreateInspectOrMSGFDBSynopsisFile = blnValue
						End If
					End If

					If .RetrieveValueForParameter("SynPvalue", strValue) Then
						If Single.TryParse(strValue, sngValue) Then
							mInspectSynopsisFilePValueThreshold = sngValue
						End If
					End If

					If .RetrieveValueForParameter("S", strValue) Then
						mRecurseFolders = True
						If Integer.TryParse(strValue, intValue) Then
							mRecurseFoldersMaxLevels = intValue
						End If
					End If
					If .RetrieveValueForParameter("A", strValue) Then mOutputFolderAlternatePath = String.Copy(strValue)
					If .RetrieveValueForParameter("R", strValue) Then mRecreateFolderHierarchyInAlternatePath = True

					If .RetrieveValueForParameter("L", strValue) Then
						mLogMessagesToFile = True

						If Not strValue Is Nothing AndAlso strValue.Length > 0 Then
							mLogFilePath = String.Copy(strValue).Trim(""""c)
						End If
					End If

					If .RetrieveValueForParameter("Q", strValue) Then mQuietMode = True
				End With

				Return True
			End If

		Catch ex As Exception
			ShowErrorMessage("Error parsing the command line parameters: " & System.Environment.NewLine & ex.Message)
		End Try

	End Function

	Private Sub ShowErrorMessage(ByVal strMessage As String)
		Dim strSeparator As String = "------------------------------------------------------------------------------"

		Console.WriteLine()
		Console.WriteLine(strSeparator)
		Console.WriteLine(strMessage)
		Console.WriteLine(strSeparator)
		Console.WriteLine()

		WriteToErrorStream(strMessage)
	End Sub

	Private Sub ShowProgramHelp()

		Try

			Console.WriteLine("This program reads in an XTandem results file (XML format), Sequest Synopsis/First Hits file, Inspect search result file, MSGF-DB search result file, or MSGF+ search result file, then creates a tab-delimited text file with the data in a standard format used at PNNL.  ")
			Console.WriteLine("It will insert modification symbols into the peptide sequences for modified peptides.  Parallel files will be created containing sequence info and modification details.  ")
			Console.WriteLine("The user can optionally provide a modification definition file which specifies the symbol to use for each modification mass.")
			Console.WriteLine()
			Console.WriteLine("Program syntax:" & ControlChars.NewLine & System.IO.Path.GetFileName(System.Reflection.Assembly.GetExecutingAssembly().Location) & _
										" InputFilePath [/O:OutputFolderPath]")
			Console.WriteLine(" [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath]")
			Console.WriteLine(" [/ProteinMods] [/F:FastaFilePath] [/ProteinModsViaPHRP] [/IgnorePepToProtMapErrors]")
			Console.WriteLine(" [/ProteinModsIncludeReversed] [/UseExistingPepToProteinMapFile]")
			Console.WriteLine(" [/T:MassCorrectionTagsFilePath] [/N:SearchToolParameterFilePath] [/SynPvalue:0.2]")
			Console.WriteLine(" [/InsFHT:True|False] [/InsSyn:True|False]")
			Console.WriteLine(" [/S:[MaxLevel]] [/A:AlternateOutputFolderPath] [/R] [/L:[LogFilePath]] [/Q]")
			Console.WriteLine()
			Console.WriteLine("The input file should be an XTandem Results file (_xt.xml), a Sequest Synopsis File (_syn.txt), a Sequest First Hits file (_fht.txt), an Inspect results file (_inspect.txt), an MSGF-DB results file (_msgfdb.txt), or an MSGF+ results file (_msgfdb.tsv or _msgfplus.tsv)")
			Console.WriteLine("The output folder switch is optional.  If omitted, the output file will be created in the same folder as the input file.")
			Console.WriteLine("The parameter file path is optional.  If included, it should point to a valid XML parameter file.")
			Console.WriteLine()
			Console.WriteLine("Use /M to specify the file containing the modification definitions.  This file should be tab delimited, with the first column containing the modification symbol, the second column containing the modification mass, plus optionally a third column listing the residues that can be modified with the given mass (1 letter residue symbols, no need to separated with commas or spaces).")
			Console.WriteLine()
			Console.WriteLine("Use /ProteinMods to indicate that the _ProteinMods.txt file should be created.  This requires that either an existing _PepToProtMapMTS.txt file exist, or that the Fasta file be defined using /F")
			Console.WriteLine("Use /ProteinModsViaPHRP to indicate that InputFilePath specifies a valid PHRP data file and thus the PHRP data files should not be re-created; only the _ProteinMods.txt file should be created.  This requires that either an existing _PepToProtMapMTS.txt file exist, or that the Fasta file be defined using /F")
			Console.WriteLine("Use /IgnorePepToProtMapErrors to ignore peptide to protein mapping errors that occur when creating a missing _PepToProtMapMTS.txt file")
			Console.WriteLine("Use /ProteinModsIncludeReversed to include Reversed proteins in the _ProteinMods.txt file")
			Console.WriteLine("Use /UseExistingPepToProteinMapFile to use an existing _PepToProtMapMTS.txt file if it exists")
			Console.WriteLine()
			Console.WriteLine("Use /T to specify the file containing the mass correction tag info.  This file should be tab delimited, with the first column containing the mass correction tag name and the second column containing the mass (the name cannot contain commas or colons and can be, at most, 8 characters long).")
			Console.WriteLine("Use /N to specify the parameter file provided to the search tool.  This is only used when processing Inspect or MSGF-DB files.")
			Console.WriteLine()
			Console.WriteLine("When processing an Inspect results file, use /SynPvalue to customize the PValue threshold used to determine which peptides are written to the the synopsis file.  The default is /SynPvalue:0.2  Note that peptides with a TotalPRMScore >= " & PeptideHitResultsProcessor.clsInSpecTResultsProcessor.TOTALPRMSCORE_THRESHOLD.ToString() & " or an FScore >= " & PeptideHitResultsProcessor.clsInSpecTResultsProcessor.FSCORE_THRESHOLD & " will also be included in the synopsis file.")
			Console.WriteLine("Use /InsFHT:True or /InsFHT:False to toggle the creation of a first-hits file (_fht.txt) when processing Inspect or MSGF-DB results (default is /InsFHT:True)")
			Console.WriteLine("Use /InsSyn:True or /InsSyn:False to toggle the creation of a synopsis file (_syn.txt) when processing Inspect or MSGF-DB results (default is /InsSyn:True)")
			Console.WriteLine()
			Console.WriteLine("Use /S to process all valid files in the input folder and subfolders. Include a number after /S (like /S:2) to limit the level of subfolders to examine.")
			Console.WriteLine("When using /S, you can redirect the output of the results using /A.")
			Console.WriteLine("When using /S, you can use /R to re-create the input folder hierarchy in the alternate output folder (if defined).")
			Console.WriteLine()
			Console.WriteLine("Use /L to specify that a log file should be created.  Use /L:LogFilePath to specify the name (or full path) for the log file.")
			Console.WriteLine("Use the optional /Q switch will suppress all error messages.")
			Console.WriteLine()

			Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006")
			Console.WriteLine("Version: " & GetAppVersion())

			Console.WriteLine()

			Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com")
			Console.WriteLine("Website: http://panomics.pnnl.gov/ or http://www.sysbio.org/resources/staff/")

			' Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
			System.Threading.Thread.Sleep(750)

		Catch ex As Exception
			ShowErrorMessage("Error displaying the program syntax: " & ex.Message)
		End Try

	End Sub

	Private Sub WriteToErrorStream(strErrorMessage As String)
		Try
			Using swErrorStream As System.IO.StreamWriter = New System.IO.StreamWriter(Console.OpenStandardError())
				swErrorStream.WriteLine(strErrorMessage)
			End Using
		Catch ex As Exception
			' Ignore errors here
		End Try
	End Sub

	Private Sub mPeptideHitResultsProcRunner_ErrorEvent(strMessage As String) Handles mPeptideHitResultsProcRunner.ErrorEvent
		WriteToErrorStream(strMessage)
	End Sub

	Private Sub mPeptideHitResultsProcRunner_MessageEvent(strMessage As String) Handles mPeptideHitResultsProcRunner.MessageEvent
		Console.WriteLine(strMessage)
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
			mLastProgressReportTime = DateTime.UtcNow
		Else
			If DateTime.UtcNow.Subtract(mLastProgressReportTime).TotalMilliseconds > PROGRESS_DOT_INTERVAL_MSEC Then
				mLastProgressReportTime = DateTime.UtcNow
				Console.Write(".")
			End If
		End If
	End Sub

	Private Sub mPeptideHitResultsProcRunner_ProgressReset() Handles mPeptideHitResultsProcRunner.ProgressReset
		mLastProgressReportTime = DateTime.UtcNow
		mLastProgressReportValue = 0
	End Sub

	Private Sub mPeptideHitResultsProcRunner_WarningEvent(strMessage As String) Handles mPeptideHitResultsProcRunner.WarningEvent
		Console.WriteLine("Warning: " & strMessage)
	End Sub
End Module
