// This program processes search results from several LC-MS/MS search engines to
// determine the modifications present, determine the cleaveage and terminus state
// of each peptide, and compute the monoisotopic mass of each peptide. See
// clsSequestSynopsisFileProcessor and clsXTandemResultsConverter for
// additional information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/ or https://panomics.pnnl.gov/
// -------------------------------------------------------------------------------
//

using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Reflection;
using System.Threading;
using PRISM;

namespace PeptideHitResultsProcRunner
{
    static class Program
    {
        public const string PROGRAM_DATE = "September 20, 2018";

        private static string mInputFilePath;
        private static string mOutputFolderPath;                         // Optional
        private static string mParameterFilePath;                        // Optional

        private static string mMassCorrectionTagsFilePath;               // Optional
        private static string mModificationDefinitionsFilePath;          // Optional
        private static string mSearchToolParameterFilePath;              // Optional

        // Note: If this is true and the _PepToProtMap.txt file isn't found then it will be created using the the Fasta file specified by mFastaFilePath
        private static bool mCreateProteinModsFile;
        private static string mFastaFilePath;
        private static bool mIgnorePeptideToProteinMapperErrors;
        private static bool mProteinModsFileIncludesReversedProteins;
        private static bool mUseExistingMTSPepToProteinMapFile;

        // Setting this to true assumes the input file is a valid PHRP data file
        // Consequently, the code will only try to create the _ProteinMods.txt file, it will not re-create the PHRP data files
        private static bool mCreateProteinModsUsingPHRPDataFile;

        private static bool mCreateInspectOrMSGFDBFirstHitsFile;
        private static bool mCreateInspectOrMSGFDBSynopsisFile;

        private static float mMsgfPlusEValueThreshold;
        private static float mMsgfPlusSpecEValueThreshold;
        private static float mInspectSynopsisFilePValueThreshold;

        private static float mMODaMODPlusSynopsisFileProbabilityThreshold;

        private static string mOutputFolderAlternatePath;                // Optional
        private static bool mRecreateFolderHierarchyInAlternatePath;     // Optional

        private static bool mRecurseFolders;
        private static int mRecurseFoldersMaxLevels;

        private static bool mLogMessagesToFile;
        private static string mLogFilePath = string.Empty;
        private static string mLogFolderPath = string.Empty;

        private static clsPeptideHitResultsProcRunner mPeptideHitResultsProcRunner;

        private static DateTime mLastProgressReportTime;
        private static int mLastProgressReportValue;
        private static DateTime mLastProgressReportValueTime;

        public static int Main()
        {
            // Returns 0 if no error, error code if an error

            int returnCode;
            var parseCommandLine = new clsParseCommandLine();

            mInputFilePath = string.Empty;
            mOutputFolderPath = string.Empty;
            mParameterFilePath = string.Empty;

            mMassCorrectionTagsFilePath = string.Empty;
            mModificationDefinitionsFilePath = string.Empty;
            mSearchToolParameterFilePath = string.Empty;

            mCreateProteinModsFile = false;
            mFastaFilePath = string.Empty;
            mIgnorePeptideToProteinMapperErrors = false;
            mProteinModsFileIncludesReversedProteins = false;
            mUseExistingMTSPepToProteinMapFile = false;

            mCreateProteinModsUsingPHRPDataFile = false;

            mMsgfPlusEValueThreshold = PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            mMsgfPlusSpecEValueThreshold = PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            // These should default to True
            mCreateInspectOrMSGFDBFirstHitsFile = true;
            mCreateInspectOrMSGFDBSynopsisFile = true;
            mInspectSynopsisFilePValueThreshold = PeptideHitResultsProcessor.clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            mMODaMODPlusSynopsisFileProbabilityThreshold = PeptideHitResultsProcessor.clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            mRecurseFoldersMaxLevels = 0;

            mLogMessagesToFile = false;
            mLogFilePath = string.Empty;
            mLogFolderPath = string.Empty;

            try
            {
                var proceed = false;
                if (parseCommandLine.ParseCommandLine())
                {
                    if (SetOptionsUsingCommandLineParameters(parseCommandLine))
                        proceed = true;
                }

                if (!proceed ||
                    parseCommandLine.NeedToShowHelp ||
                    parseCommandLine.ParameterCount + parseCommandLine.NonSwitchParameterCount == 0 ||
                    string.IsNullOrWhiteSpace(mInputFilePath))
                {
                    ShowProgramHelp();
                    returnCode = -1;
                }
                else
                {
                    // Note: Most of the options will get overridden if defined in the parameter file
                    mPeptideHitResultsProcRunner = new clsPeptideHitResultsProcRunner
                    {
                        LogMessagesToFile = mLogMessagesToFile,
                        LogFilePath = mLogFilePath,
                        LogFolderPath = mLogFolderPath,
                        MassCorrectionTagsFilePath = mMassCorrectionTagsFilePath,
                        ModificationDefinitionsFilePath = mModificationDefinitionsFilePath,
                        SearchToolParameterFilePath = mSearchToolParameterFilePath,
                        WarnMissingParameterFileSection = true,
                        CreateProteinModsFile = mCreateProteinModsFile,
                        FastaFilePath = mFastaFilePath,
                        IgnorePeptideToProteinMapperErrors = mIgnorePeptideToProteinMapperErrors,
                        ProteinModsFileIncludesReversedProteins = mProteinModsFileIncludesReversedProteins,
                        UseExistingMTSPepToProteinMapFile = mUseExistingMTSPepToProteinMapFile,
                        CreateProteinModsUsingPHRPDataFile = mCreateProteinModsUsingPHRPDataFile,
                        MsgfPlusEValueThreshold = mMsgfPlusEValueThreshold,
                        MsgfPlusSpecEValueThreshold = mMsgfPlusSpecEValueThreshold,
                        CreateInspectOrMSGFDbFirstHitsFile = mCreateInspectOrMSGFDBFirstHitsFile,
                        CreateInspectOrMSGFDbSynopsisFile = mCreateInspectOrMSGFDBSynopsisFile,
                        InspectSynopsisFilePValueThreshold = mInspectSynopsisFilePValueThreshold,
                        MODaMODPlusSynopsisFileProbabilityThreshold = mMODaMODPlusSynopsisFileProbabilityThreshold
                    };

                    mPeptideHitResultsProcRunner.ErrorEvent += mPeptideHitResultsProcRunner_ErrorEvent;
                    mPeptideHitResultsProcRunner.StatusEvent += mPeptideHitResultsProcRunner_MessageEvent;
                    mPeptideHitResultsProcRunner.ProgressUpdate += mPeptideHitResultsProcRunner_ProgressChanged;
                    mPeptideHitResultsProcRunner.ProgressReset += mPeptideHitResultsProcRunner_ProgressReset;
                    mPeptideHitResultsProcRunner.WarningEvent += mPeptideHitResultsProcRunner_WarningEvent;

                    if (mRecurseFolders)
                    {
                        if (mPeptideHitResultsProcRunner.ProcessFilesAndRecurseFolders(mInputFilePath, mOutputFolderPath, mOutputFolderAlternatePath, mRecreateFolderHierarchyInAlternatePath, mParameterFilePath, mRecurseFoldersMaxLevels))
                        {
                            returnCode = 0;
                        }
                        else
                        {
                            returnCode = (int)mPeptideHitResultsProcRunner.ErrorCode;
                        }
                    }
                    else
                    {
                        if (mPeptideHitResultsProcRunner.ProcessFilesWildcard(mInputFilePath, mOutputFolderPath, mParameterFilePath))
                        {
                            returnCode = 0;
                        }
                        else
                        {
                            returnCode = (int)mPeptideHitResultsProcRunner.ErrorCode;
                            if (returnCode == 0)
                            {
                                returnCode = -1;
                                ShowErrorMessage("ProcessFilesWildcard returned Success=False");
                            }
                            else
                            {
                                ShowErrorMessage("Error while processing: " + mPeptideHitResultsProcRunner.GetErrorMessage());
                            }
                        }
                    }

                    DisplayProgressPercent(mLastProgressReportValue, true);
                }
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in modMain->Main: " + Environment.NewLine + ex.Message);
                returnCode = -1;
            }

            return returnCode;
        }

        private static void DisplayProgressPercent(int percentComplete, bool addCarriageReturn)
        {
            if (addCarriageReturn)
            {
                Console.WriteLine();
            }
            if (percentComplete > 100)
                percentComplete = 100;
            Console.Write("Processing: " + percentComplete + "% ");

            if (addCarriageReturn)
            {
                Console.WriteLine();
            }
        }

        private static string GetAppVersion()
        {
            return PRISM.FileProcessor.ProcessFilesOrFoldersBase.GetAppVersion(PROGRAM_DATE);
        }

        /// <summary>
        /// Parse out True/False or Yes/No or T/F or Y/N or 1/0 from value
        /// </summary>
        /// <param name="valueText">Text to parse</param>
        /// <param name="value">Output parameter</param>
        /// <returns>True if successfully parsed value; the result of the parse is in value</returns>
        /// <remarks></remarks>
        private static bool ParseBoolean(string valueText, out bool value)
        {
            if (string.IsNullOrEmpty(valueText))
            {
                value = false;
                return false;
            }

            if (bool.TryParse(valueText, out value))
            {
                return true;
            }

            switch (valueText.ToUpper()[0])
            {
                case 'T':
                case 'Y':
                case '1':
                    // True or Yes or 1
                    value = true;
                    return true;
                case 'F':
                case 'N':
                case '0':
                    // False or No or 0
                    value = false;
                    return true;
            }

            return false;
        }

        private static bool SetOptionsUsingCommandLineParameters(clsParseCommandLine parseCommandLine)
        {
            // Returns True if no problems; otherwise, returns false

            var validParameters = new List<string> { "I", "O", "Folder", "P", "M", "T", "N", "ProteinMods",
                "F", "Fasta", "IgnorePepToProtMapErrors", "ProteinModsViaPHRP", "ProteinModsIncludeReversed",
                "MSGFPlusEValue", "MSGFPlusSpecEValue", "SynPvalue", "InsFHT", "InsSyn", "SynProb", "S", "A",
                "R", "L" };

            try
            {
                // Make sure no invalid parameters are present
                if (parseCommandLine.InvalidParametersPresent(validParameters))
                {
                    ShowErrorMessage("Invalid command line parameters",
                        (from item in parseCommandLine.InvalidParameters(validParameters) select "/" + item).ToList());
                    return false;
                }

                // Query parseCommandLine to see if various parameters are present
                if (parseCommandLine.RetrieveValueForParameter("I", out var value))
                {
                    mInputFilePath = string.Copy(value);
                }
                else if (parseCommandLine.NonSwitchParameterCount > 0)
                {
                    mInputFilePath = parseCommandLine.RetrieveNonSwitchParameter(0);
                }

                if (parseCommandLine.RetrieveValueForParameter("O", out value))
                    mOutputFolderPath = string.Copy(value);

                // Future
                // If .RetrieveValueForParameter("Folder", value) Then mDatasetFolderPath = String.Copy(value)
                //

                if (parseCommandLine.RetrieveValueForParameter("P", out value))
                    mParameterFilePath = string.Copy(value);
                if (parseCommandLine.RetrieveValueForParameter("M", out value))
                    mModificationDefinitionsFilePath = string.Copy(value);
                if (parseCommandLine.RetrieveValueForParameter("T", out value))
                    mMassCorrectionTagsFilePath = string.Copy(value);
                if (parseCommandLine.RetrieveValueForParameter("N", out value))
                    mSearchToolParameterFilePath = string.Copy(value);

                if (parseCommandLine.IsParameterPresent("ProteinMods"))
                {
                    mCreateProteinModsFile = true;
                }

                if (parseCommandLine.IsParameterPresent("ProteinModsViaPHRP"))
                {
                    mCreateProteinModsUsingPHRPDataFile = true;
                }

                if (parseCommandLine.RetrieveValueForParameter("F", out value))
                    mFastaFilePath = string.Copy(value);
                if (parseCommandLine.RetrieveValueForParameter("Fasta", out value))
                    mFastaFilePath = string.Copy(value);

                if (parseCommandLine.IsParameterPresent("IgnorePepToProtMapErrors"))
                    mIgnorePeptideToProteinMapperErrors = true;
                if (parseCommandLine.IsParameterPresent("ProteinModsIncludeReversed"))
                    mProteinModsFileIncludesReversedProteins = true;
                if (parseCommandLine.IsParameterPresent("UseExistingPepToProteinMapFile"))
                    mUseExistingMTSPepToProteinMapFile = true;

                if (parseCommandLine.RetrieveValueForParameter("InsFHT", out value))
                {
                    if (ParseBoolean(value, out var blnValue))
                    {
                        mCreateInspectOrMSGFDBFirstHitsFile = blnValue;
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("InsSyn", out value))
                {
                    if (ParseBoolean(value, out var blnValue ))
                    {
                        mCreateInspectOrMSGFDBSynopsisFile = blnValue;
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("MSGFPlusEValue", out value))
                {
                    if (float.TryParse(value, out var floatValue))
                    {
                        mMsgfPlusEValueThreshold = floatValue;
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("MSGFPlusSpecEValue", out value))
                {
                    if (float.TryParse(value, out var floatValue))
                    {
                        mMsgfPlusSpecEValueThreshold = floatValue;
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("SynPvalue", out value))
                {
                    if (float.TryParse(value, out var floatValue))
                    {
                        mInspectSynopsisFilePValueThreshold = floatValue;
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("SynProb", out value))
                {
                    if (float.TryParse(value, out var floatValue))
                    {
                        mMODaMODPlusSynopsisFileProbabilityThreshold = floatValue;
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("S", out value))
                {
                    mRecurseFolders = true;
                    if (int.TryParse(value, out var intValue))
                    {
                        mRecurseFoldersMaxLevels = intValue;
                    }
                }
                if (parseCommandLine.RetrieveValueForParameter("A", out value))
                    mOutputFolderAlternatePath = string.Copy(value);
                if (parseCommandLine.IsParameterPresent("R"))
                    mRecreateFolderHierarchyInAlternatePath = true;

                if (parseCommandLine.RetrieveValueForParameter("L", out value))
                {
                    mLogMessagesToFile = true;

                    if (!string.IsNullOrEmpty(value))
                    {
                        mLogFilePath = string.Copy(value).Trim('"');
                    }
                }

                if (parseCommandLine.RetrieveValueForParameter("LogFolder", out value))
                {
                    mLogMessagesToFile = true;
                    if (!string.IsNullOrEmpty(value))
                    {
                        mLogFolderPath = string.Copy(value);
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error parsing the command line parameters: " + Environment.NewLine + ex.Message);
            }

            return false;
        }

        private static void ShowErrorMessage(string message, Exception ex = null)
        {
            ConsoleMsgUtils.ShowError(message, ex);
        }

        private static void ShowErrorMessage(string title, IEnumerable<string> errorMessages)
        {
            ConsoleMsgUtils.ShowErrors(title, errorMessages);
        }

        private static void ShowProgramHelp()
        {
            try
            {
                Console.WriteLine("This program reads in an XTandem results file (XML format), Sequest Synopsis/First Hits file, Inspect search result file, MSGF+ search result file, or MSAlign results file then creates a tab-delimited text file with the data in a standard format used at PNNL.");
                Console.WriteLine("It will insert modification symbols into the peptide sequences for modified peptides. Parallel files will be created containing sequence info and modification details.");
                Console.WriteLine("The user can optionally provide a modification definition file which specifies the symbol to use for each modification mass.");
                Console.WriteLine();
                Console.WriteLine("Program syntax:" + Environment.NewLine +
                                  Path.GetFileName(Assembly.GetExecutingAssembly().Location) +
                                  " InputFilePath [/O:OutputFolderPath]");
                // Future:
                // Console.WriteLine(" [/Folder:DatasetFolderPath]")
                //
                Console.WriteLine(" [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath]");
                Console.WriteLine(" [/ProteinMods] [/F:FastaFilePath] [/ProteinModsViaPHRP] [/IgnorePepToProtMapErrors]");
                Console.WriteLine(" [/ProteinModsIncludeReversed] [/UseExistingPepToProteinMapFile]");
                Console.WriteLine(" [/T:MassCorrectionTagsFilePath] [/N:SearchToolParameterFilePath]");
                Console.WriteLine(" [/MSGFPlusSpecEValue:0.0000005] [/MSGFPlusEValue:0.75]");
                Console.WriteLine(" [/SynPvalue:0.2] [/InsFHT:True|False] [/InsSyn:True|False]");
                Console.WriteLine(" [/SynProb:0.05]");
                Console.WriteLine(" [/S:[MaxLevel]] [/A:AlternateOutputFolderPath] [/R] [/L:[LogFilePath]]");
                Console.WriteLine();
                Console.WriteLine("The input file should be an XTandem Results file (_xt.xml), a Sequest Synopsis File (_syn.txt), a Sequest First Hits file (_fht.txt), an Inspect results file (_inspect.txt), an MSGF-DB results file (_msgfdb.txt), an MSGF+ results file (_msgfdb.tsv or _msgfplus.tsv), or an MSAlign results files (_MSAlign_ResultTable.txt)");
                Console.WriteLine("The output folder switch is optional. If omitted, the output file will be created in the same folder as the input file.");
                Console.WriteLine();
                // Future:
                // Console.WriteLine("As an alternative to specifying an input file, you can specify an input folder. In this case the program will look for the best file to process from that folder, and will auto-determine /T and /N")
                //
                Console.WriteLine();
                Console.WriteLine("The parameter file path is optional. If included, it should point to a valid XML parameter file.");
                Console.WriteLine();
                Console.WriteLine("Use /M to specify the file containing the modification definitions. This file should be tab delimited, with the first column containing the modification symbol, the second column containing the modification mass, plus optionally a third column listing the residues that can be modified with the given mass (1 letter residue symbols, no need to separated with commas or spaces).");
                Console.WriteLine();
                Console.WriteLine("Use /ProteinMods to indicate that the _ProteinMods.txt file should be created. This requires that either an existing _PepToProtMapMTS.txt file exist, or that the Fasta file be defined using /F");
                Console.WriteLine("Use /ProteinModsViaPHRP to indicate that InputFilePath specifies a valid PHRP data file and thus the PHRP data files should not be re-created; only the _ProteinMods.txt file should be created. This requires that either an existing _PepToProtMapMTS.txt file exist, or that the Fasta file be defined using /F");
                Console.WriteLine("Use /F to specify the path to the fasta file. When provided, the order of the proteins in the FASTA file dictates which protein is listed for each peptide in the First Hits file");
                Console.WriteLine();
                Console.WriteLine("Use /IgnorePepToProtMapErrors to ignore peptide to protein mapping errors that occur when creating a missing _PepToProtMapMTS.txt file");
                Console.WriteLine("Use /ProteinModsIncludeReversed to include Reversed proteins in the _ProteinMods.txt file");
                Console.WriteLine("Use /UseExistingPepToProteinMapFile to use an existing _PepToProtMapMTS.txt file if it exists");
                Console.WriteLine();
                Console.WriteLine("Use /T to specify the file containing the mass correction tag info. This file should be tab delimited, with the first column containing the mass correction tag name and the second column containing the mass (the name cannot contain commas or colons and can be, at most, 8 characters long).");
                Console.WriteLine("Use /N to specify the parameter file provided to the search tool. This is only used when processing Inspect or MSGF+ files.");
                Console.WriteLine();

                Console.WriteLine("When processing an MSGF+ results file, use /MSGFPlusSpecEValue and /MSGFPlusEValue to customize the thresholds used to determine which peptides are written to the the synopsis file");
                Console.WriteLine("Defaults are /MSGFPlusSpecEValue:" + PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD +
                                  " and /MSGFPlusEValue:" + PeptideHitResultsProcessor.clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD);
                Console.WriteLine();
                Console.WriteLine("When processing an Inspect results file, use /SynPvalue to customize the PValue threshold used to determine which peptides are written to the the synopsis file. The default is /SynPvalue:0.2  Note that peptides with a TotalPRMScore >= " + PeptideHitResultsProcessor.clsInSpecTResultsProcessor.TOTALPRMSCORE_THRESHOLD + " or an FScore >= " + PeptideHitResultsProcessor.clsInSpecTResultsProcessor.FSCORE_THRESHOLD + " will also be included in the synopsis file.");
                Console.WriteLine("Use /InsFHT:True or /InsFHT:False to toggle the creation of a first-hits file (_fht.txt) when processing Inspect or MSGF+ results (default is /InsFHT:True)");
                Console.WriteLine("Use /InsSyn:True or /InsSyn:False to toggle the creation of a synopsis file (_syn.txt) when processing Inspect or MSGF+ results (default is /InsSyn:True)");
                Console.WriteLine();
                Console.WriteLine("When processing a MODPlus or MODa results file, use /SynProb to customize the Probability threshold used to determine which peptides are written to the the synopsis file. The default is /Synprob:0.05");
                Console.WriteLine();
                Console.WriteLine("Use /S to process all valid files in the input folder and subfolders. Include a number after /S (like /S:2) to limit the level of subfolders to examine.");
                Console.WriteLine("When using /S, you can redirect the output of the results using /A.");
                Console.WriteLine("When using /S, you can use /R to re-create the input folder hierarchy in the alternate output folder (if defined).");
                Console.WriteLine();
                Console.WriteLine("Use /L to specify that a log file should be created. Use /L:LogFilePath to specify the name (or full path) for the log file.");
                Console.WriteLine();

                Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006");
                Console.WriteLine("Version: " + GetAppVersion());

                Console.WriteLine();

                Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov");
                Console.WriteLine("Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/");
                Console.WriteLine();

                // Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
                Thread.Sleep(750);
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error displaying the program syntax: " + ex.Message);
            }
        }

        private static void mPeptideHitResultsProcRunner_ErrorEvent(string message, Exception ex)
        {
            ShowErrorMessage(message, ex);
        }

        private static void mPeptideHitResultsProcRunner_MessageEvent(string message)
        {
            Console.WriteLine(message);
        }

        private static void mPeptideHitResultsProcRunner_ProgressChanged(string taskDescription, float percentComplete)
        {
            const int PERCENT_REPORT_INTERVAL = 25;
            const int PROGRESS_DOT_INTERVAL_MSEC = 250;
            const int PROGRESS_VALUE_INTERVAL_SEC = 60;

            if (percentComplete >= mLastProgressReportValue)
            {
                if (mLastProgressReportValue > 0)
                {
                    Console.WriteLine();
                }
                DisplayProgressPercent(mLastProgressReportValue, false);
                mLastProgressReportValue += PERCENT_REPORT_INTERVAL;
                mLastProgressReportTime = DateTime.UtcNow;
            }
            else
            {
                if (DateTime.UtcNow.Subtract(mLastProgressReportTime).TotalMilliseconds > PROGRESS_DOT_INTERVAL_MSEC)
                {
                    mLastProgressReportTime = DateTime.UtcNow;
                    if (DateTime.UtcNow.Subtract(mLastProgressReportValueTime).TotalSeconds > PROGRESS_VALUE_INTERVAL_SEC)
                    {
                        mLastProgressReportValueTime = DateTime.UtcNow;
                        Console.WriteLine();
                        Console.Write(percentComplete.ToString("0.00") + "% complete ");
                    }
                    else
                    {
                        Console.Write(".");
                    }
                }
            }
        }

        private static void mPeptideHitResultsProcRunner_ProgressReset()
        {
            mLastProgressReportTime = DateTime.UtcNow;
            mLastProgressReportValueTime = DateTime.UtcNow;
            mLastProgressReportValue = 0;
        }

        private static void mPeptideHitResultsProcRunner_WarningEvent(string message)
        {
            if (message.StartsWith("Warning", StringComparison.InvariantCultureIgnoreCase))
                ConsoleMsgUtils.ShowWarning(message);
            else
                ConsoleMsgUtils.ShowWarning("Warning: " + message);
        }
    }
}
