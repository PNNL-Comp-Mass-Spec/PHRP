// This program converts search results from various MS/MS identification tools into a series
// of tab-delimited text files that organize the data in a similar format for each tool.
// It supports MS-GF+, MaxQuant, MSFragger, MODa, MODPlus, MSAlign, MSPathFinder,
// TopPIC, and X!Tandem, along with SEQUEST Synopsis/First Hits files.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://www.pnnl.gov/integrative-omics
// -------------------------------------------------------------------------------

using System;
using System.IO;
using System.Reflection;
using System.Threading;
using PeptideHitResultsProcessor;
using PeptideHitResultsProcessor.Processor;
using PRISM;

namespace PeptideHitResultsProcRunner
{
    internal static class Program
    {
        // Ignore Spelling: conf, Prot, MaxQuant, MODa, txt

        private static DateTime mLastProgressReportTime;
        private static int mLastProgressReportValue;
        private static DateTime mLastProgressReportValueTime;
        private static bool mSkippedInitialProgressValue;

        /// <summary>
        /// Program entry point
        /// </summary>
        /// <returns>0 if no error, error code if an error</returns>
        public static int Main(string[] args)
        {
            var exeName = Assembly.GetEntryAssembly()?.GetName().Name;

            var parser = new CommandLineParser<PHRPOptions>(exeName, PHRPBaseClass.GetAppVersion())
            {
                ProgramInfo = "This program converts search results from various MS/MS identification tools " +
                              "into a series of tab-delimited text files that organize the data in a similar format for each tool. " +
                              "It supports MS-GF+, MaxQuant, MSFragger, MODa, MODPlus, MSAlign, " +
                              "MSPathFinder, TopPIC, and X!Tandem, along with SEQUEST Synopsis/First Hits files.",
                ContactInfo = "Program written by Matthew Monroe for PNNL (Richland, WA)" + Environment.NewLine +
                              "E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov" + Environment.NewLine +
                              "Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://www.pnnl.gov/integrative-omics"
            };

            parser.UsageExamples.Add(string.Format(
                "The input file should be one of the following:\n" +
                "  MaxQuant results files (msms.txt and peptides.txt)\n" +
                "  MS-GF+ results file ({0}.tsv or {1}.tsv or .tsv)\n" +
                "  MSGF-DB results file ({1}.txt)\n" +
                "  MSAlign results file ({2}.txt)\n" +
                "  MODa results file ({3}.txt)\n" +
                "  MODPlus results file ({4}.txt)\n" +
                "  MsFragger results file (_psm.tsv)\n" +
                "  MSPathFinder results file ({5}.txt)\n" +
                "  InSpecT results file ({6}.txt)\n" +
                "  SEQUEST Synopsis File ({7}.txt)\n" +
                "  SEQUEST First Hits file ({8}.txt)\n" +
                "  TopPIC results file ({9}.txt)\n" +
                "  X!Tandem Results file (_xt.xml)",
                MSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE,
                MSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE,
                MSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE,
                MODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE,
                MODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE,
                MSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE,
                InSpecTResultsProcessor.FILENAME_SUFFIX_INSPECT_FILE,
                SequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE,
                SequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE,
                TopPICResultsProcessor.FILENAME_SUFFIX_TopPIC_PRSMs_FILE
                ));

            // The default argument name for parameter files is /ParamFile or -ParamFile
            // Also allow /Conf or /P
            parser.AddParamFileKey("Conf");
            parser.AddParamFileKey("P");

            var result = parser.ParseArgs(args);
            var options = result.ParsedResults;

            if (!result.Success || !options.Validate())
            {
                if (parser.CreateParamFileProvided)
                {
                    return 0;
                }

                // Delay for 750 msec in case the user double-clicked this file from within Windows Explorer (or started the program via a shortcut)
                Thread.Sleep(750);

                return -1;
            }

            try
            {
                if (InvalidParameterFile(parser.ParameterFilePath))
                {
                    return -1;
                }

                var peptideHitResultsProcessor = new PeptideHitResultsProcRunner(options);

                peptideHitResultsProcessor.DebugEvent += PeptideHitResultsProcRunner_DebugEvent;
                peptideHitResultsProcessor.ErrorEvent += PeptideHitResultsProcRunner_ErrorEvent;
                peptideHitResultsProcessor.StatusEvent += PeptideHitResultsProcRunner_MessageEvent;
                peptideHitResultsProcessor.ProgressUpdate += PeptideHitResultsProcRunner_ProgressChanged;
                peptideHitResultsProcessor.ProgressReset += PeptideHitResultsProcRunner_ProgressReset;
                peptideHitResultsProcessor.WarningEvent += PeptideHitResultsProcRunner_WarningEvent;

                if (options.LogMessagesToFile && args.Length > 0)
                {
                    // Append the command line arguments to the log file
                    peptideHitResultsProcessor.SkipConsoleWriteIfNoStatusListener = true;
                    peptideHitResultsProcessor.LogAdditionalMessage(string.Join(" ", args));
                    peptideHitResultsProcessor.SkipConsoleWriteIfNoStatusListener = false;
                }

                int returnCode;
                if (options.RecurseDirectories)
                {
                    Console.WriteLine("Recursively processing files in the input file's directory and below");

                    if (peptideHitResultsProcessor.ProcessFilesAndRecurseDirectories(
                        options.InputFilePath,
                        options.OutputDirectoryPath,
                        string.Empty,           // options.OutputDirectoryAlternatePath
                        false,        // options.RecreateDirectoryHierarchyInAlternatePath
                        options.XmlParameterFile,
                        options.MaxLevelsToRecurse))
                    {
                        returnCode = 0;
                    }
                    else
                    {
                        returnCode = (int)peptideHitResultsProcessor.ErrorCode;
                    }
                }
                else
                {
                    bool success;

                    if (options.InputFilePath.Contains("*") || options.InputFilePath.Contains("?"))
                    {
                        success = peptideHitResultsProcessor.ProcessFilesWildcard(
                            options.InputFilePath,
                            options.OutputDirectoryPath,
                            options.XmlParameterFile);
                    }
                    else
                    {
                        success = peptideHitResultsProcessor.ProcessFile(
                            options.InputFilePath,
                            options.OutputDirectoryPath,
                            options.XmlParameterFile);
                    }

                    if (success)
                    {
                        returnCode = 0;
                    }
                    else
                    {
                        var errorCode = (int)peptideHitResultsProcessor.ErrorCode;

                        if (errorCode == 0)
                        {
                            returnCode = -1;
                            ShowErrorMessage("ProcessFilesWildcard returned Success=False");
                        }
                        else
                        {
                            returnCode = errorCode;
                            ShowErrorMessage("Error while processing: " + peptideHitResultsProcessor.GetErrorMessage());
                        }
                    }
                }

                DisplayProgressPercent(mLastProgressReportValue, true);
                return returnCode;
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in modMain->Main: " + Environment.NewLine + ex.Message);
                return -1;
            }
        }

        private static void DisplayProgressPercent(int percentComplete, bool addCarriageReturn)
        {
            if (addCarriageReturn)
            {
                Console.WriteLine();
            }

            if (percentComplete > 100)
            {
                percentComplete = 100;
            }

            Console.Write("Processing: " + percentComplete + "% ");

            if (addCarriageReturn)
            {
                Console.WriteLine();
            }
        }

        private static bool InvalidParameterFile(string parameterFilePath)
        {
            if (string.IsNullOrWhiteSpace(parameterFilePath))
            {
                return false;
            }

            // Assure that the user did not provide an old-style XML-based parameter file
            var paramFile = new FileInfo(parameterFilePath);

            if (!paramFile.Extension.Equals(".xml", StringComparison.OrdinalIgnoreCase))
            {
                return false;
            }

            using var paramFileReader = new StreamReader(new FileStream(paramFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

            var linesRead = 0;

            while (!paramFileReader.EndOfStream && linesRead < 5)
            {
                var dataLine = paramFileReader.ReadLine();

                if (string.IsNullOrWhiteSpace(dataLine))
                {
                    continue;
                }

                linesRead++;

                var trimmedLine = dataLine.Trim();

                if (trimmedLine.StartsWith("<?xml", StringComparison.OrdinalIgnoreCase) ||
                    trimmedLine.StartsWith("<sections", StringComparison.OrdinalIgnoreCase))
                {
                    ConsoleMsgUtils.ShowWarning(
                        "PeptideHitResultsProcRunner v3.1 uses Key=Value parameter files\n" +
                        "{0} is an XML file\n" +
                        "For an example parameter file, see file PHRP_Options.conf at {1}",
                        paramFile.Name, "https://github.com/PNNL-Comp-Mass-Spec/PHRP/tree/master/Data");
                    ConsoleMsgUtils.ShowWarning("Aborting");
                    return true;
                }
            }

            return false;
        }

        private static void ShowErrorMessage(string message, Exception ex = null)
        {
            ConsoleMsgUtils.ShowError(message, ex);
        }

        private static void PeptideHitResultsProcRunner_DebugEvent(string message)
        {
            ConsoleMsgUtils.ShowDebug(message);
        }

        private static void PeptideHitResultsProcRunner_ErrorEvent(string message, Exception ex)
        {
            ShowErrorMessage(message, ex);
        }

        private static void PeptideHitResultsProcRunner_MessageEvent(string message)
        {
            Console.WriteLine(message);
        }

        private static void PeptideHitResultsProcRunner_ProgressChanged(string taskDescription, float percentComplete)
        {
            const int PERCENT_REPORT_INTERVAL = 25;
            const int PROGRESS_DOT_INTERVAL_MSEC = 250;
            const int PROGRESS_VALUE_INTERVAL_SEC = 60;

            if (percentComplete >= mLastProgressReportValue)
            {
                if (!mSkippedInitialProgressValue)
                {
                    mSkippedInitialProgressValue = true;
                    return;
                }

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

        private static void PeptideHitResultsProcRunner_ProgressReset()
        {
            mLastProgressReportTime = DateTime.UtcNow;
            mLastProgressReportValueTime = DateTime.UtcNow;
            mLastProgressReportValue = 0;
        }

        private static void PeptideHitResultsProcRunner_WarningEvent(string message)
        {
            if (message.StartsWith("Warning", StringComparison.OrdinalIgnoreCase))
            {
                ConsoleMsgUtils.ShowWarning(message);
            }
            else
            {
                ConsoleMsgUtils.ShowWarning("Warning: " + message);
            }
        }
    }
}
