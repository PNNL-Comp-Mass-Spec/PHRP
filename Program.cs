// This program converts search results from various MS/MS identification tools
// into a series of tab-delimited text files summarizing the results.
// It supports MS-GF+, MaxQuant, MSFragger, MODa, MODPlus, MSAlign,
// MSPathFinder, TopPIC, and X!Tandem, along with SEQUEST Synopsis/First Hits files.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
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

            var cmdLineParser = new CommandLineParser<PHRPOptions>(exeName, PHRPBaseClass.GetAppVersion())
            {
                ProgramInfo = "This program converts search results from various MS/MS identification tools " +
                              "into a series of tab-delimited text files summarizing the results.  " +
                              "It supports MS-GF+, MaxQuant, MSFragger, MODa, MODPlus, MSAlign,  " +
                              "MSPathFinder, TopPIC, and X!Tandem, along with SEQUEST Synopsis/First Hits files.",
                ContactInfo = "Program written by Matthew Monroe for PNNL (Richland, WA)" + Environment.NewLine +
                              "E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov" + Environment.NewLine +
                              "Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics"
            };

            cmdLineParser.UsageExamples.Add(string.Format(
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
            cmdLineParser.AddParamFileKey("Conf");
            cmdLineParser.AddParamFileKey("P");

            var result = cmdLineParser.ParseArgs(args);
            var options = result.ParsedResults;
            if (!result.Success || !options.Validate())
            {
                // Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
                Thread.Sleep(750);
                return -1;
            }

            try
            {
                if (InvalidParameterFile(cmdLineParser.ParameterFilePath))
                    return -1;

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
                        options.ParameterFilePath,
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
                    if (peptideHitResultsProcessor.ProcessFilesWildcard(
                        options.InputFilePath,
                        options.OutputDirectoryPath,
                        options.ParameterFilePath))
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
                Console.WriteLine();

            if (percentComplete > 100)
                percentComplete = 100;

            Console.Write("Processing: " + percentComplete + "% ");

            if (addCarriageReturn)
                Console.WriteLine();
        }

        private static bool InvalidParameterFile(string parameterFilePath)
        {
            if (string.IsNullOrWhiteSpace(parameterFilePath))
                return false;

            // Assure that the user did not provide an old-style XML-based parameter file
            var paramFile = new FileInfo(parameterFilePath);
            if (!paramFile.Extension.Equals(".xml", StringComparison.OrdinalIgnoreCase))
                return false;

            using var paramFileReader = new StreamReader(new FileStream(paramFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

            var linesRead = 0;
            while (!paramFileReader.EndOfStream && linesRead < 5)
            {
                var dataLine = paramFileReader.ReadLine();
                if (string.IsNullOrWhiteSpace(dataLine))
                    continue;

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

        private static void ShowProgramHelp()
        {
            try
            {
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "This program converts search results from various MS/MS identification tools into a series of tab-delimited text files summarizing the results. " +
                                      "It supports MS-GF+, MaxQuant, MSFragger, MODa, MODPlus, MSAlign, MSPathFinder, " +
                                      "TopPIC, and X!Tandem, along with SEQUEST Synopsis/First Hits files."));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "It will insert modification symbols into the peptide sequences for modified peptides. " +
                                      "Parallel files will be created containing sequence info and modification details."));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "The user can optionally provide a modification definition file " +
                                      "that specifies the symbol to use for each modification mass."));
                Console.WriteLine();
                Console.WriteLine("Program syntax:" + Environment.NewLine +
                                  Path.GetFileName(Assembly.GetExecutingAssembly().Location) +
                                  " InputFilePath [/O:OutputDirectoryPath]");
                // Future:
                // Console.WriteLine(" [/DatasetDir:DatasetDirectoryPath]")
                //
                Console.WriteLine(" [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath]");
                Console.WriteLine(" [/ProteinMods] [/F:FastaFilePath]");
                Console.WriteLine(" [/ProteinModsViaPHRP] [/IgnorePepToProtMapErrors]");
                Console.WriteLine(" [/ProteinModsIncludeReversed] [/UseExistingPepToProteinMapFile]");
                Console.WriteLine(" [/T:MassCorrectionTagsFilePath] [/N:SearchToolParameterFilePath]");
                Console.WriteLine(" [/MSGFPlusSpecEValue:0.0000005] [/MSGFPlusEValue:0.75]");
                // Console.WriteLine(" [/InsSynPValue:0.2]");
                Console.WriteLine(" [/FHT:True|False] [/Syn:True|False]");
                Console.WriteLine(" [/SynProb:0.05] [/SynPValue:0.95]");
                Console.WriteLine(" [/MaxQScore:50] [/MaxQPEP:0.01]");
                Console.WriteLine(" [/DB:DatabaseConnectionString]");
                Console.WriteLine(" [/S:[MaxLevel]] [/A:AlternateOutputDirectoryPath] [/R]");
                Console.WriteLine(" [/L:[LogFilePath]] [/LogDir:LogDirectoryPath]");
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(string.Format(
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
                                      )));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "The output directory switch is optional. " +
                                      "If omitted, the output file will be created " +
                                      "in the same directory as the input file."));
                Console.WriteLine();
                // Future:
                // Console.WriteLine("As an alternative to specifying an input file, you can specify an input directory. " +
                //                   "In this case the program will look for the best file to process from that directory, and will auto-determine /T and /N")
                //

                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "The parameter file path is optional. If included, it should point to a valid XML parameter file."));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /M to specify the file containing the modification definitions. " +
                                      "This file should be tab delimited, with the first column containing the modification symbol, " +
                                      "the second column containing the modification mass, plus optionally a third column " +
                                      "listing the residues that can be modified with the given mass " +
                                      "(1 letter residue symbols, no need to separate with commas or spaces)."));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /ProteinMods to indicate that the _ProteinMods.txt file should be created. " +
                                      "This requires that either an existing _PepToProtMapMTS.txt file exist, " +
                                      "or that the FASTA file be defined using /F"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /ProteinModsViaPHRP to indicate that InputFilePath specifies a valid PHRP data file " +
                                      "and thus the PHRP data files should not be re-created; only the _ProteinMods.txt file " +
                                      "should be created. This requires that either an existing _PepToProtMapMTS.txt file exist, " +
                                      "or that the FASTA file be defined using /F"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /F to specify the path to the FASTA file. When provided, the order of the proteins " +
                                      "in the FASTA file dictates which protein is listed for each peptide in the First Hits file"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /IgnorePepToProtMapErrors to ignore peptide to protein mapping errors " +
                                      "that occur when creating a missing _PepToProtMapMTS.txt file"));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /ProteinModsIncludeReversed to include Reversed proteins in the _ProteinMods.txt file"));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /UseExistingPepToProteinMapFile to use an existing _PepToProtMapMTS.txt file if it exists"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /T to specify the file containing the mass correction tag info. This file " +
                                      "should be tab delimited, with the first column containing the mass correction tag name " +
                                      "and the second column containing the mass (the name cannot contain commas or colons)."));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /N to specify the parameter file provided to the search tool. " +
                                      "This is used when processing results from MS-GF+, MSPathFinder, MaxQuant, MODa, MODPlus, MSAlign, TopPIC, and InSpecT. " +
                                      "For MaxQuant, provide either an XML-based parameter file (root element is <MaxQuantParams>) " +
                                      "or provide the parameters.txt file created in the txt results directory. " +
                                      "The XML-based parameter file is preferred, since it is required to allow PHRP " +
                                      "to accurately compute monoisotopic masses of peptides identified by MaxQuant."));
                Console.WriteLine();

                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "When processing an MS-GF+ results file, use /MSGFPlusSpecEValue and /MSGFPlusEValue " +
                                      "to customize the thresholds used to determine which peptides are written to the synopsis file"));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Defaults are /MSGFPlusSpecEValue:" +
                                      MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD +
                                      " and /MSGFPlusEValue:" +
                                      MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD));
                //Console.WriteLine();
                //Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                //                      "When processing an InSpecT results file, use /InsSynPValue to customize " +
                //                      "the PValue threshold used to determine which peptides are written to the synopsis file. " +
                //                      "The default is /InsSynPValue:0.2  Note that peptides with a " +
                //                      "TotalPRMScore >= " + InSpecTResultsProcessor.TOTALPRMSCORE_THRESHOLD +
                //                      " or an FScore >= " + InSpecTResultsProcessor.FSCORE_THRESHOLD +
                //                      " will also be included in the synopsis file."));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /FHT:True or /FHT:False to control the creation of a first-hits file (_fht.txt) " +
                                      "when processing results from MS-GF+, MaxQuant, MSFragger, etc. (default is /FHT:True)"));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /Syn:True or /Syn:False to control the creation of a synopsis file (_syn.txt) " +
                                      "when processing results from MS-GF+, MaxQuant, MSFragger, etc. (default is /Syn:True)"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "When processing a MODPlus or MODa results file, use /SynProb to customize " +
                                      "the probability threshold used to determine which peptides are written to the synopsis file. " +
                                      "The default is /SynProb:0.05"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                    "When processing a MODPlus or MODa results file, use /SynPValue to customize " +
                    "the p-value threshold used to determine which peptides are written to the synopsis file. " +
                    "The default is /SynPValue:0.95"));

                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                    "When processing a MaxQuant results file, use /MaxQScore to customize " +
                    "the Andromeda score threshold used to determine which peptides are written to the synopsis file. " +
                    "A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold. " +
                    "The default is /MaxQScore:50"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                    "When processing a MaxQuant results file, use /MaxQPEP to customize " +
                    "the Posterior Error Probability (PEP) score threshold used to determine which peptides are written to the synopsis file. " +
                    "The default is /MaxQPEP:0.01"));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(string.Format(
                    "When processing MaxQuant results using a computer on the pnl.gov domain, " +
                    "the DMS database is contacted to lookup dataset IDs by dataset name, " +
                    "where dataset name comes from the 'Raw file' column in the msms.txt file. " +
                    "The default is \n/DB:\"{0}\"", PHRPOptions.DMS_CONNECTION_STRING)));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /S to process all valid files in the input directory and subdirectories. " +
                                      "Include a number after /S (like /S:2) to limit the level of subdirectories to examine."));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "When using /S, you can redirect the output of the results using /A."));
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "When using /S, you can use /R to re-create the input directory hierarchy " +
                                      "in the alternate output directory (if defined)."));
                Console.WriteLine();
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                                      "Use /L to specify that a log file should be created. " +
                                      "Use /L:LogFilePath to specify the name (or full path) for the log file. " +
                                      "Use /LogDir to specify the directory to create the log file (ignored if the LogFilePath is rooted)"));
                Console.WriteLine();

                Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006");
                Console.WriteLine("Version: " + PHRPBaseClass.GetAppVersion());

                Console.WriteLine();

                Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov");
                Console.WriteLine("Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics");
                Console.WriteLine();

                // Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
                Thread.Sleep(750);
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error displaying the program syntax: " + ex.Message);
            }
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
                ConsoleMsgUtils.ShowWarning(message);
            else
                ConsoleMsgUtils.ShowWarning("Warning: " + message);
        }
    }
}
