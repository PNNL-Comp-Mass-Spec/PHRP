﻿// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://www.pnnl.gov/integrative-omics
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.IO;
using PeptideHitResultsProcessor;
using PeptideHitResultsProcessor.Processor;
using PHRPReader;
using PRISM;

namespace PeptideHitResultsProcRunner
{
    /// <summary>
    /// This class calls SequestSynopsisFileProcessor or XTandemResultsConverter
    /// to process the files to determine the modifications present for each peptide,
    /// along with other information
    /// </summary>
    public class PeptideHitResultsProcRunner : PRISM.FileProcessor.ProcessFilesBase
    {
        // ReSharper disable once CommentTypo
        // Ignore Spelling: Battelle, enums, Hyperscore, MaxQuant, MODa, parm, proc

        /// <summary>
        /// Constructor
        /// </summary>
        public PeptideHitResultsProcRunner(PHRPOptions options)
        {
            mFileDate = PHRPBaseClass.PROGRAM_DATE;
            Options = options;

            UpdateLogFileOptions();

            mPeptideHitResultsFileFormat = ResultsFileFormat.AutoDetermine;
            mUseExistingMTSPepToProteinMapFile = false;
            mLocalErrorCode = ResultsProcessorErrorCodes.NoError;
        }

        /// <summary>
        /// Error codes specialized for this class
        /// </summary>
        public enum ResultsProcessorErrorCodes
        {
            NoError = 0,
            UnspecifiedError = -1
        }

        private bool mFilePathsShown;

        private bool mLogFileOptionsUpdated;

        protected ResultsFileFormat mPeptideHitResultsFileFormat;

        private PHRPBaseClass mPeptideHitResultsProcessor;

        protected ResultsProcessorErrorCodes mLocalErrorCode;

        private bool mOptionsDisplayed;

        protected bool mUseExistingMTSPepToProteinMapFile;

        /// <summary>
        /// Processing options
        /// </summary>
        public PHRPOptions Options { get; }

        // ReSharper disable once UnusedMember.Global
        public ResultsProcessorErrorCodes LocalErrorCode => mLocalErrorCode;

        /// <summary>
        /// If filePath is an empty string, return textIfEmptyPath
        /// Otherwise, return the file path, truncated if over the specified maximum length
        /// </summary>
        /// <param name="filePath">File path</param>
        /// <param name="textIfEmptyPath">Text to return if filePath is an empty string</param>
        /// <param name="maxPathLength">Maximum path length</param>
        private string FilePathOrText(string filePath, string textIfEmptyPath, int maxPathLength = 110)
        {
            if (string.IsNullOrWhiteSpace(filePath))
                return textIfEmptyPath;

            return PathUtils.CompactPathString(filePath, maxPathLength);
        }

        /// <summary>
        /// Get the default file extensions to parse
        /// </summary>
        public override IList<string> GetDefaultExtensionsToParse()
        {
            return new List<string>
            {
                ".txt",
                ".xml",
                ".tsv"
            };
        }

        /// <summary>
        /// Get the error message, or an empty string if no error
        /// </summary>
        public override string GetErrorMessage()
        {
            if (ErrorCode is ProcessFilesErrorCodes.LocalizedError or ProcessFilesErrorCodes.NoError)
            {
                return mLocalErrorCode switch
                {
                    ResultsProcessorErrorCodes.NoError => string.Empty,
                    ResultsProcessorErrorCodes.UnspecifiedError => "Unspecified localized error",
                    _ => "Unknown error state"             // This shouldn't happen
                };
            }

            return GetBaseClassErrorMessage();
        }

        private void InitializePeptideHitResultsProcessor(string inputFilePath)
        {
            var sourceFile = new FileInfo(PathUtils.GetCleanPath(inputFilePath));

            mPeptideHitResultsProcessor.Options.MassCorrectionTagsFilePath = PHRPBaseClass.ResolveFilePath(sourceFile.DirectoryName, Options.MassCorrectionTagsFilePath);
            mPeptideHitResultsProcessor.Options.ModificationDefinitionsFilePath = PHRPBaseClass.ResolveFilePath(sourceFile.DirectoryName, Options.ModificationDefinitionsFilePath);
            mPeptideHitResultsProcessor.Options.SearchToolParameterFilePath = PHRPBaseClass.ResolveFilePath(sourceFile.DirectoryName, Options.SearchToolParameterFilePath);
        }

        /// <summary>
        /// Log an additional message to the log file
        /// </summary>
        /// <param name="message">Message</param>
        public void LogAdditionalMessage(string message)
        {
            LogMessage(message);
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">PSM tool results file</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="parameterFilePath">Parameter file path</param>
        /// <param name="resetErrorCode">When true, reset the error code</param>
        /// <returns>True if successful, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath, bool resetErrorCode)
        {
            if (!mLogFileOptionsUpdated)
            {
                UpdateLogFileOptions();
            }

            if (resetErrorCode)
            {
                SetLocalErrorCode(ResultsProcessorErrorCodes.NoError);
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    ShowErrorMessage("Input file name is empty");
                    SetBaseClassErrorCode(ProcessFilesErrorCodes.InvalidInputFilePath);
                    return false;
                }

                if (PHRPBaseClass.FindInputFile(inputFilePath, Options.AlternateBasePath, out var inputFileToUse))
                {
                    inputFilePath = inputFileToUse.FullName;
                }

                // Note that CleanupFilePaths() will update mOutputDirectoryPath, which is used by LogMessage()
                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                {
                    SetBaseClassErrorCode(ProcessFilesErrorCodes.FilePathError);
                    if (inputFilePath?.Contains("..") == true)
                    {
                        var inputFile = new FileInfo(inputFilePath);
                        OnStatusEvent("Absolute path: " + inputFile.DirectoryName);
                    }

                    return false;
                }

                UpdateProgress("Parsing " + Path.GetFileName(inputFilePath));
                ResetProgress();

                // ReSharper disable once ConvertIfStatementToReturnStatement
                if (Options.CreateProteinModsUsingPHRPDataFile)
                {
                    return StartCreateProteinModsViaPHRPData(inputFilePath, outputDirectoryPath);
                }

                return StartPHRP(inputFilePath, outputDirectoryPath, parameterFilePath);
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFile", ex);
                OnDebugEvent(StackTraceFormatter.GetExceptionStackTraceMultiLine(ex));
                return false;
            }
        }

        private void RegisterResultsProcessEvents(PHRPBaseClass resultsProcessor)
        {
            resultsProcessor.ErrorEvent += PeptideHitResultsProcessor_ErrorOccurred;
            resultsProcessor.StatusEvent += PeptideHitResultsProcessor_MessageEvent;
            resultsProcessor.ProgressUpdate += PeptideHitResultsProcessor_ProgressChanged;
            resultsProcessor.WarningEvent += PeptideHitResultsProcessor_WarningMessageEvent;
            resultsProcessor.ProgressReset += PeptideHitResultsProcessor_ProgressReset;
        }

        private void SetLocalErrorCode(ResultsProcessorErrorCodes newErrorCode, bool leaveExistingErrorCodeUnchanged = false)
        {
            if (leaveExistingErrorCodeUnchanged && mLocalErrorCode != ResultsProcessorErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mLocalErrorCode = newErrorCode;

                if (newErrorCode == ResultsProcessorErrorCodes.NoError)
                {
                    if (ErrorCode == ProcessFilesErrorCodes.LocalizedError)
                    {
                        SetBaseClassErrorCode(ProcessFilesErrorCodes.NoError);
                    }
                }
                else
                {
                    SetBaseClassErrorCode(ProcessFilesErrorCodes.LocalizedError);
                }
            }
        }

        private void ShowCurrentFilePaths(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            Console.WriteLine();

            LogMessage(string.Format("{0,-45} {1}",
                "Input file:", PathUtils.CompactPathString(inputFilePath, 110)));

            if (string.IsNullOrWhiteSpace(outputDirectoryPath))
            {
                LogMessage(string.Format("{0,-45} {1}",
                    "Output directory:", "Same directory as the input file"));
            }
            else
            {
                LogMessage(string.Format("{0,-45} {1}",
                    "Output directory:", PathUtils.CompactPathString(outputDirectoryPath, 110)));
            }

            if (!string.IsNullOrWhiteSpace(Options.OutputFileBaseName))
            {
                LogMessage(string.Format("{0,-45} {1}",
                    "Output file base name:", Options.OutputFileBaseName));
            }

            if (!string.IsNullOrWhiteSpace(parameterFilePath))
            {
                LogMessage(string.Format("{0,-45} {1}",
                    "Parameter file:", PathUtils.CompactPathString(parameterFilePath, 110)));
            }
        }

        private void ShowProcessingOptions(PHRPBaseClass resultsProcessor)
        {
            Console.WriteLine();

            // Processing options for MSGFPlus results processor
            LogMessage("Processing options for " + resultsProcessor);

            if (!string.IsNullOrWhiteSpace(resultsProcessor.Options.XmlParameterFile))
            {
                LogMessage(string.Format("{0,-45} {1}",
                    "Legacy XML parameter file:", PathUtils.CompactPathString(resultsProcessor.Options.XmlParameterFile, 110)));
            }

            LogMessage(string.Format("{0,-45} {1}",
                "Search Tool Parameter File:", FilePathOrText(resultsProcessor.Options.SearchToolParameterFilePath, "Not defined")));

            LogMessage(string.Format("{0,-45} {1}",
                "Modification Definitions File:", FilePathOrText(resultsProcessor.Options.ModificationDefinitionsFilePath, "Not defined")));

            LogMessage(string.Format("{0,-45} {1}",
                "Mass Correction Tags File:", FilePathOrText(resultsProcessor.Options.MassCorrectionTagsFilePath, "Use internally defined tags")));

            LogMessage(string.Format("{0,-45} {1}",
                "FASTA File:", FilePathOrText(resultsProcessor.Options.FastaFilePath, "Not defined")));

            Console.WriteLine();
            LogMessage(string.Format("{0,-45} {1}",
                "Create Protein Mods File:", resultsProcessor.Options.CreateProteinModsFile));

            LogMessage(string.Format("{0,-45} {1}",
                "Ignore Peptide to Protein Mapper Errors:", resultsProcessor.Options.IgnorePeptideToProteinMapperErrors));

            LogMessage(string.Format("{0,-45} {1}",
                "Protein Mods File Includes Reversed Proteins:", resultsProcessor.Options.ProteinModsFileIncludesReversedProteins));

            LogMessage(string.Format("{0,-45} {1}",
                "Use Existing MTS PepToProtein Map File:", resultsProcessor.Options.UseExistingMTSPepToProteinMapFile));

            Console.WriteLine();
            LogMessage(string.Format("{0,-45} {1}",
                "Create First Hits File:", resultsProcessor.Options.CreateFirstHitsFile));

            LogMessage(string.Format("{0,-45} {1}",
                "Create Synopsis File:", resultsProcessor.Options.CreateSynopsisFile));

            if (resultsProcessor is InSpecTResultsProcessor)
            {
                LogMessage(string.Format("{0,-45} {1:E2}",
                    "InSpecT Synopsis File PValue Threshold:", resultsProcessor.Options.InspectSynopsisFilePValueThreshold));
            }

            Console.WriteLine();

            if (resultsProcessor is MODaResultsProcessor or MODPlusResultsProcessor)
            {
                LogMessage(string.Format("{0,-49} {1:F2}",
                    "MODa/MODPlus Synopsis File Probability Threshold:", resultsProcessor.Options.MODaMODPlusSynopsisFileProbabilityThreshold));
            }

            if (resultsProcessor is MSGFPlusResultsProcessor)
            {
                LogMessage(string.Format("{0,-45} {1:F2}",
                    "MSGFPlus Synopsis File EValue Threshold:", resultsProcessor.Options.MSGFPlusSynopsisFileEValueThreshold));
            }

            if (resultsProcessor is MSGFPlusResultsProcessor or MSPathFinderResultsProcessor)
            {
                LogMessage(string.Format("{0,-45} {1:E2}",
                    "MSGFPlus Synopsis File SpecEValue Threshold:", resultsProcessor.Options.MSGFPlusSynopsisFileSpecEValueThreshold));
            }

            if (resultsProcessor is DiaNNResultsProcessor)
            {
                LogMessage(string.Format("{0,-49} {1:F3}",
                    "DIA-NN Synopsis File Q-Value Threshold:", resultsProcessor.Options.DiaNNQValueThreshold));

                LogMessage(string.Format("{0,-49} {1:F3}",
                    "DIA-NN Synopsis File Confidence Score Threshold:", resultsProcessor.Options.DiaNNConfidenceScoreThreshold));
            }

            if (resultsProcessor is MaxQuantResultsProcessor)
            {
                LogMessage(string.Format("{0,-49} {1:N0}",
                    "MaxQuant Synopsis File Andromeda Score Threshold:", resultsProcessor.Options.MaxQuantAndromedaScoreThreshold));

                LogMessage(string.Format("{0,-49} {1:F3}",
                    "MaxQuant Synopsis File PEP Threshold:", resultsProcessor.Options.MaxQuantPosteriorErrorProbabilityThreshold));
            }

            if (resultsProcessor is MSFraggerResultsProcessor)
            {
                LogMessage(string.Format("{0,-45} {1:F3}",
                    "MSFragger Synopsis File E-Value Threshold:", resultsProcessor.Options.MSGFPlusSynopsisFileEValueThreshold));

                LogMessage(string.Format("{0,-45} {1:N0}",
                    "MaxQuant Synopsis File Hyperscore Threshold:", resultsProcessor.Options.MSFraggerHyperscoreThreshold));
            }

            Console.WriteLine();
        }

        private bool StartCreateProteinModsViaPHRPData(string inputFilePath, string outputDirectoryPath)
        {
            try
            {
                if (!mFilePathsShown)
                {
                    ShowCurrentFilePaths(inputFilePath, outputDirectoryPath, string.Empty);
                    mFilePathsShown = true;
                }

                PeptideHitResultTypes PeptideHitResultType;
                switch (mPeptideHitResultsFileFormat)
                {
                    case ResultsFileFormat.SequestFirstHitsFile:
                        PeptideHitResultType = PeptideHitResultTypes.Sequest;
                        LogMessage("Detected SEQUEST First Hits file");
                        break;

                    case ResultsFileFormat.SequestSynopsisFile:
                        PeptideHitResultType = PeptideHitResultTypes.Sequest;
                        LogMessage("Detected SEQUEST Synopsis file");
                        break;

                    case ResultsFileFormat.XTandemXMLFile:
                        PeptideHitResultType = PeptideHitResultTypes.XTandem;
                        LogMessage("Detected X!Tandem XML file");
                        break;

                    case ResultsFileFormat.InspectTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.Inspect;
                        LogMessage("Detected InSpecT results file");
                        break;

                    case ResultsFileFormat.MSGFPlusTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.MSGFPlus;
                        LogMessage("Detected MS-GF+ results file");
                        break;

                    case ResultsFileFormat.MSAlignTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.MSAlign;
                        LogMessage("Detected MSAlign results file");
                        break;

                    case ResultsFileFormat.MODPlusTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.MODPlus;
                        LogMessage("Detected MODPlus results file");
                        break;

                    case ResultsFileFormat.MSPathFinderTSVFile:
                        PeptideHitResultType = PeptideHitResultTypes.MSPathFinder;
                        LogMessage("Detected MSPathfinder results file");
                        break;

                    case ResultsFileFormat.TopPICTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.TopPIC;
                        LogMessage("Detected TopPIC results file");
                        break;

                    case ResultsFileFormat.MaxQuantTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.MaxQuant;
                        LogMessage("Detected MaxQuant results file");
                        break;

                    default:
                        // Includes PeptideHitResultsFileFormatConstants.AutoDetermine
                        PeptideHitResultType = ReaderFactory.AutoDetermineResultType(inputFilePath);
                        break;
                }

                if (PeptideHitResultType == PeptideHitResultTypes.Unknown)
                {
                    ShowErrorMessage("Error: Could not determine the format of the PHRP data file: " + inputFilePath);
                    return false;
                }

                switch (PeptideHitResultType)
                {
                    case PeptideHitResultTypes.Sequest:
                        mPeptideHitResultsProcessor = new SequestResultsProcessor(Options);
                        break;

                    case PeptideHitResultTypes.XTandem:
                        mPeptideHitResultsProcessor = new XTandemResultsProcessor(Options);
                        break;

                    case PeptideHitResultTypes.Inspect:
                        mPeptideHitResultsProcessor = new InSpecTResultsProcessor(Options);
                        break;

                    case PeptideHitResultTypes.MSGFPlus:
                        mPeptideHitResultsProcessor = new MSGFPlusResultsProcessor(Options);
                        break;

                    case PeptideHitResultTypes.MSAlign:
                        mPeptideHitResultsProcessor = new MSAlignResultsProcessor(Options);
                        break;

                    case PeptideHitResultTypes.MODa:
                        mPeptideHitResultsProcessor = new MODaResultsProcessor(Options);
                        break;

                    default:
                        // Unknown format
                        ShowErrorMessage("Error: Unrecognized value for PeptideHitResultType: " + PeptideHitResultType);
                        return false;
                }

                if (mPeptideHitResultsProcessor == null)
                {
                    return false;
                }

                // Do not call RegisterEvents
                // Instead use local event handlers that optionally log to a file
                RegisterResultsProcessEvents(mPeptideHitResultsProcessor);

                InitializePeptideHitResultsProcessor(inputFilePath);

                if (!mOptionsDisplayed)
                {
                    ShowProcessingOptions(mPeptideHitResultsProcessor);
                    mOptionsDisplayed = true;
                }

                var success = mPeptideHitResultsProcessor.CreateProteinModDetailsFile(inputFilePath, outputDirectoryPath);

                if (!success)
                {
                    ShowErrorMessage(mPeptideHitResultsProcessor.ErrorMessage);
                    return false;
                }

                LogMessage("Processing Complete");
                OperationComplete();
                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error calling CreateProteinModDetailsFile in StartCreateProteinModsViaPHRPData", ex);
                return false;
            }
        }

        private bool StartPHRP(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            try
            {
                if (!mFilePathsShown)
                {
                    ShowCurrentFilePaths(inputFilePath, outputDirectoryPath, string.Empty);
                    mFilePathsShown = true;
                }

                ResultsFileFormat peptideHitResultsFormat;
                bool inputFileIsSynOrFht;

                if (mPeptideHitResultsFileFormat == ResultsFileFormat.AutoDetermine)
                {
                    inputFileIsSynOrFht = ReaderFactory.IsSynopsisOrFirstHitsFile(inputFilePath);
                    if (inputFileIsSynOrFht)
                    {
                        var fileType = inputFilePath.EndsWith("fht.txt", StringComparison.OrdinalIgnoreCase)
                            ? "first hits"
                            : "synopsis";

                        ShowWarning(string.Format(
                            "Invalid file format for PHRP; the input file is already a {0} file: {1}",
                            fileType, inputFilePath));

                        peptideHitResultsFormat = ResultsFileFormat.AutoDetermine;
                    }
                    else
                    {
                        peptideHitResultsFormat = PHRPBaseClass.DetermineResultsFileFormat(inputFilePath);
                    }
                }
                else
                {
                    peptideHitResultsFormat = mPeptideHitResultsFileFormat;
                    inputFileIsSynOrFht = false;
                }

                if (inputFileIsSynOrFht || peptideHitResultsFormat == ResultsFileFormat.AutoDetermine)
                {
                    // If peptideHitResultsFormat is still AutoDetermine that means we couldn't figure out the format

                    if (!inputFileIsSynOrFht)
                    {
                        ShowErrorMessage("Could not determine the format of the input file.");
                    }

                    Console.WriteLine();
                    ShowMessage(
                        "The filename must end in:\n" +
                        "  " + MSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE + ".tsv or " +
                        MSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE + ".tsv (for MS-GF+),\n" +
                        "  " + MSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE + ".txt (for MSAlign),\n" +
                        "  " + MODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE + ".txt (for MODA),\n" +
                        "  " + MODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE + ".txt (for MODPlus),\n" +
                        "  " + MSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE + ".tsv (for MSPathFinder),\n" +
                        "  " + SequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE + ".txt or " +
                        SequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE + ".txt (for SEQUEST),\n" +
                        "  " + TopPICResultsProcessor.FILENAME_SUFFIX_TopPIC_PRSMs_FILE + ".txt (for TopPIC), or\n" +
                        "  .xml (for X!Tandem).");

                    Console.WriteLine();
                    ShowMessage(ConsoleMsgUtils.WrapParagraph(
                                    "If the file ends in .tsv but does not match the other known .tsv suffixes shown above, " +
                                    "this program will assume the results come from MS-GF+"));

                    return false;
                }

                switch (peptideHitResultsFormat)
                {
                    case ResultsFileFormat.DiannParquetFile:
                    case ResultsFileFormat.DiannTSVFile:
                        mPeptideHitResultsProcessor = new DiaNNResultsProcessor(Options);
                        LogMessage("Detected DIA-NN results file");
                        break;

                    case ResultsFileFormat.InspectTXTFile:
                        mPeptideHitResultsProcessor = new InSpecTResultsProcessor(Options);
                        LogMessage("Detected InSpecT results file");
                        break;

                    case ResultsFileFormat.MaxQuantTXTFile:
                        mPeptideHitResultsProcessor = new MaxQuantResultsProcessor(Options);
                        LogMessage("Detected MaxQuant results file");
                        break;

                    case ResultsFileFormat.MODaTXTFile:
                        mPeptideHitResultsProcessor = new MODaResultsProcessor(Options);
                        LogMessage("Detected MODa results file");
                        break;

                    case ResultsFileFormat.MODPlusTXTFile:
                        mPeptideHitResultsProcessor = new MODPlusResultsProcessor(Options);
                        LogMessage("Detected MODPlus results file");
                        break;

                    case ResultsFileFormat.MSAlignTXTFile:
                        mPeptideHitResultsProcessor = new MSAlignResultsProcessor(Options);
                        LogMessage("Detected MSAlign results file");
                        break;

                    case ResultsFileFormat.MSFraggerTSVFile:
                        mPeptideHitResultsProcessor = new MSFraggerResultsProcessor(Options);
                        LogMessage("Detected MSFragger results file");
                        break;

                    case ResultsFileFormat.MSGFPlusTXTFile:
                        mPeptideHitResultsProcessor = new MSGFPlusResultsProcessor(Options);
                        LogMessage("Detected MS-GF+ results file");
                        break;

                    case ResultsFileFormat.MSPathFinderTSVFile:
                        mPeptideHitResultsProcessor = new MSPathFinderResultsProcessor(Options);
                        LogMessage("Detected MSPathFinder results file");
                        break;

                    case ResultsFileFormat.SequestFirstHitsFile:
                        mPeptideHitResultsProcessor = new SequestResultsProcessor(Options);
                        LogMessage("Detected SEQUEST First Hits file");
                        break;

                    case ResultsFileFormat.SequestSynopsisFile:
                        mPeptideHitResultsProcessor = new SequestResultsProcessor(Options);
                        LogMessage("Detected SEQUEST Synopsis file");
                        break;

                    case ResultsFileFormat.TopPICTXTFile:
                        mPeptideHitResultsProcessor = new TopPICResultsProcessor(Options);
                        LogMessage("Detected TopPIC results file");
                        break;

                    case ResultsFileFormat.XTandemXMLFile:
                        mPeptideHitResultsProcessor = new XTandemResultsProcessor(Options);
                        LogMessage("Detected X!Tandem XML file");
                        break;

                    case ResultsFileFormat.AutoDetermine:
                        throw new Exception("This code should not be reached; logic error in AutoDetermine: branch of switch (peptideHitResultsFormat)");

                    default:
                        // Unknown format
                        throw new Exception("This code should not be reached; logic error in default: branch of switch (peptideHitResultsFormat)");
                }

                if (mPeptideHitResultsProcessor == null)
                {
                    return false;
                }

                // Do not call RegisterEvents
                // Instead use local event handlers that optionally log to a file
                RegisterResultsProcessEvents(mPeptideHitResultsProcessor);

                InitializePeptideHitResultsProcessor(inputFilePath);

                if (!mOptionsDisplayed)
                {
                    ShowProcessingOptions(mPeptideHitResultsProcessor);
                    mOptionsDisplayed = true;
                }

                var success = mPeptideHitResultsProcessor.ProcessFile(inputFilePath, outputDirectoryPath, parameterFilePath);

                if (!success)
                {
                    ShowErrorMessage(mPeptideHitResultsProcessor.ErrorMessage);
                    return false;
                }

                LogMessage("Processing Complete");
                OperationComplete();
                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error calling ProcessFile in StartPHRP", ex);
                return false;
            }
        }

        private void UpdateLogFileOptions()
        {
            if (Options.LogFilePath != null && Options.LogFilePath.Equals("False", StringComparison.OrdinalIgnoreCase))
            {
                Options.LogMessagesToFile = false;
                Options.LogFilePath = string.Empty;
            }

            LogMessagesToFile = Options.LogMessagesToFile;

            if (Options.LogMessagesToFile)
            {
                LogFilePath = Options.LogFilePath == null ||
                              Options.LogFilePath.Equals("True", StringComparison.OrdinalIgnoreCase) ||
                              Options.LogFilePath.Equals(".")
                    ? string.Empty
                    : Options.LogFilePath;

                LogDirectoryPath = Options.LogDirectoryPath;
            }

            mLogFileOptionsUpdated = true;
        }

        private void PeptideHitResultsProcessor_ErrorOccurred(string message, Exception ex)
        {
            if (ex != null && !message.Contains(ex.Message))
            {
                LogMessage(message + ": " + ex.Message, MessageTypeConstants.ErrorMsg);
            }
            else
            {
                LogMessage(message, MessageTypeConstants.ErrorMsg);
            }

            if (ex != null)
            {
                Console.WriteLine(StackTraceFormatter.GetExceptionStackTraceMultiLine(ex));
            }
        }

        private void PeptideHitResultsProcessor_MessageEvent(string message)
        {
            LogMessage(message);
        }

        private void PeptideHitResultsProcessor_ProgressChanged(string taskDescription, float percentComplete)
        {
            UpdateProgress(taskDescription, percentComplete);
        }

        private void PeptideHitResultsProcessor_ProgressReset()
        {
            ResetProgress();
        }

        private void PeptideHitResultsProcessor_WarningMessageEvent(string warningMessage)
        {
            LogMessage(warningMessage, MessageTypeConstants.Warning);
        }
    }
}
