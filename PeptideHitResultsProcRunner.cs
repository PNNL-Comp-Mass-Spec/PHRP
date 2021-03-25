// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or https://www.pnnl.gov/sysbio/ or https://panomics.pnnl.gov/
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
        // Ignore Spelling: Battelle, parm, proc, MODa, enums, fasta

        /// <summary>
        /// Constructor
        /// </summary>
        public PeptideHitResultsProcRunner()
        {
            mFileDate = Program.PROGRAM_DATE;
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        /// <summary>
        /// Error codes specialized for this class
        /// </summary>
        public enum ResultsProcessorErrorCodes
        {
            NoError = 0,
            UnspecifiedError = -1
        }

        #endregion

        #region "Structures"

        #endregion

        #region "Class wide Variables"

        protected ResultsFileFormat mPeptideHitResultsFileFormat;

        protected bool mObtainModificationDefinitionsFromDMS;

        protected bool mUseExistingMTSPepToProteinMapFile;

        private PHRPBaseClass mPeptideHitResultsProcessor;

        protected ResultsProcessorErrorCodes mLocalErrorCode;

        private bool mOptionsDisplayed = false;

        private bool mFilePathsShown = false;

        #endregion

        #region "Properties"

        public bool CreateFirstHitsFile { get; set; }

        public bool CreateSynopsisFile { get; set; }

        public bool CreateProteinModsFile { get; set; }

        /// <summary>
        /// Create Protein Mods Using PHRP Data File
        /// </summary>
        /// <remarks>
        /// Setting this to true assumes the input file is a valid PHRP data file
        /// Consequently, the code will only try to create the _ProteinMods.txt file, it will not re-create the PHRP data files
        /// When this is True, mCreateProteinModsFile is assumed to be true
        /// </remarks>
        public bool CreateProteinModsUsingPHRPDataFile { get; set; }

        public string FastaFilePath { get; set; }

        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        public float InspectSynopsisFilePValueThreshold { get; set; }

        public ResultsProcessorErrorCodes LocalErrorCode => mLocalErrorCode;

        public string MassCorrectionTagsFilePath { get; set; }

        public string ModificationDefinitionsFilePath { get; set; }

        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        public float MsgfPlusEValueThreshold { get; set; }

        public float MsgfPlusSpecEValueThreshold { get; set; }

        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        public string SearchToolParameterFilePath { get; set; }

        /// <summary>
        /// Use Existing MTS PepToProtein Map File
        /// </summary>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the Fasta file specified by mFastaFilePath</remarks>
        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        #endregion

        /// <summary>
        /// If filePath is an empty string, return textIfEmptyPath
        /// Otherwise, return the file path, truncated if over the specified maximum length
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="textIfEmptyPath"></param>
        /// <param name="maxPathLength"></param>
        private string FilePathOrText(string filePath, string textIfEmptyPath, int maxPathLength = 110)
        {
            if (string.IsNullOrWhiteSpace(filePath))
                return textIfEmptyPath;

            return PathUtils.CompactPathString(filePath, maxPathLength);
        }

        public override IList<string> GetDefaultExtensionsToParse()
        {
            var extensionsToParse = new string[2];

            // Note: If mPeptideHitResultsFileFormat = .AutoDetermine
            //  then this class will only parse .txt files if they match
            // PeptideHitResultsProcessor.SequestResultsProcessor.SEQUEST_FIRST_HITS_FILE_SUFFIX or
            // PeptideHitResultsProcessor.SequestResultsProcessor.SEQUEST_SYNOPSIS_FILE_SUFFIX

            extensionsToParse[0] = ".txt";
            extensionsToParse[1] = ".xml";

            return extensionsToParse;
        }

        public override string GetErrorMessage()
        {
            // Returns "" if no error

            string errorMessage;

            if (ErrorCode == ProcessFilesErrorCodes.LocalizedError ||
                ErrorCode == ProcessFilesErrorCodes.NoError)
            {
                switch (mLocalErrorCode)
                {
                    case ResultsProcessorErrorCodes.NoError:
                        errorMessage = string.Empty;
                        break;
                    case ResultsProcessorErrorCodes.UnspecifiedError:
                        errorMessage = "Unspecified localized error";
                        break;
                    default:
                        // This shouldn't happen
                        errorMessage = "Unknown error state";
                        break;
                }
            }
            else
            {
                errorMessage = GetBaseClassErrorMessage();
            }

            return errorMessage;
        }

        private void InitializeLocalVariables()
        {
            mPeptideHitResultsFileFormat = ResultsFileFormat.AutoDetermine;

            MassCorrectionTagsFilePath = string.Empty;
            ModificationDefinitionsFilePath = string.Empty;
            SearchToolParameterFilePath = string.Empty;

            CreateProteinModsFile = false;
            FastaFilePath = string.Empty;
            IgnorePeptideToProteinMapperErrors = false;
            ProteinModsFileIncludesReversedProteins = false;
            mUseExistingMTSPepToProteinMapFile = false;

            CreateFirstHitsFile = false;
            CreateSynopsisFile = false;
            InspectSynopsisFilePValueThreshold = InSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MODaMODPlusSynopsisFileProbabilityThreshold = MODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            MsgfPlusEValueThreshold = MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MsgfPlusSpecEValueThreshold = MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            WarnMissingParameterFileSection = true;

            mLocalErrorCode = ResultsProcessorErrorCodes.NoError;
        }

        private void InitializePeptideHitResultsProcessor(string inputFilePath)
        {
            if (mObtainModificationDefinitionsFromDMS)
            {
                LoadModificationInfoFromDMS();
            }

            var sourceFile = new FileInfo(inputFilePath);

            mPeptideHitResultsProcessor.MassCorrectionTagsFilePath = ResolveFilePath(sourceFile.DirectoryName, MassCorrectionTagsFilePath);
            mPeptideHitResultsProcessor.ModificationDefinitionsFilePath = ResolveFilePath(sourceFile.DirectoryName, ModificationDefinitionsFilePath);
            mPeptideHitResultsProcessor.SearchToolParameterFilePath = ResolveFilePath(sourceFile.DirectoryName, SearchToolParameterFilePath);

            mPeptideHitResultsProcessor.CreateProteinModsFile = CreateProteinModsFile;
            mPeptideHitResultsProcessor.FastaFilePath = FastaFilePath;
            mPeptideHitResultsProcessor.IgnorePeptideToProteinMapperErrors = IgnorePeptideToProteinMapperErrors;
            mPeptideHitResultsProcessor.ProteinModsFileIncludesReversedProteins = ProteinModsFileIncludesReversedProteins;
            mPeptideHitResultsProcessor.UseExistingMTSPepToProteinMapFile = mUseExistingMTSPepToProteinMapFile;

            mPeptideHitResultsProcessor.CreateFirstHitsFile = CreateFirstHitsFile;
            mPeptideHitResultsProcessor.CreateSynopsisFile = CreateSynopsisFile;
            mPeptideHitResultsProcessor.InspectSynopsisFilePValueThreshold = InspectSynopsisFilePValueThreshold;

            mPeptideHitResultsProcessor.MODaMODPlusSynopsisFileProbabilityThreshold = MODaMODPlusSynopsisFileProbabilityThreshold;

            mPeptideHitResultsProcessor.MSGFPlusSynopsisFileEValueThreshold = MsgfPlusEValueThreshold;
            mPeptideHitResultsProcessor.MSGFPlusSynopsisFileSpecEValueThreshold = MsgfPlusSpecEValueThreshold;

            mPeptideHitResultsProcessor.WarnMissingParameterFileSection = WarnMissingParameterFileSection;
        }

        private void LoadModificationInfoFromDMS()
        {
            // ToDo: Contact DMS to get the modification information
            // The results of this query will need to be filtered to get the info for just this analysis job

            // ReSharper disable CommentTypo

            //SELECT D.Dataset_Num,
            //    AJ.AJ_jobID, PFMI.Local_Symbol,
            //    PFMI.Monoisotopic_Mass, PFMI.Residue_Symbol,
            //    PFMI.Mod_Type_Symbol, PFMI.Mass_Correction_Tag
            //FROM dbo.T_Analysis_Job AJ INNER JOIN
            //    dbo.T_Dataset D ON
            //    AJ.AJ_datasetID = D.Dataset_ID LEFT
            //     OUTER JOIN
            //    dbo.V_Param_File_Mass_Mod_Info PFMI ON
            //    AJ.AJ_parmFileName = PFMI.Param_File_Name
            //WHERE (D.Dataset_Num = 'QC_05_2_a_24Oct05_Doc_0508-08')
            //ORDER BY AJ.AJ_jobID, PFMI.Local_Symbol

            //SELECT D.Dataset_Num,
            //    AJ.AJ_jobID, PFMI.Local_Symbol,
            //    PFMI.Monoisotopic_Mass, PFMI.Residue_Symbol,
            //    PFMI.Mod_Type_Symbol, PFMI.Mass_Correction_Tag
            //FROM dbo.T_Analysis_Job AJ INNER JOIN
            //    dbo.T_Dataset D ON
            //    AJ.AJ_datasetID = D.Dataset_ID LEFT
            //     OUTER JOIN
            //    dbo.V_Param_File_Mass_Mod_Info PFMI ON
            //    AJ.AJ_parmFileName = PFMI.Param_File_Name
            //WHERE (AJ.AJ_jobID = 47703)
            //ORDER BY AJ.AJ_jobID, PFMI.Local_Symbol

            // ReSharper restore CommentTypo
        }

        private bool LoadParameterFileSettings(string parameterFilePath)
        {
            const string OPTIONS_SECTION = "PeptideHitResultsProcRunner";

            var settingsFile = new XmlSettingsFileAccessor();

            try
            {
                if (string.IsNullOrWhiteSpace(parameterFilePath))
                {
                    // No parameter file specified; nothing to load
                    return true;
                }

                if (!File.Exists(parameterFilePath))
                {
                    // See if parameterFilePath points to a file in the same directory as the application
                    var appDirPath = GetAppDirectoryPath();
                    if (string.IsNullOrWhiteSpace(appDirPath))
                    {
                        SetBaseClassErrorCode(ProcessFilesErrorCodes.ParameterFileNotFound);
                        return false;
                    }

                    parameterFilePath = Path.Combine(appDirPath, Path.GetFileName(parameterFilePath));
                    if (!File.Exists(parameterFilePath))
                    {
                        SetBaseClassErrorCode(ProcessFilesErrorCodes.ParameterFileNotFound);
                        return false;
                    }
                }

                if (settingsFile.LoadSettings(parameterFilePath))
                {
                    if (!settingsFile.SectionPresent(OPTIONS_SECTION))
                    {
                        ShowErrorMessage("The node '<section name=\"" + OPTIONS_SECTION + "\"> was not found in the parameter file: " + parameterFilePath);
                        SetBaseClassErrorCode(ProcessFilesErrorCodes.InvalidParameterFile);
                        return false;
                    }

                    mObtainModificationDefinitionsFromDMS = settingsFile.GetParam(OPTIONS_SECTION, "ObtainModificationDefinitionsFromDMS", mObtainModificationDefinitionsFromDMS);

                    var value = settingsFile.GetParam(OPTIONS_SECTION, "PeptideHitResultsFileFormat", Convert.ToInt32(mPeptideHitResultsFileFormat));
                    try
                    {
                        mPeptideHitResultsFileFormat = (ResultsFileFormat)value;
                    }
                    catch (Exception)
                    {
                        mPeptideHitResultsFileFormat = ResultsFileFormat.AutoDetermine;
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in LoadParameterFileSettings", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Log an additional message to the log file
        /// </summary>
        /// <param name="message"></param>
        public void LogAdditionalMessage(string message)
        {
            LogMessage(message);
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="parameterFilePath"></param>
        /// <param name="resetErrorCode"></param>
        /// <returns>True if successful, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath, bool resetErrorCode)
        {
            var success = false;

            if (resetErrorCode)
            {
                SetLocalErrorCode(ResultsProcessorErrorCodes.NoError);
            }

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                var statusMessage = "Parameter file load error: " + parameterFilePath;
                ShowErrorMessage(statusMessage);

                if (ErrorCode == ProcessFilesErrorCodes.NoError)
                {
                    SetBaseClassErrorCode(ProcessFilesErrorCodes.InvalidParameterFile);
                }
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    ShowErrorMessage("Input file name is empty");
                    SetBaseClassErrorCode(ProcessFilesErrorCodes.InvalidInputFilePath);
                }
                else
                {
                    // Note that CleanupFilePaths() will update mOutputDirectoryPath, which is used by LogMessage()
                    if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                    {
                        SetBaseClassErrorCode(ProcessFilesErrorCodes.FilePathError);
                        if (inputFilePath?.Contains("..") == true)
                        {
                            var inputFile = new FileInfo(inputFilePath);
                            OnStatusEvent("Absolute path: " + inputFile.DirectoryName);
                        }
                    }
                    else
                    {
                        UpdateProgress("Parsing " + Path.GetFileName(inputFilePath));
                        ResetProgress();

                        if (CreateProteinModsUsingPHRPDataFile)
                        {
                            success = StartCreateProteinModsViaPHRPData(inputFilePath, outputDirectoryPath);
                        }
                        else
                        {
                            success = StartPHRP(inputFilePath, outputDirectoryPath, parameterFilePath);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFile", ex);
                OnDebugEvent(StackTraceFormatter.GetExceptionStackTraceMultiLine(ex));
            }

            return success;
        }

        private void RegisterResultsProcessEvents(PHRPBaseClass resultsProcessor)
        {
            resultsProcessor.ErrorEvent += PeptideHitResultsProcessor_ErrorOccurred;
            resultsProcessor.StatusEvent += PeptideHitResultsProcessor_MessageEvent;
            resultsProcessor.ProgressUpdate += PeptideHitResultsProcessor_ProgressChanged;
            resultsProcessor.WarningEvent += PeptideHitResultsProcessor_WarningMessageEvent;
            resultsProcessor.ProgressReset += PeptideHitResultsProcessor_ProgressReset;
        }

        /// <summary>
        /// Looks for fileNameOrPath in the current working directory
        /// If not found, looks in sourceDirectoryPath
        /// </summary>
        /// <param name="sourceDirectoryPath">Path to the directory containing the input file</param>
        /// <param name="fileNameOrPath">File to find (either filename or full file path)</param>
        /// <returns>The path to the file if found, or fileNameOrPath if not found</returns>
        protected string ResolveFilePath(string sourceDirectoryPath, string fileNameOrPath)
        {
            if (File.Exists(fileNameOrPath))
            {
                return fileNameOrPath;
            }

            var fileName = Path.GetFileName(fileNameOrPath);
            if (string.IsNullOrWhiteSpace(fileName))
                return fileNameOrPath;

            var newFilePath = Path.Combine(sourceDirectoryPath, fileName);
            if (File.Exists(newFilePath))
            {
                return newFilePath;
            }

            return fileNameOrPath;
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

            if(!string.IsNullOrWhiteSpace(parameterFilePath))
            {
                LogMessage(string.Format("{0,-45} {1}",
                    "Parameter file:", PathUtils.CompactPathString(parameterFilePath, 110)));
            }
        }

        private void ShowProcessingOptions(PHRPBaseClass resultsProcessor)
        {
            Console.WriteLine();
            LogMessage("Processing options for " + resultsProcessor);

            LogMessage(string.Format("{0,-45} {1}",
                "Search Tool Parameter File:", FilePathOrText(resultsProcessor.SearchToolParameterFilePath, "Not defined")));

            LogMessage(string.Format("{0,-45} {1}",
                "Modification Definitions File:", FilePathOrText(resultsProcessor.ModificationDefinitionsFilePath, "Not defined")));

            LogMessage(string.Format("{0,-45} {1}",
                "Mass Correction Tags File:", FilePathOrText(resultsProcessor.MassCorrectionTagsFilePath, "Use internally defined tags")));

            LogMessage(string.Format("{0,-45} {1}",
                "FASTA File:", FilePathOrText(resultsProcessor.FastaFilePath, "Not defined")));

            Console.WriteLine();
            LogMessage(string.Format("{0,-45} {1}",
                "Create Protein Mods File:", resultsProcessor.CreateProteinModsFile));

            LogMessage(string.Format("{0,-45} {1}",
                "Ignore Peptide to Protein Mapper Errors:", resultsProcessor.IgnorePeptideToProteinMapperErrors));

            LogMessage(string.Format("{0,-45} {1}",
                "Protein Mods File Includes Reversed Proteins:", resultsProcessor.ProteinModsFileIncludesReversedProteins));

            LogMessage(string.Format("{0,-45} {1}",
                "Use Existing MTS PepToProtein Map File:", resultsProcessor.UseExistingMTSPepToProteinMapFile));

            Console.WriteLine();
            LogMessage(string.Format("{0,-45} {1}",
                "Create First Hits File:", resultsProcessor.CreateFirstHitsFile));

            LogMessage(string.Format("{0,-45} {1}",
                "Create Synopsis File:", resultsProcessor.CreateSynopsisFile));

            if (resultsProcessor is InSpecTResultsProcessor)
            {
                LogMessage(string.Format("{0,-45} {1:E2}",
                    "Inspect Synopsis File PValue Threshold:", resultsProcessor.InspectSynopsisFilePValueThreshold));
            }

            Console.WriteLine();
            LogMessage(string.Format("{0,-49} {1:E2}",
                "MODa/MODPlus Synopsis File Probability Threshold:", resultsProcessor.MODaMODPlusSynopsisFileProbabilityThreshold));

            LogMessage(string.Format("{0,-49} {1:E2}",
                "MSGFPlus Synopsis File EValue Threshold:", resultsProcessor.MSGFPlusSynopsisFileEValueThreshold));

            LogMessage(string.Format("{0,-49} {1:E2}",
                "MSGFPlus Synopsis File SpecEValue Threshold:", resultsProcessor.MSGFPlusSynopsisFileSpecEValueThreshold));

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
                        LogMessage("Detected Inspect results file");
                        break;

                    case ResultsFileFormat.MSGFPlusTXTFile:
                        PeptideHitResultType = PeptideHitResultTypes.MSGFPlus;
                        LogMessage("Detected MSGF+ results file");
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
                        mPeptideHitResultsProcessor = new SequestResultsProcessor();
                        break;

                    case PeptideHitResultTypes.XTandem:
                        mPeptideHitResultsProcessor = new XTandemResultsProcessor();
                        break;

                    case PeptideHitResultTypes.Inspect:
                        mPeptideHitResultsProcessor = new InSpecTResultsProcessor();
                        break;

                    case PeptideHitResultTypes.MSGFPlus:
                        mPeptideHitResultsProcessor = new MSGFPlusResultsProcessor();
                        break;

                    case PeptideHitResultTypes.MSAlign:
                        mPeptideHitResultsProcessor = new MSAlignResultsProcessor();
                        break;

                    case PeptideHitResultTypes.MODa:
                        mPeptideHitResultsProcessor = new MODaResultsProcessor();
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
                if (mPeptideHitResultsFileFormat == ResultsFileFormat.AutoDetermine)
                {
                    peptideHitResultsFormat = PHRPBaseClass.DetermineResultsFileFormat(inputFilePath);
                }
                else
                {
                    peptideHitResultsFormat = mPeptideHitResultsFileFormat;
                }

                if (peptideHitResultsFormat == ResultsFileFormat.AutoDetermine)
                {
                    // If peptideHitResultsFormat is still AutoDetermine that means we couldn't figure out the format

                    ShowErrorMessage("Could not determine the format of the input file.");

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
                    case ResultsFileFormat.SequestFirstHitsFile:
                        mPeptideHitResultsProcessor = new SequestResultsProcessor();
                        LogMessage("Detected SEQUEST First Hits file");
                        break;

                    case ResultsFileFormat.SequestSynopsisFile:
                        mPeptideHitResultsProcessor = new SequestResultsProcessor();
                        LogMessage("Detected SEQUEST Synopsis file");
                        break;

                    case ResultsFileFormat.XTandemXMLFile:
                        mPeptideHitResultsProcessor = new XTandemResultsProcessor();
                        LogMessage("Detected X!Tandem XML file");
                        break;

                    case ResultsFileFormat.InspectTXTFile:
                        mPeptideHitResultsProcessor = new InSpecTResultsProcessor();
                        LogMessage("Detected Inspect results file");
                        break;

                    case ResultsFileFormat.MSGFPlusTXTFile:
                        mPeptideHitResultsProcessor = new MSGFPlusResultsProcessor();
                        LogMessage("Detected MSGF+ results file");
                        break;

                    case ResultsFileFormat.MSAlignTXTFile:
                        mPeptideHitResultsProcessor = new MSAlignResultsProcessor();
                        LogMessage("Detected MSAlign results file");
                        break;

                    case ResultsFileFormat.MODaTXTFile:
                        mPeptideHitResultsProcessor = new MODaResultsProcessor();
                        LogMessage("Detected MODa results file");
                        break;

                    case ResultsFileFormat.MODPlusTXTFile:
                        mPeptideHitResultsProcessor = new MODPlusResultsProcessor();
                        LogMessage("Detected MODPlus results file");
                        break;

                    case ResultsFileFormat.MSPathFinderTSVFile:
                        mPeptideHitResultsProcessor = new MSPathFinderResultsProcessor();
                        LogMessage("Detected MSPathFinder results file");
                        break;

                    case ResultsFileFormat.TopPICTXTFile:
                        mPeptideHitResultsProcessor = new TopPICResultsProcessor();
                        LogMessage("Detected TopPIC results file");
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
