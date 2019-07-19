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
using PHRPReader;
using PRISM;

namespace PeptideHitResultsProcRunner
{
    /// <summary>
    /// This class calls clsSequestSynopsisFileProcessor or clsXTandemResultsConverter
    /// to process the files to determine the modifications present for each peptide,
    /// along with other information
    /// </summary>
    public class clsPeptideHitResultsProcRunner : PRISM.FileProcessor.ProcessFilesBase
    {
        public clsPeptideHitResultsProcRunner()
        {
            mFileDate = Program.PROGRAM_DATE;
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        /// <summary>
        /// Error codes specialized for this class
        /// </summary>
        public enum eResultsProcessorErrorCodes
        {
            NoError = 0,
            UnspecifiedError = -1
        }

        #endregion

        #region "Structures"

        #endregion

        #region "Classwide Variables"
        protected clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants mPeptideHitResultsFileFormat;

        protected bool mObtainModificationDefinitionsFromDMS;

        protected bool mUseExistingMTSPepToProteinMapFile;

        private clsPHRPBaseClass mPeptideHitResultsProcessor;

        protected eResultsProcessorErrorCodes mLocalErrorCode;
        #endregion

        #region "Properties"

        public bool CreateInspectOrMSGFDbFirstHitsFile { get; set; }

        public bool CreateInspectOrMSGFDbSynopsisFile { get; set; }

        public bool CreateProteinModsFile { get; set; }

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// Setting this to true assumes the input file is a valid PHRP data file
        /// Consequently, the code will only try to create the _ProteinMods.txt file, it will not re-create the PHRP data files
        /// When this is True, mCreateProteinModsFile is assumed to be true
        /// </remarks>
        public bool CreateProteinModsUsingPHRPDataFile { get; set; }

        public string FastaFilePath { get; set; }

        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        public float InspectSynopsisFilePValueThreshold { get; set; }

        public eResultsProcessorErrorCodes LocalErrorCode => mLocalErrorCode;

        public string MassCorrectionTagsFilePath { get; set; }

        public string ModificationDefinitionsFilePath { get; set; }

        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        public float MsgfPlusEValueThreshold { get; set; }

        public float MsgfPlusSpecEValueThreshold { get; set; }

        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        public string SearchToolParameterFilePath { get; set; }

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the the Fasta file specified by mFastaFilePath</remarks>
        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        #endregion

        public override IList<string> GetDefaultExtensionsToParse()
        {
            var extensionsToParse = new string[2];

            // Note: If mPeptideHitResultsFileFormat = .AutoDetermine
            //  then this class will only parse .txt files if they match
            // PeptideHitResultsProcessor.clsSequestResultsProcessor.SEQUEST_FIRST_HITS_FILE_SUFFIX or
            // PeptideHitResultsProcessor.clsSequestResultsProcessor.SEQUEST_SYNOPSIS_FILE_SUFFIX

            extensionsToParse[0] = ".txt";
            extensionsToParse[1] = ".xml";

            return extensionsToParse;
        }

        public override string GetErrorMessage()
        {
            // Returns "" if no error

            string errorMessage;

            if (ErrorCode == ProcessFilesErrorCodes.LocalizedError |
                ErrorCode == ProcessFilesErrorCodes.NoError)
            {
                switch (mLocalErrorCode)
                {
                    case eResultsProcessorErrorCodes.NoError:
                        errorMessage = string.Empty;
                        break;
                    case eResultsProcessorErrorCodes.UnspecifiedError:
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
            mPeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine;

            MassCorrectionTagsFilePath = string.Empty;
            ModificationDefinitionsFilePath = string.Empty;
            SearchToolParameterFilePath = string.Empty;

            CreateProteinModsFile = false;
            FastaFilePath = string.Empty;
            IgnorePeptideToProteinMapperErrors = false;
            ProteinModsFileIncludesReversedProteins = false;
            mUseExistingMTSPepToProteinMapFile = false;

            CreateInspectOrMSGFDbFirstHitsFile = false;
            CreateInspectOrMSGFDbSynopsisFile = false;
            InspectSynopsisFilePValueThreshold = clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MODaMODPlusSynopsisFileProbabilityThreshold = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            MsgfPlusEValueThreshold = clsMSGFPlusResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MsgfPlusSpecEValueThreshold = clsMSGFPlusResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            WarnMissingParameterFileSection = true;

            mLocalErrorCode = eResultsProcessorErrorCodes.NoError;
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

            mPeptideHitResultsProcessor.CreateInspectFirstHitsFile = CreateInspectOrMSGFDbFirstHitsFile;
            mPeptideHitResultsProcessor.CreateInspectSynopsisFile = CreateInspectOrMSGFDbSynopsisFile;
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
                        mPeptideHitResultsFileFormat = (clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants)value;
                    }
                    catch (Exception)
                    {
                        mPeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine;
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
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath, bool resetErrorCode)
        {

            var success = false;

            if (resetErrorCode)
            {
                SetLocalErrorCode(eResultsProcessorErrorCodes.NoError);
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
                        if (inputFilePath != null && inputFilePath.Contains(".."))
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
            }

            return success;
        }

        private void RegisterResultsProcessEvents(EventNotifier resultsProcessor)
        {
            resultsProcessor.ErrorEvent += PeptideHitResultsProcessor_ErrorOccurred;
            resultsProcessor.StatusEvent += PeptideHitResultsProcessor_MessageEvent;
            resultsProcessor.ProgressUpdate += PeptideHitResultsProcessor_ProgressChanged;
            resultsProcessor.WarningEvent += PeptideHitResultsProcessor_WarningMessageEvent;
        }

        /// <summary>
        /// Looks for file fileNameOrPath in the current working directory
        /// If not found, looks in sourceDirectoryPath
        /// </summary>
        /// <param name="sourceDirectoryPath">Path to the directory containing the input file</param>
        /// <param name="fileNameOrPath">File to find (either filename or full file path)</param>
        /// <returns>The path to the file if found, or fileNameOrPath if not found</returns>
        /// <remarks></remarks>
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

        private void SetLocalErrorCode(eResultsProcessorErrorCodes eNewErrorCode, bool leaveExistingErrorCodeUnchanged = false)
        {
            if (leaveExistingErrorCodeUnchanged && mLocalErrorCode != eResultsProcessorErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mLocalErrorCode = eNewErrorCode;

                if (eNewErrorCode == eResultsProcessorErrorCodes.NoError)
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

        private bool StartCreateProteinModsViaPHRPData(string inputFilePath, string outputDirectoryPath)
        {

            try
            {
                clsPHRPReader.ePeptideHitResultType ePeptideHitResultType;
                switch (mPeptideHitResultsFileFormat)
                {
                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Sequest;
                        LogMessage("Detected SEQUEST First Hits file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestSynopsisFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Sequest;
                        LogMessage("Detected SEQUEST Synopsis file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.XTandemXMLFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.XTandem;
                        LogMessage("Detected X!Tandem XML file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InspectTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Inspect;
                        LogMessage("Detected Inspect results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFPlusTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MSGFPlus;
                        LogMessage("Detected MSGF+ results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSAlignTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MSAlign;
                        LogMessage("Detected MSAlign results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODPlusTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MODPlus;
                        LogMessage("Detected MODPlus results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MSPathFinder;
                        LogMessage("Detected MSPathfinder results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.TopPICTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.TopPIC;
                        LogMessage("Detected TopPIC results file");
                        break;

                    default:
                        // Includes ePeptideHitResultsFileFormatConstants.AutoDetermine
                        ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(inputFilePath);
                        break;
                }

                if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
                {
                    ShowErrorMessage("Error: Could not determine the format of the PHRP data file: " + inputFilePath);
                    return false;
                }

                switch (ePeptideHitResultType)
                {
                    case clsPHRPReader.ePeptideHitResultType.Sequest:
                        mPeptideHitResultsProcessor = new clsSequestResultsProcessor();
                        break;

                    case clsPHRPReader.ePeptideHitResultType.XTandem:
                        mPeptideHitResultsProcessor = new clsXTandemResultsProcessor();
                        break;

                    case clsPHRPReader.ePeptideHitResultType.Inspect:
                        mPeptideHitResultsProcessor = new clsInSpecTResultsProcessor();
                        break;

                    case clsPHRPReader.ePeptideHitResultType.MSGFPlus:
                        mPeptideHitResultsProcessor = new clsMSGFPlusResultsProcessor();
                        break;

                    case clsPHRPReader.ePeptideHitResultType.MSAlign:
                        mPeptideHitResultsProcessor = new clsMSAlignResultsProcessor();
                        break;

                    case clsPHRPReader.ePeptideHitResultType.MODa:
                        mPeptideHitResultsProcessor = new clsMODaResultsProcessor();
                        break;

                    default:
                        // Unknown format
                        ShowErrorMessage("Error: Unrecognized value for ePeptideHitResultType: " + ePeptideHitResultType);
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
                clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants ePeptideHitResultsFormat;
                if (mPeptideHitResultsFileFormat == clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine)
                {
                    ePeptideHitResultsFormat = clsPHRPBaseClass.DetermineResultsFileFormat(inputFilePath);
                }
                else
                {
                    ePeptideHitResultsFormat = mPeptideHitResultsFileFormat;
                }

                if (ePeptideHitResultsFormat == clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine)
                {
                    // If ePeptideHitResultsFormat is still AutoDetermine that means we couldn't figure out the format

                    ShowErrorMessage("Could not determine the format of the input file.");

                    Console.WriteLine();
                    ShowMessage(
                        "The filename must end in:\n" +
                        "  " + clsMSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE + ".tsv or " +
                        clsMSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE + ".tsv (for MS-GF+),\n" +
                        "  " + clsMSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE + ".txt (for MSAlign),\n" +
                        "  " + clsMODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE + ".txt (for MODA),\n" +
                        "  " + clsMODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE + ".txt (for MODPlus),\n" +
                        "  " + clsMSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE + ".tsv (for MSPathFinder),\n" +
                        "  " + clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE + ".txt or " +
                        clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE + ".txt (for SEQUEST),\n" +
                        "  " + clsTopPICResultsProcessor.FILENAME_SUFFIX_TopPIC_PRSMs_FILE + ".txt (for TopPIC), or\n" +
                        "  .xml (for X!Tandem).");

                    Console.WriteLine();
                    ShowMessage(ConsoleMsgUtils.WrapParagraph(
                                    "If the file ends in .tsv but does not match the other known .tsv suffixes shown above, " +
                                    "this program will assume the results come from MS-GF+"));

                    return false;
                }

                switch (ePeptideHitResultsFormat)
                {
                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile:
                        mPeptideHitResultsProcessor = new clsSequestResultsProcessor();
                        LogMessage("Detected SEQUEST First Hits file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestSynopsisFile:
                        mPeptideHitResultsProcessor = new clsSequestResultsProcessor();
                        LogMessage("Detected SEQUEST Synopsis file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.XTandemXMLFile:
                        mPeptideHitResultsProcessor = new clsXTandemResultsProcessor();
                        LogMessage("Detected X!Tandem XML file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InspectTXTFile:
                        mPeptideHitResultsProcessor = new clsInSpecTResultsProcessor();
                        LogMessage("Detected Inspect results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFPlusTXTFile:
                        mPeptideHitResultsProcessor = new clsMSGFPlusResultsProcessor();
                        LogMessage("Detected MSGF+ results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSAlignTXTFile:
                        mPeptideHitResultsProcessor = new clsMSAlignResultsProcessor();
                        LogMessage("Detected MSAlign results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODaTXTFile:
                        mPeptideHitResultsProcessor = new clsMODaResultsProcessor();
                        LogMessage("Detected MODa results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODPlusTXTFile:
                        mPeptideHitResultsProcessor = new clsMODPlusResultsProcessor();
                        LogMessage("Detected MODPlus results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile:
                        mPeptideHitResultsProcessor = new clsMSPathFinderResultsProcessor();
                        LogMessage("Detected MSPathFinder results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.TopPICTXTFile:
                        mPeptideHitResultsProcessor = new clsTopPICResultsProcessor();
                        LogMessage("Detected TopPIC results file");
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine:
                        throw new Exception("This code should not be reached; logic error in AutoDetermine: branch of switch (ePeptideHitResultsFormat)");
                    default:
                        // Unknown format
                        throw new Exception("This code should not be reached; logic error in default: branch of switch (ePeptideHitResultsFormat)");
                }

                if (mPeptideHitResultsProcessor == null)
                {
                    return false;
                }

                // Do not call RegisterEvents
                // Instead use local event handlers that optionally log to a file
                RegisterResultsProcessEvents(mPeptideHitResultsProcessor);

                InitializePeptideHitResultsProcessor(inputFilePath);

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

        private void PeptideHitResultsProcessor_WarningMessageEvent(string warningMessage)
        {
            LogMessage(warningMessage, MessageTypeConstants.Warning);
        }

    }
}
