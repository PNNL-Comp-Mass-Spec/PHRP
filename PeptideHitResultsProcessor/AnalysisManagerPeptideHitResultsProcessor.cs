// This class implements the IPeptideHitResultsProcessor interface
// It uses class XTandemResultsProcessor to read an XTandem search results XML file
//  or class SequestResultsProcessor to read a Sequest Synopsis or First Hits file
//
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// Started January 6, 2006

using System;
using System.IO;
using System.Threading;
using PeptideHitResultsProcessor.Processor;

namespace PeptideHitResultsProcessor
{
    [Obsolete("This class is unused")]
    public class AnalysisManagerPeptideHitResultsProcessor : IPeptideHitResultsProcessor
    {
        // Ignore Spelling: msgfdb

        #region "Constants and enums"
        private const string DEFAULT_MASS_CORRECTION_TAGS_FILENAME = "Mass_Correction_Tags.txt";
        private const string MODIFICATION_DEFINITIONS_FILE_SUFFIX = "_ModDefs.txt";
        #endregion

        #region "Class wide Variables"

        private Enums.ResultsFileFormat m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.AutoDetermine;

        private string m_ParameterFilePath = string.Empty;      // Peptide search tool parameter file path
        private string m_SettingsFilePath = string.Empty;       // XML settings file with section PeptideHitResultsProcessorOptions

        private string m_PeptideHitResultsFilePath = string.Empty;
        private string m_MassCorrectionTagsFilePath = string.Empty;
        private string m_ModificationDefinitionsFilePath = string.Empty;

        private string m_ErrMsg = string.Empty;

        private ProcessStatus m_Status;
        private ProcessResults m_Results;

        private PHRPBaseClass m_PeptideHitResultsProcessor;

        private Thread m_thThread;

        #endregion

        #region "Properties"

        public override string ErrMsg => m_ErrMsg;

        public override float PercentComplete
        {
            get
            {
                if (m_PeptideHitResultsProcessor != null)
                {
                    return m_PeptideHitResultsProcessor.ProgressPercentComplete;
                }
                return 0;
            }
        }

        public override ProcessStatus Status => m_Status;

        public override ProcessResults Results => m_Results;

        #endregion

        public override ProcessStatus Abort()
        {
            if (m_PeptideHitResultsProcessor != null)
            {
                m_Status = ProcessStatus.PH_ABORTING;
                m_PeptideHitResultsProcessor.AbortProcessingNow();
            }
            return m_Status;
        }

        /// <summary>
        /// Copies all input data required for plugin operation to appropriate memory variables
        /// </summary>
        /// <param name="options"></param>
        public override void Setup(InitializationParams options)
        {
            SourceDirectoryPath = options.SourceDirectoryPath;
            OutputDirectoryPath = options.OutputDirectoryPath;

            PeptideHitResultsFileName = options.PeptideHitResultsFileName;
            MassCorrectionTagsFileName = options.MassCorrectionTagsFileName;
            ModificationDefinitionsFileName = options.ModificationDefinitionsFileName;

            // This is unused and thus obsolete
            // m_MiscParams = .MiscParams

            DebugLevel = options.DebugLevel;

            AnalysisToolName = options.AnalysisToolName;
            DatasetName = options.DatasetName;
            ParameterFileName = options.ParameterFileName;
            SettingsFileName = options.SettingsFileName;

            CreateFirstHitsFile = options.CreateFirstHitsFile;
            CreateSynopsisFile = options.CreateSynopsisFile;
        }

        public override ProcessStatus Start()
        {
            m_Status = ProcessStatus.PH_STARTING;

            // Verify necessary files are in specified locations
            if (!InitSetup())
            {
                m_Results = ProcessResults.PH_FAILURE;
                m_Status = ProcessStatus.PH_ERROR;
                return m_Status;
            }

            // Process the results file (the process runs in a separate thread)
            m_Status = ProcessPeptideHitResultsFile();

            if (m_Status == ProcessStatus.PH_ERROR)
            {
                m_Results = ProcessResults.PH_FAILURE;
            }

            return m_Status;
        }

        /// <summary>
        /// Initializes m_PeptideHitResultsProcessor then starts a separate thread to process the file
        /// </summary>
        protected virtual ProcessStatus ProcessPeptideHitResultsFile()
        {
            try
            {
                // Initialize m_PeptideHitResultsProcessor
                switch (m_PeptideHitResultsFileFormat)
                {
                    case Enums.ResultsFileFormat.XTandemXMLFile:
                        m_PeptideHitResultsProcessor = new XtandemResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.SequestFirstHitsFile:
                    case Enums.ResultsFileFormat.SequestSynopsisFile:
                        m_PeptideHitResultsProcessor = new SequestResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.InspectTXTFile:
                        m_PeptideHitResultsProcessor = new InSpecTResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.MSGFPlusTXTFile:
                        m_PeptideHitResultsProcessor = new MSGFPlusResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.MSAlignTXTFile:
                        m_PeptideHitResultsProcessor = new MSAlignResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.MODaTXTFile:
                        m_PeptideHitResultsProcessor = new MODaResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.MODPlusTXTFile:
                        m_PeptideHitResultsProcessor = new MODaResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.MSPathFinderTSVFile:
                        m_PeptideHitResultsProcessor = new MSPathFinderResultsProcessor();
                        break;

                    case Enums.ResultsFileFormat.TopPICTXTFile:
                        m_PeptideHitResultsProcessor = new TopPICResultsProcessor();
                        break;

                    default:
                        // Unknown format; cannot continue
                        LogErrors("ProcessPeptideHitResultsFile", "Unknown peptide hit results file format: " + m_PeptideHitResultsFileFormat.ToString(), null);
                        m_Status = ProcessStatus.PH_ERROR;
                        return m_Status;
                }

                RegisterEvents(m_PeptideHitResultsProcessor);

                m_PeptideHitResultsProcessor.ProgressComplete += mPeptideHitResultsProcessor_ProgressComplete;
                m_PeptideHitResultsProcessor.ProgressReset += mPeptideHitResultsProcessor_ProgressReset;

                // Define the auxiliary file paths
                m_PeptideHitResultsProcessor.MassCorrectionTagsFilePath = m_MassCorrectionTagsFilePath;
                m_PeptideHitResultsProcessor.ModificationDefinitionsFilePath = m_ModificationDefinitionsFilePath;
                m_PeptideHitResultsProcessor.SearchToolParameterFilePath = m_ParameterFilePath;
                m_PeptideHitResultsProcessor.InspectSynopsisFilePValueThreshold = InSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

                m_PeptideHitResultsProcessor.CreateFirstHitsFile = CreateFirstHitsFile;
                m_PeptideHitResultsProcessor.CreateSynopsisFile = CreateSynopsisFile;

                m_thThread = new Thread(ProcessPeptideHitResultsFileWork);
                m_thThread.Start();

                m_Status = ProcessStatus.PH_RUNNING;
            }
            catch (Exception ex)
            {
                LogErrors("ProcessPeptideHitResultsFile", "Error initializing and running m_PeptideHitResultsProcessor", ex);
                m_Status = ProcessStatus.PH_ERROR;
            }

            return m_Status;
        }

        protected virtual void ProcessPeptideHitResultsFileWork()
        {
            bool success;

            try
            {
                success = m_PeptideHitResultsProcessor.ProcessFile(m_PeptideHitResultsFilePath, OutputDirectoryPath, m_SettingsFilePath);
            }
            catch (Exception)
            {
                success = false;
            }

            if (success)
            {
                m_Status = ProcessStatus.PH_COMPLETE;
            }
            else
            {
                if (m_PeptideHitResultsProcessor.AbortProcessing)
                {
                    LogErrors("ProcessPeptideHitResultsFileWork", "Processing aborted", null);
                    m_Results = ProcessResults.PH_ABORTED;
                    m_Status = ProcessStatus.PH_ABORTING;
                }
                else
                {
                    LogErrors("ProcessPeptideHitResultsFileWork", m_PeptideHitResultsProcessor.ErrorMessage, null);
                    m_Results = ProcessResults.PH_FAILURE;
                    m_Status = ProcessStatus.PH_ERROR;
                }
            }
        }

        /// <summary>
        /// Initializes module variables and verifies that mandatory parameters have been properly specified
        /// </summary>
        /// <returns>True if successful, false if an error</returns>
        protected virtual bool InitSetup()
        {
            // Output directory name
            if (string.IsNullOrEmpty(OutputDirectoryPath))
            {
                m_ErrMsg = "Output directory path not specified";
                return false;
            }

            // Source directory name
            if (string.IsNullOrEmpty(SourceDirectoryPath))
            {
                m_ErrMsg = "Source directory not specified";
                return false;
            }

            if (DebugLevel >= 3)
            {
                OnDebugEvent("Setup params: OutputDirectoryPath = " + OutputDirectoryPath);
                OnDebugEvent("Setup params: SourceDirectoryPath = " + SourceDirectoryPath);
            }

            // Source directory exists?
            if (!VerifyDirExists(SourceDirectoryPath))
                return false; // Error msg handled by VerifyDirExists

            // Output directory exists?
            if (!VerifyDirExists(OutputDirectoryPath))
                return false; // Error msg handled by VerifyDirExists

            // Analysis tool name defined?
            if (string.IsNullOrWhiteSpace(AnalysisToolName))
            {
                m_ErrMsg = "Analysis tool name not specified";
                return false;
            }

            // Dataset name defined?
            if (string.IsNullOrWhiteSpace(DatasetName))
            {
                m_ErrMsg = "Dataset name not specified";
                return false;
            }

            // Settings file name defined?
            if (string.IsNullOrWhiteSpace(SettingsFileName))
            {
                m_ErrMsg = "Settings file name not specified";
                return false;
            }

            // Parameter file name defined?
            if (string.IsNullOrWhiteSpace(ParameterFileName))
            {
                m_ErrMsg = "Parameter file name not specified";
                return false;
            }

            // Define the parameter file path; this is passed as the search tool parameter file
            m_ParameterFilePath = Path.Combine(SourceDirectoryPath, ParameterFileName);

            // Define the settings file path; this is passed as the parameter file name to m_PeptideHitResultsProcessor
            m_SettingsFilePath = Path.Combine(SourceDirectoryPath, SettingsFileName);

            // Define the peptide hit results format based on the analysis tool name
            if (AnalysisToolName.IndexOf(XtandemResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.XTandemXMLFile;
            }
            else if (AnalysisToolName.IndexOf(SequestResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.SequestSynopsisFile;
            }
            else if (AnalysisToolName.IndexOf(InSpecTResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.InspectTXTFile;
            }
            else if (AnalysisToolName.IndexOf("msgfdb", StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.MSGFPlusTXTFile;
            }
            else if (AnalysisToolName.IndexOf(MSGFPlusResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.MSGFPlusTXTFile;
            }
            else if (AnalysisToolName.IndexOf(MSAlignResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.MSAlignTXTFile;
            }
            else if (AnalysisToolName.IndexOf(MODaResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.MODaTXTFile;
            }
            else if (AnalysisToolName.IndexOf(MODPlusResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.MODPlusTXTFile;
            }
            else if (AnalysisToolName.IndexOf(MSPathFinderResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.MSPathFinderTSVFile;
            }
            else if (AnalysisToolName.IndexOf(TopPICResultsProcessor.TOOL_NAME, StringComparison.OrdinalIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.TopPICTXTFile;
            }
            else if (AnalysisToolName.IndexOf("DataExtractor", StringComparison.OrdinalIgnoreCase) >= 0)
            {
                // Data Extractor step-tool; we'll need to auto-determine the results format
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.AutoDetermine;
            }
            else
            {
                // Unrecognized analysis tool name
                m_PeptideHitResultsFileFormat = Enums.ResultsFileFormat.AutoDetermine;
            }

            if (DebugLevel >= 3)
            {
                OnDebugEvent("Setup params: AnalysisToolName = " + AnalysisToolName);
                OnDebugEvent("Setup params: PeptideHitResultsFileFormat = " + m_PeptideHitResultsFileFormat.ToString());

                OnDebugEvent("Setup params: DSName = " + DatasetName);
                OnDebugEvent("Setup params: SettingsFilePath = " + m_SettingsFilePath);
                OnDebugEvent("Setup params: ParameterFilePath = " + m_ParameterFilePath);
            }

            // Define the peptide hit results file name
            if (string.IsNullOrWhiteSpace(PeptideHitResultsFileName))
            {
                m_PeptideHitResultsFilePath = PHRPBaseClass.AutoDefinePeptideHitResultsFilePath(m_PeptideHitResultsFileFormat, SourceDirectoryPath, DatasetName);
            }
            else
            {
                m_PeptideHitResultsFilePath = Path.Combine(SourceDirectoryPath, PeptideHitResultsFileName);
            }

            if (DebugLevel >= 3)
            {
                OnDebugEvent("Setup params: PeptideHitResultsFilePath = " + m_PeptideHitResultsFilePath);
            }

            // Now that m_PeptideHitResultsFilePath has been determined, if m_PeptideHitResultsFileFormat is .AutoDetermine, try to determine the correct format
            if (m_PeptideHitResultsFileFormat == Enums.ResultsFileFormat.AutoDetermine)
            {
                m_PeptideHitResultsFileFormat = PHRPBaseClass.DetermineResultsFileFormat(m_PeptideHitResultsFilePath);
            }

            // Define the mass correction tags file path
            if (string.IsNullOrWhiteSpace(MassCorrectionTagsFileName))
            {
                m_MassCorrectionTagsFilePath = Path.Combine(SourceDirectoryPath, DEFAULT_MASS_CORRECTION_TAGS_FILENAME);
            }
            else
            {
                m_MassCorrectionTagsFilePath = Path.Combine(SourceDirectoryPath, MassCorrectionTagsFileName);
            }

            // Define the modification definitions file path
            if (string.IsNullOrWhiteSpace(ModificationDefinitionsFileName))
            {
                m_ModificationDefinitionsFilePath = Path.Combine(SourceDirectoryPath, Path.GetFileNameWithoutExtension(ParameterFileName) + MODIFICATION_DEFINITIONS_FILE_SUFFIX);
            }
            else
            {
                m_ModificationDefinitionsFilePath = Path.Combine(SourceDirectoryPath, ModificationDefinitionsFileName);
            }

            if (DebugLevel >= 3)
            {
                OnDebugEvent("Setup params: PeptideHitResultsFileFormat = " + m_PeptideHitResultsFileFormat.ToString());
                OnDebugEvent("Setup params: MassCorrectionTagsFilePath = " + m_MassCorrectionTagsFilePath);
                OnDebugEvent("Setup params: ModificationDefinitionsFilePath = " + m_ModificationDefinitionsFilePath);
            }

            // Parameter file exists?
            if (!VerifyFileExists(m_ParameterFilePath))
                return false; // Error msg handled by VerifyFileExists

            // Settings file exists?
            if (!VerifyFileExists(m_SettingsFilePath))
                return false; // Error msg handled by VerifyFileExists

            // Peptide hit results file exists?
            if (!VerifyFileExists(m_PeptideHitResultsFilePath))
                return false; // Error msg handled by VerifyFileExists

            // Modification definitions file exists?
            if (!VerifyFileExists(m_ModificationDefinitionsFilePath))
                return false; // Error msg handled by VerifyFileExists

            // Mass correction tags file exists?
            if (!VerifyFileExists(m_MassCorrectionTagsFilePath))
                return false; // Error msg handled by VerifyFileExists

            // If we got here, everything's OK
            return true;
        }

        private void LogErrors(string source, string message, Exception ex)
        {
            m_ErrMsg = string.Copy(message).Replace("\n", "; ");

            if (ex?.Message != null && ex.Message.Length > 0)
            {
                m_ErrMsg += "; " + ex.Message;
            }

            OnErrorEvent(source + ": " + m_ErrMsg);
        }

        /// <summary>
        /// Verifies that the specified directory exists
        /// </summary>
        /// <param name="directoryPath"></param>
        /// <returns>True if the directory exists, otherwise false</returns>
        /// <remarks>Updates m_ErrMsg if the directory is not found</remarks>
        protected virtual bool VerifyDirExists(string directoryPath)
        {
            if (Directory.Exists(directoryPath))
            {
                m_ErrMsg = string.Empty;
                return true;
            }

            m_ErrMsg = "Directory " + directoryPath + " not found";
            return false;
        }

        /// <summary>
        /// Verifies that the specified file exists
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns>True if the file exists, otherwise false</returns>
        /// <remarks>Updates m_ErrMsg if the file is not found</remarks>
        protected virtual bool VerifyFileExists(string filePath)
        {
            if (File.Exists(filePath))
            {
                m_ErrMsg = string.Empty;
                return true;
            }

            m_ErrMsg = "File " + filePath + " not found";
            return false;
        }

        private void mPeptideHitResultsProcessor_ProgressComplete()
        {
            // OperationComplete()
        }

        private void mPeptideHitResultsProcessor_ProgressReset()
        {
            // ResetProgress()
        }
    }
}
