// This class implements the IPeptideHitResultsProcessor interface
// It uses class clsXTandemResultsProcessor to read an XTandem search results XML file
//  or class clsSequestResultsProcessor to read a Sequest Synopsis or First Hits file
//
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute
// Started January 6, 2006

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Threading;

namespace PeptideHitResultsProcessor
{
    public class clsAnalysisManagerPeptideHitResultsProcessor : IPeptideHitResultsProcessor
    {
        #region "Constants and enums"
        private const string DEFAULT_MASS_CORRECTION_TAGS_FILENAME = "Mass_Correction_Tags.txt";
        private const string MODIFICATION_DEFINITIONS_FILE_SUFFIX = "_ModDefs.txt";
        #endregion

        #region "Classwide variables"

        private clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine;

        [Obsolete("Unused")]
        private Dictionary<string, string> m_MiscParams;

        private string m_ParameterFilePath = string.Empty;      // Peptide search tool parameter file path
        private string m_SettingsFilePath = string.Empty;       // XML settings file with section PeptideHitResultsProcessorOptions

        private string m_PeptideHitResultsFilePath = string.Empty;
        private string m_MassCorrectionTagsFilePath = string.Empty;
        private string m_ModificationDefinitionsFilePath = string.Empty;

        private string m_ErrMsg = string.Empty;

        private IPeptideHitResultsProcessor.ProcessStatus m_Status;
        private IPeptideHitResultsProcessor.ProcessResults m_Results;

        private clsPHRPBaseClass withEventsField_m_PeptideHitResultsProcessor;
        private clsPHRPBaseClass m_PeptideHitResultsProcessor
        {
            get { return withEventsField_m_PeptideHitResultsProcessor; }
            set
            {
                if (withEventsField_m_PeptideHitResultsProcessor != null)
                {
                    withEventsField_m_PeptideHitResultsProcessor.ErrorOccurred -= m_PeptideHitResultsProcessor_ErrorOccurred;
                    withEventsField_m_PeptideHitResultsProcessor.ProgressChanged -= mPeptideHitResultsProcessor_ProgressChanged;
                    withEventsField_m_PeptideHitResultsProcessor.ProgressComplete -= mPeptideHitResultsProcessor_ProgressComplete;
                    withEventsField_m_PeptideHitResultsProcessor.ProgressReset -= mPeptideHitResultsProcessor_ProgressReset;
                }
                withEventsField_m_PeptideHitResultsProcessor = value;
                if (withEventsField_m_PeptideHitResultsProcessor != null)
                {
                    withEventsField_m_PeptideHitResultsProcessor.ErrorOccurred += m_PeptideHitResultsProcessor_ErrorOccurred;
                    withEventsField_m_PeptideHitResultsProcessor.ProgressChanged += mPeptideHitResultsProcessor_ProgressChanged;
                    withEventsField_m_PeptideHitResultsProcessor.ProgressComplete += mPeptideHitResultsProcessor_ProgressComplete;
                    withEventsField_m_PeptideHitResultsProcessor.ProgressReset += mPeptideHitResultsProcessor_ProgressReset;
                }
            }
        }

        private Thread m_thThread;

        #endregion

        #region "Events"
        public override event ErrorOccurredEventHandler ErrorOccurred;
        public override event DebugEventEventHandler DebugEvent;

        // PercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public override event ProgressChangedEventHandler ProgressChanged;

        #endregion

        #region "Properties"
        public override string AnalysisToolName { get; set; }

        public override bool CreateInspectFirstHitsFile { get; set; }

        public override bool CreateInspectSynopsisFile { get; set; }

        public override string DatasetName { get; set; }

        public override int DebugLevel { get; set; }

        public override string ErrMsg
        {
            get { return m_ErrMsg; }
        }

        public override string MassCorrectionTagsFileName { get; set; }

        [Obsolete("Unused")]
        public override Dictionary<string, string> MiscParams
        {
            get { return m_MiscParams; }
            set { m_MiscParams = value; }
        }

        public override string ModificationDefinitionsFileName { get; set; }

        public override string OutputFolderPath { get; set; }

        /// <summary>
        /// Peptide search tool parameter file name
        /// </summary>
        /// <returns></returns>
        public override string ParameterFileName { get; set; }

        public override string PeptideHitResultsFileName { get; set; }

        public override float PercentComplete
        {
            get
            {
                if ((m_PeptideHitResultsProcessor != null))
                {
                    return m_PeptideHitResultsProcessor.ProgressPercentComplete;
                }
                return 0;
            }
        }

        /// <summary>
        /// XML settings file with section PeptideHitResultsProcessorOptions
        /// </summary>
        /// <returns></returns>
        public override string SettingsFileName { get; set; }

        public override string SourceFolderPath { get; set; }

        public override IPeptideHitResultsProcessor.ProcessStatus Status
        {
            get { return m_Status; }
        }

        public override IPeptideHitResultsProcessor.ProcessResults Results
        {
            get { return m_Results; }
        }
        #endregion

        public override IPeptideHitResultsProcessor.ProcessStatus Abort()
        {
            if ((m_PeptideHitResultsProcessor != null))
            {
                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ABORTING;
                m_PeptideHitResultsProcessor.AbortProcessingNow();
            }
            return m_Status;
        }

        public override void Setup(IPeptideHitResultsProcessor.InitializationParams InitParams)
        {
            //Copies all input data required for plugin operation to appropriate memory variables
            SourceFolderPath = InitParams.SourceFolderPath;
            OutputFolderPath = InitParams.OutputFolderPath;

            PeptideHitResultsFileName = InitParams.PeptideHitResultsFileName;
            MassCorrectionTagsFileName = InitParams.MassCorrectionTagsFileName;
            ModificationDefinitionsFileName = InitParams.ModificationDefinitionsFileName;

            // This is unused and thus obsolete
            // m_MiscParams = .MiscParams

            DebugLevel = InitParams.DebugLevel;

            AnalysisToolName = InitParams.AnalysisToolName;
            DatasetName = InitParams.DatasetName;
            ParameterFileName = InitParams.ParameterFileName;
            SettingsFileName = InitParams.SettingsFileName;

            CreateInspectFirstHitsFile = InitParams.CreateInspectFirstHitsFile;
            CreateInspectSynopsisFile = InitParams.CreateInspectSynopsisFile;
        }

        public override IPeptideHitResultsProcessor.ProcessStatus Start()
        {
            m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_STARTING;

            //Verify necessary files are in specified locations
            if (!InitSetup())
            {
                m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_FAILURE;
                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR;
                return m_Status;
            }

            // Process the results file (the process runs in a separate thread)
            m_Status = ProcessPeptideHitResultsFile();

            if (m_Status == IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR)
            {
                m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_FAILURE;
            }

            return m_Status;
        }

        /// <summary>
        /// Initializes m_PeptideHitResultsProcessor then starts a separate thread to process the file
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        protected virtual IPeptideHitResultsProcessor.ProcessStatus ProcessPeptideHitResultsFile()
        {
            try
            {
                //Initialize m_PeptideHitResultsProcessor
                switch (m_PeptideHitResultsFileFormat)
                {
                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.XTandemXMLFile:
                        m_PeptideHitResultsProcessor = new clsXTandemResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile:
                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestSynopsisFile:
                        m_PeptideHitResultsProcessor = new clsSequestResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InSpectTXTFile:
                        m_PeptideHitResultsProcessor = new clsInSpecTResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile:
                        m_PeptideHitResultsProcessor = new clsMSGFDBResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSAlignTXTFile:
                        m_PeptideHitResultsProcessor = new clsMSAlignResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODaTXTFile:
                        m_PeptideHitResultsProcessor = new clsMODaResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODPlusTXTFile:
                        m_PeptideHitResultsProcessor = new clsMODaResultsProcessor();
                        break;

                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile:
                        m_PeptideHitResultsProcessor = new clsMSPathFinderResultsProcessor();
                        break;

                    default:
                        // Unknown format; cannot continue
                        LogErrors("ProcessPeptideHitResultsFile", "Unknown peptide hit results file format: " + m_PeptideHitResultsFileFormat.ToString(), null);
                        m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR;
                        return m_Status;
                }

                // Define the auxiliary file paths
                m_PeptideHitResultsProcessor.MassCorrectionTagsFilePath = m_MassCorrectionTagsFilePath;
                m_PeptideHitResultsProcessor.ModificationDefinitionsFilePath = m_ModificationDefinitionsFilePath;
                m_PeptideHitResultsProcessor.SearchToolParameterFilePath = m_ParameterFilePath;
                m_PeptideHitResultsProcessor.InspectSynopsisFilePValueThreshold = clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

                m_PeptideHitResultsProcessor.CreateInspectFirstHitsFile = CreateInspectFirstHitsFile;
                m_PeptideHitResultsProcessor.CreateInspectSynopsisFile = CreateInspectSynopsisFile;

                m_thThread = new Thread(ProcessPeptideHitResultsFileWork);
                m_thThread.Start();

                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_RUNNING;
            }
            catch (Exception ex)
            {
                LogErrors("ProcessPeptideHitResultsFile", "Error initializing and running m_PeptideHitResultsProcessor", ex);
                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR;
            }

            return m_Status;
        }

        protected virtual void ProcessPeptideHitResultsFileWork()
        {
            bool blnSuccess = false;

            try
            {
                blnSuccess = m_PeptideHitResultsProcessor.ProcessFile(m_PeptideHitResultsFilePath, OutputFolderPath, m_SettingsFilePath);
            }
            catch (Exception)
            {
                blnSuccess = false;
            }

            if (blnSuccess)
            {
                m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_COMPLETE;
            }
            else
            {
                if (m_PeptideHitResultsProcessor.AbortProcessing)
                {
                    LogErrors("ProcessPeptideHitResultsFileWork", "Processing aborted", null);
                    m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_ABORTED;
                    m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ABORTING;
                }
                else
                {
                    LogErrors("ProcessPeptideHitResultsFileWork", m_PeptideHitResultsProcessor.ErrorMessage, null);
                    m_Results = IPeptideHitResultsProcessor.ProcessResults.PH_FAILURE;
                    m_Status = IPeptideHitResultsProcessor.ProcessStatus.PH_ERROR;
                }
            }
        }

        protected virtual bool InitSetup()
        {
            //Initializes module variables and verifies mandatory parameters have been propery specified

            //Output folder name
            if (string.IsNullOrEmpty(OutputFolderPath))
            {
                m_ErrMsg = "Output folder path not specified";
                return false;
            }

            //Source folder name
            if (string.IsNullOrEmpty(SourceFolderPath))
            {
                m_ErrMsg = "Source folder not specified";
                return false;
            }

            if (this.DebugLevel >= 3)
            {
                if (DebugEvent != null)
                {
                    DebugEvent("Setup params: OutFolderPath = " + OutputFolderPath);
                }
                if (DebugEvent != null)
                {
                    DebugEvent("Setup params: SourceFolderPath = " + SourceFolderPath);
                }
            }

            //Source directory exists?
            if (!VerifyDirExists(SourceFolderPath))
                return false; //Error msg handled by VerifyDirExists

            //Output directory exists?
            if (!VerifyDirExists(OutputFolderPath))
                return false; //Error msg handled by VerifyDirExists

            //Analysis tool name defined?
            if (string.IsNullOrWhiteSpace(AnalysisToolName))
            {
                m_ErrMsg = "Analysis tool name not specified";
                return false;
            }

            //Dataset name defined?
            if (string.IsNullOrWhiteSpace(DatasetName))
            {
                m_ErrMsg = "Dataset name not specified";
                return false;
            }

            //Settings file name defined?
            if (string.IsNullOrWhiteSpace(SettingsFileName))
            {
                m_ErrMsg = "Settings file name not specified";
                return false;
            }

            //Parameter file name defined?
            if (string.IsNullOrWhiteSpace(ParameterFileName))
            {
                m_ErrMsg = "Parameter file name not specified";
                return false;
            }

            //Define the parameter file path; this is passed as the search tool parameter file
            m_ParameterFilePath = Path.Combine(SourceFolderPath, ParameterFileName);

            //Define the settings file path; this is passed as the parameter file name to m_PeptideHitResultsProcessor
            m_SettingsFilePath = Path.Combine(SourceFolderPath, SettingsFileName);

            //Define the peptide hit results format based on the analysis tool name
            if (AnalysisToolName.IndexOf("xtandem", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.XTandemXMLFile;
            }
            else if (AnalysisToolName.IndexOf("sequest", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.SequestSynopsisFile;
            }
            else if (AnalysisToolName.IndexOf("inspect", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InSpectTXTFile;
            }
            else if (AnalysisToolName.IndexOf("msgfdb", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile;
            }
            else if (AnalysisToolName.IndexOf("msalign", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSAlignTXTFile;
            }
            else if (AnalysisToolName.IndexOf("moda", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODaTXTFile;
            }
            else if (AnalysisToolName.IndexOf("modplus", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MODPlusTXTFile;
            }
            else if (AnalysisToolName.IndexOf("mspathfinder", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile;
            }
            else if (AnalysisToolName.IndexOf("dataextractor", StringComparison.InvariantCultureIgnoreCase) >= 0)
            {
                // Data Extractor step-tool; we'll need to auto-determine the results format
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine;
            }
            else
            {
                // Unrecognized analysis tool name
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine;
            }

            if (this.DebugLevel >= 3)
            {
                if (DebugEvent != null)
                {
                    DebugEvent("Setup params: AnalysisToolName = " + AnalysisToolName);
                    DebugEvent("Setup params: PeptideHitResultsFileFormat = " + m_PeptideHitResultsFileFormat.ToString());

                    DebugEvent("Setup params: DSName = " + DatasetName);
                    DebugEvent("Setup params: SettingsFilePath = " + m_SettingsFilePath);
                    DebugEvent("Setup params: ParameterFilePath = " + m_ParameterFilePath);
                }
            }

            //Define the peptide hit results file name
            if (string.IsNullOrWhiteSpace(PeptideHitResultsFileName))
            {
                m_PeptideHitResultsFilePath = clsPHRPBaseClass.AutoDefinePeptideHitResultsFilePath(m_PeptideHitResultsFileFormat, SourceFolderPath, DatasetName);
            }
            else
            {
                m_PeptideHitResultsFilePath = Path.Combine(SourceFolderPath, PeptideHitResultsFileName);
            }

            if (this.DebugLevel >= 3)
            {
                if (DebugEvent != null)
                {
                    DebugEvent("Setup params: PeptideHitResultsFilePath = " + m_PeptideHitResultsFilePath);
                }
            }

            //Now that m_PeptideHitResultsFilePath has been determined, if m_PeptideHitResultsFileFormat is .AutoDetermine then try to determine the correct format
            if (m_PeptideHitResultsFileFormat == clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine)
            {
                m_PeptideHitResultsFileFormat = clsPHRPBaseClass.DetermineResultsFileFormat(m_PeptideHitResultsFilePath);
            }

            //Define the mass correction tags file path
            if (string.IsNullOrWhiteSpace(MassCorrectionTagsFileName))
            {
                m_MassCorrectionTagsFilePath = Path.Combine(SourceFolderPath, DEFAULT_MASS_CORRECTION_TAGS_FILENAME);
            }
            else
            {
                m_MassCorrectionTagsFilePath = Path.Combine(SourceFolderPath, MassCorrectionTagsFileName);
            }

            //Define the modification definitions file path
            if (string.IsNullOrWhiteSpace(ModificationDefinitionsFileName))
            {
                m_ModificationDefinitionsFilePath = Path.Combine(SourceFolderPath, Path.GetFileNameWithoutExtension(ParameterFileName) + MODIFICATION_DEFINITIONS_FILE_SUFFIX);
            }
            else
            {
                m_ModificationDefinitionsFilePath = Path.Combine(SourceFolderPath, ModificationDefinitionsFileName);
            }

            if (this.DebugLevel >= 3)
            {
                if (DebugEvent != null)
                {
                    DebugEvent("Setup params: PeptideHitResultsFileFormat = " + m_PeptideHitResultsFileFormat.ToString());
                    DebugEvent("Setup params: MassCorrectionTagsFilePath = " + m_MassCorrectionTagsFilePath);
                    DebugEvent("Setup params: ModificationDefinitionsFilePath = " + m_ModificationDefinitionsFilePath);
                }
            }

            //Parameter file exists?
            if (!VerifyFileExists(m_ParameterFilePath))
                return false; //Error msg handled by VerifyFileExists

            //Settings file exists?
            if (!VerifyFileExists(m_SettingsFilePath))
                return false; //Error msg handled by VerifyFileExists

            //Peptide hit results file exists?
            if (!VerifyFileExists(m_PeptideHitResultsFilePath))
                return false; //Error msg handled by VerifyFileExists

            //Modification definitions file exists?
            if (!VerifyFileExists(m_ModificationDefinitionsFilePath))
                return false; //Error msg handled by VerifyFileExists

            //Mass correction tags file exists?
            if (!VerifyFileExists(m_MassCorrectionTagsFilePath))
                return false; //Error msg handled by VerifyFileExists

            //If we got here, everything's OK
            return true;
        }

        private void LogErrors(string strSource, string strMessage, Exception ex, bool blnLogLocalOnly = true)
        {
            m_ErrMsg = string.Copy(strMessage).Replace("\n", "; ");

            if (ex != null)
            {
                if ((ex.Message != null) && ex.Message.Length > 0)
                {
                    m_ErrMsg += "; " + ex.Message;
                }
            }

            Trace.WriteLine(System.DateTime.Now.ToLongTimeString() + "; " + m_ErrMsg, strSource);
            Console.WriteLine(System.DateTime.Now.ToLongTimeString() + "; " + m_ErrMsg, strSource);

            if (ErrorOccurred != null)
            {
                ErrorOccurred(m_ErrMsg);
            }
        }

        private string strProgressStepDescriptionSaved = string.Empty;
        private float sngProgressPercentComplete = 0;
        private void UpdateProgress(string strProgressStepDescription, float sngPercentComplete)
        {
            var blnDescriptionChanged = false;

            if (strProgressStepDescription != strProgressStepDescriptionSaved)
            {
                blnDescriptionChanged = true;
            }

            strProgressStepDescriptionSaved = string.Copy(strProgressStepDescription);
            if (sngPercentComplete < 0)
            {
                sngPercentComplete = 0;
            }
            else if (sngPercentComplete > 100)
            {
                sngPercentComplete = 100;
            }
            sngProgressPercentComplete = sngPercentComplete;

            if (blnDescriptionChanged & this.DebugLevel >= 2)
            {
                if (Math.Abs(sngProgressPercentComplete) < float.Epsilon)
                {
                    if (DebugEvent != null)
                    {
                        DebugEvent(strProgressStepDescriptionSaved);
                    }
                }
                else
                {
                    if (DebugEvent != null)
                    {
                        DebugEvent(strProgressStepDescriptionSaved + " (" + sngProgressPercentComplete.ToString("0.0") + "% complete)");
                    }
                }
            }

            if (ProgressChanged != null)
            {
                ProgressChanged(strProgressStepDescription, sngPercentComplete);
            }
        }

        protected virtual bool VerifyDirExists(string TestDir)
        {
            //Verifies that the specified directory exists
            if (Directory.Exists(TestDir))
            {
                m_ErrMsg = "";
                return true;
            }
            else
            {
                m_ErrMsg = "Directory " + TestDir + " not found";
                return false;
            }
        }

        protected virtual bool VerifyFileExists(string TestFile)
        {
            //Verifies specified file exists
            if (File.Exists(TestFile))
            {
                m_ErrMsg = "";
                return true;
            }
            else
            {
                m_ErrMsg = "File " + TestFile + " not found";
                return false;
            }
        }

        private void m_PeptideHitResultsProcessor_ErrorOccurred(string ErrorMessage)
        {
            LogErrors("PeptideHitResultsProcessor", ErrorMessage, null, true);
        }

        private void mPeptideHitResultsProcessor_ProgressChanged(string taskDescription, float pctComplete)
        {
            UpdateProgress(taskDescription, pctComplete);
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
