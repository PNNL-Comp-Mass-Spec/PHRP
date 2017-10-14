// This class calls clsSequestSynopsisFileProcessor or clsXTandemResultsConverter
// to process the files to determine the modifications present for each peptide,
// along with other information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.

using System;
using System.IO;
using System.Reflection;
using PeptideHitResultsProcessor;
using PHRPReader;

namespace PeptideHitResultsProcRunner
{
    public class clsPeptideHitResultsProcRunner : clsProcessFilesBaseClass
    {
        public clsPeptideHitResultsProcRunner()
        {
            mFileDate = Program.PROGRAM_DATE;
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        // Error codes specialized for this class
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
        /// When this is True, then mCreateProteinModsFile is assumed to be true
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
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found then it will be created using the the Fasta file specified by mFastaFilePath</remarks>
        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        #endregion

        public override string[] GetDefaultExtensionsToParse()
        {
            var strExtensionsToParse = new string[2];

            // Note: If mPeptideHitResultsFileFormat = .AutoDetermine
            //  then this class will only parse .txt files if they match
            // PeptideHitResultsProcessor.clsSequestResultsProcessor.SEQUEST_FIRST_HITS_FILE_SUFFIX or
            // PeptideHitResultsProcessor.clsSequestResultsProcessor.SEQUEST_SYNOPSIS_FILE_SUFFIX

            strExtensionsToParse[0] = ".txt";
            strExtensionsToParse[1] = ".xml";

            return strExtensionsToParse;
        }

        public override string GetErrorMessage()
        {
            // Returns "" if no error

            string strErrorMessage;

            if (ErrorCode == eProcessFilesErrorCodes.LocalizedError |
                ErrorCode == eProcessFilesErrorCodes.NoError)
            {
                switch (mLocalErrorCode)
                {
                    case eResultsProcessorErrorCodes.NoError:
                        strErrorMessage = "";
                        break;
                    case eResultsProcessorErrorCodes.UnspecifiedError:
                        strErrorMessage = "Unspecified localized error";
                        break;
                    default:
                        // This shouldn't happen
                        strErrorMessage = "Unknown error state";
                        break;
                }
            }
            else
            {
                strErrorMessage = GetBaseClassErrorMessage();
            }

            return strErrorMessage;
        }

        private void InitializeLocalVariables()
        {
            ShowMessages = false;

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

            MsgfPlusEValueThreshold = clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MsgfPlusSpecEValueThreshold = clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            WarnMissingParameterFileSection = true;

            mLocalErrorCode = eResultsProcessorErrorCodes.NoError;
        }

        private void InitializePeptideHitResultsProcessor(string strInputFilePath)
        {
            if (mObtainModificationDefinitionsFromDMS)
            {
                LoadModificationInfoFromDMS();
            }

            var sourceFile = new FileInfo(strInputFilePath);

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

            mPeptideHitResultsProcessor.MSGFDBSynopsisFileEValueThreshold = MsgfPlusEValueThreshold;
            mPeptideHitResultsProcessor.MSGFDBSynopsisFileSpecEValueThreshold = MsgfPlusSpecEValueThreshold;

            mPeptideHitResultsProcessor.WarnMissingParameterFileSection = WarnMissingParameterFileSection;
        }

        private void LoadModificationInfoFromDMS()
        {
            // ToDo: Contact DMS to get the modification information
            // The results of this query will need to be filtered to get the info for just this analysis job

            //SELECT D.Dataset_Num,
            //    AJ.AJ_jobID, PFMI.Local_Symbol,
            //    PFMI.Monoisotopic_Mass_Correction, PFMI.Residue_Symbol,
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
            //    PFMI.Monoisotopic_Mass_Correction, PFMI.Residue_Symbol,
            //    PFMI.Mod_Type_Symbol, PFMI.Mass_Correction_Tag
            //FROM dbo.T_Analysis_Job AJ INNER JOIN
            //    dbo.T_Dataset D ON
            //    AJ.AJ_datasetID = D.Dataset_ID LEFT
            //     OUTER JOIN
            //    dbo.V_Param_File_Mass_Mod_Info PFMI ON
            //    AJ.AJ_parmFileName = PFMI.Param_File_Name
            //WHERE (AJ.AJ_jobID = 47703)
            //ORDER BY AJ.AJ_jobID, PFMI.Local_Symbol
        }

        private bool LoadParameterFileSettings(string strParameterFilePath)
        {
            const string OPTIONS_SECTION = "PeptideHitResultsProcRunner";

            var objSettingsFile = new PRISM.XmlSettingsFileAccessor();

            try
            {
                if (string.IsNullOrWhiteSpace(strParameterFilePath))
                {
                    // No parameter file specified; nothing to load
                    return true;
                }

                if (!File.Exists(strParameterFilePath))
                {
                    // See if strParameterFilePath points to a file in the same directory as the application
                    strParameterFilePath = Path.Combine(Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location), Path.GetFileName(strParameterFilePath));
                    if (!File.Exists(strParameterFilePath))
                    {
                        SetBaseClassErrorCode(eProcessFilesErrorCodes.ParameterFileNotFound);
                        return false;
                    }
                }

                if (objSettingsFile.LoadSettings(strParameterFilePath))
                {
                    if (!objSettingsFile.SectionPresent(OPTIONS_SECTION))
                    {
                        ShowErrorMessage("The node '<section name=\"" + OPTIONS_SECTION + "\"> was not found in the parameter file: " + strParameterFilePath);
                        SetBaseClassErrorCode(eProcessFilesErrorCodes.InvalidParameterFile);
                        return false;
                    }

                    mObtainModificationDefinitionsFromDMS = objSettingsFile.GetParam(OPTIONS_SECTION, "ObtainModificationDefinitionsFromDMS", mObtainModificationDefinitionsFromDMS);

                    var intValue = objSettingsFile.GetParam(OPTIONS_SECTION, "PeptideHitResultsFileFormat", Convert.ToInt32(mPeptideHitResultsFileFormat));
                    try
                    {
                        mPeptideHitResultsFileFormat = (clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants) intValue;
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

        // Main processing function
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath, bool blnResetErrorCode)
        {
            // Returns True if success, False if failure

            var blnSuccess = false;

            if (blnResetErrorCode)
            {
                SetLocalErrorCode(eResultsProcessorErrorCodes.NoError);
            }

            if (!LoadParameterFileSettings(strParameterFilePath))
            {
                var strStatusMessage = "Parameter file load error: " + strParameterFilePath;
                ShowErrorMessage(strStatusMessage);

                if (ErrorCode == eProcessFilesErrorCodes.NoError)
                {
                    SetBaseClassErrorCode(eProcessFilesErrorCodes.InvalidParameterFile);
                }
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(strInputFilePath))
                {
                    ShowErrorMessage("Input file name is empty");
                    SetBaseClassErrorCode(eProcessFilesErrorCodes.InvalidInputFilePath);
                }
                else
                {
                    // Note that CleanupFilePaths() will update mOutputFolderPath, which is used by LogMessage()
                    if (!CleanupFilePaths(ref strInputFilePath, ref strOutputFolderPath))
                    {
                        SetBaseClassErrorCode(eProcessFilesErrorCodes.FilePathError);
                    }
                    else
                    {
                        mProgressStepDescription = "Parsing " + Path.GetFileName(strInputFilePath);
                        LogMessage(mProgressStepDescription);
                        ResetProgress();

                        if (CreateProteinModsUsingPHRPDataFile)
                        {
                            blnSuccess = StartCreateProteinModsViaPHRPData(strInputFilePath, strOutputFolderPath);
                        }
                        else
                        {
                            blnSuccess = StartPHRP(strInputFilePath, strOutputFolderPath, strParameterFilePath);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFile", ex);
            }

            return blnSuccess;
        }

        /// <summary>
        /// Looks for file strFileNameOrPath in the current working directory
        /// If not found, then looks in strSourceFolderPath
        /// </summary>
        /// <param name="strSourceFolderPath">Path to the folder containing the input file</param>
        /// <param name="strFileNameOrPath">File to find (either filename or full file path)</param>
        /// <returns>The path to the file if found, or strFileNameOrPath if not found</returns>
        /// <remarks></remarks>
        protected string ResolveFilePath(string strSourceFolderPath, string strFileNameOrPath)
        {
            if (File.Exists(strFileNameOrPath))
            {
                return strFileNameOrPath;
            }

            var strNewPath = Path.Combine(strSourceFolderPath, Path.GetFileName(strFileNameOrPath));
            if (File.Exists(strNewPath))
            {
                return strNewPath;
            }

            return strFileNameOrPath;
        }

        private void SetLocalErrorCode(eResultsProcessorErrorCodes eNewErrorCode, bool blnLeaveExistingErrorCodeUnchanged = false)
        {
            if (blnLeaveExistingErrorCodeUnchanged && mLocalErrorCode != eResultsProcessorErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mLocalErrorCode = eNewErrorCode;

                if (eNewErrorCode == eResultsProcessorErrorCodes.NoError)
                {
                    if (ErrorCode == eProcessFilesErrorCodes.LocalizedError)
                    {
                        SetBaseClassErrorCode(eProcessFilesErrorCodes.NoError);
                    }
                }
                else
                {
                    SetBaseClassErrorCode(eProcessFilesErrorCodes.LocalizedError);
                }
            }
        }

        private bool StartCreateProteinModsViaPHRPData(string strInputFilePath, string strOutputFolderPath)
        {
            var blnSuccess = false;

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
                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InSpectTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Inspect;
                        LogMessage("Detected Inspect results file");

                        break;
                    case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile:
                        ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MSGFDB;
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

                    default:
                        // Includes ePeptideHitResultsFileFormatConstants.AutoDetermine
                        ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath);
                        break;
                }

                string strMessage;
                if (ePeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
                {
                    blnSuccess = false;

                    strMessage = "Error: Could not determine the format of the PHRP data file: " + strInputFilePath;
                    ShowErrorMessage(strMessage);
                }
                else
                {
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

                        case clsPHRPReader.ePeptideHitResultType.MSGFDB:
                            mPeptideHitResultsProcessor = new clsMSGFDBResultsProcessor();
                            break;

                        case clsPHRPReader.ePeptideHitResultType.MSAlign:
                            mPeptideHitResultsProcessor = new clsMSAlignResultsProcessor();
                            break;

                        case clsPHRPReader.ePeptideHitResultType.MODa:
                            mPeptideHitResultsProcessor = new clsMODaResultsProcessor();
                            break;

                        default:
                            // Unknown format
                            strMessage = "Error: Unrecognized value for ePeptideHitResultType: " + ePeptideHitResultType.ToString();
                            ShowErrorMessage(strMessage);
                            blnSuccess = false;
                            break;
                    }

                    // Do not call RegisterEvents
                    // Instead use local event handlers that optionally log to a file
                    mPeptideHitResultsProcessor.ErrorEvent += mPeptideHitResultsProcessor_ErrorOccurred;
                    mPeptideHitResultsProcessor.StatusEvent += mPeptideHitResultsProcessor_MessageEvent;
                    mPeptideHitResultsProcessor.ProgressUpdate += mPeptideHitResultsProcessor_ProgressChanged;
                    mPeptideHitResultsProcessor.WarningEvent += mPeptideHitResultsProcessor_WarningMessageEvent;

                    if (mPeptideHitResultsProcessor != null)
                    {
                        InitializePeptideHitResultsProcessor(strInputFilePath);

                        blnSuccess = mPeptideHitResultsProcessor.CreateProteinModDetailsFile(strInputFilePath, strOutputFolderPath);

                        if (!blnSuccess)
                        {
                            ShowErrorMessage(mPeptideHitResultsProcessor.ErrorMessage);
                        }
                        else
                        {
                            LogMessage("Processing Complete");
                            OperationComplete();
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error calling CreateProteinModDetailsFile in CreateProteinModsViaPHRPData", ex);
            }

            return blnSuccess;
        }

        private bool StartPHRP(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            var blnSuccess = false;

            try
            {
                clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants ePeptideHitResultsFormat;
                if (mPeptideHitResultsFileFormat == clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine)
                {
                    ePeptideHitResultsFormat = clsPHRPBaseClass.DetermineResultsFileFormat(strInputFilePath);
                }
                else
                {
                    ePeptideHitResultsFormat = mPeptideHitResultsFileFormat;
                }

                if (ePeptideHitResultsFormat == clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.AutoDetermine)
                {
                    // If ePeptideHitResultsFormat is still AutoDetermine that means we couldn't figure out the format
                    blnSuccess = false;

                    var strMessage = "Error: Could not determine the format of the input file.  It must end in " +
                                        clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE + ".txt, " +
                                        clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE + ".txt, .xml (for X!Tandem), " +
                                        clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE + ".txt, " +
                                        clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE + ".tsv, " +
                                        clsMSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE + ".txt, " +
                                        clsMODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE + ".txt, " +
                                        clsMODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE + ".txt, or" +
                                        clsMSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE + ".tsv";

                    ShowErrorMessage(strMessage);
                }
                else
                {
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

                        case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.InSpectTXTFile:
                            mPeptideHitResultsProcessor = new clsInSpecTResultsProcessor();
                            LogMessage("Detected Inspect results file");
                            break;

                        case clsPHRPBaseClass.ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile:
                            mPeptideHitResultsProcessor = new clsMSGFDBResultsProcessor();
                            LogMessage("Detected MSGFDB (or MSGF+) results file");
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
                        default:
                            // Unknown format
                            blnSuccess = false;
                            break;
                    }

                    if ((mPeptideHitResultsProcessor != null))
                    {
                        InitializePeptideHitResultsProcessor(strInputFilePath);

                        blnSuccess = mPeptideHitResultsProcessor.ProcessFile(strInputFilePath, strOutputFolderPath, strParameterFilePath);
                        if (!blnSuccess)
                        {
                            ShowErrorMessage(mPeptideHitResultsProcessor.ErrorMessage);
                        }
                        else
                        {
                            LogMessage("Processing Complete");
                            OperationComplete();
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error calling ProcessFile in StartPHRP", ex);
            }

            return blnSuccess;
        }

        private void mPeptideHitResultsProcessor_ErrorOccurred(string message, Exception ex)
        {
            if (ex != null && !message.Contains(ex.Message))
            {
                LogMessage(message + ": " + ex.Message, eMessageTypeConstants.ErrorMsg);
            }
            else
            {
                LogMessage(message, eMessageTypeConstants.ErrorMsg);
            }

        }

        private void mPeptideHitResultsProcessor_MessageEvent(string message)
        {
            LogMessage(message);
        }

        private void mPeptideHitResultsProcessor_ProgressChanged(string taskDescription, float percentComplete)
        {
            UpdateProgress(taskDescription, percentComplete);
        }

        private void mPeptideHitResultsProcessor_WarningMessageEvent(string WarningMessage)
        {
            LogMessage(WarningMessage, eMessageTypeConstants.Warning);
        }
    }
}
