// This class can be used as a base class for peptide hit results processor classes
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
// Started January 6, 2006
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
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Threading;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public abstract class clsPHRPBaseClass : PRISM.clsEventNotifier
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        protected clsPHRPBaseClass()
        {
            mFileDate = "October 13, 2017";

            mPeptideSeqMassCalculator = new clsPeptideMassCalculator {ChargeCarrierMass = clsPeptideMassCalculator.MASS_PROTON};

            // Initialize mPeptideMods
            mPeptideMods = new clsPeptideModificationContainer();

            // Initialize mUniqueSequences
            mUniqueSequences = new clsUniqueSequencesContainer();

            // Initialize mSeqToProteinMap
            mSeqToProteinMap = new Hashtable();

            // Define a RegEx to replace all of the non-letter characters
            mReplaceSymbols = new Regex(@"[^A-Za-z]", RegexOptions.Compiled);

            mProteinNameOrder = new Dictionary<string, int>();

            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        protected const string SEP_CHAR = "\t";
        private const string UNIQUE_SEQ_TO_PROTEIN_MAP_SEP = "_";

        private const string COLUMN_NAME_UNIQUE_SEQ_ID = "Unique_Seq_ID";
        private const string COLUMN_NAME_PROTEIN_NAME = "Protein_Name";
        protected const string COLUMN_NAME_RESULTID = "ResultID";
        protected const string COLUMN_NAME_PEPTIDE = "Peptide";
        private const string COLUMN_NAME_RESIDUE = "Residue";
        private const string COLUMN_NAME_PROTEIN_RESIDUE_NUMBER = "Protein_Residue_Num";
        private const string COLUMN_NAME_RESIDUE_MOD_NAME = "Mod_Name";
        private const string COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER = "Peptide_Residue_Num";
        private const string COLUMN_NAME_MSGF_SPECPROB = "MSGF_SpecProb";

        public const string XTANDEM_RESULTS_FILE_SUFFIX = "_xt.xml";
        public const string SEQUEST_SYNOPSIS_FILE_SUFFIX = "_syn.txt";
        public const string SEQUEST_FIRST_HITS_FILE_SUFFIX = "_fht.txt";

        public const string INSPECT_RESULTS_FILE_SUFFIX = "_inspect.txt";
        public const string INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX = "_fht.txt";
        public const string INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX = "_Fscore_fht.txt";

        public const string MSGFDB_RESULTS_FILE_SUFFIX = "_msgfdb.txt";

        public const string MSALIGN_RESULTS_FILE_SUFFIX = "_MSAlign_ResultTable.txt";

        public const string MODa_RESULTS_FILE_SUFFIX = "_moda.id.txt";
        public const string MODPlus_RESULTS_FILE_SUFFIX = "_modp.id.txt";
        public const string MSPathFinder_RESULTS_FILE_SUFFIX = "_IcTda.tsv";

        public const string FILENAME_SUFFIX_RESULT_TO_SEQ_MAP = "_ResultToSeqMap.txt";
        public const string FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP = "_SeqToProteinMap.txt";

        public const string FILENAME_SUFFIX_SEQ_INFO = "_SeqInfo.txt";
        public const string FILENAME_SUFFIX_MOD_DETAILS = "_ModDetails.txt";
        public const string FILENAME_SUFFIX_MOD_SUMMARY = "_ModSummary.txt";

        public const string FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING = "_PepToProtMap";
        public const string FILENAME_SUFFIX_PROTEIN_MODS = "_ProteinMods.txt";
        public const string FILENAME_SUFFIX_MSGF = "_MSGF.txt";

        protected const float PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE = 90;
        private const float PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE = 95;

        private const string PROTEIN_NAME_NO_MATCH  = "__NoMatch__";

        public enum ePeptideHitResultsFileFormatConstants
        {
            AutoDetermine = 0,
            SequestSynopsisFile = 1,
            SequestFirstHitsFile = 2,
            XTandemXMLFile = 3,
            InSpectTXTFile = 4,
            MSGFDbTXTFile = 5,
            MSAlignTXTFile = 6,
            MODaTXTFile = 7,
            MODPlusTXTFile = 8,
            MSPathFinderTSVFile = 9
        }

        public enum ePHRPErrorCodes
        {
            NoError = 0,
            InvalidInputFilePath = 1,
            InvalidOutputFolderPath = 2,
            ParameterFileNotFound = 3,
            MassCorrectionTagsFileNotFound = 4,
            ModificationDefinitionFileNotFound = 5,

            ErrorReadingInputFile = 6,
            ErrorCreatingOutputFiles = 7,
            ErrorReadingParameterFile = 8,
            ErrorReadingMassCorrectionTagsFile = 9,
            ErrorReadingModificationDefinitionsFile = 10,

            FilePathError = 11,
            UnspecifiedError = -1
        }

        #endregion

        #region "Structures"
        protected struct udtSearchOptionModificationInfoType
        {
            public int SortOrder;
            public double ModificationMass;
            public string TargetResidues;
            public clsModificationDefinition.eModificationTypeConstants ModificationType;

            public override string ToString()
            {
                return ModificationType + ": " + ModificationMass + " @ " + TargetResidues;
            }
        }

        internal struct udtModNameAndResidueLocType
        {
            public string ModName;
            public int ResidueLocInPeptide;

            public override string ToString()
            {
                return ResidueLocInPeptide + ": " + ModName;
            }
        }

        protected struct udtPepToProteinMappingType
        {
            public string Peptide;
            public string Protein;
            public int ResidueStart;
            public int ResidueEnd;

            public override string ToString()
            {
                return Peptide + ", " + Protein;
            }
        }

        #endregion

        #region "Classwide Variables"

        protected ePHRPErrorCodes mErrorCode = ePHRPErrorCodes.NoError;
        protected string mErrorMessage = string.Empty;

        protected string mFileDate;

        protected readonly clsPeptideMassCalculator mPeptideSeqMassCalculator;

        protected readonly clsPeptideModificationContainer mPeptideMods;
        private readonly clsUniqueSequencesContainer mUniqueSequences;
        private readonly Hashtable mSeqToProteinMap;

        private StreamWriter mResultToSeqMapFile;
        private StreamWriter mSeqInfoFile;
        private StreamWriter mModDetailsFile;
        private StreamWriter mSeqToProteinMapFile;

        private int mNextPeptideToProteinMapperLevel;

        /// <summary>
        /// Tracks the protein names in the order that they are listed in the FASTA file
        /// Keys are protein Names, values are a sequentially assigned integer
        /// </summary>
        protected readonly Dictionary<string, int> mProteinNameOrder;

        private readonly Regex mReplaceSymbols;

        #endregion

        #region "Progress Events and Variables"

        public event ProgressResetEventHandler ProgressReset;
        public delegate void ProgressResetEventHandler();

        public event ProgressCompleteEventHandler ProgressComplete;
        public delegate void ProgressCompleteEventHandler();

        protected string mProgressStepDescription = string.Empty;

        /// <summary>
        /// Ranges from 0 to 100, but can contain decimal percentage values
        /// </summary>
        protected float mProgressPercentComplete;
        #endregion

        #region "Properties"

        public bool AbortProcessing { get; set; }

        public bool CreateModificationSummaryFile { get; set; }

        public bool CreateInspectFirstHitsFile { get; set; }

        public bool CreateInspectSynopsisFile { get; set; }

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found then it will be created using the the Fasta file specified by mFastaFilePath</remarks>
        public bool CreateProteinModsFile { get; set; }

        public clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType EnzymeMatchSpec { get; set; }

        public ePHRPErrorCodes ErrorCode => mErrorCode;

        public string ErrorMessage => GetErrorMessage();

        public string FastaFilePath { get; set; }

        public string FileVersion => GetVersionForExecutingAssembly();

        public string FileDate => mFileDate;

        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float InspectSynopsisFilePValueThreshold { get; set; }

        public string MassCorrectionTagsFilePath { get; set; }

        public string ModificationDefinitionsFilePath { get; set; }

        /// <summary>
        /// Used by clsMODaResultsProcessor and clsMODPlusResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Higher probability are higher confidence results</remarks>
        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float MSAlignSynopsisFilePValueThreshold { get; set; }

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower p-values are higher confidence results; renamed to EValue in MSGF+</remarks>
        [Obsolete("Use MSGFDBSynopsisFileEValueThreshold")]
        public float MSGFDBSynopsisFilePValueThreshold
        {
            get => MSGFDBSynopsisFileEValueThreshold;
            set => MSGFDBSynopsisFileEValueThreshold = value;
        }

        /// <summary>
        /// MSGFDB synopsis file specEValue threshold
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower SpecProb values are higher confidence results; renamed to SpecEValue in MSGF+</remarks>
        [Obsolete("Use MSGFDBSynopsisFileSpecEValueThreshold")]
        public float MSGFDBSynopsisFileSpecProbThreshold
        {
            get => MSGFDBSynopsisFileSpecEValueThreshold;
            set => MSGFDBSynopsisFileSpecEValueThreshold = value;
        }

        /// <summary>
        /// Used by clsMSGFDBResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower e-values are higher confidence results</remarks>
        public float MSGFDBSynopsisFileEValueThreshold { get; set; }

        /// <summary>
        /// clsMSGFDBResultsProcessor and clsMSPathFinderResultsProcessor
        /// </summary>
        /// <returns></returns>
        /// <remarks>Lower SpecEValue values are higher confidence results</remarks>
        public float MSGFDBSynopsisFileSpecEValueThreshold { get; set; }

        /// <summary>
        /// Typical non-zero value is 17.0027387
        /// </summary>
        /// <returns></returns>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideCTerminusMassChange { get; set; }

        /// <summary>
        ///  typical non-zero value is 1.0078246
        /// </summary>
        /// <returns></returns>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideNTerminusMassChange { get; set; }

        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        public string ProgressStepDescription => mProgressStepDescription;

        // ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public float ProgressPercentComplete => Convert.ToSingle(Math.Round(mProgressPercentComplete, 2));

        /// <summary>
        ///
        /// </summary>
        /// <returns></returns>
        /// <remarks>Used by clsInSpecTResultsProcessor and clsMSGFDBResultsProcessor (aka SearchEngineParamFileName)</remarks>
        public string SearchToolParameterFilePath { get; set; }

        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        #endregion

        public void AbortProcessingNow()
        {
            AbortProcessing = true;
        }

        public static string AutoDefinePeptideHitResultsFilePath(ePeptideHitResultsFileFormatConstants ePeptideHitResultFileFormat,
            string strSourceFolderPath, string strBaseName)
        {
            if (!string.IsNullOrEmpty(strBaseName))
            {
                switch (ePeptideHitResultFileFormat)
                {
                    case ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + SEQUEST_FIRST_HITS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.SequestSynopsisFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.XTandemXMLFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + XTANDEM_RESULTS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.InSpectTXTFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + INSPECT_RESULTS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + MSGFDB_RESULTS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.MSAlignTXTFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + MSALIGN_RESULTS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.MODaTXTFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + MODa_RESULTS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.MODPlusTXTFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + MODPlus_RESULTS_FILE_SUFFIX);
                    case ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile:

                        return Path.Combine(strSourceFolderPath, strBaseName + MSPathFinder_RESULTS_FILE_SUFFIX);
                    default:
                        // Includes ePeptideHitResultsFileFormatConstants.AutoDetermine
                        // Call AutoDefinePeptideHitResultsFilePath below sending it only strSourceFolderPath
                        break;
                }
            }

            return AutoDefinePeptideHitResultsFilePath(strSourceFolderPath);
        }

        public static string AutoDefinePeptideHitResultsFilePath(string strSourceFolderPath)
        {
            // Looks for a file ending in _syn.txt, _fht.txt, _xt.xml, or _inspect.txt in folder strSourceFolderPath
            // Returns the first matching file found

            var strMatchSpec = string.Empty;

            try
            {
                for (var intIndex = 0; intIndex <= 3; intIndex++)
                {
                    switch (intIndex)
                    {
                        case 0:
                            strMatchSpec = "*" + SEQUEST_SYNOPSIS_FILE_SUFFIX;
                            break;
                        case 1:
                            strMatchSpec = "*" + SEQUEST_FIRST_HITS_FILE_SUFFIX;
                            break;
                        case 2:
                            strMatchSpec = "*" + XTANDEM_RESULTS_FILE_SUFFIX;
                            break;
                        case 3:
                            strMatchSpec = "*" + INSPECT_RESULTS_FILE_SUFFIX;
                            break;
                    }

                    var sourceFolder = new DirectoryInfo(strSourceFolderPath);
                    foreach (var resultsFile in sourceFolder.GetFiles(strMatchSpec))
                    {
                        // If we get here, a match was found; return its path
                        return resultsFile.FullName;
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }

            // No match; return empty
            return string.Empty;
        }

        protected bool CacheProteinNamesFromFasta()
        {
            if (string.IsNullOrWhiteSpace(FastaFilePath))
            {
                // Nothing to do
                return true;
            }

            mProteinNameOrder.Clear();
            var reExtractProteinName = new Regex(@"^>([^ ]+)", RegexOptions.Compiled);

            ReportMessage("Caching protein names from the FASTA file");

            try
            {
                var proteinNumber = 0;

                using (var reader = new StreamReader(new FileStream(FastaFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(lineIn))
                            continue;

                        var reMatch = reExtractProteinName.Match(lineIn);

                        if (!reMatch.Success)
                            continue;
                        var proteinName = reMatch.Groups[1].Value;

                        if (mProteinNameOrder.ContainsKey(proteinName))
                            continue;

                        proteinNumber += 1;

                        mProteinNameOrder.Add(proteinName, proteinNumber);
                    }
                }

                ReportMessage("Cached " + mProteinNameOrder.Count + " proteins");

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error caching protein names from fasta file " + Path.GetFileName(FastaFilePath) + ": " + ex.Message);
                return false;
            }
        }

        protected bool CheckSeqToProteinMapDefined(int intUniqueSeqID, string strProteinName)
        {
            // Returns True if the sequence to protein map was already defined
            // Returns False if the mapping was not defined (will also update mSeqToProteinMap)

            bool blnExistingMapFound;

            try
            {
                if (strProteinName == null)
                    strProteinName = string.Empty;

                var strKey = intUniqueSeqID + UNIQUE_SEQ_TO_PROTEIN_MAP_SEP + strProteinName;

                if (mSeqToProteinMap.ContainsKey(strKey))
                {
                    blnExistingMapFound = true;
                }
                else
                {
                    mSeqToProteinMap.Add(strKey, 1);
                    blnExistingMapFound = false;
                }
            }
            catch (Exception)
            {
                blnExistingMapFound = false;
            }

            return blnExistingMapFound;
        }

        protected int CIntSafe(string strValue, int intDefaultValue)
        {
            try
            {
                // Note: Integer.Parse() fails if strValue contains a decimal point, even if it is "8.000"
                // Thus, we're converting to a double first, and then rounding
                return (int) Math.Round(Convert.ToDouble(strValue));
            }
            catch (Exception)
            {
                // Error converting strValue to a number; return the default
                return intDefaultValue;
            }
        }

        protected double CDblSafe(string strValue, double dblDefaultValue)
        {
            try
            {
                return double.Parse(strValue);
            }
            catch (Exception)
            {
                // Error converting strValue to a number; return the default
                return dblDefaultValue;
            }
        }

        protected float CSngSafe(string strValue, float sngDefaultValue)
        {
            try
            {
                return float.Parse(strValue);
            }
            catch (Exception)
            {
                // Error converting strValue to a number; return the default
                return sngDefaultValue;
            }
        }

        protected bool CleanupFilePaths(ref string strInputFilePath, ref string strOutputFolderPath)
        {
            // Returns True if success, False if failure

            try
            {
                // Make sure strInputFilePath points to a valid file
                var inputFile = new FileInfo(strInputFilePath);

                if (!inputFile.Exists)
                {
                    SetErrorMessage("Input file not found: " + strInputFilePath);
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                if (string.IsNullOrWhiteSpace(strOutputFolderPath))
                {
                    // Define strOutputFolderPath based on strInputFilePath
                    strOutputFolderPath = inputFile.DirectoryName;
                }

                // Make sure strOutputFolderPath points to a folder
                var outputFolder = new DirectoryInfo(strOutputFolderPath);

                if (!outputFolder.Exists)
                {
                    // strOutputFolderPath points to a non-existent folder; attempt to create it
                    try
                    {
                        outputFolder.Create();
                    }
                    catch (Exception)
                    {
                        SetErrorMessage("Invalid output folder: " + strOutputFolderPath);
                        SetErrorCode(ePHRPErrorCodes.InvalidOutputFolderPath);
                        return false;
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error cleaning up the file paths: " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.FilePathError);
                return false;
            }
        }

        /// <summary>
        /// Collapses a list of strings to a tab-delimited line of text
        /// </summary>
        /// <param name="lstFields"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected string CollapseList(List<string> lstFields)
        {
            return string.Join("\t", lstFields);
        }

        public static ePeptideHitResultsFileFormatConstants DetermineResultsFileFormat(string strFilePath)
        {
            // Examine the extension on strFilePath to determine the file format

            var strExtensionLCase = Path.GetExtension(strFilePath).ToLower();
            var baseFileName = Path.GetFileNameWithoutExtension(strFilePath).ToLower();

            if (strExtensionLCase == ".xml")
            {
                return ePeptideHitResultsFileFormatConstants.XTandemXMLFile;
            }

            if (baseFileName.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.SequestFirstHitsFile;
            }

            if (baseFileName.EndsWith(clsSequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.SequestSynopsisFile;
            }

            if (baseFileName.EndsWith(clsInSpecTResultsProcessor.FILENAME_SUFFIX_INSPECT_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.InSpectTXTFile;
            }

            if (baseFileName.EndsWith(clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile;
            }

            if (baseFileName.EndsWith(clsMSGFDBResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile;
            }
            if (baseFileName.EndsWith(clsMSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSAlignTXTFile;
            }

            if (baseFileName.EndsWith(clsMODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MODaTXTFile;
            }

            if (baseFileName.EndsWith(clsMODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MODPlusTXTFile;
            }

            if (baseFileName.EndsWith(clsMSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                return ePeptideHitResultsFileFormatConstants.MSPathFinderTSVFile;
            }

            if (strExtensionLCase == ".tsv")
            {
                // Assume this is an MSGF+ TSV file
                return ePeptideHitResultsFileFormatConstants.MSGFDbTXTFile;
            }

            // Unknown extension
            return ePeptideHitResultsFileFormatConstants.AutoDetermine;
        }

        protected void CloseSequenceOutputFiles()
        {
            try
            {
                if (mResultToSeqMapFile != null)
                {
                    mResultToSeqMapFile.Close();
                    mResultToSeqMapFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }

            try
            {
                if (mSeqInfoFile != null)
                {
                    mSeqInfoFile.Close();
                    mSeqInfoFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }

            try
            {
                if (mModDetailsFile != null)
                {
                    mModDetailsFile.Close();
                    mModDetailsFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }

            try
            {
                if (mSeqToProteinMapFile != null)
                {
                    mSeqToProteinMapFile.Close();
                    mSeqToProteinMapFile = null;
                }
            }
            catch (Exception)
            {
                // Ignore errors
            }
        }

        protected void ComputePseudoPeptideLocInProtein(clsSearchResultsBaseClass objSearchResult)
        {
            // Set these to 1 and 10000 since MSGFDB, Sequest, and Inspect results files do not contain protein sequence information
            // If we find later that the peptide sequence spans the length of the protein, we'll revise .ProteinSeqResidueNumberEnd as needed
            objSearchResult.ProteinSeqResidueNumberStart = 1;
            objSearchResult.ProteinSeqResidueNumberEnd = 10000;

            if (objSearchResult.PeptidePreResidues.Trim().EndsWith(clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST.ToString()))
            {
                // The peptide is at the N-Terminus of the protein
                objSearchResult.PeptideLocInProteinStart = objSearchResult.ProteinSeqResidueNumberStart;
                objSearchResult.PeptideLocInProteinEnd = objSearchResult.PeptideLocInProteinStart + objSearchResult.PeptideCleanSequence.Length - 1;

                if (objSearchResult.PeptidePostResidues.Trim()[0] == clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST)
                {
                    // The peptide spans the entire length of the protein
                    objSearchResult.ProteinSeqResidueNumberEnd = objSearchResult.PeptideLocInProteinEnd;
                }
                else
                {
                    if (objSearchResult.PeptideLocInProteinEnd > objSearchResult.ProteinSeqResidueNumberEnd)
                    {
                        // The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                        objSearchResult.ProteinSeqResidueNumberEnd = objSearchResult.PeptideLocInProteinEnd + 1;
                    }
                }
            }
            else if (objSearchResult.PeptidePostResidues.Trim().StartsWith(clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST.ToString()))
            {
                // The peptide is at the C-Terminus of the protein
                objSearchResult.PeptideLocInProteinEnd = objSearchResult.ProteinSeqResidueNumberEnd;
                objSearchResult.PeptideLocInProteinStart = objSearchResult.PeptideLocInProteinEnd - objSearchResult.PeptideCleanSequence.Length + 1;

                if (objSearchResult.PeptideLocInProteinStart < objSearchResult.ProteinSeqResidueNumberStart)
                {
                    // The peptide is more than 10000 characters long; this is highly unlikely
                    objSearchResult.ProteinSeqResidueNumberEnd = objSearchResult.ProteinSeqResidueNumberStart + 1 + objSearchResult.PeptideCleanSequence.Length;
                    objSearchResult.PeptideLocInProteinEnd = objSearchResult.ProteinSeqResidueNumberEnd;
                    objSearchResult.PeptideLocInProteinStart = objSearchResult.PeptideLocInProteinEnd - objSearchResult.PeptideCleanSequence.Length + 1;
                }
            }
            else
            {
                objSearchResult.PeptideLocInProteinStart = objSearchResult.ProteinSeqResidueNumberStart + 1;
                objSearchResult.PeptideLocInProteinEnd = objSearchResult.PeptideLocInProteinStart + objSearchResult.PeptideCleanSequence.Length - 1;

                if (objSearchResult.PeptideLocInProteinEnd > objSearchResult.ProteinSeqResidueNumberEnd)
                {
                    // The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                    objSearchResult.ProteinSeqResidueNumberEnd = objSearchResult.PeptideLocInProteinEnd + 1;
                }
            }
        }

        protected virtual string ConstructPepToProteinMapFilePath(string strInputFilePath, string strOutputFolderPath, bool MTS)
        {
            var strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);

            if (MTS)
            {
                strPepToProteinMapFilePath += FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + "MTS.txt";
            }
            else
            {
                strPepToProteinMapFilePath += FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + ".txt";
            }

            if (string.IsNullOrWhiteSpace(strOutputFolderPath))
            {
                var inputFile = new FileInfo(strInputFilePath);
                strOutputFolderPath = inputFile.DirectoryName;
            }

            strPepToProteinMapFilePath = Path.Combine(strOutputFolderPath, strPepToProteinMapFilePath);

            return strPepToProteinMapFilePath;
        }

        /// <summary>
        /// Use the PeptideToProteinMapEngine to create the Peptide to Protein map file for the file or files in lstSourcePHRPDataFiles
        /// </summary>
        /// <param name="lstSourcePHRPDataFiles"></param>
        /// <param name="strMTSPepToProteinMapFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreatePepToProteinMapFile(List<string> lstSourcePHRPDataFiles, string strMTSPepToProteinMapFilePath)
        {
            var blnSuccess = false;

            try
            {
                if (string.IsNullOrEmpty(strMTSPepToProteinMapFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because strMTSPepToProteinMapFilePath is empty; likely a programming bug");
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                if (string.IsNullOrEmpty(FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File Path is not defined");
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                if (!File.Exists(FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File was not found: " + FastaFilePath);
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                // Verify that the fasta file is not a DNA-sequence based fasta file
                blnSuccess = ValidateProteinFastaFile(FastaFilePath);
                if (!blnSuccess)
                {
                    return false;
                }

                Console.WriteLine();
                Console.WriteLine();
                UpdateProgress("Creating Peptide to Protein Map file", PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE);

                // Initialize items
                var mtsPepToProteinMapFile = new FileInfo(strMTSPepToProteinMapFilePath);
                var strOutputFolderPath = mtsPepToProteinMapFile.DirectoryName;

                var htPeptideToProteinMapResults = new Hashtable();

                var objPeptideToProteinMapper = new PeptideToProteinMapEngine.clsPeptideToProteinMapEngine
                {
                    DeleteTempFiles = true,
                    IgnoreILDifferences = false,
                    InspectParameterFilePath = string.Empty,
                    LogMessagesToFile = false,
                    MatchPeptidePrefixAndSuffixToProtein = false,
                    OutputProteinSequence = false,
                    PeptideInputFileFormat = PeptideToProteinMapEngine.clsPeptideToProteinMapEngine.ePeptideInputFileFormatConstants.PHRPFile,
                    PeptideFileSkipFirstLine = false,
                    ProteinDataRemoveSymbolCharacters = true,
                    ProteinInputFilePath = FastaFilePath,
                    SaveProteinToPeptideMappingFile = true,
                    SearchAllProteinsForPeptideSequence = true,
                    SearchAllProteinsSkipCoverageComputationSteps = true,
                    ShowMessages = true
                };

                using (var swMTSpepToProteinMapFile =
                       new StreamWriter(new FileStream(strMTSPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    foreach (var strInputFilePath in lstSourcePHRPDataFiles)
                    {
                        var strResultsFilePath = Path.GetFileNameWithoutExtension(strInputFilePath) +
                            PeptideToProteinMapEngine.clsPeptideToProteinMapEngine.FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING;

                        strResultsFilePath = Path.Combine(strOutputFolderPath, strResultsFilePath);

                        // Make sure the results file doesn't already exist
                        DeleteFileIgnoreErrors(strResultsFilePath);

                        objPeptideToProteinMapper.ProgressChanged += PeptideToProteinMapper_ProgressChanged;
                        mNextPeptideToProteinMapperLevel = 25;

                        blnSuccess = objPeptideToProteinMapper.ProcessFile(strInputFilePath, strOutputFolderPath, string.Empty, true);

                        objPeptideToProteinMapper.ProgressChanged -= PeptideToProteinMapper_ProgressChanged;

                        if (blnSuccess)
                        {
                            if (!File.Exists(strResultsFilePath))
                            {
                                SetErrorMessage("Peptide to protein mapping file was not created for " + strInputFilePath);
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                blnSuccess = false;
                                break;
                            }
                            blnSuccess = ValidatePeptideToProteinMapResults(strResultsFilePath, IgnorePeptideToProteinMapperErrors);
                        }
                        else
                        {
                            if (string.IsNullOrWhiteSpace(objPeptideToProteinMapper.GetErrorMessage()) && objPeptideToProteinMapper.StatusMessage.ToLower().Contains("error"))
                            {
                                SetErrorMessage("Error running clsPeptideToProteinMapEngine: " + objPeptideToProteinMapper.StatusMessage);
                            }
                            else
                            {
                                if (objPeptideToProteinMapper.StatusMessage.Length > 0)
                                {
                                    SetErrorMessage("clsPeptideToProteinMapEngine status: " + objPeptideToProteinMapper.StatusMessage);
                                }
                                SetErrorMessage("Error running clsPeptideToProteinMapEngine: " + objPeptideToProteinMapper.GetErrorMessage());
                            }

                            if (IgnorePeptideToProteinMapperErrors)
                            {
                                ReportWarning("Ignoring protein mapping error since 'IgnorePeptideToProteinMapperErrors' = True");

                                if (File.Exists(strResultsFilePath))
                                {
                                    blnSuccess = ValidatePeptideToProteinMapResults(strResultsFilePath, IgnorePeptideToProteinMapperErrors);
                                }
                                else
                                {
                                    mErrorMessage = string.Empty;
                                    mErrorCode = ePHRPErrorCodes.NoError;
                                    blnSuccess = true;
                                }
                            }
                            else
                            {
                                SetErrorMessage("Error in CreatePepToProteinMapFile");
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                            }
                        }

                        if (!File.Exists(strResultsFilePath))
                        {
                            continue;
                        }

                        // Read the newly created file and append new entries to strMTSPepToProteinMapFilePath
                        using (var srResultsFile = new StreamReader(new FileStream(strResultsFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                        {
                            while (!srResultsFile.EndOfStream)
                            {
                                var strLineIn = srResultsFile.ReadLine();

                                if (!string.IsNullOrWhiteSpace(strLineIn))
                                {
                                    var strSplitLine = strLineIn.Split(new [] {'\t'}, 2);
                                    if (strSplitLine.Length >= 2)
                                    {
                                        var strPeptideAndProteinKey = strSplitLine[0] + "_" + strSplitLine[1];

                                        if (!htPeptideToProteinMapResults.ContainsKey(strPeptideAndProteinKey))
                                        {
                                            htPeptideToProteinMapResults.Add(strPeptideAndProteinKey, 0);
                                            swMTSpepToProteinMapFile.WriteLine(strLineIn);
                                        }
                                    }
                                }
                            }
                        }

                        // Delete the interim results file
                        DeleteFileIgnoreErrors(strResultsFilePath);
                    }
                }

                objPeptideToProteinMapper.CloseLogFileNow();
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreatePepToProteinMapFile:" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
            }

            return blnSuccess;
        }

        /// <summary>
        /// Create the protein mod details file for the specified PHRP data file
        /// </summary>
        /// <param name="strPHRPDataFilePath"></param>
        /// <param name="strOutputFolderPath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool CreateProteinModDetailsFile(string strPHRPDataFilePath, string strOutputFolderPath)
        {
            var blnSuccess = false;

            try
            {
                var inputFile = new FileInfo(strPHRPDataFilePath);

                var lstSourcePHRPDataFiles = new List<string> {
                    inputFile.FullName
                };

                var strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, strOutputFolderPath, MTS: true);

                blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);

                if (blnSuccess)
                {
                    blnSuccess = CreateProteinModDetailsFile(strPHRPDataFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Unknown);
                }
            }
            catch (Exception ex)
            {
                ReportWarning("Error in CreateProteinModDetailsFile: " + ex.Message);
            }

            return blnSuccess;
        }

        public bool CreateProteinModDetailsFile(
            string strPHRPDataFilePath,
            string strOutputFolderPath,
            string strMTSPepToProteinMapFilePath,
            clsPHRPReader.ePeptideHitResultType ePHRPResultType)
        {
            var blnSuccess = false;

            try
            {
                Console.WriteLine();

                var sngProgressAtStart = mProgressPercentComplete;
                UpdateProgress("Creating the Protein Mod Details file", PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE);

                // Confirm that the PHRP data file exists
                var phrpDataFile = new FileInfo(strPHRPDataFilePath);
                if (!phrpDataFile.Exists)
                {
                    SetErrorMessage("PHRP data file not found in CreateProteinModDetailsFile: " + strPHRPDataFilePath);
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                // Confirm that the _PepToProtMapMTS.txt file exists
                if (string.IsNullOrEmpty(strMTSPepToProteinMapFilePath))
                {
                    SetErrorMessage("Cannot create the ProteinMods file because strMTSPepToProteinMapFilePath is empty");
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    return false;
                }

                // Initialize lstPepToProteinMapping
                var lstPepToProteinMapping = new List<udtPepToProteinMappingType>();

                // Read the _PepToProtMapMTS file
                blnSuccess = LoadPeptideToProteinMapInfo(strMTSPepToProteinMapFilePath, lstPepToProteinMapping, out _);
                if (!blnSuccess)
                {
                    return false;
                }

                // Assure that lstPepToProteinMapping is sorted on peptide
                if (lstPepToProteinMapping.Count > 1)
                {
                    lstPepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                var strProteinModsFilePath = ReplaceFilenameSuffix(phrpDataFile, FILENAME_SUFFIX_PROTEIN_MODS);
                if (!string.IsNullOrEmpty(strOutputFolderPath))
                {
                    strProteinModsFilePath = Path.Combine(strOutputFolderPath, Path.GetFileName(strProteinModsFilePath));
                }

                var intPSMCount = 0;
                var intPSMCountSkippedSinceReversedOrScrambledProtein = 0;

                // Create a ProteinMods file parallel to the PHRP file
                using (var swProteinModsFile = new StreamWriter(new FileStream(strProteinModsFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    // Write the header line
                    swProteinModsFile.WriteLine(COLUMN_NAME_RESULTID + "\t" +
                                                COLUMN_NAME_PEPTIDE + "\t" +
                                                COLUMN_NAME_UNIQUE_SEQ_ID + "\t" +
                                                COLUMN_NAME_PROTEIN_NAME + "\t" +
                                                COLUMN_NAME_RESIDUE + "\t" +
                                                COLUMN_NAME_PROTEIN_RESIDUE_NUMBER + "\t" +
                                                COLUMN_NAME_RESIDUE_MOD_NAME + "\t" +
                                                COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER + "\t" +
                                                COLUMN_NAME_MSGF_SPECPROB);

                    var blnLoadMSGFResults = ePHRPResultType != clsPHRPReader.ePeptideHitResultType.MSGFDB;

                    // Update the Mass Calculator to use the one tracked by this class
                    // (since this class's calculator knows about custom amino acids and custom charge carriers)
                    var oStartupOptions = new clsPHRPStartupOptions
                    {
                        LoadModsAndSeqInfo = true,
                        LoadMSGFResults = blnLoadMSGFResults,
                        LoadScanStatsData = false,
                        PeptideMassCalculator = mPeptideSeqMassCalculator
                    };

                    using (var objReader = new clsPHRPReader(strPHRPDataFilePath, ePHRPResultType, oStartupOptions))
                    {
                        objReader.EchoMessagesToConsole = false;
                        objReader.SkipDuplicatePSMs = true;

                        foreach (var strErrorMessage in objReader.ErrorMessages)
                        {
                            OnErrorEvent(strErrorMessage);
                        }
                        RegisterEvents(objReader);

                        foreach (var strWarningMessage in objReader.WarningMessages)
                        {
                            var msg = strWarningMessage;
                            if (strWarningMessage.StartsWith("MSGF file not found", StringComparison.InvariantCultureIgnoreCase))
                            {
                                msg = "MSGF file not found; column " + COLUMN_NAME_MSGF_SPECPROB + " will not have any data";
                            }
                            ReportWarning(msg);
                        }

                        objReader.ClearErrors();
                        objReader.ClearWarnings();

                        while (objReader.MoveNext())
                        {
                            // Use binary search to find this peptide in lstPepToProteinMapping
                            var intPepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(lstPepToProteinMapping, objReader.CurrentPSM.Peptide);

                            if (intPepToProteinMapIndex >= 0)
                            {
                                do
                                {
                                    intPSMCount += 1;

                                    var blnSkipProtein = false;
                                    if (!ProteinModsFileIncludesReversedProteins)
                                    {
                                        blnSkipProtein = IsReversedProtein(lstPepToProteinMapping[intPepToProteinMapIndex].Protein);
                                        if (blnSkipProtein)
                                        {
                                            intPSMCountSkippedSinceReversedOrScrambledProtein += 1;
                                        }
                                    }

                                    if (!blnSkipProtein)
                                    {
                                        foreach (var objMod in objReader.CurrentPSM.ModifiedResidues)
                                        {
                                            var intResidueLocInProtein = lstPepToProteinMapping[intPepToProteinMapIndex].ResidueStart + objMod.ResidueLocInPeptide - 1;
                                            string strResidue;

                                            if (IsLetterAtoZ(objMod.Residue))
                                            {
                                                strResidue = objMod.Residue.ToString();
                                            }
                                            else
                                            {
                                                var strCleanSequence = objReader.CurrentPSM.PeptideCleanSequence;

                                                if (objMod.ResidueLocInPeptide < 1)
                                                {
                                                    // This shouldn't be the case, but we'll check for it anyway
                                                    strResidue = strCleanSequence.Substring(0, 1);
                                                }
                                                else if (objMod.ResidueLocInPeptide > strCleanSequence.Length)
                                                {
                                                    // This shouldn't be the case, but we'll check for it anyway
                                                    strResidue = strCleanSequence.Substring(strCleanSequence.Length - 1, 1);
                                                }
                                                else
                                                {
                                                    strResidue = strCleanSequence.Substring(objMod.ResidueLocInPeptide - 1, 1);
                                                }
                                            }

                                            if (lstPepToProteinMapping[intPepToProteinMapIndex].Protein == PROTEIN_NAME_NO_MATCH && IsReversedProtein(objReader.CurrentPSM.ProteinFirst))
                                            {
                                                // Skip this result
                                                intPSMCountSkippedSinceReversedOrScrambledProtein += 1;
                                            }
                                            else
                                            {
                                                swProteinModsFile.WriteLine(objReader.CurrentPSM.ResultID + "\t" +
                                                                            objReader.CurrentPSM.Peptide + "\t" +
                                                                            objReader.CurrentPSM.SeqID + "\t" +
                                                                            lstPepToProteinMapping[intPepToProteinMapIndex].Protein + "\t" +
                                                                            strResidue + "\t" +
                                                                            intResidueLocInProtein + "\t" +
                                                                            objMod.ModDefinition.MassCorrectionTag + "\t" +
                                                                            objMod.ResidueLocInPeptide + "\t" +
                                                                            objReader.CurrentPSM.MSGFSpecEValue);
                                            }
                                        }
                                    }

                                    intPepToProteinMapIndex += 1;
                                } while (intPepToProteinMapIndex < lstPepToProteinMapping.Count &&
                                         objReader.CurrentPSM.Peptide == lstPepToProteinMapping[intPepToProteinMapIndex].Peptide);
                            }
                            else
                            {
                                ReportWarning("Peptide not found in lstPepToProteinMapping: " + objReader.CurrentPSM.Peptide);
                            }

                            UpdateProgress(sngProgressAtStart + objReader.PercentComplete * (100 - sngProgressAtStart) / 100);
                        }

                        if (intPSMCount > 0)
                        {
                            if (intPSMCountSkippedSinceReversedOrScrambledProtein == intPSMCount)
                            {
                                Console.WriteLine();
                                ReportWarning("All PSMs map to reversed or scrambled proteins; the _ProteinMods.txt file is empty");
                            }
                            else if (intPSMCountSkippedSinceReversedOrScrambledProtein > 0)
                            {
                                Console.WriteLine();
                                Console.WriteLine("Note: skipped " + intPSMCountSkippedSinceReversedOrScrambledProtein + " / " + intPSMCount + " PSMs that map to reversed or scrambled proteins while creating the _ProteinMods.txt file");
                            }
                        }
                    }
                }

                blnSuccess = true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateProteinModDetailsFile:" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
            }

            return blnSuccess;
        }

        protected void DeleteFileIgnoreErrors(string strFilePath)
        {
            try
            {
                if (File.Exists(strFilePath))
                {
                    Thread.Sleep(200);
                    File.Delete(strFilePath);
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        protected void ExpandListIfRequired<T>(List<T> lstItems, int countToAdd, int largeListThreshold = 1000000)
        {
            if (lstItems.Count > largeListThreshold && lstItems.Count + countToAdd > lstItems.Capacity)
            {
                // .NET by default will double the size of the list to accomodate these new items
                // Instead, expand the list by 20% of the current size
                lstItems.Capacity = lstItems.Capacity + Convert.ToInt32(lstItems.Count / 5);
            }
        }

        private readonly IComparer<udtPepToProteinMappingType> objPeptideSearchComparer = new PepToProteinMappingPeptideSearchComparer();

        protected int FindFirstMatchInPepToProteinMapping(List<udtPepToProteinMappingType> lstPepToProteinMapping, string strPeptideToFind)
        {
            // Use binary search to find this peptide in lstPepToProteinMapping
            var udtItemToFind = new udtPepToProteinMappingType {
                Peptide = strPeptideToFind
            };

            var intPepToProteinMapIndex = lstPepToProteinMapping.BinarySearch(udtItemToFind, objPeptideSearchComparer);

            if (intPepToProteinMapIndex > 0)
            {
                // Step Backward until the first match is found
                while (intPepToProteinMapIndex > 0 && lstPepToProteinMapping[intPepToProteinMapIndex - 1].Peptide == strPeptideToFind)
                {
                    intPepToProteinMapIndex -= 1;
                }
            }

            return intPepToProteinMapIndex;
        }

        protected string GetCleanSequence(string strSequenceWithMods)
        {
            return GetCleanSequence(strSequenceWithMods, out _, out _);
        }

        protected string GetCleanSequence(string strSequenceWithMods, out string strPrefix, out string strSuffix)
        {
            strPrefix = string.Empty;
            strSuffix = string.Empty;

            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, out var strPrimarySequence, out strPrefix, out strSuffix))
            {
                // Remove any non-letter characters when calling .ComputeCleavageState()

                strPrimarySequence = mReplaceSymbols.Replace(strPrimarySequence, string.Empty);
            }
            else
            {
                // Unable to determine cleavage-state
                strPrimarySequence = mReplaceSymbols.Replace(strSequenceWithMods, string.Empty);
            }

            return strPrimarySequence;
        }

        /// <summary>
        /// If intColumnIndex is >= 0 then updates strValue with the value at strSplitLine[intColumnIndex]
        /// Otherwise, updates strValue to String.Empty
        /// </summary>
        /// <returns>True if intColumnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] strSplitLine, int intColumnIndex, out string strValue)
        {
            return GetColumnValue(strSplitLine, intColumnIndex, out strValue, string.Empty);
        }

        /// <summary>
        /// If intColumnIndex is >= 0 then updates intValue with the value at strSplitLine[intColumnIndex]
        /// Otherwise, updates intValue to 0
        /// </summary>
        /// <returns>True if intColumnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] strSplitLine, int intColumnIndex, out int intValue)
        {
            return GetColumnValue(strSplitLine, intColumnIndex, out intValue, 0);
        }

        /// <summary>
        /// If intColumnIndex is >= 0 then updates strValue with the value at strSplitLine[intColumnIndex]
        /// Otherwise, updates strValue to strValueIfMissing
        /// </summary>
        /// <returns>True if intColumnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] strSplitLine, int intColumnIndex, out string strValue, string strValueIfMissing)
        {
            strValueIfMissing = "";
            if (intColumnIndex >= 0 && intColumnIndex < strSplitLine.Length)
            {
                strValue = string.Copy(strSplitLine[intColumnIndex]);
                return true;
            }

            strValue = string.Copy(strValueIfMissing);
            return false;
        }

        /// <summary>
        /// If intColumnIndex is >= 0 then updates intValue with the value at strSplitLine[intColumnIndex]
        /// Otherwise, updates strValue to intValueIfMissing
        /// </summary>
        /// <returns>True if intColumnIndex >= 0</returns>
        /// <remarks></remarks>
        protected bool GetColumnValue(string[] strSplitLine, int intColumnIndex, out int intValue, int intValueIfMissing)
        {
            if (GetColumnValue(strSplitLine, intColumnIndex, out var strValue, intValueIfMissing.ToString()))
            {
                if (int.TryParse(strValue, out intValue))
                {
                    return true;
                }

                intValue = intValueIfMissing;
                return false;
            }

            intValue = intValueIfMissing;
            return false;
        }

        protected string GetErrorMessage()
        {
            // Returns String.Empty if no error

            string strMessage;

            switch (ErrorCode)
            {
                case ePHRPErrorCodes.NoError:
                    strMessage = string.Empty;
                    break;
                case ePHRPErrorCodes.InvalidInputFilePath:
                    strMessage = "Invalid input file path";
                    break;
                case ePHRPErrorCodes.InvalidOutputFolderPath:
                    strMessage = "Invalid output folder path";
                    break;
                case ePHRPErrorCodes.ParameterFileNotFound:
                    strMessage = "Parameter file not found";
                    break;
                case ePHRPErrorCodes.MassCorrectionTagsFileNotFound:
                    strMessage = "Mass correction tags file not found";
                    break;
                case ePHRPErrorCodes.ModificationDefinitionFileNotFound:
                    strMessage = "Modification definition file not found";

                    break;
                case ePHRPErrorCodes.ErrorReadingInputFile:
                    strMessage = "Error reading input file";
                    break;
                case ePHRPErrorCodes.ErrorCreatingOutputFiles:
                    strMessage = "Error creating output files";
                    break;
                case ePHRPErrorCodes.ErrorReadingParameterFile:
                    strMessage = "Invalid parameter file";
                    break;
                case ePHRPErrorCodes.ErrorReadingMassCorrectionTagsFile:
                    strMessage = "Error reading mass correction tags file";
                    break;
                case ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile:
                    strMessage = "Error reading modification definitions file";

                    break;
                case ePHRPErrorCodes.FilePathError:
                    strMessage = "General file path error";
                    break;
                case ePHRPErrorCodes.UnspecifiedError:
                    strMessage = "Unspecified error";
                    break;
                default:
                    // This shouldn't happen
                    strMessage = "Unknown error state";
                    break;
            }

            if (mErrorMessage.Length > 0)
            {
                if (strMessage.Length > 0)
                {
                    strMessage += "; ";
                }
                strMessage += mErrorMessage;
            }

            return strMessage;
        }

        private string GetVersionForExecutingAssembly()
        {
            string strVersion;

            try
            {
                strVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            }
            catch (Exception)
            {
                strVersion = "??.??.??.??";
            }

            return strVersion;
        }

        private void InitializeLocalVariables()
        {
            mErrorCode = ePHRPErrorCodes.NoError;
            mErrorMessage = string.Empty;
            WarnMissingParameterFileSection = false;

            CreateModificationSummaryFile = true;
            CreateProteinModsFile = false;
            FastaFilePath = string.Empty;
            IgnorePeptideToProteinMapperErrors = false;
            ProteinModsFileIncludesReversedProteins = false;
            UseExistingMTSPepToProteinMapFile = false;

            CreateInspectFirstHitsFile = true;
            CreateInspectSynopsisFile = true;

            MassCorrectionTagsFilePath = string.Empty;
            ModificationDefinitionsFilePath = string.Empty;
            SearchToolParameterFilePath = string.Empty;

            InspectSynopsisFilePValueThreshold = clsInSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MODaMODPlusSynopsisFileProbabilityThreshold = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            MSAlignSynopsisFilePValueThreshold = clsMSAlignResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MSGFDBSynopsisFileEValueThreshold = clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MSGFDBSynopsisFileSpecEValueThreshold = clsMSGFDBResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            EnzymeMatchSpec = clsPeptideCleavageStateCalculator.GetDefaultEnzymeMatchSpec();

            PeptideNTerminusMassChange = clsPeptideMassCalculator.DEFAULT_N_TERMINUS_MASS_CHANGE;
            PeptideCTerminusMassChange = clsPeptideMassCalculator.DEFAULT_C_TERMINUS_MASS_CHANGE;
        }

        protected bool InitializeSequenceOutputFiles(string strBaseOutputFilePath)
        {
            // Initializes the StreamWriter objects using strBaseOutputFilePath as a base name and replacing the suffix with the default suffix names
            // Returns True if success; does not catch errors; they will be thrown to the calling function if they occur

            var outputFileInfo = new FileInfo(strBaseOutputFilePath);

            // Initialize the file paths based on strBaseOutputFilePath
            var strResultToSeqMapFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_RESULT_TO_SEQ_MAP);
            var strSeqInfoFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_SEQ_INFO);
            var strModDetailsFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_MOD_DETAILS);
            var strSeqToProteinMapFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP);

            // Clear the unique sequences container
            mUniqueSequences.Clear();

            // Clear the sequence to protein map
            mSeqToProteinMap.Clear();

            // Initialize the ResultToSeqMap file
            mResultToSeqMapFile = new StreamWriter(strResultToSeqMapFilePath);
            mResultToSeqMapFile.WriteLine("Result_ID" + SEP_CHAR + COLUMN_NAME_UNIQUE_SEQ_ID);

            // Initialize the SeqInfo file
            mSeqInfoFile = new StreamWriter(strSeqInfoFilePath, false);
            mSeqInfoFile.WriteLine(COLUMN_NAME_UNIQUE_SEQ_ID + SEP_CHAR +
                                   "Mod_Count" + SEP_CHAR +
                                   "Mod_Description" + SEP_CHAR +
                                   "Monoisotopic_Mass");

            // Initialize the ModDetails file
            mModDetailsFile = new StreamWriter(strModDetailsFilePath);
            mModDetailsFile.WriteLine(COLUMN_NAME_UNIQUE_SEQ_ID + SEP_CHAR +
                                      "Mass_Correction_Tag" + SEP_CHAR +
                                      "Position");

            // Initialize the SeqToProtein map file
            mSeqToProteinMapFile = new StreamWriter(strSeqToProteinMapFilePath, false);
            mSeqToProteinMapFile.WriteLine(COLUMN_NAME_UNIQUE_SEQ_ID + SEP_CHAR +
                                           "Cleavage_State" + SEP_CHAR +
                                           "Terminus_State" + SEP_CHAR +
                                           COLUMN_NAME_PROTEIN_NAME + SEP_CHAR +
                                           "Protein_Expectation_Value_Log(e)" + SEP_CHAR +
                                           "Protein_Intensity_Log(I)");

            return true;
        }

        private static readonly Regex RegexIsLetter = new Regex(@"[A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Returns true if the character is a letter between A and Z or a and z
        /// </summary>
        /// <param name="chChar">Character to examine</param>
        /// <returns></returns>
        /// <remarks>The Char.IsLetter() function returns True for "�" and various other Unicode ModifierLetter characters; use this function to only return True for normal letters between A and Z</remarks>
        public static bool IsLetterAtoZ(char chChar)
        {
            if (RegexIsLetter.IsMatch(chChar.ToString()))
            {
                return true;
            }

            return false;
        }

        protected bool IsReversedProtein(string strProteinName)
        {
            if (strProteinName.StartsWith("reversed_", StringComparison.InvariantCultureIgnoreCase))
            {
                // Used in DMS-generated protein collections
                return true;
            }

            if (strProteinName.StartsWith("REV_", StringComparison.InvariantCultureIgnoreCase))
            {
                // Used by MSGFDB
                return true;
            }

            if (strProteinName.StartsWith("scrambled_", StringComparison.InvariantCultureIgnoreCase))
            {
                // Used in DMS-generated protein collections
                return true;
            }
            if (strProteinName.StartsWith("xxx_", StringComparison.InvariantCultureIgnoreCase))
            {
                // Used by MSGF+
                return true;
            }

            if (strProteinName.StartsWith("xxx.", StringComparison.InvariantCultureIgnoreCase))
            {
                // Used by Inspect
                return true;
            }

            if (strProteinName.EndsWith(":reversed", StringComparison.InvariantCultureIgnoreCase))
            {
                // Used by X!Tandem
                return true;
            }

            return false;
        }

        protected virtual bool LoadParameterFileSettings(string strParameterFilePath)
        {
            const string OPTIONS_SECTION = "PeptideHitResultsProcessorOptions";

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
                        SetErrorCode(ePHRPErrorCodes.ParameterFileNotFound);
                        return false;
                    }
                }

                if (objSettingsFile.LoadSettings(strParameterFilePath))
                {
                    if (!objSettingsFile.SectionPresent(OPTIONS_SECTION))
                    {
                        // Section OPTIONS_SECTION was not found in the parameter file; warn the user if mWarnMissingParameterFileSection = True
                        if (WarnMissingParameterFileSection)
                        {
                            SetErrorMessage("The node '<section name=\"" + OPTIONS_SECTION + "\"> was not found in the parameter file: " + strParameterFilePath);
                            SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile);
                            return false;
                        }
                        return true;
                    }

                    MassCorrectionTagsFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "MassCorrectionTagsFilePath", MassCorrectionTagsFilePath);
                    ModificationDefinitionsFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "ModificationDefinitionsFilePath", ModificationDefinitionsFilePath);
                    SearchToolParameterFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "SearchToolParameterFilePath", SearchToolParameterFilePath);

                    CreateModificationSummaryFile = objSettingsFile.GetParam(OPTIONS_SECTION, "CreateModificationSummaryFile", CreateModificationSummaryFile);

                    CreateProteinModsFile = objSettingsFile.GetParam(OPTIONS_SECTION, "CreateProteinModsFile", CreateProteinModsFile);
                    FastaFilePath = objSettingsFile.GetParam(OPTIONS_SECTION, "FastaFilePath", FastaFilePath);
                    ProteinModsFileIncludesReversedProteins = objSettingsFile.GetParam(OPTIONS_SECTION, "ProteinModsFileIncludesReversedProteins", ProteinModsFileIncludesReversedProteins);
                    UseExistingMTSPepToProteinMapFile = objSettingsFile.GetParam(OPTIONS_SECTION, "UseExistingMTSPepToProteinMapFile", UseExistingMTSPepToProteinMapFile);

                    var strLeftResidueRegEx = string.Copy(EnzymeMatchSpec.LeftResidueRegEx);
                    var strRightResidueRegEx = string.Copy(EnzymeMatchSpec.RightResidueRegEx);

                    strLeftResidueRegEx = objSettingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecLeftResidue", strLeftResidueRegEx, out var blnValueNotPresent);
                    if (!blnValueNotPresent)
                    {
                        strRightResidueRegEx = objSettingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecRightResidue", strRightResidueRegEx, out blnValueNotPresent);

                        if (!blnValueNotPresent)
                        {
                            EnzymeMatchSpec = new clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType(strLeftResidueRegEx, strRightResidueRegEx);
                        }
                    }

                    PeptideNTerminusMassChange = objSettingsFile.GetParam(OPTIONS_SECTION, "PeptideNTerminusMassChange", PeptideNTerminusMassChange);
                    PeptideCTerminusMassChange = objSettingsFile.GetParam(OPTIONS_SECTION, "PeptideCTerminusMassChange", PeptideCTerminusMassChange);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadParameterFileSettings:" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Load the PeptideToProteinMap information
        /// </summary>
        /// <param name="strPepToProteinMapFilePath">File to read</param>
        /// <param name="lstPepToProteinMapping">Output parameter: peptide to protein mapping (calling function must pre-initialize the list)</param>
        /// <param name="strHeaderLine">Output parameter: Header line text</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool LoadPeptideToProteinMapInfo(
            string strPepToProteinMapFilePath,
            List<udtPepToProteinMappingType> lstPepToProteinMapping,
            out string strHeaderLine)
        {
            strHeaderLine = "";

            bool blnSuccess;

            try
            {
                // Initialize the output parameters
                lstPepToProteinMapping.Clear();
                strHeaderLine = string.Empty;

                if (string.IsNullOrWhiteSpace(strPepToProteinMapFilePath))
                {
                    SetErrorMessage("Warning: PepToProteinMap file is not defined");
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                if (!File.Exists(strPepToProteinMapFilePath))
                {
                    SetErrorMessage("Warning: PepToProteinMap file does not exist: " + strPepToProteinMapFilePath);
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                // Open strProteinToPeptideMappingFilePath for reading
                using (var srInFile = new StreamReader(new FileStream(strPepToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    var intLinesRead = 0;
                    while (!srInFile.EndOfStream)
                    {
                        var strLineIn = srInFile.ReadLine();
                        if (string.IsNullOrEmpty(strLineIn))
                            continue;

                        // Split the line on tabs
                        var strSplitLine = strLineIn.Trim().Split('\t');

                        if (strSplitLine.Length >= 4)
                        {
                            if (intLinesRead == 0 && !int.TryParse(strSplitLine[2], out _))
                            {
                                // Header line; cache it
                                strHeaderLine = string.Copy(strLineIn);
                            }
                            else
                            {
                                var udtPepToProteinMappingEntry = default(udtPepToProteinMappingType);
                                udtPepToProteinMappingEntry.Peptide = string.Copy(strSplitLine[0]);
                                udtPepToProteinMappingEntry.Protein = string.Copy(strSplitLine[1]);
                                int.TryParse(strSplitLine[2], out udtPepToProteinMappingEntry.ResidueStart);
                                int.TryParse(strSplitLine[3], out udtPepToProteinMappingEntry.ResidueEnd);

                                ExpandListIfRequired(lstPepToProteinMapping, 1);

                                lstPepToProteinMapping.Add(udtPepToProteinMappingEntry);
                            }
                        }
                    }
                }

                blnSuccess = true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the Peptide to Protein Map File (" + Path.GetFileName(strPepToProteinMapFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        protected void OperationComplete()
        {
            ProgressComplete?.Invoke();
        }

        public bool ProcessFile(string strInputFilePath, string strOutputFolderPath)
        {
            return ProcessFile(strInputFilePath, strOutputFolderPath, string.Empty);
        }

        public abstract bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath);

        protected string ReplaceFilenameSuffix(FileInfo originalFile, string strNewSuffix)
        {
            // Appends strNewSuffix to the base name of the original file, then returns a full path to the file using the folder associated with strOriginalFilePath
            // Note that strNewSuffix may contain a file extension though it does not have to
            //  If strNewSuffix does not contain an extension, then the path returned will end in the same extension as strOriginalFilePath

            // Keep track of the original extension on strOriginalFilePath
            var strOriginalExtension = originalFile.Extension;

            // Make sure strNewSuffix is not nothing
            if (strNewSuffix == null)
                strNewSuffix = string.Empty;

            // Obtain the filename, without its extension
            var strNewFileName = Path.GetFileNameWithoutExtension(originalFile.Name);

            // Append strNewSuffix to strNewFileName
            if (Path.HasExtension(strNewSuffix))
            {
                strNewFileName += strNewSuffix;
            }
            else
            {
                strNewFileName += strNewSuffix + strOriginalExtension;
            }

            strNewFileName = Path.Combine(originalFile.DirectoryName, strNewFileName);

            return strNewFileName;
        }

        protected void ReportError(string errMsg, bool throwException = false, Exception innerException = null)
        {
            SetErrorMessage(errMsg);

            if (throwException)
            {
                if (innerException == null)
                {
                    throw new Exception(errMsg);
                }

                throw new Exception(errMsg, innerException);
            }
        }

        protected void ReportMessage(string message)
        {
            OnStatusEvent(message);
        }

        protected void ReportWarning(string message)
        {
            OnWarningEvent(message);
        }

        public bool ResetMassCorrectionTagsAndModificationDefinitions()
        {
            var blnFileNotFound = false;

            // Note: If mMassCorrectionTagsFilePath is blank then the mass correction tags will be reset to the defaults and blnSuccess will be True
            var blnSuccess = mPeptideMods.ReadMassCorrectionTagsFile(MassCorrectionTagsFilePath, ref blnFileNotFound);
            if (!blnSuccess)
            {
                if (blnFileNotFound)
                {
                    SetErrorCode(ePHRPErrorCodes.MassCorrectionTagsFileNotFound);
                }
                else
                {
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingMassCorrectionTagsFile);
                }
            }

            // Note: If mModificationDefinitionsFilePath is blank, then the modifications will be cleared and blnSuccess will be True
            blnSuccess = mPeptideMods.ReadModificationDefinitionsFile(ModificationDefinitionsFilePath, ref blnFileNotFound);
            if (!blnSuccess)
            {
                if (blnFileNotFound)
                {
                    SetErrorCode(ePHRPErrorCodes.ModificationDefinitionFileNotFound);
                    ReportWarning("File not found: " + ModificationDefinitionsFilePath);
                }
                else
                {
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                }
            }

            return blnSuccess;
        }

        protected void ResetProgress()
        {
            ResetProgress(string.Empty);
        }

        protected void ResetProgress(string strProgressStepDescription, bool echoToConsole = false)
        {
            mProgressStepDescription = string.Copy(strProgressStepDescription);
            mProgressPercentComplete = 0;
            ProgressReset?.Invoke();

            if (echoToConsole)
            {
                Console.WriteLine();
                Console.WriteLine();
                Console.WriteLine(ProgressStepDescription);
            }
        }

        protected void SaveModificationSummaryFile(string strModificationSummaryFilePath)
        {
            using (var swOutFile = new StreamWriter(strModificationSummaryFilePath, false))
            {
                // Write the header line
                swOutFile.WriteLine(clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol + SEP_CHAR +
                                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass + SEP_CHAR +
                                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Target_Residues + SEP_CHAR +
                                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Type + SEP_CHAR +
                                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Mass_Correction_Tag + SEP_CHAR +
                                    clsPHRPModSummaryReader.MOD_SUMMARY_COLUMN_Occurrence_Count);

                for (var intIndex = 0; intIndex <= mPeptideMods.ModificationCount - 1; intIndex++)
                {
                    var oModInfo = mPeptideMods.GetModificationByIndex(intIndex);
                    var objSearchResult = oModInfo;
                    if (objSearchResult.OccurrenceCount > 0 || !objSearchResult.UnknownModAutoDefined)
                    {
                        swOutFile.WriteLine(objSearchResult.ModificationSymbol + SEP_CHAR +
                                            objSearchResult.ModificationMass.ToString(CultureInfo.InvariantCulture) + SEP_CHAR +
                                            objSearchResult.TargetResidues + SEP_CHAR +
                                            clsModificationDefinition.ModificationTypeToModificationSymbol(objSearchResult.ModificationType) + SEP_CHAR +
                                            objSearchResult.MassCorrectionTag + SEP_CHAR +
                                            objSearchResult.OccurrenceCount);
                    }
                }
            }
        }

        protected void SaveResultsFileEntrySeqInfo(clsSearchResultsBaseClass searchResult, bool blnUpdateResultToSeqMapFile)
        {
            // Note: Be sure to call Me.InitializeOutputFiles before calling this function
            // blnUpdateResultToSeqMapFile should be set to True only for the first protein of each peptide in each group

            // This ID is assigned using a hashtable containing mPeptideCleanSequence and mPeptideModDescription
            var intUniqueSeqID = mUniqueSequences.GetNextUniqueSequenceID(
                searchResult.PeptideCleanSequence,
                searchResult.PeptideModDescription,
                out var blnExistingSequenceFound);

            if (blnUpdateResultToSeqMapFile)
            {
                // Write a new entry to the ResultToSeqMap file
                mResultToSeqMapFile.WriteLine(searchResult.ResultID + SEP_CHAR + intUniqueSeqID);

                // Only write this entry to the SeqInfo and ModDetails files if blnExistingSequenceFound is False

                if (!blnExistingSequenceFound)
                {
                    // Write a new entry to the SeqInfo file
                    mSeqInfoFile.WriteLine(intUniqueSeqID + SEP_CHAR +
                                           searchResult.SearchResultModificationCount + SEP_CHAR +
                                           searchResult.PeptideModDescription + SEP_CHAR +
                                           searchResult.PeptideMonoisotopicMass.ToString("0.0000000"));

                    if (searchResult.SearchResultModificationCount > 0)
                    {
                        var udtModNameAndResidueLoc = new udtModNameAndResidueLocType[searchResult.SearchResultModificationCount];
                        var intPointerArray = new int[searchResult.SearchResultModificationCount];

                        if (searchResult.SearchResultModificationCount == 1)
                        {
                            intPointerArray[0] = 0;
                        }
                        else
                        {
                            // Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                            for (var intIndex = 0; intIndex <= searchResult.SearchResultModificationCount - 1; intIndex++)
                            {
                                var resultModDetails = searchResult.GetSearchResultModDetailsByIndex(intIndex);
                                udtModNameAndResidueLoc[intIndex].ResidueLocInPeptide = resultModDetails.ResidueLocInPeptide;
                                udtModNameAndResidueLoc[intIndex].ModName = resultModDetails.ModDefinition.MassCorrectionTag;
                                intPointerArray[intIndex] = intIndex;
                            }

                            Array.Sort(udtModNameAndResidueLoc, intPointerArray, new IModNameAndResidueLocComparer());
                        }

                        // Write out the modifications to the ModDetails file
                        // Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
                        for (var intIndex = 0; intIndex <= searchResult.SearchResultModificationCount - 1; intIndex++)
                        {
                            var resultModDetails = searchResult.GetSearchResultModDetailsByIndex(intPointerArray[intIndex]);
                            mModDetailsFile.WriteLine(intUniqueSeqID + SEP_CHAR +
                                                      resultModDetails.ModDefinition.MassCorrectionTag + SEP_CHAR +
                                                      resultModDetails.ResidueLocInPeptide);
                        }
                    }
                }
            }

            // Write a new entry to the SeqToProteinMap file if not yet defined
            if (!CheckSeqToProteinMapDefined(intUniqueSeqID, searchResult.ProteinName))
            {
                mSeqToProteinMapFile.WriteLine(intUniqueSeqID + SEP_CHAR +
                                               Convert.ToInt32(searchResult.PeptideCleavageState) + SEP_CHAR +
                                               Convert.ToInt32(searchResult.PeptideTerminusState) + SEP_CHAR +
                                               searchResult.ProteinName + SEP_CHAR +
                                               searchResult.ProteinExpectationValue + SEP_CHAR +
                                               searchResult.ProteinIntensity);
            }
        }

        protected void SetErrorCode(ePHRPErrorCodes eNewErrorCode)
        {
            SetErrorCode(eNewErrorCode, false);
        }

        protected void SetErrorCode(ePHRPErrorCodes eNewErrorCode, bool blnLeaveExistingErrorCodeUnchanged)
        {
            if (blnLeaveExistingErrorCodeUnchanged && mErrorCode != ePHRPErrorCodes.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mErrorCode = eNewErrorCode;
            }
        }

        protected void SetErrorMessage(string message)
        {
            if (message == null)
                message = string.Empty;

            mErrorMessage = message;
            if (message.Length > 0)
            {
                OnErrorEvent(message);
            }
        }

        /// <summary>
        /// Return the text up to (but not including) the first space in strProteinNameAndDescription
        /// </summary>
        /// <param name="strProteinNameAndDescription"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected virtual string TruncateProteinName(string strProteinNameAndDescription)
        {
            var intIndex = strProteinNameAndDescription.IndexOf(' ');
            if (intIndex > 0)
            {
                return strProteinNameAndDescription.Substring(0, intIndex);
            }

            return strProteinNameAndDescription;
        }

        protected void UpdatePepToProteinMapPeptide(List<udtPepToProteinMappingType> lstPepToProteinMapping, int intIndex, string strPeptide)
        {
            var udtItem = lstPepToProteinMapping[intIndex];
            udtItem.Peptide = strPeptide;
            lstPepToProteinMapping[intIndex] = udtItem;
        }

        protected void UpdateProgress(string strProgressStepDescription)
        {
            UpdateProgress(strProgressStepDescription, mProgressPercentComplete);
        }

        protected void UpdateProgress(float sngPercentComplete)
        {
            UpdateProgress(this.ProgressStepDescription, sngPercentComplete);
        }

        protected void UpdateProgress(string strProgressStepDescription, float sngPercentComplete)
        {
            mProgressStepDescription = string.Copy(strProgressStepDescription);
            if (sngPercentComplete < 0)
            {
                sngPercentComplete = 0;
            }
            else if (sngPercentComplete > 100)
            {
                sngPercentComplete = 100;
            }
            mProgressPercentComplete = sngPercentComplete;

            OnProgressUpdate(ProgressStepDescription, ProgressPercentComplete);
        }

        private bool ValidatePeptideToProteinMapResults(string strPeptideToProteinMapFilePath, bool blnIgnorePeptideToProteinMapperErrors)
        {
            bool blnSuccess;

            var intPeptideCount = 0;
            var intPeptideCountNoMatch = 0;
            var intLinesRead = 0;
            var chSplitChars = new[] {'\t'};

            try
            {
                // Validate that none of the results in strPeptideToProteinMapFilePath has protein name PROTEIN_NAME_NO_MATCH ( __NoMatch__ )

                var strLastPeptide = string.Empty;

                using (var srInFile = new StreamReader(new FileStream(strPeptideToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        var strLineIn = srInFile.ReadLine();
                        intLinesRead += 1;

                        if (intLinesRead <= 1 || string.IsNullOrEmpty(strLineIn))
                            continue;

                        var strSplitLine = strLineIn.Split(chSplitChars, 2);
                        if (strSplitLine.Length <= 0)
                            continue;

                        if (strSplitLine[0] != strLastPeptide)
                        {
                            intPeptideCount += 1;
                            strLastPeptide = string.Copy(strSplitLine[0]);
                        }

                        if (strLineIn.Contains(PROTEIN_NAME_NO_MATCH))
                        {
                            intPeptideCountNoMatch += 1;
                        }
                    }
                }

                if (intPeptideCount == 0)
                {
                    SetErrorMessage("Peptide to protein mapping file is empty: " + strPeptideToProteinMapFilePath);
                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                    blnSuccess = false;
                }
                else if (intPeptideCountNoMatch == 0)
                {
                    blnSuccess = true;
                }
                else
                {
                    // Value between 0 and 100
                    var dblErrorPercent = intPeptideCountNoMatch / (double)intPeptideCount * 100.0;

                    var message = string.Format("{0:0.00}% of the entries ({1} / {2}) in the peptide to protein map file ({3}) " +
                                                "did not match to a protein in the FASTA file ({4})",
                                                dblErrorPercent, intPeptideCountNoMatch, intPeptideCount,
                                                Path.GetFileName(strPeptideToProteinMapFilePath),
                                                Path.GetFileName(FastaFilePath));

                    if (blnIgnorePeptideToProteinMapperErrors || dblErrorPercent < 0.1)
                    {
                        ReportWarning(message);
                        blnSuccess = true;
                    }
                    else
                    {
                        SetErrorMessage(message);
                        blnSuccess = false;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ValidatePeptideToProteinMapResults:" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        protected void ValidatePHRPReaderSupportFiles(string strPHRPDataFilePath, string strOutputFolderPath)
        {

            try
            {
                if (!string.IsNullOrWhiteSpace(strOutputFolderPath))
                {
                    var phrpDataFile = new FileInfo(strPHRPDataFilePath);
                    var outputFolder = new DirectoryInfo(strOutputFolderPath);

                    if (!string.Equals(phrpDataFile.DirectoryName, outputFolder.FullName, StringComparison.CurrentCultureIgnoreCase))
                    {
                        var strMSGFFileName = Path.GetFileName(ReplaceFilenameSuffix(phrpDataFile, FILENAME_SUFFIX_MSGF));

                        var strSourcePath = Path.Combine(phrpDataFile.DirectoryName, strMSGFFileName);
                        var strTargetPath = Path.Combine(outputFolder.FullName, strMSGFFileName);

                        if (File.Exists(strSourcePath) & !File.Exists(strTargetPath))
                        {
                            File.Copy(strSourcePath, strTargetPath);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                ReportWarning("Error in ValidatePHRPReaderSupportFiles: " + ex.Message);
            }
        }

        protected bool ValidateProteinFastaFile(string strFastaFilePath)
        {
            var blnSuccess = ValidateProteinFastaFile(strFastaFilePath, out var strWarningMessage);

            if (!blnSuccess)
            {
                ReportWarning(strWarningMessage);
            }

            return blnSuccess;
        }

        public static bool ValidateProteinFastaFile(string strFastaFilePath, out string strWarningMessage)
        {
            // This RegEx looks for standard amino acids, skipping A, T, C, and G
            var reDefiniteAminoAcid = new Regex(@"[DEFHIKLMNPQRSVWY]", RegexOptions.IgnoreCase | RegexOptions.Compiled);

            // This RegEx looks for A, T, C, and G
            var rePotentialNucleicAcid = new Regex(@"[ATCG]", RegexOptions.IgnoreCase | RegexOptions.Compiled);

            // This matches any letter
            var reLetter = new Regex(@"[A-Z]", RegexOptions.IgnoreCase | RegexOptions.Compiled);

            var intValidProteinCount = 0;
            var intInvalidProteinCount = 0;

            try
            {
                strWarningMessage = string.Empty;

                if (string.IsNullOrEmpty(strFastaFilePath))
                {
                    Console.WriteLine();
                    strWarningMessage = "strFastaFilePath is not defined in ValidateProteinFastaFile";
                    return false;
                }

                if (!File.Exists(strFastaFilePath))
                {
                    Console.WriteLine();
                    strWarningMessage = "Fasta file not found: " + strFastaFilePath;
                    return false;
                }

                var objFastaFile = new ProteinFileReader.FastaFileReader();
                if (!objFastaFile.OpenFile(strFastaFilePath))
                {
                    Console.WriteLine();
                    strWarningMessage = "Error opening the fasta file: " + strFastaFilePath;
                    return false;
                }

                // Read the first 500 proteins and confirm that each contains amino acid residues
                while (objFastaFile.ReadNextProteinEntry())
                {
                    var intDefiniteAminoAcidCount = reDefiniteAminoAcid.Matches(objFastaFile.ProteinSequence).Count;
                    var intPotentialNucleicAcidCount = rePotentialNucleicAcid.Matches(objFastaFile.ProteinSequence).Count;
                    var intLetterCount = reLetter.Matches(objFastaFile.ProteinSequence).Count;

                    if (intDefiniteAminoAcidCount > 0.1 * intLetterCount)
                    {
                        intValidProteinCount += 1;
                    }
                    else if (intPotentialNucleicAcidCount > 0.95 * intLetterCount)
                    {
                        intInvalidProteinCount += 1;
                    }

                    if (intValidProteinCount + intInvalidProteinCount >= 500)
                    {
                        break;
                    }
                }

                if (intValidProteinCount < intInvalidProteinCount)
                {
                    Console.WriteLine();
                    strWarningMessage = "Fasta file contains Nucleic Acids, not Amino Acids: " + Path.GetFileName(strFastaFilePath);
                    return false;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine();
                strWarningMessage = "Exception in ValidateProteinFastaFile: " + ex.Message;
                return false;
            }

            return true;
        }

        #region "PeptideToProteinMapper Event Handlers"

        private void PeptideToProteinMapper_ProgressChanged(string taskDescription, float percentComplete)
        {
            if (percentComplete >= mNextPeptideToProteinMapperLevel)
            {
                mNextPeptideToProteinMapperLevel += 25;
                UpdateProgress(PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE + percentComplete * (PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE - PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE) / 100);
                Console.WriteLine(" PeptideToProteinMapper is " + percentComplete.ToString("0") + "% complete");
            }
        }

        #endregion

        #region "IComparer classes"

        protected class ISearchOptionModificationInfoComparer : IComparer<udtSearchOptionModificationInfoType>
        {
            public int Compare(udtSearchOptionModificationInfoType x, udtSearchOptionModificationInfoType y)
            {
                if (x.SortOrder > y.SortOrder)
                {
                    return 1;
                }

                if (x.SortOrder < y.SortOrder)
                {
                    return -1;
                }

                if (x.ModificationMass > y.ModificationMass)
                {
                    return 1;
                }

                if (x.ModificationMass < y.ModificationMass)
                {
                    return -1;
                }

                return 0;
            }
        }

        internal class IModNameAndResidueLocComparer : IComparer<udtModNameAndResidueLocType>
        {
            public int Compare(udtModNameAndResidueLocType x, udtModNameAndResidueLocType y)
            {
                if (x.ResidueLocInPeptide > y.ResidueLocInPeptide)
                {
                    return 1;
                }

                if (x.ResidueLocInPeptide < y.ResidueLocInPeptide)
                {
                    return -1;
                }

                if (x.ModName == null)
                    x.ModName = string.Empty;

                if (y.ModName == null)
                    y.ModName = string.Empty;

                return string.Compare(x.ModName, y.ModName, StringComparison.Ordinal);
            }
        }

        protected class PepToProteinMappingComparer : IComparer<udtPepToProteinMappingType>
        {
            public int Compare(udtPepToProteinMappingType x, udtPepToProteinMappingType y)
            {
                var result = string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
                if (result == 0)
                {
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }
                return result;
            }
        }

        protected class PepToProteinMappingPeptideSearchComparer : IComparer<udtPepToProteinMappingType>
        {
            public int Compare(udtPepToProteinMappingType x, udtPepToProteinMappingType y)
            {
                return string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
            }
        }

        #endregion
    }
}