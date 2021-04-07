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
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Threading;
using PeptideHitResultsProcessor.Data;
using PeptideToProteinMapEngine;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM.FileProcessor;
using ProteinCoverageSummarizer;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class can be used as a base class for peptide hit results processor classes
    /// </summary>
    public abstract class PHRPBaseClass : PRISM.EventNotifier
    {
        // Ignore Spelling: Da, A-Za-z, Fscore, MaxQuant, prot, mts, MSFragger, xxx

        /// <summary>
        /// Constructor
        /// </summary>
        protected PHRPBaseClass(PHRPOptions options)
        {
            FileDate = "April 6, 2021";

            Options = options;

            mPeptideSeqMassCalculator = new PeptideMassCalculator { ChargeCarrierMass = PeptideMassCalculator.MASS_PROTON };

            // Initialize mPeptideMods
            mPeptideMods = new PeptideModificationContainer();

            // Initialize mUniqueSequences
            mUniqueSequences = new UniqueSequencesContainer();

            // Initialize mSeqToProteinMap
            mSeqToProteinMap = new SortedSet<string>();

            // Define a RegEx to replace all of the non-letter characters
            mReplaceSymbols = new Regex("[^A-Za-z]", RegexOptions.Compiled);

            mProteinNameOrder = new Dictionary<string, int>();

            mErrorCode = PHRPErrorCode.NoError;
            mErrorMessage = string.Empty;
        }

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
        public const string TopPIC_RESULTS_FILE_SUFFIX = "_TopPIC_PrSMs.txt";

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

        private const string PROTEIN_NAME_NO_MATCH = "__NoMatch__";

        protected struct SearchOptionModificationInfo
        {
            public int SortOrder;
            public double ModificationMass;
            public string TargetResidues;
            public ModificationDefinition.ResidueModificationType ModificationType;

            /// <summary>
            /// Duplicate this modification via a deep copy
            /// </summary>
            public SearchOptionModificationInfo Clone()
            {
                var modInfo = new SearchOptionModificationInfo
                {
                    SortOrder = SortOrder,
                    ModificationMass = ModificationMass,
                    TargetResidues = TargetResidues,
                    ModificationType = ModificationType
                };

                return modInfo;
            }

            /// <summary>
            /// Show the modification type and mass
            /// </summary>
            public override string ToString()
            {
                return ModificationType + ": " + ModificationMass + " @ " + TargetResidues;
            }
        }

        internal struct ModNameAndResidueLoc
        {
            public string ModName;
            public int ResidueLocInPeptide;

            public override string ToString()
            {
                return ResidueLocInPeptide + ": " + ModName;
            }
        }

        protected struct PepToProteinMapping
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

        protected PHRPErrorCode mErrorCode;
        protected string mErrorMessage;

        protected readonly PeptideMassCalculator mPeptideSeqMassCalculator;

        protected readonly PeptideModificationContainer mPeptideMods;
        private readonly UniqueSequencesContainer mUniqueSequences;
        private readonly SortedSet<string> mSeqToProteinMap;

        private StreamWriter mResultToSeqMapFile;
        private StreamWriter mSeqInfoFile;
        private StreamWriter mModDetailsFile;
        private StreamWriter mSeqToProteinMapFile;

        private int mNextPeptideToProteinMapperLevel;

        /// <summary>
        /// Tracks the protein names in the order that they are listed in the FASTA file
        /// Keys are protein names, values are a sequentially assigned integer
        /// </summary>
        protected readonly Dictionary<string, int> mProteinNameOrder;

        private readonly Regex mReplaceSymbols;

        public event ProgressResetEventHandler ProgressReset;
        public delegate void ProgressResetEventHandler();

        public event ProgressCompleteEventHandler ProgressComplete;
        public delegate void ProgressCompleteEventHandler();

        protected string mProgressStepDescription = string.Empty;

        /// <summary>
        /// Ranges from 0 to 100, but can contain decimal percentage values
        /// </summary>
        protected float mProgressPercentComplete;

        public bool AbortProcessing { get; set; }

        public PHRPErrorCode ErrorCode => mErrorCode;

        public string ErrorMessage => GetErrorMessage();

        public string FileVersion => GetVersionForExecutingAssembly();

        public string FileDate { get; protected set; }

        public PHRPOptions Options { get; }

        public string ProgressStepDescription => mProgressStepDescription;

        // ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public float ProgressPercentComplete => Convert.ToSingle(Math.Round(mProgressPercentComplete, 2));

        public void AbortProcessingNow()
        {
            AbortProcessing = true;
        }

        public static string AutoDefinePeptideHitResultsFilePath(
            ResultsFileFormat peptideHitResultFileFormat,
            string sourceDirectoryPath,
            string baseName)
        {
            if (string.IsNullOrEmpty(baseName))
                return AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath);

            return peptideHitResultFileFormat switch
            {
                ResultsFileFormat.SequestFirstHitsFile => Path.Combine(sourceDirectoryPath, baseName + SEQUEST_FIRST_HITS_FILE_SUFFIX),
                ResultsFileFormat.SequestSynopsisFile => Path.Combine(sourceDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX),
                ResultsFileFormat.XTandemXMLFile => Path.Combine(sourceDirectoryPath, baseName + XTANDEM_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.InspectTXTFile => Path.Combine(sourceDirectoryPath, baseName + INSPECT_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.MSGFPlusTXTFile => Path.Combine(sourceDirectoryPath, baseName + MSGFDB_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.MSAlignTXTFile => Path.Combine(sourceDirectoryPath, baseName + MSALIGN_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.MODaTXTFile => Path.Combine(sourceDirectoryPath, baseName + MODa_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.MODPlusTXTFile => Path.Combine(sourceDirectoryPath, baseName + MODPlus_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.MSPathFinderTSVFile => Path.Combine(sourceDirectoryPath, baseName + MSPathFinder_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.TopPICTXTFile => Path.Combine(sourceDirectoryPath, baseName + TopPIC_RESULTS_FILE_SUFFIX),
                ResultsFileFormat.MaxQuantTXTFile => Path.Combine(sourceDirectoryPath, MaxQuantResultsProcessor.MSMS_FILE_NAME),
                _ => AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath)
            };
        }

        public static string AutoDefinePeptideHitResultsFilePath(string sourceDirectoryPath)
        {
            // Looks for a file ending in _syn.txt, _fht.txt, _xt.xml, or _inspect.txt in directory sourceDirectoryPath
            // Returns the first matching file found

            var matchSpec = string.Empty;

            try
            {
                for (var index = 0; index <= 3; index++)
                {
                    matchSpec = index switch
                    {
                        0 => "*" + SEQUEST_SYNOPSIS_FILE_SUFFIX,
                        1 => "*" + SEQUEST_FIRST_HITS_FILE_SUFFIX,
                        2 => "*" + XTANDEM_RESULTS_FILE_SUFFIX,
                        3 => "*" + INSPECT_RESULTS_FILE_SUFFIX,
                        _ => matchSpec
                    };

                    var sourceDirectory = new DirectoryInfo(sourceDirectoryPath);
                    foreach (var resultsFile in sourceDirectory.GetFiles(matchSpec))
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
            if (string.IsNullOrWhiteSpace(Options.FastaFilePath))
            {
                // Nothing to do
                return true;
            }

            mProteinNameOrder.Clear();
            var reExtractProteinName = new Regex("^>([^ ]+)", RegexOptions.Compiled);

            ReportMessage("Caching protein names from the FASTA file");

            try
            {
                var proteinNumber = 0;

                using (var reader = new StreamReader(new FileStream(Options.FastaFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
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

                        proteinNumber++;

                        mProteinNameOrder.Add(proteinName, proteinNumber);
                    }
                }

                ReportMessage(string.Format("Cached {0:N0} proteins", mProteinNameOrder.Count));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error caching protein names from fasta file " + Path.GetFileName(Options.FastaFilePath) + ": " + ex.Message);
                return false;
            }
        }

        protected bool CheckSeqToProteinMapDefined(int uniqueSeqID, string proteinName)
        {
            // Returns True if the sequence to protein map was already defined
            // Returns False if the mapping was not defined (will also update mSeqToProteinMap)

            bool existingMapFound;

            try
            {
                if (proteinName == null)
                    proteinName = string.Empty;

                var key = uniqueSeqID + UNIQUE_SEQ_TO_PROTEIN_MAP_SEP + proteinName;

                if (mSeqToProteinMap.Contains(key))
                {
                    existingMapFound = true;
                }
                else
                {
                    mSeqToProteinMap.Add(key);
                    existingMapFound = false;
                }
            }
            catch (Exception)
            {
                existingMapFound = false;
            }

            return existingMapFound;
        }

        protected int CIntSafe(string value, int defaultValue)
        {
            try
            {
                // Note: Integer.Parse() fails if value contains a decimal point, even if it is "8.000"
                // Thus, we're converting to a double first, and then rounding
                return (int)Math.Round(Convert.ToDouble(value));
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        protected double CDblSafe(string value, double defaultValue)
        {
            try
            {
                return double.Parse(value);
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        protected float CSngSafe(string value, float defaultValue)
        {
            try
            {
                return float.Parse(value);
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        /// <summary>
        /// Validate the input file and output directory paths
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <returns>True if successful, False if failure</returns>
        protected bool CleanupFilePaths(ref string inputFilePath, ref string outputDirectoryPath)
        {
            try
            {
                // Make sure inputFilePath points to a valid file
                var inputFile = new FileInfo(inputFilePath);

                if (!inputFile.Exists)
                {
                    SetErrorMessage("Input file not found: " + inputFilePath);
                    if (inputFilePath.Contains(".."))
                    {
                        ReportWarning("Absolute path: " + inputFile.DirectoryName);
                    }
                    SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                    return false;
                }

                if (string.IsNullOrWhiteSpace(outputDirectoryPath))
                {
                    // Define outputDirectoryPath based on inputFilePath
                    outputDirectoryPath = inputFile.DirectoryName;
                }

                // Make sure outputDirectoryPath points to a directory
                DirectoryInfo outputDirectory;

                if (string.IsNullOrWhiteSpace(outputDirectoryPath))
                {
                    outputDirectory = inputFile.Directory;
                }
                else
                {
                    outputDirectory = new DirectoryInfo(outputDirectoryPath);
                }

                if (outputDirectory == null)
                {
                    outputDirectoryPath = ".";
                    return true;
                }

                if (outputDirectory.Exists)
                    return true;

                // outputDirectoryPath points to a non-existent directory; attempt to create it
                try
                {
                    outputDirectory.Create();
                }
                catch (Exception ex2)
                {
                    SetErrorMessage("Invalid output directory: " + outputDirectoryPath, ex2);
                    SetErrorCode(PHRPErrorCode.InvalidOutputDirectoryPath);
                    return false;
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error cleaning up the file paths: " + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.FilePathError);
                return false;
            }
        }

        /// <summary>
        /// Collapses a list of strings to a tab-delimited line of text
        /// </summary>
        /// <param name="fields"></param>
        protected string CollapseList(List<string> fields)
        {
            return string.Join("\t", fields);
        }

        /// <summary>
        /// Examine the extension on filePath to determine the file format
        /// </summary>
        /// <param name="filePath"></param>
        public static ResultsFileFormat DetermineResultsFileFormat(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath))
                return ResultsFileFormat.AutoDetermine;

            var extensionLCase = Path.GetExtension(filePath).ToLower();

            var fileName = Path.GetFileName(filePath);
            var baseFileName = Path.GetFileNameWithoutExtension(filePath);

            if (extensionLCase == ".xml")
            {
                return ResultsFileFormat.XTandemXMLFile;
            }

            if (baseFileName.EndsWith(SequestResultsProcessor.FILENAME_SUFFIX_FIRST_HITS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.SequestFirstHitsFile;
            }

            if (baseFileName.EndsWith(SequestResultsProcessor.FILENAME_SUFFIX_SYNOPSIS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.SequestSynopsisFile;
            }

            if (baseFileName.EndsWith(InSpecTResultsProcessor.FILENAME_SUFFIX_INSPECT_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.InspectTXTFile;
            }

            if (baseFileName.EndsWith(MSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFDB_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MSGFPlusTXTFile;
            }

            if (baseFileName.EndsWith(MSGFPlusResultsProcessor.FILENAME_SUFFIX_MSGFPLUS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MSGFPlusTXTFile;
            }
            if (baseFileName.EndsWith(MSAlignResultsProcessor.FILENAME_SUFFIX_MSALIGN_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MSAlignTXTFile;
            }

            if (baseFileName.EndsWith(MODaResultsProcessor.FILENAME_SUFFIX_MODA_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MODaTXTFile;
            }

            if (baseFileName.EndsWith(MODPlusResultsProcessor.FILENAME_SUFFIX_MODPlus_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MODPlusTXTFile;
            }

            if (baseFileName.EndsWith(MSPathFinderResultsProcessor.FILENAME_SUFFIX_MSPathFinder_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MSPathFinderTSVFile;
            }

            if (baseFileName.EndsWith(TopPICResultsProcessor.FILENAME_SUFFIX_TopPIC_PRSMs_FILE, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.TopPICTXTFile;
            }

            if (fileName.Equals(MaxQuantResultsProcessor.MSMS_FILE_NAME, StringComparison.OrdinalIgnoreCase) ||
                fileName.EndsWith(MaxQuantResultsProcessor.MSMS_FILE_NAME, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MaxQuantTXTFile;
            }
            
            if (extensionLCase == ".tsv")
            {
                // Assume this is an MS-GF+ TSV file
                return ResultsFileFormat.MSGFPlusTXTFile;
            }

            var candidateDirectory = new DirectoryInfo(filePath);
            if (candidateDirectory.Exists && candidateDirectory.GetFiles(MaxQuantResultsProcessor.MSMS_FILE_NAME).Length > 1)
            {
                return ResultsFileFormat.MaxQuantTXTFile;
            }

            // Unknown extension
            return ResultsFileFormat.AutoDetermine;
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

        protected void ComputePseudoPeptideLocInProtein(SearchResultsBaseClass searchResult)
        {
            // Set these to 1 and 10000 since MSGFDB, Sequest, and Inspect results files do not contain protein sequence information
            // If we find later that the peptide sequence spans the length of the protein, we'll revise .ProteinSeqResidueNumberEnd as needed
            searchResult.ProteinSeqResidueNumberStart = 1;
            searchResult.ProteinSeqResidueNumberEnd = 10000;

            if (searchResult.PeptidePreResidues.Trim().EndsWith(PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST.ToString()))
            {
                // The peptide is at the N-Terminus of the protein
                searchResult.PeptideLocInProteinStart = searchResult.ProteinSeqResidueNumberStart;
                searchResult.PeptideLocInProteinEnd = searchResult.PeptideLocInProteinStart + searchResult.PeptideCleanSequence.Length - 1;

                if (searchResult.PeptidePostResidues.Trim()[0] == PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST)
                {
                    // The peptide spans the entire length of the protein
                    searchResult.ProteinSeqResidueNumberEnd = searchResult.PeptideLocInProteinEnd;
                }
                else
                {
                    if (searchResult.PeptideLocInProteinEnd > searchResult.ProteinSeqResidueNumberEnd)
                    {
                        // The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                        searchResult.ProteinSeqResidueNumberEnd = searchResult.PeptideLocInProteinEnd + 1;
                    }
                }
            }
            else if (searchResult.PeptidePostResidues.Trim().StartsWith(PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST.ToString()))
            {
                // The peptide is at the C-Terminus of the protein
                searchResult.PeptideLocInProteinEnd = searchResult.ProteinSeqResidueNumberEnd;
                searchResult.PeptideLocInProteinStart = searchResult.PeptideLocInProteinEnd - searchResult.PeptideCleanSequence.Length + 1;

                if (searchResult.PeptideLocInProteinStart < searchResult.ProteinSeqResidueNumberStart)
                {
                    // The peptide is more than 10000 characters long; this is highly unlikely
                    searchResult.ProteinSeqResidueNumberEnd = searchResult.ProteinSeqResidueNumberStart + 1 + searchResult.PeptideCleanSequence.Length;
                    searchResult.PeptideLocInProteinEnd = searchResult.ProteinSeqResidueNumberEnd;
                    searchResult.PeptideLocInProteinStart = searchResult.PeptideLocInProteinEnd - searchResult.PeptideCleanSequence.Length + 1;
                }
            }
            else
            {
                searchResult.PeptideLocInProteinStart = searchResult.ProteinSeqResidueNumberStart + 1;
                searchResult.PeptideLocInProteinEnd = searchResult.PeptideLocInProteinStart + searchResult.PeptideCleanSequence.Length - 1;

                if (searchResult.PeptideLocInProteinEnd > searchResult.ProteinSeqResidueNumberEnd)
                {
                    // The peptide is more than 10000 characters long; this is highly unlikely, but we'll update .ProteinSeqResidueNumberEnd as needed
                    searchResult.ProteinSeqResidueNumberEnd = searchResult.PeptideLocInProteinEnd + 1;
                }
            }
        }

        protected virtual string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            return ConstructPepToProteinMapFilePathWork(inputFilePath, outputDirectoryPath, mts);
        }

        private string ConstructPepToProteinMapFilePathWork(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var pepToProteinMapBaseName = Path.GetFileNameWithoutExtension(inputFilePath);

            string pepToProteinMapFileName;
            if (mts)
            {
                pepToProteinMapFileName = pepToProteinMapBaseName + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + "MTS.txt";
            }
            else
            {
                pepToProteinMapFileName = pepToProteinMapBaseName + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + ".txt";
            }

            var inputFile = new FileInfo(inputFilePath);
            var pepToProteinMapFilePath = string.IsNullOrWhiteSpace(inputFile.DirectoryName) ?
                                              pepToProteinMapFileName :
                                              Path.Combine(inputFile.DirectoryName, pepToProteinMapFileName);

            var pepToProteinMapFile = new FileInfo(pepToProteinMapFilePath);
            if (string.IsNullOrWhiteSpace(outputDirectoryPath))
            {
                return pepToProteinMapFile.FullName;
            }

            var alternatePepToProteinMapFilePath = Path.Combine(outputDirectoryPath, pepToProteinMapFileName);
            var alternatePepToProteinMapFile = new FileInfo(alternatePepToProteinMapFilePath);

            if (pepToProteinMapFile.Exists && !alternatePepToProteinMapFile.Exists)
            {
                return pepToProteinMapFile.FullName;
            }

            if (!pepToProteinMapFile.Exists && alternatePepToProteinMapFile.Exists)
            {
                return alternatePepToProteinMapFile.FullName;
            }

            if (pepToProteinMapFile.Exists && alternatePepToProteinMapFile.Exists &&
                alternatePepToProteinMapFile.LastWriteTime > pepToProteinMapFile.LastWriteTime)
            {
                return alternatePepToProteinMapFile.FullName;
            }

            return pepToProteinMapFile.FullName;
        }

        protected string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts, List<string> suffixesToFind, int charsToRemove)
        {
            var baseName = Path.GetFileNameWithoutExtension(inputFilePath);
            if (string.IsNullOrEmpty(baseName))
                return string.Empty;

            foreach (var item in suffixesToFind)
            {
                if (!baseName.EndsWith(item, StringComparison.OrdinalIgnoreCase))
                    continue;

                // baseName matches something like Dataset_msgfplus_syn
                baseName = baseName.Substring(0, baseName.Length - charsToRemove);
                break;
            }

            var inputDirectoryPath = Path.GetDirectoryName(inputFilePath);

            if (inputDirectoryPath == null)
            {
                return ConstructPepToProteinMapFilePathWork(baseName, outputDirectoryPath, mts);
            }

            return ConstructPepToProteinMapFilePathWork(Path.Combine(inputDirectoryPath, baseName), outputDirectoryPath, mts);
        }

        /// <summary>
        /// Use the PeptideToProteinMapEngine to create the Peptide to Protein map file for the file or files in sourcePHRPDataFiles
        /// </summary>
        /// <param name="sourcePHRPDataFiles"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        protected bool CreatePepToProteinMapFile(List<string> sourcePHRPDataFiles, string mtsPepToProteinMapFilePath)
        {
            var success = false;

            try
            {
                if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because mtsPepToProteinMapFilePath is empty; likely a programming bug");
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                if (string.IsNullOrEmpty(Options.FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File Path is not defined");
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                if (!File.Exists(Options.FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the Fasta File was not found: " + Options.FastaFilePath);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                // Verify that the fasta file is not a DNA-sequence based fasta file
                success = ValidateProteinFastaFile(Options.FastaFilePath);
                if (!success)
                {
                    return false;
                }

                Console.WriteLine();
                Console.WriteLine();
                UpdateProgress("Creating Peptide to Protein Map file", PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE);

                // Initialize items
                var mtsPepToProteinMapFile = new FileInfo(mtsPepToProteinMapFilePath);
                string outputDirectoryPath;

                if (string.IsNullOrWhiteSpace(mtsPepToProteinMapFile.DirectoryName))
                    outputDirectoryPath = string.Empty;
                else
                    outputDirectoryPath = mtsPepToProteinMapFile.DirectoryName;

                var peptideToProteinMapResults = new SortedSet<string>();

                var options = new ProteinCoverageSummarizerOptions()
                {
                    IgnoreILDifferences = false,
                    MatchPeptidePrefixAndSuffixToProtein = false,
                    OutputProteinSequence = false,
                    PeptideFileSkipFirstLine = false,
                    RemoveSymbolCharacters = true,
                    ProteinInputFilePath = Options.FastaFilePath,
                    SaveProteinToPeptideMappingFile = true,
                    SearchAllProteinsForPeptideSequence = true,
                    SearchAllProteinsSkipCoverageComputationSteps = true
                };

                var peptideToProteinMapper = new clsPeptideToProteinMapEngine(options)
                {
                    DeleteTempFiles = true,
                    InspectParameterFilePath = string.Empty,
                    LogMessagesToFile = false,
                    PeptideInputFileFormat = clsPeptideToProteinMapEngine.PeptideInputFileFormatConstants.PHRPFile
                };

                RegisterEvents(peptideToProteinMapper);

                // Handle progress updates using PeptideToProteinMapper_ProgressChanged instead of OnProgressUpdate
                peptideToProteinMapper.ProgressUpdate -= OnProgressUpdate;
                peptideToProteinMapper.ProgressUpdate += PeptideToProteinMapper_ProgressChanged;
                peptideToProteinMapper.SkipConsoleWriteIfNoProgressListener = true;

                using (var writer = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    foreach (var inputFilePath in sourcePHRPDataFiles)
                    {
                        var resultsFilePath = Path.GetFileNameWithoutExtension(inputFilePath) +
                            clsPeptideToProteinMapEngine.FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING;

                        resultsFilePath = Path.Combine(outputDirectoryPath, resultsFilePath);

                        // Make sure the results file doesn't already exist
                        DeleteFileIgnoreErrors(resultsFilePath);

                        peptideToProteinMapper.ProgressUpdate += PeptideToProteinMapper_ProgressChanged;
                        mNextPeptideToProteinMapperLevel = 25;

                        success = peptideToProteinMapper.ProcessFile(inputFilePath, outputDirectoryPath, string.Empty, true);

                        peptideToProteinMapper.ProgressUpdate -= PeptideToProteinMapper_ProgressChanged;

                        if (success)
                        {
                            if (!File.Exists(resultsFilePath))
                            {
                                SetErrorMessage("Peptide to protein mapping file was not created for " + inputFilePath);
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                success = false;
                                break;
                            }
                            success = ValidatePeptideToProteinMapResults(resultsFilePath, Options.IgnorePeptideToProteinMapperErrors);
                        }
                        else
                        {
                            if (string.IsNullOrWhiteSpace(peptideToProteinMapper.GetErrorMessage()) && peptideToProteinMapper.StatusMessage.IndexOf("error", StringComparison.OrdinalIgnoreCase) >= 0)
                            {
                                SetErrorMessage("Error running clsPeptideToProteinMapEngine: " + peptideToProteinMapper.StatusMessage);
                            }
                            else
                            {
                                if (peptideToProteinMapper.StatusMessage.Length > 0)
                                {
                                    SetErrorMessage("clsPeptideToProteinMapEngine status: " + peptideToProteinMapper.StatusMessage);
                                }
                                SetErrorMessage("Error running clsPeptideToProteinMapEngine: " + peptideToProteinMapper.GetErrorMessage());
                            }

                            if (Options.IgnorePeptideToProteinMapperErrors)
                            {
                                ReportWarning("Ignoring protein mapping error since 'IgnorePeptideToProteinMapperErrors' = True");

                                if (File.Exists(resultsFilePath))
                                {
                                    success = ValidatePeptideToProteinMapResults(resultsFilePath, Options.IgnorePeptideToProteinMapperErrors);
                                }
                                else
                                {
                                    mErrorMessage = string.Empty;
                                    mErrorCode = PHRPErrorCode.NoError;
                                    success = true;
                                }
                            }
                            else
                            {
                                SetErrorMessage("Error in CreatePepToProteinMapFile");
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                            }
                        }

                        if (!File.Exists(resultsFilePath))
                        {
                            continue;
                        }

                        // Read the newly created file and append new entries to mtsPepToProteinMapFilePath
                        using (var reader = new StreamReader(new FileStream(resultsFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                        {
                            while (!reader.EndOfStream)
                            {
                                var lineIn = reader.ReadLine();

                                if (string.IsNullOrWhiteSpace(lineIn))
                                    continue;

                                var splitLine = lineIn.Split(new[] { '\t' }, 2);
                                if (splitLine.Length < 2)
                                    continue;

                                var peptideAndProteinKey = splitLine[0] + "_" + splitLine[1];

                                if (!peptideToProteinMapResults.Contains(peptideAndProteinKey))
                                {
                                    peptideToProteinMapResults.Add(peptideAndProteinKey);
                                    writer.WriteLine(lineIn);
                                }
                            }
                        }

                        // Delete the interim results file
                        DeleteFileIgnoreErrors(resultsFilePath);
                    }
                }

                peptideToProteinMapper.CloseLogFileNow();
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreatePepToProteinMapFile: " + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
            }

            return success;
        }

        /// <summary>
        /// Create the protein mod details file for the specified PHRP data file
        /// </summary>
        /// <param name="phrpDataFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        public bool CreateProteinModDetailsFile(string phrpDataFilePath, string outputDirectoryPath)
        {
            var success = false;

            try
            {
                var inputFile = new FileInfo(phrpDataFilePath);

                var sourcePHRPDataFiles = new List<string> {
                    inputFile.FullName
                };

                var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, outputDirectoryPath, mts: true);

                success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);

                if (success)
                {
                    success = CreateProteinModDetailsFile(phrpDataFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath, PeptideHitResultTypes.Unknown);
                }
            }
            catch (Exception ex)
            {
                ReportWarning("Error in CreateProteinModDetailsFile: " + ex.Message);
            }

            return success;
        }

        public bool CreateProteinModDetailsFile(
            string phrpDataFilePath,
            string outputDirectoryPath,
            string mtsPepToProteinMapFilePath,
            PeptideHitResultTypes phrpResultType)
        {
            try
            {
                Console.WriteLine();

                var progressAtStart = mProgressPercentComplete;

                ResetProgress(mProgressStepDescription);
                UpdateProgress("Creating the Protein Mod Details file", PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE);

                // Confirm that the PHRP data file exists
                var phrpDataFile = new FileInfo(phrpDataFilePath);
                if (!phrpDataFile.Exists)
                {
                    SetErrorMessage("PHRP data file not found in CreateProteinModDetailsFile: " + phrpDataFilePath);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                // Confirm that the _PepToProtMapMTS.txt file exists
                if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath))
                {
                    SetErrorMessage("Cannot create the ProteinMods file because mtsPepToProteinMapFilePath is empty");
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                // Initialize pepToProteinMapping
                var pepToProteinMapping = new List<PepToProteinMapping>();

                // Read the _PepToProtMapMTS file
                var success = LoadPeptideToProteinMapInfo(mtsPepToProteinMapFilePath, pepToProteinMapping, out _);
                if (!success)
                {
                    return false;
                }

                // Assure that pepToProteinMapping is sorted on peptide
                if (pepToProteinMapping.Count > 1)
                {
                    pepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                var proteinModsFilePath = ReplaceFilenameSuffix(phrpDataFile, FILENAME_SUFFIX_PROTEIN_MODS);
                if (!string.IsNullOrEmpty(outputDirectoryPath))
                {
                    var proteinModsFile = Path.GetFileName(proteinModsFilePath);
                    if (string.IsNullOrEmpty(proteinModsFile))
                    {
                        proteinModsFile = Path.GetFileNameWithoutExtension(phrpDataFile.Name) + FILENAME_SUFFIX_PROTEIN_MODS;
                    }
                    proteinModsFilePath = Path.Combine(outputDirectoryPath, proteinModsFile);
                }

                var psmCount = 0;
                var psmCountSkippedSinceReversedOrScrambledProtein = 0;

                // Create a ProteinMods file parallel to the PHRP file
                using var writer = new StreamWriter(new FileStream(proteinModsFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                // Write the header line
                writer.WriteLine(COLUMN_NAME_RESULTID + "\t" +
                                 COLUMN_NAME_PEPTIDE + "\t" +
                                 COLUMN_NAME_UNIQUE_SEQ_ID + "\t" +
                                 COLUMN_NAME_PROTEIN_NAME + "\t" +
                                 COLUMN_NAME_RESIDUE + "\t" +
                                 COLUMN_NAME_PROTEIN_RESIDUE_NUMBER + "\t" +
                                 COLUMN_NAME_RESIDUE_MOD_NAME + "\t" +
                                 COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER + "\t" +
                                 COLUMN_NAME_MSGF_SPECPROB);

                var loadMSGFResults = phrpResultType != PeptideHitResultTypes.MSGFPlus;

                // Update the Mass Calculator to use the one tracked by this class
                // (since this class's calculator knows about custom amino acids and custom charge carriers)
                var startupOptions = new StartupOptions
                {
                    LoadModsAndSeqInfo = true,
                    LoadMSGFResults = loadMSGFResults,
                    LoadScanStatsData = false,
                    PeptideMassCalculator = mPeptideSeqMassCalculator
                };

                using var reader = new ReaderFactory(phrpDataFilePath, phrpResultType, startupOptions)
                {
                    EchoMessagesToConsole = false,
                    SkipDuplicatePSMs = true
                };

                foreach (var errorMessage in reader.ErrorMessages)
                {
                    OnErrorEvent(errorMessage);
                }
                RegisterEvents(reader);

                foreach (var warningMessage in reader.WarningMessages)
                {
                    var msg = warningMessage;
                    if (warningMessage.StartsWith("MSGF file not found", StringComparison.OrdinalIgnoreCase))
                    {
                        msg = "MSGF file not found; column " + COLUMN_NAME_MSGF_SPECPROB + " will not have any data";
                    }
                    ReportWarning(msg);
                }

                reader.ClearErrors();
                reader.ClearWarnings();

                var peptidesNotFoundInPepToProtMapping = 0;
                while (reader.MoveNext())
                {
                    // Use binary search to find this peptide in pepToProteinMapping
                    var pepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(pepToProteinMapping, reader.CurrentPSM.Peptide);

                    if (pepToProteinMapIndex >= 0)
                    {
                        do
                        {
                            psmCount++;

                            var skipProtein = false;
                            if (!Options.ProteinModsFileIncludesReversedProteins)
                            {
                                skipProtein = IsReversedProtein(pepToProteinMapping[pepToProteinMapIndex].Protein);
                                if (skipProtein)
                                {
                                    psmCountSkippedSinceReversedOrScrambledProtein++;
                                }
                            }

                            if (!skipProtein)
                            {
                                WriteModDetailsEntry(reader,
                                    writer,
                                    pepToProteinMapping,
                                    pepToProteinMapIndex,
                                    ref psmCountSkippedSinceReversedOrScrambledProtein);
                            }

                            pepToProteinMapIndex++;
                        } while (pepToProteinMapIndex < pepToProteinMapping.Count &&
                                 reader.CurrentPSM.Peptide == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                    }
                    else
                    {
                        peptidesNotFoundInPepToProtMapping++;
                        ShowPeriodicWarning(peptidesNotFoundInPepToProtMapping, 10,
                            "Peptide not found in pepToProteinMapping: " + reader.CurrentPSM.Peptide);
                    }

                    var overallProgress = ProcessFilesOrDirectoriesBase.ComputeIncrementalProgress(
                        progressAtStart, 100, reader.PercentComplete);

                    UpdateProgress(overallProgress);
                }

                if (psmCount > 0)
                {
                    if (psmCountSkippedSinceReversedOrScrambledProtein == psmCount)
                    {
                        Console.WriteLine();
                        ReportWarning("All PSMs map to reversed or scrambled proteins; the _ProteinMods.txt file is empty");
                    }
                    else if (psmCountSkippedSinceReversedOrScrambledProtein > 0)
                    {
                        Console.WriteLine();
                        Console.WriteLine("Note: skipped {0:N0} / {1:N0} PSMs that map to reversed or scrambled proteins " +
                                          "while creating the _ProteinMods.txt file",
                            psmCountSkippedSinceReversedOrScrambledProtein, psmCount);
                    }
                }

                if (peptidesNotFoundInPepToProtMapping > 10)
                {
                    Console.WriteLine("Note: {0:N0} peptides were not found in pepToProteinMapping", peptidesNotFoundInPepToProtMapping);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateProteinModDetailsFile:" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        protected void DeleteFileIgnoreErrors(string filePath)
        {
            try
            {
                if (File.Exists(filePath))
                {
                    Thread.Sleep(200);
                    File.Delete(filePath);
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        protected void ExpandListIfRequired<T>(List<T> items, int countToAdd, int largeListThreshold = 1000000)
        {
            if (items.Count > largeListThreshold && items.Count + countToAdd > items.Capacity)
            {
                // .NET by default will double the size of the list to accommodate these new items
                // Instead, expand the list by 20% of the current size
                items.Capacity += Convert.ToInt32(items.Count / 5);
            }
        }

        private readonly IComparer<PepToProteinMapping> peptideSearchComparer = new PepToProteinMappingPeptideSearchComparer();

        protected int FindFirstMatchInPepToProteinMapping(List<PepToProteinMapping> pepToProteinMapping, string peptideToFind)
        {
            // Use binary search to find this peptide in pepToProteinMapping
            var udtItemToFind = new PepToProteinMapping
            {
                Peptide = peptideToFind
            };

            var pepToProteinMapIndex = pepToProteinMapping.BinarySearch(udtItemToFind, peptideSearchComparer);

            if (pepToProteinMapIndex > 0)
            {
                // Step Backward until the first match is found
                while (pepToProteinMapIndex > 0 && pepToProteinMapping[pepToProteinMapIndex - 1].Peptide == peptideToFind)
                {
                    pepToProteinMapIndex--;
                }
            }

            return pepToProteinMapIndex;
        }

        protected string GetCleanSequence(string sequenceWithMods)
        {
            return GetCleanSequence(sequenceWithMods, out _, out _);
        }

        protected string GetCleanSequence(string sequenceWithMods, out string prefix, out string suffix)
        {
            if (PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out prefix, out suffix))
            {
                // Remove any non-letter characters
                primarySequence = mReplaceSymbols.Replace(primarySequence, string.Empty);
            }
            else
            {
                // Sequence does not have prefix or suffix letters; use sequenceWithMods
                primarySequence = mReplaceSymbols.Replace(sequenceWithMods, string.Empty);
            }

            return primarySequence;
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to string.Empty
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out string value)
        {
            return GetColumnValue(splitLine, columnIndex, out value, string.Empty);
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to 0
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out int value)
        {
            return GetColumnValue(splitLine, columnIndex, out value, 0);
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out string value, string valueIfMissing)
        {
            if (columnIndex >= 0 && columnIndex < splitLine.Length)
            {
                value = string.Copy(splitLine[columnIndex]);
                return true;
            }

            value = string.Copy(valueIfMissing);
            return false;
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        protected bool GetColumnValue(string[] splitLine, int columnIndex, out int value, int valueIfMissing)
        {
            if (GetColumnValue(splitLine, columnIndex, out var valueText, valueIfMissing.ToString()))
            {
                if (int.TryParse(valueText, out value))
                {
                    return true;
                }

                value = valueIfMissing;
                return false;
            }

            value = valueIfMissing;
            return false;
        }

        /// <summary>
        /// Get the error message, or an empty string if no error
        /// </summary>
        protected string GetErrorMessage()
        {
            var message = ErrorCode switch
            {
                PHRPErrorCode.NoError => string.Empty,
                PHRPErrorCode.InvalidInputFilePath => "Invalid input file path",
                PHRPErrorCode.InvalidOutputDirectoryPath => "Invalid output directory path",
                PHRPErrorCode.ParameterFileNotFound => "Parameter file not found",
                PHRPErrorCode.MassCorrectionTagsFileNotFound => "Mass correction tags file not found",
                PHRPErrorCode.ModificationDefinitionFileNotFound => "Modification definition file not found",
                PHRPErrorCode.ErrorReadingInputFile => "Error reading input file",
                PHRPErrorCode.ErrorCreatingOutputFiles => "Error creating output files",
                PHRPErrorCode.ErrorReadingParameterFile => "Invalid parameter file",
                PHRPErrorCode.ErrorReadingMassCorrectionTagsFile => "Error reading mass correction tags file",
                PHRPErrorCode.ErrorReadingModificationDefinitionsFile => "Error reading modification definitions file",
                PHRPErrorCode.FilePathError => "General file path error",
                PHRPErrorCode.UnspecifiedError => "Unspecified error",
                _ => "Unknown error state"
            };

            if (mErrorMessage.Length > 0)
            {
                if (message.Length > 0)
                {
                    message += "; ";
                }
                message += mErrorMessage;
            }

            return message;
        }

        private string GetVersionForExecutingAssembly()
        {
            string version;

            try
            {
                version = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            }
            catch (Exception)
            {
                version = "??.??.??.??";
            }

            return version;
        }

        /// <summary>
        /// Initializes the StreamWriter objects using baseOutputFilePath as a base name and replacing the suffix with the default suffix names
        /// </summary>
        /// <param name="baseOutputFilePath"></param>
        /// <returns>True if success; does not catch errors; they will be thrown to the calling function if they occur</returns>
        protected bool InitializeSequenceOutputFiles(string baseOutputFilePath)
        {
            var outputFileInfo = new FileInfo(baseOutputFilePath);

            // Initialize the file paths based on baseOutputFilePath
            var resultToSeqMapFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_RESULT_TO_SEQ_MAP);
            var seqInfoFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_SEQ_INFO);
            var modDetailsFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_MOD_DETAILS);
            var seqToProteinMapFilePath = ReplaceFilenameSuffix(outputFileInfo, FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP);

            // Clear the unique sequences container
            mUniqueSequences.Clear();

            // Clear the sequence to protein map
            mSeqToProteinMap.Clear();

            var resultToSeqMapHeaders = new List<string>
            {
                "Result_ID",
                COLUMN_NAME_UNIQUE_SEQ_ID
            };

            var seqInfoHeaders = new List<string>
            {
                COLUMN_NAME_UNIQUE_SEQ_ID,
                "Mod_Count",
                "Mod_Description",
                "Monoisotopic_Mass"
            };

            var modDetailsHeaders = new List<string>
            {
                COLUMN_NAME_UNIQUE_SEQ_ID,
                "Mass_Correction_Tag",
                "Position"
            };

            var seqToProteinMapHeaders = new List<string>
            {
                COLUMN_NAME_UNIQUE_SEQ_ID,
                "Cleavage_State",
                "Terminus_State",
                COLUMN_NAME_PROTEIN_NAME,
                "Protein_Expectation_Value_Log(e)",
                "Protein_Intensity_Log(I)"
            };

            // Initialize the ResultToSeqMap file
            mResultToSeqMapFile = new StreamWriter(resultToSeqMapFilePath);
            mResultToSeqMapFile.WriteLine(CollapseList(resultToSeqMapHeaders));

            // Initialize the SeqInfo file
            mSeqInfoFile = new StreamWriter(seqInfoFilePath, false);
            mSeqInfoFile.WriteLine(CollapseList(seqInfoHeaders));

            // Initialize the ModDetails file
            mModDetailsFile = new StreamWriter(modDetailsFilePath);
            mModDetailsFile.WriteLine(CollapseList(modDetailsHeaders));

            // Initialize the SeqToProtein map file
            mSeqToProteinMapFile = new StreamWriter(seqToProteinMapFilePath, false);
            mSeqToProteinMapFile.WriteLine(CollapseList(seqToProteinMapHeaders));

            return true;
        }

        private static readonly Regex RegexIsLetter = new("[A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Returns true if the character is a letter between A and Z or a and z
        /// </summary>
        /// <param name="chChar">Character to examine</param>
        /// <remarks>
        /// Note that the Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters
        /// In contrast, this method only returns True for normal letters between A and Z (case insensitive)
        /// </remarks>
        public static bool IsLetterAtoZ(char chChar)
        {
            return RegexIsLetter.IsMatch(chChar.ToString());
        }

        protected bool IsReversedProtein(string proteinName)
        {
            if (proteinName.StartsWith("reversed_", StringComparison.OrdinalIgnoreCase))
            {
                // Used in DMS-generated protein collections
                return true;
            }

            if (proteinName.StartsWith("REV_", StringComparison.OrdinalIgnoreCase))
            {
                // Used by MSGFDB
                return true;
            }

            if (proteinName.StartsWith("scrambled_", StringComparison.OrdinalIgnoreCase))
            {
                // Used in DMS-generated protein collections
                return true;
            }
            if (proteinName.StartsWith("xxx_", StringComparison.OrdinalIgnoreCase))
            {
                // Used by MS-GF+ and MSFragger
                return true;
            }

            if (proteinName.StartsWith("xxx.", StringComparison.OrdinalIgnoreCase))
            {
                // Used by Inspect
                return true;
            }

            if (proteinName.EndsWith(":reversed", StringComparison.OrdinalIgnoreCase))
            {
                // Used by X!Tandem
                return true;
            }

            return false;
        }

        protected virtual bool LoadParameterFileSettings(string parameterFilePath)
        {
            const string OPTIONS_SECTION = "PeptideHitResultsProcessorOptions";

            var settingsFile = new PRISM.XmlSettingsFileAccessor();

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
                    var appDirPath = ProcessFilesOrDirectoriesBase.GetAppDirectoryPath();
                    if (string.IsNullOrWhiteSpace(appDirPath))
                    {
                        SetErrorCode(PHRPErrorCode.ParameterFileNotFound);
                        return false;
                    }

                    parameterFilePath = Path.Combine(appDirPath, Path.GetFileName(parameterFilePath));
                    if (!File.Exists(parameterFilePath))
                    {
                        SetErrorCode(PHRPErrorCode.ParameterFileNotFound);
                        return false;
                    }
                }

                if (settingsFile.LoadSettings(parameterFilePath))
                {
                    if (!settingsFile.SectionPresent(OPTIONS_SECTION))
                    {
                        // Section OPTIONS_SECTION was not found in the parameter file; warn the user if mWarnMissingParameterFileSection = True
                        if (Options.WarnMissingParameterFileSection)
                        {
                            SetErrorMessage("The node '<section name=\"" + OPTIONS_SECTION + "\"> was not found in the parameter file: " + parameterFilePath);
                            SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile);
                            return false;
                        }
                        return true;
                    }

                    Options.MassCorrectionTagsFilePath = settingsFile.GetParam(OPTIONS_SECTION, "MassCorrectionTagsFilePath", Options.MassCorrectionTagsFilePath);
                    Options.ModificationDefinitionsFilePath = settingsFile.GetParam(OPTIONS_SECTION, "ModificationDefinitionsFilePath", Options.ModificationDefinitionsFilePath);
                    Options.SearchToolParameterFilePath = settingsFile.GetParam(OPTIONS_SECTION, "SearchToolParameterFilePath", Options.SearchToolParameterFilePath);

                    Options.CreateModificationSummaryFile = settingsFile.GetParam(OPTIONS_SECTION, "CreateModificationSummaryFile", Options.CreateModificationSummaryFile);

                    Options.CreateProteinModsFile = settingsFile.GetParam(OPTIONS_SECTION, "CreateProteinModsFile", Options.CreateProteinModsFile);
                    Options.FastaFilePath = settingsFile.GetParam(OPTIONS_SECTION, "FastaFilePath", Options.FastaFilePath);
                    Options.ProteinModsFileIncludesReversedProteins = settingsFile.GetParam(OPTIONS_SECTION, "ProteinModsFileIncludesReversedProteins", Options.ProteinModsFileIncludesReversedProteins);
                    Options.UseExistingMTSPepToProteinMapFile = settingsFile.GetParam(OPTIONS_SECTION, "UseExistingMTSPepToProteinMapFile", Options.UseExistingMTSPepToProteinMapFile);

                    var leftResidueRegEx = string.Copy(Options.EnzymeMatchSpec.LeftResidueRegEx);
                    var rightResidueRegEx = string.Copy(Options.EnzymeMatchSpec.RightResidueRegEx);

                    leftResidueRegEx = settingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecLeftResidue", leftResidueRegEx, out var valueNotPresent);
                    if (!valueNotPresent)
                    {
                        rightResidueRegEx = settingsFile.GetParam(OPTIONS_SECTION, "EnzymeMatchSpecRightResidue", rightResidueRegEx, out valueNotPresent);

                        if (!valueNotPresent)
                        {
                            Options.EnzymeMatchSpec = new PeptideCleavageStateCalculator.EnzymeMatchSpecInfo(leftResidueRegEx, rightResidueRegEx);
                        }
                    }

                    Options.PeptideNTerminusMassChange = settingsFile.GetParam(OPTIONS_SECTION, "PeptideNTerminusMassChange", Options.PeptideNTerminusMassChange);
                    Options.PeptideCTerminusMassChange = settingsFile.GetParam(OPTIONS_SECTION, "PeptideCTerminusMassChange", Options.PeptideCTerminusMassChange);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadParameterFileSettings:" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Load the PeptideToProteinMap information
        /// </summary>
        /// <param name="pepToProteinMapFilePath">File to read</param>
        /// <param name="pepToProteinMapping">Output parameter: peptide to protein mapping (calling function must pre-initialize the list)</param>
        /// <param name="headerLine">Output parameter: Header line text</param>
        /// <returns>True if successful, false if an error</returns>
        protected bool LoadPeptideToProteinMapInfo(
            string pepToProteinMapFilePath,
            List<PepToProteinMapping> pepToProteinMapping,
            out string headerLine)
        {
            headerLine = string.Empty;

            bool success;

            try
            {
                // Initialize the output parameters
                pepToProteinMapping.Clear();
                headerLine = string.Empty;

                if (string.IsNullOrWhiteSpace(pepToProteinMapFilePath))
                {
                    SetErrorMessage("Warning: PepToProteinMap file is not defined");
                    SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                    return false;
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    SetErrorMessage("Warning: PepToProteinMap file does not exist: " + pepToProteinMapFilePath);
                    SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                    return false;
                }

                // Open the peptide to protein map file for reading
                using (var reader = new StreamReader(new FileStream(pepToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    var linesRead = 0;
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(lineIn))
                            continue;

                        var dataLine = lineIn.Trim();
                        if (dataLine.Length == 0)
                            continue;

                        linesRead++;

                        // Split the line on tabs
                        var splitLine = dataLine.TrimEnd().Split('\t');

                        if (splitLine.Length >= 4)
                        {
                            if (linesRead == 1 && !int.TryParse(splitLine[2], out _))
                            {
                                // Header line; cache it
                                headerLine = string.Copy(dataLine);
                            }
                            else
                            {
                                var pepToProteinMappingEntry = new PepToProteinMapping
                                {
                                    Peptide = string.Copy(splitLine[0]),
                                    Protein = string.Copy(splitLine[1])
                                };
                                int.TryParse(splitLine[2], out pepToProteinMappingEntry.ResidueStart);
                                int.TryParse(splitLine[3], out pepToProteinMappingEntry.ResidueEnd);

                                ExpandListIfRequired(pepToProteinMapping, 1);

                                pepToProteinMapping.Add(pepToProteinMappingEntry);
                            }
                        }
                    }
                }

                success = true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the Peptide to Protein Map File (" + Path.GetFileName(pepToProteinMapFilePath) + "): " + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                success = false;
            }

            return success;
        }

        protected string MassErrorToString(double massErrorDa)
        {
            if (Math.Abs(massErrorDa) < 0.000001)
                return "0";

            return Math.Abs(massErrorDa) < 0.0001 ?
                PRISM.StringUtilities.DblToString(massErrorDa, 6, 0.0000001) :
                PRISM.StringUtilities.DblToString(massErrorDa, 5, 0.000001);
        }

        protected void OperationComplete()
        {
            ProgressComplete?.Invoke();
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">PSM tool results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <returns>True if successful, False if failure</returns>
        public bool ProcessFile(string inputFilePath, string outputDirectoryPath)
        {
            return ProcessFile(inputFilePath, outputDirectoryPath, string.Empty);
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">PSM tool results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if successful, False if failure</returns>
        public abstract bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath);

        /// <summary>
        /// Appends newSuffix to the base name of the original file, then returns the full path
        /// Note that newSuffix may contain a file extension though it does not have to
        /// If newSuffix does not contain an extension, the path returned will end in the same extension as originalFilePath
        /// </summary>
        /// <param name="originalFile"></param>
        /// <param name="newSuffix"></param>
        /// <returns>Full path to the file using the directory associated with originalFilePath</returns>
        protected string ReplaceFilenameSuffix(FileInfo originalFile, string newSuffix)
        {
            // Keep track of the original extension on originalFilePath
            var originalExtension = originalFile.Extension;

            // Make sure newSuffix is not nothing
            if (newSuffix == null)
                newSuffix = string.Empty;

            // Obtain the filename, without its extension
            var newFileBaseName = Path.GetFileNameWithoutExtension(originalFile.Name);

            // Append newSuffix to newFileBaseName
            string newFileName;
            if (Path.HasExtension(newSuffix))
            {
                newFileName = newFileBaseName + newSuffix;
            }
            else
            {
                newFileName = newFileBaseName + newSuffix + originalExtension;
            }

            if (string.IsNullOrWhiteSpace(originalFile.DirectoryName))
                return newFileName;

            var newFilePath = Path.Combine(originalFile.DirectoryName, newFileName);

            return newFilePath;
        }

        protected void ReportError(string errMsg, bool throwException = false, Exception ex = null)
        {
            SetErrorMessage(errMsg);

            if (throwException)
            {
                if (ex == null)
                {
                    throw new Exception(errMsg);
                }

                throw new Exception(errMsg, ex);
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
            var fileNotFound = false;

            // Note: If mMassCorrectionTagsFilePath is blank, the mass correction tags will be reset to the defaults and success will be True
            var success = mPeptideMods.ReadMassCorrectionTagsFile(Options.MassCorrectionTagsFilePath, ref fileNotFound);
            if (!success)
            {
                if (fileNotFound)
                {
                    SetErrorCode(PHRPErrorCode.MassCorrectionTagsFileNotFound);
                }
                else
                {
                    SetErrorCode(PHRPErrorCode.ErrorReadingMassCorrectionTagsFile);
                }
            }

            // Note: If mModificationDefinitionsFilePath is blank, the modifications will be cleared and success will be True
            success = mPeptideMods.ReadModificationDefinitionsFile(Options.ModificationDefinitionsFilePath, ref fileNotFound);
            if (!success)
            {
                if (fileNotFound)
                {
                    SetErrorCode(PHRPErrorCode.ModificationDefinitionFileNotFound);
                    ReportWarning("File not found: " + Options.ModificationDefinitionsFilePath);
                }
                else
                {
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
            }

            return success;
        }

        protected void ResetProgress()
        {
            ResetProgress(string.Empty);
        }

        protected void ResetProgress(string progressStepDescription, bool echoToConsole = false)
        {
            mProgressStepDescription = string.Copy(progressStepDescription);
            mProgressPercentComplete = 0;
            ProgressReset?.Invoke();

            if (echoToConsole)
            {
                Console.WriteLine();
                Console.WriteLine();
                Console.WriteLine(ProgressStepDescription);
            }
        }

        protected void SaveModificationSummaryFile(string modificationSummaryFilePath)
        {
            using var writer = new StreamWriter(modificationSummaryFilePath, false);

            // Write the header line
            var headerNames = new List<string>
            {
                PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol,
                PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass,
                PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Target_Residues,
                PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Type,
                PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Mass_Correction_Tag,
                PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Occurrence_Count
            };

            writer.WriteLine(CollapseList(headerNames));

            for (var index = 0; index <= mPeptideMods.ModificationCount - 1; index++)
            {
                var modInfo = mPeptideMods.GetModificationByIndex(index);
                var searchResult = modInfo;
                if (searchResult.OccurrenceCount <= 0 && searchResult.UnknownModAutoDefined)
                    continue;

                var data = new List<string>
                {
                    searchResult.ModificationSymbol.ToString(),
                    searchResult.ModificationMass.ToString(CultureInfo.InvariantCulture),
                    searchResult.TargetResidues,
                    ModificationDefinition.ModificationTypeToModificationSymbol(searchResult.ModificationType).ToString(),
                    searchResult.MassCorrectionTag,
                    searchResult.OccurrenceCount.ToString()
                };
                writer.WriteLine(CollapseList(data));
            }
        }

        protected void SaveResultsFileEntrySeqInfo(SearchResultsBaseClass searchResult, bool updateResultToSeqMapFile)
        {
            // Note: Be sure to call Me.InitializeOutputFiles before calling this function
            // updateResultToSeqMapFile should be set to True only for the first protein of each peptide in each group

            // This ID is assigned using a SortedSet containing mPeptideCleanSequence and mPeptideModDescription
            var uniqueSeqID = mUniqueSequences.GetNextUniqueSequenceID(
                searchResult.PeptideCleanSequence,
                searchResult.PeptideModDescription,
                out var existingSequenceFound);

            if (updateResultToSeqMapFile)
            {
                // Write a new entry to the ResultToSeqMap file
                var seqMapData = new List<string>
                {
                    searchResult.ResultID.ToString(),
                    uniqueSeqID.ToString()
                };
                mResultToSeqMapFile.WriteLine(CollapseList(seqMapData));

                // Only write this entry to the SeqInfo and ModDetails files if existingSequenceFound is False

                if (!existingSequenceFound)
                {
                    // Write a new entry to the SeqInfo file
                    var seqInfoData = new List<string>
                    {
                        uniqueSeqID.ToString(),
                        searchResult.SearchResultModificationCount.ToString(),
                        searchResult.PeptideModDescription,
                        PRISM.StringUtilities.DblToString(searchResult.PeptideMonoisotopicMass, 5, 0.000001)
                    };
                    mSeqInfoFile.WriteLine(CollapseList(seqInfoData));

                    if (searchResult.SearchResultModificationCount > 0)
                    {
                        var udtModNameAndResidueLoc = new ModNameAndResidueLoc[searchResult.SearchResultModificationCount];
                        var pointerArray = new int[searchResult.SearchResultModificationCount];

                        if (searchResult.SearchResultModificationCount == 1)
                        {
                            pointerArray[0] = 0;
                        }
                        else
                        {
                            // Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                            for (var index = 0; index <= searchResult.SearchResultModificationCount - 1; index++)
                            {
                                var resultModDetails = searchResult.GetSearchResultModDetailsByIndex(index);
                                udtModNameAndResidueLoc[index].ResidueLocInPeptide = resultModDetails.ResidueLocInPeptide;
                                udtModNameAndResidueLoc[index].ModName = resultModDetails.ModDefinition.MassCorrectionTag;
                                pointerArray[index] = index;
                            }

                            Array.Sort(udtModNameAndResidueLoc, pointerArray, new IModNameAndResidueLocComparer());
                        }

                        // Write out the modifications to the ModDetails file
                        // Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
                        for (var index = 0; index <= searchResult.SearchResultModificationCount - 1; index++)
                        {
                            var resultModDetails = searchResult.GetSearchResultModDetailsByIndex(pointerArray[index]);

                            var modDetailsData = new List<string>
                            {
                                uniqueSeqID.ToString(),
                                resultModDetails.ModDefinition.MassCorrectionTag,
                                resultModDetails.ResidueLocInPeptide.ToString()
                            };

                            mModDetailsFile.WriteLine(CollapseList(modDetailsData));
                        }
                    }
                }
            }

            // Write a new entry to the SeqToProteinMap file if not yet defined
            if (!CheckSeqToProteinMapDefined(uniqueSeqID, searchResult.ProteinName))
            {
                var seqToProteinData = new List<string>
                {
                    uniqueSeqID.ToString(),
                    Convert.ToInt32(searchResult.CleavageState).ToString(),
                    Convert.ToInt32(searchResult.TerminusState).ToString(),
                    searchResult.ProteinName,
                    searchResult.ProteinExpectationValue,
                    searchResult.ProteinIntensity
                };

                mSeqToProteinMapFile.WriteLine(CollapseList(seqToProteinData));
            }
        }

        protected void SetErrorCode(PHRPErrorCode newErrorCode)
        {
            SetErrorCode(newErrorCode, false);
        }

        protected void SetErrorCode(PHRPErrorCode newErrorCode, bool leaveExistingErrorCodeUnchanged)
        {
            if (leaveExistingErrorCodeUnchanged && mErrorCode != PHRPErrorCode.NoError)
            {
                // An error code is already defined; do not change it
            }
            else
            {
                mErrorCode = newErrorCode;
            }
        }

        protected void SetErrorMessage(string message, Exception ex = null)
        {
            if (message == null)
                message = string.Empty;

            mErrorMessage = message;
            if (message.Length > 0)
            {
                OnErrorEvent(message, ex);
            }
        }

        protected void ShowPeriodicWarning(int warningCount, int thresholdCountAlwaysShow, string warningMessage)
        {
            if (warningCount <= thresholdCountAlwaysShow ||
                warningCount < 1000 && warningCount % 100 == 0 ||
                warningCount < 10000 && warningCount % 1000 == 0 ||
                warningCount < 100000 && warningCount % 10000 == 0 ||
                warningCount < 1000000 && warningCount % 100000 == 0)
            {
                ReportWarning(warningMessage);
            }
        }

        /// <summary>
        /// If resultID is 0 or 1, returns valueText
        /// Otherwise, if valueText is 0.0, returns 0
        /// Otherwise, returns valueText
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="valueText"></param>
        protected string TrimZeroIfNotFirstID(int resultID, string valueText)
        {
            return resultID > 1 ? TrimZero(valueText) : valueText;
        }

        /// <summary>
        /// If valueText is 0.0, returns 0
        /// Otherwise, returns valueText
        /// </summary>
        /// <param name="valueText"></param>
        private string TrimZero(string valueText)
        {
            return valueText.Equals("0.0") ? "0" : valueText;
        }

        /// <summary>
        /// Return the text up to (but not including) the first space in proteinNameAndDescription
        /// </summary>
        /// <param name="proteinNameAndDescription"></param>
        protected virtual string TruncateProteinName(string proteinNameAndDescription)
        {
            var index = proteinNameAndDescription.IndexOf(' ');
            if (index > 0)
            {
                return proteinNameAndDescription.Substring(0, index);
            }

            return proteinNameAndDescription;
        }

        protected void UpdatePepToProteinMapPeptide(List<PepToProteinMapping> pepToProteinMapping, int index, string peptide)
        {
            var udtItem = pepToProteinMapping[index];
            udtItem.Peptide = peptide;
            pepToProteinMapping[index] = udtItem;
        }

        protected void UpdateProgress(string progressStepDescription)
        {
            UpdateProgress(progressStepDescription, mProgressPercentComplete);
        }

        protected void UpdateProgress(float percentComplete)
        {
            UpdateProgress(ProgressStepDescription, percentComplete);
        }

        protected void UpdateProgress(string progressStepDescription, float percentComplete)
        {
            mProgressStepDescription = string.Copy(progressStepDescription);
            if (percentComplete < 0)
            {
                percentComplete = 0;
            }
            else if (percentComplete > 100)
            {
                percentComplete = 100;
            }
            mProgressPercentComplete = percentComplete;

            OnProgressUpdate(ProgressStepDescription, ProgressPercentComplete);
        }

        /// <summary>
        /// Update progress while creating the synopsis file
        /// </summary>
        /// <param name="reader"></param>
        protected void UpdateSynopsisFileCreationProgress(StreamReader reader)
        {
            var percentComplete = reader.BaseStream.Position / (float)reader.BaseStream.Length * 100;

            if (!Options.CreateProteinModsFile)
            {
                UpdateProgress(percentComplete);
                return;
            }

            var overallProgress = ProcessFilesOrDirectoriesBase.ComputeIncrementalProgress(
                0, PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE, percentComplete);

            UpdateProgress(overallProgress);
        }

        /// <summary>
        /// Validate that the specified file exists and has at least one tab-delimited row with a numeric value in the first column
        /// </summary>
        /// <param name="filePath">Path to the file</param>
        /// <param name="fileDescription">File description, e.g. Synopsis</param>
        /// <param name="errorMessage"></param>
        /// <returns>True if the file has data; otherwise false</returns>
        public static bool ValidateFileHasData(string filePath, string fileDescription, out string errorMessage)
        {
            const int numericDataColIndex = 0;
            return ValidateFileHasData(filePath, fileDescription, out errorMessage, numericDataColIndex);
        }

        /// <summary>
        /// Validate that the specified file exists and has at least one tab-delimited row with a numeric value
        /// </summary>
        /// <param name="filePath">Path to the file</param>
        /// <param name="fileDescription">File description, e.g. Synopsis</param>
        /// <param name="errorMessage"></param>
        /// <param name="numericDataColIndex">Index of the numeric data column; use -1 to simply look for any text in the file</param>
        /// <returns>True if the file has data; otherwise false</returns>
        public static bool ValidateFileHasData(string filePath, string fileDescription, out string errorMessage, int numericDataColIndex)
        {
            var dataFound = false;

            errorMessage = string.Empty;

            try
            {
                var fileInfo = new FileInfo(filePath);

                if (!fileInfo.Exists)
                {
                    errorMessage = fileDescription + " file not found: " + fileInfo.Name;
                    return false;
                }

                if (fileInfo.Length == 0)
                {
                    errorMessage = fileDescription + " file is empty (zero-bytes)";
                    return false;
                }

                // Open the file and confirm it has data rows
                using (var reader = new StreamReader(new FileStream(fileInfo.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream && !dataFound)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        if (numericDataColIndex < 0)
                        {
                            dataFound = true;
                        }
                        else
                        {
                            // Split on the tab character and check if the first column is numeric
                            var splitLine = lineIn.Split('\t');

                            if (splitLine.Length <= numericDataColIndex)
                                continue;

                            if (double.TryParse(splitLine[numericDataColIndex], out _))
                            {
                                dataFound = true;
                            }
                        }
                    }
                }

                if (!dataFound)
                {
                    errorMessage = fileDescription + " is empty (no data)";
                }
            }
            catch (Exception)
            {
                errorMessage = "Exception validating " + fileDescription + " file";
                return false;
            }

            return dataFound;
        }

        /// <summary>
        /// Compare the two mass values; warn the user if more than 0.1 Da apart (slightly larger threshold if over 5000 Da)
        /// </summary>
        /// <param name="toolName">Tool name</param>
        /// <param name="peptide">Peptide sequence</param>
        /// <param name="peptideMonoMassFromPHRP">Peptide monoisotopic mass, as computed by PHRP</param>
        /// <param name="peptideMonoMassFromTool">Peptide monoisotopic mass, as computed by the search engine</param>
        /// <param name="deltaMassWarningCount">Keeps track of the number of times the monoisotopic mass values differ more than the threshold</param>
        protected void ValidateMatchingMonoisotopicMass(
            string toolName,
            string peptide,
            double peptideMonoMassFromPHRP,
            double peptideMonoMassFromTool,
            ref int deltaMassWarningCount)
        {
            var massDiffThreshold = peptideMonoMassFromTool / 5000 / 10;
            if (massDiffThreshold < 0.1)
                massDiffThreshold = 0.1;

            if (Math.Abs(peptideMonoMassFromPHRP - peptideMonoMassFromTool) <= massDiffThreshold)
                return;

            // Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da
            // (or by a slightly larger value if over 5000 Da)

            // This is unexpected

            string first30Residues;
            if (peptide.Length < 27)
            {
                first30Residues = peptide;
            }
            else
            {
                first30Residues = peptide.Substring(0, 27) + "...";
            }

            deltaMassWarningCount++;
            ShowPeriodicWarning(deltaMassWarningCount,
                                10,
                                string.Format(
                                    "The monoisotopic mass computed by PHRP is more than {0:F2} Da away from " +
                                    "the mass computed by {1}: {2:F4} vs. {3:F4}; peptide {4}",
                                    massDiffThreshold, toolName, peptideMonoMassFromPHRP, peptideMonoMassFromTool, first30Residues));
        }

        private bool ValidatePeptideToProteinMapResults(string peptideToProteinMapFilePath, bool ignorePeptideToProteinMapperErrors)
        {
            bool success;

            var peptideCount = 0;
            var peptideCountNoMatch = 0;
            var linesRead = 0;
            var chSplitChars = new[] { '\t' };

            try
            {
                // Validate that none of the results in peptideToProteinMapFilePath has protein name PROTEIN_NAME_NO_MATCH ( __NoMatch__ )

                var lastPeptide = string.Empty;

                using (var reader = new StreamReader(new FileStream(peptideToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        linesRead++;

                        if (linesRead <= 1 || string.IsNullOrEmpty(lineIn))
                            continue;

                        var splitLine = lineIn.Split(chSplitChars, 2);
                        if (splitLine.Length == 0)
                            continue;

                        if (splitLine[0] != lastPeptide)
                        {
                            peptideCount++;
                            lastPeptide = string.Copy(splitLine[0]);
                        }

                        if (lineIn.Contains(PROTEIN_NAME_NO_MATCH))
                        {
                            peptideCountNoMatch++;
                        }
                    }
                }

                if (peptideCount == 0)
                {
                    SetErrorMessage("Peptide to protein mapping file is empty: " + peptideToProteinMapFilePath);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    success = false;
                }
                else if (peptideCountNoMatch == 0)
                {
                    success = true;
                }
                else
                {
                    // Value between 0 and 100
                    var errorPercent = peptideCountNoMatch / (double)peptideCount * 100.0;

                    var message = string.Format("{0:0.00}% of the entries ({1:N0} / {2:N0}) in the peptide to protein map file ({3}) " +
                                                "did not match to a protein in the FASTA file ({4})",
                                                errorPercent, peptideCountNoMatch, peptideCount,
                                                Path.GetFileName(peptideToProteinMapFilePath),
                                                Path.GetFileName(Options.FastaFilePath));

                    if (ignorePeptideToProteinMapperErrors || errorPercent < 0.1)
                    {
                        ReportWarning(message);
                        success = true;
                    }
                    else
                    {
                        SetErrorMessage(message);
                        success = false;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ValidatePeptideToProteinMapResults:" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        protected void ValidatePHRPReaderSupportFiles(string phrpDataFilePath, string outputDirectoryPath)
        {
            try
            {
                if (string.IsNullOrWhiteSpace(outputDirectoryPath))
                    return;

                var phrpDataFile = new FileInfo(phrpDataFilePath);
                var outputDirectory = new DirectoryInfo(outputDirectoryPath);

                if (string.Equals(phrpDataFile.DirectoryName, outputDirectory.FullName, StringComparison.OrdinalIgnoreCase))
                    return;

                var msgfFileName = Path.GetFileName(ReplaceFilenameSuffix(phrpDataFile, FILENAME_SUFFIX_MSGF));

                if (string.IsNullOrWhiteSpace(msgfFileName))
                    return;

                string sourcePath;

                if (string.IsNullOrWhiteSpace(phrpDataFile.DirectoryName))
                    sourcePath = msgfFileName;
                else
                    sourcePath = Path.Combine(phrpDataFile.DirectoryName, msgfFileName);

                var targetPath = Path.Combine(outputDirectory.FullName, msgfFileName);

                if (File.Exists(sourcePath) && !File.Exists(targetPath))
                {
                    File.Copy(sourcePath, targetPath);
                }
            }
            catch (Exception ex)
            {
                ReportWarning("Error in ValidatePHRPReaderSupportFiles: " + ex.Message);
            }
        }

        protected bool ValidateProteinFastaFile(string fastaFilePath)
        {
            var success = ValidateProteinFastaFile(fastaFilePath, out var warningMessage);

            if (!success)
            {
                ReportWarning(warningMessage);
            }

            return success;
        }

        public static bool ValidateProteinFastaFile(string fastaFilePath, out string warningMessage)
        {
            // This RegEx looks for standard amino acids, skipping A, T, C, and G
            var reDefiniteAminoAcid = new Regex("[DEFHIKLMNPQRSVWY]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            // This RegEx looks for A, T, C, and G
            var rePotentialNucleicAcid = new Regex("[ATCG]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            // This matches any letter
            var reLetter = new Regex("[A-Z]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            var validProteinCount = 0;
            var invalidProteinCount = 0;

            try
            {
                warningMessage = string.Empty;

                if (string.IsNullOrEmpty(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "fastaFilePath is not defined in ValidateProteinFastaFile";
                    return false;
                }

                if (!File.Exists(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "Fasta file not found: " + fastaFilePath;
                    return false;
                }

                var fastaFile = new ProteinFileReader.FastaFileReader();
                if (!fastaFile.OpenFile(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "Error opening the fasta file: " + fastaFilePath;
                    return false;
                }

                // Read the first 500 proteins and confirm that each contains amino acid residues
                while (fastaFile.ReadNextProteinEntry())
                {
                    var definiteAminoAcidCount = reDefiniteAminoAcid.Matches(fastaFile.ProteinSequence).Count;
                    var potentialNucleicAcidCount = rePotentialNucleicAcid.Matches(fastaFile.ProteinSequence).Count;
                    var letterCount = reLetter.Matches(fastaFile.ProteinSequence).Count;

                    if (definiteAminoAcidCount > 0.1 * letterCount)
                    {
                        validProteinCount++;
                    }
                    else if (potentialNucleicAcidCount > 0.95 * letterCount)
                    {
                        invalidProteinCount++;
                    }

                    if (validProteinCount + invalidProteinCount >= 500)
                    {
                        break;
                    }
                }

                if (validProteinCount < invalidProteinCount)
                {
                    Console.WriteLine();
                    warningMessage = "Fasta file contains Nucleic Acids, not Amino Acids: " + Path.GetFileName(fastaFilePath);
                    return false;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine();
                warningMessage = "Exception in ValidateProteinFastaFile: " + ex.Message;
                return false;
            }

            return true;
        }

        private void WriteModDetailsEntry(
            ReaderFactory reader,
            TextWriter writer,
            IReadOnlyList<PepToProteinMapping> pepToProteinMapping,
            int pepToProteinMapIndex,
            ref int psmCountSkippedSinceReversedOrScrambledProtein)
        {
            foreach (var mod in reader.CurrentPSM.ModifiedResidues)
            {
                var residueLocInProtein = pepToProteinMapping[pepToProteinMapIndex].ResidueStart + mod.ResidueLocInPeptide - 1;
                string residue;

                if (IsLetterAtoZ(mod.Residue))
                {
                    residue = mod.Residue.ToString();
                }
                else
                {
                    var cleanSequence = reader.CurrentPSM.PeptideCleanSequence;

                    if (mod.ResidueLocInPeptide < 1)
                    {
                        // This shouldn't be the case, but we'll check for it anyway
                        residue = cleanSequence.Substring(0, 1);
                    }
                    else if (mod.ResidueLocInPeptide > cleanSequence.Length)
                    {
                        // This shouldn't be the case, but we'll check for it anyway
                        residue = cleanSequence.Substring(cleanSequence.Length - 1, 1);
                    }
                    else
                    {
                        residue = cleanSequence.Substring(mod.ResidueLocInPeptide - 1, 1);
                    }
                }

                if (pepToProteinMapping[pepToProteinMapIndex].Protein == PROTEIN_NAME_NO_MATCH && IsReversedProtein(reader.CurrentPSM.ProteinFirst))
                {
                    // Skip this result
                    psmCountSkippedSinceReversedOrScrambledProtein++;
                }
                else
                {
                    writer.WriteLine(reader.CurrentPSM.ResultID + "\t" +
                                     reader.CurrentPSM.Peptide + "\t" +
                                     reader.CurrentPSM.SeqID + "\t" +
                                     pepToProteinMapping[pepToProteinMapIndex].Protein + "\t" +
                                     residue + "\t" +
                                     residueLocInProtein + "\t" +
                                     mod.ModDefinition.MassCorrectionTag + "\t" +
                                     mod.ResidueLocInPeptide + "\t" +
                                     reader.CurrentPSM.MSGFSpecEValue);
                }
            }
        }

        /// <summary>
        /// Override this method to display the name of each class
        /// </summary>
        public abstract override string ToString();

        private void PeptideToProteinMapper_ProgressChanged(string taskDescription, float percentComplete)
        {
            if (percentComplete >= mNextPeptideToProteinMapperLevel)
            {
                mNextPeptideToProteinMapperLevel += 25;

                var overallProgress = ProcessFilesOrDirectoriesBase.ComputeIncrementalProgress(
                    PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE, PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE, percentComplete);

                UpdateProgress(overallProgress);

                if (!HasEventListenerProgressUpdate)
                {
                    Console.WriteLine(" PeptideToProteinMapper is {0:0}% complete", percentComplete);
                }
            }
        }

        protected class ISearchOptionModificationInfoComparer : IComparer<SearchOptionModificationInfo>
        {
            public int Compare(SearchOptionModificationInfo x, SearchOptionModificationInfo y)
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

        internal class IModNameAndResidueLocComparer : IComparer<ModNameAndResidueLoc>
        {
            public int Compare(ModNameAndResidueLoc x, ModNameAndResidueLoc y)
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

                return string.CompareOrdinal(x.ModName, y.ModName);
            }
        }

        protected class PepToProteinMappingComparer : IComparer<PepToProteinMapping>
        {
            public int Compare(PepToProteinMapping x, PepToProteinMapping y)
            {
                var result = string.CompareOrdinal(x.Peptide, y.Peptide);
                if (result == 0)
                {
                    result = string.CompareOrdinal(x.Protein, y.Protein);
                }
                return result;
            }
        }

        private class PepToProteinMappingPeptideSearchComparer : IComparer<PepToProteinMapping>
        {
            public int Compare(PepToProteinMapping x, PepToProteinMapping y)
            {
                return string.CompareOrdinal(x.Peptide, y.Peptide);
            }
        }
    }
}
