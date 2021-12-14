// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.DirectoryServices.ActiveDirectory;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using PeptideHitResultsProcessor.Data;
using PeptideHitResultsProcessor.SearchToolResults;
using PeptideToProteinMapEngine;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM;
using PRISM.FileProcessor;
using PRISMDatabaseUtils;
using ProteinCoverageSummarizer;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class can be used as a base class for peptide hit results processor classes
    /// </summary>
    public abstract class PHRPBaseClass : EventNotifier
    {
        // Ignore Spelling: A-Za-z, Da, Daltons, Fscore, MaxQuant, MSFragger, mts, pre, prot, xxx

        /// <summary>
        /// Program date
        /// </summary>
        public const string PROGRAM_DATE = "December 13, 2021";

        /// <summary>
        /// Constructor
        /// </summary>
        protected PHRPBaseClass(PHRPOptions options)
        {
            FileDate = PROGRAM_DATE;

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

        /// <summary>
        /// ResultID column
        /// </summary>
        protected const string COLUMN_NAME_RESULTID = "ResultID";

        /// <summary>
        /// Peptide column
        /// </summary>
        protected const string COLUMN_NAME_PEPTIDE = "Peptide";

        private const string COLUMN_NAME_RESIDUE = "Residue";

        private const string COLUMN_NAME_PROTEIN_RESIDUE_NUMBER = "Protein_Residue_Num";

        private const string COLUMN_NAME_RESIDUE_MOD_NAME = "Mod_Name";

        private const string COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER = "Peptide_Residue_Num";

        private const string COLUMN_NAME_MSGF_SPECPROB = "MSGF_SpecProb";

        /// <summary>
        /// X!Tandem results file suffix
        /// </summary>
        public const string XTANDEM_RESULTS_FILE_SUFFIX = "_xt.xml";

        /// <summary>
        /// Synopsis file suffix
        /// </summary>
        public const string SYNOPSIS_FILE_SUFFIX = "_syn.txt";

        /// <summary>
        /// First-hits file suffix
        /// </summary>
        public const string FIRST_HITS_FILE_SUFFIX = "_fht.txt";

        /// <summary>
        /// Inspect results file suffix
        /// </summary>
        public const string INSPECT_RESULTS_FILE_SUFFIX = "_inspect.txt";

        /// <summary>
        /// Inspect TotalPRMScore-filtered first hits file suffix
        /// </summary>
        public const string INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX = "_fht.txt";

        /// <summary>
        /// Inspect FScore filtered first hits file suffix
        /// </summary>
        public const string INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX = "_Fscore_fht.txt";

        /// <summary>
        /// MSGFDB results file suffix
        /// </summary>
        public const string MSGFDB_RESULTS_FILE_SUFFIX = "_msgfdb.txt";

        /// <summary>
        /// MSAlign results file suffix
        /// </summary>
        public const string MSALIGN_RESULTS_FILE_SUFFIX = "_MSAlign_ResultTable.txt";

        /// <summary>
        /// MODa results file suffix
        /// </summary>
        public const string MODa_RESULTS_FILE_SUFFIX = "_moda.id.txt";

        /// <summary>
        /// MODPlus results file suffix
        /// </summary>
        public const string MODPlus_RESULTS_FILE_SUFFIX = "_modp.id.txt";

        /// <summary>
        /// MSPathFinder results file suffix
        /// </summary>
        public const string MSPathFinder_RESULTS_FILE_SUFFIX = "_IcTda.tsv";

        /// <summary>
        /// TopPIC results file suffix
        /// </summary>
        public const string TopPIC_RESULTS_FILE_SUFFIX = "_TopPIC_PrSMs.txt";

        /// <summary>
        /// Result to sequence map file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_RESULT_TO_SEQ_MAP = "_ResultToSeqMap.txt";

        /// <summary>
        /// Sequence to protein map file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SEQ_TO_PROTEIN_MAP = "_SeqToProteinMap.txt";

        /// <summary>
        /// Sequence info file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SEQ_INFO = "_SeqInfo.txt";

        /// <summary>
        /// Modification details file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_MOD_DETAILS = "_ModDetails.txt";

        /// <summary>
        /// Modification summary file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_MOD_SUMMARY = "_ModSummary.txt";

        /// <summary>
        /// Peptide to protein map file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING = "_PepToProtMap";

        /// <summary>
        /// Protein modifications file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_PROTEIN_MODS = "_ProteinMods.txt";

        /// <summary>
        /// MSGF score file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_MSGF = "_MSGF.txt";

        /// <summary>
        /// Progress percent value at the start of peptide to protein mapping
        /// </summary>
        protected const float PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE = 90;

        /// <summary>
        /// Progress percent value at the start of creating the protein modifications file
        /// </summary>
        private const float PROGRESS_PERCENT_CREATING_PROTEIN_MODS_FILE = 95;

        /// <summary>
        /// Protein name to use when a peptide does not match any of the loaded proteins
        /// </summary>
        private const string PROTEIN_NAME_NO_MATCH = "__NoMatch__";

        /// <summary>
        /// ProteinMods file skipped warning message
        /// </summary>
        public const string WARNING_MESSAGE_SKIPPING_PROTEIN_MODS_FILE_CREATION = "Skipping creation of the ProteinMods file";

        /// <summary>
        /// Search options for finding modifications
        /// </summary>
        protected struct SearchOptionModificationInfo
        {
            /// <summary>
            /// Sort order
            /// </summary>
            public int SortOrder;

            /// <summary>
            /// Modification mass to find
            /// </summary>
            public double ModificationMass;

            /// <summary>
            /// Target residues to match
            /// </summary>
            public string TargetResidues;

            /// <summary>
            /// Modification type to match
            /// </summary>
            public ModificationDefinition.ResidueModificationType ModificationType;

            /// <summary>
            /// Duplicate this modification via a deep copy
            /// </summary>
            public SearchOptionModificationInfo Clone()
            {
                return new SearchOptionModificationInfo
                {
                    SortOrder = SortOrder,
                    ModificationMass = ModificationMass,
                    TargetResidues = TargetResidues,
                    ModificationType = ModificationType
                };
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

        /// <summary>
        /// Class for tracking the location of a peptide in a protein
        /// </summary>
        protected struct PepToProteinMapping
        {
            /// <summary>
            /// peptide sequence
            /// </summary>
            public string Peptide;

            /// <summary>
            /// Protein name
            /// </summary>
            public string Protein;

            /// <summary>
            /// Location in the protein sequence where the peptide starts
            /// </summary>
            public int ResidueStart;

            /// <summary>
            /// Location in the protein sequence where the peptide ends
            /// </summary>
            public int ResidueEnd;

            /// <summary>
            /// Show the peptide sequence and protein name
            /// </summary>
            public override string ToString()
            {
                return Peptide + ", " + Protein;
            }
        }

        /// <summary>
        /// Error code
        /// </summary>
        protected PHRPErrorCode mErrorCode;

        /// <summary>
        /// Error message
        /// </summary>
        protected string mErrorMessage;

        /// <summary>
        /// Peptide sequence mass calculator
        /// </summary>
        protected readonly PeptideMassCalculator mPeptideSeqMassCalculator;

        /// <summary>
        /// Tracks peptide modifications
        /// </summary>
        protected readonly PeptideModificationContainer mPeptideMods;

        /// <summary>
        /// Tracks unique sequences
        /// </summary>
        private readonly UniqueSequencesContainer mUniqueSequences;

        /// <summary>
        /// Tracks unique sequence IDs and corresponding protein names
        /// </summary>
        /// <remarks>
        /// Data is stored as UniqueSequenceID_ProteinName
        /// </remarks>
        private readonly SortedSet<string> mSeqToProteinMap;

        private StreamWriter mResultToSeqMapFile;
        private StreamWriter mSeqInfoFile;
        private StreamWriter mModDetailsFile;
        private StreamWriter mSeqToProteinMapFile;

        private int mNextPeptideToProteinMapperLevel;

        /// <summary>
        /// Tracks the protein names in the order that they are listed in the FASTA file
        /// Keys are protein names, values are a sequentially assigned integer (starting with 1)
        /// </summary>
        protected readonly Dictionary<string, int> mProteinNameOrder;

        private readonly Regex mReplaceSymbols;

        /// <summary>
        /// Progress reset event
        /// </summary>
        public event ProgressResetEventHandler ProgressReset;

        /// <summary>
        /// Progress reset event handler
        /// </summary>
        public delegate void ProgressResetEventHandler();

        /// <summary>
        /// Progress complete event
        /// </summary>
        public event ProgressCompleteEventHandler ProgressComplete;

        /// <summary>
        /// Progress complete event handler
        /// </summary>
        public delegate void ProgressCompleteEventHandler();

        /// <summary>
        /// Progress step description
        /// </summary>
        protected string mProgressStepDescription = string.Empty;

        /// <summary>
        /// Ranges from 0 to 100, but can contain decimal percentage values
        /// </summary>
        protected float mProgressPercentComplete;

        /// <summary>
        /// The calling procedure can set this to true to request that processing should be aborted
        /// </summary>
        public bool AbortProcessing { get; set; }

        /// <summary>
        /// Error code enum
        /// </summary>
        public PHRPErrorCode ErrorCode => mErrorCode;

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage => GetErrorMessage();

        /// <summary>
        /// Version of the executing assembly
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public string FileVersion => GetVersionForExecutingAssembly();

        /// <summary>
        /// File date
        /// </summary>
        public string FileDate { get; protected set; }

        /// <summary>
        /// Processing options
        /// </summary>
        public PHRPOptions Options { get; }

        /// <summary>
        /// Description of the current processing step
        /// </summary>
        public string ProgressStepDescription => mProgressStepDescription;

        /// <summary>
        /// Overall processing percent complete (value between 0 and 100)
        /// </summary>
        public float ProgressPercentComplete => Convert.ToSingle(Math.Round(mProgressPercentComplete, 2));

        /// <summary>
        /// The calling procedure can call this method to request that processing should be aborted
        /// </summary>
        /// <remarks>
        /// Sets property AbortProcessing to true
        /// </remarks>
        // ReSharper disable once UnusedMember.Global
        public void AbortProcessingNow()
        {
            AbortProcessing = true;
        }

        /// <summary>
        /// Construct the peptide hit results file path for the given results file format
        /// </summary>
        /// <param name="peptideHitResultFileFormat"></param>
        /// <param name="sourceDirectoryPath"></param>
        /// <param name="baseName"></param>
        // ReSharper disable once UnusedMember.Global
        public static string AutoDefinePeptideHitResultsFilePath(
            ResultsFileFormat peptideHitResultFileFormat,
            string sourceDirectoryPath,
            string baseName)
        {
            if (string.IsNullOrEmpty(baseName))
                return AutoDefinePeptideHitResultsFilePath(sourceDirectoryPath);

            return peptideHitResultFileFormat switch
            {
                ResultsFileFormat.SequestFirstHitsFile => Path.Combine(sourceDirectoryPath, baseName + FIRST_HITS_FILE_SUFFIX),
                ResultsFileFormat.SequestSynopsisFile => Path.Combine(sourceDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX),
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

        /// <summary>
        /// Auto-define the peptide hit results file path by looking for a file with a known suffix in the given directory
        /// </summary>
        /// <remarks>
        /// Looks for a file ending in _syn.txt, _fht.txt, _xt.xml, or _inspect.txt in directory sourceDirectoryPath
        /// </remarks>
        /// <param name="sourceDirectoryPath"></param>
        /// <returns>The first matching file found, or an empty string</returns>
        public static string AutoDefinePeptideHitResultsFilePath(string sourceDirectoryPath)
        {
            try
            {
                for (var index = 0; index <= 3; index++)
                {
                    var matchSpec = index switch
                    {
                        0 => "*" + SYNOPSIS_FILE_SUFFIX,
                        1 => "*" + FIRST_HITS_FILE_SUFFIX,
                        2 => "*" + XTANDEM_RESULTS_FILE_SUFFIX,
                        3 => "*" + INSPECT_RESULTS_FILE_SUFFIX,
                        _ => throw new ArgumentOutOfRangeException()
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

        /// <summary>
        /// Cache the protein names defined by Options.FastaFilepath
        /// </summary>
        protected bool CacheProteinNamesFromFasta()
        {
            if (string.IsNullOrWhiteSpace(Options.FastaFilePath))
            {
                // Nothing to do
                return true;
            }

            mProteinNameOrder.Clear();
            var proteinNameMatcher = new Regex("^>([^ ]+)", RegexOptions.Compiled);

            OnStatusEvent("Caching protein names from the FASTA file");

            try
            {
                var proteinNumber = 0;

                using var reader = new StreamReader(new FileStream(Options.FastaFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(lineIn))
                        continue;

                    var match = proteinNameMatcher.Match(lineIn);

                    if (!match.Success)
                        continue;

                    var proteinName = match.Groups[1].Value;

                    if (mProteinNameOrder.ContainsKey(proteinName))
                        continue;

                    proteinNumber++;

                    mProteinNameOrder.Add(proteinName, proteinNumber);
                }

                OnStatusEvent("Cached {0:N0} proteins", mProteinNameOrder.Count);

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error caching protein names from FASTA file " + Path.GetFileName(Options.FastaFilePath), false, ex);
                return false;
            }
        }

        /// <summary>
        /// Check whether the sequence to protein map is defined in mSeqToProteinMap, adding it if missing
        /// </summary>
        /// <param name="uniqueSeqID">Unique sequence ID</param>
        /// <param name="proteinName">Protein name</param>
        /// <returns>
        /// True if the sequence to protein map was already defined;
        /// False if the mapping was not defined (will also update mSeqToProteinMap)
        /// </returns>
        protected bool CheckSeqToProteinMapDefined(int uniqueSeqID, string proteinName)
        {
            bool existingMapFound;

            try
            {
                proteinName ??= string.Empty;

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
                        OnWarningEvent("Absolute path: " + inputFile.DirectoryName);
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
                SetErrorMessage("Error cleaning up the file paths", ex);
                SetErrorCode(PHRPErrorCode.FilePathError);
                return false;
            }
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

            if (fileName.EndsWith(MSFraggerResultsProcessor.PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase))
            {
                return ResultsFileFormat.MSFraggerTSVFile;
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

        /// <summary>
        /// Close the result to sequence map file and the sequence info file
        /// </summary>
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

        /// <summary>
        /// Compute a pseudo peptide location in a protein
        /// </summary>
        /// <remarks>
        /// <para>
        /// This method is used when a peptide is associated with a protein for which we do not know the protein's sequence
        /// </para>
        /// <para>
        /// Peptides at the N-terminus of a protein will have
        /// searchResult.PeptideLocInProteinStart = 1 and PeptideLocInProteinEnd updated accordingly
        /// </para>
        /// <para>
        /// Peptides at the C-terminus of a protein will have
        /// searchResult.PeptideLocInProteinEnd = 10000 and PeptideLocInProteinStart updated accordingly
        /// </para>
        /// <para>
        /// Peptides not at the N- or C- terminus of a protein will have
        /// searchResult.PeptideLocInProteinStart = 2 and PeptideLocInProteinEnd updated accordingly
        /// </para>
        /// </remarks>
        /// <param name="searchResult"></param>
        protected void ComputePseudoPeptideLocInProtein(SearchResultsBaseClass searchResult)
        {
            // Set these to 1 and 10000 since MSGFDB, SEQUEST, and InSpecT results files do not contain protein sequence information
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

        /// <summary>
        /// Compute the mass of the residues in the sequence, then add the modification mass (if non-zero) and return the result
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="totalModMass"></param>
        /// <returns>Mass, in Daltons</returns>
        protected double ComputePeptideMassForCleanSequence(string cleanSequence, double totalModMass)
        {
            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);

            if (Math.Abs(totalModMass) > double.Epsilon)
            {
                return mass + totalModMass;
            }

            return mass;
        }

        /// <summary>
        /// Compute FDR values, then assign QValues
        /// </summary>
        /// <remarks>Assumes the data is sorted by highest confidence to lowest confidence</remarks>
        /// <param name="searchResults"></param>
        protected void ComputeQValues(List<ToolResultsBaseClass> searchResults)
        {
            var forwardPeptideCount = 0;
            var reversePeptideCount = 0;

            foreach (var searchResult in searchResults)
            {
                if (searchResult.Reverse)
                {
                    reversePeptideCount++;
                }
                else
                {
                    forwardPeptideCount++;
                }

                double fdr = 1;

                if (forwardPeptideCount > 0)
                {
                    fdr = reversePeptideCount / Convert.ToDouble(forwardPeptideCount);
                }

                searchResult.FDR = fdr;
            }

            // Now compute Q-Values
            // We step through the list, from the worst scoring result to the best result
            // The first Q-Value is the FDR of the final entry
            // The next Q-Value is the minimum of (QValue, CurrentFDR)

            var qValue = searchResults.Last().FDR;
            if (qValue > 1)
                qValue = 1;

            for (var index = searchResults.Count - 1; index >= 0; index += -1)
            {
                qValue = Math.Min(qValue, searchResults[index].FDR);
                searchResults[index].QValue = qValue;
            }
        }

        /// <summary>
        /// Construct the peptide to protein map file path
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="mts">If true, the map file will end with MTS.txt; otherwise, just .txt</param>
        /// <returns>_PepToProtMap file that corresponds to the input file</returns>
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

            if (!Options.UseExistingMTSPepToProteinMapFile)
            {
                return alternatePepToProteinMapFile.FullName;
            }

            // ReSharper disable once ConvertIfStatementToSwitchStatement
            // ReSharper disable once ConvertIfStatementToSwitchExpression

            if (pepToProteinMapFile.Exists && !alternatePepToProteinMapFile.Exists)
            {
                return pepToProteinMapFile.FullName;
            }

            if (!pepToProteinMapFile.Exists && alternatePepToProteinMapFile.Exists)
            {
                return alternatePepToProteinMapFile.FullName;
            }

            if (pepToProteinMapFile.Exists && alternatePepToProteinMapFile.Exists &&
                pepToProteinMapFile.LastWriteTime > alternatePepToProteinMapFile.LastWriteTime)
            {
                return pepToProteinMapFile.FullName;
            }

            return alternatePepToProteinMapFile.FullName;
        }

        /// <summary>
        /// Construct the peptide to protein map file path
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="mts">If true, the map file will end with MTS.txt; otherwise, just .txt</param>
        /// <param name="suffixesToFind">Filename suffixes to find</param>
        /// <param name="charsToRemove">Number of characters to remove from the base name</param>
        /// <returns>_PepToProtMap file that corresponds to the input file</returns>
        protected string ConstructPepToProteinMapFilePath(
            string inputFilePath,
            string outputDirectoryPath,
            bool mts,
            List<string> suffixesToFind,
            int charsToRemove)
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
        /// <param name="maximumAllowableMatchErrorPercentThreshold">
        /// Maximum percentage of peptides in the peptide to protein map file that are allowed to have not matched a protein in the FASTA file (value between 0 and 100)
        /// This is typically 0.1, but for MS-GF+, MaxQuant, and other tools we set this to 50
        /// </param>
        /// <param name="matchErrorPercentWarningThreshold">
        /// When at least one peptide did not have a matched protein in the FASTA file, this threshold defines at what percent level a warning should be shown (value between 0 and 100)
        /// </param>
        protected bool CreatePepToProteinMapFile(
            List<string> sourcePHRPDataFiles,
            string mtsPepToProteinMapFilePath,
            double maximumAllowableMatchErrorPercentThreshold = 0.1,
            double matchErrorPercentWarningThreshold = 0)
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
                    SetErrorMessage("Cannot create the PepToProtein map file because the FASTA File Path is not defined");
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                if (!File.Exists(Options.FastaFilePath))
                {
                    SetErrorMessage("Cannot create the PepToProtein map file because the FASTA File was not found: " + Options.FastaFilePath);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                // Verify that the FASTA file is not a DNA-sequence based FASTA file
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

                using var writer = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));
                var headerWritten = false;

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

                        success = ValidatePeptideToProteinMapResults(
                            resultsFilePath,
                            Options.IgnorePeptideToProteinMapperErrors,
                            maximumAllowableMatchErrorPercentThreshold,
                            matchErrorPercentWarningThreshold);
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
                            OnWarningEvent("Ignoring protein mapping error since 'IgnorePeptideToProteinMapperErrors' = True");

                            if (File.Exists(resultsFilePath))
                            {
                                success = ValidatePeptideToProteinMapResults(
                                    resultsFilePath,
                                    Options.IgnorePeptideToProteinMapperErrors,
                                    maximumAllowableMatchErrorPercentThreshold,
                                    matchErrorPercentWarningThreshold);
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

                    // This tracks the number of non-empty lines read
                    var linesRead = 0;

                    using (var reader = new StreamReader(new FileStream(resultsFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!reader.EndOfStream)
                        {
                            var lineIn = reader.ReadLine();

                            if (string.IsNullOrWhiteSpace(lineIn))
                                continue;

                            linesRead++;

                            var splitLine = lineIn.Split(new[] { '\t' }, 3);
                            if (splitLine.Length < 2)
                                continue;

                            if (linesRead == 1)
                            {
                                if (splitLine[0].Equals("Peptide") && splitLine[1].Equals("Protein"))
                                {
                                    // Header line; write to disk if no results have been written yet
                                    if (!headerWritten)
                                    {
                                        writer.WriteLine(lineIn);
                                        headerWritten = true;
                                    }

                                    continue;
                                }

                                if (!headerWritten)
                                {
                                    // Write the default header line

                                    var headerLine = string.Format(
                                        "Peptide\tProtein{0}",
                                        splitLine.Length == 2 ? string.Empty : "\tResidue_Start\tResidue_End");

                                    OnWarningEvent(
                                        "Input file {0} does not have the expected header line; using the default:\n{1}",
                                        resultsFilePath, headerLine);

                                    writer.WriteLine(headerLine);
                                    headerWritten = true;
                                }
                            }

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

                peptideToProteinMapper.CloseLogFileNow();
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreatePepToProteinMapFile", ex);
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
                OnWarningEvent("Error in CreateProteinModDetailsFile: " + ex.Message);
            }

            return success;
        }

        /// <summary>
        /// Create the protein modification details file
        /// </summary>
        /// <param name="phrpDataFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <param name="phrpResultType"></param>
        /// <returns>True if successful, false if an error</returns>
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
                var headerNames = new List<string>
                {
                    COLUMN_NAME_RESULTID,
                    COLUMN_NAME_PEPTIDE,
                    COLUMN_NAME_UNIQUE_SEQ_ID,
                    COLUMN_NAME_PROTEIN_NAME,
                    COLUMN_NAME_RESIDUE,
                    COLUMN_NAME_PROTEIN_RESIDUE_NUMBER,
                    COLUMN_NAME_RESIDUE_MOD_NAME,
                    COLUMN_NAME_PEPTIDE_RESIDUE_NUMBER,
                    COLUMN_NAME_MSGF_SPECPROB
                };

                writer.WriteLine(string.Join("\t", headerNames));

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
                    OnWarningEvent(msg);
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
                        OnWarningEvent("All PSMs map to reversed or scrambled proteins; the _ProteinMods.txt file is empty");
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
                SetErrorMessage("Error in CreateProteinModDetailsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Create the MTSPepToProteinMap file if missing, then create the protein modification details file
        /// </summary>
        /// <param name="baseName"></param>
        /// <param name="inputFile"></param>
        /// <param name="synOutputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="phrpResultType"></param>
        /// <param name="maximumAllowableMatchErrorPercentThreshold">
        /// Maximum percentage of peptides in the peptide to protein map file that are allowed to have not matched a protein in the FASTA file (value between 0 and 100)
        /// </param>
        /// <param name="matchErrorPercentWarningThreshold">
        /// When at least one peptide did not have a matched protein in the FASTA file, this threshold defines at what percent level a warning should be shown (value between 0 and 100)
        /// </param>
        /// <param name="fhtOutputFilePath"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <returns>True if successful, false if an error</returns>
        protected bool CreateProteinModsFileWork(
            string baseName,
            FileInfo inputFile,
            string synOutputFilePath,
            string outputDirectoryPath,
            PeptideHitResultTypes phrpResultType,
            double maximumAllowableMatchErrorPercentThreshold = 0.1,
            double matchErrorPercentWarningThreshold = 0,
            string fhtOutputFilePath = "",
            string mtsPepToProteinMapFilePath = ""
            )
        {
            bool success;

            if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath))
            {
                var baseNameFilePath = Path.Combine(inputFile.DirectoryName ?? string.Empty, baseName);
                mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseNameFilePath, outputDirectoryPath, mts: true);
            }

            var sourcePHRPDataFiles = new List<string>();

            if (!string.IsNullOrEmpty(fhtOutputFilePath))
            {
                sourcePHRPDataFiles.Add(fhtOutputFilePath);
            }

            if (!string.IsNullOrEmpty(synOutputFilePath))
            {
                sourcePHRPDataFiles.Add(synOutputFilePath);
            }

            if (sourcePHRPDataFiles.Count == 0)
            {
                SetErrorMessage("Cannot call CreatePepToProteinMapFile since sourcePHRPDataFiles is empty");
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                success = false;
            }
            else
            {
                if (File.Exists(mtsPepToProteinMapFilePath) && Options.UseExistingMTSPepToProteinMapFile)
                {
                    success = true;
                }
                else
                {
                    success = CreatePepToProteinMapFile(
                        sourcePHRPDataFiles,
                        mtsPepToProteinMapFilePath,
                        maximumAllowableMatchErrorPercentThreshold,
                        matchErrorPercentWarningThreshold);

                    if (!success)
                    {
                        // Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False
                        OnWarningEvent(WARNING_MESSAGE_SKIPPING_PROTEIN_MODS_FILE_CREATION + " since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                return true;
            }

            if (inputFile.Directory == null)
            {
                OnWarningEvent("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
            }
            else if (string.IsNullOrWhiteSpace(synOutputFilePath))
            {
                OnWarningEvent("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
            }
            else
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)),
                    outputDirectoryPath);

                // Create the Protein Mods file
                var modsFileCreated = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath, phrpResultType);

                if (!modsFileCreated)
                {
                    // Do not treat this as a fatal error
                    return true;
                }
            }

            return true;
        }

        /// <summary>
        /// Delete the file, ignoring any errors
        /// </summary>
        /// <param name="filePath"></param>
        protected void DeleteFileIgnoreErrors(string filePath)
        {
            try
            {
                if (File.Exists(filePath))
                {
                    Thread.Sleep(100);
                    File.Delete(filePath);
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        /// <summary>
        /// Expand the capacity of the list if required
        /// </summary>
        /// <remarks>
        /// If the list has fewer than largeListThreshold items, the capacity is left unchanged
        /// </remarks>
        /// <typeparam name="T"></typeparam>
        /// <param name="items"></param>
        /// <param name="countToAdd"></param>
        /// <param name="largeListThreshold"></param>
        protected void ExpandListIfRequired<T>(List<T> items, int countToAdd, int largeListThreshold = 1000000)
        {
            if (items.Count > largeListThreshold && items.Count + countToAdd > items.Capacity)
            {
                // .NET by default will double the size of the list to accommodate these new items
                // Instead, expand the list by 20% of the current size
                items.Capacity += Convert.ToInt32(items.Count / 5);
            }
        }

        /// <summary>
        /// Peptide to protein map results comparer
        /// </summary>
        private readonly IComparer<PepToProteinMapping> mPeptideSearchComparer = new PepToProteinMappingPeptideSearchComparer();

        /// <summary>
        /// Look for the peptide in the peptide to protein map list
        /// </summary>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="peptideToFind"></param>
        /// <returns>Index of the match if found, otherwise a negative number (based on a bitwise complement calculation)</returns>
        protected int FindFirstMatchInPepToProteinMapping(List<PepToProteinMapping> pepToProteinMapping, string peptideToFind)
        {
            // Use binary search to find this peptide in pepToProteinMapping
            var udtItemToFind = new PepToProteinMapping
            {
                Peptide = peptideToFind
            };

            var pepToProteinMapIndex = pepToProteinMapping.BinarySearch(udtItemToFind, mPeptideSearchComparer);

            if (pepToProteinMapIndex <= 0)
            {
                return pepToProteinMapIndex;
            }

            // Step Backward until the first match is found
            while (pepToProteinMapIndex > 0 && pepToProteinMapping[pepToProteinMapIndex - 1].Peptide == peptideToFind)
            {
                pepToProteinMapIndex--;
            }

            return pepToProteinMapIndex;
        }

        /// <summary>
        /// Return the application version (including the program date)
        /// </summary>
        public static string GetAppVersion()
        {
            return ProcessFilesOrDirectoriesBase.GetAppVersion(PROGRAM_DATE);
        }

        /// <summary>
        /// Remove the prefix and suffix residues from the given peptide sequence and return the result
        /// </summary>
        /// <param name="sequenceWithMods"></param>
        protected string GetCleanSequence(string sequenceWithMods)
        {
            return GetCleanSequence(sequenceWithMods, out _, out _);
        }

        /// <summary>
        /// Remove the prefix and suffix residues from the given peptide sequence and return the result
        /// </summary>
        /// <param name="sequenceWithMods"></param>
        /// <param name="prefix">Output: prefix residue</param>
        /// <param name="suffix">Output: suffix residue</param>
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
        /// <returns>True if columnIndex >= 0 and an integer value is present</returns>
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
                value = splitLine[columnIndex];
                return true;
            }

            value = valueIfMissing;
            return false;
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0 and an integer value is present</returns>
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
        /// <returns>True if success; does not catch errors; they will be thrown to the calling method if they occur</returns>
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
            mResultToSeqMapFile.WriteLine(StringUtilities.CollapseList(resultToSeqMapHeaders));

            // Initialize the SeqInfo file
            mSeqInfoFile = new StreamWriter(seqInfoFilePath, false);
            mSeqInfoFile.WriteLine(StringUtilities.CollapseList(seqInfoHeaders));

            // Initialize the ModDetails file
            mModDetailsFile = new StreamWriter(modDetailsFilePath);
            mModDetailsFile.WriteLine(StringUtilities.CollapseList(modDetailsHeaders));

            // Initialize the SeqToProtein map file
            mSeqToProteinMapFile = new StreamWriter(seqToProteinMapFilePath, false);
            mSeqToProteinMapFile.WriteLine(StringUtilities.CollapseList(seqToProteinMapHeaders));

            return true;
        }

        /// <summary>
        /// Return true if the protein name starts with or ends with a known reversed protein name tag
        /// </summary>
        /// <param name="proteinName"></param>
        /// <returns>
        /// True if the protein starts with reversed_, REV_, scrambled_, xxx_, or REV__, or if it ends with :reversed
        /// Otherwise, false
        /// </returns>
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

            if (proteinName.StartsWith("REV__", StringComparison.OrdinalIgnoreCase))
            {
                // Used by MSFragger
                return true;
            }

            if (proteinName.StartsWith("xxx.", StringComparison.OrdinalIgnoreCase))
            {
                // Used by InSpecT
                return true;
            }

            // ReSharper disable once ConvertIfStatementToReturnStatement
            if (proteinName.EndsWith(":reversed", StringComparison.OrdinalIgnoreCase))
            {
                // Used by X!Tandem
                return true;
            }

            return false;
        }

        /// <summary>
        /// Load settings from either a Key=Value parameter file or an XML-based parameter file
        /// </summary>
        /// <param name="parameterFilePath"></param>
        /// <returns>True if parameters were loaded (or if parameterFilePath is an empty string), false if an error</returns>
        protected bool LoadParameterFileSettings(string parameterFilePath)
        {
            try
            {
                if (string.IsNullOrWhiteSpace(parameterFilePath))
                {
                    // No parameter file specified; nothing to load
                    return true;
                }

                var parameterFile = new FileInfo(parameterFilePath);

                if (parameterFile.Extension.Equals(".xml", StringComparison.OrdinalIgnoreCase))
                    return LoadParameterFileSettingsXML(parameterFilePath);

                // Read settings from a Key=Value parameter file

                var assemblyName = Assembly.GetEntryAssembly()?.GetName().Name;

                var cmdLineParser = new CommandLineParser<PHRPOptions>(assemblyName, GetAppVersion())
                {
                    ProgramInfo = "This program converts search results from various MS/MS identification tools " +
                                  "into a series of tab-delimited text files that organize the data in a similar format for each tool.",
                    ContactInfo = "Program written by Matthew Monroe for PNNL (Richland, WA)"
                };

                var args = new List<string>
                {
                    "-ParamFile:" + PathUtils.PossiblyQuotePath(parameterFile.FullName)
                };

                var result = cmdLineParser.ParseArgs(args.ToArray());
                var options = result.ParsedResults;
                if (!result.Success || !options.Validate())
                {
                    SetErrorMessage(string.Format("Error in LoadParameterFileSettings reading settings from parameter file {0}", parameterFile.FullName));
                    SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile);
                    return false;
                }

                Options.UpdateAll(options);

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadParameterFileSettings", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile);
                return false;
            }
        }

        /// <summary>
        /// Load settings from an XML-based parameter file
        /// </summary>
        /// <param name="parameterFilePath"></param>
        /// <returns>True if parameters were loaded (or if parameterFilePath is an empty string), false if an error</returns>
        private bool LoadParameterFileSettingsXML(string parameterFilePath)
        {
            const string OPTIONS_SECTION = "PeptideHitResultsProcessorOptions";

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
                        // Section OPTIONS_SECTION was not found in the parameter file
                        SetErrorMessage("The node '<section name=\"" + OPTIONS_SECTION + "\"> was not found in the parameter file: " + parameterFilePath);
                        SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile);
                        return false;
                    }

                    Options.MassCorrectionTagsFilePath = settingsFile.GetParam(OPTIONS_SECTION, "MassCorrectionTagsFilePath", Options.MassCorrectionTagsFilePath);
                    Options.ModificationDefinitionsFilePath = settingsFile.GetParam(OPTIONS_SECTION, "ModificationDefinitionsFilePath", Options.ModificationDefinitionsFilePath);
                    Options.SearchToolParameterFilePath = settingsFile.GetParam(OPTIONS_SECTION, "SearchToolParameterFilePath", Options.SearchToolParameterFilePath);

                    Options.CreateModificationSummaryFile = settingsFile.GetParam(OPTIONS_SECTION, "CreateModificationSummaryFile", Options.CreateModificationSummaryFile);

                    Options.CreateProteinModsFile = settingsFile.GetParam(OPTIONS_SECTION, "CreateProteinModsFile", Options.CreateProteinModsFile);
                    Options.FastaFilePath = settingsFile.GetParam(OPTIONS_SECTION, "FastaFilePath", Options.FastaFilePath);
                    Options.ProteinModsFileIncludesReversedProteins = settingsFile.GetParam(OPTIONS_SECTION, "ProteinModsFileIncludesReversedProteins", Options.ProteinModsFileIncludesReversedProteins);
                    Options.UseExistingMTSPepToProteinMapFile = settingsFile.GetParam(OPTIONS_SECTION, "UseExistingMTSPepToProteinMapFile", Options.UseExistingMTSPepToProteinMapFile);

                    var leftResidueRegEx = Options.EnzymeMatchSpec.LeftResidueRegEx;
                    var rightResidueRegEx = Options.EnzymeMatchSpec.RightResidueRegEx;

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

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadParameterFileSettingsXML", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile);
                return false;
            }
        }

        /// <summary>
        /// Load the PeptideToProteinMap information
        /// </summary>
        /// <param name="pepToProteinMapFilePath">File to read</param>
        /// <param name="pepToProteinMapping">Output parameter: peptide to protein mapping (calling method must pre-initialize the list)</param>
        /// <param name="headerLine">Output parameter: Header line text</param>
        /// <returns>True if successful, false if an error</returns>
        protected bool LoadPeptideToProteinMapInfo(
            string pepToProteinMapFilePath,
            List<PepToProteinMapping> pepToProteinMapping,
            out string headerLine)
        {
            headerLine = string.Empty;

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
                using var reader = new StreamReader(new FileStream(pepToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

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
                            headerLine = dataLine;
                        }
                        else
                        {
                            var pepToProteinMappingEntry = new PepToProteinMapping
                            {
                                Peptide = splitLine[0],
                                Protein = splitLine[1]
                            };
                            int.TryParse(splitLine[2], out pepToProteinMappingEntry.ResidueStart);
                            int.TryParse(splitLine[3], out pepToProteinMappingEntry.ResidueEnd);

                            ExpandListIfRequired(pepToProteinMapping, 1);

                            pepToProteinMapping.Add(pepToProteinMappingEntry);
                        }
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the Peptide to Protein Map File (" + Path.GetFileName(pepToProteinMapFilePath) + ")", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        /// <summary>
        /// Contact the database to lookup dataset IDs by dataset name
        /// </summary>
        /// <param name="datasetNames"></param>
        /// <returns>Dictionary where keys are dataset names and values are dataset IDs</returns>
        protected Dictionary<string, int> LookupDatasetIDs(IEnumerable<string> datasetNames)
        {
            var datasetIDs = new Dictionary<string, int>();

            if (string.IsNullOrWhiteSpace(Options.DMSConnectionString) ||
                Options.DMSConnectionString.Equals("false", StringComparison.OrdinalIgnoreCase))
            {
                return datasetIDs;
            }

            try
            {
                var domain = Domain.GetComputerDomain();
                if (!domain.Name.EndsWith("pnl.gov", StringComparison.OrdinalIgnoreCase))
                    return datasetIDs;
            }
            catch (Exception)
            {
                // Computer is not joined to a domain
                return datasetIDs;
            }

            try
            {
                var quotedDatasetNames = new StringBuilder();
                foreach (var item in datasetNames)
                {
                    var optionalComma = quotedDatasetNames.Length == 0 ? string.Empty : ", ";
                    quotedDatasetNames.AppendFormat("{0}'{1}'", optionalComma, item);
                }

                var sqlQuery = string.Format(
                    "SELECT Dataset, ID " +
                    "FROM V_Dataset_Export " +
                    "WHERE dataset in ({0})", quotedDatasetNames);

                var connectionStringToUse = DbToolsFactory.AddApplicationNameToConnectionString(Options.DMSConnectionString, "PeptideHitResultsProcessor");

                var dbTools = DbToolsFactory.GetDBTools(connectionStringToUse);

                var success = dbTools.GetQueryResults(sqlQuery, out var results);

                if (!success)
                {
                    OnWarningEvent("DbTools returned false querying V_Dataset_Export to look up dataset IDs by dataset name");
                    return datasetIDs;
                }

                foreach (var result in results)
                {
                    var dataset = result[0];
                    var datasetId = result[1];

                    datasetIDs.Add(dataset, int.Parse(datasetId));
                }
            }
            catch (Exception ex)
            {
                OnWarningEvent("Error looking up dataset IDs by dataset name: " + ex.Message);
            }

            return datasetIDs;
        }

        /// <summary>
        /// Raise event ProgressComplete
        /// </summary>
        protected void OperationComplete()
        {
            ProgressComplete?.Invoke();
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">PSM tool results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <returns>True if successful, False if failure</returns>
        // ReSharper disable once UnusedMember.Global
        public bool ProcessFile(string inputFilePath, string outputDirectoryPath)
        {
            return ProcessFile(inputFilePath, outputDirectoryPath, string.Empty);
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">PSM tool results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file (either legacy XML-based or Key=Value)</param>
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
            newSuffix ??= string.Empty;

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

            return string.IsNullOrWhiteSpace(originalFile.DirectoryName)
                ? newFileName
                : Path.Combine(originalFile.DirectoryName, newFileName);
        }

        /// <summary>
        /// Report an error by raising event ErrorEvent; store the error in mErrorMessage
        /// </summary>
        /// <param name="errMsg"></param>
        /// <param name="throwException">If True, throw an exception</param>
        /// <param name="ex"></param>
        protected void ReportError(string errMsg, bool throwException = false, Exception ex = null)
        {
            SetErrorMessage(errMsg, ex);

            if (!throwException)
            {
                return;
            }

            if (ex == null)
            {
                throw new Exception(errMsg);
            }

            throw new Exception(errMsg, ex);
        }

        /// <summary>
        /// Reset the mass correction tags and modification definitions
        /// </summary>
        /// <remarks>
        /// Reads data from Options.MassCorrectionTagsFilePath or Options.ModificationDefinitionsFilePath if defined,
        /// otherwise, resets to default values
        /// </remarks>
        public bool ResetMassCorrectionTagsAndModificationDefinitions()
        {
            var fileNotFound = false;

            // Note: If mMassCorrectionTagsFilePath is blank, the mass correction tags will be reset to the defaults and success will be True
            var massCorrectionTagsSuccess = mPeptideMods.ReadMassCorrectionTagsFile(massCorrectionFilePath, out var massCorrectionTagsFileNotFound);

            if (!massCorrectionTagsSuccess)
            {
                if (massCorrectionTagsFileNotFound)
                {
                    SetErrorCode(PHRPErrorCode.MassCorrectionTagsFileNotFound);
                    OnWarningEvent("Mass Correction Tags file not found: " + Options.MassCorrectionTagsFilePath);
                }
                else
                {
                    SetErrorCode(PHRPErrorCode.ErrorReadingMassCorrectionTagsFile);
                }
            }


            // Note: If mModificationDefinitionsFilePath is blank, the modifications will be cleared and success will be True
            var modDefsSuccess = mPeptideMods.ReadModificationDefinitionsFile(modificationDefinitionsFilePath, out var modDefsFileNotFound);

            {
                if (fileNotFound)
                {
                    SetErrorCode(PHRPErrorCode.ModificationDefinitionFileNotFound);
                    OnWarningEvent("File not found: " + Options.ModificationDefinitionsFilePath);
                }
                else
                {
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
            }

            return success;
        }

        /// <summary>
        /// Reset percent complete to 0 and update mProgressStepDescription to the given description
        /// </summary>
        /// <param name="progressStepDescription"></param>
        /// <param name="echoToConsole"></param>
        protected void ResetProgress(string progressStepDescription, bool echoToConsole = false)
        {
            mProgressStepDescription = progressStepDescription;
            mProgressPercentComplete = 0;
            ProgressReset?.Invoke();

            if (!echoToConsole)
            {
                return;
            }

            Console.WriteLine();
            Console.WriteLine();
            Console.WriteLine(ProgressStepDescription);
        }

        /// <summary>
        /// Looks for fileNameOrPath in the current working directory
        /// If not found, looks in sourceDirectoryPath
        /// If searchParentDirectory is true, also looks in the parent directory of both the source file directory and the working directory
        /// </summary>
        /// <param name="sourceDirectoryPath">Path to the directory containing the input file</param>
        /// <param name="fileNameOrPath">File to find (either filename or full file path)</param>
        /// <param name="searchParentDirectory">
        /// If true and the file is not found in the working directory or the source directory, also examine the parent directory of each
        /// </param>
        /// <returns>The path to the file if found, or fileNameOrPath if not found</returns>
        public static string ResolveFilePath(string sourceDirectoryPath, string fileNameOrPath, bool searchParentDirectory = true)
        {
            if (File.Exists(fileNameOrPath))
            {
                return fileNameOrPath;
            }

            var fileName = Path.GetFileName(fileNameOrPath);
            if (string.IsNullOrWhiteSpace(fileName))
                return fileNameOrPath;

            var sourceDirectoryCandidateFile = new FileInfo(Path.Combine(sourceDirectoryPath, fileName));
            if (sourceDirectoryCandidateFile.Exists)
            {
                return sourceDirectoryCandidateFile.FullName;
            }

            var workingDirectoryCandidateFile = new FileInfo(fileName);
            if (workingDirectoryCandidateFile.Exists)
            {
                return workingDirectoryCandidateFile.FullName;
            }

            if (!searchParentDirectory)
            {
                return fileNameOrPath;
            }

            var parentDirectories = new List<DirectoryInfo>();

            if (sourceDirectoryCandidateFile.Directory != null)
            {
                parentDirectories.Add(sourceDirectoryCandidateFile.Directory.Parent);
            }

            if (workingDirectoryCandidateFile.Directory != null)
            {
                parentDirectories.Add(workingDirectoryCandidateFile.Directory.Parent);
            }

            // ReSharper disable once ForeachCanBePartlyConvertedToQueryUsingAnotherGetEnumerator
            foreach (var directory in parentDirectories)
            {
                var parentCandidateFile = new FileInfo(Path.Combine(directory.FullName, fileName));

                if (parentCandidateFile.Exists)
                {
                    return parentCandidateFile.FullName;
                }
            }

            return fileNameOrPath;
        }

        /// <summary>
        /// Round a value (stored as text) to the given number of digits after the decimal point
        /// </summary>
        /// <param name="valueToRound"></param>
        /// <param name="digitsAfterDecimal"></param>
        protected string RoundValue(string valueToRound, byte digitsAfterDecimal = 4)
        {
            if (!double.TryParse(valueToRound, out var value))
            {
                return valueToRound;
            }

            return PRISM.StringUtilities.DblToString(value, digitsAfterDecimal, 0.00001);
        }

        /// <summary>
        /// Create the modification summary file
        /// </summary>
        /// <param name="modificationSummaryFilePath"></param>
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

            writer.WriteLine(StringUtilities.CollapseList(headerNames));

            for (var index = 0; index <= mPeptideMods.ModificationCount - 1; index++)
            {
                var modDef = mPeptideMods.GetModificationByIndex(index);
                if (modDef.OccurrenceCount <= 0 && modDef.UnknownModAutoDefined)
                    continue;

                var data = new List<string>
                {
                    modDef.ModificationSymbol.ToString(),
                    PRISM.StringUtilities.DblToString(modDef.ModificationMass, 6),
                    modDef.TargetResidues,
                    ModificationDefinition.ModificationTypeToModificationSymbol(modDef.ModificationType).ToString(),
                    modDef.MassCorrectionTag,
                    modDef.OccurrenceCount.ToString()
                };

                writer.WriteLine(StringUtilities.CollapseList(data));
            }
        }

        /// <summary>
        /// Append an entry to the _ResultToSeqMap.txt and _SeqInfo.txt files
        /// </summary>
        /// <remarks>Call InitializeOutputFiles once prior to calling this method</remarks>
        /// <param name="searchResult"></param>
        /// <param name="updateResultToSeqMapFile">Set to True only for the first protein of each peptide in each group</param>
        protected void SaveResultsFileEntrySeqInfo(SearchResultsBaseClass searchResult, bool updateResultToSeqMapFile)
        {
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
                mResultToSeqMapFile.WriteLine(StringUtilities.CollapseList(seqMapData));

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
                    mSeqInfoFile.WriteLine(StringUtilities.CollapseList(seqInfoData));

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

                            mModDetailsFile.WriteLine(StringUtilities.CollapseList(modDetailsData));
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

                mSeqToProteinMapFile.WriteLine(StringUtilities.CollapseList(seqToProteinData));
            }
        }

        /// <summary>
        /// Update the error code stored in mErrorCode
        /// </summary>
        /// <param name="newErrorCode"></param>
        protected void SetErrorCode(PHRPErrorCode newErrorCode)
        {
            SetErrorCode(newErrorCode, false);
        }

        /// <summary>
        /// Update the error code stored in mErrorCode
        /// </summary>
        /// <param name="newErrorCode"></param>
        /// <param name="leaveExistingErrorCodeUnchanged">If true and mErrorCode is already non-zero, leave unchanged</param>
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

        /// <summary>
        /// Update the error message stored in mErrorMessage
        /// </summary>
        /// <remarks>
        /// If ex is not null, appends the exception message to the error message (if not already present)
        /// </remarks>
        /// <param name="message"></param>
        /// <param name="ex"></param>
        protected void SetErrorMessage(string message, Exception ex = null)
        {
            message ??= string.Empty;

            if (ex != null && message.IndexOf(ex.Message, StringComparison.Ordinal) < 0)
            {
                message += ": " + ex.Message;
            }

            mErrorMessage = message;

            if (message.Length > 0)
            {
                OnErrorEvent(message, ex);
            }
        }

        /// <summary>
        /// Show a warning periodically
        /// </summary>
        /// <param name="warningCount">Number of times this warning has been encountered</param>
        /// <param name="thresholdCountAlwaysShow">Always show the warning up to this many times</param>
        /// <param name="warningMessage">Warning message</param>
        /// <remarks>As the warning count surpasses, 1000, 10000, etc., the warning is shown less frequently</remarks>
        protected void ShowPeriodicWarning(int warningCount, int thresholdCountAlwaysShow, string warningMessage)
        {
            if (warningCount <= thresholdCountAlwaysShow ||
                warningCount < 1000 && warningCount % 100 == 0 ||
                warningCount < 10000 && warningCount % 1000 == 0 ||
                warningCount < 100000 && warningCount % 10000 == 0 ||
                warningCount < 1000000 && warningCount % 100000 == 0)
            {
                OnWarningEvent(warningMessage);
            }
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

        /// <summary>
        /// Update the peptide to protein mapping at the given index
        /// </summary>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="index"></param>
        /// <param name="peptide"></param>
        protected void UpdatePepToProteinMapPeptide(List<PepToProteinMapping> pepToProteinMapping, int index, string peptide)
        {
            var udtItem = pepToProteinMapping[index];
            udtItem.Peptide = peptide;
            pepToProteinMapping[index] = udtItem;
        }

        /// <summary>
        /// Update the progress percent complete
        /// </summary>
        /// <param name="percentComplete">Value between 0 and 100</param>
        protected void UpdateProgress(float percentComplete)
        {
            UpdateProgress(ProgressStepDescription, percentComplete);
        }

        /// <summary>
        /// Update the progress percent complete
        /// </summary>
        /// <param name="progressStepDescription">Current processing step description</param>
        /// <param name="percentComplete">Value between 0 and 100</param>
        protected void UpdateProgress(string progressStepDescription, float percentComplete)
        {
            mProgressStepDescription = progressStepDescription;

            if (percentComplete < 0)
            {
                mProgressPercentComplete = 0;
            }
            else if (percentComplete > 100)
            {
                mProgressPercentComplete = 100;
            }
            else
            {
                mProgressPercentComplete = percentComplete;
            }

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
                using var reader = new StreamReader(new FileStream(fileInfo.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

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

                if (!dataFound)
                {
                    errorMessage = fileDescription + " is empty (no data)";
                }

                return dataFound;
            }
            catch (Exception)
            {
                errorMessage = "Exception validating " + fileDescription + " file";
                return false;
            }
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

        /// <summary>
        /// Determine the percentage of results in peptideToProteinMapFilePath that have protein name __NoMatch__
        /// </summary>
        /// <param name="peptideToProteinMapFilePath">Peptide to protein map file</param>
        /// <param name="ignorePeptideToProteinMapperErrors">When true, return true even if one or more peptides did not match to a known protein</param>
        /// <param name="maximumAllowableMatchErrorPercentThreshold">
        /// Maximum percentage of peptides in the peptide to protein map file that are allowed to have not matched a protein in the FASTA file (value between 0 and 100)
        /// This is typically 0.1, but for MaxQuant we set this to 50
        /// </param>
        /// <param name="matchErrorPercentWarningThreshold">
        /// When at least one peptide did not have a matched protein in the FASTA file, this threshold defines at what percent level a warning should be shown (value between 0 and 100)
        /// This is typically 0, but for MaxQuant and MSFragger we set it to 5
        /// </param>
        /// <returns>
        /// True if the required percentage of peptides matched a known protein, false if they did not and ignorePeptideToProteinMapperErrors is false
        /// </returns>
        private bool ValidatePeptideToProteinMapResults(
            string peptideToProteinMapFilePath,
            bool ignorePeptideToProteinMapperErrors,
            double maximumAllowableMatchErrorPercentThreshold,
            double matchErrorPercentWarningThreshold)
        {
            try
            {
                var peptideCount = ValidatePeptideToProteinMapResultsWork(peptideToProteinMapFilePath, out var peptideCountNoMatch);

                if (peptideCount == 0)
                {
                    SetErrorMessage("Peptide to protein mapping file is empty: " + peptideToProteinMapFilePath);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                if (peptideCountNoMatch == 0)
                {
                    return true;
                }

                // Value between 0 and 100
                var errorPercent = peptideCountNoMatch / (double)peptideCount * 100.0;

                var message = string.Format("{0:0.00}% of the entries ({1:N0} / {2:N0}) in the peptide to protein map file ({3}) " +
                                            "did not match to a protein in the FASTA file ({4})",
                    errorPercent, peptideCountNoMatch, peptideCount,
                    Path.GetFileName(peptideToProteinMapFilePath),
                    Path.GetFileName(Options.FastaFilePath));

                if (ignorePeptideToProteinMapperErrors || errorPercent < maximumAllowableMatchErrorPercentThreshold)
                {
                    if (errorPercent >= matchErrorPercentWarningThreshold)
                    {
                        OnWarningEvent(message);
                    }
                    else
                    {
                        Console.WriteLine();
                        OnStatusEvent(message);
                    }

                    return true;
                }

                SetErrorMessage(message);
                return false;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ValidatePeptideToProteinMapResults", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        private static int ValidatePeptideToProteinMapResultsWork(string peptideToProteinMapFilePath, out int peptideCountNoMatch)
        {
            var peptideCount = 0;
            var linesRead = 0;
            var splitChars = new[] { '\t' };

            var lastPeptide = string.Empty;

            peptideCountNoMatch = 0;

            using var reader = new StreamReader(new FileStream(peptideToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

            while (!reader.EndOfStream)
            {
                var lineIn = reader.ReadLine();
                linesRead++;

                if (linesRead <= 1 || string.IsNullOrEmpty(lineIn))
                    continue;

                var splitLine = lineIn.Split(splitChars, 2);
                if (splitLine.Length == 0)
                    continue;

                if (splitLine[0] != lastPeptide)
                {
                    peptideCount++;
                    lastPeptide = splitLine[0];
                }

                if (lineIn.Contains(PROTEIN_NAME_NO_MATCH))
                {
                    peptideCountNoMatch++;
                }
            }

            return peptideCount;
        }

        /// <summary>
        /// Assures that the _MSGF.txt file exists in the output directory, provided it exists in the input directory
        /// </summary>
        /// <remarks>Copies the file to the output directory if missing</remarks>
        /// <param name="phrpDataFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
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

                var sourcePath = string.IsNullOrWhiteSpace(phrpDataFile.DirectoryName)
                    ? msgfFileName
                    : Path.Combine(phrpDataFile.DirectoryName, msgfFileName);

                var targetPath = Path.Combine(outputDirectory.FullName, msgfFileName);

                if (File.Exists(sourcePath) && !File.Exists(targetPath))
                {
                    File.Copy(sourcePath, targetPath);
                }
            }
            catch (Exception ex)
            {
                OnWarningEvent("Error in ValidatePHRPReaderSupportFiles: " + ex.Message);
            }
        }

        /// <summary>
        /// Verify that the FASTA file exists and has amino acid based proteins
        /// </summary>
        /// <param name="fastaFilePath"></param>
        /// <returns>True if found and valid, otherwise false</returns>
        protected bool ValidateProteinFastaFile(string fastaFilePath)
        {
            var success = ValidateProteinFastaFile(fastaFilePath, out var warningMessage);

            if (!success)
            {
                OnWarningEvent(warningMessage);
            }

            return success;
        }

        /// <summary>
        /// Verify that the FASTA file exists and has amino acid based proteins
        /// </summary>
        /// <remarks>
        /// <para>
        /// Proteins where at least 10% of the residues are not A, T, C, or G are considered valid proteins
        /// </para>
        /// <para>
        /// Proteins where over 95% of the residues are A, T, C, or G are considered invalid proteins
        /// </para>
        /// <para>
        /// After examining the first 500 proteins, if there are more invalid proteins than valid proteins, returns false
        /// </para>
        /// </remarks>
        /// <param name="fastaFilePath"></param>
        /// <param name="warningMessage"></param>
        /// <returns>True if found and valid, otherwise false</returns>
        public static bool ValidateProteinFastaFile(string fastaFilePath, out string warningMessage)
        {
            // This RegEx looks for standard amino acids, skipping A, T, C, and G
            var definiteAminoAcidMatcher = new Regex("[DEFHIKLMNPQRSVWY]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            // This RegEx looks for A, T, C, and G
            var potentialNucleicAcidMatcher = new Regex("[ATCG]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            // This matches any letter
            var letterMatcher = new Regex("[A-Z]", RegexOptions.Compiled | RegexOptions.IgnoreCase);

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
                    warningMessage = "FASTA file not found: " + fastaFilePath;
                    return false;
                }

                var fastaFile = new ProteinFileReader.FastaFileReader();
                if (!fastaFile.OpenFile(fastaFilePath))
                {
                    Console.WriteLine();
                    warningMessage = "Error opening the FASTA file: " + fastaFilePath;
                    return false;
                }

                // Read the first 500 proteins and confirm that each contains amino acid residues
                while (fastaFile.ReadNextProteinEntry())
                {
                    var definiteAminoAcidCount = definiteAminoAcidMatcher.Matches(fastaFile.ProteinSequence).Count;
                    var potentialNucleicAcidCount = potentialNucleicAcidMatcher.Matches(fastaFile.ProteinSequence).Count;
                    var letterCount = letterMatcher.Matches(fastaFile.ProteinSequence).Count;

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
                    warningMessage = "FASTA file contains Nucleic Acids, not Amino Acids: " + Path.GetFileName(fastaFilePath);
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
            var dataValues = new List<string>();

            foreach (var mod in reader.CurrentPSM.ModifiedResidues)
            {
                var residueLocInProtein = pepToProteinMapping[pepToProteinMapIndex].ResidueStart + mod.ResidueLocInPeptide - 1;
                string residue;

                if (StringUtilities.IsLetterAtoZ(mod.Residue))
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
                    dataValues.Clear();

                    dataValues.Add(reader.CurrentPSM.ResultID.ToString());
                    dataValues.Add(reader.CurrentPSM.Peptide);
                    dataValues.Add(reader.CurrentPSM.SeqID.ToString());
                    dataValues.Add(pepToProteinMapping[pepToProteinMapIndex].Protein);
                    dataValues.Add(residue);
                    dataValues.Add(residueLocInProtein.ToString());
                    dataValues.Add(mod.ModDefinition.MassCorrectionTag);
                    dataValues.Add(mod.ResidueLocInPeptide.ToString());
                    dataValues.Add(reader.CurrentPSM.MSGFSpecEValue);

                    writer.WriteLine(string.Join("\t", dataValues));
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

                x.ModName ??= string.Empty;

                y.ModName ??= string.Empty;

                return string.CompareOrdinal(x.ModName, y.ModName);
            }
        }

        /// <summary>
        /// Peptide to protein map item comparer
        /// </summary>
        protected class PepToProteinMappingComparer : IComparer<PepToProteinMapping>
        {
            /// <summary>
            /// Compare to peptide to protein map values, sorting first on peptide sequence and then on protein name
            /// </summary>
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
