// This class reads in an DIA-NN report.tsv file and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PeptideHitResultsProcessor.Data;
using PeptideHitResultsProcessor.SearchToolResults;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads in DIA-NN report.tsv file and creates
    /// a tab-delimited text file with the data
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads DIA-NN results file report.tsv
    /// </para>
    /// <para>
    /// 2) It calls CreateSynResultsFile to create the _syn.txt file
    /// </para>
    /// <para>
    /// 3) ParseDiaNNResultsFileHeaderLine reads the header line to determine the column mapping
    ///      columnMapping = new Dictionary of DiaNNReportFileColumns, int
    /// </para>
    /// <para>
    /// 4) ParseDiaNNResultsFileEntry reads each data line and stores in an instance of DiaNNSearchResult, which is a private structure
    ///    The data is stored in a list
    ///      searchResultsUnfiltered = new List of DiaNNSearchResult
    /// </para>
    /// <para>
    /// 5) Once the entire .tsv has been read, searchResultsUnfiltered is sorted by scan, charge, and ascending QValue
    /// </para>
    /// <para>
    /// 6) StoreSynMatches stores filter-passing values in a new list
    ///      filteredSearchResults = new List of DiaNNSearchResult
    /// </para>
    /// <para>
    /// 7) SortAndWriteFilteredSearchResults performs one more sort, then writes out to disk
    ///    Sorts ascending Expectation value, Scan, Peptide, and Protein
    /// </para>
    /// </remarks>
    public class DiaNNResultsProcessor : MultiDatasetResultsProcessor
    {
        // ReSharper disable CommentTypo

        // Ignore Spelling: acetylated, Da, Carbamidomethyl, expectscore, Hyperscore, massdiff
        // Ignore Spelling: Nextscore, Prev, proline, scannum, sp, tryptic, txt

        // ReSharper restore CommentTypo

        /// <summary>
        /// Constructor
        /// </summary>
        public DiaNNResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "April 6, 2023";

            mModificationMassByName = new Dictionary<string, double>();

            mPeptideCleavageStateCalculator = new PeptideCleavageStateCalculator();
        }

        /// <summary>
        /// Default Q-Value threshold to use when creating the synopsis file
        /// </summary>
        public const double DEFAULT_QVALUE_THRESHOLD = 0.10;

        /// <summary>
        /// Default confidence score (CScore) threshold to use when creating the synopsis file
        /// </summary>
        public const double DEFAULT_CONFIDENCE_SCORE_THRESHOLD = 0.25;

        /// <summary>
        /// DIA-NN tool name
        /// </summary>
        public const string TOOL_NAME = "DIA-NN";

        /// <summary>
        /// Default Hyperscore threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Hyperscore is over the threshold, or if its expectation value is below the threshold
        /// </remarks>
        public const int DEFAULT_HYPERSCORE_THRESHOLD = 20;

        /// <summary>
        /// N-terminus symbol used by DIA-NN
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_DiaNN = "-";

        /// <summary>
        /// C-terminus symbol used by DIA-NN
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_DiaNN = "-";

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        /// <summary>
        /// These columns correspond to the DIA-NN report.tsv file
        /// </summary>
        private enum DiaNNReportFileColumns
        {
            /// <summary>
            /// Undefined
            /// </summary>
            Undefined = -1,

            /// <summary>
            /// Dataset .mzML file path
            /// </summary>
            DatasetFile = 0,

            /// <summary>
            /// Dataset Name (from the Run column)
            /// </summary>
            /// <remarks>
            /// This is an abbreviated dataset name, as assigned by the analysis manager
            /// </remarks>
            DatasetName = 1,

            /// <summary>
            /// Protein group name
            /// </summary>
            ProteinGroup = 2,

            /// <summary>
            /// Protein names (from the FASTA file)
            /// </summary>
            ProteinIDs = 3,

            /// <summary>
            /// Protein names, as determined by DIA-NN, typically corresponding to UniProt Name
            /// </summary>
            ///
            ProteinNames = 4,

            /// <summary>
            /// Gene names associated with the peptide
            /// </summary>
            GeneNames = 5,

            /// <summary>
            /// Protein Group Quantity
            /// </summary>
            ProteinGroupQuantity = 6,

            /// <summary>
            /// Protein Group Normalized
            /// </summary>
            ProteinGroupNormalized = 7,

            /// <summary>
            /// Protein Group Max LFQ
            /// </summary>
            ProteinGroupMaxLFQ = 8,

            /// <summary>
            /// Genes Quantity
            /// </summary>
            GenesQuantity = 9,

            /// <summary>
            /// Genes Normalized
            /// </summary>
            GenesNormalized = 10,

            /// <summary>
            /// Genes Max LFQ
            /// </summary>
            GenesMaxLFQ = 11,

            /// <summary>
            /// Genes Max LFQUnique
            /// </summary>
            GenesMaxLFQUnique = 12,

            // ReSharper disable CommentTypo

            /// <summary>
            /// Peptide sequence with modifications
            /// </summary>
            /// <remarks>
            /// Examples:
            ///   AAALEAM(UniMod:35)K
            ///   AAEAHVDAHYYEQNEQPTGTC(UniMod:4)AAC(UniMod:4)ITGDNR
            /// </remarks>
            ModifiedSequence = 13,

            // ReSharper restore CommentTypo

            /// <summary>
            /// Peptide sequence without modifications
            /// </summary>
            StrippedSequence = 14,

            /// <summary>
            /// Precursor sequence (unused by PHRP)
            /// </summary>
            PrecursorId = 15,

            /// <summary>
            /// Precursor charge state
            /// </summary>
            PrecursorCharge = 16,

            /// <summary>
            /// QValue
            /// </summary>
            QValue = 17,

            /// <summary>
            /// PEP (posterior error probability)
            /// </summary>
            PEP = 18,

            /// <summary>
            /// Global QValue
            /// </summary>
            GlobalQValue = 19,

            /// <summary>
            /// Protein QValue
            /// </summary>
            ProteinQValue = 20,

            /// <summary>
            /// Protein Group QValue
            /// </summary>
            ProteinGroupQValue = 21,

            /// <summary>
            /// Global Protein Group QValue
            /// </summary>
            GlobalProteinGroupQValue = 22,

            /// <summary>
            /// Gene Group QValue
            /// </summary>
            GeneGroupQValue = 23,

            /// <summary>
            /// Translated QValue
            /// </summary>
            TranslatedQValue = 24,

            /// <summary>
            /// 1 if the peptide is specific to a given protein or gene, otherwise 0
            /// </summary>
            Proteotypic = 25,

            /// <summary>
            /// Precursor Quantity
            /// </summary>
            PrecursorQuantity = 26,

            /// <summary>
            /// Precursor Normalized
            /// </summary>
            PrecursorNormalized = 27,

            /// <summary>
            /// Precursor Translated
            /// </summary>
            PrecursorTranslated = 28,

            /// <summary>
            /// Translated Quality
            /// </summary>
            TranslatedQuality = 29,

            /// <summary>
            /// MS1 Translated
            /// </summary>
            MS1Translated = 30,

            /// <summary>
            /// Quantity Quality
            /// </summary>
            QuantityQuality = 31,

            /// <summary>
            /// Elution Time (aka retention time)
            /// </summary>
            RT = 32,

            /// <summary>
            /// Elution Time Start
            /// </summary>
            RTStart = 33,

            /// <summary>
            /// Elution Time Stop
            /// </summary>
            RTStop = 34,

            /// <summary>
            /// Indexed retention time (iRT)
            /// </summary>
            IndexedRT = 35,

            /// <summary>
            /// Predicted retention time (unused by PHRP)
            /// </summary>
            PredictedRT = 36,

            /// <summary>
            /// Predicted indexed retention time (unused by PHRP)
            /// </summary>
            PredictedIndexedRT = 37,

            /// <summary>
            /// First protein description (unused by PHRP)
            /// </summary>
            FirstProteinDescription = 38,

            /// <summary>
            /// Spectral Library QValue
            /// </summary>
            LibQValue = 39,

            /// <summary>
            /// Spectral Library Protein Group QValue
            /// </summary>
            LibProteinGroupQValue = 40,

            /// <summary>
            /// MS1 Profile Correlation
            /// </summary>
            MS1ProfileCorr = 41,

            /// <summary>
            /// MS1 Area
            /// </summary>
            MS1Area = 42,

            /// <summary>
            /// Evidence (score)
            /// </summary>
            /// <remarks>Higher values are better</remarks>
            Evidence = 43,

            /// <summary>
            /// Spectrum Similarity
            /// </summary>
            SpectrumSimilarity = 44,

            /// <summary>
            /// Averagine
            /// </summary>
            Averagine = 45,

            /// <summary>
            /// Mass Evidence
            /// </summary>
            /// <remarks>Higher values are better</remarks>
            MassEvidence = 46,

            /// <summary>
            /// CScore
            /// </summary>
            CScore = 47,

            /// <summary>
            /// Decoy Evidence
            /// </summary>
            DecoyEvidence = 48,

            /// <summary>
            /// Decoy CScore
            /// </summary>
            DecoyCScore = 49,

            /// <summary>
            /// Fragmentation Ion Abundances
            /// </summary>
            FragmentQuantRaw = 50,

            /// <summary>
            /// Corrected Fragmentation Ion Abundances
            /// </summary>
            FragmentQuantCorrected = 51,

            /// <summary>
            /// Fragmentation ion correlations
            /// </summary>
            FragmentCorrelations = 52,

            /// <summary>
            /// MS/MS Scan Number
            /// </summary>
            MS2Scan = 53,

            /// <summary>
            /// Ion Mobility
            /// </summary>
            IonMobility = 54,

            /// <summary>
            /// Indexed Ion Mobility
            /// </summary>
            IndexedIonMobility = 55,

            /// <summary>
            /// Predicted Ion Mobility
            /// </summary>
            PredictedIonMobility = 56,

            /// <summary>
            /// Predicted Indexed Ion Mobility
            /// </summary>
            PredictedIndexedIonMobility = 57
        }

        private struct DiaNNModInfo
        {
            public double ModMass;
            public char ResidueSymbol;
            public int ResidueLocInPeptide;
            public AminoAcidModInfo.ResidueTerminusState TerminusState;

            /// <summary>
            /// Show the residue symbol and mod mass
            /// </summary>
            public override string ToString()
            {
                return ResidueSymbol is default(char)
                    ? string.Format("{0:F4}", ModMass)
                    : string.Format("{0}: {1:F4}", ResidueSymbol, ModMass);
            }
        }

        /// <summary>
        /// Keys in this dictionary are modification names, values are modification masses
        /// </summary>
        private readonly Dictionary <string, double> mModificationMassByName;

        private readonly PeptideCleavageStateCalculator mPeptideCleavageStateCalculator;

        /// <summary>
        /// Add modifications to a peptide read from the DIA-NN synopsis file
        /// Next, compute the monoisotopic mass
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <returns>True if success, false if an error</returns>
        private bool AddModificationsAndComputeMass(
            DiaNNResults searchResult,
            bool updateModOccurrenceCounts)
        {
            try
            {
                // Some of the other tools add IsotopicMods here
                // This is not supported for DIA-NN
                //
                // searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts)

                // Parse .Modifications to determine the modified residues present
                AddModificationsToResidues(searchResult, updateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();

                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }

        /// <summary>
        /// Add modifications to a peptide read from the DIA-NN synopsis file
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        private void AddModificationsToResidues(
            DiaNNResults searchResult,
            bool updateModOccurrenceCounts)
        {
            if (string.IsNullOrWhiteSpace(searchResult.Modifications))
            {
                return;
            }

            // Associate each static and dynamic mod with its corresponding residue
            foreach (var modEntry in GetPeptideModifications(searchResult))
            {
                searchResult.SearchResultAddModification(
                    modEntry.ModMass, modEntry.ResidueSymbol, modEntry.ResidueLocInPeptide,
                    modEntry.TerminusState, updateModOccurrenceCounts);
            }
        }

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        private void AssignRankByScore(
            IList<DiaNNSearchResult> searchResults,
            int startIndex,
            int endIndex)
        {
            if (startIndex == endIndex)
            {
                // Only one result
                var currentResult = searchResults[startIndex];
                currentResult.RankScore = 1;
                searchResults[startIndex] = currentResult;
                return;
            }

            // Duplicate a portion of searchResults so that we can sort by ascending QValue

            var resultsSubset = new Dictionary<int, DiaNNSearchResult>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsByScore = (from item in resultsSubset orderby item.Value.QValue, item.Value.ConfidenceScore descending select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByScore)
            {
                var result = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = result.QValue;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(result.QValue - lastValue) > double.Epsilon)
                    {
                        lastValue = result.QValue;
                        currentRank++;
                    }
                }

                result.RankScore = currentRank;
                searchResults[entry.Key] = result;
            }
        }

        /// <summary>
        /// Compute the monoisotopic MH value using the calculated monoisotopic mass
        /// </summary>
        /// <remarks>This is (M+H)+ when the charge carrier is a proton</remarks>
        /// <param name="searchResult"></param>
        /// <returns>(M+H)+, as a string</returns>
        private string ComputeMH(ToolResultsBaseClass searchResult)
        {
            return PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(searchResult.CalculatedMonoMassValue, 0), 6);
        }

        /// <summary>
        /// Computes the total of all modifications defined for the sequence
        /// </summary>
        /// <param name="searchResult"></param>
        private double ComputeTotalModMass(ToolResultsBaseClass searchResult)
        {
            if (string.IsNullOrWhiteSpace(searchResult.ModificationList))
            {
                return 0;
            }

            return GetPeptideModifications(searchResult).Sum(modEntry => modEntry.ModMass);
        }

        /// <summary>
        /// This routine creates a synopsis file from the output from DIA-NN (file report.tsv)
        /// The synopsis file includes every result with a probability above a set threshold
        /// </summary>
        /// <param name="inputFilePath">DIA-NN results file</param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="baseName">Output: base synopsis file name</param>
        /// <param name="synOutputFilePath">Output: synopsis file path created by this method</param>
        /// <param name="filterPassingResultCount">Output: number of filter passing results</param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputDirectoryPath,
            out string baseName,
            out string synOutputFilePath,
            out int filterPassingResultCount)
        {
            filterPassingResultCount = 0;

            try
            {
                var errorMessages = new List<string>();

                var inputFile = new FileInfo(inputFilePath);
                if (inputFile.Directory == null)
                {
                    SetErrorMessage("Unable to determine the parent directory of file " + inputFile.FullName);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);

                    baseName = string.Empty;
                    synOutputFilePath = string.Empty;

                    return false;
                }

                var success = ReadDiaNNResults(inputFile, errorMessages, out var filteredSearchResults, out var datasetNameToBaseNameMap);

                if (!success)
                {
                    baseName = string.Empty;
                    synOutputFilePath = string.Empty;

                    return false;
                }

                if (filteredSearchResults.Count == 0)
                {
                    // The report.tsv file does not have any filter passing data

                    baseName = string.Empty;
                    synOutputFilePath = string.Empty;

                    return false;
                }

                Console.WriteLine();

                var datasetNames = new SortedSet<string>();

                foreach (var datasetName in datasetNameToBaseNameMap.Keys)
                {
                    datasetNames.Add(datasetName);
                }

                GetDatasetNameMap(inputFile.Name, datasetNames, out var longestCommonBaseName);

                // The synopsis file name will be of the form DatasetName_diann_syn.txt
                // If datasetNameToBaseNameMap only has one item, will use the full dataset name
                // If datasetNameToBaseNameMap has multiple items, will use either Options.OutputFileBaseName,
                // or the longest string in common for the keys in datasetNameToBaseNameMap

                baseName = GetBaseNameForOutputFiles(datasetNameToBaseNameMap, "diann", longestCommonBaseName);

                synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);

                Console.WriteLine();
                OnStatusEvent("Creating synopsis file, " + PathUtils.CompactPathString(synOutputFilePath, 80));

                using var writer = new StreamWriter(new FileStream(synOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                // Write the header line to the output file
                WriteSynFHTFileHeader(writer, errorMessages);

                // Write the search results to disk
                WriteFilteredSearchResults(datasetNameToBaseNameMap, writer, filteredSearchResults, errorMessages);

                filterPassingResultCount = filteredSearchResults.Count;

                // Inform the user if any errors occurred
                if (errorMessages.Count > 0)
                {
                    SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateSynResultsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);

                baseName = string.Empty;
                synOutputFilePath = string.Empty;

                return false;
            }
        }

        private double GetModificationMass(string modificationName)
        {
            if (mModificationMassByName.TryGetValue(modificationName, out var modMass))
                return modMass;

            return 0;
        }

        private List<DiaNNModInfo> GetPeptideModifications(DiaNNResults searchResult)
        {
            return GetPeptideModifications(searchResult.PeptideCleanSequence, searchResult.Modifications);
        }

        private List<DiaNNModInfo> GetPeptideModifications(ToolResultsBaseClass searchResult)
        {
            return GetPeptideModifications(searchResult.Sequence, searchResult.ModificationList);
        }

        private List<DiaNNModInfo> GetPeptideModifications(string cleanSequence, string peptideWithModifications)
        {
            // ReSharper disable CommentTypo

            // peptideWithModifications should have the peptide sequence, including modification names; examples:
            //   AAAGDLGGDHLAFSC(UniMod:4)DVAK
            //   AAGVGLVDC(UniMod:4)HC(UniMod:4)HLSAPDFDR

            // ReSharper enable CommentTypo

            var mods = new List<DiaNNModInfo>();

            var finalResidueLoc = cleanSequence.Length;

            var currentIndex = -1;
            var maxIndex = peptideWithModifications.Length - 1;

            var residueNumber = 0;
            var currentResidue = '-';

            while (currentIndex < maxIndex)
            {
                currentIndex++;

                if (!peptideWithModifications[currentIndex + 1].Equals('('))
                {
                    if (char.IsLetter(peptideWithModifications[currentIndex]))
                    {
                        currentResidue = peptideWithModifications[currentIndex];
                        residueNumber++;
                    }

                    continue;
                }

                var closingParenthesisIndex = peptideWithModifications.IndexOf(')', currentIndex + 1);

                if (closingParenthesisIndex < 0)
                {
                    OnWarningEvent("Mismatched parenthesis after index {0} in {1}", currentIndex, peptideWithModifications);
                    return mods;
                }

                if (closingParenthesisIndex == currentIndex + 1)
                {
                    OnWarningEvent("Empty modification name after index {0} in {1}", currentIndex, peptideWithModifications);
                    continue;
                }

                // Parse out the modification name
                var modificationName = peptideWithModifications.Substring(currentIndex + 2, closingParenthesisIndex - currentIndex - 2);

                var modMass = GetModificationMass(modificationName);

                var currentMod = new DiaNNModInfo
                {
                    ResidueSymbol = currentResidue,
                    ResidueLocInPeptide = residueNumber,
                    ModMass = modMass
                };

                if (currentMod.ResidueLocInPeptide <= 1)
                {
                    currentMod.TerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                }
                else if (currentMod.ResidueLocInPeptide >= finalResidueLoc)
                {
                    currentMod.TerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                }
                else
                {
                    currentMod.TerminusState = AminoAcidModInfo.ResidueTerminusState.None;
                }

                // Assure that ResidueLocInPeptide is between 1 and finalResidueLoc
                if (currentMod.ResidueLocInPeptide < 1)
                {
                    currentMod.ResidueLocInPeptide = 1;
                }
                else if (currentMod.ResidueLocInPeptide > finalResidueLoc)
                {
                    currentMod.ResidueLocInPeptide = finalResidueLoc;
                }

                mods.Add(currentMod);

                currentIndex = closingParenthesisIndex;
            }

            return mods;
        }

        /// <summary>
        /// Read the static and dynamic modifications from the DIA-NN parameter file
        /// </summary>
        /// <param name="diannParamFilePath"></param>
        /// <returns>True on success, false if an error</returns>
        private bool LoadSearchEngineParamFile(string diannParamFilePath)
        {
            try
            {
                var sourceFile = new FileInfo(diannParamFilePath);
                if (!sourceFile.Exists)
                {
                    SetErrorMessage("DIA-NN parameter file not found: " + diannParamFilePath);
                    SetErrorCode(PHRPErrorCode.ParameterFileNotFound);
                    return false;
                }

                OnStatusEvent("Reading the DIA-NN parameter file: " + PathUtils.CompactPathString(sourceFile.FullName, 110));

                var startupOptions = new StartupOptions
                {
                    DisableOpeningInputFiles = true,
                    LoadModsAndSeqInfo = false
                };

                var reader = new DiaNNSynFileReader("DiaNN_ParamFile_Reader", sourceFile.FullName, startupOptions);
                RegisterEvents(reader);

                // ReSharper disable once UnusedVariable
                var success = reader.LoadSearchEngineParameters(sourceFile.FullName, out var searchEngineParams);

                if (!success)
                {
                    return false;
                }

                return ExtractModInfoFromParamFile(diannParamFilePath);
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadSearchEngineParamFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the DIA-NN parameter file
        /// </summary>
        /// <remarks>The DMS-based parameter file for DIA-NN uses the same formatting as MS-GF+</remarks>
        /// <param name="diannParamFilePath"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromParamFile(
            string diannParamFilePath)
        {
            var modFileProcessor = new MSGFPlusParamFileModExtractor(TOOL_NAME);
            RegisterEvents(modFileProcessor);

            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                diannParamFilePath,
                MSGFPlusParamFileModExtractor.ModSpecFormats.DiaNN,
                out var modList);

            if (!success || mErrorCode != PHRPErrorCode.NoError)
            {
                if (mErrorCode == PHRPErrorCode.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the DIA-NN parameter file");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            // ToDo: Verify the behavior of this call
            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modList, mPeptideMods);

            mModificationMassByName.Clear();

            // Cache the modification names and masses
            foreach (var item in modList)
            {
                if (mModificationMassByName.TryGetValue(item.ModName, out var existingModMass))
                {
                    if (Math.Abs(existingModMass - item.ModMassVal) > 0.01)
                    {
                        OnWarningEvent(
                            "The DIA-NN parameter file has two static and/or dynamic mods with the same name ({0}) but different masses ({1} and {2})",
                            item.ModName, existingModMass, item.ModMassVal);

                        return false;
                    }

                    continue;
                }

                mModificationMassByName.Add(item.ModName, item.ModMassVal);
            }

            return true;
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            var columnMapping = new Dictionary<DiaNNSynFileColumns, int>();

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (resetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                // Initialize searchResult
                var searchResult = new DiaNNResults(mPeptideMods, mPeptideSeqMassCalculator);

                var errorMessages = new List<string>();

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(Options);

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                    var resultsProcessed = 0;
                    var headerParsed = false;

                    // Create the output files
                    var baseOutputFilePath = Path.Combine(outputDirectoryPath, Path.GetFileName(inputFilePath));

                    var filesInitialized = InitializeSequenceOutputFiles(baseOutputFilePath);
                    if (!filesInitialized)
                        return false;

                    // Parse the input file
                    while (!reader.EndOfStream && !AbortProcessing)
                    {
                        var lineIn = reader.ReadLine();

                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        if (!headerParsed)
                        {
                            var validHeader = ParseDiaNNSynFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }
                            headerParsed = true;
                            continue;
                        }

                        var validSearchResult = ParseDiaNNSynFileEntry(
                            lineIn, searchResult, errorMessages,
                            resultsProcessed, columnMapping);

                        resultsProcessed++;
                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, true);
                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
                            errorMessages.Add(string.Format("Error adding modifications to sequence for ResultID '{0}'", searchResult.ResultID));
                        }

                        SaveResultsFileEntrySeqInfo(searchResult, true);

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    if (Options.CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(inputFilePath);
                        var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));

                        if (string.IsNullOrWhiteSpace(modificationSummaryFilePath))
                        {
                            OnWarningEvent("ParseDiaNNSynopsisFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
                        }
                        else
                        {
                            modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                            SaveModificationSummaryFile(modificationSummaryFilePath);
                        }
                    }

                    // Inform the user if any errors occurred
                    if (errorMessages.Count > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                    }

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in ParseDiaNNSynopsisFile", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                    return false;
                }
                finally
                {
                    CloseSequenceOutputFiles();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating the output file in ParseDiaNNSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parse a DIA-NN results line while creating the DIA-NN synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="columnMapping"></param>
        /// <param name="lineNumber">Line number in the input file (used for error reporting)</param>
        /// <param name="datasetNameToBaseNameMap">Keys are full dataset names, values are abbreviated dataset names</param>
        /// <param name="currentDatasetFile">Current dataset file path; updated by this method if the Spectrum File path is not an empty string</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNResultsFileEntry(
            string lineIn,
            out DiaNNSearchResult searchResult,
            ICollection<string> errorMessages,
            IDictionary<DiaNNReportFileColumns, int> columnMapping,
            int lineNumber,
            IDictionary<string, string> datasetNameToBaseNameMap,
            ref string currentDatasetFile)
        {
            searchResult = new DiaNNSearchResult();

            try
            {
                searchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                // The file should have over 20 columns, but we'll only require 15
                if (splitLine.Length < 15)
                {
                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.DatasetFile], out string spectrumFile);

                bool updateDatasetNameMapDictionary;
                if (!string.IsNullOrWhiteSpace(spectrumFile))
                {
                    currentDatasetFile = spectrumFile;
                    searchResult.DatasetFile = spectrumFile;
                    updateDatasetNameMapDictionary = true;
                }
                else
                {
                    searchResult.DatasetFile = currentDatasetFile;
                    updateDatasetNameMapDictionary = false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.DatasetName], out searchResult.DatasetName);

                if (updateDatasetNameMapDictionary)
                {
                    var fullDatasetName = Path.GetFileNameWithoutExtension(spectrumFile);

                    if (!baseNameByDatasetName.ContainsKey(fullDatasetName) && !string.IsNullOrWhiteSpace(searchResult.DatasetName))
                    {
                        baseNameByDatasetName.Add(fullDatasetName, searchResult.DatasetName);
                    }
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinGroup], out searchResult.ProteinGroup);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinIDs], out searchResult.ProteinIDs);

                var proteinIdList = searchResult.ProteinIDs.Split(';').ToList();

                if (proteinIdList.Count > 0)
                    searchResult.Protein = proteinIdList[0];

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinNames], out searchResult.ProteinNames);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GeneNames], out searchResult.GeneNames);

                // Cannot compute NTT since peptides in report.tsv do not have prefix and suffix residues
                // DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NTT], out searchResult.NumberOfTrypticTermini);


                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinGroupQuantity], out searchResult.ProteinGroupQuantity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinGroupNormalized], out searchResult.ProteinGroupNormalized);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinGroupMaxLFQ], out searchResult.ProteinGroupMaxLFQ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GenesQuantity], out searchResult.GenesQuantity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GenesNormalized], out searchResult.GenesNormalized);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GenesMaxLFQ], out searchResult.GenesMaxLFQ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GenesMaxLFQUnique], out searchResult.GenesMaxLFQUnique);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ModifiedSequence], out searchResult.ModifiedSequence);

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.StrippedSequence], out searchResult.Sequence))
                {
                    ReportError("Stripped.Sequence column is missing or invalid on line " + lineNumber, true);
                }

                searchResult.Length = searchResult.Sequence.Length;

                // Peptide sequences in report.tsv do not include prefix or suffix residues
                // DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrevAA], out searchResult.PrefixResidue);
                // DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NextAA], out searchResult.SuffixResidue);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorId], out searchResult.PrecursorId);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorCharge], out searchResult.Charge);
                searchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(searchResult.Charge, 0));

                if (DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.QValue], out searchResult.QValueDiaNN))
                {
                    double.TryParse(searchResult.QValueDiaNN, out searchResult.QValue);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PEP], out searchResult.PEP);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GlobalQValue], out searchResult.GlobalQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinQValue], out searchResult.ProteinQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinGroupQValue], out searchResult.ProteinGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GlobalProteinGroupQValue], out searchResult.GlobalProteinGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.GeneGroupQValue], out searchResult.GeneGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.TranslatedQValue], out searchResult.TranslatedQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Proteotypic], out searchResult.Proteotypic);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorQuantity], out searchResult.PrecursorQuantity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorNormalized], out searchResult.PrecursorNormalized);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorTranslated], out searchResult.PrecursorTranslated);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.TranslatedQuality], out searchResult.TranslatedQuality);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MS1Translated], out searchResult.MS1Translated);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.QuantityQuality], out searchResult.QuantityQuality);

                // Retention time is in minutes
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.RT], out searchResult.ElutionTime);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.RTStart], out searchResult.RTStart);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.RTStop], out searchResult.RTStop);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.IndexedRT], out searchResult.IndexedRT);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PredictedRT], out searchResult.PredictedRT);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PredictedIndexedRT], out searchResult.PredictedIndexedRT);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.FirstProteinDescription], out searchResult.FirstProteinDescription);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.LibQValue], out searchResult.LibQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.LibProteinGroupQValue], out searchResult.LibProteinGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MS1ProfileCorr], out searchResult.MS1ProfileCorr);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MS1Area], out searchResult.MS1Area);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Evidence], out searchResult.Evidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.SpectrumSimilarity], out searchResult.SpectrumSimilarity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Averagine], out searchResult.Averagine);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MassEvidence], out searchResult.MassEvidence);

                if (DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.CScore], out searchResult.CScore))
                {
                    double.TryParse(searchResult.CScore, out searchResult.ConfidenceScore);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.DecoyEvidence], out searchResult.DecoyEvidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.DecoyCScore], out searchResult.DecoyCScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.FragmentQuantRaw], out searchResult.FragmentQuantRaw);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.FragmentQuantCorrected], out searchResult.FragmentQuantCorrected);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.FragmentCorrelations], out searchResult.FragmentCorrelations);

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MS2Scan], out searchResult.Scan))
                {
                    ReportError("MS2.Scan column is missing or invalid on line " + lineNumber, true);
                }

                if (!int.TryParse(searchResult.Scan, out searchResult.ScanNum))
                {
                    ReportError("MS2.Scan column is not numeric on line " + lineNumber, true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.IonMobility], out searchResult.IonMobility);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.IndexedIonMobility], out searchResult.IndexedIonMobility);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PredictedIonMobility], out searchResult.PredictedIonMobility);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PredictedIndexedIonMobility], out searchResult.PredictedIndexedIonMobility);

                searchResult.Reverse = IsReversedProtein(searchResult.Protein);

                // Parse the modification list to determine the total mod mass
                var totalModMass = ComputeTotalModMass(searchResult);

                // Compute monoisotopic mass of the peptide
                searchResult.CalculatedMonoMassPHRP = ComputePeptideMassForCleanSequence(searchResult.Sequence, totalModMass);

                searchResult.CalculatedMonoMassValue = searchResult.CalculatedMonoMassPHRP;


                // ToDo: compute the peptide mass,then compute PrecursorMZ, MH, and Mass

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                searchResult.MH = ComputeMH(searchResult);


                // Assume the peptide is preceded by a K or R, then compute cleavage state and missed cleavages

                var peptideWithPrefixAndSuffix = "K." + searchResult.Sequence + ".-";

                // Use the peptide cleavage state calculator to compute the missed cleavage count

                var missedCleavageCount = mPeptideCleavageStateCalculator.ComputeNumberOfMissedCleavages(peptideWithPrefixAndSuffix);

                searchResult.MissedCleavageCount = missedCleavageCount.ToString();

                var cleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(peptideWithPrefixAndSuffix);

                searchResult.NumberOfTrypticTermini = cleavageState switch
                {
                    PeptideCleavageStateCalculator.PeptideCleavageState.Full => 2,
                    PeptideCleavageStateCalculator.PeptideCleavageState.Partial => 1,
                    PeptideCleavageStateCalculator.PeptideCleavageState.NonSpecific => 0,
                    PeptideCleavageStateCalculator.PeptideCleavageState.Unknown => 0,
                    _ => 0
                };

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the DIA-NN results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing DiaNN results in ParseDiaNNResultsFileEntry, line {0}: {1}", lineNumber, ex.Message));
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the DIA-NN report.tsv file header line, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseDiaNNResultsFileHeaderLine(
            string lineIn,
            IDictionary<DiaNNReportFileColumns, int> columnMapping)
        {
            // ReSharper disable StringLiteralTypo

            // Columns in report.tsv files
            var columnNames = new SortedDictionary<string, DiaNNReportFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "File.Name", DiaNNReportFileColumns.DatasetFile },
                { "Run", DiaNNReportFileColumns.DatasetName },
                { "Protein.Group", DiaNNReportFileColumns.ProteinGroup },
                { "Protein.Ids", DiaNNReportFileColumns.ProteinIDs },
                { "Protein.Names", DiaNNReportFileColumns.ProteinNames },
                { "Genes", DiaNNReportFileColumns.GeneNames },
                { "PG.Quantity", DiaNNReportFileColumns.ProteinGroupQuantity },
                { "PG.Normalised", DiaNNReportFileColumns.ProteinGroupNormalized },
                { "PG.MaxLFQ", DiaNNReportFileColumns.ProteinGroupMaxLFQ },
                { "Genes.Quantity", DiaNNReportFileColumns.GenesQuantity },
                { "Genes.Normalised", DiaNNReportFileColumns.GenesNormalized },
                { "Genes.MaxLFQ", DiaNNReportFileColumns.GenesMaxLFQ },
                { "Genes.MaxLFQ.Unique", DiaNNReportFileColumns.GenesMaxLFQUnique },
                { "Modified.Sequence", DiaNNReportFileColumns.ModifiedSequence },
                { "Stripped.Sequence", DiaNNReportFileColumns.StrippedSequence },
                { "Precursor.Id", DiaNNReportFileColumns.PrecursorId },
                { "Precursor.Charge", DiaNNReportFileColumns.PrecursorCharge },
                { "Q.Value", DiaNNReportFileColumns.QValue },
                { "PEP", DiaNNReportFileColumns.PEP },
                { "Global.Q.Value", DiaNNReportFileColumns.GlobalQValue },
                { "Protein.Q.Value", DiaNNReportFileColumns.ProteinQValue },
                { "PG.Q.Value", DiaNNReportFileColumns.ProteinGroupQValue },
                { "Global.PG.Q.Value", DiaNNReportFileColumns.GlobalProteinGroupQValue },
                { "GG.Q.Value", DiaNNReportFileColumns.GeneGroupQValue },
                { "Translated.Q.Value", DiaNNReportFileColumns.TranslatedQValue },
                { "Proteotypic", DiaNNReportFileColumns.Proteotypic },
                { "Precursor.Quantity", DiaNNReportFileColumns.PrecursorQuantity },
                { "Precursor.Normalised", DiaNNReportFileColumns.PrecursorNormalized },
                { "Precursor.Translated", DiaNNReportFileColumns.PrecursorTranslated },
                { "Translated.Quality", DiaNNReportFileColumns.TranslatedQuality },
                { "Ms1.Translated", DiaNNReportFileColumns.MS1Translated },
                { "Quantity.Quality", DiaNNReportFileColumns.QuantityQuality },
                { "RT", DiaNNReportFileColumns.RT },
                { "RT.Start", DiaNNReportFileColumns.RTStart },
                { "RT.Stop", DiaNNReportFileColumns.RTStop },
                { "iRT", DiaNNReportFileColumns.IndexedRT },
                { "Predicted.RT", DiaNNReportFileColumns.PredictedRT },
                { "Predicted.iRT", DiaNNReportFileColumns.PredictedIndexedRT },
                { "First.Protein.Description", DiaNNReportFileColumns.FirstProteinDescription },
                { "Lib.Q.Value", DiaNNReportFileColumns.LibQValue },
                { "Lib.PG.Q.Value", DiaNNReportFileColumns.LibProteinGroupQValue },
                { "Ms1.Profile.Corr", DiaNNReportFileColumns.MS1ProfileCorr },
                { "Ms1.Area", DiaNNReportFileColumns.MS1Area },
                { "Evidence", DiaNNReportFileColumns.Evidence },
                { "Spectrum.Similarity", DiaNNReportFileColumns.SpectrumSimilarity },
                { "Averagine", DiaNNReportFileColumns.Averagine },
                { "Mass.Evidence", DiaNNReportFileColumns.MassEvidence },
                { "CScore", DiaNNReportFileColumns.CScore },
                { "Decoy.Evidence", DiaNNReportFileColumns.DecoyEvidence },
                { "Decoy.CScore", DiaNNReportFileColumns.DecoyCScore },
                { "Fragment.Quant.Raw", DiaNNReportFileColumns.FragmentQuantRaw },
                { "Fragment.Quant.Corrected", DiaNNReportFileColumns.FragmentQuantCorrected },
                { "Fragment.Correlations", DiaNNReportFileColumns.FragmentCorrelations },
                { "MS2.Scan", DiaNNReportFileColumns.MS2Scan },
                { "IM", DiaNNReportFileColumns.IonMobility },
                { "iIM", DiaNNReportFileColumns.IndexedIonMobility },
                { "Predicted.IM", DiaNNReportFileColumns.PredictedIonMobility },
                { "Predicted.iIM", DiaNNReportFileColumns.PredictedIndexedIonMobility }
            };

            // ReSharper restore StringLiteralTypo

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (DiaNNReportFileColumns resultColumn in Enum.GetValues(typeof(DiaNNReportFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');

                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the DiaNN results file", ex);
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a DIA-NN _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNSynFileHeaderLine(string lineIn, IDictionary<DiaNNSynFileColumns, int> columnMapping)
        {
            var columnNames = DiaNNSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (DiaNNSynFileColumns resultColumn in Enum.GetValues(typeof(DiaNNSynFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');

                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the DiaNN synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a DIA-NN Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNSynFileEntry(
            string lineIn,
            DiaNNResults searchResult,
            ICollection<string> errorMessages,
            int resultsProcessed,
            IDictionary<DiaNNSynFileColumns, int> columnMapping)
        {
            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 15)
                {
                    return false;
                }

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ResultID], out string value))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading ResultID value from DiaNN results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Dataset], out string dataset);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DatasetID], out int datasetId);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Scan], out string scan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Charge], out string charge);

                searchResult.DatasetName = dataset;
                searchResult.DatasetID = datasetId;

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Peptide], out string peptideSequence))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading Peptide sequence value from DiaNN results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Protein], out string proteinName);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.AdditionalProteins], out string additionalProteins);

                if (!string.IsNullOrWhiteSpace(proteinName))
                {
                    var trimmedName = proteinName.Trim();

                    searchResult.ProteinName = trimmedName;
                    searchResult.Proteins.Add(trimmedName);
                }

                if (!string.IsNullOrWhiteSpace(additionalProteins))
                {
                    // Protein names in the AdditionalProteins column are a comma separated list
                    foreach (var additionalProtein in additionalProteins.Split(','))
                    {
                        var trimmedName = additionalProtein.Trim();

                        if (string.IsNullOrWhiteSpace(trimmedName))
                            continue;

                        searchResult.Proteins.Add(trimmedName);
                    }

                    if (string.IsNullOrWhiteSpace(searchResult.ProteinName) && searchResult.Proteins.Count > 0)
                    {
                        searchResult.ProteinName = searchResult.Proteins[0];
                    }
                }

                if (searchResult.Proteins.Count > 0)
                {
                    searchResult.MultipleProteinCount = (searchResult.Proteins.Count - 1).ToString();
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PrecursorMZ], out string precursorMz);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MH], out string parentIonMH);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Mass], out string monoisotopicMass);

                searchResult.PrecursorMZ = precursorMz;
                searchResult.ParentIonMH = parentIonMH;
                searchResult.CalculatedMonoMass = monoisotopicMass;

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DelM], out string phrpComputedDelM);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DelM_PPM], out string phrpComputedDelMppm);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DelM_DiaNN], out string DiaNNComputedDelM);

                searchResult.PeptideDeltaMass = phrpComputedDelM;
                searchResult.PHRPComputedDelM = phrpComputedDelM;
                searchResult.PHRPComputedDelMPPM = phrpComputedDelMppm;

                searchResult.DiaNNComputedDelM = DiaNNComputedDelM;

                // Note that DIA-NN peptides don't actually have mod symbols; that information is tracked via searchResult.Modifications
                // Thus, .PeptideSequenceWithMods will not have any mod symbols

                // Calling this method will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the cleavage state and terminus state
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Modifications], out string modifications);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.NTT], out string ntt);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.EValue], out string eValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Hyperscore], out string hyperscore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Nextscore], out string nextScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PeptideProphetProbability], out string peptideProphetProbability);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ElutionTime], out string elutionTime);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ElutionTimeAverage], out string elutionTimeAverage);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MissedCleavages], out string missedCleavages);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.NumberOfMatchedIons], out string numberOfMatchedIons);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.TotalNumberOfIons], out string totalNumberOfIons);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.QValue], out string qValue);

                // Store the data
                searchResult.Modifications = modifications;

                searchResult.NTT = ntt;
                searchResult.EValue = eValue;
                searchResult.Hyperscore = hyperscore;
                searchResult.Nextscore = nextScore;
                searchResult.PeptideProphetProbability = peptideProphetProbability;
                searchResult.ElutionTime = elutionTime;
                searchResult.ElutionTimeAverage = elutionTimeAverage;
                searchResult.MissedCleavageCount = missedCleavages;
                searchResult.NumberOfMatchedIons = numberOfMatchedIons;
                searchResult.TotalNumberOfIons = totalNumberOfIons;
                searchResult.QValue = qValue;

                if (string.IsNullOrWhiteSpace(searchResult.MultipleProteinCount))
                {
                    searchResult.MultipleProteinCount = "0";
                }

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing DiaNN results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing DiaNN Results in ParseDiaNNSynFileEntry: " + ex.Message);
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">DIA-NN results file (report.tsv)</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if successful, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            var success = false;

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                    return false;
                }

                success = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!success)
                {
                    return false;
                }

                ResetProgress("Parsing " + PathUtils.CompactPathString(inputFilePath, 110));

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath, Options.AlternateBasePath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);
                    if (inputFile.Directory == null)
                    {
                        SetErrorMessage("Unable to determine the parent directory of the input file: " + inputFilePath);
                        SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                        return false;
                    }

                    if (string.IsNullOrWhiteSpace(Options.SearchToolParameterFilePath))
                    {
                        OnWarningEvent("DIA-NN parameter file not defined; unable to determine modification masses");
                    }
                    else
                    {
                        var diannParameterFilePath = ResolveFilePath(inputFile.DirectoryName, Options.SearchToolParameterFilePath);

                        // Examine the DIA-NN parameter file to determine dynamic and static mods
                        var paramFileLoaded = LoadSearchEngineParamFile(diannParameterFilePath);

                        if (!paramFileLoaded)
                            return false;
                    }

                    // Do not create a first-hits file for DIA-NN results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    success = CreateSynResultsFile(
                        inputFilePath, outputDirectoryPath,
                        out var baseName,
                        out var synOutputFilePath,
                        out var filterPassingResultCount);

                    if (!success)
                    {
                        return false;
                    }

                    if (filterPassingResultCount == 0)
                    {
                        OnWarningEvent("Aborting processing since no filter passing results were found");
                        return true;
                    }

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to create the other PHRP files
                    success = ParseDiaNNSynopsisFile(synOutputFilePath, outputDirectoryPath, false);
                    if (!success)
                    {
                        return false;
                    }

                    if (Options.CreateProteinModsFile)
                    {
                        // Check for an empty synopsis file
                        if (!ValidateFileHasData(synOutputFilePath, "Synopsis file", out var errorMessage))
                        {
                            OnWarningEvent(errorMessage);
                        }
                        else
                        {
                            // ToDo: validate this statement, or remove
                            // ToDo: Use a higher match error threshold since some peptides reported by DIA-NN don't perfectly match the FASTA file

                            const int MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD = 15;
                            const int MATCH_ERROR_PERCENT_WARNING_THRESHOLD = 5;

                            success = CreateProteinModsFileWork(
                                baseName, inputFile,
                                synOutputFilePath, outputDirectoryPath,
                                PeptideHitResultTypes.DiaNN,
                                MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD,
                                MATCH_ERROR_PERCENT_WARNING_THRESHOLD);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in DiaNNResultsProcessor.ProcessFile (2)", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in DiaNNResultsProcessor.ProcessFile (1)", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        /// <summary>
        /// Load DIA-NN search results from a report.tsv file
        /// </summary>
        /// <remarks>If the file is found, but has no results, this method still returns true</remarks>
        /// <param name="inputFile">report.tsv file</param>
        /// <param name="errorMessages"></param>
        /// <param name="filteredSearchResults">Output: DIA-NN results</param>
        /// <param name="datasetNameToBaseNameMap">Keys are full dataset names, values are abbreviated dataset names</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ReadDiaNNResults(
            FileSystemInfo inputFile,
            ICollection<string> errorMessages,
            out List<DiaNNSearchResult> filteredSearchResults,
            out Dictionary<string, string> datasetNameToBaseNameMap)
        {
            filteredSearchResults = new List<DiaNNSearchResult>();
            datasetNameToBaseNameMap = new Dictionary<string, string>();

            try
            {
                if (!inputFile.Exists)
                {
                    OnErrorEvent("File not found: {0}", inputFile.FullName);
                    Console.WriteLine();
                    return false;
                }

                OnStatusEvent("Reading DiaNN results file, " + PathUtils.CompactPathString(inputFile.FullName, 80));

                var columnMapping = new Dictionary<DiaNNReportFileColumns, int>();

                // Open the input file and parse it
                using var reader = new StreamReader(new FileStream(inputFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var headerParsed = false;
                var lineNumber = 0;

                // Initialize the list that will hold all of the records in the DIA-NN result file
                var searchResultsUnfiltered = new List<DiaNNSearchResult>();

                var currentDatasetFile = string.Empty;

                // Parse the input file
                while (!reader.EndOfStream && !AbortProcessing)
                {
                    var lineIn = reader.ReadLine();
                    lineNumber++;

                    if (string.IsNullOrWhiteSpace(lineIn))
                    {
                        continue;
                    }

                    if (!headerParsed)
                    {
                        // Parse the header line
                        var success = ParseDiaNNResultsFileHeaderLine(lineIn, columnMapping);

                        if (!success)
                        {
                            if (string.IsNullOrEmpty(mErrorMessage))
                            {
                                SetErrorMessage("Invalid header line in " + inputFile.Name);
                            }

                            return false;
                        }

                        headerParsed = true;
                        continue;
                    }

                    var validSearchResult = ParseDiaNNResultsFileEntry(
                        lineIn,
                        out var searchResult,
                        errorMessages,
                        columnMapping,
                        lineNumber,
                        datasetNameToBaseNameMap,
                        ref currentDatasetFile);

                    if (validSearchResult)
                    {
                        searchResultsUnfiltered.Add(searchResult);
                    }

                    // Update the progress
                    UpdateSynopsisFileCreationProgress(reader);
                }

                if (searchResultsUnfiltered.Count == 0)
                {
                    OnWarningEvent("No results were found in file {0}", PathUtils.CompactPathString(inputFile.FullName, 110));

                    Console.WriteLine();
                    return true;
                }

                // Sort the SearchResults by dataset name, scan, charge, and ascending QValue
                searchResultsUnfiltered.Sort(new DiaNNSearchResultsComparerDatasetScanChargeQValue());

                // Now filter the data
                var startIndex = 0;

                while (startIndex < searchResultsUnfiltered.Count)
                {
                    // Find all of the matches for the current result's scan
                    // (we sorted by dataset, then scan, so adjacent results will be from the same dataset, except when a new dataset is encountered)
                    // DIA-NN will typically report just one match

                    var endIndex = startIndex;
                    while (endIndex + 1 < searchResultsUnfiltered.Count &&
                           searchResultsUnfiltered[endIndex + 1].ScanNum == searchResultsUnfiltered[startIndex].ScanNum)
                    {
                        endIndex++;
                    }

                    // Store the results for this scan
                    StoreSynMatches(searchResultsUnfiltered, startIndex, endIndex, filteredSearchResults);

                    startIndex = endIndex + 1;
                }

                if (filteredSearchResults.Count == 0)
                {
                    OnWarningEvent("No filter passing results were found in file {0}", PathUtils.CompactPathString(inputFile.FullName, 110));
                }

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in DiaNNResultsProcessor.ReadDiaNNResults", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
                return false;
            }
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan in a single dataset (typically there will only be one result for DIA-NN)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Output parameter: the actual filtered search results</param>
        private void StoreSynMatches(
            IList<DiaNNSearchResult> searchResults,
            int startIndex,
            int endIndex,
            List<DiaNNSearchResult> filteredSearchResults)
        {
            AssignRankByScore(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by dataset name, scan, charge, and Score; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            // Now store the matches that pass the filters
            //  Either QValue < 0.10
            //  or     ConfidenceScore > 0.25
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].QValue < Options.DiaNNQValueThreshold ||
                    searchResults[index].ConfidenceScore > Options.DiaNNConfidenceScoreThreshold)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorMessages"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ICollection<string> errorMessages)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = DiaNNSynFileReader.GetColumnHeaderNamesAndIDs();

                var headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();

                writer.WriteLine(StringUtilities.CollapseList(headerNames));
            }
            catch (Exception)
            {
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add("Error writing synopsis / first hits header");
                }
            }
        }

        /// <summary>
        /// Write search results to disk
        /// </summary>
        /// <param name="datasetNameToBaseNameMap">Keys are full dataset names, values are abbreviated dataset names</param>
        /// <param name="writer"></param>
        /// <param name="filteredSearchResults"></param>
        /// <param name="errorMessages"></param>
        private void WriteFilteredSearchResults(
            Dictionary<string, string> datasetNameToBaseNameMap,
            TextWriter writer,
            List<DiaNNSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            // Lookup the Dataset ID for each dataset (only if on the pnl.gov domain)
            var datasetIDs = LookupDatasetIDs(datasetNameToBaseNameMap.Keys.ToList());

            var index = 1;
            foreach (var result in filteredSearchResults)
            {
                GetBaseNameAndDatasetID(datasetNameToBaseNameMap, datasetIDs, result.DatasetName, out var baseDatasetName, out var datasetID);

                WriteSearchResultToFile(index, baseDatasetName, datasetID, writer, result, errorMessages);
                index++;
            }
        }

        /// <summary>
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="baseDatasetName"></param>
        /// <param name="datasetID"></param>
        /// <param name="writer"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        private void WriteSearchResultToFile(
            int resultID,
            string baseDatasetName,
            int datasetID,
            TextWriter writer,
            DiaNNSearchResult searchResult,
            ICollection<string> errorMessages)
        {
            try
            {
                var data = new List<string>
                {
                    resultID.ToString(),
                    baseDatasetName,
                    datasetID.ToString(),
                    searchResult.Scan,
                    searchResult.Charge,
                    searchResult.PrecursorMZ,
                    searchResult.MassErrorDa,
                    searchResult.MassErrorPpm,
                    searchResult.MassErrorDaDiaNN,
                    searchResult.MH,
                    searchResult.CalculatedMonoMass,
                    searchResult.Sequence,
                    searchResult.ModificationList,
                    searchResult.Protein,
                    searchResult.AdditionalProteins,
                    searchResult.NumberOfTrypticTermini.ToString(),
                    PRISM.StringUtilities.DblToString(searchResult.EValue,5, 0.00005),
                    searchResult.RankScore.ToString(),
                    searchResult.Hyperscore,
                    searchResult.Nextscore,
                    searchResult.PeptideProphetProbability,
                    searchResult.ElutionTime,
                    searchResult.ElutionTimeAverage,
                    searchResult.MissedCleavageCount,
                    searchResult.NumberOfMatchedIons,
                    searchResult.TotalNumberOfIons,
                    PRISM.StringUtilities.DblToString(searchResult.QValue, 5, 0.00005)
                };

                writer.WriteLine(StringUtilities.CollapseList(data));
            }
            catch (Exception)
            {
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add("Error writing synopsis / first hits record");
                }
            }
        }

        /// <summary>
        /// Override this method to display the name of each class
        /// </summary>
        public override string ToString()
        {
            return TOOL_NAME + " results processor";
        }

        private void ModExtractorErrorHandler(string message, Exception ex)
        {
            SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
        }

        private class DiaNNSearchResultsComparerDatasetScanChargeQValue : IComparer<DiaNNSearchResult>
        {
            public int Compare(DiaNNSearchResult x, DiaNNSearchResult y)
            {
                if (x == null || y == null)
                    return 0;

                // First sort on dataset name
                var nameComparisonResult = string.CompareOrdinal(x.DatasetName, y.DatasetName);
                if (nameComparisonResult != 0)
                {
                    return nameComparisonResult;
                }

                // Same dataset, check scan
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                // Scan is the same, check charge
                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {
                    return -1;
                }

                // Charge is the same; check QValue
                if (x.QValue < y.QValue)
                {
                    return -1;
                }

                if (x.QValue > y.QValue)
                {
                    return 1;
                }

                // QValue is the same; check sequence
                var peptideComparisonResult = string.CompareOrdinal(x.Sequence, y.Sequence);

                // ReSharper disable once ConvertIfStatementToReturnStatement
                if (peptideComparisonResult != 0)
                {
                    return peptideComparisonResult;
                }

                // Peptide is the same, check Protein
                return string.CompareOrdinal(x.Protein, y.Protein);
            }
        }

        private class DiaNNSearchResultsComparerQValueScanChargePeptide : IComparer<DiaNNSearchResult>
        {
            public int Compare(DiaNNSearchResult x, DiaNNSearchResult y)
            {
                if (x == null || y == null)
                    return 0;

                if (x.EValue < y.EValue)
                {
                    return -1;
                }

                if (x.EValue > y.EValue)
                {
                    return 1;
                }

                // QValue is the same; check scan number
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                // Scan is the same, check charge
                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {
                    return -1;
                }

                // Charge is the same; check peptide
                return string.CompareOrdinal(x.Sequence, y.Sequence);
            }
        }
    }
}
