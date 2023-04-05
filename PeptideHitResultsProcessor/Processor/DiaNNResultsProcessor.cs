// This class reads in an DiaNN PSM results file and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.
//
// The modification definition information is determined from the DiaNN parameter file
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
using System.Text.RegularExpressions;
using PeptideHitResultsProcessor.Data;
using PeptideHitResultsProcessor.SearchToolResults;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads in DiaNN report.tsv file and creates
    /// a tab-delimited text file with the data
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads DiaNN results file report.tsv
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
    /// 5) Once the entire .tsv has been read, searchResultsUnfiltered is sorted by scan, charge, and ascending expectation value (e-value)
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
            FileDate = "April 4, 2023";

            mPeptideCleavageStateCalculator = new PeptideCleavageStateCalculator();
        }

        /// <summary>
        /// DiaNN tool name
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
        /// C-terminus symbol used by DiaNN
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
            Genes = 5,

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

            /// <summary>
            /// Peptide sequence with modifications
            /// </summary>
            /// <remarks>
            /// Examples:
            ///   AAALEAM(UniMod:35)K
            ///   AAEAHVDAHYYEQNEQPTGTC(UniMod:4)AAC(UniMod:4)ITGDNR
            /// </remarks>
            ModifiedSequence = 13,

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
            /// Translated Quality (unused by PHRP)
            /// </summary>
            TranslatedQuality = 29,

            /// <summary>
            /// MS1 Translated (unused by PHRP)
            /// </summary>
            Ms1Translated = 30,

            /// <summary>
            /// Quantity Quality (unused by PHRP)
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
            PredictedIRT = 37,

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
            Ms1ProfileCorr = 41,

            /// <summary>
            /// MS1 Area
            /// </summary>
            Ms1Area = 42,

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
            IM = 54,

            /// <summary>
            /// Indexed Ion Mobility
            /// </summary>
            IndexedIM = 55,

            /// <summary>
            /// Predicted Ion Mobility
            /// </summary>
            PredictedIM = 56,

            /// <summary>
            /// Predicted Indexed Ion Mobility
            /// </summary>
            PredictedIndexedIM = 57
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

        private readonly Regex mModListResidueModMatcher = new(@"(?<ResidueNumber>\d+)(?<ResidueSymbol>[A-Z])\((?<ModMass>[0-9.-]+)\)", RegexOptions.Compiled);

        private readonly Regex mModListTerminalModMatcher = new(@"(?<TerminusName>[^ ]+-term)\((?<ModMass>[0-9.-]+)\)", RegexOptions.Compiled);

        private readonly PeptideCleavageStateCalculator mPeptideCleavageStateCalculator;

        /// <summary>
        /// This variable keeps track of the number of PSMs whose computed monoisotopic mass
        /// does not agree with the monoisotopic mass computed from the precursor m/z, within a reasonable tolerance
        /// </summary>
        private int mPrecursorMassErrorWarningCount;

        /// <summary>
        /// Precursor match tolerance read from the DiaNN parameter file
        /// </summary>
        private PrecursorMassTolerance mPrecursorMassTolerance;

        /// <summary>
        /// Add modifications to a peptide read from the DiaNN synopsis file
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
                // This is not supported for DiaNN
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
        /// Add modifications to a peptide read from the DiaNN synopsis file
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

            // Duplicate a portion of searchResults so that we can sort by descending Andromeda Score

            var resultsSubset = new Dictionary<int, DiaNNSearchResult>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsByScore = (from item in resultsSubset orderby item.Value.EValue, item.Value.HyperscoreValue descending select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByScore)
            {
                var result = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = result.EValue;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(result.EValue - lastValue) > double.Epsilon)
                    {
                        lastValue = result.EValue;
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
        /// Compute observed DelM and DelM_PPM values
        /// </summary>
        /// <param name="filteredSearchResults">Search results</param>
        private void ComputeObservedMassErrors(IEnumerable<DiaNNSearchResult> filteredSearchResults)
        {
            mPrecursorMassErrorWarningCount = 0;

            foreach (var searchResult in filteredSearchResults)
            {
                if (!double.TryParse(searchResult.PrecursorMZ, out var precursorMz))
                {
                    OnWarningEvent("Invalid Precursor m/z value for scan {0}: {1}", searchResult.Scan, searchResult.PrecursorMZ);
                    continue;
                }

                var observedPrecursorMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMz, searchResult.ChargeNum, 0);

                var deltaMassDa = observedPrecursorMass - searchResult.CalculatedMonoMassPHRP;
                var warningShown = false;

                double peptideDeltaMassPpm;
                double precursorErrorDa;

                if (Math.Abs(deltaMassDa) < 15)
                {
                    peptideDeltaMassPpm = SearchResultsBaseClass.ComputeDelMCorrectedPPM(
                        mPeptideSeqMassCalculator, deltaMassDa, precursorMz, searchResult.ChargeNum,
                        searchResult.CalculatedMonoMassPHRP, true);

                    // Note that this will be a C13-corrected precursor error; not the absolute precursor error
                    precursorErrorDa = PeptideMassCalculator.PPMToMass(peptideDeltaMassPpm, searchResult.CalculatedMonoMassPHRP);
                }
                else
                {
                    // Delta mass value is unreasonably large; do not try to correct the delta mass
                    peptideDeltaMassPpm = PeptideMassCalculator.MassToPPM(deltaMassDa, searchResult.CalculatedMonoMassPHRP);
                    precursorErrorDa = deltaMassDa;

                    mPrecursorMassErrorWarningCount++;
                    ShowPeriodicWarning(mPrecursorMassErrorWarningCount,
                        10,
                        string.Format(
                            "Peptide mass computed by PHRP differs from the precursor mass by more than 15 Da, indicating an error adding static and/or dynamic mods: {0:F2} Da for {1}, Scan {2}",
                            deltaMassDa,
                            searchResult.Sequence,
                            searchResult.Scan));

                    warningShown = true;
                }

                searchResult.MassErrorPpm = PRISM.StringUtilities.DblToString(peptideDeltaMassPpm, 5, 0.00005);
                searchResult.MassErrorDa = StringUtilities.MassErrorToString(precursorErrorDa);

                if (warningShown)
                    continue;

                // ReSharper disable once ConvertIfStatementToSwitchStatement
                if (mPrecursorMassTolerance.IsPPM &&
                    Math.Abs(peptideDeltaMassPpm) > mPrecursorMassTolerance.ToleranceLeft * 3)
                {
                    // Mass computed by PHRP differs from the precursor m/z by a larger amount than expected
                    mPrecursorMassErrorWarningCount++;
                    ShowPeriodicWarning(mPrecursorMassErrorWarningCount,
                        10,
                        string.Format(
                            "Precursor mass error computed by PHRP is more than {0:F0} ppm, indicating a possible error adding static and/or dynamic mods: {1:F2} ppm for {2}, Scan {3}",
                            mPrecursorMassTolerance.ToleranceLeft * 3,
                            peptideDeltaMassPpm,
                            searchResult.Sequence,
                            searchResult.Scan));
                }
                else if (!mPrecursorMassTolerance.IsPPM &&
                         Math.Abs(precursorErrorDa) > mPrecursorMassTolerance.ToleranceLeft * 1.5)
                {
                    // Mass computed by PHRP differs from the precursor m/z by a larger amount than expected
                    mPrecursorMassErrorWarningCount++;
                    ShowPeriodicWarning(mPrecursorMassErrorWarningCount,
                        10,
                        string.Format(
                            "Precursor mass error computed by PHRP is more than {0:F0} Da, indicating a possible error adding static and/or dynamic mods: {1:F2} Da",
                            mPrecursorMassTolerance.ToleranceLeft * 1.5,
                            precursorErrorDa));
                }
            }
        }

        private void ComputeQValues(List<DiaNNSearchResult> filteredSearchResults)
        {
            // Sort filteredSearchResults by E-value, Scan, and Peptide
            filteredSearchResults.Sort(new DiaNNSearchResultsComparerEValueScanChargePeptide());

            var listForQValue = new List<ToolResultsBaseClass>();
            listForQValue.AddRange(filteredSearchResults);

            // Compute FDR values, then assign QValues
            ComputeQValues(listForQValue);
        }

        private IDictionary<int, double> ComputeScanToElutionTimeMap(Dictionary<int, List<double>> elutionTimesByScanNumber)
        {
            var startScan = elutionTimesByScanNumber.Keys.Min();
            var endScan = elutionTimesByScanNumber.Keys.Max();

            var elutionTimeByScanNumber = new SortedDictionary<int, double>();

            foreach (var item in elutionTimesByScanNumber)
            {
                if (item.Value.Count == 0)
                    continue;

                var averageElutionTime = item.Value.Average();
                elutionTimeByScanNumber.Add(item.Key, averageElutionTime);
            }

            if (elutionTimeByScanNumber.Count < 2)
                return elutionTimeByScanNumber;

            var interpolationTool = new BinarySearchFindNearest();
            interpolationTool.AddData(elutionTimeByScanNumber);

            for (var scan = startScan; scan <= endScan; scan++)
            {
                if (elutionTimeByScanNumber.ContainsKey(scan))
                    continue;

                var elutionTime = interpolationTool.GetYForX(scan);

                elutionTimeByScanNumber.Add(scan, elutionTime);
            }

            // Apply a moving average smooth to the times in elutionTimeByScanNumber
            // Note that the keys in the SortedDictionary are sorted ascending and are in the same order as the Values

            var smoothedTimes = DataUtilities.MovingAverageSmooth(elutionTimeByScanNumber.Values);

            var elutionTimeByScanNumberSmoothed = new Dictionary<int, double>();

            var i = 0;
            foreach (var item in elutionTimeByScanNumber)
            {
                elutionTimeByScanNumberSmoothed.Add(item.Key, smoothedTimes[i]);
                i++;
            }

            return elutionTimeByScanNumberSmoothed;
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
        /// Use scan number, charge state, and clean sequence (no prefix or suffix residues) to create a key for a PSM
        /// </summary>
        /// <param name="additionalResult"></param>
        /// <returns>String in the form Scan-Charge-CleanSequence</returns>
        /// <exception cref="NotImplementedException"></exception>
        private string ConstructKeyForPSM(ToolResultsBaseClass additionalResult)
        {
            var cleanSequence = GetCleanSequence(additionalResult.Sequence);
            return string.Format("{0}-{1}-{2}", additionalResult.Scan, additionalResult.Charge, cleanSequence);
        }

        /// <summary>
        /// This routine creates a synopsis file from the output from DiaNN (file Dataset_psm.tsv or Dataset.tsv)
        /// The synopsis file includes every result with a probability above a set threshold
        /// </summary>
        /// <param name="inputFilePath">DiaNN results file</param>
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

                var success = ReadDiaNNResults(inputFile, errorMessages, out var filteredSearchResults);

                // Sort the data in filteredSearchResults and compute Q Values
                ComputeQValues(filteredSearchResults);

                if (!success)
                {
                    baseName = string.Empty;
                    synOutputFilePath = string.Empty;

                    return false;
                }

                if (filteredSearchResults.Count == 0 && inputFile.Name.EndsWith(PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase))
                {
                    // The _psm.tsv file does not have any filter passing data
                    // Instead, load data from the Dataset.tsv file (or files)

                    if (!LoadNonPsmResults(inputFile, errorMessages, filteredSearchResults))
                    {
                        baseName = string.Empty;
                        synOutputFilePath = string.Empty;

                        return false;
                    }
                }
                else
                {
                    // If we loaded results from a _psm.tsv file, merge in the ion counts from the Dataset.tsv file (if it exists)
                    // If we loaded results from a Dataset.tsv file, merge in peptide prophet values and other scores from the _psm.tsv file (if it exists)

                    if (!MergeRelatedDiaNNResults(inputFile, errorMessages, filteredSearchResults))
                    {
                        baseName = string.Empty;
                        synOutputFilePath = string.Empty;

                        return false;
                    }
                }

                Console.WriteLine();

                // Keys in this dictionary are dataset names, values are abbreviated names
                var baseNameByDatasetName = GetDatasetNameMap(inputFile.Name, filteredSearchResults, out var longestCommonBaseName);

                // Compute MassErrorPpm and MassErrorDa
                ComputeObservedMassErrors(filteredSearchResults);

                // The synopsis file name will be of the form DatasetName_DiaNN_syn.txt
                // If baseNameByDatasetName only has one item, will use the full dataset name
                // If baseNameByDatasetName has multiple items, will use either Options.OutputFileBaseName,
                // or the longest string in common for the keys in baseNameByDatasetName

                baseName = GetBaseNameForOutputFiles(baseNameByDatasetName, "DiaNN", longestCommonBaseName);

                synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);

                Console.WriteLine();
                OnStatusEvent("Creating synopsis file, " + PathUtils.CompactPathString(synOutputFilePath, 80));

                using var writer = new StreamWriter(new FileStream(synOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                // Write the header line to the output file
                WriteSynFHTFileHeader(writer, errorMessages);

                // Write the search results to disk
                WriteFilteredSearchResults(baseNameByDatasetName, writer, filteredSearchResults, errorMessages);

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

        /// <summary>
        /// Find the Dataset.tsv files that correspond to an Aggregation_psm.tsv file
        /// </summary>
        /// <param name="inputFile"></param>
        /// <param name="datasetNames"></param>
        /// <returns>Dictionary with the matching files, along with the dataset name for each file</returns>
        private Dictionary<FileInfo, string> FindAggregationPsmSourceFiles(FileInfo inputFile, IEnumerable<string> datasetNames)
        {
            var sourceFiles = new Dictionary<FileInfo, string>();

            if (inputFile.Directory == null)
            {
                return sourceFiles;
            }

            // ReSharper disable once CommentTypo
            // Look for .tsv files in this directory that have columns "expectscore", "hyperscore", and "nextscore"

            var datasetNamesFound = new SortedSet<string>(StringComparer.OrdinalIgnoreCase);

            foreach (var candidateFile in inputFile.Directory.GetFiles("*.tsv"))
            {
                if (candidateFile.FullName.Equals(inputFile.FullName))
                    continue;

                if (candidateFile.Name.EndsWith("_ion.tsv", StringComparison.OrdinalIgnoreCase) ||
                    candidateFile.Name.EndsWith("_peptide.tsv", StringComparison.OrdinalIgnoreCase) ||
                    candidateFile.Name.EndsWith("_protein.tsv", StringComparison.OrdinalIgnoreCase) ||
                    candidateFile.Name.Equals("Aggregation_psm.tsv", StringComparison.OrdinalIgnoreCase))
                {
                    continue;
                }

                var tsvFileFormat = DetermineResultsFileFormat(candidateFile.FullName);

                if (tsvFileFormat != ResultsFileFormat.DiannTSVFile)
                    continue;

                var datasetName = GetDatasetName(candidateFile.Name);

                sourceFiles.Add(candidateFile, datasetName);

                datasetNamesFound.Add(datasetName);
            }

            foreach (var expectedDataset in datasetNames.Where(expectedDataset => !datasetNamesFound.Contains(expectedDataset)))
            {
                OnWarningEvent("Did not find file {0} in {1}", expectedDataset + ".tsv", inputFile.Directory.FullName);
            }

            return sourceFiles;
        }

        private List<string> GetDatasetNames(IEnumerable<DiaNNSearchResult> filteredSearchResults)
        {
            return (from item in filteredSearchResults select item.DatasetName).Distinct().ToList();
        }

        private string GetDatasetName(string fileName)
        {
            var psmSuffixWithoutExtension = Path.GetFileNameWithoutExtension(PSM_FILE_SUFFIX);

            var baseName = Path.GetFileNameWithoutExtension(fileName);

            return baseName.EndsWith(psmSuffixWithoutExtension, StringComparison.OrdinalIgnoreCase)
                ? baseName.Substring(0, baseName.Length - psmSuffixWithoutExtension.Length)
                : baseName;
        }

        /// <summary>
        /// Examine the dataset names in filteredSearchResults
        /// Create a mapping from full name to abbreviated name
        /// </summary>
        /// <param name="inputFileName"></param>
        /// <param name="filteredSearchResults"></param>
        /// <param name="longestCommonBaseName"></param>
        /// <returns>Dictionary where keys are dataset names and values are abbreviated names</returns>
        private Dictionary<string, string> GetDatasetNameMap(
            string inputFileName,
            IEnumerable<DiaNNSearchResult> filteredSearchResults,
            out string longestCommonBaseName)
        {
            var datasetNames = new SortedSet<string>();

            foreach (var item in filteredSearchResults)
            {
                var datasetName = item.DatasetName;
                if (string.IsNullOrWhiteSpace(datasetName))
                    continue;

                // Note that .Add() calls .AddIfNotPresent() internally, so it is safe to call for dataset names already in the SortedSet
                datasetNames.Add(datasetName);
            }

            return GetDatasetNameMap(inputFileName, datasetNames, out longestCommonBaseName, PSM_FILE_SUFFIX);
        }

        private List<DiaNNModInfo> GetPeptideModifications(DiaNNResults searchResult)
        {
            return GetPeptideModifications(searchResult.PeptideCleanSequence, searchResult.Modifications);
        }

        private List<DiaNNModInfo> GetPeptideModifications(ToolResultsBaseClass searchResult)
        {
            return GetPeptideModifications(searchResult.Sequence, searchResult.ModificationList);
        }

        private List<DiaNNModInfo> GetPeptideModifications(string cleanSequence, string modificationList)
        {
            // modificationList should have both the static and dynamic mods, as a comma separated list
            // of residue number, residue symbol, and mod mass; examples:
            //   15M(15.9949)
            //   1M(15.9949), 5C(57.0215)
            //   N-term(42.0106)

            var mods = new List<DiaNNModInfo>();
            var finalResidueLoc = cleanSequence.Length;

            foreach (var modEntry in modificationList.Split(','))
            {
                // Parse out the residue and mod mass
                var residueMatch = mModListResidueModMatcher.Match(modEntry);
                var termMatch = mModListTerminalModMatcher.Match(modEntry);

                if (!(residueMatch.Success || termMatch.Success))
                {
                    ReportError("Invalid DiaNN mod entry format; must be residue number, symbol, and mod mass: " + modEntry);
                    continue;
                }

                var modMassText = residueMatch.Success
                    ? residueMatch.Groups["ModMass"].Value
                    : termMatch.Groups["ModMass"].Value;

                if (!double.TryParse(modMassText, out var modMass))
                {
                    ReportError(string.Format("Unable to parse the mod mass value from {0}; invalid number: {1}", modEntry, modMassText));
                    continue;
                }

                var currentMod = new DiaNNModInfo
                {
                    ModMass = modMass
                };

                if (residueMatch.Success)
                {
                    // Matched a modified residue
                    var residueNumber = residueMatch.Groups["ResidueNumber"].Value;
                    currentMod.ResidueSymbol = residueMatch.Groups["ResidueSymbol"].Value[0];

                    if (!int.TryParse(residueNumber, out currentMod.ResidueLocInPeptide))
                    {
                        ReportError("Unable to parse the residue number from the mod entry: " + modEntry);
                        continue;
                    }

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
                }
                else
                {
                    // Matched a terminal mod
                    switch (termMatch.Groups["TerminusName"].Value)
                    {
                        case "N-term":
                            currentMod.ResidueLocInPeptide = 1;
                            currentMod.TerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                            currentMod.ResidueSymbol = cleanSequence[0];
                            break;

                        case "C-term":
                            currentMod.ResidueLocInPeptide = finalResidueLoc;
                            currentMod.TerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                            currentMod.ResidueSymbol = cleanSequence[finalResidueLoc - 1];
                            break;

                        default:
                            ReportError("Unrecognized terminus name in the mod entry: " + modEntry);
                            continue;
                    }
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
            }

            return mods;
        }

        private string GetPeptideSequence(ToolResultsBaseClass peptideInfo, bool includePrefixAndSuffix = true)
        {
            if (!includePrefixAndSuffix)
                return peptideInfo.Sequence;

            return peptideInfo.PrefixResidue + "." + peptideInfo.Sequence + "." + peptideInfo.SuffixResidue;
        }

        private bool GetTsvColumnNameSynonym(
            IReadOnlyDictionary<string, DiaNNDatasetTsvFileColumns> tsvFileColumnNames,
            string columnName,
            out DiaNNReportFileColumns psmFileColumn)
        {
            if (!tsvFileColumnNames.TryGetValue(columnName, out var resultFileColumn))
            {
                psmFileColumn = DiaNNReportFileColumns.Undefined;
                return false;
            }

            psmFileColumn = resultFileColumn switch
            {
                DiaNNDatasetTsvFileColumns.ScanNumber => DiaNNReportFileColumns.Spectrum,
                DiaNNDatasetTsvFileColumns.Peptide => DiaNNReportFileColumns.Peptide,
                DiaNNDatasetTsvFileColumns.PeptidePrevAA => DiaNNReportFileColumns.PrevAA,
                DiaNNDatasetTsvFileColumns.PeptideNextAA => DiaNNReportFileColumns.NextAA,
                DiaNNDatasetTsvFileColumns.Charge => DiaNNReportFileColumns.Charge,
                DiaNNDatasetTsvFileColumns.RetentionTimeMinutes => DiaNNReportFileColumns.RetentionTime,
                DiaNNDatasetTsvFileColumns.PrecursorNeutralMass => DiaNNReportFileColumns.ObservedMass,
                DiaNNDatasetTsvFileColumns.CalcNeutralPepMass => DiaNNReportFileColumns.CalculatedPeptideMass,
                DiaNNDatasetTsvFileColumns.MassDiff => DiaNNReportFileColumns.DeltaMass,
                DiaNNDatasetTsvFileColumns.ExpectScore => DiaNNReportFileColumns.Expectation,
                DiaNNDatasetTsvFileColumns.Hyperscore => DiaNNReportFileColumns.Hyperscore,
                DiaNNDatasetTsvFileColumns.NextScore => DiaNNReportFileColumns.Nextscore,
                DiaNNDatasetTsvFileColumns.NumTolTerm => DiaNNReportFileColumns.NumberOfEnzymaticTermini,
                DiaNNDatasetTsvFileColumns.NumMissedCleavages => DiaNNReportFileColumns.NumberOfMissedCleavages,
                DiaNNDatasetTsvFileColumns.NumMatchedIons => DiaNNReportFileColumns.NumberOfMatchedIons,
                DiaNNDatasetTsvFileColumns.TotNumIons => DiaNNReportFileColumns.TotalNumberOfIons,
                DiaNNDatasetTsvFileColumns.ModificationInfo => DiaNNReportFileColumns.AssignedModifications,
                DiaNNDatasetTsvFileColumns.Protein => DiaNNReportFileColumns.Protein,
                DiaNNDatasetTsvFileColumns.AlternativeProteins => DiaNNReportFileColumns.MappedProteins,
                _ => DiaNNReportFileColumns.Undefined
            };

            return psmFileColumn != DiaNNReportFileColumns.Undefined;
        }

        private bool LoadNonPsmResults(FileInfo inputFile, ICollection<string> errorMessages, List<DiaNNSearchResult> filteredSearchResults)
        {
            if (!inputFile.Name.EndsWith(PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase))
            {
                throw new ArgumentOutOfRangeException(
                    nameof(inputFile),
                    string.Format("The file name does not end with {0}; invalid call to method LoadNonPsmResults", inputFile.Name));
            }

            if (filteredSearchResults.Count > 0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(filteredSearchResults), "LoadNonPsmResults should only be called if there are no search results");
            }

            try
            {
                if (inputFile.Directory == null)
                {
                    SetErrorMessage("Unable to determine the parent directory of file " + inputFile.FullName);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                Dictionary<FileInfo, string> additionalResultFiles;

                if (inputFile.Name.Equals("Aggregation_psm.tsv", StringComparison.OrdinalIgnoreCase))
                {
                    additionalResultFiles = FindAggregationPsmSourceFiles(inputFile, new List<string>());
                }
                else
                {
                    additionalResultFiles = new Dictionary<FileInfo, string>();

                    var datasetName = GetDatasetName(inputFile.Name);
                    var tsvToFind = new FileInfo(Path.Combine(inputFile.Directory.FullName, datasetName + ".tsv"));

                    if (tsvToFind.Exists)
                    {
                        additionalResultFiles.Add(tsvToFind, datasetName);
                    }
                }

                if (additionalResultFiles.Count == 0)
                {
                    OnWarningEvent("The _psm.tsv file does not have any filter passing data; unable to find any Dataset.tsv files to load instead");
                    return false;
                }

                OnWarningEvent("The _psm.tsv file does not have any filter passing data; loading data from {0} instead",
                    additionalResultFiles.Count == 1
                        ? additionalResultFiles.First().Key.Name
                        : string.Format("{0} Dataset.tsv files", additionalResultFiles.Count));

                var successCountAdditionalFiles = 0;

                foreach (var additionalFile in additionalResultFiles)
                {
                    if (!ReadDiaNNResults(additionalFile.Key, errorMessages, out var additionalSearchResults))
                        continue;

                    successCountAdditionalFiles++;
                    filteredSearchResults.AddRange(additionalSearchResults);
                }

                ComputeQValues(filteredSearchResults);

                return successCountAdditionalFiles >= additionalResultFiles.Count;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadNonPsmResults", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read the precursor match tolerance from the DiaNN parameter file
        /// </summary>
        /// <param name="DiaNNParamFilePath"></param>
        /// <returns>True on success, false if an error</returns>
        private bool LoadSearchEngineParamFile(string DiaNNParamFilePath)
        {
            try
            {
                var sourceFile = new FileInfo(DiaNNParamFilePath);
                if (!sourceFile.Exists)
                {
                    SetErrorMessage("DiaNN parameter file not found: " + DiaNNParamFilePath);
                    SetErrorCode(PHRPErrorCode.ParameterFileNotFound);
                    return false;
                }

                OnStatusEvent("Reading the DiaNN parameter file: " + PathUtils.CompactPathString(sourceFile.FullName, 110));

                var startupOptions = new StartupOptions
                {
                    DisableOpeningInputFiles = true,
                    LoadModsAndSeqInfo = false
                };

                var reader = new DiaNNSynFileReader("DiaNN_ParamFile_Reader", sourceFile.FullName, startupOptions);
                RegisterEvents(reader);

                var success = reader.LoadSearchEngineParameters(sourceFile.FullName, out var searchEngineParams);

                if (!success)
                {
                    return false;
                }

                var validTolerances = reader.GetPrecursorSearchTolerances(
                    searchEngineParams,
                    out var toleranceLower, out var toleranceUpper,
                    out var ppmBased, out var singleTolerance);

                if (!validTolerances)
                {
                    OnWarningEvent("Unable to extract the precursor ion match tolerances from the parameter file");
                    mPrecursorMassTolerance = new PrecursorMassTolerance(75, true);
                    return false;
                }

                if (singleTolerance)
                {
                    mPrecursorMassTolerance = new PrecursorMassTolerance(toleranceLower, ppmBased);
                }
                else
                {
                    mPrecursorMassTolerance = new PrecursorMassTolerance(toleranceLower, toleranceUpper, ppmBased);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadSearchEngineParamFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        private bool MergeRelatedDiaNNResults(
            FileInfo inputFile,
            ICollection<string> errorMessages,
            List<DiaNNSearchResult> filteredSearchResults)
        {
            try
            {
                if (inputFile.Directory == null)
                {
                    return false;
                }

                // Check whether the input file ends with _psm.tsv
                var targetResultsAreFromPsmFile = inputFile.Name.EndsWith(PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase);

                // Keys are FileInfo objects, values are the dataset name for the given file
                Dictionary<FileInfo, string> additionalResultFiles;

                var datasetNames = GetDatasetNames(filteredSearchResults);

                if (datasetNames.Count > 1 || inputFile.Name.Equals("Aggregation_psm.tsv", StringComparison.OrdinalIgnoreCase))
                {
                    additionalResultFiles = FindAggregationPsmSourceFiles(inputFile, datasetNames);
                }
                else
                {
                    var datasetName = GetDatasetName(inputFile.Name);

                    var nameToFind = targetResultsAreFromPsmFile
                        ? datasetName + ".tsv"
                        : datasetName + PSM_FILE_SUFFIX;

                    var tsvToFind = new FileInfo(Path.Combine(inputFile.Directory.FullName, nameToFind));

                    additionalResultFiles = new Dictionary<FileInfo, string>();

                    if (tsvToFind.Exists)
                    {
                        additionalResultFiles.Add(tsvToFind, datasetName);
                    }
                }

                if (additionalResultFiles.Count == 0)
                {
                    if (targetResultsAreFromPsmFile)
                    {
                        OnDebugEvent("Unable to find the Dataset.tsv file that corresponds to {0}; synopsis file will not have matched ion counts", inputFile.Name);
                    }
                    else
                    {
                        OnDebugEvent("Unable to find the Dataset_psm.tsv file that corresponds to {0}; synopsis file will not have peptide prophet probabilities", inputFile.Name);
                    }

                    // This is not a fatal error, so return true
                    return true;
                }

                // additionalResultFiles may have a mix of both Dataset.tsv and Dataset_psm.tsv files
                // Preferentially use the Dataset.tsv files

                var additionalResultFilesByDataset = new Dictionary<string, KeyValuePair<FileInfo, string>>(StringComparer.OrdinalIgnoreCase);

                foreach (var additionalFile in additionalResultFiles)
                {
                    if (!additionalResultFilesByDataset.TryGetValue(additionalFile.Value, out var datasetInfo))
                    {
                        additionalResultFilesByDataset.Add(additionalFile.Value, additionalFile);
                        continue;
                    }

                    if (datasetInfo.Key.Name.EndsWith(PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase) &&
                        !additionalFile.Key.Name.EndsWith(PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase))
                    {
                        additionalResultFilesByDataset[additionalFile.Value] = additionalFile;
                    }
                }

                var matchDatasetNames = additionalResultFilesByDataset.Count > 1;

                // This dictionary tracks elution times observed for each scan number, across all datasets
                var elutionTimesByScanNumber = new Dictionary<int, List<double>>();

                foreach (var additionalFile in additionalResultFilesByDataset.Values)
                {
                    if (!ReadDiaNNResults(additionalFile.Key, errorMessages, out var additionalSearchResults))
                        continue;

                    // If the input file is a _psm.tsv file and this additional file is a Dataset.tsv file, compute Q-Values
                    if (targetResultsAreFromPsmFile && !additionalFile.Key.Name.EndsWith(PSM_FILE_SUFFIX, StringComparison.OrdinalIgnoreCase))
                    {
                        // When search results are loaded from a Dataset.tsv file, it should have reverse hits,
                        // allowing for computing Q-Values using additionalSearchResults
                        ComputeQValues(additionalSearchResults);
                    }

                    // Populate a dictionary mapping scan, charge, and peptide to each row in additionalSearchResults
                    var searchResultsLookup = new Dictionary<string, DiaNNSearchResult>();

                    foreach (var additionalResult in additionalSearchResults)
                    {
                        var key = ConstructKeyForPSM(additionalResult);
                        searchResultsLookup.Add(key, additionalResult);

                        if (!double.TryParse(additionalResult.ElutionTime, out var elutionTime))
                            continue;

                        if (elutionTimesByScanNumber.TryGetValue(additionalResult.ScanNum, out var elutionTimes))
                        {
                            elutionTimes.Add(elutionTime);
                            continue;
                        }

                        elutionTimesByScanNumber.Add(additionalResult.ScanNum, new List<double> { elutionTime });
                    }

                    var datasetNameFilter = additionalFile.Value;

                    // Merge information from additionalSearchResults into filteredSearchResults by matching on scan and peptide

                    foreach (var item in filteredSearchResults)
                    {
                        if (matchDatasetNames && item.DatasetName.Length > 0 && !item.DatasetName.Equals(datasetNameFilter, StringComparison.OrdinalIgnoreCase))
                            continue;

                        var key = ConstructKeyForPSM(item);
                        if (!searchResultsLookup.TryGetValue(key, out var additionalResult))
                            continue;

                        UpdateUndefinedValue(ref item.NumberOfMatchedIons, additionalResult.NumberOfMatchedIons);
                        UpdateUndefinedValue(ref item.TotalNumberOfIons, additionalResult.TotalNumberOfIons);
                        UpdateUndefinedValue(ref item.PeptideProphetProbability, additionalResult.PeptideProphetProbability);

                        if (!targetResultsAreFromPsmFile)
                            continue;

                        if (additionalResult.QValue > 0)
                            item.QValue = additionalResult.QValue;

                        if (string.IsNullOrWhiteSpace(item.ElutionTime))
                            item.ElutionTime = additionalResult.ElutionTime;
                    }
                }

                if (!targetResultsAreFromPsmFile || datasetNames.Count < 2)
                    return true;

                // Obtain a dictionary of average elution time by scan number (a moving average smooth is applied to smooth out variations)
                var averageElutionTimeByScanNumber = ComputeScanToElutionTimeMap(elutionTimesByScanNumber);

                foreach (var item in filteredSearchResults)
                {
                    if (!averageElutionTimeByScanNumber.TryGetValue(item.ScanNum, out var elutionTime))
                        continue;

                    var elutionTimeText = PRISM.StringUtilities.DblToString(elutionTime, 4);

                    if (string.IsNullOrWhiteSpace(item.ElutionTime))
                        item.ElutionTime = elutionTimeText;

                    item.ElutionTimeAverage = elutionTimeText;
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in MergeRelatedDiaNNResults", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
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
        /// Parse a DiaNN results line while creating the DiaNN synopsis file
        /// </summary>
        /// <param name="readingPsmFile"></param>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="columnMapping"></param>
        /// <param name="lineNumber">Line number in the input file (used for error reporting)</param>
        /// <param name="currentDatasetName">Current dataset name; updated by this method if the Spectrum File name is not an empty string</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNResultsFileEntry(
            bool readingPsmFile,
            string lineIn,
            out DiaNNSearchResult searchResult,
            ICollection<string> errorMessages,
            IDictionary<DiaNNReportFileColumns, int> columnMapping,
            int lineNumber,
            ref string currentDatasetName)
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

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.SpectrumFile], out string spectrumFile);

                if (!string.IsNullOrWhiteSpace(spectrumFile))
                {
                    // Extract the dataset name from the spectrum file name
                    // Example names:
                    //   interact-Dataset.pep.xml
                    //   interact-Dataset.pin.pep.xml

                    var baseName = Path.GetFileNameWithoutExtension(spectrumFile);
                    if (baseName.StartsWith("interact-", StringComparison.OrdinalIgnoreCase))
                    {
                        currentDatasetName = baseName.Substring("interact-".Length);
                    }
                    else
                    {
                        currentDatasetName = baseName;
                    }

                    var suffixesToRemove = new List<string>
                    {
                        ".pep",
                        ".pin"
                    };

                    while (true)
                    {
                        var suffixRemoved = false;

                        foreach (var suffix in suffixesToRemove)
                        {
                            if (currentDatasetName.EndsWith(suffix))
                            {
                                currentDatasetName = currentDatasetName.Substring(0, currentDatasetName.Length - suffix.Length);
                                suffixRemoved = true;
                            }
                        }

                        if (!suffixRemoved)
                            break;
                    }
                }

                searchResult.DatasetName = currentDatasetName;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Spectrum], out searchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid on line " + lineNumber, true);
                }

                if (!int.TryParse(searchResult.Scan, out searchResult.ScanNum))
                {
                    ReportError("Scan column is not numeric on line " + lineNumber, true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Charge], out searchResult.Charge);
                searchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(searchResult.Charge, 0));

                // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by DiaNN
                if (DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.CalculatedPeptideMass], out searchResult.CalculatedMonoMass))
                {
                    double.TryParse(searchResult.CalculatedMonoMass, out searchResult.CalculatedMonoMassValue);
                }

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Peptide], out searchResult.Sequence))
                {
                    ReportError("Peptide column is missing or invalid on line " + lineNumber, true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ModifiedPeptide], out searchResult.ModifiedPeptide);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrevAA], out searchResult.PrefixResidue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NextAA], out searchResult.SuffixResidue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PeptideLength], out searchResult.Length);

                if (searchResult.Length == 0)
                {
                    searchResult.Length = searchResult.Sequence.Length;
                }

                if (readingPsmFile)
                {
                    // The aggregated results file reports retention time in seconds
                    if (DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.RetentionTime], out double retentionTimeSeconds))
                    {
                        searchResult.ElutionTime = PRISM.StringUtilities.DblToString(retentionTimeSeconds / 60.0, 4);
                    }
                }
                else
                {
                    // The single dataset results file reports retention time in minutes
                    DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.RetentionTime], out searchResult.ElutionTime);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ObservedMass], out searchResult.PrecursorMonoMass);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.CalibratedObservedMass], out searchResult.CalibratedObservedMass);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ObservedMZ], out searchResult.PrecursorMZ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.CalibratedObservedMZ], out searchResult.CalibratedObservedMZ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.CalculatedMZ], out searchResult.CalculatedMZ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.DeltaMass], out searchResult.MassErrorDaDiaNN);

                if (string.IsNullOrWhiteSpace(searchResult.PrecursorMZ) && double.TryParse(searchResult.PrecursorMonoMass, out var precursorMonoMass))
                {
                    var observedPrecursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, searchResult.ChargeNum);
                    searchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(observedPrecursorMZ, 4);
                }

                if (string.IsNullOrWhiteSpace(searchResult.CalculatedMZ) && searchResult.CalculatedMonoMassValue > 0)
                {
                    var calculatedPrecursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(searchResult.CalculatedMonoMassValue, 0, searchResult.ChargeNum);
                    searchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(calculatedPrecursorMZ, 4);
                }

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                searchResult.MH = ComputeMH(searchResult);

                if (DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Expectation], out searchResult.Expectation))
                {
                    double.TryParse(searchResult.Expectation, out searchResult.EValue);
                }

                if (DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Hyperscore], out searchResult.Hyperscore))
                {
                    double.TryParse(searchResult.Hyperscore, out searchResult.HyperscoreValue);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Nextscore], out searchResult.Nextscore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PeptideProphetProbability], out searchResult.PeptideProphetProbability);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NumberOfEnzymaticTermini], out searchResult.NumberOfTrypticTermini);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NumberOfMissedCleavages], out searchResult.MissedCleavageCount);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NumberOfMatchedIons], out searchResult.NumberOfMatchedIons);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.TotalNumberOfIons], out searchResult.TotalNumberOfIons);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinStart], out searchResult.ProteinStart);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinEnd], out searchResult.ProteinEnd);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Intensity], out searchResult.Intensity, "0");

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.AssignedModifications], out searchResult.ModificationList);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ObservedModifications], out searchResult.ObservedModifications);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.IsUnique], out searchResult.IsUnique);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Protein], out searchResult.Protein);

                // The Protein column in Dataset.tsv files has both protein name and description; extract out the protein description
                // For _psm.tsv files, .Protein should just have a single protein name, but we'll check for a space anyway

                if (SplitProteinNameAndDescription(searchResult.Protein, out var proteinName, out var proteinDescription))
                {
                    searchResult.Protein = proteinName;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinID], out searchResult.ProteinID);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.EntryName], out searchResult.EntryName);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Gene], out searchResult.Gene);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinDescription], out searchResult.ProteinDescription, proteinDescription);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MappedGenes], out searchResult.AdditionalGenes);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MappedProteins], out searchResult.AdditionalProteins);

                if (string.IsNullOrWhiteSpace(searchResult.IsUnique))
                {
                    searchResult.IsUnique = searchResult.AdditionalProteins.Length > 0 ? "FALSE" : "TRUE";
                }

                if (!readingPsmFile && searchResult.AdditionalProteins.Length > 0 ||
                    searchResult.AdditionalProteins.IndexOf("@@", StringComparison.Ordinal) > 0)
                {
                    // The list of additional proteins will have both protein names and descriptions, each separated by @@

                    var additionalProteinAndDescriptions = searchResult.AdditionalProteins.Split(new[] { "@@" }, StringSplitOptions.RemoveEmptyEntries);
                    var additionalProteins = new List<string>();

                    foreach (var item in additionalProteinAndDescriptions)
                    {
                        SplitProteinNameAndDescription(item, out var additionalProteinName, out _);

                        additionalProteins.Add(additionalProteinName);
                    }

                    searchResult.AdditionalProteins = string.Join(", ", additionalProteins);
                }

                searchResult.Reverse = IsReversedProtein(searchResult.Protein);

                // Parse the modification list to determine the total mod mass
                var totalModMass = ComputeTotalModMass(searchResult);

                // Compute monoisotopic mass of the peptide
                searchResult.CalculatedMonoMassPHRP = ComputePeptideMassForCleanSequence(searchResult.Sequence, totalModMass);

                // MassErrorPpm and MassErrorDa will be computed later by method ComputeObservedMassErrors

                // Compare the mass computed by PHRP to the one reported by DiaNN
                var deltaMassVsDiaNN = searchResult.CalculatedMonoMassPHRP - searchResult.CalculatedMonoMassValue;

                if (Math.Abs(deltaMassVsDiaNN) > 0.01)
                {
                    OnWarningEvent("Calculated monoisotopic mass differs from the value reported by DiaNN on line {0} for {1}, {2} Da; delta mass: {3:F3} Da", lineNumber, searchResult.Sequence, searchResult.CalculatedMonoMassValue, deltaMassVsDiaNN);
                }

                var peptideWithPrefixAndSuffix = GetPeptideSequence(searchResult);

                // ReSharper disable once CommentTypo

                // The missed cleavage count listed in _psm.tsv files is sometimes wrong
                // For example, peptide R.TEMENEFVLIKK.D has a missed cleavage count of 1 in the Dataset.tsv file but 0 in the Dataset_psm.tsv file
                // Use the peptide cleavage state calculator to compute the value and update if required

                var missedCleavageCount = mPeptideCleavageStateCalculator.ComputeNumberOfMissedCleavages(peptideWithPrefixAndSuffix);

                if (int.TryParse(searchResult.MissedCleavageCount, out var missedCleavageCountDiaNN) && missedCleavageCount != missedCleavageCountDiaNN)
                {
                    searchResult.MissedCleavageCount = missedCleavageCount.ToString();
                }

                var cleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(peptideWithPrefixAndSuffix);

                var ntt = cleavageState switch
                {
                    PeptideCleavageStateCalculator.PeptideCleavageState.Full => 2,
                    PeptideCleavageStateCalculator.PeptideCleavageState.Partial => 1,
                    PeptideCleavageStateCalculator.PeptideCleavageState.NonSpecific => 0,
                    PeptideCleavageStateCalculator.PeptideCleavageState.Unknown => 0,
                    _ => 0
                };

                if (cleavageState == PeptideCleavageStateCalculator.PeptideCleavageState.Unknown || searchResult.NumberOfTrypticTermini == ntt)
                {
                    return true;
                }

                if (ntt <= searchResult.NumberOfTrypticTermini)
                {
                    // If a peptide starts after the M at the start of a protein,
                    // DiaNN treats this as valid for a fully tryptic peptide

                    // DiaNN also does not use the proline rule when determining cleavage state

                    // Thus, allow NumberOfTrypticTermini to be larger than the computed cleavage state value
                    return true;
                }

                OnWarningEvent("Changing cleavage state from {0} to {1} for {2} in scan {3}",
                    searchResult.NumberOfTrypticTermini, ntt, peptideWithPrefixAndSuffix, searchResult.Scan);

                searchResult.NumberOfTrypticTermini = ntt;

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the MassDiaNN results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing DiaNN results in ParseDiaNNResultsFileEntry, line {0}: {1}", lineNumber, ex.Message));
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the DiaNN results file header line, populating columnMapping
        /// </summary>
        /// <remarks>
        /// This method supports both Dataset_psm.tsv files and Dataset.tsv (which have similar, but different column names)
        /// </remarks>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <param name="readingPsmFile"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseDiaNNResultsFileHeaderLine(
            string lineIn,
            IDictionary<DiaNNReportFileColumns, int> columnMapping,
            bool readingPsmFile)
        {
            // Columns in _psm.tsv files
            var columnNames = new SortedDictionary<string, DiaNNReportFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "Spectrum", DiaNNReportFileColumns.Spectrum },
                { "Spectrum File", DiaNNReportFileColumns.SpectrumFile },
                { "Peptide", DiaNNReportFileColumns.Peptide },
                { "Modified Peptide", DiaNNReportFileColumns.ModifiedPeptide },
                { "Prev AA", DiaNNReportFileColumns.PrevAA },
                { "Next AA", DiaNNReportFileColumns.NextAA },
                { "Peptide Length", DiaNNReportFileColumns.PeptideLength },
                { "Charge", DiaNNReportFileColumns.Charge },
                { "Retention", DiaNNReportFileColumns.RetentionTime },
                { "Observed Mass", DiaNNReportFileColumns.ObservedMass },
                { "Calibrated Observed Mass", DiaNNReportFileColumns.CalibratedObservedMass },
                { "Observed M/Z", DiaNNReportFileColumns.ObservedMZ },
                { "Calibrated Observed M/Z", DiaNNReportFileColumns.CalibratedObservedMZ },
                { "Calculated Peptide Mass", DiaNNReportFileColumns.CalculatedPeptideMass },
                { "Calculated M/Z", DiaNNReportFileColumns.CalculatedMZ },
                { "Delta Mass", DiaNNReportFileColumns.DeltaMass },
                { "Expectation", DiaNNReportFileColumns.Expectation },
                { "Hyperscore", DiaNNReportFileColumns.Hyperscore },
                { "Nextscore", DiaNNReportFileColumns.Nextscore },
                { "PeptideProphet Probability", DiaNNReportFileColumns.PeptideProphetProbability },
                { "Number of Enzymatic Termini", DiaNNReportFileColumns.NumberOfEnzymaticTermini },
                { "Number of Missed Cleavages", DiaNNReportFileColumns.NumberOfMissedCleavages },
                { "Protein Start", DiaNNReportFileColumns.ProteinStart },
                { "Protein End", DiaNNReportFileColumns.ProteinEnd },
                { "Intensity", DiaNNReportFileColumns.Intensity },
                { "Assigned Modifications", DiaNNReportFileColumns.AssignedModifications },
                { "Observed Modifications", DiaNNReportFileColumns.ObservedModifications },
                { "Is Unique", DiaNNReportFileColumns.IsUnique },
                { "Protein", DiaNNReportFileColumns.Protein },
                { "Protein ID", DiaNNReportFileColumns.ProteinID },
                { "Entry Name", DiaNNReportFileColumns.EntryName },
                { "Gene", DiaNNReportFileColumns.Gene },
                { "Protein Description", DiaNNReportFileColumns.ProteinDescription },
                { "Mapped Genes", DiaNNReportFileColumns.MappedGenes },
                { "Mapped Proteins", DiaNNReportFileColumns.MappedProteins }
            };

            // ReSharper disable StringLiteralTypo

            var tsvFileColumnNames = new SortedDictionary<string, DiaNNDatasetTsvFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "scannum", DiaNNDatasetTsvFileColumns.ScanNumber },
                { "precursor_neutral_mass", DiaNNDatasetTsvFileColumns.PrecursorNeutralMass },
                { "retention_time", DiaNNDatasetTsvFileColumns.RetentionTimeMinutes },
                { "charge", DiaNNDatasetTsvFileColumns.Charge },
                { "hit_rank", DiaNNDatasetTsvFileColumns.HitRank },
                { "peptide", DiaNNDatasetTsvFileColumns.Peptide },
                { "peptide_prev_aa", DiaNNDatasetTsvFileColumns.PeptidePrevAA },
                { "peptide_next_aa", DiaNNDatasetTsvFileColumns.PeptideNextAA },
                { "protein", DiaNNDatasetTsvFileColumns.Protein },
                { "num_matched_ions", DiaNNDatasetTsvFileColumns.NumMatchedIons },
                { "tot_num_ions", DiaNNDatasetTsvFileColumns.TotNumIons },
                { "calc_neutral_pep_mass", DiaNNDatasetTsvFileColumns.CalcNeutralPepMass },
                { "massdiff", DiaNNDatasetTsvFileColumns.MassDiff },
                { "num_tol_term", DiaNNDatasetTsvFileColumns.NumTolTerm },
                { "num_missed_cleavages", DiaNNDatasetTsvFileColumns.NumMissedCleavages },
                { "modification_info", DiaNNDatasetTsvFileColumns.ModificationInfo },
                { "hyperscore", DiaNNDatasetTsvFileColumns.Hyperscore },
                { "nextscore", DiaNNDatasetTsvFileColumns.NextScore },
                { "expectscore", DiaNNDatasetTsvFileColumns.ExpectScore },
                { "best_locs", DiaNNDatasetTsvFileColumns.BestLocs },
                { "score_without_delta_mass", DiaNNDatasetTsvFileColumns.ScoreWithoutDeltaMass },
                { "best_score_with_delta_mass", DiaNNDatasetTsvFileColumns.BestScoreWithDeltaMass },
                { "second_best_score_with_delta_mass", DiaNNDatasetTsvFileColumns.SecondBestScoreWithDeltaMass },
                { "delta_score", DiaNNDatasetTsvFileColumns.DeltaScore },
                { "alternative_proteins", DiaNNDatasetTsvFileColumns.AlternativeProteins }
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
                    if (!readingPsmFile && GetTsvColumnNameSynonym(tsvFileColumnNames, splitLine[index], out var psmFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[psmFileColumn] = index;
                        continue;
                    }

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
        /// Parse the header line of a DiaNN _syn.txt file, populating columnMapping
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
        /// Parses an entry from a DiaNN Synopsis file
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

                // Note that DiaNN peptides don't actually have mod symbols; that information is tracked via searchResult.Modifications
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
        /// <param name="inputFilePath">DiaNN results file (Dataset_psm.tsv or Dataset.tsv)</param>
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
                        OnWarningEvent("DiaNN parameter file not defined; unable to determine the precursor match tolerance");
                        mPrecursorMassTolerance = new PrecursorMassTolerance(75, true);
                    }
                    else
                    {
                        var DiaNNParameterFilePath = ResolveFilePath(inputFile.DirectoryName, Options.SearchToolParameterFilePath);

                        // Examine the DiaNN parameter file to determine the precursor match tolerance
                        var toleranceExtracted = LoadSearchEngineParamFile(DiaNNParameterFilePath);
                        if (!toleranceExtracted)
                        {
                            return false;
                        }
                    }

                    // Do not create a first-hits file for DiaNN results

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
                            // Use a higher match error threshold since some peptides reported by DiaNN don't perfectly match the FASTA file
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
        /// Load DiaNN search results
        /// </summary>
        /// <remarks>If the file is found, but has no results, this method still returns true</remarks>
        /// <param name="inputFile">report.tsv file</param>
        /// <param name="errorMessages"></param>
        /// <param name="filteredSearchResults">Output: DiaNN results</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ReadDiaNNResults(
            FileSystemInfo inputFile,
            ICollection<string> errorMessages,
            out List<DiaNNSearchResult> filteredSearchResults)
        {
            filteredSearchResults = new List<DiaNNSearchResult>();

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

                // Initialize the list that will hold all of the records in the DiaNN result file
                var searchResultsUnfiltered = new List<DiaNNSearchResult>();

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
                        lineNumber);

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

                // Sort the SearchResults by dataset name, scan, charge, and ascending E-Value
                searchResultsUnfiltered.Sort(new DiaNNSearchResultsComparerDatasetScanChargeEValuePeptide());

                // Now filter the data
                var startIndex = 0;

                while (startIndex < searchResultsUnfiltered.Count)
                {
                    // Find all of the matches for the current result's scan
                    // (we sorted by dataset, then scan, so adjacent results will be from the same dataset, except when a new dataset is encountered)
                    // DiaNN will typically report just one match

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
        /// If proteinNameToCheck contains a space, extract the protein name and description, then return true
        /// Otherwise, store the string in proteinName, set proteinDescription to an empty string, and return false
        /// </summary>
        /// <param name="proteinNameToCheck"></param>
        /// <param name="proteinName"></param>
        /// <param name="proteinDescription"></param>
        /// <returns>True if a space was found, otherwise false</returns>
        private bool SplitProteinNameAndDescription(string proteinNameToCheck, out string proteinName, out string proteinDescription)
        {
            var trimmedName = proteinNameToCheck.Trim();
            var spaceIndex = trimmedName.IndexOf(' ');

            if (spaceIndex > 0)
            {
                proteinDescription = spaceIndex < trimmedName.Length - 1
                    ? trimmedName.Substring(spaceIndex + 1).Trim()
                    : string.Empty;

                proteinName = trimmedName.Substring(0, spaceIndex);
                return true;
            }

            proteinName = trimmedName;
            proteinDescription = string.Empty;
            return false;
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan in a single dataset (typically there will only be one result for DiaNN)
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
            //  Either E-value < 0.75
            //  or     HyperscoreValue > 20
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].EValue < Options.MSGFPlusSynopsisFileEValueThreshold ||
                    searchResults[index].HyperscoreValue > Options.DiaNNHyperscoreThreshold)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Update targetValue if it is empty, but sourceValue is not empty
        /// </summary>
        /// <param name="targetValue"></param>
        /// <param name="sourceValue"></param>
        private void UpdateUndefinedValue(ref string targetValue, string sourceValue)
        {
            if (!string.IsNullOrWhiteSpace(targetValue) || string.IsNullOrWhiteSpace(sourceValue))
                return;

            targetValue = sourceValue;
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
        /// <param name="baseNameByDatasetName">Keys are dataset names, values are dataset ID (or 0 if undefined)</param>
        /// <param name="writer"></param>
        /// <param name="filteredSearchResults"></param>
        /// <param name="errorMessages"></param>
        private void WriteFilteredSearchResults(
            Dictionary<string, string> baseNameByDatasetName,
            TextWriter writer,
            List<DiaNNSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            // Lookup the Dataset ID for each dataset (only if on the pnl.gov domain)
            var datasetIDs = LookupDatasetIDs(baseNameByDatasetName.Keys.ToList());

            var index = 1;
            foreach (var result in filteredSearchResults)
            {
                GetBaseNameAndDatasetID(baseNameByDatasetName, datasetIDs, result.DatasetName, out var baseDatasetName, out var datasetID);

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
                    GetPeptideSequence(searchResult),
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

        private class DiaNNSearchResultsComparerDatasetScanChargeEValuePeptide : IComparer<DiaNNSearchResult>
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

                // Charge is the same; check E-value
                if (x.EValue < y.EValue)
                {
                    return -1;
                }

                if (x.EValue > y.EValue)
                {
                    return 1;
                }

                // E-value is the same; check sequence
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

        private class DiaNNSearchResultsComparerEValueScanChargePeptide : IComparer<DiaNNSearchResult>
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

                // E-value is the same; check scan number
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
