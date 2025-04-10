﻿// This class reads search results from a DIA-NN report.tsv file and creates
// a tab-delimited text file with the data
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using ParquetSharp;
using PeptideHitResultsProcessor.Data;
using PeptideHitResultsProcessor.SearchToolResults;
using PeptideToProteinMapEngine;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM;
using ProteinCoverageSummarizer;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads search results from a DIA-NN report.tsv file and creates
    /// a tab-delimited text file with the data
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads a DIA-NN results file (Dataset_report.tsv)
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

        // Ignore Spelling: Averagine, Dia, diann, proteotypic, Quant, Qvalue, tryptic, txt

        // ReSharper restore CommentTypo

        /// <summary>
        /// Constructor
        /// </summary>
        public DiaNNResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "March 22, 2025";

            mMissingScanInfoDatasets = new SortedSet<string>(StringComparer.OrdinalIgnoreCase);
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
        /// N-terminus symbol used by DIA-NN
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_DiaNN = "-";

        /// <summary>
        /// C-terminus symbol used by DIA-NN
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_DiaNN = "-";

        /// <summary>
        /// .parquet report file suffix
        /// </summary>
        public const string PARQUET_REPORT_FILE_SUFFIX = "_report.parquet";

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        /// <summary>
        /// These columns correspond to the DIA-NN report.tsv file
        /// </summary>
        internal enum DiaNNReportFileColumns
        {
            // ReSharper disable InconsistentNaming

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
            ///   C(UniMod:4)TSYNIPC(UniMod:4)TSDM(UniMod:35)AK
            ///   DIHFM(UniMod:35)PC(UniMod:4)SGLTGANLK
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
            /// Precursor m/z value
            /// </summary>
            PrecursorMz = 17,

            /// <summary>
            /// QValue
            /// </summary>
            QValue = 18,

            /// <summary>
            /// PEP (posterior error probability)
            /// </summary>
            PEP = 19,

            /// <summary>
            /// Global QValue
            /// </summary>
            GlobalQValue = 20,

            /// <summary>
            /// Protein QValue
            /// </summary>
            ProteinQValue = 21,

            /// <summary>
            /// Protein Group QValue
            /// </summary>
            ProteinGroupQValue = 22,

            /// <summary>
            /// Global Protein Group QValue
            /// </summary>
            GlobalProteinGroupQValue = 23,

            /// <summary>
            /// Gene Group QValue
            /// </summary>
            GeneGroupQValue = 24,

            /// <summary>
            /// Translated QValue
            /// </summary>
            TranslatedQValue = 25,

            /// <summary>
            /// 1 if the peptide is specific to a given protein or gene, otherwise 0
            /// </summary>
            Proteotypic = 26,

            /// <summary>
            /// Precursor Quantity
            /// </summary>
            PrecursorQuantity = 27,

            /// <summary>
            /// Precursor Normalized
            /// </summary>
            PrecursorNormalized = 28,

            /// <summary>
            /// Precursor Translated
            /// </summary>
            PrecursorTranslated = 29,

            /// <summary>
            /// Translated Quality
            /// </summary>
            TranslatedQuality = 30,

            /// <summary>
            /// MS1 Translated
            /// </summary>
            MS1Translated = 31,

            /// <summary>
            /// Quantity Quality
            /// </summary>
            QuantityQuality = 32,

            /// <summary>
            /// Elution Time (aka retention time)
            /// </summary>
            RT = 33,

            /// <summary>
            /// Elution Time Start
            /// </summary>
            RTStart = 34,

            /// <summary>
            /// Elution Time Stop
            /// </summary>
            RTStop = 35,

            /// <summary>
            /// Indexed retention time (iRT)
            /// </summary>
            IndexedRT = 36,

            /// <summary>
            /// Predicted retention time (unused by PHRP)
            /// </summary>
            PredictedRT = 37,

            /// <summary>
            /// Predicted indexed retention time (unused by PHRP)
            /// </summary>
            PredictedIndexedRT = 38,

            /// <summary>
            /// First protein description (unused by PHRP)
            /// </summary>
            FirstProteinDescription = 39,

            /// <summary>
            /// Spectral Library QValue
            /// </summary>
            LibQValue = 40,

            /// <summary>
            /// Spectral Library Protein Group QValue
            /// </summary>
            LibProteinGroupQValue = 41,

            /// <summary>
            /// MS1 Profile Correlation
            /// </summary>
            MS1ProfileCorr = 42,

            /// <summary>
            /// MS1 Area
            /// </summary>
            MS1Area = 43,

            /// <summary>
            /// Evidence (score)
            /// </summary>
            /// <remarks>Higher values are better</remarks>
            Evidence = 44,

            /// <summary>
            /// Spectrum Similarity
            /// </summary>
            SpectrumSimilarity = 45,

            /// <summary>
            /// Averagine
            /// </summary>
            Averagine = 46,

            /// <summary>
            /// Mass Evidence
            /// </summary>
            /// <remarks>Higher values are better</remarks>
            MassEvidence = 47,

            /// <summary>
            /// CScore
            /// </summary>
            CScore = 48,

            /// <summary>
            /// Decoy Evidence
            /// </summary>
            DecoyEvidence = 49,

            /// <summary>
            /// Decoy CScore
            /// </summary>
            DecoyCScore = 50,

            /// <summary>
            /// Fragmentation Ion Abundances
            /// </summary>
            FragmentQuantRaw = 51,

            /// <summary>
            /// Corrected Fragmentation Ion Abundances
            /// </summary>
            FragmentQuantCorrected = 52,

            /// <summary>
            /// Fragmentation ion correlations
            /// </summary>
            FragmentCorrelations = 53,

            /// <summary>
            /// MS/MS Scan Number
            /// </summary>
            MS2Scan = 54,

            /// <summary>
            /// Ion Mobility
            /// </summary>
            IonMobility = 55,

            /// <summary>
            /// Indexed Ion Mobility
            /// </summary>
            IndexedIonMobility = 56,

            /// <summary>
            /// Predicted Ion Mobility
            /// </summary>
            PredictedIonMobility = 57,

            /// <summary>
            /// Predicted Indexed Ion Mobility
            /// </summary>
            PredictedIndexedIonMobility = 58,

            /// <summary>
            /// Run Index
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            RunIndex = 59,

            /// <summary>
            /// Channel
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            Channel = 60,

            /// <summary>
            /// Precursor Library Index
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            PrecursorLibIndex = 61,

            /// <summary>
            /// Decoy (0 if not a decoy, 1 if a decoy)
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            Decoy = 62,

            /// <summary>
            /// Protein Group MaxLFQ Quality
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            ProteinGroupMaxLFQQuality = 63,

            /// <summary>
            /// Genes MaxLFQ Quality
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            GenesMaxLFQQuality = 64,

            /// <summary>
            /// Genes MaxLFQ Unique Quality
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            GenesMaxLFQUniqueQuality = 65,

            /// <summary>
            /// Peptidoform QValue
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            PeptidoformQValue = 66,

            /// <summary>
            /// Global Peptidoform QValue
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            GlobalPeptidoformQValue = 67,

            /// <summary>
            /// Lib Peptidoform QValue
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            LibPeptidoformQValue = 68,

            /// <summary>
            /// Ptm Site Confidence
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            PtmSiteConfidence = 69,

            /// <summary>
            /// Site Occupancy Probabilities
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            SiteOccupancyProbabilities = 70,

            /// <summary>
            /// Protein Sites
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            ProteinSites = 71,

            /// <summary>
            /// Lib Ptm Site Confidence
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            LibPtmSiteConfidence = 72,

            /// <summary>
            /// Channel QValue
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            ChannelQValue = 73,

            /// <summary>
            /// Protein Group PEP
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            ProteinGroupPEP = 74,

            /// <summary>
            /// MS1 Normalised
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            MS1Normalised = 75,

            /// <summary>
            /// MS1 Apex Area
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            MS1ApexArea = 76,

            /// <summary>
            /// MS1 Apex Mz Delta
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            MS1ApexMzDelta = 77,

            /// <summary>
            /// Normalisation Factor
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            NormalisationFactor = 78,

            /// <summary>
            /// Empirical Quality
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            EmpiricalQuality = 79,

            /// <summary>
            /// Normalisation Noise
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            NormalisationNoise = 80,

            /// <summary>
            /// FWHM
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            FWHM = 81,

            /// <summary>
            /// Protein Group Top N
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            ProteinGroupTopN = 82,

            /// <summary>
            /// Genes Top N
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            GenesTopN = 83,

            /// <summary>
            /// Channel Evidence
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            ChannelEvidence = 84,

            /// <summary>
            /// MS1 Total Signal Before
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            MS1TotalSignalBefore = 85,

            /// <summary>
            /// MS1 Total Signal After
            /// </summary>
            /// <remarks>Only in .parquet files</remarks>
            MS1TotalSignalAfter = 86

            // ReSharper restore InconsistentNaming
        }

        /// <summary>
        /// These columns correspond to tab-delimited _ScanInfo.txt files
        /// </summary>
        private enum ScanInfoFileColumns
        {
            /// <summary>
            /// Scan number
            /// </summary>
            Scan = 0,

            /// <summary>
            /// Scan
            /// </summary>
            ScanTime = 1,

            /// <summary>
            /// Charge
            /// </summary>
            ScanType = 2
        }

        /// <summary>
        /// Keys in this dictionary are modification names, values are modification masses
        /// </summary>
        private readonly Dictionary<string, double> mModificationMassByName;

        /// <summary>
        /// This holds the names of datasets for which a _ScanInfo.txt file was not loaded
        /// </summary>
        private readonly SortedSet<string> mMissingScanInfoDatasets;

        private readonly PeptideCleavageStateCalculator mPeptideCleavageStateCalculator;

        /// <summary>
        /// Add modifications to a peptide read from the DIA-NN synopsis file
        /// Next, compute the monoisotopic mass
        /// </summary>
        /// <param name="searchResult">Search result</param>
        /// <param name="updateModOccurrenceCounts">When true, update mod occurrence counts</param>
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
        /// <param name="searchResult">Search result</param>
        /// <param name="updateModOccurrenceCounts">When true, update mod occurrence counts</param>
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

        private bool AddPrefixAndSuffixResiduesUsingFASTA(
            List<DiaNNSearchResult> filteredSearchResults,
            string outputDirectoryPath)
        {
            const bool mapFileIncludesPrefixAndSuffixColumns = true;
            const bool ignorePeptideToProteinMapErrors = true;

            const double maximumAllowableMatchErrorPercentThreshold = 0.1;
            const double matchErrorPercentWarningThreshold = 0;

            var peptidesWithPrefixAndSuffix = 0;

            var uniquePeptides = new Dictionary<string, string>();

            foreach (var item in filteredSearchResults)
            {
                GetCleanSequence(item.Sequence, out var prefix, out var suffix);

                var peptideHasPrefixAndSuffix = !string.IsNullOrWhiteSpace(prefix) && !string.IsNullOrWhiteSpace(suffix);

                if (peptideHasPrefixAndSuffix)
                {
                    peptidesWithPrefixAndSuffix++;
                }

                if (!uniquePeptides.ContainsKey(item.Sequence))
                {
                    uniquePeptides.Add(item.Sequence, peptideHasPrefixAndSuffix ? item.Sequence : string.Empty);
                }
            }

            if (peptidesWithPrefixAndSuffix == filteredSearchResults.Count)
            {
                // All the peptides already have prefix and suffix residues
                OnStatusEvent(
                    "All the peptides loaded from the DIA-NN results file already have prefix and suffix letters; " +
                    "skipping peptide to protein mapping prior to creating the synopsis file");

                return true;
            }

            if (string.IsNullOrWhiteSpace(Options.FastaFilePath))
            {
                OnWarningEvent("Cannot add prefix and suffix residues to peptides identified by DIA-NN since a FASTA file was not specified when starting PHRP");
                return true;
            }

            var peptideFilePath = Path.Combine(outputDirectoryPath, "Temp_DiaNN_peptides.txt");
            var peptideToProteinMapFilePath = Path.Combine(outputDirectoryPath, "Temp_DiaNN_peptide_to_protein_map.txt");

            using (var writer = new StreamWriter(new FileStream(peptideFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
            {
                // Header line
                writer.WriteLine("Peptide");

                foreach (var item in uniquePeptides)
                {
                    writer.WriteLine(item.Key);
                }
            }

            var sourceDataFiles = new List<FileInfo>
            {
                new(peptideFilePath)
            };

            var success = CreatePepToProteinMapFile(
                sourceDataFiles,
                peptideToProteinMapFilePath,
                mapFileIncludesPrefixAndSuffixColumns,
                ignorePeptideToProteinMapErrors,
                maximumAllowableMatchErrorPercentThreshold,
                matchErrorPercentWarningThreshold,
                clsPeptideToProteinMapEngine.PeptideInputFileFormatConstants.TabDelimitedText);

            if (!success)
                return false;

            var trypticPeptide = new Regex(@"^([KR]\.[^P].+[KR]\.[A-O,Q-Z]|-\..+[KR]\.[A-O,Q-Z]|[KR]\.[^P].+\.-|-\..+\.-)$", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            using (var reader = new StreamReader(new FileStream(peptideToProteinMapFilePath, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                var headerMap = new Dictionary<string, int>();
                var peptideColumnIndex = -1;
                var prefixColumnIndex = -1;
                var suffixColumnIndex = -1;

                while (!reader.EndOfStream)
                {
                    var dataLine = reader.ReadLine();

                    if (string.IsNullOrWhiteSpace(dataLine))
                        continue;

                    var lineParts = dataLine.Split('\t');

                    if (headerMap.Count == 0)
                    {
                        for (var i = 0; i < lineParts.Length; i++)
                        {
                            headerMap.Add(lineParts[i], i);
                        }

                        if (!headerMap.TryGetValue(clsProteinCoverageSummarizer.PROTEIN_TO_PEPTIDE_MAP_FILE_PEPTIDE_COLUMN, out peptideColumnIndex))
                        {
                            SetErrorMessage(string.Format(
                                "The peptide to protein map file does not have column {0} in the header line: {1}",
                                clsProteinCoverageSummarizer.PROTEIN_TO_PEPTIDE_MAP_FILE_PEPTIDE_COLUMN, peptideToProteinMapFilePath));

                            return false;
                        }

                        if (!headerMap.TryGetValue(clsProteinCoverageSummarizer.PROTEIN_TO_PEPTIDE_MAP_FILE_PREFIX_RESIDUE_COLUMN, out prefixColumnIndex))
                        {
                            SetErrorMessage(string.Format(
                                "The peptide to protein map file does not have column {0} in the header line: {1}",
                                clsProteinCoverageSummarizer.PROTEIN_TO_PEPTIDE_MAP_FILE_PREFIX_RESIDUE_COLUMN, peptideToProteinMapFilePath));

                            return false;
                        }

                        if (!headerMap.TryGetValue(clsProteinCoverageSummarizer.PROTEIN_TO_PEPTIDE_MAP_FILE_SUFFIX_RESIDUE_COLUMN, out suffixColumnIndex))
                        {
                            SetErrorMessage(string.Format(
                                "The peptide to protein map file does not have column {0} in the header line: {1}",
                                clsProteinCoverageSummarizer.PROTEIN_TO_PEPTIDE_MAP_FILE_SUFFIX_RESIDUE_COLUMN, peptideToProteinMapFilePath));

                            return false;
                        }

                        continue;
                    }

                    var peptide = lineParts[peptideColumnIndex];
                    var prefixResidue = lineParts[prefixColumnIndex];
                    var suffixResidue = lineParts[suffixColumnIndex];

                    if (!uniquePeptides.TryGetValue(peptide, out var existingContextInfo))
                        continue;

                    if (string.IsNullOrWhiteSpace(prefixResidue) && string.IsNullOrWhiteSpace(suffixResidue))
                        continue;

                    var newPeptideWithContext = string.Format("{0}.{1}.{2}", prefixResidue, peptide, suffixResidue);

                    if (string.IsNullOrWhiteSpace(existingContextInfo))
                    {
                        uniquePeptides[peptide] = newPeptideWithContext;
                        continue;
                    }

                    var existingIsTryptic = trypticPeptide.IsMatch(existingContextInfo);

                    if (!existingIsTryptic && trypticPeptide.IsMatch(newPeptideWithContext))
                    {
                        uniquePeptides[peptide] = newPeptideWithContext;
                    }
                }
            }

            foreach (var item in filteredSearchResults)
            {
                if (!uniquePeptides.TryGetValue(item.Sequence, out var peptideWithContext))
                {
                    throw new KeyNotFoundException(string.Format(
                        "uniquePeptides dictionary is missing peptide {0}; programming bug", item.Sequence));
                }

                if (string.IsNullOrWhiteSpace(peptideWithContext))
                {
                    OnWarningEvent("Peptide loaded from the peptide to protein map file does not have a mapping for peptide {0}", item.Sequence);
                    continue;
                }

                if (peptideWithContext.IndexOf(item.Sequence, StringComparison.Ordinal) < 0)
                {
                    OnWarningEvent(
                        "Peptide loaded from the peptide to protein map file does not include the original peptide sequence: {0} not in {1}",
                        item.Sequence, peptideWithContext);

                    continue;
                }

                item.Sequence = peptideWithContext;
            }

            // Delete the temp files
            DeleteFileIgnoreErrors(peptideFilePath);
            DeleteFileIgnoreErrors(peptideToProteinMapFilePath);

            return true;
        }

        /// <summary>
        /// Compute the monoisotopic MH value using the calculated monoisotopic mass
        /// </summary>
        /// <remarks>This is (M+H)+ when the charge carrier is a proton</remarks>
        /// <param name="searchResult">Search result</param>
        /// <returns>(M+H)+, as a string</returns>
        private string ComputeMH(ToolResultsBaseClass searchResult)
        {
            return PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(searchResult.CalculatedMonoMassValue, 0), 6);
        }

        /// <summary>
        /// Computes the total of all modifications defined for the sequence
        /// </summary>
        /// <param name="searchResult">Search result</param>
        /// <param name="modList">Output: list of modifications</param>
        private double ComputeTotalModMass(DiaNNSearchResult searchResult, out List<MSFraggerResultsProcessor.MSFraggerModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(searchResult.ModifiedSequence))
            {
                modList = new List<MSFraggerResultsProcessor.MSFraggerModInfo>();
                return 0;
            }

            modList = GetPeptideModifications(searchResult);

            return modList.Sum(modEntry => modEntry.ModMass);
        }

        private void ComputeTrypticStateAndMissedCleavages(ToolResultsBaseClass searchResult, string peptideWithPrefixAndSuffix)
        {
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
        }

        private string ConvertModListToMSFraggerNotation(List<MSFraggerResultsProcessor.MSFraggerModInfo> modList)
        {
            var modificationList = new StringBuilder();

            foreach (var modItem in modList)
            {
                if (modificationList.Length > 0)
                    modificationList.Append(",");

                // ReSharper disable once SwitchStatementHandlesSomeKnownEnumValuesWithDefault
                switch (modItem.TerminusState)
                {
                    case AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus:
                        modificationList.AppendFormat("{0}({1:0.0###})", MSFraggerResultsProcessor.MOD_POSITION_NAME_N_TERMINUS, modItem.ModMass);
                        break;

                    case AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus:
                        modificationList.AppendFormat("{0}({1:0.0###})", MSFraggerResultsProcessor.MOD_POSITION_NAME_C_TERMINUS, modItem.ModMass);
                        break;

                    case AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus:
                    case AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus:
                    case AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus:
                        throw new Exception(string.Format(
                            "Unexpected terminus state for mod item: {0} ({1})", modItem.TerminusState, (int)modItem.TerminusState));

                    default:
                        // Includes: AminoAcidModInfo.ResidueTerminusState.None:
                        modificationList.AppendFormat("{0}{1}({2:0.0###})", modItem.ResidueLocInPeptide, modItem.ResidueSymbol, modItem.ModMass);
                        break;
                }
            }

            return modificationList.ToString();
        }

        /// <summary>
        /// This routine creates a synopsis file using search results read from the DIA-NN report.tsv file
        /// The synopsis file includes every result with a probability above a set threshold
        /// </summary>
        /// <param name="inputFilePath">DIA-NN results file (either .tsv or .parquet)</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
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

                var resultsLoaded = ReadDiaNNResults(
                    inputFile,
                    errorMessages,
                    out var filteredSearchResults,
                    out var datasetNameToBaseNameMap,
                    out var usingParquetFile);

                if (!resultsLoaded)
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

                // Write the results to a tab-delimited text file then use the peptide to protein mapper to determine the prefix and suffix letters for each peptide
                var residuesAdded = AddPrefixAndSuffixResiduesUsingFASTA(filteredSearchResults, outputDirectoryPath);

                if (!residuesAdded)
                {
                    baseName = string.Empty;
                    synOutputFilePath = string.Empty;

                    return false;
                }

                // Now that prefix and suffix residues are defined, update the cleavage state and missed cleavage values
                foreach (var searchResult in filteredSearchResults)
                {
                    ComputeTrypticStateAndMissedCleavages(searchResult, searchResult.Sequence);
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
                WriteSynFHTFileHeader(writer, usingParquetFile, errorMessages);

                // Write the search results to disk
                WriteFilteredSearchResults(datasetNameToBaseNameMap, writer, filteredSearchResults, usingParquetFile, errorMessages);

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
        /// Read mod info from the DIA-NN parameter file
        /// </summary>
        /// <remarks>The DMS-based parameter file for DIA-NN uses the same formatting as MS-GF+</remarks>
        /// <param name="diannParamFilePath">DIA-NN parameter file path</param>
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
        /// Get the list of static and dynamic modifications for a peptide read from the DIA-NN report.tsv file
        /// </summary>
        /// <param name="cleanSequence">Peptide sequence without modification symbols</param>
        /// <param name="peptideWithModifications">Peptide sequence, with embedded modification names (DIA-NN style)</param>
        /// <returns>List of modifications</returns>
        private List<MSFraggerResultsProcessor.MSFraggerModInfo> ExtractModsFromModifiedSequence(string cleanSequence, string peptideWithModifications)
        {
            // ReSharper disable CommentTypo

            // peptideWithModifications should have the peptide sequence, including modification names; examples:
            //   AAALEAM(UniMod:35)K
            //   C(UniMod:4)TSYNIPC(UniMod:4)TSDM(UniMod:35)AK
            //   DIHFM(UniMod:35)PC(UniMod:4)SGLTGANLK

            // ReSharper enable CommentTypo

            var mods = new List<MSFraggerResultsProcessor.MSFraggerModInfo>();

            var finalResidueLoc = cleanSequence.Length;

            var currentIndex = -1;
            var maxIndexToExamine = peptideWithModifications.Length - 2;

            var residueNumber = 0;
            var currentResidue = '-';

            while (currentIndex < maxIndexToExamine)
            {
                currentIndex++;

                if (char.IsLetter(peptideWithModifications[currentIndex]))
                {
                    currentResidue = peptideWithModifications[currentIndex];
                    residueNumber++;
                }

                int addon;

                if (currentIndex == 0 && peptideWithModifications[currentIndex].Equals('('))
                {
                    // N-terminal peptide mod

                    addon = 1;
                }
                else if (!peptideWithModifications[currentIndex + 1].Equals('('))
                {
                    continue;
                }
                else
                {
                    addon = 2;
                }

                var closingParenthesisIndex = peptideWithModifications.IndexOf(')', currentIndex + 1);

                if (closingParenthesisIndex < 0)
                {
                    OnWarningEvent("Mismatched parenthesis after index {0} in {1}", currentIndex, peptideWithModifications);
                    return mods;
                }

                if (closingParenthesisIndex == currentIndex + addon)
                {
                    OnWarningEvent("Empty modification name after index {0} in {1}", currentIndex, peptideWithModifications);
                    continue;
                }

                // Parse out the modification name
                var modificationName = peptideWithModifications.Substring(currentIndex + addon, closingParenthesisIndex - currentIndex - addon);

                var modMass = GetModificationMass(modificationName);

                var currentMod = new MSFraggerResultsProcessor.MSFraggerModInfo
                {
                    ResidueSymbol = currentResidue,
                    ResidueLocInPeptide = residueNumber,
                    ModMass = modMass
                };

                if (currentMod.ResidueLocInPeptide < 1)
                {
                    currentMod.TerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                }
                else if (currentMod.ResidueLocInPeptide > finalResidueLoc)
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
                    throw new Exception(string.Format(
                        "Residue number tracked by residueNumber is larger than the final residue location: {0} vs. {1} for {2}",
                        residueNumber, finalResidueLoc, peptideWithModifications));
                }

                if (residueNumber > 0 && currentResidue != cleanSequence[residueNumber - 1])
                {
                    throw new Exception(string.Format(
                        "Residue {0} tracked by currentResidue does not match the residue in the clean sequence: {1} vs. {2} for {3}",
                        residueNumber, currentResidue, cleanSequence[residueNumber - 1], peptideWithModifications));
                }

                mods.Add(currentMod);

                currentIndex = closingParenthesisIndex;
            }

            return mods;
        }

        private string GetDatasetNameFromFileName(string fileName, string fileNameSuffix)
        {
            if (fileName.EndsWith(fileNameSuffix, StringComparison.OrdinalIgnoreCase))
            {
                return fileName.Substring(0, fileName.Length - fileNameSuffix.Length);
            }

            OnWarningEvent("Filename does not end with {0}: {1}", fileNameSuffix, fileName);

            return Path.GetFileNameWithoutExtension(fileName);
        }

        private double GetModificationMass(string modificationName)
        {
            if (mModificationMassByName.TryGetValue(modificationName, out var modMass))
                return modMass;

            return 0;
        }

        /// <summary>
        /// Populate list outputData with the header names for the .tsv file created from the .parquet file
        /// </summary>
        /// <param name="columnNamesByEnum">Dictionary where keys are enum DiaNNReportFileColumns and values are the corresponding column name</param>
        private List<string> GetParquetTsvHeaderNames(IReadOnlyDictionary<DiaNNReportFileColumns, string> columnNamesByEnum)
        {
            return new List<string>
            {
                columnNamesByEnum[DiaNNReportFileColumns.RunIndex],
                columnNamesByEnum[DiaNNReportFileColumns.DatasetName],
                columnNamesByEnum[DiaNNReportFileColumns.Channel],
                columnNamesByEnum[DiaNNReportFileColumns.PrecursorId],
                columnNamesByEnum[DiaNNReportFileColumns.ModifiedSequence],
                columnNamesByEnum[DiaNNReportFileColumns.StrippedSequence],
                columnNamesByEnum[DiaNNReportFileColumns.PrecursorCharge],
                columnNamesByEnum[DiaNNReportFileColumns.PrecursorLibIndex],
                columnNamesByEnum[DiaNNReportFileColumns.Decoy],
                columnNamesByEnum[DiaNNReportFileColumns.Proteotypic],
                columnNamesByEnum[DiaNNReportFileColumns.PrecursorMz],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinIDs],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinGroup],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinNames],
                columnNamesByEnum[DiaNNReportFileColumns.GeneNames],
                columnNamesByEnum[DiaNNReportFileColumns.RT],
                columnNamesByEnum[DiaNNReportFileColumns.IndexedRT],
                columnNamesByEnum[DiaNNReportFileColumns.PredictedRT],
                columnNamesByEnum[DiaNNReportFileColumns.PredictedIndexedRT],
                columnNamesByEnum[DiaNNReportFileColumns.IonMobility],
                columnNamesByEnum[DiaNNReportFileColumns.IndexedIonMobility],
                columnNamesByEnum[DiaNNReportFileColumns.PredictedIonMobility],
                columnNamesByEnum[DiaNNReportFileColumns.PredictedIndexedIonMobility],
                columnNamesByEnum[DiaNNReportFileColumns.PrecursorQuantity],
                columnNamesByEnum[DiaNNReportFileColumns.PrecursorNormalized],
                columnNamesByEnum[DiaNNReportFileColumns.MS1Area],
                columnNamesByEnum[DiaNNReportFileColumns.MS1Normalised],
                columnNamesByEnum[DiaNNReportFileColumns.MS1ApexArea],
                columnNamesByEnum[DiaNNReportFileColumns.MS1ApexMzDelta],
                columnNamesByEnum[DiaNNReportFileColumns.NormalisationFactor],
                columnNamesByEnum[DiaNNReportFileColumns.QuantityQuality],
                columnNamesByEnum[DiaNNReportFileColumns.EmpiricalQuality],
                columnNamesByEnum[DiaNNReportFileColumns.NormalisationNoise],
                columnNamesByEnum[DiaNNReportFileColumns.MS1ProfileCorr],
                columnNamesByEnum[DiaNNReportFileColumns.Evidence],
                columnNamesByEnum[DiaNNReportFileColumns.MassEvidence],
                columnNamesByEnum[DiaNNReportFileColumns.ChannelEvidence],
                columnNamesByEnum[DiaNNReportFileColumns.MS1TotalSignalBefore],
                columnNamesByEnum[DiaNNReportFileColumns.MS1TotalSignalAfter],
                columnNamesByEnum[DiaNNReportFileColumns.RTStart],
                columnNamesByEnum[DiaNNReportFileColumns.RTStop],
                columnNamesByEnum[DiaNNReportFileColumns.FWHM],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinGroupTopN],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinGroupMaxLFQ],
                columnNamesByEnum[DiaNNReportFileColumns.GenesTopN],
                columnNamesByEnum[DiaNNReportFileColumns.GenesMaxLFQ],
                columnNamesByEnum[DiaNNReportFileColumns.GenesMaxLFQUnique],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinGroupMaxLFQQuality],
                columnNamesByEnum[DiaNNReportFileColumns.GenesMaxLFQQuality],
                columnNamesByEnum[DiaNNReportFileColumns.GenesMaxLFQUniqueQuality],
                columnNamesByEnum[DiaNNReportFileColumns.QValue],
                columnNamesByEnum[DiaNNReportFileColumns.PEP],
                columnNamesByEnum[DiaNNReportFileColumns.GlobalQValue],
                columnNamesByEnum[DiaNNReportFileColumns.LibQValue],
                columnNamesByEnum[DiaNNReportFileColumns.PeptidoformQValue],
                columnNamesByEnum[DiaNNReportFileColumns.GlobalPeptidoformQValue],
                columnNamesByEnum[DiaNNReportFileColumns.LibPeptidoformQValue],
                columnNamesByEnum[DiaNNReportFileColumns.PtmSiteConfidence],
                columnNamesByEnum[DiaNNReportFileColumns.SiteOccupancyProbabilities],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinSites],
                columnNamesByEnum[DiaNNReportFileColumns.LibPtmSiteConfidence],
                columnNamesByEnum[DiaNNReportFileColumns.TranslatedQValue],
                columnNamesByEnum[DiaNNReportFileColumns.ChannelQValue],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinGroupQValue],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinGroupPEP],
                columnNamesByEnum[DiaNNReportFileColumns.GeneGroupQValue],
                columnNamesByEnum[DiaNNReportFileColumns.ProteinQValue],
                columnNamesByEnum[DiaNNReportFileColumns.GlobalProteinGroupQValue],
                columnNamesByEnum[DiaNNReportFileColumns.LibProteinGroupQValue]
            };
        }

        /// <summary>
        /// Get the list of static and dynamic modifications for a peptide read from a DIA-NN synopsis file
        /// </summary>
        /// <param name="searchResult">Search result</param>
        /// <returns>List of modifications</returns>
        private List<MSFraggerResultsProcessor.MSFraggerModInfo> GetPeptideModifications(DiaNNResults searchResult)
        {
            return GetPeptideModifications(searchResult.PeptideCleanSequence, searchResult.Modifications);
        }

        /// <summary>
        /// Get the list of static and dynamic modifications for a peptide read from the DIA-NN report.tsv file
        /// </summary>
        /// <param name="searchResult">Search result</param>
        /// <returns>List of modifications</returns>
        private List<MSFraggerResultsProcessor.MSFraggerModInfo> GetPeptideModifications(DiaNNSearchResult searchResult)
        {
            return ExtractModsFromModifiedSequence(searchResult.Sequence, searchResult.ModifiedSequence);
        }

        /// <summary>
        /// Get the list of static and dynamic modifications for a peptide with MSFragger style modification descriptions
        /// </summary>
        /// <param name="cleanSequence">Peptide sequence without modification symbols</param>
        /// <param name="modificationList">Comma separated list of modification info, e.g. 1M(15.9949), 5C(57.0215)</param>
        /// <returns>List of modifications</returns>
        internal List<MSFraggerResultsProcessor.MSFraggerModInfo> GetPeptideModifications(string cleanSequence, string modificationList)
        {
            // modificationList should have both the static and dynamic mods,
            // as a comma separated list of residue number, residue symbol, and mod mass; examples:
            //   15M(15.9949)
            //   1M(15.9949), 5C(57.0215)
            //   N-term(42.0106)

            var mods = new List<MSFraggerResultsProcessor.MSFraggerModInfo>();
            var finalResidueLoc = cleanSequence.Length;

            foreach (var modEntry in modificationList.Split(','))
            {
                if (MSFraggerResultsProcessor.ParseModificationDescription(cleanSequence, finalResidueLoc, modEntry.Trim(), out var modInfo, out var errorMessage))
                {
                    mods.Add(modInfo);
                    continue;
                }

                ReportError(errorMessage);
            }

            return mods;
        }

        /// <summary>
        /// Read the static and dynamic modifications from the DIA-NN parameter file
        /// </summary>
        /// <param name="diannParamFilePath">DIA-NN parameter file path</param>
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
        /// Determine scan number that corresponds to the given elution time, for the given dataset
        /// </summary>
        /// <param name="scanInfoByDataset">Dictionary where keys are dataset name and values are a class mapping elution times to scan numbers</param>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="elutionTime">Elution time</param>
        /// <returns>Closest scan number</returns>
        private int LookupScanNumber(
            IReadOnlyDictionary<string, BinarySearchFindNearest> scanInfoByDataset,
            string datasetName,
            string elutionTime)
        {
            if (!scanInfoByDataset.TryGetValue(datasetName, out var elutionTimeToScanMap))
            {
                // ReSharper disable once InvertIf
                if (!mMissingScanInfoDatasets.Contains(datasetName))
                {
                    OnWarningEvent("Dataset not found in scanInfoByDataset: {0}", datasetName);
                    mMissingScanInfoDatasets.Add(datasetName);
                }

                return 0;
            }

            if (!double.TryParse(elutionTime, out var elutionTimeValue))
            {
                OnWarningEvent("Non-numeric elution time for dataset {0}: {1}", datasetName, elutionTime);
            }

            var scanNumber = elutionTimeToScanMap.GetYForX(elutionTimeValue);

            return (int)Math.Round(scanNumber);
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions">If true, reset the mass correction tags and modification definitions</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
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

                    var columnMapping = new Dictionary<DiaNNSynFileColumns, int>();

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
        /// Parse a DIA-NN results line while reading the report.tsv file
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="searchResult">Search result</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <param name="lineNumber">Line number in the input file (used for error reporting)</param>
        /// <param name="datasetNameToBaseNameMap">Keys are full dataset names, values are abbreviated dataset names</param>
        /// <param name="baseNameToFullDatasetNameMap">Keys are abbreviated dataset names, values are the full dataset names</param>
        /// <param name="currentDatasetFile">Current dataset file path; updated by this method if the Spectrum File path is not an empty string</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNResultsFileEntry(
            string lineIn,
            out DiaNNSearchResult searchResult,
            ICollection<string> errorMessages,
            IDictionary<DiaNNReportFileColumns, int> columnMapping,
            int lineNumber,
            IDictionary<string, string> datasetNameToBaseNameMap,
            IDictionary<string, string> baseNameToFullDatasetNameMap,
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

                // Read the abbreviated dataset name
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.DatasetName], out searchResult.DatasetName);

                if (!string.IsNullOrWhiteSpace(searchResult.DatasetName))
                {
                    if (updateDatasetNameMapDictionary)
                    {
                        searchResult.FullDatasetName = Path.GetFileNameWithoutExtension(spectrumFile);

                        if (!string.IsNullOrEmpty(searchResult.FullDatasetName))
                        {
                            if (!datasetNameToBaseNameMap.ContainsKey(searchResult.FullDatasetName))
                            {
                                datasetNameToBaseNameMap.Add(searchResult.FullDatasetName, searchResult.DatasetName);
                            }

                            if (!baseNameToFullDatasetNameMap.ContainsKey(searchResult.DatasetName))
                            {
                                baseNameToFullDatasetNameMap.Add(searchResult.DatasetName, searchResult.FullDatasetName);
                            }
                        }
                    }
                    else if (!baseNameToFullDatasetNameMap.TryGetValue(searchResult.DatasetName, out searchResult.FullDatasetName))
                    {
                        // This code should never be reached
                        OnWarningEvent("Abbreviated dataset name not found in baseNameToFullDatasetNameMap; this is unexpected");

                        searchResult.FullDatasetName = Path.GetFileNameWithoutExtension(spectrumFile);

                        if (!string.IsNullOrEmpty(searchResult.FullDatasetName))
                        {
                            baseNameToFullDatasetNameMap.Add(searchResult.DatasetName, searchResult.FullDatasetName);
                        }
                    }
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinGroup], out searchResult.ProteinGroup);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.ProteinIDs], out searchResult.ProteinIDs);
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

                // Peptide sequences in DIA-NN report.tsv files do not include prefix or suffix residues (but this class will add them later)
                // DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrevAA], out searchResult.PrefixResidue);
                // DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.NextAA], out searchResult.SuffixResidue);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorId], out searchResult.PrecursorId);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.PrecursorCharge], out searchResult.Charge);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.QValue], out searchResult.QValueDiaNN);
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
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MS1ProfileCorr], out searchResult.MS1ProfileCorrelation);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MS1Area], out searchResult.MS1Area);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Evidence], out searchResult.Evidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.SpectrumSimilarity], out searchResult.SpectrumSimilarity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.Averagine], out searchResult.Averagine);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.MassEvidence], out searchResult.MassEvidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNReportFileColumns.CScore], out searchResult.CScore);
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

                searchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(searchResult.CalculatedMonoMassValue, 0, searchResult.ChargeNum), 6);

                return ParseDiaNNSearchResult(searchResult, errorMessages, lineNumber, true);
            }
            catch (Exception ex)
            {
                // Error parsing this row from the DIA-NN results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing DIA-NN results in ParseDiaNNResultsFileEntry, line {0}: {1}", lineNumber, ex.Message));
                }

                return false;
            }
        }

        /// <summary>
        /// Parse a DIA-NN results line while reading the report.tsv file
        /// </summary>
        /// <param name="searchResult">Search result</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="lineOrRowNumber">Line number or row number in the input file (used for error reporting)</param>
        /// <param name="readingTsvFile">True if reading a .tsv file, false if reading a .parquet file</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNSearchResult(
            DiaNNSearchResult searchResult,
            ICollection<string> errorMessages,
            int lineOrRowNumber,
            bool readingTsvFile)
        {
            var lineOrRowDescription = readingTsvFile ? "line" : "row";

            try
            {
                var proteinIdList = searchResult.ProteinIDs.Split(';').ToList();

                if (proteinIdList.Count > 0)
                    searchResult.Protein = proteinIdList[0];

                if (string.IsNullOrWhiteSpace(searchResult.Sequence))
                {
                    // Example messages:
                    // - Missing sequence on line 5
                    // - Missing sequence on row 5
                    ReportError(string.Format("Missing sequence on {0} {1}", lineOrRowDescription, lineOrRowNumber), true);
                }

                searchResult.Length = searchResult.Sequence.Length;

                searchResult.ChargeNum = (short)StringUtilities.CIntSafe(searchResult.Charge, 0);

                if (!string.IsNullOrWhiteSpace(searchResult.QValueDiaNN))
                {
                    double.TryParse(searchResult.QValueDiaNN, out searchResult.QValue);
                }

                if (!string.IsNullOrWhiteSpace(searchResult.CScore))
                {
                    double.TryParse(searchResult.CScore, out searchResult.ConfidenceScore);
                }

                searchResult.Reverse = IsReversedProtein(searchResult.Protein);

                // Parse the modification list to determine the total mod mass
                // In addition, obtain the list of static and dynamic modifications in Modified.Sequence
                var totalModMass = ComputeTotalModMass(searchResult, out var modList);

                searchResult.ModificationList = ConvertModListToMSFraggerNotation(modList);

                // Compute monoisotopic mass of the peptide
                searchResult.CalculatedMonoMassPHRP = ComputePeptideMassForCleanSequence(searchResult.Sequence, totalModMass);

                searchResult.CalculatedMonoMass = PRISM.StringUtilities.DblToString(searchResult.CalculatedMonoMassPHRP, 6);
                searchResult.CalculatedMonoMassValue = searchResult.CalculatedMonoMassPHRP;

                searchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(searchResult.CalculatedMonoMassValue, 0, searchResult.ChargeNum), 6);

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                searchResult.MH = ComputeMH(searchResult);

                // Assume the peptide is fully tryptic (preceded by a K or R), then compute cleavage state and missed cleavages

                var peptideWithPrefixAndSuffix = string.Format("K.{0}.A", searchResult.Sequence);

                ComputeTrypticStateAndMissedCleavages(searchResult, peptideWithPrefixAndSuffix);

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this line/row from the DIA-NN results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing DIA-NN results in ParseDiaNNSearchResult, {0} {1}: {2}", lineOrRowDescription, lineOrRowNumber, ex.Message));
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the DIA-NN report.tsv file header line, populating columnMapping
        /// </summary>
        /// <param name="headerLine">Tab-delimited list of header names</param>
        /// <param name="columnMapping">Column mapping, where keys are an enum in DiaNNReportFileColumns and values are the zero-based index</param>
        /// <param name="columnNamesByEnum">Output: Dictionary where keys are enum DiaNNReportFileColumns and values are the corresponding column name</param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseDiaNNResultsFileHeaderLine(
            string headerLine,
            IDictionary<DiaNNReportFileColumns, int> columnMapping,
            out Dictionary<DiaNNReportFileColumns, string> columnNamesByEnum)
        {
            columnNamesByEnum = new Dictionary<DiaNNReportFileColumns, string>();

            // ReSharper disable StringLiteralTypo

            // Columns in report.tsv and/or report.parquet files
            var columnNames = new SortedDictionary<string, DiaNNReportFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "Run.Index", DiaNNReportFileColumns.RunIndex },                                       // Only in .parquet
                { "File.Name", DiaNNReportFileColumns.DatasetFile },
                { "Run", DiaNNReportFileColumns.DatasetName },
                { "Channel", DiaNNReportFileColumns.Channel },                                          // Only in .parquet
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
                { "Precursor.Lib.Index", DiaNNReportFileColumns.PrecursorLibIndex },                    // Only in .parquet
                { "Decoy", DiaNNReportFileColumns.Decoy },                                              // Only in .parquet
                { "Precursor.Mz", DiaNNReportFileColumns.PrecursorMz },                                 // Only in .parquet
                { "PG.MaxLFQ.Quality", DiaNNReportFileColumns.ProteinGroupMaxLFQQuality },              // Only in .parquet
                { "Genes.MaxLFQ.Quality", DiaNNReportFileColumns.GenesMaxLFQQuality },                  // Only in .parquet
                { "Genes.MaxLFQ.Unique.Quality", DiaNNReportFileColumns.GenesMaxLFQUniqueQuality },          // Only in .parquet
                { "Q.Value", DiaNNReportFileColumns.QValue },
                { "PEP", DiaNNReportFileColumns.PEP },
                { "Global.Q.Value", DiaNNReportFileColumns.GlobalQValue },
                { "Peptidoform.Q.Value", DiaNNReportFileColumns.PeptidoformQValue },                   // Only in .parquet
                { "Global.Peptidoform.Q.Value", DiaNNReportFileColumns.GlobalPeptidoformQValue },      // Only in .parquet
                { "Lib.Peptidoform.Q.Value", DiaNNReportFileColumns.LibPeptidoformQValue },            // Only in .parquet
                { "PTM.Site.Confidence", DiaNNReportFileColumns.PtmSiteConfidence },                   // Only in .parquet
                { "Site.Occupancy.Probabilities", DiaNNReportFileColumns.SiteOccupancyProbabilities }, // Only in .parquet
                { "Protein.Sites", DiaNNReportFileColumns.ProteinSites },                              // Only in .parquet
                { "Protein.Q.Value", DiaNNReportFileColumns.ProteinQValue  },
                { "Lib.PTM.Site.Confidence", DiaNNReportFileColumns.LibPtmSiteConfidence },             // Only in .parquet
                { "Channel.Q.Value", DiaNNReportFileColumns.ChannelQValue },                            // Only in .parquet
                { "PG.Q.Value", DiaNNReportFileColumns.ProteinGroupQValue },
                { "PG.PEP", DiaNNReportFileColumns.ProteinGroupPEP },                                   // Only in .parquet
                { "Global.PG.Q.Value", DiaNNReportFileColumns.GlobalProteinGroupQValue },
                { "GG.Q.Value", DiaNNReportFileColumns.GeneGroupQValue },
                { "Translated.Q.Value", DiaNNReportFileColumns.TranslatedQValue },
                { "Proteotypic", DiaNNReportFileColumns.Proteotypic },
                { "Precursor.Quantity", DiaNNReportFileColumns.PrecursorQuantity },
                { "Precursor.Normalised", DiaNNReportFileColumns.PrecursorNormalized },
                { "Precursor.Translated", DiaNNReportFileColumns.PrecursorTranslated },
                { "Translated.Quality", DiaNNReportFileColumns.TranslatedQuality },
                { "Ms1.Area", DiaNNReportFileColumns.MS1Area },
                { "Ms1.Normalised", DiaNNReportFileColumns.MS1Normalised },             // Only in .parquet
                { "Ms1.Apex.Area", DiaNNReportFileColumns.MS1ApexArea },                // Only in .parquet
                { "Ms1.Apex.Mz.Delta", DiaNNReportFileColumns.MS1ApexMzDelta },         // Only in .parquet
                { "Normalisation.Factor", DiaNNReportFileColumns.NormalisationFactor }, // Only in .parquet
                { "Ms1.Translated", DiaNNReportFileColumns.MS1Translated },
                { "Quantity.Quality", DiaNNReportFileColumns.QuantityQuality },
                { "Empirical.Quality", DiaNNReportFileColumns.EmpiricalQuality },       // Only in .parquet
                { "Normalisation.Noise", DiaNNReportFileColumns.NormalisationNoise },   // Only in .parquet
                { "RT", DiaNNReportFileColumns.RT },
                { "RT.Start", DiaNNReportFileColumns.RTStart },
                { "RT.Stop", DiaNNReportFileColumns.RTStop },
                { "FWHM", DiaNNReportFileColumns.FWHM },                                // Only in .parquet
                { "PG.TopN", DiaNNReportFileColumns.ProteinGroupTopN },                 // Only in .parquet
                { "Genes.TopN", DiaNNReportFileColumns.GenesTopN },                     // Only in .parquet
                { "iRT", DiaNNReportFileColumns.IndexedRT },
                { "Predicted.RT", DiaNNReportFileColumns.PredictedRT },
                { "Predicted.iRT", DiaNNReportFileColumns.PredictedIndexedRT },
                { "First.Protein.Description", DiaNNReportFileColumns.FirstProteinDescription },
                { "Lib.Q.Value", DiaNNReportFileColumns.LibQValue },
                { "Lib.PG.Q.Value", DiaNNReportFileColumns.LibProteinGroupQValue },
                { "Ms1.Profile.Corr", DiaNNReportFileColumns.MS1ProfileCorr },
                { "Evidence", DiaNNReportFileColumns.Evidence },
                { "Spectrum.Similarity", DiaNNReportFileColumns.SpectrumSimilarity },
                { "Averagine", DiaNNReportFileColumns.Averagine },
                { "Mass.Evidence", DiaNNReportFileColumns.MassEvidence },
                { "Channel.Evidence", DiaNNReportFileColumns.ChannelEvidence },                 // Only in .parquet
                { "Ms1.Total.Signal.Before", DiaNNReportFileColumns.MS1TotalSignalBefore },     // Only in .parquet
                { "Ms1.Total.Signal.After", DiaNNReportFileColumns.MS1TotalSignalAfter },       // Only in .parquet
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

                var splitLine = headerLine.Split('\t');

                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }

                // Populate columnNamesByEnum
                foreach (var item in columnNames)
                {
                    columnNamesByEnum.Add(item.Value, item.Key);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the DIA-NN results file", ex);
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a DIA-NN _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseDiaNNSynFileHeaderLine(string lineIn, IDictionary<DiaNNSynFileColumns, int> columnMapping)
        {
            var columnNames = DiaNNSynFileReader.GetColumnHeaderNamesAndIDs(false);

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
                SetErrorMessage("Error parsing the header line in the DIA-NN synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a DIA-NN Synopsis file
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="searchResult">Search result</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="resultsProcessed">Number of results loaded when this method is called</param>
        /// <param name="columnMapping">Column mapping</param>
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
                            "Error reading ResultID value from DIA-NN results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Dataset], out string datasetName);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DatasetID], out int datasetId);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Scan], out string scan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.IonMobility], out string ionMobility);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Charge], out string charge);

                // The dataset name will typically be an abbreviated version of the full dataset name
                searchResult.DatasetName = datasetName;
                searchResult.DatasetID = datasetId;

                searchResult.Scan = scan;
                searchResult.IonMobility = ionMobility;
                searchResult.Charge = charge;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Peptide], out string peptideSequence))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading Peptide sequence value from DIA-NN results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                // searchResult.ModificationList,
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Modifications], out string modifications);
                searchResult.Modifications = modifications;

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PrecursorMZ], out string precursorMz);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MH], out string parentIonMH);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Mass], out string monoisotopicMass);

                searchResult.PrecursorMZ = precursorMz;
                searchResult.ParentIonMH = parentIonMH;
                searchResult.CalculatedMonoMass = monoisotopicMass;

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinGroup], out string proteinGroup);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinIDs], out string proteinIDs);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinNames], out string proteinNames);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GeneNames], out string geneNames);

                searchResult.ProteinGroup = proteinGroup;
                searchResult.ProteinIDs.AddRange(ParseDelimitedList(proteinIDs, ';'));
                searchResult.ProteinNames.AddRange(ParseDelimitedList(proteinNames, ';'));
                searchResult.GeneNames.AddRange(ParseDelimitedList(geneNames, ';'));

                if (searchResult.ProteinIDs.Count > 0)
                {
                    searchResult.ProteinName = searchResult.ProteinIDs[0];
                    searchResult.MultipleProteinCount = (searchResult.ProteinIDs.Count - 1).ToString();
                }
                else
                {
                    searchResult.MultipleProteinCount = "0";
                }

                // ReSharper disable once GrammarMistakeInComment

                // Note that DIA-NN peptides don't actually have mod symbols; that information is tracked via searchResult.Modifications
                // Thus, .PeptideSequenceWithMods will not have any mod symbols

                // Calling this method will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the cleavage state and terminus state

                if (string.IsNullOrWhiteSpace(searchResult.PeptidePreResidues) && string.IsNullOrWhiteSpace(searchResult.PeptidePostResidues))
                {
                    // Peptide sequences identified by DIA-NN do not include prefix or suffix residues
                    // Assume the peptide is fully tryptic (preceded by a K or R), then compute cleavage state and missed cleavages
                    searchResult.ComputePeptideCleavageStateInProtein(searchResult.PeptideCleanSequence, "K", "A");
                }
                else
                {
                    searchResult.ComputePeptideCleavageStateInProtein();
                }

                // Read the remaining data values

                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.NTT], out string numberOfTrypticTermini);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinGroupQuantity], out string proteinGroupQuantity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinGroupNormalized], out string proteinGroupNormalized);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinGroupMaxLFQ], out string proteinGroupMaxLFQ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GenesQuantity], out string genesQuantity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GenesNormalized], out string genesNormalized);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GenesMaxLFQ], out string genesMaxLFQ);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GenesMaxLFQUnique], out string genesMaxLFQUnique);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.QValue], out string qValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PEP], out string pep);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GlobalQValue], out string globalQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinQValue], out string proteinQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ProteinGroupQValue], out string proteinGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GlobalProteinGroupQValue], out string globalProteinGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.GeneGroupQValue], out string geneGroupQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.TranslatedQValue], out string translatedQValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Proteotypic], out string proteotypic);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PrecursorQuantity], out string precursorQuantity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PrecursorNormalized], out string precursorNormalized);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PrecursorTranslated], out string precursorTranslated);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.TranslatedQuality], out string translatedQuality);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MS1Translated], out string ms1Translated);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.QuantityQuality], out string quantityQuality);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ElutionTime], out string elutionTime);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ElutionTimeStart], out string elutionTimeStart);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.ElutionTimeStop], out string elutionTimeStop);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.IndexedRT], out string indexedRT);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.IndexedIonMobility], out string indexedIonMobility);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PredictedRT], out string predictedRT);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.PredictedIndexedRT], out string predictedIndexedRT);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MS1ProfileCorrelation], out string mS1ProfileCorrelation);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MS1Area], out string ms1Area);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Evidence], out string evidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.SpectrumSimilarity], out string spectrumSimilarity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.Averagine], out string averagine);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.MassEvidence], out string massEvidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.CScore], out string cScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DecoyEvidence], out string decoyEvidence);
                DataUtilities.GetColumnValue(splitLine, columnMapping[DiaNNSynFileColumns.DecoyCScore], out string decoyCScore);

                // Store the data
                searchResult.NTT = numberOfTrypticTermini;
                searchResult.ProteinGroupQuantity = proteinGroupQuantity;
                searchResult.ProteinGroupNormalized = proteinGroupNormalized;
                searchResult.ProteinGroupMaxLFQ = proteinGroupMaxLFQ;
                searchResult.GenesQuantity = genesQuantity;
                searchResult.GenesNormalized = genesNormalized;
                searchResult.GenesMaxLFQ = genesMaxLFQ;
                searchResult.GenesMaxLFQUnique = genesMaxLFQUnique;
                searchResult.QValue = qValue;
                searchResult.PEP = pep;
                searchResult.GlobalQValue = globalQValue;
                searchResult.ProteinQValue = proteinQValue;
                searchResult.ProteinGroupQValue = proteinGroupQValue;
                searchResult.GlobalProteinGroupQValue = globalProteinGroupQValue;
                searchResult.GeneGroupQValue = geneGroupQValue;
                searchResult.TranslatedQValue = translatedQValue;
                searchResult.Proteotypic = proteotypic;
                searchResult.PrecursorQuantity = precursorQuantity;
                searchResult.PrecursorNormalized = precursorNormalized;
                searchResult.PrecursorTranslated = precursorTranslated;
                searchResult.TranslatedQuality = translatedQuality;
                searchResult.MS1Translated = ms1Translated;
                searchResult.QuantityQuality = quantityQuality;
                searchResult.ElutionTime = elutionTime;
                searchResult.ElutionTimeStart = elutionTimeStart;
                searchResult.ElutionTimeStop = elutionTimeStop;
                searchResult.IndexedRT = indexedRT;
                searchResult.IndexedIonMobility = indexedIonMobility;
                searchResult.PredictedRT = predictedRT;
                searchResult.PredictedIndexedRT = predictedIndexedRT;
                searchResult.MS1ProfileCorrelation = mS1ProfileCorrelation;
                searchResult.MS1Area = ms1Area;
                searchResult.Evidence = evidence;
                searchResult.SpectrumSimilarity = spectrumSimilarity;
                searchResult.Averagine = averagine;
                searchResult.MassEvidence = massEvidence;
                searchResult.CScore = cScore;
                searchResult.DecoyEvidence = decoyEvidence;
                searchResult.DecoyCScore = decoyCScore;

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
                            "Error parsing DIA-NN results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing DIA-NN Results in ParseDiaNNSynFileEntry: " + ex.Message);
                    }
                }
            }

            return false;
        }

        private bool ParseScanInfoFileHeaderLine(string lineIn, IDictionary<ScanInfoFileColumns, int> columnMapping)
        {
            var columnNames = new SortedDictionary<string, ScanInfoFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ScanNumber", ScanInfoFileColumns.Scan },
                { "ScanTime", ScanInfoFileColumns.ScanTime },
                { "ScanType", ScanInfoFileColumns.ScanType }
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (ScanInfoFileColumns resultColumn in Enum.GetValues(typeof(ScanInfoFileColumns)))
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
                SetErrorMessage("Error parsing the header line in a _ScanInfo.txt file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">DIA-NN results file (Dataset_report.tsv)</param>
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
                            // Uncomment to use a higher match error threshold in case some peptides reported by DIA-NN don't perfectly match the FASTA file
                            // const int MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD = 15;
                            // const int MATCH_ERROR_PERCENT_WARNING_THRESHOLD = 5;

                            success = CreateProteinModsFileWork(
                                baseName, inputFile,
                                synOutputFilePath, outputDirectoryPath,
                                PeptideHitResultTypes.DiaNN
                            // MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD,
                            // MATCH_ERROR_PERCENT_WARNING_THRESHOLD
                            );
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
        /// <param name="inputFile">report.tsv or report.parquet file</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="filteredSearchResults">Output: DIA-NN results</param>
        /// <param name="datasetNameToBaseNameMap">Output: Keys are full dataset names, values are abbreviated dataset names</param>
        /// <param name="usingParquetFile">Output: true if the input file is a .parquet file</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ReadDiaNNResults(
            FileInfo inputFile,
            ICollection<string> errorMessages,
            out List<DiaNNSearchResult> filteredSearchResults,
            out Dictionary<string, string> datasetNameToBaseNameMap,
            out bool usingParquetFile)
        {
            filteredSearchResults = new List<DiaNNSearchResult>();
            datasetNameToBaseNameMap = new Dictionary<string, string>();

            usingParquetFile = inputFile.Extension.Equals(".parquet", StringComparison.OrdinalIgnoreCase);

            try
            {
                if (!inputFile.Exists)
                {
                    OnErrorEvent("File not found: {0}", inputFile.FullName);
                    Console.WriteLine();
                    return false;
                }

                OnStatusEvent("Reading DIA-NN results file, " + PathUtils.CompactPathString(inputFile.FullName, 80));

                bool success;
                List<DiaNNSearchResult> unfilteredSearchResults;

                // ReSharper disable once ConvertIfStatementToConditionalTernaryExpression
                if (usingParquetFile)
                {
                    success = ReadReportParquetFile(inputFile, errorMessages, datasetNameToBaseNameMap, out unfilteredSearchResults);
                }
                else
                {
                    success = ReadReportTsvFile(inputFile, errorMessages, datasetNameToBaseNameMap, out unfilteredSearchResults);
                }

                if (!success)
                    return false;

                if (unfilteredSearchResults.Count == 0)
                {
                    OnWarningEvent("No results were found in file {0}", PathUtils.CompactPathString(inputFile.FullName, 110));

                    Console.WriteLine();
                    return true;
                }

                // Sort the SearchResults by dataset name, scan, charge, and ascending QValue
                unfilteredSearchResults.Sort(new DiaNNSearchResultsComparerDatasetScanChargeQValue());

                // Now filter the data
                var startIndex = 0;

                while (startIndex < unfilteredSearchResults.Count)
                {
                    // Find all the matches for the current result's scan
                    // (we sorted by dataset, then scan, so adjacent results will be from the same dataset, except when a new dataset is encountered)
                    // DIA-NN will typically report just one match

                    var endIndex = startIndex;

                    while (endIndex + 1 < unfilteredSearchResults.Count &&
                           unfilteredSearchResults[endIndex + 1].ScanNum == unfilteredSearchResults[startIndex].ScanNum)
                    {
                        endIndex++;
                    }

                    // Store the results for this scan
                    StoreSynMatches(unfilteredSearchResults, startIndex, endIndex, filteredSearchResults);

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

        private bool ReadReportParquetFile
        (
            FileInfo inputFile,
            ICollection<string> errorMessages,
            IDictionary<string, string> datasetNameToBaseNameMap,
            out List<DiaNNSearchResult> unfilteredSearchResults)
        {
            // Initialize the list that will hold all the records in the DIA-NN result file
            unfilteredSearchResults = new List<DiaNNSearchResult>();
            var baseNameToFullDatasetNameMap = new Dictionary<string, string>();

            try
            {
                // Cache the scan numbers and elution times in the _ScanInfo.txt files
                var success = ReadScanInfoFiles(inputFile.Directory, out var scanInfoByDataset);

                if (!success)
                {
                    return false;
                }

                var defaultDatasetName = GetDatasetNameFromFileName(inputFile.Name, PARQUET_REPORT_FILE_SUFFIX);

                var columnMapping = new Dictionary<DiaNNReportFileColumns, int>();

                var scanNumberOutOfRangeCount = 0;

                // Define the path to the _report.tsv file that will be created from the _report.parquet file
                var tsvFile = new FileInfo(Path.ChangeExtension(inputFile.FullName, ".tsv"));

                using var writer = new StreamWriter(new FileStream(tsvFile.FullName, FileMode.Create, FileAccess.Write, FileShare.Read));

                // Open the .parquet file and parse it
                using var reader = new ParquetFileReader(inputFile.FullName);

                for (var rowGroup = 0; rowGroup < reader.FileMetaData.NumRowGroups; ++rowGroup)
                {
                    using var rowGroupReader = reader.RowGroup(rowGroup);

                    var columnCount = rowGroupReader.MetaData.NumColumns;

                    var columnNames = new List<string>();

                    // Construct a tab-delimited list of the column names
                    for (var columnIndex = 0; columnIndex < columnCount; columnIndex++)
                    {
                        columnNames.Add(rowGroupReader.Column(columnIndex).ColumnDescriptor.Name);
                    }

                    var columnNameList = string.Join("\t", columnNames);

                    // Parse the header line
                    var headerParsed = ParseDiaNNResultsFileHeaderLine(columnNameList, columnMapping, out var columnNamesByEnum);

                    if (!headerParsed)
                    {
                        if (string.IsNullOrEmpty(mErrorMessage))
                        {
                            SetErrorMessage("Error parsing the column names in .parquet file " + inputFile.Name);
                        }

                        return false;
                    }

                    var rowCount = checked((int)rowGroupReader.MetaData.NumRows);

                    if (rowGroup > 0 && rowCount <= 0)
                        continue;

                    // Write the header names to the _report.tsv file
                    var headerNames = GetParquetTsvHeaderNames(columnNamesByEnum);
                    writer.WriteLine(string.Join("\t", headerNames));

                    var parquetData = new DiaNNParquetDataReader(columnMapping);

                    parquetData.ReadParquetData(rowGroupReader, rowCount);

                    var currentDataset = string.Empty;

                    for (var i = 0; i < rowCount; i++)
                    {
                        var searchResult = new DiaNNSearchResult
                        {
                            DatasetName = string.IsNullOrWhiteSpace(parquetData.DatasetName[i]) ? defaultDatasetName : parquetData.DatasetName[i]
                        };

                        if (!string.IsNullOrWhiteSpace(searchResult.DatasetName))
                        {
                            searchResult.FullDatasetName = Path.GetFileNameWithoutExtension(searchResult.DatasetName);

                            if (!datasetNameToBaseNameMap.ContainsKey(searchResult.FullDatasetName))
                            {
                                datasetNameToBaseNameMap.Add(searchResult.FullDatasetName, searchResult.DatasetName);
                            }

                            if (!baseNameToFullDatasetNameMap.ContainsKey(searchResult.DatasetName))
                            {
                                baseNameToFullDatasetNameMap.Add(searchResult.DatasetName, searchResult.FullDatasetName);
                            }

                            searchResult.DatasetFile = Path.HasExtension(searchResult.DatasetName)
                                ? searchResult.DatasetName
                                : string.Format("{0}.mzML", searchResult.DatasetName);
                        }

                        string datasetNameToUse;

                        if (string.IsNullOrWhiteSpace(currentDataset) || !currentDataset.Equals(searchResult.DatasetName))
                        {
                            currentDataset = searchResult.DatasetName;
                            datasetNameToUse = searchResult.DatasetName;
                        }
                        else
                        {
                            datasetNameToUse = string.Empty;
                        }

                        var outputData = parquetData.GetRowData(i, datasetNameToUse);

                        writer.WriteLine(string.Join("\t", outputData));

                        searchResult.ProteinGroup = parquetData.ProteinGroup[i];
                        searchResult.ProteinIDs = parquetData.ProteinIDs[i];
                        searchResult.ProteinNames = parquetData.ProteinNames[i];
                        searchResult.GeneNames = parquetData.GeneNames[i];

                        // searchResult.ProteinGroupQuantity = parquetData.ProteinGroupQuantity[i];
                        // searchResult.ProteinGroupNormalized = parquetData.ProteinGroupNormalized[i];

                        searchResult.ProteinGroupMaxLFQ = parquetData.ProteinGroupMaxLFQ[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.GenesQuantity = parquetData.GenesQuantity[i];
                        // searchResult.GenesNormalized = parquetData.GenesNormalized[i];
                        searchResult.GenesMaxLFQ = parquetData.GenesMaxLFQ[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.GenesMaxLFQUnique = parquetData.GenesMaxLFQUnique[i].ToString(CultureInfo.InvariantCulture);

                        searchResult.ModifiedSequence = parquetData.ModifiedSequence[i];

                        searchResult.Sequence = parquetData.StrippedSequence[i];

                        // Peptide sequences in DIA-NN report.tsv files do not include prefix or suffix residues (but this class will add them later)
                        // searchResult.PrevAA = parquetData.PrevAA[i];
                        // searchResult.NextAA = parquetData.NextAA[i];

                        searchResult.PrecursorId = parquetData.PrecursorId[i];

                        searchResult.Charge = parquetData.PrecursorCharge[i].ToString();

                        // Actual precursor m/z (only defined in .parquet files, not in .tsv files)
                        searchResult.PrecursorMzDiaNN = parquetData.PrecursorMz[i].ToString(CultureInfo.InvariantCulture);

                        searchResult.QValueDiaNN = parquetData.QValueDiaNN[i].ToString(CultureInfo.InvariantCulture);

                        searchResult.PEP = parquetData.Pep[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.GlobalQValue = parquetData.GlobalQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.ProteinQValue = parquetData.ProteinQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.ProteinGroupQValue = parquetData.ProteinGroupQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.GlobalProteinGroupQValue = parquetData.GlobalProteinGroupQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.GeneGroupQValue = parquetData.GeneGroupQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.TranslatedQValue = parquetData.TranslatedQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.Proteotypic = parquetData.Proteotypic[i].ToString();
                        searchResult.PrecursorQuantity = parquetData.PrecursorQuantity[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.PrecursorNormalized = parquetData.PrecursorNormalized[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.PrecursorTranslated = parquetData.PrecursorTranslated[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.TranslatedQuality = parquetData.TranslatedQuality[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.MS1Translated = parquetData.MS1Translated[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.QuantityQuality = parquetData.QuantityQuality[i].ToString(CultureInfo.InvariantCulture);

                        // Retention time is in minutes
                        searchResult.ElutionTime = parquetData.ElutionTime[i].ToString(CultureInfo.InvariantCulture);

                        searchResult.RTStart = parquetData.RtStart[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.RTStop = parquetData.RtStop[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.IndexedRT = parquetData.IndexedRT[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.PredictedRT = parquetData.PredictedRT[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.PredictedIndexedRT = parquetData.PredictedIndexedRT[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.FirstProteinDescription = parquetData.FirstProteinDescription[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.LibQValue = parquetData.LibQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.LibProteinGroupQValue = parquetData.LibProteinGroupQValue[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.MS1ProfileCorrelation = parquetData.Ms1ProfileCorrelation[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.MS1Area = parquetData.Ms1Area[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.Evidence = parquetData.Evidence[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.SpectrumSimilarity = parquetData.SpectrumSimilarity[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.Averagine = parquetData.Averagine[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.MassEvidence = parquetData.MassEvidence[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.CScore = parquetData.CScore[i].ToString(CultureInfo.InvariantCulture);

                        // searchResult.DecoyEvidence = decoyEvidence[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.DecoyCScore = decoyCScore[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.FragmentQuantRaw = fragmentQuantRaw[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.FragmentQuantCorrected = fragmentQuantCorrected[i].ToString(CultureInfo.InvariantCulture);
                        // searchResult.FragmentCorrelations = fragmentCorrelations[i].ToString(CultureInfo.InvariantCulture);

                        var scanNumber = LookupScanNumber(scanInfoByDataset, searchResult.DatasetName, searchResult.ElutionTime);
                        var startScan = LookupScanNumber(scanInfoByDataset, searchResult.DatasetName, searchResult.RTStart);
                        var endScan = LookupScanNumber(scanInfoByDataset, searchResult.DatasetName, searchResult.RTStop);

                        if (scanNumber < startScan || scanNumber > endScan)
                        {
                            scanNumberOutOfRangeCount++;

                            if (scanNumberOutOfRangeCount <= 5)
                            {
                                OnDebugEvent("Elution time (RT) is not between RT Start and RT End: {0} vs. time range {1} to {2}",
                                    searchResult.ElutionTime, searchResult.RTStart, searchResult.RTStop);
                            }
                        }

                        searchResult.ScanNum = scanNumber;
                        searchResult.Scan = searchResult.ScanNum.ToString();
                        searchResult.IonMobility = parquetData.IonMobility[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.IndexedIonMobility = parquetData.IndexedIonMobility[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.PredictedIonMobility = parquetData.PredictedIonMobility[i].ToString(CultureInfo.InvariantCulture);
                        searchResult.PredictedIndexedIonMobility = parquetData.PredictedIndexedIonMobility[i].ToString(CultureInfo.InvariantCulture);

                        var validSearchResult = ParseDiaNNSearchResult(
                            searchResult,
                            errorMessages,
                            i,
                            false);

                        if (validSearchResult)
                        {
                            unfilteredSearchResults.Add(searchResult);
                        }

                        // Update the progress
                        var percentComplete = i / (float)rowCount * 100;
                        UpdateProgress(percentComplete);
                    }
                }

                if (scanNumberOutOfRangeCount > 5)
                {
                    OnDebugEvent("Note: {0} PSMs in the .parquet file had an elution time not between RT Start and RT End", scanNumberOutOfRangeCount);
                }

                reader.Close();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in DiaNNResultsProcessor.ReadReportParquetFile", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
                return false;
            }
        }

        /// <summary>
        /// Load DIA-NN search results from a report.tsv file
        /// </summary>
        /// <remarks>If the file is found, but has no results, this method still returns true</remarks>
        /// <param name="inputFile">report.tsv file</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="datasetNameToBaseNameMap">Keys are full dataset names, values are abbreviated dataset names</param>
        /// <param name="unfilteredSearchResults">Output: DIA-NN results</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ReadReportTsvFile(
            FileSystemInfo inputFile,
            ICollection<string> errorMessages,
            IDictionary<string, string> datasetNameToBaseNameMap,
            out List<DiaNNSearchResult> unfilteredSearchResults)
        {
            // Initialize the list that will hold all the records in the DIA-NN result file
            unfilteredSearchResults = new List<DiaNNSearchResult>();
            var baseNameToFullDatasetNameMap = new Dictionary<string, string>();

            try
            {
                var columnMapping = new Dictionary<DiaNNReportFileColumns, int>();

                // Open the input file and parse it
                using var reader = new StreamReader(new FileStream(inputFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var headerParsed = false;
                var lineNumber = 0;

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
                        var success = ParseDiaNNResultsFileHeaderLine(lineIn, columnMapping, out _);

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
                        baseNameToFullDatasetNameMap,
                        ref currentDatasetFile);

                    if (validSearchResult)
                    {
                        unfilteredSearchResults.Add(searchResult);
                    }

                    // Update the progress
                    UpdateSynopsisFileCreationProgress(reader);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in DiaNNResultsProcessor.ReadReportTsvFile", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
                return false;
            }
        }

        private bool ReadScanInfoFiles(DirectoryInfo inputFileDirectory, out Dictionary<string, BinarySearchFindNearest> scanInfoByDataset)
        {
            scanInfoByDataset = new Dictionary<string, BinarySearchFindNearest>(StringComparer.OrdinalIgnoreCase);

            try
            {
                foreach (var scanInfoFile in inputFileDirectory.GetFiles("*_ScanInfo.txt"))
                {
                    OnStatusEvent("Reading scan info file,      " + PathUtils.CompactPathString(scanInfoFile.FullName, 80));

                    var datasetName = GetDatasetNameFromFileName(scanInfoFile.Name, "_ScanInfo.txt");
                    var scanNumbersByElutionTime = new Dictionary<double, double>();
                    var duplicateElutionTimes = new SortedSet<double>();

                    var columnMapping = new Dictionary<ScanInfoFileColumns, int>();

                    using var reader = new StreamReader(new FileStream(scanInfoFile.FullName, FileMode.Open, FileAccess.Read, FileShare.Read));

                    var headerParsed = false;

                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();

                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        if (!headerParsed)
                        {
                            // Parse the header line
                            var success = ParseScanInfoFileHeaderLine(lineIn, columnMapping);

                            if (!success)
                            {
                                if (string.IsNullOrEmpty(mErrorMessage))
                                {
                                    SetErrorMessage("Invalid header line in " + scanInfoFile.Name);
                                }

                                return false;
                            }

                            if (columnMapping[ScanInfoFileColumns.Scan] < 0)
                            {
                                SetErrorMessage(string.Format("Column 'ScanNumber' is missing from file {0}", scanInfoFile.FullName));
                                return false;
                            }

                            if (columnMapping[ScanInfoFileColumns.ScanTime] < 0)
                            {
                                SetErrorMessage(string.Format("Column 'ScanTime' is missing from file {0}", scanInfoFile.FullName));
                                return false;
                            }

                            headerParsed = true;
                            continue;
                        }

                        var splitLine = lineIn.TrimEnd().Split('\t');

                        DataUtilities.GetColumnValue(splitLine, columnMapping[ScanInfoFileColumns.Scan], out int scan);
                        DataUtilities.GetColumnValue(splitLine, columnMapping[ScanInfoFileColumns.ScanTime], out double elutionTime);
                        // DataUtilities.GetColumnValue(splitLine, columnMapping[ScanInfoFileColumns.ScanType], out int scanType);

                        if (scanNumbersByElutionTime.ContainsKey(elutionTime))
                        {
                            if (duplicateElutionTimes.Add(elutionTime) && duplicateElutionTimes.Count <= 5)
                            {
                                OnDebugEvent("Elution time {0} is listed more than once in the scan info file; ignoring subsequent occurrences", elutionTime);
                            }

                            continue;
                        }

                        scanNumbersByElutionTime.Add(elutionTime, scan);
                    }

                    if (duplicateElutionTimes.Count > 10)
                    {
                        OnDebugEvent("Note: {0} elution times were listed more than once in the _ScanInfo.txt file", duplicateElutionTimes.Count);
                    }

                    var elutionTimeToScanMap = new BinarySearchFindNearest();
                    elutionTimeToScanMap.AddData(scanNumbersByElutionTime);

                    scanInfoByDataset.Add(datasetName, elutionTimeToScanMap);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in DiaNNResultsProcessor.ReadScanInfoFiles", ex);
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
            // When processing peptides identified by isolating a specific precursor m/z, we assign the ranks using AssignRankByScore:
            // AssignRankByScore(searchResults, startIndex, endIndex);

            // Since this is DIA data, set RankScore to 1 for each result
            for (var i = startIndex; i <= endIndex; i++)
            {
                searchResults[i].RankScore = 1;
            }

            // The calling procedure already sorted by dataset name, scan, charge, and score; no need to re-sort

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
        /// <param name="writer">Writer</param>
        /// <param name="usingParquetFile">When true, the input file is a .parquet file and thus several columns should not be written to the output file</param>
        /// <param name="errorMessages">Error messages</param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            bool usingParquetFile,
            ICollection<string> errorMessages)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = DiaNNSynFileReader.GetColumnHeaderNamesAndIDs(usingParquetFile);

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
        /// <param name="writer">Writer</param>
        /// <param name="filteredSearchResults">Filtered search results</param>
        /// <param name="usingParquetFile">When true, the input file is a .parquet file and thus several columns should not be written to the output file</param>
        /// <param name="errorMessages">Error messages</param>
        private void WriteFilteredSearchResults(
            Dictionary<string, string> datasetNameToBaseNameMap,
            TextWriter writer,
            List<DiaNNSearchResult> filteredSearchResults,
            bool usingParquetFile,
            ICollection<string> errorMessages)
        {
            // Lookup the Dataset ID for each dataset (only if on the pnl.gov domain)
            var datasetIDs = LookupDatasetIDs(datasetNameToBaseNameMap.Keys.ToList());

            // Sort filteredSearchResults by ascending QValue, ascending Scan, Charge, and Peptide
            filteredSearchResults.Sort(new DiaNNSearchResultsComparerQValueScanChargePeptide());

            var index = 1;

            foreach (var result in filteredSearchResults)
            {
                GetBaseNameAndDatasetID(datasetNameToBaseNameMap, datasetIDs, result.FullDatasetName, out var baseDatasetName, out var datasetID);

                if (string.IsNullOrEmpty(baseDatasetName))
                {
                    OnWarningEvent("Method GetBaseNameAndDatasetID returned an empty string for baseDatasetName for result {0}; this is unexpected", result);
                    WriteSearchResultToFile(index, result.DatasetName, datasetID, writer, result, usingParquetFile, errorMessages);
                }
                else
                {
                    WriteSearchResultToFile(index, baseDatasetName, datasetID, writer, result, usingParquetFile, errorMessages);
                }

                index++;
            }
        }

        /// <summary>
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID">Result ID</param>
        /// <param name="baseDatasetName">Base dataset name</param>
        /// <param name="datasetID">Dataset ID</param>
        /// <param name="writer">Writer</param>
        /// <param name="searchResult">Search result</param>
        /// <param name="usingParquetFile">When true, the input file is a .parquet file and thus several columns should not be written to the output file</param>
        /// <param name="errorMessages">Error messages</param>
        private void WriteSearchResultToFile(
            int resultID,
            string baseDatasetName,
            int datasetID,
            TextWriter writer,
            DiaNNSearchResult searchResult,
            bool usingParquetFile,
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
                    searchResult.IonMobility,
                    searchResult.Charge,
                    searchResult.PrecursorMZ,
                    searchResult.MH,
                    searchResult.CalculatedMonoMass,
                    searchResult.Sequence,
                    searchResult.ModificationList,
                    searchResult.ProteinGroup,
                    searchResult.ProteinIDs,
                    searchResult.ProteinNames,
                    searchResult.GeneNames,
                    searchResult.NumberOfTrypticTermini.ToString(), // NTT
                };

                if (!usingParquetFile)
                {
                    data.Add(searchResult.ProteinGroupQuantity);
                    data.Add(searchResult.ProteinGroupNormalized);
                }

                data.Add(searchResult.ProteinGroupMaxLFQ);

                if (!usingParquetFile)
                {
                    data.Add(searchResult.GenesQuantity);
                    data.Add(searchResult.GenesNormalized);
                }

                data.Add(searchResult.GenesMaxLFQ);
                data.Add(searchResult.GenesMaxLFQUnique);
                data.Add(searchResult.QValueDiaNN);
                data.Add(searchResult.PEP);
                data.Add(searchResult.GlobalQValue);
                data.Add(searchResult.ProteinQValue);
                data.Add(searchResult.ProteinGroupQValue);
                data.Add(searchResult.GlobalProteinGroupQValue);
                data.Add(searchResult.GeneGroupQValue);
                data.Add(searchResult.TranslatedQValue);
                data.Add(searchResult.Proteotypic);
                data.Add(searchResult.PrecursorQuantity);
                data.Add(searchResult.PrecursorNormalized);

                if (!usingParquetFile)
                {
                    data.Add(searchResult.PrecursorTranslated);
                    data.Add(searchResult.TranslatedQuality);
                    data.Add(searchResult.MS1Translated);
                }

                data.Add(searchResult.QuantityQuality);
                data.Add(searchResult.ElutionTime);
                data.Add(searchResult.RTStart);
                data.Add(searchResult.RTStop);
                data.Add(searchResult.IndexedRT);
                data.Add(searchResult.IndexedIonMobility);
                data.Add(searchResult.PredictedRT);
                data.Add(searchResult.PredictedIndexedRT);
                data.Add(searchResult.MS1ProfileCorrelation);
                data.Add(searchResult.MS1Area);
                data.Add(searchResult.Evidence);

                if (!usingParquetFile)
                {
                    data.Add(searchResult.SpectrumSimilarity);
                    data.Add(searchResult.Averagine);
                }

                data.Add(searchResult.MassEvidence);

                if (!usingParquetFile)
                {
                    data.Add(searchResult.CScore);
                    data.Add(searchResult.DecoyEvidence);
                    data.Add(searchResult.DecoyCScore);
                }

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
            return string.Format("{0} results processor", TOOL_NAME);
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
                // Sort ascending
                if (x.QValue > y.QValue)
                {
                    return 1;
                }

                if (x.QValue < y.QValue)
                {
                    return -1;
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

                // Sort ascending
                if (x.QValue > y.QValue)
                {
                    return 1;
                }

                if (x.QValue < y.QValue)
                {
                    return -1;
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
