// This class reads in an MS-GF+ results file (.txt format) and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.
//
// The modification definition information is determined from the MS-GF+ parameter file
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 8/12/2011
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads a MS-GF+ results file (e.g. Dataset.tsv) and creates the first hits and synopsis files
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads MS-GF+ results file Dataset.tsv
    /// </para>
    /// <para>
    /// 2) It calls CreateFHTorSYNResultsFile to create the _fht.txt and _syn.txt files
    /// </para>
    /// <para>
    /// 3) ParseMSGFPlusResultsFileHeaderLine reads the header line to determine the column mapping
    ///      columnMapping = new Dictionary of MSGFPlusResultsFileColumns, int
    /// </para>
    /// <para>
    /// 4) ParseMSGFPlusResultsFileEntry reads each data line and stores in an instance of MSGFPlusSearchResult, which is a private struct
    ///    The data is stored in a list
    ///      searchResultsCurrentScan = new List of MSGFPlusSearchResult
    /// </para>
    /// <para>
    /// 5) Filter-passing results for a given scan are stored in another list
    ///      searchResultsPrefiltered = new List of MSGFPlusSearchResult
    /// </para>
    /// <para>
    /// 6) Once the entire .tsv has been read, searchResultsPrefiltered is sorted by score
    /// </para>
    /// <para>
    /// 7) Data to be written to the _syn.txt or _fht.txt file is stored in another list
    ///      filteredSearchResults = new List of MSGFPlusSearchResult
    ///    For syn files, all filter-passing PSMs are stored
    ///    For fht files, only the highest scoring PSM is stored for each scan/charge combo
    /// </para>
    /// <para>
    /// 8) SortAndWriteFilteredSearchResults performs one more sort, then writes out to disk
    /// </para>
    /// </remarks>
    public class MSGFPlusResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: da, fht, frag, kv, methylation, msgfdb, novo, pre, Prefiltered, struct, structs, tda, tsv, udt

        /// <summary>
        /// Constructor
        /// </summary>
        public MSGFPlusResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "January 13, 2022";

            mModMassRegEx = new Regex(MSGFPlus_MOD_MASS_REGEX, REGEX_OPTIONS);

            mPeptideCleavageStateCalculator = new PeptideCleavageStateCalculator();
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(PeptideCleavageStateCalculator.StandardCleavageAgent.Trypsin);

            mNumericModErrors = 0;
        }

        /// <summary>
        /// MS-GF+ tool name
        /// </summary>
        public const string TOOL_NAME = "MSGFPlus";

        /// <summary>
        /// MSGF-DB results file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_MSGFDB_FILE = "_msgfdb";

        /// <summary>
        /// MS-GF+ results file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_MSGFPLUS_FILE = "_msgfplus";

        /// <summary>
        /// N-terminus symbol used by MS-GF+
        /// </summary>
        public const string N_TERMINUS_SYMBOL_MSGFPlus = "_.";

        /// <summary>
        /// C-terminus symbol used by MS-GF+
        /// </summary>
        public const string C_TERMINUS_SYMBOL_MSGFPlus = "._";

        /// <summary>
        /// Default synopsis file Spec E-Value threshold
        /// </summary>
        /// <remarks>
        /// <para>
        /// Filter passing peptides have Spec_E-value less than 5E-7 or E-Value (EValue) less than 0.75 or Q-Value (QValue) less than 1%
        /// </para>
        /// <para>
        /// This filter is also used by MSPathFinder
        /// </para>
        /// </remarks>
        public const float DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD = 5E-07f;

        /// <summary>
        /// Default synopsis file E-value threshold
        /// </summary>
        /// <remarks>
        /// <para>
        /// Filter passing peptides have Spec_E-value less than 5E-7 or E-Value (EValue) less than 0.75 or Q-Value (QValue) less than 1%
        /// </para>
        /// <para>
        /// This filter is also used by MSFragger
        /// </para>
        /// </remarks>
        public const float DEFAULT_SYN_FILE_EVALUE_THRESHOLD = 0.75f;

        private const string SEARCH_ENGINE_NAME = "MS-GF+";

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        // Match mod masses (positive or negative) at start, e.g.
        // ReSharper disable CommentTypo
        // +57.021HWWTLTTDRINK         matches +57.021
        // -57.021+42.011HWWTLTTDRINK  matches -57.021+42.011 (two separate mods)
        // +42.011MDHTPQSQLK           matches +42.011
        // ReSharper restore CommentTypo
        private const string MSGFPlus_N_TERMINAL_MOD_MASS_REGEX = @"^([0-9\.\+\-]+)";

        private const string MSGFPlus_MOD_MASS_REGEX = @"([+-][0-9\.]+)";

        private const string PROTEIN_AND_TERM_SYMBOLS_REGEX = @"([^;]+)\(pre=(.),post=(.)\)";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to TSV file created by MzidToTsvConverter.exe from the MS-GF+ .mzid file
        /// </summary>
        private enum MSGFPlusResultsFileColumns
        {
            SpectrumFile = 0,
            SpecIndex = 1,               // SpecID in MS-GF+
            Scan = 2,
            ScanTimeMinutes = 3,         // Added to MzidToTsvConverter in April 2019
            FragMethod = 4,
            PrecursorMZ = 5,
            PMErrorDa = 6,               // Corresponds to PrecursorError(Da)
            PMErrorPPM = 7,              // Corresponds to PrecursorError(ppm)
            Charge = 8,
            Peptide = 9,
            Protein = 10,
            DeNovoScore = 11,
            MSGFScore = 12,
            SpecProb_EValue = 13,
            PValue_EValue = 14,
            FDR_QValue = 15,             // Only present if searched using -tda 1
            PepFDR_PepQValue = 16,       // Only present if searched using -tda 1
            EFDR = 17,                   // Only present if did not search using -tda 1
            IMSScan = 18,                // Only present for MSGFDB_IMS results
            IMSDriftTime = 19,           // Only present for MSGFDB_IMS results
            IsotopeError = 20            // Only reported by MS-GF+
        }

        private enum FilteredOutputFileTypeConstants
        {
            SynFile = 0,
            FHTFile = 1
        }

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file (Dataset.tsv) created MzidToTsvConverter
        /// </summary>
        /// <remarks>
        /// These columns hold data that this class will use when creating the synopsis file
        /// </remarks>
        private struct MSGFPlusSearchResult
        {
            // ReSharper disable once NotAccessedField.Local

            /// <summary>
            /// Spectrum file name
            /// </summary>
            public string SpectrumFileName;

            /// <summary>
            /// Spectrum index
            /// </summary>
            public string SpecIndex;

            /// <summary>
            /// Scan number
            /// </summary>
            public string Scan;

            /// <summary>
            /// Numeric value of Scan
            /// </summary>
            public int ScanNum;

            /// <summary>
            /// Fragmentation method
            /// </summary>
            public string FragMethod;

            /// <summary>
            /// Precursor ion m/z (observed value)
            /// </summary>
            public string PrecursorMZ;

            /// <summary>
            /// Precursor mass error, in Da
            /// </summary>
            /// <remarks>
            /// From column "PMError(Da)"; MS-GF+ stores this value as Observed - Theoretical
            /// </remarks>
            public string PMErrorDa;

            /// <summary>
            /// Precursor mass error, in ppm
            /// </summary>
            /// <remarks>
            /// From column "PMError(ppm)"; MS-GF+ stores this value as Observed - Theoretical
            /// </remarks>
            public string PMErrorPPM;

            /// <summary>
            /// Peptide M+H value
            /// </summary>
            public string MH;

            /// <summary>
            /// Charge state
            /// </summary>
            public string Charge;

            /// <summary>
            /// Numeric value of Charge
            /// </summary>
            public short ChargeNum;

            /// <summary>
            /// Peptide sequence, including prefix, suffix, and any mod symbols or mod masses
            /// </summary>
            public string Peptide;

            /// <summary>
            /// Protein name
            /// </summary>
            public string Protein;

            /// <summary>
            /// Number of tolerable termini
            /// </summary>
            public string NTT;

            /// <summary>
            /// De-novo score
            /// </summary>
            public string DeNovoScore;

            /// <summary>
            /// MSGF score
            /// </summary>
            public string MSGFScore;

            /// <summary>
            /// Spec E-value
            /// </summary>
            /// <remarks>
            /// <para>
            /// Smaller values are better scores (e.g. 1E-9 is better than 1E-6)
            /// </para>
            /// <para>
            /// Renamed from SpecProb to SpecEValue when MSGF-DB was renamed to MS-GF+
            /// </para>
            /// </remarks>
            public string SpecEValue;

            /// <summary>
            /// Spec E-value rank
            /// </summary>
            public double SpecEValueNum;

            /// <summary>
            /// E-value
            /// </summary>
            /// <remarks>
            /// <para>Smaller values are better scores (e.g. 1E-7 is better than 1E-3)</para>
            /// <para>
            /// Renamed from PValue to EValue when MSGF-DB was renamed to MS-GF+
            /// </para>
            /// </remarks>
            public string EValue;

            /// <summary>
            /// Numeric value of EValue
            /// </summary>
            public double EValueNum;

            /// <summary>
            /// Holds FDR when a target/decoy search was used;
            /// Holds EFDR when a non-decoy search was used;
            /// Holds QValue for MS-GF+
            /// </summary>
            public string QValue;

            /// <summary>
            /// Numeric value of QValue
            /// </summary>
            public double QValueNum;

            /// <summary>
            /// Pep Q-value
            /// </summary>
            /// <remarks>
            /// Only used when target/decoy search was used; holds PepQValue for MS-GF+
            /// </remarks>
            public string PepQValue;

            /// <summary>
            /// Rank Spec E-value
            /// </summary>
            public int RankSpecEValue;

            /// <summary>
            /// IMS Scan
            /// </summary>
            public int IMSScan;

            /// <summary>
            /// IMS drift time
            /// </summary>
            public string IMSDriftTime;

            /// <summary>
            /// Isotope error
            /// </summary>
            /// <remarks>
            /// Only used by MS-GF+
            /// </remarks>
            public int IsotopeError;

            /// <summary>
            /// Reset stored values to empty strings and zeros
            /// </summary>
            public void Clear()
            {
                SpectrumFileName = string.Empty;
                SpecIndex = string.Empty;
                ScanNum = 0;
                FragMethod = string.Empty;
                PrecursorMZ = string.Empty;
                PMErrorDa = string.Empty;
                PMErrorPPM = string.Empty;
                MH = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                Peptide = string.Empty;
                Protein = string.Empty;
                NTT = string.Empty;
                DeNovoScore = string.Empty;
                MSGFScore = string.Empty;
                SpecEValue = string.Empty;
                SpecEValueNum = 0;
                EValue = string.Empty;
                EValueNum = 0;
                QValue = string.Empty;
                QValueNum = 0;
                PepQValue = string.Empty;
                RankSpecEValue = 0;
                IMSScan = 0;
                IMSDriftTime = string.Empty;
                IsotopeError = 0;
            }

            /// <summary>
            /// Show scan, peptide, and Spec E-value
            /// </summary>
            public override string ToString()
            {
                return string.Format("Scan {0}: {1}, SpecEValue {2}", ScanNum, Peptide, SpecEValue);
            }
        }

        private struct ScanGroupInfo
        {
            public int ScanGroupID;
            public short Charge;
            public int Scan;
        }

        private struct TerminusChars
        {
            public char NTerm;
            public char CTerm;
        }

        private readonly PeptideCleavageStateCalculator mPeptideCleavageStateCalculator;

        /// <summary>
        /// This variable keeps track of the number of PSMs whose computed monoisotopic mass
        /// does not agree with the monoisotopic mass computed from the precursor m/z, within 1.5x of the match tolerance
        /// </summary>
        private int mPrecursorMassErrorWarningCount;

        /// <summary>
        /// Precursor match tolerance read from the MS-GF+ parameter file
        /// </summary>
        private PrecursorMassTolerance mPrecursorMassTolerance;

        /// <summary>
        /// Looks for numeric mods in MS-GF+ results
        /// For example, +14.016 in K.LQVPAGK+14.016ANPSPPIGPALGQR.G
        /// </summary>
        private readonly Regex mModMassRegEx;

        private int mNumericModErrors;

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod symbol, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        private void AddDynamicAndStaticResidueMods(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            var mostRecentLetter = '-';
            var residueLocInPeptide = 0;

            var sequence = searchResult.PeptideSequenceWithMods;

            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var character = sequence[index];

                if (StringUtilities.IsLetterAtoZ(character))
                {
                    mostRecentLetter = character;
                    residueLocInPeptide++;

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ResidueModificationType.StaticMod)
                        {
                            continue;
                        }

                        var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                        if (!modificationDefinition.TargetResiduesContain(character))
                        {
                            continue;
                        }

                        // Match found; add this modification
                        var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocInPeptide);

                        searchResult.SearchResultAddModification(
                            modificationDefinition, character, residueLocInPeptide,
                            residueTerminusState, updateModOccurrenceCounts);
                    }
                }
                else if (StringUtilities.IsLetterAtoZ(mostRecentLetter))
                {
                    var success = searchResult.SearchResultAddDynamicModification(character, mostRecentLetter, residueLocInPeptide, searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
                    if (success)
                    {
                        continue;
                    }

                    var errorMessage = searchResult.ErrorMessage;
                    if (string.IsNullOrEmpty(errorMessage))
                    {
                        errorMessage = "SearchResultAddDynamicModification returned false for symbol " + character;
                    }
                    SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
                }
                // ReSharper disable once RedundantIfElseBlock
                else
                {
                    // We found a modification symbol but mostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }
        }

        /// <summary>
        /// Adds or updates the prefix and suffix residues to the peptide, as defined in kvProteinInfo
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="kvProteinInfo"></param>
        /// <returns>Peptide sequence with N-terminal and C-Terminal residues</returns>
        private string AddUpdatePrefixAndSuffixResidues(string peptide, KeyValuePair<string, TerminusChars> kvProteinInfo)
        {
            if (peptide.IndexOf('.') < 0)
            {
                return kvProteinInfo.Value.NTerm + "." + peptide + "." + kvProteinInfo.Value.CTerm;
            }

            string peptideNew;

            if (peptide.Length >= 2)
            {
                if (peptide[1] == '.')
                {
                    // Peptide already has the N-terminal residue
                    // Replace it using kvProteinInfo
                    peptideNew = kvProteinInfo.Value.NTerm + "." + peptide.Substring(2);
                }
                else if (peptide[0] == '.')
                {
                    peptideNew = kvProteinInfo.Value.NTerm + peptide;
                }
                else
                {
                    peptideNew = kvProteinInfo.Value.NTerm + "." + peptide;
                }
            }
            else
            {
                peptideNew = peptide;
            }

            if (peptideNew.Length < 4)
            {
                return peptideNew;
            }

            if (peptideNew[peptideNew.Length - 2] == '.')
            {
                // Peptide already has the C-terminal residue
                // Replace it using kvProteinInfo
                peptideNew = peptideNew.Substring(0, peptideNew.Length - 2) + "." + kvProteinInfo.Value.CTerm;
            }
            else if (peptideNew[peptideNew.Length - 1] == '.')
            {
                peptideNew += kvProteinInfo.Value.CTerm;
            }
            else
            {
                peptideNew = peptideNew + "." + kvProteinInfo.Value.CTerm;
            }

            return peptideNew;
        }

        private void AppendToScanGroupDetails(
            ICollection<ScanGroupInfo> scanGroupDetails,
            IDictionary<string, bool> scanGroupCombo,
            ScanGroupInfo udtScanGroupInfo,
            ref int currentScanGroupID,
            ref int nextScanGroupID)
        {
            var chargeScanComboText = udtScanGroupInfo.Charge + "_" + udtScanGroupInfo.Scan;

            if (scanGroupCombo.ContainsKey(chargeScanComboText))
            {
                return;
            }

            if (currentScanGroupID < 0)
            {
                currentScanGroupID = nextScanGroupID;
                nextScanGroupID++;
            }

            udtScanGroupInfo.ScanGroupID = currentScanGroupID;

            scanGroupDetails.Add(udtScanGroupInfo);
            scanGroupCombo.Add(chargeScanComboText, true);
        }

        private void AppendToSearchResults(
            ICollection<MSGFPlusSearchResult> searchResults,
            MSGFPlusSearchResult udtSearchResult,
            Dictionary<string, TerminusChars> proteinInfo)
        {
            if (proteinInfo.Count == 0)
            {
                searchResults.Add(udtSearchResult);
            }
            else
            {
                foreach (var kvEntry in proteinInfo)
                {
                    udtSearchResult.Protein = kvEntry.Key;
                    udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, kvEntry);

                    searchResults.Add(udtSearchResult);
                }
            }
        }

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        private void AssignRankByScore(
            IList<MSGFPlusSearchResult> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assigned per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            if (startIndex == endIndex)
            {
                // Only one result
                var currentResult = searchResults[startIndex];
                currentResult.RankSpecEValue = 1;
                searchResults[startIndex] = currentResult;
                return;
            }

            // Duplicate a portion of searchResults so that we can sort by ascending Spec E-Value

            var resultsSubset = new Dictionary<int, MSGFPlusSearchResult>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsBySpecProb = (from item in resultsSubset orderby item.Value.SpecEValueNum select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsBySpecProb)
            {
                var currentResult = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = currentResult.SpecEValueNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(currentResult.SpecEValueNum - lastValue) > double.Epsilon)
                    {
                        lastValue = currentResult.SpecEValueNum;
                        currentRank++;
                    }
                }

                currentResult.RankSpecEValue = currentRank;

                // Because this is a list of structs, we have to copy currentResult back into the current position in searchResults
                searchResults[entry.Key] = currentResult;
            }
        }

        private bool AddModificationsAndComputeMass(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            try
            {
                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts);

                // Make sure .PeptideSequenceWithMods does not have any generic mod masses
                // It should only have mod symbols
                var match = mModMassRegEx.Match(searchResult.PeptideSequenceWithMods);
                if (match.Success)
                {
                    // Modification mass did not have a symbol associated with it in the _ModDefs.txt file
                    // We could try to handle this, listing the modification mass in place of the modification symbol in the _ModDetails.txt file,
                    // but will instead abort processing

                    mNumericModErrors++;

                    if (mNumericModErrors < 250)
                    {
                        var localErrorMessage = string.Format(
                            "Search result contains a numeric mod mass that could not be associated with a modification symbol; ResultID = {0}, ModMass = {1}",
                            searchResult.ResultID, match.Value);

                        SetErrorMessage(localErrorMessage);
                    }
                    else if (mNumericModErrors == 250)
                    {
                        SetErrorMessage("Too many numeric mod mass results have been found; suppressing further logging");
                    }

                    return false;
                }

                // Parse .PeptideSequenceWithMods to determine the modified residues present
                AddDynamicAndStaticResidueMods(searchResult, updateModOccurrenceCounts);

                // Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
                // Since InSpecT allows a terminal peptide residue to be modified twice, we'll allow that to happen,
                //  even though, biologically, that's typically not possible
                // However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus
                //  (where two COOH groups are present)
                searchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, updateModOccurrenceCounts);

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

        private short ComputeCleavageState(string sequenceWithMods)
        {
            // Remove any non-letter characters before calling .ComputeCleavageState()
            var cleanSequence = GetCleanSequence(sequenceWithMods, out var prefix, out var suffix);

            var cleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(cleanSequence, prefix, suffix);

            return Convert.ToInt16(cleavageState);
        }

        /// <summary>
        /// Compute the delta mass, in ppm, optionally correcting for C13 isotopic selection errors
        /// </summary>
        /// <remarks>This method should only be called when column PMError(Da) is present (and PMError(ppm) is not present)</remarks>
        /// <param name="precursorErrorDa">Mass error (Observed - theoretical)</param>
        /// <param name="precursorMZ">Precursor m/z</param>
        /// <param name="charge">Precursor charge</param>
        /// <param name="peptideMonoisotopicMass">Monoisotopic mass of the peptide</param>
        /// <param name="adjustPrecursorMassForC13"></param>
        /// <returns>DelM, in ppm</returns>
        private double ComputeDelMCorrectedPPM(
            double precursorErrorDa,
            double precursorMZ,
            int charge,
            double peptideMonoisotopicMass,
            bool adjustPrecursorMassForC13)
        {
            // Compute the observed precursor monoisotopic mass
            var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, charge, 0);

            return SearchResultsBaseClass.ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, peptideMonoisotopicMass, adjustPrecursorMassForC13);
        }

        /// <summary>
        /// Compute the monoisotopic mass of the peptide
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="totalModMass"></param>
        private double ComputePeptideMass(string peptide, double totalModMass)
        {
            var cleanSequence = GetCleanSequence(peptide);

            return ComputePeptideMassForCleanSequence(cleanSequence, totalModMass);
        }

        /// <summary>
        /// Construct the peptide to protein map file path
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="mts">If true, the map file will end with MTS.txt; otherwise, just .txt</param>
        /// <returns>_PepToProtMap file that corresponds to the input file</returns>
        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var suffixesToFind = new List<string> {
                "_msgfplus_syn",
                "_msgfplus_fht",
                "_msgfdb_syn",
                "_msgfdb_fht",
                "_syn",
                "_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// Parses the digits in modDigits to convert them to one or more modification symbols
        /// </summary>
        /// <param name="currentResidue"></param>
        /// <param name="modDigits">Example: +57.021 or +79.9663+14.0157 or -18.0106</param>
        /// <param name="modSymbols"></param>
        /// <param name="dynModSymbols"></param>
        /// <param name="msgfPlusModInfo"></param>
        /// <param name="nTerminalMod"></param>
        /// <param name="possibleCTerminalMod"></param>
        /// <param name="modMassFound"></param>
        /// <param name="containsStaticMod"></param>
        /// <returns>True if success; false if a problem</returns>
        private bool ConvertMSGFModMassesToSymbols(
            string currentResidue,
            string modDigits,
            out string modSymbols,
            out string dynModSymbols,
            IReadOnlyList<MSGFPlusParamFileModExtractor.ModInfo> msgfPlusModInfo,
            bool nTerminalMod,
            bool possibleCTerminalMod,
            out double modMassFound,
            out bool containsStaticMod)
        {
            double bestMassDiff = 0;
            var modSymbolsFound = 0;
            var symbolBestMatch = ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
            var residuesBestBatch = string.Empty;

            modSymbols = string.Empty;
            dynModSymbols = string.Empty;
            modMassFound = 0;
            containsStaticMod = false;

            foreach (Match match in mModMassRegEx.Matches(modDigits))
            {
                var modMassText = match.Value;

                // Convert modMass to a mass value
                var modMass = double.Parse(modMassText);
                modMassFound += modMass;

                var matchFound = false;
                int bestMatchIndex;

                while (true)
                {
                    bestMatchIndex = -1;

                    // Step through the known modifications to find the closest match
                    for (var index = 0; index <= msgfPlusModInfo.Count - 1; index++)
                    {
                        var testMod = true;

                        if (nTerminalMod)
                        {
                            // Only test N-terminal mods
                            if (!(msgfPlusModInfo[index].ModType is MSGFPlusParamFileModExtractor.MSGFPlusModType.DynNTermPeptide or MSGFPlusParamFileModExtractor.MSGFPlusModType.DynNTermProtein))
                            {
                                testMod = false;
                            }
                        }
                        else if (!possibleCTerminalMod)
                        {
                            // Skip C-terminal mods since we're not at the C-terminus
                            if (msgfPlusModInfo[index].ModType is MSGFPlusParamFileModExtractor.MSGFPlusModType.DynCTermPeptide or MSGFPlusParamFileModExtractor.MSGFPlusModType.DynCTermProtein)
                            {
                                testMod = false;
                            }
                        }

                        if (!testMod)
                            continue;

                        var candidateMassDiff = Math.Abs(msgfPlusModInfo[index].ModMassVal - modMass);
                        if (!(candidateMassDiff < 0.25))
                            continue;

                        // Possible match found
                        var updateCandidate = false;

                        if (bestMatchIndex < 0 ||
                            candidateMassDiff < bestMassDiff)
                        {
                            updateCandidate = true;
                        }
                        else if (Math.Abs(candidateMassDiff - bestMassDiff) < float.Epsilon &&
                                 symbolBestMatch == '-' &&
                                 msgfPlusModInfo[index].ModSymbol != '-')
                        {
                            // Masses are the same, but the existing candidate is a static mod

                            // If this new one is a dynamic mod, switch to it, but only if residuesBestBatch does not contain the residue while
                            // the residues for the new candidate does contain the residue

                            if (!residuesBestBatch.Contains(currentResidue) &&
                                msgfPlusModInfo[index].Residues.Contains(currentResidue))
                            {
                                updateCandidate = true;
                            }
                        }

                        if (!updateCandidate)
                        {
                            continue;
                        }

                        bestMatchIndex = index;
                        bestMassDiff = candidateMassDiff;
                        symbolBestMatch = msgfPlusModInfo[index].ModSymbol;
                        residuesBestBatch = msgfPlusModInfo[index].Residues;
                        matchFound = true;
                    }

                    if (!matchFound && nTerminalMod)
                    {
                        // Set this to false, then search again
                        nTerminalMod = false;
                    }
                    else if (!matchFound && !possibleCTerminalMod)
                    {
                        // Set this to true, then search again
                        possibleCTerminalMod = true;
                    }
                    else
                    {
                        // matchFound is true
                        break;
                    }
                }

                if (matchFound)
                {
                    // Match found; use the mod symbol
                    modSymbols += msgfPlusModInfo[bestMatchIndex].ModSymbol;
                    modSymbolsFound++;

                    if (msgfPlusModInfo[bestMatchIndex].ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod)
                    {
                        containsStaticMod = true;
                    }
                    else
                    {
                        dynModSymbols += msgfPlusModInfo[bestMatchIndex].ModSymbol;
                    }
                }
                else
                {
                    // Match not found; use the mass value
                    modSymbols += modMass;
                }
            }

            return modSymbolsFound > 0;
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MS-GF+
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="scanGroupFilePath"></param>
        /// <param name="msgfPlusModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
        /// <param name="isMsgfPlus">Output parameter: this method will set this to True if we're processing MS-GF+ results</param>
        /// <param name="specIdToIndex"></param>
        /// <param name="filteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateFHTorSYNResultsFile(
            string inputFilePath,
            string outputFilePath,
            string scanGroupFilePath,
            IReadOnlyList<MSGFPlusParamFileModExtractor.ModInfo> msgfPlusModInfo,
            out bool isMsgfPlus,
            IDictionary<string, int> specIdToIndex,
            FilteredOutputFileTypeConstants filteredOutputFileType)
        {
            var searchResultsCurrentScan = new List<MSGFPlusSearchResult>();
            var searchResultsPrefiltered = new List<MSGFPlusSearchResult>();

            isMsgfPlus = false;

            try
            {
                var scanGroupDetails = new List<ScanGroupInfo>();
                var scanGroupCombo = new Dictionary<string, bool>();

                mPrecursorMassErrorWarningCount = 0;

                // Look for custom amino acids
                var customAA =
                    (from item in msgfPlusModInfo
                     where item.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.CustomAA
                     select item).ToList();

                foreach (var customAADef in customAA)
                {
                    var aminoAcidSymbol = customAADef.Residues[0];
                    var empiricalFormulaString = customAADef.ModMass;
                    var aminoAcidMass = customAADef.ModMassVal;

                    try
                    {
                        var elementalComposition = PeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormulaString);

                        mPeptideSeqMassCalculator.SetAminoAcidMass(aminoAcidSymbol, aminoAcidMass);
                        mPeptideSeqMassCalculator.SetAminoAcidAtomCounts(aminoAcidSymbol, elementalComposition);
                    }
                    catch (Exception ex)
                    {
                        ReportError(ex.Message);
                    }
                }

                try
                {
                    // Open the input file and parse it
                    // Initialize the stream reader and the stream writer

                    using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                    using var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                    var errorMessages = new List<string>();
                    var headerParsed = false;
                    var includeFDRandPepFDR = false;
                    var includeEFDR = false;
                    var includeIMSFields = false;

                    var nextScanGroupID = 1;
                    scanGroupDetails.Clear();
                    scanGroupCombo.Clear();

                    // Initialize the array that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<MSGFPlusSearchResult>();

                    // Initialize a dictionary that tracks the peptide sequence for each combo of scan and charge
                    // Keys are Scan_Charge, values track the clean sequence, the associated protein name, and the protein number for that name
                    // Note that we can only track protein numbers if the FASTA file path was provided at the command line
                    var scanChargeFirstHit = new Dictionary<string, FirstHitInfo>();

                    var columnMapping = new Dictionary<MSGFPlusResultsFileColumns, int>();

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
                            var validHeader = ParseMSGFPlusResultsFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }

                            headerParsed = true;

                            if (columnMapping[MSGFPlusResultsFileColumns.FDR_QValue] >= 0 ||
                                columnMapping[MSGFPlusResultsFileColumns.PepFDR_PepQValue] >= 0)
                            {
                                includeFDRandPepFDR = true;
                            }
                            else if (columnMapping[MSGFPlusResultsFileColumns.EFDR] >= 0)
                            {
                                includeEFDR = true;
                            }

                            if (columnMapping[MSGFPlusResultsFileColumns.IMSDriftTime] >= 0)
                            {
                                includeIMSFields = true;
                            }

                            if (columnMapping[MSGFPlusResultsFileColumns.IsotopeError] >= 0)
                            {
                                isMsgfPlus = true;
                            }

                            // Write the header line
                            WriteSynFHTFileHeader(writer, errorMessages, includeFDRandPepFDR, includeEFDR, includeIMSFields, isMsgfPlus);

                            continue;
                        }

                        var validSearchResult = ParseMSGFPlusResultsFileEntry(
                            lineIn, isMsgfPlus, msgfPlusModInfo,
                            searchResultsCurrentScan, errorMessages,
                            columnMapping, ref nextScanGroupID, scanGroupDetails,
                            scanGroupCombo, specIdToIndex);

                        if (!validSearchResult || searchResultsCurrentScan.Count == 0)
                        {
                            continue;
                        }

                        if (filteredOutputFileType == FilteredOutputFileTypeConstants.SynFile)
                        {
                            // Synopsis file
                            validSearchResult = MSGFPlusResultPassesSynFilter(searchResultsCurrentScan[0]);
                        }
                        else
                        {
                            // First Hits file
                            var scanChargeKey = searchResultsCurrentScan[0].Scan + "_" + searchResultsCurrentScan[0].Charge;

                            if (scanChargeFirstHit.TryGetValue(scanChargeKey, out var firstHitPeptide))
                            {
                                // A result has already been stored for this scan/charge combo
                                validSearchResult = false;

                                // Possibly update the associated protein name
                                if (firstHitPeptide.CleanSequence.Equals(
                                    GetCleanSequence(searchResultsCurrentScan[0].Peptide, out var newPrefix, out var newSuffix)))
                                {
                                    var bestProtein = GetBestProteinName(
                                        firstHitPeptide.ProteinName,
                                        firstHitPeptide.ProteinNumber,
                                        searchResultsCurrentScan[0].Protein);

                                    if (bestProtein.Value < firstHitPeptide.ProteinNumber)
                                    {
                                        firstHitPeptide.ProteinName = bestProtein.Key;
                                        firstHitPeptide.ProteinNumber = bestProtein.Value;
                                        firstHitPeptide.UpdatePrefixAndSuffix(newPrefix, newSuffix);
                                    }
                                }
                            }
                            else
                            {
                                firstHitPeptide = new FirstHitInfo(searchResultsCurrentScan[0].Peptide,
                                    GetCleanSequence(searchResultsCurrentScan[0].Peptide))
                                {
                                    ProteinName = searchResultsCurrentScan[0].Protein,
                                    ProteinNumber = int.MaxValue
                                };

                                if (mProteinNameOrder.TryGetValue(searchResultsCurrentScan[0].Protein, out var proteinNumber))
                                {
                                    firstHitPeptide.ProteinNumber = proteinNumber;
                                }

                                scanChargeFirstHit.Add(scanChargeKey, firstHitPeptide);
                            }
                        }

                        if (validSearchResult)
                        {
                            ExpandListIfRequired(searchResultsPrefiltered, searchResultsCurrentScan.Count);
                            searchResultsPrefiltered.AddRange(searchResultsCurrentScan);
                        }

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    searchResultsPrefiltered.TrimExcess();

                    // Sort the SearchResults by scan, charge, and ascending SpecEValue
                    searchResultsPrefiltered.Sort(new MSGFPlusSearchResultsComparerScanChargeSpecEValuePeptide());

                    if (filteredOutputFileType == FilteredOutputFileTypeConstants.FHTFile)
                    {
                        // Update the protein names in searchResultsPrefiltered using scanChargeFirstHit
                        // This step is typically not necessary, but is often required for SplitFasta results

                        for (var index = 0; index <= searchResultsPrefiltered.Count - 1; index++)
                        {
                            var scanChargeKey = searchResultsPrefiltered[index].Scan + "_" + searchResultsPrefiltered[index].Charge;

                            if (!scanChargeFirstHit.TryGetValue(scanChargeKey, out var firstHitPeptide))
                                continue;

                            if (searchResultsPrefiltered[index].Protein.Equals(firstHitPeptide.ProteinName))
                                continue;

                            // Protein name doesn't match the expected name
                            if (firstHitPeptide.CleanSequence.Equals(GetCleanSequence(searchResultsPrefiltered[index].Peptide)))
                            {
                                // Update the protein name
                                var updatedSearchResult = searchResultsPrefiltered[index];
                                updatedSearchResult.Peptide = firstHitPeptide.SequenceWithModsAndContext;
                                updatedSearchResult.Protein = firstHitPeptide.ProteinName;
                                searchResultsPrefiltered[index] = updatedSearchResult;
                            }
                            else
                            {
                                Console.WriteLine("Possible programming bug; " +
                                                  "mix of peptides tracked for a given scan/charge combo when caching data for First Hits files; " +
                                                  $"see scan_charge {scanChargeKey}");
                            }
                        }
                    }

                    // Now filter the data and store in filteredSearchResults
                    // Due to code updates in October 2016, searchResultsPrefiltered already has filtered data
                    var startIndex = 0;

                    while (startIndex < searchResultsPrefiltered.Count)
                    {
                        var endIndex = startIndex;
                        // Find all of the peptides with the same scan number
                        while (endIndex + 1 < searchResultsPrefiltered.Count &&
                               searchResultsPrefiltered[endIndex + 1].ScanNum == searchResultsPrefiltered[startIndex].ScanNum)
                        {
                            endIndex++;
                        }

                        // Store the results for this scan
                        if (filteredOutputFileType == FilteredOutputFileTypeConstants.SynFile)
                        {
                            StoreSynMatches(searchResultsPrefiltered, startIndex, endIndex, filteredSearchResults);
                        }
                        else
                        {
                            StoreTopFHTMatch(searchResultsPrefiltered, startIndex, endIndex, filteredSearchResults);
                        }

                        startIndex = endIndex + 1;
                    }

                    // Sort the data in filteredSearchResults then write out to disk
                    SortAndWriteFilteredSearchResults(writer, filteredSearchResults, errorMessages,
                        includeFDRandPepFDR, includeEFDR, includeIMSFields, isMsgfPlus);

                    // Write out the scan group info
                    if (!string.IsNullOrEmpty(scanGroupFilePath))
                    {
                        StoreScanGroupInfo(scanGroupFilePath, scanGroupDetails);
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
                    SetErrorMessage("Error reading input file in CreateFHTorSYNResultsFile", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                    return false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating the output file in CreateFHTorSYNResultsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Extracts mod info from either a MS-GF+ param file or from a MSGFPlus_Mods.txt file (previously MSGFDB_Mods.txt)
        /// </summary>
        /// <param name="msgfPlusParamFilePath"></param>
        /// <param name="modList"></param>
        /// <returns>True if success; false if a problem</returns>
        private bool ExtractModInfoFromParamFile(
            string msgfPlusParamFilePath,
            out List<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            var modFileProcessor = new MSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                msgfPlusParamFilePath,
                MSGFPlusParamFileModExtractor.ModSpecFormats.MSGFPlusAndMSPathFinder,
                out modList);

            if (!success || mErrorCode != PHRPErrorCode.NoError)
            {
                if (mErrorCode == PHRPErrorCode.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MS-GF+ parameter file");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modList, mPeptideMods);

            return true;
        }

        /// <summary>
        /// Extracts parent mass tolerance from the parameters loaded from an MS-GF+ parameter file
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <returns>Parent mass tolerance info.  Tolerances will be 0 if an error occurs</returns>
        private PrecursorMassTolerance ExtractParentMassToleranceFromParamFile(SearchEngineParameters searchEngineParams)
        {
            try
            {
                if (!searchEngineParams.Parameters.TryGetValue(MSGFPlusSynFileReader.PRECURSOR_TOLERANCE_PARAM_NAME, out var value) &&
                    !searchEngineParams.Parameters.TryGetValue(MSGFPlusSynFileReader.PRECURSOR_TOLERANCE_PARAM_NAME_SYNONYM, out value))
                {
                    OnWarningEvent(
                        "Could not find parameter {0} or {1} in parameter file {2}; cannot determine the precursor mass tolerance",
                        MSGFPlusSynFileReader.PRECURSOR_TOLERANCE_PARAM_NAME,
                        MSGFPlusSynFileReader.PRECURSOR_TOLERANCE_PARAM_NAME_SYNONYM,
                        Path.GetFileName(searchEngineParams.SearchEngineParamFilePath));

                    return new PrecursorMassTolerance(75, true);
                }

                // Parent ion tolerance line found

                // Split the line on commas
                var splitLine = value.Split(',');

                // ReSharper disable once ConvertIfStatementToSwitchStatement
                if (splitLine.Length == 1 && ParseParentMassTolerance(splitLine[0], out var tolerance, out var isPPM))
                {
                    // Tolerance does not contain a comma
                    return new PrecursorMassTolerance(tolerance, isPPM);
                }

                if (splitLine.Length <= 1)
                {
                    // Tolerance does not contain a comma and does not contain a valid number
                    OnWarningEvent("Invalid parent ion tolerance format: {0}; should be 20ppm or 1.5Da or 20ppm,30ppm", value);
                    return new PrecursorMassTolerance(75, true);
                }

                // Tolerance line does contain a comma

                if (ParseParentMassTolerance(splitLine[0], out var leftTolerance, out var leftToleranceIsPPM))
                {
                    // ReSharper disable once ConvertIfStatementToReturnStatement
                    if (ParseParentMassTolerance(splitLine[1], out var rightTolerance, out _))
                    {
                        return new PrecursorMassTolerance(leftTolerance, rightTolerance, leftToleranceIsPPM);
                    }

                    return new PrecursorMassTolerance(leftTolerance, leftTolerance, leftToleranceIsPPM);
                }

                OnWarningEvent("Invalid parent ion tolerance format: {0}; should be 0.5Da,2.5Da", value);
                return new PrecursorMassTolerance(75, true);
            }
            catch (Exception ex)
            {
                SetErrorMessage(string.Format("Error parsing the ParentMass tolerance from the MS-GF+ parameter file ({0}): {1}",
                    Path.GetFileName(searchEngineParams.SearchEngineParamFilePath), ex.Message), ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                return new PrecursorMassTolerance(75, true);
            }
        }

        /// <summary>
        /// Look for candidateProteinName in mProteinNameOrder
        /// If found, and if the proteinNumber for that protein is less than currentProteinNumber,
        /// return a KeyValuePair with candidateProteinName and the proteinNumber for that protein
        /// Otherwise, return a KeyValuePair with currentProteinName and currentProteinNumber
        /// </summary>
        /// <param name="currentProteinName"></param>
        /// <param name="currentProteinNumber"></param>
        /// <param name="candidateProteinName"></param>
        /// <returns>KeyValuePair of the best protein name</returns>
        private KeyValuePair<string, int> GetBestProteinName(string currentProteinName, int currentProteinNumber, string candidateProteinName)
        {
            if (mProteinNameOrder.Count == 0)
            {
                return new KeyValuePair<string, int>(currentProteinName, currentProteinNumber);
            }

            // Lookup the protein number (to make sure we use the protein name that occurs first in the FASTA file)
            // Only possible if the user provided the path to the FASTA file

            if (mProteinNameOrder.TryGetValue(candidateProteinName, out var proteinNumber))
            {
                if (proteinNumber < currentProteinNumber)
                {
                    // A better protein name has been found (or this is the first protein and we just determined the protein number to associate with it)
                    return new KeyValuePair<string, int>(candidateProteinName, proteinNumber);
                }
            }
            // ReSharper disable once RedundantIfElseBlock
            else
            {
                // Protein not found in mProteinNameOrder
                // It's likely a reverse-hit protein
            }

            // A better protein name was not found; return the current info
            return new KeyValuePair<string, int>(currentProteinName, currentProteinNumber);
        }

        /// <summary>
        /// Load the PeptideToProteinMap information; in addition, creates the _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected termini symbols
        /// </summary>
        /// <param name="pepToProteinMapFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="msgfPlusModInfo"></param>
        /// <param name="isMsgfPlus">Should be set to True if processing MS-GF+ results</param>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <returns>True if successful, false if an error</returns>
        private void LoadPeptideToProteinMapInfoMSGFPlus(
            string pepToProteinMapFilePath,
            string outputDirectoryPath,
            IReadOnlyList<MSGFPlusParamFileModExtractor.ModInfo> msgfPlusModInfo,
            bool isMsgfPlus,
            List<PepToProteinMapping> pepToProteinMapping,
            out string mtsPepToProteinMapFilePath)
        {
            mtsPepToProteinMapFilePath = string.Empty;

            try
            {
                if (string.IsNullOrWhiteSpace(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    OnWarningEvent("PepToProteinMap file is not defined");
                    return;
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    var pepToProteinMapAlternate = ReaderFactory.AutoSwitchToLegacyMSGFDBIfRequired(pepToProteinMapFilePath, "Dataset_msgfdb.txt");
                    if (File.Exists(pepToProteinMapAlternate))
                    {
                        pepToProteinMapFilePath = pepToProteinMapAlternate;
                    }
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    // The analysis manager creates file Dataset_msgfplus_PepToProtMap.txt after running MS-GF+
                    // See MSGFPlusUtils.CreatePeptideToProteinMapping

                    Console.WriteLine();

                    var fastaFilePath = FindInputFile(Options.FastaFilePath, Options.AlternateBasePath, out var fastaFile, true)
                        ? fastaFile.FullName
                        : Options.FastaFilePath;

                    // If the .FASTA file was defined (and exists), this program will use it to define the peptide to protein mapping
                    if (!string.IsNullOrWhiteSpace(fastaFilePath) && File.Exists(fastaFilePath))
                    {
                        OnDebugEvent(
                            "The PepToProteinMap file does not exist ({0}), but the FASTA file does exist ({1}); " +
                            "it will be used to determine the peptide to protein mapping",
                            Path.GetFileName(pepToProteinMapFilePath),
                            Path.GetFileName(Options.FastaFilePath));
                    }
                    else
                    {
                        OnWarningEvent("PepToProteinMap file does not exist: " + pepToProteinMapFilePath);
                        OnStatusEvent(
                            "The PepToProteinMap file is typically created by the Protein Coverage Summarizer. " +
                            "Proteins associated with each peptide will be based on protein names in the input file");
                    }

                    return;
                }

                // Initialize pepToProteinMapping
                pepToProteinMapping.Clear();

                // Read the data in the peptide to protein map file
                var success = LoadPeptideToProteinMapInfo(pepToProteinMapFilePath, pepToProteinMapping, out var headerLine);

                if (!success)
                {
                    return;
                }

                mtsPepToProteinMapFilePath =
                    Path.Combine(outputDirectoryPath, Path.GetFileNameWithoutExtension(pepToProteinMapFilePath) + "MTS.txt");

                using var writer = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                if (!string.IsNullOrEmpty(headerLine))
                {
                    // Header line
                    writer.WriteLine(headerLine);
                }

                for (var index = 0; index <= pepToProteinMapping.Count - 1; index++)
                {
                    // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                    // In addition, replace the * terminus symbols with dashes
                    var mtsCompatiblePeptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(pepToProteinMapping[index].Peptide),
                        msgfPlusModInfo, isMsgfPlus, out _);

                    if (pepToProteinMapping[index].Peptide != mtsCompatiblePeptide)
                    {
                        UpdatePepToProteinMapPeptide(pepToProteinMapping, index, mtsCompatiblePeptide);
                    }

                    writer.WriteLine(pepToProteinMapping[index].Peptide + "\t" +
                                     pepToProteinMapping[index].Protein + "\t" +
                                     pepToProteinMapping[index].ResidueStart + "\t" +
                                     pepToProteinMapping[index].ResidueEnd);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" + Path.GetFileName(mtsPepToProteinMapFilePath) + "): " + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
            }
        }

        /// <summary>
        /// Load the MS-GF+ parameter file and updates settings
        /// </summary>
        /// <param name="msgfPlusParamFilePath"></param>
        /// <returns>
        /// True if success, false if an error.
        /// Returns True if msgfPlusParamFilePath is empty
        /// Returns False if the paramFilePath is defined but the file is not found or cannot be parsed</returns>
        private bool LoadSearchEngineParamFile(string msgfPlusParamFilePath)
        {
            if (string.IsNullOrWhiteSpace(msgfPlusParamFilePath))
            {
                OnWarningEvent("MS-GF+ parameter file is not defined. Unable to extract parent mass tolerance info or custom charge carrier masses");
                return true;
            }

            var searchEngineParams = new SearchEngineParameters(SEARCH_ENGINE_NAME);

            var success = SynFileReaderBaseClass.ReadKeyValuePairSearchEngineParamFile(
                SEARCH_ENGINE_NAME, msgfPlusParamFilePath, PeptideHitResultTypes.MSGFPlus,
                searchEngineParams, out var localErrorMessage, out var localWarningMessage);

            if (!string.IsNullOrWhiteSpace(localErrorMessage))
            {
                ReportError(localErrorMessage);
                return false;
            }

            if (!string.IsNullOrWhiteSpace(localWarningMessage))
            {
                OnWarningEvent(localWarningMessage);
            }

            if (searchEngineParams.Parameters.Count == 0)
            {
                SetErrorMessage("MS-GF+ parameter file is empty; unable to extract parent mass tolerance info");
                SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                return false;
            }

            // Parse the PrecursorMassTolerance setting
            mPrecursorMassTolerance = ExtractParentMassToleranceFromParamFile(searchEngineParams);

            // Parse the ChargeCarrierMass setting
            if (MSGFPlusSynFileReader.GetCustomChargeCarrierMass(searchEngineParams, out var customChargeCarrierMass))
            {
                OnStatusEvent("Using a charge carrier mass of {0:F3} Da", customChargeCarrierMass);
                mPeptideSeqMassCalculator.ChargeCarrierMass = customChargeCarrierMass;
            }

            return success;
        }

        private void ModExtractorErrorHandler(string message, Exception ex)
        {
            SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
        }

        private bool MSGFPlusResultPassesSynFilter(MSGFPlusSearchResult msgfPlusSearchResultType)
        {
            return msgfPlusSearchResultType.EValueNum <= Options.MSGFPlusSynopsisFileEValueThreshold ||
                   msgfPlusSearchResultType.SpecEValueNum <= Options.MSGFPlusSynopsisFileSpecEValueThreshold ||
                   msgfPlusSearchResultType.QValueNum is > 0 and < 0.01;
        }

        private bool ParseMSGFPlusSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            List<PepToProteinMapping> pepToProteinMapping,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            var successOverall = true;

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (resetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                mNumericModErrors = 0;

                // Initialize searchResult
                var searchResult = new MSGFPlusResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Note that MS-GF+ synopsis files are normally sorted on SpecEValue, ascending
                // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
                // we will keep track of the scan, charge, and peptide information parsed for each unique SpecEValue encountered
                // (see peptidesFoundForSpecEValueLevel below)

                // This is required since a PSM with multiple proteins will be listed on multiple lines in the synopsis file
                // Values are PeptideSequenceWithMods_Scan_Charge

                var peptidesFoundForSpecEValueLevel = new SortedSet<string>();

                var previousSpecEValue = string.Empty;

                // Assure that pepToProteinMapping is sorted on peptide
                if (pepToProteinMapping.Count > 1)
                {
                    pepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(Options);

                    var errorMessages = new List<string>();

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

                    var columnMapping = new Dictionary<MSGFPlusSynFileColumns, int>();

                    var peptidesNotFoundInPepToProtMapping = 0;

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
                            var validHeader = ParseMSGFPlusSynFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }
                            headerParsed = true;
                            continue;
                        }

                        var validSearchResult = ParseMSGFPlusSynFileEntry(
                            lineIn, searchResult, errorMessages,
                            resultsProcessed, columnMapping,
                            out var currentPeptideWithMods);

                        resultsProcessed++;
                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                        bool firstMatchForGroup;

                        if (searchResult.SpecEValue == previousSpecEValue)
                        {
                            // New result has the same SpecEValue as the previous result
                            // See if peptidesFoundForSpecEValueLevel contains the peptide, scan and charge

                            if (peptidesFoundForSpecEValueLevel.Contains(key))
                            {
                                firstMatchForGroup = false;
                            }
                            else
                            {
                                peptidesFoundForSpecEValueLevel.Add(key);
                                firstMatchForGroup = true;
                            }
                        }
                        else
                        {
                            // New SpecEValue
                            peptidesFoundForSpecEValueLevel.Clear();

                            // Update previousSpecEValue
                            previousSpecEValue = searchResult.SpecEValue;

                            // Append a new entry
                            peptidesFoundForSpecEValueLevel.Add(key);
                            firstMatchForGroup = true;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
                            successOverall = false;
                            errorMessages.Add(string.Format("Error adding modifications to sequence for ResultID '{0}'", searchResult.ResultID));
                        }

                        SaveResultsFileEntrySeqInfo(searchResult, firstMatchForGroup);

                        if (pepToProteinMapping.Count > 0)
                        {
                            // Add the additional proteins for this peptide

                            // Use binary search to find this peptide in pepToProteinMapping
                            var pepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(pepToProteinMapping, currentPeptideWithMods);

                            if (pepToProteinMapIndex >= 0)
                            {
                                // Call MyBase.SaveResultsFileEntrySeqInfo for each entry in pepToProteinMapping() for peptide , skipping searchResult.ProteinName
                                var currentProtein = searchResult.ProteinName;
                                do
                                {
                                    if (pepToProteinMapping[pepToProteinMapIndex].Protein != currentProtein)
                                    {
                                        searchResult.ProteinName = pepToProteinMapping[pepToProteinMapIndex].Protein;
                                        SaveResultsFileEntrySeqInfo(searchResult, false);
                                    }

                                    pepToProteinMapIndex++;
                                } while (pepToProteinMapIndex < pepToProteinMapping.Count &&
                                         currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                            }
                            else
                            {
                                // Match not found; this is unexpected
                                peptidesNotFoundInPepToProtMapping++;
                                ShowPeriodicWarning(peptidesNotFoundInPepToProtMapping,
                                    10,
                                    "no match for '" + currentPeptideWithMods + "' in pepToProteinMapping");
                            }
                        }

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    if (Options.CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(inputFilePath);
                        var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                        SaveModificationSummaryFile(modificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (errorMessages.Count > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                    }

                    return successOverall;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in ParseMSGFPlusSynopsisFile", ex);
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
                SetErrorMessage("Error creating the output file in ParseMSGFPlusSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        private bool ParseMSGFPlusResultsFileEntry(
            string lineIn,
            bool isMsgfPlus,
            IReadOnlyList<MSGFPlusParamFileModExtractor.ModInfo> msgfPlusModInfo,
            ICollection<MSGFPlusSearchResult> searchResultsCurrentScan,
            ICollection<string> errorMessages,
            IDictionary<MSGFPlusResultsFileColumns, int> columnMapping,
            ref int nextScanGroupID,
            ICollection<ScanGroupInfo> scanGroupDetails,
            IDictionary<string, bool> scanGroupCombo,
            IDictionary<string, int> specIdToIndex)
        {
            // Parses an entry from the MS-GF+ results file

            var udtSearchResult = new MSGFPlusSearchResult();
            string rowIndex = null;

            MSGFPlusSearchResult[] udtMergedScanInfo = null;

            try
            {
                var proteinInfo = new Dictionary<string, TerminusChars>();

                // Reset searchResults
                searchResultsCurrentScan.Clear();

                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 13)
                {
                    // Not a valid result
                    return false;
                }

                rowIndex = splitLine[0];

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.SpectrumFile], out udtSearchResult.SpectrumFileName))
                {
                    ReportError("SpectrumFile column is missing or invalid", true);
                }
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.SpecIndex], out udtSearchResult.SpecIndex);

                if (isMsgfPlus)
                {
                    var generateSpecIndex = true;

                    if (!int.TryParse(udtSearchResult.SpecIndex, out var specIndex))
                    {
                        // MS-GF+ includes text in the SpecID column, for example: "controllerType=0 controllerNumber=1 scan=6390" or "index=4323"
                        // Need to convert these to an integer

                        if (udtSearchResult.SpecIndex.StartsWith("index="))
                        {
                            udtSearchResult.SpecIndex = udtSearchResult.SpecIndex.Substring("index=".Length);
                            if (int.TryParse(udtSearchResult.SpecIndex, out specIndex))
                            {
                                generateSpecIndex = false;
                            }
                        }

                        if (generateSpecIndex)
                        {
                            if (!specIdToIndex.TryGetValue(udtSearchResult.SpecIndex, out specIndex))
                            {
                                specIndex = specIdToIndex.Count + 1;
                                specIdToIndex.Add(udtSearchResult.SpecIndex, specIndex);
                            }

                            udtSearchResult.SpecIndex = specIndex.ToString();
                        }
                    }
                }

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.Scan], out udtSearchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid", true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                var slashIndex = udtSearchResult.Scan.IndexOf('/');
                int scanCount;
                if (slashIndex > 0)
                {
                    // This is a merged spectrum and thus scan number looks like: 3010/3011/3012
                    // Split the Scan list on the slash
                    // Later in this method, we'll append searchResults with this scan plus the other scans

                    var splitResult = udtSearchResult.Scan.Split('/');
                    scanCount = splitResult.Length;
                    udtMergedScanInfo = new MSGFPlusSearchResult[scanCount];

                    for (var index = 0; index <= scanCount - 1; index++)
                    {
                        udtMergedScanInfo[index] = new MSGFPlusSearchResult();
                        udtMergedScanInfo[index].Clear();
                        udtMergedScanInfo[index].Scan = splitResult[index];
                        udtMergedScanInfo[index].ScanNum = StringUtilities.CIntSafe(splitResult[index], 0);
                    }

                    // Now split SpecIndex and store in udtMergedScanInfo
                    splitResult = udtSearchResult.SpecIndex.Split('/');

                    for (var index = 0; index <= splitResult.Length - 1; index++)
                    {
                        if (index >= udtMergedScanInfo.Length)
                        {
                            // There are more entries for SpecIndex than there are for Scan#; this is unexpected
                            break;
                        }
                        udtMergedScanInfo[index].SpecIndex = splitResult[index];
                    }

                    // Now split FragMethod and store in udtMergedScanInfo
                    splitResult = udtSearchResult.FragMethod.Split('/');

                    for (var index = 0; index <= splitResult.Length - 1; index++)
                    {
                        if (index >= udtMergedScanInfo.Length)
                        {
                            // There are more entries for FragMethod than there are for Scan#; this is unexpected
                            break;
                        }
                        udtMergedScanInfo[index].FragMethod = splitResult[index];
                    }
                }
                else
                {
                    udtSearchResult.ScanNum = StringUtilities.CIntSafe(udtSearchResult.Scan, 0);
                    scanCount = 1;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.PrecursorMZ], out udtSearchResult.PrecursorMZ);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(udtSearchResult.Charge, 0));

                // Precursor mass error could be in PPM or Da
                //   In MSGFDB, the header line will have PMError(ppm)        or PMError(Da)
                //   In MS-GF+, the header line will have PrecursorError(ppm) or PrecursorError(Da)
                double precursorErrorDa = 0;

                if (columnMapping[MSGFPlusResultsFileColumns.PMErrorPPM] >= 0)
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.PMErrorPPM], out udtSearchResult.PMErrorPPM);
                }
                else
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.PMErrorDa], out udtSearchResult.PMErrorDa);
                    precursorErrorDa = StringUtilities.CDblSafe(udtSearchResult.PMErrorDa, 0);
                    udtSearchResult.PMErrorPPM = string.Empty; // We'll populate this column later in this method
                }

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.Protein], out udtSearchResult.Protein);

                // MS-GF+ .tsv files may have a semicolon separated list of protein names; check for this
                udtSearchResult.Protein = SplitProteinList(udtSearchResult.Protein, proteinInfo);

                if (proteinInfo.Count > 0)
                {
                    // Need to add the prefix and suffix residues
                    udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, proteinInfo.First());
                }

                // Replace any mod text values in the peptide sequence with the appropriate mod symbols
                // In addition, replace the terminus symbols with dashes
                udtSearchResult.Peptide =
                    ReplaceMSGFModTextWithSymbol(ReplaceTerminus(udtSearchResult.Peptide), msgfPlusModInfo, isMsgfPlus,
                                                 out var totalModMass);

                // Compute monoisotopic mass of the peptide
                var peptideMonoisotopicMass = ComputePeptideMass(udtSearchResult.Peptide, totalModMass);

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoisotopicMass, 0), 6);

                if (!string.IsNullOrEmpty(udtSearchResult.PMErrorPPM))
                {
                    // Convert the ppm-based PM Error to Da-based

                    if (double.TryParse(udtSearchResult.PrecursorMZ, out var precursorMZ))
                    {
                        // Note that since .PMErrorPPM is present, the Precursor m/z is a C13-corrected m/z value
                        // In other words, it may not be the actual m/z selected for fragmentation.

                        if (double.TryParse(udtSearchResult.PMErrorPPM, out var parentMassErrorPPM))
                        {
                            if (mPrecursorMassTolerance.IsPPM &&
                                (parentMassErrorPPM < -mPrecursorMassTolerance.ToleranceLeft * 1.5 ||
                                 parentMassErrorPPM > mPrecursorMassTolerance.ToleranceRight * 1.5))
                            {
                                // PPM error computed by MS-GF+ is more than 1.5-fold larger than the ppm-based parent ion tolerance; don't trust the value computed by MS-GF+
                                mPrecursorMassErrorWarningCount++;
                                ShowPeriodicWarning(mPrecursorMassErrorWarningCount,
                                                    10,
                                                    string.Format("Precursor mass error computed by MS-GF+ is 1.5-fold larger than the search tolerance: {0} vs. {1:F0}ppm,{2:F0}ppm",
                                                    udtSearchResult.PMErrorPPM,
                                                    mPrecursorMassTolerance.ToleranceLeft,
                                                    mPrecursorMassTolerance.ToleranceRight));

                                var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, udtSearchResult.ChargeNum, 0);

                                precursorErrorDa = precursorMonoMass - peptideMonoisotopicMass;

                                udtSearchResult.PMErrorPPM = string.Empty;
                            }
                            else
                            {
                                precursorErrorDa = PeptideMassCalculator.PPMToMass(parentMassErrorPPM, peptideMonoisotopicMass);

                                // Note that this will be a C13-corrected precursor error; not the absolute precursor error
                                udtSearchResult.PMErrorDa = StringUtilities.MassErrorToString(precursorErrorDa);
                            }
                        }
                    }
                }

                if (string.IsNullOrEmpty(udtSearchResult.PMErrorPPM))
                {
                    if (double.TryParse(udtSearchResult.PrecursorMZ, out var precursorMZ))
                    {
                        var peptideDeltaMassCorrectedPpm =
                            ComputeDelMCorrectedPPM(precursorErrorDa, precursorMZ, udtSearchResult.ChargeNum, peptideMonoisotopicMass, true);

                        udtSearchResult.PMErrorPPM = PRISM.StringUtilities.DblToString(peptideDeltaMassCorrectedPpm, 5, 0.00005);

                        if (string.IsNullOrEmpty(udtSearchResult.PMErrorDa))
                        {
                            precursorErrorDa = PeptideMassCalculator.PPMToMass(peptideDeltaMassCorrectedPpm, peptideMonoisotopicMass);

                            // Note that this will be a C13-corrected precursor error; not the absolute precursor error
                            udtSearchResult.PMErrorDa = StringUtilities.MassErrorToString(precursorErrorDa);
                        }
                    }
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.DeNovoScore], out udtSearchResult.DeNovoScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.MSGFScore], out udtSearchResult.MSGFScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.SpecProb_EValue], out udtSearchResult.SpecEValue);
                if (!double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum))
                    udtSearchResult.SpecEValueNum = 0;

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.PValue_EValue], out udtSearchResult.EValue);
                if (!double.TryParse(udtSearchResult.EValue, out udtSearchResult.EValueNum))
                    udtSearchResult.EValueNum = 0;

                var targetDecoyFDRValid = DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.FDR_QValue], out udtSearchResult.QValue);
                if (!double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum))
                    udtSearchResult.QValueNum = 0;

                if (targetDecoyFDRValid)
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.PepFDR_PepQValue], out udtSearchResult.PepQValue);
                }
                else
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.EFDR], out udtSearchResult.QValue);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.IsotopeError], out udtSearchResult.IsotopeError);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.IMSScan], out udtSearchResult.IMSScan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusResultsFileColumns.IMSDriftTime], out udtSearchResult.IMSDriftTime);

                udtSearchResult.NTT = ComputeCleavageState(udtSearchResult.Peptide).ToString();

                var udtScanGroupInfo = new ScanGroupInfo();
                var currentScanGroupID = -1;

                udtScanGroupInfo.Charge = udtSearchResult.ChargeNum;

                if (scanCount > 1)
                {
                    // This result came from a merged spectrum and thus has a scan number formatted like: 3010/3011/3012
                    // Append one entry to searchResults for each item in udtMergedScanInfo()

                    for (var index = 0; index <= scanCount - 1; index++)
                    {
                        udtSearchResult.Scan = udtMergedScanInfo[index].Scan;
                        udtSearchResult.ScanNum = udtMergedScanInfo[index].ScanNum;

                        udtSearchResult.SpecIndex = udtMergedScanInfo[index].SpecIndex;
                        udtSearchResult.FragMethod = udtMergedScanInfo[index].FragMethod;

                        AppendToSearchResults(searchResultsCurrentScan, udtSearchResult, proteinInfo);

                        // Append an entry to scanGroupDetails
                        udtScanGroupInfo.Scan = udtSearchResult.ScanNum;
                        AppendToScanGroupDetails(scanGroupDetails, scanGroupCombo, udtScanGroupInfo, ref currentScanGroupID,
                                                 ref nextScanGroupID);
                    }
                }
                else
                {
                    // This is not a merged result; simply append udtSearchResult to searchResults
                    AppendToSearchResults(searchResultsCurrentScan, udtSearchResult, proteinInfo);

                    // Also append an entry to scanGroupDetails
                    udtScanGroupInfo.Scan = udtSearchResult.ScanNum;
                    AppendToScanGroupDetails(scanGroupDetails, scanGroupCombo, udtScanGroupInfo, ref currentScanGroupID,
                                             ref nextScanGroupID);
                }

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the MS-GF+ results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing MS-GF+ results in ParseMSGFPlusResultsFileEntry for RowIndex '{0}': {1}", rowIndex, ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing MS-GF+ Results in ParseMSGFPlusResultsFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the MS-GF+ results file header line, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSGFPlusResultsFileHeaderLine(string lineIn, IDictionary<MSGFPlusResultsFileColumns, int> columnMapping)
        {
            // The expected header from MSGFDB is:
            // #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(Da)           Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR
            // or
            // #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(ppm)          Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR

            // The expected header from MS-GF+ is:
            // #SpecFile    SpecID       ScanNum   ScanTime(Min)    FragMethod    Precursor    IsotopeError    PrecursorError(Da)    Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue
            // or
            // #SpecFile    SpecID       ScanNum   ScanTime(Min)    FragMethod    Precursor    IsotopeError    PrecursorError(ppm)   Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue

            var columnNames = new SortedDictionary<string, MSGFPlusResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"#SpecFile", MSGFPlusResultsFileColumns.SpectrumFile},
                {"SpecIndex", MSGFPlusResultsFileColumns.SpecIndex},
                {"SpecID", MSGFPlusResultsFileColumns.SpecIndex},
                {"Scan#", MSGFPlusResultsFileColumns.Scan},
                {"ScanNum", MSGFPlusResultsFileColumns.Scan},
                {"ScanTime(Min)", MSGFPlusResultsFileColumns.ScanTimeMinutes},
                {"FragMethod", MSGFPlusResultsFileColumns.FragMethod},
                {"Precursor", MSGFPlusResultsFileColumns.PrecursorMZ},
                {"IsotopeError", MSGFPlusResultsFileColumns.IsotopeError},
                {"PMError(Da)", MSGFPlusResultsFileColumns.PMErrorDa},
                {"PrecursorError(Da)", MSGFPlusResultsFileColumns.PMErrorDa},
                {"PMError(ppm)", MSGFPlusResultsFileColumns.PMErrorPPM},
                {"PrecursorError(ppm)", MSGFPlusResultsFileColumns.PMErrorPPM},
                {"Charge", MSGFPlusResultsFileColumns.Charge},
                {"Peptide", MSGFPlusResultsFileColumns.Peptide},
                {"Protein", MSGFPlusResultsFileColumns.Protein},
                {"DeNovoScore", MSGFPlusResultsFileColumns.DeNovoScore},
                {"MSGFScore", MSGFPlusResultsFileColumns.MSGFScore},
                {"SpecProb", MSGFPlusResultsFileColumns.SpecProb_EValue},
                {"SpecEValue", MSGFPlusResultsFileColumns.SpecProb_EValue},
                {"P-value", MSGFPlusResultsFileColumns.PValue_EValue},
                {"EValue", MSGFPlusResultsFileColumns.PValue_EValue},
                {"FDR", MSGFPlusResultsFileColumns.FDR_QValue},
                {"QValue", MSGFPlusResultsFileColumns.FDR_QValue},
                {"PepFDR", MSGFPlusResultsFileColumns.PepFDR_PepQValue},
                {"PepQValue", MSGFPlusResultsFileColumns.PepFDR_PepQValue},
                {"EFDR", MSGFPlusResultsFileColumns.EFDR},
                {"IMS_Scan", MSGFPlusResultsFileColumns.IMSScan},
                {"IMS_Drift_Time", MSGFPlusResultsFileColumns.IMSDriftTime}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MSGFPlusResultsFileColumns resultColumn in Enum.GetValues(typeof(MSGFPlusResultsFileColumns)))
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
                    else
                    {
                        // Unrecognized column name
                        Console.WriteLine("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMSGFPlusResultsFileHeaderLine");
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the MSGFPlus results file: " + ex.Message, ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse the header line of a MS-GF+ _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSGFPlusSynFileHeaderLine(string lineIn, IDictionary<MSGFPlusSynFileColumns, int> columnMapping)
        {
            var columnNames = MSGFPlusSynFileReader.GetColumnHeaderNamesAndIDs(true, true);

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MSGFPlusSynFileColumns resultColumn in Enum.GetValues(typeof(MSGFPlusSynFileColumns)))
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
                SetErrorMessage("Error parsing the header line in the MSGFPlus synopsis file: " + ex.Message, ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse an entry from a MS-GF+ synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequenceWithMods"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSGFPlusSynFileEntry(
            string lineIn,
            MSGFPlusResults searchResult,
            ICollection<string> errorMessages,
            int resultsProcessed,
            IDictionary<MSGFPlusSynFileColumns, int> columnMapping,
            out string peptideSequenceWithMods)
        {
            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();
            peptideSequenceWithMods = string.Empty;

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 15)
                {
                    return false;
                }

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.ResultID], out string value))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading ResultID value from MS-GF+ results, line {0}", resultsProcessed + 1));
                    }

                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.Scan], out string scan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading peptide sequence from MS-GF+ results, line {0}", resultsProcessed + 1));
                    }

                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.DelM], out string msgfPlusComputedDelM);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.DelMPPM], out string msgfPlusComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.MSGFPlusComputedDelM = msgfPlusComputedDelM;
                searchResult.MSGFPlusComputedDelMPPM = msgfPlusComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.MSGFPlusComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MS-GF+ results file as "Observed_Mass - Theoretical_Mass"
                // However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                // Therefore, we will negate .peptideDeltaMass
                try
                {
                    searchResult.PeptideDeltaMass = (-double.Parse(searchResult.PeptideDeltaMass)).ToString(CultureInfo.InvariantCulture);
                }
                catch (Exception)
                {
                    // Error; Leave .peptideDeltaMass unchanged
                }

                // Calling this method will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the cleavage state and terminus state
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since InSpecT only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.FragMethod], out string fragMethod);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.SpecIndex], out string specIndex);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.PrecursorMZ], out string precursorMz);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.MH], out string peptideMh);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.NTT], out string ntt);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.DeNovoScore], out string deNovoScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.MSGFScore], out string msgfScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.SpecEValue], out string specEValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.RankSpecEValue], out string rankSpecEValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.EValue], out string eValue);

                var targetDecoyFDRValid = DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.QValue], out string qValue);

                searchResult.FragMethod = fragMethod;
                searchResult.SpecIndex = specIndex;
                searchResult.PrecursorMZ = precursorMz;
                searchResult.PeptideMH = peptideMh;
                searchResult.NTT = ntt;
                searchResult.DeNovoScore = deNovoScore;
                searchResult.MSGFScore = msgfScore;
                searchResult.SpecEValue = specEValue;
                searchResult.RankSpecEValue = rankSpecEValue;
                searchResult.EValue = eValue;
                searchResult.QValue = qValue;

                if (targetDecoyFDRValid)
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.PepQValue], out string pepQValue);
                    searchResult.PepQValue = pepQValue;
                }
                else
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.EFDR], out string efdr);
                    searchResult.QValue = efdr;
                }

                if (columnMapping[MSGFPlusSynFileColumns.IsotopeError] >= 0)
                {
                    DataUtilities.GetColumnValue(splitLine, columnMapping[MSGFPlusSynFileColumns.IsotopeError], out string isotopeError);
                    searchResult.IsotopeError = isotopeError;
                    searchResult.UsedMSGFPlus = true;
                }
                else
                {
                    searchResult.UsedMSGFPlus = false;
                }

                // Compute PrecursorMH using PrecursorMZ and charge
                if (double.TryParse(searchResult.PrecursorMZ, out var precursorMZ))
                {
                    if (int.TryParse(searchResult.Charge, out var chargeValue))
                    {
                        searchResult.ParentIonMH =
                            PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, chargeValue), 6);
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count >= MAX_ERROR_MESSAGE_COUNT)
                {
                    return false;
                }

                if (splitLine?.Length > 0)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing MS-GF+ results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                }
                else
                {
                    errorMessages.Add("Error parsing MS-GF+ Results in ParseMSGFPlusSynFileEntry: " + ex.Message);
                }

                return false;
            }
        }

        private bool ParseParentMassTolerance(string toleranceText, out double tolerance, out bool isPPM)
        {
            tolerance = 0;

            toleranceText = toleranceText.ToLower().Trim();

            if (toleranceText.EndsWith("da"))
            {
                toleranceText = toleranceText.Substring(0, toleranceText.Length - 2);
                isPPM = false;
            }
            else if (toleranceText.EndsWith("ppm"))
            {
                toleranceText = toleranceText.Substring(0, toleranceText.Length - 3);
                isPPM = true;
            }
            else
            {
                isPPM = false;
                return false;
            }

            return double.TryParse(toleranceText, out tolerance);
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <remarks>Use SearchToolParameterFilePath to define the search engine parameter file</remarks>
        /// <param name="inputFilePath">MS-GF+ results file (Dataset.tsv)</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">
        /// Parameter file for data processing
        /// <para>
        /// This is an empty string if being called from the .exe and no parameter file was used or a Key=Value parameter file was provided
        /// </para>
        /// <para>
        /// If hidden command line argument /XmlParamFile was used, this will be an XML-based parameter file
        /// </para>
        /// <para>
        /// If calling this method from the DLL, this can be an empty string, a Key=Value parameter file, or an XML-based parameter file
        /// </para>
        /// </param>
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

                mPeptideSeqMassCalculator.ResetAminoAcidMasses();

                var specIdToIndex = new Dictionary<string, int>();

                ResetProgress("Parsing " + Path.GetFileName(inputFilePath));

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath, Options.AlternateBasePath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);
                    if (inputFile.DirectoryName == null)
                    {
                        OnWarningEvent("Unable to determine the parent directory of the input file: " + inputFile.FullName);
                        return false;
                    }

                    var pepToProteinMapping = new List<PepToProteinMapping>();

                    var msgfPlusParameterFilePath = ResolveFilePath(inputFile.DirectoryName, Options.SearchToolParameterFilePath);

                    // Load the MS-GF+ Parameter File so that we can determine the modification names and masses
                    // If the MSGFPlus_Mods.txt or MSGFDB_Mods.txt file was defined, the mod symbols in that file will be used to define the mod symbols in msgfPlusModInfo
                    var modInfoExtracted = ExtractModInfoFromParamFile(msgfPlusParameterFilePath, out var msgfPlusModInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    if (!LoadSearchEngineParamFile(msgfPlusParameterFilePath))
                    {
                        return false;
                    }

                    var query =
                        from item in msgfPlusModInfo
                        where item.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.CustomAA
                        select item;

                    if (query.Any())
                    {
                        // Custom amino acids are defined; read their values and update the mass calculator

                        var modFileProcessor = new MSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME);
                        RegisterEvents(modFileProcessor);

                        modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

                        MSGFPlusSynFileReader.UpdateMassCalculatorMasses(
                            msgfPlusParameterFilePath,
                            modFileProcessor,
                            mPeptideSeqMassCalculator,
                            out var localErrorMsg);

                        if (!string.IsNullOrWhiteSpace(localErrorMsg) && string.IsNullOrWhiteSpace(mErrorMessage))
                        {
                            ReportError(localErrorMsg);
                        }
                    }

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_msgfdb" with "_msgfplus"
                    if (baseName.EndsWith("_msgfdb", StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - "_msgfdb".Length) + "_msgfplus";
                    }

                    string fhtOutputFilePath;

                    if (!Options.CreateFirstHitsFile && !Options.CreateSynopsisFile)
                    {
                        OnWarningEvent("Both 'CreateFirstHitsFile' and 'CreateSynopsisFile' are false; aborting since nothing to do");
                        return true;
                    }

                    if (Options.CreateFirstHitsFile)
                    {
                        // Read the FASTA file to cache the protein names in memory
                        // These will be used when creating the first hits file
                        if (!CacheProteinNamesFromFasta())
                        {
                            return false;
                        }

                        // Create the first hits output file
                        ResetProgress("Creating the FHT file", true);

                        fhtOutputFilePath = Path.Combine(outputDirectoryPath, baseName + FIRST_HITS_FILE_SUFFIX);

                        var scanGroupFilePath = string.Empty;

                        success = CreateFHTorSYNResultsFile(
                            inputFilePath, fhtOutputFilePath, scanGroupFilePath, msgfPlusModInfo,
                            out _, specIdToIndex, FilteredOutputFileTypeConstants.FHTFile);

                        if (!success)
                        {
                            return false;
                        }
                    }
                    else
                    {
                        fhtOutputFilePath = string.Empty;
                    }

                    if (Options.CreateSynopsisFile)
                    {
                        // Create the synopsis output file
                        ResetProgress("Creating the SYN file", true);

                        // The synopsis file name will be of the form BasePath_msgfplus_syn.txt
                        var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);

                        var scanGroupFilePath = Path.Combine(outputDirectoryPath, baseName + "_ScanGroupInfo.txt");

                        success = CreateFHTorSYNResultsFile(
                            inputFilePath, synOutputFilePath, scanGroupFilePath, msgfPlusModInfo,
                            out var isMsgfPlus, specIdToIndex, FilteredOutputFileTypeConstants.SynFile);

                        if (!success)
                        {
                            return false;
                        }

                        // Load the PeptideToProteinMap information; if the file doesn't exist, a warning will be displayed, but processing will continue
                        // LoadPeptideToProteinMapInfoMSGFPlus also creates _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected termini symbols

                        var baseNameFilePath = Path.Combine(inputFile.DirectoryName ?? string.Empty, baseName);

                        var pepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseNameFilePath, outputDirectoryPath, mts: false);

                        ResetProgress("Loading the PepToProtein map file (if it exists): " + Path.GetFileName(pepToProteinMapFilePath), true);

                        LoadPeptideToProteinMapInfoMSGFPlus(
                            pepToProteinMapFilePath,
                            outputDirectoryPath,
                            msgfPlusModInfo,
                            isMsgfPlus,
                            pepToProteinMapping,
                            out var mtsPepToProteinMapFilePath);

                        // Create the other PHRP-specific files
                        ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                        // Now parse the _syn.txt file that we just created to create the other PHRP files
                        success = ParseMSGFPlusSynopsisFile(synOutputFilePath, outputDirectoryPath, pepToProteinMapping, false);
                        if (!success)
                        {
                            return false;
                        }

                        // Remove all items from pepToProteinMapping to reduce memory overhead
                        pepToProteinMapping.Clear();
                        pepToProteinMapping.TrimExcess();

                        if (Options.CreateProteinModsFile)
                        {
                            // Use a higher match error threshold since MS-GF+ often includes reverse protein peptides in the results
                            // even though the FASTA typically does not have reverse proteins
                            const int MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD = 50;

                            success = CreateProteinModsFileWork(
                                baseName, inputFile,
                                synOutputFilePath, outputDirectoryPath,
                                PeptideHitResultTypes.MSGFPlus,
                                MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD,
                                0,
                                fhtOutputFilePath,
                                mtsPepToProteinMapFilePath);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in MSGFPlusResultsProcessor.ProcessFile (2)", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in MSGFPlusResultsProcessor.ProcessFile (1)", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        private static readonly Regex NTerminalModMassMatcher = new(MSGFPlus_N_TERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);
        private static readonly Regex ModMassMatcher = new(MSGFPlus_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Replaces modification masses in peptide sequences with modification symbols (uses case-sensitive comparisons)
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="msgfPlusModInfo">This method assumes that each entry in msgfPlusModInfo has both .ModName and .ModSymbol defined</param>
        /// <param name="isMsgfPlus">Should be set to True if processing MS-GF+ results</param>
        /// <param name="totalModMass">Output parameter: total mass of all modifications</param>
        /// <returns>Updated peptide sequence</returns>
        public string ReplaceMSGFModTextWithSymbol(
            string peptide,
            IReadOnlyList<MSGFPlusParamFileModExtractor.ModInfo> msgfPlusModInfo,
            bool isMsgfPlus,
            out double totalModMass)
        {
            var prefix = string.Empty;
            var suffix = string.Empty;

            var possibleCTerminalMod = false;
            bool containsStaticMod;

            // Reset the total mod mass
            totalModMass = 0;

            // Remove the prefix and suffix residues
            if (peptide.Length >= 4 && peptide[1] == '.' && peptide[peptide.Length - 2] == '.')
            {
                prefix = peptide.Substring(0, 2);
                suffix = peptide.Substring(peptide.Length - 2, 2);

                peptide = peptide.Substring(2, peptide.Length - 4);
            }

            // Peptide should now be the primary peptide sequence, without the prefix or suffix residues

            // First look for dynamic N-terminal mods (NTermPeptide or NTermProtein)
            // This RegEx will match one or more mods, all at the N-terminus
            var match = NTerminalModMassMatcher.Match(peptide);

            if (match.Success)
            {
                // Convert the mod mass (or masses) to one or more mod symbols

                if (ConvertMSGFModMassesToSymbols("-", match.Groups[1].Value,
                    out var modSymbols, out var dynModSymbols, msgfPlusModInfo,
                    true, false,
                    out var modMassFound, out containsStaticMod))
                {
                    // Replace the mod digits with the mod symbols

                    peptide = ReplaceMSGFModTextWithMatchedSymbol(peptide, match.Groups[1], modSymbols, dynModSymbols, isMsgfPlus, containsStaticMod);
                    totalModMass += modMassFound;
                }
            }

            // Next, step through the peptide and parse each mod mass that follows a residue
            // Any mod mass at the end must be considered a C-terminal mod

            // Need to start at the first letter
            // If we had N-terminal mods, they're currently notated like this: _.+42.011MDHTPQSQLK.L or _.+42.011+57.021MNDR.Q
            // We want things to look like this: -.#MDHTPQSQLK.L or -.#*MNDRQLNHR.S

            // In MSGFDB, static mods do not have a mod mass listed
            // In MS-GF+, static mods do have a mod mass listed
            // Regardless, we do not add mod symbols for static mods, but we do increment totalModMass

            // Find the index of the last residue
            var index = peptide.Length - 1;
            while (index > 0 && !StringUtilities.IsLetterAtoZ(peptide[index]))
            {
                index--;
            }
            var indexLastResidue = index;

            // Find the index of the first residue
            index = 0;
            while (index < peptide.Length && !StringUtilities.IsLetterAtoZ(peptide[index]))
            {
                index++;
            }
            var indexFirstResidue = index;

            var currentResidue = "-";

            while (index < peptide.Length)
            {
                if (StringUtilities.IsLetterAtoZ(peptide[index]))
                {
                    currentResidue = peptide[index].ToString();

                    if (!isMsgfPlus)
                    {
                        // Look for static mods that should be applied to this residue (only applies to MSGFDB, not MS-GF+)
                        for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                        {
                            var modificationType = mPeptideMods.GetModificationTypeByIndex(modIndex);

                            ModificationDefinition modificationDefinition;
                            if (modificationType == ModificationDefinition.ResidueModificationType.StaticMod)
                            {
                                modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                                if (modificationDefinition.TargetResiduesContain(peptide[index]))
                                {
                                    // Match found; update totalModMass but do not add a static mod symbol
                                    totalModMass += modificationDefinition.ModificationMass;
                                }
                            }
                            else if (index == indexFirstResidue)
                            {
                                if (modificationType == ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod && prefix == "_")
                                {
                                    // N-terminal protein static mod
                                    modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                                    totalModMass += modificationDefinition.ModificationMass;
                                }
                                else if (modificationType == ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod)
                                {
                                    // N-terminal peptide static mod
                                    modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                                    totalModMass += modificationDefinition.ModificationMass;
                                }
                            }
                        }
                    }

                    index++;

                    if (index == indexLastResidue)
                        possibleCTerminalMod = true;
                }
                else
                {
                    // Found a mod; find the extent of the mod digits
                    match = ModMassMatcher.Match(peptide, index);

                    // Note that possibleCTerminalMod will be set to True once we hit the last residue

                    // Convert the mod mass (or masses) to one or more mod symbols

                    if (ConvertMSGFModMassesToSymbols(currentResidue, match.Groups[1].Value,
                        out var modSymbols, out var dynModSymbols, msgfPlusModInfo,
                        false, possibleCTerminalMod, out var modMassFound, out containsStaticMod))
                    {
                        peptide = ReplaceMSGFModTextWithMatchedSymbol(peptide, match.Groups[1], modSymbols, dynModSymbols, isMsgfPlus, containsStaticMod);
                        totalModMass += modMassFound;

                        if (isMsgfPlus && containsStaticMod)
                        {
                            // MS-GF+ shows mod masses for static mods
                            // Thus, we have removed the static mod mass and did not add a mod symbol
                            // Therefore, leave index unchanged
                        }
                        else
                        {
                            index += modSymbols.Length;
                        }
                    }
                    else
                    {
                        var addOn = match.Groups[1].Value.Length;
                        if (addOn == 0)
                            index++;
                        else
                            index += addOn;
                    }
                }
            }

            // If any N-terminal mods were present, we need to move them to after the first residue
            // in other words, switch from #MDHTPQSQLK to M#DHTPQSQLK
            //                          or #*MNDRQLNHR to M#*NDRQLNHR

            // Update indexFirstResidue
            indexFirstResidue = 0;
            while (indexFirstResidue < peptide.Length && !StringUtilities.IsLetterAtoZ(peptide[indexFirstResidue]))
            {
                indexFirstResidue++;
            }

            if (indexFirstResidue > 0 && indexFirstResidue < peptide.Length)
            {
                var peptideNew = peptide[indexFirstResidue] + peptide.Substring(0, indexFirstResidue);
                if (indexFirstResidue < peptide.Length - 1)
                {
                    peptideNew += peptide.Substring(indexFirstResidue + 1);
                }
                peptide = peptideNew;
            }

            return prefix + peptide + suffix;
        }

        private string ReplaceMSGFModTextWithMatchedSymbol(
            string peptide,
            Capture captureGroup,
            string modSymbols,
            string dynModSymbols,
            bool isMsgfPlus,
            bool containsStaticMod)
        {
            string peptideNew;

            if (captureGroup.Index > 0)
            {
                peptideNew = peptide.Substring(0, captureGroup.Index);
            }
            else
            {
                peptideNew = string.Empty;
            }

            if (isMsgfPlus && containsStaticMod)
            {
                // MS-GF+ shows mod masses for static mods
                // However, for consistency with other PHRP results, we do not add a symbol to the peptide for this static mod
                // Catch: If we have a peptide/terminus affected by both a static and a dynamic mod, we still want the dynamic mod.
                if (!string.IsNullOrWhiteSpace(dynModSymbols))
                {
                    peptideNew += dynModSymbols;
                }
            }
            else
            {
                peptideNew += modSymbols;
            }

            if (captureGroup.Index + captureGroup.Length < peptide.Length)
            {
                peptideNew += peptide.Substring(captureGroup.Index + captureGroup.Length);
            }

            return peptideNew;
        }

        private string ReplaceTerminus(string peptide)
        {
            if (peptide.StartsWith(N_TERMINUS_SYMBOL_MSGFPlus))
            {
                peptide = PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_MSGFPlus.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_MSGFPlus))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_MSGFPlus.Length) + "." + PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        private static readonly Regex ProteinInfoMatcher = new(PROTEIN_AND_TERM_SYMBOLS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Examines proteinList to look for a semicolon separated list of proteins and terminus symbols, for example
        /// AT1G26570.1(pre=K,post=N);AT3G29360.1(pre=K,post=N);AT3G29360.2(pre=K,post=N)
        /// </summary>
        /// <param name="proteinList">Protein list to examine</param>
        /// <param name="proteinInfo">Protein information, if it is of the form ProteinName(pre=X,post=Y)</param>
        /// <returns>The name of the first protein</returns>
        private string SplitProteinList(string proteinList, IDictionary<string, TerminusChars> proteinInfo)
        {
            proteinInfo.Clear();

            var matches = ProteinInfoMatcher.Matches(proteinList);

            if (matches.Count == 0)
            {
                // No match; likely just one protein
                return TruncateProteinName(proteinList);
            }

            foreach (Match match in matches)
            {
                var proteinName = TruncateProteinName(match.Groups[1].Value);

                if (proteinInfo.ContainsKey(proteinName))
                {
                    // Skip this protein since it's already present
                }
                else
                {
                    var terminusChars = new TerminusChars
                    {
                        NTerm = match.Groups[2].Value[0],
                        CTerm = match.Groups[3].Value[0]
                    };

                    proteinInfo.Add(proteinName, terminusChars);
                }
            }

            return proteinInfo.First().Key;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<MSGFPlusSearchResult> filteredSearchResults,
            ICollection<string> errorMessages,
            bool includeFDRandPepFDR,
            bool includeEFDR,
            bool includeIMSFields,
            bool isMsgfPlus)
        {
            // Sort filteredSearchResults by ascending SpecEValue, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.SpecEValueNum, item.ScanNum, item.ChargeNum, item.Peptide, item.Protein select item;

            var resultID = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(resultID, writer, result, errorMessages, includeFDRandPepFDR, includeEFDR, includeIMSFields, isMsgfPlus);
                resultID++;
            }
        }

        private void StoreScanGroupInfo(string scanGroupFilePath, IReadOnlyCollection<ScanGroupInfo> scanGroupDetails)
        {
            try
            {
                // Only create the ScanGroup file if one or more scan groups exist
                // Step through scanGroupDetails to check for this
                var scanGroupIDPrevious = -1;
                var createFile = false;
                foreach (var udtScanGroupInfo in scanGroupDetails)
                {
                    if (udtScanGroupInfo.ScanGroupID == scanGroupIDPrevious)
                    {
                        createFile = true;
                        break;
                    }
                    scanGroupIDPrevious = udtScanGroupInfo.ScanGroupID;
                }

                if (createFile)
                {
                    using var writer = new StreamWriter(new FileStream(scanGroupFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                    writer.WriteLine("Scan_Group_ID\tCharge\tScan");

                    foreach (var udtScanGroupInfo in scanGroupDetails)
                    {
                        writer.WriteLine(udtScanGroupInfo.ScanGroupID + "\t" + udtScanGroupInfo.Charge + "\t" + udtScanGroupInfo.Scan);
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating ScanGroupInfo file", ex);
            }
        }

        /// <summary>
        /// Stores the first hits file matches for a single scan
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Filtered search results</param>
        private void StoreTopFHTMatch(
            IList<MSGFPlusSearchResult> searchResults,
            int startIndex,
            int endIndex,
            List<MSGFPlusSearchResult> filteredSearchResults)
        {
            AssignRankByScore(searchResults, startIndex, endIndex);

            // The calling procedure should have already sorted by scan, charge, and SpecEValue; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            // Now store the first match for each charge for this scan
            // When storing, we use the protein name that occurred first in the FASTA file

            var udtCurrentResult = searchResults[startIndex];
            var currentCharge = udtCurrentResult.ChargeNum;
            var currentProteinNumber = int.MaxValue;
            var currentPeptide = GetCleanSequence(udtCurrentResult.Peptide);

            for (var index = startIndex; index <= endIndex; index++)
            {
                if (currentCharge != searchResults[index].ChargeNum)
                {
                    // New charge state
                    // Store udtCurrentResult (from the previous charge state)
                    filteredSearchResults.Add(udtCurrentResult);

                    udtCurrentResult = searchResults[index];
                    currentCharge = udtCurrentResult.ChargeNum;
                    currentProteinNumber = int.MaxValue;
                    currentPeptide = GetCleanSequence(udtCurrentResult.Peptide);
                }

                var newPeptide = GetCleanSequence(searchResults[index].Peptide);
                if (currentPeptide.Equals(newPeptide))
                {
                    var bestProtein = GetBestProteinName(udtCurrentResult.Protein, currentProteinNumber, searchResults[index].Protein);
                    if (bestProtein.Value < currentProteinNumber)
                    {
                        currentProteinNumber = bestProtein.Value;
                        if (!udtCurrentResult.Protein.Equals(bestProtein.Key))
                        {
                            udtCurrentResult.Protein = bestProtein.Key;
                        }
                    }
                }
            }

            // Store udtCurrentResult (from the previous charge state)
            filteredSearchResults.Add(udtCurrentResult);
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Filtered search results</param>
        private void StoreSynMatches(
            IList<MSGFPlusSearchResult> searchResults,
            int startIndex,
            int endIndex,
            List<MSGFPlusSearchResult> filteredSearchResults)
        {
            AssignRankByScore(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            var results = new SortedSet<string>();

            // Now store or write out the matches that pass the filters
            // By default, filter passing peptides have MSGFDB_SpecEValue <= 5E-7 Or EValue less than 0.75 or QValue less than 1% (but not 0)
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (MSGFPlusResultPassesSynFilter(searchResults[index]))
                {
                    // Check for identical results
                    // This can happen if the search used an n-terminal dynamic mod that also could apply to a specific residue and that residue is at the N-terminus
                    // For example:
                    //   R.S+229.163IGLPDVHSGYGFAIGNMAAFDMNDPEAVVSPGGVGFDINC+57.021GVR.L
                    //   R.+229.163SIGLPDVHSGYGFAIGNMAAFDMNDPEAVVSPGGVGFDINC+57.021GVR.L
                    //
                    // This software will have turned both of these results info:
                    //  R.S#IGLPDVHSGYGFAIGNMAAFDMNDPEAVVSPGGVGFDINCGVR.L
                    // We only want to include the result once in the _syn.txt file
                    // (though if the peptide maps to multiple proteins it will be listed multiple times; one line per protein)

                    var resultKey = searchResults[index].Peptide + "_" +
                                    searchResults[index].Protein + "_" +
                                    searchResults[index].MH + "_" +
                                    searchResults[index].SpecEValue;

                    if (results.Contains(resultKey))
                    {
                        continue;
                    }

                    results.Add(resultKey);
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorMessages"></param>
        /// <param name="includeFDRandPepFDR"></param>
        /// <param name="includeEFDR"></param>
        /// <param name="includeIMSFields"></param>
        /// <param name="isMsgfPlus"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ICollection<string> errorMessages,
            bool includeFDRandPepFDR,
            bool includeEFDR,
            bool includeIMSFields,
            bool isMsgfPlus)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var knownHeaderColumns = MSGFPlusSynFileReader.GetColumnHeaderNamesAndIDs(true, true);

                var data = new List<string>
                {
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.ResultID),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.Scan),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.FragMethod),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.SpecIndex),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.Charge),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.PrecursorMZ),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.DelM),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.DelMPPM),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.MH),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.Peptide),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.Protein),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.NTT),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.DeNovoScore),
                    MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.MSGFScore)
                };

                if (isMsgfPlus)
                {
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.SpecEValue));
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.RankSpecEValue));
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.EValue));
                }
                else
                {
                    data.Add(MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.SpecProb));
                    data.Add(MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.RankSpecProb));
                    data.Add(MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PValue));
                }

                if (includeFDRandPepFDR)
                {
                    if (isMsgfPlus)
                    {
                        data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.QValue));
                        data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.PepQValue));
                    }
                    else
                    {
                        data.Add(MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.FDR));
                        data.Add(MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PepFDR));
                    }
                }
                else if (includeEFDR)
                {
                    // Note that we'll write out a "1" for "PepFDR" for every result
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.EFDR));
                    data.Add(MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PepFDR));
                }

                if (isMsgfPlus)
                {
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.IsotopeError));
                }

                if (includeIMSFields)
                {
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.IMSScan));
                    data.Add(MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.IMSDriftTime));
                }

                foreach (var headerName in data)
                {
                    if (!knownHeaderColumns.ContainsKey(headerName))
                    {
                        errorMessages.Add(string.Format(
                            "Unrecognized header name for the synopsis / first hits file: {0}", headerName));
                    }
                }

                writer.WriteLine(StringUtilities.CollapseList(data));
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
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="includeFDRandPepFDR"></param>
        /// <param name="includeEFDR"></param>
        /// <param name="includeIMSFields"></param>
        /// <param name="isMsgfPlus"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            MSGFPlusSearchResult udtSearchResult,
            ICollection<string> errorMessages,
            bool includeFDRandPepFDR,
            bool includeEFDR,
            bool includeIMSFields,
            bool isMsgfPlus)
        {
            try
            {
                // Primary Columns (other columns are added in certain circumstances):
                //
                // MSGFDB
                // ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecProb    Rank_MSGFDB_SpecProb    PValue  FDR     PepFDR

                // MS-GF+
                // ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecEValue  Rank_MSGFDB_SpecEValue  EValue  QValue  PepQValue  IsotopeError

                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.Scan,
                    udtSearchResult.FragMethod,
                    udtSearchResult.SpecIndex,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.PMErrorDa,
                    udtSearchResult.PMErrorPPM,
                    udtSearchResult.MH,
                    udtSearchResult.Peptide,
                    udtSearchResult.Protein,
                    udtSearchResult.NTT,
                    udtSearchResult.DeNovoScore,
                    udtSearchResult.MSGFScore,
                    udtSearchResult.SpecEValue,
                    udtSearchResult.RankSpecEValue.ToString(),
                    udtSearchResult.EValue
                };

                if (includeFDRandPepFDR)
                {
                    data.Add(StringUtilities.TrimZeroIfNotFirstID(resultID, udtSearchResult.QValue));
                    data.Add(StringUtilities.TrimZeroIfNotFirstID(resultID, udtSearchResult.PepQValue));
                }
                else if (includeEFDR)
                {
                    data.Add(StringUtilities.TrimZeroIfNotFirstID(resultID, udtSearchResult.QValue));
                    data.Add("1");
                }

                if (isMsgfPlus)
                {
                    data.Add(udtSearchResult.IsotopeError.ToString());
                }

                if (includeIMSFields)
                {
                    data.Add(udtSearchResult.IMSScan.ToString());
                    data.Add(udtSearchResult.IMSDriftTime);
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
            return TOOL_NAME + " results processor";
        }

        private class MSGFPlusSearchResultsComparerScanChargeSpecEValuePeptide : IComparer<MSGFPlusSearchResult>
        {
            public int Compare(MSGFPlusSearchResult x, MSGFPlusSearchResult y)
            {
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

                // Charge is the same; check SpecEValue
                if (x.SpecEValueNum > y.SpecEValueNum)
                {
                    return 1;
                }

                if (x.SpecEValueNum < y.SpecEValueNum)
                {
                    return -1;
                }

                // SpecEValue is the same; check peptide
                var result = string.CompareOrdinal(x.Peptide, y.Peptide);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.CompareOrdinal(x.Protein, y.Protein);
                }
                return result;
            }
        }
    }
}
