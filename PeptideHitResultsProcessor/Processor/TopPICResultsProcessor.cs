// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
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
    /// This class reads a TopPIC results file (txt format) and creates
    /// a tab-delimited text file with the data
    /// </summary>
    public class TopPICResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: acetyl, Cysteine, Da, enums, fht, Frag, methylation, monoisotopic
        // Ignore Spelling: proteoform, proteoforms, Prsm, Pvalue, syn, toppic, txt

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="options"></param>
        public TopPICResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "September 12, 2022";
        }

        /// <summary>
        /// TopPIC tool name
        /// </summary>
        public const string TOOL_NAME = "TopPIC";

        /// <summary>
        /// TopPIC proteoforms file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_TopPIC_PROTEOFORMS_FILE = "_TopPIC_Proteoforms";

        /// <summary>
        /// TopPIC results file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_TopPIC_PRSMs_FILE = "_TopPIC_PrSMs";

        /// <summary>
        /// N-terminus symbol used by TopPIC
        /// </summary>
        public const string N_TERMINUS_SYMBOL_TopPIC = ".";

        /// <summary>
        /// C-terminus symbol used by TopPIC
        /// </summary>
        public const string C_TERMINUS_SYMBOL_TopPIC = ".";

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        /// <summary>
        /// RegEx to match mods in square brackets, examples:
        /// [11.94403]
        /// [15.98154]
        /// [-109.08458]
        /// [4.52e-003]
        /// [Acetyl]
        /// </summary>
        private const string TopPIC_MOD_MASS_OR_NAME_REGEX = @"\[(?<ModMass>[+-]*[0-9\.e-]+)\]|\[(?<NamedMod>[^\]]+)\]";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to the tab-delimited file created directly by TopPIC
        /// </summary>
        private enum TopPICResultsFileColumns
        {
            SpectrumFileName = 0,
            Prsm_ID = 1,
            Spectrum_ID = 2,
            FragMethod = 3,
            Scans = 4,
            RetentionTime = 5,
            Peaks = 6,
            Charge = 7,
            Precursor_mass = 8,             // Monoisotopic mass value of the observed precursor_mz
            Adjusted_precursor_mass = 9,    // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC; typically identical to Proteoform_mass
            Proteoform_ID = 10,
            Feature_intensity = 11,
            Feature_score = 12,
            Feature_apex_time = 13,         // Feature apex in v1.4; Feature apex time in v1.5
            Protein_hits = 14,
            Protein_accession = 15,
            Protein_description = 16,
            First_residue = 17,
            Last_residue = 18,
            Special_amino_acids = 19,
            Proteoform = 20,
            Proteoform_mass = 21,           // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC; typically identical to Adjusted_precursor_mass
            Protein_Nterminal_form = 22,
            Unexpected_modifications = 23,
            Variable_PTMs = 24,
            MIScore = 25,
            Matched_peaks = 26,
            Matched_fragment_ions = 27,
            Pvalue = 28,                    // Deprecated with 1.5
            Evalue = 29,
            Qvalue = 30,                    // Spectral FDR, or PepFDR, or Spectrum-level Q-value
            Proteoform_QValue = 31          // Proteoform-level Q-value
        }

        private int mDeltaMassWarningCount;

        private readonly SortedSet<string> mUnknownNamedMods = new();

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        private void AddDynamicAndStaticResidueMods(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const char NO_RESIDUE = '-';

            var parsingModInfo = false;
            var modMassOrName = string.Empty;

            var mostRecentResidue = NO_RESIDUE;
            var residueLocInPeptide = 0;

            var ambiguousResidue = NO_RESIDUE;
            var ambiguousResidueLocInPeptide = 0;

            var clearAmbiguousResidue = false;
            var storeAmbiguousResidue = false;

            var sequence = searchResult.PeptideSequenceWithMods;
            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var character = sequence[index];

                if (!parsingModInfo && StringUtilities.IsLetterAtoZ(character))
                {
                    mostRecentResidue = character;
                    residueLocInPeptide++;

                    if (storeAmbiguousResidue)
                    {
                        ambiguousResidue = character;
                        ambiguousResidueLocInPeptide = residueLocInPeptide;
                        storeAmbiguousResidue = false;
                    }
                    else if (clearAmbiguousResidue)
                    {
                        ambiguousResidue = NO_RESIDUE;
                        clearAmbiguousResidue = false;
                    }

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == ModificationDefinition.ResidueModificationType.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                            if (modificationDefinition.TargetResiduesContain(character))
                            {
                                // Match found; add this modification
                                var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocInPeptide);

                                searchResult.SearchResultAddModification(
                                    modificationDefinition, character, residueLocInPeptide,
                                    residueTerminusState, updateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (character == '(')
                {
                    // Start of a mod group
                    storeAmbiguousResidue = true;
                }
                else if (character == ')')
                {
                    // End of a mod group
                    clearAmbiguousResidue = true;
                }
                else if (character == '[')
                {
                    // Mod Info Start
                    modMassOrName = string.Empty;
                    parsingModInfo = true;
                }
                else if (character == ']')
                {
                    // Mod Info End

                    if (!parsingModInfo)
                        continue;

                    char residueForMod;
                    int residueLocForMod;

                    if (ambiguousResidue == NO_RESIDUE)
                    {
                        residueForMod = mostRecentResidue;
                        residueLocForMod = residueLocInPeptide;
                    }
                    else
                    {
                        // Ambiguous mod
                        // We'll associate it with the first residue of the mod group
                        residueForMod = ambiguousResidue;
                        residueLocForMod = ambiguousResidueLocInPeptide;
                    }

                    if (residueLocForMod == 0)
                    {
                        // Modification is at the peptide N-terminus
                        residueLocForMod = 1;
                    }

                    var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocForMod);

                    if (double.TryParse(modMassOrName, out var modMass))
                    {
                        var success = searchResult.SearchResultAddModification(
                            modMass, residueForMod, residueLocForMod,
                            residueTerminusState, updateModOccurrenceCounts);

                        if (!success)
                        {
                            var errorMessage = searchResult.ErrorMessage;
                            if (string.IsNullOrEmpty(errorMessage))
                            {
                                errorMessage = "SearchResultAddDynamicModification returned false for mod mass " + modMassOrName;
                            }
                            SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
                        }
                    }
                    else
                    {
                        // Named mod
                        if (LookupModificationMassByName(modMassOrName, out var modMassFromUniMod))
                        {
                            var success = searchResult.SearchResultAddModification(
                                modMassFromUniMod, residueForMod, residueLocForMod,
                                residueTerminusState, updateModOccurrenceCounts);

                            if (!success)
                            {
                                var errorMessage = searchResult.ErrorMessage;
                                if (string.IsNullOrEmpty(errorMessage))
                                {
                                    errorMessage = string.Format("SearchResultAddDynamicModification returned false for mod named {0} with mass {1} ",
                                                                 modMassOrName, modMassFromUniMod);
                                }
                                SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
                            }
                        }
                        else
                        {
                            if (string.IsNullOrEmpty(mErrorMessage))
                            {
                                var errorMessage = "Unknown mod name: " + modMassOrName;
                                SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
                            }
                            else
                            {
                                OnDebugEvent(ErrorMessage);
                            }
                        }
                    }

                    parsingModInfo = false;
                }
                else if (parsingModInfo)
                {
                    modMassOrName += character;
                }
                // ReSharper disable once RedundantIfElseBlock
                else
                {
                    // Unrecognized symbol; ignore it
                }
            }
        }

        private bool AddModificationsAndComputeMass(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool success;

            try
            {
                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts);

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

                // Populate .ProteoformModDescription
                searchResult.UpdateModDescription();

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
        }

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        /// <param name="sortOnPValue"></param>
        private void AssignRankByScore(
            IList<TopPICPrSMs> searchResults,
            int startIndex,
            int endIndex,
            bool sortOnPValue)
        {
            if (startIndex == endIndex)
            {
                // Only one result
                searchResults[startIndex].RankPValue = 1;
                return;
            }

            // Duplicate a portion of searchResults so that we can sort by PValue

            var resultsSubset = new Dictionary<int, TopPICPrSMs>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            List<KeyValuePair<int, TopPICPrSMs>> resultsByProbability;

            if (sortOnPValue)
            {
                resultsByProbability = (from item in resultsSubset orderby item.Value.PValueNum select item).ToList();
            }
            else
            {
                resultsByProbability = (from item in resultsSubset orderby item.Value.EValueNum select item).ToList();
            }

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByProbability)
            {
                var result = searchResults[entry.Key];

                var currentValue = sortOnPValue ? result.PValueNum : result.EValueNum;

                if (currentRank < 0)
                {
                    lastValue = currentValue;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(currentValue - lastValue) > double.Epsilon)
                    {
                        lastValue = currentValue;
                        currentRank++;
                    }
                }

                result.RankPValue = currentRank;
            }
        }

        private string AssureInteger(string integer, int defaultValue)
        {
            if (integer.EndsWith(".0"))
                integer = integer.Substring(0, integer.Length - 2);

            if (int.TryParse(integer, out var value))
            {
                return value.ToString();
            }

            if (double.TryParse(integer, out var doubleValue))
            {
                return doubleValue.ToString("0");
            }

            return defaultValue.ToString();
        }

        /// <summary>
        /// Copy values from sourceResult to targetResult, but only if both have the same Prsm_ID value
        /// </summary>
        /// <remarks>Does not copy protein name or description</remarks>
        /// <param name="sourceResult"></param>
        /// <param name="targetResult"></param>
        /// <returns>True if the results both have the same Prsm_ID and scores were thus cloned, otherwise false</returns>
        private bool CloneScores(TopPICPrSMs sourceResult, TopPICPrSMs targetResult)
        {
            if (sourceResult.Prsm_ID != targetResult.Prsm_ID)
            {
                OnWarningEvent(
                    "Source and target TopPIC search results do not have the same Prsm_ID: {0} vs. {1}; will not clone scores",
                    sourceResult.Prsm_ID, targetResult.Prsm_ID);

                return false;
            }

            if (string.IsNullOrWhiteSpace(targetResult.Spectrum_ID))
                targetResult.Spectrum_ID = sourceResult.Spectrum_ID;

            if (string.IsNullOrWhiteSpace(targetResult.FragMethod))
                targetResult.FragMethod = sourceResult.FragMethod;

            if (string.IsNullOrWhiteSpace(targetResult.Scans))
                targetResult.Scans = sourceResult.Scans;

            targetResult.ScanNum = sourceResult.ScanNum;

            if (string.IsNullOrWhiteSpace(targetResult.RetentionTime))
                targetResult.RetentionTime = sourceResult.RetentionTime;

            if (string.IsNullOrWhiteSpace(targetResult.Peaks))
                targetResult.Peaks = sourceResult.Peaks;

            if (string.IsNullOrWhiteSpace(targetResult.Charge))
                targetResult.Charge = sourceResult.Charge;

            targetResult.ChargeNum = sourceResult.ChargeNum;

            if (string.IsNullOrWhiteSpace(targetResult.Precursor_mass))
                targetResult.Precursor_mass = sourceResult.Precursor_mass;

            if (string.IsNullOrWhiteSpace(targetResult.PrecursorMZ))
                targetResult.PrecursorMZ = sourceResult.PrecursorMZ;

            if (string.IsNullOrWhiteSpace(targetResult.Adjusted_precursor_mass))
                targetResult.Adjusted_precursor_mass = sourceResult.Adjusted_precursor_mass;

            if (string.IsNullOrWhiteSpace(targetResult.MH))
                targetResult.MH = sourceResult.MH;

            if (string.IsNullOrWhiteSpace(targetResult.DelM))
                targetResult.DelM = sourceResult.DelM;

            if (string.IsNullOrWhiteSpace(targetResult.DelM_PPM))
                targetResult.DelM_PPM = sourceResult.DelM_PPM;

            if (string.IsNullOrWhiteSpace(targetResult.Proteoform_ID))
                targetResult.Proteoform_ID = sourceResult.Proteoform_ID;

            if (string.IsNullOrWhiteSpace(targetResult.Feature_Intensity))
                targetResult.Feature_Intensity = sourceResult.Feature_Intensity;

            if (string.IsNullOrWhiteSpace(targetResult.Feature_Score))
                targetResult.Feature_Score = sourceResult.Feature_Score;

            if (string.IsNullOrWhiteSpace(targetResult.Feature_Apex_Time))
                targetResult.Feature_Apex_Time = sourceResult.Feature_Apex_Time;

            if (string.IsNullOrWhiteSpace(targetResult.Protein_Hits))
                targetResult.Protein_Hits = sourceResult.Protein_Hits;

            // Skip targetResult.Protein
            // Skip targetResult.ProteinDescription

            if (string.IsNullOrWhiteSpace(targetResult.ResidueStart))
                targetResult.ResidueStart = sourceResult.ResidueStart;

            if (string.IsNullOrWhiteSpace(targetResult.ResidueEnd))
                targetResult.ResidueEnd = sourceResult.ResidueEnd;

            if (string.IsNullOrWhiteSpace(targetResult.Special_amino_acids))
                targetResult.Special_amino_acids = sourceResult.Special_amino_acids;

            // The prefix and suffix residues in the proteoform sequence correspond to the first protein listed in the TopPIC_PrSMs.txt file
            if (string.IsNullOrWhiteSpace(targetResult.Proteoform))
                targetResult.Proteoform = sourceResult.Proteoform;

            if (string.IsNullOrWhiteSpace(targetResult.Protein_Nterminal_Form))
                targetResult.Protein_Nterminal_Form = sourceResult.Protein_Nterminal_Form;

            if (string.IsNullOrWhiteSpace(targetResult.Proteoform_mass))
                targetResult.Proteoform_mass = sourceResult.Proteoform_mass;

            if (string.IsNullOrWhiteSpace(targetResult.Unexpected_Mod_Count))
                targetResult.Unexpected_Mod_Count = sourceResult.Unexpected_Mod_Count;

            if (string.IsNullOrWhiteSpace(targetResult.MIScore))
                targetResult.MIScore = sourceResult.MIScore;

            if (string.IsNullOrWhiteSpace(targetResult.VariablePTMs))
                targetResult.VariablePTMs = sourceResult.VariablePTMs;

            if (string.IsNullOrWhiteSpace(targetResult.Matched_peaks))
                targetResult.Matched_peaks = sourceResult.Matched_peaks;

            if (string.IsNullOrWhiteSpace(targetResult.Matched_fragment_ions))
                targetResult.Matched_fragment_ions = sourceResult.Matched_fragment_ions;

            if (string.IsNullOrWhiteSpace(targetResult.PValue))
                targetResult.PValue = sourceResult.PValue;

            targetResult.PValueNum = sourceResult.PValueNum;

            targetResult.RankPValue = sourceResult.RankPValue;

            if (string.IsNullOrWhiteSpace(targetResult.Evalue))
                targetResult.Evalue = sourceResult.Evalue;

            targetResult.EValueNum = sourceResult.EValueNum;

            if (string.IsNullOrWhiteSpace(targetResult.Qvalue))
                targetResult.Qvalue = sourceResult.Qvalue;

            if (string.IsNullOrWhiteSpace(targetResult.Proteoform_QValue))
                targetResult.Proteoform_QValue = sourceResult.Proteoform_QValue;

            return true;
        }

        /// <summary>
        /// Compute the peptide mass
        /// </summary>
        /// <param name="peptide">Sequence with mods, which can be either numeric or named ([15.98154] or [Acetyl])</param>
        /// <param name="totalModMass"></param>
        /// <returns>Mass of the peptide plus totalModMass</returns>
        private double ComputePeptideMass(string peptide, double totalModMass)
        {
            var cleanSequence = GetCleanSequence(peptide);

            return ComputePeptideMassForCleanSequence(cleanSequence, totalModMass);
        }

        private static readonly Regex ModMatcher = new(TopPIC_MOD_MASS_OR_NAME_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form [23.5432] or [Acetyl]</param>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            foreach (Match match in ModMatcher.Matches(peptide))
            {
                if (double.TryParse(match.Groups["ModMass"].Value, out var modMassFound))
                {
                    totalModMass += modMassFound;
                }
                else
                {
                    if (LookupModificationMassByName(match.Groups["NamedMod"].Value, out var modMass))
                        totalModMass += modMass;
                }
            }

            return totalModMass;
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
                "_toppic_syn",
                "_toppic_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// This routine creates a first hits file and/or a synopsis file from the output from TopPIC
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath">TopPIC results file (Dataset_TopPIC_PrSMs.txt)</param>
        /// <param name="fhtFilePath">First-hits file path (Dataset_toppic_fht.txt)</param>
        /// <param name="synFilePath">Synopsis file path (Dataset_toppic_syn.txt)</param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateFHTorSYNResultsFile(
            string inputFilePath,
            string fhtFilePath,
            string synFilePath)
        {
            var columnMapping = new Dictionary<TopPICResultsFileColumns, int>();

            try
            {
                var errorMessages = new List<string>();

                // Initialize a SortedSet that tracks each combo of scan and charge
                var scanChargeFirstHit = new SortedSet<string>();

                // Open the input file and parse it
                // Initialize the stream reader and the stream writer
                using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var firstHitsFileWriter = string.IsNullOrWhiteSpace(fhtFilePath)
                    ? null
                    : new StreamWriter(new FileStream(fhtFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                var synopsisFileWriter = string.IsNullOrWhiteSpace(synFilePath)
                    ? null
                    : new StreamWriter(new FileStream(synFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                var headerParsed = false;
                var dataHasPValues = false;

                mDeltaMassWarningCount = 0;

                // Initialize array that will hold all of the records in the TopPIC result file
                var searchResultsUnfiltered = new List<TopPICPrSMs>();

                // Initialize the array that will hold all of the records that will ultimately be written out to disk
                var filteredSearchResults = new List<TopPICPrSMs>();

                var previousSearchResult = new TopPICPrSMs();

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
                        var validHeader = ParseTopPICResultsFileHeaderLine(lineIn, columnMapping);
                        if (!validHeader)
                        {
                            // Error parsing header
                            SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                            return false;
                        }

                        headerParsed = true;

                        dataHasPValues = columnMapping[TopPICResultsFileColumns.Pvalue] >= 0;

                        // Write the header line to each file
                        WriteSynFHTFileHeader(firstHitsFileWriter, dataHasPValues, errorMessages);

                        WriteSynFHTFileHeader(synopsisFileWriter, dataHasPValues, errorMessages);

                        continue;
                    }

                    var validSearchResult = ParseTopPICResultsFileEntry(lineIn, out var searchResult, out var isAdditionalProtein, errorMessages, columnMapping);

                    if (validSearchResult)
                    {
                        if (isAdditionalProtein)
                        {
                            // If creating the synopsis file, copy scores from the previous search result to searchResult
                            // (both should have the same Prsm_ID value; CloneScores will return false if they do not)

                            if (synopsisFileWriter != null && CloneScores(previousSearchResult, searchResult))
                            {
                                searchResultsUnfiltered.Add(searchResult);
                            }
                        }
                        else
                        {
                            var scanChargeKey = searchResult.Scans + "_" + searchResult.Charge;

                            if (!scanChargeFirstHit.Contains(scanChargeKey))
                            {
                                scanChargeFirstHit.Add(scanChargeKey);
                                searchResult.StoreInFirstHitsFile = true;
                            }

                            searchResultsUnfiltered.Add(searchResult);
                            previousSearchResult = searchResult;
                        }
                    }

                    // Update the progress
                    UpdateSynopsisFileCreationProgress(reader);
                }

                // Sort the SearchResults by scan, charge, and ascending PValue
                searchResultsUnfiltered.Sort(new TopPICSearchResultsComparerScanChargePValuePeptide());

                // Now filter the data

                // Initialize variables
                var startIndex = 0;

                while (startIndex < searchResultsUnfiltered.Count)
                {
                    var endIndex = startIndex;
                    while (endIndex + 1 < searchResultsUnfiltered.Count &&
                           searchResultsUnfiltered[endIndex + 1].ScanNum == searchResultsUnfiltered[startIndex].ScanNum)
                    {
                        endIndex++;
                    }

                    // Store the results for this scan
                    StoreSynMatches(searchResultsUnfiltered, startIndex, endIndex, filteredSearchResults, dataHasPValues);

                    startIndex = endIndex + 1;
                }

                // Sort the data in filteredSearchResults then write out to disk
                SortAndWriteFilteredSearchResults(firstHitsFileWriter, synopsisFileWriter, filteredSearchResults, errorMessages, dataHasPValues);

                // Inform the user if any errors occurred
                if (errorMessages.Count > 0)
                {
                    SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateFHTorSYNResultsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the TopPIC parameter file
        /// </summary>
        /// <remarks>The DMS-based parameter file for TopPIC uses the same formatting as MS-GF+</remarks>
        /// <param name="topPICParamFilePath"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromParamFile(string topPICParamFilePath)
        {
            var modFileProcessor = new MSGFPlusParamFileModExtractor(TOOL_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                topPICParamFilePath,
                MSGFPlusParamFileModExtractor.ModSpecFormats.TopPIC,
                out var modList);

            if (!success || mErrorCode != PHRPErrorCode.NoError)
            {
                if (mErrorCode == PHRPErrorCode.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the TopPIC parameter file");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modList, mPeptideMods);

            return true;
        }

        private new string GetCleanSequence(string sequenceWithMods)
        {
            return GetCleanSequence(sequenceWithMods, out _, out _, out _);
        }

        private string GetCleanSequence(string sequenceWithMods, out string prefix, out string suffix, out string primarySequenceWithMods)
        {
            string primarySequence;

            if (PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out primarySequenceWithMods, out prefix, out suffix))
            {
                // Remove all mods
                primarySequence = ModMatcher.Replace(primarySequenceWithMods, string.Empty);
            }
            else
            {
                // Sequence does not have prefix or suffix letters; use sequenceWithMods
                primarySequenceWithMods = sequenceWithMods;
                primarySequence = ModMatcher.Replace(sequenceWithMods, string.Empty);
            }

            // The primary sequence may still contain parentheses; remove them
            return base.GetCleanSequence(primarySequence);
        }

        private bool LookupModificationMassByName(string modName, out double modMass)
        {
            if (mPeptideMods.LookupModificationMassByName(modName, out modMass))
                return true;

            if (!mUnknownNamedMods.Contains(modName))
            {
                mUnknownNamedMods.Add(modName);
                OnErrorEvent("Unrecognized named mod: " + modName);
            }

            return false;
        }

        private bool ParseTopPICSynopsisFile(string inputFilePath, string outputDirectoryPath, ref List<PepToProteinMapping> pepToProteinMapping, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Note that TopPIC synopsis files are normally sorted on PValue, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            // we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
            // (see peptidesFoundForPValueLevel below)

            // Although this was a possibility with InSpecT, it likely never occurs for TopPIC
            // But, we'll keep the check in place just in case

            var columnMapping = new Dictionary<TopPICSynFileColumns, int>();

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
                var searchResult = new TopPICResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize a SortedSet that will be used to avoid double-counting the same PSM in the same scan
                // This is required since a PSM with multiple proteins will be listed on multiple lines in the synopsis file
                // Values are PeptideSequenceWithMods_Scan_Charge

                var peptidesFoundForPValueLevel = new SortedSet<string>();
                var previousPValue = string.Empty;

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
                            var validHeader = ParseTopPICSynFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }
                            headerParsed = true;
                            continue;
                        }

                        var validSearchResult = ParseTopPICSynFileEntry(
                            lineIn, searchResult, errorMessages,
                            resultsProcessed, columnMapping,
                            out var currentPeptideWithMods);

                        resultsProcessed++;
                        if (!validSearchResult)
                            continue;

                        var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                        bool firstMatchForGroup;
                        if (searchResult.PValue == previousPValue)
                        {
                            // New result has the same PValue as the previous result
                            // See if peptidesFoundForPValueLevel contains the peptide, scan and charge

                            if (peptidesFoundForPValueLevel.Contains(key))
                            {
                                firstMatchForGroup = false;
                            }
                            else
                            {
                                peptidesFoundForPValueLevel.Add(key);
                                firstMatchForGroup = true;
                            }
                        }
                        else
                        {
                            // New PValue
                            // Reset peptidesFoundForPValueLevel
                            peptidesFoundForPValueLevel.Clear();

                            // Update previousPValue
                            previousPValue = searchResult.PValue;

                            // Append a new entry to peptidesFoundForPValueLevel
                            peptidesFoundForPValueLevel.Add(key);
                            firstMatchForGroup = true;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
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
                                OnWarningEvent("no match for '" + currentPeptideWithMods + "' in pepToProteinMapping");
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

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in ParseTopPICSynopsisFile", ex);
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
                SetErrorMessage("Error creating the output file in ParseTopPICSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parses an entry from the TopPIC results file (Dataset_TopPIC_PrSMs.txt)
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="isAdditionalProtein">Output: true if this is an additional protein for the current PSM</param>
        /// <param name="errorMessages"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICResultsFileEntry(
            string lineIn,
            out TopPICPrSMs searchResult,
            out bool isAdditionalProtein,
            ICollection<string> errorMessages,
            IDictionary<TopPICResultsFileColumns, int> columnMapping)
        {
            searchResult = new TopPICPrSMs();
            isAdditionalProtein = false;

            string[] splitLine = null;

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                // The file should have over 20 columns, but we'll only require 15
                if (splitLine.Length < 15)
                {
                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.SpectrumFileName], out searchResult.SpectrumFileName);
                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Prsm_ID], out searchResult.Prsm_ID))
                {
                    ReportError("Prsm_ID column is missing or invalid", true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Spectrum_ID], out searchResult.Spectrum_ID);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.FragMethod], out searchResult.FragMethod);

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Scans], out searchResult.Scans))
                {
                    ReportError("Scan(s) column is missing or invalid", true);
                }

                if (!string.IsNullOrWhiteSpace(searchResult.Scans) && !int.TryParse(searchResult.Scans, out searchResult.ScanNum))
                {
                    // .Scans may have a list of scan numbers; extract the first scan number from .scans
                    var scanNumberDigits = string.Empty;
                    foreach (var character in searchResult.Scans)
                    {
                        if (char.IsDigit(character))
                        {
                            scanNumberDigits += character;
                        }
                    }

                    if (!int.TryParse(scanNumberDigits, out searchResult.ScanNum))
                    {
                        OnWarningEvent("Error parsing out the scan number from the scan list; could not find an integer: " +
                                      searchResult.Scans);
                        searchResult.ScanNum = 0;
                    }
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.RetentionTime], out searchResult.RetentionTime);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Peaks], out searchResult.Peaks);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Charge], out searchResult.Charge);
                searchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(searchResult.Charge, 0));

                // Monoisotopic mass value of the observed precursor_mz
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Precursor_mass], out searchResult.Precursor_mass);

                var precursorMZ = 0.0;

                // precursorMonoMass is Observed m/z, converted to monoisotopic mass
                if (double.TryParse(searchResult.Precursor_mass, out var precursorMonoMass))
                {
                    if (searchResult.ChargeNum > 0)
                    {
                        precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, searchResult.ChargeNum);
                        searchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                    }
                }

                // peptideMonoMassTopPIC is theoretical peptide monoisotopic mass, including mods, as computed by TopPIC
                double peptideMonoMassTopPIC;

                if (columnMapping[TopPICResultsFileColumns.Adjusted_precursor_mass] >= 0)
                {
                    // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC; typically identical to Proteoform_Mass
                    DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Adjusted_precursor_mass], out searchResult.Adjusted_precursor_mass);

                    double.TryParse(searchResult.Adjusted_precursor_mass, out peptideMonoMassTopPIC);
                }
                else
                {
                    peptideMonoMassTopPIC = 0;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform_ID], out searchResult.Proteoform_ID);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Feature_intensity], out searchResult.Feature_Intensity);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Feature_score], out searchResult.Feature_Score);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Feature_apex_time], out searchResult.Feature_Apex_Time);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Protein_hits], out searchResult.Protein_Hits);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Protein_accession], out searchResult.Protein);
                searchResult.Protein = TruncateProteinName(searchResult.Protein);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Protein_description], out searchResult.ProteinDescription);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.First_residue], out searchResult.ResidueStart);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Last_residue], out searchResult.ResidueEnd);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Special_amino_acids], out searchResult.Special_amino_acids);

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform], out searchResult.Proteoform))
                {
                    // Starting with TopPIC 1.5, proteoforms that map to multiple proteins will be listed multiple times in the _TopPIC_PrSMs.txt file
                    // The first protein will be in a result line with all of the scores
                    // Subsequent proteins will be in a line with Data File Name, Prism ID, Protein accession, and Protein description, while all of the other columns are blank

                    var emptyColumns = 0;

                    if (string.IsNullOrWhiteSpace(searchResult.Spectrum_ID))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.FragMethod))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Scans))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.RetentionTime))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Peaks))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Charge))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Precursor_mass))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Proteoform_ID))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Feature_Intensity))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Feature_Score))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Feature_Apex_Time))
                        emptyColumns++;
                    if (string.IsNullOrWhiteSpace(searchResult.Protein_Hits))
                        emptyColumns++;

                    if (emptyColumns > 5)
                        isAdditionalProtein = true;
                    else
                        ReportError("Proteoform column is missing or invalid", true);
                }

                if (!isAdditionalProtein)
                {
                    // Add the standard terminus symbols to the peptide sequence
                    searchResult.Proteoform = ReplaceTerminus(searchResult.Proteoform);

                    // Parse the sequence to determine the total mod mass
                    // Note that we do not remove any of the mod symbols since TopPIC since mods can ambiguously apply to residues
                    var totalModMass = ComputeTotalModMass(searchResult.Proteoform);

                    // Compute theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                    var peptideMonoMassPHRP = ComputePeptideMass(searchResult.Proteoform, totalModMass);

                    if (Math.Abs(peptideMonoMassTopPIC) < double.Epsilon)
                    {
                        peptideMonoMassTopPIC = peptideMonoMassPHRP;
                    }

                    // Warn the user if the monoisotopic mass values differ by more than 0.1 Da
                    ValidateMatchingMonoisotopicMass(TOOL_NAME, searchResult.Proteoform, peptideMonoMassPHRP, peptideMonoMassTopPIC,
                        ref mDeltaMassWarningCount);

                    if (peptideMonoMassTopPIC > 0)
                    {
                        // Compute DelM and DelM_PPM
                        var delM = precursorMonoMass - peptideMonoMassTopPIC;
                        searchResult.DelM = StringUtilities.MassErrorToString(delM);

                        if (precursorMZ > 0)
                        {
                            searchResult.DelM_PPM =
                                PRISM.StringUtilities.DblToString(PeptideMassCalculator.MassToPPM(delM, precursorMZ), 5, 0.00005);
                        }
                        else
                        {
                            searchResult.DelM_PPM = PRISM.StringUtilities.DblToString(PeptideMassCalculator.MassToPPM(delM, 1000), 5, 0.00005);
                        }
                    }

                    // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    searchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoMassPHRP, 0), 6);
                }

                // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC; typically identical to Adjusted_precursor_mass
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform_mass], out searchResult.Proteoform_mass);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Protein_Nterminal_form], out searchResult.Protein_Nterminal_Form);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Unexpected_modifications], out searchResult.Unexpected_Mod_Count);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Variable_PTMs], out searchResult.VariablePTMs);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.MIScore], out searchResult.MIScore);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Matched_peaks], out searchResult.Matched_peaks);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Matched_fragment_ions], out searchResult.Matched_fragment_ions);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Pvalue], out searchResult.PValue);
                if (!double.TryParse(searchResult.PValue, out searchResult.PValueNum))
                    searchResult.PValueNum = 0;

                if (!isAdditionalProtein)
                {
                    // Assure that the following are truly integers (Matched_peaks and Matched_fragment_ions are often of the form 8.0)
                    searchResult.Unexpected_Mod_Count = AssureInteger(searchResult.Unexpected_Mod_Count, 0);   // Unexpected_Mod_Count
                    searchResult.Peaks = AssureInteger(searchResult.Peaks, 0);                                 // Peak_count
                    searchResult.Matched_peaks = AssureInteger(searchResult.Matched_peaks, 0);                 // Matched_Peak_Count
                    searchResult.Matched_fragment_ions = AssureInteger(searchResult.Matched_fragment_ions, 0); // Matched_Fragment_Ion_Count
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Evalue], out searchResult.Evalue);
                if (!double.TryParse(searchResult.Evalue, out searchResult.EValueNum))
                    searchResult.EValueNum = 0;

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Qvalue], out searchResult.Qvalue);

                if (string.Equals(searchResult.Qvalue, "infinity", StringComparison.OrdinalIgnoreCase))
                {
                    searchResult.Qvalue = "10";
                }
                else if (!string.IsNullOrEmpty(searchResult.Qvalue) && !double.TryParse(searchResult.Qvalue, out _))
                {
                    searchResult.Qvalue = string.Empty;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform_QValue], out searchResult.Proteoform_QValue);

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the TopPIC results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing TopPIC Results in ParseTopPICResultsFileEntry for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing TopPIC Results in ParseTopPICResultsFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a TopPIC results file (Dataset_TopPIC_PrSMs.txt), populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICResultsFileHeaderLine(string lineIn, IDictionary<TopPICResultsFileColumns, int> columnMapping)
        {
            // Headers prior to November 2018:
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Protein name    First residue    Last residue                                                                                                                                Proteoform                                                  #unexpected modifications                                 #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR    #Variable PTMs

            // Headers for TopPIC 1.2 and 1.3
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity                                                           Protein accession    Protein description    First residue    Last residue                           Proteoform                                                  #unexpected modifications    MIScore    #variable PTMs    #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR

            // Headers for TopPIC 1.4
            // The P-Value column has been removed and the Q-Value and "Proteoform FDR" columns have been renamed
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Feature score                                          Protein accession    Protein description    First residue    Last residue                           Proteoform                                                  #unexpected modifications    MIScore    #variable PTMs    #matched peaks    #matched fragment ions    E-value    Spectrum-level Q-value    Proteoform-level Q-value

            // Headers for TopPIC 1.4.13
            // Column Feature apex was added
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Feature score    Feature apex                          Protein accession    Protein description    First residue    Last residue                           Proteoform                                                  #unexpected modifications    MIScore    #variable PTMs    #matched peaks    #matched fragment ions    E-value    Spectrum-level Q-value    Proteoform-level Q-value

            // Headers for TopPIC 1.5.4
            // "Feature apex" was renamed to "Feature apex time", columns "MIScore" and "#variable PTMs" swapped places, and four columns were added
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Feature score    Feature apex time    #Protein hits    Protein accession    Protein description    First residue    Last residue    Special amino acids    Proteoform    Proteoform mass    Protein N-terminal form    #unexpected modifications    #variable PTMs    MIScore    #matched peaks    #matched fragment ions    E-value    Spectrum-level Q-value    Proteoform-level Q-value

            var columnNames = new SortedDictionary<string, TopPICResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Data file name", TopPICResultsFileColumns.SpectrumFileName},
                {"Prsm ID", TopPICResultsFileColumns.Prsm_ID},
                {"Spectrum ID", TopPICResultsFileColumns.Spectrum_ID},
                {"Fragmentation", TopPICResultsFileColumns.FragMethod},
                {"Scan(s)", TopPICResultsFileColumns.Scans},
                {"Retention time", TopPICResultsFileColumns.RetentionTime},
                {"#peaks", TopPICResultsFileColumns.Peaks},
                {"Charge", TopPICResultsFileColumns.Charge},
                {"Precursor mass", TopPICResultsFileColumns.Precursor_mass},
                {"Adjusted precursor mass", TopPICResultsFileColumns.Adjusted_precursor_mass},
                {"Proteoform ID", TopPICResultsFileColumns.Proteoform_ID},
                {"Feature intensity", TopPICResultsFileColumns.Feature_intensity},
                {"Feature score", TopPICResultsFileColumns.Feature_score},
                {"Feature apex", TopPICResultsFileColumns.Feature_apex_time},
                {"Feature apex time", TopPICResultsFileColumns.Feature_apex_time},
                {"#Protein hits", TopPICResultsFileColumns.Protein_hits},
                {"Protein name", TopPICResultsFileColumns.Protein_accession},
                {"Protein accession", TopPICResultsFileColumns.Protein_accession},
                {"Protein description", TopPICResultsFileColumns.Protein_description},
                {"First residue", TopPICResultsFileColumns.First_residue},
                {"Last residue", TopPICResultsFileColumns.Last_residue},
                {"Special amino acids", TopPICResultsFileColumns.Special_amino_acids},
                {"Proteoform", TopPICResultsFileColumns.Proteoform},
                {"Proteoform mass", TopPICResultsFileColumns.Proteoform_mass},
                {"Protein N-terminal form", TopPICResultsFileColumns.Protein_Nterminal_form},
                {"#unexpected modifications", TopPICResultsFileColumns.Unexpected_modifications},
                {"#variable PTMs", TopPICResultsFileColumns.Variable_PTMs},
                {"MIScore", TopPICResultsFileColumns.MIScore},
                {"#matched peaks", TopPICResultsFileColumns.Matched_peaks},
                {"#matched fragment ions", TopPICResultsFileColumns.Matched_fragment_ions},
                {"P-value", TopPICResultsFileColumns.Pvalue},
                {"E-value", TopPICResultsFileColumns.Evalue},
                {"Q-value (spectral FDR)", TopPICResultsFileColumns.Qvalue},
                {"Spectrum-level Q-value", TopPICResultsFileColumns.Qvalue},
                {"Proteoform FDR", TopPICResultsFileColumns.Proteoform_QValue},
                {"Proteoform-level Q-value", TopPICResultsFileColumns.Proteoform_QValue}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (TopPICResultsFileColumns resultColumn in Enum.GetValues(typeof(TopPICResultsFileColumns)))
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
                        Console.WriteLine(
                            "Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseTopPICResultsFileHeaderLine");
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the TopPIC results file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse the header line of a TopPIC _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICSynFileHeaderLine(string lineIn, IDictionary<TopPICSynFileColumns, int> columnMapping)
        {
            var columnNames = TopPICSynFileReader.GetColumnHeaderNamesAndIDs(true);

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (TopPICSynFileColumns resultColumn in Enum.GetValues(typeof(TopPICSynFileColumns)))
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
                SetErrorMessage("Error parsing the header line in the TopPIC synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from the TopPIC Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequenceWithMods"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICSynFileEntry(
            string lineIn,
            TopPICResults searchResult,
            ICollection<string> errorMessages,
            int resultsProcessed,
            IDictionary<TopPICSynFileColumns, int> columnMapping,
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

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.ResultID], out string value))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading ResultID value from TopPIC results line '{0}'", resultsProcessed + 1));
                    }

                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Scan], out string scan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading peptide sequence from TopPIC results, line {0}", resultsProcessed + 1));
                    }

                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.DelM], out string msAlignComputedDelM);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.DelMPPM], out string msAlignComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.TopPICComputedDelM = msAlignComputedDelM;
                searchResult.TopPICComputedDelMPPM = msAlignComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.TopPICComputedDelM;

                // Note: .PeptideDeltaMass is stored in the TopPIC results file as "Observed_Mass - Theoretical_Mass"
                // However, in MTS .PeptideDeltaMass is "Theoretical - Observed"
                // Therefore, we will negate .PeptideDeltaMass
                try
                {
                    searchResult.PeptideDeltaMass = (-double.Parse(searchResult.PeptideDeltaMass)).ToString(CultureInfo.InvariantCulture);
                }
                catch (Exception)
                {
                    // Error; Leave .PeptideDeltaMass unchanged
                }

                // Since peptideSequenceWithMods may contain named mods in square brackets,
                // do not call searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods)
                // and instead emulate what SetPeptideSequenceWithMods does

                var primarySequenceWithoutMods = GetCleanSequence(peptideSequenceWithMods, out var prefix, out var suffix, out var primarySequenceWithMods);

                searchResult.PeptidePreResidues = prefix;
                searchResult.PeptidePostResidues = suffix;
                searchResult.PeptideCleanSequence = primarySequenceWithoutMods;
                searchResult.PeptideSequenceWithMods = primarySequenceWithMods;

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the cleavage state and terminus state
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since InSpecT only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Prsm_ID], out string prsmId);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Spectrum_ID], out string spectrumId);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.PrecursorMZ], out string precursorMz);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.MH], out string parentIonMH);

                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Unexpected_Mod_Count], out string unexpectedModCount);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Peak_Count], out string peakCount);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Matched_Peak_Count], out string matchedPeakCount);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Matched_Fragment_Ion_Count], out string matchedFragmentIonCount);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.PValue], out string pValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Rank_PValue], out string rankPValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.EValue], out string eValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.QValue], out string qValue);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.FragMethod], out string fragMethod);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Proteoform_QValue], out string proteoformFDR);
                DataUtilities.GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.VariablePTMs], out string variablePTMs);

                searchResult.Prsm_ID = prsmId;
                searchResult.Spectrum_ID = spectrumId;
                searchResult.Precursor_mz = precursorMz;
                searchResult.ParentIonMH = parentIonMH;
                searchResult.Unexpected_Mod_Count = unexpectedModCount;
                searchResult.Peak_Count = peakCount;
                searchResult.Matched_Peak_Count = matchedPeakCount;
                searchResult.Matched_Fragment_Ion_Count = matchedFragmentIonCount;
                searchResult.PValue = pValue;
                searchResult.Rank_PValue = rankPValue;
                searchResult.EValue = eValue;
                searchResult.QValue = qValue;
                searchResult.FragMethod = fragMethod;
                searchResult.ProteoformFDR = proteoformFDR;
                searchResult.VariablePTMs = variablePTMs;

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
                            "Error parsing TopPIC results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing TopPIC Results in ParseTopPICSynFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">TopPIC results file (Dataset_TopPIC_PrSMs.txt)</param>
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

                ResetProgress("Parsing " + Path.GetFileName(inputFilePath));

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath, Options.AlternateBasePath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);

                    var pepToProteinMapping = new List<PepToProteinMapping>();

                    var topPICParameterFilePath = ResolveFilePath(inputFile.DirectoryName, Options.SearchToolParameterFilePath);

                    // Load the TopPIC Parameter File so that we can determine whether Cysteine residues are statically modified
                    var modInfoExtracted = ExtractModInfoFromParamFile(topPICParameterFilePath);

                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_TopPIC_PrSMs" or "_TopPIC_Proteoforms" with "_toppic"
                    foreach (var suffix in new List<string> { FILENAME_SUFFIX_TopPIC_PRSMs_FILE, FILENAME_SUFFIX_TopPIC_PROTEOFORMS_FILE })
                    {
                        if (!baseName.EndsWith(suffix, StringComparison.OrdinalIgnoreCase))
                            continue;
                        baseName = baseName.Substring(0, baseName.Length - suffix.Length) + "_toppic";
                        break;
                    }

                    if (!Options.CreateFirstHitsFile && !Options.CreateSynopsisFile)
                    {
                        OnWarningEvent("Both 'CreateFirstHitsFile' and 'CreateSynopsisFile' are false; aborting since nothing to do");
                        return true;
                    }

                    string fhtOutputFilePath;
                    string synOutputFilePath;

                    if (Options.CreateFirstHitsFile)
                    {
                        fhtOutputFilePath = Path.Combine(outputDirectoryPath, baseName + FIRST_HITS_FILE_SUFFIX);
                    }
                    else
                    {
                        fhtOutputFilePath = string.Empty;
                    }

                    if (Options.CreateSynopsisFile)
                    {
                        // The synopsis file name will be of the form BasePath_msalign_syn.txt
                        synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);
                    }
                    else
                    {
                        synOutputFilePath = string.Empty;
                    }

                    // Create the first hits output file
                    ResetProgress("Creating the FHT and SYN files", true);

                    success = CreateFHTorSYNResultsFile(inputFilePath, fhtOutputFilePath, synOutputFilePath);

                    if (!success)
                    {
                        return false;
                    }

                    if (Options.CreateSynopsisFile)
                    {
                        // Create the other PHRP-specific files
                        ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                        // Now parse the _syn.txt file that we just created to next create the other PHRP files
                        success = ParseTopPICSynopsisFile(synOutputFilePath, outputDirectoryPath, ref pepToProteinMapping, false);

                        if (!success)
                        {
                            return false;
                        }

                        // Remove all items from pepToProteinMapping to reduce memory overhead
                        pepToProteinMapping.Clear();
                        pepToProteinMapping.TrimExcess();

                        if (Options.CreateProteinModsFile)
                        {
                            // Use a higher match error threshold because some peptides reported by TopPIC don't perfectly match the FASTA file
                            const int MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD = 50;

                            success = CreateProteinModsFileWork(
                                baseName, inputFile,
                                synOutputFilePath, outputDirectoryPath,
                                PeptideHitResultTypes.TopPIC,
                                MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD,
                                0,
                                fhtOutputFilePath);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in TopPICResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in TopPICResultsProcessor.ProcessFile (1):" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        private string ReplaceTerminus(string peptide)
        {
            if (peptide.StartsWith(N_TERMINUS_SYMBOL_TopPIC))
            {
                peptide = PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_TopPIC.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_TopPIC))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_TopPIC.Length) + "." + PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter firstHitsFileWriter,
            TextWriter synopsisFileWriter,
            IEnumerable<TopPICPrSMs> filteredSearchResults,
            ICollection<string> errorMessages,
            bool dataHasPValues)
        {
            IOrderedEnumerable<TopPICPrSMs> query;

            if (dataHasPValues)
            {
                // Sort filteredSearchResults by ascending PValue, ascending scan, ascending charge, ascending peptide, and ascending protein
                query = from item in filteredSearchResults orderby item.PValueNum, item.ScanNum, item.ChargeNum, item.Proteoform, item.Protein select item;
            }
            else
            {
                // Sort filteredSearchResults by ascending EValue, ascending scan, ascending charge, ascending peptide, and ascending protein
                query = from item in filteredSearchResults orderby item.EValueNum, item.ScanNum, item.ChargeNum, item.Proteoform, item.Protein select item;
            }

            var index = 1;
            foreach (var result in query)
            {
                if (firstHitsFileWriter != null && result.StoreInFirstHitsFile)
                {
                    WriteSearchResultToFile(index, firstHitsFileWriter, result, dataHasPValues, errorMessages);
                }

                if (synopsisFileWriter != null)
                {
                    WriteSearchResultToFile(index, synopsisFileWriter, result, dataHasPValues, errorMessages);
                }

                index++;
            }
        }

        private void StoreSynMatches(
            IList<TopPICPrSMs> searchResults,
            int startIndex,
            int endIndex,
            ICollection<TopPICPrSMs> filteredSearchResults,
            bool dataHasPValues)
        {
            AssignRankByScore(searchResults, startIndex, endIndex, dataHasPValues);

            // The calling procedure already sorted by scan, charge, and PValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (dataHasPValues && searchResults[index].PValueNum <= Options.MSAlignAndTopPICSynopsisFilePValueThreshold ||
                    !dataHasPValues && searchResults[index].EValueNum <= Options.MSAlignAndTopPICSynopsisFilePValueThreshold)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="dataHasPValues"></param>
        /// <param name="errorMessages"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            bool dataHasPValues,
            ICollection<string> errorMessages)
        {
            if (writer == null)
                return;

            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = TopPICSynFileReader.GetColumnHeaderNamesAndIDs(false);

                if (!dataHasPValues)
                {
                    headerColumns.Remove(TopPICSynFileReader.GetColumnNameByID(TopPICSynFileColumns.PValue));
                    headerColumns.Remove(TopPICSynFileReader.GetColumnNameByID(TopPICSynFileColumns.Rank_PValue));
                }

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
        /// Writes an entry to the synopsis file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="searchResult"></param>
        /// <param name="dataHasPValues"></param>
        /// <param name="errorMessages"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            TopPICPrSMs searchResult,
            bool dataHasPValues,
            ICollection<string> errorMessages)
        {
            try
            {
                // Primary Columns
                //
                // TopPIC
                // ResultID  Scan  Prsm_ID  Spectrum_ID  FragMethod  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Proteoform_ID  Feature_Intensity  Feature_Score  Feature_Apex_Time  Protein_Count  Protein  Protein_N-terminal_Form  ResidueStart  ResidueEnd  Unexpected_Mod_Count  Peak_Count  Matched_Peak_Count  Matched_Fragment_Ion_Count  PValue  Rank_PValue  EValue  QValue  ProteoformFDR  VariablePTMs

                var data = new List<string>
                {
                    resultID.ToString(),
                    searchResult.ScanNum.ToString(),
                    searchResult.Prsm_ID,
                    searchResult.Spectrum_ID,
                    searchResult.FragMethod,
                    searchResult.Charge,
                    searchResult.PrecursorMZ,
                    searchResult.DelM,
                    searchResult.DelM_PPM,
                    searchResult.MH,                        // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
                    searchResult.Proteoform,                // aka Peptide
                    searchResult.Proteoform_ID,
                    searchResult.Feature_Intensity,
                    searchResult.Feature_Score,
                    searchResult.Feature_Apex_Time,
                    searchResult.Protein_Hits,
                    searchResult.Protein,
                    searchResult.Protein_Nterminal_Form,
                    // Skip Proteoform_mass since nearly the same as MH
                    searchResult.ResidueStart,
                    searchResult.ResidueEnd,
                    searchResult.Unexpected_Mod_Count,       // Unexpected_Mod_Count
                    searchResult.Peaks,                      // Peak_count
                    searchResult.Matched_peaks,              // Matched_Peak_Count
                    searchResult.Matched_fragment_ions       // Matched_Fragment_Ion_Count
                };

                if (dataHasPValues)
                {
                    data.Add(searchResult.PValue);
                    data.Add(searchResult.RankPValue.ToString());
                }

                data.Add(searchResult.Evalue);
                data.Add(searchResult.Qvalue);
                data.Add(searchResult.Proteoform_QValue);
                data.Add(searchResult.VariablePTMs);

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

        private class TopPICSearchResultsComparerScanChargePValuePeptide : IComparer<TopPICPrSMs>
        {
            public int Compare(TopPICPrSMs x, TopPICPrSMs y)
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

                // Charge is the same; check PValue
                var result = string.CompareOrdinal(x.PValue, y.PValue);
                if (result == 0)
                {
                    // PValue is the same; check peptide
                    result = string.CompareOrdinal(x.Proteoform, y.Proteoform);
                    if (result == 0)
                    {
                        // Peptide is the same, check Protein
                        result = string.CompareOrdinal(x.Protein, y.Protein);
                    }
                }
                return result;
            }
        }
    }
}
