// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
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
using System.Linq;
using System.Text.RegularExpressions;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads in a TopPIC results file (txt format) and creates
    /// a tab-delimited text file with the data.
    /// </summary>
    public class TopPICResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: monoisotopic, Battelle, txt, enums, acetyl, methylation, toppic, syn, fht, Prsm, Da, Pvalue, Frag

        public TopPICResultsProcessor()
        {
            FileDate = "July 10, 2019";
            InitializeLocalVariables();
        }

        public const string TOOL_NAME = "TopPIC";

        public const string FILENAME_SUFFIX_TopPIC_PROTEOFORMS_FILE = "_TopPIC_Proteoforms";

        public const string FILENAME_SUFFIX_TopPIC_PRSMs_FILE = "_TopPIC_PrSMs";

        public const string N_TERMINUS_SYMBOL_TopPIC = ".";

        public const string C_TERMINUS_SYMBOL_TopPIC = ".";

        private const int MAX_ERROR_LOG_LENGTH = 4096;

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
            Precursor_mass = 8,              // Monoisotopic mass value of the observed precursor_mz
            Adjusted_precursor_mass = 9,     // Theoretical monoisotopic mass of the peptide (including mods)
            Proteoform_ID = 10,
            Feature_intensity = 11,
            Feature_score = 12,
            Protein_accession = 13,
            Protein_description = 14,
            First_residue = 15,
            Last_residue = 16,
            Proteoform = 17,
            Unexpected_modifications = 18,
            MIScore = 19,
            Variable_PTMs = 20,
            Matched_peaks = 21,
            Matched_fragment_ions = 22,
            Pvalue = 23,
            Evalue = 24,
            Qvalue = 25,                     //  Spectral FDR, or PepFDR
            Proteoform_QValue = 26
        }

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file created directly by TopPIC
        /// </summary>
        private struct TopPICSearchResult
        {
            public string SpectrumFileName;
            public string Prsm_ID;
            public string Spectrum_ID;
            public string FragMethod;
            public string Scans;
            public int ScanNum;
            public string Peaks;
            public string Charge;
            public short ChargeNum;
            public string Precursor_mass;               // Monoisotopic mass value of the observed precursor_mz
            public string PrecursorMZ;                  // Computed from Precursor_mass
            public string Adjusted_precursor_mass;      // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC
            public string MH;                           // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                         // Computed using Precursor_mass - Adjusted_precursor_mass
            public string DelM_PPM;                     // Computed using DelM and Adjusted_precursor_mass
            public string Proteoform_ID;
            public string Feature_Intensity;
            public string Feature_Score;
            public string Protein;
            public string ProteinDescription;
            public string ResidueStart;                 // First_residue
            public string ResidueEnd;                   // Last_residue
            public string Proteoform;
            public string Unexpected_Mod_Count;         // unexpected modifications
            public string MIScore;
            public string VariablePTMs;
            public string Matched_peaks;
            public string Matched_fragment_ions;
            public string Pvalue;
            public double PValueNum;
            public int RankPValue;
            public string Evalue;
            public double EValueNum;
            public string Qvalue;
            public string Proteoform_QValue;

            public void Clear()
            {
                SpectrumFileName = string.Empty;
                Prsm_ID = string.Empty;
                Spectrum_ID = string.Empty;
                FragMethod = string.Empty;
                Scans = string.Empty;
                ScanNum = 0;
                Peaks = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                Precursor_mass = string.Empty;
                PrecursorMZ = string.Empty;
                Adjusted_precursor_mass = string.Empty;
                MH = string.Empty;
                DelM = string.Empty;
                DelM_PPM = string.Empty;
                Proteoform_ID = string.Empty;
                Feature_Intensity = string.Empty;
                Feature_Score = string.Empty;
                Protein = string.Empty;
                ProteinDescription = string.Empty;
                ResidueStart = string.Empty;
                ResidueEnd = string.Empty;
                Proteoform = string.Empty;
                Unexpected_Mod_Count = string.Empty;
                MIScore = string.Empty;
                VariablePTMs = string.Empty;
                Matched_peaks = string.Empty;
                Matched_fragment_ions = string.Empty;
                Pvalue = string.Empty;
                PValueNum = 0;
                RankPValue = 0;
                Evalue = string.Empty;
                EValueNum = 0;
                Qvalue = string.Empty;
                Proteoform_QValue = string.Empty;
            }
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

            var chMostRecentResidue = NO_RESIDUE;
            var residueLocInPeptide = 0;

            var chAmbiguousResidue = NO_RESIDUE;
            var ambiguousResidueLocInPeptide = 0;

            var clearAmbiguousResidue = false;
            var storeAmbiguousResidue = false;

            var sequence = searchResult.PeptideSequenceWithMods;
            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var chChar = sequence[index];

                if (!parsingModInfo && IsLetterAtoZ(chChar))
                {
                    chMostRecentResidue = chChar;
                    residueLocInPeptide++;

                    if (storeAmbiguousResidue)
                    {
                        chAmbiguousResidue = chChar;
                        ambiguousResidueLocInPeptide = residueLocInPeptide;
                        storeAmbiguousResidue = false;
                    }
                    else if (clearAmbiguousResidue)
                    {
                        chAmbiguousResidue = NO_RESIDUE;
                        clearAmbiguousResidue = false;
                    }

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == ModificationDefinition.ResidueModificationType.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                            if (modificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocInPeptide);

                                searchResult.SearchResultAddModification(
                                    modificationDefinition, chChar, residueLocInPeptide,
                                    residueTerminusState, updateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (chChar == '(')
                {
                    // Start of a mod group
                    storeAmbiguousResidue = true;
                }
                else if (chChar == ')')
                {
                    // End of a mod group
                    clearAmbiguousResidue = true;
                }
                else if (chChar == '[')
                {
                    // Mod Info Start
                    modMassOrName = string.Empty;
                    parsingModInfo = true;
                }
                else if (chChar == ']')
                {
                    // Mod Info End

                    if (!parsingModInfo)
                        continue;

                    char residueForMod;
                    int residueLocForMod;

                    if (chAmbiguousResidue == NO_RESIDUE)
                    {
                        residueForMod = chMostRecentResidue;
                        residueLocForMod = residueLocInPeptide;
                    }
                    else
                    {
                        // Ambiguous mod
                        // We'll associate it with the first residue of the mod group
                        residueForMod = chAmbiguousResidue;
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
                    modMassOrName += chChar;
                }
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
                // Since Inspect allows a terminal peptide residue to be modified twice, we'll allow that to happen,
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
        private void AssignRankAndDeltaNormValues(
            IList<TopPICSearchResult> searchResults,
            int startIndex,
            int endIndex,
            bool sortOnPValue)
        {
            // Duplicate a portion of searchResults so that we can sort by PValue

            var resultsSubset = new Dictionary<int, TopPICSearchResult>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            List<KeyValuePair<int, TopPICSearchResult>> resultsByProbability;

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
                searchResults[entry.Key] = result;
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
        /// Compute the peptide mass
        /// </summary>
        /// <param name="peptide">Sequence with mods, which can be either numeric or named ([15.98154] or [Acetyl])</param>
        /// <param name="totalModMass"></param>
        /// <returns>Mass of the peptide plus totalModMass</returns>
        private double ComputePeptideMass(string peptide, double totalModMass)
        {
            var cleanSequence = GetCleanSequence(peptide);

            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);

            if (Math.Abs(totalModMass) > double.Epsilon)
            {
                mass += totalModMass;
            }

            return mass;
        }

        private static readonly Regex ModMatcher = new(TopPIC_MOD_MASS_OR_NAME_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form [23.5432] or [Acetyl]</param>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            foreach (Match reMatch in ModMatcher.Matches(peptide))
            {
                if (double.TryParse(reMatch.Groups["ModMass"].Value, out var modMassFound))
                {
                    totalModMass += modMassFound;
                }
                else
                {
                    if (LookupModificationMassByName(reMatch.Groups["NamedMod"].Value, out var modMass))
                        totalModMass += modMass;
                }
            }

            return totalModMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var suffixesToFind = new List<string> {
                "_toppic_syn",
                "_toppic_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from TopPIC
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath)
        {
            var columnMapping = new Dictionary<TopPICResultsFileColumns, int>();

            try
            {
                var errorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var headerParsed = false;
                    var dataHasPValues = false;

                    mDeltaMassWarningCount = 0;

                    // Initialize array that will hold all of the records in the TopPIC result file
                    var searchResultsUnfiltered = new List<TopPICSearchResult>();

                    // Initialize the array that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<TopPICSearchResult>();

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

                            // Write the header line
                            WriteSynFHTFileHeader(writer, dataHasPValues, ref errorLog);

                            continue;
                        }

                        var udtSearchResult = new TopPICSearchResult();
                        var validSearchResult = ParseTopPICResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping);

                        if (validSearchResult)
                        {
                            searchResultsUnfiltered.Add(udtSearchResult);
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
                    SortAndWriteFilteredSearchResults(writer, filteredSearchResults, ref errorLog, dataHasPValues);
                }

                // Inform the user if any errors occurred
                if (errorLog.Length > 0)
                {
                    SetErrorMessage("Invalid Lines: \n" + errorLog);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the TopPIC parameter file
        /// </summary>
        /// <param name="topPICParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True on success, false if an error</returns>
        /// <remarks>The DMS-based parameter file for TopPIC uses the same formatting as MS-GF+</remarks>
        private bool ExtractModInfoFromParamFile(
            string topPICParamFilePath,
            out List<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            var modFileProcessor = new MSGFPlusParamFileModExtractor(TOOL_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                topPICParamFilePath,
                MSGFPlusParamFileModExtractor.ModSpecFormats.TopPIC,
                out modInfo);

            if (!success || mErrorCode != PHRPErrorCode.NoError)
            {
                if (mErrorCode == PHRPErrorCode.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the TopPIC parameter file");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modInfo, mPeptideMods);

            return true;
        }

        private new string GetCleanSequence(string sequenceWithMods)
        {
            var primarySequence = GetCleanSequence(sequenceWithMods, out _, out _, out _);
            return primarySequence;
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
                primarySequenceWithMods = string.Copy(sequenceWithMods);
                primarySequence = ModMatcher.Replace(sequenceWithMods, string.Empty);
            }

            // The primary sequence may still contain parentheses; remove them
            return base.GetCleanSequence(primarySequence);
        }

        private void InitializeLocalVariables()
        {
            // Nothing to do at present
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
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that TopPIC synopsis files are normally sorted on PValue, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
            // Although this was a possibility with Inspect, it likely never occurs for TopPIC
            //  But, we'll keep the check in place just in case

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

                // Initialize peptidesFoundForPValueLevel
                var peptidesFoundForPValueLevel = new SortedSet<string>();
                var previousPValue = string.Empty;

                // Assure that pepToProteinMapping is sorted on peptide
                if (pepToProteinMapping.Count > 1)
                {
                    pepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    var errorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
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

                            var validSearchResult = ParseTopPICSynFileEntry(lineIn, searchResult, ref errorLog,
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
                            if (!modsAdded)
                            {
                                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                {
                                    errorLog += "Error adding modifications to sequence at RowIndex '" +
                                                searchResult.ResultID + "'\n";
                                }
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
                                    var currentProtein = string.Copy(searchResult.ProteinName);
                                    do
                                    {
                                        if (pepToProteinMapping[pepToProteinMapIndex].Protein != currentProtein)
                                        {
                                            searchResult.ProteinName = string.Copy(pepToProteinMapping[pepToProteinMapIndex].Protein);
                                            SaveResultsFileEntrySeqInfo(searchResult, false);
                                        }

                                        pepToProteinMapIndex++;
                                    } while (pepToProteinMapIndex < pepToProteinMapping.Count &&
                                             currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                                }
                                else
                                {
                                    // Match not found; this is unexpected
                                    ReportWarning("no match for '" + currentPeptideWithMods + "' in pepToProteinMapping");
                                }
                            }

                            // Update the progress
                            UpdateSynopsisFileCreationProgress(reader);
                        }
                    }

                    if (CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(inputFilePath);
                        var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                        SaveModificationSummaryFile(modificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (errorLog.Length > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + errorLog);
                    }

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
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
                SetErrorMessage(ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parses an entry from the TopPIC results file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICResultsFileEntry(
            string lineIn,
            ref TopPICSearchResult udtSearchResult,
            ref string errorLog,
            IDictionary<TopPICResultsFileColumns, int> columnMapping)
        {
            string[] splitLine = null;

            double precursorMZ = 0;

            try
            {
                udtSearchResult.Clear();
                splitLine = lineIn.TrimEnd().Split('\t');

                // The file should have over 20 columns, but we'll only require 15
                if (splitLine.Length < 15)
                {
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);
                if (!GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Prsm_ID], out udtSearchResult.Prsm_ID))
                {
                    ReportError("Prsm_ID column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Spectrum_ID], out udtSearchResult.Spectrum_ID);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                if (!GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Scans], out udtSearchResult.Scans))
                {
                    ReportError("Scan(s) column is missing or invalid", true);
                }

                if (!int.TryParse(udtSearchResult.Scans, out udtSearchResult.ScanNum))
                {
                    // .Scans may have a list of scan numbers; extract the first scan number from .scans
                    var scanNumberDigits = string.Empty;
                    foreach (var chChar in udtSearchResult.Scans)
                    {
                        if (char.IsDigit(chChar))
                        {
                            scanNumberDigits += chChar;
                        }
                    }

                    if (!int.TryParse(scanNumberDigits, out udtSearchResult.ScanNum))
                    {
                        ReportWarning("Error parsing out the scan number from the scan list; could not find an integer: " +
                                      udtSearchResult.Scans);
                        udtSearchResult.ScanNum = 0;
                    }
                }

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Peaks], out udtSearchResult.Peaks);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Monoisotopic mass value of the observed precursor_mz
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Precursor_mass], out udtSearchResult.Precursor_mass);

                // precursorMonoMass is Observed m/z, converted to monoisotopic mass
                if (double.TryParse(udtSearchResult.Precursor_mass, out var precursorMonoMass))
                {
                    if (udtSearchResult.ChargeNum > 0)
                    {
                        precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, udtSearchResult.ChargeNum);
                        udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                    }
                }

                // peptideMonoMassTopPIC is Theoretical peptide monoisotopic mass, including mods, as computed by TopPIC
                double peptideMonoMassTopPIC;

                if (columnMapping[TopPICResultsFileColumns.Adjusted_precursor_mass] >= 0)
                {
                    // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC
                    GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Adjusted_precursor_mass], out udtSearchResult.Adjusted_precursor_mass);

                    double.TryParse(udtSearchResult.Adjusted_precursor_mass, out peptideMonoMassTopPIC);
                }
                else
                {
                    peptideMonoMassTopPIC = 0;
                }

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform_ID], out udtSearchResult.Proteoform_ID);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Feature_intensity], out udtSearchResult.Feature_Intensity);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Feature_score], out udtSearchResult.Feature_Score);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Protein_accession], out udtSearchResult.Protein);
                udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Protein_description], out udtSearchResult.ProteinDescription);

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.First_residue], out udtSearchResult.ResidueStart);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Last_residue], out udtSearchResult.ResidueEnd);

                if (!GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform], out udtSearchResult.Proteoform))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                // Add the standard terminus symbols to the peptide sequence
                udtSearchResult.Proteoform = ReplaceTerminus(udtSearchResult.Proteoform);

                // Parse the sequence to determine the total mod mass
                // Note that we do not remove any of the mod symbols since TopPIC since mods can ambiguously apply to residues
                var totalModMass = ComputeTotalModMass(udtSearchResult.Proteoform);

                // Compute theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                var peptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Proteoform, totalModMass);

                if (Math.Abs(peptideMonoMassTopPIC) < double.Epsilon)
                {
                    peptideMonoMassTopPIC = peptideMonoMassPHRP;
                }

                // Warn the user if the monoisotopic mass values differ by more than 0.1 Da
                ValidateMatchingMonoisotopicMass(TOOL_NAME, udtSearchResult.Proteoform, peptideMonoMassPHRP, peptideMonoMassTopPIC, ref mDeltaMassWarningCount);

                if (peptideMonoMassTopPIC > 0)
                {
                    // Compute DelM and DelM_PPM
                    var delM = precursorMonoMass - peptideMonoMassTopPIC;
                    udtSearchResult.DelM = MassErrorToString(delM);

                    if (precursorMZ > 0)
                    {
                        udtSearchResult.DelM_PPM =
                            PRISM.StringUtilities.DblToString(PeptideMassCalculator.MassToPPM(delM, precursorMZ), 5, 0.00005);
                    }
                    else
                    {
                        udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(PeptideMassCalculator.MassToPPM(delM, 1000), 5, 0.00005);
                    }
                }

                // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoMassPHRP, 0), 6);

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Unexpected_modifications], out udtSearchResult.Unexpected_Mod_Count);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.MIScore], out udtSearchResult.MIScore);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Variable_PTMs], out udtSearchResult.VariablePTMs);

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Matched_peaks], out udtSearchResult.Matched_peaks);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Matched_fragment_ions], out udtSearchResult.Matched_fragment_ions);
                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Pvalue], out udtSearchResult.Pvalue);
                if (!double.TryParse(udtSearchResult.Pvalue, out udtSearchResult.PValueNum))
                    udtSearchResult.PValueNum = 0;

                // Assure that the following are truly integers (Matched_peaks and Matched_fragment_ions are often of the form 8.0)
                udtSearchResult.Unexpected_Mod_Count = AssureInteger(udtSearchResult.Unexpected_Mod_Count, 0);      // Unexpected_Mod_Count
                udtSearchResult.Peaks = AssureInteger(udtSearchResult.Peaks, 0);                                    // Peak_count
                udtSearchResult.Matched_peaks = AssureInteger(udtSearchResult.Matched_peaks, 0);                    // Matched_Peak_Count
                udtSearchResult.Matched_fragment_ions = AssureInteger(udtSearchResult.Matched_fragment_ions, 0);    // Matched_Fragment_Ion_Count

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Evalue], out udtSearchResult.Evalue);
                if (!double.TryParse(udtSearchResult.Evalue, out udtSearchResult.EValueNum))
                    udtSearchResult.EValueNum = 0;

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Qvalue], out udtSearchResult.Qvalue);

                if (string.Equals(udtSearchResult.Qvalue, "infinity", StringComparison.OrdinalIgnoreCase))
                {
                    udtSearchResult.Qvalue = "10";
                }
                else if (!string.IsNullOrEmpty(udtSearchResult.Qvalue) && !double.TryParse(udtSearchResult.Qvalue, out _))
                {
                    udtSearchResult.Qvalue = string.Empty;
                }

                GetColumnValue(splitLine, columnMapping[TopPICResultsFileColumns.Proteoform_QValue], out udtSearchResult.Proteoform_QValue);

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the TopPIC results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorLog += "Error parsing TopPIC Results in ParseTopPICResultsFileEntry for RowIndex '" + splitLine[0] + "'" +
                                       "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing TopPIC Results in ParseTopPICResultsFileEntry\n";
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a TopPIC results file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICResultsFileHeaderLine(string lineIn, IDictionary<TopPICResultsFileColumns, int> columnMapping)
        {
            // Header prior to November 2018:
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Protein name    First residue    Last residue    Proteoform    #unexpected modifications    #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR    #Variable PTMs

            // Header for TopPIC 1.2 and 1.3
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Protein accession    Protein description    First residue    Last residue    Proteoform    #unexpected modifications    MIScore    #variable PTMs    #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR

            // Header for TopPIC 1.4
            // The P-Value column has been removed and the Q-Value and "Proteoform FDR" columns have been renamed
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Feature score    Protein accession    Protein description    First residue    Last residue    Proteoform    #unexpected modifications    MIScore    #variable PTMs    #matched peaks    #matched fragment ions    E-value    Spectrum-level Q-value    Proteoform-level Q-value

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
                {"Protein name", TopPICResultsFileColumns.Protein_accession},
                {"Protein accession", TopPICResultsFileColumns.Protein_accession},
                {"Protein description", TopPICResultsFileColumns.Protein_description},
                {"First residue", TopPICResultsFileColumns.First_residue},
                {"Last residue", TopPICResultsFileColumns.Last_residue},
                {"Proteoform", TopPICResultsFileColumns.Proteoform},
                {"#unexpected modifications", TopPICResultsFileColumns.Unexpected_modifications},
                {"MIScore", TopPICResultsFileColumns.MIScore},
                {"#variable PTMs", TopPICResultsFileColumns.Variable_PTMs},
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
                for (var index = 0; index <= splitLine.Length - 1; index++)
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
                SetErrorMessage("Error parsing header in TopPIC results file: " + ex.Message);
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
                for (var index = 0; index <= splitLine.Length - 1; index++)
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
                SetErrorMessage("Error parsing header in TopPIC synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from the TopPIC Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequenceWithMods"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseTopPICSynFileEntry(
            string lineIn,
            TopPICResults searchResult,
            ref string errorLog,
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

                if (!GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from TopPIC Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from TopPIC Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.DelM], out string msAlignComputedDelM);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.DelMPPM], out string msAlignComputedDelMppm);

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

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Prsm_ID], out string prsmId);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Spectrum_ID], out string spectrumId);

                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.MH], out string parentIonMH);

                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Unexpected_Mod_Count], out string unexpectedModCount);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Peak_Count], out string peakCount);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Matched_Peak_Count], out string matchedPeakCount);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Matched_Fragment_Ion_Count], out string matchedFragmentIonCount);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.PValue], out string pValue);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Rank_PValue], out string rankPValue);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.EValue], out string eValue);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.QValue], out string qValue);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.FragMethod], out string fragMethod);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.Proteoform_QValue], out string proteoformFDR);
                GetColumnValue(splitLine, columnMapping[TopPICSynFileColumns.VariablePTMs], out string variablePTMs);

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
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorLog += "Error parsing TopPIC Results for RowIndex '" + splitLine[0] + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing TopPIC Results in ParseTopPICSynFileEntry\n";
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Main processing function
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

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);

                    var pepToProteinMapping = new List<PepToProteinMapping>();

                    // Load the TopPIC Parameter File so that we can determine whether Cysteine residues are statically modified
                    var modInfoExtracted = ExtractModInfoFromParamFile(SearchToolParameterFilePath, out var topPicModInfo);
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

                    // Do not create a first-hits file for TopPIC results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_msalign_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    success = ParseTopPICSynopsisFile(synOutputFilePath, outputDirectoryPath, ref pepToProteinMapping, false);

                    // Remove all items from pepToProteinMapping to reduce memory overhead
                    pepToProteinMapping.Clear();
                    pepToProteinMapping.TrimExcess();

                    if (success && CreateProteinModsFile)
                    {
                        success = CreateProteinModsFileWork(baseName, inputFile, synOutputFilePath, outputDirectoryPath);
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

        private bool CreateProteinModsFileWork(string baseName, FileInfo inputFile, string synOutputFilePath, string outputDirectoryPath)
        {
            bool success;

            // Create the MTSPepToProteinMap file

            var baseNameFilePath = Path.Combine(inputFile.DirectoryName ?? string.Empty, baseName);
            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseNameFilePath, outputDirectoryPath, mts: true);

            var sourcePHRPDataFiles = new List<string>();

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
                if (File.Exists(mtsPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
                {
                    success = true;
                }
                else
                {
                    // Auto-change mIgnorePeptideToProteinMapperErrors to True
                    // We do this since a small number of peptides reported by TopPIC don't perfectly match the fasta file
                    IgnorePeptideToProteinMapperErrors = true;
                    success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                    if (!success)
                    {
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (success)
            {
                if (inputFile.Directory == null)
                {
                    ReportWarning("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                }
                else if (string.IsNullOrWhiteSpace(synOutputFilePath))
                {
                    ReportWarning("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputDirectoryPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          PeptideHitResultTypes.TopPIC);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                success = true;
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
            TextWriter writer,
            IEnumerable<TopPICSearchResult> filteredSearchResults,
            ref string errorLog,
            bool dataHasPValues)
        {
            IOrderedEnumerable<TopPICSearchResult> query;

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
                WriteSearchResultToFile(index, writer, result, dataHasPValues, ref errorLog);
                index++;
            }
        }

        private void StoreSynMatches(
            IList<TopPICSearchResult> searchResults,
            int startIndex,
            int endIndex,
            ICollection<TopPICSearchResult> filteredSearchResults,
            bool dataHasPValues)
        {
            AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex, dataHasPValues);

            // The calling procedure already sorted by scan, charge, and PValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (dataHasPValues && searchResults[index].PValueNum <= MSAlignAndTopPICSynopsisFilePValueThreshold ||
                    !dataHasPValues && searchResults[index].EValueNum <= MSAlignAndTopPICSynopsisFilePValueThreshold)
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
        /// <param name="errorLog"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            bool dataHasPValues,
            ref string errorLog)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = TopPICSynFileReader.GetColumnHeaderNamesAndIDs(false);

                if (!dataHasPValues)
                {
                    headerColumns.Remove(TopPICSynFileReader.DATA_COLUMN_PValue);
                    headerColumns.Remove(TopPICSynFileReader.DATA_COLUMN_Rank_PValue);
                }

                var headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();

                writer.WriteLine(CollapseList(headerNames));
            }
            catch (Exception)
            {
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error writing synopsis / first hits header\n";
                }
            }
        }

        /// <summary>
        /// Writes an entry to the synopsis file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="dataHasPValues"></param>
        /// <param name="errorLog"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            TopPICSearchResult udtSearchResult,
            bool dataHasPValues,
            ref string errorLog)
        {
            try
            {
                // Primary Columns
                //
                // TopPIC
                // ResultID  Scan  Prsm_ID  Spectrum_ID  FragMethod  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Proteoform_ID  Feature_Intensity  Feature_Score  Protein  ResidueStart  ResidueEnd  Unexpected_Mod_Count  Peak_Count  Matched_Peak_Count  Matched_Fragment_Ion_Count  PValue  Rank_PValue  EValue  QValue  ProteoformFDR  VariablePTMs

                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.ScanNum.ToString(),
                    udtSearchResult.Prsm_ID,
                    udtSearchResult.Spectrum_ID,
                    udtSearchResult.FragMethod,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.DelM,
                    udtSearchResult.DelM_PPM,
                    udtSearchResult.MH,
                    udtSearchResult.Proteoform,                // aka Peptide
                    udtSearchResult.Proteoform_ID,
                    udtSearchResult.Feature_Intensity,
                    udtSearchResult.Feature_Score,
                    udtSearchResult.Protein,
                    udtSearchResult.ResidueStart,
                    udtSearchResult.ResidueEnd,
                    udtSearchResult.Unexpected_Mod_Count,       // Unexpected_Mod_Count
                    udtSearchResult.Peaks,                      // Peak_count
                    udtSearchResult.Matched_peaks,              // Matched_Peak_Count
                    udtSearchResult.Matched_fragment_ions       // Matched_Fragment_Ion_Count
                };

                if (dataHasPValues)
                {
                    data.Add(udtSearchResult.Pvalue);
                    data.Add(udtSearchResult.RankPValue.ToString());
                }

                data.Add(udtSearchResult.Evalue);
                data.Add(udtSearchResult.Qvalue);
                data.Add(udtSearchResult.Proteoform_QValue);
                data.Add(udtSearchResult.VariablePTMs);

                writer.WriteLine(CollapseList(data));
            }
            catch (Exception)
            {
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error writing synopsis / first hits record\n";
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

        private class TopPICSearchResultsComparerScanChargePValuePeptide : IComparer<TopPICSearchResult>
        {
            public int Compare(TopPICSearchResult x, TopPICSearchResult y)
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

                // Charge is the same; check Pvalue
                var result = string.CompareOrdinal(x.Pvalue, y.Pvalue);
                if (result == 0)
                {
                    // Pvalue is the same; check peptide
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
