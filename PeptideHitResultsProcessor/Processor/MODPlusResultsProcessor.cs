// This class reads in an MODPlus results file (.txt format) and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 05/15/2015
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Processor
{
    public class MODPlusResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: MODa, tda, methylation, udt, fht, modp, ModDefs, massdiff

        public MODPlusResultsProcessor()
        {
            FileDate = "July 10, 2019";
        }

        public const string TOOL_NAME = "ModPlus";

        public const string FILENAME_SUFFIX_MODPlus_FILE = "_modp.id";

        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_MODPlus = "-";

        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_MODPlus = "-";

        // This is used for filtering both MODa and MODPlus results
        public const float DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD = 0.05f;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string MODPlus_MOD_MASS_REGEX = "([+-][0-9.]+)";

        private const byte MODPlus_MASS_DIGITS_OF_PRECISION = 3;

        private const byte MODPlus_MASS_DIGITS_OF_PRECISION_LOOSE = 2;

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to the tab-delimited file (_MODPlus.id.txt) created by MODPlus's tda_plus.jar file
        /// </summary>
        private enum MODPlusResultsFileColumns
        {
            SpectrumFileName = 0,
            SpectrumIndex = 1,
            ScanNumber = 2,
            ObservedMonoMass = 3,
            Charge = 4,
            CalculatedMonoMass = 5,
            DeltaMass = 6,
            Score = 7,
            Probability = 8,
            Peptide = 9,
            NTT = 10,
            ProteinAndPeptidePositionList = 11,
            ModificationAnnotation = 12
        }

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file (_MODPlus.id.txt) created by MODPlus's tda_plus.jar file
        /// </summary>
        private struct MODPlusSearchResult
        {
            public string SpectrumFileName;
            public string SpectrumIndex;
            public int ScanNum;
            public string Precursor_mass;           // Uncharged monoisotopic mass value of the observed precursor_mz, reported as ObservedMW by MODPlus
            public string PrecursorMZ;              // Computed by this class from ObservedMonoMass
            public string Charge;
            public short ChargeNum;
            public string CalculatedMonoMass;       // Theoretical monoisotopic mass of the peptide (including mods), as computed by MODPlus
            public string DeltaMass;                // Computed by MODPlus
            public string MH;                       // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                     // Computed by this class using Precursor_mass - CalculatedMonoMass
            public string DelM_PPM;                 // Computed by this class using DelM and CalculatedMonoMass
            public string Score;
            public double ScoreNum;
            public string Probability;              // Higher values are better
            public double ProbabilityNum;           // Higher values are better
            public int RankScore;
            public string Peptide;
            public string NTT;                      // Number of Tryptic Termini
            public string ProteinList;              // One or more protein entries of the form ref|YP_003651515.1[K.196~206.Q(2)] where the text in brackets is the start/stop residues of the peptide; Multiple entries will be separated by semicolons, e.g. ref|YP_003651515.1[K.196~206.Q(2)];ref|YP_003201491.1[K.223~233.Q(2)];ref|YP_003313784.1[K.266~276.Q(2)]
            public string ModificationAnnotation;
            public double FDR;                      // Computed by this class
            public double QValue;                   // Computed by this class

            public void Clear()
            {
                SpectrumFileName = string.Empty;
                SpectrumIndex = string.Empty;
                ScanNum = 0;
                Precursor_mass = string.Empty;
                PrecursorMZ = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                CalculatedMonoMass = string.Empty;
                DeltaMass = string.Empty;
                MH = string.Empty;
                DelM = string.Empty;
                DelM_PPM = string.Empty;
                Score = string.Empty;
                ScoreNum = 0;
                Probability = string.Empty;
                ProbabilityNum = 0;
                RankScore = 0;
                Peptide = string.Empty;
                NTT = string.Empty;
                ProteinList = string.Empty;
                ModificationAnnotation = string.Empty;
                FDR = 0;
                QValue = 0;
            }

            public override string ToString()
            {
                return Probability + ": " + Peptide;
            }
        }

        private int mDeltaMassWarningCount;

        private Regex mProteinNamePositionSplit;

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

            var parsingModMass = false;

            var modMassDigits = string.Empty;

            var mostRecentResidue = NO_RESIDUE;
            var residueLocInPeptide = 0;

            var sequence = searchResult.PeptideSequenceWithMods;
            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var chChar = sequence[index];

                if (IsLetterAtoZ(chChar))
                {
                    if (parsingModMass)
                    {
                        // Associate the mod mass in modMassDigits with the previous residue
                        AssociateDynamicModWithResidue(searchResult, mostRecentResidue, residueLocInPeptide, modMassDigits, updateModOccurrenceCounts);
                        parsingModMass = false;
                    }

                    mostRecentResidue = chChar;
                    residueLocInPeptide++;

                    // Look for static mods to associate with this residue
                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ResidueModificationType.StaticMod)
                            continue;

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
                else
                {
                    var isNumberChar = chChar == '+' || chChar == '-' || char.IsDigit(chChar);

                    if (parsingModMass)
                    {
                        if (isNumberChar || chChar == '.')
                        {
                            modMassDigits += chChar;
                        }
                    }
                    else if (isNumberChar)
                    {
                        // Mod Mass Start
                        modMassDigits = chChar.ToString();
                        parsingModMass = true;
                    }
                    else
                    {
                        // Unrecognized symbol; ignore it
                    }
                }
            }

            if (parsingModMass)
            {
                // Associate the mod mass in modMassDigits with the previous residue
                AssociateDynamicModWithResidue(searchResult, mostRecentResidue, residueLocInPeptide, modMassDigits, updateModOccurrenceCounts);
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

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
        }

        private void AssociateDynamicModWithResidue(
            SearchResultsBaseClass searchResult,
            char chMostRecentResidue,
            int residueLocInPeptide,
            string modMassDigits,
            bool updateModOccurrenceCounts)
        {
            var residueForMod = chMostRecentResidue;
            var residueLocForMod = residueLocInPeptide;

            if (double.TryParse(modMassDigits, out var modMass))
            {
                if (residueLocForMod == 0)
                {
                    // Modification is at the peptide N-terminus
                    residueLocForMod = 1;
                }

                var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocForMod);

                var success = searchResult.SearchResultAddModification(
                    modMass, residueForMod, residueLocForMod,
                    residueTerminusState, updateModOccurrenceCounts,
                    MODPlus_MASS_DIGITS_OF_PRECISION, MODPlus_MASS_DIGITS_OF_PRECISION_LOOSE);

                if (!success)
                {
                    var errorMessage = searchResult.ErrorMessage;
                    if (string.IsNullOrEmpty(errorMessage))
                    {
                        errorMessage = "SearchResultAddDynamicModification returned false for mod mass " + modMassDigits;
                    }
                    SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
                }
            }
        }

        /// <summary>
        /// Ranks each entry assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        private void AssignRankAndDeltaNormValues(
            IList<MODPlusSearchResult> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of searchResults so that we can sort by descending Probability

            var resultsSubset = new Dictionary<int, MODPlusSearchResult>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsByScore = (from item in resultsSubset orderby item.Value.ScoreNum descending select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByScore)
            {
                var result = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = result.ScoreNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(result.ScoreNum - lastValue) > double.Epsilon)
                    {
                        lastValue = result.ScoreNum;
                        currentRank++;
                    }
                }

                result.RankScore = currentRank;
                searchResults[entry.Key] = result;
            }
        }

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

        private static readonly Regex ModMassMatcher = new(MODPlus_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptide, out var primarySequence, out var _, out var _);

            // Parse the dynamic mods reported by MODPlus
            foreach (Match reMatch in ModMassMatcher.Matches(primarySequence))
            {
                // We use .TrimEnd() because the matched mod mass will end in a period if this mod applies to the final residue in a peptide
                if (double.TryParse(reMatch.Groups[1].Value.TrimEnd('.'), out var modMassFound))
                {
                    totalModMass += modMassFound;
                }
            }

            // Now look for static mods
            // First determine the index of the last residue in primarySequence
            var indexLastChar = primarySequence.Length;

            for (var index = primarySequence.Length - 1; index >= 0; index += -1)
            {
                if (IsLetterAtoZ(primarySequence[index]))
                {
                    indexLastChar = index;
                    break;
                }
            }

            for (var index = 0; index <= primarySequence.Length - 1; index++)
            {
                var chChar = primarySequence[index];

                if (!IsLetterAtoZ(chChar))
                    continue;

                // Look for static mods to associate with this residue
                for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                {
                    if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ResidueModificationType.StaticMod)
                        continue;

                    var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                    var matchFound = modificationDefinition.TargetResiduesContain(chChar);

                    if (!matchFound && index == 0)
                    {
                        matchFound = modificationDefinition.TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
                    }

                    if (!matchFound && index == indexLastChar)
                    {
                        matchFound = modificationDefinition.TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);
                    }

                    if (matchFound)
                    {
                        totalModMass += modificationDefinition.ModificationMass;
                    }
                }
            }

            return totalModMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var suffixesToFind = new List<string> {
                "_MODPlus_syn",
                "_MODPlus_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MODPlus
        /// The synopsis file includes every result with a probability above a set threshold
        /// The first-hits file includes the result with the highest probability (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath)
        {
            try
            {
                var columnMapping = new Dictionary<MODPlusResultsFileColumns, int>();
                var errorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var headerParsed = false;
                    mDeltaMassWarningCount = 0;

                    // Initialize the list that will hold all of the records in the MODPlus result file
                    var searchResultsUnfiltered = new List<MODPlusSearchResult>();

                    // Initialize the list that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<MODPlusSearchResult>();

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
                            // Parse the header line

                            var success = ParseMODPlusResultsFileHeaderLine(lineIn, columnMapping);
                            if (!success)
                            {
                                if (string.IsNullOrEmpty(mErrorMessage))
                                {
                                    SetErrorMessage("Invalid header line in " + Path.GetFileName(inputFilePath));
                                }

                                return false;
                            }

                            // Write the header line to the output file
                            WriteSynFHTFileHeader(writer, ref errorLog);

                            headerParsed = true;
                            continue;
                        }

                        var udtSearchResult = new MODPlusSearchResult();

                        var validSearchResult = ParseMODPlusResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping);

                        if (validSearchResult)
                        {
                            searchResultsUnfiltered.Add(udtSearchResult);
                        }

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    // Sort the SearchResults by scan, charge, and descending score
                    searchResultsUnfiltered.Sort(new MODPlusSearchResultsComparerScanChargeScorePeptide());

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
                        StoreSynMatches(searchResultsUnfiltered, startIndex, endIndex, filteredSearchResults);

                        startIndex = endIndex + 1;
                    }

                    // Sort the data in filteredSearchResults then write out to disk
                    SortAndWriteFilteredSearchResults(writer, filteredSearchResults, ref errorLog);
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
                SetErrorMessage("Error in CreateSynResultsFile: " + ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Load the static mods defined in the MODPlus parameter file
        /// </summary>
        /// <param name="modPlusParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>
        /// We don't care about the dynamic mods because there are so many possible mods.
        /// We'll add each dynamic mod as we encounter it in the results
        /// </remarks>
        private bool ExtractModInfoFromMODPlusParamFile(string modPlusParamFilePath, out List<ModificationDefinition> modInfo)
        {
            modInfo = new List<ModificationDefinition>();

            try
            {
                if (string.IsNullOrEmpty(modPlusParamFilePath))
                {
                    SetErrorMessage("MODPlus Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                var paramFile = new FileInfo(modPlusParamFilePath);
                if (!paramFile.Exists)
                {
                    SetErrorMessage("MODPlus param file not found: " + modPlusParamFilePath);
                    return false;
                }

                // Read the contents of the parameter file
                var doc = new XmlDocument();
                doc.Load(paramFile.FullName);

                var nodeList = doc.SelectNodes("/search/modifications/fixed/mod");
                if (nodeList?.Count > 0)
                {
                    // Store the fixed mods

                    foreach (XmlNode node in nodeList)
                    {
                        var modName = node.Attributes?["name"].Value;
                        var residue = node.Attributes?["site"].Value.Trim();
                        var modPosition = node.Attributes?["position"].Value;
                        var modMass = node.Attributes?["massdiff"].Value;

                        // Replace N-Term or C-Term with < or >
                        if (string.Equals(residue, "n-term", StringComparison.OrdinalIgnoreCase))
                            residue = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        if (string.Equals(residue, "c-term", StringComparison.OrdinalIgnoreCase))
                            residue = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                        if (double.TryParse(modMass, out var modMassDa))
                        {
                            if (Math.Abs(modMassDa - 0) > float.Epsilon)
                            {
                                var massCorrectionTag = mPeptideMods.LookupMassCorrectionTagByMass(modMassDa);

                                var modType = ModificationDefinition.ResidueModificationType.StaticMod;
                                if (residue == AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || residue == AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                                {
                                    modType = ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod;
                                }

                                var modDef = new ModificationDefinition(ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, residue, modType, massCorrectionTag);
                                modInfo.Add(modDef);
                            }
                        }
                    }
                }

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MODPlus parameter file (" + Path.GetFileName(modPlusParamFilePath) + "): " + ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
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
        private bool ParseMODPlusSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that MODPlus synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            var columnMapping = new Dictionary<MODPlusSynFileColumns, int>();

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
                var searchResult = new MODPlusResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize peptidesFoundForProbabilityLevel
                var peptidesFoundForProbabilityLevel = new SortedSet<string>();

                var previousProbability = string.Empty;

                var errorLog = string.Empty;

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

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
                                var validHeader = ParseMODPlusSynFileHeaderLine(lineIn, columnMapping);
                                if (!validHeader)
                                {
                                    // Error parsing header
                                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMODPlusSynFileEntry(
                                lineIn, searchResult, ref errorLog, resultsProcessed, columnMapping, out _);

                            resultsProcessed++;
                            if (!validSearchResult)
                            {
                                continue;
                            }

                            var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                            bool firstMatchForGroup;
                            if (searchResult.Probability == previousProbability)
                            {
                                // New result has the same Probability as the previous result
                                // See if peptidesFoundForProbabilityLevel contains the peptide, scan and charge

                                if (peptidesFoundForProbabilityLevel.Contains(key))
                                {
                                    firstMatchForGroup = false;
                                }
                                else
                                {
                                    peptidesFoundForProbabilityLevel.Add(key);
                                    firstMatchForGroup = true;
                                }
                            }
                            else
                            {
                                // New Probability
                                // Reset peptidesFoundForProbabilityLevel
                                peptidesFoundForProbabilityLevel.Clear();

                                // Update previousProbability
                                previousProbability = searchResult.Probability;

                                // Append a new entry to peptidesFoundForProbabilityLevel
                                peptidesFoundForProbabilityLevel.Add(key);
                                firstMatchForGroup = true;
                            }

                            var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                            if (!modsAdded)
                            {
                                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                {
                                    errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'" +
                                                   "\n";
                                }
                            }

                            SaveResultsFileEntrySeqInfo(searchResult, firstMatchForGroup);

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
                    SetErrorMessage("Error reading input file in ParseMODPlusSynopsisFile: " + ex.Message);
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
                SetErrorMessage("Error creating the output file in ParseMODPlusSynopsisFile: " + ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parse a MODPlus results line while creating the MODPlus synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODPlusResultsFileEntry(
            string lineIn,
            ref MODPlusSearchResult udtSearchResult,
            ref string errorLog,
            IDictionary<MODPlusResultsFileColumns, int> columnMapping)
        {
            // Parses an entry from the MODPlus results file

            var rowIndex = "?";

            try
            {
                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 11)
                {
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);

                if (!GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.SpectrumIndex], out udtSearchResult.SpectrumIndex))
                {
                    ReportError("Index column is missing or invalid", true);
                }
                else
                {
                    rowIndex = udtSearchResult.SpectrumIndex;
                }

                if (!int.TryParse(udtSearchResult.SpectrumIndex, out _))
                {
                    ReportError("Index column is not numeric", true);
                }

                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.ScanNumber], out udtSearchResult.ScanNum);

                // Monoisotopic mass value of the observed precursor_mz
                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.ObservedMonoMass], out udtSearchResult.Precursor_mass);
                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // precursorMonoMass is Observed m/z, converted to monoisotopic mass
                if (double.TryParse(udtSearchResult.Precursor_mass, out var precursorMonoMass))
                {
                    if (udtSearchResult.ChargeNum > 0)
                    {
                        var precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, udtSearchResult.ChargeNum);
                        udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                    }
                }

                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                // Theoretical peptide monoisotopic mass, including mods, as computed by MODPlus
                double.TryParse(udtSearchResult.CalculatedMonoMass, out var peptideMonoMassMODPlus);

                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.DeltaMass], out udtSearchResult.DeltaMass);
                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.Score], out udtSearchResult.Score);
                if (!double.TryParse(udtSearchResult.Score, out udtSearchResult.ScoreNum))
                    udtSearchResult.ScoreNum = 0;

                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.Probability], out udtSearchResult.Probability);
                if (!double.TryParse(udtSearchResult.Probability, out udtSearchResult.ProbabilityNum))
                    udtSearchResult.ProbabilityNum = 0;

                if (!GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.NTT], out udtSearchResult.NTT);

                if (splitLine.Length > (int)MODPlusResultsFileColumns.ProteinAndPeptidePositionList)
                {
                    GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.ProteinAndPeptidePositionList], out udtSearchResult.ProteinList);

                    // The protein column will have both the protein name and the peptide position
                    // For example, ref|YP_001038741.1[R.67~78.L(2)]
                    // It may have multiple proteins listed, separated by semicolons
                    // We will split the list on semicolons in function ParseMODPlusSynFileEntry

                    if (!udtSearchResult.ProteinList.Contains('['))
                    {
                        // This is likely a reverse-hit protein
                        udtSearchResult.ModificationAnnotation = string.Copy(udtSearchResult.ProteinList);
                        udtSearchResult.ProteinList = string.Empty;
                    }
                    else
                    {
                        if (splitLine.Length > (int)MODPlusResultsFileColumns.ModificationAnnotation)
                        {
                            GetColumnValue(splitLine, columnMapping[MODPlusResultsFileColumns.ModificationAnnotation], out udtSearchResult.ModificationAnnotation);
                        }
                    }
                }

                // Parse the sequence to determine the total mod mass
                // Note that we do not remove any of the mod symbols since MODPlus identifies mods by mass alone
                // Note that static mods are implied (thus are not explicitly displayed by MODPlus)
                var totalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                // Compute the theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                var peptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, totalModMass);

                // Only override peptideMonoMassMODPlus if it is 0
                if (Math.Abs(peptideMonoMassMODPlus) < double.Epsilon)
                {
                    peptideMonoMassMODPlus = peptideMonoMassPHRP;
                }

                // Warn the user if the monoisotopic mass values differ by more than 0.1 Da
                ValidateMatchingMonoisotopicMass(TOOL_NAME, udtSearchResult.Peptide, peptideMonoMassPHRP, peptideMonoMassMODPlus, ref mDeltaMassWarningCount);

                if (peptideMonoMassMODPlus > 0)
                {
                    // Compute DelM and DelM_PPM
                    var delM = precursorMonoMass - peptideMonoMassMODPlus;
                    udtSearchResult.DelM = MassErrorToString(delM);

                    var peptideDeltaMassCorrectedPpm =
                        SearchResultsBaseClass.ComputeDelMCorrectedPPM(delM, precursorMonoMass, true, peptideMonoMassMODPlus);

                    udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(peptideDeltaMassCorrectedPpm, 5, 0.00005);
                }

                // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoMassPHRP, 0), 6);

                if (string.Equals(udtSearchResult.Probability, "infinity", StringComparison.OrdinalIgnoreCase))
                {
                    udtSearchResult.Probability = "0";
                }
                else if (!string.IsNullOrEmpty(udtSearchResult.Probability) && !double.TryParse(udtSearchResult.Probability, out _))
                {
                    udtSearchResult.Probability = string.Empty;
                }

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the MODPlus results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        errorLog += "Error parsing MODPlus Results in ParseMODPlusResultsFileEntry for RowIndex '" + rowIndex + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MODPlus Results in ParseMODPlusResultsFileEntry\n";
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a MOD+ results file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseMODPlusResultsFileHeaderLine(string lineIn, IDictionary<MODPlusResultsFileColumns, int> columnMapping)
        {
            // The expected column order from MODPlus:
            //   SpectrumFile   Index   ScanNo   ObservedMW   Charge   CalculatedMW   DeltaMass   Score   Probability   Peptide   NTT    Protein   ModificationAnnotation

            var columnNames = new SortedDictionary<string, MODPlusResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"SpectrumFile", MODPlusResultsFileColumns.SpectrumFileName},
                {"Index", MODPlusResultsFileColumns.SpectrumIndex},
                {"ScanNo", MODPlusResultsFileColumns.ScanNumber},
                {"ObservedMW", MODPlusResultsFileColumns.ObservedMonoMass},
                {"Charge", MODPlusResultsFileColumns.Charge},
                {"CalculatedMW", MODPlusResultsFileColumns.CalculatedMonoMass},
                {"DeltaMass", MODPlusResultsFileColumns.DeltaMass},
                {"Score", MODPlusResultsFileColumns.Score},
                {"Probability", MODPlusResultsFileColumns.Probability},
                {"Peptide", MODPlusResultsFileColumns.Peptide},
                {"NTT", MODPlusResultsFileColumns.NTT},
                {"Protein", MODPlusResultsFileColumns.ProteinAndPeptidePositionList},
                {"ModificationAnnotation", MODPlusResultsFileColumns.ModificationAnnotation}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MODPlusResultsFileColumns resultColumn in Enum.GetValues(typeof(MODPlusResultsFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');
                var useDefaultHeaders = false;

                if (splitLine.Length >= 2)
                {
                    if (int.TryParse(splitLine[1], out _))
                    {
                        // Second column has a number; this is not a header line
                        useDefaultHeaders = true;
                    }
                    else
                    {
                        for (var index = 0; index < splitLine.Length; index++)
                        {
                            if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                            {
                                // Recognized column name; update columnMapping
                                columnMapping[resultFileColumn] = index;
                                useDefaultHeaders = false;
                            }
                            else
                            {
                                // Unrecognized column name
                                Console.WriteLine("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMODPlusResultsFileHeaderLine");
                            }
                        }
                    }
                }

                if (useDefaultHeaders)
                {
                    // Use default column mappings
                    foreach (MODPlusResultsFileColumns resultColumn in Enum.GetValues(typeof(MODPlusResultsFileColumns)))
                    {
                        columnMapping[resultColumn] = (int)resultColumn;
                    }

                    // This is not a header line; return false
                    return false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MODPlus results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        /// <summary>
        /// Parse the header line of a MOD+ _syn.txt file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODPlusSynFileHeaderLine(string lineIn, IDictionary<MODPlusSynFileColumns, int> columnMapping)
        {
            var columnNames = MODPlusSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MODPlusSynFileColumns resultColumn in Enum.GetValues(typeof(MODPlusSynFileColumns)))
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
                SetErrorMessage("Error parsing header in MODPlus synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMODPlusSynFileEntry(
            string lineIn,
            MODPlusResults searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<MODPlusSynFileColumns, int> columnMapping,
            out string peptideSequenceWithMods)
        {
            // Parses an entry from the MODPlus Synopsis file

            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();
            peptideSequenceWithMods = string.Empty;

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 13)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MODPlus Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from MODPlus Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.DelM], out string modPlusComputedDelM);
                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.DelM_PPM], out string modPlusComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.MODPlusComputedDelM = modPlusComputedDelM;
                searchResult.MODPlusComputedDelMPPM = modPlusComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.MODPlusComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MODPlus results file as "Observed_Mass - Theoretical_Mass"
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

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (SearchResultsBaseClass) searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Spectrum_Index], out string spectrumIndex);

                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.MH], out string parentIonMh);

                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Score], out string modPlusScore);
                GetColumnValue(splitLine, columnMapping[MODPlusSynFileColumns.Probability], out string probability);

                searchResult.Spectrum_Index = spectrumIndex;
                searchResult.Precursor_mz = precursorMz;
                searchResult.ParentIonMH = parentIonMh;
                searchResult.MODPlusScore = modPlusScore;
                searchResult.Probability = probability;

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorLog += "Error parsing MODPlus Results for RowIndex '" + splitLine[0] + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MODPlus Results in ParseMODPlusSynFileEntry\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MODPlus results file (Dataset_MODPlus.id.txt)</param>
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

                    // Load the MODPlus Parameter File to look for any static mods
                    var modInfoExtracted = ExtractModInfoFromMODPlusParamFile(SearchToolParameterFilePath, out var modPlusModInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Resolve the mods in modPlusModInfo with the ModDefs mods
                    ResolveMODPlusModsWithModDefinitions(modPlusModInfo);

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "modp.id" with "_modp"
                    if (baseName.EndsWith(FILENAME_SUFFIX_MODPlus_FILE, StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - FILENAME_SUFFIX_MODPlus_FILE.Length) + "_modp";
                    }

                    // Do not create a first-hits file for MODPlus results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_modp_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    success = ParseMODPlusSynopsisFile(synOutputFilePath, outputDirectoryPath, false);

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
                    SetErrorMessage("Error in MODPlusResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile (1):" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        private bool CreateProteinModsFileWork(
            string baseName,
            FileInfo inputFile,
            string synOutputFilePath,
            string outputDirectoryPath)
        {
            bool success;

            if (inputFile.Directory == null)
            {
                ReportWarning("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                return false;
            }

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
                    // We do this because some peptides reported by MODPlus may not match the fasta file (due to amino acid substitutions)
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
                if (string.IsNullOrWhiteSpace(synOutputFilePath))
                {
                    ReportWarning("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputDirectoryPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          PeptideHitResultTypes.MODPlus);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                return true;
            }

            return true;
        }

        private void ResolveMODPlusModsWithModDefinitions(IReadOnlyCollection<ModificationDefinition> modPlusModInfo)
        {
            if (modPlusModInfo == null)
            {
                return;
            }

            // Call .LookupModificationDefinitionByMass for each entry in modPlusModInfo
            foreach (var modInfo in modPlusModInfo)
            {
                if (string.IsNullOrEmpty(modInfo.TargetResidues))
                {
                    mPeptideMods.LookupModificationDefinitionByMassAndModType(
                        modInfo.ModificationMass, modInfo.ModificationType, default,
                        AminoAcidModInfo.ResidueTerminusState.None, out _, true);
                }
                else
                {
                    foreach (var chTargetResidue in modInfo.TargetResidues)
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(
                            modInfo.ModificationMass, modInfo.ModificationType, chTargetResidue,
                            AminoAcidModInfo.ResidueTerminusState.None, out _, true);
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            List<MODPlusSearchResult> filteredSearchResults,
            ref string errorLog)
        {
            // Sort filteredSearchResults by descending score, ascending scan, ascending charge, ascending peptide, and ascending protein
            filteredSearchResults.Sort(new MODPlusSearchResultsComparerScoreScanChargePeptide());

            // Compute FDR values then assign QValues
            ComputeQValues(filteredSearchResults);

            if (mProteinNamePositionSplit == null)
            {
                mProteinNamePositionSplit = new Regex(@"(.+)\[([^\]]+)\]", RegexOptions.Compiled);
            }

            var resultID = 1;
            foreach (var result in filteredSearchResults)
            {
                var proteinList = result.ProteinList.Split(';');

                if (proteinList.Length == 0)
                {
                    // This code should not be reached
                    WriteSearchResultToFile(resultID, writer, result, "Unknown_Protein", string.Empty, ref errorLog);
                    resultID++;
                }

                foreach (var proteinEntry in proteinList)
                {
                    string proteinName;
                    string peptidePosition;

                    var reMatch = mProteinNamePositionSplit.Match(proteinEntry);
                    if (reMatch.Success)
                    {
                        proteinName = reMatch.Groups[1].Value;
                        peptidePosition = reMatch.Groups[2].Value;
                    }
                    else
                    {
                        proteinName = string.Copy(proteinEntry);
                        peptidePosition = string.Empty;
                    }

                    WriteSearchResultToFile(resultID, writer, result, proteinName, peptidePosition, ref errorLog);
                    resultID++;
                }
            }
        }

        /// <summary>
        /// Compute FDR values then assign QValues
        /// </summary>
        /// <param name="searchResults"></param>
        /// <remarks>Assumes the data is sorted by descending score using MODPlusSearchResultsComparerScoreScanChargePeptide</remarks>
        private void ComputeQValues(IList<MODPlusSearchResult> searchResults)
        {
            var forwardPeptideCount = 0;
            var reversePeptideCount = 0;

            for (var index = 0; index < searchResults.Count;)
            {
                // Check for entries with multiple proteins listed
                var indexEnd = index;
                while (indexEnd + 1 < searchResults.Count)
                {
                    if (searchResults[index].ScanNum == searchResults[indexEnd + 1].ScanNum &&
                        searchResults[index].ChargeNum == searchResults[indexEnd + 1].ChargeNum &&
                        searchResults[index].Peptide == searchResults[indexEnd + 1].Peptide)
                    {
                        indexEnd++;
                    }
                    else
                    {
                        break;
                    }
                }

                var isReverse = true;

                // Look for non-reverse proteins
                for (var indexCheck = index; indexCheck <= indexEnd; indexCheck++)
                {
                    var proteinList = searchResults[indexCheck].ProteinList.Split(';');

                    foreach (var proteinEntry in proteinList)
                    {
                        if (!IsReversedProtein(proteinEntry))
                        {
                            isReverse = false;
                            break;
                        }
                    }
                }

                if (isReverse)
                {
                    reversePeptideCount++;
                }
                else
                {
                    forwardPeptideCount++;
                }

                double fDR = 1;

                if (forwardPeptideCount > 0)
                {
                    fDR = reversePeptideCount / Convert.ToDouble(forwardPeptideCount);
                }

                // Store the FDR values
                for (var indexStore = index; indexStore <= indexEnd; indexStore++)
                {
                    var udtResult = searchResults[indexStore];
                    udtResult.FDR = fDR;

                    searchResults[indexStore] = udtResult;
                }

                index = indexEnd + 1;
            }

            // Now compute QValues
            // We step through the list, from the worst scoring result to the best result
            // The first QValue is the FDR of the final entry
            // The next QValue is the minimum of (QValue, CurrentFDR)

            var qValue = searchResults.Last().FDR;
            if (qValue > 1)
                qValue = 1;

            for (var index = searchResults.Count - 1; index >= 0; index += -1)
            {
                var udtResult = searchResults[index];

                qValue = Math.Min(qValue, udtResult.FDR);
                udtResult.QValue = qValue;

                searchResults[index] = udtResult;
            }
        }

        private void StoreSynMatches(
            IList<MODPlusSearchResult> searchResults,
            int startIndex,
            int endIndex,
            ICollection<MODPlusSearchResult> filteredSearchResults)
        {
            AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].ProbabilityNum >= MODaMODPlusSynopsisFileProbabilityThreshold)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorLog"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ref string errorLog)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = MODPlusSynFileReader.GetColumnHeaderNamesAndIDs();

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
        /// <param name="proteinName"></param>
        /// <param name="peptidePosition"></param>
        /// <param name="errorLog"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            MODPlusSearchResult udtSearchResult,
            string proteinName,
            string peptidePosition,
            ref string errorLog)
        {
            try
            {
                // Primary Columns
                //
                // MODPlus
                // ResultID	Scan	Spectrum_Index	Charge	PrecursorMZ	DelM	DelM_PPM	MH	Peptide	NTT	ModificationAnnotation	Protein	Peptide_Position	Score	Probability	Rank_Probability   QValue

                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.ScanNum.ToString(),
                    udtSearchResult.SpectrumIndex,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.DelM,
                    udtSearchResult.DelM_PPM,
                    udtSearchResult.MH,
                    udtSearchResult.Peptide,
                    udtSearchResult.NTT,
                    udtSearchResult.ModificationAnnotation,
                    proteinName,
                    peptidePosition,
                    udtSearchResult.Score,
                    udtSearchResult.Probability,
                    udtSearchResult.RankScore.ToString(),
                    PRISM.StringUtilities.DblToString(udtSearchResult.QValue, 5, 0.00005)
                };

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

        private class MODPlusSearchResultsComparerScanChargeScorePeptide : IComparer<MODPlusSearchResult>
        {
            public int Compare(MODPlusSearchResult x, MODPlusSearchResult y)
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

                // Charge is the same; check ScoreNum
                if (x.ScoreNum < y.ScoreNum)
                {
                    return 1;
                }

                if (x.ScoreNum > y.ScoreNum)
                {
                    return -1;
                }

                // Probability is the same; check peptide
                var result = string.CompareOrdinal(x.Peptide, y.Peptide);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.CompareOrdinal(x.ProteinList, y.ProteinList);
                }
                return result;
            }
        }

        private class MODPlusSearchResultsComparerScoreScanChargePeptide : IComparer<MODPlusSearchResult>
        {
            public int Compare(MODPlusSearchResult x, MODPlusSearchResult y)
            {
                if (x.ScoreNum < y.ScoreNum)
                {
                    return 1;
                }

                if (x.ScoreNum > y.ScoreNum)
                {
                    return -1;
                }

                // P-value is the same; check scan number
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
                var result = string.CompareOrdinal(x.Peptide, y.Peptide);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.CompareOrdinal(x.ProteinList, y.ProteinList);
                }
                return result;
            }
        }
    }
}
