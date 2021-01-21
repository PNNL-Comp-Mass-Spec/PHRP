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
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class reads in a TopPIC results file (txt format) and creates
    /// a tab-delimited text file with the data.
    /// </summary>
    public class clsTopPICResultsProcessor : clsPHRPBaseClass
    {
        public clsTopPICResultsProcessor()
        {
            FileDate = "July 10, 2019";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string TOOL_NAME = "TopPIC";

        public const string FILENAME_SUFFIX_TopPIC_PROTEOFORMS_FILE = "_TopPIC_Proteoforms";

        public const string FILENAME_SUFFIX_TopPIC_PRSMs_FILE = "_TopPIC_PrSMs";

        public const string N_TERMINUS_SYMBOL_TopPIC = ".";

        public const string C_TERMINUS_SYMBOL_TopPIC = ".";

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        /// <summary>
        /// Regex to match mods in square brackets, examples:
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
        private enum eTopPICResultsFileColumns
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
            Protein_accession = 12,
            Protein_description = 13,
            First_residue = 14,
            Last_residue = 15,
            Proteoform = 16,
            Unexpected_modifications = 17,
            MIScore = 18,
            Variable_PTMs = 19,
            Matched_peaks = 20,
            Matched_fragment_ions = 21,
            Pvalue = 22,
            Evalue = 23,
            Qvalue = 24,                     //  Spectral FDR, or PepFDR
            Proteoform_FDR = 25
        }

        #endregion

        #region "Structures"

        // This data structure holds rows read from the tab-delimited file created directly by TopPIC
        private struct udtTopPICSearchResultType
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
            public string Qvalue;
            public string Proteoform_FDR;

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
                Qvalue = string.Empty;
                Proteoform_FDR = string.Empty;
            }
        }

        #endregion

        #region "Class wide Variables"

        private int mDeltaMassWarningCount;

        private readonly SortedSet<string> mUnknownNamedMods = new SortedSet<string>();

        #endregion

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
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
                    residueLocInPeptide += 1;

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
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                            if (modificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                searchResult.SearchResultAddModification(
                                    modificationDefinition, chChar, residueLocInPeptide,
                                    searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
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

                    if (!parsingModInfo) continue;

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

                    if (double.TryParse(modMassOrName, out var modMass))
                    {
                        var success = searchResult.SearchResultAddModification(modMass, residueForMod, residueLocForMod, searchResult.DetermineResidueTerminusState(residueLocForMod), updateModOccurrenceCounts);
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
                            var success = searchResult.SearchResultAddModification(modMassFromUniMod, residueForMod, residueLocForMod, searchResult.DetermineResidueTerminusState(residueLocForMod), updateModOccurrenceCounts);
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

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
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
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            IList<udtTopPICSearchResultType> searchResults,
            int startIndex,
            int endIndex)
        {
            // Duplicate a portion of searchResults so that we can sort by PValue

            var resultsSubset = new Dictionary<int, udtTopPICSearchResultType>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsByProbability = (from item in resultsSubset orderby item.Value.PValueNum select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByProbability)
            {
                var result = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = result.PValueNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(result.PValueNum - lastValue) > double.Epsilon)
                    {
                        lastValue = result.PValueNum;
                        currentRank += 1;
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
        /// <returns></returns>
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

        private static readonly Regex ModMatcher = new Regex(TopPIC_MOD_MASS_OR_NAME_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form [23.5432] or [Acetyl]</param>
        /// <returns></returns>
        /// <remarks></remarks>
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
        /// <returns></returns>
        /// <remarks></remarks>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath)
        {

            var columnMapping = new Dictionary<eTopPICResultsFileColumns, int>();

            try
            {
                var errorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var headerParsed = false;

                    mDeltaMassWarningCount = 0;

                    // Initialize array that will hold all of the records in the TopPIC result file
                    var searchResultsUnfiltered = new List<udtTopPICSearchResultType>();

                    // Initialize the array that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<udtTopPICSearchResultType>();

                    // Parse the input file
                    while (!reader.EndOfStream & !AbortProcessing)
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
                                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                return false;
                            }

                            headerParsed = true;

                            // Write the header line
                            WriteSynFHTFileHeader(writer, ref errorLog);

                            continue;
                        }

                        var udtSearchResult = new udtTopPICSearchResultType();
                        var validSearchResult = ParseTopPICResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping);

                        if (validSearchResult)
                        {
                            searchResultsUnfiltered.Add(udtSearchResult);
                        }

                        // Update the progress
                        var percentComplete = Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100);
                        if (CreateProteinModsFile)
                        {
                            percentComplete = percentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                        }

                        UpdateProgress(percentComplete);
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
                            endIndex += 1;
                        }

                        // Store the results for this scan
                        StoreSynMatches(searchResultsUnfiltered, startIndex, endIndex, filteredSearchResults);

                        startIndex = endIndex + 1;
                    }

                    // Sort the data in udtFilteredSearchResults then write out to disk
                    SortAndWriteFilteredSearchResults(writer, filteredSearchResults, ref errorLog);
                }

                // Inform the user if any errors occurred
                if (errorLog.Length > 0)
                {
                    SetErrorMessage("Invalid Lines: " + "\n" + errorLog);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
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
            out List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor(TOOL_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                topPICParamFilePath,
                clsMSGFPlusParamFileModExtractor.ModSpecFormats.TopPIC,
                out modInfo);

            if (!success || mErrorCode != ePHRPErrorCodes.NoError)
            {
                if (mErrorCode == ePHRPErrorCodes.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the TopPIC parameter file");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
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

            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out primarySequenceWithMods, out prefix, out suffix))
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

        private bool ParseTopPICSynopsisFile(string inputFilePath, string outputDirectoryPath, ref List<udtPepToProteinMappingType> pepToProteinMapping, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that TopPIC synopsis files are normally sorted on PValue value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
            // Although this was a possibility with Inspect, it likely never occurs for TopPIC
            //  But, we'll keep the check in place just in case

            var columnMapping = new Dictionary<clsPHRPParserTopPIC.TopPICSynFileColumns, int>();

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
                var searchResult = new clsSearchResultsTopPIC(mPeptideMods, mPeptideSeqMassCalculator);

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
                        while (!reader.EndOfStream & !AbortProcessing)
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
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseTopPICSynFileEntry(lineIn, searchResult, ref errorLog,
                                                                            resultsProcessed, columnMapping,
                                                                            out var currentPeptideWithMods);

                            resultsProcessed += 1;
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
                                                searchResult.ResultID + "'" + "\n";
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

                                        pepToProteinMapIndex += 1;
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
                            var percentComplete = Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100);
                            if (CreateProteinModsFile)
                            {
                                percentComplete = percentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                            }
                            UpdateProgress(percentComplete);
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
                        SetErrorMessage("Invalid Lines: " + "\n" + errorLog);
                    }

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
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
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
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
        /// <returns></returns>
        private bool ParseTopPICResultsFileEntry(
            string lineIn,
            ref udtTopPICSearchResultType udtSearchResult,
            ref string errorLog,
            IDictionary<eTopPICResultsFileColumns, int> columnMapping)
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

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);
                if (!GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Prsm_ID], out udtSearchResult.Prsm_ID))
                {
                    ReportError("Prsm_ID column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Spectrum_ID], out udtSearchResult.Spectrum_ID);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                if (!GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Scans], out udtSearchResult.Scans))
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

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Peaks], out udtSearchResult.Peaks);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Monoisotopic mass value of the observed precursor_mz
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Precursor_mass], out udtSearchResult.Precursor_mass);

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

                if (columnMapping[eTopPICResultsFileColumns.Adjusted_precursor_mass] >= 0)
                {
                    // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC
                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Adjusted_precursor_mass], out udtSearchResult.Adjusted_precursor_mass);

                    double.TryParse(udtSearchResult.Adjusted_precursor_mass, out peptideMonoMassTopPIC);
                }
                else
                {
                    peptideMonoMassTopPIC = 0;
                }

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Proteoform_ID], out udtSearchResult.Proteoform_ID);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Feature_intensity], out udtSearchResult.Feature_Intensity);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Protein_accession], out udtSearchResult.Protein);
                udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Protein_description], out udtSearchResult.ProteinDescription);

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.First_residue], out udtSearchResult.ResidueStart);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Last_residue], out udtSearchResult.ResidueEnd);

                if (!GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Proteoform], out udtSearchResult.Proteoform))
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
                            PRISM.StringUtilities.DblToString(clsPeptideMassCalculator.MassToPPM(delM, precursorMZ), 5, 0.00005);
                    }
                    else
                    {
                        udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(clsPeptideMassCalculator.MassToPPM(delM, 1000), 5, 0.00005);
                    }
                }

                // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoMassPHRP, 0), 6);

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Unexpected_modifications], out udtSearchResult.Unexpected_Mod_Count);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.MIScore], out udtSearchResult.MIScore);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Variable_PTMs], out udtSearchResult.VariablePTMs);

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Matched_peaks], out udtSearchResult.Matched_peaks);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Matched_fragment_ions], out udtSearchResult.Matched_fragment_ions);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Pvalue], out udtSearchResult.Pvalue);
                if (!double.TryParse(udtSearchResult.Pvalue, out udtSearchResult.PValueNum))
                    udtSearchResult.PValueNum = 0;

                // Assure that the following are truly integers (Matched_peaks and Matched_fragment_ions are often of the form 8.0)
                udtSearchResult.Unexpected_Mod_Count = AssureInteger(udtSearchResult.Unexpected_Mod_Count, 0);      // Unexpected_Mod_Count
                udtSearchResult.Peaks = AssureInteger(udtSearchResult.Peaks, 0);                                    // Peak_count
                udtSearchResult.Matched_peaks = AssureInteger(udtSearchResult.Matched_peaks, 0);                    // Matched_Peak_Count
                udtSearchResult.Matched_fragment_ions = AssureInteger(udtSearchResult.Matched_fragment_ions, 0);    // Matched_Fragment_Ion_Count

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Evalue], out udtSearchResult.Evalue);
                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Qvalue], out udtSearchResult.Qvalue);

                if (string.Equals(udtSearchResult.Qvalue, "infinity", StringComparison.OrdinalIgnoreCase))
                {
                    udtSearchResult.Qvalue = "10";
                }
                else if (!string.IsNullOrEmpty(udtSearchResult.Qvalue) & !double.TryParse(udtSearchResult.Qvalue, out _))
                {
                    udtSearchResult.Qvalue = string.Empty;
                }

                GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Proteoform_FDR], out udtSearchResult.Proteoform_FDR);

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the TopPIC results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing TopPIC Results in ParseTopPICResultsFileEntry for RowIndex '" + splitLine[0] + "'" +
                                       "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing TopPIC Results in ParseTopPICResultsFileEntry" + "\n";
                    }
                }
                return false;
            }

        }

        /// <summary>
        /// Parse the header line of a TopPIC results file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns></returns>
        private bool ParseTopPICResultsFileHeaderLine(string lineIn, IDictionary<eTopPICResultsFileColumns, int> columnMapping)
        {
            // Header prior to November 2018:
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Protein name    First residue    Last residue    Proteoform    #unexpected modifications    #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR    #Variable PTMs

            // Header for TopPIC 1.2 and newer
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    Retention time    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Protein accession    Protein description    First residue    Last residue    Proteoform    #unexpected modifications    MIScore    #variable PTMs    #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR

            var columnNames = new SortedDictionary<string, eTopPICResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Data file name", eTopPICResultsFileColumns.SpectrumFileName},
                {"Prsm ID", eTopPICResultsFileColumns.Prsm_ID},
                {"Spectrum ID", eTopPICResultsFileColumns.Spectrum_ID},
                {"Fragmentation", eTopPICResultsFileColumns.FragMethod},
                {"Scan(s)", eTopPICResultsFileColumns.Scans},
                {"Retention time", eTopPICResultsFileColumns.RetentionTime},
                {"#peaks", eTopPICResultsFileColumns.Peaks},
                {"Charge", eTopPICResultsFileColumns.Charge},
                {"Precursor mass", eTopPICResultsFileColumns.Precursor_mass},
                {"Adjusted precursor mass", eTopPICResultsFileColumns.Adjusted_precursor_mass},
                {"Proteoform ID", eTopPICResultsFileColumns.Proteoform_ID},
                {"Feature intensity", eTopPICResultsFileColumns.Feature_intensity},
                {"Protein name", eTopPICResultsFileColumns.Protein_accession},
                {"Protein accession", eTopPICResultsFileColumns.Protein_accession},
                {"Protein description", eTopPICResultsFileColumns.Protein_description},
                {"First residue", eTopPICResultsFileColumns.First_residue},
                {"Last residue", eTopPICResultsFileColumns.Last_residue},
                {"Proteoform", eTopPICResultsFileColumns.Proteoform},
                {"#unexpected modifications", eTopPICResultsFileColumns.Unexpected_modifications},
                {"MIScore", eTopPICResultsFileColumns.MIScore},
                {"#variable PTMs", eTopPICResultsFileColumns.Variable_PTMs},
                {"#matched peaks", eTopPICResultsFileColumns.Matched_peaks},
                {"#matched fragment ions", eTopPICResultsFileColumns.Matched_fragment_ions},
                {"P-value", eTopPICResultsFileColumns.Pvalue},
                {"E-value", eTopPICResultsFileColumns.Evalue},
                {"Q-value (spectral FDR)", eTopPICResultsFileColumns.Qvalue},
                {"Proteoform FDR", eTopPICResultsFileColumns.Proteoform_FDR}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (eTopPICResultsFileColumns resultColumn in Enum.GetValues(typeof(eTopPICResultsFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');
                for (var index = 0; index <= splitLine.Length - 1; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var eResultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[eResultFileColumn] = index;
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
        /// Parse the header line of a TopPIC _syn.txt file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns></returns>
        private bool ParseTopPICSynFileHeaderLine(string lineIn, IDictionary<clsPHRPParserTopPIC.TopPICSynFileColumns, int> columnMapping)
        {
            var columnNames = clsPHRPParserTopPIC.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (clsPHRPParserTopPIC.TopPICSynFileColumns resultColumn in Enum.GetValues(typeof(clsPHRPParserTopPIC.TopPICSynFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');
                for (var index = 0; index <= splitLine.Length - 1; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var eResultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[eResultFileColumn] = index;
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
        /// <returns></returns>
        private bool ParseTopPICSynFileEntry(
            string lineIn,
            clsSearchResultsTopPIC searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<clsPHRPParserTopPIC.TopPICSynFileColumns, int> columnMapping,
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

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from TopPIC Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from TopPIC Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.DelM], out string msAlignComputedDelM);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.DelMPPM], out string msAlignComputedDelMppm);

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

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Prsm_ID], out string prsmId);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Spectrum_ID], out string spectrumId);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.MH], out string parentIonMH);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Unexpected_Mod_Count], out string unexpectedModCount);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Peak_Count], out string peakCount);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Matched_Peak_Count], out string matchedPeakCount);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Matched_Fragment_Ion_Count], out string matchedFragmentIonCount);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.PValue], out string pValue);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Rank_PValue], out string rankPValue);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.EValue], out string eValue);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.QValue], out string qValue);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.FragMethod], out string fragMethod);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.Proteoform_FDR], out string proteoformFDR);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserTopPIC.TopPICSynFileColumns.VariablePTMs], out string variablePTMs);

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
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing TopPIC Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing TopPIC Results in ParseTopPICSynFileEntry" + "\n";
                    }
                }
                return false;
            }

        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">TopPIC results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            var success = false;

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
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

                    var pepToProteinMapping = new List<udtPepToProteinMappingType>();

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
                        if (!baseName.EndsWith(suffix, StringComparison.OrdinalIgnoreCase)) continue;
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
                    SetErrorMessage("Error in clsTopPICResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in clsTopPICResultsProcessor.ProcessFile (1):" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return success;
        }

        private bool CreateProteinModsFileWork(string baseName, FileInfo inputFile, string synOutputFilePath, string outputDirectoryPath)
        {
            bool success;

            // Create the MTSPepToProteinMap file

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseName, outputDirectoryPath, mts: true);

            var sourcePHRPDataFiles = new List<string>();

            if (!string.IsNullOrEmpty(synOutputFilePath))
            {
                sourcePHRPDataFiles.Add(synOutputFilePath);
            }

            if (sourcePHRPDataFiles.Count == 0)
            {
                SetErrorMessage("Cannot call CreatePepToProteinMapFile since sourcePHRPDataFiles is empty");
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
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
                                                          clsPHRPReader.ePeptideHitResultType.TopPIC);
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
                peptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_TopPIC.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_TopPIC))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_TopPIC.Length) + "." + clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<udtTopPICSearchResultType> filteredSearchResults,
            ref string errorLog)
        {
            // Sort filteredSearchResults by ascending PValue, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.PValueNum, item.ScanNum, item.ChargeNum, item.Proteoform, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, ref errorLog);
                index += 1;
            }
        }

        private void StoreSynMatches(
            IList<udtTopPICSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            ICollection<udtTopPICSearchResultType> filteredSearchResults)
        {
            AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by scan, charge, and PValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].PValueNum <= MSAlignAndTopPICSynopsisFilePValueThreshold)
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
                var headerColumns = clsPHRPParserTopPIC.GetColumnHeaderNamesAndIDs();

                var headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();

                writer.WriteLine(CollapseList(headerNames));
            }
            catch (Exception)
            {
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error writing synopsis / first hits header" + "\n";
                }
            }
        }

        /// <summary>
        /// Writes an entry to the synopsis file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtTopPICSearchResultType udtSearchResult,
            ref string errorLog)
        {
            try
            {
                // Primary Columns
                //
                // TopPIC
                // ResultID  Scan  Prsm_ID  Spectrum_ID  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  Protein_Mass  Unexpected_Mod_Count  Peak_Count  Matched_Peak_Count  Matched_Fragment_Ion_Count  PValue  Rank_PValue  EValue  QValue  ProteoformFDR  FragMethod  VariablePTMs

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
                    udtSearchResult.Proteoform,
                    udtSearchResult.Proteoform_ID,
                    udtSearchResult.Feature_Intensity,
                    udtSearchResult.Protein,
                    udtSearchResult.ProteinDescription,
                    udtSearchResult.ResidueStart,
                    udtSearchResult.ResidueEnd,
                    udtSearchResult.Unexpected_Mod_Count,       // Unexpected_Mod_Count
                    udtSearchResult.Peaks,                      // Peak_count
                    udtSearchResult.Matched_peaks,              // Matched_Peak_Count
                    udtSearchResult.Matched_fragment_ions,      // Matched_Fragment_Ion_Count
                    udtSearchResult.Pvalue,
                    udtSearchResult.RankPValue.ToString(),
                    udtSearchResult.Evalue,
                    udtSearchResult.Qvalue,
                    udtSearchResult.Proteoform_FDR,
                    udtSearchResult.VariablePTMs
                };

                writer.WriteLine(CollapseList(data));
            }
            catch (Exception)
            {
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error writing synopsis / first hits record" + "\n";
                }
            }
        }

        #region "Event Handlers"

        private void ModExtractorErrorHandler(string message, Exception ex)
        {
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
        }

        #endregion

        #region "IComparer Classes"

        private class TopPICSearchResultsComparerScanChargePValuePeptide : IComparer<udtTopPICSearchResultType>
        {
            public int Compare(udtTopPICSearchResultType x, udtTopPICSearchResultType y)
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
                var result = string.Compare(x.Pvalue, y.Pvalue, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Pvalue is the same; check peptide
                    result = string.Compare(x.Proteoform, y.Proteoform, StringComparison.Ordinal);
                    if (result == 0)
                    {
                        // Peptide is the same, check Protein
                        result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                    }
                }
                return result;
            }
        }

        #endregion
    }
}
