﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class reads in an MSAlign results file (txt format) and creates
    /// a tab-delimited text file with the data.
    /// </summary>
    /// <remarks>
    /// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
    /// Started 11/28/2012
    /// </remarks>
    public class clsMSAlignResultsProcessor : clsPHRPBaseClass
    {
        // Ignore Spelling: monoisotopic, enums, txt, histone, methylation, fht, Da, Prsm, IodoAcet, IodoAcid, Frag, pvalue, Defs

        public clsMSAlignResultsProcessor()
        {
            FileDate = "January 27, 2021";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string TOOL_NAME = "MSAlign";

        public const string FILENAME_SUFFIX_MSALIGN_FILE = "_MSAlign_ResultTable";

        public const string N_TERMINUS_SYMBOL_MSALIGN = ".";

        public const string C_TERMINUS_SYMBOL_MSALIGN = ".";

        public const float DEFAULT_SYN_FILE_PVALUE_THRESHOLD = 0.95f;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string MSALIGN_MOD_MASS_REGEX = @"\[([+-]*[0-9\.]+)\]";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to the tab-delimited file created directly by MSAlign
        /// </summary>
        private enum MSAlignResultsFileColumns
        {
            SpectrumFileName = 0,
            Prsm_ID = 1,
            Spectrum_ID = 2,
            Protein_Sequence_ID = 3,         // Used by MSAlign v0.5, but not by v0.6
            Scans = 4,
            Peaks = 5,
            Charge = 6,
            Precursor_mass = 7,              // Monoisotopic mass value of the observed precursor_mz
            Adjusted_precursor_mass = 8,     // Theoretical monoisotopic mass of the peptide (including mods)
            Protein_ID = 9,
            Protein_name = 10,               // Protein name and description
            Protein_mass = 11,
            First_residue = 12,
            Last_residue = 13,
            Peptide = 14,
            Unexpected_modifications = 15,
            Matched_peaks = 16,
            Matched_fragment_ions = 17,
            Pvalue = 18,
            Evalue = 19,
            FDR = 20,
            Species_ID = 21,             // Present between Protein_ID and Protein_name in MSAlign_Histone result files
            FragMethod = 22              // Present as the last column in MSAlign_Histone result files
        }

        #endregion

        #region "Structures"

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file created directly by MSAlign
        /// </summary>
        private struct udtMSAlignSearchResultType
        {
            public string SpectrumFileName;
            public string Scans;
            public int ScanNum;
            public string Prsm_ID;
            public string Spectrum_ID;
            public string Peaks;
            public string Charge;
            public short ChargeNum;
            public string Precursor_mass;               // Monoisotopic mass value of the observed precursor_mz
            public string PrecursorMZ;                  // Computed from Precursor_mass
            public string Adjusted_precursor_mass;      // Theoretical monoisotopic mass of the peptide (including mods), as computed by MSAlign
            public string MH;                           // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                         // Computed using Precursor_mass - Adjusted_precursor_mass
            public string DelM_PPM;                     // Computed using DelM and Adjusted_precursor_mass
            public string Protein_ID;
            public string Species_ID;                   // Only present in MSAlign_Histone results
            public string Protein;
            public string Protein_mass;
            public string First_residue;
            public string Last_residue;
            public string Peptide;
            public string Unexpected_modifications;
            public string Matched_peaks;
            public string Matched_fragment_ions;
            public string Pvalue;
            public double PValueNum;
            public string Evalue;
            public string FDR;
            public string FragMethod;                   // Only present in MSAlign_Histone results
            public int RankPValue;

            public void Clear()
            {
                SpectrumFileName = string.Empty;
                Prsm_ID = string.Empty;
                Spectrum_ID = string.Empty;
                Scans = string.Empty;
                ScanNum = 0;
                Peaks = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                Precursor_mass = string.Empty;
                Adjusted_precursor_mass = string.Empty;
                MH = string.Empty;
                DelM = string.Empty;
                DelM_PPM = string.Empty;
                Protein_ID = string.Empty;
                Species_ID = string.Empty;
                Protein = string.Empty;
                Protein_mass = string.Empty;
                First_residue = string.Empty;
                Last_residue = string.Empty;
                Peptide = string.Empty;
                Unexpected_modifications = string.Empty;
                Matched_peaks = string.Empty;
                Matched_fragment_ions = string.Empty;
                Pvalue = string.Empty;
                PValueNum = 0;
                Evalue = string.Empty;
                FDR = string.Empty;
                FragMethod = string.Empty;
                RankPValue = 0;
            }
        }

        #endregion

        #region "Class wide Variables"

        private int mDeltaMassWarningCount;

        #endregion

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const char NO_RESIDUE = '-';

            // Preferably use three digits of precision, but allow two digits of precision
            // since the MSAlign_ResultTable.txt file lists modifications with just two digits after the decimal, e.g. (DIQM)[16.00]
            const int MSALIGN_MASS_DIGITS_OF_PRECISION = 3;

            const int MSALIGN_MASS_DIGITS_OF_PRECISION_LOOSE = 2;

            var parsingModMass = false;
            var modMassDigits = string.Empty;

            var chMostRecentResidue = NO_RESIDUE;
            var residueLocInPeptide = 0;

            var chAmbiguousResidue = NO_RESIDUE;
            var ambiguousResidueLocInPeptide = 0;

            var clearAmbiguousResidue = false;
            var storeAmbiguousResidue = false;

            // ReSharper disable CommentTypo

            // Examine sequence, which will look something like:
            // (ST)[-1.02]YSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
            // (VHTFPAVLQSSGLYSLSSV)[-.99]VTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTH

            // ReSharper restore CommentTypo

            // When a modification mass follows a series of residues in parentheses, that means an ambiguous mod,
            // i.e. we don't know which residue to associate the mod with

            var sequence = searchResult.PeptideSequenceWithMods;

            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var chChar = sequence[index];

                if (IsLetterAtoZ(chChar))
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
                        // Only add this modification if it is a static mod; dynamic mods are handled later in this method when a ']' is found
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ModificationTypeConstants.StaticMod)
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
                    // Mod Mass Start
                    modMassDigits = string.Empty;
                    parsingModMass = true;
                }
                else if (chChar == ']')
                {
                    // Mod Mass End

                    if (!parsingModMass)
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
                            MSALIGN_MASS_DIGITS_OF_PRECISION, MSALIGN_MASS_DIGITS_OF_PRECISION_LOOSE);

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

                    parsingModMass = false;
                }
                else if (parsingModMass)
                {
                    modMassDigits += chChar;
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

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        private void AssignRankAndDeltaNormValues(
            IList<udtMSAlignSearchResultType> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of searchResults so that we can sort by PValue

            var resultsSubset = new Dictionary<int, udtMSAlignSearchResultType>();
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

        private static readonly Regex RegexModMassRegEx = new Regex(MSALIGN_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form [23.5432]</param>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            foreach (Match reMatch in RegexModMassRegEx.Matches(peptide))
            {
                if (double.TryParse(reMatch.Groups[1].Value, out var modMassFound))
                {
                    totalModMass += modMassFound;
                }
            }

            return totalModMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var suffixesToFind = new List<string> {
                "_msalign_syn",
                "_msalign_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MSAlign
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath)
        {
            var columnMapping = new Dictionary<MSAlignResultsFileColumns, int>();

            try
            {
                var errorLog = string.Empty;
                var includeSpeciesAndFragMethod = false;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var headerParsed = false;

                    mDeltaMassWarningCount = 0;

                    // Initialize array that will hold all of the records in the MSAlign result file
                    var searchResultsUnfiltered = new List<udtMSAlignSearchResultType>();

                    // Initialize the array that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<udtMSAlignSearchResultType>();

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
                            var validHeader = ParseMSAlignResultsFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
                                return false;
                            }

                            headerParsed = true;

                            // Write the header line
                            WriteSynFHTFileHeader(writer, columnMapping, out includeSpeciesAndFragMethod, ref errorLog);

                            continue;
                        }

                        var udtSearchResult = new udtMSAlignSearchResultType();
                        var validSearchResult = ParseMSAlignResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping);

                        if (validSearchResult)
                        {
                            searchResultsUnfiltered.Add(udtSearchResult);
                        }

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    // Sort the SearchResults by scan, charge, and ascending PValue
                    searchResultsUnfiltered.Sort(new MSAlignSearchResultsComparerScanChargePValuePeptide());

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
                    SortAndWriteFilteredSearchResults(writer, filteredSearchResults, includeSpeciesAndFragMethod, ref errorLog);
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
                SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the MSAlign parameter file
        /// </summary>
        /// <param name="msAlignParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromMSAlignParamFile(string msAlignParamFilePath, out List<ModificationDefinition> modInfo)
        {
            var success = false;

            modInfo = new List<ModificationDefinition>();

            try
            {
                if (string.IsNullOrEmpty(msAlignParamFilePath))
                {
                    SetErrorMessage("MSAlign Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(PHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(msAlignParamFilePath))
                {
                    SetErrorMessage("MSAlign param file not found: " + msAlignParamFilePath);
                }
                else
                {
                    // Read the contents of the parameter (or mods) file
                    using (var reader = new StreamReader(new FileStream(msAlignParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!reader.EndOfStream)
                        {
                            var lineIn = reader.ReadLine();
                            if (string.IsNullOrWhiteSpace(lineIn))
                                continue;

                            var dataLine = lineIn.Trim();
                            if (dataLine.Length == 0)
                                continue;

                            if (dataLine.StartsWith("#"))
                            {
                                // Comment line; skip it
                            }
                            else
                            {
                                // Split the line on the equals sign
                                var kvSetting = SynFileReaderBaseClass.ParseKeyValueSetting(dataLine, '=', "#");

                                if (string.Equals(kvSetting.Key, "cysteineProtection", StringComparison.OrdinalIgnoreCase))
                                {
                                    ModificationDefinition modDef;
                                    switch (kvSetting.Value.ToUpper())
                                    {
                                        case "C57":
                                            modDef = new ModificationDefinition(ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                                                                   57.0215, "C",
                                                                                   ModificationDefinition.ModificationTypeConstants.StaticMod,
                                                                                   "IodoAcet");
                                            modInfo.Add(modDef);
                                            break;

                                        case "C58":
                                            modDef = new ModificationDefinition(ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                                                                   58.0055, "C",
                                                                                   ModificationDefinition.ModificationTypeConstants.StaticMod,
                                                                                   "IodoAcid");
                                            modInfo.Add(modDef);
                                            break;
                                    }
                                }
                            }
                        }
                    }

                    Console.WriteLine();

                    success = true;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MSAlign parameter file (" + Path.GetFileName(msAlignParamFilePath) + "): " + ex.Message);
                SetErrorCode(PHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                success = false;
            }

            return success;
        }

        private void InitializeLocalVariables()
        {
            // Nothing to do at present
        }

        private bool ParseMSAlignSynopsisFile(string inputFilePath, string outputDirectoryPath, ref List<udtPepToProteinMappingType> pepToProteinMapping, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that MSAlign synopsis files are normally sorted on PValue, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
            // Although this was a possibility with Inspect, it likely never occurs for MSAlign
            //  But, we'll keep the check in place just in case

            var columnMapping = new Dictionary<MSAlignSynFileReader.MSAlignSynFileColumns, int>();

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
                var searchResult = new clsSearchResultsMSAlign(mPeptideMods, mPeptideSeqMassCalculator);

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
                                var validHeader = ParseMSAlignSynFileHeaderLine(lineIn, columnMapping);
                                if (!validHeader)
                                {
                                    // Error parsing header
                                    SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMSAlignSynFileEntry(lineIn, searchResult, ref errorLog,
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
                                    errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'" +
                                                "\n";
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
                    SetErrorCode(PHRPErrorCodes.ErrorReadingInputFile);
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
                SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }
        }

        private bool ParseMSAlignResultsFileEntry(
            string lineIn,
            ref udtMSAlignSearchResultType udtSearchResult,
            ref string errorLog,
            IDictionary<MSAlignResultsFileColumns, int> columnMapping)
        {
            // Parses an entry from the MSAlign results file

            string[] splitLine = null;

            double precursorMZ = 0;

            try
            {
                udtSearchResult.Clear();
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 13)
                {
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);
                if (!GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Prsm_ID], out udtSearchResult.Prsm_ID))
                {
                    ReportError("Prsm_ID column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Spectrum_ID], out udtSearchResult.Spectrum_ID);

                if (!GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Scans], out udtSearchResult.Scans))
                {
                    ReportError("Scan(s) column is missing or invalid", true);
                }

                if (!int.TryParse(udtSearchResult.Scans, out udtSearchResult.ScanNum))
                {
                    // .Scans likely has a list of scan numbers; extract the first scan number from .scans
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

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Peaks], out udtSearchResult.Peaks);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Monoisotopic mass value of the observed precursor_mz
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Precursor_mass], out udtSearchResult.Precursor_mass);

                // precursorMonoMass is Observed m/z, converted to monoisotopic mass
                if (double.TryParse(udtSearchResult.Precursor_mass, out var precursorMonoMass))
                {
                    if (udtSearchResult.ChargeNum > 0)
                    {
                        precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, udtSearchResult.ChargeNum);
                        udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                    }
                }

                // peptideMonoMassMSAlign is Theoretical peptide monoisotopic mass, including mods, as computed by MSAlign
                double peptideMonoMassMSAlign;

                if (columnMapping[MSAlignResultsFileColumns.Adjusted_precursor_mass] >= 0)
                {
                    // Theoretical monoisotopic mass of the peptide (including mods), as computed by MSAlign
                    GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Adjusted_precursor_mass], out udtSearchResult.Adjusted_precursor_mass);

                    double.TryParse(udtSearchResult.Adjusted_precursor_mass, out peptideMonoMassMSAlign);
                }
                else
                {
                    peptideMonoMassMSAlign = 0;
                }

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Protein_ID], out udtSearchResult.Protein_ID);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Species_ID], out udtSearchResult.Species_ID);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Protein_name], out udtSearchResult.Protein);
                udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Protein_mass], out udtSearchResult.Protein_mass);

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.First_residue], out udtSearchResult.First_residue);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Last_residue], out udtSearchResult.Last_residue);

                if (!GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                // Add the standard terminus symbols to the peptide sequence
                udtSearchResult.Peptide = ReplaceTerminus(udtSearchResult.Peptide);

                // Parse the sequence to determine the total mod mass
                // Note that we do not remove any of the mod symbols since MSAlign identifies mods by mass alone, and since mods can ambiguously apply to residues
                var totalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                // Compute theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                var peptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, totalModMass);

                if (Math.Abs(peptideMonoMassMSAlign) < double.Epsilon)
                {
                    peptideMonoMassMSAlign = peptideMonoMassPHRP;
                }

                // Warn the user if the monoisotopic mass values differ by more than 0.1 Da
                ValidateMatchingMonoisotopicMass(TOOL_NAME, udtSearchResult.Peptide, peptideMonoMassPHRP, peptideMonoMassMSAlign, ref mDeltaMassWarningCount);

                if (peptideMonoMassMSAlign > 0)
                {
                    // Compute DelM and DelM_PPM
                    var delM = precursorMonoMass - peptideMonoMassMSAlign;
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

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Unexpected_modifications], out udtSearchResult.Unexpected_modifications);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Matched_peaks], out udtSearchResult.Matched_peaks);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Matched_fragment_ions], out udtSearchResult.Matched_fragment_ions);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Pvalue], out udtSearchResult.Pvalue);
                if (!double.TryParse(udtSearchResult.Pvalue, out udtSearchResult.PValueNum))
                    udtSearchResult.PValueNum = 0;

                // Assure that the following are truly integers (Matched_peaks and Matched_fragment_ions are often of the form 8.0)
                udtSearchResult.Unexpected_modifications = AssureInteger(udtSearchResult.Unexpected_modifications, 0); // Unexpected_Mod_Count
                udtSearchResult.Peaks = AssureInteger(udtSearchResult.Peaks, 0); // Peak_count
                udtSearchResult.Matched_peaks = AssureInteger(udtSearchResult.Matched_peaks, 0); // Matched_Peak_Count
                udtSearchResult.Matched_fragment_ions = AssureInteger(udtSearchResult.Matched_fragment_ions, 0); // Matched_Fragment_Ion_Count

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.Evalue], out udtSearchResult.Evalue);
                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.FDR], out udtSearchResult.FDR);

                if (string.Equals(udtSearchResult.FDR, "infinity", StringComparison.OrdinalIgnoreCase))
                {
                    udtSearchResult.FDR = "10";
                }
                else if (!string.IsNullOrEmpty(udtSearchResult.FDR) && !double.TryParse(udtSearchResult.FDR, out _))
                {
                    udtSearchResult.FDR = string.Empty;
                }

                GetColumnValue(splitLine, columnMapping[MSAlignResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the MSAlign results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorLog += "Error parsing MSAlign Results in ParseMSAlignResultsFileEntry for RowIndex '" + splitLine[0] + "'" +
                                       "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MSAlign Results in ParseMSAlignResultsFileEntry\n";
                    }
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a MSAlign results file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSAlignResultsFileHeaderLine(string lineIn, IDictionary<MSAlignResultsFileColumns, int> columnMapping)
        {
            // The expected header from MSAlign v0.5 is:
            //                   Prsm_ID    Spectrum_ID    Protein_Sequence_ID    Spectrum_ID    Scan(s)    #peaks    Charge    Precursor_mass                                                           Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions               E-value

            // The expected header from MSAlign v0.6 is:
            // Data_file_name    Prsm_ID    Spectrum_ID                                          Scan(s)    #peaks    Charge    Precursor_mass    Adjusted_precursor_mass    Protein_ID                  Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions    P-value    E-value    FDR

            // The expected header from MSAlign_Histone v0.9 is
            // Data_file_name    Prsm_ID    Spectrum_ID                                          Scan(s)    #peaks    Charge    Precursor_mass    Adjusted_precursor_mass    Protein_ID    Species_ID    Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions    P-value    E-value    FDR    FragMethod

            var columnNames = new SortedDictionary<string, MSAlignResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Data_file_name", MSAlignResultsFileColumns.SpectrumFileName},
                {"Prsm_ID", MSAlignResultsFileColumns.Prsm_ID},
                {"Spectrum_ID", MSAlignResultsFileColumns.Spectrum_ID},
                {"Protein_Sequence_ID", MSAlignResultsFileColumns.Protein_Sequence_ID},
                {"Scan(s)", MSAlignResultsFileColumns.Scans},
                {"#peaks", MSAlignResultsFileColumns.Peaks},
                {"Charge", MSAlignResultsFileColumns.Charge},
                {"Precursor_mass", MSAlignResultsFileColumns.Precursor_mass},
                {"Adjusted_precursor_mass", MSAlignResultsFileColumns.Adjusted_precursor_mass},
                {"Protein_ID", MSAlignResultsFileColumns.Protein_ID},
                {"Species_ID", MSAlignResultsFileColumns.Species_ID},
                {"Protein_name", MSAlignResultsFileColumns.Protein_name},
                {"Protein_mass", MSAlignResultsFileColumns.Protein_mass},
                {"First_residue", MSAlignResultsFileColumns.First_residue},
                {"Last_residue", MSAlignResultsFileColumns.Last_residue},
                {"Peptide", MSAlignResultsFileColumns.Peptide},
                {"#unexpected_modifications", MSAlignResultsFileColumns.Unexpected_modifications},
                {"#matched_peaks", MSAlignResultsFileColumns.Matched_peaks},
                {"#matched_fragment_ions", MSAlignResultsFileColumns.Matched_fragment_ions},
                {"P-value", MSAlignResultsFileColumns.Pvalue},
                {"E-value", MSAlignResultsFileColumns.Evalue},
                {"FDR", MSAlignResultsFileColumns.FDR},
                {"FragMethod", MSAlignResultsFileColumns.FragMethod}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MSAlignResultsFileColumns resultColumn in Enum.GetValues(typeof(MSAlignResultsFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');
                for (var index = 0; index <= splitLine.Length - 1; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultColumn] = index;
                    }
                    else
                    {
                        // Unrecognized column name
                        Console.WriteLine(
                            "Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMSAlignResultsFileHeaderLine");
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MSAlign results file: " + ex.Message);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse the header line of a MSAlign _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSAlignSynFileHeaderLine(string lineIn, IDictionary<MSAlignSynFileReader.MSAlignSynFileColumns, int> columnMapping)
        {
            var columnNames = MSAlignSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MSAlignSynFileReader.MSAlignSynFileColumns resultColumn in Enum.GetValues(typeof(MSAlignSynFileReader.MSAlignSynFileColumns)))
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
                SetErrorMessage("Error parsing header in MSAlign synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMSAlignSynFileEntry(
            string lineIn,
            clsSearchResultsMSAlign searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<MSAlignSynFileReader.MSAlignSynFileColumns, int> columnMapping,
            out string peptideSequenceWithMods)
        {
            // Parses an entry from the MSAlign Synopsis file

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

                if (!GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MSAlign Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from MSAlign Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.DelM], out string msAlignComputedDelM);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.DelMPPM], out string msAlignComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.MSAlignComputedDelM = msAlignComputedDelM;
                searchResult.MSAlignComputedDelMPPM = msAlignComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.MSAlignComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MSAlign results file as "Observed_Mass - Theoretical_Mass"
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

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Prsm_ID], out string prsmId);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Spectrum_ID], out string spectrumId);

                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.MH], out string parentIonMH);

                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Protein_Mass], out string proteinMass);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Unexpected_Mod_Count], out string unexpectedModCount);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Peak_Count], out string peakCount);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Matched_Peak_Count], out string matchedPeakCount);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Matched_Fragment_Ion_Count], out string matchedFragmentIonCount);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.PValue], out string pValue);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Rank_PValue], out string rankPValue);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.EValue], out string eValue);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.FDR], out string fdr);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.Species_ID], out string speciesId);
                GetColumnValue(splitLine, columnMapping[MSAlignSynFileReader.MSAlignSynFileColumns.FragMethod], out string fragMethod);

                searchResult.Prsm_ID = prsmId;
                searchResult.Spectrum_ID = spectrumId;
                searchResult.Precursor_mz = precursorMz;
                searchResult.ParentIonMH = parentIonMH;
                searchResult.Protein_Mass = proteinMass;
                searchResult.Unexpected_Mod_Count = unexpectedModCount;
                searchResult.Peak_Count = peakCount;
                searchResult.Matched_Peak_Count = matchedPeakCount;
                searchResult.Matched_Fragment_Ion_Count = matchedFragmentIonCount;
                searchResult.PValue = pValue;
                searchResult.Rank_PValue = rankPValue;
                searchResult.EValue = eValue;
                searchResult.FDR = fdr;
                searchResult.Species_ID = speciesId;
                searchResult.FragMethod = fragMethod;

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorLog += "Error parsing MSAlign Results for RowIndex '" + splitLine[0] + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MSAlign Results in ParseMSAlignSynFileEntry\n";
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MSAlign results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if successful, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            var success = false;

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                SetErrorCode(PHRPErrorCodes.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(PHRPErrorCodes.InvalidInputFilePath);
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

                    // Load the MSAlign Parameter File so that we can determine whether Cysteine residues are statically modified
                    var modInfoExtracted = ExtractModInfoFromMSAlignParamFile(SearchToolParameterFilePath, out var msAlignModInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Resolve the mods in msAlignModInfo with the ModDefs mods
                    ResolveMSAlignModsWithModDefinitions(msAlignModInfo);

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_MSAlign_ResultTable" with "_msalign"
                    if (baseName.EndsWith(FILENAME_SUFFIX_MSALIGN_FILE, StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - FILENAME_SUFFIX_MSALIGN_FILE.Length) + "_msalign";
                    }

                    // Do not create a first-hits file for MSAlign results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_msalign_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    success = ParseMSAlignSynopsisFile(synOutputFilePath, outputDirectoryPath, ref pepToProteinMapping, false);

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
                    SetErrorMessage("Error in clsMSAlignResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
                    SetErrorCode(PHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile (1):" + ex.Message, ex);
                SetErrorCode(PHRPErrorCodes.UnspecifiedError);
            }

            return success;
        }

        /// <summary>
        /// Create the MTSPepToProteinMap file
        /// </summary>
        /// <param name="baseName"></param>
        /// <param name="inputFile"></param>
        /// <param name="synOutputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateProteinModsFileWork(string baseName, FileInfo inputFile, string synOutputFilePath, string outputDirectoryPath)
        {
            bool success;

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
                SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
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
                    // We do this since a small number of peptides reported by MSAlign don't perfectly match the fasta file
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
                                                          PHRPReader.PHRPReader.PeptideHitResultTypes.MSAlign);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                success = true;
            }

            // ReSharper disable once ConditionIsAlwaysTrueOrFalse
            return success;
        }

        private string ReplaceTerminus(string peptide)
        {
            if (peptide.StartsWith(N_TERMINUS_SYMBOL_MSALIGN))
            {
                peptide = PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_MSALIGN.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_MSALIGN))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_MSALIGN.Length) + "." + PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        /// <summary>
        /// Call .LookupModificationDefinitionByMass for each entry in msAlignModInfo
        /// </summary>
        /// <param name="msAlignModInfo"></param>
        private void ResolveMSAlignModsWithModDefinitions(IReadOnlyCollection<ModificationDefinition> msAlignModInfo)
        {
            if (msAlignModInfo == null)
                return;

            foreach (var modInfo in msAlignModInfo)
            {
                if (string.IsNullOrEmpty(modInfo.TargetResidues))
                {
                    mPeptideMods.LookupModificationDefinitionByMassAndModType(
                        modInfo.ModificationMass, modInfo.ModificationType, default,
                        AminoAcidModInfo.ResidueTerminusStateConstants.None, out _, true);
                }
                else
                {
                    foreach (var chTargetResidue in modInfo.TargetResidues)
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(
                            modInfo.ModificationMass, modInfo.ModificationType, chTargetResidue,
                            AminoAcidModInfo.ResidueTerminusStateConstants.None, out _, true);
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<udtMSAlignSearchResultType> filteredSearchResults,
            bool includeSpeciesAndFragMethod,
            ref string errorLog)
        {
            // Sort filteredSearchResults by ascending PValue, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.PValueNum, item.ScanNum, item.ChargeNum, item.Peptide, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, includeSpeciesAndFragMethod, ref errorLog);
                index++;
            }
        }

        private void StoreSynMatches(
            IList<udtMSAlignSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            ICollection<udtMSAlignSearchResultType> filteredSearchResults)
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
        /// <param name="columnMapping"></param>
        /// <param name="includeSpeciesAndFragMethod"></param>
        /// <param name="errorLog"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            IDictionary<MSAlignResultsFileColumns, int> columnMapping,
            out bool includeSpeciesAndFragMethod,
            ref string errorLog)
        {
            includeSpeciesAndFragMethod = false;

            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = MSAlignSynFileReader.GetColumnHeaderNamesAndIDs();

                List<string> headerNames;

                // Only include Species_ID and FragMethod if defined in the column mapping
                columnMapping.TryGetValue(MSAlignResultsFileColumns.Species_ID, out var colIndexSpeciesId);
                columnMapping.TryGetValue(MSAlignResultsFileColumns.FragMethod, out var colIndexFragMethod);

                if (colIndexSpeciesId >= 0 || colIndexFragMethod >= 0)
                {
                    includeSpeciesAndFragMethod = true;
                    headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();
                }
                else
                {
                    headerNames = (from item in headerColumns
                                   where item.Key != MSAlignSynFileReader.DATA_COLUMN_Species_ID &&
                                         item.Key != MSAlignSynFileReader.DATA_COLUMN_FragMethod
                                   orderby item.Value
                                   select item.Key).ToList();
                }

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
        /// <param name="includeSpeciesAndFragMethod"></param>
        /// <param name="errorLog"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtMSAlignSearchResultType udtSearchResult,
            bool includeSpeciesAndFragMethod,
            ref string errorLog)
        {
            try
            {
                // Primary Columns
                //
                // MSAlign
                // ResultID  Scan  Prsm_ID  Spectrum_ID  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  Protein_Mass  Unexpected_Mod_Count  Peak_Count  Matched_Peak_Count  Matched_Fragment_Ion_Count  PValue  Rank_PValue  EValue  FDR

                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.ScanNum.ToString(),
                    udtSearchResult.Prsm_ID,
                    udtSearchResult.Spectrum_ID,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.DelM,
                    udtSearchResult.DelM_PPM,
                    udtSearchResult.MH,
                    udtSearchResult.Peptide,
                    udtSearchResult.Protein,
                    udtSearchResult.Protein_mass,
                    udtSearchResult.Unexpected_modifications,     // Unexpected_Mod_Count
                    udtSearchResult.Peaks,                        // Peak_count
                    udtSearchResult.Matched_peaks,                // Matched_Peak_Count
                    udtSearchResult.Matched_fragment_ions,        // Matched_Fragment_Ion_Count
                    udtSearchResult.Pvalue,
                    udtSearchResult.RankPValue.ToString(),
                    udtSearchResult.Evalue,
                    udtSearchResult.FDR
                };

                if (includeSpeciesAndFragMethod)
                {
                    data.Add(udtSearchResult.Species_ID);
                    data.Add(udtSearchResult.FragMethod);
                }

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

        #region "IComparer Classes"

        private class MSAlignSearchResultsComparerScanChargePValuePeptide : IComparer<udtMSAlignSearchResultType>
        {
            public int Compare(udtMSAlignSearchResultType x, udtMSAlignSearchResultType y)
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
                    result = string.CompareOrdinal(x.Peptide, y.Peptide);
                    if (result == 0)
                    {
                        // Peptide is the same, check Protein
                        result = string.CompareOrdinal(x.Protein, y.Protein);
                    }
                }
                return result;
            }
        }

        #endregion
    }
}