﻿// -------------------------------------------------------------------------------
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
    class clsTopPICResultsProcessor : clsPHRPBaseClass
    {
        public clsTopPICResultsProcessor()
        {
            mFileDate = "March 28, 2019";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_TopPIC_PROTEOFORMS_FILE = "_TopPIC_Proteoforms";

        public const string FILENAME_SUFFIX_TopPIC_PRSMs_FILE = "_TopPIC_PrSMs";

        public const string N_TERMINUS_SYMBOL_TopPIC = ".";

        public const string C_TERMINUS_SYMBOL_TopPIC = ".";

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string TopPIC_MOD_MASS_REGEX = @"\[([+-]*[0-9\.]+)\]";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to the tab-delimited file created directly by TopPIC
        /// </summary>
        public enum eTopPICResultsFileColumns
        {
            SpectrumFileName = 0,
            Prsm_ID = 1,
            Spectrum_ID = 2,
            FragMethod = 3,
            Scans = 4,
            Peaks = 5,
            Charge = 6,
            Precursor_mass = 7,              // Monoisotopic mass value of the observed precursor_mz
            Adjusted_precursor_mass = 8,     // Theoretical monoisotopic mass of the peptide (including mods)
            Proteoform_ID = 9,
            Feature_intensity = 10,
            Protein_name = 11,               // Protein name and description
            First_residue = 12,
            Last_residue = 13,
            Proteoform = 14,
            Unexpected_modifications = 15,
            Matched_peaks = 16,
            Matched_fragment_ions = 17,
            Pvalue = 18,
            Evalue = 19,
            Qvalue = 20,                     //  Spectral FDR, or PepFDR
            Proteoform_FDR = 21,
            Variable_PTMs = 22
        }

        #endregion

        #region "Structures"

        // This data structure holds rows read from the tab-delimited file created directly by TopPIC
        protected struct udtTopPICSearchResultType
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
            public string ResidueStart;                 // First_residue
            public string ResidueEnd;                   // Last_residue
            public string Proteoform;
            public string Unexpected_Mod_Count;         // unexpected modifications
            public string Matched_peaks;
            public string Matched_fragment_ions;
            public string Pvalue;
            public double PValueNum;
            public int RankPValue;
            public string Evalue;
            public string Qvalue;
            public string Proteoform_FDR;
            public string VariablePTMs;

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
                ResidueStart = string.Empty;
                ResidueEnd = string.Empty;
                Proteoform = string.Empty;
                Unexpected_Mod_Count = string.Empty;
                Matched_peaks = string.Empty;
                Matched_fragment_ions = string.Empty;
                Pvalue = string.Empty;
                PValueNum = 0;
                RankPValue = 0;
                Evalue = string.Empty;
                Qvalue = string.Empty;
                Proteoform_FDR = string.Empty;
                VariablePTMs = string.Empty;
            }
        }

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

            var parsingModMass = false;
            var modMassDigits = string.Empty;

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

                if (IsLetterAtoZ(chChar))
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
                    // Mod Mass Start
                    modMassDigits = string.Empty;
                    parsingModMass = true;
                }
                else if (chChar == ']')
                {
                    // Mod Mass End

                    if (parsingModMass)
                    {
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

                            var success = searchResult.SearchResultAddModification(modMass, residueForMod, residueLocForMod, searchResult.DetermineResidueTerminusState(residueLocForMod), updateModOccurrenceCounts);

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
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of searchResults so that we can sort by PValue

            var dctResultsSubset = new Dictionary<int, udtTopPICSearchResultType>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                dctResultsSubset.Add(index, searchResults[index]);
            }

            var resultsByProbability = (from item in dctResultsSubset orderby item.Value.PValueNum select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByProbability)
            {
                var oResult = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = oResult.PValueNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(oResult.PValueNum - lastValue) > double.Epsilon)
                    {
                        lastValue = oResult.PValueNum;
                        currentRank += 1;
                    }
                }

                oResult.RankPValue = currentRank;
                searchResults[entry.Key] = oResult;
            }
        }

        protected string AssureInteger(string integer, int defaultValue)
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

        protected double ComputePeptideMass(string peptide, double totalModMass)
        {
            var cleanSequence = GetCleanSequence(peptide);

            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);

            if (Math.Abs(totalModMass) > double.Epsilon)
            {
                mass += totalModMass;
            }

            return mass;
        }

        private static readonly Regex RegexModMassRegEx = new Regex(TopPIC_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form [23.5432]</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected double ComputeTotalModMass(string peptide)
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

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool MTS)
        {
            var pepToProteinMapFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
            if (pepToProteinMapFilePath.EndsWith("_toppic_syn", StringComparison.OrdinalIgnoreCase) ||
                pepToProteinMapFilePath.EndsWith("_toppic_fht", StringComparison.OrdinalIgnoreCase))
            {
                // Remove _syn or _fht
                pepToProteinMapFilePath = pepToProteinMapFilePath.Substring(0, pepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(pepToProteinMapFilePath, outputDirectoryPath, MTS);
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
        protected bool CreateSynResultsFile(
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
                {
                    using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        var headerParsed = false;

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
                                var success = ParseTopPICResultsFileHeaderLine(lineIn, columnMapping);
                                if (!success)
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
                            while (endIndex + 1 < searchResultsUnfiltered.Count && searchResultsUnfiltered[endIndex + 1].ScanNum == searchResultsUnfiltered[startIndex].ScanNum)
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

        protected bool ExtractModInfoFromTopPICParamFile(string topPICParamFilePath, ref List<clsModificationDefinition> modInfo)
        {
            var success = false;

            try
            {
                // Initialize the modification list
                if (modInfo == null)
                {
                    modInfo = new List<clsModificationDefinition>();
                }
                else
                {
                    modInfo.Clear();
                }

                if (string.IsNullOrEmpty(topPICParamFilePath))
                {
                    SetErrorMessage("TopPIC Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(topPICParamFilePath))
                {
                    SetErrorMessage("TopPIC param file not found: " + topPICParamFilePath);
                }
                else
                {
                    // Read the contents of the parameter (or mods) file
                    using (var reader = new StreamReader(new FileStream(topPICParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!reader.EndOfStream)
                        {
                            var lineIn = reader.ReadLine();
                            if (string.IsNullOrWhiteSpace(lineIn))
                                continue;

                            var dataLine = lineIn.Trim();
                            if (dataLine.Length <= 0)
                                continue;

                            if (dataLine.StartsWith("#"))
                            {
                                // Comment line; skip it
                            }
                            else
                            {
                                // Split the line on the equals sign
                                var kvSetting = clsPHRPParser.ParseKeyValueSetting(dataLine, '=', "#");

                                if (string.Equals(kvSetting.Key, "cysteineProtection", StringComparison.OrdinalIgnoreCase))
                                {
                                    clsModificationDefinition modDef;
                                    switch (kvSetting.Value.ToUpper())
                                    {
                                        case "C57":
                                            modDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                                                                   57.0215, "C",
                                                                                   clsModificationDefinition.eModificationTypeConstants.StaticMod,
                                                                                   "IodoAcet");
                                            modInfo.Add(modDef);
                                            break;

                                        case "C58":
                                            modDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                                                                   58.0055, "C",
                                                                                   clsModificationDefinition.eModificationTypeConstants.StaticMod,
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
                SetErrorMessage("Error reading the TopPIC parameter file (" + Path.GetFileName(topPICParamFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                success = false;
            }

            return success;
        }

        private void InitializeLocalVariables()
        {
            // Nothing to do at present
        }

        protected bool ParseTopPICSynopsisFile(string inputFilePath, string outputDirectoryPath, ref List<udtPepToProteinMappingType> pepToProteinMapping, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that TopPIC synopsis files are normally sorted on PValue value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
            // Although this was a possibility with Inspect, it likely never occurs for TopPIC
            //  But, we'll keep the check in place just in case

            var columnMapping = new Dictionary<clsPHRPParserTopPIC.TopPICSynFileColumns, int>();

            bool success;

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

                // Initialize htPeptidesFoundForPValueLevel
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
                        success = InitializeSequenceOutputFiles(baseOutputFilePath);

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
                                success = ParseTopPICSynFileHeaderLine(lineIn, columnMapping);
                                if (!success)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return success;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseTopPICSynFileEntry(lineIn, searchResult, ref errorLog,
                                                                            resultsProcessed, columnMapping,
                                                                            out var currentPeptideWithMods);

                            if (validSearchResult)
                            {
                                var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                                bool firstMatchForGroup;
                                if (searchResult.PValue == previousPValue)
                                {
                                    // New result has the same PValue as the previous result
                                    // See if htPeptidesFoundForPValueLevel contains the peptide, scan and charge

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
                                    // Reset htPeptidesFoundForPValueLevel
                                    peptidesFoundForPValueLevel.Clear();

                                    // Update previousPValue
                                    previousPValue = searchResult.PValue;

                                    // Append a new entry to htPeptidesFoundForPValueLevel
                                    peptidesFoundForPValueLevel.Add(key);
                                    firstMatchForGroup = true;
                                }

                                success = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                                if (!success)
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
                                        } while (pepToProteinMapIndex < pepToProteinMapping.Count && currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                                    }
                                    else
                                    {
                                        // Match not found; this is unexpected
                                        ReportWarning("no match for '" + currentPeptideWithMods + "' in pepToProteinMapping");
                                    }
                                }
                            }

                            // Update the progress
                            var percentComplete = Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100);
                            if (CreateProteinModsFile)
                            {
                                percentComplete = percentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                            }
                            UpdateProgress(percentComplete);

                            resultsProcessed += 1;
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

                    success = true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    success = false;
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
                success = false;
            }

            return success;
        }

        private bool ParseTopPICResultsFileEntry(
            string lineIn,
            ref udtTopPICSearchResultType udtSearchResult,
            ref string errorLog,
            IDictionary<eTopPICResultsFileColumns, int> columnMapping)
        {
            // Parses an entry from the TopPIC results file

            string[] splitLine = null;

            double precursorMZ = 0;

            bool validSearchResult;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                udtSearchResult.Clear();
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length >= 13)
                {
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
                            ReportWarning("Error parsing out the scan number from the scan list; could not find an integer: " + udtSearchResult.Scans);
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
                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Protein_name], out udtSearchResult.Protein);
                    udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);

                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.First_residue], out udtSearchResult.ResidueStart);
                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Last_residue], out udtSearchResult.ResidueEnd);

                    if (!GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Proteoform], out udtSearchResult.Proteoform))
                    {
                        ReportError("Peptide column is missing or invalid", true);
                    }

                    // Add the standard terminus symbols to the peptide sequence
                    udtSearchResult.Proteoform = ReplaceTerminus(udtSearchResult.Proteoform);

                    // Parse the sequence to determine the total mod mass
                    // Note that we do not remove any of the mod symbols since TopPIC identifies mods by mass alone, and since mods can ambiguously apply to residues
                    var totalModMass = ComputeTotalModMass(udtSearchResult.Proteoform);

                    // Compute theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                    var peptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Proteoform, totalModMass);

                    if (Math.Abs(peptideMonoMassTopPIC) < double.Epsilon)
                    {
                        peptideMonoMassTopPIC = peptideMonoMassPHRP;
                    }

                    var massDiffThreshold = peptideMonoMassTopPIC / 50000;
                    if (massDiffThreshold < 0.1)
                        massDiffThreshold = 0.1;

                    if (Math.Abs(peptideMonoMassPHRP - peptideMonoMassTopPIC) > massDiffThreshold)
                    {
                        // Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da or by a slightly larger value if over 5000 Da; this is unexpected
                        string first30Residues;
                        if (udtSearchResult.Proteoform.Length < 27)
                        {
                            first30Residues = udtSearchResult.Proteoform;
                        }
                        else
                        {
                            first30Residues = udtSearchResult.Proteoform.Substring(0, 27) + "...";
                        }
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " + massDiffThreshold.ToString("0.00") + " Da away from the mass computed by TopPIC: " + peptideMonoMassPHRP.ToString("0.0000") + " vs. " + peptideMonoMassTopPIC.ToString("0.0000") + "; peptide " + first30Residues);
                    }

                    if (peptideMonoMassTopPIC > 0)
                    {
                        // Compute DelM and DelM_PPM
                        var delM = precursorMonoMass - peptideMonoMassTopPIC;
                        udtSearchResult.DelM = MassErrorToString(delM);

                        if (precursorMZ > 0)
                        {
                            udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(clsPeptideMassCalculator.MassToPPM(delM, precursorMZ), 5, 0.00005);
                        }
                        else
                        {
                            udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(clsPeptideMassCalculator.MassToPPM(delM, 1000), 5, 0.00005);
                        }
                    }

                    // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoMassPHRP, 0), 6);

                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Unexpected_modifications], out udtSearchResult.Unexpected_Mod_Count);
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

                    if (udtSearchResult.Qvalue.ToLower() == "infinity")
                    {
                        udtSearchResult.Qvalue = "10";
                    }
                    else if (!string.IsNullOrEmpty(udtSearchResult.Qvalue) & !double.TryParse(udtSearchResult.Qvalue, out _))
                    {
                        udtSearchResult.Qvalue = string.Empty;
                    }

                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Proteoform_FDR], out udtSearchResult.Proteoform_FDR);
                    GetColumnValue(splitLine, columnMapping[eTopPICResultsFileColumns.Variable_PTMs], out udtSearchResult.VariablePTMs);

                    validSearchResult = true;
                }
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
                validSearchResult = false;
            }

            return validSearchResult;
        }

        private bool ParseTopPICResultsFileHeaderLine(string lineIn, IDictionary<eTopPICResultsFileColumns, int> columnMapping)
        {
            // Parse the header line

            // The expected header:
            // Data file name    Prsm ID    Spectrum ID    Fragmentation    Scan(s)    #peaks    Charge    Precursor mass    Adjusted precursor mass    Proteoform ID    Feature intensity    Protein name    First residue    Last residue    Proteoform    #unexpected modifications    #matched peaks    #matched fragment ions    P-value    E-value    Q-value (spectral FDR)    Proteoform FDR    #Variable PTMs

            var columnNames = new SortedDictionary<string, eTopPICResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Data file name", eTopPICResultsFileColumns.SpectrumFileName},
                {"Prsm ID", eTopPICResultsFileColumns.Prsm_ID},
                {"Spectrum ID", eTopPICResultsFileColumns.Spectrum_ID},
                {"Fragmentation", eTopPICResultsFileColumns.FragMethod},
                {"Scan(s)", eTopPICResultsFileColumns.Scans},
                {"#peaks", eTopPICResultsFileColumns.Peaks},
                {"Charge", eTopPICResultsFileColumns.Charge},
                {"Precursor mass", eTopPICResultsFileColumns.Precursor_mass},
                {"Adjusted precursor mass", eTopPICResultsFileColumns.Adjusted_precursor_mass},
                {"Proteoform ID", eTopPICResultsFileColumns.Proteoform_ID},
                {"Feature intensity", eTopPICResultsFileColumns.Feature_intensity},
                {"Protein name", eTopPICResultsFileColumns.Protein_name},
                {"First residue", eTopPICResultsFileColumns.First_residue},
                {"Last residue", eTopPICResultsFileColumns.Last_residue},
                {"Proteoform", eTopPICResultsFileColumns.Proteoform},
                {"#unexpected modifications", eTopPICResultsFileColumns.Unexpected_modifications},
                {"#matched peaks", eTopPICResultsFileColumns.Matched_peaks},
                {"#matched fragment ions", eTopPICResultsFileColumns.Matched_fragment_ions},
                {"P-value", eTopPICResultsFileColumns.Pvalue},
                {"E-value", eTopPICResultsFileColumns.Evalue},
                {"Q-value (spectral FDR)", eTopPICResultsFileColumns.Qvalue},
                {"Proteoform FDR", eTopPICResultsFileColumns.Proteoform_FDR},
                {"#Variable PTMs", eTopPICResultsFileColumns.Variable_PTMs}
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

        private bool ParseTopPICSynFileHeaderLine(string lineIn, IDictionary<clsPHRPParserTopPIC.TopPICSynFileColumns, int> columnMapping)
        {
            // Parse the header line

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

        private bool ParseTopPICSynFileEntry(
            string lineIn,
            clsSearchResultsTopPIC searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<clsPHRPParserTopPIC.TopPICSynFileColumns, int> columnMapping,
            out string peptideSequenceWithMods)
        {
            // Parses an entry from the TopPIC Synopsis file

            string[] splitLine = null;

            bool validSearchResult;

            // Reset searchResult
            searchResult.Clear();
            peptideSequenceWithMods = string.Empty;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length >= 15)
                {
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

                    // Calling this function will set .ProteoformPreResidues, .ProteoformPostResidues, .PeptideSequenceWithMods, and .ProteoformCleanSequence
                    searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

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

                    validSearchResult = true;
                }
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
                validSearchResult = false;
            }

            return validSearchResult;
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

                    var msAlignModInfo = new List<clsModificationDefinition>();
                    var pepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the TopPIC Parameter File so that we can determine whether Cysteine residues are statically modified
                    ExtractModInfoFromTopPICParamFile(SearchToolParameterFilePath, ref msAlignModInfo);

                    // Resolve the mods in mSAlignModInfo with the ModDefs mods
                    ResolveTopPICModsWithModDefinitions(ref msAlignModInfo);

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_TopPIC_ResultTable" with "_msalign"
                    if (baseName.EndsWith("_TopPIC_ResultTable", StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - "_TopPIC_ResultTable".Length) + "_msalign";
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

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseName, outputDirectoryPath, MTS: true);

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
                    // We only do this since a small number of peptides reported by TopPIC don't perfectly match the fasta file
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

        protected string ReplaceTerminus(string peptide)
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

        protected void ResolveTopPICModsWithModDefinitions(ref List<clsModificationDefinition> mSAlignModInfo)
        {
            if (mSAlignModInfo != null)
            {
                // Call .LookupModificationDefinitionByMass for each entry in mSAlignModInfo
                foreach (var modInfo in mSAlignModInfo)
                {
                    if (string.IsNullOrEmpty(modInfo.TargetResidues))
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(
                            modInfo.ModificationMass, modInfo.ModificationType, default(char),
                            clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out _, true);
                    }
                    else
                    {
                        foreach (var chTargetResidue in modInfo.TargetResidues)
                        {
                            mPeptideMods.LookupModificationDefinitionByMassAndModType(
                                modInfo.ModificationMass, modInfo.ModificationType, chTargetResidue,
                                clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out _, true);
                        }
                    }
                }
            }
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
                    udtSearchResult.FragMethod,
                    udtSearchResult.VariablePTMs,
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

        #region "IComparer Classes"

        protected class TopPICSearchResultsComparerScanChargePValuePeptide : IComparer<udtTopPICSearchResultType>
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