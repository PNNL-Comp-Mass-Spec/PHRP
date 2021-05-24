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
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM.AppSettings;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads in a MSPathFinder results file (_IcTda.tsv) and creates
    /// a tab-delimited text file with the data.
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads MSPathFinder results file Dataset_IcTda.tsv
    /// </para>
    /// <para>
    /// 2) It calls CreateSynResultsFile to create the _syn.txt file
    /// </para>
    /// <para>
    /// 3) ParseMSPathFinderResultsFileHeaderLine reads the header line to determine the column mapping
    ///      columnMapping = new Dictionary of MSPathFinderResultsFileColumns, int
    /// </para>
    /// <para>
    /// 4) ParseMSPathFinderResultsFileEntry reads each data line and stores in an instance of MSPathFinderSearchResult, which is a private struct
    ///    The data is stored in a list
    ///      searchResultsUnfiltered = new List of MSPathFinderSearchResult
    /// </para>
    /// <para>
    /// 5) Once the entire .tsv has been read, searchResultsUnfiltered is sorted by scan, charge, and ascending SpecEValue
    /// </para>
    /// <para>
    /// 6) StoreSynMatches stores filter-passing values in a new list
    ///      filteredSearchResults = new List of MSPathFinderSearchResult
    /// </para>
    /// <para>
    /// 7) SortAndWriteFilteredSearchResults performs one more sort, then writes out to disk
    ///    Sorts ascending by SpecEValue, QValue, Scan, Peptide, and Protein
    /// </para>
    /// </remarks>
    public class MSPathFinderResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: Dehydro, Desc, IcTda, mspath, Pre, struct

        /// <summary>
        /// Constructor
        /// </summary>
        public MSPathFinderResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "September 28, 2020";
        }

        public const string TOOL_NAME = "MSPathFinder";

        public const string FILENAME_SUFFIX_MSPathFinder_FILE = "_IcTda";

        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_MSPATHFINDER = "-";

        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_MSPATHFINDER = "-";

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        /// <summary>
        /// These columns correspond to the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
        /// </summary>
        private enum MSPathFinderResultsFileColumns
        {
            Scan = 0,
            PrefixResidue = 1,
            Sequence = 2,
            SuffixResidue = 3,
            Modifications = 4,
            Composition = 5,
            Protein = 6,
            ProteinDesc = 7,
            ProteinLength = 8,
            ResidueStart = 9,
            ResidueEnd = 10,
            Charge = 11,
            MostAbundantIsotopeMz = 12,
            CalculatedMonoMass = 13,
            MS1Features = 14,               // Column added 2017-10-14
            NumMatchedFragments = 15,
            Probability = 16,               // Column added 2015-11-18
            SpecEValue = 17,                // Column added 2015-08-25
            EValue = 18,                    // Column added 2015-08-25
            QValue = 19,
            PepQValue = 20
        }

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
        /// </summary>
        /// <remarks>
        /// These columns hold data that this class will use when creating the synopsis file
        /// </remarks>
        private struct MSPathFinderSearchResult
        {
            public string Scan;
            public int ScanNum;
            public string PrefixResidue;
            public string Sequence;
            public string SuffixResidue;
            public string Modifications;
            public string Composition;
            public string Protein;
            public string ProteinDesc;
            public string ProteinLength;
            public string ResidueStart;
            public string ResidueEnd;
            public string Charge;
            public short ChargeNum;
            public string MostAbundantIsotopeMz;        // As reported by MSPathfinder
            public string CalculatedMonoMass;           // Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MSPathFinder
            public double CalculatedMonoMassPHRP;       // Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by PHRP
            public string NumMatchedFragments;
            public string SpecEValue;                   // EValue, at the scan level
            public double SpecEValueNum;
            public string EValue;                       // EValue, at the peptide level
            public string QValue;                       // FDR, at the scan level
            public double QValueNum;
            public string PepQValue;                    // FDR, at the peptide level

            // The following are typically defined for other search engines, but are not used for MSPathFinder
            //   Public DelM As String                   ' Computed using Precursor_mass - CalculatedMonoMass
            //   Public DelM_PPM As String               ' Computed using DelM and CalculatedMonoMass

            public void Clear()
            {
                Scan = string.Empty;
                ScanNum = 0;
                PrefixResidue = string.Empty;
                Sequence = string.Empty;
                SuffixResidue = string.Empty;
                Modifications = string.Empty;
                Composition = string.Empty;
                Protein = string.Empty;
                ProteinDesc = string.Empty;
                ProteinLength = string.Empty;
                ResidueStart = string.Empty;
                ResidueEnd = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                MostAbundantIsotopeMz = string.Empty;
                CalculatedMonoMass = string.Empty;
                CalculatedMonoMassPHRP = 0;
                NumMatchedFragments = string.Empty;
                SpecEValue = string.Empty;
                SpecEValueNum = 0;
                EValue = string.Empty;
                QValue = string.Empty;
                PepQValue = string.Empty;

                // Unused at present: MH = string.Empty
                // Unused at present: DelM = string.Empty
                // Unused at present: DelM_PPM = string.Empty
            }

            public override string ToString()
            {
                return string.Format("Scan {0}: {1}, SpecEValue {2}", Scan, Sequence, SpecEValue);
            }
        }

        private readonly Regex mModListModNameMatcher = new(@"(?<ModName>.+) (?<ResidueNumber>\d+)", RegexOptions.Compiled);

        /// <summary>
        /// Step through the Modifications and associate each modification with the residues
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modList"></param>
        private void AddModificationsToResidues(
            MSPathFinderResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(searchResult.Modifications))
            {
                return;
            }

            // Modifications are listed as a comma separated list for Mod name and residue number; examples:
            // Oxidation 21
            // Oxidation 11,Dehydro 12
            // Dehydro 1,Dehydro 4,Dehydro 7

            var mods = searchResult.Modifications.Split(',');
            var finalResidueLoc = searchResult.PeptideCleanSequence.Length;

            foreach (var modEntry in mods)
            {
                // Find the mod in the list of known modifications (loaded from the MSPathFinder parameter file)
                var matchFound = false;

                // Obtain the mod name, for example "Dehydro" from "Dehydro 52"
                var match = mModListModNameMatcher.Match(modEntry);

                if (!match.Success)
                {
                    ReportError("Invalid MSPathFinder mod entry format; must be a name then a space then a number: " + modEntry);
                    continue;
                }

                var modName = match.Groups["ModName"].Value;
                var residueNumber = match.Groups["ResidueNumber"].Value;

                foreach (var modDef in modList)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.OrdinalIgnoreCase))
                    {
                        if (!int.TryParse(residueNumber, out var residueLocInPeptide))
                        {
                            ReportError("Mod entry does not have a number after the name: " + modEntry);
                            continue;
                        }

                        var residueTerminusState = AminoAcidModInfo.ResidueTerminusState.None;
                        if (residueLocInPeptide <= 1)
                        {
                            residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                        }
                        else if (residueLocInPeptide >= finalResidueLoc)
                        {
                            residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                        }

                        // Now that we know the terminus position, assure that residueLocInPeptide is 1 not 0
                        if (residueLocInPeptide < 1)
                        {
                            residueLocInPeptide = 1;
                        }
                        else if (residueLocInPeptide > finalResidueLoc)
                        {
                            residueLocInPeptide = finalResidueLoc;
                        }

                        var chMostRecentResidue = '-';

                        if (residueLocInPeptide >= 1 && residueLocInPeptide <= finalResidueLoc)
                        {
                            chMostRecentResidue = searchResult.PeptideCleanSequence[residueLocInPeptide - 1];
                        }

                        // Associate the mod with the given residue
                        searchResult.SearchResultAddModification(
                            modDef.ModMassVal, chMostRecentResidue, residueLocInPeptide,
                            residueTerminusState, updateModOccurrenceCounts);

                        matchFound = true;
                        break;
                    }
                }

                if (!matchFound)
                {
                    ReportError("Mod name " + modName + " was not defined in the MSPathFinder parameter file; cannot determine mod mass");
                }
            }
        }

        private bool AddModificationsAndComputeMass(
            MSPathFinderResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            bool success;

            try
            {
                // For other tools, we would add IsotopicMods here
                // This is not supported for MSPathFinder
                //
                // searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts)

                // Parse .Modifications to determine the modified residues present
                AddModificationsToResidues(searchResult, updateModOccurrenceCounts, modList);

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

        private double ComputePeptideMass(string cleanSequence, double totalModMass)
        {
            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);
            mass += totalModMass;

            return mass;
        }

        /// <summary>
        /// Computes the total of all modifications defined for the sequence
        /// </summary>
        /// <param name="modificationList">Comma separated list of modifications, e.g. Dehydro 52,Dehydro 63</param>
        /// <param name="modList"></param>
        private double ComputeTotalModMass(
            string modificationList,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(modificationList))
            {
                return 0;
            }

            double totalModMass = 0;
            var mods = modificationList.Split(',');

            foreach (var modEntry in mods)
            {
                // Convert the mod name to a mass value
                var matchFound = false;

                // Obtain the mod name, for example "Dehydro" from "Dehydro 52"
                var match = mModListModNameMatcher.Match(modEntry);

                if (!match.Success)
                {
                    ReportError("Mod entry does not have a name separated by a number: " + modEntry, true);
                }

                var modName = match.Groups["ModName"].Value;

                foreach (var modDef in modList)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.OrdinalIgnoreCase))
                    {
                        totalModMass += modDef.ModMassVal;
                        matchFound = true;
                        break;
                    }
                }

                if (!matchFound)
                {
                    ReportError("Mod name " + modName + " was not defined in the MSPathFinder parameter file; cannot determine the mod mass", true);
                }
            }

            return totalModMass;
        }

        private bool CreateProteinModsFileWork(
            string baseName,
            FileInfo inputFile,
            string synOutputFilePath,
            string outputDirectoryPath)
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
                if (File.Exists(mtsPepToProteinMapFilePath) && Options.UseExistingMTSPepToProteinMapFile)
                {
                    success = true;
                }
                else
                {
                    success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                    if (!success)
                    {
                        OnWarningEvent("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (success)
            {
                if (inputFile.Directory == null)
                {
                    OnWarningEvent("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                }
                else if (string.IsNullOrWhiteSpace(synOutputFilePath))
                {
                    OnWarningEvent("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputDirectoryPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          PeptideHitResultTypes.MSPathFinder);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                return true;
            }

            return true;
        }

        /// <summary>
        /// This routine creates a synopsis file from the output from MSPathFinder
        /// The synopsis file includes every result with a probability above a set threshold
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="modList"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            try
            {
                var columnMapping = new Dictionary<MSPathFinderResultsFileColumns, int>();
                var errorMessages = new List<string>();

                // Open the input file and parse it
                // Initialize the stream reader and the stream writer
                using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                using var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                var headerParsed = false;
                var lineNumber = 0;

                // Initialize the list that will hold all of the records in the MSPathFinder result file
                var searchResultsUnfiltered = new List<MSPathFinderSearchResult>();

                // Initialize the list that will hold all of the records that will ultimately be written out to disk
                var filteredSearchResults = new List<MSPathFinderSearchResult>();

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
                        var success = ParseMSPathFinderResultsFileHeaderLine(lineIn, columnMapping);
                        if (!success)
                        {
                            if (string.IsNullOrEmpty(mErrorMessage))
                            {
                                SetErrorMessage("Invalid header line in " + Path.GetFileName(inputFilePath));
                            }

                            return false;
                        }

                        // Write the header line to the output file
                        WriteSynFHTFileHeader(writer, errorMessages);

                        headerParsed = true;
                        continue;
                    }

                    var validSearchResult = ParseMSPathFinderResultsFileEntry(
                        lineIn, out var udtSearchResult, errorMessages,
                        columnMapping, modList, lineNumber);

                    if (validSearchResult)
                    {
                        searchResultsUnfiltered.Add(udtSearchResult);
                    }

                    // Update the progress
                    UpdateSynopsisFileCreationProgress(reader);
                }

                // Sort the SearchResults by scan, ascending SpecEValue, and peptide
                searchResultsUnfiltered.Sort(new MSPathFinderSearchResultsComparerScanScorePeptide());

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
                SortAndWriteFilteredSearchResults(writer, filteredSearchResults, errorMessages);

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
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the MSPathFinder parameter file
        /// </summary>
        /// <param name="msPathFinderParamFilePath"></param>
        /// <param name="modList"></param>
        /// <returns>True on success, false if an error</returns>
        /// <remarks>The DMS-based parameter file for MSPathFinder uses the same formatting as MS-GF+</remarks>
        private bool ExtractModInfoFromParamFile(
            string msPathFinderParamFilePath,
            out List<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            var modFileProcessor = new MSGFPlusParamFileModExtractor(TOOL_NAME);
            RegisterEvents(modFileProcessor);

            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                msPathFinderParamFilePath,
                MSGFPlusParamFileModExtractor.ModSpecFormats.MSGFPlusAndMSPathFinder,
                out modList);

            if (!success || mErrorCode != PHRPErrorCode.NoError)
            {
                if (mErrorCode == PHRPErrorCode.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MSPathFinder parameter file");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modList, mPeptideMods);

            return true;
        }

        private int GetNumMatchesPerSpectrumToReport(string searchToolParameterFilePath)
        {
            try
            {
                if (string.IsNullOrEmpty(searchToolParameterFilePath))
                {
                    ReportError("MSPathFinder parameter file name not defined; cannot determine NumMatchesPerSpec");
                    return 0;
                }

                var paramFile = new FileInfo(searchToolParameterFilePath);
                if (!paramFile.Exists)
                {
                    ReportError("MSPathFinder parameter file not found: " + searchToolParameterFilePath);
                    return 0;
                }

                var paramFileReader = new KeyValueParamFileReader("MSPathFinder", searchToolParameterFilePath);
                RegisterEvents(paramFileReader);

                var success = paramFileReader.ParseKeyValueParameterFile(out var paramFileEntries);
                if (!success)
                {
                    ReportError("Error reading MSPathFinder parameter file in GetNumMatchesPerSpectrumToReport: " + paramFileReader.ErrorMessage);
                    return 0;
                }

                return KeyValueParamFileReader.GetParameterValue(paramFileEntries, "NumMatchesPerSpec", 1);
            }
            catch (Exception ex)
            {
                ReportError(string.Format(
                                "Error reading the MSPathFinder parameter file ({0}): {1}",
                                Path.GetFileName(searchToolParameterFilePath), ex.Message));
                return 0;
            }
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <param name="modList"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSPathfinderSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            bool resetMassCorrectionTagsAndModificationDefinitions,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            // Note that MSPathfinder synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            var columnMapping = new Dictionary<MSPathFinderSynFileColumns, int>();

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
                var searchResult = new MSPathFinderResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize two SortedSets that will be used to avoid double-counting the same PSM in the same scan
                // This is required since a PSM with multiple proteins will be listed on multiple lines in the synopsis file
                // Values are PeptideSequenceWithMods_Scan_Charge

                var peptidesFoundForSpecEValue = new SortedSet<string>();
                var peptidesFoundForQValue = new SortedSet<string>();

                var firstMatchForGroup = false;

                var previousSpecEValue = string.Empty;
                var previousQValue = string.Empty;

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
                            var validHeader = ParseMSPathFinderSynFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }
                            headerParsed = true;
                            continue;
                        }

                        var validSearchResult = ParseMSPathFinderSynFileEntry(
                            lineIn, searchResult, errorMessages,
                            resultsProcessed, columnMapping,
                            out _);

                        resultsProcessed++;
                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                        var newValue = true;

                        if (string.IsNullOrEmpty(searchResult.SpecEValue))
                        {
                            if (searchResult.QValue == previousQValue)
                            {
                                // New result has the same QValue as the previous result
                                // See if peptidesFoundForQValue contains the peptide, scan and charge

                                if (peptidesFoundForQValue.Contains(key))
                                {
                                    firstMatchForGroup = false;
                                }
                                else
                                {
                                    peptidesFoundForQValue.Add(key);
                                    firstMatchForGroup = true;
                                }

                                newValue = false;
                            }
                        }
                        else if (searchResult.SpecEValue == previousSpecEValue)
                        {
                            // New result has the same SpecEValue as the previous result
                            // See if peptidesFoundForSpecEValue contains the peptide, scan and charge

                            if (peptidesFoundForSpecEValue.Contains(key))
                            {
                                firstMatchForGroup = false;
                            }
                            else
                            {
                                peptidesFoundForSpecEValue.Add(key);
                                firstMatchForGroup = true;
                            }

                            newValue = false;
                        }

                        if (newValue)
                        {
                            // New SpecEValue or new QValue
                            // Reset the SortedSets
                            peptidesFoundForSpecEValue.Clear();
                            peptidesFoundForQValue.Clear();

                            // Update the cached values
                            previousSpecEValue = searchResult.SpecEValue;
                            previousQValue = searchResult.QValue;

                            // Append a new entry to the SortedSets
                            peptidesFoundForSpecEValue.Add(key);
                            peptidesFoundForQValue.Add(key);

                            firstMatchForGroup = true;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup, modList);
                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
                            errorMessages.Add(string.Format("Error adding modifications to sequence for ResultID '{0}'", searchResult.ResultID));
                        }

                        SaveResultsFileEntrySeqInfo(searchResult, firstMatchForGroup);

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
                            OnWarningEvent("ParseMSPathfinderSynopsisFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
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
                    SetErrorMessage("Error reading input file in ParseMSPathfinderSynopsisFile", ex);
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
                SetErrorMessage("Error creating the output file in ParseMSPathfinderSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parse a MSPathFinder results line while creating the MSPathFinder synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="columnMapping"></param>
        /// <param name="modList"></param>
        /// <param name="lineNumber">Line number in the input file (used for error reporting)</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSPathFinderResultsFileEntry(
            string lineIn,
            out MSPathFinderSearchResult udtSearchResult,
            ICollection<string> errorMessages,
            IDictionary<MSPathFinderResultsFileColumns, int> columnMapping,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList,
            int lineNumber)
        {
            udtSearchResult = new MSPathFinderSearchResult();

            try
            {
                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 11)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.Scan], out udtSearchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid on line " + lineNumber, true);
                }

                if (!int.TryParse(udtSearchResult.Scan, out udtSearchResult.ScanNum))
                {
                    ReportError("Scan column is not numeric on line " + lineNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(udtSearchResult.Charge, 0));

                // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MSPathFinder
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                // double.TryParse(udtSearchResult.CalculatedMonoMass, out var sequenceMonoMassMSPathFinder);

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.PrefixResidue], out udtSearchResult.PrefixResidue);

                if (!GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.Sequence], out udtSearchResult.Sequence))
                {
                    ReportError("Sequence column is missing or invalid on line " + lineNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.SuffixResidue], out udtSearchResult.SuffixResidue);

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.Modifications], out udtSearchResult.Modifications);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.Composition], out udtSearchResult.Composition);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.Protein], out udtSearchResult.Protein);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.ProteinDesc], out udtSearchResult.ProteinDesc);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.ProteinLength], out udtSearchResult.ProteinLength);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.ResidueEnd], out udtSearchResult.ResidueEnd);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.ResidueStart], out udtSearchResult.ResidueStart);

                // Parse the list of modified residues to determine the total mod mass
                var totalModMass = ComputeTotalModMass(udtSearchResult.Modifications, modList);

                // Compute monoisotopic mass of the peptide
                udtSearchResult.CalculatedMonoMassPHRP = ComputePeptideMass(udtSearchResult.Sequence, totalModMass);

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.MostAbundantIsotopeMz], out udtSearchResult.MostAbundantIsotopeMz);
                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.NumMatchedFragments], out udtSearchResult.NumMatchedFragments);

                if (GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.SpecEValue], out udtSearchResult.SpecEValue))
                {
                    double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum);
                }

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.EValue], out udtSearchResult.EValue);

                if (GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.QValue], out udtSearchResult.QValue))
                {
                    double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum);
                }

                GetColumnValue(splitLine, columnMapping[MSPathFinderResultsFileColumns.PepQValue], out udtSearchResult.PepQValue);

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the MSPathFinder results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing MSPathFinder results in ParseMSPathFinderResultsFileEntry, line {0}", lineNumber));
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the MSPathFinder results file header line, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseMSPathFinderResultsFileHeaderLine(string lineIn, IDictionary<MSPathFinderResultsFileColumns, int> columnMapping)
        {
            // The expected column order from MSPathFinder:
            //   Scan	Pre	Sequence	Post	Modifications	Composition	ProteinName	ProteinDesc	ProteinLength	Start	End	Charge	MostAbundantIsotopeMz	Mass	#MatchedFragments	Probability SpecEValue    EValue    QValue    PepQValue

            var columnNames = new SortedDictionary<string, MSPathFinderResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Scan", MSPathFinderResultsFileColumns.Scan},
                {"Pre", MSPathFinderResultsFileColumns.PrefixResidue},
                {"Sequence", MSPathFinderResultsFileColumns.Sequence},
                {"Post", MSPathFinderResultsFileColumns.SuffixResidue},
                {"Modifications", MSPathFinderResultsFileColumns.Modifications},
                {"Composition", MSPathFinderResultsFileColumns.Composition},
                {"ProteinName", MSPathFinderResultsFileColumns.Protein},
                {"ProteinDesc", MSPathFinderResultsFileColumns.ProteinDesc},
                {"ProteinLength", MSPathFinderResultsFileColumns.ProteinLength},
                {"Start", MSPathFinderResultsFileColumns.ResidueStart},
                {"End", MSPathFinderResultsFileColumns.ResidueEnd},
                {"Charge", MSPathFinderResultsFileColumns.Charge},
                {"MostAbundantIsotopeMz", MSPathFinderResultsFileColumns.MostAbundantIsotopeMz},
                {"Mass", MSPathFinderResultsFileColumns.CalculatedMonoMass},
                {"MS1Features", MSPathFinderResultsFileColumns.MS1Features},
                {"#MatchedFragments", MSPathFinderResultsFileColumns.NumMatchedFragments},
                {"Probability", MSPathFinderResultsFileColumns.Probability},
                {"SpecEValue", MSPathFinderResultsFileColumns.SpecEValue},
                {"EValue", MSPathFinderResultsFileColumns.EValue},
                {"QValue", MSPathFinderResultsFileColumns.QValue},
                {"PepQValue", MSPathFinderResultsFileColumns.PepQValue}
            };

            // These columns were in older versions of MSPathFinder
            // Ignore them
            var legacyColumnNames = new SortedSet<string>(StringComparer.OrdinalIgnoreCase)
            {
                "IsotopeCorrPrevMs1",
                "IsotopeCorrNextMs1",
                "CorrMostAbundantPlusOneIsotope",
                "ChargeCorrMinusOne",
                "ChargeCorrPlusOne"
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MSPathFinderResultsFileColumns resultColumn in Enum.GetValues(typeof(MSPathFinderResultsFileColumns)))
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
                            }
                            else if (!legacyColumnNames.Contains(splitLine[index]))
                            {
                                // Unrecognized column name
                                OnWarningEvent("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMSPathFinderResultsFileHeaderLine");
                            }
                        }
                    }
                }

                if (useDefaultHeaders)
                {
                    // Use default column mappings
                    foreach (MSPathFinderResultsFileColumns resultColumn in Enum.GetValues(typeof(MSPathFinderResultsFileColumns)))
                    {
                        columnMapping[resultColumn] = (int)resultColumn;
                    }

                    // This is not a header line; return false
                    return false;
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the MSPathFinder results file", ex);
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a MSPathFinder _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSPathFinderSynFileHeaderLine(string lineIn, IDictionary<MSPathFinderSynFileColumns, int> columnMapping)
        {
            var columnNames = MSPathFinderSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MSPathFinderSynFileColumns resultColumn in Enum.GetValues(typeof(MSPathFinderSynFileColumns)))
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
                SetErrorMessage("Error parsing the header line in the MSPathFinder synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a MSPathFinder Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequence"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMSPathFinderSynFileEntry(
            string lineIn,
            MSPathFinderResults searchResult,
            ICollection<string> errorMessages,
            int resultsProcessed,
            IDictionary<MSPathFinderSynFileColumns, int> columnMapping,
            out string peptideSequence)
        {
            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();
            peptideSequence = string.Empty;

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 15)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.ResultID], out string value))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading ResultID value from MSPathFinder results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.Sequence], out peptideSequence))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading Peptide sequence value from MSPathFinder results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.Protein], out string proteinName);
                searchResult.ProteinName = proteinName;
                searchResult.MultipleProteinCount = "0";

                searchResult.PeptideDeltaMass = "0";

                // Note that MSPathFinder sequences don't actually have mod symbols; that information is tracked via searchResult.Modifications

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.MostAbundantIsotopeMz], out string mostAbundantIsotopeMz);

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.Modifications], out string modifications);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.Composition], out string composition);

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.ProteinDesc], out string proteinDesc);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.ProteinLength], out string proteinLength);

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.ResidueStart], out string residueStart);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.ResidueEnd], out string residueEnd);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.MatchedFragments], out string matchedFragments);

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.SpecEValue], out string specEValue);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.EValue], out string eValue);

                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.QValue], out string qValue);
                GetColumnValue(splitLine, columnMapping[MSPathFinderSynFileColumns.PepQValue], out string pepQValue);

                searchResult.MostAbundantIsotopeMz = mostAbundantIsotopeMz;
                searchResult.Modifications = modifications;
                searchResult.Composition = composition;
                searchResult.ProteinDesc = proteinDesc;
                searchResult.ProteinLength = proteinLength;
                searchResult.ResidueStart = residueStart;
                searchResult.ResidueEnd = residueEnd;
                searchResult.MatchedFragments = matchedFragments;
                searchResult.SpecEValue = specEValue;
                searchResult.EValue = eValue;
                searchResult.QValue = qValue;
                searchResult.PepQValue = pepQValue;

                // Update the peptide location in the protein
                if (!string.IsNullOrEmpty(searchResult.ResidueStart))
                {
                    int.TryParse(searchResult.ResidueStart, out var peptideLocInProteinStart);
                    searchResult.PeptideLocInProteinStart = peptideLocInProteinStart;
                }

                if (!string.IsNullOrEmpty(searchResult.ResidueEnd))
                {
                    int.TryParse(searchResult.ResidueEnd, out var peptideLocInProteinEnd);
                    searchResult.PeptideLocInProteinEnd = peptideLocInProteinEnd;
                }

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing MSPathFinder results for RowIndex '{0}'", splitLine[0]));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing MSPathFinder Results in ParseMSPathFinderSynFileEntry");
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MSPathFinder results file (Dataset_IcTda.tsv)</param>
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

                    // Load the MSPathFinder Parameter File so that we can determine the modification names and masses
                    var modInfoExtracted = ExtractModInfoFromParamFile(Options.SearchToolParameterFilePath, out var modList);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Re-parse the MSPathFinder parameter file to look for NumMatchesPerSpec
                    var numMatchesPerSpec = GetNumMatchesPerSpectrumToReport(Options.SearchToolParameterFilePath);
                    if (numMatchesPerSpec > 1)
                    {
                        // Auto-change IgnorePeptideToProteinMapperErrors to True
                        // since the MSPathFinder parameter file has NumMatchesPerSpec of 2 or higher
                        // (which results in PSMs associated with decoy proteins)
                        Options.IgnorePeptideToProteinMapperErrors = true;

                        OnDebugEvent(string.Format(
                            "Set IgnorePeptideToProteinMapperErrors to true since NumMatchesPerSpec is {0} in the MSPathFinder parameter file",
                            numMatchesPerSpec));
                    }

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_IcTda" with "_mspath"
                    if (baseName.EndsWith(FILENAME_SUFFIX_MSPathFinder_FILE, StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - FILENAME_SUFFIX_MSPathFinder_FILE.Length) + "_mspath";
                    }

                    // Do not create a first-hits file for MSPathFinder results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_mspath_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath, modList);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to create the other PHRP files
                    success = ParseMSPathfinderSynopsisFile(synOutputFilePath, outputDirectoryPath, false, modList);

                    if (success && Options.CreateProteinModsFile)
                    {
                        // Check for an empty synopsis file
                        if (!ValidateFileHasData(synOutputFilePath, "Synopsis file", out var errorMessage))
                        {
                            OnWarningEvent(errorMessage);
                        }
                        else
                        {
                            success = CreateProteinModsFileWork(baseName, inputFile, synOutputFilePath, outputDirectoryPath);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in MSPathFinderResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<MSPathFinderSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            // Sort filteredSearchResults by ascending SpecEValue, QValue, Scan, Peptide, and Protein
            var query = from item in filteredSearchResults orderby item.SpecEValueNum, item.QValueNum, item.ScanNum, item.Sequence, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, errorMessages);
                index++;
            }
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan (typically there will only be one result for MSPathFinder)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Output parameter: the actual filtered search results</param>
        private void StoreSynMatches(
            IList<MSPathFinderSearchResult> searchResults,
            int startIndex,
            int endIndex,
            List<MSPathFinderSearchResult> filteredSearchResults)
        {
            // If there was more than one result, we could rank them by score
            // AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex)

            // The calling procedure already sorted by scan, charge, and Score; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            // Now store the matches that pass the filters
            //  Either SpecEValue < 5E-07 (0.0000005)
            //  or     QValue < 10
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].SpecEValueNum <= Options.MSGFPlusSynopsisFileSpecEValueThreshold || searchResults[index].QValueNum < 0.1)
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
                var headerColumns = MSPathFinderSynFileReader.GetColumnHeaderNamesAndIDs();

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
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorMessages"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            MSPathFinderSearchResult udtSearchResult,
            ICollection<string> errorMessages)
        {
            try
            {
                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.Scan,
                    udtSearchResult.Charge,
                    udtSearchResult.MostAbundantIsotopeMz,
                    udtSearchResult.CalculatedMonoMass,
                    udtSearchResult.PrefixResidue + "." + udtSearchResult.Sequence + "." + udtSearchResult.SuffixResidue,
                    udtSearchResult.Modifications,
                    udtSearchResult.Composition,
                    udtSearchResult.Protein,
                    udtSearchResult.ProteinDesc,
                    udtSearchResult.ProteinLength,
                    udtSearchResult.ResidueStart,
                    udtSearchResult.ResidueEnd,
                    udtSearchResult.NumMatchedFragments,
                    udtSearchResult.SpecEValue,
                    udtSearchResult.EValue,
                    udtSearchResult.QValue,
                    udtSearchResult.PepQValue
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

        private class MSPathFinderSearchResultsComparerScanScorePeptide : IComparer<MSPathFinderSearchResult>
        {
            public int Compare(MSPathFinderSearchResult x, MSPathFinderSearchResult y)
            {
                // First sort on Scan
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                // Scan is the same; check SpecEValue
                if (x.SpecEValueNum > y.SpecEValueNum)
                {
                    return 1;
                }

                if (x.SpecEValueNum < y.SpecEValueNum)
                {
                    return -1;
                }

                // SpecEValue is the same; check QValue
                if (x.QValueNum > y.QValueNum)
                {
                    return 1;
                }

                if (x.QValueNum < y.QValueNum)
                {
                    return -1;
                }

                // SpecEValue is the same; check sequence
                var result = string.CompareOrdinal(x.Sequence, y.Sequence);
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
