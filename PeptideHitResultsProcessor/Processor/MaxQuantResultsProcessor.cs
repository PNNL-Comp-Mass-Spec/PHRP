// This class reads in a MaxQuant msms.txt results file and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.
//
// The modification definition information is determined from the MaxQuant parameter file
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
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM.AppSettings;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads in MaxQuant results file msms.txt and creates
    /// a tab-delimited text file with the data.
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads MaxQuant results file Dataset_IcTda.tsv
    /// </para>
    /// <para>
    /// 2) It calls CreateSynResultsFile to create the _syn.txt file
    /// </para>
    /// <para>
    /// 3) ParseMaxQuantResultsFileHeaderLine reads the header line to determine the column mapping
    ///      columnMapping = new Dictionary of MaxQuantResultsFileColumns, int
    /// </para>
    /// <para>
    /// 4) ParseMaxQuantResultsFileEntry reads each data line and stores in an instance of MaxQuantSearchResult, which is a private struct
    ///    The data is stored in a list
    ///      searchResultsUnfiltered = new List of MaxQuantSearchResult
    /// </para>
    /// <para>
    /// 5) Once the entire .tsv has been read, searchResultsUnfiltered is sorted by scan, charge, and ascending SpecEValue
    /// </para>
    /// <para>
    /// 6) StoreSynMatches stores filter-passing values in a new list
    ///      filteredSearchResults = new List of MaxQuantSearchResult
    /// </para>
    /// <para>
    /// 7) SortAndWriteFilteredSearchResults performs one more sort, then writes out to disk
    ///    Sorts ascending by SpecEValue, QValue, Scan, Peptide, and Protein
    /// </para>
    /// </remarks>
    public class MaxQuantResultsProcessor : PHRPBaseClass
    {
        // ReSharper disable CommentTypo

        // Ignore Spelling: acetyl, Carbamidomethyl, Desc, DimethNter, Glu, Gln, MaxQuant, plex, pyro, struct

        // ReSharper restore CommentTypo

        /// <summary>
        /// Constructor
        /// </summary>
        public MaxQuantResultsProcessor()
        {
            FileDate = "April 3, 2021";

            mModCountMatcher = new Regex(@"^(?<ModCount>\d+) (?<ModName>.+)", RegexOptions.Compiled);

            // ReSharper disable CommentTypo

            // The following RegEx matches mod names like these:
            // DimethNter0
            // DimethNter4
            // DimethNter8
            // ICPL-Nter10
            // mTRAQ-Nter8
            // DimethNter2
            // DimethNter6
            // iTRAQ4plex-Nter114
            // iTRAQ8plex-Nter121
            // TMT2plex-Nter126
            // TMT8plex-Nter129N
            // TMT8plex-Nter129C
            // TMT11plex-Nter131C
            mNTermModMatcher = new Regex(@"Nter\d+[NC]*$", RegexOptions.Compiled);

            // ReSharper restore CommentTypo
        }

        public const string TOOL_NAME = "MaxQuant";

        public const string MSMS_FILE_NAME = "msms.txt";

        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_MaxQuant = "_";

        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_MaxQuant = "_";

        private const string SEARCH_ENGINE_NAME = "MaxQuant";

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        /// <summary>
        /// These columns correspond to MaxQuant file msms.txt
        /// </summary>
        private enum MaxQuantResultsFileColumns
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
            MS1Features = 14,
            NumMatchedFragments = 15,
            Probability = 16,
            SpecEValue = 17,
            EValue = 18,
            QValue = 19,
            PepQValue = 20
        }

        /// <summary>
        /// This data structure holds rows read from MaxQuant file msms.txt
        /// </summary>
        /// <remarks>
        /// These columns correspond to the Synopsis file created by this class
        /// </remarks>
        private struct MaxQuantSearchResult
        {
            // ToDo: Customize these for MaxQuant
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
            public string MostAbundantIsotopeMz;        // As reported by MaxQuant
            public string CalculatedMonoMass;           // Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MaxQuant
            public double CalculatedMonoMassPHRP;       // Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by PHRP
            public string NumMatchedFragments;
            public string SpecEValue;                   // EValue, at the scan level
            public double SpecEValueNum;
            public string EValue;                       // EValue, at the peptide level
            public string QValue;                       // FDR, at the scan level
            public double QValueNum;
            public string PepQValue;                    // FDR, at the peptide level

            public string DelM;                         // Computed using Precursor_mass - CalculatedMonoMass
            public string DelM_PPM;                     // Computed using DelM and CalculatedMonoMass

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

                // Unused at present: MH = string.Empty;
                DelM = string.Empty;
                DelM_PPM = string.Empty;
            }
        }

        private readonly Regex mModCountMatcher;

        private readonly Regex mNTermModMatcher;

        /// <summary>
        /// Step through the Modifications and associate each modification with the residues
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modInfo"></param>
        private void AddModificationsToResidues(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            if (string.IsNullOrWhiteSpace(searchResult.Modifications) || searchResult.Modifications.Equals("Unmodified"))
            {
                return;
            }

            var modsWithCounts = searchResult.Modifications.Split(',');
            var finalResidueLoc = searchResult.PeptideCleanSequence.Length;

            var modNames = new SortedSet<string>();

            foreach (var modItem in modsWithCounts)
            {
                // Check whether the mod name is preceded by an integer (which indicates the number of times this peptide has this modification)
                var match = mModCountMatcher.Match(modItem);

                var modName = match.Success ? match.Groups["ModName"].Value : modItem;

                if (!modNames.Contains(modName))
                    modNames.Add(modName);
            }

            // Keys are mod names, values are index in searchResult.ModifiedSequence
            var modPositions = new Dictionary<string, int>();

            while (true)
            {
                // Look for the earliest occurrence of any of the mods in searchResult.ModifiedSequence

                foreach (var mod in modNames)
                {
                    var charIndex = searchResult.ModifiedSequence.IndexOf(string.Format("({0})", mod), StringComparison.Ordinal);
                    modPositions.Add(mod, charIndex);
                }

                var minIndex = int.MaxValue;
                var firstModName = string.Empty;

                foreach (var modName in modPositions.Keys)
                {
                    if (modPositions[modName] >= 0 && modPositions[modName] < minIndex)
                    {
                        minIndex = modPositions[modName];
                        firstModName = string.Copy(modName);
                    }
                }

                if (string.IsNullOrEmpty(firstModName))
                {
                    // No more matches
                    break;
                }

                // ToDo: possibly get this info from modifications.xml
                //       Look for any mod with <position>anyNterm</position> or <position>proteinNterm</position>

                var nTerminalMod = firstModName.EndsWith("N-term)") ||
                                   firstModName.Equals("Glu-&gt;pyro-Glu") ||
                                   firstModName.Equals("Gln-&gt;pyro-Glu") ||
                                   mNTermModMatcher.IsMatch(firstModName);

                if (nTerminalMod)
                {
                    // The mod applies to the residue just after the mod name (which is surrounded by parentheses)
                }
                else
                {
                    // The mod applies to the residue just before the mod name (which is surrounded by parentheses)
                }

                // ToDo: Associate the mod with the given residue
                // searchResult.SearchResultAddModification(
                //     modDef.ModMassVal, chMostRecentResidue, residueLocInPeptide,
                //     residueTerminusState, updateModOccurrenceCounts);

            }
        }

        private bool AddModificationsAndComputeMass(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            bool success;

            try
            {
                // Some of the other tools add IsotopicMods here
                // This is not supported for MaxQuant
                //
                // searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts)

                // Parse .Modifications to determine the modified residues present
                AddModificationsToResidues(searchResult, updateModOccurrenceCounts, modInfo);

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
        /// <param name="modificationList">Comma separated list of modifications, e.g. "Acetyl (Protein N-term),Oxidation (M)" or "2 Oxidation (M) "</param>
        /// <param name="modInfo"></param>
        private double ComputeTotalModMass(
            string modificationList,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            if (string.IsNullOrWhiteSpace(modificationList))
            {
                return 0;
            }

            double totalModMass = 0;
            var mods = modificationList.Split(',');

            foreach (var modItem in mods)
            {
                // Convert the mod name to a mass value
                var matchFound = false;

                // Obtain the mod name, for example "Oxidation (M)" from "2 Oxidation (M)"

                var match = mModCountMatcher.Match(modItem);

                var modName = match.Success ? match.Groups["ModName"].Value : modItem;

                foreach (var modDef in modInfo)
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
                    ReportError("Mod name " + modName + " was not defined in the MaxQuant parameter file; cannot determine the mod mass", true);
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
                if (File.Exists(mtsPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
                {
                    success = true;
                }
                else
                {
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
                    // If necessary, copy various PHRPReader support files to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputDirectoryPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          PeptideHitResultTypes.MaxQuant);
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
        /// This routine creates a synopsis file from the output from MaxQuant
        /// The synopsis file includes every result with a probability above a set threshold
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            try
            {
                var columnMapping = new Dictionary<MaxQuantResultsFileColumns, int>();
                var errorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var headerParsed = false;
                    var rowNumber = 0;

                    // Initialize the list that will hold all of the records in the MaxQuant result file
                    var searchResultsUnfiltered = new List<MaxQuantSearchResult>();

                    // Initialize the list that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<MaxQuantSearchResult>();

                    // Parse the input file
                    while (!reader.EndOfStream && !AbortProcessing)
                    {
                        var lineIn = reader.ReadLine();
                        rowNumber++;

                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        if (!headerParsed)
                        {
                            // Parse the header line
                            var success = ParseMaxQuantResultsFileHeaderLine(lineIn, columnMapping);
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

                        var udtSearchResult = new MaxQuantSearchResult();

                        var validSearchResult =
                            ParseMaxQuantResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping, modInfo, rowNumber);

                        if (validSearchResult)
                        {
                            searchResultsUnfiltered.Add(udtSearchResult);
                        }

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    // Sort the SearchResults by scan, charge, and ascending SpecEValue
                    searchResultsUnfiltered.Sort(new MaxQuantSearchResultsComparerScanChargeScorePeptide());

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
                SetErrorMessage(ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the MaxQuant parameter file
        /// </summary>
        /// <param name="maxQuantParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromParamFile(
            string maxQuantParamFilePath,
            out List<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            var success = false;

            // ToDo: open the MaxQuant parameter file with an XML reader and look for static and dynamic mod names

            // <fixedModifications>
            //    <string>Carbamidomethyl (C)</string>
            // </fixedModifications>

            // <variableModifications>
            //    <string>Oxidation (M)</string>
            //    <string>Acetyl (Protein N-term)</string>
            // </variableModifications>

            modInfo = new List<MSGFPlusParamFileModExtractor.ModInfo>();

            if (!success || mErrorCode != PHRPErrorCode.NoError)
            {
                if (mErrorCode == PHRPErrorCode.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MaxQuant parameter file");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            var modFileProcessor = new MSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modInfo, mPeptideMods);

            return true;
        }

        private int GetNumMatchesPerSpectrumToReport(string searchToolParameterFilePath)
        {
            try
            {
                if (string.IsNullOrEmpty(searchToolParameterFilePath))
                {
                    ReportError("MaxQuant parameter file name not defined; cannot determine NumMatchesPerSpec");
                    return 0;
                }

                var paramFile = new FileInfo(searchToolParameterFilePath);
                if (!paramFile.Exists)
                {
                    ReportError("MaxQuant parameter file not found: " + searchToolParameterFilePath);
                    return 0;
                }

                var paramFileReader = new KeyValueParamFileReader("MaxQuant", searchToolParameterFilePath);
                RegisterEvents(paramFileReader);

                var success = paramFileReader.ParseKeyValueParameterFile(out var paramFileEntries);
                if (!success)
                {
                    ReportError("Error reading MaxQuant parameter file in GetNumMatchesPerSpectrumToReport: " + paramFileReader.ErrorMessage);
                    return 0;
                }

                return KeyValueParamFileReader.GetParameterValue(paramFileEntries, "NumMatchesPerSpec", 1);
            }
            catch (Exception ex)
            {
                ReportError(string.Format(
                                "Error reading the MaxQuant parameter file ({0}): {1}",
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
        /// <param name="modInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            bool resetMassCorrectionTagsAndModificationDefinitions,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modInfo)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that MaxQuant synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            var columnMapping = new Dictionary<MaxQuantSynFileColumns, int>();

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
                var searchResult = new MaxQuantResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize two SortedSets
                var peptidesFoundForSpecEValue = new SortedSet<string>();
                var peptidesFoundForQValue = new SortedSet<string>();

                var firstMatchForGroup = false;

                var previousSpecEValue = string.Empty;
                var previousQValue = string.Empty;

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
                                var validHeader = ParseMaxQuantSynFileHeaderLine(lineIn, columnMapping);
                                if (!validHeader)
                                {
                                    // Error parsing header
                                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMaxQuantSynFileEntry(lineIn, searchResult, ref errorLog,
                                                                                     resultsProcessed, columnMapping,
                                                                                     out _);

                            resultsProcessed++;
                            if (!validSearchResult)
                            {
                                continue;
                            }

                            var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                            var newValue = true;

                            // ToDo:
                            //if (string.IsNullOrEmpty(searchResult.SpecEValue))
                            //{
                            //    if (searchResult.QValue == previousQValue)
                            //    {
                            //        // New result has the same QValue as the previous result
                            //        // See if peptidesFoundForQValue contains the peptide, scan and charge

                            //        if (peptidesFoundForQValue.Contains(key))
                            //        {
                            //            firstMatchForGroup = false;
                            //        }
                            //        else
                            //        {
                            //            peptidesFoundForQValue.Add(key);
                            //            firstMatchForGroup = true;
                            //        }

                            //        newValue = false;
                            //    }
                            //}
                            //else if (searchResult.SpecEValue == previousSpecEValue)
                            //{
                            //    // New result has the same SpecEValue as the previous result
                            //    // See if peptidesFoundForSpecEValue contains the peptide, scan and charge

                            //    if (peptidesFoundForSpecEValue.Contains(key))
                            //    {
                            //        firstMatchForGroup = false;
                            //    }
                            //    else
                            //    {
                            //        peptidesFoundForSpecEValue.Add(key);
                            //        firstMatchForGroup = true;
                            //    }

                            //    newValue = false;
                            //}

                            if (newValue)
                            {
                                // New SpecEValue or new QValue
                                // Reset the SortedSets
                                peptidesFoundForSpecEValue.Clear();
                                peptidesFoundForQValue.Clear();

                                // ToDo: Customize

                                // Update the cached values
                                //previousSpecEValue = searchResult.SpecEValue;
                                //previousQValue = searchResult.QValue;

                                // Append a new entry to the SortedSets
                                peptidesFoundForSpecEValue.Add(key);
                                peptidesFoundForQValue.Add(key);

                                firstMatchForGroup = true;
                            }

                            var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup, modInfo);
                            if (!modsAdded)
                            {
                                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                {
                                    errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'\n";
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

                        if (string.IsNullOrWhiteSpace(modificationSummaryFilePath))
                        {
                            ReportWarning("ParseMaxQuantSynopsisFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
                        }
                        else
                        {
                            modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                            SaveModificationSummaryFile(modificationSummaryFilePath);
                        }
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
        /// Parse a MaxQuant results line while creating the MaxQuant synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="columnMapping"></param>
        /// <param name="modInfo"></param>
        /// <param name="rowNumber">Row number (used for error reporting)</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantResultsFileEntry(
            string lineIn,
            ref MaxQuantSearchResult udtSearchResult,
            ref string errorLog,
            IDictionary<MaxQuantResultsFileColumns, int> columnMapping,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modInfo,
            int rowNumber)
        {
            try
            {
                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 11)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Scan], out udtSearchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid in row " + rowNumber, true);
                }

                if (!int.TryParse(udtSearchResult.Scan, out udtSearchResult.ScanNum))
                {
                    ReportError("Scan column is not numeric in row " + rowNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MaxQuant
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                // double.TryParse(udtSearchResult.CalculatedMonoMass, out var sequenceMonoMassMaxQuant);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrefixResidue], out udtSearchResult.PrefixResidue);

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Sequence], out udtSearchResult.Sequence))
                {
                    ReportError("Sequence column is missing or invalid in row " + rowNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.SuffixResidue], out udtSearchResult.SuffixResidue);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Modifications], out udtSearchResult.Modifications);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Composition], out udtSearchResult.Composition);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Protein], out udtSearchResult.Protein);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ProteinDesc], out udtSearchResult.ProteinDesc);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ProteinLength], out udtSearchResult.ProteinLength);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ResidueEnd], out udtSearchResult.ResidueEnd);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ResidueStart], out udtSearchResult.ResidueStart);

                // Parse the modification list to determine the total mod mass
                var totalModMass = ComputeTotalModMass(udtSearchResult.Modifications, modInfo);

                // Compute monoisotopic mass of the peptide
                udtSearchResult.CalculatedMonoMassPHRP = ComputePeptideMass(udtSearchResult.Sequence, totalModMass);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MostAbundantIsotopeMz], out udtSearchResult.MostAbundantIsotopeMz);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.NumMatchedFragments], out udtSearchResult.NumMatchedFragments);

                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.SpecEValue], out udtSearchResult.SpecEValue))
                {
                    double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.EValue], out udtSearchResult.EValue);

                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.QValue], out udtSearchResult.QValue))
                {
                    double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PepQValue], out udtSearchResult.PepQValue);

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the MassMaxQuant results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error parsing MassMaxQuant Results in ParseMaxQuantResultsFileEntry for Row " + rowNumber + "\n";
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the MaxQuant results file header line, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseMaxQuantResultsFileHeaderLine(string lineIn, IDictionary<MaxQuantResultsFileColumns, int> columnMapping)
        {
            // ToDo: Customize this

            // The expected column order from MaxQuant:
            //   Scan	Pre	Sequence	Post	Modifications	Composition	ProteinName	ProteinDesc	ProteinLength	Start	End	Charge	MostAbundantIsotopeMz	Mass	#MatchedFragments	Probability SpecEValue    EValue    QValue    PepQValue

            var columnNames = new SortedDictionary<string, MaxQuantResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Scan", MaxQuantResultsFileColumns.Scan},
                {"Pre", MaxQuantResultsFileColumns.PrefixResidue},
                {"Sequence", MaxQuantResultsFileColumns.Sequence},
                {"Post", MaxQuantResultsFileColumns.SuffixResidue},
                {"Modifications", MaxQuantResultsFileColumns.Modifications},
                {"Composition", MaxQuantResultsFileColumns.Composition},
                {"ProteinName", MaxQuantResultsFileColumns.Protein},
                {"ProteinDesc", MaxQuantResultsFileColumns.ProteinDesc},
                {"ProteinLength", MaxQuantResultsFileColumns.ProteinLength},
                {"Start", MaxQuantResultsFileColumns.ResidueStart},
                {"End", MaxQuantResultsFileColumns.ResidueEnd},
                {"Charge", MaxQuantResultsFileColumns.Charge},
                {"MostAbundantIsotopeMz", MaxQuantResultsFileColumns.MostAbundantIsotopeMz},
                {"Mass", MaxQuantResultsFileColumns.CalculatedMonoMass},
                {"MS1Features", MaxQuantResultsFileColumns.MS1Features},
                {"#MatchedFragments", MaxQuantResultsFileColumns.NumMatchedFragments},
                {"Probability", MaxQuantResultsFileColumns.Probability},
                {"SpecEValue", MaxQuantResultsFileColumns.SpecEValue},
                {"EValue", MaxQuantResultsFileColumns.EValue},
                {"QValue", MaxQuantResultsFileColumns.QValue},
                {"PepQValue", MaxQuantResultsFileColumns.PepQValue}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MaxQuantResultsFileColumns resultColumn in Enum.GetValues(typeof(MaxQuantResultsFileColumns)))
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
                        for (var index = 0; index <= splitLine.Length - 1; index++)
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
                                OnWarningEvent("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMaxQuantResultsFileHeaderLine");
                            }
                        }
                    }
                }

                if (useDefaultHeaders)
                {
                    // Use default column mappings
                    foreach (MaxQuantResultsFileColumns resultColumn in Enum.GetValues(typeof(MaxQuantResultsFileColumns)))
                    {
                        columnMapping[resultColumn] = (int)resultColumn;
                    }

                    // This is not a header line; return false
                    return false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MaxQuant results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        /// <summary>
        /// Parse the header line of a MaxQuant _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantSynFileHeaderLine(string lineIn, IDictionary<MaxQuantSynFileColumns, int> columnMapping)
        {
            var columnNames = MaxQuantSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MaxQuantSynFileColumns resultColumn in Enum.GetValues(typeof(MaxQuantSynFileColumns)))
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
                SetErrorMessage("Error parsing header in MaxQuant synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a MaxQuant Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequence"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantSynFileEntry(
            string lineIn,
            MaxQuantResults searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<MaxQuantSynFileColumns, int> columnMapping,
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

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MaxQuant Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Scan], out string scan);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Charge], out string charge);

                //searchResult.Scan = scan;
                //searchResult.Charge = charge;

                //if (!GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Sequence], out peptideSequence))
                //{
                //    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                //    {
                //        errorLog += "Error reading Peptide sequence value from MaxQuant Results line " +
                //                    (resultsProcessed + 1) + "\n";
                //    }
                //    return false;
                //}

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Protein], out string proteinName);
                //searchResult.ProteinName = proteinName;
                //searchResult.MultipleProteinCount = "0";

                //searchResult.PeptideDeltaMass = "0";

                //// Note that MaxQuant sequences don't actually have mod symbols; that information is tracked via searchResult.Modifications

                //// Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                //searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                //var searchResultBase = (SearchResultsBaseClass)searchResult;

                //ComputePseudoPeptideLocInProtein(searchResultBase);

                //// Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                //searchResult.ComputePeptideCleavageStateInProtein();

                //// Read the remaining data values
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MostAbundantIsotopeMz], out string mostAbundantIsotopeMz);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Modifications], out string modifications);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Composition], out string composition);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ProteinDesc], out string proteinDesc);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ProteinLength], out string proteinLength);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ResidueStart], out string residueStart);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ResidueEnd], out string residueEnd);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MatchedFragments], out string matchedFragments);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.SpecEValue], out string specEValue);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.EValue], out string eValue);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.QValue], out string qValue);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PepQValue], out string pepQValue);

                //searchResult.MostAbundantIsotopeMz = mostAbundantIsotopeMz;
                //searchResult.Modifications = modifications;
                //searchResult.Composition = composition;
                //searchResult.ProteinDesc = proteinDesc;
                //searchResult.ProteinLength = proteinLength;
                //searchResult.ResidueStart = residueStart;
                //searchResult.ResidueEnd = residueEnd;
                //searchResult.MatchedFragments = matchedFragments;
                //searchResult.SpecEValue = specEValue;
                //searchResult.EValue = eValue;
                //searchResult.QValue = qValue;
                //searchResult.PepQValue = pepQValue;

                //// Update the peptide location in the protein
                //if (!string.IsNullOrEmpty(searchResult.ResidueStart))
                //{
                //    int.TryParse(searchResult.ResidueStart, out var peptideLocInProteinStart);
                //    searchResult.PeptideLocInProteinStart = peptideLocInProteinStart;
                //}

                //if (!string.IsNullOrEmpty(searchResult.ResidueEnd))
                //{
                //    int.TryParse(searchResult.ResidueEnd, out var peptideLocInProteinEnd);
                //    searchResult.PeptideLocInProteinEnd = peptideLocInProteinEnd;
                //}

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorLog += "Error parsing MaxQuant Results for RowIndex '" + splitLine[0] + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MaxQuant Results in ParseMaxQuantSynFileEntry\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MaxQuant results file (Dataset_IcTda.tsv)</param>
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

                    // Load the MaxQuant Parameter File so that we can determine the modification names and masses
                    var modInfoExtracted = ExtractModInfoFromParamFile(SearchToolParameterFilePath, out var modInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Re-parse the MaxQuant parameter file to look for NumMatchesPerSpec
                    var numMatchesPerSpec = GetNumMatchesPerSpectrumToReport(SearchToolParameterFilePath);
                    if (numMatchesPerSpec > 1)
                    {
                        // Auto-change IgnorePeptideToProteinMapperErrors to True
                        // since the MaxQuant parameter file has NumMatchesPerSpec of 2 or higher
                        // (which results in PSMs associated with decoy proteins)
                        IgnorePeptideToProteinMapperErrors = true;

                        OnDebugEvent(string.Format(
                            "Set IgnorePeptideToProteinMapperErrors to true since NumMatchesPerSpec is {0} in the MaxQuant parameter file",
                            numMatchesPerSpec));
                    }

                    var baseName = "Dataset_maxq";

                    // Do not create a first-hits file for MaxQuant results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form Dataset_maxq_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath, modInfo);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to create the other PHRP files
                    success = ParseMaxQuantSynopsisFile(synOutputFilePath, outputDirectoryPath, false, modInfo);

                    if (success && CreateProteinModsFile)
                    {
                        // Check for an empty synopsis file
                        if (!ValidateFileHasData(synOutputFilePath, "Synopsis file", out var errorMessage))
                        {
                            ReportWarning(errorMessage);
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
                    SetErrorMessage("Error in MaxQuantResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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
            IEnumerable<MaxQuantSearchResult> filteredSearchResults,
            ref string errorLog)
        {
            // Sort filteredSearchResults by ascending SpecEValue, QValue, Scan, Peptide, and Protein
            var query = from item in filteredSearchResults orderby item.SpecEValueNum, item.QValueNum, item.ScanNum, item.Sequence, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, ref errorLog);
                index++;
            }
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan (typically there will only be one result for MaxQuant)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Output parameter: the actual filtered search results</param>
        private void StoreSynMatches(
            IList<MaxQuantSearchResult> searchResults,
            int startIndex,
            int endIndex,
            List<MaxQuantSearchResult> filteredSearchResults)
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
                if (searchResults[index].SpecEValueNum <= MSGFPlusSynopsisFileSpecEValueThreshold || searchResults[index].QValueNum < 0.1)
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
                var headerColumns = MaxQuantSynFileReader.GetColumnHeaderNamesAndIDs();

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
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            MaxQuantSearchResult udtSearchResult,
            ref string errorLog)
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

        private class MaxQuantSearchResultsComparerScanChargeScorePeptide : IComparer<MaxQuantSearchResult>
        {
            public int Compare(MaxQuantSearchResult x, MaxQuantSearchResult y)
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
