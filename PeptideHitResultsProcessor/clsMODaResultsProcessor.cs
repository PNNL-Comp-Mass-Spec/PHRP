// This class reads in an MODa results file (txt format) and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 04/01/2014
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsMODaResultsProcessor : clsPHRPBaseClass
    {
        // Ignore Spelling: MODa

        public clsMODaResultsProcessor()
        {
            FileDate = "July 10, 2019";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string TOOL_NAME = "MODa";

        public const string FILENAME_SUFFIX_MODA_FILE = "_moda.id";

        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_MODA = "-";

        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_MODA = "-";

        // ReSharper disable once UnusedMember.Global
        public const float DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        // Note that as of April 2014, all mod masses reported by MODa are simply integers, meaning matching a trailing period is not necessary
        private const string MODA_MOD_MASS_REGEX = "([+-][0-9.]+)";

        private const byte MODA_MASS_DIGITS_OF_PRECISION = 0;

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to the tab-delimited file (_moda.id.txt) created by MODa's anal_moda.jar file
        /// </summary>
        private enum MODaResultsFileColumns
        {
            SpectrumFileName = 0,
            SpectrumIndex = 1,
            ObservedMonoMass = 2,
            Charge = 3,
            CalculatedMonoMass = 4,
            DeltaMass = 5,
            Score = 6,
            Probability = 7,
            Peptide = 8,
            Protein = 9,
            PeptidePosition = 10
        }

        #endregion

        #region "Structures"

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file (_moda.id.txt) created by MODa's anal_moda.jar file
        /// </summary>
        private struct udtMODaSearchResultType
        {
            public string SpectrumFileName;
            public string SpectrumIndex;
            public int ScanNum;                     // Determined by looking for SpectrumIndex in the _mgf_IndexToScanMap.txt file
            public string Precursor_mass;           // Uncharged monoisotopic mass value of the observed precursor_mz, reported as ObservedMonoMass by MODa
            public string PrecursorMZ;              // Computed from ObservedMonoMass
            public string Charge;
            public short ChargeNum;
            public string CalculatedMonoMass;       // Theoretical monoisotopic mass of the peptide (including mods), as computed by MODa
            public string DeltaMass;                // Computed by MODa
            public string MH;                       // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                     // Computed using Precursor_mass - CalculatedMonoMass
            public string DelM_PPM;                 // Computed using DelM and CalculatedMonoMass
            public string Score;
            public string Probability;              // Higher values are better
            public double ProbabilityNum;           // Higher values are better
            public int RankProbability;
            public string Peptide;
            public string Protein;
            public string PeptidePosition;          // Protein start/stop residues of the peptide, e.g. 108~115
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
                Probability = string.Empty;
                ProbabilityNum = 0;
                RankProbability = 0;
                Peptide = string.Empty;
                Protein = string.Empty;
                PeptidePosition = string.Empty;
                FDR = 0;
                QValue = 0;
            }
        }

        #endregion

        #region "Class wide Variables"

        private int mDeltaMassWarningCount;

        private Dictionary<int, int> mSpectrumIndexToScanMap;

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
                        AssociateDynamicModWithResidue(searchResult, mostRecentResidue, residueLocInPeptide,
                            modMassDigits, updateModOccurrenceCounts);
                        parsingModMass = false;
                    }

                    mostRecentResidue = chChar;
                    residueLocInPeptide++;

                    // Look for static mods to associate with this residue
                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) != clsModificationDefinition.ModificationTypeConstants.StaticMod)
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

        private void AssociateDynamicModWithResidue(
            clsSearchResultsBaseClass searchResult,
            char mostRecentResidue,
            int residueLocInPeptide,
            string modMassDigits,
            bool updateModOccurrenceCounts)
        {
            var residueForMod = mostRecentResidue;
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
                    MODA_MASS_DIGITS_OF_PRECISION, MODA_MASS_DIGITS_OF_PRECISION);

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
        /// Ranks each entry  assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        private void AssignRankAndDeltaNormValues(
            IList<udtMODaSearchResultType> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of searchResults so that we can sort by descending Probability

            var resultsSubset = new Dictionary<int, udtMODaSearchResultType>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsByProbability = (from item in resultsSubset orderby item.Value.ProbabilityNum descending select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByProbability)
            {
                var result = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = result.ProbabilityNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(result.ProbabilityNum - lastValue) > double.Epsilon)
                    {
                        lastValue = result.ProbabilityNum;
                        currentRank++;
                    }
                }

                result.RankProbability = currentRank;
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

        private static readonly Regex RegexModMassRegEx = new Regex(MODA_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptide, out var primarySequence, out _, out _);

            // Parse the dynamic mods reported by MODa
            foreach (Match reMatch in RegexModMassRegEx.Matches(primarySequence))
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
                    if (mPeptideMods.GetModificationTypeByIndex(modIndex) != clsModificationDefinition.ModificationTypeConstants.StaticMod)
                        continue;

                    var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                    var matchFound = modificationDefinition.TargetResiduesContain(chChar);

                    if (!matchFound && index == 0)
                    {
                        matchFound = modificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
                    }

                    if (!matchFound && index == indexLastChar)
                    {
                        matchFound = modificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);
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
                "_MODa_syn",
                "_MODa_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MODa
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
                var columnMapping = new Dictionary<MODaResultsFileColumns, int>();
                var errorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var resultsProcessed = 0;
                    mDeltaMassWarningCount = 0;

                    // Initialize the list that will hold all of the records in the MODa result file
                    var searchResultsUnfiltered = new List<udtMODaSearchResultType>();

                    // Initialize the list that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<udtMODaSearchResultType>();

                    // Parse the input file
                    while (!reader.EndOfStream && !AbortProcessing)
                    {
                        var lineIn = reader.ReadLine();
                        var skipLine = false;

                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        if (resultsProcessed == 0)
                        {
                            // The first line might be a header line
                            // However, as of April 2014, MODa id.txt files do not have a header line

                            skipLine = ParseMODaResultsFileHeaderLine(lineIn, columnMapping);

                            // Write the header line to the output file
                            WriteSynFHTFileHeader(writer, ref errorLog);
                        }

                        if (!skipLine)
                        {
                            var udtSearchResult = new udtMODaSearchResultType();

                            var validSearchResult = ParseMODaResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping);

                            if (validSearchResult)
                            {
                                searchResultsUnfiltered.Add(udtSearchResult);
                            }

                            // Update the progress
                            UpdateSynopsisFileCreationProgress(reader);
                        }

                        resultsProcessed++;
                    }

                    // Sort the SearchResults by scan, charge, and Descending probability
                    searchResultsUnfiltered.Sort(new MODaSearchResultsComparerScanChargeProbabilityPeptide());

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

                    // Sort the data in udtFilteredSearchResults then write out to disk
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
                SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Load the static mods defined in the MODa parameter file
        /// </summary>
        /// <param name="modaParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ExtractModInfoFromMODaParamFile(string modaParamFilePath, out List<clsModificationDefinition> modInfo)
        {
            var success = false;

            modInfo = new List<clsModificationDefinition>();

            try
            {
                if (string.IsNullOrEmpty(modaParamFilePath))
                {
                    SetErrorMessage("MODa Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(PHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(modaParamFilePath))
                {
                    SetErrorMessage("MODa param file not found: " + modaParamFilePath);
                }
                else
                {
                    // Read the contents of the parameter (or mods) file
                    using (var reader = new StreamReader(new FileStream(modaParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
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
                                continue;
                            }

                            // Split the line on the equals sign
                            var kvSetting = clsPHRPParser.ParseKeyValueSetting(dataLine, '=', "#");

                            if (string.IsNullOrEmpty(kvSetting.Key))
                            {
                                continue;
                            }

                            if (string.Equals(kvSetting.Key, "add", StringComparison.OrdinalIgnoreCase))
                            {
                                // ModA defines all of its static modifications with the ADD keyword
                                // Split the value at the comma and create a new setting entry with the residue name

                                var value = kvSetting.Value;
                                var commaIndex = value.IndexOf(',');

                                var residue = value.Substring(0, commaIndex).Trim();
                                value = value.Substring(commaIndex + 1).Trim();

                                // Replace NTerm or CTerm with < or >
                                if (residue.Equals("NTerm", StringComparison.OrdinalIgnoreCase))
                                    residue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                                if (residue.Equals("CTerm", StringComparison.OrdinalIgnoreCase))
                                    residue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                                if (double.TryParse(value, out var modMass))
                                {
                                    if (Math.Abs(modMass - 0) > float.Epsilon)
                                    {
                                        var massCorrectionTag = mPeptideMods.LookupMassCorrectionTagByMass(modMass);

                                        var modDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                                                                   modMass,
                                                                                   residue,
                                                                                   clsModificationDefinition.ModificationTypeConstants.StaticMod,
                                                                                   massCorrectionTag);
                                        modInfo.Add(modDef);
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
                SetErrorMessage("Error reading the MODa parameter file (" + Path.GetFileName(modaParamFilePath) + "): " + ex.Message);
                SetErrorCode(PHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                success = false;
            }

            return success;
        }

        private void InitializeLocalVariables()
        {
            mSpectrumIndexToScanMap = new Dictionary<int, int>();
        }

        private bool LoadMGFIndexToScanMapFile(FileInfo inputFile)
        {
            var indexToScanMapFilePath = string.Empty;

            if (inputFile.Directory == null)
            {
                ReportWarning("LoadMGFIndexToScanMapFile: Could not determine the parent directory of " + inputFile.FullName);
                return false;
            }

            try
            {
                mSpectrumIndexToScanMap.Clear();

                // Look for the IndexToScanMap file that corresponds to inputFile
                List<FileInfo> scanMapFiles;
                var matchIndex = inputFile.Name.LastIndexOf("_moda", StringComparison.OrdinalIgnoreCase);
                string sourceFileDescription;

                if (matchIndex > 0)
                {
                    var datasetName = inputFile.Name.Substring(0, matchIndex);
                    scanMapFiles = inputFile.Directory.GetFiles(datasetName + "*mgf_IndexToScanMap*").ToList();
                    sourceFileDescription = " dataset " + datasetName;
                }
                else
                {
                    // Source file does not have "_moda" in the name
                    // Look for any mgf_IndexToScanMap file
                    scanMapFiles = inputFile.Directory.GetFiles("*mgf_IndexToScanMap*").ToList();
                    sourceFileDescription = inputFile.Name;
                }

                if (scanMapFiles.Count == 1)
                {
                    indexToScanMapFilePath = scanMapFiles.First().FullName;
                }
                else if (scanMapFiles.Count == 0)
                {
                    ReportWarning("Did not find a mgf_IndexToScanMap file for " + sourceFileDescription + " in directory " + inputFile.Directory.FullName + "; scan numbers will be 0 in the synopsis file");
                    return false;
                }
                else
                {
                    ReportWarning("Found more than one potential mgf_IndexToScanMap file for " + sourceFileDescription + " in directory " + inputFile.Directory.FullName + " scan numbers will be 0 in the synopsis file");
                    return false;
                }

                var sourceFile = new FileInfo(indexToScanMapFilePath);

                if (!sourceFile.Exists)
                {
                    ReportWarning("MGF Index to Scan Map file not found; scan numbers will be 0 in the synopsis file: " + indexToScanMapFilePath);
                    return false;
                }

                using (var reader = new StreamReader(new FileStream(sourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        var splitLine = lineIn.Split('\t');
                        if (splitLine.Length >= 3)
                        {
                            if (int.TryParse(splitLine[0], out var spectrumIndex))
                            {
                                if (int.TryParse(splitLine[1], out var scanStart))
                                {
                                    // int.TryParse(splitLine[2], out scanEnd);

                                    mSpectrumIndexToScanMap.Add(spectrumIndex, scanStart);
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MGF Index to Scan Map file (" + Path.GetFileName(indexToScanMapFilePath) + "): " + ex.Message);
                SetErrorCode(PHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }

            return true;
        }

        private int LookupScanBySpectrumIndex(int spectrumIndex)
        {
            if (mSpectrumIndexToScanMap.TryGetValue(spectrumIndex, out var scanNumber))
            {
                return scanNumber;
            }

            return 0;
        }

        /// <summary>
        /// Parse a MODa synopsis file
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function</remarks>
        private bool ParseMODaSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            List<udtPepToProteinMappingType> pepToProteinMapping,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Note that MODa synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            var columnMapping = new Dictionary<clsPHRPParserMODa.MODaSynFileColumns, int>();

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
                var searchResult = new clsSearchResultsMODa(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize peptidesFoundForProbabilityLevel
                var peptidesFoundForProbabilityLevel = new SortedSet<string>();
                var previousProbability = string.Empty;

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
                                var validHeader = ParseMODaSynFileHeaderLine(lineIn, columnMapping);
                                if (!validHeader)
                                {
                                    // Error parsing header
                                    SetErrorCode(PHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMODaSynFileEntry(lineIn,
                                                                          searchResult,
                                                                          ref errorLog,
                                                                          resultsProcessed,
                                                                          columnMapping,
                                                                          out var currentPeptideWithMods);

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
                                    } while (pepToProteinMapIndex < pepToProteinMapping.Count && currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
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

        /// <summary>
        /// Parses an entry from the MODa results file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaResultsFileEntry(
            string lineIn,
            ref udtMODaSearchResultType udtSearchResult,
            ref string errorLog,
            IDictionary<MODaResultsFileColumns, int> columnMapping)
        {
            var rowIndex = "?";

            try
            {
                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 11)
                {
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);

                if (!GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.SpectrumIndex], out udtSearchResult.SpectrumIndex))
                {
                    ReportError("Index column is missing or invalid", true);
                }
                else
                {
                    rowIndex = udtSearchResult.SpectrumIndex;
                }

                if (!int.TryParse(udtSearchResult.SpectrumIndex, out var spectrumIndex))
                {
                    ReportError("Index column is not numeric", true);
                }

                udtSearchResult.ScanNum = LookupScanBySpectrumIndex(spectrumIndex);
                if (udtSearchResult.ScanNum == 0)
                {
                    ReportWarning("Error, could not resolve spectrumIndex to Scan Number: " + spectrumIndex);
                }

                // Monoisotopic mass value of the observed precursor_mz
                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.ObservedMonoMass], out udtSearchResult.Precursor_mass);
                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Observed m/z, converted to monoisotopic mass
                if (double.TryParse(udtSearchResult.Precursor_mass, out var precursorMonoMass))
                {
                    if (udtSearchResult.ChargeNum > 0)
                    {
                        var precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, udtSearchResult.ChargeNum);
                        udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                    }
                }

                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                // Theoretical peptide monoisotopic mass, including mods, as computed by MODa
                double.TryParse(udtSearchResult.CalculatedMonoMass, out var peptideMonoMassMODa);

                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.DeltaMass], out udtSearchResult.DeltaMass);
                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Score], out udtSearchResult.Score);

                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Probability], out udtSearchResult.Probability);
                if (!double.TryParse(udtSearchResult.Probability, out udtSearchResult.ProbabilityNum))
                    udtSearchResult.ProbabilityNum = 0;

                if (!GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Protein], out udtSearchResult.Protein);
                udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);
                GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.PeptidePosition], out udtSearchResult.PeptidePosition);

                // Parse the sequence to determine the total mod mass
                // Note that we do not remove any of the mod symbols since MODa identifies mods by mass alone
                // Note that static mods are implied (thus are not explicitly displayed by MODa)
                var totalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                // Compute monoisotopic mass of the peptide
                // Theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                var peptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, totalModMass);

                if (Math.Abs(peptideMonoMassMODa) < double.Epsilon)
                {
                    peptideMonoMassMODa = peptideMonoMassPHRP;
                }

                ValidateMatchingMonoisotopicMass(TOOL_NAME, udtSearchResult.Peptide, peptideMonoMassPHRP, peptideMonoMassMODa, ref mDeltaMassWarningCount);

                if (peptideMonoMassMODa > 0)
                {
                    // Compute DelM and DelM_PPM
                    var delM = precursorMonoMass - peptideMonoMassMODa;
                    udtSearchResult.DelM = MassErrorToString(delM);

                    var peptideDeltaMassCorrectedPpm =
                        clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(delM, precursorMonoMass, true, peptideMonoMassMODa);

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
                // Error parsing this row from the MODa results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        errorLog += "Error parsing MODa Results in ParseMODaResultsFileEntry for RowIndex '" + rowIndex + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MODa Results in ParseMODaResultsFileEntry\n";
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the header line from a MODa results file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseMODaResultsFileHeaderLine(string lineIn, IDictionary<MODaResultsFileColumns, int> columnMapping)
        {
            // The expected column order from MODa:
            //   SpectrumFile	Index	ObservedMonoMass	Charge	CalculatedMonoMass	DeltaMass	Score	Probability	Peptide	Protein	PeptidePosition

            var columnNames = new SortedDictionary<string, MODaResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"SpectrumFile", MODaResultsFileColumns.SpectrumFileName},
                {"Index", MODaResultsFileColumns.SpectrumIndex},
                {"ObservedMW", MODaResultsFileColumns.ObservedMonoMass},
                {"Charge", MODaResultsFileColumns.Charge},
                {"CalculatedMW", MODaResultsFileColumns.CalculatedMonoMass},
                {"DeltaMass", MODaResultsFileColumns.DeltaMass},
                {"Score", MODaResultsFileColumns.Score},
                {"Probability", MODaResultsFileColumns.Probability},
                {"Peptide", MODaResultsFileColumns.Peptide},
                {"Protein", MODaResultsFileColumns.Protein},
                {"PeptidePosition", MODaResultsFileColumns.PeptidePosition}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MODaResultsFileColumns resultColumn in Enum.GetValues(typeof(MODaResultsFileColumns)))
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
                            }
                            else
                            {
                                // Unrecognized column name
                                Console.WriteLine("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMODaResultsFileHeaderLine");
                            }
                        }
                    }
                }

                if (useDefaultHeaders)
                {
                    // Use default column mappings
                    foreach (MODaResultsFileColumns resultColumn in Enum.GetValues(typeof(MODaResultsFileColumns)))
                    {
                        columnMapping[resultColumn] = (int)resultColumn;
                    }

                    // This is not a header line; return false
                    return false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MODa results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        /// <summary>
        /// Parse the header line from a MODa _syn.txt file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaSynFileHeaderLine(string lineIn, IDictionary<clsPHRPParserMODa.MODaSynFileColumns, int> columnMapping)
        {
            var columnNames = clsPHRPParserMODa.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (clsPHRPParserMODa.MODaSynFileColumns resultColumn in Enum.GetValues(typeof(clsPHRPParserMODa.MODaSynFileColumns)))
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
                SetErrorMessage("Error parsing header in MODa synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a MODa Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequenceWithMods"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaSynFileEntry(
            string lineIn,
            clsSearchResultsMODa searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<clsPHRPParserMODa.MODaSynFileColumns, int> columnMapping,
            out string peptideSequenceWithMods)
        {
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

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MODa Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from MODa Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.DelM], out string moDaComputedDelM);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.DelM_PPM], out string moDaComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.MODaComputedDelM = moDaComputedDelM;
                searchResult.MODaComputedDelMPPM = moDaComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.MODaComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MODa results file as "Observed_Mass - Theoretical_Mass"
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

                var searchResultBase = searchResult as clsSearchResultsBaseClass;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Spectrum_Index], out string spectrumIndex);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.MH], out string parentIonMh);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Score], out string moDaScore);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMODa.MODaSynFileColumns.Probability], out string probability);

                searchResult.Spectrum_Index = spectrumIndex;
                searchResult.Precursor_mz = precursorMz;
                searchResult.ParentIonMH = parentIonMh;
                searchResult.MODaScore = moDaScore;
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
                        errorLog += "Error parsing MODa Results for RowIndex '" + splitLine[0] + "'\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MODa Results in ParseMODaSynFileEntry\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MODa results file (Dataset_moda.id.txt)</param>
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

                    // Load the MODa Parameter File to look for any static mods
                    var modInfoExtracted = ExtractModInfoFromMODaParamFile(SearchToolParameterFilePath, out var modaModInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Resolve the mods in modaModInfo with the ModDefs mods
                    ResolveMODaModsWithModDefinitions(modaModInfo);

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_moda.id" with "_moda"
                    if (baseName.EndsWith(FILENAME_SUFFIX_MODA_FILE, StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - FILENAME_SUFFIX_MODA_FILE.Length) + "_moda";
                    }

                    // Load the MSG IndexToScanMap file (if it exists)
                    LoadMGFIndexToScanMapFile(inputFile);

                    // Do not create a first-hits file for MODa results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_moda_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    success = ParseMODaSynopsisFile(synOutputFilePath, outputDirectoryPath, pepToProteinMapping, false);

                    // This step is not necessary:
                    // if (success)
                    //	success = AppendDelMPPMRefinedToSynFile(synOutputFilePath);

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
                    SetErrorMessage("Error in clsMODaResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseName, outputDirectoryPath, mts: true);

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
                    // We do this since a small number of peptides reported by MODa don't perfectly match the fasta file
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
                                                          clsPHRPReader.PeptideHitResultTypes.MODa);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
            }

            return true;
        }

        private void ResolveMODaModsWithModDefinitions(IReadOnlyCollection<clsModificationDefinition> modaModInfo)
        {
            if (modaModInfo == null)
            {
                return;
            }

            // Call .LookupModificationDefinitionByMass for each entry in modaModInfo
            foreach (var modInfo in modaModInfo)
            {
                if (string.IsNullOrEmpty(modInfo.TargetResidues))
                {
                    mPeptideMods.LookupModificationDefinitionByMassAndModType(modInfo.ModificationMass, modInfo.ModificationType, default, clsAminoAcidModInfo.ResidueTerminusStateConstants.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                }
                else
                {
                    foreach (var targetResidue in modInfo.TargetResidues)
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(modInfo.ModificationMass, modInfo.ModificationType, targetResidue, clsAminoAcidModInfo.ResidueTerminusStateConstants.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<udtMODaSearchResultType> filteredSearchResults,
            ref string errorLog)
        {
            // Sort filteredSearchResults by descending probability, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.ProbabilityNum descending, item.ScanNum, item.ChargeNum, item.Peptide, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, ref errorLog);
                index++;
            }
        }

        private void StoreSynMatches(
            IList<udtMODaSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            ICollection<udtMODaSearchResultType> filteredSearchResults)
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
                var headerColumns = clsPHRPParserMODa.GetColumnHeaderNamesAndIDs();

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
        /// <param name="errorLog"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtMODaSearchResultType udtSearchResult,
            ref string errorLog)
        {
            try
            {
                // Primary Columns
                //
                // MODa
                // ResultID   Scan   Spectrum_Index   Charge   PrecursorMZ   DelM   DelM_PPM   MH   Peptide   Protein   Score   Probability   Rank_Probability   PeptidePosition      QValue

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
                    udtSearchResult.Protein,
                    udtSearchResult.Score,
                    udtSearchResult.Probability,
                    udtSearchResult.RankProbability.ToString(),
                    udtSearchResult.PeptidePosition,
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

        #region "IComparer Classes"

        private class MODaSearchResultsComparerScanChargeProbabilityPeptide : IComparer<udtMODaSearchResultType>
        {
            public int Compare(udtMODaSearchResultType x, udtMODaSearchResultType y)
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

                // Charge is the same; check ProbabilityNum
                if (x.ProbabilityNum < y.ProbabilityNum)
                {
                    return 1;
                }

                if (x.ProbabilityNum > y.ProbabilityNum)
                {
                    return -1;
                }

                // Probability is the same; check peptide
                var result = string.CompareOrdinal(x.Peptide, y.Peptide);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.CompareOrdinal(x.Protein, y.Protein);
                }
                return result;
            }
        }

        #endregion
    }
}
