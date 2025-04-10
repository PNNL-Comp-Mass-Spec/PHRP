﻿// This class reads in an MODa results file (txt format) and creates
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
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads a MODa results file (_moda.id) and creates
    /// a tab-delimited text file with the data
    /// </summary>
    public class MODaResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: fht, methylation, mgf, moda, MODa, ModDefs, txt

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="options">Options</param>
        public MODaResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "January 14, 2022";

            InitializeLocalVariables();
        }

        /// <summary>
        /// MODa tool name
        /// </summary>
        public const string TOOL_NAME = "MODa";

        /// <summary>
        /// MODa results file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_MODA_FILE = "_moda.id";

        /// <summary>
        /// N-terminus symbol used by MODa
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_MODA = "-";

        /// <summary>
        /// C-terminus symbol used by MODa
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_MODA = "-";

        /// <summary>
        /// Default synopsis file probability threshold
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const float DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD = MODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        // Note that as of April 2014, all mod masses reported by MODa are simply integers, meaning matching a trailing period is not necessary
        private const string MODA_MOD_MASS_REGEX = "([+-][0-9.]+)";

        private const byte MODA_MASS_DIGITS_OF_PRECISION = 0;

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to the tab-delimited file (_moda.id.txt) created by Java file anal_moda.jar
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

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file (_moda.id.txt) created by Java file anal_moda.jar
        /// </summary>
        /// <remarks>
        /// These columns hold data that this class will use when creating the synopsis file
        /// </remarks>
        private struct MODaSearchResult
        {
            // ReSharper disable NotAccessedField.Local

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
            public string Probability;              // Larger values are better
            public double ProbabilityNum;           // Larger values are better
            public int RankProbability;
            public string Peptide;
            public string Protein;
            public string PeptidePosition;          // Protein start/stop residues of the peptide, e.g. 108~115
            public double QValue;                   // Computed by this class

            // ReSharper restore NotAccessedField.Local

            /// <summary>
            /// Reset stored values to empty strings and zeros
            /// </summary>
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
                QValue = 0;
            }

            /// <summary>
            /// Show scan, peptide, and probability
            /// </summary>
            public readonly override string ToString()
            {
                return string.Format("Scan {0}: {1}, Probability {2}", ScanNum, Peptide, Probability);
            }
        }

        private int mDeltaMassWarningCount;

        private Dictionary<int, int> mSpectrumIndexToScanMap;

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult">Search result</param>
        /// <param name="updateModOccurrenceCounts">When true, update mod occurrence counts</param>
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
                var character = sequence[index];

                if (StringUtilities.IsLetterAtoZ(character))
                {
                    if (parsingModMass)
                    {
                        // Associate the mod mass in modMassDigits with the previous residue
                        AssociateDynamicModWithResidue(searchResult, mostRecentResidue, residueLocInPeptide,
                            modMassDigits, updateModOccurrenceCounts);
                        parsingModMass = false;
                    }

                    mostRecentResidue = character;
                    residueLocInPeptide++;

                    // Look for static mods to associate with this residue
                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ResidueModificationType.StaticMod)
                            continue;

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
                else
                {
                    var isNumberChar = character is '+' or '-' || char.IsDigit(character);

                    if (parsingModMass)
                    {
                        if (isNumberChar || character == '.')
                        {
                            modMassDigits += character;
                        }
                    }
                    else if (isNumberChar)
                    {
                        // Mod Mass Start
                        modMassDigits = character.ToString();
                        parsingModMass = true;
                    }
                    // ReSharper disable once RedundantIfElseBlock
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
                // Since InSpecT allows a terminal peptide residue to be modified twice, we'll allow that to happen,
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
        /// Ranks each entry (assumes all the data is from the same scan)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index</param>
        /// <param name="endIndex">End index</param>
        private void AssignRankByScore(
            IList<MODaSearchResult> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            if (startIndex == endIndex)
            {
                // Only one result
                var currentResult = searchResults[startIndex];
                currentResult.RankProbability = 1;
                searchResults[startIndex] = currentResult;
                return;
            }

            // Duplicate a portion of searchResults so that we can sort by descending Probability

            var resultsSubset = new Dictionary<int, MODaSearchResult>();

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

            return ComputePeptideMassForCleanSequence(cleanSequence, totalModMass);
        }

        private static readonly Regex RegexModMassRegEx = new(MODA_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptide, out var primarySequence, out _, out _);

            // Parse the dynamic mods reported by MODa
            foreach (Match match in RegexModMassRegEx.Matches(primarySequence))
            {
                // We use .TrimEnd() because the matched mod mass will end in a period if this mod applies to the final residue in a peptide
                if (double.TryParse(match.Groups[1].Value.TrimEnd('.'), out var modMassFound))
                {
                    totalModMass += modMassFound;
                }
            }

            // Now look for static mods
            // First determine the index of the last residue in primarySequence
            var indexLastChar = primarySequence.Length;

            for (var index = primarySequence.Length - 1; index >= 0; index += -1)
            {
                if (StringUtilities.IsLetterAtoZ(primarySequence[index]))
                {
                    indexLastChar = index;
                    break;
                }
            }

            for (var index = 0; index <= primarySequence.Length - 1; index++)
            {
                var character = primarySequence[index];

                if (!StringUtilities.IsLetterAtoZ(character))
                    continue;

                // Look for static mods to associate with this residue
                for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                {
                    if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ResidueModificationType.StaticMod)
                        continue;

                    var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                    var matchFound = modificationDefinition.TargetResiduesContain(character);

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
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputFilePath">Output file path</param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath)
        {
            try
            {
                var columnMapping = new Dictionary<MODaResultsFileColumns, int>();
                var errorMessages = new List<string>();

                // Open the input file and parse it
                // Initialize the stream reader and the stream writer
                using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                using var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                var resultsProcessed = 0;
                mDeltaMassWarningCount = 0;

                // Initialize the list that will hold all the records in the MODa result file
                var searchResultsUnfiltered = new List<MODaSearchResult>();

                // Initialize the list that will hold all the records that will ultimately be written out to disk
                var filteredSearchResults = new List<MODaSearchResult>();

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
                        WriteSynFHTFileHeader(writer, errorMessages);
                    }

                    if (!skipLine)
                    {
                        var validSearchResult = ParseMODaResultsFileEntry(lineIn, out var udtSearchResult, errorMessages, columnMapping);

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
        /// Load the static mods defined in the MODa parameter file
        /// </summary>
        /// <param name="modaParamFilePath">MODa parameter file path</param>
        /// <param name="modList">List of modifications</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ExtractModInfoFromMODaParamFile(string modaParamFilePath, out List<ModificationDefinition> modList)
        {
            modList = new List<ModificationDefinition>();

            try
            {
                if (string.IsNullOrEmpty(modaParamFilePath))
                {
                    SetErrorMessage("MODa Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(modaParamFilePath))
                {
                    SetErrorMessage("MODa param file not found: " + modaParamFilePath);
                    return false;
                }

                // Read the contents of the parameter (or mods) file
                using var reader = new StreamReader(new FileStream(modaParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

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
                    var kvSetting = SynFileReaderBaseClass.ParseKeyValueSetting(dataLine, '=', "#");

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
                            residue = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        if (residue.Equals("CTerm", StringComparison.OrdinalIgnoreCase))
                            residue = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                        if (double.TryParse(value, out var modMass))
                        {
                            if (Math.Abs(modMass - 0) > float.Epsilon)
                            {
                                var massCorrectionTag = mPeptideMods.LookupMassCorrectionTagByMass(modMass);

                                var modDef = new ModificationDefinition(ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                    modMass,
                                    residue,
                                    ModificationDefinition.ResidueModificationType.StaticMod,
                                    massCorrectionTag);
                                modList.Add(modDef);
                            }
                        }
                    }
                }

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MODa parameter file (" + Path.GetFileName(modaParamFilePath) + ")", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                return false;
            }
        }

        private void InitializeLocalVariables()
        {
            mSpectrumIndexToScanMap = new Dictionary<int, int>();
        }

        // ReSharper disable once UnusedMethodReturnValue.Local
        private bool LoadMGFIndexToScanMapFile(FileInfo inputFile)
        {
            var indexToScanMapFilePath = string.Empty;

            if (inputFile.Directory == null)
            {
                OnWarningEvent("LoadMGFIndexToScanMapFile: Could not determine the parent directory of " + inputFile.FullName);
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
                    indexToScanMapFilePath = scanMapFiles[0].FullName;
                }
                else if (scanMapFiles.Count == 0)
                {
                    OnWarningEvent("Did not find a mgf_IndexToScanMap file for " + sourceFileDescription + " in directory " + inputFile.Directory.FullName + "; scan numbers will be 0 in the synopsis file");
                    return false;
                }
                else
                {
                    OnWarningEvent("Found more than one potential mgf_IndexToScanMap file for " + sourceFileDescription + " in directory " + inputFile.Directory.FullName + " scan numbers will be 0 in the synopsis file");
                    return false;
                }

                var sourceFile = new FileInfo(indexToScanMapFilePath);

                if (!sourceFile.Exists)
                {
                    OnWarningEvent("MGF Index to Scan Map file not found; scan numbers will be 0 in the synopsis file: " + indexToScanMapFilePath);
                    return false;
                }

                using var reader = new StreamReader(new FileStream(sourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

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
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MGF Index to Scan Map file (" + Path.GetFileName(indexToScanMapFilePath) + ")", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
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
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="pepToProteinMapping">List of peptide to protein map values</param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions">If true, reset the mass correction tags and modification definitions</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            List<PepToProteinMapping> pepToProteinMapping,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Note that MODa synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            // we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered
            // (see peptidesFoundForProbabilityLevel below)

            var columnMapping = new Dictionary<MODaSynFileColumns, int>();

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
                var searchResult = new MODaResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize a SortedSet that will be used to avoid double-counting the same PSM in the same scan
                // This is required since a PSM with multiple proteins will be listed on multiple lines in the synopsis file
                // Values are PeptideSequenceWithMods_Scan_Charge

                var peptidesFoundForProbabilityLevel = new SortedSet<string>();

                var previousProbability = string.Empty;

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
                            var validHeader = ParseMODaSynFileHeaderLine(lineIn, columnMapping);

                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }
                            headerParsed = true;
                            continue;
                        }

                        var validSearchResult = ParseMODaSynFileEntry(
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
                                } while (pepToProteinMapIndex < pepToProteinMapping.Count && currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
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
                    SetErrorMessage("Error reading input file in ParseMODaSynopsisFile", ex);
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
                SetErrorMessage("Error creating the output file in ParseMODaSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parses an entry from the MODa results file
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="udtSearchResult">Search result</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaResultsFileEntry(
            string lineIn,
            out MODaSearchResult udtSearchResult,
            ICollection<string> errorMessages,
            IDictionary<MODaResultsFileColumns, int> columnMapping)
        {
            udtSearchResult = new MODaSearchResult();

            var rowIndex = "?";

            try
            {
                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 11)
                {
                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.SpectrumIndex], out udtSearchResult.SpectrumIndex))
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
                    OnWarningEvent("Error, could not resolve spectrumIndex to Scan Number: " + spectrumIndex);
                }

                // Monoisotopic mass value of the observed precursor_mz
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.ObservedMonoMass], out udtSearchResult.Precursor_mass);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = (short)StringUtilities.CIntSafe(udtSearchResult.Charge, 0);

                // Observed m/z, converted to monoisotopic mass
                if (double.TryParse(udtSearchResult.Precursor_mass, out var precursorMonoMass))
                {
                    if (udtSearchResult.ChargeNum > 0)
                    {
                        var precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, udtSearchResult.ChargeNum);
                        udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                    }
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                // Theoretical peptide monoisotopic mass, including mods, as computed by MODa
                double.TryParse(udtSearchResult.CalculatedMonoMass, out var peptideMonoMassMODa);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.DeltaMass], out udtSearchResult.DeltaMass);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Score], out udtSearchResult.Score);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Probability], out udtSearchResult.Probability);
                if (!double.TryParse(udtSearchResult.Probability, out udtSearchResult.ProbabilityNum))
                    udtSearchResult.ProbabilityNum = 0;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.Protein], out udtSearchResult.Protein);
                udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaResultsFileColumns.PeptidePosition], out udtSearchResult.PeptidePosition);

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
                    udtSearchResult.DelM = StringUtilities.MassErrorToString(delM);

                    var peptideDeltaMassCorrectedPpm =
                        SearchResultsBaseClass.ComputeDelMCorrectedPPM(delM, precursorMonoMass, peptideMonoMassMODa, true);

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
            catch (Exception ex)
            {
                // Error parsing this row from the MODa results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing MODa results in ParseMODaResultsFileEntry for RowIndex '{0}': {1}", rowIndex, ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing MODa Results in ParseMODaResultsFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the header line from a MODa results file
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseMODaResultsFileHeaderLine(string lineIn, IDictionary<MODaResultsFileColumns, int> columnMapping)
        {
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

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the MODa results file", ex);
                return false;
            }
        }

        /// <summary>
        /// Parse the header line from a MODa _syn.txt file
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaSynFileHeaderLine(string lineIn, IDictionary<MODaSynFileColumns, int> columnMapping)
        {
            var columnNames = MODaSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MODaSynFileColumns resultColumn in Enum.GetValues(typeof(MODaSynFileColumns)))
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
                SetErrorMessage("Error parsing the header line in the MODa synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a MODa Synopsis file
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="searchResult">Search result</param>
        /// <param name="errorMessages">Error messages</param>
        /// <param name="resultsProcessed">Number of results loaded when this method is called</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <param name="peptideSequenceWithMods">Peptide sequence, with modification symbols</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMODaSynFileEntry(
            string lineIn,
            MODaResults searchResult,
            ICollection<string> errorMessages,
            int resultsProcessed,
            IDictionary<MODaSynFileColumns, int> columnMapping,
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

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.ResultID], out string value))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading ResultID value from MODa results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Scan], out string scan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading peptide sequence from MODa results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.DelM], out string moDaComputedDelM);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.DelM_PPM], out string moDaComputedDelMppm);

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

                // Calling this method will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = searchResult as SearchResultsBaseClass;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the cleavage state and terminus state
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since InSpecT only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Spectrum_Index], out string spectrumIndex);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.PrecursorMZ], out string precursorMz);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.MH], out string parentIonMh);

                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Score], out string moDaScore);
                DataUtilities.GetColumnValue(splitLine, columnMapping[MODaSynFileColumns.Probability], out string probability);

                searchResult.Spectrum_Index = spectrumIndex;
                searchResult.Precursor_mz = precursorMz;
                searchResult.ParentIonMH = parentIonMh;
                searchResult.MODaScore = moDaScore;
                searchResult.Probability = probability;

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
                            "Error parsing MODa results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing MODa Results in ParseMODaSynFileEntry: " + ex.Message);
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing routine
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

                    var modaParameterFilePath = ResolveFilePath(inputFile.DirectoryName, Options.SearchToolParameterFilePath);

                    // Load the MODa Parameter File to look for any static mods
                    var modInfoExtracted = ExtractModInfoFromMODaParamFile(modaParameterFilePath, out var modaModInfo);

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
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath);

                    if (!success)
                    {
                        return false;
                    }

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    success = ParseMODaSynopsisFile(synOutputFilePath, outputDirectoryPath, pepToProteinMapping, false);

                    if (!success)
                    {
                        return false;
                    }

                    // This step is not necessary:
                    // if (success)
                    //	success = AppendDelMPPMRefinedToSynFile(synOutputFilePath);

                    // Remove all items from pepToProteinMapping to reduce memory overhead
                    pepToProteinMapping.Clear();
                    pepToProteinMapping.TrimExcess();

                    if (Options.CreateProteinModsFile)
                    {
                        // Use a higher match error threshold since some peptides reported by MODa don't perfectly match the FASTA file
                        const int MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD = 50;

                        success = CreateProteinModsFileWork(
                            baseName, inputFile,
                            synOutputFilePath, outputDirectoryPath,
                            PeptideHitResultTypes.MODa,
                            MAXIMUM_ALLOWABLE_MATCH_ERROR_PERCENT_THRESHOLD);
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in MODaResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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

        private void ResolveMODaModsWithModDefinitions(IReadOnlyCollection<ModificationDefinition> modaModInfo)
        {
            if (modaModInfo == null)
            {
                return;
            }

            // Call .LookupModificationDefinitionByMass for each entry in modaModInfo
            foreach (var modDef in modaModInfo)
            {
                if (string.IsNullOrEmpty(modDef.TargetResidues))
                {
                    mPeptideMods.LookupModificationDefinitionByMassAndModType(modDef.ModificationMass, modDef.ModificationType, default, AminoAcidModInfo.ResidueTerminusState.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                }
                else
                {
                    foreach (var targetResidue in modDef.TargetResidues)
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(modDef.ModificationMass, modDef.ModificationType, targetResidue, AminoAcidModInfo.ResidueTerminusState.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<MODaSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            // Sort filteredSearchResults by descending probability, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.ProbabilityNum descending, item.ScanNum, item.ChargeNum, item.Peptide, item.Protein select item;

            var index = 1;

            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, errorMessages);
                index++;
            }
        }

        private void StoreSynMatches(
            IList<MODaSearchResult> searchResults,
            int startIndex,
            int endIndex,
            ICollection<MODaSearchResult> filteredSearchResults)
        {
            AssignRankByScore(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].ProbabilityNum >= Options.MODaMODPlusSynopsisFileProbabilityThreshold)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer">Writer</param>
        /// <param name="errorMessages">Error messages</param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ICollection<string> errorMessages)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = MODaSynFileReader.GetColumnHeaderNamesAndIDs();

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
        /// <param name="resultID">Result ID</param>
        /// <param name="writer">Writer</param>
        /// <param name="udtSearchResult">Search result</param>
        /// <param name="errorMessages">Error messages</param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            MODaSearchResult udtSearchResult,
            ICollection<string> errorMessages)
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
            return string.Format("{0} results processor", TOOL_NAME);
        }

        private class MODaSearchResultsComparerScanChargeProbabilityPeptide : IComparer<MODaSearchResult>
        {
            public int Compare(MODaSearchResult x, MODaSearchResult y)
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
    }
}
