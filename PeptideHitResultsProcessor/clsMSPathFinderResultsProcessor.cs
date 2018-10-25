// This class reads in a MSPathFinder results file (_IcTda.tsv) and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 05/01/2015
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsMSPathFinderResultsProcessor : clsPHRPBaseClass
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        public clsMSPathFinderResultsProcessor()
        {
            mFileDate = "October 15, 2017";

            mGetModName = new Regex(@"(?<ModName>.+) (?<ResidueNumber>\d+)", RegexOptions.Compiled);
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_MSPathFinder_FILE = "_IcTda";

        public const string N_TERMINUS_SYMBOL_MSPATHFINDER = "-";
        public const string C_TERMINUS_SYMBOL_MSPATHFINDER = "-";

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        // These columns correspond to the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
        protected const int MSPathFinderResultsFileColCount = 21;
        public enum eMSPathFinderResultsFileColumns
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

        // These columns correspond to the Synopsis file created by this class
        protected const int MSPathFinderSynFileColCount = 18;
        public enum eMSPathFinderSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Charge = 2,
            MostAbundantIsotopeMz = 3,
            Mass = 4,
            Sequence = 5,                // PrefixLetter.Sequence.SuffixLetter
            Modifications = 6,
            Composition = 7,
            Protein = 8,
            ProteinDesc = 9,
            ProteinLength = 10,
            ResidueStart = 11,
            ResidueEnd = 12,
            MatchedFragments = 13,
            SpecEValue = 14,             // Column added 2015-08-25
            EValue = 15,                 // Column added 2015-08-25
            QValue = 16,
            PepQValue = 17
        }

        #endregion

        #region "Structures"
        // This data structure holds rows read from the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
        protected struct udtMSPathFinderSearchResultType
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

                // Unused at present: MH = String.Empty
                // Unused at present: DelM = String.Empty
                // Unused at present: DelM_PPM = String.Empty
            }
        }

        #endregion

        #region "Classwide Variables"

        #endregion
        private readonly Regex mGetModName;

        /// <summary>
        /// Step through the Modifications and associate each modification with the residues
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modInfo"></param>
        /// <remarks></remarks>
        private void AddModificationsToResidues(
            clsSearchResultsMSPathFinder searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            if (string.IsNullOrWhiteSpace(searchResult.Modifications))
            {
                return;
            }

            var mods = searchResult.Modifications.Split(',');
            var finalResidueLoc = searchResult.PeptideCleanSequence.Length;

            foreach (var modEntry in mods)
            {
                // Find the mod in the list of known modifications (loaded from the MSPathFinder parameter file)
                var matchFound = false;

                // Obtain the mod name, for example "Dehydro" from "Dehydro 52"
                var reMatch = mGetModName.Match(modEntry);

                if (!reMatch.Success)
                {
                    ReportError("Invalid MSPathFinder mod entry format; must be a name then a space then a number: " + modEntry);
                    continue;
                }

                var modName = reMatch.Groups["ModName"].Value;
                var residueNumber = reMatch.Groups["ResidueNumber"].Value;

                foreach (var modDef in modInfo)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.InvariantCultureIgnoreCase))
                    {
                        if (!int.TryParse(residueNumber, out var residueLocInPeptide))
                        {
                            ReportError("Mod entry does not have a number after the name: " + modEntry);
                            continue;
                        }

                        var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                        if (residueLocInPeptide <= 1)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                        }
                        else if (residueLocInPeptide >= finalResidueLoc)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
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
                        searchResult.SearchResultAddModification(modDef.ModMassVal, chMostRecentResidue, residueLocInPeptide, eResidueTerminusState, updateModOccurrenceCounts);

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
            clsSearchResultsMSPathFinder searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            bool success;

            try
            {
                // For other tools, we would add IsotopicMods here
                // This is not supported for MSPathFinder
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

        protected double ComputePeptideMass(string cleanSequence, double totalModMass)
        {
            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);
            mass += totalModMass;

            return mass;
        }

        /// <summary>
        /// Computes the total of all modifications defined for the sequence
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="modificationList">Comma separated list of modifications, e.g. Dehydro 52,Dehydro 63</param>
        /// <param name="modInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected double ComputeTotalModMass(
            string cleanSequence,
            string modificationList,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
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
                var reMatch = mGetModName.Match(modEntry);

                if (!reMatch.Success)
                {
                    ReportError("Mod entry does not have a name separated by a number: " + modEntry, true);
                }

                var modName = reMatch.Groups["ModName"].Value;

                foreach (var modDef in modInfo)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.InvariantCultureIgnoreCase))
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
            string outputFolderPath)
        {
            bool success;

            // Create the MTSPepToProteinMap file

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseName, outputFolderPath, mts: true);

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
                    success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                    if (!success)
                    {
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (success)
            {
                // If necessary, copy various PHRPReader support files to the output folder
                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.DirectoryName, Path.GetFileName(synOutputFilePath)), outputFolderPath);

                // Create the Protein Mods file
                success = CreateProteinModDetailsFile(synOutputFilePath, outputFolderPath, mtsPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSPathFinder);
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
        /// <param name="modInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            try
            {
                int[] columnMapping = null;
                var errorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        var headerParsed = false;
                        var rowNumber = 0;

                        // Initialize the list that will hold all of the records in the MSPathFinder result file
                        var searchResultsUnfiltered = new List<udtMSPathFinderSearchResultType>();

                        // Initialize the list that will hold all of the records that will ultimately be written out to disk
                        var filteredSearchResults = new List<udtMSPathFinderSearchResultType>();

                        // Parse the input file
                        while (!reader.EndOfStream & !AbortProcessing)
                        {
                            var lineIn = reader.ReadLine();
                            rowNumber += 1;

                            if (string.IsNullOrWhiteSpace(lineIn))
                            {
                                continue;
                            }

                            if (!headerParsed)
                            {
                                // Parse the header line

                                var success = ParseMSPathFinderResultsFileHeaderLine(lineIn, out columnMapping);
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

                            var udtSearchResult = new udtMSPathFinderSearchResultType();

                            var validSearchResult = ParseMSPathFinderResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping, modInfo, rowNumber);

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

                        // Sort the SearchResults by scan, charge, and ascending SpecEValue
                        searchResultsUnfiltered.Sort(new MSPathFinderSearchResultsComparerScanChargeScorePeptide());

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

        private bool ExtractModInfoFromParamFile(string mSGFDBParamFilePath,
            out List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            // The DMS-based parameter file for MSPathFinder uses the same formatting as MSGF+

            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MSPathFinder");
            RegisterEvents(modFileProcessor);

            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            // Note that this call will initialize modInfo
            var success = modFileProcessor.ExtractModInfoFromParamFile(mSGFDBParamFilePath, out modInfo);

            if (!success || mErrorCode != ePHRPErrorCodes.NoError)
            {
                if (mErrorCode == ePHRPErrorCodes.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MSPathFinder parameter file");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFDBModsWithModDefinitions(modInfo, mPeptideMods);

            return true;
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFolderPath"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <param name="modInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool ParseMSPathfinderSynopsisFile(
            string inputFilePath,
            string outputFolderPath,
            bool resetMassCorrectionTagsAndModificationDefinitions,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that ParseMSPathfinderSynopsisFile synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            int[] columnMapping = null;
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
                var searchResult = new clsSearchResultsMSPathFinder(mPeptideMods, mPeptideSeqMassCalculator);

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
                        var baseOutputFilePath = Path.Combine(outputFolderPath, Path.GetFileName(inputFilePath));
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
                                success = ParseMSPathFinderSynFileHeaderLine(lineIn, out columnMapping);
                                if (!success)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMSPathFinderSynFileEntry(lineIn, searchResult, ref errorLog,
                                                                                     resultsProcessed, columnMapping,
                                                                                     out _);

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

                            success = AddModificationsAndComputeMass(searchResult, firstMatchForGroup, modInfo);
                            if (!success)
                            {
                                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                {
                                    errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'" + "\n";
                                }
                            }

                            SaveResultsFileEntrySeqInfo(searchResult, firstMatchForGroup);

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
                        modificationSummaryFilePath = Path.Combine(outputFolderPath, modificationSummaryFilePath);

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

        /// <summary>
        /// Parse a MSPathFinder results line while creating the MSPathFinder synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="columnMapping"></param>
        /// <param name="modInfo"></param>
        /// <param name="rowNumber">Row number (used for error reporting)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool ParseMSPathFinderResultsFileEntry(
            string lineIn,
            ref udtMSPathFinderSearchResultType udtSearchResult,
            ref string errorLog,
            IList<int> columnMapping,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo,
            int rowNumber)
        {
            // Parses an entry from the MSPathFinder results file

            double sequenceMonoMassMSPathFinder = 0;    // Theoretical peptide monoisotopic mass, including mods, as computed by MSPathFinder

            bool validSearchResult;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length >= 11)
                {
                    if (!GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.Scan], out udtSearchResult.Scan))
                    {
                        ReportError("Scan column is missing or invalid in row " + rowNumber, true);
                    }

                    if (!int.TryParse(udtSearchResult.Scan, out udtSearchResult.ScanNum))
                    {
                        ReportError("Scan column is not numeric in row " + rowNumber, true);
                    }

                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MSPathFinder
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);
                    double.TryParse(udtSearchResult.CalculatedMonoMass, out sequenceMonoMassMSPathFinder);

                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.PrefixResidue], out udtSearchResult.PrefixResidue);

                    if (!GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.Sequence], out udtSearchResult.Sequence))
                    {
                        ReportError("Sequence column is missing or invalid in row " + rowNumber, true);
                    }

                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.SuffixResidue], out udtSearchResult.SuffixResidue);

                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.Modifications], out udtSearchResult.Modifications);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.Composition], out udtSearchResult.Composition);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.Protein], out udtSearchResult.Protein);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.ProteinDesc], out udtSearchResult.ProteinDesc);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.ProteinLength], out udtSearchResult.ProteinLength);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.ResidueEnd], out udtSearchResult.ResidueEnd);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.ResidueStart], out udtSearchResult.ResidueStart);

                    // Parse the list of modified residues to determine the total mod mass
                    var totalModMass = ComputeTotalModMass(udtSearchResult.Sequence, udtSearchResult.Modifications, modInfo);

                    // Compute monoisotopic mass of the peptide
                    udtSearchResult.CalculatedMonoMassPHRP = ComputePeptideMass(udtSearchResult.Sequence, totalModMass);

                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.MostAbundantIsotopeMz], out udtSearchResult.MostAbundantIsotopeMz);
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.NumMatchedFragments], out udtSearchResult.NumMatchedFragments);

                    if (GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.SpecEValue], out udtSearchResult.SpecEValue))
                    {
                        double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum);
                    }
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.EValue], out udtSearchResult.EValue);

                    if (GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.QValue], out udtSearchResult.QValue))
                    {
                        double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum);
                    }
                    GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderResultsFileColumns.PepQValue], out udtSearchResult.PepQValue);

                    validSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MassMSPathFinder results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error parsing MassMSPathFinder Results in ParseMSPathFinderResultsFileEntry for Row " + rowNumber + "\n";
                }

                validSearchResult = false;
            }

            return validSearchResult;
        }

        /// <summary>
        /// Parse the MSPathFinder results file header line
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        /// <remarks></remarks>
        private bool ParseMSPathFinderResultsFileHeaderLine(string lineIn, out int[] columnMapping)
        {
            // Parse the header line

            // The expected column order from MassMSPathFinder:
            //   Scan	Pre	Sequence	Post	Modifications	Composition	ProteinName	ProteinDesc	ProteinLength	Start	End	Charge	MostAbundantIsotopeMz	Mass	#MatchedFragments	Probability SpecEValue    EValue    QValue    PepQValue

            var columnNames = new SortedDictionary<string, eMSPathFinderResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {"Scan", eMSPathFinderResultsFileColumns.Scan},
                {"Pre", eMSPathFinderResultsFileColumns.PrefixResidue},
                {"Sequence", eMSPathFinderResultsFileColumns.Sequence},
                {"Post", eMSPathFinderResultsFileColumns.SuffixResidue},
                {"Modifications", eMSPathFinderResultsFileColumns.Modifications},
                {"Composition", eMSPathFinderResultsFileColumns.Composition},
                {"ProteinName", eMSPathFinderResultsFileColumns.Protein},
                {"ProteinDesc", eMSPathFinderResultsFileColumns.ProteinDesc},
                {"ProteinLength", eMSPathFinderResultsFileColumns.ProteinLength},
                {"Start", eMSPathFinderResultsFileColumns.ResidueStart},
                {"End", eMSPathFinderResultsFileColumns.ResidueEnd},
                {"Charge", eMSPathFinderResultsFileColumns.Charge},
                {"MostAbundantIsotopeMz", eMSPathFinderResultsFileColumns.MostAbundantIsotopeMz},
                {"Mass", eMSPathFinderResultsFileColumns.CalculatedMonoMass},
                {"MS1Features", eMSPathFinderResultsFileColumns.MS1Features},
                {"#MatchedFragments", eMSPathFinderResultsFileColumns.NumMatchedFragments},
                {"Probability", eMSPathFinderResultsFileColumns.Probability},
                {"SpecEValue", eMSPathFinderResultsFileColumns.SpecEValue},
                {"EValue", eMSPathFinderResultsFileColumns.EValue},
                {"QValue", eMSPathFinderResultsFileColumns.QValue},
                {"PepQValue", eMSPathFinderResultsFileColumns.PepQValue}
            };

            columnMapping = new int[MSPathFinderResultsFileColCount];

            try
            {
                // Initialize each entry in columnMapping to -1
                for (var index = 0; index <= columnMapping.Length - 1; index++)
                {
                    columnMapping[index] = -1;
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
                            if (columnNames.TryGetValue(splitLine[index], out var eResultFileColumn))
                            {
                                // Recognized column name; update columnMapping
                                columnMapping[(int)eResultFileColumn] = index;
                                useDefaultHeaders = false;
                            }
                            else
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
                    for (var index = 0; index <= columnMapping.Length - 1; index++)
                    {
                        columnMapping[index] = index;
                    }

                    // This is not a header line; return false
                    return false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MSPathFinder results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        private bool ParseMSPathFinderSynFileHeaderLine(string lineIn, out int[] columnMapping)
        {
            // Parse the header line

            var columnNames = new SortedDictionary<string, eMSPathFinderSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {clsPHRPParserMSPathFinder.DATA_COLUMN_ResultID, eMSPathFinderSynFileColumns.ResultID},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Scan, eMSPathFinderSynFileColumns.Scan},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Charge, eMSPathFinderSynFileColumns.Charge},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_MostAbundantIsotopeMz, eMSPathFinderSynFileColumns.MostAbundantIsotopeMz},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Mass, eMSPathFinderSynFileColumns.Mass},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Sequence, eMSPathFinderSynFileColumns.Sequence},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Modifications, eMSPathFinderSynFileColumns.Modifications},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Composition, eMSPathFinderSynFileColumns.Composition},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_Protein, eMSPathFinderSynFileColumns.Protein},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinDesc, eMSPathFinderSynFileColumns.ProteinDesc},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinLength, eMSPathFinderSynFileColumns.ProteinLength},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueStart, eMSPathFinderSynFileColumns.ResidueStart},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueEnd, eMSPathFinderSynFileColumns.ResidueEnd},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_MatchedFragments, eMSPathFinderSynFileColumns.MatchedFragments},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_SpecEValue, eMSPathFinderSynFileColumns.SpecEValue},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_EValue, eMSPathFinderSynFileColumns.EValue},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_QValue, eMSPathFinderSynFileColumns.QValue},
                {clsPHRPParserMSPathFinder.DATA_COLUMN_PepQValue, eMSPathFinderSynFileColumns.PepQValue}
            };

            columnMapping = new int[MSPathFinderSynFileColCount];

            try
            {
                // Initialize each entry in columnMapping to -1
                for (var index = 0; index <= columnMapping.Length - 1; index++)
                {
                    columnMapping[index] = -1;
                }

                var splitLine = lineIn.Split('\t');
                for (var index = 0; index <= splitLine.Length - 1; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var eResultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[(int)eResultFileColumn] = index;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MSPathFinder synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMSPathFinderSynFileEntry(
            string lineIn,
            clsSearchResultsMSPathFinder searchResult,
            ref string errorLog,
            int resultsProcessed,
            IReadOnlyList<int> columnMapping,
            out string peptideSequence)
        {
            // Parses an entry from the MSPathFinder Synopsis file

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

                if (!GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MSPathFinder Results line " + (resultsProcessed + 1).ToString() + "\n";
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.Sequence], out peptideSequence))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from MSPathFinder Results line " + (resultsProcessed + 1).ToString() + "\n";
                    }
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.Protein], out string proteinName);
                searchResult.ProteinName = proteinName;
                searchResult.MultipleProteinCount = "0";

                searchResult.PeptideDeltaMass = "0";

                // Note that MSPathFinder sequences don't actually have mod symbols; that information is tracked via searchResult.Modifications

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.MostAbundantIsotopeMz], out string mostAbundantIsotopeMz);

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.Modifications], out string modifications);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.Composition], out string composition);

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.ProteinDesc], out string proteinDesc);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.ProteinLength], out string proteinLength);

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.ResidueStart], out string residueStart);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.ResidueEnd], out string residueEnd);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.MatchedFragments], out string matchedFragments);

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.SpecEValue], out string specEValue);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.EValue], out string eValue);

                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.QValue], out string qValue);
                GetColumnValue(splitLine, columnMapping[(int)eMSPathFinderSynFileColumns.PepQValue], out string pepQValue);

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
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing MSPathFinder Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MSPathFinder Results in ParseMSPathFinderSynFileEntry" + "\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MSPathFinder results file (Dataset_IcTda.tsv)</param>
        /// <param name="outputFolderPath">Output folder</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputFolderPath, string parameterFilePath)
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

                if (!CleanupFilePaths(ref inputFilePath, ref outputFolderPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);

                    // Load the MSPathFinder Parameter File so that we can determine the modification names and masses
                    // Note that this call will initialize modInfo
                    var modInfoExtracted = ExtractModInfoFromParamFile(SearchToolParameterFilePath, out var modInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_IcTda.tsv" with "_mspath"
                    if (baseName.EndsWith("_IcTda", StringComparison.InvariantCultureIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - "_IcTda".Length) + "_mspath";
                    }

                    // Do not create a first-hits file for MSPathFinder results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_mspath_syn.txt
                    var synOutputFilePath = Path.Combine(outputFolderPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath, modInfo);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to create the other PHRP files
                    success = ParseMSPathfinderSynopsisFile(synOutputFilePath, outputFolderPath, false, modInfo);

                    if (success && CreateProteinModsFile)
                    {
                        // Check for an empty synopsis file
                        if (!ValidateFileHasData(synOutputFilePath, "Synopsis file", out var errorMessage))
                        {
                            ReportWarning(errorMessage);
                        }
                        else
                        {
                            success = CreateProteinModsFileWork(baseName, inputFile, synOutputFilePath, outputFolderPath);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in clsMSPathFinderResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile (1):" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return success;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            List<udtMSPathFinderSearchResultType> filteredSearchResults,
            ref string errorLog)
        {
            // Sort udtFilteredSearchResults by ascending SpecEValue, QValue, Scan, Peptide, and Protein
            var query = from item in filteredSearchResults orderby item.SpecEValueNum, item.QValueNum, item.ScanNum, item.Sequence, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, ref errorLog);
                index += 1;
            }
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan (typically there will only be one result for MSPathFinder)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Output parmaeter: the actual filtered search results</param>
        /// <remarks></remarks>
        private void StoreSynMatches(
            IList<udtMSPathFinderSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            List<udtMSPathFinderSearchResultType> filteredSearchResults)
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
                if (searchResults[index].SpecEValueNum <= MSGFDBSynopsisFileSpecEValueThreshold || searchResults[index].QValueNum < 0.1)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ref string errorLog)
        {
            // Write out the header line for synopsis / first hits files
            try
            {
                var data = new List<string>
                {
                    clsPHRPParserMSPathFinder.DATA_COLUMN_ResultID,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Scan,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Charge,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_MostAbundantIsotopeMz,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Mass,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Sequence,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Modifications,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Composition,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_Protein,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinDesc,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_ProteinLength,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueStart,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_ResidueEnd,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_MatchedFragments,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_SpecEValue,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_EValue,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_QValue,
                    clsPHRPParserMSPathFinder.DATA_COLUMN_PepQValue
                };

                writer.WriteLine(CollapseList(data));
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
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtMSPathFinderSearchResultType udtSearchResult,
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

        protected class MSPathFinderSearchResultsComparerScanChargeScorePeptide : IComparer<udtMSPathFinderSearchResultType>
        {
            public int Compare(udtMSPathFinderSearchResultType x, udtMSPathFinderSearchResultType y)
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

                // SpecEValue is the same; check qvalue
                if (x.QValueNum > y.QValueNum)
                {
                    return 1;
                }

                if (x.QValueNum < y.QValueNum)
                {
                    return -1;
                }

                // SpecEValue is the same; check sequence
                var result = string.Compare(x.Sequence, y.Sequence, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }
                return result;
            }
        }

        #endregion
    }
}
