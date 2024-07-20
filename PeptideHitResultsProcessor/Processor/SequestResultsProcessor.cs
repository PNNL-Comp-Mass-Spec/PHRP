// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://www.pnnl.gov/integrative-omics
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
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads a SEQUEST synopsis or first hits file (tab-delimited representation
    /// of the out file data, created from STARSuite Extractor) and creates a new file
    /// containing columns for cleavage and terminus state, modification information, and
    /// the monoisotopic mass of each peptide.  The data in the new file is linked to the
    /// original file by the Row ID number in the original file.  The user must provide a
    /// modification definition file that specifies the dynamic and/or static modifications
    /// used in the search.
    /// </summary>
    public class SequestResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: fht, methylation, sequest

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="options">Options</param>
        public SequestResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "January 13, 2022";
        }

        /// <summary>
        /// SEQUEST tool name
        /// </summary>
        public const string TOOL_NAME = "SEQUEST";

        /// <summary>
        /// First hits file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_FIRST_HITS_FILE = "_fht";

        /// <summary>
        /// Synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYNOPSIS_FILE = "_syn";

        private const int SEQUEST_SYN_FILE_MIN_COL_COUNT = 5;
        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        private bool AddDynamicAndStaticResidueMods(
            SearchResultsBaseClass searchResult,
            bool updateModOccurrenceCounts)
        {
            // Step through .PeptideSequenceWithMods
            // For each residue, check if a static mod is defined that affects that residue
            // For each mod symbol, determine the modification and add to searchResult

            var mostRecentLetter = '-';
            var residueLocInPeptide = 0;

            // Assume success for now
            var success = true;

            var sequenceWithMods = searchResult.PeptideSequenceWithMods;

            for (var index = 0; index <= sequenceWithMods.Length - 1; index++)
            {
                var character = sequenceWithMods[index];

                if (StringUtilities.IsLetterAtoZ(character))
                {
                    mostRecentLetter = character;
                    residueLocInPeptide++;

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == ModificationDefinition.ResidueModificationType.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                            if (modificationDefinition.TargetResiduesContain(character))
                            {
                                // Match found; add this modification
                                var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocInPeptide);

                                success = searchResult.SearchResultAddModification(
                                    modificationDefinition, character, residueLocInPeptide,
                                    residueTerminusState, updateModOccurrenceCounts);

                                if (!success)
                                {
                                    // Error adding this static mod
                                    SetErrorCode(PHRPErrorCode.UnspecifiedError);
                                    mErrorMessage = "Error calling searchResult.SearchResultAddModification for peptide '" + sequenceWithMods + "': " + searchResult.ErrorMessage;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (StringUtilities.IsLetterAtoZ(mostRecentLetter))
                {
                    success = searchResult.SearchResultAddDynamicModification(character, mostRecentLetter, residueLocInPeptide, searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);

                    if (!success)
                    {
                        // Error adding this dynamic mod
                        SetErrorCode(PHRPErrorCode.UnspecifiedError);
                        mErrorMessage = "Error calling searchResult.SearchResultAddDynamicModification for peptide '" + sequenceWithMods + "': " + searchResult.ErrorMessage;
                        break;
                    }
                }
                else
                {
                    // We found a modification symbol but mostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }

            return success;
        }

        private bool AddModificationsAndComputeMass(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool success;

            try
            {
                // Assume success for now
                success = true;

                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts);

                // Parse .PeptideSequenceWithMods to determine the modified residues present
                if (!AddDynamicAndStaticResidueMods(searchResult, updateModOccurrenceCounts))
                {
                    success = false;
                }

                // Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
                // Since Sequest allows a terminal peptide residue to be modified twice, we'll allow that to happen,
                //  even though, biologically, that's typically not possible
                // However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus
                //  (where two COOH groups are present)
                searchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, updateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
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
                "_syn",
                "_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// Parse a synopsis or first hits file
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions">If true, reset the mass correction tags and modification definitions</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseSynopsisOrFirstHitsFile(string inputFilePath, string outputDirectoryPath, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Note that synopsis files are normally sorted on XCorr descending, with lines
            // duplicated when peptide search results are mapped to multiple proteins

            // In order to prevent duplicate entries from being made to the ResultToSeqMap file,
            // we will keep track of the scan, charge, and peptide information parsed for each unique XCorr encountered
            // (see peptidesFoundForXCorrLevel below)

            var columnMapping = new Dictionary<SequestSynopsisFileColumns, int>();

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
                var searchResult = new SequestResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize a SortedSet that will be used to avoid double-counting the same PSM in the same scan
                // This is required since a PSM with multiple proteins will be listed on multiple lines in the synopsis file
                // Values are PeptideSequenceWithMods_Scan_NumScans_Charge_MH

                var peptidesFoundForXCorrLevel = new SortedSet<string>();

                var previousXCorr = string.Empty;

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(Options);

                    var errorMessages = new List<string>();

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using var reader = new StreamReader(inputFilePath);

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

                        var dataLine = true;

                        if (!headerParsed)
                        {
                            var validHeader = ParseSequestSynFileHeaderLine(lineIn, columnMapping);

                            if (validHeader)
                            {
                                dataLine = false;
                            }
                            else
                            {
                                // Error parsing header; assume this is a data line
                                dataLine = true;
                            }
                            headerParsed = true;
                        }

                        bool validSearchResult;
                        if (dataLine)
                        {
                            validSearchResult = ParseSequestResultsFileEntry(lineIn, columnMapping, searchResult, errorMessages);
                        }
                        else
                        {
                            validSearchResult = false;
                        }

                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.NumScans + "_" +
                                  searchResult.Charge + "_" + searchResult.PeptideMH;

                        bool firstMatchForGroup;

                        if (searchResult.PeptideXCorr == previousXCorr)
                        {
                            // New result has the same XCorr as the previous results
                            // See if peptidesFoundForXCorrLevel contains the peptide, scan, charge, and MH

                            if (peptidesFoundForXCorrLevel.Contains(key))
                            {
                                firstMatchForGroup = false;
                            }
                            else
                            {
                                peptidesFoundForXCorrLevel.Add(key);
                                firstMatchForGroup = true;
                            }
                        }
                        else
                        {
                            // New XCorr
                            // Reset peptidesFoundForXCorrLevel
                            peptidesFoundForXCorrLevel.Clear();

                            // Update previousXCorr
                            previousXCorr = searchResult.PeptideXCorr;

                            // Append a new entry to peptidesFoundForXCorrLevel
                            peptidesFoundForXCorrLevel.Add(key);
                            firstMatchForGroup = true;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);

                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
                            errorMessages.Add(string.Format(
                                "Error adding modifications to sequence for ResultID '{0}'{1}",
                                searchResult.ResultID, string.IsNullOrEmpty(mErrorMessage) ? string.Empty : ": " + mErrorMessage));

                            if (!string.IsNullOrEmpty(mErrorMessage))
                            {
                                mErrorMessage = string.Empty;
                            }
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
                            OnWarningEvent("ParseSynopsisOrFirstHitsFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
                        }
                        else
                        {
                            SaveModificationSummaryFile(Path.Combine(outputDirectoryPath, modificationSummaryFilePath));
                        }
                    }

                    // Inform the user if any errors occurred
                    if (errorMessages.Count > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                        return false;
                    }

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in ParseSynopsisOrFirstHitsFile", ex);
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
                SetErrorMessage("Error creating the output file in ParseSynopsisOrFirstHitsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        private bool ParseSequestResultsFileEntry(
            string lineIn,
            IDictionary<SequestSynopsisFileColumns, int> columnMapping,
            SequestResults searchResult,
            ICollection<string> errorMessages)
        {
            string[] splitLine = null;

            try
            {
                // Reset searchResult
                searchResult.Clear();

                splitLine = lineIn.TrimEnd().Split('\t');
                if (splitLine.Length < SEQUEST_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.RowIndex], out int resultId))
                {
                    ReportError("RowIndex column is missing or invalid", true);
                }

                searchResult.ResultID = resultId;

                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.Scan], out string scan);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.NumScans], out string numScans);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.Charge], out string charge);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.PeptideMH], out string peptideMh);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.XCorr], out string peptideXCorr);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.DeltaCn], out string peptideDeltaCn);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.Sp], out string peptideSp);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.ProteinName], out string proteinName);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.MultipleProteinCount], out string multipleProteinCount);

                searchResult.Scan = scan;
                searchResult.NumScans = numScans;
                searchResult.Charge = charge;
                searchResult.PeptideMH = peptideMh;
                searchResult.PeptideXCorr = peptideXCorr;
                searchResult.PeptideDeltaCn = peptideDeltaCn;
                searchResult.PeptideSp = peptideSp;
                searchResult.ProteinName = proteinName;
                searchResult.MultipleProteinCount = multipleProteinCount;

                if (!DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.PeptideSequence], out string peptideSequenceWithMods))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                // Calling this method will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the cleavage state and terminus state
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Sequest only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.DeltaCn2], out string peptideDeltaCn2);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.RankSP], out string peptideRankSp);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.RankXC], out string peptideRankXc);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.DelM], out string peptideDeltaMass);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.XcRatio], out string peptideXcRatio);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.PassFilt], out string peptidePassFilt);           // Legacy/Unused
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.MScore], out string peptideMScore);               // Legacy/Unused
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.NTT], out string peptideNtt);

                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.IonsObserved], out string ionsObserved);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.IonsExpected], out string ionsExpected);
                DataUtilities.GetColumnValue(splitLine, columnMapping[SequestSynopsisFileColumns.DelMPPM], out string delMppm);

                searchResult.PeptideDeltaCn2 = peptideDeltaCn2;
                searchResult.PeptideRankSP = peptideRankSp;
                searchResult.PeptideRankXC = peptideRankXc;
                searchResult.PeptideDeltaMass = peptideDeltaMass;
                searchResult.PeptideXcRatio = peptideXcRatio;
                searchResult.PeptidePassFilt = peptidePassFilt;
                searchResult.PeptideMScore = peptideMScore;
                searchResult.PeptideNTT = peptideNtt;
                searchResult.IonsObserved = ionsObserved;
                searchResult.IonsExpected = ionsExpected;
                searchResult.DelMPPM = delMppm;

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format("Error parsing Sequest Results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing Sequest Results in ParseSequestResultsFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Main processing routine
        /// </summary>
        /// <param name="inputFilePath">SEQUEST Synopsis or First-hits file</param>
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

                // Set this to true since SEQUEST param files can have the same mod mass on different residues, and those residues may use different symbols
                mPeptideMods.ConsiderModSymbolWhenFindingIdenticalMods = true;

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

                // Obtain the full path to the input file
                var inputFile = new FileInfo(inputFilePath);

                try
                {
                    success = ParseSynopsisOrFirstHitsFile(inputFile.FullName, outputDirectoryPath, false);
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error calling ParseSynopsisOrFirstHitsFile" + ex.Message);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                    success = false;
                }

                if (!success)
                {
                    return false;
                }

                if (Options.CreateProteinModsFile)
                {
                    success = CreateProteinModsFileWork(inputFile, outputDirectoryPath);
                }

                if (success)
                {
                    OperationComplete();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in SequestResultsProcessor.ProcessFile", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        private bool CreateProteinModsFileWork(FileInfo inputFile, string outputDirectoryPath)
        {
            bool success;

            // First create the MTS PepToProteinMap file using inputFile
            // Will also look for the first hits or synopsis file and use that too if it is present

            var sourcePHRPDataFiles = new List<string>();

            string additionalFile;
            var inputFileBaseName = Path.GetFileNameWithoutExtension(inputFile.Name);

            sourcePHRPDataFiles.Add(inputFile.FullName);
            if (inputFileBaseName.EndsWith(FILENAME_SUFFIX_SYNOPSIS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                additionalFile = ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_FIRST_HITS_FILE);
                if (File.Exists(additionalFile))
                {
                    sourcePHRPDataFiles.Add(additionalFile);
                }
            }
            else if (inputFileBaseName.EndsWith(FILENAME_SUFFIX_FIRST_HITS_FILE, StringComparison.OrdinalIgnoreCase))
            {
                additionalFile = ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_SYNOPSIS_FILE);
                if (File.Exists(additionalFile))
                {
                    sourcePHRPDataFiles.Add(additionalFile);
                }
            }

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, outputDirectoryPath, mts: true);

            if (File.Exists(mtsPepToProteinMapFilePath) && Options.UseExistingMTSPepToProteinMapFile)
            {
                success = true;
            }
            else
            {
                // Mapping file does not exist
                success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath, false);

                if (!success)
                {
                    OnWarningEvent(WARNING_MESSAGE_SKIPPING_PROTEIN_MODS_FILE_CREATION + " since CreatePepToProteinMapFile returned False");
                }
            }

            if (success)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                ValidatePHRPReaderSupportFiles(inputFile.FullName, outputDirectoryPath);

                // Now create the Protein Mods file
                success = CreateProteinModDetailsFile(inputFile.FullName, outputDirectoryPath, mtsPepToProteinMapFilePath, PeptideHitResultTypes.Sequest);
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                return true;
            }

            return true;
        }

        /// <summary>
        /// Parse the header line of a SEQUEST _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn">Data line</param>
        /// <param name="columnMapping">Column mapping</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseSequestSynFileHeaderLine(string lineIn, IDictionary<SequestSynopsisFileColumns, int> columnMapping)
        {
            var columnNames = SequestSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (SequestSynopsisFileColumns resultColumn in Enum.GetValues(typeof(SequestSynopsisFileColumns)))
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
                SetErrorMessage("Error parsing the header line in the Sequest synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Override this method to display the name of each class
        /// </summary>
        public override string ToString()
        {
            return string.Format("{0} results processor", TOOL_NAME);
        }
    }
}
