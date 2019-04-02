// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
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
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class reads in a synopsis or first hits file (tab-delimited representation
    /// of the out file data, created from STARSuite Extractor) and creates a new file
    /// containing columns for cleavage and terminus state, modification information, and
    /// the monoisotopic mass of each peptide.  The data in the new file is linked to the
    /// original file by the Row ID number in the original file.  The user must provide a
    /// modification definition file that specifies the dynamic and/or static modifications
    /// used in the search.
    /// </summary>
    public class clsSequestResultsProcessor : clsPHRPBaseClass
    {
        public clsSequestResultsProcessor()
        {
            mFileDate = "March 28, 2019";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"
        public const string FILENAME_SUFFIX_FIRST_HITS_FILE = "_fht";
        public const string FILENAME_SUFFIX_SYNOPSIS_FILE = "_syn";

        private const int SEQUEST_SYN_FILE_MIN_COL_COUNT = 5;
        private const int MAX_ERROR_LOG_LENGTH = 4096;

        #endregion

        #region "Structures"
        #endregion

        #region "Classwide Variables"

        #endregion

        #region "Properties"
        #endregion

        private bool AddDynamicAndStaticResidueMods(
            clsSearchResultsBaseClass searchResult,
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
                var chChar = sequenceWithMods[index];

                if (IsLetterAtoZ(chChar))
                {
                    mostRecentLetter = chChar;
                    residueLocInPeptide += 1;

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                            if (modificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                success = searchResult.SearchResultAddModification(
                                    modificationDefinition, chChar, residueLocInPeptide,
                                    searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);

                                if (!success)
                                {
                                    // Error adding this static mod
                                    SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
                                    mErrorMessage = "Error calling searchResult.SearchResultAddModification for peptide '" + sequenceWithMods + "': " + searchResult.ErrorMessage;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(mostRecentLetter))
                {
                    success = searchResult.SearchResultAddDynamicModification(chChar, mostRecentLetter, residueLocInPeptide, searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);

                    if (!success)
                    {
                        // Error adding this dynamic mod
                        SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
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

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
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

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool MTS)
        {
            var pepToProteinMapFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
            if (pepToProteinMapFilePath.EndsWith("_syn", StringComparison.OrdinalIgnoreCase) ||
                pepToProteinMapFilePath.EndsWith("_fht", StringComparison.OrdinalIgnoreCase))
            {
                // Remove _syn or _fht
                pepToProteinMapFilePath = pepToProteinMapFilePath.Substring(0, pepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(pepToProteinMapFilePath, outputDirectoryPath, MTS);
        }

        private void InitializeLocalVariables()
        {
            // Reserved for future use
        }

        /// <summary>
        /// Parse a synopsis or first hits file
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <returns></returns>
        /// <remarks>Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function</remarks>
        protected bool ParseSynopsisOrFirstHitsFile(string inputFilePath, string outputDirectoryPath, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Note that synopsis files are normally sorted on XCorr descending, with lines
            //  duplicated when peptide search results are mapped to multiple proteins
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file,
            //  we will keep track of the scan, charge, and peptide information parsed for each unique XCorr encountered

            var resultsProcessed = 0;

            var columnMapping = new Dictionary<clsPHRPParserSequest.SequestSynopsisFileColumns, int>();

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
                var searchResult = new clsSearchResultsSequest(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForXCorrLevel
                var peptidesFoundForXCorrLevel = new SortedSet<string>();
                var previousXCorr = string.Empty;

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    var errorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var reader = new StreamReader(inputFilePath))
                    {
                        resultsProcessed = 0;
                        var headerParsed = false;

                        // Create the output files
                        var baseOutputFilePath = Path.Combine(outputDirectoryPath, Path.GetFileName(inputFilePath));
                        var success = InitializeSequenceOutputFiles(baseOutputFilePath);

                        // Parse the input file
                        while (!reader.EndOfStream & !AbortProcessing)
                        {
                            var lineIn = reader.ReadLine();
                            if (string.IsNullOrWhiteSpace(lineIn))
                            {
                                continue;
                            }

                            var dataLine = true;

                            if (!headerParsed)
                            {
                                success = ParseSequestSynFileHeaderLine(lineIn, columnMapping);
                                if (success)
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
                                validSearchResult = ParseSequestResultsFileEntry(lineIn, columnMapping, searchResult, ref errorLog);
                            }
                            else
                            {
                                validSearchResult = false;
                            }

                            if (validSearchResult)
                            {
                                var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.NumScans + "_" + searchResult.Charge + "_" + searchResult.PeptideMH;
                                bool firstMatchForGroup;

                                if (searchResult.PeptideXCorr == previousXCorr)
                                {
                                    // New result has the same XCorr as the previous results
                                    // See if htPeptidesFoundForXCorrLevel contains the peptide, scan, charge, and MH

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
                                    // Reset htPeptidesFoundForXCorrLevel
                                    peptidesFoundForXCorrLevel.Clear();

                                    // Update previousXCorr
                                    previousXCorr = searchResult.PeptideXCorr;

                                    // Append a new entry to htPeptidesFoundForXCorrLevel
                                    peptidesFoundForXCorrLevel.Add(key);
                                    firstMatchForGroup = true;
                                }

                                success = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                                if (!success)
                                {
                                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'";
                                        if (!string.IsNullOrEmpty(mErrorMessage))
                                        {
                                            errorLog += ": " + mErrorMessage;
                                            mErrorMessage = string.Empty;
                                        }
                                        errorLog += "\n";
                                    }
                                }
                                SaveResultsFileEntrySeqInfo(searchResult, firstMatchForGroup);
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

                        if (string.IsNullOrWhiteSpace(modificationSummaryFilePath))
                        {
                            ReportWarning("ParseSynopsisOrFirstHitsFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
                        }
                        else
                        {
                            SaveModificationSummaryFile(Path.Combine(outputDirectoryPath, modificationSummaryFilePath));
                        }

                    }

                    // Inform the user if any errors occurred
                    if (errorLog.Length > 0)
                    {
                        SetErrorMessage("Invalid Lines: " + "\n" + errorLog);
                        return false;
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

        private bool ParseSequestResultsFileEntry(
            string lineIn,
            IDictionary<clsPHRPParserSequest.SequestSynopsisFileColumns, int> columnMapping,
            clsSearchResultsSequest searchResult,
            ref string errorLog)
        {
            string[] splitLine = null;

            bool validSearchResult;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                // Reset searchResult
                searchResult.Clear();

                splitLine = lineIn.TrimEnd().Split('\t');
                if (splitLine.Length < SEQUEST_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.RowIndex], out int resultId))
                {
                    ReportError("RowIndex column is missing or invalid", true);
                }

                searchResult.ResultID = resultId;

                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.NumScans], out string numScans);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.Charge], out string charge);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.PeptideMH], out string peptideMh);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.XCorr], out string peptideXCorr);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.DeltaCn], out string peptideDeltaCn);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.Sp], out string peptideSp);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.ProteinName], out string proteinName);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.MultipleProteinCount], out string multipleProteinCount);

                searchResult.Scan = scan;
                searchResult.NumScans = numScans;
                searchResult.Charge = charge;
                searchResult.PeptideMH = peptideMh;
                searchResult.PeptideXCorr = peptideXCorr;
                searchResult.PeptideDeltaCn = peptideDeltaCn;
                searchResult.PeptideSp = peptideSp;
                searchResult.ProteinName = proteinName;
                searchResult.MultipleProteinCount = multipleProteinCount;

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.PeptideSequence], out string peptideSequenceWithMods))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Sequest only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.DeltaCn2], out string peptideDeltaCn2);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.RankSP], out string peptideRankSp);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.RankXC], out string peptideRankXc);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.DelM], out string peptideDeltaMass);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.XcRatio], out string peptideXcRatio);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.PassFilt], out string peptidePassFilt);           // Legacy/Unused
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.MScore], out string peptideMScore);               // Legacy/Unused
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.NTT], out string peptideNtt);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.IonsObserved], out string ionsObserved);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.IonsExpected], out string ionsExpected);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserSequest.SequestSynopsisFileColumns.DelMPPM], out string delMppm);

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

                validSearchResult = true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing Sequest Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing Sequest Results in ParseSequestResultsFileEntry" + "\n";
                    }
                }
                validSearchResult = false;
            }

            return validSearchResult;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">Sequest Synopsis or First-hits file</param>
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

                // Set this to true since Sequest param files can have the same mod mass on different residues, and those residues may use different symbols
                mPeptideMods.ConsiderModSymbolWhenFindingIdenticalMods = true;

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

                // Obtain the full path to the input file
                var inputFile = new FileInfo(inputFilePath);

                try
                {
                    success = ParseSynopsisOrFirstHitsFile(inputFile.FullName, outputDirectoryPath, false);
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error calling ParseSynopsisOrFirstHitsFile" + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    success = false;
                }

                if (success && CreateProteinModsFile)
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
                SetErrorMessage("Error in clsSequestResultsProcessor.ProcessFile:  " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
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

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, outputDirectoryPath, MTS: true);

            if (File.Exists(mtsPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
            {
                success = true;
            }
            else
            {
                // Mapping file does not exist
                success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                if (!success)
                {
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                }
            }

            if (success)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                ValidatePHRPReaderSupportFiles(inputFile.FullName, outputDirectoryPath);

                // Now create the Protein Mods file
                success = CreateProteinModDetailsFile(inputFile.FullName, outputDirectoryPath, mtsPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Sequest);
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                success = true;
            }

            return true;
        }

        private bool ParseSequestSynFileHeaderLine(
            string lineIn,
            IDictionary<clsPHRPParserSequest.SequestSynopsisFileColumns, int> columnMapping)
        {
            // Parse the header line

            var columnNames = clsPHRPParserSequest.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (clsPHRPParserSequest.SequestSynopsisFileColumns resultColumn in Enum.GetValues(typeof(clsPHRPParserSequest.SequestSynopsisFileColumns)))
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
                SetErrorMessage("Error parsing header in Sequest synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }
    }
}
