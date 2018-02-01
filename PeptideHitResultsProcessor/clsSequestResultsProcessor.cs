// This class reads in a synopsis or first hits file (tab-delimited representation
// of the out file data, created from STARSuite Extractor) and creates a new file
// containing columns for cleavage and terminus state, modification information, and
// the monoisotopic mass of each peptide.  The data in the new file is linked to the
// original file by the Row ID number in the original file.  The user must provide a
// modification definition file that specifies the dynamic and/or static modifications
// used in the search.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsSequestResultsProcessor : clsPHRPBaseClass
    {
        public clsSequestResultsProcessor()
        {
            mFileDate = "October 13, 2017";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"
        public const string FILENAME_SUFFIX_FIRST_HITS_FILE = "_fht";
        public const string FILENAME_SUFFIX_SYNOPSIS_FILE = "_syn";

        private const int SEQUEST_SYN_FILE_MIN_COL_COUNT = 5;
        private const int MAX_ERROR_LOG_LENGTH = 4096;

        // These columns correspond to the tab-delimited file created directly by SEQUEST
        protected const int SequestSynopsisFileColCount = 27;
        public enum eSequestSynopsisFileColumns
        {
            RowIndex = 0,
            Scan = 1,
            NumScans = 2,
            Charge = 3,
            PeptideMH = 4,
            XCorr = 5,
            DeltaCn = 6,
            Sp = 7,
            ProteinName = 8,                 // Aka Reference
            MultipleProteinCount = 9,        // Aka MO = MultipleORFCount; this is 0 if the peptide is in just one protein; 1 if in 2 proteins, etc.
            PeptideSequence = 10,            // This is the sequence with prefix and suffix residues and also with modification symbols
            DeltaCn2 = 11,
            RankSP = 12,
            RankXC = 13,
            DelM = 14,
            XcRatio = 15,
            PassFilt = 16,                   // Legacy/unused
            MScore = 17,                     // Legacy/unused
            NTT = 18,                        // Number of tryptic terminii
            IonsObserved = 19,               // Added in August 2011
            IonsExpected = 20,               // Added in August 2011
            DelMPPM = 21,                    // Added in August 2011
            Cleavage_State = 22,             // This column and the ones after it are computed by this program and appended to the input file or saved in a new file
            Terminus_State = 23,
            Mod_Count = 24,
            Mod_Description = 25,
            Monoisotopic_Mass = 26
        }

        #endregion

        #region "Structures"
        #endregion

        #region "Classwide Variables"

        #endregion

        #region "Properties"
        #endregion

        private bool AddDynamicAndStaticResidueMods(
            clsSearchResultsBaseClass searchResult,
            bool blnUpdateModOccurrenceCounts)
        {
            // Step through .PeptideSequenceWithMods
            // For each residue, check if a static mod is defined that affects that residue
            // For each mod symbol, determine the modification and add to searchResult

            var chMostRecentLetter = '-';
            var intResidueLocInPeptide = 0;

            // Assume success for now
            var blnSuccess = true;

            var strSequence = searchResult.PeptideSequenceWithMods;
            for (var intIndex = 0; intIndex <= strSequence.Length - 1; intIndex++)
            {
                var chChar = strSequence[intIndex];

                if (IsLetterAtoZ(chChar))
                {
                    chMostRecentLetter = chChar;
                    intResidueLocInPeptide += 1;

                    for (var intModIndex = 0; intModIndex <= mPeptideMods.ModificationCount - 1; intModIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(intModIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);

                            if (objModificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                blnSuccess = searchResult.SearchResultAddModification(
                                    objModificationDefinition, chChar, intResidueLocInPeptide,
                                    searchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);

                                if (!blnSuccess)
                                {
                                    // Error adding this static mod
                                    SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
                                    mErrorMessage = "Error calling searchResult.SearchResultAddModification for peptide '" + strSequence + "': " + searchResult.ErrorMessage;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(chMostRecentLetter))
                {
                    blnSuccess = searchResult.SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, searchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);

                    if (!blnSuccess)
                    {
                        // Error adding this dynamic mod
                        SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
                        mErrorMessage = "Error calling searchResult.SearchResultAddDynamicModification for peptide '" + strSequence + "': " + searchResult.ErrorMessage;
                        break;
                    }
                }
                else
                {
                    // We found a modification symbol but chMostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }

            return blnSuccess;
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass searchResult, bool blnUpdateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool blnSuccess;

            try
            {
                // Assume success for now
                blnSuccess = true;

                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                searchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts);

                // Parse .PeptideSequenceWithMods to determine the modified residues present
                if (!AddDynamicAndStaticResidueMods(searchResult, blnUpdateModOccurrenceCounts))
                {
                    blnSuccess = false;
                }

                // Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
                // Since Sequest allows a terminal peptide residue to be modified twice, we'll allow that to happen,
                //  even though, biologically, that's typically not possible
                // However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus
                //  (where two COOH groups are present)
                searchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();
            }
            catch (Exception)
            {
                blnSuccess = false;
            }

            return blnSuccess;
        }

        protected override string ConstructPepToProteinMapFilePath(string strInputFilePath, string strOutputFolderPath, bool MTS)
        {
            var strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
            if (strPepToProteinMapFilePath.EndsWith("_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS);
        }

        private void InitializeLocalVariables()
        {
            // Reserved for future use
        }

        protected bool ParseSynopsisOrFirstHitsFile(string strInputFilePath, string strOutputFolderPath, bool blnResetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that synopsis files are normally sorted on XCorr descending, with lines
            //  duplicated when peptide search results are mapped to multiple proteins
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file,
            //  we will keep track of the scan, charge, and peptide information parsed for each unique XCorr encountered

            var intResultsProcessed = 0;

            int[] intColumnMapping = null;

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (blnResetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                // Initialize searchResult
                var searchResult = new clsSearchResultsSequest(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForXCorrLevel
                var htPeptidesFoundForXCorrLevel = new Hashtable();
                var strPreviousXCorr = string.Empty;

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    var strErrorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var srDataFile = new StreamReader(strInputFilePath))
                    {
                        intResultsProcessed = 0;
                        var blnHeaderParsed = false;

                        // Create the output files
                        var strBaseOutputFilePath = Path.Combine(strOutputFolderPath, Path.GetFileName(strInputFilePath));
                        var blnSuccess = InitializeSequenceOutputFiles(strBaseOutputFilePath);

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !AbortProcessing)
                        {
                            var strLineIn = srDataFile.ReadLine();
                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            var blnDataLine = true;

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseSequestSynFileHeaderLine(strLineIn, out intColumnMapping);
                                if (blnSuccess)
                                {
                                    blnDataLine = false;
                                }
                                else
                                {
                                    // Error parsing header; assume this is a data line
                                    blnDataLine = true;
                                }
                                blnHeaderParsed = true;
                            }

                            bool blnValidSearchResult;
                            if (blnDataLine)
                            {
                                blnValidSearchResult = ParseSequestResultsFileEntry(strLineIn, intColumnMapping, searchResult, ref strErrorLog);
                            }
                            else
                            {
                                blnValidSearchResult = false;
                            }

                            if (blnValidSearchResult)
                            {
                                var strKey = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.NumScans + "_" + searchResult.Charge + "_" + searchResult.PeptideMH;
                                bool blnFirstMatchForGroup;

                                if (searchResult.PeptideXCorr == strPreviousXCorr)
                                {
                                    // New result has the same XCorr as the previous results
                                    // See if htPeptidesFoundForXCorrLevel contains the peptide, scan, charge, and MH

                                    if (htPeptidesFoundForXCorrLevel.ContainsKey(strKey))
                                    {
                                        blnFirstMatchForGroup = false;
                                    }
                                    else
                                    {
                                        htPeptidesFoundForXCorrLevel.Add(strKey, 1);
                                        blnFirstMatchForGroup = true;
                                    }
                                }
                                else
                                {
                                    // New XCorr
                                    // Reset htPeptidesFoundForXCorrLevel
                                    htPeptidesFoundForXCorrLevel.Clear();

                                    // Update strPreviousXCorr
                                    strPreviousXCorr = searchResult.PeptideXCorr;

                                    // Append a new entry to htPeptidesFoundForXCorrLevel
                                    htPeptidesFoundForXCorrLevel.Add(strKey, 1);
                                    blnFirstMatchForGroup = true;
                                }

                                blnSuccess = AddModificationsAndComputeMass(searchResult, blnFirstMatchForGroup);
                                if (!blnSuccess)
                                {
                                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        strErrorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'";
                                        if (!string.IsNullOrEmpty(mErrorMessage))
                                        {
                                            strErrorLog += ": " + mErrorMessage;
                                            mErrorMessage = string.Empty;
                                        }
                                        strErrorLog += "\n";
                                    }
                                }
                                SaveResultsFileEntrySeqInfo(searchResult, blnFirstMatchForGroup);
                            }

                            // Update the progress
                            var sngPercentComplete = Convert.ToSingle(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100);
                            if (CreateProteinModsFile)
                            {
                                sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                            }
                            UpdateProgress(sngPercentComplete);

                            intResultsProcessed += 1;
                        }
                    }

                    if (CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(strInputFilePath);
                        var strModificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        strModificationSummaryFilePath = Path.Combine(strOutputFolderPath, strModificationSummaryFilePath);

                        SaveModificationSummaryFile(strModificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (strErrorLog.Length > 0)
                    {
                        SetErrorMessage("Invalid Lines: " + "\n" + strErrorLog);
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
            string strLineIn,
            IReadOnlyList<int> intColumnMapping,
            clsSearchResultsSequest searchResult,
            ref string strErrorLog)
        {
            string[] strSplitLine = null;

            bool blnValidSearchResult;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                // Reset searchResult
                searchResult.Clear();

                strSplitLine = strLineIn.TrimEnd().Split('\t');
                if (strSplitLine.Length < SEQUEST_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.RowIndex], out int resultId))
                {
                    ReportError("RowIndex column is missing or invalid", true);
                }

                searchResult.ResultID = resultId;

                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.Scan], out string scan);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.NumScans], out string numScans);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.Charge], out string charge);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.PeptideMH], out string peptideMh);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.XCorr], out string peptideXCorr);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.DeltaCn], out string peptideDeltaCn);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.Sp], out string peptideSp);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.ProteinName], out string proteinName);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.MultipleProteinCount], out string multipleProteinCount);

                searchResult.Scan = scan;
                searchResult.NumScans = numScans;
                searchResult.Charge = charge;
                searchResult.PeptideMH = peptideMh;
                searchResult.PeptideXCorr = peptideXCorr;
                searchResult.PeptideDeltaCn = peptideDeltaCn;
                searchResult.PeptideSp = peptideSp;
                searchResult.ProteinName = proteinName;
                searchResult.MultipleProteinCount = multipleProteinCount;

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.PeptideSequence], out string strPeptideSequenceWithMods))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(strPeptideSequenceWithMods, true, true);

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Sequest only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.DeltaCn2], out string peptideDeltaCn2);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.RankSP], out string peptideRankSp);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.RankXC], out string peptideRankXc);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.DelM], out string peptideDeltaMass);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.XcRatio], out string peptideXcRatio);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.PassFilt], out string peptidePassFilt);           // Legacy/Unused
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.MScore], out string peptideMScore);               // Legacy/Unused
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.NTT], out string peptideNtt);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.IonsObserved], out string ionsObserved);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.IonsExpected], out string ionsExpected);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.DelMPPM], out string delMppm);

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

                blnValidSearchResult = true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (strSplitLine != null && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing Sequest Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing Sequest Results in ParseSequestResultsFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">Sequest Synopsis or First-hits file</param>
        /// <param name="strOutputFolderPath">Output folder</param>
        /// <param name="strParameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            var blnSuccess = false;

            if (!LoadParameterFileSettings(strParameterFilePath))
            {
                SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(strInputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                // Set this to true since Sequest param files can have the same mod mass on different residues, and those residues may use different symbols
                mPeptideMods.ConsiderModSymbolWhenFindingIdenticalMods = true;

                blnSuccess = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!blnSuccess)
                {
                    return false;
                }

                ResetProgress("Parsing " + Path.GetFileName(strInputFilePath));

                if (!CleanupFilePaths(ref strInputFilePath, ref strOutputFolderPath))
                {
                    return false;
                }

                // Obtain the full path to the input file
                var inputFile = new FileInfo(strInputFilePath);

                try
                {
                    blnSuccess = ParseSynopsisOrFirstHitsFile(inputFile.FullName, strOutputFolderPath, false);
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error calling ParseSynopsisOrFirstHitsFile" + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    blnSuccess = false;
                }

                if (blnSuccess && CreateProteinModsFile)
                {
                    blnSuccess = CreateProteinModsFileWork(inputFile, strOutputFolderPath);
                }

                if (blnSuccess)
                {
                    OperationComplete();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in clsSequestResultsProcessor.ProcessFile:  " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return blnSuccess;
        }

        private bool CreateProteinModsFileWork(FileInfo inputFile, string strOutputFolderPath)
        {
            bool blnSuccess;

            // First create the MTS PepToProteinMap file using inputFile
            // Will also look for the first hits or synopsis file and use that too if it is present

            var lstSourcePHRPDataFiles = new List<string>();

            string strAdditionalFile;
            var strInputFileBaseName = Path.GetFileNameWithoutExtension(inputFile.Name);

            lstSourcePHRPDataFiles.Add(inputFile.FullName);
            if (strInputFileBaseName.EndsWith(FILENAME_SUFFIX_SYNOPSIS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                strAdditionalFile = ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_FIRST_HITS_FILE);
                if (File.Exists(strAdditionalFile))
                {
                    lstSourcePHRPDataFiles.Add(strAdditionalFile);
                }
            }
            else if (strInputFileBaseName.EndsWith(FILENAME_SUFFIX_FIRST_HITS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                strAdditionalFile = ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_SYNOPSIS_FILE);
                if (File.Exists(strAdditionalFile))
                {
                    lstSourcePHRPDataFiles.Add(strAdditionalFile);
                }
            }

            var strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, strOutputFolderPath, MTS: true);

            if (File.Exists(strMTSPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
            {
                blnSuccess = true;
            }
            else
            {
                // Mapping file does not exist
                blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                if (!blnSuccess)
                {
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                ValidatePHRPReaderSupportFiles(inputFile.FullName, strOutputFolderPath);

                // Now create the Protein Mods file
                blnSuccess = CreateProteinModDetailsFile(inputFile.FullName, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Sequest);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                blnSuccess = true;
            }

            return true;
        }

        private bool ParseSequestSynFileHeaderLine(
            string strLineIn,
            out int[] intColumnMapping)
        {
            // Parse the header line



            var lstColumnNames = new SortedDictionary<string, eSequestSynopsisFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {"HitNum", eSequestSynopsisFileColumns.RowIndex},
                {"ScanNum", eSequestSynopsisFileColumns.Scan},
                {"ScanCount", eSequestSynopsisFileColumns.NumScans},
                {"ChargeState", eSequestSynopsisFileColumns.Charge},
                {"MH", eSequestSynopsisFileColumns.PeptideMH},
                {"XCorr", eSequestSynopsisFileColumns.XCorr},
                {"DelCn", eSequestSynopsisFileColumns.DeltaCn},
                {"Sp", eSequestSynopsisFileColumns.Sp},
                {"Reference", eSequestSynopsisFileColumns.ProteinName},
                {"MultiProtein", eSequestSynopsisFileColumns.MultipleProteinCount},     // Multiple protein count: 0 if the peptide is in 1 protein, 1 if the peptide is in 2 proteins, etc.
                {"Peptide", eSequestSynopsisFileColumns.PeptideSequence},
                {"DelCn2", eSequestSynopsisFileColumns.DeltaCn2},
                {"RankSP", eSequestSynopsisFileColumns.RankSP},
                {"RankXC", eSequestSynopsisFileColumns.RankXC},
                {"DelM", eSequestSynopsisFileColumns.DelM},
                {"XcRatio", eSequestSynopsisFileColumns.XcRatio},
                {"PassFilt", eSequestSynopsisFileColumns.PassFilt},                     // Legacy/unused
                {"MScore", eSequestSynopsisFileColumns.MScore},                         // Legacy/unused
                {"NumTrypticEnds", eSequestSynopsisFileColumns.NTT},
                {"Ions_Observed", eSequestSynopsisFileColumns.IonsObserved},
                {"Ions_Expected", eSequestSynopsisFileColumns.IonsExpected},
                {"DelM_PPM", eSequestSynopsisFileColumns.DelMPPM},
                {"Cleavage_State", eSequestSynopsisFileColumns.Cleavage_State},         // Computed by this program and appended to the input file or saved in a new file
                {"Terminus_State", eSequestSynopsisFileColumns.Terminus_State},         // Computed by this program
                {"Mod_Count", eSequestSynopsisFileColumns.Mod_Count},                   // Computed by this program
                {"Mod_Description", eSequestSynopsisFileColumns.Mod_Description},       // Computed by this program
                {"Monoisotopic_Mass", eSequestSynopsisFileColumns.Monoisotopic_Mass}    // Computed by this program
            };



            intColumnMapping = new int[SequestSynopsisFileColCount];

            try
            {
                // Initialize each entry in intColumnMapping to -1
                for (var intIndex = 0; intIndex <= intColumnMapping.Length - 1; intIndex++)
                {
                    intColumnMapping[intIndex] = -1;
                }

                var strSplitLine = strLineIn.Split('\t');
                for (var intIndex = 0; intIndex <= strSplitLine.Length - 1; intIndex++)
                {
                    if (lstColumnNames.TryGetValue(strSplitLine[intIndex], out var eResultFileColumn))
                    {
                        // Recognized column name; update intColumnMapping
                        intColumnMapping[(int)eResultFileColumn] = intIndex;
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
