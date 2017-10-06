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
// E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
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
        public clsSequestResultsProcessor() : base()
        {
            base.mFileDate = "June 28, 2013";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"
        public const string FILENAME_SUFFIX_FIRST_HITS_FILE = "_fht";
        public const string FILENAME_SUFFIX_SYNOPSIS_FILE = "_syn";

        private const int SEQUEST_SYN_FILE_MIN_COL_COUNT = 5;
        private const int MAX_ERROR_LOG_LENGTH = 4096;

        // These columns correspond to the tab-delimited file created directly by SEQUEST
        protected const int SequestSynopsisFileColCount = 27;
        public enum eSequestSynopsisFileColumns : int
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
            clsSearchResultsSequest objSearchResult,
            bool blnUpdateModOccurrenceCounts)
        {
            // Step through .PeptideSequenceWithMods
            // For each residue, check if a static mod is defined that affects that residue
            // For each mod symbol, determine the modification and add to objSearchResult

            int intModIndex = 0;
            char chChar = default(char);
            clsModificationDefinition objModificationDefinition = default(clsModificationDefinition);

            string strSequence = null;
            char chMostRecentLetter = default(char);
            int intResidueLocInPeptide = 0;

            bool blnSuccess = false;

            chMostRecentLetter = '-';
            intResidueLocInPeptide = 0;

            // Assume success for now
            blnSuccess = true;

            strSequence = objSearchResult.PeptideSequenceWithMods;
            for (var intIndex = 0; intIndex <= strSequence.Length - 1; intIndex++)
            {
                chChar = strSequence[intIndex];

                if (IsLetterAtoZ(chChar))
                {
                    chMostRecentLetter = chChar;
                    intResidueLocInPeptide += 1;

                    for (intModIndex = 0; intModIndex <= mPeptideMods.ModificationCount - 1; intModIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(intModIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);

                            if (objModificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                blnSuccess = objSearchResult.SearchResultAddModification(ref objModificationDefinition, chChar, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);

                                if (!blnSuccess)
                                {
                                    // Error adding this static mod
                                    SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
                                    mErrorMessage = "Error calling objSearchResult.SearchResultAddModification for peptide '" + strSequence + "': " + objSearchResult.ErrorMessage;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(chMostRecentLetter))
                {
                    blnSuccess = objSearchResult.SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);

                    if (!blnSuccess)
                    {
                        // Error adding this dynamic mod
                        SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
                        mErrorMessage = "Error calling objSearchResult.SearchResultAddDynamicModification for peptide '" + strSequence + "': " + objSearchResult.ErrorMessage;
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

        private bool AddModificationsAndComputeMass(clsSearchResultsSequest objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool blnSuccess = false;

            try
            {
                // Assume success for now
                blnSuccess = true;

                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts);

                // Parse .PeptideSequenceWithMods to determine the modified residues present
                if (!AddDynamicAndStaticResidueMods(objSearchResult, blnUpdateModOccurrenceCounts))
                {
                    blnSuccess = false;
                }

                // Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
                // Since Sequest allows a terminal peptide residue to be modified twice, we'll allow that to happen,
                //  even though, biologically, that's typically not possible
                // However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus
                //  (where two COOH groups are present)
                objSearchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                objSearchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                objSearchResult.UpdateModDescription();
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

            string strPreviousXCorr = null;

            // Note that synopsis files are normally sorted on XCorr descending, with lines
            //  duplicated when peptide search results are mapped to multiple proteins
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file,
            //  we will keep track of the scan, charge, and peptide information parsed for each unique XCorr encountered

            Hashtable htPeptidesFoundForXCorrLevel = default(Hashtable);

            string strKey = null;

            string strLineIn = null;
            string strModificationSummaryFilePath = null;

            clsSearchResultsSequest objSearchResult = default(clsSearchResultsSequest);

            int intResultsProcessed = 0;
            float sngPercentComplete = 0;

            bool blnSuccess = false;
            bool blnDataLine = false;
            bool blnValidSearchResult = false;
            bool blnFirstMatchForGroup = false;

            bool blnHeaderParsed = false;
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

                // Initialize objSearchResult
                objSearchResult = new clsSearchResultsSequest(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForXCorrLevel
                htPeptidesFoundForXCorrLevel = new Hashtable();
                strPreviousXCorr = string.Empty;

                try
                {
                    objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    string strErrorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var srDataFile = new StreamReader(strInputFilePath))
                    {
                        intResultsProcessed = 0;
                        blnHeaderParsed = false;

                        // Create the output files
                        string strBaseOutputFilePath = Path.Combine(strOutputFolderPath, Path.GetFileName(strInputFilePath));
                        blnSuccess = base.InitializeSequenceOutputFiles(strBaseOutputFilePath);

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !base.AbortProcessing)
                        {
                            strLineIn = srDataFile.ReadLine();
                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            blnDataLine = true;

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseSequestSynFileHeaderLine(strLineIn, ref intColumnMapping);
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

                            if (blnDataLine)
                            {
                                blnValidSearchResult = ParseSequestResultsFileEntry(ref strLineIn, ref intColumnMapping, objSearchResult, ref strErrorLog);
                            }
                            else
                            {
                                blnValidSearchResult = false;
                            }

                            if (blnValidSearchResult)
                            {
                                strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.NumScans + "_" + objSearchResult.Charge + "_" + objSearchResult.PeptideMH;

                                if (objSearchResult.PeptideXCorr == strPreviousXCorr)
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
                                    strPreviousXCorr = objSearchResult.PeptideXCorr;

                                    // Append a new entry to htPeptidesFoundForXCorrLevel
                                    htPeptidesFoundForXCorrLevel.Add(strKey, 1);
                                    blnFirstMatchForGroup = true;
                                }

                                blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup);
                                if (!blnSuccess)
                                {
                                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        strErrorLog += "Error adding modifications to sequence at RowIndex '" + objSearchResult.ResultID + "'";
                                        if ((mErrorMessage != null) && mErrorMessage.Length > 0)
                                        {
                                            strErrorLog += ": " + mErrorMessage;
                                            mErrorMessage = string.Empty;
                                        }
                                        strErrorLog += "\n";
                                    }
                                }
                                base.SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass)objSearchResult, blnFirstMatchForGroup);
                            }

                            // Update the progress
                            sngPercentComplete = Convert.ToSingle(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100);
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
                        var fiInputFile = new FileInfo(strInputFilePath);
                        strModificationSummaryFilePath = Path.GetFileName(base.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        strModificationSummaryFilePath = Path.Combine(strOutputFolderPath, strModificationSummaryFilePath);

                        SaveModificationSummaryFile(strModificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (strErrorLog.Length > 0)
                    {
                        SetErrorMessage("Invalid Lines: " + "\n" + strErrorLog);
                        blnSuccess = false;
                    }
                    else
                    {
                        blnSuccess = true;
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    blnSuccess = false;
                }
                finally
                {
                    base.CloseSequenceOutputFiles();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private bool ParseSequestResultsFileEntry(
            ref string strLineIn,
            ref int[] intColumnMapping,
            clsSearchResultsSequest objSearchResult,
            ref string strErrorLog)
        {
            string[] strSplitLine = null;
            string strPeptideSequenceWithMods = string.Empty;

            bool blnValidSearchResult = false;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                // Reset objSearchResult
                objSearchResult.Clear();

                strSplitLine = strLineIn.Trim().Split('\t');
                if (strSplitLine.Length < SEQUEST_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.RowIndex], out int resultId))
                {
                    ReportError("RowIndex column is missing or invalid", true);
                }

                objSearchResult.ResultID = resultId;

                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.Scan], out string scan);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.NumScans], out string numScans);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.Charge], out string charge);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.PeptideMH], out string peptideMh);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.XCorr], out string peptideXCorr);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.DeltaCn], out string peptideDeltaCn);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.Sp], out string peptideSp);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.ProteinName], out string proteinName);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.MultipleProteinCount], out string multipleProteinCount);

                objSearchResult.Scan = scan;
                objSearchResult.NumScans = numScans;
                objSearchResult.Charge = charge;
                objSearchResult.PeptideMH = peptideMh;
                objSearchResult.PeptideXCorr = peptideXCorr;
                objSearchResult.PeptideDeltaCn = peptideDeltaCn;
                objSearchResult.PeptideSp = peptideSp;
                objSearchResult.ProteinName = proteinName;
                objSearchResult.MultipleProteinCount = multipleProteinCount;

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eSequestSynopsisFileColumns.PeptideSequence], out strPeptideSequenceWithMods))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                objSearchResult.SetPeptideSequenceWithMods(strPeptideSequenceWithMods, true, true);

                clsSearchResultsBaseClass objSearchResultBase = default(clsSearchResultsBaseClass);
                objSearchResultBase = (clsSearchResultsBaseClass)objSearchResult;

                base.ComputePseudoPeptideLocInProtein(objSearchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Sequest only outputs the prefix and suffix letters for the first protein
                objSearchResult.ComputePeptideCleavageStateInProtein();

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

                objSearchResult.PeptideDeltaCn2 = peptideDeltaCn2;
                objSearchResult.PeptideRankSP = peptideRankSp;
                objSearchResult.PeptideRankXC = peptideRankXc;
                objSearchResult.PeptideDeltaMass = peptideDeltaMass;
                objSearchResult.PeptideXcRatio = peptideXcRatio;
                objSearchResult.PeptidePassFilt = peptidePassFilt;
                objSearchResult.PeptideMScore = peptideMScore;
                objSearchResult.PeptideNTT = peptideNtt;
                objSearchResult.IonsObserved = ionsObserved;
                objSearchResult.IonsExpected = ionsExpected;
                objSearchResult.DelMPPM = delMppm;

                blnValidSearchResult = true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if ((strSplitLine != null) && strSplitLine.Length > 0)
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
            bool blnSuccess = false;

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

                base.ResetProgress("Parsing " + Path.GetFileName(strInputFilePath));

                if (!CleanupFilePaths(ref strInputFilePath, ref strOutputFolderPath))
                {
                    return false;
                }

                // Obtain the full path to the input file
                var fiInputFile = new FileInfo(strInputFilePath);

                try
                {
                    blnSuccess = ParseSynopsisOrFirstHitsFile(fiInputFile.FullName, strOutputFolderPath, false);
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error calling ParseSynopsisOrFirstHitsFile" + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    blnSuccess = false;
                }

                if (blnSuccess && CreateProteinModsFile)
                {
                    blnSuccess = CreateProteinModsFileWork(fiInputFile, strOutputFolderPath);
                }

                if (blnSuccess)
                {
                    base.OperationComplete();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in clsSequestResultsProcessor.ProcessFile:  " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return blnSuccess;
        }

        private bool CreateProteinModsFileWork(FileInfo fiInputFile, string strOutputFolderPath)
        {
            bool blnSuccess = false;

            // First create the MTS PepToProteinMap file using fiInputFile
            // Will also look for the first hits or synopsis file and use that too if it is present

            var lstSourcePHRPDataFiles = new List<string>();

            string strAdditionalFile = null;
            string strInputFileBaseName = Path.GetFileNameWithoutExtension(fiInputFile.Name);

            lstSourcePHRPDataFiles.Add(fiInputFile.FullName);
            if (strInputFileBaseName.EndsWith(FILENAME_SUFFIX_SYNOPSIS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                strAdditionalFile = base.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_FIRST_HITS_FILE);
                if (File.Exists(strAdditionalFile))
                {
                    lstSourcePHRPDataFiles.Add(strAdditionalFile);
                }
            }
            else if (strInputFileBaseName.EndsWith(FILENAME_SUFFIX_FIRST_HITS_FILE, StringComparison.InvariantCultureIgnoreCase))
            {
                strAdditionalFile = base.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_SYNOPSIS_FILE);
                if (File.Exists(strAdditionalFile))
                {
                    lstSourcePHRPDataFiles.Add(strAdditionalFile);
                }
            }

            var strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(fiInputFile.FullName, strOutputFolderPath, MTS: true);

            if (File.Exists(strMTSPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
            {
                blnSuccess = true;
            }
            else
            {
                // Mapping file does not exist
                blnSuccess = base.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                if (!blnSuccess)
                {
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                base.ValidatePHRPReaderSupportFiles(fiInputFile.FullName, strOutputFolderPath);

                // Now create the Protein Mods file
                blnSuccess = base.CreateProteinModDetailsFile(fiInputFile.FullName, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Sequest);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                blnSuccess = true;
            }

            return blnSuccess;
        }

        private bool ParseSequestSynFileHeaderLine(
            string strLineIn,
            ref int[] intColumnMapping)
        {
            // Parse the header line

            string[] strSplitLine = null;
            eSequestSynopsisFileColumns eResultFileColumn = default(eSequestSynopsisFileColumns);
            SortedDictionary<string, eSequestSynopsisFileColumns> lstColumnNames = default(SortedDictionary<string, eSequestSynopsisFileColumns>);
            lstColumnNames = new SortedDictionary<string, eSequestSynopsisFileColumns>(StringComparer.InvariantCultureIgnoreCase);

            intColumnMapping = new int[SequestSynopsisFileColCount];

            lstColumnNames.Add("HitNum", eSequestSynopsisFileColumns.RowIndex);
            lstColumnNames.Add("ScanNum", eSequestSynopsisFileColumns.Scan);
            lstColumnNames.Add("ScanCount", eSequestSynopsisFileColumns.NumScans);
            lstColumnNames.Add("ChargeState", eSequestSynopsisFileColumns.Charge);
            lstColumnNames.Add("MH", eSequestSynopsisFileColumns.PeptideMH);
            lstColumnNames.Add("XCorr", eSequestSynopsisFileColumns.XCorr);
            lstColumnNames.Add("DelCn", eSequestSynopsisFileColumns.DeltaCn);
            lstColumnNames.Add("Sp", eSequestSynopsisFileColumns.Sp);
            lstColumnNames.Add("Reference", eSequestSynopsisFileColumns.ProteinName);
            lstColumnNames.Add("MultiProtein", eSequestSynopsisFileColumns.MultipleProteinCount);                     // Multiple protein count: 0 if the peptide is in 1 protein, 1 if the peptide is in 2 proteins, etc.
            lstColumnNames.Add("Peptide", eSequestSynopsisFileColumns.PeptideSequence);
            lstColumnNames.Add("DelCn2", eSequestSynopsisFileColumns.DeltaCn2);
            lstColumnNames.Add("RankSP", eSequestSynopsisFileColumns.RankSP);
            lstColumnNames.Add("RankXC", eSequestSynopsisFileColumns.RankXC);
            lstColumnNames.Add("DelM", eSequestSynopsisFileColumns.DelM);
            lstColumnNames.Add("XcRatio", eSequestSynopsisFileColumns.XcRatio);
            lstColumnNames.Add("PassFilt", eSequestSynopsisFileColumns.PassFilt);                // Legacy/unused
            lstColumnNames.Add("MScore", eSequestSynopsisFileColumns.MScore);                    // Legacy/unused
            lstColumnNames.Add("NumTrypticEnds", eSequestSynopsisFileColumns.NTT);
            lstColumnNames.Add("Ions_Observed", eSequestSynopsisFileColumns.IonsObserved);
            lstColumnNames.Add("Ions_Expected", eSequestSynopsisFileColumns.IonsExpected);
            lstColumnNames.Add("DelM_PPM", eSequestSynopsisFileColumns.DelMPPM);

            // The following columns are computed by this program and appended to the input file or saved in a new file
            lstColumnNames.Add("Cleavage_State", eSequestSynopsisFileColumns.Cleavage_State);
            lstColumnNames.Add("Terminus_State", eSequestSynopsisFileColumns.Terminus_State);
            lstColumnNames.Add("Mod_Count", eSequestSynopsisFileColumns.Mod_Count);
            lstColumnNames.Add("Mod_Description", eSequestSynopsisFileColumns.Mod_Description);
            lstColumnNames.Add("Monoisotopic_Mass", eSequestSynopsisFileColumns.Monoisotopic_Mass);

            try
            {
                // Initialize each entry in intColumnMapping to -1
                for (var intIndex = 0; intIndex <= intColumnMapping.Length - 1; intIndex++)
                {
                    intColumnMapping[intIndex] = -1;
                }

                strSplitLine = strLineIn.Split('\t');
                for (var intIndex = 0; intIndex <= strSplitLine.Length - 1; intIndex++)
                {
                    if (lstColumnNames.TryGetValue(strSplitLine[intIndex], out eResultFileColumn))
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
