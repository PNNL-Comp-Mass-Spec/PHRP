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
using System.Collections;
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
        /// For each mod mass, determine the modification and add to objSearchResult
        /// </summary>
        /// <param name="objSearchResult"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <param name="lstModInfo"></param>
        /// <remarks></remarks>
        private void AddModificationsToResidues(
            clsSearchResultsMSPathFinder objSearchResult,
            bool blnUpdateModOccurrenceCounts,
            IReadOnlyCollection<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            if (string.IsNullOrWhiteSpace(objSearchResult.Modifications))
            {
                return;
            }

            var lstMods = objSearchResult.Modifications.Split(',');
            var finalResidueLoc = objSearchResult.PeptideCleanSequence.Length;

            foreach (var modEntry in lstMods)
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

                foreach (var modDef in lstModInfo)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.InvariantCultureIgnoreCase))
                    {
                        if (!int.TryParse(residueNumber, out var intResidueLocInPeptide))
                        {
                            ReportError("Mod entry does not have a number after the name: " + modEntry);
                            continue;
                        }

                        var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                        if (intResidueLocInPeptide <= 1)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                        }
                        else if (intResidueLocInPeptide >= finalResidueLoc)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                        }

                        // Now that we know the terminus position, assure that intResidueLocInPeptide is 1 not 0
                        if (intResidueLocInPeptide < 1)
                        {
                            intResidueLocInPeptide = 1;
                        }
                        else if (intResidueLocInPeptide > finalResidueLoc)
                        {
                            intResidueLocInPeptide = finalResidueLoc;
                        }

                        var chMostRecentResidue = '-';

                        if (intResidueLocInPeptide >= 1 && intResidueLocInPeptide <= finalResidueLoc)
                        {
                            chMostRecentResidue = objSearchResult.PeptideCleanSequence[intResidueLocInPeptide - 1];
                        }

                        // Associate the mod with the given residue
                        objSearchResult.SearchResultAddModification(modDef.ModMassVal, chMostRecentResidue, intResidueLocInPeptide, eResidueTerminusState, blnUpdateModOccurrenceCounts);

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
            clsSearchResultsMSPathFinder objSearchResult,
            bool blnUpdateModOccurrenceCounts,
            IReadOnlyCollection<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            bool blnSuccess;

            try
            {
                // For other tools, we would add IsotopicMods here
                // This is not supported for MSPathFinder
                //
                // objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

                // Parse .Modifications to determine the modified residues present
                AddModificationsToResidues(objSearchResult, blnUpdateModOccurrenceCounts, lstModInfo);

                // Compute the monoisotopic mass for this peptide
                objSearchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                objSearchResult.UpdateModDescription();

                blnSuccess = true;
            }
            catch (Exception)
            {
                blnSuccess = false;
            }

            return blnSuccess;
        }

        protected double ComputePeptideMass(string strCleanSequence, double dblTotalModMass)
        {
            var dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence);
            dblMass += dblTotalModMass;

            return dblMass;
        }

        /// <summary>
        /// Computes the total of all modifications defined for the sequence
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="modificationList">Comma separated list of modifications, e.g. Dehydro 52,Dehydro 63</param>
        /// <param name="lstModInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected double ComputeTotalModMass(
            string cleanSequence,
            string modificationList,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            if (string.IsNullOrWhiteSpace(modificationList))
            {
                return 0;
            }

            double dblTotalModMass = 0;
            var lstMods = modificationList.Split(',');

            foreach (var modEntry in lstMods)
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

                foreach (var modDef in lstModInfo)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.InvariantCultureIgnoreCase))
                    {
                        dblTotalModMass += modDef.ModMassVal;
                        matchFound = true;
                        break;
                    }
                }

                if (!matchFound)
                {
                    ReportError("Mod name " + modName + " was not defined in the MSPathFinder parameter file; cannot determine the mod mass", true);
                }
            }

            return dblTotalModMass;
        }

        private bool CreateProteinModsFileWork(
            string strBaseName,
            FileInfo inputFile,
            string strSynOutputFilePath,
            string strOutputFolderPath)
        {
            bool blnSuccess;

            // Create the MTSPepToProteinMap file

            var strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS: true);

            var lstSourcePHRPDataFiles = new List<string>();

            if (!string.IsNullOrEmpty(strSynOutputFilePath))
            {
                lstSourcePHRPDataFiles.Add(strSynOutputFilePath);
            }

            if (lstSourcePHRPDataFiles.Count == 0)
            {
                SetErrorMessage("Cannot call CreatePepToProteinMapFile since lstSourcePHRPDataFiles is empty");
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }
            else
            {
                if (File.Exists(strMTSPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
                {
                    blnSuccess = true;
                }
                else
                {
                    blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                    if (!blnSuccess)
                    {
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files to the output folder
                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath);

                // Create the Protein Mods file
                blnSuccess = CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSPathFinder);
            }

            if (!blnSuccess)
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
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFilePath"></param>
        /// <param name="lstModInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreateSynResultsFile(
            string strInputFilePath,
            string strOutputFilePath,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            try
            {
                int[] intColumnMapping = null;
                var strErrorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var srDataFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    using (var swResultFile = new StreamWriter(new FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        var headerParsed = false;
                        var rowNumber = 0;

                        // Initialize the list that will hold all of the records in the MSPathFinder result file
                        var lstSearchResultsUnfiltered = new List<udtMSPathFinderSearchResultType>();

                        // Initialize the list that will hold all of the records that will ultimately be written out to disk
                        var lstFilteredSearchResults = new List<udtMSPathFinderSearchResultType>();

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !AbortProcessing)
                        {
                            var strLineIn = srDataFile.ReadLine();
                            rowNumber += 1;

                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            if (!headerParsed)
                            {
                                // Parse the header line

                                var blnSuccess = ParseMSPathFinderResultsFileHeaderLine(strLineIn, out intColumnMapping);
                                if (!blnSuccess)
                                {
                                    if (string.IsNullOrEmpty(mErrorMessage))
                                    {
                                        SetErrorMessage("Invalid header line in " + Path.GetFileName(strInputFilePath));
                                    }
                                    return false;
                                }

                                // Write the header line to the output file
                                WriteSynFHTFileHeader(swResultFile, ref strErrorLog);

                                headerParsed = true;
                                continue;
                            }

                            var udtSearchResult = new udtMSPathFinderSearchResultType();

                            var blnValidSearchResult = ParseMSPathFinderResultsFileEntry(strLineIn, ref udtSearchResult, ref strErrorLog, intColumnMapping, lstModInfo, rowNumber);

                            if (blnValidSearchResult)
                            {
                                lstSearchResultsUnfiltered.Add(udtSearchResult);
                            }

                            // Update the progress
                            var sngPercentComplete = Convert.ToSingle(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100);
                            if (CreateProteinModsFile)
                            {
                                sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                            }
                            UpdateProgress(sngPercentComplete);
                        }

                        // Sort the SearchResults by scan, charge, and ascending SpecEValue
                        lstSearchResultsUnfiltered.Sort(new MSPathFinderSearchResultsComparerScanChargeScorePeptide());

                        // Now filter the data

                        // Initialize variables
                        var intStartIndex = 0;

                        while (intStartIndex < lstSearchResultsUnfiltered.Count)
                        {
                            var intEndIndex = intStartIndex;
                            while (intEndIndex + 1 < lstSearchResultsUnfiltered.Count &&
                                   lstSearchResultsUnfiltered[intEndIndex + 1].ScanNum == lstSearchResultsUnfiltered[intStartIndex].ScanNum)
                            {
                                intEndIndex += 1;
                            }

                            // Store the results for this scan
                            StoreSynMatches(lstSearchResultsUnfiltered, intStartIndex, intEndIndex, lstFilteredSearchResults);

                            intStartIndex = intEndIndex + 1;
                        }

                        // Sort the data in udtFilteredSearchResults then write out to disk
                        SortAndWriteFilteredSearchResults(swResultFile, lstFilteredSearchResults, ref strErrorLog);
                    }
                }

                // Inform the user if any errors occurred
                if (strErrorLog.Length > 0)
                {
                    SetErrorMessage("Invalid Lines: " + "\n" + strErrorLog);
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

        private bool ExtractModInfoFromParamFile(string strMSGFDBParamFilePath,
            out List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            // The DMS-based parameter file for MSPathFinder uses the same formatting as MSGF+

            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MSPathFinder");
            RegisterEvents(modFileProcessor);

            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            // Note that this call will initialize lstModInfo
            var success = modFileProcessor.ExtractModInfoFromParamFile(strMSGFDBParamFilePath, out lstModInfo);

            if (!success || mErrorCode != ePHRPErrorCodes.NoError)
            {
                if (mErrorCode == ePHRPErrorCodes.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MSPathFinder parameter file");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFDBModsWithModDefinitions(lstModInfo, mPeptideMods);

            return true;
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFolderPath"></param>
        /// <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <param name="lstModInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool ParseMSPathfinderSynopsisFile(
            string strInputFilePath,
            string strOutputFolderPath,
            bool blnResetMassCorrectionTagsAndModificationDefinitions,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that ParseMSPathfinderSynopsisFile synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            int[] intColumnMapping = null;
            bool blnSuccess;

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
                var objSearchResult = new clsSearchResultsMSPathFinder(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize two hashtables
                var htPeptidesFoundForSpecEValue = new Hashtable();
                var htPeptidesFoundForQValue = new Hashtable();

                var blnFirstMatchForGroup = false;

                var strPreviousSpecEValue = string.Empty;
                var strPreviousQValue = string.Empty;

                var strErrorLog = string.Empty;

                try
                {
                    objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var srDataFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        var intResultsProcessed = 0;
                        var blnHeaderParsed = false;

                        // Create the output files
                        var strBaseOutputFilePath = Path.Combine(strOutputFolderPath, Path.GetFileName(strInputFilePath));
                        blnSuccess = InitializeSequenceOutputFiles(strBaseOutputFilePath);

                        // Parse the input file

                        while (!srDataFile.EndOfStream & !AbortProcessing)
                        {
                            var strLineIn = srDataFile.ReadLine();

                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseMSPathFinderSynFileHeaderLine(strLineIn, out intColumnMapping);
                                if (!blnSuccess)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                blnHeaderParsed = true;
                                continue;
                            }

                            var blnValidSearchResult = ParseMSPathFinderSynFileEntry(strLineIn, objSearchResult, ref strErrorLog,
                                                                                     intResultsProcessed, intColumnMapping,
                                                                                     out _);

                            if (!blnValidSearchResult)
                            {
                                continue;
                            }

                            var strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.Charge;

                            var blnNewValue = true;

                            if (string.IsNullOrEmpty(objSearchResult.SpecEValue))
                            {
                                if (objSearchResult.QValue == strPreviousQValue)
                                {
                                    // New result has the same QValue as the previous result
                                    // See if htPeptidesFoundForQValue contains the peptide, scan and charge

                                    if (htPeptidesFoundForQValue.ContainsKey(strKey))
                                    {
                                        blnFirstMatchForGroup = false;
                                    }
                                    else
                                    {
                                        htPeptidesFoundForQValue.Add(strKey, 1);
                                        blnFirstMatchForGroup = true;
                                    }

                                    blnNewValue = false;
                                }
                            }
                            else if (objSearchResult.SpecEValue == strPreviousSpecEValue)
                            {
                                // New result has the same SpecEValue as the previous result
                                // See if htPeptidesFoundForSpecEValue contains the peptide, scan and charge

                                if (htPeptidesFoundForSpecEValue.ContainsKey(strKey))
                                {
                                    blnFirstMatchForGroup = false;
                                }
                                else
                                {
                                    htPeptidesFoundForSpecEValue.Add(strKey, 1);
                                    blnFirstMatchForGroup = true;
                                }

                                blnNewValue = false;
                            }

                            if (blnNewValue)
                            {
                                // New SpecEValue or new QValue
                                // Reset the hashtables
                                htPeptidesFoundForSpecEValue.Clear();
                                htPeptidesFoundForQValue.Clear();

                                // Update the cached values
                                strPreviousSpecEValue = objSearchResult.SpecEValue;
                                strPreviousQValue = objSearchResult.QValue;

                                // Append a new entry to the hashtables
                                htPeptidesFoundForSpecEValue.Add(strKey, 1);
                                htPeptidesFoundForQValue.Add(strKey, 1);

                                blnFirstMatchForGroup = true;
                            }

                            blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup, lstModInfo);
                            if (!blnSuccess)
                            {
                                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                {
                                    strErrorLog += "Error adding modifications to sequence at RowIndex '" + objSearchResult.ResultID + "'" + "\n";
                                }
                            }

                            SaveResultsFileEntrySeqInfo(objSearchResult, blnFirstMatchForGroup);

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
                    }

                    blnSuccess = true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    blnSuccess = false;
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
                blnSuccess = false;
            }

            return blnSuccess;
        }

        /// <summary>
        /// Parse a MSPathFinder results line while creating the MSPathFinder synopsis file
        /// </summary>
        /// <param name="strLineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="strErrorLog"></param>
        /// <param name="intColumnMapping"></param>
        /// <param name="lstModInfo"></param>
        /// <param name="rowNumber">Row number (used for error reporting)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool ParseMSPathFinderResultsFileEntry(
            string strLineIn,
            ref udtMSPathFinderSearchResultType udtSearchResult,
            ref string strErrorLog,
            IList<int> intColumnMapping,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo,
            int rowNumber)
        {
            // Parses an entry from the MSPathFinder results file

            double dblSequenceMonoMassMSPathFinder = 0;    // Theoretical peptide monoisotopic mass, including mods, as computed by MSPathFinder

            bool blnValidSearchResult;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                udtSearchResult.Clear();
                var strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 11)
                {
                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.Scan], out udtSearchResult.Scan))
                    {
                        ReportError("Scan column is missing or invalid in row " + rowNumber, true);
                    }

                    if (!int.TryParse(udtSearchResult.Scan, out udtSearchResult.ScanNum))
                    {
                        ReportError("Scan column is not numeric in row " + rowNumber, true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MSPathFinder
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);
                    double.TryParse(udtSearchResult.CalculatedMonoMass, out dblSequenceMonoMassMSPathFinder);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.PrefixResidue], out udtSearchResult.PrefixResidue);

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.Sequence], out udtSearchResult.Sequence))
                    {
                        ReportError("Sequence column is missing or invalid in row " + rowNumber, true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.SuffixResidue], out udtSearchResult.SuffixResidue);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.Modifications], out udtSearchResult.Modifications);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.Composition], out udtSearchResult.Composition);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.Protein], out udtSearchResult.Protein);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.ProteinDesc], out udtSearchResult.ProteinDesc);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.ProteinLength], out udtSearchResult.ProteinLength);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.ResidueEnd], out udtSearchResult.ResidueEnd);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.ResidueStart], out udtSearchResult.ResidueStart);

                    // Parse the list of modified residues to determine the total mod mass
                    var dblTotalModMass = ComputeTotalModMass(udtSearchResult.Sequence, udtSearchResult.Modifications, lstModInfo);

                    // Compute monoisotopic mass of the peptide
                    udtSearchResult.CalculatedMonoMassPHRP = ComputePeptideMass(udtSearchResult.Sequence, dblTotalModMass);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.MostAbundantIsotopeMz], out udtSearchResult.MostAbundantIsotopeMz);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.NumMatchedFragments], out udtSearchResult.NumMatchedFragments);

                    if (GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.SpecEValue], out udtSearchResult.SpecEValue))
                    {
                        double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum);
                    }
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.EValue], out udtSearchResult.EValue);

                    if (GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.QValue], out udtSearchResult.QValue))
                    {
                        double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum);
                    }
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderResultsFileColumns.PepQValue], out udtSearchResult.PepQValue);

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MassMSPathFinder results file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    strErrorLog += "Error parsing MassMSPathFinder Results in ParseMSPathFinderResultsFileEntry for Row " + rowNumber + "\n";
                }

                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        /// <summary>
        /// Parse the MSPathFinder results file header line
        /// </summary>
        /// <param name="strLineIn"></param>
        /// <param name="intColumnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        /// <remarks></remarks>
        private bool ParseMSPathFinderResultsFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            // The expected column order from MassMSPathFinder:
            //   Scan	Pre	Sequence	Post	Modifications	Composition	ProteinName	ProteinDesc	ProteinLength	Start	End	Charge	MostAbundantIsotopeMz	Mass	#MatchedFragments	Probability SpecEValue    EValue    QValue    PepQValue

            var lstColumnNames = new SortedDictionary<string, eMSPathFinderResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase)
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

            intColumnMapping = new int[MSPathFinderResultsFileColCount];

            try
            {
                // Initialize each entry in intColumnMapping to -1
                for (var intIndex = 0; intIndex <= intColumnMapping.Length - 1; intIndex++)
                {
                    intColumnMapping[intIndex] = -1;
                }

                var strSplitLine = strLineIn.Split('\t');
                var blnUseDefaultHeaders = false;

                if (strSplitLine.Length >= 2)
                {
                    if (int.TryParse(strSplitLine[1], out _))
                    {
                        // Second column has a number; this is not a header line
                        blnUseDefaultHeaders = true;
                    }
                    else
                    {
                        for (var intIndex = 0; intIndex <= strSplitLine.Length - 1; intIndex++)
                        {
                            if (lstColumnNames.TryGetValue(strSplitLine[intIndex], out var eResultFileColumn))
                            {
                                // Recognized column name; update intColumnMapping
                                intColumnMapping[(int)eResultFileColumn] = intIndex;
                                blnUseDefaultHeaders = false;
                            }
                            else
                            {
                                // Unrecognized column name
                                OnWarningEvent("Warning: Unrecognized column header name '" + strSplitLine[intIndex] + "' in ParseMSPathFinderResultsFileHeaderLine");
                            }
                        }
                    }
                }

                if (blnUseDefaultHeaders)
                {
                    // Use default column mappings
                    for (var intIndex = 0; intIndex <= intColumnMapping.Length - 1; intIndex++)
                    {
                        intColumnMapping[intIndex] = intIndex;
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

        private bool ParseMSPathFinderSynFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            var lstColumnNames = new SortedDictionary<string, eMSPathFinderSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
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

            intColumnMapping = new int[MSPathFinderSynFileColCount];

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
                SetErrorMessage("Error parsing header in MSPathFinder synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMSPathFinderSynFileEntry(
            string strLineIn,
            clsSearchResultsMSPathFinder objSearchResult,
            ref string strErrorLog,
            int intResultsProcessed,
            IReadOnlyList<int> intColumnMapping,
            out string strPeptideSequence)
        {
            // Parses an entry from the MSPathFinder Synopsis file

            string[] strSplitLine = null;

            // Reset objSearchResult
            objSearchResult.Clear();
            strPeptideSequence = string.Empty;

            try
            {

                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length < 15)
                {
                    return false;
                }

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.ResultID], out string strValue))
                {
                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        strErrorLog += "Error reading ResultID value from MSPathFinder Results line " + (intResultsProcessed + 1).ToString() + "\n";
                    }
                    return false;
                }

                objSearchResult.ResultID = int.Parse(strValue);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.Scan], out string scan);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.Charge], out string charge);

                objSearchResult.Scan = scan;
                objSearchResult.Charge = charge;

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.Sequence], out strPeptideSequence))
                {
                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        strErrorLog += "Error reading Peptide sequence value from MSPathFinder Results line " + (intResultsProcessed + 1).ToString() + "\n";
                    }
                    return false;
                }

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.Protein], out string proteinName);
                objSearchResult.ProteinName = proteinName;
                objSearchResult.MultipleProteinCount = "0";

                objSearchResult.PeptideDeltaMass = "0";

                // Note that MSPathFinder sequences don't actually have mod symbols; that information is tracked via objSearchResult.Modifications

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                objSearchResult.SetPeptideSequenceWithMods(strPeptideSequence, true, true);

                var objSearchResultBase = (clsSearchResultsBaseClass)objSearchResult;

                ComputePseudoPeptideLocInProtein(objSearchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                objSearchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.MostAbundantIsotopeMz], out string mostAbundantIsotopeMz);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.Modifications], out string modifications);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.Composition], out string composition);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.ProteinDesc], out string proteinDesc);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.ProteinLength], out string proteinLength);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.ResidueStart], out string residueStart);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.ResidueEnd], out string residueEnd);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.MatchedFragments], out string matchedFragments);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.SpecEValue], out string specEValue);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.EValue], out string eValue);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.QValue], out string qValue);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMSPathFinderSynFileColumns.PepQValue], out string pepQValue);

                objSearchResult.MostAbundantIsotopeMz = mostAbundantIsotopeMz;
                objSearchResult.Modifications = modifications;
                objSearchResult.Composition = composition;
                objSearchResult.ProteinDesc = proteinDesc;
                objSearchResult.ProteinLength = proteinLength;
                objSearchResult.ResidueStart = residueStart;
                objSearchResult.ResidueEnd = residueEnd;
                objSearchResult.MatchedFragments = matchedFragments;
                objSearchResult.SpecEValue = specEValue;
                objSearchResult.EValue = eValue;
                objSearchResult.QValue = qValue;
                objSearchResult.PepQValue = pepQValue;

                // Update the peptide location in the protein
                if (!string.IsNullOrEmpty(objSearchResult.ResidueStart))
                {
                    int.TryParse(objSearchResult.ResidueStart, out var peptideLocInProteinStart);
                    objSearchResult.PeptideLocInProteinStart = peptideLocInProteinStart;
                }

                if (!string.IsNullOrEmpty(objSearchResult.ResidueEnd))
                {
                    int.TryParse(objSearchResult.ResidueEnd, out var peptideLocInProteinEnd);
                    objSearchResult.PeptideLocInProteinEnd = peptideLocInProteinEnd;
                }

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (strSplitLine != null && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing MSPathFinder Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MSPathFinder Results in ParseMSPathFinderSynFileEntry" + "\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">MSPathFinder results file (Dataset_IcTda.tsv)</param>
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

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(strInputFilePath);

                    // Load the MSPathFinder Parameter File so that we can determine the modification names and masses
                    // Note that this call will initialize lstModInfo
                    var success = ExtractModInfoFromParamFile(SearchToolParameterFilePath, out var lstModInfo);
                    if (!success)
                    {
                        return false;
                    }

                    // Define the base output filename using strInputFilePath
                    var strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath);

                    // Auto-replace "_IcTda.tsv" with "_mspath"
                    if (strBaseName.EndsWith("_IcTda", StringComparison.InvariantCultureIgnoreCase))
                    {
                        strBaseName = strBaseName.Substring(0, strBaseName.Length - "_IcTda".Length) + "_mspath";
                    }

                    // Do not create a first-hits file for MSPathFinder results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_mspath_syn.txt
                    var strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath, lstModInfo);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(strSynOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to create the other PHRP files
                    blnSuccess = ParseMSPathfinderSynopsisFile(strSynOutputFilePath, strOutputFolderPath, false, lstModInfo);

                    if (blnSuccess && CreateProteinModsFile)
                    {
                        // Check for an empty synopsis file
                        if (!ValidateFileHasData(strSynOutputFilePath, "Synopsis file", out var errorMessage))
                        {
                            ReportWarning(errorMessage);
                        }
                        else
                        {
                            blnSuccess = CreateProteinModsFileWork(strBaseName, inputFile, strSynOutputFilePath, strOutputFolderPath);
                        }
                    }

                    if (blnSuccess)
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

            return blnSuccess;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter swResultFile,
            List<udtMSPathFinderSearchResultType> lstFilteredSearchResults,
            ref string strErrorLog)
        {
            // Sort udtFilteredSearchResults by ascending SpecEValue, QValue, Scan, Peptide, and Protein
            var query = from item in lstFilteredSearchResults orderby item.SpecEValueNum, item.QValueNum, item.ScanNum, item.Sequence, item.Protein select item;

            var intIndex = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(intIndex, swResultFile, result, ref strErrorLog);
                intIndex += 1;
            }
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan (typically there will only be one result for MSPathFinder)
        /// </summary>
        /// <param name="lstSearchResults">Search results</param>
        /// <param name="intStartIndex">Start index for data in this scan</param>
        /// <param name="intEndIndex">End index for data in this scan</param>
        /// <param name="lstFilteredSearchResults">Output parmaeter: the actual filtered search results</param>
        /// <remarks></remarks>
        private void StoreSynMatches(
            IList<udtMSPathFinderSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex,
            List<udtMSPathFinderSearchResultType> lstFilteredSearchResults)
        {
            // If there was more than one result, we could rank them by score
            // AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex)

            // The calling procedure already sorted by scan, charge, and Score; no need to re-sort

            ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1);

            // Now store the matches that pass the filters
            //  Either SpecEValue < 5E-07 (0.0000005)
            //  or     QValue < 10
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                if (lstSearchResults[intIndex].SpecEValueNum <= MSGFDBSynopsisFileSpecEValueThreshold || lstSearchResults[intIndex].QValueNum < 0.1)
                {
                    lstFilteredSearchResults.Add(lstSearchResults[intIndex]);
                }
            }
        }

        private void WriteSynFHTFileHeader(
            TextWriter swResultFile,
            ref string strErrorLog)
        {
            // Write out the header line for synopsis / first hits files
            try
            {
                var lstData = new List<string>
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

                swResultFile.WriteLine(CollapseList(lstData));
            }
            catch (Exception)
            {
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    strErrorLog += "Error writing synopsis / first hits header" + "\n";
                }
            }
        }

        /// <summary>
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="intResultID"></param>
        /// <param name="swResultFile"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="strErrorLog"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int intResultID,
            TextWriter swResultFile,
            udtMSPathFinderSearchResultType udtSearchResult,
            ref string strErrorLog)
        {
            try
            {
                var lstData = new List<string>
                {
                    intResultID.ToString(),
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

                swResultFile.WriteLine(CollapseList(lstData));
            }
            catch (Exception)
            {
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    strErrorLog += "Error writing synopsis / first hits record" + "\n";
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
