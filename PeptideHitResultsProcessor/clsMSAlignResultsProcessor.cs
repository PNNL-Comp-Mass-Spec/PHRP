// This class reads in an MSalign results file (txt format) and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 11/28/2012
//
// E-mail: matthew.monroe@pnnl.gov
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
    public class clsMSAlignResultsProcessor : clsPHRPBaseClass
    {
        public clsMSAlignResultsProcessor() : base()
        {
            base.mFileDate = "September 3, 2014";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_MSALIGN_FILE = "_MSAlign_ResultTable";

        public const string N_TERMINUS_SYMBOL_MSALIGN = ".";
        public const string C_TERMINUS_SYMBOL_MSALIGN = ".";

        public const float DEFAULT_SYN_FILE_PVALUE_THRESHOLD = 0.95f;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string MSALIGN_MOD_MASS_REGEX = @"\[([+-]*[0-9\.]+)\]";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        // These columns correspond to the tab-delimited file created directly by MSAlign
        protected const int MSAlignResultsFileColCount = 23;
        public enum eMSAlignResultsFileColumns : int
        {
            SpectrumFileName = 0,
            Prsm_ID = 1,
            Spectrum_ID = 2,
            Protein_Sequence_ID = 3,         // Used by MSAlign v0.5, but not by v0.6
            Scans = 4,
            Peaks = 5,
            Charge = 6,
            Precursor_mass = 7,              // Monoisotopic mass value of the observed precursor_mz
            Adjusted_precursor_mass = 8,     // Theoretical monoisotopic mass of the peptide (including mods)
            Protein_ID = 9,
            Protein_name = 10,               // Protein name and description
            Protein_mass = 11,
            First_residue = 12,
            Last_residue = 13,
            Peptide = 14,
            Unexpected_modifications = 15,
            Matched_peaks = 16,
            Matched_fragment_ions = 17,
            Pvalue = 18,
            Evalue = 19,
            FDR = 20,
            Species_ID = 21,             // Present between Protein_ID and Protein_name in MSAlign_Histone result files
            FragMethod = 22              // Present as the last column in MSAlign_Histone result files
        }

        // These columns correspond to the Synopsis file created by this class
        protected const int MSAlignSynFileColCount = 22;
        public enum eMSAlignSynFileColumns : int
        {
            ResultID = 0,
            Scan = 1,
            Prsm_ID = 2,
            Spectrum_ID = 3,
            Charge = 4,
            PrecursorMZ = 5,
            DelM = 6,                            // Precursor error, in Da
            DelMPPM = 7,                         // Precursor error, in ppm
            MH = 8,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
            Peptide = 9,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
            Protein = 10,                        // Protein Name
            Protein_Mass = 11,
            Unexpected_Mod_Count = 12,
            Peak_Count = 13,
            Matched_Peak_Count = 14,
            Matched_Fragment_Ion_Count = 15,
            PValue = 16,
            Rank_PValue = 17,
            EValue = 18,
            FDR = 19,
            Species_ID = 20,
            FragMethod = 21
        }

        #endregion

        #region "Structures"
        // This data structure holds rows read from the tab-delimited file created directly by MSAlign
        protected struct udtMSAlignSearchResultType
        {
            public string SpectrumFileName;
            public string Scans;
            public int ScanNum;
            public string Prsm_ID;
            public string Spectrum_ID;
            public string Peaks;
            public string Charge;
            public short ChargeNum;
            public string Precursor_mass;               // Monoisotopic mass value of the observed precursor_mz
            public string PrecursorMZ;                  // Computed from Precursor_mass
            public string Adjusted_precursor_mass;      // Theoretical monoisotopic mass of the peptide (including mods), as computed by MSAlign
            public string MH;                           // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                         // Computed using Precursor_mass - Adjusted_precursor_mass
            public string DelM_PPM;                     // Computed using DelM and Adjusted_precursor_mass
            public string Protein_ID;
            public string Species_ID;                   // Only present in MSAlign_Histone results
            public string Protein;
            public string Protein_mass;
            public string First_residue;
            public string Last_residue;
            public string Peptide;
            public string Unexpected_modifications;
            public string Matched_peaks;
            public string Matched_fragment_ions;
            public string Pvalue;
            public double PValueNum;
            public string Evalue;
            public string FDR;
            public string FragMethod;                   // Only present in MSAlign_Histone results
            public int RankPValue;

            public void Clear()
            {
                SpectrumFileName = string.Empty;
                Prsm_ID = string.Empty;
                Spectrum_ID = string.Empty;
                Scans = string.Empty;
                ScanNum = 0;
                Peaks = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                Precursor_mass = string.Empty;
                Adjusted_precursor_mass = string.Empty;
                MH = string.Empty;
                DelM = string.Empty;
                DelM_PPM = string.Empty;
                Protein_ID = string.Empty;
                Species_ID = string.Empty;
                Protein = string.Empty;
                Protein_mass = string.Empty;
                First_residue = string.Empty;
                Last_residue = string.Empty;
                Peptide = string.Empty;
                Unexpected_modifications = string.Empty;
                Matched_peaks = string.Empty;
                Matched_fragment_ions = string.Empty;
                Pvalue = string.Empty;
                PValueNum = 0;
                Evalue = string.Empty;
                FDR = string.Empty;
                FragMethod = string.Empty;
                RankPValue = 0;
            }
        }

        #endregion

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to objSearchResult
        /// </summary>
        /// <param name="objSearchResult"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsMSAlign objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const char NO_RESIDUE = '-';

            int intModIndex = 0;
            char chChar = default(char);
            clsModificationDefinition objModificationDefinition = default(clsModificationDefinition);

            string strSequence = null;

            bool blnParsingModMass = false;

            char chMostRecentResidue = default(char);
            int intResidueLocInPeptide = 0;

            char chAmbiguousResidue = default(char);
            int intAmbiguousResidueLocInPeptide = 0;

            bool blnClearAmbiguousResidue = false;
            bool blnStoreAmbiguousResidue = false;

            bool blnSuccess = false;

            blnParsingModMass = false;
            var strModMassDigits = string.Empty;

            chMostRecentResidue = NO_RESIDUE;
            intResidueLocInPeptide = 0;

            chAmbiguousResidue = NO_RESIDUE;
            intAmbiguousResidueLocInPeptide = 0;

            blnClearAmbiguousResidue = false;
            blnStoreAmbiguousResidue = false;

            strSequence = objSearchResult.PeptideSequenceWithMods;
            for (var intIndex = 0; intIndex <= strSequence.Length - 1; intIndex++)
            {
                chChar = strSequence[intIndex];

                if (IsLetterAtoZ(chChar))
                {
                    chMostRecentResidue = chChar;
                    intResidueLocInPeptide += 1;

                    if (blnStoreAmbiguousResidue)
                    {
                        chAmbiguousResidue = chChar;
                        intAmbiguousResidueLocInPeptide = intResidueLocInPeptide;
                        blnStoreAmbiguousResidue = false;
                    }
                    else if (blnClearAmbiguousResidue)
                    {
                        chAmbiguousResidue = NO_RESIDUE;
                        blnClearAmbiguousResidue = false;
                    }

                    for (intModIndex = 0; intModIndex <= mPeptideMods.ModificationCount - 1; intModIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(intModIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);

                            if (objModificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                objSearchResult.SearchResultAddModification(ref objModificationDefinition, chChar, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (chChar == '(')
                {
                    // Start of a mod group
                    blnStoreAmbiguousResidue = true;
                }
                else if (chChar == ')')
                {
                    // End of a mod group
                    blnClearAmbiguousResidue = true;
                }
                else if (chChar == '[')
                {
                    // Mod Mass Start
                    strModMassDigits = string.Empty;
                    blnParsingModMass = true;
                }
                else if (chChar == ']')
                {
                    // Mod Mass End

                    if (blnParsingModMass)
                    {
                        char chResidueForMod = default(char);
                        int intResidueLocForMod = 0;
                        double dblModMass = 0;

                        if (chAmbiguousResidue == NO_RESIDUE)
                        {
                            chResidueForMod = chMostRecentResidue;
                            intResidueLocForMod = intResidueLocInPeptide;
                        }
                        else
                        {
                            // Ambigous mod
                            // We'll associate it with the first residue of the mod group
                            chResidueForMod = chAmbiguousResidue;
                            intResidueLocForMod = intAmbiguousResidueLocInPeptide;
                        }

                        if (double.TryParse(strModMassDigits, out dblModMass))
                        {
                            if (intResidueLocForMod == 0)
                            {
                                // Modification is at the peptide N-terminus
                                intResidueLocForMod = 1;
                            }

                            blnSuccess = objSearchResult.SearchResultAddModification(dblModMass, chResidueForMod, intResidueLocForMod, objSearchResult.DetermineResidueTerminusState(intResidueLocForMod), blnUpdateModOccurrenceCounts);

                            if (!blnSuccess)
                            {
                                string strErrorMessage = objSearchResult.ErrorMessage;
                                if (string.IsNullOrEmpty(strErrorMessage))
                                {
                                    strErrorMessage = "SearchResultAddDynamicModification returned false for mod mass " + strModMassDigits;
                                }
                                SetErrorMessage(strErrorMessage + "; ResultID = " + objSearchResult.ResultID);
                            }
                        }

                        blnParsingModMass = false;
                    }
                }
                else if (blnParsingModMass)
                {
                    strModMassDigits += chChar;
                }
                else
                {
                    // Unrecognized symbol; ignore it
                }
            }
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsMSAlign objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool blnSuccess = false;

            try
            {
                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts);

                // Parse .PeptideSequenceWithMods to determine the modified residues present
                AddDynamicAndStaticResidueMods(objSearchResult, blnUpdateModOccurrenceCounts);

                // Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
                // Since Inspect allows a terminal peptide residue to be modified twice, we'll allow that to happen,
                //  even though, biologically, that's typically not possible
                // However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus
                //  (where two COOH groups are present)
                objSearchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts);

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

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="lstSearchResults"></param>
        /// <param name="intStartIndex"></param>
        /// <param name="intEndIndex"></param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            List<udtMSAlignSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of lstSearchResults so that we can sort by PValue

            var dctResultsSubset = new Dictionary<int, udtMSAlignSearchResultType>();
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                dctResultsSubset.Add(intIndex, lstSearchResults[intIndex]);
            }

            var lstResultsByProbability = (from item in dctResultsSubset orderby item.Value.PValueNum select item).ToList();

            double dblLastValue = 0;
            int intCurrentRank = -1;

            foreach (var entry in lstResultsByProbability)
            {
                var oResult = lstSearchResults[entry.Key];

                if (intCurrentRank < 0)
                {
                    dblLastValue = oResult.PValueNum;
                    intCurrentRank = 1;
                }
                else
                {
                    if (Math.Abs(oResult.PValueNum - dblLastValue) > double.Epsilon)
                    {
                        dblLastValue = oResult.PValueNum;
                        intCurrentRank += 1;
                    }
                }

                oResult.RankPValue = intCurrentRank;
                lstSearchResults[entry.Key] = oResult;
            }
        }

        protected string AssureInteger(string strInteger, int intDefaultValue)
        {
            int intValue = 0;
            double dblValue = 0;

            if (strInteger.EndsWith(".0"))
                strInteger = strInteger.Substring(0, strInteger.Length - 2);

            if (int.TryParse(strInteger, out intValue))
            {
                return intValue.ToString();
            }
            else if (double.TryParse(strInteger, out dblValue))
            {
                return dblValue.ToString("0");
            }
            else
            {
                return intDefaultValue.ToString();
            }
        }

        protected double ComputePeptideMass(string strPeptide, double dblTotalModMass)
        {
            string strCleanSequence = null;
            double dblMass = 0;

            strCleanSequence = GetCleanSequence(strPeptide);

            dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence);

            if (Math.Abs(dblTotalModMass) > double.Epsilon)
            {
                dblMass += dblTotalModMass;
            }

            return dblMass;
        }

        private static readonly Regex RegexModMassRegEx = new Regex(MSALIGN_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="strPeptide">Peptide sequence, with mod masses in the form [23.5432]</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected double ComputeTotalModMass(string strPeptide)
        {
            double dblTotalModMass = 0;
            double dblModMassFound = 0;

            dblTotalModMass = 0;

            foreach (Match reMatch in RegexModMassRegEx.Matches(strPeptide))
            {
                if (double.TryParse(reMatch.Groups[1].Value, out dblModMassFound))
                {
                    dblTotalModMass += dblModMassFound;
                }
            }

            return dblTotalModMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string strInputFilePath, string strOutputFolderPath, bool MTS)
        {
            var strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
            if (strPepToProteinMapFilePath.EndsWith("_msalign_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_msalign_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MSAlign
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreateSynResultsFile(
            string strInputFilePath,
            string strOutputFilePath)
        {
            string strLineIn = null;

            float sngPercentComplete = 0;

            bool blnHeaderParsed = false;

            int[] intColumnMapping = null;

            int intResultsProcessed = 0;

            bool blnSuccess = false;
            bool blnValidSearchResult = false;

            try
            {
                string strErrorLog = string.Empty;

                // Open the input file and parse it
                // Initialize the stream reader and the stream Text writer
                using (var srDataFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    using (var swResultFile = new StreamWriter(new FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        intResultsProcessed = 0;
                        blnHeaderParsed = false;

                        // Initialize array that will hold all of the records in the MSAlign result file
                        var lstSearchResultsUnfiltered = new List<udtMSAlignSearchResultType>();

                        // Initialize the array that will hold all of the records that will ultimately be written out to disk
                        var lstFilteredSearchResults = new List<udtMSAlignSearchResultType>();

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !base.AbortProcessing)
                        {
                            strLineIn = srDataFile.ReadLine();
                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseMSAlignResultsFileHeaderLine(strLineIn, ref intColumnMapping);
                                if (!blnSuccess)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                blnHeaderParsed = true;

                                // Write the header line
                                WriteSynFHTFileHeader(swResultFile, ref strErrorLog);

                                continue;
                            }

                            var udtSearchResult = new udtMSAlignSearchResultType();
                            blnValidSearchResult = ParseMSAlignResultsFileEntry(ref strLineIn, ref udtSearchResult, ref strErrorLog, ref intColumnMapping);

                            if (blnValidSearchResult)
                            {
                                lstSearchResultsUnfiltered.Add(udtSearchResult);
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

                        // Sort the SearchResults by scan, charge, and ascending PValue
                        lstSearchResultsUnfiltered.Sort(new MSAlignSearchResultsComparerScanChargePValuePeptide());

                        // Now filter the data

                        // Initialize variables
                        var intStartIndex = 0;

                        while (intStartIndex < lstSearchResultsUnfiltered.Count)
                        {
                            var intEndIndex = intStartIndex;
                            while (intEndIndex + 1 < lstSearchResultsUnfiltered.Count && lstSearchResultsUnfiltered[intEndIndex + 1].ScanNum == lstSearchResultsUnfiltered[intStartIndex].ScanNum)
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

        protected bool ExtractModInfoFromMSAlignParamFile(string strMSAlignParamFilePath, ref List<clsModificationDefinition> lstModInfo)
        {
            string strLineIn = null;
            KeyValuePair<string, string> kvSetting = default(KeyValuePair<string, string>);

            clsModificationDefinition objModDef = default(clsModificationDefinition);

            var blnSuccess = false;

            try
            {
                // Initialize the modification list
                if (lstModInfo == null)
                {
                    lstModInfo = new List<clsModificationDefinition>();
                }
                else
                {
                    lstModInfo.Clear();
                }

                if (string.IsNullOrEmpty(strMSAlignParamFilePath))
                {
                    SetErrorMessage("MSAlign Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(strMSAlignParamFilePath))
                {
                    SetErrorMessage("MSAlign param file not found: " + strMSAlignParamFilePath);
                }
                else
                {
                    // Read the contents of the parameter (or mods) file
                    using (var srInFile = new StreamReader(new FileStream(strMSAlignParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!srInFile.EndOfStream)
                        {
                            strLineIn = srInFile.ReadLine().Trim();

                            if (strLineIn.Length > 0)
                            {
                                if (strLineIn.StartsWith("#"))
                                {
                                    // Comment line; skip it
                                }
                                else
                                {
                                    // Split the line on the equals sign
                                    kvSetting = clsPHRPParser.ParseKeyValueSetting(strLineIn, '=', "#");

                                    if (kvSetting.Key.ToLower() == "cysteineProtection".ToLower())
                                    {
                                        switch (kvSetting.Value.ToUpper())
                                        {
                                            case "C57":
                                                objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 57.0215, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcet");
                                                lstModInfo.Add(objModDef);
                                                break;

                                            case "C58":
                                                objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 58.0055, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcid");
                                                lstModInfo.Add(objModDef);
                                                break;

                                        }
                                    }
                                }
                            }
                        }
                    }

                    Console.WriteLine();

                    blnSuccess = true;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MSAlign parameter file (" + Path.GetFileName(strMSAlignParamFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private void InitializeLocalVariables()
        {
            // Nothing to do at present
        }

        protected bool ParseMSAlignSynopsisFile(string strInputFilePath, string strOutputFolderPath, ref List<udtPepToProteinMappingType> lstPepToProteinMapping, bool blnResetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            string strPreviousPValue = null;

            // Note that MSAlign synopsis files are normally sorted on PValue value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique PValue encountered
            // Although this was a possiblity with Inspect, it likely never occurs for MSAlign
            //  But, we'll keep the check in place just in case

            Hashtable htPeptidesFoundForPValueLevel = default(Hashtable);

            string strKey = null;

            string strLineIn = null;
            string strModificationSummaryFilePath = null;

            clsSearchResultsMSAlign objSearchResult = default(clsSearchResultsMSAlign);

            int intResultsProcessed = 0;
            int intPepToProteinMapIndex = 0;
            float sngPercentComplete = 0;

            bool blnHeaderParsed = false;
            int[] intColumnMapping = null;

            string strCurrentPeptideWithMods = string.Empty;
            string strCurrentProtein = null;

            bool blnSuccess = false;
            bool blnValidSearchResult = false;
            bool blnFirstMatchForGroup = false;

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
                objSearchResult = new clsSearchResultsMSAlign(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForPValueLevel
                htPeptidesFoundForPValueLevel = new Hashtable();
                strPreviousPValue = string.Empty;

                // Assure that lstPepToProteinMapping is sorted on peptide
                if (lstPepToProteinMapping.Count > 1)
                {
                    lstPepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                try
                {
                    objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    string strErrorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var srDataFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
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

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseMSAlignSynFileHeaderLine(strLineIn, ref intColumnMapping);
                                if (!blnSuccess)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return blnSuccess;
                                }
                                blnHeaderParsed = true;
                                continue;
                            }

                            blnValidSearchResult = ParseMSAlignSynFileEntry(ref strLineIn, objSearchResult, ref strErrorLog,
                                                                            intResultsProcessed, ref intColumnMapping,
                                                                            ref strCurrentPeptideWithMods);

                            if (blnValidSearchResult)
                            {
                                strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.Charge;

                                if (objSearchResult.PValue == strPreviousPValue)
                                {
                                    // New result has the same PValue as the previous result
                                    // See if htPeptidesFoundForPValueLevel contains the peptide, scan and charge

                                    if (htPeptidesFoundForPValueLevel.ContainsKey(strKey))
                                    {
                                        blnFirstMatchForGroup = false;
                                    }
                                    else
                                    {
                                        htPeptidesFoundForPValueLevel.Add(strKey, 1);
                                        blnFirstMatchForGroup = true;
                                    }
                                }
                                else
                                {
                                    // New PValue
                                    // Reset htPeptidesFoundForPValueLevel
                                    htPeptidesFoundForPValueLevel.Clear();

                                    // Update strPreviousPValue
                                    strPreviousPValue = objSearchResult.PValue;

                                    // Append a new entry to htPeptidesFoundForPValueLevel
                                    htPeptidesFoundForPValueLevel.Add(strKey, 1);
                                    blnFirstMatchForGroup = true;
                                }

                                blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup);
                                if (!blnSuccess)
                                {
                                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        strErrorLog += "Error adding modifications to sequence at RowIndex '" + objSearchResult.ResultID + "'" +
                                                       "\n";
                                    }
                                }

                                base.SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass) objSearchResult, blnFirstMatchForGroup);

                                if (lstPepToProteinMapping.Count > 0)
                                {
                                    // Add the additional proteins for this peptide

                                    // Use binary search to find this peptide in lstPepToProteinMapping
                                    intPepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(lstPepToProteinMapping, strCurrentPeptideWithMods);

                                    if (intPepToProteinMapIndex >= 0)
                                    {
                                        // Call MyBase.SaveResultsFileEntrySeqInfo for each entry in lstPepToProteinMapping() for peptide , skipping objSearchResult.ProteinName
                                        strCurrentProtein = string.Copy(objSearchResult.ProteinName);
                                        do
                                        {
                                            if (lstPepToProteinMapping[intPepToProteinMapIndex].Protein != strCurrentProtein)
                                            {
                                                objSearchResult.ProteinName = string.Copy(lstPepToProteinMapping[intPepToProteinMapIndex].Protein);
                                                base.SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass) objSearchResult, false);
                                            }

                                            intPepToProteinMapIndex += 1;
                                        } while (intPepToProteinMapIndex < lstPepToProteinMapping.Count && strCurrentPeptideWithMods == lstPepToProteinMapping[intPepToProteinMapIndex].Peptide);
                                    }
                                    else
                                    {
                                        // Match not found; this is unexpected
                                        ReportWarning("no match for '" + strCurrentPeptideWithMods + "' in lstPepToProteinMapping");
                                    }
                                }
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

        private bool ParseMSAlignResultsFileEntry(
            ref string strLineIn,
            ref udtMSAlignSearchResultType udtSearchResult,
            ref string strErrorLog,
            ref int[] intColumnMapping)
        {
            // Parses an entry from the MSAlign results file

            string[] strSplitLine = null;

            double dblPrecursorMonoMass = 0;           // Observed m/z, converted to monoisotopic mass
            double dblPeptideMonoMassMSAlign = 0;      // Theoretical peptide monoisotopic mass, including mods, as computed by MSAlign
            double dblPeptideMonoMassPHRP = 0;         // Theoretical peptide monoisotopic mass, including mods, as computed by PHRP

            double dblPrecursorMZ = 0;
            double dblDelM = 0;
            double dblFDR = 0;

            double dblTotalModMass = 0;

            bool blnValidSearchResult = false;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                udtSearchResult.Clear();
                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 13)
                {
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);
                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Prsm_ID], out udtSearchResult.Prsm_ID))
                    {
                        ReportError("Prsm_ID column is missing or invalid", true);
                    }
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Spectrum_ID], out udtSearchResult.Spectrum_ID);

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Scans], out udtSearchResult.Scans))
                    {
                        ReportError("Scan(s) column is missing or invalid", true);
                    }

                    if (!int.TryParse(udtSearchResult.Scans, out udtSearchResult.ScanNum))
                    {
                        // .Scans likely has a list of scan numbers; extract the first scan number from .scans
                        string strScanNumberDigits = string.Empty;
                        foreach (char chChar in udtSearchResult.Scans)
                        {
                            if (char.IsDigit(chChar))
                            {
                                strScanNumberDigits += chChar;
                            }
                        }
                        if (!int.TryParse(strScanNumberDigits, out udtSearchResult.ScanNum))
                        {
                            ReportWarning("Error parsing out the scan number from the scan list; could not find an integer: " + udtSearchResult.Scans);
                            udtSearchResult.ScanNum = 0;
                        }
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Peaks], out udtSearchResult.Peaks);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    // Monoisotopic mass value of the observed precursor_mz
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Precursor_mass], out udtSearchResult.Precursor_mass);

                    if (double.TryParse(udtSearchResult.Precursor_mass, out dblPrecursorMonoMass))
                    {
                        if (udtSearchResult.ChargeNum > 0)
                        {
                            dblPrecursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMonoMass, 0, udtSearchResult.ChargeNum);
                            udtSearchResult.PrecursorMZ = NumToString(dblPrecursorMZ, 6, true);
                        }
                    }

                    if (intColumnMapping[(int)eMSAlignResultsFileColumns.Adjusted_precursor_mass] >= 0)
                    {
                        // Theoretical monoisotopic mass of the peptide (including mods), as computed by MSAlign
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Adjusted_precursor_mass], out udtSearchResult.Adjusted_precursor_mass);

                        double.TryParse(udtSearchResult.Adjusted_precursor_mass, out dblPeptideMonoMassMSAlign);
                    }
                    else
                    {
                        dblPeptideMonoMassMSAlign = 0;
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Protein_ID], out udtSearchResult.Protein_ID);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Species_ID], out udtSearchResult.Species_ID);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Protein_name], out udtSearchResult.Protein);
                    udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Protein_mass], out udtSearchResult.Protein_mass);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.First_residue], out udtSearchResult.First_residue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Last_residue], out udtSearchResult.Last_residue);

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                    {
                        ReportError("Peptide column is missing or invalid", true);
                    }

                    // Add the standard terminus symbols to the peptide sequence
                    udtSearchResult.Peptide = ReplaceTerminus(udtSearchResult.Peptide);

                    // Parse the sequence to determine the total mod mass
                    // Note that we do not remove any of the mod symbols since MSAlign identifies mods by mass alone, and since mods can ambiguously apply to residues
                    dblTotalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                    // Compute monoisotopic mass of the peptide
                    dblPeptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, dblTotalModMass);

                    if (Math.Abs(dblPeptideMonoMassMSAlign) < double.Epsilon)
                    {
                        dblPeptideMonoMassMSAlign = dblPeptideMonoMassPHRP;
                    }

                    double dblMassDiffThreshold = dblPeptideMonoMassMSAlign / 50000;
                    if (dblMassDiffThreshold < 0.1)
                        dblMassDiffThreshold = 0.1;

                    if (Math.Abs(dblPeptideMonoMassPHRP - dblPeptideMonoMassMSAlign) > dblMassDiffThreshold)
                    {
                        // Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da or by a slightly larger value if over 5000 Da; this is unexpected
                        string strFirst30Residues = null;
                        if (udtSearchResult.Peptide.Length < 27)
                        {
                            strFirst30Residues = udtSearchResult.Peptide;
                        }
                        else
                        {
                            strFirst30Residues = udtSearchResult.Peptide.Substring(0, 27) + "...";
                        }
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " + dblMassDiffThreshold.ToString("0.00") + " Da away from the mass computed by MSAlign: " + dblPeptideMonoMassPHRP.ToString("0.0000") + " vs. " + dblPeptideMonoMassMSAlign.ToString("0.0000") + "; peptide " + strFirst30Residues);
                    }

                    if (dblPeptideMonoMassMSAlign > 0)
                    {
                        // Compute DelM and DelM_PPM
                        dblDelM = dblPrecursorMonoMass - dblPeptideMonoMassMSAlign;
                        udtSearchResult.DelM = NumToString(dblDelM, 6, true);

                        if (dblPrecursorMZ > 0)
                        {
                            udtSearchResult.DelM_PPM = NumToString(clsPeptideMassCalculator.MassToPPM(dblDelM, dblPrecursorMZ), 5, true);
                        }
                        else
                        {
                            udtSearchResult.DelM_PPM = NumToString(clsPeptideMassCalculator.MassToPPM(dblDelM, 1000), 5, true);
                        }
                    }

                    // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    udtSearchResult.MH = NumToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPeptideMonoMassPHRP, 0, 1), 6, true);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Unexpected_modifications], out udtSearchResult.Unexpected_modifications);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Matched_peaks], out udtSearchResult.Matched_peaks);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Matched_fragment_ions], out udtSearchResult.Matched_fragment_ions);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Pvalue], out udtSearchResult.Pvalue);
                    if (!double.TryParse(udtSearchResult.Pvalue, out udtSearchResult.PValueNum))
                        udtSearchResult.PValueNum = 0;

                    // Assure that the following are truly integers (Matched_peaks and Matched_fragment_ions are often of the form 8.0)
                    udtSearchResult.Unexpected_modifications = AssureInteger(udtSearchResult.Unexpected_modifications, 0); // Unexpected_Mod_Count
                    udtSearchResult.Peaks = AssureInteger(udtSearchResult.Peaks, 0);                                       // Peak_count
                    udtSearchResult.Matched_peaks = AssureInteger(udtSearchResult.Matched_peaks, 0);                       // Matched_Peak_Count
                    udtSearchResult.Matched_fragment_ions = AssureInteger(udtSearchResult.Matched_fragment_ions, 0);       // Matched_Fragment_Ion_Count

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.Evalue], out udtSearchResult.Evalue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.FDR], out udtSearchResult.FDR);

                    if (udtSearchResult.FDR.ToLower() == "infinity")
                    {
                        udtSearchResult.FDR = "10";
                    }
                    else if (!string.IsNullOrEmpty(udtSearchResult.FDR) & !double.TryParse(udtSearchResult.FDR, out dblFDR))
                    {
                        udtSearchResult.FDR = "";
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MSAlign results file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if ((strSplitLine != null) && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing MSAlign Results in ParseMSAlignResultsFileEntry for RowIndex '" + strSplitLine[0] + "'" +
                                       "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MSAlign Results in ParseMSAlignResultsFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        private bool ParseMSAlignResultsFileHeaderLine(string strLineIn, ref int[] intColumnMapping)
        {
            // Parse the header line

            // The expected header from MSAlign v0.5 is:
            //                   Prsm_ID    Spectrum_ID    Protein_Sequence_ID    Spectrum_ID    Scan(s)    #peaks    Charge    Precursor_mass                                                           Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions               E-value

            // The expected header from MSAlign v0.6 is:
            // Data_file_name    Prsm_ID    Spectrum_ID                                          Scan(s)    #peaks    Charge    Precursor_mass    Adjusted_precursor_mass    Protein_ID                  Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions    P-value    E-value    FDR

            // The expected header from MSAlign_Histone v0.9 is
            // Data_file_name    Prsm_ID    Spectrum_ID                                          Scan(s)    #peaks    Charge    Precursor_mass    Adjusted_precursor_mass    Protein_ID    Species_ID    Protein_name    Protein_mass    First_residue    Last_residue    Peptide    #unexpected_modifications    #matched_peaks    #matched_fragment_ions    P-value    E-value    FDR    FragMethod

            string[] strSplitLine = null;
            eMSAlignResultsFileColumns eResultFileColumn = default(eMSAlignResultsFileColumns);
            SortedDictionary<string, eMSAlignResultsFileColumns> lstColumnNames = default(SortedDictionary<string, eMSAlignResultsFileColumns>);
            lstColumnNames = new SortedDictionary<string, eMSAlignResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase);

            intColumnMapping = new int[MSAlignResultsFileColCount];

            lstColumnNames.Add("Data_file_name", eMSAlignResultsFileColumns.SpectrumFileName);
            lstColumnNames.Add("Prsm_ID", eMSAlignResultsFileColumns.Prsm_ID);
            lstColumnNames.Add("Spectrum_ID", eMSAlignResultsFileColumns.Spectrum_ID);
            lstColumnNames.Add("Protein_Sequence_ID", eMSAlignResultsFileColumns.Protein_Sequence_ID);
            lstColumnNames.Add("Scan(s)", eMSAlignResultsFileColumns.Scans);
            lstColumnNames.Add("#peaks", eMSAlignResultsFileColumns.Peaks);
            lstColumnNames.Add("Charge", eMSAlignResultsFileColumns.Charge);
            lstColumnNames.Add("Precursor_mass", eMSAlignResultsFileColumns.Precursor_mass);
            lstColumnNames.Add("Adjusted_precursor_mass", eMSAlignResultsFileColumns.Adjusted_precursor_mass);
            lstColumnNames.Add("Protein_ID", eMSAlignResultsFileColumns.Protein_ID);
            lstColumnNames.Add("Species_ID", eMSAlignResultsFileColumns.Species_ID);
            lstColumnNames.Add("Protein_name", eMSAlignResultsFileColumns.Protein_name);
            lstColumnNames.Add("Protein_mass", eMSAlignResultsFileColumns.Protein_mass);
            lstColumnNames.Add("First_residue", eMSAlignResultsFileColumns.First_residue);
            lstColumnNames.Add("Last_residue", eMSAlignResultsFileColumns.Last_residue);
            lstColumnNames.Add("Peptide", eMSAlignResultsFileColumns.Peptide);
            lstColumnNames.Add("#unexpected_modifications", eMSAlignResultsFileColumns.Unexpected_modifications);
            lstColumnNames.Add("#matched_peaks", eMSAlignResultsFileColumns.Matched_peaks);
            lstColumnNames.Add("#matched_fragment_ions", eMSAlignResultsFileColumns.Matched_fragment_ions);
            lstColumnNames.Add("P-value", eMSAlignResultsFileColumns.Pvalue);
            lstColumnNames.Add("E-value", eMSAlignResultsFileColumns.Evalue);
            lstColumnNames.Add("FDR", eMSAlignResultsFileColumns.FDR);
            lstColumnNames.Add("FragMethod", eMSAlignResultsFileColumns.FragMethod);

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
                    else
                    {
                        // Unrecognized column name
                        Console.WriteLine(
                            "Warning: Unrecognized column header name '" + strSplitLine[intIndex] + "' in ParseMSAlignResultsFileHeaderLine");
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MSAlign results file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMSAlignSynFileHeaderLine(string strLineIn, ref int[] intColumnMapping)
        {
            // Parse the header line

            string[] strSplitLine = null;
            eMSAlignSynFileColumns eResultFileColumn = default(eMSAlignSynFileColumns);
            SortedDictionary<string, eMSAlignSynFileColumns> lstColumnNames = default(SortedDictionary<string, eMSAlignSynFileColumns>);
            lstColumnNames = new SortedDictionary<string, eMSAlignSynFileColumns>(StringComparer.InvariantCultureIgnoreCase);

            intColumnMapping = new int[MSAlignSynFileColCount];

            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_ResultID, eMSAlignSynFileColumns.ResultID);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Scan, eMSAlignSynFileColumns.Scan);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Prsm_ID, eMSAlignSynFileColumns.Prsm_ID);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Spectrum_ID, eMSAlignSynFileColumns.Spectrum_ID);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Charge, eMSAlignSynFileColumns.Charge);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_PrecursorMZ, eMSAlignSynFileColumns.PrecursorMZ);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_DelM, eMSAlignSynFileColumns.DelM);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_DelM_PPM, eMSAlignSynFileColumns.DelMPPM);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_MH, eMSAlignSynFileColumns.MH);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Peptide, eMSAlignSynFileColumns.Peptide);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Protein, eMSAlignSynFileColumns.Protein);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Protein_Mass, eMSAlignSynFileColumns.Protein_Mass);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Unexpected_Mod_Count, eMSAlignSynFileColumns.Unexpected_Mod_Count);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Peak_Count, eMSAlignSynFileColumns.Peak_Count);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Matched_Peak_Count, eMSAlignSynFileColumns.Matched_Peak_Count);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Matched_Fragment_Ion_Count, eMSAlignSynFileColumns.Matched_Fragment_Ion_Count);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_PValue, eMSAlignSynFileColumns.PValue);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Rank_PValue, eMSAlignSynFileColumns.Rank_PValue);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_EValue, eMSAlignSynFileColumns.EValue);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_FDR, eMSAlignSynFileColumns.FDR);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_Species_ID, eMSAlignSynFileColumns.Species_ID);
            lstColumnNames.Add(clsPHRPParserMSAlign.DATA_COLUMN_FragMethod, eMSAlignSynFileColumns.FragMethod);

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
                SetErrorMessage("Error parsing header in MSAlign synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMSAlignSynFileEntry(
            ref string strLineIn,
            clsSearchResultsMSAlign objSearchResult,
            ref string strErrorLog,
            int intResultsProcessed,
            ref int[] intColumnMapping,
            ref string strPeptideSequenceWithMods)
        {
            // Parses an entry from the MSAlign Synopsis file

            string[] strSplitLine = null;

            bool blnValidSearchResult = false;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                // Reset objSearchResult
                objSearchResult.Clear();
                strPeptideSequenceWithMods = string.Empty;

                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 15)
                {
                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.ResultID], out string strValue))
                    {
                        if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                        {
                            strErrorLog += "Error reading ResultID value from MSAlign Results line " + (intResultsProcessed + 1).ToString() +
                                           "\n";
                        }
                        return false;
                    }

                    objSearchResult.ResultID = int.Parse(strValue);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Scan], out string scan);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Charge], out string charge);

                    objSearchResult.Scan = scan;
                    objSearchResult.Charge = charge;

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Peptide], out strPeptideSequenceWithMods))
                    {
                        if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                        {
                            strErrorLog += "Error reading Peptide sequence value from MSAlign Results line " + (intResultsProcessed + 1).ToString() +
                                           "\n";
                        }
                        return false;
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Protein], out string proteinName);
                    objSearchResult.MultipleProteinCount = "0";

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.DelM], out string msAlignComputedDelM);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.DelMPPM], out string msAlignComputedDelMppm);

                    objSearchResult.ProteinName = proteinName;
                    objSearchResult.MSAlignComputedDelM = msAlignComputedDelM;
                    objSearchResult.MSAlignComputedDelMPPM = msAlignComputedDelMppm;

                    objSearchResult.PeptideDeltaMass = objSearchResult.MSAlignComputedDelM;

                    // Note: .PeptideDeltaMass is stored in the MSAlign results file as "Observed_Mass - Theoretical_Mass"
                    // However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                    // Therefore, we will negate .peptideDeltaMass
                    try
                    {
                        objSearchResult.PeptideDeltaMass = (-double.Parse(objSearchResult.PeptideDeltaMass)).ToString();
                    }
                    catch (Exception)
                    {
                        // Error; Leave .peptideDeltaMass unchanged
                    }

                    // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                    objSearchResult.SetPeptideSequenceWithMods(strPeptideSequenceWithMods, true, true);

                    clsSearchResultsBaseClass objSearchResultBase = default(clsSearchResultsBaseClass);
                    objSearchResultBase = (clsSearchResultsBaseClass) objSearchResult;

                    base.ComputePseudoPeptideLocInProtein(objSearchResultBase);

                    // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                    // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                    // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                    objSearchResult.ComputePeptideCleavageStateInProtein();

                    // Read the remaining data values
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Prsm_ID], out string prsmId);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Spectrum_ID], out string spectrumId);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.PrecursorMZ], out string precursorMz);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.MH], out string parentIonMH);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Protein_Mass], out string proteinMass);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Unexpected_Mod_Count], out string unexpectedModCount);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Peak_Count], out string peakCount);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Matched_Peak_Count], out string matchedPeakCount);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Matched_Fragment_Ion_Count], out string matchedFragmentIonCount);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.PValue], out string pValue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Rank_PValue], out string rankPValue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.EValue], out string eValue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.FDR], out string fdr);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.Species_ID], out string speciesId);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSAlignSynFileColumns.FragMethod], out string fragMethod);

                    objSearchResult.Prsm_ID = prsmId;
                    objSearchResult.Spectrum_ID = spectrumId;
                    objSearchResult.Precursor_mz = precursorMz;
                    objSearchResult.ParentIonMH = parentIonMH;
                    objSearchResult.Protein_Mass = proteinMass;
                    objSearchResult.Unexpected_Mod_Count = unexpectedModCount;
                    objSearchResult.Peak_Count = peakCount;
                    objSearchResult.Matched_Peak_Count = matchedPeakCount;
                    objSearchResult.Matched_Fragment_Ion_Count = matchedFragmentIonCount;
                    objSearchResult.PValue = pValue;
                    objSearchResult.Rank_PValue = rankPValue;
                    objSearchResult.EValue = eValue;
                    objSearchResult.FDR = fdr;
                    objSearchResult.Species_ID = speciesId;
                    objSearchResult.FragMethod = fragMethod;

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if ((strSplitLine != null) && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing MSAlign Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MSAlign Results in ParseMSAlignSynFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">MSAlign results file</param>
        /// <param name="strOutputFolderPath">Output folder</param>
        /// <param name="strParameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            string strBaseName = string.Empty;
            string strSynOutputFilePath = string.Empty;

            List<clsModificationDefinition> lstMSAlignModInfo = default(List<clsModificationDefinition>);
            List<udtPepToProteinMappingType> lstPepToProteinMapping = default(List<udtPepToProteinMappingType>);

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

                try
                {
                    // Obtain the full path to the input file
                    var fiInputFile = new FileInfo(strInputFilePath);

                    lstMSAlignModInfo = new List<clsModificationDefinition>();
                    lstPepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the MSAlign Parameter File so that we can determine whether Cysteine residues are statically modified
                    ExtractModInfoFromMSAlignParamFile(SearchToolParameterFilePath, ref lstMSAlignModInfo);

                    // Resolve the mods in lstMSAlignModInfo with the ModDefs mods
                    ResolveMSAlignModsWithModDefinitions(ref lstMSAlignModInfo);

                    // Define the base output filename using strInputFilePath
                    strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath);

                    // Auto-replace "_MSAlign_ResultTable" with "_msalign"
                    if (strBaseName.EndsWith("_MSAlign_ResultTable", StringComparison.InvariantCultureIgnoreCase))
                    {
                        strBaseName = strBaseName.Substring(0, strBaseName.Length - "_MSAlign_ResultTable".Length) + "_msalign";
                    }

                    // Do not create a first-hits file for MSAlign results

                    // Create the synopsis output file
                    base.ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_msalign_syn.txt
                    strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath);

                    // Create the other PHRP-specific files
                    base.ResetProgress("Creating the PHRP files for " + Path.GetFileName(strSynOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    blnSuccess = ParseMSAlignSynopsisFile(strSynOutputFilePath, strOutputFolderPath, ref lstPepToProteinMapping, false);

                    // Remove all items from lstPepToProteinMapping to reduce memory overhead
                    lstPepToProteinMapping.Clear();
                    lstPepToProteinMapping.TrimExcess();

                    if (blnSuccess && CreateProteinModsFile)
                    {
                        blnSuccess = CreateProteinModsFileWork(strBaseName, fiInputFile, strSynOutputFilePath, strOutputFolderPath);
                    }

                    if (blnSuccess)
                    {
                        base.OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in clsMSAlignResultsProcessor.ProcessFile (2):  " + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile (1):" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return blnSuccess;
        }

        private bool CreateProteinModsFileWork(string strBaseName, FileInfo fiInputFile, string strSynOutputFilePath, string strOutputFolderPath)
        {
            bool blnSuccess = false;
            string strMTSPepToProteinMapFilePath = null;

            // Create the MTSPepToProteinMap file

            strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS: true);

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
                    // Auto-change mIgnorePeptideToProteinMapperErrors to True
                    // We only do this since a small number of peptides reported by MSAlign don't perfectly match the fasta file
                    IgnorePeptideToProteinMapperErrors = true;
                    blnSuccess = base.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                    if (!blnSuccess)
                    {
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                base.ValidatePHRPReaderSupportFiles(Path.Combine(fiInputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath);

                // Create the Protein Mods file
                blnSuccess = base.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSAlign);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                blnSuccess = true;
            }
            return blnSuccess;
        }

        protected string ReplaceTerminus(string strPeptide)
        {
            if (strPeptide.StartsWith(N_TERMINUS_SYMBOL_MSALIGN))
            {
                strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + strPeptide.Substring(N_TERMINUS_SYMBOL_MSALIGN.Length);
            }

            if (strPeptide.EndsWith(C_TERMINUS_SYMBOL_MSALIGN))
            {
                strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_MSALIGN.Length) + "." + clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return strPeptide;
        }

        protected void ResolveMSAlignModsWithModDefinitions(ref List<clsModificationDefinition> lstMSAlignModInfo)
        {
            bool blnExistingModFound = false;
            clsModificationDefinition objModDef = default(clsModificationDefinition);

            if ((lstMSAlignModInfo != null))
            {
                // Call .LookupModificationDefinitionByMass for each entry in lstMSAlignModInfo
                foreach (clsModificationDefinition objModInfo in lstMSAlignModInfo)
                {
                    if (string.IsNullOrEmpty(objModInfo.TargetResidues))
                    {
                        objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, default(char), clsAminoAcidModInfo.eResidueTerminusStateConstants.None,  ref blnExistingModFound, true, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION);
                    }
                    else
                    {
                        foreach (char chTargetResidue in objModInfo.TargetResidues)
                        {
                            objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, chTargetResidue, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, ref blnExistingModFound, true, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION);
                        }
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            StreamWriter swResultFile,
            List<udtMSAlignSearchResultType> lstFilteredSearchResults,
            ref string strErrorLog)
        {
            // Sort udtFilteredSearchResults by ascending PVAlue, ascending scan, ascending charge, ascending peptide, and ascending protein
            lstFilteredSearchResults.Sort(new MSAlignSearchResultsComparerPValueScanChargePeptide());

            for (var intIndex = 0; intIndex <= lstFilteredSearchResults.Count - 1; intIndex++)
            {
                WriteSearchResultToFile(intIndex + 1, swResultFile, lstFilteredSearchResults[intIndex], ref strErrorLog);
            }
        }

        private void StoreSynMatches(
            List<udtMSAlignSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex,
            List<udtMSAlignSearchResultType> lstFilteredSearchResults)
        {
            AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex);

            // The calling procedure already sorted by scan, charge, and PValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                if (lstSearchResults[intIndex].PValueNum <= MSAlignSynopsisFilePValueThreshold)
                {
                    lstFilteredSearchResults.Add(lstSearchResults[intIndex]);
                }
            }
        }

        private void WriteSynFHTFileHeader(
            StreamWriter swResultFile,
            ref string strErrorLog)
        {
            // Write out the header line for synopsis / first hits files
            try
            {
                List<string> lstData = new List<string>();
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_ResultID);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Scan);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Prsm_ID);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Spectrum_ID);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Charge);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_PrecursorMZ);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_DelM);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_DelM_PPM);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_MH);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Peptide);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Protein);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Protein_Mass);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Unexpected_Mod_Count);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Peak_Count);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Matched_Peak_Count);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Matched_Fragment_Ion_Count);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_PValue);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Rank_PValue);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_EValue);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_FDR);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_Species_ID);
                lstData.Add(clsPHRPParserMSAlign.DATA_COLUMN_FragMethod);

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
        /// Writes an entry to the synopsis file
        /// </summary>
        /// <param name="intResultID"></param>
        /// <param name="swResultFile"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="strErrorLog"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int intResultID,
            StreamWriter swResultFile,
            udtMSAlignSearchResultType udtSearchResult,
            ref string strErrorLog)
        {
            try
            {
                // Primary Columns
                //
                // MSAlign
                // ResultID  Scan  Prsm_ID  Spectrum_ID  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  Protein_Mass  Unexpected_Mod_Count  Peak_Count  Matched_Peak_Count  Matched_Fragment_Ion_Count  PValue  Rank_PValue  EValue  FDR

                List<string> lstData = new List<string>();
                lstData.Add(intResultID.ToString());
                lstData.Add(udtSearchResult.ScanNum.ToString());
                lstData.Add(udtSearchResult.Prsm_ID);
                lstData.Add(udtSearchResult.Spectrum_ID);
                lstData.Add(udtSearchResult.Charge);
                lstData.Add(udtSearchResult.PrecursorMZ);
                lstData.Add(udtSearchResult.DelM);
                lstData.Add(udtSearchResult.DelM_PPM);
                lstData.Add(udtSearchResult.MH);
                lstData.Add(udtSearchResult.Peptide);
                lstData.Add(udtSearchResult.Protein);
                lstData.Add(udtSearchResult.Protein_mass);
                lstData.Add(udtSearchResult.Unexpected_modifications);   // Unexpected_Mod_Count
                lstData.Add(udtSearchResult.Peaks);                      // Peak_count
                lstData.Add(udtSearchResult.Matched_peaks);              // Matched_Peak_Count
                lstData.Add(udtSearchResult.Matched_fragment_ions);      // Matched_Fragment_Ion_Count
                lstData.Add(udtSearchResult.Pvalue);
                lstData.Add(udtSearchResult.RankPValue.ToString());
                lstData.Add(udtSearchResult.Evalue);
                lstData.Add(udtSearchResult.FDR);
                lstData.Add(udtSearchResult.Species_ID);
                lstData.Add(udtSearchResult.FragMethod);

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

        #region "IComparer Classes"

        protected class MSAlignSearchResultsComparerScanChargePValuePeptide : IComparer<udtMSAlignSearchResultType>
        {
            public int Compare(udtMSAlignSearchResultType x, udtMSAlignSearchResultType y)
            {
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }
                else if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }
                else
                {
                    // Scan is the same, check charge
                    if (x.ChargeNum > y.ChargeNum)
                    {
                        return 1;
                    }
                    else if (x.ChargeNum < y.ChargeNum)
                    {
                        return -1;
                    }
                    else
                    {
                        // Charge is the same; check Pvalue
                        var result = x.Pvalue.CompareTo(y.Pvalue);
                        if (result == 0)
                        {
                            // Pvalue is the same; check peptide
                            result = x.Peptide.CompareTo(y.Peptide);
                            if (result == 0)
                            {
                                // Peptide is the same, check Protein
                                result = x.Protein.CompareTo(y.Protein);
                            }
                        }
                        return result;
                    }
                }
            }
        }

        protected class MSAlignSearchResultsComparerPValueScanChargePeptide : IComparer<udtMSAlignSearchResultType>
        {
            public int Compare(udtMSAlignSearchResultType x, udtMSAlignSearchResultType y)
            {
                var result1 = x.Pvalue.CompareTo(y.Pvalue);
                if (result1 == 0)
                {
                    // Pvalue is the same; check scan number
                    if (x.ScanNum > y.ScanNum)
                    {
                        return 1;
                    }
                    else if (x.ScanNum < y.ScanNum)
                    {
                        return -1;
                    }
                    else
                    {
                        // Scan is the same, check charge
                        if (x.ChargeNum > y.ChargeNum)
                        {
                            return 1;
                        }
                        else if (x.ChargeNum < y.ChargeNum)
                        {
                            return -1;
                        }
                        else
                        {
                            // Charge is the same; check peptide
                            var result = x.Peptide.CompareTo(y.Peptide);
                            if (result == 0)
                            {
                                // Peptide is the same, check Protein
                                result = x.Protein.CompareTo(y.Protein);
                            }
                            return result;
                        }
                    }
                }
                return result1;
            }
        }

        #endregion
    }
}
