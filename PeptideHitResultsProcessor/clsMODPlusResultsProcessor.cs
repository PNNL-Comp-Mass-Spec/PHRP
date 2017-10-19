// This class reads in an MODPlus results file (txt format) and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 05/15/2015
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------
using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsMODPlusResultsProcessor : clsPHRPBaseClass
    {
        public clsMODPlusResultsProcessor()
        {
            mFileDate = "October 13, 2017";
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_MODPlus_FILE = "_modp.id";

        public const string N_TERMINUS_SYMBOL_MODPlus = "-";
        public const string C_TERMINUS_SYMBOL_MODPlus = "-";

        // This is used for filtering both MODa and MODPlus results
        public const float DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD = 0.05f;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string MODPlus_MOD_MASS_REGEX = @"([+-][0-9.]+)";

        private const byte MODPlus_MASS_DIGITS_OF_PRECISION = 3;

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        // These columns correspond to the tab-delimited file (_MODPlus.id.txt) created by MODPlus's tda_plus.jar file
        protected const int MODPlusResultsFileColCount = 13;
        public enum eMODPlusResultsFileColumns
        {
            SpectrumFileName = 0,
            SpectrumIndex = 1,
            ScanNumber = 2,
            ObservedMonoMass = 3,
            Charge = 4,
            CalculatedMonoMass = 5,
            DeltaMass = 6,
            Score = 7,
            Probability = 8,
            Peptide = 9,
            NTT = 10,
            ProteinAndPeptidePositionList = 11,
            ModificationAnnotation = 12
        }

        // These columns correspond to the Synopsis file created by this class
        protected const int MODPlusSynFileColCount = 17;
        public enum eMODPlusSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Spectrum_Index = 2,
            Charge = 3,
            PrecursorMZ = 4,
            DelM = 5,                            // Precursor error, in Da
            DelM_PPM = 6,                        // Precursor error, in ppm
            MH = 7,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
            Peptide = 8,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
            NTT = 9,
            ModificationAnnotation = 10,
            Protein = 11,
            Peptide_Position = 12,
            Score = 13,
            Probability = 14,
            Rank_Score = 15,
            QValue = 16
        }

        #endregion

        #region "Structures"
        // This data structure holds rows read from the tab-delimited file (_MODPlus.id.txt) created by MODPlus's tda_plus.jar file
        protected struct udtMODPlusSearchResultType
        {
            public string SpectrumFileName;
            public string SpectrumIndex;
            public int ScanNum;
            public string Precursor_mass;           // Uncharged monoisotopic mass value of the observed precursor_mz, reported as ObservedMW by MODPlus
            public string PrecursorMZ;              // Computed by this class from ObservedMonoMass
            public string Charge;
            public short ChargeNum;
            public string CalculatedMonoMass;       // Theoretical monoisotopic mass of the peptide (including mods), as computed by MODPlus
            public string DeltaMass;                // Computed by MODPlus
            public string MH;                       // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                     // Computed by this class using Precursor_mass - CalculatedMonoMass
            public string DelM_PPM;                 // Computed by this class using DelM and CalculatedMonoMass
            public string Score;
            public double ScoreNum;
            public string Probability;              // Higher values are better
            public double ProbabilityNum;           // Higher values are better
            public int RankScore;
            public string Peptide;
            public string NTT;                      // Number of Tryptic Terminii
            public string ProteinList;              // One or more protein entries of the form ref|YP_003651515.1[K.196~206.Q(2)] where the text in brackets is the start/stop residues of the peptide; Multiple entries will be separated by semicolons, e.g. ref|YP_003651515.1[K.196~206.Q(2)];ref|YP_003201491.1[K.223~233.Q(2)];ref|YP_003313784.1[K.266~276.Q(2)]
            public string ModificationAnnotation;
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
                ScoreNum = 0;
                Probability = string.Empty;
                ProbabilityNum = 0;
                RankScore = 0;
                Peptide = string.Empty;
                NTT = string.Empty;
                ProteinList = string.Empty;
                ModificationAnnotation = string.Empty;
                FDR = 0;
                QValue = 0;
            }

            public override string ToString()
            {
                return Probability + ": " + Peptide;
            }
        }

        #endregion

        #region "Classwide Variables"
        protected Regex mProteinNamePositionSplit;
        #endregion

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to objSearchResult
        /// </summary>
        /// <param name="objSearchResult"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const char NO_RESIDUE = '-';

            var blnParsingModMass = false;

            var strModMassDigits = string.Empty;

            var chMostRecentResidue = NO_RESIDUE;
            var intResidueLocInPeptide = 0;

            var strSequence = objSearchResult.PeptideSequenceWithMods;
            for (var intIndex = 0; intIndex <= strSequence.Length - 1; intIndex++)
            {
                var chChar = strSequence[intIndex];

                if (IsLetterAtoZ(chChar))
                {
                    if (blnParsingModMass)
                    {
                        // Associate the mod mass in strModMassDigits with the previous residue
                        AssociateDynamicModWithResidue(objSearchResult, chMostRecentResidue, intResidueLocInPeptide, strModMassDigits, blnUpdateModOccurrenceCounts);
                        blnParsingModMass = false;
                    }

                    chMostRecentResidue = chChar;
                    intResidueLocInPeptide += 1;

                    // Look for static mods to associate with this residue
                    for (var intModIndex = 0; intModIndex <= mPeptideMods.ModificationCount - 1; intModIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(intModIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);

                            if (objModificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                objSearchResult.SearchResultAddModification(
                                    objModificationDefinition, chChar, intResidueLocInPeptide,
                                    objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else
                {
                    var blnIsNumberChar = chChar == '+' || chChar == '-' || char.IsDigit(chChar);

                    if (blnParsingModMass)
                    {
                        if (blnIsNumberChar || chChar == '.')
                        {
                            strModMassDigits += chChar;
                        }
                    }
                    else if (blnIsNumberChar)
                    {
                        // Mod Mass Start
                        strModMassDigits = chChar.ToString();
                        blnParsingModMass = true;
                    }
                    else
                    {
                        // Unrecognized symbol; ignore it
                    }
                }
            }

            if (blnParsingModMass)
            {
                // Associate the mod mass in strModMassDigits with the previous residue
                AssociateDynamicModWithResidue(objSearchResult, chMostRecentResidue, intResidueLocInPeptide, strModMassDigits, blnUpdateModOccurrenceCounts);
            }
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool blnSuccess;

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

        private void AssociateDynamicModWithResidue(
            clsSearchResultsBaseClass objSearchResult,
            char chMostRecentResidue,
            int intResidueLocInPeptide,
            string strModMassDigits,
            bool blnUpdateModOccurrenceCounts)
        {
            var chResidueForMod = chMostRecentResidue;
            var intResidueLocForMod = intResidueLocInPeptide;

            if (double.TryParse(strModMassDigits, out var dblModMass))
            {
                if (intResidueLocForMod == 0)
                {
                    // Modification is at the peptide N-terminus
                    intResidueLocForMod = 1;
                }

                var blnSuccess = objSearchResult.SearchResultAddModification(
                    dblModMass, chResidueForMod, intResidueLocForMod,
                    objSearchResult.DetermineResidueTerminusState(intResidueLocForMod),
                    blnUpdateModOccurrenceCounts, MODPlus_MASS_DIGITS_OF_PRECISION);

                if (!blnSuccess)
                {
                    var strErrorMessage = objSearchResult.ErrorMessage;
                    if (string.IsNullOrEmpty(strErrorMessage))
                    {
                        strErrorMessage = "SearchResultAddDynamicModification returned false for mod mass " + strModMassDigits;
                    }
                    SetErrorMessage(strErrorMessage + "; ResultID = " + objSearchResult.ResultID);
                }
            }
        }

        /// <summary>
        /// Ranks each entry assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="lstSearchResults"></param>
        /// <param name="intStartIndex"></param>
        /// <param name="intEndIndex"></param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            IList<udtMODPlusSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of lstSearchResults so that we can sort by descending Probability

            var dctResultsSubset = new Dictionary<int, udtMODPlusSearchResultType>();
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                dctResultsSubset.Add(intIndex, lstSearchResults[intIndex]);
            }

            var lstResultsByScore = (from item in dctResultsSubset orderby item.Value.ScoreNum descending select item).ToList();

            double dblLastValue = 0;
            var intCurrentRank = -1;

            foreach (var entry in lstResultsByScore)
            {
                var oResult = lstSearchResults[entry.Key];

                if (intCurrentRank < 0)
                {
                    dblLastValue = oResult.ScoreNum;
                    intCurrentRank = 1;
                }
                else
                {
                    if (Math.Abs(oResult.ScoreNum - dblLastValue) > double.Epsilon)
                    {
                        dblLastValue = oResult.ScoreNum;
                        intCurrentRank += 1;
                    }
                }

                oResult.RankScore = intCurrentRank;
                lstSearchResults[entry.Key] = oResult;
            }
        }

        protected string AssureInteger(string strInteger, int intDefaultValue)
        {
            if (strInteger.EndsWith(".0"))
                strInteger = strInteger.Substring(0, strInteger.Length - 2);

            if (int.TryParse(strInteger, out var intValue))
            {
                return intValue.ToString();
            }

            if (double.TryParse(strInteger, out var dblValue))
            {
                return dblValue.ToString("0");
            }

            return intDefaultValue.ToString();
        }

        protected double ComputePeptideMass(string strPeptide, double dblTotalModMass)
        {
            var strCleanSequence = GetCleanSequence(strPeptide);

            var dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence);

            if (Math.Abs(dblTotalModMass) > double.Epsilon)
            {
                dblMass += dblTotalModMass;
            }

            return dblMass;
        }

        private static readonly Regex RegexModMassRegEx = new Regex(MODPlus_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="strPeptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected double ComputeTotalModMass(string strPeptide)
        {
            double dblTotalModMass = 0;

            clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, out var strPrimarySequence, out var _, out var _);

            // Parse the dynamic mods reported by MODPlus
            foreach (Match reMatch in RegexModMassRegEx.Matches(strPrimarySequence))
            {
                // We use .TrimEnd() because the matched mod mass will end in a period if this mod applies to the final residue in a peptide
                if (double.TryParse(reMatch.Groups[1].Value.TrimEnd('.'), out var dblModMassFound))
                {
                    dblTotalModMass += dblModMassFound;
                }
            }

            // Now look for static mods
            // First determine the index of the last residue in strPrimarySequence
            var intIndexLastChar = strPrimarySequence.Length;

            for (var intIndex = strPrimarySequence.Length - 1; intIndex >= 0; intIndex += -1)
            {
                if (IsLetterAtoZ(strPrimarySequence[intIndex]))
                {
                    intIndexLastChar = intIndex;
                    break;
                }
            }

            for (var intIndex = 0; intIndex <= strPrimarySequence.Length - 1; intIndex++)
            {
                var chChar = strPrimarySequence[intIndex];

                if (IsLetterAtoZ(chChar))
                {
                    // Look for static mods to associate with this residue
                    for (var intModIndex = 0; intModIndex <= mPeptideMods.ModificationCount - 1; intModIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(intModIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);
                            var blnMatchFound = objModificationDefinition.TargetResiduesContain(chChar);

                            if (!blnMatchFound && intIndex == 0)
                            {
                                blnMatchFound = objModificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
                            }

                            if (!blnMatchFound && intIndex == intIndexLastChar)
                            {
                                blnMatchFound = objModificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);
                            }

                            if (blnMatchFound)
                            {
                                dblTotalModMass += objModificationDefinition.ModificationMass;
                            }
                        }
                    }
                }
            }

            return dblTotalModMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string strInputFilePath, string strOutputFolderPath, bool MTS)
        {
            var strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
            if (strPepToProteinMapFilePath.EndsWith("_MODPlus_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_MODPlus_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MODPlus
        /// The synopsis file includes every result with a probability above a set threshold
        /// The first-hits file includes the result with the highest probability (for each scan and charge)
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreateSynResultsFile(
            string strInputFilePath,
            string strOutputFilePath)
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

                        // Initialize the list that will hold all of the records in the MODPlus result file
                        var lstSearchResultsUnfiltered = new List<udtMODPlusSearchResultType>();

                        // Initialize the list that will hold all of the records that will ultimately be written out to disk
                        var lstFilteredSearchResults = new List<udtMODPlusSearchResultType>();

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !AbortProcessing)
                        {
                            var strLineIn = srDataFile.ReadLine();

                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            if (!headerParsed)
                            {
                                // Parse the header line

                                var blnsuccess = ParseMODPlusResultsFileHeaderLine(strLineIn, out intColumnMapping);
                                if (!blnsuccess)
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

                            var udtSearchResult = new udtMODPlusSearchResultType();

                            var blnValidSearchResult = ParseMODPlusResultsFileEntry(strLineIn, ref udtSearchResult, ref strErrorLog, intColumnMapping);

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

                        // Sort the SearchResults by scan, charge, and descending score
                        lstSearchResultsUnfiltered.Sort(new MODPlusSearchResultsComparerScanChargeScorePeptide());

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

        /// <summary>
        /// Load the static mods defined in the MODPlus parameter file
        /// </summary>
        /// <param name="strMODPlusParamFilePath"></param>
        /// <param name="lstModInfo"></param>
        /// <returns></returns>
        /// <remarks>We don't care about the dynamic mods because there are so many possible mods.  We'll add each dynamic mod as we encounter it in the results</remarks>
        protected bool ExtractModInfoFromMODPlusParamFile(string strMODPlusParamFilePath, ref List<clsModificationDefinition> lstModInfo)
        {
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

                if (string.IsNullOrEmpty(strMODPlusParamFilePath))
                {
                    SetErrorMessage("MODPlus Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                var paramFile = new FileInfo(strMODPlusParamFilePath);
                if (!paramFile.Exists)
                {
                    SetErrorMessage("MODPlus param file not found: " + strMODPlusParamFilePath);
                    return false;
                }

                // Read the contents of the parameter file
                var doc = new XmlDocument();
                doc.Load(paramFile.FullName);

                var nodeList = doc.SelectNodes("/search/modifications/fixed/mod");
                if (nodeList != null && nodeList.Count > 0)
                {
                    // Store the fixed mods

                    foreach (XmlNode node in nodeList)
                    {
                        var modName = node.Attributes["name"].Value;
                        var strResidue = node.Attributes["site"].Value.Trim();
                        var modPosition = node.Attributes["position"].Value;
                        var modMass = node.Attributes["massdiff"].Value;

                        // Replace N-Term or C-Term with < or >
                        if (strResidue.ToLower() == "n-term")
                            strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        if (strResidue.ToLower() == "c-term")
                            strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                        if (double.TryParse(modMass, out var modMassDa))
                        {
                            if (Math.Abs(modMassDa - 0) > float.Epsilon)
                            {
                                var strMassCorrectionTag = mPeptideMods.LookupMassCorrectionTagByMass(modMassDa);

                                var eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod;
                                if (strResidue == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || strResidue == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                                {
                                    eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                                }

                                var objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, strResidue, eModType, strMassCorrectionTag);
                                lstModInfo.Add(objModDef);
                            }
                        }
                    }
                }

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MODPlus parameter file (" + Path.GetFileName(strMODPlusParamFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFolderPath"></param>
        /// <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool ParseMODPlusSynopsisFile(
            string strInputFilePath,
            string strOutputFolderPath,
            bool blnResetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that MODPlus synopsis files are normally sorted on Probability value, ascending
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
                var objSearchResult = new clsSearchResultsMODPlus(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForProbabilityLevel
                var htPeptidesFoundForProbabilityLevel = new Hashtable();

                var strPreviousProbability = string.Empty;

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
                                blnSuccess = ParseMODPlusSynFileHeaderLine(strLineIn, out intColumnMapping);
                                if (!blnSuccess)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return blnSuccess;
                                }
                                blnHeaderParsed = true;
                                continue;
                            }

                            var strCurrentPeptideWithMods = string.Empty;

                            var blnValidSearchResult = ParseMODPlusSynFileEntry(strLineIn, objSearchResult, ref strErrorLog,
                                                                                intResultsProcessed, intColumnMapping,
                                                                                out strCurrentPeptideWithMods);

                            if (!blnValidSearchResult)
                            {
                                continue;
                            }

                            var strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.Charge;

                            bool blnFirstMatchForGroup;
                            if (objSearchResult.Probability == strPreviousProbability)
                            {
                                // New result has the same Probability as the previous result
                                // See if htPeptidesFoundForProbabilityLevel contains the peptide, scan and charge

                                if (htPeptidesFoundForProbabilityLevel.ContainsKey(strKey))
                                {
                                    blnFirstMatchForGroup = false;
                                }
                                else
                                {
                                    htPeptidesFoundForProbabilityLevel.Add(strKey, 1);
                                    blnFirstMatchForGroup = true;
                                }
                            }
                            else
                            {
                                // New Probability
                                // Reset htPeptidesFoundForProbabilityLevel
                                htPeptidesFoundForProbabilityLevel.Clear();

                                // Update strPreviousProbability
                                strPreviousProbability = objSearchResult.Probability;

                                // Append a new entry to htPeptidesFoundForProbabilityLevel
                                htPeptidesFoundForProbabilityLevel.Add(strKey, 1);
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
        /// Parse a MODPlus results line while creating the MODPlus synopsis file
        /// </summary>
        /// <param name="strLineIn"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="strErrorLog"></param>
        /// <param name="intColumnMapping"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool ParseMODPlusResultsFileEntry(
            string strLineIn,
            ref udtMODPlusSearchResultType udtSearchResult,
            ref string strErrorLog,
            IReadOnlyList<int> intColumnMapping)
        {
            // Parses an entry from the MODPlus results file

            var rowIndex = "?";

            bool blnValidSearchResult;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                udtSearchResult.Clear();
                var strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 11)
                {
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.SpectrumIndex], out udtSearchResult.SpectrumIndex))
                    {
                        ReportError("Index column is missing or invalid", true);
                    }
                    else
                    {
                        rowIndex = udtSearchResult.SpectrumIndex;
                    }

                    if (!int.TryParse(udtSearchResult.SpectrumIndex, out _))
                    {
                        ReportError("Index column is not numeric", true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.ScanNumber], out udtSearchResult.ScanNum);

                    // Monoisotopic mass value of the observed precursor_mz
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.ObservedMonoMass], out udtSearchResult.Precursor_mass);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    // dblPrecursorMonoMass is Observed m/z, converted to monoisotopic mass
                    if (double.TryParse(udtSearchResult.Precursor_mass, out var dblPrecursorMonoMass))
                    {
                        if (udtSearchResult.ChargeNum > 0)
                        {
                            var dblPrecursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMonoMass, 0, udtSearchResult.ChargeNum);
                            udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(dblPrecursorMZ, 6);
                        }
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                    // Theoretical peptide monoisotopic mass, including mods, as computed by MODPlus
                    double.TryParse(udtSearchResult.CalculatedMonoMass, out var dblPeptideMonoMassMODPlus);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.DeltaMass], out udtSearchResult.DeltaMass);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.Score], out udtSearchResult.Score);
                    if (!double.TryParse(udtSearchResult.Score, out udtSearchResult.ScoreNum))
                        udtSearchResult.ScoreNum = 0;

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.Probability], out udtSearchResult.Probability);
                    if (!double.TryParse(udtSearchResult.Probability, out udtSearchResult.ProbabilityNum))
                        udtSearchResult.ProbabilityNum = 0;

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                    {
                        ReportError("Peptide column is missing or invalid", true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.NTT], out udtSearchResult.NTT);

                    if (strSplitLine.Length > (int)eMODPlusResultsFileColumns.ProteinAndPeptidePositionList)
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.ProteinAndPeptidePositionList], out udtSearchResult.ProteinList);

                        // The protein column will have both the protein name and the peptide position
                        // For example, ref|YP_001038741.1[R.67~78.L(2)]
                        // It may have multiple proteins listed, separated by semicolons
                        // We will split the list on semicolons in function ParseMODPlusSynFileEntry

                        if (!udtSearchResult.ProteinList.Contains('['))
                        {
                            // This is likely a reverse-hit protein
                            udtSearchResult.ModificationAnnotation = string.Copy(udtSearchResult.ProteinList);
                            udtSearchResult.ProteinList = string.Empty;
                        }
                        else
                        {
                            if (strSplitLine.Length > (int)eMODPlusResultsFileColumns.ModificationAnnotation)
                            {
                                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusResultsFileColumns.ModificationAnnotation],
                                    out udtSearchResult.ModificationAnnotation);
                            }
                        }
                    }

                    // Parse the sequence to determine the total mod mass
                    // Note that we do not remove any of the mod symbols since MODPlus identifies mods by mass alone
                    // Note that static mods are implied (thus are not explicitly displayed by MODPlus)
                    var dblTotalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                    // Compute the theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                    var dblPeptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, dblTotalModMass);

                    // Only override dblPeptideMonoMassMODPlus if it is 0
                    if (Math.Abs(dblPeptideMonoMassMODPlus) < double.Epsilon)
                    {
                        dblPeptideMonoMassMODPlus = dblPeptideMonoMassPHRP;
                    }

                    var dblMassDiffThreshold = dblPeptideMonoMassMODPlus / 50000;
                    if (dblMassDiffThreshold < 0.1)
                        dblMassDiffThreshold = 0.1;

                    if (Math.Abs(dblPeptideMonoMassPHRP - dblPeptideMonoMassMODPlus) > dblMassDiffThreshold)
                    {
                        // Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da or by a slightly larger value if over 5000 Da; this is unexpected
                        string strFirst30Residues;
                        if (udtSearchResult.Peptide.Length < 27)
                        {
                            strFirst30Residues = udtSearchResult.Peptide;
                        }
                        else
                        {
                            strFirst30Residues = udtSearchResult.Peptide.Substring(0, 27) + "...";
                        }
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " + dblMassDiffThreshold.ToString("0.00") + " Da away from the mass computed by MODPlus: " + dblPeptideMonoMassPHRP.ToString("0.0000") + " vs. " + dblPeptideMonoMassMODPlus.ToString("0.0000") + "; peptide " + strFirst30Residues);
                    }

                    if (dblPeptideMonoMassMODPlus > 0)
                    {
                        // Compute DelM and DelM_PPM
                        var dblDelM = dblPrecursorMonoMass - dblPeptideMonoMassMODPlus;
                        udtSearchResult.DelM = MassErrorToString(dblDelM);

                        var dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, dblPrecursorMonoMass, true, dblPeptideMonoMassMODPlus);

                        udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(dblPeptideDeltaMassCorrectedPpm, 5, 0.00005);
                    }

                    // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPeptideMonoMassPHRP, 0), 6);

                    if (udtSearchResult.Probability.ToLower() == "infinity")
                    {
                        udtSearchResult.Probability = "0";
                    }
                    else if (!string.IsNullOrEmpty(udtSearchResult.Probability) & !double.TryParse(udtSearchResult.Probability, out _))
                    {
                        udtSearchResult.Probability = "";
                    }

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MODPlus results file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        strErrorLog += "Error parsing MODPlus Results in ParseMODPlusResultsFileEntry for RowIndex '" + rowIndex + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MODPlus Results in ParseMODPlusResultsFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="strLineIn"></param>
        /// <param name="intColumnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        /// <remarks></remarks>
        private bool ParseMODPlusResultsFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            // The expected column order from MODPlus:
            //   SpectrumFile   Index   ScanNo   ObservedMW   Charge   CalculatedMW   DeltaMass   Score   Probability   Peptide   NTT    Protein   ModificationAnnotation

            var lstColumnNames = new SortedDictionary<string, eMODPlusResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {"SpectrumFile", eMODPlusResultsFileColumns.SpectrumFileName},
                {"Index", eMODPlusResultsFileColumns.SpectrumIndex},
                {"ScanNo", eMODPlusResultsFileColumns.ScanNumber},
                {"ObservedMW", eMODPlusResultsFileColumns.ObservedMonoMass},
                {"Charge", eMODPlusResultsFileColumns.Charge},
                {"CalculatedMW", eMODPlusResultsFileColumns.CalculatedMonoMass},
                {"DeltaMass", eMODPlusResultsFileColumns.DeltaMass},
                {"Score", eMODPlusResultsFileColumns.Score},
                {"Probability", eMODPlusResultsFileColumns.Probability},
                {"Peptide", eMODPlusResultsFileColumns.Peptide},
                {"NTT", eMODPlusResultsFileColumns.NTT},
                {"Protein", eMODPlusResultsFileColumns.ProteinAndPeptidePositionList},
                {"ModificationAnnotation", eMODPlusResultsFileColumns.ModificationAnnotation}
            };

            intColumnMapping = new int[MODPlusResultsFileColCount];

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
                                Console.WriteLine("Warning: Unrecognized column header name '" + strSplitLine[intIndex] + "' in ParseMODPlusResultsFileHeaderLine");
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
                SetErrorMessage("Error parsing header in MODPlus results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        private bool ParseMODPlusSynFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            var lstColumnNames = new SortedDictionary<string, eMODPlusSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {clsPHRPParserMODPlus.DATA_COLUMN_ResultID, eMODPlusSynFileColumns.ResultID},
                {clsPHRPParserMODPlus.DATA_COLUMN_Scan, eMODPlusSynFileColumns.Scan},
                {clsPHRPParserMODPlus.DATA_COLUMN_Spectrum_Index, eMODPlusSynFileColumns.Spectrum_Index},
                {clsPHRPParserMODPlus.DATA_COLUMN_Charge, eMODPlusSynFileColumns.Charge},
                {clsPHRPParserMODPlus.DATA_COLUMN_PrecursorMZ, eMODPlusSynFileColumns.PrecursorMZ},
                {clsPHRPParserMODPlus.DATA_COLUMN_DelM, eMODPlusSynFileColumns.DelM},
                {clsPHRPParserMODPlus.DATA_COLUMN_DelM_PPM, eMODPlusSynFileColumns.DelM_PPM},
                {clsPHRPParserMODPlus.DATA_COLUMN_MH, eMODPlusSynFileColumns.MH},
                {clsPHRPParserMODPlus.DATA_COLUMN_Peptide, eMODPlusSynFileColumns.Peptide},
                {clsPHRPParserMODPlus.DATA_COLUMN_NTT, eMODPlusSynFileColumns.NTT},
                {clsPHRPParserMODPlus.DATA_COLUMN_Modification_Annotation, eMODPlusSynFileColumns.ModificationAnnotation},
                {clsPHRPParserMODPlus.DATA_COLUMN_Protein, eMODPlusSynFileColumns.Protein},
                {clsPHRPParserMODPlus.DATA_COLUMN_Peptide_Position, eMODPlusSynFileColumns.Peptide_Position},
                {clsPHRPParserMODPlus.DATA_COLUMN_Score, eMODPlusSynFileColumns.Score},
                {clsPHRPParserMODPlus.DATA_COLUMN_Probability, eMODPlusSynFileColumns.Probability},
                {clsPHRPParserMODPlus.DATA_COLUMN_Rank_Score, eMODPlusSynFileColumns.Rank_Score},
                {clsPHRPParserMODPlus.DATA_COLUMN_QValue, eMODPlusSynFileColumns.QValue}
            };

            intColumnMapping = new int[MODPlusSynFileColCount];

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
                SetErrorMessage("Error parsing header in MODPlus synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMODPlusSynFileEntry(
            string strLineIn,
            clsSearchResultsMODPlus objSearchResult,
            ref string strErrorLog,
            int intResultsProcessed,
            IReadOnlyList<int> intColumnMapping,
            out string strPeptideSequenceWithMods)
        {
            // Parses an entry from the MODPlus Synopsis file

            string[] strSplitLine = null;

            // Reset objSearchResult
            objSearchResult.Clear();
            strPeptideSequenceWithMods = string.Empty;

            try
            {

                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length < 13)
                {
                    return false;
                }

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.ResultID], out string strValue))
                {
                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        strErrorLog += "Error reading ResultID value from MODPlus Results line " + (intResultsProcessed + 1).ToString() +
                                       "\n";
                    }
                    return false;
                }

                objSearchResult.ResultID = int.Parse(strValue);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Scan], out string scan);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Charge], out string charge);

                objSearchResult.Scan = scan;
                objSearchResult.Charge = charge;

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Peptide], out strPeptideSequenceWithMods))
                {
                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        strErrorLog += "Error reading Peptide sequence value from MODPlus Results line " + (intResultsProcessed + 1).ToString() +
                                       "\n";
                    }
                    return false;
                }

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Protein], out string proteinName);
                objSearchResult.MultipleProteinCount = "0";

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.DelM], out string modPlusComputedDelM);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.DelM_PPM], out string modPlusComputedDelMppm);

                objSearchResult.ProteinName = proteinName;
                objSearchResult.MODPlusComputedDelM = modPlusComputedDelM;
                objSearchResult.MODPlusComputedDelMPPM = modPlusComputedDelMppm;

                objSearchResult.PeptideDeltaMass = objSearchResult.MODPlusComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MODPlus results file as "Observed_Mass - Theoretical_Mass"
                // However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                // Therefore, we will negate .peptideDeltaMass
                try
                {
                    objSearchResult.PeptideDeltaMass = (-double.Parse(objSearchResult.PeptideDeltaMass)).ToString(CultureInfo.InvariantCulture);
                }
                catch (Exception)
                {
                    // Error; Leave .peptideDeltaMass unchanged
                }

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                objSearchResult.SetPeptideSequenceWithMods(strPeptideSequenceWithMods, true, true);

                var objSearchResultBase = (clsSearchResultsBaseClass) objSearchResult;

                ComputePseudoPeptideLocInProtein(objSearchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                objSearchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Spectrum_Index], out string spectrumIndex);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.MH], out string parentIonMh);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Score], out string modPlusScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODPlusSynFileColumns.Probability], out string probability);

                objSearchResult.Spectrum_Index = spectrumIndex;
                objSearchResult.Precursor_mz = precursorMz;
                objSearchResult.ParentIonMH = parentIonMh;
                objSearchResult.MODPlusScore = modPlusScore;
                objSearchResult.Probability = probability;

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (strSplitLine != null && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing MODPlus Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MODPlus Results in ParseMODPlusSynFileEntry" + "\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">MODPlus results file (Dataset_MODPlus.id.txt)</param>
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

                    var lstMODPlusModInfo = new List<clsModificationDefinition>();

                    // Load the MODPlus Parameter File to look for any static mods
                    ExtractModInfoFromMODPlusParamFile(SearchToolParameterFilePath, ref lstMODPlusModInfo);

                    // Resolve the mods in lstMODPlusModInfo with the ModDefs mods
                    ResolveMODPlusModsWithModDefinitions(ref lstMODPlusModInfo);

                    // Define the base output filename using strInputFilePath
                    var strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath);

                    // Auto-replace "modp.id" with "_modp"
                    if (strBaseName.EndsWith("_modp.id", StringComparison.InvariantCultureIgnoreCase))
                    {
                        strBaseName = strBaseName.Substring(0, strBaseName.Length - "_modp.id".Length) + "_modp";
                    }

                    // Do not create a first-hits file for MODPlus results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_modp_syn.txt
                    var strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(strSynOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    blnSuccess = ParseMODPlusSynopsisFile(strSynOutputFilePath, strOutputFolderPath, false);

                    if (blnSuccess && CreateProteinModsFile)
                    {
                        blnSuccess = CreateProteinModsFileWork(strBaseName, inputFile, strSynOutputFilePath, strOutputFolderPath);
                    }

                    if (blnSuccess)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in clsMODPlusResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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
                    // Auto-change mIgnorePeptideToProteinMapperErrors to True
                    // We only do this because some peptides reported by MODPlus may not match the fasta file (due to amino acid substitutions)
                    IgnorePeptideToProteinMapperErrors = true;
                    blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                    if (!blnSuccess)
                    {
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath);

                // Create the Protein Mods file
                blnSuccess = CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MODPlus);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                blnSuccess = true;
            }

            return true;
        }

        protected void ResolveMODPlusModsWithModDefinitions(ref List<clsModificationDefinition> lstMODPlusModInfo)
        {
            var blnExistingModFound = false;
            clsModificationDefinition objModDef;

            if (lstMODPlusModInfo != null)
            {
                // Call .LookupModificationDefinitionByMass for each entry in lstMODPlusModInfo
                foreach (var objModInfo in lstMODPlusModInfo)
                {
                    if (string.IsNullOrEmpty(objModInfo.TargetResidues))
                    {
                        objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(
                            objModInfo.ModificationMass, objModInfo.ModificationType, default(char),
                            clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out blnExistingModFound, true);
                    }
                    else
                    {
                        foreach (var chTargetResidue in objModInfo.TargetResidues)
                        {
                            objModDef = mPeptideMods.LookupModificationDefinitionByMassAndModType(
                                objModInfo.ModificationMass, objModInfo.ModificationType, chTargetResidue,
                                clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out blnExistingModFound, true);
                        }
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter swResultFile,
            List<udtMODPlusSearchResultType> lstFilteredSearchResults,
            ref string strErrorLog)
        {
            // Sort udtFilteredSearchResults by descending score, ascending scan, ascending charge, ascending peptide, and ascending protein
            lstFilteredSearchResults.Sort(new MODPlusSearchResultsComparerScoreScanChargePeptide());

            // Compute FDR values then assign QValues
            ComputeQValues(lstFilteredSearchResults);

            if (mProteinNamePositionSplit == null)
            {
                mProteinNamePositionSplit = new Regex(@"(.+)\[([^\]]+)\]", RegexOptions.Compiled);
            }

            var resultID = 1;
            foreach (var result in lstFilteredSearchResults)
            {
                var proteinList = result.ProteinList.Split(';');

                if (proteinList.Length == 0)
                {
                    // This code should not be reached
                    WriteSearchResultToFile(resultID, swResultFile, result, "Unknown_Protein", string.Empty, ref strErrorLog);
                    resultID += 1;
                }

                foreach (var proteinEntry in proteinList)
                {
                    string proteinName;
                    string peptidePosition;

                    var reMatch = mProteinNamePositionSplit.Match(proteinEntry);
                    if (reMatch.Success)
                    {
                        proteinName = reMatch.Groups[1].Value;
                        peptidePosition = reMatch.Groups[2].Value;
                    }
                    else
                    {
                        proteinName = string.Copy(proteinEntry);
                        peptidePosition = string.Empty;
                    }

                    WriteSearchResultToFile(resultID, swResultFile, result, proteinName, peptidePosition, ref strErrorLog);
                    resultID += 1;
                }
            }
        }

        /// <summary>
        /// Compute FDR values then assign QValues
        /// </summary>
        /// <param name="lstSearchResults"></param>
        /// <remarks>Assumes the data is sorted by descending score using MODPlusSearchResultsComparerScoreScanChargePeptide</remarks>
        private void ComputeQValues(IList<udtMODPlusSearchResultType> lstSearchResults)
        {
            var forwardPeptideCount = 0;
            var reversePeptideCount = 0;

            for (var intIndex = 0; intIndex < lstSearchResults.Count;)
            {
                // Check for entries with multiple proteins listed
                var intIndexEnd = intIndex;
                while (intIndexEnd + 1 < lstSearchResults.Count)
                {
                    if (lstSearchResults[intIndex].ScanNum == lstSearchResults[intIndexEnd + 1].ScanNum &&
                        lstSearchResults[intIndex].ChargeNum == lstSearchResults[intIndexEnd + 1].ChargeNum &&
                        lstSearchResults[intIndex].Peptide == lstSearchResults[intIndexEnd + 1].Peptide)
                    {
                        intIndexEnd += 1;
                    }
                    else
                    {
                        break;
                    }
                }

                var isReverse = true;

                // Look for non-reverse proteins
                for (var intIndexCheck = intIndex; intIndexCheck <= intIndexEnd; intIndexCheck++)
                {
                    var proteinList = lstSearchResults[intIndexCheck].ProteinList.Split(';');

                    foreach (var proteinEntry in proteinList)
                    {
                        if (!IsReversedProtein(proteinEntry))
                        {
                            isReverse = false;
                            break;
                        }
                    }
                }

                if (isReverse)
                {
                    reversePeptideCount += 1;
                }
                else
                {
                    forwardPeptideCount += 1;
                }

                double dblFDR = 1;

                if (forwardPeptideCount > 0)
                {
                    dblFDR = reversePeptideCount / Convert.ToDouble(forwardPeptideCount);
                }

                // Store the FDR values
                for (var intIndexStore = intIndex; intIndexStore <= intIndexEnd; intIndexStore++)
                {
                    var udtResult = lstSearchResults[intIndexStore];
                    udtResult.FDR = dblFDR;

                    lstSearchResults[intIndexStore] = udtResult;
                }

                intIndex = intIndexEnd + 1;
            }

            // Now compute QValues
            // We step through the list, from the worst scoring result to the best result
            // The first QValue is the FDR of the final entry
            // The next QValue is the minimum of (QValue, CurrentFDR)

            var dblQValue = lstSearchResults.Last().FDR;
            if (dblQValue > 1)
                dblQValue = 1;

            for (var intIndex = lstSearchResults.Count - 1; intIndex >= 0; intIndex += -1)
            {
                var udtResult = lstSearchResults[intIndex];

                dblQValue = Math.Min(dblQValue, udtResult.FDR);
                udtResult.QValue = dblQValue;

                lstSearchResults[intIndex] = udtResult;
            }
        }

        private void StoreSynMatches(
            IList<udtMODPlusSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex,
            ICollection<udtMODPlusSearchResultType> lstFilteredSearchResults)
        {
            AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                if (lstSearchResults[intIndex].ProbabilityNum >= MODaMODPlusSynopsisFileProbabilityThreshold)
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
                    clsPHRPParserMODPlus.DATA_COLUMN_ResultID,
                    clsPHRPParserMODPlus.DATA_COLUMN_Scan,
                    clsPHRPParserMODPlus.DATA_COLUMN_Spectrum_Index,
                    clsPHRPParserMODPlus.DATA_COLUMN_Charge,
                    clsPHRPParserMODPlus.DATA_COLUMN_PrecursorMZ,
                    clsPHRPParserMODPlus.DATA_COLUMN_DelM,
                    clsPHRPParserMODPlus.DATA_COLUMN_DelM_PPM,
                    clsPHRPParserMODPlus.DATA_COLUMN_MH,
                    clsPHRPParserMODPlus.DATA_COLUMN_Peptide,
                    clsPHRPParserMODPlus.DATA_COLUMN_NTT,
                    clsPHRPParserMODPlus.DATA_COLUMN_Modification_Annotation,
                    clsPHRPParserMODPlus.DATA_COLUMN_Protein,
                    clsPHRPParserMODPlus.DATA_COLUMN_Peptide_Position,
                    clsPHRPParserMODPlus.DATA_COLUMN_Score,
                    clsPHRPParserMODPlus.DATA_COLUMN_Probability,
                    clsPHRPParserMODPlus.DATA_COLUMN_Rank_Score,
                    clsPHRPParserMODPlus.DATA_COLUMN_QValue
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
        /// Writes an entry to the synopsis file
        /// </summary>
        /// <param name="intResultID"></param>
        /// <param name="swResultFile"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="proteinName"></param>
        /// <param name="peptidePosition"></param>
        /// <param name="strErrorLog"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int intResultID,
            TextWriter swResultFile,
            udtMODPlusSearchResultType udtSearchResult,
            string proteinName,
            string peptidePosition,
            ref string strErrorLog)
        {
            try
            {
                // Primary Columns
                //
                // MODPlus
                // ResultID	Scan	Spectrum_Index	Charge	PrecursorMZ	DelM	DelM_PPM	MH	Peptide	NTT	ModificationAnnotation	Protein	Peptide_Position	Score	Probability	Rank_Probability   QValue

                var lstData = new List<string>
                {
                    intResultID.ToString(),
                    udtSearchResult.ScanNum.ToString(),
                    udtSearchResult.SpectrumIndex,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.DelM,
                    udtSearchResult.DelM_PPM,
                    udtSearchResult.MH,
                    udtSearchResult.Peptide,
                    udtSearchResult.NTT,
                    udtSearchResult.ModificationAnnotation,
                    proteinName,
                    peptidePosition,
                    udtSearchResult.Score,
                    udtSearchResult.Probability,
                    udtSearchResult.RankScore.ToString(),
                    PRISM.StringUtilities.DblToString(udtSearchResult.QValue, 5, 0.00005)
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

        #region "IComparer Classes"

        protected class MODPlusSearchResultsComparerScanChargeScorePeptide : IComparer<udtMODPlusSearchResultType>
        {
            public int Compare(udtMODPlusSearchResultType x, udtMODPlusSearchResultType y)
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

                // Charge is the same; check ScoreNum
                if (x.ScoreNum < y.ScoreNum)
                {
                    return 1;
                }

                if (x.ScoreNum > y.ScoreNum)
                {
                    return -1;
                }

                // Probability is the same; check peptide
                var result = string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.Compare(x.ProteinList, y.ProteinList, StringComparison.Ordinal);
                }
                return result;
            }
        }

        protected class MODPlusSearchResultsComparerScoreScanChargePeptide : IComparer<udtMODPlusSearchResultType>
        {
            public int Compare(udtMODPlusSearchResultType x, udtMODPlusSearchResultType y)
            {
                if (x.ScoreNum < y.ScoreNum)
                {
                    return 1;
                }

                if (x.ScoreNum > y.ScoreNum)
                {
                    return -1;
                }

                // Pvalue is the same; check scan number
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

                // Charge is the same; check peptide
                var result = string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.Compare(x.ProteinList, y.ProteinList, StringComparison.Ordinal);
                }
                return result;
            }
        }

        #endregion
    }
}
