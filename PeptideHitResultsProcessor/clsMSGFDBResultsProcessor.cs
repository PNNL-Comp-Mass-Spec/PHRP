// This class reads in an MSGF_DB results file (txt format) and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.
//
// The modification definition information is determined from the MSGF+ parameter file
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 8/12/2011
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
    public class clsMSGFDBResultsProcessor : clsPHRPBaseClass
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        public clsMSGFDBResultsProcessor()
        {
            mFileDate = "October 13, 2017";
            mModMassRegEx = new Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS);

            mPeptideCleavageStateCalculator = new clsPeptideCleavageStateCalculator();
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin);

            mNumericModErrors = 0;
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_MSGFDB_FILE = "_msgfdb";
        public const string FILENAME_SUFFIX_MSGFPLUS_FILE = "_msgfplus";

        public const string N_TERMINUS_SYMBOL_MSGFDB = "_.";
        public const string C_TERMINUS_SYMBOL_MSGFDB = "._";

        [Obsolete("Used by MSGF-DB; renamed to SpecEValue in MSGF+")]
        public const float DEFAULT_SYN_FILE_MSGF_SPECPROB_THRESHOLD = 5E-07f;

        [Obsolete("Used by MSGF-DB; renamed to EValue in MSGF+")]
        public const float DEFAULT_SYN_FILE_PVALUE_THRESHOLD = 0.75f;

        /// <summary>
        /// Filter passing peptides have MSGFDB_SpecEValue less than 5E-7 Or EValue less than 0.75 or QValue less than 10%
        /// This filter is also used by MSPathFinder
        /// </summary>
        public const float DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD = 5E-07f;

        /// <summary>
        /// Filter passing peptides have MSGFDB_SpecEValue less than 5E-7 Or EValue less than 0.75 or QValue less than 10%
        /// This filter is also used by MSPathFinder
        /// </summary>
        public const float DEFAULT_SYN_FILE_EVALUE_THRESHOLD = 0.75f;

        private const string SEARCH_ENGINE_NAME = "MSGF+";

        private const int MAX_ERROR_LOG_LENGTH = 4096;
        // Match mod masses (positive or negative) at start, e.g.
        // +57.021HWWTLTTDRINK         matches +57.021
        // -57.021+42.011HWWTLTTDRINK  matches -57.021+42.011 (two separate mods)
        // +42.011MDHTPQSQLK           matches +42.011

        private const string MSGFDB_NTERMINAL_MOD_MASS_REGEX = @"^([0-9\.\+\-]+)";
        // Match mod masses (positive or negative) at end, e.g.
        // FAACPLTCE+14.0157VS+79.9663+14.0157   matches +79.9663+14.0157

        private const string MSGFDB_CTERMINAL_MOD_MASS_REGEX = @"([0-9\.\+\-]+)$";
        private const string MSGFDB_MOD_MASS_REGEX = @"([+-][0-9\.]+)";

        private const string PROTEIN_AND_TERM_SYMBOLS_REGEX = @"([^;]+)\(pre=(.),post=(.)\)";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        // These columns correspond to the tab-delimited file created directly by MSGF+
        private const int MSGFDBResultsFileColCount = 20;
        public enum eMSGFDBResultsFileColumns
        {
            SpectrumFile = 0,
            SpecIndex = 1,               // SpecID in MSGF+
            Scan = 2,
            FragMethod = 3,
            PrecursorMZ = 4,
            PMErrorDa = 5,               // Corresponds to PMError(Da)
            PMErrorPPM = 6,              // Corresponds to PMError(ppm)
            Charge = 7,
            Peptide = 8,
            Protein = 9,
            DeNovoScore = 10,
            MSGFScore = 11,
            SpecProb_EValue = 12,
            PValue_EValue = 13,
            FDR_QValue = 14,             // Only present if searched using -tda 1
            PepFDR_PepQValue = 15,       // Only present if searched using -tda 1
            EFDR = 16,                   // Only present if did not search using -tda 1
            IMSScan = 17,                // Only present for MSGFDB_IMS results
            IMSDriftTime = 18,           // Only present for MSGFDB_IMS results
            IsotopeError = 19            // Only reported by MSGF+
        }

        // These columns correspond to the Synopsis and First-Hits files created by this class
        private const int MSGFDBSynFileColCount = 23;
        public enum eMSFDBSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            FragMethod = 2,
            SpecIndex = 3,
            Charge = 4,
            PrecursorMZ = 5,
            DelM = 6,                            // Precursor error, in Da; if the search used a tolerance less than 0.5 Da or less than 500 ppm, this value is computed from the DelMPPM value
            DelMPPM = 7,                         // Precursor error, in ppm; corrected for isotope selection errors
            MH = 8,                              // Theoretical monoisotopic peptide mass (computed by PHRP)
            Peptide = 9,                         // This is the sequence with prefix and suffix residues and also with modification symbols
            Protein = 10,                        // Protein Name (remove description)
            NTT = 11,                            // Number of tryptic terminii
            DeNovoScore = 12,
            MSGFScore = 13,
            SpecProb_EValue = 14,
            RankSpecProb = 15,                   // Rank 1 means lowest SpecEValue, 2 means next higher score, etc. (ties get the same rank)
            PValue_EValue = 16,
            FDR_QValue = 17,                     // Only present if searched using -tda 1
            PepFDR_PepQValue = 18,               // Only present if searched using -tda 1
            EFDR = 19,                           // Only present if did not search using -tda 1
            IMSScan = 20,                        // Only present for MSGFDB_IMS results
            IMSDriftTime = 21,                   // Only present for MSGFDB_IMS results
            IsotopeError = 22
        }

        private enum eFilteredOutputFileTypeConstants
        {
            SynFile = 0,
            FHTFile = 1
        }

        #endregion

        #region "Structures"
        private struct udtMSGFDBSearchResultType
        {
            public string SpectrumFileName;
            public string SpecIndex;
            public string Scan;
            public int ScanNum;
            public string FragMethod;
            public string PrecursorMZ;
            public string PMErrorDa;                // Corresponds to PMError(Da); MSGFDB stores this value as Observed - Theoretical
            public string PMErrorPPM;               // Corresponds to PMError(ppm); MSGFDB stores this value as Observed - Theoretical
            public string MH;
            public string Charge;
            public short ChargeNum;
            public string Peptide;                  // Peptide sequence, including prefix, suffix, and any mod symbols or mod masses
            public string Protein;
            public string NTT;
            public string DeNovoScore;
            public string MSGFScore;
            public string SpecEValue;               // Smaller values are better scores (e.g. 1E-9 is better than 1E-6); MSGF+ renamed this from SpecProb to SpecEValue
            public double SpecEValueNum;
            public string EValue;                   // Smaller values are better scores (e.g. 1E-7 is better than 1E-3); MSGF+ renamed this from PValue to EValue
            public double EValueNum;
            public string QValue;                   // Holds FDR when a target/decoy search was used; holds EFDR when a non-decoy search was used; holds QValue for MSGF+
            public double QValueNum;                // Numeric equivalent of QValue
            public string PepQValue;                // Only used when target/decoy search was used; holds PepQValue for MSGF+
            public int RankSpecProb;
            public int IMSScan;
            public string IMSDriftTime;
            public int IsotopeError;                // Only used by MSGF+

            public void Clear()
            {
                SpectrumFileName = string.Empty;
                SpecIndex = string.Empty;
                ScanNum = 0;
                FragMethod = string.Empty;
                PrecursorMZ = string.Empty;
                PMErrorDa = string.Empty;
                PMErrorPPM = string.Empty;
                MH = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                Peptide = string.Empty;
                Protein = string.Empty;
                NTT = string.Empty;
                DeNovoScore = string.Empty;
                MSGFScore = string.Empty;
                SpecEValue = string.Empty;
                SpecEValueNum = 0;
                EValue = string.Empty;
                EValueNum = 0;
                QValue = string.Empty;
                QValueNum = 0;
                PepQValue = string.Empty;
                RankSpecProb = 0;
                IMSScan = 0;
                IMSDriftTime = string.Empty;
                IsotopeError = 0;
            }
        }

        private struct udtScanGroupInfoType
        {
            public int ScanGroupID;
            public short Charge;
            public int Scan;
        }

        private struct udtTerminusCharsType
        {
            public char NTerm;
            public char CTerm;
        }

        private struct udtParentMassToleranceType
        {
            // Given a tolerance of 20ppm, we would have ToleranceLeft=20, ToleranceRight=20, and ToleranceIsPPM=True
            // Given a tolerance of 0.5Da,2.5Da, we would have ToleranceLeft=0.5, ToleranceRight=2.5, and ToleranceIsPPM=False
            public double ToleranceLeft;
            public double ToleranceRight;
            public bool IsPPM;

            public void Clear()
            {
                ToleranceLeft = 0;
                ToleranceRight = 0;
                IsPPM = false;
            }

            public override string ToString()
            {
                string units = null;
                double equivalenceThreshold = 0;

                if (IsPPM)
                {
                    units = "ppm";
                    equivalenceThreshold = 0.01;
                }
                else
                {
                    units = "Da";
                    equivalenceThreshold = 0.0001;
                }

                if (Math.Abs(ToleranceLeft - ToleranceRight) < equivalenceThreshold)
                {
                    return "+/-" + ToleranceLeft + " " + units;
                }
                else
                {
                    return "-" + ToleranceRight + ", +" + ToleranceLeft + " " + units;
                }
            }
        }

        #endregion

        #region "Classwide Variables"
        private readonly clsPeptideCleavageStateCalculator mPeptideCleavageStateCalculator;

        private udtParentMassToleranceType mParentMassToleranceInfo;

        private int mPrecursorMassErrorWarningCount;

        /// <summary>
        /// Looks for numeric mods in MSGF+ results
        /// For example, +14.016 in K.LQVPAGK+14.016ANPSPPIGPALGQR.G
        /// </summary>
        /// <remarks></remarks>
        private readonly Regex mModMassRegEx;

        private int mNumericModErrors;

        #endregion

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod symbol, determine the modification and add to objSearchResult
        /// </summary>
        /// <param name="objSearchResult"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            var intModIndex = 0;
            var chChar = default(char);
            var objModificationDefinition = default(clsModificationDefinition);

            string strSequence = null;
            var chMostRecentLetter = default(char);
            var intResidueLocInPeptide = 0;
            var blnSuccess = false;

            chMostRecentLetter = '-';
            intResidueLocInPeptide = 0;

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
                                objSearchResult.SearchResultAddModification(
                                    objModificationDefinition, chChar, intResidueLocInPeptide,
                                    objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(chMostRecentLetter))
                {
                    blnSuccess = objSearchResult.SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);
                    if (!blnSuccess)
                    {
                        var strErrorMessage = objSearchResult.ErrorMessage;
                        if (string.IsNullOrEmpty(strErrorMessage))
                        {
                            strErrorMessage = "SearchResultAddDynamicModification returned false for symbol " + chChar;
                        }
                        SetErrorMessage(strErrorMessage + "; ResultID = " + objSearchResult.ResultID);
                    }
                }
                else
                {
                    // We found a modification symbol but chMostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }
        }

        /// <summary>
        /// Adds or updates the prefix and suffix residues to the peptide, as defined in kvProteinInfo
        /// </summary>
        /// <param name="strPeptide"></param>
        /// <param name="kvProteinInfo"></param>
        /// <returns>Peptide sequence with N-terminal and C-Terminal residues</returns>
        /// <remarks></remarks>
        private string AddUpdatePrefixAndSuffixResidues(string strPeptide, KeyValuePair<string, udtTerminusCharsType> kvProteinInfo)
        {
            if (strPeptide.IndexOf('.') < 0)
            {
                return kvProteinInfo.Value.NTerm + "." + strPeptide + "." + kvProteinInfo.Value.CTerm;
            }

            string strPeptideNew;

            if (strPeptide.Length >= 2)
            {
                if (strPeptide[1] == '.')
                {
                    // Peptide already has the N-terminal residue
                    // Replace it using kvProteinInfo
                    strPeptideNew = kvProteinInfo.Value.NTerm + "." + strPeptide.Substring(2);
                }
                else if (strPeptide[0] == '.')
                {
                    strPeptideNew = kvProteinInfo.Value.NTerm + strPeptide;
                }
                else
                {
                    strPeptideNew = kvProteinInfo.Value.NTerm + "." + strPeptide;
                }
            }
            else
            {
                strPeptideNew = string.Copy(strPeptide);
            }

            if (strPeptideNew.Length >= 4)
            {
                if (strPeptideNew[strPeptideNew.Length - 2] == '.')
                {
                    // Peptide already has the C-terminal residue
                    // Replace it using kvProteinInfo
                    strPeptideNew = strPeptideNew.Substring(0, strPeptideNew.Length - 2) + "." + kvProteinInfo.Value.CTerm;
                }
                else if (strPeptideNew[strPeptideNew.Length - 1] == '.')
                {
                    strPeptideNew = strPeptideNew + kvProteinInfo.Value.CTerm;
                }
                else
                {
                    strPeptideNew = strPeptideNew + "." + kvProteinInfo.Value.CTerm;
                }
            }

            return strPeptideNew;
        }

        private void AppendToScanGroupDetails(
            ICollection<udtScanGroupInfoType> lstScanGroupDetails,
            IDictionary<string, bool> htScanGroupCombo,
            udtScanGroupInfoType udtScanGroupInfo,
            ref int intCurrentScanGroupID,
            ref int intNextScanGroupID)
        {
            var strChargeScanComboText = udtScanGroupInfo.Charge + "_" + udtScanGroupInfo.Scan;

            if (!htScanGroupCombo.ContainsKey(strChargeScanComboText))
            {
                if (intCurrentScanGroupID < 0)
                {
                    intCurrentScanGroupID = intNextScanGroupID;
                    intNextScanGroupID += 1;
                }

                udtScanGroupInfo.ScanGroupID = intCurrentScanGroupID;

                lstScanGroupDetails.Add(udtScanGroupInfo);
                htScanGroupCombo.Add(strChargeScanComboText, true);
            }
        }

        private void AppendToSearchResults(
            ICollection<udtMSGFDBSearchResultType> lstSearchResults,
            udtMSGFDBSearchResultType udtSearchResult,
            Dictionary<string, udtTerminusCharsType> lstProteinInfo)
        {
            if (lstProteinInfo.Count == 0)
            {
                lstSearchResults.Add(udtSearchResult);
            }
            else
            {
                foreach (var kvEntry in lstProteinInfo)
                {
                    udtSearchResult.Protein = kvEntry.Key;
                    udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, kvEntry);

                    lstSearchResults.Add(udtSearchResult);
                }
            }
        }

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="lstSearchResults">Search results</param>
        /// <param name="intStartIndex">Start index for data in this scan</param>
        /// <param name="intEndIndex">End index for data in this scan</param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            IList<udtMSGFDBSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex)
        {
            // Prior to September 2014 ranks were assigned per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            if (intStartIndex == intEndIndex)
            {
                // Only one result
                var currentResult = lstSearchResults[intStartIndex];
                currentResult.RankSpecProb = 1;
                lstSearchResults[intStartIndex] = currentResult;
                return;
            }

            // Duplicate a portion of udtSearchResults so that we can sort by ascending Spectral Probability

            var dctResultsSubset = new Dictionary<int, udtMSGFDBSearchResultType>();
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                dctResultsSubset.Add(intIndex, lstSearchResults[intIndex]);
            }

            var lstResultsBySpecProb = (from item in dctResultsSubset orderby item.Value.SpecEValueNum select item).ToList();

            double dblLastValue = 0;
            var intCurrentRank = -1;

            foreach (var entry in lstResultsBySpecProb)
            {
                var currentResult = lstSearchResults[entry.Key];

                if (intCurrentRank < 0)
                {
                    dblLastValue = currentResult.SpecEValueNum;
                    intCurrentRank = 1;
                }
                else
                {
                    if (Math.Abs(currentResult.SpecEValueNum - dblLastValue) > double.Epsilon)
                    {
                        dblLastValue = currentResult.SpecEValueNum;
                        intCurrentRank += 1;
                    }
                }

                currentResult.RankSpecProb = intCurrentRank;

                // Because this is a list of structs, we have to copy currentResult back into the current position in lstSearchResults
                lstSearchResults[entry.Key] = currentResult;
            }
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            var blnSuccess = false;

            try
            {
                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts);

                // Make sure .PeptideSequenceWithMods does not have any generic mod masses
                // It should only have mod symbols
                var reMatch = mModMassRegEx.Match(objSearchResult.PeptideSequenceWithMods);
                if (reMatch.Success)
                {
                    // Modification mass did not have a symbol associated with it in the _ModDefs.txt file
                    // We could try to handle this, listing the modification mass in place of the modification symbol in the _ModDetails.txt file, but will
                    // instead abort processing

                    mNumericModErrors += 1;

                    if (mNumericModErrors < 250)
                    {
                        var localErrorMessage = "Search result contains a numeric mod mass that could not be associated with a modification symbol; ResultID = " + objSearchResult.ResultID + ", ModMass = " + reMatch.Value;
                        SetErrorMessage(localErrorMessage);
                    }
                    else if (mNumericModErrors == 250)
                    {
                        SetErrorMessage("Too many numeric mod mass results have been found; suppressing further logging");
                    }

                    return false;
                }

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

        private short ComputeCleaveageState(string strSequenceWithMods)
        {
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;

            var eCleavageState = default(clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants);

            string strCleanSequence = null;

            // Remove any non-letter characters when before .ComputeCleavageState()
            strCleanSequence = GetCleanSequence(strSequenceWithMods, out strPrefix, out strSuffix);

            eCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(strCleanSequence, strPrefix, strSuffix);

            return Convert.ToInt16(eCleavageState);
        }

        /// <summary>
        /// This function should only be called when column PMError(Da) is present (and PMError(ppm) is not present)
        /// </summary>
        /// <param name="dblPrecursorErrorDa">Mass error (Observed - theoretical)</param>
        /// <param name="dblPrecursorMZ"></param>
        /// <param name="intCharge"></param>
        /// <param name="dblPeptideMonoisotopicMass"></param>
        /// <param name="blnAdjustPrecursorMassForC13"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private double ComputeDelMCorrectedPPM(
            double dblPrecursorErrorDa,
            double dblPrecursorMZ,
            int intCharge,
            double dblPeptideMonoisotopicMass,
            bool blnAdjustPrecursorMassForC13)
        {
            double dblPeptideDeltaMassCorrectedPpm = 0;

            double dblPrecursorMonoMass = 0;

            // Compute the original value for the precursor monoisotopic mass
            dblPrecursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, intCharge, 0);

            dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, dblPeptideMonoisotopicMass);

            return dblPeptideDeltaMassCorrectedPpm;
        }

        /// <summary>
        /// Compute the monoisotopic mass of the peptide
        /// </summary>
        /// <param name="strPeptide"></param>
        /// <param name="dblTotalModMass"></param>
        /// <returns></returns>
        private double ComputePeptideMass(string strPeptide, double dblTotalModMass)
        {
            string strCleanSequence = null;
            double dblMass = 0;

            strCleanSequence = GetCleanSequence(strPeptide);

            dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence);
            dblMass += dblTotalModMass;

            return dblMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string strInputFilePath, string strOutputFolderPath, bool MTS)
        {
            var strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);

            if (strPepToProteinMapFilePath.EndsWith("_msgfplus_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_msgfplus_fht", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_msgfdb_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_msgfdb_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS);
        }

        /// <summary>
        /// Parses the digits in strModDigits to convert them to one or more modification symbols
        /// </summary>
        /// <param name="currentResidue"></param>
        /// <param name="strModDigits">Example: +57.021 or +79.9663+14.0157 or -18.0106</param>
        /// <param name="strModSymbols"></param>
        /// <param name="strDynModSymbols"></param>
        /// <param name="lstMSGFDBModInfo"></param>
        /// <param name="blnNterminalMod"></param>
        /// <param name="blnPossibleCTerminalMod"></param>
        /// <param name="dblModMassFound"></param>
        /// <param name="blnContainsStaticMod"></param>
        /// <returns>True if success; false if a problem</returns>
        /// <remarks></remarks>
        private bool ConvertMGSFModMassesToSymbols(
            string currentResidue,
            string strModDigits,
            out string strModSymbols,
            out string strDynModSymbols,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstMSGFDBModInfo,
            bool blnNterminalMod,
            bool blnPossibleCTerminalMod,
            out double dblModMassFound,
            out bool blnContainsStaticMod)
        {
            var reMatches = default(MatchCollection);

            var blnMatchFound = false;
            var blnTestMod = false;

            var intBestMatchIndex = 0;
            double dblBestMassDiff = 0;
            var intModSymbolsFound = 0;
            var chSymbolBestMatch = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
            var residuesBestBatch = string.Empty;

            strModSymbols = string.Empty;
            strDynModSymbols = string.Empty;
            dblModMassFound = 0;
            blnContainsStaticMod = false;

            reMatches = mModMassRegEx.Matches(strModDigits);

            foreach (Match reMatch in reMatches)
            {
                string strModMass = null;
                strModMass = reMatch.Value;

                // Convert strModMass to a mass value
                double dblModMass = 0;
                dblModMass = double.Parse(strModMass);
                dblModMassFound += dblModMass;

                blnMatchFound = false;
                while (!blnMatchFound)
                {
                    intBestMatchIndex = -1;

                    // Step through the known modifications to find the closest match
                    for (var intIndex = 0; intIndex <= lstMSGFDBModInfo.Count - 1; intIndex++)
                    {
                        blnTestMod = true;

                        if (blnNterminalMod)
                        {
                            // Only test N-terminal mods
                            if (!(lstMSGFDBModInfo[intIndex].ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynNTermPeptide ||
                                  lstMSGFDBModInfo[intIndex].ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynNTermProtein))
                            {
                                blnTestMod = false;
                            }
                        }
                        else if (!blnPossibleCTerminalMod)
                        {
                            // Skip C-terminal mods since we're not at the C-terminus
                            if (lstMSGFDBModInfo[intIndex].ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynCTermPeptide ||
                                lstMSGFDBModInfo[intIndex].ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynCTermProtein)
                            {
                                blnTestMod = false;
                            }
                        }

                        if (blnTestMod)
                        {
                            var dblCandidateMassDiff = Math.Abs(lstMSGFDBModInfo[intIndex].ModMassVal - dblModMass);
                            if (dblCandidateMassDiff < 0.25)
                            {
                                // Possible match found
                                var updateCandidate = false;

                                if (intBestMatchIndex < 0 ||
                                    dblCandidateMassDiff < dblBestMassDiff)
                                {
                                    updateCandidate = true;
                                }
                                else if (Math.Abs(dblCandidateMassDiff - dblBestMassDiff) < float.Epsilon &&
                                         chSymbolBestMatch == '-' &&
                                         lstMSGFDBModInfo[intIndex].ModSymbol != '-')
                                {
                                    // Masses are the same, but the existing candidate is a static mod

                                    // If this new one is a dynamic mod, switch to it, but only if residuesBestBatch does not contain the residue while
                                    // the residues for the new candidate does contain the residue

                                    if (!residuesBestBatch.Contains(currentResidue) &&
                                        lstMSGFDBModInfo[intIndex].Residues.Contains(currentResidue))
                                    {
                                        updateCandidate = true;
                                    }
                                }

                                if (updateCandidate)
                                {
                                    intBestMatchIndex = intIndex;
                                    dblBestMassDiff = dblCandidateMassDiff;
                                    chSymbolBestMatch = lstMSGFDBModInfo[intIndex].ModSymbol;
                                    residuesBestBatch = string.Copy(lstMSGFDBModInfo[intIndex].Residues);
                                    blnMatchFound = true;
                                }
                            }
                        }
                    }

                    if (!blnMatchFound && blnNterminalMod)
                    {
                        // Set this to false, then search again
                        blnNterminalMod = false;
                    }
                    else if (!blnMatchFound && !blnPossibleCTerminalMod)
                    {
                        // Set this to true, then search again
                        blnPossibleCTerminalMod = true;
                    }
                    else
                    {
                        break;
                    }
                }

                if (blnMatchFound)
                {
                    // Match found; use the mod symbol
                    strModSymbols += lstMSGFDBModInfo[intBestMatchIndex].ModSymbol;
                    intModSymbolsFound += 1;

                    if (lstMSGFDBModInfo[intBestMatchIndex].ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.StaticMod)
                    {
                        blnContainsStaticMod = true;
                    }
                    else
                    {
                        strDynModSymbols += lstMSGFDBModInfo[intBestMatchIndex].ModSymbol;
                    }
                }
                else
                {
                    // Match not found; use the mass value
                    strModSymbols += strModMass;
                }
            }

            if (intModSymbolsFound > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MSGF+
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFilePath"></param>
        /// <param name="strScanGroupFilePath"></param>
        /// <param name="lstMSGFDBModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
        /// <param name="blnMSGFPlus">Output parameter: this function will set this to True if we're processing MSGF+ results</param>
        /// <param name="lstSpecIdToIndex"></param>
        /// <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool CreateFHTorSYNResultsFile(
            string strInputFilePath,
            string strOutputFilePath,
            string strScanGroupFilePath,
            List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstMSGFDBModInfo,
            out bool blnMSGFPlus,
            IDictionary<string, int> lstSpecIdToIndex,
            eFilteredOutputFileTypeConstants eFilteredOutputFileType)
        {
            string strLineIn = null;

            var lstSearchResultsCurrentScan = new List<udtMSGFDBSearchResultType>();
            var lstSearchResultsPrefiltered = new List<udtMSGFDBSearchResultType>();

            float sngPercentComplete = 0;

            var blnHeaderParsed = false;
            var blnIncludeFDRandPepFDR = false;
            var blnIncludeEFDR = false;
            var blnIncludeIMSFields = false;

            int[] intColumnMapping = null;

            var intNextScanGroupID = 0;
            var lstScanGroupDetails = default(List<udtScanGroupInfoType>);
            var htScanGroupCombo = default(Dictionary<string, bool>);

            var blnSuccess = false;
            var blnValidSearchResult = false;

            string strErrorLog = null;

            blnMSGFPlus = false;

            try
            {
                lstScanGroupDetails = new List<udtScanGroupInfoType>();
                htScanGroupCombo = new Dictionary<string, bool>();

                mPrecursorMassErrorWarningCount = 0;

                // Look for custom amino acids
                var lstCustomAA = (from item in lstMSGFDBModInfo where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA select item)
                    .ToList();
                foreach (var customAADef in lstCustomAA)
                {
                    var aminoAcidSymbol = customAADef.Residues[0];
                    var empiricalFormulaString = customAADef.ModMass;
                    var aminoAcidMass = customAADef.ModMassVal;

                    try
                    {
                        var elementalComposition = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormulaString);

                        mPeptideSeqMassCalculator.SetAminoAcidMass(aminoAcidSymbol, aminoAcidMass);
                        mPeptideSeqMassCalculator.SetAminoAcidAtomCounts(aminoAcidSymbol, elementalComposition);
                    }
                    catch (Exception ex)
                    {
                        ReportError(ex.Message);
                    }
                }

                try
                {
                    // Open the input file and parse it
                    // Initialize the stream reader and the stream Text writer
                    using (var srDataFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        using (var swResultFile = new StreamWriter(new FileStream(strOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                        {
                            strErrorLog = string.Empty;
                            blnHeaderParsed = false;

                            intNextScanGroupID = 1;
                            lstScanGroupDetails.Clear();
                            htScanGroupCombo.Clear();

                            // Initialize the array that will hold all of the records that will ultimately be written out to disk
                            var lstFilteredSearchResults = new List<udtMSGFDBSearchResultType>();

                            // Initialize a dictionary that tracks the peptide sequence for each combo of scan and charge
                            // Keys are Scan_Charge, values track the clean sequence, the associated protein name, and the protein number for that name
                            // Note that we can only track protein numbers if the FASTA file path was provided at the command line
                            var scanChargeFirstHit = new Dictionary<string, clsFirstHitInfo>();

                            // Parse the input file
                            while (!srDataFile.EndOfStream & !AbortProcessing)
                            {
                                strLineIn = srDataFile.ReadLine();
                                if (string.IsNullOrWhiteSpace(strLineIn))
                                {
                                    continue;
                                }

                                if (!blnHeaderParsed)
                                {
                                    blnSuccess = ParseMSGFDBResultsFileHeaderLine(strLineIn, out intColumnMapping);
                                    if (!blnSuccess)
                                    {
                                        // Error parsing header
                                        SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                        return false;
                                    }
                                    blnHeaderParsed = true;

                                    blnIncludeFDRandPepFDR = false;
                                    blnIncludeEFDR = false;

                                    if (intColumnMapping[(int)eMSGFDBResultsFileColumns.FDR_QValue] >= 0 || intColumnMapping[(int)eMSGFDBResultsFileColumns.PepFDR_PepQValue] >= 0)
                                    {
                                        blnIncludeFDRandPepFDR = true;
                                    }
                                    else if (intColumnMapping[(int)eMSGFDBResultsFileColumns.EFDR] >= 0)
                                    {
                                        blnIncludeEFDR = true;
                                    }

                                    if (intColumnMapping[(int)eMSGFDBResultsFileColumns.IMSDriftTime] >= 0)
                                    {
                                        blnIncludeIMSFields = true;
                                    }

                                    if (intColumnMapping[(int)eMSGFDBResultsFileColumns.IsotopeError] >= 0)
                                    {
                                        blnMSGFPlus = true;
                                    }

                                    // Write the header line
                                    WriteSynFHTFileHeader(swResultFile, ref strErrorLog, blnIncludeFDRandPepFDR, blnIncludeEFDR, blnIncludeIMSFields,
                                        blnMSGFPlus);
                                    continue;
                                }

                                blnValidSearchResult = ParseMSGFDBResultsFileEntry(strLineIn, blnMSGFPlus, lstMSGFDBModInfo,
                                                                                   lstSearchResultsCurrentScan, ref strErrorLog,
                                                                                   intColumnMapping, ref intNextScanGroupID, lstScanGroupDetails, htScanGroupCombo, lstSpecIdToIndex);

                                if (blnValidSearchResult && lstSearchResultsCurrentScan.Count > 0)
                                {
                                    if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                                    {
                                        // Synopsis file
                                        blnValidSearchResult = MSGFPlusResultPassesSynFilter(lstSearchResultsCurrentScan[0]);
                                    }
                                    else
                                    {
                                        // First Hits file
                                        var scanChargeKey = lstSearchResultsCurrentScan[0].Scan + "_" + lstSearchResultsCurrentScan[0].Charge;

                                        if (scanChargeFirstHit.TryGetValue(scanChargeKey, out var firstHitPeptide))
                                        {
                                            // A result has already been stored for this scan/charge combo
                                            blnValidSearchResult = false;

                                            // Possibly update the associated protein name
                                            if (firstHitPeptide.CleanSequence.Equals(GetCleanSequence(lstSearchResultsCurrentScan[0].Peptide, out var strNewPrefix, out var strNewSuffix)))
                                            {
                                                var bestProtein = GetBestProteinName(firstHitPeptide.ProteinName, firstHitPeptide.ProteinNumber, lstSearchResultsCurrentScan[0].Protein);
                                                if (bestProtein.Value < firstHitPeptide.ProteinNumber)
                                                {
                                                    firstHitPeptide.ProteinName = bestProtein.Key;
                                                    firstHitPeptide.ProteinNumber = bestProtein.Value;
                                                    firstHitPeptide.UpdatePrefixAndSuffix(strNewPrefix, strNewSuffix);
                                                }
                                            }
                                        }
                                        else
                                        {
                                            firstHitPeptide = new clsFirstHitInfo(lstSearchResultsCurrentScan[0].Peptide, GetCleanSequence(lstSearchResultsCurrentScan[0].Peptide))
                                            {
                                                ProteinName = lstSearchResultsCurrentScan[0].Protein,
                                                ProteinNumber = int.MaxValue
                                            };

                                            if (mProteinNameOrder.TryGetValue(lstSearchResultsCurrentScan[0].Protein, out var proteinNumber))
                                            {
                                                firstHitPeptide.ProteinNumber = proteinNumber;
                                            }

                                            scanChargeFirstHit.Add(scanChargeKey, firstHitPeptide);
                                        }
                                    }

                                    if (blnValidSearchResult)
                                    {
                                        ExpandListIfRequired(lstSearchResultsPrefiltered, lstSearchResultsCurrentScan.Count);
                                        lstSearchResultsPrefiltered.AddRange(lstSearchResultsCurrentScan);
                                    }
                                }

                                // Update the progress
                                sngPercentComplete = Convert.ToSingle(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100);
                                if (CreateProteinModsFile)
                                {
                                    sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                                }
                                UpdateProgress(sngPercentComplete);
                            }

                            lstSearchResultsPrefiltered.TrimExcess();

                            // Sort the SearchResults by scan, charge, and ascending SpecEValue
                            lstSearchResultsPrefiltered.Sort(new MSGFDBSearchResultsComparerScanChargeSpecEValuePeptide());

                            if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.FHTFile)
                            {
                                // Update the protein names in lstSearchResultsPrefiltered using scanChargeFirstHit
                                // This step is likely unnecessary, thus the "Unexpected code reached" message below

                                for (var intIndex = 0; intIndex <= lstSearchResultsPrefiltered.Count - 1; intIndex++)
                                {
                                    var scanChargeKey = lstSearchResultsPrefiltered[intIndex].Scan + "_" + lstSearchResultsPrefiltered[intIndex].Charge;
                                    clsFirstHitInfo firstHitPeptide = null;

                                    if (scanChargeFirstHit.TryGetValue(scanChargeKey, out firstHitPeptide))
                                    {
                                        if (!lstSearchResultsPrefiltered[intIndex].Protein.Equals(firstHitPeptide.ProteinName))
                                        {
                                            Console.WriteLine("Unexpected code reached; possible logic error in clsMSGFDBResultsProcessor.CreateFHTorSYNCResultsFile");

                                            if (firstHitPeptide.CleanSequence.Equals(GetCleanSequence(lstSearchResultsPrefiltered[intIndex].Peptide)))
                                            {
                                                var updatedSearchResult = lstSearchResultsPrefiltered[intIndex];
                                                updatedSearchResult.Peptide = firstHitPeptide.SequenceWithModsAndContext;
                                                updatedSearchResult.Protein = string.Copy(firstHitPeptide.ProteinName);
                                                lstSearchResultsPrefiltered[intIndex] = updatedSearchResult;
                                            }
                                            else
                                            {
                                                Console.WriteLine(string.Format("Possible programming bug; " +
                                                                                "mix of peptides tracked for a given scan/charge combo when caching data for First Hits files; " +
                                                                                "see scan_charge {0}", scanChargeKey));
                                            }
                                        }
                                    }
                                }
                            }

                            // Now filter the data and store in lstFilteredSearchResults
                            // Due to code updates in October 2016, lstSearchResultsPrefiltered already has filtered data
                            var intStartIndex = 0;
                            var intEndIndex = 0;

                            while (intStartIndex < lstSearchResultsPrefiltered.Count)
                            {
                                intEndIndex = intStartIndex;
                                // Find all of the peptides with the same scan number
                                while (intEndIndex + 1 < lstSearchResultsPrefiltered.Count &&
                                       lstSearchResultsPrefiltered[intEndIndex + 1].ScanNum == lstSearchResultsPrefiltered[intStartIndex].ScanNum)
                                {
                                    intEndIndex += 1;
                                }

                                // Store the results for this scan
                                if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                                {
                                    StoreSynMatches(lstSearchResultsPrefiltered, intStartIndex, intEndIndex, lstFilteredSearchResults);
                                }
                                else
                                {
                                    StoreTopFHTMatch(lstSearchResultsPrefiltered, intStartIndex, intEndIndex, lstFilteredSearchResults);
                                }

                                intStartIndex = intEndIndex + 1;
                            }

                            // Sort the data in udtFilteredSearchResults then write out to disk
                            SortAndWriteFilteredSearchResults(swResultFile, lstFilteredSearchResults, ref strErrorLog,
                                                              blnIncludeFDRandPepFDR, blnIncludeEFDR, blnIncludeIMSFields, blnMSGFPlus);
                        }
                    }

                    // Write out the scan group info
                    if (!string.IsNullOrEmpty(strScanGroupFilePath))
                    {
                        StoreScanGroupInfo(strScanGroupFilePath, lstScanGroupDetails);
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
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        /// <summary>
        /// Extracts mod info from either a MSGF+ param file or from a MSGFPlus_Mods.txt file (previously MSGFDB_Mods.txt)
        /// </summary>
        /// <param name="strMSGFDBParamFilePath"></param>
        /// <param name="lstModInfo"></param>
        /// <returns>True if success; false if a problem</returns>
        /// <remarks></remarks>
        private bool ExtractModInfoFromParamFile(
            string strMSGFDBParamFilePath,
            out List<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstModInfo)
        {
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            // Note that this call will initialize lstModInfo
            var success = modFileProcessor.ExtractModInfoFromParamFile(strMSGFDBParamFilePath, out lstModInfo);

            if (!success || mErrorCode != ePHRPErrorCodes.NoError)
            {
                if (mErrorCode == ePHRPErrorCodes.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MSGF+ parameter file");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFDBModsWithModDefinitions(lstModInfo, mPeptideMods);

            return true;
        }

        /// <summary>
        /// Extracts parent mass tolerance from the parameters loaded from an MSGF+ parameter file
        /// </summary>
        /// <param name="objSearchEngineParams"></param>
        /// <returns>Parent mass tolerance info.  Tolerances will be 0 if an error occurs</returns>
        /// <remarks></remarks>
        private udtParentMassToleranceType ExtractParentMassToleranceFromParamFile(clsSearchEngineParameters objSearchEngineParams)
        {
            const string PM_TOLERANCE_TAG = "PMTolerance";

            var udtParentMassToleranceInfo = new udtParentMassToleranceType();

            try
            {
                udtParentMassToleranceInfo.Clear();

                if (objSearchEngineParams.Parameters.TryGetValue(PM_TOLERANCE_TAG, out var strValue))
                {
                    // Parent ion tolerance line found

                    // Split the line on commas
                    var strSplitLine = strValue.Split(',');

                    double dblTolerance;
                    bool blnIsPPM;

                    if (strSplitLine.Length == 1)
                    {
                        if (ParseParentMassTolerance(strSplitLine[0], out dblTolerance, out blnIsPPM))
                        {
                            udtParentMassToleranceInfo.ToleranceLeft = dblTolerance;
                            udtParentMassToleranceInfo.ToleranceRight = dblTolerance;
                            udtParentMassToleranceInfo.IsPPM = blnIsPPM;
                        }
                    }
                    else if (strSplitLine.Length > 1)
                    {
                        if (ParseParentMassTolerance(strSplitLine[0], out dblTolerance, out blnIsPPM))
                        {
                            udtParentMassToleranceInfo.ToleranceLeft = dblTolerance;
                            udtParentMassToleranceInfo.IsPPM = blnIsPPM;

                            if (ParseParentMassTolerance(strSplitLine[1], out dblTolerance, out blnIsPPM))
                            {
                                udtParentMassToleranceInfo.ToleranceRight = dblTolerance;
                            }
                            else
                            {
                                udtParentMassToleranceInfo.ToleranceRight = udtParentMassToleranceInfo.ToleranceLeft;
                            }
                        }
                    }
                }

                Console.WriteLine();
            }
            catch (Exception ex)
            {
                SetErrorMessage(string.Format("Error parsing the ParentMass tolerance from the MSGF+ parameter file ({0}): {1}",
                    Path.GetFileName(objSearchEngineParams.SearchEngineParamFilePath), ex.Message), ex);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
            }

            return udtParentMassToleranceInfo;
        }

        /// <summary>
        /// Look for candidateProteinName in mProteinNameOrder
        /// If found, and if the proteinNumber for that protein is less than currentProteinNumber,
        /// return a KeyValuePair with candidateProteinName and the proteinNumber for that protein
        /// Otherwise, return a KeyValuePair with currentProteinName and currentProteinNumber
        /// </summary>
        /// <param name="currentProteinName"></param>
        /// <param name="currentProteinNumber"></param>
        /// <param name="candidateProteinName"></param>
        /// <returns></returns>
        private KeyValuePair<string, int> GetBestProteinName(string currentProteinName, int currentProteinNumber, string candidateProteinName)
        {
            if (mProteinNameOrder.Count > 0)
            {
                // Lookup the protein number (to make sure we use the protein name that occurs first in the FASTA file)
                // Only possible if the user provided the path to the FASTA file

                if (mProteinNameOrder.TryGetValue(candidateProteinName, out var proteinNumber))
                {
                    if (proteinNumber < currentProteinNumber)
                    {
                        // A better protein name has been found (or this is the first protein and we just determined the protein number to associate with it)
                        return new KeyValuePair<string, int>(candidateProteinName, proteinNumber);
                    }
                }
                else
                {
                    // Protein not found in mProteinNameOrder
                    // It's likely a reverse-hit protein
                }
            }

            // A better protein name was not found; return the current info
            return new KeyValuePair<string, int>(currentProteinName, currentProteinNumber);
        }

        /// <summary>
        /// Load the PeptideToProteinMap information; in addition, creates the _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
        /// </summary>
        /// <param name="strPepToProteinMapFilePath"></param>
        /// <param name="strOutputFolderPath"></param>
        /// <param name="lstMSGFDBModInfo"></param>
        /// <param name="blnMSGFPlus">Should be set to True if processing MSGF+ results</param>
        /// <param name="lstPepToProteinMapping"></param>
        /// <param name="strMTSPepToProteinMapFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadPeptideToProteinMapInfoMSGFDB(
            string strPepToProteinMapFilePath,
            string strOutputFolderPath,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstMSGFDBModInfo,
            bool blnMSGFPlus,
            List<udtPepToProteinMappingType> lstPepToProteinMapping,
            out string strMTSPepToProteinMapFilePath)
        {
            var strHeaderLine = string.Empty;
            string strMTSCompatiblePeptide = null;

            // Not used by this function but required for the call to ReplaceMSGFModTextWithSymbol
            double dblTotalModMass = 0;

            var blnSuccess = false;

            strMTSPepToProteinMapFilePath = string.Empty;

            try
            {
                if (string.IsNullOrWhiteSpace(strPepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file is not defined");
                    return false;
                }

                if (!File.Exists(strPepToProteinMapFilePath))
                {
                    var strPepToProteinMapAlternate = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strPepToProteinMapFilePath, "Dataset_msgfdb.txt");
                    if (File.Exists(strPepToProteinMapAlternate))
                    {
                        strPepToProteinMapFilePath = strPepToProteinMapAlternate;
                    }
                }

                if (!File.Exists(strPepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file does not exist: " + strPepToProteinMapFilePath);
                    return false;
                }

                // Initialize lstPepToProteinMapping
                lstPepToProteinMapping.Clear();

                // Read the data in strProteinToPeptideMappingFilePath
                blnSuccess = LoadPeptideToProteinMapInfo(strPepToProteinMapFilePath, lstPepToProteinMapping, out strHeaderLine);

                if (blnSuccess)
                {
                    strMTSPepToProteinMapFilePath = Path.Combine(strOutputFolderPath, Path.GetFileNameWithoutExtension(strPepToProteinMapFilePath) + "MTS.txt");

                    using (var swOutFile = new StreamWriter(new FileStream(strMTSPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
                    )
                    {
                        if (!string.IsNullOrEmpty(strHeaderLine))
                        {
                            // Header line
                            swOutFile.WriteLine(strHeaderLine);
                        }

                        for (var intIndex = 0; intIndex <= lstPepToProteinMapping.Count - 1; intIndex++)
                        {
                            // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                            // In addition, replace the * terminus symbols with dashes
                            strMTSCompatiblePeptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(lstPepToProteinMapping[intIndex].Peptide), lstMSGFDBModInfo, blnMSGFPlus, out dblTotalModMass);

                            if (lstPepToProteinMapping[intIndex].Peptide != strMTSCompatiblePeptide)
                            {
                                UpdatePepToProteinMapPeptide(lstPepToProteinMapping, intIndex, strMTSCompatiblePeptide);
                            }

                            swOutFile.WriteLine(lstPepToProteinMapping[intIndex].Peptide + "\t" +
                                                lstPepToProteinMapping[intIndex].Protein + "\t" +
                                                lstPepToProteinMapping[intIndex].ResidueStart + "\t" +
                                                lstPepToProteinMapping[intIndex].ResidueEnd);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" + Path.GetFileName(strMTSPepToProteinMapFilePath) + "): " + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        /// <summary>
        /// Load the MSGF+ parameter file and updates settings
        /// </summary>
        /// <param name="msgfPlusParamFilePath"></param>
        /// <returns>
        /// True if success, false if an error.
        /// Returns True if msgfPlusParamFilePath is empty
        /// Returns False if the paramFilePath is defined but the file is not found or cannot be parsed</returns>
        private bool LoadSearchEngineParamFile(string msgfPlusParamFilePath)
        {
            if (string.IsNullOrWhiteSpace(msgfPlusParamFilePath))
            {
                ReportWarning("MSGF+ parameter file is not defined. Unable to extract parent mass tolerance info or custom charge carrier masses");
                return true;
            }

            var objSearchEngineParams = new clsSearchEngineParameters(SEARCH_ENGINE_NAME);

            var success = clsPHRPParser.ReadKeyValuePairSearchEngineParamFile(SEARCH_ENGINE_NAME, msgfPlusParamFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB,
                                                                              objSearchEngineParams, out var localErrorMessage, out var localWarningMessage);

            if (!string.IsNullOrWhiteSpace(localErrorMessage))
            {
                ReportError(localErrorMessage);
                return false;
            }

            if (!string.IsNullOrWhiteSpace(localWarningMessage))
            {
                ReportWarning(localWarningMessage);
            }

            if (objSearchEngineParams == null || objSearchEngineParams.Parameters.Count == 0)
            {
                SetErrorMessage("MSGF+ parameter file is empty; unable to extract parent mass tolerance info");
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }

            // Parse the PMTolerance setting
            mParentMassToleranceInfo = ExtractParentMassToleranceFromParamFile(objSearchEngineParams);

            // Parse the ChargeCarrierMass setting
            if (clsPHRPParserMSGFDB.GetCustomChargeCarrierMass(objSearchEngineParams, out var customChargeCarrierMass))
            {
                ReportMessage(string.Format("Using a charge carrier mass of {0:F3} Da", customChargeCarrierMass));
                mPeptideSeqMassCalculator.ChargeCarrierMass = customChargeCarrierMass;
            }

            return success;
        }

        private void ModExtractorErrorHandler(string message, Exception ex)
        {
            SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
        }

        private bool MSGFPlusResultPassesSynFilter(udtMSGFDBSearchResultType udtMsgfdbSearchResultType)
        {
            if (udtMsgfdbSearchResultType.EValueNum <= MSGFDBSynopsisFileEValueThreshold ||
                udtMsgfdbSearchResultType.SpecEValueNum <= MSGFDBSynopsisFileSpecEValueThreshold ||
                udtMsgfdbSearchResultType.QValueNum > 0 && udtMsgfdbSearchResultType.QValueNum < 0.01)
            {
                return true;
            }

            return false;
        }

        private bool ParseMSGFDBSynopsisFile(
            string strInputFilePath,
            string strOutputFolderPath,
            List<udtPepToProteinMappingType> lstPepToProteinMapping,
            bool blnResetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            string strPreviousSpecEValue = null;
            string strKey = null;

            string strLineIn = null;
            string strModificationSummaryFilePath = null;

            var objSearchResult = default(clsSearchResultsMSGFDB);

            var intResultsProcessed = 0;
            var intPepToProteinMapIndex = 0;
            float sngPercentComplete = 0;

            var blnHeaderParsed = false;
            int[] intColumnMapping = null;

            var strCurrentPeptideWithMods = string.Empty;
            string strCurrentProtein = null;

            var blnSuccess = false;
            var blnValidSearchResult = false;
            var blnFirstMatchForGroup = false;

            var successOverall = true;

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (blnResetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                mNumericModErrors = 0;

                // Initialize objSearchResult
                objSearchResult = new clsSearchResultsMSGFDB(mPeptideMods, mPeptideSeqMassCalculator);

                // Note that MSGF+ synopsis files are normally sorted on SpecEValue value, ascending
                // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
                //  we will keep track of the scan, charge, and peptide information parsed for each unique SpecEValue encountered
                // Although this was a possiblity with Inspect, it likely never occurs for MSGF+
                //  But, we'll keep the check in place just in case

                var peptidesFoundForSpecEValueLevel = new SortedSet<string>();

                strPreviousSpecEValue = string.Empty;

                // Assure that lstPepToProteinMapping is sorted on peptide
                if (lstPepToProteinMapping.Count > 1)
                {
                    lstPepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                try
                {
                    objSearchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    var strErrorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var srDataFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        intResultsProcessed = 0;
                        blnHeaderParsed = false;

                        // Create the output files
                        var strBaseOutputFilePath = Path.Combine(strOutputFolderPath, Path.GetFileName(strInputFilePath));
                        blnSuccess = InitializeSequenceOutputFiles(strBaseOutputFilePath);

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !AbortProcessing)
                        {
                            strLineIn = srDataFile.ReadLine();
                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseMSGFDBSynFileHeaderLine(strLineIn, out intColumnMapping);
                                if (!blnSuccess)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return blnSuccess;
                                }
                                blnHeaderParsed = true;
                                continue;
                            }

                            blnValidSearchResult = ParseMSGFDBSynFileEntry(strLineIn, objSearchResult, ref strErrorLog,
                                                                           intResultsProcessed, intColumnMapping,
                                                                           out strCurrentPeptideWithMods);

                            if (blnValidSearchResult)
                            {
                                strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.Charge;

                                if (objSearchResult.SpecEValue == strPreviousSpecEValue)
                                {
                                    // New result has the same SpecEValue as the previous result
                                    // See if peptidesFoundForSpecEValueLevel contains the peptide, scan and charge

                                    if (peptidesFoundForSpecEValueLevel.Contains(strKey))
                                    {
                                        blnFirstMatchForGroup = false;
                                    }
                                    else
                                    {
                                        peptidesFoundForSpecEValueLevel.Add(strKey);
                                        blnFirstMatchForGroup = true;
                                    }
                                }
                                else
                                {
                                    // New SpecEValue
                                    peptidesFoundForSpecEValueLevel.Clear();

                                    // Update strPreviousSpecEValue
                                    strPreviousSpecEValue = objSearchResult.SpecEValue;

                                    // Append a new entry
                                    peptidesFoundForSpecEValueLevel.Add(strKey);
                                    blnFirstMatchForGroup = true;
                                }

                                blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup);
                                if (!blnSuccess)
                                {
                                    successOverall = false;
                                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        strErrorLog += "Error adding modifications to sequence at RowIndex '" + objSearchResult.ResultID + "'" +
                                                       "\n";
                                    }
                                }

                                SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass) objSearchResult, blnFirstMatchForGroup);

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
                                                SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass) objSearchResult, false);
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
                        var inputFile = new FileInfo(strInputFilePath);
                        strModificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        strModificationSummaryFilePath = Path.Combine(strOutputFolderPath, strModificationSummaryFilePath);

                        SaveModificationSummaryFile(strModificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (strErrorLog.Length > 0)
                    {
                        SetErrorMessage("Invalid Lines: " + "\n" + strErrorLog);
                    }

                    blnSuccess = successOverall;
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
                SetErrorMessage(ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private bool ParseMSGFDBResultsFileEntry(
            string strLineIn,
            bool blnMSGFPlus,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstMSGFDBModInfo,
            ICollection<udtMSGFDBSearchResultType> lstSearchResultsCurrentScan,
            ref string strErrorLog,
            IList<int> intColumnMapping,
            ref int intNextScanGroupID,
            ICollection<udtScanGroupInfoType> lstScanGroupDetails,
            IDictionary<string, bool> htScanGroupCombo,
            IDictionary<string, int> lstSpecIdToIndex)
        {
            // Parses an entry from the MSGF+ results file

            var udtSearchResult = new udtMSGFDBSearchResultType();
            string[] strSplitLine = null;

            var intScanCount = 0;
            string[] strSplitResult = null;

            udtMSGFDBSearchResultType[] udtMergedScanInfo = null;

            var lstProteinInfo = default(Dictionary<string, udtTerminusCharsType>);

            var blnValidSearchResult = false;
            var intSlashIndex = 0;
            var blnTargetDecoyFDRValid = false;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                lstProteinInfo = new Dictionary<string, udtTerminusCharsType>();

                // Reset lstSearchResults
                lstSearchResultsCurrentScan.Clear();

                udtSearchResult.Clear();
                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 13)
                {
                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.SpectrumFile], out udtSearchResult.SpectrumFileName))
                    {
                        ReportError("SpectrumFile column is missing or invalid", true);
                    }
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.SpecIndex], out udtSearchResult.SpecIndex);

                    if (blnMSGFPlus)
                    {
                        var intSpecIndex = 0;
                        var blnGenerateSpecIndex = true;

                        if (!int.TryParse(udtSearchResult.SpecIndex, out intSpecIndex))
                        {
                            // MSGF+ includes text in the SpecID column, for example: "controllerType=0 controllerNumber=1 scan=6390" or "index=4323"
                            // Need to convert these to an integer

                            if (udtSearchResult.SpecIndex.StartsWith("index="))
                            {
                                udtSearchResult.SpecIndex = udtSearchResult.SpecIndex.Substring("index=".Length);
                                if (int.TryParse(udtSearchResult.SpecIndex, out intSpecIndex))
                                {
                                    blnGenerateSpecIndex = false;
                                }
                            }

                            if (blnGenerateSpecIndex)
                            {
                                if (!lstSpecIdToIndex.TryGetValue(udtSearchResult.SpecIndex, out intSpecIndex))
                                {
                                    intSpecIndex = lstSpecIdToIndex.Count + 1;
                                    lstSpecIdToIndex.Add(udtSearchResult.SpecIndex, intSpecIndex);
                                }

                                udtSearchResult.SpecIndex = intSpecIndex.ToString();
                            }
                        }
                    }

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.Scan], out udtSearchResult.Scan))
                    {
                        ReportError("Scan column is missing or invalid", true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                    intSlashIndex = udtSearchResult.Scan.IndexOf('/');
                    if (intSlashIndex > 0)
                    {
                        // This is a merged spectrum and thus scan number looks like: 3010/3011/3012
                        // Split the Scan list on the slash
                        // Later in this function, we'll append lstSearchResults with this scan plus the other scans

                        strSplitResult = udtSearchResult.Scan.Split('/');
                        intScanCount = strSplitResult.Length;
                        udtMergedScanInfo = new udtMSGFDBSearchResultType[intScanCount];

                        for (var intIndex = 0; intIndex <= intScanCount - 1; intIndex++)
                        {
                            udtMergedScanInfo[intIndex] = new udtMSGFDBSearchResultType();
                            udtMergedScanInfo[intIndex].Clear();
                            udtMergedScanInfo[intIndex].Scan = strSplitResult[intIndex];
                            udtMergedScanInfo[intIndex].ScanNum = CIntSafe(strSplitResult[intIndex], 0);
                        }

                        // Now split SpecIndex and store in udtMergedScanInfo
                        strSplitResult = udtSearchResult.SpecIndex.Split('/');

                        for (var intIndex = 0; intIndex <= strSplitResult.Length - 1; intIndex++)
                        {
                            if (intIndex >= udtMergedScanInfo.Length)
                            {
                                // There are more entries for SpecIndex than there are for Scan#; this is unexpected
                                break;
                            }
                            udtMergedScanInfo[intIndex].SpecIndex = strSplitResult[intIndex];
                        }

                        // Now split FragMethod and store in udtMergedScanInfo
                        strSplitResult = udtSearchResult.FragMethod.Split('/');

                        for (var intIndex = 0; intIndex <= strSplitResult.Length - 1; intIndex++)
                        {
                            if (intIndex >= udtMergedScanInfo.Length)
                            {
                                // There are more entries for FragMethod than there are for Scan#; this is unexpected
                                break;
                            }
                            udtMergedScanInfo[intIndex].FragMethod = strSplitResult[intIndex];
                        }
                    }
                    else
                    {
                        udtSearchResult.ScanNum = CIntSafe(udtSearchResult.Scan, 0);
                        intScanCount = 1;
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.PrecursorMZ], out udtSearchResult.PrecursorMZ);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    // Precursor mass error could be in PPM or Da
                    //   In MSGFDB, the header line will have PMError(ppm)        or PMError(Da)
                    //   In MSGF+,  the header line will have PrecursorError(ppm) or PrecursorError(Da)
                    double dblPrecursorErrorDa = 0;

                    if (intColumnMapping[(int)eMSGFDBResultsFileColumns.PMErrorPPM] >= 0)
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.PMErrorPPM], out udtSearchResult.PMErrorPPM);
                    }
                    else
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.PMErrorDa], out udtSearchResult.PMErrorDa);
                        dblPrecursorErrorDa = CDblSafe(udtSearchResult.PMErrorDa, 0);
                        udtSearchResult.PMErrorPPM = string.Empty;              // We'll populate this column later in this function
                    }

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                    {
                        ReportError("Peptide column is missing or invalid", true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.Protein], out udtSearchResult.Protein);

                    // MSGF+ .tsv files may have a semicolon separated list of protein names; check for this
                    udtSearchResult.Protein = SplitProteinList(udtSearchResult.Protein, lstProteinInfo);

                    if (lstProteinInfo.Count > 0)
                    {
                        // Need to add the prefix and suffix residues
                        udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, lstProteinInfo.First());
                    }

                    // Replace any mod text values in the peptide sequence with the appropriate mod symbols
                    // In addition, replace the terminus symbols with dashes
                    double dblTotalModMass = 0;
                    udtSearchResult.Peptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(udtSearchResult.Peptide), lstMSGFDBModInfo, blnMSGFPlus, out dblTotalModMass);

                    // Compute monoisotopic mass of the peptide
                    var dblPeptideMonoisotopicMass = ComputePeptideMass(udtSearchResult.Peptide, dblTotalModMass);

                    // Store the monoisotopic MH value in .MH
                    // This is (M+H)+ when the charge carrier is a proton
                    udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPeptideMonoisotopicMass, 0), 6);

                    if (!string.IsNullOrEmpty(udtSearchResult.PMErrorPPM))
                    {
                        // Convert the ppm-based PM Error to Da-based

                        double dblPMErrorPPM = 0;
                        double dblPrecursorMZ = 0;

                        if (double.TryParse(udtSearchResult.PrecursorMZ, out dblPrecursorMZ))
                        {
                            // Note that since .PMErrorPPM is present, the Precursor m/z is a C13-corrected m/z value
                            // In other words, it may not be the actual m/z selected for fragmentation.

                            if (double.TryParse(udtSearchResult.PMErrorPPM, out dblPMErrorPPM))
                            {
                                if (mParentMassToleranceInfo.IsPPM &&
                                    (dblPMErrorPPM < -mParentMassToleranceInfo.ToleranceLeft * 1.5 ||
                                     dblPMErrorPPM > mParentMassToleranceInfo.ToleranceRight * 1.5))
                                {
                                    // PPM error computed by MSGF+ is more than 1.5-fold larger than the ppm-based parent ion tolerance; don't trust the value computed by MSGF+

                                    mPrecursorMassErrorWarningCount += 1;
                                    if (mPrecursorMassErrorWarningCount <= 10)
                                    {
                                        ReportWarning("Precursor mass error computed by MSGF+ is 1.5-fold larger than search tolerance: " + udtSearchResult.PMErrorPPM + " vs. " + mParentMassToleranceInfo.ToleranceLeft.ToString("0") + "ppm," + mParentMassToleranceInfo.ToleranceRight.ToString("0") + "ppm");
                                        if (mPrecursorMassErrorWarningCount == 10)
                                        {
                                            ReportWarning("Additional mass errors will not be reported");
                                        }
                                    }

                                    double dblPrecursorMonoMass = 0;
                                    dblPrecursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, udtSearchResult.ChargeNum, 0);

                                    dblPrecursorErrorDa = dblPrecursorMonoMass - dblPeptideMonoisotopicMass;

                                    udtSearchResult.PMErrorPPM = string.Empty;
                                }
                                else
                                {
                                    dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPMErrorPPM, dblPeptideMonoisotopicMass);

                                    // Note that this will be a C13-corrected precursor error; not the true precursor error
                                    udtSearchResult.PMErrorDa = PRISM.StringUtilities.DblToString(dblPrecursorErrorDa, 6);
                                }
                            }
                        }
                    }

                    if (string.IsNullOrEmpty(udtSearchResult.PMErrorPPM))
                    {
                        double dblPrecursorMZ = 0;
                        if (double.TryParse(udtSearchResult.PrecursorMZ, out dblPrecursorMZ))
                        {
                            double dblPeptideDeltaMassCorrectedPpm = 0;

                            dblPeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMZ,
                                                                                      udtSearchResult.ChargeNum, dblPeptideMonoisotopicMass, true);

                            udtSearchResult.PMErrorPPM = PRISM.StringUtilities.DblToString(dblPeptideDeltaMassCorrectedPpm, 5);

                            if (string.IsNullOrEmpty(udtSearchResult.PMErrorDa))
                            {
                                dblPrecursorErrorDa = clsPeptideMassCalculator.PPMToMass(dblPeptideDeltaMassCorrectedPpm, dblPeptideMonoisotopicMass);

                                // Note that this will be a C13-corrected precursor error; not the true precursor error
                                udtSearchResult.PMErrorDa = PRISM.StringUtilities.DblToString(dblPrecursorErrorDa, 6);
                            }
                        }
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.DeNovoScore], out udtSearchResult.DeNovoScore);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.MSGFScore], out udtSearchResult.MSGFScore);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.SpecProb_EValue], out udtSearchResult.SpecEValue);
                    if (!double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum))
                        udtSearchResult.SpecEValueNum = 0;

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.PValue_EValue], out udtSearchResult.EValue);
                    if (!double.TryParse(udtSearchResult.EValue, out udtSearchResult.EValueNum))
                        udtSearchResult.EValueNum = 0;

                    blnTargetDecoyFDRValid = GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.FDR_QValue], out udtSearchResult.QValue);
                    if (!double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum))
                        udtSearchResult.QValueNum = 0;

                    if (blnTargetDecoyFDRValid)
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.PepFDR_PepQValue], out udtSearchResult.PepQValue);
                    }
                    else
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.EFDR], out udtSearchResult.QValue);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.IsotopeError], out udtSearchResult.IsotopeError);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.IMSScan], out udtSearchResult.IMSScan);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSGFDBResultsFileColumns.IMSDriftTime], out udtSearchResult.IMSDriftTime);

                    udtSearchResult.NTT = ComputeCleaveageState(udtSearchResult.Peptide).ToString();

                    var udtScanGroupInfo = default(udtScanGroupInfoType);
                    var intCurrentScanGroupID = -1;

                    udtScanGroupInfo.Charge = udtSearchResult.ChargeNum;

                    if (intScanCount > 1)
                    {
                        // Append one entry to lstSearchResults for each item in udtMergedScanInfo()

                        for (var intIndex = 0; intIndex <= intScanCount - 1; intIndex++)
                        {
                            udtSearchResult.Scan = udtMergedScanInfo[intIndex].Scan;
                            udtSearchResult.ScanNum = udtMergedScanInfo[intIndex].ScanNum;

                            udtSearchResult.SpecIndex = udtMergedScanInfo[intIndex].SpecIndex;
                            udtSearchResult.FragMethod = udtMergedScanInfo[intIndex].FragMethod;

                            AppendToSearchResults(lstSearchResultsCurrentScan, udtSearchResult, lstProteinInfo);

                            // Append an entry to lstScanGroupDetails
                            udtScanGroupInfo.Scan = udtSearchResult.ScanNum;
                            AppendToScanGroupDetails(lstScanGroupDetails, htScanGroupCombo, udtScanGroupInfo, ref intCurrentScanGroupID, ref intNextScanGroupID);
                        }
                    }
                    else
                    {
                        // This is not a merged result; simply append udtSearchResult to lstSearchResults
                        AppendToSearchResults(lstSearchResultsCurrentScan, udtSearchResult, lstProteinInfo);

                        // Also append an entry to lstScanGroupDetails
                        udtScanGroupInfo.Scan = udtSearchResult.ScanNum;
                        AppendToScanGroupDetails(lstScanGroupDetails, htScanGroupCombo, udtScanGroupInfo, ref intCurrentScanGroupID, ref intNextScanGroupID);
                    }

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MSGF+ results file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (strSplitLine != null && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing MSGF+ Results in ParseMSGFDBResultsFileEntry for RowIndex '" + strSplitLine[0] + "'" +
                                       "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MSGF+ Results in ParseMSGFDBResultsFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        private bool ParseMSGFDBResultsFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            // The expected header from MSGFDB is:
            // #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(Da)           Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR
            // or
            // #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(ppm)          Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR

            // The expected header from MSGF+ is:
            // #SpecFile    SpecID       ScanNum   FragMethod    Precursor    IsotopeError    PrecursorError(Da)    Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue
            // or
            // #SpecFile    SpecID       ScanNum   FragMethod    Precursor    IsotopeError    PrecursorError(ppm)   Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue

            var lstColumnNames = new SortedDictionary<string, eMSGFDBResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {"#SpecFile", eMSGFDBResultsFileColumns.SpectrumFile},
                {"SpecIndex", eMSGFDBResultsFileColumns.SpecIndex},
                {"SpecID", eMSGFDBResultsFileColumns.SpecIndex},
                {"Scan#", eMSGFDBResultsFileColumns.Scan},
                {"ScanNum", eMSGFDBResultsFileColumns.Scan},
                {"FragMethod", eMSGFDBResultsFileColumns.FragMethod},
                {"Precursor", eMSGFDBResultsFileColumns.PrecursorMZ},
                {"IsotopeError", eMSGFDBResultsFileColumns.IsotopeError},
                {"PMError(Da)", eMSGFDBResultsFileColumns.PMErrorDa},
                {"PrecursorError(Da)", eMSGFDBResultsFileColumns.PMErrorDa},
                {"PMError(ppm)", eMSGFDBResultsFileColumns.PMErrorPPM},
                {"PrecursorError(ppm)", eMSGFDBResultsFileColumns.PMErrorPPM},
                {"Charge", eMSGFDBResultsFileColumns.Charge},
                {"Peptide", eMSGFDBResultsFileColumns.Peptide},
                {"Protein", eMSGFDBResultsFileColumns.Protein},
                {"DeNovoScore", eMSGFDBResultsFileColumns.DeNovoScore},
                {"MSGFScore", eMSGFDBResultsFileColumns.MSGFScore},
                {"SpecProb", eMSGFDBResultsFileColumns.SpecProb_EValue},
                {"SpecEValue", eMSGFDBResultsFileColumns.SpecProb_EValue},
                {"P-value", eMSGFDBResultsFileColumns.PValue_EValue},
                {"EValue", eMSGFDBResultsFileColumns.PValue_EValue},
                {"FDR", eMSGFDBResultsFileColumns.FDR_QValue},
                {"QValue", eMSGFDBResultsFileColumns.FDR_QValue},
                {"PepFDR", eMSGFDBResultsFileColumns.PepFDR_PepQValue},
                {"PepQValue", eMSGFDBResultsFileColumns.PepFDR_PepQValue},
                {"EFDR", eMSGFDBResultsFileColumns.EFDR},
                {clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Scan, eMSGFDBResultsFileColumns.IMSScan},
                {clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Drift_Time, eMSGFDBResultsFileColumns.IMSDriftTime}
            };

            intColumnMapping = new int[MSGFDBResultsFileColCount];

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
                    else
                    {
                        // Unrecognized column name
                        Console.WriteLine("Warning: Unrecognized column header name '" + strSplitLine[intIndex] + "' in ParseMSGFDBResultsFileHeaderLine");
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MSGFDB results file: " + ex.Message, ex);
                return false;
            }

            return true;
        }

        private bool ParseMSGFDBSynFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            var lstColumnNames = new SortedDictionary<string, eMSFDBSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {clsPHRPParserMSGFDB.DATA_COLUMN_ResultID, eMSFDBSynFileColumns.ResultID},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Scan, eMSFDBSynFileColumns.Scan},
                {clsPHRPParserMSGFDB.DATA_COLUMN_FragMethod, eMSFDBSynFileColumns.FragMethod},
                {clsPHRPParserMSGFDB.DATA_COLUMN_SpecIndex, eMSFDBSynFileColumns.SpecIndex},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Charge, eMSFDBSynFileColumns.Charge},
                {clsPHRPParserMSGFDB.DATA_COLUMN_PrecursorMZ, eMSFDBSynFileColumns.PrecursorMZ},
                {clsPHRPParserMSGFDB.DATA_COLUMN_DelM, eMSFDBSynFileColumns.DelM},
                {clsPHRPParserMSGFDB.DATA_COLUMN_DelM_PPM, eMSFDBSynFileColumns.DelMPPM},
                {clsPHRPParserMSGFDB.DATA_COLUMN_MH, eMSFDBSynFileColumns.MH},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Peptide, eMSFDBSynFileColumns.Peptide},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Protein, eMSFDBSynFileColumns.Protein},
                {clsPHRPParserMSGFDB.DATA_COLUMN_NTT, eMSFDBSynFileColumns.NTT},
                {clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore, eMSFDBSynFileColumns.DeNovoScore},
                {clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore, eMSFDBSynFileColumns.MSGFScore},
                {clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb, eMSFDBSynFileColumns.SpecProb_EValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue, eMSFDBSynFileColumns.SpecProb_EValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecProb, eMSFDBSynFileColumns.RankSpecProb},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue, eMSFDBSynFileColumns.RankSpecProb},
                {clsPHRPParserMSGFDB.DATA_COLUMN_PValue, eMSFDBSynFileColumns.PValue_EValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_EValue, eMSFDBSynFileColumns.PValue_EValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_FDR, eMSFDBSynFileColumns.FDR_QValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_QValue, eMSFDBSynFileColumns.FDR_QValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR, eMSFDBSynFileColumns.PepFDR_PepQValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_PepQValue, eMSFDBSynFileColumns.PepFDR_PepQValue},
                {clsPHRPParserMSGFDB.DATA_COLUMN_EFDR, eMSFDBSynFileColumns.EFDR},
                {clsPHRPParserMSGFDB.DATA_COLUMN_Isotope_Error, eMSFDBSynFileColumns.IsotopeError}
            };

            intColumnMapping = new int[MSGFDBSynFileColCount];

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
                SetErrorMessage("Error parsing header in MSGFDB synopsis file: " + ex.Message, ex);
                return false;
            }

            return true;
        }

        private bool ParseMSGFDBSynFileEntry(
            string strLineIn,
            clsSearchResultsMSGFDB objSearchResult,
            ref string strErrorLog,
            int intResultsProcessed,
            IReadOnlyList<int> intColumnMapping,
            out string strPeptideSequenceWithMods)
        {
            // Parses an entry from the MSGFDB Synopsis file

            string[] strSplitLine = null;

            bool blnValidSearchResult;

            // Reset objSearchResult
            objSearchResult.Clear();
            strPeptideSequenceWithMods = string.Empty;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 15)
                {
                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.ResultID], out string strValue))
                    {
                        if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                        {
                            strErrorLog += "Error reading ResultID value from MSGFDB Results line " + (intResultsProcessed + 1) +
                                           "\n";
                        }
                        return false;
                    }

                    objSearchResult.ResultID = int.Parse(strValue);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.Scan], out string scan);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.Charge], out string charge);

                    objSearchResult.Scan = scan;
                    objSearchResult.Charge = charge;

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.Peptide], out strPeptideSequenceWithMods))
                    {
                        if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                        {
                            strErrorLog += "Error reading Peptide sequence value from MSGFDB Results line " + (intResultsProcessed + 1) +
                                           "\n";
                        }
                        return false;
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.Protein], out string proteinName);
                    objSearchResult.MultipleProteinCount = "0";

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.DelM], out string msgfPlusComputedDelM);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.DelMPPM], out string msgfPlusComputedDelMppm);

                    objSearchResult.ProteinName = proteinName;
                    objSearchResult.MSGFPlusComputedDelM = msgfPlusComputedDelM;
                    objSearchResult.MSGFPlusComputedDelMPPM = msgfPlusComputedDelMppm;

                    objSearchResult.PeptideDeltaMass = objSearchResult.MSGFPlusComputedDelM;

                    // Note: .PeptideDeltaMass is stored in the MSGF+ results file as "Observed_Mass - Theoretical_Mass"
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
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.FragMethod], out string fragMethod);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.PrecursorMZ], out string precursorMz);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.MH], out string peptideMh);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.NTT], out string ntt);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.DeNovoScore], out string deNovoScore);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.MSGFScore], out string msgfScore);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.SpecProb_EValue], out string specEValue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.RankSpecProb], out string rankSpecEValue);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.PValue_EValue], out string eValue);

                    var blnTargetDecoyFDRValid = GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.FDR_QValue], out string qValue);

                    objSearchResult.FragMethod = fragMethod;
                    objSearchResult.PrecursorMZ = precursorMz;
                    objSearchResult.PeptideMH = peptideMh;
                    objSearchResult.NTT = ntt;
                    objSearchResult.DeNovoScore = deNovoScore;
                    objSearchResult.MSGFScore = msgfScore;
                    objSearchResult.SpecEValue = specEValue;
                    objSearchResult.RankSpecEValue = rankSpecEValue;
                    objSearchResult.EValue = eValue;
                    objSearchResult.QValue = qValue;

                    if (blnTargetDecoyFDRValid)
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.PepFDR_PepQValue], out string pepQValue);
                        objSearchResult.PepQValue = pepQValue;
                    }
                    else
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.EFDR], out string efdr);
                        objSearchResult.QValue = efdr;
                    }

                    if (intColumnMapping[(int)eMSFDBSynFileColumns.IsotopeError] >= 0)
                    {
                        GetColumnValue(strSplitLine, intColumnMapping[(int)eMSFDBSynFileColumns.IsotopeError], out string isotopeError);
                        objSearchResult.IsotopeError = isotopeError;
                        objSearchResult.MSGFPlusResults = true;
                    }
                    else
                    {
                        objSearchResult.MSGFPlusResults = false;
                    }

                    // Compute PrecursorMH using PrecursorMZ and charge
                    if (double.TryParse(objSearchResult.PrecursorMZ, out var dblPrecursorMZ))
                    {
                        if (int.TryParse(objSearchResult.Charge, out var intCharge))
                        {
                            objSearchResult.ParentIonMH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, intCharge), 6);
                        }
                    }

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (strSplitLine != null && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing MSGFDB Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MSGFDB Results in ParseMSGFDBSynFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        private bool ParseParentMassTolerance(string strToleranceText, out double dblTolerance, out bool blnIsPPM)
        {
            dblTolerance = 0;
            blnIsPPM = false;

            strToleranceText = strToleranceText.ToLower().Trim();

            if (strToleranceText.EndsWith("da"))
            {
                strToleranceText = strToleranceText.Substring(0, strToleranceText.Length - 2);
                blnIsPPM = false;
            }
            else if (strToleranceText.EndsWith("ppm"))
            {
                strToleranceText = strToleranceText.Substring(0, strToleranceText.Length - 3);
                blnIsPPM = true;
            }
            else
            {
                return false;
            }

            if (double.TryParse(strToleranceText, out dblTolerance))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">MSGFDB results file</param>
        /// <param name="strOutputFolderPath">Output folder</param>
        /// <param name="strParameterFilePath">Parameter file for data processing</param>
        /// <returns>True if success, False if failure</returns>
        /// <remarks>Use SearchToolParameterFilePath to define the search engine parameter file</remarks>
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            var strFhtOutputFilePath = string.Empty;
            string strSynOutputFilePath = null;
            string strPepToProteinMapFilePath = null;
            string strScanGroupFilePath = null;

            var lstMSGFDBModInfo = default(List<clsMSGFPlusParamFileModExtractor.udtModInfoType>);
            var lstPepToProteinMapping = default(List<udtPepToProteinMappingType>);
            var strMTSPepToProteinMapFilePath = string.Empty;

            var blnMSGFPlus = false;
            var lstSpecIdToIndex = default(Dictionary<string, int>);

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

                mPeptideSeqMassCalculator.ResetAminoAcidMasses();

                lstSpecIdToIndex = new Dictionary<string, int>();

                ResetProgress("Parsing " + Path.GetFileName(strInputFilePath));

                if (!CleanupFilePaths(ref strInputFilePath, ref strOutputFolderPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(strInputFilePath);

                    lstMSGFDBModInfo = new List<clsMSGFPlusParamFileModExtractor.udtModInfoType>();
                    lstPepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the MSGF+ Parameter File so that we can determine the modification names and masses
                    // If the MSGFPlus_Mods.txt or MSGFDB_Mods.txt file was defined, the mod symbols in that file will be used to define the mod symbols in lstMSGFDBModInfo
                    var success = ExtractModInfoFromParamFile(SearchToolParameterFilePath, out lstMSGFDBModInfo);
                    if (!success)
                    {
                        return false;
                    }

                    if (!LoadSearchEngineParamFile(SearchToolParameterFilePath))
                    {
                        return false;
                    }

                    var query = from item in lstMSGFDBModInfo where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFDBModType.CustomAA select item;
                    if (query.Any())
                    {
                        // Custom amino acids are defined; read their values and update the mass calculator

                        var modFileProcessor = new clsMSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME);
                        RegisterEvents(modFileProcessor);

                        modFileProcessor.ErrorEvent += ModExtractorErrorHandler;


                        clsPHRPParserMSGFDB.UpdateMassCalculatorMasses(SearchToolParameterFilePath, modFileProcessor, mPeptideSeqMassCalculator,
                            out var localErrorMsg);

                        if (!string.IsNullOrWhiteSpace(localErrorMsg) && string.IsNullOrWhiteSpace(mErrorMessage))
                        {
                            ReportError(localErrorMsg);
                        }
                    }

                    // Define the base output filename using strInputFilePath
                    var strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath);

                    // Auto-replace "_msgfdb" with "_msgfplus"
                    if (strBaseName.EndsWith("_msgfdb", StringComparison.InvariantCultureIgnoreCase))
                    {
                        strBaseName = strBaseName.Substring(0, strBaseName.Length - "_msgfdb".Length) + "_msgfplus";
                    }

                    if (CreateInspectFirstHitsFile)
                    {
                        // Read the FASTA file to cache the protein names in memory
                        // These will be used when creating the first hits file
                        if (!CacheProteinNamesFromFasta())
                        {
                            return false;
                        }

                        // Create the first hits output file
                        ResetProgress("Creating the FHT file", true);

                        strFhtOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName + SEQUEST_FIRST_HITS_FILE_SUFFIX);

                        strScanGroupFilePath = string.Empty;

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strFhtOutputFilePath, strScanGroupFilePath, lstMSGFDBModInfo, out blnMSGFPlus, lstSpecIdToIndex, eFilteredOutputFileTypeConstants.FHTFile);
                    }

                    if (CreateInspectSynopsisFile)
                    {
                        // Create the synopsis output file
                        ResetProgress("Creating the SYN file", true);

                        // The synopsis file name will be of the form BasePath_msgfplus_syn.txt
                        strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                        strScanGroupFilePath = Path.Combine(strOutputFolderPath, strBaseName + "_ScanGroupInfo.txt");

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, strScanGroupFilePath, lstMSGFDBModInfo, out blnMSGFPlus, lstSpecIdToIndex, eFilteredOutputFileTypeConstants.SynFile);

                        // Load the PeptideToProteinMap information; if the file doesn't exist, a warning will be displayed, but processing will continue
                        // LoadPeptideToProteinMapInfoMSGFDB also creates _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
                        strPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(Path.Combine(strOutputFolderPath, strBaseName) + ".txt", strOutputFolderPath, MTS: false);

                        ResetProgress("Loading the PepToProtein map file: " + Path.GetFileName(strPepToProteinMapFilePath), true);

                        LoadPeptideToProteinMapInfoMSGFDB(strPepToProteinMapFilePath, strOutputFolderPath, lstMSGFDBModInfo, blnMSGFPlus, lstPepToProteinMapping, out strMTSPepToProteinMapFilePath);

                        // Create the other PHRP-specific files
                        ResetProgress("Creating the PHRP files for " + Path.GetFileName(strSynOutputFilePath), true);

                        // Now parse the _syn.txt file that we just created to create the other PHRP files
                        blnSuccess = ParseMSGFDBSynopsisFile(strSynOutputFilePath, strOutputFolderPath, lstPepToProteinMapping, false);

                        // Remove all items from lstPepToProteinMapping to reduce memory overhead
                        lstPepToProteinMapping.Clear();
                        lstPepToProteinMapping.TrimExcess();

                        if (blnSuccess && CreateProteinModsFile)
                        {
                            blnSuccess = CreateProteinModsFileWork(strBaseName, inputFile, strFhtOutputFilePath, strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath);
                        }
                    }

                    if (blnSuccess)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in clsMSGFDBResultsProcessor.ProcessFile (2):  " + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in clsMSGFDBResultsProcessor.ProcessFile (1):" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return blnSuccess;
        }

        private bool CreateProteinModsFileWork(
            string strBaseName,
            FileInfo inputFile,
            string strFhtOutputFilePath,
            string strSynOutputFilePath,
            string strOutputFolderPath,
            string strMTSPepToProteinMapFilePath)
        {
            var blnSuccess = true;

            if (string.IsNullOrEmpty(strMTSPepToProteinMapFilePath) || !File.Exists(strMTSPepToProteinMapFilePath))
            {
                // MTSPepToProteinMap file not found; auto-create it

                if (string.IsNullOrEmpty(strMTSPepToProteinMapFilePath))
                {
                    strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strBaseName, strOutputFolderPath, MTS: true);
                }

                var lstSourcePHRPDataFiles = new List<string>();

                if (!string.IsNullOrEmpty(strFhtOutputFilePath))
                {
                    lstSourcePHRPDataFiles.Add(strFhtOutputFilePath);
                }

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
                        // We only do this for MSGFDB since it often includes reverse protein peptides in the results even though the FASTA file often does not have reverse proteins
                        IgnorePeptideToProteinMapperErrors = true;
                        blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                        if (!blnSuccess)
                        {
                            ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                        }
                    }
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath);

                // Create the Protein Mods file
                blnSuccess = CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MSGFDB);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                blnSuccess = true;
            }

            return blnSuccess;
        }

        private static readonly Regex RegexNTerminalModMassRegEx = new Regex(MSGFDB_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);
        private static readonly Regex RegexModMassRegEx = new Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Replaces modification masses in peptide sequences with modification symbols (uses case-sensitive comparisons)
        /// </summary>
        /// <param name="strPeptide"></param>
        /// <param name="lstMSGFDBModInfo">This function assumes that each entry in lstMSGFDBModInfo has both .ModName and .ModSymbol defined</param>
        /// <param name="blnMSGFPlus">Should be set to True if processing MSGF+ results</param>
        /// <param name="dblTotalModMass">Output parameter: total mass of all modifications</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ReplaceMSGFModTextWithSymbol(
            string strPeptide,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> lstMSGFDBModInfo,
            bool blnMSGFPlus,
            out double dblTotalModMass)
        {
            var intIndex = 0;
            var intIndexFirstResidue = 0;
            var intIndexLastResidue = 0;

            var strPrefix = string.Empty;
            var strSuffix = string.Empty;

            var strModSymbols = string.Empty;
            var strDynModSymbols = string.Empty;

            var blnNterminalMod = false;
            var blnPossibleCTerminalMod = false;
            var blnContainsStaticMod = false;

            var reMatch = default(Match);

            double dblModMassFound = 0;

            // Reset the total mod mass
            dblTotalModMass = 0;

            // Remove the prefix and suffix residues
            if (strPeptide.Length >= 4)
            {
                if (strPeptide[1] == '.' &&
                    strPeptide[strPeptide.Length - 2] == '.')
                {
                    strPrefix = strPeptide.Substring(0, 2);
                    strSuffix = strPeptide.Substring(strPeptide.Length - 2, 2);

                    strPeptide = strPeptide.Substring(2, strPeptide.Length - 4);
                }
            }

            // strPeptide should now be the primary peptide sequence, without the prefix or suffix residues

            // First look for dynamic N-terminal mods (NTermPeptide or NTermProtein)
            // This RegEx will match one or more mods, all at the N-terminus
            reMatch = RegexNTerminalModMassRegEx.Match(strPeptide);

            if (reMatch != null && reMatch.Success)
            {
                blnNterminalMod = true;
                blnPossibleCTerminalMod = false;
                blnContainsStaticMod = false;

                // Convert the mod mass (or masses) to one or more mod symbols

                if (ConvertMGSFModMassesToSymbols("-", reMatch.Groups[1].Value, out strModSymbols, out strDynModSymbols, lstMSGFDBModInfo, blnNterminalMod, blnPossibleCTerminalMod, out dblModMassFound, out blnContainsStaticMod))
                {
                    // Replace the mod digits with the mod symbols

                    strPeptide = ReplaceMSGFModTextWithMatchedSymbol(strPeptide, reMatch.Groups[1], strModSymbols, strDynModSymbols, blnMSGFPlus, blnContainsStaticMod);
                    dblTotalModMass += dblModMassFound;
                }
            }

            // Next, step through the peptide and parse each mod mass that follows a residue
            // Any mod mass at the end must be considered a C-terminal mod

            // Need to start at the first letter
            // If we had N-terminal mods, they're currently notated like this: _.+42.011MDHTPQSQLK.L or _.+42.011+57.021MNDR.Q
            // We want things to look like this: -.#MDHTPQSQLK.L or -.#*MNDRQLNHR.S

            // In MSGFDB, static mods do not have a mod mass listed
            // In MSGF+,  static mods do have a mod mass listed
            // Regardless, we do not add mod symbols for static mods, but we do increment dblTotalModMass

            // Find the index of the last residue
            intIndex = strPeptide.Length - 1;
            while (intIndex > 0 && !IsLetterAtoZ(strPeptide[intIndex]))
            {
                intIndex -= 1;
            }
            intIndexLastResidue = intIndex;

            // Find the index of the first residue
            intIndex = 0;
            while (intIndex < strPeptide.Length && !IsLetterAtoZ(strPeptide[intIndex]))
            {
                intIndex += 1;
            }
            intIndexFirstResidue = intIndex;

            var currentResidue = "-";

            while (intIndex < strPeptide.Length)
            {
                if (IsLetterAtoZ(strPeptide[intIndex]))
                {
                    currentResidue = strPeptide[intIndex].ToString();

                    var objModificationDefinition = default(clsModificationDefinition);

                    if (!blnMSGFPlus)
                    {
                        // Look for static mods that should be applied to this residue (only applies to MSGFDB, not MSGF+)
                        for (var intModIndex = 0; intModIndex <= mPeptideMods.ModificationCount - 1; intModIndex++)
                        {
                            var eModificationType = default(clsModificationDefinition.eModificationTypeConstants);
                            eModificationType = mPeptideMods.GetModificationTypeByIndex(intModIndex);

                            if (eModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                            {
                                objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);

                                if (objModificationDefinition.TargetResiduesContain(strPeptide[intIndex]))
                                {
                                    // Match found; update dblTotalModMass but do not add a static mod symbol
                                    dblTotalModMass += objModificationDefinition.ModificationMass;
                                }
                            }
                            else if (intIndex == intIndexFirstResidue)
                            {
                                if (eModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod && strPrefix == "_")
                                {
                                    // N-terminal protein static mod
                                    objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);
                                    dblTotalModMass += objModificationDefinition.ModificationMass;
                                }
                                else if (eModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod)
                                {
                                    // N-terminal peptide static mod
                                    objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex);
                                    dblTotalModMass += objModificationDefinition.ModificationMass;
                                }
                            }
                        }
                    }

                    intIndex += 1;

                    if (intIndex == intIndexLastResidue)
                        blnPossibleCTerminalMod = true;
                }
                else
                {
                    // Found a mod; find the extent of the mod digits
                    reMatch = RegexModMassRegEx.Match(strPeptide, intIndex);

                    // Note that blnPossibleCTerminalMod will be set to True once we hit the last residue
                    // Assure blnNterminalMod is false
                    blnNterminalMod = false;

                    // Convert the mod mass (or masses) to one or more mod symbols

                    if (ConvertMGSFModMassesToSymbols(currentResidue, reMatch.Groups[1].Value, out strModSymbols, out strDynModSymbols, lstMSGFDBModInfo,
                        blnNterminalMod, blnPossibleCTerminalMod, out dblModMassFound, out blnContainsStaticMod))
                    {
                        strPeptide = ReplaceMSGFModTextWithMatchedSymbol(strPeptide, reMatch.Groups[1], strModSymbols, strDynModSymbols, blnMSGFPlus, blnContainsStaticMod);
                        dblTotalModMass += dblModMassFound;

                        if (blnMSGFPlus && blnContainsStaticMod)
                        {
                            // MSGF+ shows mod masses for static mods
                            // Thus, we have removed the static mod mass and did not add a mod symbol
                            // Therefore, leave intIndex unchanged
                        }
                        else
                        {
                            intIndex += strModSymbols.Length;
                        }
                    }
                    else
                    {
                        intIndex += reMatch.Groups[1].Value.Length;
                    }
                }
            }

            // If any N-terminal mods were present, we need to move them to after the first residue
            // in other words, switch from #MDHTPQSQLK to M#DHTPQSQLK
            //                          or #*MNDRQLNHR to M#*NDRQLNHR

            // Update intIndexFirstResidue
            intIndexFirstResidue = 0;
            while (intIndexFirstResidue < strPeptide.Length && !IsLetterAtoZ(strPeptide[intIndexFirstResidue]))
            {
                intIndexFirstResidue += 1;
            }

            if (intIndexFirstResidue > 0 && intIndexFirstResidue < strPeptide.Length)
            {
                string strPeptideNew = null;
                strPeptideNew = strPeptide[intIndexFirstResidue] + strPeptide.Substring(0, intIndexFirstResidue);
                if (intIndexFirstResidue < strPeptide.Length - 1)
                {
                    strPeptideNew += strPeptide.Substring(intIndexFirstResidue + 1);
                }
                strPeptide = string.Copy(strPeptideNew);
            }

            return strPrefix + strPeptide + strSuffix;
        }

        private string ReplaceMSGFModTextWithMatchedSymbol(
            string strPeptide,
            Capture reGroup,
            string strModSymbols,
            string strDynModSymbols,
            bool blnMSGFPlus,
            bool blnContainsStaticMod)
        {
            string strPeptideNew = null;

            if (reGroup.Index > 0)
            {
                strPeptideNew = strPeptide.Substring(0, reGroup.Index);
            }
            else
            {
                strPeptideNew = string.Empty;
            }

            if (blnMSGFPlus && blnContainsStaticMod)
            {
                // MSGF+ shows mod masses for static mods
                // However, for consistency with other PHRP results, we do not add a symbol to the peptide for this static mod
                // Catch: If we have a peptide/terminus affected by both a static and a dynamic mod, we still want the dynamic mod.
                if (!string.IsNullOrWhiteSpace(strDynModSymbols))
                {
                    strPeptideNew += strDynModSymbols;
                }
            }
            else
            {
                strPeptideNew += strModSymbols;
            }

            if (reGroup.Index + reGroup.Length < strPeptide.Length)
            {
                strPeptideNew += strPeptide.Substring(reGroup.Index + reGroup.Length);
            }

            return strPeptideNew;
        }

        private string ReplaceTerminus(string strPeptide)
        {
            if (strPeptide.StartsWith(N_TERMINUS_SYMBOL_MSGFDB))
            {
                strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + strPeptide.Substring(N_TERMINUS_SYMBOL_MSGFDB.Length);
            }

            if (strPeptide.EndsWith(C_TERMINUS_SYMBOL_MSGFDB))
            {
                strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_MSGFDB.Length) + "." + clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return strPeptide;
        }

        private static readonly Regex RegexProteinInfo = new Regex(PROTEIN_AND_TERM_SYMBOLS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Examines strProteinList to look for a semi-colon separated list of proteins and terminus symbols, for example
        /// AT1G26570.1(pre=K,post=N);AT3G29360.1(pre=K,post=N);AT3G29360.2(pre=K,post=N)
        /// </summary>
        /// <param name="strProteinList">Protein list to examine</param>
        /// <param name="lstProteinInfo">Protein information, if it is of the form ProteinName(pre=X,post=Y)</param>
        /// <returns>The name of the first protein</returns>
        /// <remarks></remarks>
        private string SplitProteinList(string strProteinList, IDictionary<string, udtTerminusCharsType> lstProteinInfo)
        {
            var reMatches = default(MatchCollection);

            lstProteinInfo.Clear();

            reMatches = RegexProteinInfo.Matches(strProteinList);

            if (reMatches.Count == 0)
            {
                // No match; likely just one protein
                return TruncateProteinName(strProteinList);
            }
            else
            {
                foreach (Match reMatch in reMatches)
                {
                    string strProteinName = null;
                    var udtTerminusChars = default(udtTerminusCharsType);

                    strProteinName = TruncateProteinName(reMatch.Groups[1].Value);
                    udtTerminusChars.NTerm = reMatch.Groups[2].Value[0];
                    udtTerminusChars.CTerm = reMatch.Groups[3].Value[0];

                    if (lstProteinInfo.ContainsKey(strProteinName))
                    {
                        // Skip this protein since it's already present
                    }
                    else
                    {
                        lstProteinInfo.Add(strProteinName, udtTerminusChars);
                    }
                }

                return lstProteinInfo.First().Key;
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter swResultFile,
            List<udtMSGFDBSearchResultType> lstFilteredSearchResults,
            ref string strErrorLog,
            bool blnIncludeFDRandPepFDR,
            bool blnIncludeEFDR,
            bool blnIncludeIMSFields,
            bool blnMSGFPlus)
        {
            // Sort udtFilteredSearchResults by ascending SpecEValue, ascending scan, ascending charge, ascending peptide, and ascending protein
            lstFilteredSearchResults.Sort(new MSGFDBSearchResultsComparerSpecEValueScanChargePeptide());

            for (var intIndex = 0; intIndex <= lstFilteredSearchResults.Count - 1; intIndex++)
            {
                WriteSearchResultToFile(intIndex + 1, swResultFile, lstFilteredSearchResults[intIndex], ref strErrorLog, blnIncludeFDRandPepFDR, blnIncludeEFDR, blnIncludeIMSFields, blnMSGFPlus);
            }
        }

        private void StoreScanGroupInfo(string strScanGroupFilePath, IReadOnlyCollection<udtScanGroupInfoType> lstScanGroupDetails)
        {
            var intScanGroupIDPrevious = 0;
            var blnCreateFile = false;

            try
            {
                // Only create the ScanGroup file if one or more scan groups exist
                // Step through lstScanGroupDetails to check for this
                intScanGroupIDPrevious = -1;
                blnCreateFile = false;
                foreach (var udtScanGroupInfo in lstScanGroupDetails)
                {
                    if (udtScanGroupInfo.ScanGroupID == intScanGroupIDPrevious)
                    {
                        blnCreateFile = true;
                        break;
                    }
                    intScanGroupIDPrevious = udtScanGroupInfo.ScanGroupID;
                }

                if (blnCreateFile)
                {
                    using (var swScanGroupFile = new StreamWriter(new FileStream(strScanGroupFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        swScanGroupFile.WriteLine("Scan_Group_ID" + "\t" + "Charge" + "\t" + "Scan");

                        foreach (var udtScanGroupInfo in lstScanGroupDetails)
                        {
                            swScanGroupFile.WriteLine(udtScanGroupInfo.ScanGroupID + "\t" + udtScanGroupInfo.Charge + "\t" + udtScanGroupInfo.Scan);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating ScanGroupInfo file: " + ex.Message);
            }
        }

        /// <summary>
        /// Stores the first hits file matches for a single scan
        /// </summary>
        /// <param name="lstSearchResults">Search results</param>
        /// <param name="intStartIndex">Start index for data in this scan</param>
        /// <param name="intEndIndex">End index for data in this scan</param>
        /// <param name="lstFilteredSearchResults">Filtered search results</param>
        /// <remarks></remarks>
        private void StoreTopFHTMatch(
            IList<udtMSGFDBSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex,
            List<udtMSGFDBSearchResultType> lstFilteredSearchResults)
        {
            AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex);

            // The calling procedure should have already sorted by scan, charge, and SpecEValue; no need to re-sort

            ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1);

            // Now store the first match for each charge for this scan
            // When storing, we use the protein name that occurred first in the FASTA file

            var udtCurrentResult = lstSearchResults[intStartIndex];
            var intCurrentCharge = udtCurrentResult.ChargeNum;
            var currentProteinNumber = int.MaxValue;
            var currentPeptide = GetCleanSequence(udtCurrentResult.Peptide);

            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                if (intCurrentCharge != lstSearchResults[intIndex].ChargeNum)
                {
                    // New charge state
                    // Store udtCurrentResult (from the previous charge state)
                    lstFilteredSearchResults.Add(udtCurrentResult);

                    udtCurrentResult = lstSearchResults[intIndex];
                    intCurrentCharge = udtCurrentResult.ChargeNum;
                    currentProteinNumber = int.MaxValue;
                    currentPeptide = GetCleanSequence(udtCurrentResult.Peptide);
                }

                var newPeptide = GetCleanSequence(lstSearchResults[intIndex].Peptide);
                if (currentPeptide.Equals(newPeptide))
                {
                    var bestProtein = GetBestProteinName(udtCurrentResult.Protein, currentProteinNumber, lstSearchResults[intIndex].Protein);
                    if (bestProtein.Value < currentProteinNumber)
                    {
                        currentProteinNumber = bestProtein.Value;
                        if (!udtCurrentResult.Protein.Equals(bestProtein.Key))
                        {
                            udtCurrentResult.Protein = bestProtein.Key;
                        }
                    }
                }
            }

            // Store udtCurrentResult (from the previous charge state)
            lstFilteredSearchResults.Add(udtCurrentResult);
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan
        /// </summary>
        /// <param name="lstSearchResults">Search results</param>
        /// <param name="intStartIndex">Start index for data in this scan</param>
        /// <param name="intEndIndex">End index for data in this scan</param>
        /// <param name="lstFilteredSearchResults">Filtered search results</param>
        /// <remarks></remarks>
        private void StoreSynMatches(
            IList<udtMSGFDBSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex,
            List<udtMSGFDBSearchResultType> lstFilteredSearchResults)
        {
            AssignRankAndDeltaNormValues(lstSearchResults, intStartIndex, intEndIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            ExpandListIfRequired(lstFilteredSearchResults, intEndIndex - intStartIndex + 1);

            var results = new SortedSet<string>();

            // Now store or write out the matches that pass the filters
            // By default, filter passing peptides have MSGFDB_SpecEValue <= 5E-7 Or EValue less than 0.75 or QValue less than 1% (but not 0)
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                if (MSGFPlusResultPassesSynFilter(lstSearchResults[intIndex]))
                {
                    // Check for identical results
                    // This can happen if the search used an n-terminal dynamic mod that also could apply to a specific residue and that residue is at the N-terminus
                    // For example:
                    //   R.S+229.163IGLPDVHSGYGFAIGNMAAFDMNDPEAVVSPGGVGFDINC+57.021GVR.L
                    //   R.+229.163SIGLPDVHSGYGFAIGNMAAFDMNDPEAVVSPGGVGFDINC+57.021GVR.L
                    //
                    // This software will have turned both of these results info:
                    //  R.S#IGLPDVHSGYGFAIGNMAAFDMNDPEAVVSPGGVGFDINCGVR.L
                    // We only want to include the result once in the _syn.txt file
                    // (though if the peptide maps to multiple proteins it will be listed multiple times; one line per protein)

                    var resultKey = lstSearchResults[intIndex].Peptide + "_" +
                                    lstSearchResults[intIndex].Protein + "_" +
                                    lstSearchResults[intIndex].MH + "_" +
                                    lstSearchResults[intIndex].SpecEValue;

                    if (results.Contains(resultKey))
                    {
                        continue;
                    }

                    results.Add(resultKey);
                    lstFilteredSearchResults.Add(lstSearchResults[intIndex]);
                }
            }
        }

        private void WriteSynFHTFileHeader(
            TextWriter swResultFile,
            ref string strErrorLog,
            bool blnIncludeFDRandPepFDR,
            bool blnIncludeEFDR,
            bool blnIncludeIMSFields,
            bool blnMSGFPlus)
        {
            // Write out the header line for synopsis / first hits files
            try
            {
                var lstData = new List<string>
                {
                    clsPHRPParserMSGFDB.DATA_COLUMN_ResultID,
                    clsPHRPParserMSGFDB.DATA_COLUMN_Scan,
                    clsPHRPParserMSGFDB.DATA_COLUMN_FragMethod,
                    clsPHRPParserMSGFDB.DATA_COLUMN_SpecIndex,
                    clsPHRPParserMSGFDB.DATA_COLUMN_Charge,
                    clsPHRPParserMSGFDB.DATA_COLUMN_PrecursorMZ,
                    clsPHRPParserMSGFDB.DATA_COLUMN_DelM,
                    clsPHRPParserMSGFDB.DATA_COLUMN_DelM_PPM,
                    clsPHRPParserMSGFDB.DATA_COLUMN_MH,
                    clsPHRPParserMSGFDB.DATA_COLUMN_Peptide,
                    clsPHRPParserMSGFDB.DATA_COLUMN_Protein,
                    clsPHRPParserMSGFDB.DATA_COLUMN_NTT,
                    clsPHRPParserMSGFDB.DATA_COLUMN_DeNovoScore,
                    clsPHRPParserMSGFDB.DATA_COLUMN_MSGFScore
                };

                if (blnMSGFPlus)
                {
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue);
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue);
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EValue);
                }
                else
                {
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb);
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecProb);
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PValue);
                }

                if (blnIncludeFDRandPepFDR)
                {
                    if (blnMSGFPlus)
                    {
                        lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_QValue);
                        lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepQValue);
                    }
                    else
                    {
                        lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_FDR);
                        lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR);
                    }
                }
                else if (blnIncludeEFDR)
                {
                    // Note that we'll write out a "1" for "PepFDR" for every result
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EFDR);
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR);
                }

                if (blnMSGFPlus)
                {
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Isotope_Error);
                }

                if (blnIncludeIMSFields)
                {
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Scan);
                    lstData.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Drift_Time);
                }

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
        /// <param name="blnIncludeFDRandPepFDR"></param>
        /// <param name="blnIncludeEFDR"></param>
        /// <param name="blnIncludeIMSFields"></param>
        /// <param name="blnMSGFPlus"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int intResultID,
            TextWriter swResultFile,
            udtMSGFDBSearchResultType udtSearchResult,
            ref string strErrorLog,
            bool blnIncludeFDRandPepFDR,
            bool blnIncludeEFDR,
            bool blnIncludeIMSFields,
            bool blnMSGFPlus)
        {
            try
            {
                // Primary Columns (other columns are added in certain circumstances):
                //
                // MSGFDB
                // ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecProb    Rank_MSGFDB_SpecProb    PValue  FDR     PepFDR

                // MSGF+
                // ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecEValue  Rank_MSGFDB_SpecEValue  EValue  QValue  PepQValue  IsotopeError

                var lstData = new List<string>
                {
                    intResultID.ToString(),
                    udtSearchResult.Scan,
                    udtSearchResult.FragMethod,
                    udtSearchResult.SpecIndex,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.PMErrorDa,
                    udtSearchResult.PMErrorPPM,
                    udtSearchResult.MH,
                    udtSearchResult.Peptide,
                    udtSearchResult.Protein,
                    udtSearchResult.NTT,
                    udtSearchResult.DeNovoScore,
                    udtSearchResult.MSGFScore,
                    udtSearchResult.SpecEValue,
                    udtSearchResult.RankSpecProb.ToString(),
                    udtSearchResult.EValue
                };

                if (blnIncludeFDRandPepFDR)
                {
                    lstData.Add(udtSearchResult.QValue);
                    lstData.Add(udtSearchResult.PepQValue);
                }
                else if (blnIncludeEFDR)
                {
                    lstData.Add(udtSearchResult.QValue);
                    lstData.Add("1");
                }

                if (blnMSGFPlus)
                {
                    lstData.Add(udtSearchResult.IsotopeError.ToString());
                }

                if (blnIncludeIMSFields)
                {
                    lstData.Add(udtSearchResult.IMSScan.ToString());
                    lstData.Add(udtSearchResult.IMSDriftTime);
                }

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

        private class MSGFDBSearchResultsComparerScanChargeSpecEValuePeptide : IComparer<udtMSGFDBSearchResultType>
        {
            public int Compare(udtMSGFDBSearchResultType x, udtMSGFDBSearchResultType y)
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

                // Charge is the same; check SpecEValue
                if (x.SpecEValueNum > y.SpecEValueNum)
                {
                    return 1;
                }

                if (x.SpecEValueNum < y.SpecEValueNum)
                {
                    return -1;
                }

                // SpecEValue is the same; check peptide
                var result = string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }
                return result;
            }
        }

        private class MSGFDBSearchResultsComparerSpecEValueScanChargePeptide : IComparer<udtMSGFDBSearchResultType>
        {
            public int Compare(udtMSGFDBSearchResultType x, udtMSGFDBSearchResultType y)
            {
                if (x.SpecEValueNum > y.SpecEValueNum)
                {
                    return 1;
                }

                if (x.SpecEValueNum < y.SpecEValueNum)
                {
                    return -1;
                }

                // SpecEValueNum is the same; check scan number
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
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }

                return result;
            }
        }

        #endregion
    }
}
