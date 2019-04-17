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
            mFileDate = "April 4, 2019";
            mModMassRegEx = new Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS);

            mPeptideCleavageStateCalculator = new clsPeptideCleavageStateCalculator();
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin);

            mNumericModErrors = 0;
        }

        #region "Constants and Enums"

        public const string TOOL_NAME = "MSGFPlus";

        public const string FILENAME_SUFFIX_MSGFDB_FILE = "_msgfdb";
        public const string FILENAME_SUFFIX_MSGFPLUS_FILE = "_msgfplus";

        public const string N_TERMINUS_SYMBOL_MSGFDB = "_.";

        public const string C_TERMINUS_SYMBOL_MSGFDB = "._";

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
        // ReSharper disable CommentTypo
        // +57.021HWWTLTTDRINK         matches +57.021
        // -57.021+42.011HWWTLTTDRINK  matches -57.021+42.011 (two separate mods)
        // +42.011MDHTPQSQLK           matches +42.011
        // ReSharper restore CommentTypo
        private const string MSGFDB_N_TERMINAL_MOD_MASS_REGEX = @"^([0-9\.\+\-]+)";

        private const string MSGFDB_MOD_MASS_REGEX = @"([+-][0-9\.]+)";

        private const string PROTEIN_AND_TERM_SYMBOLS_REGEX = @"([^;]+)\(pre=(.),post=(.)\)";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        /// <summary>
        /// These columns correspond to TSV file created by MzidToTsvConverter.exe from the MSGF+ .mzid file
        /// </summary>
        private enum eMSGFPlusResultsFileColumns
        {
            SpectrumFile = 0,
            SpecIndex = 1,               // SpecID in MSGF+
            Scan = 2,
            FragMethod = 3,
            PrecursorMZ = 4,
            PMErrorDa = 5,               // Corresponds to PrecursorError(Da)
            PMErrorPPM = 6,              // Corresponds to PrecursorError(ppm)
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

        private enum eFilteredOutputFileTypeConstants
        {
            SynFile = 0,
            FHTFile = 1
        }

        #endregion

        #region "Structures"
        private struct udtMSGFPlusSearchResultType
        {
            // ReSharper disable once NotAccessedField.Local
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
                string units;
                double equivalenceThreshold;

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
        /// For each mod symbol, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {

            var chMostRecentLetter = '-';
            var residueLocInPeptide = 0;

            var sequence = searchResult.PeptideSequenceWithMods;

            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var chChar = sequence[index];

                if (IsLetterAtoZ(chChar))
                {
                    chMostRecentLetter = chChar;
                    residueLocInPeptide += 1;

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                            if (modificationDefinition.TargetResiduesContain(chChar))
                            {
                                // Match found; add this modification
                                searchResult.SearchResultAddModification(
                                    modificationDefinition, chChar, residueLocInPeptide,
                                    searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(chMostRecentLetter))
                {
                    var success = searchResult.SearchResultAddDynamicModification(chChar, chMostRecentLetter, residueLocInPeptide, searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
                    if (!success)
                    {
                        var errorMessage = searchResult.ErrorMessage;
                        if (string.IsNullOrEmpty(errorMessage))
                        {
                            errorMessage = "SearchResultAddDynamicModification returned false for symbol " + chChar;
                        }
                        SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
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
        /// <param name="peptide"></param>
        /// <param name="kvProteinInfo"></param>
        /// <returns>Peptide sequence with N-terminal and C-Terminal residues</returns>
        /// <remarks></remarks>
        private string AddUpdatePrefixAndSuffixResidues(string peptide, KeyValuePair<string, udtTerminusCharsType> kvProteinInfo)
        {
            if (peptide.IndexOf('.') < 0)
            {
                return kvProteinInfo.Value.NTerm + "." + peptide + "." + kvProteinInfo.Value.CTerm;
            }

            string peptideNew;

            if (peptide.Length >= 2)
            {
                if (peptide[1] == '.')
                {
                    // Peptide already has the N-terminal residue
                    // Replace it using kvProteinInfo
                    peptideNew = kvProteinInfo.Value.NTerm + "." + peptide.Substring(2);
                }
                else if (peptide[0] == '.')
                {
                    peptideNew = kvProteinInfo.Value.NTerm + peptide;
                }
                else
                {
                    peptideNew = kvProteinInfo.Value.NTerm + "." + peptide;
                }
            }
            else
            {
                peptideNew = string.Copy(peptide);
            }

            if (peptideNew.Length >= 4)
            {
                if (peptideNew[peptideNew.Length - 2] == '.')
                {
                    // Peptide already has the C-terminal residue
                    // Replace it using kvProteinInfo
                    peptideNew = peptideNew.Substring(0, peptideNew.Length - 2) + "." + kvProteinInfo.Value.CTerm;
                }
                else if (peptideNew[peptideNew.Length - 1] == '.')
                {
                    peptideNew = peptideNew + kvProteinInfo.Value.CTerm;
                }
                else
                {
                    peptideNew = peptideNew + "." + kvProteinInfo.Value.CTerm;
                }
            }

            return peptideNew;
        }

        private void AppendToScanGroupDetails(
            ICollection<udtScanGroupInfoType> scanGroupDetails,
            IDictionary<string, bool> scanGroupCombo,
            udtScanGroupInfoType udtScanGroupInfo,
            ref int currentScanGroupID,
            ref int nextScanGroupID)
        {
            var chargeScanComboText = udtScanGroupInfo.Charge + "_" + udtScanGroupInfo.Scan;

            if (!scanGroupCombo.ContainsKey(chargeScanComboText))
            {
                if (currentScanGroupID < 0)
                {
                    currentScanGroupID = nextScanGroupID;
                    nextScanGroupID += 1;
                }

                udtScanGroupInfo.ScanGroupID = currentScanGroupID;

                scanGroupDetails.Add(udtScanGroupInfo);
                scanGroupCombo.Add(chargeScanComboText, true);
            }
        }

        private void AppendToSearchResults(
            ICollection<udtMSGFPlusSearchResultType> searchResults,
            udtMSGFPlusSearchResultType udtSearchResult,
            Dictionary<string, udtTerminusCharsType> proteinInfo)
        {
            if (proteinInfo.Count == 0)
            {
                searchResults.Add(udtSearchResult);
            }
            else
            {
                foreach (var kvEntry in proteinInfo)
                {
                    udtSearchResult.Protein = kvEntry.Key;
                    udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, kvEntry);

                    searchResults.Add(udtSearchResult);
                }
            }
        }

        /// <summary>
        /// Ranks each entry (assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            IList<udtMSGFPlusSearchResultType> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assigned per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            if (startIndex == endIndex)
            {
                // Only one result
                var currentResult = searchResults[startIndex];
                currentResult.RankSpecProb = 1;
                searchResults[startIndex] = currentResult;
                return;
            }

            // Duplicate a portion of searchResults so that we can sort by ascending Spectral Probability

            var resultsSubset = new Dictionary<int, udtMSGFPlusSearchResultType>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsBySpecProb = (from item in resultsSubset orderby item.Value.SpecEValueNum select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsBySpecProb)
            {
                var currentResult = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = currentResult.SpecEValueNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(currentResult.SpecEValueNum - lastValue) > double.Epsilon)
                    {
                        lastValue = currentResult.SpecEValueNum;
                        currentRank += 1;
                    }
                }

                currentResult.RankSpecProb = currentRank;

                // Because this is a list of structs, we have to copy currentResult back into the current position in searchResults
                searchResults[entry.Key] = currentResult;
            }
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            try
            {
                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts);

                // Make sure .PeptideSequenceWithMods does not have any generic mod masses
                // It should only have mod symbols
                var reMatch = mModMassRegEx.Match(searchResult.PeptideSequenceWithMods);
                if (reMatch.Success)
                {
                    // Modification mass did not have a symbol associated with it in the _ModDefs.txt file
                    // We could try to handle this, listing the modification mass in place of the modification symbol in the _ModDetails.txt file, but will
                    // instead abort processing

                    mNumericModErrors += 1;

                    if (mNumericModErrors < 250)
                    {
                        var localErrorMessage = "Search result contains a numeric mod mass that could not be associated with a modification symbol; ResultID = " + searchResult.ResultID + ", ModMass = " + reMatch.Value;
                        SetErrorMessage(localErrorMessage);
                    }
                    else if (mNumericModErrors == 250)
                    {
                        SetErrorMessage("Too many numeric mod mass results have been found; suppressing further logging");
                    }

                    return false;
                }

                // Parse .PeptideSequenceWithMods to determine the modified residues present
                AddDynamicAndStaticResidueMods(searchResult, updateModOccurrenceCounts);

                // Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
                // Since Inspect allows a terminal peptide residue to be modified twice, we'll allow that to happen,
                //  even though, biologically, that's typically not possible
                // However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus
                //  (where two COOH groups are present)
                searchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, updateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();

                return true;
            }
            catch (Exception)
            {
                return false;
            }

        }

        private short ComputeCleavageState(string sequenceWithMods)
        {
            // Remove any non-letter characters before calling .ComputeCleavageState()
            var cleanSequence = GetCleanSequence(sequenceWithMods, out var prefix, out var suffix);

            var eCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(cleanSequence, prefix, suffix);

            return Convert.ToInt16(eCleavageState);
        }

        /// <summary>
        /// This function should only be called when column PMError(Da) is present (and PMError(ppm) is not present)
        /// </summary>
        /// <param name="precursorErrorDa">Mass error (Observed - theoretical)</param>
        /// <param name="precursorMZ"></param>
        /// <param name="charge"></param>
        /// <param name="peptideMonoisotopicMass"></param>
        /// <param name="adjustPrecursorMassForC13"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private double ComputeDelMCorrectedPPM(
            double precursorErrorDa,
            double precursorMZ,
            int charge,
            double peptideMonoisotopicMass,
            bool adjustPrecursorMassForC13)
        {
            // Compute the original value for the precursor monoisotopic mass
            var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, charge, 0);

            var peptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, adjustPrecursorMassForC13, peptideMonoisotopicMass);

            return peptideDeltaMassCorrectedPpm;
        }

        /// <summary>
        /// Compute the monoisotopic mass of the peptide
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="totalModMass"></param>
        /// <returns></returns>
        private double ComputePeptideMass(string peptide, double totalModMass)
        {
            var cleanSequence = GetCleanSequence(peptide);

            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);
            mass += totalModMass;

            return mass;
        }

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var suffixesToFind = new List<string> {
                "_msgfplus_syn",
                "_msgfplus_fht",
                "_msgfdb_syn",
                "_msgfdb_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// Parses the digits in modDigits to convert them to one or more modification symbols
        /// </summary>
        /// <param name="currentResidue"></param>
        /// <param name="modDigits">Example: +57.021 or +79.9663+14.0157 or -18.0106</param>
        /// <param name="modSymbols"></param>
        /// <param name="dynModSymbols"></param>
        /// <param name="msgfPlusModInfo"></param>
        /// <param name="nTerminalMod"></param>
        /// <param name="possibleCTerminalMod"></param>
        /// <param name="modMassFound"></param>
        /// <param name="containsStaticMod"></param>
        /// <returns>True if success; false if a problem</returns>
        /// <remarks></remarks>
        private bool ConvertMSGFModMassesToSymbols(
            string currentResidue,
            string modDigits,
            out string modSymbols,
            out string dynModSymbols,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> msgfPlusModInfo,
            bool nTerminalMod,
            bool possibleCTerminalMod,
            out double modMassFound,
            out bool containsStaticMod)
        {
            double bestMassDiff = 0;
            var modSymbolsFound = 0;
            var symbolBestMatch = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
            var residuesBestBatch = string.Empty;

            modSymbols = string.Empty;
            dynModSymbols = string.Empty;
            modMassFound = 0;
            containsStaticMod = false;

            var reMatches = mModMassRegEx.Matches(modDigits);

            foreach (Match reMatch in reMatches)
            {
                var modMassText = reMatch.Value;

                // Convert modMass to a mass value
                var modMass = double.Parse(modMassText);
                modMassFound += modMass;

                var matchFound = false;
                int bestMatchIndex;

                while (true)
                {
                    bestMatchIndex = -1;

                    // Step through the known modifications to find the closest match
                    for (var index = 0; index <= msgfPlusModInfo.Count - 1; index++)
                    {
                        var testMod = true;

                        if (nTerminalMod)
                        {
                            // Only test N-terminal mods
                            if (!(msgfPlusModInfo[index].ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.DynNTermPeptide ||
                                  msgfPlusModInfo[index].ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.DynNTermProtein))
                            {
                                testMod = false;
                            }
                        }
                        else if (!possibleCTerminalMod)
                        {
                            // Skip C-terminal mods since we're not at the C-terminus
                            if (msgfPlusModInfo[index].ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.DynCTermPeptide ||
                                msgfPlusModInfo[index].ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.DynCTermProtein)
                            {
                                testMod = false;
                            }
                        }

                        if (!testMod)
                            continue;

                        var candidateMassDiff = Math.Abs(msgfPlusModInfo[index].ModMassVal - modMass);
                        if (!(candidateMassDiff < 0.25))
                            continue;

                        // Possible match found
                        var updateCandidate = false;

                        if (bestMatchIndex < 0 ||
                            candidateMassDiff < bestMassDiff)
                        {
                            updateCandidate = true;
                        }
                        else if (Math.Abs(candidateMassDiff - bestMassDiff) < float.Epsilon &&
                                 symbolBestMatch == '-' &&
                                 msgfPlusModInfo[index].ModSymbol != '-')
                        {
                            // Masses are the same, but the existing candidate is a static mod

                            // If this new one is a dynamic mod, switch to it, but only if residuesBestBatch does not contain the residue while
                            // the residues for the new candidate does contain the residue

                            if (!residuesBestBatch.Contains(currentResidue) &&
                                msgfPlusModInfo[index].Residues.Contains(currentResidue))
                            {
                                updateCandidate = true;
                            }
                        }

                        if (updateCandidate)
                        {
                            bestMatchIndex = index;
                            bestMassDiff = candidateMassDiff;
                            symbolBestMatch = msgfPlusModInfo[index].ModSymbol;
                            residuesBestBatch = string.Copy(msgfPlusModInfo[index].Residues);
                            matchFound = true;
                        }
                    }

                    if (!matchFound && nTerminalMod)
                    {
                        // Set this to false, then search again
                        nTerminalMod = false;
                    }
                    else if (!matchFound && !possibleCTerminalMod)
                    {
                        // Set this to true, then search again
                        possibleCTerminalMod = true;
                    }
                    else
                    {
                        // matchFound is true
                        break;
                    }
                }

                if (matchFound)
                {
                    // Match found; use the mod symbol
                    modSymbols += msgfPlusModInfo[bestMatchIndex].ModSymbol;
                    modSymbolsFound += 1;

                    if (msgfPlusModInfo[bestMatchIndex].ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.StaticMod)
                    {
                        containsStaticMod = true;
                    }
                    else
                    {
                        dynModSymbols += msgfPlusModInfo[bestMatchIndex].ModSymbol;
                    }
                }
                else
                {
                    // Match not found; use the mass value
                    modSymbols += modMass;
                }
            }

            return modSymbolsFound > 0;
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MSGF+
        /// The synopsis file includes every result with a p-value below a set threshold or a SpecEValue below a certain threshold
        /// The first-hits file includes the results with the lowest SpecEValue (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="scanGroupFilePath"></param>
        /// <param name="msgfPlusModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
        /// <param name="isMsgfPlus">Output parameter: this function will set this to True if we're processing MSGF+ results</param>
        /// <param name="specIdToIndex"></param>
        /// <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool CreateFHTorSYNResultsFile(
            string inputFilePath,
            string outputFilePath,
            string scanGroupFilePath,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> msgfPlusModInfo,
            out bool isMsgfPlus,
            IDictionary<string, int> specIdToIndex,
            eFilteredOutputFileTypeConstants eFilteredOutputFileType)
        {

            var searchResultsCurrentScan = new List<udtMSGFPlusSearchResultType>();
            var searchResultsPrefiltered = new List<udtMSGFPlusSearchResultType>();

            isMsgfPlus = false;

            try
            {
                var scanGroupDetails = new List<udtScanGroupInfoType>();
                var scanGroupCombo = new Dictionary<string, bool>();

                mPrecursorMassErrorWarningCount = 0;

                // Look for custom amino acids
                var customAA = (from item in msgfPlusModInfo
                                where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.CustomAA
                                select item).ToList();

                foreach (var customAADef in customAA)
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
                    string errorLog;

                    using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        errorLog = string.Empty;
                        var headerParsed = false;
                        var includeFDRandPepFDR = false;
                        var includeEFDR = false;
                        var includeIMSFields = false;

                        var nextScanGroupID = 1;
                        scanGroupDetails.Clear();
                        scanGroupCombo.Clear();

                        // Initialize the array that will hold all of the records that will ultimately be written out to disk
                        var filteredSearchResults = new List<udtMSGFPlusSearchResultType>();

                        // Initialize a dictionary that tracks the peptide sequence for each combo of scan and charge
                        // Keys are Scan_Charge, values track the clean sequence, the associated protein name, and the protein number for that name
                        // Note that we can only track protein numbers if the FASTA file path was provided at the command line
                        var scanChargeFirstHit = new Dictionary<string, clsFirstHitInfo>();

                        var columnMapping = new Dictionary<eMSGFPlusResultsFileColumns, int>();

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
                                var validHeader = ParseMSGFPlusResultsFileHeaderLine(lineIn, columnMapping);
                                if (!validHeader)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }

                                headerParsed = true;

                                if (columnMapping[eMSGFPlusResultsFileColumns.FDR_QValue] >= 0 ||
                                    columnMapping[eMSGFPlusResultsFileColumns.PepFDR_PepQValue] >= 0)
                                {
                                    includeFDRandPepFDR = true;
                                }
                                else if (columnMapping[eMSGFPlusResultsFileColumns.EFDR] >= 0)
                                {
                                    includeEFDR = true;
                                }

                                if (columnMapping[eMSGFPlusResultsFileColumns.IMSDriftTime] >= 0)
                                {
                                    includeIMSFields = true;
                                }

                                if (columnMapping[eMSGFPlusResultsFileColumns.IsotopeError] >= 0)
                                {
                                    isMsgfPlus = true;
                                }

                                // Write the header line
                                WriteSynFHTFileHeader(writer, ref errorLog, includeFDRandPepFDR, includeEFDR, includeIMSFields, isMsgfPlus);

                                continue;
                            }

                            var validSearchResult = ParseMSGFPlusResultsFileEntry(lineIn, isMsgfPlus, msgfPlusModInfo,
                                                                                  searchResultsCurrentScan, ref errorLog,
                                                                                  columnMapping, ref nextScanGroupID, scanGroupDetails,
                                                                                  scanGroupCombo, specIdToIndex);

                            if (!validSearchResult || searchResultsCurrentScan.Count <= 0)
                            {
                                continue;
                            }

                            if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                            {
                                // Synopsis file
                                validSearchResult = MSGFPlusResultPassesSynFilter(searchResultsCurrentScan[0]);
                            }
                            else
                            {
                                // First Hits file
                                var scanChargeKey = searchResultsCurrentScan[0].Scan + "_" + searchResultsCurrentScan[0].Charge;

                                if (scanChargeFirstHit.TryGetValue(scanChargeKey, out var firstHitPeptide))
                                {
                                    // A result has already been stored for this scan/charge combo
                                    validSearchResult = false;

                                    // Possibly update the associated protein name
                                    if (firstHitPeptide.CleanSequence.Equals(
                                        GetCleanSequence(searchResultsCurrentScan[0].Peptide, out var newPrefix, out var newSuffix)))
                                    {
                                        var bestProtein = GetBestProteinName(firstHitPeptide.ProteinName, firstHitPeptide.ProteinNumber,
                                                                             searchResultsCurrentScan[0].Protein);
                                        if (bestProtein.Value < firstHitPeptide.ProteinNumber)
                                        {
                                            firstHitPeptide.ProteinName = bestProtein.Key;
                                            firstHitPeptide.ProteinNumber = bestProtein.Value;
                                            firstHitPeptide.UpdatePrefixAndSuffix(newPrefix, newSuffix);
                                        }
                                    }
                                }
                                else
                                {
                                    firstHitPeptide = new clsFirstHitInfo(searchResultsCurrentScan[0].Peptide,
                                                                          GetCleanSequence(searchResultsCurrentScan[0].Peptide))
                                    {
                                        ProteinName = searchResultsCurrentScan[0].Protein,
                                        ProteinNumber = int.MaxValue
                                    };

                                    if (mProteinNameOrder.TryGetValue(searchResultsCurrentScan[0].Protein, out var proteinNumber))
                                    {
                                        firstHitPeptide.ProteinNumber = proteinNumber;
                                    }

                                    scanChargeFirstHit.Add(scanChargeKey, firstHitPeptide);
                                }
                            }

                            if (validSearchResult)
                            {
                                ExpandListIfRequired(searchResultsPrefiltered, searchResultsCurrentScan.Count);
                                searchResultsPrefiltered.AddRange(searchResultsCurrentScan);
                            }

                            // Update the progress
                            var percentComplete = Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100);
                            if (CreateProteinModsFile)
                            {
                                percentComplete = percentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                            }

                            UpdateProgress(percentComplete);
                        }

                        searchResultsPrefiltered.TrimExcess();

                        // Sort the SearchResults by scan, charge, and ascending SpecEValue
                        searchResultsPrefiltered.Sort(new MSGFDBSearchResultsComparerScanChargeSpecEValuePeptide());

                        if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.FHTFile)
                        {
                            // Update the protein names in searchResultsPrefiltered using scanChargeFirstHit
                            // This step is typically not necessary, but is often required for SplitFasta results

                            for (var index = 0; index <= searchResultsPrefiltered.Count - 1; index++)
                            {
                                var scanChargeKey = searchResultsPrefiltered[index].Scan + "_" + searchResultsPrefiltered[index].Charge;

                                if (!scanChargeFirstHit.TryGetValue(scanChargeKey, out var firstHitPeptide))
                                    continue;

                                if (searchResultsPrefiltered[index].Protein.Equals(firstHitPeptide.ProteinName))
                                    continue;

                                // Protein name doesn't match the expected name
                                if (firstHitPeptide.CleanSequence.Equals(GetCleanSequence(searchResultsPrefiltered[index].Peptide)))
                                {
                                    // Update the protein name
                                    var updatedSearchResult = searchResultsPrefiltered[index];
                                    updatedSearchResult.Peptide = firstHitPeptide.SequenceWithModsAndContext;
                                    updatedSearchResult.Protein = string.Copy(firstHitPeptide.ProteinName);
                                    searchResultsPrefiltered[index] = updatedSearchResult;
                                }
                                else
                                {
                                    Console.WriteLine("Possible programming bug; " +
                                                      "mix of peptides tracked for a given scan/charge combo when caching data for First Hits files; " +
                                                      $"see scan_charge {scanChargeKey}");
                                }
                            }
                        }

                        // Now filter the data and store in filteredSearchResults
                        // Due to code updates in October 2016, searchResultsPrefiltered already has filtered data
                        var startIndex = 0;

                        while (startIndex < searchResultsPrefiltered.Count)
                        {
                            var endIndex = startIndex;
                            // Find all of the peptides with the same scan number
                            while (endIndex + 1 < searchResultsPrefiltered.Count &&
                                   searchResultsPrefiltered[endIndex + 1].ScanNum == searchResultsPrefiltered[startIndex].ScanNum)
                            {
                                endIndex += 1;
                            }

                            // Store the results for this scan
                            if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                            {
                                StoreSynMatches(searchResultsPrefiltered, startIndex, endIndex, filteredSearchResults);
                            }
                            else
                            {
                                StoreTopFHTMatch(searchResultsPrefiltered, startIndex, endIndex, filteredSearchResults);
                            }

                            startIndex = endIndex + 1;
                        }

                        // Sort the data in udtFilteredSearchResults then write out to disk
                        SortAndWriteFilteredSearchResults(writer, filteredSearchResults, ref errorLog,
                                                          includeFDRandPepFDR, includeEFDR, includeIMSFields, isMsgfPlus);
                    }

                    // Write out the scan group info
                    if (!string.IsNullOrEmpty(scanGroupFilePath))
                    {
                        StoreScanGroupInfo(scanGroupFilePath, scanGroupDetails);
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
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    return false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }

        }

        /// <summary>
        /// Extracts mod info from either a MSGF+ param file or from a MSGFPlus_Mods.txt file (previously MSGFDB_Mods.txt)
        /// </summary>
        /// <param name="msgfPlusParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns>True if success; false if a problem</returns>
        /// <remarks></remarks>
        private bool ExtractModInfoFromParamFile(
            string msgfPlusParamFilePath,
            out List<clsMSGFPlusParamFileModExtractor.udtModInfoType> modInfo)
        {
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor(SEARCH_ENGINE_NAME);

            RegisterEvents(modFileProcessor);
            modFileProcessor.ErrorEvent += ModExtractorErrorHandler;

            var success = modFileProcessor.ExtractModInfoFromParamFile(
                msgfPlusParamFilePath,
                clsMSGFPlusParamFileModExtractor.ModSpecFormats.MSGFPlusAndMSPathFinder,
                out modInfo);

            if (!success || mErrorCode != ePHRPErrorCodes.NoError)
            {
                if (mErrorCode == ePHRPErrorCodes.NoError)
                {
                    SetErrorMessage("Unknown error extracting the modification definitions from the MSGF+ parameter file");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                }
                return false;
            }

            modFileProcessor.ResolveMSGFPlusModsWithModDefinitions(modInfo, mPeptideMods);

            return true;
        }

        /// <summary>
        /// Extracts parent mass tolerance from the parameters loaded from an MSGF+ parameter file
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <returns>Parent mass tolerance info.  Tolerances will be 0 if an error occurs</returns>
        /// <remarks></remarks>
        private udtParentMassToleranceType ExtractParentMassToleranceFromParamFile(clsSearchEngineParameters searchEngineParams)
        {
            const string PM_TOLERANCE_TAG = "PMTolerance";

            var udtParentMassToleranceInfo = new udtParentMassToleranceType();

            try
            {
                udtParentMassToleranceInfo.Clear();

                if (searchEngineParams.Parameters.TryGetValue(PM_TOLERANCE_TAG, out var value))
                {
                    // Parent ion tolerance line found

                    // Split the line on commas
                    var splitLine = value.Split(',');

                    double tolerance;
                    bool isPPM;

                    if (splitLine.Length == 1)
                    {
                        if (ParseParentMassTolerance(splitLine[0], out tolerance, out isPPM))
                        {
                            udtParentMassToleranceInfo.ToleranceLeft = tolerance;
                            udtParentMassToleranceInfo.ToleranceRight = tolerance;
                            udtParentMassToleranceInfo.IsPPM = isPPM;
                        }
                    }
                    else if (splitLine.Length > 1)
                    {
                        if (ParseParentMassTolerance(splitLine[0], out tolerance, out isPPM))
                        {
                            udtParentMassToleranceInfo.ToleranceLeft = tolerance;
                            udtParentMassToleranceInfo.IsPPM = isPPM;

                            if (ParseParentMassTolerance(splitLine[1], out tolerance, out isPPM))
                            {
                                udtParentMassToleranceInfo.ToleranceRight = tolerance;
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
                    Path.GetFileName(searchEngineParams.SearchEngineParamFilePath), ex.Message), ex);
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
        /// <param name="pepToProteinMapFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="msgfPlusModInfo"></param>
        /// <param name="isMsgfPlus">Should be set to True if processing MSGF+ results</param>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadPeptideToProteinMapInfoMSGFDB(
            string pepToProteinMapFilePath,
            string outputDirectoryPath,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> msgfPlusModInfo,
            bool isMsgfPlus,
            List<udtPepToProteinMappingType> pepToProteinMapping,
            out string mtsPepToProteinMapFilePath)
        {


            mtsPepToProteinMapFilePath = string.Empty;

            try
            {
                if (string.IsNullOrWhiteSpace(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file is not defined");
                    return false;
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    var pepToProteinMapAlternate = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(pepToProteinMapFilePath, "Dataset_msgfdb.txt");
                    if (File.Exists(pepToProteinMapAlternate))
                    {
                        pepToProteinMapFilePath = pepToProteinMapAlternate;
                    }
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file does not exist: " + pepToProteinMapFilePath);
                    return false;
                }

                // Initialize pepToProteinMapping
                pepToProteinMapping.Clear();

                // Read the data in proteinToPeptideMappingFilePath
                var success = LoadPeptideToProteinMapInfo(pepToProteinMapFilePath, pepToProteinMapping, out var headerLine);

                if (!success)
                {
                    return false;
                }

                mtsPepToProteinMapFilePath =
                    Path.Combine(outputDirectoryPath, Path.GetFileNameWithoutExtension(pepToProteinMapFilePath) + "MTS.txt");

                using (var writer = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read))
                )
                {
                    if (!string.IsNullOrEmpty(headerLine))
                    {
                        // Header line
                        writer.WriteLine(headerLine);
                    }

                    for (var index = 0; index <= pepToProteinMapping.Count - 1; index++)
                    {
                        // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                        // In addition, replace the * terminus symbols with dashes
                        var mtsCompatiblePeptide = ReplaceMSGFModTextWithSymbol(ReplaceTerminus(pepToProteinMapping[index].Peptide),
                                                                               msgfPlusModInfo, isMsgfPlus, out _);

                        if (pepToProteinMapping[index].Peptide != mtsCompatiblePeptide)
                        {
                            UpdatePepToProteinMapPeptide(pepToProteinMapping, index, mtsCompatiblePeptide);
                        }

                        writer.WriteLine(pepToProteinMapping[index].Peptide + "\t" +
                                         pepToProteinMapping[index].Protein + "\t" +
                                         pepToProteinMapping[index].ResidueStart + "\t" +
                                         pepToProteinMapping[index].ResidueEnd);
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" + Path.GetFileName(mtsPepToProteinMapFilePath) + "): " + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }

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

            var searchEngineParams = new clsSearchEngineParameters(SEARCH_ENGINE_NAME);

            var success = clsPHRPParser.ReadKeyValuePairSearchEngineParamFile(SEARCH_ENGINE_NAME, msgfPlusParamFilePath, clsPHRPReader.ePeptideHitResultType.MSGFPlus,
                                                                              searchEngineParams, out var localErrorMessage, out var localWarningMessage);

            if (!string.IsNullOrWhiteSpace(localErrorMessage))
            {
                ReportError(localErrorMessage);
                return false;
            }

            if (!string.IsNullOrWhiteSpace(localWarningMessage))
            {
                ReportWarning(localWarningMessage);
            }

            if (searchEngineParams.Parameters.Count == 0)
            {
                SetErrorMessage("MSGF+ parameter file is empty; unable to extract parent mass tolerance info");
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }

            // Parse the PMTolerance setting
            mParentMassToleranceInfo = ExtractParentMassToleranceFromParamFile(searchEngineParams);

            // Parse the ChargeCarrierMass setting
            if (clsPHRPParserMSGFDB.GetCustomChargeCarrierMass(searchEngineParams, out var customChargeCarrierMass))
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

        private bool MSGFPlusResultPassesSynFilter(udtMSGFPlusSearchResultType udtMSGFPlusSearchResultType)
        {
            if (udtMSGFPlusSearchResultType.EValueNum <= MSGFPlusSynopsisFileEValueThreshold ||
                udtMSGFPlusSearchResultType.SpecEValueNum <= MSGFPlusSynopsisFileSpecEValueThreshold ||
                udtMSGFPlusSearchResultType.QValueNum > 0 && udtMSGFPlusSearchResultType.QValueNum < 0.01)
            {
                return true;
            }

            return false;
        }

        private bool ParseMSGFPlusSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            List<udtPepToProteinMappingType> pepToProteinMapping,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            var successOverall = true;

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (resetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                mNumericModErrors = 0;

                // Initialize searchResult
                var searchResult = new clsSearchResultsMSGFDB(mPeptideMods, mPeptideSeqMassCalculator);

                // Note that MSGF+ synopsis files are normally sorted on SpecEValue value, ascending
                // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
                //  we will keep track of the scan, charge, and peptide information parsed for each unique SpecEValue encountered
                // Although this was a possibility with Inspect, it likely never occurs for MSGF+
                //  But, we'll keep the check in place just in case

                var peptidesFoundForSpecEValueLevel = new SortedSet<string>();

                var previousSpecEValue = string.Empty;

                // Assure that pepToProteinMapping is sorted on peptide
                if (pepToProteinMapping.Count > 1)
                {
                    pepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);

                    var errorLog = string.Empty;

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        var resultsProcessed = 0;
                        var headerParsed = false;

                        // Create the output files
                        var baseOutputFilePath = Path.Combine(outputDirectoryPath, Path.GetFileName(inputFilePath));
                        var filesInitialized = InitializeSequenceOutputFiles(baseOutputFilePath);
                        if (!filesInitialized)
                            return false;

                        var columnMapping = new Dictionary<clsPHRPParserMSGFDB.MSGFPlusSynFileColumns, int>();

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
                                var validHeader = ParseMSGFPlusSynFileHeaderLine(lineIn, columnMapping);
                                if (!validHeader)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return false;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMSGFPlusSynFileEntry(lineIn, searchResult, ref errorLog,
                                                                           resultsProcessed, columnMapping,
                                                                           out var currentPeptideWithMods);

                            resultsProcessed += 1;
                            if (!validSearchResult)
                            {
                                continue;
                            }

                            var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                            bool firstMatchForGroup;

                            if (searchResult.SpecEValue == previousSpecEValue)
                            {
                                // New result has the same SpecEValue as the previous result
                                // See if peptidesFoundForSpecEValueLevel contains the peptide, scan and charge

                                if (peptidesFoundForSpecEValueLevel.Contains(key))
                                {
                                    firstMatchForGroup = false;
                                }
                                else
                                {
                                    peptidesFoundForSpecEValueLevel.Add(key);
                                    firstMatchForGroup = true;
                                }
                            }
                            else
                            {
                                // New SpecEValue
                                peptidesFoundForSpecEValueLevel.Clear();

                                // Update previousSpecEValue
                                previousSpecEValue = searchResult.SpecEValue;

                                // Append a new entry
                                peptidesFoundForSpecEValueLevel.Add(key);
                                firstMatchForGroup = true;
                            }

                            var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                            if (!modsAdded)
                            {
                                successOverall = false;
                                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                {
                                    errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'" +
                                                "\n";
                                }
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
                                    var currentProtein = string.Copy(searchResult.ProteinName);
                                    do
                                    {
                                        if (pepToProteinMapping[pepToProteinMapIndex].Protein != currentProtein)
                                        {
                                            searchResult.ProteinName = string.Copy(pepToProteinMapping[pepToProteinMapIndex].Protein);
                                            SaveResultsFileEntrySeqInfo(searchResult, false);
                                        }

                                        pepToProteinMapIndex += 1;
                                    } while (pepToProteinMapIndex < pepToProteinMapping.Count &&
                                             currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                                }
                                else
                                {
                                    // Match not found; this is unexpected
                                    ReportWarning("no match for '" + currentPeptideWithMods + "' in pepToProteinMapping");
                                }
                            }

                            // Update the progress
                            var percentComplete = Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100);
                            if (CreateProteinModsFile)
                            {
                                percentComplete = percentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                            }
                            UpdateProgress(percentComplete);
                        }
                    }

                    if (CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(inputFilePath);
                        var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                        SaveModificationSummaryFile(modificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (errorLog.Length > 0)
                    {
                        SetErrorMessage("Invalid Lines: " + "\n" + errorLog);
                    }

                    return successOverall;
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
                SetErrorMessage(ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                return false;
            }

        }

        private bool ParseMSGFPlusResultsFileEntry(
            string lineIn,
            bool isMsgfPlus,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> msgfPlusModInfo,
            ICollection<udtMSGFPlusSearchResultType> searchResultsCurrentScan,
            ref string errorLog,
            IDictionary<eMSGFPlusResultsFileColumns, int> columnMapping,
            ref int nextScanGroupID,
            ICollection<udtScanGroupInfoType> scanGroupDetails,
            IDictionary<string, bool> scanGroupCombo,
            IDictionary<string, int> specIdToIndex)
        {
            // Parses an entry from the MSGF+ results file

            var udtSearchResult = new udtMSGFPlusSearchResultType();
            string rowIndex = null;

            udtMSGFPlusSearchResultType[] udtMergedScanInfo = null;

            try
            {

                var proteinInfo = new Dictionary<string, udtTerminusCharsType>();

                // Reset searchResults
                searchResultsCurrentScan.Clear();

                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 13)
                {
                    // Not a valid result
                    return false;
                }

                rowIndex = splitLine[0];

                if (!GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.SpectrumFile], out udtSearchResult.SpectrumFileName))
                {
                    ReportError("SpectrumFile column is missing or invalid", true);
                }
                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.SpecIndex], out udtSearchResult.SpecIndex);

                if (isMsgfPlus)
                {
                    var generateSpecIndex = true;

                    if (!int.TryParse(udtSearchResult.SpecIndex, out var specIndex))
                    {
                        // MSGF+ includes text in the SpecID column, for example: "controllerType=0 controllerNumber=1 scan=6390" or "index=4323"
                        // Need to convert these to an integer

                        if (udtSearchResult.SpecIndex.StartsWith("index="))
                        {
                            udtSearchResult.SpecIndex = udtSearchResult.SpecIndex.Substring("index=".Length);
                            if (int.TryParse(udtSearchResult.SpecIndex, out specIndex))
                            {
                                generateSpecIndex = false;
                            }
                        }

                        if (generateSpecIndex)
                        {
                            if (!specIdToIndex.TryGetValue(udtSearchResult.SpecIndex, out specIndex))
                            {
                                specIndex = specIdToIndex.Count + 1;
                                specIdToIndex.Add(udtSearchResult.SpecIndex, specIndex);
                            }

                            udtSearchResult.SpecIndex = specIndex.ToString();
                        }
                    }
                }

                if (!GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.Scan], out udtSearchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.FragMethod], out udtSearchResult.FragMethod);

                var slashIndex = udtSearchResult.Scan.IndexOf('/');
                int scanCount;
                if (slashIndex > 0)
                {
                    // This is a merged spectrum and thus scan number looks like: 3010/3011/3012
                    // Split the Scan list on the slash
                    // Later in this function, we'll append searchResults with this scan plus the other scans

                    var splitResult = udtSearchResult.Scan.Split('/');
                    scanCount = splitResult.Length;
                    udtMergedScanInfo = new udtMSGFPlusSearchResultType[scanCount];

                    for (var index = 0; index <= scanCount - 1; index++)
                    {
                        udtMergedScanInfo[index] = new udtMSGFPlusSearchResultType();
                        udtMergedScanInfo[index].Clear();
                        udtMergedScanInfo[index].Scan = splitResult[index];
                        udtMergedScanInfo[index].ScanNum = CIntSafe(splitResult[index], 0);
                    }

                    // Now split SpecIndex and store in udtMergedScanInfo
                    splitResult = udtSearchResult.SpecIndex.Split('/');

                    for (var index = 0; index <= splitResult.Length - 1; index++)
                    {
                        if (index >= udtMergedScanInfo.Length)
                        {
                            // There are more entries for SpecIndex than there are for Scan#; this is unexpected
                            break;
                        }
                        udtMergedScanInfo[index].SpecIndex = splitResult[index];
                    }

                    // Now split FragMethod and store in udtMergedScanInfo
                    splitResult = udtSearchResult.FragMethod.Split('/');

                    for (var index = 0; index <= splitResult.Length - 1; index++)
                    {
                        if (index >= udtMergedScanInfo.Length)
                        {
                            // There are more entries for FragMethod than there are for Scan#; this is unexpected
                            break;
                        }
                        udtMergedScanInfo[index].FragMethod = splitResult[index];
                    }
                }
                else
                {
                    udtSearchResult.ScanNum = CIntSafe(udtSearchResult.Scan, 0);
                    scanCount = 1;
                }

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.PrecursorMZ], out udtSearchResult.PrecursorMZ);

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Precursor mass error could be in PPM or Da
                //   In MSGFDB, the header line will have PMError(ppm)        or PMError(Da)
                //   In MSGF+,  the header line will have PrecursorError(ppm) or PrecursorError(Da)
                double precursorErrorDa = 0;

                if (columnMapping[eMSGFPlusResultsFileColumns.PMErrorPPM] >= 0)
                {
                    GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.PMErrorPPM], out udtSearchResult.PMErrorPPM);
                }
                else
                {
                    GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.PMErrorDa], out udtSearchResult.PMErrorDa);
                    precursorErrorDa = CDblSafe(udtSearchResult.PMErrorDa, 0);
                    udtSearchResult.PMErrorPPM = string.Empty; // We'll populate this column later in this function
                }

                if (!GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                {
                    ReportError("Peptide column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.Protein], out udtSearchResult.Protein);

                // MSGF+ .tsv files may have a semicolon separated list of protein names; check for this
                udtSearchResult.Protein = SplitProteinList(udtSearchResult.Protein, proteinInfo);

                if (proteinInfo.Count > 0)
                {
                    // Need to add the prefix and suffix residues
                    udtSearchResult.Peptide = AddUpdatePrefixAndSuffixResidues(udtSearchResult.Peptide, proteinInfo.First());
                }

                // Replace any mod text values in the peptide sequence with the appropriate mod symbols
                // In addition, replace the terminus symbols with dashes
                udtSearchResult.Peptide =
                    ReplaceMSGFModTextWithSymbol(ReplaceTerminus(udtSearchResult.Peptide), msgfPlusModInfo, isMsgfPlus,
                                                 out var totalModMass);

                // Compute monoisotopic mass of the peptide
                var peptideMonoisotopicMass = ComputePeptideMass(udtSearchResult.Peptide, totalModMass);

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoisotopicMass, 0), 6);

                if (!string.IsNullOrEmpty(udtSearchResult.PMErrorPPM))
                {
                    // Convert the ppm-based PM Error to Da-based

                    if (double.TryParse(udtSearchResult.PrecursorMZ, out var precursorMZ))
                    {
                        // Note that since .PMErrorPPM is present, the Precursor m/z is a C13-corrected m/z value
                        // In other words, it may not be the actual m/z selected for fragmentation.

                        if (double.TryParse(udtSearchResult.PMErrorPPM, out var pMErrorPPM))
                        {
                            if (mParentMassToleranceInfo.IsPPM &&
                                (pMErrorPPM < -mParentMassToleranceInfo.ToleranceLeft * 1.5 ||
                                 pMErrorPPM > mParentMassToleranceInfo.ToleranceRight * 1.5))
                            {
                                // PPM error computed by MSGF+ is more than 1.5-fold larger than the ppm-based parent ion tolerance; don't trust the value computed by MSGF+

                                mPrecursorMassErrorWarningCount += 1;
                                if (mPrecursorMassErrorWarningCount <= 10)
                                {
                                    ReportWarning("Precursor mass error computed by MSGF+ is 1.5-fold larger than search tolerance: " +
                                                  udtSearchResult.PMErrorPPM + " vs. " + mParentMassToleranceInfo.ToleranceLeft.ToString("0") +
                                                  "ppm," + mParentMassToleranceInfo.ToleranceRight.ToString("0") + "ppm");
                                    if (mPrecursorMassErrorWarningCount == 10)
                                    {
                                        ReportWarning("Additional mass errors will not be reported");
                                    }
                                }

                                var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, udtSearchResult.ChargeNum, 0);

                                precursorErrorDa = precursorMonoMass - peptideMonoisotopicMass;

                                udtSearchResult.PMErrorPPM = string.Empty;
                            }
                            else
                            {
                                precursorErrorDa = clsPeptideMassCalculator.PPMToMass(pMErrorPPM, peptideMonoisotopicMass);

                                // Note that this will be a C13-corrected precursor error; not the true precursor error
                                udtSearchResult.PMErrorDa = MassErrorToString(precursorErrorDa);

                            }
                        }
                    }
                }

                if (string.IsNullOrEmpty(udtSearchResult.PMErrorPPM))
                {
                    if (double.TryParse(udtSearchResult.PrecursorMZ, out var precursorMZ))
                    {
                        var peptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(precursorErrorDa, precursorMZ,
                                                                                      udtSearchResult.ChargeNum, peptideMonoisotopicMass,
                                                                                      true);

                        udtSearchResult.PMErrorPPM = PRISM.StringUtilities.DblToString(peptideDeltaMassCorrectedPpm, 5, 0.00005);

                        if (string.IsNullOrEmpty(udtSearchResult.PMErrorDa))
                        {
                            precursorErrorDa = clsPeptideMassCalculator.PPMToMass(peptideDeltaMassCorrectedPpm, peptideMonoisotopicMass);

                            // Note that this will be a C13-corrected precursor error; not the true precursor error
                            udtSearchResult.PMErrorDa = MassErrorToString(precursorErrorDa);
                        }
                    }
                }

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.DeNovoScore], out udtSearchResult.DeNovoScore);
                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.MSGFScore], out udtSearchResult.MSGFScore);
                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.SpecProb_EValue], out udtSearchResult.SpecEValue);
                if (!double.TryParse(udtSearchResult.SpecEValue, out udtSearchResult.SpecEValueNum))
                    udtSearchResult.SpecEValueNum = 0;

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.PValue_EValue], out udtSearchResult.EValue);
                if (!double.TryParse(udtSearchResult.EValue, out udtSearchResult.EValueNum))
                    udtSearchResult.EValueNum = 0;

                var targetDecoyFDRValid = GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.FDR_QValue], out udtSearchResult.QValue);
                if (!double.TryParse(udtSearchResult.QValue, out udtSearchResult.QValueNum))
                    udtSearchResult.QValueNum = 0;

                if (targetDecoyFDRValid)
                {
                    GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.PepFDR_PepQValue], out udtSearchResult.PepQValue);
                }
                else
                {
                    GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.EFDR], out udtSearchResult.QValue);
                }

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.IsotopeError], out udtSearchResult.IsotopeError);

                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.IMSScan], out udtSearchResult.IMSScan);
                GetColumnValue(splitLine, columnMapping[eMSGFPlusResultsFileColumns.IMSDriftTime], out udtSearchResult.IMSDriftTime);

                udtSearchResult.NTT = ComputeCleavageState(udtSearchResult.Peptide).ToString();

                var udtScanGroupInfo = new udtScanGroupInfoType();
                var currentScanGroupID = -1;

                udtScanGroupInfo.Charge = udtSearchResult.ChargeNum;

                if (scanCount > 1)
                {
                    // Append one entry to searchResults for each item in udtMergedScanInfo()

                    for (var index = 0; index <= scanCount - 1; index++)
                    {
                        udtSearchResult.Scan = udtMergedScanInfo[index].Scan;
                        udtSearchResult.ScanNum = udtMergedScanInfo[index].ScanNum;

                        udtSearchResult.SpecIndex = udtMergedScanInfo[index].SpecIndex;
                        udtSearchResult.FragMethod = udtMergedScanInfo[index].FragMethod;

                        AppendToSearchResults(searchResultsCurrentScan, udtSearchResult, proteinInfo);

                        // Append an entry to scanGroupDetails
                        udtScanGroupInfo.Scan = udtSearchResult.ScanNum;
                        AppendToScanGroupDetails(scanGroupDetails, scanGroupCombo, udtScanGroupInfo, ref currentScanGroupID,
                                                 ref nextScanGroupID);
                    }
                }
                else
                {
                    // This is not a merged result; simply append udtSearchResult to searchResults
                    AppendToSearchResults(searchResultsCurrentScan, udtSearchResult, proteinInfo);

                    // Also append an entry to scanGroupDetails
                    udtScanGroupInfo.Scan = udtSearchResult.ScanNum;
                    AppendToScanGroupDetails(scanGroupDetails, scanGroupCombo, udtScanGroupInfo, ref currentScanGroupID,
                                             ref nextScanGroupID);
                }

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the MSGF+ results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        errorLog += "Error parsing MSGF+ Results in ParseMSGFPlusResultsFileEntry for RowIndex '" + rowIndex + "'" +
                                       "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MSGF+ Results in ParseMSGFPlusResultsFileEntry" + "\n";
                    }
                }
                return false;
            }

        }

        /// <summary>
        /// Parse the MSGF+ results file header line
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns></returns>
        private bool ParseMSGFPlusResultsFileHeaderLine(string lineIn, IDictionary<eMSGFPlusResultsFileColumns, int> columnMapping)
        {

            // The expected header from MSGFDB is:
            // #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(Da)           Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR
            // or
            // #SpecFile    SpecIndex    Scan#     FragMethod    Precursor                    PMError(ppm)          Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecProb      P-value   FDR       PepFDR

            // The expected header from MSGF+ is:
            // #SpecFile    SpecID       ScanNum   FragMethod    Precursor    IsotopeError    PrecursorError(Da)    Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue
            // or
            // #SpecFile    SpecID       ScanNum   FragMethod    Precursor    IsotopeError    PrecursorError(ppm)   Charge    Peptide    Protein    DeNovoScore    MSGFScore    SpecEValue    EValue    QValue    PepQValue

            var columnNames = new SortedDictionary<string, eMSGFPlusResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"#SpecFile", eMSGFPlusResultsFileColumns.SpectrumFile},
                {"SpecIndex", eMSGFPlusResultsFileColumns.SpecIndex},
                {"SpecID", eMSGFPlusResultsFileColumns.SpecIndex},
                {"Scan#", eMSGFPlusResultsFileColumns.Scan},
                {"ScanNum", eMSGFPlusResultsFileColumns.Scan},
                {"FragMethod", eMSGFPlusResultsFileColumns.FragMethod},
                {"Precursor", eMSGFPlusResultsFileColumns.PrecursorMZ},
                {"IsotopeError", eMSGFPlusResultsFileColumns.IsotopeError},
                {"PMError(Da)", eMSGFPlusResultsFileColumns.PMErrorDa},
                {"PrecursorError(Da)", eMSGFPlusResultsFileColumns.PMErrorDa},
                {"PMError(ppm)", eMSGFPlusResultsFileColumns.PMErrorPPM},
                {"PrecursorError(ppm)", eMSGFPlusResultsFileColumns.PMErrorPPM},
                {"Charge", eMSGFPlusResultsFileColumns.Charge},
                {"Peptide", eMSGFPlusResultsFileColumns.Peptide},
                {"Protein", eMSGFPlusResultsFileColumns.Protein},
                {"DeNovoScore", eMSGFPlusResultsFileColumns.DeNovoScore},
                {"MSGFScore", eMSGFPlusResultsFileColumns.MSGFScore},
                {"SpecProb", eMSGFPlusResultsFileColumns.SpecProb_EValue},
                {"SpecEValue", eMSGFPlusResultsFileColumns.SpecProb_EValue},
                {"P-value", eMSGFPlusResultsFileColumns.PValue_EValue},
                {"EValue", eMSGFPlusResultsFileColumns.PValue_EValue},
                {"FDR", eMSGFPlusResultsFileColumns.FDR_QValue},
                {"QValue", eMSGFPlusResultsFileColumns.FDR_QValue},
                {"PepFDR", eMSGFPlusResultsFileColumns.PepFDR_PepQValue},
                {"PepQValue", eMSGFPlusResultsFileColumns.PepFDR_PepQValue},
                {"EFDR", eMSGFPlusResultsFileColumns.EFDR},
                {clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Scan, eMSGFPlusResultsFileColumns.IMSScan},
                {clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Drift_Time, eMSGFPlusResultsFileColumns.IMSDriftTime}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (eMSGFPlusResultsFileColumns resultColumn in Enum.GetValues(typeof(eMSGFPlusResultsFileColumns)))
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
                    else
                    {
                        // Unrecognized column name
                        Console.WriteLine("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMSGFPlusResultsFileHeaderLine");
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing header in MSGFPlus results file: " + ex.Message, ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse the header line of a MSGF+ _syn.txt file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns></returns>
        private bool ParseMSGFPlusSynFileHeaderLine(string lineIn, IDictionary<clsPHRPParserMSGFDB.MSGFPlusSynFileColumns, int> columnMapping)
        {
            var columnNames = clsPHRPParserMSGFDB.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (clsPHRPParserMSGFDB.MSGFPlusSynFileColumns resultColumn in Enum.GetValues(typeof(clsPHRPParserMSGFDB.MSGFPlusSynFileColumns)))
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
                SetErrorMessage("Error parsing header in MSGFPlus synopsis file: " + ex.Message, ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse an entry from a MSGF+ Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorLog"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequenceWithMods"></param>
        /// <returns></returns>
        private bool ParseMSGFPlusSynFileEntry(
            string lineIn,
            clsSearchResultsMSGFDB searchResult,
            ref string errorLog,
            int resultsProcessed,
            IDictionary<clsPHRPParserMSGFDB.MSGFPlusSynFileColumns, int> columnMapping,
            out string peptideSequenceWithMods)
        {

            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();
            peptideSequenceWithMods = string.Empty;

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 15)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MSGF+ Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from MSGF+ Results line " +
                                    (resultsProcessed + 1) + "\n";
                    }

                    return false;
                }

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.DelM], out string msgfPlusComputedDelM);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.DelMPPM], out string msgfPlusComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.MSGFPlusComputedDelM = msgfPlusComputedDelM;
                searchResult.MSGFPlusComputedDelMPPM = msgfPlusComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.MSGFPlusComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MSGF+ results file as "Observed_Mass - Theoretical_Mass"
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

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.FragMethod], out string fragMethod);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.MH], out string peptideMh);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.NTT], out string ntt);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.DeNovoScore], out string deNovoScore);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.MSGFScore], out string msgfScore);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.SpecProb_EValue], out string specEValue);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.RankSpecProb], out string rankSpecEValue);
                GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.PValue_EValue], out string eValue);

                var targetDecoyFDRValid = GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.FDR_QValue], out string qValue);

                searchResult.FragMethod = fragMethod;
                searchResult.PrecursorMZ = precursorMz;
                searchResult.PeptideMH = peptideMh;
                searchResult.NTT = ntt;
                searchResult.DeNovoScore = deNovoScore;
                searchResult.MSGFScore = msgfScore;
                searchResult.SpecEValue = specEValue;
                searchResult.RankSpecEValue = rankSpecEValue;
                searchResult.EValue = eValue;
                searchResult.QValue = qValue;

                if (targetDecoyFDRValid)
                {
                    GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.PepFDR_PepQValue], out string pepQValue);
                    searchResult.PepQValue = pepQValue;
                }
                else
                {
                    GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.EFDR], out string efdr);
                    searchResult.QValue = efdr;
                }

                if (columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.IsotopeError] >= 0)
                {
                    GetColumnValue(splitLine, columnMapping[clsPHRPParserMSGFDB.MSGFPlusSynFileColumns.IsotopeError], out string isotopeError);
                    searchResult.IsotopeError = isotopeError;
                    searchResult.MSGFPlusResults = true;
                }
                else
                {
                    searchResult.MSGFPlusResults = false;
                }

                // Compute PrecursorMH using PrecursorMZ and charge
                if (double.TryParse(searchResult.PrecursorMZ, out var precursorMZ))
                {
                    if (int.TryParse(searchResult.Charge, out var chargeValue))
                    {
                        searchResult.ParentIonMH =
                            PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, chargeValue), 6);
                    }
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
                        errorLog += "Error parsing MSGFDB Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MSGFDB Results in ParseMSGFPlusSynFileEntry" + "\n";
                    }
                }
                return false;
            }

        }

        private bool ParseParentMassTolerance(string toleranceText, out double tolerance, out bool isPPM)
        {
            tolerance = 0;

            toleranceText = toleranceText.ToLower().Trim();

            if (toleranceText.EndsWith("da"))
            {
                toleranceText = toleranceText.Substring(0, toleranceText.Length - 2);
                isPPM = false;
            }
            else if (toleranceText.EndsWith("ppm"))
            {
                toleranceText = toleranceText.Substring(0, toleranceText.Length - 3);
                isPPM = true;
            }
            else
            {
                isPPM = false;
                return false;
            }

            if (double.TryParse(toleranceText, out tolerance))
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MSGFDB results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file for data processing</param>
        /// <returns>True if success, False if failure</returns>
        /// <remarks>Use SearchToolParameterFilePath to define the search engine parameter file</remarks>
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

                success = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!success)
                {
                    return false;
                }

                mPeptideSeqMassCalculator.ResetAminoAcidMasses();

                var specIdToIndex = new Dictionary<string, int>();

                ResetProgress("Parsing " + Path.GetFileName(inputFilePath));

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);

                    var pepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the MSGF+ Parameter File so that we can determine the modification names and masses
                    // If the MSGFPlus_Mods.txt or MSGFDB_Mods.txt file was defined, the mod symbols in that file will be used to define the mod symbols in msgfPlusModInfo
                    var modInfoExtracted = ExtractModInfoFromParamFile(SearchToolParameterFilePath, out var msgfPlusModInfo);
                    if (!modInfoExtracted)
                    {
                        return false;
                    }

                    if (!LoadSearchEngineParamFile(SearchToolParameterFilePath))
                    {
                        return false;
                    }

                    var query = from item in msgfPlusModInfo where item.ModType == clsMSGFPlusParamFileModExtractor.eMSGFPlusModType.CustomAA select item;
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

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_msgfdb" with "_msgfplus"
                    if (baseName.EndsWith("_msgfdb", StringComparison.OrdinalIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - "_msgfdb".Length) + "_msgfplus";
                    }

                    string fhtOutputFilePath;

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

                        fhtOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_FIRST_HITS_FILE_SUFFIX);

                        var scanGroupFilePath = string.Empty;

                        success = CreateFHTorSYNResultsFile(
                            inputFilePath, fhtOutputFilePath, scanGroupFilePath, msgfPlusModInfo,
                            out _, specIdToIndex, eFilteredOutputFileTypeConstants.FHTFile);
                    }
                    else
                    {
                        fhtOutputFilePath = string.Empty;
                    }

                    if (CreateInspectSynopsisFile)
                    {
                        // Create the synopsis output file
                        ResetProgress("Creating the SYN file", true);

                        // The synopsis file name will be of the form BasePath_msgfplus_syn.txt
                        var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                        var scanGroupFilePath = Path.Combine(outputDirectoryPath, baseName + "_ScanGroupInfo.txt");

                        success = CreateFHTorSYNResultsFile(
                            inputFilePath, synOutputFilePath, scanGroupFilePath, msgfPlusModInfo,
                            out var isMsgfPlus, specIdToIndex, eFilteredOutputFileTypeConstants.SynFile);

                        // Load the PeptideToProteinMap information; if the file doesn't exist, a warning will be displayed, but processing will continue
                        // LoadPeptideToProteinMapInfoMSGFDB also creates _msgfplus_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
                        var pepToProteinMapFilePath = ConstructPepToProteinMapFilePath(Path.Combine(outputDirectoryPath, baseName) + ".txt", outputDirectoryPath, mts: false);

                        ResetProgress("Loading the PepToProtein map file: " + Path.GetFileName(pepToProteinMapFilePath), true);

                        LoadPeptideToProteinMapInfoMSGFDB(pepToProteinMapFilePath, outputDirectoryPath, msgfPlusModInfo, isMsgfPlus, pepToProteinMapping, out var mtsPepToProteinMapFilePath);

                        // Create the other PHRP-specific files
                        ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                        // Now parse the _syn.txt file that we just created to create the other PHRP files
                        success = ParseMSGFPlusSynopsisFile(synOutputFilePath, outputDirectoryPath, pepToProteinMapping, false);

                        // Remove all items from pepToProteinMapping to reduce memory overhead
                        pepToProteinMapping.Clear();
                        pepToProteinMapping.TrimExcess();

                        if (success && CreateProteinModsFile)
                        {
                            success = CreateProteinModsFileWork(baseName, inputFile, fhtOutputFilePath, synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath);
                        }
                    }

                    if (success)
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

            return success;
        }

        private bool CreateProteinModsFileWork(
            string baseName,
            FileInfo inputFile,
            string fhtOutputFilePath,
            string synOutputFilePath,
            string outputDirectoryPath,
            string mtsPepToProteinMapFilePath)
        {
            var success = true;

            if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath) || !File.Exists(mtsPepToProteinMapFilePath))
            {
                // MTSPepToProteinMap file not found; auto-create it

                if (string.IsNullOrEmpty(mtsPepToProteinMapFilePath))
                {
                    mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseName, outputDirectoryPath, mts: true);
                }

                var sourcePHRPDataFiles = new List<string>();

                if (!string.IsNullOrEmpty(fhtOutputFilePath))
                {
                    sourcePHRPDataFiles.Add(fhtOutputFilePath);
                }

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
                        // Auto-change mIgnorePeptideToProteinMapperErrors to True
                        // We only do this for MSGFDB since it often includes reverse protein peptides in the results even though the FASTA file often does not have reverse proteins
                        IgnorePeptideToProteinMapperErrors = true;
                        success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                        if (!success)
                        {
                            ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                        }
                    }
                }
            }

            if (success)
            {
                if (inputFile.Directory == null)
                {
                    ReportWarning("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                }
                else if (string.IsNullOrWhiteSpace(synOutputFilePath))
                {
                    ReportWarning("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputDirectoryPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          clsPHRPReader.ePeptideHitResultType.MSGFPlus);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                success = true;
            }

            return true;
        }

        private static readonly Regex NTerminalModMassMatcher = new Regex(MSGFDB_N_TERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);
        private static readonly Regex ModMassMatcher = new Regex(MSGFDB_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Replaces modification masses in peptide sequences with modification symbols (uses case-sensitive comparisons)
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="msgfPlusModInfo">This function assumes that each entry in msgfPlusModInfo has both .ModName and .ModSymbol defined</param>
        /// <param name="isMsgfPlus">Should be set to True if processing MSGF+ results</param>
        /// <param name="totalModMass">Output parameter: total mass of all modifications</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ReplaceMSGFModTextWithSymbol(
            string peptide,
            IReadOnlyList<clsMSGFPlusParamFileModExtractor.udtModInfoType> msgfPlusModInfo,
            bool isMsgfPlus,
            out double totalModMass)
        {
            var prefix = string.Empty;
            var suffix = string.Empty;

            var possibleCTerminalMod = false;
            bool containsStaticMod;

            // Reset the total mod mass
            totalModMass = 0;

            // Remove the prefix and suffix residues
            if (peptide.Length >= 4)
            {
                if (peptide[1] == '.' &&
                    peptide[peptide.Length - 2] == '.')
                {
                    prefix = peptide.Substring(0, 2);
                    suffix = peptide.Substring(peptide.Length - 2, 2);

                    peptide = peptide.Substring(2, peptide.Length - 4);
                }
            }

            // peptide should now be the primary peptide sequence, without the prefix or suffix residues

            // First look for dynamic N-terminal mods (NTermPeptide or NTermProtein)
            // This RegEx will match one or more mods, all at the N-terminus
            var reMatch = NTerminalModMassMatcher.Match(peptide);

            if (reMatch.Success)
            {
                // Convert the mod mass (or masses) to one or more mod symbols

                if (ConvertMSGFModMassesToSymbols("-", reMatch.Groups[1].Value,
                    out var modSymbols, out var dynModSymbols, msgfPlusModInfo,
                    true, false,
                    out var modMassFound, out containsStaticMod))
                {
                    // Replace the mod digits with the mod symbols

                    peptide = ReplaceMSGFModTextWithMatchedSymbol(peptide, reMatch.Groups[1], modSymbols, dynModSymbols, isMsgfPlus, containsStaticMod);
                    totalModMass += modMassFound;
                }
            }

            // Next, step through the peptide and parse each mod mass that follows a residue
            // Any mod mass at the end must be considered a C-terminal mod

            // Need to start at the first letter
            // If we had N-terminal mods, they're currently notated like this: _.+42.011MDHTPQSQLK.L or _.+42.011+57.021MNDR.Q
            // We want things to look like this: -.#MDHTPQSQLK.L or -.#*MNDRQLNHR.S

            // In MSGFDB, static mods do not have a mod mass listed
            // In MSGF+,  static mods do have a mod mass listed
            // Regardless, we do not add mod symbols for static mods, but we do increment totalModMass

            // Find the index of the last residue
            var index = peptide.Length - 1;
            while (index > 0 && !IsLetterAtoZ(peptide[index]))
            {
                index -= 1;
            }
            var indexLastResidue = index;

            // Find the index of the first residue
            index = 0;
            while (index < peptide.Length && !IsLetterAtoZ(peptide[index]))
            {
                index += 1;
            }
            var indexFirstResidue = index;

            var currentResidue = "-";

            while (index < peptide.Length)
            {
                if (IsLetterAtoZ(peptide[index]))
                {
                    currentResidue = peptide[index].ToString();

                    if (!isMsgfPlus)
                    {
                        // Look for static mods that should be applied to this residue (only applies to MSGFDB, not MSGF+)
                        for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                        {
                            var eModificationType = mPeptideMods.GetModificationTypeByIndex(modIndex);

                            clsModificationDefinition modificationDefinition;
                            if (eModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                            {
                                modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                                if (modificationDefinition.TargetResiduesContain(peptide[index]))
                                {
                                    // Match found; update totalModMass but do not add a static mod symbol
                                    totalModMass += modificationDefinition.ModificationMass;
                                }
                            }
                            else if (index == indexFirstResidue)
                            {
                                if (eModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod && prefix == "_")
                                {
                                    // N-terminal protein static mod
                                    modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                                    totalModMass += modificationDefinition.ModificationMass;
                                }
                                else if (eModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod)
                                {
                                    // N-terminal peptide static mod
                                    modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                                    totalModMass += modificationDefinition.ModificationMass;
                                }
                            }
                        }
                    }

                    index += 1;

                    if (index == indexLastResidue)
                        possibleCTerminalMod = true;
                }
                else
                {
                    // Found a mod; find the extent of the mod digits
                    reMatch = ModMassMatcher.Match(peptide, index);

                    // Note that possibleCTerminalMod will be set to True once we hit the last residue

                    // Convert the mod mass (or masses) to one or more mod symbols

                    if (ConvertMSGFModMassesToSymbols(currentResidue, reMatch.Groups[1].Value,
                        out var modSymbols, out var dynModSymbols, msgfPlusModInfo,
                        false, possibleCTerminalMod, out var modMassFound, out containsStaticMod))
                    {
                        peptide = ReplaceMSGFModTextWithMatchedSymbol(peptide, reMatch.Groups[1], modSymbols, dynModSymbols, isMsgfPlus, containsStaticMod);
                        totalModMass += modMassFound;

                        if (isMsgfPlus && containsStaticMod)
                        {
                            // MSGF+ shows mod masses for static mods
                            // Thus, we have removed the static mod mass and did not add a mod symbol
                            // Therefore, leave index unchanged
                        }
                        else
                        {
                            index += modSymbols.Length;
                        }
                    }
                    else
                    {
                        var addOn = reMatch.Groups[1].Value.Length;
                        if (addOn == 0)
                            index += 1;
                        else
                            index += addOn;
                    }
                }
            }

            // If any N-terminal mods were present, we need to move them to after the first residue
            // in other words, switch from #MDHTPQSQLK to M#DHTPQSQLK
            //                          or #*MNDRQLNHR to M#*NDRQLNHR

            // Update indexFirstResidue
            indexFirstResidue = 0;
            while (indexFirstResidue < peptide.Length && !IsLetterAtoZ(peptide[indexFirstResidue]))
            {
                indexFirstResidue += 1;
            }

            if (indexFirstResidue > 0 && indexFirstResidue < peptide.Length)
            {
                var peptideNew = peptide[indexFirstResidue] + peptide.Substring(0, indexFirstResidue);
                if (indexFirstResidue < peptide.Length - 1)
                {
                    peptideNew += peptide.Substring(indexFirstResidue + 1);
                }
                peptide = string.Copy(peptideNew);
            }

            return prefix + peptide + suffix;
        }

        private string ReplaceMSGFModTextWithMatchedSymbol(
            string peptide,
            Capture reGroup,
            string modSymbols,
            string dynModSymbols,
            bool isMsgfPlus,
            bool containsStaticMod)
        {
            string peptideNew;

            if (reGroup.Index > 0)
            {
                peptideNew = peptide.Substring(0, reGroup.Index);
            }
            else
            {
                peptideNew = string.Empty;
            }

            if (isMsgfPlus && containsStaticMod)
            {
                // MSGF+ shows mod masses for static mods
                // However, for consistency with other PHRP results, we do not add a symbol to the peptide for this static mod
                // Catch: If we have a peptide/terminus affected by both a static and a dynamic mod, we still want the dynamic mod.
                if (!string.IsNullOrWhiteSpace(dynModSymbols))
                {
                    peptideNew += dynModSymbols;
                }
            }
            else
            {
                peptideNew += modSymbols;
            }

            if (reGroup.Index + reGroup.Length < peptide.Length)
            {
                peptideNew += peptide.Substring(reGroup.Index + reGroup.Length);
            }

            return peptideNew;
        }

        private string ReplaceTerminus(string peptide)
        {
            if (peptide.StartsWith(N_TERMINUS_SYMBOL_MSGFDB))
            {
                peptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_MSGFDB.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_MSGFDB))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_MSGFDB.Length) + "." + clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        private static readonly Regex ProteinInfoMatcher = new Regex(PROTEIN_AND_TERM_SYMBOLS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Examines proteinList to look for a semi-colon separated list of proteins and terminus symbols, for example
        /// AT1G26570.1(pre=K,post=N);AT3G29360.1(pre=K,post=N);AT3G29360.2(pre=K,post=N)
        /// </summary>
        /// <param name="proteinList">Protein list to examine</param>
        /// <param name="proteinInfo">Protein information, if it is of the form ProteinName(pre=X,post=Y)</param>
        /// <returns>The name of the first protein</returns>
        /// <remarks></remarks>
        private string SplitProteinList(string proteinList, IDictionary<string, udtTerminusCharsType> proteinInfo)
        {
            proteinInfo.Clear();

            var reMatches = ProteinInfoMatcher.Matches(proteinList);

            if (reMatches.Count == 0)
            {
                // No match; likely just one protein
                return TruncateProteinName(proteinList);
            }

            foreach (Match reMatch in reMatches)
            {
                var proteinName = TruncateProteinName(reMatch.Groups[1].Value);

                if (proteinInfo.ContainsKey(proteinName))
                {
                    // Skip this protein since it's already present
                }
                else
                {
                    var terminusChars = new udtTerminusCharsType
                    {
                        NTerm = reMatch.Groups[2].Value[0],
                        CTerm = reMatch.Groups[3].Value[0]
                    };

                    proteinInfo.Add(proteinName, terminusChars);
                }
            }

            return proteinInfo.First().Key;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<udtMSGFPlusSearchResultType> filteredSearchResults,
            ref string errorLog,
            bool includeFDRandPepFDR,
            bool includeEFDR,
            bool includeIMSFields,
            bool isMsgfPlus)
        {
            // Sort filteredSearchResults by ascending SpecEValue, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.SpecEValueNum, item.ScanNum, item.ChargeNum, item.Peptide, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, ref errorLog, includeFDRandPepFDR, includeEFDR, includeIMSFields, isMsgfPlus);
                index += 1;
            }
        }

        private void StoreScanGroupInfo(string scanGroupFilePath, IReadOnlyCollection<udtScanGroupInfoType> scanGroupDetails)
        {
            try
            {
                // Only create the ScanGroup file if one or more scan groups exist
                // Step through scanGroupDetails to check for this
                var scanGroupIDPrevious = -1;
                var createFile = false;
                foreach (var udtScanGroupInfo in scanGroupDetails)
                {
                    if (udtScanGroupInfo.ScanGroupID == scanGroupIDPrevious)
                    {
                        createFile = true;
                        break;
                    }
                    scanGroupIDPrevious = udtScanGroupInfo.ScanGroupID;
                }

                if (createFile)
                {
                    using (var writer = new StreamWriter(new FileStream(scanGroupFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        writer.WriteLine("Scan_Group_ID" + "\t" + "Charge" + "\t" + "Scan");

                        foreach (var udtScanGroupInfo in scanGroupDetails)
                        {
                            writer.WriteLine(udtScanGroupInfo.ScanGroupID + "\t" + udtScanGroupInfo.Charge + "\t" + udtScanGroupInfo.Scan);
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
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Filtered search results</param>
        /// <remarks></remarks>
        private void StoreTopFHTMatch(
            IList<udtMSGFPlusSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            List<udtMSGFPlusSearchResultType> filteredSearchResults)
        {
            AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex);

            // The calling procedure should have already sorted by scan, charge, and SpecEValue; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            // Now store the first match for each charge for this scan
            // When storing, we use the protein name that occurred first in the FASTA file

            var udtCurrentResult = searchResults[startIndex];
            var currentCharge = udtCurrentResult.ChargeNum;
            var currentProteinNumber = int.MaxValue;
            var currentPeptide = GetCleanSequence(udtCurrentResult.Peptide);

            for (var index = startIndex; index <= endIndex; index++)
            {
                if (currentCharge != searchResults[index].ChargeNum)
                {
                    // New charge state
                    // Store udtCurrentResult (from the previous charge state)
                    filteredSearchResults.Add(udtCurrentResult);

                    udtCurrentResult = searchResults[index];
                    currentCharge = udtCurrentResult.ChargeNum;
                    currentProteinNumber = int.MaxValue;
                    currentPeptide = GetCleanSequence(udtCurrentResult.Peptide);
                }

                var newPeptide = GetCleanSequence(searchResults[index].Peptide);
                if (currentPeptide.Equals(newPeptide))
                {
                    var bestProtein = GetBestProteinName(udtCurrentResult.Protein, currentProteinNumber, searchResults[index].Protein);
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
            filteredSearchResults.Add(udtCurrentResult);
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Filtered search results</param>
        /// <remarks></remarks>
        private void StoreSynMatches(
            IList<udtMSGFPlusSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            List<udtMSGFPlusSearchResultType> filteredSearchResults)
        {
            AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            var results = new SortedSet<string>();

            // Now store or write out the matches that pass the filters
            // By default, filter passing peptides have MSGFDB_SpecEValue <= 5E-7 Or EValue less than 0.75 or QValue less than 1% (but not 0)
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (MSGFPlusResultPassesSynFilter(searchResults[index]))
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

                    var resultKey = searchResults[index].Peptide + "_" +
                                    searchResults[index].Protein + "_" +
                                    searchResults[index].MH + "_" +
                                    searchResults[index].SpecEValue;

                    if (results.Contains(resultKey))
                    {
                        continue;
                    }

                    results.Add(resultKey);
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorLog"></param>
        /// <param name="includeFDRandPepFDR"></param>
        /// <param name="includeEFDR"></param>
        /// <param name="includeIMSFields"></param>
        /// <param name="isMsgfPlus"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ref string errorLog,
            bool includeFDRandPepFDR,
            bool includeEFDR,
            bool includeIMSFields,
            bool isMsgfPlus)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var knownHeaderColumns = clsPHRPParserMSGFDB.GetColumnHeaderNamesAndIDs();

                var data = new List<string>
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

                if (isMsgfPlus)
                {
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFPlus_SpecEValue);
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue);
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EValue);
                }
                else
                {
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_MSGFDB_SpecProb);
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFDB_SpecProb);
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PValue);
                }

                if (includeFDRandPepFDR)
                {
                    if (isMsgfPlus)
                    {
                        data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_QValue);
                        data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepQValue);
                    }
                    else
                    {
                        data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_FDR);
                        data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR);
                    }
                }
                else if (includeEFDR)
                {
                    // Note that we'll write out a "1" for "PepFDR" for every result
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_EFDR);
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_PepFDR);
                }

                if (isMsgfPlus)
                {
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_Isotope_Error);
                }

                if (includeIMSFields)
                {
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Scan);
                    data.Add(clsPHRPParserMSGFDB.DATA_COLUMN_IMS_Drift_Time);
                }

                foreach (var headerName in data)
                {
                    if (!knownHeaderColumns.ContainsKey(headerName))
                    {
                        errorLog += "Unrecognized header name for the synopsis / first hits file: " + headerName + "\n";
                    }
                }

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
        /// <param name="includeFDRandPepFDR"></param>
        /// <param name="includeEFDR"></param>
        /// <param name="includeIMSFields"></param>
        /// <param name="isMsgfPlus"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtMSGFPlusSearchResultType udtSearchResult,
            ref string errorLog,
            bool includeFDRandPepFDR,
            bool includeEFDR,
            bool includeIMSFields,
            bool isMsgfPlus)
        {
            try
            {
                // Primary Columns (other columns are added in certain circumstances):
                //
                // MSGFDB
                // ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecProb    Rank_MSGFDB_SpecProb    PValue  FDR     PepFDR

                // MSGF+
                // ResultID  Scan FragMethod  SpecIndex  Charge  PrecursorMZ  DelM  DelM_PPM  MH  Peptide  Protein  NTT  DeNovoScore  MSGFScore  MSGFDB_SpecEValue  Rank_MSGFDB_SpecEValue  EValue  QValue  PepQValue  IsotopeError

                var data = new List<string>
                {
                    resultID.ToString(),
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

                if (includeFDRandPepFDR)
                {
                    data.Add(udtSearchResult.QValue);
                    data.Add(udtSearchResult.PepQValue);
                }
                else if (includeEFDR)
                {
                    data.Add(udtSearchResult.QValue);
                    data.Add("1");
                }

                if (isMsgfPlus)
                {
                    data.Add(udtSearchResult.IsotopeError.ToString());
                }

                if (includeIMSFields)
                {
                    data.Add(udtSearchResult.IMSScan.ToString());
                    data.Add(udtSearchResult.IMSDriftTime);
                }

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

        #region "IComparer Classes"

        private class MSGFDBSearchResultsComparerScanChargeSpecEValuePeptide : IComparer<udtMSGFPlusSearchResultType>
        {
            public int Compare(udtMSGFPlusSearchResultType x, udtMSGFPlusSearchResultType y)
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

        #endregion
    }
}
