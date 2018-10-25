// This class reads in an MODa results file (txt format) and creates
// a tab-delimited text file with the data.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started 04/01/2014
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
    public class clsMODaResultsProcessor : clsPHRPBaseClass
    {
        public clsMODaResultsProcessor()
        {
            mFileDate = "October 13, 2017";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_MODA_FILE = "_moda.id";

        public const string N_TERMINUS_SYMBOL_MODA = "-";
        public const string C_TERMINUS_SYMBOL_MODA = "-";

        public const float DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD = clsMODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        // Note that as of April 2014, all mod masses reported by MODa are simply integers, meaning matching a trailing period is not necessary
        private const string MODA_MOD_MASS_REGEX = @"([+-][0-9.]+)";

        private const byte MODA_MASS_DIGITS_OF_PRECISION = 0;

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        // These columns correspond to the tab-delimited file (_moda.id.txt) created by MODa's anal_moda.jar file
        private const int MODaResultsFileColCount = 11;
        public enum eMODaResultsFileColumns
        {
            SpectrumFileName = 0,
            SpectrumIndex = 1,
            ObservedMonoMass = 2,
            Charge = 3,
            CalculatedMonoMass = 4,
            DeltaMass = 5,
            Score = 6,
            Probability = 7,
            Peptide = 8,
            Protein = 9,
            PeptidePosition = 10
        }

        // These columns correspond to the Synopsis file created by this class
        private const int MODaSynFileColCount = 15;
        public enum eMODaSynFileColumns
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
            Protein = 9,                         // Protein Name
            Score = 10,
            Probability = 11,
            Rank_Probability = 12,
            Peptide_Position = 13,
            QValue = 14
        }

        #endregion

        #region "Structures"
        // This data structure holds rows read from the tab-delimited file (_moda.id.txt) created by MODa's anal_moda.jar file
        private struct udtMODaSearchResultType
        {
            public string SpectrumFileName;
            public string SpectrumIndex;
            public int ScanNum;                     // Determined by looking for SpectrumIndex in the _mgf_IndexToScanMap.txt file
            public string Precursor_mass;           // Uncharged monoisotopic mass value of the observed precursor_mz, reported as ObservedMonoMass by MODa
            public string PrecursorMZ;              // Computed from ObservedMonoMass
            public string Charge;
            public short ChargeNum;
            public string CalculatedMonoMass;       // Theoretical monoisotopic mass of the peptide (including mods), as computed by MODa
            public string DeltaMass;                // Computed by MODa
            public string MH;                       // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
            public string DelM;                     // Computed using Precursor_mass - CalculatedMonoMass
            public string DelM_PPM;                 // Computed using DelM and CalculatedMonoMass
            public string Score;
            public string Probability;              // Higher values are better
            public double ProbabilityNum;           // Higher values are better
            public int RankProbability;
            public string Peptide;
            public string Protein;
            public string PeptidePosition;          // Protein start/stop residues of the peptide, e.g. 108~115
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
                Probability = string.Empty;
                ProbabilityNum = 0;
                RankProbability = 0;
                Peptide = string.Empty;
                Protein = string.Empty;
                PeptidePosition = string.Empty;
                FDR = 0;
                QValue = 0;
            }
        }

        #endregion

        #region "Classwide Variables"
        private Dictionary<int, int> mSpectrumIndexToScanMap;
        #endregion

        /// <summary>
        /// Step through .PeptideSequenceWithMods
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const char NO_RESIDUE = '-';

            var parsingModMass = false;
            var modMassDigits = string.Empty;

            var mostRecentResidue = NO_RESIDUE;
            var residueLocInPeptide = 0;

            var sequence = searchResult.PeptideSequenceWithMods;
            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var chChar = sequence[index];

                if (IsLetterAtoZ(chChar))
                {
                    if (parsingModMass)
                    {
                        // Associate the mod mass in modMassDigits with the previous residue
                        AssociateDynamicModWithResidue(searchResult, mostRecentResidue, residueLocInPeptide,
                            modMassDigits, updateModOccurrenceCounts);
                        parsingModMass = false;
                    }

                    mostRecentResidue = chChar;
                    residueLocInPeptide += 1;

                    // Look for static mods to associate with this residue
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
                else
                {
                    var isNumberChar = chChar == '+' || chChar == '-' || char.IsDigit(chChar);

                    if (parsingModMass)
                    {
                        if (isNumberChar || chChar == '.')
                        {
                            modMassDigits += chChar;
                        }
                    }
                    else if (isNumberChar)
                    {
                        // Mod Mass Start
                        modMassDigits = chChar.ToString();
                        parsingModMass = true;
                    }
                    else
                    {
                        // Unrecognized symbol; ignore it
                    }
                }
            }

            if (parsingModMass)
            {
                // Associate the mod mass in modMassDigits with the previous residue
                AssociateDynamicModWithResidue(searchResult, mostRecentResidue, residueLocInPeptide, modMassDigits, updateModOccurrenceCounts);
            }
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = true;

            bool success;

            try
            {
                // If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
                searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts);

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

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
        }

        //// This function was an experiment to compute better DelM_PPM values
        //// by reading the synopsis file with PHRPReader and re-computing the DelM_PPM values based on the monoisotopic mass values computed for the sequences
        //// It turned out to not be required, since the DelM_PPM values reported by MODa are quite accurate (despite the fact that it reports integer mod mass values)
        //private bool AppendDelMPPMRefinedToSynFile(string synOutputFilePath)
        //{
        //    const string SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED = "DelM_PPM_Refined";
        //
        //    bool success = false;
        //
        //    try
        //    {
        //        // Keys in this dictionary are ResultID values from the synopsis file
        //        // Values are refined DelM_PPM values
        //        var dctRefinedDelMPPMErrors = new Dictionary<int, double>();
        //
        //        using (clsPHRPReader reader = new clsPHRPReader(synOutputFilePath, clsPHRPReader.ePeptideHitResultType.MODa, loadModsAndSeqInfo: true, loadMSGFResults: false))
        //        {
        //            reader.EchoMessagesToConsole = true;
        //            reader.SkipDuplicatePSMs = true;
        //            reader.SkipDuplicatePSMs = false;
        //
        //            foreach (string errorMessage in reader.ErrorMessages)
        //            {
        //                SetErrorMessage(errorMessage);
        //            }
        //
        //            foreach (string warningMessage in reader.WarningMessages)
        //            {
        //                ReportWarning(warningMessage);
        //            }
        //
        //            reader.ClearErrors();
        //            reader.ClearWarnings();
        //
        //            reader.ErrorEvent += PHRPReader_ErrorEvent;
        //            reader.WarningEvent += PHRPReader_WarningEvent;
        //
        //            while (reader.MoveNext())
        //            {
        //                var oPSM = reader.CurrentPSM;
        //                var oSeqInfo = reader.CurrentPSMSeqInfo;
        //
        //                if (oSeqInfo != null)
        //                {
        //                    var delM = oPSM.PrecursorNeutralMass - oSeqInfo.MonoisotopicMass;
        //
        //                    var peptideDeltaMassRefinedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(delM, oPSM.PrecursorNeutralMass, true, oSeqInfo.MonoisotopicMass);
        //
        //                    double originalDelMPPM;
        //                    if (double.TryParse(oPSM.MassErrorPPM, out originalDelMPPM))
        //                    {
        //                        if (Math.Abs(peptideDeltaMassRefinedPpm - originalDelMPPM) > 2)
        //                        {
        //                            Console.WriteLine("Computed a refined DelMPPM value: " + peptideDeltaMassRefinedPpm.ToString("0.0") + " vs. " + originalDelMPPM.ToString("0.0"));
        //                        }
        //                    }
        //
        //                    dctRefinedDelMPPMErrors.Add(oPSM.ResultID, peptideDeltaMassRefinedPpm);
        //                }
        //            }
        //
        //            reader.ErrorEvent -= PHRPReader_ErrorEvent;
        //            reader.WarningEvent -= PHRPReader_WarningEvent;
        //
        //        }
        //
        //        var synOutputFilePathNew = synOutputFilePath + ".refinedDelMPPM";
        //        bool headersParsed = false;
        //        bool swapFiles = true;
        //
        //        using (var reader = new StreamReader(new FileStream(synOutputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
        //        {
        //            using (var swOutfile = new StreamWriter(new FileStream(synOutputFilePathNew, FileMode.Create, FileAccess.Write, FileShare.Read)))
        //            {
        //                while (!reader.EndOfStream)
        //                {
        //                    var lineIn = reader.ReadLine();
        //
        //                    var splitLine = lineIn.Split('\t');
        //
        //                    if (!headersParsed)
        //                    {
        //                        if (splitLine.Contains(SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED))
        //                        {
        //                            // This file already has the refined DelM_PPM column
        //                            // Do not update it
        //                            swapFiles = false;
        //                            break;
        //                        }
        //
        //                        swOutfile.WriteLine(lineIn + "\t" + SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED);
        //                        headersParsed = true;
        //                    }
        //                    else
        //                    {
        //                        int resultID = 0;
        //                        string delMPPMRefined = string.Empty;
        //
        //                        if (int.TryParse(splitLine[0], out resultID))
        //            {
        //                            double delMPPMRefined;
        //                            if (dctRefinedDelMPPMErrors.TryGetValue(resultID, out delMPPMRefined))
        //                            {
        //                                delMPPMRefined = PRISM.StringUtilities.DblToString(delMPPMRefined, 5, 0.00005);
        //                            }
        //                        }
        //
        //                        swOutfile.WriteLine(lineIn + "\t" + delMPPMRefined);
        //                    }
        //                }
        //            }
        //        }
        //
        //        if (swapFiles)
        //        {
        //            System.Threading.Thread.Sleep(150);
        //
        //            try
        //            {
        //                // Replace the original synopsis file with the updated one
        //
        //                File.Delete(synOutputFilePath);
        //                System.Threading.Thread.Sleep(150);
        //
        //                File.Move(synOutputFilePathNew, synOutputFilePath);
        //
        //                success = true;
        //
        //            }
        //            catch (Exception ex)
        //            {
        //                SetErrorMessage("Exception adding column " + SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED + " to the synopsis file: " + ex.Message);
        //                success = false;
        //            }
        //        }
        //    }
        //    catch (Exception ex)
        //    {
        //        SetErrorMessage(ex.Message);
        //        SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
        //        success = false;
        //    }
        //
        //    return success;
        //}

        private void AssociateDynamicModWithResidue(
            clsSearchResultsBaseClass searchResult,
            char mostRecentResidue,
            int residueLocInPeptide,
            string modMassDigits,
            bool updateModOccurrenceCounts)
        {
            var residueForMod = mostRecentResidue;
            var residueLocForMod = residueLocInPeptide;

            if (double.TryParse(modMassDigits, out var modMass))
            {
                if (residueLocForMod == 0)
                {
                    // Modification is at the peptide N-terminus
                    residueLocForMod = 1;
                }

                var success = searchResult.SearchResultAddModification(modMass, residueForMod, residueLocForMod,
                                                                       searchResult.DetermineResidueTerminusState(residueLocForMod),
                                                                       updateModOccurrenceCounts,
                                                                       MODA_MASS_DIGITS_OF_PRECISION);

                if (!success)
                {
                    var errorMessage = searchResult.ErrorMessage;
                    if (string.IsNullOrEmpty(errorMessage))
                    {
                        errorMessage = "SearchResultAddDynamicModification returned false for mod mass " + modMassDigits;
                    }
                    SetErrorMessage(errorMessage + "; ResultID = " + searchResult.ResultID);
                }
            }
        }

        /// <summary>
        /// Ranks each entry  assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            IList<udtMODaSearchResultType> searchResults,
            int startIndex,
            int endIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of searchResults so that we can sort by descending Probability

            var dctResultsSubset = new Dictionary<int, udtMODaSearchResultType>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                dctResultsSubset.Add(index, searchResults[index]);
            }

            var resultsByProbability = (from item in dctResultsSubset orderby item.Value.ProbabilityNum descending select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByProbability)
            {
                var oResult = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = oResult.ProbabilityNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(oResult.ProbabilityNum - lastValue) > double.Epsilon)
                    {
                        lastValue = oResult.ProbabilityNum;
                        currentRank += 1;
                    }
                }

                oResult.RankProbability = currentRank;
                searchResults[entry.Key] = oResult;
            }
        }

        private string AssureInteger(string integerText, int defaultValue)
        {
            if (integerText.EndsWith(".0"))
                integerText = integerText.Substring(0, integerText.Length - 2);

            if (int.TryParse(integerText, out var intValue))
            {
                return intValue.ToString();
            }

            if (double.TryParse(integerText, out var dblValue))
            {
                return dblValue.ToString("0");
            }

            return defaultValue.ToString();
        }

        private double ComputePeptideMass(string peptide, double totalModMass)
        {
            var cleanSequence = GetCleanSequence(peptide);

            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);

            if (Math.Abs(totalModMass) > double.Epsilon)
            {
                mass += totalModMass;
            }

            return mass;
        }

        private static readonly Regex RegexModMassRegEx = new Regex(MODA_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="peptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private double ComputeTotalModMass(string peptide)
        {
            double totalModMass = 0;

            clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptide, out var primarySequence, out _, out _);

            // Parse the dynamic mods reported by MODa
            foreach (Match reMatch in RegexModMassRegEx.Matches(primarySequence))
            {
                // We use .TrimEnd() because the matched mod mass will end in a period if this mod applies to the final residue in a peptide
                if (double.TryParse(reMatch.Groups[1].Value.TrimEnd('.'), out var modMassFound))
                {
                    totalModMass += modMassFound;
                }
            }

            // Now look for static mods
            // First determine the index of the last residue in primarySequence
            var indexLastChar = primarySequence.Length;

            for (var index = primarySequence.Length - 1; index >= 0; index += -1)
            {
                if (IsLetterAtoZ(primarySequence[index]))
                {
                    indexLastChar = index;
                    break;
                }
            }

            for (var index = 0; index <= primarySequence.Length - 1; index++)
            {
                var chChar = primarySequence[index];

                if (IsLetterAtoZ(chChar))
                {
                    // Look for static mods to associate with this residue
                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);
                            var matchFound = modificationDefinition.TargetResiduesContain(chChar);

                            if (!matchFound && index == 0)
                            {
                                matchFound = modificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
                            }

                            if (!matchFound && index == indexLastChar)
                            {
                                matchFound = modificationDefinition.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);
                            }

                            if (matchFound)
                            {
                                totalModMass += modificationDefinition.ModificationMass;
                            }
                        }
                    }
                }
            }

            return totalModMass;
        }

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputFolderPath, bool mts)
        {
            var pepToProteinMapFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
            if (pepToProteinMapFilePath != null &&
                  (pepToProteinMapFilePath.EndsWith("_MODa_syn", StringComparison.InvariantCultureIgnoreCase) ||
                   pepToProteinMapFilePath.EndsWith("_MODa_fht", StringComparison.InvariantCultureIgnoreCase)))
            {
                // Remove _syn or _fht
                pepToProteinMapFilePath = pepToProteinMapFilePath.Substring(0, pepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(pepToProteinMapFilePath, outputFolderPath, mts);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MODa
        /// The synopsis file includes every result with a probability above a set threshold
        /// The first-hits file includes the result with the highest probability (for each scan and charge)
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool CreateSynResultsFile(
            string inputFilePath,
            string outputFilePath)
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
                        var resultsProcessed = 0;

                        // Initialize the list that will hold all of the records in the MODa result file
                        var searchResultsUnfiltered = new List<udtMODaSearchResultType>();

                        // Initialize the list that will hold all of the records that will ultimately be written out to disk
                        var filteredSearchResults = new List<udtMODaSearchResultType>();

                        // Parse the input file
                        while (!reader.EndOfStream & !AbortProcessing)
                        {
                            var lineIn = reader.ReadLine();
                            var skipLine = false;

                            if (string.IsNullOrWhiteSpace(lineIn))
                            {
                                continue;
                            }

                            if (resultsProcessed == 0)
                            {
                                // The first line might be a header line
                                // However, as of April 2014, MODa id.txt files do not have a header line

                                skipLine = ParseMODaResultsFileHeaderLine(lineIn, out columnMapping);

                                // Write the header line to the output file
                                WriteSynFHTFileHeader(writer, ref errorLog);
                            }

                            if (!skipLine)
                            {
                                var udtSearchResult = new udtMODaSearchResultType();

                                var validSearchResult = ParseMODaResultsFileEntry(lineIn, ref udtSearchResult, ref errorLog, columnMapping);

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

                            resultsProcessed += 1;
                        }

                        // Sort the SearchResults by scan, charge, and Descending probability
                        searchResultsUnfiltered.Sort(new MODaSearchResultsComparerScanChargeProbabilityPeptide());

                        // Now filter the data

                        // Initialize variables
                        var startIndex = 0;

                        while (startIndex < searchResultsUnfiltered.Count)
                        {
                            var endIndex = startIndex;
                            while (endIndex + 1 < searchResultsUnfiltered.Count && searchResultsUnfiltered[endIndex + 1].ScanNum == searchResultsUnfiltered[startIndex].ScanNum)
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

        /// <summary>
        /// Load the static mods defined in the MODa parameter file
        /// </summary>
        /// <param name="mODaParamFilePath"></param>
        /// <param name="modInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool ExtractModInfoFromMODaParamFile(string mODaParamFilePath, ref List<clsModificationDefinition> modInfo)
        {
            var success = false;

            try
            {
                // Initialize the modification list
                if (modInfo == null)
                {
                    modInfo = new List<clsModificationDefinition>();
                }
                else
                {
                    modInfo.Clear();
                }

                if (string.IsNullOrEmpty(mODaParamFilePath))
                {
                    SetErrorMessage("MODa Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(mODaParamFilePath))
                {
                    SetErrorMessage("MODa param file not found: " + mODaParamFilePath);
                }
                else
                {
                    // Read the contents of the parameter (or mods) file
                    using (var reader = new StreamReader(new FileStream(mODaParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!reader.EndOfStream)
                        {
                            var lineIn = reader.ReadLine();
                            if (string.IsNullOrWhiteSpace(lineIn))
                                continue;

                            var dataLine = lineIn.Trim();
                            if (dataLine.Length <= 0)
                                continue;

                            if (dataLine.StartsWith("#"))
                            {
                                // Comment line; skip it
                                continue;
                            }

                            // Split the line on the equals sign
                            var kvSetting = clsPHRPParser.ParseKeyValueSetting(dataLine, '=', "#");

                            if (string.IsNullOrEmpty(kvSetting.Key))
                            {
                                continue;
                            }

                            if (string.Equals(kvSetting.Key, "add", StringComparison.InvariantCultureIgnoreCase))
                            {
                                // ModA defines all of its static modifications with the ADD keyword
                                // Split the value at the comma and create a new setting entry with the residue name

                                var value = kvSetting.Value;
                                var commaIndex = value.IndexOf(',');

                                var residue = value.Substring(0, commaIndex).Trim();
                                value = value.Substring(commaIndex + 1).Trim();

                                // Replace NTerm or CTerm with < or >
                                if (residue.Equals("NTerm", StringComparison.OrdinalIgnoreCase))
                                    residue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                                if (residue.Equals("CTerm", StringComparison.OrdinalIgnoreCase))
                                    residue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                                if (double.TryParse(value, out var modMass))
                                {
                                    if (Math.Abs(modMass - 0) > float.Epsilon)
                                    {
                                        var massCorrectionTag = mPeptideMods.LookupMassCorrectionTagByMass(modMass);

                                        var modDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                                                                   modMass,
                                                                                   residue,
                                                                                   clsModificationDefinition.eModificationTypeConstants.StaticMod,
                                                                                   massCorrectionTag);
                                        modInfo.Add(modDef);
                                    }
                                }
                            }
                        }
                    }

                    Console.WriteLine();

                    success = true;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MODa parameter file (" + Path.GetFileName(mODaParamFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                success = false;
            }

            return success;
        }

        private void InitializeLocalVariables()
        {
            mSpectrumIndexToScanMap = new Dictionary<int, int>();
        }

        private bool LoadMGFIndexToScanMapFile(FileInfo inputFile)
        {
            var indexToScanMapFilePath = string.Empty;

            if (inputFile.Directory == null)
            {
                ReportWarning("LoadMGFIndexToScanMapFile: Could not determine the parent directory of " + inputFile.FullName);
                return false;
            }

            try
            {
                mSpectrumIndexToScanMap.Clear();

                // Look for the IndexToScanMap file that corresponds to inputFile
                List<FileInfo> scanMapFiles;
                var matchIndex = inputFile.Name.LastIndexOf("_moda", StringComparison.InvariantCultureIgnoreCase);
                string sourceFileDescription;

                if (matchIndex > 0)
                {
                    var datasetName = inputFile.Name.Substring(0, matchIndex);
                    scanMapFiles = inputFile.Directory.GetFiles(datasetName + "*mgf_IndexToScanMap*").ToList();
                    sourceFileDescription = " dataset " + datasetName;
                }
                else
                {
                    // Source file does not have "_moda" in the name
                    // Look for any mgf_IndexToScanMap file
                    scanMapFiles = inputFile.Directory.GetFiles("*mgf_IndexToScanMap*").ToList();
                    sourceFileDescription = inputFile.Name;
                }

                if (scanMapFiles.Count == 1)
                {
                    indexToScanMapFilePath = scanMapFiles.First().FullName;
                }
                else if (scanMapFiles.Count == 0)
                {
                    ReportWarning("Did not find a mgf_IndexToScanMap file for " + sourceFileDescription + " in folder " + inputFile.Directory.FullName + "; scan numbers will be 0 in the synopsis file");
                    return false;
                }
                else
                {
                    ReportWarning("Found more than one potential mgf_IndexToScanMap file for " + sourceFileDescription + " in folder " + inputFile.Directory.FullName + " scan numbers will be 0 in the synopsis file");
                    return false;
                }

                var sourceFile = new FileInfo(indexToScanMapFilePath);

                if (!sourceFile.Exists)
                {
                    ReportWarning("MGF Index to Scan Map file not found; scan numbers will be 0 in the synopsis file: " + indexToScanMapFilePath);
                    return false;
                }

                using (var srMapFile = new StreamReader(new FileStream(sourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srMapFile.EndOfStream)
                    {
                        var lineIn = srMapFile.ReadLine();
                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        var splitLine = lineIn.Split('\t');
                        if (splitLine.Length >= 3)
                        {
                            if (int.TryParse(splitLine[0], out var spectrumIndex))
                            {
                                if (int.TryParse(splitLine[1], out var scanStart))
                                {
                                    // int.TryParse(splitLine[2], out scanEnd);

                                    mSpectrumIndexToScanMap.Add(spectrumIndex, scanStart);
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the MGF Index to Scan Map file (" + Path.GetFileName(indexToScanMapFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }

            return true;
        }

        private int LookupScanBySpectrumIndex(int spectrumIndex)
        {
            if (mSpectrumIndexToScanMap.TryGetValue(spectrumIndex, out var scanNumber))
            {
                return scanNumber;
            }

            return 0;
        }

        private bool ParseMODaSynopsisFile(
            string inputFilePath,
            string outputFolderPath,
            List<udtPepToProteinMappingType> pepToProteinMapping,
            bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that MODa synopsis files are normally sorted on Probability value, ascending
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
                var searchResult = new clsSearchResultsMODa(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize peptidesFoundForProbabilityLevel
                var peptidesFoundForProbabilityLevel = new SortedSet<string>();
                var previousProbability = string.Empty;

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
                                success = ParseMODaSynFileHeaderLine(lineIn, out columnMapping);
                                if (!success)
                                {
                                    // Error parsing header
                                    SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                                    return success;
                                }
                                headerParsed = true;
                                continue;
                            }

                            var validSearchResult = ParseMODaSynFileEntry(lineIn,
                                                                          searchResult,
                                                                          ref errorLog,
                                                                          resultsProcessed,
                                                                          columnMapping,
                                                                          out var currentPeptideWithMods);

                            if (!validSearchResult)
                            {
                                continue;
                            }

                            var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;
                            bool firstMatchForGroup;

                            if (searchResult.Probability == previousProbability)
                            {
                                // New result has the same Probability as the previous result
                                // See if peptidesFoundForProbabilityLevel contains the peptide, scan and charge

                                if (peptidesFoundForProbabilityLevel.Contains(key))
                                {
                                    firstMatchForGroup = false;
                                }
                                else
                                {
                                    peptidesFoundForProbabilityLevel.Add(key);
                                    firstMatchForGroup = true;
                                }
                            }
                            else
                            {
                                // New Probability
                                // Reset peptidesFoundForProbabilityLevel
                                peptidesFoundForProbabilityLevel.Clear();

                                // Update previousProbability
                                previousProbability = searchResult.Probability;

                                // Append a new entry to peptidesFoundForProbabilityLevel
                                peptidesFoundForProbabilityLevel.Add(key);
                                firstMatchForGroup = true;
                            }

                            success = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                            if (!success)
                            {
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
                                    } while (pepToProteinMapIndex < pepToProteinMapping.Count && currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
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

        private bool ParseMODaResultsFileEntry(
            string lineIn,
            ref udtMODaSearchResultType udtSearchResult,
            ref string errorLog,
            IList<int> columnMapping)
        {
            // Parses an entry from the MODa results file

            var rowIndex = "?";

            var validSearchResult = false;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length >= 11)
                {
                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);

                    if (!GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.SpectrumIndex], out udtSearchResult.SpectrumIndex))
                    {
                        ReportError("Index column is missing or invalid", true);
                    }
                    else
                    {
                        rowIndex = udtSearchResult.SpectrumIndex;
                    }

                    if (!int.TryParse(udtSearchResult.SpectrumIndex, out var spectrumIndex))
                    {
                        ReportError("Index column is not numeric", true);
                    }
                    udtSearchResult.ScanNum = LookupScanBySpectrumIndex(spectrumIndex);
                    if (udtSearchResult.ScanNum == 0)
                    {
                        ReportWarning("Error, could not resolve spectrumIndex to Scan Number: " + spectrumIndex);
                    }

                    // Monoisotopic mass value of the observed precursor_mz
                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.ObservedMonoMass], out udtSearchResult.Precursor_mass);
                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    // Observed m/z, converted to monoisotopic mass
                    if (double.TryParse(udtSearchResult.Precursor_mass, out var precursorMonoMass))
                    {
                        if (udtSearchResult.ChargeNum > 0)
                        {
                            var precursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(precursorMonoMass, 0, udtSearchResult.ChargeNum);
                            udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(precursorMZ, 6);
                        }
                    }

                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);

                    // Theoretical peptide monoisotopic mass, including mods, as computed by MODa
                    double.TryParse(udtSearchResult.CalculatedMonoMass, out var peptideMonoMassMODa);

                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.DeltaMass], out udtSearchResult.DeltaMass);
                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.Score], out udtSearchResult.Score);

                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.Probability], out udtSearchResult.Probability);
                    if (!double.TryParse(udtSearchResult.Probability, out udtSearchResult.ProbabilityNum))
                        udtSearchResult.ProbabilityNum = 0;

                    if (!GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                    {
                        ReportError("Peptide column is missing or invalid", true);
                    }

                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.Protein], out udtSearchResult.Protein);
                    udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);
                    GetColumnValue(splitLine, columnMapping[(int)eMODaResultsFileColumns.PeptidePosition], out udtSearchResult.PeptidePosition);

                    // Parse the sequence to determine the total mod mass
                    // Note that we do not remove any of the mod symbols since MODa identifies mods by mass alone
                    // Note that static mods are implied (thus are not explicitly displayed by MODa)
                    var totalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                    // Compute monoisotopic mass of the peptide
                    // Theoretical peptide monoisotopic mass, including mods, as computed by PHRP
                    var peptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, totalModMass);

                    if (Math.Abs(peptideMonoMassMODa) < double.Epsilon)
                    {
                        peptideMonoMassMODa = peptideMonoMassPHRP;
                    }

                    var massDiffThreshold = peptideMonoMassMODa / 50000;
                    if (massDiffThreshold < 0.1)
                        massDiffThreshold = 0.1;

                    if (Math.Abs(peptideMonoMassPHRP - peptideMonoMassMODa) > massDiffThreshold)
                    {
                        // Computed monoisotopic mass values differ by more than 0.1 Da if less than 5000 Da or by a slightly larger value if over 5000 Da; this is unexpected
                        string first30Residues;
                        if (udtSearchResult.Peptide.Length < 27)
                        {
                            first30Residues = udtSearchResult.Peptide;
                        }
                        else
                        {
                            first30Residues = udtSearchResult.Peptide.Substring(0, 27) + "...";
                        }
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " + massDiffThreshold.ToString("0.00") + " Da away from the mass computed by MODa: " + peptideMonoMassPHRP.ToString("0.0000") + " vs. " + peptideMonoMassMODa.ToString("0.0000") + "; peptide " + first30Residues);
                    }

                    if (peptideMonoMassMODa > 0)
                    {
                        // Compute DelM and DelM_PPM
                        var delM = precursorMonoMass - peptideMonoMassMODa;
                        udtSearchResult.DelM = MassErrorToString(delM);

                        var peptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(delM, precursorMonoMass, true, peptideMonoMassMODa);

                        udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(peptideDeltaMassCorrectedPpm, 5, 0.00005);
                    }

                    // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(peptideMonoMassPHRP, 0), 6);

                    if (udtSearchResult.Probability.ToLower() == "infinity")
                    {
                        udtSearchResult.Probability = "0";
                    }
                    else if (!string.IsNullOrEmpty(udtSearchResult.Probability) & !double.TryParse(udtSearchResult.Probability, out _))
                    {
                        udtSearchResult.Probability = "";
                    }

                    validSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MODa results file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        errorLog += "Error parsing MODa Results in ParseMODaResultsFileEntry for RowIndex '" + rowIndex + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MODa Results in ParseMODaResultsFileEntry" + "\n";
                    }
                }
                validSearchResult = false;
            }

            return validSearchResult;
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        /// <remarks></remarks>
        private bool ParseMODaResultsFileHeaderLine(string lineIn, out int[] columnMapping)
        {
            // Parse the header line

            // The expected column order from MODa:
            //   SpectrumFile	Index	ObservedMonoMass	Charge	CalculatedMonoMass	DeltaMass	Score	Probability	Peptide	Protein	PeptidePosition

            var columnNames = new SortedDictionary<string, eMODaResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {"SpectrumFile", eMODaResultsFileColumns.SpectrumFileName},
                {"Index", eMODaResultsFileColumns.SpectrumIndex},
                {"ObservedMW", eMODaResultsFileColumns.ObservedMonoMass},
                {"Charge", eMODaResultsFileColumns.Charge},
                {"CalculatedMW", eMODaResultsFileColumns.CalculatedMonoMass},
                {"DeltaMass", eMODaResultsFileColumns.DeltaMass},
                {"Score", eMODaResultsFileColumns.Score},
                {"Probability", eMODaResultsFileColumns.Probability},
                {"Peptide", eMODaResultsFileColumns.Peptide},
                {"Protein", eMODaResultsFileColumns.Protein},
                {"PeptidePosition", eMODaResultsFileColumns.PeptidePosition}
            };

            columnMapping = new int[MODaResultsFileColCount];

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
                            }
                            else
                            {
                                // Unrecognized column name
                                Console.WriteLine("Warning: Unrecognized column header name '" + splitLine[index] + "' in ParseMODaResultsFileHeaderLine");
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
                SetErrorMessage("Error parsing header in MODa results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        private bool ParseMODaSynFileHeaderLine(string lineIn, out int[] columnMapping)
        {
            // Parse the header line

            var columnNames = new SortedDictionary<string, eMODaSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {clsPHRPParserMODa.DATA_COLUMN_ResultID, eMODaSynFileColumns.ResultID},
                {clsPHRPParserMODa.DATA_COLUMN_Scan, eMODaSynFileColumns.Scan},
                {clsPHRPParserMODa.DATA_COLUMN_Spectrum_Index, eMODaSynFileColumns.Spectrum_Index},
                {clsPHRPParserMODa.DATA_COLUMN_Charge, eMODaSynFileColumns.Charge},
                {clsPHRPParserMODa.DATA_COLUMN_PrecursorMZ, eMODaSynFileColumns.PrecursorMZ},
                {clsPHRPParserMODa.DATA_COLUMN_DelM, eMODaSynFileColumns.DelM},
                {clsPHRPParserMODa.DATA_COLUMN_DelM_PPM, eMODaSynFileColumns.DelM_PPM},
                {clsPHRPParserMODa.DATA_COLUMN_MH, eMODaSynFileColumns.MH},
                {clsPHRPParserMODa.DATA_COLUMN_Peptide, eMODaSynFileColumns.Peptide},
                {clsPHRPParserMODa.DATA_COLUMN_Protein, eMODaSynFileColumns.Protein},
                {clsPHRPParserMODa.DATA_COLUMN_Score, eMODaSynFileColumns.Score},
                {clsPHRPParserMODa.DATA_COLUMN_Probability, eMODaSynFileColumns.Probability},
                {clsPHRPParserMODa.DATA_COLUMN_Rank_Probability, eMODaSynFileColumns.Rank_Probability},
                {clsPHRPParserMODa.DATA_COLUMN_Peptide_Position, eMODaSynFileColumns.Peptide_Position},
                {clsPHRPParserMODa.DATA_COLUMN_QValue, eMODaSynFileColumns.QValue}
            };

            columnMapping = new int[MODaSynFileColCount];

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
                SetErrorMessage("Error parsing header in MODa synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMODaSynFileEntry(
            string lineIn,
            clsSearchResultsMODa searchResult,
            ref string errorLog,
            int resultsProcessed,
            IReadOnlyList<int> columnMapping,
            out string peptideSequenceWithMods)
        {
            // Parses an entry from the MODa Synopsis file

            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();
            peptideSequenceWithMods = string.Empty;

            try
            {

                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 13)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.ResultID], out string value))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading ResultID value from MODa Results line " + (resultsProcessed + 1) +
                                       "\n";
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Charge], out string charge);

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Peptide], out peptideSequenceWithMods))
                {
                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        errorLog += "Error reading Peptide sequence value from MODa Results line " + (resultsProcessed + 1) +
                                       "\n";
                    }
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Protein], out string proteinName);
                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.DelM], out string moDaComputedDelM);
                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.DelM_PPM], out string moDaComputedDelMppm);

                searchResult.ProteinName = proteinName;
                searchResult.MODaComputedDelM = moDaComputedDelM;
                searchResult.MODaComputedDelMPPM = moDaComputedDelMppm;

                searchResult.PeptideDeltaMass = searchResult.MODaComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MODa results file as "Observed_Mass - Theoretical_Mass"
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

                var searchResultBase = searchResult as clsSearchResultsBaseClass;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Spectrum_Index], out string spectrumIndex);

                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.MH], out string parentIonMh);

                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Score], out string moDaScore);
                GetColumnValue(splitLine, columnMapping[(int)eMODaSynFileColumns.Probability], out string probability);

                searchResult.Spectrum_Index = spectrumIndex;
                searchResult.Precursor_mz = precursorMz;
                searchResult.ParentIonMH = parentIonMh;
                searchResult.MODaScore = moDaScore;
                searchResult.Probability = probability;

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing MODa Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing MODa Results in ParseMODaSynFileEntry" + "\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MODa results file (Dataset_moda.id.txt)</param>
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

                    var mODaModInfo = new List<clsModificationDefinition>();
                    var pepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the MODa Parameter File to look for any static mods
                    ExtractModInfoFromMODaParamFile(SearchToolParameterFilePath, ref mODaModInfo);

                    // Resolve the mods in mODaModInfo with the ModDefs mods
                    ResolveMODaModsWithModDefinitions(mODaModInfo);

                    // Define the base output filename using inputFilePath
                    var baseName = Path.GetFileNameWithoutExtension(inputFilePath);

                    // Auto-replace "_moda.id" with "_moda"
                    if (baseName.EndsWith("_moda.id", StringComparison.InvariantCultureIgnoreCase))
                    {
                        baseName = baseName.Substring(0, baseName.Length - "_moda.id".Length) + "_moda";
                    }

                    // Load the MSG IndexToScanMap file (if it exists)
                    LoadMGFIndexToScanMapFile(inputFile);

                    // Do not create a first-hits file for MODa results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_moda_syn.txt
                    var synOutputFilePath = Path.Combine(outputFolderPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(inputFilePath, synOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    success = ParseMODaSynopsisFile(synOutputFilePath, outputFolderPath, pepToProteinMapping, false);

                    // This step is not necessary
                    //If success Then
                    //	success = AppendDelMPPMRefinedToSynFile(synOutputFilePath)
                    //End If

                    // Remove all items from pepToProteinMapping to reduce memory overhead
                    pepToProteinMapping.Clear();
                    pepToProteinMapping.TrimExcess();

                    if (success && CreateProteinModsFile)
                    {
                        success = CreateProteinModsFileWork(baseName, inputFile, synOutputFilePath, outputFolderPath);
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in clsMODaResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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

        private bool CreateProteinModsFileWork(
            string baseName,
            FileInfo inputFile,
            string synOutputFilePath,
            string outputFolderPath)
        {
            bool success;

            if (inputFile.Directory == null)
            {
                ReportWarning("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                return false;
            }

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
                    // Auto-change mIgnorePeptideToProteinMapperErrors to True
                    // We only do this since a small number of peptides reported by MODa don't perfectly match the fasta file
                    IgnorePeptideToProteinMapperErrors = true;
                    success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                    if (!success)
                    {
                        ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (success)
            {
                if (string.IsNullOrWhiteSpace(synOutputFilePath))
                {
                    ReportWarning("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputFolderPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputFolderPath, mtsPepToProteinMapFilePath,
                                                          clsPHRPReader.ePeptideHitResultType.MODa);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
            }

            return true;
        }

        private void ResolveMODaModsWithModDefinitions(IReadOnlyCollection<clsModificationDefinition> mODaModInfo)
        {
            if (mODaModInfo != null)
            {
                // Call .LookupModificationDefinitionByMass for each entry in mODaModInfo
                foreach (var modInfo in mODaModInfo)
                {
                    if (string.IsNullOrEmpty(modInfo.TargetResidues))
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(modInfo.ModificationMass, modInfo.ModificationType, default(char), clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                    }
                    else
                    {
                        foreach (var targetResidue in modInfo.TargetResidues)
                        {
                            mPeptideMods.LookupModificationDefinitionByMassAndModType(modInfo.ModificationMass, modInfo.ModificationType, targetResidue, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                        }
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            List<udtMODaSearchResultType> filteredSearchResults,
            ref string errorLog)
        {
            // Sort udtFilteredSearchResults by descending probability, ascending scan, ascending charge, ascending peptide, and ascending protein
            filteredSearchResults.Sort(new MODaSearchResultsComparerProbabilityScanChargePeptide());

            // Compute FDR values then assign QValues
            ComputeQValues(filteredSearchResults);

            for (var index = 0; index <= filteredSearchResults.Count - 1; index++)
            {
                WriteSearchResultToFile(index + 1, writer, filteredSearchResults[index], ref errorLog);
            }
        }

        /// <summary>
        /// Compute FDR values then assign QValues
        /// </summary>
        /// <param name="searchResults"></param>
        /// <remarks>Assumes the data is sorted by descending probability using MODaSearchResultsComparerProbabilityScanChargePeptide</remarks>
        private void ComputeQValues(IList<udtMODaSearchResultType> searchResults)
        {
            var forwardPeptideCount = 0;
            var reversePeptideCount = 0;

            for (var index = 0; index < searchResults.Count;)
            {
                // Check for entries with multiple proteins listed
                var indexEnd = index;
                while (indexEnd + 1 < searchResults.Count)
                {
                    if (searchResults[index].ScanNum == searchResults[indexEnd + 1].ScanNum &&
                        searchResults[index].ChargeNum == searchResults[indexEnd + 1].ChargeNum &&
                        searchResults[index].Peptide == searchResults[indexEnd + 1].Peptide)
                    {
                        indexEnd += 1;
                    }
                    else
                    {
                        break;
                    }
                }

                var isReverse = true;

                // Look for non-reverse proteins
                for (var indexCheck = index; indexCheck <= indexEnd; indexCheck++)
                {
                    if (!IsReversedProtein(searchResults[indexCheck].Protein))
                    {
                        isReverse = false;
                        break;
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

                double fDR = 1;

                if (forwardPeptideCount > 0)
                {
                    fDR = reversePeptideCount / Convert.ToDouble(forwardPeptideCount);
                }

                // Store the FDR values
                for (var indexStore = index; indexStore <= indexEnd; indexStore++)
                {
                    var udtResult = searchResults[indexStore];
                    udtResult.FDR = fDR;

                    searchResults[indexStore] = udtResult;
                }

                index = indexEnd + 1;
            }

            // Now compute QValues
            // We step through the list, from the worst scoring result to the best result
            // The first QValue is the FDR of the final entry
            // The next QValue is the minimum of (QValue, CurrentFDR)

            var qValue = searchResults.Last().FDR;
            if (qValue > 1)
                qValue = 1;

            for (var index = searchResults.Count - 1; index >= 0; index += -1)
            {
                var udtResult = searchResults[index];

                qValue = Math.Min(qValue, udtResult.FDR);
                udtResult.QValue = qValue;

                searchResults[index] = udtResult;
            }
        }

        private void StoreSynMatches(
            IList<udtMODaSearchResultType> searchResults,
            int startIndex,
            int endIndex,
            ICollection<udtMODaSearchResultType> filteredSearchResults)
        {
            AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex);

            // The calling procedure already sorted by scan, charge, and SpecEValue; no need to re-sort

            // Now store or write out the matches that pass the filters
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].ProbabilityNum >= MODaMODPlusSynopsisFileProbabilityThreshold)
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
                    clsPHRPParserMODa.DATA_COLUMN_ResultID,
                    clsPHRPParserMODa.DATA_COLUMN_Scan,
                    clsPHRPParserMODa.DATA_COLUMN_Spectrum_Index,
                    clsPHRPParserMODa.DATA_COLUMN_Charge,
                    clsPHRPParserMODa.DATA_COLUMN_PrecursorMZ,
                    clsPHRPParserMODa.DATA_COLUMN_DelM,
                    clsPHRPParserMODa.DATA_COLUMN_DelM_PPM,
                    clsPHRPParserMODa.DATA_COLUMN_MH,
                    clsPHRPParserMODa.DATA_COLUMN_Peptide,
                    clsPHRPParserMODa.DATA_COLUMN_Protein,
                    clsPHRPParserMODa.DATA_COLUMN_Score,
                    clsPHRPParserMODa.DATA_COLUMN_Probability,
                    clsPHRPParserMODa.DATA_COLUMN_Rank_Probability,
                    clsPHRPParserMODa.DATA_COLUMN_Peptide_Position,
                    clsPHRPParserMODa.DATA_COLUMN_QValue
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
        /// Writes an entry to the synopsis file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorLog"></param>
        /// <remarks></remarks>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtMODaSearchResultType udtSearchResult,
            ref string errorLog)
        {
            try
            {
                // Primary Columns
                //
                // MODa
                // ResultID   Scan   Spectrum_Index   Charge   PrecursorMZ   DelM   DelM_PPM   MH   Peptide   Protein   Score   Probability   Rank_Probability   PeptidePosition      QValue

                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.ScanNum.ToString(),
                    udtSearchResult.SpectrumIndex,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.DelM,
                    udtSearchResult.DelM_PPM,
                    udtSearchResult.MH,
                    udtSearchResult.Peptide,
                    udtSearchResult.Protein,
                    udtSearchResult.Score,
                    udtSearchResult.Probability,
                    udtSearchResult.RankProbability.ToString(),
                    udtSearchResult.PeptidePosition,
                    PRISM.StringUtilities.DblToString(udtSearchResult.QValue, 5, 0.00005)
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

        #region "IComparer Classes"

        private class MODaSearchResultsComparerScanChargeProbabilityPeptide : IComparer<udtMODaSearchResultType>
        {
            public int Compare(udtMODaSearchResultType x, udtMODaSearchResultType y)
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

                // Charge is the same; check ProbabilityNum
                if (x.ProbabilityNum < y.ProbabilityNum)
                {
                    return 1;
                }

                if (x.ProbabilityNum > y.ProbabilityNum)
                {
                    return -1;
                }

                // Probability is the same; check peptide
                var result = string.Compare(x.Peptide, y.Peptide, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }
                return result;
            }
        }

        private class MODaSearchResultsComparerProbabilityScanChargePeptide : IComparer<udtMODaSearchResultType>
        {
            public int Compare(udtMODaSearchResultType x, udtMODaSearchResultType y)
            {
                if (x.ProbabilityNum < y.ProbabilityNum)
                {
                    return 1;
                }

                if (x.ProbabilityNum > y.ProbabilityNum)
                {
                    return -1;
                }

                // PValue is the same; check scan number
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
