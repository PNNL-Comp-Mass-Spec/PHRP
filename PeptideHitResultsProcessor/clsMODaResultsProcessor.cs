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
using System.Collections;
using System.Collections.Generic;
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
        /// For each mod mass, determine the modification and add to objSearchResult
        /// </summary>
        /// <param name="objSearchResult"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <remarks></remarks>
        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const char NO_RESIDUE = '-';

            var intModIndex = 0;
            var chChar = default(char);
            var objModificationDefinition = default(clsModificationDefinition);

            string strSequence = null;

            var blnParsingModMass = false;

            var chMostRecentResidue = default(char);
            var intResidueLocInPeptide = 0;

            blnParsingModMass = false;
            var strModMassDigits = string.Empty;

            chMostRecentResidue = NO_RESIDUE;
            intResidueLocInPeptide = 0;

            strSequence = objSearchResult.PeptideSequenceWithMods;
            for (var intIndex = 0; intIndex <= strSequence.Length - 1; intIndex++)
            {
                chChar = strSequence[intIndex];

                if (IsLetterAtoZ(chChar))
                {
                    if (blnParsingModMass)
                    {
                        // Associate the mod mass in strModMassDigits with the previous residue
                        AssociateDynamicModWithResidue(objSearchResult, chMostRecentResidue, intResidueLocInPeptide, strModMassDigits,
                            blnUpdateModOccurrenceCounts);
                        blnParsingModMass = false;
                    }

                    chMostRecentResidue = chChar;
                    intResidueLocInPeptide += 1;

                    // Look for static mods to associate with this residue
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

            var blnSuccess = false;

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

        //// This function was an experiment to compute better DelM_PPM values
        //// by reading the synopsis file with PHRPReader and re-computing the DelM_PPM values based on the monoisotopic mass values computed for the sequences
        //// It turned out to not be required, since the DelM_PPM values reported by MODa are quite accurate (despite the fact that it reports integer mod mass values)
        //private bool AppendDelMPPMRefinedToSynFile(string strSynOutputFilePath)
        //{
        //    const string SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED = "DelM_PPM_Refined";
        //
        //    bool blnSuccess = false;
        //
        //    try
        //    {
        //        // Keys in this dictionary are ResultID values from the synopsis file
        //        // Values are refined DelM_PPM values
        //        var dctRefinedDelMPPMErrors = new Dictionary<int, double>();
        //
        //        using (clsPHRPReader objReader = new clsPHRPReader(strSynOutputFilePath, clsPHRPReader.ePeptideHitResultType.MODa, blnLoadModsAndSeqInfo: true, blnLoadMSGFResults: false))
        //        {
        //            objReader.EchoMessagesToConsole = true;
        //            objReader.SkipDuplicatePSMs = true;
        //            objReader.SkipDuplicatePSMs = false;
        //
        //            foreach (string strErrorMessage in objReader.ErrorMessages)
        //            {
        //                SetErrorMessage(strErrorMessage);
        //            }
        //
        //            foreach (string strWarningMessage in objReader.WarningMessages)
        //            {
        //                ReportWarning(strWarningMessage);
        //            }
        //
        //            objReader.ClearErrors();
        //            objReader.ClearWarnings();
        //
        //            objReader.ErrorEvent += PHRPReader_ErrorEvent;
        //            objReader.WarningEvent += PHRPReader_WarningEvent;
        //
        //            while (objReader.MoveNext())
        //            {
        //                var oPSM = objReader.CurrentPSM;
        //                var oSeqInfo = objReader.CurrentPSMSeqInfo;
        //
        //                if (oSeqInfo != null)
        //                {
        //                    var dblDelM = oPSM.PrecursorNeutralMass - oSeqInfo.MonoisotopicMass;
        //
        //                    var dblPeptideDeltaMassRefinedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, oPSM.PrecursorNeutralMass, true, oSeqInfo.MonoisotopicMass);
        //
        //                    double dblOriginalDelMPPM;
        //                    if (double.TryParse(oPSM.MassErrorPPM, out dblOriginalDelMPPM))
        //                    {
        //                        if (Math.Abs(dblPeptideDeltaMassRefinedPpm - dblOriginalDelMPPM) > 2)
        //                        {
        //                            Console.WriteLine("Computed a refined DelMPPM value: " + dblPeptideDeltaMassRefinedPpm.ToString("0.0") + " vs. " + dblOriginalDelMPPM.ToString("0.0"));
        //                        }
        //                    }
        //
        //                    dctRefinedDelMPPMErrors.Add(oPSM.ResultID, dblPeptideDeltaMassRefinedPpm);
        //                }
        //            }
        //
        //            objReader.ErrorEvent -= PHRPReader_ErrorEvent;
        //            objReader.WarningEvent -= PHRPReader_WarningEvent;
        //
        //        }
        //
        //        var strSynOutputFilePathNew = strSynOutputFilePath + ".refinedDelMPPM";
        //        bool blnHeadersParsed = false;
        //        bool blnSwapFiles = true;
        //
        //        using (var srDataFile = new StreamReader(new FileStream(strSynOutputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
        //        {
        //            using (var swOutfile = new StreamWriter(new FileStream(strSynOutputFilePathNew, FileMode.Create, FileAccess.Write, FileShare.Read)))
        //            {
        //                while (!srDataFile.EndOfStream)
        //                {
        //                    var strLineIn = srDataFile.ReadLine();
        //
        //                    var strSplitLine = strLineIn.Split('\t');
        //
        //                    if (!blnHeadersParsed)
        //                    {
        //                        if (strSplitLine.Contains(SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED))
        //                        {
        //                            // This file already has the refined DelM_PPM column
        //                            // Do not update it
        //                            blnSwapFiles = false;
        //                            break;
        //                        }
        //
        //                        swOutfile.WriteLine(strLineIn + "\t" + SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED);
        //                        blnHeadersParsed = true;
        //                    }
        //                    else
        //                    {
        //                        int resultID = 0;
        //                        string strDelMPPMRefined = string.Empty;
        //
        //                        if (int.TryParse(strSplitLine[0], out resultID))
        //            {
        //                            double delMPPMRefined;
        //                            if (dctRefinedDelMPPMErrors.TryGetValue(resultID, out delMPPMRefined))
        //                            {
        //                                strDelMPPMRefined = PRISM.StringUtilities.DblToString(delMPPMRefined, 5, 0.00005);
        //                            }
        //                        }
        //
        //                        swOutfile.WriteLine(strLineIn + "\t" + strDelMPPMRefined);
        //                    }
        //                }
        //            }
        //        }
        //
        //        if (blnSwapFiles)
        //        {
        //            System.Threading.Thread.Sleep(150);
        //
        //            try
        //            {
        //                // Replace the original synopsis file with the updated one
        //
        //                File.Delete(strSynOutputFilePath);
        //                System.Threading.Thread.Sleep(150);
        //
        //                File.Move(strSynOutputFilePathNew, strSynOutputFilePath);
        //
        //                blnSuccess = true;
        //
        //            }
        //            catch (Exception ex)
        //            {
        //                SetErrorMessage("Exception adding column " + SYNOPSIS_FILE_COLUMN_DELM_PPM_REFINED + " to the synopsis file: " + ex.Message);
        //                blnSuccess = false;
        //            }
        //        }
        //    }
        //    catch (Exception ex)
        //    {
        //        SetErrorMessage(ex.Message);
        //        SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
        //        blnSuccess = false;
        //    }
        //
        //    return blnSuccess;
        //}

        private void AssociateDynamicModWithResidue(
            clsSearchResultsBaseClass objSearchResult,
            char chMostRecentResidue,
            int intResidueLocInPeptide,
            string strModMassDigits,
            bool blnUpdateModOccurrenceCounts)
        {
            var blnSuccess = false;

            var chResidueForMod = default(char);
            var intResidueLocForMod = 0;
            double dblModMass = 0;

            chResidueForMod = chMostRecentResidue;
            intResidueLocForMod = intResidueLocInPeptide;

            if (double.TryParse(strModMassDigits, out dblModMass))
            {
                if (intResidueLocForMod == 0)
                {
                    // Modification is at the peptide N-terminus
                    intResidueLocForMod = 1;
                }

                blnSuccess = objSearchResult.SearchResultAddModification(dblModMass, chResidueForMod, intResidueLocForMod, objSearchResult.DetermineResidueTerminusState(intResidueLocForMod), blnUpdateModOccurrenceCounts, MODA_MASS_DIGITS_OF_PRECISION);

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
        /// Ranks each entry  assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="lstSearchResults"></param>
        /// <param name="intStartIndex"></param>
        /// <param name="intEndIndex"></param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(
            IList<udtMODaSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex)
        {
            // Prior to September 2014 ranks were assign per charge state per scan;
            // Ranks are now assigned per scan (across all charge states)

            // Duplicate a portion of lstSearchResults so that we can sort by descending Probability

            var dctResultsSubset = new Dictionary<int, udtMODaSearchResultType>();
            for (var intIndex = intStartIndex; intIndex <= intEndIndex; intIndex++)
            {
                dctResultsSubset.Add(intIndex, lstSearchResults[intIndex]);
            }

            var lstResultsByProbability = (from item in dctResultsSubset orderby item.Value.ProbabilityNum descending select item).ToList();

            double dblLastValue = 0;
            var intCurrentRank = -1;

            foreach (var entry in lstResultsByProbability)
            {
                var oResult = lstSearchResults[entry.Key];

                if (intCurrentRank < 0)
                {
                    dblLastValue = oResult.ProbabilityNum;
                    intCurrentRank = 1;
                }
                else
                {
                    if (Math.Abs(oResult.ProbabilityNum - dblLastValue) > double.Epsilon)
                    {
                        dblLastValue = oResult.ProbabilityNum;
                        intCurrentRank += 1;
                    }
                }

                oResult.RankProbability = intCurrentRank;
                lstSearchResults[entry.Key] = oResult;
            }
        }

        private string AssureInteger(string strInteger, int intDefaultValue)
        {
            var intValue = 0;
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

        private double ComputePeptideMass(string strPeptide, double dblTotalModMass)
        {
            var strCleanSequence = GetCleanSequence(strPeptide);

            var dblMass = mPeptideSeqMassCalculator.ComputeSequenceMass(strCleanSequence);

            if (Math.Abs(dblTotalModMass) > double.Epsilon)
            {
                dblMass += dblTotalModMass;
            }

            return dblMass;
        }

        private static readonly Regex RegexModMassRegEx = new Regex(MODA_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Computes the total of all modification masses defined for the peptide
        /// </summary>
        /// <param name="strPeptide">Peptide sequence, with mod masses in the form +53.8 or -23</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private double ComputeTotalModMass(string strPeptide)
        {
            double dblTotalModMass = 0;

            var strPrimarySequence = string.Empty;
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;

            clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, out strPrimarySequence, out strPrefix, out strSuffix);

            // Parse the dynamic mods reported by MODa
            foreach (Match reMatch in RegexModMassRegEx.Matches(strPrimarySequence))
            {
                double dblModMassFound = 0;
                // We use .TrimEnd() because the matched mod mass will end in a period if this mod applies to the final residue in a peptide
                if (double.TryParse(reMatch.Groups[1].Value.TrimEnd('.'), out dblModMassFound))
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
            if (strPepToProteinMapFilePath.EndsWith("_MODa_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_MODa_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS);
        }

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from MODa
        /// The synopsis file includes every result with a probability above a set threshold
        /// The first-hits file includes the result with the highest probability (for each scan and charge)
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool CreateSynResultsFile(
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
                        var intResultsProcessed = 0;

                        // Initialize the list that will hold all of the records in the MODa result file
                        var lstSearchResultsUnfiltered = new List<udtMODaSearchResultType>();

                        // Initialize the list that will hold all of the records that will ultimately be written out to disk
                        var lstFilteredSearchResults = new List<udtMODaSearchResultType>();

                        // Parse the input file
                        while (!srDataFile.EndOfStream & !AbortProcessing)
                        {
                            var strLineIn = srDataFile.ReadLine();
                            var blnSkipLine = false;

                            if (string.IsNullOrWhiteSpace(strLineIn))
                            {
                                continue;
                            }

                            if (intResultsProcessed == 0)
                            {
                                // The first line might be a header line
                                // However, as of April 2014, MODa id.txt files do not have a header line

                                blnSkipLine = ParseMODaResultsFileHeaderLine(strLineIn, out intColumnMapping);

                                // Write the header line to the output file
                                WriteSynFHTFileHeader(swResultFile, ref strErrorLog);
                            }

                            if (!blnSkipLine)
                            {
                                var udtSearchResult = new udtMODaSearchResultType();

                                var blnValidSearchResult = ParseMODaResultsFileEntry(strLineIn, ref udtSearchResult, ref strErrorLog, intColumnMapping);

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

                            intResultsProcessed += 1;
                        }

                        // Sort the SearchResults by scan, charge, and Descending probability
                        lstSearchResultsUnfiltered.Sort(new MODaSearchResultsComparerScanChargeProbabilityPeptide());

                        // Now filter the data

                        // Initialize variables
                        var intStartIndex = 0;
                        var intEndIndex = 0;

                        while (intStartIndex < lstSearchResultsUnfiltered.Count)
                        {
                            intEndIndex = intStartIndex;
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
        /// Load the static mods defined in the MODa parameter file
        /// </summary>
        /// <param name="strMODaParamFilePath"></param>
        /// <param name="lstModInfo"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool ExtractModInfoFromMODaParamFile(string strMODaParamFilePath, ref List<clsModificationDefinition> lstModInfo)
        {
            var kvSetting = default(KeyValuePair<string, string>);

            var objModDef = default(clsModificationDefinition);

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

                if (string.IsNullOrEmpty(strMODaParamFilePath))
                {
                    SetErrorMessage("MODa Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                if (!File.Exists(strMODaParamFilePath))
                {
                    SetErrorMessage("MODa param file not found: " + strMODaParamFilePath);
                }
                else
                {
                    // Read the contents of the parameter (or mods) file
                    using (var srInFile = new StreamReader(new FileStream(strMODaParamFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        while (!srInFile.EndOfStream)
                        {
                            var lineIn = srInFile.ReadLine();
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
                            kvSetting = clsPHRPParser.ParseKeyValueSetting(dataLine, '=', "#");

                            if (string.IsNullOrEmpty(kvSetting.Key))
                            {
                                continue;
                            }

                            if (string.Equals(kvSetting.Key, "add", StringComparison.InvariantCultureIgnoreCase))
                            {
                                // ModA defines all of its static modifications with the ADD keyword
                                // Split the value at the comma and create a new setting entry with the residue name

                                var strValue = kvSetting.Value;
                                var commaIndex = strValue.IndexOf(',');

                                var strResidue = strValue.Substring(0, commaIndex).Trim();
                                strValue = strValue.Substring(commaIndex + 1).Trim();

                                // Replace NTerm or CTerm with < or >
                                if (strResidue.ToLower() == "nterm")
                                    strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                                if (strResidue.ToLower() == "cterm")
                                    strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                                double modMass = 0;
                                if (double.TryParse(strValue, out modMass))
                                {
                                    if (Math.Abs(modMass - 0) > float.Epsilon)
                                    {
                                        var strMassCorrectionTag = mPeptideMods.LookupMassCorrectionTagByMass(modMass);

                                        objModDef = new clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMass, strResidue, clsModificationDefinition.eModificationTypeConstants.StaticMod, strMassCorrectionTag);
                                        lstModInfo.Add(objModDef);
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
                SetErrorMessage("Error reading the MODa parameter file (" + Path.GetFileName(strMODaParamFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private void InitializeLocalVariables()
        {
            mSpectrumIndexToScanMap = new Dictionary<int, int>();
        }

        private bool LoadMGFIndexToScanMapFile(FileInfo inputFile)
        {
            var indexToScanMapFilePath = string.Empty;

            try
            {
                mSpectrumIndexToScanMap.Clear();

                // Look for the IndexToScanMap file that corresponds to inputFile
                List<FileInfo> lstScanMapFiles;
                var matchIndex = inputFile.Name.LastIndexOf("_moda", StringComparison.InvariantCultureIgnoreCase);
                string sourceFileDescription;

                if (matchIndex > 0)
                {
                    var datasetName = inputFile.Name.Substring(0, matchIndex);
                    lstScanMapFiles = inputFile.Directory.GetFiles(datasetName + "*mgf_IndexToScanMap*").ToList();
                    sourceFileDescription = " dataset " + datasetName;
                }
                else
                {
                    // Source file does not have "_moda" in the name
                    // Look for any mgf_IndexToScanMap file
                    lstScanMapFiles = inputFile.Directory.GetFiles("*mgf_IndexToScanMap*").ToList();
                    sourceFileDescription = inputFile.Name;
                }

                if (lstScanMapFiles.Count == 1)
                {
                    indexToScanMapFilePath = lstScanMapFiles.First().FullName;
                }
                else if (lstScanMapFiles.Count == 0)
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
                        var strLineIn = srMapFile.ReadLine();
                        if (string.IsNullOrEmpty(strLineIn))
                            continue;

                        var strSplitLine = strLineIn.Split('\t');
                        if (strSplitLine.Length >= 3)
                        {
                            var spectrumIndex = 0;
                            var scanStart = 0;
                            var scanEnd = 0;

                            if (int.TryParse(strSplitLine[0], out spectrumIndex))
                            {
                                if (int.TryParse(strSplitLine[1], out scanStart))
                                {
                                    int.TryParse(strSplitLine[2], out scanEnd);

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
            var scanNumber = 0;
            if (mSpectrumIndexToScanMap.TryGetValue(spectrumIndex, out scanNumber))
            {
                return scanNumber;
            }

            return 0;
        }

        private bool ParseMODaSynopsisFile(
            string strInputFilePath,
            string strOutputFolderPath,
            List<udtPepToProteinMappingType> lstPepToProteinMapping,
            bool blnResetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            string strPreviousProbability = null;

            // Note that MODa synopsis files are normally sorted on Probability value, ascending
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
                var objSearchResult = new clsSearchResultsMODa(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForProbabilityLevel
                var htPeptidesFoundForProbabilityLevel = new Hashtable();
                strPreviousProbability = string.Empty;

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
                                blnSuccess = ParseMODaSynFileHeaderLine(strLineIn, out intColumnMapping);
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

                            var blnValidSearchResult = ParseMODaSynFileEntry(strLineIn,
                                                                             objSearchResult,
                                                                             ref strErrorLog,
                                                                             intResultsProcessed,
                                                                             intColumnMapping,
                                                                             out strCurrentPeptideWithMods);

                            if (!blnValidSearchResult)
                            {
                                continue;
                            }

                            var strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.Charge;
                            var blnFirstMatchForGroup = false;

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

                            SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass) objSearchResult, blnFirstMatchForGroup);

                            if (lstPepToProteinMapping.Count > 0)
                            {
                                // Add the additional proteins for this peptide

                                // Use binary search to find this peptide in lstPepToProteinMapping
                                var intPepToProteinMapIndex = FindFirstMatchInPepToProteinMapping(lstPepToProteinMapping, strCurrentPeptideWithMods);

                                if (intPepToProteinMapIndex >= 0)
                                {
                                    // Call MyBase.SaveResultsFileEntrySeqInfo for each entry in lstPepToProteinMapping() for peptide , skipping objSearchResult.ProteinName
                                    var strCurrentProtein = string.Copy(objSearchResult.ProteinName);
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

        private bool ParseMODaResultsFileEntry(
            string strLineIn,
            ref udtMODaSearchResultType udtSearchResult,
            ref string strErrorLog,
            IList<int> intColumnMapping)
        {
            // Parses an entry from the MODa results file

            var rowIndex = "?";
            double dblPrecursorMonoMass = 0;           // Observed m/z, converted to monoisotopic mass
            double dblPeptideMonoMassMODa = 0;         // Theoretical peptide monoisotopic mass, including mods, as computed by MODa
            double dblPeptideMonoMassPHRP = 0;         // Theoretical peptide monoisotopic mass, including mods, as computed by PHRP

            double dblPrecursorMZ = 0;
            double dblDelM = 0;

            double dblTotalModMass = 0;

            var blnValidSearchResult = false;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                udtSearchResult.Clear();
                var strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 11)
                {
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.SpectrumFileName], out udtSearchResult.SpectrumFileName);

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.SpectrumIndex], out udtSearchResult.SpectrumIndex))
                    {
                        ReportError("Index column is missing or invalid", true);
                    }
                    else
                    {
                        rowIndex = udtSearchResult.SpectrumIndex;
                    }

                    var spectrumIndex = 0;
                    if (!int.TryParse(udtSearchResult.SpectrumIndex, out spectrumIndex))
                    {
                        ReportError("Index column is not numeric", true);
                    }
                    udtSearchResult.ScanNum = LookupScanBySpectrumIndex(spectrumIndex);
                    if (udtSearchResult.ScanNum == 0)
                    {
                        ReportWarning("Error, could not resolve spectrumIndex to Scan Number: " + spectrumIndex);
                    }

                    // Monoisotopic mass value of the observed precursor_mz
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.ObservedMonoMass], out udtSearchResult.Precursor_mass);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.Charge], out udtSearchResult.Charge);
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    if (double.TryParse(udtSearchResult.Precursor_mass, out dblPrecursorMonoMass))
                    {
                        if (udtSearchResult.ChargeNum > 0)
                        {
                            dblPrecursorMZ = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMonoMass, 0, udtSearchResult.ChargeNum);
                            udtSearchResult.PrecursorMZ = PRISM.StringUtilities.DblToString(dblPrecursorMZ, 6);
                        }
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass);
                    double.TryParse(udtSearchResult.CalculatedMonoMass, out dblPeptideMonoMassMODa);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.DeltaMass], out udtSearchResult.DeltaMass);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.Score], out udtSearchResult.Score);

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.Probability], out udtSearchResult.Probability);
                    if (!double.TryParse(udtSearchResult.Probability, out udtSearchResult.ProbabilityNum))
                        udtSearchResult.ProbabilityNum = 0;

                    if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.Peptide], out udtSearchResult.Peptide))
                    {
                        ReportError("Peptide column is missing or invalid", true);
                    }

                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.Protein], out udtSearchResult.Protein);
                    udtSearchResult.Protein = TruncateProteinName(udtSearchResult.Protein);
                    GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaResultsFileColumns.PeptidePosition], out udtSearchResult.PeptidePosition);

                    // Parse the sequence to determine the total mod mass
                    // Note that we do not remove any of the mod symbols since MODa identifies mods by mass alone
                    // Note that static mods are implied (thus are not explicitly displayed by MODa)
                    dblTotalModMass = ComputeTotalModMass(udtSearchResult.Peptide);

                    // Compute monoisotopic mass of the peptide
                    dblPeptideMonoMassPHRP = ComputePeptideMass(udtSearchResult.Peptide, dblTotalModMass);

                    if (Math.Abs(dblPeptideMonoMassMODa) < double.Epsilon)
                    {
                        dblPeptideMonoMassMODa = dblPeptideMonoMassPHRP;
                    }

                    var dblMassDiffThreshold = dblPeptideMonoMassMODa / 50000;
                    if (dblMassDiffThreshold < 0.1)
                        dblMassDiffThreshold = 0.1;

                    if (Math.Abs(dblPeptideMonoMassPHRP - dblPeptideMonoMassMODa) > dblMassDiffThreshold)
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
                        ReportWarning("The monoisotopic mass computed by PHRP is more than " + dblMassDiffThreshold.ToString("0.00") + " Da away from the mass computed by MODa: " + dblPeptideMonoMassPHRP.ToString("0.0000") + " vs. " + dblPeptideMonoMassMODa.ToString("0.0000") + "; peptide " + strFirst30Residues);
                    }

                    if (dblPeptideMonoMassMODa > 0)
                    {
                        // Compute DelM and DelM_PPM
                        dblDelM = dblPrecursorMonoMass - dblPeptideMonoMassMODa;
                        udtSearchResult.DelM = MassErrorToString(dblDelM);

                        var dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, dblPrecursorMonoMass, true, dblPeptideMonoMassMODa);

                        udtSearchResult.DelM_PPM = PRISM.StringUtilities.DblToString(dblPeptideDeltaMassCorrectedPpm, 5, 0.00005);
                    }

                    // Store the monoisotopic MH value in .MH; note that this is (M+H)+
                    udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(dblPeptideMonoMassPHRP, 0), 6);

                    double dblProbability = 0;

                    if (udtSearchResult.Probability.ToLower() == "infinity")
                    {
                        udtSearchResult.Probability = "0";
                    }
                    else if (!string.IsNullOrEmpty(udtSearchResult.Probability) & !double.TryParse(udtSearchResult.Probability, out dblProbability))
                    {
                        udtSearchResult.Probability = "";
                    }

                    blnValidSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the MODa results file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (!string.IsNullOrEmpty(rowIndex))
                    {
                        strErrorLog += "Error parsing MODa Results in ParseMODaResultsFileEntry for RowIndex '" + rowIndex + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MODa Results in ParseMODaResultsFileEntry" + "\n";
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
        private bool ParseMODaResultsFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            // The expected column order from MODa:
            //   SpectrumFile	Index	ObservedMonoMass	Charge	CalculatedMonoMass	DeltaMass	Score	Probability	Peptide	Protein	PeptidePosition

            var lstColumnNames = new SortedDictionary<string, eMODaResultsFileColumns>(StringComparer.InvariantCultureIgnoreCase)
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

            intColumnMapping = new int[MODaResultsFileColCount];

            try
            {
                // Initialize each entry in intColumnMapping to -1
                for (var intIndex = 0; intIndex <= intColumnMapping.Length - 1; intIndex++)
                {
                    intColumnMapping[intIndex] = -1;
                }

                var strSplitLine = strLineIn.Split('\t');
                var blnUseDefaultHeaders = false;

                var value = 0;
                if (strSplitLine.Length >= 2)
                {
                    if (int.TryParse(strSplitLine[1], out value))
                    {
                        // Second column has a number; this is not a header line
                        blnUseDefaultHeaders = true;
                    }
                    else
                    {
                        for (var intIndex = 0; intIndex <= strSplitLine.Length - 1; intIndex++)
                        {
                            var eResultFileColumn = default(eMODaResultsFileColumns);

                            if (lstColumnNames.TryGetValue(strSplitLine[intIndex], out eResultFileColumn))
                            {
                                // Recognized column name; update intColumnMapping
                                intColumnMapping[(int)eResultFileColumn] = intIndex;
                                blnUseDefaultHeaders = false;
                            }
                            else
                            {
                                // Unrecognized column name
                                Console.WriteLine("Warning: Unrecognized column header name '" + strSplitLine[intIndex] + "' in ParseMODaResultsFileHeaderLine");
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
                SetErrorMessage("Error parsing header in MODa results file: " + ex.Message);
                return false;
            }

            // Header line found and parsed; return true
            return true;
        }

        private bool ParseMODaSynFileHeaderLine(string strLineIn, out int[] intColumnMapping)
        {
            // Parse the header line

            var lstColumnNames = new SortedDictionary<string, eMODaSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
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

            intColumnMapping = new int[MODaSynFileColCount];

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
                SetErrorMessage("Error parsing header in MODa synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        private bool ParseMODaSynFileEntry(
            string strLineIn,
            clsSearchResultsMODa objSearchResult,
            ref string strErrorLog,
            int intResultsProcessed,
            IReadOnlyList<int> intColumnMapping,
            out string strPeptideSequenceWithMods)
        {
            // Parses an entry from the MODa Synopsis file

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

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.ResultID], out string strValue))
                {
                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        strErrorLog += "Error reading ResultID value from MODa Results line " + (intResultsProcessed + 1) +
                                       "\n";
                    }
                    return false;
                }

                objSearchResult.ResultID = int.Parse(strValue);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Scan], out string scan);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Charge], out string charge);

                objSearchResult.Scan = scan;
                objSearchResult.Charge = charge;

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Peptide], out strPeptideSequenceWithMods))
                {
                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                    {
                        strErrorLog += "Error reading Peptide sequence value from MODa Results line " + (intResultsProcessed + 1) +
                                       "\n";
                    }
                    return false;
                }

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Protein], out string proteinName);
                objSearchResult.MultipleProteinCount = "0";

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.DelM], out string moDaComputedDelM);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.DelM_PPM], out string moDaComputedDelMppm);

                objSearchResult.ProteinName = proteinName;
                objSearchResult.MODaComputedDelM = moDaComputedDelM;
                objSearchResult.MODaComputedDelMPPM = moDaComputedDelMppm;

                objSearchResult.PeptideDeltaMass = objSearchResult.MODaComputedDelM;

                // Note: .PeptideDeltaMass is stored in the MODa results file as "Observed_Mass - Theoretical_Mass"
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

                var objSearchResultBase = default(clsSearchResultsBaseClass);
                objSearchResultBase = (clsSearchResultsBaseClass) objSearchResult;

                ComputePseudoPeptideLocInProtein(objSearchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                objSearchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Spectrum_Index], out string spectrumIndex);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.PrecursorMZ], out string precursorMz);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.MH], out string parentIonMh);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Score], out string moDaScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eMODaSynFileColumns.Probability], out string probability);

                objSearchResult.Spectrum_Index = spectrumIndex;
                objSearchResult.Precursor_mz = precursorMz;
                objSearchResult.ParentIonMH = parentIonMh;
                objSearchResult.MODaScore = moDaScore;
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
                        strErrorLog += "Error parsing MODa Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing MODa Results in ParseMODaSynFileEntry" + "\n";
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">MODa results file (Dataset_moda.id.txt)</param>
        /// <param name="strOutputFolderPath">Output folder</param>
        /// <param name="strParameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            var strBaseName = string.Empty;
            var strSynOutputFilePath = string.Empty;

            var lstMODaModInfo = default(List<clsModificationDefinition>);
            var lstPepToProteinMapping = default(List<udtPepToProteinMappingType>);

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

                    lstMODaModInfo = new List<clsModificationDefinition>();
                    lstPepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the MODa Parameter File to look for any static mods
                    ExtractModInfoFromMODaParamFile(SearchToolParameterFilePath, ref lstMODaModInfo);

                    // Resolve the mods in lstMODaModInfo with the ModDefs mods
                    ResolveMODaModsWithModDefinitions(lstMODaModInfo);

                    // Define the base output filename using strInputFilePath
                    strBaseName = Path.GetFileNameWithoutExtension(strInputFilePath);

                    // Auto-replace "_moda.id" with "_moda"
                    if (strBaseName.EndsWith("_moda.id", StringComparison.InvariantCultureIgnoreCase))
                    {
                        strBaseName = strBaseName.Substring(0, strBaseName.Length - "_moda.id".Length) + "_moda";
                    }

                    // Load the MSG IndexToScanMap file (if it exists)
                    LoadMGFIndexToScanMapFile(inputFile);

                    // Do not create a first-hits file for MODa results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form BasePath_moda_syn.txt
                    strSynOutputFilePath = Path.Combine(strOutputFolderPath, strBaseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    blnSuccess = CreateSynResultsFile(strInputFilePath, strSynOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(strSynOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to next create the other PHRP files
                    blnSuccess = ParseMODaSynopsisFile(strSynOutputFilePath, strOutputFolderPath, lstPepToProteinMapping, false);

                    // This step is not necessary
                    //If blnSuccess Then
                    //	blnSuccess = AppendDelMPPMRefinedToSynFile(strSynOutputFilePath)
                    //End If

                    // Remove all items from lstPepToProteinMapping to reduce memory overhead
                    lstPepToProteinMapping.Clear();
                    lstPepToProteinMapping.TrimExcess();

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
                    SetErrorMessage("Error in clsMODaResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
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
            var blnSuccess = false;
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
                    // We only do this since a small number of peptides reported by MODa don't perfectly match the fasta file
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
                blnSuccess = CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.MODa);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                blnSuccess = true;
            }
            return blnSuccess;
        }

        private void ResolveMODaModsWithModDefinitions(IReadOnlyCollection<clsModificationDefinition> lstMODaModInfo)
        {
            if (lstMODaModInfo != null)
            {
                // Call .LookupModificationDefinitionByMass for each entry in lstMODaModInfo
                foreach (var objModInfo in lstMODaModInfo)
                {
                    if (string.IsNullOrEmpty(objModInfo.TargetResidues))
                    {
                        mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, default(char), clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                    }
                    else
                    {
                        foreach (var chTargetResidue in objModInfo.TargetResidues)
                        {
                            mPeptideMods.LookupModificationDefinitionByMassAndModType(objModInfo.ModificationMass, objModInfo.ModificationType, chTargetResidue, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, out _, true, MODA_MASS_DIGITS_OF_PRECISION);
                        }
                    }
                }
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter swResultFile,
            List<udtMODaSearchResultType> lstFilteredSearchResults,
            ref string strErrorLog)
        {
            // Sort udtFilteredSearchResults by descending probability, ascending scan, ascending charge, ascending peptide, and ascending protein
            lstFilteredSearchResults.Sort(new MODaSearchResultsComparerProbabilityScanChargePeptide());

            // Compute FDR values then assign QValues
            ComputeQValues(lstFilteredSearchResults);

            for (var intIndex = 0; intIndex <= lstFilteredSearchResults.Count - 1; intIndex++)
            {
                WriteSearchResultToFile(intIndex + 1, swResultFile, lstFilteredSearchResults[intIndex], ref strErrorLog);
            }
        }

        /// <summary>
        /// Compute FDR values then assign QValues
        /// </summary>
        /// <param name="lstSearchResults"></param>
        /// <remarks>Assumes the data is sorted by descending probability using MODaSearchResultsComparerProbabilityScanChargePeptide</remarks>
        private void ComputeQValues(IList<udtMODaSearchResultType> lstSearchResults)
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
                    if (!IsReversedProtein(lstSearchResults[intIndexCheck].Protein))
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
            IList<udtMODaSearchResultType> lstSearchResults,
            int intStartIndex,
            int intEndIndex,
            IList<udtMODaSearchResultType> lstFilteredSearchResults)
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
            TextWriter swResultFile,
            udtMODaSearchResultType udtSearchResult,
            ref string strErrorLog)
        {
            try
            {
                // Primary Columns
                //
                // MODa
                // ResultID   Scan   Spectrum_Index   Charge   PrecursorMZ   DelM   DelM_PPM   MH   Peptide   Protein   Score   Probability   Rank_Probability   PeptidePosition      QValue

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
                    udtSearchResult.Protein,
                    udtSearchResult.Score,
                    udtSearchResult.Probability,
                    udtSearchResult.RankProbability.ToString(),
                    udtSearchResult.PeptidePosition,
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

        private class MODaSearchResultsComparerScanChargeProbabilityPeptide : IComparer<udtMODaSearchResultType>
        {
            public int Compare(udtMODaSearchResultType x, udtMODaSearchResultType y)
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
                        // Charge is the same; check ProbabilityNum
                        if (x.ProbabilityNum < y.ProbabilityNum)
                        {
                            return 1;
                        }
                        else if (x.ProbabilityNum > y.ProbabilityNum)
                        {
                            return -1;
                        }
                        else
                        {
                            // Probability is the same; check peptide
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
                else if (x.ProbabilityNum > y.ProbabilityNum)
                {
                    return -1;
                }
                else
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
            }
        }

        #endregion
    }
}
