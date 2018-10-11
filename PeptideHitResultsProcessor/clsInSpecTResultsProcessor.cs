// -------------------------------------------------------------------------------
// Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
// Program started August 12, 2008
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
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text.RegularExpressions;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class reads in an InSpecT results file (txt format) and creates
    /// a tab-delimited text file with the data.  It will insert modification symbols
    /// into the peptide sequences for modified peptides.
    /// </summary>
    /// <remarks>The modification definition information is determined from the InSpecT parameter file</remarks>
    public class clsInSpecTResultsProcessor : clsPHRPBaseClass
    {
        public clsInSpecTResultsProcessor()
        {
            mFileDate = "October 13, 2017";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"

        public const string FILENAME_SUFFIX_INSPECT_FILE = "_inspect";

        private const int INSPECT_SYN_FILE_MIN_COL_COUNT = 5;

        public const string N_TERMINUS_SYMBOL_INSPECT = "*.";
        public const string C_TERMINUS_SYMBOL_INSPECT = ".*";

        private const char UNKNOWN_INSPECT_MOD_SYMBOL = '?';

        // When writing the synopsis file, we keep data that passes any of these thresholds (thus, it's an OR comparison, not an AND comparison)
        // pValue <= 0.2 Or TotalPRMScore >= 50 or FScore >= 0
        public const float DEFAULT_SYN_FILE_PVALUE_THRESHOLD = 0.2f;
        public const float TOTALPRMSCORE_THRESHOLD = 50;
        public const float FSCORE_THRESHOLD = 0;

        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string DTA_FILENAME_SCAN_NUMBER_REGEX = @"(\d+)\.\d+\.\d+\.dta";
        private const string INSPECT_NTERMINAL_MOD_MASS_REGEX = @"^\+(\d+)";
        private const string INSPECT_CTERMINAL_MOD_MASS_REGEX = @"\+(\d+)$";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        private const string PHOS_MOD_NAME = "phos";
        private const string PHOS_MOD_MASS = "79.9663";
        private const string PHOS_MOD_RESIDUES = "STY";

        // These columns correspond to the tab-delimited file created directly by Inspect
        public enum eInspectResultsFileColumns
        {
            SpectrumFile = 0,
            Scan = 1,
            Annotation = 2,
            Protein = 3,
            Charge = 4,
            MQScore = 5,
            Length = 6,
            TotalPRMScore = 7,
            MedianPRMScore = 8,
            FractionY = 9,
            FractionB = 10,
            Intensity = 11,
            NTT = 12,
            PValue = 13,
            FScore = 14,
            DeltaScore = 15,
            DeltaScoreOther = 16,
            RecordNumber = 17,
            DBFilePos = 18,
            SpecFilePos = 19,
            PrecursorMZ = 20,
            PrecursorError = 21
        }

        // These columns correspond to the Synopsis and First-Hits files created by this class
        protected const int InspectSynopsisFileColCount = 27;
        public enum eInspectSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Peptide = 2,
            Protein = 3,
            Charge = 4,
            MQScore = 5,
            Length = 6,
            TotalPRMScore = 7,
            MedianPRMScore = 8,
            FractionY = 9,
            FractionB = 10,
            Intensity = 11,
            NTT = 12,
            PValue = 13,
            FScore = 14,
            DeltaScore = 15,
            DeltaScoreOther = 16,
            DeltaNormMQScore = 17,                   // Computed as Abs((MQScore(n) - MQScore(n+1)) / MQScore(n)); storing 0 for the lowest scoring result in each set. If MQScore(n) is 0, then also storing 0.   This value is not usable when MQScore(n) is <= 0, and should generally not be used when MQScore(n) is < 0.5
            DeltaNormTotalPRMScore = 18,             // Computed as Abs((TotalPRMScore(n) - TotalPRMScore(n+1)) / TotalPRMScore(n)); storing 0 for the lowest scoring result in each set.  If TotalPRMScore(n) is 0, then also storing 0.  This value is not usable when TotalPRMScore(n) is <= 0, and should generally not be used when TotalPRMScore(n) is < 0.5
            RankTotalPRMScore = 19,                  // Rank 1 means highest TotalPRMScore, 2 means next lower score, etc. (ties get the same rank)
            RankFScore = 20,                         // Rank 1 means highest FScore, 2 means next lower, etc. (ties get the same rank)
            MH = 21,                                 // Theoretical monoisotopic peptide mass (computed by PHRP); note that this is (M+H)+
            RecordNumber = 22,
            DBFilePos = 23,
            SpecFilePos = 24,
            PrecursorMZ = 25,
            PrecursorError = 26
        }

        protected enum eInspectModType
        {
            Unknown = 0,
            DynamicMod = 1,
            StaticMod = 2,
            DynNTermPeptide = 3,
            DynCTermPeptide = 4
        }

        protected enum eFilteredOutputFileTypeConstants
        {
            SynFile = 0,
            FHTbyFScore = 1,
            FHTbyTotalPRM = 2
        }
        #endregion

        #region "Structures"
        protected struct udtInspectSearchResultType
        {
            public string SpectrumFileName;
            public string Scan;
            public int ScanNum;
            public string PeptideAnnotation;
            public string Protein;
            public string Charge;
            public short ChargeNum;
            public string MQScore;                  // Higher values are better scores; note that MQScore can be negative
            public float MQScoreNum;                // Store the value of the string for quick reference when sorting
            public int Length;
            public string TotalPRMScore;            // Higher values are better scores
            public float TotalPRMScoreNum;          // We store the value of the string for quick reference when sorting
            public string MedianPRMScore;
            public string FractionY;
            public string FractionB;
            public string Intensity;
            public int NTT;
            public string pValue;                   // Lower values are better scores
            public float PValueNum;                 // Store the value of the string for quick reference when sorting
            public string FScore;                   // Higher values are better scores
            public float FScoreNum;                 // Store the value of the string for quick reference when sorting
            public string DeltaScore;
            public string DeltaScoreOther;
            public float DeltaNormMQScore;
            public float DeltaNormTotalPRMScore;
            public int RankTotalPRMScore;
            public int RankFScore;
            public double MH;
            public string RecordNumber;
            public string DBFilePos;
            public string SpecFilePos;
            public string PrecursorMZ;
            public string PrecursorError;           // Precursor error in; units are m/z (NOT Daltons)
            public string DelMPPM;                  // Computed by this application

            public void Clear()
            {
                SpectrumFileName = string.Empty;
                Scan = string.Empty;
                ScanNum = 0;
                PeptideAnnotation = string.Empty;
                Protein = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                MQScore = string.Empty;
                MQScoreNum = 0;
                Length = 0;
                TotalPRMScore = string.Empty;
                TotalPRMScoreNum = 0;
                MedianPRMScore = string.Empty;
                FractionY = string.Empty;
                FractionB = string.Empty;
                Intensity = string.Empty;
                NTT = 0;
                pValue = string.Empty;
                PValueNum = 0;
                FScore = string.Empty;
                FScoreNum = 0;
                DeltaScore = string.Empty;
                DeltaScoreOther = string.Empty;
                DeltaNormMQScore = 0;
                DeltaNormTotalPRMScore = 0;
                RankTotalPRMScore = 0;
                RankFScore = 0;
                MH = 0;
                RecordNumber = string.Empty;
                DBFilePos = string.Empty;
                SpecFilePos = string.Empty;
                PrecursorMZ = string.Empty;
                PrecursorError = string.Empty;
                DelMPPM = string.Empty;
            }
        }

        protected struct udtModInfoType
        {
            public string ModName;              // Mod names must be lower case, and 4 characters in length (or shorter)
            public string ModMass;              // Storing as a string since reading from a text file and writing to a text file
            public string Residues;
            public eInspectModType ModType;
            public string ModSymbol;
        }

        #endregion

        #region "Classwide Variables"

        #endregion

        #region "Properties"

        public bool SortFHTAndSynFiles { get; set; }

        #endregion

        private void AddCurrentRecordToSearchResults(ref int currentScanResultsCount,
            udtInspectSearchResultType[] udtSearchResultsCurrentScan,
            udtInspectSearchResultType udtSearchResult)
        {
            if (currentScanResultsCount >= udtSearchResultsCurrentScan.Length)
            {
                Array.Resize(ref udtSearchResultsCurrentScan, udtSearchResultsCurrentScan.Length * 2);
            }

            udtSearchResultsCurrentScan[currentScanResultsCount] = udtSearchResult;
            currentScanResultsCount += 1;
        }

        private void AddDynamicAndStaticResidueMods(clsSearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            // Step through .PeptideSequenceWithMods
            // For each residue, check if a static mod is defined that affects that residue
            // For each mod symbol, determine the modification and add to searchResult

            var mostRecentLetter = '-';
            var residueLocInPeptide = 0;

            var sequence = searchResult.PeptideSequenceWithMods;
            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var chChar = sequence[index];

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
                                searchResult.SearchResultAddModification(
                                    modificationDefinition, chChar, residueLocInPeptide,
                                    searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(mostRecentLetter))
                {
                    searchResult.SearchResultAddDynamicModification(chChar, mostRecentLetter, residueLocInPeptide, searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
                }
                else
                {
                    // We found a modification symbol but mostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }
        }

        private readonly InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc sortScanChargeFScore = new InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc();
        private readonly InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc sortScanChargeMQScore = new InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc();
        private readonly InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc sortScanChargeTotalPRMDesc = new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc();

        /// <summary>
        /// Sorts the data by descending TotalPRMScore, than ranks each entry; in addition, computes normalized delta score (DeltaNorm) values
        /// </summary>
        /// <param name="udtSearchResultsCurrentScan"></param>
        /// <param name="currentScanResultsCount"></param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(ref udtInspectSearchResultType[] udtSearchResultsCurrentScan, int currentScanResultsCount)
        {
            const float DeltaNormMQScore_If_Undefined = 0;
            const float DeltaNormTotalPRMScore_If_Undefined = 0;

            var lastCharge = 0;
            double lastValue = 0;

            var currentRank = 0;

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, and descending RankFScore
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, currentScanResultsCount, sortScanChargeFScore);

            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index == 0 || udtSearchResultsCurrentScan[index].ChargeNum != lastCharge)
                {
                    lastCharge = udtSearchResultsCurrentScan[index].ChargeNum;
                    lastValue = udtSearchResultsCurrentScan[index].FScoreNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(udtSearchResultsCurrentScan[index].FScoreNum - lastValue) > float.Epsilon)
                    {
                        lastValue = udtSearchResultsCurrentScan[index].FScoreNum;
                        currentRank += 1;
                    }
                }

                udtSearchResultsCurrentScan[index].RankFScore = currentRank;
            }

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, and descending MQScore (note that MQScore can be negative)
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, currentScanResultsCount, sortScanChargeMQScore);

            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index < currentScanResultsCount - 1 && udtSearchResultsCurrentScan[index].ChargeNum == udtSearchResultsCurrentScan[index + 1].ChargeNum)
                {
                    udtSearchResultsCurrentScan[index].DeltaNormMQScore = ComputeDeltaNormScore(udtSearchResultsCurrentScan[index].MQScoreNum, udtSearchResultsCurrentScan[index + 1].MQScoreNum, DeltaNormMQScore_If_Undefined);
                }
                else
                {
                    udtSearchResultsCurrentScan[index].DeltaNormMQScore = 0;
                }
            }

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending PValue
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, currentScanResultsCount, sortScanChargeTotalPRMDesc);

            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index == 0 || udtSearchResultsCurrentScan[index].ChargeNum != lastCharge)
                {
                    lastCharge = udtSearchResultsCurrentScan[index].ChargeNum;
                    lastValue = udtSearchResultsCurrentScan[index].TotalPRMScoreNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(udtSearchResultsCurrentScan[index].TotalPRMScoreNum - lastValue) > float.Epsilon)
                    {
                        lastValue = udtSearchResultsCurrentScan[index].TotalPRMScoreNum;
                        currentRank += 1;
                    }
                }

                udtSearchResultsCurrentScan[index].RankTotalPRMScore = currentRank;

                if (index < currentScanResultsCount - 1 && udtSearchResultsCurrentScan[index].ChargeNum == udtSearchResultsCurrentScan[index + 1].ChargeNum)
                {
                    udtSearchResultsCurrentScan[index].DeltaNormTotalPRMScore = ComputeDeltaNormScore(udtSearchResultsCurrentScan[index].TotalPRMScoreNum, udtSearchResultsCurrentScan[index + 1].TotalPRMScoreNum, DeltaNormTotalPRMScore_If_Undefined);
                }
                else
                {
                    udtSearchResultsCurrentScan[index].DeltaNormTotalPRMScore = 0;
                }
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

        /// <summary>
        /// This routine creates a first hits file or synopsis file from the output from InSpecT
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="udtInspectModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
        /// <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreateFHTorSYNResultsFile(
            string inputFilePath,
            string outputFilePath,
            ref udtModInfoType[] udtInspectModInfo,
            eFilteredOutputFileTypeConstants eFilteredOutputFileType)
        {

            var udtSearchResult = new udtInspectSearchResultType();

            var resultID = 0;

            bool success;

            var errorLog = string.Empty;

            try
            {
                // Initialize variables
                var previousScan = Int32.MinValue;
                IComparer<udtInspectSearchResultType> sortComparer;

                if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                {
                    // Writes the synopsis file, which writes every record with a p-value below a set threshold or a TotalPRMScore above a certain threshold
                    sortComparer = new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc();
                }
                else
                {
                    if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.FHTbyTotalPRM)
                    {
                        // Write the PRM first-hits file, which writes the record with the highest TotalPRMScore
                        sortComparer = new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc();
                    }
                    else
                    {
                        // eFilteredOutputFileTypeConstants.FHTbyFScore
                        // Write the FScore first-hits file, which writes the record with the highest FScore
                        sortComparer = new InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc();
                    }
                }

                try
                {
                    // Open the input file and parse it
                    // Initialize the stream reader and the stream Text writer
                    using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
                        using (var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                        {
                            // Write the header line
                            WriteSynFHTFileHeader(writer, ref errorLog);

                            errorLog = string.Empty;
                            var resultsProcessed = 0;

                            // Initialize array that will hold all of the records for a given scan
                            var currentScanResultsCount = 0;
                            var udtSearchResultsCurrentScan = new udtInspectSearchResultType[10];

                            // Initialize the array that will hold all of the records that will ultimately be written out to disk
                            var filteredSearchResultCount = 0;
                            var udtFilteredSearchResults = new udtInspectSearchResultType[1000];

                            // Parse the input file
                            while (!reader.EndOfStream & !AbortProcessing)
                            {
                                var lineIn = reader.ReadLine();

                                if (string.IsNullOrWhiteSpace(lineIn))
                                {
                                    continue;
                                }

                                // Initialize udtSearchResult
                                udtSearchResult.Clear();

                                var validSearchResult = ParseInspectResultsFileEntry(ref lineIn, ref udtInspectModInfo, ref udtSearchResult, ref errorLog, resultsProcessed);

                                if (validSearchResult)
                                {
                                    if (previousScan != Int32.MinValue && previousScan != udtSearchResult.ScanNum)
                                    {
                                        // New scan encountered; sort and filter the data in udtSearchResultsCurrentScan, then call StoreTopFHTMatch or StoreSynMatches
                                        if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                                        {
                                            StoreSynMatches(writer, ref resultID, currentScanResultsCount, udtSearchResultsCurrentScan, ref filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog, ref sortComparer);
                                        }
                                        else
                                        {
                                            StoreTopFHTMatch(writer, ref resultID, currentScanResultsCount, udtSearchResultsCurrentScan, ref filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog, ref sortComparer);
                                        }
                                        currentScanResultsCount = 0;
                                    }

                                    AddCurrentRecordToSearchResults(ref currentScanResultsCount, udtSearchResultsCurrentScan, udtSearchResult);

                                    previousScan = udtSearchResult.ScanNum;
                                }

                                // Update the progress
                                UpdateProgress(Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100));

                                resultsProcessed += 1;
                            }

                            // Store the last record
                            if (currentScanResultsCount > 0)
                            {
                                if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                                {
                                    StoreSynMatches(writer, ref resultID, currentScanResultsCount, udtSearchResultsCurrentScan, ref filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog, ref sortComparer);
                                }
                                else
                                {
                                    StoreTopFHTMatch(writer, ref resultID, currentScanResultsCount, udtSearchResultsCurrentScan, ref filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog, ref sortComparer);
                                }
                                currentScanResultsCount = 0;
                            }

                            if (SortFHTAndSynFiles)
                            {
                                // Sort the data in udtFilteredSearchResults then write out to disk
                                SortAndWriteFilteredSearchResults(writer, filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog);
                            }
                        }
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
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        protected double ComputeDelMCorrectedPPM(
            double precursorErrorDa,
            double precursorMonoMass,
            double peptideMonoisotopicMass,
            bool adjustPrecursorMassForC13)
        {
            var peptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, adjustPrecursorMassForC13, peptideMonoisotopicMass);

            return peptideDeltaMassCorrectedPpm;
        }

        private float ComputeDeltaNormScore(float currentScore, float nextScore, float valueIfCurrentScoreZero)
        {
            try
            {
                if (Math.Abs(currentScore) > float.Epsilon)
                {
                    return Math.Abs((currentScore - nextScore) / currentScore);
                }

                return valueIfCurrentScoreZero;
            }
            catch (Exception)
            {
                return valueIfCurrentScoreZero;
            }
        }

        private double ComputePeptideMHFromPrecursorInfo(string precursorMZText, string precursorErrorText, string chargeText)
        {
            // Compute the theoretical peptide MH using the precursor m/z value and the precursor error values

            double peptideMH = 0;

            if (string.IsNullOrWhiteSpace(precursorMZText) || precursorMZText == "0")
            {
                // Precursor m/z is undefined; cannot continue
                peptideMH = 0;
            }
            else if (string.IsNullOrWhiteSpace(precursorErrorText))
            {
                // Precursor error is undefined; cannot continue
                peptideMH = 0;
            }
            else
            {
                var charge = CIntSafe(chargeText, -1);

                if (charge >= 1)
                {
                    if (double.TryParse(precursorMZText, out var precursorMZ))
                    {
                        if (double.TryParse(precursorErrorText, out var precursorError))
                        {
                            // Note: the October 2008 version of Inspect uses an Absolute Value function when computing the PrecursorError; the version used by PNNL does not use Absolute Value
                            // Note: switched to compute (M+H)+ in August 2011; prior to this, we were computing uncharged monoisotopic mass
                            peptideMH = (precursorMZ - precursorError) * charge - (charge - 1) * clsPeptideMassCalculator.MASS_PROTON;
                        }
                    }
                }
            }

            return peptideMH;
        }

        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputFolderPath, bool mts)
        {
            var pepToProteinMapFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
            if (pepToProteinMapFilePath.EndsWith("_inspect_syn", StringComparison.InvariantCultureIgnoreCase) ||
                pepToProteinMapFilePath.EndsWith("_inspect_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                pepToProteinMapFilePath = pepToProteinMapFilePath.Substring(0, pepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(pepToProteinMapFilePath, outputFolderPath, mts);
        }

        private bool ExtractModInfoFromInspectParamFile(string inspectParameterFilePath, ref udtModInfoType[] udtModList)
        {

            try
            {
                // Initialize udtModList and unnamedModID
                var modCount = 0;
                udtModList = new udtModInfoType[0];

                var unnamedModID = 0;

                if (string.IsNullOrWhiteSpace(inspectParameterFilePath))
                {
                    SetErrorMessage("Inspect Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                // Read the contents of the inspect parameter file
                using (var srInFile = new StreamReader(new FileStream(inspectParameterFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        var lineIn = srInFile.ReadLine();
                        if (string.IsNullOrWhiteSpace(lineIn))
                            continue;

                        var dataLine = lineIn.Trim();
                        if (dataLine.Length <= 0)
                            continue;

                        if (dataLine[0] == '#')
                        {
                            // Comment line; skip it
                            continue;
                        }

                        if (!dataLine.StartsWith("mod", StringComparison.InvariantCultureIgnoreCase))
                            continue;

                        // Modification definition line

                        // Split the line on commas
                        var splitLine = dataLine.Split(',');

                        if (splitLine.Length < 3 || splitLine[0].ToLower().Trim() != "mod")
                            continue;

                        if (udtModList == null || udtModList.Length == 0)
                        {
                            udtModList = new udtModInfoType[1];
                        }
                        else if (modCount >= udtModList.Length)
                        {
                            Array.Resize(ref udtModList, udtModList.Length * 2);
                        }

                        var udtMod = new udtModInfoType
                        {
                            ModMass = splitLine[1],
                            Residues = splitLine[2],
                            ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL.ToString()
                        };

                        if (splitLine.Length >= 4)
                        {
                            switch (splitLine[3].ToLower())
                            {
                                case "opt":
                                    udtMod.ModType = eInspectModType.DynamicMod;
                                    break;
                                case "fix":
                                    udtMod.ModType = eInspectModType.StaticMod;
                                    break;
                                case "nterminal":
                                    udtMod.ModType = eInspectModType.DynNTermPeptide;
                                    break;
                                case "cterminal":
                                    udtMod.ModType = eInspectModType.DynCTermPeptide;
                                    break;
                                default:
                                    ReportWarning("Unrecognized Mod Type in the Inspect parameter file");
                                    udtMod.ModType = eInspectModType.DynamicMod;
                                    break;
                            }
                        }
                        else
                        {
                            // Assume dynamic if not specified
                            udtMod.ModType = eInspectModType.DynamicMod;
                        }

                        if (splitLine.Length >= 5)
                        {
                            udtMod.ModName = splitLine[4].ToLower();
                            if (udtMod.ModName.Length > 4)
                            {
                                // Only keep the first 4 characters of the modification name
                                udtMod.ModName = udtMod.ModName.Substring(0, 4);
                            }
                        }
                        else
                        {
                            unnamedModID += 1;
                            udtMod.ModName = "UnnamedMod" + unnamedModID.ToString();
                        }

                        // Check for phosphorylation
                        // Inspect requires that it be defined in the parameter file as: mod,80,STY,opt,phosphorylation
                        //  However, we want to use the more precise mass of 79.9663
                        if (udtMod.ModName == PHOS_MOD_NAME.ToLower() & udtMod.ModMass == "80")
                        {
                            udtMod.ModMass = PHOS_MOD_MASS;
                        }
                        udtModList[modCount] = udtMod;

                        modCount += 1;
                    }
                }

                // Shrink udtModList to the appropriate length
                Array.Resize(ref udtModList, modCount);

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the Inspect parameter file (" + Path.GetFileName(inspectParameterFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }
        }

        private static readonly Regex RegexScanNumberRegEx = new Regex(DTA_FILENAME_SCAN_NUMBER_REGEX, REGEX_OPTIONS);

        private string ExtractScanNumFromDTAName(string spectrumFile)
        {
            var scanNum = string.Empty;

            // See if value resembles a .Dta file name
            // For example, "MyDataset.300.300.2.dta"

            try
            {
                var match = RegexScanNumberRegEx.Match(spectrumFile);
                if (match.Success && match.Groups.Count > 1)
                {
                    scanNum = match.Groups[1].Value;
                }
            }
            catch (Exception)
            {
                // Ignore errors here
                scanNum = "0";
            }

            return scanNum;
        }

        private void InitializeLocalVariables()
        {
            SortFHTAndSynFiles = true;
        }

        /// <summary>
        /// Load the PeptideToProteinMap information; in addition, creates the _inspect_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
        /// </summary>
        /// <param name="pepToProteinMapFilePath"></param>
        /// <param name="outputFolderPath"></param>
        /// <param name="udtInspectModInfo"></param>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool LoadPeptideToProteinMapInfoInspect(
            string pepToProteinMapFilePath,
            string outputFolderPath,
            ref udtModInfoType[] udtInspectModInfo,
            ref List<udtPepToProteinMappingType> pepToProteinMapping,
            ref string mtsPepToProteinMapFilePath)
        {
            bool success;

            try
            {
                mtsPepToProteinMapFilePath = string.Empty;

                if (string.IsNullOrWhiteSpace(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file is not defined");
                    return false;
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file does not exist: " + pepToProteinMapFilePath);
                    return false;
                }

                // Initialize pepToProteinMapping
                pepToProteinMapping = new List<udtPepToProteinMappingType>();

                // Read the data in proteinToPeptideMappingFilePath
                success = LoadPeptideToProteinMapInfo(pepToProteinMapFilePath, pepToProteinMapping, out var headerLine);

                if (success)
                {
                    mtsPepToProteinMapFilePath = Path.Combine(outputFolderPath, Path.GetFileNameWithoutExtension(pepToProteinMapFilePath) + "MTS.txt");

                    using (var swOutFile = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                    {
                        if (!string.IsNullOrEmpty(headerLine))
                        {
                            // Header line
                            swOutFile.WriteLine(headerLine);
                        }

                        for (var index = 0; index <= pepToProteinMapping.Count - 1; index++)
                        {
                            // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                            // In addition, replace the * terminus symbols with dashes
                            var mtsCompatiblePeptide = ReplaceInspectModTextWithSymbol(ReplaceTerminus(pepToProteinMapping[index].Peptide), ref udtInspectModInfo);

                            if (pepToProteinMapping[index].Peptide != mtsCompatiblePeptide)
                            {
                                UpdatePepToProteinMapPeptide(pepToProteinMapping, index, mtsCompatiblePeptide);
                            }

                            swOutFile.WriteLine(pepToProteinMapping[index].Peptide + "\t" +
                                                pepToProteinMapping[index].Protein + "\t" +
                                                pepToProteinMapping[index].ResidueStart + "\t" +
                                                pepToProteinMapping[index].ResidueEnd);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" + Path.GetFileName(mtsPepToProteinMapFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        private bool ParseInspectSynFileHeaderLine(string lineIn, out int[] columnMapping)
        {
            // Parse the header line

            var columnNames = new SortedDictionary<string, eInspectSynFileColumns>(StringComparer.InvariantCultureIgnoreCase)
            {
                {"ResultID", eInspectSynFileColumns.ResultID},
                {"Scan", eInspectSynFileColumns.Scan},
                {"Peptide", eInspectSynFileColumns.Peptide},
                {"Protein", eInspectSynFileColumns.Protein},
                {"Charge", eInspectSynFileColumns.Charge},
                {"MQScore", eInspectSynFileColumns.MQScore},
                {"Length", eInspectSynFileColumns.Length},
                {"TotalPRMScore", eInspectSynFileColumns.TotalPRMScore},
                {"MedianPRMScore", eInspectSynFileColumns.MedianPRMScore},
                {"FractionY", eInspectSynFileColumns.FractionY},
                {"FractionB", eInspectSynFileColumns.FractionB},
                {"Intensity", eInspectSynFileColumns.Intensity},
                {"NTT", eInspectSynFileColumns.NTT},
                {"PValue", eInspectSynFileColumns.PValue},
                {"FScore", eInspectSynFileColumns.FScore},
                {"DeltaScore", eInspectSynFileColumns.DeltaScore},
                {"DeltaScoreOther", eInspectSynFileColumns.DeltaScoreOther},
                {"DeltaNormMQScore", eInspectSynFileColumns.DeltaNormMQScore},
                {"DeltaNormTotalPRMScore", eInspectSynFileColumns.DeltaNormTotalPRMScore},
                {"RankTotalPRMScore", eInspectSynFileColumns.RankTotalPRMScore},
                {"RankFScore", eInspectSynFileColumns.RankFScore},
                {"MH", eInspectSynFileColumns.MH},
                {"RecordNumber", eInspectSynFileColumns.RecordNumber},
                {"DBFilePos", eInspectSynFileColumns.DBFilePos},
                {"SpecFilePos", eInspectSynFileColumns.SpecFilePos},
                {"PrecursorMZ", eInspectSynFileColumns.PrecursorMZ},
                {"PrecursorError", eInspectSynFileColumns.PrecursorError}
            };

            columnMapping = new int[InspectSynopsisFileColCount];

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
                SetErrorMessage("Error parsing header in Inspect synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        protected bool ParseInspectSynopsisFile(string inputFilePath, string outputFolderPath, ref List<udtPepToProteinMappingType> pepToProteinMapping, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that Inspect synopsis files are normally sorted on TotalPRMScore descending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique TotalPRMScore encountered

            var currentPeptideWithMods = string.Empty;

            int[] columnMapping = null;

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
                var searchResult = new clsSearchResultsInSpecT(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForTotalPRMScoreLevel
                var htPeptidesFoundForTotalPRMScoreLevel = new Hashtable();
                var previousTotalPRMScore = string.Empty;

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
                                success = ParseInspectSynFileHeaderLine(lineIn, out columnMapping);
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
                                validSearchResult = ParseInspectSynFileEntry(lineIn, columnMapping, searchResult, ref errorLog, out currentPeptideWithMods);
                            }
                            else
                            {
                                validSearchResult = false;
                            }

                            if (validSearchResult)
                            {
                                var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;
                                bool firstMatchForGroup;

                                if (searchResult.TotalPRMScore == previousTotalPRMScore)
                                {
                                    // New result has the same TotalPRMScore as the previous result
                                    // See if htPeptidesFoundForTotalPRMScoreLevel contains the peptide, scan and charge

                                    if (htPeptidesFoundForTotalPRMScoreLevel.ContainsKey(key))
                                    {
                                        firstMatchForGroup = false;
                                    }
                                    else
                                    {
                                        htPeptidesFoundForTotalPRMScoreLevel.Add(key, 1);
                                        firstMatchForGroup = true;
                                    }
                                }
                                else
                                {
                                    // New TotalPRMScore
                                    // Reset htPeptidesFoundForTotalPRMScoreLevel
                                    htPeptidesFoundForTotalPRMScoreLevel.Clear();

                                    // Update previousTotalPRMScore
                                    previousTotalPRMScore = searchResult.TotalPRMScore;

                                    // Append a new entry to htPeptidesFoundForTotalPRMScoreLevel
                                    htPeptidesFoundForTotalPRMScoreLevel.Add(key, 1);
                                    firstMatchForGroup = true;
                                }

                                success = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                                if (!success)
                                {
                                    if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        errorLog += "Error adding modifications to sequence at RowIndex '" + searchResult.ResultID + "'" + "\n";
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

        private bool ParseInspectResultsFileEntry(
            ref string lineIn,
            ref udtModInfoType[] udtInspectModInfo,
            ref udtInspectSearchResultType udtSearchResult,
            ref string errorLog,
            int resultsProcessed)
        {
            // Parses an entry from the Inspect results file
            // The expected header line is:
            // #SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos PrecursorMZ	PrecursorError DelM_PPM

            string[] splitLine = null;

            bool validSearchResult;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                // Reset udtSearchResult
                udtSearchResult.Clear();

                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length >= 15)
                {
                    if (resultsProcessed == 0)
                    {
                        // This is the first line of the file; it may be a header row
                        // Determine this by seeing if any of the first three columns contains a number
                        if (!(clsPHRPParser.IsNumber(splitLine[0]) ||
                              clsPHRPParser.IsNumber(splitLine[1]) ||
                              clsPHRPParser.IsNumber(splitLine[2])))
                        {
                            // This is a header line; ignore it
                            validSearchResult = false;
                            return false;
                        }
                    }

                    udtSearchResult.SpectrumFileName = splitLine[(int)eInspectResultsFileColumns.SpectrumFile];
                    if (splitLine[(int)eInspectResultsFileColumns.Scan] == "0")
                    {
                        udtSearchResult.Scan = ExtractScanNumFromDTAName(udtSearchResult.SpectrumFileName);
                    }
                    else
                    {
                        udtSearchResult.Scan = splitLine[(int)eInspectResultsFileColumns.Scan];
                    }
                    udtSearchResult.ScanNum = CIntSafe(udtSearchResult.Scan, 0);

                    // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                    // In addition, replace the * terminus symbols with dashes
                    udtSearchResult.PeptideAnnotation = ReplaceInspectModTextWithSymbol(ReplaceTerminus(splitLine[(int)eInspectResultsFileColumns.Annotation]), ref udtInspectModInfo);
                    udtSearchResult.Protein = TruncateProteinName(splitLine[(int)eInspectResultsFileColumns.Protein]);

                    udtSearchResult.Charge = splitLine[(int)eInspectResultsFileColumns.Charge];
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    udtSearchResult.MQScore = splitLine[(int)eInspectResultsFileColumns.MQScore];
                    udtSearchResult.MQScoreNum = CSngSafe(udtSearchResult.MQScore, 0);

                    udtSearchResult.Length = CIntSafe(splitLine[(int)eInspectResultsFileColumns.Length], 0);

                    udtSearchResult.TotalPRMScore = splitLine[(int)eInspectResultsFileColumns.TotalPRMScore];
                    udtSearchResult.TotalPRMScoreNum = CSngSafe(udtSearchResult.TotalPRMScore, 0);

                    udtSearchResult.MedianPRMScore = splitLine[(int)eInspectResultsFileColumns.MedianPRMScore];
                    udtSearchResult.FractionY = RemoveExtraneousDigits(splitLine[(int)eInspectResultsFileColumns.FractionY]);
                    udtSearchResult.FractionB = RemoveExtraneousDigits(splitLine[(int)eInspectResultsFileColumns.FractionB]);
                    udtSearchResult.Intensity = splitLine[(int)eInspectResultsFileColumns.Intensity];
                    udtSearchResult.NTT = CIntSafe(splitLine[(int)eInspectResultsFileColumns.NTT], 0);

                    udtSearchResult.pValue = RemoveExtraneousDigits(splitLine[(int)eInspectResultsFileColumns.PValue]);
                    udtSearchResult.PValueNum = CSngSafe(udtSearchResult.pValue, 0);

                    udtSearchResult.FScore = splitLine[(int)eInspectResultsFileColumns.FScore];
                    udtSearchResult.FScoreNum = CSngSafe(udtSearchResult.FScore, 0);

                    udtSearchResult.DeltaScore = splitLine[(int)eInspectResultsFileColumns.DeltaScore];
                    udtSearchResult.DeltaScoreOther = splitLine[(int)eInspectResultsFileColumns.DeltaScoreOther];

                    udtSearchResult.RecordNumber = splitLine[(int)eInspectResultsFileColumns.RecordNumber];
                    udtSearchResult.DBFilePos = splitLine[(int)eInspectResultsFileColumns.DBFilePos];
                    udtSearchResult.SpecFilePos = splitLine[(int)eInspectResultsFileColumns.SpecFilePos];

                    if (splitLine.Length >= (int)eInspectResultsFileColumns.PrecursorError + 1)
                    {
                        // Inspect version 2008-10-14 added these two Precursor mass columns
                        udtSearchResult.PrecursorMZ = splitLine[(int)eInspectResultsFileColumns.PrecursorMZ];
                        udtSearchResult.PrecursorError = splitLine[(int)eInspectResultsFileColumns.PrecursorError];

                        udtSearchResult.MH = ComputePeptideMHFromPrecursorInfo(udtSearchResult.PrecursorMZ, udtSearchResult.PrecursorError, udtSearchResult.Charge);

                        if (double.TryParse(udtSearchResult.PrecursorMZ, out var precursorMZ))
                        {
                            var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, udtSearchResult.ChargeNum, 0);
                            var peptideMonoisotopicMass = udtSearchResult.MH - clsPeptideMassCalculator.MASS_PROTON;

                            var precursorErrorDa = precursorMonoMass - peptideMonoisotopicMass;

                            var peptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, peptideMonoisotopicMass, true);

                            udtSearchResult.DelMPPM = PRISM.StringUtilities.DblToString(peptideDeltaMassCorrectedPpm, 5, 0.00005);
                        }
                    }
                    else
                    {
                        udtSearchResult.PrecursorMZ = "0";
                        udtSearchResult.PrecursorError = "0";
                        udtSearchResult.MH = 0;
                        udtSearchResult.DelMPPM = "0";
                    }

                    validSearchResult = true;
                }
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing InSpecT Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing InSpecT Results in ParseInspectResultsFileEntry" + "\n";
                    }
                }
                validSearchResult = false;
            }

            return validSearchResult;
        }

        private bool ParseInspectSynFileEntry(
            string lineIn,
            IReadOnlyList<int> columnMapping,
            clsSearchResultsInSpecT searchResult,
            ref string errorLog,
            out string peptideSequenceWithMods)
        {
            // Parses an entry from the Inspect Synopsis file

            string[] splitLine = null;

            bool validSearchResult;

            peptideSequenceWithMods = string.Empty;

            try
            {
                // Set this to False for now
                validSearchResult = false;

                // Reset searchResult
                searchResult.Clear();

                splitLine = lineIn.TrimEnd().Split('\t');
                if (splitLine.Length < INSPECT_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.ResultID], out int resultId))
                {
                    ReportError("ResultID column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.Charge], out string charge);

                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.Protein], out string proteinName);

                searchResult.ResultID = resultId;
                searchResult.Scan = scan;
                searchResult.Charge = charge;
                searchResult.ProteinName = proteinName;

                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.Peptide], out peptideSequenceWithMods);

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (clsSearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                if (GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.PrecursorError], out string peptideDeltaMass))
                {
                    searchResult.PeptideDeltaMass = peptideDeltaMass;
                    // Note: .peptideDeltaMass is stored in the Inspect results file as "Observed_Mass - Theoretical_Mass"
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
                }
                else
                {
                    searchResult.PeptideDeltaMass = "0";
                }

                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.MQScore], out string mqScore);

                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.Length], out string length);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.TotalPRMScore], out string totalPrmScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.MedianPRMScore], out string medianPrmScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.FractionY], out string fractionY);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.FractionB], out string fractionB);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.Intensity], out string intensity);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.NTT], out string ntt);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.PValue], out string pValue);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.FScore], out string fScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.DeltaScore], out string deltaScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.DeltaScoreOther], out string deltaScoreOther);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.DeltaNormMQScore], out string deltaNormMqScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.DeltaNormTotalPRMScore], out string deltaNormTotalPrmScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.RankTotalPRMScore], out string rankTotalPrmScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.RankFScore], out string rankFScore);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.MH], out string peptideMh);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.RecordNumber], out string recordNumber);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.DBFilePos], out string dbFilePos);
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.SpecFilePos], out string specFilePos);

                // Note: .PrecursorError was processed earlier in this function
                GetColumnValue(splitLine, columnMapping[(int)eInspectSynFileColumns.PrecursorMZ], out string precursorMz);

                searchResult.MQScore = mqScore;
                searchResult.Length = length;
                searchResult.TotalPRMScore = totalPrmScore;
                searchResult.MedianPRMScore = medianPrmScore;
                searchResult.FractionY = fractionY;
                searchResult.FractionB = fractionB;
                searchResult.Intensity = intensity;
                searchResult.NTT = ntt;
                searchResult.pValue = pValue;
                searchResult.FScore = fScore;
                searchResult.DeltaScore = deltaScore;
                searchResult.DeltaScoreOther = deltaScoreOther;
                searchResult.DeltaNormMQScore = deltaNormMqScore;
                searchResult.DeltaNormTotalPRMScore = deltaNormTotalPrmScore;
                searchResult.RankTotalPRMScore = rankTotalPrmScore;
                searchResult.RankFScore = rankFScore;
                searchResult.PeptideMH = peptideMh;
                searchResult.RecordNumber = recordNumber;
                searchResult.DBFilePos = dbFilePos;
                searchResult.SpecFilePos = specFilePos;
                searchResult.PrecursorMz = precursorMz;

                validSearchResult = true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if (splitLine != null && splitLine.Length > 0)
                    {
                        errorLog += "Error parsing InSpecT Results for RowIndex '" + splitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        errorLog += "Error parsing InSpecT Results in ParseInspectSynFileEntry" + "\n";
                    }
                }
                validSearchResult = false;
            }

            return validSearchResult;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">Inspect results file</param>
        /// <param name="outputFolderPath">Output folder</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputFolderPath, string parameterFilePath)
        {
            var mtsPepToProteinMapFilePath = string.Empty;

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

                    var udtInspectModInfo = new udtModInfoType[0];
                    var pepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the Inspect Parameter File so that we can determine the modification names and masses
                    if (!ExtractModInfoFromInspectParamFile(SearchToolParameterFilePath, ref udtInspectModInfo))
                    {
                        if (udtInspectModInfo == null || udtInspectModInfo.Length == 0)
                        {
                            udtInspectModInfo = new udtModInfoType[1];
                            var modInfo = new udtModInfoType
                            {
                                ModName = PHOS_MOD_NAME.ToLower(),
                                ModMass = PHOS_MOD_MASS,
                                Residues = PHOS_MOD_RESIDUES,
                                ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL.ToString()
                            };
                            udtInspectModInfo[1] = modInfo;
                        }
                    }

                    // Resolve the mods in mInspectModInfo with the ModDefs mods
                    ResolveInspectModsWithModDefinitions(ref udtInspectModInfo);

                    if (CreateInspectFirstHitsFile)
                    {
                        // Create the first hits output file
                        ResetProgress("Creating the FHT file (top TotalPRMScore)", true);

                        var outputFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
                        outputFilePath = Path.Combine(outputFolderPath, outputFilePath + INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX);

                        success = CreateFHTorSYNResultsFile(inputFilePath, outputFilePath, ref udtInspectModInfo, eFilteredOutputFileTypeConstants.FHTbyTotalPRM);

                        // Create the first hits output file
                        ResetProgress("Creating the FHT file (top FScore)", true);

                        outputFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
                        outputFilePath = Path.Combine(outputFolderPath, outputFilePath + INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX);

                        success = CreateFHTorSYNResultsFile(inputFilePath, outputFilePath, ref udtInspectModInfo, eFilteredOutputFileTypeConstants.FHTbyFScore);
                    }

                    if (CreateInspectSynopsisFile)
                    {
                        // Create the synopsis output file
                        ResetProgress("Creating the SYN file", true);

                        // Define the synopsis output file name based on inputFilePath
                        var synOutputFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
                        synOutputFilePath = Path.Combine(outputFolderPath, synOutputFilePath + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                        success = CreateFHTorSYNResultsFile(inputFilePath, synOutputFilePath, ref udtInspectModInfo, eFilteredOutputFileTypeConstants.SynFile);

                        // Load the PeptideToProteinMap information; if the file doesn't exist, then a warning will be displayed, but processing will continue
                        // LoadPeptideToProteinMapInfoInspect also creates _inspect_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
                        var pepToProteinMapFilePath = Path.Combine(Path.GetDirectoryName(inputFile.FullName), Path.GetFileNameWithoutExtension(inputFile.FullName) + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + ".txt");

                        ResetProgress("Loading the PepToProtein map file: " + Path.GetFileName(pepToProteinMapFilePath), true);

                        LoadPeptideToProteinMapInfoInspect(pepToProteinMapFilePath, outputFolderPath, ref udtInspectModInfo, ref pepToProteinMapping, ref mtsPepToProteinMapFilePath);

                        // Create the other PHRP-specific files
                        ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                        success = ParseInspectSynopsisFile(synOutputFilePath, outputFolderPath, ref pepToProteinMapping, false);

                        // Remove all items from pepToProteinMapping to reduce memory overhead
                        pepToProteinMapping.Clear();
                        pepToProteinMapping.TrimExcess();

                        if (success && CreateProteinModsFile)
                        {
                            // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                            ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.DirectoryName, Path.GetFileName(synOutputFilePath)), outputFolderPath);

                            // Create the Protein Mods file
                            success = CreateProteinModDetailsFile(synOutputFilePath, outputFolderPath, mtsPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Inspect);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error calling CreateFHTorSYNResultsFile: " + ex.Message, ex);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile:" + ex.Message, ex);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return success;
        }

        private static readonly Regex RegexNumPlusZeroes = new Regex(@"(\.\d*[1-9])0+$", RegexOptions.Compiled);
        private static readonly Regex RegexAllZeroes = new Regex(@"\.0+$", RegexOptions.Compiled);

        // If value ends in .0000, then remove the .0000 portion
        private string RemoveExtraneousDigits(string value)
        {
            if (string.IsNullOrWhiteSpace(value))
            {
                return string.Empty;
            }

            var reMatch = RegexAllZeroes.Match(value);
            if (reMatch.Success && reMatch.Index > 0)
            {
                value = value.Substring(0, reMatch.Index);
            }
            else
            {
                reMatch = RegexNumPlusZeroes.Match(value);
                if (reMatch.Success && reMatch.Index > 0)
                {
                    if (reMatch.Groups.Count > 1)
                    {
                        // Number is of the form 1.0030 or 1.300 or 1.030
                        value = value.Substring(0, reMatch.Index) + reMatch.Groups[1].Value;
                    }
                }
            }

            return value;
        }

        private static readonly Regex RegexNTerminalModMassRegEx = new Regex(INSPECT_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);
        private static readonly Regex RegexCTerminalModMassRegEx = new Regex(INSPECT_CTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Replaces modification name text in peptide sequences with modification symbols (uses case-sensitive comparisons)
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="udtInspectModInfo">This function assumes that each entry in udtInspectModInfo() has both .ModName and .ModSymbol defined</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected string ReplaceInspectModTextWithSymbol(string peptide, ref udtModInfoType[] udtInspectModInfo)
        {
            var prefix = string.Empty;
            var suffix = string.Empty;

            if (peptide.Length >= 4)
            {
                if (peptide[1] == '.' && peptide[peptide.Length - 2] == '.')
                {
                    prefix = peptide.Substring(0, 2);
                    suffix = peptide.Substring(peptide.Length - 2, 2);

                    peptide = peptide.Substring(2, peptide.Length - 4);
                }
            }

            // peptide should now be the clean peptide, without the prefix or suffix residues
            for (var index = 0; index <= udtInspectModInfo.Length - 1; index++)
            {
                if (udtInspectModInfo[index].ModType != eInspectModType.StaticMod)
                {
                    peptide = peptide.Replace(udtInspectModInfo[index].ModName, udtInspectModInfo[index].ModSymbol);

                    if (udtInspectModInfo[index].ModType == eInspectModType.DynNTermPeptide |
                        udtInspectModInfo[index].ModType == eInspectModType.DynCTermPeptide)
                    {
                        var reMatch = default(Match);
                        if (udtInspectModInfo[index].ModType == eInspectModType.DynNTermPeptide)
                        {
                            // Inspect notates N-terminal mods like this: R.+14HVIFLAER.R   (Note: This behavior is not yet confirmed)
                            // Look for this using reNTerminalModMassRegEx
                            reMatch = RegexNTerminalModMassRegEx.Match(peptide);
                        }
                        else if (udtInspectModInfo[index].ModType == eInspectModType.DynCTermPeptide)
                        {
                            // Inspect notates C-terminal mods like this: R.HVIFLAER+14.R
                            // Look for this using reCTerminalModMassRegEx
                            reMatch = RegexCTerminalModMassRegEx.Match(peptide);
                        }
                        else
                        {
                            // This code should never be reached
                            reMatch = null;
                        }

                        if (reMatch != null)
                        {
                            if (reMatch.Success && reMatch.Groups.Count > 1)
                            {
                                // Match found
                                try
                                {
                                    var modMass = Convert.ToInt32(reMatch.Groups[1].Value);

                                    // Compare the mod mass in the specification to this Mod's mod mass
                                    // If they are less than 0.5 Da apart, then assume we have a match; yes, this assumption is a bit flaky
                                    if (Math.Abs(modMass - Convert.ToDouble(udtInspectModInfo[index].ModMass)) <= 0.5)
                                    {
                                        // Match found
                                        // Replace the matched region with .ModSymbol

                                        string peptideNew;

                                        if (reMatch.Groups[0].Index > 0)
                                        {
                                            peptideNew = peptide.Substring(0, reMatch.Groups[0].Index);
                                        }
                                        else
                                        {
                                            peptideNew = string.Empty;
                                        }

                                        peptideNew += udtInspectModInfo[index].ModSymbol;

                                        if (reMatch.Groups[0].Index + reMatch.Groups[0].Length < peptide.Length)
                                        {
                                            peptideNew += peptide.Substring(reMatch.Groups[0].Index + reMatch.Groups[0].Length);
                                        }

                                        peptide = string.Copy(peptideNew);
                                    }
                                }
                                catch (Exception ex)
                                {
                                    ReportError("Error comparing mod mass in peptide to mod mass in udtInspectModInfo", true, ex);
                                }
                            }
                        }
                    }
                }
            }

            return prefix + peptide + suffix;
        }

        protected string ReplaceTerminus(string peptide)
        {
            if (peptide.StartsWith(N_TERMINUS_SYMBOL_INSPECT))
            {
                peptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_INSPECT.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_INSPECT))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_INSPECT.Length) + "." + clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        protected void ResolveInspectModsWithModDefinitions(ref udtModInfoType[] udtInspectModInfo)
        {
            var existingModFound = false;

            if (udtInspectModInfo == null)
                return;

            // Call .LookupModificationDefinitionByMass for each entry in udtInspectModInfo
            for (var index = 0; index <= udtInspectModInfo.Length - 1; index++)
            {
                var modInfo = udtInspectModInfo[index];
                if (double.TryParse(modInfo.ModMass, out var modMass))
                {
                    int resIndexStart;
                    int resIndexEnd;

                    if (modInfo.Residues.Length > 0)
                    {
                        resIndexStart = 0;
                        resIndexEnd = modInfo.Residues.Length - 1;
                    }
                    else
                    {
                        resIndexStart = -1;
                        resIndexEnd = -1;
                    }

                    for (var residueIndex = resIndexStart; residueIndex <= resIndexEnd; residueIndex++)
                    {
                        char targetResidue;
                        if (residueIndex >= 0)
                        {
                            targetResidue = modInfo.Residues[residueIndex];
                        }
                        else
                        {
                            targetResidue = default(char);
                        }

                        clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState;
                        if (modInfo.ModType == eInspectModType.DynNTermPeptide)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                        }
                        else if (modInfo.ModType == eInspectModType.DynCTermPeptide)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                        }
                        else
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
                        }

                        var modificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(modMass, targetResidue, eResidueTerminusState, out existingModFound, true);

                        if (residueIndex == resIndexStart)
                        {
                            modInfo.ModSymbol = modificationDefinition.ModificationSymbol.ToString();
                        }
                    }
                }
                udtInspectModInfo[index] = modInfo;
            }
        }

        protected void StoreOrWriteSearchResult(
            StreamWriter writer,
            ref int resultID,
            udtInspectSearchResultType udtSearchResult,
            ref int filteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string errorLog)
        {
            if (SortFHTAndSynFiles)
            {
                if (filteredSearchResultCount == udtFilteredSearchResults.Length)
                {
                    Array.Resize(ref udtFilteredSearchResults, udtFilteredSearchResults.Length * 2);
                }

                udtFilteredSearchResults[filteredSearchResultCount] = udtSearchResult;
                filteredSearchResultCount += 1;
            }
            else
            {
                resultID += 1;
                WriteSearchResultToFile(resultID, writer, udtSearchResult, ref errorLog);
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            int filteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string errorLog)
        {
            // Sort udtFilteredSearchResults by descending TotalPRMScore, ascending scan, ascending charge, ascending peptide, and ascending protein
            Array.Sort(udtFilteredSearchResults, 0, filteredSearchResultCount, new InspectSearchResultsComparerTotalPRMDescScanChargePeptide());

            for (var index = 0; index <= filteredSearchResultCount - 1; index++)
            {
                WriteSearchResultToFile(index + 1, writer, udtFilteredSearchResults[index], ref errorLog);
            }
        }

        private void StoreTopFHTMatch(
            StreamWriter writer,
            ref int resultID,
            int currentScanResultsCount,
            udtInspectSearchResultType[] udtSearchResultsCurrentScan,
            ref int filteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string errorLog,
            ref IComparer<udtInspectSearchResultType> sortComparer)
        {
            var currentCharge = short.MinValue;

            AssignRankAndDeltaNormValues(ref udtSearchResultsCurrentScan, currentScanResultsCount);

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, then descending TotalPRMScore or descending FScore (depending on sortComparer)
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, currentScanResultsCount, sortComparer);

            // Now store or write out the first match for each charge for this scan
            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index == 0 || currentCharge != udtSearchResultsCurrentScan[index].ChargeNum)
                {
                    StoreOrWriteSearchResult(writer, ref resultID, udtSearchResultsCurrentScan[index], ref filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog);
                    currentCharge = udtSearchResultsCurrentScan[index].ChargeNum;
                }
            }
        }

        private void StoreSynMatches(
            StreamWriter writer,
            ref int resultID,
            int currentScanResultsCount,
            udtInspectSearchResultType[] udtSearchResultsCurrentScan,
            ref int filteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string errorLog,
            ref IComparer<udtInspectSearchResultType> sortComparer)
        {
            AssignRankAndDeltaNormValues(ref udtSearchResultsCurrentScan, currentScanResultsCount);

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending FScore
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, currentScanResultsCount, sortComparer);

            // Now store or write out the matches that pass the filters
            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (udtSearchResultsCurrentScan[index].PValueNum <= InspectSynopsisFilePValueThreshold ||
                    udtSearchResultsCurrentScan[index].TotalPRMScoreNum >= TOTALPRMSCORE_THRESHOLD ||
                    udtSearchResultsCurrentScan[index].FScoreNum >= FSCORE_THRESHOLD)
                {
                    StoreOrWriteSearchResult(writer, ref resultID, udtSearchResultsCurrentScan[index], ref filteredSearchResultCount, ref udtFilteredSearchResults, ref errorLog);
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
                    COLUMN_NAME_RESULTID,
                    "Scan",
                    COLUMN_NAME_PEPTIDE,
                    "Protein",
                    "Charge",
                    "MQScore",
                    "Length",
                    "TotalPRMScore",
                    "MedianPRMScore",
                    "FractionY",
                    "FractionB",
                    "Intensity",
                    "NTT",
                    "PValue",
                    "FScore",
                    "DeltaScore",
                    "DeltaScoreOther",
                    "DeltaNormMQScore",
                    "DeltaNormTotalPRMScore",
                    "RankTotalPRMScore",
                    "RankFScore",
                    "MH",
                    "RecordNumber",
                    "DBFilePos",
                    "SpecFilePos",
                    "PrecursorMZ",
                    "PrecursorError",
                    "DelM_PPM"
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

        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            udtInspectSearchResultType udtSearchResult,
            ref string errorLog)
        {
            // Writes an entry to a synopsis or first hits file
            try
            {
                writer.WriteLine(resultID + "\t" +
                                       udtSearchResult.Scan + "\t" +
                                       udtSearchResult.PeptideAnnotation + "\t" +
                                       udtSearchResult.Protein + "\t" +
                                       udtSearchResult.Charge + "\t" +
                                       udtSearchResult.MQScore + "\t" +
                                       udtSearchResult.Length + "\t" +
                                       udtSearchResult.TotalPRMScore + "\t" +
                                       udtSearchResult.MedianPRMScore + "\t" +
                                       udtSearchResult.FractionY + "\t" +
                                       udtSearchResult.FractionB + "\t" +
                                       udtSearchResult.Intensity + "\t" +
                                       udtSearchResult.NTT + "\t" +
                                       udtSearchResult.pValue + "\t" +
                                       udtSearchResult.FScore + "\t" +
                                       udtSearchResult.DeltaScore + "\t" +
                                       udtSearchResult.DeltaScoreOther + "\t" +
                                       PRISM.StringUtilities.DblToString(udtSearchResult.DeltaNormMQScore, 5) + "\t" +
                                       PRISM.StringUtilities.DblToString(udtSearchResult.DeltaNormTotalPRMScore, 5) + "\t" +
                                       udtSearchResult.RankTotalPRMScore + "\t" +
                                       udtSearchResult.RankFScore + "\t" +
                                       PRISM.StringUtilities.DblToString(udtSearchResult.MH, 6) + "\t" +
                                       udtSearchResult.RecordNumber + "\t" +
                                       udtSearchResult.DBFilePos + "\t" +
                                       udtSearchResult.SpecFilePos + "\t" +
                                       udtSearchResult.PrecursorMZ + "\t" +
                                       udtSearchResult.PrecursorError + "\t" +
                                       udtSearchResult.DelMPPM);
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

        protected class InspectSearchResultsComparerTotalPRMDescScanChargePeptide : IComparer<udtInspectSearchResultType>
        {
            public int Compare(udtInspectSearchResultType x, udtInspectSearchResultType y)
            {
                if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                {
                    return -1;
                }

                if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                {
                    return 1;
                }

                // TotalPRMScore is the same; check scan number
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
                var result = string.Compare(x.PeptideAnnotation, y.PeptideAnnotation, StringComparison.Ordinal);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.Compare(x.Protein, y.Protein, StringComparison.Ordinal);
                }
                return result;
            }
        }

        protected class InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc : IComparer<udtInspectSearchResultType>
        {
            public int Compare(udtInspectSearchResultType x, udtInspectSearchResultType y)
            {
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {
                    return -1;
                }

                // Charge is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                {
                    return -1;
                }

                if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                {
                    return 1;
                }

                // TotalPRMScore is the same; check FScore (sort on descending FScore)
                if (x.FScoreNum > y.FScoreNum)
                {
                    return -1;
                }

                if (x.FScoreNum < y.FScoreNum)
                {
                    return 1;
                }

                return 0;
            }
        }

        protected class InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc : IComparer<udtInspectSearchResultType>
        {
            public int Compare(udtInspectSearchResultType x, udtInspectSearchResultType y)
            {
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {

                    return -1;
                }
                // Charge is the same; check FScore (sort on descending FScore)
                if (x.FScoreNum > y.FScoreNum)
                {
                    return -1;
                }

                if (x.FScoreNum < y.FScoreNum)
                {
                    return 1;
                }

                // FScore is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                {
                    return -1;
                }

                if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                {
                    return 1;
                }

                return 0;
            }
        }

        protected class InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc : IComparer<udtInspectSearchResultType>
        {
            public int Compare(udtInspectSearchResultType x, udtInspectSearchResultType y)
            {
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {
                    return -1;
                }

                // Charge is the same; check MQScore (sort on descending MQScore)
                if (x.MQScoreNum > y.MQScoreNum)
                {
                    return -1;
                }

                if (x.MQScoreNum < y.MQScoreNum)
                {
                    return 1;
                }

                // MQScore is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                {
                    return -1;
                }

                if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                {
                    return 1;
                }

                return 0;
            }
        }

        #endregion
    }
}
