// -------------------------------------------------------------------------------
// Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
// Program started August 12, 2008
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads an InSpecT results file (txt format) and creates
    /// a tab-delimited text file with the data.  It will insert modification symbols
    /// into the peptide sequences for modified peptides.
    /// </summary>
    /// <remarks>The modification definition information is determined from the InSpecT parameter file</remarks>
    public class InSpecTResultsProcessor : PHRPBaseClass
    {
        // Ignore Spelling: cterminal, Da, Daltons, fht, methylation, ModDefs, nterminal, phos, phosphorylation, Pos, txt

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="options"></param>
        public InSpecTResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "April 17, 2019";
            InitializeLocalVariables();
        }

        /// <summary>
        /// Tool name
        /// </summary>
        public const string TOOL_NAME = "InSpecT";

        /// <summary>
        /// InSpecT results file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_INSPECT_FILE = "_inspect";

        private const int INSPECT_SYN_FILE_MIN_COL_COUNT = 5;

        /// <summary>
        /// N-terminus symbol used by InSpecT
        /// </summary>
        public const string N_TERMINUS_SYMBOL_INSPECT = "*.";

        /// <summary>
        /// C-terminus symbol used by InSpecT
        /// </summary>
        public const string C_TERMINUS_SYMBOL_INSPECT = ".*";

        private const char UNKNOWN_INSPECT_MOD_SYMBOL = '?';

        // When writing the synopsis file, we keep data that passes any of these thresholds (thus, it's an OR comparison, not an AND comparison)
        // pValue <= 0.2 Or TotalPRMScore >= 50 or FScore >= 0

        /// <summary>
        /// Default synopsis file p-value threshold
        /// </summary>
        public const float DEFAULT_SYN_FILE_PVALUE_THRESHOLD = 0.2f;

        /// <summary>
        /// Default synopsis file TotalPRMScore threshold
        /// </summary>
        public const float TOTALPRMSCORE_THRESHOLD = 50;

        /// <summary>
        /// Default synopsis file FScore threshold
        /// </summary>
        public const float FSCORE_THRESHOLD = 0;

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        private const string DTA_FILENAME_SCAN_NUMBER_REGEX = @"(\d+)\.\d+\.\d+\.dta";
        private const string INSPECT_NTERMINAL_MOD_MASS_REGEX = @"^\+(\d+)";
        private const string INSPECT_CTERMINAL_MOD_MASS_REGEX = @"\+(\d+)$";

        private const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

        private const string PHOS_MOD_NAME = "phos";
        private const string PHOS_MOD_MASS = "79.9663";
        private const string PHOS_MOD_RESIDUES = "STY";

        /// <summary>
        /// These columns correspond to the tab-delimited file created directly by InSpecT
        /// </summary>
        private enum InspectResultsFileColumns
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

        private enum InspectModType
        {
            Unknown = 0,
            DynamicMod = 1,
            StaticMod = 2,
            DynNTermPeptide = 3,
            DynCTermPeptide = 4
        }

        private enum FilteredOutputFileTypeConstants
        {
            SynFile = 0,
            FHTbyFScore = 1,
            FHTbyTotalPRM = 2
        }

        /// <summary>
        /// This data structure holds rows read from the tab-delimited file (Dataset.tsv) created InSpecT
        /// </summary>
        /// <remarks>
        /// These columns hold data that this class will use when creating the synopsis file
        /// </remarks>
        private struct InspectSearchResult
        {
            public string SpectrumFileName;
            public string Scan;
            public int ScanNum;
            public string PeptideAnnotation;
            public string Protein;
            public string Charge;
            public short ChargeNum;
            public string MQScore;                  // Larger values are better scores; note that MQScore can be negative
            public float MQScoreNum;                // Store the value of the string for quick reference when sorting
            public int Length;
            public string TotalPRMScore;            // Larger values are better scores
            public float TotalPRMScoreNum;          // We store the value of the string for quick reference when sorting
            public string MedianPRMScore;
            public string FractionY;
            public string FractionB;
            public string Intensity;
            public int NTT;
            public string PValue;                   // Smaller values are better scores
            public float PValueNum;                 // Store the value of the string for quick reference when sorting
            public string FScore;                   // Larger values are better scores
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

            /// <summary>
            /// Reset stored values to empty strings and zeros
            /// </summary>
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
                PValue = string.Empty;
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

            /// <summary>
            /// Show scan, peptide, and p-value
            /// </summary>
            public override string ToString()
            {
                return string.Format("Scan {0}: {1}, PValue {2}", ScanNum, PeptideAnnotation, PValue);
            }
        }

        private struct ModInfo
        {
            public string ModName;              // Mod names must be lower case, and 4 characters in length (or shorter)
            public string ModMass;              // Storing as a string since reading from a text file and writing to a text file
            public string Residues;
            public InspectModType ModType;
            public string ModSymbol;
        }

        /// <summary>
        /// When true, sort identifications by score when writing to first hits and synopsis files
        /// </summary>
        public bool SortFHTAndSynFiles { get; set; }

        private void AddCurrentRecordToSearchResults(ref int currentScanResultsCount,
            InspectSearchResult[] searchResultsCurrentScan,
            InspectSearchResult udtSearchResult)
        {
            if (currentScanResultsCount >= searchResultsCurrentScan.Length)
            {
                Array.Resize(ref searchResultsCurrentScan, searchResultsCurrentScan.Length * 2);
            }

            searchResultsCurrentScan[currentScanResultsCount] = udtSearchResult;
            currentScanResultsCount++;
        }

        private void AddDynamicAndStaticResidueMods(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
        {
            // Step through .PeptideSequenceWithMods
            // For each residue, check if a static mod is defined that affects that residue
            // For each mod symbol, determine the modification and add to searchResult

            var mostRecentLetter = '-';
            var residueLocInPeptide = 0;

            var sequence = searchResult.PeptideSequenceWithMods;
            for (var index = 0; index <= sequence.Length - 1; index++)
            {
                var character = sequence[index];

                if (StringUtilities.IsLetterAtoZ(character))
                {
                    mostRecentLetter = character;
                    residueLocInPeptide++;

                    for (var modIndex = 0; modIndex <= mPeptideMods.ModificationCount - 1; modIndex++)
                    {
                        if (mPeptideMods.GetModificationTypeByIndex(modIndex) != ModificationDefinition.ResidueModificationType.StaticMod)
                            continue;

                        var modificationDefinition = mPeptideMods.GetModificationByIndex(modIndex);

                        if (modificationDefinition.TargetResiduesContain(character))
                        {
                            // Match found; add this modification
                            var residueTerminusState = searchResult.DetermineResidueTerminusState(residueLocInPeptide);

                            searchResult.SearchResultAddModification(
                                modificationDefinition, character, residueLocInPeptide,
                                residueTerminusState, updateModOccurrenceCounts);
                        }
                    }
                }
                else if (StringUtilities.IsLetterAtoZ(mostRecentLetter))
                {
                    searchResult.SearchResultAddDynamicModification(character, mostRecentLetter, residueLocInPeptide, searchResult.DetermineResidueTerminusState(residueLocInPeptide), updateModOccurrenceCounts);
                }
                else
                {
                    // We found a modification symbol but mostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }
        }

        private readonly InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc sortScanChargeFScore = new();
        private readonly InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc sortScanChargeMQScore = new();
        private readonly InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc sortScanChargeTotalPRMDesc = new();

        /// <summary>
        /// Sorts the data by descending TotalPRMScore, than ranks each entry; in addition, computes normalized delta score (DeltaNorm) values
        /// </summary>
        /// <param name="searchResultsCurrentScan"></param>
        /// <param name="currentScanResultsCount"></param>
        private void AssignRankAndDeltaNormValues(ref InspectSearchResult[] searchResultsCurrentScan, int currentScanResultsCount)
        {
            const float DeltaNormMQScore_If_Undefined = 0;
            const float DeltaNormTotalPRMScore_If_Undefined = 0;

            var lastCharge = 0;
            double lastValue = 0;

            var currentRank = 0;

            // Sort searchResultsCurrentScan by ascending scan, ascending charge, and descending RankFScore
            // All of the data in searchResultsCurrentScan should have the same scan number
            Array.Sort(searchResultsCurrentScan, 0, currentScanResultsCount, sortScanChargeFScore);

            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index == 0 || searchResultsCurrentScan[index].ChargeNum != lastCharge)
                {
                    lastCharge = searchResultsCurrentScan[index].ChargeNum;
                    lastValue = searchResultsCurrentScan[index].FScoreNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(searchResultsCurrentScan[index].FScoreNum - lastValue) > float.Epsilon)
                    {
                        lastValue = searchResultsCurrentScan[index].FScoreNum;
                        currentRank++;
                    }
                }

                searchResultsCurrentScan[index].RankFScore = currentRank;
            }

            // Sort searchResultsCurrentScan by ascending scan, ascending charge, and descending MQScore (note that MQScore can be negative)
            // All of the data in searchResultsCurrentScan should have the same scan number
            Array.Sort(searchResultsCurrentScan, 0, currentScanResultsCount, sortScanChargeMQScore);

            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index < currentScanResultsCount - 1 && searchResultsCurrentScan[index].ChargeNum == searchResultsCurrentScan[index + 1].ChargeNum)
                {
                    searchResultsCurrentScan[index].DeltaNormMQScore = ComputeDeltaNormScore(searchResultsCurrentScan[index].MQScoreNum, searchResultsCurrentScan[index + 1].MQScoreNum, DeltaNormMQScore_If_Undefined);
                }
                else
                {
                    searchResultsCurrentScan[index].DeltaNormMQScore = 0;
                }
            }

            // Sort searchResultsCurrentScan by ascending scan, ascending charge, descending TotalPRMScore, and descending PValue
            // All of the data in searchResultsCurrentScan should have the same scan number
            Array.Sort(searchResultsCurrentScan, 0, currentScanResultsCount, sortScanChargeTotalPRMDesc);

            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index == 0 || searchResultsCurrentScan[index].ChargeNum != lastCharge)
                {
                    lastCharge = searchResultsCurrentScan[index].ChargeNum;
                    lastValue = searchResultsCurrentScan[index].TotalPRMScoreNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(searchResultsCurrentScan[index].TotalPRMScoreNum - lastValue) > float.Epsilon)
                    {
                        lastValue = searchResultsCurrentScan[index].TotalPRMScoreNum;
                        currentRank++;
                    }
                }

                searchResultsCurrentScan[index].RankTotalPRMScore = currentRank;

                if (index < currentScanResultsCount - 1 && searchResultsCurrentScan[index].ChargeNum == searchResultsCurrentScan[index + 1].ChargeNum)
                {
                    searchResultsCurrentScan[index].DeltaNormTotalPRMScore = ComputeDeltaNormScore(searchResultsCurrentScan[index].TotalPRMScoreNum, searchResultsCurrentScan[index + 1].TotalPRMScoreNum, DeltaNormTotalPRMScore_If_Undefined);
                }
                else
                {
                    searchResultsCurrentScan[index].DeltaNormTotalPRMScore = 0;
                }
            }
        }

        private bool AddModificationsAndComputeMass(SearchResultsBaseClass searchResult, bool updateModOccurrenceCounts)
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
                // Since InSpecT allows a terminal peptide residue to be modified twice, we'll allow that to happen,
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
        /// <param name="inspectModInfo">Used to replace Mod text entries in the peptides with Mod Symbols</param>
        /// <param name="filteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateFHTorSYNResultsFile(
            string inputFilePath,
            string outputFilePath,
            IReadOnlyList<ModInfo> inspectModInfo,
            FilteredOutputFileTypeConstants filteredOutputFileType)
        {
            var resultID = 0;

            bool success;

            var errorMessages = new List<string>();

            try
            {
                // Initialize variables
                var previousScan = int.MinValue;

                IComparer<InspectSearchResult> sortComparer = filteredOutputFileType switch
                {
                    // Writes the synopsis file, which writes every record with a p-value below a set threshold or a TotalPRMScore above a certain threshold
                    FilteredOutputFileTypeConstants.SynFile => new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc(),

                    // Write the PRM first-hits file, which writes the record with the highest TotalPRMScore
                    FilteredOutputFileTypeConstants.FHTbyTotalPRM => new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc(),

                    // Writes the synopsis file, which writes every record with a p-value below a set threshold or a TotalPRMScore above a certain threshold
                    FilteredOutputFileTypeConstants.FHTbyFScore => new InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc(),

                    _ => throw new ArgumentOutOfRangeException(nameof(filteredOutputFileType), filteredOutputFileType, null)
                };

                try
                {
                    // Open the input file and parse it
                    // Initialize the stream reader and the stream writer
                    using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                    using var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                    // Write the header line
                    WriteSynFHTFileHeader(writer, errorMessages);

                    errorMessages.Clear();
                    var resultsProcessed = 0;

                    // Initialize array that will hold all of the records for a given scan
                    var currentScanResultsCount = 0;
                    var searchResultsCurrentScan = new InspectSearchResult[10];

                    // Initialize the list that will hold all of the records that will ultimately be written out to disk
                    var filteredSearchResults = new List<InspectSearchResult>();

                    // Parse the input file
                    while (!reader.EndOfStream && !AbortProcessing)
                    {
                        var lineIn = reader.ReadLine();

                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        var validSearchResult = ParseInspectResultsFileEntry(lineIn, inspectModInfo, out var udtSearchResult, errorMessages, resultsProcessed);

                        resultsProcessed++;
                        if (!validSearchResult)
                        {
                            continue;
                        }

                        if (previousScan != int.MinValue && previousScan != udtSearchResult.ScanNum)
                        {
                            // New scan encountered; sort and filter the data in searchResultsCurrentScan, then call StoreTopFHTMatch or StoreSynMatches
                            if (filteredOutputFileType == FilteredOutputFileTypeConstants.SynFile)
                            {
                                StoreSynMatches(writer, ref resultID, currentScanResultsCount, searchResultsCurrentScan,
                                    filteredSearchResults, errorMessages, ref sortComparer);
                            }
                            else
                            {
                                StoreTopFHTMatch(writer, ref resultID, currentScanResultsCount, searchResultsCurrentScan,
                                    filteredSearchResults, errorMessages, ref sortComparer);
                            }

                            currentScanResultsCount = 0;
                        }

                        AddCurrentRecordToSearchResults(ref currentScanResultsCount, searchResultsCurrentScan, udtSearchResult);

                        previousScan = udtSearchResult.ScanNum;

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    // Store the last record
                    if (currentScanResultsCount > 0)
                    {
                        if (filteredOutputFileType == FilteredOutputFileTypeConstants.SynFile)
                        {
                            StoreSynMatches(writer, ref resultID, currentScanResultsCount, searchResultsCurrentScan, filteredSearchResults,
                                errorMessages, ref sortComparer);
                        }
                        else
                        {
                            StoreTopFHTMatch(writer, ref resultID, currentScanResultsCount, searchResultsCurrentScan, filteredSearchResults,
                                errorMessages, ref sortComparer);
                        }

                        currentScanResultsCount = 0;
                    }

                    if (SortFHTAndSynFiles)
                    {
                        // Sort the data in filteredSearchResults then write out to disk
                        SortAndWriteFilteredSearchResults(writer, filteredSearchResults, errorMessages);
                    }

                    // Inform the user if any errors occurred
                    if (errorMessages.Count > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                    }

                    success = true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in CreateFHTorSYNResultsFile", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                    success = false;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating the output file in CreateFHTorSYNResultsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        private double ComputeDelMCorrectedPPM(
            double precursorErrorDa,
            double precursorMonoMass,
            double peptideMonoisotopicMass,
            bool adjustPrecursorMassForC13)
        {
            return SearchResultsBaseClass.ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, peptideMonoisotopicMass, adjustPrecursorMassForC13);
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

        /// <summary>
        /// Compute the theoretical peptide MH using the precursor m/z value and the precursor error values
        /// </summary>
        /// <param name="precursorMZText"></param>
        /// <param name="precursorErrorText"></param>
        /// <param name="chargeText"></param>
        private double ComputePeptideMHFromPrecursorInfo(string precursorMZText, string precursorErrorText, string chargeText)
        {
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
                var charge = StringUtilities.CIntSafe(chargeText, -1);

                if (charge >= 1)
                {
                    if (double.TryParse(precursorMZText, out var precursorMZ))
                    {
                        if (double.TryParse(precursorErrorText, out var precursorError))
                        {
                            // Note: the October 2008 version of InSpecT uses an Absolute Value function when computing the PrecursorError; the version used by PNNL does not use Absolute Value
                            // Note: switched to compute (M+H)+ in August 2011; prior to this, we were computing uncharged monoisotopic mass
                            peptideMH = (precursorMZ - precursorError) * charge - (charge - 1) * PeptideMassCalculator.MASS_PROTON;
                        }
                    }
                }
            }

            return peptideMH;
        }

        /// <summary>
        /// Construct the peptide to protein map file path
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="outputDirectoryPath">Output directory path</param>
        /// <param name="mts">If true, the map file will end with MTS.txt; otherwise, just .txt</param>
        /// <returns>_PepToProtMap file that corresponds to the input file</returns>
        protected override string ConstructPepToProteinMapFilePath(string inputFilePath, string outputDirectoryPath, bool mts)
        {
            var suffixesToFind = new List<string> {
                "_inspect_syn",
                "_inspect_fht"
            };

            return ConstructPepToProteinMapFilePath(inputFilePath, outputDirectoryPath, mts, suffixesToFind, 4);
        }

        /// <summary>
        /// Read mod info from the InSpecT parameter file
        /// </summary>
        /// <param name="inspectParameterFilePath"></param>
        /// <param name="modList"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromInspectParamFile(string inspectParameterFilePath, out List<ModInfo> modList)
        {
            modList = new List<ModInfo>();

            try
            {
                var unnamedModID = 0;

                if (string.IsNullOrWhiteSpace(inspectParameterFilePath))
                {
                    SetErrorMessage("InSpecT Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                // Read the contents of the inspect parameter file
                using var reader = new StreamReader(new FileStream(inspectParameterFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(lineIn))
                        continue;

                    var dataLine = lineIn.Trim();
                    if (dataLine.Length == 0)
                        continue;

                    if (dataLine[0] == '#')
                    {
                        // Comment line; skip it
                        continue;
                    }

                    if (!dataLine.StartsWith("mod", StringComparison.OrdinalIgnoreCase))
                        continue;

                    // Modification definition line

                    // Split the line on commas
                    var splitLine = dataLine.Split(',');

                    if (splitLine.Length < 3 || splitLine[0].ToLower().Trim() != "mod")
                        continue;

                    var modDef = new ModInfo()
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
                                modDef.ModType = InspectModType.DynamicMod;
                                break;

                            case "fix":
                                modDef.ModType = InspectModType.StaticMod;
                                break;

                            case "nterminal":
                                modDef.ModType = InspectModType.DynNTermPeptide;
                                break;

                            case "cterminal":
                                modDef.ModType = InspectModType.DynCTermPeptide;
                                break;

                            default:
                                OnWarningEvent("Unrecognized Mod Type in the InSpecT parameter file");
                                modDef.ModType = InspectModType.DynamicMod;
                                break;
                        }
                    }
                    else
                    {
                        // Assume dynamic if not specified
                        modDef.ModType = InspectModType.DynamicMod;
                    }

                    if (splitLine.Length >= 5)
                    {
                        modDef.ModName = splitLine[4].ToLower();
                        if (modDef.ModName.Length > 4)
                        {
                            // Only keep the first 4 characters of the modification name
                            modDef.ModName = modDef.ModName.Substring(0, 4);
                        }
                    }
                    else
                    {
                        unnamedModID++;
                        modDef.ModName = "UnnamedMod" + unnamedModID.ToString();
                    }

                    // Check for phosphorylation
                    // InSpecT requires that it be defined in the parameter file as: mod,80,STY,opt,phosphorylation
                    //  However, we want to use the more precise mass of 79.9663
                    if (modDef.ModName == PHOS_MOD_NAME.ToLower() && modDef.ModMass == "80")
                    {
                        modDef.ModMass = PHOS_MOD_MASS;
                    }

                    modList.Add(modDef);
                }

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the InSpecT parameter file (" + Path.GetFileName(inspectParameterFilePath) + ")", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingModificationDefinitionsFile);
                return false;
            }
        }

        private static readonly Regex RegexScanNumberRegEx = new(DTA_FILENAME_SCAN_NUMBER_REGEX, REGEX_OPTIONS);

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
        /// Load the PeptideToProteinMap information
        /// In addition, creates the _inspect_PepToProtMapMTS.txt file with the new mod symbols and corrected termini symbols
        /// </summary>
        /// <param name="pepToProteinMapFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="inspectModInfo"></param>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="mtsPepToProteinMapFilePath"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool LoadPeptideToProteinMapInfoInspect(
            string pepToProteinMapFilePath,
            string outputDirectoryPath,
            IReadOnlyList<ModInfo> inspectModInfo,
            ref List<PepToProteinMapping> pepToProteinMapping,
            ref string mtsPepToProteinMapFilePath)
        {
            bool success;

            try
            {
                mtsPepToProteinMapFilePath = string.Empty;

                if (string.IsNullOrWhiteSpace(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    OnWarningEvent("PepToProteinMap file is not defined");
                    return false;
                }

                if (!File.Exists(pepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    OnWarningEvent("PepToProteinMap file does not exist: " + pepToProteinMapFilePath);
                    return false;
                }

                // Initialize pepToProteinMapping
                pepToProteinMapping = new List<PepToProteinMapping>();

                // Read the data in the peptide to protein map file
                success = LoadPeptideToProteinMapInfo(pepToProteinMapFilePath, pepToProteinMapping, out var headerLine);

                if (success)
                {
                    mtsPepToProteinMapFilePath = Path.Combine(outputDirectoryPath, Path.GetFileNameWithoutExtension(pepToProteinMapFilePath) + "MTS.txt");

                    using var writer = new StreamWriter(new FileStream(mtsPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                    if (!string.IsNullOrEmpty(headerLine))
                    {
                        // Header line
                        writer.WriteLine(headerLine);
                    }

                    for (var index = 0; index <= pepToProteinMapping.Count - 1; index++)
                    {
                        // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                        // In addition, replace the * terminus symbols with dashes
                        var mtsCompatiblePeptide = ReplaceInspectModTextWithSymbol(ReplaceTerminus(pepToProteinMapping[index].Peptide), inspectModInfo);

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
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" + Path.GetFileName(mtsPepToProteinMapFilePath) + ")", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        /// <summary>
        /// Parse the header line of an Inspect _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseInspectSynFileHeaderLine(string lineIn, IDictionary<InspectSynFileColumns, int> columnMapping)
        {
            var columnNames = InspectSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (InspectSynFileColumns resultColumn in Enum.GetValues(typeof(InspectSynFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');
                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the InSpecT synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse and InSpecT synopsis file
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="pepToProteinMapping"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseInspectSynopsisFile(string inputFilePath, string outputDirectoryPath, ref List<PepToProteinMapping> pepToProteinMapping, bool resetMassCorrectionTagsAndModificationDefinitions)
        {
            // Note that InSpecT synopsis files are normally sorted on TotalPRMScore descending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            // we will keep track of the scan, charge, and peptide information parsed for each unique TotalPRMScore encountered
            // (see peptidesFoundForTotalPRMScoreLevel below)

            var currentPeptideWithMods = string.Empty;

            var columnMapping = new Dictionary<InspectSynFileColumns, int>();

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
                var searchResult = new InSpecTResults(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize a SortedSet that will be used to avoid double-counting the same PSM in the same scan
                // This is required since a PSM with multiple proteins will be listed on multiple lines in the synopsis file
                // Values are PeptideSequenceWithMods_Scan_Charge

                var peptidesFoundForTotalPRMScoreLevel = new SortedSet<string>();

                var previousTotalPRMScore = string.Empty;

                // Assure that pepToProteinMapping is sorted on peptide
                if (pepToProteinMapping.Count > 1)
                {
                    pepToProteinMapping.Sort(new PepToProteinMappingComparer());
                }

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(Options);

                    var errorMessages = new List<string>();

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                    var headerParsed = false;

                    // Create the output files
                    var baseOutputFilePath = Path.Combine(outputDirectoryPath, Path.GetFileName(inputFilePath));
                    var filesInitialized = InitializeSequenceOutputFiles(baseOutputFilePath);
                    if (!filesInitialized)
                        return false;

                    // Parse the input file
                    while (!reader.EndOfStream && !AbortProcessing)
                    {
                        var lineIn = reader.ReadLine();
                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        var dataLine = true;

                        if (!headerParsed)
                        {
                            var validHeader = ParseInspectSynFileHeaderLine(lineIn, columnMapping);
                            if (validHeader)
                            {
                                dataLine = false;
                            }
                            else
                            {
                                // Error parsing header; assume this is a data line
                            }
                            headerParsed = true;
                        }

                        bool validSearchResult;
                        if (dataLine)
                        {
                            validSearchResult = ParseInspectSynFileEntry(lineIn, columnMapping, searchResult, errorMessages, out currentPeptideWithMods);
                        }
                        else
                        {
                            validSearchResult = false;
                        }

                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;
                        bool firstMatchForGroup;

                        if (searchResult.TotalPRMScore == previousTotalPRMScore)
                        {
                            // New result has the same TotalPRMScore as the previous result
                            // See if peptidesFoundForTotalPRMScoreLevel contains the peptide, scan and charge

                            if (peptidesFoundForTotalPRMScoreLevel.Contains(key))
                            {
                                firstMatchForGroup = false;
                            }
                            else
                            {
                                peptidesFoundForTotalPRMScoreLevel.Add(key);
                                firstMatchForGroup = true;
                            }
                        }
                        else
                        {
                            // New TotalPRMScore
                            // Reset peptidesFoundForTotalPRMScoreLevel
                            peptidesFoundForTotalPRMScoreLevel.Clear();

                            // Update previousTotalPRMScore
                            previousTotalPRMScore = searchResult.TotalPRMScore;

                            // Append a new entry to peptidesFoundForTotalPRMScoreLevel
                            peptidesFoundForTotalPRMScoreLevel.Add(key);
                            firstMatchForGroup = true;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup);
                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
                            errorMessages.Add(string.Format("Error adding modifications to sequence for ResultID '{0}'", searchResult.ResultID));
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
                                var currentProtein = searchResult.ProteinName;
                                do
                                {
                                    if (pepToProteinMapping[pepToProteinMapIndex].Protein != currentProtein)
                                    {
                                        searchResult.ProteinName = pepToProteinMapping[pepToProteinMapIndex].Protein;
                                        SaveResultsFileEntrySeqInfo(searchResult, false);
                                    }

                                    pepToProteinMapIndex++;
                                } while (pepToProteinMapIndex < pepToProteinMapping.Count &&
                                         currentPeptideWithMods == pepToProteinMapping[pepToProteinMapIndex].Peptide);
                            }
                            else
                            {
                                // Match not found; this is unexpected
                                OnWarningEvent("no match for '" + currentPeptideWithMods + "' in pepToProteinMapping");
                            }
                        }

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    if (Options.CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(inputFilePath);
                        var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                        modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                        SaveModificationSummaryFile(modificationSummaryFilePath);
                    }

                    // Inform the user if any errors occurred
                    if (errorMessages.Count > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                    }

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in ParseInspectSynopsisFile", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                    return false;
                }
                finally
                {
                    CloseSequenceOutputFiles();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating the output file in ParseInspectSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of an InSpecT results file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="inspectModInfo"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="resultsProcessed"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseInspectResultsFileEntry(
            string lineIn,
            IReadOnlyList<ModInfo> inspectModInfo,
            out InspectSearchResult udtSearchResult,
            ICollection<string> errorMessages,
            int resultsProcessed)
        {
            // Parses an entry from the InSpecT results file
            // The expected header line is:
            // #SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos PrecursorMZ	PrecursorError DelM_PPM

            udtSearchResult = new InspectSearchResult();

            string[] splitLine = null;

            try
            {
                udtSearchResult.Clear();

                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length >= 15)
                {
                    if (resultsProcessed == 0)
                    {
                        // This is the first line of the file; it may be a header row
                        // Determine this by seeing if any of the first three columns contains a number
                        if (!(SynFileReaderBaseClass.IsNumber(splitLine[0]) ||
                              SynFileReaderBaseClass.IsNumber(splitLine[1]) ||
                              SynFileReaderBaseClass.IsNumber(splitLine[2])))
                        {
                            // This is a header line; ignore it
                            return false;
                        }
                    }

                    udtSearchResult.SpectrumFileName = splitLine[(int)InspectResultsFileColumns.SpectrumFile];
                    if (splitLine[(int)InspectResultsFileColumns.Scan] == "0")
                    {
                        udtSearchResult.Scan = ExtractScanNumFromDTAName(udtSearchResult.SpectrumFileName);
                    }
                    else
                    {
                        udtSearchResult.Scan = splitLine[(int)InspectResultsFileColumns.Scan];
                    }
                    udtSearchResult.ScanNum = StringUtilities.CIntSafe(udtSearchResult.Scan, 0);

                    // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                    // In addition, replace the * terminus symbols with dashes
                    udtSearchResult.PeptideAnnotation = ReplaceInspectModTextWithSymbol(ReplaceTerminus(splitLine[(int)InspectResultsFileColumns.Annotation]), inspectModInfo);
                    udtSearchResult.Protein = TruncateProteinName(splitLine[(int)InspectResultsFileColumns.Protein]);

                    udtSearchResult.Charge = splitLine[(int)InspectResultsFileColumns.Charge];
                    udtSearchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(udtSearchResult.Charge, 0));

                    udtSearchResult.MQScore = splitLine[(int)InspectResultsFileColumns.MQScore];
                    udtSearchResult.MQScoreNum = StringUtilities.CSngSafe(udtSearchResult.MQScore, 0);

                    udtSearchResult.Length = StringUtilities.CIntSafe(splitLine[(int)InspectResultsFileColumns.Length], 0);

                    udtSearchResult.TotalPRMScore = splitLine[(int)InspectResultsFileColumns.TotalPRMScore];
                    udtSearchResult.TotalPRMScoreNum = StringUtilities.CSngSafe(udtSearchResult.TotalPRMScore, 0);

                    udtSearchResult.MedianPRMScore = splitLine[(int)InspectResultsFileColumns.MedianPRMScore];
                    udtSearchResult.FractionY = RemoveExtraneousDigits(splitLine[(int)InspectResultsFileColumns.FractionY]);
                    udtSearchResult.FractionB = RemoveExtraneousDigits(splitLine[(int)InspectResultsFileColumns.FractionB]);
                    udtSearchResult.Intensity = splitLine[(int)InspectResultsFileColumns.Intensity];
                    udtSearchResult.NTT = StringUtilities.CIntSafe(splitLine[(int)InspectResultsFileColumns.NTT], 0);

                    udtSearchResult.PValue = RemoveExtraneousDigits(splitLine[(int)InspectResultsFileColumns.PValue]);
                    udtSearchResult.PValueNum = StringUtilities.CSngSafe(udtSearchResult.PValue, 0);

                    udtSearchResult.FScore = splitLine[(int)InspectResultsFileColumns.FScore];
                    udtSearchResult.FScoreNum = StringUtilities.CSngSafe(udtSearchResult.FScore, 0);

                    udtSearchResult.DeltaScore = splitLine[(int)InspectResultsFileColumns.DeltaScore];
                    udtSearchResult.DeltaScoreOther = splitLine[(int)InspectResultsFileColumns.DeltaScoreOther];

                    udtSearchResult.RecordNumber = splitLine[(int)InspectResultsFileColumns.RecordNumber];
                    udtSearchResult.DBFilePos = splitLine[(int)InspectResultsFileColumns.DBFilePos];
                    udtSearchResult.SpecFilePos = splitLine[(int)InspectResultsFileColumns.SpecFilePos];

                    if (splitLine.Length >= (int)InspectResultsFileColumns.PrecursorError + 1)
                    {
                        // InSpecT version 2008-10-14 added these two Precursor mass columns
                        udtSearchResult.PrecursorMZ = splitLine[(int)InspectResultsFileColumns.PrecursorMZ];
                        udtSearchResult.PrecursorError = splitLine[(int)InspectResultsFileColumns.PrecursorError];

                        udtSearchResult.MH = ComputePeptideMHFromPrecursorInfo(udtSearchResult.PrecursorMZ, udtSearchResult.PrecursorError, udtSearchResult.Charge);

                        if (double.TryParse(udtSearchResult.PrecursorMZ, out var precursorMZ))
                        {
                            var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, udtSearchResult.ChargeNum, 0);
                            var peptideMonoisotopicMass = udtSearchResult.MH - PeptideMassCalculator.MASS_PROTON;

                            var precursorErrorDa = precursorMonoMass - peptideMonoisotopicMass;

                            var peptideDeltaMassCorrectedPpm =
                                ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, peptideMonoisotopicMass, true);

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

                    return true;
                }

                return false;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing InSpecT results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing InSpecT Results in ParseInspectResultsFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of an InSpecT _syn.txt file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="peptideSequenceWithMods"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseInspectSynFileEntry(
            string lineIn,
            IDictionary<InspectSynFileColumns, int> columnMapping,
            InSpecTResults searchResult,
            ICollection<string> errorMessages,
            out string peptideSequenceWithMods)
        {
            string[] splitLine = null;

            peptideSequenceWithMods = string.Empty;

            try
            {
                // Reset searchResult
                searchResult.Clear();

                splitLine = lineIn.TrimEnd().Split('\t');
                if (splitLine.Length < INSPECT_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.ResultID], out int resultId))
                {
                    ReportError("ResultID column is missing or invalid", true);
                }

                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.Charge], out string charge);

                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.Protein], out string proteinName);

                searchResult.ResultID = resultId;
                searchResult.Scan = scan;
                searchResult.Charge = charge;
                searchResult.ProteinName = proteinName;

                searchResult.MultipleProteinCount = "0";

                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.Peptide], out peptideSequenceWithMods);

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequenceWithMods, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since InSpecT only outputs the prefix and suffix letters for the first protein
                searchResult.ComputePeptideCleavageStateInProtein();

                if (GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.PrecursorError], out string peptideDeltaMass))
                {
                    searchResult.PeptideDeltaMass = peptideDeltaMass;
                    // Note: .peptideDeltaMass is stored in the InSpecT results file as "Observed_Mass - Theoretical_Mass"
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

                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.MQScore], out string mqScore);

                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.Length], out string length);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.TotalPRMScore], out string totalPrmScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.MedianPRMScore], out string medianPrmScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.FractionY], out string fractionY);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.FractionB], out string fractionB);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.Intensity], out string intensity);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.NTT], out string ntt);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.PValue], out string pValue);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.FScore], out string fScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.DeltaScore], out string deltaScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.DeltaScoreOther], out string deltaScoreOther);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.DeltaNormMQScore], out string deltaNormMqScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.DeltaNormTotalPRMScore], out string deltaNormTotalPrmScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.RankTotalPRMScore], out string rankTotalPrmScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.RankFScore], out string rankFScore);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.MH], out string peptideMh);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.RecordNumber], out string recordNumber);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.DBFilePos], out string dbFilePos);
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.SpecFilePos], out string specFilePos);

                // Note: .PrecursorError was processed earlier in this function
                GetColumnValue(splitLine, columnMapping[InspectSynFileColumns.PrecursorMZ], out string precursorMz);

                searchResult.MQScore = mqScore;
                searchResult.Length = length;
                searchResult.TotalPRMScore = totalPrmScore;
                searchResult.MedianPRMScore = medianPrmScore;
                searchResult.FractionY = fractionY;
                searchResult.FractionB = fractionB;
                searchResult.Intensity = intensity;
                searchResult.NTT = ntt;
                searchResult.PValue = pValue;
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

                return true;
            }
            catch (Exception ex)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing InSpecT results for RowIndex '{0}': {1}", splitLine[0], ex.Message));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing InSpecT Results in ParseInspectSynFileEntry: " + ex.Message);
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">InSpecT results file (Dataset_inspect.txt)</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if successful, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            var mtsPepToProteinMapFilePath = string.Empty;

            var success = false;

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                    return false;
                }

                success = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!success)
                {
                    return false;
                }

                ResetProgress("Parsing " + Path.GetFileName(inputFilePath));

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);

                    if (inputFile.Directory == null)
                    {
                        OnWarningEvent("InSpecTResultsProcessor.ProcessFile: Could not determine the parent directory of " + inputFile.FullName);
                        return false;
                    }

                    var pepToProteinMapping = new List<PepToProteinMapping>();

                    // Load the InSpecT Parameter File so that we can determine the modification names and masses
                    if (!ExtractModInfoFromInspectParamFile(Options.SearchToolParameterFilePath, out var inspectModInfo))
                    {
                        if (inspectModInfo == null || inspectModInfo.Count == 0)
                        {
                            inspectModInfo = new List<ModInfo>();
                            var modDef = new ModInfo()
                            {
                                ModName = PHOS_MOD_NAME.ToLower(),
                                ModMass = PHOS_MOD_MASS,
                                Residues = PHOS_MOD_RESIDUES,
                                ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL.ToString()
                            };
                            inspectModInfo.Add(modDef);
                        }
                    }

                    // Resolve the mods in mInspectModInfo with the ModDefs mods
                    ResolveInspectModsWithModDefinitions(inspectModInfo);

                    if (Options.CreateFirstHitsFile)
                    {
                        // Create the first hits output file
                        ResetProgress("Creating the FHT file (top TotalPRMScore)", true);

                        var outputFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
                        outputFilePath = Path.Combine(outputDirectoryPath, outputFilePath + INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX);

                        success = CreateFHTorSYNResultsFile(inputFilePath, outputFilePath, inspectModInfo, FilteredOutputFileTypeConstants.FHTbyTotalPRM);

                        // Create the first hits output file
                        ResetProgress("Creating the FHT file (top FScore)", true);

                        outputFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
                        outputFilePath = Path.Combine(outputDirectoryPath, outputFilePath + INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX);

                        success = CreateFHTorSYNResultsFile(inputFilePath, outputFilePath, inspectModInfo, FilteredOutputFileTypeConstants.FHTbyFScore);
                    }

                    if (Options.CreateSynopsisFile)
                    {
                        // Create the synopsis output file
                        ResetProgress("Creating the SYN file", true);

                        // Define the synopsis output file name based on inputFilePath
                        var synOutputFilePath = Path.GetFileNameWithoutExtension(inputFilePath);
                        synOutputFilePath = Path.Combine(outputDirectoryPath, synOutputFilePath + SYNOPSIS_FILE_SUFFIX);

                        success = CreateFHTorSYNResultsFile(inputFilePath, synOutputFilePath, inspectModInfo, FilteredOutputFileTypeConstants.SynFile);
                        if (!success)
                        {
                            return false;
                        }

                        // Load the PeptideToProteinMap information; if the file doesn't exist, a warning will be displayed, but processing will continue
                        // LoadPeptideToProteinMapInfoInspect also creates _inspect_PepToProtMapMTS.txt file with the new mod symbols and corrected termini symbols
                        var pepToProteinMapFilePath = Path.Combine(inputFile.Directory.FullName, Path.GetFileNameWithoutExtension(inputFile.Name) + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + ".txt");

                        ResetProgress("Loading the PepToProtein map file: " + Path.GetFileName(pepToProteinMapFilePath), true);

                        LoadPeptideToProteinMapInfoInspect(pepToProteinMapFilePath, outputDirectoryPath, inspectModInfo, ref pepToProteinMapping, ref mtsPepToProteinMapFilePath);

                        // Create the other PHRP-specific files
                        ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                        success = ParseInspectSynopsisFile(synOutputFilePath, outputDirectoryPath, ref pepToProteinMapping, false);
                        if (!success)
                        {
                            return false;
                        }

                        // Remove all items from pepToProteinMapping to reduce memory overhead
                        pepToProteinMapping.Clear();
                        pepToProteinMapping.TrimExcess();

                        if (Options.CreateProteinModsFile)
                        {
                            if (string.IsNullOrWhiteSpace(synOutputFilePath))
                            {
                                OnWarningEvent("InSpecTResultsProcessor.ProcessFile: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                            }
                            else
                            {
                                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)),
                                                               outputDirectoryPath);

                                // Create the Protein Mods file
                                success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                                      PeptideHitResultTypes.Inspect);
                            }
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
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile:" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        private static readonly Regex RegexNumPlusZeroes = new(@"(\.\d*[1-9])0+$", RegexOptions.Compiled);
        private static readonly Regex RegexAllZeroes = new(@"\.0+$", RegexOptions.Compiled);

        /// <summary>
        /// If value ends in .0000, remove the .0000 portion
        /// </summary>
        /// <param name="value"></param>
        private string RemoveExtraneousDigits(string value)
        {
            if (string.IsNullOrWhiteSpace(value))
            {
                return string.Empty;
            }

            var match = RegexAllZeroes.Match(value);
            if (match.Success && match.Index > 0)
            {
                value = value.Substring(0, match.Index);
            }
            else
            {
                match = RegexNumPlusZeroes.Match(value);
                if (match.Success && match.Index > 0)
                {
                    if (match.Groups.Count > 1)
                    {
                        // Number is of the form 1.0030 or 1.300 or 1.030
                        value = value.Substring(0, match.Index) + match.Groups[1].Value;
                    }
                }
            }

            return value;
        }

        private static readonly Regex NTerminalModMassMatcher = new(INSPECT_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);
        private static readonly Regex CTerminalModMassMatcher = new(INSPECT_CTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Replaces modification name text in peptide sequences with modification symbols (uses case-sensitive comparisons)
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="inspectModInfo">This function assumes that each entry in inspectModInfo has both .ModName and .ModSymbol defined</param>
        private string ReplaceInspectModTextWithSymbol(string peptide, IReadOnlyList<ModInfo> inspectModInfo)
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
            for (var index = 0; index <= inspectModInfo.Count - 1; index++)
            {
                if (inspectModInfo[index].ModType == InspectModType.StaticMod)
                {
                    continue;
                }

                peptide = peptide.Replace(inspectModInfo[index].ModName, inspectModInfo[index].ModSymbol);

                if (inspectModInfo[index].ModType != InspectModType.DynNTermPeptide &&
                    inspectModInfo[index].ModType != InspectModType.DynCTermPeptide)
                {
                    continue;
                }

                Match match;
                if (inspectModInfo[index].ModType == InspectModType.DynNTermPeptide)
                {
                    // InSpecT notates N-terminal mods like this: R.+14HVIFLAER.R   (Note: This behavior is not yet confirmed)
                    // Look for this using a RegEx
                    match = NTerminalModMassMatcher.Match(peptide);
                }
                else if (inspectModInfo[index].ModType == InspectModType.DynCTermPeptide)
                {
                    // InSpecT notates C-terminal mods like this: R.HVIFLAER+14.R
                    // Look for this using a RegEx
                    match = CTerminalModMassMatcher.Match(peptide);
                }
                else
                {
                    // This code should never be reached
                    match = null;
                }

                if (match == null)
                {
                    continue;
                }

                if (match.Success && match.Groups.Count > 1)
                {
                    // Match found
                    try
                    {
                        var modMass = Convert.ToInt32(match.Groups[1].Value);

                        // Compare the mod mass in the specification to this Mod's mod mass
                        // If they are less than 0.5 Da apart, assume we have a match; yes, this assumption is a bit flaky
                        if (Math.Abs(modMass - Convert.ToDouble(inspectModInfo[index].ModMass)) <= 0.5)
                        {
                            // Match found
                            // Replace the matched region with .ModSymbol

                            string peptideNew;

                            if (match.Groups[0].Index > 0)
                            {
                                peptideNew = peptide.Substring(0, match.Groups[0].Index);
                            }
                            else
                            {
                                peptideNew = string.Empty;
                            }

                            peptideNew += inspectModInfo[index].ModSymbol;

                            if (match.Groups[0].Index + match.Groups[0].Length < peptide.Length)
                            {
                                peptideNew += peptide.Substring(match.Groups[0].Index + match.Groups[0].Length);
                            }

                            peptide = peptideNew;
                        }
                    }
                    catch (Exception ex)
                    {
                        ReportError("Error comparing mod mass in peptide to mod mass in inspectModInfo", true, ex);
                    }
                }
            }

            return prefix + peptide + suffix;
        }

        private string ReplaceTerminus(string peptide)
        {
            if (peptide.StartsWith(N_TERMINUS_SYMBOL_INSPECT))
            {
                peptide = PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + peptide.Substring(N_TERMINUS_SYMBOL_INSPECT.Length);
            }

            if (peptide.EndsWith(C_TERMINUS_SYMBOL_INSPECT))
            {
                peptide = peptide.Substring(0, peptide.Length - C_TERMINUS_SYMBOL_INSPECT.Length) + "." + PeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return peptide;
        }

        private void ResolveInspectModsWithModDefinitions(IList<ModInfo> inspectModInfo)
        {
            if (inspectModInfo == null)
            {
                return;
            }

            // Call .LookupModificationDefinitionByMass for each entry in inspectModInfo
            for (var index = 0; index <= inspectModInfo.Count - 1; index++)
            {
                var modDef = inspectModInfo[index];

                if (!double.TryParse(modDef.ModMass, out var modMass))
                {
                    continue;
                }

                int resIndexStart;
                int resIndexEnd;

                if (modDef.Residues.Length > 0)
                {
                    resIndexStart = 0;
                    resIndexEnd = modDef.Residues.Length - 1;
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
                        targetResidue = modDef.Residues[residueIndex];
                    }
                    else
                    {
                        targetResidue = default;
                    }

                    AminoAcidModInfo.ResidueTerminusState residueTerminusState;
                    if (modDef.ModType == InspectModType.DynNTermPeptide)
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                    }
                    else if (modDef.ModType == InspectModType.DynCTermPeptide)
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                    }
                    else
                    {
                        residueTerminusState = AminoAcidModInfo.ResidueTerminusState.None;
                    }

                    var modificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(modMass, targetResidue, residueTerminusState, out _, true);

                    if (residueIndex == resIndexStart)
                    {
                        modDef.ModSymbol = modificationDefinition.ModificationSymbol.ToString();
                    }
                }

                inspectModInfo[index] = modDef;
            }
        }

        private void StoreOrWriteSearchResult(
            TextWriter writer,
            ref int resultID,
            InspectSearchResult udtSearchResult,
            ICollection<InspectSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            if (SortFHTAndSynFiles)
            {
                filteredSearchResults.Add(udtSearchResult);
            }
            else
            {
                resultID++;
                WriteSearchResultToFile(resultID, writer, udtSearchResult, errorMessages);
            }
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<InspectSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            // Sort filteredSearchResults by descending TotalPRMScore, ascending scan, ascending charge, ascending peptide, and ascending protein
            var query = from item in filteredSearchResults orderby item.TotalPRMScoreNum descending, item.ScanNum, item.ChargeNum, item.PeptideAnnotation, item.Protein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, errorMessages);
                index++;
            }
        }

        private void StoreTopFHTMatch(
            TextWriter writer,
            ref int resultID,
            int currentScanResultsCount,
            InspectSearchResult[] searchResultsCurrentScan,
            ICollection<InspectSearchResult> filteredSearchResults,
            ICollection<string> errorMessages,
            ref IComparer<InspectSearchResult> sortComparer)
        {
            var currentCharge = short.MinValue;

            AssignRankAndDeltaNormValues(ref searchResultsCurrentScan, currentScanResultsCount);

            // Sort searchResultsCurrentScan by ascending scan, ascending charge, then descending TotalPRMScore or descending FScore (depending on sortComparer)
            // All of the data in searchResultsCurrentScan should have the same scan number
            Array.Sort(searchResultsCurrentScan, 0, currentScanResultsCount, sortComparer);

            // Now store or write out the first match for each charge for this scan
            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (index == 0 || currentCharge != searchResultsCurrentScan[index].ChargeNum)
                {
                    StoreOrWriteSearchResult(writer, ref resultID, searchResultsCurrentScan[index], filteredSearchResults, errorMessages);
                    currentCharge = searchResultsCurrentScan[index].ChargeNum;
                }
            }
        }

        private void StoreSynMatches(
            TextWriter writer,
            ref int resultID,
            int currentScanResultsCount,
            InspectSearchResult[] searchResultsCurrentScan,
            ICollection<InspectSearchResult> filteredSearchResults,
            ICollection<string> errorMessages,
            ref IComparer<InspectSearchResult> sortComparer)
        {
            AssignRankAndDeltaNormValues(ref searchResultsCurrentScan, currentScanResultsCount);

            // Sort searchResultsCurrentScan by ascending scan, ascending charge, descending TotalPRMScore, and descending FScore
            // All of the data in searchResultsCurrentScan should have the same scan number
            Array.Sort(searchResultsCurrentScan, 0, currentScanResultsCount, sortComparer);

            // Now store or write out the matches that pass the filters
            for (var index = 0; index <= currentScanResultsCount - 1; index++)
            {
                if (searchResultsCurrentScan[index].PValueNum <= Options.InspectSynopsisFilePValueThreshold ||
                    searchResultsCurrentScan[index].TotalPRMScoreNum >= TOTALPRMSCORE_THRESHOLD ||
                    searchResultsCurrentScan[index].FScoreNum >= FSCORE_THRESHOLD)
                {
                    StoreOrWriteSearchResult(writer, ref resultID, searchResultsCurrentScan[index], filteredSearchResults, errorMessages);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorMessages"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ICollection<string> errorMessages)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = InspectSynFileReader.GetColumnHeaderNamesAndIDs();

                var headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();

                writer.WriteLine(StringUtilities.CollapseList(headerNames));
            }
            catch (Exception)
            {
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add("Error writing synopsis / first hits header");
                }
            }
        }

        /// <summary>
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorMessages"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            InspectSearchResult udtSearchResult,
            ICollection<string> errorMessages)
        {
            try
            {
                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.Scan,
                    udtSearchResult.PeptideAnnotation,
                    udtSearchResult.Protein,
                    udtSearchResult.Charge,
                    udtSearchResult.MQScore,
                    udtSearchResult.Length.ToString(),
                    udtSearchResult.TotalPRMScore,
                    udtSearchResult.MedianPRMScore,
                    udtSearchResult.FractionY,
                    udtSearchResult.FractionB,
                    udtSearchResult.Intensity,
                    udtSearchResult.NTT.ToString(),
                    udtSearchResult.PValue,
                    udtSearchResult.FScore,
                    udtSearchResult.DeltaScore,
                    udtSearchResult.DeltaScoreOther,
                    PRISM.StringUtilities.DblToString(udtSearchResult.DeltaNormMQScore, 5),
                    PRISM.StringUtilities.DblToString(udtSearchResult.DeltaNormTotalPRMScore, 5),
                    udtSearchResult.RankTotalPRMScore.ToString(),
                    udtSearchResult.RankFScore.ToString(),
                    PRISM.StringUtilities.DblToString(udtSearchResult.MH, 6),
                    udtSearchResult.RecordNumber,
                    udtSearchResult.DBFilePos,
                    udtSearchResult.SpecFilePos,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.PrecursorError,
                    udtSearchResult.DelMPPM
                };

                writer.WriteLine(StringUtilities.CollapseList(data));
            }
            catch (Exception)
            {
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add("Error writing synopsis / first hits record");
                }
            }
        }

        /// <summary>
        /// Override this method to display the name of each class
        /// </summary>
        public override string ToString()
        {
            return TOOL_NAME + " results processor";
        }

        private class InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc : IComparer<InspectSearchResult>
        {
            public int Compare(InspectSearchResult x, InspectSearchResult y)
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

        private class InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc : IComparer<InspectSearchResult>
        {
            public int Compare(InspectSearchResult x, InspectSearchResult y)
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

        private class InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc : IComparer<InspectSearchResult>
        {
            public int Compare(InspectSearchResult x, InspectSearchResult y)
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
    }
}
