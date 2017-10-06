// This class reads in an InSpecT results file (txt format) and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.
//
// The modification definition information is determined from the InSpecT parameter file
//
// -------------------------------------------------------------------------------
// Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
// Program started August 12, 2008
//
// E-mail: john.sandoval@pnl.gov
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
using System.Text.RegularExpressions;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsInSpecTResultsProcessor : clsPHRPBaseClass
    {
        public clsInSpecTResultsProcessor() : base()
        {
            base.mFileDate = "June 28, 2013";
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
        public enum eInspectResultsFileColumns : int
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
            pvalue = 13,
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
        public enum eInspectSynFileColumns : int
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

        protected enum eInspectModType : int
        {
            Unknown = 0,
            DynamicMod = 1,
            StaticMod = 2,
            DynNTermPeptide = 3,
            DynCTermPeptide = 4
        }

        protected enum eFilteredOutputFileTypeConstants : int
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
        public bool SortFHTandSynFiles { get; set; }

        #endregion

        private void AddCurrentRecordToSearchResults(ref int intCurrentScanResultsCount,
            udtInspectSearchResultType[] udtSearchResultsCurrentScan,
            udtInspectSearchResultType udtSearchResult)
        {
            if (intCurrentScanResultsCount >= udtSearchResultsCurrentScan.Length)
            {
                Array.Resize(ref udtSearchResultsCurrentScan, udtSearchResultsCurrentScan.Length * 2);
            }

            udtSearchResultsCurrentScan[intCurrentScanResultsCount] = udtSearchResult;
            intCurrentScanResultsCount += 1;
        }

        private void AddDynamicAndStaticResidueMods(clsSearchResultsInSpecT objSearchResult, bool blnUpdateModOccurrenceCounts)
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
                                objSearchResult.SearchResultAddModification(ref objModificationDefinition, chChar, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);
                            }
                        }
                    }
                }
                else if (IsLetterAtoZ(chMostRecentLetter))
                {
                    objSearchResult.SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts);
                }
                else
                {
                    // We found a modification symbol but chMostRecentLetter is not a letter
                    // Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                }
            }
        }

        private readonly InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc objSortScanChargeFScore = new InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc();
        private readonly InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc objSortScanChargeMQScore = new InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc();
        private readonly InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc objSortScanChargeTotalPRMDesc = new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc();

        /// <summary>
        /// Sorts the data by descending TotalPRMScore, than ranks each entry; in addition, computes normalized delta score (DeltaNorm) values
        /// </summary>
        /// <param name="udtSearchResultsCurrentScan"></param>
        /// <param name="intCurrentScanResultsCount"></param>
        /// <remarks></remarks>
        private void AssignRankAndDeltaNormValues(ref udtInspectSearchResultType[] udtSearchResultsCurrentScan, int intCurrentScanResultsCount)
        {
            const float DeltaNormMQScore_If_Undefined = 0;
            const float DeltaNormTotalPRMScore_If_Undefined = 0;

            int intLastCharge = 0;
            double dblLastValue = 0;

            int intCurrentRank = 0;

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, and descending RankFScore
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortScanChargeFScore);

            for (var intIndex = 0; intIndex <= intCurrentScanResultsCount - 1; intIndex++)
            {
                if (intIndex == 0 || udtSearchResultsCurrentScan[intIndex].ChargeNum != intLastCharge)
                {
                    intLastCharge = udtSearchResultsCurrentScan[intIndex].ChargeNum;
                    dblLastValue = udtSearchResultsCurrentScan[intIndex].FScoreNum;
                    intCurrentRank = 1;
                }
                else
                {
                    if (Math.Abs(udtSearchResultsCurrentScan[intIndex].FScoreNum - dblLastValue) > float.Epsilon)
                    {
                        dblLastValue = udtSearchResultsCurrentScan[intIndex].FScoreNum;
                        intCurrentRank += 1;
                    }
                }

                udtSearchResultsCurrentScan[intIndex].RankFScore = intCurrentRank;
            }

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, and descending MQScore (note that MQScore can be negative)
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortScanChargeMQScore);

            for (var intIndex = 0; intIndex <= intCurrentScanResultsCount - 1; intIndex++)
            {
                if (intIndex < intCurrentScanResultsCount - 1 && udtSearchResultsCurrentScan[intIndex].ChargeNum == udtSearchResultsCurrentScan[intIndex + 1].ChargeNum)
                {
                    udtSearchResultsCurrentScan[intIndex].DeltaNormMQScore = ComputeDeltaNormScore(udtSearchResultsCurrentScan[intIndex].MQScoreNum, udtSearchResultsCurrentScan[intIndex + 1].MQScoreNum, DeltaNormMQScore_If_Undefined);
                }
                else
                {
                    udtSearchResultsCurrentScan[intIndex].DeltaNormMQScore = 0;
                }
            }

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending PValue
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortScanChargeTotalPRMDesc);

            for (var intIndex = 0; intIndex <= intCurrentScanResultsCount - 1; intIndex++)
            {
                if (intIndex == 0 || udtSearchResultsCurrentScan[intIndex].ChargeNum != intLastCharge)
                {
                    intLastCharge = udtSearchResultsCurrentScan[intIndex].ChargeNum;
                    dblLastValue = udtSearchResultsCurrentScan[intIndex].TotalPRMScoreNum;
                    intCurrentRank = 1;
                }
                else
                {
                    if (Math.Abs(udtSearchResultsCurrentScan[intIndex].TotalPRMScoreNum - dblLastValue) > float.Epsilon)
                    {
                        dblLastValue = udtSearchResultsCurrentScan[intIndex].TotalPRMScoreNum;
                        intCurrentRank += 1;
                    }
                }

                udtSearchResultsCurrentScan[intIndex].RankTotalPRMScore = intCurrentRank;

                if (intIndex < intCurrentScanResultsCount - 1 && udtSearchResultsCurrentScan[intIndex].ChargeNum == udtSearchResultsCurrentScan[intIndex + 1].ChargeNum)
                {
                    udtSearchResultsCurrentScan[intIndex].DeltaNormTotalPRMScore = ComputeDeltaNormScore(udtSearchResultsCurrentScan[intIndex].TotalPRMScoreNum, udtSearchResultsCurrentScan[intIndex + 1].TotalPRMScoreNum, DeltaNormTotalPRMScore_If_Undefined);
                }
                else
                {
                    udtSearchResultsCurrentScan[intIndex].DeltaNormTotalPRMScore = 0;
                }
            }
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsInSpecT objSearchResult, bool blnUpdateModOccurrenceCounts)
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
        /// This routine creates a first hits file or synopsis file from the output from InSpecT
        /// </summary>
        /// <param name="strInputFilePath"></param>
        /// <param name="strOutputFilePath"></param>
        /// <param name="udtInspectModInfo">Used to replace Mod text entries in the peptides with Mod Symbols; assumes each entryin </param>
        /// <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool CreateFHTorSYNResultsFile(
            string strInputFilePath,
            string strOutputFilePath,
            ref udtModInfoType[] udtInspectModInfo,
            eFilteredOutputFileTypeConstants eFilteredOutputFileType)
        {
            int intPreviousScan = 0;

            string strLineIn = null;

            var udtSearchResult = new udtInspectSearchResultType();

            int intCurrentScanResultsCount = 0;
            udtInspectSearchResultType[] udtSearchResultsCurrentScan = null;

            int intFilteredSearchResultCount = 0;
            udtInspectSearchResultType[] udtFilteredSearchResults = null;

            int intResultsProcessed = 0;
            var intResultID = 0;

            bool blnSuccess = false;
            bool blnValidSearchResult = false;

            string strErrorLog = string.Empty;
            IComparer<udtInspectSearchResultType> objSortComparer = default(IComparer<udtInspectSearchResultType>);

            try
            {
                // Initialize variables
                intPreviousScan = Int32.MinValue;

                if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                {
                    // Writes the synopsis file, which writes every record with a p-value below a set threshold or a TotalPRMScore above a certain threshold
                    objSortComparer = new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc();
                }
                else
                {
                    if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.FHTbyTotalPRM)
                    {
                        // Write the PRM first-hits file, which writes the record with the highest TotalPRMScore
                        objSortComparer = new InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc();
                    }
                    else
                    {
                        // eFilteredOutputFileTypeConstants.FHTbyFScore
                        // Write the FScore first-hits file, which writes the record with the highest FScore
                        objSortComparer = new InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc();
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
                            // Write the header line
                            WriteSynFHTFileHeader(swResultFile, ref strErrorLog);

                            strErrorLog = string.Empty;
                            intResultsProcessed = 0;

                            // Initialize array that will hold all of the records for a given scan
                            intCurrentScanResultsCount = 0;
                            udtSearchResultsCurrentScan = new udtInspectSearchResultType[10];

                            // Initialize the array that will hold all of the records that will ultimately be written out to disk
                            intFilteredSearchResultCount = 0;
                            udtFilteredSearchResults = new udtInspectSearchResultType[1000];

                            // Parse the input file
                            while (!srDataFile.EndOfStream & !base.AbortProcessing)
                            {
                                strLineIn = srDataFile.ReadLine();

                                if (string.IsNullOrWhiteSpace(strLineIn))
                                {
                                    continue;
                                }

                                // Initialize udtSearchResult
                                udtSearchResult.Clear();

                                blnValidSearchResult = ParseInspectResultsFileEntry(ref strLineIn, ref udtInspectModInfo, ref udtSearchResult, ref strErrorLog, intResultsProcessed);

                                if (blnValidSearchResult)
                                {
                                    if (intPreviousScan != Int32.MinValue && intPreviousScan != udtSearchResult.ScanNum)
                                    {
                                        // New scan encountered; sort and filter the data in udtSearchResultsCurrentScan, then call StoreTopFHTMatch or StoreSynMatches
                                        if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                                        {
                                            StoreSynMatches(swResultFile, ref intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, ref intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog, ref objSortComparer);
                                        }
                                        else
                                        {
                                            StoreTopFHTMatch(swResultFile, ref intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, ref intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog, ref objSortComparer);
                                        }
                                        intCurrentScanResultsCount = 0;
                                    }

                                    AddCurrentRecordToSearchResults(ref intCurrentScanResultsCount, udtSearchResultsCurrentScan, udtSearchResult);

                                    intPreviousScan = udtSearchResult.ScanNum;
                                }

                                // Update the progress
                                UpdateProgress(Convert.ToSingle(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100));

                                intResultsProcessed += 1;
                            }

                            // Store the last record
                            if (intCurrentScanResultsCount > 0)
                            {
                                if (eFilteredOutputFileType == eFilteredOutputFileTypeConstants.SynFile)
                                {
                                    StoreSynMatches(swResultFile, ref intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, ref intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog, ref objSortComparer);
                                }
                                else
                                {
                                    StoreTopFHTMatch(swResultFile, ref intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, ref intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog, ref objSortComparer);
                                }
                                intCurrentScanResultsCount = 0;
                            }

                            if (SortFHTandSynFiles)
                            {
                                // Sort the data in udtFilteredSearchResults then write out to disk
                                SortAndWriteFilteredSearchResults(swResultFile, intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog);
                            }
                        }
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
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        protected double ComputeDelMCorrectedPPM(
            double dblPrecursorErrorDa,
            double dblPrecursorMonoMass,
            double dblPeptideMonoisotopicMass,
            bool blnAdjustPrecursorMassForC13)
        {
            double dblPeptideDeltaMassCorrectedPpm = 0;

            dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, dblPeptideMonoisotopicMass);

            return dblPeptideDeltaMassCorrectedPpm;
        }

        private float ComputeDeltaNormScore(float sngCurrentScore, float sngNextScore, float sngValueIfCurrentScoreZero)
        {
            try
            {
                if (Math.Abs(sngCurrentScore) > float.Epsilon)
                {
                    return Math.Abs((sngCurrentScore - sngNextScore) / sngCurrentScore);
                }
                else
                {
                    return sngValueIfCurrentScoreZero;
                }
            }
            catch (Exception)
            {
                return sngValueIfCurrentScoreZero;
            }
        }

        private double ComputePeptideMHFromPrecursorInfo(string strPrecursorMZ, string strPrecursorError, string strCharge)
        {
            // Compute the theoretical peptide MH using the precursor m/z value and the precursor error values

            double dblPrecursorMZ = 0;
            double dblPrecursorError = 0;
            int intCharge = 0;
            double dblPeptideMH = 0;

            if (string.IsNullOrWhiteSpace(strPrecursorMZ) || strPrecursorMZ == "0")
            {
                // Precursor m/z is undefined; cannot continue
                dblPeptideMH = 0;
            }
            else if (string.IsNullOrWhiteSpace(strPrecursorError))
            {
                // Precursor error is undefined; cannot continue
                dblPeptideMH = 0;
            }
            else
            {
                intCharge = CIntSafe(strCharge, -1);

                if (intCharge >= 1)
                {
                    if (double.TryParse(strPrecursorMZ, out dblPrecursorMZ))
                    {
                        if (double.TryParse(strPrecursorError, out dblPrecursorError))
                        {
                            // Note: the October 2008 version of Inspect uses an Absolute Value function when computing the PrecursorError; the version used by PNNL does not use Absolute Value
                            // Note: switched to compute (M+H)+ in August 2011; prior to this, we were computing uncharged monoisotopic mass
                            dblPeptideMH = (dblPrecursorMZ - dblPrecursorError) * intCharge - (intCharge - 1) * clsPeptideMassCalculator.MASS_PROTON;
                        }
                    }
                }
            }

            return dblPeptideMH;
        }

        protected override string ConstructPepToProteinMapFilePath(string strInputFilePath, string strOutputFolderPath, bool MTS)
        {
            var strPepToProteinMapFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
            if (strPepToProteinMapFilePath.EndsWith("_inspect_syn", StringComparison.InvariantCultureIgnoreCase) ||
                strPepToProteinMapFilePath.EndsWith("_inspect_fht", StringComparison.InvariantCultureIgnoreCase))
            {
                // Remove _syn or _fht
                strPepToProteinMapFilePath = strPepToProteinMapFilePath.Substring(0, strPepToProteinMapFilePath.Length - 4);
            }

            return base.ConstructPepToProteinMapFilePath(strPepToProteinMapFilePath, strOutputFolderPath, MTS);
        }

        private bool ExtractModInfoFromInspectParamFile(string strInspectParameterFilePath, ref udtModInfoType[] udtModList)
        {
            string strLineIn = null;
            string[] strSplitLine = null;

            int intModCount = 0;
            int intUnnamedModID = 0;

            try
            {
                // Initialize udtModList and intUnnamedModID
                intModCount = 0;
                udtModList = new udtModInfoType[0];

                intUnnamedModID = 0;

                if (string.IsNullOrWhiteSpace(strInspectParameterFilePath))
                {
                    SetErrorMessage("Inspect Parameter File name not defined; unable to extract mod info");
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                    return false;
                }

                // Read the contents of the inspect parameter file
                using (var srInFile = new StreamReader(new FileStream(strInspectParameterFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        strLineIn = srInFile.ReadLine().Trim();

                        if (strLineIn.Length > 0)
                        {
                            if (strLineIn[0] == '#')
                            {
                                // Comment line; skip it
                            }
                            else if (strLineIn.StartsWith("mod", StringComparison.InvariantCultureIgnoreCase))
                            {
                                // Modification definition line

                                // Split the line on commas
                                strSplitLine = strLineIn.Split(',');

                                if (strSplitLine.Length >= 3 && strSplitLine[0].ToLower().Trim() == "mod")
                                {
                                    if (udtModList == null || udtModList.Length == 0)
                                    {
                                        udtModList = new udtModInfoType[1];
                                    }
                                    else if (intModCount >= udtModList.Length)
                                    {
                                        Array.Resize(ref udtModList, udtModList.Length * 2);
                                    }

                                    var udtMod = udtModList[intModCount];
                                    udtMod.ModMass = strSplitLine[1];
                                    udtMod.Residues = strSplitLine[2];
                                    udtMod.ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL.ToString();

                                    if (strSplitLine.Length >= 4)
                                    {
                                        switch (strSplitLine[3].ToLower())
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
                                        // Assume dynamic if not specifed
                                        udtMod.ModType = eInspectModType.DynamicMod;
                                    }

                                    if (strSplitLine.Length >= 5)
                                    {
                                        udtMod.ModName = strSplitLine[4].ToLower();
                                        if (udtMod.ModName.Length > 4)
                                        {
                                            // Only keep the first 4 characters of the modification name
                                            udtMod.ModName = udtMod.ModName.Substring(0, 4);
                                        }
                                    }
                                    else
                                    {
                                        intUnnamedModID += 1;
                                        udtMod.ModName = "UnnamedMod" + intUnnamedModID.ToString();
                                    }

                                    // Check for phosphorylation
                                    // Inspect requires that it be defined in the parameter file as: mod,80,STY,opt,phosphorylation
                                    //  However, we want to use the more precise mass of 79.9663
                                    if (udtMod.ModName == PHOS_MOD_NAME.ToLower() & udtMod.ModMass == "80")
                                    {
                                        udtMod.ModMass = PHOS_MOD_MASS;
                                    }

                                    intModCount += 1;
                                }
                            }
                        }
                    }
                }

                // Shrink udtModList to the appropriate length
                Array.Resize(ref udtModList, intModCount);

                Console.WriteLine();

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error reading the Inspect parameter file (" + Path.GetFileName(strInspectParameterFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile);
                return false;
            }
        }

        private static readonly Regex RegexScanNumberRegEx = new Regex(DTA_FILENAME_SCAN_NUMBER_REGEX, REGEX_OPTIONS);

        private string ExtractScanNumFromDTAName(string spectrumFile)
        {
            string scanNum = string.Empty;

            // See if strValue resembles a .Dta file name
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
            SortFHTandSynFiles = true;
        }

        /// <summary>
        /// Load the PeptideToProteinMap information; in addition, creates the _inspect_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
        /// </summary>
        /// <param name="strPepToProteinMapFilePath"></param>
        /// <param name="strOutputFolderPath"></param>
        /// <param name="udtInspectModInfo"></param>
        /// <param name="lstPepToProteinMapping"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected bool LoadPeptideToProteinMapInfoInspect(
            string strPepToProteinMapFilePath,
            string strOutputFolderPath,
            ref udtModInfoType[] udtInspectModInfo,
            ref List<udtPepToProteinMappingType> lstPepToProteinMapping,
            ref string strMTSPepToProteinMapFilePath)
        {
            string strHeaderLine = string.Empty;
            string strMTSCompatiblePeptide = null;

            var blnSuccess = false;

            try
            {
                strMTSPepToProteinMapFilePath = string.Empty;

                if (string.IsNullOrWhiteSpace(strPepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file is not defined");
                    return false;
                }
                else if (!File.Exists(strPepToProteinMapFilePath))
                {
                    Console.WriteLine();
                    ReportWarning("PepToProteinMap file does not exist: " + strPepToProteinMapFilePath);
                    return false;
                }

                // Initialize lstPepToProteinMapping
                lstPepToProteinMapping = new List<udtPepToProteinMappingType>();

                // Read the data in strProteinToPeptideMappingFilePath
                blnSuccess = LoadPeptideToProteinMapInfo(strPepToProteinMapFilePath, lstPepToProteinMapping, out strHeaderLine);

                if (blnSuccess)
                {
                    strMTSPepToProteinMapFilePath = Path.Combine(strOutputFolderPath, Path.GetFileNameWithoutExtension(strPepToProteinMapFilePath) + "MTS.txt");

                    using (var swOutFile = new StreamWriter(new FileStream(strMTSPepToProteinMapFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
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
                            strMTSCompatiblePeptide = ReplaceInspectModTextWithSymbol(ReplaceTerminus(lstPepToProteinMapping[intIndex].Peptide), ref udtInspectModInfo);

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
                SetErrorMessage("Error writing MTS-compatible Peptide to Protein Map File (" + Path.GetFileName(strMTSPepToProteinMapFilePath) + "): " + ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private bool ParseInspectSynFileHeaderLine(string strLineIn, ref int[] intColumnMapping)
        {
            // Parse the header line

            string[] strSplitLine = null;
            eInspectSynFileColumns eResultFileColumn = default(eInspectSynFileColumns);
            SortedDictionary<string, eInspectSynFileColumns> lstColumnNames = default(SortedDictionary<string, eInspectSynFileColumns>);
            lstColumnNames = new SortedDictionary<string, eInspectSynFileColumns>(StringComparer.InvariantCultureIgnoreCase);

            intColumnMapping = new int[InspectSynopsisFileColCount];

            lstColumnNames.Add("ResultID", eInspectSynFileColumns.ResultID);
            lstColumnNames.Add("Scan", eInspectSynFileColumns.Scan);
            lstColumnNames.Add("Peptide", eInspectSynFileColumns.Peptide);
            lstColumnNames.Add("Protein", eInspectSynFileColumns.Protein);
            lstColumnNames.Add("Charge", eInspectSynFileColumns.Charge);
            lstColumnNames.Add("MQScore", eInspectSynFileColumns.MQScore);
            lstColumnNames.Add("Length", eInspectSynFileColumns.Length);
            lstColumnNames.Add("TotalPRMScore", eInspectSynFileColumns.TotalPRMScore);
            lstColumnNames.Add("MedianPRMScore", eInspectSynFileColumns.MedianPRMScore);
            lstColumnNames.Add("FractionY", eInspectSynFileColumns.FractionY);
            lstColumnNames.Add("FractionB", eInspectSynFileColumns.FractionB);
            lstColumnNames.Add("Intensity", eInspectSynFileColumns.Intensity);
            lstColumnNames.Add("NTT", eInspectSynFileColumns.NTT);
            lstColumnNames.Add("PValue", eInspectSynFileColumns.PValue);
            lstColumnNames.Add("FScore", eInspectSynFileColumns.FScore);
            lstColumnNames.Add("DeltaScore", eInspectSynFileColumns.DeltaScore);
            lstColumnNames.Add("DeltaScoreOther", eInspectSynFileColumns.DeltaScoreOther);
            lstColumnNames.Add("DeltaNormMQScore", eInspectSynFileColumns.DeltaNormMQScore);
            lstColumnNames.Add("DeltaNormTotalPRMScore", eInspectSynFileColumns.DeltaNormTotalPRMScore);
            lstColumnNames.Add("RankTotalPRMScore", eInspectSynFileColumns.RankTotalPRMScore);
            lstColumnNames.Add("RankFScore", eInspectSynFileColumns.RankFScore);
            lstColumnNames.Add("MH", eInspectSynFileColumns.MH);
            lstColumnNames.Add("RecordNumber", eInspectSynFileColumns.RecordNumber);
            lstColumnNames.Add("DBFilePos", eInspectSynFileColumns.DBFilePos);
            lstColumnNames.Add("SpecFilePos", eInspectSynFileColumns.SpecFilePos);
            lstColumnNames.Add("PrecursorMZ", eInspectSynFileColumns.PrecursorMZ);
            lstColumnNames.Add("PrecursorError", eInspectSynFileColumns.PrecursorError);

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
                SetErrorMessage("Error parsing header in Inspect synopsis file: " + ex.Message);
                return false;
            }

            return true;
        }

        protected bool ParseInSpectSynopsisFile(string strInputFilePath, string strOutputFolderPath, ref List<udtPepToProteinMappingType> lstPepToProteinMapping, bool blnResetMassCorrectionTagsAndModificationDefinitions)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            string strPreviousTotalPRMScore = null;

            // Note that Inspect synopsis files are normally sorted on TotalPRMScore descending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique TotalPRMScore encountered

            Hashtable htPeptidesFoundForTotalPRMScoreLevel = default(Hashtable);

            string strKey = null;

            string strLineIn = null;
            string strModificationSummaryFilePath = null;

            clsSearchResultsInSpecT objSearchResult = default(clsSearchResultsInSpecT);

            int intResultsProcessed = 0;
            int intPepToProteinMapIndex = 0;
            float sngPercentComplete = 0;

            string strCurrentPeptideWithMods = string.Empty;
            string strCurrentProtein = null;

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
                objSearchResult = new clsSearchResultsInSpecT(mPeptideMods, mPeptideSeqMassCalculator);

                // Initialize htPeptidesFoundForTotalPRMScoreLevel
                htPeptidesFoundForTotalPRMScoreLevel = new Hashtable();
                strPreviousTotalPRMScore = string.Empty;

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

                            blnDataLine = true;

                            if (!blnHeaderParsed)
                            {
                                blnSuccess = ParseInspectSynFileHeaderLine(strLineIn, ref intColumnMapping);
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
                                blnValidSearchResult = ParseInSpectSynFileEntry(ref strLineIn, ref intColumnMapping, objSearchResult, ref strErrorLog, ref strCurrentPeptideWithMods);
                            }
                            else
                            {
                                blnValidSearchResult = false;
                            }

                            if (blnValidSearchResult)
                            {
                                strKey = objSearchResult.PeptideSequenceWithMods + "_" + objSearchResult.Scan + "_" + objSearchResult.Charge;

                                if (objSearchResult.TotalPRMScore == strPreviousTotalPRMScore)
                                {
                                    // New result has the same TotalPRMScore as the previous result
                                    // See if htPeptidesFoundForTotalPRMScoreLevel contains the peptide, scan and charge

                                    if (htPeptidesFoundForTotalPRMScoreLevel.ContainsKey(strKey))
                                    {
                                        blnFirstMatchForGroup = false;
                                    }
                                    else
                                    {
                                        htPeptidesFoundForTotalPRMScoreLevel.Add(strKey, 1);
                                        blnFirstMatchForGroup = true;
                                    }
                                }
                                else
                                {
                                    // New TotalPRMScore
                                    // Reset htPeptidesFoundForTotalPRMScoreLevel
                                    htPeptidesFoundForTotalPRMScoreLevel.Clear();

                                    // Update strPreviousTotalPRMScore
                                    strPreviousTotalPRMScore = objSearchResult.TotalPRMScore;

                                    // Append a new entry to htPeptidesFoundForTotalPRMScoreLevel
                                    htPeptidesFoundForTotalPRMScoreLevel.Add(strKey, 1);
                                    blnFirstMatchForGroup = true;
                                }

                                blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup);
                                if (!blnSuccess)
                                {
                                    if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                    {
                                        strErrorLog += "Error adding modifications to sequence at RowIndex '" + objSearchResult.ResultID + "'" + "\n";
                                    }
                                }
                                base.SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass)objSearchResult, blnFirstMatchForGroup);

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
                                                base.SaveResultsFileEntrySeqInfo((clsSearchResultsBaseClass)objSearchResult, false);
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

        private bool ParseInspectResultsFileEntry(
            ref string strLineIn,
            ref udtModInfoType[] udtInspectModInfo,
            ref udtInspectSearchResultType udtSearchResult,
            ref string strErrorLog,
            int intResultsProcessed)
        {
            // Parses an entry from the Inspect results file
            // The expected header line is:
            // #SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos PrecursorMZ	PrecursorError DelM_PPM

            string[] strSplitLine = null;

            bool blnValidSearchResult = false;

            try
            {
                // Set this to False for now
                blnValidSearchResult = false;

                // Reset udtSearchResult
                udtSearchResult.Clear();

                strSplitLine = strLineIn.Trim().Split('\t');

                if (strSplitLine.Length >= 15)
                {
                    if (intResultsProcessed == 0)
                    {
                        // This is the first line of the file; it may be a header row
                        // Determine this by seeing if any of the first three columns contains a number
                        if (!(clsPHRPParser.IsNumber(strSplitLine[0]) ||
                              clsPHRPParser.IsNumber(strSplitLine[1]) ||
                              clsPHRPParser.IsNumber(strSplitLine[2])))
                        {
                            // This is a header line; ignore it
                            blnValidSearchResult = false;
                            return false;
                        }
                    }

                    udtSearchResult.SpectrumFileName = strSplitLine[(int)eInspectResultsFileColumns.SpectrumFile];
                    if (strSplitLine[(int)eInspectResultsFileColumns.Scan] == "0")
                    {
                        udtSearchResult.Scan = ExtractScanNumFromDTAName(udtSearchResult.SpectrumFileName);
                    }
                    else
                    {
                        udtSearchResult.Scan = strSplitLine[(int)eInspectResultsFileColumns.Scan];
                    }
                    udtSearchResult.ScanNum = CIntSafe(udtSearchResult.Scan, 0);

                    // Replace any mod text names in the peptide sequence with the appropriate mod symbols
                    // In addition, replace the * terminus symbols with dashes
                    udtSearchResult.PeptideAnnotation = ReplaceInspectModTextWithSymbol(ReplaceTerminus(strSplitLine[(int)eInspectResultsFileColumns.Annotation]), ref udtInspectModInfo);
                    udtSearchResult.Protein = TruncateProteinName(strSplitLine[(int)eInspectResultsFileColumns.Protein]);

                    udtSearchResult.Charge = strSplitLine[(int)eInspectResultsFileColumns.Charge];
                    udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                    udtSearchResult.MQScore = strSplitLine[(int)eInspectResultsFileColumns.MQScore];
                    udtSearchResult.MQScoreNum = CSngSafe(udtSearchResult.MQScore, 0);

                    udtSearchResult.Length = CIntSafe(strSplitLine[(int)eInspectResultsFileColumns.Length], 0);

                    udtSearchResult.TotalPRMScore = strSplitLine[(int)eInspectResultsFileColumns.TotalPRMScore];
                    udtSearchResult.TotalPRMScoreNum = CSngSafe(udtSearchResult.TotalPRMScore, 0);

                    udtSearchResult.MedianPRMScore = strSplitLine[(int)eInspectResultsFileColumns.MedianPRMScore];
                    udtSearchResult.FractionY = RemoveExtraneousDigits(strSplitLine[(int)eInspectResultsFileColumns.FractionY]);
                    udtSearchResult.FractionB = RemoveExtraneousDigits(strSplitLine[(int)eInspectResultsFileColumns.FractionB]);
                    udtSearchResult.Intensity = strSplitLine[(int)eInspectResultsFileColumns.Intensity];
                    udtSearchResult.NTT = CIntSafe(strSplitLine[(int)eInspectResultsFileColumns.NTT], 0);

                    udtSearchResult.pValue = RemoveExtraneousDigits(strSplitLine[(int)eInspectResultsFileColumns.pvalue]);
                    udtSearchResult.PValueNum = CSngSafe(udtSearchResult.pValue, 0);

                    udtSearchResult.FScore = strSplitLine[(int)eInspectResultsFileColumns.FScore];
                    udtSearchResult.FScoreNum = CSngSafe(udtSearchResult.FScore, 0);

                    udtSearchResult.DeltaScore = strSplitLine[(int)eInspectResultsFileColumns.DeltaScore];
                    udtSearchResult.DeltaScoreOther = strSplitLine[(int)eInspectResultsFileColumns.DeltaScoreOther];

                    udtSearchResult.RecordNumber = strSplitLine[(int)eInspectResultsFileColumns.RecordNumber];
                    udtSearchResult.DBFilePos = strSplitLine[(int)eInspectResultsFileColumns.DBFilePos];
                    udtSearchResult.SpecFilePos = strSplitLine[(int)eInspectResultsFileColumns.SpecFilePos];

                    if (strSplitLine.Length >= (int)eInspectResultsFileColumns.PrecursorError + 1)
                    {
                        // Inspect version 2008-10-14 added these two Precursor mass columns
                        udtSearchResult.PrecursorMZ = strSplitLine[(int)eInspectResultsFileColumns.PrecursorMZ];
                        udtSearchResult.PrecursorError = strSplitLine[(int)eInspectResultsFileColumns.PrecursorError];

                        udtSearchResult.MH = ComputePeptideMHFromPrecursorInfo(udtSearchResult.PrecursorMZ, udtSearchResult.PrecursorError, udtSearchResult.Charge);

                        double dblPrecursorMZ = 0;
                        double dblPrecursorMonoMass = 0;
                        double dblPrecursorErrorDa = 0;

                        double dblPeptideMonoisotopicMass = 0;
                        double dblPeptideDeltaMassCorrectedPpm = 0;

                        if (double.TryParse(udtSearchResult.PrecursorMZ, out dblPrecursorMZ))
                        {
                            dblPrecursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(dblPrecursorMZ, udtSearchResult.ChargeNum, 0);
                            dblPeptideMonoisotopicMass = udtSearchResult.MH - clsPeptideMassCalculator.MASS_PROTON;

                            dblPrecursorErrorDa = dblPrecursorMonoMass - dblPeptideMonoisotopicMass;

                            dblPeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, dblPeptideMonoisotopicMass, true);

                            udtSearchResult.DelMPPM = NumToString(dblPeptideDeltaMassCorrectedPpm, 4, true);
                        }
                    }
                    else
                    {
                        udtSearchResult.PrecursorMZ = "0";
                        udtSearchResult.PrecursorError = "0";
                        udtSearchResult.MH = 0;
                        udtSearchResult.DelMPPM = "0";
                    }

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
                        strErrorLog += "Error parsing InSpecT Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing InSpecT Results in ParseInspectResultsFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        private bool ParseInSpectSynFileEntry(
            ref string strLineIn,
            ref int[] intColumnMapping,
            clsSearchResultsInSpecT objSearchResult,
            ref string strErrorLog,
            ref string strPeptideSequenceWithMods)
        {
            // Parses an entry from the Inspect Synopsis file

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
                if (strSplitLine.Length < INSPECT_SYN_FILE_MIN_COL_COUNT)
                {
                    return false;
                }

                if (!GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.ResultID], out int resultId))
                {
                    ReportError("ResultID column is missing or invalid", true);
                }

                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.Scan], out string scan);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.Charge], out string charge);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.Protein], out string proteinName);

                objSearchResult.ResultID = resultId;
                objSearchResult.Scan = scan;
                objSearchResult.Charge = charge;
                objSearchResult.ProteinName = proteinName;

                objSearchResult.MultipleProteinCount = "0";

                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.Peptide], out strPeptideSequenceWithMods);

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                objSearchResult.SetPeptideSequenceWithMods(strPeptideSequenceWithMods, true, true);

                clsSearchResultsBaseClass objSearchResultBase = default(clsSearchResultsBaseClass);
                objSearchResultBase = (clsSearchResultsBaseClass)objSearchResult;

                base.ComputePseudoPeptideLocInProtein(objSearchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                // If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide
                // will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                objSearchResult.ComputePeptideCleavageStateInProtein();

                if (GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.PrecursorError], out string peptideDeltaMass))
                {
                    objSearchResult.PeptideDeltaMass = peptideDeltaMass;
                    // Note: .peptideDeltaMass is stored in the Inspect results file as "Observed_Mass - Theoretical_Mass"
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
                }
                else
                {
                    objSearchResult.PeptideDeltaMass = "0";
                }

                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.MQScore], out string mqScore);

                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.Length], out string length);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.TotalPRMScore], out string totalPrmScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.MedianPRMScore], out string medianPrmScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.FractionY], out string fractionY);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.FractionB], out string fractionB);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.Intensity], out string intensity);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.NTT], out string ntt);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.PValue], out string pValue);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.FScore], out string fScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.DeltaScore], out string deltaScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.DeltaScoreOther], out string deltaScoreOther);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.DeltaNormMQScore], out string deltaNormMqScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.DeltaNormTotalPRMScore], out string deltaNormTotalPrmScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.RankTotalPRMScore], out string rankTotalPrmScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.RankFScore], out string rankFScore);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.MH], out string peptideMh);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.RecordNumber], out string recordNumber);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.DBFilePos], out string dbFilePos);
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.SpecFilePos], out string specFilePos);

                // Note: .PrecursorError was processed earlier in this function
                GetColumnValue(strSplitLine, intColumnMapping[(int)eInspectSynFileColumns.PrecursorMZ], out string precursorMz);

                objSearchResult.MQScore = mqScore;
                objSearchResult.Length = length;
                objSearchResult.TotalPRMScore = totalPrmScore;
                objSearchResult.MedianPRMScore = medianPrmScore;
                objSearchResult.FractionY = fractionY;
                objSearchResult.FractionB = fractionB;
                objSearchResult.Intensity = intensity;
                objSearchResult.NTT = ntt;
                objSearchResult.pValue = pValue;
                objSearchResult.FScore = fScore;
                objSearchResult.DeltaScore = deltaScore;
                objSearchResult.DeltaScoreOther = deltaScoreOther;
                objSearchResult.DeltaNormMQScore = deltaNormMqScore;
                objSearchResult.DeltaNormTotalPRMScore = deltaNormTotalPrmScore;
                objSearchResult.RankTotalPRMScore = rankTotalPrmScore;
                objSearchResult.RankFScore = rankFScore;
                objSearchResult.PeptideMH = peptideMh;
                objSearchResult.RecordNumber = recordNumber;
                objSearchResult.DBFilePos = dbFilePos;
                objSearchResult.SpecFilePos = specFilePos;
                objSearchResult.PrecursorMz = precursorMz;

                blnValidSearchResult = true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    if ((strSplitLine != null) && strSplitLine.Length > 0)
                    {
                        strErrorLog += "Error parsing InSpecT Results for RowIndex '" + strSplitLine[0] + "'" + "\n";
                    }
                    else
                    {
                        strErrorLog += "Error parsing InSpecT Results in ParseInspectSynFileEntry" + "\n";
                    }
                }
                blnValidSearchResult = false;
            }

            return blnValidSearchResult;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">Inspect results file</param>
        /// <param name="strOutputFolderPath">Output folder</param>
        /// <param name="strParameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            string strOutputFilePath = null;
            string strSynOutputFilePath = null;
            string strPepToProteinMapFilePath = null;

            udtModInfoType[] udtInspectModInfo = null;
            List<udtPepToProteinMappingType> lstPepToProteinMapping = default(List<udtPepToProteinMappingType>);
            string strMTSPepToProteinMapFilePath = string.Empty;

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

                    udtInspectModInfo = new udtModInfoType[-1 + 1];
                    lstPepToProteinMapping = new List<udtPepToProteinMappingType>();

                    // Load the Inspect Parameter File so that we can determine the modification names and masses
                    if (!ExtractModInfoFromInspectParamFile(SearchToolParameterFilePath, ref udtInspectModInfo))
                    {
                        if (udtInspectModInfo == null || udtInspectModInfo.Length == 0)
                        {
                            udtInspectModInfo = new udtModInfoType[1];
                            var modInfo = udtInspectModInfo[0];
                            modInfo.ModName = PHOS_MOD_NAME.ToLower();
                            modInfo.ModMass = PHOS_MOD_MASS;
                            modInfo.Residues = PHOS_MOD_RESIDUES;
                            modInfo.ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL.ToString();
                        }
                    }

                    // Resolve the mods in mInspectModInfo with the ModDefs mods
                    ResolveInspectModsWithModDefinitions(ref udtInspectModInfo);

                    if (base.CreateInspectFirstHitsFile)
                    {
                        // Create the first hits output file
                        base.ResetProgress("Creating the FHT file (top TotalPRMScore)", true);

                        strOutputFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
                        strOutputFilePath = Path.Combine(strOutputFolderPath, strOutputFilePath + INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX);

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strOutputFilePath, ref udtInspectModInfo, eFilteredOutputFileTypeConstants.FHTbyTotalPRM);

                        // Create the first hits output file
                        base.ResetProgress("Creating the FHT file (top FScore)", true);

                        strOutputFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
                        strOutputFilePath = Path.Combine(strOutputFolderPath, strOutputFilePath + INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX);

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strOutputFilePath, ref udtInspectModInfo, eFilteredOutputFileTypeConstants.FHTbyFScore);
                    }

                    if (base.CreateInspectSynopsisFile)
                    {
                        // Create the synopsis output file
                        base.ResetProgress("Creating the SYN file", true);

                        //Define the synopsis output file name based on strInputFilePath
                        strSynOutputFilePath = Path.GetFileNameWithoutExtension(strInputFilePath);
                        strSynOutputFilePath = Path.Combine(strOutputFolderPath, strSynOutputFilePath + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                        blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, ref udtInspectModInfo, eFilteredOutputFileTypeConstants.SynFile);

                        // Load the PeptideToProteinMap information; if the file doesn't exist, then a warning will be displayed, but processing will continue
                        // LoadPeptideToProteinMapInfoInspect also creates _inspect_PepToProtMapMTS.txt file with the new mod symbols and corrected terminii symbols
                        strPepToProteinMapFilePath = Path.Combine(Path.GetDirectoryName(fiInputFile.FullName), Path.GetFileNameWithoutExtension(fiInputFile.FullName) + FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING + ".txt");

                        base.ResetProgress("Loading the PepToProtein map file: " + Path.GetFileName(strPepToProteinMapFilePath), true);

                        LoadPeptideToProteinMapInfoInspect(strPepToProteinMapFilePath, strOutputFolderPath, ref udtInspectModInfo, ref lstPepToProteinMapping, ref strMTSPepToProteinMapFilePath);

                        // Create the other PHRP-specific files
                        base.ResetProgress("Creating the PHRP files for " + Path.GetFileName(strSynOutputFilePath), true);

                        blnSuccess = ParseInSpectSynopsisFile(strSynOutputFilePath, strOutputFolderPath, ref lstPepToProteinMapping, false);

                        // Remove all items from lstPepToProteinMapping to reduce memory overhead
                        lstPepToProteinMapping.Clear();
                        lstPepToProteinMapping.TrimExcess();

                        if (blnSuccess && CreateProteinModsFile)
                        {
                            // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                            base.ValidatePHRPReaderSupportFiles(Path.Combine(fiInputFile.DirectoryName, Path.GetFileName(strSynOutputFilePath)), strOutputFolderPath);

                            // Create the Protein Mods file
                            blnSuccess = base.CreateProteinModDetailsFile(strSynOutputFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.Inspect);
                        }
                    }

                    if (blnSuccess)
                    {
                        base.OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error calling CreateFHTorSYNResultsFile: " + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile:" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return blnSuccess;
        }

        private static readonly Regex RegexNumPlusZeroes = new Regex(@"(\.\d*[1-9])0+$", RegexOptions.Compiled);
        private static readonly Regex RegexAllZeroes = new Regex(@"\.0+$", RegexOptions.Compiled);

        // If strValue ends in .0000, then remove the .0000 portion
        private string RemoveExtraneousDigits(string strValue)
        {
            Match reMatch = default(Match);

            if (string.IsNullOrWhiteSpace(strValue))
            {
                return string.Empty;
            }
            else
            {
                reMatch = RegexAllZeroes.Match(strValue);
                if (reMatch.Success && reMatch.Index > 0)
                {
                    strValue = strValue.Substring(0, reMatch.Index);
                }
                else
                {
                    reMatch = RegexNumPlusZeroes.Match(strValue);
                    if (reMatch.Success && reMatch.Index > 0)
                    {
                        if (reMatch.Groups.Count > 1)
                        {
                            // Number is of the form 1.0030 or 1.300 or 1.030
                            strValue = strValue.Substring(0, reMatch.Index) + reMatch.Groups[1].Value;
                        }
                    }
                }

                return strValue;
            }
        }

        private static readonly Regex RegexNTerminalModMassRegEx = new Regex(INSPECT_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);
        private static readonly Regex RegexCTerminalModMassRegEx = new Regex(INSPECT_CTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS);

        /// <summary>
        /// Replaces modification name text in peptide sequences with modification symbols (uses case-sensitive comparisons)
        /// </summary>
        /// <param name="strPeptide"></param>
        /// <param name="udtInspectModInfo">This function assumes that each entry in udtInspectModInfo() has both .ModName and .ModSymbol defined</param>
        /// <returns></returns>
        /// <remarks></remarks>
        protected string ReplaceInspectModTextWithSymbol(string strPeptide, ref udtModInfoType[] udtInspectModInfo)
        {
            string strPrefix = string.Empty;
            string strSuffix = string.Empty;

            int intModMass = 0;
            Match reMatch = default(Match);

            if (strPeptide.Length >= 4)
            {
                if (strPeptide[1] == '.' && strPeptide[strPeptide.Length - 2] == '.')
                {
                    strPrefix = strPeptide.Substring(0, 2);
                    strSuffix = strPeptide.Substring(strPeptide.Length - 2, 2);

                    strPeptide = strPeptide.Substring(2, strPeptide.Length - 4);
                }
            }

            // strPeptide should now be the clean peptide, without the prefix or suffix residues
            for (var intIndex = 0; intIndex <= udtInspectModInfo.Length - 1; intIndex++)
            {
                if (udtInspectModInfo[intIndex].ModType != eInspectModType.StaticMod)
                {
                    strPeptide = strPeptide.Replace(udtInspectModInfo[intIndex].ModName, udtInspectModInfo[intIndex].ModSymbol);

                    if (udtInspectModInfo[intIndex].ModType == eInspectModType.DynNTermPeptide |
                        udtInspectModInfo[intIndex].ModType == eInspectModType.DynCTermPeptide)
                    {
                        if (udtInspectModInfo[intIndex].ModType == eInspectModType.DynNTermPeptide)
                        {
                            // Inspect notates N-terminal mods like this: R.+14HVIFLAER.R   (Note: This behavior is not yet confirmed)
                            // Look for this using reNTerminalModMassRegEx
                            reMatch = RegexNTerminalModMassRegEx.Match(strPeptide);
                        }
                        else if (udtInspectModInfo[intIndex].ModType == eInspectModType.DynCTermPeptide)
                        {
                            // Inspect notates C-terminal mods like this: R.HVIFLAER+14.R
                            // Look for this using reCTerminalModMassRegEx
                            reMatch = RegexCTerminalModMassRegEx.Match(strPeptide);
                        }
                        else
                        {
                            // This code should never be reached
                            reMatch = null;
                        }

                        if ((reMatch != null))
                        {
                            if (reMatch.Success && reMatch.Groups.Count > 1)
                            {
                                // Match found
                                try
                                {
                                    intModMass = Convert.ToInt32(reMatch.Groups[1].Value);

                                    // Compare the mod mass in the specification to this Mod's mod mass
                                    // If they are less than 0.5 Da apart, then assume we have a match; yes, this assumption is a bit flaky
                                    if (Math.Abs(intModMass - Convert.ToDouble(udtInspectModInfo[intIndex].ModMass)) <= 0.5)
                                    {
                                        // Match found
                                        // Replace the matched region with .ModSymbol

                                        string strPeptideNew = null;

                                        if (reMatch.Groups[0].Index > 0)
                                        {
                                            strPeptideNew = strPeptide.Substring(0, reMatch.Groups[0].Index);
                                        }
                                        else
                                        {
                                            strPeptideNew = string.Empty;
                                        }

                                        strPeptideNew += udtInspectModInfo[intIndex].ModSymbol;

                                        if (reMatch.Groups[0].Index + reMatch.Groups[0].Length < strPeptide.Length)
                                        {
                                            strPeptideNew += strPeptide.Substring(reMatch.Groups[0].Index + reMatch.Groups[0].Length);
                                        }

                                        strPeptide = string.Copy(strPeptideNew);
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

            return strPrefix + strPeptide + strSuffix;
        }

        protected string ReplaceTerminus(string strPeptide)
        {
            if (strPeptide.StartsWith(N_TERMINUS_SYMBOL_INSPECT))
            {
                strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST + "." + strPeptide.Substring(N_TERMINUS_SYMBOL_INSPECT.Length);
            }

            if (strPeptide.EndsWith(C_TERMINUS_SYMBOL_INSPECT))
            {
                strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_INSPECT.Length) + "." + clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
            }

            return strPeptide;
        }

        protected void ResolveInspectModsWithModDefinitions(ref udtModInfoType[] udtInspectModInfo)
        {
            int intResidueIndex = 0;
            int intResIndexStart = 0;
            int intResIndexEnd = 0;

            double dblModMass = 0;
            char chTargetResidue = default(char);
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState = default(clsAminoAcidModInfo.eResidueTerminusStateConstants);
            bool blnExistingModFound = false;

            clsModificationDefinition objModificationDefinition = default(clsModificationDefinition);

            if ((udtInspectModInfo != null))
            {
                // Call .LookupModificationDefinitionByMass for each entry in udtInspectModInfo
                for (var intIndex = 0; intIndex <= udtInspectModInfo.Length - 1; intIndex++)
                {
                    var modInfo = udtInspectModInfo[intIndex];
                    if (double.TryParse(modInfo.ModMass, out dblModMass))
                    {
                        if (modInfo.Residues.Length > 0)
                        {
                            intResIndexStart = 0;
                            intResIndexEnd = modInfo.Residues.Length - 1;
                        }
                        else
                        {
                            intResIndexStart = -1;
                            intResIndexEnd = -1;
                        }

                        for (intResidueIndex = intResIndexStart; intResidueIndex <= intResIndexEnd; intResidueIndex++)
                        {
                            if (intResidueIndex >= 0)
                            {
                                chTargetResidue = modInfo.Residues[intResidueIndex];
                            }
                            else
                            {
                                chTargetResidue = default(char);
                            }

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
                            blnExistingModFound = false;

                            objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(dblModMass, chTargetResidue, eResidueTerminusState, ref blnExistingModFound, true, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION);

                            if (intResidueIndex == intResIndexStart)
                            {
                                modInfo.ModSymbol = objModificationDefinition.ModificationSymbol.ToString();
                            }
                        }
                    }
                }
            }
        }

        protected void StoreOrWriteSearchResult(
            StreamWriter swResultFile,
            ref int intResultID,
            udtInspectSearchResultType udtSearchResult,
            ref int intFilteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string strErrorLog)
        {
            if (SortFHTandSynFiles)
            {
                if (intFilteredSearchResultCount == udtFilteredSearchResults.Length)
                {
                    Array.Resize(ref udtFilteredSearchResults, udtFilteredSearchResults.Length * 2);
                }

                udtFilteredSearchResults[intFilteredSearchResultCount] = udtSearchResult;
                intFilteredSearchResultCount += 1;
            }
            else
            {
                intResultID += 1;
                WriteSearchResultToFile(intResultID, swResultFile, udtSearchResult, ref strErrorLog);
            }
        }

        private void SortAndWriteFilteredSearchResults(
            StreamWriter swResultFile,
            int intFilteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string strErrorLog)
        {
            // Sort udtFilteredSearchResults by descending TotalPRMScore, ascending scan, ascending charge, ascending peptide, and ascending protein
            Array.Sort(udtFilteredSearchResults, 0, intFilteredSearchResultCount, new InspectSearchResultsComparerTotalPRMDescScanChargePeptide());

            for (var intIndex = 0; intIndex <= intFilteredSearchResultCount - 1; intIndex++)
            {
                WriteSearchResultToFile(intIndex + 1, swResultFile, udtFilteredSearchResults[intIndex], ref strErrorLog);
            }
        }

        private void StoreTopFHTMatch(
            StreamWriter swResultFile,
            ref int intResultID,
            int intCurrentScanResultsCount,
            udtInspectSearchResultType[] udtSearchResultsCurrentScan,
            ref int intFilteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string strErrorLog,
            ref IComparer<udtInspectSearchResultType> objSortComparer)
        {
            short intCurrentCharge = short.MinValue;

            AssignRankAndDeltaNormValues(ref udtSearchResultsCurrentScan, intCurrentScanResultsCount);

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, then descending TotalPRMScore or descending FScore (depending on objSortComparer)
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortComparer);

            // Now store or write out the first match for each charge for this scan
            for (var intIndex = 0; intIndex <= intCurrentScanResultsCount - 1; intIndex++)
            {
                if (intIndex == 0 || intCurrentCharge != udtSearchResultsCurrentScan[intIndex].ChargeNum)
                {
                    StoreOrWriteSearchResult(swResultFile, ref intResultID, udtSearchResultsCurrentScan[intIndex], ref intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog);
                    intCurrentCharge = udtSearchResultsCurrentScan[intIndex].ChargeNum;
                }
            }
        }

        private void StoreSynMatches(
            StreamWriter swResultFile,
            ref int intResultID,
            int intCurrentScanResultsCount,
            udtInspectSearchResultType[] udtSearchResultsCurrentScan,
            ref int intFilteredSearchResultCount,
            ref udtInspectSearchResultType[] udtFilteredSearchResults,
            ref string strErrorLog,
            ref IComparer<udtInspectSearchResultType> objSortComparer)
        {
            AssignRankAndDeltaNormValues(ref udtSearchResultsCurrentScan, intCurrentScanResultsCount);

            // Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending FScore
            // All of the data in udtSearchResultsCurrentScan should have the same scan number
            Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortComparer);

            // Now store or write out the matches that pass the filters
            for (var intIndex = 0; intIndex <= intCurrentScanResultsCount - 1; intIndex++)
            {
                if (udtSearchResultsCurrentScan[intIndex].PValueNum <= InspectSynopsisFilePValueThreshold ||
                    udtSearchResultsCurrentScan[intIndex].TotalPRMScoreNum >= TOTALPRMSCORE_THRESHOLD ||
                    udtSearchResultsCurrentScan[intIndex].FScoreNum >= FSCORE_THRESHOLD)
                {
                    StoreOrWriteSearchResult(swResultFile, ref intResultID, udtSearchResultsCurrentScan[intIndex], ref intFilteredSearchResultCount, ref udtFilteredSearchResults, ref strErrorLog);
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
                lstData.Add(COLUMN_NAME_RESULTID);
                lstData.Add("Scan");
                lstData.Add(COLUMN_NAME_PEPTIDE);
                lstData.Add("Protein");
                lstData.Add("Charge");
                lstData.Add("MQScore");
                lstData.Add("Length");
                lstData.Add("TotalPRMScore");
                lstData.Add("MedianPRMScore");
                lstData.Add("FractionY");
                lstData.Add("FractionB");
                lstData.Add("Intensity");
                lstData.Add("NTT");
                lstData.Add("PValue");
                lstData.Add("FScore");
                lstData.Add("DeltaScore");
                lstData.Add("DeltaScoreOther");
                lstData.Add("DeltaNormMQScore");
                lstData.Add("DeltaNormTotalPRMScore");
                lstData.Add("RankTotalPRMScore");
                lstData.Add("RankFScore");
                lstData.Add("MH");
                lstData.Add("RecordNumber");
                lstData.Add("DBFilePos");
                lstData.Add("SpecFilePos");
                lstData.Add("PrecursorMZ");
                lstData.Add("PrecursorError");
                lstData.Add("DelM_PPM");

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

        private void WriteSearchResultToFile(
            int intResultID,
            StreamWriter swResultFile,
            udtInspectSearchResultType udtSearchResult,
            ref string strErrorLog)
        {
            // Writes an entry to a synopsis or first hits file
            try
            {
                swResultFile.WriteLine(intResultID.ToString() + "\t" +
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
                                       NumToString(udtSearchResult.DeltaNormMQScore, 5, true) + "\t" +
                                       NumToString(udtSearchResult.DeltaNormTotalPRMScore, 5, true) + "\t" +
                                       udtSearchResult.RankTotalPRMScore + "\t" +
                                       udtSearchResult.RankFScore + "\t" +
                                       NumToString(udtSearchResult.MH, 6, true) + "\t" +
                                       udtSearchResult.RecordNumber + "\t" +
                                       udtSearchResult.DBFilePos + "\t" +
                                       udtSearchResult.SpecFilePos + "\t" +
                                       udtSearchResult.PrecursorMZ + "\t" +
                                       udtSearchResult.PrecursorError + "\t" +
                                       udtSearchResult.DelMPPM);
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

        protected class InspectSearchResultsComparerTotalPRMDescScanChargePeptide : IComparer<udtInspectSearchResultType>
        {
            public int Compare(udtInspectSearchResultType x, udtInspectSearchResultType y)
            {
                if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                {
                    return -1;
                }
                else if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                {
                    return 1;
                }
                else
                {
                    // TotalPRMScore is the same; check scan number
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
                            var result = x.PeptideAnnotation.CompareTo(y.PeptideAnnotation);
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

        protected class InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc : IComparer<udtInspectSearchResultType>
        {
            public int Compare(udtInspectSearchResultType x, udtInspectSearchResultType y)
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
                        // Charge is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                        if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                        {
                            return -1;
                        }
                        else if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                        {
                            return 1;
                        }
                        else
                        {
                            // TotalPRMScore is the same; check FScore (sort on descending FScore)
                            if (x.FScoreNum > y.FScoreNum)
                            {
                                return -1;
                            }
                            else if (x.FScoreNum < y.FScoreNum)
                            {
                                return 1;
                            }
                            else
                            {
                                return 0;
                            }
                        }
                    }
                }
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
                else if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }
                else
                {
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
                        // Charge is the same; check FScore (sort on descending FScore)
                        if (x.FScoreNum > y.FScoreNum)
                        {
                            return -1;
                        }
                        else if (x.FScoreNum < y.FScoreNum)
                        {
                            return 1;
                        }
                        else
                        {
                            // FScore is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                            if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                            {
                                return -1;
                            }
                            else if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                            {
                                return 1;
                            }
                            else
                            {
                                return 0;
                            }
                        }
                    }
                }
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
                else if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }
                else
                {
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
                        // Charge is the same; check MQScore (sort on descending MQScore)
                        if (x.MQScoreNum > y.MQScoreNum)
                        {
                            return -1;
                        }
                        else if (x.MQScoreNum < y.MQScoreNum)
                        {
                            return 1;
                        }
                        else
                        {
                            // MQScore is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                            if (x.TotalPRMScoreNum > y.TotalPRMScoreNum)
                            {
                                return -1;
                            }
                            else if (x.TotalPRMScoreNum < y.TotalPRMScoreNum)
                            {
                                return 1;
                            }
                            else
                            {
                                return 0;
                            }
                        }
                    }
                }
            }
        }

        #endregion
    }
}
