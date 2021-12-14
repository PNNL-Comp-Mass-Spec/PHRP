//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 07/17/2015
//
// This class parses data lines from mspath_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for MSPathfinder
    /// </summary>
    public class MSPathFinderSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: mspath, ProteinDesc, PepToProt

        /// <summary>
        /// MSPathFinder synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_mspath_syn.txt";

        /// <summary>
        /// MSPathFinder first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for MaxQuant results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_mspath_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string MSPathFinder_SEARCH_ENGINE_NAME = "MSPathFinder";

        /// <summary>
        /// Mapping from enum to synopsis file column name for MSPathFinder
        /// </summary>
        private static readonly Dictionary<MSPathFinderSynFileColumns, string> mSynopsisFileColumn = new();

        /// <summary>
        /// First hits file
        /// </summary>
        public override string PHRPFirstHitsFileName => GetPHRPFirstHitsFileName(mDatasetName);

        /// <summary>
        /// Mod summary file
        /// </summary>
        public override string PHRPModSummaryFileName => GetPHRPModSummaryFileName(mDatasetName);

        /// <summary>
        /// Peptide to protein map file
        /// </summary>
        public override string PHRPPepToProteinMapFileName => GetPHRPPepToProteinMapFileName(mDatasetName);

        /// <summary>
        /// Protein mods file
        /// </summary>
        public override string PHRPProteinModsFileName => GetPHRPProteinModsFileName(mDatasetName);

        /// <summary>
        /// Synopsis file
        /// </summary>
        public override string PHRPSynopsisFileName => GetPHRPSynopsisFileName(mDatasetName);

        /// <summary>
        /// Result to sequence map file
        /// </summary>
        public override string PHRPResultToSeqMapFileName => GetPHRPResultToSeqMapFileName(mDatasetName);

        /// <summary>
        /// Sequence info file
        /// </summary>
        public override string PHRPSeqInfoFileName => GetPHRPSeqInfoFileName(mDatasetName);

        /// <summary>
        /// Sequence to protein map file
        /// </summary>
        public override string PHRPSeqToProteinMapFileName => GetPHRPSeqToProteinMapFileName(mDatasetName);

        /// <summary>
        /// Search engine name
        /// </summary>
        public override string SearchEngineName => GetSearchEngineName();

        /// <summary>
        /// Constructor; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        public MSPathFinderSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MSPathFinderSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MSPathFinder, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MSPathFinderSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MSPathFinder, startupOptions)
        {
        }

        /// <summary>
        /// Get the header names in the PHRP synopsis or first hits file for this tool
        /// </summary>
        /// <returns>List of header names</returns>
        protected override List<string> GetColumnHeaderNames()
        {
            var headerNames = new List<string>();
            headerNames.AddRange(GetColumnHeaderNamesAndIDs().Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, MSPathFinderSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, MSPathFinderSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", MSPathFinderSynFileColumns.ResultID },
                { "Scan", MSPathFinderSynFileColumns.Scan },
                { "Charge", MSPathFinderSynFileColumns.Charge },
                { "MostAbundantIsotopeMz", MSPathFinderSynFileColumns.MostAbundantIsotopeMz },
                { "Mass", MSPathFinderSynFileColumns.Mass },
                { "Sequence", MSPathFinderSynFileColumns.Sequence },
                { "Modifications", MSPathFinderSynFileColumns.Modifications },
                { "Composition", MSPathFinderSynFileColumns.Composition },
                { "ProteinName", MSPathFinderSynFileColumns.Protein },
                { "ProteinDesc", MSPathFinderSynFileColumns.ProteinDesc },
                { "ProteinLength", MSPathFinderSynFileColumns.ProteinLength },
                { "ResidueStart", MSPathFinderSynFileColumns.ResidueStart },
                { "ResidueEnd", MSPathFinderSynFileColumns.ResidueEnd },
                { "MatchedFragments", MSPathFinderSynFileColumns.MatchedFragments },
                { "SpecEValue", MSPathFinderSynFileColumns.SpecEValue },
                { "EValue", MSPathFinderSynFileColumns.EValue },
                { "QValue", MSPathFinderSynFileColumns.QValue },
                { "PepQValue", MSPathFinderSynFileColumns.PepQValue },
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSPathFinderSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MSPathFinderSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MSPathFinderSynFileColumns column)
        {
            if (mSynopsisFileColumn.Count > 0)
            {
                return mSynopsisFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mSynopsisFileColumn.Add(item.Value, item.Key);
            }

            return mSynopsisFileColumn[column];
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Empty string, since MSPathFinder does not have a first-hits file; just the _syn.txt file</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_mspath_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_mspath_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_mspath_syn_ProteinMods.txt";
        }

        /// <summary>
        /// Default Synopsis file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_SYN;
        }

        /// <summary>
        /// Default ResultToSeq map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(string datasetName)
        {
            return datasetName + "_mspath_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_mspath_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_mspath_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MSPathFinder_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MSPathFinder parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MSPathFinder_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                const PeptideHitResultTypes resultType = PeptideHitResultTypes.MSPathFinder;
                var success = ReadKeyValuePairSearchEngineParamFile(MSPathFinder_SEARCH_ENGINE_NAME, searchEngineParamFileName, resultType, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                searchEngineParams.Enzyme = "no_enzyme";
                searchEngineParams.MinNumberTermini = 0;

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                searchEngineParams.PrecursorMassToleranceDa = MSGFPlusSynFileReader.DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM, resultType);
                searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM details</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;

            psm = new PSM();

            try
            {
                var columns = line.Split('\t');

                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSPathFinderSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSPathFinderSynFileColumns.ResultID), mColumnHeaders, 0);
                psm.ScoreRank = 1;

                var sequence = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSPathFinderSynFileColumns.Sequence), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(sequence, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(sequence, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSPathFinderSynFileColumns.Charge), mColumnHeaders, 0);

                var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSPathFinderSynFileColumns.Protein), mColumnHeaders);
                psm.AddProtein(protein);

                // Store the sequence mass as the "precursor" mass, though MSPathFinderT results are from MS1 spectra, and thus we didn't do MS/MS on a precursor
                psm.PrecursorNeutralMass = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSPathFinderSynFileColumns.Mass), mColumnHeaders, 0.0);

                // Collision mode, precursor neutral mass, etc. are not applicable
                // psm.CollisionMode =
                // psm.MassErrorDa =
                // psm.MassErrorPPM =

                // psm.MSGFSpecEValue =

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                // Store the remaining data

                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.MostAbundantIsotopeMz));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.Modifications));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.Composition));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.ProteinDesc));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.ProteinLength));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.ResidueStart));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.ResidueEnd));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.MatchedFragments));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.SpecEValue));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.EValue));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.QValue));
                AddScore(psm, columns, GetColumnNameByID(MSPathFinderSynFileColumns.PepQValue));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MSGFDB data file: " + ex.Message);
                return false;
            }
        }
    }
}
