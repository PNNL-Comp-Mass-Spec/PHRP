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

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for MSPathfinder
    /// </summary>
    public class clsPHRPParserMSPathFinder : clsPHRPParser
    {
        #region "Constants"

#pragma warning disable 1591
        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_MostAbundantIsotopeMz = "MostAbundantIsotopeMz";
        public const string DATA_COLUMN_Mass = "Mass";
        public const string DATA_COLUMN_Sequence = "Sequence";
        public const string DATA_COLUMN_Modifications = "Modifications";
        public const string DATA_COLUMN_Composition = "Composition";
        public const string DATA_COLUMN_Protein = "ProteinName";
        public const string DATA_COLUMN_ProteinDesc = "ProteinDesc";
        public const string DATA_COLUMN_ProteinLength = "ProteinLength";
        public const string DATA_COLUMN_ResidueStart = "ResidueStart";
        public const string DATA_COLUMN_ResidueEnd = "ResidueEnd";
        public const string DATA_COLUMN_MatchedFragments = "MatchedFragments";
        public const string DATA_COLUMN_SpecEValue = "SpecEValue";
        public const string DATA_COLUMN_EValue = "EValue";
        public const string DATA_COLUMN_QValue = "QValue";
        public const string DATA_COLUMN_PepQValue = "PepQValue";

        public const string FILENAME_SUFFIX_SYN = "_mspath_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_mspath_fht.txt";

        private const string MSPathFinder_SEARCH_ENGINE_NAME = "MSPathFinder";
#pragma warning restore 1591

        #endregion

        #region "Properties"

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

        #endregion

        /// <summary>
        /// Constructor; assumes loadModsAndSeqInfo=True
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <remarks></remarks>
        public clsPHRPParserMSPathFinder(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, then load the ModSummary file and SeqInfo files</param>
        /// <remarks></remarks>
        public clsPHRPParserMSPathFinder(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MSPathFinder, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserMSPathFinder(string datasetName, string inputFilePath, clsPHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.MSPathFinder, startupOptions)
        {
        }

        /// <summary>
        /// Define column header names
        /// </summary>
        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_ResultID);

            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_MostAbundantIsotopeMz);
            AddHeaderColumn(DATA_COLUMN_Mass);
            AddHeaderColumn(DATA_COLUMN_Sequence);
            AddHeaderColumn(DATA_COLUMN_Modifications);
            AddHeaderColumn(DATA_COLUMN_Composition);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_ProteinDesc);
            AddHeaderColumn(DATA_COLUMN_ProteinLength);
            AddHeaderColumn(DATA_COLUMN_ResidueStart);
            AddHeaderColumn(DATA_COLUMN_ResidueEnd);
            AddHeaderColumn(DATA_COLUMN_MatchedFragments);
            AddHeaderColumn(DATA_COLUMN_SpecEValue);
            AddHeaderColumn(DATA_COLUMN_EValue);
            AddHeaderColumn(DATA_COLUMN_QValue);
            AddHeaderColumn(DATA_COLUMN_PepQValue);
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // MSPathFinder does not have a first-hits file; just the _syn.txt file
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
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out clsSearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new clsSearchEngineParameters(MSPathFinder_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams, clsPHRPReader.ePeptideHitResultType.MSPathFinder);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, clsSearchEngineParameters searchEngineParams, clsPHRPReader.ePeptideHitResultType resultType)
        {
            try
            {
                var success = ReadKeyValuePairSearchEngineParamFile(MSPathFinder_SEARCH_ENGINE_NAME, searchEngineParamFileName, clsPHRPReader.ePeptideHitResultType.MSGFPlus, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                searchEngineParams.Enzyme = "no_enzyme";
                searchEngineParams.MinNumberTermini = 0;

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                searchEngineParams.PrecursorMassToleranceDa = clsPHRPParserMSGFDB.DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM, resultType);
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
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">clsPSM object (output)</param>
        /// <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out clsPSM psm, bool fastReadMode)
        {
            psm = new clsPSM();

            try
            {
                var columns = line.Split('\t');
                var success = false;

                psm.DataLineText = line;
                psm.ScanNumber = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, -100);
                if (psm.ScanNumber == -100)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_ResultID, mColumnHeaders, 0);
                    psm.ScoreRank = 1;

                    var sequence = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Sequence, mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(sequence, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(sequence, mCleavageStateCalculator);
                    }

                    psm.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                    var protein = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Protein, mColumnHeaders);
                    psm.AddProtein(protein);

                    // Store the sequence mass as the "precursor" mass, though MSPathFinderT results are from MS1 spectra, and thus we didn't do MS/MS on a precursor
                    psm.PrecursorNeutralMass = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Mass, mColumnHeaders, 0.0);

                    // Thus collision mode, precursor neutral mass, etc. are not applicable
                    // psm.CollisionMode =
                    // psm.MassErrorDa =
                    // psm.MassErrorPPM =

                    // psm.MSGFSpecEValue =

                    success = true;
                }

                if (!success)
                {
                    return false;
                }

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                // Store the remaining data

                AddScore(psm, columns, DATA_COLUMN_MostAbundantIsotopeMz);
                AddScore(psm, columns, DATA_COLUMN_Modifications);
                AddScore(psm, columns, DATA_COLUMN_Composition);
                AddScore(psm, columns, DATA_COLUMN_ProteinDesc);
                AddScore(psm, columns, DATA_COLUMN_ProteinLength);
                AddScore(psm, columns, DATA_COLUMN_ResidueStart);
                AddScore(psm, columns, DATA_COLUMN_ResidueEnd);
                AddScore(psm, columns, DATA_COLUMN_MatchedFragments);
                AddScore(psm, columns, DATA_COLUMN_SpecEValue);
                AddScore(psm, columns, DATA_COLUMN_EValue);
                AddScore(psm, columns, DATA_COLUMN_QValue);
                AddScore(psm, columns, DATA_COLUMN_PepQValue);

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
