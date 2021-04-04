using PHRPReader.Data;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for MaxQuant
    /// </summary>
    public class MaxQuantSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: maxq, MaxQuant

#pragma warning disable 1591

        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Spectrum_Index = "Spectrum_Index";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Score = "Score";
        public const string DATA_COLUMN_Probability = "Probability";
        public const string DATA_COLUMN_Rank_Probability = "Rank_Probability";
        public const string DATA_COLUMN_Peptide_Position = "Peptide_Position";
        public const string DATA_COLUMN_QValue = "QValue";

        public const string FILENAME_SUFFIX_SYN = "_maxq_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_maxq_fht.txt";

        private const string MAXQUANT_SEARCH_ENGINE_NAME = "MaxQuant";

#pragma warning restore 1591

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
        public MaxQuantSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MaxQuantSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MaxQuant, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MaxQuantSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MaxQuant, startupOptions)
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
        public static SortedDictionary<string, MaxQuantSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            var headerColumns = new SortedDictionary<string, MaxQuantSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {DATA_COLUMN_ResultID, MaxQuantSynFileColumns.ResultID},
                {DATA_COLUMN_Scan, MaxQuantSynFileColumns.Scan},
                //{DATA_COLUMN_Spectrum_Index, MaxQuantSynFileColumns.Spectrum_Index},
                //{DATA_COLUMN_Charge, MaxQuantSynFileColumns.Charge},
                //{DATA_COLUMN_PrecursorMZ, MaxQuantSynFileColumns.PrecursorMZ},
                //{DATA_COLUMN_DelM, MaxQuantSynFileColumns.DelM},
                //{DATA_COLUMN_DelM_PPM, MaxQuantSynFileColumns.DelM_PPM},
                //{DATA_COLUMN_MH, MaxQuantSynFileColumns.MH},
                //{DATA_COLUMN_Peptide, MaxQuantSynFileColumns.Peptide},
                //{DATA_COLUMN_Protein, MaxQuantSynFileColumns.Protein},
                //{DATA_COLUMN_Score, MaxQuantSynFileColumns.Score},
                //{DATA_COLUMN_Probability, MaxQuantSynFileColumns.Probability},
                //{DATA_COLUMN_Rank_Probability, MaxQuantSynFileColumns.Rank_Probability},
                //{DATA_COLUMN_Peptide_Position, MaxQuantSynFileColumns.Peptide_Position},
                //{DATA_COLUMN_QValue, MaxQuantSynFileColumns.QValue}
            };

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MaxQuantSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MaxQuantSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // MaxQuant does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_maxq_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_maxq_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_maxq_syn_ProteinMods.txt";
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
            return datasetName + "_maxq_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_maxq_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_maxq_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MAXQUANT_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MaxQuant parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MAXQUANT_SEARCH_ENGINE_NAME);

            //var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            // ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            var success = false;
            return success;
        }

        /// <summary>
        /// Parse the data line read from a PHRP results file
        /// </summary>
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM details</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode)
        {
            throw new NotImplementedException();
        }
    }
}
