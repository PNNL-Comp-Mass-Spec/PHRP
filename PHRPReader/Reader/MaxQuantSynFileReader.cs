using PHRPReader.Data;
using System;
using System.Collections.Generic;
using System.IO;
using System.Xml.Linq;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for MaxQuant
    /// </summary>
    public class MaxQuantSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: Chymotrypsin, maxq, MaxQuant, PepToProtMap, tryptic

        /// <summary>
        /// MaxQuant synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_maxq_syn.txt";

        /// <summary>
        /// MaxQuant first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for MaxQuant results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_maxq_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        // ReSharper disable once IdentifierTypo
        private const string MAXQUANT_SEARCH_ENGINE_NAME = "MaxQuant";

        /// <summary>
        /// Mapping from enum to synopsis file column name for MaxQuant
        /// </summary>
        private static readonly Dictionary<MaxQuantSynFileColumns, string> mSynopsisFileColumn = new();

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
            return new(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", MaxQuantSynFileColumns.ResultID },
                { "Dataset", MaxQuantSynFileColumns.Dataset },
                { "DatasetID", MaxQuantSynFileColumns.DatasetID },
                { "Scan", MaxQuantSynFileColumns.Scan },
                { "FragMethod", MaxQuantSynFileColumns.FragMethod },
                { "SpecIndex", MaxQuantSynFileColumns.SpecIndex },
                { "Charge", MaxQuantSynFileColumns.Charge },
                { "PrecursorMZ", MaxQuantSynFileColumns.PrecursorMZ },
                { "PrecursorMZ_MaxQuant", MaxQuantSynFileColumns.PrecursorMZ_MaxQuant },
                { "DelM", MaxQuantSynFileColumns.DelM },
                { "DelM_PPM", MaxQuantSynFileColumns.DelM_PPM },
                { "DelM_MaxQuant", MaxQuantSynFileColumns.DelM_MaxQuant },
                { "DelM_PPM_MaxQuant", MaxQuantSynFileColumns.DelM_PPM_MaxQuant },
                { "MH", MaxQuantSynFileColumns.MH },
                { "Mass", MaxQuantSynFileColumns.Mass },
                { "Peptide", MaxQuantSynFileColumns.Peptide },
                { "Modifications", MaxQuantSynFileColumns.DynamicModifications },
                { "Proteins", MaxQuantSynFileColumns.Proteins },
                { "LeadingRazorProtein", MaxQuantSynFileColumns.LeadingRazorProtein },
                { "NTT", MaxQuantSynFileColumns.NTT },
                { "PEP", MaxQuantSynFileColumns.PEP },
                { "Score", MaxQuantSynFileColumns.Score },
                { "DeltaScore", MaxQuantSynFileColumns.DeltaScore },
                { "RankScore", MaxQuantSynFileColumns.RankScore },
                { "TotalPeptideIntensity", MaxQuantSynFileColumns.TotalPeptideIntensity },
                { "MassAnalyzer", MaxQuantSynFileColumns.MassAnalyzer },
                { "PrecursorType", MaxQuantSynFileColumns.PrecursorType },
                { "RetentionTime", MaxQuantSynFileColumns.RetentionTime },
                { "PrecursorScan", MaxQuantSynFileColumns.PrecursorScan },
                { "PrecursorIntensity", MaxQuantSynFileColumns.PrecursorIntensity },
                { "NumberOfMatches", MaxQuantSynFileColumns.NumberOfMatches },
                { "IntensityCoverage", MaxQuantSynFileColumns.IntensityCoverage },
                { "MissedCleavages", MaxQuantSynFileColumns.MissedCleavages },
                { "MsMsID", MaxQuantSynFileColumns.MsMsID },
                { "ProteinGroupIDs", MaxQuantSynFileColumns.ProteinGroupIDs },
                { "PeptideID", MaxQuantSynFileColumns.PeptideID },
                { "ModPeptideID", MaxQuantSynFileColumns.ModPeptideID },
                { "EvidenceID", MaxQuantSynFileColumns.EvidenceID },
                { "QValue", MaxQuantSynFileColumns.QValue }
            };
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
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MaxQuantSynFileColumns column)
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

            return ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);
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

            var columns = line.Split('\t');

            var success = false;

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.ResultID), mColumnHeaders, 0);
                    psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.RankScore), mColumnHeaders, 0);

                    var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.Peptide), mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.Charge), mColumnHeaders, 0);

                    var proteinNames = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.Proteins), mColumnHeaders);
                    if (!string.IsNullOrWhiteSpace(proteinNames))
                    {
                        foreach (var protein in proteinNames.Split(';'))
                        {
                            if (string.IsNullOrWhiteSpace(protein))
                                continue;

                            psm.AddProtein(protein.Trim());
                        }
                    }

                    var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);

                    if (Math.Abs(precursorMZ) > float.Epsilon)
                    {
                        psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);
                    }
                    else
                    {
                        var precursorMZ_MaxQuant = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.PrecursorMZ_MaxQuant), mColumnHeaders, 0.0);
                        psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ_MaxQuant, psm.Charge, 0);
                    }

                    psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.DelM), mColumnHeaders);
                    psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MaxQuantSynFileColumns.DelM_PPM), mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining data
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.Dataset));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.DatasetID));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.FragMethod));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.SpecIndex));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.DelM_MaxQuant));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.DelM_PPM_MaxQuant));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.MH));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.Mass));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.DynamicModifications));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.LeadingRazorProtein));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.NTT));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.PEP));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.Score));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.DeltaScore));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.TotalPeptideIntensity));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.MassAnalyzer));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.PrecursorType));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.RetentionTime));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.PrecursorScan));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.PrecursorIntensity));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.NumberOfMatches));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.IntensityCoverage));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.MissedCleavages));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.MsMsID));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.ProteinGroupIDs));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.PeptideID));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.ModPeptideID));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.EvidenceID));
                    AddScore(psm, columns, GetColumnNameByID(MaxQuantSynFileColumns.QValue));

                    if (psm.Proteins.Count == 0)
                    {
                        // Add the leading razor protein to the protein list
                        if (psm.TryGetScore(GetColumnNameByID(MaxQuantSynFileColumns.LeadingRazorProtein), out var leadingRazorProtein) &&
                            !string.IsNullOrWhiteSpace(leadingRazorProtein))
                        {
                            psm.AddProtein(leadingRazorProtein.Trim());
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MSAlign data file: " + ex.Message);
            }

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);
                var paramFile = new FileInfo(paramFilePath);

                if (!paramFile.Exists)
                {
                    ReportWarning("MaxQuant param file not found: " + paramFilePath);
                    return false;
                }

                searchEngineParams.UpdateSearchEngineParamFilePath(paramFilePath);

                using var reader = new StreamReader(new FileStream(paramFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                // Note that XDocument supersedes XmlDocument and XPathDocument
                // XDocument can often be easier to use since XDocument is LINQ-based

                var doc = XDocument.Parse(reader.ReadToEnd());

                // Determine the enzyme name

                var enzymeName = XmlReaderUtilities.GetElementValueOrDefault(
                    doc.Root,
                    "parameterGroups", "parameterGroup", "enzymes", "string");

                if (string.IsNullOrWhiteSpace(enzymeName))
                {
                    ReportWarning("'enzymes' node not found in the MaxQuant parameter file");
                }
                else
                {
                    searchEngineParams.Enzyme = enzymeName;

                    switch (searchEngineParams.Enzyme)
                    {
                        case "Trypsin":
                        case "Trypsin/P":
                        case "LysC":
                        case "LysC/P":
                        case "D.P":
                        case "ArgC":
                        case "AspC":
                        case "GluC":
                        case "GluN":
                        case "AspN":
                        case "LysN":
                        case "Chymotrypsin+":
                        case "Chymotrypsin":
                            break;

                        default:
                            ReportWarning(string.Format("Unrecognized enzyme '{0}' in the MaxQuant parameter file", searchEngineParams.Enzyme));
                            break;
                    }
                }

                // Determine the cleavage specificity
                var enzymeMode = XmlReaderUtilities.GetElementValueOrDefault(
                    doc.Root,
                    "parameterGroups", "parameterGroup", "enzymeMode");

                if (string.IsNullOrWhiteSpace(enzymeMode))
                {
                    ReportWarning("'enzymeMode' node not found in the MaxQuant parameter file");
                }
                else
                {
                    switch (enzymeMode)
                    {
                        case "0":
                            // Fully tryptic
                            searchEngineParams.MinNumberTermini = 2;
                            break;

                        case "1":
                            // Partially tryptic, free N-terminus
                            searchEngineParams.MinNumberTermini = 1;
                            break;

                        case "2":
                            // Partially tryptic, free C-terminus
                            searchEngineParams.MinNumberTermini = 1;
                            break;

                        case "3":
                            // Partially tryptic
                            searchEngineParams.MinNumberTermini = 1;
                            break;

                        case "4":
                            // Non-tryptic
                            searchEngineParams.MinNumberTermini = 0;
                            break;

                        default:
                            ReportWarning(string.Format("Unrecognized enzymeMode '{0}'' in the MaxQuant parameter file", enzymeMode));
                            break;
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                //var firstSearchPrecursorTolerancePPM = XmlReaderUtilities.GetElementValueOrDefault(
                //    doc.Root,
                //    "parameterGroups", "parameterGroup", "firstSearchTol");

                var mainSearchPrecursorTolerancePPM = XmlReaderUtilities.GetElementValueOrDefault(
                    doc.Root,
                    "parameterGroups", "parameterGroup", "mainSearchTol");

                if (double.TryParse(mainSearchPrecursorTolerancePPM, out var tolerancePPM))
                {
                    searchEngineParams.PrecursorMassToleranceDa = PeptideMassCalculator.PPMToMass(tolerancePPM, 2000);
                    searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;
                }
                else
                {
                    searchEngineParams.PrecursorMassToleranceDa = 0;
                    searchEngineParams.PrecursorMassTolerancePpm = 0;
                }

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }
    }
}
