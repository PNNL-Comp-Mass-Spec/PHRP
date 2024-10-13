using PHRPReader.Data;
using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for MSFragger
    /// </summary>
    public class MSFraggerSynFileReader : SynFileReaderBaseClass
    {
        // ReSharper disable CommentTypo

        // Ignore Spelling: Chymotrypsin, Daltons, FHT, Fragger, Hyperscore, MSFragger, NextScore, PepToProtMap, psm, tryptic
        // Ignore Spelling: argc, aspn, cnbr, elastase, formicacid, gluc, lysc, lysn, stricttrypsin, thermolysin

        // ReSharper restore CommentTypo

        /// <summary>
        /// MSFragger synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_msfragger_syn.txt";

        /// <summary>
        /// MSFragger first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for MSFragger results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_msfragger_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        // ReSharper disable once IdentifierTypo
        private const string MSFRAGGER_SEARCH_ENGINE_NAME = "MSFragger";

        /// <summary>
        /// Mapping from enum to synopsis file column name for MSFragger
        /// </summary>
        private static readonly Dictionary<MSFraggerSynFileColumns, string> mSynopsisFileColumn = new();

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
        public MSFraggerSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MSFraggerSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MSFragger, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MSFraggerSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MSFragger, startupOptions)
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
        public static SortedDictionary<string, MSFraggerSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, MSFraggerSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", MSFraggerSynFileColumns.ResultID },
                { "Dataset", MSFraggerSynFileColumns.Dataset },
                { "DatasetID", MSFraggerSynFileColumns.DatasetID },
                { "Scan", MSFraggerSynFileColumns.Scan },
                { "Charge", MSFraggerSynFileColumns.Charge },
                { "PrecursorMZ", MSFraggerSynFileColumns.PrecursorMZ },
                { "DelM", MSFraggerSynFileColumns.DelM },
                { "DelM_PPM", MSFraggerSynFileColumns.DelM_PPM },
                { "DelM_MSFragger", MSFraggerSynFileColumns.DelM_MSFragger },
                { "MH", MSFraggerSynFileColumns.MH },
                { "Mass", MSFraggerSynFileColumns.Mass },
                { "Peptide", MSFraggerSynFileColumns.Peptide },
                { "Modifications", MSFraggerSynFileColumns.Modifications },
                { "Protein", MSFraggerSynFileColumns.Protein },
                { "AdditionalProteins", MSFraggerSynFileColumns.AdditionalProteins },
                { "NTT", MSFraggerSynFileColumns.NTT },
                { "EValue", MSFraggerSynFileColumns.EValue },
                { "Rank_EValue", MSFraggerSynFileColumns.RankEValue },
                { "Hyperscore", MSFraggerSynFileColumns.Hyperscore },
                { "Nextscore", MSFraggerSynFileColumns.Nextscore },
                { "PeptideProphetProbability", MSFraggerSynFileColumns.PeptideProphetProbability },
                { "ElutionTime", MSFraggerSynFileColumns.ElutionTime },
                { "ElutionTimeAverage", MSFraggerSynFileColumns.ElutionTimeAverage },
                { "MissedCleavages", MSFraggerSynFileColumns.MissedCleavages },
                { "MatchedIons", MSFraggerSynFileColumns.NumberOfMatchedIons },
                { "TotalIons", MSFraggerSynFileColumns.TotalNumberOfIons },
                { "QValue", MSFraggerSynFileColumns.QValue }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSFraggerSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames">List of header names</param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MSFraggerSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column">Column enum</param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MSFraggerSynFileColumns column)
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

        private bool GetParameterValue(
            SearchEngineParameters searchEngineParams,
            string parameterPrefix,
            string parameterName,
            out string parameterValue)
        {
            return searchEngineParams.Parameters.TryGetValue(string.Format("{0}{1}", parameterPrefix, parameterName), out parameterValue);
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // MSFragger does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_msfragger_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msfragger_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_msfragger_syn_ProteinMods.txt";
        }

        /// <summary>
        /// Default Synopsis file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSynopsisFileName(string datasetName)
        {
            return datasetName + FILENAME_SUFFIX_SYN;
        }

        /// <summary>
        /// Default ResultToSeq map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPResultToSeqMapFileName(string datasetName)
        {
            return datasetName + "_msfragger_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_msfragger_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msfragger_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Examine the MSFragger or FragPipe parameters to determine the precursor mass tolerance(s)
        /// </summary>
        /// <param name="searchEngineParams">MSFragger search engine parameters loaded from a Key=Value parameter file</param>
        /// <param name="parameterPrefix">Parameter name prefix: empty string for MSFragger, "msfragger." for FragPipe</param>
        /// <param name="toleranceLower">Output: Tolerance to the left, e.g. -20</param>
        /// <param name="toleranceUpper">Output: Tolerance to the right, e.g. 20</param>
        /// <param name="ppmBased">Output: True if ppm-based tolerances</param>
        /// <param name="singleTolerance">Output: true if a single tolerance is defined, false if two tolerances are defined</param>
        /// <returns>True if the tolerance parameters were found, false if not found or an error</returns>
        public bool GetPrecursorSearchTolerances(
            SearchEngineParameters searchEngineParams,
            string parameterPrefix,
            out double toleranceLower,
            out double toleranceUpper,
            out bool ppmBased,
            out bool singleTolerance)
        {
            // First look for the legacy symmetric tolerance parameter
            GetParameterValue(searchEngineParams, parameterPrefix, "precursor_true_tolerance", out var precursorTrueTolerance);
            GetParameterValue(searchEngineParams, parameterPrefix, "precursor_true_units", out var precursorTrueUnits);

            // Next look for the newer tolerance parameters
            GetParameterValue(searchEngineParams, parameterPrefix, "precursor_mass_lower", out var precursorMassLower);
            GetParameterValue(searchEngineParams, parameterPrefix, "precursor_mass_upper", out var precursorMassUpper);
            GetParameterValue(searchEngineParams, parameterPrefix, "precursor_mass_units", out var precursorMassUnits);

            if (int.TryParse(precursorMassUnits, out var precursorMassUnitsVal) &&
                double.TryParse(precursorMassLower, out var precursorMassLowerVal) &&
                double.TryParse(precursorMassUpper, out var precursorMassUpperVal))
            {
                // Two tolerances are defined (though they may be equivalent)
                ppmBased = MassToleranceUnitsArePPM(precursorMassUnitsVal, "precursor_mass_units");
                singleTolerance = false;
                toleranceLower = precursorMassLowerVal;
                toleranceUpper = precursorMassUpperVal;

                return true;
            }

            if (int.TryParse(precursorTrueUnits, out var precursorTrueUnitsVal) &&
                double.TryParse(precursorTrueTolerance, out var precursorTrueToleranceVal))
            {
                // A single, symmetric tolerance is defined
                ppmBased = MassToleranceUnitsArePPM(precursorTrueUnitsVal, "precursor_true_units");
                singleTolerance = true;
                toleranceLower = precursorTrueToleranceVal;
                toleranceUpper = precursorTrueToleranceVal;

                return true;
            }

            ppmBased = false;
            singleTolerance = true;
            toleranceLower = 0;
            toleranceUpper = 0;
            return false;
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MSFRAGGER_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Determine if the search engine parameters were loaded from a MSFragger parameter file, or from a FragPipe workflow file
        /// </summary>
        /// <param name="searchEngineParams">MSFragger or FragPipe parameters</param>
        /// <returns>Empty string if MSFragger parameters, "msfragger." if FragPipe workflow parameters</returns>
        public static string GetSearchEngineParameterPrefix(SearchEngineParameters searchEngineParams)
        {
            if (searchEngineParams.Parameters.TryGetValue("msfragger.precursor_mass_lower", out _) ||
                searchEngineParams.Parameters.TryGetValue("msfragger.precursor_mass_units", out _))
            {
                // The search engine parameter file is a FragPipe workflow file
                return "msfragger.";
            }

            // The search engine parameter file is a MSFragger parameter file
            return string.Empty;
        }

        /// <summary>
        /// Parses the specified MSFragger parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName">MSFragger parameter file name</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MSFRAGGER_SEARCH_ENGINE_NAME);

            return ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);
        }

        private bool MassToleranceUnitsArePPM(int massToleranceUnits, string parameterName)
        {
            switch (massToleranceUnits)
            {
                case 0:
                    // Dalton based
                    return false;

                case 1:
                    // ppm based
                    return true;

                default:
                    ReportWarning(string.Format(
                        "Unrecognized value for parameter {0} in the MSFragger parameter file: {1}",
                        parameterName, massToleranceUnits));

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

            var columns = line.Split('\t');

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.ResultID), mColumnHeaders, 0);
                psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.RankEValue), mColumnHeaders, 0);

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.Charge), mColumnHeaders, 0);

                var proteinName = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.Protein), mColumnHeaders);

                if (!string.IsNullOrWhiteSpace(proteinName))
                {
                    psm.AddProtein(proteinName.Trim());
                }

                var additionalProteins = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.AdditionalProteins), mColumnHeaders);

                if (!string.IsNullOrWhiteSpace(additionalProteins))
                {
                    foreach (var protein in additionalProteins.Split(';'))
                    {
                        if (string.IsNullOrWhiteSpace(protein))
                            continue;

                        psm.AddProtein(protein.Trim());
                    }
                }

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);

                if (Math.Abs(precursorMZ) > float.Epsilon)
                {
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);
                }

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.DelM), mColumnHeaders);
                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.DelM_PPM), mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                var elutionTime = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.ElutionTime), mColumnHeaders, string.Empty);

                if (float.TryParse(elutionTime, out var elutionTimeMinutes))
                {
                    // Update the elution time
                    psm.ElutionTimeMinutes = elutionTimeMinutes;
                }

                // Store the remaining data
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.Dataset));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.DatasetID));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.DelM_MSFragger));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.MH));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.Mass));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.Modifications));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.NTT));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.EValue));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.Hyperscore));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.Nextscore));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.PeptideProphetProbability));
                AddScore(psm, columns, elutionTime);
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.ElutionTimeAverage));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.MissedCleavages));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.NumberOfMatchedIons));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.TotalNumberOfIons));
                AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.QValue));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MSFragger data file: " + ex.Message);
                return false;
            }
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

                var success = ReadKeyValuePairSearchEngineParamFile(MSFRAGGER_SEARCH_ENGINE_NAME, paramFilePath, PeptideHitResultTypes.MSFragger, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                var parameterPrefix = GetSearchEngineParameterPrefix(searchEngineParams);

                string enzymeName;

                // Starting with MSFragger v3.4, enzyme name is specified by parameter search_enzyme_name_1
                // Previous versions used search_enzyme_name
                if (GetParameterValue(searchEngineParams, parameterPrefix, "search_enzyme_name_1", out var enzymeNameA))
                {
                    enzymeName = enzymeNameA;
                }
                else if (GetParameterValue(searchEngineParams, parameterPrefix, "search_enzyme_name", out var enzymeNameB))
                {
                    enzymeName = enzymeNameB;
                }
                else
                {
                    enzymeName = string.Empty;
                    ReportWarning("The MSFragger parameter file does not have parameter 'search_enzyme_name' or 'search_enzyme_name_1'");
                }

                // Determine the enzyme name
                if (!string.IsNullOrWhiteSpace(enzymeName))
                {
                    searchEngineParams.Enzyme = enzymeName;

                    // ReSharper disable StringLiteralTypo

                    switch (searchEngineParams.Enzyme)
                    {
                        case "argc":
                        case "aspn":
                        case "chymotrypsin":
                        case "cnbr":
                        case "elastase":
                        case "formicacid":
                        case "gluc":
                        case "gluc_bicarb":
                        case "lysc":
                        case "lysc-p":
                        case "lysn":
                        case "lysn_promisc":
                        case "nonspecific":
                        case "null":
                        case "stricttrypsin":
                        case "thermolysin":
                        case "trypsin":
                        case "trypsin/chymotrypsin":
                        case "trypsin/cnbr":
                        case "trypsin_gluc":
                        case "trypsin_k":
                        case "trypsin_r":
                            break;

                        default:
                            ReportWarning(string.Format("Unrecognized enzyme '{0}' in the MSFragger parameter file", searchEngineParams.Enzyme));
                            break;
                    }

                    // ReSharper restore StringLiteralTypo
                }

                // Determine the cleavage specificity

                if (!GetParameterValue(searchEngineParams, parameterPrefix, "num_enzyme_termini", out var numTermini))
                {
                    ReportWarning("'num_enzyme_termini' parameter not found in the MSFragger parameter file");
                }
                else
                {
                    switch (numTermini)
                    {
                        case "2":
                            // Fully-enzymatic (fully tryptic)
                            searchEngineParams.MinNumberTermini = 2;
                            break;

                        case "1":
                            // Semi-enzymatic
                            searchEngineParams.MinNumberTermini = 1;
                            break;

                        case "0":
                            // non-enzymatic (non-tryptic)
                            searchEngineParams.MinNumberTermini = 0;
                            break;

                        default:
                            ReportWarning(string.Format("Unrecognized value for num_enzyme_termini in the MSFragger parameter file: {0}", numTermini));
                            break;
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)

                var validTolerance = GetPrecursorSearchTolerances(
                    searchEngineParams,
                    parameterPrefix,
                    out var toleranceLower, out var toleranceUpper,
                    out var ppmBased, out _);

                if (!validTolerance)
                {
                    searchEngineParams.PrecursorMassToleranceDa = 0;
                    searchEngineParams.PrecursorMassTolerancePpm = 0;

                    return true;
                }

                UpdatePrecursorMassTolerance(searchEngineParams, toleranceLower, toleranceUpper, ppmBased);

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }

        /// <summary>
        /// Update PrecursorMassToleranceDa and PrecursorMassTolerancePpm in searchEngineParams
        /// </summary>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <param name="massToleranceLower">Tolerance to the left, e.g. -20</param>
        /// <param name="massToleranceUpper">Tolerance to the right, e.g. 20</param>
        /// <param name="ppmBased">True if ppm-based tolerances</param>
        private void UpdatePrecursorMassTolerance(
            SearchEngineParameters searchEngineParams,
            double massToleranceLower,
            double massToleranceUpper,
            bool ppmBased)
        {
            double toleranceToStore;

            if (Math.Abs(Math.Abs(massToleranceLower) - Math.Abs(massToleranceUpper)) < float.Epsilon)
            {
                toleranceToStore = Math.Abs(massToleranceUpper);
            }
            else
            {
                // Tolerances are not symmetric; store the average value
                toleranceToStore = (Math.Abs(massToleranceLower) + Math.Abs(massToleranceUpper)) / 2.0;
            }

            if (ppmBased)
            {
                searchEngineParams.PrecursorMassToleranceDa = PeptideMassCalculator.PPMToMass(toleranceToStore, 2000);
                searchEngineParams.PrecursorMassTolerancePpm = toleranceToStore;
            }
            else
            {
                searchEngineParams.PrecursorMassToleranceDa = toleranceToStore;
                searchEngineParams.PrecursorMassTolerancePpm = PeptideMassCalculator.MassToPPM(toleranceToStore, 1000);
            }
        }
    }
}
