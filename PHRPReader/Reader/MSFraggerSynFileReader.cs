﻿using PHRPReader.Data;
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

        // Ignore Spelling: Chymotrypsin, Daltons, Hyperscore, MSFragger, NextScore, PepToProtMap, tryptic
        // Ignore Spelling: argc, aspn, cnbr, elastase, formicacid, gluc, lysc, lysn, stricttrypsin, thermolysin

        // ReSharper restore CommentTypo

        /// <summary>
        /// MSFragger synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_msf_syn.txt";

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
        public const string FILENAME_SUFFIX_FHT = "_msf_fht.txt";

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
                { "RetentionTime", MSFraggerSynFileColumns.RetentionTime },
                { "MissedCleavages", MSFraggerSynFileColumns.MissedCleavages },
                { "QValue", MSFraggerSynFileColumns.QValue }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSFraggerSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
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
        /// <param name="column"></param>
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

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // MSFragger does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_msf_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msf_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_msf_syn_ProteinMods.txt";
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
            return datasetName + "_msf_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_msf_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msf_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MSFRAGGER_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MSFragger parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MSFRAGGER_SEARCH_ENGINE_NAME);

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
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MSFraggerSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
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

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
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
                    AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.RetentionTime));
                    AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.MissedCleavages));
                    AddScore(psm, columns, GetColumnNameByID(MSFraggerSynFileColumns.QValue));
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MSFragger data file: " + ex.Message);
            }

            return success;
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

                // Determine the enzyme name
                if (!searchEngineParams.Parameters.TryGetValue("search_enzyme_name", out var enzymeName))
                {
                    ReportWarning("'search_enzyme_name' parameter not found in the MSFragger parameter file");
                }
                else
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

                if (!searchEngineParams.Parameters.TryGetValue("num_enzyme_termini", out var numTermini))
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

                // First look for the legacy symmetric tolerance parameter
                searchEngineParams.Parameters.TryGetValue("precursor_true_tolerance", out var precursorTrueTolerance);
                searchEngineParams.Parameters.TryGetValue("precursor_true_units", out var precursorTrueUnits);

                // Next look for the newer tolerance parameters
                searchEngineParams.Parameters.TryGetValue("precursor_mass_lower", out var precursorMassLower);
                searchEngineParams.Parameters.TryGetValue("precursor_mass_upper", out var precursorMassUpper);
                searchEngineParams.Parameters.TryGetValue("precursor_mass_units", out var precursorMassUnits);

                if (int.TryParse(precursorMassUnits, out var precursorMassUnitsVal) &&
                    double.TryParse(precursorMassLower, out var precursorMassLowerVal) &&
                    double.TryParse(precursorMassUpper, out var precursorMassUpperVal))
                {
                    UpdatePrecursorMassTolerance(
                        searchEngineParams, precursorMassUnitsVal,
                        precursorMassLowerVal, precursorMassUpperVal);
                }
                else if (int.TryParse(precursorTrueUnits, out var precursorTrueUnitsVal) &&
                         double.TryParse(precursorTrueTolerance, out var precursorTrueToleranceVal))
                {
                    UpdatePrecursorMassTolerance(
                        searchEngineParams, precursorTrueUnitsVal,
                        precursorTrueToleranceVal, precursorTrueToleranceVal);
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

        /// <summary>
        /// Update PrecursorMassToleranceDa and PrecursorMassTolerancePpm in searchEngineParams
        /// </summary>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <param name="massToleranceUnits">Precursor mass units: 0 means Daltons, 1 means ppm</param>
        /// <param name="massToleranceLower">Tolerance to the left, e.g. -20</param>
        /// <param name="massToleranceUpper">Tolerance to the right, e.g. 20</param>
        private void UpdatePrecursorMassTolerance(
            SearchEngineParameters searchEngineParams,
            int massToleranceUnits,
            double massToleranceLower,
            double massToleranceUpper)
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

            switch (massToleranceUnits)
            {
                case 0:
                    // Dalton based
                    searchEngineParams.PrecursorMassToleranceDa = toleranceToStore;
                    searchEngineParams.PrecursorMassTolerancePpm = PeptideMassCalculator.MassToPPM(toleranceToStore, 1000);
                    break;

                case 1:
                    // ppm based
                    searchEngineParams.PrecursorMassToleranceDa = PeptideMassCalculator.PPMToMass(toleranceToStore, 2000);
                    searchEngineParams.PrecursorMassTolerancePpm = toleranceToStore;
                    break;

                default:
                    ReportWarning(string.Format("Unrecognized value for precursor_mass_units in the MSFragger parameter file: {0}", massToleranceUnits));
                    break;
            }
        }
    }
}