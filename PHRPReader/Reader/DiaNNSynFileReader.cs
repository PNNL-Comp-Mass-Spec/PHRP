using PHRPReader.Data;
using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for DIA-NN
    /// </summary>
    public class DiaNNSynFileReader : SynFileReaderBaseClass
    {
        // ReSharper disable CommentTypo

        // Ignore Spelling: Chymotrypsin, Daltons, Hyperscore, DiaNN, NextScore, PepToProtMap, tryptic
        // Ignore Spelling: argc, aspn, cnbr, elastase, formicacid, gluc, lysc, lysn, stricttrypsin, thermolysin

        // ReSharper restore CommentTypo

        /// <summary>
        /// DIA-NN synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_diann_syn.txt";

        /// <summary>
        /// DIA-NN first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for DIA-NN results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_diann_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        // ReSharper disable once IdentifierTypo
        private const string DiaNN_SEARCH_ENGINE_NAME = "DIA-NN";

        /// <summary>
        /// Mapping from enum to synopsis file column name for DIA-NN
        /// </summary>
        private static readonly Dictionary<DiaNNSynFileColumns, string> mSynopsisFileColumn = new();

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
        public DiaNNSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public DiaNNSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.DiaNN, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public DiaNNSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.DiaNN, startupOptions)
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
        public static SortedDictionary<string, DiaNNSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, DiaNNSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", DiaNNSynFileColumns.ResultID },
                { "Dataset", DiaNNSynFileColumns.Dataset },
                { "DatasetID", DiaNNSynFileColumns.DatasetID },
                { "Scan", DiaNNSynFileColumns.Scan },
                { "IonMobility", DiaNNSynFileColumns.IonMobility },
                { "Charge", DiaNNSynFileColumns.Charge },
                { "PrecursorMZ", DiaNNSynFileColumns.PrecursorMZ },
                { "MH", DiaNNSynFileColumns.MH },
                { "Mass", DiaNNSynFileColumns.Mass },
                { "Peptide", DiaNNSynFileColumns.Peptide },
                { "Modifications", DiaNNSynFileColumns.Modifications },
                { "Peptide", DiaNNSynFileColumns.Peptide },
                { "Modifications", DiaNNSynFileColumns.Modifications },
                { "ProteinGroup", DiaNNSynFileColumns.ProteinGroup },
                { "ProteinIDs", DiaNNSynFileColumns.ProteinIDs },
                { "ProteinNames", DiaNNSynFileColumns.ProteinNames },
                { "Genes", DiaNNSynFileColumns.GeneNames },
                { "NTT", DiaNNSynFileColumns.NTT },
                { "ProteinGroupQuantity", DiaNNSynFileColumns.ProteinGroupQuantity },
                { "ProteinGroupNormalized", DiaNNSynFileColumns.ProteinGroupNormalized },
                { "ProteinGroupMaxLFQ", DiaNNSynFileColumns.ProteinGroupMaxLFQ },
                { "GenesQuantity", DiaNNSynFileColumns.GenesQuantity },
                { "GenesNormalized", DiaNNSynFileColumns.GenesNormalized },
                { "GenesMaxLFQ", DiaNNSynFileColumns.GenesMaxLFQ },
                { "GenesMaxLFQUnique", DiaNNSynFileColumns.GenesMaxLFQUnique },
                { "QValue", DiaNNSynFileColumns.QValue },
                { "PEP", DiaNNSynFileColumns.PEP },
                { "GlobalQValue", DiaNNSynFileColumns.GlobalQValue },
                { "ProteinQValue", DiaNNSynFileColumns.ProteinQValue },
                { "ProteinGroupQValue", DiaNNSynFileColumns.ProteinGroupQValue },
                { "GlobalProteinGroupQValue", DiaNNSynFileColumns.GlobalProteinGroupQValue },
                { "GeneGroupQValue", DiaNNSynFileColumns.GeneGroupQValue },
                { "TranslatedQValue", DiaNNSynFileColumns.TranslatedQValue },
                { "PrecursorQuantity", DiaNNSynFileColumns.PrecursorQuantity },
                { "PrecursorNormalized", DiaNNSynFileColumns.PrecursorNormalized },
                { "PrecursorTranslated", DiaNNSynFileColumns.PrecursorTranslated },
                { "TranslatedQuality", DiaNNSynFileColumns.TranslatedQuality },
                { "MS1Translated", DiaNNSynFileColumns.MS1Translated },
                { "QuantityQuality", DiaNNSynFileColumns.QuantityQuality },
                { "ElutionTime", DiaNNSynFileColumns.ElutionTime },
                { "ElutionTimeStart", DiaNNSynFileColumns.ElutionTimeStart },
                { "ElutionTimeStop", DiaNNSynFileColumns.ElutionTimeStop },
                { "IndexedRT", DiaNNSynFileColumns.IndexedRT },
                { "PredictedRT", DiaNNSynFileColumns.PredictedRT },
                { "PredictedIndexedRT", DiaNNSynFileColumns.PredictedIndexedRT },
                { "MS1ProfileCorrelation", DiaNNSynFileColumns.MS1ProfileCorrelation },
                { "MS1Area", DiaNNSynFileColumns.MS1Area },
                { "Evidence", DiaNNSynFileColumns.Evidence },
                { "SpectrumSimilarity", DiaNNSynFileColumns.SpectrumSimilarity },
                { "Averagine", DiaNNSynFileColumns.Averagine },
                { "MassEvidence", DiaNNSynFileColumns.MassEvidence },
                { "CScore", DiaNNSynFileColumns.CScore },
                { "DecoyEvidence", DiaNNSynFileColumns.DecoyEvidence },
                { "DecoyCScore", DiaNNSynFileColumns.DecoyCScore },
                { "Scan", DiaNNSynFileColumns.Scan },
                { "IonMobility", DiaNNSynFileColumns.IonMobility },
                { "IndexedIonMobility", DiaNNSynFileColumns.IndexedIonMobility }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum DiaNNSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<DiaNNSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(DiaNNSynFileColumns column)
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
            // DIA-NN does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_diann_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_diann_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_diann_syn_ProteinMods.txt";
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
            return datasetName + "_diann_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_diann_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_diann_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Examine the DiaNN parameters to determine the precursor mass tolerance(s)
        /// </summary>
        /// <param name="searchEngineParams">DiaNN search engine parameters loaded from a Key=Value parameter file</param>
        /// <param name="toleranceLower">Output: Tolerance to the left, e.g. -20</param>
        /// <param name="toleranceUpper">Output: Tolerance to the right, e.g. 20</param>
        /// <param name="ppmBased">Output: True if ppm-based tolerances</param>
        /// <param name="singleTolerance">Output: true if a single tolerance is defined, false if two tolerances are defined</param>
        /// <returns>True if the tolerance parameters were found, false if not found or an error</returns>
        public bool GetPrecursorSearchTolerances(
            SearchEngineParameters searchEngineParams,
            out double toleranceLower,
            out double toleranceUpper,
            out bool ppmBased,
            out bool singleTolerance)
        {
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
            return DiaNN_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified DIA-NN parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(DiaNN_SEARCH_ENGINE_NAME);

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
                        "Unrecognized value for parameter {0} in the DiaNN parameter file: {1}",
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
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.ResultID), mColumnHeaders, 0);
                psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.RankEValue), mColumnHeaders, 0);

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.Charge), mColumnHeaders, 0);

                var proteinName = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.Protein), mColumnHeaders);

                if (!string.IsNullOrWhiteSpace(proteinName))
                {
                    psm.AddProtein(proteinName.Trim());
                }

                var additionalProteins = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.AdditionalProteins), mColumnHeaders);
                if (!string.IsNullOrWhiteSpace(additionalProteins))
                {
                    foreach (var protein in additionalProteins.Split(';'))
                    {
                        if (string.IsNullOrWhiteSpace(protein))
                            continue;

                        psm.AddProtein(protein.Trim());
                    }
                }

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);

                if (Math.Abs(precursorMZ) > float.Epsilon)
                {
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);
                }

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.DelM), mColumnHeaders);
                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.DelM_PPM), mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                // Store the remaining data
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Dataset));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.DatasetID));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.DelM_DiaNN));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MH));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Mass));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Modifications));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.NTT));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.EValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Hyperscore));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Nextscore));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PeptideProphetProbability));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ElutionTime));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ElutionTimeAverage));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MissedCleavages));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.NumberOfMatchedIons));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.TotalNumberOfIons));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.QValue));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the DIA-NN data file: " + ex.Message);
                return false;
            }
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

                var success = ReadKeyValuePairSearchEngineParamFile(DiaNN_SEARCH_ENGINE_NAME, paramFilePath, PeptideHitResultTypes.DiaNN, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                if (!searchEngineParams.Parameters.TryGetValue("CleavageSpecificity", out var cleavageSpecificity))
                {
                    cleavageSpecificity = string.Empty;
                    ReportWarning("The DIA-NN parameter file does not have parameter 'CleavageSpecificity'");
                }

                // Determine the enzyme name
                if (!string.IsNullOrWhiteSpace(cleavageSpecificity))
                {
                    // ReSharper disable StringLiteralTypo

                    switch (cleavageSpecificity)
                    {
                        case "K*,R*":           // Trypsin (ignore the proline rule)
                        case "K*,R*,!*P":       // Strict trypsin
                            searchEngineParams.Enzyme = "trypsin";
                            break;

                        case "K*":              // Lys/C
                            searchEngineParams.Enzyme = "lysc";
                            break;

                        case "F*,W*,Y*,L*":     // Chymotrypsin
                            searchEngineParams.Enzyme = "chymotrypsin";
                            break;

                        case "D*":              // Asp-N
                            searchEngineParams.Enzyme = "aspn";
                            break;

                        case "E*,D*":           // Glu-C
                            searchEngineParams.Enzyme = "gluc";
                            break;

                        default:
                            ReportWarning(string.Format("Unrecognized cleavage specificity '{0}' in the DIA-NN parameter file", cleavageSpecificity));
                            searchEngineParams.Enzyme = cleavageSpecificity;
                            break;
                    }

                    // ReSharper restore StringLiteralTypo
                }

                // Assume fully-tryptic
                searchEngineParams.MinNumberTermini = 2;

                if (searchEngineParams.Parameters.TryGetValue("MissedCleavages", out var missedCleavages) && int.TryParse(missedCleavages, out var cleavageCount))
                {
                    searchEngineParams.MaxNumberInternalCleavages = cleavageCount;
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
