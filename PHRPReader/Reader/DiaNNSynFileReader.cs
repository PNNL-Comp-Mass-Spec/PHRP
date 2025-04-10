﻿using PHRPReader.Data;
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

        // Ignore Spelling: aspn, Averagine, chymotrypsin, Dia, diann, FHT, Glu, gluc, Lys, lysc
        // Ignore Spelling: proline, proteotypic, psm, tryptic

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
            headerNames.AddRange(GetColumnHeaderNamesAndIDs(false).Keys);
            return headerNames;
        }

        /// <summary>
        /// Header names and enums for the PHRP synopsis file for this tool
        /// </summary>
        /// <param name="excludeColumnsNotInParquetFile">When true, do not include column names that are not in .parquet files (used by DIA-NN 2.x)</param>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, DiaNNSynFileColumns> GetColumnHeaderNamesAndIDs(bool excludeColumnsNotInParquetFile)
        {
            var headerColumns = new SortedDictionary<string, DiaNNSynFileColumns>(StringComparer.OrdinalIgnoreCase)
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
                { "ProteinGroup", DiaNNSynFileColumns.ProteinGroup },
                { "ProteinIDs", DiaNNSynFileColumns.ProteinIDs },
                { "ProteinNames", DiaNNSynFileColumns.ProteinNames },
                { "GeneNames", DiaNNSynFileColumns.GeneNames },
                { "NTT", DiaNNSynFileColumns.NTT }
            };

            if (!excludeColumnsNotInParquetFile)
            {
                headerColumns.Add("ProteinGroupQuantity", DiaNNSynFileColumns.ProteinGroupQuantity);
                headerColumns.Add("ProteinGroupNormalized", DiaNNSynFileColumns.ProteinGroupNormalized);
            }

            headerColumns.Add("ProteinGroupMaxLFQ", DiaNNSynFileColumns.ProteinGroupMaxLFQ);

            if (!excludeColumnsNotInParquetFile)
            {
                headerColumns.Add("GenesQuantity", DiaNNSynFileColumns.GenesQuantity);
                headerColumns.Add("GenesNormalized", DiaNNSynFileColumns.GenesNormalized);
            }

            headerColumns.Add("GenesMaxLFQ", DiaNNSynFileColumns.GenesMaxLFQ);
            headerColumns.Add("GenesMaxLFQUnique", DiaNNSynFileColumns.GenesMaxLFQUnique);
            headerColumns.Add("QValue", DiaNNSynFileColumns.QValue);
            headerColumns.Add("PEP", DiaNNSynFileColumns.PEP);
            headerColumns.Add("GlobalQValue", DiaNNSynFileColumns.GlobalQValue);
            headerColumns.Add("ProteinQValue", DiaNNSynFileColumns.ProteinQValue);
            headerColumns.Add("ProteinGroupQValue", DiaNNSynFileColumns.ProteinGroupQValue);
            headerColumns.Add("GlobalProteinGroupQValue", DiaNNSynFileColumns.GlobalProteinGroupQValue);
            headerColumns.Add("GeneGroupQValue", DiaNNSynFileColumns.GeneGroupQValue);
            headerColumns.Add("TranslatedQValue", DiaNNSynFileColumns.TranslatedQValue);
            headerColumns.Add("Proteotypic", DiaNNSynFileColumns.Proteotypic);
            headerColumns.Add("PrecursorQuantity", DiaNNSynFileColumns.PrecursorQuantity);
            headerColumns.Add("PrecursorNormalized", DiaNNSynFileColumns.PrecursorNormalized);

            if (!excludeColumnsNotInParquetFile)
            {
                headerColumns.Add("PrecursorTranslated", DiaNNSynFileColumns.PrecursorTranslated);
                headerColumns.Add("TranslatedQuality", DiaNNSynFileColumns.TranslatedQuality);
                headerColumns.Add("MS1Translated", DiaNNSynFileColumns.MS1Translated);
            }

            headerColumns.Add("QuantityQuality", DiaNNSynFileColumns.QuantityQuality);
            headerColumns.Add("ElutionTime", DiaNNSynFileColumns.ElutionTime);
            headerColumns.Add("ElutionTimeStart", DiaNNSynFileColumns.ElutionTimeStart);
            headerColumns.Add("ElutionTimeStop", DiaNNSynFileColumns.ElutionTimeStop);
            headerColumns.Add("IndexedRT", DiaNNSynFileColumns.IndexedRT);
            headerColumns.Add("IndexedIonMobility", DiaNNSynFileColumns.IndexedIonMobility);
            headerColumns.Add("PredictedRT", DiaNNSynFileColumns.PredictedRT);
            headerColumns.Add("PredictedIndexedRT", DiaNNSynFileColumns.PredictedIndexedRT);
            headerColumns.Add("MS1ProfileCorrelation", DiaNNSynFileColumns.MS1ProfileCorrelation);
            headerColumns.Add("MS1Area", DiaNNSynFileColumns.MS1Area);
            headerColumns.Add("Evidence", DiaNNSynFileColumns.Evidence);

            if (!excludeColumnsNotInParquetFile)
            {
                headerColumns.Add("SpectrumSimilarity", DiaNNSynFileColumns.SpectrumSimilarity);
                headerColumns.Add("Averagine", DiaNNSynFileColumns.Averagine);
            }

            headerColumns.Add("MassEvidence", DiaNNSynFileColumns.MassEvidence);

            // ReSharper disable once InvertIf
            if (!excludeColumnsNotInParquetFile)
            {
                headerColumns.Add("CScore", DiaNNSynFileColumns.CScore);
                headerColumns.Add("DecoyEvidence", DiaNNSynFileColumns.DecoyEvidence);
                headerColumns.Add("DecoyCScore", DiaNNSynFileColumns.DecoyCScore);
            }

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum DiaNNSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames">List of header names</param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<DiaNNSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs(false);
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column">Column enum</param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(DiaNNSynFileColumns column)
        {
            if (mSynopsisFileColumn.Count > 0)
            {
                return mSynopsisFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs(false))
            {
                mSynopsisFileColumn.Add(item.Value, item.Key);
            }

            return mSynopsisFileColumn[column];
        }

        /// <summary>
        /// Default first hits file for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // DIA-NN does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_diann_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_diann_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_diann_syn_ProteinMods.txt";
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
            return datasetName + "_diann_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_diann_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_diann_syn_SeqToProteinMap.txt";
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
        /// <param name="searchEngineParamFileName">DIA-NN parameter file name</param>
        /// <param name="searchEngineParams">Search engine parameters</param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(DiaNN_SEARCH_ENGINE_NAME);

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

                // Since this is DIA data, always use 1 for score rank
                psm.ScoreRank = 1;

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    if (peptide.IndexOf('.') < 0)
                    {
                        // Assume the peptide is tryptic
                        psm.SetPeptide(peptide, mCleavageStateCalculator, true);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.Charge), mColumnHeaders, 0);

                var proteinIDs = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinIDs), mColumnHeaders);

                foreach (var proteinID in ParseDelimitedList(proteinIDs, ';'))
                {
                    psm.AddProtein(proteinID);
                }

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);

                if (Math.Abs(precursorMZ) > float.Epsilon)
                {
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);
                }

                // DIA-NN results do not have mass error information
                // psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.DelM), mColumnHeaders);
                // psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.DelM_PPM), mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                var elutionTime = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(DiaNNSynFileColumns.ElutionTime), mColumnHeaders, string.Empty);

                if (float.TryParse(elutionTime, out var elutionTimeMinutes))
                {
                    // Update the elution time
                    psm.ElutionTimeMinutes = elutionTimeMinutes;
                }

                // Store the remaining data
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Dataset));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.DatasetID));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.IonMobility));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MH));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Mass));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Peptide));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Modifications));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinGroup));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinNames));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GeneNames));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.NTT));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinGroupQuantity));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinGroupNormalized));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinGroupMaxLFQ));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GenesQuantity));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GenesNormalized));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GenesMaxLFQ));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GenesMaxLFQUnique));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.QValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PEP));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GlobalQValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinQValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ProteinGroupQValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GlobalProteinGroupQValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.GeneGroupQValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.TranslatedQValue));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Proteotypic));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PrecursorQuantity));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PrecursorNormalized));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PrecursorTranslated));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.TranslatedQuality));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MS1Translated));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.QuantityQuality));
                AddScore(psm, columns, elutionTime);
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ElutionTimeStart));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.ElutionTimeStop));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.IndexedRT));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.IndexedIonMobility));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PredictedRT));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.PredictedIndexedRT));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MS1ProfileCorrelation));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MS1Area));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Evidence));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.SpectrumSimilarity));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.Averagine));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.MassEvidence));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.CScore));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.DecoyEvidence));
                AddScore(psm, columns, GetColumnNameByID(DiaNNSynFileColumns.DecoyCScore));

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
