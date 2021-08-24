//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from SEQUEST _syn.txt and _fht.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for SEQUEST
    /// </summary>
    public class SequestSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: Chymotrypsin, Cn, Da, Daltons, GluC, LysC, milli, mmu, PassFilt, Prot, tryptic, Xc

        /// <summary>
        /// SEQUEST synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_syn.txt";

        /// <summary>
        /// SEQUEST first hits file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_FHT = "_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string SEQ_SEARCH_ENGINE_NAME = "SEQUEST";

        /// <summary>
        /// Mapping from enum to synopsis file column name for SEQUEST
        /// </summary>
        private static readonly Dictionary<SequestSynopsisFileColumns, string> mSynopsisFileColumn = new();

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
        public SequestSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public SequestSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.Sequest, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public SequestSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.Sequest, startupOptions)
        {
        }

        /// <summary>
        /// Define column header names for SEQUEST synopsis and first hits files
        /// </summary>
        protected override void DefineColumnHeaders()
        {
            base.DefineColumnHeaders();

            // These are legacy column names that aren't always present, so change their mapping to -1
            mColumnHeaders["PassFilt"] = -1;
            mColumnHeaders["MScore"] = -1;
        }

        /// <summary>
        /// Determines the precursor mass tolerance
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <returns>Precursor tolerance, in Da</returns>
        private double DeterminePrecursorMassTolerance(SearchEngineParameters searchEngineParams, out double tolerancePPM)
        {
            double toleranceDa = 0;

            tolerancePPM = 0;

            if (searchEngineParams.Parameters.TryGetValue("peptide_mass_tolerance", out var peptideMassTolerance))
            {
                if (double.TryParse(peptideMassTolerance, out var value))
                {
                    // Determine the mass units
                    // 0 means Da, 1 means mmu, 2 means ppm
                    var units = 0;

                    if (searchEngineParams.Parameters.TryGetValue("peptide_mass_units", out var peptideMassUnits))
                    {
                        if (!string.IsNullOrEmpty(peptideMassUnits))
                        {
                            int.TryParse(peptideMassUnits, out units);
                        }
                    }

                    if (units == 2)
                    {
                        // Tolerance is in ppm; convert to Da at 2000 m/z
                        tolerancePPM = value;

                        toleranceDa = PeptideMassCalculator.PPMToMass(value, 2000);
                    }
                    else if (units == 1)
                    {
                        // Tolerance is in mmu (milli mass units)
                        toleranceDa = value / 1000.0;

                        // Convert from Dalton to PPM (assuming a mass of 2000 m/z)
                        tolerancePPM = PeptideMassCalculator.MassToPPM(toleranceDa, 2000);
                    }
                    else
                    {
                        // Tolerance is in Daltons
                        toleranceDa = value;

                        // Convert from Dalton to PPM (assuming a mass of 2000 m/z)
                        tolerancePPM = PeptideMassCalculator.MassToPPM(toleranceDa, 2000);
                    }
                }
            }

            return toleranceDa;
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
        public static SortedDictionary<string, SequestSynopsisFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new(StringComparer.OrdinalIgnoreCase)
            {
                { "HitNum", SequestSynopsisFileColumns.RowIndex },
                { "ScanNum", SequestSynopsisFileColumns.Scan },
                { "ScanCount", SequestSynopsisFileColumns.NumScans },
                { "ChargeState", SequestSynopsisFileColumns.Charge },
                { "MH", SequestSynopsisFileColumns.PeptideMH },
                { "XCorr", SequestSynopsisFileColumns.XCorr },
                { "DelCn", SequestSynopsisFileColumns.DeltaCn },
                { "Sp", SequestSynopsisFileColumns.Sp },
                { "Reference", SequestSynopsisFileColumns.ProteinName },
                { "MultiProtein", SequestSynopsisFileColumns.MultipleProteinCount },     // Multiple protein count: 0 if the peptide is in 1 protein, 1 if the peptide is in 2 proteins, etc.
                { "Peptide", SequestSynopsisFileColumns.PeptideSequence },
                { "DelCn2", SequestSynopsisFileColumns.DeltaCn2 },
                { "RankSp", SequestSynopsisFileColumns.RankSP },
                { "RankXc", SequestSynopsisFileColumns.RankXC },
                { "DelM", SequestSynopsisFileColumns.DelM },
                { "XcRatio", SequestSynopsisFileColumns.XcRatio },
                { "PassFilt", SequestSynopsisFileColumns.PassFilt },                    // Legacy/unused
                { "MScore", SequestSynopsisFileColumns.MScore },                        // Legacy/unused
                { "NumTrypticEnds", SequestSynopsisFileColumns.NTT },
                { "Ions_Observed", SequestSynopsisFileColumns.IonsObserved },
                { "Ions_Expected", SequestSynopsisFileColumns.IonsExpected },
                { "DelM_PPM", SequestSynopsisFileColumns.DelMPPM },
                { "Cleavage_State", SequestSynopsisFileColumns.Cleavage_State },         // Computed by this program and appended to the input file or saved in a new file
                { "Terminus_State", SequestSynopsisFileColumns.Terminus_State },         // Computed by this program
                { "Mod_Count", SequestSynopsisFileColumns.Mod_Count },                   // Computed by this program
                { "Mod_Description", SequestSynopsisFileColumns.Mod_Description },       // Computed by this program
                { "Monoisotopic_Mass", SequestSynopsisFileColumns.Monoisotopic_Mass }    // Computed by this program
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum SequestSynopsisFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<SequestSynopsisFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(SequestSynopsisFileColumns column)
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
            return datasetName + FILENAME_SUFFIX_FHT;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_syn_ProteinMods.txt";
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
            return datasetName + "_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return SEQ_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified SEQUEST parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(SEQ_SEARCH_ENGINE_NAME, mModInfo)
            {
                Enzyme = "trypsin"
            };

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            var enzymeSpecificityMatcher = new Regex(@"^\S+\s(\d)\s\d\s.+", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            try
            {
                var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

                if (!File.Exists(paramFilePath))
                {
                    ReportError("SEQUEST param file not found: " + paramFilePath);
                    return false;
                }

                using var reader = new StreamReader(new FileStream(paramFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(lineIn))
                        continue;

                    var dataLine = lineIn.TrimStart();

                    if (dataLine.StartsWith(";") || dataLine.StartsWith("[") || !dataLine.Contains("="))
                        continue;

                    // Split the line on the equals sign
                    var kvSetting = ParseKeyValueSetting(dataLine, '=');

                    if (string.IsNullOrEmpty(kvSetting.Key))
                        continue;

                    // Trim off any text that occurs after a semicolon in kvSetting.Value
                    var settingValue = kvSetting.Value;
                    var charIndex = settingValue.IndexOf(';');
                    if (charIndex > 0)
                    {
                        settingValue = settingValue.Substring(charIndex).Trim();
                    }

                    searchEngineParams.AddUpdateParameter(kvSetting.Key, settingValue);

                    int value;
                    switch (kvSetting.Key.ToLower())
                    {
                        case "first_database_name":
                        case "database_name":
                            string fastaFilePath;
                            try
                            {
                                fastaFilePath = Path.Combine(@"C:\Database", Path.GetFileName(settingValue));
                            }
                            catch (Exception)
                            {
                                fastaFilePath = settingValue;
                            }
                            searchEngineParams.FastaFilePath = fastaFilePath;

                            break;
                        case "mass_type_parent":
                            if (settingValue == "0")
                            {
                                // Average mass
                                searchEngineParams.PrecursorMassType = SearchEngineParameters.MASS_TYPE_AVERAGE;
                            }
                            else
                            {
                                // Monoisotopic mass
                                searchEngineParams.PrecursorMassType = SearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
                            }

                            break;
                        case "mass_type_fragment":
                            if (settingValue == "0")
                            {
                                // Average mass
                                searchEngineParams.PrecursorMassType = SearchEngineParameters.MASS_TYPE_AVERAGE;
                            }
                            else
                            {
                                // Monoisotopic mass
                                searchEngineParams.PrecursorMassType = SearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
                            }

                            break;
                        case "max_num_internal_cleavage_sites":
                            if (int.TryParse(settingValue, out value))
                            {
                                searchEngineParams.MaxNumberInternalCleavages = value;
                            }

                            break;
                        case "enzyme_info":
                            // Used in new-style SEQUEST parameter files

                            // Examples include:
                            // Fully-tryptic:     Trypsin(KR) 1 1 KR -
                            // Partially-tryptic: Trypsin(KR) 2 1 KR -
                            // No-enzyme:         No_Enzyme(-) 0 0 - -
                            //
                            searchEngineParams.Enzyme = "trypsin";

                            if (settingValue.StartsWith("no_enzyme", StringComparison.OrdinalIgnoreCase))
                            {
                                searchEngineParams.MinNumberTermini = 0;
                            }
                            else
                            {
                                // Parse out the cleavage specificity number
                                // This is the first number after the closing parenthesis in the above examples
                                var match = enzymeSpecificityMatcher.Match(settingValue);
                                if (match.Success)
                                {
                                    if (int.TryParse(match.Groups[1].Value, out value))
                                    {
                                        searchEngineParams.MinNumberTermini = value;
                                    }
                                }
                            }

                            break;
                        case "enzyme_number":
                            // Used in old-style SEQUEST parameter files
                            if (int.TryParse(settingValue, out value))
                            {
                                if (value == 0)
                                {
                                    // No-enzyme
                                    searchEngineParams.Enzyme = "trypsin";
                                    searchEngineParams.MinNumberTermini = 0;
                                }
                                else
                                {
                                    searchEngineParams.Enzyme = value switch
                                    {
                                        1 => "trypsin",
                                        2 => "trypsin_modified",
                                        3 => "Chymotrypsin",
                                        4 => "Chymotrypsin_modified",
                                        5 => "Clostripain",
                                        6 => "Cyanogen_Bromide",
                                        7 => "IodosoBenzoate",
                                        8 => "Proline_Endopept",
                                        9 => "Staph_Protease",
                                        10 => "Trypsin_K",
                                        11 => "Trypsin_R",
                                        12 => "GluC",
                                        13 => "LysC",
                                        14 => "AspN",
                                        15 => "Elastase",
                                        16 => "Elastase/Tryp/Chymo",
                                        _ => "Unknown"
                                    };
                                    searchEngineParams.MinNumberTermini = 2;
                                }
                            }
                            break;
                    }
                }

                // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM);
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

            var columns = line.Split('\t');

            var success = false;

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.RowIndex), mColumnHeaders, 0);
                    psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.RankXC), mColumnHeaders, 1);

                    var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.PeptideSequence), mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.Charge), mColumnHeaders, 0);

                    var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.ProteinName), mColumnHeaders);
                    psm.AddProtein(protein);

                    // Note that the MH value listed in SEQUEST files is not the precursor MH but is instead the theoretical (computed) MH of the peptide
                    // We'll update this value below using massErrorDa
                    // We'll further update this value using the ScanStatsEx data
                    var precursorMH = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.PeptideMH), mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMH, 1, 0);

                    psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.DelM), mColumnHeaders);
                    if (double.TryParse(psm.MassErrorDa, out var massErrorDa))
                    {
                        // Adjust the precursor mass
                        psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMH - massErrorDa, 1, 0);
                    }

                    psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(SequestSynopsisFileColumns.DelMPPM), mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.XCorr));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.DeltaCn));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.Sp));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.DeltaCn2));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.RankSP));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.RankXC));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.XcRatio));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.IonsObserved));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.IonsExpected));
                    AddScore(psm, columns, GetColumnNameByID(SequestSynopsisFileColumns.NTT));
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the SEQUEST data file: " + ex.Message);
            }

            return success;
        }
    }
}
