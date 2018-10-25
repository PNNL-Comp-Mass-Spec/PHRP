//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
// This class parses data lines from Sequest _syn.txt and _fht.txt ifles
//
//*********************************************************************************************************
using System;
using System.IO;
using System.Text.RegularExpressions;

namespace PHRPReader
{
    /// <summary>
    /// PHRP parser for SEQUEST
    /// </summary>
    public class clsPHRPParserSequest : clsPHRPParser
    {
        #region "Constants"

#pragma warning disable 1591
        public const string DATA_COLUMN_HitNum = "HitNum";
        public const string DATA_COLUMN_ScanNum = "ScanNum";
        public const string DATA_COLUMN_ScanCount = "ScanCount";
        public const string DATA_COLUMN_ChargeState = "ChargeState";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_XCorr = "XCorr";
        public const string DATA_COLUMN_DelCn = "DelCn";
        public const string DATA_COLUMN_Sp = "Sp";
        public const string DATA_COLUMN_Reference = "Reference";
        public const string DATA_COLUMN_MultiProtein = "MultiProtein";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_DelCn2 = "DelCn2";
        public const string DATA_COLUMN_RankSp = "RankSp";
        public const string DATA_COLUMN_RankXc = "RankXc";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_XcRatio = "XcRatio";
        public const string DATA_COLUMN_PassFilt = "PassFilt";
        public const string DATA_COLUMN_MScore = "MScore";
        public const string DATA_COLUMN_Ions_Observed = "Ions_Observed";
        public const string DATA_COLUMN_Ions_Expected = "Ions_Expected";
        public const string DATA_COLUMN_NumTrypticEnds = "NumTrypticEnds";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";

        public const string FILENAME_SUFFIX_SYN = "_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_fht.txt";

        private const string SEQ_SEARCH_ENGINE_NAME = "SEQUEST";
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
        public clsPHRPParserSequest(string datasetName, string inputFilePath)
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
        public clsPHRPParserSequest(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.Sequest, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        /// <remarks></remarks>
        public clsPHRPParserSequest(string datasetName, string inputFilePath, clsPHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, clsPHRPReader.ePeptideHitResultType.Sequest, startupOptions)
        {
        }

        /// <summary>
        /// Define column header names
        /// </summary>
        protected override void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_HitNum);
            AddHeaderColumn(DATA_COLUMN_ScanNum);
            AddHeaderColumn(DATA_COLUMN_ScanCount);
            AddHeaderColumn(DATA_COLUMN_ChargeState);
            AddHeaderColumn(DATA_COLUMN_MH);
            AddHeaderColumn(DATA_COLUMN_XCorr);
            AddHeaderColumn(DATA_COLUMN_DelCn);
            AddHeaderColumn(DATA_COLUMN_Sp);
            AddHeaderColumn(DATA_COLUMN_Reference);
            AddHeaderColumn(DATA_COLUMN_MultiProtein);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_DelCn2);
            AddHeaderColumn(DATA_COLUMN_RankSp);
            AddHeaderColumn(DATA_COLUMN_RankXc);
            AddHeaderColumn(DATA_COLUMN_DelM);
            AddHeaderColumn(DATA_COLUMN_XcRatio);

            mColumnHeaders.Add(DATA_COLUMN_PassFilt, -1);
            mColumnHeaders.Add(DATA_COLUMN_MScore, -1);

            AddHeaderColumn(DATA_COLUMN_Ions_Observed);
            AddHeaderColumn(DATA_COLUMN_Ions_Expected);

            AddHeaderColumn(DATA_COLUMN_NumTrypticEnds);
            AddHeaderColumn(DATA_COLUMN_DelM_PPM);
        }

        /// <summary>
        /// Determines the precursor mass tolerance
        /// </summary>
        /// <param name="searchEngineParams"></param>
        /// <param name="tolerancePPM">Precursor mass tolerance, in ppm</param>
        /// <returns>Precursor tolerance, in Da</returns>
        /// <remarks></remarks>
        private double DeterminePrecursorMassTolerance(clsSearchEngineParameters searchEngineParams, out double tolerancePPM)
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

                        toleranceDa = clsPeptideMassCalculator.PPMToMass(value, 2000);
                    }
                    else if (units == 1)
                    {
                        // Tolerance is in milli mass units
                        toleranceDa = value / 1000.0;

                        // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                        tolerancePPM = clsPeptideMassCalculator.MassToPPM(toleranceDa, 2000);
                    }
                    else
                    {
                        // Tolerance is in daltons
                        toleranceDa = value;

                        // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                        tolerancePPM = clsPeptideMassCalculator.MassToPPM(toleranceDa, 2000);
                    }
                }
            }

            return toleranceDa;
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
        /// Defeault ModSummary file name for the given dataset
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
        /// Parses the specified Sequest parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out clsSearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new clsSearchEngineParameters(SEQ_SEARCH_ENGINE_NAME, mModInfo) {
                Enzyme = "trypsin"
            };

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, clsSearchEngineParameters searchEngineParams)
        {
            var reEnzymeSpecificity = new Regex(@"^\S+\s(\d)\s\d\s.+", RegexOptions.Compiled | RegexOptions.IgnoreCase);

            var success = false;

            try
            {
                var paramFilePath = Path.Combine(mInputDirectoryPath, searchEngineParamFileName);

                if (!File.Exists(paramFilePath))
                {
                    ReportError("Sequest param file not found: " + paramFilePath);
                }
                else
                {
                    using (var reader = new StreamReader(new FileStream(paramFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                    {
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
                                        fastaFilePath = Path.Combine("C:\\Database", Path.GetFileName(settingValue));
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
                                        searchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_AVERAGE;
                                    }
                                    else
                                    {
                                        // Monoisotopic mass
                                        searchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
                                    }

                                    break;
                                case "mass_type_fragment":
                                    if (settingValue == "0")
                                    {
                                        // Average mass
                                        searchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_AVERAGE;
                                    }
                                    else
                                    {
                                        // Monoisotopic mass
                                        searchEngineParams.PrecursorMassType = clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
                                    }

                                    break;
                                case "max_num_internal_cleavage_sites":
                                    if (int.TryParse(settingValue, out value))
                                    {
                                        searchEngineParams.MaxNumberInternalCleavages = value;
                                    }

                                    break;
                                case "enzyme_info":
                                    // Used in new-style sequest parameter files

                                    // Examples include:
                                    // Fully-tryptic:     Trypsin(KR) 1 1 KR -
                                    // Partially-tryptic: Trypsin(KR) 2 1 KR -
                                    // No-enzyme:         No_Enzyme(-) 0 0 - -
                                    //
                                    searchEngineParams.Enzyme = "trypsin";

                                    if (settingValue.StartsWith("no_enzyme", StringComparison.InvariantCultureIgnoreCase))
                                    {
                                        searchEngineParams.MinNumberTermini = 0;
                                    }
                                    else
                                    {
                                        // Parse out the cleavage specificity number
                                        // This is the first number after the closing parenthesis in the above examples
                                        var reMatch = reEnzymeSpecificity.Match(settingValue);
                                        if (reMatch.Success)
                                        {
                                            if (int.TryParse(reMatch.Groups[1].Value, out value))
                                            {
                                                searchEngineParams.MinNumberTermini = value;
                                            }
                                        }
                                    }

                                    break;
                                case "enzyme_number":
                                    // Used in old-style sequest parameter files
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
                                            switch (value)
                                            {
                                                case 1:
                                                    searchEngineParams.Enzyme = "trypsin";
                                                    break;
                                                case 2:
                                                    searchEngineParams.Enzyme = "trypsin_modified";
                                                    break;
                                                case 3:
                                                    searchEngineParams.Enzyme = "Chymotrypsin";
                                                    break;
                                                case 4:
                                                    searchEngineParams.Enzyme = "Chymotrypsin_modified";
                                                    break;
                                                case 5:
                                                    searchEngineParams.Enzyme = "Clostripain";
                                                    break;
                                                case 6:
                                                    searchEngineParams.Enzyme = "Cyanogen_Bromide";
                                                    break;
                                                case 7:
                                                    searchEngineParams.Enzyme = "IodosoBenzoate";
                                                    break;
                                                case 8:
                                                    searchEngineParams.Enzyme = "Proline_Endopept";
                                                    break;
                                                case 9:
                                                    searchEngineParams.Enzyme = "Staph_Protease";
                                                    break;
                                                case 10:
                                                    searchEngineParams.Enzyme = "Trypsin_K";
                                                    break;
                                                case 11:
                                                    searchEngineParams.Enzyme = "Trypsin_R";
                                                    break;
                                                case 12:
                                                    searchEngineParams.Enzyme = "GluC";
                                                    break;
                                                case 13:
                                                    searchEngineParams.Enzyme = "LysC";
                                                    break;
                                                case 14:
                                                    searchEngineParams.Enzyme = "AspN";
                                                    break;
                                                case 15:
                                                    searchEngineParams.Enzyme = "Elastase";
                                                    break;
                                                case 16:
                                                    searchEngineParams.Enzyme = "Elastase/Tryp/Chymo";
                                                    break;
                                                default:
                                                    searchEngineParams.Enzyme = "Unknown";
                                                    break;
                                            }
                                            searchEngineParams.MinNumberTermini = 2;
                                        }
                                    }
                                    break;
                            }
                        }
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM);
                    searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;

                    success = true;
                }
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
            }

            return success;
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
            var columns = line.Split('\t');

            var success = false;

            psm = new clsPSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_ScanNum, mColumnHeaders, -100);
                if (psm.ScanNumber == -100)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_HitNum, mColumnHeaders, 0);
                    psm.ScoreRank = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_RankXc, mColumnHeaders, 1);

                    var peptide = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Peptide, mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = Convert.ToInt16(clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_ChargeState, mColumnHeaders, 0));

                    var protein = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_Reference, mColumnHeaders);
                    psm.AddProtein(protein);

                    // Note that the MH value listed in Sequest files is not the precursor MH but is instead the theoretical (computed) MH of the peptide
                    // We'll update this value below using massErrorDa
                    // We'll further update this value using the ScanStatsEx data
                    var precursorMH = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_MH, mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMH, 1, 0);

                    psm.MassErrorDa = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM, mColumnHeaders);
                    if (double.TryParse(psm.MassErrorDa, out var massErrorDa))
                    {
                        // Adjust the precursor mass
                        psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMH - massErrorDa, 1, 0);
                    }

                    psm.MassErrorPPM = clsPHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, DATA_COLUMN_XCorr);
                    AddScore(psm, columns, DATA_COLUMN_DelCn);
                    AddScore(psm, columns, DATA_COLUMN_Sp);
                    AddScore(psm, columns, DATA_COLUMN_DelCn2);
                    AddScore(psm, columns, DATA_COLUMN_RankSp);
                    AddScore(psm, columns, DATA_COLUMN_RankXc);
                    AddScore(psm, columns, DATA_COLUMN_XcRatio);
                    AddScore(psm, columns, DATA_COLUMN_Ions_Observed);
                    AddScore(psm, columns, DATA_COLUMN_Ions_Expected);
                    AddScore(psm, columns, DATA_COLUMN_NumTrypticEnds);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the Sequest data file: " + ex.Message);
            }

            return success;
        }
    }
}
