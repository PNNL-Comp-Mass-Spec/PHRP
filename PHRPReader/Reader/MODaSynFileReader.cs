//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/02/2014
//
// This class parses data lines from moda_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP parser for MODa
    /// </summary>
    public class MODaSynFileReader : SynFileReaderBaseClass
    {
        // Ignore Spelling: MODa

        #region "Constants"

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

        public const string FILENAME_SUFFIX_SYN = "_moda_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_moda_fht.txt";

        private const string MODa_SEARCH_ENGINE_NAME = "MODa";

        /// <summary>
        /// These columns correspond to the Synopsis file created by MODaResultsProcessor
        /// </summary>
        public enum MODaSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Spectrum_Index = 2,
            Charge = 3,
            PrecursorMZ = 4,
            DelM = 5,                            // Precursor error, in Da
            DelM_PPM = 6,                        // Precursor error, in ppm
            MH = 7,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
            Peptide = 8,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
            Protein = 9,                         // Protein Name
            Score = 10,
            Probability = 11,
            Rank_Probability = 12,
            Peptide_Position = 13,
            QValue = 14
        }

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
        public MODaSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MODaSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PHRPReader.PeptideHitResultTypes.MODa, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MODaSynFileReader(string datasetName, string inputFilePath, PHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, PHRPReader.PeptideHitResultTypes.MODa, startupOptions)
        {
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

            if (searchEngineParams.Parameters.TryGetValue("PPMTolerance", out var tolerance))
            {
                // Parent mass tolerance, in ppm
                if (double.TryParse(tolerance, out tolerancePPM))
                {
                    toleranceDa = PeptideMassCalculator.PPMToMass(tolerancePPM, 2000);
                }
            }

            // ReSharper disable once StringLiteralTypo
            else if (searchEngineParams.Parameters.TryGetValue("PeptTolerance", out tolerance))
            {
                // Parent mass tolerance, in Da
                double.TryParse(tolerance, out toleranceDa);

                // Convert from dalton to PPM (assuming a mass of 2000 m/z)
                tolerancePPM = PeptideMassCalculator.MassToPPM(toleranceDa, 2000);
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
        public static SortedDictionary<string, MODaSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            var headerColumns = new SortedDictionary<string, MODaSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {DATA_COLUMN_ResultID, MODaSynFileColumns.ResultID},
                {DATA_COLUMN_Scan, MODaSynFileColumns.Scan},
                {DATA_COLUMN_Spectrum_Index, MODaSynFileColumns.Spectrum_Index},
                {DATA_COLUMN_Charge, MODaSynFileColumns.Charge},
                {DATA_COLUMN_PrecursorMZ, MODaSynFileColumns.PrecursorMZ},
                {DATA_COLUMN_DelM, MODaSynFileColumns.DelM},
                {DATA_COLUMN_DelM_PPM, MODaSynFileColumns.DelM_PPM},
                {DATA_COLUMN_MH, MODaSynFileColumns.MH},
                {DATA_COLUMN_Peptide, MODaSynFileColumns.Peptide},
                {DATA_COLUMN_Protein, MODaSynFileColumns.Protein},
                {DATA_COLUMN_Score, MODaSynFileColumns.Score},
                {DATA_COLUMN_Probability, MODaSynFileColumns.Probability},
                {DATA_COLUMN_Rank_Probability, MODaSynFileColumns.Rank_Probability},
                {DATA_COLUMN_Peptide_Position, MODaSynFileColumns.Peptide_Position},
                {DATA_COLUMN_QValue, MODaSynFileColumns.QValue}
            };

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MODaSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MODaSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
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
            // MODa does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_moda_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_moda_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_moda_syn_ProteinMods.txt";
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
            return datasetName + "_moda_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_moda_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_moda_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MODa_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MODa parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MODa_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                var success = ReadKeyValuePairSearchEngineParamFile(MODa_SEARCH_ENGINE_NAME, searchEngineParamFileName, PHRPReader.PeptideHitResultTypes.MODa, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                // For MS-GF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                // But MODa does not have a _ModDefs.txt file because it performs a blind search
                // The user can define static mods on any of the residues, plus the peptide termini; check for these now

                var residuesToFind = new List<string> {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

                // This dictionary tracks the static mod names we will look for
                // It is populated using the amino acid letters in residuesToFind, plus also the N and T terminus tags
                var dctResiduesAndSymbols = new Dictionary<string, string>();

                foreach (var residueSymbol in residuesToFind)
                {
                    dctResiduesAndSymbols.Add(residueSymbol, residueSymbol);
                }

                dctResiduesAndSymbols.Add("ADD_NTerm", AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString());
                dctResiduesAndSymbols.Add("ADD_CTerm", AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString());

                foreach (var residueSpec in dctResiduesAndSymbols)
                {
                    var key = "ADD_" + residueSpec.Key;

                    if (!searchEngineParams.Parameters.TryGetValue(key, out var settingValue))
                        continue;

                    if (!double.TryParse(settingValue, out var modMassDa))
                        continue;

                    if (Math.Abs(modMassDa) < float.Epsilon)
                        continue;

                    var modType = Enums.ResidueModificationType.StaticMod;
                    if (residueSpec.Value == AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || residueSpec.Value == AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        modType = Enums.ResidueModificationType.TerminalPeptideStaticMod;
                    }

                    var modDef = new ModificationDefinition(
                        ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                        modMassDa,
                        residueSpec.Value,
                        modType,
                        "Mod" + modMassDa.ToString("0"));

                    searchEngineParams.AddModification(modDef);
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
        /// <param name="line">Data line</param>
        /// <param name="linesRead">Number of lines read so far (used for error reporting)</param>
        /// <param name="psm">Output: PSM details</param>
        /// <param name="fastReadMode">When set to true, reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
        public override bool ParsePHRPDataLine(string line, int linesRead, out PSM psm, bool fastReadMode)
        {
            const int SCAN_NOT_FOUND_FLAG = -100;

            var columns = line.Split('\t');

            var success = false;

            psm = new PSM();

            try
            {
                psm.DataLineText = line;
                psm.ScanNumber = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_Scan, mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_ResultID, mColumnHeaders, 0);
                    psm.ScoreRank = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_Rank_Probability, mColumnHeaders, 1);

                    var peptide = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_Peptide, mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = Convert.ToInt16(PHRPReader.LookupColumnValue(columns, DATA_COLUMN_Charge, mColumnHeaders, 0));

                    var protein = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_Protein, mColumnHeaders);
                    psm.AddProtein(protein);

                    var precursorMZ = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                    psm.MassErrorDa = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM, mColumnHeaders);
                    psm.MassErrorPPM = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_DelM_PPM, mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, DATA_COLUMN_Spectrum_Index);

                    AddScore(psm, columns, DATA_COLUMN_MH);

                    AddScore(psm, columns, DATA_COLUMN_Score);
                    AddScore(psm, columns, DATA_COLUMN_Probability);
                    AddScore(psm, columns, DATA_COLUMN_Peptide_Position);
                    AddScore(psm, columns, DATA_COLUMN_QValue);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MODa data file: " + ex.Message);
            }

            return success;
        }
    }
}
