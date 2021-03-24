//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 11/28/2012
//
// This class parses data lines from MSAlign msalign_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP parser for MSAlign
    /// </summary>
    public class MSAlignSynFileReader : SynFileReaderBaseClass
    {
        #region "Constants"

#pragma warning disable 1591

        public const string DATA_COLUMN_ResultID = "ResultID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Prsm_ID = "Prsm_ID";
        public const string DATA_COLUMN_Spectrum_ID = "Spectrum_ID";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_PrecursorMZ = "PrecursorMZ";
        public const string DATA_COLUMN_DelM = "DelM";
        public const string DATA_COLUMN_DelM_PPM = "DelM_PPM";
        public const string DATA_COLUMN_MH = "MH";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Protein_Mass = "Protein_Mass";
        public const string DATA_COLUMN_Unexpected_Mod_Count = "Unexpected_Mod_Count";
        public const string DATA_COLUMN_Peak_Count = "Peak_Count";
        public const string DATA_COLUMN_Matched_Peak_Count = "Matched_Peak_Count";
        public const string DATA_COLUMN_Matched_Fragment_Ion_Count = "Matched_Fragment_Ion_Count";
        public const string DATA_COLUMN_PValue = "PValue";
        public const string DATA_COLUMN_Rank_PValue = "Rank_PValue";
        public const string DATA_COLUMN_EValue = "EValue";
        public const string DATA_COLUMN_FDR = "FDR";
        public const string DATA_COLUMN_Species_ID = "Species_ID";      // Only in MSAlign_Histone results
        public const string DATA_COLUMN_FragMethod = "FragMethod";      // Only in MSAlign_Histone results

        public const string FILENAME_SUFFIX_SYN = "_msalign_syn.txt";
        public const string FILENAME_SUFFIX_FHT = "_msalign_fht.txt";

        private const string MSAlign_SEARCH_ENGINE_NAME = "MSAlign";

        /// <summary>
        /// These columns correspond to the Synopsis file created by MSAlignResultsProcessor
        /// </summary>
        public enum MSAlignSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            Prsm_ID = 2,
            Spectrum_ID = 3,
            Charge = 4,
            PrecursorMZ = 5,
            DelM = 6,                            // Precursor error, in Da
            DelMPPM = 7,                         // Precursor error, in ppm
            MH = 8,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
            Peptide = 9,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
            Protein = 10,                        // Protein Name
            Protein_Mass = 11,
            Unexpected_Mod_Count = 12,
            Peak_Count = 13,
            Matched_Peak_Count = 14,
            Matched_Fragment_Ion_Count = 15,
            PValue = 16,
            Rank_PValue = 17,
            EValue = 18,
            FDR = 19,
            Species_ID = 20,
            FragMethod = 21
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
        public MSAlignSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MSAlignSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PHRPReader.PeptideHitResultTypes.MSAlign, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MSAlignSynFileReader(string datasetName, string inputFilePath, PHRPStartupOptions startupOptions)
            : base(datasetName, inputFilePath, PHRPReader.PeptideHitResultTypes.MSAlign, startupOptions)
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

            if (searchEngineParams.Parameters.TryGetValue("errorTolerance", out var tolerance))
            {
                // Parent mass tolerance, in ppm
                if (Double.TryParse(tolerance, out tolerancePPM))
                {
                    toleranceDa = PeptideMassCalculator.PPMToMass(tolerancePPM, 2000);
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
        public static SortedDictionary<string, MSAlignSynFileColumns> GetColumnHeaderNamesAndIDs()
        {
            var headerColumns = new SortedDictionary<string, MSAlignSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {DATA_COLUMN_ResultID, MSAlignSynFileColumns.ResultID},
                {DATA_COLUMN_Scan, MSAlignSynFileColumns.Scan},
                {DATA_COLUMN_Prsm_ID, MSAlignSynFileColumns.Prsm_ID},
                {DATA_COLUMN_Spectrum_ID, MSAlignSynFileColumns.Spectrum_ID},
                {DATA_COLUMN_Charge, MSAlignSynFileColumns.Charge},
                {DATA_COLUMN_PrecursorMZ, MSAlignSynFileColumns.PrecursorMZ},
                {DATA_COLUMN_DelM, MSAlignSynFileColumns.DelM},
                {DATA_COLUMN_DelM_PPM, MSAlignSynFileColumns.DelMPPM},
                {DATA_COLUMN_MH, MSAlignSynFileColumns.MH},
                {DATA_COLUMN_Peptide, MSAlignSynFileColumns.Peptide},
                {DATA_COLUMN_Protein, MSAlignSynFileColumns.Protein},
                {DATA_COLUMN_Protein_Mass, MSAlignSynFileColumns.Protein_Mass},
                {DATA_COLUMN_Unexpected_Mod_Count, MSAlignSynFileColumns.Unexpected_Mod_Count},
                {DATA_COLUMN_Peak_Count, MSAlignSynFileColumns.Peak_Count},
                {DATA_COLUMN_Matched_Peak_Count, MSAlignSynFileColumns.Matched_Peak_Count},
                {DATA_COLUMN_Matched_Fragment_Ion_Count, MSAlignSynFileColumns.Matched_Fragment_Ion_Count},
                {DATA_COLUMN_PValue, MSAlignSynFileColumns.PValue},
                {DATA_COLUMN_Rank_PValue, MSAlignSynFileColumns.Rank_PValue},
                {DATA_COLUMN_EValue, MSAlignSynFileColumns.EValue},
                {DATA_COLUMN_FDR, MSAlignSynFileColumns.FDR},
                {DATA_COLUMN_Species_ID, MSAlignSynFileColumns.Species_ID},         // Only in MSAlign_Histone results
                {DATA_COLUMN_FragMethod, MSAlignSynFileColumns.FragMethod}          // Only in MSAlign_Histone results
            };

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSAlignSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MSAlignSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
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
            return datasetName + FILENAME_SUFFIX_FHT;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_msalign_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msalign_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_msalign_syn_ProteinMods.txt";
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
            return datasetName + "_msalign_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_msalign_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_msalign_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MSAlign_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MSAlign parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MSAlign_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            var success = false;

            try
            {
                success = ReadKeyValuePairSearchEngineParamFile(MSAlign_SEARCH_ENGINE_NAME, searchEngineParamFileName, PHRPReader.PeptideHitResultTypes.MSAlign, searchEngineParams);

                if (success)
                {
                    // For MS-GF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                    // But MSAlign does not have a _ModDefs.txt file because it performs a blind search
                    // The user can define static mods on cysteine; check for these now

                    if (searchEngineParams.Parameters.TryGetValue("cysteineProtection", out var settingValue))
                    {
                        ModificationDefinition modDef;
                        switch (settingValue)
                        {
                            case "C57":
                                modDef = new ModificationDefinition(
                                    ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                    57.0215, "C",
                                    ModificationDefinition.ModificationTypeConstants.StaticMod,
                                    "IodoAcet");

                                searchEngineParams.AddModification(modDef);
                                break;

                            case "C58":
                                modDef = new ModificationDefinition(
                                    ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                                    58.0055, "C",
                                    ModificationDefinition.ModificationTypeConstants.StaticMod,
                                    "IodoAcid");

                                searchEngineParams.AddModification(modDef);
                                break;
                        }
                    }

                    // Determine the precursor mass tolerance (will store 0 if a problem or not found)
                    searchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(searchEngineParams, out var tolerancePPM);
                    searchEngineParams.PrecursorMassTolerancePpm = tolerancePPM;
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
                    psm.ScoreRank = PHRPReader.LookupColumnValue(columns, DATA_COLUMN_Rank_PValue, mColumnHeaders, 1);

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
                    AddScore(psm, columns, DATA_COLUMN_Prsm_ID);
                    AddScore(psm, columns, DATA_COLUMN_Spectrum_ID);

                    AddScore(psm, columns, DATA_COLUMN_MH);

                    AddScore(psm, columns, DATA_COLUMN_Protein_Mass);
                    AddScore(psm, columns, DATA_COLUMN_Unexpected_Mod_Count);
                    AddScore(psm, columns, DATA_COLUMN_Peak_Count);
                    AddScore(psm, columns, DATA_COLUMN_Matched_Peak_Count);
                    AddScore(psm, columns, DATA_COLUMN_Matched_Fragment_Ion_Count);

                    AddScore(psm, columns, DATA_COLUMN_PValue);
                    AddScore(psm, columns, DATA_COLUMN_EValue);
                    AddScore(psm, columns, DATA_COLUMN_FDR);

                    AddScore(psm, columns, DATA_COLUMN_Species_ID);
                    AddScore(psm, columns, DATA_COLUMN_FragMethod);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MSAlign data file: " + ex.Message);
            }

            return success;
        }
    }
}
