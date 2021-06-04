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
    /// PHRP SynFileReader for MODa
    /// </summary>
    public class MODaSynFileReader : SynFileReaderBaseClass
    {
        // ReSharper disable once CommentTypo
        // Ignore Spelling: Da, moda, MODa, PeptTolerance, PepToProtMap

        /// <summary>
        /// MODa synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_moda_syn.txt";

        /// <summary>
        /// MODa first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for MaxQuant results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_moda_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string MODa_SEARCH_ENGINE_NAME = "MODa";

        /// <summary>
        /// Mapping from enum to synopsis file column name for MODa
        /// </summary>
        private static readonly Dictionary<MODaSynFileColumns, string> mSynopsisFileColumn = new();

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
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MODa, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MODaSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MODa, startupOptions)
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

                // Convert from Dalton to PPM (assuming a mass of 2000 m/z)
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
            return new(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", MODaSynFileColumns.ResultID },
                { "Scan", MODaSynFileColumns.Scan },
                { "Spectrum_Index", MODaSynFileColumns.Spectrum_Index },
                { "Charge", MODaSynFileColumns.Charge },
                { "PrecursorMZ", MODaSynFileColumns.PrecursorMZ },
                { "DelM", MODaSynFileColumns.DelM },
                { "DelM_PPM", MODaSynFileColumns.DelM_PPM },
                { "MH", MODaSynFileColumns.MH },
                { "Peptide", MODaSynFileColumns.Peptide },
                { "Protein", MODaSynFileColumns.Protein },
                { "Score", MODaSynFileColumns.Score },
                { "Probability", MODaSynFileColumns.Probability },
                { "Rank_Probability", MODaSynFileColumns.Rank_Probability },
                { "Peptide_Position", MODaSynFileColumns.Peptide_Position },
                { "QValue", MODaSynFileColumns.QValue }
            };
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
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MODaSynFileColumns column)
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
                var success = ReadKeyValuePairSearchEngineParamFile(MODa_SEARCH_ENGINE_NAME, searchEngineParamFileName, PeptideHitResultTypes.MODa, searchEngineParams);

                if (!success)
                {
                    return false;
                }

                // For MS-GF+ or SEQUEST we load mod info from the _ModDefs.txt file for the parameter file
                // But MODa does not have a _ModDefs.txt file because it performs a blind search
                // The user can define static mods on any of the residues, plus the peptide termini; check for these now

                var residuesToFind = new List<string> { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" };

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

                    var modType = ModificationDefinition.ResidueModificationType.StaticMod;
                    if (residueSpec.Value == AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || residueSpec.Value == AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        modType = ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod;
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
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                }
                else
                {
                    psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.ResultID), mColumnHeaders, 0);
                    psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.Rank_Probability), mColumnHeaders, 1);

                    var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.Peptide), mColumnHeaders);

                    if (fastReadMode)
                    {
                        psm.SetPeptide(peptide, updateCleanSequence: false);
                    }
                    else
                    {
                        psm.SetPeptide(peptide, mCleavageStateCalculator);
                    }

                    psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.Charge), mColumnHeaders, 0);

                    var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.Protein), mColumnHeaders);
                    psm.AddProtein(protein);

                    var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);
                    psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                    psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.DelM), mColumnHeaders);
                    psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODaSynFileColumns.DelM_PPM), mColumnHeaders);

                    success = true;
                }

                if (success)
                {
                    if (!fastReadMode)
                    {
                        UpdatePSMUsingSeqInfo(psm);
                    }

                    // Store the remaining scores
                    AddScore(psm, columns, GetColumnNameByID(MODaSynFileColumns.Spectrum_Index));

                    AddScore(psm, columns, GetColumnNameByID(MODaSynFileColumns.MH));

                    AddScore(psm, columns, GetColumnNameByID(MODaSynFileColumns.Score));
                    AddScore(psm, columns, GetColumnNameByID(MODaSynFileColumns.Probability));
                    AddScore(psm, columns, GetColumnNameByID(MODaSynFileColumns.Peptide_Position));
                    AddScore(psm, columns, GetColumnNameByID(MODaSynFileColumns.QValue));
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
