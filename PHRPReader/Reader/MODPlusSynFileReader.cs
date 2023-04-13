//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 05/18/2015
//
// This class parses data lines from modp_syn.txt files
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Xml;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SynFileReader for MOD+
    /// </summary>
    public class MODPlusSynFileReader : SynFileReaderBaseClass
    {
        // ReSharper disable once CommentTypo
        // Ignore Spelling: da, Daltons, massdiff, msms, MODp, modp, prot, tol

        /// <summary>
        /// MODPlus synopsis file suffix
        /// </summary>
        public const string FILENAME_SUFFIX_SYN = "_modp_syn.txt";

        /// <summary>
        /// MODPlus first hits file suffix
        /// </summary>
        /// <remarks>
        /// <para>
        /// This public constant is defined for compatibility with other classes
        /// </para>
        /// <para>
        /// We don't actually create first-hits files for MaxQuant results
        /// </para>
        /// </remarks>
        public const string FILENAME_SUFFIX_FHT = "_modp_fht.txt";

        /// <summary>
        /// Search engine name
        /// </summary>
        private const string MODPlus_SEARCH_ENGINE_NAME = "MODPlus";

        /// <summary>
        /// Mapping from enum to synopsis file column name for MODPlus
        /// </summary>
        private static readonly Dictionary<MODPlusSynFileColumns, string> mSynopsisFileColumn = new();

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
        public MODPlusSynFileReader(string datasetName, string inputFilePath)
            : this(datasetName, inputFilePath, loadModsAndSeqInfo: true)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="loadModsAndSeqInfo">If True, load the ModSummary file and SeqInfo files</param>
        public MODPlusSynFileReader(string datasetName, string inputFilePath, bool loadModsAndSeqInfo)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MODPlus, loadModsAndSeqInfo)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputFilePath">Input file path</param>
        /// <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and MaxProteinsPerPSM</param>
        public MODPlusSynFileReader(string datasetName, string inputFilePath, StartupOptions startupOptions)
            : base(datasetName, inputFilePath, PeptideHitResultTypes.MODPlus, startupOptions)
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
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, MODPlusSynFileColumns> GetColumnHeaderNamesAndIDs(bool includeLegacyNames)
        {
            var headerColumns = new SortedDictionary<string, MODPlusSynFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ResultID", MODPlusSynFileColumns.ResultID },
                { "Scan", MODPlusSynFileColumns.Scan },
                { "Spectrum_Index", MODPlusSynFileColumns.Spectrum_Index },
                { "Charge", MODPlusSynFileColumns.Charge },
                { "PrecursorMZ", MODPlusSynFileColumns.PrecursorMZ },
                { "DelM", MODPlusSynFileColumns.DelM },
                { "DelM_PPM", MODPlusSynFileColumns.DelM_PPM },
                { "MH", MODPlusSynFileColumns.MH },
                { "Peptide", MODPlusSynFileColumns.Peptide },
                { "NTT", MODPlusSynFileColumns.NTT },
                { "ModificationAnnotation", MODPlusSynFileColumns.ModificationAnnotation },
                { "Protein", MODPlusSynFileColumns.Protein },
                { "Peptide_Position", MODPlusSynFileColumns.Peptide_Position },
                { "Score", MODPlusSynFileColumns.Score },
                { "Probability", MODPlusSynFileColumns.Probability },
                { "Rank_Score", MODPlusSynFileColumns.Rank_Score },
                { "QValue", MODPlusSynFileColumns.QValue }
            };

            if (includeLegacyNames)
            {
                headerColumns.Add("Modification_Annotation", MODPlusSynFileColumns.ModificationAnnotation);
            }

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MODPlusSynFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MODPlusSynFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs(true);
            return GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the synopsis file column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MODPlusSynFileColumns column)
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
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPFirstHitsFileName(string datasetName)
        {
            // MODPlus does not have a first-hits file; just the _syn.txt file
            return string.Empty;
        }

        /// <summary>
        /// Default ModSummary file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPModSummaryFileName(string datasetName)
        {
            return datasetName + "_modp_syn_ModSummary.txt";
        }

        /// <summary>
        /// Default PepToProtMap file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPPepToProteinMapFileName(string datasetName)
        {
            return datasetName + "_modp_PepToProtMapMTS.txt";
        }

        /// <summary>
        /// Default ProteinMods file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPProteinModsFileName(string datasetName)
        {
            return datasetName + "_modp_syn_ProteinMods.txt";
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
            return datasetName + "_modp_syn_ResultToSeqMap.txt";
        }

        /// <summary>
        /// Default SeqInfo map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqInfoFileName(string datasetName)
        {
            return datasetName + "_modp_syn_SeqInfo.txt";
        }

        /// <summary>
        /// Default SeqToProtein map file name for the given dataset
        /// </summary>
        /// <param name="datasetName"></param>
        /// <returns>Filename</returns>
        public static string GetPHRPSeqToProteinMapFileName(string datasetName)
        {
            return datasetName + "_modp_syn_SeqToProteinMap.txt";
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public static string GetSearchEngineName()
        {
            return MODPlus_SEARCH_ENGINE_NAME;
        }

        /// <summary>
        /// Parses the specified MODp parameter file
        /// </summary>
        /// <param name="searchEngineParamFileName"></param>
        /// <param name="searchEngineParams"></param>
        /// <returns>True if successful, false if an error</returns>
        public override bool LoadSearchEngineParameters(string searchEngineParamFileName, out SearchEngineParameters searchEngineParams)
        {
            searchEngineParams = new SearchEngineParameters(MODPlus_SEARCH_ENGINE_NAME);

            var success = ReadSearchEngineParamFile(searchEngineParamFileName, searchEngineParams);

            ReadSearchEngineVersion(mPeptideHitResultType, searchEngineParams);

            return success;
        }

        private bool ReadSearchEngineParamFile(string searchEngineParamFileName, SearchEngineParameters searchEngineParams)
        {
            try
            {
                // For MS-GF+ or SEQUEST we load mod info from the _ModDefs.txt file for the parameter file
                // But MODPlus does not have a _ModDefs.txt file because it performs a blind search
                // The user can define static mods on any of the residues, plus the peptide termini; check for these now

                var paramFilePath = Path.Combine(InputDirectoryPath, searchEngineParamFileName);

                // Read the contents of the parameter file
                var doc = new XmlDocument();
                doc.Load(paramFilePath);

                var dbNodes = doc.SelectNodes("/search/database");

                if (dbNodes?.Count > 0)
                {
                    searchEngineParams.FastaFilePath = GetAttribute(dbNodes[0], "local_path");
                }

                var enzymeNodes = doc.SelectNodes("/search/enzyme_rule");

                if (enzymeNodes?.Count > 0)
                {
                    searchEngineParams.Enzyme = GetAttribute(enzymeNodes[0], "name");
                }

                var instrumentResolutionNodes = doc.SelectNodes("/search/instrument_resolution");

                if (instrumentResolutionNodes?.Count > 0)
                {
                    searchEngineParams.PrecursorMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes[0], "ms"));

                    // ReSharper disable once StringLiteralTypo
                    searchEngineParams.FragmentMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes[0], "msms"));
                }

                var enzymeConstraintNodes = doc.SelectNodes("/search/parameters/enzyme_constraint");

                if (enzymeConstraintNodes?.Count > 0)
                {
                    if (int.TryParse(GetAttribute(enzymeConstraintNodes[0], "max_miss_cleavages"), out var maxNumberInternalCleavages))
                    {
                        searchEngineParams.MaxNumberInternalCleavages = maxNumberInternalCleavages;
                    }
                    if (int.TryParse(GetAttribute(enzymeConstraintNodes[0], "min_number_termini"), out var minNumberTermini))
                    {
                        searchEngineParams.MinNumberTermini = minNumberTermini;
                    }
                }

                var massTolParamNode = doc.SelectNodes("/search/parameters/peptide_mass_tol");

                if (massTolParamNode?.Count > 0)
                {
                    var tolerance = GetAttribute(massTolParamNode[0], "value");
                    var massUnits = GetAttribute(massTolParamNode[0], "unit");

                    if (string.Equals(massUnits, "ppm", StringComparison.OrdinalIgnoreCase))
                    {
                        // Parent mass tolerance, in ppm

                        if (double.TryParse(tolerance, out var precursorMassTolerancePpm))
                        {
                            searchEngineParams.PrecursorMassTolerancePpm = precursorMassTolerancePpm;
                            searchEngineParams.PrecursorMassToleranceDa = PeptideMassCalculator.PPMToMass(searchEngineParams.PrecursorMassTolerancePpm, 2000);
                        }
                    }
                    else if (string.Equals(massUnits, "da", StringComparison.OrdinalIgnoreCase))
                    {
                        if (double.TryParse(tolerance, out var precursorMassToleranceDa))
                        {
                            searchEngineParams.PrecursorMassToleranceDa = precursorMassToleranceDa;
                        }

                        // Convert from Daltons to PPM (assuming a mass of 2000 m/z)
                        searchEngineParams.PrecursorMassToleranceDa = PeptideMassCalculator.MassToPPM(searchEngineParams.PrecursorMassToleranceDa, 2000);
                    }
                }

                var modNodes = doc.SelectNodes("/search/modifications/fixed/mod");

                if (modNodes == null || modNodes.Count <= 0)
                    return true;

                // Store the fixed mods

                foreach (XmlNode node in modNodes)
                {
                    // var modName = GetAttribute(node, "name");
                    var residue = GetAttribute(node, "site").Trim();
                    // var modPosition = GetAttribute(node, "position");

                    // ReSharper disable once StringLiteralTypo
                    var modMass = GetAttribute(node, "massdiff");

                    // Replace N-Term or C-Term with < or >
                    if (string.Equals(residue, "n-term", StringComparison.OrdinalIgnoreCase))
                        residue = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                    if (string.Equals(residue, "c-term", StringComparison.OrdinalIgnoreCase))
                        residue = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();

                    if (!double.TryParse(modMass, out var modMassDa))
                        continue;

                    if (Math.Abs(modMassDa - 0) < float.Epsilon)
                        continue;

                    var modType = ModificationDefinition.ResidueModificationType.StaticMod;

                    if (residue == AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() || residue == AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        modType = ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod;
                    }

                    var modDef = new ModificationDefinition(
                        ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL,
                        modMassDa,
                        residue,
                        modType,
                        "Mod" + modMassDa.ToString("0"));

                    searchEngineParams.AddModification(modDef);
                }

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error in ReadSearchEngineParamFile: " + ex.Message);
                return false;
            }
        }

        private string ConvertResolutionModeToMassType(string resolutionType)
        {
            if (resolutionType == "high")
            {
                return SearchEngineParameters.MASS_TYPE_MONOISOTOPIC;
            }

            return resolutionType == "low" ? SearchEngineParameters.MASS_TYPE_AVERAGE : "unknown";
        }

        private string GetAttribute(XmlNode node, string attributeName)
        {
            if (node.Attributes?.Count > 0)
            {
                try
                {
                    var attribute = node.Attributes[attributeName];

                    if (attribute != null)
                    {
                        return attribute.Value;
                    }
                }
                catch (Exception ex)
                {
                    // Ignore errors
                    Console.WriteLine("Attribute lookup error: " + ex.Message);
                }
            }

            return string.Empty;
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
                psm.ScanNumber = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.Scan), mColumnHeaders, SCAN_NOT_FOUND_FLAG);
                if (psm.ScanNumber == SCAN_NOT_FOUND_FLAG)
                {
                    // Data line is not valid
                    return false;
                }

                psm.ResultID = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.ResultID), mColumnHeaders, 0);
                psm.ScoreRank = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.Rank_Score), mColumnHeaders, 1);

                var peptide = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.Peptide), mColumnHeaders);

                if (fastReadMode)
                {
                    psm.SetPeptide(peptide, updateCleanSequence: false);
                }
                else
                {
                    psm.SetPeptide(peptide, mCleavageStateCalculator);
                }

                psm.Charge = (short)ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.Charge), mColumnHeaders, 0);

                var protein = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.Protein), mColumnHeaders);
                psm.AddProtein(protein);

                var precursorMZ = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.PrecursorMZ), mColumnHeaders, 0.0);
                psm.PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(precursorMZ, psm.Charge, 0);

                psm.MassErrorDa = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.DelM), mColumnHeaders);
                psm.MassErrorPPM = ReaderFactory.LookupColumnValue(columns, GetColumnNameByID(MODPlusSynFileColumns.DelM_PPM), mColumnHeaders);

                if (!fastReadMode)
                {
                    UpdatePSMUsingSeqInfo(psm);
                }

                // Store the remaining scores
                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.Spectrum_Index));

                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.MH));

                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.ModificationAnnotation));
                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.Peptide_Position));

                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.Score));
                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.Probability));
                AddScore(psm, columns, GetColumnNameByID(MODPlusSynFileColumns.QValue));

                return true;
            }
            catch (Exception ex)
            {
                ReportError("Error parsing line " + linesRead + " in the MODPlus data file: " + ex.Message);
                return false;
            }
        }
    }
}
