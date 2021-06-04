using System;
using System.Collections.Generic;
using System.IO;
using PHRPReader.Data;

// ReSharper disable UnusedMember.Global

namespace PHRPReader.Reader
{
    /// <summary>
    /// PHRP SeqMap reader
    /// </summary>
    public class PHRPSeqMapReader
    {
        // Ignore Spelling: PepToProtMap

#pragma warning disable 1591

        public const string SEQ_PROT_MAP_COLUMN_Unique_Seq_ID = "Unique_Seq_ID";
        public const string SEQ_PROT_MAP_COLUMN_Cleavage_State = "Cleavage_State";
        public const string SEQ_PROT_MAP_COLUMN_Terminus_State = "Terminus_State";
        public const string SEQ_PROT_MAP_COLUMN_Protein_Name = "Protein_Name";

        public const string SEQ_PROT_MAP_COLUMN_Protein_EValue = "Protein_Expectation_Value_Log(e)";      // Only used by X!Tandem
        public const string SEQ_PROT_MAP_COLUMN_Protein_Intensity = "Protein_Intensity_Log(I)";           // Only used by X!Tandem

        public const string SEQ_INFO_COLUMN_Unique_Seq_ID = "Unique_Seq_ID";
        public const string SEQ_INFO_COLUMN_Mod_Count = "Mod_Count";
        public const string SEQ_INFO_COLUMN_Mod_Description = "Mod_Description";
        public const string SEQ_INFO_COLUMN_Monoisotopic_Mass = "Monoisotopic_Mass";

#pragma warning restore 1591

        private readonly string mSeqInfoFilename;

        /// <summary>
        /// Dataset name
        /// </summary>
        public string DatasetName { get; }

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage { get; private set; } = string.Empty;

        /// <summary>
        /// Input directory path
        /// </summary>
        public string InputDirectoryPath { get; }

        /// <summary>
        /// Input directory path
        /// </summary>
        [Obsolete("Use InputDirectoryPath")]
        public string InputFolderPath => InputDirectoryPath;

        /// <summary>
        /// Max proteins to track for each SeqID
        /// </summary>
        public int MaxProteinsPerSeqID { get; set; }

        /// <summary>
        /// PHRP result type
        /// </summary>
        public PeptideHitResultTypes PeptideHitResultType { get; }

        /// <summary>
        /// PepToProtMap filename
        /// </summary>
        public string PepToProteinMapFilename { get; }

        /// <summary>
        /// ResultToSeqMap filename
        /// </summary>
        public string ResultToSeqMapFilename { get; }

        /// <summary>
        /// SeqToProteinMap filename
        /// </summary>
        public string SeqToProteinMapFilename { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputDirectoryPath">Input file path</param>
        /// <param name="peptideHitResultType">Peptide Hit result type</param>
        public PHRPSeqMapReader(string datasetName, string inputDirectoryPath, PeptideHitResultTypes peptideHitResultType)
            : this(datasetName, inputDirectoryPath, peptideHitResultType, ReaderFactory.GetPHRPSynopsisFileName(peptideHitResultType, datasetName))
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputDirectoryPath">Input file path</param>
        /// <param name="peptideHitResultType">Peptide Hit result type</param>
        /// <param name="phrpDataFileName">The base PHRP data file name; used when calling AutoSwitchToLegacyMSGFDBIfRequired and AutoSwitchToFHTIfRequired</param>
        public PHRPSeqMapReader(string datasetName, string inputDirectoryPath, PeptideHitResultTypes peptideHitResultType, string phrpDataFileName)
        {
            DatasetName = datasetName;

            if (string.IsNullOrEmpty(DatasetName))
            {
                ErrorMessage = "Dataset name cannot be empty";
                throw new Exception(ErrorMessage);
            }

            InputDirectoryPath = inputDirectoryPath;
            if (string.IsNullOrEmpty(InputDirectoryPath))
            {
                InputDirectoryPath = string.Empty;
            }

            PeptideHitResultType = peptideHitResultType;

            ResultToSeqMapFilename = ReaderFactory.GetPHRPResultToSeqMapFileName(peptideHitResultType, DatasetName);
            if (string.IsNullOrEmpty(ResultToSeqMapFilename))
            {
                ErrorMessage = "Unable to determine ResultToSeqMap filename for peptideHitResultType: " + peptideHitResultType.ToString();
                throw new Exception(ErrorMessage);
            }

            ResultToSeqMapFilename = ReaderFactory.FindPHRPFile(InputDirectoryPath, phrpDataFileName, ResultToSeqMapFilename, out _);

            SeqToProteinMapFilename = ReaderFactory.GetPHRPSeqToProteinMapFileName(peptideHitResultType, DatasetName);
            if (string.IsNullOrEmpty(SeqToProteinMapFilename))
            {
                ErrorMessage = "Unable to determine SeqToProteinMap filename for peptideHitResultType: " + peptideHitResultType.ToString();
                throw new Exception(ErrorMessage);
            }

            SeqToProteinMapFilename = ReaderFactory.FindPHRPFile(InputDirectoryPath, phrpDataFileName, SeqToProteinMapFilename, out _);

            mSeqInfoFilename = ReaderFactory.GetPHRPSeqInfoFileName(peptideHitResultType, DatasetName);
            if (string.IsNullOrEmpty(mSeqInfoFilename))
            {
                ErrorMessage = "Unable to determine SeqInfo filename for peptideHitResultType: " + peptideHitResultType.ToString();
                throw new Exception(ErrorMessage);
            }

            mSeqInfoFilename = ReaderFactory.FindPHRPFile(InputDirectoryPath, phrpDataFileName, mSeqInfoFilename, out _);

            PepToProteinMapFilename = ReaderFactory.GetPHRPPepToProteinMapFileName(peptideHitResultType, DatasetName);
            if (string.IsNullOrEmpty(PepToProteinMapFilename))
            {
                ErrorMessage = "Unable to determine PepToProtMap filename for PeptideHitResultType: " + peptideHitResultType.ToString();
                throw new Exception(ErrorMessage);
            }

            PepToProteinMapFilename = ReaderFactory.FindPHRPFile(InputDirectoryPath, phrpDataFileName, PepToProteinMapFilename, out _);

            MaxProteinsPerSeqID = 0;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="resultToSeqMapFilename">ResultToSeqMap filename</param>
        /// <param name="seqToProteinMapFilename"></param>
        /// <param name="seqInfoFilename">SeqInfo filename</param>
        public PHRPSeqMapReader(string inputDirectoryPath, string resultToSeqMapFilename, string seqToProteinMapFilename, string seqInfoFilename)
        {
            InputDirectoryPath = inputDirectoryPath;

            if (string.IsNullOrEmpty(InputDirectoryPath))
            {
                InputDirectoryPath = string.Empty;
            }

            PeptideHitResultType = ReaderFactory.AutoDetermineResultType(resultToSeqMapFilename);
            if (PeptideHitResultType == PeptideHitResultTypes.Unknown)
            {
                ErrorMessage = "Unable to auto-determine the PeptideHit result type based on filename " + resultToSeqMapFilename;
                throw new Exception(ErrorMessage);
            }

            DatasetName = ReaderFactory.AutoDetermineDatasetName(resultToSeqMapFilename);
            if (string.IsNullOrEmpty(DatasetName))
            {
                ErrorMessage = "Unable to auto-determine the dataset name using filename '" + resultToSeqMapFilename + "'";
                throw new Exception(ErrorMessage);
            }

            ResultToSeqMapFilename = resultToSeqMapFilename;
            SeqToProteinMapFilename = seqToProteinMapFilename;
            mSeqInfoFilename = seqInfoFilename;

            MaxProteinsPerSeqID = 0;
        }

        /// <summary>
        /// Load the mapping between ResultID and Protein Name
        /// </summary>
        /// <param name="resultToSeqMap">Output: ResultToSeqMap list; keys are ResultID, Values as SeqID</param>
        /// <param name="seqToProteinMap">Output: SeqToProteinMap list; keys are SeqID, Values are list of ProteinInfo objects</param>
        /// <param name="seqInfo">Output: SeqInfo list; keys are SeqID, Values are seq details stored in SeqInfo objects</param>
        /// <returns>True if successful, false if an error</returns>
        public bool GetProteinMapping(
            SortedList<int, int> resultToSeqMap,
            SortedList<int, List<ProteinInfo>> seqToProteinMap,
            SortedList<int, SequenceInfo> seqInfo)
        {
            var pepToProteinMap = new Dictionary<string, PepToProteinMapInfo>();
            return GetProteinMapping(resultToSeqMap, seqToProteinMap, seqInfo, pepToProteinMap);
        }

        /// <summary>
        /// Load the mapping between ResultID and Protein Name
        /// </summary>
        /// <param name="resultToSeqMap">Output: ResultToSeqMap list; keys are ResultID, Values as SeqID</param>
        /// <param name="seqToProteinMap">Output: SeqToProteinMap list; keys are SeqID, Values are list of ProteinInfo objects</param>
        /// <param name="seqInfo">Output: SeqInfo list; keys are SeqID, Values are seq details stored in SeqInfo objects</param>
        /// <param name="pepToProteinMap">Output: PepToProteinMap list; keys are clean peptide sequences (no mods), Values are Protein name and residue start/end locations for the peptide</param>
        /// <returns>True if successful, false if an error</returns>
        public bool GetProteinMapping(
            SortedList<int, int> resultToSeqMap,
            SortedList<int, List<ProteinInfo>> seqToProteinMap,
            SortedList<int, SequenceInfo> seqInfo,
            Dictionary<string, PepToProteinMapInfo> pepToProteinMap)
        {
            bool success;
            string filePath;

            resultToSeqMap.Clear();
            seqToProteinMap.Clear();
            seqInfo.Clear();
            pepToProteinMap.Clear();

            // Note: do not put a Try/Catch handler in this function
            //       Instead, allow LoadResultToSeqMapping or LoadSeqToProteinMapping to raise exceptions

            if (string.IsNullOrEmpty(ResultToSeqMapFilename))
            {
                success = false;
            }
            else
            {
                filePath = Path.Combine(InputDirectoryPath, ResultToSeqMapFilename);
                if (!File.Exists(filePath))
                {
                    ErrorMessage = "SeqInfo file not found: " + filePath;
                    success = false;
                }
                else
                {
                    success = LoadResultToSeqMapping(filePath, resultToSeqMap);
                }
            }

            if (!success)
                return false;

            if (!string.IsNullOrEmpty(mSeqInfoFilename))
            {
                filePath = Path.Combine(InputDirectoryPath, mSeqInfoFilename);
                if (File.Exists(filePath))
                {
                    LoadSeqInfo(filePath, seqInfo);
                }
            }

            if (!string.IsNullOrEmpty(SeqToProteinMapFilename))
            {
                filePath = Path.Combine(InputDirectoryPath, SeqToProteinMapFilename);
                if (!File.Exists(filePath))
                {
                    ErrorMessage = "SeqInfo file not found: " + filePath;
                    success = false;
                }
                else
                {
                    success = LoadSeqToProteinMapping(filePath, seqToProteinMap);
                }
            }
            else
            {
                success = false;
            }

            if (!success || string.IsNullOrEmpty(PepToProteinMapFilename))
                return success;

            filePath = Path.Combine(InputDirectoryPath, PepToProteinMapFilename);
            if (!File.Exists(filePath))
            {
                var filePathAlternate = ReaderFactory.AutoSwitchToLegacyMSGFDBIfRequired(filePath, "Dataset_msgfdb.txt");
                if (File.Exists(filePathAlternate))
                {
                    filePath = filePathAlternate;
                }
            }

            if (!File.Exists(filePath))
            {
                Console.WriteLine("Warning: PepToProtMap file not found; protein residue start/end values will be zero");
                Console.WriteLine("         " + filePath);
            }
            else
            {
                success = LoadPepToProtMapData(filePath, pepToProteinMap);
            }

            return success;
        }

        /// <summary>
        /// Load the Peptide to Protein mapping using the specified PHRP result file
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="pepToProteinMap">Peptide to protein mapping</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>The PepToProtMap file contains Residue_Start and Residue_End columns</remarks>
        private bool LoadPepToProtMapData(string filePath, IDictionary<string, PepToProteinMapInfo> pepToProteinMap)
        {
            var linesRead = 0;
            var lastProgress = DateTime.UtcNow;
            var notifyComplete = false;

            try
            {
                // Read the data from the PepToProtMap file
                using var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    linesRead++;

                    if (!string.IsNullOrEmpty(lineIn))
                    {
                        var splitLine = lineIn.Split('\t');

                        if (splitLine.Length >= 4)
                        {
                            // Parse out the numbers from the last two columns
                            // (the first line of the file is the header line, and it will get skipped)
                            if (int.TryParse(splitLine[2], out var residueStart))
                            {
                                if (int.TryParse(splitLine[3], out var residueEnd))
                                {
                                    var peptide = PeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(splitLine[0], true);

                                    if (pepToProteinMap.TryGetValue(peptide, out var pepToProtMapInfo))
                                    {
                                        if (MaxProteinsPerSeqID == 0 || pepToProtMapInfo.ProteinCount < MaxProteinsPerSeqID)
                                        {
                                            pepToProtMapInfo.AddProtein(splitLine[1], residueStart, residueEnd);
                                        }
                                    }
                                    else
                                    {
                                        pepToProtMapInfo = new PepToProteinMapInfo(splitLine[1], residueStart, residueEnd);

                                        pepToProteinMap.Add(peptide, pepToProtMapInfo);
                                    }
                                }
                            }
                        }

                        if (linesRead % 100 != 0)
                            continue;

                        if (DateTime.UtcNow.Subtract(lastProgress).TotalSeconds < 5)
                        {
                            continue;
                        }

                        var percentComplete = reader.BaseStream.Position / (float)reader.BaseStream.Length * 100;
                        Console.WriteLine(" ... caching PepToProtMapData: {0:0.0}% complete", percentComplete);
                        lastProgress = DateTime.UtcNow;
                        notifyComplete = true;
                    }
                }

                if (notifyComplete)
                {
                    Console.WriteLine(" ... caching PepToProtMapData: 100% complete");
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Peptide to Protein Map data from " + Path.GetFileName(filePath) + ": " + ex.Message);
            }

            return true;
        }

        /// <summary>
        /// Load the Result to Seq mapping using the specified PHRP result file
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="resultToSeqMap">Result to sequence mapping</param>
        /// <returns>True if successful, false if an error</returns>
        private bool LoadResultToSeqMapping(string filePath, IDictionary<int, int> resultToSeqMap)
        {
            try
            {
                // Read the data from the result to sequence map file
                using var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    var splitLine = lineIn.Split('\t');

                    if (splitLine.Length < 2)
                        continue;

                    // Parse out the numbers from the first two columns
                    // (the first line of the file is the header line, and it will get skipped)
                    if (int.TryParse(splitLine[0], out var resultID))
                    {
                        if (int.TryParse(splitLine[1], out var seqID))
                        {
                            if (!resultToSeqMap.ContainsKey(resultID))
                            {
                                resultToSeqMap.Add(resultID, seqID);
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Result to Seq Mapping from " + Path.GetFileName(filePath) + ": " + ex.Message);
            }

            return true;
        }

        /// <summary>
        /// Load the sequence info
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="seqInfo">Sequences</param>
        private void LoadSeqInfo(string filePath, IDictionary<int, SequenceInfo> seqInfo)
        {
            var headerLineParsed = false;

            try
            {
                // Initialize the column mapping
                // Using a case-insensitive comparer
                var columnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase)
                {
                    {SEQ_INFO_COLUMN_Unique_Seq_ID, 0},
                    {SEQ_INFO_COLUMN_Mod_Count, 1},
                    {SEQ_INFO_COLUMN_Mod_Description, 2},
                    {SEQ_INFO_COLUMN_Monoisotopic_Mass, 3}
                };

                // Read the data from the sequence info file
                using var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    var skipLine = false;

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    var splitLine = lineIn.Split('\t');

                    if (!headerLineParsed)
                    {
                        if (string.Equals(splitLine[0], SEQ_INFO_COLUMN_Unique_Seq_ID, StringComparison.OrdinalIgnoreCase))
                        {
                            // Parse the header line to confirm the column ordering
                            ReaderFactory.ParseColumnHeaders(splitLine, columnHeaders);
                            skipLine = true;
                        }

                        headerLineParsed = true;
                    }

                    if (!skipLine && splitLine.Length >= 3)
                    {
                        if (int.TryParse(splitLine[0], out var seqID))
                        {
                            var modCount = ReaderFactory.LookupColumnValue(splitLine, SEQ_INFO_COLUMN_Mod_Count, columnHeaders, 0);
                            var modDescription = ReaderFactory.LookupColumnValue(splitLine, SEQ_INFO_COLUMN_Mod_Description, columnHeaders, string.Empty);
                            var monoisotopicMass = ReaderFactory.LookupColumnValue(splitLine, SEQ_INFO_COLUMN_Monoisotopic_Mass, columnHeaders, 0.0);

                            if (!seqInfo.ContainsKey(seqID))
                            {
                                seqInfo.Add(seqID, new SequenceInfo(seqID, monoisotopicMass, modCount, modDescription));
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Seq Info from " + Path.GetFileName(filePath) + ": " + ex.Message);
            }
        }

        /// <summary>
        /// Load the Sequence to Protein mapping using the specified PHRP result file
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="seqToProteinMap">Sequence to protein map</param>
        /// <returns>True if successful, false if an error</returns>
        private bool LoadSeqToProteinMapping(string filePath, IDictionary<int, List<ProteinInfo>> seqToProteinMap)
        {
            var headerLineParsed = false;

            try
            {
                // Initialize the column mapping
                // Using a case-insensitive comparer
                var columnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase)
                {
                    {SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, 0},
                    {SEQ_PROT_MAP_COLUMN_Cleavage_State, 1},
                    {SEQ_PROT_MAP_COLUMN_Terminus_State, 2},
                    {SEQ_PROT_MAP_COLUMN_Protein_Name, 3},
                    {SEQ_PROT_MAP_COLUMN_Protein_EValue, 4},
                    {SEQ_PROT_MAP_COLUMN_Protein_Intensity, 5}
                };

                // Read the data from the sequence to protein map file
                using var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    var skipLine = false;

                    if (string.IsNullOrEmpty(lineIn))
                    {
                        continue;
                    }

                    var splitLine = lineIn.Split('\t');

                    if (!headerLineParsed)
                    {
                        if (string.Equals(splitLine[0], SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, StringComparison.OrdinalIgnoreCase))
                        {
                            // Parse the header line to confirm the column ordering
                            ReaderFactory.ParseColumnHeaders(splitLine, columnHeaders);
                            skipLine = true;
                        }

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 3)
                    {
                        continue;
                    }

                    if (!int.TryParse(splitLine[0], out var seqID))
                    {
                        continue;
                    }

                    var proteinName = ReaderFactory.LookupColumnValue(splitLine, SEQ_PROT_MAP_COLUMN_Protein_Name, columnHeaders, string.Empty);

                    if (string.IsNullOrEmpty(proteinName))
                    {
                        continue;
                    }

                    var cleavageState = (PeptideCleavageStateCalculator.PeptideCleavageState)ReaderFactory.LookupColumnValue(splitLine, SEQ_PROT_MAP_COLUMN_Cleavage_State, columnHeaders, 0);
                    var terminusState = (PeptideCleavageStateCalculator.PeptideTerminusState)ReaderFactory.LookupColumnValue(splitLine, SEQ_PROT_MAP_COLUMN_Terminus_State, columnHeaders, 0);

                    var proteinInfo = new ProteinInfo(proteinName, seqID, cleavageState, terminusState);

                    if (seqToProteinMap.TryGetValue(seqID, out var proteins))
                    {
                        // Sequence already exists in seqToProteinMap; add the new protein info
                        if (MaxProteinsPerSeqID == 0 || proteins.Count < MaxProteinsPerSeqID)
                        {
                            proteins.Add(proteinInfo);
                        }
                    }
                    else
                    {
                        // New Sequence ID
                        proteins = new List<ProteinInfo> {
                            proteinInfo
                        };
                        seqToProteinMap.Add(seqID, proteins);
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Seq to Protein Mapping from " + Path.GetFileName(filePath) + ": " + ex.Message);
            }

            return true;
        }
    }
}
