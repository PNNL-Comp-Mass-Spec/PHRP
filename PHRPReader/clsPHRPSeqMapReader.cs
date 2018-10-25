using System;
using System.Collections.Generic;
using System.IO;

// ReSharper disable UnusedMember.Global

namespace PHRPReader
{
    /// <summary>
    /// PHRP SeqMap reader
    /// </summary>
    public class clsPHRPSeqMapReader
    {
        #region "Constants"

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

        #endregion

        #region "Module-wide variables"
        private readonly string mDatasetName;
        private readonly string mInputDirectoryPath;

        private readonly string mResultToSeqMapFilename;
        private readonly string mSeqToProteinMapFilename;
        private readonly string mSeqInfoFilename;
        private readonly string mPepToProteinMapFilename;

        private readonly clsPHRPReader.ePeptideHitResultType mPeptideHitResultType;

        private string mErrorMessage = string.Empty;
        #endregion

        #region "Properties"

        /// <summary>
        /// Dataset name
        /// </summary>
        public string DatasetName => mDatasetName;

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage => mErrorMessage;

        /// <summary>
        /// Input directory path
        /// </summary>
        public string InputDirectoryPath => mInputDirectoryPath;

        /// <summary>
        /// Input directory path
        /// </summary>
        [Obsolete("Use InputDirectoryPath")]
        public string InputFolderPath => mInputDirectoryPath;

        /// <summary>
        /// Max proteins to track for each SeqID
        /// </summary>
        public int MaxProteinsPerSeqID { get; set; }

        /// <summary>
        /// PHRP result type
        /// </summary>
        public clsPHRPReader.ePeptideHitResultType PeptideHitResultType => mPeptideHitResultType;

        /// <summary>
        /// PepToProtMap filename
        /// </summary>
        public string PepToProteinMapFilename => mPepToProteinMapFilename;

        /// <summary>
        /// ResultToSeqMap filename
        /// </summary>
        public string ResultToSeqMapFilename => mResultToSeqMapFilename;

        /// <summary>
        /// SeqToProteinMap filename
        /// </summary>
        public string SeqToProteinMapFilename => mSeqToProteinMapFilename;

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputDirectoryPath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit result type</param>
        /// <remarks></remarks>
        public clsPHRPSeqMapReader(string datasetName, string inputDirectoryPath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType) : this(datasetName, inputDirectoryPath, ePeptideHitResultType, clsPHRPReader.GetPHRPSynopsisFileName(ePeptideHitResultType, datasetName))
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="datasetName">Dataset name</param>
        /// <param name="inputDirectoryPath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit result type</param>
        /// <param name="phrpDataFileName">The base PHRP data file name; used when calling AutoSwitchToLegacyMSGFDBIfRequired and AutoSwitchToFHTIfRequired</param>
        /// <remarks></remarks>
        public clsPHRPSeqMapReader(string datasetName, string inputDirectoryPath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType, string phrpDataFileName)
        {
            mDatasetName = datasetName;

            if (string.IsNullOrEmpty(mDatasetName))
            {
                mErrorMessage = "Dataset name cannot be empty";
                throw new Exception(mErrorMessage);
            }

            mInputDirectoryPath = inputDirectoryPath;
            if (string.IsNullOrEmpty(mInputDirectoryPath))
            {
                mInputDirectoryPath = string.Empty;
            }

            mPeptideHitResultType = ePeptideHitResultType;

            mResultToSeqMapFilename = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mResultToSeqMapFilename))
            {
                mErrorMessage = "Unable to determine ResultToSeqMap filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }

            mResultToSeqMapFilename = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(mResultToSeqMapFilename, phrpDataFileName);
            mResultToSeqMapFilename = Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mResultToSeqMapFilename, phrpDataFileName));

            mSeqToProteinMapFilename = clsPHRPReader.GetPHRPSeqToProteinMapFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mSeqToProteinMapFilename))
            {
                mErrorMessage = "Unable to determine SeqToProteinMap filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }

            mSeqToProteinMapFilename = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(mSeqToProteinMapFilename, phrpDataFileName);
            mSeqToProteinMapFilename = Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mSeqToProteinMapFilename, phrpDataFileName));

            mSeqInfoFilename = clsPHRPReader.GetPHRPSeqInfoFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mSeqInfoFilename))
            {
                mErrorMessage = "Unable to determine SeqInfo filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }

            mSeqInfoFilename = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(mSeqInfoFilename, phrpDataFileName);
            mSeqInfoFilename = Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mSeqInfoFilename, phrpDataFileName));

            mPepToProteinMapFilename = clsPHRPReader.GetPHRPPepToProteinMapFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mPepToProteinMapFilename))
            {
                mErrorMessage = "Unable to determine PepToProtMap filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }

            MaxProteinsPerSeqID = 0;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="inputDirectoryPath">Input directory path</param>
        /// <param name="resultToSeqMapFilename">ResultToSeqMap filename</param>
        /// <param name="seqToProteinMapFilename"></param>
        /// <param name="seqInfoFilename">SeqInfo filename</param>
        /// <remarks></remarks>
        public clsPHRPSeqMapReader(string inputDirectoryPath, string resultToSeqMapFilename, string seqToProteinMapFilename, string seqInfoFilename)
        {
            mInputDirectoryPath = inputDirectoryPath;

            if (string.IsNullOrEmpty(mInputDirectoryPath))
            {
                mInputDirectoryPath = string.Empty;
            }

            mPeptideHitResultType = clsPHRPReader.AutoDetermineResultType(resultToSeqMapFilename);
            if (mPeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
            {
                mErrorMessage = "Unable to auto-determine the PepthideHit result type based on filename " + resultToSeqMapFilename;
                throw new Exception(mErrorMessage);
            }

            mDatasetName = clsPHRPReader.AutoDetermineDatasetName(resultToSeqMapFilename);
            if (string.IsNullOrEmpty(mDatasetName))
            {
                mErrorMessage = "Unable to auto-determine the dataset name using filename '" + resultToSeqMapFilename + "'";
                throw new Exception(mErrorMessage);
            }

            mResultToSeqMapFilename = resultToSeqMapFilename;
            mSeqToProteinMapFilename = seqToProteinMapFilename;
            mSeqInfoFilename = seqInfoFilename;

            MaxProteinsPerSeqID = 0;
        }

        /// <summary>
        /// Load the mapping between ResultID and Protein Name
        /// </summary>
        /// <param name="resultToSeqMap">ResultToSeqMap list (output); keys are ResultID, Values as SeqID</param>
        /// <param name="seqToProteinMap">SeqToProteinMap list (output); keys are SeqID, Values are list of clsProteinInfo objects</param>
        /// <param name="seqInfo">SeqInfo list (output); keys are SeqID, Values are seq details stored in clsSeqInfo objects</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks></remarks>
        public bool GetProteinMapping(
            out SortedList<int, int> resultToSeqMap,
            out SortedList<int, List<clsProteinInfo>> seqToProteinMap,
            out SortedList<int, clsSeqInfo> seqInfo)
        {
            return GetProteinMapping(out resultToSeqMap, out seqToProteinMap, out seqInfo, out _);
        }

        /// <summary>
        /// Load the mapping between ResultID and Protein Name
        /// </summary>
        /// <param name="resultToSeqMap">ResultToSeqMap list (output); keys are ResultID, Values as SeqID</param>
        /// <param name="seqToProteinMap">SeqToProteinMap list (output); keys are SeqID, Values are list of clsProteinInfo objects</param>
        /// <param name="seqInfo">SeqInfo list (output); keys are SeqID, Values are seq details stored in clsSeqInfo objects</param>
        /// <param name="pepToProteinMap">PepToProteinMap list (ouput); keys are clean peptide sequences (no mods), Values are Protein name and residue start/end locations for the peptide</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks></remarks>
        public bool GetProteinMapping(
            out SortedList<int, int> resultToSeqMap,
            out SortedList<int, List<clsProteinInfo>> seqToProteinMap,
            out SortedList<int, clsSeqInfo> seqInfo,
            out Dictionary<string, clsPepToProteinMapInfo> pepToProteinMap)
        {
            bool success;
            string filePath;

            // Note: do not put a Try/Catch handler in this function
            //       Instead, allow LoadResultToSeqMapping or LoadSeqToProteinMapping to raise exceptions

            resultToSeqMap = new SortedList<int, int>();
            seqToProteinMap = new SortedList<int, List<clsProteinInfo>>();
            seqInfo = new SortedList<int, clsSeqInfo>();
            pepToProteinMap = new Dictionary<string, clsPepToProteinMapInfo>();

            if (string.IsNullOrEmpty(mResultToSeqMapFilename))
            {
                success = false;
            }
            else
            {
                filePath = Path.Combine(mInputDirectoryPath, mResultToSeqMapFilename);
                if (!File.Exists(filePath))
                {
                    mErrorMessage = "SeqInfo file not found: " + filePath;
                    success = false;
                }
                else
                {
                    success = LoadResultToSeqMapping(filePath, resultToSeqMap);
                }
            }

            if (success)
            {
                if (!string.IsNullOrEmpty(mSeqInfoFilename))
                {
                    filePath = Path.Combine(mInputDirectoryPath, mSeqInfoFilename);
                    if (File.Exists(filePath))
                    {
                        LoadSeqInfo(filePath, seqInfo);
                    }
                }

                if (!string.IsNullOrEmpty(mSeqToProteinMapFilename))
                {
                    filePath = Path.Combine(mInputDirectoryPath, mSeqToProteinMapFilename);
                    if (!File.Exists(filePath))
                    {
                        mErrorMessage = "SeqInfo file not found: " + filePath;
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

                if (success && !string.IsNullOrEmpty(mPepToProteinMapFilename))
                {
                    filePath = Path.Combine(mInputDirectoryPath, mPepToProteinMapFilename);
                    if (!File.Exists(filePath))
                    {
                        var filePathAlternate = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(filePath, "Dataset_msgfdb.txt");
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
                }
            }

            return success;
        }

        /// <summary>
        /// Load the Peptide to Protein mapping using the specified PHRP result file
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="pepToProteinMap">Peptide to protein mapping</param>
        /// <returns></returns>
        /// <remarks>The PepToProtMap file contains Residue_Start and Residue_End columns</remarks>
        private bool LoadPepToProtMapData(string filePath, IDictionary<string, clsPepToProteinMapInfo> pepToProteinMap)
        {
            var linesRead = 0;
            var lastProgress = DateTime.UtcNow;
            var notifyComplete = false;

            try
            {
                // Read the data from the PepToProtMap file
                using (var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        linesRead += 1;

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
                                        var peptide = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(splitLine[0], true);

                                        if (pepToProteinMap.TryGetValue(peptide, out var pepToProtMapInfo))
                                        {
                                            if (MaxProteinsPerSeqID == 0 || pepToProtMapInfo.ProteinCount < MaxProteinsPerSeqID)
                                            {
                                                pepToProtMapInfo.AddProtein(splitLine[1], residueStart, residueEnd);
                                            }
                                        }
                                        else
                                        {
                                            pepToProtMapInfo = new clsPepToProteinMapInfo(splitLine[1], residueStart, residueEnd);

                                            pepToProteinMap.Add(peptide, pepToProtMapInfo);
                                        }
                                    }
                                }
                            }

                            if (linesRead % 100 == 0)
                            {
                                if (DateTime.UtcNow.Subtract(lastProgress).TotalSeconds >= 5)
                                {
                                    var pctComplete = reader.BaseStream.Position / Convert.ToDouble(reader.BaseStream.Length) * 100;
                                    Console.WriteLine(" ... caching PepToProtMapData: " + pctComplete.ToString("0.0") + "% complete");
                                    lastProgress = DateTime.UtcNow;
                                    notifyComplete = true;
                                }
                            }
                        }
                    }
                }

                if (notifyComplete)
                {
                    Console.WriteLine(" ... caching PepToProtMapData: 100% complete");
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Pep to Prot Map data from " + Path.GetFileName(filePath) + ": " + ex.Message);
            }

            return true;
        }

        /// <summary>
        /// Load the Result to Seq mapping using the specified PHRP result file
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="resultToSeqMap">Result to sequence mapping</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadResultToSeqMapping(string filePath, IDictionary<int, int> resultToSeqMap)
        {

            try
            {
                // Read the data from the result to sequence map file
                using (var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
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
        /// <returns></returns>
        /// <remarks></remarks>
        private void LoadSeqInfo(string filePath, IDictionary<int, clsSeqInfo> seqInfo)
        {

            var headerLineParsed = false;

            try
            {
                // Initialize the column mapping
                // Using a case-insensitive comparer
                var columnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase)
                {
                    {SEQ_INFO_COLUMN_Unique_Seq_ID, 0},
                    {SEQ_INFO_COLUMN_Mod_Count, 1},
                    {SEQ_INFO_COLUMN_Mod_Description, 2},
                    {SEQ_INFO_COLUMN_Monoisotopic_Mass, 3}
                };

                // Read the data from the sequence info file
                using (var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!reader.EndOfStream)
                    {
                        var lineIn = reader.ReadLine();
                        var skipLine = false;

                        if (string.IsNullOrEmpty(lineIn))
                            continue;

                        var splitLine = lineIn.Split('\t');

                        if (!headerLineParsed)
                        {
                            if (splitLine[0].ToLower() == SEQ_INFO_COLUMN_Unique_Seq_ID.ToLower())
                            {
                                // Parse the header line to confirm the column ordering
                                clsPHRPReader.ParseColumnHeaders(splitLine, columnHeaders);
                                skipLine = true;
                            }

                            headerLineParsed = true;
                        }

                        if (!skipLine && splitLine.Length >= 3)
                        {
                            if (int.TryParse(splitLine[0], out var seqID))
                            {
                                var modCount = clsPHRPReader.LookupColumnValue(splitLine, SEQ_INFO_COLUMN_Mod_Count, columnHeaders, 0);
                                var modDescription = clsPHRPReader.LookupColumnValue(splitLine, SEQ_INFO_COLUMN_Mod_Description, columnHeaders, string.Empty);
                                var monoisotopicMass = clsPHRPReader.LookupColumnValue(splitLine, SEQ_INFO_COLUMN_Monoisotopic_Mass, columnHeaders, 0.0);

                                if (!seqInfo.ContainsKey(seqID))
                                {
                                    seqInfo.Add(seqID, new clsSeqInfo(seqID, monoisotopicMass, modCount, modDescription));
                                }
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
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadSeqToProteinMapping(string filePath, IDictionary<int, List<clsProteinInfo>> seqToProteinMap)
        {
            var headerLineParsed = false;

            try
            {
                // Initialize the column mapping
                // Using a case-insensitive comparer
                var columnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase)
                {
                    {SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, 0},
                    {SEQ_PROT_MAP_COLUMN_Cleavage_State, 1},
                    {SEQ_PROT_MAP_COLUMN_Terminus_State, 2},
                    {SEQ_PROT_MAP_COLUMN_Protein_Name, 3},
                    {SEQ_PROT_MAP_COLUMN_Protein_EValue, 4},
                    {SEQ_PROT_MAP_COLUMN_Protein_Intensity, 5}
                };


                // Read the data from the sequence to protein map file
                using (var reader = new StreamReader(new FileStream(filePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
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
                            if (splitLine[0].ToLower() == SEQ_PROT_MAP_COLUMN_Unique_Seq_ID.ToLower())
                            {
                                // Parse the header line to confirm the column ordering
                                clsPHRPReader.ParseColumnHeaders(splitLine, columnHeaders);
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

                        var proteinName = clsPHRPReader.LookupColumnValue(splitLine, SEQ_PROT_MAP_COLUMN_Protein_Name, columnHeaders, string.Empty);

                        if (string.IsNullOrEmpty(proteinName))
                        {
                            continue;
                        }

                        var eCleavageState = (clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants)clsPHRPReader.LookupColumnValue(splitLine, SEQ_PROT_MAP_COLUMN_Cleavage_State, columnHeaders, 0);
                        var eTerminusState = (clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants)clsPHRPReader.LookupColumnValue(splitLine, SEQ_PROT_MAP_COLUMN_Terminus_State, columnHeaders, 0);

                        var proteinInfo = new clsProteinInfo(proteinName, seqID, eCleavageState, eTerminusState);

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
                            proteins = new List<clsProteinInfo> {
                                proteinInfo
                            };
                            seqToProteinMap.Add(seqID, proteins);
                        }
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
