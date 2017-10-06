using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader
{
    public class clsPHRPSeqMapReader
    {
        #region "Constants"
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

        #endregion

        #region "Module-wide variables"
        private readonly string mDatasetName;
        private readonly string mInputFolderPath;

        private readonly string mResultToSeqMapFilename;
        private readonly string mSeqToProteinMapFilename;
        private readonly string mSeqInfoFilename;
        private readonly string mPepToProteinMapFilename;

        private readonly clsPHRPReader.ePeptideHitResultType mPeptideHitResultType;

        private int mMaxProteinsPerSeqID;

        private string mErrorMessage = string.Empty;
        #endregion

        #region "Properties"

        public string DatasetName
        {
            get { return mDatasetName; }
        }

        public string ErrorMessage
        {
            get { return mErrorMessage; }
        }

        public string InputFolderPath
        {
            get { return mInputFolderPath; }
        }

        public int MaxProteinsPerSeqID
        {
            get { return mMaxProteinsPerSeqID; }
            set { mMaxProteinsPerSeqID = value; }
        }

        public clsPHRPReader.ePeptideHitResultType PeptideHitResultType
        {
            get { return mPeptideHitResultType; }
        }

        public string PepToProteinMapFilename
        {
            get { return mPepToProteinMapFilename; }
        }

        public string ResultToSeqMapFilename
        {
            get { return mResultToSeqMapFilename; }
        }

        public string SeqToProteinMapFilename
        {
            get { return mSeqToProteinMapFilename; }
        }

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFolderPath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit result type</param>
        /// <remarks></remarks>
        public clsPHRPSeqMapReader(string strDatasetName, string strInputFolderPath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType) : this(strDatasetName, strInputFolderPath, ePeptideHitResultType, clsPHRPReader.GetPHRPSynopsisFileName(ePeptideHitResultType, strDatasetName))
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strDatasetName">Dataset name</param>
        /// <param name="strInputFolderPath">Input file path</param>
        /// <param name="ePeptideHitResultType">Peptide Hit result type</param>
        /// <param name="strPHRPDataFileName">The base PHRP data file name; used when calling AutoSwitchToLegacyMSGFDBIfRequired and AutoSwitchToFHTIfRequired</param>
        /// <remarks></remarks>
        public clsPHRPSeqMapReader(string strDatasetName, string strInputFolderPath, clsPHRPReader.ePeptideHitResultType ePeptideHitResultType, string strPHRPDataFileName)
        {
            mDatasetName = strDatasetName;

            if (string.IsNullOrEmpty(mDatasetName))
            {
                mErrorMessage = "Dataset name cannot be empty";
                throw new Exception(mErrorMessage);
            }

            mInputFolderPath = strInputFolderPath;
            if (string.IsNullOrEmpty(mInputFolderPath))
            {
                mInputFolderPath = string.Empty;
            }

            mPeptideHitResultType = ePeptideHitResultType;

            mResultToSeqMapFilename = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mResultToSeqMapFilename))
            {
                mErrorMessage = "Unable to determine ResultToSeqMap filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }
            else
            {
                mResultToSeqMapFilename = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(mResultToSeqMapFilename, strPHRPDataFileName);
                mResultToSeqMapFilename = Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mResultToSeqMapFilename, strPHRPDataFileName));
            }

            mSeqToProteinMapFilename = clsPHRPReader.GetPHRPSeqToProteinMapFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mSeqToProteinMapFilename))
            {
                mErrorMessage = "Unable to determine SeqToProteinMap filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }
            else
            {
                mSeqToProteinMapFilename = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(mSeqToProteinMapFilename, strPHRPDataFileName);
                mSeqToProteinMapFilename = Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mSeqToProteinMapFilename, strPHRPDataFileName));
            }

            mSeqInfoFilename = clsPHRPReader.GetPHRPSeqInfoFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mSeqInfoFilename))
            {
                mErrorMessage = "Unable to determine SeqInfo filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }
            else
            {
                mSeqInfoFilename = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(mSeqInfoFilename, strPHRPDataFileName);
                mSeqInfoFilename = Path.GetFileName(clsPHRPReader.AutoSwitchToFHTIfRequired(mSeqInfoFilename, strPHRPDataFileName));
            }

            mPepToProteinMapFilename = clsPHRPReader.GetPHRPPepToProteinMapFileName(mPeptideHitResultType, mDatasetName);
            if (string.IsNullOrEmpty(mPepToProteinMapFilename))
            {
                mErrorMessage = "Unable to determine PepToProtMap filename for PeptideHitResultType: " + mPeptideHitResultType.ToString();
                throw new Exception(mErrorMessage);
            }

            mMaxProteinsPerSeqID = 0;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strInputFolderPath">Input folder path</param>
        /// <param name="strResultToSeqMapFilename">ResultToSeqMap filename</param>
        /// <param name="strSeqToProteinMapFilename"></param>
        /// <param name="strSeqInfoFilename">SeqInfo filename</param>
        /// <remarks></remarks>
        public clsPHRPSeqMapReader(string strInputFolderPath, string strResultToSeqMapFilename, string strSeqToProteinMapFilename, string strSeqInfoFilename)
        {
            mInputFolderPath = strInputFolderPath;

            if (string.IsNullOrEmpty(mInputFolderPath))
            {
                mInputFolderPath = string.Empty;
            }

            mPeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strResultToSeqMapFilename);
            if (mPeptideHitResultType == clsPHRPReader.ePeptideHitResultType.Unknown)
            {
                mErrorMessage = "Unable to auto-determine the PepthideHit result type based on filename " + strResultToSeqMapFilename;
                throw new Exception(mErrorMessage);
            }

            mDatasetName = clsPHRPReader.AutoDetermineDatasetName(strResultToSeqMapFilename);
            if (string.IsNullOrEmpty(mDatasetName))
            {
                mErrorMessage = "Unable to auto-determine the dataset name using filename '" + strResultToSeqMapFilename + "'";
                throw new Exception(mErrorMessage);
            }

            mResultToSeqMapFilename = strResultToSeqMapFilename;
            mSeqToProteinMapFilename = strSeqToProteinMapFilename;
            mSeqInfoFilename = strSeqInfoFilename;

            mMaxProteinsPerSeqID = 0;
        }

        /// <summary>
        /// Load the mapping between ResultID and Protein Name
        /// </summary>
        /// <param name="lstResultToSeqMap">ResultToSeqMap list (output); keys are ResultID, Values as SeqID</param>
        /// <param name="lstSeqToProteinMap">SeqToProteinMap list (output); keys are SeqID, Values are list of clsProteinInfo objects</param>
        /// <param name="lstSeqInfo">SeqInfo list (output); keys are SeqID, Values are seq details stored in clsSeqInfo objects</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks></remarks>
        public bool GetProteinMapping(
            out SortedList<int, int> lstResultToSeqMap,
            out SortedList<int, List<clsProteinInfo>> lstSeqToProteinMap,
            out SortedList<int, clsSeqInfo> lstSeqInfo)
        {
            Dictionary<string, clsPepToProteinMapInfo> lstPepToProteinMap = null;

            return GetProteinMapping(out lstResultToSeqMap, out lstSeqToProteinMap, out lstSeqInfo, out lstPepToProteinMap);
        }

        /// <summary>
        /// Load the mapping between ResultID and Protein Name
        /// </summary>
        /// <param name="lstResultToSeqMap">ResultToSeqMap list (output); keys are ResultID, Values as SeqID</param>
        /// <param name="lstSeqToProteinMap">SeqToProteinMap list (output); keys are SeqID, Values are list of clsProteinInfo objects</param>
        /// <param name="lstSeqInfo">SeqInfo list (output); keys are SeqID, Values are seq details stored in clsSeqInfo objects</param>
        /// <param name="lstPepToProteinMap">PepToProteinMap list (ouput); keys are clean peptide sequences (no mods), Values are Protein name and residue start/end locations for the peptide</param>
        /// <returns>True if success, false if an error</returns>
        /// <remarks></remarks>
        public bool GetProteinMapping(
            out SortedList<int, int> lstResultToSeqMap,
            out SortedList<int, List<clsProteinInfo>> lstSeqToProteinMap,
            out SortedList<int, clsSeqInfo> lstSeqInfo,
            out Dictionary<string, clsPepToProteinMapInfo> lstPepToProteinMap)
        {
            bool blnSuccess = false;
            string strFilePath = null;

            // Note: do not put a Try/Catch handler in this function
            //       Instead, allow LoadResultToSeqMapping or LoadSeqToProteinMapping to raise exceptions

            lstResultToSeqMap = new SortedList<int, int>();
            lstSeqToProteinMap = new SortedList<int, List<clsProteinInfo>>();
            lstSeqInfo = new SortedList<int, clsSeqInfo>();
            lstPepToProteinMap = new Dictionary<string, clsPepToProteinMapInfo>();

            if (string.IsNullOrEmpty(mResultToSeqMapFilename))
            {
                blnSuccess = false;
            }
            else
            {
                strFilePath = Path.Combine(mInputFolderPath, mResultToSeqMapFilename);
                if (!File.Exists(strFilePath))
                {
                    mErrorMessage = "SeqInfo file not found: " + strFilePath;
                    blnSuccess = false;
                }
                else
                {
                    blnSuccess = LoadResultToSeqMapping(strFilePath, lstResultToSeqMap);
                }
            }

            if (blnSuccess)
            {
                if (!string.IsNullOrEmpty(mSeqInfoFilename))
                {
                    strFilePath = Path.Combine(mInputFolderPath, mSeqInfoFilename);
                    if (File.Exists(strFilePath))
                    {
                        LoadSeqInfo(strFilePath, lstSeqInfo);
                    }
                }

                if (!string.IsNullOrEmpty(mSeqToProteinMapFilename))
                {
                    strFilePath = Path.Combine(mInputFolderPath, mSeqToProteinMapFilename);
                    if (!File.Exists(strFilePath))
                    {
                        mErrorMessage = "SeqInfo file not found: " + strFilePath;
                        blnSuccess = false;
                    }
                    else
                    {
                        blnSuccess = LoadSeqToProteinMapping(strFilePath, lstSeqToProteinMap);
                    }
                }
                else
                {
                    blnSuccess = false;
                }

                if (blnSuccess && !string.IsNullOrEmpty(mPepToProteinMapFilename))
                {
                    strFilePath = Path.Combine(mInputFolderPath, mPepToProteinMapFilename);
                    if (!File.Exists(strFilePath))
                    {
                        var strFilePathAlternate = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strFilePath, "Dataset_msgfdb.txt");
                        if (File.Exists(strFilePathAlternate))
                        {
                            strFilePath = strFilePathAlternate;
                        }
                    }

                    if (!File.Exists(strFilePath))
                    {
                        Console.WriteLine("Warning: PepToProtMap file not found; protein residue start/end values will be zero");
                        Console.WriteLine("         " + strFilePath);
                    }
                    else
                    {
                        blnSuccess = LoadPepToProtMapData(strFilePath, lstPepToProteinMap);
                    }
                }
            }

            return blnSuccess;
        }

        /// <summary>
        /// Load the Peptide to Protein mapping using the specified PHRP result file
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="lstPepToProteinMap">Peptide to protein mapping</param>
        /// <returns></returns>
        /// <remarks>The PepToProtMap file contains Residue_Start and Residue_End columns</remarks>
        private bool LoadPepToProtMapData(string strFilePath, Dictionary<string, clsPepToProteinMapInfo> lstPepToProteinMap)
        {
            var linesRead = 0;
            DateTime dtLastProgress = DateTime.UtcNow;
            bool blnNotifyComplete = false;

            try
            {
                // Read the data from the PepToProtMap file
                using (var srInFile = new StreamReader(new FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        var strLineIn = srInFile.ReadLine();
                        linesRead += 1;

                        if (!string.IsNullOrEmpty(strLineIn))
                        {
                            var strSplitLine = strLineIn.Split('\t');

                            if (strSplitLine.Length >= 4)
                            {
                                int residueStart = 0;
                                int residueEnd = 0;

                                // Parse out the numbers from the last two columns
                                // (the first line of the file is the header line, and it will get skipped)
                                if (int.TryParse(strSplitLine[2], out residueStart))
                                {
                                    if (int.TryParse(strSplitLine[3], out residueEnd))
                                    {
                                        var strPeptide = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strSplitLine[0], true);

                                        clsPepToProteinMapInfo oPepToProtMapInfo = null;

                                        if (lstPepToProteinMap.TryGetValue(strPeptide, out oPepToProtMapInfo))
                                        {
                                            if (mMaxProteinsPerSeqID == 0 || oPepToProtMapInfo.ProteinCount < mMaxProteinsPerSeqID)
                                            {
                                                oPepToProtMapInfo.AddProtein(strSplitLine[1], residueStart, residueEnd);
                                            }
                                        }
                                        else
                                        {
                                            oPepToProtMapInfo = new clsPepToProteinMapInfo(strSplitLine[1], residueStart, residueEnd);

                                            lstPepToProteinMap.Add(strPeptide, oPepToProtMapInfo);
                                        }
                                    }
                                }
                            }

                            if (linesRead % 100 == 0)
                            {
                                if (DateTime.UtcNow.Subtract(dtLastProgress).TotalSeconds >= 5)
                                {
                                    var pctComplete = srInFile.BaseStream.Position / Convert.ToDouble(srInFile.BaseStream.Length) * 100;
                                    Console.WriteLine(" ... caching PepToProtMapData: " + pctComplete.ToString("0.0") + "% complete");
                                    dtLastProgress = DateTime.UtcNow;
                                    blnNotifyComplete = true;
                                }
                            }
                        }
                    }
                }

                if (blnNotifyComplete)
                {
                    Console.WriteLine(" ... caching PepToProtMapData: 100% complete");
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Pep to Prot Map data from " + Path.GetFileName(strFilePath) + ": " + ex.Message);
            }

            return true;
        }

        /// <summary>
        /// Load the Result to Seq mapping using the specified PHRP result file
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="lstResultToSeqMap">Result to sequence mapping</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadResultToSeqMapping(string strFilePath, SortedList<int, int> lstResultToSeqMap)
        {
            string strLineIn = null;
            string[] strSplitLine = null;

            int intResultID = 0;
            int intSeqID = 0;

            try
            {
                // Read the data from the result to sequence map file
                using (var srInFile = new StreamReader(new FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        strLineIn = srInFile.ReadLine();

                        if (!string.IsNullOrEmpty(strLineIn))
                        {
                            strSplitLine = strLineIn.Split('\t');

                            if (strSplitLine.Length >= 2)
                            {
                                // Parse out the numbers from the first two columns
                                // (the first line of the file is the header line, and it will get skipped)
                                if (int.TryParse(strSplitLine[0], out intResultID))
                                {
                                    if (int.TryParse(strSplitLine[1], out intSeqID))
                                    {
                                        if (!lstResultToSeqMap.ContainsKey(intResultID))
                                        {
                                            lstResultToSeqMap.Add(intResultID, intSeqID);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Result to Seq Mapping from " + Path.GetFileName(strFilePath) + ": " + ex.Message);
            }

            return true;
        }

        /// <summary>
        /// Load the sequence info
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="lstSeqInfo">Sequences</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadSeqInfo(string strFilePath, SortedList<int, clsSeqInfo> lstSeqInfo)
        {
            SortedDictionary<string, int> objColumnHeaders = default(SortedDictionary<string, int>);

            string strLineIn = null;
            string[] strSplitLine = null;

            int intSeqID = 0;
            int intModCount = 0;
            string strModDescription = null;
            double dblMonoisotopicMass = 0;

            bool blnHeaderLineParsed = false;
            bool blnSkipLine = false;

            try
            {
                // Initialize the column mapping
                // Using a case-insensitive comparer
                objColumnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase);

                // Define the default column mapping
                objColumnHeaders.Add(SEQ_INFO_COLUMN_Unique_Seq_ID, 0);
                objColumnHeaders.Add(SEQ_INFO_COLUMN_Mod_Count, 1);
                objColumnHeaders.Add(SEQ_INFO_COLUMN_Mod_Description, 2);
                objColumnHeaders.Add(SEQ_INFO_COLUMN_Monoisotopic_Mass, 3);

                // Read the data from the sequence info file
                using (var srInFile = new StreamReader(new FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        strLineIn = srInFile.ReadLine();
                        blnSkipLine = false;

                        if (!string.IsNullOrEmpty(strLineIn))
                        {
                            strSplitLine = strLineIn.Split('\t');

                            if (!blnHeaderLineParsed)
                            {
                                if (strSplitLine[0].ToLower() == SEQ_INFO_COLUMN_Unique_Seq_ID.ToLower())
                                {
                                    // Parse the header line to confirm the column ordering
                                    clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders);
                                    blnSkipLine = true;
                                }

                                blnHeaderLineParsed = true;
                            }

                            if (!blnSkipLine && strSplitLine.Length >= 3)
                            {
                                if (int.TryParse(strSplitLine[0], out intSeqID))
                                {
                                    intModCount = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_INFO_COLUMN_Mod_Count, objColumnHeaders, 0);
                                    strModDescription = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_INFO_COLUMN_Mod_Description, objColumnHeaders, string.Empty);
                                    dblMonoisotopicMass = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_INFO_COLUMN_Monoisotopic_Mass, objColumnHeaders, 0.0);

                                    if (!lstSeqInfo.ContainsKey(intSeqID))
                                    {
                                        lstSeqInfo.Add(intSeqID, new clsSeqInfo(intSeqID, dblMonoisotopicMass, intModCount, strModDescription));
                                    }
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Seq Info from " + Path.GetFileName(strFilePath) + ": " + ex.Message);
            }

            return true;
        }

        /// <summary>
        /// Load the Sequence to Protein mapping using the specified PHRP result file
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="lstSeqToProteinMap">Sequence to protein map</param>
        /// <returns></returns>
        /// <remarks></remarks>
        private bool LoadSeqToProteinMapping(string strFilePath, SortedList<int, List<clsProteinInfo>> lstSeqToProteinMap)
        {
            List<clsProteinInfo> lstProteins = null;

            SortedDictionary<string, int> objColumnHeaders = default(SortedDictionary<string, int>);

            string strLineIn = null;
            string[] strSplitLine = null;

            clsProteinInfo objProteinInfo = default(clsProteinInfo);

            string strProteinName = null;
            int intSeqID = 0;
            clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants eCleavageState = default(clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants);
            clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants eTerminusState = default(clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants);

            bool blnHeaderLineParsed = false;
            bool blnSkipLine = false;

            try
            {
                // Initialize the column mapping
                // Using a case-insensitive comparer
                objColumnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase);

                // Define the default column mapping
                objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Unique_Seq_ID, 0);
                objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Cleavage_State, 1);
                objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Terminus_State, 2);
                objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_Name, 3);
                objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_EValue, 4);
                objColumnHeaders.Add(SEQ_PROT_MAP_COLUMN_Protein_Intensity, 5);

                // Read the data from the sequence to protein map file
                using (var srInFile = new StreamReader(new FileStream(strFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        strLineIn = srInFile.ReadLine();
                        blnSkipLine = false;

                        if (string.IsNullOrEmpty(strLineIn))
                        {
                            continue;
                        }

                        strSplitLine = strLineIn.Split('\t');

                        if (!blnHeaderLineParsed)
                        {
                            if (strSplitLine[0].ToLower() == SEQ_PROT_MAP_COLUMN_Unique_Seq_ID.ToLower())
                            {
                                // Parse the header line to confirm the column ordering
                                clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders);
                                blnSkipLine = true;
                            }

                            blnHeaderLineParsed = true;
                        }

                        if (blnSkipLine || strSplitLine.Length < 3)
                        {
                            continue;
                        }

                        if (!int.TryParse(strSplitLine[0], out intSeqID))
                        {
                            continue;
                        }

                        strProteinName = clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Protein_Name, objColumnHeaders, string.Empty);

                        if (string.IsNullOrEmpty(strProteinName))
                        {
                            continue;
                        }

                        eCleavageState = (clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants)clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Cleavage_State, objColumnHeaders, 0);
                        eTerminusState = (clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants)clsPHRPReader.LookupColumnValue(strSplitLine, SEQ_PROT_MAP_COLUMN_Terminus_State, objColumnHeaders, 0);

                        objProteinInfo = new clsProteinInfo(strProteinName, intSeqID, eCleavageState, eTerminusState);

                        if (lstSeqToProteinMap.TryGetValue(intSeqID, out lstProteins))
                        {
                            // Sequence already exists in lstSeqToProteinMap; add the new protein info
                            if (mMaxProteinsPerSeqID == 0 || lstProteins.Count < mMaxProteinsPerSeqID)
                            {
                                lstProteins.Add(objProteinInfo);
                            }
                        }
                        else
                        {
                            // New Sequence ID
                            lstProteins = new List<clsProteinInfo>();
                            lstProteins.Add(objProteinInfo);
                            lstSeqToProteinMap.Add(intSeqID, lstProteins);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Exception loading Seq to Protein Mapping from " + Path.GetFileName(strFilePath) + ": " + ex.Message);
            }

            return true;
        }
    }
}
