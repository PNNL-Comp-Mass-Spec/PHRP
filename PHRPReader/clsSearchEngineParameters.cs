using System;
using System.Collections.Generic;

namespace PHRPReader
{
    /// <summary>
    /// Search engine parameters container
    /// </summary>
    public class clsSearchEngineParameters
    {
        /// <summary>
        /// Monoisotopic mass
        /// </summary>
        public const string MASS_TYPE_MONOISOTOPIC = "monoisotopic";

        /// <summary>
        /// Average mass
        /// </summary>
        public const string MASS_TYPE_AVERAGE = "average";

        private string mSearchEngineName;
        private string mSearchEngineVersion;
        private DateTime mSearchDate;

        private string mPrecursorMassType;             // Typically "monoisotopic" or "average"
        private string mFragmentMassType;

        private string mSearchEngineParamFilePath;

        private string mEnzyme;

        private readonly List<clsModificationDefinition> mModInfo;

        private readonly Dictionary<string, string> mParameters;

        #region "Properties"
        /// <summary>
        /// Enzyme name
        /// </summary>
        public string Enzyme
        {
            get => mEnzyme;
            set
            {
                if (string.IsNullOrEmpty(value))
                    value = "none";
                mEnzyme = value;
            }
        }

        /// <summary>
        /// FASTA file path
        /// </summary>
        public string FastaFilePath { get; set; }

        /// <summary>
        /// Fragment mass type
        /// </summary>
        /// <remarks>Typically "monoisotopic" or "average"</remarks>
        public string FragmentMassType
        {
            get => mFragmentMassType;
            set
            {
                if (string.IsNullOrEmpty(value))
                    value = MASS_TYPE_MONOISOTOPIC;
                mFragmentMassType = value;
            }
        }

        /// <summary>
        /// Maximum number of internal cleavages (missed cleavage points)
        /// </summary>
        public int MaxNumberInternalCleavages { get; set; }

        /// <summary>
        /// 0 means no-enzyme, 1 means partially tryptic, 2 means fully tryptic
        /// </summary>
        /// <remarks>For trypsin, this is NTT or Number of Tryptic Termini</remarks>
        public int MinNumberTermini { get; set; }

        /// <summary>
        /// Dynamic and static mods to search for
        /// </summary>
        public List<clsModificationDefinition> ModInfo => mModInfo;

        /// <summary>
        /// Parameter dictionary (key/value pairs)
        /// </summary>
        public Dictionary<string, string> Parameters => mParameters;

        /// <summary>
        /// Precursor mass tolerance, in Daltons; 0 if unknown
        /// </summary>
        public double PrecursorMassToleranceDa { get; set; }

        /// <summary>
        /// Precursor mass tolerance, in ppm; 0 if unknown
        /// </summary>
        public double PrecursorMassTolerancePpm { get; set; }

        /// <summary>
        /// Precursor mass type
        /// </summary>
        /// <remarks>Typically "monoisotopic" or "average"</remarks>
        public string PrecursorMassType
        {
            get => mPrecursorMassType;
            set
            {
                if (string.IsNullOrWhiteSpace(value))
                    value = MASS_TYPE_MONOISOTOPIC;
                mPrecursorMassType = value;
            }
        }

        /// <summary>
        /// Search engine name
        /// </summary>
        public string SearchEngineName => mSearchEngineName;

        /// <summary>
        /// Search engine parameter file path
        /// </summary>
        public string SearchEngineParamFilePath => string.IsNullOrWhiteSpace(mSearchEngineParamFilePath) ? string.Empty : mSearchEngineParamFilePath;

        /// <summary>
        /// Search engine version
        /// </summary>
        public string SearchEngineVersion => mSearchEngineVersion;

        /// <summary>
        /// Search date
        /// </summary>
        public DateTime SearchDate => mSearchDate;

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="searchEngineName"></param>
        public clsSearchEngineParameters(string searchEngineName) : this(searchEngineName, new List<clsModificationDefinition>(), null)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="searchEngineName"></param>
        /// <param name="modInfo"></param>
        public clsSearchEngineParameters(string searchEngineName, List<clsModificationDefinition> modInfo) : this(searchEngineName, modInfo, null)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="searchEngineName"></param>
        /// <param name="modInfo"></param>
        /// <param name="parameters"></param>
        public clsSearchEngineParameters(string searchEngineName, List<clsModificationDefinition> modInfo, Dictionary<string, string> parameters)
        {
            InitializeDefaults();

            mSearchEngineName = searchEngineName;
            mSearchEngineParamFilePath = string.Empty;

            mModInfo = modInfo;

            mModInfo = modInfo ?? new List<clsModificationDefinition>();

            if (mParameters == null)
            {
                mParameters = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            }
            else
            {
                mParameters = parameters;
            }
        }

        /// <summary>
        /// Add a new dynamic or static modification
        /// </summary>
        /// <param name="modInfo"></param>
        public void AddModification(clsModificationDefinition modInfo)
        {
            mModInfo.Add(modInfo);
        }

        /// <summary>
        /// Add/update a parameter
        /// </summary>
        /// <param name="kvSetting"></param>
        public void AddUpdateParameter(KeyValuePair<string, string> kvSetting)
        {
            AddUpdateParameter(kvSetting.Key, kvSetting.Value);
        }

        /// <summary>
        /// Add/update a parameter
        /// </summary>
        /// <param name="paramName"></param>
        /// <param name="paramValue"></param>
        public void AddUpdateParameter(string paramName, string paramValue)
        {
            if (mParameters.ContainsKey(paramName))
            {
                mParameters[paramName] = paramValue;
            }
            else
            {
                mParameters.Add(paramName, paramValue);
            }
        }

        /// <summary>
        /// Clear stored dynamic and static modifications
        /// </summary>
        public void ClearModifications()
        {
            mModInfo.Clear();
        }

        /// <summary>
        /// Clear stored key/value parameters
        /// </summary>
        public void ClearParameters()
        {
            mParameters.Clear();
        }

        private void InitializeDefaults()
        {
            mSearchEngineName = "Unknown";
            mSearchEngineVersion = "Unknown";
            mSearchDate = new DateTime(1980, 1, 1);

            FastaFilePath = string.Empty;

            PrecursorMassToleranceDa = 0;
            PrecursorMassTolerancePpm = 0;

            mPrecursorMassType = MASS_TYPE_MONOISOTOPIC;
            mFragmentMassType = MASS_TYPE_MONOISOTOPIC;

            mEnzyme = "trypsin";
            MaxNumberInternalCleavages = 4;
            MinNumberTermini = 0;
        }

        /// <summary>
        /// Update the search engine parameter file path
        /// </summary>
        /// <param name="paramFilePath"></param>
        public void UpdateSearchEngineParamFilePath(string paramFilePath)
        {
            mSearchEngineParamFilePath = paramFilePath;
        }

        /// <summary>
        /// Update the search engine version
        /// </summary>
        /// <param name="searchEngineVersion"></param>
        public void UpdateSearchEngineVersion(string searchEngineVersion)
        {
            mSearchEngineVersion = string.Copy(searchEngineVersion);
        }

        /// <summary>
        /// Update the search date
        /// </summary>
        /// <param name="searchDate"></param>
        public void UpdateSearchDate(DateTime searchDate)
        {
            mSearchDate = searchDate;
        }
    }
}
