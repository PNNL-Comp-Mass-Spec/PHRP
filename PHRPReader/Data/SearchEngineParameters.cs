using System;
using System.Collections.Generic;

namespace PHRPReader.Data
{
    /// <summary>
    /// Search engine parameters container
    /// </summary>
    public class SearchEngineParameters
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
        public List<ModificationDefinition> ModList { get; }

        /// <summary>
        /// Parameter dictionary (key/value pairs)
        /// </summary>
        public Dictionary<string, string> Parameters { get; }

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

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="searchEngineName"></param>
        public SearchEngineParameters(string searchEngineName) : this(searchEngineName, new List<ModificationDefinition>(), null)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="searchEngineName"></param>
        /// <param name="modList"></param>
        public SearchEngineParameters(string searchEngineName, List<ModificationDefinition> modList) : this(searchEngineName, modList, null)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="searchEngineName"></param>
        /// <param name="modList"></param>
        /// <param name="parameters"></param>
        public SearchEngineParameters(string searchEngineName, List<ModificationDefinition> modList, Dictionary<string, string> parameters)
        {
            InitializeDefaults();

            mSearchEngineName = searchEngineName;
            mSearchEngineParamFilePath = string.Empty;

            ModList = modList;

            ModList = modList ?? new List<ModificationDefinition>();

            if (Parameters == null)
            {
                Parameters = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            }
            else
            {
                Parameters = parameters;
            }
        }

        /// <summary>
        /// Add a new dynamic or static modification
        /// </summary>
        /// <param name="modDef"></param>
        public void AddModification(ModificationDefinition modDef)
        {
            ModList.Add(modDef);
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
            // Add/update the dictionary
            Parameters[paramName] = paramValue;
        }

        /// <summary>
        /// Clear stored dynamic and static modifications
        /// </summary>
        public void ClearModifications()
        {
            ModList.Clear();
        }

        /// <summary>
        /// Clear stored key/value parameters
        /// </summary>
        public void ClearParameters()
        {
            Parameters.Clear();
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
