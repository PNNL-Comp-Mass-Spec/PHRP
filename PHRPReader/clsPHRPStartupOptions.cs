namespace PHRPReader
{
    public class clsPHRPStartupOptions
    {
        private int mMaxProteinsPerPSM;

        public bool LoadModsAndSeqInfo { get; set; }
        public bool LoadMSGFResults { get; set; }
        public bool LoadScanStatsData { get; set; }

        /// <summary>
        /// Maximum number of proteins to associate with each PSM
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>0 means to load all proteins</remarks>
        public int MaxProteinsPerPSM
        {
            get
            {
                if (mMaxProteinsPerPSM <= 0)
                {
                    return int.MaxValue;
                }
                return mMaxProteinsPerPSM;
            }
            set { mMaxProteinsPerPSM = value; }
        }

        /// <summary>
        /// Use this to override the default peptide mass calculator class;
        /// this is useful if custom amino acids are in use
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPeptideMassCalculator PeptideMassCalculator { get; set; }

        public clsPHRPStartupOptions()
        {
            LoadModsAndSeqInfo = true;
            LoadMSGFResults = true;
            LoadScanStatsData = false;
            mMaxProteinsPerPSM = 0;      // 0 means to load all proteins
        }
    }
}
