﻿namespace PHRPReader
{
    /// <summary>
    /// PHRP Startup options
    /// </summary>
    /// <remarks>Use these options to define load behavior to be used when instantiating PHRP reader</remarks>
    public class PHRPStartupOptions
    {
        /// <summary>
        /// If true, load the modification and SeqInfo data
        /// </summary>
        public bool LoadModsAndSeqInfo { get; set; }

        /// <summary>
        /// If true, load MSGF results (not MS-GF+)
        /// </summary>
        public bool LoadMSGFResults { get; set; }

        /// <summary>
        /// If true, load ScanStats data
        /// </summary>
        public bool LoadScanStatsData { get; set; }

        /// <summary>
        /// Maximum number of proteins to associate with each PSM
        /// </summary>
        /// <remarks>Set to 0 to load all proteins</remarks>
        public int MaxProteinsPerPSM { get; set; }

        /// <summary>
        /// Use this to override the default peptide mass calculator class;
        /// this is useful if custom amino acids are in use
        /// </summary>
        public PeptideMassCalculator PeptideMassCalculator { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public PHRPStartupOptions()
        {
            LoadModsAndSeqInfo = true;
            LoadMSGFResults = true;
            LoadScanStatsData = false;
            MaxProteinsPerPSM = 0;      // 0 means to load all proteins
        }
    }
}