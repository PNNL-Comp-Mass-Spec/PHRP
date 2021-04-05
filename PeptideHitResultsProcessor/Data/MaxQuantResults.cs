using PHRPReader;
using PHRPReader.Data;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for a MaxQuant search result
    /// See SearchResultsBaseClass for additional information
    /// </summary>
    public class MaxQuantResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Acetyl, MaxQuant

        public string DatasetName { get; set; }
        public int DatasetID { get; set; }
        public string MissedCleavageCount { get; set; }

        /// <summary>
        /// Comma-separated list of modification names that are in the ModifiedSequence field
        /// </summary>
        /// <remarks>
        /// Examples:
        ///   Acetyl (Protein N-term)
        ///   Oxidation (M)
        ///   Acetyl (Protein N-term),Oxidation (M)
        ///   2 Oxidation (M)
        ///   4 Oxidation (M)
        /// </remarks>
        public string Modifications { get; set; }

        public string ModifiedSequence { get; set; }
        public string CollisionMode { get; set; }

        /// <summary>
        /// Mass analyzer
        /// </summary>
        /// <remarks>
        /// FTMS, ITMS
        /// </remarks>
        public string MassAnalyzer { get; set; }

        /// <summary>
        /// Identification type
        /// </summary>
        /// <remarks>
        /// MULTI-MSMS, MULTI-SECPEP, or MSMS
        /// </remarks>
        public string IdType { get; set; }
        public string ScanEvent { get; set; }
        public string IsotopeIndex { get; set; }
        public string Precursor_mz { get; set; }
        public string MaxQuantComputedDelMPPM { get; set; }
        public string MaxQuantComputedDelM { get; set; }
        public string SingleMassErrorPPM { get; set; }
        public string ScanTimeMinutes { get; set; }

        /// <summary>
        /// Posterior error probability
        /// </summary>
        /// <remarks>
        /// Lower is better
        /// </remarks>
        public string PEP { get; set; }

        /// <summary>
        /// Confidence score; higher is better
        /// </summary>
        public string Score { get; set; }
        public string DeltaScore { get; set; }
        public string ScoreDiff { get; set; }
        public string LocalizationProb { get; set; }
        public string Combinatorics { get; set; }
        public string PIF { get; set; }
        public string PrecursorScan { get; set; }
        public string PrecursorIntensity { get; set; }
        public string MatchedFragmentIonCount { get; set; }
        public string IntensityCoverage { get; set; }
        public string PeakCoverage { get; set; }
        public string IsReversedProtein { get; set; }

        /// <summary>
        /// Semicolon separated list of protein group IDs
        /// </summary>
        public string ProteinGroupIDs { get; set; }

        public string PeptideID { get; set; }
        public string ModPeptideID { get; set; }
        public string EvidenceID { get; set; }
        public string OxidationSiteIDs { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        public MaxQuantResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            DatasetName = string.Empty;
            DatasetID = 0;
            MissedCleavageCount = string.Empty;
            CollisionMode = string.Empty;
            MassAnalyzer = string.Empty;
            IdType = string.Empty;
            ScanEvent = string.Empty;
            IsotopeIndex = string.Empty;
            Precursor_mz = string.Empty;
            MaxQuantComputedDelMPPM = string.Empty;
            MaxQuantComputedDelM = string.Empty;
            SingleMassErrorPPM = string.Empty;
            ScanTimeMinutes = string.Empty;
            PEP = string.Empty;
            Score = string.Empty;
            DeltaScore = string.Empty;
            ScoreDiff = string.Empty;
            LocalizationProb = string.Empty;
            Combinatorics = string.Empty;
            PIF = string.Empty;
            PrecursorScan = string.Empty;
            PrecursorIntensity = string.Empty;
            MatchedFragmentIonCount = string.Empty;
            IntensityCoverage = string.Empty;
            PeakCoverage = string.Empty;
            IsReversedProtein = string.Empty;
            ProteinGroupIDs = string.Empty;
            PeptideID = string.Empty;
            ModPeptideID = string.Empty;
            EvidenceID = string.Empty;
            OxidationSiteIDs = string.Empty;
        }
    }
}
