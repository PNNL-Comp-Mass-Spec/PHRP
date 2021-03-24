using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PHRPReader;
using PHRPReader.Data;

namespace PeptideHitResultsProcessor
{
    class clsSearchResultsMaxQuant : clsSearchResultsBaseClass
    {
        #region "Properties"

        public string DatasetName { get; set; }
        public int DatasetID { get; set; }
        public string MissedCleavageCount { get; set; }
        public string CollisionMode { get; set; }
        public string MassAnalyzer { get; set; }
        public string IdType { get; set; }
        public string ScanEvent { get; set; }
        public string IsotopeIndex { get; set; }
        public string Precursor_mz { get; set; }
        public string MaxQuantComputedDelMPPM { get; set; }
        public string MaxQuantComputedDelM { get; set; }
        public string SingleMassErrorPPM { get; set; }
        public string ScanTimeMinutes { get; set; }
        public string PEP { get; set; }
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

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        public clsSearchResultsMaxQuant(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
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
