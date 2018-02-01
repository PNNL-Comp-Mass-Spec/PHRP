// This class is used to track the peptide details for a MSAlign search result
// See clsSearchResultsBaseClass for additional information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Created 11/27/2012
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsSearchResultsMSAlign : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables
        #endregion

        #region "Properties"

        // Auto-Properties
        public string Prsm_ID { get; set; }
        public string Spectrum_ID { get; set; }

        public string Protein_Mass { get; set; }
        public string Unexpected_Mod_Count { get; set; }

        public string Peak_Count { get; set; }
        public string Matched_Peak_Count { get; set; }
        public string Matched_Fragment_Ion_Count { get; set; }

        public string PValue { get; set; }
        public string Rank_PValue { get; set; }

        public string EValue { get; set; }
        public string FDR { get; set; }

        public string Species_ID { get; set; }
        public string FragMethod { get; set; }

        public string Precursor_mz { get; set; }            // Observed precursor_mz
        public string MSAlignComputedDelM { get; set; }
        public string MSAlignComputedDelMPPM { get; set; }

        #endregion

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        public clsSearchResultsMSAlign(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            Prsm_ID = string.Empty;
            Spectrum_ID = string.Empty;

            Protein_Mass = string.Empty;
            Unexpected_Mod_Count = string.Empty;

            Peak_Count = string.Empty;
            Matched_Peak_Count = string.Empty;
            Matched_Fragment_Ion_Count = string.Empty;

            PValue = string.Empty;
            Rank_PValue = string.Empty;

            EValue = string.Empty;
            FDR = string.Empty;

            Species_ID = string.Empty;
            FragMethod = string.Empty;

            Precursor_mz = string.Empty;
            MSAlignComputedDelM = string.Empty;
            MSAlignComputedDelMPPM = string.Empty;
        }
    }
}
