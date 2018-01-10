// This class is used to track the peptide details for a MODPlus search result
// See clsSearchResultsBaseClass for additional information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Created 5/15/2015
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsSearchResultsMODPlus : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables
        #endregion

        #region "Properties"

        // Auto-Properties
        public string Spectrum_Index { get; set; }

        public string Precursor_mz { get; set; }            // Observed precursor m/z, converted to monoisotopic mass by MODa

        public string MODPlusComputedDelM { get; set; }
        public string MODPlusComputedDelMPPM { get; set; }

        public string MODPlusScore { get; set; }

        public string Probability { get; set; }

        public string PeptidePosition { get; set; }

        public string ModificationAnnotation { get; set; }

        #endregion

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        public clsSearchResultsMODPlus(clsPeptideModificationContainer objPeptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(objPeptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();
            Spectrum_Index = string.Empty;

            Precursor_mz = string.Empty;

            MODPlusComputedDelM = string.Empty;
            MODPlusComputedDelMPPM = string.Empty;

            MODPlusScore = string.Empty;

            Probability = string.Empty;

            PeptidePosition = string.Empty;
        }
    }
}
