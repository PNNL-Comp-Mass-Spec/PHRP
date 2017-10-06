// This class is used to track the peptide details for a MODa search result
// See clsSearchResultsBaseClass for additional information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Created 4/01/2014
//
// E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsSearchResultsMODa : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables
        #endregion

        #region "Properties"

        // Auto-Properties
        public string Spectrum_Index { get; set; }

        public string Precursor_mz { get; set; }            // Observed precursor m/z, converted to monoisotopic mass by MODa

        public string MODaComputedDelM { get; set; }
        public string MODaComputedDelMPPM { get; set; }

        public string MODaScore { get; set; }

        public string Probability { get; set; }

        #endregion

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        public clsSearchResultsMODa(clsPeptideModificationContainer objPeptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(objPeptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();
            Spectrum_Index = string.Empty;

            Precursor_mz = string.Empty;

            MODaComputedDelM = string.Empty;
            MODaComputedDelMPPM = string.Empty;

            MODaScore = string.Empty;

            Probability = string.Empty;
        }
    }
}
