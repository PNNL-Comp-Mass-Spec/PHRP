// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or https://www.pnnl.gov/sysbio/ or https://panomics.pnnl.gov/
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class is used to track the peptide details for a MODPlus search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
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

        // Note that the following call will call both the base class's Clear method and this class's Clear method
        public clsSearchResultsMODPlus(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
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
