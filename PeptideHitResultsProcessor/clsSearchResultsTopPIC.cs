﻿// -------------------------------------------------------------------------------
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
    /// This class is used to track the peptide details for a TopPIC search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    class clsSearchResultsTopPIC : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables
        #endregion

        #region "Properties"

        // Auto-Properties
        public string Prsm_ID { get; set; }
        public string Spectrum_ID { get; set; }

        public string Unexpected_Mod_Count { get; set; }

        public string Peak_Count { get; set; }
        public string Matched_Peak_Count { get; set; }
        public string Matched_Fragment_Ion_Count { get; set; }

        public string PValue { get; set; }
        public string Rank_PValue { get; set; }

        public string EValue { get; set; }
        public string QValue { get; set; }

        public string FragMethod { get; set; }

        public string ProteoformFDR { get; set; }
        public string VariablePTMs { get; set; }

        public string Precursor_mz { get; set; }            // Observed precursor_mz
        public string TopPICComputedDelM { get; set; }
        public string TopPICComputedDelMPPM { get; set; }

        #endregion

        // Note that the following call will call both the base class's Clear method and this class's Clear method
        public clsSearchResultsTopPIC(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            Prsm_ID = string.Empty;
            Spectrum_ID = string.Empty;

            Unexpected_Mod_Count = string.Empty;

            Peak_Count = string.Empty;
            Matched_Peak_Count = string.Empty;
            Matched_Fragment_Ion_Count = string.Empty;

            PValue = string.Empty;
            Rank_PValue = string.Empty;

            EValue = string.Empty;
            QValue = string.Empty;

            FragMethod = string.Empty;
            ProteoformFDR = string.Empty;
            VariablePTMs = string.Empty;

            Precursor_mz = string.Empty;
            TopPICComputedDelM = string.Empty;
            TopPICComputedDelMPPM = string.Empty;
        }
    }
}