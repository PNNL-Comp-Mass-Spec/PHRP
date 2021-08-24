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
using PHRPReader.Data;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for a MSAlign search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class MSAlignResults : SearchResultsBaseClass
    {
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

        /// <summary>
        /// Observed precursor_mz
        /// </summary>
        public string Precursor_mz { get; set; }
        public string MSAlignComputedDelM { get; set; }
        public string MSAlignComputedDelMPPM { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        public MSAlignResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
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
