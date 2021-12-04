// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
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
        /// <summary>
        /// Result ID
        /// </summary>
        public string Prsm_ID { get; set; }

        /// <summary>
        /// Spectrum ID
        /// </summary>
        public string Spectrum_ID { get; set; }

        /// <summary>
        /// Protein Mass
        /// </summary>
        public string Protein_Mass { get; set; }

        /// <summary>
        /// Unexpected Mod Count
        /// </summary>
        public string Unexpected_Mod_Count { get; set; }

        /// <summary>
        /// Peak Count
        /// </summary>
        public string Peak_Count { get; set; }

        /// <summary>
        /// Matched Peak Count
        /// </summary>
        public string Matched_Peak_Count { get; set; }

        /// <summary>
        /// Matched Fragment Ion Count
        /// </summary>
        public string Matched_Fragment_Ion_Count { get; set; }

        /// <summary>
        /// p-value
        /// </summary>
        public string PValue { get; set; }

        /// <summary>
        /// Rank p-value
        /// </summary>
        public string Rank_PValue { get; set; }

        /// <summary>
        /// E-value
        /// </summary>
        public string EValue { get; set; }

        /// <summary>
        /// False discovery rate
        /// </summary>
        public string FDR { get; set; }

        /// <summary>
        /// Species ID
        /// </summary>
        public string Species_ID { get; set; }

        /// <summary>
        /// Fragmentation Method
        /// </summary>
        public string FragMethod { get; set; }

        /// <summary>
        /// Observed precursor m/z
        /// </summary>
        public string Precursor_mz { get; set; }

        /// <summary>
        /// Mass error, as computed by MSAlign
        /// </summary>
        public string MSAlignComputedDelM { get; set; }

        /// <summary>
        /// Mass error (in ppm), as computed by MSAlign
        /// </summary>
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

        /// <summary>
        /// Clear stored values
        /// </summary>
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
