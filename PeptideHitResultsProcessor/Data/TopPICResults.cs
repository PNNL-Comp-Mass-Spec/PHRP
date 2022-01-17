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
    /// This class is used to track the peptide details for a TopPIC search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    internal class TopPICResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Proteoform

        /// <summary>
        ///  Result ID
        /// </summary>
        public string Prsm_ID { get; set; }

        /// <summary>
        /// Spectrum ID
        /// </summary>
        public string Spectrum_ID { get; set; }

        /// <summary>
        /// Unexpected modification count
        /// </summary>
        public string Unexpected_Mod_Count { get; set; }

        /// <summary>
        /// Number of peaks (m/z values) in the MS/MS spectrum
        /// </summary>
        public string Peak_Count { get; set; }

        /// <summary>
        /// Matched peak count
        /// </summary>
        public string Matched_Peak_Count { get; set; }

        /// <summary>
        /// Matched fragment ion count
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
        /// Q-value
        /// </summary>
        public string QValue { get; set; }

        /// <summary>
        /// Fragmentation method
        /// </summary>
        public string FragMethod { get; set; }

        /// <summary>
        /// Proteoform level FDR
        /// </summary>
        public string ProteoformFDR { get; set; }

        /// <summary>
        /// Dynamic modifications
        /// </summary>
        public string VariablePTMs { get; set; }

        /// <summary>
        /// Observed precursor_mz
        /// </summary>
        public string Precursor_mz { get; set; }

        /// <summary>
        /// Mass error, as computed by TopPIC
        /// </summary>
        public string TopPICComputedDelM { get; set; }

        /// <summary>
        /// Mass error (in ppm), as computed by TopPIC
        /// </summary>
        public string TopPICComputedDelMPPM { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        public TopPICResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
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
