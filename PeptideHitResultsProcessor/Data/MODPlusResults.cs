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
    /// This class is used to track the peptide details for a MODPlus search result
    /// See SearchResultsBaseClass for additional information
    /// </summary>
    public class MODPlusResults : SearchResultsBaseClass
    {
        public string Spectrum_Index { get; set; }

        /// <summary>
        /// Observed precursor m/z
        /// </summary>
        public string Precursor_mz { get; set; }

        public string MODPlusComputedDelM { get; set; }

        public string MODPlusComputedDelMPPM { get; set; }

        public string MODPlusScore { get; set; }

        public string Probability { get; set; }

        public string PeptidePosition { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        public MODPlusResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
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
