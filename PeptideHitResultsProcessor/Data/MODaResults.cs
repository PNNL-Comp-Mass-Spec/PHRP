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
using PHRPReader.Data;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class is used to track the peptide details for a MODa search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    public class clsSearchResultsMODa : clsSearchResultsBaseClass
    {
        // Ignore Spelling: MODa

        #region "Properties"

        public string Spectrum_Index { get; set; }

        /// <summary>
        /// Observed precursor m/z, converted to monoisotopic mass by MODa
        /// </summary>
        public string Precursor_mz { get; set; }

        public string MODaComputedDelM { get; set; }
        public string MODaComputedDelMPPM { get; set; }

        public string MODaScore { get; set; }

        public string Probability { get; set; }

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
        public clsSearchResultsMODa(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
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