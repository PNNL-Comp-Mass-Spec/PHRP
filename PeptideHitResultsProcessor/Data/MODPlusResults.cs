﻿// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://www.pnnl.gov/integrative-omics
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
    /// This class is used to track the peptide details for a MODPlus search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class MODPlusResults : SearchResultsBaseClass
    {
        /// <summary>
        /// Spectrum index
        /// </summary>
        public string Spectrum_Index { get; set; }

        /// <summary>
        /// Observed precursor m/z
        /// </summary>
        public string Precursor_mz { get; set; }

        /// <summary>
        /// Mass error, as computed by MODPlus
        /// </summary>
        public string MODPlusComputedDelM { get; set; }

        /// <summary>
        /// Mass error (in ppm), as computed by MODPlus
        /// </summary>
        public string MODPlusComputedDelMPPM { get; set; }

        /// <summary>
        /// MODPlus score
        /// </summary>
        public string MODPlusScore { get; set; }

        /// <summary>
        /// Probability
        /// </summary>
        public string Probability { get; set; }

        /// <summary>
        /// Peptide position
        /// </summary>
        public string PeptidePosition { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public MODPlusResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        /// <summary>
        /// Clear stored values
        /// </summary>
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
