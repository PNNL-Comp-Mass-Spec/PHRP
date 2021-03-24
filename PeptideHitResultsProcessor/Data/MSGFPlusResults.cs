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
    /// This class is used to track the peptide details for an MS-GF+ search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    public class clsSearchResultsMSGFPlus : clsSearchResultsBaseClass
    {
        // Ignore Spelling: tda

        #region "Properties"

        public string FragMethod { get; set; }
        public string NTT { get; set; }

        public string DeNovoScore { get; set; }
        public string MSGFScore { get; set; }

        /// <summary>
        /// SpecProb in MSGFDB; SpecEValue in MS-GF+
        /// </summary>
        public string SpecEValue { get; set; }
        public string RankSpecEValue { get; set; }

        /// <summary>
        /// PValue in MSGFDB; EValue in MS-GF+
        /// </summary>
        public string EValue { get; set; }

        /// <summary>
        /// QValue
        /// </summary>
        /// <remarks>
        /// Will contain target/decoy FDR when -tda 1 was used; will contain EFDR when -tda 1 was not used; FDR in MSGFDB; QValue in MS-GF+
        /// </remarks>
        public string QValue { get; set; }

        /// <summary>
        /// Pep QValue
        /// </summary>
        /// <remarks>Only present if searched using -tda 1; PepFDR in MSGFDB; PepQValue in MS-GF+</remarks>
        public string PepQValue { get; set; }

        public string PrecursorMZ { get; set; }
        public string MSGFPlusComputedDelM { get; set; }
        public string MSGFPlusComputedDelMPPM { get; set; }

        /// <summary>
        /// Isotope error
        /// </summary>
        /// <remarks>Only reported by MS-GF+</remarks>
        public string IsotopeError { get; set; }

        public bool MSGFPlusResults { get; set; }

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
        public clsSearchResultsMSGFPlus(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            FragMethod = string.Empty;
            NTT = string.Empty;

            DeNovoScore = string.Empty;
            MSGFScore = string.Empty;
            SpecEValue = string.Empty;
            RankSpecEValue = string.Empty;
            EValue = string.Empty;

            QValue = string.Empty;
            PepQValue = string.Empty;

            PrecursorMZ = string.Empty;
            MSGFPlusComputedDelM = string.Empty;
            MSGFPlusComputedDelMPPM = string.Empty;

            IsotopeError = string.Empty;
            MSGFPlusResults = false;
        }
    }
}
