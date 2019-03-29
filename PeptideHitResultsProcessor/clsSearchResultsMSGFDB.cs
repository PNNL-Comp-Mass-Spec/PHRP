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
    /// This class is used to track the peptide details for an MSGF+ search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    public class clsSearchResultsMSGFDB : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables
        #endregion

        #region "Properties"

        public string FragMethod { get; set; }
        public string NTT { get; set; }

        public string DeNovoScore { get; set; }
        public string MSGFScore { get; set; }
        public string SpecEValue { get; set; }      // SpecProb in MSGFDB; SpecEValue in MSGF+
        public string RankSpecEValue { get; set; }
        public string EValue { get; set; }          // PValue in MSGFDB; EValue in MSGF+

        // ReSharper disable once CommentTypo
        public string QValue { get; set; }          // Will contain target/decoy FDR when -tda 1 was used; will contain EFDR when -tda 1 was not used; FDR in MSGFDB; QValue in MSGF+
        public string PepQValue { get; set; }       // Only present if searched using -tda 1; PepFDR in MSGFDB; PepQValue in MSGF+

        public string PrecursorMZ { get; set; }
        public string MSGFPlusComputedDelM { get; set; }
        public string MSGFPlusComputedDelMPPM { get; set; }

        public string IsotopeError { get; set; }    // Only reported by MSGF+

        public bool MSGFPlusResults { get; set; }
        #endregion

        // Note that the following call will call both the base class's Clear method and this class's Clear method
        public clsSearchResultsMSGFDB(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
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
