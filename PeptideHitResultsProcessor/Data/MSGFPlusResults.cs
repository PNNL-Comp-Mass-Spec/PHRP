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
    /// This class is used to track the peptide details for an MS-GF+ search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class MSGFPlusResults : SearchResultsBaseClass
    {
        // Ignore Spelling: De novo, tda, tryptic

        /// <summary>
        /// Fragmentation method for the given MS/MS spectrum
        /// </summary>
        public string FragMethod { get; set; }

        /// <summary>
        /// Spectrum index
        /// </summary>
        public string SpecIndex { get; set; }

        /// <summary>
        /// Precursor ion m/z (observed value), as reported by MS-GF+
        /// </summary>
        public string PrecursorMZ { get; set; }

        /// <summary>
        /// Mass error value computed by MS-GF+
        /// </summary>
        public string MSGFPlusComputedDelM { get; set; }

        /// <summary>
        /// Ppm-based mass error value computed by MS-GF+
        /// </summary>
        public string MSGFPlusComputedDelMPPM { get; set; }

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        public string NTT { get; set; }

        /// <summary>
        /// De novo score
        /// </summary>
        public string DeNovoScore { get; set; }

        /// <summary>
        /// MSGF Score
        /// </summary>
        public string MSGFScore { get; set; }

        /// <summary>
        /// SpecEValue in MS-GF+
        /// SpecProb in MSGFDB
        /// </summary>
        public string SpecEValue { get; set; }

        /// <summary>
        /// Rank_MSGFDB_SpecEValue in MS-GF+
        /// Rank_MSGFDB_SpecProb in MSGFDB
        /// </summary>
        public string RankSpecEValue { get; set; }

        /// <summary>
        /// EValue in MS-GF+
        /// PValue in MSGFDB
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
        /// <remarks>
        /// Only present if searched using -tda 1
        /// PepQValue in MS-GF+
        /// PepFDR in MSGFDB
        /// </remarks>
        public string PepQValue { get; set; }

        /// <summary>
        /// Isotope error (integer value)
        /// </summary>
        /// <remarks>Only reported by MS-GF+</remarks>
        public string IsotopeError { get; set; }

        /// <summary>
        /// True if the loaded data is from MS-GF+
        /// False if from MSGFDB
        /// </summary>
        public bool UsedMSGFPlus { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public MSGFPlusResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        /// <summary>
        /// Clear stored values
        /// </summary>
        public override void Clear()
        {
            base.Clear();

            FragMethod = string.Empty;
            SpecIndex = string.Empty;

            PrecursorMZ = string.Empty;
            MSGFPlusComputedDelM = string.Empty;
            MSGFPlusComputedDelMPPM = string.Empty;

            NTT = string.Empty;

            DeNovoScore = string.Empty;
            MSGFScore = string.Empty;

            SpecEValue = string.Empty;
            RankSpecEValue = string.Empty;

            EValue = string.Empty;
            QValue = string.Empty;
            PepQValue = string.Empty;

            IsotopeError = string.Empty;
            UsedMSGFPlus = false;
        }
    }
}
