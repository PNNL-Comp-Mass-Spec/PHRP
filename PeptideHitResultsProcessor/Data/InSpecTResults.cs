// -------------------------------------------------------------------------------
// Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
// Program started August 19, 2008
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
    /// This class is used to track the peptide details for an InSpecT search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class InSpecTResults : SearchResultsBaseClass
    {
        // Ignore Spelling: NTT, tryptic

        // Scan: tracked by the base class
        // Annotation: aka Peptide, which is tracked by the base class
        // Protein: tracked by the base class

        /// <summary>
        /// Spectrum file name
        /// </summary>
        public string SpectrumFile { get; set; }

        /// <summary>
        /// MQScore
        /// </summary>
        public string MQScore { get; set; }

        /// <summary>
        /// Length
        /// </summary>
        public string Length { get; set; }

        /// <summary>
        /// TotalPRMScore
        /// </summary>
        public string TotalPRMScore { get; set; }

        /// <summary>
        /// MedianPRMScore
        /// </summary>
        public string MedianPRMScore { get; set; }

        /// <summary>
        /// FractionY
        /// </summary>
        public string FractionY { get; set; }

        /// <summary>
        /// FractionB
        /// </summary>
        public string FractionB { get; set; }

        /// <summary>
        /// Intensity
        /// </summary>
        public string Intensity { get; set; }

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        public string NTT { get; set; }

        /// <summary>
        /// PValue
        /// </summary>
        public string PValue { get; set; }

        /// <summary>
        /// FScore
        /// </summary>
        public string FScore { get; set; }

        /// <summary>
        /// DeltaScore
        /// </summary>
        public string DeltaScore { get; set; }

        /// <summary>
        /// DeltaScoreOther
        /// </summary>
        public string DeltaScoreOther { get; set; }

        /// <summary>
        /// DeltaNormMQScore
        /// </summary>
        public string DeltaNormMQScore { get; set; }

        /// <summary>
        /// DeltaNormTotalPRMScore
        /// </summary>
        public string DeltaNormTotalPRMScore { get; set; }

        /// <summary>
        /// RankTotalPRMScore
        /// </summary>
        public string RankTotalPRMScore { get; set; }

        /// <summary>
        /// RankFScore
        /// </summary>
        public string RankFScore { get; set; }

        /// <summary>
        /// RecordNumber
        /// </summary>
        public string RecordNumber { get; set; }

        /// <summary>
        /// DBFilePos
        /// </summary>
        public string DBFilePos { get; set; }

        /// <summary>
        /// SpecFilePos
        /// </summary>
        public string SpecFilePos { get; set; }

        /// <summary>
        /// PrecursorMz
        /// </summary>
        public string PrecursorMz { get; set; }

        /// <summary>
        /// PrecursorError
        /// </summary>
        public string PrecursorError { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public InSpecTResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        /// <summary>
        /// Clear stored values
        /// </summary>
        public override void Clear()
        {
            base.Clear();

            SpectrumFile = string.Empty;
            MQScore = string.Empty;
            Length = string.Empty;
            TotalPRMScore = string.Empty;
            MedianPRMScore = string.Empty;
            FractionY = string.Empty;
            FractionB = string.Empty;
            Intensity = string.Empty;
            NTT = string.Empty;
            PValue = string.Empty;
            FScore = string.Empty;
            DeltaScore = string.Empty;
            DeltaScoreOther = string.Empty;
            DeltaNormMQScore = string.Empty;
            DeltaNormTotalPRMScore = string.Empty;
            RankTotalPRMScore = string.Empty;
            RankFScore = string.Empty;
            RecordNumber = string.Empty;
            DBFilePos = string.Empty;
            SpecFilePos = string.Empty;
            PrecursorMz = string.Empty;
            PrecursorError = string.Empty;
        }
    }
}
