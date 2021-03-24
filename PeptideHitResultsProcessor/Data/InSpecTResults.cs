// -------------------------------------------------------------------------------
// Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
// Program started August 19, 2008
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
    /// This class is used to track the peptide details for an InSpecT search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    public class clsSearchResultsInSpecT : clsSearchResultsBaseClass
    {
        #region "Class wide Variables"

        // Scan: tracked by the base class
        // Annotation: aka Peptide, which is tracked by the base class
        // Protein: tracked by the base class

        #endregion

        #region "Properties"

        public string SpectrumFile { get; set; }
        public string MQScore { get; set; }
        public string Length { get; set; }
        public string TotalPRMScore { get; set; }
        public string MedianPRMScore { get; set; }
        public string FractionY { get; set; }
        public string FractionB { get; set; }
        public string Intensity { get; set; }
        public string NTT { get; set; }
        public string pValue { get; set; }
        public string FScore { get; set; }
        public string DeltaScore { get; set; }
        public string DeltaScoreOther { get; set; }
        public string DeltaNormMQScore { get; set; }
        public string DeltaNormTotalPRMScore { get; set; }
        public string RankTotalPRMScore { get; set; }
        public string RankFScore { get; set; }
        public string RecordNumber { get; set; }
        public string DBFilePos { get; set; }
        public string SpecFilePos { get; set; }
        public string PrecursorMz { get; set; }
        public string PrecursorError { get; set; }

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
        public clsSearchResultsInSpecT(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

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
            pValue = string.Empty;
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
