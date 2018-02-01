// This class is used to track the peptide details for an InSpecT search result
// See clsSearchResultsBaseClass for additional information
//
// -------------------------------------------------------------------------------
// Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
// Program started August 19, 2008
//
// E-mail: john.sandoval@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsSearchResultsInSpecT : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"

        // Scan: tracked by the base class
        // Annotation: aka Peptide, which is tracked by the base class
        // Protein: tracked by the base class

        #endregion

        #region "Properties"

        // Auto-Properties
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

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
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
