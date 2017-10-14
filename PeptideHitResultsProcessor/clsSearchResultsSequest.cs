// This class is used to track the peptide details for a Sequest search result
// See clsSearchResultsBaseClass for additional information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 7, 2006
//
// E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
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
using System;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsSearchResultsSequest : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"

        protected string mPeptideXCorrNext;             // Not stored in the synopsis or first hits file; it can be computed from XCorr and DeltaCn2
        #endregion

        #region "Properties"
        // Auto-properties (only work with Visual Studio 2010)
        public string NumScans { get; set; }
        public string PeptideDeltaCn { get; set; }
        public string PeptideDeltaCn2 { get; set; }
        public string PeptideMScore { get; set; }
        public string PeptideNTT { get; set; }
        public string PeptidePassFilt { get; set; }
        public string PeptideRankSP { get; set; }
        public string PeptideRankXC { get; set; }
        public string PeptideSp { get; set; }
        public string PeptideXCorr { get; set; }
        public string PeptideXcRatio { get; set; }
        public string IonsObserved { get; set; }
        public string IonsExpected { get; set; }
        public string DelMPPM { get; set; }

        public string PeptideXCorrNext => mPeptideXCorrNext;

        #endregion

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        public clsSearchResultsSequest(clsPeptideModificationContainer objPeptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(objPeptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            NumScans = string.Empty;
            PeptideDeltaCn = string.Empty;
            PeptideDeltaCn2 = string.Empty;
            PeptideMScore = string.Empty;
            PeptideNTT = string.Empty;
            PeptidePassFilt = string.Empty;
            PeptideRankSP = string.Empty;
            PeptideRankXC = string.Empty;
            PeptideSp = string.Empty;
            PeptideXCorr = string.Empty;
            PeptideXcRatio = string.Empty;
            IonsObserved = string.Empty;
            IonsExpected = string.Empty;
            DelMPPM = string.Empty;

            mPeptideXCorrNext = string.Empty;
        }

        protected void ComputePeptideXCorrNext()
        {
            float sngXCorr = 0;
            float sngDelCN2 = 0;

            try
            {
                if (float.TryParse(PeptideXCorr, out sngXCorr) &&
                    float.TryParse(PeptideDeltaCn2, out sngDelCN2))
                {
                    mPeptideXCorrNext = (sngXCorr - sngDelCN2 * sngXCorr).ToString();
                }
                else
                {
                    mPeptideXCorrNext = "0";
                }
            }
            catch (Exception)
            {
                mPeptideXCorrNext = "0";
            }
        }
    }
}
