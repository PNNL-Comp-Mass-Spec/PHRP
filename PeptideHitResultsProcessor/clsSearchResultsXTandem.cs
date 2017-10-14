// This class is used to track the peptide details for an XTandem search result
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
    public class clsSearchResultsXTandem : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note: ProteinExpectationValue and ProteinIntensity are defined in clsSearchResultsBaseClass
        // The raw expectation value from the results file is converted to the Base-10 Log form when read into this program
        protected string mPeptideNextScore;

        #endregion

        #region "Properties"

        public string fI { get; set; }

        // The raw expectation value from the results file is converted to the Base-10 Log form when read into this program
        public string PeptideExpectationValue { get; set; }

        public string PeptideHyperscore { get; set; }

        public string PeptideNextScore
        {
            get => mPeptideNextScore;
            set
            {
                mPeptideNextScore = value;
                ComputePeptideDeltaCn2();
            }
        }

        public float PeptideDeltaCn2 { get; set; }

        public string PeptideYScore { get; set; }

        public string PeptideYIons { get; set; }

        public string PeptideBScore { get; set; }

        public string PeptideBIons { get; set; }

        public string PeptideIntensity { get; set; }

        public string PeptideIntensityMax { get; set; }

        public double PeptideDeltaMassCorrectedPpm { get; set; }

        #endregion

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        public clsSearchResultsXTandem(clsPeptideModificationContainer objPeptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(objPeptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            fI = string.Empty;

            PeptideExpectationValue = string.Empty;
            PeptideHyperscore = string.Empty;
            mPeptideNextScore = string.Empty;
            PeptideDeltaCn2 = 0;

            PeptideYScore = string.Empty;
            PeptideYIons = string.Empty;
            PeptideBScore = string.Empty;
            PeptideBIons = string.Empty;

            PeptideIntensity = string.Empty;
            PeptideIntensityMax = string.Empty;

            PeptideDeltaMassCorrectedPpm = 0;
        }

        public void ComputeDelMCorrectedXT()
        {
            double dblDelM = 0;

            double dblPrecursorMonoMass = 0;

            var blnParseError = false;

            // Note that mPeptideDeltaMass is the DeltaMass value reported by X!Tandem
            // (though clsXtandemResultsProcessor took the negative of the value in the results file so it currently represents "theoretical - observed")
            if (double.TryParse(PeptideDeltaMass, out dblDelM))
            {
                // Negate dblDelM so that it represents observed - theoretical
                dblDelM = -dblDelM;

                // Compute the original value for the precursor monoisotopic mass
                double dblParentIonMH = 0;
                if (double.TryParse(base.ParentIonMH, out dblParentIonMH))
                {
                    dblPrecursorMonoMass = dblParentIonMH - PHRPReader.clsPeptideMassCalculator.MASS_PROTON;
                }
                else
                {
                    blnParseError = true;
                }

                if (blnParseError)
                {
                    dblPrecursorMonoMass = PeptideMonoisotopicMass + dblDelM;
                }

                const bool blnAdjustPrecursorMassForC13 = true;
                PeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, PeptideMonoisotopicMass);
            }
            else
            {
                PeptideDeltaMassCorrectedPpm = 0;
            }
        }

        protected void ComputePeptideDeltaCn2()
        {
            try
            {
                if (clsPHRPParser.IsNumber(PeptideHyperscore) & clsPHRPParser.IsNumber(mPeptideNextScore))
                {
                    PeptideDeltaCn2 = (Convert.ToSingle(PeptideHyperscore) - Convert.ToSingle(mPeptideNextScore)) / Convert.ToSingle(PeptideHyperscore);
                }
                else
                {
                    PeptideDeltaCn2 = 0;
                }
            }
            catch (Exception)
            {
                PeptideDeltaCn2 = 0;
            }
        }
    }
}
