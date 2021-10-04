// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 7, 2006
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

using System;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for an XTandem search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class XTandemResults : SearchResultsBaseClass
    {
        // Ignore Spelling: delM

        // Note: ProteinExpectationValue and ProteinIntensity are defined in SearchResultsBaseClass
        // The raw expectation value from the results file is converted to the Base-10 Log form when read into this program
        private string mPeptideNextScore;

        /// <summary>
        /// A multiplier to convert the normalized spectrum contained in a group back to the original intensity values
        /// </summary>
        public string IntensityMultiplier { get; set; }

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

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        public XTandemResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            IntensityMultiplier = string.Empty;

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
            double precursorMonoMass = 0;

            var parseError = false;

            // Note that mPeptideDeltaMass is the DeltaMass value reported by X!Tandem
            // (though XtandemResultsProcessor took the negative of the value in the results file so it currently represents "theoretical - observed")
            if (double.TryParse(PeptideDeltaMass, out var delM))
            {
                // Negate delM so that it represents observed - theoretical
                delM = -delM;

                // Compute the observed precursor monoisotopic mass
                if (double.TryParse(ParentIonMH, out var parentIonMH))
                {
                    precursorMonoMass = parentIonMH - PeptideMassCalculator.MASS_PROTON;
                }
                else
                {
                    parseError = true;
                }

                if (parseError)
                {
                    precursorMonoMass = PeptideMonoisotopicMass + delM;
                }

                const bool adjustPrecursorMassForC13 = true;
                PeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(delM, precursorMonoMass, PeptideMonoisotopicMass, adjustPrecursorMassForC13);
            }
            else
            {
                PeptideDeltaMassCorrectedPpm = 0;
            }
        }

        private void ComputePeptideDeltaCn2()
        {
            try
            {
                if (SynFileReaderBaseClass.IsNumber(PeptideHyperscore) && SynFileReaderBaseClass.IsNumber(mPeptideNextScore))
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
