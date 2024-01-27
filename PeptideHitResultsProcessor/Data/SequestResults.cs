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
using System.Globalization;
using PHRPReader;
using PHRPReader.Data;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for a SEQUEST search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class SequestResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Cn, Filt, Sequest, tryptic, Xc

        /// <summary>
        /// NumScans
        /// </summary>
        public string NumScans { get; set; }

        /// <summary>
        /// DeltaCn (score difference)
        /// </summary>
        public string PeptideDeltaCn { get; set; }

        /// <summary>
        /// DeltaCn2
        /// </summary>
        public string PeptideDeltaCn2 { get; set; }

        /// <summary>
        /// MScore
        /// </summary>
        public string PeptideMScore { get; set; }

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        public string PeptideNTT { get; set; }

        /// <summary>
        /// PassFilt score
        /// </summary>
        public string PeptidePassFilt { get; set; }

        /// <summary>
        /// Rank SP
        /// </summary>
        public string PeptideRankSP { get; set; }

        /// <summary>
        /// Rank XC
        /// </summary>
        public string PeptideRankXC { get; set; }

        /// <summary>
        /// Sp
        /// </summary>
        public string PeptideSp { get; set; }

        /// <summary>
        /// XCorr
        /// </summary>
        public string PeptideXCorr { get; set; }

        /// <summary>
        /// XcRatio
        /// </summary>
        public string PeptideXcRatio { get; set; }

        /// <summary>
        /// Ions Observed
        /// </summary>
        public string IonsObserved { get; set; }

        /// <summary>
        /// Ions Expected
        /// </summary>
        public string IonsExpected { get; set; }

        /// <summary>
        /// DelM PPM
        /// </summary>
        public string DelMPPM { get; set; }

        /// <summary>
        /// Not stored in the synopsis or first hits file; it can be computed from XCorr and DeltaCn2
        /// </summary>
        public string PeptideXCorrNext { get; private set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public SequestResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        /// <summary>
        /// Clear stored values
        /// </summary>
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

            PeptideXCorrNext = string.Empty;
        }

        // ReSharper disable once UnusedMember.Local
        private void ComputePeptideXCorrNext()
        {
            try
            {
                if (float.TryParse(PeptideXCorr, out var xCorr) &&
                    float.TryParse(PeptideDeltaCn2, out var delCN2))
                {
                    PeptideXCorrNext = (xCorr - delCN2 * xCorr).ToString(CultureInfo.InvariantCulture);
                }
                else
                {
                    PeptideXCorrNext = "0";
                }
            }
            catch (Exception)
            {
                PeptideXCorrNext = "0";
            }
        }
    }
}
