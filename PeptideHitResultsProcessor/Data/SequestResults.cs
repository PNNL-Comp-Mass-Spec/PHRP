// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 7, 2006
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

using System;
using PHRPReader;
using PHRPReader.Data;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class is used to track the peptide details for a Sequest search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    public class clsSearchResultsSequest : clsSearchResultsBaseClass
    {
        #region "Properties"

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

        /// <summary>
        /// Not stored in the synopsis or first hits file; it can be computed from XCorr and DeltaCn2
        /// </summary>
        public string PeptideXCorrNext { get; private set; }

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
        public clsSearchResultsSequest(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
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

            PeptideXCorrNext = string.Empty;
        }

        private void ComputePeptideXCorrNext()
        {
            try
            {
                if (float.TryParse(PeptideXCorr, out var xCorr) &&
                    float.TryParse(PeptideDeltaCn2, out var delCN2))
                {
                    PeptideXCorrNext = (xCorr - delCN2 * xCorr).ToString();
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