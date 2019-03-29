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
    /// This class is used to track the peptide details for a MSPathFinder search result
    /// See clsSearchResultsBaseClass for additional information
    /// </summary>
    public class clsSearchResultsMSPathFinder : clsSearchResultsBaseClass
    {
        #region "Classwide Variables"
        // Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables
        #endregion

        #region "Properties"

        // Auto-Properties
        public string MostAbundantIsotopeMz { get; set; }
        public string Modifications { get; set; }
        public string Composition { get; set; }
        public string ProteinDesc { get; set; }
        public string ProteinLength { get; set; }
        public string ResidueStart { get; set; }
        public string ResidueEnd { get; set; }
        public string MatchedFragments { get; set; }
        public string SpecEValue { get; set; }
        public string EValue { get; set; }
        public string QValue { get; set; }
        public string PepQValue { get; set; }

        #endregion

        // Note that the following call will call both the base class's Clear method and this class's Clear method
        public clsSearchResultsMSPathFinder(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        public override void Clear()
        {
            base.Clear();

            MostAbundantIsotopeMz = string.Empty;
            Modifications = string.Empty;
            Composition = string.Empty;
            ProteinDesc = string.Empty;
            ProteinLength = string.Empty;
            ResidueStart = string.Empty;
            ResidueEnd = string.Empty;
            MatchedFragments = string.Empty;
            SpecEValue = string.Empty;
            EValue = string.Empty;
            QValue = string.Empty;
            PepQValue = string.Empty;
        }
    }
}
