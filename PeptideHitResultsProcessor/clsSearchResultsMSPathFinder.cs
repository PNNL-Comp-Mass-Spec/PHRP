// This class is used to track the peptide details for a MSPathFinder search result
// See clsSearchResultsBaseClass for additional information
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Created 7/16/2015
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
using PHRPReader;

namespace PeptideHitResultsProcessor
{
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

        // Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        public clsSearchResultsMSPathFinder(clsPeptideModificationContainer objPeptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
            : base(objPeptideMods, peptideSeqMassCalculator)
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
