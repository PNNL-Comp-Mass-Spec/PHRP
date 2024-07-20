// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
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
    /// This class is used to track the peptide details for a MSPathFinder search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class MSPathFinderResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Dehydro

        /// <summary>
        /// Most abundant isotope m/z
        /// </summary>
        public string MostAbundantIsotopeMz { get; set; }

        /// <summary>
        /// Comma-separated list of modification names and affected residue number
        /// </summary>
        /// <remarks>
        /// Example values:
        ///   Oxidation 7,Dehydro 16
        ///   Dehydro 1,Dehydro 4,Dehydro 7
        /// </remarks>
        public string Modifications { get; set; }

        /// <summary>
        /// Empirical formula
        /// </summary>
        /// <remarks>
        /// Example: C(74) H(119) N(19) O(26) S(0)
        /// </remarks>
        public string Composition { get; set; }

        /// <summary>
        /// Protein description
        /// </summary>
        public string ProteinDesc { get; set; }

        /// <summary>
        /// Number of residues in the full protein
        /// </summary>
        public string ProteinLength { get; set; }

        /// <summary>
        /// Residue number in the protein where this PSM starts
        /// </summary>
        public string ResidueStart { get; set; }

        /// <summary>
        /// Residue number in the protein where this PSM ends
        /// </summary>
        public string ResidueEnd { get; set; }

        /// <summary>
        /// Number of matched fragment ions
        /// </summary>
        public string MatchedFragments { get; set; }

        /// <summary>
        /// Spectrum-level E-value for this PSM
        /// </summary>
        public string SpecEValue { get; set; }

        /// <summary>
        /// Dataset-wide E-value for this PSM
        /// </summary>
        public string EValue { get; set; }

        /// <summary>
        /// Q-value (FDR) for this PSM
        /// </summary>
        /// <remarks>
        /// Minimum false discovery rate (FDR) at which a test may be called significant
        /// </remarks>
        public string QValue { get; set; }

        /// <summary>
        /// Peptide-level QValue (FDR) estimated using the target-decoy approach
        /// </summary>
        public string PepQValue { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public MSPathFinderResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        /// <summary>
        /// Clear stored values
        /// </summary>
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
