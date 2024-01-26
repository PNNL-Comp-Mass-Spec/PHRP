using System;
using System.Text;
using System.Text.RegularExpressions;
using MolecularWeightCalculator.Formula;
using PHRPReader;

namespace PeptideHitResultsProcessor.Data
{
    internal class MaxQuantModInfo
    {
        // Ignore Spelling: Cx, Da, Hx, Nx, MaxQuant, UniMod

        private static readonly MolecularWeightCalculator.MolecularWeightTool mMolecularWeightCalculator = new(ElementMassMode.Isotopic);

        private static readonly PeptideMassCalculator mPeptideSeqMassCalculator = new()
        {
            ChargeCarrierMass = PeptideMassCalculator.MASS_PROTON,
            RemovePrefixAndSuffixIfPresent = false
        };

        private static readonly Regex mNegativeCountMatcher = new(@"(?<Element>[a-z]+)\((?<ElementCount>-\d+)\)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private static readonly Regex mNonAminoAcidCharacterMatcher = new("[^ACDEFGHIKLMNOPQRSTUVWY]+", RegexOptions.Compiled);

        /// <summary>
        /// Empirical formula
        /// </summary>
        public string Composition { get; }

        /// <summary>
        /// Description
        /// </summary>
        public string Description { get; set; }

        /// <summary>
        /// Monoisotopic mass
        /// </summary>
        public double MonoisotopicMass { get; }

        /// <summary>
        /// Residues that this modification can apply to, tracked as a
        /// space-free, comma-free list of one letter amino acid residue symbols
        /// </summary>
        public string Residues { get; set; }

        /// <summary>
        /// Position in the peptide where this mod can occur
        /// </summary>
        public MaxQuantModPosition Position { get; set; }

        /// <summary>
        /// Mod title, which traditionally has the mod name and affected residue, e.g. "Oxidation (M)"
        /// </summary>
        public string Title { get; }

        /// <summary>
        /// Modification type
        /// </summary>
        public MaxQuantModType ModType { get; set; }

        /// <summary>
        /// Terminus type
        /// </summary>
        [Obsolete("Deprecated with MaxQuant v2.4.0")]
        public MaxQuantTerminusType TerminusType { get; set; }

        /// <summary>
        /// Constructor that takes a mod name and empirical formula
        /// </summary>
        /// <param name="title">Modification name</param>
        /// <param name="composition">Empirical formula</param>
        [Obsolete("Use the constructor that accepts Mod Type")]
        public MaxQuantModInfo(string title, string composition)
        {
            Title = title;
            Composition = composition;
            Position = MaxQuantModPosition.Anywhere;
            ModType = MaxQuantModType.Standard;
            TerminusType = MaxQuantTerminusType.None;
            Residues = string.Empty;
            MonoisotopicMass = ComputeMass(title, composition, ModType);
        }

        /// <summary>
        /// Constructor that takes a mod name and empirical formula
        /// </summary>
        /// <param name="title">Modification name</param>
        /// <param name="composition">Empirical formula</param>
        /// <param name="modType">MaxQuant mod type</param>
        public MaxQuantModInfo(string title, string composition, MaxQuantModType modType)
        {
            Title = title;
            Composition = composition.Trim();
            Position = MaxQuantModPosition.Anywhere;
            ModType = MaxQuantModType.Standard;
            Residues = string.Empty;
            MonoisotopicMass = ComputeMass(title, composition, modType);
        }

        /// <summary>
        /// Constructor that takes a mod name and mass
        /// </summary>
        /// <param name="title">Modification name</param>
        /// <param name="monoisotopicMass">Monoisotopic mass, in Da</param>
        public MaxQuantModInfo(string title, double monoisotopicMass)
        {
            Title = title;
            Composition = string.Empty;
            Position = MaxQuantModPosition.Anywhere;
            ModType = MaxQuantModType.Standard;
            Residues = string.Empty;
            MonoisotopicMass = monoisotopicMass;
        }

        /// <summary>
        /// Parse the empirical formula to compute the modification mass
        /// </summary>
        /// <param name="title">Modification name</param>
        /// <param name="composition">Empirical formula</param>
        /// <param name="modType">MaxQuant mod type</param>
        /// <returns>Monoisotopic mass, in Da</returns>
        public static double ComputeMass(string title, string composition, MaxQuantModType modType)
        {
            // Examples formulas for composition:

            // C(2) H(2) O
            // H(3) O(4) P
            // Cx(5) Nx C(-5) N(-1)
            // H(-1) N(-1) Ox
            // C(-6) Cx(6) N(-1) Nx

            // MaxQuant v2.4.13 introduced modType SequenceBasedModifier, which indicates that the composition is amino acid based, e.g.
            // GG (two glycines)
            // RGG (arginine and two glycines)
            // QQTGG
            // DVFQQQTGG

            if (modType == MaxQuantModType.SequenceBasedModifier)
            {
                // Amino acid based formula
                // Check for any unrecognized characters

                var unexpectedCharacter = mNonAminoAcidCharacterMatcher.Match(composition);

                if (unexpectedCharacter.Success)
                {
                    throw new Exception(string.Format(
                        "Error computing modification mass for '{0}', composition {1}, unexpected character found: {2}",
                        title, composition, unexpectedCharacter.Value));
                }

                // By default, the peptide sequence mass calculator adds the mass of H and OH when computing the mass of a peptide
                // We don't want that for MaxQuant mods, so assure that the N and C terminus mass values are zero

                if (mPeptideSeqMassCalculator.PeptideNTerminusMass != 0)
                {
                    mPeptideSeqMassCalculator.PeptideNTerminusMass = 0;
                }

                if (mPeptideSeqMassCalculator.PeptideCTerminusMass != 0)
                {
                    mPeptideSeqMassCalculator.PeptideCTerminusMass = 0;
                }

                return mPeptideSeqMassCalculator.ComputeSequenceMass(composition);
            }

            // Look for elements (or groups) with a negative element count
            // The Molecular Weight Calculator does not support negative counts, so we'll need to process these elements separately
            var negativeCountMatches = mNegativeCountMatcher.Matches(composition);

            string updatedFormula1;
            double massToSubtract;
            if (negativeCountMatches.Count > 0)
            {
                var formulaToSubtract = new StringBuilder();

                foreach (Match item in negativeCountMatches)
                {
                    var elementCount = int.Parse(item.Groups["ElementCount"].Value);

                    formulaToSubtract.AppendFormat("{0}{1}", item.Groups["Element"].Value, Math.Abs(elementCount));
                }

                var updatedFormulaToSubtract = ConvertFormulaNotation(formulaToSubtract.ToString());

                massToSubtract = mMolecularWeightCalculator.ComputeMass(updatedFormulaToSubtract);

                if (massToSubtract == 0)
                {
                    throw new Exception(string.Format(
                        "Error computing modification mass for '{0}', formula {1}, subtracting {2}: {3}",
                        title, composition, formulaToSubtract, mMolecularWeightCalculator.ErrorDescription));
                }

                updatedFormula1 = mNegativeCountMatcher.Replace(composition, string.Empty).Trim();

                if (string.IsNullOrWhiteSpace(updatedFormula1))
                {
                    // This modification only subtracts a mass
                    return -massToSubtract;
                }
            }
            else
            {
                massToSubtract = 0;
                updatedFormula1 = composition;
            }

            var updatedFormula2 = ConvertFormulaNotation(updatedFormula1);

            var computedMass = mMolecularWeightCalculator.ComputeMass(updatedFormula2);

            if (computedMass != 0)
            {
                return computedMass - massToSubtract;
            }

            throw new Exception(string.Format(
                "Error computing modification mass for '{0}', formula {1}: {2}",
                title, composition, mMolecularWeightCalculator.ErrorDescription));
        }

        /// <summary>
        /// The Molecular Weight Calculator can properly parse UniMod style formulas, e.g. H(25) C(8) 13C(7) N 15N(2) O(3)
        /// However, MaxQuant uses Cx, Nx, Ox, and Hx for heavy isotopes of elements
        /// This method replaces those symbols with the format supported by the Molecular Weight Calculator
        /// </summary>
        /// <param name="formula">Empirical formula</param>
        private static string ConvertFormulaNotation(string formula)
        {
            return formula.Replace("Cx", "^13.003355C").Replace("Nx", "^15.000109N").Replace("Ox", "^17.999161O").Replace("Hx", "D");
        }

        /// <summary>
        /// Show the mod name and mass
        /// </summary>
        public override string ToString()
        {
            return string.Format("{0}: {1} Da", Title, MonoisotopicMass);
        }
    }
}
