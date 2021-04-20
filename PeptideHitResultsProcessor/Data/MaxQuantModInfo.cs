using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;
using MolecularWeightCalculator.Formula;

namespace PeptideHitResultsProcessor.Data
{
    internal class MaxQuantModInfo
    {
        // Ignore Spelling: Cx, Da, Hx, Nx, UniMod

        private static readonly MolecularWeightCalculator.MolecularWeightTool mMolecularWeightCalculator = new(ElementMassMode.Isotopic);

        private static readonly Regex mNegativeCountMatcher = new(@"(?<Element>[a-z]+)\((?<ElementCount>-\d+)\)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        public string Composition { get; }

        public string Description { get; set; }

        public double MonoisotopicMass { get; }

        public string Residues { get; set; }

        public MaxQuantModPosition Position { get; set; }

        public string Title { get; }

        public MaxQuantModType ModType { get; set; }

        public MaxQuantTerminusType TerminusType { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="title"></param>
        /// <param name="composition"></param>
        public MaxQuantModInfo(string title, string composition)
        {
            Title = title;
            Composition = composition;
            Position = MaxQuantModPosition.Anywhere;
            ModType = MaxQuantModType.Standard;
            TerminusType = MaxQuantTerminusType.None;

            Residues = string.Empty;

            MonoisotopicMass = ComputeMass(title, composition);
        }

        public static double ComputeMass(string title, string composition)
        {
            // Examples formulas for composition:

            // C(2) H(2) O
            // H(3) O(4) P
            // Cx(5) Nx C(-5) N(-1)
            // H(-1) N(-1) Ox
            // C(-6) Cx(6) N(-1) Nx

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

                massToSubtract = mMolecularWeightCalculator.ComputeMass(formulaToSubtract.ToString());

                if (massToSubtract == 0)
                {
                    throw new Exception(string.Format(
                        "Error computing modification mass for '{0}', formula {1}, subtracting {2}: {3}",
                        title, composition, formulaToSubtract, mMolecularWeightCalculator.ErrorDescription));
                }

                updatedFormula1 = mNegativeCountMatcher.Replace(composition, string.Empty);

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
        /// This method replaces those symbols with the format support by the Molecular Weight Calculator
        /// </summary>
        /// <param name="formula"></param>
        /// <returns></returns>
        private static string ConvertFormulaNotation(string formula)
        {
            return formula.Replace("Cx", "^13.003355C").Replace("Nx", "^15.000109N").Replace("Ox", "^17.999161O").Replace("Hx", "D");
        }

        public override string ToString()
        {
            return string.Format("{0}: {1} Da", Title, MonoisotopicMass);
        }
    }
}
