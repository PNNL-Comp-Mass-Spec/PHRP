﻿using System.Collections.Generic;

namespace PHRPReader.Data
{
    /// <summary>
    /// Empirical formula
    /// </summary>
    public class EmpiricalFormula
    {
        /// <summary>
        /// Elements in the empirical formula
        /// Keys are element symbols, values are element counts
        /// </summary>
        public Dictionary<string, int> ElementCounts { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        public EmpiricalFormula()
        {
            ElementCounts = new Dictionary<string, int>();
        }

        /// <summary>
        /// Constructor, initialized with an existing dictionary of element symbols and counts
        /// </summary>
        public EmpiricalFormula(Dictionary<string, int> elementInfo)
        {
            ElementCounts = elementInfo;
        }

        /// <summary>
        /// Constructor, initialized with a list of element symbols
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public EmpiricalFormula(IEnumerable<string> elementInfo)
        {
            ElementCounts = new Dictionary<string, int>();
            foreach (var element in elementInfo)
            {
                AddElement(element, 1);
            }
        }

        /// <summary>
        /// Constructor, initialized with a list of KeyValuePairs of element symbol and element count
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public EmpiricalFormula(IEnumerable<KeyValuePair<string, int>> elementInfo)
        {
            ElementCounts = new Dictionary<string, int>();
            foreach (var element in elementInfo)
            {
                AddElement(element.Key, element.Value);
            }
        }

        /// <summary>
        /// Add a new element to the empirical formula
        /// </summary>
        /// <param name="elementSymbol">Element symbol</param>
        /// <param name="elementCount">Element count</param>
        public void AddElement(string elementSymbol, int elementCount)
        {
            if (ElementCounts.TryGetValue(elementSymbol, out var existingCount))
            {
                ElementCounts[elementSymbol] = existingCount + elementCount;
            }
            else
            {
                ElementCounts.Add(elementSymbol, elementCount);
            }
        }

        /// <summary>
        /// Adds all the elements from the given empirical formula
        /// </summary>
        /// <param name="empiricalFormula">Empirical formula</param>
        public void AddElements(EmpiricalFormula empiricalFormula)
        {
            foreach (var element in empiricalFormula.ElementCounts)
            {
                AddElement(element.Key, element.Value);
            }
        }

        /// <summary>
        /// Return the number of atoms of the given element in the empirical formula
        /// </summary>
        /// <param name="elementSymbol">Element symbol</param>
        /// <returns>Element Count, or 0 if the element is not in ElementCounts</returns>
        public int GetElementCount(char elementSymbol)
        {
            if (ElementCounts.TryGetValue(elementSymbol.ToString(), out var elementCount))
            {
                return elementCount;
            }

            return 0;
        }

        /// <summary>
        /// String representation of the empirical formula
        /// </summary>
        public override string ToString()
        {
            if (ElementCounts.Count == 0)
            {
                return "<Undefined>";
            }

            var formulaDescription = string.Empty;

            foreach (var element in ElementCounts)
            {
                if (element.Value == 1)
                {
                    formulaDescription += element.Key;
                }
                else
                {
                    formulaDescription += element.Key + element.Value;
                }
            }

            return formulaDescription;
        }
    }
}
