using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
namespace PeptideHitResultsProcessor.Processor
{
    public static class StringUtilities
    {
        // Ignore Spelling: Da, A-Za-z, Fscore, MaxQuant, prot, mts, MSFragger, xxx
        private static readonly Regex RegexIsLetter = new("[A-Za-z]", RegexOptions.Compiled);

        public static int CIntSafe(string value, int defaultValue)
        {
            try
            {
                // Note: Integer.Parse() fails if value contains a decimal point, even if it is "8.000"
                // Thus, we're converting to a double first, and then rounding
                return (int)Math.Round(Convert.ToDouble(value));
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        public static double CDblSafe(string value, double defaultValue)
        {
            try
            {
                return double.Parse(value);
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        public static float CSngSafe(string value, float defaultValue)
        {
            try
            {
                return float.Parse(value);
            }
            catch (Exception)
            {
                // Error converting value to a number; return the default
                return defaultValue;
            }
        }

        /// <summary>
        /// Collapses a list of strings to a tab-delimited line of text
        /// </summary>
        /// <param name="fields"></param>
        public static string CollapseList(List<string> fields)
        {
            return string.Join("\t", fields);
        }

        /// <summary>
        /// Returns true if the character is a letter between A and Z or a and z
        /// </summary>
        /// <param name="character">Character to examine</param>
        /// <remarks>
        /// Note that the Char.IsLetter() function returns True for "º" and various other Unicode ModifierLetter characters
        /// In contrast, this method only returns True for normal letters between A and Z (case insensitive)
        /// </remarks>
        public static bool IsLetterAtoZ(char character)
        {
            return RegexIsLetter.IsMatch(character.ToString());
        }

        public static string LongestCommonStringFromStart(List<string> items, bool caseSensitive = false)
        {
            if (items.Count == 0)
                return string.Empty;

            if (items.Count == 1)
                return items.First();

            var comparisonType = caseSensitive ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase;

            var longestCommonString = items[0] ?? string.Empty;

            foreach (var item in items.Skip(1))
            {
                longestCommonString = LongestCommonStringFromStart(longestCommonString, item, comparisonType);

                if (longestCommonString.Length == 0)
                {
                    // The items do not all start with the same characters
                    return string.Empty;
                }
            }

            return longestCommonString;
        }

        public static string LongestCommonStringFromStart(string string1, string string2, StringComparison comparisonType = StringComparison.OrdinalIgnoreCase)
        {
            if (string2.Length < string1.Length)
            {
                // Swap strings so that string2 has the longer string
                var temp = string1;
                string1 = string2;
                string2 = temp;
            }

            for (var length = string1.Length; length > 0; length--)
            {
                var startOfString = string1.Substring(0, length);
                if (string2.StartsWith(startOfString, comparisonType))
                {
                    return string1.Substring(0, length);
                }
            }

            return string.Empty;
        }

        public static string LongestCommonStringFromStartVerbose(List<string> items, bool caseSensitive = false)
        {
            if (items.Count == 0)
                return string.Empty;

            if (items.Count == 1)
                return items.First();

            var shortestItem = items.First();

            foreach (var item in items.Skip(1))
            {
                if (item.Length < shortestItem.Length)
                    shortestItem = item;
            }

            var comparisonType = caseSensitive ? StringComparison.Ordinal : StringComparison.OrdinalIgnoreCase;

            var charCount = 1;
            while (charCount <= shortestItem.Length)
            {
                var matchingItems = 0;
                foreach (var item in items)
                {
                    if (item.Equals(shortestItem))
                    {
                        matchingItems++;
                        continue;
                    }

                    if (item.StartsWith(shortestItem.Substring(0, charCount), comparisonType))
                        matchingItems++;
                    else
                        break;
                }

                if (matchingItems < items.Count)
                {
                    break;
                }

                charCount++;
            }

            if (charCount == 1)
            {
                // No text in common
                return string.Empty;
            }

            return shortestItem.Substring(0, charCount - 1);
        }

        public static string MassErrorToString(double massErrorDa)
        {
            if (Math.Abs(massErrorDa) < 0.000001)
                return "0";

            return Math.Abs(massErrorDa) < 0.0001 ?
                       PRISM.StringUtilities.DblToString(massErrorDa, 6, 0.0000001) :
                       PRISM.StringUtilities.DblToString(massErrorDa, 5, 0.000001);
        }

        /// <summary>
        /// If resultID is 0 or 1, returns valueText
        /// Otherwise, if valueText is 0.0, returns 0
        /// Otherwise, returns valueText
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="valueText"></param>
        public static string TrimZeroIfNotFirstID(int resultID, string valueText)
        {
            return resultID > 1 ? TrimZero(valueText) : valueText;
        }

        /// <summary>
        /// If valueText is 0.0, returns 0
        /// Otherwise, returns valueText
        /// </summary>
        /// <param name="valueText"></param>
        public static string TrimZero(string valueText)
        {
            return valueText.Equals("0.0") ? "0" : valueText;
        }
    }
}
