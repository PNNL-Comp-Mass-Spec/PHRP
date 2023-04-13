using System;

namespace PeptideHitResultsProcessor.Data
{
    internal class MaxQuantModifiedSequenceModInfo
    {
        // Ignore Spelling: Acetyl

        /// <summary>
        /// MaxQuant modification name, as read from the MaxQuant results file (msms.txt)
        /// </summary>
        /// <remarks>
        /// <para>
        /// Example names:
        ///   Oxidation (M)
        ///   Acetyl (Protein N-term)
        /// </para>
        /// <para>
        /// If the results file has abbreviated modification names, PHRP tries to determine the official modification names, storing the result in MatchedModName
        /// Example abbreviated names:
        ///   ac
        ///   ox
        /// </para>
        /// </remarks>
        public string MaxQuantModName { get; }

        /// <summary>
        /// Official MaxQuant modification name
        /// </summary>
        public string MatchedModName { get; set; }

        /// <summary>
        /// Modification name, surrounded by parentheses
        /// </summary>
        /// <remarks>
        /// Empty string if MatchedModName is not defined
        /// </remarks>
        public string MatchedModNameWithParentheses => string.IsNullOrWhiteSpace(MatchedModName) ? string.Empty : string.Format("({0})", MatchedModName);

        /// <summary>
        /// For mods with names like "Oxidation (M)" or "Acetyl (Protein N-term)", return "Oxidation" or "Acetyl"
        /// Otherwise, return the full name
        /// </summary>
        public string GetModNameWithoutResidues()
        {
            return GetModNameWithoutResidues(MaxQuantModName);
        }

        /// <summary>
        /// For mods with names like "Oxidation (M)" or "Acetyl (Protein N-term)", return "Oxidation" or "Acetyl"
        /// Otherwise, return the full name
        /// </summary>
        public static string GetModNameWithoutResidues(string maxQuantModName)
        {
            // Look for a space followed by an open parenthesis
            var charIndex = maxQuantModName.IndexOf(" (", StringComparison.Ordinal);

            if (charIndex < 0)
                return maxQuantModName;

            return maxQuantModName.Substring(0, charIndex);
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="maxQuantModName"></param>
        public MaxQuantModifiedSequenceModInfo(string maxQuantModName)
        {
            MaxQuantModName = maxQuantModName;
            MatchedModName = string.Empty;
        }

        /// <summary>
        /// Show the modification name
        /// </summary>
        public override string ToString()
        {
            return string.Format("{0}", MaxQuantModName);
        }
    }
}
