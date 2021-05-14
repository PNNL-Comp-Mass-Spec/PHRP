using System;

namespace PeptideHitResultsProcessor.Data
{
    internal class MaxQuantModifiedSequenceModInfo
    {
        // Ignore Spelling: Acetyl

        public string MaxQuantModName { get; }

        public string MatchedModName { get; set; }

        public string MatchedModNameWithParentheses => string.IsNullOrWhiteSpace(MatchedModName) ? string.Empty : string.Format("({0})", MatchedModName);

        /// <summary>
        /// For mods with names like "Oxidation (M)" or "Acetyl (Protein N-term)", return "Oxidation" or "Acetyl"
        /// Otherwise, return the full name
        /// </summary>
        public string GetModNameWithoutResidues()
        {
            // Look for a space followed by an open parenthesis
            var charIndex = MaxQuantModName.IndexOf(" (", StringComparison.Ordinal);
            if (charIndex < 0)
                return MaxQuantModName;

            return MaxQuantModName.Substring(0, charIndex);
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
    }
}
