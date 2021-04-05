using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class holds rows read from MaxQuant file peptides.txt
    /// </summary>
    internal class MaxQuantPeptideInfo
    {
        // Ignore Spelling: MaxQuant

        public string Sequence;
        public string Prefix;
        public string Suffix;
        public readonly List<string> Proteins = new();
        public string LeadingRazorProtein;
        public string Intensity;
        public readonly Dictionary<string, string> IntensityByExperiment = new();

        /// <summary>
        /// Constructor
        /// </summary>
        public MaxQuantPeptideInfo()
        {
            Clear();
        }

        public void Clear()
        {
            Sequence = string.Empty;
            Prefix = string.Empty;
            Suffix = string.Empty;
            Proteins.Clear();
            LeadingRazorProtein = string.Empty;
            Intensity = string.Empty;
            IntensityByExperiment.Clear();
        }
    }
}
