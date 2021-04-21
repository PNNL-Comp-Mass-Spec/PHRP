using System.Collections.Generic;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class holds rows read from MaxQuant file peptides.txt
    /// </summary>
    internal class MaxQuantPeptideInfo
    {
        // Ignore Spelling: MaxQuant

        /// <summary>
        /// Id number for a given row in the peptides.txt file
        /// </summary>
        public int PeptideId { get; }

        /// <summary>
        /// Amino acid sequence
        /// </summary>
        public string Sequence;

        /// <summary>
        /// Amino acid in the protein sequence before the peptide
        /// </summary>
        /// <remarks>
        /// Dash if at the start of a protein
        /// Empty string if from a reversed protein
        /// </remarks>
        public string Prefix;

        /// <summary>
        /// Amino acid in the protein sequence after the peptide
        /// </summary>
        /// <remarks>
        /// Dash if at the end of a protein
        /// Empty string if from a reversed protein
        /// </remarks>
        public string Suffix;

        /// <summary>
        /// Protein names that have this peptide sequence
        /// </summary>
        public List<string> Proteins { get; } = new();

        /// <summary>
        /// Name of the best scoring protein this peptide is associated with
        /// </summary>
        /// <remarks>
        /// Typically there is only one protein name here
        /// However, in cases of a tied score, will be a semicolon separated list
        /// </remarks>
        public string LeadingRazorProtein;

        /// <summary>
        /// Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
        /// </summary>
        public string Intensity;

        /// <summary>
        /// Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide, by experiment
        /// </summary>
        /// <remarks>Keys are intensity name, values are the intensity value (an integer)</remarks>
        public Dictionary<string, string> IntensityByExperiment { get; } = new();

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideId">Peptide ID</param>
        public MaxQuantPeptideInfo(int peptideId)
        {
            PeptideId = peptideId;
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
