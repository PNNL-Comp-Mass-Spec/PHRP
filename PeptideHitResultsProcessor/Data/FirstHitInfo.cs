using System;
using PHRPReader;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// First hits file entry
    /// </summary>
    public class FirstHitInfo
    {
        // Ignore Spelling: Readonly

        /// <summary>
        /// Primary sequence
        /// </summary>
        /// <remarks>
        /// Includes modification symbols, but not the prefix or suffix residue
        /// </remarks>
        // ReSharper disable once UnassignedReadonlyField
        private readonly string mPrimarySequence;

        /// <summary>
        /// Prefix residue
        /// </summary>
        private string mPrefix;

        /// <summary>
        /// Suffix residue
        /// </summary>
        private string mSuffix;

        /// <summary>
        /// Clean sequence (no modification symbols)
        /// </summary>
        public string CleanSequence { get; }

        /// <summary>
        /// Sequence with prefix and suffix residues, plus any modification symbols
        /// </summary>
        public string SequenceWithModsAndContext => mPrefix + "." + mPrimarySequence + "." + mSuffix;

        /// <summary>
        /// Protein name associated with this peptide
        /// </summary>
        public string ProteinName { get; set; }

        /// <summary>
        /// Protein number (unique integer for assigned to each protein)
        /// </summary>
        public int ProteinNumber { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideSeqWithModsAndContext"></param>
        /// <param name="peptideCleanSeq"></param>
        public FirstHitInfo(string peptideSeqWithModsAndContext, string peptideCleanSeq)
        {
            if (!PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptideSeqWithModsAndContext, out mPrimarySequence, out mPrefix, out mSuffix))
            {
                throw new Exception("Unable to split the prefix and suffix from peptide " + peptideSeqWithModsAndContext);
            }

            CleanSequence = peptideCleanSeq;
        }

        /// <summary>
        /// Show clean sequence, protein name, and protein number
        /// </summary>
        public override string ToString()
        {
            return CleanSequence + ", " + ProteinName + ", " + ProteinNumber;
        }

        /// <summary>
        /// Update the prefix and suffix residues for this peptide
        /// </summary>
        /// <param name="newPrefix"></param>
        /// <param name="newSuffix"></param>
        public void UpdatePrefixAndSuffix(string newPrefix, string newSuffix)
        {
            mPrefix = newPrefix;
            mSuffix = newSuffix;
        }
    }
}
