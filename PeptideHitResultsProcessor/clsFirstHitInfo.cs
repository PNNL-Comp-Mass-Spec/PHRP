using System;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsFirstHitInfo
    {
        // ReSharper disable once UnassignedReadonlyField
        private readonly string mPrimarySequence;

        private string mPrefix;
        private string mSuffix;

        public string CleanSequence { get; }

        public string SequenceWithModsAndContext
        {
            get { return mPrefix + "." + mPrimarySequence + "." + mSuffix; }
        }

        public string ProteinName { get; set; }

        public int ProteinNumber { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideSeqWithModsAndContext"></param>
        /// <param name="peptideCleanSeq"></param>
        public clsFirstHitInfo(string peptideSeqWithModsAndContext, string peptideCleanSeq)
        {
            if (!clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptideSeqWithModsAndContext, out mPrimarySequence, out mPrefix, out mSuffix))
            {
                throw new Exception("Unable to split the prefix and suffix from peptide " + peptideSeqWithModsAndContext);
            }

            CleanSequence = peptideCleanSeq;
        }

        public override string ToString()
        {
            return CleanSequence + ", " + ProteinName + ", " + ProteinNumber;
        }

        public void UpdatePrefixAndSuffix(string newPrefix, string newSuffix)
        {
            mPrefix = newPrefix;
            mSuffix = newSuffix;
        }
    }
}
