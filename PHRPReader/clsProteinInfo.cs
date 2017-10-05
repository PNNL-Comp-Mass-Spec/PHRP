namespace PHRPReader
{
    public class clsProteinInfo
    {
        private readonly string mProteinName;
        private string mProteinDescription;
        private readonly int mSeqID;
        private readonly clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants mCleavageState;
        private readonly clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants mTerminusState;
        private int mResidueStart;                         // Residue number in the protein at which this sequence starts
        private int mResidueEnd;                           // Residue number in the protein at which this sequence ends

        public string ProteinName
        {
            get { return mProteinName; }
        }

        // ReSharper disable once ConvertToAutoProperty
        public clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants CleavageState
        {
            get { return mCleavageState; }
        }

        public int ResidueStart
        {
            get { return mResidueStart; }
        }

        public int ResidueEnd
        {
            get { return mResidueEnd; }
        }

        // ReSharper disable once ConvertToAutoProperty
        public int SeqID
        {
            get { return mSeqID; }
        }

        // ReSharper disable once ConvertToAutoProperty
        public clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants TerminusState
        {
            get { return mTerminusState; }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="ProteinName"></param>
        /// <param name="SeqID"></param>
        /// <param name="CleavageState"></param>
        /// <param name="TerminusState"></param>
        /// <remarks></remarks>
        public clsProteinInfo(string ProteinName, int SeqID, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants CleavageState, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants TerminusState) : this(ProteinName, string.Empty, SeqID, CleavageState, TerminusState)
        {
        }

        public clsProteinInfo(string ProteinName, string ProteinDescription, int SeqID, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants CleavageState, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants TerminusState) : this(ProteinName, string.Empty, SeqID, CleavageState, TerminusState, 0, 0)
        {
        }

        public clsProteinInfo(string ProteinName, string ProteinDescription, int SeqID, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants CleavageState, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants TerminusState, int ProteinResidueStart, int ProteinResidueEnd)
        {
            mProteinName = ProteinName;
            if (string.IsNullOrEmpty(ProteinDescription))
            {
                mProteinDescription = string.Empty;
            }
            else
            {
                mProteinDescription = ProteinDescription;
            }
            mSeqID = SeqID;
            mCleavageState = CleavageState;
            mTerminusState = TerminusState;

            UpdateLocationInProtein(ProteinResidueStart, ProteinResidueEnd);
        }

        public void UpdateLocationInProtein(int ProteinResidueStart, int ProteinResidueEnd)
        {
            mResidueStart = ProteinResidueStart;
            mResidueEnd = ProteinResidueEnd;
        }

        public override string ToString()
        {
            return mProteinName;
        }
    }
}
