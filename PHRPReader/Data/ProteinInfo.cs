namespace PHRPReader.Data
{
    /// <summary>
    /// Protein (or peptide) metadata
    /// </summary>
    public class ProteinInfo
    {
        /// <summary>
        /// Protein name
        /// </summary>
        public string ProteinName { get; }

        /// <summary>
        /// Protein description
        /// </summary>
        public string Description { get; }

        /// <summary>
        /// Cleavage state of a protein fragment
        /// </summary>
        public PeptideCleavageStateCalculator.PeptideCleavageStateConstants CleavageState { get; }

        /// <summary>
        /// Residue number in the protein at which this sequence starts
        /// </summary>
        public int ResidueStart { get; private set; }

        /// <summary>
        /// Residue number in the protein at which this sequence ends
        /// </summary>
        public int ResidueEnd { get; private set; }

        /// <summary>
        /// Sequence ID
        /// </summary>
        public int SeqID { get; }

        /// <summary>
        /// Terminus state of a protein fragment
        /// </summary>
        public PeptideCleavageStateCalculator.PeptideTerminusStateConstants TerminusState { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName"></param>
        /// <param name="seqID"></param>
        /// <param name="cleavageState"></param>
        /// <param name="terminusState"></param>
        public ProteinInfo(
            string proteinName,
            int seqID,
            PeptideCleavageStateCalculator.PeptideCleavageStateConstants cleavageState,
            PeptideCleavageStateCalculator.PeptideTerminusStateConstants terminusState) : this(proteinName, string.Empty, seqID, cleavageState, terminusState)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName"></param>
        /// <param name="proteinDescription"></param>
        /// <param name="seqID"></param>
        /// <param name="cleavageState"></param>
        /// <param name="terminusState"></param>
        public ProteinInfo(
            string proteinName,
            string proteinDescription,
            int seqID,
            PeptideCleavageStateCalculator.PeptideCleavageStateConstants cleavageState,
            PeptideCleavageStateCalculator.PeptideTerminusStateConstants terminusState) : this(proteinName, proteinDescription, seqID, cleavageState, terminusState, 0, 0)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName"></param>
        /// <param name="proteinDescription"></param>
        /// <param name="seqID"></param>
        /// <param name="cleavageState"></param>
        /// <param name="terminusState"></param>
        /// <param name="proteinResidueStart"></param>
        /// <param name="proteinResidueEnd"></param>
        public ProteinInfo(
            string proteinName,
            string proteinDescription,
            int seqID,
            PeptideCleavageStateCalculator.PeptideCleavageStateConstants cleavageState,
            PeptideCleavageStateCalculator.PeptideTerminusStateConstants terminusState,
            int proteinResidueStart,
            int proteinResidueEnd)
        {
            ProteinName = proteinName;
            if (string.IsNullOrEmpty(proteinDescription))
            {
                Description = string.Empty;
            }
            else
            {
                Description = proteinDescription;
            }
            SeqID = seqID;
            CleavageState = cleavageState;
            TerminusState = terminusState;

            UpdateLocationInProtein(proteinResidueStart, proteinResidueEnd);
        }

        /// <summary>
        /// Update the start/end residues for this protein (or peptide)
        /// </summary>
        /// <param name="proteinResidueStart"></param>
        /// <param name="proteinResidueEnd"></param>
        public void UpdateLocationInProtein(int proteinResidueStart, int proteinResidueEnd)
        {
            ResidueStart = proteinResidueStart;
            ResidueEnd = proteinResidueEnd;
        }

        /// <summary>
        /// Protein name
        /// </summary>
        public override string ToString()
        {
            return ProteinName;
        }
    }
}
