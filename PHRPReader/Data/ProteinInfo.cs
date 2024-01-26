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
        public PeptideCleavageStateCalculator.PeptideCleavageState CleavageState { get; }

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
        public PeptideCleavageStateCalculator.PeptideTerminusState TerminusState { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName">Protein name</param>
        /// <param name="seqID">Protein sequence ID</param>
        /// <param name="cleavageState">Cleavage state</param>
        /// <param name="terminusState">Terminus state</param>
        public ProteinInfo(
            string proteinName,
            int seqID,
            PeptideCleavageStateCalculator.PeptideCleavageState cleavageState,
            PeptideCleavageStateCalculator.PeptideTerminusState terminusState) : this(proteinName, string.Empty, seqID, cleavageState, terminusState)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName">Protein name</param>
        /// <param name="proteinDescription">Protein description</param>
        /// <param name="seqID">Protein sequence ID</param>
        /// <param name="cleavageState">Cleavage state</param>
        /// <param name="terminusState">Terminus state</param>
        public ProteinInfo(
            string proteinName,
            string proteinDescription,
            int seqID,
            PeptideCleavageStateCalculator.PeptideCleavageState cleavageState,
            PeptideCleavageStateCalculator.PeptideTerminusState terminusState) : this(proteinName, proteinDescription, seqID, cleavageState, terminusState, 0, 0)
        {
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName">Protein name</param>
        /// <param name="proteinDescription">Protein description</param>
        /// <param name="seqID">Protein sequence ID</param>
        /// <param name="cleavageState">Cleavage state</param>
        /// <param name="terminusState">Terminus state</param>
        /// <param name="proteinResidueStart">Residue number in the protein at which this sequence starts</param>
        /// <param name="proteinResidueEnd">Residue number in the protein at which this sequence ends</param>
        public ProteinInfo(
            string proteinName,
            string proteinDescription,
            int seqID,
            PeptideCleavageStateCalculator.PeptideCleavageState cleavageState,
            PeptideCleavageStateCalculator.PeptideTerminusState terminusState,
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
        /// <param name="proteinResidueStart">Residue number in the protein at which this sequence starts</param>
        /// <param name="proteinResidueEnd">Residue number in the protein at which this sequence ends</param>
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
