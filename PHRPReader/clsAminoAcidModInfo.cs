namespace PHRPReader
{
    public class clsAminoAcidModInfo
    {
        #region "Constants and Enums"
        public const char N_TERMINAL_PEPTIDE_SYMBOL_DMS = '<';
        public const char C_TERMINAL_PEPTIDE_SYMBOL_DMS = '>';
        public const char N_TERMINAL_PROTEIN_SYMBOL_DMS = '[';
        public const char C_TERMINAL_PROTEIN_SYMBOL_DMS = ']';

        public enum eResidueTerminusStateConstants : int
        {
            None = 0,
            // The residue is in the middle of the peptide
            PeptideNTerminus = 1,
            // The residue is located at the peptide's N-terminus; superseded by ProteinNTerminus if applicable
            PeptideCTerminus = 2,
            // The residue is located at the peptide's C-terminus; superseded by ProteinCTerminus if applicable
            ProteinNTerminus = 3,
            // The residue is located at the protein's N-terminus
            ProteinCTerminus = 4,
            // The residue is located at the protein's C-terminus
            ProteinNandCCTerminus = 5
            // The protein only has one residue
        }
        #endregion

        private readonly clsModificationDefinition mModDefinition;
        private readonly char mResidue;
        private readonly int mResidueLocInPeptide;
        private readonly int mEndResidueLocInPeptide;
        private readonly eResidueTerminusStateConstants mResidueTerminusState;

        /// <summary>
        /// True if the location of the modification is ambiguous
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool AmbiguousMod
        {
            get { return mEndResidueLocInPeptide > mResidueLocInPeptide; }
        }

        /// <summary>
        /// For ambiguous mods, indicates the last residue on which the mod could appear.  For non-ambiguous mods, whill be the same as ResidueLocInPeptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int EndResidueLocInPeptide
        {
            get { return mEndResidueLocInPeptide; }
        }

        public clsModificationDefinition ModDefinition
        {
            get { return mModDefinition; }
        }

        public char Residue
        {
            get { return mResidue; }
        }

        /// <summary>
        /// Indicates the residue number modified; the first residue is at position 1
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>For ambiguous mods, indicates the first residue on which the mod could appear</remarks>
        public int ResidueLocInPeptide
        {
            get { return mResidueLocInPeptide; }
        }

        public eResidueTerminusStateConstants ResidueTerminusState
        {
            get { return mResidueTerminusState; }
        }

        public clsAminoAcidModInfo(char chResidue, int intResidueLocInPeptide, eResidueTerminusStateConstants eResidueTerminusState, clsModificationDefinition objModDefinition)
        {
            mModDefinition = objModDefinition;
            mResidue = chResidue;
            mResidueLocInPeptide = intResidueLocInPeptide;
            mEndResidueLocInPeptide = mResidueLocInPeptide;
            mResidueTerminusState = eResidueTerminusState;
        }

        public clsAminoAcidModInfo(char chResidue, int intResidueLocInPeptide, eResidueTerminusStateConstants eResidueTerminusState, clsModificationDefinition objModDefinition, int intEndResidueLocInPeptide)
        {
            mModDefinition = objModDefinition;
            mResidue = chResidue;
            mResidueLocInPeptide = intResidueLocInPeptide;
            mEndResidueLocInPeptide = intEndResidueLocInPeptide;
            mResidueTerminusState = eResidueTerminusState;
        }
    }
}
