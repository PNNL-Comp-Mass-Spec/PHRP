namespace PHRPReader
{
    public class clsAminoAcidModInfo
    {
        #region "Constants and Enums"
        public const char N_TERMINAL_PEPTIDE_SYMBOL_DMS = '<';
        public const char C_TERMINAL_PEPTIDE_SYMBOL_DMS = '>';
        public const char N_TERMINAL_PROTEIN_SYMBOL_DMS = '[';
        public const char C_TERMINAL_PROTEIN_SYMBOL_DMS = ']';

        public enum eResidueTerminusStateConstants
        {
            // The residue is in the middle of the peptide
            None = 0,

            // The residue is located at the peptide's N-terminus; superseded by ProteinNTerminus if applicable
            PeptideNTerminus = 1,

            // The residue is located at the peptide's C-terminus; superseded by ProteinCTerminus if applicable
            PeptideCTerminus = 2,

            // The residue is located at the protein's N-terminus
            ProteinNTerminus = 3,

            // The residue is located at the protein's C-terminus
            ProteinCTerminus = 4,

            // The protein only has one residue
            ProteinNandCCTerminus = 5

        }
        #endregion

        /// <summary>
        /// True if the location of the modification is ambiguous
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool AmbiguousMod => EndResidueLocInPeptide > ResidueLocInPeptide;

        /// <summary>
        /// For ambiguous mods, indicates the last residue on which the mod could appear.  For non-ambiguous mods, whill be the same as ResidueLocInPeptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int EndResidueLocInPeptide { get; }

        public clsModificationDefinition ModDefinition { get; }

        public char Residue { get; }

        /// <summary>
        /// Indicates the residue number modified; the first residue is at position 1
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>For ambiguous mods, indicates the first residue on which the mod could appear</remarks>
        public int ResidueLocInPeptide { get; }

        public eResidueTerminusStateConstants ResidueTerminusState { get; }

        public clsAminoAcidModInfo(char chResidue, int intResidueLocInPeptide, eResidueTerminusStateConstants eResidueTerminusState, clsModificationDefinition objModDefinition)
        {
            ModDefinition = objModDefinition;
            Residue = chResidue;
            ResidueLocInPeptide = intResidueLocInPeptide;
            EndResidueLocInPeptide = ResidueLocInPeptide;
            ResidueTerminusState = eResidueTerminusState;
        }

        public clsAminoAcidModInfo(char chResidue, int intResidueLocInPeptide, eResidueTerminusStateConstants eResidueTerminusState, clsModificationDefinition objModDefinition, int intEndResidueLocInPeptide)
        {
            ModDefinition = objModDefinition;
            Residue = chResidue;
            ResidueLocInPeptide = intResidueLocInPeptide;
            EndResidueLocInPeptide = intEndResidueLocInPeptide;
            ResidueTerminusState = eResidueTerminusState;
        }
    }
}
