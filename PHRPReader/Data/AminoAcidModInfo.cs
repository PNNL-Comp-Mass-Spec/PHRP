namespace PHRPReader.Data
{
    /// <summary>
    /// Tracks modifications on a residue
    /// </summary>
    public class AminoAcidModInfo
    {
        #region "Constants and Enums"

        /// <summary>
        /// Symbol used by DMS for tracking the N-terminus of a peptide
        /// </summary>
        public const char N_TERMINAL_PEPTIDE_SYMBOL_DMS = '<';

        /// <summary>
        /// Symbol used by DMS for tracking the C-terminus of a peptide
        /// </summary>
        public const char C_TERMINAL_PEPTIDE_SYMBOL_DMS = '>';

        /// <summary>
        /// Symbol used by DMS for tracking the N-terminus of a protein
        /// </summary>
        public const char N_TERMINAL_PROTEIN_SYMBOL_DMS = '[';

        /// <summary>
        /// Symbol used by DMS for tracking the C-terminus of a protein
        /// </summary>
        public const char C_TERMINAL_PROTEIN_SYMBOL_DMS = ']';

        /// <summary>
        /// Terminus state enum
        /// </summary>
        public enum ResidueTerminusStateConstants
        {
            /// <summary>
            /// The residue is in the middle of the peptide
            /// </summary>
            None = 0,

            /// <summary>
            /// The residue is located at the peptide's N-terminus; superseded by ProteinNTerminus if applicable
            /// </summary>
            PeptideNTerminus = 1,

            /// <summary>
            /// The residue is located at the peptide's C-terminus; superseded by ProteinCTerminus if applicable
            /// </summary>
            PeptideCTerminus = 2,

            /// <summary>
            /// The residue is located at the protein's N-terminus
            /// </summary>
            ProteinNTerminus = 3,

            /// <summary>
            /// The residue is located at the protein's C-terminus
            /// </summary>
            ProteinCTerminus = 4,

            /// <summary>
            /// The protein only has one residue
            /// </summary>
            ProteinNandCCTerminus = 5

        }
        #endregion

        /// <summary>
        /// True if the location of the modification is ambiguous
        /// </summary>
        public bool AmbiguousMod => EndResidueLocInPeptide > ResidueLocInPeptide;

        /// <summary>
        /// For ambiguous mods, indicates the last residue on which the mod could appear.  For non-ambiguous mods, will be the same as ResidueLocInPeptide
        /// </summary>
        public int EndResidueLocInPeptide { get; }

        /// <summary>
        /// Modification definition
        /// </summary>
        public ModificationDefinition ModDefinition { get; }

        /// <summary>
        /// Amino acid residue symbol
        /// </summary>
        public char Residue { get; }

        /// <summary>
        /// Indicates the residue number modified; the first residue is at position 1
        /// </summary>
        /// <remarks>For ambiguous mods, indicates the first residue on which the mod could appear</remarks>
        public int ResidueLocInPeptide { get; }

        /// <summary>
        /// Residue terminus state
        /// </summary>
        public ResidueTerminusStateConstants ResidueTerminusState { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="residue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="residueTerminusState"></param>
        /// <param name="modDefinition"></param>
        public AminoAcidModInfo(char residue, int residueLocInPeptide, ResidueTerminusStateConstants residueTerminusState, ModificationDefinition modDefinition)
        {
            ModDefinition = modDefinition;
            Residue = residue;
            ResidueLocInPeptide = residueLocInPeptide;
            EndResidueLocInPeptide = ResidueLocInPeptide;
            ResidueTerminusState = residueTerminusState;
        }

        /// <summary>
        /// Constructor with endResidueLocInPeptide
        /// </summary>
        /// <param name="residue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="residueTerminusState"></param>
        /// <param name="modDefinition"></param>
        /// <param name="endResidueLocInPeptide"></param>
        public AminoAcidModInfo(char residue, int residueLocInPeptide, ResidueTerminusStateConstants residueTerminusState, ModificationDefinition modDefinition, int endResidueLocInPeptide)
        {
            ModDefinition = modDefinition;
            Residue = residue;
            ResidueLocInPeptide = residueLocInPeptide;
            EndResidueLocInPeptide = endResidueLocInPeptide;
            ResidueTerminusState = residueTerminusState;
        }
    }
}
