using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PHRPReader
{
    public class Enums
    {
        /// <summary>
        /// Modification types
        /// </summary>
        public enum ResidueModificationType
        {
            /// <summary>
            /// Unknown mod type on a residue; essentially treated as a dynamic mod
            /// </summary>
            UnknownType = 0,

            /// <summary>
            /// Dynamic mod on a residue or peptide terminus; supported by Sequest and notated via a modification symbol; this mod is explicitly notated by X!Tandem; if a terminus mod, the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            DynamicMod = 1,

            /// <summary>
            /// Static mod on a residue or peptide terminus; supported by Sequest but not explicitly notated; this mod is explicitly notated by X!Tandem; if a terminus mod, the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            StaticMod = 2,

            /// <summary>
            /// Peptide terminus static mod (DMS Symbol is T); used by Sequest and MSGFDB; note that terminal mods are always dynamic in X!Tandem
            /// </summary>
            TerminalPeptideStaticMod = 3,

            /// <summary>
            /// Isotopic mod, e.g. N15, or C13; supported by Sequest; most likely not supported by XTandem
            /// </summary>
            IsotopicMod = 4,

            /// <summary>
            /// Protein terminus static mod; supported by Sequest; this mod is also supported by X!Tandem but modified residues are not explicitly notated; instead, all peptides have their mass implicitly modified by this amount
            /// </summary>
            ProteinTerminusStaticMod = 5
        }
    }
}
