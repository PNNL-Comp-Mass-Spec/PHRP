// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or https://www.pnnl.gov/sysbio/ or https://panomics.pnnl.gov/
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.Globalization;

namespace PHRPReader.Data
{
    /// <summary>
    /// This class describes an amino acid modification
    /// </summary>
    public class ModificationDefinition
    {
        // Ignore Spelling: UnkMod, Phosph, IodoAcet, Plus1Oxy

        /// <summary>
        /// Modification symbol used after all of the DEFAULT_MODIFICATION_SYMBOLS have been used
        /// </summary>
        public const char LAST_RESORT_MODIFICATION_SYMBOL = '_';

        /// <summary>
        /// Symbol to indicate a modification does not have a mod symbol
        /// </summary>
        /// <remarks>Used with isotopic mods and protein terminus static mods</remarks>
        public const char NO_SYMBOL_MODIFICATION_SYMBOL = '-';

        /// <summary>
        /// Unknown mod base name
        /// </summary>
        public const string UNKNOWN_MOD_BASE_NAME = "UnkMod";

        /// <summary>
        /// Initial unknown mass correction tag name
        /// </summary>
        public const string INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME = UNKNOWN_MOD_BASE_NAME + "00";

        /// <summary>
        /// Legacy modification type constants
        /// </summary>
        [Obsolete("Superseded by ResidueModificationType")]
        public enum ModificationTypeConstants
        {
            /// <summary>
            /// Unknown mod type on a residue; essentially treated as a dynamic mod
            /// </summary>
            UnknownType = 0,

            /// <summary>
            /// Dynamic mod on a residue or peptide terminus; supported by SEQUEST and notated via a modification symbol; this mod is explicitly notated by X!Tandem; if a terminus mod, the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            DynamicMod = 1,

            /// <summary>
            /// Static mod on a residue or peptide terminus; supported by SEQUEST but not explicitly notated; this mod is explicitly notated by X!Tandem; if a terminus mod, the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            StaticMod = 2,

            /// <summary>
            /// Peptide terminus static mod (DMS Symbol is T); used by SEQUEST and MSGFDB; note that terminal mods are always dynamic in X!Tandem
            /// </summary>
            TerminalPeptideStaticMod = 3,

            /// <summary>
            /// Isotopic mod, e.g. N15, or C13; supported by SEQUEST; most likely not supported by XTandem
            /// </summary>
            IsotopicMod = 4,

            /// <summary>
            /// Protein terminus static mod; supported by SEQUEST; this mod is also supported by X!Tandem but modified residues are not explicitly notated; instead, all peptides have their mass implicitly modified by this amount
            /// </summary>
            ProteinTerminusStaticMod = 5
        }

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
            /// Dynamic mod on a residue or peptide terminus; supported by SEQUEST and notated via a modification symbol; this mod is explicitly notated by X!Tandem; if a terminus mod, the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            DynamicMod = 1,

            /// <summary>
            /// Static mod on a residue or peptide terminus; supported by SEQUEST but not explicitly notated; this mod is explicitly notated by X!Tandem; if a terminus mod, the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            StaticMod = 2,

            /// <summary>
            /// Peptide terminus static mod (DMS Symbol is T); used by SEQUEST and MSGFDB; note that terminal mods are always dynamic in X!Tandem
            /// </summary>
            TerminalPeptideStaticMod = 3,

            /// <summary>
            /// Isotopic mod, e.g. N15, or C13; supported by SEQUEST; most likely not supported by XTandem
            /// </summary>
            IsotopicMod = 4,

            /// <summary>
            /// Protein terminus static mod; supported by SEQUEST; this mod is also supported by X!Tandem but modified residues are not explicitly notated; instead, all peptides have their mass implicitly modified by this amount
            /// </summary>
            ProteinTerminusStaticMod = 5
        }

        /// <summary>
        /// One letter symbol for this modification; use NO_SYMBOL_MODIFICATION_SYMBOL if no symbol (necessary for isotopic mods or protein terminus static mods)
        /// </summary>
        private char mModificationSymbol;

        /// <summary>
        /// Monoisotopic modification mass
        /// </summary>
        private double mModificationMass;

        /// <summary>
        /// Target residues, tracked as a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to
        /// Use the *_SYMBOL_DMS constants for the peptide and protein termini symbols (&lt; and &gt; for the peptide termini; [ and ] for the protein termini)
        /// </summary>
        /// <remarks>
        /// If this is empty, the given modification can apply to any residue or terminus
        /// </remarks>
        private string mTargetResidues;

        /// <summary>
        /// Name associated with the given ModificationMass; maximum length is 8 characters
        /// Cannot contain a colon, comma, or space
        /// </summary>
        private string mMassCorrectionTag;

        /// <summary>
        /// Set to Nothing or to PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications)
        /// For Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
        /// </summary>
        private char mAffectedAtom;

        /// <summary>
        /// One letter symbol for this modification
        /// </summary>
        /// <remarks>
        /// Use NO_SYMBOL_MODIFICATION_SYMBOL (a dash) if no symbol
        /// (necessary for isotopic mods or protein terminus static mods)
        /// </remarks>
        public char ModificationSymbol
        {
            get => mModificationSymbol;
            set
            {
                if (value == default(char))
                {
                    mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL;
                }
                else
                {
                    mModificationSymbol = value;
                }
            }
        }

        /// <summary>
        /// Monoisotopic modification mass
        /// </summary>
        public double ModificationMass
        {
            get => mModificationMass;
            set
            {
                mModificationMass = value;
                ModificationMassAsText = mModificationMass.ToString(CultureInfo.InvariantCulture);
            }
        }

        /// <summary>
        /// Modification mass, stored as text
        /// </summary>
        /// <remarks>Represents the original string value read from the data file</remarks>
        public string ModificationMassAsText { get; set; }

        /// <summary>
        /// Residues that this modification can apply to
        /// </summary>
        /// <remarks>
        /// If an empty string, the modification can apply to any residue or terminus;
        /// Otherwise, should contain a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to.
        /// Use the *_SYMBOL_DMS constants for the peptide and protein termini symbols
        /// (less than and greater than signs for the peptide termini; [ and ] for the protein termini)
        /// </remarks>
        public string TargetResidues
        {
            get => mTargetResidues;
            set => mTargetResidues = value ?? string.Empty;
        }

        /// <summary>
        /// Modification type
        /// </summary>
        public ResidueModificationType ModificationType { get; set; }

        /// <summary>
        /// Modification name, for example Phosph, IodoAcet, Plus1Oxy, or Methyl
        /// </summary>
        /// <remarks>Maximum length is 8 characters; cannot contain a colon, comma, or space</remarks>
        public string MassCorrectionTag
        {
            get => mMassCorrectionTag;
            set => mMassCorrectionTag = value ?? string.Empty;
        }

        /// <summary>
        /// Only used with Isotopic modifications, indicating the atom affected (e.g. C, H, N, O, or S)
        /// </summary>
        /// <remarks>
        /// Set to Nothing or to PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL (a dash) for positional modifications
        /// (including terminus modifications)
        /// </remarks>
        public char AffectedAtom
        {
            get => mAffectedAtom;
            set
            {
                if (value == default(char))
                {
                    mAffectedAtom = PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
                }
                else
                {
                    mAffectedAtom = value;
                }
            }
        }

        /// <summary>
        /// Number of times this modification was observed in the given dataset
        /// </summary>
        public int OccurrenceCount { get; set; }

        /// <summary>
        /// True if this was an unknown mass that was auto defined
        /// </summary>
        public bool UnknownModAutoDefined { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public ModificationDefinition()
        {
            Clear();
        }

        /// <summary>
        /// Constructor that takes a mod symbol and mod mass
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="modificationMass"></param>
        public ModificationDefinition(char modificationSymbol, double modificationMass)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
        }

        /// <summary>
        /// Constructor that takes a mod mass, target residues, and modification type
        /// </summary>
        /// <param name="modificationMass"></param>
        /// <param name="targetResidues"></param>
        /// <param name="modificationType"></param>
        public ModificationDefinition(double modificationMass, string targetResidues, ModificationDefinition.ResidueModificationType modificationType)
        {
            Clear();

            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = modificationType;
        }

        /// <summary>
        /// Constructor that takes a mod symbol, mod mass, target residues, modification type, and mass correction tag
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="modificationMass"></param>
        /// <param name="targetResidues"></param>
        /// <param name="modificationType"></param>
        /// <param name="massCorrectionTag"></param>
        public ModificationDefinition(char modificationSymbol, double modificationMass, string targetResidues, ModificationDefinition.ResidueModificationType modificationType, string massCorrectionTag)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = modificationType;
            MassCorrectionTag = massCorrectionTag;
        }

        /// <summary>
        /// Constructor that takes a mod symbol, mod mass, target residues, modification type, mass correction tag, and affected atom
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="modificationMass"></param>
        /// <param name="targetResidues"></param>
        /// <param name="modificationType"></param>
        /// <param name="massCorrectionTag"></param>
        /// <param name="chAffectedAtom"></param>
        /// <param name="unknownModAutoDefined"></param>
        public ModificationDefinition(char modificationSymbol, double modificationMass, string targetResidues, ModificationDefinition.ResidueModificationType modificationType, string massCorrectionTag, char chAffectedAtom, bool unknownModAutoDefined)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = modificationType;
            MassCorrectionTag = massCorrectionTag;
            AffectedAtom = chAffectedAtom;
            UnknownModAutoDefined = unknownModAutoDefined;
        }

        /// <summary>
        /// Initialize the modification definition
        /// </summary>
        public void Clear()
        {
            mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL;
            mModificationMass = 0;
            ModificationMassAsText = "0";
            mTargetResidues = string.Empty;
            ModificationType = ModificationDefinition.ResidueModificationType.UnknownType;
            mMassCorrectionTag = INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME;
            mAffectedAtom = PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
            OccurrenceCount = 0;
            UnknownModAutoDefined = false;
        }

        /// <summary>
        /// Compares b to this object, ignoring .ModificationSymbol and ignoring .AffectedResidues
        /// </summary>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        public bool EquivalentMassTypeTagAndAtom(ModificationDefinition b)
        {
            return EquivalentMassTypeTagAndAtom(this, b);
        }

        /// <summary>
        /// Compare a to b but ignore .ModificationSymbol and .AffectedResidues
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        public bool EquivalentMassTypeTagAndAtom(ModificationDefinition a, ModificationDefinition b)
        {
            var equivalent =
                Math.Abs(Math.Round(a.ModificationMass - b.ModificationMass, PeptideModificationContainer.MASS_DIGITS_OF_PRECISION) - 0) < float.Epsilon &&
                a.ModificationType == b.ModificationType &&
                a.MassCorrectionTag == b.MassCorrectionTag &&
                a.AffectedAtom == b.AffectedAtom;

            return equivalent;
        }

        /// <summary>
        /// Compares b to this object, ignoring .ModificationSymbol
        /// </summary>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        public bool EquivalentMassTypeTagAtomAndResidues(ModificationDefinition b)
        {
            return EquivalentMassTypeTagAtomAndResidues(this, b);
        }

        /// <summary>
        /// Compares b to this object
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        public bool EquivalentMassTypeTagAtomAndResidues(ModificationDefinition a, ModificationDefinition b)
        {
            // First compare a to b but ignore .ModificationSymbol and .AffectedResidues
            var equivalent = EquivalentMassTypeTagAndAtom(a, b);

            if (equivalent)
            {
                // Mass, ModificationType, MassCorrectionTag, and AffectedAtom are identical
                // What about the residues?

                equivalent = false;
                if (a.TargetResidues == null && b.TargetResidues == null)
                {
                    equivalent = true;
                }
                else if (a.TargetResidues != null && b.TargetResidues != null)
                {
                    if (a.ModificationType is ModificationDefinition.ResidueModificationType.DynamicMod or ModificationDefinition.ResidueModificationType.StaticMod)
                    {
                        // Matching dynamic or static modification definitions
                        // Make sure each of the residues in b.TargetResidues is present in .TargetResidues
                        if (EquivalentTargetResidues(a.TargetResidues, b.TargetResidues, false))
                        {
                            equivalent = true;
                        }
                    }
                    else
                    {
                        // Not a dynamic or static mod; require identical target residue lists in order to flag them as identical
                        if (a.TargetResidues == b.TargetResidues)
                        {
                            equivalent = true;
                        }
                    }
                }
            }

            return equivalent;
        }

        /// <summary>
        /// Compare the residue lists (ignoring order)
        /// </summary>
        /// <param name="residues1"></param>
        /// <param name="residues2"></param>
        /// <param name="allowResidues2ToBeSubsetOfResidues1"></param>
        /// <returns>True if they contain the same residues</returns>
        public static bool EquivalentTargetResidues(string residues1, string residues2, bool allowResidues2ToBeSubsetOfResidues1)
        {
            var equivalent = false;

            if (residues1 == null && residues2 == null)
            {
                // Both residues lists are blank
                equivalent = true;
            }
            else if (!(residues1 == null || residues2 == null))
            {
                if (residues1.Length == 0 && residues2.Length == 0)
                {
                    equivalent = true;
                }
                else if (residues1 == residues2)
                {
                    equivalent = true;
                }
                else if (residues1.Length >= residues2.Length)
                {
                    // See if each of the residues in residues2 is in residues1
                    var matchCount = 0;
                    foreach (var chChar in residues2)
                    {
                        if (residues1.IndexOf(chChar) >= 0)
                        {
                            matchCount++;
                        }
                        else
                        {
                            break;
                        }
                    }

                    if (matchCount == residues1.Length)
                    {
                        equivalent = true;
                    }
                    else if (allowResidues2ToBeSubsetOfResidues1 && matchCount > 0)
                    {
                        equivalent = true;
                    }
                }
            }

            return equivalent;
        }

        /// <summary>
        /// Returns True if this modification can affect the peptide or protein terminus
        /// </summary>
        /// <remarks>Note that some modifications can affect either peptide termini or internal residues</remarks>
        public bool CanAffectPeptideOrProteinTerminus()
        {
            var terminalSymbols = GetTerminalSymbols();

            if (ModificationType is ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod or ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod)
            {
                return true;
            }

            foreach (var chChar in mTargetResidues)
            {
                if (terminalSymbols.Contains(chChar))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// Returns true if this modification can affect peptide residues
        /// </summary>
        public bool CanAffectPeptideResidues()
        {
            var terminalSymbols = GetTerminalSymbols();

            if (ModificationType is ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod or ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod)
            {
                return false;
            }

            if (string.IsNullOrEmpty(mTargetResidues))
            {
                return true;
            }

            foreach (var chChar in mTargetResidues)
            {
                if (!terminalSymbols.Contains(chChar))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Retrieve the protein and peptide terminus symbols
        /// </summary>
        public static SortedSet<char> GetTerminalSymbols()
        {
            var terminalSymbols = new SortedSet<char>
            {
                AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS,
                AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS,
                AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS,
                AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS
            };

            return terminalSymbols;
        }

        /// <summary>
        /// Retrieve the modification type for the given modification type symbol
        /// </summary>
        /// <param name="modificationTypeSymbol">D, S, T, I, or P</param>
        public static ModificationDefinition.ResidueModificationType ModificationSymbolToModificationType(char modificationTypeSymbol)
        {
            if (modificationTypeSymbol == default(char))
            {
                return ModificationDefinition.ResidueModificationType.UnknownType;
            }

            return modificationTypeSymbol switch
            {
                'D' => ModificationDefinition.ResidueModificationType.DynamicMod,
                'S' => ModificationDefinition.ResidueModificationType.StaticMod,
                'T' => ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod,
                'I' => ModificationDefinition.ResidueModificationType.IsotopicMod,
                'P' => ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod,
                _ => ModificationDefinition.ResidueModificationType.UnknownType
            };
        }

        /// <summary>
        /// Retrieve the modification type symbol for the given modification Type
        /// </summary>
        /// <param name="modificationType"></param>
        /// <returns>D, S, T, I, or P</returns>
        public static char ModificationTypeToModificationSymbol(ModificationDefinition.ResidueModificationType modificationType)
        {
            return modificationType switch
            {
                ResidueModificationType.DynamicMod => 'D',
                ResidueModificationType.StaticMod => 'S',
                ResidueModificationType.TerminalPeptideStaticMod => 'T',
                ResidueModificationType.IsotopicMod => 'I',
                ResidueModificationType.ProteinTerminusStaticMod => 'P',
                _ => '?'
            };
        }

        /// <summary>
        /// Check whether the target residues contain the given residue
        /// </summary>
        /// <param name="chComparisonResidue"></param>
        /// <returns>True if successful, false if an error</returns>
        public bool TargetResiduesContain(char chComparisonResidue)
        {
            if (chComparisonResidue == default(char))
            {
                return false;
            }

            return TargetResidues.IndexOf(chComparisonResidue) >= 0;
        }

        /// <summary>
        /// Description of this modification definition
        /// </summary>
        public override string ToString()
        {
            return string.Format("{0} {1}, {2:F4}; {3}", ModificationType, MassCorrectionTag, ModificationMass, TargetResidues);
        }
    }
}
