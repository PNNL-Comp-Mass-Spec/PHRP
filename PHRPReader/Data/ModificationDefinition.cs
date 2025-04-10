// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://www.pnnl.gov/integrative-omics
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
        // Ignore Spelling: UnkMod, Phosph, IodoAcet, Plus1Oxy, Quant, Uni

        /// <summary>
        /// Modification symbol used after all the DEFAULT_MODIFICATION_SYMBOLS have been used
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
        /// For Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
        /// Otherwise, set to PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications)
        /// </summary>
        private char mAffectedAtom;

        /// <summary>
        /// Monoisotopic modification mass
        /// </summary>
        private double mModificationMass;

        /// <summary>
        /// One letter symbol for this modification; use NO_SYMBOL_MODIFICATION_SYMBOL if no symbol (necessary for isotopic mods or protein terminus static mods)
        /// </summary>
        private char mModificationSymbol;

        /// <summary>
        /// Only used with Isotopic modifications, indicating the atom affected (e.g. C, H, N, O, or S)
        /// </summary>
        /// <remarks>
        /// Set to PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL (a dash) for positional modifications
        /// (including terminus modifications)
        /// </remarks>
        public char AffectedAtom
        {
            get => mAffectedAtom;
            set
            {
                if (value == '\0')
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
        /// Modification name, for example Phosph, IodoAcet, Plus1Oxy, or Methyl
        /// </summary>
        /// <remarks>
        /// <para>
        /// Cannot contain a colon, comma, or space
        /// </para>
        /// <para>
        /// Prior to 2020, was limited to 8 characters in length, but that restriction has been removed
        /// </para>
        /// </remarks>
        public string MassCorrectionTag { get; set; }

        /// <summary>
        /// MaxQuant modification title
        /// </summary>
        public string MaxQuantModName { get; set; }

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
                if (value == '\0')
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
        /// Modification type
        /// </summary>
        public ResidueModificationType ModificationType { get; set; }

        /// <summary>
        /// Number of times this modification was observed in the given dataset
        /// </summary>
        public int OccurrenceCount { get; set; }

        /// <summary>
        /// Residues that this modification can apply to, tracked as a
        /// space-free, comma-free list of one letter amino acid residue symbols
        /// </summary>
        /// <remarks>
        /// <para>
        /// Use the *_SYMBOL_DMS constants for the peptide and protein termini symbols
        /// (less than and greater than signs for the peptide termini; [ and ] for the protein termini)
        /// </para>
        /// <para>
        /// If an empty string, the modification can apply to any residue or terminus
        /// </para>
        /// </remarks>
        public string TargetResidues { get; set; }

        /// <summary>
        /// True if this was an unknown mass that was auto defined
        /// </summary>
        public bool UnknownModAutoDefined { get; set; }

        /// <summary>
        /// UniMod name
        /// </summary>
        /// <remarks>
        /// To see standard mod names, login as a guest to http://www.unimod.org/modifications_list.php?
        /// </remarks>
        public string UniModName { get; set; }

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
        /// <param name="modificationSymbol">Modification symbol</param>
        /// <param name="modificationMass">Modification mass</param>
        public ModificationDefinition(char modificationSymbol, double modificationMass)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
        }

        /// <summary>
        /// Constructor that takes a mod mass, target residues, and modification type
        /// </summary>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="targetResidues">Target residues</param>
        /// <param name="modificationType">Modification type</param>
        public ModificationDefinition(double modificationMass, string targetResidues, ResidueModificationType modificationType)
        {
            Clear();

            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = modificationType;
        }

        /// <summary>
        /// Constructor that takes a mod symbol, mod mass, target residues, modification type, and mass correction tag
        /// </summary>
        /// <param name="modificationSymbol">Modification symbol</param>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="targetResidues">Target residues</param>
        /// <param name="modificationType">Modification type</param>
        /// <param name="massCorrectionTag">Mass correction tag</param>
        public ModificationDefinition(char modificationSymbol, double modificationMass, string targetResidues, ResidueModificationType modificationType, string massCorrectionTag)
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
        /// <param name="modificationSymbol">Modification symbol</param>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="targetResidues">Target residues</param>
        /// <param name="modificationType">Modification type</param>
        /// <param name="massCorrectionTag">Mass correction tag</param>
        /// <param name="chAffectedAtom">Affected atom</param>
        /// <param name="unknownModAutoDefined">True if this is an auto-defined unknown modification</param>
        public ModificationDefinition(char modificationSymbol, double modificationMass, string targetResidues, ResidueModificationType modificationType, string massCorrectionTag, char chAffectedAtom, bool unknownModAutoDefined)
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
            AffectedAtom = PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
            MassCorrectionTag = INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME;
            MaxQuantModName = string.Empty;
            ModificationMass = 0;
            ModificationMassAsText = "0";
            ModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL;
            ModificationType = ResidueModificationType.UnknownType;
            OccurrenceCount = 0;
            TargetResidues = string.Empty;
            UnknownModAutoDefined = false;
            UniModName = string.Empty;
        }

        /// <summary>
        /// Compares b to this object, ignoring .ModificationSymbol and ignoring .AffectedResidues
        /// </summary>
        /// <param name="b">Modification definition to compare</param>
        /// <returns>True if the items are equivalent</returns>
        public bool EquivalentMassTypeTagAndAtom(ModificationDefinition b)
        {
            return EquivalentMassTypeTagAndAtom(this, b);
        }

        /// <summary>
        /// Compare a to b but ignore .ModificationSymbol and .AffectedResidues
        /// </summary>
        /// <param name="a">First modification definition</param>
        /// <param name="b">Second modification definition</param>
        /// <returns>True if the items are equivalent</returns>
        public bool EquivalentMassTypeTagAndAtom(ModificationDefinition a, ModificationDefinition b)
        {
            return
                Math.Abs(Math.Round(a.ModificationMass - b.ModificationMass, PeptideModificationContainer.MASS_DIGITS_OF_PRECISION) - 0) < float.Epsilon &&
                a.ModificationType == b.ModificationType &&
                a.MassCorrectionTag == b.MassCorrectionTag &&
                a.AffectedAtom == b.AffectedAtom;
        }

        /// <summary>
        /// Compares b to this object, ignoring .ModificationSymbol
        /// </summary>
        /// <param name="b">Modification definition to compare</param>
        /// <returns>True if the items are equivalent</returns>
        // ReSharper disable once UnusedMember.Global
        public bool EquivalentMassTypeTagAtomAndResidues(ModificationDefinition b)
        {
            return EquivalentMassTypeTagAtomAndResidues(this, b);
        }

        /// <summary>
        /// Compares b to this object
        /// </summary>
        /// <param name="a">First modification definition</param>
        /// <param name="b">Second modification definition</param>
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
                    if (a.ModificationType is ResidueModificationType.DynamicMod or ResidueModificationType.StaticMod)
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
        /// <param name="residues1">First list of residues</param>
        /// <param name="residues2">Second list of residues</param>
        /// <param name="allowResidues2ToBeSubsetOfResidues1">True if residues2 is allowed to be a subset of residues1</param>
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
        // ReSharper disable once UnusedMember.Global
        public bool CanAffectPeptideOrProteinTerminus()
        {
            var terminalSymbols = GetTerminalSymbols();

            if (ModificationType is ResidueModificationType.ProteinTerminusStaticMod or ResidueModificationType.TerminalPeptideStaticMod)
            {
                return true;
            }

            foreach (var chChar in TargetResidues)
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
        // ReSharper disable once UnusedMember.Global
        public bool CanAffectPeptideResidues()
        {
            var terminalSymbols = GetTerminalSymbols();

            if (ModificationType is ResidueModificationType.ProteinTerminusStaticMod or ResidueModificationType.TerminalPeptideStaticMod)
            {
                return false;
            }

            if (string.IsNullOrEmpty(TargetResidues))
            {
                return true;
            }

            foreach (var chChar in TargetResidues)
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
            return new SortedSet<char>
            {
                AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS,
                AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS,
                AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS,
                AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS
            };
        }

        /// <summary>
        /// Retrieve the modification type for the given modification type symbol
        /// </summary>
        /// <param name="modificationTypeSymbol">D, S, T, I, or P</param>
        public static ResidueModificationType ModificationSymbolToModificationType(char modificationTypeSymbol)
        {
            if (modificationTypeSymbol == '\0')
            {
                return ResidueModificationType.UnknownType;
            }

            return modificationTypeSymbol switch
            {
                'D' => ResidueModificationType.DynamicMod,
                'S' => ResidueModificationType.StaticMod,
                'T' => ResidueModificationType.TerminalPeptideStaticMod,
                'I' => ResidueModificationType.IsotopicMod,
                'P' => ResidueModificationType.ProteinTerminusStaticMod,
                _ => ResidueModificationType.UnknownType
            };
        }

        /// <summary>
        /// Retrieve the modification type symbol for the given modification Type
        /// </summary>
        /// <param name="modificationType">Modification type</param>
        /// <returns>D, S, T, I, or P</returns>
        public static char ModificationTypeToModificationSymbol(ResidueModificationType modificationType)
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
        /// <param name="chComparisonResidue">Comparison residue</param>
        /// <returns>True if successful, false if an error</returns>
        public bool TargetResiduesContain(char chComparisonResidue)
        {
            if (chComparisonResidue == '\0')
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
