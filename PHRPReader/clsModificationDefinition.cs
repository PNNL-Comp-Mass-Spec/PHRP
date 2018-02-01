// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
// Started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//

using System;
using System.Collections.Generic;
using System.Globalization;

namespace PHRPReader
{
    /// <summary>
    /// This class describes an amino acid modification
    /// </summary>
    public class clsModificationDefinition
    {
        #region "Constants and Enums"

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
        /// Modification types
        /// </summary>
        public enum eModificationTypeConstants
        {
            /// <summary>
            /// Unknown mod type on a residue; essentially treated as a dynamic mod
            /// </summary>
            /// <remarks></remarks>
            UnknownType = 0,

            /// <summary>
            /// Dynamic mod on a residue or peptide terminus; supported by Sequest and notated via a modification symbol; this mod is explicitly notated by X!Tandem; if a terminus mod, then the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            /// <remarks></remarks>
            DynamicMod = 1,

            /// <summary>
            /// Static mod on a residue or peptide terminus; supported by Sequest but not explicitly notated; this mod is explicitly notated by X!Tandem; if a terminus mod, then the mod symbol is associated with the first or last residue in the peptide
            /// </summary>
            /// <remarks></remarks>
            StaticMod = 2,

            /// <summary>
            /// Peptide terminus static mod (DMS Symbol is T); used by Sequest and MSGFDB; note that terminal mods are always dynamic in X!Tandem
            /// </summary>
            /// <remarks></remarks>
            TerminalPeptideStaticMod = 3,

            /// <summary>
            /// Isotopic mod, e.g. N15, or C13; supported by Sequest; most likely not supported by XTandem
            /// </summary>
            /// <remarks></remarks>
            IsotopicMod = 4,

            /// <summary>
            /// Protein terminus static mod; supported by Sequest; this mod is also supported by X!Tandem but modified residues are not explicitly notated; instead, all peptides have their mass implicitly modified by this amount
            /// </summary>
            /// <remarks></remarks>
            ProteinTerminusStaticMod = 5
        }

        #endregion

        #region "Classwide Variables"
        /// <summary>
        /// One letter symbol for this modification; use NO_SYMBOL_MODIFICATION_SYMBOL if no symbol (necessary for isotopic mods or protein terminus static mods)
        /// </summary>
        private char mModificationSymbol;

        /// <summary>
        /// Monoisotopic modification mass
        /// </summary>
        private double mModificationMass;

        /// <summary>
        /// Modification mass, stored as text
        /// </summary>
        private string mModificationMassAsText;

        /// <summary>
        /// Target residues, tracked as a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to
        /// Use the *_SYMBOL_DMS constants for the peptide and protein terminii symbols (&lt; and &gt; for the peptide terminii; [ and ] for the protein terminii)
        /// </summary>
        /// <remarks>
        /// If this is empty, the given modification can apply to any residue or terminus
        /// </remarks>
        private string mTargetResidues;

        /// <summary>
        /// Modification type
        /// </summary>
        private eModificationTypeConstants mModificationType;

        /// <summary>
        /// Name associated with the given ModificationMass; maximum length is 8 characters
        /// Cannot contain a colon, comma, or space
        /// </summary>
        private string mMassCorrectionTag;

        /// <summary>
        /// Set to Nothing or to clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications)
        /// For Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
        /// </summary>
        private char mAffectedAtom;

        /// <summary>
        /// Number of times this modification was observed in the given XML file
        /// </summary>
        private int mOccurrenceCount;

        /// <summary>
        /// True if this was an unknown mass that was auto defined
        /// </summary>
        private bool mUnknownModAutoDefined;

        #endregion

        #region "Properties"

        /// <summary>
        /// One letter symbol for this modification
        /// </summary>
        /// <value></value>
        /// <returns></returns>
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
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public double ModificationMass
        {
            get => mModificationMass;
            set
            {
                mModificationMass = value;
                mModificationMassAsText = mModificationMass.ToString(CultureInfo.InvariantCulture);
            }
        }

        /// <summary>
        /// Modification mass, stored as text
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Represents the original string value read from the data file</remarks>
        public string ModificationMassAsText
        {
            get => mModificationMassAsText;
            set => mModificationMassAsText = value;
        }

        /// <summary>
        /// Residues that this modification can apply to
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>
        /// If an empty string, then the modification can apply to any residue or terminus;
        /// Otherwise, should contain a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to.
        /// Use the *_SYMBOL_DMS constants for the peptide and protein terminii symbols
        /// (less than and greater than signs for the peptide terminii; [ and ] for the protein terminii)
        /// </remarks>
        public string TargetResidues
        {
            get => mTargetResidues;
            set
            {
                if (value == null)
                {
                    mTargetResidues = string.Empty;
                }
                else
                {
                    mTargetResidues = string.Copy(value);
                }
            }
        }

        /// <summary>
        /// Modification type
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public eModificationTypeConstants ModificationType
        {
            get => mModificationType;
            set => mModificationType = value;
        }

        /// <summary>
        /// Modification name, for example Phosph, IodoAcet, Plus1Oxy, or Methyl
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Maximum length is 8 characters; cannot contain a colon, comma, or space</remarks>
        public string MassCorrectionTag
        {
            get => mMassCorrectionTag;
            set
            {
                if (value == null)
                {
                    mMassCorrectionTag = string.Empty;
                }
                else
                {
                    mMassCorrectionTag = string.Copy(value);
                }
            }
        }

        /// <summary>
        /// Only used with Isotopic modifications, indicating the atom affected (e.g. C, H, N, O, or S)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>
        /// Set to Nothing or to clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL (a dash) for positional modifications
        /// (including terminus modifications)
        /// </remarks>
        public char AffectedAtom
        {
            get => mAffectedAtom;
            set
            {
                if (value == default(char))
                {
                    mAffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
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
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int OccurrenceCount
        {
            get => mOccurrenceCount;
            set => mOccurrenceCount = value;
        }

        /// <summary>
        /// True if this was an unknown mass that was auto defined
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool UnknownModAutoDefined
        {
            get => mUnknownModAutoDefined;
            set => mUnknownModAutoDefined = value;
        }
        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        public clsModificationDefinition()
        {
            Clear();
        }

        /// <summary>
        /// Constructor that takes a mod symbol and mod mass
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="modificationMass"></param>
        public clsModificationDefinition(char modificationSymbol, double modificationMass)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
        }

        /// <summary>
        /// Constructor that takes a mod mass, target residues, and modification tye
        /// </summary>
        /// <param name="modificationMass"></param>
        /// <param name="targetResidues"></param>
        /// <param name="eModificationType"></param>
        public clsModificationDefinition(double modificationMass, string targetResidues, eModificationTypeConstants eModificationType)
        {
            Clear();

            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = eModificationType;
        }

        /// <summary>
        /// Constructor that takes a mod symbol, mod mass, target residues, modification type, and mass correction tag
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="modificationMass"></param>
        /// <param name="targetResidues"></param>
        /// <param name="eModificationType"></param>
        /// <param name="massCorrectionTag"></param>
        public clsModificationDefinition(char modificationSymbol, double modificationMass, string targetResidues, eModificationTypeConstants eModificationType, string massCorrectionTag)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = eModificationType;
            MassCorrectionTag = massCorrectionTag;
        }

        /// <summary>
        /// Constructor that takes a mod symbol, mod mass, target residues, modification type, mass correction tag, and affected atom
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="modificationMass"></param>
        /// <param name="targetResidues"></param>
        /// <param name="eModificationType"></param>
        /// <param name="massCorrectionTag"></param>
        /// <param name="chAffectedAtom"></param>
        /// <param name="unknownModAutoDefined"></param>
        public clsModificationDefinition(char modificationSymbol, double modificationMass, string targetResidues, eModificationTypeConstants eModificationType, string massCorrectionTag, char chAffectedAtom, bool unknownModAutoDefined)
        {
            Clear();

            ModificationSymbol = modificationSymbol;
            ModificationMass = modificationMass;
            TargetResidues = targetResidues;
            ModificationType = eModificationType;
            MassCorrectionTag = massCorrectionTag;
            AffectedAtom = chAffectedAtom;
            UnknownModAutoDefined = unknownModAutoDefined;
        }

        /// <summary>
        /// Initialize the modification definition
        /// </summary>
        /// <remarks></remarks>
        public void Clear()
        {
            mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL;
            mModificationMass = 0;
            mModificationMassAsText = "0";
            mTargetResidues = string.Empty;
            mModificationType = eModificationTypeConstants.UnknownType;
            mMassCorrectionTag = INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME;
            mAffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
            mOccurrenceCount = 0;
            mUnknownModAutoDefined = false;
        }

        /// <summary>
        /// Compares b to this object, ignoring .ModificationSymbol and ignoring .AffectedResidues
        /// </summary>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAndAtom(clsModificationDefinition b)
        {
            return EquivalentMassTypeTagAndAtom(this, b);
        }

        /// <summary>
        /// Compare a to b but ignore .ModificationSymbol and .AffectedResidues
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAndAtom(clsModificationDefinition a, clsModificationDefinition b)
        {
            var equivalent =
                Math.Abs(Math.Round(a.ModificationMass - b.ModificationMass, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) - 0) < float.Epsilon &&
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
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAtomAndResidues(clsModificationDefinition b)
        {
            return EquivalentMassTypeTagAtomAndResidues(this, b);
        }

        /// <summary>
        /// Compares b to this object
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAtomAndResidues(clsModificationDefinition a, clsModificationDefinition b)
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
                    if (a.ModificationType == eModificationTypeConstants.DynamicMod ||
                        a.ModificationType == eModificationTypeConstants.StaticMod)
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
        /// <remarks></remarks>
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
                            matchCount += 1;
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
        /// Note that some modifications can affect either peptide teriminii or internal residues
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool CanAffectPeptideOrProteinTerminus()
        {
            var terminalSymbols = GetTerminalSymbols();

            if (mModificationType == eModificationTypeConstants.ProteinTerminusStaticMod || mModificationType == eModificationTypeConstants.TerminalPeptideStaticMod)
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
        /// <returns></returns>
        /// <remarks></remarks>
        public bool CanAffectPeptideResidues()
        {
            var terminalSymbols = GetTerminalSymbols();

            if (mModificationType == eModificationTypeConstants.ProteinTerminusStaticMod || mModificationType == eModificationTypeConstants.TerminalPeptideStaticMod)
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
        /// <returns></returns>
        /// <remarks></remarks>
        public static SortedSet<char> GetTerminalSymbols()
        {
            var terminalSymbols = new SortedSet<char>
            {
                clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS,
                clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS,
                clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS,
                clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS
            };

            return terminalSymbols;
        }

        /// <summary>
        /// Retrieve the modification type for the given modification type symbol
        /// </summary>
        /// <param name="modificationTypeSymbol">D, S, T, I, or P</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static eModificationTypeConstants ModificationSymbolToModificationType(char modificationTypeSymbol)
        {
            if (modificationTypeSymbol == default(char))
            {
                return eModificationTypeConstants.UnknownType;
            }

            switch (modificationTypeSymbol)
            {
                case 'D':
                    return eModificationTypeConstants.DynamicMod;
                case 'S':
                    return eModificationTypeConstants.StaticMod;
                case 'T':
                    return eModificationTypeConstants.TerminalPeptideStaticMod;
                case 'I':
                    return eModificationTypeConstants.IsotopicMod;
                case 'P':
                    return eModificationTypeConstants.ProteinTerminusStaticMod;
                default:
                    return eModificationTypeConstants.UnknownType;
            }
        }

        /// <summary>
        /// Retrieve the modification type symbol for the given modification Type
        /// </summary>
        /// <param name="eModificationType"></param>
        /// <returns>D, S, T, I, or P</returns>
        /// <remarks></remarks>
        public static char ModificationTypeToModificationSymbol(eModificationTypeConstants eModificationType)
        {
            switch (eModificationType)
            {
                case eModificationTypeConstants.DynamicMod:
                    return 'D';
                case eModificationTypeConstants.StaticMod:
                    return 'S';
                case eModificationTypeConstants.TerminalPeptideStaticMod:
                    return 'T';
                case eModificationTypeConstants.IsotopicMod:
                    return 'I';
                case eModificationTypeConstants.ProteinTerminusStaticMod:
                    return 'P';
                default:
                    return '?';
            }
        }

        /// <summary>
        /// Check whether the target residues contain the given residue
        /// </summary>
        /// <param name="chComparisonResidue"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool TargetResiduesContain(char chComparisonResidue)
        {
            if (chComparisonResidue == default(char))
            {
                return false;
            }

            if (TargetResidues.IndexOf(chComparisonResidue) >= 0)
            {
                return true;
            }

            return false;
        }
    }
}
