// This class describes an amino acid modification
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
// Started January 6, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.
using System;
using System.Collections.Generic;

namespace PHRPReader
{
    public class clsModificationDefinition
    {
        #region "Constants and Enums"
        public const char LAST_RESORT_MODIFICATION_SYMBOL = '_';                        // The underscore is used if all of the DEFAULT_MODIFICATION_SYMBOLS are used up
        public const char NO_SYMBOL_MODIFICATION_SYMBOL = '-';
        public const string UNKNOWN_MOD_BASE_NAME = "UnkMod";
        public const string INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME = UNKNOWN_MOD_BASE_NAME + "00";

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
        private char mModificationSymbol;              // One letter symbol for this modification; use NO_SYMBOL_MODIFICATION_SYMBOL if no symbol (necessary for isotopic mods or protein terminus static mods)
        private double mModificationMass;              // Monoisotopic modification mass
        private string mModificationMassAsText;        // Modification mass, stored as text
        private string mTargetResidues;                // If this string is empty, then the given modification can apply to any residue or terminus; Otherwise, should contain a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to; Use the *_SYMBOL_DMS constants for the peptide and protein terminii symbols (< and > for the peptide terminii; [ and ] for the protein terminii)
        private eModificationTypeConstants mModificationType;
        private string mMassCorrectionTag;             // Name associated with the given ModificationMass; maximum length is 8 characters; cannot contain a colon, comma, or space
        private char mAffectedAtom;                    // Set to Nothing or to clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications); for Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
        private int mOccurrenceCount;                  // Number of times this modification was observed in the given XML file
        private bool mUnknownModAutoDefined;           // True if this was an unknown mass that was auto defined
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
                mModificationMassAsText = mModificationMass.ToString();
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

        public clsModificationDefinition()
        {
            this.Clear();
        }

        public clsModificationDefinition(char chModificationSymbol, double dblModificationMass)
        {
            this.Clear();

            this.ModificationSymbol = chModificationSymbol;
            this.ModificationMass = dblModificationMass;
        }

        public clsModificationDefinition(double dblModificationMass, string strTargetResidues, eModificationTypeConstants eModificationType)
        {
            this.Clear();

            this.ModificationMass = dblModificationMass;
            this.TargetResidues = strTargetResidues;
            this.ModificationType = eModificationType;
        }

        public clsModificationDefinition(char chModificationSymbol, double dblModificationMass, string strTargetResidues, eModificationTypeConstants eModificationType, string strMassCorrectionTag)
        {
            this.Clear();

            this.ModificationSymbol = chModificationSymbol;
            this.ModificationMass = dblModificationMass;
            this.TargetResidues = strTargetResidues;
            this.ModificationType = eModificationType;
            this.MassCorrectionTag = strMassCorrectionTag;
        }

        public clsModificationDefinition(char chModificationSymbol, double dblModificationMass, string strTargetResidues, eModificationTypeConstants eModificationType, string strMassCorrectionTag, char chAffectedAtom, bool blnUnknownModAutoDefined)
        {
            this.Clear();

            this.ModificationSymbol = chModificationSymbol;
            this.ModificationMass = dblModificationMass;
            this.TargetResidues = strTargetResidues;
            this.ModificationType = eModificationType;
            this.MassCorrectionTag = strMassCorrectionTag;
            this.AffectedAtom = chAffectedAtom;
            this.UnknownModAutoDefined = blnUnknownModAutoDefined;
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
        /// Compares objB to this object, ignoring .ModificationSymbol and ignoring .AffectedResidues
        /// </summary>
        /// <param name="objB"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAndAtom(clsModificationDefinition objB)
        {
            return EquivalentMassTypeTagAndAtom(this, objB);
        }

        /// <summary>
        /// Compare objA to objB but ignore .ModificationSymbol and .AffectedResidues
        /// </summary>
        /// <param name="objA"></param>
        /// <param name="objB"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAndAtom(clsModificationDefinition objA, clsModificationDefinition objB)
        {
            var blnEquivalent = false;

            blnEquivalent = false;
            if (Math.Abs(Math.Round(objA.ModificationMass - objB.ModificationMass, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) - 0) < float.Epsilon &&
                objA.ModificationType == objB.ModificationType &&
                objA.MassCorrectionTag == objB.MassCorrectionTag &&
                objA.AffectedAtom == objB.AffectedAtom)
            {
                blnEquivalent = true;
            }

            return blnEquivalent;
        }

        /// <summary>
        /// Compares objB to this object, ignoring .ModificationSymbol
        /// </summary>
        /// <param name="objB"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAtomAndResidues(clsModificationDefinition objB)
        {
            return EquivalentMassTypeTagAtomAndResidues(this, objB);
        }

        /// <summary>
        /// Compares objB to this object
        /// </summary>
        /// <param name="objA"></param>
        /// <param name="objB"></param>
        /// <returns>True if the items are equivalent</returns>
        /// <remarks></remarks>
        public bool EquivalentMassTypeTagAtomAndResidues(clsModificationDefinition objA, clsModificationDefinition objB)
        {
            var blnEquivalent = false;

            // First compare objA to objB but ignore .ModificationSymbol and .AffectedResidues
            blnEquivalent = EquivalentMassTypeTagAndAtom(objA, objB);

            if (blnEquivalent)
            {
                // Mass, ModificationType, MassCorrectionTag, and AffectedAtom are identical
                // What about the residues?

                blnEquivalent = false;
                if (objA.TargetResidues == null && objB.TargetResidues == null)
                {
                    blnEquivalent = true;
                }
                else if (objA.TargetResidues != null && objB.TargetResidues != null)
                {
                    if (objA.ModificationType == eModificationTypeConstants.DynamicMod ||
                        objA.ModificationType == eModificationTypeConstants.StaticMod)
                    {
                        // Matching dynamic or static modification definitions
                        // Make sure each of the residues in objB.TargetResidues is present in .TargetResidues
                        if (EquivalentTargetResidues(objA.TargetResidues, objB.TargetResidues, false))
                        {
                            blnEquivalent = true;
                        }
                    }
                    else
                    {
                        // Not a dynamic or static mod; require identical target residue lists in order to flag them as identical
                        if (objA.TargetResidues == objB.TargetResidues)
                        {
                            blnEquivalent = true;
                        }
                    }
                }
            }

            return blnEquivalent;
        }

        /// <summary>
        /// Compare the residue lists (ignoring order)
        /// </summary>
        /// <param name="strResidues1"></param>
        /// <param name="strResidues2"></param>
        /// <param name="blnAllowResidues2ToBeSubsetOfResidues1"></param>
        /// <returns>True if they contain the same residues</returns>
        /// <remarks></remarks>
        public static bool EquivalentTargetResidues(string strResidues1, string strResidues2, bool blnAllowResidues2ToBeSubsetOfResidues1)
        {
            var intMatchCount = 0;
            var blnEquivalent = false;

            blnEquivalent = false;

            if (strResidues1 == null && strResidues2 == null)
            {
                // Both residues lists are blank
                blnEquivalent = true;
            }
            else if (!(strResidues1 == null || strResidues2 == null))
            {
                if (strResidues1.Length == 0 && strResidues2.Length == 0)
                {
                    blnEquivalent = true;
                }
                else if (strResidues1 == strResidues2)
                {
                    blnEquivalent = true;
                }
                else if (strResidues1.Length >= strResidues2.Length)
                {
                    // See if each of the residues in strResidues2 is in strResidues1
                    intMatchCount = 0;
                    foreach (var chChar in strResidues2)
                    {
                        if (strResidues1.IndexOf(chChar) >= 0)
                        {
                            intMatchCount += 1;
                        }
                        else
                        {
                            break;
                        }
                    }

                    if (intMatchCount == strResidues1.Length)
                    {
                        blnEquivalent = true;
                    }
                    else if (blnAllowResidues2ToBeSubsetOfResidues1 && intMatchCount > 0)
                    {
                        blnEquivalent = true;
                    }
                }
            }

            return blnEquivalent;
        }

        /// <summary>
        /// Returns True if this modification can affect the peptide or protein terminus
        /// Note that some modifications can affect either peptide teriminii or internal residues
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool CanAffectPeptideOrProteinTerminus()
        {
            var lstTerminalSymbols = GetTerminalSymbols();

            if (mModificationType == eModificationTypeConstants.ProteinTerminusStaticMod || mModificationType == eModificationTypeConstants.TerminalPeptideStaticMod)
            {
                return true;
            }
            else
            {
                foreach (var chChar in mTargetResidues)
                {
                    if (lstTerminalSymbols.Contains(chChar))
                    {
                        return true;
                    }
                }
                return false;
            }
        }

        /// <summary>
        /// Returns true if this modification can affect peptide residues
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool CanAffectPeptideResidues()
        {
            var lstTerminalSymbols = GetTerminalSymbols();

            if (mModificationType == eModificationTypeConstants.ProteinTerminusStaticMod || mModificationType == eModificationTypeConstants.TerminalPeptideStaticMod)
            {
                return false;
            }
            else
            {
                if (string.IsNullOrEmpty(mTargetResidues))
                {
                    return true;
                }
                else
                {
                    foreach (var chChar in mTargetResidues)
                    {
                        if (!lstTerminalSymbols.Contains(chChar))
                        {
                            return true;
                        }
                    }
                }

                return false;
            }
        }

        /// <summary>
        /// Retrieve the protein and peptide terminus symbols
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public static SortedSet<char> GetTerminalSymbols()
        {
            var lstTerminalSymbols = default(SortedSet<char>);

            lstTerminalSymbols = new SortedSet<char>();

            lstTerminalSymbols.Add(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS);
            lstTerminalSymbols.Add(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
            lstTerminalSymbols.Add(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS);
            lstTerminalSymbols.Add(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);

            return lstTerminalSymbols;
        }

        /// <summary>
        /// Retrieve the modification type for the given modification type symbol
        /// </summary>
        /// <param name="chModificationTypeSymbol">D, S, T, I, or P</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static clsModificationDefinition.eModificationTypeConstants ModificationSymbolToModificationType(char chModificationTypeSymbol)
        {
            if (chModificationTypeSymbol == default(char))
            {
                return clsModificationDefinition.eModificationTypeConstants.UnknownType;
            }
            else
            {
                switch (chModificationTypeSymbol)
                {
                    case 'D':
                        return clsModificationDefinition.eModificationTypeConstants.DynamicMod;
                    case 'S':
                        return clsModificationDefinition.eModificationTypeConstants.StaticMod;
                    case 'T':
                        return clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                    case 'I':
                        return clsModificationDefinition.eModificationTypeConstants.IsotopicMod;
                    case 'P':
                        return clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod;
                    default:
                        return clsModificationDefinition.eModificationTypeConstants.UnknownType;
                }
            }
        }

        /// <summary>
        /// Retrieve the modification type symbol for the given modification Type
        /// </summary>
        /// <param name="eModificationType"></param>
        /// <returns>D, S, T, I, or P</returns>
        /// <remarks></remarks>
        public static char ModificationTypeToModificationSymbol(clsModificationDefinition.eModificationTypeConstants eModificationType)
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
            else if (this.TargetResidues.IndexOf(chComparisonResidue) >= 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
}
