// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 3, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;
using PHRPReader.Data;

// ReSharper disable UnusedMember.Global

namespace PHRPReader
{
    /// <summary>
    /// This class will compute the mass of a given peptide sequence.  The sequence
    /// must consist of only capital letters, though if RemovePrefixAndSuffixIfPresent = True,
    /// characters up to the first . and after the last . in the sequence will be removed.
    /// Residue modification information can be supplied by passing an array of modifications
    /// using the structure PeptideSequenceModInfo
    /// </summary>
    public class PeptideMassCalculator
    {
        // Ignore Spelling: Asn, Da, Gln, Glu, Ile, Leu, Selenocysteine

        /// <summary>
        /// Symbol used when the modification is not an isotopic modification
        /// </summary>
        public const char NO_AFFECTED_ATOM_SYMBOL = '-';

        /// <summary>
        /// Monoisotopic mass of hydrogen
        /// </summary>
        public const double MASS_HYDROGEN = 1.0078246;

        /// <summary>
        /// Monoisotopic mass of oxygen
        /// </summary>
        public const double MASS_OXYGEN = 15.9949141;

        /// <summary>
        /// Monoisotopic mass of a proton
        /// </summary>
        /// <remarks>This is the mass of hydrogen minus the mass of one electron</remarks>
        public const double MASS_PROTON = 1.00727649;

        /// <summary>
        /// Monoisotopic mass of an electron
        /// </summary>
        public const double MASS_ELECTRON = 0.00054811;

        /// <summary>
        /// Default N-terminal mass change (+1.007276)
        /// </summary>
        public const double DEFAULT_N_TERMINUS_MASS_CHANGE = MASS_HYDROGEN;

        /// <summary>
        /// Default C-terminal mass change (+15.9949)
        /// </summary>
        public const double DEFAULT_C_TERMINUS_MASS_CHANGE = MASS_OXYGEN + MASS_HYDROGEN;

        private const byte ASCII_VAL_LETTER_A = 65;

        /// <summary>
        /// Peptide sequence mod info
        /// </summary>
        public struct PeptideSequenceModInfo
        {
            /// <summary>
            /// Position that the modification occurs; not used by PeptideMassCalculator
            /// </summary>
            public int ResidueLocInPeptide;

            /// <summary>
            /// Modification mass
            /// </summary>
            public double ModificationMass;

            /// <summary>
            /// Affected atom
            /// </summary>
            /// <remarks>
            /// Set to Nothing or to NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications)
            /// For Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
            /// </remarks>
            public char AffectedAtom;

            /// <summary>
            /// Modification mass and residue number
            /// </summary>
            public override string ToString()
            {
                return string.Format("{0:F4} Da on residue {1}", ModificationMass, ResidueLocInPeptide);
            }
        }

        // The Amino Acid arrays contain 26 entries, corresponding to A through Z
        // Invalid/Undefined letters (J and U) have values of 0 for the mass and atom counts
        // The values can be customized using SetAminoAcidMass and SetAminoAcidAtomCounts

        private const byte AMINO_ACID_LIST_MAX_INDEX = 25;
        private double[] mAminoAcidMasses;
        private EmpiricalFormula[] mAminoAcidEmpiricalFormulas;

        // Typically mPeptideNTerminusMass + mPeptideCTerminusMass = 18.0105633 (the mass of water)

        /// <summary>
        /// Regular expression for parsing an empirical formula
        /// </summary>
        private static readonly Regex mAtomicFormulaRegEx;

        /// <summary>
        /// This dictionary tracks element symbols and monoisotopic masses
        /// </summary>
        private static readonly Dictionary<string, double> mElementMonoMasses;

        /// <summary>
        /// Charge carrier mass
        /// </summary>
        public double ChargeCarrierMass { get; set; }

        /// <summary>
        /// Most recent error message
        /// </summary>
        public string ErrorMessage { get; private set; }

        /// <summary>
        /// Peptide C-terminus mass
        /// </summary>
        public double PeptideCTerminusMass { get; set; }

        /// <summary>
        /// Peptide N-terminus mass
        /// </summary>
        public double PeptideNTerminusMass { get; set; }

        /// <summary>
        /// If true, look for and remove prefix and suffix residues
        /// </summary>
        public bool RemovePrefixAndSuffixIfPresent { get; set; }

        /// <summary>
        /// Constructor for shared (static) variables
        /// </summary>
        static PeptideMassCalculator()
        {
            mElementMonoMasses = GetElementMonoMasses();

            mAtomicFormulaRegEx = GetAtomicFormulaRegEx(mElementMonoMasses);
        }

        /// <summary>
        /// Constructor
        /// </summary>
        public PeptideMassCalculator()
        {
            ChargeCarrierMass = MASS_PROTON;
            ErrorMessage = string.Empty;
            RemovePrefixAndSuffixIfPresent = true;
            InitializeAminoAcidData();
        }

        /// <summary>
        /// Compute the monoisotopic mass of the given empirical formula
        /// </summary>
        /// <remarks>Throws an exception if an unknown symbol is encountered</remarks>
        /// <param name="empiricalFormula"></param>
        public static double ComputeMonoisotopicMass(EmpiricalFormula empiricalFormula)
        {
            double monoisotopicMass = 0;

            foreach (var element in empiricalFormula.ElementCounts)
            {
                if (mElementMonoMasses.TryGetValue(element.Key, out var elementMass))
                {
                    monoisotopicMass += element.Value * elementMass;
                }
                else
                {
                    throw new Exception("Unrecognized symbol " + element.Key);
                }
            }

            return monoisotopicMass;
        }

        /// <summary>
        /// Compute the monoisotopic mass of the compound represented by elementalComposition
        /// </summary>
        /// <param name="elementalComposition"></param>
        /// <param name="unknownSymbols"></param>
        public static double ComputeMonoisotopicMass(Dictionary<string, int> elementalComposition, out List<string> unknownSymbols)
        {
            double monoisotopicMass = 0;

            unknownSymbols = new List<string>();

            foreach (var elementItem in elementalComposition)
            {
                if (mElementMonoMasses.TryGetValue(elementItem.Key, out var elementMass))
                {
                    monoisotopicMass += elementItem.Value * elementMass;
                }
                else
                {
                    unknownSymbols.Add(elementItem.Key);
                }
            }

            return monoisotopicMass;
        }

        /// <summary>
        /// Compute the mass of the peptide sequence (it cannot contain modification symbols)
        /// </summary>
        /// <remarks>
        /// Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True
        /// If modification symbols are present, returns -1</remarks>
        /// <param name="sequence">One letter amino acid symbols (no modification symbols or numbers); can have prefix and suffix letters</param>
        /// <returns>Monoisotopic mass, or -1 if an error</returns>
        public double ComputeSequenceMass(string sequence)
        {
            var primarySequence = sequence;
            double mass = 0;
            short validResidueCount = 0;

            if (RemovePrefixAndSuffixIfPresent)
            {
                if (!PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequence, out primarySequence, out _, out _))
                {
                    // Prefix and suffix residues not present; simply copy sequence to primarySequence
                    primarySequence = sequence;
                }
            }

            if (string.IsNullOrWhiteSpace(primarySequence))
            {
                // This code should never be reached; including this as a fail-safe
                primarySequence = sequence;
            }

            ErrorMessage = string.Empty;
            foreach (var chChar in primarySequence)
            {
                // Use Convert.ToInt16 to convert to the ASCII value, then subtract 65
                var aminoAcidIndex = ConvertAminoAcidCharToIndex(chChar);

                try
                {
                    if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
                    {
                        ErrorMessage = "Unknown symbol " + chChar + " in sequence " + primarySequence;
                        validResidueCount = 0;
                        mass = -1;
                        break;
                    }

                    mass += mAminoAcidMasses[aminoAcidIndex];
                    validResidueCount++;
                }
                catch (Exception)
                {
                    // Invalid value; ignore
                }
            }

            if (validResidueCount > 0)
            {
                mass += PeptideNTerminusMass + PeptideCTerminusMass;
            }

            return mass;
        }

        /// <summary>
        /// Compute the mass of the peptide sequence; uses the information in modifiedResidues() to determine modification masses
        /// </summary>
        /// <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
        /// <param name="sequence">One letter amino acid symbols (no modification symbols or numbers)</param>
        /// <param name="modifiedResidues">List of modified residues</param>
        /// <returns>The computed mass, or -1 if an error</returns>
        public double ComputeSequenceMass(string sequence, List<PeptideSequenceModInfo> modifiedResidues)
        {
            // Note that this call to ComputeSequenceMass will reset mErrorMessage
            var mass = ComputeSequenceMass(sequence);

            if (mass < 0 || modifiedResidues == null || modifiedResidues.Count == 0)
            {
                return mass;
            }

            var empiricalFormula = new EmpiricalFormula();

            foreach (var modifiedResidue in modifiedResidues)
            {
                // Note: do not use string.IsNullOrWhiteSpace(modifiedResidue.AffectedAtom) since that does not work on a char

                if (modifiedResidue.AffectedAtom is default(char) or NO_AFFECTED_ATOM_SYMBOL)
                {
                    // Positional modification (static or dynamic mod)
                    // Simply add the modification mass to mass
                    mass += modifiedResidue.ModificationMass;
                    continue;
                }

                // Isotopic modification
                if (empiricalFormula.ElementCounts.Count == 0)
                {
                    // Initialize empiricalFormula using the amino acid sequence
                    var empiricalFormulaToAdd = ConvertAminoAcidSequenceToEmpiricalFormula(sequence);
                    empiricalFormula.AddElements(empiricalFormulaToAdd);
                }

                if (!mElementMonoMasses.ContainsKey(modifiedResidue.AffectedAtom.ToString()))
                {
                    ErrorMessage = "Unknown Affected Atom '" + modifiedResidue.AffectedAtom + "'";
                    mass = -1;
                    break;
                }

                var elementCount = empiricalFormula.GetElementCount(modifiedResidue.AffectedAtom);

                if (elementCount == 0)
                {
                    Console.WriteLine("Warning: no amino acids in {0} contain element {1}", sequence, modifiedResidue.AffectedAtom);
                }
                else
                {
                    mass += elementCount * modifiedResidue.ModificationMass;
                }
            }

            return mass;
        }

        private static readonly Regex RegexModMasses = new("[+-][0-9.]+", RegexOptions.Compiled);

        /// <summary>
        /// Compute the mass of the peptide sequence.  Supports peptide sequences with numeric mod masses
        /// Examples of numeric mods:
        ///  R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A
        ///  K.Q-17.0265QIEESTSDYDKEK.L
        /// </summary>
        /// <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
        /// <param name="sequence"></param>
        public double ComputeSequenceMassNumericMods(string sequence)
        {
            string primarySequence;

            if (RemovePrefixAndSuffixIfPresent)
            {
                if (!PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequence, out primarySequence, out _, out _))
                {
                    // Prefix and suffix residues not present; simply copy sequence to primarySequence
                    primarySequence = sequence;
                }
            }
            else
            {
                primarySequence = sequence;
            }

            var match = RegexModMasses.Match(primarySequence);

            var sequenceWithoutMods = new StringBuilder();
            var startIndex = 0;
            double modMassTotal = 0;

            while (match.Success)
            {
                if (match.Index > startIndex)
                {
                    sequenceWithoutMods.Append(primarySequence, startIndex, match.Index - startIndex);
                }

                var modMassText = match.ToString();

                if (double.TryParse(modMassText, out var modMass))
                {
                    modMassTotal += modMass;
                }

                startIndex = match.Index + modMassText.Length;
                match = match.NextMatch();
            }

            if (startIndex < primarySequence.Length)
            {
                sequenceWithoutMods.Append(primarySequence, startIndex, primarySequence.Length - startIndex);
            }

            var peptideMass = ComputeSequenceMass(sequenceWithoutMods.ToString());

            return peptideMass < 0 ? -1 : peptideMass + modMassTotal;
        }

        private short ConvertAminoAcidCharToIndex(char aminoAcidSymbol)
        {
            return (short)(Convert.ToInt16(Convert.ToByte(aminoAcidSymbol)) - ASCII_VAL_LETTER_A);
        }

        private char ConvertAminoAcidIndexToChar(byte aminoAcidIndex)
        {
            return Convert.ToChar(aminoAcidIndex + ASCII_VAL_LETTER_A);
        }

        /// <summary>
        /// Converts the m/z value from one charge state to another charge state.  Either charge state can be 0, which means an uncharged peptide
        /// </summary>
        /// <remarks>Uses the charge carrier mass defined by ChargeCarrierMass</remarks>
        /// <param name="massMZ"></param>
        /// <param name="currentCharge"></param>
        /// <param name="desiredCharge"></param>
        public double ConvoluteMass(double massMZ, int currentCharge, int desiredCharge = 1)
        {
            return ConvoluteMass(massMZ, currentCharge, desiredCharge, ChargeCarrierMass);
        }

        /// <summary>
        /// Converts the m/z value from one charge state to another charge state.  Either charge state can be 0, which means an uncharged peptide
        /// </summary>
        /// <remarks>To return the neutral mass, set desiredCharge to 0</remarks>
        /// <param name="massMZ">m/z</param>
        /// <param name="currentCharge">Current charge; if 0, assumes massMZ is the neutral, monoisotopic mass</param>
        /// <param name="desiredCharge">Desired charge</param>
        /// <param name="chargeCarrierMass">Charge carrier mass (Default is the mass of a proton)</param>
        public double ConvoluteMass(double massMZ, int currentCharge, int desiredCharge, double chargeCarrierMass)
        {
            double newMZ;

            if (Math.Abs(chargeCarrierMass) < float.Epsilon)
            {
                chargeCarrierMass = MASS_PROTON;
            }

            try
            {
                if (currentCharge == desiredCharge)
                {
                    newMZ = massMZ;
                }
                else
                {
                    if (currentCharge == 1)
                    {
                        newMZ = massMZ;
                    }
                    else if (currentCharge > 1)
                    {
                        // Convert massMZ to M+H
                        newMZ = massMZ * currentCharge - chargeCarrierMass * (currentCharge - 1);
                    }
                    else if (currentCharge == 0)
                    {
                        // Convert massMZ (which is neutral) to M+H and store in newMZ
                        newMZ = massMZ + chargeCarrierMass;
                    }
                    else
                    {
                        // Negative charges are not supported; return 0
                        return 0;
                    }

                    if (desiredCharge > 1)
                    {
                        newMZ = (newMZ + chargeCarrierMass * (desiredCharge - 1)) / desiredCharge;
                    }
                    else if (desiredCharge == 1)
                    {
                        // Return M+H, which is currently stored in newMZ
                    }
                    else if (desiredCharge == 0)
                    {
                        // Return the neutral mass
                        newMZ -= chargeCarrierMass;
                    }
                    else
                    {
                        // Negative charges are not supported; return 0
                        newMZ = 0;
                    }
                }
            }
            catch (Exception)
            {
                // Error occurred
                newMZ = 0;
            }

            return newMZ;
        }

        /// <summary>
        /// Convert an amino acid sequence into an empirical formula
        /// </summary>
        /// <param name="sequence">One letter amino acid symbols (no modification symbols or numbers)</param>
        private EmpiricalFormula ConvertAminoAcidSequenceToEmpiricalFormula(string sequence)
        {
            var empiricalFormula = new EmpiricalFormula();

            foreach (var chAminoAcidSymbol in sequence)
            {
                // Use Convert.ToInt32 to convert to the ASCII value, then subtract 65
                var aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol);

                try
                {
                    if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
                    {
                        ErrorMessage = "Unknown symbol " + chAminoAcidSymbol + " in sequence " + sequence;
                        break;
                    }

                    empiricalFormula.AddElements(mAminoAcidEmpiricalFormulas[aminoAcidIndex]);
                }
                catch (Exception)
                {
                    // Invalid value; ignore
                }
            }

            return empiricalFormula;
        }

        /// <summary>
        /// Returns the mass of the specified amino acid
        /// </summary>
        /// <param name="chAminoAcidSymbol"></param>
        public double GetAminoAcidMass(char chAminoAcidSymbol)
        {
            // Returns the mass if success, 0 if an error

            if (chAminoAcidSymbol != default(char))
            {
                var aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol);

                if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
                {
                    // Invalid Index
                    return 0;
                }

                return mAminoAcidMasses[aminoAcidIndex];
            }

            return 0;
        }

        /// <summary>
        /// Returns a List with the number of atoms of C, H, N, O, and S in the specified amino acid
        /// </summary>
        /// <param name="chAminoAcidSymbol"></param>
        public EmpiricalFormula GetAminoAcidEmpiricalFormula(char chAminoAcidSymbol)
        {
            // Returns the atom counts if success, 0 if an error

            if (chAminoAcidSymbol == default(char))
            {
                return new EmpiricalFormula();
            }

            var aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol);

            if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
            {
                // Invalid Index
                return new EmpiricalFormula();
            }

            return mAminoAcidEmpiricalFormulas[aminoAcidIndex];
        }

        /// <summary>
        /// Create a new EmpiricalFormula instance with the specified number of atoms
        /// </summary>
        /// <param name="countC"></param>
        /// <param name="countH"></param>
        /// <param name="countN"></param>
        /// <param name="countO"></param>
        /// <param name="countS"></param>
        private EmpiricalFormula GetAminoAcidEmpiricalFormula(int countC, int countH, int countN, int countO, int countS)
        {
            var empiricalFormula = new EmpiricalFormula();

            empiricalFormula.AddElement("C", countC);
            empiricalFormula.AddElement("H", countH);
            empiricalFormula.AddElement("N", countN);
            empiricalFormula.AddElement("O", countO);
            empiricalFormula.AddElement("S", countS);

            return empiricalFormula;
        }

        /// <summary>
        /// Create a RegEx for parsing an empirical formula that optionally contains element counts and optionally contains plus or minus signs
        /// Examples of supported empirical formulas:
        ///  CHNOS
        ///  C3H3NOS4
        ///  CH23NO-5S+4
        /// </summary>
        /// <param name="elementMonoMasses"></param>
        /// <returns>RegEx with named capture groups ElementSymbol and ElementCount</returns>
        private static Regex GetAtomicFormulaRegEx(Dictionary<string, double> elementMonoMasses)
        {
            const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline;

            var regExSpec = new StringBuilder();

            regExSpec.Append("(?<ElementSymbol>");

            foreach (var element in elementMonoMasses)
            {
                regExSpec.Append(element.Key).Append('|');
            }

            // Remove the trailing vertical bar
            regExSpec.Remove(regExSpec.Length - 1, 1);

            regExSpec.Append(")");

            // RegEx will be of the form: (?<ElementSymbol>H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al)(?<ElementCount>[+-]?\d*)
            return new Regex(regExSpec + @"(?<ElementCount>[+-]?\d*)", REGEX_OPTIONS);
        }

        private double GetDefaultAminoAcidMass(char aminoAcidSymbol, out EmpiricalFormula empiricalFormula)
        {
            // These monoisotopic masses come from those traditionally used in DMS
            // They were originally assembled by Gordon Anderson for use in ICR-2LS

            double monoMass;

            switch (aminoAcidSymbol)
            {
                case 'A':
                    monoMass = 71.0371100902557;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 1, 0);
                    break;
                case 'B':
                    // Use N or D (aka Asn/Asp)
                    monoMass = 114.042921543121;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(4, 6, 2, 2, 0);
                    break;
                case 'C':
                    monoMass = 103.009180784225;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 1, 1);
                    break;
                case 'D':
                    monoMass = 115.026938199997;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(4, 5, 1, 3, 0);
                    break;
                case 'E':
                    monoMass = 129.042587518692;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 7, 1, 3, 0);
                    break;
                case 'F':
                    monoMass = 147.068408727646;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(9, 9, 1, 1, 0);
                    break;
                case 'G':
                    monoMass = 57.0214607715607;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(2, 3, 1, 1, 0);
                    break;
                case 'H':
                    monoMass = 137.058904886246;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(6, 7, 3, 1, 0);
                    break;
                case 'I':
                    monoMass = 113.084058046341;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(6, 11, 1, 1, 0);
                    break;
                case 'J':
                    // Could use mass of Ile/Leu, but we're instead treating this as an invalid, zero-mass amino acid
                    monoMass = 0;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(0, 0, 0, 0, 0);
                    break;
                case 'K':
                    monoMass = 128.094955444336;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(6, 12, 2, 1, 0);
                    break;
                case 'L':
                    monoMass = 113.084058046341;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(6, 11, 1, 1, 0);
                    break;
                case 'M':
                    monoMass = 131.040479421616;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 9, 1, 1, 1);
                    break;
                case 'N':
                    monoMass = 114.042921543121;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(4, 6, 2, 2, 0);
                    break;
                case 'O':
                    monoMass = 114.079306125641;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 10, 2, 1, 0);
                    break;
                case 'P':
                    monoMass = 97.0527594089508;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 7, 1, 1, 0);
                    break;
                case 'Q':
                    monoMass = 128.058570861816;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 8, 2, 2, 0);
                    break;
                case 'R':
                    monoMass = 156.101100921631;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(6, 12, 4, 1, 0);
                    break;
                case 'S':
                    monoMass = 87.0320241451263;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 2, 0);
                    break;
                case 'T':
                    monoMass = 101.047673463821;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(4, 7, 1, 2, 0);
                    break;
                case 'U':
                    // Corresponds to Sec = Selenocysteine (C3H5NOSe)
                    monoMass = 150.95363;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 1, 0);
                    empiricalFormula.AddElement("Se", 1);
                    break;
                case 'V':
                    monoMass = 99.0684087276459;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 9, 1, 1, 0);
                    break;
                case 'W':
                    monoMass = 186.079306125641;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(11, 10, 2, 1, 0);
                    break;
                case 'X':
                    // Unknown; use mass of Ile/Leu
                    monoMass = 113.084058046341;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(6, 11, 1, 1, 0);
                    break;
                case 'Y':
                    monoMass = 163.063322782516;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(9, 9, 1, 2, 0);
                    break;
                case 'Z':
                    // use Q or E (aka Gln/Glu); note that these are 0.984 Da apart
                    monoMass = 128.058570861816;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(5, 8, 2, 2, 0);
                    break;
                default:
                    monoMass = 0;
                    empiricalFormula = GetAminoAcidEmpiricalFormula(0, 0, 0, 0, 0);
                    break;
            }

            var computedMass = ComputeMonoisotopicMass(empiricalFormula);

            if (Math.Abs(computedMass - monoMass) > 0.00001)
            {
                Console.WriteLine("Mass discrepancy for amino acid {0}. DMS uses {1:F4} but this class computed {2:F4}", aminoAcidSymbol, monoMass, computedMass);
            }

            return monoMass;
        }

        /// <summary>
        /// Return a dictionary of element symbols and element masses
        /// </summary>
        private static Dictionary<string, double> GetElementMonoMasses()
        {
            return new Dictionary<string, double> {
                {"H", MASS_HYDROGEN}, {"He", 4.0026029}, {"Li", 7.016005}, {"Be", 9.012183},
                {"B", 11.009305}, {"C", 12}, {"N", 14.003074}, {"O", MASS_OXYGEN},
                {"F", 18.9984032}, {"Ne", 19.992439}, {"Na", 22.98977}, {"Mg", 23.98505},
                {"Al", 26.981541}, {"Si", 27.976928}, {"P", 30.973763}, {"S", 31.972072},
                {"Cl", 34.968853}, {"Ar", 39.962383}, {"K", 38.963708}, {"Ca", 39.962591},
                {"Sc", 44.955914}, {"Ti", 47.947947}, {"V", 50.943963}, {"Cr", 51.94051},
                {"Mn", 54.938046}, {"Fe", 55.934939}, {"Co", 58.933198}, {"Ni", 57.935347},
                {"Cu", 62.929599}, {"Zn", 63.929145}, {"Ga", 68.925581}, {"Ge", 71.92208},
                {"As", 74.921596}, {"Se", 79.916521}, {"Br", 78.918336}, {"Kr", 83.911506},
                {"Rb", 84.9118}, {"Sr", 87.905625}, {"Y", 88.905856}, {"Zr", 89.904708},
                {"Nb", 92.906378}, {"Mo", 97.905405}, {"Tc", 98}, {"Ru", 101.90434},
                {"Rh", 102.905503}, {"Pd", 105.903475}, {"Ag", 106.905095}, {"Cd", 113.903361},
                {"In", 114.903875}, {"Sn", 119.902199}, {"Sb", 120.903824}, {"Te", 129.906229},
                {"I", 126.904477}, {"Xe", 131.904148}, {"Cs", 132.905433}, {"Ba", 137.905236},
                {"La", 138.906355}, {"Ce", 139.905442}, {"Pr", 140.907657}, {"Nd", 141.907731},
                {"Pm", 145}, {"Sm", 151.919741}, {"Eu", 152.921243}, {"Gd", 157.924111},
                {"Tb", 158.92535}, {"Dy", 163.929183}, {"Ho", 164.930332}, {"Er", 165.930305},
                {"Tm", 168.934225}, {"Yb", 173.938873}, {"Lu", 174.940785}, {"Hf", 179.946561},
                {"Ta", 180.948014}, {"W", 183.950953}, {"Re", 186.955765}, {"Os", 191.960603},
                {"Ir", 192.962942}, {"Pt", 194.964785}, {"Au", 196.96656}, {"Hg", 201.970632},
                {"Tl", 204.97441}, {"Pb", 207.976641}, {"Bi", 208.980388}, {"Po", 209},
                {"At", 210}, {"Rn", 222}, {"Fr", 223}, {"Ra", 227},
                {"Ac", 227}, {"Th", 232.038054}, {"Pa", 231}, {"U", 238.050786},
                {"Np", 237}, {"Pu", 244}, {"Am", 243}, {"Cm", 247},
                {"Bk", 247}, {"Cf", 251}, {"Es", 252}, {"Fm", 257},
                {"Md", 258}, {"No", 269}, {"Lr", 260}
            };
        }

        /// <summary>
        /// Parse the given empirical formula to return a dictionary of the elements
        /// Examples of supported empirical formulas:
        ///  CHNOS
        ///  C3H3NOS4
        ///  CH23NO-5S+4
        /// </summary>
        /// <param name="empiricalFormula"></param>
        /// <returns>EmpiricalFormula instance tracking the element symbols and counts</returns>
        public static EmpiricalFormula GetEmpiricalFormulaComponents(string empiricalFormula)
        {
            // Originally MS-GF+ only allowed for elements C, H, N, O, S, and P in a dynamic or static mod definition
            // It now allows for any element

            var matches = mAtomicFormulaRegEx.Matches(empiricalFormula);

            var empiricalFormulaInstance = new EmpiricalFormula();

            if (matches.Count > 0)
            {
                foreach (Match match in matches)
                {
                    var elementSymbol = match.Groups["ElementSymbol"].ToString();
                    var elementCountText = match.Groups["ElementCount"].ToString();

                    var elementCount = 1;

                    if (!string.IsNullOrEmpty(elementCountText) && elementCountText.Length > 0)
                    {
                        if (!int.TryParse(elementCountText, out elementCount))
                        {
                            throw new Exception("Error parsing empirical formula '" + empiricalFormula + "', number not found in " + elementCountText);
                        }
                    }

                    empiricalFormulaInstance.AddElement(elementSymbol, elementCount);
                }
            }

            return empiricalFormulaInstance;
        }

        private void InitializeAminoAcidData()
        {
            mAminoAcidMasses = new double[AMINO_ACID_LIST_MAX_INDEX + 1];
            mAminoAcidEmpiricalFormulas = new EmpiricalFormula[AMINO_ACID_LIST_MAX_INDEX + 1];

            for (byte index = 0; index <= AMINO_ACID_LIST_MAX_INDEX; index++)
            {
                UpdateAminoAcidStatEntry(index);
            }

            ResetTerminusMasses();
        }

        /// <summary>
        /// Converts massToConvert to ppm, based on the value of currentMZ
        /// </summary>
        /// <param name="massToConvert"></param>
        /// <param name="currentMZ"></param>
        public static double MassToPPM(double massToConvert, double currentMZ)
        {
            return massToConvert * 1000000.0 / currentMZ;
        }

        /// <summary>
        /// Converts and MH mass to the uncharged (neutral) mass
        /// </summary>
        /// <remarks>Equivalent to ConvoluteMass(mH, 1, 0)</remarks>
        /// <param name="mH"></param>
        public double MHToMonoisotopicMass(double mH)
        {
            return ConvoluteMass(mH, 1, 0);
        }

        /// <summary>
        /// Converts an uncharged (neutral) mass to the m/z value for the specified charge
        /// </summary>
        /// <remarks>Equivalent to ConvoluteMass(monoisotopicMass, 0, desiredCharge)</remarks>
        /// <param name="monoisotopicMass"></param>
        /// <param name="desiredCharge"></param>
        public double MonoisotopicMassToMZ(double monoisotopicMass, int desiredCharge)
        {
            return ConvoluteMass(monoisotopicMass, 0, desiredCharge);
        }

        /// <summary>
        /// Converts from a ppm value to a mass value, using the specified mass value as a reference point
        /// </summary>
        /// <param name="ppmToConvert"></param>
        /// <param name="currentMass"></param>
        public static double PPMToMass(double ppmToConvert, double currentMass)
        {
            return ppmToConvert / 1000000.0 * currentMass;
        }

        /// <summary>
        /// Reset all of the amino acid masses and atom counts to default values
        /// </summary>
        public void ResetAminoAcidMasses()
        {
            for (byte aminoAcidIndex = 0; aminoAcidIndex <= AMINO_ACID_LIST_MAX_INDEX; aminoAcidIndex++)
            {
                var aminoAcidSymbol = ConvertAminoAcidIndexToChar(aminoAcidIndex);
                ResetAminoAcidToDefault(aminoAcidSymbol);
            }
        }

        /// <summary>
        /// Reset the mass and atom counts of the given amino acid to use default values
        /// </summary>
        /// <param name="aminoAcidSymbol">Letter between A and Z</param>
        public void ResetAminoAcidToDefault(char aminoAcidSymbol)
        {
            var aminoAcidIndex = ConvertAminoAcidCharToIndex(aminoAcidSymbol);

            if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
            {
                // Invalid Index
                return;
            }

            UpdateAminoAcidStatEntry(Convert.ToByte(aminoAcidIndex));
        }

        /// <summary>
        /// Reset the N and C terminus default mass values
        /// </summary>
        public void ResetTerminusMasses()
        {
            // See comment in method InitializeAminoAcidData concerning these masses

            PeptideNTerminusMass = DEFAULT_N_TERMINUS_MASS_CHANGE;
            PeptideCTerminusMass = DEFAULT_C_TERMINUS_MASS_CHANGE;
        }

        /// <summary>
        /// Defines the number of C, H, N, O, S, etc. elements in an amino acid
        /// </summary>
        /// <param name="chAminoAcidSymbol">Amino acid symbol</param>
        /// <param name="elementalComposition">Dictionary where keys are element symbols and values are the element counts</param>
        /// <returns>True if successful, False if an invalid amino acid symbol</returns>
        public bool SetAminoAcidAtomCounts(char chAminoAcidSymbol, Dictionary<string, int> elementalComposition)
        {
            var empiricalFormula = new EmpiricalFormula(elementalComposition);
            return SetAminoAcidAtomCounts(chAminoAcidSymbol, empiricalFormula);
        }

        /// <summary>
        /// Defines the number of C, H, N, O, S, etc. elements in an amino acid
        /// </summary>
        /// <param name="chAminoAcidSymbol">>Amino acid symbol</param>
        /// <param name="empiricalFormula">Empirical formula class</param>
        /// <returns>True if successful, False if an invalid amino acid symbol</returns>
        public bool SetAminoAcidAtomCounts(char chAminoAcidSymbol, EmpiricalFormula empiricalFormula)
        {
            if (chAminoAcidSymbol == default(char))
            {
                return false;
            }

            var aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol);

            if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
            {
                // Invalid Index
                return false;
            }

            mAminoAcidEmpiricalFormulas[aminoAcidIndex] = empiricalFormula;
            return true;
        }

        /// <summary>
        /// Defines a custom mass for an amino acid
        /// </summary>
        /// <param name="chAminoAcidSymbol"></param>
        /// <param name="mass"></param>
        /// <returns>True if successful, False if an invalid amino acid symbol</returns>
        public bool SetAminoAcidMass(char chAminoAcidSymbol, double mass)
        {
            if (chAminoAcidSymbol != default(char))
            {
                var aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol);

                if (aminoAcidIndex is < 0 or > AMINO_ACID_LIST_MAX_INDEX)
                {
                    // Invalid Index
                    return false;
                }

                mAminoAcidMasses[aminoAcidIndex] = mass;
                return true;
            }

            return false;
        }

        /// <summary>
        /// Updates an entry in parallel arrays AminoAcidMasses and AminoAcidSymbols
        /// </summary>
        /// <param name="aminoAcidIndex"></param>
        private void UpdateAminoAcidStatEntry(byte aminoAcidIndex)
        {
            // Use Convert.ToChar to convert from ASCII code to the letter
            var aminoAcidSymbol = ConvertAminoAcidIndexToChar(aminoAcidIndex);

            var monoMass = GetDefaultAminoAcidMass(aminoAcidSymbol, out var empiricalFormula);

            mAminoAcidMasses[aminoAcidIndex] = monoMass;
            mAminoAcidEmpiricalFormulas[aminoAcidIndex] = empiricalFormula;
        }
    }
}
