// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 4, 2006
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
using System.Text.RegularExpressions;

namespace PHRPReader
{
    /// <summary>
    /// This class will compute the cleavage state and terminus state of a given peptide sequence.
    /// It can also be used to remove modification symbols from a sequence using ExtractCleanSequenceFromSequenceWithMods
    /// </summary>
    /// <remarks>
    /// The sequence can simply contain single-letter amino acid symbols (capital letters) or a mix
    /// of amino acid symbols and modification symbols, for example:
    ///   A.BCDEFGHIJK.L
    ///   A.B*CDEFGHIJK.L
    ///   A.BCDEFGHIJK*.L
    ///   A.BCDEFGHIJK.L
    ///
    /// Function ComputeCleavageState is overloaded to either except the peptide sequence with
    /// prefix and suffix letters (e.g. A.BCDEFGHIJK.L) or accept the primary peptide sequence,
    /// the prefix residue(s), and the suffix residue(s).
    ///
    /// Use EnzymeMatchSpec to specify the residues to match for cleavage
    ///
    /// The default cleavage specification is for trypsin: [KR]|[^P]
    ///
    /// Note: Function SplitPrefixAndSuffixFromSequence will change peptides that look like:
    ///      E.TGMLTQKFARSLGMLAVDNQARV..   to   E.TGMLTQKFARSLGMLAVDNQARV.
    ///   or ..TGMLTQKFARSLGMLAVDNQARV.R   to   .TGMLTQKFARSLGMLAVDNQARV.R
    /// </remarks>
    public class clsPeptideCleavageStateCalculator
    {
        #region "Constants and Enums"

        /// <summary>
        /// Generic residue symbol
        /// </summary>
        public const char GENERIC_RESIDUE_SYMBOL = 'X';

        /// <summary>
        /// Peptide terminus symbol for SEQUEST
        /// </summary>
        public const char TERMINUS_SYMBOL_SEQUEST = '-';

        /// <summary>
        /// Peptide N-terminus symbol for X!Tandem
        /// </summary>
        public const char TERMINUS_SYMBOL_XTANDEM_NTerminus = '[';

        /// <summary>
        /// /// Peptide C-terminus symbol for X!Tandem
        /// </summary>
        public const char TERMINUS_SYMBOL_XTANDEM_CTerminus = ']';

        private const string TRYPSIN_LEFT_RESIDUE_REGEX = "[KR]";
        private const string TRYPSIN_RIGHT_RESIDUE_REGEX = "[^P]";

        /// <summary>
        /// Peptide cleavage state
        /// </summary>
        public enum PeptideCleavageStateConstants
        {
            /// <summary>
            /// Unknown cleavage specificity
            /// </summary>
            Unknown = -1,

            /// <summary>
            /// E.g., non-tryptic
            /// </summary>
            NonSpecific = 0,

            /// <summary>
            /// E.g., partially tryptic
            /// </summary>
            Partial = 1,

            /// <summary>
            /// E.g., fully tryptic
            /// </summary>
            Full = 2
        }

        /// <summary>
        /// Peptide terminus state
        /// </summary>
        public enum PeptideTerminusStateConstants
        {
            /// <summary>
            /// The peptide is located in the middle of the protein
            /// </summary>
            None = 0,

            /// <summary>
            /// The peptide is located at the protein's N-terminus
            /// </summary>
            ProteinNTerminus = 1,

            /// <summary>
            /// The peptide is located at the protein's C-terminus
            /// </summary>
            ProteinCTerminus = 2,

            /// <summary>
            /// The peptide spans the entire length of the protein
            /// </summary>
            ProteinNandCCTerminus = 3
        }

        /// <summary>
        /// Standard enzymes
        /// </summary>
        public enum StandardCleavageAgentConstants
        {
#pragma warning disable 1591
            Trypsin = 0,
            TrypsinWithoutProlineRule = 1,
            TrypsinPlusFVLEY = 2,
            Chymotrypsin = 3,
            ChymotrypsinAndTrypsin = 4,
            V8_aka_GluC = 5,
            CyanBr = 6,
            EndoArgC = 7,
            EndoLysC = 8,
            EndoAspN = 9,
            V8 = 10
#pragma warning restore 1591
        }
        #endregion

        #region "Structures"

        /// <summary>
        /// Example RegEx match strings for udtEnzymeMatchSpecType:
        /// [KR] means to match K or R
        /// [^P] means the residue cannot be P
        /// [A-Z] means to match anything; empty string also means match anything
        /// </summary>
        /// <remarks>Note, this class will automatically change [X] to [A-Z] (provided GENERIC_RESIDUE_SYMBOL = "X")</remarks>
        public struct udtEnzymeMatchSpecType
        {
            /// <summary>
            /// RegEx match string for matching the residue to the left of the cleavage point
            /// </summary>
            public string LeftResidueRegEx;

            /// <summary>
            /// RegEx match string for matching the residue to the right of the cleavage point
            /// </summary>
            public string RightResidueRegEx;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="leftResidueRegEx"></param>
            /// <param name="rightResidueRegEx"></param>
            public udtEnzymeMatchSpecType(string leftResidueRegEx, string rightResidueRegEx)
            {
                LeftResidueRegEx = leftResidueRegEx;
                RightResidueRegEx = rightResidueRegEx;
            }
        }

        #endregion

        #region "Class wide Variables"
        private udtEnzymeMatchSpecType mEnzymeMatchSpec;
        private Regex mLeftRegEx;
        private Regex mRightRegEx;
        private bool mUsingStandardTrypsinRules;

        /// <summary>
        /// This array holds TERMINUS_SYMBOL_SEQUEST, TERMINUS_SYMBOL_XTANDEM_NTerminus, and TERMINUS_SYMBOL_XTANDEM_CTerminus
        /// and is useful for quickly checking for the presence of a terminus symbol using a binary search
        /// </summary>
        private readonly SortedSet<char> mTerminusSymbols;

        #endregion

        #region "Properties"

        /// <summary>
        /// RegEx patterns for matching cleavage site residues
        /// </summary>
        public udtEnzymeMatchSpecType EnzymeMatchSpec
        {
            get => mEnzymeMatchSpec;
            set => SetEnzymeMatchSpec(value.LeftResidueRegEx, value.RightResidueRegEx);
        }

        /// <summary>
        /// Array of peptide terminus symbols
        /// </summary>
        public SortedSet<char> TerminusSymbols => mTerminusSymbols;

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        public clsPeptideCleavageStateCalculator()
        {
            mTerminusSymbols = new SortedSet<char> {
                TERMINUS_SYMBOL_SEQUEST,
                TERMINUS_SYMBOL_XTANDEM_NTerminus,
                TERMINUS_SYMBOL_XTANDEM_CTerminus };

            SetStandardEnzymeMatchSpec(StandardCleavageAgentConstants.Trypsin);
        }

        /// <summary>
        /// Converts Cleavage State to 0, 1, or 2
        /// </summary>
        /// <param name="cleavageState"></param>
        public static short CleavageStateToShort(PeptideCleavageStateConstants cleavageState)
        {
            return Convert.ToInt16(cleavageState);
        }

        /// <summary>
        /// Determines the cleavage state of the specified peptide
        /// </summary>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public PeptideCleavageStateConstants ComputeCleavageState(string sequenceWithPrefixAndSuffix)
        {
            if (SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out var prefix, out var suffix))
            {
                return ComputeCleavageState(primarySequence, prefix, suffix);
            }

            return PeptideCleavageStateConstants.NonSpecific;
        }

        /// <summary>
        /// Determine the cleavage state of cleanSequence utilizing the rules specified in mEnzymeMatchSpec
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="prefixResidues"></param>
        /// <param name="suffixResidues"></param>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        public PeptideCleavageStateConstants ComputeCleavageState(string cleanSequence, string prefixResidues, string suffixResidues)
        {
            if (string.IsNullOrEmpty(cleanSequence))
                return PeptideCleavageStateConstants.NonSpecific;

            // Find the letter closest to the end of prefixResidues
            var prefix = FindLetterNearestEnd(prefixResidues);

            // Find the letter closest to the start of suffixResidues
            var suffix = FindLetterNearestStart(suffixResidues);

            // Find the letter closest to the start of cleanSequence
            var chSequenceStart = FindLetterNearestStart(cleanSequence);

            // Find the letter closest to the end of cleanSequence
            var chSequenceEnd = FindLetterNearestEnd(cleanSequence);

            // Determine the terminus state of this peptide
            var peptideTerminusState = ComputeTerminusState(prefix, suffix);

            PeptideCleavageStateConstants peptideCleavageState;

            if (peptideTerminusState == PeptideTerminusStateConstants.ProteinNandCCTerminus)
            {
                // The peptide spans the entire length of the protein; mark it as fully tryptic
                peptideCleavageState = PeptideCleavageStateConstants.Full;
            }
            else if (peptideTerminusState == PeptideTerminusStateConstants.ProteinNTerminus)
            {
                // Peptides at the N-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                if (TestCleavageRule(chSequenceEnd, suffix))
                {
                    peptideCleavageState = PeptideCleavageStateConstants.Full;
                }
                else
                {
                    peptideCleavageState = PeptideCleavageStateConstants.NonSpecific;
                }
            }
            else if (peptideTerminusState == PeptideTerminusStateConstants.ProteinCTerminus)
            {
                // Peptides at the C-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                if (TestCleavageRule(prefix, chSequenceStart))
                {
                    peptideCleavageState = PeptideCleavageStateConstants.Full;
                }
                else
                {
                    peptideCleavageState = PeptideCleavageStateConstants.NonSpecific;
                }
            }
            else
            {
                // Check whether prefix matches mLeftRegEx and chSequenceStart matches mRightRegEx
                var ruleMatchStart = TestCleavageRule(prefix, chSequenceStart);
                var ruleMatchEnd = TestCleavageRule(chSequenceEnd, suffix);

                if (ruleMatchStart && ruleMatchEnd)
                {
                    peptideCleavageState = PeptideCleavageStateConstants.Full;
                }
                else if (ruleMatchStart || ruleMatchEnd)
                {
                    peptideCleavageState = PeptideCleavageStateConstants.Partial;
                }
                else
                {
                    peptideCleavageState = PeptideCleavageStateConstants.NonSpecific;
                }
            }

            return peptideCleavageState;
        }

        /// <summary>
        /// Count the number of missed cleavages in the peptide
        /// </summary>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public short ComputeNumberOfMissedCleavages(string sequenceWithPrefixAndSuffix)
        {
            short numMissedCleavages = 0;

            if (!SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out _, out _))
                return numMissedCleavages;

            if (string.IsNullOrWhiteSpace(primarySequence))
                return numMissedCleavages;

            var previousLetter = string.Empty;
            for (var index = 0; index <= primarySequence.Length - 1; index++)
            {
                var chCurrent = primarySequence[index];

                if (!clsPHRPReader.IsLetterAtoZ(chCurrent))
                    continue;

                if (!string.IsNullOrEmpty(previousLetter))
                {
                    if (TestCleavageRule(previousLetter[0], chCurrent))
                    {
                        numMissedCleavages++;
                    }
                }

                previousLetter = chCurrent.ToString();
            }

            return numMissedCleavages;
        }

        /// <summary>
        /// Determine the terminus state of the peptide
        /// </summary>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        /// <remarks>Peptide must have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public PeptideTerminusStateConstants ComputeTerminusState(string sequenceWithPrefixAndSuffix)
        {
            if (SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out var prefix, out var suffix))
            {
                return ComputeTerminusState(primarySequence, prefix, suffix);
            }

            return PeptideTerminusStateConstants.None;
        }

        /// <summary>
        /// Determine the terminus state given the prefix and suffix characters
        /// </summary>
        /// <param name="prefix"></param>
        /// <param name="suffix"></param>
        /// <remarks>For example, if the peptide is -.PEPTIDE.G, pass prefix="-" and suffix="G"</remarks>
        public PeptideTerminusStateConstants ComputeTerminusState(char prefix, char suffix)
        {
            PeptideTerminusStateConstants peptideTerminusState;

            if (mTerminusSymbols.Contains(prefix))
            {
                // Prefix character matches a terminus symbol
                if (mTerminusSymbols.Contains(suffix))
                {
                    // The peptide spans the entire length of the protein
                    peptideTerminusState = PeptideTerminusStateConstants.ProteinNandCCTerminus;
                }
                else
                {
                    // The peptide is located at the protein's N-terminus
                    peptideTerminusState = PeptideTerminusStateConstants.ProteinNTerminus;
                }
            }
            else if (mTerminusSymbols.Contains(suffix))
            {
                // Suffix character matches a terminus symbol
                // The peptide is located at the protein's C-terminus
                peptideTerminusState = PeptideTerminusStateConstants.ProteinCTerminus;
            }
            else
            {
                peptideTerminusState = PeptideTerminusStateConstants.None;
            }

            return peptideTerminusState;
        }

        /// <summary>
        /// Determine the terminus state of the peptide
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="prefixResidues"></param>
        /// <param name="suffixResidues"></param>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        public PeptideTerminusStateConstants ComputeTerminusState(string cleanSequence, string prefixResidues, string suffixResidues)
        {
            // Determine the terminus state of cleanSequence

            if (string.IsNullOrEmpty(cleanSequence))
            {
                return PeptideTerminusStateConstants.None;
            }

            // Find the letter closest to the end of prefixResidues
            var prefix = FindLetterNearestEnd(prefixResidues);

            // Find the letter closest to the start of suffixResidues
            var suffix = FindLetterNearestStart(suffixResidues);

            return ComputeTerminusState(prefix, suffix);
        }

        private static readonly Regex RegexNotLetter = new Regex("[^A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Removes all modification symbols (*, #, +, 8, etc.) from the peptide; optionally removes prefix and suffix letters
        /// </summary>
        /// <param name="sequenceWithMods"></param>
        /// <param name="checkForPrefixAndSuffixResidues"></param>
        /// <returns>Clean peptide sequence</returns>
        public static string ExtractCleanSequenceFromSequenceWithMods(string sequenceWithMods, bool checkForPrefixAndSuffixResidues)
        {
            if (sequenceWithMods == null)
            {
                return string.Empty;
            }

            // Use a RegEx to remove any characters that are not letters, then return the result
            // This method of string parsing is 4x faster than using a StringBuilder object

            if (checkForPrefixAndSuffixResidues)
            {
                if (SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out _, out _))
                {
                    return RegexNotLetter.Replace(primarySequence, string.Empty);
                }
            }

            return RegexNotLetter.Replace(sequenceWithMods, string.Empty);
        }

        private char FindLetterNearestEnd(string text)
        {
            char chMatch;

            if (string.IsNullOrEmpty(text))
            {
                chMatch = TERMINUS_SYMBOL_SEQUEST;
            }
            else
            {
                var index = text.Length - 1;
                chMatch = text[index];
                while (!(clsPHRPReader.IsLetterAtoZ(chMatch) || mTerminusSymbols.Contains(chMatch)) && index > 0)
                {
                    index--;
                    chMatch = text[index];
                }
            }

            return chMatch;
        }

        private char FindLetterNearestStart(string text)
        {
            char chMatch;

            if (string.IsNullOrEmpty(text))
            {
                chMatch = TERMINUS_SYMBOL_SEQUEST;
            }
            else
            {
                var index = 0;
                chMatch = text[index];
                while (!(clsPHRPReader.IsLetterAtoZ(chMatch) || mTerminusSymbols.Contains(chMatch)) && index < text.Length - 1)
                {
                    index++;
                    chMatch = text[index];
                }
            }

            return chMatch;
        }

        /// <summary>
        /// Returns the default enzyme RegEx match specifications
        /// </summary>
        public static udtEnzymeMatchSpecType GetDefaultEnzymeMatchSpec()
        {
            var udtEnzymeMatchSpec = default(udtEnzymeMatchSpecType);
            udtEnzymeMatchSpec.LeftResidueRegEx = TRYPSIN_LEFT_RESIDUE_REGEX;
            udtEnzymeMatchSpec.RightResidueRegEx = TRYPSIN_RIGHT_RESIDUE_REGEX;
            return udtEnzymeMatchSpec;
        }

        private void InitializeRegExObjects()
        {
            const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

            if (mEnzymeMatchSpec.LeftResidueRegEx == null || mEnzymeMatchSpec.RightResidueRegEx == null)
            {
                // Note that calling SetStandardEnzymeMatchSpec will cause this sub to also be called
                SetStandardEnzymeMatchSpec(StandardCleavageAgentConstants.Trypsin);
                return;
            }

            try
            {
                mLeftRegEx = new Regex(mEnzymeMatchSpec.LeftResidueRegEx, REGEX_OPTIONS);
                mRightRegEx = new Regex(mEnzymeMatchSpec.RightResidueRegEx, REGEX_OPTIONS);

                mUsingStandardTrypsinRules =
                    mEnzymeMatchSpec.LeftResidueRegEx == TRYPSIN_LEFT_RESIDUE_REGEX &&
                    mEnzymeMatchSpec.RightResidueRegEx == TRYPSIN_RIGHT_RESIDUE_REGEX;
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        /// <summary>
        /// Define custom enzyme match rules using RegEx strings
        /// </summary>
        /// <param name="leftResidueRegEx"></param>
        /// <param name="rightResidueRegEx"></param>
        public void SetEnzymeMatchSpec(string leftResidueRegEx, string rightResidueRegEx)
        {
            if (leftResidueRegEx != null && rightResidueRegEx != null)
            {
                if (leftResidueRegEx.Length == 0)
                    leftResidueRegEx = "[A-Z]";
                if (rightResidueRegEx.Length == 0)
                    rightResidueRegEx = "[A-Z]";

                if (leftResidueRegEx == GENERIC_RESIDUE_SYMBOL.ToString() || leftResidueRegEx == "[" + GENERIC_RESIDUE_SYMBOL + "]")
                {
                    leftResidueRegEx = "[A-Z]";
                }

                if (rightResidueRegEx == GENERIC_RESIDUE_SYMBOL.ToString() || rightResidueRegEx == "[" + GENERIC_RESIDUE_SYMBOL + "]")
                {
                    rightResidueRegEx = "[A-Z]";
                }

                if (leftResidueRegEx == "[^" + GENERIC_RESIDUE_SYMBOL + "]")
                {
                    leftResidueRegEx = "[^A-Z]";
                }

                if (rightResidueRegEx == "[^" + GENERIC_RESIDUE_SYMBOL + "]")
                {
                    rightResidueRegEx = "[^A-Z]";
                }

                mEnzymeMatchSpec.LeftResidueRegEx = leftResidueRegEx;
                mEnzymeMatchSpec.RightResidueRegEx = rightResidueRegEx;
            }

            InitializeRegExObjects();
        }

        /// <summary>
        /// Select a standard enzyme match rule
        /// </summary>
        /// <param name="standardCleavageAgent"></param>
        public void SetStandardEnzymeMatchSpec(StandardCleavageAgentConstants standardCleavageAgent)
        {
            switch (standardCleavageAgent)
            {
                case StandardCleavageAgentConstants.Trypsin:
                    SetEnzymeMatchSpec(TRYPSIN_LEFT_RESIDUE_REGEX, TRYPSIN_RIGHT_RESIDUE_REGEX);
                    break;
                case StandardCleavageAgentConstants.TrypsinWithoutProlineRule:
                    SetEnzymeMatchSpec("[KR]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.TrypsinPlusFVLEY:
                    SetEnzymeMatchSpec("[KRFYVEL]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.Chymotrypsin:
                    SetEnzymeMatchSpec("[FWYL]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.ChymotrypsinAndTrypsin:
                    SetEnzymeMatchSpec("[FWYLKR]", "[A-Z]");

                    break;
                case StandardCleavageAgentConstants.V8_aka_GluC:
                    SetEnzymeMatchSpec("[ED]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.CyanBr:
                    SetEnzymeMatchSpec("[M]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.EndoArgC:
                    SetEnzymeMatchSpec("[R]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.EndoLysC:
                    SetEnzymeMatchSpec("[K]", "[A-Z]");
                    break;
                case StandardCleavageAgentConstants.EndoAspN:
                    SetEnzymeMatchSpec("[A-Z]", "[D]");
                    break;
                default:
                    break;
                // Unknown agent; leave unchanged
            }
        }

        /// <summary>
        /// Examines sequenceIn and splits apart into prefix, primary sequence, and suffix
        /// </summary>
        /// <param name="sequenceIn">Peptide sequence to examine</param>
        /// <param name="primarySequence">Primary sequence (output)</param>
        /// <param name="prefix">Prefix residue (output)</param>
        /// <param name="suffix">Suffix residue (output)</param>
        /// <returns> Returns True if success, False if prefix and suffix residues were not found</returns>
        /// <remarks>If more than one character is present before the first period or after the last period, all characters are returned
        /// If the peptide starts with ".." it is auto-changed to start with "."
        /// If the peptide ends with ".." it is auto-changed to end with "."
        /// </remarks>
        public static bool SplitPrefixAndSuffixFromSequence(string sequenceIn,
            out string primarySequence,
            out string prefix,
            out string suffix)
        {
            var success = false;

            prefix = string.Empty;
            suffix = string.Empty;
            primarySequence = string.Empty;

            if (string.IsNullOrEmpty(sequenceIn))
            {
                return false;
            }

            if (sequenceIn.StartsWith("..") && sequenceIn.Length > 2)
            {
                sequenceIn = "." + sequenceIn.Substring(2);
            }

            if (sequenceIn.EndsWith("..") && sequenceIn.Length > 2)
            {
                sequenceIn = sequenceIn.Substring(0, sequenceIn.Length - 2) + ".";
            }

            primarySequence = string.Copy(sequenceIn);

            // See if sequenceIn contains two periods
            var periodLoc1 = sequenceIn.IndexOf('.');
            if (periodLoc1 >= 0)
            {
                var periodLoc2 = sequenceIn.LastIndexOf('.');

                if (periodLoc2 > periodLoc1 + 1)
                {
                    // Sequence contains two periods with letters between the periods,
                    // For example, A.BCDEFGHIJK.L or ABCD.BCDEFGHIJK.L
                    // Extract out the text between the periods
                    primarySequence = sequenceIn.Substring(periodLoc1 + 1, periodLoc2 - periodLoc1 - 1);
                    if (periodLoc1 > 0)
                    {
                        prefix = sequenceIn.Substring(0, periodLoc1);
                    }
                    suffix = sequenceIn.Substring(periodLoc2 + 1);

                    success = true;
                }
                else if (periodLoc2 == periodLoc1 + 1)
                {
                    // Peptide contains two periods in a row
                    if (periodLoc1 <= 1)
                    {
                        primarySequence = string.Empty;

                        if (periodLoc1 > 0)
                        {
                            prefix = sequenceIn.Substring(0, periodLoc1);
                        }
                        suffix = sequenceIn.Substring(periodLoc2 + 1);

                        success = true;
                    }
                    else
                    {
                        // Leave the sequence unchanged
                        primarySequence = string.Copy(sequenceIn);
                        success = false;
                    }
                }
                else if (periodLoc1 == periodLoc2)
                {
                    // Peptide only contains one period
                    if (periodLoc1 == 0)
                    {
                        primarySequence = sequenceIn.Substring(1);
                        success = true;
                    }
                    else if (periodLoc1 == sequenceIn.Length - 1)
                    {
                        primarySequence = sequenceIn.Substring(0, periodLoc1);
                        success = true;
                    }
                    else if (periodLoc1 == 1 && sequenceIn.Length > 2)
                    {
                        primarySequence = sequenceIn.Substring(periodLoc1 + 1);
                        prefix = sequenceIn.Substring(0, periodLoc1);
                        success = true;
                    }
                    else if (periodLoc1 == sequenceIn.Length - 2)
                    {
                        primarySequence = sequenceIn.Substring(0, periodLoc1);
                        suffix = sequenceIn.Substring(periodLoc1 + 1);
                        success = true;
                    }
                    else
                    {
                        // Leave the sequence unchanged
                        primarySequence = string.Copy(sequenceIn);
                    }
                }
            }

            return success;
        }

        /// <summary>
        /// Examines the two residues to see if they represent an expected cleavage point
        /// </summary>
        /// <param name="chLeftChar"></param>
        /// <param name="chRightChar"></param>
        /// <returns>True if the characters match the currently defined cleavage rule</returns>
        public bool TestCleavageRule(char chLeftChar, char chRightChar)
        {
            if (mUsingStandardTrypsinRules)
            {
                if (chLeftChar == 'K' || chLeftChar == 'R')
                {
                    if (chRightChar != 'P')
                    {
                        return true;
                    }
                }
                return false;
            }

            if (mLeftRegEx.Match(chLeftChar.ToString()).Success)
            {
                if (mRightRegEx.Match(chRightChar.ToString()).Success)
                {
                    return true;
                }
            }

            return false;
        }
    }
}
