// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 4, 2006
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
using System.Text.RegularExpressions;

namespace PHRPReader
{
    /// <summary>
    /// This class will compute the cleavage state and terminus state of a given peptide sequence.
    /// It can also be used to remove modification symbols from a sequence using ExtractCleanSequenceFromSequenceWithMods
    /// </summary>
    /// <remarks>
    /// <para>
    /// The sequence can simply contain single-letter amino acid symbols (capital letters) or a mix
    /// of amino acid symbols and modification symbols, for example:
    ///   A.BCDEFGHIJK.L
    ///   A.B*CDEFGHIJK.L
    ///   A.BCDEFGHIJK*.L
    ///   A.BCDEFGHIJK.L
    /// </para>
    /// <para>
    /// Method ComputeCleavageState is overloaded to either except the peptide sequence with
    /// prefix and suffix letters (e.g. A.BCDEFGHIJK.L) or accept the primary peptide sequence,
    /// the prefix residue(s), and the suffix residue(s).
    /// </para>
    /// <para>Use EnzymeMatchSpec to specify the residues to match for cleavage</para>
    /// <para>The default cleavage specification is for trypsin: [KR]|[^P]</para>
    /// <para>
    /// Note: Method SplitPrefixAndSuffixFromSequence will change peptides that look like:
    ///      E.TGMLTQKFARSLGMLAVDNQARV..   to   E.TGMLTQKFARSLGMLAVDNQARV.
    ///   or ..TGMLTQKFARSLGMLAVDNQARV.R   to   .TGMLTQKFARSLGMLAVDNQARV.R
    /// </para>
    /// </remarks>
    public class PeptideCleavageStateCalculator
    {
        // Ignore Spelling: A-Za-z, tryptic, udt

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
        /// Peptide C-terminus symbol for X!Tandem
        /// </summary>
        public const char TERMINUS_SYMBOL_XTANDEM_CTerminus = ']';

        private const string TRYPSIN_LEFT_RESIDUE_REGEX = "[KR]";
        private const string TRYPSIN_RIGHT_RESIDUE_REGEX = "[^P]";

        /// <summary>
        /// Peptide cleavage state
        /// </summary>
        public enum PeptideCleavageState
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
        public enum PeptideTerminusState
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
        public enum StandardCleavageAgent
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

        /// <summary>
        /// Example RegEx match strings for EnzymeMatchSpecInfo:
        /// [KR] means to match K or R
        /// [^P] means the residue cannot be P
        /// [A-Z] means to match anything; empty string also means match anything
        /// </summary>
        /// <remarks>Note, this class will automatically change [X] to [A-Z] (provided GENERIC_RESIDUE_SYMBOL = "X")</remarks>
        public struct EnzymeMatchSpecInfo
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
            public EnzymeMatchSpecInfo(string leftResidueRegEx, string rightResidueRegEx)
            {
                LeftResidueRegEx = leftResidueRegEx;
                RightResidueRegEx = rightResidueRegEx;
            }
        }

        private EnzymeMatchSpecInfo mEnzymeMatchSpec;
        private Regex mLeftRegEx;
        private Regex mRightRegEx;
        private bool mUsingStandardTrypsinRules;

        /// <summary>
        /// RegEx patterns for matching cleavage site residues
        /// </summary>
        public EnzymeMatchSpecInfo EnzymeMatchSpec
        {
            get => mEnzymeMatchSpec;
            set => SetEnzymeMatchSpec(value.LeftResidueRegEx, value.RightResidueRegEx);
        }

        /// <summary>
        /// Array of peptide terminus symbols
        /// </summary>
        /// <remarks>
        /// This holds TERMINUS_SYMBOL_SEQUEST, TERMINUS_SYMBOL_XTANDEM_NTerminus, and TERMINUS_SYMBOL_XTANDEM_CTerminus
        /// and is useful for quickly checking for the presence of a terminus symbol using a binary search
        /// </remarks>
        public SortedSet<char> TerminusSymbols { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        public PeptideCleavageStateCalculator()
        {
            TerminusSymbols = new SortedSet<char> {
                TERMINUS_SYMBOL_SEQUEST,
                TERMINUS_SYMBOL_XTANDEM_NTerminus,
                TERMINUS_SYMBOL_XTANDEM_CTerminus };

            SetStandardEnzymeMatchSpec(StandardCleavageAgent.Trypsin);
        }

        /// <summary>
        /// Converts Cleavage State to 0, 1, or 2
        /// </summary>
        /// <param name="cleavageState"></param>
        public static short CleavageStateToShort(PeptideCleavageState cleavageState)
        {
            return (short)cleavageState;
        }

        /// <summary>
        /// Determines the cleavage state of the specified peptide
        /// </summary>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        public PeptideCleavageState ComputeCleavageState(string sequenceWithPrefixAndSuffix)
        {
            if (SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out var prefix, out var suffix))
            {
                return ComputeCleavageState(primarySequence, prefix, suffix);
            }

            return PeptideCleavageState.NonSpecific;
        }

        /// <summary>
        /// Determine the cleavage state of cleanSequence utilizing the rules specified in mEnzymeMatchSpec
        /// </summary>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        /// <param name="cleanSequence"></param>
        /// <param name="prefixResidues"></param>
        /// <param name="suffixResidues"></param>
        public PeptideCleavageState ComputeCleavageState(string cleanSequence, string prefixResidues, string suffixResidues)
        {
            if (string.IsNullOrEmpty(cleanSequence))
                return PeptideCleavageState.NonSpecific;

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

            PeptideCleavageState peptideCleavageState;

            if (peptideTerminusState == PeptideTerminusState.ProteinNandCCTerminus)
            {
                // The peptide spans the entire length of the protein; mark it as fully tryptic
                peptideCleavageState = PeptideCleavageState.Full;
            }
            else if (peptideTerminusState == PeptideTerminusState.ProteinNTerminus)
            {
                // Peptides at the N-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                if (TestCleavageRule(chSequenceEnd, suffix))
                {
                    peptideCleavageState = PeptideCleavageState.Full;
                }
                else
                {
                    peptideCleavageState = PeptideCleavageState.NonSpecific;
                }
            }
            else if (peptideTerminusState == PeptideTerminusState.ProteinCTerminus)
            {
                // Peptides at the C-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                if (TestCleavageRule(prefix, chSequenceStart))
                {
                    peptideCleavageState = PeptideCleavageState.Full;
                }
                else
                {
                    peptideCleavageState = PeptideCleavageState.NonSpecific;
                }
            }
            else
            {
                // Check whether prefix matches mLeftRegEx and chSequenceStart matches mRightRegEx
                var ruleMatchStart = TestCleavageRule(prefix, chSequenceStart);
                var ruleMatchEnd = TestCleavageRule(chSequenceEnd, suffix);

                if (ruleMatchStart && ruleMatchEnd)
                {
                    peptideCleavageState = PeptideCleavageState.Full;
                }
                else if (ruleMatchStart || ruleMatchEnd)
                {
                    peptideCleavageState = PeptideCleavageState.Partial;
                }
                else
                {
                    peptideCleavageState = PeptideCleavageState.NonSpecific;
                }
            }

            return peptideCleavageState;
        }

        /// <summary>
        /// Count the number of missed cleavages in the peptide
        /// </summary>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
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

                if (!ReaderFactory.IsLetterAtoZ(chCurrent))
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
        /// <remarks>Peptide must have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        public PeptideTerminusState ComputeTerminusState(string sequenceWithPrefixAndSuffix)
        {
            if (SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out var prefix, out var suffix))
            {
                return ComputeTerminusState(primarySequence, prefix, suffix);
            }

            return PeptideTerminusState.None;
        }

        /// <summary>
        /// Determine the terminus state given the prefix and suffix characters
        /// </summary>
        /// <remarks>For example, if the peptide is -.PEPTIDE.G, pass prefix="-" and suffix="G"</remarks>
        /// <param name="prefix"></param>
        /// <param name="suffix"></param>
        public PeptideTerminusState ComputeTerminusState(char prefix, char suffix)
        {
            PeptideTerminusState peptideTerminusState;

            if (TerminusSymbols.Contains(prefix))
            {
                // Prefix character matches a terminus symbol
                if (TerminusSymbols.Contains(suffix))
                {
                    // The peptide spans the entire length of the protein
                    peptideTerminusState = PeptideTerminusState.ProteinNandCCTerminus;
                }
                else
                {
                    // The peptide is located at the protein's N-terminus
                    peptideTerminusState = PeptideTerminusState.ProteinNTerminus;
                }
            }
            else if (TerminusSymbols.Contains(suffix))
            {
                // Suffix character matches a terminus symbol
                // The peptide is located at the protein's C-terminus
                peptideTerminusState = PeptideTerminusState.ProteinCTerminus;
            }
            else
            {
                peptideTerminusState = PeptideTerminusState.None;
            }

            return peptideTerminusState;
        }

        /// <summary>
        /// Determine the terminus state of the peptide
        /// </summary>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        /// <param name="cleanSequence"></param>
        /// <param name="prefixResidues"></param>
        /// <param name="suffixResidues"></param>
        public PeptideTerminusState ComputeTerminusState(string cleanSequence, string prefixResidues, string suffixResidues)
        {
            // Determine the terminus state of cleanSequence

            if (string.IsNullOrEmpty(cleanSequence))
            {
                return PeptideTerminusState.None;
            }

            // Find the letter closest to the end of prefixResidues
            var prefix = FindLetterNearestEnd(prefixResidues);

            // Find the letter closest to the start of suffixResidues
            var suffix = FindLetterNearestStart(suffixResidues);

            return ComputeTerminusState(prefix, suffix);
        }

        private static readonly Regex RegexNotLetter = new("[^A-Za-z]", RegexOptions.Compiled);

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
                while (!(ReaderFactory.IsLetterAtoZ(chMatch) || TerminusSymbols.Contains(chMatch)) && index > 0)
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
                while (!(ReaderFactory.IsLetterAtoZ(chMatch) || TerminusSymbols.Contains(chMatch)) && index < text.Length - 1)
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
        public static EnzymeMatchSpecInfo GetDefaultEnzymeMatchSpec()
        {
            var udtEnzymeMatchSpec = default(EnzymeMatchSpecInfo);
            udtEnzymeMatchSpec.LeftResidueRegEx = TRYPSIN_LEFT_RESIDUE_REGEX;
            udtEnzymeMatchSpec.RightResidueRegEx = TRYPSIN_RIGHT_RESIDUE_REGEX;
            return udtEnzymeMatchSpec;
        }

        private void InitializeRegExObjects()
        {
            const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

            if (mEnzymeMatchSpec.LeftResidueRegEx == null || mEnzymeMatchSpec.RightResidueRegEx == null)
            {
                // Note that calling SetStandardEnzymeMatchSpec will cause this method to also be called
                SetStandardEnzymeMatchSpec(StandardCleavageAgent.Trypsin);
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
        public void SetStandardEnzymeMatchSpec(StandardCleavageAgent standardCleavageAgent)
        {
            switch (standardCleavageAgent)
            {
                case StandardCleavageAgent.Trypsin:
                    SetEnzymeMatchSpec(TRYPSIN_LEFT_RESIDUE_REGEX, TRYPSIN_RIGHT_RESIDUE_REGEX);
                    break;
                case StandardCleavageAgent.TrypsinWithoutProlineRule:
                    SetEnzymeMatchSpec("[KR]", "[A-Z]");
                    break;
                case StandardCleavageAgent.TrypsinPlusFVLEY:
                    SetEnzymeMatchSpec("[KRFYVEL]", "[A-Z]");
                    break;
                case StandardCleavageAgent.Chymotrypsin:
                    SetEnzymeMatchSpec("[FWYL]", "[A-Z]");
                    break;
                case StandardCleavageAgent.ChymotrypsinAndTrypsin:
                    SetEnzymeMatchSpec("[FWYLKR]", "[A-Z]");

                    break;
                case StandardCleavageAgent.V8_aka_GluC:
                    SetEnzymeMatchSpec("[ED]", "[A-Z]");
                    break;
                case StandardCleavageAgent.CyanBr:
                    SetEnzymeMatchSpec("[M]", "[A-Z]");
                    break;
                case StandardCleavageAgent.EndoArgC:
                    SetEnzymeMatchSpec("[R]", "[A-Z]");
                    break;
                case StandardCleavageAgent.EndoLysC:
                    SetEnzymeMatchSpec("[K]", "[A-Z]");
                    break;
                case StandardCleavageAgent.EndoAspN:
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
        /// <remarks>If more than one character is present before the first period or after the last period, all characters are returned
        /// If the peptide starts with ".." it is auto-changed to start with "."
        /// If the peptide ends with ".." it is auto-changed to end with "."
        /// </remarks>
        /// <param name="sequenceIn">Peptide sequence to examine</param>
        /// <param name="primarySequence">Output: Primary sequence</param>
        /// <param name="prefix">Output: Prefix residue</param>
        /// <param name="suffix">Output: Suffix residue</param>
        /// <returns> Returns True if success, False if prefix and suffix residues were not found</returns>
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

            primarySequence = sequenceIn;

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
                        primarySequence = sequenceIn;
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
                        primarySequence = sequenceIn;
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
                if (chLeftChar is 'K' or 'R')
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
