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
using System.ComponentModel;
using System.Text.RegularExpressions;

namespace PHRPReader
{
    // ReSharper disable CommentTypo

    /// <summary>
    /// This class will compute the cleavage state and terminus state of a given peptide sequence.
    /// It can also be used to remove modification symbols from a sequence using ExtractCleanSequenceFromSequenceWithMods
    /// </summary>
    /// <remarks>
    /// <para>
    /// The sequence can simply contain single-letter amino acid symbols (capital letters) or a mix
    /// of amino acid symbols and modification symbols, for example:
    ///   R.PEPTIDEK.L
    ///   R.P*EPTIDEK.L
    ///   R.PEPTIDEK*.L
    /// </para>
    /// <para>
    /// Method ComputeCleavageState is overloaded to either except the peptide sequence with
    /// prefix and suffix letters (e.g. R.PEPTIDEK.L) or accept the primary peptide sequence,
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
    // ReSharper restore CommentTypo
    public class PeptideCleavageStateCalculator
    {
        // Ignore Spelling: A-Za-z, Arg, chymotrypsin, Glu, Lys, tryptic, udt

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

        /// <summary>
        /// Regular expression for the residue to the left of cleavage points for the given enzyme
        /// </summary>
        public const string TRYPSIN_LEFT_RESIDUE_REGEX = "[KR]";

        /// <summary>
        /// Regular expression for the residue to the right of cleavage points for the given enzyme
        /// </summary>
        public const string TRYPSIN_RIGHT_RESIDUE_REGEX = "[^P]";

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

        // ReSharper disable IdentifierTypo

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
            /// <summary>
            /// Trypsin: cleave after K or R, unless followed by P
            /// </summary>
            Trypsin = 0,

            /// <summary>
            /// Trypsin: cleave after any K or R
            /// </summary>
            TrypsinWithoutProlineRule = 1,

            /// <summary>
            /// Cleave after K, R, F, Y, V, E, or L
            /// </summary>
            TrypsinPlusFVLEY = 2,

            /// <summary>
            /// Chymotrypsin: cleave after F, W, Y, or L
            /// </summary>
            Chymotrypsin = 3,

            /// <summary>
            /// Both trypsin and chymotrypsin: cleave after F, W, Y, L, K, or R
            /// </summary>
            ChymotrypsinAndTrypsin = 4,

            /// <summary>
            /// Glu-C: cleave after E or D
            /// </summary>
            GluC = 5,

            /// <summary>
            /// Cyanogen bromide: cleave after M
            /// </summary>
            CyanBr = 6,

            /// <summary>
            /// Arg-C: cleave after R
            /// </summary>
            EndoArgC = 7,

            /// <summary>
            /// Lys-C: cleave after K
            /// </summary>
            EndoLysC = 8,

            /// <summary>
            /// Asp-N: cleave before D
            /// </summary>
            EndoAspN = 9
        }

        // ReSharper restore IdentifierTypo

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
        // ReSharper disable once UnusedMember.Global
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

            if (checkForPrefixAndSuffixResidues && SplitPrefixAndSuffixFromSequence(sequenceWithMods, out var primarySequence, out _, out _))
            {
                return RegexNotLetter.Replace(primarySequence, string.Empty);
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
            return new EnzymeMatchSpecInfo
            {
                LeftResidueRegEx = TRYPSIN_LEFT_RESIDUE_REGEX,
                RightResidueRegEx = TRYPSIN_RIGHT_RESIDUE_REGEX
            };
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
            // ReSharper disable StringLiteralTypo

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

                case StandardCleavageAgent.GluC:
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
                    throw new InvalidEnumArgumentException(nameof(standardCleavageAgent));
            }

            // ReSharper restore StringLiteralTypo
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

            if (periodLoc1 < 0)
            {
                return false;
            }

            var periodLoc2 = sequenceIn.LastIndexOf('.');

            if (periodLoc2 > periodLoc1 + 1)
            {
                // ReSharper disable CommentTypo

                // Sequence contains two periods with letters between the periods,
                // For example, R.PEPTIDEK.L or RPEP.TIDESEQK.L
                // Extract out the text between the periods

                primarySequence = sequenceIn.Substring(periodLoc1 + 1, periodLoc2 - periodLoc1 - 1);
                if (periodLoc1 > 0)
                {
                    prefix = sequenceIn.Substring(0, periodLoc1);
                }
                suffix = sequenceIn.Substring(periodLoc2 + 1);

                // ReSharper restore CommentTypo

                return true;
            }

            if (periodLoc2 == periodLoc1 + 1)
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

                    return true;
                }

                // Leave the sequence unchanged
                primarySequence = sequenceIn;
                return false;
            }

            if (periodLoc1 != periodLoc2)
            {
                return false;
            }

            // Peptide only contains one period
            if (periodLoc1 == 0)
            {
                primarySequence = sequenceIn.Substring(1);
                return true;
            }

            if (periodLoc1 == sequenceIn.Length - 1)
            {
                primarySequence = sequenceIn.Substring(0, periodLoc1);
                return true;
            }

            if (periodLoc1 == 1 && sequenceIn.Length > 2)
            {
                primarySequence = sequenceIn.Substring(periodLoc1 + 1);
                prefix = sequenceIn.Substring(0, periodLoc1);
                return true;
            }

            if (periodLoc1 == sequenceIn.Length - 2)
            {
                primarySequence = sequenceIn.Substring(0, periodLoc1);
                suffix = sequenceIn.Substring(periodLoc1 + 1);
                return true;
            }

            // Leave the sequence unchanged
            primarySequence = sequenceIn;

            return false;
        }

        /// <summary>
        /// Examines the two residues to see if they represent an expected cleavage point
        /// </summary>
        /// <param name="leftChar"></param>
        /// <param name="rightChar"></param>
        /// <returns>True if the characters match the currently defined cleavage rule</returns>
        public bool TestCleavageRule(char leftChar, char rightChar)
        {
            if (mUsingStandardTrypsinRules)
            {
                return leftChar is 'K' or 'R' && rightChar != 'P';
            }

            return mLeftRegEx.Match(leftChar.ToString()).Success && mRightRegEx.Match(rightChar.ToString()).Success;
        }
    }
}
