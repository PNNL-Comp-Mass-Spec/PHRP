// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 4, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/ or https://panomics.pnnl.gov/' -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
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

        private const string TRYPSIN_LEFT_RESIDUE_REGEX = @"[KR]";
        private const string TRYPSIN_RIGHT_RESIDUE_REGEX = @"[^P]";

        /// <summary>
        /// Peptide cleavage state
        /// </summary>
        public enum ePeptideCleavageStateConstants
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
        public enum ePeptideTerminusStateConstants
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
        public enum eStandardCleavageAgentConstants
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

        #region "Classwide Variables"
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
        /// <remarks></remarks>
        public clsPeptideCleavageStateCalculator()
        {
            mTerminusSymbols = new SortedSet<char> {
                TERMINUS_SYMBOL_SEQUEST,
                TERMINUS_SYMBOL_XTANDEM_NTerminus,
                TERMINUS_SYMBOL_XTANDEM_CTerminus };

            SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants.Trypsin);
        }

        /// <summary>
        /// Converts Cleavage State to 0, 1, or 2
        /// </summary>
        /// <param name="eCleavageState"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static short CleavageStateToShort(ePeptideCleavageStateConstants eCleavageState)
        {
            return Convert.ToInt16(eCleavageState);
        }

        /// <summary>
        /// Determines the cleavage state of the specified peptide
        /// </summary>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        /// <returns></returns>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public ePeptideCleavageStateConstants ComputeCleavageState(string sequenceWithPrefixAndSuffix)
        {
            if (SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out var prefix, out var suffix))
            {
                return ComputeCleavageState(primarySequence, prefix, suffix);
            }

            return ePeptideCleavageStateConstants.NonSpecific;
        }

        /// <summary>
        /// Determine the cleavage state of cleanSequence utilizing the rules specified in mEnzymeMatchSpec
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="prefixResidues"></param>
        /// <param name="suffixResidues"></param>
        /// <returns></returns>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        public ePeptideCleavageStateConstants ComputeCleavageState(string cleanSequence, string prefixResidues, string suffixResidues)
        {
            if (string.IsNullOrEmpty(cleanSequence))
                return ePeptideCleavageStateConstants.NonSpecific;

            // Find the letter closest to the end of prefixResidues
            var prefix = FindLetterNearestEnd(prefixResidues);

            // Find the letter closest to the start of suffixResidues
            var suffix = FindLetterNearestStart(suffixResidues);

            // Find the letter closest to the start of cleanSequence
            var chSequenceStart = FindLetterNearestStart(cleanSequence);

            // Find the letter closest to the end of cleanSequence
            var chSequenceEnd = FindLetterNearestEnd(cleanSequence);

            // Determine the terminus state of this peptide
            var ePeptideTerminusState = ComputeTerminusState(prefix, suffix);

            ePeptideCleavageStateConstants ePeptideCleavageState;

            if (ePeptideTerminusState == ePeptideTerminusStateConstants.ProteinNandCCTerminus)
            {
                // The peptide spans the entire length of the protein; mark it as fully tryptic
                ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
            }
            else if (ePeptideTerminusState == ePeptideTerminusStateConstants.ProteinNTerminus)
            {
                // Peptides at the N-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                if (TestCleavageRule(chSequenceEnd, suffix))
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                }
                else
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific;
                }
            }
            else if (ePeptideTerminusState == ePeptideTerminusStateConstants.ProteinCTerminus)
            {
                // Peptides at the C-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                if (TestCleavageRule(prefix, chSequenceStart))
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                }
                else
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific;
                }
            }
            else
            {
                // Check whether prefix matches mLeftRegEx and chSequenceStart matches mRightRegEx
                var ruleMatchStart = TestCleavageRule(prefix, chSequenceStart);
                var ruleMatchEnd = TestCleavageRule(chSequenceEnd, suffix);

                if (ruleMatchStart && ruleMatchEnd)
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                }
                else if (ruleMatchStart || ruleMatchEnd)
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.Partial;
                }
                else
                {
                    ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific;
                }
            }

            return ePeptideCleavageState;
        }

        /// <summary>
        /// Count the number of missed cleavages in the peptide
        /// </summary>
        /// <param name="sequenceWithPrefixAndSuffix"></param>
        /// <returns></returns>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public short ComputeNumberOfMissedCleavages(string sequenceWithPrefixAndSuffix)
        {
            short numMissedCleavages = 0;

            if (!SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out _, out _))
                return numMissedCleavages;

            if (string.IsNullOrWhiteSpace(primarySequence))
                return numMissedCleavages;

            var previousLetter = "";
            for (var index = 0; index <= primarySequence.Length - 1; index++)
            {
                var chCurrent = primarySequence[index];

                if (!clsPHRPReader.IsLetterAtoZ(chCurrent))
                    continue;

                if (!string.IsNullOrEmpty(previousLetter))
                {
                    if (TestCleavageRule(previousLetter[0], chCurrent))
                    {
                        numMissedCleavages += 1;
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
        /// <returns></returns>
        /// <remarks>Peptide must have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public ePeptideTerminusStateConstants ComputeTerminusState(string sequenceWithPrefixAndSuffix)
        {
            if (SplitPrefixAndSuffixFromSequence(sequenceWithPrefixAndSuffix, out var primarySequence, out var prefix, out var suffix))
            {
                return ComputeTerminusState(primarySequence, prefix, suffix);
            }

            return ePeptideTerminusStateConstants.None;
        }

        /// <summary>
        /// Determine the terminus state given the prefix and suffix characters
        /// </summary>
        /// <param name="prefix"></param>
        /// <param name="suffix"></param>
        /// <returns></returns>
        /// <remarks>For example, if the peptide is -.PEPTIDE.G then pass prefix="-" and suffix="G"</remarks>
        public ePeptideTerminusStateConstants ComputeTerminusState(char prefix, char suffix)
        {
            ePeptideTerminusStateConstants ePeptideTerminusState;

            if (mTerminusSymbols.Contains(prefix))
            {
                // Prefix character matches a terminus symbol
                if (mTerminusSymbols.Contains(suffix))
                {
                    // The peptide spans the entire length of the protein
                    ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus;
                }
                else
                {
                    // The peptide is located at the protein's N-terminus
                    ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus;
                }
            }
            else if (mTerminusSymbols.Contains(suffix))
            {
                // Suffix character matches a terminus symbol
                // The peptide is located at the protein's C-terminus
                ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus;
            }
            else
            {
                ePeptideTerminusState = ePeptideTerminusStateConstants.None;
            }

            return ePeptideTerminusState;
        }

        /// <summary>
        /// Determine the terminus state of the peptide
        /// </summary>
        /// <param name="cleanSequence"></param>
        /// <param name="prefixResidues"></param>
        /// <param name="suffixResidues"></param>
        /// <returns></returns>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        public ePeptideTerminusStateConstants ComputeTerminusState(string cleanSequence, string prefixResidues, string suffixResidues)
        {
            // Determine the terminus state of cleanSequence

            ePeptideTerminusStateConstants ePeptideTerminusState;

            if (string.IsNullOrEmpty(cleanSequence))
            {
                ePeptideTerminusState = ePeptideTerminusStateConstants.None;
            }
            else
            {
                // Find the letter closest to the end of prefixResidues
                var prefix = FindLetterNearestEnd(prefixResidues);

                // Find the letter closest to the start of suffixResidues
                var suffix = FindLetterNearestStart(suffixResidues);

                ePeptideTerminusState = ComputeTerminusState(prefix, suffix);
            }

            return ePeptideTerminusState;
        }

        private static readonly Regex RegexNotLetter = new Regex(@"[^A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Removes all modification symbols (*, #, +, 8, etc.) from the peptide; optionally removes prefix and suffix letters
        /// </summary>
        /// <param name="sequenceWithMods"></param>
        /// <param name="checkForPrefixAndSuffixResidues"></param>
        /// <returns>Clean peptide sequence</returns>
        /// <remarks></remarks>
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
                    index -= 1;
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
                    index += 1;
                    chMatch = text[index];
                }
            }

            return chMatch;
        }

        /// <summary>
        /// Returns the default enzyme RegEx match specifications
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
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
                SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants.Trypsin);
                return;
            }

            try
            {
                mLeftRegEx = new Regex(mEnzymeMatchSpec.LeftResidueRegEx, REGEX_OPTIONS);
                mRightRegEx = new Regex(mEnzymeMatchSpec.RightResidueRegEx, REGEX_OPTIONS);

                if (mEnzymeMatchSpec.LeftResidueRegEx == TRYPSIN_LEFT_RESIDUE_REGEX &&
                    mEnzymeMatchSpec.RightResidueRegEx == TRYPSIN_RIGHT_RESIDUE_REGEX)
                {
                    mUsingStandardTrypsinRules = true;
                }
                else
                {
                    mUsingStandardTrypsinRules = false;
                }
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
        /// <remarks></remarks>
        public void SetEnzymeMatchSpec(string leftResidueRegEx, string rightResidueRegEx)
        {
            if (leftResidueRegEx != null && rightResidueRegEx != null)
            {
                if (leftResidueRegEx.Length == 0)
                    leftResidueRegEx = @"[A-Z]";
                if (rightResidueRegEx.Length == 0)
                    rightResidueRegEx = @"[A-Z]";

                if (leftResidueRegEx == GENERIC_RESIDUE_SYMBOL.ToString() || leftResidueRegEx == @"[" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    leftResidueRegEx = @"[A-Z]";
                }

                if (rightResidueRegEx == GENERIC_RESIDUE_SYMBOL.ToString() || rightResidueRegEx == @"[" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    rightResidueRegEx = @"[A-Z]";
                }

                if (leftResidueRegEx == @"[^" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    leftResidueRegEx = @"[^A-Z]";
                }

                if (rightResidueRegEx == @"[^" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    rightResidueRegEx = @"[^A-Z]";
                }

                mEnzymeMatchSpec.LeftResidueRegEx = leftResidueRegEx;
                mEnzymeMatchSpec.RightResidueRegEx = rightResidueRegEx;
            }

            InitializeRegExObjects();
        }

        /// <summary>
        /// Select a standard enzyme match rule
        /// </summary>
        /// <param name="eStandardCleavageAgent"></param>
        /// <remarks></remarks>
        public void SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants eStandardCleavageAgent)
        {
            switch (eStandardCleavageAgent)
            {
                case eStandardCleavageAgentConstants.Trypsin:
                    SetEnzymeMatchSpec(TRYPSIN_LEFT_RESIDUE_REGEX, TRYPSIN_RIGHT_RESIDUE_REGEX);
                    break;
                case eStandardCleavageAgentConstants.TrypsinWithoutProlineRule:
                    SetEnzymeMatchSpec(@"[KR]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.TrypsinPlusFVLEY:
                    SetEnzymeMatchSpec(@"[KRFYVEL]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.Chymotrypsin:
                    SetEnzymeMatchSpec(@"[FWYL]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.ChymotrypsinAndTrypsin:
                    SetEnzymeMatchSpec(@"[FWYLKR]", @"[A-Z]");

                    break;
                case eStandardCleavageAgentConstants.V8_aka_GluC:
                    SetEnzymeMatchSpec(@"[ED]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.CyanBr:
                    SetEnzymeMatchSpec(@"[M]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.EndoArgC:
                    SetEnzymeMatchSpec(@"[R]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.EndoLysC:
                    SetEnzymeMatchSpec(@"[K]", @"[A-Z]");
                    break;
                case eStandardCleavageAgentConstants.EndoAspN:
                    SetEnzymeMatchSpec(@"[A-Z]", @"[D]");
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
        /// <remarks>If more than one character is present before the first period or after the last period, then all characters are returned
        /// If the peptide starts with ".." then it is auto-changed to start with "."
        /// If the peptide ends with ".." then it is auto-changed to end with "."
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
        /// <remarks></remarks>
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
