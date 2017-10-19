// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 4, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/ or http://panomics.pnnl.gov/' -------------------------------------------------------------------------------
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
            /// <param name="strLeftResidueRegEx"></param>
            /// <param name="strRightResidueRegEx"></param>
            public udtEnzymeMatchSpecType(string strLeftResidueRegEx, string strRightResidueRegEx)
            {
                LeftResidueRegEx = strLeftResidueRegEx;
                RightResidueRegEx = strRightResidueRegEx;
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
            mTerminusSymbols = new char[3];
            mTerminusSymbols[0] = TERMINUS_SYMBOL_SEQUEST;
            mTerminusSymbols[1] = TERMINUS_SYMBOL_XTANDEM_NTerminus;
            mTerminusSymbols[2] = TERMINUS_SYMBOL_XTANDEM_CTerminus;
            Array.Sort(mTerminusSymbols);

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
        /// <param name="strSequenceWithPrefixAndSuffix"></param>
        /// <returns></returns>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public ePeptideCleavageStateConstants ComputeCleavageState(string strSequenceWithPrefixAndSuffix)
        {
            var strPrimarySequence = string.Empty;
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;

            if (SplitPrefixAndSuffixFromSequence(strSequenceWithPrefixAndSuffix, out strPrimarySequence, out strPrefix, out strSuffix))
            {
                return ComputeCleavageState(strPrimarySequence, strPrefix, strSuffix);
            }
            else
            {
                return ePeptideCleavageStateConstants.NonSpecific;
            }
        }

        /// <summary>
        /// Determines the cleavage state of the specified peptide
        /// </summary>
        /// <param name="strCleanSequence"></param>
        /// <param name="strPrefixResidues"></param>
        /// <param name="strSuffixResidues"></param>
        /// <returns></returns>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        public ePeptideCleavageStateConstants ComputeCleavageState(string strCleanSequence, string strPrefixResidues, string strSuffixResidues)
        {
            // Determine the cleavage state of strCleanSequence utilizing the rules specified in mEnzymeMatchSpec

            var chSequenceStart = default(char);
            var chSequenceEnd = default(char);
            var chPrefix = default(char);
            var chSuffix = default(char);

            var ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific;
            var ePeptideTerminusState = ePeptideTerminusStateConstants.None;
            var blnRuleMatchStart = false;
            var blnRuleMatchEnd = false;

            if (strCleanSequence != null && strCleanSequence.Length > 0)
            {
                // Find the letter closest to the end of strPrefixResidues
                chPrefix = FindLetterNearestEnd(strPrefixResidues);

                // Find the letter closest to the start of strSuffixResidues
                chSuffix = FindLetterNearestStart(strSuffixResidues);

                // Find the letter closest to the start of strCleanSequence
                chSequenceStart = FindLetterNearestStart(strCleanSequence);

                // Find the letter closest to the end of strCleanSequence
                chSequenceEnd = FindLetterNearestEnd(strCleanSequence);

                // Determine the terminus state of this peptide
                ePeptideTerminusState = ComputeTerminusState(chPrefix, chSuffix);

                if (ePeptideTerminusState == ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                {
                    // The peptide spans the entire length of the protein; mark it as fully tryptic
                    ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                }
                else if (ePeptideTerminusState == ePeptideTerminusStateConstants.ProteinNTerminus)
                {
                    // Peptides at the N-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                    if (TestCleavageRule(chSequenceEnd, chSuffix))
                    {
                        ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                    }
                    else
                    {
                        // Leave ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
                    }
                }
                else if (ePeptideTerminusState == ePeptideTerminusStateConstants.ProteinCTerminus)
                {
                    // Peptides at the C-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
                    if (TestCleavageRule(chPrefix, chSequenceStart))
                    {
                        ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                    }
                    else
                    {
                        // Leave ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
                    }
                }
                else
                {
                    // Check whether chPrefix matches mLeftRegEx and chSequenceStart matches mRightRegEx
                    blnRuleMatchStart = TestCleavageRule(chPrefix, chSequenceStart);
                    blnRuleMatchEnd = TestCleavageRule(chSequenceEnd, chSuffix);

                    if (blnRuleMatchStart && blnRuleMatchEnd)
                    {
                        ePeptideCleavageState = ePeptideCleavageStateConstants.Full;
                    }
                    else if (blnRuleMatchStart || blnRuleMatchEnd)
                    {
                        ePeptideCleavageState = ePeptideCleavageStateConstants.Partial;
                    }
                    else
                    {
                        // Leave ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
                    }
                }
            }

            return ePeptideCleavageState;
        }

        /// <summary>
        /// Count the number of missed cleavages in the peptide
        /// </summary>
        /// <param name="strSequenceWithPrefixAndSuffix"></param>
        /// <returns></returns>
        /// <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public short ComputeNumberOfMissedCleavages(string strSequenceWithPrefixAndSuffix)
        {
            var strPrimarySequence = string.Empty;
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;
            string strPreviousLetter = null;

            short intNumMissedCleavages = 0;

            if (SplitPrefixAndSuffixFromSequence(strSequenceWithPrefixAndSuffix, out strPrimarySequence, out strPrefix, out strSuffix))
            {
                if (!string.IsNullOrWhiteSpace(strPrimarySequence))
                {
                    strPreviousLetter = "";
                    for (var intIndex = 0; intIndex <= strPrimarySequence.Length - 1; intIndex++)
                    {
                        var chCurrent = strPrimarySequence[intIndex];

                        if (clsPHRPReader.IsLetterAtoZ(chCurrent))
                        {
                            if (!string.IsNullOrEmpty(strPreviousLetter))
                            {
                                if (TestCleavageRule(strPreviousLetter[0], chCurrent))
                                {
                                    intNumMissedCleavages += 1;
                                }
                            }

                            strPreviousLetter = chCurrent.ToString();
                        }
                    }
                }
            }

            return intNumMissedCleavages;
        }

        /// <summary>
        /// Determine the terminus state of the peptide
        /// </summary>
        /// <param name="strSequenceWithPrefixAndSuffix"></param>
        /// <returns></returns>
        /// <remarks>Peptide must have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
        public ePeptideTerminusStateConstants ComputeTerminusState(string strSequenceWithPrefixAndSuffix)
        {
            var strPrimarySequence = string.Empty;
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;

            if (SplitPrefixAndSuffixFromSequence(strSequenceWithPrefixAndSuffix, out strPrimarySequence, out strPrefix, out strSuffix))
            {
                return ComputeTerminusState(strPrimarySequence, strPrefix, strSuffix);
            }
            else
            {
                return ePeptideTerminusStateConstants.None;
            }
        }

        /// <summary>
        /// Determine the terminus state given the prefix and suffix characters
        /// </summary>
        /// <param name="chPrefix"></param>
        /// <param name="chSuffix"></param>
        /// <returns></returns>
        /// <remarks>For example, if the peptide is -.PEPTIDE.G then pass chPrefix="-" and chSuffix="G"</remarks>
        public ePeptideTerminusStateConstants ComputeTerminusState(char chPrefix, char chSuffix)
        {
            var ePeptideTerminusState = ePeptideTerminusStateConstants.None;

            if (Array.BinarySearch(mTerminusSymbols, chPrefix) >= 0)
            {
                // Prefix character matches a terminus symbol
                if (Array.BinarySearch(mTerminusSymbols, chSuffix) >= 0)
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
            else if (Array.BinarySearch(mTerminusSymbols, chSuffix) >= 0)
            {
                // Suffix character matches a terminus symbol
                // The peptide is located at the protein's C-terminus
                ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus;
            }
            else
            {
                // Leave ePeptideTerminusState = ePeptideTerminusStateConstants.None
            }

            return ePeptideTerminusState;
        }

        /// <summary>
        /// Determine the terminus state of the peptide
        /// </summary>
        /// <param name="strCleanSequence"></param>
        /// <param name="strPrefixResidues"></param>
        /// <param name="strSuffixResidues"></param>
        /// <returns></returns>
        /// <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
        public ePeptideTerminusStateConstants ComputeTerminusState(string strCleanSequence, string strPrefixResidues, string strSuffixResidues)
        {
            // Determine the terminus state of strCleanSequence

            var chPrefix = default(char);
            var chSuffix = default(char);
            var ePeptideTerminusState = ePeptideTerminusStateConstants.None;

            if (strCleanSequence == null || strCleanSequence.Length == 0)
            {
                ePeptideTerminusState = ePeptideTerminusStateConstants.None;
            }
            else
            {
                // Find the letter closest to the end of strPrefixResidues
                chPrefix = FindLetterNearestEnd(strPrefixResidues);

                // Find the letter closest to the start of strSuffixResidues
                chSuffix = FindLetterNearestStart(strSuffixResidues);

                ePeptideTerminusState = ComputeTerminusState(chPrefix, chSuffix);
            }

            return ePeptideTerminusState;
        }

        private static readonly Regex RegexNotLetter = new Regex(@"[^A-Za-z]", RegexOptions.Compiled);

        /// <summary>
        /// Removes all modification symbols (*, #, +, 8, etc.) from the peptide; optionally removes prefix and suffix letters
        /// </summary>
        /// <param name="strSequenceWithMods"></param>
        /// <param name="blnCheckForPrefixAndSuffixResidues"></param>
        /// <returns>Clean peptide sequence</returns>
        /// <remarks></remarks>
        public static string ExtractCleanSequenceFromSequenceWithMods(string strSequenceWithMods, bool blnCheckForPrefixAndSuffixResidues)
        {
            if (strSequenceWithMods == null)
            {
                return string.Empty;
            }

            // Use a RegEx to remove any characters that are not letters, then return the result
            // This method of string parsing is 4x faster than using a StringBuilder object

            if (blnCheckForPrefixAndSuffixResidues)
            {
                var strPrimarySequence = string.Empty;
                var strPrefix = string.Empty;
                var strSuffix = string.Empty;

                if (SplitPrefixAndSuffixFromSequence(strSequenceWithMods, out strPrimarySequence, out strPrefix, out strSuffix))
                {
                    return RegexNotLetter.Replace(strPrimarySequence, string.Empty);
                }
            }

            return RegexNotLetter.Replace(strSequenceWithMods, string.Empty);
        }

        private char FindLetterNearestEnd(string strText)
        {
            var intIndex = 0;
            var chMatch = default(char);

            if (strText == null || strText.Length == 0)
            {
                chMatch = TERMINUS_SYMBOL_SEQUEST;
            }
            else
            {
                intIndex = strText.Length - 1;
                chMatch = strText[intIndex];
                while (!(clsPHRPReader.IsLetterAtoZ(chMatch) || Array.BinarySearch(mTerminusSymbols, chMatch) >= 0) && intIndex > 0)
                {
                    intIndex -= 1;
                    chMatch = strText[intIndex];
                }
            }

            return chMatch;
        }

        private char FindLetterNearestStart(string strText)
        {
            var intIndex = 0;
            var chMatch = default(char);

            if (strText == null || strText.Length == 0)
            {
                chMatch = TERMINUS_SYMBOL_SEQUEST;
            }
            else
            {
                intIndex = 0;
                chMatch = strText[intIndex];
                while (!(clsPHRPReader.IsLetterAtoZ(chMatch) || Array.BinarySearch(mTerminusSymbols, chMatch) >= 0) && intIndex < strText.Length - 1)
                {
                    intIndex += 1;
                    chMatch = strText[intIndex];
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
        /// <param name="strLeftResidueRegEx"></param>
        /// <param name="strRightResidueRegEx"></param>
        /// <remarks></remarks>
        public void SetEnzymeMatchSpec(string strLeftResidueRegEx, string strRightResidueRegEx)
        {
            if (strLeftResidueRegEx != null && strRightResidueRegEx != null)
            {
                if (strLeftResidueRegEx.Length == 0)
                    strLeftResidueRegEx = @"[A-Z]";
                if (strRightResidueRegEx.Length == 0)
                    strRightResidueRegEx = @"[A-Z]";

                if (strLeftResidueRegEx == GENERIC_RESIDUE_SYMBOL.ToString() || strLeftResidueRegEx == @"[" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    strLeftResidueRegEx = @"[A-Z]";
                }

                if (strRightResidueRegEx == GENERIC_RESIDUE_SYMBOL.ToString() || strRightResidueRegEx == @"[" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    strRightResidueRegEx = @"[A-Z]";
                }

                if (strLeftResidueRegEx == @"[^" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    strLeftResidueRegEx = @"[^A-Z]";
                }

                if (strRightResidueRegEx == @"[^" + GENERIC_RESIDUE_SYMBOL + @"]")
                {
                    strRightResidueRegEx = @"[^A-Z]";
                }

                mEnzymeMatchSpec.LeftResidueRegEx = strLeftResidueRegEx;
                mEnzymeMatchSpec.RightResidueRegEx = strRightResidueRegEx;
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
        /// Examines strSequenceIn and splits apart into prefix, primary sequence, and suffix
        /// </summary>
        /// <param name="strSequenceIn">Peptide sequence to examine</param>
        /// <param name="strPrimarySequence">Primary sequence (output)</param>
        /// <param name="strPrefix">Prefix residue (output)</param>
        /// <param name="strSuffix">Suffix residue (output)</param>
        /// <returns> Returns True if success, False if prefix and suffix residues were not found</returns>
        /// <remarks>If more than one character is present before the first period or after the last period, then all characters are returned
        /// If the peptide starts with ".." then it is auto-changed to start with "."
        /// If the peptide ends with ".." then it is auto-changed to end with "."
        /// </remarks>
        public static bool SplitPrefixAndSuffixFromSequence(string strSequenceIn,
            out string strPrimarySequence,
            out string strPrefix,
            out string strSuffix)
        {
            var intPeriodLoc1 = 0;
            var intPeriodLoc2 = 0;
            var blnSuccess = false;

            strPrefix = string.Empty;
            strSuffix = string.Empty;
            strPrimarySequence = string.Empty;

            if (strSequenceIn == null || strSequenceIn.Length == 0)
            {
                return false;
            }
            else
            {
                if (strSequenceIn.StartsWith("..") && strSequenceIn.Length > 2)
                {
                    strSequenceIn = "." + strSequenceIn.Substring(2);
                }

                if (strSequenceIn.EndsWith("..") && strSequenceIn.Length > 2)
                {
                    strSequenceIn = strSequenceIn.Substring(0, strSequenceIn.Length - 2) + ".";
                }

                strPrimarySequence = string.Copy(strSequenceIn);

                // See if strSequenceIn contains two periods
                intPeriodLoc1 = strSequenceIn.IndexOf('.');
                if (intPeriodLoc1 >= 0)
                {
                    intPeriodLoc2 = strSequenceIn.LastIndexOf('.');

                    if (intPeriodLoc2 > intPeriodLoc1 + 1)
                    {
                        // Sequence contains two periods with letters between the periods,
                        // For example, A.BCDEFGHIJK.L or ABCD.BCDEFGHIJK.L
                        // Extract out the text between the periods
                        strPrimarySequence = strSequenceIn.Substring(intPeriodLoc1 + 1, intPeriodLoc2 - intPeriodLoc1 - 1);
                        if (intPeriodLoc1 > 0)
                        {
                            strPrefix = strSequenceIn.Substring(0, intPeriodLoc1);
                        }
                        strSuffix = strSequenceIn.Substring(intPeriodLoc2 + 1);

                        blnSuccess = true;
                    }
                    else if (intPeriodLoc2 == intPeriodLoc1 + 1)
                    {
                        // Peptide contains two periods in a row
                        if (intPeriodLoc1 <= 1)
                        {
                            strPrimarySequence = string.Empty;

                            if (intPeriodLoc1 > 0)
                            {
                                strPrefix = strSequenceIn.Substring(0, intPeriodLoc1);
                            }
                            strSuffix = strSequenceIn.Substring(intPeriodLoc2 + 1);

                            blnSuccess = true;
                        }
                        else
                        {
                            // Leave the sequence unchanged
                            strPrimarySequence = string.Copy(strSequenceIn);
                            blnSuccess = false;
                        }
                    }
                    else if (intPeriodLoc1 == intPeriodLoc2)
                    {
                        // Peptide only contains one period
                        if (intPeriodLoc1 == 0)
                        {
                            strPrimarySequence = strSequenceIn.Substring(1);
                            blnSuccess = true;
                        }
                        else if (intPeriodLoc1 == strSequenceIn.Length - 1)
                        {
                            strPrimarySequence = strSequenceIn.Substring(0, intPeriodLoc1);
                            blnSuccess = true;
                        }
                        else if (intPeriodLoc1 == 1 && strSequenceIn.Length > 2)
                        {
                            strPrimarySequence = strSequenceIn.Substring(intPeriodLoc1 + 1);
                            strPrefix = strSequenceIn.Substring(0, intPeriodLoc1);
                            blnSuccess = true;
                        }
                        else if (intPeriodLoc1 == strSequenceIn.Length - 2)
                        {
                            strPrimarySequence = strSequenceIn.Substring(0, intPeriodLoc1);
                            strSuffix = strSequenceIn.Substring(intPeriodLoc1 + 1);
                            blnSuccess = true;
                        }
                        else
                        {
                            // Leave the sequence unchanged
                            strPrimarySequence = string.Copy(strSequenceIn);
                        }
                    }
                }
            }

            return blnSuccess;
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
