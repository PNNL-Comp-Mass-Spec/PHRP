// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 5, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//

using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader
{
    /// <summary>
    /// This class is used to track modifications that can be applied to peptides
    /// It handles both residue level modifications and static, peptide-wide modifications
    /// </summary>
    /// <remarks>
    /// Use ReadMassCorrectionTagsFile() and ReadModificationDefinitionsFile() to customize
    /// the default mass correction tag and modification definition lists
    /// </remarks>
    public class clsPeptideModificationContainer
    {
        #region "Constants and Enums"

        /// <summary>
        /// Default modification symbols
        /// </summary>
        public const string DEFAULT_MODIFICATION_SYMBOLS = "*#@$&!%~†‡¤º^`×÷+=ø¢";         // A few other possibilities: €£¥§

        /// <summary>
        /// Digits of precision to round masses to when finding mass correction tags by mass
        /// </summary>
        public const byte MASS_DIGITS_OF_PRECISION = 3;

        /// <summary>
        /// Symbol used by X!Tandem for tracking a modification at the peptide N-terminus
        /// </summary>
        public const char N_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM = '[';

        /// <summary>
        /// Symbol used by X!Tandem for tracking a modification at the peptide C-terminus
        /// </summary>
        public const char C_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM = ']';

        /// <summary>
        /// Symbol used by Inspect for tracking a modification at the peptide N-terminus
        /// </summary>
        public const char N_TERMINAL_PEPTIDE_MOD_SYMBOL_INSPECT = '[';


        /// <summary>
        /// Symbol used by Inspect for tracking a modification at the peptide C-terminus
        /// </summary>
        public const char C_TERMINAL_PEPTIDE_MOD_SYMBOL_INSPECT = ']';

        #endregion

        #region "Structures"
        #endregion

        #region "Classwide Variables"
        // List of available modification symbols
        private Queue mDefaultModificationSymbols;

        // List of known mass correction tags
        private readonly Dictionary<string, double> mMassCorrectionTags;

        // List of known modifications
        private readonly List<clsModificationDefinition> mModifications;

        private string mErrorMessage;

        // This array holds modifications that Sequest or XTandem will often use but for
        // which the auto-addition method sometimes incorrectly notes
        private List<clsModificationDefinition> mStandardRefinementModifications;

        private Dictionary<int, string> mIntegerMassCorrectionTagLookup;

        #endregion

        #region "Properties"

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage => mErrorMessage;

        /// <summary>
        /// Modification count
        /// </summary>
        public int ModificationCount => mModifications.Count;

        /// <summary>
        /// Modification list
        /// </summary>
        public List<clsModificationDefinition> Modifications => mModifications;

        /// <summary>
        /// When true, take the mod symbol into account when finding identical mods
        /// </summary>
        public bool ConsiderModSymbolWhenFindingIdenticalMods { get; set; }

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        public clsPeptideModificationContainer()
        {
            mMassCorrectionTags = new Dictionary<string, double>();
            mModifications = new List<clsModificationDefinition>();

            InitializeLocalVariables();
        }

        /// <summary>
        /// Add objModificationDefinition to mModifications
        /// However, do not add if a duplicate modification
        /// Furthermore, if everything matches except for .TargetResidues, then add the new target residues to the existing, matching mod
        /// </summary>
        /// <param name="objModificationDefinition"></param>
        /// <param name="blnUseNextAvailableModificationSymbol"></param>
        /// <returns>The index of the newly added modification, or the the index of the modification that objModificationDefinition matches </returns>
        /// <remarks></remarks>
        private int AddModification(clsModificationDefinition objModificationDefinition, bool blnUseNextAvailableModificationSymbol)
        {
            int intModificationIndex;

            var blnMatchFound = false;

            // See if any of the existing modifications match objModificationDefinition, ignoring .TargetResidues and possibly ignoring .ModificationSymbol
            for (intModificationIndex = 0; intModificationIndex <= mModifications.Count - 1; intModificationIndex++)
            {
                if (mModifications[intModificationIndex].EquivalentMassTypeTagAndAtom(objModificationDefinition))
                {
                    blnMatchFound = true;

                    if (ConsiderModSymbolWhenFindingIdenticalMods)
                    {
                        if (objModificationDefinition.ModificationSymbol != mModifications[intModificationIndex].ModificationSymbol)
                        {
                            // Symbols differ; add this as a new modification definition
                            blnMatchFound = false;
                        }
                    }

                    if (blnMatchFound)
                    {
                        var mod = mModifications[intModificationIndex];
                        if (mod.ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mod.ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                        {
                            // Matching dynamic or static modification definitions
                            // Merge the two modifications by making sure each of the residues in objModificationDefinition.TargetResidues is present in .TargetResidues
                            foreach (var chChar in objModificationDefinition.TargetResidues)
                            {
                                if (!mod.TargetResiduesContain(chChar))
                                {
                                    mod.TargetResidues += chChar;
                                }
                            }

                            if (!blnUseNextAvailableModificationSymbol)
                            {
                                // See if the new modification symbol is different than the already-defined symbol

                                if (objModificationDefinition.ModificationSymbol != mod.ModificationSymbol)
                                {
                                }
                            }
                        }
                    }
                }
                if (blnMatchFound)
                    break;
            }

            if (blnMatchFound)
                return intModificationIndex;

            if (blnUseNextAvailableModificationSymbol && mDefaultModificationSymbols.Count > 0)
            {
                // Add objModificationDefinition to the list, using the next available default modification symbol
                objModificationDefinition.ModificationSymbol = Convert.ToChar(mDefaultModificationSymbols.Dequeue());
            }

            mModifications.Add(objModificationDefinition);

            return mModifications.Count - 1;
        }

        private clsModificationDefinition AddUnknownModification(
            double dblModificationMass,
            clsModificationDefinition.eModificationTypeConstants eModType,
            char chTargetResidue,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool addToModificationListIfUnknown,
            bool blnUseNextAvailableModificationSymbol,
            char chModSymbol,
            byte massDigitsOfPrecision)
        {
            string strTargetResidues;

            if (chTargetResidue == default(char))
            {
                strTargetResidues = string.Empty;
            }
            else
            {
                strTargetResidues = chTargetResidue.ToString();
            }

            if (eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
            {
                // Assume this is a terminus mod
                switch (eResidueTerminusState)
                {
                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus:
                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus:
                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus:
                        strTargetResidues = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        break;
                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus:
                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus:
                        strTargetResidues = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        break;
                    default:
                        throw new Exception("Unrecognized eResidueTerminusStateConstants enum: " + eResidueTerminusState);
                }
            }

            if (!blnUseNextAvailableModificationSymbol)
            {
                chModSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
            }

            var strMassCorrectionTag = LookupMassCorrectionTagByMass(dblModificationMass, massDigitsOfPrecision, true, massDigitsOfPrecision);

            var objModificationDefinition = new clsModificationDefinition(chModSymbol, dblModificationMass, strTargetResidues, eModType, strMassCorrectionTag, clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, true);

            if (!addToModificationListIfUnknown)
                return objModificationDefinition;

            // Append objModificationDefinition to mModifications()
            int intNewModIndex;
            if (mDefaultModificationSymbols.Count > 0 && blnUseNextAvailableModificationSymbol)
            {
                intNewModIndex = AddModification(objModificationDefinition, blnUseNextAvailableModificationSymbol: true);
            }
            else
            {
                intNewModIndex = AddModification(objModificationDefinition, blnUseNextAvailableModificationSymbol: false);
            }

            if (intNewModIndex >= 0)
            {
                return mModifications[intNewModIndex];
            }

            return objModificationDefinition;

            // Either addToModificationListIfUnknown = False or no more default modification symbols
            // Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
        }

        /// <summary>
        /// Append standard refinement modifications
        /// </summary>
        public void AppendStandardRefinmentModifications()
        {
            for (var intIndex = 0; intIndex <= mStandardRefinementModifications.Count - 1; intIndex++)
            {
                var mod = mStandardRefinementModifications[intIndex];
                VerifyModificationPresent(mod.ModificationMass, mod.TargetResidues, mod.ModificationType);
            }
        }

        /// <summary>
        /// Clear modifications
        /// </summary>
        public void ClearModifications()
        {
            UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS);
            mModifications.Clear();
        }

        /// <summary>
        /// Converts a modification mass to a generic 8 character name
        /// The name will always start with + or - then will have the modification mass, rounded as necessary to give an 8 character name
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        private string GenerateGenericModMassName(double dblModificationMass)
        {
            int intFormatDigits;
            int intMaxLength;

            if (Math.Abs(dblModificationMass) < float.Epsilon)
            {
                return "+0.00000";
            }

            if (dblModificationMass < -9999999)
            {
                // Modification mass is too negative; always return -9999999
                return "-9999999";
            }

            if (dblModificationMass > 9999999)
            {
                // Modification mass is too positve; always return +9999999
                return "+9999999";
            }

            // Determine the number of digits that we will display to the left of the decimal point
            var dblFormatDigits = Math.Log10(Math.Abs(dblModificationMass));
            if (Math.Abs(dblFormatDigits - Convert.ToInt32(dblFormatDigits)) < float.Epsilon)
            {
                // ModMass is a power of 10
                intFormatDigits = Convert.ToInt32(dblFormatDigits) + 1;
            }
            else
            {
                intFormatDigits = Convert.ToInt32(Math.Ceiling(dblFormatDigits));
            }

            if (intFormatDigits < 1)
                intFormatDigits = 1;

            // Generate the format string
            // For example, a modification mass of 15.9994 will have strFormatString = "+00.0000"
            // Negative modification masses do not start with a minus sign in the format string since Visual Studio auto-adds it
            // Thus, a modification mass of -15.9994 will have strFormatString = "00.0000"
            var strFormatString = new string('0', intFormatDigits);

            if (dblModificationMass > 0)
            {
                strFormatString = "+" + strFormatString;
                intMaxLength = 8;
            }
            else
            {
                intMaxLength = 7;
            }

            if (strFormatString.Length < intMaxLength)
            {
                strFormatString += ".";
            }

            while (strFormatString.Length < intMaxLength)
            {
                strFormatString += "0";
            }

            var strModMassName = dblModificationMass.ToString(strFormatString);

            if (strModMassName.Length < 8 && strModMassName.IndexOf('.') < 0)
            {
                strModMassName += ".";
            }

            while (strModMassName.Length < 8)
            {
                // Append extra zeroes (this code will likely never be reached)
                strModMassName += "0";
            }

            if (strModMassName.Length > 8)
            {
                throw new ArgumentOutOfRangeException("Generated Mod Name is longer than 8 characters: " + strModMassName);
            }

            return strModMassName;
        }

        /// <summary>
        /// Looks for the best match in mIntegerMassCorrectionTagLookup for dblModificationMass (which should be close to a integer value)
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <returns>The mass correction tag name if a match, otherwise nothing</returns>
        /// <remarks></remarks>
        private string GetBestIntegerBasedMassCorrectionTag(double dblModificationMass)
        {
            var strClosestMassCorrectionTag = string.Empty;

            foreach (var tagOverride in mIntegerMassCorrectionTagLookup)
            {
                if (Math.Abs(dblModificationMass - tagOverride.Key) < 0.0001)
                {
                    strClosestMassCorrectionTag = tagOverride.Value;
                    break;
                }
            }

            return strClosestMassCorrectionTag;
        }

        /// <summary>
        /// Get a modification, by index
        /// </summary>
        /// <param name="intIndex"></param>
        /// <returns></returns>
        public clsModificationDefinition GetModificationByIndex(int intIndex)
        {
            if (intIndex >= 0 & intIndex < mModifications.Count)
            {
                return mModifications[intIndex];
            }

            return new clsModificationDefinition();
        }

        /// <summary>
        /// Get the modification type, by modification index
        /// </summary>
        /// <param name="intIndex"></param>
        /// <returns></returns>
        public clsModificationDefinition.eModificationTypeConstants GetModificationTypeByIndex(int intIndex)
        {
            if (intIndex >= 0 & intIndex < mModifications.Count)
            {
                return mModifications[intIndex].ModificationType;
            }

            return clsModificationDefinition.eModificationTypeConstants.UnknownType;
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not found
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <returns>Mod name, or empty string if no match</returns>
        /// <remarks>
        /// Searches known mods using 3 digits of precision, then 2 digits, then 1 digit
        /// If no match, adds as a new, unknown modification
        /// </remarks>
        public string LookupMassCorrectionTagByMass(double dblModificationMass)
        {
            const byte massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION;
            const bool addToModificationListIfUnknown = true;

            return LookupMassCorrectionTagByMass(dblModificationMass, massDigitsOfPrecision, addToModificationListIfUnknown);
        }

        /// <summary>
        ///  Find the mass correction tag with the given mass, adding to the unknown modification list if not  found
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="massDigitsOfPrecision"></param>
        /// <returns>Mod name, or empty string if no match</returns>
        /// <remarks>
        /// Searches known mods using massDigitsOfPrecision digits of precision, then massDigitsOfPrecision-1 digits, ... 1 digit
        /// If no match, adds as a new, unknown modification
        /// </remarks>
        public string LookupMassCorrectionTagByMass(double dblModificationMass, byte massDigitsOfPrecision)
        {
            const bool addToModificationListIfUnknown = true;

            return LookupMassCorrectionTagByMass(dblModificationMass, massDigitsOfPrecision, addToModificationListIfUnknown);
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not found and addToModificationListIfUnknown is true
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="massDigitsOfPrecision"></param>
        /// <param name="addToModificationListIfUnknown"></param>
        /// <returns>Mod name, or empty string if no match</returns>
        /// <remarks>
        /// Searches known mods using massDigitsOfPrecision digits of precision, then massDigitsOfPrecision-1 digits, ... 1 digit
        /// If no match, adds as a new, unknown modification
        /// </remarks>
        public string LookupMassCorrectionTagByMass(double dblModificationMass, byte massDigitsOfPrecision, bool addToModificationListIfUnknown)
        {
            const byte massDigitsOfPrecisionLoose = 1;

            return LookupMassCorrectionTagByMass(dblModificationMass, massDigitsOfPrecision, addToModificationListIfUnknown, massDigitsOfPrecisionLoose);
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not found and addToModificationListIfUnknown is true
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="massDigitsOfPrecision"></param>
        /// <param name="addToModificationListIfUnknown"></param>
        /// <param name="massDigitsOfPrecisionLoose"></param>
        /// <returns>Mod name, or empty string if no match</returns>
        public string LookupMassCorrectionTagByMass(
            double dblModificationMass,
            byte massDigitsOfPrecision,
            bool addToModificationListIfUnknown,
            byte massDigitsOfPrecisionLoose)
        {
            int intMassDigitsOfPrecisionStop;

            if (massDigitsOfPrecision < massDigitsOfPrecisionLoose)
                massDigitsOfPrecisionLoose = massDigitsOfPrecision;

            if (massDigitsOfPrecision >= massDigitsOfPrecisionLoose)
            {
                intMassDigitsOfPrecisionStop = massDigitsOfPrecisionLoose;
            }
            else
            {
                intMassDigitsOfPrecisionStop = massDigitsOfPrecision;
            }

            for (var currentPrecision = (int)massDigitsOfPrecision; currentPrecision >= intMassDigitsOfPrecisionStop; currentPrecision += -1)
            {
                var strClosestMassCorrectionTag = string.Empty;
                var dblClosestMassCorrectionTagMassDiff = double.MaxValue;

                try
                {
                    if (intMassDigitsOfPrecisionStop == 0)
                    {
                        strClosestMassCorrectionTag = GetBestIntegerBasedMassCorrectionTag(dblModificationMass);
                        if (!string.IsNullOrEmpty(strClosestMassCorrectionTag))
                        {
                            dblClosestMassCorrectionTagMassDiff = 0;
                            break;
                        }
                    }

                    // First look for an exact match in mMassCorrectionTags
                    foreach (var massCorrectionTag in mMassCorrectionTags)
                    {
                        var massDiff = Math.Abs(dblModificationMass - massCorrectionTag.Value);
                        if (massDiff < closestMassCorrectionTagMassDiff)
                        {
                            closestMassCorrectionTag = massCorrectionTag.Key;
                            closestMassCorrectionTagMassDiff = massDiff;
                        }
                    }
                }
                catch (Exception)
                {
                    // Error enumerating through mMassCorrectionTags
                }

                if (Math.Abs(Math.Round(dblClosestMassCorrectionTagMassDiff, currentPrecision)) < float.Epsilon)
                {
                    // Match found
                    return strClosestMassCorrectionTag;
                }

                if (currentPrecision > intMassDigitsOfPrecisionStop)
                {
                    // Let the For loop go through another iteration to see if we find a match
                }
                else
                {
                    // Match not found
                    // Name the modification based on the mod mass
                    strClosestMassCorrectionTag = GenerateGenericModMassName(dblModificationMass);

                    if (addToModificationListIfUnknown)
                    {
                        if (mMassCorrectionTags.ContainsKey(closestMassCorrectionTag))
                        {
                            Console.WriteLine("Warning: Ignoring duplicate mass correction tag: {0}, mass {1:F3}",
                                              closestMassCorrectionTag, dblModificationMass);
                        }
                        else
                        {
                            mMassCorrectionTags.Add(closestMassCorrectionTag, dblModificationMass);
                        }
                    }

                    return strClosestMassCorrectionTag;
                }
            }

            return string.Empty;
        }

        /// <summary>
        /// Looks for a modification of type .DynamicMod or type .UnknownType in mModifications having .ModificationSymbol = chModificationSymbol and chTargetResidue in .TargetResidues
        /// </summary>
        /// <param name="chModificationSymbol"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnExistingModFound"></param>
        /// <returns>Modification details</returns>
        /// <remarks>If chModificationSymbol does not match any of the mods, a modification with a mass of 0 is returned</remarks>
        public clsModificationDefinition LookupDynamicModificationDefinitionByTargetInfo(
            char chModificationSymbol,
            char chTargetResidue,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            out bool blnExistingModFound)
        {

            blnExistingModFound = false;
            if (chTargetResidue != default(char) || eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
            {
                // The residue was provided and/or the residue is located at a peptide or protein terminus
                // First compare against modifications with 1 or more residues in .TargetResidues
                for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
                {
                    if ((mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType) && mModifications[intIndex].TargetResidues.Length > 0)
                    {
                        if (mModifications[intIndex].ModificationSymbol == chModificationSymbol)
                        {
                            // Matching modification symbol found
                            // Now see if .TargetResidues contains chTargetResidue
                            if (mModifications[intIndex].TargetResiduesContain(chTargetResidue))
                            {
                                blnExistingModFound = true;
                            }

                            if (!blnExistingModFound && eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
                            {
                                switch (eResidueTerminusState)
                                {
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus:
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus:
                                        if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS))
                                        {
                                            blnExistingModFound = true;
                                        }
                                        break;
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus:
                                        if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                        {
                                            blnExistingModFound = true;
                                        }
                                        break;
                                }

                                if (!blnExistingModFound)
                                {
                                    switch (eResidueTerminusState)
                                    {
                                        case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus:
                                        case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus:
                                            if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS))
                                            {
                                                blnExistingModFound = true;
                                            }
                                            break;
                                        case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus:
                                            if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                            {
                                                blnExistingModFound = true;
                                            }
                                            break;
                                    }
                                }

                                if (!blnExistingModFound && (eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus || eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus))
                                {
                                    // Protein N-Terminus residue could also match a Peptide N-terminal mod; check for this
                                    if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        blnExistingModFound = true;
                                    }
                                }

                                if (!blnExistingModFound && (eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus || eResidueTerminusState == clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus))
                                {
                                    // Protein C-Terminus residue could also match a Peptide C-terminal mod; check for this
                                    if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        blnExistingModFound = true;
                                    }
                                }
                            }

                            if (blnExistingModFound)
                            {
                                // Match found
                                return mModifications[intIndex];
                            }
                        }
                    }
                }
            }

            // No match was found
            // First compare against modifications, only considering those with empty .TargetResidues
            // If still no match, then we'll try again but ignore .TargetResidues
            var blnConsiderTargetResidues = true;

            while (true)
            {
                for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
                {
                    if (mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType)
                    {
                        if (blnConsiderTargetResidues && string.IsNullOrWhiteSpace(mModifications[intIndex].TargetResidues) || !blnConsiderTargetResidues)
                        {
                            if (mModifications[intIndex].ModificationSymbol == chModificationSymbol)
                            {
                                // Matching mass found
                                blnExistingModFound = true;
                                return mModifications[intIndex];
                            }
                        }
                    }
                }

                if (blnConsiderTargetResidues)
                {
                    // No match; try again, but ignore .TargetResidues
                    blnConsiderTargetResidues = false;
                }
                else
                {
                    break;
                }
            }

            // Still no match; return a default modification with a mass of 0
            var objModificationDefinition = new clsModificationDefinition(chModificationSymbol, 0) {
                MassCorrectionTag = LookupMassCorrectionTagByMass(0)
            };

            return objModificationDefinition;
        }

        /// <summary>
        /// Looks for an existing modification with the given modification mass and target residues
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="chTargetResidue">If defined, then returns the first modification with the given mass and containing the residue in .TargetResidues; if no match, then looks for the first modification with the given mass and no defined .TargetResidues</param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnExistingModFound"></param>
        /// <param name="addToModificationListIfUnknown"></param>
        /// <param name="massDigitsOfPrecision"></param>
        /// <returns>The best matched modification; if no match is found, then returns a newly created modification definition, adding it to mModifications if addToModificationListIfUnknown is True</returns>
        /// <remarks>If chTargetResidue is nothing, then follows similar matching logic, but skips defined modifications with defined .TargetResidues</remarks>
        public clsModificationDefinition LookupModificationDefinitionByMass(double dblModificationMass,
            char chTargetResidue,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            out bool blnExistingModFound,
            bool addToModificationListIfUnknown,
            byte massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION)
        {
            clsModificationDefinition objModificationDefinition;

            blnExistingModFound = false;
            if (chTargetResidue != default(char) || eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
            {
                // The residue was provided and/or the residue is located at a peptide or protein terminus
                // First compare against modifications with 1 or more residues in .TargetResidues
                for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
                {
                    if (mModifications[intIndex].ModificationType != clsModificationDefinition.eModificationTypeConstants.DynamicMod &&
                        mModifications[intIndex].ModificationType != clsModificationDefinition.eModificationTypeConstants.StaticMod &&
                        mModifications[intIndex].ModificationType != clsModificationDefinition.eModificationTypeConstants.UnknownType ||
                        mModifications[intIndex].TargetResidues.Length <= 0)
                    {
                        continue;
                    }

                    if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) > float.Epsilon)
                    {
                        continue;
                    }

                    // Matching mass found
                    // Now see if .TargetResidues contains chTargetResidue
                    if (mModifications[intIndex].TargetResiduesContain(chTargetResidue))
                    {
                        blnExistingModFound = true;
                    }

                    if (!blnExistingModFound && eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
                    {
                        switch (eResidueTerminusState)
                        {
                            case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus:
                            case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus:
                            case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus:
                                if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                {
                                    blnExistingModFound = true;
                                }
                                break;
                            case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus:
                            case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus:
                                if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                {
                                    blnExistingModFound = true;
                                }
                                break;
                        }
                    }

                    if (blnExistingModFound)
                    {
                        // Match found
                        return mModifications[intIndex];
                    }
                }
            }

            // No match was found
            // Compare against modifications with empty .TargetResidues
            for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
            {
                if ((mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType) && string.IsNullOrWhiteSpace(mModifications[intIndex].TargetResidues))
                {
                    if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                    {
                        // Matching mass found
                        return mModifications[intIndex];
                    }
                }
            }

            // Still no match; look for the modification mass and residue in mStandardRefinementModifications
            // Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing
            if (chTargetResidue != default(char))
            {
                for (var intIndex = 0; intIndex <= mStandardRefinementModifications.Count - 1; intIndex++)
                {
                    if (Math.Abs(Math.Round(Math.Abs(mStandardRefinementModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                    {
                        // Matching mass found
                        // Now see if .TargetResidues contains chTargetResidue
                        if (mStandardRefinementModifications[intIndex].TargetResiduesContain(chTargetResidue))
                        {
                            blnExistingModFound = true;

                            objModificationDefinition = mStandardRefinementModifications[intIndex];
                            objModificationDefinition.ModificationSymbol = clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL;

                            if (addToModificationListIfUnknown && mDefaultModificationSymbols.Count > 0)
                            {
                                // Append objModificationDefinition to mModifications()
                                var intNewModIndex = AddModification(objModificationDefinition, true);
                                if (intNewModIndex >= 0)
                                {
                                    return mModifications[intNewModIndex];
                                }

                                return objModificationDefinition;
                            }

                            // Either addToModificationListIfUnknown = False or no more default modification symbols
                            // Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
                            return objModificationDefinition;
                        }
                    }
                }
            }

            // Still no match
            // Compare against dynamic and unknown-type modifications, but ignore .TargetResidues
            for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
            {
                if (mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType)
                {
                    if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                    {
                        // Matching mass found
                        // Assure that the target residues contain chTargetResidue
                        if (chTargetResidue != default(char) && !mModifications[intIndex].TargetResiduesContain(chTargetResidue))
                        {
                            mModifications[intIndex].TargetResidues += chTargetResidue;
                        }

                        return mModifications[intIndex];
                    }
                }
            }

            // Still no match; define a new custom modification
            const clsModificationDefinition.eModificationTypeConstants eModType = clsModificationDefinition.eModificationTypeConstants.DynamicMod;
            const char chModSymbol = clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL;
            const bool blnUseNextAvailableModificationSymbol = true;

            objModificationDefinition = AddUnknownModification(dblModificationMass, eModType, chTargetResidue, eResidueTerminusState, addToModificationListIfUnknown, blnUseNextAvailableModificationSymbol, chModSymbol, massDigitsOfPrecision);

            return objModificationDefinition;
        }

        /// <summary>
        /// Looks for an existing modification with the given modification mass, modification type, and target residues
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="eModType"></param>
        /// <param name="chTargetResidue">If defined, then returns the first modification with the given mass and containing the residue in .TargetResidues; if no match, then looks for the first modification with the given mass and no defined .TargetResidues</param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnExistingModFound"></param>
        /// <param name="addToModificationListIfUnknown"></param>
        /// <param name="massDigitsOfPrecision"></param>
        /// <returns>The best matched modification; if no match is found, then returns a newly created modification definition, adding it to mModifications if addToModificationListIfUnknown = True</returns>
        /// <remarks>If chTargetResidue is nothing, then follows similar matching logic, but skips defined modifications with defined .TargetResidues</remarks>
        public clsModificationDefinition LookupModificationDefinitionByMassAndModType(
            double dblModificationMass,
            clsModificationDefinition.eModificationTypeConstants eModType,
            char chTargetResidue,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            out bool blnExistingModFound,
            bool addToModificationListIfUnknown,
            byte massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION)
        {
            // If chTargetResidue is defined, then returns the first modification with the given mass and containing the residue in .TargetResidues
            //  If no match is found, then looks for the first modification with the given mass and no defined .TargetResidues
            //  If no match is found, then looks for the first dynamic modification with the given mass, regardless of .TargetResidues
            //  If no match is found, then returns a newly created modification definition, adding it to mModifications if addToModificationListIfUnknown = True
            // If chTargetResidue is nothing, then follows similar logic, but skips defined modifications with defined .TargetResidues

            clsModificationDefinition objModificationDefinition;

            char chModSymbol;
            bool blnUseNextAvailableModificationSymbol;

            switch (eModType)
            {
                case clsModificationDefinition.eModificationTypeConstants.StaticMod:
                case clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod:
                case clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod:
                    chModSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                    blnUseNextAvailableModificationSymbol = false;
                    break;
                default:
                    chModSymbol = clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL;
                    blnUseNextAvailableModificationSymbol = true;
                    break;
            }

            blnExistingModFound = false;
            if (chTargetResidue != default(char) || eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
            {
                // The residue was provided and/or the residue is located at a peptide or protein terminus
                // First compare against modifications with 1 or more residues in .TargetResidues
                for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
                {
                    if (mModifications[intIndex].ModificationType == eModType && (mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType) && mModifications[intIndex].TargetResidues.Length > 0)
                    {
                        if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                        {
                            // Matching mass found
                            // Now see if .TargetResidues contains chTargetResidue
                            if (chTargetResidue != default(char) && mModifications[intIndex].TargetResiduesContain(chTargetResidue))
                            {
                                blnExistingModFound = true;
                            }

                            if (!blnExistingModFound && eResidueTerminusState != clsAminoAcidModInfo.eResidueTerminusStateConstants.None)
                            {
                                switch (eResidueTerminusState)
                                {
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus:
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus:
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus:
                                        if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                        {
                                            blnExistingModFound = true;
                                        }
                                        break;
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus:
                                    case clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus:
                                        if (mModifications[intIndex].TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                        {
                                            blnExistingModFound = true;
                                        }
                                        break;
                                }
                            }

                            if (blnExistingModFound)
                            {
                                // Match found
                                return mModifications[intIndex];
                            }
                        }
                    }
                }
            }

            // No match was found
            // Compare against modifications with empty .TargetResidues
            for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
            {
                if (mModifications[intIndex].ModificationType == eModType && (mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod || mModifications[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType) && string.IsNullOrWhiteSpace(mModifications[intIndex].TargetResidues))
                {
                    if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                    {
                        // Matching mass found
                        return mModifications[intIndex];
                    }
                }
            }

            // Still no match; look for the modification mass and residue in mStandardRefinementModifications
            // Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing or chTargetResidue = '<' or chTargetResidue = '>'
            if (chTargetResidue != default(char))
            {
                for (var intIndex = 0; intIndex <= mStandardRefinementModifications.Count - 1; intIndex++)
                {
                    if (Math.Abs(Math.Round(Math.Abs(mStandardRefinementModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                    {
                        // Matching mass found
                        // Now see if .TargetResidues contains chTargetResidue
                        if (mStandardRefinementModifications[intIndex].TargetResiduesContain(chTargetResidue))
                        {
                            blnExistingModFound = true;

                            objModificationDefinition = mStandardRefinementModifications[intIndex];
                            objModificationDefinition.ModificationSymbol = chModSymbol;
                            objModificationDefinition.ModificationType = eModType;

                            if (addToModificationListIfUnknown && mDefaultModificationSymbols.Count > 0)
                            {
                                // Append objModificationDefinition to mModifications()
                                var intNewModIndex = AddModification(objModificationDefinition, true);
                                if (intNewModIndex >= 0)
                                {
                                    return mModifications[intNewModIndex];
                                }

                                return objModificationDefinition;
                            }

                            // Either addToModificationListIfUnknown = False or no more default modification symbols
                            // Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
                            return objModificationDefinition;
                        }
                    }
                }
            }

            // No match was found
            // Compare against modifications of the same type, but ignore .TargetResidues
            for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
            {
                if (mModifications[intIndex].ModificationType == eModType)
                {
                    if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) < float.Epsilon)
                    {
                        // Matching mass found
                        // Assure that the target residues contain chTargetResidue
                        if (chTargetResidue != default(char) && !mModifications[intIndex].TargetResiduesContain(chTargetResidue))
                        {
                            mModifications[intIndex].TargetResidues += chTargetResidue;
                        }

                        return mModifications[intIndex];
                    }
                }
            }

            // Still no match; define a new custom modification
            objModificationDefinition = AddUnknownModification(dblModificationMass, eModType, chTargetResidue, eResidueTerminusState, addToModificationListIfUnknown, blnUseNextAvailableModificationSymbol, chModSymbol, massDigitsOfPrecision);

            return objModificationDefinition;
        }

        private void InitializeLocalVariables()
        {
            mErrorMessage = string.Empty;
            SetDefaultMassCorrectionTags();

            // Note that this sub will call UpdateDefaultModificationSymbols()
            ClearModifications();

            UpdateStandardRefinementModifications();

            UpdateIntegerBasedModificationMap();
        }

        /// <summary>
        /// Load the mass correction tags file
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="blnFileNotFound"></param>
        /// <returns></returns>
        public bool ReadMassCorrectionTagsFile(string strFilePath, ref bool blnFileNotFound)
        {
            bool blnSuccess;

            try
            {
                // Open the mass correction tags file
                // It should have 2 columns, separated by tabs
                // Column 1 is the mass correction tag name
                // Column 2 is the monoisotopic mass for the mass correction (positive or negative number)

                if (string.IsNullOrWhiteSpace(strFilePath))
                {
                    SetDefaultMassCorrectionTags();
                    blnSuccess = true;
                }
                else if (!File.Exists(strFilePath))
                {
                    mErrorMessage = "Mass CorrectionTags File Not Found: " + strFilePath;
                    SetDefaultMassCorrectionTags();
                    blnFileNotFound = true;
                    blnSuccess = false;
                }
                else
                {
                    using (var srMassCorrectionTagsFile = new StreamReader(strFilePath))
                    {
                        mMassCorrectionTags.Clear();

                        while (!srMassCorrectionTagsFile.EndOfStream)
                        {
                            var strLineIn = srMassCorrectionTagsFile.ReadLine();
                            if (string.IsNullOrEmpty(strLineIn))
                                continue;

                            var strSplitLine = strLineIn.Split('\t');

                            if (strSplitLine.Length >= 2)
                            {
                                // See if the first column contains 1 or more characters and if the second column contains a number
                                // Note that StoreMassCorrectionTag() will trim spaces from the end of the mass correction tag names
                                if (strSplitLine[0].Trim().Length >= 1 && clsPHRPParser.IsNumber(strSplitLine[1]))
                                {
                                    StoreMassCorrectionTag(strSplitLine[0], double.Parse(strSplitLine[1]));
                                }
                            }
                        }
                    }

                    blnSuccess = true;

                    if (mMassCorrectionTags.Count == 0)
                    {
                        SetDefaultMassCorrectionTags();
                    }
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading Mass Correction Tags file (" + strFilePath + "): " + ex.Message;
                blnSuccess = false;
            }

            return blnSuccess;
        }

        /// <summary>
        /// Read a modification definitions file (_ModDefs.txt)
        /// </summary>
        /// <param name="strFilePath"></param>
        /// <param name="blnFileNotFound"></param>
        /// <returns></returns>
        public bool ReadModificationDefinitionsFile(string strFilePath, ref bool blnFileNotFound)
        {

            bool blnSuccess;

            try
            {
                // Open the modification file
                // It should have 2 or more columns, separated by tabs
                // Column 1 is the modification symbol
                // Column 2 is the modification mass
                // Column 3, which is optional, is the residues and/or terminii that can be modified; if omitted, then the modification can apply to any residues or terminii
                //   For column 3, use 1 letter amino acid abbreviations; the residues can be a continous string, or can be separated by commas and/or spaces
                //   For column 3, use the *_SYMBOL_DMS constants for the terminii (< and > for the peptide terminii; [ and ] for the protein terminii)
                // Column 4, which is optional, specifies the type of modification: D, S, T, I, or P (corresponding to clsModificationDefinition.eModificationTypeConstants)
                // Column 5, which is optional, specifies the mass correction tag associated with the given modification

                if (string.IsNullOrWhiteSpace(strFilePath))
                {
                    ClearModifications();
                    blnSuccess = true;
                }
                else if (!File.Exists(strFilePath))
                {
                    mErrorMessage = "Modification Definition File Not Found: " + strFilePath;
                    ClearModifications();
                    blnFileNotFound = true;
                    blnSuccess = false;
                }
                else
                {
                    using (var srModificationFile = new StreamReader(strFilePath))
                    {
                        ClearModifications();

                        while (!srModificationFile.EndOfStream)
                        {
                            var strLineIn = srModificationFile.ReadLine();
                            if (string.IsNullOrEmpty(strLineIn))
                                continue;

                            var strSplitLine = strLineIn.Split('\t');

                            if (strSplitLine.Length < 2)
                                continue;

                            // See if the first column contains a single character and if the second column contains a number

                            if (strSplitLine[0].Trim().Length != 1 || !clsPHRPParser.IsNumber(strSplitLine[1]))
                                continue;

                            var objModificationDefinition = new clsModificationDefinition(strSplitLine[0].Trim()[0], double.Parse(strSplitLine[1]));

                            if (strSplitLine.Length >= 3)
                            {
                                // Parse the target residues list
                                var strResidues = strSplitLine[2].Trim().ToUpper();

                                var strResiduesClean = string.Empty;
                                foreach (var chChar in strResidues)
                                {
                                    if (char.IsUpper(chChar))
                                    {
                                        strResiduesClean += chChar;
                                    }
                                    else if (chChar == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS ||
                                             chChar == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS ||
                                             chChar == clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS ||
                                             chChar == clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)
                                    {
                                        strResiduesClean += chChar;
                                    }
                                }

                                if (strResiduesClean.Length > 0)
                                {
                                    objModificationDefinition.TargetResidues = string.Copy(strResiduesClean);
                                }

                                if (strSplitLine.Length >= 4)
                                {
                                    // Store the modification type
                                    if (strSplitLine[3].Trim().Length == 1)
                                    {
                                        objModificationDefinition.ModificationType = clsModificationDefinition.ModificationSymbolToModificationType(strSplitLine[3].ToUpper().Trim()[0]);
                                    }

                                    // If the .ModificationType is unknown, then change it to Dynamic
                                    if (objModificationDefinition.ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType)
                                    {
                                        objModificationDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod;
                                    }

                                    if (strSplitLine.Length >= 5)
                                    {
                                        objModificationDefinition.MassCorrectionTag = strSplitLine[4].Trim();

                                        if (strSplitLine.Length >= 6)
                                        {
                                            strSplitLine[5] = strSplitLine[5].Trim();
                                            if (strSplitLine[5].Length > 0)
                                            {
                                                objModificationDefinition.AffectedAtom = strSplitLine[5][0];
                                            }
                                            else
                                            {
                                                objModificationDefinition.AffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
                                            }
                                        }
                                    }
                                }
                            }

                            // Check whether the modification type is Static and the .TargetResidues are one of: <>[]
                            // If so, update the modification type as needed
                            if (objModificationDefinition.TargetResidues != null && objModificationDefinition.TargetResidues.Trim().Length == 1 && objModificationDefinition.ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                            {
                                if (objModificationDefinition.TargetResidues[0] == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS |
                                    objModificationDefinition.TargetResidues[0] == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                                {
                                    objModificationDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod;
                                }
                                else if (objModificationDefinition.TargetResidues[0] == clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS |
                                         objModificationDefinition.TargetResidues[0] == clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)
                                {
                                    objModificationDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod;
                                }
                            }

                            // Validate some of the settings if the modification type is IsotopicMod or TerminalPeptideStaticMod or ProteinTerminusStaticMod
                            var blnValidMod = true;
                            switch (objModificationDefinition.ModificationType)
                            {
                                case clsModificationDefinition.eModificationTypeConstants.IsotopicMod:
                                    objModificationDefinition.ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                                    if (objModificationDefinition.AffectedAtom == clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL)
                                    {
                                        blnValidMod = false;
                                    }
                                    break;
                                case clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod:
                                    objModificationDefinition.ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                                    if (objModificationDefinition.TargetResidues != clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() &&
                                        objModificationDefinition.TargetResidues != clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                                    {
                                        blnValidMod = false;
                                    }
                                    break;
                                case clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod:
                                    objModificationDefinition.ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                                    if (objModificationDefinition.TargetResidues != clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString() &&
                                        objModificationDefinition.TargetResidues != clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString())
                                    {
                                        blnValidMod = false;
                                    }
                                    break;
                                case clsModificationDefinition.eModificationTypeConstants.UnknownType:
                                    objModificationDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod;
                                    break;
                            }

                            if (objModificationDefinition.MassCorrectionTag == clsModificationDefinition.INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME)
                            {
                                // Try to determine the mass correction name
                                objModificationDefinition.MassCorrectionTag = LookupMassCorrectionTagByMass(objModificationDefinition.ModificationMass);
                            }

                            if (blnValidMod)
                            {
                                AddModification(objModificationDefinition, false);
                            }
                        }
                    }

                    // Note that this sub will call UpdateDefaultModificationSymbols()
                    ValidateModificationsVsDefaultModificationSymbols();
                    blnSuccess = true;
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading Modification Definition file (" + strFilePath + "): " + ex.Message;
                blnSuccess = false;
            }

            return blnSuccess;
        }

        /// <summary>
        /// Reset mod occurrence count stats
        /// </summary>
        public void ResetOccurrenceCountStats()
        {
            for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
            {
                mModifications[intIndex].OccurrenceCount = 0;
            }
        }

        /// <summary>
        /// Define the default mass correction tags
        /// </summary>
        public void SetDefaultMassCorrectionTags()
        {
            try
            {
                mMassCorrectionTags.Clear();

                // Note: Function StoreMassCorrectionTag will remove spaces
                // from the beginning or end of the mass correction tag names
                StoreMassCorrectionTag("4xDeut  ", 4.025107);
                StoreMassCorrectionTag("6C134N15", 10.008269);
                StoreMassCorrectionTag("6xC13N15", 7.017164);
                StoreMassCorrectionTag("AcetAmid", 41.02655);
                StoreMassCorrectionTag("Acetyl  ", 42.010567);
                StoreMassCorrectionTag("Acrylmid", 71.037117);
                StoreMassCorrectionTag("ADPRibos", 541.061096);
                StoreMassCorrectionTag("AlkSulf ", -25.0316);
                StoreMassCorrectionTag("Aminaton", 15.010899);
                StoreMassCorrectionTag("AmOxButa", -2.01565);
                StoreMassCorrectionTag("Bromo   ", 77.910507);
                StoreMassCorrectionTag("BS3Olnk ", 156.078644);
                StoreMassCorrectionTag("C13DtFrm", 36.07567);
                StoreMassCorrectionTag("Carbamyl", 43.005814);
                StoreMassCorrectionTag("Cyano   ", 24.995249);
                StoreMassCorrectionTag("Cys-Dha ", -33.98772);
                StoreMassCorrectionTag("Cystnyl ", 119.004097);
                StoreMassCorrectionTag("Deamide ", 0.984016);
                StoreMassCorrectionTag("DeutForm", 32.056407);
                StoreMassCorrectionTag("DeutMeth", 17.034479);
                StoreMassCorrectionTag("Dimethyl", 28.0313);
                StoreMassCorrectionTag("DTBP_Alk", 144.03573);
                StoreMassCorrectionTag("Formyl  ", 27.994915);
                StoreMassCorrectionTag("GalNAFuc", 648.2603);
                StoreMassCorrectionTag("GalNAMan", 664.2551);
                StoreMassCorrectionTag("Gluthone", 305.068146);
                StoreMassCorrectionTag("Guanid  ", 42.021797);
                StoreMassCorrectionTag("Heme_615", 615.169458);
                StoreMassCorrectionTag("Hexosam ", 203.079376);
                StoreMassCorrectionTag("Hexose  ", 162.052826);
                StoreMassCorrectionTag("ICAT_D0 ", 442.225006);
                StoreMassCorrectionTag("ICAT_D8 ", 450.275208);
                StoreMassCorrectionTag("IodoAcet", 57.021465);
                StoreMassCorrectionTag("IodoAcid", 58.005478);
                StoreMassCorrectionTag("Iso_N15 ", 0.997035);
                StoreMassCorrectionTag("itrac   ", 144.102066);
                StoreMassCorrectionTag("iTRAQ8  ", 304.205353);
                StoreMassCorrectionTag("LeuToMet", 17.956421);
                StoreMassCorrectionTag("Lipid2  ", 576.51178);
                StoreMassCorrectionTag("Mercury ", 199.9549);
                StoreMassCorrectionTag("Met_O18 ", 16.028204);
                StoreMassCorrectionTag("Methyl  ", 14.01565);
                StoreMassCorrectionTag("Methylmn", 13.031634);
                StoreMassCorrectionTag("MinusH2O", -18.010565);
                StoreMassCorrectionTag("NEM     ", 125.047676);
                StoreMassCorrectionTag("NH3_Loss", -17.026548);
                StoreMassCorrectionTag("NHS_SS  ", 87.998283);
                StoreMassCorrectionTag("NO2_Addn", 44.985077);
                StoreMassCorrectionTag("None    ", 0);
                StoreMassCorrectionTag("OMinus2H", 13.979265);
                StoreMassCorrectionTag("One_C12 ", 12);
                StoreMassCorrectionTag("One_O18 ", 2.004246);
                StoreMassCorrectionTag("OxoAla  ", -17.992805);
                StoreMassCorrectionTag("palmtlic", 236.21402);
                StoreMassCorrectionTag("PCGalNAz", 502.202332);
                StoreMassCorrectionTag("PEO     ", 414.193695);
                StoreMassCorrectionTag("PhosAden", 329.052521);
                StoreMassCorrectionTag("Phosph  ", 79.966331);
                StoreMassCorrectionTag("PhosUrid", 306.025299);
                StoreMassCorrectionTag("Plus1Oxy", 15.994915);
                StoreMassCorrectionTag("Plus2Oxy", 31.989828);
                StoreMassCorrectionTag("Plus3Oxy", 47.984745);
                StoreMassCorrectionTag("Propnyl ", 56.026215);
                StoreMassCorrectionTag("Pyro-cmC", 39.994915);
                StoreMassCorrectionTag("SATA_Alk", 131.0041);
                StoreMassCorrectionTag("SATA_Lgt", 115.9932);
                StoreMassCorrectionTag("Sucinate", 116.010956);
                StoreMassCorrectionTag("SulfoNHS", 226.077591);
                StoreMassCorrectionTag("Sumoylat", 484.228149);
                StoreMassCorrectionTag("TMT0Tag ", 224.152481);
                StoreMassCorrectionTag("TMT6Tag ", 229.162933);
                StoreMassCorrectionTag("TriMeth ", 42.046951);
                StoreMassCorrectionTag("Two_O18 ", 4.008491);
                StoreMassCorrectionTag("Ubiq_02 ", 114.042931);
                StoreMassCorrectionTag("Ubiq_L  ", 100.016045);
                StoreMassCorrectionTag("ValToMet", 31.972071);
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        private void StoreMassCorrectionTag(string strTagName, double dblMass)
        {
            try
            {
                mMassCorrectionTags.Add(strTagName.Trim(), dblMass);
            }
            catch (Exception)
            {
                // If a duplicate is tag is entered into the mMassCorrectionTags hashtable, an error will occur; we'll ignore the error
                // Ignore errors here
            }
        }

        private void UpdateDefaultModificationSymbols(string strModificationChars)
        {
            try
            {
                if (string.IsNullOrEmpty(strModificationChars))
                    return;

                if (mDefaultModificationSymbols == null)
                {
                    mDefaultModificationSymbols = new Queue();
                }
                else
                {
                    mDefaultModificationSymbols.Clear();
                }

                // Populate mDefaultModificationSymbols, making sure no characters are duplicated
                // In addition, do not allow LAST_RESORT_MODIFICATION_SYMBOL or NO_SYMBOL_MODIFICATION_SYMBOL to be used
                foreach (var chChar in strModificationChars)
                {
                    if (chChar != clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL &&
                        chChar != clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL)
                    {
                        if (!mDefaultModificationSymbols.Contains(chChar))
                        {
                            mDefaultModificationSymbols.Enqueue(chChar);
                        }
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        private void UpdateIntegerBasedModificationMap()
        {
            mIntegerMassCorrectionTagLookup = new Dictionary<int, string>
            {
                {-18, "MinusH2O"},
                {-17, "NH3_Loss"},
                {-11, "AsnToCys"},
                {-8, "HisToGlu"},
                {-7, "TyrToArg"},
                {-4, "ThrToPro"},
                {-3, "MetToLys"},
                {-1, "Dehydro"},
                {1, "Deamide"},
                {2, "GluToMet"},
                {4, "TrypOxy"},
                {5, "5C13"},
                {6, "6C13"},
                {10, "D10-Leu"},
                {13, "Methylmn"},
                {14, "Methyl"},
                {16, "Plus1Oxy"},
                {18, "LeuToMet"},
                {25, "Cyano"},
                {28, "Dimethyl"},
                {32, "Plus2Oxy"},
                {42, "Acetyl"},
                {43, "Carbamyl"},
                {45, "NO2_Addn"},
                {48, "Plus3Oxy"},
                {56, "Propnyl"},
                {58, "IodoAcid"},
                {80, "Phosph"},
                {89, "Biotinyl"},
                {96, "PhosphH"},
                {104, "Ubiq_H"},
                {116, "Sucinate"},
                {119, "Cystnyl"},
                {125, "NEM"},
                {144, "itrac"},
                {215, "MethylHg"},
                {236, "ICAT_C13"},
                {442, "ICAT_D0"}
            };

        }

        private void UpdateStandardRefinementModifications()
        {
            mStandardRefinementModifications = new List<clsModificationDefinition>();

            var dblModificationMass = -17.026549;
            mStandardRefinementModifications.Add(new clsModificationDefinition(
                clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL,
                dblModificationMass,
                "Q",
                clsModificationDefinition.eModificationTypeConstants.DynamicMod,
                LookupMassCorrectionTagByMass(dblModificationMass)));

            dblModificationMass = -18.0106;
            mStandardRefinementModifications.Add(new clsModificationDefinition(
                clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL,
                dblModificationMass,
                "E",
                clsModificationDefinition.eModificationTypeConstants.DynamicMod,
                LookupMassCorrectionTagByMass(dblModificationMass)));
        }

        private void ValidateModificationsVsDefaultModificationSymbols()
        {
            try
            {
                // Reset the default modification symbols list
                UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS);

                var chDefaultModificationSymbols = new char[mDefaultModificationSymbols.Count];
                mDefaultModificationSymbols.ToArray().CopyTo(chDefaultModificationSymbols, 0);
                var intDefaultModificationSymbolCount = chDefaultModificationSymbols.Length;

                // Step through mModifications and make sure each of the modification symbols is not present in mDefaultModificationChars
                for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
                {
                    var intIndexCompare = 0;
                    while (intIndexCompare < intDefaultModificationSymbolCount)
                    {
                        if (mModifications[intIndex].ModificationSymbol == chDefaultModificationSymbols[intIndexCompare])
                        {
                            // Remove this symbol from chDefaultModificationSymbols
                            for (var intIndexCopy = intIndexCompare; intIndexCopy <= intDefaultModificationSymbolCount - 2; intIndexCopy++)
                            {
                                chDefaultModificationSymbols[intIndexCopy] = chDefaultModificationSymbols[intIndexCopy + 1];
                            }
                            intDefaultModificationSymbolCount -= 1;
                        }
                        else
                        {
                            intIndexCompare += 1;
                        }
                    }
                }

                if (intDefaultModificationSymbolCount < mDefaultModificationSymbols.Count)
                {
                    mDefaultModificationSymbols.Clear();
                    for (var intIndex = 0; intIndex <= intDefaultModificationSymbolCount - 1; intIndex++)
                    {
                        mDefaultModificationSymbols.Enqueue(chDefaultModificationSymbols[intIndex]);
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        /// <summary>
        /// Verify that a modification is present, adding it if missing
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="strTargetResidues"></param>
        /// <param name="eModificationType"></param>
        /// <param name="massDigitsOfPrecision"></param>
        /// <returns>True if the modification was matched or was added; false if an error</returns>
        public bool VerifyModificationPresent(
            double dblModificationMass,
            string strTargetResidues,
            clsModificationDefinition.eModificationTypeConstants eModificationType,
            int massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION)
        {
            // Returns True if the modification was matched or was added
            // Returns False if an error

            // Look for mods in mModifications with matching .ModificationType, .ModificationMass (within tolerance),
            //   and .TargetResidues vs. udtModDefintion
            // If not found, add a new entry to mModifications

            var blnMatchFound = false;

            if (massDigitsOfPrecision < 0)
                massDigitsOfPrecision = 0;

            try
            {
                for (var intIndex = 0; intIndex <= mModifications.Count - 1; intIndex++)
                {
                    if (mModifications[intIndex].ModificationType != eModificationType)
                        continue;

                    // Matching modification type
                    if (Math.Abs(Math.Round(Math.Abs(mModifications[intIndex].ModificationMass - dblModificationMass), massDigitsOfPrecision)) > float.Epsilon)
                        continue;

                    // Matching mass
                    // Compare .TargetResidues
                    blnMatchFound = clsModificationDefinition.EquivalentTargetResidues(mModifications[intIndex].TargetResidues, strTargetResidues, true);
                    if (blnMatchFound)
                    {
                        break;
                    }
                }

                if (!blnMatchFound)
                {
                    var objModificationDefinition = new clsModificationDefinition(dblModificationMass, strTargetResidues, eModificationType)
                    {
                        MassCorrectionTag = LookupMassCorrectionTagByMass(dblModificationMass)
                    };

                    // Append objModificationDefinition to mModifications()
                    AddModification(objModificationDefinition, true);

                    blnMatchFound = true;
                }
            }
            catch (Exception)
            {
                blnMatchFound = false;
            }

            return blnMatchFound;
        }
    }
}
