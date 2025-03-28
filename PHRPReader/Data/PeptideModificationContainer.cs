﻿// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 5, 2006
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
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PHRPReader.Reader;

namespace PHRPReader.Data
{
    /// <summary>
    /// This class is used to track modifications that can be applied to peptides
    /// It handles both residue level modifications and static, peptide-wide modifications
    /// </summary>
    /// <remarks>
    /// Use ReadMassCorrectionTagsFile() and ReadModificationDefinitionsFile() to customize
    /// the default mass correction tag and modification definition lists
    /// </remarks>
    public class PeptideModificationContainer
    {
        // Ignore Spelling: acetyl, carbamyl, deamidated, didehydro, dimethylation, methylation
        // Ignore Spelling: phospho, phosphorylation, Sequest, UniMod, XTandem

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
        /// Symbol used by InSpecT for tracking a modification at the peptide N-terminus
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const char N_TERMINAL_PEPTIDE_MOD_SYMBOL_INSPECT = '[';

        /// <summary>
        /// Symbol used by InSpecT for tracking a modification at the peptide C-terminus
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public const char C_TERMINAL_PEPTIDE_MOD_SYMBOL_INSPECT = ']';

        /// <summary>
        /// List of available modification symbols
        /// </summary>
        private Queue mDefaultModificationSymbols;

        /// <summary>
        /// List of known mass correction tags
        /// </summary>
        /// <remarks>Keys are mod names, values are mod masses</remarks>
        private readonly Dictionary<string, double> mMassCorrectionTags;

        /// <summary>
        /// This array holds modifications that Sequest or XTandem will often use but for
        /// which the auto-addition method sometimes incorrectly notes
        /// </summary>
        private List<ModificationDefinition> mStandardRefinementModifications;

        private Dictionary<int, string> mIntegerMassCorrectionTagLookup;

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage { get; private set; }

        /// <summary>
        /// Modification count
        /// </summary>
        public int ModificationCount => Modifications.Count;

        /// <summary>
        /// Modification list
        /// </summary>
        public List<ModificationDefinition> Modifications { get; }

        /// <summary>
        /// When true, take the mod symbol into account when finding identical mods
        /// </summary>
        public bool ConsiderModSymbolWhenFindingIdenticalMods { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public PeptideModificationContainer()
        {
            mMassCorrectionTags = new Dictionary<string, double>();
            Modifications = new List<ModificationDefinition>();

            InitializeLocalVariables();
        }

        /// <summary>
        /// Add modificationDefinition to mModifications
        /// However, do not add if a duplicate modification
        /// Furthermore, if everything matches except for .TargetResidues, add the new target residues to the existing, matching mod
        /// </summary>
        /// <param name="modificationDefinition">Modification definition</param>
        /// <param name="useNextAvailableModificationSymbol">When true, use the next available modification symbol</param>
        /// <returns>The index of the newly added modification, or the index of the modification that modificationDefinition matches </returns>
        private int AddModification(ModificationDefinition modificationDefinition, bool useNextAvailableModificationSymbol)
        {
            int modificationIndex;

            var matchFound = false;

            // See if any of the existing modifications match modificationDefinition, ignoring .TargetResidues and possibly ignoring .ModificationSymbol
            for (modificationIndex = 0; modificationIndex <= Modifications.Count - 1; modificationIndex++)
            {
                if (Modifications[modificationIndex].EquivalentMassTypeTagAndAtom(modificationDefinition))
                {
                    matchFound = true;

                    if (ConsiderModSymbolWhenFindingIdenticalMods)
                    {
                        if (modificationDefinition.ModificationSymbol != Modifications[modificationIndex].ModificationSymbol)
                        {
                            // Symbols differ; add this as a new modification definition
                            matchFound = false;
                        }
                    }

                    if (matchFound)
                    {
                        var mod = Modifications[modificationIndex];

                        if (mod.ModificationType is ModificationDefinition.ResidueModificationType.DynamicMod or ModificationDefinition.ResidueModificationType.StaticMod)
                        {
                            // Matching dynamic or static modification definitions
                            // Merge the two modifications by making sure each of the residues in modificationDefinition.TargetResidues is present in .TargetResidues
                            foreach (var chChar in modificationDefinition.TargetResidues)
                            {
                                if (!mod.TargetResiduesContain(chChar))
                                {
                                    mod.TargetResidues += chChar;
                                }
                            }

                            if (!useNextAvailableModificationSymbol)
                            {
                                // See if the new modification symbol is different from the already-defined symbol

                                if (modificationDefinition.ModificationSymbol != mod.ModificationSymbol)
                                {
                                }
                            }
                        }
                    }
                }
                if (matchFound)
                    break;
            }

            if (matchFound)
                return modificationIndex;

            if (useNextAvailableModificationSymbol && mDefaultModificationSymbols.Count > 0)
            {
                // Add modificationDefinition to the list, using the next available default modification symbol
                modificationDefinition.ModificationSymbol = (char)mDefaultModificationSymbols.Dequeue();
            }

            Modifications.Add(modificationDefinition);

            return Modifications.Count - 1;
        }

        private ModificationDefinition AddUnknownModification(
            double modificationMass,
            ModificationDefinition.ResidueModificationType modType,
            char chTargetResidue,
            AminoAcidModInfo.ResidueTerminusState residueTerminusState,
            bool addToModificationListIfUnknown,
            bool useNextAvailableModificationSymbol,
            char modSymbol,
            byte massDigitsOfPrecision,
            byte massDigitsOfPrecisionLoose)
        {
            string targetResidues;

            if (chTargetResidue == '\0')
            {
                targetResidues = string.Empty;
            }
            else
            {
                targetResidues = chTargetResidue.ToString();
            }

            if (residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
            {
                // Assume this is a terminus mod
                // Override the target residues
                switch (residueTerminusState)
                {
                    case AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus:
                    case AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus:
                    case AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus:
                        targetResidues = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        break;
                    case AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus:
                    case AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus:
                        targetResidues = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                        break;
                    default:
                        throw new Exception("Unrecognized ResidueTerminusStateConstants enum: " + residueTerminusState);
                }
            }

            if (!useNextAvailableModificationSymbol)
            {
                modSymbol = ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
            }

            var massCorrectionTag = LookupMassCorrectionTagByMass(modificationMass, massDigitsOfPrecision, true, massDigitsOfPrecisionLoose);

            var modificationDefinition = new ModificationDefinition(
                modSymbol,
                modificationMass,
                targetResidues,
                modType,
                massCorrectionTag,
                PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, true);

            if (!addToModificationListIfUnknown)
                return modificationDefinition;

            // Append modificationDefinition to mModifications()
            int newModIndex;
            if (mDefaultModificationSymbols.Count > 0 && useNextAvailableModificationSymbol)
            {
                newModIndex = AddModification(modificationDefinition, useNextAvailableModificationSymbol: true);
            }
            else
            {
                newModIndex = AddModification(modificationDefinition, useNextAvailableModificationSymbol: false);
            }

            return newModIndex >= 0 ? Modifications[newModIndex] : modificationDefinition;

            // Either addToModificationListIfUnknown = False or no more default modification symbols
            // Return modificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
        }

        /// <summary>
        /// Append standard refinement modifications
        /// </summary>
        public void AppendStandardRefinementModifications()
        {
            for (var index = 0; index <= mStandardRefinementModifications.Count - 1; index++)
            {
                var mod = mStandardRefinementModifications[index];
                VerifyModificationPresent(mod.ModificationMass, mod.TargetResidues, mod.ModificationType);
            }
        }

        /// <summary>
        /// Clear modifications
        /// </summary>
        public void ClearModifications()
        {
            UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS);
            Modifications.Clear();
        }

        // ReSharper disable once ParameterTypeCanBeEnumerable.Local

        /// <summary>
        /// Find the modification in matchedMods with the smallest delta mass value
        /// </summary>
        /// <param name="matchedMods">List of modifications where keys are absolute value of delta mass and values are the modification definition</param>
        /// <param name="chTargetResidue">Target residue</param>
        /// <param name="addTargetResidue">If true, add the target residue to the modification definition's target residues list (if missing)</param>
        private ModificationDefinition FindClosestMatchedMod(
            IReadOnlyCollection<KeyValuePair<double, ModificationDefinition>> matchedMods,
            char chTargetResidue,
            bool addTargetResidue = false)
        {
            var closestMatch = matchedMods.First();

            foreach (var mod in matchedMods.Skip(1))
            {
                if (mod.Key < closestMatch.Key)
                {
                    closestMatch = mod;
                }
            }

            if (!addTargetResidue)
                return closestMatch.Value;

            // Assure that the target residues contain chTargetResidue
            if (chTargetResidue != '\0' && !closestMatch.Value.TargetResiduesContain(chTargetResidue))
            {
                closestMatch.Value.TargetResidues += chTargetResidue;
            }

            return closestMatch.Value;
        }

        /// <summary>
        /// Converts a modification mass to a generic 8 character name
        /// The name will always start with + or - then will have the modification mass, rounded as necessary to give an 8 character name
        /// </summary>
        /// <param name="modificationMass">Modification mass</param>
        private string GenerateGenericModMassName(double modificationMass)
        {
            int formatDigits;
            int maxLength;

            if (Math.Abs(modificationMass) < float.Epsilon)
            {
                return "+0.00000";
            }

            if (modificationMass < -9999999)
            {
                // Modification mass is too negative; always return -9999999
                return "-9999999";
            }

            if (modificationMass > 9999999)
            {
                // Modification mass is too positive; always return +9999999
                return "+9999999";
            }

            // Determine the number of digits that we will display to the left of the decimal point
            var formatDigitsDouble = Math.Log10(Math.Abs(modificationMass));

            if (Math.Abs(formatDigitsDouble - (int)formatDigitsDouble) < float.Epsilon)
            {
                // ModMass is a power of 10
                formatDigits = (int)formatDigitsDouble + 1;
            }
            else
            {
                formatDigits = (int)Math.Ceiling(formatDigitsDouble);
            }

            if (formatDigits < 1)
                formatDigits = 1;

            // Generate the format string
            // For example, a modification mass of 15.9994 will have formatString = "+00.0000"
            // Negative modification masses do not start with a minus sign in the format string since Visual Studio auto-adds it
            // Thus, a modification mass of -15.9994 will have formatString = "00.0000"
            var formatString = new string('0', formatDigits);

            if (modificationMass > 0)
            {
                formatString = "+" + formatString;
                maxLength = 8;
            }
            else
            {
                maxLength = 7;
            }

            if (formatString.Length < maxLength)
            {
                formatString += ".";
            }

            while (formatString.Length < maxLength)
            {
                formatString += "0";
            }

            var modMassName = modificationMass.ToString(formatString);

            if (modMassName.Length < 8 && modMassName.IndexOf('.') < 0)
            {
                modMassName += ".";
            }

            while (modMassName.Length < 8)
            {
                // Append extra zeros (this code will likely never be reached)
                modMassName += "0";
            }

            return modMassName.Length > 8
                ? throw new ArgumentOutOfRangeException("Generated Mod Name is longer than 8 characters: " + modMassName)
                : modMassName;
        }

        /// <summary>
        /// Looks for the best match in mIntegerMassCorrectionTagLookup for modificationMass (which should be close to an integer value)
        /// </summary>
        /// <param name="modificationMass">Modification mass</param>
        /// <returns>The mass correction tag name if a match, otherwise an empty string</returns>
        private string GetBestIntegerBasedMassCorrectionTag(double modificationMass)
        {
            var closestMassCorrectionTag = string.Empty;

            foreach (var tagOverride in mIntegerMassCorrectionTagLookup)
            {
                if (Math.Abs(modificationMass - tagOverride.Key) < 0.0001)
                {
                    closestMassCorrectionTag = tagOverride.Value;
                    break;
                }
            }

            return closestMassCorrectionTag;
        }

        /// <summary>
        /// Get a modification, by index
        /// </summary>
        /// <param name="index">Index</param>
        public ModificationDefinition GetModificationByIndex(int index)
        {
            if (index >= 0 && index < Modifications.Count)
            {
                return Modifications[index];
            }

            return new ModificationDefinition();
        }

        /// <summary>
        /// Get the modification type, by modification index
        /// </summary>
        /// <param name="index">Index</param>
        public ModificationDefinition.ResidueModificationType GetModificationTypeByIndex(int index)
        {
            if (index >= 0 && index < Modifications.Count)
            {
                return Modifications[index].ModificationType;
            }

            return ModificationDefinition.ResidueModificationType.UnknownType;
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not found
        /// </summary>
        /// <remarks>
        /// Searches known mods using 3 digits of precision, then 2 digits, then 1 digit
        /// If no match, adds as a new, unknown modification
        /// </remarks>
        /// <param name="modificationMass">Modification mass</param>
        /// <returns>Mod name, or empty string if no match</returns>
        public string LookupMassCorrectionTagByMass(double modificationMass)
        {
            const byte massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION;
            const bool addToModificationListIfUnknown = true;

            return LookupMassCorrectionTagByMass(modificationMass, massDigitsOfPrecision, addToModificationListIfUnknown);
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not  found
        /// </summary>
        /// <remarks>
        /// Searches known mods using massDigitsOfPrecision digits of precision, then massDigitsOfPrecision-1 digits, ... 1 digit
        /// If no match, adds as a new, unknown modification
        /// </remarks>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="massDigitsOfPrecision">Mass digits of precision</param>
        /// <returns>Mod name, or empty string if no match</returns>
        // ReSharper disable once UnusedMember.Global
        public string LookupMassCorrectionTagByMass(double modificationMass, byte massDigitsOfPrecision)
        {
            const bool addToModificationListIfUnknown = true;

            return LookupMassCorrectionTagByMass(modificationMass, massDigitsOfPrecision, addToModificationListIfUnknown);
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not found and addToModificationListIfUnknown is true
        /// </summary>
        /// <remarks>
        /// Searches known mods using massDigitsOfPrecision digits of precision, then massDigitsOfPrecision-1 digits, ... 1 digit
        /// If no match, adds as a new, unknown modification
        /// </remarks>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="massDigitsOfPrecision">Mass digits of precision</param>
        /// <param name="addToModificationListIfUnknown">When true, add to the modification list if unknown</param>
        /// <returns>Mod name, or empty string if no match</returns>
        public string LookupMassCorrectionTagByMass(double modificationMass, byte massDigitsOfPrecision, bool addToModificationListIfUnknown)
        {
            const byte massDigitsOfPrecisionLoose = 1;

            return LookupMassCorrectionTagByMass(modificationMass, massDigitsOfPrecision, addToModificationListIfUnknown, massDigitsOfPrecisionLoose);
        }

        /// <summary>
        /// Find the mass correction tag with the given mass, adding to the unknown modification list if not found and addToModificationListIfUnknown is true
        /// </summary>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="massDigitsOfPrecision">Number of digits after the decimal point to round to when comparing mod masses</param>
        /// <param name="addToModificationListIfUnknown">When true, add to the modification list if unknown</param>
        /// <param name="massDigitsOfPrecisionLoose">Number of digits after the decimal point to round to, for a more lenient match (if no match found using massDigitsOfPrecision)</param>
        /// <returns>Mod name, or empty string if no match</returns>
        public string LookupMassCorrectionTagByMass(
            double modificationMass,
            byte massDigitsOfPrecision,
            bool addToModificationListIfUnknown,
            byte massDigitsOfPrecisionLoose)
        {
            int massDigitsOfPrecisionStop;

            if (massDigitsOfPrecision < massDigitsOfPrecisionLoose)
                massDigitsOfPrecisionLoose = massDigitsOfPrecision;

            if (massDigitsOfPrecision >= massDigitsOfPrecisionLoose)
            {
                massDigitsOfPrecisionStop = massDigitsOfPrecisionLoose;
            }
            else
            {
                massDigitsOfPrecisionStop = massDigitsOfPrecision;
            }

            for (var currentPrecision = (int)massDigitsOfPrecision; currentPrecision >= massDigitsOfPrecisionStop; currentPrecision += -1)
            {
                var closestMassCorrectionTag = string.Empty;
                var closestMassCorrectionTagMassDiff = double.MaxValue;

                try
                {
                    if (massDigitsOfPrecisionStop == 0)
                    {
                        closestMassCorrectionTag = GetBestIntegerBasedMassCorrectionTag(modificationMass);
                        if (!string.IsNullOrEmpty(closestMassCorrectionTag))
                        {
                            return closestMassCorrectionTag;
                        }
                    }

                    // First look for an exact match in mMassCorrectionTags
                    foreach (var massCorrectionTag in mMassCorrectionTags)
                    {
                        var massDiff = Math.Abs(modificationMass - massCorrectionTag.Value);

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

                if (Math.Abs(Math.Round(closestMassCorrectionTagMassDiff, currentPrecision)) < float.Epsilon)
                {
                    // Match found
                    return closestMassCorrectionTag;
                }

                if (currentPrecision > massDigitsOfPrecisionStop)
                {
                    // Let the For loop go through another iteration to see if we find a match
                }
                else
                {
                    // Match not found
                    // Name the modification based on the mod mass
                    closestMassCorrectionTag = GenerateGenericModMassName(modificationMass);

                    if (addToModificationListIfUnknown)
                    {
                        if (mMassCorrectionTags.ContainsKey(closestMassCorrectionTag))
                        {
                            Console.WriteLine("Warning: Ignoring duplicate mass correction tag: {0}, mass {1:F3}",
                                              closestMassCorrectionTag, modificationMass);
                        }
                        else
                        {
                            mMassCorrectionTags.Add(closestMassCorrectionTag, modificationMass);
                        }
                    }

                    return closestMassCorrectionTag;
                }
            }

            return string.Empty;
        }

        /// <summary>
        /// Looks for a modification of type .DynamicMod or type .UnknownType in mModifications having .ModificationSymbol = modificationSymbol and chTargetResidue in .TargetResidues
        /// </summary>
        /// <remarks>If modificationSymbol does not match any of the mods, a modification with a mass of 0 is returned</remarks>
        /// <param name="modificationSymbol">Modification symbol</param>
        /// <param name="chTargetResidue">Target residue</param>
        /// <param name="residueTerminusState">Residue terminus state</param>
        /// <param name="existingModFound">Output: true if an existing modification was found</param>
        /// <returns>Modification details</returns>
        public ModificationDefinition LookupDynamicModificationDefinitionByTargetInfo(
            char modificationSymbol,
            char chTargetResidue,
            AminoAcidModInfo.ResidueTerminusState residueTerminusState,
            out bool existingModFound)
        {
            existingModFound = false;
            if (chTargetResidue != '\0' || residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
            {
                // The residue was provided and/or the residue is located at a peptide or protein terminus
                // First compare against modifications with 1 or more residues in .TargetResidues
                for (var index = 0; index <= Modifications.Count - 1; index++)
                {
                    if (Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.DynamicMod &&
                        Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.UnknownType ||
                        Modifications[index].TargetResidues.Length == 0)
                    {
                        continue;
                    }

                    if (Modifications[index].ModificationSymbol != modificationSymbol)
                    {
                        continue;
                    }

                    // Matching modification symbol found
                    // Now see if .TargetResidues contains chTargetResidue
                    if (Modifications[index].TargetResiduesContain(chTargetResidue))
                    {
                        existingModFound = true;
                    }

                    if (!existingModFound && residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
                    {
                        switch (residueTerminusState)
                        {
                            case AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus:
                            case AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus:
                                if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS))
                                {
                                    existingModFound = true;
                                }
                                break;
                            case AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus:
                                if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                {
                                    existingModFound = true;
                                }
                                break;
                        }

                        if (!existingModFound)
                        {
                            switch (residueTerminusState)
                            {
                                case AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus:
                                    if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS))
                                    {
                                        existingModFound = true;
                                    }
                                    break;
                                case AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus:
                                    if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        existingModFound = true;
                                    }
                                    break;
                            }
                        }

                        if (!existingModFound && residueTerminusState is AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus or AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus)
                        {
                            // Protein N-Terminus residue could also match a Peptide N-terminal mod; check for this
                            if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                            {
                                existingModFound = true;
                            }
                        }

                        if (!existingModFound && residueTerminusState is AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus or AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus)
                        {
                            // Protein C-Terminus residue could also match a Peptide C-terminal mod; check for this
                            if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                            {
                                existingModFound = true;
                            }
                        }
                    }

                    if (existingModFound)
                    {
                        // Match found
                        return Modifications[index];
                    }
                }
            }

            // No match was found
            // First compare against modifications, only considering those with empty .TargetResidues
            // If still no match, we'll try again but ignore .TargetResidues
            var considerTargetResidues = true;

            while (true)
            {
                for (var index = 0; index <= Modifications.Count - 1; index++)
                {
                    if (Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.DynamicMod &&
                        Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.UnknownType)
                    {
                        continue;
                    }

                    if (considerTargetResidues && string.IsNullOrWhiteSpace(Modifications[index].TargetResidues) || !considerTargetResidues)
                    {
                        if (Modifications[index].ModificationSymbol == modificationSymbol)
                        {
                            // Matching mass found
                            existingModFound = true;
                            return Modifications[index];
                        }
                    }
                }

                if (considerTargetResidues)
                {
                    // No match; try again, but ignore .TargetResidues
                    considerTargetResidues = false;
                }
                else
                {
                    break;
                }
            }

            // Still no match; return a default modification with a mass of 0
            return new ModificationDefinition(modificationSymbol, 0)
            {
                MassCorrectionTag = LookupMassCorrectionTagByMass(0)
            };
        }

        /// <summary>
        /// Looks for an existing modification with the given modification mass and target residues
        /// </summary>
        /// <remarks>If chTargetResidue is '\0' (which is default(char)) or AminoAcidModInfo.ResidueTerminusState.None, follows similar matching logic, but skips defined modifications with defined .TargetResidues</remarks>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="chTargetResidue">
        /// If defined, returns the first modification with the given mass and containing the residue in .TargetResidues;
        /// if no match, looks for the first modification with the given mass and no defined .TargetResidues
        /// </param>
        /// <param name="residueTerminusState">Residue terminus state</param>
        /// <param name="existingModFound">Output: true if an existing modification was found</param>
        /// <param name="addToModificationListIfUnknown">When true, add to the modification list if unknown</param>
        /// <param name="massDigitsOfPrecision">Number of digits after the decimal point to round to when comparing mod masses</param>
        /// <param name="massDigitsOfPrecisionLoose">Number of digits after the decimal point to round to, for a more lenient match (if no match found using massDigitsOfPrecision)</param>
        /// <returns>The best matched modification; if no match is found, returns a newly created modification definition, adding it to mModifications if addToModificationListIfUnknown is True</returns>
        public ModificationDefinition LookupModificationDefinitionByMass(
            double modificationMass,
            char chTargetResidue,
            AminoAcidModInfo.ResidueTerminusState residueTerminusState,
            out bool existingModFound,
            bool addToModificationListIfUnknown,
            byte massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION,
            byte massDigitsOfPrecisionLoose = MASS_DIGITS_OF_PRECISION)
        {
            existingModFound = false;

            var matchedMods = new List<KeyValuePair<double, ModificationDefinition>>();

            if (chTargetResidue != '\0' || residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
            {
                // The residue was provided and/or the residue is located at a peptide or protein terminus
                // First compare against modifications with 1 or more residues in .TargetResidues
                // On the first iteration, ignore terminal peptide and terminal protein mods
                // On the second iteration, match terminal peptide and terminal protein mods

                for (var modTypeIteration = 0; modTypeIteration < 2; modTypeIteration++)
                {
                    for (var index = 0; index <= Modifications.Count - 1; index++)
                    {
                        if (Modifications[index].TargetResidues.Length == 0)
                            continue;

                        // ReSharper disable once ConvertIfStatementToSwitchStatement
                        if (modTypeIteration == 0 &&
                            Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.DynamicMod &&
                            Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.StaticMod &&
                            Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.UnknownType)
                        {
                            continue;
                        }

                        if (modTypeIteration == 1 &&
                            Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod &&
                            Modifications[index].ModificationType != ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod)
                        {
                            continue;
                        }

                        // Absolute value of the mass difference
                        var deltaMassAbs = Math.Abs(Modifications[index].ModificationMass - modificationMass);

                        if (Math.Abs(Math.Round(deltaMassAbs, massDigitsOfPrecision)) > float.Epsilon)
                            continue;

                        // Matching mass found
                        // Now see if .TargetResidues contains chTargetResidue
                        if (Modifications[index].TargetResiduesContain(chTargetResidue))
                        {
                            existingModFound = true;
                        }

                        if (!existingModFound && residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
                        {
                            // ReSharper disable once CommentTypo

                            // Note that MaxQuant allows N-terminal protein mods to occur at either the true N-terminal residue, or one residue to the right
                            // For this reason, the _ModSummary.txt file will have entries like this:
                            // Symbol  Mass       TargetResidues  ModType  MassCorrectionTag  Occurrence_Count
                            // *       42.010567  [AMTSC          D        Acetyl             70

                            switch (residueTerminusState)
                            {
                                case AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus:
                                    if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        existingModFound = true;
                                    }

                                    break;

                                case AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus:
                                    if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        existingModFound = true;
                                    }

                                    break;
                            }
                        }

                        if (existingModFound)
                        {
                            // Match found; add to the list of candidate mods
                            matchedMods.Add(new KeyValuePair<double, ModificationDefinition>(deltaMassAbs, Modifications[index]));
                        }
                    }

                    if (matchedMods.Count > 0)
                    {
                        return FindClosestMatchedMod(matchedMods, chTargetResidue);
                    }
                }
            }

            // No match was found
            // Compare against modifications with empty .TargetResidues
            for (var index = 0; index <= Modifications.Count - 1; index++)
            {
                // ReSharper disable once InvertIf
                if (Modifications[index].ModificationType is ModificationDefinition.ResidueModificationType.DynamicMod or ModificationDefinition.ResidueModificationType.StaticMod or ModificationDefinition.ResidueModificationType.UnknownType &&
                    string.IsNullOrWhiteSpace(Modifications[index].TargetResidues))
                {
                    // Absolute value of the mass difference
                    var deltaMassAbs = Math.Abs(Modifications[index].ModificationMass - modificationMass);

                    if (Math.Abs(Math.Round(deltaMassAbs, massDigitsOfPrecision)) > float.Epsilon)
                        continue;

                    // Match found; add to the list of candidate mods
                    matchedMods.Add(new KeyValuePair<double, ModificationDefinition>(deltaMassAbs, Modifications[index]));
                }
            }

            if (matchedMods.Count > 0)
            {
                return FindClosestMatchedMod(matchedMods, chTargetResidue);
            }

            // Still no match; look for the modification mass and residue in mStandardRefinementModifications
            // Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing
            if (chTargetResidue != '\0')
            {
                for (var index = 0; index <= mStandardRefinementModifications.Count - 1; index++)
                {
                    // Absolute value of the mass difference
                    var deltaMassAbs = Math.Abs(mStandardRefinementModifications[index].ModificationMass - modificationMass);

                    if (Math.Abs(Math.Round(deltaMassAbs, massDigitsOfPrecision)) > float.Epsilon)
                        continue;

                    // Matching mass found
                    // Now see if .TargetResidues contains chTargetResidue
                    if (!mStandardRefinementModifications[index].TargetResiduesContain(chTargetResidue))
                        continue;

                    existingModFound = true;

                    var modificationDefinition = mStandardRefinementModifications[index];
                    modificationDefinition.ModificationSymbol = ModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL;

                    if (addToModificationListIfUnknown && mDefaultModificationSymbols.Count > 0)
                    {
                        // Append modificationDefinition to mModifications()
                        var newModIndex = AddModification(modificationDefinition, true);
                        return newModIndex >= 0 ? Modifications[newModIndex] : modificationDefinition;
                    }

                    // Either addToModificationListIfUnknown = False or no more default modification symbols
                    // Return modificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
                    return modificationDefinition;
                }
            }

            // Still no match
            // Compare against dynamic and unknown-type modifications, but ignore .TargetResidues
            for (var index = 0; index <= Modifications.Count - 1; index++)
            {
                if (Modifications[index].ModificationType is ModificationDefinition.ResidueModificationType.DynamicMod or ModificationDefinition.ResidueModificationType.UnknownType)
                {
                    var deltaMassAbs = Math.Abs(Modifications[index].ModificationMass - modificationMass);

                    if (Math.Abs(Math.Round(deltaMassAbs, massDigitsOfPrecision)) > float.Epsilon)
                        continue;

                    // Match found; add to the list of candidate mods
                    matchedMods.Add(new KeyValuePair<double, ModificationDefinition>(deltaMassAbs, Modifications[index]));
                }
            }

            if (matchedMods.Count > 0)
            {
                return FindClosestMatchedMod(matchedMods, chTargetResidue, true);
            }

            // Still no match; define a new custom modification
            const ModificationDefinition.ResidueModificationType modType = ModificationDefinition.ResidueModificationType.DynamicMod;
            const char modSymbol = ModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL;
            const bool useNextAvailableModificationSymbol = true;

            return AddUnknownModification(
                modificationMass, modType, chTargetResidue, residueTerminusState,
                addToModificationListIfUnknown, useNextAvailableModificationSymbol, modSymbol,
                massDigitsOfPrecision, massDigitsOfPrecisionLoose);
        }

        /// <summary>
        /// Resolve a modification name to a modification mass
        /// Checks both names in mMassCorrectionTags and UniMod names
        /// </summary>
        /// <param name="modName">Modification name</param>
        /// <param name="modMass">Monoisotopic mass</param>
        /// <returns>True if found, otherwise false</returns>
        public bool LookupModificationMassByName(string modName, out double modMass)
        {
            foreach (var item in mMassCorrectionTags)
            {
                if (string.Equals(item.Key, modName, StringComparison.OrdinalIgnoreCase))
                {
                    modMass = item.Value;
                    return true;
                }
            }

            foreach (var item in Modifications)
            {
                if (item.UniModName?.Equals(modName, StringComparison.OrdinalIgnoreCase) == true ||
                    item.MassCorrectionTag?.Equals(modName, StringComparison.OrdinalIgnoreCase) == true)
                {
                    modMass = item.ModificationMass;
                    return true;
                }
            }

            switch (modName.ToLower())
            {
                case "carbamyl-n":          // Used in TopPIC parameter file TopPIC_18Mods_15ppmParTol_0.8DaClusters_NumShift1_FASTA_has_decoy_EValue_0.5_2022-09-08.txt
                case "carbamyl_n":          // Used in TopPIC parameter file TopPIC_12Mods_15ppmParTol_0.8DaClusters_NumShift1_FASTA_has_decoy_EValue_0.5_2022-09-07.txt
                    modMass = 43.005814;
                    return true;

                case "deamidated":
                    modMass = 0.984016;
                    return true;

                case "dimethylation":
                    modMass = 28.0313;
                    return true;

                case "disulfide":           // Didehydro
                    modMass = -2.01565;
                    return true;

                case "methyl":
                case "methylation":
                    modMass = 14.01565;
                    return true;

                case "oxidation":
                    modMass = 15.994915;
                    return true;

                case "acetyl":
                    modMass = 42.010567;
                    return true;

                case "phospho":
                case "phosphorylation":
                    modMass = 79.966331;
                    return true;

                default:
                    modMass = 0;
                    return false;
            }
        }

        /// <summary>
        /// Looks for an existing modification with the given modification mass, modification type, and target residues
        /// </summary>
        /// <remarks>If chTargetResidue is '\0' (which is default(char)) or AminoAcidModInfo.ResidueTerminusState.None, follows similar matching logic, but skips defined modifications with defined .TargetResidues</remarks>
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="modType">Modification type</param>
        /// <param name="chTargetResidue">
        /// If defined, returns the first modification with the given mass and containing the residue in .TargetResidues;
        /// if no match, looks for the first modification with the given mass and no defined .TargetResidues</param>
        /// <param name="residueTerminusState">Residue terminus state</param>
        /// <param name="existingModFound">Output: true if an existing modification was found</param>
        /// <param name="addToModificationListIfUnknown">When true, add to the modification list if unknown</param>
        /// <param name="massDigitsOfPrecision">Number of digits after the decimal point to round to when comparing mod masses</param>
        /// <param name="massDigitsOfPrecisionLoose">Number of digits after the decimal point to round to, for a more lenient match (if no match found using massDigitsOfPrecision)</param>
        /// <returns>The best matched modification; if no match is found, returns a newly created modification definition, adding it to mModifications if addToModificationListIfUnknown = True</returns>
        public ModificationDefinition LookupModificationDefinitionByMassAndModType(
            double modificationMass,
            ModificationDefinition.ResidueModificationType modType,
            char chTargetResidue,
            AminoAcidModInfo.ResidueTerminusState residueTerminusState,
            out bool existingModFound,
            bool addToModificationListIfUnknown,
            byte massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION,
            byte massDigitsOfPrecisionLoose = MASS_DIGITS_OF_PRECISION)
        {
            // If chTargetResidue is defined, returns the first modification with the given mass and containing the residue in .TargetResidues
            // If no match is found, looks for the first modification with the given mass and no defined .TargetResidues
            // If no match is found, looks for the first dynamic modification with the given mass, regardless of .TargetResidues
            // If no match is found, returns a newly created modification definition, adding it to mModifications if addToModificationListIfUnknown = True
            // If chTargetResidue is '\0' or AminoAcidModInfo.ResidueTerminusState.None, follows similar logic, but skips defined modifications with defined .TargetResidues

            char modSymbol;
            bool useNextAvailableModificationSymbol;

            switch (modType)
            {
                case ModificationDefinition.ResidueModificationType.StaticMod:
                case ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod:
                case ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod:
                    modSymbol = ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                    useNextAvailableModificationSymbol = false;
                    break;

                default:
                    modSymbol = ModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL;
                    useNextAvailableModificationSymbol = true;
                    break;
            }

            existingModFound = false;
            if (chTargetResidue != '\0' || residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
            {
                // The residue was provided and/or the residue is located at a peptide or protein terminus
                // First compare against modifications with 1 or more residues in .TargetResidues
                for (var index = 0; index <= Modifications.Count - 1; index++)
                {
                    if (Modifications[index].ModificationType == modType &&
                        Modifications[index].ModificationType is ModificationDefinition.ResidueModificationType.DynamicMod or ModificationDefinition.ResidueModificationType.StaticMod or ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod or ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod or ModificationDefinition.ResidueModificationType.UnknownType &&
                        Modifications[index].TargetResidues.Length > 0)
                    {
                        if (Math.Round(Math.Abs(Modifications[index].ModificationMass - modificationMass), massDigitsOfPrecision) > float.Epsilon)
                            continue;

                        // Matching mass found
                        // Now see if .TargetResidues contains chTargetResidue
                        if (chTargetResidue != '\0' && Modifications[index].TargetResiduesContain(chTargetResidue))
                        {
                            existingModFound = true;
                        }

                        if (!existingModFound && residueTerminusState != AminoAcidModInfo.ResidueTerminusState.None)
                        {
                            switch (residueTerminusState)
                            {
                                case AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinNandCCTerminus:
                                    if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        existingModFound = true;
                                    }
                                    break;

                                case AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus:
                                case AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus:
                                    if (Modifications[index].TargetResiduesContain(AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS))
                                    {
                                        existingModFound = true;
                                    }
                                    break;
                            }
                        }

                        if (existingModFound)
                        {
                            // Match found
                            return Modifications[index];
                        }
                    }
                }
            }

            // No match was found
            // Compare against modifications with empty .TargetResidues
            for (var index = 0; index <= Modifications.Count - 1; index++)
            {
                if (Modifications[index].ModificationType == modType &&
                    Modifications[index].ModificationType is ModificationDefinition.ResidueModificationType.DynamicMod or ModificationDefinition.ResidueModificationType.StaticMod or ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod or ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod or ModificationDefinition.ResidueModificationType.UnknownType &&
                    string.IsNullOrWhiteSpace(Modifications[index].TargetResidues))
                {
                    if (Math.Round(Math.Abs(Modifications[index].ModificationMass - modificationMass), massDigitsOfPrecision) < float.Epsilon)
                    {
                        // Matching mass found
                        return Modifications[index];
                    }
                }
            }

            // Still no match; look for the modification mass and residue in mStandardRefinementModifications
            // Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing or chTargetResidue = '<' or chTargetResidue = '>'
            if (chTargetResidue != '\0')
            {
                for (var index = 0; index <= mStandardRefinementModifications.Count - 1; index++)
                {
                    if (Math.Round(Math.Abs(mStandardRefinementModifications[index].ModificationMass - modificationMass), massDigitsOfPrecision) > float.Epsilon)
                        continue;

                    // Matching mass found
                    // Now see if .TargetResidues contains chTargetResidue
                    if (!mStandardRefinementModifications[index].TargetResiduesContain(chTargetResidue))
                        continue;

                    existingModFound = true;

                    var modificationDefinition = mStandardRefinementModifications[index];
                    modificationDefinition.ModificationSymbol = modSymbol;
                    modificationDefinition.ModificationType = modType;

                    if (!addToModificationListIfUnknown || mDefaultModificationSymbols.Count <= 0)
                        return modificationDefinition;

                    // Append modificationDefinition to mModifications()
                    var newModIndex = AddModification(modificationDefinition, true);

                    if (newModIndex >= 0)
                    {
                        return Modifications[newModIndex];
                    }

                    // Either addToModificationListIfUnknown = False or no more default modification symbols
                    // Return modificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
                    return modificationDefinition;
                }
            }

            // No match was found
            // Compare against modifications of the same type, but ignore .TargetResidues
            for (var index = 0; index <= Modifications.Count - 1; index++)
            {
                if (Modifications[index].ModificationType != modType)
                    continue;

                if (Math.Round(Math.Abs(Modifications[index].ModificationMass - modificationMass), massDigitsOfPrecision) > float.Epsilon)
                    continue;

                // Matching mass found
                // Assure that the target residues contain chTargetResidue
                if (chTargetResidue != '\0' && !Modifications[index].TargetResiduesContain(chTargetResidue))
                {
                    Modifications[index].TargetResidues += chTargetResidue;
                }

                return Modifications[index];
            }

            // Still no match; define a new custom modification
            return AddUnknownModification(
                modificationMass, modType, chTargetResidue, residueTerminusState,
                addToModificationListIfUnknown, useNextAvailableModificationSymbol, modSymbol,
                massDigitsOfPrecision, massDigitsOfPrecisionLoose);
        }

        private void InitializeLocalVariables()
        {
            ErrorMessage = string.Empty;
            SetDefaultMassCorrectionTags();

            // Note that this method will call UpdateDefaultModificationSymbols()
            ClearModifications();

            UpdateStandardRefinementModifications();

            UpdateIntegerBasedModificationMap();
        }

        /// <summary>
        /// Load the mass correction tags file
        /// </summary>
        /// <param name="filePath">Mass correction tags file path</param>
        /// <param name="fileNotFound">Output: true if the file was not found</param>
        /// <returns>True if successful, false if an error</returns>
        public bool ReadMassCorrectionTagsFile(string filePath, out bool fileNotFound)
        {
            fileNotFound = false;

            try
            {
                // Open the mass correction tags file
                // Prior to September 2021, it did not have a header line

                // It should have 3 columns, separated by tabs
                // Column 1 is the mass correction tag name
                // Column 2 is the monoisotopic mass for the mass correction (positive or negative number)
                // Column 3 is the affected atom (only used for isotopic mods)

                if (string.IsNullOrWhiteSpace(filePath))
                {
                    SetDefaultMassCorrectionTags();
                    return true;
                }

                if (!File.Exists(filePath))
                {
                    ErrorMessage = "Mass CorrectionTags File Not Found: " + filePath;
                    SetDefaultMassCorrectionTags();
                    fileNotFound = true;
                    return false;
                }

                using var massCorrectionTagsReader = new StreamReader(filePath);

                mMassCorrectionTags.Clear();

                while (!massCorrectionTagsReader.EndOfStream)
                {
                    var dataLine = massCorrectionTagsReader.ReadLine();

                    if (string.IsNullOrEmpty(dataLine) || dataLine.StartsWith("Mass_Correction_Tag"))
                        continue;

                    var splitLine = dataLine.Split('\t');

                    if (splitLine.Length < 2)
                    {
                        continue;
                    }

                    // See if the first column contains 1 or more characters and if the second column contains a number
                    // Note that StoreMassCorrectionTag() will trim spaces from the end of the mass correction tag names
                    if (splitLine[0].Trim().Length >= 1 && SynFileReaderBaseClass.IsNumber(splitLine[1]))
                    {
                        StoreMassCorrectionTag(splitLine[0], double.Parse(splitLine[1]));
                    }
                }

                if (mMassCorrectionTags.Count == 0)
                {
                    SetDefaultMassCorrectionTags();
                }

                return true;
            }
            catch (Exception ex)
            {
                ErrorMessage = "Error reading Mass Correction Tags file (" + filePath + "): " + ex.Message;
                return false;
            }
        }

        /// <summary>
        /// Read a modification definitions file (_ModDefs.txt)
        /// </summary>
        /// <param name="filePath">Modification definitions file path</param>
        /// <param name="fileNotFound">Output: true if the file was not found</param>
        /// <returns>True if successful, false if an error</returns>
        public bool ReadModificationDefinitionsFile(string filePath, out bool fileNotFound)
        {
            fileNotFound = false;

            try
            {
                // Open the modification file

                // It should have 2 or more columns, separated by tabs
                // It may or may not have a header line

                // When no header line is defined, assume a fixed column order:
                // Column 1: modification symbol
                // Column 2: modification mass
                // Column 3: residues and/or termini that can be modified; if omitted, the modification can apply to any residues or termini
                //   For column 3, use 1 letter amino acid abbreviations; the residues can be a continuous string, or can be separated by commas and/or spaces
                //   For column 3, use the *_SYMBOL_DMS constants for the termini (< and > for the peptide termini; [ and ] for the protein termini)
                // Column 4: type of modification: D, S, T, I, or P (corresponding to ModificationDefinition.ResidueModificationType)
                // Column 5: DMS-specific mass correction tag associated with the given modification
                // Column 6: affected atom

                // If a header line is present, column 6 could be UniMod_Mod_Name
                // For MaxQuant parameter files, column 7 will be MaxQuant_Mod_Name

                if (string.IsNullOrWhiteSpace(filePath))
                {
                    ClearModifications();
                    return true;
                }

                if (!File.Exists(filePath))
                {
                    ErrorMessage = "Modification Definition File Not Found: " + filePath;
                    ClearModifications();
                    fileNotFound = true;
                    return false;
                }

                ClearModifications();

                var headersDefined = false;
                var headerNamesFound = false;

                // Initialize the column mapping object
                var columnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);

                const string MOD_SUMMARY_COLUMN_Monoisotopic_Mass = "Monoisotopic_Mass";
                const string MOD_SUMMARY_COLUMN_Affected_Atom = "Affected_Atom";
                const string MOD_SUMMARY_COLUMN_UniMod_Mod_Name = "UniMod_Mod_Name";
                const string MOD_SUMMARY_COLUMN_MaxQuant_Mod_Name = "MaxQuant_Mod_Name";

                var headerNames = new List<string>
                {
                    PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol,
                    MOD_SUMMARY_COLUMN_Monoisotopic_Mass,
                    PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Target_Residues,
                    PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Type,
                    PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Mass_Correction_Tag,
                    MOD_SUMMARY_COLUMN_Affected_Atom
                };

                using var modFileReader = new StreamReader(filePath);

                while (!modFileReader.EndOfStream)
                {
                    var dataLine = modFileReader.ReadLine();

                    if (string.IsNullOrEmpty(dataLine))
                        continue;

                    var columns = dataLine.Split('\t');

                    if (!headersDefined)
                    {
                        headersDefined = true;

                        if (dataLine.StartsWith(PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol))
                        {
                            // Add the optional header names
                            headerNames.Add(MOD_SUMMARY_COLUMN_UniMod_Mod_Name);
                            headerNames.Add(MOD_SUMMARY_COLUMN_MaxQuant_Mod_Name);
                            headerNames.Add(PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Occurrence_Count);

                            // _ModDefs.txt files    have Monoisotopic_Mass as the second column
                            // _ModSummary.txt files have Modification_Mass as the second column
                            // Add "Modification_Mass" just in case somebody manually copied headers over
                            headerNames.Add(PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass);

                            SynFileReaderBaseClass.DefineColumnHeaders(columnHeaders, headerNames);

                            ReaderFactory.ParseColumnHeaders(columns, columnHeaders);

                            if (columnHeaders[MOD_SUMMARY_COLUMN_Monoisotopic_Mass] < 0 &&
                                columnHeaders[PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass] >= 0)
                            {
                                // Update the column index for the monoisotopic mass column
                                columnHeaders[MOD_SUMMARY_COLUMN_Monoisotopic_Mass] = columnHeaders[PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Mass];
                            }

                            headerNamesFound = true;
                            continue;
                        }

                        // Assume the default header order
                        SynFileReaderBaseClass.DefineColumnHeaders(columnHeaders, headerNames);
                    }

                    if (!headerNamesFound)
                    {
                        if (columns.Length < 2)
                            continue;

                        // See if the first column contains a single character and if the second column contains a number
                        if (columns[0].Trim().Length != 1 || !SynFileReaderBaseClass.IsNumber(columns[1]))
                            continue;
                    }

                    var modificationSymbol = ReaderFactory.LookupColumnValue(columns, PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Symbol, columnHeaders);
                    var modificationMass = ReaderFactory.LookupColumnValue(columns, MOD_SUMMARY_COLUMN_Monoisotopic_Mass, columnHeaders, 0.0);

                    var modificationDefinition = new ModificationDefinition(modificationSymbol.Trim()[0], modificationMass);

                    if (ReaderFactory.TryGetColumnValue(
                        columns, PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Target_Residues, columnHeaders,
                        out var targetResidues))
                    {
                        // Parse the target residues list
                        var residues = targetResidues.Trim().ToUpper();

                        var residuesClean = string.Empty;

                        foreach (var residue in residues)
                        {
                            if (char.IsUpper(residue))
                            {
                                residuesClean += residue;
                            }
                            else if (residue is AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS or AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS or
                                AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS or AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)
                            {
                                residuesClean += residue;
                            }
                        }

                        if (residuesClean.Length > 0)
                        {
                            modificationDefinition.TargetResidues = residuesClean;
                        }
                    }

                    if (ReaderFactory.TryGetColumnValue(
                        columns, PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Modification_Type, columnHeaders,
                        out var modificationType))
                    {
                        // Store the modification type
                        if (modificationType.Trim().Length == 1)
                        {
                            var modTypeSymbol = modificationType.ToUpper().Trim()[0];
                            modificationDefinition.ModificationType = ModificationDefinition.ModificationSymbolToModificationType(modTypeSymbol);
                        }

                        // If the .ModificationType is unknown, change it to Dynamic
                        if (modificationDefinition.ModificationType == ModificationDefinition.ResidueModificationType.UnknownType)
                        {
                            modificationDefinition.ModificationType = ModificationDefinition.ResidueModificationType.DynamicMod;
                        }
                    }

                    if (ReaderFactory.TryGetColumnValue(
                        columns, PHRPModSummaryReader.MOD_SUMMARY_COLUMN_Mass_Correction_Tag, columnHeaders,
                        out var massCorrectionTag))
                    {
                        modificationDefinition.MassCorrectionTag = massCorrectionTag.Trim();
                    }

                    if (ReaderFactory.TryGetColumnValue(
                        columns, MOD_SUMMARY_COLUMN_Affected_Atom, columnHeaders,
                        out var affectedAtom))
                    {
                        if (affectedAtom.Trim().Length > 0)
                        {
                            modificationDefinition.AffectedAtom = affectedAtom[0];
                        }
                        else
                        {
                            modificationDefinition.AffectedAtom = PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL;
                        }
                    }

                    if (ReaderFactory.TryGetColumnValue(
                        columns, MOD_SUMMARY_COLUMN_UniMod_Mod_Name, columnHeaders,
                        out var uniModName))
                    {
                        modificationDefinition.UniModName = uniModName.Trim();
                    }

                    if (ReaderFactory.TryGetColumnValue(
                        columns, MOD_SUMMARY_COLUMN_MaxQuant_Mod_Name, columnHeaders,
                        out var maxQuantMod))
                    {
                        modificationDefinition.MaxQuantModName= maxQuantMod.Trim();
                    }

                    // Check whether the modification type is Static and the .TargetResidues are one of: <>[]
                    // If so, update the modification type as needed
                    if (modificationDefinition.TargetResidues?.Trim().Length == 1 && modificationDefinition.ModificationType == ModificationDefinition.ResidueModificationType.StaticMod)
                    {
                        if (modificationDefinition.TargetResidues[0] == AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS ||
                            modificationDefinition.TargetResidues[0] == AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                        {
                            modificationDefinition.ModificationType = ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod;
                        }
                        else if (modificationDefinition.TargetResidues[0] == AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS ||
                                 modificationDefinition.TargetResidues[0] == AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)
                        {
                            modificationDefinition.ModificationType = ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod;
                        }
                    }

                    // Validate some of the settings if the modification type is IsotopicMod or TerminalPeptideStaticMod or ProteinTerminusStaticMod
                    var validMod = true;

                    switch (modificationDefinition.ModificationType)
                    {
                        case ModificationDefinition.ResidueModificationType.IsotopicMod:
                            modificationDefinition.ModificationSymbol = ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                            if (modificationDefinition.AffectedAtom == PeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL)
                            {
                                validMod = false;
                            }
                            break;

                        case ModificationDefinition.ResidueModificationType.TerminalPeptideStaticMod:
                            modificationDefinition.ModificationSymbol = ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                            if (modificationDefinition.TargetResidues != AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString() &&
                                modificationDefinition.TargetResidues != AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                            {
                                validMod = false;
                            }
                            break;

                        case ModificationDefinition.ResidueModificationType.ProteinTerminusStaticMod:
                            modificationDefinition.ModificationSymbol = ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL;
                            if (modificationDefinition.TargetResidues != AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString() &&
                                modificationDefinition.TargetResidues != AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString())
                            {
                                validMod = false;
                            }
                            break;

                        case ModificationDefinition.ResidueModificationType.UnknownType:
                            modificationDefinition.ModificationType = ModificationDefinition.ResidueModificationType.DynamicMod;
                            break;
                    }

                    if (modificationDefinition.MassCorrectionTag == ModificationDefinition.INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME)
                    {
                        // Try to determine the mass correction name
                        modificationDefinition.MassCorrectionTag = LookupMassCorrectionTagByMass(modificationDefinition.ModificationMass);
                    }

                    if (validMod)
                    {
                        AddModification(modificationDefinition, false);
                    }
                }

                // Note that this method will call UpdateDefaultModificationSymbols()
                ValidateModificationsVsDefaultModificationSymbols();
                return true;
            }
            catch (Exception ex)
            {
                ErrorMessage = "Error reading Modification Definition file (" + filePath + "): " + ex.Message;
                return false;
            }
        }

        /// <summary>
        /// Reset mod occurrence count stats
        /// </summary>
        public void ResetOccurrenceCountStats()
        {
            for (var index = 0; index <= Modifications.Count - 1; index++)
            {
                Modifications[index].OccurrenceCount = 0;
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

                // ReSharper disable StringLiteralTypo
                StoreMassCorrectionTag("4xDeut  ", 4.025107);
                StoreMassCorrectionTag("6C134N15", 10.008269);      // HeavyR
                StoreMassCorrectionTag("6xC13N15", 7.017164);
                StoreMassCorrectionTag("6C132N15", 8.014199);       // HeavyK
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
                StoreMassCorrectionTag("HeavyK  ", 8.014199);
                StoreMassCorrectionTag("HeavyR  ", 10.008269);
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

                // ReSharper restore StringLiteralTypo
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        private void StoreMassCorrectionTag(string tagName, double mass)
        {
            try
            {
                mMassCorrectionTags.Add(tagName.Trim(), mass);
            }
            catch (Exception)
            {
                // If a duplicate is tag is entered into the mMassCorrectionTags Dictionary, an error will occur; we'll ignore the error
                // Ignore errors here
            }
        }

        private void UpdateDefaultModificationSymbols(string modificationChars)
        {
            try
            {
                if (string.IsNullOrEmpty(modificationChars))
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
                foreach (var chChar in modificationChars)
                {
                    if (chChar != ModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL &&
                        chChar != ModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL)
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
            // ReSharper disable StringLiteralTypo

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
            // ReSharper restore StringLiteralTypo

        }

        private void UpdateStandardRefinementModifications()
        {
            mStandardRefinementModifications = new List<ModificationDefinition>();

            const double NH3_LOSS = -17.026549;
            mStandardRefinementModifications.Add(new ModificationDefinition(
                ModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL,
                NH3_LOSS,
                "Q",
                ModificationDefinition.ResidueModificationType.DynamicMod,
                LookupMassCorrectionTagByMass(NH3_LOSS)));

            const double H2O_LOSS = -18.0106;
            mStandardRefinementModifications.Add(new ModificationDefinition(
                ModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL,
                H2O_LOSS,
                "E",
                ModificationDefinition.ResidueModificationType.DynamicMod,
                LookupMassCorrectionTagByMass(H2O_LOSS)));
        }

        private void ValidateModificationsVsDefaultModificationSymbols()
        {
            try
            {
                // Reset the default modification symbols list
                UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS);

                var chDefaultModificationSymbols = new char[mDefaultModificationSymbols.Count];
                mDefaultModificationSymbols.ToArray().CopyTo(chDefaultModificationSymbols, 0);
                var defaultModificationSymbolCount = chDefaultModificationSymbols.Length;

                // Step through mModifications and make sure each of the modification symbols is not present in mDefaultModificationChars
                for (var index = 0; index <= Modifications.Count - 1; index++)
                {
                    var indexCompare = 0;

                    while (indexCompare < defaultModificationSymbolCount)
                    {
                        if (Modifications[index].ModificationSymbol == chDefaultModificationSymbols[indexCompare])
                        {
                            // Remove this symbol from chDefaultModificationSymbols
                            for (var indexCopy = indexCompare; indexCopy <= defaultModificationSymbolCount - 2; indexCopy++)
                            {
                                chDefaultModificationSymbols[indexCopy] = chDefaultModificationSymbols[indexCopy + 1];
                            }
                            defaultModificationSymbolCount--;
                        }
                        else
                        {
                            indexCompare++;
                        }
                    }
                }

                if (defaultModificationSymbolCount < mDefaultModificationSymbols.Count)
                {
                    mDefaultModificationSymbols.Clear();
                    for (var index = 0; index <= defaultModificationSymbolCount - 1; index++)
                    {
                        mDefaultModificationSymbols.Enqueue(chDefaultModificationSymbols[index]);
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
        /// <param name="modificationMass">Modification mass</param>
        /// <param name="targetResidues">Target residues</param>
        /// <param name="modificationType">Modification type</param>
        /// <param name="massDigitsOfPrecision">Mass digits of precision</param>
        /// <returns>True if the modification was matched or was added; false if an error</returns>
        public bool VerifyModificationPresent(
            double modificationMass,
            string targetResidues,
            ModificationDefinition.ResidueModificationType modificationType,
            int massDigitsOfPrecision = MASS_DIGITS_OF_PRECISION)
        {
            // ReSharper disable once GrammarMistakeInComment

            // Look for mods in mModifications with matching .ModificationType, .ModificationMass (within tolerance), and .TargetResidues
            // If not found, add a new entry to mModifications

            var matchFound = false;

            if (massDigitsOfPrecision < 0)
                massDigitsOfPrecision = 0;

            try
            {
                for (var index = 0; index <= Modifications.Count - 1; index++)
                {
                    if (Modifications[index].ModificationType != modificationType)
                        continue;

                    // Matching modification type
                    if (Math.Abs(Math.Round(Math.Abs(Modifications[index].ModificationMass - modificationMass), massDigitsOfPrecision)) > float.Epsilon)
                        continue;

                    // Matching mass
                    // Compare .TargetResidues
                    matchFound = ModificationDefinition.EquivalentTargetResidues(Modifications[index].TargetResidues, targetResidues, true);
                    if (matchFound)
                    {
                        break;
                    }
                }

                if (!matchFound)
                {
                    var modificationDefinition = new ModificationDefinition(modificationMass, targetResidues, modificationType)
                    {
                        MassCorrectionTag = LookupMassCorrectionTagByMass(modificationMass)
                    };

                    // Append modificationDefinition to mModifications()
                    AddModification(modificationDefinition, true);

                    matchFound = true;
                }
            }
            catch (Exception)
            {
                matchFound = false;
            }

            return matchFound;
        }
    }
}
