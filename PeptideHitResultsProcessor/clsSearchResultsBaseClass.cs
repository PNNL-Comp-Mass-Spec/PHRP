// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 7, 2006
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
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class is used to track the peptide details for a given MS/MS search result.
    /// It can track peptide residue level modifications and static, peptide-wide modifications.
    /// </summary>
    /// <remarks>
    /// Use SearchResultClearModifications() and SearchResultAddModification()
    /// to track specific modifications for a given peptide
    /// </remarks>
    public abstract class clsSearchResultsBaseClass
    {
        #region "Constants and Enums"
        /// <summary>
        /// Mass difference, in Da, between ^12C and ^13C
        /// </summary>
        /// <remarks></remarks>
        public const double MASS_C13 = 1.00335483;
        #endregion

        #region "Class wide Variables"
        // Note: Many of these variables typically hold numbers but we're storing the numbers as strings
        //       to prevent the numeric representation from changing when converting to a number then back to a string

        /// <summary>
        /// Residue or residues before the start of the peptide sequence
        /// </summary>
        /// <remarks></remarks>
        private string mPeptidePreResidues;

        /// <summary>
        /// Residue or residues after the end of the peptide sequence
        /// </summary>
        /// <remarks></remarks>
        private string mPeptidePostResidues;

        /// <summary>
        /// Peptide sequence without any modification symbols
        /// </summary>
        /// <remarks></remarks>
        private string mPeptideCleanSequence;

        private clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants mPeptideCleavageState;

        /// <summary>
        /// Terminus state of the peptide
        /// </summary>
        /// <remarks></remarks>
        private clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants mPeptideTerminusState;

        /// <summary>
        /// List of modifications present in the current peptide
        /// </summary>
        /// <remarks></remarks>
        private readonly List<clsAminoAcidModInfo> mSearchResultModifications;

        /// <summary>
        /// Possible modifications that the peptide could have
        /// </summary>
        /// <remarks></remarks>
        private readonly clsPeptideModificationContainer mPeptideMods;

        private readonly clsPeptideCleavageStateCalculator mPeptideCleavageStateCalculator;
        private readonly clsPeptideMassCalculator mPeptideSeqMassCalculator;

        private string mErrorMessage = string.Empty;
        #endregion

        #region "Properties"

        /// <summary>
        /// Most recent error message
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ErrorMessage => mErrorMessage;

        /// <summary>
        /// RowIndex for Synopsis/First Hits files; auto-assigned for XTandem, Inspect, MSGFDB, and MODa
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int ResultID { get; set; }

        /// <summary>
        /// Group ID assigned by XTandem
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int GroupID { get; set; }

        /// <summary>
        /// Scan number
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string Scan { get; set; }

        /// <summary>
        /// Charge state
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string Charge { get; set; }

        /// <summary>
        /// Observed precursor m/z value converted to M+H
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ParentIonMH { get; set; }

        /// <summary>
        /// Multiple protein count: 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string MultipleProteinCount { get; set; }

        /// <summary>
        /// First protein for this PSM
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ProteinName { get; set; }

        /// <summary>
        /// Typically only used by XTandem; actually holds the Log of the expectation value
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ProteinExpectationValue { get; set; }

        /// <summary>
        /// Typically only used by XTandem; actually holds the Log of the intensity
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ProteinIntensity { get; set; }

        /// <summary>
        /// Number of the first residue of a protein; typically always 1
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Initialized to 0 (unknown) then later changed 1 when ProteinSeqResidueNumberEnd is defined </remarks>
        public int ProteinSeqResidueNumberStart { get; set; }

        /// <summary>
        /// The residue number of the last residue in the protein's sequence
        /// For example, 100 if the protein has 100 residues total
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>>Initialized to 0 (unknown)</remarks>
        public int ProteinSeqResidueNumberEnd { get; set; }

        /// <summary>
        /// Position in the protein's residues of the first residue in the peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int PeptideLocInProteinStart { get; set; }

        /// <summary>
        /// Position in the protein's residues of the last residue in the peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int PeptideLocInProteinEnd { get; set; }

        /// <summary>
        ///  Residue or residues before the start of the peptide sequence
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Calls ComputePeptideCleavageStateInProtein</remarks>
        public string PeptidePreResidues
        {
            get => mPeptidePreResidues;
            set
            {
                if (value == null)
                    value = string.Empty;
                mPeptidePreResidues = value;
                ComputePeptideCleavageStateInProtein();
            }
        }

        /// <summary>
        /// Residue or residues after the end of the peptide sequence
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Calls ComputePeptideCleavageStateInProtein</remarks>
        public string PeptidePostResidues
        {
            get => mPeptidePostResidues;
            set
            {
                if (value == null)
                    value = string.Empty;
                mPeptidePostResidues = value;
                ComputePeptideCleavageStateInProtein();
            }
        }

        /// <summary>
        /// Peptide sequence without any modification symbols
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Calls ComputePeptideCleavageStateInProtein</remarks>
        public string PeptideCleanSequence
        {
            get => mPeptideCleanSequence;
            set
            {
                mPeptideCleanSequence = value;
                ComputePeptideCleavageStateInProtein();
            }
        }

        /// <summary>
        /// Peptide sequence with modification symbols
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Mod symbols are single characters, like *, #, @, etc.</remarks>
        public string PeptideSequenceWithMods { get; set; }

        /// <summary>
        /// Cleavage state of the peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants PeptideCleavageState => mPeptideCleavageState;

        /// <summary>
        /// Terminus state of the peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants PeptideTerminusState => mPeptideTerminusState;

        /// <summary>
        /// In XTandem this is the theoretical monoisotopic MH
        /// In Sequest it was historically the average mass MH, though when a monoisotopic mass parent tolerance is specified, this is a monoisotopic mass
        /// In Inspect, MS-GF+, and MSAlign, this is the theoretical monoisotopic MH; note that this is (M+H)+
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string PeptideMH { get; set; }

        /// <summary>
        /// Difference in mass between the peptide's computed mass and the parent ion mass (i.e. the mass chosen for fragmentation)
        /// In Sequest this is Theoretical Mass - Observed Mass
        /// In XTandem, Inspect, MS-GF+, and MSAlign the DelM value is listed as Observed - Theoretical,
        /// however, PHRP negates that value while reading the synopsis file to match Sequest
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string PeptideDeltaMass { get; set; }

        /// <summary>
        /// Comma separated list of the modifications present and the residue modified, with a colon separating the mod name and residue location
        /// </summary>
        /// <remarks>
        /// Examples:
        ///   Acetyl:1
        ///   MinusH2O:1,Plus1Oxy:4,Plus1Oxy:19
        /// </remarks>
        public string PeptideModDescription { get; set; }

        /// <summary>
        /// Theoretical (computed) monoisotopic mass for a given peptide sequence, including any modified residues
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public double PeptideMonoisotopicMass { get; set; }

        /// <summary>
        /// Number of peptide residue modifications (dynamic or static)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int SearchResultModificationCount => mSearchResultModifications.Count;

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        /// <remarks></remarks>
        protected clsSearchResultsBaseClass(clsPeptideModificationContainer peptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
        {
            mSearchResultModifications = new List<clsAminoAcidModInfo>();

            mPeptideMods = peptideMods;

            mPeptideCleavageStateCalculator = new clsPeptideCleavageStateCalculator();
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin);

            mPeptideSeqMassCalculator = peptideSeqMassCalculator ?? throw new Exception("peptideSeqMassCalculator instance cannot be null");

            InitializeLocalVariables();
        }

        /// <summary>
        /// Add the modification symbols for the current peptide to the clean sequence
        /// </summary>
        /// <returns>Sequence with mod symbols</returns>
        /// <remarks></remarks>
        public string AddSearchResultModificationsToCleanSequence()
        {
            // Initialize sequenceWithMods to the clean sequence; we'll insert the mod symbols below if mSearchResultModifications.Count > 0
            var sequenceWithMods = string.Copy(mPeptideCleanSequence);

            if (mSearchResultModifications.Count > 0)
            {
                // Insert the modification symbols into sequenceWithMods
                // First, sort mSearchResultModifications on .ResidueLocInPeptide and .MassCorrectionTag

                if (mSearchResultModifications.Count > 1)
                {
                    mSearchResultModifications.Sort(new IGenericResidueModificationInfoComparer());
                }

                // Now step backward through residueModificationPositions and add the symbols to sequenceWithMods
                for (var index = mSearchResultModifications.Count - 1; index >= 0; index += -1)
                {
                    var resultMod = mSearchResultModifications[index];
                    if (resultMod.ModDefinition.ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod ||
                        resultMod.ModDefinition.ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType)
                    {
                        sequenceWithMods = sequenceWithMods.Insert(resultMod.ResidueLocInPeptide, resultMod.ModDefinition.ModificationSymbol.ToString());
                    }
                }
            }

            return sequenceWithMods;
        }

        /// <summary>
        /// Update .PeptideSequenceWithMods and .PeptideModDescription using the modifications defined for this peptide
        /// </summary>
        /// <remarks></remarks>
        public void ApplyModificationInformation()
        {
            PeptideSequenceWithMods = AddSearchResultModificationsToCleanSequence();
            UpdateModDescription();
        }

        /// <summary>
        /// Clear all values
        /// </summary>
        /// <remarks></remarks>
        public virtual void Clear()
        {
            ResultID = 0;
            GroupID = 0;

            Scan = string.Empty;
            Charge = string.Empty;
            ParentIonMH = string.Empty;

            MultipleProteinCount = string.Empty;
            ProteinName = string.Empty;
            ProteinExpectationValue = string.Empty;
            ProteinIntensity = string.Empty;

            ClearProteinSequenceInfo();

            // Note that this calls ClearSearchResultModifications()
            ClearPeptideDetailsInfo();
        }

        /// <summary>
        /// Clear cached peptide information
        /// </summary>
        /// <remarks></remarks>
        public void ClearPeptideDetailsInfo()
        {
            PeptideLocInProteinStart = 0;
            PeptideLocInProteinEnd = 0;

            mPeptidePreResidues = string.Empty;
            mPeptidePostResidues = string.Empty;
            mPeptideCleanSequence = string.Empty;
            PeptideSequenceWithMods = string.Empty;

            mPeptideCleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific;
            mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None;

            PeptideMH = string.Empty;
            PeptideDeltaMass = string.Empty;

            PeptideModDescription = string.Empty;
            PeptideMonoisotopicMass = 0;

            ClearSearchResultModifications();
        }

        /// <summary>
        /// Clear the Protein sequence start residue and end residue values (ProteinSeqResidueNumberStart and ProteinSeqResidueNumberEnd)
        /// </summary>
        /// <remarks></remarks>
        public void ClearProteinSequenceInfo()
        {
            ProteinSeqResidueNumberStart = 0;
            ProteinSeqResidueNumberEnd = 0;
        }

        /// <summary>
        /// Clear the modifications for this peptide
        /// </summary>
        /// <remarks></remarks>
        public void ClearSearchResultModifications()
        {
            mSearchResultModifications.Clear();
        }

        /// <summary>
        /// Compute the delta mass, in ppm, optionally correcting for C13 isotopic selection errors
        /// </summary>
        /// <param name="delM">Delta mass, in Da</param>
        /// <param name="precursorMonoMass">Precursor monoisotopic mass</param>
        /// <param name="adjustPrecursorMassForC13">True to correct for C13 isotopic selection errors</param>
        /// <param name="peptideMonoisotopicMass"></param>
        /// <returns>Delta mass, in ppm</returns>
        /// <remarks></remarks>
        public static double ComputeDelMCorrectedPPM(
            double delM,
            double precursorMonoMass,
            bool adjustPrecursorMassForC13,
            double peptideMonoisotopicMass)
        {
            var correctionCount = 0;

            // Examine delM to determine which isotope was chosen
            if (delM >= -0.5)
            {
                // This is the typical case
                while (delM > 0.5)
                {
                    delM -= MASS_C13;
                    correctionCount += 1;
                }
            }
            else
            {
                // This happens less often; but we'll still account for it
                // In this case, correctionCount will be negative
                while (delM < -0.5)
                {
                    delM += MASS_C13;
                    correctionCount -= 1;
                }
            }

            if (correctionCount != 0)
            {
                if (adjustPrecursorMassForC13)
                {
                    // Adjust the precursor mono mass based on correctionCount
                    precursorMonoMass -= correctionCount * MASS_C13;
                }

                // Compute a new DelM value
                delM = precursorMonoMass - peptideMonoisotopicMass;
            }

            return clsPeptideMassCalculator.MassToPPM(delM, peptideMonoisotopicMass);
        }

        /// <summary>
        /// Computes the theoretical monoisotopic mass for the given peptide
        /// Also updates mPeptideDeltaMassCorrectedPpm
        /// </summary>
        /// <remarks></remarks>
        public void ComputeMonoisotopicMass()
        {
            var modifiedResidues = new List<clsPeptideMassCalculator.udtPeptideSequenceModInfoType>();

            // Copy the mod info from mSearchResultModifications to list modifiedResidues
            foreach (var searchResultMod in mSearchResultModifications)
            {
                var modifiedResidue =
                    new clsPeptideMassCalculator.udtPeptideSequenceModInfoType
                    {
                        ResidueLocInPeptide = searchResultMod.ResidueLocInPeptide,
                        ModificationMass = searchResultMod.ModDefinition.ModificationMass,
                        AffectedAtom = searchResultMod.ModDefinition.AffectedAtom
                    };

                modifiedResidues.Add(modifiedResidue);
            }

            PeptideMonoisotopicMass = mPeptideSeqMassCalculator.ComputeSequenceMass(mPeptideCleanSequence, modifiedResidues);
        }

        /// <summary>
        /// Update PeptideCleavageState and PeptideTerminusState based on PeptideCleanSequence, PeptidePreResidues and PeptidePostResidues
        /// </summary>
        /// <remarks></remarks>
        public void ComputePeptideCleavageStateInProtein()
        {
            // Determine the peptide's terminus state and cleavage state within the protein
            mPeptideCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues);
            mPeptideTerminusState = mPeptideCleavageStateCalculator.ComputeTerminusState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues);
        }

        /// <summary>
        /// Determine the terminus state for the given residue in the peptide
        /// </summary>
        /// <param name="residueLocInPeptide">Residue number (1-based)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsAminoAcidModInfo.eResidueTerminusStateConstants DetermineResidueTerminusState(int residueLocInPeptide)
        {
            var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
            if (residueLocInPeptide == 1)
            {
                // Residue is at the peptide's N-terminus
                if (PeptideLocInProteinStart == ProteinSeqResidueNumberStart)
                {
                    // Residue is at the protein's N-terminus
                    if (PeptideLocInProteinEnd == ProteinSeqResidueNumberEnd)
                    {
                        // The protein only contains one Residue, and we're currently examining it
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus;
                    }
                    else
                    {
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                    }
                }
                else
                {
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                }
            }
            else
            {
                if (residueLocInPeptide == PeptideLocInProteinEnd - PeptideLocInProteinStart + 1)
                {
                    // Residue is at the peptide's C-terminus
                    if (PeptideLocInProteinEnd == ProteinSeqResidueNumberEnd)
                    {
                        // Residue is at the protein's C-terminus
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                    }
                    else
                    {
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                    }
                }
            }

            return eResidueTerminusState;
        }

        /// <summary>
        /// Get the modification info, by index
        /// </summary>
        /// <param name="index">Modification entry index (0-based)</param>
        /// <param name="residue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="modificationMass"></param>
        /// <param name="chAffectedAtom"></param>
        /// <returns></returns>
        /// <remarks>Returns True if index is valid; otherwise, returns false</remarks>
        public bool GetSearchResultModDetailsByIndex(
            int index,
            out char residue,
            out int residueLocInPeptide,
            out double modificationMass,
            out char chAffectedAtom)
        {
            if (index >= 0 && index < mSearchResultModifications.Count)
            {
                var resultMod = mSearchResultModifications[index];
                residue = resultMod.Residue;
                residueLocInPeptide = resultMod.ResidueLocInPeptide;
                modificationMass = resultMod.ModDefinition.ModificationMass;
                chAffectedAtom = resultMod.ModDefinition.AffectedAtom;
                return true;
            }

            residue = Convert.ToChar(0);
            residueLocInPeptide = 0;
            modificationMass = 0;
            chAffectedAtom = Convert.ToChar(0);
            return false;
        }

        /// <summary>
        /// Get the modification info, by index
        /// </summary>
        /// <param name="index">Modification entry index (0-based)</param>
        /// <returns>Modification info details if a valid index, otherwise nothing</returns>
        /// <remarks></remarks>
        public clsAminoAcidModInfo GetSearchResultModDetailsByIndex(int index)
        {
            if (index >= 0 && index < mSearchResultModifications.Count)
            {
                return mSearchResultModifications[index];
            }

            return null;
        }

        private void InitializeLocalVariables()
        {
            mErrorMessage = string.Empty;
            Clear();
        }

        /// <summary>
        /// Associates the given modification symbol with the given residue
        /// Returns false if the modification symbol is unknown
        /// </summary>
        /// <param name="modificationSymbol"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SearchResultAddDynamicModification(
            char modificationSymbol,
            char chTargetResidue,
            int residueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool updateModOccurrenceCounts)
        {
            var success = false;

            // Find the modification that uses this modification symbol and applies to this target residue
            var modificationDefinition = mPeptideMods.LookupDynamicModificationDefinitionByTargetInfo(
                modificationSymbol, chTargetResidue, eResidueTerminusState, out var existingModFound);

            if (existingModFound)
            {
                if (residueLocInPeptide < 1)
                {
                    // Invalid position; ignore this modification
                    mErrorMessage = "Invalid value for residueLocInPeptide: " + residueLocInPeptide.ToString();
                }
                else
                {
                    success = SearchResultAddModification(
                        modificationDefinition,
                        chTargetResidue,
                        residueLocInPeptide,
                        eResidueTerminusState,
                        updateModOccurrenceCounts);
                }

                return success;
            }

            // Modification not found
            mErrorMessage = "Modification symbol not found: " + modificationSymbol + "; TerminusState = " + eResidueTerminusState.ToString();
            return false;

        }

        /// <summary>
        /// Associates the given modification mass with the given residue
        /// If the modification mass is unknown, will auto-add it to the list of known modifications
        /// </summary>
        /// <param name="modificationMass"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <returns>True if mod successfully added</returns>
        /// <remarks></remarks>
        public bool SearchResultAddModification(
            double modificationMass,
            char chTargetResidue,
            int residueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool updateModOccurrenceCounts)
        {
            return SearchResultAddModification(modificationMass, chTargetResidue, residueLocInPeptide, eResidueTerminusState, updateModOccurrenceCounts, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION);
        }

        /// <summary>
        /// Associates the given modification mass with the given residue
        /// If the modification mass is unknown, will auto-add it to the list of known modifications
        /// </summary>
        /// <param name="modificationMass"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modMassDigitsOfPrecision">Digits of precision to use when comparing modificationMass to the masses of known mods</param>
        /// <returns>True if mod successfully added</returns>
        /// <remarks></remarks>
        public bool SearchResultAddModification(
            double modificationMass,
            char chTargetResidue,
            int residueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool updateModOccurrenceCounts, byte modMassDigitsOfPrecision)
        {
            var success = false;

            if (residueLocInPeptide < 1)
            {
                // Invalid position; ignore this modification
                mErrorMessage = "Invalid value for residueLocInPeptide: " + residueLocInPeptide.ToString();
            }
            else
            {
                // Lookup the modification definition given the modification information
                // If the modification mass is unknown, will auto-add it to the list of known modifications
                var modificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(
                    modificationMass,
                    chTargetResidue,
                    eResidueTerminusState,
                    out _,
                    true, modMassDigitsOfPrecision);

                success = SearchResultAddModification(
                    modificationDefinition,
                    chTargetResidue,
                    residueLocInPeptide,
                    eResidueTerminusState,
                    updateModOccurrenceCounts);
            }

            return success;
        }

        /// <summary>
        /// Associates the given modification with the given residue
        /// </summary>
        /// <param name="modificationDefinition"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="residueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SearchResultAddModification(
            clsModificationDefinition modificationDefinition,
            char chTargetResidue,
            int residueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool updateModOccurrenceCounts)
        {
            var success = false;

            if (residueLocInPeptide < 1 && modificationDefinition.ModificationType != clsModificationDefinition.eModificationTypeConstants.IsotopicMod)
            {
                // Invalid position; ignore this modification
                mErrorMessage = "Invalid value for residueLocInPeptide: " + residueLocInPeptide + " (modificationDefinition.ModificationType = " + modificationDefinition.ModificationType + ")";
            }
            else
            {
                if (updateModOccurrenceCounts)
                {
                    // Increment OccurrenceCount
                    modificationDefinition.OccurrenceCount += 1;
                }

                mSearchResultModifications.Add(new clsAminoAcidModInfo(chTargetResidue, residueLocInPeptide, eResidueTerminusState, modificationDefinition));
                success = true;
            }

            return success;
        }

        /// <summary>
        /// Adds any defined isotopic modifications to the peptide
        /// </summary>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SearchResultAddIsotopicModifications(bool updateModOccurrenceCounts)
        {
            var success = false;

            for (var index = 0; index <= mPeptideMods.ModificationCount - 1; index++)
            {
                if (mPeptideMods.GetModificationTypeByIndex(index) == clsModificationDefinition.eModificationTypeConstants.IsotopicMod)
                {
                    var residueLocInPeptide = 0;

                    var mod = mPeptideMods.GetModificationByIndex(index);
                    // TODO: VERIFY: success = SearchResultAddModification(ref mPeptideMods.GetModificationByIndex(index), clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, residueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, updateModOccurrenceCounts);
                    success = SearchResultAddModification(mod, clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, residueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, updateModOccurrenceCounts);
                }
            }

            return success;
        }

        /// <summary>
        /// Add any protein or peptide terminus static mods that are defined
        /// </summary>
        /// <param name="allowDuplicateModOnTerminus">When false, only add the modification if the terminal residue does not already have the given modification associated with it</param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <remarks>
        /// Peptide terminus static mods are always considered for a given peptide
        /// Protein terminus static mods are only considered if the peptide is at the appropriate terminus for the modification
        /// </remarks>
        public void SearchResultAddStaticTerminusMods(bool allowDuplicateModOnTerminus, bool updateModOccurrenceCounts)
        {
            var residueLocInPeptide = 0;
            clsModificationDefinition modificationDefinition = null;

            for (var modificationIndex = 0; modificationIndex <= mPeptideMods.ModificationCount - 1; modificationIndex++)
            {
                var addModification = false;

                var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;

                if (mPeptideMods.GetModificationTypeByIndex(modificationIndex) == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod)
                {
                    modificationDefinition = mPeptideMods.GetModificationByIndex(modificationIndex);
                    if (modificationDefinition.TargetResidues == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        residueLocInPeptide = 1;
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus ||
                            mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                        }
                        else
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                        }
                        addModification = true;
                    }
                    else if (modificationDefinition.TargetResidues == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        residueLocInPeptide = mPeptideCleanSequence.Length;
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus ||
                            mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                        }
                        else
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                        }
                        addModification = true;
                    }
                    else
                    {
                        // Invalid target residue for a peptide terminus static mod; do not add the modification
                    }
                }
                else if (mPeptideMods.GetModificationTypeByIndex(modificationIndex) == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod)
                {
                    modificationDefinition = mPeptideMods.GetModificationByIndex(modificationIndex);
                    if (modificationDefinition.TargetResidues == clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString())
                    {
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus ||
                            mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            residueLocInPeptide = 1;
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                            addModification = true;
                        }
                    }
                    else if (modificationDefinition.TargetResidues == clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString())
                    {
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus ||
                            mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            residueLocInPeptide = mPeptideCleanSequence.Length;
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                            addModification = true;
                        }
                    }
                    else
                    {
                        // Invalid target residue for a protein terminus static mod; do not add the modification
                    }
                }

                if (!addModification)
                    continue;

                if (!allowDuplicateModOnTerminus)
                {
                    // Look through mSearchResultModifications to see if this residue already has a modification with this modification's mass
                    // If it does, do not re-add the modification
                    for (var indexCompare = 0; indexCompare <= mSearchResultModifications.Count - 1; indexCompare++)
                    {
                        if (ReferenceEquals(modificationDefinition, mSearchResultModifications[indexCompare].ModDefinition))
                        {
                            addModification = false;
                            break;
                        }

                        if (mSearchResultModifications[indexCompare].ResidueLocInPeptide == residueLocInPeptide)
                        {
                            // First compare the MassCorrectionTag names
                            if (mSearchResultModifications[indexCompare].ModDefinition.MassCorrectionTag == modificationDefinition.MassCorrectionTag)
                            {
                                addModification = false;
                                break;
                            }

                            var massDifference = Math.Abs(mSearchResultModifications[indexCompare].ModDefinition.ModificationMass - modificationDefinition.ModificationMass);
                            if (Math.Abs(Math.Round(massDifference, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)) < float.Epsilon)
                            {
                                addModification = false;
                                break;
                            }
                        }
                    }
                }

                if (addModification)
                {
                    SearchResultAddModification(
                        modificationDefinition,
                        mPeptideCleanSequence[residueLocInPeptide - 1],
                        residueLocInPeptide,
                        eResidueTerminusState,
                        updateModOccurrenceCounts);
                }
            }
        }

        /// <summary>
        /// Updates the N-Terminal mass applied to peptides when computing their mass if it is significantly different than the currently defined N-terminal peptide mass
        /// </summary>
        /// <param name="nTerminalMassChange"></param>
        /// <remarks></remarks>
        public void UpdatePeptideNTerminusMass(double nTerminalMassChange)
        {
            if (Math.Round(Math.Abs(nTerminalMassChange - mPeptideSeqMassCalculator.PeptideNTerminusMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) > 0)
            {
                mPeptideSeqMassCalculator.PeptideNTerminusMass = nTerminalMassChange;
            }
        }

        /// <summary>
        /// Updates the C-Terminal mass applied to peptides when computing their mass if significantly different than the currently defined C-terminal peptide mass
        /// </summary>
        /// <param name="cTerminalMassChange"></param>
        /// <remarks></remarks>
        public void UpdatePeptideCTerminusMass(double cTerminalMassChange)
        {
            if (Math.Round(Math.Abs(cTerminalMassChange - mPeptideSeqMassCalculator.PeptideCTerminusMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) > 0)
            {
                mPeptideSeqMassCalculator.PeptideCTerminusMass = cTerminalMassChange;
            }
        }

        public void UpdateSearchResultEnzymeAndTerminusInfo(clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType udtEnzymeMatchSpec, double peptideNTerminusMassChange, double peptideCTerminusMassChange)
        {
            SetEnzymeMatchSpec(udtEnzymeMatchSpec);

            // Update the N-Terminus and/or C-Terminus masses if those in the XML file are significantly different than the defaults
            if (Math.Abs(peptideNTerminusMassChange - 0) > float.Epsilon)
            {
                UpdatePeptideNTerminusMass(peptideNTerminusMassChange);
            }

            if (Math.Abs(peptideCTerminusMassChange - 0) > float.Epsilon)
            {
                UpdatePeptideCTerminusMass(peptideCTerminusMassChange);
            }
        }

        /// <summary>
        /// Obtain the peptide sequence
        /// </summary>
        /// <param name="returnSequenceWithMods">When true, include the mod symbols in the sequence</param>
        /// <returns></returns>
        /// <remarks>
        /// If you want to guarantee that mod symbols are included in the peptide sequence,
        /// You must call ApplyModificationInformation before using this function
        /// </remarks>
        public string SequenceWithPrefixAndSuffix(bool returnSequenceWithMods)
        {
            string work;
            char chPrefix = default;
            char chSuffix = default;

            try
            {
                chPrefix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;

                if (mPeptidePreResidues != null)
                {
                    work = mPeptidePreResidues.Trim();
                    if (work.Length > 0)
                    {
                        chPrefix = work[work.Length - 1];
                        if (chPrefix == clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_XTANDEM_NTerminus)
                        {
                            chPrefix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
                        }
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }

            try
            {
                chSuffix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
                if (mPeptidePostResidues != null)
                {
                    work = mPeptidePostResidues.Trim();
                    if (work.Length > 0)
                    {
                        chSuffix = work[0];
                        if (chSuffix == clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_XTANDEM_CTerminus)
                        {
                            chSuffix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;
                        }
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }

            if (returnSequenceWithMods && PeptideSequenceWithMods != null)
            {
                return chPrefix + "." + PeptideSequenceWithMods + "." + chSuffix;
            }

            if (mPeptideCleanSequence == null)
            {
                return string.Empty;
            }

            return chPrefix + "." + mPeptideCleanSequence + "." + chSuffix;
        }

        /// <summary>
        /// Define custom RegEx specs used to find enzyme cleavage sites
        /// </summary>
        /// <param name="udtEnzymeMatchSpec"></param>
        /// <remarks>Define standard RegEx values using SetStandardEnzymeMatchSpec</remarks>
        public void SetEnzymeMatchSpec(clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType udtEnzymeMatchSpec)
        {
            mPeptideCleavageStateCalculator.SetEnzymeMatchSpec(udtEnzymeMatchSpec.LeftResidueRegEx, udtEnzymeMatchSpec.RightResidueRegEx);
        }

        /// <summary>
        /// Stores sequenceWithMods in PeptideSequenceWithMods
        /// </summary>
        /// <param name="sequenceWithMods"></param>
        /// <param name="checkForPrefixAndSuffixResidues"></param>
        /// <param name="autoPopulateCleanSequence">When true, populates PeptideCleanSequence, which automatically calls ComputePeptideCleavageStateInProtein</param>
        /// <remarks></remarks>
        public void SetPeptideSequenceWithMods(string sequenceWithMods, bool checkForPrefixAndSuffixResidues, bool autoPopulateCleanSequence)
        {
            // Sequence with mods, but without the prefix or suffix residues
            string primarySequence;

            if (checkForPrefixAndSuffixResidues)
            {
                if (!clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequenceWithMods, out primarySequence, out var prefix, out var suffix))
                {
                    primarySequence = string.Copy(sequenceWithMods);
                }
                else
                {
                    if (autoPopulateCleanSequence)
                    {
                        PeptidePreResidues = prefix;
                        PeptidePostResidues = suffix;
                    }
                }
            }
            else
            {
                primarySequence = string.Copy(sequenceWithMods);
            }

            if (autoPopulateCleanSequence)
            {
                // Note: Property PeptideCleanSequence will call ComputePeptideCleavageStateInProtein()
                PeptideCleanSequence = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(primarySequence, false);
            }

            PeptideSequenceWithMods = string.Copy(primarySequence);
        }

        /// <summary>
        /// Define standard RegEx values for finding enzyme cleavage sites
        /// </summary>
        /// <param name="eStandardCleavageAgent"></param>
        /// <remarks>Define custom RegEx values using SetEnzymeMatchSpec</remarks>
        public void SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants eStandardCleavageAgent)
        {
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(eStandardCleavageAgent);
        }

        /// <summary>
        /// Generate a comma separated list of the modifications present and the residue modified, with a colon separating the mod name and residue location
        /// Accessible via PeptideModDescription
        /// </summary>
        /// <remarks>
        /// Examples:
        ///   Acetyl:1
        ///   MinusH2O:1,Plus1Oxy:4,Plus1Oxy:19
        /// </remarks>
        public void UpdateModDescription()
        {
            const char MOD_LIST_SEP_CHAR = ',';

            PeptideModDescription = string.Empty;

            if (mSearchResultModifications.Count > 0)
            {
                var udtModNameAndResidueLoc = new clsPHRPBaseClass.udtModNameAndResidueLocType[mSearchResultModifications.Count];
                var pointerArray = new int[mSearchResultModifications.Count];

                if (mSearchResultModifications.Count == 1)
                {
                    pointerArray[0] = 0;
                }
                else
                {
                    // Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                    for (var index = 0; index <= mSearchResultModifications.Count - 1; index++)
                    {
                        udtModNameAndResidueLoc[index].ResidueLocInPeptide = mSearchResultModifications[index].ResidueLocInPeptide;
                        udtModNameAndResidueLoc[index].ModName = mSearchResultModifications[index].ModDefinition.MassCorrectionTag;
                        pointerArray[index] = index;
                    }

                    Array.Sort(udtModNameAndResidueLoc, pointerArray, new clsPHRPBaseClass.IModNameAndResidueLocComparer());
                }

                // Step through the modifications and add the modification name and residue position to mPeptideModDescription
                // Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
                for (var index = 0; index <= mSearchResultModifications.Count - 1; index++)
                {
                    var resultMods = mSearchResultModifications[pointerArray[index]];
                    if (index > 0)
                        PeptideModDescription += MOD_LIST_SEP_CHAR;
                    PeptideModDescription += resultMods.ModDefinition.MassCorrectionTag.Trim() + ':' + resultMods.ResidueLocInPeptide;
                }
            }
        }

        private class IGenericResidueModificationInfoComparer
            : IComparer<clsAminoAcidModInfo>
        {
            /// <summary>
            /// Comparer
            /// </summary>
            /// <param name="x"></param>
            /// <param name="y"></param>
            /// <returns></returns>
            /// <remarks></remarks>
            public int Compare(clsAminoAcidModInfo x, clsAminoAcidModInfo y)
            {
                if (x == null)
                    throw new ArgumentNullException(nameof(x));

                if (y == null)
                    throw new ArgumentNullException(nameof(y));

                if (x.ResidueLocInPeptide > y.ResidueLocInPeptide)
                {
                    return 1;
                }

                if (x.ResidueLocInPeptide < y.ResidueLocInPeptide)
                {
                    return -1;
                }

                return string.Compare(x.ModDefinition.MassCorrectionTag, y.ModDefinition.MassCorrectionTag, StringComparison.Ordinal);
            }
        }
    }
}
