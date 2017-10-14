// This class is used to track the peptide details for a given MS/MS search result
// It can track peptide residue level modifications and static, peptide-wide modifications
//
// Use SearchResultClearModifications() and SearchResultAddModification() to track
//  specific modifications for a given peptide
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 7, 2006
// Last updated June 26, 2006
//
// E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
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
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public abstract class clsSearchResultsBaseClass
    {
        #region "Constants and Enums"
        /// <summary>
        /// Mass difference, in Da, between ^12C and ^13C
        /// </summary>
        /// <remarks></remarks>
        public const double MASS_C13 = 1.00335483;
        #endregion

        #region "Classwide Variables"
        // Note: Many of these variables typically hold numbers but we're storing the numbers as strings
        //       to prevent the numeric representation from changing when converting to a number then back to a string

        /// <summary>
        /// Residue or residues before the start of the peptide sequence
        /// </summary>
        /// <remarks></remarks>
        protected string mPeptidePreResidues;

        /// <summary>
        /// Residue or residues after the end of the peptide sequence
        /// </summary>
        /// <remarks></remarks>
        protected string mPeptidePostResidues;

        /// <summary>
        /// Peptide sequence without any modification symbols
        /// </summary>
        /// <remarks></remarks>
        protected string mPeptideCleanSequence;

        /// <summary>
        /// Cleavage state of the peptide
        /// </summary>
        /// <remarks></remarks>
        protected clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants mPeptideCleavageState;

        /// <summary>
        /// Terminus state of the peptide
        /// </summary>
        /// <remarks></remarks>
        protected clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants mPeptideTerminusState;

        /// <summary>
        /// List of modifications present in the current peptide
        /// </summary>
        /// <remarks></remarks>
        protected readonly List<clsAminoAcidModInfo> mSearchResultModifications;

        /// <summary>
        /// Possible modifications that the peptide could have
        /// </summary>
        /// <remarks></remarks>
        protected readonly clsPeptideModificationContainer mPeptideMods;

        protected readonly clsPeptideCleavageStateCalculator mPeptideCleavageStateCalculator;
        protected readonly clsPeptideMassCalculator mPeptideSeqMassCalculator;

        protected string mErrorMessage = "";
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
        /// In Sequest it was historically the average mass MH, though when a monoisotopic mass parent tolerance is specified, then this is a monoisotopic mass
        /// In Inspect, MSGF+, and MSAlign, this is the theoretical monoisotopic MH; note that this is (M+H)+
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string PeptideMH { get; set; }

        /// <summary>
        /// Difference in mass between the peptide's computed mass and the parent ion mass (i.e. the mass chosen for fragmentation)
        /// In Sequest this is Theoretical Mass - Observed Mass
        /// In XTandem, Inspect, MSGF+, and MSAlign the DelM value is listed as Observed - Theoretical,
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
        /// <param name="objPeptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        /// <remarks></remarks>
        protected clsSearchResultsBaseClass(clsPeptideModificationContainer objPeptideMods, clsPeptideMassCalculator peptideSeqMassCalculator)
        {
            mSearchResultModifications = new List<clsAminoAcidModInfo>();

            mPeptideMods = objPeptideMods;

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
            // Initialize strSequenceWithMods to the clean sequence; we'll insert the mod symbols below if mSearchResultModifications.Count > 0
            var strSequenceWithMods = string.Copy(mPeptideCleanSequence);

            if (mSearchResultModifications.Count > 0)
            {
                // Insert the modification symbols into strSequenceWithMods
                // First, sort mSearchResultModifications on .ResidueLocInPeptide and .MassCorrectionTag

                if (mSearchResultModifications.Count > 1)
                {
                    mSearchResultModifications.Sort(new IGenericResidueModificationInfoComparer());
                }

                // Now step backward through intResidueModificationPositions and add the symbols to strSequenceWithMods
                for (var intIndex = mSearchResultModifications.Count - 1; intIndex >= 0; intIndex += -1)
                {
                    var resultMod = mSearchResultModifications[intIndex];
                    if (resultMod.ModDefinition.ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod ||
                        resultMod.ModDefinition.ModificationType == clsModificationDefinition.eModificationTypeConstants.UnknownType)
                    {
                        strSequenceWithMods = strSequenceWithMods.Insert(resultMod.ResidueLocInPeptide, resultMod.ModDefinition.ModificationSymbol.ToString());
                    }
                }
            }

            return strSequenceWithMods;
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
        /// <param name="dblDelM">Delta mass, in Da</param>
        /// <param name="dblPrecursorMonoMass">Precursor monoisotopic mass</param>
        /// <param name="blnAdjustPrecursorMassForC13">True to correct for C13 isotopic selection errors</param>
        /// <param name="dblPeptideMonoisotopicMass"></param>
        /// <returns>Delta mass, in ppm</returns>
        /// <remarks></remarks>
        public static double ComputeDelMCorrectedPPM(
            double dblDelM,
            double dblPrecursorMonoMass,
            bool blnAdjustPrecursorMassForC13,
            double dblPeptideMonoisotopicMass)
        {
            var intCorrectionCount = 0;

            // Examine dblDelM to determine which isotope was chosen
            if (dblDelM >= -0.5)
            {
                // This is the typical case
                while (dblDelM > 0.5)
                {
                    dblDelM -= MASS_C13;
                    intCorrectionCount += 1;
                }
            }
            else
            {
                // This happens less often; but we'll still account for it
                // In this case, intCorrectionCount will be negative
                while (dblDelM < -0.5)
                {
                    dblDelM += MASS_C13;
                    intCorrectionCount -= 1;
                }
            }

            if (intCorrectionCount != 0)
            {
                if (blnAdjustPrecursorMassForC13)
                {
                    // Adjust the precursor mono mass based on intCorrectionCount
                    dblPrecursorMonoMass -= intCorrectionCount * MASS_C13;
                }

                // Compute a new DelM value
                dblDelM = dblPrecursorMonoMass - dblPeptideMonoisotopicMass;
            }

            return clsPeptideMassCalculator.MassToPPM(dblDelM, dblPeptideMonoisotopicMass);
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
        /// <param name="intResidueLocInPeptide">Residue number (1-based)</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsAminoAcidModInfo.eResidueTerminusStateConstants DetermineResidueTerminusState(int intResidueLocInPeptide)
        {
            var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;
            if (intResidueLocInPeptide == 1)
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
                if (intResidueLocInPeptide == PeptideLocInProteinEnd - PeptideLocInProteinStart + 1)
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
        /// <param name="intIndex">Modification entry index (0-based)</param>
        /// <param name="chResidue"></param>
        /// <param name="intResidueLocInPeptide"></param>
        /// <param name="dblModificationMass"></param>
        /// <param name="chAffectedAtom"></param>
        /// <returns></returns>
        /// <remarks>Returns True if intIndex is valid; otherwise, returns false</remarks>
        public bool GetSearchResultModDetailsByIndex(
            int intIndex,
            out char chResidue,
            out int intResidueLocInPeptide,
            out double dblModificationMass,
            out char chAffectedAtom)
        {
            if (intIndex >= 0 & intIndex < mSearchResultModifications.Count)
            {
                var resultMod = mSearchResultModifications[intIndex];
                chResidue = resultMod.Residue;
                intResidueLocInPeptide = resultMod.ResidueLocInPeptide;
                dblModificationMass = resultMod.ModDefinition.ModificationMass;
                chAffectedAtom = resultMod.ModDefinition.AffectedAtom;
                return true;
            }

            chResidue = Convert.ToChar(0);
            intResidueLocInPeptide = 0;
            dblModificationMass = 0;
            chAffectedAtom = Convert.ToChar(0);
            return false;
        }

        /// <summary>
        /// Get the modification info, by index
        /// </summary>
        /// <param name="intIndex">Modification entry index (0-based)</param>
        /// <returns>Modification info details if a valid index, otherwise nothing</returns>
        /// <remarks></remarks>
        public clsAminoAcidModInfo GetSearchResultModDetailsByIndex(int intIndex)
        {
            if (intIndex >= 0 & intIndex < mSearchResultModifications.Count)
            {
                return mSearchResultModifications[intIndex];
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
        /// If the modification symbol is unknown, then will return False
        /// </summary>
        /// <param name="chModificationSymbol"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="intResidueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SearchResultAddDynamicModification(
            char chModificationSymbol,
            char chTargetResidue,
            int intResidueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool blnUpdateModOccurrenceCounts)
        {
            var blnSuccess = false;

            // Find the modification that uses this modification symbol and applies to this target residue
            var objModificationDefinition = mPeptideMods.LookupDynamicModificationDefinitionByTargetInfo(
                chModificationSymbol, chTargetResidue, eResidueTerminusState, out var blnExistingModFound);

            if (blnExistingModFound)
            {
                if (intResidueLocInPeptide < 1)
                {
                    // Invalid position; ignore this modification
                    mErrorMessage = "Invalid value for intResidueLocInPeptide: " + intResidueLocInPeptide.ToString();
                }
                else
                {
                    blnSuccess = SearchResultAddModification(
                        objModificationDefinition,
                        chTargetResidue,
                        intResidueLocInPeptide,
                        eResidueTerminusState,
                        blnUpdateModOccurrenceCounts);
                }

                return blnSuccess;
            }

            // Modification not found
            mErrorMessage = "Modification symbol not found: " + chModificationSymbol + "; TerminusState = " + eResidueTerminusState.ToString();
            return false;

        }

        /// <summary>
        /// Associates the given modification mass with the given residue
        /// If the modification mass is unknown, then will auto-add it to the list of known modifications
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="intResidueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <returns>True if mod successfully added</returns>
        /// <remarks></remarks>
        public bool SearchResultAddModification(
            double dblModificationMass,
            char chTargetResidue,
            int intResidueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool blnUpdateModOccurrenceCounts)
        {
            return SearchResultAddModification(dblModificationMass, chTargetResidue, intResidueLocInPeptide, eResidueTerminusState, blnUpdateModOccurrenceCounts, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION);
        }

        /// <summary>
        /// Associates the given modification mass with the given residue
        /// If the modification mass is unknown, then will auto-add it to the list of known modifications
        /// </summary>
        /// <param name="dblModificationMass"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="intResidueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <param name="modMassDigitsOfPrecision">Digits of precision to use when comparinig dblModificationMass to the masses of known mods</param>
        /// <returns>True if mod successfully added</returns>
        /// <remarks></remarks>
        public bool SearchResultAddModification(
            double dblModificationMass,
            char chTargetResidue,
            int intResidueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool blnUpdateModOccurrenceCounts, byte modMassDigitsOfPrecision)
        {
            var blnSuccess = false;

            if (intResidueLocInPeptide < 1)
            {
                // Invalid position; ignore this modification
                mErrorMessage = "Invalid value for intResidueLocInPeptide: " + intResidueLocInPeptide.ToString();
            }
            else
            {
                // Lookup the modification definition given the modification information
                // If the modification mass is unknown, then will auto-add it to the list of known modifications
                var objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(
                    dblModificationMass,
                    chTargetResidue,
                    eResidueTerminusState,
                    out _,
                    true, modMassDigitsOfPrecision);

                blnSuccess = SearchResultAddModification(
                    objModificationDefinition,
                    chTargetResidue,
                    intResidueLocInPeptide,
                    eResidueTerminusState,
                    blnUpdateModOccurrenceCounts);
            }

            return blnSuccess;
        }

        /// <summary>
        /// Associates the given modification with the given residue
        /// </summary>
        /// <param name="objModificationDefinition"></param>
        /// <param name="chTargetResidue"></param>
        /// <param name="intResidueLocInPeptide"></param>
        /// <param name="eResidueTerminusState"></param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SearchResultAddModification(
            clsModificationDefinition objModificationDefinition,
            char chTargetResidue,
            int intResidueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants eResidueTerminusState,
            bool blnUpdateModOccurrenceCounts)
        {
            var blnSuccess = false;

            if (intResidueLocInPeptide < 1 & objModificationDefinition.ModificationType != clsModificationDefinition.eModificationTypeConstants.IsotopicMod)
            {
                // Invalid position; ignore this modification
                mErrorMessage = "Invalid value for intResidueLocInPeptide: " + intResidueLocInPeptide + " (objModificationDefinition.ModificationType = " + objModificationDefinition.ModificationType + ")";
            }
            else
            {
                if (blnUpdateModOccurrenceCounts)
                {
                    // Increment OccurrenceCount
                    objModificationDefinition.OccurrenceCount += 1;
                }

                mSearchResultModifications.Add(new clsAminoAcidModInfo(chTargetResidue, intResidueLocInPeptide, eResidueTerminusState, objModificationDefinition));
                blnSuccess = true;
            }

            return blnSuccess;
        }

        /// <summary>
        /// Adds any defined isotopic modifications to the peptide
        /// </summary>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool SearchResultAddIsotopicModifications(bool blnUpdateModOccurrenceCounts)
        {
            var blnSuccess = false;

            for (var intIndex = 0; intIndex <= mPeptideMods.ModificationCount - 1; intIndex++)
            {
                if (mPeptideMods.GetModificationTypeByIndex(intIndex) == clsModificationDefinition.eModificationTypeConstants.IsotopicMod)
                {
                    var intResidueLocInPeptide = 0;

                    var mod = mPeptideMods.GetModificationByIndex(intIndex);
                    // TODO: VERIFY: blnSuccess = SearchResultAddModification(ref mPeptideMods.GetModificationByIndex(intIndex), clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnUpdateModOccurrenceCounts);
                    blnSuccess = SearchResultAddModification(mod, clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, intResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants.None, blnUpdateModOccurrenceCounts);
                }
            }

            return blnSuccess;
        }

        /// <summary>
        /// Add any protein or peptide terminus static mods that are defined
        /// </summary>
        /// <param name="blnAllowDuplicateModOnTerminus">When false, only add the modification if the terminal residue does not already have the given modification associated with it</param>
        /// <param name="blnUpdateModOccurrenceCounts"></param>
        /// <remarks>
        /// Peptide terminus static mods are always considered for a given peptide
        /// Protein terminus static mods are only considered if the peptide is at the appropriate terminus for the modification
        /// </remarks>
        public void SearchResultAddStaticTerminusMods(bool blnAllowDuplicateModOnTerminus, bool blnUpdateModOccurrenceCounts)
        {
            var intResidueLocInPeptide = 0;
            clsModificationDefinition objModificationDefinition = null;

            for (var intModificationIndex = 0; intModificationIndex <= mPeptideMods.ModificationCount - 1; intModificationIndex++)
            {
                var blnAddModification = false;

                var eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None;

                if (mPeptideMods.GetModificationTypeByIndex(intModificationIndex) == clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod)
                {
                    objModificationDefinition = mPeptideMods.GetModificationByIndex(intModificationIndex);
                    if (objModificationDefinition.TargetResidues == clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        intResidueLocInPeptide = 1;
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus | mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                        }
                        else
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus;
                        }
                        blnAddModification = true;
                    }
                    else if (objModificationDefinition.TargetResidues == clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString())
                    {
                        intResidueLocInPeptide = mPeptideCleanSequence.Length;
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus | mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                        }
                        else
                        {
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus;
                        }
                        blnAddModification = true;
                    }
                    else
                    {
                        // Invalid target residue for a peptide terminus static mod; do not add the modification
                    }
                }
                else if (mPeptideMods.GetModificationTypeByIndex(intModificationIndex) == clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod)
                {
                    objModificationDefinition = mPeptideMods.GetModificationByIndex(intModificationIndex);
                    if (objModificationDefinition.TargetResidues == clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString())
                    {
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus | mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            intResidueLocInPeptide = 1;
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus;
                            blnAddModification = true;
                        }
                    }
                    else if (objModificationDefinition.TargetResidues == clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString())
                    {
                        if (mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus | mPeptideTerminusState == clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)
                        {
                            intResidueLocInPeptide = mPeptideCleanSequence.Length;
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus;
                            blnAddModification = true;
                        }
                    }
                    else
                    {
                        // Invalid target residue for a protein terminus static mod; do not add the modification
                    }
                }

                if (blnAddModification & (objModificationDefinition != null))
                {
                    if (!blnAllowDuplicateModOnTerminus)
                    {
                        // Look through udtResidueModificationInfo to see if this residue already has a modification with this modification's mass
                        // If it does, do not re-add the modification
                        for (var intIndexCompare = 0; intIndexCompare <= mSearchResultModifications.Count - 1; intIndexCompare++)
                        {
                            if (object.ReferenceEquals(objModificationDefinition, mSearchResultModifications[intIndexCompare].ModDefinition))
                            {
                                blnAddModification = false;
                                break;
                            }

                            if (mSearchResultModifications[intIndexCompare].ResidueLocInPeptide == intResidueLocInPeptide)
                            {
                                // First compare the MassCorrectionTag names
                                if (mSearchResultModifications[intIndexCompare].ModDefinition.MassCorrectionTag == objModificationDefinition.MassCorrectionTag)
                                {
                                    blnAddModification = false;
                                    break;
                                }

                                var dblMassDifference = Math.Abs(mSearchResultModifications[intIndexCompare].ModDefinition.ModificationMass - objModificationDefinition.ModificationMass);
                                if (Math.Abs(Math.Round(dblMassDifference, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)) < float.Epsilon)
                                {
                                    blnAddModification = false;
                                    break;
                                }
                            }
                        }
                    }

                    if (blnAddModification)
                    {
                        SearchResultAddModification(
                            objModificationDefinition,
                            mPeptideCleanSequence[intResidueLocInPeptide - 1],
                            intResidueLocInPeptide,
                            eResidueTerminusState,
                            blnUpdateModOccurrenceCounts);
                    }
                }
            }
        }

        /// <summary>
        /// Updates the N-Terminal mass applied to peptides when computing their mass if it is significantly different than the currently defined N-terminal peptide mass
        /// </summary>
        /// <param name="dblNTerminalMassChange"></param>
        /// <remarks></remarks>
        public void UpdatePeptideNTerminusMass(double dblNTerminalMassChange)
        {
            if (Math.Round(Math.Abs(dblNTerminalMassChange - mPeptideSeqMassCalculator.PeptideNTerminusMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) > 0)
            {
                mPeptideSeqMassCalculator.PeptideNTerminusMass = dblNTerminalMassChange;
            }
        }

        /// <summary>
        /// Updates the C-Terminal mass applied to peptides when computing their mass if significantly different than the currently defined C-terminal peptide mass
        /// </summary>
        /// <param name="dblCTerminalMassChange"></param>
        /// <remarks></remarks>
        public void UpdatePeptideCTerminusMass(double dblCTerminalMassChange)
        {
            if (Math.Round(Math.Abs(dblCTerminalMassChange - mPeptideSeqMassCalculator.PeptideCTerminusMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) > 0)
            {
                mPeptideSeqMassCalculator.PeptideCTerminusMass = dblCTerminalMassChange;
            }
        }

        public void UpdateSearchResultEnzymeAndTerminusInfo(clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType udtEnzymeMatchSpec, double dblPeptideNTerminusMassChange, double dblPeptideCTerminusMassChange)
        {
            SetEnzymeMatchSpec(udtEnzymeMatchSpec);

            // Update the N-Terminus and/or C-Terminus masses if those in the XML file are significantly different than the defaults
            if (Math.Abs(dblPeptideNTerminusMassChange - 0) > float.Epsilon)
            {
                UpdatePeptideNTerminusMass(dblPeptideNTerminusMassChange);
            }

            if (Math.Abs(dblPeptideCTerminusMassChange - 0) > float.Epsilon)
            {
                UpdatePeptideCTerminusMass(dblPeptideCTerminusMassChange);
            }
        }

        /// <summary>
        /// Obtain the peptide sequence
        /// </summary>
        /// <param name="blnReturnSequenceWithMods">When true, then include the mod symbols in the sequence</param>
        /// <returns></returns>
        /// <remarks>
        /// If you want to guarantee that mod symbols are included in the peptide sequence,
        /// You must call ApplyModificationInformation before using this function
        /// </remarks>
        public string SequenceWithPrefixAndSuffix(bool blnReturnSequenceWithMods)
        {
            string strWork;
            var chPrefix = default(char);
            var chSuffix = default(char);

            try
            {
                chPrefix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST;

                if (mPeptidePreResidues != null)
                {
                    strWork = mPeptidePreResidues.Trim();
                    if (strWork.Length > 0)
                    {
                        chPrefix = strWork[strWork.Length - 1];
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
                    strWork = mPeptidePostResidues.Trim();
                    if (strWork.Length > 0)
                    {
                        chSuffix = strWork[0];
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

            if (blnReturnSequenceWithMods && PeptideSequenceWithMods != null)
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
        /// Stores strSequenceWithMods in PeptideSequenceWithMods
        /// </summary>
        /// <param name="strSequenceWithMods"></param>
        /// <param name="blnCheckForPrefixAndSuffixResidues"></param>
        /// <param name="blnAutoPopulateCleanSequence">When true, populates PeptideCleanSequence, which automatically calls ComputePeptideCleavageStateInProtein</param>
        /// <remarks></remarks>
        public void SetPeptideSequenceWithMods(string strSequenceWithMods, bool blnCheckForPrefixAndSuffixResidues, bool blnAutoPopulateCleanSequence)
        {
            // Sequence with mods, but without the prefix or suffix residues
            string strPrimarySequence;

            if (blnCheckForPrefixAndSuffixResidues)
            {
                if (!clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, out strPrimarySequence, out var strPrefix, out var strSuffix))
                {
                    strPrimarySequence = string.Copy(strSequenceWithMods);
                }
                else
                {
                    if (blnAutoPopulateCleanSequence)
                    {
                        PeptidePreResidues = strPrefix;
                        PeptidePostResidues = strSuffix;
                    }
                }
            }
            else
            {
                strPrimarySequence = string.Copy(strSequenceWithMods);
            }

            if (blnAutoPopulateCleanSequence)
            {
                // Note: Property PeptideCleanSequence will call ComputePeptideCleavageStateInProtein()
                PeptideCleanSequence = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strPrimarySequence, false);
            }

            PeptideSequenceWithMods = string.Copy(strPrimarySequence);
        }

        /// <summary>
        /// Define standard RegEx values for finding enzyming cleavage sites
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
                var intPointerArray = new int[mSearchResultModifications.Count];

                if (mSearchResultModifications.Count == 1)
                {
                    intPointerArray[0] = 0;
                }
                else
                {
                    // Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                    for (var intIndex = 0; intIndex <= mSearchResultModifications.Count - 1; intIndex++)
                    {
                        udtModNameAndResidueLoc[intIndex].ResidueLocInPeptide = mSearchResultModifications[intIndex].ResidueLocInPeptide;
                        udtModNameAndResidueLoc[intIndex].ModName = mSearchResultModifications[intIndex].ModDefinition.MassCorrectionTag;
                        intPointerArray[intIndex] = intIndex;
                    }

                    Array.Sort(udtModNameAndResidueLoc, intPointerArray, new clsPHRPBaseClass.IModNameAndResidueLocComparer());
                }

                // Step through the modifications and add the modification name and residue position to mPeptideModDescription
                // Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
                for (var intIndex = 0; intIndex <= mSearchResultModifications.Count - 1; intIndex++)
                {
                    var resultMods = mSearchResultModifications[intPointerArray[intIndex]];
                    if (intIndex > 0)
                        PeptideModDescription += MOD_LIST_SEP_CHAR;
                    PeptideModDescription += resultMods.ModDefinition.MassCorrectionTag.Trim() + ':' + resultMods.ResidueLocInPeptide;
                }
            }
        }

        protected class IGenericResidueModificationInfoComparer
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
