//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;

namespace PHRPReader.Data
{
    /// <summary>
    /// This class tracks the details for a peptide hit search result
    /// (typically loaded from a tab-delimited text file created by the Peptide File Extractor or by PHRP)
    /// </summary>
    public class PSM
    {
        /// <summary>
        /// Unknown collision mode
        /// </summary>
        public const string UNKNOWN_COLLISION_MODE = "n/a";

        private string mDataLineText = string.Empty;

        // Note: Be sure to update the Clone() function if you add new class-wide variables or properties

        private int mScanNumber;

        /// <summary>
        /// Peptide sequence where modified residues have the modification mass indicated as a number
        /// Example: R.N+144.102063SNPVIAELSQAINSGTLLSK+144.102063PS+79.9663PPLPPK+144.102063.R
        /// </summary>
        private string mPeptideWithNumericMods;

        /// <summary>
        /// Protein names
        /// </summary>
        /// <remarks>Note that names are case-sensitive</remarks>
        private readonly List<string> mProteins;

        /// <summary>
        /// Dictionary with info on each protein, including name, description, cleavage state, terminus state, residue start, and residue end
        /// </summary>
        private readonly Dictionary<string, ProteinInfo> mProteinDetails;

        /// <summary>
        /// Dictionary tracking additional, tool-specific scores
        /// </summary>
        private readonly Dictionary<string, string> mAdditionalScores;

        #region "Properties"

        /// <summary>
        /// Returns a dictionary with additional search engine scores stored as key/value pairs
        /// </summary>
        /// <remarks>Update scores using SetScore</remarks>
        // ReSharper disable once UnusedMember.Global
        public IReadOnlyDictionary<string, string> AdditionalScores => mAdditionalScores;

        /// <summary>
        /// Assumed charge of the spectrum in which this peptide was identified
        /// </summary>
        public short Charge { get; set; }

        /// <summary>
        /// Peptide cleavage state (with regards to ProteinFirst)
        /// </summary>
        /// <remarks>
        /// CleavageState, NumMissedCleavages, and NumTrypticTermini are typically populated using UpdateCleavageInfo
        /// </remarks>
        public PeptideCleavageStateCalculator.PeptideCleavageStateConstants CleavageState { get; set; }

        /// <summary>
        /// Collision mode (CID, ETD, HCD)
        /// PepXML allows this to be CID, ETD, ECD, ETD/CID, or HCD
        /// </summary>
        public string CollisionMode { get; set; }

        /// <summary>
        /// Single line of data read from a PHRP data file
        /// </summary>
        public string DataLineText
        {
            get => mDataLineText;
            set
            {
                if (string.IsNullOrEmpty(value))
                {
                    mDataLineText = string.Empty;
                }
                else
                {
                    mDataLineText = value;
                }
            }
        }

        /// <summary>
        /// Elution time (in minutes) of the spectrum
        /// </summary>
        public float ElutionTimeMinutes { get; set; }

        /// <summary>
        /// Mass difference, in daltons, between the monoisotopic mass of the precursor ion and the calculated (theoretical) monoisotopic mass of the peptide
        /// </summary>
        public string MassErrorDa { get; set; }

        /// <summary>
        /// Mass difference, in ppm, between the monoisotopic mass of the precursor ion and the calculated (theoretical) monoisotopic mass of the peptide
        /// </summary>
        public string MassErrorPPM { get; set; }

        /// <summary>
        /// List of modified residues
        /// </summary>
        /// <remarks>A given residue is allowed to have more than one modification</remarks>
        public List<AminoAcidModInfo> ModifiedResidues { get; }

        /// <summary>
        /// MSGF Spectral E-Value associated with this peptide (aka SpecEValue or SpecProb)
        /// </summary>
        /// <remarks>
        /// Ranges from 0 to 1, where 0 is the best score and 1 is the worse score
        /// Stored as a string to preserve formatting
        /// </remarks>
        public string MSGFSpecEValue { get; set; }

        /// <summary>
        /// MSGF Spectral E-Value associated with this peptide (aka SpecEValue or SpecProb)
        /// </summary>
        /// <remarks>
        /// Ranges from 0 to 1, where 0 is the best score and 1 is the worse score
        /// Stored as a string to preserve formatting
        /// </remarks>
        [Obsolete("Use MSGFSpecEValue")]
        public string MSGFSpecProb
        {
            get => MSGFSpecEValue;
            set => MSGFSpecEValue = value;
        }

        /// <summary>
        /// Number of missed cleavages (internal K or R)
        /// </summary>
        /// <remarks>
        /// CleavageState, NumMissedCleavages, and NumTrypticTermini are typically populated using UpdateCleavageInfo
        /// </remarks>
        public short NumMissedCleavages { get; set; }

        /// <summary>
        /// Number of tryptic termini (or similar if not using trypsin)
        /// </summary>
        /// <remarks>
        /// 2 means fully tryptic, 1 means partially tryptic, 0 means non-tryptic
        /// CleavageState, NumMissedCleavages, and NumTrypticTermini are typically populated using UpdateCleavageInfo
        /// </remarks>
        public short NumTrypticTermini { get; set; }

        /// <summary>
        /// Peptide sequence, including any modification symbols that were assigned by the search engine
        /// Peptide sequence, with or without prefix and suffix residues; may contain mod symbols
        /// Example, R.AAS*PQDLAGGYTSSLACHR.A
        /// </summary>
        public string Peptide { get; private set; }

        /// <summary>
        /// Peptide residues without any modification symbols or flanking residues
        /// For example, AASPQDLAGGYTSSLACHR
        /// </summary>
        public string PeptideCleanSequence { get; private set; }

        /// <summary>
        /// Computed monoisotopic mass (uncharged, theoretical mass, including mods)
        /// </summary>
        /// <remarks>This mass is computed by PHRP using the PrecursorNeutralMass plus any modification masses associated with the peptide's residues</remarks>
        public double PeptideMonoisotopicMass { get; set; }

        /// <summary>
        /// Peptide sequence where all modified residues have the modification masses displayed as numeric values
        /// For example, R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A
        /// </summary>
        public string PeptideWithNumericMods
        {
            get => mPeptideWithNumericMods;
            set
            {
                if (string.IsNullOrEmpty(value))
                {
                    mPeptideWithNumericMods = string.Empty;
                }
                else
                {
                    mPeptideWithNumericMods = value;
                }
            }
        }

        /// <summary>
        /// First protein associated with this peptide
        /// </summary>
        /// <remarks>Retrieve full list of proteins using the Proteins property</remarks>
        // ReSharper disable once UnusedMember.Global
        public string ProteinFirst
        {
            get
            {
                return mProteins.Count == 0 ? string.Empty : mProteins[0];
            }
        }

        /// <summary>
        /// Uncharged monoisotopic mass of the precursor (observed mass based on m/z and charge)
        /// </summary>
        /// <remarks>This mass is based on the mass or m/z value reported by the search engine</remarks>
        public double PrecursorNeutralMass { get; set; }

        /// <summary>
        /// List of proteins associated with this peptide
        /// </summary>
        public IReadOnlyList<string> Proteins => mProteins;

        /// <summary>
        /// Dictionary with info on each protein, including name, description, cleavage state, terminus state, residue start, and residue end
        /// </summary>
        public IReadOnlyDictionary<string, ProteinInfo> ProteinDetails => mProteinDetails;

        /// <summary>
        /// ResultID of this peptide (typically assigned by the search engine)
        /// </summary>
        public int ResultID { get; set; }

        /// <summary>
        /// List of scans that were combined prior to identifying this peptide
        /// </summary>
        public SortedSet<int> ScanList { get; }

        /// <summary>
        /// Scan number of the mass spectrum in which this peptide was identified
        /// Will automatically update ScanList if it does not yet contain this scan number
        /// </summary>
        public int ScanNumber
        {
            get => mScanNumber;
            set
            {
                mScanNumber = value;
                if (!ScanList.Contains(value))
                {
                    ScanList.Add(value);
                }
            }
        }

        /// <summary>
        /// First scan number
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public int ScanNumberStart => ScanList.Min;

        /// <summary>
        /// Last scan number
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public int ScanNumberEnd => ScanList.Max;

        /// <summary>
        /// Rank of this peptide in the given spectrum
        /// </summary>
        /// <remarks>Top scoring peptide is rank 1, next lowest score is rank 2, etc.</remarks>
        public int ScoreRank { get; set; }

        /// <summary>
        /// Sequence ID value assigned by PHRP
        /// Required for looking up information from the SeqInfo files
        /// </summary>
        public int SeqID { get; set; }

        #endregion

        /// <summary>
        /// Constructor; auto-calls Clear()
        /// </summary>
        public PSM()
        {
            ScanList = new SortedSet<int>();
            mProteins = new List<string>();
            mProteinDetails = new Dictionary<string, ProteinInfo>(StringComparer.OrdinalIgnoreCase);
            ModifiedResidues = new List<AminoAcidModInfo>();
            mAdditionalScores = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            Clear();
        }

        /// <summary>
        /// Add an additional scan number to associate with this PSM
        /// </summary>
        /// <param name="scanNumber"></param>
        public void AddCombinedScan(int scanNumber)
        {
            if (!ScanList.Contains(scanNumber))
            {
                ScanList.Add(scanNumber);
            }
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="modInfo">Modification info class</param>
        public void AddModifiedResidue(AminoAcidModInfo modInfo)
        {
            ModifiedResidues.Add(modInfo);
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="residue">Amino acid letter; use angle brackets or square brackets for peptide or protein termini (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
        /// <param name="residueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
        /// <param name="residueTerminusState">Terminus state of residue</param>
        /// <param name="modDefinition">Modification details</param>
        public void AddModifiedResidue(
            char residue,
            int residueLocInPeptide,
            AminoAcidModInfo.ResidueTerminusStateConstants residueTerminusState,
            ModificationDefinition modDefinition)
        {
            ModifiedResidues.Add(new AminoAcidModInfo(residue, residueLocInPeptide, residueTerminusState, modDefinition));
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="residue">Amino acid letter; use angle brackets or square brackets for peptide or protein termini (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
        /// <param name="residueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
        /// <param name="residueTerminusState">Terminus state of residue</param>
        /// <param name="modDefinition">Modification details</param>
        /// <param name="endResidueLocInPeptide">For ambiguous mods, the residue number of the last residue that could have this modification</param>
        public void AddModifiedResidue(
            char residue,
            int residueLocInPeptide,
            AminoAcidModInfo.ResidueTerminusStateConstants residueTerminusState,
            ModificationDefinition modDefinition,
            int endResidueLocInPeptide)
        {
            ModifiedResidues.Add(new AminoAcidModInfo(residue, residueLocInPeptide, residueTerminusState, modDefinition, endResidueLocInPeptide));
        }

        /// <summary>
        /// Add a new protein to associate with this peptide
        /// </summary>
        /// <param name="proteinName"></param>
        /// <remarks>Does not update the ProteinDetails dictionary</remarks>
        public void AddProtein(string proteinName)
        {
            if (!string.IsNullOrWhiteSpace(proteinName) && !mProteins.Contains(proteinName))
            {
                mProteins.Add(proteinName);
            }
        }

        /// <summary>
        /// Add detailed info of a protein associated with this peptide
        /// </summary>
        /// <param name="proteinInfo"></param>
        /// <remarks>Updates both the Protein list and the ProteinDetails dictionary</remarks>
        public void AddProtein(ProteinInfo proteinInfo)
        {
            AddProteinDetail(proteinInfo);
        }

        /// <summary>
        /// Add detailed info of a protein associated with this peptide
        /// </summary>
        /// <param name="proteinInfo"></param>
        /// <remarks>Updates both the Protein list and the ProteinDetails dictionary</remarks>
        public void AddProteinDetail(ProteinInfo proteinInfo)
        {
            var proteinName = proteinInfo.ProteinName;

            // Add/update the dictionary
            mProteinDetails[proteinName] = proteinInfo;

            if (!mProteins.Contains(proteinName))
            {
                mProteins.Add(proteinName);
            }
        }

        /// <summary>
        /// Reset the peptide to default values (and empty strings)
        /// </summary>
        public void Clear()
        {
            mDataLineText = string.Empty;
            mScanNumber = 0;
            ElutionTimeMinutes = 0;

            ScanList.Clear();

            Peptide = string.Empty;
            mPeptideWithNumericMods = string.Empty;
            PeptideCleanSequence = string.Empty;
            Charge = 0;
            ResultID = 0;
            ScoreRank = 0;

            CollisionMode = UNKNOWN_COLLISION_MODE;
            MSGFSpecEValue = string.Empty;

            CleavageState = PeptideCleavageStateCalculator.PeptideCleavageStateConstants.Unknown;
            NumMissedCleavages = 0;
            NumTrypticTermini = 0;

            PrecursorNeutralMass = 0;
            MassErrorDa = string.Empty;
            MassErrorPPM = string.Empty;

            PeptideMonoisotopicMass = 0;

            mProteins.Clear();
            mProteinDetails.Clear();

            ModifiedResidues.Clear();
            mAdditionalScores.Clear();
        }

        /// <summary>
        /// Clear any residue modifications
        /// </summary>
        public void ClearModifiedResidues()
        {
            ModifiedResidues.Clear();
        }

        /// <summary>
        /// Duplicate this PSM object and return a new one
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public PSM Clone()
        {
            var newPSM = new PSM
            {
                ResultID = ResultID,
                ScoreRank = ScoreRank,
                ScanNumber = mScanNumber,
                ElutionTimeMinutes = ElutionTimeMinutes
            };

            foreach (var scanNumber in ScanList)
            {
                newPSM.AddCombinedScan(scanNumber);
            }

            // Note: this call will auto-update mPeptideCleanSequence in newPSM
            newPSM.SetPeptide(Peptide);

            newPSM.PeptideWithNumericMods = mPeptideWithNumericMods;
            newPSM.Charge = Charge;
            newPSM.CollisionMode = CollisionMode;
            newPSM.MSGFSpecEValue = MSGFSpecEValue;

            newPSM.CleavageState = CleavageState;
            newPSM.NumMissedCleavages = NumMissedCleavages;
            newPSM.NumTrypticTermini = NumTrypticTermini;

            newPSM.PrecursorNeutralMass = PrecursorNeutralMass;
            newPSM.MassErrorDa = MassErrorDa;
            newPSM.MassErrorPPM = MassErrorPPM;

            foreach (var protein in mProteins)
            {
                newPSM.AddProtein(protein);
            }

            foreach (var item in mProteinDetails.Values)
            {
                newPSM.AddProteinDetail(item);
            }

            foreach (var item in ModifiedResidues)
            {
                newPSM.AddModifiedResidue(item.Residue, item.ResidueLocInPeptide, item.ResidueTerminusState, item.ModDefinition);
            }

            foreach (var score in mAdditionalScores)
            {
                newPSM.SetScore(score.Key, score.Value);
            }

            return newPSM;
        }

        /// <summary>
        /// Update the clean sequence
        /// </summary>
        public void UpdateCleanSequence()
        {
            if (string.IsNullOrEmpty(Peptide))
            {
                PeptideCleanSequence = string.Empty;
            }
            else
            {
                PeptideCleanSequence = PeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(Peptide, true);
            }
        }

        /// <summary>
        /// Returns the value stored for the specified score
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <returns>Score if defined, otherwise an empty string</returns>
        public string GetScore(string scoreName)
        {
            if (mAdditionalScores.TryGetValue(scoreName, out var scoreValue))
            {
                return scoreValue;
            }

            return string.Empty;
        }

        /// <summary>
        /// Returns the value stored for the specified score (as a double)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <returns>Score if defined, otherwise 0</returns>
        // ReSharper disable once UnusedMember.Global
        public double GetScoreDbl(string scoreName)
        {
            return GetScoreDbl(scoreName, 0);
        }

        /// <summary>
        /// Returns the value stored for the specified score (as a double)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <param name="valueIfMissing">Value to return if the score is not defined</param>
        /// <returns>Score if defined, otherwise valueIfMissing</returns>
        public double GetScoreDbl(string scoreName, double valueIfMissing)
        {
            var scoreValue = GetScore(scoreName);
            if (!string.IsNullOrEmpty(scoreValue) && double.TryParse(scoreValue, out var score))
            {
                return score;
            }

            return valueIfMissing;
        }

        /// <summary>
        /// Returns the value stored for the specified score (as an integer)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <returns>Score if defined, otherwise 0</returns>
        // ReSharper disable once UnusedMember.Global
        public int GetScoreInt(string scoreName)
        {
            return GetScoreInt(scoreName, 0);
        }

        /// <summary>
        /// Returns the value stored for the specified score (as an integer)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <param name="valueIfMissing">Value to return if the score is not defined</param>
        /// <returns>Score if defined, otherwise valueIfMissing</returns>
        public int GetScoreInt(string scoreName, int valueIfMissing)
        {
            var scoreValue = GetScore(scoreName);
            if (!string.IsNullOrEmpty(scoreValue) && int.TryParse(scoreValue, out var score))
            {
                return score;
            }

            return valueIfMissing;
        }

        /// <summary>
        /// Update the peptide sequence, auto-determining the clean sequence if updateCleanSequence is true
        /// </summary>
        /// <param name="peptideSequence">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
        /// <param name="updateCleanSequence"></param>
        /// <remarks>Does not update the cleavage state info.  If updateCleanSequence is false, call UpdateCleanSequence at a later time to populate mPeptideCleanSequence</remarks>
        public void SetPeptide(string peptideSequence, bool updateCleanSequence = true)
        {
            if (string.IsNullOrEmpty(peptideSequence))
            {
                Peptide = string.Empty;
            }
            else
            {
                Peptide = peptideSequence;
            }

            if (updateCleanSequence)
            {
                UpdateCleanSequence();
            }
        }

        /// <summary>
        /// Update the peptide sequence (auto-determines the clean sequence); also auto-update the cleavage state info
        /// </summary>
        /// <param name="peptide">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
        /// <param name="cleavageStateCalculator">Cleavage state calculator object</param>
        public void SetPeptide(string peptide, PeptideCleavageStateCalculator cleavageStateCalculator)
        {
            SetPeptide(peptide);
            UpdateCleavageInfo(cleavageStateCalculator);
        }

        /// <summary>
        /// Add/update an additional score to associate with this peptide
        /// </summary>
        /// <param name="scoreName"></param>
        /// <param name="scoreValue"></param>
        public void SetScore(string scoreName, string scoreValue)
        {
            // Add/update the dictionary
            mAdditionalScores[scoreName] = scoreValue;
        }

        /// <summary>
        /// Returns the value stored for the specified score
        /// </summary>
        /// <param name="scoreName"></param>
        /// <param name="scoreValue"></param>
        /// <returns>True if the score is defined, otherwise false</returns>
        public bool TryGetScore(string scoreName, out string scoreValue)
        {
            return mAdditionalScores.TryGetValue(scoreName, out scoreValue);
        }

        /// <summary>
        /// Auto-determine the number of missed cleavages, cleavage state, and number of tryptic termini based on the peptide sequence
        /// </summary>
        /// <param name="cleavageStateCalculator"></param>
        public void UpdateCleavageInfo(PeptideCleavageStateCalculator cleavageStateCalculator)
        {
            NumMissedCleavages = cleavageStateCalculator.ComputeNumberOfMissedCleavages(Peptide);

            CleavageState = cleavageStateCalculator.ComputeCleavageState(Peptide);

            if (CleavageState == PeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)
            {
                NumTrypticTermini = 2;
            }
            else if (CleavageState == PeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)
            {
                NumTrypticTermini = 1;
            }
            else
            {
                NumTrypticTermini = 0;
            }
        }
    }
}
