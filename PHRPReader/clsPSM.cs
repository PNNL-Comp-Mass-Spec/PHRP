//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/04/2012
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;

namespace PHRPReader
{
    /// <summary>
    /// This class tracks the details for a peptide hit search result
    /// (typically loaded from a tab-delimited text file created by the Peptide File Extractor or by PHRP)
    /// </summary>
    public class clsPSM
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
        private readonly Dictionary<string, clsProteinInfo> mProteinDetails;

        /// <summary>
        /// Dictionary tracking additional, tool-specific scores
        /// </summary>
        private readonly Dictionary<string, string> mAdditionalScores;

        #region "Properties"

        /// <summary>
        /// Returns a dictionary with additional search engine scores stored as key/value pairs
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Update scores using SetScore</remarks>
        // ReSharper disable once UnusedMember.Global
        public IReadOnlyDictionary<string, string> AdditionalScores => mAdditionalScores;

        /// <summary>
        /// Assumed charge of the spectrum in which this peptide was identified
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public short Charge { get; set; }

        /// <summary>
        /// Peptide cleavage state (with regards to ProteinFirst)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>
        /// CleavageState, NumMissedCleavages, and NumTrypticTerminii are typically populated using UpdateCleavageInfo
        /// </remarks>
        public clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants CleavageState { get; set; }

        /// <summary>
        /// Collision mode (CID, ETD, HCD)
        /// PepXML allows this to be CID, ETD, ECD, ETD/CID, or HCD
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
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
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public float ElutionTimeMinutes { get; set; }

        /// <summary>
        /// Mass difference, in daltons, between the monoisotopic mass of the precursor ion and the calculated (theoretical) monoisotopic mass of the peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string MassErrorDa { get; set; }

        /// <summary>
        /// Mass difference, in ppm, between the monoisotopic mass of the precursor ion and the calculated (theoretical) monoisotopic mass of the peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string MassErrorPPM { get; set; }

        /// <summary>
        /// List of modified residues
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>A given residue is allowed to have more than one modification</remarks>
        public List<clsAminoAcidModInfo> ModifiedResidues { get; }

        /// <summary>
        /// MSGF Spectral E-Value associated with this peptide (aka SpecEValue or SpecProb)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>
        /// Ranges from 0 to 1, where 0 is the best score and 1 is the worse score
        /// Stored as a string to preserve formatting
        /// </remarks>
        public string MSGFSpecEValue { get; set; }

        /// <summary>
        /// MSGF Spectral E-Value associated with this peptide (aka SpecEValue or SpecProb)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
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
        /// <value></value>
        /// <returns></returns>
        /// <remarks>
        /// CleavageState, NumMissedCleavages, and NumTrypticTerminii are typically populated using UpdateCleavageInfo
        /// </remarks>
        public short NumMissedCleavages { get; set; }

        /// <summary>
        /// Number of tryptic terminii (or similar if not using trypsin)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>
        /// 2 means fully tryptic, 1 means partially tryptic, 0 means non-tryptic
        /// CleavageState, NumMissedCleavages, and NumTrypticTerminii are typically populated using UpdateCleavageInfo
        /// </remarks>
        public short NumTrypticTerminii { get; set; }

        /// <summary>
        /// Peptide sequence, including any modification symbols that were assigned by the search engine
        /// Peptide sequence, with or without prefix and suffix residues; may contain mod symbols
        /// Example, R.AAS*PQDLAGGYTSSLACHR.A
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string Peptide { get; private set; }

        /// <summary>
        /// Peptide residues without any modification symbols or flanking residues
        /// For example, AASPQDLAGGYTSSLACHR
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string PeptideCleanSequence { get; private set; }

        /// <summary>
        /// Computed monoisotopic mass (uncharged, theoretical mass, including mods)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>This mass is computed by PHRP using the PrecursorNeutralMass plus any modification masses associated with the peptide's residues</remarks>
        public double PeptideMonoisotopicMass { get; set; }

        /// <summary>
        /// Peptide sequence where all modified residues have the modification masses displayed as numeric values
        /// For example, R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
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
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Retrieve full list of proteins using the Proteins property</remarks>
        // ReSharper disable once UnusedMember.Global
        public string ProteinFirst
        {
            get
            {
                if (mProteins.Count == 0)
                {
                    return string.Empty;
                }

                return mProteins[0];
            }
        }

        /// <summary>
        /// Uncharged monoisotopic mass of the precursor (observed mass based on m/z and charge)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks>This mass is based on the mass or m/z value reported by the search engine</remarks>
        public double PrecursorNeutralMass { get; set; }

        /// <summary>
        /// List of proteins associated with this peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public IReadOnlyList<string> Proteins => mProteins;

        /// <summary>
        /// Dictionary with info on each protein, including name, description, cleavage state, terminus state, residue start, and residue end
        /// </summary>
        public IReadOnlyDictionary<string, clsProteinInfo> ProteinDetails => mProteinDetails;

        /// <summary>
        /// ResultID of this peptide (typically assigned by the search engine)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int ResultID { get; set; }

        /// <summary>
        /// List of scans that were combined prior to identifying this peptide
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public SortedSet<int> ScanList { get; }

        /// <summary>
        /// Scan number of the mass spectrum in which this peptide was identified
        /// Will automatically update ScanList if it does not yet contain this scan number
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
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
        /// <value></value>
        /// <returns></returns>
        /// <remarks>Top scoring peptide is rank 1, next lowest score is rank 2, etc.</remarks>
        public int ScoreRank { get; set; }

        /// <summary>
        /// Sequence ID value assigned by PHRP
        /// Required for looking up information from the SeqInfo files
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int SeqID { get; set; }

        #endregion

        /// <summary>
        /// Constructor; auto-calls Clear()
        /// </summary>
        /// <remarks></remarks>
        public clsPSM()
        {
            ScanList = new SortedSet<int>();
            mProteins = new List<string>();
            mProteinDetails = new Dictionary<string, clsProteinInfo>(StringComparer.OrdinalIgnoreCase);
            ModifiedResidues = new List<clsAminoAcidModInfo>();
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
        /// <remarks></remarks>
        public void AddModifiedResidue(clsAminoAcidModInfo modInfo)
        {
            ModifiedResidues.Add(modInfo);
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="residue">Amino acid letter; use angle brackets or square brackets for peptide or protein terminii (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
        /// <param name="residueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
        /// <param name="residueTerminusState">Terminus state of residue</param>
        /// <param name="modDefinition">Modification details</param>
        /// <remarks></remarks>
        public void AddModifiedResidue(
            char residue,
            int residueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants residueTerminusState,
            clsModificationDefinition modDefinition)
        {
            ModifiedResidues.Add(new clsAminoAcidModInfo(residue, residueLocInPeptide, residueTerminusState, modDefinition));
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="residue">Amino acid letter; use angle brackets or square brackets for peptide or protein terminii (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
        /// <param name="residueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
        /// <param name="residueTerminusState">Terminus state of residue</param>
        /// <param name="modDefinition">Modification details</param>
        /// <param name="endResidueLocInPeptide">For ambiguous mods, the residue number of the last residue that could have this modification</param>
        /// <remarks></remarks>
        public void AddModifiedResidue(
            char residue,
            int residueLocInPeptide,
            clsAminoAcidModInfo.eResidueTerminusStateConstants residueTerminusState,
            clsModificationDefinition modDefinition,
            int endResidueLocInPeptide)
        {
            ModifiedResidues.Add(new clsAminoAcidModInfo(residue, residueLocInPeptide, residueTerminusState, modDefinition, endResidueLocInPeptide));
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
        public void AddProtein(clsProteinInfo proteinInfo)
        {
            AddProteinDetail(proteinInfo);
        }

        /// <summary>
        /// Add detailed info of a protein associated with this peptide
        /// </summary>
        /// <param name="proteinInfo"></param>
        /// <remarks>Updates both the Protein list and the ProteinDetails dictionary</remarks>
        public void AddProteinDetail(clsProteinInfo proteinInfo)
        {

            var proteinName = proteinInfo.ProteinName;

            if (mProteinDetails.ContainsKey(proteinName))
            {
                mProteinDetails[proteinName] = proteinInfo;
            }
            else
            {
                mProteinDetails.Add(proteinName, proteinInfo);
            }

            if (!mProteins.Contains(proteinName))
            {
                mProteins.Add(proteinName);
            }

        }

        /// <summary>
        /// Reset the peptide to default values (and empty strings)
        /// </summary>
        /// <remarks></remarks>
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

            CleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Unknown;
            NumMissedCleavages = 0;
            NumTrypticTerminii = 0;

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
        /// <returns></returns>
        /// <remarks></remarks>
        // ReSharper disable once UnusedMember.Global
        public clsPSM Clone()
        {
            var newPSM = new clsPSM
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
            newPSM.NumTrypticTerminii = NumTrypticTerminii;

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
                PeptideCleanSequence = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(Peptide, true);
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
        ///  Returns the value stored for the specified score (as a double)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <returns>Score if defined, otherwise 0</returns>
        /// <remarks></remarks>
        // ReSharper disable once UnusedMember.Global
        public double GetScoreDbl(string scoreName)
        {
            return GetScoreDbl(scoreName, 0);
        }

        /// <summary>
        ///  Returns the value stored for the specified score (as a double)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <param name="valueIfMissing">Value to return if the score is not defined</param>
        /// <returns>Score if defined, otherwise valueIfMissing</returns>
        /// <remarks></remarks>
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
        ///  Returns the value stored for the specified score (as an integer)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <returns>Score if defined, otherwise 0</returns>
        /// <remarks></remarks>
        // ReSharper disable once UnusedMember.Global
        public int GetScoreInt(string scoreName)
        {
            return GetScoreInt(scoreName, 0);
        }

        /// <summary>
        ///  Returns the value stored for the specified score (as an integer)
        /// </summary>
        /// <param name="scoreName">Score name</param>
        /// <param name="valueIfMissing">Value to return if the score is not defined</param>
        /// <returns>Score if defined, otherwise valueIfMissing</returns>
        /// <remarks></remarks>
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
        /// <remarks>Does not update the cleavage state info.  If updateCleanSequence is false, then call UpdateCleanSequence at a later time to populate mPeptideCleanSequence</remarks>
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
        /// Update the peptide sequence (auto-determines the clean sequence); also auto-update the the cleavage state info
        /// </summary>
        /// <param name="peptide">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
        /// <param name="cleavageStateCalculator">Cleavage state calculator object</param>
        /// <remarks></remarks>
        public void SetPeptide(string peptide, clsPeptideCleavageStateCalculator cleavageStateCalculator)
        {
            SetPeptide(peptide);
            UpdateCleavageInfo(cleavageStateCalculator);
        }

        /// <summary>
        /// Add/update an additional score to associate with this peptide
        /// </summary>
        /// <param name="scoreName"></param>
        /// <param name="scoreValue"></param>
        /// <remarks></remarks>
        public void SetScore(string scoreName, string scoreValue)
        {
            if (mAdditionalScores.ContainsKey(scoreName))
            {
                mAdditionalScores[scoreName] = scoreValue;
            }
            else
            {
                mAdditionalScores.Add(scoreName, scoreValue);
            }
        }

        /// <summary>
        /// Returns the value stored for the specified score
        /// </summary>
        /// <param name="scoreName"></param>
        /// <param name="scoreValue"></param>
        /// <returns>True if the score is defined, otherwise false</returns>
        /// <remarks></remarks>
        public bool TryGetScore(string scoreName, out string scoreValue)
        {
            scoreValue = string.Empty;
            if (mAdditionalScores.TryGetValue(scoreName, out scoreValue))
            {
                return true;
            }

            return false;
        }

        /// <summary>
        /// Auto-determine the number of missed cleavages, cleavage state, and number of tryptic terminii based on the peptide sequence
        /// </summary>
        /// <param name="cleavageStateCalculator"></param>
        /// <remarks></remarks>
        public void UpdateCleavageInfo(clsPeptideCleavageStateCalculator cleavageStateCalculator)
        {
            NumMissedCleavages = cleavageStateCalculator.ComputeNumberOfMissedCleavages(Peptide);

            CleavageState = cleavageStateCalculator.ComputeCleavageState(Peptide);

            if (CleavageState == clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full)
            {
                NumTrypticTerminii = 2;
            }
            else if (CleavageState == clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)
            {
                NumTrypticTerminii = 1;
            }
            else
            {
                NumTrypticTerminii = 0;
            }
        }
    }
}
