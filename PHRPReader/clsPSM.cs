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
        /// List of scans that were combined prior to identifying this peptide
        /// </summary>
        private readonly SortedSet<int> mScanList;

        /// <summary>
        /// Peptide sequence, with or without prefix and suffix residues; may contain mod symbols; example: R.RM*VNSGSGADSAVDLNSIPVAMIAR.V
        /// </summary>
        private string mPeptide;

        /// <summary>
        /// Peptide sequence where modified residues have the modification mass indicated as a number, example: R.N+144.102063SNPVIAELSQAINSGTLLSK+144.102063PS+79.9663PPLPPK+144.102063.R
        /// </summary>
        private string mPeptideWithNumericMods;

        /// <summary>
        /// Pepide sequence without any mod symbols
        /// </summary>
        private string mPeptideCleanSequence;

        /// <summary>
        /// Modified residues
        /// </summary>
        private readonly List<clsAminoAcidModInfo> mModifiedPeptideResidues;

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
        public List<clsAminoAcidModInfo> ModifiedResidues => mModifiedPeptideResidues;

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
        /// For example, R.AAS*PQDLAGGYTSSLACHR.A
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string Peptide => mPeptide;

        /// <summary>
        /// Peptide residues without any modification symbols or flanking residues
        /// For example, AASPQDLAGGYTSSLACHR
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string PeptideCleanSequence => mPeptideCleanSequence;

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
        public List<string> Proteins => mProteins;

        public Dictionary<string, clsProteinInfo> ProteinDetails => mProteinDetails;

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
        public SortedSet<int> ScanList => mScanList;

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
                if (!mScanList.Contains(value))
                {
                    mScanList.Add(value);
                }
            }
        }

        public int ScanNumberStart => mScanList.Min;

        public int ScanNumberEnd => mScanList.Max;

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
            mScanList = new SortedSet<int>();
            mProteins = new List<string>();
            mProteinDetails = new Dictionary<string, clsProteinInfo>(StringComparer.CurrentCultureIgnoreCase);
            mModifiedPeptideResidues = new List<clsAminoAcidModInfo>();
            mAdditionalScores = new Dictionary<string, string>(StringComparer.CurrentCultureIgnoreCase);
            this.Clear();
        }

        public void AddCombinedScan(int intScanNumber)
        {
            if (!mScanList.Contains(intScanNumber))
            {
                mScanList.Add(intScanNumber);
            }
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="objModInfo">Modification info class</param>
        /// <remarks></remarks>
        public void AddModifiedResidue(clsAminoAcidModInfo objModInfo)
        {
            mModifiedPeptideResidues.Add(objModInfo);
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="Residue">Amino acid letter; use angle brackets or square brackes for peptide or protein terminii (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
        /// <param name="ResidueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
        /// <param name="ResidueTerminusState">Terminus state of residue</param>
        /// <param name="ModDefinition">Modification details</param>
        /// <remarks></remarks>
        public void AddModifiedResidue(char Residue, int ResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants ResidueTerminusState, clsModificationDefinition ModDefinition)
        {
            mModifiedPeptideResidues.Add(new clsAminoAcidModInfo(Residue, ResidueLocInPeptide, ResidueTerminusState, ModDefinition));
        }

        /// <summary>
        /// Add the details for a modified residue
        /// </summary>
        /// <param name="Residue">Amino acid letter; use angle brackets or square brackes for peptide or protein terminii (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
        /// <param name="ResidueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
        /// <param name="ResidueTerminusState">Terminus state of residue</param>
        /// <param name="ModDefinition">Modification details</param>
        /// <param name="EndResidueLocInPeptide">For ambiguous mods, the residue number of the last residue that could have this modification</param>
        /// <remarks></remarks>
        public void AddModifiedResidue(char Residue, int ResidueLocInPeptide, clsAminoAcidModInfo.eResidueTerminusStateConstants ResidueTerminusState, clsModificationDefinition ModDefinition, int EndResidueLocInPeptide)
        {
            mModifiedPeptideResidues.Add(new clsAminoAcidModInfo(Residue, ResidueLocInPeptide, ResidueTerminusState, ModDefinition, EndResidueLocInPeptide));
        }

        /// <summary>
        /// Add a new protein to associate with this peptide
        /// </summary>
        /// <param name="strProteinName"></param>
        /// <remarks></remarks>
        public void AddProtein(string strProteinName)
        {
            if (!string.IsNullOrWhiteSpace(strProteinName) && !mProteins.Contains(strProteinName))
            {
                mProteins.Add(strProteinName);
            }
        }

        /// <summary>
        /// Add new detailed protein info for this peptide
        /// </summary>
        /// <param name="oProteinInfo"></param>
        /// <remarks></remarks>
        public void AddProteinDetail(clsProteinInfo oProteinInfo)
        {
            clsProteinInfo oCachedInfo = null;
            if (mProteinDetails.TryGetValue(oProteinInfo.ProteinName, out oCachedInfo))
            {
                mProteinDetails[oProteinInfo.ProteinName] = oProteinInfo;
            }
            else
            {
                mProteinDetails.Add(oProteinInfo.ProteinName, oProteinInfo);
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

            mScanList.Clear();

            mPeptide = string.Empty;
            mPeptideWithNumericMods = string.Empty;
            mPeptideCleanSequence = string.Empty;
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

            mModifiedPeptideResidues.Clear();
            mAdditionalScores.Clear();
        }

        public void ClearModifiedResidues()
        {
            mModifiedPeptideResidues.Clear();
        }

        /// <summary>
        /// Duplicate this PSM object and return a new one
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public clsPSM Clone()
        {
            clsPSM objNew = null;

            objNew = new clsPSM();

            objNew.ResultID = ResultID;
            objNew.ScoreRank = ScoreRank;
            objNew.ScanNumber = mScanNumber;
            objNew.ElutionTimeMinutes = ElutionTimeMinutes;

            foreach (var intScanNumber in mScanList)
            {
                objNew.AddCombinedScan(intScanNumber);
            }

            objNew.SetPeptide(mPeptide);               // Note: this will auto-update mPeptideCleanSequence in objNew
            objNew.PeptideWithNumericMods = mPeptideWithNumericMods;
            objNew.Charge = Charge;
            objNew.CollisionMode = CollisionMode;
            objNew.MSGFSpecEValue = MSGFSpecEValue;

            objNew.CleavageState = CleavageState;
            objNew.NumMissedCleavages = NumMissedCleavages;
            objNew.NumTrypticTerminii = NumTrypticTerminii;

            objNew.PrecursorNeutralMass = PrecursorNeutralMass;
            objNew.MassErrorDa = MassErrorDa;
            objNew.MassErrorPPM = MassErrorPPM;

            foreach (var strProtein in mProteins)
            {
                objNew.AddProtein(strProtein);
            }

            foreach (var item in mProteinDetails.Values)
            {
                objNew.AddProteinDetail(item);
            }

            foreach (var objItem in mModifiedPeptideResidues)
            {
                objNew.AddModifiedResidue(objItem.Residue, objItem.ResidueLocInPeptide, objItem.ResidueTerminusState, objItem.ModDefinition);
            }

            foreach (var objScore in mAdditionalScores)
            {
                objNew.SetScore(objScore.Key, objScore.Value);
            }

            return objNew;
        }

        public void UpdateCleanSequence()
        {
            UpdateCleanSequence(mPeptide);
        }

        private void UpdateCleanSequence(string strPeptide)
        {
            if (string.IsNullOrEmpty(strPeptide))
            {
                mPeptideCleanSequence = string.Empty;
            }
            else
            {
                mPeptideCleanSequence = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strPeptide, true);
            }
        }

        /// <summary>
        /// Returns the value stored for the specified score
        /// </summary>
        /// <param name="strScoreName">Score name</param>
        /// <returns>Score if defined, otherwise an empty string</returns>
        public string GetScore(string strScoreName)
        {
            var strScoreValue = string.Empty;
            if (mAdditionalScores.TryGetValue(strScoreName, out strScoreValue))
            {
                return strScoreValue;
            }
            else
            {
                return string.Empty;
            }
        }

        /// <summary>
        ///  Returns the value stored for the specified score (as a double)
        /// </summary>
        /// <param name="strScoreName">Score name</param>
        /// <returns>Score if defined, otherwise 0</returns>
        /// <remarks></remarks>
        public double GetScoreDbl(string strScoreName)
        {
            return GetScoreDbl(strScoreName, 0);
        }

        /// <summary>
        ///  Returns the value stored for the specified score (as a double)
        /// </summary>
        /// <param name="strScoreName">Score name</param>
        /// <param name="dblValueIfMissing">Value to return if the score is not defined</param>
        /// <returns>Score if defined, otherwise dblValueIfMissing</returns>
        /// <remarks></remarks>
        public double GetScoreDbl(string strScoreName, double dblValueIfMissing)
        {
            string strScoreValue = null;
            double dblScore = 0;

            strScoreValue = GetScore(strScoreName);
            if (!string.IsNullOrEmpty(strScoreValue) && double.TryParse(strScoreValue, out dblScore))
            {
                return dblScore;
            }
            else
            {
                return dblValueIfMissing;
            }
        }

        /// <summary>
        ///  Returns the value stored for the specified score (as an integer)
        /// </summary>
        /// <param name="strScoreName">Score name</param>
        /// <returns>Score if defined, otherwise 0</returns>
        /// <remarks></remarks>
        public int GetScoreInt(string strScoreName)
        {
            return GetScoreInt(strScoreName, 0);
        }

        /// <summary>
        ///  Returns the value stored for the specified score (as an integer)
        /// </summary>
        /// <param name="strScoreName">Score name</param>
        /// <param name="intValueIfMissing">Value to return if the score is not defined</param>
        /// <returns>Score if defined, otherwise intValueIfMissing</returns>
        /// <remarks></remarks>
        public int GetScoreInt(string strScoreName, int intValueIfMissing)
        {
            string strScoreValue = null;
            var intScore = 0;

            strScoreValue = GetScore(strScoreName);
            if (!string.IsNullOrEmpty(strScoreValue) && int.TryParse(strScoreValue, out intScore))
            {
                return intScore;
            }
            else
            {
                return intValueIfMissing;
            }
        }

        public void SetPeptide(string strPeptide)
        {
            SetPeptide(strPeptide, blnUpdateCleanSequence: true);
        }

        /// <summary>
        /// Update the peptide sequence, auto-determining the clean sequence if blnUpdateCleanSequence is true
        /// </summary>
        /// <param name="strPeptide">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
        /// <param name="blnUpdateCleanSequence"></param>
        /// <remarks>Does not update the cleavage state info.  If blnUpdateCleanSequence is false, then call UpdateCleanSequence at a later time to populate mPeptideCleanSequence</remarks>
        public void SetPeptide(string strPeptide, bool blnUpdateCleanSequence)
        {
            if (string.IsNullOrEmpty(strPeptide))
            {
                mPeptide = string.Empty;
            }
            else
            {
                mPeptide = strPeptide;
            }

            if (blnUpdateCleanSequence)
            {
                UpdateCleanSequence(mPeptide);
            }
        }

        /// <summary>
        /// Update the peptide sequence (auto-determines the clean sequence); also auto-update the the cleavage state info
        /// </summary>
        /// <param name="strPeptide">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
        /// <param name="objCleavageStateCalculator">Cleavage state calculator object</param>
        /// <remarks></remarks>
        public void SetPeptide(string strPeptide, clsPeptideCleavageStateCalculator objCleavageStateCalculator)
        {
            SetPeptide(strPeptide);
            UpdateCleavageInfo(objCleavageStateCalculator);
        }

        /// <summary>
        /// Add/update an additional score to associate with this peptide
        /// </summary>
        /// <param name="strScoreName"></param>
        /// <param name="strScoreValue"></param>
        /// <remarks></remarks>
        public void SetScore(string strScoreName, string strScoreValue)
        {
            if (mAdditionalScores.ContainsKey(strScoreName))
            {
                mAdditionalScores[strScoreName] = strScoreValue;
            }
            else
            {
                mAdditionalScores.Add(strScoreName, strScoreValue);
            }
        }

        /// <summary>
        /// Returns the value stored for the specified score
        /// </summary>
        /// <param name="strScoreName"></param>
        /// <param name="strScoreValue"></param>
        /// <returns>True if the score is defined, otherwise false</returns>
        /// <remarks></remarks>
        public bool TryGetScore(string strScoreName, out string strScoreValue)
        {
            strScoreValue = string.Empty;
            if (mAdditionalScores.TryGetValue(strScoreName, out strScoreValue))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// Auto-determine the number of missed cleavages, cleavage state, and number of tryptic terminii based on the peptide sequence
        /// </summary>
        /// <param name="objCleavageStateCalculator"></param>
        /// <remarks></remarks>
        public void UpdateCleavageInfo(clsPeptideCleavageStateCalculator objCleavageStateCalculator)
        {
            NumMissedCleavages = objCleavageStateCalculator.ComputeNumberOfMissedCleavages(mPeptide);

            CleavageState = objCleavageStateCalculator.ComputeCleavageState(mPeptide);

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
