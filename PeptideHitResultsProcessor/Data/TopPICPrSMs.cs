namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class holds rows read from the tab-delimited file (_TopPIC_PrSMs.txt) created directly by TopPIC
    /// </summary>
    /// <remarks>
    /// This class is only used by TopPICResultsProcessor
    /// </remarks>
    internal class TopPICPrSMs
    {
        // Ignore Spelling: proteoform

        // ReSharper disable NotAccessedField.Local

        public string SpectrumFileName;
        public string Prsm_ID;
        public string Spectrum_ID;
        public string FragMethod;
        public string Scans;
        public int ScanNum;
        public string RetentionTime;
        public string Peaks;
        public string Charge;
        public short ChargeNum;
        public string Precursor_mass;               // Monoisotopic mass value of the observed precursor_mz
        public string PrecursorMZ;                  // Computed from Precursor_mass
        public string Adjusted_precursor_mass;      // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC; typically identical to Proteoform_mass
        public string MH;                           // Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
        public string DelM;                         // Computed using Precursor_mass - Adjusted_precursor_mass
        public string DelM_PPM;                     // Computed using DelM and Adjusted_precursor_mass
        public string Proteoform_ID;
        public string Feature_Intensity;
        public string Feature_Score;
        public string Feature_Apex_Time;
        public string Protein_Hits;
        public string Protein;
        public string ProteinDescription;
        public string ResidueStart;                 // First_residue
        public string ResidueEnd;                   // Last_residue
        public string Special_amino_acids;
        public string Proteoform;
        public string Protein_Nterminal_Form;
        public string Proteoform_mass;              // Theoretical monoisotopic mass of the peptide (including mods), as computed by TopPIC; typically identical to Adjusted_precursor_mass
        public string Unexpected_Mod_Count;         // unexpected modifications
        public string MIScore;
        public string VariablePTMs;
        public string Matched_peaks;
        public string Matched_fragment_ions;
        public string PValue;
        public double PValueNum;
        public int RankPValue;
        public string Evalue;
        public double EValueNum;
        public string Qvalue;
        public string Proteoform_QValue;

        public bool StoreInFirstHitsFile;

        // ReSharper restore NotAccessedField.Local

        /// <summary>
        /// Constructor
        /// </summary>
        public TopPICPrSMs()
        {
            Clear();
        }

        /// <summary>
        /// Reset stored values to empty strings and zeros
        /// </summary>
        // ReSharper disable once UnusedMember.Local
        public void Clear()
        {
            SpectrumFileName = string.Empty;
            Prsm_ID = string.Empty;
            Spectrum_ID = string.Empty;
            FragMethod = string.Empty;
            Scans = string.Empty;
            ScanNum = 0;
            RetentionTime = string.Empty;
            Peaks = string.Empty;
            Charge = string.Empty;
            ChargeNum = 0;
            Precursor_mass = string.Empty;
            PrecursorMZ = string.Empty;
            Adjusted_precursor_mass = string.Empty;
            MH = string.Empty;
            DelM = string.Empty;
            DelM_PPM = string.Empty;
            Proteoform_ID = string.Empty;
            Feature_Intensity = string.Empty;
            Feature_Score = string.Empty;
            Feature_Apex_Time = string.Empty;
            Protein_Hits = string.Empty;
            Protein = string.Empty;
            ProteinDescription = string.Empty;
            ResidueStart = string.Empty;
            ResidueEnd = string.Empty;
            Special_amino_acids = string.Empty;
            Proteoform = string.Empty;
            Protein_Nterminal_Form = string.Empty;
            Proteoform_mass = string.Empty;
            Unexpected_Mod_Count = string.Empty;
            MIScore = string.Empty;
            VariablePTMs = string.Empty;
            Matched_peaks = string.Empty;
            Matched_fragment_ions = string.Empty;
            PValue = string.Empty;
            PValueNum = 0;
            RankPValue = 0;
            Evalue = string.Empty;
            EValueNum = 0;
            Qvalue = string.Empty;
            Proteoform_QValue = string.Empty;
        }

        /// <summary>
        /// Show scan, proteoform, and p-value
        /// </summary>
        public override string ToString()
        {
            return string.Format("Scan {0}: {1}, PValue {2}", ScanNum, Proteoform, PValue);
        }
    }
}
