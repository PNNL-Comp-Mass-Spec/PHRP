﻿namespace PeptideHitResultsProcessor.SearchToolResults
{
    /// <summary>
    /// This class is used to track search results loaded from a MaxQuant, MSFragger, or DIA-NN result file
    /// </summary>
    public abstract class ToolResultsBaseClass
    {
        // Ignore Spelling: Acetyl, Da, prev, tryptic

        /// <summary>
        /// Dataset name
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, the .raw file name, without the extension
        /// </para>
        /// <para>
        /// For MSFragger, comes from column "Spectrum File" in Dataset_psm.tsv (typically extracted from interact-Dataset.pep.xml)
        /// </para>
        /// <para>
        /// For MSFragger, the _psm.tsv file will typically only show this for the first result for each dataset,
        /// but MSFraggerResultsProcessor stores the correct name for each PSM
        /// </para>
        /// <para>
        /// For DIA-NN, comes from column Run in the Dataset_report.tsv file
        /// This is an abbreviated dataset name, as assigned by the analysis manager
        /// </para>
        /// </remarks>
        public string DatasetName;

        /// <summary>
        /// MS/MS scan number in the dataset
        /// </summary>
        /// <remarks>
        /// For MSFragger, comes from column "Spectrum"
        /// </remarks>
        public string Scan;

        /// <summary>
        /// Integer value of Scan
        /// </summary>
        public int ScanNum;

        /// <summary>
        /// Identified peptide
        /// </summary>
        /// <remarks>
        /// <para>
        /// Amino acid symbols only; no modification info
        /// </para>
        /// <para>
        /// For MSFragger, from column "Peptide"
        /// </para>
        /// <para>
        /// For DIA-NN, originally from column "Stripped.Sequence", and thus no prefix or suffix residues
        /// However, method AddPrefixAndSuffixResiduesUsingFASTA uses clsPeptideToProteinMapEngine to determine the prefix and suffix residues for each peptide
        /// </para>
        /// </remarks>
        public string Sequence;

        /// <summary>
        /// Residue in the protein before this peptide
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, from the peptides.txt file, column "Amino acid before"
        /// </para>
        /// <para>
        /// For MSFragger, from column "Prev AA"
        /// </para>
        /// </remarks>
        public string PrefixResidue;

        /// <summary>
        /// Residue in the protein after this peptide
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, from the peptides.txt file, column "Amino acid after"
        /// </para>
        /// <para>
        /// For MSFragger, from column "Next AA"
        /// </para>
        /// </remarks>
        public string SuffixResidue;

        /// <summary>
        /// Number of residues in the peptide
        /// </summary>
        public int Length;

        // ReSharper disable CommentTypo

        /// <summary>
        /// Post-translational modifications contained within the identified peptide
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, this is read from the msms.txt file and includes dynamic modifications, but not static modifications
        /// </para>
        /// <para>
        /// For MSFragger, this includes both static and dynamic modifications
        /// and comes from column "Assigned Modifications"
        /// </para>
        /// <para>
        /// For DIA-NN, this is read from column Modified.Sequence in Dataset_report.tsv
        /// </para>
        /// <para>
        /// MaxQuant Examples:
        /// One modified residue:
        ///   Oxidation (M)
        /// Two residues with the same modification
        ///   2 Oxidation (M)
        /// Two separate modifications (comma-separated list)
        ///   Acetyl (Protein N-term),Oxidation (M)
        /// </para>
        /// <para>
        /// MSFragger Examples:
        ///   15M(15.9949)
        ///   8C(57.0215)
        ///   11C(57.0215), 4C(57.0215)
        ///   5M(15.9949), 8M(15.9949)
        ///   1M(15.9949), 5C(57.0215)
        ///   12M(15.9949), 2M(15.9949), 8M(15.9949)
        ///   N-term(42.0106)
        ///   6M(15.9949), N-term(42.0106)
        /// </para>
        /// <para>
        /// DIA-NN Examples:
        ///   PHRP converts "AAAGDLGGDHLAFSC(UniMod:4)DVAK"            to "15C(57.0215)"
        ///   PHRP converts "AAGVGLVDC(UniMod:4)HC(UniMod:4)HLSAPDFDR" to "9C(57.0215), 11C(57.0215)"
        ///   PHRP converts "DIHFM(UniMod:35)PC(UniMod:4)SGLTGANLK"    to "5M(15.9949), 7C(57.0215)"
        /// </para>
        /// </remarks>
        public string ModificationList;

        // ReSharper restore CommentTypo

        /// <summary>
        /// Charge state of the precursor ion
        /// </summary>
        public string Charge;

        /// <summary>
        /// Numeric value of Charge
        /// </summary>
        public short ChargeNum;

        /// <summary>
        /// Elution time of the MS/MS spectrum, in minutes
        /// </summary>
        /// <remarks>
        /// For MaxQuant, in the msms.txt file, column "Retention time" is in minutes
        /// For MSFragger, in the Dataset_psm.tsv file, column "Retention" is in seconds
        /// For MSFragger, in the MSFragger results file, column "retention_time" is in minutes
        /// For DIA-NN, in the Dataset_report.tsv file, column "RT" is in minutes
        /// </remarks>
        public string ElutionTime;

        /// <summary>
        /// Monoisotopic (M+H)+ value, computed from PrecursorMZ and Charge
        /// </summary>
        public string MH;

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as reported by MaxQuant or MSFragger
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, this mass will be overridden by PHRP if isobaric mods were used
        /// </para>
        /// <para>
        /// For DIA-NN, this value is computed by PHRP
        /// </para>
        /// </remarks>
        public string CalculatedMonoMass;

        /// <summary>
        /// Numeric value of CalculatedMonoMass
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, this mass will be overridden by PHRP if isobaric mods were used
        /// </para>
        /// <para>
        /// For MSFragger, from column "Calculated Peptide Mass"
        /// </para>
        /// <para>
        /// For DIA-NN, this value is computed by PHRP and is identical to CalculatedMonoMassPHRP
        /// </para>
        /// </remarks>
        public double CalculatedMonoMassValue;

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by PHRP
        /// </summary>
        public double CalculatedMonoMassPHRP;

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence,
        /// as computed by PHRP
        /// </summary>
        /// <remarks>
        /// This is a C13-corrected precursor error
        /// </remarks>
        public string MassErrorPpm;

        /// <summary>
        /// Mass error, in Da, as computed by PHRP
        /// </summary>
        /// <remarks>
        /// This is a C13-corrected precursor error
        /// </remarks>
        public string MassErrorDa;

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        /// <remarks>
        /// For MaxQuant, initially read from column "Number of Enzymatic Termini", but may be updated by PHRP
        /// </remarks>
        public int NumberOfTrypticTermini;

        /// <summary>
        /// Number of missed enzymatic cleavages
        /// </summary>
        /// <remarks>
        /// For MSFragger, from column "Number of Missed Cleavages", but may be updated by PHRP (since sometimes wrong for _psm.tsv files)
        /// </remarks>
        public string MissedCleavageCount;

        /// <summary>
        /// True for peptides that are associated with a decoy protein
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, the msms.txt file has '+' in the Reverse column for decoy peptides
        /// </para>
        /// <para>
        /// For MSFragger, determined based on protein name
        /// </para>
        /// <para>
        /// Not stored in the synopsis file
        /// </para>
        /// </remarks>
        public bool Reverse;

        /// <summary>
        /// Protein name (best name, if more than one protein)
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, this is the Leading Razor Protein
        /// </para>
        /// <para>
        /// For MSFragger, comes from the Protein column
        /// </para>
        /// <para>
        /// For DIA-NN, the first protein in the Protein.Ids column
        /// </para>
        /// </remarks>
        public string Protein;

        /// <summary>
        /// Residue number in the protein where this peptide starts
        /// </summary>
        public string ProteinStart;

        /// <summary>
        /// Residue number in the protein where this peptide ends
        /// </summary>
        public string ProteinEnd;

        /// <summary>
        /// Peptide intensity
        /// </summary>
        /// <remarks>
        /// <para>
        /// For MaxQuant, from the peptides.txt file, column "Intensity", representing the
        /// summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
        /// </para>
        /// <para>
        /// For MSFragger, Often "0"
        /// </para>
        /// </remarks>
        public string Intensity;

        /// <summary>
        /// Score rank (1 is best match, 2 is second best, etc.)
        /// </summary>
        public int RankScore;

        /// <summary>
        /// FDR
        /// </summary>
        /// <remarks>
        /// Computed by PHRP for MSFragger, MaxQuant, and DIA-NN
        /// </remarks>
        public double FDR;

        /// <summary>
        /// Q-Value
        /// </summary>
        /// <remarks>
        /// <para>
        /// Smaller values (closer to zero) are higher confidence
        /// </para>
        /// <para>
        /// Computed by PHRP for MSFragger and MaxQuant
        /// </para>
        /// <para>
        /// From column Q.Value for DIA-NN
        /// </para>
        /// <para>
        /// The MSGF Results Summarizer filters on QValue when computing PSM stats, which are stored in the database using procedure store_job_psm_stats
        /// Stats are displayed on the Analysis Job Detail Report page, for example:
        /// Spectra Searched: 23949, Total PSMs: 12317, Unique Peptides: 10439, Unique Proteins: 2247 (FDR &lt; 1.00%)
        /// </para>
        /// </remarks>
        public double QValue;
    }
}
