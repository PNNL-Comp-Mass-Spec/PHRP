using PHRPReader;
using PHRPReader.Data;
using System.Collections.Generic;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for a MaxQuant search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class MaxQuantResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Acetyl, Carbamidomethyl, Da, MaxQuant, Orbitrap, plex, tryptic

        /// <summary>
        /// Dataset name
        /// </summary>
        /// <remarks>This is an abbreviated form of the full dataset name</remarks>
        public string DatasetName { get; set; }

        /// <summary>
        /// Dataset ID
        /// </summary>
        public int DatasetID { get; set; }

        /// <summary>
        /// Type of fragmentation used to create the MS/MS spectrum
        /// </summary>
        /// <remarks>
        /// Types:
        ///   CID - Collision Induced Dissociation
        ///   HCD - High energy Collision induced Dissociation
        ///   ETD - Electron Transfer Dissociation
        /// </remarks>
        public string FragMethod { get; set; }

        /// <summary>
        /// Index of the spectrum in the dataset (1-based, consecutive integer)
        /// </summary>
        public string SpecIndex { get; set; }

        /// <summary>
        /// Precursor ion m/z (observed value, read from a _PrecursorInfo.txt file by PHRP)
        /// </summary>
        public string PrecursorMZ { get; set; }

        // ReSharper disable CommentTypo

        /// <summary>
        /// Precursor ion m/z (theoretical value, not observed value), as reported by MaxQuant
        /// </summary>
        /// <remarks>
        /// <para>
        /// This is the theoretical m/z of the first isotope of the isotopic distribution of the parent ion
        /// It does not account for isobaric mods
        /// </para>
        /// <para>
        /// For example, given a 3+ parent ion whose isotopic distribution in the MS1 spectrum
        /// has ions at 827.769, 828.105, 828.438, and 828.772 m/z, the most intense ion is the 828.434 ion.
        /// The instrument might then isolate that ion for fragmentation, giving a spectrum label for the MS2 spectrum of
        /// "Full ms2 828.44@cid35.00"
        /// </para>
        /// <para>
        /// MS-GF+ reports the precursor m/z as 827.76965, which is the observed m/z value
        /// MaxQuant reports the PrecursorMZ as 827.76176, which is the theoretical value based on the identified peptide (TAHEVRPGNVIMFEGSPWVVQK)
        /// </para>
        /// <para>
        /// Example 2: given a 3+ parent ion whose isotopic distribution in the MS1 spectrum
        /// has ions at 755.402, 755.737, 756.071, and 756.407 m/z, the most intense ion is the 755.737 ion.
        /// The instrument might then isolate that ion for fragmentation, giving a spectrum label for the MS2 spectrum of
        /// "Full ms2 755.74@cid35.00"
        /// </para>
        /// <para>
        /// MS-GF+ reports the precursor m/z as 755.40515, which is the observed m/z value (after correction by MZ Refinery)
        /// MaxQuant reports the PrecursorMZ as 755.39497, which is the theoretical value based on the identified peptide (IINIGSVVGTMGNAGQVNYSAAK)
        /// </para>
        /// <para>
        /// Isobaric Label Example: given a 3+ parent ion whose isotopic distribution in the MS1 spectrum
        /// has ions at 902.844, 903.176, 903.512, and 903.843 m/z, the most intense ion is the 903.17 ion.
        /// The instrument might then isolate that ion for fragmentation, giving a spectrum label for the MS2 spectrum of
        /// "Full ms2 903.18@hcd28.00"
        /// </para>
        /// <para>
        /// MS-GF+ reports the precursor m/z as 903.17218, which is the observed m/z value (after correction by MZ Refinery)
        /// MaxQuant reports the PrecursorMZ as 750.0655, which is the theoretical value based on the monoisotopic mass of the peptide (ILSLLEGQKIAFGGETDEATR),
        /// not counting the isobaric 6-plex TMT mod that is present at both the N-terminus and on the internal K residue
        /// Mono mass (without actual mods) = 2247.174619
        /// Charge = 3
        /// Reported m/z = (2247.174619 / 3) + 1.007276467, where 1.00727 is the charge carrier mass (one proton)
        /// </para>
        /// <para>
        /// Correct mono mass = 2705.5004
        /// Correct m/z calculation = (2705.5004 + 3 * 1.007276467) / 3 = 902.8407495
        /// </para>
        /// </remarks>
        public string PrecursorMZ_MaxQuant { get; set; }

        // ReSharper restore CommentTypo

        /// <summary>
        /// Mass error, in Da, as computed by PHRP
        /// </summary>
        public string PHRPComputedDelM { get; set; }

        /// <summary>
        /// Mass error, in ppm, as computed by PHRP
        /// </summary>
        public string PHRPComputedDelMPPM { get; set; }

        /// <summary>
        /// Mass error, in Da, as computed by MaxQuant
        /// </summary>
        /// <remarks>
        /// This is an empty string if a peptide has Type=MSMS
        /// In contrast, if a peptide has Type=MULTI-MSMS, a value will be defined
        /// </remarks>
        public string MaxQuantComputedDelM { get; set; }

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence,
        /// as computed by MaxQuant
        /// </summary>
        /// <remarks>
        /// This is an empty string if a peptide has Type=MSMS
        /// In contrast, if a peptide has Type=MULTI-MSMS, a value will be defined
        /// </remarks>
        public string MaxQuantComputedDelMPPM { get; set; }

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MaxQuant
        /// </summary>
        public string CalculatedMonoMass { get; set; }

        /// <summary>
        /// Comma-separated list of dynamic modification names and affected residue number
        /// </summary>
        /// <remarks>
        /// Example values:
        ///   Oxidation 7
        ///   Acetyl 1,Oxidation 4
        /// </remarks>
        public string Modifications { get; set; }

        /// <summary>
        /// Protein(s) associated with this peptide
        /// </summary>
        /// <remarks>
        /// Empty string for peptides resulting from a reverse hit protein;
        /// the reverse-hit protein name will be in LeadingRazorProtein
        /// </remarks>
        public List<string> Proteins { get; } = new();

        /// <summary>
        /// Name of the best scoring protein this peptide is associated with
        /// </summary>
        /// <remarks>
        /// Typically there is only one protein name here
        /// However, in cases of a tied score, will be a semicolon separated list
        /// </remarks>
        public string LeadingRazorProtein { get; set; }

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        public string NTT { get; set; }

        /// <summary>
        /// Posterior error probability
        /// </summary>
        /// <remarks>
        /// Similar to p-value
        /// Smaller values (closer to zero) are higher confidence
        /// </remarks>
        public string PEP { get; set; }

        /// <summary>
        /// Andromeda score for the best MS/MS spectrum with this peptide
        /// </summary>
        /// <remarks>
        /// Larger values are better
        /// </remarks>
        public string Score { get; set; }

        /// <summary>
        /// Score difference to the second best identified peptide with a different amino acid sequence
        /// </summary>
        public string DeltaScore { get; set; }

        /// <summary>
        /// Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
        /// </summary>
        /// <remarks>
        /// <para>
        /// All peptides with the same sequence will have the same total peptide intensity value for a given dataset
        /// </para>
        /// <para>
        /// From the peptides.txt file, column "Intensity"
        /// </para>
        /// </remarks>
        public string TotalPeptideIntensity { get; set; }

        /// <summary>
        /// Mass Analyzer of the instrument
        /// </summary>
        /// <remarks>
        /// Types:
        ///   ITMS - Ion trap
        ///   FTMS - Fourier transform ICR or Orbitrap
        ///   TOF - Time of flight
        /// </remarks>
        public string MassAnalyzer { get; set; }

        /// <summary>
        /// Type of precursor ion as identified by MaxQuant
        /// </summary>
        /// <remarks>
        /// <para>
        /// Prefixes:
        ///   ISO - isotopic cluster.
        ///   PEAK - single peak.
        ///   MULTI - labeling cluster.
        /// </para>
        /// <para>
        /// Example values:
        ///   MULTI-MSMS
        ///   MULTI-SECPEP
        ///   MSMS
        /// </para>
        /// </remarks>
        public string PrecursorType { get; set; }

        /// <summary>
        /// Elution time of the MS/MS spectrum
        /// </summary>
        public string RetentionTime { get; set; }

        /// <summary>
        /// Scan number where the precursor ion was observed
        /// </summary>
        public string PrecursorScanNumber { get; set; }

        /// <summary>
        /// Intensity of the precursor ion in the scan that it was observed
        /// </summary>
        public string PrecursorIntensity { get; set; }

        /// <summary>
        /// Number of peaks (MS/MS ions) matching to the predicted fragmentation spectrum
        /// </summary>
        public string NumberOfMatches { get; set; }

        /// <summary>
        /// Fraction of intensity in the MS/MS spectrum that is annotated
        /// </summary>
        public string IntensityCoverage { get; set; }

        /// <summary>
        /// Number of missed enzymatic cleavages
        /// </summary>
        public string MissedCleavageCount { get; set; }

        /// <summary>
        /// Unique (consecutive) identifier for each row in msms.txt
        /// </summary>
        /// <remarks>
        /// Used to cross-link the information in msms.txt with information stored in other files
        /// </remarks>
        public string MsMsID { get; set; }

        /// <summary>
        /// Identifier of the protein-group this redundant peptide sequence is associated with
        /// </summary>
        /// <remarks>
        /// <para>
        /// Typically a single number, but could be a semicolon separated list
        /// if a peptide is associated with multiple protein groups
        /// </para>
        /// <para>
        /// Can be used to look up the extended protein information in the proteinGroups.txt file
        /// </para>
        /// </remarks>
        public string ProteinGroupIDs { get; set; }

        /// <summary>
        /// The identifier of the non-redundant peptide sequence
        /// </summary>
        /// <remarks>
        /// Corresponds to the id column in the peptides.txt file
        /// </remarks>
        public string PeptideID { get; set; }

        /// <summary>
        /// Identifier referencing a row in the modificationSpecificPeptides.txt file
        /// </summary>
        public string ModPeptideID { get; set; }

        /// <summary>
        /// Identifier referencing a row in the evidence.txt file
        /// </summary>
        public string EvidenceID { get; set; }

        /// <summary>
        /// Q-Value
        /// </summary>
        /// <remarks>
        /// Computed by PHRP
        /// </remarks>
        public string QValue { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods"></param>
        /// <param name="peptideSeqMassCalculator"></param>
        public MaxQuantResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
            : base(peptideMods, peptideSeqMassCalculator)
        {
        }

        /// <summary>
        /// Clear stored values
        /// </summary>
        public override void Clear()
        {
            base.Clear();

            DatasetName = string.Empty;
            DatasetID = 0;
            FragMethod = string.Empty;
            SpecIndex = string.Empty;
            PrecursorMZ = string.Empty;
            PrecursorMZ_MaxQuant = string.Empty;
            PHRPComputedDelM = string.Empty;
            PHRPComputedDelMPPM = string.Empty;
            MaxQuantComputedDelM = string.Empty;
            MaxQuantComputedDelMPPM = string.Empty;
            CalculatedMonoMass = string.Empty;
            Modifications = string.Empty;
            Proteins.Clear();
            LeadingRazorProtein = string.Empty;
            NTT = string.Empty;
            PEP = string.Empty;
            Score = string.Empty;
            DeltaScore = string.Empty;
            TotalPeptideIntensity = string.Empty;
            MassAnalyzer = string.Empty;
            PrecursorType = string.Empty;
            RetentionTime = string.Empty;
            PrecursorScanNumber = string.Empty;
            PrecursorIntensity = string.Empty;
            NumberOfMatches = string.Empty;
            IntensityCoverage = string.Empty;
            MissedCleavageCount = string.Empty;
            MsMsID = string.Empty;
            ProteinGroupIDs = string.Empty;
            PeptideID = string.Empty;
            ModPeptideID = string.Empty;
            EvidenceID = string.Empty;
            QValue = string.Empty;
        }
    }
}
