using PHRPReader;
using PHRPReader.Data;
using System.Collections.Generic;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for a MSFragger search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class MSFraggerResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Da, Fragger, Hyperscore, Nextscore, NTT, tryptic

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
        /// Precursor ion m/z (observed)
        /// </summary>
        public string PrecursorMZ { get; set; }

        /// <summary>
        /// Mass error, in Da, as computed by PHRP
        /// </summary>
        public string PHRPComputedDelM { get; set; }

        /// <summary>
        /// Mass error, in ppm, as computed by PHRP
        /// </summary>
        public string PHRPComputedDelMPPM { get; set; }

        /// <summary>
        /// Mass error, in Da, as computed by MSFragger
        /// </summary>
        /// <remarks>
        /// <para>
        /// Mass difference between "Calibrated Observed Mass" and "Calculated Peptide Mass", as reported by MSFragger
        /// </para>
        /// <para>
        /// From column "Delta Mass"
        /// </para>
        /// </remarks>
        public string MSFraggerComputedDelM { get; set; }

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MSFragger
        /// </summary>
        public string CalculatedMonoMass { get; set; }

        /// <summary>
        /// Comma-separated list of static and dynamic modification names and affected residue number
        /// </summary>
        /// <remarks>
        /// Example values:
        ///   15M(15.9949)
        ///   8C(57.0215)
        ///   5M(15.9949), 8M(15.9949)
        ///   1M(15.9949), 5C(57.0215)
        ///   12M(15.9949), 2M(15.9949), 8M(15.9949)
        ///   N-term(42.0106)
        /// </remarks>
        public string Modifications { get; set; }

        /// <summary>
        /// Name of the best scoring protein this peptide is associated with
        /// </summary>
        public string Protein { get; set; }

        /// <summary>
        /// Proteins associated with this peptide
        /// </summary>
        /// <remarks>
        /// This list is populated using the protein name in the Protein column,
        /// plus any other proteins in the AdditionalProteins column
        /// </remarks>
        public List<string> Proteins { get; } = new();

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        public string NTT { get; set; }

        /// <summary>
        /// Expectation score (E-value)
        /// </summary>
        /// <remarks>
        /// Smaller values (closer to zero) are higher confidence
        /// </remarks>
        public string EValue { get; set; }

        /// <summary>
        /// Hyperscore of the top match
        /// </summary>
        /// <remarks>
        /// Larger scores are better
        /// </remarks>
        public string Hyperscore { get; set; }

        /// <summary>
        /// Hyperscore of the next best match
        /// </summary>
        /// <remarks>
        /// Larger values are better
        /// </remarks>
        public string Nextscore { get; set; }

        /// <summary>
        /// Peptide prophet probability
        /// </summary>
        /// <remarks>
        /// Values closer to 1 are better
        /// </remarks>
        public string PeptideProphetProbability { get; set; }

        /// <summary>
        /// Elution time of the MS/MS spectrum
        /// </summary>
        public string ElutionTime { get; set; }

        /// <summary>
        /// Average elution time of MS/MS spectra (only applicable for a multi-dataset based search)
        /// </summary>
        public string ElutionTimeAverage { get; set; }

        /// <summary>
        /// Number of missed enzymatic cleavages
        /// </summary>
        public string MissedCleavageCount { get; set; }

        /// <summary>
        /// Number of matched ions in the fragmentation spectrum
        /// </summary>
        /// <remarks>
        /// Read from the Dataset.tsv file since not in the _psm.tsv file
        /// </remarks>
        public string NumberOfMatchedIons { get; set; }

        /// <summary>
        /// Total number of ions in the fragmentation spectrum
        /// </summary>
        /// <remarks>
        /// Read from the Dataset.tsv file since not in the _psm.tsv file
        /// </remarks>
        public string TotalNumberOfIons { get; set; }

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
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public MSFraggerResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
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
            PrecursorMZ = string.Empty;
            PHRPComputedDelM = string.Empty;
            PHRPComputedDelMPPM = string.Empty;
            MSFraggerComputedDelM = string.Empty;
            CalculatedMonoMass = string.Empty;
            Modifications = string.Empty;
            Protein = string.Empty;
            Proteins.Clear();
            NTT = string.Empty;
            EValue = string.Empty;
            Hyperscore = string.Empty;
            Nextscore = string.Empty;
            PeptideProphetProbability = string.Empty;
            ElutionTime = string.Empty;
            ElutionTimeAverage = string.Empty;
            MissedCleavageCount = string.Empty;
            NumberOfMatchedIons = string.Empty;
            TotalNumberOfIons = string.Empty;
            QValue = string.Empty;
        }
    }
}
