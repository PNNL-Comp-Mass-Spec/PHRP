namespace PeptideHitResultsProcessor.SearchToolResults
{
    /// <summary>
    /// This class holds rows read from the MSFragger PSM results file
    /// </summary>
    /// <remarks>
    /// It is used when creating the synopsis file
    /// </remarks>
    internal class MSFraggerSearchResult : ToolResultsBaseClass
    {
        // Ignore Spelling: Hyperscore, sp

        /// <summary>
        /// Average elution time of MS/MS spectra when reading results from a multi-dataset search, in minutes
        /// </summary>
        public string ElutionTimeAverage;

        // ReSharper disable CommentTypo

        /// <summary>
        /// Peptide sequence with integer values representing the total mass of the any modified residue
        /// </summary>
        /// <remarks>
        /// <para>
        /// Example values:
        ///   GASQAGM[147]TAPGTK
        ///   n[43]GEDAAQAEK
        ///   HQGVM[147]VGM[147]GQK
        /// </para>
        /// <para>
        /// Empty string if no dynamic modifications
        /// </para>
        /// </remarks>
        public string ModifiedPeptide;

        // ReSharper restore CommentTypo

        /// <summary>
        /// Monoisotopic mass of the peptide, as reported by MSFragger
        /// </summary>
        /// <remarks>
        /// <para>
        /// Computed from PrecursorMZ and charge
        /// </para>
        /// <para>
        /// From column "Observed Mass" in _psm.tsv files
        /// From column precursor_neutral_mass in Dataset.tsv files
        /// </para>
        /// </remarks>
        public string PrecursorMonoMass;

        /// <summary>
        /// Calibrated observed monoisotopic mass of the peptide, as reported by MSFragger
        /// </summary>
        /// <remarks>
        /// From column "Calibrated Observed Mass"
        /// </remarks>
        public string CalibratedObservedMass;

        /// <summary>
        /// Precursor ion m/z (observed)
        /// </summary>
        /// <remarks>
        /// From column "Observed M/Z" in _psm.tsv files
        /// Computed using precursor_neutral_mass and charge for Dataset.tsv files
        /// </remarks>
        public string PrecursorMZ;

        /// <summary>
        /// Calibrated observed precursor m/z, as reported by MSFragger
        /// </summary>
        /// <remarks>
        /// From column "Calibrated Observed M/Z"
        /// </remarks>
        public string CalibratedObservedMZ;

        /// <summary>
        /// Expected m/z of the identified peptide, as reported by MSFragger
        /// </summary>
        /// <remarks>
        /// From column "Calculated M/Z"
        /// </remarks>
        public string CalculatedMZ;

        /// <summary>
        /// Mass difference between "Calibrated Observed Mass" and "Calculated Peptide Mass", as reported by MSFragger
        /// </summary>
        /// <remarks>
        /// From column "Delta Mass"
        /// </remarks>
        public string MassErrorDaMSFragger;

        /// <summary>
        /// Expectation value (E-value)
        /// </summary>
        /// <remarks>
        /// Smaller values (closer to zero) are better
        /// </remarks>
        public string Expectation;

        /// <summary>
        /// Numeric value of Expectation
        /// </summary>
        public double EValue;

        /// <summary>
        /// Hyperscore of the top match
        /// </summary>
        /// <remarks>
        /// Larger values are better
        /// </remarks>
        public string Hyperscore;

        /// <summary>
        /// Numeric value of Hyperscore
        /// </summary>
        public double HyperscoreValue;

        /// <summary>
        /// Hyperscore of the next best match
        /// </summary>
        public string Nextscore;

        /// <summary>
        /// Peptide prophet probability
        /// </summary>
        /// <remarks>
        /// Values closer to 1 are better
        /// </remarks>
        public string PeptideProphetProbability;

        /// <summary>
        /// Observed modifications (typically blank)
        /// </summary>
        /// <remarks>
        /// From column "Observed Modifications"
        /// </remarks>
        public string ObservedModifications;

        /// <summary>
        /// Protein uniqueness flag
        /// </summary>
        /// <remarks>
        /// <para>
        /// "true" if the peptide only maps to one protein, otherwise "false"
        /// </para>
        /// <para>
        /// From column "Is Unique"
        /// </para>
        /// </remarks>
        public string IsUnique;

        /// <summary>
        /// Protein ID (typically blank)
        /// </summary>
        /// <remarks>
        /// From column "Protein ID"
        /// </remarks>
        public string ProteinID;

        // ReSharper disable CommentTypo

        /// <summary>
        /// Protein entry name parsed from protein name
        /// </summary>
        /// <remarks>
        /// <para>
        /// For example, if Protein is "sp|P30043|BLVRB_HUMAN", EntryName will be "BLVRB_HUMAN"
        /// </para>
        /// <para>
        /// From column "Entry Name"
        /// </para>
        /// </remarks>
        public string EntryName;

        /// <summary>
        /// Gene name, parsed from protein name
        /// </summary>
        /// <remarks>
        /// <para>
        /// For example, if Protein is "sp|P30043|BLVRB_HUMAN", Gene will be "BLVRB"
        /// </para>
        /// <para>
        /// From column "Gene"
        /// </para>
        /// </remarks>
        public string Gene;

        // ReSharper restore CommentTypo

        /// <summary>
        /// Protein description
        /// </summary>
        public string ProteinDescription;

        /// <summary>
        /// Additional genes (comma-separated list)
        /// </summary>
        /// <remarks>
        /// From column "Mapped Genes"
        /// </remarks>
        public string AdditionalGenes;

        /// <summary>
        /// Additional protein names (comma-separated list)
        /// </summary>
        /// <remarks>
        /// From column "Mapped Proteins"
        /// </remarks>
        public string AdditionalProteins;

        /// <summary>
        /// Number of matched ions in the fragmentation spectrum
        /// </summary>
        /// <remarks>
        /// Read from the Dataset.tsv file since not in the _psm.tsv file
        /// </remarks>
        public string NumberOfMatchedIons;

        /// <summary>
        /// Total number of ions in the fragmentation spectrum
        /// </summary>
        /// <remarks>
        /// Read from the Dataset.tsv file since not in the _psm.tsv file
        /// </remarks>
        public string TotalNumberOfIons;

        /// <summary>
        /// Reset stored values to empty strings and zeros
        /// </summary>
        public void Clear()
        {
            Scan = string.Empty;
            ScanNum = 0;
            DatasetName = string.Empty;
            Sequence = string.Empty;
            ModifiedPeptide = string.Empty;
            PrefixResidue = string.Empty;
            SuffixResidue = string.Empty;
            Length = 0;
            Charge = string.Empty;
            ChargeNum = 0;
            ElutionTime = string.Empty;
            ElutionTimeAverage = string.Empty;
            PrecursorMonoMass = string.Empty;
            CalibratedObservedMass = string.Empty;
            PrecursorMZ = string.Empty;
            CalibratedObservedMZ = string.Empty;
            MH = string.Empty;
            CalculatedMonoMass = string.Empty;
            CalculatedMonoMassValue = 0;
            CalculatedMonoMassPHRP = 0;
            CalculatedMZ = string.Empty;
            MassErrorDaMSFragger = string.Empty;
            Expectation = string.Empty;
            EValue = 0;
            Hyperscore = string.Empty;
            Nextscore = string.Empty;
            PeptideProphetProbability = string.Empty;
            NumberOfTrypticTermini = 0;
            MissedCleavageCount = string.Empty;
            ProteinStart = string.Empty;
            ProteinEnd = string.Empty;
            Intensity = string.Empty;
            ModificationList = string.Empty;
            ObservedModifications = string.Empty;
            IsUnique = string.Empty;
            Protein = string.Empty;
            ProteinID = string.Empty;
            EntryName = string.Empty;
            Gene = string.Empty;
            ProteinDescription = string.Empty;
            AdditionalGenes = string.Empty;
            AdditionalProteins = string.Empty;
            NumberOfMatchedIons = string.Empty;
            TotalNumberOfIons = string.Empty;
        }

        /// <summary>
        /// Show scan, peptide, and E-value
        /// </summary>
        public override string ToString()
        {
            return string.Format("Scan {0}: {1}, E-value {2}", Scan, Sequence, Expectation);
        }
    }
}
