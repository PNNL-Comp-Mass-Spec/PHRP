using PHRPReader;
using PHRPReader.Data;
using System.Collections.Generic;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track the peptide details for a DIA-NN search result loaded from a synopsis file
    /// </summary>
    /// <remarks>
    /// See SearchResultsBaseClass for additional information
    /// </remarks>
    public class DiaNNResults : SearchResultsBaseClass
    {
        // Ignore Spelling: Averagine, Dia, NTT, Proteotypic, tryptic

        /// <summary>
        /// Dataset name
        /// </summary>
        /// <remarks>
        /// This is an abbreviated form of the full dataset name, as assigned by the analysis manager
        /// </remarks>
        public string DatasetName { get; set; }

        /// <summary>
        /// Dataset ID
        /// </summary>
        public int DatasetID { get; set; }

        /// <summary>
        /// Ion Mobility
        /// </summary>
        public string IonMobility { get; set; }

        /// <summary>
        /// Precursor ion m/z (computed from peptide mass and charge by PHRP)
        /// </summary>
        public string PrecursorMZ { get; set; }

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by PHRP
        /// </summary>
        public string CalculatedMonoMass { get; set; }

        /// <summary>
        /// Comma-separated list of static and dynamic modification names and affected residue number
        /// </summary>
        /// <remarks>
        /// <para>
        /// This is read from column Modified.Sequence in Dataset_report.tsv, then converted into MSFragger style modification names by PHRP
        /// </para>
        /// <para>
        /// Example values:
        ///   15M(15.9949)
        ///   8C(57.0215)
        ///   5M(15.9949), 8M(15.9949)
        ///   1M(15.9949), 5C(57.0215)
        ///   12M(15.9949), 2M(15.9949), 8M(15.9949)
        ///   N-term(42.0106)
        /// </para>
        /// </remarks>
        public string Modifications { get; set; }

        /// <summary>
        /// Protein group name
        /// </summary>
        public string ProteinGroup { get; set; }

        /// <summary>
        /// Proteins associated with this peptide
        /// </summary>
        /// <remarks>
        /// Corresponds to the protein names in the FASTA file
        /// </remarks>
        public List<string> ProteinIDs { get; } = new();

        /// <summary>
        /// UniProt protein names associated with this peptide
        /// </summary>
        /// <remarks>
        /// Determined by DIA-NN, typically corresponding to UniProt Name
        /// </remarks>
        public List<string> ProteinNames { get; } = new();

        /// <summary>
        /// Gene names associated with this peptide
        /// </summary>
        /// <remarks>
        /// Determined by DIA-NN
        /// </remarks>
        public List<string> GeneNames { get; } = new();

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        public string NTT { get; set; }

        /// <summary>
        /// Protein Group Quantity
        /// </summary>
        public string ProteinGroupQuantity { get; set; }

        /// <summary>
        /// Protein Group Normalized
        /// </summary>
        public string ProteinGroupNormalized { get; set; }

        /// <summary>
        /// Protein Group Max LFQ
        /// </summary>
        public string ProteinGroupMaxLFQ { get; set; }

        /// <summary>
        /// Genes Quantity
        /// </summary>
        public string GenesQuantity { get; set; }

        /// <summary>
        /// Genes Normalized
        /// </summary>
        public string GenesNormalized { get; set; }

        /// <summary>
        /// Genes Max LFQ
        /// </summary>
        public string GenesMaxLFQ { get; set; }

        /// <summary>
        /// Genes Max LFQ Unique
        /// </summary>
        public string GenesMaxLFQUnique { get; set; }

        /// <summary>
        /// QValue computed by DIA-NN
        /// </summary>
        public string QValue { get; set; }

        /// <summary>
        /// PEP (posterior error probability)
        /// </summary>
        /// <remarks>
        /// Similar to p-value
        /// Smaller values (closer to zero) are higher confidence
        /// </remarks>
        public string PEP { get; set; }

        /// <summary>
        /// Global QValue
        /// </summary>
        public string GlobalQValue { get; set; }

        /// <summary>
        /// Protein QValue
        /// </summary>
        public string ProteinQValue { get; set; }

        /// <summary>
        /// Protein Group QValue
        /// </summary>
        public string ProteinGroupQValue { get; set; }

        /// <summary>
        /// Global Protein Group QValue
        /// </summary>
        public string GlobalProteinGroupQValue { get; set; }

        /// <summary>
        /// Gene Group QValue
        /// </summary>
        public string GeneGroupQValue { get; set; }

        /// <summary>
        /// Translated QValue
        /// </summary>
        /// <remarks>
        /// Translation involves sequence propagation
        /// </remarks>
        public string TranslatedQValue;

        /// <summary>
        /// Proteotypic
        /// </summary>
        /// <remarks>
        /// 1 if the peptide is specific to a given protein or gene, otherwise 0
        /// </remarks>
        public string Proteotypic { get; set; }

        /// <summary>
        /// Precursor Quantity
        /// </summary>
        public string PrecursorQuantity { get; set; }

        /// <summary>
        /// Precursor Normalized
        /// </summary>
        public string PrecursorNormalized { get; set; }

        /// <summary>
        /// Precursor Translated
        /// </summary>
        public string PrecursorTranslated { get; set; }

        /// <summary>
        /// Translated Quality
        /// </summary>
        public string TranslatedQuality;

        /// <summary>
        /// MS1 Translated
        /// </summary>
        public string MS1Translated;

        /// <summary>
        /// Quantity Quality
        /// </summary>
        public string QuantityQuality;

        /// <summary>
        /// Elution Time (aka retention time)
        /// </summary>
        public string ElutionTime { get; set; }

        /// <summary>
        /// Elution Time Start
        /// </summary>
        public string ElutionTimeStart { get; set; }

        /// <summary>
        /// Elution Time Stop
        /// </summary>
        public string ElutionTimeStop { get; set; }

        /// <summary>
        /// Indexed Retention Time (iRT)
        /// </summary>
        public string IndexedRT { get; set; }

        /// <summary>
        /// Indexed Ion Mobility
        /// </summary>
        public string IndexedIonMobility { get; set; }

        /// <summary>
        /// Predicted Retention Time
        /// </summary>
        public string PredictedRT { get; set; }

        /// <summary>
        /// Predicted Indexed Retention Time
        /// </summary>
        public string PredictedIndexedRT { get; set; }

        /// <summary>
        /// MS1 Profile Correlation
        /// </summary>
        public string MS1ProfileCorrelation { get; set; }

        /// <summary>
        /// MS1 Area
        /// </summary>
        public string MS1Area { get; set; }

        /// <summary>
        /// Evidence (score)
        /// </summary>
        /// <remarks>Higher values are better</remarks>
        public string Evidence { get; set; }

        /// <summary>
        /// Spectrum Similarity
        /// </summary>
        public string SpectrumSimilarity { get; set; }

        /// <summary>
        /// Averagine
        /// </summary>
        public string Averagine { get; set; }

        /// <summary>
        /// Mass Evidence
        /// </summary>
        /// <remarks>Higher values are better</remarks>
        public string MassEvidence { get; set; }

        /// <summary>
        /// CScore
        /// </summary>
        public string CScore { get; set; }

        /// <summary>
        /// Decoy Evidence
        /// </summary>
        public string DecoyEvidence { get; set; }

        /// <summary>
        /// Decoy CScore
        /// </summary>
        public string DecoyCScore { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks>
        /// The base class constructor calls InitializeLocalVariables,
        /// which calls both the base class's Clear method and this class's Clear method
        /// </remarks>
        /// <param name="peptideMods">Peptide modifications</param>
        /// <param name="peptideSeqMassCalculator">Peptide sequence mass calculator</param>
        public DiaNNResults(PeptideModificationContainer peptideMods, PeptideMassCalculator peptideSeqMassCalculator)
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
            IonMobility = string.Empty;
            PrecursorMZ = string.Empty;
            CalculatedMonoMass = string.Empty;
            Modifications = string.Empty;
            ProteinGroup = string.Empty;
            ProteinIDs.Clear();
            ProteinNames.Clear();
            GeneNames.Clear();
            NTT = string.Empty;
            ProteinGroupQuantity = string.Empty;
            ProteinGroupNormalized = string.Empty;
            ProteinGroupMaxLFQ = string.Empty;
            GenesQuantity = string.Empty;
            GenesNormalized = string.Empty;
            GenesMaxLFQ = string.Empty;
            GenesMaxLFQUnique = string.Empty;
            QValue = string.Empty;
            PEP = string.Empty;
            GlobalQValue = string.Empty;
            ProteinQValue = string.Empty;
            ProteinGroupQValue = string.Empty;
            GlobalProteinGroupQValue = string.Empty;
            GeneGroupQValue = string.Empty;
            TranslatedQValue = string.Empty;
            Proteotypic = string.Empty;
            PrecursorQuantity = string.Empty;
            PrecursorNormalized = string.Empty;
            PrecursorTranslated = string.Empty;
            TranslatedQuality = string.Empty;
            MS1Translated = string.Empty;
            QuantityQuality = string.Empty;
            ElutionTime = string.Empty;
            ElutionTimeStart = string.Empty;
            ElutionTimeStop = string.Empty;
            IndexedRT = string.Empty;
            IndexedIonMobility = string.Empty;
            PredictedRT = string.Empty;
            PredictedIndexedRT = string.Empty;
            MS1ProfileCorrelation = string.Empty;
            MS1Area = string.Empty;
            Evidence = string.Empty;
            SpectrumSimilarity = string.Empty;
            Averagine = string.Empty;
            MassEvidence = string.Empty;
            CScore = string.Empty;
            DecoyEvidence = string.Empty;
            DecoyCScore = string.Empty;
        }
    }
}
