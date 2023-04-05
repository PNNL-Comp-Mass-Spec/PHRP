namespace PeptideHitResultsProcessor.SearchToolResults
{
    /// <summary>
    /// This class holds rows read from the DIA-NN report.tsv file
    /// </summary>
    /// <remarks>
    /// It is used when creating the synopsis file
    /// </remarks>
    internal class DiaNNSearchResult : ToolResultsBaseClass
    {
        // Ignore Spelling: Averagine, Hyperscore

        /// <summary>
        /// Dataset .mzML file path
        /// </summary>
        public string DatasetFile;

        /// <summary>
        /// Protein group name
        /// </summary>
        public string ProteinGroup;

        /// <summary>
        /// Protein names (from the FASTA file)
        /// </summary>
        public string ProteinIDs;

        /// <summary>
        /// Protein names, as determined by DIA-NN, typically corresponding to UniProt Name
        /// </summary>
        ///
        public string ProteinNames;

        /// <summary>
        /// Gene names associated with the peptide
        /// </summary>
        public string Genes;

        /// <summary>
        /// Protein Group Quantity
        /// </summary>
        public string ProteinGroupQuantity;

        /// <summary>
        /// Protein Group Normalized
        /// </summary>
        public string ProteinGroupNormalized;

        /// <summary>
        /// Protein Group Max LFQ
        /// </summary>
        public string ProteinGroupMaxLFQ;

        /// <summary>
        /// Genes Quantity
        /// </summary>
        public string GenesQuantity;

        /// <summary>
        /// Genes Normalized
        /// </summary>
        public string GenesNormalized;

        /// <summary>
        /// Genes Max LFQ
        /// </summary>
        public string GenesMaxLFQ;

        /// <summary>
        /// Genes Max LFQ Unique
        /// </summary>
        public string GenesMaxLFQUnique;

        // ReSharper disable CommentTypo

        /// <summary>
        /// Peptide sequence with modifications
        /// </summary>
        /// <remarks>
        /// Examples:
        ///   AAALEAM(UniMod:35)K
        ///   AAEAHVDAHYYEQNEQPTGTC(UniMod:4)AAC(UniMod:4)ITGDNR
        /// </remarks>
        public string ModifiedSequence;

        // ReSharper restore CommentTypo

        // Stripped.Sequence is tracked by Sequence in the base class
        // public string StrippedSequence;

        /// <summary>
        /// Precursor sequence (unused by PHRP)
        /// </summary>
        public string PrecursorId;

        /// <summary>
        /// QValue computed by DIA-NN
        /// </summary>
        public string QValueDiaNN;

        /// <summary>
        /// PEP (posterior error probability)
        /// </summary>
        public string PEP;

        /// <summary>
        /// Global QValue
        /// </summary>
        public string GlobalQValue;

        /// <summary>
        /// Protein QValue
        /// </summary>
        public string ProteinQValue;

        /// <summary>
        /// Protein Group QValue
        /// </summary>
        public string ProteinGroupQValue;

        /// <summary>
        /// Global Protein Group QValue
        /// </summary>
        public string GlobalProteinGroupQValue;

        /// <summary>
        /// Gene Group QValue
        /// </summary>
        public string GeneGroupQValue;

        /// <summary>
        /// Translated QValue
        /// </summary>
        public string TranslatedQValue;

        /// <summary>
        /// 1 if the peptide is specific to a given protein or gene, otherwise 0
        /// </summary>
        public string Proteotypic;

        /// <summary>
        /// Precursor Quantity
        /// </summary>
        public string PrecursorQuantity;

        /// <summary>
        /// Precursor Normalized
        /// </summary>
        public string PrecursorNormalized;

        /// <summary>
        /// Precursor Translated
        /// </summary>
        public string PrecursorTranslated;

        /// <summary>
        /// Translated Quality (unused by PHRP)
        /// </summary>
        public string TranslatedQuality;

        /// <summary>
        /// MS1 Translated (unused by PHRP)
        /// </summary>
        public string Ms1Translated;

        /// <summary>
        /// Quantity Quality (unused by PHRP)
        /// </summary>
        public string QuantityQuality;

        // Elution time (aka retention time) is tracked by ElutionTime in the base class
        // public string RT;

        /// <summary>
        /// Elution Time Start
        /// </summary>
        public string RTStart;

        /// <summary>
        /// Elution Time Stop
        /// </summary>
        public string RTStop;

        /// <summary>
        /// Indexed retention time (iRT)
        /// </summary>
        public string IndexedRT;

        /// <summary>
        /// Predicted retention time (unused by PHRP)
        /// </summary>
        public string PredictedRT;

        /// <summary>
        /// Predicted indexed retention time (unused by PHRP)
        /// </summary>
        public string PredictedIRT;

        /// <summary>
        /// First protein description (unused by PHRP)
        /// </summary>
        public string FirstProteinDescription;

        /// <summary>
        /// Spectral Library QValue
        /// </summary>
        public string LibQValue;

        /// <summary>
        /// Spectral Library Protein Group QValue
        /// </summary>
        public string LibProteinGroupQValue;

        /// <summary>
        /// MS1 Profile Correlation
        /// </summary>
        public string Ms1ProfileCorr;

        /// <summary>
        /// MS1 Area
        /// </summary>
        public string Ms1Area;

        /// <summary>
        /// Evidence (score)
        /// </summary>
        /// <remarks>Higher values are better</remarks>
        public string Evidence;

        /// <summary>
        /// Spectrum Similarity
        /// </summary>
        public string SpectrumSimilarity;

        /// <summary>
        /// Averagine
        /// </summary>
        public string Averagine;

        /// <summary>
        /// Mass Evidence
        /// </summary>
        /// <remarks>Higher values are better</remarks>
        public string MassEvidence;

        /// <summary>
        /// CScore
        /// </summary>
        public string CScore;

        /// <summary>
        /// Decoy Evidence
        /// </summary>
        public string DecoyEvidence;

        /// <summary>
        /// Decoy CScore
        /// </summary>
        public string DecoyCScore;

        /// <summary>
        /// Fragmentation Ion Abundances
        /// </summary>
        public string FragmentQuantRaw;

        /// <summary>
        /// Corrected Fragmentation Ion Abundances
        /// </summary>
        public string FragmentQuantCorrected;

        /// <summary>
        /// Fragmentation ion correlations
        /// </summary>
        public string FragmentCorrelations;

        /// <summary>
        /// Ion Mobility
        /// </summary>
        public string IM;

        /// <summary>
        /// Indexed Ion Mobility
        /// </summary>
        public string IndexedIM;

        /// <summary>
        /// Predicted Ion Mobility
        /// </summary>
        public string PredictedIM;

        /// <summary>
        /// Predicted Indexed Ion Mobility
        /// </summary>
        public string PredictedIndexedIM;

        /// <summary>
        /// Reset stored values to empty strings and zeros
        /// </summary>
        public void Clear()
        {
            Scan = string.Empty;
            ScanNum = 0;
            DatasetFile = string.Empty;
            DatasetName = string.Empty;

            ProteinGroup = string.Empty;
            ProteinIDs = string.Empty;
            ProteinNames = string.Empty;
            Genes = string.Empty;
            ProteinGroupQuantity = string.Empty;
            ProteinGroupNormalized = string.Empty;
            ProteinGroupMaxLFQ = string.Empty;
            GenesQuantity = string.Empty;
            GenesNormalized = string.Empty;
            GenesMaxLFQ = string.Empty;
            GenesMaxLFQUnique = string.Empty;
            ModifiedSequence = string.Empty;
            Sequence = string.Empty;
            PrecursorId = string.Empty;
            Charge = string.Empty;
            ChargeNum = 0;
            QValueDiaNN = string.Empty;
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
            Ms1Translated = string.Empty;
            QuantityQuality = string.Empty;
            ElutionTime = string.Empty;
            RTStart = string.Empty;
            RTStop = string.Empty;
            IndexedRT = string.Empty;
            PredictedRT = string.Empty;
            PredictedIRT = string.Empty;
            FirstProteinDescription = string.Empty;
            LibQValue = string.Empty;
            LibProteinGroupQValue = string.Empty;
            Ms1ProfileCorr = string.Empty;
            Ms1Area = string.Empty;
            Evidence = string.Empty;
            SpectrumSimilarity = string.Empty;
            Averagine = string.Empty;
            MassEvidence = string.Empty;
            CScore = string.Empty;
            DecoyEvidence = string.Empty;
            DecoyCScore = string.Empty;
            FragmentQuantRaw = string.Empty;
            FragmentQuantCorrected = string.Empty;
            FragmentCorrelations = string.Empty;
            IM = string.Empty;
            IndexedIM = string.Empty;
            PredictedIM = string.Empty;
            PredictedIndexedIM = string.Empty;

            PrefixResidue = string.Empty;
            SuffixResidue = string.Empty;
            Length = 0;
            MH = string.Empty;
            CalculatedMonoMass = string.Empty;
            CalculatedMonoMassValue = 0;
            CalculatedMonoMassPHRP = 0;
            NumberOfTrypticTermini = 0;
            MissedCleavageCount = string.Empty;
            ProteinStart = string.Empty;
            ProteinEnd = string.Empty;
            Intensity = string.Empty;
            ModificationList = string.Empty;
            Protein = string.Empty;
        }

        /// <summary>
        /// Show scan, peptide, and Q Value
        /// </summary>
        public override string ToString()
        {
            return string.Format("Scan {0}: {1}, E-value {2}", Scan, Sequence, QValueDiaNN);
        }
    }
}
