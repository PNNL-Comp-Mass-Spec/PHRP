// ReSharper disable UnusedMember.Global

namespace PHRPReader.Data
{
    // Ignore Spelling: Da, Daltons, Hyperscore, MaxQuant, proteoform, novo, tda, tryptic

    /// <summary>
    /// These columns correspond to the tab-delimited _ScanStatsEx.txt file created by MASIC or MSFileInfoScanner
    /// </summary>
    public enum ExtendedScanStatsFileColumns
    {
        /// <summary>
        /// Dataset Id
        /// </summary>
        DatasetId = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        ScanNumber = 1,

        /// <summary>
        /// Ion injection time
        /// </summary>
        IonInjectionTime = 2,

        /// <summary>
        /// Scan event
        /// </summary>
        ScanEvent = 3,

        /// <summary>
        /// Master index
        /// </summary>
        MasterIndex = 4,

        /// <summary>
        /// Elapsed scan time
        /// </summary>
        ElapsedScanTime = 5,

        /// <summary>
        /// Charge state
        /// </summary>
        ChargeState = 6,

        /// <summary>
        /// Monoisotopic m/z
        /// </summary>
        MonoisotopicMZ = 7,

        /// <summary>
        /// MS2 isolation width
        /// </summary>
        MS2IsolationWidth = 8,

        /// <summary>
        /// FT analyzer settings
        /// </summary>
        FTAnalyzerSettings = 9,

        /// <summary>
        /// FT analyzer message
        /// </summary>
        FTAnalyzerMessage = 10,

        /// <summary>
        /// FTResolution
        /// </summary>
        FTResolution = 11,

        /// <summary>
        /// Conversion parameter B
        /// </summary>
        ConversionParameterB = 12,

        /// <summary>
        /// Conversion parameter C
        /// </summary>
        ConversionParameterC = 13,

        /// <summary>
        /// Conversion parameter D
        /// </summary>
        ConversionParameterD = 14,

        /// <summary>
        /// Conversion parameter E
        /// </summary>
        ConversionParameterE = 15,

        /// <summary>
        /// Collision mode
        /// </summary>
        CollisionMode = 16,

        /// <summary>
        /// Scan filter text
        /// </summary>
        ScanFilterText = 17,

        /// <summary>
        /// Source voltage
        /// </summary>
        SourceVoltage = 18,

        /// <summary>
        /// Source current
        /// </summary>
        Source_Current = 19
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by InSpecTResultsProcessor
    /// </summary>
    public enum InspectSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Peptide
        /// </summary>
        Peptide = 2,

        /// <summary>
        /// Protein
        /// </summary>
        Protein = 3,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 4,

        /// <summary>
        /// MQScore
        /// </summary>
        MQScore = 5,

        /// <summary>
        /// Length
        /// </summary>
        Length = 6,

        /// <summary>
        /// Total PRM score
        /// </summary>
        TotalPRMScore = 7,

        /// <summary>
        /// Median PRM score
        /// </summary>
        MedianPRMScore = 8,

        /// <summary>
        /// Fraction Y
        /// </summary>
        FractionY = 9,

        /// <summary>
        /// Fraction B
        /// </summary>
        FractionB = 10,

        /// <summary>
        /// Intensity
        /// </summary>
        Intensity = 11,

        /// <summary>
        /// NTT
        /// </summary>
        NTT = 12,

        /// <summary>
        /// PValue
        /// </summary>
        PValue = 13,

        /// <summary>
        /// FScore
        /// </summary>
        FScore = 14,

        /// <summary>
        /// Delta score
        /// </summary>
        DeltaScore = 15,

        /// <summary>
        /// Delta score other
        /// </summary>
        DeltaScoreOther = 16,

        /// <summary>
        /// Delta norm MQ Score
        /// </summary>
        /// <remarks>
        /// <para>
        /// Computed as Abs((MQScore(n) - MQScore(n+1)) / MQScore(n)); storing 0 for the lowest scoring result in each set.
        /// If MQScore(n) is 0, stores 0.
        /// </para>
        /// <para>
        /// This value is not usable when MQScore(n) is 0 or negative, and should generally not be used when MQScore(n) is less than 0.5
        /// </para>
        /// </remarks>
        DeltaNormMQScore = 17,

        /// <summary>
        /// Delta norm total PRM score
        /// </summary>
        /// <remarks>
        /// <para>
        /// Computed as Abs((TotalPRMScore(n) - TotalPRMScore(n+1)) / TotalPRMScore(n)); storing 0 for the lowest scoring result in each set.
        /// If TotalPRMScore(n) is 0, stores 0.
        /// </para>
        /// <para>
        /// This value is not usable when TotalPRMScore(n) is 0 or negative, and should generally not be used when TotalPRMScore(n) is less than 0.5
        /// </para>
        /// </remarks>
        DeltaNormTotalPRMScore = 18,

        /// <summary>
        /// Rank total PRM score
        /// </summary>
        /// <remarks>
        /// Rank 1 means highest TotalPRMScore, 2 means next lower score, etc. (ties get the same rank)
        /// </remarks>
        RankTotalPRMScore = 19,

        /// <summary>
        /// Rank F score
        /// </summary>
        /// <remarks>
        /// Rank 1 means highest FScore, 2 means next lower, etc. (ties get the same rank)
        /// </remarks>
        RankFScore = 20,

        /// <summary>
        /// Theoretical monoisotopic peptide mass (computed by PHRP)
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 21,

        /// <summary>
        /// Record number
        /// </summary>
        RecordNumber = 22,

        /// <summary>
        /// DB file position
        /// </summary>
        DBFilePos = 23,

        /// <summary>
        /// Spec file position
        /// </summary>
        SpecFilePos = 24,

        /// <summary>
        /// Precursor m/z
        /// </summary>
        PrecursorMZ = 25,

        /// <summary>
        /// Precursor error
        /// </summary>
        PrecursorError = 26,

        /// <summary>
        /// Precursor error, in ppm
        /// </summary>
        DelM_PPM = 27
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MaxQuantResultsProcessor
    /// </summary>
    public enum MaxQuantSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Dataset name (typically the .raw file name, without the extension)
        /// </summary>
        Dataset = 1,

        /// <summary>
        /// Dataset ID
        /// </summary>
        DatasetID = 2,

        /// <summary>
        /// MS/MS scan number in the dataset
        /// </summary>
        Scan = 3,

        /// <summary>
        /// Type of fragmentation used to create the MS/MS spectrum
        /// </summary>
        /// <remarks>
        /// Types:
        ///   CID - Collision Induced Dissociation
        ///   HCD - High energy Collision induced Dissociation
        ///   ETD - Electron Transfer Dissociation
        /// </remarks>
        FragMethod = 4,

        /// <summary>
        /// Index of the spectrum in the dataset (1-based, consecutive integer)
        /// </summary>
        SpecIndex = 5,

        /// <summary>
        /// Charge state of the precursor ion
        /// </summary>
        Charge = 6,

        /// <summary>
        /// Precursor ion m/z (observed value, read from a _PrecursorInfo.txt file by PHRP)
        /// </summary>
        PrecursorMZ = 7,

        /// <summary>
        /// Precursor ion m/z (theoretical value, not observed value)
        /// </summary>
        /// <remarks>
        /// This is the theoretical m/z of the first isotope of the isotopic distribution of the parent ion
        /// </remarks>
        PrecursorMZ_MaxQuant = 8,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence
        /// </summary>
        DelM = 9,

        /// <summary>
        /// Mass error, in ppm
        /// </summary>
        DelM_PPM = 10,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence,
        /// as computed by MaxQuant
        /// </summary>
        DelM_MaxQuant = 11,

        /// <summary>
        /// Mass error, in ppm, as computed by MaxQuant
        /// </summary>
        DelM_PPM_MaxQuant = 12,

        /// <summary>
        /// Monoisotopic (M+H)+ value, computed from PrecursorMZ and Charge by PHRP
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 13,

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MaxQuant
        /// </summary>
        Mass = 14,

        /// <summary>
        /// Peptide sequence, including modification symbols
        /// </summary>
        Peptide = 15,

        /// <summary>
        /// Comma separated list of dynamic modifications in the peptide
        /// </summary>
        DynamicModifications = 16,

        /// <summary>
        /// Protein names
        /// </summary>
        /// <remarks>
        /// <para>
        /// Semicolon separated list
        /// </para>
        /// <para>
        /// If the peptide is from a reverse-hit protein, the Proteins column will be empty in the _syn.txt file,
        /// and the decoy protein name will be in the LeadingRazorProtein column.
        /// When the MaxQuantSynFileReader class reads the _syn.txt file, it adds the LeadingRazorProtein name
        /// to the list of protein names, if the Proteins column is empty.
        /// </para>
        /// </remarks>
        Proteins = 17,

        /// <summary>
        /// Name of the best scoring protein this peptide is associated with
        /// </summary>
        /// <remarks>
        /// Typically there is only one protein name here
        /// However, in cases of a tied score, will be a semicolon separated list
        /// </remarks>
        LeadingRazorProtein = 18,

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        NTT = 19,

        /// <summary>
        /// Posterior error probability
        /// </summary>
        /// <remarks>
        /// Similar to p-value
        /// Smaller values (closer to zero) are higher confidence
        /// </remarks>
        PEP = 20,

        /// <summary>
        /// Andromeda score for the best MS/MS spectrum with this peptide
        /// </summary>
        /// <remarks>
        /// Larger values are better
        /// </remarks>
        Score = 21,

        /// <summary>
        /// Score difference to the second best identified peptide with a different amino acid sequence
        /// </summary>
        DeltaScore = 22,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means highest Andromeda score, 2 means next lower score, etc. (ties get the same rank)
        /// </remarks>
        RankScore = 23,

        /// <summary>
        /// Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
        /// </summary>
        /// <remarks>
        /// All peptides with the same sequence will have the same total peptide intensity value for a given dataset
        /// </remarks>
        TotalPeptideIntensity = 24,

        /// <summary>
        /// Mass Analyzer of the instrument
        /// </summary>
        MassAnalyzer = 25,

        /// <summary>
        /// Type of precursor ion as identified by MaxQuant
        /// </summary>
        PrecursorType = 26,

        /// <summary>
        /// Elution time of the MS/MS spectrum
        /// </summary>
        ElutionTime = 27,

        /// <summary>
        /// Scan number where the precursor ion was observed
        /// </summary>
        PrecursorScan = 28,

        /// <summary>
        /// Intensity of the precursor ion in the scan that it was observed
        /// </summary>
        PrecursorIntensity = 29,

        /// <summary>
        /// Number of peaks (MS/MS ions) matching to the predicted fragmentation spectrum
        /// </summary>
        NumberOfMatches = 30,

        /// <summary>
        /// Fraction of intensity in the MS/MS spectrum that is annotated
        /// </summary>
        IntensityCoverage = 31,

        /// <summary>
        /// Number of missed enzymatic cleavages
        /// </summary>
        MissedCleavages = 32,

        /// <summary>
        /// Unique (consecutive) identifier for each row in msms.txt
        /// </summary>
        MsMsID = 33,

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
        ProteinGroupIDs = 34,

        /// <summary>
        /// The identifier of the non-redundant peptide sequence
        /// </summary>
        /// <remarks>
        /// Corresponds to the id column in the peptides.txt file
        /// </remarks>
        PeptideID = 35,

        /// <summary>
        /// Identifier referencing a row in the modificationSpecificPeptides.txt file
        /// </summary>
        ModPeptideID = 36,

        /// <summary>
        /// Identifier referencing a row in the evidence.txt file
        /// </summary>
        EvidenceID = 37,

        /// <summary>
        /// Q-Value (computed by PHRP)
        /// </summary>
        QValue = 38
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MODPlusResultsProcessor
    /// </summary>
    public enum MODPlusSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Spectrum index
        /// </summary>
        Spectrum_Index = 2,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 3,

        /// <summary>
        /// Precursor m/z
        /// </summary>
        PrecursorMZ = 4,

        /// <summary>
        /// Precursor error, in Da
        /// </summary>
        DelM = 5,

        /// <summary>
        /// Precursor error, in ppm; corrected for isotope selection errors
        /// </summary>
        DelM_PPM = 6,

        /// <summary>
        /// Theoretical monoisotopic peptide MH (computed by PHRP)
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 7,

        /// <summary>
        /// This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
        /// </summary>
        Peptide = 8,

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        NTT = 9,

        /// <summary>
        /// Modification annotation
        /// </summary>
        ModificationAnnotation = 10,

        /// <summary>
        /// Protein
        /// </summary>
        Protein = 11,

        /// <summary>
        /// Peptide position
        /// </summary>
        Peptide_Position = 12,

        /// <summary>
        /// Score
        /// </summary>
        Score = 13,

        /// <summary>
        /// Probability
        /// </summary>
        Probability = 14,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means highest MODPlus score, 2 means next lower score, etc. (ties get the same rank)
        /// </remarks>
        Rank_Score = 15,

        /// <summary>
        /// Q-Value (FDR)
        /// </summary>
        QValue = 16
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MSAlignResultsProcessor
    /// </summary>
    public enum MSAlignSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Proteoform spectrum match ID
        /// </summary>
        Prsm_ID = 2,

        /// <summary>
        /// Spectrum ID
        /// </summary>
        Spectrum_ID = 3,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 4,

        /// <summary>
        /// Precursor m/z
        /// </summary>
        PrecursorMZ = 5,

        /// <summary>
        /// Precursor error, in Da
        /// </summary>
        DelM = 6,

        /// <summary>
        /// Precursor error, in ppm; corrected for isotope selection errors
        /// </summary>
        DelMPPM = 7,

        /// <summary>
        /// Theoretical monoisotopic peptide MH (computed by PHRP)
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 8,

        /// <summary>
        /// This is the sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
        /// </summary>
        Peptide = 9,

        /// <summary>
        /// Protein name
        /// </summary>
        Protein = 10,

        /// <summary>
        /// Protein mass
        /// </summary>
        Protein_Mass = 11,

        /// <summary>
        /// Unexpected mod count
        /// </summary>
        Unexpected_Mod_Count = 12,

        /// <summary>
        /// Peak count
        /// </summary>
        Peak_Count = 13,

        /// <summary>
        /// Matched peak count
        /// </summary>
        Matched_Peak_Count = 14,

        /// <summary>
        /// Matched fragment ion count
        /// </summary>
        Matched_Fragment_Ion_Count = 15,

        /// <summary>
        /// P-value
        /// </summary>
        PValue = 16,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means lowest (best) p-value, 2 means next higher p-value, etc. (ties get the same rank)
        /// </remarks>
        Rank_PValue = 17,

        /// <summary>
        /// E-value
        /// </summary>
        EValue = 18,

        /// <summary>
        /// FDR
        /// </summary>
        FDR = 19,

        /// <summary>
        /// Species ID
        /// </summary>
        Species_ID = 20,

        /// <summary>
        /// Fragmentation method
        /// </summary>
        FragMethod = 21
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _MSGF.txt file created by the analysis manager
    /// </summary>
    public enum MSGFFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 2,

        /// <summary>
        /// Protein
        /// </summary>
        Protein = 3,

        /// <summary>
        /// Peptide
        /// </summary>
        Peptide = 4,

        /// <summary>
        /// Spectral probability
        /// </summary>
        SpecProb = 5,

        /// <summary>
        /// Notes
        /// </summary>
        Notes = 6
    }

    /// <summary>
    /// DIA-NN synopsis file columns
    /// </summary>
    public enum DiaNNSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Dataset name (abbreviated dataset name, as assigned by the analysis manager)
        /// </summary>
        Dataset = 1,

        /// <summary>
        /// Dataset ID
        /// </summary>
        DatasetID = 2,

        /// <summary>
        /// MS/MS scan number in the dataset
        /// </summary>
        Scan = 3,

        /// <summary>
        /// Ion Mobility
        /// </summary>
        IonMobility = 4,

        /// <summary>
        /// Charge state of the precursor ion
        /// </summary>
        Charge = 5,

        /// <summary>
        /// Precursor ion m/z (computed from peptide mass and charge by PHRP)
        /// </summary>
        PrecursorMZ = 6,

        /// <summary>
        /// Monoisotopic (M+H)+ value, computed from peptide mass by PHRP
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 9,

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by PHRP
        /// </summary>
        Mass = 10,

        /// <summary>
        /// Peptide sequence, without any modifications
        /// </summary>
        Peptide = 11,

        /// <summary>
        /// ost-translational modifications contained within the identified peptide
        /// </summary>
        /// <remarks>
        /// <para>
        /// For DIA-NN, this is read from column Modified.Sequence in report.tsv, then converted into MSFragger style modification names
        /// </para>
        /// <para>
        /// DIA-NN Examples:
        ///   15M(15.9949)
        ///   1M(15.9949), 5C(57.0215)
        ///   N-term(42.0106)
        /// </para>
        /// </remarks>
        Modifications = 12,

        /// <summary>
        /// Protein group name
        /// </summary>
        ProteinGroup = 13,

        /// <summary>
        /// Protein names (from the FASTA file)
        /// </summary>
        ProteinIDs = 14,

        /// <summary>
        /// Protein names, as determined by DIA-NN, typically corresponding to UniProt Name
        /// </summary>
        ProteinNames = 15,

        /// <summary>
        /// Gene names associated with the peptide
        /// </summary>
        GeneNames = 16,

        /// <summary>
        /// Number of tryptic terminii
        /// </summary>
        NTT = 17,

        /// <summary>
        /// Protein Group Quantity
        /// </summary>
        ProteinGroupQuantity = 18,

        /// <summary>
        /// Protein Group Normalized
        /// </summary>
        ProteinGroupNormalized = 19,

        /// <summary>
        /// Protein Group Max LFQ
        /// </summary>
        ProteinGroupMaxLFQ = 20,

        /// <summary>
        /// Genes Quantity
        /// </summary>
        GenesQuantity = 21,

         /// <summary>
        /// Genes Normalized
        /// </summary>
        GenesNormalized = 22,

        /// <summary>
        /// Genes Max LFQ
        /// </summary>
        GenesMaxLFQ = 23,

        /// <summary>
        /// Genes Max LFQ Unique
        /// </summary>
        GenesMaxLFQUnique = 24,

        /// <summary>
        /// QValue computed by DIA-NN
        /// </summary>
        QValue = 25,
        /// <summary>
        /// PEP (posterior error probability)
        /// </summary>
        PEP = 26,
        /// <summary>
        /// Global QValue
        /// </summary>
        GlobalQValue = 27,

        /// <summary>
        /// Protein QValue
        /// </summary>
        ProteinQValue = 28,

        /// <summary>
        /// Protein Group QValue
        /// </summary>
        ProteinGroupQValue = 29,

        /// <summary>
        /// Global Protein Group QValue
        /// </summary>
        GlobalProteinGroupQValue = 30,

        /// <summary>
        /// Gene Group QValue
        /// </summary>
        GeneGroupQValue = 31,

        /// <summary>
        /// Precursor Quantity
        /// </summary>
        PrecursorQuantity = 32,

        /// <summary>
        /// Precursor Normalized
        /// </summary>
        PrecursorNormalized = 33,

        /// <summary>
        /// Precursor Translated
        /// </summary>
        PrecursorTranslated = 34,

        /// <summary>
        /// Elution time
        /// </summary>
        ElutionTime = 35,

        /// <summary>
        /// Elution Time Start
        /// </summary>
        ElutionTimeStart = 36,

        /// <summary>
        /// Elution Time Stop
        /// </summary>
        ElutionTimeStop = 37,
        /// <summary>
        /// MS1 Profile Correlation
        /// </summary>
        MS1ProfileCorrelation = 38,

        /// <summary>
        /// MS1 Area
        /// </summary>
        MS1Area = 39,

        /// <summary>
        /// Evidence (score)
        /// </summary>
        /// <remarks>Higher values are better</remarks>
        Evidence = 40,

        /// <summary>
        /// Spectrum Similarity
        /// </summary>
        SpectrumSimilarity = 41,

        /// <summary>
        /// Averagine
        /// </summary>
        Averagine = 42,

        /// <summary>
        /// Mass Evidence
        /// </summary>
        /// <remarks>Higher values are better</remarks>
        MassEvidence = 43,

        /// <summary>
        /// CScore
        /// </summary>
        CScore = 44,

        /// <summary>
        /// Decoy Evidence
        /// </summary>
        DecoyEvidence = 45,

        /// <summary>
        /// Decoy CScore
        /// </summary>
        DecoyCScore = 46
    }

    /// <summary>
    /// MSFragger synopsis file columns
    /// </summary>
    public enum MSFraggerSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Dataset name (typically the .raw file name, without the extension)
        /// </summary>
        Dataset = 1,

        /// <summary>
        /// Dataset ID
        /// </summary>
        DatasetID = 2,

        /// <summary>
        /// MS/MS scan number in the dataset
        /// </summary>
        Scan = 3,

        /// <summary>
        /// Charge state of the precursor ion
        /// </summary>
        Charge = 4,

        /// <summary>
        /// Precursor ion m/z (observed value)
        /// </summary>
        PrecursorMZ = 5,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence
        /// </summary>
        DelM = 6,

        /// <summary>
        /// Mass Error, in ppm
        /// </summary>
        DelM_PPM = 7,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence,
        /// as computed by MSFragger
        /// </summary>
        DelM_MSFragger = 8,

        /// <summary>
        /// Monoisotopic (M+H)+ value, computed from PrecursorMZ and Charge by PHRP
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 9,

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MSFragger
        /// </summary>
        Mass = 10,

        /// <summary>
        /// Peptide sequence, including modification symbols
        /// </summary>
        Peptide = 11,

        /// <summary>
        /// Comma separated list of static and dynamic modifications in the peptide
        /// </summary>
        Modifications = 12,

        /// <summary>
        /// Primary protein name
        /// </summary>
        Protein = 13,

        /// <summary>
        /// Additional protein names (comma separated list)
        /// </summary>
        AdditionalProteins = 14,

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        NTT = 15,

        /// <summary>
        /// Expectation value (E-value)
        /// </summary>
        EValue = 16,

        /// <summary>
        /// Ranked E-value
        /// </summary>
        /// <remarks>
        /// Rank 1 means lowest (best) E-value, 2 means next higher E-value, etc. (ties get the same rank)
        /// </remarks>
        RankEValue = 17,

        /// <summary>
        /// Hyperscore of the top match
        /// </summary>
        /// <remarks>
        /// Higher scores are better
        /// </remarks>
        Hyperscore = 18,

        /// <summary>
        /// Hyperscore of the next best match
        /// </summary>
        /// <remarks>
        /// Larger values are better
        /// </remarks>
        Nextscore = 19,

        /// <summary>
        /// Peptide prophet probability
        /// </summary>
        /// <remarks>
        /// Values closer to 1 are better
        /// </remarks>
        PeptideProphetProbability = 20,

        /// <summary>
        /// Elution time of the MS/MS spectrum
        /// </summary>
        ElutionTime = 21,

        /// <summary>
        /// Average elution time of MS/MS spectra (only applicable for a multi-dataset based search)
        /// </summary>
        ElutionTimeAverage = 22,

        /// <summary>
        /// Number of missed enzymatic cleavages
        /// </summary>
        MissedCleavages = 23,

        /// <summary>
        /// Number of matched ions in the fragmentation spectrum
        /// </summary>
        /// <remarks>
        /// Read from the Dataset.tsv file since not in the _psm.tsv file
        /// </remarks>
        NumberOfMatchedIons = 24,

        /// <summary>
        /// Total number of ions in the fragmentation spectrum
        /// </summary>
        /// <remarks>
        /// Read from the Dataset.tsv file since not in the _psm.tsv file
        /// </remarks>
        TotalNumberOfIons = 25,

        /// <summary>
        /// Q-Value (computed by PHRP)
        /// </summary>
        QValue = 26
    }

    /// <summary>
    /// These columns correspond Synopsis files created for legacy tool MSGFDB
    /// </summary>
    public enum MSGFDBSynFileColumns
    {
        /// <summary>
        /// Spectral probability
        /// </summary>
        SpecProb = 0,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means lowest (best) SpecProb, 2 means next higher SpecProb, etc. (ties get the same rank)
        /// </remarks>
        RankSpecProb = 1,

        /// <summary>
        /// PValue
        /// </summary>
        PValue = 2,

        /// <summary>
        /// FDR
        /// </summary>
        FDR = 3,

        /// <summary>
        /// Peptide level FDR
        /// </summary>
        PepFDR = 4
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MSGFPlusResultsProcessor
    /// </summary>
    public enum MSGFPlusSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Fragmentation method
        /// </summary>
        FragMethod = 2,

        /// <summary>
        /// Spectrum index
        /// </summary>
        SpecIndex = 3,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 4,

        /// <summary>
        /// Precursor ion m/z
        /// </summary>
        PrecursorMZ = 5,

        /// <summary>
        /// Precursor error, in Daltons
        /// </summary>
        /// <remarks>
        /// If the search used a tolerance less than 0.5 Da or less than 500 ppm,
        /// this value is computed from the DelMPPM value
        /// </remarks>
        DelM = 6,

        /// <summary>
        /// Precursor error, in ppm; corrected for isotope selection errors
        /// </summary>
        DelMPPM = 7,

        /// <summary>
        /// Theoretical monoisotopic peptide mass (computed by PHRP)
        /// </summary>
        MH = 8,

        /// <summary>
        /// Peptide sequence with prefix and suffix residues, plus also with modification symbols
        /// </summary>
        Peptide = 9,

        /// <summary>
        /// Protein Name (no description)
        /// </summary>
        Protein = 10,

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        NTT = 11,

        /// <summary>
        /// De-novo score
        /// </summary>
        DeNovoScore = 12,

        /// <summary>
        /// MSGF Score
        /// </summary>
        MSGFScore = 13,

        /// <summary>
        /// Spec E-Value
        /// </summary>
        SpecEValue = 14,

        // <summary>
        // Obsolete name for Spec E-Value (MS-GF+) or SpecProb (MSGFDB)
        // </summary>
        // [Obsolete("Use .SpecEValue")]
        // SpecProb_EValue = SpecEValue,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means lowest (best) Spec E-value, 2 means next higher E-value, etc. (ties get the same rank)
        /// </remarks>
        RankSpecEValue = 15,

        // <summary>
        // Obsolete name for Rank Spec E-Value (MS-GF+) or Rank SpecProb (MSGFDB)
        // </summary>
        // [Obsolete("Use .RankSpecEValue")]
        // RankSpecProb = RankSpecEValue,

        /// <summary>
        /// E-Value
        /// </summary>
        EValue = 16,

        // <summary>
        // Obsolete name for E-Value (MS-GF+) or P-value (MSGFDB)
        // </summary>
        // [Obsolete("Use .EValue")]
        // PValue_EValue = EValue,

        /// <summary>
        /// Q-Value
        /// </summary>
        /// <remarks>
        /// Only present if searched using -tda 1
        /// </remarks>
        QValue = 17,

        // <summary>
        // Obsolete name for Q-Value (MS-GF+) or FDR (MSGFDB)
        // </summary>
        // [Obsolete("Use .QValue")]
        // FDR_QValue = QValue,

        /// <summary>
        /// Peptide Q-Value
        /// </summary>
        /// <remarks>
        /// Only present if searched using -tda 1
        /// </remarks>
        PepQValue = 18,

        // <summary>
        // Obsolete name for PepQValue (MS-GF+) or PepFDR (MSGFDB)
        // </summary>
        // [Obsolete("Use .PepQValue")]
        // PepFDR_PepQValue = PepQValue,

        /// <summary>
        /// EFDR
        /// </summary>
        /// <remarks>Only present if did not search using -tda 1</remarks>
        EFDR = 19,

        // ReSharper disable UnusedMember.Global

        /// <summary>
        /// IMS scan
        /// </summary>
        /// <remarks>Only present for MSGFDB_IMS results</remarks>
        IMSScan = 20,

        /// <summary>
        /// IMS drift time
        /// </summary>
        /// <remarks>Only present for MSGFDB_IMS results</remarks>
        IMSDriftTime = 21,

        // ReSharper restore UnusedMember.Global

        /// <summary>
        /// Isotope error
        /// </summary>
        IsotopeError = 22
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MSPathFinderResultsProcessor
    /// </summary>
    public enum MSPathFinderSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 2,

        /// <summary>
        /// Most abundant isotope m/z
        /// </summary>
        MostAbundantIsotopeMz = 3,

        /// <summary>
        /// Mass
        /// </summary>
        Mass = 4,

        /// <summary>
        /// Sequence
        /// </summary>
        /// <remarks>PrefixLetter.Sequence.SuffixLetter</remarks>
        Sequence = 5,

        /// <summary>
        /// Modifications
        /// </summary>
        Modifications = 6,

        /// <summary>
        /// Composition
        /// </summary>
        Composition = 7,

        /// <summary>
        /// Protein
        /// </summary>
        Protein = 8,

        /// <summary>
        /// Protein description
        /// </summary>
        ProteinDesc = 9,

        /// <summary>
        /// Protein length
        /// </summary>
        ProteinLength = 10,

        /// <summary>
        /// Protein residue number where the sequence starts
        /// </summary>
        ResidueStart = 11,

        /// <summary>
        /// Protein residue number where the sequence ends
        /// </summary>
        ResidueEnd = 12,

        /// <summary>
        /// Count of matched fragments
        /// </summary>
        MatchedFragments = 13,

        /// <summary>
        /// Spec E-Value
        /// </summary>
        /// <remarks>Column added 2015-08-25</remarks>
        SpecEValue = 14,

        /// <summary>
        /// E-Value
        /// </summary>
        /// <remarks>Column added 2015-08-25</remarks>
        EValue = 15,

        /// <summary>
        /// Q-Value
        /// </summary>
        QValue = 16,

        /// <summary>
        /// Peptide Q-Value
        /// </summary>
        PepQValue = 17
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MODaResultsProcessor
    /// </summary>
    public enum MODaSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Spectrum index
        /// </summary>
        Spectrum_Index = 2,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 3,

        /// <summary>
        /// Precursor m/z
        /// </summary>
        PrecursorMZ = 4,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence
        /// </summary>
        DelM = 5,

        /// <summary>
        /// Precursor error, in ppm; corrected for isotope selection errors
        /// </summary>
        DelM_PPM = 6,

        /// <summary>
        /// Theoretical monoisotopic peptide MH (computed by PHRP)
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 7,

        /// <summary>
        /// This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
        /// </summary>
        Peptide = 8,

        /// <summary>
        /// Protein Name
        /// </summary>
        Protein = 9,

        /// <summary>
        /// Score
        /// </summary>
        Score = 10,

        /// <summary>
        /// Probability
        /// </summary>
        Probability = 11,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means highest probability, 2 means next lower probability, etc. (ties get the same rank)
        /// </remarks>
        Rank_Probability = 12,

        /// <summary>
        /// Peptide position
        /// </summary>
        Peptide_Position = 13,

        /// <summary>
        /// Q-Value
        /// </summary>
        QValue = 14
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _PrecursorInfo.txt file created by the Analysis Manager
    /// </summary>
    public enum PrecursorInfoFileColumns
    {
        /// <summary>
        /// Scan number
        /// </summary>
        ScanNumber = 0,

        /// <summary>
        /// Scan time
        /// </summary>
        ScanTime = 1,

        /// <summary>
        /// Scan type
        /// </summary>
        ScanType = 2,

        /// <summary>
        /// Scan type name
        /// </summary>
        ScanTypeName = 3,

        /// <summary>
        /// Precursor m/z
        /// </summary>
        PrecursorMz = 4,

        /// <summary>
        /// Scan filter text
        /// </summary>
        ScanFilterText = 5
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _ReporterIons.txt file (created by MASIC)
    /// </summary>
    /// <remarks>
    /// Note that the reporter ion intensity columns appear between the ReporterIonIntensityMax column and the WeightedAvgPctIntensityCorrection column
    /// </remarks>
    public enum ReporterIonsFileColumns
    {
        /// <summary>
        /// Dataset ID
        /// </summary>
        Dataset = 0,

        /// <summary>
        /// Scan Number
        /// </summary>
        ScanNumber = 1,

        /// <summary>
        /// Collision mode
        /// </summary>
        CollisionMode = 2,

        /// <summary>
        /// Parent ion m/z
        /// </summary>
        ParentIonMZ = 3,

        /// <summary>
        /// Base peak intensity
        /// </summary>
        BasePeakIntensity = 4,

        /// <summary>
        /// Base peak m/z
        /// </summary>
        BasePeakMZ = 5,

        /// <summary>
        /// Parent scan
        /// </summary>
        ParentScan = 6,

        /// <summary>
        /// Reporter ion intensity max
        /// </summary>
        ReporterIonIntensityMax = 7,

        /// <summary>
        /// Weighted average percent intensity correction
        /// </summary>
        WeightedAvgPctIntensityCorrection = 8
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _ScanStats.txt file created by MASIC or MSFileInfoScanner
    /// </summary>
    public enum ScanStatsFileColumns
    {
        /// <summary>
        /// Dataset ID
        /// </summary>
        DatasetId = 0,

        /// <summary>
        /// Scan number
        /// </summary>
        ScanNumber = 1,

        /// <summary>
        /// Scan time
        /// </summary>
        ScanTime = 2,

        /// <summary>
        /// Scan type
        /// </summary>
        ScanType = 3,

        /// <summary>
        /// Total ion intensity
        /// </summary>
        TotalIonIntensity = 4,

        /// <summary>
        /// Base peak intensity
        /// </summary>
        BasePeakIntensity = 5,

        /// <summary>
        /// Base peak m/z
        /// </summary>
        BasePeakMZ = 6,

        /// <summary>
        /// Base peak S/N ratio
        /// </summary>
        BasePeakSignalToNoiseRatio = 7,

        /// <summary>
        /// Ion count
        /// </summary>
        IonCount = 8,

        /// <summary>
        /// Ion count raw
        /// </summary>
        IonCountRaw = 9,

        /// <summary>
        /// Scan type name
        /// </summary>
        ScanTypeName = 10
    }

    /// <summary>
    /// These columns correspond to the tab-delimited file created from SEQUEST .out files
    /// </summary>
    public enum SequestSynopsisFileColumns
    {
        /// <summary>
        /// Row index
        /// </summary>
        RowIndex = 0,

        /// <summary>
        /// Scan
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Number of merged scans
        /// </summary>
        NumScans = 2,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 3,

        /// <summary>
        /// Peptide MH
        /// </summary>
        PeptideMH = 4,

        /// <summary>
        /// XCorr
        /// </summary>
        XCorr = 5,

        /// <summary>
        /// DeltaCn
        /// </summary>
        DeltaCn = 6,

        /// <summary>
        /// Sp
        /// </summary>
        Sp = 7,

        /// <summary>
        /// Protein name, aka Reference
        /// </summary>
        ProteinName = 8,

        /// <summary>
        /// Multiple protein count
        /// </summary>
        /// <remarks>
        /// Comes from MO = MultipleORFCount; this is 0 if the peptide is in just one protein; 1 if in 2 proteins, etc.
        /// </remarks>
        MultipleProteinCount = 9,

        /// <summary>
        /// Peptide sequence
        /// </summary>
        /// <remarks>
        /// This is the sequence with prefix and suffix residues and also with modification symbols
        /// </remarks>
        PeptideSequence = 10,

        /// <summary>
        /// DeltaCn
        /// </summary>
        DeltaCn2 = 11,

        /// <summary>
        /// RankSP
        /// </summary>
        RankSP = 12,

        /// <summary>
        /// RankXC
        /// </summary>
        RankXC = 13,

        /// <summary>
        /// Precursor error, in Daltons
        /// </summary>
        DelM = 14,

        /// <summary>
        /// XcRatio
        /// </summary>
        XcRatio = 15,

        /// <summary>
        /// Pass filter flag
        /// </summary>
        /// <remarks>
        /// Legacy/unused
        /// </remarks>
        PassFilt = 16,

        /// <summary>
        /// MScore
        /// </summary>
        /// <remarks>
        /// Legacy/unused
        /// </remarks>
        MScore = 17,

        /// <summary>
        /// Number of tryptic termini
        /// </summary>
        NTT = 18,

        /// <summary>
        /// Ions observed
        /// </summary>
        /// <remarks>
        /// Added in August 2011
        /// </remarks>
        IonsObserved = 19,

        /// <summary>
        /// Ions expected
        /// </summary>
        /// <remarks>
        /// Added in August 2011
        /// </remarks>
        IonsExpected = 20,

        /// <summary>
        /// Precursor error, in ppm
        /// </summary>
        /// <remarks>
        /// Added in August 2011
        /// </remarks>
        DelMPPM = 21,

        // The following column are computed by this program and appended to the input file or saved in a new file

        /// <summary>
        /// Cleavage state
        /// </summary>
        Cleavage_State = 22,

        /// <summary>
        /// Terminus state
        /// </summary>
        Terminus_State = 23,

        /// <summary>
        /// Mod count
        /// </summary>
        Mod_Count = 24,

        /// <summary>
        /// Mod description
        /// </summary>
        Mod_Description = 25,

        /// <summary>
        /// Monoisotopic mass
        /// </summary>
        Monoisotopic_Mass = 26
    }

    /// <summary>
    /// These columns correspond to the tab-delimited SICStats.txt file created by MASIC
    /// </summary>
    public enum SICStatsFileColumns
    {
        /// <summary>
        /// Dataset
        /// </summary>
        Dataset = 0,

        /// <summary>
        /// Parent ion index
        /// </summary>
        ParentIonIndex = 1,

        /// <summary>
        /// MZ
        /// </summary>
        MZ = 2,

        /// <summary>
        /// Survey scan number
        /// </summary>
        SurveyScanNumber = 3,

        /// <summary>
        /// Fragmentation scan number
        /// </summary>
        FragScanNumber = 4,

        /// <summary>
        /// Optimal peak apex scan number
        /// </summary>
        OptimalPeakApexScanNumber = 5,

        /// <summary>
        /// Peak apex override parent ion index
        /// </summary>
        PeakApexOverrideParentIonIndex = 6,

        /// <summary>
        /// CustomSICPeak
        /// </summary>
        CustomSICPeak = 7,

        /// <summary>
        /// Peak scan start
        /// </summary>
        PeakScanStart = 8,

        /// <summary>
        /// Peak scan end
        /// </summary>
        PeakScanEnd = 9,

        /// <summary>
        /// Peak scan max intensity
        /// </summary>
        PeakScanMaxIntensity = 10,

        /// <summary>
        /// Peak max intensity
        /// </summary>
        PeakMaxIntensity = 11,

        /// <summary>
        /// Peak S/N ratio
        /// </summary>
        PeakSignalToNoiseRatio = 12,

        /// <summary>
        /// Full-width at half-maximum, in scans
        /// </summary>
        FWHMInScans = 13,

        /// <summary>
        /// Peak area
        /// </summary>
        PeakArea = 14,

        /// <summary>
        /// Parent ion intensity
        /// </summary>
        ParentIonIntensity = 15,

        /// <summary>
        /// Peak baseline noise level
        /// </summary>
        PeakBaselineNoiseLevel = 16,

        /// <summary>
        /// Peak baseline noise StDev
        /// </summary>
        PeakBaselineNoiseStDev = 17,

        /// <summary>
        /// Peak baseline points used
        /// </summary>
        PeakBaselinePointsUsed = 18,

        /// <summary>
        /// Stat moments area
        /// </summary>
        StatMomentsArea = 19,

        /// <summary>
        /// Center of mass scan
        /// </summary>
        CenterOfMassScan = 20,

        /// <summary>
        /// Peak StDev
        /// </summary>
        PeakStDev = 21,

        /// <summary>
        /// Peak skew
        /// </summary>
        PeakSkew = 22,

        /// <summary>
        /// Peak KSStat
        /// </summary>
        PeakKSStat = 23,

        /// <summary>
        /// Stat moments data count used
        /// </summary>
        StatMomentsDataCountUsed = 24
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by TopPICResultsProcessor
    /// </summary>
    public enum TopPICSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Scan
        /// </summary>
        Scan = 1,

        /// <summary>
        /// Proteoform spectrum match ID
        /// </summary>
        Prsm_ID = 2,

        /// <summary>
        /// Spectrum ID
        /// </summary>
        Spectrum_ID = 3,

        /// <summary>
        /// Fragmentation method
        /// </summary>
        FragMethod = 4,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 5,

        /// <summary>
        /// Precursor m/z
        /// </summary>
        PrecursorMZ = 6,

        /// <summary>
        /// Precursor error, in Daltons
        /// </summary>
        DelM = 7,

        /// <summary>
        /// Precursor error, in ppm
        /// </summary>
        DelMPPM = 8,

        /// <summary>
        /// Theoretical monoisotopic peptide MH (computed by PHRP)
        /// </summary>
        /// <remarks>
        /// This is (M+H)+
        /// </remarks>
        MH = 9,

        /// <summary>
        /// Peptide sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
        /// </summary>
        Peptide = 10,

        /// <summary>
        /// Proteoform ID
        /// </summary>
        Proteoform_ID = 11,

        /// <summary>
        /// Feature intensity
        /// </summary>
        Feature_Intensity = 12,

        /// <summary>
        /// Feature score
        /// </summary>
        Feature_Score = 13,

        /// <summary>
        /// Feature apex time (in minutes)
        /// </summary>
        Feature_Apex_Time = 14,

        /// <summary>
        /// Number of proteins for this peptide (aka proteoform)
        /// </summary>
        Protein_Count = 15,

        /// <summary>
        /// Protein name
        /// </summary>
        Protein = 16,

        /// <summary>
        /// Protein N-terminus state
        /// </summary>
        /// <remarks>
        /// Examples:
        ///   NONE
        ///   M_ACETYLATION
        ///   NME
        ///   NME_ACETYLATION
        /// </remarks>
        Protein_Nterminal_Form = 17,

        /// <summary>
        /// Residue start
        /// </summary>
        ResidueStart = 18,

        /// <summary>
        /// Residue end
        /// </summary>
        ResidueEnd = 19,

        /// <summary>
        /// Unexpected Mod Count
        /// </summary>
        Unexpected_Mod_Count = 20,

        /// <summary>
        /// Peak Count
        /// </summary>
        Peak_Count = 21,

        /// <summary>
        /// Matched peak count
        /// </summary>
        Matched_Peak_Count = 22,

        /// <summary>
        /// Matched fragment ion count
        /// </summary>
        Matched_Fragment_Ion_Count = 23,

        /// <summary>
        /// P-value
        /// </summary>
        PValue = 24,

        /// <summary>
        /// Ranked score
        /// </summary>
        /// <remarks>
        /// Rank 1 means lowest (best) p-value, 2 means next higher p-value, etc. (ties get the same rank)
        /// </remarks>
        Rank_PValue = 25,

        /// <summary>
        /// E-Value
        /// </summary>
        EValue = 26,

        /// <summary>
        /// Q Value
        /// </summary>
        QValue = 27,

        /// <summary>
        /// Proteoform Q Value
        /// </summary>
        Proteoform_QValue = 28,

        /// <summary>
        /// Variable PTMs
        /// </summary>
        VariablePTMs = 29
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by XTandemResultsProcessor
    /// </summary>
    public enum XTandemSynFileColumns
    {
        /// <summary>
        /// Result ID
        /// </summary>
        ResultID = 0,

        /// <summary>
        /// Group ID
        /// </summary>
        GroupID = 1,

        /// <summary>
        /// Scan number
        /// </summary>
        Scan = 2,

        /// <summary>
        /// Charge
        /// </summary>
        Charge = 3,

        /// <summary>
        /// Monoisotopic (M+H)+ value, computed from PrecursorMZ and Charge
        /// </summary>
        MH = 4,

        /// <summary>
        /// Hyperscore
        /// </summary>
        Hyperscore = 5,

        /// <summary>
        /// Peptide expectation value LogE
        /// </summary>
        EValue = 6,

        /// <summary>
        /// Multiple protein count
        /// </summary>
        ProteinCount = 7,

        /// <summary>
        /// Peptide
        /// </summary>
        Peptide = 8,

        /// <summary>
        /// DeltaCn2
        /// </summary>
        DeltaCn2 = 9,

        /// <summary>
        /// Y ions score
        /// </summary>
        YScore = 10,

        /// <summary>
        /// Y ion count
        /// </summary>
        YIons = 11,

        /// <summary>
        /// B ions score
        /// </summary>
        BScore = 12,

        /// <summary>
        /// B ion count
        /// </summary>
        BIons = 13,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence
        /// </summary>
        DelM = 14,

        /// <summary>
        /// Peptide intensity LogI
        /// </summary>
        Intensity = 15,

        /// <summary>
        /// Precursor error, in ppm; corrected for isotope selection errors
        /// </summary>
        DelMPPM = 16
    }
}