using System;

// ReSharper disable UnusedMember.Global

namespace PHRPReader.Data
{
    // Ignore Spelling: Da, Daltons, MaxQuant, novo, tda, terminii, tryptic

#pragma warning disable 1591

    /// <summary>
    /// These columns correspond to the tab-delimited _ScanStatsEx.txt file created by MASIC or MSFileInfoScanner
    /// </summary>
    public enum ExtendedScanStatsFileColumns
    {
        DatasetId = 0,
        ScanNumber = 1,
        IonInjectionTime = 2,
        ScanEvent = 3,
        MasterIndex = 4,
        ElapsedScanTime = 5,
        ChargeState = 6,
        MonoisotopicMZ = 7,
        MS2IsolationWidth = 8,
        FTAnalyzerSettings = 9,
        FTAnalyzerMessage = 10,
        FTResolution = 11,
        ConversionParameterB = 12,
        ConversionParameterC = 13,
        ConversionParameterD = 14,
        ConversionParameterE = 15,
        CollisionMode = 16,
        ScanFilterText = 17,
        SourceVoltage = 18,
        Source_Current = 19
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by InSpecTResultsProcessor
    /// </summary>
    public enum InspectSynFileColumns
    {
        ResultID = 0,
        Scan = 1,
        Peptide = 2,
        Protein = 3,
        Charge = 4,
        MQScore = 5,
        Length = 6,
        TotalPRMScore = 7,
        MedianPRMScore = 8,
        FractionY = 9,
        FractionB = 10,
        Intensity = 11,
        NTT = 12,
        PValue = 13,
        FScore = 14,
        DeltaScore = 15,
        DeltaScoreOther = 16,
        DeltaNormMQScore = 17,                   // Computed as Abs((MQScore(n) - MQScore(n+1)) / MQScore(n)); storing 0 for the lowest scoring result in each set. If MQScore(n) is 0, stores 0.  This value is not usable when MQScore(n) is <= 0, and should generally not be used when MQScore(n) is < 0.5
        DeltaNormTotalPRMScore = 18,             // Computed as Abs((TotalPRMScore(n) - TotalPRMScore(n+1)) / TotalPRMScore(n)); storing 0 for the lowest scoring result in each set.  If TotalPRMScore(n) is 0, stores 0.  This value is not usable when TotalPRMScore(n) is <= 0, and should generally not be used when TotalPRMScore(n) is < 0.5
        RankTotalPRMScore = 19,                  // Rank 1 means highest TotalPRMScore, 2 means next lower score, etc. (ties get the same rank)
        RankFScore = 20,                         // Rank 1 means highest FScore, 2 means next lower, etc. (ties get the same rank)
        MH = 21,                                 // Theoretical monoisotopic peptide mass (computed by PHRP); note that this is (M+H)+
        RecordNumber = 22,
        DBFilePos = 23,
        SpecFilePos = 24,
        PrecursorMZ = 25,
        PrecursorError = 26,
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
        /// Precursor ion m/z (theoretical value, not observed value)
        /// </summary>
        /// <remarks>
        /// This is the theoretical m/z of the first isotope of the isotopic distribution of the parent ion
        /// </remarks>
        PrecursorMZ = 7,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence
        /// </summary>
        DelM = 8,

        /// <summary>
        /// Mass Error, in ppm
        /// </summary>
        DelM_PPM = 9,

        /// <summary>
        /// Mass error of the precursor ion equivalent monoisotopic mass value
        /// vs. the predicted monoisotopic mass of the identified peptide sequence,
        /// as computed by MaxQuant
        /// </summary>
        DelM_MaxQuant = 10,

        /// <summary>
        /// Mass Error, in ppm, as computed by MaxQuant
        /// </summary>
        DelM_PPM_MaxQuant = 11,

        /// <summary>
        /// Monoisotopic (M+H)+ value, computed from PrecursorMZ and Charge
        /// </summary>
        MH = 12,

        /// <summary>
        /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MaxQuant
        /// </summary>
        Mass = 13,

        /// <summary>
        /// Peptide sequence, including modification symbols
        /// </summary>
        Peptide = 14,

        /// <summary>
        /// Proteins
        /// </summary>
        /// <remarks>
        /// Semicolon separated list
        /// </remarks>
        Proteins = 15,

        /// <summary>
        /// Name of the best scoring protein this peptide is associated with
        /// </summary>
        /// <remarks>
        /// Typically there is only one protein name here
        /// However, in cases of a tied score, will be a semicolon separated list
        /// </remarks>
        LeadingRazorProtein = 16,

        /// <summary>
        /// Number of tryptic terminii
        /// </summary>
        NTT = 17,

        /// <summary>
        /// Posterior error probability
        /// </summary>
        /// <remarks>
        /// Similar to p-value
        /// Smaller values (closer to zero) are higher confidence
        /// </remarks>
        PEP = 18,

        /// <summary>
        ///  Andromeda score for the best MS/MS spectrum with this peptide
        /// </summary>
        /// <remarks>
        /// Higher scores are better
        /// </remarks>
        Score = 19,

        /// <summary>
        /// Score difference to the second best identified peptide with a different amino acid sequence
        /// </summary>
        DeltaScore = 20,

        /// <summary>
        /// Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
        /// </summary>
        Intensity = 21,

        /// <summary>
        /// Mass Analyzer of the instrument
        /// </summary>
        MassAnalyzer = 22,

        /// <summary>
        /// Type of precursor ion as identified by MaxQuant
        /// </summary>
        PrecursorType = 23,

        /// <summary>
        /// Elution time of the MS/MS spectrum
        /// </summary>
        RetentionTime = 24,

        /// <summary>
        /// Scan number where the precursor ion was observed
        /// </summary>
        PrecursorScan = 25,

        /// <summary>
        /// Intensity of the precursor ion in the scan that it was observed
        /// </summary>
        PrecursorIntensity = 26,

        /// <summary>
        /// Number of peaks (MS/MS ions) matching to the predicted fragmentation spectrum
        /// </summary>
        NumberOfMatches = 27,

        /// <summary>
        /// Fraction of intensity in the MS/MS spectrum that is annotated
        /// </summary>
        IntensityCoverage = 28,

        /// <summary>
        /// Number of missed enzymatic cleavages
        /// </summary>
        MissedCleavages = 29,

        /// <summary>
        /// Unique (consecutive) identifier for each row in msms.txt
        /// </summary>
        MsMsID = 30,

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
        ProteinGroupIDs = 31,

        /// <summary>
        /// The identifier of the non-redundant peptide sequence
        /// </summary>
        /// <remarks>
        /// Corresponds to the id column in the peptides.txt file
        /// </remarks>
        PeptideID = 32,

        /// <summary>
        /// Identifier referencing a row in the modificationSpecificPeptides.txt file
        /// </summary>
        ModPeptideID = 33,

        /// <summary>
        /// Identifier referencing a row in the evidence.txt file
        /// </summary>
        EvidenceID = 34
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MODPlusResultsProcessor
    /// </summary>
    public enum MODPlusSynFileColumns
    {
        ResultID = 0,
        Scan = 1,
        Spectrum_Index = 2,
        Charge = 3,
        PrecursorMZ = 4,
        DelM = 5,                            // Precursor error, in Daltons
        DelM_PPM = 6,                        // Precursor error, in ppm
        MH = 7,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
        Peptide = 8,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
        NTT = 9,
        ModificationAnnotation = 10,
        Protein = 11,
        Peptide_Position = 12,
        Score = 13,
        Probability = 14,
        Rank_Score = 15,
        QValue = 16
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by MSAlignResultsProcessor
    /// </summary>
    public enum MSAlignSynFileColumns
    {
        ResultID = 0,
        Scan = 1,
        Prsm_ID = 2,
        Spectrum_ID = 3,
        Charge = 4,
        PrecursorMZ = 5,
        DelM = 6,                            // Precursor error, in Da
        DelMPPM = 7,                         // Precursor error, in ppm
        MH = 8,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
        Peptide = 9,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
        Protein = 10,                        // Protein Name
        Protein_Mass = 11,
        Unexpected_Mod_Count = 12,
        Peak_Count = 13,
        Matched_Peak_Count = 14,
        Matched_Fragment_Ion_Count = 15,
        PValue = 16,
        Rank_PValue = 17,
        EValue = 18,
        FDR = 19,
        Species_ID = 20,
        FragMethod = 21
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _MSGF.txt file created by the analysis manager
    /// </summary>
    public enum MSGFFileColumns
    {
        ResultID = 0,
        Scan = 1,
        Charge = 2,
        Protein = 3,
        Peptide = 4,
        SpecProb = 5,
        Notes = 6
    }

#pragma warning restore 1591

    /// <summary>
    /// These columns correspond Synopsis files created for legacy tool MSGFDB
    /// </summary>
    public enum MSGFDBSynFileColumns
    {
        /// <summary>
        /// SpecProb
        /// </summary>
        SpecProb = 0,

        /// <summary>
        /// RankSpecProb
        /// </summary>
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
        /// PepFDR
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
        /// Precursor m/z
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

        /// <summary>
        /// Obsolete name for Spec E-Value (MS-GF+) or SpecProb (MSGFDB)
        /// </summary>
        [Obsolete("Use .SpecEValue")]
        SpecProb_EValue = SpecEValue,

        /// <summary>
        /// Rank SpecEValue
        /// </summary>
        RankSpecEValue = 15,

        /// <summary>
        /// Obsolete name for Rank Spec E-Value (MS-GF+) or Rank SpecProb (MSGFDB)
        /// </summary>
        [Obsolete("Use .RankSpecEValue")]
        RankSpecProb = RankSpecEValue,

        /// <summary>
        /// E-Value
        /// </summary>
        EValue = 16,

        /// <summary>
        /// Obsolete name for E-Value (MS-GF+) or P-Value (MSGFDB)
        /// </summary>
        [Obsolete("Use .EValue")]
        PValue_EValue = EValue,

        /// <summary>
        /// Q-Value
        /// </summary>
        /// <remarks>
        /// Only present if searched using -tda 1
        /// </remarks>
        QValue = 17,

        /// <summary>
        /// Obsolete name for Q-Value (MS-GF+) or FDR (MSGFDB)
        /// </summary>
        [Obsolete("Use .QValue")]
        FDR_QValue = QValue,

        /// <summary>
        /// Peptide QValue
        /// </summary>
        /// <remarks>
        /// Only present if searched using -tda 1
        /// </remarks>
        PepQValue = 18,

        /// <summary>
        /// Obsolete name for PepQValue (MS-GF+) or PepFDR (MSGFDB)
        /// </summary>
        [Obsolete("Use .PepQValue")]
        PepFDR_PepQValue = PepQValue,

        /// <summary>
        /// EFDR
        /// </summary>
        /// <remarks>Only present if did not search using -tda 1</remarks>
        EFDR = 19,

        // ReSharper disable UnusedMember.Global
        /// <summary>
        /// IMSScan
        /// </summary>
        /// <remarks>Only present for MSGFDB_IMS results</remarks>
        IMSScan = 20,

        /// <summary>
        /// IMSDriftTime
        /// </summary>
        /// <remarks>Only present for MSGFDB_IMS results</remarks>
        IMSDriftTime = 21,

        // ReSharper restore UnusedMember.Global
        /// <summary>
        /// IsotopeError
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
        /// Scan
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

#pragma warning disable 1591

    /// <summary>
    /// These columns correspond to the Synopsis file created by MODaResultsProcessor
    /// </summary>
    public enum MODaSynFileColumns
    {
        ResultID = 0,
        Scan = 1,
        Spectrum_Index = 2,
        Charge = 3,
        PrecursorMZ = 4,
        DelM = 5,                            // Precursor error, in Da
        DelM_PPM = 6,                        // Precursor error, in ppm
        MH = 7,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
        Peptide = 8,                         // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. +42
        Protein = 9,                         // Protein Name
        Score = 10,
        Probability = 11,
        Rank_Probability = 12,
        Peptide_Position = 13,
        QValue = 14
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _PrecursorInfo.txt file created by the Analysis Manager
    /// </summary>
    public enum PrecursorInfoFileColumns
    {
        ScanNumber = 0,
        ScanTime = 1,
        ScanType = 2,
        ScanTypeName = 3,
        PrecursorMz = 4,
        ScanFilterText = 5
    }

    /// <summary>
    /// These columns correspond to the tab-delimited _ScanStats.txt file created by MASIC or MSFileInfoScanner
    /// </summary>
    public enum ScanStatsFileColumns
    {
        DatasetId = 0,
        ScanNumber = 1,
        ScanTime = 2,
        ScanType = 3,
        TotalIonIntensity = 4,
        BasePeakIntensity = 5,
        BasePeakMZ = 6,
        BasePeakSignalToNoiseRatio = 7,
        IonCount = 8,
        IonCountRaw = 9,
        ScanTypeName = 10
    }

    /// <summary>
    /// These columns correspond to the tab-delimited file created from SEQUEST .out files
    /// </summary>
    public enum SequestSynopsisFileColumns
    {
        RowIndex = 0,
        Scan = 1,
        NumScans = 2,
        Charge = 3,
        PeptideMH = 4,
        XCorr = 5,
        DeltaCn = 6,
        Sp = 7,
        ProteinName = 8,                 // Aka Reference
        MultipleProteinCount = 9,        // Aka MO = MultipleORFCount; this is 0 if the peptide is in just one protein; 1 if in 2 proteins, etc.
        PeptideSequence = 10,            // This is the sequence with prefix and suffix residues and also with modification symbols
        DeltaCn2 = 11,
        RankSP = 12,
        RankXC = 13,
        DelM = 14,
        XcRatio = 15,
        PassFilt = 16,                   // Legacy/unused
        MScore = 17,                     // Legacy/unused
        NTT = 18,                        // Number of tryptic termini
        IonsObserved = 19,               // Added in August 2011
        IonsExpected = 20,               // Added in August 2011
        DelMPPM = 21,                    // Added in August 2011
        Cleavage_State = 22,             // This column and the ones after it are computed by this program and appended to the input file or saved in a new file
        Terminus_State = 23,
        Mod_Count = 24,
        Mod_Description = 25,
        Monoisotopic_Mass = 26
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by TopPICResultsProcessor
    /// </summary>
    public enum TopPICSynFileColumns
    {
        ResultID = 0,
        Scan = 1,
        Prsm_ID = 2,
        Spectrum_ID = 3,
        FragMethod = 4,
        Charge = 5,
        PrecursorMZ = 6,
        DelM = 7,                            // Precursor error, in Da
        DelMPPM = 8,                         // Precursor error, in ppm
        MH = 9,                              // Theoretical monoisotopic peptide MH (computed by PHRP); note that this is (M+H)+
        Peptide = 10,                        // This is the sequence with prefix and suffix residues and also with modification mass values, e.g. [42.01]
        Proteoform_ID = 11,
        Feature_Intensity = 12,
        Feature_Score = 13,
        Protein = 14,                        // Protein Name
        ResidueStart = 15,
        ResidueEnd = 16,
        Unexpected_Mod_Count = 17,
        Peak_Count = 18,
        Matched_Peak_Count = 19,
        Matched_Fragment_Ion_Count = 20,
        PValue = 21,
        Rank_PValue = 22,
        EValue = 23,
        QValue = 24,
        Proteoform_QValue = 25,
        VariablePTMs = 26
    }

    /// <summary>
    /// These columns correspond to the Synopsis file created by XTandemResultsProcessor
    /// </summary>
    public enum XTandemSynFileColumns
    {
        ResultID = 0,
        GroupID = 1,
        Scan = 2,
        Charge = 3,
        MH = 4,
        Hyperscore = 5,
        EValue = 6,                 // Peptide_Expectation_Value_LogE
        ProteinCount = 7,           // Multiple_Protein_Count
        Peptide = 8,
        DeltaCn2 = 9,
        YScore = 10,
        YIons = 11,
        BScore = 12,
        BIons = 13,
        DelM = 14,
        Intensity = 15,             // Peptide_Intensity_LogI
        DelMPPM = 16
    }

#pragma warning restore 1591

}