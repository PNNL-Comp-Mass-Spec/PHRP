using System;

namespace PHRPReader.Data
{
    // Ignore Spelling: Da, novo, tda

#pragma warning disable 1591

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
        /// Scan
        /// </summary>
        DatasetID = 1,

        /// <summary>
        /// Scan
        /// </summary>
        Scan = 2
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
        SpecProb_EValue = 14,

        /// <summary>
        /// Rank SpecEValue
        /// </summary>
        RankSpecEValue = 15,

        /// <summary>
        /// Obsolete name for Rank Spec E-Value (MS-GF+) or Rank SpecProb (MSGFDB)
        /// </summary>
        [Obsolete("Use .RankSpecEValue")]
        RankSpecProb = 15,

        /// <summary>
        /// E-Value
        /// </summary>
        EValue = 16,

        /// <summary>
        /// Obsolete name for E-Value (MS-GF+) or P-Value (MSGFDB)
        /// </summary>
        [Obsolete("Use .EValue")]
        PValue_EValue = 16,

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
        FDR_QValue = 17,

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
        PepFDR_PepQValue = 18,

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