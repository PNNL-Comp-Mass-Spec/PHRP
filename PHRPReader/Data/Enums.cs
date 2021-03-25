namespace PHRPReader.Data
{
    // Ignore Spelling: Da, novo, tda

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
        /// SpecProb E-Value
        /// </summary>
        SpecProb_EValue = 14,

        /// <summary>
        /// Rank SpecProb
        /// </summary>
        /// <remarks>Rank 1 means lowest SpecEValue, 2 means next higher score, etc. (ties get the same rank)</remarks>
        RankSpecProb = 15,

        /// <summary>
        /// P-Value or E-Value
        /// </summary>
        PValue_EValue = 16,

        /// <summary>
        /// FDR or Q-Value
        /// </summary>
        /// <remarks>Only present if searched using -tda 1</remarks>
        FDR_QValue = 17,

        /// <summary>
        /// Peptide FDR or Peptide QValue
        /// </summary>
        /// <remarks>Only present if searched using -tda 1</remarks>
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
    public enum MSPathFinderSynFile
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
}