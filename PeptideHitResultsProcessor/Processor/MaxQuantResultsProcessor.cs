// This class reads in a MaxQuant msms.txt results file and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.
//
// The modification definition information is determined from the MaxQuant parameter file
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// -------------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// This class reads in MaxQuant results files msms.txt and peptides.txt and creates
    /// a tab-delimited text file with the data.
    /// </summary>
    /// <remarks>
    /// <para>
    /// 1) ProcessFile reads MaxQuant results file msms.txt
    /// </para>
    /// <para>
    /// 2) It calls CreateSynResultsFile to create the _syn.txt file
    /// </para>
    /// <para>
    /// 3) ParseMaxQuantResultsFileHeaderLine reads the header line to determine the column mapping
    ///      columnMapping = new Dictionary of MaxQuantResultsFileColumns, int
    /// </para>
    /// <para>
    /// 4) ParseMaxQuantResultsFileEntry reads each data line and stores in an instance of MaxQuantSearchResult, which is a private structure
    ///    The data is stored in a list
    ///      searchResultsUnfiltered = new List of MaxQuantSearchResult
    /// </para>
    /// <para>
    /// 5) Once the entire .tsv has been read, searchResultsUnfiltered is sorted by scan, charge, and descending Andromeda score
    /// </para>
    /// <para>
    /// 6) StoreSynMatches stores filter-passing values in a new list
    ///      filteredSearchResults = new List of MaxQuantSearchResult
    /// </para>
    /// <para>
    /// 7) SortAndWriteFilteredSearchResults performs one more sort, then writes out to disk
    ///    Sorts descending Andromeda score, Scan, Peptide, and Razor Protein
    /// </para>
    /// </remarks>
    public class MaxQuantResultsProcessor : PHRPBaseClass
    {
        // ReSharper disable CommentTypo

        // Ignore Spelling: AaSubstitution, acetyl, Carbamidomethyl, conf, Cterm, Crosslink, Da, Dehydro, Desc, diff, diffs, DimethNter
        // Ignore Spelling: Glu, Gln, Glycan, maxq, MaxQuant, NeuCode, Nterm, Orbitrap, plex, pyro, struct, terminii, tryptic, txt

        // ReSharper restore CommentTypo

        /// <summary>
        /// Constructor
        /// </summary>
        public MaxQuantResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "May 21, 2021";

            MaxQuantMods = new Dictionary<string, MaxQuantModInfo>(StringComparer.OrdinalIgnoreCase);

            mPeptideCleavageStateCalculator = new PeptideCleavageStateCalculator();
        }

        public const string TOOL_NAME = "MaxQuant";

        public const string MSMS_FILE_NAME = "msms.txt";

        public const string PEPTIDES_FILE_NAME = "peptides.txt";

        /// <summary>
        /// Andromeda score threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
        /// </remarks>
        public const int DEFAULT_ANDROMEDA_SCORE_THRESHOLD = 50;

        /// <summary>
        /// PEP score threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
        /// </remarks>
        public const float DEFAULT_PEP_THRESHOLD = 0.01f;

        // ReSharper disable once UnusedMember.Global
        public const string N_TERMINUS_SYMBOL_MaxQuant = "_";

        // ReSharper disable once UnusedMember.Global
        public const string C_TERMINUS_SYMBOL_MaxQuant = "_";

        private const string SEARCH_ENGINE_NAME = "MaxQuant";

        private const int MAX_ERROR_MESSAGE_COUNT = 255;

        /// <summary>
        /// These columns correspond to MaxQuant file peptides.txt
        /// </summary>
        private enum MaxQuantPeptidesFileColumns
        {
            Sequence = 0,
            Prefix = 1,
            Suffix = 2,
            Proteins = 3,
            LeadingRazorProtein = 4,
            TotalPeptideIntensity = 5,
            Id = 6
        }

        /// <summary>
        /// These columns correspond to MaxQuant file msms.txt
        /// </summary>
        private enum MaxQuantResultsFileColumns
        {
            RawFile = 0,
            Scan = 1,
            ScanIndex = 2,
            Sequence = 3,
            Length = 4,
            MissedCleavageCount = 5,
            Modifications = 6,
            ModifiedSequence = 7,
            Proteins = 8,
            Charge = 9,
            Fragmentation = 10,
            MassAnalyzer = 11,
            PrecursorType = 12,
            ScanEventNumber = 13,
            IsotopeIndex = 14,
            MZ = 15,
            CalculatedMonoMass = 16,
            MassErrorPPM = 17,
            MassErrorDa = 18,
            SimpleMassErrorPPM = 19,
            RetentionTime = 20,
            PEP = 21,
            Score = 22,
            DeltaScore = 23,
            ScoreDiff = 24,
            LocalizationProb = 25,
            Combinatorics = 26,
            PIF = 27,
            FractionOfTotalSpectrum = 28,
            BasePeakFraction = 29,
            PrecursorScanNumber = 30,
            PrecursorIntensity = 31,
            PrecursorApexFraction = 32,
            PrecursorApexOffset = 33,
            PrecursorApexOffsetTime = 34,
            NumberOfMatches = 35,
            IntensityCoverage = 36,
            PeakCoverage = 37,
            Reverse = 38,
            ID = 39,
            ProteinGroupIDs = 40,
            PeptideID = 41,
            ModPeptideID = 42,
            EvidenceID = 43
        }

        private struct ResidueModificationInfo
        {
            /// <summary>
            /// Amino acid symbol
            /// </summary>
            public char Residue;

            /// <summary>
            /// Residue number in the peptide
            /// </summary>
            public int ResidueNumber;

            /// <summary>
            /// Tracks whether the residue is at a peptide or protein N or C terminus
            /// </summary>
            public AminoAcidModInfo.ResidueTerminusState ResidueTerminusState;

            /// <summary>
            /// Modification info
            /// </summary>
            public MSGFPlusParamFileModExtractor.ModInfo ModInfo;

            /// <summary>
            /// Show the residue symbol and mod mass
            /// </summary>
            /// <returns></returns>
            public override string ToString()
            {
                return string.Format("{0}: {1:F2}", Residue, ModInfo.ModMassVal);
            }
        }

        /// <summary>
        /// This data structure holds rows read from MaxQuant file msms.txt
        /// </summary>
        /// <remarks>
        /// These columns hold data that this class will use when creating the synopsis file
        /// </remarks>
        private struct MaxQuantSearchResult
        {
            /// <summary>
            /// Dataset name (typically the .raw file name, without the extension)
            /// </summary>
            public string DatasetName;

            /// <summary>
            /// MS/MS scan number in the dataset
            /// </summary>
            public string Scan;

            /// <summary>
            /// Integer value of Scan
            /// </summary>
            public int ScanNum;

            /// <summary>
            /// Index of the spectrum in the dataset (1-based, consecutive integer)
            /// </summary>
            public string ScanIndex;

            /// <summary>
            /// Identified peptide
            /// </summary>
            /// <remarks>
            /// Amino acid symbols only; no modification info
            /// </remarks>
            public string Sequence;

            /// <summary>
            /// Residue in the protein before this peptide
            /// </summary>
            /// <remarks>
            /// From the peptides.txt file, column "Amino acid before"
            /// </remarks>
            public string PrefixResidue;

            /// <summary>
            /// Residue in the protein after this peptide
            /// </summary>
            /// <remarks>
            /// From the peptides.txt file, column "Amino acid after"
            /// </remarks>
            public string SuffixResidue;

            /// <summary>
            /// Number of tryptic terminii
            /// </summary>
            /// <remarks>
            /// Computed by this class
            /// </remarks>
            public int NumberOfTrypticTerminii;

            /// <summary>
            /// Number of residues in the peptide
            /// </summary>
            public int Length;

            /// <summary>
            /// Number of missed enzymatic cleavages
            /// </summary>
            public string MissedCleavageCount;

            /// <summary>
            /// Post-translational modifications contained within the identified peptide
            /// This is read from the msms.txt file and includes dynamic modifications, but not static modifications
            /// </summary>
            /// <remarks>
            /// Examples
            /// <para>
            /// One modified residue:
            /// Oxidation (M)
            /// </para>
            /// <para>
            /// Two residues with the same modification
            /// 2 Oxidation (M)
            /// </para>
            /// <para>
            /// Two separate modifications (comma-separated list)
            /// Acetyl (Protein N-term),Oxidation (M)
            /// </para>
            /// </remarks>
            public string ModificationSummary;

            /// <summary>
            /// Comma-separated list of dynamic modification names and affected residue number
            /// This list is generated by this class
            /// </summary>
            /// <remarks>
            /// Examples:
            /// Oxidation 7
            /// Acetyl 1,Oxidation 4
            /// </remarks>
            public string Modifications { get; set; }

            /// <summary>
            /// Peptide sequence with embedded dynamic modifications
            /// </summary>
            /// <remarks>
            /// Examples:
            /// _YAEGYPGKR_
            /// _M(Oxidation (M))TSVGSQDTTGPMTR_
            /// _M(Oxidation (M))TSVGSQDTTGPM(Oxidation (M))TR_
            /// _(Acetyl (Protein N-term))M(Oxidation (M))GSHVAPTALTCAR_
            /// </remarks>
            public string ModifiedSequence;

            /// <summary>
            /// Protein associated with this peptide
            /// </summary>
            /// <remarks>
            /// Semicolon separated list
            /// </remarks>
            public string Proteins;

            /// <summary>
            /// Name of the best scoring protein this peptide is associated with
            /// </summary>
            /// <remarks>
            /// Typically there is only one protein name here
            /// However, in cases of a tied score, will be a semicolon separated list
            /// </remarks>
            public string LeadingRazorProtein;

            /// <summary>
            /// Charge state of the precursor ion
            /// </summary>
            public string Charge;

            /// <summary>
            /// Numeric value of Charge
            /// </summary>
            public short ChargeNum;

            /// <summary>
            /// Type of fragmentation used to create the MS/MS spectrum
            /// </summary>
            /// <remarks>
            /// Types:
            ///   CID - Collision Induced Dissociation
            ///   HCD - High energy Collision induced Dissociation
            ///   ETD - Electron Transfer Dissociation
            /// </remarks>
            public string Fragmentation;

            /// <summary>
            /// Mass Analyzer of the instrument
            /// </summary>
            /// <remarks>
            /// Types:
            ///   ITMS - Ion trap
            ///   FTMS - Fourier transform ICR or Orbitrap
            ///   TOF - Time of flight
            /// </remarks>
            public string MassAnalyzer;

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
            /// Examples:
            /// MULTI-MSMS
            /// MULTI-SECPEP
            /// MSMS
            /// </para>
            /// </remarks>
            public string PrecursorType;

            /// <summary>
            /// Scan event number
            /// </summary>
            /// <remarks>
            /// The first MS2 spectrum after a MS1 spectrum has event number 1
            /// The second MS2 spectrum has event number 2
            /// etc.
            /// Once the next MS1 spectrum is reached, scan event number resets to 1 for the next MS2 spectrum
            /// </remarks>
            public string ScanEventNumber;

            /// <summary>
            /// Isotope index
            /// </summary>
            public string IsotopeIndex;

            /// <summary>
            /// Precursor ion m/z (theoretical value, not observed value)
            /// </summary>
            /// <remarks>
            /// This is the theoretical m/z of the first isotope of the isotopic distribution of the parent ion
            /// For more information, see class <see cref="MaxQuantResults"/>
            /// </remarks>
            public string PrecursorMZ;

            /// <summary>
            /// Monoisotopic (M+H)+ value, computed from PrecursorMZ and Charge
            /// </summary>
            public string MH;

            /// <summary>
            /// Theoretical monoisotopic mass of the identified sequence (uncharged, including mods), as computed by MaxQuant
            /// </summary>
            public string CalculatedMonoMass;

            /// <summary>
            /// Numeric value of CalculatedMonoMass
            /// </summary>
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
            /// Mass error of the precursor ion equivalent monoisotopic mass value
            /// vs. the predicted monoisotopic mass of the identified peptide sequence,
            /// as computed by MaxQuant
            /// </summary>
            /// <remarks>
            /// This is NaN if a peptide has Type=MSMS
            /// In contrast, if a peptide has Type=MULTI-MSMS, a value will be defined
            /// </remarks>
            public string MassErrorPpmMaxQuant;

            /// <summary>
            /// Mass error, in Da, as computed by MaxQuant
            /// </summary>
            /// <remarks>
            /// This is NaN if a peptide has Type=MSMS
            /// In contrast, if a peptide has Type=MULTI-MSMS, a value will be defined
            /// </remarks>
            public string MassErrorDaMaxQuant;

            /// <summary>
            /// Simple mass error (ppm)
            /// </summary>
            /// <remarks>
            /// Definition not known
            /// </remarks>
            public string SimpleMassErrorPPM;

            /// <summary>
            /// Elution time of the MS/MS spectrum
            /// </summary>
            public string RetentionTime;

            /// <summary>
            /// Posterior error probability
            /// </summary>
            /// <remarks>
            /// Similar to p-value
            /// Smaller values (closer to zero) are higher confidence
            /// </remarks>
            public string PEP;

            /// <summary>
            /// Numeric value of PEP
            /// </summary>
            public double PEPValue;

            /// <summary>
            /// Andromeda score for the best MS/MS spectrum with this peptide
            /// </summary>
            /// <remarks>
            /// Higher scores are better
            /// </remarks>
            public string Score;

            /// <summary>
            /// Numeric value of Score
            /// </summary>
            public double ScoreValue;

            /// <summary>
            /// Score difference to the second best identified peptide with a different amino acid sequence
            /// </summary>
            public string DeltaScore;

            /// <summary>
            /// Score difference to the second best positioning of modifications identified peptide with the same amino acid sequence
            /// </summary>
            public string ScoreDiff;

            public int RankScore;

            /// <summary>
            /// PTM localization score
            /// </summary>
            public string LocalizationProb;

            /// <summary>
            /// Number of possible distributions of the modifications over the peptide sequence
            /// </summary>
            public string Combinatorics;

            /// <summary>
            /// Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
            /// </summary>
            /// <remarks>
            /// From the peptides.txt file, column "Intensity"
            /// </remarks>
            public string TotalPeptideIntensity;

            /// <summary>
            /// Parent Ion Fraction: the fraction of the target peak that makes up of the total intensity in the inclusion window
            /// </summary>
            public string PIF;

            /// <summary>
            /// Percentage the parent ion intensity makes up of the total intensity of the whole spectrum.
            /// </summary>
            public string FractionOfTotalSpectrum;

            /// <summary>
            /// Percentage the parent ion intensity in comparison to the highest peak in he MS spectrum
            /// </summary>
            public string BasePeakFraction;

            /// <summary>
            /// Scan number where the precursor ion was observed
            /// </summary>
            public string PrecursorScanNumber;

            /// <summary>
            /// Intensity of the precursor ion in the scan that it was observed
            /// </summary>
            public string PrecursorIntensity;

            /// <summary>
            /// Fraction the intensity of the precursor ion makes up of the peak (apex) intensity
            /// </summary>
            public string PrecursorApexFraction;

            /// <summary>
            /// The number of scans the precursor ion is offset from the peak (apex) position
            /// </summary>
            /// <remarks>
            /// For example, if the precursor scan is 3220 and the offset is -45, the peak apex is in scan 3265
            /// </remarks>
            public string PrecursorApexOffset;

            /// <summary>
            /// How much time the precursor ion is offset from the peak (apex) position
            /// </summary>
            public string PrecursorApexOffsetTime;

            /// <summary>
            /// Number of peaks (MS/MS ions) matching to the predicted fragmentation spectrum
            /// </summary>
            public string NumberOfMatches;

            /// <summary>
            /// Fraction of intensity in the MS/MS spectrum that is annotated
            /// </summary>
            public string IntensityCoverage;

            /// <summary>
            /// Fraction of peaks in the MS/MS spectrum that are annotated
            /// </summary>
            public string PeakCoverage;

            /// <summary>
            /// This will have '+' for peptides that are associated with a decoy protein
            /// </summary>
            public string Reverse;

            /// <summary>
            /// Unique (consecutive) identifier for each row in msms.txt
            /// </summary>
            /// <remarks>
            /// Used to cross-link the information in msms.txt with information stored in other files
            /// </remarks>
            public string MsMsID;

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
            public string ProteinGroupIDs;

            /// <summary>
            /// The identifier of the non-redundant peptide sequence
            /// </summary>
            /// <remarks>
            /// Corresponds to the id column in the peptides.txt file
            /// </remarks>
            public string PeptideID;

            /// <summary>
            /// Identifier referencing a row in the modificationSpecificPeptides.txt file
            /// </summary>
            public string ModPeptideID;

            /// <summary>
            /// Identifier referencing a row in the evidence.txt file
            /// </summary>
            public string EvidenceID;

            /// <summary>
            /// FDR
            /// </summary>
            /// <remarks>
            /// Computed by this class
            /// </remarks>
            public double FDR;

            /// <summary>
            /// Q-Value
            /// </summary>
            /// <remarks>
            /// Computed by this class
            /// </remarks>
            public double QValue;

            /// <summary>
            /// Reset all fields to default
            /// </summary>
            public void Clear()
            {
                DatasetName = string.Empty;
                Scan = string.Empty;
                ScanNum = 0;
                ScanIndex = string.Empty;
                Sequence = string.Empty;
                PrefixResidue = string.Empty;
                SuffixResidue = string.Empty;
                NumberOfTrypticTerminii = 0;
                Length = 0;
                MissedCleavageCount = string.Empty;
                ModificationSummary = string.Empty;
                Modifications = string.Empty;
                ModifiedSequence = string.Empty;
                Proteins = string.Empty;
                LeadingRazorProtein = string.Empty;
                Charge = string.Empty;
                ChargeNum = 0;
                Fragmentation = string.Empty;
                MassAnalyzer = string.Empty;
                PrecursorType = string.Empty;
                ScanEventNumber = string.Empty;
                IsotopeIndex = string.Empty;
                PrecursorMZ = string.Empty;
                MH = string.Empty;
                CalculatedMonoMass = string.Empty;
                CalculatedMonoMassValue = 0;
                CalculatedMonoMassPHRP = 0;
                MassErrorPpm = string.Empty;
                MassErrorDa = string.Empty;
                MassErrorPpmMaxQuant = string.Empty;
                MassErrorDaMaxQuant = string.Empty;
                SimpleMassErrorPPM = string.Empty;
                RetentionTime = string.Empty;
                PEP = string.Empty;
                PEPValue = 0;
                Score = string.Empty;
                ScoreValue = 0;
                DeltaScore = string.Empty;
                ScoreDiff = string.Empty;
                RankScore = 0;
                LocalizationProb = string.Empty;
                Combinatorics = string.Empty;
                TotalPeptideIntensity = string.Empty;
                PIF = string.Empty;
                FractionOfTotalSpectrum = string.Empty;
                BasePeakFraction = string.Empty;
                PrecursorScanNumber = string.Empty;
                PrecursorIntensity = string.Empty;
                PrecursorApexFraction = string.Empty;
                PrecursorApexOffset = string.Empty;
                PrecursorApexOffsetTime = string.Empty;
                NumberOfMatches = string.Empty;
                IntensityCoverage = string.Empty;
                PeakCoverage = string.Empty;
                Reverse = string.Empty;
                MsMsID = string.Empty;
                ProteinGroupIDs = string.Empty;
                PeptideID = string.Empty;
                ModPeptideID = string.Empty;
                EvidenceID = string.Empty;
                FDR = 0;
                QValue = 0;
            }

            public override string ToString()
            {
                return string.Format("Scan {0}: {1}, PEP {2}", Scan, Sequence, PEP);
            }
        }

        private readonly Regex mModListModNameMatcher = new(@"(?<ModName>.+) (?<ResidueNumber>\d+)", RegexOptions.Compiled);

        /// <summary>
        /// This RegEx matches dynamic modifications in ModificationSummary, e.g.
        /// Oxidation (M)
        /// 2 Oxidation (M)
        /// Acetyl (Protein N-term)
        /// </summary>
        private readonly Regex mModCountMatcher = new(@"^(?<ModCount>\d+) (?<ModName>.+)", RegexOptions.Compiled);

        /// <summary>
        /// Dictionary of MaxQuant modifications, loaded from modifications.xml
        /// </summary>
        /// <remarks>Keys are mod names, values are the details of each modification</remarks>
        private Dictionary<string, MaxQuantModInfo> MaxQuantMods { get; }

        private readonly PeptideCleavageStateCalculator mPeptideCleavageStateCalculator;

        /// <summary>
        /// Keys in this dictionary are amino acid symbol
        /// Values are the residue numbers of this amino acid in the peptide sequence
        /// </summary>
        private static readonly Dictionary<char, List<int>> mResiduePositions = new();

        /// <summary>
        /// Add dynamic modifications to a peptide read from the MaxQuant synopsis file
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modList"></param>
        private void AddDynamicModificationsToResidues(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(searchResult.Modifications))
            {
                return;
            }

            // Dynamic modifications in .Modifications are listed as a comma separated list of Mod name and residue number; examples:
            // Oxidation 21
            // Oxidation 11,Dehydro 12
            // Dehydro 1,Dehydro 4,Dehydro 7

            var mods = searchResult.Modifications.Split(',');
            var finalResidueLoc = searchResult.PeptideCleanSequence.Length;

            foreach (var modEntry in mods)
            {
                // Obtain the mod name, for example "Dehydro" from "Dehydro 52"
                var match = mModListModNameMatcher.Match(modEntry);

                if (!match.Success)
                {
                    ReportError("Invalid MaxQuant mod entry format; must be a name then a space then a number: " + modEntry);
                    continue;
                }

                var modName = match.Groups["ModName"].Value;
                var residueNumber = match.Groups["ResidueNumber"].Value;

                if (!GetMaxQuantModMass(modList, modName, true, out var modMass))
                {
                    ReportError(string.Format(
                        "Mod name '{0}' was not defined in the MaxQuant parameter file or in modifications.xml; " +
                        "cannot determine mod mass", modName));

                    continue;
                }

                if (!int.TryParse(residueNumber, out var residueLocInPeptide))
                {
                    ReportError("Mod entry does not have a number after the name: " + modEntry);
                    continue;
                }

                var residueTerminusState = AminoAcidModInfo.ResidueTerminusState.None;
                if (residueLocInPeptide <= 1)
                {
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
                }
                else if (residueLocInPeptide >= finalResidueLoc)
                {
                    residueTerminusState = AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
                }

                // Now that we know the terminus position, assure that residueLocInPeptide is 1 not 0
                if (residueLocInPeptide < 1)
                {
                    residueLocInPeptide = 1;
                }
                else if (residueLocInPeptide > finalResidueLoc)
                {
                    residueLocInPeptide = finalResidueLoc;
                }

                var chMostRecentResidue = '-';

                if (residueLocInPeptide >= 1 && residueLocInPeptide <= finalResidueLoc)
                {
                    chMostRecentResidue = searchResult.PeptideCleanSequence[residueLocInPeptide - 1];
                }

                // Associate the mod with the given residue
                searchResult.SearchResultAddModification(
                    modMass, chMostRecentResidue, residueLocInPeptide,
                    residueTerminusState, updateModOccurrenceCounts);
            }
        }

        /// <summary>
        /// Add modifications to a peptide read from the MaxQuant synopsis file
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modList"></param>
        /// <param name="staticModPresent">The calling method should set this to true if modList has static mods</param>
        private void AddModificationsToResidues(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList,
            bool staticModPresent)
        {
            // Add dynamic mods
            AddDynamicModificationsToResidues(searchResult, updateModOccurrenceCounts, modList);

            if (!staticModPresent)
                return;

            // Add static mods
            var staticModResidues = GetResiduesWithStaticMods(
                modList, searchResult.PeptideCleanSequence,
                searchResult.PeptidePreResidues, searchResult.PeptidePostResidues);

            foreach (var modifiedResidue in staticModResidues)
            {
                searchResult.SearchResultAddModification(
                    modifiedResidue.ModInfo.ModMassVal, modifiedResidue.Residue, modifiedResidue.ResidueNumber,
                    modifiedResidue.ResidueTerminusState, updateModOccurrenceCounts);
            }
        }

        /// <summary>
        /// Add modifications to a peptide read from the MaxQuant synopsis file
        /// Next, compute the monoisotopic mass
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modList"></param>
        /// <param name="staticModPresent">The calling method should set this to true if modList has static mods</param>
        /// <returns>True if success, false if an error</returns>
        private bool AddModificationsAndComputeMass(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList,
            bool staticModPresent)
        {
            try
            {
                // Some of the other tools add IsotopicMods here
                // This is not supported for MaxQuant
                //
                // searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts)

                // Parse .Modifications to determine the modified residues present
                AddModificationsToResidues(searchResult, updateModOccurrenceCounts, modList, staticModPresent);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();

                return true;
            }
            catch (Exception)
            {
                return false;
            }
        }

        /// <summary>
        /// Ranks each entry assumes all of the data is from the same scan)
        /// </summary>
        /// <param name="searchResults"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        private void AssignRankByScore(
            IList<MaxQuantSearchResult> searchResults,
            int startIndex,
            int endIndex)
        {
            if (startIndex == endIndex)
            {
                // Only one result
                var currentResult = searchResults[startIndex];
                currentResult.RankScore = 1;
                searchResults[startIndex] = currentResult;
                return;
            }

            // Duplicate a portion of searchResults so that we can sort by descending Andromeda Score

            var resultsSubset = new Dictionary<int, MaxQuantSearchResult>();
            for (var index = startIndex; index <= endIndex; index++)
            {
                resultsSubset.Add(index, searchResults[index]);
            }

            var resultsByScore = (from item in resultsSubset orderby item.Value.ScoreNum descending select item).ToList();

            double lastValue = 0;
            var currentRank = -1;

            foreach (var entry in resultsByScore)
            {
                var result = searchResults[entry.Key];

                if (currentRank < 0)
                {
                    lastValue = result.ScoreNum;
                    currentRank = 1;
                }
                else
                {
                    if (Math.Abs(result.ScoreNum - lastValue) > double.Epsilon)
                    {
                        lastValue = result.ScoreNum;
                        currentRank++;
                    }
                }

                result.RankScore = currentRank;
                searchResults[entry.Key] = result;
            }
        }

        private static string CombineDatasetNameParts(
            string datasetName,
            IReadOnlyList<string> datasetNameParts,
            int partCountToUse,
            int minimumLength = 0,
            int maximumLengthThreshold = 0)
        {
            if (partCountToUse >= datasetNameParts.Count)
                return datasetName;

            var nextIndex = 0;
            var combinedName = new StringBuilder();

            while (nextIndex < datasetNameParts.Count)
            {
                if (nextIndex >= partCountToUse)
                {
                    if (maximumLengthThreshold == 0 && minimumLength == 0)
                    {
                        // Minimum or maximum length not defined, and we've used the required number of parts
                        break;
                    }

                    var tooShort = minimumLength > 0 && combinedName.ToString().Length < minimumLength;

                    if (maximumLengthThreshold > 0)
                    {
                        // Maximum length defined
                        if (!tooShort && combinedName.ToString().Length + datasetNameParts[nextIndex].Length > maximumLengthThreshold)
                        {
                            // Adding the next part will result in the total length exceeding the maximum
                            // Do not add any more name parts
                            break;
                        }
                    }

                    if (!tooShort)
                    {
                        // Minimum length defined and the name is now long enough
                        break;
                    }
                }

                combinedName.Append(datasetNameParts[nextIndex]);
                nextIndex++;
            }

            return combinedName.ToString();
        }

        /// <summary>
        /// Compute the delta mass, in ppm, optionally correcting for C13 isotopic selection errors
        /// </summary>
        /// <param name="precursorErrorDa">Mass error (Observed - theoretical)</param>
        /// <param name="precursorMZ">Precursor m/z</param>
        /// <param name="charge">Precursor charge</param>
        /// <param name="peptideMonoisotopicMass"></param>
        /// <param name="adjustPrecursorMassForC13">Peptide's monoisotopic mass</param>
        /// <returns>DelM, in ppm</returns>
        /// <remarks>This function should only be called when column PMError(Da) is present (and PMError(ppm) is not present)</remarks>
        private double ComputeDelMCorrectedPPM(
            double precursorErrorDa,
            double precursorMZ,
            int charge,
            double peptideMonoisotopicMass,
            bool adjustPrecursorMassForC13)
        {
            // Compute the original value for the precursor monoisotopic mass
            var precursorMonoMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMZ, charge, 0);

            return SearchResultsBaseClass.ComputeDelMCorrectedPPM(precursorErrorDa, precursorMonoMass, adjustPrecursorMassForC13, peptideMonoisotopicMass);
        }

        /// <summary>
        /// Compute observed DelM and DelM_PPM values
        /// </summary>
        /// <param name="filteredSearchResults">Search results</param>
        /// <param name="precursorsByDataset">Keys are dataset names, values are dictionaries of precursor m/z by scan number</param>
        private void ComputeObservedMassErrors(
            IList<MaxQuantSearchResult> filteredSearchResults,
            IReadOnlyDictionary<string, Dictionary<int, double>> precursorsByDataset)
        {
            for (var i = 0; i < filteredSearchResults.Count; i++)
            {
                var searchResult = filteredSearchResults[i];

                if (!precursorsByDataset.TryGetValue(searchResult.DatasetName, out var precursorsByScan))
                    continue;

                if (!precursorsByScan.TryGetValue(searchResult.ScanNum, out var precursorMz))
                    continue;

                var observedPrecursorMass = mPeptideSeqMassCalculator.ConvoluteMass(precursorMz, searchResult.ChargeNum, 0);

                var deltaMassDa = observedPrecursorMass - searchResult.CalculatedMonoMassPHRP;

                var peptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(deltaMassDa, precursorMz,
                    searchResult.ChargeNum, searchResult.CalculatedMonoMassPHRP,
                    true);

                searchResult.MassErrorPpm = PRISM.StringUtilities.DblToString(peptideDeltaMassCorrectedPpm, 5, 0.00005);

                var precursorErrorDa = PeptideMassCalculator.PPMToMass(peptideDeltaMassCorrectedPpm, searchResult.CalculatedMonoMassPHRP);

                // Note that this will be a C13-corrected precursor error; not the absolute precursor error
                searchResult.MassErrorDa = StringUtilities.MassErrorToString(precursorErrorDa);

                // Update the value in the list of structs
                filteredSearchResults[i] = searchResult;
            }
        }

        private double ComputePeptideMass(string cleanSequence, double totalModMass)
        {
            var mass = mPeptideSeqMassCalculator.ComputeSequenceMass(cleanSequence);
            mass += totalModMass;

            return mass;
        }

        /// <summary>
        /// Computes the total of all modifications defined for the sequence
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="modList"></param>
        /// <param name="staticModPresent">The calling method should set this to true if modList has static mods</param>
        private double ComputeTotalModMass(
            MaxQuantSearchResult searchResult,
            IList<MSGFPlusParamFileModExtractor.ModInfo> modList,
            bool staticModPresent)
        {
            if (string.IsNullOrWhiteSpace(searchResult.ModificationSummary))
            {
                return 0;
            }

            double totalModMass = 0;

            // Parse the dynamic modification info stored in ModificationSummary, e.g.
            // "Acetyl (Protein N-term),Oxidation (M)" or "2 Oxidation (M)"

            foreach (var modItem in searchResult.ModificationSummary.Split(','))
            {
                if (modItem.Equals("Unmodified", StringComparison.OrdinalIgnoreCase))
                    continue;

                // Convert the mod name to a mass value

                // Check whether the mod name is preceded by an integer,
                // for example the 2 in "2 Oxidation (M)"

                var match = mModCountMatcher.Match(modItem);

                string modName;
                int modCount;

                if (match.Success)
                {
                    modName = match.Groups["ModName"].Value;
                    modCount = int.Parse(match.Groups["ModCount"].Value);
                }
                else
                {
                    modName = modItem;
                    modCount = 1;
                }

                if (GetMaxQuantModMass(modList, modName, false, out var modMass))
                {
                    totalModMass += modCount * modMass;
                    continue;
                }

                ReportError(string.Format(
                    "Mod name '{0}' was not defined in the MaxQuant parameter file or in modifications.xml; " +
                    "cannot determine mod mass", modName));
            }

            if (!staticModPresent)
                return totalModMass;

            // Add static mods
            var staticModResidues = GetResiduesWithStaticMods(
                modList, searchResult.Sequence,
                searchResult.PrefixResidue, searchResult.SuffixResidue);

            foreach (var modifiedResidue in staticModResidues)
            {
                totalModMass += modifiedResidue.ModInfo.ModMassVal;
            }

            return totalModMass;
        }

        /// <summary>
        /// Parse out dynamic modifications from ModifiedSequence (read from msms.txt)
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="modList"></param>
        /// <returns>Comma separated list of dynamic modification names and affected residues</returns>
        private string ConstructDynamicModificationList(
            MaxQuantSearchResult searchResult,
            IList<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(searchResult.ModificationSummary) || searchResult.ModificationSummary.Equals("Unmodified", StringComparison.OrdinalIgnoreCase))
            {
                return string.Empty;
            }

            var modsWithCounts = searchResult.ModificationSummary.Split(',');

            // This dictionary tracks information about mod names in ModificationSummary and ModifiedSequence
            // Keys are MaxQuant mod name, values are instances of MaxQuantModifiedSequenceModInfo
            var modifiedSeqMods = new Dictionary<string, MaxQuantModifiedSequenceModInfo>();

            foreach (var modItem in modsWithCounts)
            {
                // Check whether the mod name is preceded by an integer (which indicates the number of times this peptide has this modification)
                var match = mModCountMatcher.Match(modItem);

                var modName = match.Success ? match.Groups["ModName"].Value : modItem;

                if (!modifiedSeqMods.ContainsKey(modName))
                {
                    modifiedSeqMods.Add(modName, new MaxQuantModifiedSequenceModInfo(modName));
                }
            }

            // ReSharper disable CommentTypo

            // Determine the character index and residue number of each residue in searchResult.ModifiedSequence
            // Example strings
            // _YAEGYPGKR_
            // _M(Oxidation (M))TSVGSQDTTGPMTR_
            // _M(Oxidation (M))TSVGSQDTTGPM(Oxidation (M))TR_
            // _(Acetyl (Protein N-term))M(Oxidation (M))GSHVAPTALTCAR_

            // ReSharper restore CommentTypo

            // To do this, copy ModifiedSequence to a new variable and replace the modification descriptions (inside parentheses) with underscores
            var maskedModifiedSequence = searchResult.ModifiedSequence;

            foreach (var modItem in modifiedSeqMods)
            {
                var modWithParentheses = string.Format("({0})", modItem.Value.MaxQuantModName);

                var charIndex = maskedModifiedSequence.IndexOf(modWithParentheses, StringComparison.Ordinal);
                if (charIndex >= 0)
                {
                    maskedModifiedSequence = maskedModifiedSequence.Replace(modWithParentheses, new string('_', modWithParentheses.Length));
                    modItem.Value.MatchedModName = modItem.Value.MaxQuantModName;
                    continue;
                }

                // ReSharper disable CommentTypo

                // Mod not found
                // The results might be from an older version of MaxQuant that used abbreviated names
                // For example: _(ac)AAAAAAGAGPEM(ox)VR_

                // ReSharper restore CommentTypo

                var altModName = modItem.Value.MaxQuantModName.Substring(0, 2).ToLower();
                var altModNameWithParentheses = string.Format("({0})", altModName);

                var charIndex2 = maskedModifiedSequence.IndexOf(altModNameWithParentheses, StringComparison.Ordinal);
                if (charIndex2 >= 0)
                {
                    maskedModifiedSequence = maskedModifiedSequence.Replace(altModNameWithParentheses, new string('_', altModNameWithParentheses.Length));
                    modItem.Value.MatchedModName = altModName;
                    continue;
                }

                OnWarningEvent(string.Format(
                    "Did not find modification name {0} in modified sequence {1}; also tried {2}",
                    modWithParentheses, searchResult.ModifiedSequence, altModNameWithParentheses));
            }

            // Assure that there are no parentheses remaining in the masked sequence
            if (maskedModifiedSequence.IndexOfAny(new[] { '(', ')' }) >= 0)
            {
                OnWarningEvent(string.Format(
                    "Masked modified sequence still has parentheses: {0} for {1}",
                    maskedModifiedSequence, searchResult.ModifiedSequence));
            }

            // Populate a dictionary mapping character index in maskedModifiedSequence to residue number
            var charIndexResidueNumberMap = new Dictionary<int, int>();
            var currentResidueNumber = 0;

            for (var i = 0; i < maskedModifiedSequence.Length; i++)
            {
                if (char.IsLetter(maskedModifiedSequence, i))
                    currentResidueNumber++;

                charIndexResidueNumberMap.Add(i, currentResidueNumber);
            }

            if (currentResidueNumber != searchResult.Length)
            {
                OnWarningEvent(string.Format(
                    "Final residue number determined for peptide does not match amino acid count: {0} instead of {1} for {2}, originally {3}",
                    currentResidueNumber, searchResult.Length, maskedModifiedSequence, searchResult.ModifiedSequence));
            }

            // This list holds modification names and affected residue number, separated by a space
            // Modification names do not include affected residues
            // For example, if the MaxQuant modification name is "Methyl (KR)" and the affected residue is the 5th amino acid in the protein,
            // we append "Methyl 5" to modificationList
            var dynamicModifications = new List<string>();

            // This dictionary maps from the abbreviated modification name stored in the synopsis file to the the full modification name (as used by MaxQuant)
            var modificationNameMap = new Dictionary<string, string>();

            // Look for each dynamic modification and store it in dynamicModifications
            foreach (var modItem in modifiedSeqMods)
            {
                var matchedModNameWithParentheses = modItem.Value.MatchedModNameWithParentheses;

                var startIndex = 0;
                while (startIndex < matchedModNameWithParentheses.Length)
                {
                    var charIndex = searchResult.ModifiedSequence.IndexOf(matchedModNameWithParentheses, startIndex, StringComparison.Ordinal);

                    if (charIndex < 0)
                        break;

                    // Examine this modification's metadata to see if it is an N-terminal mod
                    var nTerminalMod = MaxQuantMods.TryGetValue(modItem.Value.MaxQuantModName, out var modInfo) &&
                                       modInfo.Position is MaxQuantModPosition.AnyNterm or MaxQuantModPosition.ProteinNterm;

                    // ReSharper disable CommentTypo

                    int residueNumber;

                    if (nTerminalMod)
                    {
                        // The mod applies to the residue just after the mod name (which is surrounded by parentheses)
                        // For example Acetyl in:
                        // _(Acetyl (Protein N-term))AAAAAAGAGPEM(Oxidation (M))VR_
                        // or
                        // _(ac)AAAAAAGAGPEM(ox)VR_

                        if (charIndex + matchedModNameWithParentheses.Length >= charIndexResidueNumberMap.Count)
                        {
                            OnWarningEvent(string.Format(
                                "{0} found at the end of ModifiedSequence {1}; this is unexpected",
                                matchedModNameWithParentheses, searchResult.ModifiedSequence));
                            residueNumber = searchResult.Length;
                        }
                        else
                        {
                            residueNumber = charIndexResidueNumberMap[charIndex + matchedModNameWithParentheses.Length];
                        }
                    }
                    else if (charIndex == 0)
                    {
                        // This shouldn't happen, but we'll check for it
                        OnWarningEvent(string.Format(
                            "{0} found at the start of ModifiedSequence {1}; this is unexpected",
                            matchedModNameWithParentheses, searchResult.ModifiedSequence));
                        residueNumber = 0;
                    }
                    else
                    {
                        // The mod applies to the residue just before the mod name (which is surrounded by parentheses)
                        // For example, Oxidation in
                        // _HAM(Oxidation (M))LLQTLQLM(Oxidation (M))HK_

                        residueNumber = charIndexResidueNumberMap[charIndex - 1];
                    }

                    // ReSharper restore CommentTypo

                    // Append the modification to searchResult.Modifications

                    var maxQuantModName = modItem.Key;
                    var shortModName = modItem.Value.GetModNameWithoutResidues();

                    dynamicModifications.Add(string.Format("{0} {1}", shortModName, residueNumber));

                    if (modificationNameMap.TryGetValue(shortModName, out var modificationNameToCompare))
                    {
                        if (!maxQuantModName.Equals(modificationNameToCompare))
                        {
                            OnWarningEvent(string.Format(
                                "Multiple MaxQuant modifications have the same short modification name (with the residue removed); " +
                                "short name {0} resolves to both {1} and {2}",
                                shortModName, modificationNameToCompare, maxQuantModName));
                        }
                    }
                    else
                    {
                        modificationNameMap.Add(shortModName, maxQuantModName);
                    }

                    var existingModFound = false;
                    for (var i = 0; i < modList.Count; i++)
                    {
                        var paramFileMod = modList[i];
                        if (!paramFileMod.ModName.Equals(maxQuantModName, StringComparison.OrdinalIgnoreCase))
                        {
                            continue;
                        }

                        if (string.IsNullOrWhiteSpace(modList[i].ShortName))
                        {
                            paramFileMod.ShortName = shortModName;

                            // Update the value in the list of structs
                            modList[i] = paramFileMod;
                        }

                        existingModFound = true;
                        break;
                    }

                    if (!existingModFound)
                    {
                        OnWarningEvent(string.Format("Modification {0} not found in modList; this is unexpected", maxQuantModName));
                        var modDef = GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod, maxQuantModName);
                        modList.Add(modDef);
                    }

                    startIndex = charIndex + 1;
                }
            }

            return string.Join(",", dynamicModifications);
        }

        private bool CreateProteinModsFileWork(
            string baseName,
            FileInfo inputFile,
            string synOutputFilePath,
            string outputDirectoryPath)
        {
            bool success;

            // Create the MTSPepToProteinMap file

            var baseNameFilePath = Path.Combine(inputFile.DirectoryName ?? string.Empty, baseName);
            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(baseNameFilePath, outputDirectoryPath, mts: true);

            var sourcePHRPDataFiles = new List<string>();

            if (!string.IsNullOrEmpty(synOutputFilePath))
            {
                sourcePHRPDataFiles.Add(synOutputFilePath);
            }

            if (sourcePHRPDataFiles.Count == 0)
            {
                SetErrorMessage("Cannot call CreatePepToProteinMapFile since sourcePHRPDataFiles is empty");
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                success = false;
            }
            else
            {
                if (File.Exists(mtsPepToProteinMapFilePath) && Options.UseExistingMTSPepToProteinMapFile)
                {
                    success = true;
                }
                else
                {
                    success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                    if (!success)
                    {
                        OnWarningEvent("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                    }
                }
            }

            if (success)
            {
                if (inputFile.Directory == null)
                {
                    OnWarningEvent("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                }
                else if (string.IsNullOrWhiteSpace(synOutputFilePath))
                {
                    OnWarningEvent("CreateProteinModsFileWork: synOutputFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(synOutputFilePath)), outputDirectoryPath);

                    // Create the Protein Mods file
                    success = CreateProteinModDetailsFile(synOutputFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          PeptideHitResultTypes.MaxQuant);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                return true;
            }

            return true;
        }

        /// <summary>
        /// This routine creates a synopsis file from the output from MaxQuant (file msms.txt)
        /// The synopsis file includes every result with a probability above a set threshold
        /// </summary>
        /// <param name="maxQuantPeptides"></param>
        /// <param name="inputFilePath">MaxQuant results file (msms.txt)</param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="modList"></param>
        /// <param name="baseName">Output: base synopsis file name</param>
        /// <param name="synOutputFilePath">Output: synopsis file path created by this method</param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            IReadOnlyDictionary<string, MaxQuantPeptideInfo> maxQuantPeptides,
            string inputFilePath,
            string outputDirectoryPath,
            IList<MSGFPlusParamFileModExtractor.ModInfo> modList,
            out string baseName,
            out string synOutputFilePath)
        {
            baseName = string.Empty;
            synOutputFilePath = string.Empty;

            try
            {
                var columnMapping = new Dictionary<MaxQuantResultsFileColumns, int>();
                var errorMessages = new List<string>();

                var inputFile = new FileInfo(inputFilePath);
                if (inputFile.Directory == null)
                {
                    SetErrorMessage("Unable to determine the parent directory of file " + inputFile.FullName);
                    SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                    return false;
                }

                OnStatusEvent("Reading MaxQuant results file, " + PathUtils.CompactPathString(inputFile.FullName, 80));

                // Open the input file and parse it
                using var reader = new StreamReader(new FileStream(inputFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var headerParsed = false;
                var lineNumber = 0;

                // Initialize the list that will hold all of the records in the MaxQuant result file
                var searchResultsUnfiltered = new List<MaxQuantSearchResult>();

                // Initialize the list that will hold all of the records that will ultimately be written out to disk
                var filteredSearchResults = new List<MaxQuantSearchResult>();

                // Parse the input file
                while (!reader.EndOfStream && !AbortProcessing)
                {
                    var lineIn = reader.ReadLine();
                    lineNumber++;

                    if (string.IsNullOrWhiteSpace(lineIn))
                    {
                        continue;
                    }

                    if (!headerParsed)
                    {
                        // Parse the header line
                        var success = ParseMaxQuantResultsFileHeaderLine(lineIn, columnMapping);
                        if (!success)
                        {
                            if (string.IsNullOrEmpty(mErrorMessage))
                            {
                                SetErrorMessage("Invalid header line in " + inputFile.Name);
                            }

                            return false;
                        }

                        headerParsed = true;
                        continue;
                    }

                    var validSearchResult = ParseMaxQuantResultsFileEntry(
                        maxQuantPeptides, lineIn, out var searchResult,
                        errorMessages, columnMapping, modList, lineNumber);

                    if (validSearchResult)
                    {
                        searchResultsUnfiltered.Add(searchResult);
                    }

                    // Update the progress
                    UpdateSynopsisFileCreationProgress(reader);
                }

                // Sort the SearchResults by scan, charge, and descending Andromeda score
                searchResultsUnfiltered.Sort(new MaxQuantSearchResultsComparerScanChargeScorePeptide());

                // Now filter the data
                var startIndex = 0;

                while (startIndex < searchResultsUnfiltered.Count)
                {
                    // Find all of the matches for the current result's scan
                    // MaxQuant will typically report just one match
                    var endIndex = startIndex;
                    while (endIndex + 1 < searchResultsUnfiltered.Count &&
                           searchResultsUnfiltered[endIndex + 1].ScanNum == searchResultsUnfiltered[startIndex].ScanNum)
                    {
                        endIndex++;
                    }

                    // Store the results for this scan
                    StoreSynMatches(searchResultsUnfiltered, startIndex, endIndex, filteredSearchResults);

                    startIndex = endIndex + 1;
                }

                Console.WriteLine();

                // Keys in this dictionary are dataset names, values are abbreviated names
                var baseNameByDatasetName = GetDatasetNameMap(filteredSearchResults, out var longestCommonBaseName);

                // Load Precursor m/z masses (if an appropriate file can be found)
                // If found, compute MassErrorPpm and MassErrorDa
                var precursorMzLoader = new ScanPrecursorMzLoader(inputFile.Directory.FullName);
                RegisterEvents(precursorMzLoader);

                var precursorsByDataset = precursorMzLoader.LoadPrecursorIons(baseNameByDatasetName.Keys.ToList());

                if (precursorsByDataset.Count > 0)
                {
                    ComputeObservedMassErrors(filteredSearchResults, precursorsByDataset);
                }

                // The synopsis file name will be of the form DatasetName_maxq_syn.txt
                // If baseDatasetNames only has one item, will use the full dataset name
                // If baseDatasetNames has multiple items, will use the longest string in common for the keys in baseDatasetNames

                baseName = baseNameByDatasetName.Count switch
                {
                    0 => "Dataset_maxq",
                    1 => baseNameByDatasetName.First().Key + "_maxq",
                    _ => longestCommonBaseName + "_maxq"
                };

                synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SYNOPSIS_FILE_SUFFIX);

                OnStatusEvent("Creating synopsis file, " + PathUtils.CompactPathString(synOutputFilePath, 80));

                using var writer = new StreamWriter(new FileStream(synOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                // Write the header line to the output file
                WriteSynFHTFileHeader(writer, errorMessages);

                // Sort the data in filteredSearchResults then write out to disk
                SortAndWriteFilteredSearchResults(baseNameByDatasetName, writer, filteredSearchResults, errorMessages);

                // Inform the user if any errors occurred
                if (errorMessages.Count > 0)
                {
                    SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateSynResultsFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the MaxQuant parameter file
        /// </summary>
        /// <param name="maxQuantParamFilePath">
        /// This is typically the XML-based MaxQuant parameter file, but it can alternatively be
        /// the parameters.txt file created in the txt output directory</param>
        /// <param name="modList"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromParamFile(
            string maxQuantParamFilePath,
            out List<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            try
            {
                var sourceFile = new FileInfo(maxQuantParamFilePath);
                if (!sourceFile.Exists)
                {
                    SetErrorMessage("MaxQuant parameter file not found: " + maxQuantParamFilePath);
                    SetErrorCode(PHRPErrorCode.ParameterFileNotFound);
                    modList = new List<MSGFPlusParamFileModExtractor.ModInfo>();
                    return false;
                }

                if (sourceFile.Name.Equals("parameters.txt", StringComparison.OrdinalIgnoreCase))
                    return ExtractModInfoFromTxtParameterFile(sourceFile, out modList);

                if (sourceFile.Extension.Equals(".xml", StringComparison.OrdinalIgnoreCase))
                    return ExtractModInfoFromXmlParameterFile(sourceFile, out modList);

                SetErrorMessage("MaxQuant parameter should either have a .xml extension or be named parameters.txt; " +
                                "unsupported file: " + maxQuantParamFilePath);

                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                modList = new List<MSGFPlusParamFileModExtractor.ModInfo>();
                return false;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ExtractModInfoFromParamFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                modList = new List<MSGFPlusParamFileModExtractor.ModInfo>();
                return false;
            }
        }

        private bool ExtractModInfoFromTxtParameterFile(FileInfo sourceFile, out List<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            modList = new List<MSGFPlusParamFileModExtractor.ModInfo>();

            try
            {
                // Open the parameters.txt file and look for static mod names
                // The parameters.txt file does not list dynamic mods
                // It also does not list reporter ion based mods (e.g. TMT)

                OnStatusEvent("Reading the MaxQuant parameters.txt file in " + PathUtils.CompactPathString(sourceFile.Directory?.FullName, 80));

                using var reader = new StreamReader(new FileStream(sourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var lineNumber = 0;

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    lineNumber++;

                    if (string.IsNullOrWhiteSpace(lineIn) || !lineIn.StartsWith("Fixed modifications", StringComparison.OrdinalIgnoreCase))
                        continue;

                    // Found the static mods line, e.g.
                    // Fixed modifications	Carbamidomethyl (C)

                    var lineParts = lineIn.Split('\t');

                    if (lineParts.Length < 2)
                    {
                        OnWarningEvent(string.Format(
                            "Line number {0} in the parameters.txt file has 'Fixed modifications' but does not have a tab character; this is unexpected",
                            lineNumber));

                        continue;
                    }

                    var staticMods = lineParts[1];

                    // ToDo: determine how multiple static mods are listed

                    foreach (var maxQuantModName in staticMods.Split(','))
                    {
                        var modDef = GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod, maxQuantModName);
                        modList.Add(modDef);
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ExtractModInfoFromTxtParameterFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        private bool ExtractModInfoFromXmlParameterFile(FileSystemInfo sourceFile, out List<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            modList = new List<MSGFPlusParamFileModExtractor.ModInfo>();

            try
            {
                // Open the MaxQuant parameter file with an XML reader and look for static and dynamic mod names
                OnStatusEvent("Reading the MaxQuant parameter file, " + PathUtils.CompactPathString(sourceFile.FullName, 80));

                using var reader = new StreamReader(new FileStream(sourceFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                // Note that XDocument supersedes XmlDocument and XPathDocument
                // XDocument can often be easier to use since XDocument is LINQ-based

                var doc = XDocument.Parse(reader.ReadToEnd());

                // <fixedModifications>
                //    <string>Carbamidomethyl (C)</string>
                // </fixedModifications>

                // <variableModifications>
                //    <string>Oxidation (M)</string>
                //    <string>Acetyl (Protein N-term)</string>
                // </variableModifications>

                // <isobaricLabels>
                //   <IsobaricLabelInfo>
                //     <internalLabel>TMT6plex-Lys126</internalLabel>
                //     <terminalLabel>TMT6plex-Nter126</terminalLabel>
                //     <correctionFactorM2>0</correctionFactorM2>
                //     <correctionFactorM1>0</correctionFactorM1>
                //     <correctionFactorP1>0</correctionFactorP1>
                //     <correctionFactorP2>0</correctionFactorP2>
                //     <tmtLike>True</tmtLike>
                //   </IsobaricLabelInfo>
                //   ...
                // </isobaricLabels>

                var parameterGroupNodes = doc.Elements("MaxQuantParams").Elements("parameterGroups").Elements("parameterGroup").ToList();

                if (parameterGroupNodes.Count == 0)
                {
                    OnWarningEvent("MaxQuant parameter file is missing the <parameterGroup> element; cannot add modification symbols to peptides in the synopsis file");
                    return false;
                }

                if (parameterGroupNodes.Count > 1)
                {
                    OnWarningEvent("MaxQuant parameter file has more than one <parameterGroup> element; this is unexpected");
                }

                var fixedModNodes = parameterGroupNodes.Elements("fixedModifications").Elements("string").ToList();
                var dynamicModNodes = parameterGroupNodes.Elements("variableModifications").Elements("string").ToList();

                // Check for reporter ion mods, e.g. 6-plex or 10-plex TMT
                var internalIsobaricLabelNodes = parameterGroupNodes.Elements("isobaricLabels").Elements("IsobaricLabelInfo").Elements("internalLabel").ToList();
                var terminalIsobaricLabelNodes = parameterGroupNodes.Elements("isobaricLabels").Elements("IsobaricLabelInfo").Elements("terminalLabel").ToList();

                foreach (var fixedMod in fixedModNodes)
                {
                    var maxQuantModName = fixedMod.Value;
                    var modDef = GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod, maxQuantModName);
                    modList.Add(modDef);
                }

                foreach (var dynamicMod in dynamicModNodes)
                {
                    var maxQuantModName = dynamicMod.Value;
                    var modDef = GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod, maxQuantModName);
                    modList.Add(modDef);
                }

                if (internalIsobaricLabelNodes.Count > 0)
                {
                    var maxQuantModName = internalIsobaricLabelNodes[0].Value;
                    var modDef = GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod, maxQuantModName);
                    modList.Add(modDef);
                }

                if (terminalIsobaricLabelNodes.Count > 0)
                {
                    var maxQuantModName = terminalIsobaricLabelNodes[0].Value;
                    var modDef = GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod, maxQuantModName);
                    modList.Add(modDef);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ExtractModInfoFromXmlParameterFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to string.Empty
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        protected bool GetColumnValueCheckNaN(string[] splitLine, int columnIndex, out string value)
        {
            var valueDefined = GetColumnValue(splitLine, columnIndex, out value, string.Empty);

            if (valueDefined && value.Equals("NaN", StringComparison.OrdinalIgnoreCase))
                value = string.Empty;

            return valueDefined;
        }

        /// <summary>
        /// Examine the Raw File names in filteredSearchResults
        /// Create a mapping from full name to abbreviated name
        /// </summary>
        /// <param name="filteredSearchResults"></param>
        /// <param name="longestCommonBaseName"></param>
        /// <returns>Dictionary where keys are dataset names and values are abbreviated names</returns>
        private Dictionary<string, string> GetDatasetNameMap(IEnumerable<MaxQuantSearchResult> filteredSearchResults, out string longestCommonBaseName)
        {
            var datasetNames = new SortedSet<string>();

            foreach (var item in filteredSearchResults)
            {
                var datasetName = item.DatasetName;

                if (datasetNames.Contains(datasetName))
                    continue;

                datasetNames.Add(datasetName);
            }

            var baseNameByDatasetName = GetDatasetNameMap(datasetNames, out longestCommonBaseName, out var warnings);

            foreach (var warning in warnings)
            {
                OnWarningEvent(warning);
            }

            return baseNameByDatasetName;
        }

        /// <summary>
        /// Examine the names in datasetNames
        /// Create a mapping from full name to abbreviated name
        /// </summary>
        /// <param name="datasetNames"></param>
        /// <param name="longestCommonBaseName">Output: longest common base name</param>
        /// <param name="warnings">Output: warning messages</param>
        /// <returns>Dictionary where keys are dataset names and values are abbreviated names</returns>
        public static Dictionary<string, string> GetDatasetNameMap(
            SortedSet<string> datasetNames,
            out string longestCommonBaseName,
            out List<string> warnings)
        {
            warnings = new List<string>();

            var datasetNameParts = new Dictionary<string, List<string>>();
            var maxPartCount = 0;
            var splitChars = new[] { '_', '-' };

            foreach (var datasetName in datasetNames)
            {
                var nameParts = new List<string>();
                var startIndex = 0;
                while (startIndex < datasetName.Length)
                {
                    var matchIndex = datasetName.IndexOfAny(splitChars, startIndex + 1);
                    if (matchIndex < 0)
                    {
                        nameParts.Add(datasetName.Substring(startIndex));
                        break;
                    }

                    nameParts.Add(datasetName.Substring(startIndex, matchIndex - startIndex));
                    startIndex = matchIndex;
                }

                datasetNameParts.Add(datasetName, nameParts);

                maxPartCount = Math.Max(maxPartCount, nameParts.Count);
            }

            if (datasetNameParts.Count == 0)
            {
                longestCommonBaseName = string.Empty;
                return new Dictionary<string, string>();
            }

            var candidateBaseNames = new SortedSet<string>();

            var partCountToUse = 1;

            var datasetNameKeys = datasetNameParts.Keys.ToList();

            while (partCountToUse <= maxPartCount)
            {
                candidateBaseNames.Clear();
                candidateBaseNames.Add(CombineDatasetNameParts(datasetNameKeys[0], datasetNameParts[datasetNameKeys[0]], partCountToUse));

                for (var i = 1; i < datasetNameKeys.Count; i++)
                {
                    var baseNameToAdd = CombineDatasetNameParts(datasetNameKeys[i], datasetNameParts[datasetNameKeys[i]], partCountToUse);
                    if (candidateBaseNames.Contains(baseNameToAdd))
                    {
                        // Name collision found
                        break;
                    }

                    candidateBaseNames.Add(baseNameToAdd);
                }

                if (candidateBaseNames.Count == datasetNameKeys.Count)
                    break;

                partCountToUse++;
            }

            var baseDatasetNames = new SortedSet<string>();

            // Dictionary where keys are dataset names and values are abbreviated names
            var baseNameByDatasetName = new Dictionary<string, string>();

            if (candidateBaseNames.Count == datasetNameKeys.Count)
            {
                // Can use a subsection of the dataset name(s)
                // Combine subsections to create the base name for each dataset
                foreach (var item in datasetNameParts)
                {
                    var baseNameToAdd = CombineDatasetNameParts(item.Key, item.Value, partCountToUse, 12, 25);
                    baseNameByDatasetName.Add(item.Key, baseNameToAdd);

                    if (baseDatasetNames.Contains(baseNameToAdd))
                    {
                        warnings.Add(string.Format(
                            "Warning: baseDatasetNames already contains: {0}\nLogic error for dataset {1}",
                            baseNameToAdd, item.Key));

                        continue;
                    }

                    baseDatasetNames.Add(baseNameToAdd);
                }
            }
            else
            {
                // Not able to shorten the dataset names since they are too similar
                // Use full dataset names
                foreach (var item in datasetNameParts)
                {
                    baseNameByDatasetName.Add(item.Key, item.Key);
                    baseDatasetNames.Add(item.Key);
                }
            }

            longestCommonBaseName = StringUtilities.LongestCommonStringFromStart(baseNameByDatasetName.Values.ToList());
            longestCommonBaseName = longestCommonBaseName.TrimEnd('_', '-');

            if (longestCommonBaseName.Length > 7 && longestCommonBaseName.EndsWith("_0"))
            {
                longestCommonBaseName = longestCommonBaseName.Substring(0, longestCommonBaseName.Length - 2);
            }

            return baseNameByDatasetName;
        }

        /// <summary>
        /// Look for the modification, by name, in modList, which has modification info loaded from the MaxQuant parameter file
        /// If not found, look in MaxQuantMods, which has modification info loaded from modifications.xml
        /// </summary>
        /// <param name="modList">List of modifications loaded from the MaxQuant parameter file (or from parameters.txt)</param>
        /// <param name="modName">Modification name to find</param>
        /// <param name="matchShortName">When true, modName is the short modification name</param>
        /// <param name="modMass">Output: modification mass</param>
        /// <returns>True if found, otherwise false</returns>
        private bool GetMaxQuantModMass(
            IEnumerable<MSGFPlusParamFileModExtractor.ModInfo> modList,
            string modName,
            bool matchShortName,
            out double modMass)
        {
            foreach (var modItem in modList)
            {
                if (matchShortName)
                {
                    if (!string.Equals(modItem.ShortName, modName, StringComparison.OrdinalIgnoreCase))
                    {
                        continue;
                    }
                }
                else
                {
                    if (!string.Equals(modItem.ModName, modName, StringComparison.OrdinalIgnoreCase))
                    {
                        continue;
                    }
                }

                modMass = modItem.ModMassVal;
                return true;
            }

            if (MaxQuantMods.TryGetValue(modName, out var modInfo))
            {
                modMass = modInfo.MonoisotopicMass;
                return true;
            }

            if (matchShortName)
            {
                foreach (var modItem in MaxQuantMods)
                {
                    var shortModName = MaxQuantModifiedSequenceModInfo.GetModNameWithoutResidues(modItem.Key);
                    if (string.Equals(shortModName, modName, StringComparison.OrdinalIgnoreCase))
                    {
                        modMass = modItem.Value.MonoisotopicMass;
                        return true;
                    }
                }
            }

            OnWarningEvent(string.Format("Modification {0} not found in modList or MaxQuantMods", modName));
            modMass = 0;
            return false;
        }

        private MSGFPlusParamFileModExtractor.ModInfo GetModDetails(MSGFPlusParamFileModExtractor.MSGFPlusModType modType, string maxQuantModName)
        {
            var modDef = new MSGFPlusParamFileModExtractor.ModInfo
            {
                ModType = modType,
                ModTypeInParameterFile = modType
            };

            if (!MaxQuantMods.TryGetValue(maxQuantModName, out var modInfo))
            {
                OnWarningEvent(string.Format(
                    "The MaxQuant modifications.xml file does not have a mod named '{0}'; unrecognized mod name in the MaxQuant parameter file",
                    maxQuantModName));

                return modDef;
            }

            modDef.ModName = modInfo.Title;
            modDef.ModMass = modInfo.MonoisotopicMass.ToString(CultureInfo.InvariantCulture);
            modDef.ModMassVal = modInfo.MonoisotopicMass;
            modDef.Residues = modInfo.Residues;
            modDef.ModSymbol = MSGFPlusParamFileModExtractor.UNKNOWN_MSGFPlus_MOD_SYMBOL;

            switch (modInfo.Position)
            {
                case MaxQuantModPosition.Anywhere:
                case MaxQuantModPosition.NotCterm:
                    // Leave .ModType unchanged; this is a static or dynamic mod
                    break;

                case MaxQuantModPosition.AnyNterm:
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod;
                    }

                    modDef.Residues = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynNTermPeptide;

                    break;

                case MaxQuantModPosition.AnyCterm:
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod;
                    }

                    modDef.Residues = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynCTermPeptide;

                    break;

                case MaxQuantModPosition.ProteinNterm:
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod;
                    }

                    modDef.Residues = AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynNTermProtein;

                    break;

                case MaxQuantModPosition.ProteinCterm:
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod;
                    }

                    modDef.Residues = AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynCTermProtein;

                    break;

                default:
                    throw new ArgumentOutOfRangeException();
            }

            return modDef;
        }

        private string GetPeptideSequence(MaxQuantSearchResult peptideInfo, bool includePrefixAndSuffix = true)
        {
            if (!includePrefixAndSuffix)
                return peptideInfo.Sequence;

            return peptideInfo.PrefixResidue + "." + peptideInfo.Sequence + "." + peptideInfo.SuffixResidue;
        }

        private static IEnumerable<ResidueModificationInfo> GetResiduesWithStaticMods(
            IEnumerable<MSGFPlusParamFileModExtractor.ModInfo> modList,
            string cleanSequence,
            string prefixResidue,
            string suffixResidue)
        {
            // Construct a list of residues with static modification
            // If a residue has multiple static mods, it will be added multiple times

            mResiduePositions.Clear();

            for (var residueNumber = 1; residueNumber <= cleanSequence.Length; residueNumber++)
            {
                if (mResiduePositions.TryGetValue(cleanSequence[residueNumber - 1], out var residuePositions))
                {
                    residuePositions.Add(residueNumber);
                    continue;
                }

                mResiduePositions.Add(cleanSequence[residueNumber - 1], new List<int> { residueNumber });
            }

            var staticModResidues = new List<ResidueModificationInfo>();

            foreach (var modItem in modList)
            {
                if (modItem.ModTypeInParameterFile != MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod)
                    continue;

                foreach (var residueSymbol in modItem.Residues)
                {
                    if (residueSymbol == '-')
                    {
                        // Apply this static mod to all residues
                        // This is most likely invalid
                        for (var residueNumber = 1; residueNumber <= cleanSequence.Length; residueNumber++)
                        {
                            foreach (var item in mResiduePositions)
                            {
                                foreach (var residuePosition in item.Value)
                                {
                                    var terminusState = GetResidueTerminusState(cleanSequence.Length, prefixResidue, suffixResidue, residuePosition);

                                    staticModResidues.Add(new ResidueModificationInfo
                                    {
                                        Residue = item.Key,
                                        ResidueNumber = residuePosition,
                                        ResidueTerminusState = terminusState,
                                        ModInfo = modItem
                                    });
                                }
                            }
                        }

                        continue;
                    }

                    if (!mResiduePositions.TryGetValue(residueSymbol, out var residuePositions))
                        continue;

                    // ReSharper disable once ForeachCanBeConvertedToQueryUsingAnotherGetEnumerator
                    foreach (var residuePosition in residuePositions)
                    {
                        var terminusState = GetResidueTerminusState(cleanSequence.Length, prefixResidue, suffixResidue, residuePosition);

                        staticModResidues.Add(new ResidueModificationInfo
                        {
                            Residue = residueSymbol,
                            ResidueNumber = residuePosition,
                            ResidueTerminusState = terminusState,
                            ModInfo = modItem
                        });
                    }
                }
            }

            return staticModResidues;
        }

        private static AminoAcidModInfo.ResidueTerminusState GetResidueTerminusState(
            int cleanSequenceLength, string prefixResidue, string suffixResidue, int residuePosition)
        {
            if (residuePosition == 1)
            {
                if (!string.IsNullOrWhiteSpace(prefixResidue) && prefixResidue == "-")
                    return AminoAcidModInfo.ResidueTerminusState.ProteinNTerminus;

                return AminoAcidModInfo.ResidueTerminusState.PeptideNTerminus;
            }

            if (residuePosition == cleanSequenceLength)
            {
                if (!string.IsNullOrWhiteSpace(suffixResidue) && suffixResidue == "-")
                    return AminoAcidModInfo.ResidueTerminusState.ProteinCTerminus;

                return AminoAcidModInfo.ResidueTerminusState.PeptideCTerminus;
            }

            return AminoAcidModInfo.ResidueTerminusState.None;
        }

        /// <summary>
        /// Look for the MaxQuant modifications.xml file, both in the input directory and in other predefined paths
        /// </summary>
        /// <param name="inputDirectory"></param>
        /// <returns>True if no error (including file not found; which is a warning); false if an exception</returns>
        private bool LoadMaxQuantModificationDefinitions(DirectoryInfo inputDirectory)
        {
            const string DEFAULT_MODIFICATION_FILE_PATH = @"C:\MaxQuant\bin\conf\modifications.xml";
            const string DEFAULT_MODIFICATION_FILE_PATH_DMS = @"C:\DMS_Programs\MaxQuant\bin\conf\modifications.xml";

            try
            {
                FileInfo modificationDefinitionFile;
                var candidateFiles = inputDirectory.GetFiles("modifications.xml").ToList();
                if (candidateFiles.Count > 0)
                {
                    modificationDefinitionFile = candidateFiles[0];
                }
                else
                {
                    var candidateFile1 = new FileInfo(DEFAULT_MODIFICATION_FILE_PATH);
                    var candidateFile2 = new FileInfo(DEFAULT_MODIFICATION_FILE_PATH_DMS);

                    if (candidateFile1.Exists)
                    {
                        modificationDefinitionFile = candidateFile1;
                    }
                    else if (candidateFile2.Exists)
                    {
                        modificationDefinitionFile = candidateFile2;
                    }
                    else
                    {
                        OnWarningEvent(string.Format(
                            "MaxQuant modifications file not found; cannot add modification symbols to peptides in the synopsis file. " +
                            "File modifications.xml should either be in the input directory or in the default MaxQuant directory: \n" +
                            "{0} or\n{1} or \n{2}",
                            inputDirectory.FullName, candidateFile1.FullName, candidateFile2.FullName));

                        OnWarningEvent("Download the default modifications.xml file from ");
                        return true;
                    }
                }

                OnStatusEvent("Loading MaxQuant modifications from " + modificationDefinitionFile.FullName);

                MaxQuantMods.Clear();

                using var reader = new StreamReader(new FileStream(modificationDefinitionFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                // Note that XDocument supersedes XmlDocument and XPathDocument
                // XDocument can often be easier to use since XDocument is LINQ-based

                var doc = XDocument.Parse(reader.ReadToEnd());

                // ReSharper disable CommentTypo
                //   <modification title="Carbamidomethyl (C)" description="Iodoacetamide derivative" create_date="2010-01-19T15:54:02.1976309+01:00" last_modified_date="2010-01-19T15:55:13.1353165+01:00" user="neuhause" reporterCorrectionM2="0" reporterCorrectionM1="0" reporterCorrectionP1="0" reporterCorrectionP2="0" reporterCorrectionType="false" composition="C(2) H(3) N O">
                //      <position>anywhere</position>
                //      <modification_site site="C">
                //         <neutralloss_collection />
                //         <diagnostic_collection />
                //      </modification_site>
                //      <type>Standard</type>
                //      <terminus_type>none</terminus_type>
                //   </modification>

                //  <modification title="Acetyl (K)" description="Acetylation" create_date="2010-01-19T13:25:24.82357+01:00" last_modified_date="2010-01-19T14:15:40.2445414+01:00" user="neuhause" reporterCorrectionM2="0" reporterCorrectionM1="0" reporterCorrectionP1="0" reporterCorrectionP2="0" reporterCorrectionType="false" composition="C(2) H(2) O">
                //     <position>notCterm</position>
                //     <modification_site site="K">
                //        <neutralloss_collection />
                //        <diagnostic_collection>
                //           <diagnostic name="acK*" shortname="acK*" composition="C(7) H(11) O N" />
                //           <diagnostic name="acK" shortname="acK" composition="C(7) H(14) O N(2)" />
                //        </diagnostic_collection>
                //     </modification_site>
                //     <type>Standard</type>
                //     <terminus_type>none</terminus_type>
                //  </modification>

                // ReSharper restore CommentTypo

                foreach (var modificationNode in doc.Elements("modifications").Elements("modification"))
                {
                    if (!XmlReaderUtilities.TryGetAttribute(modificationNode, "title", out var modTitle))
                    {
                        OnWarningEvent("Modification node in the MaxQuant modifications file is missing the title attribute");
                        continue;
                    }

                    if (!XmlReaderUtilities.TryGetAttribute(modificationNode, "composition", out var composition))
                    {
                        OnWarningEvent("Modification node in the MaxQuant modifications file is missing the composition attribute");
                        continue;
                    }

                    XmlReaderUtilities.TryGetAttribute(modificationNode, "description", out var modDescription);

                    if (MaxQuantMods.ContainsKey(modTitle))
                    {
                        OnWarningEvent(string.Format(
                            "MaxQuant modifications file has a duplicate definition for the modification titled '{0}'",
                            modTitle));

                        continue;
                    }

                    var maxQuantMod = new MaxQuantModInfo(modTitle, composition)
                    {
                        Description = modDescription
                    };

                    MaxQuantMods.Add(modTitle, maxQuantMod);

                    if (XmlReaderUtilities.TryGetElementValue(modificationNode, "position", out var positionText))
                    {
                        maxQuantMod.Position = positionText switch
                        {
                            "anywhere" => MaxQuantModPosition.Anywhere,
                            "anyNterm" => MaxQuantModPosition.AnyNterm,
                            "anyCterm" => MaxQuantModPosition.AnyCterm,
                            "notCterm" => MaxQuantModPosition.NotCterm,
                            "proteinNterm" => MaxQuantModPosition.ProteinNterm,
                            "proteinCterm" => MaxQuantModPosition.ProteinCterm,
                            _ => throw new ArgumentOutOfRangeException()
                        };
                    }
                    else
                    {
                        OnWarningEvent(string.Format(
                            "Modification node '{0}' in the MaxQuant modifications file is missing the 'position' element",
                            maxQuantMod.Title));
                    }

                    if (XmlReaderUtilities.TryGetElementValue(modificationNode, "type", out var modTypeText))
                    {
                        maxQuantMod.ModType = modTypeText switch
                        {
                            "Standard" => MaxQuantModType.Standard,
                            "IsobaricLabel" => MaxQuantModType.IsobaricLabel,
                            "Label" => MaxQuantModType.Label,
                            "NeuCodeLabel" => MaxQuantModType.NeuCodeLabel,
                            "Glycan" => MaxQuantModType.Glycan,
                            "AaSubstitution" => MaxQuantModType.AaSubstitution,
                            "CleavedCrosslink" => MaxQuantModType.CleavedCrosslink,
                            _ => throw new ArgumentOutOfRangeException()
                        };
                    }
                    else
                    {
                        OnWarningEvent(string.Format(
                            "Modification node '{0}' in the MaxQuant modifications file is missing the 'type' element",
                            maxQuantMod.Title));
                    }

                    if (XmlReaderUtilities.TryGetElementValue(modificationNode, "terminus_type", out var terminusTypeText))
                    {
                        maxQuantMod.TerminusType = terminusTypeText switch
                        {
                            "none" => MaxQuantTerminusType.None,
                            "nterm" => MaxQuantTerminusType.NTerminus,
                            _ => throw new ArgumentOutOfRangeException()
                        };
                    }
                    else
                    {
                        OnWarningEvent(string.Format(
                            "Modification node '{0}' in the MaxQuant modifications file is missing the 'terminus_type' element",
                            maxQuantMod.Title));
                    }

                    var modificationSiteNodes = modificationNode.Elements("modification_site").ToList();
                    if (modificationSiteNodes.Count == 0)
                    {
                        OnWarningEvent(string.Format(
                            "Modification node '{0}' in the MaxQuant modifications file 'modification_site' elements",
                            maxQuantMod.Title));

                        continue;
                    }

                    foreach (var modificationSiteNode in modificationSiteNodes)
                    {
                        if (XmlReaderUtilities.TryGetAttribute(modificationSiteNode, "site", out var modificationSite))
                        {
                            if (string.IsNullOrWhiteSpace(modificationSite))
                            {
                                OnWarningEvent(string.Format(
                                    "Modification node '{0}' in the MaxQuant modifications file has a modification_site element with an empty string for the 'site' attribute",
                                    maxQuantMod.Title));
                            }
                            else
                            {
                                var residue = modificationSite[0];
                                if (maxQuantMod.Residues.Contains(residue))
                                {
                                    OnWarningEvent(string.Format(
                                        "Modification node '{0}' in the MaxQuant modifications file has multiple modification_site elements with the same 'site' attribute value",
                                        maxQuantMod.Title));
                                    continue;
                                }

                                maxQuantMod.Residues += residue;
                            }
                        }
                        else
                        {
                            OnWarningEvent(string.Format(
                                "Modification node '{0}' in the MaxQuant modifications file is missing the 'site' attribute for the modification_site element",
                                maxQuantMod.Title));
                        }
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadMaxQuantModificationDefinitions", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        /// <summary>
        /// Read data from the peptides.txt file
        /// </summary>
        /// <param name="inputDirectory"></param>
        /// <param name="maxQuantPeptides">Keys are peptide sequence, values are the metadata for the peptide</param>
        /// <remarks>
        /// For column details, see http://www.coxdocs.org/doku.php?id=maxquant:table:peptidetable
        /// </remarks>
        private void LoadPeptideInfo(FileSystemInfo inputDirectory, out Dictionary<string, MaxQuantPeptideInfo> maxQuantPeptides)
        {
            maxQuantPeptides = new Dictionary<string, MaxQuantPeptideInfo>();

            try
            {
                var inputFile = new FileInfo(Path.Combine(inputDirectory.FullName, PEPTIDES_FILE_NAME));
                if (!inputFile.Exists)
                {
                    OnWarningEvent("MaxQuant peptides.txt file not found: " + inputFile.FullName);
                    return;
                }

                OnStatusEvent("Reading peptides from file " + PathUtils.CompactPathString(inputFile.FullName, 80));

                using var reader = new StreamReader(new FileStream(inputFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var headerParsed = false;

                var columnMapping = new Dictionary<MaxQuantPeptidesFileColumns, int>();

                // Keys are column indices, values are experiment names
                var intensityByExperimentColumns = new Dictionary<int, string>();
                var lineNumber = 0;

                while (!reader.EndOfStream && !AbortProcessing)
                {
                    var lineIn = reader.ReadLine();
                    lineNumber++;

                    if (string.IsNullOrWhiteSpace(lineIn))
                    {
                        continue;
                    }

                    if (!headerParsed)
                    {
                        var validHeader = ParseMaxQuantPeptidesFileHeaderLine(lineIn, columnMapping, intensityByExperimentColumns);
                        if (!validHeader)
                            return;

                        headerParsed = true;
                        continue;
                    }

                    var splitLine = lineIn.Split('\t');

                    if (splitLine.Length < 30)
                        continue;

                    if (!GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.Id], out var peptideId, -1))
                    {
                        OnWarningEvent(string.Format(
                            "Line {0} in file {1} does not have an integer in the id column", lineNumber, inputFile.Name));
                        continue;
                    }

                    var peptideInfo = new MaxQuantPeptideInfo(peptideId);

                    GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.Sequence], out peptideInfo.Sequence);

                    if (string.IsNullOrWhiteSpace(peptideInfo.Sequence))
                    {
                        OnWarningEvent(string.Format(
                            "Line {0} in file {1} does not have a peptide in the Sequence column", lineNumber, inputFile.Name));
                        continue;
                    }

                    maxQuantPeptides.Add(peptideInfo.Sequence, peptideInfo);

                    GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.Prefix], out peptideInfo.Prefix);
                    GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.Suffix], out peptideInfo.Suffix);

                    if (GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.Proteins], out string proteinList))
                    {
                        foreach (var protein in proteinList.Split(';'))
                        {
                            peptideInfo.Proteins.Add(protein);
                        }
                    }

                    GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.LeadingRazorProtein], out peptideInfo.LeadingRazorProtein);
                    GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.TotalPeptideIntensity], out peptideInfo.TotalPeptideIntensity);

                    foreach (var item in intensityByExperimentColumns)
                    {
                        GetColumnValue(splitLine, item.Key, out string experimentIntensity);
                        peptideInfo.IntensityByExperiment.Add(item.Value, experimentIntensity);
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadPeptideInfo", ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
            }
        }

        private bool LookupMaxQuantPeptideInfo(
            IReadOnlyDictionary<string, MaxQuantPeptideInfo> maxQuantPeptides,
            string sequence,
            out MaxQuantPeptideInfo peptideInfo)
        {
            if (!maxQuantPeptides.TryGetValue(sequence, out peptideInfo))
            {
                OnWarningEvent("Peptide from msms.txt file was not present in the peptides.txt file: " + sequence);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parse the Synopsis file to create the other PHRP-compatible files
        /// </summary>
        /// <param name="inputFilePath"></param>
        /// <param name="outputDirectoryPath"></param>
        /// <param name="resetMassCorrectionTagsAndModificationDefinitions"></param>
        /// <param name="modList"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantSynopsisFile(
            string inputFilePath,
            string outputDirectoryPath,
            bool resetMassCorrectionTagsAndModificationDefinitions,
            IReadOnlyCollection<MSGFPlusParamFileModExtractor.ModInfo> modList)
        {
            var columnMapping = new Dictionary<MaxQuantSynFileColumns, int>();

            var staticModPresent = modList.Any(modItem => modItem.ModTypeInParameterFile == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod);

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (resetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                // Initialize searchResult
                var searchResult = new MaxQuantResults(mPeptideMods, mPeptideSeqMassCalculator);

                var errorMessages = new List<string>();

                try
                {
                    searchResult.UpdateSearchResultEnzymeAndTerminusInfo(Options);

                    // Open the input file and parse it
                    // Initialize the stream reader
                    using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                    var resultsProcessed = 0;
                    var headerParsed = false;

                    // Create the output files
                    var baseOutputFilePath = Path.Combine(outputDirectoryPath, Path.GetFileName(inputFilePath));
                    var filesInitialized = InitializeSequenceOutputFiles(baseOutputFilePath);
                    if (!filesInitialized)
                        return false;

                    // Parse the input file
                    while (!reader.EndOfStream && !AbortProcessing)
                    {
                        var lineIn = reader.ReadLine();

                        if (string.IsNullOrWhiteSpace(lineIn))
                        {
                            continue;
                        }

                        if (!headerParsed)
                        {
                            var validHeader = ParseMaxQuantSynFileHeaderLine(lineIn, columnMapping);
                            if (!validHeader)
                            {
                                // Error parsing header
                                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                                return false;
                            }
                            headerParsed = true;
                            continue;
                        }

                        var validSearchResult = ParseMaxQuantSynFileEntry(
                            lineIn, searchResult, errorMessages,
                            resultsProcessed, columnMapping, out _);

                        resultsProcessed++;
                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, true, modList, staticModPresent);
                        if (!modsAdded && errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                        {
                            errorMessages.Add(string.Format("Error adding modifications to sequence for ResultID '{0}'", searchResult.ResultID));
                        }

                        SaveResultsFileEntrySeqInfo(searchResult, true);

                        // Update the progress
                        UpdateSynopsisFileCreationProgress(reader);
                    }

                    if (Options.CreateModificationSummaryFile)
                    {
                        // Create the modification summary file
                        var inputFile = new FileInfo(inputFilePath);
                        var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));

                        if (string.IsNullOrWhiteSpace(modificationSummaryFilePath))
                        {
                            OnWarningEvent("ParseMaxQuantSynopsisFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
                        }
                        else
                        {
                            modificationSummaryFilePath = Path.Combine(outputDirectoryPath, modificationSummaryFilePath);

                            SaveModificationSummaryFile(modificationSummaryFilePath);
                        }
                    }

                    // Inform the user if any errors occurred
                    if (errorMessages.Count > 0)
                    {
                        SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                    }

                    return true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error reading input file in ParseMaxQuantSynopsisFile", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                    return false;
                }
                finally
                {
                    CloseSequenceOutputFiles();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error creating the output file in ParseMaxQuantSynopsisFile", ex);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parse a MaxQuant results line while creating the MaxQuant synopsis file
        /// </summary>
        /// <param name="maxQuantPeptides"></param>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="columnMapping"></param>
        /// <param name="modList"></param>
        /// <param name="lineNumber">Line number in the input file (used for error reporting)</param>
        /// <returns>True if successful, false if an error</returns>
        /// <remarks>
        /// This method is called while reading the msms.txt file
        /// For column details, see http://www.coxdocs.org/doku.php?id=maxquant:table:msmstable
        /// </remarks>
        private bool ParseMaxQuantResultsFileEntry(
            IReadOnlyDictionary<string, MaxQuantPeptideInfo> maxQuantPeptides,
            string lineIn,
            out MaxQuantSearchResult searchResult,
            ICollection<string> errorMessages,
            IDictionary<MaxQuantResultsFileColumns, int> columnMapping,
            IList<MSGFPlusParamFileModExtractor.ModInfo> modList,
            int lineNumber)
        {
            searchResult = new MaxQuantSearchResult();

            var staticModPresent = modList.Any(modItem => modItem.ModTypeInParameterFile == MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod);

            try
            {
                searchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 20)
                {
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.RawFile], out searchResult.DatasetName);

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Scan], out searchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid on line " + lineNumber, true);
                }

                if (!int.TryParse(searchResult.Scan, out searchResult.ScanNum))
                {
                    ReportError("Scan column is not numeric on line " + lineNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ScanIndex], out searchResult.ScanIndex);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Charge], out searchResult.Charge);
                searchResult.ChargeNum = Convert.ToInt16(StringUtilities.CIntSafe(searchResult.Charge, 0));

                // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MaxQuant
                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.CalculatedMonoMass], out searchResult.CalculatedMonoMass))
                {
                    double.TryParse(searchResult.CalculatedMonoMass, out searchResult.CalculatedMonoMassValue);
                }

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Sequence], out searchResult.Sequence))
                {
                    ReportError("Sequence column is missing or invalid on line " + lineNumber, true);
                }

                // Note that we replace NaN (not-a-number) with an empty string for several of the values

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Length], out searchResult.Length);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MissedCleavageCount], out searchResult.MissedCleavageCount);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Modifications], out searchResult.ModificationSummary);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ModifiedSequence], out searchResult.ModifiedSequence);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Proteins], out searchResult.Proteins);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Fragmentation], out searchResult.Fragmentation);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MassAnalyzer], out searchResult.MassAnalyzer);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorType], out searchResult.PrecursorType);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ScanEventNumber], out searchResult.ScanEventNumber);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.IsotopeIndex], out searchResult.IsotopeIndex);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MZ], out searchResult.PrecursorMZ);

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                searchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(searchResult.CalculatedMonoMassValue, 0), 6);

                GetColumnValueCheckNaN(splitLine, columnMapping[MaxQuantResultsFileColumns.MassErrorPPM], out searchResult.MassErrorPpmMaxQuant);
                GetColumnValueCheckNaN(splitLine, columnMapping[MaxQuantResultsFileColumns.MassErrorDa], out searchResult.MassErrorDaMaxQuant);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.SimpleMassErrorPPM], out searchResult.SimpleMassErrorPPM);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.RetentionTime], out searchResult.RetentionTime);

                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PEP], out searchResult.PEP))
                {
                    double.TryParse(searchResult.PEP, out searchResult.PEPValue);
                }

                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Score], out searchResult.Score))
                {
                    double.TryParse(searchResult.Score, out searchResult.ScoreValue);
                }

                GetColumnValueCheckNaN(splitLine, columnMapping[MaxQuantResultsFileColumns.DeltaScore], out searchResult.DeltaScore);
                GetColumnValueCheckNaN(splitLine, columnMapping[MaxQuantResultsFileColumns.ScoreDiff], out searchResult.ScoreDiff);
                GetColumnValueCheckNaN(splitLine, columnMapping[MaxQuantResultsFileColumns.LocalizationProb], out searchResult.LocalizationProb);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Combinatorics], out searchResult.Combinatorics);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PIF], out searchResult.PIF);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.FractionOfTotalSpectrum], out searchResult.FractionOfTotalSpectrum);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.BasePeakFraction], out searchResult.BasePeakFraction);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorScanNumber], out searchResult.PrecursorScanNumber);
                GetColumnValueCheckNaN(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorIntensity], out searchResult.PrecursorIntensity);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorApexFraction], out searchResult.PrecursorApexFraction);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorApexOffset], out searchResult.PrecursorApexOffset);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorApexOffsetTime], out searchResult.PrecursorApexOffsetTime);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.NumberOfMatches], out searchResult.NumberOfMatches);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.IntensityCoverage], out searchResult.IntensityCoverage);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PeakCoverage], out searchResult.PeakCoverage);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Reverse], out searchResult.Reverse);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ID], out searchResult.MsMsID);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ProteinGroupIDs], out searchResult.ProteinGroupIDs);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PeptideID], out searchResult.PeptideID);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ModPeptideID], out searchResult.ModPeptideID);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.EvidenceID], out searchResult.EvidenceID);

                // Parse the modification list to determine the total mod mass
                var totalModMass = ComputeTotalModMass(searchResult, modList, staticModPresent);

                // Construct the list of dynamic modifications
                searchResult.Modifications = ConstructDynamicModificationList(searchResult, modList);

                // Compute monoisotopic mass of the peptide
                searchResult.CalculatedMonoMassPHRP = ComputePeptideMass(searchResult.Sequence, totalModMass);

                // For other tools, we manually compute the mass error using the observed precursor ion m/z value
                // This cannot be done with the information in the msms.txt file, but it can be done if a PrecursorInfo file exists
                // Method ComputeObservedMassErrors will do this

                // Compare the mass computed by PHRP to the one reported by MaxQuant
                var deltaMassVsMaxQuant = searchResult.CalculatedMonoMassValue - searchResult.CalculatedMonoMassPHRP;
                if (Math.Abs(deltaMassVsMaxQuant) > 0.01)
                {
                    OnWarningEvent(string.Format(
                        "Calculated monoisotopic mass differs from the value reported by MaxQuant; delta mass: {0:F3} Da", deltaMassVsMaxQuant));
                }

                // Lookup the Prefix and Suffix residues
                var peptideInfoFound = LookupMaxQuantPeptideInfo(maxQuantPeptides, searchResult.Sequence, out var peptideInfo);
                if (peptideInfoFound)
                {
                    searchResult.PrefixResidue = peptideInfo.Prefix;
                    searchResult.SuffixResidue = peptideInfo.Suffix;
                    searchResult.LeadingRazorProtein = peptideInfo.LeadingRazorProtein;
                    searchResult.TotalPeptideIntensity = peptideInfo.TotalPeptideIntensity;

                    var cleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(GetPeptideSequence(searchResult));

                    searchResult.NumberOfTrypticTerminii = cleavageState switch
                    {
                        PeptideCleavageStateCalculator.PeptideCleavageState.Full => 2,
                        PeptideCleavageStateCalculator.PeptideCleavageState.Partial => 1,
                        PeptideCleavageStateCalculator.PeptideCleavageState.NonSpecific => 0,
                        PeptideCleavageStateCalculator.PeptideCleavageState.Unknown => 0,
                        _ => 0
                    };
                }
                else
                {
                    searchResult.PrefixResidue = "_";
                    searchResult.SuffixResidue = "_";
                }

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the MassMaxQuant results file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add(string.Format(
                        "Error parsing MaxQuant results in ParseMaxQuantResultsFileEntry, line {0}", lineNumber));
                }

                return false;
            }
        }

        /// <summary>
        /// Parse the header line in the MaxQuant peptides.txt file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <param name="intensityByExperimentColumns">Keys are column index, values are experiment name</param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantPeptidesFileHeaderLine(
            string lineIn,
            IDictionary<MaxQuantPeptidesFileColumns, int> columnMapping,
            IDictionary<int, string> intensityByExperimentColumns)
        {
            var columnNames = new SortedDictionary<string, MaxQuantPeptidesFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Sequence", MaxQuantPeptidesFileColumns.Sequence},
                {"Amino acid before", MaxQuantPeptidesFileColumns.Prefix},
                {"Amino acid after", MaxQuantPeptidesFileColumns.Suffix},
                {"Proteins", MaxQuantPeptidesFileColumns.Proteins},
                {"Leading razor protein", MaxQuantPeptidesFileColumns.LeadingRazorProtein},
                {"Intensity", MaxQuantPeptidesFileColumns.TotalPeptideIntensity},
                {"id", MaxQuantPeptidesFileColumns.Id}
            };

            columnMapping.Clear();
            intensityByExperimentColumns.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MaxQuantPeptidesFileColumns resultColumn in Enum.GetValues(typeof(MaxQuantPeptidesFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');

                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }

                var intensityColumnIndex = columnMapping[MaxQuantPeptidesFileColumns.TotalPeptideIntensity];
                if (intensityColumnIndex < 0)
                {
                    OnWarningEvent("Intensity column not found in the MaxQuant peptides.txt file");
                    return true;
                }

                // The columns after the "Intensity" column list peptide intensity, by experiment
                // If only a single dataset was searched, there will be only one column
                // For more info, see http://www.coxdocs.org/doku.php?id=maxquant:table:peptidetable

                for (var index = intensityColumnIndex + 1; index < splitLine.Length; index++)
                {
                    if (splitLine[index].StartsWith("Intensity L ") ||
                        splitLine[index].StartsWith("Intensity M ") ||
                        splitLine[index].StartsWith("Intensity H "))
                    {
                        // These represent intensity from a light, medium, or heavy label partner
                        // Ignore theme
                        continue;
                    }

                    if (!splitLine[index].StartsWith("Intensity "))
                    {
                        break;
                    }

                    var experimentName = splitLine[index].Substring("Intensity ".Length);
                    intensityByExperimentColumns.Add(index, experimentName);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the MaxQuant peptides file", ex);
                return false;
            }
        }

        /// <summary>
        /// Parse the MaxQuant results file header line (file msms.txt), populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if this is a valid header line, otherwise false (meaning it is a data line)</returns>
        private bool ParseMaxQuantResultsFileHeaderLine(string lineIn, IDictionary<MaxQuantResultsFileColumns, int> columnMapping)
        {
            // The expected column order from MaxQuant:
            //   Raw file	Dataset_ID	Scan number	Scan index	Sequence	Length	Missed cleavages	Modifications	Modified sequence	Oxidation (M) Probabilities	Oxidation (M) Score diffs	Acetyl (Protein N-term)	Oxidation (M)	Proteins	Charge	Fragmentation	Mass analyzer	Type	Scan event number	Isotope index	m/z	Mass	etc.

            var columnNames = new SortedDictionary<string, MaxQuantResultsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Raw file", MaxQuantResultsFileColumns.RawFile},
                {"Scan number", MaxQuantResultsFileColumns.Scan},
                {"Scan index", MaxQuantResultsFileColumns.ScanIndex},
                {"Sequence", MaxQuantResultsFileColumns.Sequence},
                {"Length", MaxQuantResultsFileColumns.Length},
                {"Missed cleavages", MaxQuantResultsFileColumns.MissedCleavageCount},
                {"Modifications", MaxQuantResultsFileColumns.Modifications},
                {"Modified sequence", MaxQuantResultsFileColumns.ModifiedSequence},
                {"Proteins", MaxQuantResultsFileColumns.Proteins},
                {"Charge", MaxQuantResultsFileColumns.Charge},
                {"Fragmentation", MaxQuantResultsFileColumns.Fragmentation},
                {"Mass analyzer", MaxQuantResultsFileColumns.MassAnalyzer},
                {"Type", MaxQuantResultsFileColumns.PrecursorType},
                {"Scan event number", MaxQuantResultsFileColumns.ScanEventNumber},
                {"Isotope index", MaxQuantResultsFileColumns.IsotopeIndex},
                {"m/z", MaxQuantResultsFileColumns.MZ},
                {"Mass", MaxQuantResultsFileColumns.CalculatedMonoMass},
                {"Mass error [ppm]", MaxQuantResultsFileColumns.MassErrorPPM},
                {"Mass error [Da]", MaxQuantResultsFileColumns.MassErrorDa},
                {"Simple mass error [ppm]", MaxQuantResultsFileColumns.SimpleMassErrorPPM},
                {"Retention time", MaxQuantResultsFileColumns.RetentionTime},
                {"PEP", MaxQuantResultsFileColumns.PEP},
                {"Score", MaxQuantResultsFileColumns.Score},
                {"Delta score", MaxQuantResultsFileColumns.DeltaScore},
                {"Score diff", MaxQuantResultsFileColumns.ScoreDiff},
                {"Localization prob", MaxQuantResultsFileColumns.LocalizationProb},
                {"Combinatorics", MaxQuantResultsFileColumns.Combinatorics},
                {"PIF", MaxQuantResultsFileColumns.PIF},
                {"Fraction of total spectrum", MaxQuantResultsFileColumns.FractionOfTotalSpectrum},
                {"Base peak fraction", MaxQuantResultsFileColumns.BasePeakFraction},
                {"Precursor full scan number", MaxQuantResultsFileColumns.PrecursorScanNumber},
                {"Precursor Intensity", MaxQuantResultsFileColumns.PrecursorIntensity},
                {"Precursor apex fraction", MaxQuantResultsFileColumns.PrecursorApexFraction},
                {"Precursor apex offset", MaxQuantResultsFileColumns.PrecursorApexOffset},
                {"Precursor apex offset time", MaxQuantResultsFileColumns.PrecursorApexOffsetTime},
                {"Number of matches", MaxQuantResultsFileColumns.NumberOfMatches},
                {"Intensity coverage", MaxQuantResultsFileColumns.IntensityCoverage},
                {"Peak coverage", MaxQuantResultsFileColumns.PeakCoverage},
                {"Reverse", MaxQuantResultsFileColumns.Reverse},
                {"id", MaxQuantResultsFileColumns.ID},
                {"Protein group IDs", MaxQuantResultsFileColumns.ProteinGroupIDs},
                {"Peptide ID", MaxQuantResultsFileColumns.PeptideID},
                {"Mod. peptide ID", MaxQuantResultsFileColumns.ModPeptideID},
                {"Evidence ID", MaxQuantResultsFileColumns.EvidenceID}
            };

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MaxQuantResultsFileColumns resultColumn in Enum.GetValues(typeof(MaxQuantResultsFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');

                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the MaxQuant results file", ex);
                return false;
            }
        }

        /// <summary>
        /// Parse the header line of a MaxQuant _syn.txt file, populating columnMapping
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="columnMapping"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantSynFileHeaderLine(string lineIn, IDictionary<MaxQuantSynFileColumns, int> columnMapping)
        {
            var columnNames = MaxQuantSynFileReader.GetColumnHeaderNamesAndIDs();

            columnMapping.Clear();

            try
            {
                // Initialize each entry in columnMapping to -1
                foreach (MaxQuantSynFileColumns resultColumn in Enum.GetValues(typeof(MaxQuantSynFileColumns)))
                {
                    columnMapping.Add(resultColumn, -1);
                }

                var splitLine = lineIn.Split('\t');

                for (var index = 0; index < splitLine.Length; index++)
                {
                    if (columnNames.TryGetValue(splitLine[index], out var resultFileColumn))
                    {
                        // Recognized column name; update columnMapping
                        columnMapping[resultFileColumn] = index;
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error parsing the header line in the MaxQuant synopsis file", ex);
                return false;
            }

            return true;
        }

        /// <summary>
        /// Parses an entry from a MaxQuant Synopsis file
        /// </summary>
        /// <param name="lineIn"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        /// <param name="resultsProcessed"></param>
        /// <param name="columnMapping"></param>
        /// <param name="peptideSequence"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool ParseMaxQuantSynFileEntry(
            string lineIn,
            MaxQuantResults searchResult,
            ICollection<string> errorMessages,
            int resultsProcessed,
            IDictionary<MaxQuantSynFileColumns, int> columnMapping,
            out string peptideSequence)
        {
            string[] splitLine = null;

            // Reset searchResult
            searchResult.Clear();
            peptideSequence = string.Empty;

            try
            {
                splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 15)
                {
                    return false;
                }

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ResultID], out string value))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading ResultID value from MaxQuant results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Dataset], out string dataset);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DatasetID], out int datasetId);

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Scan], out string scan);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Charge], out string charge);

                searchResult.DatasetName = dataset;
                searchResult.DatasetID = datasetId;

                searchResult.Scan = scan;
                searchResult.Charge = charge;

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Peptide], out peptideSequence))
                {
                    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                    {
                        errorMessages.Add(string.Format(
                            "Error reading Peptide sequence value from MaxQuant results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Proteins], out string proteinNames);

                if (!string.IsNullOrWhiteSpace(proteinNames))
                {
                    // Protein names in the synopsis file are a semicolon separated list
                    foreach (var proteinName in proteinNames.Split(';'))
                    {
                        if (string.IsNullOrWhiteSpace(proteinName))
                            continue;

                        searchResult.Proteins.Add(proteinName);
                    }

                    if (searchResult.Proteins.Count > 0)
                    {
                        searchResult.ProteinName = searchResult.Proteins[0];
                        searchResult.MultipleProteinCount = (searchResult.Proteins.Count - 1).ToString();
                    }
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PrecursorMZ], out string precursorMz);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MH], out string parentIonMH);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Mass], out string monoisotopicMass);

                searchResult.PrecursorMZ = precursorMz;
                searchResult.ParentIonMH = parentIonMH;
                searchResult.CalculatedMonoMass = monoisotopicMass;

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DelM], out string phrpComputedDelM);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DelM_PPM], out string phrpComputedDelMppm);

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DelM_MaxQuant], out string maxquantComputedDelM);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DelM_PPM_MaxQuant], out string maxquantComputedDelMppm);

                searchResult.PeptideDeltaMass = phrpComputedDelM;
                searchResult.PHRPComputedDelM = phrpComputedDelM;
                searchResult.PHRPComputedDelMPPM = phrpComputedDelMppm;

                searchResult.MaxQuantComputedDelM = maxquantComputedDelM;
                searchResult.MaxQuantComputedDelMPPM = maxquantComputedDelMppm;

                // Note that MaxQuant peptides don't actually have mod symbols; that information is tracked via searchResult.Modifications
                // Thus, .PeptideSequenceWithMods will not have any mod symbols

                // Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                var searchResultBase = (SearchResultsBaseClass)searchResult;

                ComputePseudoPeptideLocInProtein(searchResultBase);

                // Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                searchResult.ComputePeptideCleavageStateInProtein();

                // Read the remaining data values

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.FragMethod], out string fragMethod);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.SpecIndex], out string specIndex);

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DynamicModifications], out string dynamicModifications);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.LeadingRazorProtein], out string leadingRazorProtein);

                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.NTT], out string ntt);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PEP], out string posteriorErrorProbability);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Score], out string andromedaScore);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.DeltaScore], out string deltaScore);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.TotalPeptideIntensity], out string totalPeptideIntensity);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MassAnalyzer], out string massAnalyzer);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PrecursorType], out string precursorType);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.RetentionTime], out string retentionTime);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PrecursorScan], out string precursorScan);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PrecursorIntensity], out string precursorIntensity);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.NumberOfMatches], out string numberOfMatches);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.IntensityCoverage], out string intensityCoverage);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MissedCleavages], out string missedCleavages);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MsMsID], out string msMsID);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ProteinGroupIDs], out string proteinGroupIDs);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PeptideID], out string peptideID);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ModPeptideID], out string modPeptideID);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.EvidenceID], out string evidenceID);
                GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.QValue], out string qValue);

                // Store the data
                searchResult.FragMethod = fragMethod;
                searchResult.SpecIndex = specIndex;
                searchResult.Modifications = dynamicModifications;

                searchResult.LeadingRazorProtein = leadingRazorProtein;

                if (leadingRazorProtein.Length > 0 && searchResult.Proteins.Count == 0)
                {
                    // This is most likely a peptide from a decoy protein (whose name should start with REV__)
                    // Add it to the protein list
                    searchResult.Proteins.Add(leadingRazorProtein);
                    searchResult.ProteinName = leadingRazorProtein;
                }

                searchResult.NTT = ntt;
                searchResult.PEP = posteriorErrorProbability;
                searchResult.Score = andromedaScore;
                searchResult.DeltaScore = deltaScore;
                searchResult.TotalPeptideIntensity = totalPeptideIntensity;
                searchResult.MassAnalyzer = massAnalyzer;
                searchResult.PrecursorType = precursorType;
                searchResult.RetentionTime = retentionTime;
                searchResult.PrecursorScanNumber = precursorScan;
                searchResult.PrecursorIntensity = precursorIntensity;
                searchResult.NumberOfMatches = numberOfMatches;
                searchResult.IntensityCoverage = intensityCoverage;
                searchResult.MissedCleavageCount = missedCleavages;
                searchResult.MsMsID = msMsID;
                searchResult.ProteinGroupIDs = proteinGroupIDs;
                searchResult.PeptideID = peptideID;
                searchResult.ModPeptideID = modPeptideID;
                searchResult.EvidenceID = evidenceID;
                searchResult.QValue = qValue;

                if (string.IsNullOrWhiteSpace(searchResult.MultipleProteinCount))
                {
                    searchResult.MultipleProteinCount = "0";
                }

                return true;
            }
            catch (Exception)
            {
                // Error parsing this row from the synopsis or first hits file
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    if (splitLine?.Length > 0)
                    {
                        errorMessages.Add(string.Format(
                            "Error parsing MaxQuant results for RowIndex '{0}'", splitLine[0]));
                    }
                    else
                    {
                        errorMessages.Add("Error parsing MaxQuant Results in ParseMaxQuantSynFileEntry");
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">MaxQuant results file (msms.txt); alternatively, a directory with files msms.txt and peptides.txt</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if successful, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            var success = false;

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                SetErrorCode(PHRPErrorCode.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                    return false;
                }

                success = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!success)
                {
                    return false;
                }

                // Check whether inputFilePath is a directory path
                var candidateDirectory = new DirectoryInfo(inputFilePath);
                if (candidateDirectory.Exists)
                {
                    ResetProgress("Parsing " + candidateDirectory.FullName);
                    inputFilePath = Path.Combine(candidateDirectory.FullName, MSMS_FILE_NAME);
                }
                else
                {
                    ResetProgress("Parsing " + PathUtils.CompactPathString(inputFilePath, 100));
                }

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);
                    if (inputFile.Directory == null)
                    {
                        SetErrorMessage("Unable to determine the parent directory of the input file: " + inputFilePath);
                        SetErrorCode(PHRPErrorCode.InvalidInputFilePath);
                        return false;
                    }

                    List<MSGFPlusParamFileModExtractor.ModInfo> modList;
                    if (string.IsNullOrWhiteSpace(Options.SearchToolParameterFilePath))
                    {
                        OnWarningEvent(
                            "MaxQuant parameter file not defined; cannot generate the _ModSummary.txt and _ModDetails.txt files. " +
                            "To define, use /N at the command line or SearchToolParameterFilePath in a parameter file");

                        modList = new List<MSGFPlusParamFileModExtractor.ModInfo>();
                    }
                    else
                    {
                        var modificationDefinitionSuccess = LoadMaxQuantModificationDefinitions(inputFile.Directory);
                        if (!modificationDefinitionSuccess)
                            return false;

                        // Load the MaxQuant Parameter File so that we can determine the modification names
                        var modInfoExtracted = ExtractModInfoFromParamFile(Options.SearchToolParameterFilePath, out modList);
                        if (!modInfoExtracted)
                        {
                            return false;
                        }
                    }

                    LoadPeptideInfo(inputFile.Directory, out var maxQuantPeptides);

                    // Do not create a first-hits file for MaxQuant results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    success = CreateSynResultsFile(maxQuantPeptides, inputFilePath, outputDirectoryPath, modList, out var baseName, out var synOutputFilePath);

                    // Create the other PHRP-specific files
                    ResetProgress("Creating the PHRP files for " + Path.GetFileName(synOutputFilePath), true);

                    // Now parse the _syn.txt file that we just created to create the other PHRP files
                    success = ParseMaxQuantSynopsisFile(synOutputFilePath, outputDirectoryPath, false, modList);

                    if (success && Options.CreateProteinModsFile)
                    {
                        // Check for an empty synopsis file
                        if (!ValidateFileHasData(synOutputFilePath, "Synopsis file", out var errorMessage))
                        {
                            OnWarningEvent(errorMessage);
                        }
                        else
                        {
                            success = CreateProteinModsFileWork(baseName, inputFile, synOutputFilePath, outputDirectoryPath);
                        }
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in MaxQuantResultsProcessor.ProcessFile (2)", ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile (1)", ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        /// <summary>
        /// Sort filteredSearchResults and write to disk
        /// </summary>
        /// <param name="baseNameByDatasetName">Keys are dataset names, values are dataset ID (or 0 if undefined)</param>
        /// <param name="writer"></param>
        /// <param name="filteredSearchResults"></param>
        /// <param name="errorMessages"></param>
        private void SortAndWriteFilteredSearchResults(
            Dictionary<string, string> baseNameByDatasetName,
            TextWriter writer,
            List<MaxQuantSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            var datasetIDs = LookupDatasetIDs(baseNameByDatasetName.Keys.ToList());

            // Sort filteredSearchResults by descending Andromeda score, Scan, Peptide, and Razor Protein
            filteredSearchResults.Sort(new MaxQuantSearchResultsComparerScoreScanChargePeptide());

            // Compute FDR values, then assign QValues
            ComputeQValues(filteredSearchResults);

            var index = 1;
            foreach (var result in filteredSearchResults)
            {
                var baseDatasetName = baseNameByDatasetName[result.DatasetName];

                if (!datasetIDs.TryGetValue(result.DatasetName, out var datasetID))
                {
                    datasetID = 0;
                }

                WriteSearchResultToFile(index, baseDatasetName, datasetID, writer, result, errorMessages);
                index++;
            }
        }

        /// <summary>
        /// Compute FDR values, then assign QValues
        /// </summary>
        /// <param name="searchResults"></param>
        /// <remarks>Assumes the data is sorted by descending score using MaxQuantSearchResultsComparerScoreScanChargePeptide</remarks>
        private void ComputeQValues(IList<MaxQuantSearchResult> searchResults)
        {
            var forwardPeptideCount = 0;
            var reversePeptideCount = 0;

            for (var index = 0; index < searchResults.Count; index++)
            {
                // Reverse should be true for reverse-hit peptides
                // Alternative, the LeadingRazor protein field will start with REV__ for reverse-hit peptides

                if (searchResults[index].Reverse || searchResults[index].LeadingRazorProtein.StartsWith("REV__"))
                {
                    reversePeptideCount++;
                }
                else
                {
                    forwardPeptideCount++;
                }

                double fdr = 1;

                if (forwardPeptideCount > 0)
                {
                    fdr = reversePeptideCount / Convert.ToDouble(forwardPeptideCount);
                }

                var udtResult = searchResults[index];
                udtResult.FDR = fdr;

                searchResults[index] = udtResult;
            }

            // Now compute Q-Values
            // We step through the list, from the worst scoring result to the best result
            // The first Q-Value is the FDR of the final entry
            // The next Q-Value is the minimum of (QValue, CurrentFDR)

            var qValue = searchResults.Last().FDR;
            if (qValue > 1)
                qValue = 1;

            for (var index = searchResults.Count - 1; index >= 0; index += -1)
            {
                var udtResult = searchResults[index];

                qValue = Math.Min(qValue, udtResult.FDR);
                udtResult.QValue = qValue;

                searchResults[index] = udtResult;
            }
        }

        /// <summary>
        /// Stores the synopsis file matches for a single scan (typically there will only be one result for MaxQuant)
        /// </summary>
        /// <param name="searchResults">Search results</param>
        /// <param name="startIndex">Start index for data in this scan</param>
        /// <param name="endIndex">End index for data in this scan</param>
        /// <param name="filteredSearchResults">Output parameter: the actual filtered search results</param>
        private void StoreSynMatches(
            IList<MaxQuantSearchResult> searchResults,
            int startIndex,
            int endIndex,
            List<MaxQuantSearchResult> filteredSearchResults)
        {
            AssignRankByScore(searchResults, startIndex, endIndex);


            // If the MaxQuant results file has multiple datasets, there will often be multiple results with the same scan number

            // The calling procedure already sorted by scan, charge, and Score; no need to re-sort

            ExpandListIfRequired(filteredSearchResults, endIndex - startIndex + 1);

            // Now store the matches that pass the filters
            //  Either Andromeda Score > AndromedaScoreThreshold
            //  or     pep < PosteriorErrorProbabilityThreshold
            for (var index = startIndex; index <= endIndex; index++)
            {
                if (searchResults[index].ScoreValue >= Options.MaxQuantAndromedaScoreThreshold ||
                    searchResults[index].PEPValue < Options.MaxQuantPosteriorErrorProbabilityThreshold)
                {
                    filteredSearchResults.Add(searchResults[index]);
                }
            }
        }

        /// <summary>
        /// Write out the header line for synopsis / first hits files
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorMessages"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ICollection<string> errorMessages)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = MaxQuantSynFileReader.GetColumnHeaderNamesAndIDs();

                var headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();

                writer.WriteLine(StringUtilities.CollapseList(headerNames));
            }
            catch (Exception)
            {
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add("Error writing synopsis / first hits header");
                }
            }
        }

        /// <summary>
        /// Writes an entry to a synopsis or first hits file
        /// </summary>
        /// <param name="resultID"></param>
        /// <param name="baseDatasetName"></param>
        /// <param name="datasetID"></param>
        /// <param name="writer"></param>
        /// <param name="searchResult"></param>
        /// <param name="errorMessages"></param>
        private void WriteSearchResultToFile(
            int resultID,
            string baseDatasetName,
            int datasetID,
            TextWriter writer,
            MaxQuantSearchResult searchResult,
            ICollection<string> errorMessages)
        {
            try
            {
                var data = new List<string>
                {
                    resultID.ToString(),
                    baseDatasetName,
                    datasetID.ToString(),
                    searchResult.Scan,
                    searchResult.Fragmentation,
                    searchResult.ScanIndex,
                    searchResult.Charge,
                    searchResult.PrecursorMZ,
                    searchResult.MassErrorDa,
                    searchResult.MassErrorPpm,
                    searchResult.MassErrorDaMaxQuant,
                    searchResult.MassErrorPpmMaxQuant,
                    searchResult.MH,
                    searchResult.CalculatedMonoMass,
                    GetPeptideSequence(searchResult),
                    searchResult.Modifications,
                    searchResult.Proteins,
                    searchResult.LeadingRazorProtein,
                    searchResult.NumberOfTrypticTerminii.ToString(),
                    searchResult.PEP,
                    searchResult.Score,
                    searchResult.DeltaScore,
                    searchResult.RankScore.ToString(),
                    searchResult.TotalPeptideIntensity,
                    searchResult.MassAnalyzer,
                    searchResult.PrecursorType,
                    searchResult.RetentionTime,
                    searchResult.PrecursorScanNumber,
                    searchResult.PrecursorIntensity,
                    searchResult.NumberOfMatches,
                    searchResult.IntensityCoverage,
                    searchResult.MissedCleavageCount,
                    searchResult.MsMsID,
                    searchResult.ProteinGroupIDs,
                    searchResult.PeptideID,
                    searchResult.ModPeptideID,
                    searchResult.EvidenceID,
                    PRISM.StringUtilities.DblToString(searchResult.QValue, 5, 0.00005)
                };

                writer.WriteLine(StringUtilities.CollapseList(data));
            }
            catch (Exception)
            {
                if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                {
                    errorMessages.Add("Error writing synopsis / first hits record");
                }
            }
        }

        /// <summary>
        /// Override this method to display the name of each class
        /// </summary>
        public override string ToString()
        {
            return TOOL_NAME + " results processor";
        }

        private class MaxQuantSearchResultsComparerScanChargeScorePeptide : IComparer<MaxQuantSearchResult>
        {
            public int Compare(MaxQuantSearchResult x, MaxQuantSearchResult y)
            {
                // First sort on Scan
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                // Scan is the same, check charge
                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {
                    return -1;
                }

                // Charge is the same; check Score
                if (x.ScoreValue < y.ScoreValue)
                {
                    return 1;
                }

                if (x.ScoreValue > y.ScoreValue)
                {
                    return -1;
                }

                // Score is the same; check sequence
                var result = string.CompareOrdinal(x.Sequence, y.Sequence);
                if (result == 0)
                {
                    // Peptide is the same, check Protein
                    result = string.CompareOrdinal(x.LeadingRazorProtein, y.LeadingRazorProtein);
                }
                return result;
            }
        }


        private class MaxQuantSearchResultsComparerScoreScanChargePeptide : IComparer<MaxQuantSearchResult>
        {
            public int Compare(MaxQuantSearchResult x, MaxQuantSearchResult y)
            {
                if (x.ScoreValue < y.ScoreValue)
                {
                    return 1;
                }

                if (x.ScoreValue > y.ScoreValue)
                {
                    return -1;
                }

                // P-value is the same; check scan number
                if (x.ScanNum > y.ScanNum)
                {
                    return 1;
                }

                if (x.ScanNum < y.ScanNum)
                {
                    return -1;
                }

                // Scan is the same, check charge
                if (x.ChargeNum > y.ChargeNum)
                {
                    return 1;
                }

                if (x.ChargeNum < y.ChargeNum)
                {
                    return -1;
                }

                // Charge is the same; check peptide
                return string.CompareOrdinal(x.Sequence, y.Sequence);
            }
        }
    }
}
