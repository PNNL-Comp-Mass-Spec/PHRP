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
using System.Text.RegularExpressions;
using System.Xml.Linq;
using PeptideHitResultsProcessor.Data;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;
using PRISM;
using static PHRPReader.MSGFPlusParamFileModExtractor;

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
    /// 4) ParseMaxQuantResultsFileEntry reads each data line and stores in an instance of MaxQuantSearchResult, which is a private struct
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

        // Ignore Spelling: AaSubstitution, acetyl, Carbamidomethyl, conf, Cterm, Crosslink, Da, Desc, diff, diffs, DimethNter
        // Ignore Spelling: Glu, Gln, Glycan, maxq, MaxQuant, NeuCode, Nterm, plex, pyro, struct, terminii

        // ReSharper restore CommentTypo

        /// <summary>
        /// Constructor
        /// </summary>
        public MaxQuantResultsProcessor(PHRPOptions options) : base(options)
        {
            FileDate = "April 20, 2021";

            mModCountMatcher = new Regex(@"^(?<ModCount>\d+) (?<ModName>.+)", RegexOptions.Compiled);

            MaxQuantMods = new Dictionary<string, MaxQuantModInfo>(StringComparer.OrdinalIgnoreCase);
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
            Intensity = 5,
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
            MissedCleavages = 5,
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
            public string Length;

            /// <summary>
            /// Number of missed enzymatic cleavages
            /// </summary>
            public string MissedCleavages;

            /// <summary>
            /// Post-translational modifications contained within the identified peptide
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
            public string Modifications;

            /// <summary>
            /// Peptide sequence with embedded modifications
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
            /// Charge state of the precursor ion.
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
            /// Precursor ion m/z
            /// </summary>
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
            /// vs. the predicted monoisotopic mass of the identified peptide sequence
            /// </summary>
            public string MassErrorPPM;

            /// <summary>
            /// Mass error, in Da
            /// </summary>
            public string MassErrorDa;

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

            /// <summary>
            /// PTM localization score
            /// </summary>
            public string LocalizationProb;

            /// <summary>
            /// Number of possible distributions of the modifications over the peptide sequence
            /// </summary>
            public string Combinatorics;

            /// <summary>
            ///Summed up extracted ion current (XIC) of all isotopic clusters associated with this peptide
            /// </summary>
            /// <remarks>
            /// From the peptides.txt file, column "Intensity"
            /// </remarks>
            public string PeptideIntensity;

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
                Length = string.Empty;
                MissedCleavages = string.Empty;
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
                MassErrorPPM = string.Empty;
                MassErrorDa = string.Empty;
                SimpleMassErrorPPM = string.Empty;
                RetentionTime = string.Empty;
                PEP = string.Empty;
                PEPValue = 0;
                Score = string.Empty;
                ScoreValue = 0;
                DeltaScore = string.Empty;
                ScoreDiff = string.Empty;
                LocalizationProb = string.Empty;
                Combinatorics = string.Empty;
                PeptideIntensity = string.Empty;
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
            }
        }

        private readonly Regex mModCountMatcher;

        private Dictionary<string, MaxQuantModInfo> MaxQuantMods { get; }

        /// <summary>
        /// Step through the Modifications and associate each modification with the residues
        /// For each residue, check if a static mod is defined that affects that residue
        /// For each mod mass, determine the modification and add to searchResult
        /// </summary>
        /// <param name="searchResult"></param>
        /// <param name="updateModOccurrenceCounts"></param>
        /// <param name="modList"></param>
        private void AddModificationsToResidues(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<ModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(searchResult.Modifications) || searchResult.Modifications.Equals("Unmodified"))
            {
                return;
            }

            var modsWithCounts = searchResult.Modifications.Split(',');
            var finalResidueLoc = searchResult.PeptideCleanSequence.Length;

            var modNames = new SortedSet<string>();

            foreach (var modItem in modsWithCounts)
            {
                // Check whether the mod name is preceded by an integer (which indicates the number of times this peptide has this modification)
                var match = mModCountMatcher.Match(modItem);

                var modName = match.Success ? match.Groups["ModName"].Value : modItem;

                if (!modNames.Contains(modName))
                    modNames.Add(modName);
            }

            // Keys are mod names, values are index in searchResult.ModifiedSequence
            var modPositions = new Dictionary<string, int>();

            while (true)
            {
                // Look for the earliest occurrence of any of the mods in searchResult.ModifiedSequence

                foreach (var mod in modNames)
                {
                    var charIndex = searchResult.ModifiedSequence.IndexOf(string.Format("({0})", mod), StringComparison.Ordinal);
                    modPositions.Add(mod, charIndex);
                }

                var minIndex = int.MaxValue;
                var firstModName = string.Empty;

                foreach (var modName in modPositions.Keys)
                {
                    if (modPositions[modName] >= 0 && modPositions[modName] < minIndex)
                    {
                        minIndex = modPositions[modName];
                        firstModName = modName;
                    }
                }

                if (string.IsNullOrEmpty(firstModName))
                {
                    // No more matches
                    break;
                }

                // Examine this modification's position to see if it is an N-terminal mod
                var nTerminalMod = MaxQuantMods.TryGetValue(firstModName, out var modInfo) &&
                                   modInfo.Position is MaxQuantModPosition.AnyNterm or MaxQuantModPosition.ProteinNterm;

                if (nTerminalMod)
                {
                    // The mod applies to the residue just after the mod name (which is surrounded by parentheses)
                }
                else
                {
                    // The mod applies to the residue just before the mod name (which is surrounded by parentheses)
                }

                // ToDo: Associate the mod with the given residue
                // searchResult.SearchResultAddModification(
                //     modDef.ModMassVal, chMostRecentResidue, residueLocInPeptide,
                //     residueTerminusState, updateModOccurrenceCounts);

            }
        }

        private bool AddModificationsAndComputeMass(
            MaxQuantResults searchResult,
            bool updateModOccurrenceCounts,
            IReadOnlyCollection<ModInfo> modList)
        {
            bool success;

            try
            {
                // Some of the other tools add IsotopicMods here
                // This is not supported for MaxQuant
                //
                // searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts)

                // Parse .Modifications to determine the modified residues present
                AddModificationsToResidues(searchResult, updateModOccurrenceCounts, modList);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Populate .PeptideModDescription
                searchResult.UpdateModDescription();

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
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
        /// <param name="modificationList">Comma separated list of modifications, e.g. "Acetyl (Protein N-term),Oxidation (M)" or "2 Oxidation (M) "</param>
        /// <param name="modList"></param>
        private double ComputeTotalModMass(
            string modificationList,
            IReadOnlyCollection<ModInfo> modList)
        {
            if (string.IsNullOrWhiteSpace(modificationList))
            {
                return 0;
            }

            double totalModMass = 0;
            var mods = modificationList.Split(',');

            foreach (var modItem in mods)
            {
                // Convert the mod name to a mass value
                var matchFound = false;

                // Obtain the mod name, for example "Oxidation (M)" from "2 Oxidation (M)"

                var match = mModCountMatcher.Match(modItem);

                var modName = match.Success ? match.Groups["ModName"].Value : modItem;

                foreach (var modDef in modList)
                {
                    if (string.Equals(modDef.ModName, modName, StringComparison.OrdinalIgnoreCase))
                    {
                        totalModMass += modDef.ModMassVal;
                        matchFound = true;
                        break;
                    }
                }

                if (!matchFound)
                {
                    ReportError("Mod name " + modName + " was not defined in the MaxQuant parameter file; cannot determine the mod mass", true);
                }
            }

            return totalModMass;
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
        /// <param name="inputFilePath"></param>
        /// <param name="outputFilePath"></param>
        /// <param name="modList"></param>
        /// <returns>True if successful, false if an error</returns>
        private bool CreateSynResultsFile(
            IReadOnlyDictionary<string, MaxQuantPeptideInfo> maxQuantPeptides,
            string inputFilePath,
            string outputFilePath,
            IReadOnlyCollection<ModInfo> modList)
        {
            try
            {
                var columnMapping = new Dictionary<MaxQuantResultsFileColumns, int>();
                var errorMessages = new List<string>();

                OnStatusEvent("Reading MaxQuant results file, " + inputFilePath);

                // Open the input file and parse it
                // Initialize the stream reader and the stream writer
                using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));
                using var writer = new StreamWriter(new FileStream(outputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

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
                                SetErrorMessage("Invalid header line in " + Path.GetFileName(inputFilePath));
                            }

                            return false;
                        }

                        // Write the header line to the output file
                        WriteSynFHTFileHeader(writer, errorMessages);

                        headerParsed = true;
                        continue;
                    }

                    var validSearchResult =
                        ParseMaxQuantResultsFileEntry(maxQuantPeptides, lineIn, out var udtSearchResult, errorMessages, columnMapping, modList, lineNumber);

                    if (validSearchResult)
                    {
                        searchResultsUnfiltered.Add(udtSearchResult);
                    }

                    // Update the progress
                    UpdateSynopsisFileCreationProgress(reader);
                }

                // Sort the SearchResults by scan, charge, and descending Andromeda score
                searchResultsUnfiltered.Sort(new MaxQuantSearchResultsComparerScanChargeScorePeptide());

                // Now filter the data

                // Initialize variables
                var startIndex = 0;

                while (startIndex < searchResultsUnfiltered.Count)
                {
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

                // Sort the data in filteredSearchResults then write out to disk
                SortAndWriteFilteredSearchResults(writer, filteredSearchResults, errorMessages);

                // Inform the user if any errors occurred
                if (errorMessages.Count > 0)
                {
                    SetErrorMessage("Invalid Lines: \n" + string.Join("\n", errorMessages));
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in CreateSynResultsFile: " + ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Read mod info from the MaxQuant parameter file
        /// </summary>
        /// <param name="maxQuantParamFilePath"></param>
        /// <param name="modList"></param>
        /// <returns>True on success, false if an error</returns>
        private bool ExtractModInfoFromParamFile(
            string maxQuantParamFilePath,
            out List<ModInfo> modList)
        {
            modList = new List<ModInfo>();

            try
            {
                var sourceFile = new FileInfo(maxQuantParamFilePath);
                if (!sourceFile.Exists)
                {
                    SetErrorMessage("MaxQuant parameter file not found: " + maxQuantParamFilePath);
                    SetErrorCode(PHRPErrorCode.ParameterFileNotFound);
                    return false;
                }

                // Open the MaxQuant parameter file with an XML reader and look for static and dynamic mod names
                OnStatusEvent("Reading the MaxQuant parameter file, " + sourceFile.Name);

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

                foreach (var fixedMod in fixedModNodes)
                {
                    var maxQuantModName = fixedMod.Value;
                    var modDef = GetModDetails(MSGFPlusModType.StaticMod, maxQuantModName);
                    modList.Add(modDef);
                }

                foreach (var dynamicMod in dynamicModNodes)
                {
                    var maxQuantModName = dynamicMod.Value;
                    var modDef = GetModDetails(MSGFPlusModType.DynamicMod, maxQuantModName);
                    modList.Add(modDef);
                }

                return true;
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ExtractModInfoFromParamFile:  " + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                return false;
            }
        }

        private ModInfo GetModDetails(MSGFPlusModType modType, string maxQuantModName)
        {
            var modDef = new ModInfo
            {
                ModType = modType
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
            modDef.ModSymbol = UNKNOWN_MSGFPlus_MOD_SYMBOL;

            switch (modInfo.Position)
            {
                case MaxQuantModPosition.Anywhere:
                case MaxQuantModPosition.NotCterm:
                    // Leave .ModType unchanged; this is a static or dynamic mod
                    break;

                case MaxQuantModPosition.AnyNterm:
                    if (modDef.ModType == MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusModType.DynamicMod;
                    }
                    modDef.Residues = AminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusModType.DynNTermPeptide;

                    break;

                case MaxQuantModPosition.AnyCterm:
                    if (modDef.ModType == MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusModType.DynamicMod;
                    }
                    modDef.Residues = AminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusModType.DynCTermPeptide;

                    break;

                case MaxQuantModPosition.ProteinNterm:
                    if (modDef.ModType == MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusModType.DynamicMod;
                    }
                    modDef.Residues = AminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusModType.DynNTermProtein;

                    break;

                case MaxQuantModPosition.ProteinCterm:
                    if (modDef.ModType == MSGFPlusModType.StaticMod && modDef.Residues != "-")
                    {
                        // This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                        modDef.ModType = MSGFPlusModType.DynamicMod;
                    }
                    modDef.Residues = AminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                    if (modDef.ModType == MSGFPlusModType.DynamicMod)
                        modDef.ModType = MSGFPlusModType.DynCTermProtein;

                    break;

                default:
                    throw new ArgumentOutOfRangeException();
            }

            return modDef;
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
                    if (!TryGetAttribute(modificationNode, "title", out var modTitle))
                    {
                        OnWarningEvent("Modification node in the MaxQuant modifications file is missing the title attribute");
                        continue;
                    }

                    if (!TryGetAttribute(modificationNode, "composition", out var composition))
                    {
                        OnWarningEvent("Modification node in the MaxQuant modifications file is missing the composition attribute");
                        continue;
                    }

                    TryGetAttribute(modificationNode, "description", out var modDescription);

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

                    if (TryGetElementValue(modificationNode, "position", out var positionText))
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

                    if (TryGetElementValue(modificationNode, "type", out var modTypeText))
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

                    if (TryGetElementValue(modificationNode, "terminus_type", out var terminusTypeText))
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
                        if (TryGetAttribute(modificationSiteNode, "site", out var modificationSite))
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
                SetErrorMessage("Error in LoadMaxQuantModificationDefinitions: " + ex.Message);
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

                OnStatusEvent("Reading peptides from file " + inputFile.FullName);

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
                    GetColumnValue(splitLine, columnMapping[MaxQuantPeptidesFileColumns.Intensity], out peptideInfo.Intensity);

                    foreach (var item in intensityByExperimentColumns)
                    {
                        GetColumnValue(splitLine, item.Key, out string experimentIntensity);
                        peptideInfo.IntensityByExperiment.Add(item.Value, experimentIntensity);
                    }
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in LoadPeptideInfo: " + ex.Message);
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
            IReadOnlyCollection<ModInfo> modList)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

            // Note that MaxQuant synopsis files are normally sorted on Probability value, ascending
            // In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
            //  we will keep track of the scan, charge, and peptide information parsed for each unique Probability encountered

            var columnMapping = new Dictionary<MaxQuantSynFileColumns, int>();

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

                // Initialize two SortedSets
                var peptidesFoundForSpecEValue = new SortedSet<string>();
                var peptidesFoundForQValue = new SortedSet<string>();

                var firstMatchForGroup = false;

                var previousSpecEValue = string.Empty;
                var previousQValue = string.Empty;

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

                        var validSearchResult = ParseMaxQuantSynFileEntry(lineIn, searchResult, errorMessages,
                            resultsProcessed, columnMapping,
                            out _);

                        resultsProcessed++;
                        if (!validSearchResult)
                        {
                            continue;
                        }

                        var key = searchResult.PeptideSequenceWithMods + "_" + searchResult.Scan + "_" + searchResult.Charge;

                        var newValue = true;

                        // ToDo:
                        //if (string.IsNullOrEmpty(searchResult.SpecEValue))
                        //{
                        //    if (searchResult.QValue == previousQValue)
                        //    {
                        //        // New result has the same QValue as the previous result
                        //        // See if peptidesFoundForQValue contains the peptide, scan and charge

                        //        if (peptidesFoundForQValue.Contains(key))
                        //        {
                        //            firstMatchForGroup = false;
                        //        }
                        //        else
                        //        {
                        //            peptidesFoundForQValue.Add(key);
                        //            firstMatchForGroup = true;
                        //        }

                        //        newValue = false;
                        //    }
                        //}
                        //else if (searchResult.SpecEValue == previousSpecEValue)
                        //{
                        //    // New result has the same SpecEValue as the previous result
                        //    // See if peptidesFoundForSpecEValue contains the peptide, scan and charge

                        //    if (peptidesFoundForSpecEValue.Contains(key))
                        //    {
                        //        firstMatchForGroup = false;
                        //    }
                        //    else
                        //    {
                        //        peptidesFoundForSpecEValue.Add(key);
                        //        firstMatchForGroup = true;
                        //    }

                        //    newValue = false;
                        //}

                        if (newValue)
                        {
                            // New SpecEValue or new QValue
                            // Reset the SortedSets
                            peptidesFoundForSpecEValue.Clear();
                            peptidesFoundForQValue.Clear();

                            // ToDo: Customize

                            // Update the cached values
                            //previousSpecEValue = searchResult.SpecEValue;
                            //previousQValue = searchResult.QValue;

                            // Append a new entry to the SortedSets
                            peptidesFoundForSpecEValue.Add(key);
                            peptidesFoundForQValue.Add(key);

                            firstMatchForGroup = true;
                        }

                        var modsAdded = AddModificationsAndComputeMass(searchResult, firstMatchForGroup, modList);
                        if (!modsAdded)
                        {
                            if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                            {
                                errorMessages.Add(string.Format("Error adding modifications to sequence for ResultID '{0}'", searchResult.ResultID));
                            }
                        }

                        SaveResultsFileEntrySeqInfo(searchResult, firstMatchForGroup);

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
                    SetErrorMessage("Error reading input file in ParseMaxQuantSynopsisFile: " + ex.Message);
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
                SetErrorMessage("Error creating the output file in ParseMaxQuantSynopsisFile: " + ex.Message);
                SetErrorCode(PHRPErrorCode.ErrorCreatingOutputFiles);
                return false;
            }
        }

        /// <summary>
        /// Parse a MaxQuant results line while creating the MaxQuant synopsis file
        /// </summary>
        /// <param name="maxQuantPeptides"></param>
        /// <param name="lineIn"></param>
        /// <param name="udtSearchResult"></param>
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
            out MaxQuantSearchResult udtSearchResult,
            ICollection<string> errorMessages,
            IDictionary<MaxQuantResultsFileColumns, int> columnMapping,
            IReadOnlyCollection<ModInfo> modList,
            int lineNumber)
        {
            udtSearchResult = new MaxQuantSearchResult();

            try
            {
                udtSearchResult.Clear();
                var splitLine = lineIn.TrimEnd().Split('\t');

                if (splitLine.Length < 20)
                {
                    return false;
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.RawFile], out udtSearchResult.DatasetName);

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Scan], out udtSearchResult.Scan))
                {
                    ReportError("Scan column is missing or invalid on line " + lineNumber, true);
                }

                if (!int.TryParse(udtSearchResult.Scan, out udtSearchResult.ScanNum))
                {
                    ReportError("Scan column is not numeric on line " + lineNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ScanIndex], out udtSearchResult.ScanIndex);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Charge], out udtSearchResult.Charge);
                udtSearchResult.ChargeNum = Convert.ToInt16(CIntSafe(udtSearchResult.Charge, 0));

                // Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MaxQuant
                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.CalculatedMonoMass], out udtSearchResult.CalculatedMonoMass))
                {
                    double.TryParse(udtSearchResult.CalculatedMonoMass, out udtSearchResult.CalculatedMonoMassValue);
                }

                // double.TryParse(udtSearchResult.CalculatedMonoMass, out var sequenceMonoMassMaxQuant);

                if (!GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Sequence], out udtSearchResult.Sequence))
                {
                    ReportError("Sequence column is missing or invalid on line " + lineNumber, true);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Length], out udtSearchResult.Length);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MissedCleavages], out udtSearchResult.MissedCleavages);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Modifications], out udtSearchResult.Modifications);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ModifiedSequence], out udtSearchResult.ModifiedSequence);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Proteins], out string proteins);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Fragmentation], out udtSearchResult.Fragmentation);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MassAnalyzer], out udtSearchResult.MassAnalyzer);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorType], out udtSearchResult.PrecursorType);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ScanEventNumber], out udtSearchResult.ScanEventNumber);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.IsotopeIndex], out udtSearchResult.IsotopeIndex);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MZ], out udtSearchResult.PrecursorMZ);

                // Store the monoisotopic MH value in .MH
                // This is (M+H)+ when the charge carrier is a proton
                udtSearchResult.MH = PRISM.StringUtilities.DblToString(mPeptideSeqMassCalculator.ConvoluteMass(udtSearchResult.CalculatedMonoMassValue, 0), 6);

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MassErrorPPM], out udtSearchResult.MassErrorPPM);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.MassErrorDa], out udtSearchResult.MassErrorDa);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.SimpleMassErrorPPM], out udtSearchResult.SimpleMassErrorPPM);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.RetentionTime], out udtSearchResult.RetentionTime);

                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PEP], out udtSearchResult.PEP))
                {
                    double.TryParse(udtSearchResult.PEP, out udtSearchResult.PEPValue);
                }

                if (GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Score], out udtSearchResult.Score))
                {
                    double.TryParse(udtSearchResult.Score, out udtSearchResult.ScoreValue);
                }

                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.DeltaScore], out udtSearchResult.DeltaScore);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ScoreDiff], out udtSearchResult.ScoreDiff);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.LocalizationProb], out udtSearchResult.LocalizationProb);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Combinatorics], out udtSearchResult.Combinatorics);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PIF], out udtSearchResult.PIF);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.FractionOfTotalSpectrum], out udtSearchResult.FractionOfTotalSpectrum);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.BasePeakFraction], out udtSearchResult.BasePeakFraction);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorScanNumber], out udtSearchResult.PrecursorScanNumber);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorIntensity], out udtSearchResult.PrecursorIntensity);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorApexFraction], out udtSearchResult.PrecursorApexFraction);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorApexOffset], out udtSearchResult.PrecursorApexOffset);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PrecursorApexOffsetTime], out udtSearchResult.PrecursorApexOffsetTime);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.NumberOfMatches], out udtSearchResult.NumberOfMatches);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.IntensityCoverage], out udtSearchResult.IntensityCoverage);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PeakCoverage], out udtSearchResult.PeakCoverage);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.Reverse], out udtSearchResult.Reverse);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ID], out udtSearchResult.MsMsID);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ProteinGroupIDs], out udtSearchResult.ProteinGroupIDs);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.PeptideID], out udtSearchResult.PeptideID);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.ModPeptideID], out udtSearchResult.ModPeptideID);
                GetColumnValue(splitLine, columnMapping[MaxQuantResultsFileColumns.EvidenceID], out udtSearchResult.EvidenceID);

                // Parse the modification list to determine the total mod mass
                var totalModMass = ComputeTotalModMass(udtSearchResult.Modifications, modList);

                // Compute monoisotopic mass of the peptide
                udtSearchResult.CalculatedMonoMassPHRP = ComputePeptideMass(udtSearchResult.Sequence, totalModMass);

                // Lookup the Prefix and Suffix residues
                var peptideInfoFound = LookupMaxQuantPeptideInfo(maxQuantPeptides, udtSearchResult.Sequence, out var peptideInfo);
                if (peptideInfoFound)
                {
                    udtSearchResult.PrefixResidue = peptideInfo.Prefix;
                    udtSearchResult.SuffixResidue = peptideInfo.Suffix;
                    udtSearchResult.LeadingRazorProtein = peptideInfo.LeadingRazorProtein;
                    udtSearchResult.PeptideIntensity = peptideInfo.Intensity;
                }
                else
                {
                    udtSearchResult.PrefixResidue = "_";
                    udtSearchResult.SuffixResidue = "_";
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
                {"Intensity", MaxQuantPeptidesFileColumns.Intensity},
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

                var intensityColumnIndex = columnMapping[MaxQuantPeptidesFileColumns.Intensity];
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
                SetErrorMessage("Error parsing the header line in the MaxQuant peptides file: " + ex.Message);
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
                {"Missed cleavages", MaxQuantResultsFileColumns.MissedCleavages},
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
                SetErrorMessage("Error parsing the header line in the MaxQuant results file: " + ex.Message);
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
                SetErrorMessage("Error parsing the header line in the MaxQuant synopsis file: " + ex.Message);
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
                            "Error reading ResultID value fromMaxQuant results, line {0}", resultsProcessed + 1));
                    }
                    return false;
                }

                searchResult.ResultID = int.Parse(value);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Scan], out string scan);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Charge], out string charge);

                //searchResult.Scan = scan;
                //searchResult.Charge = charge;

                //if (!GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Sequence], out peptideSequence))
                //{
                //    if (errorMessages.Count < MAX_ERROR_MESSAGE_COUNT)
                //    {
                //        errorMessages.Add(string.Format(
                //            "Error reading Peptide sequence value from MaxQuant results, line {0}", resultsProcessed + 1));
                //    }
                //    return false;
                //}

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Protein], out string proteinName);
                //searchResult.ProteinName = proteinName;
                //searchResult.MultipleProteinCount = "0";

                //searchResult.PeptideDeltaMass = "0";

                //// Note that MaxQuant sequences don't actually have mod symbols; that information is tracked via searchResult.Modifications

                //// Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                //searchResult.SetPeptideSequenceWithMods(peptideSequence, true, true);

                //var searchResultBase = (SearchResultsBaseClass)searchResult;

                //ComputePseudoPeptideLocInProtein(searchResultBase);

                //// Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                //searchResult.ComputePeptideCleavageStateInProtein();

                //// Read the remaining data values
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MostAbundantIsotopeMz], out string mostAbundantIsotopeMz);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Modifications], out string modifications);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.Composition], out string composition);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ProteinDesc], out string proteinDesc);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ProteinLength], out string proteinLength);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ResidueStart], out string residueStart);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.ResidueEnd], out string residueEnd);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.MatchedFragments], out string matchedFragments);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.SpecEValue], out string specEValue);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.EValue], out string eValue);

                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.QValue], out string qValue);
                //GetColumnValue(splitLine, columnMapping[MaxQuantSynFileColumns.PepQValue], out string pepQValue);

                //searchResult.MostAbundantIsotopeMz = mostAbundantIsotopeMz;
                //searchResult.Modifications = modifications;
                //searchResult.Composition = composition;
                //searchResult.ProteinDesc = proteinDesc;
                //searchResult.ProteinLength = proteinLength;
                //searchResult.ResidueStart = residueStart;
                //searchResult.ResidueEnd = residueEnd;
                //searchResult.MatchedFragments = matchedFragments;
                //searchResult.SpecEValue = specEValue;
                //searchResult.EValue = eValue;
                //searchResult.QValue = qValue;
                //searchResult.PepQValue = pepQValue;

                //// Update the peptide location in the protein
                //if (!string.IsNullOrEmpty(searchResult.ResidueStart))
                //{
                //    int.TryParse(searchResult.ResidueStart, out var peptideLocInProteinStart);
                //    searchResult.PeptideLocInProteinStart = peptideLocInProteinStart;
                //}

                //if (!string.IsNullOrEmpty(searchResult.ResidueEnd))
                //{
                //    int.TryParse(searchResult.ResidueEnd, out var peptideLocInProteinEnd);
                //    searchResult.PeptideLocInProteinEnd = peptideLocInProteinEnd;
                //}

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

                    List<ModInfo> modList;
                    if (string.IsNullOrWhiteSpace(Options.SearchToolParameterFilePath))
                    {
                        OnWarningEvent(
                            "MaxQuant parameter file not defined; cannot add modification symbols to peptides in the synopsis file. " +
                            "To define, use /N at the command line or SearchToolParameterFilePath in a parameter file");

                        modList = new List<ModInfo>();
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

                    var baseName = "Dataset_maxq";

                    // Do not create a first-hits file for MaxQuant results

                    // Create the synopsis output file
                    ResetProgress("Creating the SYN file", true);

                    // The synopsis file name will be of the form Dataset_maxq_syn.txt
                    var synOutputFilePath = Path.Combine(outputDirectoryPath, baseName + SEQUEST_SYNOPSIS_FILE_SUFFIX);

                    success = CreateSynResultsFile(maxQuantPeptides, inputFilePath, synOutputFilePath, modList);

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
                    SetErrorMessage("Error in MaxQuantResultsProcessor.ProcessFile (2):  " + ex.Message, ex);
                    SetErrorCode(PHRPErrorCode.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in ProcessFile (1):" + ex.Message, ex);
                SetErrorCode(PHRPErrorCode.UnspecifiedError);
            }

            return success;
        }

        private void SortAndWriteFilteredSearchResults(
            TextWriter writer,
            IEnumerable<MaxQuantSearchResult> filteredSearchResults,
            ICollection<string> errorMessages)
        {
            // Sort filteredSearchResults by descending Andromeda score, Scan, Peptide, and Razor Protein
            var query = from item in filteredSearchResults orderby item.ScoreValue descending, item.ScanNum, item.Sequence, item.LeadingRazorProtein select item;

            var index = 1;
            foreach (var result in query)
            {
                WriteSearchResultToFile(index, writer, result, errorMessages);
                index++;
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
            // If there was more than one result, we could rank them by score
            // AssignRankAndDeltaNormValues(searchResults, startIndex, endIndex)

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

                writer.WriteLine(CollapseList(headerNames));
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
        /// <param name="writer"></param>
        /// <param name="udtSearchResult"></param>
        /// <param name="errorMessages"></param>
        private void WriteSearchResultToFile(
            int resultID,
            TextWriter writer,
            MaxQuantSearchResult udtSearchResult,
            ICollection<string> errorMessages)
        {
            try
            {
                var data = new List<string>
                {
                    resultID.ToString(),
                    udtSearchResult.Scan,
                    udtSearchResult.Fragmentation,
                    udtSearchResult.ScanIndex,
                    udtSearchResult.Charge,
                    udtSearchResult.PrecursorMZ,
                    udtSearchResult.MassErrorDa,
                    udtSearchResult.MassErrorPPM,
                    udtSearchResult.MH,
                    udtSearchResult.PrefixResidue + "." + udtSearchResult.Sequence + udtSearchResult.SuffixResidue,
                    udtSearchResult.Proteins,
                    udtSearchResult.NumberOfTrypticTerminii.ToString(),
                    udtSearchResult.PEP,
                    udtSearchResult.Score,
                    udtSearchResult.DeltaScore,
                    udtSearchResult.PrecursorIntensity,
                    udtSearchResult.PeptideIntensity,
                    udtSearchResult.LeadingRazorProtein,
                    udtSearchResult.ProteinGroupIDs,
                    udtSearchResult.PeptideID,
                    udtSearchResult.ModPeptideID,
                    udtSearchResult.EvidenceID
                };

                writer.WriteLine(CollapseList(data));
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

                // Scan is the same; check Score
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
    }
}
