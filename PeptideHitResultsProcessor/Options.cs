﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PeptideHitResultsProcessor.Processor;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class PHRPOptions
    {
        // Ignore Spelling: MaxQuant, MODa

        public bool CreateModificationSummaryFile { get; set; }

        public bool CreateFirstHitsFile { get; set; }

        public bool CreateSynopsisFile { get; set; }

        /// <summary>
        /// Create protein mods file
        /// </summary>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the Fasta file specified by mFastaFilePath</remarks>
        public bool CreateProteinModsFile { get; set; }

        public PeptideCleavageStateCalculator.EnzymeMatchSpecInfo EnzymeMatchSpec { get; set; }

        public string FastaFilePath { get; set; }

        public bool IgnorePeptideToProteinMapperErrors { get; set; }

        /// <summary>
        /// Inspect synopsis file p-value threshold
        /// </summary>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float InspectSynopsisFilePValueThreshold { get; set; }

        public string MassCorrectionTagsFilePath { get; set; }

        public string ModificationDefinitionsFilePath { get; set; }

        /// <summary>
        /// MaxQuant Andromeda score threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
        /// </remarks>
        public int MaxQuantAndromedaScoreThreshold { get; set; }

        /// <summary>
        /// MaxQuant PEP score threshold to use when creating the synopsis file
        /// </summary>
        /// <remarks>
        /// A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
        /// </remarks>
        public float MaxQuantPosteriorErrorProbabilityThreshold { get; set; }

        /// <summary>
        /// Used by MODaResultsProcessor and MODPlusResultsProcessor
        /// </summary>
        /// <remarks>Higher probability are higher confidence results</remarks>
        public float MODaMODPlusSynopsisFileProbabilityThreshold { get; set; }

        /// <summary>
        /// Used by MSAlign and TopPIC
        /// </summary>
        /// <remarks>Lower p-values are higher confidence results</remarks>
        public float MSAlignAndTopPICSynopsisFilePValueThreshold { get; set; }

        /// <summary>
        /// Used by MSGFPlusResultsProcessor and MSPathFinderResultsProcessor
        /// </summary>
        /// <remarks>Lower E-values are higher confidence results</remarks>
        public float MSGFPlusSynopsisFileEValueThreshold { get; set; }

        /// <summary>
        /// MSGFPlusResultsProcessor and MSPathFinderResultsProcessor
        /// </summary>
        /// <remarks>Lower SpecEValue values are higher confidence results</remarks>
        public float MSGFPlusSynopsisFileSpecEValueThreshold { get; set; }

        /// <summary>
        /// Typical non-zero value is 17.0027387
        /// </summary>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideCTerminusMassChange { get; set; }

        /// <summary>
        /// Typical non-zero value is 1.0078246
        /// </summary>
        /// <remarks>Ignored if equal to 0</remarks>
        public double PeptideNTerminusMassChange { get; set; }

        public bool ProteinModsFileIncludesReversedProteins { get; set; }

        /// <summary>
        /// Search tool parameter file path (aka SearchEngineParamFileName)
        /// </summary>
        /// <remarks>Used by MSGFPlusResultsProcessor, MaxQuantResultsProcessor, and others</remarks>
        public string SearchToolParameterFilePath { get; set; }

        /// <summary>
        /// Use Existing MTS PepToProtein Map File
        /// </summary>
        /// <remarks>If this is true and the _PepToProtMap.txt file isn't found, it will be created using the Fasta file specified by mFastaFilePath</remarks>
        public bool UseExistingMTSPepToProteinMapFile { get; set; }

        public bool WarnMissingParameterFileSection { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public PHRPOptions()
        {
            ResetToDefaults();
        }

        /// <summary>
        /// Reset options to defaults
        /// </summary>
        public void ResetToDefaults()
        {
            WarnMissingParameterFileSection = false;

            CreateModificationSummaryFile = true;

            CreateProteinModsFile = false;
            ProteinModsFileIncludesReversedProteins = false;
            UseExistingMTSPepToProteinMapFile = false;

            CreateFirstHitsFile = true;
            CreateSynopsisFile = true;

            FastaFilePath = string.Empty;
            IgnorePeptideToProteinMapperErrors = false;

            MassCorrectionTagsFilePath = string.Empty;
            ModificationDefinitionsFilePath = string.Empty;
            SearchToolParameterFilePath = string.Empty;

            InspectSynopsisFilePValueThreshold = InSpecTResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MaxQuantAndromedaScoreThreshold = MaxQuantResultsProcessor.DEFAULT_ANDROMEDA_SCORE_THRESHOLD;
            MaxQuantPosteriorErrorProbabilityThreshold = MaxQuantResultsProcessor.DEFAULT_PEP_THRESHOLD;

            MODaMODPlusSynopsisFileProbabilityThreshold = MODPlusResultsProcessor.DEFAULT_SYN_FILE_PROBABILITY_THRESHOLD;

            MSAlignAndTopPICSynopsisFilePValueThreshold = MSAlignResultsProcessor.DEFAULT_SYN_FILE_PVALUE_THRESHOLD;

            MSGFPlusSynopsisFileEValueThreshold = MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_EVALUE_THRESHOLD;
            MSGFPlusSynopsisFileSpecEValueThreshold = MSGFPlusResultsProcessor.DEFAULT_SYN_FILE_MSGF_SPEC_EVALUE_THRESHOLD;

            EnzymeMatchSpec = PeptideCleavageStateCalculator.GetDefaultEnzymeMatchSpec();

            PeptideNTerminusMassChange = PeptideMassCalculator.DEFAULT_N_TERMINUS_MASS_CHANGE;
            PeptideCTerminusMassChange = PeptideMassCalculator.DEFAULT_C_TERMINUS_MASS_CHANGE;
        }
    }
}