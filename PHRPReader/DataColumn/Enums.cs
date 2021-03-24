using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PHRPReader.DataColumn
{
    public class Enums
    {
        /// <summary>
        /// These columns correspond to the Synopsis file created by MSPathFinderResultsProcessor
        /// </summary>
        public enum MSPathFinderSynFile
        {
            ResultID = 0,
            Scan = 1,
            Charge = 2,
            MostAbundantIsotopeMz = 3,
            Mass = 4,
            Sequence = 5,                // PrefixLetter.Sequence.SuffixLetter
            Modifications = 6,
            Composition = 7,
            Protein = 8,
            ProteinDesc = 9,
            ProteinLength = 10,
            ResidueStart = 11,
            ResidueEnd = 12,
            MatchedFragments = 13,
            SpecEValue = 14,             // Column added 2015-08-25
            EValue = 15,                 // Column added 2015-08-25
            QValue = 16,
            PepQValue = 17
        }


        /// <summary>
        /// These columns correspond to the Synopsis file created by MSGFPlusResultsProcessor
        /// </summary>
        public enum MSGFPlusSynFileColumns
        {
            ResultID = 0,
            Scan = 1,
            FragMethod = 2,
            SpecIndex = 3,
            Charge = 4,
            PrecursorMZ = 5,
            DelM = 6,                            // Precursor error, in Daltons; if the search used a tolerance less than 0.5 Da or less than 500 ppm, this value is computed from the DelMPPM value
            DelMPPM = 7,                         // Precursor error, in ppm; corrected for isotope selection errors
            MH = 8,                              // Theoretical monoisotopic peptide mass (computed by PHRP)
            Peptide = 9,                         // This is the sequence with prefix and suffix residues and also with modification symbols
            Protein = 10,                        // Protein Name (remove description)
            NTT = 11,                            // Number of tryptic termini
            DeNovoScore = 12,
            MSGFScore = 13,
            SpecProb_EValue = 14,
            RankSpecProb = 15,                   // Rank 1 means lowest SpecEValue, 2 means next higher score, etc. (ties get the same rank)
            PValue_EValue = 16,
            FDR_QValue = 17,                     // Only present if searched using -tda 1
            PepFDR_PepQValue = 18,               // Only present if searched using -tda 1
            EFDR = 19,                           // Only present if did not search using -tda 1
            // ReSharper disable UnusedMember.Global
            IMSScan = 20,                        // Only present for MSGFDB_IMS results
            IMSDriftTime = 21,                   // Only present for MSGFDB_IMS results
            // ReSharper restore UnusedMember.Global
            IsotopeError = 22
        }
    }
}
