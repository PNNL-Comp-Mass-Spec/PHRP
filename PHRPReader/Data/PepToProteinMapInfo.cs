using System;
using System.Collections.Generic;

namespace PHRPReader.Data
{
    /// <summary>
    /// Track location of a given peptide in one or more proteins
    /// </summary>
    public class PepToProteinMapInfo
    {
        /// <summary>
        /// Start and end residue locations in a protein
        /// </summary>
        public class ProteinLocationInfo
        {
            /// <summary>
            /// Start residue (first residue is 1)
            /// </summary>
            public int ResidueStart;

            /// <summary>
            /// End residue
            /// </summary>
            public int ResidueEnd;
        }

        /// <summary>
        /// Dictionary of protein names and residue start/end positions for a given peptide
        /// </summary>
        private readonly Dictionary<string, List<ProteinLocationInfo>> mProteinMapInfo;

        /// <summary>
        /// Number of proteins that contain a given peptide
        /// </summary>
        public int ProteinCount => ProteinMapInfo.Count;

        /// <summary>
        /// Dictionary of protein names and residue start/end positions for a given peptide
        /// </summary>
        public IReadOnlyDictionary<string, List<ProteinLocationInfo>> ProteinMapInfo => mProteinMapInfo;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName">Protein name</param>
        /// <param name="residueStart">Location that a peptide starts in the protein</param>
        /// <param name="residueEnd">Location that a peptide ends in the protein</param>
        public PepToProteinMapInfo(string proteinName, int residueStart, int residueEnd)
        {
            mProteinMapInfo = new Dictionary<string, List<ProteinLocationInfo>>(StringComparer.OrdinalIgnoreCase);

            AddProtein(proteinName, residueStart, residueEnd);
        }

        /// <summary>
        /// Add another peptide to protein mapping for a given peptide
        /// </summary>
        /// <remarks>If an entry already exists for a protein at a given start position, the end position will be updated</remarks>
        /// <param name="proteinName">Protein name</param>
        /// <param name="residueStart">Location that a peptide starts in the protein</param>
        /// <param name="residueEnd">Location that a peptide ends in the protein</param>
        public void AddProtein(string proteinName, int residueStart, int residueEnd)
        {
            if (ProteinMapInfo.TryGetValue(proteinName, out var locations))
            {
                // Protein mapping already exists; check residueStart
                foreach (var udtLoc in locations)
                {
                    if (udtLoc.ResidueStart == residueStart)
                    {
                        // Update this entry
                        udtLoc.ResidueEnd = residueEnd;
                        return;
                    }
                }

                var udtLocInfoAdditional = new ProteinLocationInfo
                {
                    ResidueStart = residueStart,
                    ResidueEnd = residueEnd
                };

                locations.Add(udtLocInfoAdditional);
                return;
            }

            // Protein not found
            var udtLocInfo = new ProteinLocationInfo
            {
                ResidueStart = residueStart,
                ResidueEnd = residueEnd
            };

            mProteinMapInfo.Add(proteinName, new List<ProteinLocationInfo> { udtLocInfo });
        }
    }
}
