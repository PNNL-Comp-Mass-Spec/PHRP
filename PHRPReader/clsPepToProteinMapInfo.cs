using System;
using System.Collections.Generic;

namespace PHRPReader
{
    /// <summary>
    /// Track location of a given peptide in one or more proteins
    /// </summary>
    public class clsPepToProteinMapInfo
    {
        /// <summary>
        /// Start and end residue locations in a protein
        /// </summary>
        public class udtProteinLocationInfo
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
        private readonly Dictionary<string, List<udtProteinLocationInfo>> mProteinMapInfo;

        /// <summary>
        /// Number of proteins that contain a given peptide
        /// </summary>
        public int ProteinCount => ProteinMapInfo.Count;

        /// <summary>
        /// Dictionary of protein names and residue start/end positions for a given peptide
        /// </summary>
        public IReadOnlyDictionary<string, List<udtProteinLocationInfo>> ProteinMapInfo => mProteinMapInfo;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="proteinName">Protein name</param>
        /// <param name="residueStart">Location that a peptide starts in the protein</param>
        /// <param name="residueEnd">Location that a peptide ends in the protein</param>
        public clsPepToProteinMapInfo(string proteinName, int residueStart, int residueEnd)
        {
            mProteinMapInfo = new Dictionary<string, List<udtProteinLocationInfo>>(StringComparer.CurrentCultureIgnoreCase);

            AddProtein(proteinName, residueStart, residueEnd);
        }

        public void AddProtein(string proteinName, int residueStart, int residueEnd)
        {
            List<udtProteinLocationInfo> lstLocations = null;

            if (mProteinMapInfo.TryGetValue(proteinName, out lstLocations))
            {
                // Protein mapping already exists; check residueStart
                foreach (var udtLoc in lstLocations)
                {
                    if (udtLoc.ResidueStart == residueStart)
                    {
                        // Update this entry
                        if (udtLoc.ResidueEnd != residueEnd)
                        {
                            udtLoc.ResidueEnd = residueEnd;
                        }
                        return;
                    }
                }

                var udtLocInfoAddnl = new udtProteinLocationInfo();
                udtLocInfoAddnl.ResidueStart = residueStart;
                udtLocInfoAddnl.ResidueEnd = residueEnd;

                lstLocations.Add(udtLocInfoAddnl);
                return;
            }

            var udtLocInfo = new udtProteinLocationInfo();
            udtLocInfo.ResidueStart = residueStart;
            udtLocInfo.ResidueEnd = residueEnd;

            mProteinMapInfo.Add(proteinName, new List<udtProteinLocationInfo> { udtLocInfo });
        }
    }
}
