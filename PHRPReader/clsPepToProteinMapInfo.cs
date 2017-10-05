using System;
using System.Collections.Generic;

namespace PHRPReader
{
    public class clsPepToProteinMapInfo
    {
        public class udtProteinLocationInfo
        {
            public int ResidueStart;
            public int ResidueEnd;
        }

        /// <summary>
        /// Dictionary of protein names and residue start/end positions
        /// </summary>
        /// <remarks></remarks>
        private readonly Dictionary<string, List<udtProteinLocationInfo>> mProteinMapInfo;

        public int ProteinCount
        {
            get { return mProteinMapInfo.Count; }
        }

        public Dictionary<string, List<udtProteinLocationInfo>> ProteinMapInfo
        {
            get { return mProteinMapInfo; }
        }

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
