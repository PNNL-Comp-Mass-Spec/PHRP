//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/20/2012
//
//*********************************************************************************************************

namespace PHRPReader
{
    /// <summary>
    /// This class tracks the sequence information determined by PHRP and stored in a _SeqInfo.txt file
    /// </summary>
    public class clsSeqInfo
    {
        /// <summary>
        /// Sequence ID
        /// </summary>
        public int SeqID { get; }

        /// <summary>
        /// Number of modifications
        /// </summary>
        public int ModCount { get; }

        /// <summary>
        /// Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"
        /// </summary>
        public string ModDescription { get; }

        /// <summary>
        /// Theoretical, monoisotopic mass (including the modified residues)
        /// </summary>
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Constructor using Sequence ID and mass
        /// </summary>
        public clsSeqInfo(int seqID, double monoisotopicMass) : this(seqID, monoisotopicMass, 0, string.Empty)
        {
        }

        /// <summary>
        /// Constructor using Sequence ID, mass, mod count, and list of modifications
        /// </summary>
        /// <param name="seqID">Sequence ID</param>
        /// <param name="monoisotopicMass">Theoretical, monoisotopic mass (including the modified residues)</param>
        /// <param name="modCount">Number of modifications</param>
        /// <param name="modDescription">Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"</param>
        public clsSeqInfo(int seqID, double monoisotopicMass, int modCount, string modDescription)
        {
            SeqID = seqID;
            MonoisotopicMass = monoisotopicMass;
            ModCount = modCount;

            ModDescription = string.IsNullOrEmpty(modDescription) ? string.Empty : modDescription;
        }

        /// <summary>
        /// Update the monoisotopic mass for this sequence
        /// </summary>
        /// <param name="monoMass"></param>
        public void UpdateMonoisotopicMass(double monoMass)
        {
            MonoisotopicMass = monoMass;
        }

        /// <summary>
        /// Show the sequence ID and monoisotopic mass
        /// </summary>
        public override string ToString()
        {
            return string.Format("ID {0}: {1:F3} Da", SeqID, MonoisotopicMass);
        }
    }
}
