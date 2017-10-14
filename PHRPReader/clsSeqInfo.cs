//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/20/2012
//
// This class tracks the sequence information determined by PHRP and stored in a _SeqInfo.txt file
//
//*********************************************************************************************************

namespace PHRPReader
{
    public class clsSeqInfo
    {
        /// <summary>
        /// Sequence ID
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int SeqID { get; }

        /// <summary>
        /// Number of modifications
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int ModCount { get; }

        /// <summary>
        /// Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ModDescription { get; }

        /// <summary>
        /// Theoretical, monoisotopic mass (including the modified residues)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Constructor using Sequence ID and mass
        /// </summary>
        /// <remarks></remarks>
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
        /// <remarks></remarks>
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
        /// <remarks></remarks>
        public void UpdateMonoisotopicMass(double monoMass)
        {
            MonoisotopicMass = monoMass;
        }
    }
}
