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
        private readonly int mSeqID;
        private readonly int mModCount;
        private readonly string mModDescription;
        private double mMonoisotopicMass;

        /// <summary>
        /// Sequence ID
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int SeqID
        {
            get { return mSeqID; }
        }

        /// <summary>
        /// Number of modifications
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public int ModCount
        {
            get { return mModCount; }
        }

        /// <summary>
        /// Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public string ModDescription
        {
            get { return mModDescription; }
        }

        /// <summary>
        /// Theoretical, monoisotopic mass (including the modified residues)
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public double MonoisotopicMass
        {
            get { return mMonoisotopicMass; }
        }

        /// <summary>
        /// Constructor using Sequence ID and mass
        /// </summary>
        /// <remarks></remarks>
        public clsSeqInfo(int SeqID, double MonoisotopicMass) : this(SeqID, MonoisotopicMass, 0, string.Empty)
        {
        }

        /// <summary>
        /// Constructor using Sequence ID, mass, mod count, and list of modifications
        /// </summary>
        /// <param name="SeqID">Sequence ID</param>
        /// <param name="MonoisotopicMass">Theoretical, monoisotopic mass (including the modified residues)</param>
        /// <param name="ModCount">Number of modifications</param>
        /// <param name="ModDescription">Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"</param>
        /// <remarks></remarks>
        public clsSeqInfo(int SeqID, double MonoisotopicMass, int ModCount, string ModDescription)
        {
            mSeqID = SeqID;
            mMonoisotopicMass = MonoisotopicMass;
            mModCount = ModCount;

            if (string.IsNullOrEmpty(ModDescription))
            {
                mModDescription = string.Empty;
            }
            else
            {
                mModDescription = ModDescription;
            }
        }

        /// <summary>
        /// Update the monoisotopic mass for this sequence
        /// </summary>
        /// <param name="monoMass"></param>
        /// <remarks></remarks>
        public void UpdateMonoisotopicMass(double monoMass)
        {
            mMonoisotopicMass = monoMass;
        }
    }
}
