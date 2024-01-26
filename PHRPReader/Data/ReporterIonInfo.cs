using System.Collections.Generic;

namespace PHRPReader.Data
{
    /// <summary>
    /// This class tracks data read from a _ReporterIons.txt file
    /// </summary>
    public class ReporterIonInfo
    {
        /// <summary>
        /// Scan number
        /// </summary>
        public int ScanNumber { get; }

        /// <summary>
        /// Collision mode
        /// </summary>
        public string CollisionMode { get; set; }

        /// <summary>
        /// Parent ion m/z
        /// </summary>
        public double ParentIonMZ { get; set; }

        /// <summary>
        /// Parent ion m/z
        /// </summary>
        public string ParentIonMzText { get; set; }

        /// <summary>
        /// Base peak intensity
        /// </summary>
        public double BasePeakIntensity { get; set; }

        /// <summary>
        /// Base peak m/z
        /// </summary>
        public double BasePeakMZ { get; set; }

        /// <summary>
        /// Parent scan number
        /// </summary>
        public int ParentScan { get; set; }

        /// <summary>
        /// Maximum reporter ion intensity
        /// </summary>
        public double ReporterIonIntensityMax { get; set; }

        /// <summary>
        /// Weighted average intensity correction percentage
        /// </summary>
        public double WeightedAvgPctIntensityCorrection { get; set; }

        /// <summary>
        /// Reporter ion data for this scan
        /// </summary>
        public List<ReporterIonStats> ReporterIons { get; }

        /// <summary>
        /// Tab-delimited list of data in the reporter-ion related columns, starting with the ReporterIonIntensityMax column
        /// </summary>
        /// <remarks>
        /// The MASICResultsMerger program uses this property
        /// </remarks>
        public string ReporterIonColumnData { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="scanNumber">Scan number</param>
        public ReporterIonInfo(int scanNumber)
        {
            ScanNumber = scanNumber;

            ReporterIons = new List<ReporterIonStats>();
        }

        /// <summary>
        /// Show the scan number
        /// </summary>
        public override string ToString()
        {
            return string.Format("Scan {0}", ScanNumber);
        }
    }
}
