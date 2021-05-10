namespace PHRPReader.Data
{
    /// <summary>
    /// Data loaded from a SICStats file
    /// </summary>
    public class SICStatsInfo
    {
        /// <summary>
        /// Dataset ID
        /// </summary>
        public int Dataset { get; set; }

        /// <summary>
        /// Parent Ion Index
        /// </summary>
        public int ParentIonIndex { get; set; }

        /// <summary>
        /// Precursor m/z
        /// </summary>
        public double MZ { get; set; }

        /// <summary>
        /// Precursor m/z
        /// </summary>
        public string MzText { get; set; }

        /// <summary>
        /// Survey Scan Number
        /// </summary>
        public int SurveyScanNumber { get; set; }

        /// <summary>
        /// Fragmentation Scan Number
        /// </summary>
        public int FragScanNumber { get; set; }

        /// <summary>
        /// Optimal Peak Apex Scan Number
        /// </summary>
        public int OptimalPeakApexScanNumber { get; set; }

        /// <summary>
        /// Peak Apex Override Parent Ion Index
        /// </summary>
        public int PeakApexOverrideParentIonIndex { get; set; }

        /// <summary>
        /// Custom SIC Peak
        /// </summary>
        public int CustomSICPeak { get; set; }

        /// <summary>
        /// Peak Scan Start
        /// </summary>
        public int PeakScanStart { get; set; }

        /// <summary>
        /// Peak Scan End
        /// </summary>
        public int PeakScanEnd { get; set; }

        /// <summary>
        /// Peak Scan Max Intensity
        /// </summary>
        public int PeakScanMaxIntensity { get; set; }

        /// <summary>
        /// Peak Max Intensity
        /// </summary>
        public double PeakMaxIntensity { get; set; }

        /// <summary>
        /// Peak Max Intensity
        /// </summary>
        public string PeakMaxIntensityText { get; set; }

        /// <summary>
        /// Peak Signal To Noise Ratio
        /// </summary>
        public double PeakSignalToNoiseRatio { get; set; }

        /// <summary>
        /// Peak Signal To Noise Ratio
        /// </summary>
        public string PeakSignalToNoiseRatioText { get; set; }

        /// <summary>
        /// FWHM In Scans
        /// </summary>
        public int FWHMInScans { get; set; }

        /// <summary>
        /// Peak Area
        /// </summary>
        public double PeakArea { get; set; }

        /// <summary>
        /// Peak Area
        /// </summary>
        public string PeakAreaText { get; set; }

        /// <summary>
        /// Parent Ion Intensity
        /// </summary>
        public double ParentIonIntensity { get; set; }

        /// <summary>
        /// Parent Ion Intensity
        /// </summary>
        public string ParentIonIntensityText { get; set; }

        /// <summary>
        /// Peak Baseline Noise Level
        /// </summary>
        public double PeakBaselineNoiseLevel { get; set; }

        /// <summary>
        /// Peak Baseline Noise Level
        /// </summary>
        public string PeakBaselineNoiseLevelText { get; set; }

        /// <summary>
        /// Peak Baseline Noise StDev
        /// </summary>
        public double PeakBaselineNoiseStDev { get; set; }

        /// <summary>
        /// Peak Baseline Noise StDev
        /// </summary>
        public string PeakBaselineNoiseStDevText { get; set; }

        /// <summary>
        /// Peak Baseline Points Used
        /// </summary>
        public int PeakBaselinePointsUsed { get; set; }

        /// <summary>
        /// Stat Moments Area
        /// </summary>
        public double StatMomentsArea { get; set; }

        /// <summary>
        /// Stat Moments Area
        /// </summary>
        public string StatMomentsAreaText { get; set; }

        /// <summary>
        /// Center Of Mass Scan
        /// </summary>
        public int CenterOfMassScan { get; set; }

        /// <summary>
        /// Peak StDev
        /// </summary>
        public double PeakStDev { get; set; }

        /// <summary>
        /// Peak StDev
        /// </summary>
        public string PeakStDevText { get; set; }

        /// <summary>
        /// Peak Skew
        /// </summary>
        public double PeakSkew { get; set; }

        /// <summary>
        /// Peak Skew
        /// </summary>
        public string PeakSkewText { get; set; }

        /// <summary>
        /// Peak KSStat
        /// </summary>
        public double PeakKSStat { get; set; }

        /// <summary>
        /// Peak KSStat
        /// </summary>
        public string PeakKSStatText { get; set; }

        /// <summary>
        /// Stat Moments Data Count Used
        /// </summary>
        public int StatMomentsDataCountUsed { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="parentIonIndex"></param>
        public SICStatsInfo(int parentIonIndex)
        {
            ParentIonIndex = parentIonIndex;
        }

        /// <summary>
        /// Show the scan number and scan filter
        /// </summary>
        public override string ToString()
        {
            return string.Format("Parent Ion Index {0}, {1} m/z", ParentIonIndex, MZ);
        }
    }
}
