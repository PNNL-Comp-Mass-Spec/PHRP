namespace PHRPReader.Data
{
    /// <summary>
    /// Data loaded from a ScanStats file
    /// </summary>
    public class ScanStatsInfo
    {
        private string mScanTypeName;

        /// <summary>
        /// Scan number
        /// </summary>
        public int ScanNumber { get; }

        /// <summary>
        /// Scan time, in minutes
        /// </summary>
        public float ScanTimeMinutes { get; set; }

        /// <summary>
        /// Scan time, in minutes
        /// </summary>
        public string ScanTimeMinutesText { get; set; }

        /// <summary>
        /// Scan type (aka scan level)
        /// </summary>
        /// <remarks>1 for MS1, 2 for MS2</remarks>
        public int ScanType { get; set; }

        /// <summary>
        /// Total ion intensity (TIC)
        /// </summary>
        public double TotalIonIntensity { get; set; }

        /// <summary>
        /// Total ion intensity (TIC)
        /// </summary>
        public string TotalIonIntensityText { get; set; }

        /// <summary>
        /// Base peak intensity (BPI)
        /// </summary>
        public double BasePeakIntensity { get; set; }

        /// <summary>
        /// Base peak intensity (BPI)
        /// </summary>
        public string BasePeakIntensityText { get; set; }

        /// <summary>
        /// Base peak m/z
        /// </summary>
        public double BasePeakMZ { get; set; }

        /// <summary>
        /// Base peak m/z
        /// </summary>
        public string BasePeakMzText { get; set; }

        /// <summary>
        /// Base peak signal to noise ratio (S/N)
        /// </summary>
        public double BasePeakSignalToNoiseRatio { get; set; }

        /// <summary>
        /// Base peak signal to noise ratio (S/N)
        /// </summary>
        public string BasePeakSignalToNoiseRatioText { get; set; }

        /// <summary>
        /// Ion count (after filters)
        /// </summary>
        public int IonCount { get; set; }

        /// <summary>
        /// Ion count (before filtering)
        /// </summary>
        public int IonCountRaw { get; set; }

        /// <summary>
        /// Scan type name
        /// </summary>
        public string ScanTypeName
        {
            get => mScanTypeName;
            set => mScanTypeName = string.IsNullOrEmpty(value) ? string.Empty : value;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="scanNumber"></param>
        /// <param name="scanTimeMinutes"></param>
        /// <param name="scanType"></param>
        public ScanStatsInfo(int scanNumber, float scanTimeMinutes, int scanType)
        {
            ScanNumber = scanNumber;
            ScanTimeMinutes = scanTimeMinutes;
            ScanType = scanType;
        }
    }
}
