namespace PHRPReader.Data
{
    public class PrecursorInfo
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
        /// Scan type (aka scan level)
        /// </summary>
        /// <remarks>1 for MS1, 2 for MS2</remarks>
        public int ScanType { get; set; }

        /// <summary>
        /// Scan type name
        /// </summary>
        public string ScanTypeName
        {
            get => mScanTypeName;
            set => mScanTypeName = string.IsNullOrEmpty(value) ? string.Empty : value;
        }

        /// <summary>
        /// Precursor m/z
        /// </summary>
        /// <remarks>
        /// For MS3 spectra, this is the m/z of the preceding MS2 scan
        /// </remarks>
        public double PrecursorMz { get; set; }

        /// <summary>
        /// Scan filter
        /// </summary>
        public string ScanFilterText { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="scanNumber"></param>
        /// <param name="scanTimeMinutes"></param>
        /// <param name="scanType"></param>
        public PrecursorInfo(int scanNumber, float scanTimeMinutes, int scanType)
        {
            ScanNumber = scanNumber;
            ScanTimeMinutes = scanTimeMinutes;
            ScanType = scanType;
        }
    }
}
