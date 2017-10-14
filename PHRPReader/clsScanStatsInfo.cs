namespace PHRPReader
{
    public class clsScanStatsInfo
    {
        private readonly int mScanNumber;

        private string mScanTypeName;
        // ReSharper disable once ConvertToAutoProperty
        public int ScanNumber => mScanNumber;

        public float ScanTimeMinutes { get; set; }

        public int ScanType { get; set; }

        public double TotalIonIntensity { get; set; }

        public double BasePeakIntensity { get; set; }

        public double BasePeakMZ { get; set; }

        public double BasePeakSignalToNoiseRatio { get; set; }

        public int IonCount { get; set; }

        public int IonCountRaw { get; set; }

        public string ScanTypeName
        {
            get => mScanTypeName;
            set => mScanTypeName = string.IsNullOrEmpty(value) ? string.Empty : value;
        }

        public clsScanStatsInfo(int intScanNumber, float sngScanTimeMinutes, int intScanType)
        {
            mScanNumber = intScanNumber;
            ScanTimeMinutes = sngScanTimeMinutes;
            ScanType = intScanType;
        }
    }
}
