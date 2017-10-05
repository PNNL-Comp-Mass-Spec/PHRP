namespace PHRPReader
{
    public class clsScanStatsInfo
    {
        private readonly int mScanNumber;

        private string mScanTypeName;
        // ReSharper disable once ConvertToAutoProperty
        public int ScanNumber
        {
            get { return mScanNumber; }
        }

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
            get { return mScanTypeName; }
            set
            {
                if (string.IsNullOrEmpty(value))
                {
                    mScanTypeName = string.Empty;
                }
                else
                {
                    mScanTypeName = value;
                }
            }
        }

        public clsScanStatsInfo(int intScanNumber, float sngScanTimeMinutes, int intScanType)
        {
            mScanNumber = intScanNumber;
            ScanTimeMinutes = sngScanTimeMinutes;
            ScanType = intScanType;
        }
    }
}
