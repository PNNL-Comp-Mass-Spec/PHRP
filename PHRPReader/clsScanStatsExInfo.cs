namespace PHRPReader
{
    public class clsScanStatsExInfo
    {
        private readonly int mScanNumber;

        // ReSharper disable once ConvertToAutoProperty
        public int ScanNumber
        {
            get { return mScanNumber; }
        }

        public double IonInjectionTime { get; set; }

        public int ScanEvent { get; set; }

        public int MasterIndex { get; set; }

        public double ElapsedScanTime { get; set; }

        public int ChargeState { get; set; }

        public double MonoisotopicMZ { get; set; }

        public double MS2IsolationWidth { get; set; }

        public string FTAnalyzerSettings { get; set; }

        public string FTAnalyzerMessage { get; set; }

        public double FTResolution { get; set; }

        public double ConversionParameterB { get; set; }

        public double ConversionParameterC { get; set; }

        public double ConversionParameterD { get; set; }

        public double ConversionParameterE { get; set; }

        public string CollisionMode { get; set; }

        public string ScanFilterText { get; set; }

        public double SourceVoltage { get; set; }

        public double Source_Current { get; set; }

        public clsScanStatsExInfo(int intScanNumber)
        {
            mScanNumber = intScanNumber;
        }
    }
}
