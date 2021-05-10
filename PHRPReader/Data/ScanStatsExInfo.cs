namespace PHRPReader.Data
{
    /// <summary>
    /// Data loaded from a ScanStatsEx file
    /// </summary>
    public class ScanStatsExInfo
    {
        /// <summary>
        /// Scan number
        /// </summary>
        public int ScanNumber { get; }

        /// <summary>
        /// Ion injection time
        /// </summary>
        public double IonInjectionTime { get; set; }

        /// <summary>
        /// scan event
        /// </summary>
        public int ScanEvent { get; set; }

        /// <summary>
        /// Master index
        /// </summary>
        public int MasterIndex { get; set; }

        /// <summary>
        /// Elapsed scan time (in seconds)
        /// </summary>
        /// <remarks>This is the time required to acquire the mass spectrum</remarks>
        public double ElapsedScanTime { get; set; }

        /// <summary>
        /// charge state
        /// </summary>
        public int ChargeState { get; set; }

        /// <summary>
        /// Monoisotopic m/z of the parent ion
        /// </summary>
        public double MonoisotopicMZ { get; set; }

        /// <summary>
        /// MS2 isolation width
        /// </summary>
        public double MS2IsolationWidth { get; set; }

        /// <summary>
        /// FT analyzer settings
        /// </summary>
        public string FTAnalyzerSettings { get; set; }

        /// <summary>
        /// FT analyzer message
        /// </summary>
        public string FTAnalyzerMessage { get; set; }

        /// <summary>
        /// FT resolution
        /// </summary>
        public double FTResolution { get; set; }

        /// <summary>
        /// Conversion parameter B
        /// </summary>
        public double ConversionParameterB { get; set; }

        /// <summary>
        /// Conversion parameter C
        /// </summary>
        public double ConversionParameterC { get; set; }

        /// <summary>
        /// Conversion parameter D
        /// </summary>
        public double ConversionParameterD { get; set; }

        /// <summary>
        /// Conversion parameter E
        /// </summary>
        public double ConversionParameterE { get; set; }

        /// <summary>
        /// Collision mode
        /// </summary>
        public string CollisionMode { get; set; }

        /// <summary>
        /// Scan filter
        /// </summary>
        public string ScanFilterText { get; set; }

        /// <summary>
        /// Source voltage
        /// </summary>
        public double SourceVoltage { get; set; }

        /// <summary>
        /// Source current
        /// </summary>
        public double Source_Current { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="scanNumber"></param>
        public ScanStatsExInfo(int scanNumber)
        {
            ScanNumber = scanNumber;
        }

        /// <summary>
        /// Show the scan number and scan filter
        /// </summary>
        public override string ToString()
        {
            return string.Format("Scan {0}, {1}", ScanNumber, ScanFilterText);
        }
    }
}
