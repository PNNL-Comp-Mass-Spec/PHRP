namespace PHRPReader.Data
{
    /// <summary>
    /// Details of an observed reporter ion
    /// </summary>
    public class ReporterIonStats
    {
        /// <summary>
        /// Reporter ion m/z
        /// </summary>
        public double MZ { get; }

        /// <summary>
        /// Reporter ion intensity
        /// </summary>
        public double Intensity { get; set; }

        /// <summary>
        /// Reporter ion intensity, before interference correction
        /// </summary>
        public double? OriginalIntensity { get; set; }

        /// <summary>
        /// Signal to noise ratio
        /// </summary>
        public double? SignalToNoise { get; set; }

        /// <summary>
        /// Resolution
        /// </summary>
        public double? Resolution { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        public ReporterIonStats(double reporterIonMz)
        {
            MZ = reporterIonMz;
        }
    }
}
