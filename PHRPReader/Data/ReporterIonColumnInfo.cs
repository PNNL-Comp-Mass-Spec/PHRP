namespace PHRPReader.Data
{
    /// <summary>
    /// Information about reporter ion related columns in a _ReporterIons.txt file
    /// </summary>
    public class ReporterIonColumnInfo
    {
        /// <summary>
        /// Index of the column with reporter ion intensity values
        /// </summary>
        public int IntensityColumnIndex { get; }

        /// <summary>
        /// Name of the column with reporter ion intensity values
        /// </summary>
        public string IntensityColumnName { get; }

        /// <summary>
        /// Reporter ion m/z
        /// </summary>
        public double MZ { get; }

        /// <summary>
        /// Index of the column with reporter ion original intensity values
        /// </summary>
        /// <remarks>
        /// These columns will only be present if the MASIC parameter file had ReporterIonSaveUncorrectedIntensities=True
        /// </remarks>
        public int OriginalIntensityColumnIndex { get; set; }

        /// <summary>
        /// Index of the column with reporter ion signal/noise values
        /// </summary>
        public int SignalToNoiseColumnIndex { get; set; }

        /// <summary>
        /// Index of the column with reporter ion resolution values
        /// </summary>
        public int ResolutionColumnIndex { get; set; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="intensityColumnIndex">Intensity column index</param>
        /// <param name="intensityColumnName">Intensity column name</param>
        /// <param name="reporterIonMz">Reporter ion m/z</param>
        public ReporterIonColumnInfo(int intensityColumnIndex, string intensityColumnName, double reporterIonMz)
        {
            IntensityColumnIndex = intensityColumnIndex;
            IntensityColumnName = intensityColumnName;
            MZ = reporterIonMz;
        }

        /// <summary>
        /// Show the intensity column name
        /// </summary>
        public override string ToString()
        {
            return IntensityColumnName;
        }
    }
}
