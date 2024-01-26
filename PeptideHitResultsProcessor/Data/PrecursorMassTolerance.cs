using System;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// Tracks precursor mass tolerance (aka parent mass tolerance),
    /// used by the search tool when finding candidate peptides for a given MS/MS spectrum
    /// </summary>
    /// <remarks>
    /// <para>
    /// If a parameter file specifies a tolerance of 20ppm, this class will have
    /// ToleranceLeft=20, ToleranceRight=20, and ToleranceIsPPM=True
    /// </para>
    /// <para>
    /// If a MS-GF+ parameter file specifies a tolerance of 0.5Da,2.5Da, this class will have
    /// ToleranceLeft=0.5, ToleranceRight=2.5, and ToleranceIsPPM=False
    /// </para>
    /// </remarks>
    internal class PrecursorMassTolerance
    {
        // Ignore Spelling: Da

        private double mToleranceLeft;

        private double mToleranceRight;

        /// <summary>
        /// Mass tolerance to the left of a m/z value (i.e., more negative)
        /// </summary>
        /// <remarks>
        /// MSFragger tracks this as a negative value, but this class stores it as a positive number
        /// </remarks>
        public double ToleranceLeft
        {
            get => mToleranceLeft;
            set => mToleranceLeft = Math.Abs(value);
        }

        /// <summary>
        /// Mass tolerance to the right of a m/z value (i.e., more positive)
        /// </summary>
        /// <remarks>
        /// Using a backing variable here for symmetry with property ToleranceLeft
        /// </remarks>
        public double ToleranceRight
        {
            get => mToleranceRight;
            set => mToleranceRight = Math.Abs(value);
        }

        /// <summary>
        /// True if the tolerance is specified in parts per million
        /// </summary>
        public bool IsPPM { get; set; }

        /// <summary>
        /// Constructor for a symmetric tolerance
        /// </summary>
        /// <param name="tolerance">Tolerance</param>
        /// <param name="isPPM">True if the tolerance is in PPM</param>
        public PrecursorMassTolerance(double tolerance, bool isPPM) : this(tolerance, tolerance, isPPM)
        {
        }

        /// <summary>
        /// Constructor for an asymmetric tolerance
        /// </summary>
        /// <param name="toleranceLeft">Tolerance to the left</param>
        /// <param name="toleranceRight">Tolerance to the right</param>
        /// <param name="isPPM">True if the tolerance is in PPM</param>
        public PrecursorMassTolerance(double toleranceLeft, double toleranceRight, bool isPPM)
        {
            ToleranceLeft = toleranceLeft;
            ToleranceRight = toleranceRight;
            IsPPM = isPPM;
        }

        /// <summary>
        /// Clear tolerance values
        /// </summary>
        public void Clear()
        {
            ToleranceLeft = 0;
            ToleranceRight = 0;
            IsPPM = false;
        }

        /// <summary>
        /// Show the tolerance and units
        /// </summary>
        public override string ToString()
        {
            string units;
            double equivalenceThreshold;

            if (IsPPM)
            {
                units = "ppm";
                equivalenceThreshold = 0.01;
            }
            else
            {
                units = "Da";
                equivalenceThreshold = 0.0001;
            }

            if (Math.Abs(ToleranceLeft - ToleranceRight) < equivalenceThreshold)
            {
                return "+/-" + ToleranceLeft + " " + units;
            }

            return "-" + ToleranceRight + ", +" + ToleranceLeft + " " + units;
        }
    }
}
