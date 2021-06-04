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

        /// <summary>
        /// Mass tolerance to the left of a m/z value (i.e., more negative)
        /// </summary>
        public double ToleranceLeft { get; set; }

        /// <summary>
        /// Mass tolerance to the right of a m/z value (i.e., more positive)
        /// </summary>
        public double ToleranceRight { get; set; }

        /// <summary>
        /// True if the tolerance is specified in parts per million
        /// </summary>
        public bool IsPPM { get; set; }

        /// <summary>
        /// Parameter-less constructor
        /// </summary>
        public PrecursorMassTolerance() : this(0, 0, false)
        {
        }

        /// <summary>
        /// Constructor for a symmetric tolerance
        /// </summary>
        /// <param name="tolerance"></param>
        /// <param name="isPPM"></param>
        public PrecursorMassTolerance(double tolerance, bool isPPM) : this(tolerance, tolerance, isPPM)
        {
        }

        /// <summary>
        /// Constructor for an asymmetric tolerance
        /// </summary>
        /// <param name="toleranceLeft"></param>
        /// <param name="toleranceRight"></param>
        /// <param name="isPPM"></param>
        public PrecursorMassTolerance(double toleranceLeft, double toleranceRight, bool isPPM)
        {
            ToleranceLeft = toleranceLeft;
            ToleranceLeft = toleranceRight;
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
        /// Return the tolerance, as a string
        /// </summary>
        /// <returns></returns>
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
