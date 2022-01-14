using System.Collections.Generic;

namespace PeptideHitResultsProcessor.Processor
{
    internal static class DataUtilities
    {
        /// <summary>
        /// Smooth data using a moving average window
        /// </summary>
        /// <remarks>
        /// <para>
        /// Based on code by J.L. Haynes, https://stackoverflow.com/a/44318312/1179467
        /// </para>
        /// <para>
        /// Math.DotNet alternatives
        /// https://numerics.mathdotnet.com/api/MathNet.Numerics.Statistics/MovingStatistics.htm
        /// https://numerics.mathdotnet.com/api/MathNet.Numerics.Statistics/RunningStatistics.htm
        /// </para>
        /// </remarks>
        /// <param name="values"></param>
        /// <param name="windowSize"></param>
        /// <returns>List of smoothed values</returns>
        public static List<double> MovingAverageSmooth(IEnumerable<double> values, int windowSize = 10)
        {
            var currentValues = new Queue<double>();
            double runningSum = 0;

            var smoothedData = new List<double>();

            foreach (var value in values)
            {
                runningSum += value;
                currentValues.Enqueue(value);

                if (currentValues.Count > windowSize)
                {
                    runningSum -= currentValues.Dequeue();
                }

                smoothedData.Add(runningSum / currentValues.Count);
            }

            return smoothedData;
        }
    }
}
