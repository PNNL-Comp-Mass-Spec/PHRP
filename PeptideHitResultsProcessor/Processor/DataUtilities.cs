using System;
using System.Collections.Generic;
using System.Globalization;

namespace PeptideHitResultsProcessor.Processor
{
    internal static class DataUtilities
    {
        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to string.Empty
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        public static bool GetColumnValue(string[] splitLine, int columnIndex, out string value)
        {
            return GetColumnValue(splitLine, columnIndex, out value, string.Empty);
        }
        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to 0
        /// </summary>
        /// <returns>True if columnIndex >= 0 and an integer value is present</returns>
        public static bool GetColumnValue(string[] splitLine, int columnIndex, out int value)
        {
            return GetColumnValue(splitLine, columnIndex, out value, 0);
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        public static bool GetColumnValue(string[] splitLine, int columnIndex, out string value, string valueIfMissing)
        {
            if (columnIndex >= 0 && columnIndex < splitLine.Length)
            {
                value = splitLine[columnIndex];
                return true;
            }

            value = valueIfMissing;
            return false;
        }
        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to valueIfMissing
        /// </summary>
        /// <returns>True if columnIndex >= 0 and an integer value is present</returns>
        public static bool GetColumnValue(string[] splitLine, int columnIndex, out int value, int valueIfMissing)
        {
            if (GetColumnValue(splitLine, columnIndex, out var valueText, valueIfMissing.ToString()))
            {
                if (int.TryParse(valueText, out value))
                {
                    return true;
                }

                value = valueIfMissing;
                return false;
            }

            value = valueIfMissing;
            return false;
        }

        /// <summary>
        /// If columnIndex is >= 0, updates value with the value at splitLine[columnIndex]
        /// Otherwise, updates value to string.Empty
        /// </summary>
        /// <returns>True if columnIndex >= 0</returns>
        public static bool GetColumnValueCheckNaN(string[] splitLine, int columnIndex, out string value)
        {
            var valueDefined = GetColumnValue(splitLine, columnIndex, out value, string.Empty);

            if (valueDefined && value.Equals("NaN", StringComparison.OrdinalIgnoreCase))
                value = string.Empty;

            return valueDefined;
        }

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
