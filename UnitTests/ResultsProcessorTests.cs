using System;
using System.Collections.Generic;
using NUnit.Framework;
using PeptideHitResultsProcessor.Processor;

namespace PHRP_UnitTests
{
    /// <summary>
    /// Unit tests for MultiDatasetResultsProcessor.GetDatasetNameMap and various StringUtilities methods
    /// </summary>
    public class ResultsProcessorTests
    {
        // ReSharper disable once CommentTypo

        // Ignore Spelling: Acetamiprid_rep2_1000V, Da, Mam

        // ReSharper disable StringLiteralTypo

        [Test]
        [TestCase("", "")]
        [TestCase(" ", " ")]
        [TestCase("QC_Shew_16_01", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14")]
        [TestCase("QC_Shew_16_01-15f", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14", "QC_Shew_16_01-15f_09_15Nov16_Tiger_16-02-15")]
        [TestCase("QC_Shew_16_01-15f", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14", "QC_Shew_16_01-15f_10_15Nov16_Tiger_16-02-14")]
        [TestCase("QC_Shew_16_01-15f", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14", "QC_Shew_16_01-15f_09_15Nov16_Tiger_16-02-15", "QC_Shew_16_01-15f_10_15Nov16_Tiger_16-02-14")]
        [TestCase("QC", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14", "QC_PP_MCF-7_21_01_d_22Apr21_Rage_Rep-21-02-06")]
        [TestCase("QC", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14", "QC_Shew_16_01-15f_09_15Nov16_Tiger_16-02-15", "QC_PP_MCF-7_21_01_d_22Apr21_Rage_Rep-21-02-06")]
        [TestCase("QC", "QC_Shew_16_01-15f_09_15Nov16_Tiger_16-02-15", "QC_Mam_19_01_a_22Apr21_Rage_Rep-21-02-06")]
        [TestCase("QC", "QC_PP_MCF-7_21_01_d_22Apr21_Rage_Rep-21-02-06", "QC_Mam_19_01_a_22Apr21_Rage_Rep-21-02-06")]
        [TestCase("", "QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14", "AALEDTLAETEAR_MSonly")]
        [TestCase("AALEDTLAETEAR_MS", "AALEDTLAETEAR_MSMSonly", "AALEDTLAETEAR_MSonly")]
        [TestCase("Acetamiprid_rep1_1000V", "Acetamiprid_rep1_1000V", "Acetamiprid_rep1_1000V_neg")]
        [TestCase("Acetamiprid_rep", "Acetamiprid_rep1_1000V", "Acetamiprid_rep1_1000V_neg", "Acetamiprid_rep2_1000V", "Acetamiprid_rep2_1000V_neg")]
        [TestCase("Acetamiprid_rep2_1", "Acetamiprid_rep2_1000V", "Acetamiprid_rep2_1000V_neg", "Acetamiprid_rep2_1100V", "Acetamiprid_rep2_1100V_neg")]
        // ReSharper restore StringLiteralTypo
        public void TestGetDatasetNameMap(string expectedLongestCommonString, params string[] datasets)
        {
            var datasetNames = new SortedSet<string>();

            foreach (var item in datasets)
            {
                if (item.Length > 0)
                    datasetNames.Add(item);
            }

            var baseDatasetNames = MultiDatasetResultsProcessor.GetDatasetNameMap(datasetNames, out var longestCommonString, out var warnings);

            foreach (var warning in warnings)
            {
                Console.WriteLine(warning);
                Console.WriteLine();
            }

            foreach (var item in baseDatasetNames)
            {
                Console.WriteLine("{0} =>", item.Key);
                Console.WriteLine(item.Value);
                Console.WriteLine();
            }

            Console.WriteLine("Longest common string:");
            Console.WriteLine(longestCommonString);

            Assert.AreEqual(expectedLongestCommonString, longestCommonString);

            Assert.AreEqual(0, warnings.Count, "Warnings were reported");
        }

        [Test]
        [TestCase("0", 0)]
        [TestCase("5", 5)]
        [TestCase("5.678", 6)]
        [TestCase("5.6789012345", 6)]
        [TestCase("5.67890123456789012345", 6)]
        [TestCase("5.67x", 0)]
        [TestCase("5.67 Da", 0)]
        [TestCase("-0", 0)]
        [TestCase("-5", -5)]
        [TestCase("-5.678", -6)]
        [TestCase("-5.67x", 0)]
        [TestCase("-5.67 Da", 0)]
        [TestCase("Text 5.67", 0)]
        public void TestCIntSafe(string value, int expectedResult)
        {
            var result = StringUtilities.CIntSafe(value, 0);

            Console.WriteLine("{0,-7} =>  {1}", value, result);

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase("0", 0)]
        [TestCase("5", 5)]
        [TestCase("5.678", 5.678)]
        [TestCase("5.6789012345", 5.6789012345)]
        [TestCase("5.67890123456789012345", 5.6789012345678902)]
        [TestCase("5.67x", 0)]
        [TestCase("5.67 Da", 0)]
        [TestCase("-0", 0)]
        [TestCase("-5", -5)]
        [TestCase("-5.678", -5.678)]
        [TestCase("-5.67x", 0)]
        [TestCase("-5.67 Da", 0)]
        [TestCase("X 5.67", 0)]
        public void TestCDblSafe(string value, double expectedResult)
        {
            var result = StringUtilities.CDblSafe(value, 0);

            Console.WriteLine("{0,-7} =>  {1}", value, result);

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase("0", 0)]
        [TestCase("5", 5)]
        [TestCase("5.678", 5.678f)]
        [TestCase("5.6789012345", 5.6789012f)]
        [TestCase("5.67890123456789012345", 5.6789012f)]
        [TestCase("5.67x", 0)]
        [TestCase("5.67 Da", 0)]
        [TestCase("-0", 0)]
        [TestCase("-5", -5)]
        [TestCase("-5.678", -5.678f)]
        [TestCase("-5.67x", 0)]
        [TestCase("-5.67 Da", 0)]
        [TestCase("X 5.67", 0)]
        public void TestCSngSafe(string value, float expectedResult)
        {
            var result = StringUtilities.CSngSafe(value, 0);

            Console.WriteLine("{0,-7} =>  {1}", value, result);

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase('a', true)]
        [TestCase('m', true)]
        [TestCase('z', true)]
        [TestCase('A', true)]
        [TestCase('M', true)]
        [TestCase('Z', true)]
        [TestCase('.', false)]
        [TestCase('&', false)]
        [TestCase('#', false)]
        [TestCase('\n', false)]
        [TestCase('\r', false)]
        [TestCase('\0', false)]
        public void TestIsLetterAtoZ(char value, bool expectedResult)
        {
            var result = StringUtilities.IsLetterAtoZ(value);

            Console.WriteLine("{0,-2} =>  {1}", value, result);

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase("", "")]
        [TestCase(" ", " ")]
        [TestCase("dataset1", "dataset1")]
        [TestCase("dataset", "dataset1", "dataset2")]
        [TestCase("dataset", "dataset1", "dataset2", "dataset3")]
        [TestCase("data", "dataset1", "data2", "dataset3")]
        [TestCase("data", "dataset1", "data2", "dataset3", "data2", "dataset1")]
        [TestCase("da", "dataset", "data", "daily")]
        [TestCase("da", "data", "dataset", "daily")]
        [TestCase("da", "daily", "data", "dataset")]
        [TestCase("da", "daily", "dataset", "data")]
        [TestCase("d", "dataset", "data", "diamond")]
        [TestCase("", "dataset", "crystal", "daily")]
        public void TestLongestCommonStringFromStart(string expectedResult, params string[] itemNames)
        {
            var items = new List<string>();

            foreach (var item in itemNames)
            {
                if (item.Length > 0)
                    items.Add(item);
            }

            var result = StringUtilities.LongestCommonStringFromStart(items);

            Console.WriteLine("Longest common string: ");
            Console.WriteLine(result);
            Console.WriteLine();

            foreach (var item in items)
            {
                Console.WriteLine(item);
            }

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase("", "")]
        [TestCase(" ", " ")]
        [TestCase("dataset1", "dataset1")]
        [TestCase("dataset", "dataset1", "dataset2")]
        [TestCase("dataset", "dataset1", "dataset2", "dataset3")]
        [TestCase("data", "dataset1", "data2", "dataset3")]
        [TestCase("data", "dataset1", "data2", "dataset3", "data2", "dataset1")]
        [TestCase("da", "dataset", "data", "daily")]
        [TestCase("da", "data", "dataset", "daily")]
        [TestCase("da", "daily", "data", "dataset")]
        [TestCase("da", "daily", "dataset", "data")]
        [TestCase("d", "dataset", "data", "diamond")]
        [TestCase("", "dataset", "crystal", "daily")]
        public void TestLongestCommonStringFromStartVerbose(string expectedResult, params string[] itemNames)
        {
            var items = new List<string>();

            foreach (var item in itemNames)
            {
                if (item.Length > 0)
                    items.Add(item);
            }

            var result = StringUtilities.LongestCommonStringFromStartVerbose(items);

            Console.WriteLine("Longest common string: ");
            Console.WriteLine(result);
            Console.WriteLine();

            foreach (var item in items)
            {
                Console.WriteLine(item);
            }

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase(0, "0")]
        [TestCase(0.1, "0.1")]
        [TestCase(0.01, "0.01")]
        [TestCase(0.001, "0.001")]
        [TestCase(0.0001, "0.0001")]
        [TestCase(0.00001, "0.00001")]
        [TestCase(0.000001, "0.000001")]
        [TestCase(0.0000001, "0")]
        [TestCase(0.00000001, "0")]
        [TestCase(0.83, "0.83")]
        [TestCase(3.83, "3.83")]
        [TestCase(-0.83, "-0.83")]
        [TestCase(-4.83, "-4.83")]
        public void TestMassErrorToString(double massErrorDa, string expectedResult)
        {
            var result = StringUtilities.MassErrorToString(massErrorDa);

            Console.WriteLine("{0,-10} =>  {1}", massErrorDa, result);

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase(0, "Dataset", "Dataset")]
        [TestCase(1, "Dataset", "Dataset")]
        [TestCase(2, "Dataset", "Dataset")]
        [TestCase(3, "Dataset", "Dataset")]
        [TestCase(0, "0", "0")]
        [TestCase(1, "0", "0")]
        [TestCase(2, "0", "0")]
        [TestCase(3, "0", "0")]
        [TestCase(0, "0.0", "0.0")]
        [TestCase(1, "0.0", "0.0")]
        [TestCase(2, "0.0", "0")]
        [TestCase(3, "0.0", "0")]
        [TestCase(0, "5.0", "5.0")]
        [TestCase(1, "5.0", "5.0")]
        [TestCase(2, "5.0", "5.0")]
        [TestCase(3, "5.0", "5.0")]
        [TestCase(0, "5.1", "5.1")]
        [TestCase(1, "5.1", "5.1")]
        [TestCase(2, "5.1", "5.1")]
        [TestCase(3, "5.1", "5.1")]
        public void TestTrimZeroIfNotFirstID(int resultID, string value, string expectedResult)
        {
            var result = StringUtilities.TrimZeroIfNotFirstID(resultID, value);

            Console.WriteLine("{0}: {1,-7} =>  {2}", resultID, value, result);

            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        [TestCase("Dataset", "Dataset")]
        [TestCase("5", "5")]
        [TestCase("5.0", "5.0")]
        [TestCase("5.00", "5.00")]
        [TestCase("-5.00", "-5.00")]
        [TestCase("0", "0")]
        [TestCase("0.0", "0")]
        [TestCase("0.00", "0.00")]
        [TestCase("-0", "-0")]
        [TestCase("-0.0", "-0.0")]
        [TestCase("-0.00", "-0.00")]
        public void TestTrimZero(string value, string expectedResult)
        {
            var result = StringUtilities.TrimZero(value);

            Console.WriteLine("{0,-7} =>  {1}", value, result);

            Assert.AreEqual(expectedResult, result);
        }
    }
}
