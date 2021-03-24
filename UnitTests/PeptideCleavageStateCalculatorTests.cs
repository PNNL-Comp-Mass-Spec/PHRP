using NUnit.Framework;
using PHRPReader;
using static PHRPReader.PeptideCleavageStateCalculator;

// ReSharper disable StringLiteralTypo
namespace PHRP_UnitTests
{
    [TestFixture]
    public class PeptideCleavageStateCalculatorTests
    {
        private PeptideCleavageStateCalculator mCleavageStateCalculator;

        [SetUp]
        public void Init()
        {
            mCleavageStateCalculator = new PeptideCleavageStateCalculator();
        }

        [Test]
        [TestCase("A.BCDE.F", "BCDE", "A", "F")]
        [TestCase("-.BCDE.F", "BCDE", "-", "F")]
        [TestCase("A.BCDE.-", "BCDE", "A", "-")]
        [TestCase("A.B.F", "B", "A", "F")]
        [TestCase("A.BCDE", "BCDE", "A", "")]
        [TestCase("BCDE.F", "BCDE", "", "F")]
        [TestCase("FA.BCDE.FG", "BCDE", "FA", "FG")]
        [TestCase("BCDE", "BCDE", "", "")]
        [TestCase("BCDE.", "BCDE", "", "")]
        [TestCase(".BCDE", "BCDE", "", "")]
        [TestCase(".F.", "F", "", "")]
        [TestCase("F..E", "", "F", "E")]
        [TestCase("AF..EF", "AF..EF", "", "")]
        [TestCase("AFF..EF", "AFF..EF", "", "")]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", "TGMLTQKFARSLGMLAVDNQARV", "E", "")]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", "TGMLTQKFARSLGMLAVDNQARV", "", "R")]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV..", "TGMLTQKFARSLGMLAVDNQARV", "", "")]
        [TestCase("A.BCDEFGHIJK.L", "BCDEFGHIJK", "A", "L")]
        [TestCase("A.B*CDEFGHIJK.L", "B*CDEFGHIJK", "A", "L")]
        [TestCase("A.BCDEFGHIJK*.L", "BCDEFGHIJK*", "A", "L")]
        [TestCase("A.BCDEFGHIJK.L", "BCDEFGHIJK", "A", "L")]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", "TGMLTQKFARSLGMLAVDNQARV", "E", "")]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", "TGMLTQKFARSLGMLAVDNQARV", "", "R")]
        public void TestSplitPrefixAndSuffix(string sequence, string expectedPrimarySeq, string expectedPrefix, string expectedSuffix)
        {
            PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequence, out var primarySequence, out var prefix, out var suffix);

            Assert.AreEqual(expectedPrimarySeq, primarySequence);
            Assert.AreEqual(expectedPrefix, prefix);
            Assert.AreEqual(expectedSuffix, suffix);
        }

        [Test]
        [TestCase("A.BCDE.F", PeptideCleavageState.NonSpecific)]
        [TestCase("-.BCDE.F", PeptideCleavageState.NonSpecific)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", PeptideCleavageState.NonSpecific)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", PeptideCleavageState.NonSpecific)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQAR.F", PeptideCleavageState.Full)]
        [TestCase("A.BCDEFGHIJK.L", PeptideCleavageState.Partial)]
        [TestCase("A.B*CDEFGHIJK.L", PeptideCleavageState.Partial)]
        [TestCase("A.BCDEFGHIJK.L", PeptideCleavageState.Partial)]
        [TestCase("-.GLMVPVIR.A", PeptideCleavageState.Full)]
        [TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", PeptideCleavageState.Full)]
        [TestCase("K.MSSTFIGNSTAIQELFK.R", PeptideCleavageState.Full)]
        [TestCase("K.MSSTFIGNSTAIQELFK.P", PeptideCleavageState.Partial)]
        [TestCase("K.MSSTFIGNSTAIQELFR.P", PeptideCleavageState.Partial)]
        [TestCase("K.PSSTFIGNSTAIQELFR.D", PeptideCleavageState.Partial)]
        [TestCase("K.PSSTFIGNSTAIQELFR.P", PeptideCleavageState.NonSpecific)]
        [TestCase("K.RSSTFIGNSTAIQELFK.R", PeptideCleavageState.Full)]
        [TestCase("K.MSSTFIGNSTAIQELFD.R", PeptideCleavageState.Partial)]
        [TestCase("F.MSSTFIGNSTAIQELFK.R", PeptideCleavageState.Partial)]
        [TestCase("K.ACDEFGR.S", PeptideCleavageState.Full)]
        [TestCase("R.ACDEFGR.S", PeptideCleavageState.Full)]
        [TestCase("-.ACDEFGR.S", PeptideCleavageState.Full)]
        [TestCase("R.ACDEFGH.-", PeptideCleavageState.Full)]
        [TestCase("-.ACDEFG.-", PeptideCleavageState.Full)]
        [TestCase("-.ACDEFG.-", PeptideCleavageState.Full)]
        [TestCase("K.ACDEFGR*.S", PeptideCleavageState.Full)]
        [TestCase("-.ACDEFGR*.S", PeptideCleavageState.Full)]
        [TestCase("K.ACDEFGH.S", PeptideCleavageState.Partial)]
        [TestCase("L.ACDEFGR.S", PeptideCleavageState.Partial)]
        [TestCase("K.ACDEFGR.P", PeptideCleavageState.Partial)]
        [TestCase("K.PCDEFGR.S", PeptideCleavageState.Partial)]
        [TestCase("L.ACDEFGH.S", PeptideCleavageState.NonSpecific)]
        [TestCase("-.ACDEFGH.S", PeptideCleavageState.NonSpecific)]
        [TestCase("L.ACDEFGH.-", PeptideCleavageState.NonSpecific)]
        [TestCase("L.ACDEFGR.P", PeptideCleavageState.NonSpecific)]
        [TestCase("K.PCDEFGR.P", PeptideCleavageState.NonSpecific)]
        public void TestComputeCleavageState(string sequence, PeptideCleavageState expectedCleavageState)
        {
            var cleavageState = mCleavageStateCalculator.ComputeCleavageState(sequence);

            Assert.AreEqual(expectedCleavageState, cleavageState);
        }

        [Test]
        [TestCase("A.BCDE.F", 0)]
        [TestCase("-.BCDE.F", 0)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", 3)]
        [TestCase("F.TGMLTQKFARSLGMLAVDNQARV.R", 3)]
        [TestCase("G.TGMLTQKFARSLGMLAVDNQAR.F", 2)]
        [TestCase("K.TGMLTQKFARSLGMKLAVDNQARV.R", 4)]
        [TestCase("..TGMLTQKFARSLGMKLAVDNQARV.R", 4)]
        [TestCase("-.TGMLTQKFARSLGMKLAVDNQARV.R", 4)]
        [TestCase("M.TGMLTQKFARSLGMKPLAVDNQARV.R", 3)]
        [TestCase("A.BCDEFGHIJK.L", 0)]
        [TestCase("A.B*CDEFGHIJK.L", 0)]
        [TestCase("A.BCDEFGHIJK.L", 0)]
        [TestCase("-.GLMVPVIR.A", 0)]
        [TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", 1)]
        [TestCase("K.MSSTFIGNSTAIQELFK.R", 0)]
        public void TestComputeNumberOfMissedCleavages(string sequence, short expectedMissedCleavages)
        {
            var missedCleavages = mCleavageStateCalculator.ComputeNumberOfMissedCleavages(sequence);

            Assert.AreEqual(expectedMissedCleavages, missedCleavages);
        }

        [Test]
        [TestCase("A.BCDE.F", PeptideTerminusState.None)]
        [TestCase("-.BCDE.F", PeptideTerminusState.ProteinNTerminus)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", PeptideTerminusState.ProteinCTerminus)]
        [TestCase("K.TGMLTQKFARSLGMKLAVDNQARV.R", PeptideTerminusState.None)]
        [TestCase("..TGMLTQKFARSLGMKLAVDNQARV.R", PeptideTerminusState.ProteinNTerminus)]
        [TestCase("-.TGMLTQKFARSLGMKLAVDNQARV.R", PeptideTerminusState.ProteinNTerminus)]
        [TestCase("M.TGMLTQKFARSLGMKPLAVDNQARV.R", PeptideTerminusState.None)]
        [TestCase("A.BCDEFGHIJK.L", PeptideTerminusState.None)]
        [TestCase("A.B*CDEFGHIJK.L", PeptideTerminusState.None)]
        [TestCase("-.GLMVPVIR.A", PeptideTerminusState.ProteinNTerminus)]
        [TestCase("F.GLMVPVIR.-", PeptideTerminusState.ProteinCTerminus)]
        [TestCase("-.GLMVPVIR.-", PeptideTerminusState.ProteinNandCCTerminus)]
        [TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", PeptideTerminusState.None)]
        public void TestComputeTerminusState(string sequence, PeptideTerminusState expectedTerminusState)
        {
            var terminusState = mCleavageStateCalculator.ComputeTerminusState(sequence);

            Assert.AreEqual(expectedTerminusState, terminusState);
        }
    }
}
