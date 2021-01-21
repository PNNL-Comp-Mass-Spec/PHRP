using NUnit.Framework;
using PHRPReader;

// ReSharper disable StringLiteralTypo
namespace PHRP_UnitTests
{
    [TestFixture]
    public class PeptideCleavageStateCalculatorTests
    {
        private clsPeptideCleavageStateCalculator mCleavageStateCalculator;

        [SetUp]
        public void Init()
        {
            mCleavageStateCalculator = new clsPeptideCleavageStateCalculator();
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
            clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequence, out var primarySequence, out var prefix, out var suffix);

            Assert.AreEqual(expectedPrimarySeq, primarySequence);
            Assert.AreEqual(expectedPrefix, prefix);
            Assert.AreEqual(expectedSuffix, suffix);
        }

        [Test]
        [TestCase("A.BCDE.F", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("-.BCDE.F", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQAR.F", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("A.BCDEFGHIJK.L", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("A.B*CDEFGHIJK.L", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("A.BCDEFGHIJK.L", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("-.GLMVPVIR.A", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("K.MSSTFIGNSTAIQELFK.R", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("K.MSSTFIGNSTAIQELFK.P", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("K.MSSTFIGNSTAIQELFR.P", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("K.PSSTFIGNSTAIQELFR.D", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("K.PSSTFIGNSTAIQELFR.P", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("K.RSSTFIGNSTAIQELFK.R", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("K.MSSTFIGNSTAIQELFD.R", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("F.MSSTFIGNSTAIQELFK.R", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("K.ACDEFGR.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("R.ACDEFGR.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("-.ACDEFGR.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("R.ACDEFGH.-", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("-.ACDEFG.-", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("-.ACDEFG.-", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("K.ACDEFGR*.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("-.ACDEFGR*.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Full)]
        [TestCase("K.ACDEFGH.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("L.ACDEFGR.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("K.ACDEFGR.P", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("K.PCDEFGR.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.Partial)]
        [TestCase("L.ACDEFGH.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("-.ACDEFGH.S", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("L.ACDEFGH.-", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("L.ACDEFGR.P", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        [TestCase("K.PCDEFGR.P", clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants.NonSpecific)]
        public void TestComputeCleavageState(string sequence, clsPeptideCleavageStateCalculator.PeptideCleavageStateConstants expectedCleavageState)
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
        [TestCase("A.BCDE.F", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.None)]
        [TestCase("-.BCDE.F", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinNTerminus)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinCTerminus)]
        [TestCase("K.TGMLTQKFARSLGMKLAVDNQARV.R", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.None)]
        [TestCase("..TGMLTQKFARSLGMKLAVDNQARV.R", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinNTerminus)]
        [TestCase("-.TGMLTQKFARSLGMKLAVDNQARV.R", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinNTerminus)]
        [TestCase("M.TGMLTQKFARSLGMKPLAVDNQARV.R", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.None)]
        [TestCase("A.BCDEFGHIJK.L", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.None)]
        [TestCase("A.B*CDEFGHIJK.L", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.None)]
        [TestCase("-.GLMVPVIR.A", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinNTerminus)]
        [TestCase("F.GLMVPVIR.-", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinCTerminus)]
        [TestCase("-.GLMVPVIR.-", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.ProteinNandCCTerminus)]
        [TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants.None)]
        public void TestComputeTerminusState(string sequence, clsPeptideCleavageStateCalculator.PeptideTerminusStateConstants expectedTerminusState)
        {
            var terminusState = mCleavageStateCalculator.ComputeTerminusState(sequence);

            Assert.AreEqual(expectedTerminusState, terminusState);
        }
    }
}
