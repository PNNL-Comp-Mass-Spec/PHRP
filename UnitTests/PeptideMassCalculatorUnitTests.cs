using System;
using System.Collections.Generic;
using NUnit.Framework;
using PHRPReader;

namespace PHRP_UnitTests
{

    [TestFixture()]
    public class PeptideMassCalculatorUnitTests
    {
        private clsPeptideMassCalculator mPeptideMassCalculator;
        [SetUp()]
        public void Init()
        {
            mPeptideMassCalculator = new clsPeptideMassCalculator();
        }

        [Test()]
        [TestCase("A.LCDE.F", 478.17333)]
        [TestCase("-.LCDE.F", 478.17333)]
        [TestCase("A.LCDE.-", 478.17333)]
        [TestCase("A.LCDE", 478.17333)]
        [TestCase("LCDE.F", 478.17333)]
        [TestCase("FA.LCDE.FG", 478.17333)]
        [TestCase("LCDE", 478.17333)]
        [TestCase("LCDE.", 478.17333)]
        [TestCase(".LCDE", 478.17333)]
        [TestCase("A.L.F", 131.09462)]
        [TestCase(".F.", 165.07897)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", 2506.3147)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", 2506.3147)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV..", 2506.3147)]
        [TestCase("A.LCDEFGHIJK.L", 1060.50112)]
        [TestCase("A.VCDEFGHIJK.L", 1046.48547)]
        [TestCase("A.M*CDEFGHIJK.L", -1)]
        [TestCase("A.SCDEFGHIJK*.L", -1)]
        public void TestComputeAminoAcidMass(string strSequence, double expectedMass)
        {
            var computedMass = mPeptideMassCalculator.ComputeSequenceMass(strSequence);

            Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strSequence, computedMass, expectedMass);

            Assert.AreEqual(expectedMass, computedMass, 0.001, "Unexpected mass for the amino acid sequence");
        }

        [Test()]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV.R", "3:15.994; 11:79.996", 2602.3047)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", "3:15.994; 11:79.996", 2602.3047)]
        [TestCase("..TGMLTQKFARSLGMLAVDNQARV..", "3:15.994; 11:79.996", 2602.3047)]
        [TestCase("E.TGMLTQKFARSLGMLAVDNQARV.R.", "3:15.994; 11:79.996", -1)]
        [TestCase("A.LCDEFGHIJK.L", "6:32.5", 1093.00112)]
        [TestCase("A.MCDEFGHIJK.L", "1:15.9994", 1094.45694)]
        [TestCase("A.M*CDEFGHIJK.L", "1:15.9994; 9:138.5", -1)]
        [TestCase("A.SCDEFGHIJK.L", "3:24.3", 1058.74908)]
        public void TestComputeAminoAcidMassModifiedResidues(string strSequence, string residueModList, double expectedMass)
        {
            var residueMods = residueModList.Split(';');

            var modifiedResidues = new List<clsPeptideMassCalculator.udtPeptideSequenceModInfoType>();

            foreach (var modifiedResidue in residueMods)
            {
                var modParts = modifiedResidue.Split(':');

                var residueLocation = int.Parse(modParts[0]);
                var modMass = double.Parse(modParts[1]);

                var modInfo = new clsPeptideMassCalculator.udtPeptideSequenceModInfoType
                {
                    ResidueLocInPeptide = residueLocation,
                    ModificationMass = modMass,
                    AffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
                };

                modifiedResidues.Add(modInfo);
            }

            var computedMass = mPeptideMassCalculator.ComputeSequenceMass(strSequence, modifiedResidues);

            Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strSequence, computedMass, expectedMass);

            Assert.AreEqual(expectedMass, computedMass, 0.0001, "Unexpected mass for the amino acid sequence");
        }

        [Test()]
        [TestCase("R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A", 2184.97443)]
        [TestCase("K.Q-17.0265QIEESTSDYDKEK.L", 1681.73186)]
        [TestCase("Q-17.0265QIEESTSDYDKEK", 1681.73186)]
        public void TestComputeAminoAcidMassNumericMods(string strSequence, double expectedMass)
        {
            var computedMass = mPeptideMassCalculator.ComputeSequenceMassNumericMods(strSequence);

            Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strSequence, computedMass, expectedMass);

            Assert.AreEqual(expectedMass, computedMass, 1E-05, "Unexpected mass for the empirical formula");
        }

        [Test()]
        [TestCase("CHNOSP", 105.9516)]
        [TestCase("C2H3N4OS", 131.0027)]
        [TestCase("C2H3N-2OS3N+3S-2", 88.99353)]
        [TestCase("CH53N-3Se4P2Xe3", 513.03752)]
        public void TestComputeEmpiricalFormulaMassMethod1(string strEmpiricalFormula, double expectedMass)
        {
            var empiricalFormula = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(strEmpiricalFormula);
            var computedMass = clsPeptideMassCalculator.ComputeMonoistopicMass(empiricalFormula);

            List<string> unknownSymbols = null;
            var computedMassAlt = clsPeptideMassCalculator.ComputeMonoistopicMass(empiricalFormula.ElementCounts, out unknownSymbols);

            Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strEmpiricalFormula, computedMass, expectedMass);

            Assert.AreEqual(computedMass, computedMassAlt, 1E-05, "The two overloads ComputeMonoistopicMass reported conflicting mass values");

            Assert.AreEqual(expectedMass, computedMass, 0.0001, "Unexpected mass for the empirical formula");
        }

        [Test()]
        [TestCase(1000, 1, 0, 998.99272)]
        [TestCase(1000, 3, 1, 2997.98545)]
        [TestCase(700, 0, 2, 351.00728)]
        [TestCase(2400, 2, 5, 960.60437)]
        public void TestConvoluteMass(double massMz, int currentCharge, int newCharge, double expectedMz)
        {
            var newMz = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge);
            var newMzAlt = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge, clsPeptideMassCalculator.MASS_PROTON);

            Console.WriteLine("{0} from {1}+ to {2}+ is {3:F5}; expected {4:F5}", massMz, currentCharge, newCharge, newMz, expectedMz);

            Assert.AreEqual(newMz, newMzAlt, 1E-05, "The two overloads of ConvoluteMass reported conflicting mass values");

            Assert.AreEqual(expectedMz, newMz, 0.0001, "Unexpected convoluted m/z");

        }

        [Test()]
        [TestCase(1000, 1, 0, 977.01078)]
        [TestCase(1000, 3, 1, 2954.02156)]
        [TestCase(700, 0, 2, 372.9892)]
        [TestCase(2400, 2, 5, 973.79353)]
        public void TestConvoluteMassNa(double massMz, int currentCharge, int newCharge, double expectedMz)
        {
            const double MASS_SODIUM = 22.98922189;

            mPeptideMassCalculator.ChargeCarrierMass = MASS_SODIUM;

            var newMz = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge);
            var newMzAlt = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge, MASS_SODIUM);

            Console.WriteLine("{0} from {1}+ to {2}+ is {3:F5}; expected {4:F5}", massMz, currentCharge, newCharge, newMz, expectedMz);

            Assert.AreEqual(newMz, newMzAlt, 1E-05, "The two overloads of ConvoluteMass reported conflicting mass values");

            Assert.AreEqual(expectedMz, newMz, 0.0001, "Unexpected convoluted m/z");
        }

        [Test()]
        [TestCase('A', 71.03711)]
        [TestCase('B', 114.04292)]
        [TestCase('C', 103.00918)]
        [TestCase('D', 115.02694)]
        [TestCase('E', 129.04259)]
        [TestCase('F', 147.06841)]
        [TestCase('G', 57.02146)]
        [TestCase('H', 137.0589)]
        [TestCase('I', 113.08406)]
        [TestCase('J', 0)]
        [TestCase('K', 128.09496)]
        [TestCase('L', 113.08406)]
        [TestCase('M', 131.04048)]
        [TestCase('N', 114.04292)]
        [TestCase('O', 114.07931)]
        [TestCase('P', 97.05276)]
        [TestCase('Q', 128.05857)]
        [TestCase('R', 156.1011)]
        [TestCase('S', 87.03202)]
        [TestCase('T', 101.04767)]
        [TestCase('U', 150.95363)]
        [TestCase('V', 99.06841)]
        [TestCase('W', 186.07931)]
        [TestCase('X', 113.08406)]
        [TestCase('Y', 163.06332)]
        [TestCase('Z', 128.05857)]
        public void TestGetAminoAcidMass(char aminoAcidSymbol, double expectedMass)
        {
            var computedMass = mPeptideMassCalculator.GetAminoAcidMass(aminoAcidSymbol);

            var empiricalFormula = new clsEmpiricalFormula();
            empiricalFormula.AddElements(mPeptideMassCalculator.GetAminoAcidEmpiricalFormula(aminoAcidSymbol));

            var computedMassAlt = clsPeptideMassCalculator.ComputeMonoistopicMass(empiricalFormula);

            Assert.AreEqual(expectedMass, computedMass, 0.0001, "Amino acid does not match the expected value");

            Assert.AreEqual(computedMass, computedMassAlt, 1E-05, "GetAminoAcidMass and ComputeMonoistopicMass do not agree on the mass for the amino acid");
        }

        [Test()]
        [TestCase(0.005, 1000, 5)]
        [TestCase(0.015, 1000, 15)]
        [TestCase(0.005, 500, 10)]
        [TestCase(0.001, 500, 2)]
        [TestCase(0.005, 2300, 2.17391)]
        [TestCase(0.0012, 2300, 0.52174)]
        public void TestPPMConversion(double massToConvert, double currentMz, double expectedPPM)
        {
            var convertedPPM = clsPeptideMassCalculator.MassToPPM(massToConvert, currentMz);

            var reconvertedMass = clsPeptideMassCalculator.PPMToMass(convertedPPM, currentMz);

            Console.WriteLine("DelM of {0} Da converts to {1:F5} ppm at {2:F5} m/z ", massToConvert, convertedPPM, currentMz);

            Assert.AreEqual(expectedPPM, convertedPPM, 0.0001, "PPM Conversion error");

            Assert.AreEqual(massToConvert, reconvertedMass, 1E-05, "Da to PPM to Da round trip error");
        }

        [Test()]
        [TestCase(1200, 2, 1198.99272, 600.5036)]
        [TestCase(1200, 3, 1198.99272, 400.67152)]
        public void TestMassConversions(double dblMH, int chargeState, double expectedMonoMass, double expectedMz)
        {
            var computedMonoMass = mPeptideMassCalculator.MHToMonoisotopicMass(dblMH);

            var computedMz = mPeptideMassCalculator.MonoisotopicMassToMZ(computedMonoMass, chargeState);

            Console.WriteLine("{0} MH converts to {1:F5} Da and {2:F5} m/z at charge {3}", dblMH, computedMonoMass, computedMz, chargeState);

            Assert.AreEqual(expectedMonoMass, computedMonoMass, 0.0001, "Monoisotopic mass mismatch");
            Assert.AreEqual(expectedMz, computedMz, 0.0001, "M/Z mismatch");
        }
    }
}
