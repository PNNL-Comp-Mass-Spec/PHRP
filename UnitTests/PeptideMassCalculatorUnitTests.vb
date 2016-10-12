Imports NUnit.Framework
Imports PHRPReader

<TestFixture()>
Public Class PeptideMassCalculatorUnitTests

    Private mPeptideMassCalculator As clsPeptideMassCalculator

    <SetUp()>
    Public Sub Init()
        mPeptideMassCalculator = New clsPeptideMassCalculator()
    End Sub

    <Test()>
    <TestCase("A.LCDE.F", 478.17333)>
    <TestCase("-.LCDE.F", 478.17333)>
    <TestCase("A.LCDE.-", 478.17333)>
    <TestCase("A.L.F", 131.09462)>
    <TestCase("A.LCDE", 478.17333)>
    <TestCase("LCDE.F", 478.17333)>
    <TestCase("FA.LCDE.FG", 478.17333)>
    <TestCase("LCDE", 478.17333)>
    <TestCase("LCDE.", 478.17333)>
    <TestCase(".LCDE", 478.17333)>
    <TestCase(".F.", 165.07897)>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", 2506.3147)>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", 2506.3147)>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV..", 2506.3147)>
    <TestCase("A.LCDEFGHIJK.L", 1060.50112)>
    <TestCase("A.M*CDEFGHIJK.L", -1)>
    <TestCase("A.SCDEFGHIJK*.L", -1)>
    <TestCase("A.VCDEFGHIJK.L", 1046.48547)>
    Public Sub TestComputeAminoAcidMass(strSequence As String, expectedMass As Double)
        Dim computedMass = mPeptideMassCalculator.ComputeSequenceMass(strSequence)

        Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strSequence, computedMass, expectedMass)

        Assert.AreEqual(expectedMass, computedMass, 0.001, "Unexpected mass for the amino acid sequence")
    End Sub

    <Test()>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV.R", "3:15.994; 11:79.996", 2602.3047)>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV.R.", "3:15.994; 11:79.996", -1)>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", "3:15.994; 11:79.996", 2602.3047)>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV..", "3:15.994; 11:79.996", 2602.3047)>
    <TestCase("A.LCDEFGHIJK.L", "6:32.5", 1093.00112)>
    <TestCase("A.MCDEFGHIJK.L", "1:15.9994", 1094.45694)>
    <TestCase("A.M*CDEFGHIJK.L", "1:15.9994; 9:138.5", -1)>
    <TestCase("A.SCDEFGHIJK.L", "3:24.3", 1058.74908)>
    Public Sub TestComputeAminoAcidMassModifiedResidues(
      strSequence As String, residueModList As String, expectedMass As Double)

        Dim residueMods = residueModList.Split(";"c)

        Dim modifiedResidues = New List(Of clsPeptideMassCalculator.udtPeptideSequenceModInfoType)

        For Each modifiedResidue In residueMods
            Dim modParts = modifiedResidue.Split(":"c)

            Dim residueLocation = Integer.Parse(modParts(0))
            Dim modMass = Double.Parse(modParts(1))

            Dim modInfo = New clsPeptideMassCalculator.udtPeptideSequenceModInfoType() With {
                .ResidueLocInPeptide = residueLocation,
                .ModificationMass = modMass,
                .AffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
                }

            modifiedResidues.Add(modInfo)
        Next

        Dim computedMass = mPeptideMassCalculator.ComputeSequenceMass(strSequence, modifiedResidues)

        Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strSequence, computedMass, expectedMass)

        Assert.AreEqual(expectedMass, computedMass, 0.001, "Unexpected mass for the amino acid sequence")
    End Sub

    <Test()>
    <TestCase("R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A", 2184.97443)>
    <TestCase("K.Q-17.0265QIEESTSDYDKEK.L", 1681.73186)>
    <TestCase("Q-17.0265QIEESTSDYDKEK", 1681.73186)>
    Public Sub TestComputeAminoAcidMassNumericMods(strSequence As String, expectedMass As Double)

        Dim computedMass = mPeptideMassCalculator.ComputeSequenceMassNumericMods(strSequence)

        Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strSequence, computedMass, expectedMass)

        Assert.AreEqual(expectedMass, computedMass, 0.00001, "Unexpected mass for the empirical formula")
    End Sub

    <Test()>
    <TestCase("CHNOSP", 105.9516)>
    <TestCase("C2H3N4OS", 131.0027)>
    <TestCase("C2H3N-2OS3N+3S-2", 88.99353)>
    <TestCase("CH53N-3Se4P2Xe3", 513.03752)>
    Public Sub TestComputeEmpiricalFormulaMassMethod1(strEmpiricalFormula As String, expectedMass As Double)

        Dim empiricalFormula = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(strEmpiricalFormula)
        Dim computedMass = clsPeptideMassCalculator.ComputeMonoistopicMass(empiricalFormula)

        Dim unknownSymbols As List(Of String) = Nothing
        Dim computedMassAlt = clsPeptideMassCalculator.ComputeMonoistopicMass(empiricalFormula.ElementCounts, unknownSymbols)

        Console.WriteLine("{0,-30} is {1:F5}; expected {2:F5}", strEmpiricalFormula, computedMass, expectedMass)

        Assert.AreEqual(computedMass, computedMassAlt, 0.0001, "The two overloads ComputeMonoistopicMass reported conflicting mass values")

        Assert.AreEqual(expectedMass, computedMass, 0.001, "Unexpected mass for the empirical formula")
    End Sub

    <Test()>
    <TestCase(1000, 1, 0, 998.99272)>
    <TestCase(1000, 3, 1, 2997.98545)>
    <TestCase(700, 0, 2, 351.00728)>
    <TestCase(2400, 2, 5, 960.60437)>
    Public Sub TestConvoluteMass(massMz As Double, currentCharge As Integer, newCharge As Integer, expectedMz As Double)

        Dim newMz = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge)
        Dim newMzAlt = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge, clsPeptideMassCalculator.MASS_PROTON)

        Console.WriteLine("{0} from {1}+ to {2}+ is {3:F5}; expected {4:F5}", massMz, currentCharge, newCharge, newMz, expectedMz)

        Assert.AreEqual(newMz, newMzAlt, 0.0001, "The two overloads of ConvoluteMass reported conflicting mass values")

        Assert.AreEqual(expectedMz, newMz, 0.001, "Unexpected convoluted m/z")

    End Sub

    <Test()>
    <TestCase(1000, 1, 0, 977.01078)>
    <TestCase(1000, 3, 1, 2954.02156)>
    <TestCase(700, 0, 2, 372.9892)>
    <TestCase(2400, 2, 5, 973.79353)>
    Public Sub TestConvoluteMassNa(massMz As Double, currentCharge As Integer, newCharge As Integer, expectedMz As Double)
        Const MASS_SODIUM = 22.98922189

        mPeptideMassCalculator.ChargeCarrierMass = MASS_SODIUM

        Dim newMz = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge)
        Dim newMzAlt = mPeptideMassCalculator.ConvoluteMass(massMz, currentCharge, newCharge, MASS_SODIUM)

        Console.WriteLine("{0} from {1}+ to {2}+ is {3:F5}; expected {4:F5}", massMz, currentCharge, newCharge, newMz, expectedMz)

        Assert.AreEqual(newMz, newMzAlt, 0.0001, "The two overloads of ConvoluteMass reported conflicting mass values")

        Assert.AreEqual(expectedMz, newMz, 0.001, "Unexpected convoluted m/z")

    End Sub

    <Test()>
    <TestCase("A", 71.0371100902557)>
    <TestCase("B", 114.042921543121)>
    <TestCase("C", 103.009180784225)>
    <TestCase("D", 115.026938199997)>
    <TestCase("E", 129.042587518692)>
    <TestCase("F", 147.068408727646)>
    <TestCase("G", 57.0214607715607)>
    <TestCase("H", 137.058904886246)>
    <TestCase("I", 113.084058046341)>
    <TestCase("J", 0)>
    <TestCase("K", 128.094955444336)>
    <TestCase("L", 113.084058046341)>
    <TestCase("M", 131.040479421616)>
    <TestCase("N", 114.042921543121)>
    <TestCase("O", 114.079306125641)>
    <TestCase("P", 97.0527594089508)>
    <TestCase("Q", 128.058570861816)>
    <TestCase("R", 156.101100921631)>
    <TestCase("S", 87.0320241451263)>
    <TestCase("T", 101.047673463821)>
    <TestCase("U", 150.95363)>
    <TestCase("V", 99.0684087276459)>
    <TestCase("W", 186.079306125641)>
    <TestCase("X", 113.084058046341)>
    <TestCase("Y", 163.063322782516)>
    <TestCase("Z", 128.058570861816)>
    Public Sub TestGetAminoAcidMass(aminoAcidSymbol As String, expectedMass As Double)

        Dim computedMass = mPeptideMassCalculator.GetAminoAcidMass(aminoAcidSymbol)

        Dim empiricalFormula = New clsEmpiricalFormula()
        empiricalFormula.AddElements(mPeptideMassCalculator.GetAminoAcidEmpiricalFormula(aminoAcidSymbol))

        Dim computedMassAlt = clsPeptideMassCalculator.ComputeMonoistopicMass(empiricalFormula)

        Assert.AreEqual(expectedMass, computedMass, 0.01, "Amino acid does not match the expected value")

        Assert.AreEqual(computedMass, computedMassAlt, 0.01, "GetAminoAcidMass and ComputeMonoistopicMass do not agree on the mass for the amino acid")

    End Sub

    <Test()>
    <TestCase(0.005, 1000, 5)>
    <TestCase(0.015, 1000, 15)>
    <TestCase(0.005, 500, 10)>
    <TestCase(0.001, 500, 2)>
    <TestCase(0.005, 2300, 2.17391)>
    <TestCase(0.0012, 2300, 0.52174)>
    Public Sub TestPPMConversion(massToConvert As Double, currentMz As Double, expectedPPM As Double)

        Dim convertedPPM = clsPeptideMassCalculator.MassToPPM(massToConvert, currentMz)

        Dim reconvertedMass = clsPeptideMassCalculator.PPMToMass(convertedPPM, currentMz)

        Console.WriteLine("DelM of {0} Da converts to {1:F5} ppm at {2:F5} m/z ", massToConvert, convertedPPM, currentMz)

        Assert.AreEqual(expectedPPM, convertedPPM, 0.001, "PPM Conversion error")

        Assert.AreEqual(massToConvert, reconvertedMass, 0.0001, "Da to PPM to Da round trip error")

    End Sub

    <Test()>
    <TestCase(1200, 2, 1198.99272, 600.5036)>
    <TestCase(1200, 3, 1198.99272, 400.67152)>
    Public Sub TestMassConversions(dblMH As Double, chargeState As Double, expectedMonoMass As Double, expectedMz As Double)

        Dim computedMonoMass = mPeptideMassCalculator.MHToMonoisotopicMass(dblMH)

        Dim computedMz = mPeptideMassCalculator.MonoisotopicMassToMZ(computedMonoMass, chargeState)

        Console.WriteLine("{0} MH converts to {1:F5} Da and {2:F5} m/z at charge {3}", dblMH, computedMonoMass, computedMz, chargeState)

        Assert.AreEqual(expectedMonoMass, computedMonoMass, 0.001, "Monoisotopic mass mismatch")
        Assert.AreEqual(expectedMz, computedMz, 0.001, "M/Z mismatch")

    End Sub
End Class