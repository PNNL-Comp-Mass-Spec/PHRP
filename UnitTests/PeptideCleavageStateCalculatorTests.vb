﻿Imports NUnit.Framework
Imports PHRPReader

<TestFixture()>
Public Class PeptideCleavageStateCalculatorTests

    Private mCleavageStateCalculator As clsPeptideCleavageStateCalculator

    <SetUp()>
    Public Sub Init()
        mCleavageStateCalculator = New clsPeptideCleavageStateCalculator()
    End Sub

    <Test()>
    <TestCase("A.BCDE.F", "BCDE", "A", "F")>
    <TestCase("-.BCDE.F", "BCDE", "-", "F")>
    <TestCase("A.BCDE.-", "BCDE", "A", "-")>
    <TestCase("A.B.F", "B", "A", "F")>
    <TestCase("A.BCDE", "BCDE", "A", "")>
    <TestCase("BCDE.F", "BCDE", "", "F")>
    <TestCase("FA.BCDE.FG", "BCDE", "FA", "FG")>
    <TestCase("BCDE", "BCDE", "", "")>
    <TestCase("BCDE.", "BCDE", "", "")>
    <TestCase(".BCDE", "BCDE", "", "")>
    <TestCase(".F.", "F", "", "")>
    <TestCase("F..E", "", "F", "E")>
    <TestCase("AF..EF", "AF..EF", "", "")>
    <TestCase("AFF..EF", "AFF..EF", "", "")>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", "TGMLTQKFARSLGMLAVDNQARV", "E", "")>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", "TGMLTQKFARSLGMLAVDNQARV", "", "R")>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV..", "TGMLTQKFARSLGMLAVDNQARV", "", "")>
    <TestCase("A.BCDEFGHIJK.L", "BCDEFGHIJK", "A", "L")>
    <TestCase("A.B*CDEFGHIJK.L", "B*CDEFGHIJK", "A", "L")>
    <TestCase("A.BCDEFGHIJK*.L", "BCDEFGHIJK*", "A", "L")>
    <TestCase("A.BCDEFGHIJK.L", "BCDEFGHIJK", "A", "L")>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", "TGMLTQKFARSLGMLAVDNQARV", "E", "")>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", "TGMLTQKFARSLGMLAVDNQARV", "", "R")>
    Public Sub TestSplitPrefixAndSuffix(sequence As String, expectedPrimarySeq As String, expectedPrefix As String, expectedSuffix As String)

        Dim strPrimarySequence As String = String.Empty
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(sequence, strPrimarySequence, strPrefix, strSuffix)

        Assert.AreEqual(expectedPrimarySeq, strPrimarySequence)
        Assert.AreEqual(expectedPrefix, strPrefix)
        Assert.AreEqual(expectedSuffix, strSuffix)

    End Sub

    <Test()>
    <TestCase("A.BCDE.F", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific)>
    <TestCase("-.BCDE.F", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific)>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific)>
    <TestCase("..TGMLTQKFARSLGMLAVDNQARV.R", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific)>
    <TestCase("..TGMLTQKFARSLGMLAVDNQAR.F", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full)>
    <TestCase("A.BCDEFGHIJK.L", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("A.B*CDEFGHIJK.L", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("A.BCDEFGHIJK.L", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("-.GLMVPVIR.A", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full)>
    <TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full)>
    <TestCase("K.MSSTFIGNSTAIQELFK.R", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full)>
    <TestCase("K.MSSTFIGNSTAIQELFK.P", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("K.MSSTFIGNSTAIQELFR.P", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("K.PSSTFIGNSTAIQELFR.D", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("K.PSSTFIGNSTAIQELFR.P", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific)>
    <TestCase("K.RSSTFIGNSTAIQELFK.R", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full)>
    <TestCase("K.MSSTFIGNSTAIQELFD.R", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    <TestCase("F.MSSTFIGNSTAIQELFK.R", clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial)>
    Public Sub TestComputeCleavageState(sequence As String, expectedCleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants)

        Dim cleavageState = mCleavageStateCalculator.ComputeCleavageState(sequence)

        Assert.AreEqual(expectedCleavageState, cleavageState)

    End Sub

    <Test()>
    <TestCase("A.BCDE.F", 0)>
    <TestCase("-.BCDE.F", 0)>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", 3)>
    <TestCase("F.TGMLTQKFARSLGMLAVDNQARV.R", 3)>
    <TestCase("G.TGMLTQKFARSLGMLAVDNQAR.F", 2)>
    <TestCase("K.TGMLTQKFARSLGMKLAVDNQARV.R", 4)>
    <TestCase("..TGMLTQKFARSLGMKLAVDNQARV.R", 4)>
    <TestCase("-.TGMLTQKFARSLGMKLAVDNQARV.R", 4)>
    <TestCase("M.TGMLTQKFARSLGMKPLAVDNQARV.R", 3)>
    <TestCase("A.BCDEFGHIJK.L", 0)>
    <TestCase("A.B*CDEFGHIJK.L", 0)>
    <TestCase("A.BCDEFGHIJK.L", 0)>
    <TestCase("-.GLMVPVIR.A", 0)>
    <TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", 1)>
    <TestCase("K.MSSTFIGNSTAIQELFK.R", 0)>
    Public Sub TestComputeNumberOfMissedCleavages(sequence As String, expectedMissedCleavages As Short)

        Dim missedCleavages = mCleavageStateCalculator.ComputeNumberOfMissedCleavages(sequence)

        Assert.AreEqual(expectedMissedCleavages, missedCleavages)

    End Sub

    <Test()>
    <TestCase("A.BCDE.F", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)>
    <TestCase("-.BCDE.F", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus)>
    <TestCase("E.TGMLTQKFARSLGMLAVDNQARV..", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus)>
    <TestCase("K.TGMLTQKFARSLGMKLAVDNQARV.R", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)>
    <TestCase("..TGMLTQKFARSLGMKLAVDNQARV.R", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus)>
    <TestCase("-.TGMLTQKFARSLGMKLAVDNQARV.R", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus)>
    <TestCase("M.TGMLTQKFARSLGMKPLAVDNQARV.R", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)>
    <TestCase("A.BCDEFGHIJK.L", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)>
    <TestCase("A.B*CDEFGHIJK.L", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)>
    <TestCase("-.GLMVPVIR.A", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus)>
    <TestCase("F.GLMVPVIR.-", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus)>
    <TestCase("-.GLMVPVIR.-", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus)>
    <TestCase("R.EVSSRPS+79.966T+79.966PGLSVVSGISATSEDIPNKIEDLR.S", clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)>
    Public Sub TestComputeTerminusState(sequence As String, expectedTerminusState As clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants)

        Dim terminusState = mCleavageStateCalculator.ComputeTerminusState(sequence)

        Assert.AreEqual(expectedTerminusState, terminusState)
    End Sub
End Class