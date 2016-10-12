' This class will compute the mass of a given peptide sequence.  The sequence
'  must consist of only capital letters, though if RemovePrefixAndSuffixIfPresent = True, then
'  characters up to the first . and after the last . in the sequence will be removed
' Residue modification information can be supplied by passing an array of modifications using
'  the structure udtPeptideSequenceModInfoType
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 3, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'

Option Strict On

Imports System.Runtime.InteropServices
Imports System.Text.RegularExpressions
Imports System.Text

Public Class clsPeptideMassCalculator

#Region "Constants and Enums"

    Public Const NO_AFFECTED_ATOM_SYMBOL As Char = "-"c

    Public Const MASS_HYDROGEN As Double = 1.0078246
    Public Const MASS_OXYGEN As Double = 15.9949141
    Public Const MASS_PROTON As Double = 1.00727649               ' Note that this is the mass of hydrogen minus the mass of one electron
    Public Const MASS_ELECTRON As Double = 0.00054811

    Public Const DEFAULT_N_TERMINUS_MASS_CHANGE As Double = MASS_HYDROGEN
    Public Const DEFAULT_C_TERMINUS_MASS_CHANGE As Double = MASS_OXYGEN + MASS_HYDROGEN

    Private Const ASCII_VAL_LETTER_A As Byte = 65
#End Region

#Region "Structures"
    Public Structure udtPeptideSequenceModInfoType
        ''' <summary>
        ''' Position that the modification occurs; not used by clsPeptideMassCalculator
        ''' </summary>
        Public ResidueLocInPeptide As Integer

        ''' <summary>
        ''' Modification mass
        ''' </summary>
        Public ModificationMass As Double

        ''' <summary>
        ''' Affected atom
        ''' </summary>
        ''' <remarks>
        ''' Set to Nothing or to NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications)
        ''' For Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
        ''' </remarks>
        Public AffectedAtom As Char
    End Structure

#End Region

#Region "Classwide Variables"
    ' The Amino Acid arrays contain 26 entries, corresponding to A through Z
    ' Invalid/Undefined letters (J and U) have values of 0 for the mass and atom counts
    ' The values can be customized using SetAminoAcidMass and SetAminoAcidAtomCounts

    Private Const AMINO_ACID_LIST_MAX_INDEX As Byte = 25
    Private mAminoAcidMasses() As Double
    Private mAminoAcidEmpiricalFormulas() As clsEmpiricalFormula

    ' Typically mPeptideNTerminusMass + mPeptideCTerminusMass = 18.0105633 (the mass of water)
    Private mPeptideNTerminusMass As Double
    Private mPeptideCTerminusMass As Double

    Private mRemovePrefixAndSuffixIfPresent As Boolean

    Private mErrorMessage As String

    ''' <summary>
    ''' Regular expression for parsing an empirical formula
    ''' </summary>
    Private Shared ReadOnly mAtomicFormulaRegEx As Regex

    ''' <summary>
    ''' This dictionary tracks element symbols and monoisotopic masses
    ''' </summary>
    Private Shared ReadOnly mElementMonoMasses As Dictionary(Of String, Double)

#End Region

#Region "Properties"

    Public Property ChargeCarrierMass As Double

    ''' <summary>
    ''' Most recent error message
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ErrorMessage() As String
        Get
            Return mErrorMessage
        End Get
    End Property

    Public Property PeptideCTerminusMass() As Double
        Get
            Return mPeptideCTerminusMass
        End Get
        Set(Value As Double)
            mPeptideCTerminusMass = Value
        End Set
    End Property

    Public Property PeptideNTerminusMass() As Double
        Get
            Return mPeptideNTerminusMass
        End Get
        Set(Value As Double)
            mPeptideNTerminusMass = Value
        End Set
    End Property

    Public Property RemovePrefixAndSuffixIfPresent() As Boolean
        Get
            Return mRemovePrefixAndSuffixIfPresent
        End Get
        Set(Value As Boolean)
            mRemovePrefixAndSuffixIfPresent = Value
        End Set
    End Property
#End Region

    ''' <summary>
    ''' Constructor for shared (static) variables
    ''' </summary>
    Shared Sub New()

        mElementMonoMasses = GetElementMonoMasses()

        mAtomicFormulaRegEx = GetAtomicFormulaRegEx(mElementMonoMasses)

    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub New()
        ChargeCarrierMass = MASS_PROTON
        mErrorMessage = String.Empty
        mRemovePrefixAndSuffixIfPresent = True
        InitializeAminoAcidData()
    End Sub

    ''' <summary>
    ''' Compute the monoisotopic mass of the given empirical formula
    ''' </summary>
    ''' <param name="empiricalFormula"></param>
    ''' <returns></returns>
    ''' <remarks>Throws an exception if an unknown symbol is encountered</remarks>
    Public Shared Function ComputeMonoistopicMass(empiricalFormula As clsEmpiricalFormula) As Double

        Dim monoisotopicMass As Double = 0

        For Each element In empiricalFormula.ElementCounts

            Dim elementMass As Double
            If mElementMonoMasses.TryGetValue(element.Key, elementMass) Then
                monoisotopicMass += element.Value * elementMass
            Else
                Throw New Exception("Unrecognized symbol " & element.Key)
            End If

        Next

        Return monoisotopicMass

    End Function

    ''' <summary>
    ''' Compute the monoisotopic mass of the compound represented by elementalComposition
    ''' </summary>
    ''' <param name="elementalComposition"></param>
    ''' <param name="unknownSymbols"></param>
    ''' <returns></returns>
    Public Shared Function ComputeMonoistopicMass(elementalComposition As Dictionary(Of String, Integer), <Out()> ByRef unknownSymbols As List(Of String)) As Double

        Dim monoisotopicMass As Double = 0

        unknownSymbols = New List(Of String)

        For Each elementItem In elementalComposition

            Dim elementMass As Double
            If mElementMonoMasses.TryGetValue(elementItem.Key, elementMass) Then
                monoisotopicMass += elementItem.Value * elementMass
            Else
                unknownSymbols.Add(elementItem.Key)
            End If

        Next

        Return monoisotopicMass

    End Function

    ''' <summary>
    ''' Compute the mass of peptide sequence strSequence (it cannot contain modification symbols)
    ''' </summary>
    ''' <param name="strSequence">One letter amino acid symbols (no modification symbols or numbers); can have prefix and suffix letters</param>
    ''' <returns>Monoisotopic mass, or -1 if an error</returns>
    ''' <remarks>
    ''' Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True
    ''' If modification symbols are present, returns -1</remarks>
    Public Function ComputeSequenceMass(strSequence As String) As Double

        Dim strPrimarySequence As String = strSequence
        Dim dblMass As Double = 0
        Dim intValidResidueCount As Short = 0

        If mRemovePrefixAndSuffixIfPresent Then
            If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequence, strPrimarySequence, Nothing, Nothing) Then
                ' Prefix and suffix residues not present; simply copy strSequence to strPrimarySequence
                strPrimarySequence = strSequence
            End If
        End If

        If String.IsNullOrWhiteSpace(strPrimarySequence) Then
            ' This code should never be reached; including this as a fail-safe
            strPrimarySequence = strSequence
        End If

        mErrorMessage = String.Empty
        For Each chChar As Char In strPrimarySequence
            ' Use Convert.ToInt32 to convert to the Ascii value, then subtract 65
            Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(chChar)

            Try
                If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                    mErrorMessage = "Unknown symbol " & chChar & " in sequence " & strPrimarySequence
                    intValidResidueCount = 0
                    dblMass = -1
                    Exit For
                Else
                    dblMass += mAminoAcidMasses(aminoAcidIndex)
                    intValidResidueCount += 1S
                End If
            Catch ex As Exception
                ' Invalid value; ignore
            End Try
        Next chChar

        If intValidResidueCount > 0 Then
            dblMass += (mPeptideNTerminusMass + mPeptideCTerminusMass)
        End If

        Return dblMass

    End Function


    ''' <summary>
    ''' Compute the mass of peptide sequence strSequence; uses the information in udtResidueModificationInfo() to determine modification masses
    ''' </summary>
    ''' <param name="strSequence"></param>
    ''' <param name="intModCount"></param>
    ''' <param name="udtResidueModificationInfo">Array of modified residues; index 0 to intModCount-1</param>
    ''' <returns>The computed mass, or -1 if an error</returns>
    ''' <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
    <Obsolete("This version uses an array for modified residues; use the version that takes a list")>
    Public Function ComputeSequenceMass(strSequence As String, intModCount As Integer, ByRef udtResidueModificationInfo() As udtPeptideSequenceModInfoType) As Double

        Dim modifiedResidues = New List(Of udtPeptideSequenceModInfoType)

        If intModCount > 0 Then
            For intIndex = 0 To intModCount - 1
                modifiedResidues.Add(udtResidueModificationInfo(intIndex))
            Next
        End If

        Return ComputeSequenceMass(strSequence, modifiedResidues)

    End Function

    ''' <summary>
    ''' Compute the mass of peptide sequence strSequence; uses the information in udtResidueModificationInfo() to determine modification masses
    ''' </summary>
    ''' <param name="strSequence">One letter amino acid symbols (no modification symbols or numbers)</param>
    ''' <param name="modifiedResidues">List of modified residues</param>
    ''' <returns>The computed mass, or -1 if an error</returns>
    ''' <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
    Public Function ComputeSequenceMass(strSequence As String, modifiedResidues As List(Of udtPeptideSequenceModInfoType)) As Double

        ' Note that this call to ComputeSequenceMass will reset mErrorMessage
        Dim dblMass = ComputeSequenceMass(strSequence)

        If dblMass >= 0 AndAlso Not modifiedResidues Is Nothing AndAlso modifiedResidues.Count > 0 Then

            Dim empiricalFormula = New clsEmpiricalFormula()

            For Each modifiedResidue In modifiedResidues

                ' Note: do not use String.IsNullOrWhiteSpace(modifiedResidue.AffectedAtom) since that does not work on a char

                If modifiedResidue.AffectedAtom = Nothing OrElse modifiedResidue.AffectedAtom = NO_AFFECTED_ATOM_SYMBOL Then
                    ' Positional modification (static or dynamic mod)
                    ' Simply add the modification mass to dblMass
                    dblMass += modifiedResidue.ModificationMass
                    Continue For
                End If

                ' Isotopic modification
                If empiricalFormula.ElementCounts.Count = 0 Then
                    ' Initialize empiricalFormula using the amino acid sequence
                    Dim empiricalFormulaToAdd = ConvertAminoAcidSequenceToEmpiricalFormula(strSequence)
                    empiricalFormula.AddElements(empiricalFormulaToAdd)
                End If

                If Not mElementMonoMasses.ContainsKey(modifiedResidue.AffectedAtom) Then
                    mErrorMessage = "Unknown Affected Atom '" & modifiedResidue.AffectedAtom & "'"
                    dblMass = -1
                    Exit For
                End If

                Dim elementCount = empiricalFormula.GetElementCount(modifiedResidue.AffectedAtom)
                If elementCount = 0 Then
                    Console.WriteLine("Warning: no amino acids in {0} contain element {1}", strSequence,
                                        modifiedResidue.AffectedAtom)
                Else
                    dblMass += elementCount * modifiedResidue.ModificationMass
                End If

            Next
        End If

        Return dblMass

    End Function

    ''' <summary>
    ''' Compute the mass of peptide sequence strSequence.  Supports peptide sequences with with numeric mod masses
    ''' Examples of numeric mods:
    '''  R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A
    '''  K.Q-17.0265QIEESTSDYDKEK.L
    ''' </summary>
    ''' <param name="strSequence"></param>
    ''' <returns></returns>
    ''' <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
    Public Function ComputeSequenceMassNumericMods(strSequence As String) As Double

        Static reModMasses As Regex = New Regex("[+-][0-9.]+", RegexOptions.Compiled)

        Dim strPrimarySequence As String = String.Empty

        If mRemovePrefixAndSuffixIfPresent Then
            If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequence, strPrimarySequence, Nothing, Nothing) Then
                ' Prefix and suffix residues not present; simply copy strSequence to strPrimarySequence
                strPrimarySequence = String.Copy(strSequence)
            End If
        Else
            strPrimarySequence = String.Copy(strSequence)
        End If

        Dim reMatch = reModMasses.Match(strPrimarySequence)

        Dim sbSequenceWithoutMods = New StringBuilder
        Dim intStartIndex = 0
        Dim dblModMassTotal As Double = 0

        Do While reMatch.Success
            If reMatch.Index > intStartIndex Then
                sbSequenceWithoutMods.Append(strPrimarySequence.Substring(intStartIndex, reMatch.Index - intStartIndex))
            End If

            Dim strModMass = reMatch.ToString()
            Dim dblModMass As Double
            If Double.TryParse(strModMass, dblModMass) Then
                dblModMassTotal += dblModMass
            End If

            intStartIndex = reMatch.Index + strModMass.Length
            reMatch = reMatch.NextMatch()

        Loop

        If intStartIndex < strPrimarySequence.Length Then
            sbSequenceWithoutMods.Append(strPrimarySequence.Substring(intStartIndex, strPrimarySequence.Length - intStartIndex))
        End If

        Dim dblPeptideMass = ComputeSequenceMass(sbSequenceWithoutMods.ToString())

        If dblPeptideMass < 0 Then
            Return -1
        Else
            Return dblPeptideMass + dblModMassTotal
        End If

    End Function

    Private Function ConvertAminoAcidCharToIndex(aminoAcidSymbol As Char) As Short
        Return CShort(Convert.ToByte(aminoAcidSymbol)) - ASCII_VAL_LETTER_A
    End Function

    Private Function ConvertAminoAcidIndexToChar(aminoAcidIndex As Byte) As Char
        Return Convert.ToChar(aminoAcidIndex + ASCII_VAL_LETTER_A)
    End Function

    ''' <summary>
    ''' Converts the m/z value from one charge state to another charge state.  Either charge state can be 0, which means an uncharged peptide
    ''' </summary>
    ''' <param name="dblMassMZ"></param>
    ''' <param name="intCurrentCharge"></param>
    ''' <param name="intDesiredCharge"></param>
    ''' <returns></returns>
    ''' <remarks>Uses the charge carrier mass defined by ChargeCarrierMass</remarks>
    Public Function ConvoluteMass(
      dblMassMZ As Double,
      intCurrentCharge As Integer,
      Optional intDesiredCharge As Integer = 1
      ) As Double

        Return ConvoluteMass(dblMassMZ, intCurrentCharge, intDesiredCharge, ChargeCarrierMass)
    End Function

    ''' <summary>
    ''' Converts the m/z value from one charge state to another charge state.  Either charge state can be 0, which means an uncharged peptide
    ''' </summary>
    ''' <param name="dblMassMZ">m/z</param>
    ''' <param name="intCurrentCharge">Current charge; if 0, assumes dblMassMZ is the neutral, monoisotopic mass</param>
    ''' <param name="intDesiredCharge">Desired charge</param>
    ''' <param name="dblChargeCarrierMass">Charge carrier mass (Default is the mass of a proton)</param>
    ''' <returns></returns>
    ''' <remarks>To return the neutral mass, set intDesiredCharge to 0</remarks>
    Public Function ConvoluteMass(
      dblMassMZ As Double,
      intCurrentCharge As Integer,
      intDesiredCharge As Integer,
      dblChargeCarrierMass As Double) As Double

        Dim dblNewMZ As Double

        If Math.Abs(dblChargeCarrierMass) < Single.Epsilon Then
            dblChargeCarrierMass = MASS_PROTON
        End If

        Try
            If intCurrentCharge = intDesiredCharge Then
                dblNewMZ = dblMassMZ
            Else
                If intCurrentCharge = 1 Then
                    dblNewMZ = dblMassMZ
                ElseIf intCurrentCharge > 1 Then
                    ' Convert dblMassMZ to M+H
                    dblNewMZ = (dblMassMZ * intCurrentCharge) - dblChargeCarrierMass * (intCurrentCharge - 1)
                ElseIf intCurrentCharge = 0 Then
                    ' Convert dblMassMZ (which is neutral) to M+H and store in dblNewMZ
                    dblNewMZ = dblMassMZ + dblChargeCarrierMass
                Else
                    ' Negative charges are not supported; return 0
                    Return 0
                End If

                If intDesiredCharge > 1 Then
                    dblNewMZ = (dblNewMZ + dblChargeCarrierMass * (intDesiredCharge - 1)) / intDesiredCharge
                ElseIf intDesiredCharge = 1 Then
                    ' Return M+H, which is currently stored in dblNewMZ
                ElseIf intDesiredCharge = 0 Then
                    ' Return the neutral mass
                    dblNewMZ -= dblChargeCarrierMass
                Else
                    ' Negative charges are not supported; return 0
                    dblNewMZ = 0
                End If
            End If
        Catch ex As Exception
            ' Error occurred
            dblNewMZ = 0
        End Try

        Return dblNewMZ

    End Function

    ''' <summary>
    ''' Convert an amino acid sequence into an empirical formula
    ''' </summary>
    ''' <param name="strSequence">One letter amino acid symbols (no modification symbols or numbers)</param>
    ''' <returns></returns>
    Private Function ConvertAminoAcidSequenceToEmpiricalFormula(strSequence As String) As clsEmpiricalFormula

        Dim empiricalFormula = New clsEmpiricalFormula()

        For Each chAminoAcidSymbol As Char In strSequence
            ' Use Convert.ToInt32 to convert to the Ascii value, then subtract 65
            Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol)

            Try
                If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                    mErrorMessage = "Unknown symbol " & chAminoAcidSymbol & " in sequence " & strSequence
                    Exit For
                Else
                    empiricalFormula.AddElements(mAminoAcidEmpiricalFormulas(aminoAcidIndex))
                End If
            Catch ex As Exception
                ' Invalid value; ignore
            End Try
        Next

        Return empiricalFormula

    End Function

    ''' <summary>
    ''' Returns the mass of the specified amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function GetAminoAcidMass(chAminoAcidSymbol As Char) As Double
        ' Returns the mass if success, 0 if an error

        If Not chAminoAcidSymbol = Nothing Then
            Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol)
            If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return 0
            Else
                Return mAminoAcidMasses(aminoAcidIndex)
            End If
        End If

        Return 0

    End Function

    ''' <summary>
    ''' Returns a List with the number of atoms of C, H, N, O, and S in the specified amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function GetAminoAcidEmpiricalFormula(chAminoAcidSymbol As Char) As clsEmpiricalFormula
        ' Returns the atom counts if success, 0 if an error

        If chAminoAcidSymbol = Nothing Then
            Return New clsEmpiricalFormula()
        End If

        Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol)
        If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
            ' Invalid Index
            Return New clsEmpiricalFormula()
        Else
            Return mAminoAcidEmpiricalFormulas(aminoAcidIndex)
        End If

    End Function

    ''' <summary>
    ''' Create a new clsEmpiricalFormula instnance with the specified number of atoms
    ''' </summary>
    ''' <param name="countC"></param>
    ''' <param name="countH"></param>
    ''' <param name="countN"></param>
    ''' <param name="countO"></param>
    ''' <param name="countS"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function GetAminoAcidEmpiricalFormula(
      countC As Integer,
      countH As Integer,
      countN As Integer,
      countO As Integer,
      countS As Integer) As clsEmpiricalFormula

        Dim empiricalFormula = New clsEmpiricalFormula()

        empiricalFormula.AddElement("C", countC)
        empiricalFormula.AddElement("H", countH)
        empiricalFormula.AddElement("N", countN)
        empiricalFormula.AddElement("O", countO)
        empiricalFormula.AddElement("S", countS)

        Return empiricalFormula

    End Function

    ''' <summary>
    ''' Create a regex for parsing an empirical formula that optionally contains element counts and optionally contains plus or minus signs
    ''' Examples of supported empirical formulas:
    '''  CHNOS
    '''  C3H3NOS4
    '''  CH23NO-5S+4
    ''' </summary>
    ''' <param name="elementMonoMasses"></param>
    ''' <returns>RegEx with named capture groups ElementSymbol and ElementCount</returns>
    Private Shared Function GetAtomicFormulaRegEx(elementMonoMasses As Dictionary(Of String, Double)) As Regex

        Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline

        Dim sbRegEx = New StringBuilder()

        sbRegEx.Append("(?<ElementSymbol>")

        For Each element In elementMonoMasses
            sbRegEx.Append(element.Key & "|")
        Next

        ' Remove the trailing vertical bar
        sbRegEx.Remove(sbRegEx.Length - 1, 1)

        sbRegEx.Append(")")


        ' RegEx will be of the form: (?<ElementSymbol>H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al)(?<ElementCount>[+-]?\d*)
        Dim reAtomicFormulaRegEx = New Regex(sbRegEx.ToString() & "(?<ElementCount>[+-]?\d*)", REGEX_OPTIONS)

        Return reAtomicFormulaRegEx

    End Function

    Private Function GetDefaultAminoAcidMass(aminoAcidSymbol As Char, <Out> ByRef empiricalFormula As clsEmpiricalFormula) As Double

        ' These monoisotopic masses come from those traditionally used in DMS
        ' They were originally assembled by Gordon Anderson for use in ICR-2LS

        Dim monoMass As Double

        Select Case aminoAcidSymbol

            Case "A"c
                monoMass = 71.0371100902557
                empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 1, 0)
            Case "B"c
                ' Use N or D (aka Asn/Asp)
                monoMass = 114.042921543121
                empiricalFormula = GetAminoAcidEmpiricalFormula(4, 6, 2, 2, 0)
            Case "C"c
                monoMass = 103.009180784225
                empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 1, 1)
            Case "D"c
                monoMass = 115.026938199997
                empiricalFormula = GetAminoAcidEmpiricalFormula(4, 5, 1, 3, 0)
            Case "E"c
                monoMass = 129.042587518692
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 7, 1, 3, 0)
            Case "F"c
                monoMass = 147.068408727646
                empiricalFormula = GetAminoAcidEmpiricalFormula(9, 9, 1, 1, 0)
            Case "G"c
                monoMass = 57.0214607715607
                empiricalFormula = GetAminoAcidEmpiricalFormula(2, 3, 1, 1, 0)
            Case "H"c
                monoMass = 137.058904886246
                empiricalFormula = GetAminoAcidEmpiricalFormula(6, 7, 3, 1, 0)
            Case "I"c
                monoMass = 113.084058046341
                empiricalFormula = GetAminoAcidEmpiricalFormula(6, 11, 1, 1, 0)
            Case "J"c
                ' Could use mass of Ile/Leu, but we're instead treating this as an invalid, massless amino acid
                monoMass = 0
                empiricalFormula = GetAminoAcidEmpiricalFormula(0, 0, 0, 0, 0)
            Case "K"c
                monoMass = 128.094955444336
                empiricalFormula = GetAminoAcidEmpiricalFormula(6, 12, 2, 1, 0)
            Case "L"c
                monoMass = 113.084058046341
                empiricalFormula = GetAminoAcidEmpiricalFormula(6, 11, 1, 1, 0)
            Case "M"c
                monoMass = 131.040479421616
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 9, 1, 1, 1)
            Case "N"c
                monoMass = 114.042921543121
                empiricalFormula = GetAminoAcidEmpiricalFormula(4, 6, 2, 2, 0)
            Case "O"c
                monoMass = 114.079306125641
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 10, 2, 1, 0)
            Case "P"c
                monoMass = 97.0527594089508
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 7, 1, 1, 0)
            Case "Q"c
                monoMass = 128.058570861816
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 8, 2, 2, 0)
            Case "R"c
                monoMass = 156.101100921631
                empiricalFormula = GetAminoAcidEmpiricalFormula(6, 12, 4, 1, 0)
            Case "S"c
                monoMass = 87.0320241451263
                empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 2, 0)
            Case "T"c
                monoMass = 101.047673463821
                empiricalFormula = GetAminoAcidEmpiricalFormula(4, 7, 1, 2, 0)
            Case "U"c
                ' Corresponds to Sec = Selenocysteine (C3H5NOSe)
                monoMass = 150.95363
                empiricalFormula = GetAminoAcidEmpiricalFormula(3, 5, 1, 1, 0)
                empiricalFormula.AddElement("Se", 1)
            Case "V"c
                monoMass = 99.0684087276459
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 9, 1, 1, 0)
            Case "W"c
                monoMass = 186.079306125641
                empiricalFormula = GetAminoAcidEmpiricalFormula(11, 10, 2, 1, 0)
            Case "X"c
                ' Unknown; use mass of Ile/Leu
                monoMass = 113.084058046341
                empiricalFormula = GetAminoAcidEmpiricalFormula(6, 11, 1, 1, 0)
            Case "Y"c
                monoMass = 163.063322782516
                empiricalFormula = GetAminoAcidEmpiricalFormula(9, 9, 1, 2, 0)
            Case "Z"c
                ' use Q or E (aka Gln/Glu); note that these are 0.984 Da apart
                monoMass = 128.058570861816
                empiricalFormula = GetAminoAcidEmpiricalFormula(5, 8, 2, 2, 0)
            Case Else
                monoMass = 0
                empiricalFormula = GetAminoAcidEmpiricalFormula(0, 0, 0, 0, 0)
        End Select

        Dim computedMass = ComputeMonoistopicMass(empiricalFormula)
        If Math.Abs(computedMass - monoMass) > 0.01 Then
            Console.WriteLine("Mass discrepancy for amino acid {0}. DMS uses {1:F4} but this class computed {2:F4}", aminoAcidSymbol, monoMass, computedMass)
        End If

        Return monoMass

    End Function

    ''' <summary>
    ''' Return a dictionary of element symbols and element masses
    ''' </summary>
    ''' <returns></returns>
    Private Shared Function GetElementMonoMasses() As Dictionary(Of String, Double)
        Dim elementMonoMasses = New Dictionary(Of String, Double) From {
            {"H", MASS_HYDROGEN}, {"He", 4.0026029}, {"Li", 7.016005}, {"Be", 9.012183},
            {"B", 11.009305}, {"C", 12}, {"N", 14.003074}, {"O", MASS_OXYGEN},
            {"F", 18.9984032}, {"Ne", 19.992439}, {"Na", 22.98977}, {"Mg", 23.98505},
            {"Al", 26.981541}, {"Si", 27.976928}, {"P", 30.973763}, {"S", 31.972072},
            {"Cl", 34.968853}, {"Ar", 39.962383}, {"K", 38.963708}, {"Ca", 39.962591},
            {"Sc", 44.955914}, {"Ti", 47.947947}, {"V", 50.943963}, {"Cr", 51.94051},
            {"Mn", 54.938046}, {"Fe", 55.934939}, {"Co", 58.933198}, {"Ni", 57.935347},
            {"Cu", 62.929599}, {"Zn", 63.929145}, {"Ga", 68.925581}, {"Ge", 71.92208},
            {"As", 74.921596}, {"Se", 79.916521}, {"Br", 78.918336}, {"Kr", 83.911506},
            {"Rb", 84.9118}, {"Sr", 87.905625}, {"Y", 88.905856}, {"Zr", 89.904708},
            {"Nb", 92.906378}, {"Mo", 97.905405}, {"Tc", 98}, {"Ru", 101.90434},
            {"Rh", 102.905503}, {"Pd", 105.903475}, {"Ag", 106.905095}, {"Cd", 113.903361},
            {"In", 114.903875}, {"Sn", 119.902199}, {"Sb", 120.903824}, {"Te", 129.906229},
            {"I", 126.904477}, {"Xe", 131.904148}, {"Cs", 132.905433}, {"Ba", 137.905236},
            {"La", 138.906355}, {"Ce", 139.905442}, {"Pr", 140.907657}, {"Nd", 141.907731},
            {"Pm", 145}, {"Sm", 151.919741}, {"Eu", 152.921243}, {"Gd", 157.924111},
            {"Tb", 158.92535}, {"Dy", 163.929183}, {"Ho", 164.930332}, {"Er", 165.930305},
            {"Tm", 168.934225}, {"Yb", 173.938873}, {"Lu", 174.940785}, {"Hf", 179.946561},
            {"Ta", 180.948014}, {"W", 183.950953}, {"Re", 186.955765}, {"Os", 191.960603},
            {"Ir", 192.962942}, {"Pt", 194.964785}, {"Au", 196.96656}, {"Hg", 201.970632},
            {"Tl", 204.97441}, {"Pb", 207.976641}, {"Bi", 208.980388}, {"Po", 209},
            {"At", 210}, {"Rn", 222}, {"Fr", 223}, {"Ra", 227},
            {"Ac", 227}, {"Th", 232.038054}, {"Pa", 231}, {"U", 238.050786},
            {"Np", 237}, {"Pu", 244}, {"Am", 243}, {"Cm", 247},
            {"Bk", 247}, {"Cf", 251}, {"Es", 252}, {"Fm", 257},
            {"Md", 258}, {"No", 269}, {"Lr", 260}
        }

        Return elementMonoMasses
    End Function

    ''' <summary>
    ''' Parse the given empirical formula to return a dictionary of the elements
    ''' Examples of supported empirical formulas:
    '''  CHNOS
    '''  C3H3NOS4
    '''  CH23NO-5S+4
    ''' </summary>
    ''' <param name="strEmpiricalformula"></param>
    ''' <returns>EmpiricalFormula instnance tracking the element symbols and counts</returns>
    Public Shared Function GetEmpiricalFormulaComponents(strEmpiricalformula As String) As clsEmpiricalFormula

        ' Originally MSGF+ only allowed for elements C, H, N, O, S, and P in a dynamic or static mod definition
        ' It now allows for any element

        Dim reMatches As MatchCollection = mAtomicFormulaRegEx.Matches(strEmpiricalformula)

        Dim empiricalFormula = New clsEmpiricalFormula()

        If reMatches.Count > 0 Then

            For Each reMatch As Match In reMatches

                Dim elementSymbol = reMatch.Groups("ElementSymbol").ToString()
                Dim elementCountText = reMatch.Groups("ElementCount").ToString()

                Dim elementCount = 1
                If Not String.IsNullOrEmpty(elementCountText) AndAlso elementCountText.Length > 0 Then
                    If Not Integer.TryParse(elementCountText, elementCount) Then
                        Throw New Exception("Error parsing empirical formula '" & strEmpiricalformula & "', number not found in " & elementCountText)
                    End If
                End If

                empiricalFormula.AddElement(elementSymbol, elementCount)

            Next
        End If

        Return empiricalFormula

    End Function

    Private Sub InitializeAminoAcidData()

        ReDim mAminoAcidMasses(AMINO_ACID_LIST_MAX_INDEX)
        ReDim mAminoAcidEmpiricalFormulas(AMINO_ACID_LIST_MAX_INDEX)

        For index As Byte = 0 To AMINO_ACID_LIST_MAX_INDEX
            UpdateAminoAcidStatEntry(index)
        Next

        ResetTerminusMasses()
    End Sub

    ''' <summary>
    ''' Converts dblMassToConvert to ppm, based on the value of dblCurrentMZ
    ''' </summary>
    ''' <param name="dblMassToConvert"></param>
    ''' <param name="dblCurrentMZ"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function MassToPPM(dblMassToConvert As Double, dblCurrentMZ As Double) As Double
        Return dblMassToConvert * 1000000.0 / dblCurrentMZ
    End Function

    ''' <summary>
    ''' Converts and MH mass to the uncharged (neutral) mass
    ''' </summary>
    ''' <param name="dblMH"></param>
    ''' <returns></returns>
    ''' <remarks>Equivalent to ConvoluteMass(dblMH, 1, 0)</remarks>
    Public Function MHToMonoisotopicMass(dblMH As Double) As Double
        Return ConvoluteMass(dblMH, 1, 0)
    End Function

    ''' <summary>
    ''' Converts an uncharged (neutral) mass to the m/z value for the specified charge
    ''' </summary>
    ''' <param name="dblMonoisotopicMass"></param>
    ''' <param name="intDesiredCharge"></param>
    ''' <returns></returns>
    ''' <remarks>Equivalent to ConvoluteMass(dblMonoisotopicMass, 0, intDesiredCharge)</remarks>
    Public Function MonoisotopicMassToMZ(dblMonoisotopicMass As Double, intDesiredCharge As Integer) As Double
        Return ConvoluteMass(dblMonoisotopicMass, 0, intDesiredCharge)
    End Function

    ''' <summary>
    ''' Converts from a ppm value to a mass value, using the specified m/z as a reference point
    ''' </summary>
    ''' <param name="dblPPMToConvert"></param>
    ''' <param name="dblCurrentMZ"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function PPMToMass(dblPPMToConvert As Double, dblCurrentMZ As Double) As Double
        ' Converts dblPPMToConvert to a mass value, which is dependent on dblCurrentMZ

        Return dblPPMToConvert / 1000000.0 * dblCurrentMZ
    End Function

    ''' <summary>
    ''' Reset all of the amino acid masses and atom counts to default values
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ResetAminoAcidMasses()

        For aminoAcidIndex As Byte = 0 To AMINO_ACID_LIST_MAX_INDEX
            Dim aminoAcidSymbol As Char = ConvertAminoAcidIndexToChar(aminoAcidIndex)
            ResetAminoAcidToDefault(aminoAcidSymbol)
        Next
    End Sub

    ''' <summary>
    ''' Reset the mass and atom counts of the given amino acid to use default values
    ''' </summary>
    ''' <param name="aminoAcidSymbol">Letter between A and Z</param>
    ''' <remarks></remarks>
    Public Sub ResetAminoAcidToDefault(aminoAcidSymbol As Char)

        Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(aminoAcidSymbol)
        If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
            ' Invalid Index
            Return
        End If

        UpdateAminoAcidStatEntry(CByte(aminoAcidIndex))

    End Sub

    ''' <summary>
    ''' Reset the N and C terminus default mass values
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ResetTerminusMasses()
        ' See comment in Sub InitializeAminoAcidData concerning these masses

        mPeptideNTerminusMass = DEFAULT_N_TERMINUS_MASS_CHANGE
        mPeptideCTerminusMass = DEFAULT_C_TERMINUS_MASS_CHANGE
    End Sub

    ''' <summary>
    ''' Defines the number of C, H, N, O, S, etc. elements in an amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol">Amino acid symbol</param>
    ''' <param name="elementalComposition">Dictionary where keys are element symbols and values are the element counts</param>
    ''' <returns>True if success, False if an invalid amino acid symbol</returns>
    ''' <remarks></remarks>
    Public Function SetAminoAcidAtomCounts(chAminoAcidSymbol As Char, elementalComposition As Dictionary(Of String, Integer)) As Boolean
        Dim empiricalFormula = New clsEmpiricalFormula(elementalComposition)
        Dim success = SetAminoAcidAtomCounts(chAminoAcidSymbol, empiricalFormula)
        Return success
    End Function

    ''' <summary>
    ''' Defines the number of C, H, N, O, S, etc. elements in an amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol">>Amino acid symbol</param>
    ''' <param name="empiricalFormula">Empirical formula class</param>
    ''' <returns>True if success, False if an invalid amino acid symbol</returns>
    ''' <remarks></remarks>
    Public Function SetAminoAcidAtomCounts(chAminoAcidSymbol As Char, empiricalFormula As clsEmpiricalFormula) As Boolean

        If chAminoAcidSymbol = Nothing Then
            Return False
        End If

        Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol)
        If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
            ' Invalid Index
            Return False
        Else
            mAminoAcidEmpiricalFormulas(aminoAcidIndex) = empiricalFormula
            Return True
        End If

    End Function

    ''' <summary>
    ''' Defines a custom mass for an amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <param name="dblMass"></param>
    ''' <returns>True if success, False if an invalid amino acid symbol</returns>
    ''' <remarks></remarks>
    Public Function SetAminoAcidMass(chAminoAcidSymbol As Char, dblMass As Double) As Boolean

        If Not chAminoAcidSymbol = Nothing Then
            Dim aminoAcidIndex = ConvertAminoAcidCharToIndex(chAminoAcidSymbol)
            If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return False
            Else
                mAminoAcidMasses(aminoAcidIndex) = dblMass
                Return True
            End If
        End If

        Return False

    End Function

    ''' <summary>
    ''' Updates an entry in parallel arrays AminoAcidMasses and AminoAcidSymbols
    ''' </summary>
    ''' <param name="aminoAcidIndex"></param>
    ''' <remarks></remarks>
    Private Sub UpdateAminoAcidStatEntry(aminoAcidIndex As Byte)

        ' Use Convert.ToChar to convert from Ascii code to the letter
        Dim aminoAcidSymbol As Char = ConvertAminoAcidIndexToChar(aminoAcidIndex)

        Dim empiricalFormula As clsEmpiricalFormula = Nothing
        Dim monoMass = GetDefaultAminoAcidMass(aminoAcidSymbol, empiricalFormula)

        mAminoAcidMasses(aminoAcidIndex) = monoMass
        mAminoAcidEmpiricalFormulas(aminoAcidIndex) = empiricalFormula

    End Sub

End Class

