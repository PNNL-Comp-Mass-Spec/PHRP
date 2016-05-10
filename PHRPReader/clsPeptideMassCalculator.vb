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
        Public ResidueLocInPeptide As Integer           ' Position that the modification occurs; not used by this class
        Public ModificationMass As Double
        Public AffectedAtom As Char                     ' Set to Nothing or to NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications); for Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
    End Structure

    ''' <summary>
    ''' Tracks elemental composition of each amino acid
    ''' </summary>
    ''' <remarks>In MGSF+, amino acids can only be C, H, N, O, and S</remarks>
    Public Structure udtAtomCountsType
        Public CountC As Integer
        Public CountH As Integer
        Public CountN As Integer
        Public CountO As Integer
        Public CountS As Integer

        Public Overrides Function ToString() As String
            Return "C" & CountC & "H" & CountH & "N" & CountN & "O" & CountO & "S" & CountS
        End Function
    End Structure

#End Region

#Region "Classwide Variables"
    ' The Amino Acid arrays contain 26 entries, corresponding to A through Z
    ' Invalid/Undefined letters (J and U) have values of 0 for the mass and atom counts
    ' The values can be customized using SetAminoAcidMass and SetAminoAcidAtomCounts

    Private Const AMINO_ACID_LIST_MAX_INDEX As Integer = 25
    Private mAminoAcidMasses() As Double
    Private mAminoAcidAtomCounts() As udtAtomCountsType

    ' typically mPeptideNTerminusMass + mPeptideCTerminusMass = 18.0105633 (the mass of water)
    Private mPeptideNTerminusMass As Double
    Private mPeptideCTerminusMass As Double

    Private mRemovePrefixAndSuffixIfPresent As Boolean

    Private mErrorMessage As String
#End Region

#Region "Properties"

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
    ''' Constructor
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub New()
        mErrorMessage = String.Empty
        mRemovePrefixAndSuffixIfPresent = True
        InitializeAminoAcidData()
    End Sub

    ''' <summary>
    ''' Adds a new entry to parallel arrays AminoAcidMasses and AminoAcidSymbols
    ''' </summary>
    ''' <param name="intIndex"></param>
    ''' <remarks></remarks>
    Private Sub AddAminoAcidStatEntry(intIndex As Integer)

        ' Use Convert.ToChar to convert from Ascii code to the letter
        Dim aminoAcidSymbol As Char = Convert.ToChar(intIndex + ASCII_VAL_LETTER_A)

        Dim udtAtomCounts As udtAtomCountsType = Nothing
        Dim monoMass = GetDefaultAminoAcidMass(aminoAcidSymbol, udtAtomCounts)

        mAminoAcidMasses(intIndex) = monoMass
        mAminoAcidAtomCounts(intIndex) = udtAtomCounts

    End Sub

    ''' <summary>
    ''' Compute the mass of peptide sequence strSequence.  If modification symbols are present, returns -1
    ''' </summary>
    ''' <param name="strSequence"></param>
    ''' <returns></returns>
    ''' <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
    Public Function ComputeSequenceMass(strSequence As String) As Double
        ' Computes the mass for sequence strSequence
        ' Returns -1 if an error

        Dim strPrimarySequence As String = String.Empty
        Dim dblMass As Double = 0
        Dim intValidResidueCount As Short = 0

        If mRemovePrefixAndSuffixIfPresent Then
            If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequence, strPrimarySequence, Nothing, Nothing) Then
                ' Prefix and suffix residues not present; simply copy strSequence to strPrimarySequence
                strPrimarySequence = String.Copy(strSequence)
            End If
        End If

        mErrorMessage = String.Empty
        For Each chChar As Char In strSequence
            ' Use Convert.ToInt32 to convert to the Ascii value, then subtract 65
            Dim aminoAcidIndex = Convert.ToInt32(chChar) - ASCII_VAL_LETTER_A

            Try
                If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                    mErrorMessage = "Unknown symbol " & chChar & " in sequence " & strSequence
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
    ''' <returns></returns>
    ''' <remarks>Looks for and removes prefix and suffix letters if .RemovePrefixAndSuffixIfPresent = True</remarks>
    Public Function ComputeSequenceMass(strSequence As String, intModCount As Integer, ByRef udtResidueModificationInfo() As udtPeptideSequenceModInfoType) As Double
        ' Computes the mass for sequence strSequence using the mods in udtResidueModificationInfo()
        ' Returns -1 if an error

        ' Note that ComputeSequenceMass will reset mErorMessage
        Dim dblMass = ComputeSequenceMass(strSequence)

        If intModCount > 0 AndAlso dblMass >= 0 Then

            Dim blnAtomCountsDefined = False
            Dim udtAtomCounts As udtAtomCountsType

            For intIndex = 0 To intModCount - 1
                With udtResidueModificationInfo(intIndex)

                    If .AffectedAtom = Nothing OrElse .AffectedAtom = NO_AFFECTED_ATOM_SYMBOL Then
                        ' Positional modification (static or dynamic mod)
                        ' Simply add the modification mass to dblMass
                        dblMass += .ModificationMass
                    Else
                        
                        ' Isotopic modification
                        If Not blnAtomCountsDefined Then
                            udtAtomCounts = CountAtoms(strSequence)
                            blnAtomCountsDefined = True
                        End If

                        Select Case Char.ToUpper(.AffectedAtom)
                            Case "C"c
                                dblMass += (udtAtomCounts.CountC * .ModificationMass)
                            Case "H"c
                                dblMass += (udtAtomCounts.CountH * .ModificationMass)
                            Case "N"c
                                dblMass += (udtAtomCounts.CountN * .ModificationMass)
                            Case "O"c
                                dblMass += (udtAtomCounts.CountO * .ModificationMass)
                            Case "S"c
                                dblMass += (udtAtomCounts.CountS * .ModificationMass)
                            Case Else
                                mErrorMessage = "Unknown Affected Atom '" & .AffectedAtom & "'"
                                dblMass = -1
                                Exit For
                        End Select
                    End If

                End With
            Next intIndex
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

    ''' <summary>
    ''' Populate a udtAtomCountsType instance using the information in elementalComposition
    ''' </summary>
    ''' <param name="elementalComposition"></param>
    ''' <returns></returns>
    ''' <remarks>Only valid for amino acids, since udtAtomCountsType only supports C, H, N, O, and S</remarks>
    Public Shared Function ConvertElementalCompositionToAtomCounts(elementalComposition As Dictionary(Of Char, Integer)) As udtAtomCountsType

        Dim atomCounts = New udtAtomCountsType()

        For Each elementItem In elementalComposition

            Dim elementSymbol As Char = elementItem.Key

            Select Case Char.ToUpper(elementSymbol)
                Case "C"c : atomCounts.CountC = elementItem.Value
                Case "H"c : atomCounts.CountH = elementItem.Value
                Case "N"c : atomCounts.CountN = elementItem.Value
                Case "O"c : atomCounts.CountO = elementItem.Value
                Case "S"c : atomCounts.CountS = elementItem.Value
                Case Else
                    ' Unknown element
                    Throw New Exception("Error parsing items in elementalComposition, unknown element " & elementItem.Key & "; must be C, H, N, O, or S")
            End Select
        Next

        Return atomCounts

    End Function

    ''' <summary>
    ''' Converts the m/z value from one charge state to another charge state.  Either charge state can be 0, which means an uncharged peptide
    ''' </summary>
    ''' <param name="dblMassMZ"></param>
    ''' <param name="intCurrentCharge"></param>
    ''' <param name="intDesiredCharge"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function ConvoluteMass(dblMassMZ As Double, intCurrentCharge As Integer, Optional intDesiredCharge As Integer = 1) As Double
        ' Converts dblMassMZ to the MZ that would appear at the given intDesiredCharge
        ' If intCurrentCharge = 0, then assumes dblMassMZ is the neutral, monoisotopic mass
        ' To return the neutral mass, set intDesiredCharge to 0

        Dim dblNewMZ As Double

        Try
            If intCurrentCharge = intDesiredCharge Then
                dblNewMZ = dblMassMZ
            Else
                If intCurrentCharge = 1 Then
                    dblNewMZ = dblMassMZ
                ElseIf intCurrentCharge > 1 Then
                    ' Convert dblMassMZ to M+H
                    dblNewMZ = (dblMassMZ * intCurrentCharge) - MASS_PROTON * (intCurrentCharge - 1)
                ElseIf intCurrentCharge = 0 Then
                    ' Convert dblMassMZ (which is neutral) to M+H and store in dblNewMZ
                    dblNewMZ = dblMassMZ + MASS_PROTON
                Else
                    ' Negative charges are not supported; return 0
                    Return 0
                End If

                If intDesiredCharge > 1 Then
                    dblNewMZ = (dblNewMZ + MASS_PROTON * (intDesiredCharge - 1)) / intDesiredCharge
                ElseIf intDesiredCharge = 1 Then
                    ' Return M+H, which is currently stored in dblNewMZ
                ElseIf intDesiredCharge = 0 Then
                    ' Return the neutral mass
                    dblNewMZ -= MASS_PROTON
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

    Private Function CountAtoms(ByRef strSequence As String) As udtAtomCountsType

        Dim atomCounts = New udtAtomCountsType()

        For Each chChar As Char In strSequence
            ' Use Convert.ToInt32 to convert to the Ascii value, then subtract 65
            Dim aminoAcidIndex = Convert.ToInt32(chChar) - ASCII_VAL_LETTER_A

            Try
                If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                    mErrorMessage = "Unknown symbol " & chChar & " in sequence " & strSequence
                    Exit For
                Else
                    With atomCounts
                        .CountC += mAminoAcidAtomCounts(aminoAcidIndex).CountC
                        .CountH += mAminoAcidAtomCounts(aminoAcidIndex).CountH
                        .CountN += mAminoAcidAtomCounts(aminoAcidIndex).CountN
                        .CountO += mAminoAcidAtomCounts(aminoAcidIndex).CountO
                        .CountS += mAminoAcidAtomCounts(aminoAcidIndex).CountS
                    End With
                End If
            Catch ex As Exception
                ' Invalid value; ignore
            End Try
        Next chChar

        Return atomCounts

    End Function

    ''' <summary>
    ''' Returns a structure with the number of atoms of C, H, N, O, and S in the specified amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function GetAminoAcidAtomCounts(chAminoAcidSymbol As Char) As udtAtomCountsType
        ' Returns the atom counts if success, 0 if an error

        Dim aminoAcidIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            aminoAcidIndex = Convert.ToInt32(chAminoAcidSymbol) - ASCII_VAL_LETTER_A
            If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return New udtAtomCountsType
            Else
                Return mAminoAcidAtomCounts(aminoAcidIndex)
            End If
        End If

        Return New udtAtomCountsType

    End Function

    ''' <summary>
    ''' Returns the mass of the specified amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function GetAminoAcidMass(chAminoAcidSymbol As Char) As Double
        ' Returns the mass if success, 0 if an error

        Dim aminoAcidIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            aminoAcidIndex = Convert.ToInt32(chAminoAcidSymbol) - ASCII_VAL_LETTER_A
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
    ''' Create a new udtAtomCountsType struct with the specified number of atoms
    ''' </summary>
    ''' <param name="countC"></param>
    ''' <param name="countH"></param>
    ''' <param name="countN"></param>
    ''' <param name="countO"></param>
    ''' <param name="countS"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function GetAtomCountsStruct(
      countC As Integer,
      countH As Integer,
      countN As Integer,
      countO As Integer,
      countS As Integer) As udtAtomCountsType

        Dim atomCounts = New udtAtomCountsType()

        atomCounts.CountC = countC
        atomCounts.CountH = countH
        atomCounts.CountN = countN
        atomCounts.CountO = countO
        atomCounts.CountS = countS

        Return atomCounts

    End Function

    Private Function GetDefaultAminoAcidMass(aminoAcidSymbol As Char, <Out> ByRef udtAtomCounts As udtAtomCountsType) As Double

        ' These monoisotopic masses come from those traditionally used in DMS
        ' They were originally assembled by Gordon Anderson for use in ICR-2LS
        
        Dim monoMass As Double

        Select Case aminoAcidSymbol

            Case "A"c
                monoMass = 71.0371100902557 : udtAtomCounts = GetAtomCountsStruct(3, 5, 1, 1, 0)
            Case "B"c
                ' Use N or D (aka Asn/Asp)
                monoMass = 114.042921543121 : udtAtomCounts = GetAtomCountsStruct(4, 6, 2, 2, 0)
            Case "C"c
                monoMass = 103.009180784225 : udtAtomCounts = GetAtomCountsStruct(3, 5, 1, 1, 1)
            Case "D"c
                monoMass = 115.026938199997 : udtAtomCounts = GetAtomCountsStruct(4, 5, 1, 3, 0)
            Case "E"c
                monoMass = 129.042587518692 : udtAtomCounts = GetAtomCountsStruct(5, 7, 1, 3, 0)
            Case "F"c
                monoMass = 147.068408727646 : udtAtomCounts = GetAtomCountsStruct(9, 9, 1, 1, 0)
            Case "G"c
                monoMass = 57.0214607715607 : udtAtomCounts = GetAtomCountsStruct(2, 3, 1, 1, 0)
            Case "H"c
                monoMass = 137.058904886246 : udtAtomCounts = GetAtomCountsStruct(6, 7, 3, 1, 0)
            Case "I"c
                monoMass = 113.084058046341 : udtAtomCounts = GetAtomCountsStruct(6, 11, 1, 1, 0)
            Case "J"c
                ' Could use mass of Ile/Leu, but we're using ""
                monoMass = 0 : udtAtomCounts = GetAtomCountsStruct(0, 0, 0, 0, 0)
            Case "K"c
                monoMass = 128.094955444336 : udtAtomCounts = GetAtomCountsStruct(6, 12, 2, 1, 0)
            Case "L"c
                monoMass = 113.084058046341 : udtAtomCounts = GetAtomCountsStruct(6, 11, 1, 1, 0)
            Case "M"c
                monoMass = 131.040479421616 : udtAtomCounts = GetAtomCountsStruct(5, 9, 1, 1, 1)
            Case "N"c
                monoMass = 114.042921543121 : udtAtomCounts = GetAtomCountsStruct(4, 6, 2, 2, 0)
            Case "O"c
                monoMass = 114.079306125641 : udtAtomCounts = GetAtomCountsStruct(5, 10, 2, 1, 0)
            Case "P"c
                monoMass = 97.0527594089508 : udtAtomCounts = GetAtomCountsStruct(5, 7, 1, 1, 0)
            Case "Q"c
                monoMass = 128.058570861816 : udtAtomCounts = GetAtomCountsStruct(5, 8, 2, 2, 0)
            Case "R"c
                monoMass = 156.101100921631 : udtAtomCounts = GetAtomCountsStruct(6, 12, 4, 1, 0)
            Case "S"c
                monoMass = 87.0320241451263 : udtAtomCounts = GetAtomCountsStruct(3, 5, 1, 2, 0)
            Case "T"c
                monoMass = 101.047673463821 : udtAtomCounts = GetAtomCountsStruct(4, 7, 1, 2, 0)
            Case "U"c
                ' Corresponds to Sec = Selenocysteine (C3H5NOSe)
                monoMass = 150.95363 : udtAtomCounts = GetAtomCountsStruct(3, 5, 1, 1, 0)
            Case "V"c
                monoMass = 99.0684087276459 : udtAtomCounts = GetAtomCountsStruct(5, 9, 1, 1, 0)
            Case "W"c
                monoMass = 186.079306125641 : udtAtomCounts = GetAtomCountsStruct(11, 10, 2, 1, 0)
            Case "X"c
                ' Unknown; use mass of Ile/Leu
                monoMass = 113.084058046341 : udtAtomCounts = GetAtomCountsStruct(6, 11, 1, 1, 0)
            Case "Y"c
                monoMass = 163.063322782516 : udtAtomCounts = GetAtomCountsStruct(9, 9, 1, 2, 0)
            Case "Z"c
                ' use Q or E (aka Gln/Glu); note that these are 0.984 Da apart
                monoMass = 128.058570861816 : udtAtomCounts = GetAtomCountsStruct(5, 8, 2, 2, 0)
            Case Else
                monoMass = 0 : udtAtomCounts = GetAtomCountsStruct(0, 0, 0, 0, 0)
        End Select

        Return monoMass

    End Function
    
    Public Shared Function GetEmpiricalFormulaComponents(strEmpiricalformula As String) As Dictionary(Of Char, Integer)

        Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline Or RegexOptions.IgnoreCase

        Static reAtomicFormulaRegEx As New Regex("[CHNOSP][+-]?\d*", REGEX_OPTIONS)
        Dim reMatches As MatchCollection = reAtomicFormulaRegEx.Matches(strEmpiricalformula)

        Dim elementalComposition = New Dictionary(Of Char, Integer)

        If reMatches.Count > 0 Then

            For Each reMatch As Match In reMatches
                Dim strElement = reMatch.Value.Chars(0)
                Dim intCount As Integer

                If reMatch.Value.Length > 1 Then

                    If Not Integer.TryParse(reMatch.Value.Substring(1), intCount) Then
                        Throw New Exception("Error parsing empirical formula '" & strEmpiricalformula & "', number not found in " & reMatch.Value)
                    End If
                Else
                    intCount = 1
                End If

                elementalComposition.Add(strElement, intCount)
            Next
        End If

        Return elementalComposition

    End Function

    Private Sub InitializeAminoAcidData()

        ReDim mAminoAcidMasses(AMINO_ACID_LIST_MAX_INDEX)
        ReDim mAminoAcidAtomCounts(AMINO_ACID_LIST_MAX_INDEX)

        For index = 0 To AMINO_ACID_LIST_MAX_INDEX
            AddAminoAcidStatEntry(index)
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
    Public Shared Function MHToMonoisotopicMass(dblMH As Double) As Double
        Return ConvoluteMass(dblMH, 1, 0)
    End Function

    ''' <summary>
    ''' Converts an uncharged (neutral) mass to the m/z value for the specified charge
    ''' </summary>
    ''' <param name="dblMonoisotopicMass"></param>
    ''' <param name="intDesiredCharge"></param>
    ''' <returns></returns>
    ''' <remarks>Equivalent to ConvoluteMass(dblMonoisotopicMass, 0, intDesiredCharge)</remarks>
    Public Shared Function MonoisotopicMassToMZ(dblMonoisotopicMass As Double, intDesiredCharge As Integer) As Double
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
    ''' Reset the mass and atom counts of the given amino acid to use default values
    ''' </summary>
    ''' <param name="aminoAcidSymbol">Letter between A and Z</param>
    ''' <remarks></remarks>
    Public Sub ResetAminoAcidToDefault(aminoAcidSymbol As Char)

        Dim aminoAcidIndex = Convert.ToInt32(aminoAcidSymbol) - ASCII_VAL_LETTER_A
        If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
            ' Invalid Index
            Return
        End If

        Dim udtAtomCounts As udtAtomCountsType = Nothing
        Dim monoMass = GetDefaultAminoAcidMass(aminoAcidSymbol, udtAtomCounts)

        mAminoAcidMasses(aminoAcidIndex) = monoMass
        mAminoAcidAtomCounts(aminoAcidIndex) = udtAtomCounts

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
    ''' Defines the number of C, H, N, O, and S atoms in an amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <param name="udtAtomCounts"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function SetAminoAcidAtomCounts(chAminoAcidSymbol As Char, udtAtomCounts As udtAtomCountsType) As Boolean
        ' Returns True if success, False if an invalid amino acid symbol

        Dim aminoAcidIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            aminoAcidIndex = Convert.ToInt32(chAminoAcidSymbol) - ASCII_VAL_LETTER_A
            If aminoAcidIndex < 0 OrElse aminoAcidIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return False
            Else
                mAminoAcidAtomCounts(aminoAcidIndex) = udtAtomCounts
                Return True
            End If
        End If

        Return False

    End Function

    ''' <summary>
    ''' Defines a custom mass for an amino acid
    ''' </summary>
    ''' <param name="chAminoAcidSymbol"></param>
    ''' <param name="dblMass"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function SetAminoAcidMass(chAminoAcidSymbol As Char, dblMass As Double) As Boolean
        ' Returns True if success, False if an invalid amino acid symbol

        Dim aminoAcidIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            aminoAcidIndex = Convert.ToInt32(chAminoAcidSymbol) - ASCII_VAL_LETTER_A
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
End Class

