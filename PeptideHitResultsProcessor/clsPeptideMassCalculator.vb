Option Strict On

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
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'
' Notice: This computer software was prepared by Battelle Memorial Institute, 
' hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
' Department of Energy (DOE).  All rights in the computer software are reserved 
' by DOE on behalf of the United States Government and the Contractor as 
' provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
' SOFTWARE.  This notice including this sentence must appear on any copies of 
' this computer software.

Public Class clsPeptideMassCalculator

#Region "Constants and Enums"
    Public Const NO_AFFECTED_ATOM_SYMBOL As Char = "-"c

    Public Const MASS_HYDROGEN As Double = 1.0078246
    Public Const MASS_OXYGEN As Double = 15.9949141
    Public Const MASS_PROTON As Double = 1.00727649               ' Note that this is the mass of hydrogen minus the mass of one electron
    Public Const MASS_ELECTRON As Double = 0.00054811

    Public Const DEFAULT_N_TERMINUS_MASS_CHANGE As Double = MASS_HYDROGEN
    Public Const DEFAULT_C_TERMINUS_MASS_CHANGE As Double = MASS_OXYGEN + MASS_HYDROGEN

#End Region

#Region "Structures"
    Public Structure udtPeptideSequenceModInfoType
        Public ResidueLocInPeptide As Integer           ' Position that the modification occurs; not used by this class
        Public ModificationMass As Double
        Public AffectedAtom As Char                     ' Set to Nothing or to NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications); for Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
    End Structure

    Public Structure udtAtomCountsType
        Public CountC As Integer
        Public CountH As Integer
        Public CountN As Integer
        Public CountO As Integer
        Public CountS As Integer
    End Structure

#End Region

#Region "Classwide Variables"
    ' The Amino Acid arrays contain 26 entries, corresponding to A through Z
    ' Invalid/Undefined letters (J and U) have values of 0 for the mass and atom counts

    Private AMINO_ACID_LIST_MAX_INDEX As Integer = 25
    Private mAminoAcidMasses() As Double
    Private mAminoAcidAtomCounts() As udtAtomCountsType

    ' typically mPeptideNTerminusMass + mPeptideCTerminusMass = 18.0105633 (the mass of water)
    Private mPeptideNTerminusMass As Double
    Private mPeptideCTerminusMass As Double

    Private mRemovePrefixAndSuffixIfPresent As Boolean

    Private mErrorMessage As String
#End Region

#Region "Properties"
    Public ReadOnly Property ErrorMessage() As String
        Get
            Return mErrorMessage
        End Get
    End Property

    Public Property PeptideCTerminusMass() As Double
        Get
            Return mPeptideCTerminusMass
        End Get
        Set(ByVal Value As Double)
            mPeptideCTerminusMass = Value
        End Set
    End Property

    Public Property PeptideNTerminusMass() As Double
        Get
            Return mPeptideNTerminusMass
        End Get
        Set(ByVal Value As Double)
            mPeptideNTerminusMass = Value
        End Set
    End Property

    Public Property RemovePrefixAndSuffixIfPresent() As Boolean
        Get
            Return mRemovePrefixAndSuffixIfPresent
        End Get
        Set(ByVal Value As Boolean)
            mRemovePrefixAndSuffixIfPresent = Value
        End Set
    End Property
#End Region

    Public Sub New()
        mErrorMessage = String.Empty
        mRemovePrefixAndSuffixIfPresent = True
        InitializeAminoAcidData()
    End Sub

    Private Sub AddAminoAcidStatEntry(ByVal intIndex As Integer, ByVal dblMonoisotopicMass As Double, ByVal CountC As Integer, ByVal CountH As Integer, ByVal CountN As Integer, ByVal CountO As Integer, ByVal CountS As Integer)
        ' Adds new entry to AminoAcidMasses and AminoAcidSymbols

        mAminoAcidMasses(intIndex) = dblMonoisotopicMass

        With mAminoAcidAtomCounts(intIndex)
            .CountC = CountC
            .CountH = CountH
            .CountN = CountN
            .CountO = CountO
            .CountS = CountS
        End With
    End Sub

    Public Function ComputeSequenceMass(ByVal strSequence As String) As Double
        ' Computes the mass for sequence strSequence
        ' Returns -1 if an error

        Dim chChar As Char
        Dim intAAIndex As Integer
        Dim dblMass As Double = 0
        Dim intValidResidueCount As Short = 0

        If mRemovePrefixAndSuffixIfPresent Then
            If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequence, strSequence, Nothing, Nothing) Then
                mErrorMessage = "Error calling SplitPrefixAndSuffixFromSequence for " & strSequence
                Return -1
            End If
        End If

        mErrorMessage = String.Empty
        For Each chChar In strSequence
            ' Use Convert.ToInt32 to convert to the Ascii value, then subtract 65
            intAAIndex = System.Convert.ToInt32(chChar) - 65

            Try
                If intAAIndex < 0 OrElse intAAIndex > AMINO_ACID_LIST_MAX_INDEX Then
                    mErrorMessage = "Unknown symbol " & chChar & " in sequence " & strSequence
                    dblMass = -1
                    Exit For
                Else
                    dblMass += mAminoAcidMasses(intAAIndex)
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

    Public Function ComputeSequenceMass(ByVal strSequence As String, ByVal intModCount As Integer, ByRef udtResidueModificationInfo() As udtPeptideSequenceModInfoType) As Double
        ' Computes the mass for sequence strSequence using the mods in udtResidueModificationInfo()
        ' Returns -1 if an error

        Dim intIndex As Integer
        Dim dblMass As Double

        Dim udtAtomCounts As udtAtomCountsType
        Dim blnAtomCountsDefined As Boolean

        ' Note that ComputeSequenceMass will reset mErorMessage
        dblMass = ComputeSequenceMass(strSequence)

        If intModCount > 0 AndAlso dblMass >= 0 Then

            blnAtomCountsDefined = False
            For intIndex = 0 To intModCount - 1
                With udtResidueModificationInfo(intIndex)

                    If .AffectedAtom = Nothing OrElse .AffectedAtom = NO_AFFECTED_ATOM_SYMBOL Then
                        ' Positional modification (static or dynamic mod)
                        ' Simply add the modification mass to dblMass
                        dblMass += .ModificationMass
                    Else
                        ' Isotopic modification
                        If Not blnAtomCountsDefined Then
                            CountAtoms(strSequence, udtAtomCounts)
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

    Public Shared Function ConvoluteMass(ByVal dblMassMZ As Double, ByVal intCurrentCharge As Integer, Optional ByVal intDesiredCharge As Integer = 1) As Double
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
                ElseIf intCurrentCharge <= 0 Then
                    dblNewMZ = dblMassMZ + MASS_PROTON
                Else
                    dblNewMZ = (dblMassMZ * intCurrentCharge) - MASS_PROTON * (intCurrentCharge - 1)
                End If

                If intDesiredCharge > 1 Then
                    dblNewMZ = (dblNewMZ + MASS_PROTON * (intDesiredCharge - 1)) / intDesiredCharge
                ElseIf intDesiredCharge = 0 Then
                    dblNewMZ -= MASS_PROTON
                End If
            End If
        Catch ex As Exception
            ' Error occurred
            dblNewMZ = 0
        End Try

        Return dblNewMZ

    End Function

    Private Sub CountAtoms(ByRef strSequence As String, ByRef udtAtomCounts As udtAtomCountsType)

        Dim chChar As Char
        Dim intAAIndex As Integer

        With udtAtomCounts
            .CountC = 0
            .CountH = 0
            .CountN = 0
            .CountO = 0
            .CountS = 0
        End With

        For Each chChar In strSequence
            ' Use Convert.ToInt32 to convert to the Ascii value, then subtract 65
            intAAIndex = System.Convert.ToInt32(chChar) - 65

            Try
                If intAAIndex < 0 OrElse intAAIndex > AMINO_ACID_LIST_MAX_INDEX Then
                    mErrorMessage = "Unknown symbol " & chChar & " in sequence " & strSequence
                    Exit For
                Else
                    With udtAtomCounts
                        .CountC += mAminoAcidAtomCounts(intAAIndex).CountC
                        .CountH += mAminoAcidAtomCounts(intAAIndex).CountH
                        .CountN += mAminoAcidAtomCounts(intAAIndex).CountN
                        .CountO += mAminoAcidAtomCounts(intAAIndex).CountO
                        .CountS += mAminoAcidAtomCounts(intAAIndex).CountS
                    End With
                End If
            Catch ex As Exception
                ' Invalid value; ignore
            End Try
        Next chChar
    End Sub

    Public Function GetAminoAcidAtomCounts(ByVal chAminoAcidSymbol As Char) As udtAtomCountsType
        ' Returns the atom counts if success, 0 if an error

        Dim intAAIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            intAAIndex = System.Convert.ToInt32(chAminoAcidSymbol) - 65
            If intAAIndex < 0 OrElse intAAIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return New udtAtomCountsType
            Else
                Return mAminoAcidAtomCounts(intAAIndex)
            End If
        End If

    End Function

    Public Function GetAminoAcidMass(ByVal chAminoAcidSymbol As Char) As Double
        ' Returns the mass if success, 0 if an error

        Dim intAAIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            intAAIndex = System.Convert.ToInt32(chAminoAcidSymbol) - 65
            If intAAIndex < 0 OrElse intAAIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return 0
            Else
                Return mAminoAcidMasses(intAAIndex)
            End If
        End If
    End Function

    Private Sub InitializeAminoAcidData()
        ' These monoisotopic masses come from those traditionally used in DMS
        ' They were originally assembled by Gordon Anderson for use in ICR-2LS

        ReDim mAminoAcidMasses(AMINO_ACID_LIST_MAX_INDEX)
        ReDim mAminoAcidAtomCounts(AMINO_ACID_LIST_MAX_INDEX)

        AddAminoAcidStatEntry(0, 71.0371100902557, 3, 5, 1, 1, 0)           ' A
        AddAminoAcidStatEntry(1, 114.042921543121, 4, 6, 2, 2, 0)           ' B: use N or D (aka Asn/Asp)
        AddAminoAcidStatEntry(2, 103.009180784225, 3, 5, 1, 1, 1)           ' C
        AddAminoAcidStatEntry(3, 115.026938199997, 4, 5, 1, 3, 0)           ' D
        AddAminoAcidStatEntry(4, 129.042587518692, 5, 7, 1, 3, 0)           ' E
        AddAminoAcidStatEntry(5, 147.068408727646, 9, 9, 1, 1, 0)           ' F
        AddAminoAcidStatEntry(6, 57.0214607715607, 2, 3, 1, 1, 0)           ' G
        AddAminoAcidStatEntry(7, 137.058904886246, 6, 7, 3, 1, 0)           ' H
        AddAminoAcidStatEntry(8, 113.084058046341, 6, 11, 1, 1, 0)          ' I
        AddAminoAcidStatEntry(9, 0, 0, 0, 0, 0, 0)                          ' J: Could use mass of Ile/Leu, but we're using ""
        AddAminoAcidStatEntry(10, 128.094955444336, 6, 12, 2, 1, 0)         ' K
        AddAminoAcidStatEntry(11, 113.084058046341, 6, 11, 1, 1, 0)         ' L
        AddAminoAcidStatEntry(12, 131.040479421616, 5, 9, 1, 1, 1)          ' M
        AddAminoAcidStatEntry(13, 114.042921543121, 4, 6, 2, 2, 0)          ' N
        AddAminoAcidStatEntry(14, 114.079306125641, 5, 10, 2, 1, 0)         ' O
        AddAminoAcidStatEntry(15, 97.0527594089508, 5, 7, 1, 1, 0)          ' P
        AddAminoAcidStatEntry(16, 128.058570861816, 5, 8, 2, 2, 0)          ' Q
        AddAminoAcidStatEntry(17, 156.101100921631, 6, 12, 4, 1, 0)         ' R
        AddAminoAcidStatEntry(18, 87.0320241451263, 3, 5, 1, 2, 0)          ' S
        AddAminoAcidStatEntry(19, 101.047673463821, 4, 7, 1, 2, 0)          ' T
        AddAminoAcidStatEntry(20, 150.95363, 3, 5, 1, 1, 0)                 ' U: Corresponds to Sec = Selenocysteine (C3H5NOSe)
        AddAminoAcidStatEntry(21, 99.0684087276459, 5, 9, 1, 1, 0)          ' V
        AddAminoAcidStatEntry(22, 186.079306125641, 11, 10, 2, 1, 0)        ' W
        AddAminoAcidStatEntry(23, 113.084058046341, 6, 11, 1, 1, 0)         ' X: Unknown; use mass of Ile/Leu
        AddAminoAcidStatEntry(24, 163.063322782516, 9, 9, 1, 2, 0)          ' Y
        AddAminoAcidStatEntry(25, 128.058570861816, 5, 8, 2, 2, 0)          ' Z: use Q or E (aka Gln/Glu); note that these are 0.984 Da apart

        ResetTerminusMasses()
    End Sub

    Public Shared Function MassToPPM(ByVal dblMassToConvert As Double, ByVal dblCurrentMZ As Double) As Double
        ' Converts dblMassToConvert to ppm, based on the value of dblCurrentMZ

        Return dblMassToConvert * 1000000.0 / dblCurrentMZ
    End Function

    Public Shared Function MHToMonoisotopicMass(ByVal dblMH As Double) As Double
        Return ConvoluteMass(dblMH, 1, 0)
    End Function

    Public Shared Function MonoisotopicMassToMZ(ByVal dblMonoisotopicMass As Double, ByVal intDesiredCharge As Integer) As Double
        Return ConvoluteMass(dblMonoisotopicMass, 0, intDesiredCharge)
    End Function

    Public Shared Function PPMToMass(ByVal dblPPMToConvert As Double, ByVal dblCurrentMZ As Double) As Double
        ' Converts dblPPMToConvert to a mass value, which is dependent on dblCurrentMZ

        Return dblPPMToConvert / 1000000.0 * dblCurrentMZ
    End Function

    Public Sub ResetTerminusMasses()
        ' See comment in Sub InitializeAminoAcidData concerning these masses

        mPeptideNTerminusMass = DEFAULT_N_TERMINUS_MASS_CHANGE
        mPeptideCTerminusMass = DEFAULT_C_TERMINUS_MASS_CHANGE
    End Sub

    Public Function SetAminoAcidAtomCounts(ByVal chAminoAcidSymbol As Char, ByVal udtAtomCounts As udtAtomCountsType) As Boolean
        ' Returns True if success, False if an invalid amino acid symbol

        Dim intAAIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            intAAIndex = System.Convert.ToInt32(chAminoAcidSymbol) - 65
            If intAAIndex < 0 OrElse intAAIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return False
            Else
                mAminoAcidAtomCounts(intAAIndex) = udtAtomCounts
                Return True
            End If
        End If

    End Function

    Public Function SetAminoAcidMass(ByVal chAminoAcidSymbol As Char, ByVal dblMass As Double) As Boolean
        ' Returns True if success, False if an invalid amino acid symbol

        Dim intAAIndex As Integer

        If Not chAminoAcidSymbol = Nothing Then
            intAAIndex = System.Convert.ToInt32(chAminoAcidSymbol) - 65
            If intAAIndex < 0 OrElse intAAIndex > AMINO_ACID_LIST_MAX_INDEX Then
                ' Invalid Index
                Return False
            Else
                mAminoAcidMasses(intAAIndex) = dblMass
                Return True
            End If
        End If
    End Function
End Class

