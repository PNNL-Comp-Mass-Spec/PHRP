Option Explicit On 

' This class can be used as a base class for peptide hit results processor classes
'
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
' Started January 6, 2006
' Last updated March 6, 2006

Public Class clsModificationDefinition

#Region "Constants and Enums"
    Public Const LAST_RESORT_MODIFICATION_SYMBOL As Char = "_"c                          ' The underscore is used if all of the DEFAULT_MODIFICATION_SYMBOLS are used up
    Public Const NO_SYMBOL_MODIFICATION_SYMBOL As Char = "-"c
    Public Const UNKNOWN_MOD_BASE_NAME As String = "UnkMod"
    Public Const INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME As String = UNKNOWN_MOD_BASE_NAME & "00"

    Public Enum eModificationTypeConstants As Integer
        UnknownType = 0                         ' Unknown mod type on a residue; essentially treated as a dynamic mod
        DynamicMod = 1                          ' Dynamic mod on a residue or peptide terminus; this mod is explicitly notated by XTandem; if a terminus mod, then the mod symbol is associated with the first or last residue in the peptide
        StaticMod = 2                           ' Static mod on a residue or peptide terminus;  this mod is explicitly notated by XTandem; if a terminus mod, then the mod symbol is associated with the first or last residue in the peptide
        TerminalPeptideStaticMod = 3            ' Peptide terminus static mod (DMS Symbol is T); only used by Sequest since terminus mods are always dynamic in XTandem
        IsotopicMod = 4                         ' e.g. N15, or C13; most likely not supported by XTandem
        ProteinTerminusStaticMod = 5            ' Protein terminus static mod; this mod is supported by XTandem but modified residues are not explicitly notated; instead, all peptides have their mass implicitly modified by this amount
    End Enum

    Protected Const MASS_DIGITS_OF_PRECISION As Integer = 3

#End Region

#Region "Classwide Variables"
    Protected mModificationSymbol As Char               ' One letter symbol for this modification; use NO_SYMBOL_MODIFICATION_SYMBOL if no symbol (necessary for isotopic mods or protein terminus static mods)
    Protected mModificationMass As Double               ' Monoisotopic modification mass
    Protected mTargetResidues As String                 ' If this string is empty, then the given modification can apply to any residue or terminus; Otherwise, should contain a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to; Use the *_SYMBOL_DMS constants for the peptide and protein terminii symbols (< and > for the peptide terminii; [ and ] for the protein terminii)
    Protected mModificationType As eModificationTypeConstants
    Protected mMassCorrectionTag As String              ' Name associated with the given ModificationMass; maximum length is 8 characters; cannot contain a colon, comma, or space
    Protected mAffectedAtom As Char                     ' Set to Nothing or to clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications); for Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
    Protected mOccurrenceCount As Integer               ' Number of times this modification was observed in the given XML file
    Protected mUnknownModAutoDefined As Boolean         ' True if this was an unknown mass that was auto defined
#End Region

#Region "Properties"
    Public Property ModificationSymbol() As Char
        Get
            Return mModificationSymbol
        End Get
        Set(ByVal Value As Char)
            If Value = Nothing Then
                mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL
            Else
                mModificationSymbol = Value
            End If
        End Set
    End Property
    Public Property ModificationMass() As Double
        Get
            Return mModificationMass
        End Get
        Set(ByVal Value As Double)
            mModificationMass = Value
        End Set
    End Property
    Public Property TargetResidues() As String
        Get
            Return mTargetResidues
        End Get
        Set(ByVal Value As String)
            If Value Is Nothing Then
                mTargetResidues = String.Empty
            Else
                mTargetResidues = String.Copy(Value)
            End If
        End Set
    End Property
    Public Property ModificationType() As eModificationTypeConstants
        Get
            Return mModificationType
        End Get
        Set(ByVal Value As eModificationTypeConstants)
            mModificationType = Value
        End Set
    End Property
    Public Property MassCorrectionTag() As String
        Get
            Return mMassCorrectionTag
        End Get
        Set(ByVal Value As String)
            If Value Is Nothing Then
                mMassCorrectionTag = String.Empty
            Else
                mMassCorrectionTag = String.Copy(Value)
            End If
        End Set
    End Property
    Public Property AffectedAtom() As Char
        Get
            Return mAffectedAtom
        End Get
        Set(ByVal Value As Char)
            If Value = Nothing Then
                mAffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
            Else
                mAffectedAtom = Value
            End If
        End Set
    End Property
    Public Property OccurrenceCount() As Integer
        Get
            Return mOccurrenceCount
        End Get
        Set(ByVal Value As Integer)
            mOccurrenceCount = Value
        End Set
    End Property
    Public Property UnknownModAutoDefined() As Boolean
        Get
            Return mUnknownModAutoDefined
        End Get
        Set(ByVal Value As Boolean)
            mUnknownModAutoDefined = Value
        End Set
    End Property
#End Region

    Public Sub New()
        Me.Clear()
    End Sub

    Public Sub New(ByVal chModificationSymbol As Char, ByVal dblModificationMass As Double)
        Me.Clear()

        Me.ModificationSymbol = chModificationSymbol
        Me.ModificationMass = dblModificationMass
    End Sub
    Public Sub New(ByVal dblModificationMass As Double, ByVal strTargetResidues As String, ByVal eModificationType As eModificationTypeConstants)
        Me.Clear()

        Me.ModificationMass = dblModificationMass
        Me.TargetResidues = strTargetResidues
        Me.ModificationType = eModificationType
    End Sub
    Public Sub New(ByVal chModificationSymbol As Char, ByVal dblModificationMass As Double, ByVal strTargetResidues As String, ByVal eModificationType As eModificationTypeConstants, ByVal strMassCorrectionTag As String)
        Me.Clear()

        Me.ModificationSymbol = chModificationSymbol
        Me.ModificationMass = dblModificationMass
        Me.TargetResidues = strTargetResidues
        Me.ModificationType = eModificationType
        Me.MassCorrectionTag = strMassCorrectionTag
    End Sub
    Public Sub New(ByVal chModificationSymbol As Char, ByVal dblModificationMass As Double, ByVal strTargetResidues As String, ByVal eModificationType As eModificationTypeConstants, ByVal strMassCorrectionTag As String, ByVal chAffectedAtom As Char, ByVal blnUnknownModAutoDefined As Boolean)
        Me.Clear()

        Me.ModificationSymbol = chModificationSymbol
        Me.ModificationMass = dblModificationMass
        Me.TargetResidues = strTargetResidues
        Me.ModificationType = eModificationType
        Me.MassCorrectionTag = strMassCorrectionTag
        Me.AffectedAtom = chAffectedAtom
        Me.UnknownModAutoDefined = blnUnknownModAutoDefined
    End Sub

    Public Sub Clear()
        mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL
        mModificationMass = 0
        mTargetResidues = String.Empty
        mModificationType = eModificationTypeConstants.UnknownType
        mMassCorrectionTag = INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME
        mAffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
        mOccurrenceCount = 0
        mUnknownModAutoDefined = False
    End Sub

    Public Function EquivalentMassTypeTagAndAtom(ByVal objB As clsModificationDefinition) As Boolean
        ' Compares objB to this object, ignoring .ModificationSymbol and ignoring .AffectedResidues
        Return EquivalentMassTypeTagAndAtom(Me, objB)
    End Function

    Public Function EquivalentMassTypeTagAndAtom(ByVal objA As clsModificationDefinition, ByVal objB As clsModificationDefinition) As Boolean
        Dim blnEquivalent As Boolean

        ' Compare objA to objB but ignore .ModificationSymbol and .AffectedResidues
        blnEquivalent = False
        With objA
            If Math.Round(.ModificationMass - objB.ModificationMass, MASS_DIGITS_OF_PRECISION) = 0 AndAlso _
               .ModificationType = objB.ModificationType AndAlso _
               .MassCorrectionTag = objB.MassCorrectionTag AndAlso _
               .AffectedAtom = objB.AffectedAtom Then
                blnEquivalent = True
            End If
        End With

        Return blnEquivalent
    End Function

    Public Function EquivalentMassTypeTagAtomAndResidues(ByVal objB As clsModificationDefinition) As Boolean
        ' Compares objB to this object, ignoring .ModificationSymbol
        Return EquivalentMassTypeTagAtomAndResidues(Me, objB)
    End Function

    Public Function EquivalentMassTypeTagAtomAndResidues(ByVal objA As clsModificationDefinition, ByVal objB As clsModificationDefinition) As Boolean
        Dim blnEquivalent As Boolean

        ' First compare objA to objB but ignore .ModificationSymbol and .AffectedResidues
        blnEquivalent = EquivalentMassTypeTagAndAtom(objA, objB)

        If blnEquivalent Then
            ' Mass, ModificationType, MassCorrectionTag, and AffectedAtom are identical
            ' What about the residues?

            blnEquivalent = False
            With objA
                If .TargetResidues Is Nothing AndAlso objB.TargetResidues Is Nothing Then
                    blnEquivalent = True
                ElseIf Not .TargetResidues Is Nothing AndAlso Not objB.TargetResidues Is Nothing Then
                    If .ModificationType = eModificationTypeConstants.DynamicMod OrElse _
                       .ModificationType = eModificationTypeConstants.StaticMod Then
                        ' Matching dynamic or static modification definitions
                        ' Make sure each of the residues in objB.TargetResidues is present in .TargetResidues
                        If EquivalentTargetResidues(.TargetResidues, objB.TargetResidues, False) Then
                            blnEquivalent = True
                        End If
                    Else
                        ' Not a dynamic or static mod; require identical target residue lists in order to flag them as identical
                        If .TargetResidues = objB.TargetResidues Then
                            blnEquivalent = True
                        End If
                    End If
                End If
            End With
        End If

        Return blnEquivalent
    End Function

    Public Shared Function EquivalentTargetResidues(ByVal strResidues1 As String, ByVal strResidues2 As String, ByVal blnAllowResidues2ToBeSubsetOfResidues1 As Boolean) As Boolean
        Dim chChar As Char
        Dim intMatchCount As Integer
        Dim blnEquivalent As Boolean

        blnEquivalent = False

        If strResidues1 Is Nothing AndAlso strResidues2 Is Nothing Then
            ' Both residues lists are blank
            blnEquivalent = True
        ElseIf Not (strResidues1 Is Nothing OrElse strResidues2 Is Nothing) Then
            If strResidues1.Length = 0 AndAlso strResidues2.Length = 0 Then
                blnEquivalent = True
            ElseIf strResidues1 = strResidues2 Then
                blnEquivalent = True
            ElseIf strResidues1.Length >= strResidues2.Length Then
                ' See if each of the residues in strResidues2 is in strResidues1
                intMatchCount = 0
                For Each chChar In strResidues2
                    If strResidues1.IndexOf(chChar) >= 0 Then
                        intMatchCount += 1
                    Else
                        Exit For
                    End If
                Next chChar
                If intMatchCount = strResidues1.Length Then
                    blnEquivalent = True
                ElseIf blnAllowResidues2ToBeSubsetOfResidues1 AndAlso intMatchCount > 0 Then
                    blnEquivalent = True
                End If
            End If
        End If

        Return blnEquivalent

    End Function

    Public Function TargetResiduesContain(ByVal chComparisonResidue As Char) As Boolean
        If chComparisonResidue = Nothing Then
            Return False
        ElseIf Me.TargetResidues.IndexOf(chComparisonResidue) >= 0 Then
            Return True
        Else
            Return False
        End If
    End Function
End Class
