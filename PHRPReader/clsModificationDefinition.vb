' This class describes an amino acid modification
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.
' Started January 6, 2006
' Last updated June 26, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

Option Strict On

<Assembly: CLSCompliant(True)>

<CLSCompliant(True)>
Public Class clsModificationDefinition

#Region "Constants and Enums"
    Public Const LAST_RESORT_MODIFICATION_SYMBOL As Char = "_"c                          ' The underscore is used if all of the DEFAULT_MODIFICATION_SYMBOLS are used up
    Public Const NO_SYMBOL_MODIFICATION_SYMBOL As Char = "-"c
    Public Const UNKNOWN_MOD_BASE_NAME As String = "UnkMod"
    Public Const INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME As String = UNKNOWN_MOD_BASE_NAME & "00"

    Public Enum eModificationTypeConstants As Integer
        ''' <summary>
        ''' Unknown mod type on a residue; essentially treated as a dynamic mod
        ''' </summary>
        ''' <remarks></remarks>
        UnknownType = 0

        ''' <summary>
        ''' Dynamic mod on a residue or peptide terminus; supported by Sequest and notated via a modification symbol; this mod is explicitly notated by X!Tandem; if a terminus mod, then the mod symbol is associated with the first or last residue in the peptide
        ''' </summary>
        ''' <remarks></remarks>
        DynamicMod = 1

        ''' <summary>
        ''' Static mod on a residue or peptide terminus; supported by Sequest but not explicitly notated; this mod is explicitly notated by X!Tandem; if a terminus mod, then the mod symbol is associated with the first or last residue in the peptide
        ''' </summary>
        ''' <remarks></remarks>
        StaticMod = 2

        ''' <summary>
        ''' Peptide terminus static mod (DMS Symbol is T); used by Sequest and MSGFDB; note that terminal mods are always dynamic in X!Tandem
        ''' </summary>
        ''' <remarks></remarks>
        TerminalPeptideStaticMod = 3

        ''' <summary>
        ''' Isotopic mod, e.g. N15, or C13; supported by Sequest; most likely not supported by XTandem
        ''' </summary>
        ''' <remarks></remarks>
        IsotopicMod = 4

        ''' <summary>
        ''' Protein terminus static mod; supported by Sequest; this mod is also supported by X!Tandem but modified residues are not explicitly notated; instead, all peptides have their mass implicitly modified by this amount
        ''' </summary>
        ''' <remarks></remarks>
        ProteinTerminusStaticMod = 5
    End Enum

#End Region

#Region "Classwide Variables"
    Private mModificationSymbol As Char                 ' One letter symbol for this modification; use NO_SYMBOL_MODIFICATION_SYMBOL if no symbol (necessary for isotopic mods or protein terminus static mods)
    Private mModificationMass As Double                 ' Monoisotopic modification mass
    Private mModificationMassAsText As String           ' Modification mass, stored as text
    Private mTargetResidues As String                   ' If this string is empty, then the given modification can apply to any residue or terminus; Otherwise, should contain a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to; Use the *_SYMBOL_DMS constants for the peptide and protein terminii symbols (< and > for the peptide terminii; [ and ] for the protein terminii)
    Private mModificationType As eModificationTypeConstants
    Private mMassCorrectionTag As String                ' Name associated with the given ModificationMass; maximum length is 8 characters; cannot contain a colon, comma, or space
    Private mAffectedAtom As Char                       ' Set to Nothing or to clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL for positional modifications (including terminus modifications); for Isotopic modifications, indicate the atom affected (e.g. C, H, N, O, or S)
    Private mOccurrenceCount As Integer                 ' Number of times this modification was observed in the given XML file
    Private mUnknownModAutoDefined As Boolean           ' True if this was an unknown mass that was auto defined
#End Region

#Region "Properties"

    ''' <summary>
    ''' One letter symbol for this modification
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>
    ''' Use NO_SYMBOL_MODIFICATION_SYMBOL (a dash) if no symbol
    ''' (necessary for isotopic mods or protein terminus static mods)
    ''' </remarks>
    Public Property ModificationSymbol As Char
        Get
            Return mModificationSymbol
        End Get
        Set
            If Value = Nothing Then
                mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL
            Else
                mModificationSymbol = Value
            End If
        End Set
    End Property

    ''' <summary>
    ''' Monoisotopic modification mass
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ModificationMass As Double
        Get
            Return mModificationMass
        End Get
        Set
            mModificationMass = Value
            mModificationMassAsText = mModificationMass.ToString()
        End Set
    End Property

    ''' <summary>
    ''' Modification mass, stored as text
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Represents the original string value read from the data file</remarks>
    Public Property ModificationMassAsText As String
        Get
            Return mModificationMassAsText
        End Get
        Set
            mModificationMassAsText = Value
        End Set
    End Property

    ''' <summary>
    ''' Residues that this modification can apply to
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>
    ''' If an empty string, then the modification can apply to any residue or terminus;
    ''' Otherwise, should contain a space-free, comma-free list of one letter amino acid residue symbols that this mod can apply to.
    ''' Use the *_SYMBOL_DMS constants for the peptide and protein terminii symbols
    ''' (less than and greater than signs for the peptide terminii; [ and ] for the protein terminii)
    ''' </remarks>
    Public Property TargetResidues As String
        Get
            Return mTargetResidues
        End Get
        Set
            If Value Is Nothing Then
                mTargetResidues = String.Empty
            Else
                mTargetResidues = String.Copy(Value)
            End If
        End Set
    End Property

    ''' <summary>
    ''' Modification type
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ModificationType As eModificationTypeConstants
        Get
            Return mModificationType
        End Get
        Set
            mModificationType = Value
        End Set
    End Property

    ''' <summary>
    ''' Modification name, for example Phosph, IodoAcet, Plus1Oxy, or Methyl
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Maximum length is 8 characters; cannot contain a colon, comma, or space</remarks>
    Public Property MassCorrectionTag As String
        Get
            Return mMassCorrectionTag
        End Get
        Set
            If Value Is Nothing Then
                mMassCorrectionTag = String.Empty
            Else
                mMassCorrectionTag = String.Copy(Value)
            End If
        End Set
    End Property

    ''' <summary>
    ''' Only used with Isotopic modifications, indicating the atom affected (e.g. C, H, N, O, or S)
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>
    ''' Set to Nothing or to clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL (a dash) for positional modifications
    ''' (including terminus modifications)
    ''' </remarks>
    Public Property AffectedAtom As Char
        Get
            Return mAffectedAtom
        End Get
        Set
            If Value = Nothing Then
                mAffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
            Else
                mAffectedAtom = Value
            End If
        End Set
    End Property

    ''' <summary>
    ''' Number of times this modification was observed in the given dataset
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property OccurrenceCount As Integer
        Get
            Return mOccurrenceCount
        End Get
        Set
            mOccurrenceCount = Value
        End Set
    End Property

    ''' <summary>
    ''' True if this was an unknown mass that was auto defined
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property UnknownModAutoDefined As Boolean
        Get
            Return mUnknownModAutoDefined
        End Get
        Set
            mUnknownModAutoDefined = Value
        End Set
    End Property
#End Region

    Public Sub New()
        Me.Clear()
    End Sub

    Public Sub New(chModificationSymbol As Char, dblModificationMass As Double)
        Me.Clear()

        Me.ModificationSymbol = chModificationSymbol
        Me.ModificationMass = dblModificationMass
    End Sub

    Public Sub New(dblModificationMass As Double, strTargetResidues As String, eModificationType As eModificationTypeConstants)
        Me.Clear()

        Me.ModificationMass = dblModificationMass
        Me.TargetResidues = strTargetResidues
        Me.ModificationType = eModificationType
    End Sub

    Public Sub New(chModificationSymbol As Char, dblModificationMass As Double, strTargetResidues As String, eModificationType As eModificationTypeConstants, strMassCorrectionTag As String)
        Me.Clear()

        Me.ModificationSymbol = chModificationSymbol
        Me.ModificationMass = dblModificationMass
        Me.TargetResidues = strTargetResidues
        Me.ModificationType = eModificationType
        Me.MassCorrectionTag = strMassCorrectionTag
    End Sub

    Public Sub New(chModificationSymbol As Char, dblModificationMass As Double, strTargetResidues As String, eModificationType As eModificationTypeConstants, strMassCorrectionTag As String, chAffectedAtom As Char, blnUnknownModAutoDefined As Boolean)
        Me.Clear()

        Me.ModificationSymbol = chModificationSymbol
        Me.ModificationMass = dblModificationMass
        Me.TargetResidues = strTargetResidues
        Me.ModificationType = eModificationType
        Me.MassCorrectionTag = strMassCorrectionTag
        Me.AffectedAtom = chAffectedAtom
        Me.UnknownModAutoDefined = blnUnknownModAutoDefined
    End Sub

    ''' <summary>
    ''' Initialize the modification definition
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub Clear()
        mModificationSymbol = NO_SYMBOL_MODIFICATION_SYMBOL
        mModificationMass = 0
        mModificationMassAsText = "0"
        mTargetResidues = String.Empty
        mModificationType = eModificationTypeConstants.UnknownType
        mMassCorrectionTag = INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME
        mAffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
        mOccurrenceCount = 0
        mUnknownModAutoDefined = False
    End Sub

    ''' <summary>
    ''' Compares objB to this object, ignoring .ModificationSymbol and ignoring .AffectedResidues
    ''' </summary>
    ''' <param name="objB"></param>
    ''' <returns>True if the items are equivalent</returns>
    ''' <remarks></remarks>
    Public Function EquivalentMassTypeTagAndAtom(objB As clsModificationDefinition) As Boolean
        Return EquivalentMassTypeTagAndAtom(Me, objB)
    End Function

    ''' <summary>
    ''' Compare objA to objB but ignore .ModificationSymbol and .AffectedResidues
    ''' </summary>
    ''' <param name="objA"></param>
    ''' <param name="objB"></param>
    ''' <returns>True if the items are equivalent</returns>
    ''' <remarks></remarks>
    Public Function EquivalentMassTypeTagAndAtom(objA As clsModificationDefinition, objB As clsModificationDefinition) As Boolean
        Dim blnEquivalent As Boolean

        blnEquivalent = False
        With objA
            If Math.Abs(Math.Round(.ModificationMass - objB.ModificationMass, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) - 0) < Single.Epsilon AndAlso
               .ModificationType = objB.ModificationType AndAlso
               .MassCorrectionTag = objB.MassCorrectionTag AndAlso
               .AffectedAtom = objB.AffectedAtom Then
                blnEquivalent = True
            End If
        End With

        Return blnEquivalent
    End Function

    ''' <summary>
    ''' Compares objB to this object, ignoring .ModificationSymbol
    ''' </summary>
    ''' <param name="objB"></param>
    ''' <returns>True if the items are equivalent</returns>
    ''' <remarks></remarks>
    Public Function EquivalentMassTypeTagAtomAndResidues(objB As clsModificationDefinition) As Boolean
        Return EquivalentMassTypeTagAtomAndResidues(Me, objB)
    End Function

    ''' <summary>
    ''' Compares objB to this object
    ''' </summary>
    ''' <param name="objA"></param>
    ''' <param name="objB"></param>
    ''' <returns>True if the items are equivalent</returns>
    ''' <remarks></remarks>
    Public Function EquivalentMassTypeTagAtomAndResidues(objA As clsModificationDefinition, objB As clsModificationDefinition) As Boolean
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

    ''' <summary>
    ''' Compare the residue lists (ignoring order)
    ''' </summary>
    ''' <param name="strResidues1"></param>
    ''' <param name="strResidues2"></param>
    ''' <param name="blnAllowResidues2ToBeSubsetOfResidues1"></param>
    ''' <returns>True if they contain the same residues</returns>
    ''' <remarks></remarks>
    Public Shared Function EquivalentTargetResidues(strResidues1 As String, strResidues2 As String, blnAllowResidues2ToBeSubsetOfResidues1 As Boolean) As Boolean
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

    ''' <summary>
    ''' Returns True if this modification can affect the peptide or protein terminus
    ''' Note that some modifications can affect either peptide teriminii or internal residues
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function CanAffectPeptideOrProteinTerminus() As Boolean
        Dim lstTerminalSymbols As SortedSet(Of Char) = GetTerminalSymbols()

        If mModificationType = eModificationTypeConstants.ProteinTerminusStaticMod OrElse mModificationType = eModificationTypeConstants.TerminalPeptideStaticMod Then
            Return True
        Else
            For Each chChar As Char In mTargetResidues
                If lstTerminalSymbols.Contains(chChar) Then
                    Return True
                End If
            Next
            Return False
        End If

    End Function

    ''' <summary>
    ''' Returns true if this modification can affect peptide residues
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function CanAffectPeptideResidues() As Boolean
        Dim lstTerminalSymbols As SortedSet(Of Char) = GetTerminalSymbols()

        If mModificationType = eModificationTypeConstants.ProteinTerminusStaticMod OrElse mModificationType = eModificationTypeConstants.TerminalPeptideStaticMod Then
            Return False
        Else
            If String.IsNullOrEmpty(mTargetResidues) Then
                Return True
            Else
                For Each chChar As Char In mTargetResidues
                    If Not lstTerminalSymbols.Contains(chChar) Then
                        Return True
                    End If
                Next
            End If

            Return False
        End If

    End Function

    ''' <summary>
    ''' Retrieve the protein and peptide terminus symbols
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function GetTerminalSymbols() As SortedSet(Of Char)
        Dim lstTerminalSymbols As SortedSet(Of Char)

        lstTerminalSymbols = New SortedSet(Of Char)

        lstTerminalSymbols.Add(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)
        lstTerminalSymbols.Add(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS)
        lstTerminalSymbols.Add(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)
        lstTerminalSymbols.Add(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)

        Return lstTerminalSymbols

    End Function

    ''' <summary>
    ''' Retrieve the modification type for the given modification type symbol
    ''' </summary>
    ''' <param name="chModificationTypeSymbol">D, S, T, I, or P</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function ModificationSymbolToModificationType(chModificationTypeSymbol As Char) As clsModificationDefinition.eModificationTypeConstants
        If chModificationTypeSymbol = Nothing Then
            Return clsModificationDefinition.eModificationTypeConstants.UnknownType
        Else
            Select Case chModificationTypeSymbol
                Case "D"c
                    Return clsModificationDefinition.eModificationTypeConstants.DynamicMod
                Case "S"c
                    Return clsModificationDefinition.eModificationTypeConstants.StaticMod
                Case "T"c
                    Return clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                Case "I"c
                    Return clsModificationDefinition.eModificationTypeConstants.IsotopicMod
                Case "P"c
                    Return clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
                Case Else
                    Return clsModificationDefinition.eModificationTypeConstants.UnknownType
            End Select
        End If
    End Function

    ''' <summary>
    ''' Retrieve the modification type symbol for the given modification Type
    ''' </summary>
    ''' <param name="eModificationType"></param>
    ''' <returns>D, S, T, I, or P</returns>
    ''' <remarks></remarks>
    Public Shared Function ModificationTypeToModificationSymbol(eModificationType As clsModificationDefinition.eModificationTypeConstants) As Char
        Select Case eModificationType
            Case eModificationTypeConstants.DynamicMod
                Return "D"c
            Case eModificationTypeConstants.StaticMod
                Return "S"c
            Case eModificationTypeConstants.TerminalPeptideStaticMod
                Return "T"c
            Case eModificationTypeConstants.IsotopicMod
                Return "I"c
            Case eModificationTypeConstants.ProteinTerminusStaticMod
                Return "P"c
            Case Else
                Return "?"c
        End Select
    End Function

    ''' <summary>
    ''' Check whether the target residues contain the given residue
    ''' </summary>
    ''' <param name="chComparisonResidue"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function TargetResiduesContain(chComparisonResidue As Char) As Boolean
        If chComparisonResidue = Nothing Then
            Return False
        ElseIf Me.TargetResidues.IndexOf(chComparisonResidue) >= 0 Then
            Return True
        Else
            Return False
        End If
    End Function
End Class
