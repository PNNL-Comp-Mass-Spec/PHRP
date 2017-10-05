Option Strict On

<CLSCompliant(True)>
Public Class clsAminoAcidModInfo

#Region "Constants and Enums"
    Public Const N_TERMINAL_PEPTIDE_SYMBOL_DMS As Char = "<"c
    Public Const C_TERMINAL_PEPTIDE_SYMBOL_DMS As Char = ">"c
    Public Const N_TERMINAL_PROTEIN_SYMBOL_DMS As Char = "["c
    Public Const C_TERMINAL_PROTEIN_SYMBOL_DMS As Char = "]"c

    Public Enum eResidueTerminusStateConstants As Integer
        None = 0                        ' The residue is in the middle of the peptide
        PeptideNTerminus = 1            ' The residue is located at the peptide's N-terminus; superseded by ProteinNTerminus if applicable
        PeptideCTerminus = 2            ' The residue is located at the peptide's C-terminus; superseded by ProteinCTerminus if applicable
        ProteinNTerminus = 3            ' The residue is located at the protein's N-terminus
        ProteinCTerminus = 4            ' The residue is located at the protein's C-terminus
        ProteinNandCCTerminus = 5       ' The protein only has one residue 
    End Enum
#End Region

    Private ReadOnly mModDefinition As clsModificationDefinition
    Private ReadOnly mResidue As Char
    Private ReadOnly mResidueLocInPeptide As Integer
    Private ReadOnly mEndResidueLocInPeptide As Integer
    Private ReadOnly mResidueTerminusState As eResidueTerminusStateConstants

    ''' <summary>
    ''' True if the location of the modification is ambiguous
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property AmbiguousMod As Boolean
        Get
            If mEndResidueLocInPeptide > mResidueLocInPeptide Then
                Return True
            Else
                Return False
            End If
        End Get
    End Property

    ''' <summary>
    ''' For ambiguous mods, indicates the last residue on which the mod could appear.  For non-ambiguous mods, whill be the same as ResidueLocInPeptide
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property EndResidueLocInPeptide As Integer
        Get
            Return mEndResidueLocInPeptide
        End Get
    End Property

    Public ReadOnly Property ModDefinition As clsModificationDefinition
        Get
            Return mModDefinition
        End Get
    End Property

    Public ReadOnly Property Residue As Char
        Get
            Return mResidue
        End Get
    End Property

    ''' <summary>
    ''' Indicates the residue number modified; the first residue is at position 1
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>For ambiguous mods, indicates the first residue on which the mod could appear</remarks>
    Public ReadOnly Property ResidueLocInPeptide As Integer
        Get
            Return mResidueLocInPeptide
        End Get
    End Property

    Public ReadOnly Property ResidueTerminusState As eResidueTerminusStateConstants
        Get
            Return mResidueTerminusState
        End Get
    End Property

    Public Sub New(
      chResidue As Char,
      intResidueLocInPeptide As Integer,
      eResidueTerminusState As eResidueTerminusStateConstants,
      objModDefinition As clsModificationDefinition)

        mModDefinition = objModDefinition
        mResidue = chResidue
        mResidueLocInPeptide = intResidueLocInPeptide
        mEndResidueLocInPeptide = mResidueLocInPeptide
        mResidueTerminusState = eResidueTerminusState

    End Sub

    Public Sub New(
      chResidue As Char,
      intResidueLocInPeptide As Integer,
      eResidueTerminusState As eResidueTerminusStateConstants,
      objModDefinition As clsModificationDefinition,
      intEndResidueLocInPeptide As Integer)

        mModDefinition = objModDefinition
        mResidue = chResidue
        mResidueLocInPeptide = intResidueLocInPeptide
        mEndResidueLocInPeptide = intEndResidueLocInPeptide
        mResidueTerminusState = eResidueTerminusState

    End Sub
End Class
