Option Strict On

<CLSCompliant(True)>
Public Class clsAminoAcidModInfo

#Region "Constants and Enums"
	Public Const N_TERMINAL_PEPTIDE_SYMBOL_DMS As Char = "<"c
	Public Const C_TERMINAL_PEPTIDE_SYMBOL_DMS As Char = ">"c
	Public Const N_TERMINAL_PROTEIN_SYMBOL_DMS As Char = "["c
	Public Const C_TERMINAL_PROTEIN_SYMBOL_DMS As Char = "]"c

	Public Enum eResidueTerminusStateConstants As Integer
		None = 0						' The residue is in the middle of the peptide
		PeptideNTerminus = 1			' The residue is located at the peptide's N-terminus; superseded by ProteinNTerminus if applicable
		PeptideCTerminus = 2			' The residue is located at the peptide's C-terminus; superseded by ProteinCTerminus if applicable
		ProteinNTerminus = 3			' The residue is located at the protein's N-terminus
		ProteinCTerminus = 4			' The residue is located at the protein's C-terminus
		ProteinNandCCTerminus = 5		' The protein only has one residue 
	End Enum
#End Region

	Protected mModDefinition As clsModificationDefinition
	Protected mResidue As Char
	Protected mResidueLocInPeptide As Integer							   ' Indicates the residue number modified; the first residue is at position 1
	Protected mResidueTerminusState As eResidueTerminusStateConstants

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

	Public Sub New(Residue As Char, ResidueLocInPeptide As Integer, ResidueTerminusState As eResidueTerminusStateConstants, ModDefinition As clsModificationDefinition)
		mModDefinition = ModDefinition
		mResidue = Residue
		mResidueLocInPeptide = ResidueLocInPeptide
		mResidueTerminusState = ResidueTerminusState
	End Sub
End Class
