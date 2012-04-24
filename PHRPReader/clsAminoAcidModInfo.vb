Option Strict On

Imports PeptideHitResultsProcessor.clsPeptideModificationContainer

Public Class clsAminoAcidModInfo
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
