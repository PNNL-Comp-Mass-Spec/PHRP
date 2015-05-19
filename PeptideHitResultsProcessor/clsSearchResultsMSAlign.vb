Option Strict On

' This class is used to track the peptide details for a MSAlign search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Created 11/27/2012
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------

Public Class clsSearchResultsMSAlign
	Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"
	' Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables 
#End Region

#Region "Properties"

	' Auto-Properties
	Public Property Prsm_ID As String
	Public Property Spectrum_ID As String

	Public Property Protein_Mass As String
	Public Property Unexpected_Mod_Count As String

	Public Property Peak_Count As String
	Public Property Matched_Peak_Count As String
	Public Property Matched_Fragment_Ion_Count As String

	Public Property PValue As String
	Public Property Rank_PValue As String

	Public Property EValue As String
	Public Property FDR As String

	Public Property Species_ID As String
	Public Property FragMethod As String

	Public Property Precursor_mz As String				' Observed precursor_mz
	Public Property MSAlignComputedDelM As String
	Public Property MSAlignComputedDelMPPM As String

#End Region

	Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
		' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
		MyBase.New(objPeptideMods)
	End Sub

	Public Overrides Sub Clear()
		MyBase.Clear()

		Prsm_ID = String.Empty
		Spectrum_ID = String.Empty

		Protein_Mass = String.Empty
		Unexpected_Mod_Count = String.Empty

		Peak_Count = String.Empty
		Matched_Peak_Count = String.Empty
		Matched_Fragment_Ion_Count = String.Empty

		PValue = String.Empty
		Rank_PValue = String.Empty

		EValue = String.Empty
		FDR = String.Empty

		Species_ID = String.Empty
		FragMethod = String.Empty

		Precursor_mz = String.Empty
		MSAlignComputedDelM = String.Empty
		MSAlignComputedDelMPPM = String.Empty

	End Sub
End Class
