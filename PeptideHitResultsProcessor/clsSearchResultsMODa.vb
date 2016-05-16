Option Strict On

' This class is used to track the peptide details for a MODa search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Created 4/01/2014
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------

Imports PHRPReader

Public Class clsSearchResultsMODa
	Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"
	' Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables 
#End Region

#Region "Properties"

	' Auto-Properties	
	Public Property Spectrum_Index As String

	Public Property Precursor_mz As String				' Observed precursor m/z, converted to monoisotopic mass by MODa

	Public Property MODaComputedDelM As String
	Public Property MODaComputedDelMPPM As String

	Public Property MODaScore As String

	Public Property Probability As String

#End Region

    Public Sub New(objPeptideMods As clsPeptideModificationContainer, peptideSeqMassCalculator As clsPeptideMassCalculator)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods, peptideSeqMassCalculator)
    End Sub

	Public Overrides Sub Clear()
		MyBase.Clear()
		Spectrum_Index = String.Empty

		Precursor_mz = String.Empty

		MODaComputedDelM = String.Empty
		MODaComputedDelMPPM = String.Empty

		MODaScore = String.Empty

		Probability = String.Empty
	End Sub
End Class
