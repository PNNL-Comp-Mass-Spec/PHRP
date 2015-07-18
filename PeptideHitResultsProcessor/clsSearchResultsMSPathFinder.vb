Option Strict On

' This class is used to track the peptide details for a MSPathFinder search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Created 7/16/2015
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------

Public Class clsSearchResultsMSPathFinder
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"
    ' Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables 
#End Region

#Region "Properties"

    ' Auto-Properties
    Public Property MostAbundantIsotopeMz As String
    Public Property Composition As String
    Public Property ProteinDesc As String
    Public Property ProteinLength As String
    Public Property ResidueStart As String
    Public Property ResidueEnd As String
    Public Property MatchedFragments As String
    Public Property QValue As String
    Public Property PepQValue As String

#End Region

    Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        MostAbundantIsotopeMz = String.Empty
        Composition = String.Empty
        ProteinDesc = String.Empty
        ProteinLength = String.Empty
        ResidueStart = String.Empty
        ResidueEnd = String.Empty
        MatchedFragments = String.Empty
        QValue = String.Empty
        PepQValue = String.Empty

    End Sub



End Class
