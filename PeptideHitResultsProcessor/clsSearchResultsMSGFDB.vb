Option Strict On

' This class can be used to track the peptide details for an MSGF-DB search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Created 8/12/2011
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------

Public Class clsSearchResultsMSGFDB
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"
    ' Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables 
#End Region

#Region "Properties"

    ' Auto-Properties (only work with Visual Studio 2010)
    Public Property FragMethod As String
    Public Property NTT As String

    Public Property DeNovoScore As String
    Public Property MSGFScore As String
    Public Property SpecProb As String
    Public Property RankSpecProb As String
    Public Property PValue As String

    Public Property FDR As String           ' Only present if searched using -tda 1
    Public Property PepFDR As String        ' Only present if searched using -tda 1

    Public Property PrecursorMZ As String
    Public Property MSGFDbComputedDelM As String
    Public Property MSGFDbComputedDelMPPM As String


#End Region

    Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        FragMethod = String.Empty
        NTT = String.Empty

        DeNovoScore = String.Empty
        MSGFScore = String.Empty
        SpecProb = String.Empty
        RankSpecProb = String.Empty
        PValue = String.Empty

        FDR = String.Empty
        PepFDR = String.Empty

        PrecursorMZ = String.Empty
        MSGFDbComputedDelM = String.Empty
        MSGFDbComputedDelMPPM = String.Empty
    End Sub

End Class
