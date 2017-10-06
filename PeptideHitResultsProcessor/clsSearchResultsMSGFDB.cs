Option Strict On

' This class is used to track the peptide details for an MSGF+ search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Created 8/12/2011
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------

Imports PHRPReader

Public Class clsSearchResultsMSGFDB
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"
    ' Note that "Automatic properties" are being used; thus, we don't need to explicitly define class variables 
#End Region

#Region "Properties"

    Public Property FragMethod As String
    Public Property NTT As String

    Public Property DeNovoScore As String
    Public Property MSGFScore As String
    Public Property SpecEValue As String        ' SpecProb in MSGFDB; SpecEValue in MSGF+
    Public Property RankSpecEValue As String
    Public Property EValue As String            ' PValue in MSGFDB; EValue in MSGF+

    Public Property QValue As String            ' Will contain target/decoy FDR when -tda 1 was used; will contain EFDR when -tda 1 was not used; FDR in MSGFDB; QValue in MSGF+
    Public Property PepQValue As String         ' Only present if searched using -tda 1; PepFDR in MSGFDB; PepQValue in MSGF+

    Public Property PrecursorMZ As String
    Public Property MSGFPlusComputedDelM As String
    Public Property MSGFPlusComputedDelMPPM As String

    Public Property IsotopeError As String      ' Only reported by MSGF+

    Public Property MSGFPlusResults As Boolean
#End Region

    Public Sub New(objPeptideMods As clsPeptideModificationContainer, peptideSeqMassCalculator As clsPeptideMassCalculator)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods, peptideSeqMassCalculator)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        FragMethod = String.Empty
        NTT = String.Empty

        DeNovoScore = String.Empty
        MSGFScore = String.Empty
        SpecEValue = String.Empty
        RankSpecEValue = String.Empty
        EValue = String.Empty

        QValue = String.Empty
        PepQValue = String.Empty

        PrecursorMZ = String.Empty
        MSGFPlusComputedDelM = String.Empty
        MSGFPlusComputedDelMPPM = String.Empty

        IsotopeError = String.Empty
        MSGFPlusResults = False
    End Sub

End Class
