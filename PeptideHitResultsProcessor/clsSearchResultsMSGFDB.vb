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

        PrecursorMZ = String.Empty
        MSGFDbComputedDelM = String.Empty
        MSGFDbComputedDelMPPM = String.Empty
    End Sub

    '' Unused
    ''Protected Overrides Sub ComputeDelMCorrected()


    ''    If Not String.IsNullOrEmpty(MSGFDbComputedDelMPPM) AndAlso Double.TryParse(MSGFDbComputedDelMPPM, mPeptideDeltaMassCorrectedPpm) Then
    ''        ' Simply use the DelM_PPM value computed by MSGF; column PMError(ppm)

    ''    Else
    ''        Dim dblDelM As Double
    ''        Dim intCorrectionCount As Integer = 0

    ''        Dim dblPrecursorMZ As Double
    ''        Dim intCharge As Integer
    ''        Dim dblPrecursorMonoMass As Double

    ''        Dim blnParseError As Boolean = False

    ''        ' Note that mPeptideDeltaMass is the DeltaMass value reported by MSGFDB 
    ''        ' (though clsMSGFDBResultsProcessor took the negative of the value in the results file so it currently represents "theoretical - observed")
    ''        If Double.TryParse(mPeptideDeltaMass, dblDelM) Then

    ''            ' Negate dblDelM so that it represents observed - theoretical
    ''            dblDelM = -dblDelM

    ''            ' Compute the original value for the precursor monoisotopic mass
    ''            If Double.TryParse(PrecursorMZ, dblPrecursorMZ) Then
    ''                If Integer.TryParse(mCharge, intCharge) Then
    ''                    dblPrecursorMonoMass = dblPrecursorMZ * intCharge - intCharge * clsPeptideMassCalculator.MASS_PROTON
    ''                Else
    ''                    blnParseError = True
    ''                End If
    ''            End If

    ''            If blnParseError Then
    ''                dblPrecursorMonoMass = mPeptideMonoisotopicMass + dblDelM
    ''            End If

    ''            MyBase.ComputeDelMCorrectedWork(dblDelM, dblPrecursorMonoMass, True)

    ''        Else
    ''            mPeptideDeltaMassCorrectedPpm = 0
    ''        End If

    ''    End If

    ''End Sub

End Class
