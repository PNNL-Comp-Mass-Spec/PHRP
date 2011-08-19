Option Strict On

' This class can be used to track the peptide details for a Sequest search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 7, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnnl.gov/ or http://www.sysbio.org/resources/staff/
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

Public Class clsSearchResultsSequest
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"

    Protected mPeptideXCorrNext As String           ' Not stored in the synopsis or first hits file; it can be computed from XCorr and DeltaCn2
#End Region

#Region "Properties"
    ' Auto-properties (only work with Visual Studio 2010)
    Public Property NumScans() As String
    Public Property PeptideDeltaCn() As String
    Public Property PeptideDeltaCn2() As String
    Public Property PeptideMScore() As String
    Public Property PeptideNTT() As String
    Public Property PeptidePassFilt() As String
    Public Property PeptideRankSP() As String
    Public Property PeptideRankXC() As String
    Public Property PeptideSp() As String
    Public Property PeptideXCorr() As String
    Public Property PeptideXcRatio() As String
    Public Property IonsObserved() As String
    Public Property IonsExpected() As String
    Public Property DelMPPM() As String

    Public ReadOnly Property PeptideXCorrNext() As String
        Get
            Return mPeptideXCorrNext
        End Get
    End Property

#End Region

    Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        NumScans = String.Empty
        PeptideDeltaCn = String.Empty
        PeptideDeltaCn2 = String.Empty
        PeptideMScore = String.Empty
        PeptideNTT = String.Empty
        PeptidePassFilt = String.Empty
        PeptideRankSP = String.Empty
        PeptideRankXC = String.Empty
        PeptideSp = String.Empty
        PeptideXCorr = String.Empty
        PeptideXcRatio = String.Empty
        IonsObserved = String.Empty
        IonsExpected = String.Empty
        DelMPPM = String.Empty

        mPeptideXCorrNext = String.Empty

    End Sub

    Protected Sub ComputePeptideXCorrNext()
        Dim sngXCorr As Single
        Dim sngDelCN2 As Single

        Try
            If Single.TryParse(PeptideXCorr, sngXCorr) AndAlso _
               Single.TryParse(PeptideDeltaCn2, sngDelCN2) Then
                mPeptideXCorrNext = (sngXCorr - sngDelCN2 * sngXCorr).ToString
            Else
                mPeptideXCorrNext = "0"
            End If

        Catch ex As Exception
            mPeptideXCorrNext = "0"
        End Try
    End Sub

    '' Unused
    ''Protected Overrides Sub ComputeDelMCorrected()

    ''    Dim dblDelM As Double
    ''    Dim intCorrectionCount As Integer = 0

    ''    Dim dblPrecursorMonoMass As Double

    ''    ' Note that mPeptideDeltaMass is the DeltaMass value reported by Sequest
    ''    ' (it currently represents "theoretical - observed")
    ''    If Double.TryParse(mPeptideDeltaMass, dblDelM) Then

    ''        ' Negate dblDelM so that it represents observed - theoretical
    ''        dblDelM = -dblDelM

    ''        ' Compute the original value for the precursor monoisotopic mass
    ''        dblPrecursorMonoMass = mPeptideMonoisotopicMass + dblDelM

    ''        MyBase.ComputeDelMCorrectedWork(dblDelM, dblPrecursorMonoMass, True)

    ''    Else
    ''        mPeptideDeltaMassCorrectedPpm = 0
    ''    End If

    ''End Sub

End Class
