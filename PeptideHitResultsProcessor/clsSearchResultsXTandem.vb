Option Strict On

' This class is used to track the peptide details for an XTandem search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 7, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

Imports PHRPReader

Public Class clsSearchResultsXTandem
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"
    ' Note: ProteinExpectationValue and ProteinIntensity are defined in clsSearchResultsBaseClass
    Protected mfI As String

    Protected mPeptideExpectationValue As String            ' The raw expectation value from the results file is converted to the Base-10 Log form when read into this program
    Protected mPeptideHyperscore As String
    Protected mPeptideNextScore As String
    Protected mPeptideDeltaCn2 As Single

    Protected mPeptideYScore As String
    Protected mPeptideYIons As String
    Protected mPeptideBScore As String
    Protected mPeptideBIons As String

    Protected mPeptideIntensity As String
    Protected mPeptideIntensityMax As String

    Protected mPeptideDeltaMassCorrectedPpm As Double

#End Region

#Region "Properties"
    Public Property fI() As String
        Get
            Return mfI
        End Get
        Set(Value As String)
            mfI = Value
        End Set
    End Property
    Public Property PeptideExpectationValue() As String
        Get
            Return mPeptideExpectationValue
        End Get
        Set(Value As String)
            mPeptideExpectationValue = Value
        End Set
    End Property
    Public Property PeptideHyperscore() As String
        Get
            Return mPeptideHyperscore
        End Get
        Set(Value As String)
            mPeptideHyperscore = Value
        End Set
    End Property
    Public Property PeptideNextScore() As String
        Get
            Return mPeptideNextScore
        End Get
        Set(Value As String)
            mPeptideNextScore = Value
            ComputePeptideDeltaCn2()
        End Set
    End Property
    Public ReadOnly Property PeptideDeltaCn2() As Single
        Get
            Return mPeptideDeltaCn2
        End Get
    End Property
    Public Property PeptideYScore() As String
        Get
            Return mPeptideYScore
        End Get
        Set(Value As String)
            mPeptideYScore = Value
        End Set
    End Property
    Public Property PeptideYIons() As String
        Get
            Return mPeptideYIons
        End Get
        Set(Value As String)
            mPeptideYIons = Value
        End Set
    End Property
    Public Property PeptideBScore() As String
        Get
            Return mPeptideBScore
        End Get
        Set(Value As String)
            mPeptideBScore = Value
        End Set
    End Property
    Public Property PeptideBIons() As String
        Get
            Return mPeptideBIons
        End Get
        Set(Value As String)
            mPeptideBIons = Value
        End Set
    End Property
    Public Property PeptideIntensity() As String
        Get
            Return mPeptideIntensity
        End Get
        Set(Value As String)
            mPeptideIntensity = Value
        End Set
    End Property
    Public Property PeptideIntensityMax() As String
        Get
            Return mPeptideIntensityMax
        End Get
        Set(Value As String)
            mPeptideIntensityMax = Value
        End Set
    End Property
    Public Property PeptideDeltaMassCorrectedPpm() As Double
        Get
            Return mPeptideDeltaMassCorrectedPpm
        End Get
        Set(value As Double)
            mPeptideDeltaMassCorrectedPpm = value
        End Set
    End Property
#End Region

    Public Sub New(objPeptideMods As clsPeptideModificationContainer, peptideSeqMassCalculator As clsPeptideMassCalculator)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods, peptideSeqMassCalculator)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        mfI = String.Empty

        mPeptideExpectationValue = String.Empty
        mPeptideHyperscore = String.Empty
        mPeptideNextScore = String.Empty
        mPeptideDeltaCn2 = 0

        mPeptideYScore = String.Empty
        mPeptideYIons = String.Empty
        mPeptideBScore = String.Empty
        mPeptideBIons = String.Empty

        mPeptideIntensity = String.Empty
        mPeptideIntensityMax = String.Empty

        mPeptideDeltaMassCorrectedPpm = 0
    End Sub

    Public Sub ComputeDelMCorrectedXT()

        Dim dblDelM As Double
        Dim intCorrectionCount = 0

        Dim dblPrecursorMonoMass As Double

        Dim blnParseError = False

        ' Note that mPeptideDeltaMass is the DeltaMass value reported by X!Tandem
        ' (though clsXtandemResultsProcessor took the negative of the value in the results file so it currently represents "theoretical - observed")
        If Double.TryParse(PeptideDeltaMass, dblDelM) Then

            ' Negate dblDelM so that it represents observed - theoretical
            dblDelM = -dblDelM

            ' Compute the original value for the precursor monoisotopic mass
            Dim dblParentIonMH As Double
            If Double.TryParse(MyBase.ParentIonMH, dblParentIonMH) Then
                dblPrecursorMonoMass = dblParentIonMH - PHRPReader.clsPeptideMassCalculator.MASS_PROTON
            Else
                blnParseError = True
            End If

            If blnParseError Then
                dblPrecursorMonoMass = PeptideMonoisotopicMass + dblDelM
            End If

            Const blnAdjustPrecursorMassForC13 = True
            mPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblDelM, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, PeptideMonoisotopicMass)

        Else
            mPeptideDeltaMassCorrectedPpm = 0
        End If

    End Sub

    Protected Sub ComputePeptideDeltaCn2()
        Try
            If clsPHRPParser.IsNumber(mPeptideHyperscore) And clsPHRPParser.IsNumber(mPeptideNextScore) Then
                mPeptideDeltaCn2 = (CSng(mPeptideHyperscore) - CSng(mPeptideNextScore)) / CSng(mPeptideHyperscore)
            Else
                mPeptideDeltaCn2 = 0
            End If
        Catch ex As Exception
            mPeptideDeltaCn2 = 0
        End Try
    End Sub

End Class
