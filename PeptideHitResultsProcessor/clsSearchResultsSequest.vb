Option Strict On

' This class can be used to track the peptide details for a Sequest search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 7, 2006
'
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

    Protected mNumScans As String

    Protected mPeptideXCorr As String
    Protected mPeptideXCorrNext As String           ' Not stored in the synopsis or first hits file; it can be computed from XCorr and DeltaCn2
    Protected mPeptideDeltaCn As String
    Protected mPeptideSp As String
    Protected mPeptideDeltaCn2 As String
    Protected mPeptideRankSP As String
    Protected mPeptideRankXC As String
    Protected mPeptideXcRatio As String
    Protected mPeptidePassFilt As String
    Protected mPeptideMScore As String
    Protected mPeptideNTT As String
#End Region

#Region "Properties"
    Public Property NumScans() As String
        Get
            Return mNumScans
        End Get
        Set(ByVal Value As String)
            mNumScans = Value
        End Set
    End Property

    Public Property PeptideDeltaCn() As String
        Get
            Return mPeptideDeltaCn
        End Get
        Set(ByVal Value As String)
            mPeptideDeltaCn = Value
        End Set
    End Property
    Public Property PeptideDeltaCn2() As String
        Get
            Return mPeptideDeltaCn2
        End Get
        Set(ByVal Value As String)
            mPeptideDeltaCn2 = Value
            ComputePeptideXCorrNext()
        End Set
    End Property
    Public Property PeptideMScore() As String
        Get
            Return mPeptideMScore
        End Get
        Set(ByVal Value As String)
            mPeptideMScore = Value
        End Set
    End Property
    Public Property PeptideNTT() As String
        Get
            Return mPeptideNTT
        End Get
        Set(ByVal Value As String)
            mPeptideNTT = Value
        End Set
    End Property
    Public Property PeptidePassFilt() As String
        Get
            Return mPeptidePassFilt
        End Get
        Set(ByVal Value As String)
            mPeptidePassFilt = Value
        End Set
    End Property
    Public Property PeptideRankSP() As String
        Get
            Return mPeptideRankSP
        End Get
        Set(ByVal Value As String)
            mPeptideRankSP = Value
        End Set
    End Property
    Public Property PeptideRankXC() As String
        Get
            Return mPeptideRankXC
        End Get
        Set(ByVal Value As String)
            mPeptideRankXC = Value
        End Set
    End Property
    Public Property PeptideSp() As String
        Get
            Return mPeptideSp
        End Get
        Set(ByVal Value As String)
            mPeptideSp = Value
        End Set
    End Property
    Public Property PeptideXCorr() As String
        Get
            Return mPeptideXCorr
        End Get
        Set(ByVal Value As String)
            mPeptideXCorr = Value
        End Set
    End Property
    Public ReadOnly Property PeptideXCorrNext() As String
        Get
            Return mPeptideXCorrNext
        End Get
    End Property
    Public Property PeptideXcRatio() As String
        Get
            Return mPeptideXcRatio
        End Get
        Set(ByVal Value As String)
            mPeptideXcRatio = Value
        End Set
    End Property

#End Region

    Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        mNumScans = String.Empty

        mPeptideXCorr = String.Empty
        mPeptideXCorrNext = String.Empty
        mPeptideDeltaCn = String.Empty
        mPeptideSp = String.Empty
        mPeptideDeltaCn2 = String.Empty
        mPeptideRankSP = String.Empty
        mPeptideRankXC = String.Empty
        mPeptideXcRatio = String.Empty
        mPeptidePassFilt = String.Empty
        mPeptideMScore = String.Empty
        mPeptideNTT = String.Empty

    End Sub

    Protected Sub ComputePeptideXCorrNext()
        Try
            If clsPHRPBaseClass.IsNumber(mPeptideXCorr) And clsPHRPBaseClass.IsNumber(mPeptideDeltaCn2) Then
                mPeptideXCorrNext = (CSng(mPeptideXCorr) - CSng(mPeptideDeltaCn2) * CSng(mPeptideXCorr)).ToString
            Else
                mPeptideXCorrNext = "0"
            End If
        Catch ex As Exception
            mPeptideXCorrNext = "0"
        End Try
    End Sub

End Class
