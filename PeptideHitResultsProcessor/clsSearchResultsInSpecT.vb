Option Strict On

' This class is used to track the peptide details for an InSpecT search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
' Program started August 19, 2008
'
' E-mail: john.sandoval@pnnl.gov
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

Public Class clsSearchResultsInSpecT
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"

    Protected mSpectrumFile As String
    ' Scan: tracked by the base class
    ' Annotation: aka Peptide, which is tracked by the base class
    ' Protein: tracked by the base class
    Protected mMQScore As String
    Protected mLength As String
    Protected mTotalPRMScore As String
    Protected mMedianPRMScore As String
    Protected mFractionY As String
    Protected mFractionB As String
    Protected mIntensity As String
    Protected mNTT As String
    Protected mpvalue As String
    Protected mFScore As String
    Protected mDeltaScore As String
    Protected mDeltaScoreOther As String
    Protected mDeltaNormMQScore As String
    Protected mDeltaNormTotalPRMScore As String
    Protected mRankTotalPRMScore As String
    Protected mRankFScore As String
    Protected mRecordNumber As String
    Protected mDBFilePos As String
    Protected mSpecFilePos As String
    Protected mPrecursorMZ As String
    Protected mPrecursorError As String

#End Region

#Region "Properties"
    Public Property SpectrumFile() As String
        Get
            Return mSpectrumFile
        End Get
        Set(ByVal value As String)
            mSpectrumFile = Value
        End Set
    End Property
    Public Property MQScore() As String
        Get
            Return mMQScore
        End Get
        Set(ByVal value As String)
            mMQScore = Value
        End Set
    End Property
    Public Property Length() As String
        Get
            Return mLength
        End Get
        Set(ByVal value As String)
            mLength = Value
        End Set
    End Property
    Public Property TotalPRMScore() As String
        Get
            Return mTotalPRMScore
        End Get
        Set(ByVal value As String)
            mTotalPRMScore = Value
        End Set
    End Property
    Public Property MedianPRMScore() As String
        Get
            Return mMedianPRMScore
        End Get
        Set(ByVal value As String)
            mMedianPRMScore = Value
        End Set
    End Property
    Public Property FractionY() As String
        Get
            Return mFractionY
        End Get
        Set(ByVal value As String)
            mFractionY = Value
        End Set
    End Property
    Public Property FractionB() As String
        Get
            Return mFractionB
        End Get
        Set(ByVal value As String)
            mFractionB = Value
        End Set
    End Property
    Public Property Intensity() As String
        Get
            Return mIntensity
        End Get
        Set(ByVal value As String)
            mIntensity = Value
        End Set
    End Property
    Public Property NTT() As String
        Get
            Return mNTT
        End Get
        Set(ByVal value As String)
            mNTT = Value
        End Set
    End Property
    Public Property pValue() As String
        Get
            Return mpValue
        End Get
        Set(ByVal value As String)
            mpValue = Value
        End Set
    End Property
    Public Property FScore() As String
        Get
            Return mFScore
        End Get
        Set(ByVal value As String)
            mFScore = Value
        End Set
    End Property
    Public Property DeltaScore() As String
        Get
            Return mDeltaScore
        End Get
        Set(ByVal value As String)
            mDeltaScore = value
        End Set
    End Property
    Public Property DeltaScoreOther() As String
        Get
            Return mDeltaScoreOther
        End Get
        Set(ByVal value As String)
            mDeltaScoreOther = Value
        End Set
    End Property
    Public Property DeltaNormMQScore() As String
        Get
            Return mDeltaNormMQScore
        End Get
        Set(ByVal value As String)
            mDeltaNormMQScore = value
        End Set
    End Property
    Public Property DeltaNormTotalPRMScore() As String
        Get
            Return mDeltaNormTotalPRMScore
        End Get
        Set(ByVal value As String)
            mDeltaNormTotalPRMScore = value
        End Set
    End Property

    Public Property RankTotalPRMScore() As String
        Get
            Return mRankTotalPRMScore
        End Get
        Set(ByVal value As String)
            mRankTotalPRMScore = value
        End Set
    End Property
    Public Property RankFScore() As String
        Get
            Return mRankFScore
        End Get
        Set(ByVal value As String)
            mRankFScore = value
        End Set
    End Property
    Public Property RecordNumber() As String
        Get
            Return mRecordNumber
        End Get
        Set(ByVal value As String)
            mRecordNumber = Value
        End Set
    End Property
    Public Property DBFilePos() As String
        Get
            Return mDBFilePos
        End Get
        Set(ByVal value As String)
            mDBFilePos = Value
        End Set
    End Property
    Public Property SpecFilePos() As String
        Get
            Return mSpecFilePos
        End Get
        Set(ByVal value As String)
            mSpecFilePos = Value
        End Set
    End Property
    Public Property PrecursorMZ() As String
        Get
            Return mPrecursorMZ
        End Get
        Set(ByVal value As String)
            mPrecursorMZ = value
        End Set
    End Property
    Public Property PrecursorError() As String
        Get
            Return mPrecursorError
        End Get
        Set(ByVal value As String)
            mPrecursorError = value
        End Set
    End Property
#End Region

    Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        mSpectrumFile = String.Empty
        mMQScore = String.Empty
        mLength = String.Empty
        mTotalPRMScore = String.Empty
        mMedianPRMScore = String.Empty
        mFractionY = String.Empty
        mFractionB = String.Empty
        mIntensity = String.Empty
        mNTT = String.Empty
        mpValue = String.Empty
        mFScore = String.Empty
        mDeltaScore = String.Empty
        mDeltaScoreOther = String.Empty
        mDeltaNormMQScore = String.Empty
        mDeltaNormTotalPRMScore = String.Empty
        mRankTotalPRMScore = String.Empty
        mRankFScore = String.Empty
        mRecordNumber = String.Empty
        mDBFilePos = String.Empty
        mSpecFilePos = String.Empty
        mPrecursorMZ = String.Empty
        mPrecursorError = String.Empty
    End Sub

End Class
