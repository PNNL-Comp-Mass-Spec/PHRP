Imports PHRPReader

Public Class clsFirstHitInfo

    ' ReSharper disable once UnassignedReadonlyField
    Private ReadOnly mPrimarySequence As String

    Private mPrefix As String
    Private mSuffix As String

    Public ReadOnly Property CleanSequence As String

    Public ReadOnly Property SequenceWithModsAndContext As String
        Get
            Return mPrefix & "." & mPrimarySequence & "." & mSuffix
        End Get
    End Property

    Public Property ProteinName As String

    Public Property ProteinNumber As Integer

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="peptideSeqWithModsAndContext"></param>
    ''' <param name="peptideCleanSeq"></param>
    Public Sub New(peptideSeqWithModsAndContext As String, peptideCleanSeq As String)

        If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptideSeqWithModsAndContext, mPrimarySequence, mPrefix, mSuffix) Then
            Throw New Exception("Unable to split the prefix and suffix from peptide " & peptideSeqWithModsAndContext)
        End If

        CleanSequence = peptideCleanSeq
    End Sub

    Public Overrides Function ToString() As String
        Return CleanSequence & ", " & ProteinName & ", " & ProteinNumber
    End Function

    Public Sub UpdatePrefixAndSuffix(newPrefix As String, newSuffix As String)
        mPrefix = newPrefix
        mSuffix = newSuffix
    End Sub
End Class
