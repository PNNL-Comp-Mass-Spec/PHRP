Option Strict On

Imports PHRPReader.clsPeptideCleavageStateCalculator

Public Class clsProteinInfo
    Private ReadOnly mProteinName As String
    Private mProteinDescription As String
    Private ReadOnly mSeqID As Integer
    Private ReadOnly mCleavageState As ePeptideCleavageStateConstants
    Private ReadOnly mTerminusState As ePeptideTerminusStateConstants
    Private mResidueStart As Integer                        ' Residue number in the protein at which this sequence starts
    Private mResidueEnd As Integer                          ' Residue number in the protein at which this sequence ends

    Public ReadOnly Property ProteinName As String
        Get
            Return mProteinName
        End Get
    End Property

    ' ReSharper disable once ConvertToVbAutoProperty
    Public ReadOnly Property CleavageState As ePeptideCleavageStateConstants
        Get
            Return mCleavageState
        End Get
    End Property

    Public ReadOnly Property ResidueStart As Integer
        Get
            Return mResidueStart
        End Get
    End Property

    Public ReadOnly Property ResidueEnd As Integer
        Get
            Return mResidueEnd
        End Get
    End Property

    ' ReSharper disable once ConvertToVbAutoProperty
    Public ReadOnly Property SeqID As Integer
        Get
            Return mSeqID
        End Get
    End Property

    ' ReSharper disable once ConvertToVbAutoProperty
    Public ReadOnly Property TerminusState As ePeptideTerminusStateConstants
        Get
            Return mTerminusState
        End Get
    End Property

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="ProteinName"></param>
    ''' <param name="SeqID"></param>
    ''' <param name="CleavageState"></param>
    ''' <param name="TerminusState"></param>
    ''' <remarks></remarks>
    Public Sub New(
      ProteinName As String,
      SeqID As Integer,
      CleavageState As ePeptideCleavageStateConstants,
      TerminusState As ePeptideTerminusStateConstants)

        Me.New(ProteinName, String.Empty, SeqID, CleavageState, TerminusState)

    End Sub

    Public Sub New(
     ProteinName As String,
     ProteinDescription As String,
     SeqID As Integer,
     CleavageState As ePeptideCleavageStateConstants,
     TerminusState As ePeptideTerminusStateConstants)

        Me.New(ProteinName, String.Empty, SeqID, CleavageState, TerminusState, 0, 0)

    End Sub

    Public Sub New(
      ProteinName As String,
      ProteinDescription As String,
      SeqID As Integer,
      CleavageState As ePeptideCleavageStateConstants,
      TerminusState As ePeptideTerminusStateConstants,
      ProteinResidueStart As Integer,
      ProteinResidueEnd As Integer)

        mProteinName = ProteinName
        If String.IsNullOrEmpty(ProteinDescription) Then
            mProteinDescription = String.Empty
        Else
            mProteinDescription = ProteinDescription
        End If
        mSeqID = SeqID
        mCleavageState = CleavageState
        mTerminusState = TerminusState

        UpdateLocationInProtein(ProteinResidueStart, ProteinResidueEnd)

    End Sub

    Public Sub UpdateLocationInProtein(ProteinResidueStart As Integer, ProteinResidueEnd As Integer)

        mResidueStart = ProteinResidueStart
        mResidueEnd = ProteinResidueEnd

    End Sub

    Public Overloads Function ToString() As String
        Return mProteinName
    End Function

End Class
