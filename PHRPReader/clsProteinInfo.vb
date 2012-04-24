Option Strict On

Imports PHRPReader.clsPeptideCleavageStateCalculator

Public Class clsProteinInfo
	Protected mProteinName As String
	Protected mProteinDescription As String
	Protected mSeqID As Integer
	Protected mCleavageState As ePeptideCleavageStateConstants
	Protected mTerminusState As ePeptideTerminusStateConstants

	Public ReadOnly Property ProteinName As String
		Get
			Return mProteinName
		End Get
	End Property

	Public ReadOnly Property CleavageState As ePeptideCleavageStateConstants
		Get
			Return mCleavageState
		End Get
	End Property

	Public ReadOnly Property SeqID As Integer
		Get
			Return mSeqID
		End Get
	End Property

	Public ReadOnly Property TerminusState As ePeptideTerminusStateConstants
		Get
			Return mTerminusState
		End Get
	End Property

	Public Sub New(ByVal ProteinName As String, ByVal SeqID As Integer, ByVal CleavageState As ePeptideCleavageStateConstants, ByVal TerminusState As ePeptideTerminusStateConstants)
		Me.New(ProteinName, String.Empty, SeqID, CleavageState, TerminusState)
	End Sub

	Public Sub New(ByVal ProteinName As String, ByVal ProteinDescription As String, ByVal SeqID As Integer, ByVal CleavageState As ePeptideCleavageStateConstants, ByVal TerminusState As ePeptideTerminusStateConstants)
		mProteinName = ProteinName
		If String.IsNullOrEmpty(ProteinDescription) Then
			mProteinDescription = String.Empty
		Else
			mProteinDescription = ProteinDescription
		End If
		mSeqID = SeqID
		mCleavageState = CleavageState
		mTerminusState = TerminusState
	End Sub

End Class
