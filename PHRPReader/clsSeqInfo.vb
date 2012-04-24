Option Strict On

Public Class clsSeqInfo
	Protected mSeqID As Integer
	Protected mModCount As Integer
	Protected mModDescription As String
	Protected mMonoisotopicMass As Double

	Public ReadOnly Property SeqID As Integer
		Get
			Return mSeqID
		End Get
	End Property

	Public ReadOnly Property ModCount As Integer
		Get
			Return mModCount
		End Get
	End Property

	Public ReadOnly Property ModDescription As String
		Get
			Return mModDescription
		End Get
	End Property

	Public ReadOnly Property MonoisotopicMass As Double
		Get
			Return mMonoisotopicMass
		End Get
	End Property

	Public Sub New(ByVal SeqID As Integer, ByVal MonoisotopicMass As Double)
		Me.New(SeqID, MonoisotopicMass, 0, String.Empty)
	End Sub

	Public Sub New(ByVal SeqID As Integer, ByVal MonoisotopicMass As Double, ByVal ModCount As Integer, ByVal ModDescription As String)
		mSeqID = SeqID
		mMonoisotopicMass = MonoisotopicMass
		mModCount = ModCount

		If String.IsNullOrEmpty(ModDescription) Then
			mModDescription = String.Empty
		Else
			mModDescription = ModDescription
		End If
	End Sub

	Public Sub UpdateMonoisotopicMass(ByVal MonoisotopicMass As Double)
		mMonoisotopicMass = MonoisotopicMass
	End Sub
End Class
