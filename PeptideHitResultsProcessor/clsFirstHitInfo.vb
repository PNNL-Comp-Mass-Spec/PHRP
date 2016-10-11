Public Class clsFirstHitInfo

    Public ReadOnly Property CleanSequence As String

    Public Property ProteinName As String

    Public Property ProteinNumber As Integer

    Public Sub New(peptideCleanSequence As String)
        CleanSequence = peptideCleanSequence
    End Sub

    Public Overrides Function ToString() As String
        Return CleanSequence & ", " & ProteinName & ", " & ProteinNumber
    End Function
End Class
