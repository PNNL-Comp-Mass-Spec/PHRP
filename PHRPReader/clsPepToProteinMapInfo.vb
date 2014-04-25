Option Strict On

Public Class clsPepToProteinMapInfo

    Public Structure udtProteinLocationInfo
        Public ResidueStart As Integer
        Public ResidueEnd As Integer
    End Structure

    ''' <summary>
    ''' Dictionary of protein names and residue start/end positions
    ''' </summary>
    ''' <remarks></remarks>
    Protected mProteinMapInfo As Dictionary(Of String, List(Of udtProteinLocationInfo))

    Public ReadOnly Property ProteinMapInfo As Dictionary(Of String, List(Of udtProteinLocationInfo))
        Get
            Return mProteinMapInfo
        End Get
    End Property

    Public Sub New(ByVal proteinName As String, ByVal residueStart As Integer, ByVal residueEnd As Integer)

        mProteinMapInfo = New Dictionary(Of String, List(Of udtProteinLocationInfo))(StringComparer.CurrentCultureIgnoreCase)

        AddProtein(proteinName, residueStart, residueEnd)
    End Sub

    Public Sub AddProtein(ByVal proteinName As String, ByVal residueStart As Integer, ByVal residueEnd As Integer)

        Dim lstLocations As List(Of udtProteinLocationInfo) = Nothing

        If mProteinMapInfo.TryGetValue(proteinName, lstLocations) Then
            ' Protein mapping already exists; check residueStart
            For Each udtLoc In lstLocations
                If udtLoc.ResidueStart = residueStart Then
                    ' Update this entry
                    If udtLoc.ResidueEnd <> residueEnd Then
                        udtLoc.ResidueEnd = residueEnd
                    End If
                    Exit Sub
                End If
            Next

            Dim udtLocInfoAddnl = New udtProteinLocationInfo
            udtLocInfoAddnl.ResidueStart = residueStart
            udtLocInfoAddnl.ResidueEnd = residueEnd

            lstLocations.Add(udtLocInfoAddnl)
            Exit Sub
        End If

        Dim udtLocInfo = New udtProteinLocationInfo
        udtLocInfo.ResidueStart = residueStart
        udtLocInfo.ResidueEnd = residueEnd

        mProteinMapInfo.Add(proteinName, New List(Of udtProteinLocationInfo) From {udtLocInfo})

    End Sub

End Class
