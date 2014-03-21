Option Strict On

Public Class clsPepToProteinMapInfo

	Public Structure udtProteinMapInfoType
		Public Protein As String
		Public ResidueStart As Integer
		Public ResidueEnd As Integer
	End Structure

	''' <summary>
	''' List of proteins and residue start/end positions
	''' </summary>
	''' <remarks></remarks>
	Protected mProteinMapInfo As List(Of udtProteinMapInfoType)

	Public ReadOnly Property ProteinMapInfo As List(Of udtProteinMapInfoType)
		Get
			Return mProteinMapInfo
		End Get
	End Property

	Public Sub New(ByVal proteinName As String, ByVal residueStart As Integer, ByVal residueEnd As Integer)

		mProteinMapInfo = New List(Of udtProteinMapInfoType)

		AddProtein(proteinName, residueStart, residueEnd)
	End Sub

	Public Sub AddProtein(ByVal proteinName As String, ByVal residueStart As Integer, ByVal residueEnd As Integer)
		Dim udtProteinMapInfo = New udtProteinMapInfoType
		udtProteinMapInfo.Protein = String.Copy(proteinName)
		udtProteinMapInfo.ResidueStart = residueStart
		udtProteinMapInfo.ResidueEnd = residueEnd

		For i As Integer = 0 To mProteinMapInfo.Count - 1
			If String.Equals(mProteinMapInfo(i).Protein, proteinName, StringComparison.CurrentCultureIgnoreCase) Then
				' Protein mapping already exists; check residueStart
				If mProteinMapInfo(i).ResidueStart = residueStart Then
					If mProteinMapInfo(i).ResidueEnd <> residueEnd Then
						' Update this entry
						mProteinMapInfo(i) = udtProteinMapInfo
					End If
					Exit Sub
				End If
			End If
		Next

		mProteinMapInfo.Add(udtProteinMapInfo)
	End Sub
End Class
