Public Class clsEmpiricalFormula

    ''' <summary>
    ''' Elements in the empirical formula
    ''' Keys are element symbols, values are element counts
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property ElementCounts As Dictionary(Of String, Integer)

    ''' <summary>
    ''' Constructor
    ''' </summary>
    Public Sub New()
        ElementCounts() = New Dictionary(Of String, Integer)
    End Sub

    ''' <summary>
    ''' Constructor, initialized with an existing dictionary of element symbols and counts
    ''' </summary>
    Public Sub New(elementInfo As Dictionary(Of String, Integer))
        ElementCounts() = elementInfo
    End Sub

    ''' <summary>
    ''' Constructor, initialized with a list of element symbols
    ''' </summary>
    Public Sub New(elementInfo As IEnumerable(Of String))
        ElementCounts() = New Dictionary(Of String, Integer)
        For Each element In elementInfo
            AddElement(element, 1)
        Next
    End Sub

    ''' <summary>
    ''' Constructor, initialized with a list of KeyValuePairs of element symbol and element count
    ''' </summary>
    Public Sub New(elementInfo As IEnumerable(Of KeyValuePair(Of String, Integer)))
        ElementCounts() = New Dictionary(Of String, Integer)
        For Each element In elementInfo
            AddElement(element.Key, element.Value)
        Next
    End Sub

    ''' <summary>
    ''' Add a new element to the empirical formula
    ''' </summary>
    ''' <param name="elementSymbol"></param>
    ''' <param name="elementCount"></param>
    Public Sub AddElement(elementSymbol As String, elementCount As Integer)
        Dim existingCount As Integer
        If ElementCounts.TryGetValue(elementSymbol, existingCount) Then
            ElementCounts(elementSymbol) = existingCount + elementCount
        Else
            ElementCounts.Add(elementSymbol, elementCount)
        End If
    End Sub

    ''' <summary>
    ''' Adds all of the elements from the given empirical formula
    ''' </summary>
    ''' <param name="empiricalFormula"></param>
    Public Sub AddElements(empiricalFormula As clsEmpiricalFormula)
        For Each element In empiricalFormula.ElementCounts
            AddElement(element.Key, element.Value)
        Next
    End Sub

    ''' <summary>
    ''' Return the number of atoms of the given element in the empirical formula
    ''' </summary>
    ''' <param name="elementSymbol"></param>
    ''' <returns>Element Count, or 0 if the element is not in ElementCounts</returns>
    Public Function GetElementCount(elementSymbol As Char) As Integer
        Dim elementCount As Integer
        If ElementCounts.TryGetValue(elementSymbol, elementCount) Then
            Return elementCount
        End If

        Return 0
    End Function

    Public Overrides Function ToString() As String
        If ElementCounts.Count = 0 Then
            Return "<Undefined>"
        End If

        Dim formulaDescription As String = String.Empty
        For Each element In ElementCounts
            If element.Value = 1 Then
                formulaDescription &= element.Key
            Else
                formulaDescription &= element.Key & element.Value
            End If

        Next

        Return formulaDescription

    End Function

End Class
