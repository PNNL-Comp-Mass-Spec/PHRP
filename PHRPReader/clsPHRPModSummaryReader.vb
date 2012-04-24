Option Strict On

Public Class clsPHRPModSummaryReader

	Protected Const MOD_SUMMARY_COLUMN_Modification_Symbol As String = "Modification_Symbol"
	Protected Const MOD_SUMMARY_COLUMN_Modification_Mass As String = "Modification_Mass"
	Protected Const MOD_SUMMARY_COLUMN_Target_Residues As String = "Target_Residues"
	Protected Const MOD_SUMMARY_COLUMN_Modification_Type As String = "Modification_Type"
	Protected Const MOD_SUMMARY_COLUMN_Mass_Correction_Tag As String = "Mass_Correction_Tag"
	Protected Const MOD_SUMMARY_COLUMN_Occurence_Count As String = "Occurence_Count"

	Protected mModificationDefs As System.Collections.Generic.List(Of clsModificationDefinition)

	' The keys in this dictionary are MassCorrectionTag names and the values are the modification mass, stored as text (as it appears in the _ModSummary file)
	Protected mModDefMassesAsText As System.Collections.Generic.Dictionary(Of String, String)

	Protected mSuccess As Boolean

	Public ReadOnly Property ModificationDefs() As System.Collections.Generic.List(Of clsModificationDefinition)
		Get
			Return mModificationDefs
		End Get
	End Property

	Public ReadOnly Property Success() As Boolean
		Get
			Return mSuccess
		End Get
	End Property

	Public Sub New(ByVal strModSummaryFilePath As String)

		mModificationDefs = New System.Collections.Generic.List(Of clsModificationDefinition)
		mModDefMassesAsText = New System.Collections.Generic.Dictionary(Of String, String)

		mSuccess = False

		If String.IsNullOrEmpty(strModSummaryFilePath) Then
			Throw New Exception("ModSummaryFilePath is empty; unable to continue")
		ElseIf Not System.IO.File.Exists(strModSummaryFilePath) Then
			Throw New System.IO.FileNotFoundException("File not found: " & strModSummaryFilePath)
		End If

		mSuccess = ReadModSummaryFile(strModSummaryFilePath, mModificationDefs)

	End Sub

	''' <summary>
	''' Returns the mass value associated with the given mass correction tag
	''' </summary>
	''' <param name="strMassCorrectionTag"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Function GetModificationMassAsText(ByVal strMassCorrectionTag As String) As String
		Dim strModMass As String = String.Empty

		If mModDefMassesAsText.TryGetValue(strMassCorrectionTag, strModMass) Then
			Return strModMass
		Else
			Return String.Empty
		End If
	End Function

	Protected Function ReadModSummaryFile(ByVal strModSummaryFilePath As String, ByRef lstModInfo As System.Collections.Generic.List(Of clsModificationDefinition)) As Boolean

		Dim strLineIn As String
		Dim strSplitLine() As String

		Dim objColumnHeaders As SortedDictionary(Of String, Integer)

		Dim strModSymbol As String
		Dim strModMass As String
		Dim strTargetResidues As String
		Dim strModType As String
		Dim strMassCorrectionTag As String

		Dim chModSymbol As Char
		Dim dblModificationMass As Double
		Dim eModificationType As clsModificationDefinition.eModificationTypeConstants

		Dim lstModMasses As List(Of String) = Nothing

		Dim blnSkipLine As Boolean
		Dim blnHeaderLineParsed As Boolean


		If lstModInfo Is Nothing Then
			lstModInfo = New System.Collections.Generic.List(Of clsModificationDefinition)
		Else
			lstModInfo.Clear()
		End If

		If String.IsNullOrEmpty(strModSummaryFilePath) Then
			Return False
		End If

		' Initialize the column mapping
		' Using a case-insensitive comparer
		objColumnHeaders = New SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

		' Define the default column mapping
		objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Symbol, 0)
		objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Mass, 1)
		objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Target_Residues, 2)
		objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Type, 3)
		objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Mass_Correction_Tag, 4)
		objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Occurence_Count, 5)

		' Read the data from the ModSummary.txt file
		' The first line is typically a header line:
		' Modification_Symbol	Modification_Mass	Target_Residues	Modification_Type	Mass_Correction_Tag	Occurence_Count

		Using srModSummaryFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strModSummaryFilePath, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read))

			blnHeaderLineParsed = False

			Do While srModSummaryFile.Peek >= 0
				strLineIn = srModSummaryFile.ReadLine
				blnSkipLine = False

				If Not String.IsNullOrEmpty(strLineIn) Then
					strSplitLine = strLineIn.Split(ControlChars.Tab)

					If Not blnHeaderLineParsed Then
						If strSplitLine(0).ToLower() = MOD_SUMMARY_COLUMN_Modification_Symbol.ToLower() Then
							' Parse the header line to confirm the column ordering
							clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders)
							blnSkipLine = True
						End If

						blnHeaderLineParsed = True
					End If

					If Not blnSkipLine AndAlso strSplitLine.Length >= 4 Then
						strModSymbol = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Symbol, objColumnHeaders)
						strModMass = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Mass, objColumnHeaders)
						strTargetResidues = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Target_Residues, objColumnHeaders)
						strModType = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Type, objColumnHeaders)
						strMassCorrectionTag = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Mass_Correction_Tag, objColumnHeaders)

						If String.IsNullOrWhiteSpace(strModSymbol) Then
							strModSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
						End If
						chModSymbol = strModSymbol.Chars(0)

						If Not Double.TryParse(strModMass, dblModificationMass) Then
							Throw New Exception("Modification mass is not numeric for MassCorrectionTag: " & strMassCorrectionTag & ": " & strModMass)
						End If

						If String.IsNullOrWhiteSpace(strModType) Then
							eModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType
						Else
							eModificationType = clsModificationDefinition.ModificationSymbolToModificationType(strModType.Chars(0))
						End If

						Dim objModDef As clsModificationDefinition
						objModDef = New clsModificationDefinition(chModSymbol, dblModificationMass, strTargetResidues, eModificationType, strMassCorrectionTag)

						lstModInfo.Add(objModDef)

						If Not mModDefMassesAsText.ContainsKey(strMassCorrectionTag) Then
							mModDefMassesAsText.Add(strMassCorrectionTag, strModMass)
						End If

					End If
				End If

			Loop

		End Using

		Return True

	End Function


End Class
