' This class will compute the cleavage state and terminus state of a given peptide sequence.  
' It can also be used to remove modification symbols from a sequence using ExtractCleanSequenceFromSequenceWithMods

' The sequence can simply contain single-letter amino acid symbols (capital letters) or a mix 
'  of amino acid symbols and modification symbols, for example:
'   A.BCDEFGHIJK.L
'   A.B*CDEFGHIJK.L
'   A.BCDEFGHIJK*.L
'   A.BCDEFGHIJK.L

' Function ComputeCleavageState is overloaded to either except the peptide sequence with
'  prefix and suffix letters (e.g. A.BCDEFGHIJK.L) or accept the primary peptide sequence,
'  the prefix residue(s), and the suffix residue(s).

' Use EnzymeMatchSpec to specify the residues to match for cleavage

' The default cleavage specification is for trypsin: [KR]|[^P]

' Note: Function SplitPrefixAndSuffixFromSequence will change peptides that look like:
'      E.TGMLTQKFARSLGMLAVDNQARV..   to   E.TGMLTQKFARSLGMLAVDNQARV.
'   or ..TGMLTQKFARSLGMLAVDNQARV.R   to   .TGMLTQKFARSLGMLAVDNQARV.R
'

' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 4, 2006
'
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'

Option Strict On

Public Class clsPeptideCleavageStateCalculator

#Region "Constants and Enums"

    Public Const GENERIC_RESIDUE_SYMBOL As Char = "X"c
    Public Const TERMINUS_SYMBOL_SEQUEST As Char = "-"c
    Public Const TERMINUS_SYMBOL_XTANDEM_NTerminus As Char = "["c
    Public Const TERMINUS_SYMBOL_XTANDEM_CTerminus As Char = "]"c

    Protected Const TRYPSIN_LEFT_RESIDUE_REGEX As String = "[KR]"
    Protected Const TRYPSIN_RIGHT_RESIDUE_REGEX As String = "[^P]"

    Public Enum ePeptideCleavageStateConstants As Integer
        NonSpecific = 0                         ' e.g. Non-tryptic
		[Partial] = 1							' e.g. Partially tryptic
        Full = 2                                ' e.g. Fully tryptic
    End Enum

    Public Enum ePeptideTerminusStateConstants As Integer
        None = 0                        ' The peptide is located in the middle of the protein
        ProteinNTerminus = 1            ' The peptide is located at the protein's N-terminus
        ProteinCTerminus = 2            ' The peptide is located at the protein's C-terminus
        ProteinNandCCTerminus = 3       ' The peptide spans the entire length of the protein
    End Enum

    Public Enum eStandardCleavageAgentConstants As Integer
        Trypsin = 0
        TrypsinWithoutProlineRule = 1
        TrypsinPlusFVLEY = 2
        Chymotrypsin = 3
        ChymotrypsinAndTrypsin = 4
        V8_aka_GluC = 5
        CyanBr = 6
        EndoArgC = 7
        EndoLysC = 8
        EndoAspN = 9
        V8 = 10
    End Enum
#End Region

#Region "Structures"
    ' Example RegEx match strings for udtEnzymeMatchSpecType:
    ' [KR] means to match K or R
    ' [^P] means the residue cannot be P
    ' [A-Z] means to match anything; empty string also means match anything
    ' Note, this function will automatically change [X] to [A-Z] (provided GENERIC_RESIDUE_SYMBOL = "X")
    Public Structure udtEnzymeMatchSpecType
        Public LeftResidueRegEx As String           ' RegEx match string for matching the residue to the left of the cleavage point
        Public RightResidueRegEx As String          ' RegEx match string for matching the residue to the right of the cleavage point
    End Structure

#End Region

#Region "Classwide Variables"
    Private mEnzymeMatchSpec As udtEnzymeMatchSpecType
    Private mLeftRegEx As System.Text.RegularExpressions.Regex
    Private mRightRegEx As System.Text.RegularExpressions.Regex

    ' This array holds TERMINUS_SYMBOL_SEQUEST, TERMINUS_SYMBOL_XTANDEM_NTerminus, and TERMINUS_SYMBOL_XTANDEM_CTerminus
    '  and is useful for quickly checking for the presence of a terminus symbol using a binary search
    Private mTerminusSymbols() As Char
#End Region

#Region "Properties"
    Public Property EnzymeMatchSpec() As udtEnzymeMatchSpecType
        Get
            Return mEnzymeMatchSpec
        End Get
        Set(ByVal Value As udtEnzymeMatchSpecType)
            SetEnzymeMatchSpec(Value.LeftResidueRegEx, Value.RightResidueRegEx)
        End Set
    End Property

    Public ReadOnly Property TerminusSymbols() As Char()
        Get
            Return mTerminusSymbols
        End Get
    End Property
#End Region

	''' <summary>
	''' Constructor
	''' </summary>
	''' <remarks></remarks>
    Public Sub New()
        ReDim mTerminusSymbols(2)
        mTerminusSymbols(0) = TERMINUS_SYMBOL_SEQUEST
        mTerminusSymbols(1) = TERMINUS_SYMBOL_XTANDEM_NTerminus
        mTerminusSymbols(2) = TERMINUS_SYMBOL_XTANDEM_CTerminus
        Array.Sort(mTerminusSymbols)

        SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants.Trypsin)
    End Sub

	''' <summary>
	''' Converts Cleavage State to 0, 1, or 2
	''' </summary>
	''' <param name="eCleavageState"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function CleavageStateToShort(ByVal eCleavageState As ePeptideCleavageStateConstants) As Short
		Return CShort(eCleavageState)
	End Function

	''' <summary>
	''' Determines the cleavage state of the specified peptide
	''' </summary>
	''' <param name="strSequenceWithPrefixAndSuffix"></param>
	''' <returns></returns>
	''' <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
	Public Function ComputeCleavageState(ByVal strSequenceWithPrefixAndSuffix As String) As ePeptideCleavageStateConstants
		Dim strPrimarySequence As String = String.Empty
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty

		If SplitPrefixAndSuffixFromSequence(strSequenceWithPrefixAndSuffix, strPrimarySequence, strPrefix, strSuffix) Then
			Return ComputeCleavageState(strPrimarySequence, strPrefix, strSuffix)
		Else
			Return ePeptideCleavageStateConstants.NonSpecific
		End If
	End Function

	''' <summary>
	''' Determines the cleavage state of the specified peptide
	''' </summary>
	''' <param name="strCleanSequence"></param>
	''' <param name="strPrefixResidues"></param>
	''' <param name="strSuffixResidues"></param>
	''' <returns></returns>
	''' <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
	Public Function ComputeCleavageState(ByVal strCleanSequence As String, ByVal strPrefixResidues As String, ByVal strSuffixResidues As String) As ePeptideCleavageStateConstants
		' Determine the cleavage state of strCleanSequence utilizing the rules specified in mEnzymeMatchSpec

		Dim chSequenceStart As Char, chSequenceEnd As Char
		Dim chPrefix As Char
		Dim chSuffix As Char

		Dim ePeptideCleavageState As ePeptideCleavageStateConstants = ePeptideCleavageStateConstants.NonSpecific
		Dim ePeptideTerminusState As ePeptideTerminusStateConstants = ePeptideTerminusStateConstants.None
		Dim blnRuleMatchStart As Boolean
		Dim blnRuleMatchEnd As Boolean

		If Not strCleanSequence Is Nothing AndAlso strCleanSequence.Length > 0 Then
			' Find the letter closest to the end of strPrefixResidues
			chPrefix = FindLetterNearestEnd(strPrefixResidues)

			' Find the letter closest to the start of strSuffixResidues
			chSuffix = FindLetterNearestStart(strSuffixResidues)

			' Find the letter closest to the start of strCleanSequence
			chSequenceStart = FindLetterNearestStart(strCleanSequence)

			' Find the letter closest to the end of strCleanSequence
			chSequenceEnd = FindLetterNearestEnd(strCleanSequence)

			' Determine the terminus state of this peptide
			ePeptideTerminusState = ComputeTerminusState(chPrefix, chSuffix)

			If ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
				' The peptide spans the entire length of the protein; mark it as fully tryptic
				ePeptideCleavageState = ePeptideCleavageStateConstants.Full
			ElseIf ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus Then
				' Peptides at the N-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
				If TestCleavageRule(chSequenceEnd, chSuffix) Then
					ePeptideCleavageState = ePeptideCleavageStateConstants.Full
				Else
					' Leave ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
				End If
			ElseIf ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus Then
				' Peptides at the C-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic
				If TestCleavageRule(chPrefix, chSequenceStart) Then
					ePeptideCleavageState = ePeptideCleavageStateConstants.Full
				Else
					' Leave ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
				End If
			Else
				' Check whether chPrefix matches mLeftRegEx and chSequenceStart matches mRightRegEx
				blnRuleMatchStart = TestCleavageRule(chPrefix, chSequenceStart)
				blnRuleMatchEnd = TestCleavageRule(chSequenceEnd, chSuffix)

				If blnRuleMatchStart AndAlso blnRuleMatchEnd Then
					ePeptideCleavageState = ePeptideCleavageStateConstants.Full
				ElseIf blnRuleMatchStart OrElse blnRuleMatchEnd Then
					ePeptideCleavageState = ePeptideCleavageStateConstants.[Partial]
				Else
					' Leave ePeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
				End If
			End If
		End If

		Return ePeptideCleavageState

	End Function

	''' <summary>
	''' Count the number of missed cleavages in the peptide
	''' </summary>
	''' <param name="strSequenceWithPrefixAndSuffix"></param>
	''' <returns></returns>
	''' <remarks>Peptide can have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
	Public Function ComputeNumberOfMissedCleavages(ByVal strSequenceWithPrefixAndSuffix As String) As Short
		Dim strPrimarySequence As String = String.Empty
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty
		Dim strPreviousLetter As String

		Dim intNumMissedCleavages As Short = 0

		If SplitPrefixAndSuffixFromSequence(strSequenceWithPrefixAndSuffix, strPrimarySequence, strPrefix, strSuffix) Then
			If Not String.IsNullOrWhiteSpace(strPrimarySequence) Then

				strPreviousLetter = ""
				For intIndex As Integer = 0 To strPrimarySequence.Length - 1
					Dim chCurrent As Char = strPrimarySequence.Chars(intIndex)

					If clsPHRPReader.IsLetterAtoZ(chCurrent) Then

						If Not String.IsNullOrEmpty(strPreviousLetter) Then
							If TestCleavageRule(strPreviousLetter.Chars(0), chCurrent) Then
								intNumMissedCleavages += 1S
							End If
						End If

						strPreviousLetter = chCurrent
					End If
				Next
			End If
		End If

		Return intNumMissedCleavages

	End Function

	''' <summary>
	''' Determine the terminus state of the peptide
	''' </summary>
	''' <param name="strSequenceWithPrefixAndSuffix"></param>
	''' <returns></returns>
	''' <remarks>Peptide must have prefix and suffix letters, for example K.PEPTIDE.G</remarks>
	Public Function ComputeTerminusState(ByVal strSequenceWithPrefixAndSuffix As String) As ePeptideTerminusStateConstants
		Dim strPrimarySequence As String = String.Empty
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty

		If SplitPrefixAndSuffixFromSequence(strSequenceWithPrefixAndSuffix, strPrimarySequence, strPrefix, strSuffix) Then
			Return ComputeTerminusState(strPrimarySequence, strPrefix, strSuffix)
		Else
			Return ePeptideTerminusStateConstants.None
		End If
	End Function

	''' <summary>
	''' Determine the terminus state given the prefix and suffix characters
	''' </summary>
	''' <param name="chPrefix"></param>
	''' <param name="chSuffix"></param>
	''' <returns></returns>
	''' <remarks>For example, if the peptide is -.PEPTIDE.G then pass chPrefix="-" and chSuffix="G"</remarks>
	Public Function ComputeTerminusState(ByVal chPrefix As Char, ByVal chSuffix As Char) As ePeptideTerminusStateConstants

		Dim ePeptideTerminusState As ePeptideTerminusStateConstants = ePeptideTerminusStateConstants.None

		If Array.BinarySearch(mTerminusSymbols, chPrefix) >= 0 Then
			' Prefix character matches a terminus symbol
			If Array.BinarySearch(mTerminusSymbols, chSuffix) >= 0 Then
				' The peptide spans the entire length of the protein
				ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus
			Else
				' The peptide is located at the protein's N-terminus
				ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus
			End If
		ElseIf Array.BinarySearch(mTerminusSymbols, chSuffix) >= 0 Then
			' Suffix character matches a terminus symbol
			' The peptide is located at the protein's C-terminus
			ePeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus
		Else
			' Leave ePeptideTerminusState = ePeptideTerminusStateConstants.None
		End If

		Return ePeptideTerminusState

	End Function

	''' <summary>
	''' Determine the terminus state of the peptide
	''' </summary>
	''' <param name="strCleanSequence"></param>
	''' <param name="strPrefixResidues"></param>
	''' <param name="strSuffixResidues"></param>
	''' <returns></returns>
	''' <remarks>Peptide cannot have prefix and suffix letters, and thus must be in the form PEPTIDE</remarks>
	Public Function ComputeTerminusState(ByVal strCleanSequence As String, ByVal strPrefixResidues As String, ByVal strSuffixResidues As String) As ePeptideTerminusStateConstants
		' Determine the terminus state of strCleanSequence

		Dim chPrefix As Char
		Dim chSuffix As Char
		Dim ePeptideTerminusState As ePeptideTerminusStateConstants = ePeptideTerminusStateConstants.None

		If strCleanSequence Is Nothing OrElse strCleanSequence.Length = 0 Then
			ePeptideTerminusState = ePeptideTerminusStateConstants.None
		Else
			' Find the letter closest to the end of strPrefixResidues
			chPrefix = FindLetterNearestEnd(strPrefixResidues)

			' Find the letter closest to the start of strSuffixResidues
			chSuffix = FindLetterNearestStart(strSuffixResidues)

			ePeptideTerminusState = ComputeTerminusState(chPrefix, chSuffix)
		End If

		Return ePeptideTerminusState

	End Function

	''' <summary>
	''' Removes all modification symbols (*, #, +, 8, etc.) from the peptide; optionally removes prefix and suffix letters
	''' </summary>
	''' <param name="strSequenceWithMods"></param>
	''' <param name="blnCheckForPrefixAndSuffixResidues"></param>
	''' <returns>Clean peptide sequence</returns>
	''' <remarks></remarks>
	Public Shared Function ExtractCleanSequenceFromSequenceWithMods(ByVal strSequenceWithMods As String, ByVal blnCheckForPrefixAndSuffixResidues As Boolean) As String

		Dim strPrimarySequence As String = String.Empty
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty

		Dim chChar As Char
		Dim strCleanSequence As String

		strCleanSequence = String.Empty
		If Not strSequenceWithMods Is Nothing Then
			If blnCheckForPrefixAndSuffixResidues Then
				If Not SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then
					strPrimarySequence = String.Copy(strSequenceWithMods)
				End If
			Else
				strPrimarySequence = String.Copy(strSequenceWithMods)
			End If

			For Each chChar In strPrimarySequence
				If clsPHRPReader.IsLetterAtoZ(chChar) Then
					strCleanSequence &= chChar
				End If
			Next chChar
		End If

		Return strCleanSequence
	End Function

	Private Function FindLetterNearestEnd(ByVal strText As String) As Char
		Dim intIndex As Integer
		Dim chMatch As Char

		If strText Is Nothing OrElse strText.Length = 0 Then
			chMatch = TERMINUS_SYMBOL_SEQUEST
		Else
			intIndex = strText.Length - 1
			chMatch = strText.Chars(intIndex)
			Do While Not (clsPHRPReader.IsLetterAtoZ(chMatch) OrElse Array.BinarySearch(mTerminusSymbols, chMatch) >= 0) AndAlso intIndex > 0
				intIndex -= 1
				chMatch = strText.Chars(intIndex)
			Loop
		End If

		Return chMatch

	End Function

	Private Function FindLetterNearestStart(ByVal strText As String) As Char
		Dim intIndex As Integer
		Dim chMatch As Char

		If strText Is Nothing OrElse strText.Length = 0 Then
			chMatch = TERMINUS_SYMBOL_SEQUEST
		Else
			intIndex = 0
			chMatch = strText.Chars(intIndex)
			Do While Not (clsPHRPReader.IsLetterAtoZ(chMatch) OrElse Array.BinarySearch(mTerminusSymbols, chMatch) >= 0) AndAlso intIndex < strText.Length - 1
				intIndex += 1
				chMatch = strText.Chars(intIndex)
			Loop
		End If

		Return chMatch

	End Function

	''' <summary>
	''' Returns the default enzyme RegEx match specifications
	''' </summary>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Shared Function GetDefaultEnzymeMatchSpec() As udtEnzymeMatchSpecType
		Dim udtEnzymeMatchSpec As udtEnzymeMatchSpecType
		With udtEnzymeMatchSpec
			.LeftResidueRegEx = TRYPSIN_LEFT_RESIDUE_REGEX
			.RightResidueRegEx = TRYPSIN_RIGHT_RESIDUE_REGEX
		End With
		Return udtEnzymeMatchSpec
	End Function

	Private Sub InitializeRegExObjects()
		Const REGEX_OPTIONS As Text.RegularExpressions.RegexOptions = Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.Singleline Or Text.RegularExpressions.RegexOptions.IgnoreCase

		If mEnzymeMatchSpec.LeftResidueRegEx Is Nothing OrElse mEnzymeMatchSpec.RightResidueRegEx Is Nothing Then
			' Note that calling SetStandardEnzymeMatchSpec will cause this sub to also be called
			SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants.Trypsin)
			Exit Sub
		End If

		Try
			mLeftRegEx = New System.Text.RegularExpressions.Regex(mEnzymeMatchSpec.LeftResidueRegEx, REGEX_OPTIONS)
			mRightRegEx = New System.Text.RegularExpressions.Regex(mEnzymeMatchSpec.RightResidueRegEx, REGEX_OPTIONS)
		Catch ex As Exception
			' Ignore errors here
		End Try

	End Sub

	''' <summary>
	''' Define custom enzyme match rules using RegEx strings
	''' </summary>
	''' <param name="strLeftResidueRegEx"></param>
	''' <param name="strRightResidueRegEx"></param>
	''' <remarks></remarks>
	Public Sub SetEnzymeMatchSpec(ByVal strLeftResidueRegEx As String, ByVal strRightResidueRegEx As String)
		If Not strLeftResidueRegEx Is Nothing AndAlso Not strRightResidueRegEx Is Nothing Then
			If strLeftResidueRegEx.Length = 0 Then strLeftResidueRegEx = "[A-Z]"
			If strRightResidueRegEx.Length = 0 Then strRightResidueRegEx = "[A-Z]"

			If strLeftResidueRegEx = GENERIC_RESIDUE_SYMBOL OrElse strLeftResidueRegEx = "[" & GENERIC_RESIDUE_SYMBOL & "]" Then
				strLeftResidueRegEx = "[A-Z]"
			End If

			If strRightResidueRegEx = GENERIC_RESIDUE_SYMBOL OrElse strRightResidueRegEx = "[" & GENERIC_RESIDUE_SYMBOL & "]" Then
				strRightResidueRegEx = "[A-Z]"
			End If

			If strLeftResidueRegEx = "[^" & GENERIC_RESIDUE_SYMBOL & "]" Then
				strLeftResidueRegEx = "[^A-Z]"
			End If

			If strRightResidueRegEx = "[^" & GENERIC_RESIDUE_SYMBOL & "]" Then
				strRightResidueRegEx = "[^A-Z]"
			End If

			mEnzymeMatchSpec.LeftResidueRegEx = strLeftResidueRegEx
			mEnzymeMatchSpec.RightResidueRegEx = strRightResidueRegEx
		End If

		InitializeRegExObjects()
	End Sub

	''' <summary>
	''' Select a standard enzyme match rule
	''' </summary>
	''' <param name="eStandardCleavageAgent"></param>
	''' <remarks></remarks>
	Public Sub SetStandardEnzymeMatchSpec(ByVal eStandardCleavageAgent As eStandardCleavageAgentConstants)

		Select Case eStandardCleavageAgent
			Case eStandardCleavageAgentConstants.Trypsin
				SetEnzymeMatchSpec(TRYPSIN_LEFT_RESIDUE_REGEX, TRYPSIN_RIGHT_RESIDUE_REGEX)
			Case eStandardCleavageAgentConstants.TrypsinWithoutProlineRule
				SetEnzymeMatchSpec("[KR]", "[A-Z]")
			Case eStandardCleavageAgentConstants.TrypsinPlusFVLEY
				SetEnzymeMatchSpec("[KRFYVEL]", "[A-Z]")
			Case eStandardCleavageAgentConstants.Chymotrypsin
				SetEnzymeMatchSpec("[FWYL]", "[A-Z]")
			Case eStandardCleavageAgentConstants.ChymotrypsinAndTrypsin
				SetEnzymeMatchSpec("[FWYLKR]", "[A-Z]")

			Case eStandardCleavageAgentConstants.V8_aka_GluC
				SetEnzymeMatchSpec("[ED]", "[A-Z]")
			Case eStandardCleavageAgentConstants.CyanBr
				SetEnzymeMatchSpec("[M]", "[A-Z]")
			Case eStandardCleavageAgentConstants.EndoArgC
				SetEnzymeMatchSpec("[R]", "[A-Z]")
			Case eStandardCleavageAgentConstants.EndoLysC
				SetEnzymeMatchSpec("[K]", "[A-Z]")
			Case eStandardCleavageAgentConstants.EndoAspN
				SetEnzymeMatchSpec("[A-Z]", "[D]")
			Case Else
				' Unknown agent; leave unchanged
		End Select
	End Sub

	''' <summary>
	''' Examines strSequenceIn and splits apart into prefix, primary sequence, and suffix
	''' </summary>
	''' <param name="strSequenceIn">Peptide sequence to examine</param>
	''' <param name="strPrimarySequence">Primary sequence (output)</param>
	''' <param name="strPrefix">Prefix residue (output)</param>
	''' <param name="strSuffix">Suffix residue (output)</param>
	''' <returns> Returns True if success, False if prefix and suffix residues were not found</returns>
	''' <remarks>If more than one character is present before the first period or after the last period, then all characters are returned
	''' If the peptide starts with ".." then it is auto-changed to start with "."
	''' If the peptide ends with ".." then it is auto-changed to end with "."
	''' </remarks>
	Public Shared Function SplitPrefixAndSuffixFromSequence(ByVal strSequenceIn As String, ByRef strPrimarySequence As String, ByRef strPrefix As String, ByRef strSuffix As String) As Boolean

		Dim intPeriodLoc1 As Integer
		Dim intPeriodLoc2 As Integer
		Dim blnSuccess As Boolean

		strPrefix = String.Empty
		strSuffix = String.Empty
		strPrimarySequence = String.Empty

		blnSuccess = False

		If strSequenceIn Is Nothing OrElse strSequenceIn.Length = 0 Then
			Return False
		Else
			If strSequenceIn.StartsWith("..") AndAlso strSequenceIn.Length > 2 Then
				strSequenceIn = "." & strSequenceIn.Substring(2)
			End If

			If strSequenceIn.EndsWith("..") AndAlso strSequenceIn.Length > 2 Then
				strSequenceIn = strSequenceIn.Substring(0, strSequenceIn.Length - 2) & "."
			End If

			strPrimarySequence = String.Copy(strSequenceIn)

			' See if strSequenceIn contains two periods
			intPeriodLoc1 = strSequenceIn.IndexOf("."c)
			If intPeriodLoc1 >= 0 Then
				intPeriodLoc2 = strSequenceIn.LastIndexOf("."c)

				If intPeriodLoc2 > intPeriodLoc1 + 1 Then
					' Sequence contains two periods with letters between the periods, 
					' For example, A.BCDEFGHIJK.L or ABCD.BCDEFGHIJK.L
					' Extract out the text between the periods
					strPrimarySequence = strSequenceIn.Substring(intPeriodLoc1 + 1, intPeriodLoc2 - intPeriodLoc1 - 1)
					If intPeriodLoc1 > 0 Then
						strPrefix = strSequenceIn.Substring(0, intPeriodLoc1)
					End If
					strSuffix = strSequenceIn.Substring(intPeriodLoc2 + 1)

					blnSuccess = True
				ElseIf intPeriodLoc2 = intPeriodLoc1 + 1 Then
					' Peptide contains two periods in a row
					If intPeriodLoc1 <= 1 Then
						strPrimarySequence = String.Empty

						If intPeriodLoc1 > 0 Then
							strPrefix = strSequenceIn.Substring(0, intPeriodLoc1)
						End If
						strSuffix = strSequenceIn.Substring(intPeriodLoc2 + 1)

						blnSuccess = True
					Else
						' Leave the sequence unchanged
						strPrimarySequence = String.Copy(strSequenceIn)
						blnSuccess = False
					End If
				ElseIf intPeriodLoc1 = intPeriodLoc2 Then
					' Peptide only contains one period
					If intPeriodLoc1 = 0 Then
						strPrimarySequence = strSequenceIn.Substring(1)
						blnSuccess = True
					ElseIf intPeriodLoc1 = strSequenceIn.Length - 1 Then
						strPrimarySequence = strSequenceIn.Substring(0, intPeriodLoc1)
						blnSuccess = True
					ElseIf intPeriodLoc1 = 1 AndAlso strSequenceIn.Length > 2 Then
						strPrimarySequence = strSequenceIn.Substring(intPeriodLoc1 + 1)
						strPrefix = strSequenceIn.Substring(0, intPeriodLoc1)
						blnSuccess = True
					ElseIf intPeriodLoc1 = strSequenceIn.Length - 2 Then
						strPrimarySequence = strSequenceIn.Substring(0, intPeriodLoc1)
						strSuffix = strSequenceIn.Substring(intPeriodLoc1 + 1)
						blnSuccess = True
					Else
						' Leave the sequence unchanged
						strPrimarySequence = String.Copy(strSequenceIn)
					End If
				End If
			End If
		End If

		Return blnSuccess
	End Function

	''' <summary>
	''' Examines the two residues to see if they represent an expected cleavage point
	''' </summary>
	''' <param name="chLeftChar"></param>
	''' <param name="chRightChar"></param>
	''' <returns>True if the characters match the currently defined cleavage rule</returns>
	''' <remarks></remarks>
	Public Function TestCleavageRule(ByVal chLeftChar As Char, ByVal chRightChar As Char) As Boolean

		With mLeftRegEx.Match(chLeftChar)
			If .Success Then
				With mRightRegEx.Match(chRightChar)
					If .Success Then
						Return True
					End If
				End With
			End If
		End With

		Return False

	End Function

End Class
