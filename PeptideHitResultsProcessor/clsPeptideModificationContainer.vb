Option Strict On

' This class can be used to track modifications that can be applied to peptides
' It handles both residue level modifications and static, peptide-wide modifications
'
' Use ReadMassCorrectionTagsFile() and ReadModificationDefinitionsFile() to customize
'  the default mass correction tag and modification definition lists
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 5, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'
' Notice: This computer software was prepared by Battelle Memorial Institute, 
' hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
' Department of Energy (DOE).  All rights in the computer software are reserved 
' by DOE on behalf of the United States Government and the Contractor as 
' provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
' SOFTWARE.  This notice including this sentence must appear on any copies of 
' this computer software.

Imports PHRPReader
Imports PHRPReader.clsAminoAcidModInfo

Public Class clsPeptideModificationContainer

#Region "Constants and Enums"
    Public Const DEFAULT_MODIFICATION_SYMBOLS As String = "*#@$&!%~†‡¤º^`×÷+=ø¢"         ' A few other possibilities: €£¥§

    Public Const MASS_DIGITS_OF_PRECISION As Byte = 3

    Public Const N_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM As Char = "["c
    Public Const C_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM As Char = "]"c

    Public Const N_TERMINAL_PEPTIDE_MOD_SYMBOL_INSPECT As Char = "["c
    Public Const C_TERMINAL_PEPTIDE_MOD_SYMBOL_INSPECT As Char = "]"c

#End Region

#Region "Structures"
#End Region

#Region "Classwide Variables"
    ' List of available modification symbols
    Protected mDefaultModificationSymbols As System.Collections.Queue

    ' List of known mass correction tags
    Protected mMassCorrectionTags As Hashtable

    ' List of known modifications
	Protected mModifications As System.Collections.Generic.List(Of clsModificationDefinition)

    Protected mErrorMessage As String

    ' This array holds modifications that Sequest or XTandem will often use but for 
    ' which the auto-addition method sometimes incorrectly notes
    Protected mStandardRefinementModifications() As clsModificationDefinition

    Protected mConsiderModSymbolWhenFindingIdenticalMods As Boolean
#End Region

#Region "Properties"
    Public ReadOnly Property ErrorMessage() As String
        Get
            Return mErrorMessage
        End Get
    End Property

    Public ReadOnly Property ModificationCount() As Integer
        Get
			Return mModifications.Count
		End Get
	End Property

	Public ReadOnly Property Modifications As System.Collections.Generic.List(Of clsModificationDefinition)
		Get
			Return mModifications
		End Get
	End Property

	Public Property ConsiderModSymbolWhenFindingIdenticalMods() As Boolean
		Get
			Return mConsiderModSymbolWhenFindingIdenticalMods
		End Get
		Set(ByVal value As Boolean)
			mConsiderModSymbolWhenFindingIdenticalMods = value
		End Set
	End Property
#End Region

	Public Sub New()
		InitializeLocalVariables()
	End Sub

	''' <summary>
	''' Add objModificationDefinition to mModifications
	''' However, do not add if a duplicate modification
	''' Furthermore, if everything matches except for .TargetResidues, then add the new target residues to the existing, matching mod
	''' </summary>
	''' <param name="objModificationDefinition"></param>
	''' <param name="blnUseNextAvailableModificationSymbol"></param>
	''' <returns>The index of the newly added modification, or the the index of the modification that objModificationDefinition matches </returns>
	''' <remarks></remarks>
	Private Function AddModification(ByVal objModificationDefinition As clsModificationDefinition, ByVal blnUseNextAvailableModificationSymbol As Boolean) As Integer

		Dim chChar As Char
		Dim intModificationIndex As Integer
		Dim blnMatchFound As Boolean

		blnMatchFound = False

		' See if any of the existing modifications match objModificationDefinition, ignoring .TargetResidues and possibly ignoring .ModificationSymbol
		For intModificationIndex = 0 To mModifications.Count - 1
			If mModifications(intModificationIndex).EquivalentMassTypeTagAndAtom(objModificationDefinition) Then

				blnMatchFound = True

				If mConsiderModSymbolWhenFindingIdenticalMods Then
					If objModificationDefinition.ModificationSymbol <> mModifications(intModificationIndex).ModificationSymbol Then
						' Symbols differ; add this as a new modification definition
						blnMatchFound = False
					End If
				End If

				If blnMatchFound Then
					With mModifications(intModificationIndex)
						If .ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
						   .ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
							' Matching dynamic or static modification definitions
							' Merge the two modifications by making sure each of the residues in objModificationDefinition.TargetResidues is present in .TargetResidues
							For Each chChar In objModificationDefinition.TargetResidues
								If Not .TargetResiduesContain(chChar) Then
									.TargetResidues &= chChar
								End If
							Next chChar

							If Not blnUseNextAvailableModificationSymbol Then
								' See if the new modification symbol is different than the already-defined symbol
								If objModificationDefinition.ModificationSymbol <> .ModificationSymbol Then

								End If
							End If

						End If
					End With
				End If

			End If
			If blnMatchFound Then Exit For
		Next intModificationIndex

		If Not blnMatchFound Then
			If blnUseNextAvailableModificationSymbol AndAlso mDefaultModificationSymbols.Count > 0 Then
				' Add objModificationDefinition to the list, using the next available default modification symbol
				objModificationDefinition.ModificationSymbol = CChar(mDefaultModificationSymbols.Dequeue)
			Else
				' Leave .ModificationSymbol as-is
			End If

			mModifications.Add(objModificationDefinition)

			Return mModifications.Count - 1
		Else
			Return intModificationIndex
		End If

	End Function

	Public Sub AppendStandardRefinmentModifications()
		Dim intIndex As Integer

		For intIndex = 0 To mStandardRefinementModifications.Length - 1
			With mStandardRefinementModifications(intIndex)
				VerifyModificationPresent(.ModificationMass, .TargetResidues, .ModificationType)
			End With
		Next intIndex
	End Sub

	Public Sub ClearModifications()
		UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS)
		mModifications.Clear()
	End Sub

	Public Function GetModificationByIndex(ByVal intIndex As Integer) As clsModificationDefinition
		If intIndex >= 0 And intIndex < mModifications.Count Then
			Return mModifications(intIndex)
		Else
			Return New clsModificationDefinition
		End If
	End Function

	Public Function GetModificationTypeByIndex(ByVal intIndex As Integer) As clsModificationDefinition.eModificationTypeConstants
		If intIndex >= 0 And intIndex < mModifications.Count Then
			Return mModifications(intIndex).ModificationType
		Else
			Return clsModificationDefinition.eModificationTypeConstants.UnknownType
		End If
	End Function

	Public Function LookupMassCorrectionTagByMass(ByVal dblModificationMass As Double, Optional ByVal MassDigitsOfPrecision As Byte = MASS_DIGITS_OF_PRECISION, Optional ByVal blnAddToModificationListIfUnknown As Boolean = True) As String

		Dim objEnum As System.Collections.IDictionaryEnumerator

		Dim intMassDigitsOfPrecisionCurrent As Integer
		Dim intMassDigitsOfPrecisionStop As Integer

		Dim strMassCorrectionTag As String
		Dim dblMassDiff As Double

		Dim strClosestMassCorrectionTag As String
		Dim dblClosestMassCorrectionTagMassDiff As Double

		Dim intUnknownModValue As Integer
		Dim intLargestUnknownModValue As Integer

		If MassDigitsOfPrecision < 0 Then MassDigitsOfPrecision = 0

		If MassDigitsOfPrecision >= 1 Then
			intMassDigitsOfPrecisionStop = 1
		Else
			intMassDigitsOfPrecisionStop = MassDigitsOfPrecision
		End If

		For intMassDigitsOfPrecisionCurrent = MassDigitsOfPrecision To intMassDigitsOfPrecisionStop Step -1
			strClosestMassCorrectionTag = String.Empty
			dblClosestMassCorrectionTagMassDiff = Double.MaxValue
			intLargestUnknownModValue = 0

			Try
				' First look for an exact match in mMassCorrectionTags
				' At the same time, look for entries starting with clsModificationDefinition.UNKNOWN_MOD_BASE_NAME
				objEnum = mMassCorrectionTags.GetEnumerator
				Do While objEnum.MoveNext
					strMassCorrectionTag = CStr(objEnum.Key)
					dblMassDiff = Math.Abs(dblModificationMass - CDbl(objEnum.Value))
					If dblMassDiff < dblClosestMassCorrectionTagMassDiff Then
						strClosestMassCorrectionTag = CStr(objEnum.Key)
						dblClosestMassCorrectionTagMassDiff = dblMassDiff
					End If

					If strMassCorrectionTag.StartsWith(clsModificationDefinition.UNKNOWN_MOD_BASE_NAME) Then
						Try
							intUnknownModValue = CInt(strMassCorrectionTag.Substring(clsModificationDefinition.UNKNOWN_MOD_BASE_NAME.Length))
						Catch ex As Exception
							intUnknownModValue = 0
						End Try
						If intUnknownModValue > intLargestUnknownModValue Then
							intLargestUnknownModValue = intUnknownModValue
						End If
					End If
				Loop
			Catch ex As Exception
				' Error enumerating through mMassCorrectionTags
			End Try

			If Math.Round(dblClosestMassCorrectionTagMassDiff, intMassDigitsOfPrecisionCurrent) = 0 Then
				' Match found
				Return strClosestMassCorrectionTag
			Else
				If intMassDigitsOfPrecisionCurrent > intMassDigitsOfPrecisionStop Then
					' Let the For loop go through another iteration to see if we find a match
				Else
					If intLargestUnknownModValue < 99 Then
						strClosestMassCorrectionTag = clsModificationDefinition.UNKNOWN_MOD_BASE_NAME & (intLargestUnknownModValue + 1).ToString("00")
						If blnAddToModificationListIfUnknown Then
							Try
								mMassCorrectionTags.Add(strClosestMassCorrectionTag, dblModificationMass)
							Catch ex As Exception
								' This shouldn't happen, due to searching about for mass correction tags starting with clsModificationDefinition.UNKNOWN_MOD_BASE_NAME
								' Ignore the error
							End Try
						End If
					Else
						strClosestMassCorrectionTag = clsModificationDefinition.UNKNOWN_MOD_BASE_NAME & "99"
					End If

					Return strClosestMassCorrectionTag
				End If
			End If
		Next intMassDigitsOfPrecisionCurrent

		Return String.Empty

	End Function

	Public Function LookupDynamicModificationDefinitionByTargetInfo(ByVal chModificationSymbol As Char, ByVal chTargetResidue As Char, ByVal eResidueTerminusState As eResidueTerminusStateConstants, ByRef blnExistingModFound As Boolean) As clsModificationDefinition
		' Looks for a modification of type .DynamicMod or type .UnknownType in mModifications having .ModificationSymbol = chModificationSymbol and chTargetResidue in .TargetResidues
		' Note: If chModificationSymbol does not match any of the mods, then a modification with a mass of 0 is returned

		Dim intIndex As Integer

		Dim objModificationDefinition As clsModificationDefinition

		blnExistingModFound = False
		If Not chTargetResidue = Nothing OrElse eResidueTerminusState <> eResidueTerminusStateConstants.None Then
			' The residue was provided and/or the residue is located at a peptide or protein terminus
			' First compare against modifications with 1 or more residues in .TargetResidues
			For intIndex = 0 To mModifications.Count - 1
				If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
				   mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso _
				   mModifications(intIndex).TargetResidues.Length > 0 Then
					If mModifications(intIndex).ModificationSymbol = chModificationSymbol Then
						' Matching modification symbol found
						' Now see if .TargetResidues contains chTargetResidue
						If Not chTargetResidue = Nothing AndAlso mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
							blnExistingModFound = True
						End If

						If Not blnExistingModFound AndAlso eResidueTerminusState <> eResidueTerminusStateConstants.None Then

							Select Case eResidueTerminusState
								Case eResidueTerminusStateConstants.ProteinNTerminus, eResidueTerminusStateConstants.ProteinNandCCTerminus
									If mModifications(intIndex).TargetResiduesContain(N_TERMINAL_PROTEIN_SYMBOL_DMS) Then
										blnExistingModFound = True
									End If
								Case eResidueTerminusStateConstants.PeptideNTerminus
									If mModifications(intIndex).TargetResiduesContain(N_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
										blnExistingModFound = True
									End If
							End Select

							If Not blnExistingModFound Then
								Select Case eResidueTerminusState
									Case eResidueTerminusStateConstants.ProteinCTerminus, eResidueTerminusStateConstants.ProteinNandCCTerminus
										If mModifications(intIndex).TargetResiduesContain(C_TERMINAL_PROTEIN_SYMBOL_DMS) Then
											blnExistingModFound = True
										End If
									Case eResidueTerminusStateConstants.PeptideCTerminus
										If mModifications(intIndex).TargetResiduesContain(C_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
											blnExistingModFound = True
										End If
								End Select
							End If

							If Not blnExistingModFound AndAlso _
							 (eResidueTerminusState = eResidueTerminusStateConstants.ProteinNTerminus OrElse _
							  eResidueTerminusState = eResidueTerminusStateConstants.ProteinNandCCTerminus) Then

								' Protein N-Terminus residue could also match a Peptide N-terminal mod; check for this
								If mModifications(intIndex).TargetResiduesContain(N_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
									blnExistingModFound = True
								End If
							End If

							If Not blnExistingModFound AndAlso _
							 (eResidueTerminusState = eResidueTerminusStateConstants.ProteinCTerminus OrElse _
							  eResidueTerminusState = eResidueTerminusStateConstants.ProteinNandCCTerminus) Then

								' Protein C-Terminus residue could also match a Peptide C-terminal mod; check for this
								If mModifications(intIndex).TargetResiduesContain(C_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
									blnExistingModFound = True
								End If
							End If

						End If

						If blnExistingModFound Then
							' Match found
							Return mModifications(intIndex)
						End If
					End If
				End If
			Next intIndex
		End If

		' No match was found
		' First compare against modifications, only considering those with empty .TargetResidues
		' If still not match, then we'll try again but ignore .TargetResidues
		Dim blnConsiderTargetResidues As Boolean = True

		Do
			For intIndex = 0 To mModifications.Count - 1
				If mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
				   mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType Then

					If (blnConsiderTargetResidues And mModifications(intIndex).TargetResidues.Length = 0) OrElse Not blnConsiderTargetResidues Then
						If mModifications(intIndex).ModificationSymbol = chModificationSymbol Then
							' Matching mass found
							blnExistingModFound = True
							Return mModifications(intIndex)
						End If
					End If

				End If
			Next intIndex

			If blnConsiderTargetResidues Then
				' No match; try again, but ignore .TargetResidues
				blnConsiderTargetResidues = False
			Else				
				Exit Do
			End If
		Loop	

		' Still no match; return a default modification with a mass of 0
		objModificationDefinition = New clsModificationDefinition(chModificationSymbol, 0)
		objModificationDefinition.MassCorrectionTag = LookupMassCorrectionTagByMass(0)
		Return objModificationDefinition

	End Function

	''' <summary>
	''' Looks for an existing modification with the given modification mass and target residues
	''' </summary>
	''' <param name="dblModificationMass"></param>
	''' <param name="chTargetResidue">If defined, then returns the first modification with the given mass and containing the residue in .TargetResidues; if no match, then looks for the first modification with the given mass and no defined .TargetResidues</param>
	''' <param name="eResidueTerminusState"></param>
	''' <param name="blnExistingModFound"></param>
	''' <param name="blnAddToModificationListIfUnknown"></param>
	''' <param name="MassDigitsOfPrecision"></param>
	''' <returns>The best matched modification; if no match is found, then returns a newly created modification definition, adding it to mModifications if blnAddToModificationListIfUnknown = True</returns>
	''' <remarks>If chTargetResidue is nothing, then follows similar matching logic, but skips defined modifications with defined .TargetResidues</remarks>
	Public Function LookupModificationDefinitionByMass(ByVal dblModificationMass As Double, _
													   ByVal chTargetResidue As Char, _
													   ByVal eResidueTerminusState As eResidueTerminusStateConstants, _
													   ByRef blnExistingModFound As Boolean, _
													   ByVal blnAddToModificationListIfUnknown As Boolean, _
													   Optional ByVal MassDigitsOfPrecision As Byte = MASS_DIGITS_OF_PRECISION) As clsModificationDefinition
		Dim intIndex As Integer
		Dim intNewModIndex As Integer

		Dim objModificationDefinition As clsModificationDefinition
		Dim strTargetResidues As String

		If MassDigitsOfPrecision < 0 Then MassDigitsOfPrecision = 0

		blnExistingModFound = False
		If Not chTargetResidue = Nothing OrElse eResidueTerminusState <> eResidueTerminusStateConstants.None Then
			' The residue was provided and/or the residue is located at a peptide or protein terminus
			' First compare against modifications with 1 or more residues in .TargetResidues
			For intIndex = 0 To mModifications.Count - 1
				If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
					mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse _
					mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso _
					mModifications(intIndex).TargetResidues.Length > 0 Then
					If Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
						' Matching mass found
						' Now see if .TargetResidues contains chTargetResidue
						If Not chTargetResidue = Nothing AndAlso mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
							blnExistingModFound = True
						End If

						If Not blnExistingModFound AndAlso eResidueTerminusState <> eResidueTerminusStateConstants.None Then
							Select Case eResidueTerminusState
								Case eResidueTerminusStateConstants.PeptideNTerminus, eResidueTerminusStateConstants.ProteinNTerminus, eResidueTerminusStateConstants.ProteinNandCCTerminus
									If mModifications(intIndex).TargetResiduesContain(N_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
										blnExistingModFound = True
									End If
								Case eResidueTerminusStateConstants.PeptideCTerminus, eResidueTerminusStateConstants.ProteinCTerminus
									If mModifications(intIndex).TargetResiduesContain(C_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
										blnExistingModFound = True
									End If
							End Select
						End If

						If blnExistingModFound Then
							' Match found
							Return mModifications(intIndex)
						End If
					End If
				End If
			Next intIndex
		End If

		' No match was found
		' Compare against modifications with empty .TargetResidues
		For intIndex = 0 To mModifications.Count - 1
			If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
				mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse _
				mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso _
				mModifications(intIndex).TargetResidues.Length = 0 Then

				If Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
					' Matching mass found
					Return mModifications(intIndex)
				End If

			End If
		Next intIndex

		' Still no match; look for the modification mass and residue in mStandardRefinementModifications
		' Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing
		If Not chTargetResidue = Nothing Then
			For intIndex = 0 To mStandardRefinementModifications.Length - 1
				If Math.Round(Math.Abs(mStandardRefinementModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
					' Matching mass found
					' Now see if .TargetResidues contains chTargetResidue
					If mStandardRefinementModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
						blnExistingModFound = True

						objModificationDefinition = mStandardRefinementModifications(intIndex)
						objModificationDefinition.ModificationSymbol = clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL

						If blnAddToModificationListIfUnknown AndAlso mDefaultModificationSymbols.Count > 0 Then
							' Append objModificationDefinition to mModifications()
							intNewModIndex = AddModification(objModificationDefinition, True)
							If intNewModIndex >= 0 Then
								Return mModifications(intNewModIndex)
							Else
								Return objModificationDefinition
							End If
						End If

						' Either blnAddToModificationListIfUnknown = False or no more default modification symbols
						' Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
						Return objModificationDefinition
					End If
				End If
			Next intIndex
		End If

		' Still no match; define a new custom modification
		' First, need to populate strTargetResidues
		If chTargetResidue = Nothing Then
			strTargetResidues = String.Empty
		Else
			strTargetResidues = chTargetResidue
		End If

		If eResidueTerminusState <> eResidueTerminusStateConstants.None Then
			' Assume this is a terminus mod
			Select Case eResidueTerminusState
				Case eResidueTerminusStateConstants.PeptideNTerminus, eResidueTerminusStateConstants.ProteinNTerminus, eResidueTerminusStateConstants.ProteinNandCCTerminus
					strTargetResidues = N_TERMINAL_PEPTIDE_SYMBOL_DMS
				Case eResidueTerminusStateConstants.PeptideNTerminus, eResidueTerminusStateConstants.PeptideCTerminus
					strTargetResidues = C_TERMINAL_PEPTIDE_SYMBOL_DMS
				Case Else
					' This shouldn't occur
			End Select
		End If

		objModificationDefinition = New clsModificationDefinition( _
					clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL, _
					dblModificationMass, _
					strTargetResidues, _
					clsModificationDefinition.eModificationTypeConstants.DynamicMod, _
					LookupMassCorrectionTagByMass(dblModificationMass), _
					clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, _
					True)

		If blnAddToModificationListIfUnknown AndAlso mDefaultModificationSymbols.Count > 0 Then
			' Append objModificationDefinition to mModifications()
			intNewModIndex = AddModification(objModificationDefinition, True)
			If intNewModIndex >= 0 Then
				Return mModifications(intNewModIndex)
			Else
				Return objModificationDefinition
			End If
		End If

		' Either blnAddToModificationListIfUnknown = False or no more default modification symbols
		' Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
		Return objModificationDefinition

	End Function

	''' <summary>
	''' Looks for an existing modification with the given modification mass, modification type, and target residues
	''' </summary>
	''' <param name="dblModificationMass"></param>
	''' <param name="eModType"></param>
	''' <param name="chTargetResidue">If defined, then returns the first modification with the given mass and containing the residue in .TargetResidues; if no match, then looks for the first modification with the given mass and no defined .TargetResidues</param>
	''' <param name="eResidueTerminusState"></param>
	''' <param name="blnExistingModFound"></param>
	''' <param name="blnAddToModificationListIfUnknown"></param>
	''' <param name="MassDigitsOfPrecision"></param>
	''' <returns>The best matched modification; if no match is found, then returns a newly created modification definition, adding it to mModifications if blnAddToModificationListIfUnknown = True</returns>
	''' <remarks>If chTargetResidue is nothing, then follows similar matching logic, but skips defined modifications with defined .TargetResidues</remarks>
	Public Function LookupModificationDefinitionByMassAndModType( _
	  ByVal dblModificationMass As Double, _
	  ByVal eModType As clsModificationDefinition.eModificationTypeConstants, _
	  ByVal chTargetResidue As Char, _
	  ByVal eResidueTerminusState As eResidueTerminusStateConstants, _
	  ByRef blnExistingModFound As Boolean, _
	  ByVal blnAddToModificationListIfUnknown As Boolean, _
	  Optional ByVal MassDigitsOfPrecision As Byte = MASS_DIGITS_OF_PRECISION) As clsModificationDefinition

		' If chTargetResidue is defined, then returns the first modification with the given mass and containing the residue in .TargetResidues
		'  If no match is found, then looks for the first modification with the given mass and no defined .TargetResidues
		'  If no match is found, then returns a newly created modification definition, adding it to mModifications if blnAddToModificationListIfUnknown = True
		' If chTargetResidue is nothing, then follows similar logic, but skips defined modifications with defined .TargetResidues

		Dim intIndex As Integer
		Dim intNewModIndex As Integer

		Dim objModificationDefinition As clsModificationDefinition
		Dim strTargetResidues As String

		Dim chModSymbol As Char
		Dim blnUseNextAvailableModificationSymbol As Boolean

		If MassDigitsOfPrecision < 0 Then MassDigitsOfPrecision = 0

		Select Case eModType
			Case clsModificationDefinition.eModificationTypeConstants.StaticMod, clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod, clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
				chModSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
				blnUseNextAvailableModificationSymbol = False
			Case Else
				chModSymbol = clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL
				blnUseNextAvailableModificationSymbol = True
		End Select

		blnExistingModFound = False
		If Not chTargetResidue = Nothing OrElse eResidueTerminusState <> eResidueTerminusStateConstants.None Then
			' The residue was provided and/or the residue is located at a peptide or protein terminus
			' First compare against modifications with 1 or more residues in .TargetResidues
			For intIndex = 0 To mModifications.Count - 1
				If mModifications(intIndex).ModificationType = eModType AndAlso _
				 (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
				  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse _
				  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso _
				  mModifications(intIndex).TargetResidues.Length > 0 Then

					If Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
						' Matching mass found
						' Now see if .TargetResidues contains chTargetResidue
						If Not chTargetResidue = Nothing AndAlso mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
							blnExistingModFound = True
						End If

						If Not blnExistingModFound AndAlso eResidueTerminusState <> eResidueTerminusStateConstants.None Then
							Select Case eResidueTerminusState
								Case eResidueTerminusStateConstants.PeptideNTerminus, eResidueTerminusStateConstants.ProteinNTerminus, eResidueTerminusStateConstants.ProteinNandCCTerminus
									If mModifications(intIndex).TargetResiduesContain(N_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
										blnExistingModFound = True
									End If
								Case eResidueTerminusStateConstants.PeptideCTerminus, eResidueTerminusStateConstants.ProteinCTerminus
									If mModifications(intIndex).TargetResiduesContain(C_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
										blnExistingModFound = True
									End If
							End Select
						End If

						If blnExistingModFound Then
							' Match found
							Return mModifications(intIndex)
						End If
					End If
				End If
			Next intIndex
		End If

		' No match was found
		' Compare against modifications with empty .TargetResidues
		For intIndex = 0 To mModifications.Count - 1
			If mModifications(intIndex).ModificationType = eModType AndAlso _
			 (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
			  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse _
			  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso _
			  mModifications(intIndex).TargetResidues.Length = 0 Then

				If Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
					' Matching mass found
					Return mModifications(intIndex)
				End If

			End If
		Next intIndex

		' Still no match; look for the modification mass and residue in mStandardRefinementModifications
		' Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing
		If Not chTargetResidue = Nothing Then
			For intIndex = 0 To mStandardRefinementModifications.Length - 1
				If Math.Round(Math.Abs(mStandardRefinementModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
					' Matching mass found
					' Now see if .TargetResidues contains chTargetResidue
					If mStandardRefinementModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
						blnExistingModFound = True

						objModificationDefinition = mStandardRefinementModifications(intIndex)
						objModificationDefinition.ModificationSymbol = chModSymbol
						objModificationDefinition.ModificationType = eModType

						If blnAddToModificationListIfUnknown AndAlso mDefaultModificationSymbols.Count > 0 Then
							' Append objModificationDefinition to mModifications()
							intNewModIndex = AddModification(objModificationDefinition, True)
							If intNewModIndex >= 0 Then
								Return mModifications(intNewModIndex)
							Else
								Return objModificationDefinition
							End If
						End If

						' Either blnAddToModificationListIfUnknown = False or no more default modification symbols
						' Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
						Return objModificationDefinition
					End If
				End If
			Next intIndex
		End If

		' Still no match; define a new custom modification
		' First, need to populate strTargetResidues
		If chTargetResidue = Nothing Then
			strTargetResidues = String.Empty
		Else
			strTargetResidues = chTargetResidue
		End If

		If eResidueTerminusState <> eResidueTerminusStateConstants.None Then
			' Assume this is a terminus mod
			Select Case eResidueTerminusState
				Case eResidueTerminusStateConstants.PeptideNTerminus, eResidueTerminusStateConstants.ProteinNTerminus, eResidueTerminusStateConstants.ProteinNandCCTerminus
					strTargetResidues = N_TERMINAL_PEPTIDE_SYMBOL_DMS
				Case eResidueTerminusStateConstants.PeptideNTerminus, eResidueTerminusStateConstants.PeptideCTerminus
					strTargetResidues = C_TERMINAL_PEPTIDE_SYMBOL_DMS
				Case Else
					' This shouldn't occur
			End Select
		End If

		objModificationDefinition = New clsModificationDefinition( _
		   chModSymbol, _
		   dblModificationMass, _
		   strTargetResidues, _
		   eModType, _
		   LookupMassCorrectionTagByMass(dblModificationMass), _
		   clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, _
		   True)

		If blnAddToModificationListIfUnknown AndAlso mDefaultModificationSymbols.Count > 0 Then
			' Append objModificationDefinition to mModifications()

			intNewModIndex = AddModification(objModificationDefinition, blnUseNextAvailableModificationSymbol)
			If intNewModIndex >= 0 Then
				Return mModifications(intNewModIndex)
			Else
				Return objModificationDefinition
			End If
		End If

		' Either blnAddToModificationListIfUnknown = False or no more default modification symbols
		' Return objModificationDefinition, which has .ModificationSymbol = LAST_RESORT_MODIFICATION_SYMBOL
		Return objModificationDefinition

	End Function

	Private Sub InitializeLocalVariables()
		mErrorMessage = String.Empty
		SetDefaultMassCorrectionTags()

		mModifications = New System.Collections.Generic.List(Of clsModificationDefinition)

		' Note that this sub will call UpdateDefaultModificationSymbols()
		ClearModifications()

		UpdateStandardRefinementModifications()

	End Sub

	Public Function ReadMassCorrectionTagsFile(ByVal strFilePath As String, ByRef blnFileNotFound As Boolean) As Boolean
		Dim strLineIn As String
		Dim strSplitLine As String()

		Dim blnSuccess As Boolean

		Try
			' Open the mass correction tags file
			' It should have 2 columns, separated by tabs
			' Column 1 is the mass correction tag name
			' Column 2 is the monoisotopic mass for the mass correction (positive or negative number)

			If strFilePath Is Nothing OrElse strFilePath.Length = 0 Then
				SetDefaultMassCorrectionTags()
				blnSuccess = True
			ElseIf Not System.IO.File.Exists(strFilePath) Then
				mErrorMessage = "Mass CorrectionTags File Not Found: " & strFilePath
				SetDefaultMassCorrectionTags()
				blnFileNotFound = True
				blnSuccess = False
			Else
				Using srMassCorrectionTagsFile As System.IO.StreamReader = New System.IO.StreamReader(strFilePath)

					If mMassCorrectionTags Is Nothing Then
						mMassCorrectionTags = New Hashtable
					Else
						mMassCorrectionTags.Clear()
					End If

					Do While srMassCorrectionTagsFile.Peek() >= 0
						strLineIn = srMassCorrectionTagsFile.ReadLine()
						If Not strLineIn Is Nothing AndAlso strLineIn.Length > 0 Then
							strSplitLine = strLineIn.Split(ControlChars.Tab)

							If Not strSplitLine Is Nothing AndAlso strSplitLine.Length >= 2 Then
								' See if the first column contains 1 or more characters and if the second column contains a number
								' Note that StoreMassCorrectionTag() will trim spaces from the end of the mass correction tag names
								If strSplitLine(0).Trim.Length >= 1 AndAlso clsPHRPBaseClass.IsNumber(strSplitLine(1)) Then
									StoreMassCorrectionTag(strSplitLine(0), Double.Parse(strSplitLine(1)))
								End If
							End If
						End If
					Loop
				End Using

				blnSuccess = True

				If mMassCorrectionTags.Count = 0 Then
					SetDefaultMassCorrectionTags()
				End If
			End If
		Catch ex As Exception
			mErrorMessage = "Error reading Mass Correction Tags file (" & strFilePath & "): " & ex.Message
			blnSuccess = False
		End Try

		Return blnSuccess

	End Function

	Public Function ReadModificationDefinitionsFile(ByVal strFilePath As String, ByRef blnFileNotFound As Boolean) As Boolean

		Dim objModificationDefinition As clsModificationDefinition

		Dim strLineIn As String
		Dim strSplitLine As String()
		Dim strResidues As String
		Dim strResiduesClean As String

		Dim chChar As Char

		Dim blnValidMod As Boolean

		Dim blnSuccess As Boolean

		Try
			' Open the modification file
			' It should have 2 or more columns, separated by tabs
			' Column 1 is the modification symbol
			' Column 2 is the modification mass
			' Column 3, which is optional, is the residues and/or terminii that can be modified; if omitted, then the modification can apply to any residues or terminii
			'   For column 3, use 1 letter amino acid abbreviations; the residues can be a continous string, or can be separated by commas and/or spaces
			'   For column 3, use the *_SYMBOL_DMS constants for the terminii (< and > for the peptide terminii; [ and ] for the protein terminii)
			' Column 4, which is optional, specifies the type of modification: D, S, T, I, or P (corresponding to clsModificationDefinition.eModificationTypeConstants)
			' Column 5, which is optional, specifies the mass correction tag associated with the given modification

			If strFilePath Is Nothing OrElse strFilePath.Length = 0 Then
				ClearModifications()
				blnSuccess = True
			ElseIf Not System.IO.File.Exists(strFilePath) Then
				mErrorMessage = "Modification Definition File Not Found: " & strFilePath
				ClearModifications()
				blnFileNotFound = True
				blnSuccess = False
			Else
				Using srModificationFile As System.IO.StreamReader = New System.IO.StreamReader(strFilePath)

					ClearModifications()

					Do While srModificationFile.Peek() >= 0
						strLineIn = srModificationFile.ReadLine()
						If Not strLineIn Is Nothing AndAlso strLineIn.Length > 0 Then
							strSplitLine = strLineIn.Split(ControlChars.Tab)

							If Not strSplitLine Is Nothing AndAlso strSplitLine.Length >= 2 Then
								' See if the first column contains a single character and if the second column contains a number
								If strSplitLine(0).Trim.Length = 1 AndAlso clsPHRPBaseClass.IsNumber(strSplitLine(1)) Then

									objModificationDefinition = New clsModificationDefinition( _
									   strSplitLine(0).Trim.Chars(0), _
									   Double.Parse(strSplitLine(1)))

									With objModificationDefinition
										If strSplitLine.Length >= 3 Then
											' Parse the target residues list
											strResidues = strSplitLine(2).Trim.ToUpper

											strResiduesClean = String.Empty
											For Each chChar In strResidues
												If Char.IsUpper(chChar) Then
													strResiduesClean &= chChar
												ElseIf chChar = N_TERMINAL_PEPTIDE_SYMBOL_DMS Or _
												 chChar = C_TERMINAL_PEPTIDE_SYMBOL_DMS Or _
												 chChar = N_TERMINAL_PROTEIN_SYMBOL_DMS Or _
												 chChar = C_TERMINAL_PROTEIN_SYMBOL_DMS Then
													strResiduesClean &= chChar
												End If
											Next chChar

											If strResiduesClean.Length > 0 Then
												.TargetResidues = String.Copy(strResiduesClean)
											End If

											If strSplitLine.Length >= 4 Then
												' Store the modification type
												If strSplitLine(3).Trim.Length = 1 Then
													.ModificationType = clsModificationDefinition.ModificationSymbolToModificationType(strSplitLine(3).ToUpper.Trim.Chars(0))
												End If

												' If the .ModificationType is unknown, then change it to Dynamic
												If .ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType Then
													.ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod
												End If

												If strSplitLine.Length >= 5 Then
													.MassCorrectionTag = strSplitLine(4).Trim

													If strSplitLine.Length >= 6 Then
														strSplitLine(5) = strSplitLine(5).Trim
														If strSplitLine(5).Length > 0 Then
															.AffectedAtom = strSplitLine(5).Chars(0)
														Else
															.AffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
														End If
													End If
												End If
											End If
										End If

										' Check whether the modification type is Static and the .TargetResidues are one of: <>[]
										' If so, update the modification type as needed
										If Not .TargetResidues Is Nothing AndAlso .TargetResidues.Trim.Length = 1 AndAlso _
										   .ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
											If .TargetResidues.Chars(0) = N_TERMINAL_PEPTIDE_SYMBOL_DMS Or _
											   .TargetResidues.Chars(0) = C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
												.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
											ElseIf .TargetResidues.Chars(0) = N_TERMINAL_PROTEIN_SYMBOL_DMS Or _
											 .TargetResidues.Chars(0) = C_TERMINAL_PROTEIN_SYMBOL_DMS Then
												.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
											End If
										End If

										' Validate some of the settings if the modification type is IsotopicMod or TerminalPeptideStaticMod or ProteinTerminusStaticMod
										blnValidMod = True
										Select Case .ModificationType
											Case clsModificationDefinition.eModificationTypeConstants.IsotopicMod
												.ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
												If .AffectedAtom = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL Then
													blnValidMod = False
												End If
											Case clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
												.ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
												If .TargetResidues <> N_TERMINAL_PEPTIDE_SYMBOL_DMS AndAlso _
												   .TargetResidues <> C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
													blnValidMod = False
												End If
											Case clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
												.ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
												If .TargetResidues <> N_TERMINAL_PROTEIN_SYMBOL_DMS AndAlso _
												   .TargetResidues <> C_TERMINAL_PROTEIN_SYMBOL_DMS Then
													blnValidMod = False
												End If
											Case clsModificationDefinition.eModificationTypeConstants.UnknownType
												.ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod
										End Select

										If .MassCorrectionTag = clsModificationDefinition.INITIAL_UNKNOWN_MASS_CORRECTION_TAG_NAME Then
											' Try to determine the mass correction name
											.MassCorrectionTag = LookupMassCorrectionTagByMass(.ModificationMass)
										End If

									End With

									If blnValidMod Then
										AddModification(objModificationDefinition, False)
									End If
								End If
							End If
						End If
					Loop

				End Using

				' Note that this sub will call UpdateDefaultModificationSymbols()
				ValidateModificationsVsDefaultModificationSymbols()
				blnSuccess = True
			End If
		Catch ex As Exception
			mErrorMessage = "Error reading Modification Definition file (" & strFilePath & "): " & ex.Message
			blnSuccess = False
		End Try

		Return blnSuccess

	End Function

	Public Sub ResetOccurrenceCountStats()
		Dim intIndex As Integer

		For intIndex = 0 To mModifications.Count - 1
			mModifications(intIndex).OccurrenceCount = 0
		Next intIndex

	End Sub

	Public Sub SetDefaultMassCorrectionTags()

		Try
			If mMassCorrectionTags Is Nothing Then
				mMassCorrectionTags = New Hashtable
			Else
				mMassCorrectionTags.Clear()
			End If

			' Note: If the mass correction tag names in this list contain spaces at the 
			'       beginning or end, then function StoreMassCorrectionTag will remove them

			StoreMassCorrectionTag("4xDeut  ", 4.0251)
			StoreMassCorrectionTag("6C132N15", 8.0143)
			StoreMassCorrectionTag("6C134N15", 10.0085)
			StoreMassCorrectionTag("6xC13N15", 7.0172)
			StoreMassCorrectionTag("9xC13N15", 10.0273)
			StoreMassCorrectionTag("AcetAmid", 41.0266)
			StoreMassCorrectionTag("Acetyl  ", 42.0106)
			StoreMassCorrectionTag("AlkSulf ", -25.0316)
			StoreMassCorrectionTag("AlkSulf2", -9.0367)
			StoreMassCorrectionTag("AmAdipic", -1.031634)
			StoreMassCorrectionTag("Aminaton", 15.0109)
			StoreMassCorrectionTag("AmOxButa", -2.01565)
			StoreMassCorrectionTag("BioPeoAm", 356.1882)
			StoreMassCorrectionTag("Biotinyl", 89.0061)
			StoreMassCorrectionTag("Bromo   ", 77.9105)
			StoreMassCorrectionTag("C12_PIC ", 119.038)
			StoreMassCorrectionTag("C13_PIC ", 125.0581)
			StoreMassCorrectionTag("Carbamyl", 43.0059)
			StoreMassCorrectionTag("Chloro  ", 33.96103)
			StoreMassCorrectionTag("Cys_EDTI", 117.0248)
			StoreMassCorrectionTag("Cystnyl ", 119.0041)
			StoreMassCorrectionTag("DCAT_D0 ", 42.0375)
			StoreMassCorrectionTag("DCAT_D3 ", 44.9957)
			StoreMassCorrectionTag("Deamide ", 0.984)
			StoreMassCorrectionTag("DeutMeth", 17.0345)
			StoreMassCorrectionTag("DiAcet_K", 84.0212)
			StoreMassCorrectionTag("DiffDeut", 1.0063)
			StoreMassCorrectionTag("Dimethyl", 28.0314)
			StoreMassCorrectionTag("DuMtO18 ", 19.0387)
			StoreMassCorrectionTag("EDPHeavy", 253.1032)
			StoreMassCorrectionTag("EDPLight", 246.0819)
			StoreMassCorrectionTag("EDT_Addn", 75.9805)
			StoreMassCorrectionTag("EDT_D0  ", 177.1055)
			StoreMassCorrectionTag("EDT_D0C7", 253.0819)
			StoreMassCorrectionTag("EDT_D4C7", 257.0819)
			StoreMassCorrectionTag("EDT+Iodo", 133.002)
			StoreMassCorrectionTag("EtShD0  ", 44.0085)
			StoreMassCorrectionTag("FarnesC ", 204.3511)
			StoreMassCorrectionTag("Formyl  ", 27.9949)
			StoreMassCorrectionTag("Furylium", 78.01057)
			StoreMassCorrectionTag("GamGluAl", -43.053433)
			StoreMassCorrectionTag("GluCPD4 ", 550.1953)
			StoreMassCorrectionTag("Gluthone", 305.0682)
			StoreMassCorrectionTag("Guanid  ", 42.0218)
			StoreMassCorrectionTag("Heme_615", 615.1694)
			StoreMassCorrectionTag("Heme_617", 617)
			StoreMassCorrectionTag("Heme_Alt", 616.1772)
			StoreMassCorrectionTag("Heme_Sas", 614.1819)
			StoreMassCorrectionTag("HemeAddn", 616.18)
			StoreMassCorrectionTag("Hexosam ", 203.0794)
			StoreMassCorrectionTag("Hexose  ", 162.0528)
			StoreMassCorrectionTag("HMOSERLC", -48.1086)
			StoreMassCorrectionTag("ICAT_C12", 227.127)
			StoreMassCorrectionTag("ICAT_C13", 236.1572)
			StoreMassCorrectionTag("ICAT_D0 ", 442.225)
			StoreMassCorrectionTag("ICAT_D8 ", 450.2752)
			StoreMassCorrectionTag("IodoAcet", 57.0215)
			StoreMassCorrectionTag("IodoAcid", 58.0055)
			StoreMassCorrectionTag("Iso_C13 ", 1.00335)
			StoreMassCorrectionTag("Iso_N15 ", 0.9971)
			StoreMassCorrectionTag("Iso_O18 ", 2.0042)
			StoreMassCorrectionTag("itrac   ", 144.102063)
			StoreMassCorrectionTag("iTRAQ8  ", 304.2022)
			StoreMassCorrectionTag("LeuToMet", 17.9564)
			StoreMassCorrectionTag("Lipid2  ", 576.51)
			StoreMassCorrectionTag("Lipoyl  ", 188.033)
			StoreMassCorrectionTag("Mercury ", 199.9549)
			StoreMassCorrectionTag("Met_2O18", 18.0243)
			StoreMassCorrectionTag("Met_O18 ", 16.02)
			StoreMassCorrectionTag("Methyl  ", 14.0157)
			StoreMassCorrectionTag("MinusH2O", -18.0106)
			StoreMassCorrectionTag("MTSLAddn", 184.076)
			StoreMassCorrectionTag("NEM     ", 125.047679)
			StoreMassCorrectionTag("NH+10Da ", 25.0109)
			StoreMassCorrectionTag("NH3_Loss", -17.026549)
			StoreMassCorrectionTag("NHS_SS  ", 87.9983)
			StoreMassCorrectionTag("NHSLCBio", 339.1617)
			StoreMassCorrectionTag("NHSPEO4 ", 588.2465)
			StoreMassCorrectionTag("Nitrosyl", 28.9902)
			StoreMassCorrectionTag("NO2_Addn", 44.9851)
			StoreMassCorrectionTag("NO2_Alt ", 44.9975)
			StoreMassCorrectionTag("NO2+10Da", 54.9851)
			StoreMassCorrectionTag("None    ", 0)
			StoreMassCorrectionTag("OMinus2H", 13.9793)
			StoreMassCorrectionTag("One_O18 ", 2.0042)
			StoreMassCorrectionTag("Oxid_PEO", 430.1886)
			StoreMassCorrectionTag("PCGalNAz", 502.2023)
			StoreMassCorrectionTag("PEO     ", 414.1937)
			StoreMassCorrectionTag("PEO4Addn", 442.2553)
			StoreMassCorrectionTag("PhIATD0 ", 490.1742)
			StoreMassCorrectionTag("PhIATD4 ", 494.1993)
			StoreMassCorrectionTag("PhIATMod", 136.0017)
			StoreMassCorrectionTag("Phosph  ", 79.9663)
			StoreMassCorrectionTag("PhosphH ", 95.9612)
			StoreMassCorrectionTag("Plus1Oxy", 15.9949)
			StoreMassCorrectionTag("Plus2Oxy", 31.9898)
			StoreMassCorrectionTag("Plus3Oxy", 47.9847)
			StoreMassCorrectionTag("PMA     ", 278.0019)
			StoreMassCorrectionTag("Pro2Azet", -14.01565)
			StoreMassCorrectionTag("Propnyl ", 56.02621)
			StoreMassCorrectionTag("ROBLOSS ", -11.876)
			StoreMassCorrectionTag("SATA_Alk", 131.0041)
			StoreMassCorrectionTag("SATA_Hvy", 129.9963)
			StoreMassCorrectionTag("SATA_Lgt", 115.9932)
			StoreMassCorrectionTag("SATAIodo", 146.015)
			StoreMassCorrectionTag("SBEDBait", 88.01)
			StoreMassCorrectionTag("SBEDCapt", 547.22)
			StoreMassCorrectionTag("SelCandM", 47.9444)
			StoreMassCorrectionTag("SP_Heavy", 177.1583)
			StoreMassCorrectionTag("SP_Light", 170.1055)
			StoreMassCorrectionTag("Sucinate", 116.011)
			StoreMassCorrectionTag("Sulf-10 ", -35.0316)
			StoreMassCorrectionTag("Sulf2-10", -19.0367)
			StoreMassCorrectionTag("SulfoNHS", 226.0776)
			StoreMassCorrectionTag("SumoEstr", 498.2417)
			StoreMassCorrectionTag("Sumoylat", 484.226)
			StoreMassCorrectionTag("TriAcetK", 126.0318)
			StoreMassCorrectionTag("TriMeth ", 42.0471)
			StoreMassCorrectionTag("TrypOxy ", 3.9949)
			StoreMassCorrectionTag("TrypPD4 ", 494.74)
			StoreMassCorrectionTag("Two_O18 ", 4.0085)
			StoreMassCorrectionTag("Ubiq_02 ", 114.1)
			StoreMassCorrectionTag("Ubiq_H  ", 104.0473)
			StoreMassCorrectionTag("Ubiq_H02", 233.1583)
			StoreMassCorrectionTag("Ubiq_H03", 176.1255)
			StoreMassCorrectionTag("Ubiq_L  ", 100.016)
			StoreMassCorrectionTag("Ubiq_L02", 229.127)
			StoreMassCorrectionTag("Ubiq_L03", 172.0942)
			StoreMassCorrectionTag("UbiqLRGG", 383.228103)
			StoreMassCorrectionTag("ValToMet", 31.9721)

		Catch ex As Exception
			' Ignore errors here
		End Try

	End Sub

	Private Sub StoreMassCorrectionTag(ByVal strTagName As String, ByVal dblMass As Double)

		Try
			mMassCorrectionTags.Add(strTagName.Trim, dblMass)
		Catch ex As Exception
			' If a duplicate is tag is entered into the mMassCorrectionTags hashtable, an error will occur; we'll ignore the error
			' Ignore errors here
		End Try

	End Sub

	Private Sub UpdateDefaultModificationSymbols(ByVal strModificationChars As String)
		Dim chChar As Char

		Try
			If Not strModificationChars Is Nothing AndAlso strModificationChars.Length > 0 Then
				If mDefaultModificationSymbols Is Nothing Then
					mDefaultModificationSymbols = New System.Collections.Queue
				Else
					mDefaultModificationSymbols.Clear()
				End If

				' Populate mDefaultModificationSymbols, making sure no characters are duplicated
				' In addition, do not allow LAST_RESORT_MODIFICATION_SYMBOL or NO_SYMBOL_MODIFICATION_SYMBOL to be used
				For Each chChar In strModificationChars
					If chChar <> clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL AndAlso _
					   chChar <> clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL Then
						If Not mDefaultModificationSymbols.Contains(chChar) Then
							mDefaultModificationSymbols.Enqueue(chChar)
						End If
					End If
				Next chChar

			End If
		Catch ex As Exception
			' Ignore errors here
		End Try

	End Sub

	Private Sub UpdateStandardRefinementModifications()
		Dim dblModificationMass As Double

		ReDim mStandardRefinementModifications(1)

		dblModificationMass = -17.026549
		mStandardRefinementModifications(0) = New clsModificationDefinition( _
		  clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL, _
		  dblModificationMass, _
		  "Q", _
		  clsModificationDefinition.eModificationTypeConstants.DynamicMod, _
		  LookupMassCorrectionTagByMass(dblModificationMass))

		dblModificationMass = -18.0106
		mStandardRefinementModifications(1) = New clsModificationDefinition( _
		  clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL, _
		  dblModificationMass, _
		  "E", _
		  clsModificationDefinition.eModificationTypeConstants.DynamicMod, _
		  LookupMassCorrectionTagByMass(dblModificationMass))

	End Sub

	Private Sub ValidateModificationsVsDefaultModificationSymbols()

		Dim intIndex As Integer
		Dim intIndexCompare As Integer
		Dim intIndexCopy As Integer

		Dim intDefaultModificationSymbolCount As Integer
		Dim chDefaultModificationSymbols() As Char

		Try
			' Reset the default modification symbols list
			UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS)

			ReDim chDefaultModificationSymbols(mDefaultModificationSymbols.Count - 1)
			mDefaultModificationSymbols.ToArray.CopyTo(chDefaultModificationSymbols, 0)
			intDefaultModificationSymbolCount = chDefaultModificationSymbols.Length

			' Step through mModifications and make sure each of the modification symbols is not present in mDefaultModificationChars
			For intIndex = 0 To mModifications.Count - 1
				intIndexCompare = 0
				Do While intIndexCompare < intDefaultModificationSymbolCount
					If mModifications(intIndex).ModificationSymbol = chDefaultModificationSymbols(intIndexCompare) Then
						' Remove this symbol from chDefaultModificationSymbols
						For intIndexCopy = intIndexCompare To intDefaultModificationSymbolCount - 2
							chDefaultModificationSymbols(intIndexCopy) = chDefaultModificationSymbols(intIndexCopy + 1)
						Next intIndexCopy
						intDefaultModificationSymbolCount -= 1
					Else
						intIndexCompare += 1
					End If
				Loop
			Next intIndex

			If intDefaultModificationSymbolCount < mDefaultModificationSymbols.Count Then
				mDefaultModificationSymbols.Clear()
				For intIndex = 0 To intDefaultModificationSymbolCount - 1
					mDefaultModificationSymbols.Enqueue(chDefaultModificationSymbols(intIndex))
				Next intIndex
			End If
		Catch ex As Exception
			' Ignore errors here
		End Try

	End Sub

	Public Function VerifyModificationPresent(ByVal dblModificationMass As Double, ByVal strTargetResidues As String, ByVal eModificationType As clsModificationDefinition.eModificationTypeConstants, Optional ByVal MassDigitsOfPrecision As Integer = MASS_DIGITS_OF_PRECISION) As Boolean
		' Returns True if the modification was matched or was added
		' Returns False if an error

		' Look for mods in mModifications with matching .ModificationType, .ModificationMass (within tolerance),
		'   and .TargetResidues vs. udtModDefintion 
		' If not found, add a new entry to mModifications

		Dim intIndex As Integer

		Dim objModificationDefinition As clsModificationDefinition

		Dim blnMatchFound As Boolean

		If MassDigitsOfPrecision < 0 Then MassDigitsOfPrecision = 0

		Try
			blnMatchFound = False

			For intIndex = 0 To mModifications.Count - 1
				If mModifications(intIndex).ModificationType = eModificationType Then
					' Matching modification type
					If Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision) = 0 Then
						' Matching mass
						' Compare .TargetResidues
						blnMatchFound = clsModificationDefinition.EquivalentTargetResidues(mModifications(intIndex).TargetResidues, strTargetResidues, True)
						If blnMatchFound Then
							Exit For
						End If
					End If
				End If
			Next intIndex


			If Not blnMatchFound Then
				objModificationDefinition = New clsModificationDefinition(dblModificationMass, strTargetResidues, eModificationType)
				objModificationDefinition.MassCorrectionTag = LookupMassCorrectionTagByMass(dblModificationMass)

				' Append objModificationDefinition to mModifications()
				AddModification(objModificationDefinition, True)

				blnMatchFound = True
			End If

		Catch ex As Exception
			blnMatchFound = False
		End Try

		Return blnMatchFound

	End Function
End Class
