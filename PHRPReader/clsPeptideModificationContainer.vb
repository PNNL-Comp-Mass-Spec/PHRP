Option Strict On

' This class is used to track modifications that can be applied to peptides
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
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

Imports PHRPReader.clsAminoAcidModInfo
Imports System.IO

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
    Private mDefaultModificationSymbols As Queue

    ' List of known mass correction tags
    Private mMassCorrectionTags As Hashtable

    ' List of known modifications
    Private mModifications As List(Of clsModificationDefinition)

    Private mErrorMessage As String

    ' This array holds modifications that Sequest or XTandem will often use but for 
    ' which the auto-addition method sometimes incorrectly notes
    Private mStandardRefinementModifications As List(Of clsModificationDefinition)

    Private mConsiderModSymbolWhenFindingIdenticalMods As Boolean

    Private mIntegerMassCorrectionTagLookup As Dictionary(Of Integer, String)

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

    Public ReadOnly Property Modifications As List(Of clsModificationDefinition)
        Get
            Return mModifications
        End Get
    End Property

    Public Property ConsiderModSymbolWhenFindingIdenticalMods() As Boolean
        Get
            Return mConsiderModSymbolWhenFindingIdenticalMods
        End Get
        Set(value As Boolean)
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
    Private Function AddModification(objModificationDefinition As clsModificationDefinition, blnUseNextAvailableModificationSymbol As Boolean) As Integer

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
                        If .ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
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

    Private Function AddUnknownModification(
      dblModificationMass As Double,
      eModType As clsModificationDefinition.eModificationTypeConstants,
      chTargetResidue As Char,
      eResidueTerminusState As eResidueTerminusStateConstants,
      blnAddToModificationListIfUnknown As Boolean,
      blnUseNextAvailableModificationSymbol As Boolean,
      chModSymbol As Char,
      MassDigitsOfPrecision As Byte) As clsModificationDefinition

        Dim strTargetResidues As String
        Dim objModificationDefinition As clsModificationDefinition

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
                Case eResidueTerminusStateConstants.PeptideCTerminus, eResidueTerminusStateConstants.ProteinCTerminus
                    strTargetResidues = C_TERMINAL_PEPTIDE_SYMBOL_DMS
                Case Else
                    Throw New Exception("Unrecognized eResidueTerminusStateConstants enum: " & eResidueTerminusState)
            End Select
        End If

        If Not blnUseNextAvailableModificationSymbol Then
            chModSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
        End If

        Dim strMassCorrectionTag As String = LookupMassCorrectionTagByMass(dblModificationMass, MassDigitsOfPrecision, True, MassDigitsOfPrecision)

        objModificationDefinition = New clsModificationDefinition(
           chModSymbol,
           dblModificationMass,
           strTargetResidues,
           eModType,
           strMassCorrectionTag,
           clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL,
           True)

        If blnAddToModificationListIfUnknown Then

            ' Append objModificationDefinition to mModifications()
            Dim intNewModIndex As Integer
            If mDefaultModificationSymbols.Count > 0 AndAlso blnUseNextAvailableModificationSymbol Then
                intNewModIndex = AddModification(objModificationDefinition, blnUseNextAvailableModificationSymbol:=True)
            Else
                intNewModIndex = AddModification(objModificationDefinition, blnUseNextAvailableModificationSymbol:=False)
            End If

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

    Public Sub AppendStandardRefinmentModifications()
        Dim intIndex As Integer

        For intIndex = 0 To mStandardRefinementModifications.Count - 1
            With mStandardRefinementModifications(intIndex)
                VerifyModificationPresent(.ModificationMass, .TargetResidues, .ModificationType)
            End With
        Next intIndex
    End Sub

    Public Sub ClearModifications()
        UpdateDefaultModificationSymbols(DEFAULT_MODIFICATION_SYMBOLS)
        mModifications.Clear()
    End Sub

    ''' <summary>
    ''' Converts a modification mass to a generic 8 character name
    ''' The name will always start with + or - then will have the modification mass, rounded as necessary to give an 8 character name
    ''' </summary>
    ''' <param name="dblModificationMass"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function GenerateGenericModMassName(dblModificationMass As Double) As String

        Dim strModMassName As String
        Dim strFormatString As String
        Dim dblFormatDigits As Double
        Dim intFormatDigits As Integer
        Dim intMaxLength As Integer

        If Math.Abs(dblModificationMass) < Single.Epsilon Then
            Return "+0.00000"
        ElseIf dblModificationMass < -9999999 Then
            ' Modification mass is too negative; always return -9999999
            Return "-9999999"
        ElseIf dblModificationMass > 9999999 Then
            ' Modification mass is too positve; always return +9999999
            Return "+9999999"
        End If

        ' Determine the number of digits that we will display to the left of the decimal point
        dblFormatDigits = Math.Log10(Math.Abs(dblModificationMass))
        If Math.Abs(dblFormatDigits - CInt(dblFormatDigits)) < Single.Epsilon Then
            ' ModMass is a power of 10
            intFormatDigits = CInt(dblFormatDigits) + 1
        Else
            intFormatDigits = CInt(Math.Ceiling(dblFormatDigits))
        End If

        If intFormatDigits < 1 Then intFormatDigits = 1

        ' Generate the format string
        ' For example, a modification mass of 15.9994 will have strFormatString = "+00.0000"
        ' Negative modification masses do not start with a minus sign in the format string since Visual Studio auto-adds it
        ' Thus, a modification mass of -15.9994 will have strFormatString = "00.0000"
        strFormatString = New String("0"c, intFormatDigits)

        If dblModificationMass > 0 Then
            strFormatString = "+" & strFormatString
            intMaxLength = 8
        Else
            intMaxLength = 7
        End If

        If strFormatString.Length < intMaxLength Then
            strFormatString &= "."
        End If

        Do While strFormatString.Length() < intMaxLength
            strFormatString &= "0"
        Loop

        strModMassName = dblModificationMass.ToString(strFormatString)

        If strModMassName.Length < 8 AndAlso strModMassName.IndexOf("."c) < 0 Then
            strModMassName &= "."
        End If

        Do While strModMassName.Length < 8
            ' Append extra zeroes (this code will likely never be reached)
            strModMassName &= "0"
        Loop

        If strModMassName.Length > 8 Then
            Throw New ArgumentOutOfRangeException("Generated Mod Name is longer than 8 characters: " & strModMassName)
        End If

        Return strModMassName

    End Function

    ''' <summary>
    ''' Looks for the best match in mIntegerMassCorrectionTagLookup for dblModificationMass (which should be close to a integer value)
    ''' </summary>
    ''' <param name="dblModificationMass"></param>
    ''' <returns>The mass correction tag name if a match, otherwise nothing</returns>
    ''' <remarks></remarks>
    Private Function GetBestIntegerBasedMassCorrectionTag(dblModificationMass As Double) As String

        Dim strClosestMassCorrectionTag As String = String.Empty

        For Each tagOverride In mIntegerMassCorrectionTagLookup
            If Math.Abs(dblModificationMass - tagOverride.Key) < 0.0001 Then
                strClosestMassCorrectionTag = tagOverride.Value
                Exit For
            End If
        Next

        Return strClosestMassCorrectionTag

    End Function

    Public Function GetModificationByIndex(intIndex As Integer) As clsModificationDefinition
        If intIndex >= 0 And intIndex < mModifications.Count Then
            Return mModifications(intIndex)
        Else
            Return New clsModificationDefinition
        End If
    End Function

    Public Function GetModificationTypeByIndex(intIndex As Integer) As clsModificationDefinition.eModificationTypeConstants
        If intIndex >= 0 And intIndex < mModifications.Count Then
            Return mModifications(intIndex).ModificationType
        Else
            Return clsModificationDefinition.eModificationTypeConstants.UnknownType
        End If
    End Function

    Public Function LookupMassCorrectionTagByMass(dblModificationMass As Double) As String
        Const MassDigitsOfPrecision As Byte = MASS_DIGITS_OF_PRECISION
        Const blnAddToModificationListIfUnknown = True

        Return LookupMassCorrectionTagByMass(dblModificationMass, MassDigitsOfPrecision, blnAddToModificationListIfUnknown)
    End Function

    Public Function LookupMassCorrectionTagByMass(dblModificationMass As Double, MassDigitsOfPrecision As Byte) As String
        Const blnAddToModificationListIfUnknown = True

        Return LookupMassCorrectionTagByMass(dblModificationMass, MassDigitsOfPrecision, blnAddToModificationListIfUnknown)
    End Function

    Public Function LookupMassCorrectionTagByMass(dblModificationMass As Double, MassDigitsOfPrecision As Byte, blnAddToModificationListIfUnknown As Boolean) As String
        Const MassDigitsOfPrecisionLoose As Byte = 1

        Return LookupMassCorrectionTagByMass(dblModificationMass, MassDigitsOfPrecision, blnAddToModificationListIfUnknown, MassDigitsOfPrecisionLoose)
    End Function

    Public Function LookupMassCorrectionTagByMass(dblModificationMass As Double, MassDigitsOfPrecision As Byte, blnAddToModificationListIfUnknown As Boolean, MassDigitsOfPrecisionLoose As Byte) As String

        Dim objEnum As IDictionaryEnumerator

        Dim intMassDigitsOfPrecisionCurrent As Integer
        Dim intMassDigitsOfPrecisionStop As Integer

        Dim dblMassDiff As Double

        Dim strClosestMassCorrectionTag As String
        Dim dblClosestMassCorrectionTagMassDiff As Double

        If MassDigitsOfPrecision < 0 Then MassDigitsOfPrecision = 0

        If MassDigitsOfPrecision < MassDigitsOfPrecisionLoose Then MassDigitsOfPrecisionLoose = MassDigitsOfPrecision

        If MassDigitsOfPrecision >= MassDigitsOfPrecisionLoose Then
            intMassDigitsOfPrecisionStop = MassDigitsOfPrecisionLoose
        Else
            intMassDigitsOfPrecisionStop = MassDigitsOfPrecision
        End If

        For intMassDigitsOfPrecisionCurrent = MassDigitsOfPrecision To intMassDigitsOfPrecisionStop Step -1
            strClosestMassCorrectionTag = String.Empty
            dblClosestMassCorrectionTagMassDiff = Double.MaxValue

            Try

                If intMassDigitsOfPrecisionStop = 0 Then
                    strClosestMassCorrectionTag = GetBestIntegerBasedMassCorrectionTag(dblModificationMass)
                    If Not String.IsNullOrEmpty(strClosestMassCorrectionTag) Then
                        dblClosestMassCorrectionTagMassDiff = 0
                        Exit Try
                    End If
                End If

                ' First look for an exact match in mMassCorrectionTags
                objEnum = mMassCorrectionTags.GetEnumerator
                Do While objEnum.MoveNext
                    ' strMassCorrectionTag = CStr(objEnum.Key)
                    dblMassDiff = Math.Abs(dblModificationMass - CDbl(objEnum.Value))
                    If dblMassDiff < dblClosestMassCorrectionTagMassDiff Then
                        strClosestMassCorrectionTag = CStr(objEnum.Key)
                        dblClosestMassCorrectionTagMassDiff = dblMassDiff
                    End If
                Loop
            Catch ex As Exception
                ' Error enumerating through mMassCorrectionTags
            End Try

            If Math.Abs(Math.Round(dblClosestMassCorrectionTagMassDiff, intMassDigitsOfPrecisionCurrent)) < Single.Epsilon Then
                ' Match found
                Return strClosestMassCorrectionTag
            Else
                If intMassDigitsOfPrecisionCurrent > intMassDigitsOfPrecisionStop Then
                    ' Let the For loop go through another iteration to see if we find a match
                Else
                    ' Match not found
                    ' Name the modification based on the mod mass
                    strClosestMassCorrectionTag = GenerateGenericModMassName(dblModificationMass)

                    If blnAddToModificationListIfUnknown Then
                        Try
                            mMassCorrectionTags.Add(strClosestMassCorrectionTag, dblModificationMass)
                        Catch ex As Exception
                            ' This shouldn't happen; a match should have been found earlier in this function
                            ' Ignore the error
                        End Try
                    End If

                    Return strClosestMassCorrectionTag
                End If
            End If
        Next intMassDigitsOfPrecisionCurrent

        Return String.Empty

    End Function

    Public Function LookupDynamicModificationDefinitionByTargetInfo(chModificationSymbol As Char, chTargetResidue As Char, eResidueTerminusState As eResidueTerminusStateConstants, ByRef blnExistingModFound As Boolean) As clsModificationDefinition
        ' Looks for a modification of type .DynamicMod or type .UnknownType in mModifications having .ModificationSymbol = chModificationSymbol and chTargetResidue in .TargetResidues
        ' Note: If chModificationSymbol does not match any of the mods, then a modification with a mass of 0 is returned

        Dim intIndex As Integer

        Dim objModificationDefinition As clsModificationDefinition

        blnExistingModFound = False
        If Not chTargetResidue = Nothing OrElse eResidueTerminusState <> eResidueTerminusStateConstants.None Then
            ' The residue was provided and/or the residue is located at a peptide or protein terminus
            ' First compare against modifications with 1 or more residues in .TargetResidues
            For intIndex = 0 To mModifications.Count - 1
                If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
                   mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso
                   mModifications(intIndex).TargetResidues.Length > 0 Then
                    If mModifications(intIndex).ModificationSymbol = chModificationSymbol Then
                        ' Matching modification symbol found
                        ' Now see if .TargetResidues contains chTargetResidue
                        If mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
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

                            If Not blnExistingModFound AndAlso
                             (eResidueTerminusState = eResidueTerminusStateConstants.ProteinNTerminus OrElse
                              eResidueTerminusState = eResidueTerminusStateConstants.ProteinNandCCTerminus) Then

                                ' Protein N-Terminus residue could also match a Peptide N-terminal mod; check for this
                                If mModifications(intIndex).TargetResiduesContain(N_TERMINAL_PEPTIDE_SYMBOL_DMS) Then
                                    blnExistingModFound = True
                                End If
                            End If

                            If Not blnExistingModFound AndAlso
                             (eResidueTerminusState = eResidueTerminusStateConstants.ProteinCTerminus OrElse
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
        ' If still no match, then we'll try again but ignore .TargetResidues
        Dim blnConsiderTargetResidues = True

        Do
            For intIndex = 0 To mModifications.Count - 1
                If mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
                   mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType Then

                    If (blnConsiderTargetResidues AndAlso String.IsNullOrWhiteSpace(mModifications(intIndex).TargetResidues)) OrElse Not blnConsiderTargetResidues Then
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
    Public Function LookupModificationDefinitionByMass(dblModificationMass As Double,
       chTargetResidue As Char,
       eResidueTerminusState As eResidueTerminusStateConstants,
       ByRef blnExistingModFound As Boolean,
       blnAddToModificationListIfUnknown As Boolean,
       Optional MassDigitsOfPrecision As Byte = MASS_DIGITS_OF_PRECISION) As clsModificationDefinition

        Dim intIndex As Integer
        Dim intNewModIndex As Integer

        Dim objModificationDefinition As clsModificationDefinition

        If MassDigitsOfPrecision < 0 Then MassDigitsOfPrecision = 0

        blnExistingModFound = False
        If Not chTargetResidue = Nothing OrElse eResidueTerminusState <> eResidueTerminusStateConstants.None Then
            ' The residue was provided and/or the residue is located at a peptide or protein terminus
            ' First compare against modifications with 1 or more residues in .TargetResidues
            For intIndex = 0 To mModifications.Count - 1
                If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
                 mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse
                 mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso
                 mModifications(intIndex).TargetResidues.Length > 0 Then
                    If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
                        ' Matching mass found
                        ' Now see if .TargetResidues contains chTargetResidue
                        If mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
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
            If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
             mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse
             mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso
             String.IsNullOrWhiteSpace(mModifications(intIndex).TargetResidues) Then

                If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
                    ' Matching mass found
                    Return mModifications(intIndex)
                End If

            End If
        Next intIndex

        ' Still no match; look for the modification mass and residue in mStandardRefinementModifications
        ' Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing
        If Not chTargetResidue = Nothing Then
            For intIndex = 0 To mStandardRefinementModifications.Count - 1
                If Math.Abs(Math.Round(Math.Abs(mStandardRefinementModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
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

        ' Still no match
        ' Compare against dynamic and unknown-type modifications, but ignore .TargetResidues
        For intIndex = 0 To mModifications.Count - 1
            If (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
             mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) Then

                If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
                    ' Matching mass found
                    ' Assure that the target residues contain chTargetResidue
                    If Not chTargetResidue = Nothing AndAlso Not mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
                        mModifications(intIndex).TargetResidues &= chTargetResidue
                    End If

                    Return mModifications(intIndex)
                End If

            End If
        Next intIndex

        ' Still no match; define a new custom modification
        Const eModType As clsModificationDefinition.eModificationTypeConstants = clsModificationDefinition.eModificationTypeConstants.DynamicMod
        Const chModSymbol As Char = clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL
        Const blnUseNextAvailableModificationSymbol = True

        objModificationDefinition = AddUnknownModification(dblModificationMass, eModType, chTargetResidue, eResidueTerminusState, blnAddToModificationListIfUnknown, blnUseNextAvailableModificationSymbol, chModSymbol, MassDigitsOfPrecision)

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
    Public Function LookupModificationDefinitionByMassAndModType(
      dblModificationMass As Double,
      eModType As clsModificationDefinition.eModificationTypeConstants,
      chTargetResidue As Char,
      eResidueTerminusState As eResidueTerminusStateConstants,
      ByRef blnExistingModFound As Boolean,
      blnAddToModificationListIfUnknown As Boolean,
      Optional MassDigitsOfPrecision As Byte = MASS_DIGITS_OF_PRECISION) As clsModificationDefinition

        ' If chTargetResidue is defined, then returns the first modification with the given mass and containing the residue in .TargetResidues
        '  If no match is found, then looks for the first modification with the given mass and no defined .TargetResidues
        '  If no match is found, then looks for the first dynamic modification with the given mass, regardless of .TargetResidues
        '  If no match is found, then returns a newly created modification definition, adding it to mModifications if blnAddToModificationListIfUnknown = True
        ' If chTargetResidue is nothing, then follows similar logic, but skips defined modifications with defined .TargetResidues

        Dim intIndex As Integer

        Dim objModificationDefinition As clsModificationDefinition

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
                If mModifications(intIndex).ModificationType = eModType AndAlso
                 (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
                  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse
                  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse
                  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod OrElse
                  mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso
                  mModifications(intIndex).TargetResidues.Length > 0 Then

                    If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
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
            If mModifications(intIndex).ModificationType = eModType AndAlso
             (mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
              mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod OrElse
              mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse
              mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod OrElse
              mModifications(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType) AndAlso
              String.IsNullOrWhiteSpace(mModifications(intIndex).TargetResidues) Then

                If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
                    ' Matching mass found
                    Return mModifications(intIndex)
                End If

            End If
        Next intIndex

        ' Still no match; look for the modification mass and residue in mStandardRefinementModifications
        ' Note that N-Terminal or C-Terminal mods will have chTargetResidue = Nothing or chTargetResidue = '<' or chTargetResidue = '>'
        If Not chTargetResidue = Nothing Then
            For intIndex = 0 To mStandardRefinementModifications.Count - 1
                If Math.Abs(Math.Round(Math.Abs(mStandardRefinementModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
                    ' Matching mass found
                    ' Now see if .TargetResidues contains chTargetResidue
                    If mStandardRefinementModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
                        blnExistingModFound = True

                        objModificationDefinition = mStandardRefinementModifications(intIndex)
                        objModificationDefinition.ModificationSymbol = chModSymbol
                        objModificationDefinition.ModificationType = eModType

                        If blnAddToModificationListIfUnknown AndAlso mDefaultModificationSymbols.Count > 0 Then
                            ' Append objModificationDefinition to mModifications()
                            Dim intNewModIndex As Integer
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


        ' No match was found
        ' Compare against modifications of the same type, but ignore .TargetResidues
        For intIndex = 0 To mModifications.Count - 1
            If mModifications(intIndex).ModificationType = eModType Then

                If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
                    ' Matching mass found
                    ' Assure that the target residues contain chTargetResidue
                    If Not chTargetResidue = Nothing AndAlso Not mModifications(intIndex).TargetResiduesContain(chTargetResidue) Then
                        mModifications(intIndex).TargetResidues &= chTargetResidue
                    End If

                    Return mModifications(intIndex)
                End If

            End If
        Next intIndex

        ' Still no match; define a new custom modification
        objModificationDefinition = AddUnknownModification(dblModificationMass, eModType, chTargetResidue, eResidueTerminusState, blnAddToModificationListIfUnknown, blnUseNextAvailableModificationSymbol, chModSymbol, MassDigitsOfPrecision)

        Return objModificationDefinition

    End Function

    Private Sub InitializeLocalVariables()
        mErrorMessage = String.Empty
        SetDefaultMassCorrectionTags()

        mModifications = New List(Of clsModificationDefinition)

        ' Note that this sub will call UpdateDefaultModificationSymbols()
        ClearModifications()

        UpdateStandardRefinementModifications()

        UpdateIntegerBasedModificationMap()

    End Sub

    Public Function ReadMassCorrectionTagsFile(strFilePath As String, ByRef blnFileNotFound As Boolean) As Boolean
        Dim strLineIn As String
        Dim strSplitLine As String()

        Dim blnSuccess As Boolean

        Try
            ' Open the mass correction tags file
            ' It should have 2 columns, separated by tabs
            ' Column 1 is the mass correction tag name
            ' Column 2 is the monoisotopic mass for the mass correction (positive or negative number)

            If String.IsNullOrWhiteSpace(strFilePath) Then
                SetDefaultMassCorrectionTags()
                blnSuccess = True
            ElseIf Not File.Exists(strFilePath) Then
                mErrorMessage = "Mass CorrectionTags File Not Found: " & strFilePath
                SetDefaultMassCorrectionTags()
                blnFileNotFound = True
                blnSuccess = False
            Else
                Using srMassCorrectionTagsFile = New StreamReader(strFilePath)

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
                                If strSplitLine(0).Trim.Length >= 1 AndAlso clsPHRPParser.IsNumber(strSplitLine(1)) Then
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

    Public Function ReadModificationDefinitionsFile(strFilePath As String, ByRef blnFileNotFound As Boolean) As Boolean

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

            If String.IsNullOrWhiteSpace(strFilePath) Then
                ClearModifications()
                blnSuccess = True
            ElseIf Not File.Exists(strFilePath) Then
                mErrorMessage = "Modification Definition File Not Found: " & strFilePath
                ClearModifications()
                blnFileNotFound = True
                blnSuccess = False
            Else
                Using srModificationFile = New StreamReader(strFilePath)

                    ClearModifications()

                    Do While srModificationFile.Peek() >= 0
                        strLineIn = srModificationFile.ReadLine()
                        If Not strLineIn Is Nothing AndAlso strLineIn.Length > 0 Then
                            strSplitLine = strLineIn.Split(ControlChars.Tab)

                            If Not strSplitLine Is Nothing AndAlso strSplitLine.Length >= 2 Then
                                ' See if the first column contains a single character and if the second column contains a number
                                If strSplitLine(0).Trim.Length = 1 AndAlso clsPHRPParser.IsNumber(strSplitLine(1)) Then

                                    objModificationDefinition = New clsModificationDefinition(
                                       strSplitLine(0).Trim.Chars(0),
                                       Double.Parse(strSplitLine(1)))

                                    With objModificationDefinition
                                        If strSplitLine.Length >= 3 Then
                                            ' Parse the target residues list
                                            strResidues = strSplitLine(2).Trim.ToUpper

                                            strResiduesClean = String.Empty
                                            For Each chChar In strResidues
                                                If Char.IsUpper(chChar) Then
                                                    strResiduesClean &= chChar
                                                ElseIf chChar = N_TERMINAL_PEPTIDE_SYMBOL_DMS OrElse
                                                       chChar = C_TERMINAL_PEPTIDE_SYMBOL_DMS OrElse
                                                       chChar = N_TERMINAL_PROTEIN_SYMBOL_DMS OrElse
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
                                        If Not .TargetResidues Is Nothing AndAlso .TargetResidues.Trim.Length = 1 AndAlso
                                           .ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                                            If .TargetResidues.Chars(0) = N_TERMINAL_PEPTIDE_SYMBOL_DMS Or
                                               .TargetResidues.Chars(0) = C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                                                .ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                                            ElseIf .TargetResidues.Chars(0) = N_TERMINAL_PROTEIN_SYMBOL_DMS Or
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
                                                If .TargetResidues <> N_TERMINAL_PEPTIDE_SYMBOL_DMS AndAlso
                                                   .TargetResidues <> C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                                                    blnValidMod = False
                                                End If
                                            Case clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
                                                .ModificationSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL
                                                If .TargetResidues <> N_TERMINAL_PROTEIN_SYMBOL_DMS AndAlso
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

            ' Note: Function StoreMassCorrectionTag will remove spaces 
            ' from the beginning or end of the mass correction tag names
            StoreMassCorrectionTag("4xDeut  ", 4.025107)
            StoreMassCorrectionTag("6C134N15", 10.008269)
            StoreMassCorrectionTag("6xC13N15", 7.017164)
            StoreMassCorrectionTag("AcetAmid", 41.02655)
            StoreMassCorrectionTag("Acetyl  ", 42.010567)
            StoreMassCorrectionTag("Acrylmid", 71.037117)
            StoreMassCorrectionTag("ADPRibos", 541.061096)
            StoreMassCorrectionTag("AlkSulf ", -25.0316)
            StoreMassCorrectionTag("Aminaton", 15.010899)
            StoreMassCorrectionTag("AmOxButa", -2.01565)
            StoreMassCorrectionTag("Bromo   ", 77.910507)
            StoreMassCorrectionTag("BS3Olnk ", 156.078644)
            StoreMassCorrectionTag("C13DtFrm", 36.07567)
            StoreMassCorrectionTag("Carbamyl", 43.005814)
            StoreMassCorrectionTag("Cyano   ", 24.995249)
            StoreMassCorrectionTag("Cys-Dha ", -33.98772)
            StoreMassCorrectionTag("Cystnyl ", 119.004097)
            StoreMassCorrectionTag("Deamide ", 0.984016)
            StoreMassCorrectionTag("DeutForm", 32.056407)
            StoreMassCorrectionTag("DeutMeth", 17.034479)
            StoreMassCorrectionTag("Dimethyl", 28.0313)
            StoreMassCorrectionTag("DTBP_Alk", 144.03573)
            StoreMassCorrectionTag("Formyl  ", 27.994915)
            StoreMassCorrectionTag("GalNAFuc", 648.2603)
            StoreMassCorrectionTag("GalNAMan", 664.2551)
            StoreMassCorrectionTag("Gluthone", 305.068146)
            StoreMassCorrectionTag("Guanid  ", 42.021797)
            StoreMassCorrectionTag("Heme_615", 615.169458)
            StoreMassCorrectionTag("Hexosam ", 203.079376)
            StoreMassCorrectionTag("Hexose  ", 162.052826)
            StoreMassCorrectionTag("ICAT_D0 ", 442.225006)
            StoreMassCorrectionTag("ICAT_D8 ", 450.275208)
            StoreMassCorrectionTag("IodoAcet", 57.021465)
            StoreMassCorrectionTag("IodoAcid", 58.005478)
            StoreMassCorrectionTag("Iso_N15 ", 0.997035)
            StoreMassCorrectionTag("itrac   ", 144.102066)
            StoreMassCorrectionTag("iTRAQ8  ", 304.205353)
            StoreMassCorrectionTag("LeuToMet", 17.956421)
            StoreMassCorrectionTag("Lipid2  ", 576.51178)
            StoreMassCorrectionTag("Mercury ", 199.9549)
            StoreMassCorrectionTag("Met_O18 ", 16.028204)
            StoreMassCorrectionTag("Methyl  ", 14.01565)
            StoreMassCorrectionTag("Methylmn", 13.031634)
            StoreMassCorrectionTag("MinusH2O", -18.010565)
            StoreMassCorrectionTag("NEM     ", 125.047676)
            StoreMassCorrectionTag("NH3_Loss", -17.026548)
            StoreMassCorrectionTag("NHS_SS  ", 87.998283)
            StoreMassCorrectionTag("NO2_Addn", 44.985077)
            StoreMassCorrectionTag("None    ", 0)
            StoreMassCorrectionTag("OMinus2H", 13.979265)
            StoreMassCorrectionTag("One_C12 ", 12)
            StoreMassCorrectionTag("One_O18 ", 2.004246)
            StoreMassCorrectionTag("OxoAla  ", -17.992805)
            StoreMassCorrectionTag("palmtlic", 236.21402)
            StoreMassCorrectionTag("PCGalNAz", 502.202332)
            StoreMassCorrectionTag("PEO     ", 414.193695)
            StoreMassCorrectionTag("PhosAden", 329.052521)
            StoreMassCorrectionTag("Phosph  ", 79.966331)
            StoreMassCorrectionTag("PhosUrid", 306.025299)
            StoreMassCorrectionTag("Plus1Oxy", 15.994915)
            StoreMassCorrectionTag("Plus2Oxy", 31.989828)
            StoreMassCorrectionTag("Plus3Oxy", 47.984745)
            StoreMassCorrectionTag("Propnyl ", 56.026215)
            StoreMassCorrectionTag("Pyro-cmC", 39.994915)
            StoreMassCorrectionTag("SATA_Alk", 131.0041)
            StoreMassCorrectionTag("SATA_Lgt", 115.9932)
            StoreMassCorrectionTag("Sucinate", 116.010956)
            StoreMassCorrectionTag("SulfoNHS", 226.077591)
            StoreMassCorrectionTag("Sumoylat", 484.228149)
            StoreMassCorrectionTag("TMT0Tag ", 224.152481)
            StoreMassCorrectionTag("TMT6Tag ", 229.162933)
            StoreMassCorrectionTag("TriMeth ", 42.046951)
            StoreMassCorrectionTag("Two_O18 ", 4.008491)
            StoreMassCorrectionTag("Ubiq_02 ", 114.042931)
            StoreMassCorrectionTag("Ubiq_L  ", 100.016045)
            StoreMassCorrectionTag("ValToMet", 31.972071)

        Catch ex As Exception
            ' Ignore errors here
        End Try

    End Sub

    Private Sub StoreMassCorrectionTag(strTagName As String, dblMass As Double)

        Try
            mMassCorrectionTags.Add(strTagName.Trim, dblMass)
        Catch ex As Exception
            ' If a duplicate is tag is entered into the mMassCorrectionTags hashtable, an error will occur; we'll ignore the error
            ' Ignore errors here
        End Try

    End Sub

    Private Sub UpdateDefaultModificationSymbols(strModificationChars As String)
        Dim chChar As Char

        Try
            If Not strModificationChars Is Nothing AndAlso strModificationChars.Length > 0 Then
                If mDefaultModificationSymbols Is Nothing Then
                    mDefaultModificationSymbols = New Queue
                Else
                    mDefaultModificationSymbols.Clear()
                End If

                ' Populate mDefaultModificationSymbols, making sure no characters are duplicated
                ' In addition, do not allow LAST_RESORT_MODIFICATION_SYMBOL or NO_SYMBOL_MODIFICATION_SYMBOL to be used
                For Each chChar In strModificationChars
                    If chChar <> clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL AndAlso
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

    Private Sub UpdateIntegerBasedModificationMap()
        mIntegerMassCorrectionTagLookup = New Dictionary(Of Integer, String)

        mIntegerMassCorrectionTagLookup.Add(-18, "MinusH2O")
        mIntegerMassCorrectionTagLookup.Add(-17, "NH3_Loss")
        mIntegerMassCorrectionTagLookup.Add(-11, "AsnToCys")
        mIntegerMassCorrectionTagLookup.Add(-8, "HisToGlu")
        mIntegerMassCorrectionTagLookup.Add(-7, "TyrToArg")
        mIntegerMassCorrectionTagLookup.Add(-4, "ThrToPro")
        mIntegerMassCorrectionTagLookup.Add(-3, "MetToLys")
        mIntegerMassCorrectionTagLookup.Add(-1, "Dehydro")
        mIntegerMassCorrectionTagLookup.Add(1, "Deamide")
        mIntegerMassCorrectionTagLookup.Add(2, "GluToMet")
        mIntegerMassCorrectionTagLookup.Add(4, "TrypOxy")
        mIntegerMassCorrectionTagLookup.Add(5, "5C13")
        mIntegerMassCorrectionTagLookup.Add(6, "6C13")
        mIntegerMassCorrectionTagLookup.Add(10, "D10-Leu")
        mIntegerMassCorrectionTagLookup.Add(13, "Methylmn")
        mIntegerMassCorrectionTagLookup.Add(14, "Methyl")
        mIntegerMassCorrectionTagLookup.Add(16, "Plus1Oxy")
        mIntegerMassCorrectionTagLookup.Add(18, "LeuToMet")
        mIntegerMassCorrectionTagLookup.Add(25, "Cyano")
        mIntegerMassCorrectionTagLookup.Add(28, "Dimethyl")
        mIntegerMassCorrectionTagLookup.Add(32, "Plus2Oxy")
        mIntegerMassCorrectionTagLookup.Add(42, "Acetyl")
        mIntegerMassCorrectionTagLookup.Add(43, "Carbamyl")
        mIntegerMassCorrectionTagLookup.Add(45, "NO2_Addn")
        mIntegerMassCorrectionTagLookup.Add(48, "Plus3Oxy")
        mIntegerMassCorrectionTagLookup.Add(56, "Propnyl")
        mIntegerMassCorrectionTagLookup.Add(58, "IodoAcid")
        mIntegerMassCorrectionTagLookup.Add(80, "Phosph")
        mIntegerMassCorrectionTagLookup.Add(89, "Biotinyl")
        mIntegerMassCorrectionTagLookup.Add(96, "PhosphH")
        mIntegerMassCorrectionTagLookup.Add(104, "Ubiq_H")
        mIntegerMassCorrectionTagLookup.Add(116, "Sucinate")
        mIntegerMassCorrectionTagLookup.Add(119, "Cystnyl")
        mIntegerMassCorrectionTagLookup.Add(125, "NEM")
        mIntegerMassCorrectionTagLookup.Add(144, "itrac")
        mIntegerMassCorrectionTagLookup.Add(215, "MethylHg")
        mIntegerMassCorrectionTagLookup.Add(236, "ICAT_C13")
        mIntegerMassCorrectionTagLookup.Add(442, "ICAT_D0")

    End Sub

    Private Sub UpdateStandardRefinementModifications()
        Dim dblModificationMass As Double

        mStandardRefinementModifications = New List(Of clsModificationDefinition)

        dblModificationMass = -17.026549
        mStandardRefinementModifications.Add(New clsModificationDefinition(
          clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL,
          dblModificationMass,
          "Q",
          clsModificationDefinition.eModificationTypeConstants.DynamicMod,
          LookupMassCorrectionTagByMass(dblModificationMass)))

        dblModificationMass = -18.0106
        mStandardRefinementModifications.Add(New clsModificationDefinition(
          clsModificationDefinition.LAST_RESORT_MODIFICATION_SYMBOL,
          dblModificationMass,
          "E",
          clsModificationDefinition.eModificationTypeConstants.DynamicMod,
          LookupMassCorrectionTagByMass(dblModificationMass)))

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

    Public Function VerifyModificationPresent(dblModificationMass As Double, strTargetResidues As String, eModificationType As clsModificationDefinition.eModificationTypeConstants, Optional MassDigitsOfPrecision As Integer = MASS_DIGITS_OF_PRECISION) As Boolean
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
                    If Math.Abs(Math.Round(Math.Abs(mModifications(intIndex).ModificationMass - dblModificationMass), MassDigitsOfPrecision)) < Single.Epsilon Then
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
