Option Strict On

' This class can be used to track the peptide details for a given MS/MS search result
' It can track peptide residue level modifications and static, peptide-wide modifications
'
' Use SearchResultClearModifications() and SearchResultAddModification() to track
'  specific modifications for a given peptide
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 7, 2006
' Last updated June 26, 2006
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
Imports PHRPReader.clsPeptideCleavageStateCalculator

Public MustInherit Class clsSearchResultsBaseClass

#Region "Constants and Enums"
    Protected Const MASS_DIGITS_OF_PRECISION As Integer = 3
    Public Const MASS_C13 As Double = 1.00335483
#End Region

#Region "Structures"
	' Unused structure; deprecated in April 2012
	'Public Structure udtSearchResultModificationsType
	'    Public ModDefinition As clsModificationDefinition
	'    Public Residue As Char
	'    Public ResidueLocInPeptide As Integer                               ' Indicates the residue number modified; the first residue is at position 1
	'    Public ResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants
	'End Structure
#End Region

#Region "Classwide Variables"
	' Note: Many of these variables typically hold numbers but we're storing the numbers as strings
	'       to prevent the numeric representation from changing when converting to a number then back to a string
	Protected mResultID As Integer								' RowIndex for Synopsis/First Hits files; auto-assigned for XTandem, Inspect, and MSGFDB
	Protected mGroupID As Integer								' Group ID assigned by XTandem
	Protected mScan As String
	Protected mCharge As String
	Protected mParentIonMH As String

	Protected mMultipleProteinCount As String					' Multiple protein count: 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
	Protected mProteinName As String
	Protected mProteinSeqResidueNumberStart As Integer			' Typically always 1
	Protected mProteinSeqResidueNumberEnd As Integer			' The residue number of the last residue in the protein's sequence; e.g. 100 if the protein has 100 residues total

	Protected mProteinExpectationValue As String				 ' Typically only used by XTandem; actually holds the Log of the expectation value
	Protected mProteinIntensity As String						 ' Typically only used by XTandem; actually holds the Log of the intensity

	Protected mPeptideLocInProteinStart As Integer				' Position in the protein's residues of the first residue in the peptide
	Protected mPeptideLocInProteinEnd As Integer				' Position in the protein's residues of the last residue in the peptide

	Protected mPeptidePreResidues As String						' Residue or residues before the start of the peptide sequence
	Protected mPeptidePostResidues As String					' Residue or residues after the end of the peptide sequence
	Protected mPeptideCleanSequence As String					' Peptide sequence without any modification symbols
	Protected mPeptideSequenceWithMods As String				' Peptide sequence with modification symbols

	Protected mPeptideCleavageState As ePeptideCleavageStateConstants
	Protected mPeptideTerminusState As ePeptideTerminusStateConstants

	Protected mPeptideMH As String					' In XTandem this is the theoretical monoisotopic MH; in Sequest it was historically the average mass MH, though when a monoisotopic mass parent tolerance is specified, then this is a monoisotopic mass; in Inspect and MSGFDB, this is the theoretical monoisotopic MH; note that this is (M+H)+
	Protected mPeptideDeltaMass As String			' Difference in mass between the peptide's computed mass and the parent ion mass (i.e. the mass chosen for fragmentation); in Sequest this is Theoretical Mass - Observed Mass; The XTandem XML file stores DelM as Observed - Theoretical, but PHRP negates this to match Sequest; Inspect stores this value as Observed - Theoretical, but PHRP negates this to match Sequest; MSGFDB stores this value as Observed - Theoretical, but PHRP negates this to match Sequest

	'Protected mPeptideDeltaMassCorrectedPpm As Double         ' Computed using either mPeptideDeltaMass (negating to bring back to Observed minus Theoretical) or using PrecursorMass - mPeptideMonoisotopicMass; In either case, we must add/subtract 1 until value is between -0.5 and 0.5, then convert to ppm (using mPeptideMonoisotopicMass for ppm basis)

	Protected mPeptideModDescription As String
	Protected mPeptideMonoisotopicMass As Double				' Theoretical (computed) monoisotopic mass for a given peptide sequence, including any modified residues

	' List of modifications present in the current peptide
	Protected mSearchResultModifications As System.Collections.Generic.List(Of clsAminoAcidModInfo)

	' Possible modifications that the peptide could have
	Protected mPeptideMods As clsPeptideModificationContainer

	Protected mPeptideCleavageStateCalculator As clsPeptideCleavageStateCalculator
	Protected mPeptideSeqMassCalculator As clsPeptideMassCalculator

	Protected mErrorMessage As String = ""
#End Region

#Region "Properties"

	Public ReadOnly Property ErrorMessage() As String
		Get
			Return mErrorMessage
		End Get
	End Property
	Public Property ResultID() As Integer
		Get
			Return mResultID
		End Get
		Set(ByVal Value As Integer)
			mResultID = Value
		End Set
	End Property
	Public Property GroupID() As Integer
		Get
			Return mGroupID
		End Get
		Set(ByVal Value As Integer)
			mGroupID = Value
		End Set
	End Property
	Public Property Scan() As String
		Get
			Return mScan
		End Get
		Set(ByVal Value As String)
			mScan = Value
		End Set
	End Property
	Public Property Charge() As String
		Get
			Return mCharge
		End Get
		Set(ByVal Value As String)
			mCharge = Value
		End Set
	End Property
	Public Property ParentIonMH() As String
		Get
			Return mParentIonMH
		End Get
		Set(ByVal Value As String)
			mParentIonMH = Value
		End Set
	End Property
	Public Property MultipleProteinCount() As String
		Get
			Return mMultipleProteinCount
		End Get
		Set(ByVal Value As String)
			mMultipleProteinCount = Value
		End Set
	End Property
	Public Property ProteinName() As String
		Get
			Return mProteinName
		End Get
		Set(ByVal Value As String)
			mProteinName = Value
		End Set
	End Property
	Public Property ProteinExpectationValue() As String
		Get
			Return mProteinExpectationValue
		End Get
		Set(ByVal Value As String)
			mProteinExpectationValue = Value
		End Set
	End Property
	Public Property ProteinIntensity() As String
		Get
			Return mProteinIntensity
		End Get
		Set(ByVal Value As String)
			mProteinIntensity = Value
		End Set
	End Property
	Public Property ProteinSeqResidueNumberStart() As Integer
		Get
			Return mProteinSeqResidueNumberStart
		End Get
		Set(ByVal Value As Integer)
			mProteinSeqResidueNumberStart = Value
		End Set
	End Property
	Public Property ProteinSeqResidueNumberEnd() As Integer
		Get
			Return mProteinSeqResidueNumberEnd
		End Get
		Set(ByVal Value As Integer)
			mProteinSeqResidueNumberEnd = Value
		End Set
	End Property
	Public Property PeptideLocInProteinStart() As Integer
		Get
			Return mPeptideLocInProteinStart
		End Get
		Set(ByVal Value As Integer)
			mPeptideLocInProteinStart = Value
		End Set
	End Property
	Public Property PeptideLocInProteinEnd() As Integer
		Get
			Return mPeptideLocInProteinEnd
		End Get
		Set(ByVal Value As Integer)
			mPeptideLocInProteinEnd = Value
		End Set
	End Property
	Public Property PeptidePreResidues() As String
		Get
			Return mPeptidePreResidues
		End Get
		Set(ByVal Value As String)
			If Value Is Nothing Then Value = String.Empty
			mPeptidePreResidues = Value
			ComputePeptideCleavageStateInProtein()
		End Set
	End Property
	Public Property PeptidePostResidues() As String
		Get
			Return mPeptidePostResidues
		End Get
		Set(ByVal Value As String)
			If Value Is Nothing Then Value = String.Empty
			mPeptidePostResidues = Value
			ComputePeptideCleavageStateInProtein()
		End Set
	End Property
	Public Property PeptideCleanSequence() As String
		Get
			Return mPeptideCleanSequence
		End Get
		Set(ByVal Value As String)
			mPeptideCleanSequence = Value
			ComputePeptideCleavageStateInProtein()
		End Set
	End Property
	Public Property PeptideSequenceWithMods() As String
		Get
			Return mPeptideSequenceWithMods
		End Get
		Set(ByVal Value As String)
			mPeptideSequenceWithMods = Value
		End Set
	End Property
	Public ReadOnly Property PeptideCleavageState() As ePeptideCleavageStateConstants
		Get
			Return mPeptideCleavageState
		End Get
	End Property
	Public ReadOnly Property PeptideTerminusState() As ePeptideTerminusStateConstants
		Get
			Return mPeptideTerminusState
		End Get
	End Property
	Public Property PeptideMH() As String
		Get
			Return mPeptideMH
		End Get
		Set(ByVal Value As String)
			mPeptideMH = Value
		End Set
	End Property
	Public Property PeptideDeltaMass() As String
		Get
			Return mPeptideDeltaMass
		End Get
		Set(ByVal Value As String)
			mPeptideDeltaMass = Value
		End Set
	End Property
	'' Unused
	''Public Property PeptideDeltaMassCorrectedPpm() As Double
	''    Get
	''        Return mPeptideDeltaMassCorrectedPpm
	''    End Get
	''    Set(ByVal Value As Double)
	''        mPeptideDeltaMassCorrectedPpm = Value
	''    End Set
	''End Property
	Public Property PeptideModDescription() As String
		Get
			Return mPeptideModDescription
		End Get
		Set(ByVal Value As String)
			mPeptideModDescription = Value
		End Set
	End Property
	Public Property PeptideMonoisotopicMass() As Double
		Get
			Return mPeptideMonoisotopicMass
		End Get
		Set(ByVal Value As Double)
			mPeptideMonoisotopicMass = Value
		End Set
	End Property

	Public ReadOnly Property SearchResultModificationCount() As Integer
		Get
			Return mSearchResultModifications.Count
		End Get
	End Property
#End Region

	Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer)
		mPeptideMods = objPeptideMods

		InitializeLocalVariables()
	End Sub

	Public Function AddSearchResultModificationsToCleanSequence(ByVal strCleanSequence As String) As String
		' Generate the sequence with the mod symbols; returns the sequence

		Dim intIndex As Integer
		Dim strSequenceWithMods As String

		' Initialize strSequenceWithMods to strCleanSequence; we'll insert the mod symbols below if mSearchResultModifications.Count > 0
		strSequenceWithMods = String.Copy(strCleanSequence)

		If mSearchResultModifications.Count > 0 Then
			' Insert the modification symbols into strSequenceWithMods
			' First, sort mSearchResultModifications on .ResidueLocInPeptide and .MassCorrectionTag

			If mSearchResultModifications.Count > 1 Then
				mSearchResultModifications.Sort(New IGenericResidueModificationInfoComparer)
			End If

			' Now step backward through intResidueModificationPositions and add the symbols to strSequenceWithMods
			For intIndex = mSearchResultModifications.Count - 1 To 0 Step -1
				With mSearchResultModifications(intIndex)
					If .ModDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse _
					   .ModDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType Then
						strSequenceWithMods = strSequenceWithMods.Insert(.ResidueLocInPeptide, .ModDefinition.ModificationSymbol)
					End If
				End With
			Next intIndex
		End If

		Return strSequenceWithMods

	End Function

	Public Sub ApplyModificationInformation()
		' Populate mPeptideSequenceWithMods and .PeptideModDescription

		mPeptideSequenceWithMods = AddSearchResultModificationsToCleanSequence(mPeptideCleanSequence)
		UpdateModDescription()
	End Sub

	Public Overridable Sub Clear()
		mResultID = 0
		mGroupID = 0

		mScan = String.Empty
		mCharge = String.Empty
		mParentIonMH = String.Empty

		mMultipleProteinCount = String.Empty
		mProteinName = String.Empty
		mProteinExpectationValue = String.Empty
		mProteinIntensity = String.Empty

		ClearProteinSequenceInfo()

		ClearPeptideDetailsInfo()
	End Sub

	Public Sub ClearPeptideDetailsInfo()
		mPeptideLocInProteinStart = 0
		mPeptideLocInProteinEnd = 0

		mPeptidePreResidues = String.Empty
		mPeptidePostResidues = String.Empty
		mPeptideCleanSequence = String.Empty
		mPeptideSequenceWithMods = String.Empty

		mPeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
		mPeptideTerminusState = ePeptideTerminusStateConstants.None

		mPeptideMH = String.Empty
		mPeptideDeltaMass = String.Empty

		mPeptideModDescription = String.Empty
		mPeptideMonoisotopicMass = 0

		ClearSearchResultModifications()
	End Sub

	Public Sub ClearProteinSequenceInfo()
		mProteinSeqResidueNumberStart = 0
		mProteinSeqResidueNumberEnd = 0
	End Sub

	Public Sub ClearSearchResultModifications()

		mSearchResultModifications.Clear()

	End Sub

	Public Shared Function ComputeDelMCorrectedPPM( _
	  ByVal dblDelM As Double, _
	  ByVal dblPrecursorMonoMass As Double, _
	  ByVal blnAdjustPrecursorMassForC13 As Boolean, _
	  ByVal dblPeptideMonoisotopicMass As Double) As Double

		Dim intCorrectionCount As Integer = 0

		' Examine dblDelM to determine which isotope was chosen
		If dblDelM >= -0.5 Then
			' This is the typical case
			Do While dblDelM > 0.5
				dblDelM -= MASS_C13
				intCorrectionCount += 1
			Loop
		Else
			' This happens less often; but we'll still account for it
			' In this case, intCorrectionCount will be negative
			Do While dblDelM < -0.5
				dblDelM += MASS_C13
				intCorrectionCount -= 1
			Loop
		End If

		If intCorrectionCount <> 0 Then
			If blnAdjustPrecursorMassForC13 Then
				' Adjust the precursor mono mass based on intCorrectionCount
				dblPrecursorMonoMass -= intCorrectionCount * MASS_C13
			End If

			' Compute a new DelM value
			dblDelM = dblPrecursorMonoMass - dblPeptideMonoisotopicMass
		End If

		Return clsPeptideMassCalculator.MassToPPM(dblDelM, dblPeptideMonoisotopicMass)

	End Function

	''' <summary>
	''' Computes the theoretical monoisotopic mass for the given peptide
	''' Also updates mPeptideDeltaMassCorrectedPpm
	''' </summary>
	''' <remarks></remarks>
	Public Sub ComputeMonoisotopicMass()
		Dim intIndex As Integer

		' Ths array is static to avoid re-reserving memory for it on every function call
		Static udtPeptideSequenceModInfo() As clsPeptideMassCalculator.udtPeptideSequenceModInfoType
		If udtPeptideSequenceModInfo Is Nothing Then
			' Initially reserve space for 50 modifications
			ReDim udtPeptideSequenceModInfo(49)
		End If

		If mSearchResultModifications.Count >= udtPeptideSequenceModInfo.Length Then
			ReDim udtPeptideSequenceModInfo(mSearchResultModifications.Count - 1)
		End If

		' Copy the mod info from mPeptideMods to udtPeptideSequenceModInfo
		For intIndex = 0 To mSearchResultModifications.Count - 1
			With udtPeptideSequenceModInfo(intIndex)
				.ResidueLocInPeptide = mSearchResultModifications(intIndex).ResidueLocInPeptide
				.ModificationMass = mSearchResultModifications(intIndex).ModDefinition.ModificationMass
				.AffectedAtom = mSearchResultModifications(intIndex).ModDefinition.AffectedAtom
			End With
		Next intIndex

		mPeptideMonoisotopicMass = mPeptideSeqMassCalculator.ComputeSequenceMass(mPeptideCleanSequence, mSearchResultModifications.Count, udtPeptideSequenceModInfo)

		'' Unused
		' Update mPeptideDeltaMassCorrectedPpm
		''ComputeDelMCorrected()

	End Sub

	Public Sub ComputePeptideCleavageStateInProtein()
		' Determine the peptide's terminus state and cleavage state within the protein
		mPeptideCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues)
		mPeptideTerminusState = mPeptideCleavageStateCalculator.ComputeTerminusState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues)
	End Sub

	Public Function DetermineResidueTerminusState(ByVal intResidueLocInPeptide As Integer) As clsAminoAcidModInfo.eResidueTerminusStateConstants

		Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

		eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
		If intResidueLocInPeptide = 1 Then
			' Residue is at the peptide's N-terminus
			If mPeptideLocInProteinStart = mProteinSeqResidueNumberStart Then
				' Residue is at the protein's N-terminus
				If mPeptideLocInProteinEnd = mProteinSeqResidueNumberEnd Then
					' The protein only contains one Residue, and we're currently examining it
					eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus
				Else
					eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
				End If
			Else
				eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
			End If
		Else
			If intResidueLocInPeptide = mPeptideLocInProteinEnd - mPeptideLocInProteinStart + 1 Then
				' Residue is at the peptide's C-terminus
				If mPeptideLocInProteinEnd = mProteinSeqResidueNumberEnd Then
					' Residue is at the protein's C-terminus
					eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus
				Else
					eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
				End If
			End If
		End If

		Return eResidueTerminusState
	End Function

	Public Function GetSearchResultModDetailsByIndex(ByVal intIndex As Integer, ByRef chResidue As Char, ByRef intResidueLocInPeptide As Integer, ByRef dblModificationMass As Double, ByRef chAffectedAtom As Char) As Boolean
		' Returns True if intIndex is valid; otherwise, returns false

		If intIndex >= 0 And intIndex < mSearchResultModifications.Count Then
			With mSearchResultModifications(intIndex)
				chResidue = .Residue
				intResidueLocInPeptide = .ResidueLocInPeptide
				dblModificationMass = .ModDefinition.ModificationMass
				chAffectedAtom = .ModDefinition.AffectedAtom
			End With
			Return True
		Else
			Return False
		End If
	End Function

	Public Function GetSearchResultModDetailsByIndex(ByVal intIndex As Integer) As clsAminoAcidModInfo
		If intIndex >= 0 And intIndex < mSearchResultModifications.Count Then
			Return mSearchResultModifications(intIndex)
		Else
			Return Nothing
		End If
	End Function

	Private Sub InitializeLocalVariables()
		mSearchResultModifications = New System.Collections.Generic.List(Of clsAminoAcidModInfo)

		' Initialize mPeptideCleavageStateCalculator
		If mPeptideCleavageStateCalculator Is Nothing Then
			mPeptideCleavageStateCalculator = New clsPeptideCleavageStateCalculator
			mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants.Trypsin)
		End If

		' Initialize mPeptideSeqMassCalculator
		If mPeptideSeqMassCalculator Is Nothing Then
			mPeptideSeqMassCalculator = New clsPeptideMassCalculator

			' Set this to false to speed up the mass calculation speed
			' It's OK to do this since this class always passes the clean sequence to .ComputeSequenceMass
			mPeptideSeqMassCalculator.RemovePrefixAndSuffixIfPresent = False
		End If

		mErrorMessage = String.Empty

		Me.Clear()
	End Sub

	''' <summary>
	''' Associates the given modification symbol with the given residue
	''' If the modification symbol is unknown, then will return False
	''' </summary>
	''' <param name="chModificationSymbol"></param>
	''' <param name="chTargetResidue"></param>
	''' <param name="intResidueLocInPeptide"></param>
	''' <param name="eResidueTerminusState"></param>
	''' <param name="blnUpdateModOccurrenceCounts"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Function SearchResultAddDynamicModification(ByVal chModificationSymbol As Char, _
				   ByVal chTargetResidue As Char, _
				   ByVal intResidueLocInPeptide As Integer, _
				   ByVal eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
				   ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean

		Dim objModificationDefinition As clsModificationDefinition
		Dim blnExistingModFound As Boolean
		Dim blnSuccess As Boolean = False

		blnExistingModFound = False

		' Find the modification that uses this modification symbol and applies to this target residue
		objModificationDefinition = mPeptideMods.LookupDynamicModificationDefinitionByTargetInfo(chModificationSymbol, chTargetResidue, eResidueTerminusState, blnExistingModFound)

		If blnExistingModFound Then
			If intResidueLocInPeptide < 1 Then
				' Invalid position; ignore this modification
				mErrorMessage = "Invalid value for intResidueLocInPeptide: " & intResidueLocInPeptide.ToString
			Else
				blnSuccess = SearchResultAddModification( _
					 objModificationDefinition, _
					 chTargetResidue, _
					 intResidueLocInPeptide, _
					 eResidueTerminusState, _
					 blnUpdateModOccurrenceCounts)
			End If
		Else
			' Modification not found
			mErrorMessage = "Modification symbol not found: " & chModificationSymbol & "; TerminusState = " & eResidueTerminusState.ToString()
			blnSuccess = False
		End If

		Return blnSuccess

	End Function

	''' <summary>
	''' Associates the given modification mass with the given residue
	''' If the modification mass is unknown, then will auto-add it to the list of known modifications
	''' </summary>
	''' <param name="dblModificationMass"></param>
	''' <param name="chTargetResidue"></param>
	''' <param name="intResidueLocInPeptide"></param>
	''' <param name="eResidueTerminusState"></param>
	''' <param name="blnUpdateModOccurrenceCounts"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Function SearchResultAddModification(ByVal dblModificationMass As Double, _
			   ByVal chTargetResidue As Char, _
			   ByVal intResidueLocInPeptide As Integer, _
			   ByVal eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
			   ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean

		Dim objModificationDefinition As clsModificationDefinition
		Dim blnExistingModFound As Boolean
		Dim blnSuccess As Boolean = False

		If intResidueLocInPeptide < 1 Then
			' Invalid position; ignore this modification
			mErrorMessage = "Invalid value for intResidueLocInPeptide: " & intResidueLocInPeptide.ToString
		Else
			' Lookup the modification definition given the modification information
			' If the modification mass is unknown, then will auto-add it to the list of known modifications
			objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass( _
					 dblModificationMass, _
					 chTargetResidue, _
					 eResidueTerminusState, _
					 blnExistingModFound, _
					 True)

			blnSuccess = SearchResultAddModification( _
				 objModificationDefinition, _
				 chTargetResidue, _
				 intResidueLocInPeptide, _
				 eResidueTerminusState, _
				 blnUpdateModOccurrenceCounts)
		End If

		Return blnSuccess

	End Function

	''' <summary>
	''' Associates the given modification with the given residue
	''' </summary>
	''' <param name="objModificationDefinition"></param>
	''' <param name="chTargetResidue"></param>
	''' <param name="intResidueLocInPeptide"></param>
	''' <param name="eResidueTerminusState"></param>
	''' <param name="blnUpdateModOccurrenceCounts"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Function SearchResultAddModification(ByRef objModificationDefinition As clsModificationDefinition, _
			   ByVal chTargetResidue As Char, _
			   ByVal intResidueLocInPeptide As Integer, _
			   ByVal eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, _
			   ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean

		Dim blnSuccess As Boolean = False

		If intResidueLocInPeptide < 1 And objModificationDefinition.ModificationType <> clsModificationDefinition.eModificationTypeConstants.IsotopicMod Then
			' Invalid position; ignore this modification
			mErrorMessage = "Invalid value for intResidueLocInPeptide: " & intResidueLocInPeptide.ToString & " (objModificationDefinition.ModificationType = " & objModificationDefinition.ModificationType.ToString & ")"
		Else
			If blnUpdateModOccurrenceCounts Then
				' Increment OccurenceCount
				objModificationDefinition.OccurrenceCount += 1
			End If

			mSearchResultModifications.Add(New clsAminoAcidModInfo(chTargetResidue, intResidueLocInPeptide, eResidueTerminusState, objModificationDefinition))
			blnSuccess = True
		End If

		Return blnSuccess

	End Function

	''' <summary>
	''' Adds any defined isotopic modifications to the peptide
	''' </summary>
	''' <param name="blnUpdateModOccurrenceCounts"></param>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Function SearchResultAddIsotopicModifications(ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
		Dim intIndex As Integer
		Dim blnSuccess As Boolean = False

		For intIndex = 0 To mPeptideMods.ModificationCount - 1
			If mPeptideMods.GetModificationTypeByIndex(intIndex) = clsModificationDefinition.eModificationTypeConstants.IsotopicMod Then
				Dim intResidueLocInPeptide As Integer = 0

				blnSuccess = SearchResultAddModification( _
				  mPeptideMods.GetModificationByIndex(intIndex), _
				  clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, _
				  intResidueLocInPeptide, _
				  clsAminoAcidModInfo.eResidueTerminusStateConstants.None, _
				  blnUpdateModOccurrenceCounts)

			End If
		Next intIndex

		Return blnSuccess

	End Function

	Public Sub SearchResultAddStaticTerminusMods(ByVal blnAllowDuplicateModOnTerminus As Boolean, ByVal blnUpdateModOccurrenceCounts As Boolean)
		' See if any protein or peptide terminus static mods are defined
		' Peptide terminus static mods are always considered for a given peptide
		' Protein terminus static mods are only considered if the peptide is at the appropriate terminus for the modification
		' If blnAllowDuplicateModOnTerminus = False, then we only add the modification if the terminal residue does not already have the given modification associated with it

		Dim intModificationIndex As Integer
		Dim intIndexCompare As Integer
		Dim intResidueLocInPeptide As Integer
		Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

		Dim objModificationDefinition As clsModificationDefinition = Nothing
		Dim dblMassDifference As Double

		Dim blnAddModification As Boolean

		For intModificationIndex = 0 To mPeptideMods.ModificationCount - 1
			blnAddModification = False
			If mPeptideMods.GetModificationTypeByIndex(intModificationIndex) = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod Then
				objModificationDefinition = mPeptideMods.GetModificationByIndex(intModificationIndex)
				If objModificationDefinition.TargetResidues = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS Then
					intResidueLocInPeptide = 1
					If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus Or _
					mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
						eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
					Else
						eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
					End If
					blnAddModification = True
				ElseIf objModificationDefinition.TargetResidues = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
					intResidueLocInPeptide = mPeptideCleanSequence.Length
					If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus Or _
					mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
						eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus
					Else
						eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
					End If
					blnAddModification = True
				Else
					' Invalid target residue for a peptide terminus static mod; do not add the modification
				End If

			ElseIf mPeptideMods.GetModificationTypeByIndex(intModificationIndex) = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod Then
				objModificationDefinition = mPeptideMods.GetModificationByIndex(intModificationIndex)
				If objModificationDefinition.TargetResidues = clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS Then
					If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus Or _
					mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
						intResidueLocInPeptide = 1
						eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
						blnAddModification = True
					End If
				ElseIf objModificationDefinition.TargetResidues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS Then
					If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus Or _
					mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
						intResidueLocInPeptide = mPeptideCleanSequence.Length
						eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus
						blnAddModification = True
					End If
				Else
					' Invalid target residue for a protein terminus static mod; do not add the modification
				End If
			End If

			If blnAddModification And Not objModificationDefinition Is Nothing Then
				If Not blnAllowDuplicateModOnTerminus Then
					' Look through udtResidueModificationInfo to see if this residue already has a modification with this modification's mass
					' If it does, do not re-add the modification
					For intIndexCompare = 0 To mSearchResultModifications.Count - 1
						If objModificationDefinition Is mSearchResultModifications(intIndexCompare).ModDefinition Then
							blnAddModification = False
							Exit For
						ElseIf mSearchResultModifications(intIndexCompare).ResidueLocInPeptide = intResidueLocInPeptide Then
							' First compare the MassCorrectionTag names
							If mSearchResultModifications(intIndexCompare).ModDefinition.MassCorrectionTag = objModificationDefinition.MassCorrectionTag Then
								blnAddModification = False
								Exit For
							Else
								dblMassDifference = Math.Abs(mSearchResultModifications(intIndexCompare).ModDefinition.ModificationMass - objModificationDefinition.ModificationMass)
								If Math.Round(dblMassDifference, MASS_DIGITS_OF_PRECISION) = 0 Then
									blnAddModification = False
									Exit For
								End If
							End If

						End If
					Next intIndexCompare
				End If

				If blnAddModification Then
					SearchResultAddModification( _
					  objModificationDefinition, _
					  mPeptideCleanSequence.Chars(intResidueLocInPeptide - 1), _
					  intResidueLocInPeptide, _
					  eResidueTerminusState, _
					  blnUpdateModOccurrenceCounts)
				End If
			End If
		Next intModificationIndex

	End Sub

	Public Sub UpdatePeptideNTerminusMass(ByVal dblNTerminalMassChange As Double)
		' Updates the N-Terminal mass applied to peptides when computing their mass if it is significantly different than the
		'  currently defined N-terminal peptide mass
		If Math.Round(Math.Abs(dblNTerminalMassChange - mPeptideSeqMassCalculator.PeptideNTerminusMass), MASS_DIGITS_OF_PRECISION) > 0 Then
			mPeptideSeqMassCalculator.PeptideNTerminusMass = dblNTerminalMassChange
		End If
	End Sub

	Public Sub UpdatePeptideCTerminusMass(ByVal dblCTerminalMassChange As Double)
		' Updates the C-Terminal mass applied to peptides when computing their mass if significantly different than the
		'  currently defined C-terminal peptide mass
		If Math.Round(Math.Abs(dblCTerminalMassChange - mPeptideSeqMassCalculator.PeptideCTerminusMass), MASS_DIGITS_OF_PRECISION) > 0 Then
			mPeptideSeqMassCalculator.PeptideCTerminusMass = dblCTerminalMassChange
		End If
	End Sub

	Public Function SequenceWithPrefixAndSuffix(ByVal blnReturnSequenceWithMods As Boolean) As String
		' Note: Be sure to call ApplyModificationInformation before calling this function

		Dim strWork As String
		Dim chPrefix As Char
		Dim chSuffix As Char

		Try
			chPrefix = TERMINUS_SYMBOL_SEQUEST

			If Not mPeptidePreResidues Is Nothing Then
				strWork = mPeptidePreResidues.Trim
				If strWork.Length > 0 Then
					chPrefix = strWork.Chars(strWork.Length - 1)
					If chPrefix = TERMINUS_SYMBOL_XTANDEM_NTerminus Then
						chPrefix = TERMINUS_SYMBOL_SEQUEST
					End If
				End If
			End If
		Catch ex As Exception
			' Ignore errors here
		End Try

		Try
			chSuffix = TERMINUS_SYMBOL_SEQUEST
			If Not mPeptidePostResidues Is Nothing Then
				strWork = mPeptidePostResidues.Trim
				If strWork.Length > 0 Then
					chSuffix = strWork.Chars(0)
					If chSuffix = TERMINUS_SYMBOL_XTANDEM_CTerminus Then
						chSuffix = TERMINUS_SYMBOL_SEQUEST
					End If
				End If
			End If
		Catch ex As Exception
			' Ignore errors here
		End Try

		If blnReturnSequenceWithMods AndAlso Not mPeptideSequenceWithMods Is Nothing Then
			Return chPrefix & "." & mPeptideSequenceWithMods & "." & chSuffix
		Else
			If mPeptideCleanSequence Is Nothing Then
				Return String.Empty
			Else
				Return chPrefix & "." & mPeptideCleanSequence & "." & chSuffix
			End If
		End If
	End Function

	Public Sub SetEnzymeMatchSpec(ByVal udtEnzymeMatchSpec As udtEnzymeMatchSpecType)
		mPeptideCleavageStateCalculator.SetEnzymeMatchSpec(udtEnzymeMatchSpec.LeftResidueRegEx, udtEnzymeMatchSpec.RightResidueRegEx)
	End Sub

	Public Sub SetPeptideSequenceWithMods(ByVal strSequenceWithMods As String, ByVal blnCheckForPrefixAndSuffixResidues As Boolean, ByVal blnAutoPopulateCleanSequence As Boolean)
		' Stores strSequenceWithMods in mPeptideSequenceWithMods
		' If blnAutoPopulateCleanSequence = True, then populates mPeptideCleanSequence, 
		'  which automatically calls ComputePeptideCleavageStateInProtein

		Dim strPrimarySequence As String = String.Empty		  ' Sequence with mods, but without the prefix or suffix residues
		Dim strPrefix As String = String.Empty
		Dim strSuffix As String = String.Empty

		If blnCheckForPrefixAndSuffixResidues Then
			If Not SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then
				strPrimarySequence = String.Copy(strSequenceWithMods)
			Else
				If blnAutoPopulateCleanSequence Then
					Me.PeptidePreResidues = strPrefix
					Me.PeptidePostResidues = strSuffix
				End If
			End If
		Else
			strPrimarySequence = String.Copy(strSequenceWithMods)
		End If

		If blnAutoPopulateCleanSequence Then
			' Note: Property PeptideCleanSequence will call ComputePeptideCleavageStateInProtein()
			Me.PeptideCleanSequence = ExtractCleanSequenceFromSequenceWithMods(strPrimarySequence, False)
		End If

		mPeptideSequenceWithMods = String.Copy(strPrimarySequence)

	End Sub

	Public Sub SetStandardEnzymeMatchSpec(ByVal eStandardCleavageAgent As eStandardCleavageAgentConstants)
		mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(eStandardCleavageAgent)
	End Sub

	Public Sub UpdateModDescription()
		' Generate a comma separated list of the modifications present and the residue modified, with a colon separating the mod name and residue location
		' For example:
		'  Acetyl:1
		'  MinusH2O:1,Plus1Oxy:4,Plus1Oxy:19
		' The description is stored in mPeptideModDescription

		Const MOD_LIST_SEP_CHAR As Char = ","c
		Dim intIndex As Integer

		Dim udtModNameAndResidueLoc() As clsPHRPBaseClass.udtModNameAndResidueLocType
		Dim intPointerArray() As Integer

		mPeptideModDescription = String.Empty

		If mSearchResultModifications.Count > 0 Then
			ReDim udtModNameAndResidueLoc(mSearchResultModifications.Count - 1)
			ReDim intPointerArray(mSearchResultModifications.Count - 1)

			If mSearchResultModifications.Count = 1 Then
				intPointerArray(0) = 0
			Else
				' Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
				For intIndex = 0 To mSearchResultModifications.Count - 1
					udtModNameAndResidueLoc(intIndex).ResidueLocInPeptide = mSearchResultModifications(intIndex).ResidueLocInPeptide
					udtModNameAndResidueLoc(intIndex).ModName = mSearchResultModifications(intIndex).ModDefinition.MassCorrectionTag
					intPointerArray(intIndex) = intIndex
				Next intIndex

				Array.Sort(udtModNameAndResidueLoc, intPointerArray, New clsPHRPBaseClass.IModNameAndResidueLocComparer)
			End If

			' Step through the modifications and add the modification name and residue position to mPeptideModDescription
			' Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
			For intIndex = 0 To mSearchResultModifications.Count - 1
				With mSearchResultModifications(intPointerArray(intIndex))
					If intIndex > 0 Then mPeptideModDescription &= MOD_LIST_SEP_CHAR
					mPeptideModDescription &= .ModDefinition.MassCorrectionTag.Trim & ":"c & .ResidueLocInPeptide
				End With
			Next intIndex

		End If

	End Sub

	Protected Class IGenericResidueModificationInfoComparer
		Implements System.Collections.Generic.IComparer(Of clsAminoAcidModInfo)

		Public Function Compare(x As clsAminoAcidModInfo, y As clsAminoAcidModInfo) As Integer Implements System.Collections.Generic.IComparer(Of clsAminoAcidModInfo).Compare
			If x.ResidueLocInPeptide > y.ResidueLocInPeptide Then
				Return 1
			ElseIf x.ResidueLocInPeptide < y.ResidueLocInPeptide Then
				Return -1
			Else
				If x.ModDefinition.MassCorrectionTag > y.ModDefinition.MassCorrectionTag Then
					Return 1
				ElseIf x.ModDefinition.MassCorrectionTag < y.ModDefinition.MassCorrectionTag Then
					Return -1
				Else
					Return 0
				End If
			End If
		End Function

	End Class
End Class
