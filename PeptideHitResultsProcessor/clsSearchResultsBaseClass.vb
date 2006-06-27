Option Explicit On 

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
' E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
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

Public Class clsSearchResultsBaseClass

#Region "Constants and Enums"
    Protected Const MASS_DIGITS_OF_PRECISION As Integer = 3
#End Region

#Region "Structures"
    Public Structure udtSearchResultModificationsType
        Public ModDefinition As clsModificationDefinition
        Public Residue As Char
        Public ResidueLocInPeptide As Integer                               ' Indicates the residue number modified; the first residue is at position 1
        Public ResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants
    End Structure
#End Region

#Region "Classwide Variables"
    ' Note: Many of these variables typically hold numbers but we're storing the numbers as strings
    '       to prevent the numeric representation from changing when converting to a number then back to a string
    Protected mResultID As Integer                              ' RowIndex for Synopsis/First Hits files; auto-assigned for XTandem
    Protected mGroupID As Integer                               ' Group ID assigned by XTandem
    Protected mScan As String
    Protected mCharge As String
    Protected mParentIonMH As String

    Protected mMultipleProteinCount As String                   ' Multiple protein count: 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
    Protected mProteinName As String
    Protected mProteinSeqResidueNumberStart As Integer          ' Typically always 1
    Protected mProteinSeqResidueNumberEnd As Integer            ' The residue number of the last residue in the protein's sequence; e.g. 100 if the protein has 100 residues total

    Protected mProteinExpectationValue As String                 ' Typically only used by XTandem; actually holds the Log of the expectation value
    Protected mProteinIntensity As String                        ' Typically only used by XTandem; actually holds the Log of the intensity

    Protected mPeptideLocInProteinStart As Integer              ' Position in the protein's residues of the first residue in the peptide
    Protected mPeptideLocInProteinEnd As Integer                ' Position in the protein's residues of the last residue in the peptide

    Protected mPeptidePreResidues As String                     ' Residue or residues before the start of the peptide sequence
    Protected mPeptidePostResidues As String                    ' Residue or residues after the end of the peptide sequence
    Protected mPeptideCleanSequence As String                   ' Peptide sequence without any modification symbols
    Protected mPeptideSequenceWithMods As String                ' Peptide sequence with modification symbols

    Protected mPeptideCleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants
    Protected mPeptideTerminusState As clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants

    Protected mPeptideMH As String                  ' In XTandem this is the monoisotopic MH; in Sequest is is the average mass MH
    Protected mPeptideDeltaMass As String           ' Difference in mass between the peptide's mass and the parent ion mass (i.e. the mass chosen for fragmentation)

    Protected mPeptideModDescription As String
    Protected mPeptideMonoisotopicMass As Double

    ' List of modifications present in the current peptide
    Protected mSearchResultModificationCount As Integer
    Protected mSearchResultModifications() As udtSearchResultModificationsType

    Protected mPeptideMods As clsPeptideModificationContainer

    Protected mPeptideCleavageStateCalculator As clsPeptideCleavageStateCalculator
    Protected mPeptideSeqMassCalculator As clsPeptideMassCalculator
#End Region

#Region "Properties"

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
    Public ReadOnly Property PeptideCleavageState() As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants
        Get
            Return mPeptideCleavageState
        End Get
    End Property
    Public ReadOnly Property PeptideTerminusState() As clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants
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
            Return mSearchResultModificationCount
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

        ' Initialize strSequenceWithMods to strCleanSequence; we'll insert the mod symbols below if mSearchResultModificationCount > 0
        strSequenceWithMods = String.Copy(strCleanSequence)

        If mSearchResultModificationCount > 0 Then
            ' Insert the modification symbols into strSequenceWithMods
            ' First, sort mSearchResultModifications on .ResidueLocInPeptide
            ' Limit the sort to the first mSearchResultModificationCount items
            If mSearchResultModificationCount > 1 Then
                Array.Sort(mSearchResultModifications, 0, mSearchResultModificationCount, New IResidueModificationInfoComparer)
            End If

            ' Now step backward through intResidueModificationPositions and add the symbols to strSequenceWithMods
            For intIndex = mSearchResultModificationCount - 1 To 0 Step -1
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

        mPeptideCleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific
        mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None

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
        Dim intIndex As Integer

        mSearchResultModificationCount = 0
        For intIndex = 0 To mSearchResultModifications.Length - 1
            With mSearchResultModifications(intIndex)
                .ModDefinition = Nothing
                .Residue = clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL
                .ResidueLocInPeptide = 0
                .ResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.None
            End With
        Next intIndex

    End Sub

    Public Sub ComputeMonoisotopicMass()
        Dim intIndex As Integer

        Dim chTargetResidue As Char
        Dim eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants

        ' Ths array is static to avoid re-reserving memory for it on every function call
        Static udtPeptideSequenceModInfo() As clsPeptideMassCalculator.udtPeptideSequenceModInfoType
        If udtPeptideSequenceModInfo Is Nothing Then
            ' Initially reserve space for 50 modifications
            ReDim udtPeptideSequenceModInfo(49)
        End If

        If mSearchResultModificationCount >= udtPeptideSequenceModInfo.Length Then
            ReDim udtPeptideSequenceModInfo(mSearchResultModificationCount - 1)
        End If

        ' Copy the mod info from mPeptideMods to udtPeptideSequenceModInfo
        For intIndex = 0 To mSearchResultModificationCount - 1
            With udtPeptideSequenceModInfo(intIndex)
                .ResidueLocInPeptide = mSearchResultModifications(intIndex).ResidueLocInPeptide
                .ModificationMass = mSearchResultModifications(intIndex).ModDefinition.ModificationMass
                .AffectedAtom = mSearchResultModifications(intIndex).ModDefinition.AffectedAtom
            End With
        Next intIndex

        mPeptideMonoisotopicMass = mPeptideSeqMassCalculator.ComputeSequenceMass(mPeptideCleanSequence, mSearchResultModificationCount, udtPeptideSequenceModInfo)

    End Sub

    Public Sub ComputePeptideCleavageStateInProtein()
        ' Determine the peptide's terminus state and cleavage state within the protein
        mPeptideCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues)
        mPeptideTerminusState = mPeptideCleavageStateCalculator.ComputeTerminusState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues)
    End Sub

    Public Function DetermineResidueTerminusState(ByVal intResidueLocInPeptide As Integer) As clsPeptideModificationContainer.eResidueTerminusStateConstants

        Dim eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants

        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.None
        If intResidueLocInPeptide = 1 Then
            ' Residue is at the peptide's N-terminus
            If mPeptideLocInProteinStart = mProteinSeqResidueNumberStart Then
                ' Residue is at the protein's N-terminus
                If mPeptideLocInProteinEnd = mProteinSeqResidueNumberEnd Then
                    ' The protein only contains one Residue, and we're currently examining it
                    eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinNandCCTerminus
                Else
                    eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinNTerminus
                End If
            Else
                eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideNTerminus
            End If
        Else
            If intResidueLocInPeptide = mPeptideLocInProteinEnd - mPeptideLocInProteinStart + 1 Then
                ' Residue is at the peptide's C-terminus
                If mPeptideLocInProteinEnd = mProteinSeqResidueNumberEnd Then
                    ' Residue is at the protein's C-terminus
                    eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinCTerminus
                Else
                    eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideCTerminus
                End If
            End If
        End If

        Return eResidueTerminusState
    End Function

    Public Function ExtractCleanSequenceFromSequenceWithMods(ByVal strSequenceWithMods As String, ByVal blnCheckForPrefixAndSuffixResidues As Boolean) As String
        Dim strPrimarySequence As String
        Dim strPrefix As String
        Dim strSuffix As String

        Dim chChar As Char
        Dim strCleanSequence As String

        strCleanSequence = String.Empty
        If Not strSequenceWithMods Is Nothing Then
            If blnCheckForPrefixAndSuffixResidues Then
                If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then
                    strPrimarySequence = String.Copy(strSequenceWithMods)
                End If
            Else
                strPrimarySequence = String.Copy(strSequenceWithMods)
            End If

            For Each chChar In strPrimarySequence
                If Char.IsLetter(chChar) Then
                    strCleanSequence &= chChar
                End If
            Next chChar
        End If

        Return strCleanSequence
    End Function

    Public Function GetSearchResultModDetailsByIndex(ByVal intIndex As Integer, ByRef chResidue As Char, ByRef intResidueLocInPeptide As Integer, ByRef dblModificationMass As Double, ByRef chAffectedAtom As Char) As Boolean
        ' Returns True if intIndex is valid; otherwise, returns false

        If intIndex >= 0 And intIndex < mSearchResultModificationCount Then
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

    Public Function GetSearchResultModDetailsByIndex(ByVal intIndex As Integer) As udtSearchResultModificationsType
        If intIndex >= 0 And intIndex < mSearchResultModificationCount Then
            Return mSearchResultModifications(intIndex)
        Else
            Return Nothing
        End If
    End Function

    Private Sub InitializeLocalVariables()
        ' Note: We initially reserve space for 10 search result modifications in mSearchResultModifications
        ' This array will be expanded if needed, but is never shrunk
        ' Use mSearchResultModificationCount to determine the number of modifications present in mSearchResultModifications
        mSearchResultModificationCount = 0
        ReDim mSearchResultModifications(9)

        ' Initialize mPeptideCleavageStateCalculator
        If mPeptideCleavageStateCalculator Is Nothing Then
            mPeptideCleavageStateCalculator = New clsPeptideCleavageStateCalculator
            mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants.Trypsin)
        End If

        ' Initialize mPeptideSeqMassCalculator
        If mPeptideSeqMassCalculator Is Nothing Then
            mPeptideSeqMassCalculator = New clsPeptideMassCalculator

            ' Set this to false to speed up the mass calculation speed
            ' It's OK to do this since this class always passes the clean sequence to .ComputeSequenceMass
            mPeptideSeqMassCalculator.RemovePrefixAndSuffixIfPresent = False
        End If

        Me.Clear()
    End Sub

    Public Sub SearchResultAddDynamicModification(ByVal chModificationSymbol As Char, ByVal chTargetResidue As Char, ByVal intResidueLocInPeptide As Integer, ByVal eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants, ByVal blnUpdateModOccurrenceCounts As Boolean)

        Dim objModificationDefinition As clsModificationDefinition
        Dim blnExistingModFound As Boolean

        blnExistingModFound = False

        ' Find the modification that uses this modification symbol and applies to this target residue
        objModificationDefinition = mPeptideMods.LookupDynamicModificationDefinitionByTargetInfo(chModificationSymbol, chTargetResidue, eResidueTerminusState, blnExistingModFound)

        If blnExistingModFound Then
            If intResidueLocInPeptide < 1 Then
                ' Invalid position; ignore this modification
            Else
                SearchResultAddModification( _
                                    objModificationDefinition, _
                                    chTargetResidue, _
                                    intResidueLocInPeptide, _
                                    eResidueTerminusState, _
                                    blnUpdateModOccurrenceCounts)
            End If
        End If

    End Sub

    Public Sub SearchResultAddModification(ByVal dblModificationMass As Double, ByVal chTargetResidue As Char, ByVal intResidueLocInPeptide As Integer, ByVal eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants, ByVal blnUpdateModOccurrenceCounts As Boolean)

        Dim objModificationDefinition As clsModificationDefinition
        Dim blnExistingModFound As Boolean

        If intResidueLocInPeptide < 1 Then
            ' Invalid position; ignore this modification
        Else
            ' Lookup the modification definition given the modification information
            objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(dblModificationMass, chTargetResidue, eResidueTerminusState, blnExistingModFound, True)

            SearchResultAddModification( _
                                objModificationDefinition, _
                                chTargetResidue, _
                                intResidueLocInPeptide, _
                                eResidueTerminusState, _
                                blnUpdateModOccurrenceCounts)
        End If

    End Sub

    Public Sub SearchResultAddModification(ByRef objModificationDefinition As clsModificationDefinition, ByVal chTargetResidue As Char, ByVal intResidueLocInPeptide As Integer, ByVal eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants, ByVal blnUpdateModOccurrenceCounts As Boolean)

        If intResidueLocInPeptide < 1 And objModificationDefinition.ModificationType <> clsModificationDefinition.eModificationTypeConstants.IsotopicMod Then
            ' Invalid position; ignore this modification
        Else
            ' Possibly expand mSearchResultModifications
            If mSearchResultModificationCount >= mSearchResultModifications.Length Then
                ReDim Preserve mSearchResultModifications(mSearchResultModifications.Length * 2 - 1)
            End If

            If blnUpdateModOccurrenceCounts Then
                ' Increment OccurenceCount
                objModificationDefinition.OccurrenceCount += 1
            End If

            With mSearchResultModifications(mSearchResultModificationCount)
                .ModDefinition = objModificationDefinition
                .Residue = chTargetResidue
                .ResidueLocInPeptide = intResidueLocInPeptide
                .ResidueTerminusState = eResidueTerminusState
            End With

            mSearchResultModificationCount += 1
        End If

    End Sub

    Public Sub SearchResultAddIsotopicModifications(ByVal blnUpdateModOccurrenceCounts As Boolean)
        Dim intIndex As Integer

        For intIndex = 0 To mPeptideMods.ModificationCount - 1
            If mPeptideMods.GetModificationTypeByIndex(intIndex) = clsModificationDefinition.eModificationTypeConstants.IsotopicMod Then
                SearchResultAddModification( _
                                    mPeptideMods.GetModificationByIndex(intIndex), _
                                    clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL, _
                                    0, _
                                    clsPeptideModificationContainer.eResidueTerminusStateConstants.None, _
                                    blnUpdateModOccurrenceCounts)
            End If
        Next intIndex

    End Sub

    Public Sub SearchResultAddStaticTerminusMods(ByVal blnAllowDuplicateModOnTerminus As Boolean, ByVal blnUpdateModOccurrenceCounts As Boolean)
        ' See if any protein or peptide terminus static mods are defined
        ' Peptide terminus static mods are always considered for a given peptide
        ' Protein terminus static mods are only considered if the peptide is at the appropriate terminus for the modification
        ' If blnAllowDuplicateModOnTerminus = False, then we only add the modification if the terminal residue does not already have the given modification associated with it

        Dim intModificationIndex As Integer
        Dim intIndexCompare As Integer
        Dim intResidueLocInPeptide As Integer
        Dim eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants

        Dim objModificationDefinition As clsModificationDefinition
        Dim dblMassDifference As Double

        Dim blnAddModification As Boolean

        For intModificationIndex = 0 To mPeptideMods.ModificationCount - 1
            blnAddModification = False
            If mPeptideMods.GetModificationTypeByIndex(intModificationIndex) = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod Then
                objModificationDefinition = mPeptideMods.GetModificationByIndex(intModificationIndex)
                If objModificationDefinition.TargetResidues = clsPeptideModificationContainer.N_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                    intResidueLocInPeptide = 1
                    If mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus Or _
                       mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
                        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinNTerminus
                    Else
                        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideNTerminus
                    End If
                    blnAddModification = True
                ElseIf objModificationDefinition.TargetResidues = clsPeptideModificationContainer.C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                    intResidueLocInPeptide = mPeptideCleanSequence.Length
                    If mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus Or _
                       mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
                        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinCTerminus
                    Else
                        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideCTerminus
                    End If
                    blnAddModification = True
                Else
                    ' Invalid target residue for a peptide terminus static mod; do not add the modification
                End If

            ElseIf mPeptideMods.GetModificationTypeByIndex(intModificationIndex) = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod Then
                objModificationDefinition = mPeptideMods.GetModificationByIndex(intModificationIndex)
                If objModificationDefinition.TargetResidues = clsPeptideModificationContainer.N_TERMINAL_PROTEIN_SYMBOL_DMS Then
                    If mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNTerminus Or _
                       mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
                        intResidueLocInPeptide = 1
                        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinNTerminus
                        blnAddModification = True
                    End If
                ElseIf objModificationDefinition.TargetResidues = clsPeptideModificationContainer.C_TERMINAL_PROTEIN_SYMBOL_DMS Then
                    If mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinCTerminus Or _
                       mPeptideTerminusState = clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
                        intResidueLocInPeptide = mPeptideCleanSequence.Length
                        eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.ProteinCTerminus
                        blnAddModification = True
                    End If
                Else
                    ' Invalid target residue for a protein terminus static mod; do not add the modification
                End If
            End If

            If blnAddModification Then
                If Not blnAllowDuplicateModOnTerminus Then
                    ' Look through udtResidueModificationInfo to see if this residue already has a modification with this modification's mass
                    ' If it does, do not re-add the modification
                    For intIndexCompare = 0 To mSearchResultModificationCount - 1
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
            chPrefix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST

            If Not mPeptidePreResidues Is Nothing Then
                strWork = mPeptidePreResidues.Trim
                If strWork.Length > 0 Then
                    chPrefix = strWork.Chars(strWork.Length - 1)
                    If chPrefix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_XTANDEM_NTerminus Then
                        chPrefix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
                    End If
                End If
            End If
        Catch ex As Exception
            ' Ignore errors here
        End Try

        Try
            chSuffix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
            If Not mPeptidePostResidues Is Nothing Then
                strWork = mPeptidePostResidues.Trim
                If strWork.Length > 0 Then
                    chSuffix = strWork.Chars(0)
                    If chSuffix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_XTANDEM_CTerminus Then
                        chSuffix = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
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

    Public Sub SetEnzymeMatchSpec(ByVal udtEnzymeMatchSpec As clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType)
        mPeptideCleavageStateCalculator.SetEnzymeMatchSpec(udtEnzymeMatchSpec.LeftResidueRegEx, udtEnzymeMatchSpec.RightResidueRegEx)
    End Sub

    Public Sub SetPeptideSequenceWithMods(ByVal strSequenceWithMods As String, ByVal blnCheckForPrefixAndSuffixResidues As Boolean, ByVal blnAutoPopulateCleanSequence As Boolean)
        ' Stores strSequenceWithMods in mPeptideSequenceWithMods
        ' If blnAutoPopulateCleanSequence = True, then populates mPeptideCleanSequence, 
        '  which automatically calls ComputePeptideCleavageStateInProtein

        Dim strPrimarySequence As String        ' Sequence with mods, but without the prefix or suffix residues
        Dim strPrefix As String
        Dim strSuffix As String

        If blnCheckForPrefixAndSuffixResidues Then
            If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then
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
            ' Note: This will call ComputePeptideCleavageStateInProtein()
            Me.PeptideCleanSequence = ExtractCleanSequenceFromSequenceWithMods(strPrimarySequence, False)
        End If

        mPeptideSequenceWithMods = String.Copy(strPrimarySequence)

    End Sub

    Public Sub SetStandardEnzymeMatchSpec(ByVal eStandardCleavageAgent As clsPeptideCleavageStateCalculator.eStandardCleavageAgentConstants)
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

        If mSearchResultModificationCount > 0 Then
            ReDim udtModNameAndResidueLoc(mSearchResultModificationCount - 1)
            ReDim intPointerArray(mSearchResultModificationCount - 1)

            If mSearchResultModificationCount = 1 Then
                intPointerArray(0) = 0
            Else
                ' Construct a pointer array so that we can search the modifications by .ResidueLocInPeptide
                For intIndex = 0 To mSearchResultModificationCount - 1
                    udtModNameAndResidueLoc(intIndex).ResidueLocInPeptide = mSearchResultModifications(intIndex).ResidueLocInPeptide
                    udtModNameAndResidueLoc(intIndex).ModName = mSearchResultModifications(intIndex).ModDefinition.MassCorrectionTag
                    intPointerArray(intIndex) = intIndex
                Next intIndex

                Array.Sort(udtModNameAndResidueLoc, intPointerArray, New clsPHRPBaseClass.IModNameAndResidueLocComparer)
            End If

            ' Step through the modifications and add the modification name and residue position to mPeptideModDescription
            ' Note that mods of type IsotopicMod will have .ResidueLocInPeptide = 0; other mods will have positive .ResidueLocInPeptide values
            For intIndex = 0 To mSearchResultModificationCount - 1
                With mSearchResultModifications(intPointerArray(intIndex))
                    If intIndex > 0 Then mPeptideModDescription &= MOD_LIST_SEP_CHAR
                    mPeptideModDescription &= .ModDefinition.MassCorrectionTag.Trim & ":"c & .ResidueLocInPeptide
                End With
            Next intIndex

        End If

    End Sub

    Protected Class IResidueModificationInfoComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtSearchResultModificationsType = CType(x, udtSearchResultModificationsType)
            Dim yData As udtSearchResultModificationsType = CType(y, udtSearchResultModificationsType)

            If xData.ResidueLocInPeptide > yData.ResidueLocInPeptide Then
                Return 1
            ElseIf xData.ResidueLocInPeptide < yData.ResidueLocInPeptide Then
                Return -1
            Else
                If xData.ModDefinition.MassCorrectionTag > yData.ModDefinition.MassCorrectionTag Then
                    Return 1
                ElseIf xData.ModDefinition.MassCorrectionTag < yData.ModDefinition.MassCorrectionTag Then
                    Return -1
                Else
                    Return 0
                End If

            End If
        End Function
    End Class
End Class
