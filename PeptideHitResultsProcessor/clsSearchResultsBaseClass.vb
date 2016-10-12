Option Strict On

' This class is used to track the peptide details for a given MS/MS search result
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

Imports System.Runtime.InteropServices
Imports PHRPReader
Imports PHRPReader.clsPeptideCleavageStateCalculator

Public MustInherit Class clsSearchResultsBaseClass

#Region "Constants and Enums"
    ''' <summary>
    ''' Mass difference, in Da, between ^12C and ^13C
    ''' </summary>
    ''' <remarks></remarks>
    Public Const MASS_C13 As Double = 1.00335483
#End Region

#Region "Classwide Variables"
    ' Note: Many of these variables typically hold numbers but we're storing the numbers as strings
    '       to prevent the numeric representation from changing when converting to a number then back to a string

    ''' <summary>
    ''' Residue or residues before the start of the peptide sequence
    ''' </summary>
    ''' <remarks></remarks>
    Protected mPeptidePreResidues As String

    ''' <summary>
    ''' Residue or residues after the end of the peptide sequence
    ''' </summary>
    ''' <remarks></remarks>
    Protected mPeptidePostResidues As String

    ''' <summary>
    ''' Peptide sequence without any modification symbols
    ''' </summary>
    ''' <remarks></remarks>
    Protected mPeptideCleanSequence As String

    ''' <summary>
    ''' Cleavage state of the peptide
    ''' </summary>
    ''' <remarks></remarks>
    Protected mPeptideCleavageState As ePeptideCleavageStateConstants

    ''' <summary>
    ''' Terminus state of the peptide
    ''' </summary>
    ''' <remarks></remarks>
    Protected mPeptideTerminusState As ePeptideTerminusStateConstants

    ''' <summary>
    ''' List of modifications present in the current peptide
    ''' </summary>
    ''' <remarks></remarks>
    Protected ReadOnly mSearchResultModifications As List(Of clsAminoAcidModInfo)

    ''' <summary>
    ''' Possible modifications that the peptide could have
    ''' </summary>
    ''' <remarks></remarks>
    Protected ReadOnly mPeptideMods As clsPeptideModificationContainer

    Protected ReadOnly mPeptideCleavageStateCalculator As clsPeptideCleavageStateCalculator
    Protected ReadOnly mPeptideSeqMassCalculator As clsPeptideMassCalculator

    Protected mErrorMessage As String = ""
#End Region

#Region "Properties"

    ''' <summary>
    ''' Most recent error message
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ErrorMessage() As String
        Get
            Return mErrorMessage
        End Get
    End Property

    ''' <summary>
    ''' RowIndex for Synopsis/First Hits files; auto-assigned for XTandem, Inspect, MSGFDB, and MODa
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ResultID As Integer

    ''' <summary>
    ''' Group ID assigned by XTandem
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property GroupID As Integer

    ''' <summary>
    ''' Scan number
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property Scan As String

    ''' <summary>
    ''' Charge state
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property Charge As String

    ''' <summary>
    ''' Observed precursor m/z value converted to M+H
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ParentIonMH As String

    ''' <summary>
    ''' Multiple protein count: 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property MultipleProteinCount As String

    ''' <summary>
    ''' First protein for this PSM
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ProteinName As String

    ''' <summary>
    ''' Typically only used by XTandem; actually holds the Log of the expectation value
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ProteinExpectationValue As String

    ''' <summary>
    ''' Typically only used by XTandem; actually holds the Log of the intensity
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property ProteinIntensity As String

    ''' <summary>
    ''' Number of the first residue of a protein; typically always 1
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Initialized to 0 (unknown) then later changed 1 when ProteinSeqResidueNumberEnd is defined </remarks>
    Public Property ProteinSeqResidueNumberStart As Integer

    ''' <summary>
    ''' The residue number of the last residue in the protein's sequence
    ''' For example, 100 if the protein has 100 residues total
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>>Initialized to 0 (unknown)</remarks>
    Public Property ProteinSeqResidueNumberEnd As Integer

    ''' <summary>
    ''' Position in the protein's residues of the first residue in the peptide
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property PeptideLocInProteinStart As Integer

    ''' <summary>
    ''' Position in the protein's residues of the last residue in the peptide
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property PeptideLocInProteinEnd As Integer

    ''' <summary>
    '''  Residue or residues before the start of the peptide sequence
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Calls ComputePeptideCleavageStateInProtein</remarks>
    Public Property PeptidePreResidues() As String
        Get
            Return mPeptidePreResidues
        End Get
        Set(Value As String)
            If Value Is Nothing Then Value = String.Empty
            mPeptidePreResidues = Value
            ComputePeptideCleavageStateInProtein()
        End Set
    End Property

    ''' <summary>
    ''' Residue or residues after the end of the peptide sequence
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Calls ComputePeptideCleavageStateInProtein</remarks>
    Public Property PeptidePostResidues() As String
        Get
            Return mPeptidePostResidues
        End Get
        Set(Value As String)
            If Value Is Nothing Then Value = String.Empty
            mPeptidePostResidues = Value
            ComputePeptideCleavageStateInProtein()
        End Set
    End Property

    ''' <summary>
    ''' Peptide sequence without any modification symbols
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Calls ComputePeptideCleavageStateInProtein</remarks>
    Public Property PeptideCleanSequence() As String
        Get
            Return mPeptideCleanSequence
        End Get
        Set(Value As String)
            mPeptideCleanSequence = Value
            ComputePeptideCleavageStateInProtein()
        End Set
    End Property

    ''' <summary>
    ''' Peptide sequence with modification symbols
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>Mod symbols are single characters, like *, #, @, etc.</remarks>
    Public Property PeptideSequenceWithMods As String

    ''' <summary>
    ''' Cleavage state of the peptide
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property PeptideCleavageState() As ePeptideCleavageStateConstants
        Get
            Return mPeptideCleavageState
        End Get
    End Property

    ''' <summary>
    ''' Terminus state of the peptide
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property PeptideTerminusState() As ePeptideTerminusStateConstants
        Get
            Return mPeptideTerminusState
        End Get
    End Property

    ''' <summary>
    ''' In XTandem this is the theoretical monoisotopic MH
    ''' In Sequest it was historically the average mass MH, though when a monoisotopic mass parent tolerance is specified, then this is a monoisotopic mass
    ''' In Inspect, MSGF+, and MSAlign, this is the theoretical monoisotopic MH; note that this is (M+H)+
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property PeptideMH As String

    ''' <summary>
    ''' Difference in mass between the peptide's computed mass and the parent ion mass (i.e. the mass chosen for fragmentation)
    ''' In Sequest this is Theoretical Mass - Observed Mass
    ''' In XTandem, Inspect, MSGF+, and MSAlign the DelM value is listed as Observed - Theoretical, 
    ''' however, PHRP negates that value while reading the synopsis file to match Sequest
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property PeptideDeltaMass As String

    ''' <summary>
    ''' Comma separated list of the modifications present and the residue modified, with a colon separating the mod name and residue location
    ''' </summary>
    ''' <remarks>
    ''' Examples:
    '''   Acetyl:1
    '''   MinusH2O:1,Plus1Oxy:4,Plus1Oxy:19
    ''' </remarks>
    Public Property PeptideModDescription As String

    ''' <summary>
    ''' Theoretical (computed) monoisotopic mass for a given peptide sequence, including any modified residues
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Property PeptideMonoisotopicMass As Double

    ''' <summary>
    ''' Number of peptide residue modifications (dynamic or static)
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property SearchResultModificationCount() As Integer
        Get
            Return mSearchResultModifications.Count
        End Get
    End Property
#End Region

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="objPeptideMods"></param>
    ''' <param name="peptideSeqMassCalculator"></param>
    ''' <remarks></remarks>
    Public Sub New(objPeptideMods As clsPeptideModificationContainer, peptideSeqMassCalculator As clsPeptideMassCalculator)

        mSearchResultModifications = New List(Of clsAminoAcidModInfo)

        mPeptideMods = objPeptideMods

        mPeptideCleavageStateCalculator = New clsPeptideCleavageStateCalculator()
        mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(eStandardCleavageAgentConstants.Trypsin)

        If peptideSeqMassCalculator Is Nothing Then
            Throw New Exception("peptideSeqMassCalculator instance cannot be null")
        End If
        mPeptideSeqMassCalculator = peptideSeqMassCalculator

        InitializeLocalVariables()
    End Sub

    ''' <summary>
    ''' Add the modification symbols for the current peptide to the clean sequence
    ''' </summary>
    ''' <returns>Sequence with mod symbols</returns>
    ''' <remarks></remarks>
    Public Function AddSearchResultModificationsToCleanSequence() As String

        Dim intIndex As Integer
        Dim strSequenceWithMods As String

        ' Initialize strSequenceWithMods to the clean sequence; we'll insert the mod symbols below if mSearchResultModifications.Count > 0
        strSequenceWithMods = String.Copy(mPeptideCleanSequence)

        If mSearchResultModifications.Count > 0 Then
            ' Insert the modification symbols into strSequenceWithMods
            ' First, sort mSearchResultModifications on .ResidueLocInPeptide and .MassCorrectionTag

            If mSearchResultModifications.Count > 1 Then
                mSearchResultModifications.Sort(New IGenericResidueModificationInfoComparer)
            End If

            ' Now step backward through intResidueModificationPositions and add the symbols to strSequenceWithMods
            For intIndex = mSearchResultModifications.Count - 1 To 0 Step -1
                With mSearchResultModifications(intIndex)
                    If .ModDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod OrElse
                       .ModDefinition.ModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType Then
                        strSequenceWithMods = strSequenceWithMods.Insert(.ResidueLocInPeptide, .ModDefinition.ModificationSymbol)
                    End If
                End With
            Next intIndex
        End If

        Return strSequenceWithMods

    End Function

    ''' <summary>
    ''' Update .PeptideSequenceWithMods and .PeptideModDescription using the modifications defined for this peptide
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ApplyModificationInformation()

        PeptideSequenceWithMods = AddSearchResultModificationsToCleanSequence()
        UpdateModDescription()
    End Sub

    ''' <summary>
    ''' Clear all values
    ''' </summary>
    ''' <remarks></remarks>
    Public Overridable Sub Clear()
        ResultID = 0
        GroupID = 0

        Scan = String.Empty
        Charge = String.Empty
        ParentIonMH = String.Empty

        MultipleProteinCount = String.Empty
        ProteinName = String.Empty
        ProteinExpectationValue = String.Empty
        ProteinIntensity = String.Empty

        ClearProteinSequenceInfo()

        ' Note that this calls ClearSearchResultModifications()
        ClearPeptideDetailsInfo()

    End Sub

    ''' <summary>
    ''' Clear cached peptide information
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearPeptideDetailsInfo()
        PeptideLocInProteinStart = 0
        PeptideLocInProteinEnd = 0

        mPeptidePreResidues = String.Empty
        mPeptidePostResidues = String.Empty
        mPeptideCleanSequence = String.Empty
        PeptideSequenceWithMods = String.Empty

        mPeptideCleavageState = ePeptideCleavageStateConstants.NonSpecific
        mPeptideTerminusState = ePeptideTerminusStateConstants.None

        PeptideMH = String.Empty
        PeptideDeltaMass = String.Empty

        PeptideModDescription = String.Empty
        PeptideMonoisotopicMass = 0

        ClearSearchResultModifications()
    End Sub

    ''' <summary>
    ''' Clear the Protein sequence start residue and end residue values (ProteinSeqResidueNumberStart and ProteinSeqResidueNumberEnd)
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearProteinSequenceInfo()
        ProteinSeqResidueNumberStart = 0
        ProteinSeqResidueNumberEnd = 0
    End Sub

    ''' <summary>
    ''' Clear the modifications for this peptide
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearSearchResultModifications()
        mSearchResultModifications.Clear()
    End Sub

    ''' <summary>
    ''' Compute the delta mass, in ppm, optionally correcting for C13 isotopic selection errors
    ''' </summary>
    ''' <param name="dblDelM">Delta mass, in Da</param>
    ''' <param name="dblPrecursorMonoMass">Precursor monoisotopic mass</param>
    ''' <param name="blnAdjustPrecursorMassForC13">True to correct for C13 isotopic selection errors</param>
    ''' <param name="dblPeptideMonoisotopicMass"></param>
    ''' <returns>Delta mass, in ppm</returns>
    ''' <remarks></remarks>
    Public Shared Function ComputeDelMCorrectedPPM(
      dblDelM As Double,
      dblPrecursorMonoMass As Double,
      blnAdjustPrecursorMassForC13 As Boolean,
      dblPeptideMonoisotopicMass As Double) As Double

        Dim intCorrectionCount = 0

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

        Dim modifiedResidues = New List(Of clsPeptideMassCalculator.udtPeptideSequenceModInfoType)

        ' Copy the mod info from mSearchResultModifications to list modifiedResidues
        For Each searchResultMod In mSearchResultModifications
            Dim modifiedResidue = New clsPeptideMassCalculator.udtPeptideSequenceModInfoType
            modifiedResidue.ResidueLocInPeptide = searchResultMod.ResidueLocInPeptide
            modifiedResidue.ModificationMass = searchResultMod.ModDefinition.ModificationMass
            modifiedResidue.AffectedAtom = searchResultMod.ModDefinition.AffectedAtom

            modifiedResidues.Add(modifiedResidue)
        Next

        PeptideMonoisotopicMass = mPeptideSeqMassCalculator.ComputeSequenceMass(mPeptideCleanSequence, modifiedResidues)

    End Sub

    ''' <summary>
    ''' Update PeptideCleavageState and PeptideTerminusState based on PeptideCleanSequence, PeptidePreResidues and PeptidePostResidues
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ComputePeptideCleavageStateInProtein()
        ' Determine the peptide's terminus state and cleavage state within the protein
        mPeptideCleavageState = mPeptideCleavageStateCalculator.ComputeCleavageState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues)
        mPeptideTerminusState = mPeptideCleavageStateCalculator.ComputeTerminusState(mPeptideCleanSequence, mPeptidePreResidues, mPeptidePostResidues)
    End Sub

    ''' <summary>
    ''' Determine the terminus state for the given residue in the peptide
    ''' </summary>
    ''' <param name="intResidueLocInPeptide">Residue number (1-based)</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Function DetermineResidueTerminusState(intResidueLocInPeptide As Integer) As clsAminoAcidModInfo.eResidueTerminusStateConstants

        Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
        If intResidueLocInPeptide = 1 Then
            ' Residue is at the peptide's N-terminus
            If PeptideLocInProteinStart = ProteinSeqResidueNumberStart Then
                ' Residue is at the protein's N-terminus
                If PeptideLocInProteinEnd = ProteinSeqResidueNumberEnd Then
                    ' The protein only contains one Residue, and we're currently examining it
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNandCCTerminus
                Else
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
                End If
            Else
                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
            End If
        Else
            If intResidueLocInPeptide = PeptideLocInProteinEnd - PeptideLocInProteinStart + 1 Then
                ' Residue is at the peptide's C-terminus
                If PeptideLocInProteinEnd = ProteinSeqResidueNumberEnd Then
                    ' Residue is at the protein's C-terminus
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus
                Else
                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
                End If
            End If
        End If

        Return eResidueTerminusState
    End Function

    ''' <summary>
    ''' Get the modification info, by index
    ''' </summary>
    ''' <param name="intIndex">Modification entry index (0-based)</param>
    ''' <param name="chResidue"></param>
    ''' <param name="intResidueLocInPeptide"></param>
    ''' <param name="dblModificationMass"></param>
    ''' <param name="chAffectedAtom"></param>
    ''' <returns></returns>
    ''' <remarks>Returns True if intIndex is valid; otherwise, returns false</remarks>
    Public Function GetSearchResultModDetailsByIndex(
      intIndex As Integer,
      <Out()> ByRef chResidue As Char,
      <Out()> ByRef intResidueLocInPeptide As Integer,
      <Out()> ByRef dblModificationMass As Double,
      <Out()> ByRef chAffectedAtom As Char) As Boolean

        If intIndex >= 0 And intIndex < mSearchResultModifications.Count Then
            With mSearchResultModifications(intIndex)
                chResidue = .Residue
                intResidueLocInPeptide = .ResidueLocInPeptide
                dblModificationMass = .ModDefinition.ModificationMass
                chAffectedAtom = .ModDefinition.AffectedAtom
            End With
            Return True
        End If

        chResidue = Convert.ToChar(0)
        intResidueLocInPeptide = 0
        dblModificationMass = 0
        chAffectedAtom = Convert.ToChar(0)
        Return False

    End Function

    ''' <summary>
    ''' Get the modification info, by index
    ''' </summary>
    ''' <param name="intIndex">Modification entry index (0-based)</param>
    ''' <returns>Modification info details if a valid index, otherwise nothing</returns>
    ''' <remarks></remarks>
    Public Function GetSearchResultModDetailsByIndex(intIndex As Integer) As clsAminoAcidModInfo
        If intIndex >= 0 And intIndex < mSearchResultModifications.Count Then
            Return mSearchResultModifications(intIndex)
        Else
            Return Nothing
        End If
    End Function

    Private Sub InitializeLocalVariables()
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
    Public Function SearchResultAddDynamicModification(
       chModificationSymbol As Char,
       chTargetResidue As Char,
       intResidueLocInPeptide As Integer,
       eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
       blnUpdateModOccurrenceCounts As Boolean) As Boolean

        Dim objModificationDefinition As clsModificationDefinition
        Dim blnExistingModFound As Boolean
        Dim blnSuccess = False

        blnExistingModFound = False

        ' Find the modification that uses this modification symbol and applies to this target residue
        objModificationDefinition = mPeptideMods.LookupDynamicModificationDefinitionByTargetInfo(chModificationSymbol, chTargetResidue, eResidueTerminusState, blnExistingModFound)

        If blnExistingModFound Then
            If intResidueLocInPeptide < 1 Then
                ' Invalid position; ignore this modification
                mErrorMessage = "Invalid value for intResidueLocInPeptide: " & intResidueLocInPeptide.ToString
            Else
                blnSuccess = SearchResultAddModification(
                     objModificationDefinition,
                     chTargetResidue,
                     intResidueLocInPeptide,
                     eResidueTerminusState,
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
    ''' <returns>True if mod successfully added</returns>
    ''' <remarks></remarks>
    Public Function SearchResultAddModification(
      dblModificationMass As Double,
      chTargetResidue As Char,
      intResidueLocInPeptide As Integer,
      eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
      blnUpdateModOccurrenceCounts As Boolean) As Boolean

        Return SearchResultAddModification(dblModificationMass, chTargetResidue, intResidueLocInPeptide, eResidueTerminusState, blnUpdateModOccurrenceCounts, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)

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
    ''' <param name="modMassDigitsOfPrecision">Digits of precision to use when comparinig dblModificationMass to the masses of known mods</param>
    ''' <returns>True if mod successfully added</returns>
    ''' <remarks></remarks>
    Public Function SearchResultAddModification(
      dblModificationMass As Double,
      chTargetResidue As Char,
      intResidueLocInPeptide As Integer,
      eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
      blnUpdateModOccurrenceCounts As Boolean,
      modMassDigitsOfPrecision As Byte) As Boolean

        Dim objModificationDefinition As clsModificationDefinition
        Dim blnExistingModFound As Boolean
        Dim blnSuccess = False

        If intResidueLocInPeptide < 1 Then
            ' Invalid position; ignore this modification
            mErrorMessage = "Invalid value for intResidueLocInPeptide: " & intResidueLocInPeptide.ToString
        Else
            ' Lookup the modification definition given the modification information
            ' If the modification mass is unknown, then will auto-add it to the list of known modifications
            objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(
              dblModificationMass,
              chTargetResidue,
              eResidueTerminusState,
              blnExistingModFound,
              True, modMassDigitsOfPrecision)

            blnSuccess = SearchResultAddModification(
              objModificationDefinition,
              chTargetResidue,
              intResidueLocInPeptide,
              eResidueTerminusState,
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
    Public Function SearchResultAddModification(
       ByRef objModificationDefinition As clsModificationDefinition,
       chTargetResidue As Char,
       intResidueLocInPeptide As Integer,
       eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
       blnUpdateModOccurrenceCounts As Boolean) As Boolean

        Dim blnSuccess = False

        If intResidueLocInPeptide < 1 And objModificationDefinition.ModificationType <> clsModificationDefinition.eModificationTypeConstants.IsotopicMod Then
            ' Invalid position; ignore this modification
            mErrorMessage = "Invalid value for intResidueLocInPeptide: " & intResidueLocInPeptide.ToString & " (objModificationDefinition.ModificationType = " & objModificationDefinition.ModificationType.ToString & ")"
        Else
            If blnUpdateModOccurrenceCounts Then
                ' Increment OccurrenceCount
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
    Public Function SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts As Boolean) As Boolean
        Dim intIndex As Integer
        Dim blnSuccess = False

        For intIndex = 0 To mPeptideMods.ModificationCount - 1
            If mPeptideMods.GetModificationTypeByIndex(intIndex) = clsModificationDefinition.eModificationTypeConstants.IsotopicMod Then
                Dim intResidueLocInPeptide = 0

                blnSuccess = SearchResultAddModification(
                  mPeptideMods.GetModificationByIndex(intIndex),
                  clsPeptideMassCalculator.NO_AFFECTED_ATOM_SYMBOL,
                  intResidueLocInPeptide,
                  clsAminoAcidModInfo.eResidueTerminusStateConstants.None,
                  blnUpdateModOccurrenceCounts)

            End If
        Next intIndex

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Add any protein or peptide terminus static mods that are defined
    ''' </summary>
    ''' <param name="blnAllowDuplicateModOnTerminus">When false, only add the modification if the terminal residue does not already have the given modification associated with it</param>
    ''' <param name="blnUpdateModOccurrenceCounts"></param>
    ''' <remarks>
    ''' Peptide terminus static mods are always considered for a given peptide
    ''' Protein terminus static mods are only considered if the peptide is at the appropriate terminus for the modification
    ''' </remarks>
    Public Sub SearchResultAddStaticTerminusMods(blnAllowDuplicateModOnTerminus As Boolean, blnUpdateModOccurrenceCounts As Boolean)
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
                    If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus Or
                       mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
                    Else
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
                    End If
                    blnAddModification = True
                ElseIf objModificationDefinition.TargetResidues = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                    intResidueLocInPeptide = mPeptideCleanSequence.Length
                    If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus Or
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
                    If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNTerminus Or
                       mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinNandCCTerminus Then
                        intResidueLocInPeptide = 1
                        eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
                        blnAddModification = True
                    End If
                ElseIf objModificationDefinition.TargetResidues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS Then
                    If mPeptideTerminusState = ePeptideTerminusStateConstants.ProteinCTerminus Or
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
                                If Math.Abs(Math.Round(dblMassDifference, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)) < Single.Epsilon Then
                                    blnAddModification = False
                                    Exit For
                                End If
                            End If

                        End If
                    Next intIndexCompare
                End If

                If blnAddModification Then
                    SearchResultAddModification(
                      objModificationDefinition,
                      mPeptideCleanSequence.Chars(intResidueLocInPeptide - 1),
                      intResidueLocInPeptide,
                      eResidueTerminusState,
                      blnUpdateModOccurrenceCounts)
                End If
            End If
        Next intModificationIndex

    End Sub

    ''' <summary>
    ''' Updates the N-Terminal mass applied to peptides when computing their mass if it is significantly different than the currently defined N-terminal peptide mass
    ''' </summary>
    ''' <param name="dblNTerminalMassChange"></param>
    ''' <remarks></remarks>
    Public Sub UpdatePeptideNTerminusMass(dblNTerminalMassChange As Double)
        If Math.Round(Math.Abs(dblNTerminalMassChange - mPeptideSeqMassCalculator.PeptideNTerminusMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) > 0 Then
            mPeptideSeqMassCalculator.PeptideNTerminusMass = dblNTerminalMassChange
        End If
    End Sub

    ''' <summary>
    ''' Updates the C-Terminal mass applied to peptides when computing their mass if significantly different than the currently defined C-terminal peptide mass
    ''' </summary>
    ''' <param name="dblCTerminalMassChange"></param>
    ''' <remarks></remarks>
    Public Sub UpdatePeptideCTerminusMass(dblCTerminalMassChange As Double)
        If Math.Round(Math.Abs(dblCTerminalMassChange - mPeptideSeqMassCalculator.PeptideCTerminusMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) > 0 Then
            mPeptideSeqMassCalculator.PeptideCTerminusMass = dblCTerminalMassChange
        End If
    End Sub

    Public Sub UpdateSearchResultEnzymeAndTerminusInfo(udtEnzymeMatchSpec As udtEnzymeMatchSpecType, dblPeptideNTerminusMassChange As Double, dblPeptideCTerminusMassChange As Double)

        SetEnzymeMatchSpec(udtEnzymeMatchSpec)

        ' Update the N-Terminus and/or C-Terminus masses if those in the XML file are significantly different than the defaults
        If Math.Abs(dblPeptideNTerminusMassChange - 0) > Single.Epsilon Then
            UpdatePeptideNTerminusMass(dblPeptideNTerminusMassChange)
        End If

        If Math.Abs(dblPeptideCTerminusMassChange - 0) > Single.Epsilon Then
            UpdatePeptideCTerminusMass(dblPeptideCTerminusMassChange)
        End If

    End Sub

    ''' <summary>
    ''' Obtain the peptide sequence
    ''' </summary>
    ''' <param name="blnReturnSequenceWithMods">When true, then include the mod symbols in the sequence</param>
    ''' <returns></returns>
    ''' <remarks>
    ''' If you want to guarantee that mod symbols are included in the peptide sequence, 
    ''' You must call ApplyModificationInformation before using this function
    ''' </remarks>
    Public Function SequenceWithPrefixAndSuffix(blnReturnSequenceWithMods As Boolean) As String

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

        If blnReturnSequenceWithMods AndAlso Not PeptideSequenceWithMods Is Nothing Then
            Return chPrefix & "." & PeptideSequenceWithMods & "." & chSuffix
        Else
            If mPeptideCleanSequence Is Nothing Then
                Return String.Empty
            Else
                Return chPrefix & "." & mPeptideCleanSequence & "." & chSuffix
            End If
        End If
    End Function

    ''' <summary>
    ''' Define custom RegEx specs used to find enzyme cleavage sites
    ''' </summary>
    ''' <param name="udtEnzymeMatchSpec"></param>
    ''' <remarks>Define standard RegEx values using SetStandardEnzymeMatchSpec</remarks>
    Public Sub SetEnzymeMatchSpec(udtEnzymeMatchSpec As udtEnzymeMatchSpecType)
        mPeptideCleavageStateCalculator.SetEnzymeMatchSpec(udtEnzymeMatchSpec.LeftResidueRegEx, udtEnzymeMatchSpec.RightResidueRegEx)
    End Sub

    ''' <summary>
    ''' Stores strSequenceWithMods in PeptideSequenceWithMods
    ''' </summary>
    ''' <param name="strSequenceWithMods"></param>
    ''' <param name="blnCheckForPrefixAndSuffixResidues"></param>
    ''' <param name="blnAutoPopulateCleanSequence">When true, populates PeptideCleanSequence, which automatically calls ComputePeptideCleavageStateInProtein</param>
    ''' <remarks></remarks>
    Public Sub SetPeptideSequenceWithMods(strSequenceWithMods As String, blnCheckForPrefixAndSuffixResidues As Boolean, blnAutoPopulateCleanSequence As Boolean)

        Dim strPrimarySequence As String = String.Empty       ' Sequence with mods, but without the prefix or suffix residues
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

        PeptideSequenceWithMods = String.Copy(strPrimarySequence)

    End Sub

    ''' <summary>
    ''' Define standard RegEx values for finding enzyming cleavage sites
    ''' </summary>
    ''' <param name="eStandardCleavageAgent"></param>
    ''' <remarks>Define custom RegEx values using SetEnzymeMatchSpec</remarks>
    Public Sub SetStandardEnzymeMatchSpec(eStandardCleavageAgent As eStandardCleavageAgentConstants)
        mPeptideCleavageStateCalculator.SetStandardEnzymeMatchSpec(eStandardCleavageAgent)
    End Sub

    ''' <summary>
    ''' Generate a comma separated list of the modifications present and the residue modified, with a colon separating the mod name and residue location
    ''' Accessible via PeptideModDescription
    ''' </summary>
    ''' <remarks>
    ''' Examples:
    '''   Acetyl:1
    '''   MinusH2O:1,Plus1Oxy:4,Plus1Oxy:19
    ''' </remarks>
    Public Sub UpdateModDescription()

        Const MOD_LIST_SEP_CHAR = ","c
        Dim intIndex As Integer

        Dim udtModNameAndResidueLoc() As clsPHRPBaseClass.udtModNameAndResidueLocType
        Dim intPointerArray() As Integer

        PeptideModDescription = String.Empty

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
                    If intIndex > 0 Then PeptideModDescription &= MOD_LIST_SEP_CHAR
                    PeptideModDescription &= .ModDefinition.MassCorrectionTag.Trim & ":"c & .ResidueLocInPeptide
                End With
            Next intIndex

        End If

    End Sub

    Protected Class IGenericResidueModificationInfoComparer
        Implements IComparer(Of clsAminoAcidModInfo)

        ''' <summary>
        ''' Comparer
        ''' </summary>
        ''' <param name="x"></param>
        ''' <param name="y"></param>
        ''' <returns></returns>
        ''' <remarks></remarks>
        Public Function Compare(x As clsAminoAcidModInfo, y As clsAminoAcidModInfo) As Integer Implements IComparer(Of clsAminoAcidModInfo).Compare
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
