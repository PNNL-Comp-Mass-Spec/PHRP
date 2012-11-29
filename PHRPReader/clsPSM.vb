'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class tracks the details for a peptide hit search result (typically loaded from a tab-delimited text file created by the Peptide File Extractor or by PHRP)
'
'*********************************************************************************************************

Option Strict On

Public Class clsPSM

	Public Const UNKNOWN_COLLISION_MODE As String = "n/a"

	' Note: Be sure to update the Clone() function if you add new class-wide variables
	Protected mResultID As Integer
	Protected mScoreRank As Integer					' Top scoring peptide is rank 1, next lowest score is rank 2, etc.

	Protected mScanNumber As Integer
	Protected mElutionTimeMinutes As Single
	Protected mScanList As System.Collections.Generic.SortedSet(Of Integer)			' List of scans that were combined prior to identifying this peptide

	Protected mPeptide As String					' Peptide Sequence, with or without prefix & suffix residues; may contain mod symbols; example: R.RM*VNSGSGADSAVDLNSIPVAMIAR.V
	Protected mPeptideWithNumericMods As String		' Peptide Sequence where modified residues have the modification mass indicated as a number, example: R.N+144.102063SNPVIAELSQAINSGTLLSK+144.102063PS+79.9663PPLPPK+144.102063.R
	Protected mPeptideCleanSequence As String
	Protected mSeqID As Integer

	Protected mModifiedPeptideResidues As System.Collections.Generic.List(Of clsAminoAcidModInfo)

	Protected mCharge As Short
	Protected mCollisionMode As String		' CID, ETD, HCD, or n/a
	Protected mMSGFSpecProb As String		' MSGF SpecProb; stored as string to preserve formatting

	' The following are typically populated using UpdateCleavageInfo
	Protected mCleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants
	Protected mNumMissedCleavages As Short
	Protected mNumTrypticTerminii As Short

	Protected mPrecursorNeutralMass As Double		' Uncharged monoisotopic mass of the observed precursor
	Protected mMassErrorDa As String
	Protected mMassErrorPPM As String

	Protected mPeptideMonoisotopicMass As Double	' Theoretical (computed) monoisotopic mass

	' Note that protein names are case-sensitive
	Protected mProteins As System.Collections.Generic.List(Of String)

	' This dictionary tracks additional, tool-specific scores
	Protected mAdditionalScores As System.Collections.Generic.Dictionary(Of String, String)

#Region "Properties"

	''' <summary>
	''' Returns a dictionary with additional search engine scores stored as key/value pairs
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>Update scores using SetScore</remarks>
	Public ReadOnly Property AdditionalScores() As System.Collections.Generic.Dictionary(Of String, String)
		Get
			Return mAdditionalScores
		End Get
	End Property

	''' <summary>
	''' Assumed charge of the spectrum in which this peptide was identified
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property Charge As Short
		Get
			Return mCharge
		End Get
		Set(value As Short)
			mCharge = value
		End Set
	End Property

	''' <summary>
	''' Peptide cleavage state (with regards to ProteinFirst)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property CleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants
		Get
			Return mCleavageState
		End Get
		Set(value As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants)
			mCleavageState = value
		End Set
	End Property

	''' <summary>
	''' Collision mode (CID, ETD, HCD)
	''' PepXML allows this to be CID, ETD, ECD, ETD/CID, or HCD
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property CollisionMode As String
		Get
			Return mCollisionMode
		End Get
		Set(value As String)
			mCollisionMode = value
		End Set
	End Property

	''' <summary>
	''' Elution time (in minutes) of the spectrum
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property ElutionTimeMinutes As Single
		Get
			Return mElutionTimeMinutes
		End Get
		Set(value As Single)
			mElutionTimeMinutes = value
		End Set
	End Property

	''' <summary>
	''' Mass difference, in daltons, between the monoisotopic mass of the precursor ion and the calculated (theoretical) monoisotopic mass of the peptide
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property MassErrorDa As String
		Get
			Return mMassErrorDa
		End Get
		Set(value As String)
			mMassErrorDa = value
		End Set
	End Property

	''' <summary>
	''' Mass difference, in ppm, between the monoisotopic mass of the precursor ion and the calculated (theoretical) monoisotopic mass of the peptide
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property MassErrorPPM As String
		Get
			Return mMassErrorPPM
		End Get
		Set(value As String)
			mMassErrorPPM = value
		End Set
	End Property

	''' <summary>
	''' List of modified residues
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>A given residue is allowed to have more than one modification</remarks>
	Public ReadOnly Property ModifiedResidues As System.Collections.Generic.List(Of clsAminoAcidModInfo)
		Get
			Return mModifiedPeptideResidues
		End Get
	End Property

	''' <summary>
	''' MSGF Spectral Probability associated with this peptide
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>Ranges from 0 to 1, where 0 is the best score and 1 is the worse score</remarks>
	Public Property MSGFSpecProb As String
		Get
			Return mMSGFSpecProb
		End Get
		Set(value As String)
			mMSGFSpecProb = value
		End Set
	End Property

	''' <summary>
	''' Number of missed cleavages (internal K or R)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property NumMissedCleavages As Short
		Get
			Return mNumMissedCleavages
		End Get
		Set(value As Short)
			mNumMissedCleavages = value
		End Set
	End Property

	''' <summary>
	''' Number of tryptic terminii (or similar if not using trypsin)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>2 means fully tryptic, 1 means partially tryptic, 0 means non-tryptic</remarks>
	Public Property NumTrypticTerminii As Short
		Get
			Return mNumTrypticTerminii
		End Get
		Set(value As Short)
			mNumTrypticTerminii = value
		End Set
	End Property

	''' <summary>
	''' Peptide sequence, including any modification symbols that were assigned by the search engine
	''' For example, R.AAS*PQDLAGGYTSSLACHR.A
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property Peptide() As String
		Get
			Return mPeptide
		End Get
	End Property

	''' <summary>
	''' Peptide residues without any modification symbols or flanking residues
	''' For example, AASPQDLAGGYTSSLACHR
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property PeptideCleanSequence() As String
		Get
			Return mPeptideCleanSequence
		End Get
	End Property

	''' <summary>
	''' Computed monoisotopic mass (uncharged, theoretical mass, including mods)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property PeptideMonoisotopicMass As Double
		Get
			Return mPeptideMonoisotopicMass
		End Get
		Set(value As Double)
			mPeptideMonoisotopicMass = value
		End Set
	End Property

	''' <summary>
	''' Peptide sequence where all modified residues have the modification masses displayed as numeric values
	''' For example, R.A+144.102063AS+79.9663PQDLAGGYTSSLAC+57.0215HR.A
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property PeptideWithNumericMods() As String
		Get
			Return mPeptideWithNumericMods
		End Get
		Set(value As String)
			If String.IsNullOrEmpty(value) Then
				mPeptideWithNumericMods = String.Empty
			Else
				mPeptideWithNumericMods = value
			End If
		End Set
	End Property

	''' <summary>
	''' First protein associated with this peptide
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>Retrieve full list of proteins using the Proteins property</remarks>
	Public ReadOnly Property ProteinFirst() As String
		Get
			If mProteins.Count = 0 Then
				Return String.Empty
			Else
				Return mProteins(0)
			End If
		End Get
	End Property

	''' <summary>
	''' Uncharged monoisotopic mass of the precursor (observed mass based on m/z and charge)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property PrecursorNeutralMass As Double
		Get
			Return mPrecursorNeutralMass
		End Get
		Set(value As Double)
			mPrecursorNeutralMass = value
		End Set
	End Property

	''' <summary>
	''' List of proteins associated with this peptide
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property Proteins() As System.Collections.Generic.List(Of String)
		Get
			Return mProteins
		End Get
	End Property

	''' <summary>
	''' ResultID of this peptide (typically assigned by the search engine)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property ResultID As Integer
		Get
			Return mResultID
		End Get
		Set(value As Integer)
			mResultID = value
		End Set
	End Property

	''' <summary>
	''' List of scans that were combined prior to identifying this peptide
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property ScanList() As System.Collections.Generic.SortedSet(Of Integer)
		Get
			Return mScanList
		End Get
	End Property

	''' <summary>
	''' Scan number of the mass spectrum in which this peptide was identified
	''' Will automatically update ScanList if it does not yet contain this scan number
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property ScanNumber As Integer
		Get
			Return mScanNumber
		End Get
		Set(value As Integer)
			mScanNumber = value
			If Not mScanList.Contains(value) Then
				mScanList.Add(value)
			End If
		End Set
	End Property

	Public ReadOnly Property ScanNumberStart As Integer
		Get
			Return mScanList.Min()
		End Get
	End Property

	Public ReadOnly Property ScanNumberEnd As Integer
		Get
			Return mScanList.Max()
		End Get
	End Property

	''' <summary>
	''' Rank of this peptide in the given spectrum
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>Top scoring peptide is rank 1, next lowest score is rank 2, etc.</remarks>
	Public Property ScoreRank As Integer
		Get
			Return mScoreRank
		End Get
		Set(value As Integer)
			mScoreRank = value
		End Set
	End Property

	''' <summary>
	''' Sequence ID value assigned by PHRP
	''' Required for looking up information from the SeqInfo files
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Property SeqID As Integer
		Get
			Return mSeqID
		End Get
		Set(value As Integer)
			mSeqID = value
		End Set
	End Property

#End Region

	''' <summary>
	''' Constructor; auto-calls Clear()
	''' </summary>
	''' <remarks></remarks>
	Public Sub New()
		mScanList = New System.Collections.Generic.SortedSet(Of Integer)
		mProteins = New System.Collections.Generic.List(Of String)
		mModifiedPeptideResidues = New System.Collections.Generic.List(Of clsAminoAcidModInfo)
		mAdditionalScores = New System.Collections.Generic.Dictionary(Of String, String)(StringComparer.CurrentCultureIgnoreCase)
		Me.Clear()
	End Sub

	Public Sub AddCombinedScan(intScanNumber As Integer)
		If Not mScanList.Contains(intScanNumber) Then
			mScanList.Add(intScanNumber)
		End If
	End Sub
	''' <summary>
	''' Add the details for a modified residue
	''' </summary>
	''' <param name="objModInfo">Modification info class</param>
	''' <remarks></remarks>
	Public Sub AddModifiedResidue(objModInfo As clsAminoAcidModInfo)
		mModifiedPeptideResidues.Add(objModInfo)
	End Sub

	''' <summary>
	''' Add the details for a modified residue
	''' </summary>
	''' <param name="Residue">Amino acid letter; use angle brackets or square brackes for peptide or protein terminii (see the SYMBOL_DMS constants in clsAminoAcidModInfo)</param>
	''' <param name="ResidueLocInPeptide">Location of the residue in the peptide; use 1 for an N-terminal mod</param>
	''' <param name="ResidueTerminusState">Terminus state of residue</param>
	''' <param name="ModDefinition">Modification details</param>
	''' <remarks></remarks>
	Public Sub AddModifiedResidue(Residue As Char, ResidueLocInPeptide As Integer, ResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants, ModDefinition As clsModificationDefinition)
		mModifiedPeptideResidues.Add(New clsAminoAcidModInfo(Residue, ResidueLocInPeptide, ResidueTerminusState, ModDefinition))
	End Sub

	''' <summary>
	''' Add a new protein to associate with this peptide
	''' </summary>
	''' <param name="strProteinName"></param>
	''' <remarks></remarks>
	Public Sub AddProtein(strProteinName As String)
		If Not String.IsNullOrWhiteSpace(strProteinName) AndAlso Not mProteins.Contains(strProteinName) Then
			mProteins.Add(strProteinName)
		End If
	End Sub

	''' <summary>
	''' Reset the peptide to default values (and empty strings)
	''' </summary>
	''' <remarks></remarks>
	Public Sub Clear()
		mScanNumber = 0
		mElutionTimeMinutes = 0

		mScanList.Clear()

		mPeptide = String.Empty
		mPeptideWithNumericMods = String.Empty
		mPeptideCleanSequence = String.Empty
		mCharge = 0
		mResultID = 0
		mScoreRank = 0

		mCollisionMode = UNKNOWN_COLLISION_MODE
		mMSGFSpecProb = String.Empty

		mCleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific
		mNumMissedCleavages = 0
		mNumTrypticTerminii = 0

		mPrecursorNeutralMass = 0
		mMassErrorDa = String.Empty
		mMassErrorPPM = String.Empty

		mPeptideMonoisotopicMass = 0

		mProteins.Clear()
		mModifiedPeptideResidues.Clear()
		mAdditionalScores.Clear()
	End Sub

	Public Sub ClearModifiedResidues()
		mModifiedPeptideResidues.Clear()
	End Sub

	''' <summary>
	''' Duplicate this PSM object and return a new one
	''' </summary>
	''' <returns></returns>
	''' <remarks></remarks>
	Public Function Clone() As clsPSM
		Dim objNew As clsPSM

		objNew = New clsPSM()

		With objNew
			.ResultID = mResultID
			.ScoreRank = mScoreRank
			.ScanNumber = mScanNumber
			.ElutionTimeMinutes = mElutionTimeMinutes

			For Each intScanNumber In mScanList
				.AddCombinedScan(intScanNumber)
			Next

			.SetPeptide(mPeptide)				' Note: this will auto-update mPeptideCleanSequence in objNew
			.PeptideWithNumericMods = mPeptideWithNumericMods
			.Charge = mCharge
			.CollisionMode = mCollisionMode
			.MSGFSpecProb = mMSGFSpecProb

			.CleavageState = mCleavageState
			.NumMissedCleavages = mNumMissedCleavages
			.NumTrypticTerminii = mNumTrypticTerminii

			.PrecursorNeutralMass = mPrecursorNeutralMass
			.MassErrorDa = mMassErrorDa
			.MassErrorPPM = mMassErrorPPM

			For Each strProtein In mProteins
				.AddProtein(strProtein)
			Next

			For Each objItem In mModifiedPeptideResidues
				.AddModifiedResidue(objItem.Residue, objItem.ResidueLocInPeptide, objItem.ResidueTerminusState, objItem.ModDefinition)
			Next

			For Each objScore As System.Collections.Generic.KeyValuePair(Of String, String) In mAdditionalScores
				.SetScore(objScore.Key, objScore.Value)
			Next

		End With

		Return objNew
	End Function

	Protected Sub UpdateCleanSequence(strPeptide As String)
		If String.IsNullOrEmpty(strPeptide) Then
			mPeptideCleanSequence = String.Empty
		Else
			mPeptideCleanSequence = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strPeptide, True)
		End If
	End Sub

	''' <summary>
	''' Returns the value stored for the specified score
	''' </summary>
	''' <param name="strScoreName">Score name</param>
	''' <returns>Score if defined, otherwise an empty string</returns>
	Public Function GetScore(ByVal strScoreName As String) As String
		Dim strScoreValue As String = String.Empty
		If mAdditionalScores.TryGetValue(strScoreName, strScoreValue) Then
			Return strScoreValue
		Else
			Return String.Empty
		End If
	End Function

	''' <summary>
	'''  Returns the value stored for the specified score (as a double)
	''' </summary>
	''' <param name="strScoreName">Score name</param>
	''' <returns>Score if defined, otherwise 0</returns>
	''' <remarks></remarks>
	Public Function GetScoreDbl(ByVal strScoreName As String) As Double
		Return GetScoreDbl(strScoreName, 0)
	End Function

	''' <summary>
	'''  Returns the value stored for the specified score (as a double)
	''' </summary>
	''' <param name="strScoreName">Score name</param>
	''' <param name="dblValueIfMissing">Value to return if the score is not defined</param>
	''' <returns>Score if defined, otherwise dblValueIfMissing</returns>
	''' <remarks></remarks>
	Public Function GetScoreDbl(ByVal strScoreName As String, ByVal dblValueIfMissing As Double) As Double
		Dim strScoreValue As String
		Dim dblScore As Double

		strScoreValue = GetScore(strScoreName)
		If Not String.IsNullOrEmpty(strScoreValue) AndAlso Double.TryParse(strScoreValue, dblScore) Then
			Return dblScore
		Else
			Return dblValueIfMissing
		End If

	End Function

	''' <summary>
	'''  Returns the value stored for the specified score (as an integer)
	''' </summary>
	''' <param name="strScoreName">Score name</param>
	''' <returns>Score if defined, otherwise 0</returns>
	''' <remarks></remarks>
	Public Function GetScoreInt(ByVal strScoreName As String) As Integer
		Return GetScoreInt(strScoreName, 0)
	End Function

	''' <summary>
	'''  Returns the value stored for the specified score (as an integer)
	''' </summary>
	''' <param name="strScoreName">Score name</param>
	''' <param name="intValueIfMissing">Value to return if the score is not defined</param>
	''' <returns>Score if defined, otherwise intValueIfMissing</returns>
	''' <remarks></remarks>
	Public Function GetScoreInt(ByVal strScoreName As String, ByVal intValueIfMissing As Integer) As Integer
		Dim strScoreValue As String
		Dim intScore As Integer

		strScoreValue = GetScore(strScoreName)
		If Not String.IsNullOrEmpty(strScoreValue) AndAlso Integer.TryParse(strScoreValue, intScore) Then
			Return intScore
		Else
			Return intValueIfMissing
		End If

	End Function

	''' <summary>
	''' Update the peptide sequence (auto-determines the clean sequence)
	''' </summary>
	''' <param name="strPeptide">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
	''' <remarks>Does not update the cleavage state info</remarks>
	Public Sub SetPeptide(ByVal strPeptide As String)
		If String.IsNullOrEmpty(strPeptide) Then
			mPeptide = strPeptide
		Else
			mPeptide = strPeptide
		End If
		UpdateCleanSequence(mPeptide)
	End Sub

	''' <summary>
	''' Update the peptide sequence (auto-determines the clean sequence); also auto-update the the cleavage state info
	''' </summary>
	''' <param name="strPeptide">Peptide sequence (can optionally contain modification symbols; can optionally contain prefix and suffix residues)</param>
	''' <param name="objCleavageStateCalculator">Cleavage state calculator object</param>
	''' <remarks></remarks>
	Public Sub SetPeptide(ByVal strPeptide As String, objCleavageStateCalculator As clsPeptideCleavageStateCalculator)
		SetPeptide(strPeptide)
		UpdateCleavageInfo(objCleavageStateCalculator)
	End Sub

	''' <summary>
	''' Add/update an additional score to associate with this peptide
	''' </summary>
	''' <param name="strScoreName"></param>
	''' <param name="strScoreValue"></param>
	''' <remarks></remarks>
	Public Sub SetScore(strScoreName As String, strScoreValue As String)
		If mAdditionalScores.ContainsKey(strScoreName) Then
			mAdditionalScores(strScoreName) = strScoreValue
		Else
			mAdditionalScores.Add(strScoreName, strScoreValue)
		End If
	End Sub
	
	''' <summary>
	''' Returns the value stored for the specified score
	''' </summary>
	''' <param name="strScoreName"></param>
	''' <param name="strScoreValue"></param>
	''' <returns>True if the score is defined, otherwise false</returns>
	''' <remarks></remarks>
	Public Function TryGetScore(ByVal strScoreName As String, ByRef strScoreValue As String) As Boolean

		strScoreValue = String.Empty
		If mAdditionalScores.TryGetValue(strScoreName, strScoreValue) Then
			Return True
		Else
			Return False
		End If

	End Function

	''' <summary>
	''' Auto-determine the number of missed cleavages, cleavage state, and number of tryptic terminii based on the peptide sequence
	''' </summary>
	''' <param name="objCleavageStateCalculator"></param>
	''' <remarks></remarks>
	Public Sub UpdateCleavageInfo(objCleavageStateCalculator As clsPeptideCleavageStateCalculator)

		mNumMissedCleavages = objCleavageStateCalculator.ComputeNumberOfMissedCleavages(mPeptide)

		mCleavageState = objCleavageStateCalculator.ComputeCleavageState(mPeptide)

		If mCleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Full Then
			mNumTrypticTerminii = 2
		ElseIf mCleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.Partial Then
			mNumTrypticTerminii = 1
		Else
			mNumTrypticTerminii = 0
		End If

	End Sub

End Class
