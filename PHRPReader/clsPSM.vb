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

	' Note: Be sure to update the Clone() function if you add new class-wide variables
	Protected mResultID As Integer
	Protected mScanNumber As Integer
	Protected mPeptide As String					' Peptide Sequence, with or without prefix & suffix residues; may contain mod symbols; example: R.RM*VNSGSGADSAVDLNSIPVAMIAR.V
	Protected mPeptideWithNumericMods As String		' Peptide Sequence where modified residues have the modification mass indicated as a number, example: R.N+144.102063SNPVIAELSQAINSGTLLSK+144.102063PS+79.9663PPLPPK+144.102063.R
	Protected mPeptideCleanSequence As String

	Protected mCharge As Short
	Protected mCollisionMode As String		' CID, ETD, HCD, or n/a
	Protected mMSGFSpecProb As String		' MSGF SpecProb; stored as string to preserve formatting

	' The following are typically populated using UpdateCleavageInfo
	Protected mCleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants
	Protected mNumMissedCleavages As Short
	Protected mNumTrypticTerminii As Short

	Protected mPrecursorNeutralMass As Double		' Uncharged monoisotopic mass of the precursor
	Protected mMassErrorDa As String
	Protected mMassErrorPPM As String

	' Note that protein names are case-sensitive
	Protected mProteins As System.Collections.Generic.List(Of String)

	' This dictionary tracks additional, tool-specific scores
	Protected mAdditionalScores As System.Collections.Generic.Dictionary(Of String, String)

#Region "Properties"

	Public ReadOnly Property AdditionalScores() As System.Collections.Generic.Dictionary(Of String, String)
		Get
			Return mAdditionalScores
		End Get
	End Property

	Public Property Charge As Short
		Get
			Return mCharge
		End Get
		Set(value As Short)
			mCharge = value
		End Set
	End Property

	Public Property CleavageState As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants
		Get
			Return mCleavageState
		End Get
		Set(value As clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants)
			mCleavageState = value
		End Set
	End Property

	Public Property CollisionMode As String
		Get
			Return mCollisionMode
		End Get
		Set(value As String)
			mCollisionMode = value
		End Set
	End Property

	Public Property MassErrorDa As String
		Get
			Return mMassErrorDa
		End Get
		Set(value As String)
			mMassErrorDa = value
		End Set
	End Property

	Public Property MassErrorPPM As String
		Get
			Return mMassErrorPPM
		End Get
		Set(value As String)
			mMassErrorPPM = value
		End Set
	End Property

	Public Property MSGFSpecProb As String
		Get
			Return mMSGFSpecProb
		End Get
		Set(value As String)
			mMSGFSpecProb = value
		End Set
	End Property

	Public Property NumMissedCleavages As Short
		Get
			Return mNumMissedCleavages
		End Get
		Set(value As Short)
			mNumMissedCleavages = value
		End Set
	End Property

	Public Property NumTrypticTerminii As Short
		Get
			Return mNumTrypticTerminii
		End Get
		Set(value As Short)
			mNumTrypticTerminii = value
		End Set
	End Property

	Public Property Peptide() As String
		Get
			Return mPeptide
		End Get
		Set(value As String)
			mPeptide = value
			UpdateCleanSequence(mPeptide)
		End Set
	End Property

	Public ReadOnly Property PeptideCleanSequence() As String
		Get
			Return mPeptideCleanSequence
		End Get
	End Property

	Public Property PeptideWithNumericMods() As String
		Get
			Return mPeptideWithNumericMods
		End Get
		Set(value As String)
			mPeptideWithNumericMods = value
		End Set
	End Property

	Public ReadOnly Property ProteinFirst() As String
		Get
			If mProteins.Count = 0 Then
				Return String.Empty
			Else
				Return mProteins(0)
			End If
		End Get
	End Property

	Public Property PrecursorNeutralMass As Double
		Get
			Return mPrecursorNeutralMass
		End Get
		Set(value As Double)
			mPrecursorNeutralMass = value
		End Set
	End Property

	Public ReadOnly Property Proteins() As System.Collections.Generic.List(Of String)
		Get
			Return mProteins
		End Get
	End Property

	Public Property ResultID As Integer
		Get
			Return mResultID
		End Get
		Set(value As Integer)
			mResultID = value
		End Set
	End Property

	Public Property ScanNumber As Integer
		Get
			Return mScanNumber
		End Get
		Set(value As Integer)
			mScanNumber = value
		End Set
	End Property

#End Region

	Public Sub AddProtein(strProteinName As String)
		If Not String.IsNullOrWhiteSpace(strProteinName) AndAlso Not mProteins.Contains(strProteinName) Then
			mProteins.Add(strProteinName)
		End If
	End Sub

	Public Sub Clear()
		mScanNumber = 0
		mPeptide = String.Empty
		mPeptideWithNumericMods = String.Empty
		mPeptideCleanSequence = String.Empty

		mCharge = 0
		mResultID = 0
		mCollisionMode = "n/a"
		mMSGFSpecProb = String.Empty

		mCleavageState = clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific
		mNumMissedCleavages = 0
		mNumTrypticTerminii = 0

		mPrecursorNeutralMass = 0
		mMassErrorDa = String.Empty
		mMassErrorPPM = String.Empty

		mProteins.Clear()
		mAdditionalScores.Clear()
	End Sub

	Public Function Clone() As clsPSM
		Dim objNew As clsPSM

		objNew = New clsPSM()

		With objNew
			.ResultID = mResultID
			.ScanNumber = mScanNumber
			.Peptide = mPeptide				' Note: this will auto-update mPeptideCleanSequence in objNew
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

			For Each objScore As System.Collections.Generic.KeyValuePair(Of String, String) In mAdditionalScores
				.SetScore(objScore.Key, objScore.Value)
			Next

		End With

		Return objNew
	End Function

	Protected Sub UpdateCleanSequence(strPeptide As String)
		mPeptideCleanSequence = clsPeptideCleavageStateCalculator.ExtractCleanSequenceFromSequenceWithMods(strPeptide, True)
	End Sub
	
	Public Function GetScore(ByVal strScoreName As String) As String
		Dim strScoreValue As String = String.Empty
		If mAdditionalScores.TryGetValue(strScoreName, strScoreValue) Then
			Return strScoreValue
		Else
			Return String.Empty
		End If
	End Function

	Public Function GetScoreDbl(ByVal strScoreName As String) As Double
		Return GetScoreDbl(strScoreName, 0)
	End Function

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

	Public Function GetScoreInt(ByVal strScoreName As String) As Integer
		Return GetScoreInt(strScoreName, 0)
	End Function

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

	Public Sub SetScore(strScoreName As String, strScoreValue As String)
		If mAdditionalScores.ContainsKey(strScoreName) Then
			mAdditionalScores(strScoreName) = strScoreValue
		Else
			mAdditionalScores.Add(strScoreName, strScoreValue)
		End If
	End Sub

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

	Public Sub New()
		mProteins = New System.Collections.Generic.List(Of String)
		mAdditionalScores = New System.Collections.Generic.Dictionary(Of String, String)(StringComparer.CurrentCultureIgnoreCase)
		Me.Clear()
	End Sub
End Class
