Option Strict On

Imports PeptideHitResultsProcessor

Public Class clsSearchEngineParameters

	Public Const MASS_TYPE_MONOISOTOPIC As String = "monoisotopic"
	Public Const MASS_TYPE_AVERAGE As String = "average"

	Protected mSearchEngineName As String
	Protected mPrecursorMassType As String
	Protected mFragmentMassType As String

	Protected mEnzyme As String
	Protected mMaxNumberInternalCleavages As Integer
	Protected mMinNumberTermini As Integer						' 0 means no-enzyme, 1 means partially tryptic, 2 means fully tryptic

	Protected mModInfo As System.Collections.Generic.List(Of clsModificationDefinition)

	Protected mParameters As System.Collections.Generic.Dictionary(Of String, String)

#Region "Properties"
	Public Property Enzyme As String
		Get
			Return mEnzyme
		End Get
		Set(value As String)
			If value = "" Then value = "none"
			mEnzyme = value
		End Set
	End Property


	Public Property FragmentMassType As String
		Get
			Return mFragmentMassType
		End Get
		Set(value As String)
			If String.IsNullOrEmpty(value) Then value = MASS_TYPE_MONOISOTOPIC
			mFragmentMassType = value
		End Set
	End Property

	Public Property MaxNumberInternalCleavages As Integer
		Get
			Return mMaxNumberInternalCleavages
		End Get
		Set(value As Integer)
			mMaxNumberInternalCleavages = value
		End Set
	End Property

	Public Property MinNumberTermini As Integer
		Get
			Return mMinNumberTermini
		End Get
		Set(value As Integer)
			mMinNumberTermini = value
		End Set
	End Property

	Public ReadOnly Property ModInfo As System.Collections.Generic.List(Of clsModificationDefinition)
		Get
			Return mModInfo
		End Get
	End Property

	Public ReadOnly Property Parameters As System.Collections.Generic.Dictionary(Of String, String)
		Get
			Return mParameters
		End Get
	End Property

	Public Property PrecursorMassType As String
		Get
			Return mPrecursorMassType
		End Get
		Set(value As String)
			If String.IsNullOrEmpty(value) Then value = MASS_TYPE_MONOISOTOPIC
			mPrecursorMassType = value
		End Set
	End Property

	Public ReadOnly Property SearchEngineName As String
		Get
			Return mSearchEngineName
		End Get
	End Property
#End Region

	Public Sub New(ByVal SearchEngineName As String)
		Me.New(SearchEngineName, New System.Collections.Generic.List(Of clsModificationDefinition), New System.Collections.Generic.Dictionary(Of String, String))
	End Sub

	Public Sub New(ByVal SearchEngineName As String, ByVal objModInfo As System.Collections.Generic.List(Of clsModificationDefinition))
		Me.New(SearchEngineName, objModInfo, New System.Collections.Generic.Dictionary(Of String, String))
	End Sub

	Public Sub New(ByVal SearchEngineName As String, ByVal objModInfo As System.Collections.Generic.List(Of clsModificationDefinition), ByVal Parameters As System.Collections.Generic.Dictionary(Of String, String))
		Me.InitializeDefaults()

		mSearchEngineName = SearchEngineName
		mModInfo = objModInfo


		If objModInfo Is Nothing Then
			mModInfo = New System.Collections.Generic.List(Of clsModificationDefinition)
		Else
			mModInfo = objModInfo
		End If

		If mParameters Is Nothing Then
			mParameters = New System.Collections.Generic.Dictionary(Of String, String)
		Else
			mParameters = Parameters
		End If

	End Sub

	Public Sub AddUpdateParameter(ByVal ParamName As String, ParamValue As String)
		If mParameters.ContainsKey(ParamName) Then
			mParameters(ParamName) = ParamValue
		Else
			mParameters.Add(ParamName, ParamValue)
		End If
	End Sub

	Public Sub ClearParameters()
		mParameters.Clear()
	End Sub

	Protected Sub InitializeDefaults()
		mSearchEngineName = "Unknown"
		mPrecursorMassType = MASS_TYPE_MONOISOTOPIC
		mFragmentMassType = MASS_TYPE_MONOISOTOPIC

		mEnzyme = "trypsin"
		mMaxNumberInternalCleavages = 3
		mMinNumberTermini = 1

	End Sub

End Class

