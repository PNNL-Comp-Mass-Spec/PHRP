Option Strict On

Public Class clsSearchEngineParameters

    Public Const MASS_TYPE_MONOISOTOPIC As String = "monoisotopic"
    Public Const MASS_TYPE_AVERAGE As String = "average"

    Private mSearchEngineName As String
    Private mSearchEngineVersion As String
    Private mSearchDate As DateTime

    Private mFastaFilePath As String

    Private mPrecursorMassToleranceDa As Double     ' Precursor mass tolerance, in Da; 0 if unknown
    Private mPrecursorMassTolerancePpm As Double    ' Precursor mass tolerance, in ppm; 0 if unknown

    Private mPrecursorMassType As String            ' Typically "monoisotopic" or "average"
    Private mFragmentMassType As String

    Private mEnzyme As String
    Private mMaxNumberInternalCleavages As Integer
    Private mMinNumberTermini As Integer                        ' 0 means no-enzyme, 1 means partially tryptic, 2 means fully tryptic

    Private ReadOnly mModInfo As List(Of clsModificationDefinition)

    Private ReadOnly mParameters As Dictionary(Of String, String)

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

    Public Property FastaFilePath As String
        Get
            Return mFastaFilePath
        End Get
        Set(value As String)
            mFastaFilePath = value
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

    Public ReadOnly Property ModInfo As List(Of clsModificationDefinition)
        Get
            Return mModInfo
        End Get
    End Property

    Public ReadOnly Property Parameters As Dictionary(Of String, String)
        Get
            Return mParameters
        End Get
    End Property

    Public Property PrecursorMassToleranceDa As Double
        Get
            Return mPrecursorMassToleranceDa
        End Get
        Set(value As Double)
            mPrecursorMassToleranceDa = value
        End Set
    End Property

    Public Property PrecursorMassTolerancePpm As Double
        Get
            Return mPrecursorMassTolerancePpm
        End Get
        Set(value As Double)
            mPrecursorMassTolerancePpm = value
        End Set
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

    Public ReadOnly Property SearchEngineVersion As String
        Get
            Return mSearchEngineVersion
        End Get
    End Property

    Public ReadOnly Property SearchDate As DateTime
        Get
            Return mSearchDate
        End Get
    End Property

#End Region

    Public Sub New(searchEngineName As String)
        Me.New(searchEngineName, New List(Of clsModificationDefinition), Nothing)
    End Sub

    Public Sub New(searchEngineName As String, objModInfo As List(Of clsModificationDefinition))
        Me.New(searchEngineName, objModInfo, Nothing)
    End Sub

    Public Sub New(searchEngineName As String, objModInfo As List(Of clsModificationDefinition), Parameters As Dictionary(Of String, String))
        Me.InitializeDefaults()

        mSearchEngineName = searchEngineName

        mModInfo = objModInfo

        If objModInfo Is Nothing Then
            mModInfo = New List(Of clsModificationDefinition)
        Else
            mModInfo = objModInfo
        End If

        If mParameters Is Nothing Then
            mParameters = New Dictionary(Of String, String)(StringComparer.CurrentCultureIgnoreCase)
        Else
            mParameters = Parameters
        End If

    End Sub

    Public Sub AddModification(objModInfo As clsModificationDefinition)
        mModInfo.Add(objModInfo)
    End Sub

    Public Sub AddUpdateParameter(kvSetting As KeyValuePair(Of String, String))
        AddUpdateParameter(kvSetting.Key, kvSetting.Value)
    End Sub

    Public Sub AddUpdateParameter(ParamName As String, ParamValue As String)
        If mParameters.ContainsKey(ParamName) Then
            mParameters(ParamName) = ParamValue
        Else
            mParameters.Add(ParamName, ParamValue)
        End If
    End Sub

    Public Sub ClearModifications()
        mModInfo.Clear()
    End Sub

    Public Sub ClearParameters()
        mParameters.Clear()
    End Sub

    Private Sub InitializeDefaults()
        mSearchEngineName = "Unknown"
        mSearchEngineVersion = "Unknown"
        mSearchDate = New DateTime(1980, 1, 1)

        mFastaFilePath = String.Empty

        mPrecursorMassToleranceDa = 0
        mPrecursorMassTolerancePpm = 0

        mPrecursorMassType = MASS_TYPE_MONOISOTOPIC
        mFragmentMassType = MASS_TYPE_MONOISOTOPIC

        mEnzyme = "trypsin"
        mMaxNumberInternalCleavages = 4
        mMinNumberTermini = 0

    End Sub

    Public Sub UpdateSearchEngineVersion(strSearchEngineVersion As String)
        mSearchEngineVersion = String.Copy(strSearchEngineVersion)
    End Sub

    Public Sub UpdateSearchDate(dtSearchDate As DateTime)
        mSearchDate = dtSearchDate
    End Sub

End Class

