Option Strict On

Public Class clsSearchEngineParameters

    Public Const MASS_TYPE_MONOISOTOPIC As String = "monoisotopic"
    Public Const MASS_TYPE_AVERAGE As String = "average"

    Private mSearchEngineName As String
    Private mSearchEngineVersion As String
    Private mSearchDate As DateTime

    Private mPrecursorMassType As String            ' Typically "monoisotopic" or "average"
    Private mFragmentMassType As String

    Private mSearchEngineParamFilePath As String

    Private mEnzyme As String

    Private ReadOnly mModInfo As List(Of clsModificationDefinition)

    Private ReadOnly mParameters As Dictionary(Of String, String)

#Region "Properties"
    ''' <summary>
    ''' Enzyme name
    ''' </summary>
    ''' <returns></returns>
    Public Property Enzyme As String
        Get
            Return mEnzyme
        End Get
        Set
            If Value = "" Then Value = "none"
            mEnzyme = Value
        End Set
    End Property

    ''' <summary>
    ''' FASTA file path
    ''' </summary>
    ''' <returns></returns>
    Public Property FastaFilePath As String

    ''' <summary>
    ''' Fragment mass type
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks>Typically "monoisotopic" or "average"</remarks>
    Public Property FragmentMassType As String
        Get
            Return mFragmentMassType
        End Get
        Set
            If String.IsNullOrEmpty(Value) Then Value = MASS_TYPE_MONOISOTOPIC
            mFragmentMassType = Value
        End Set
    End Property

    ''' <summary>
    ''' Maximum number of internal cleavages (missed cleavage points)
    ''' </summary>
    ''' <returns></returns>
    Public Property MaxNumberInternalCleavages As Integer

    ''' <summary>
    ''' 0 means no-enzyme, 1 means partially tryptic, 2 means fully tryptic
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks>For trypsin, this is NTT or Number of Tryptic Terminii</remarks>
    Public Property MinNumberTermini As Integer

    ''' <summary>
    ''' Dynamic and static mods to search for
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property ModInfo As List(Of clsModificationDefinition)
        Get
            Return mModInfo
        End Get
    End Property

    ''' <summary>
    ''' Parameter dictionary (key/value pairs)
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property Parameters As Dictionary(Of String, String)
        Get
            Return mParameters
        End Get
    End Property

    ''' <summary>
    ''' Precursor mass tolerance, in Da; 0 if unknown
    ''' </summary>
    ''' <returns></returns>
    Public Property PrecursorMassToleranceDa As Double

    ''' <summary>
    ''' Precursor mass tolerance, in ppm; 0 if unknown
    ''' </summary>
    ''' <returns></returns>
    Public Property PrecursorMassTolerancePpm As Double

    ''' <summary>
    ''' Precursor mass type
    ''' </summary>
    ''' <returns></returns>
    ''' <remarks>Typically "monoisotopic" or "average"</remarks>
    Public Property PrecursorMassType As String
        Get
            Return mPrecursorMassType
        End Get
        Set
            If String.IsNullOrWhiteSpace(Value) Then Value = MASS_TYPE_MONOISOTOPIC
            mPrecursorMassType = Value
        End Set
    End Property

    ''' <summary>
    ''' Search engine name
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property SearchEngineName As String
        Get
            Return mSearchEngineName
        End Get
    End Property

    ''' <summary>
    ''' Search engine parameter file path
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property SearchEngineParamFilePath As String
        Get
            If String.IsNullOrWhiteSpace(mSearchEngineParamFilePath) Then
                Return String.Empty
            Else
                Return mSearchEngineParamFilePath
            End If
        End Get
    End Property

    ''' <summary>
    ''' Search engine version
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property SearchEngineVersion As String
        Get
            Return mSearchEngineVersion
        End Get
    End Property

    ''' <summary>
    ''' Search date
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property SearchDate As DateTime
        Get
            Return mSearchDate
        End Get
    End Property

#End Region

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="searchEngineName"></param>
    Public Sub New(searchEngineName As String)
        Me.New(searchEngineName, New List(Of clsModificationDefinition), Nothing)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="searchEngineName"></param>
    ''' <param name="objModInfo"></param>
    Public Sub New(searchEngineName As String, objModInfo As List(Of clsModificationDefinition))
        Me.New(searchEngineName, objModInfo, Nothing)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="searchEngineName"></param>
    ''' <param name="objModInfo"></param>
    ''' <param name="Parameters"></param>
    Public Sub New(searchEngineName As String, objModInfo As List(Of clsModificationDefinition), Parameters As Dictionary(Of String, String))
        Me.InitializeDefaults()

        mSearchEngineName = searchEngineName
        mSearchEngineParamFilePath = String.Empty

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

    ''' <summary>
    ''' Add a new dynamic or static modification
    ''' </summary>
    ''' <param name="objModInfo"></param>
    Public Sub AddModification(objModInfo As clsModificationDefinition)
        mModInfo.Add(objModInfo)
    End Sub

    ''' <summary>
    ''' Add/update a parameter
    ''' </summary>
    ''' <param name="kvSetting"></param>
    Public Sub AddUpdateParameter(kvSetting As KeyValuePair(Of String, String))
        AddUpdateParameter(kvSetting.Key, kvSetting.Value)
    End Sub

    ''' <summary>
    ''' Add/update a parameter
    ''' </summary>
    ''' <param name="ParamName"></param>
    ''' <param name="ParamValue"></param>
    Public Sub AddUpdateParameter(ParamName As String, ParamValue As String)
        If mParameters.ContainsKey(ParamName) Then
            mParameters(ParamName) = ParamValue
        Else
            mParameters.Add(ParamName, ParamValue)
        End If
    End Sub

    ''' <summary>
    ''' Clear stored dynamic and static modifications
    ''' </summary>
    Public Sub ClearModifications()
        mModInfo.Clear()
    End Sub

    ''' <summary>
    ''' Clear stored key/value parameters
    ''' </summary>
    Public Sub ClearParameters()
        mParameters.Clear()
    End Sub

    Private Sub InitializeDefaults()
        mSearchEngineName = "Unknown"
        mSearchEngineVersion = "Unknown"
        mSearchDate = New DateTime(1980, 1, 1)

        FastaFilePath = String.Empty

        PrecursorMassToleranceDa = 0
        PrecursorMassTolerancePpm = 0

        mPrecursorMassType = MASS_TYPE_MONOISOTOPIC
        mFragmentMassType = MASS_TYPE_MONOISOTOPIC

        mEnzyme = "trypsin"
        MaxNumberInternalCleavages = 4
        MinNumberTermini = 0

    End Sub

    ''' <summary>
    ''' Update the search engine parameter file path
    ''' </summary>
    ''' <param name="paramFilePath"></param>
    Public Sub UpdateSearchEngineParamFilePath(paramFilePath As String)
        mSearchEngineParamFilePath = paramFilePath
    End Sub

    ''' <summary>
    ''' Update the search engine version
    ''' </summary>
    ''' <param name="strSearchEngineVersion"></param>
    Public Sub UpdateSearchEngineVersion(strSearchEngineVersion As String)
        mSearchEngineVersion = String.Copy(strSearchEngineVersion)
    End Sub

    ''' <summary>
    ''' Update the search date
    ''' </summary>
    ''' <param name="dtSearchDate"></param>
    Public Sub UpdateSearchDate(dtSearchDate As DateTime)
        mSearchDate = dtSearchDate
    End Sub

End Class

