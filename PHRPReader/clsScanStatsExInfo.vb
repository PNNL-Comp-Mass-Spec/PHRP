Option Strict On

Public Class clsScanStatsExInfo

    Private ReadOnly mScanNumber As Integer

    ' ReSharper disable once ConvertToVbAutoProperty
    Public ReadOnly Property ScanNumber As Integer
        Get
            Return mScanNumber
        End Get
    End Property

    Public Property IonInjectionTime As Double

    Public Property ScanEvent As Integer

    Public Property MasterIndex As Integer

    Public Property ElapsedScanTime As Double

    Public Property ChargeState As Integer

    Public Property MonoisotopicMZ As Double

    Public Property MS2IsolationWidth As Double

    Public Property FTAnalyzerSettings As String

    Public Property FTAnalyzerMessage As String

    Public Property FTResolution As Double

    Public Property ConversionParameterB As Double

    Public Property ConversionParameterC As Double

    Public Property ConversionParameterD As Double

    Public Property ConversionParameterE As Double

    Public Property CollisionMode As String

    Public Property ScanFilterText As String

    Public Property SourceVoltage As Double

    Public Property Source_Current As Double

    Public Sub New(intScanNumber As Integer)
        mScanNumber = intScanNumber
    End Sub
End Class
