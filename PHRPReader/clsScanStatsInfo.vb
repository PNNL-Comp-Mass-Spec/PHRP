Option Strict On

Public Class clsScanStatsInfo

    Private ReadOnly mScanNumber As Integer
    Private mScanTypeName As String

    ' ReSharper disable once ConvertToVbAutoProperty
    Public ReadOnly Property ScanNumber As Integer
        Get
            Return mScanNumber
        End Get
    End Property

    Public Property ScanTimeMinutes As Single

    Public Property ScanType As Integer

    Public Property TotalIonIntensity As Double

    Public Property BasePeakIntensity As Double

    Public Property BasePeakMZ As Double

    Public Property BasePeakSignalToNoiseRatio As Double

    Public Property IonCount As Integer

    Public Property IonCountRaw As Integer

    Public Property ScanTypeName As String
        Get
            Return mScanTypeName
        End Get
        Set(value As String)
            If String.IsNullOrEmpty(value) Then
                mScanTypeName = String.Empty
            Else
                mScanTypeName = value
            End If
        End Set
    End Property

    Public Sub New(intScanNumber As Integer, sngScanTimeMinutes As Single, intScanType As Integer)
        mScanNumber = intScanNumber
        ScanTimeMinutes = sngScanTimeMinutes
        ScanType = intScanType
    End Sub
End Class
