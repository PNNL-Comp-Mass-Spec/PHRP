Option Strict On

Public Class clsScanStatsInfo

	Protected mScanNumber As Integer
	Protected mScanTimeMinutes As Single
	Protected mScanType As Integer
	Protected mTotalIonIntensity As Double
	Protected mBasePeakIntensity As Double
	Protected mBasePeakMZ As Double
	Protected mBasePeakSignalToNoiseRatio As Double
	Protected mIonCount As Integer
	Protected mIonCountRaw As Integer
	Protected mScanTypeName As String

	Public ReadOnly Property ScanNumber As Integer
		Get
			Return mScanNumber
		End Get
	End Property

	Public Property ScanTimeMinutes As Single
		Get
			Return mScanTimeMinutes
		End Get
		Set(value As Single)
			mScanTimeMinutes = value
		End Set
	End Property

	Public Property ScanType As Integer
		Get
			Return mScanType
		End Get
		Set(value As Integer)
			mScanType = value
		End Set
	End Property

	Public Property TotalIonIntensity As Double
		Get
			Return mTotalIonIntensity
		End Get
		Set(value As Double)
			mTotalIonIntensity = value
		End Set
	End Property

	Public Property BasePeakIntensity As Double
		Get
			Return mBasePeakIntensity
		End Get
		Set(value As Double)
			mBasePeakIntensity = value
		End Set
	End Property

	Public Property BasePeakMZ As Double
		Get
			Return mBasePeakMZ
		End Get
		Set(value As Double)
			mBasePeakMZ = value
		End Set
	End Property

	Public Property BasePeakSignalToNoiseRatio As Double
		Get
			Return mBasePeakSignalToNoiseRatio
		End Get
		Set(value As Double)
			mBasePeakSignalToNoiseRatio = value
		End Set
	End Property

	Public Property IonCount As Integer
		Get
			Return mIonCount
		End Get
		Set(value As Integer)
			mIonCount = value
		End Set
	End Property

	Public Property IonCountRaw As Integer
		Get
			Return mIonCountRaw
		End Get
		Set(value As Integer)
			mIonCountRaw = value
		End Set
	End Property

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
        mScanTimeMinutes = sngScanTimeMinutes
        mScanType = intScanType
    End Sub
End Class
