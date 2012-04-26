Option Strict On

Public Class clsScanStatsExInfo

	Protected mScanNumber As Integer
	Protected mIonInjectionTime As Double
	Protected mScanEvent As Integer
	Protected mMasterIndex As Integer
	Protected mElapsedScanTime As Double
	Protected mChargeState As Integer
	Protected mMonoisotopicMZ As Double
	Protected mMS2IsolationWidth As Double
	Protected mFTAnalyzerSettings As String
	Protected mFTAnalyzerMessage As String
	Protected mFTResolution As Double
	Protected mConversionParameterB As Double
	Protected mConversionParameterC As Double
	Protected mConversionParameterD As Double
	Protected mConversionParameterE As Double
	Protected mCollisionMode As String
	Protected mScanFilterText As String
	Protected mSourceVoltage As Double
	Protected mSource_Current As Double

	Public ReadOnly Property ScanNumber As Integer
		Get
			Return mScanNumber
		End Get
	End Property

	Public Property IonInjectionTime As Double
		Get
			Return mIonInjectionTime
		End Get
		Set(value As Double)
			mIonInjectionTime = value
		End Set
	End Property

	Public Property ScanEvent As Integer
		Get
			Return mScanEvent
		End Get
		Set(value As Integer)
			mScanEvent = value
		End Set
	End Property

	Public Property MasterIndex As Integer
		Get
			Return mMasterIndex
		End Get
		Set(value As Integer)
			mMasterIndex = value
		End Set
	End Property

	Public Property ElapsedScanTime As Double
		Get
			Return mElapsedScanTime
		End Get
		Set(value As Double)
			mElapsedScanTime = value
		End Set
	End Property

	Public Property ChargeState As Integer
		Get
			Return mChargeState
		End Get
		Set(value As Integer)
			mChargeState = value
		End Set
	End Property

	Public Property MonoisotopicMZ As Double
		Get
			Return mMonoisotopicMZ
		End Get
		Set(value As Double)
			mMonoisotopicMZ = value
		End Set
	End Property

	Public Property MS2IsolationWidth As Double
		Get
			Return mMS2IsolationWidth
		End Get
		Set(value As Double)
			mMS2IsolationWidth = value
		End Set
	End Property

	Public Property FTAnalyzerSettings As String
		Get
			Return mFTAnalyzerSettings
		End Get
		Set(value As String)
			mFTAnalyzerSettings = value
		End Set
	End Property

	Public Property FTAnalyzerMessage As String
		Get
			Return mFTAnalyzerMessage
		End Get
		Set(value As String)
			mFTAnalyzerMessage = value
		End Set
	End Property

	Public Property FTResolution As Double
		Get
			Return mFTResolution
		End Get
		Set(value As Double)
			mFTResolution = value
		End Set
	End Property

	Public Property ConversionParameterB As Double
		Get
			Return mConversionParameterB
		End Get
		Set(value As Double)
			mConversionParameterB = value
		End Set
	End Property

	Public Property ConversionParameterC As Double
		Get
			Return mConversionParameterC
		End Get
		Set(value As Double)
			mConversionParameterC = value
		End Set
	End Property

	Public Property ConversionParameterD As Double
		Get
			Return mConversionParameterD
		End Get
		Set(value As Double)
			mConversionParameterD = value
		End Set
	End Property

	Public Property ConversionParameterE As Double
		Get
			Return mConversionParameterE
		End Get
		Set(value As Double)
			mConversionParameterE = value
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

	Public Property ScanFilterText As String
		Get
			Return mScanFilterText
		End Get
		Set(value As String)
			mScanFilterText = value
		End Set
	End Property

	Public Property SourceVoltage As Double
		Get
			Return mSourceVoltage
		End Get
		Set(value As Double)
			mSourceVoltage = value
		End Set
	End Property

	Public Property Source_Current As Double
		Get
			Return mSource_Current
		End Get
		Set(value As Double)
			mSource_Current = value
		End Set
	End Property

	Public Sub New(ByVal intScanNumber As Integer)
		mScanNumber = intScanNumber
	End Sub
End Class
