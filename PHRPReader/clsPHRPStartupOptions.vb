﻿Option Strict On

Public Class clsPHRPStartupOptions

	Protected mMaxProteinsPerPSM As Integer

	Public Property LoadModsAndSeqInfo As Boolean
	Public Property LoadMSGFResults As Boolean
	Public Property LoadScanStatsData As Boolean

	''' <summary>
	''' Maximum number of proteins to associate with each PSM
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks>0 means to load all proteins</remarks>
	Public Property MaxProteinsPerPSM As Integer
		Get
			If mMaxProteinsPerPSM <= 0 Then
				Return Integer.MaxValue
			End If
			Return mMaxProteinsPerPSM
		End Get
		Set(value As Integer)
			mMaxProteinsPerPSM = value
		End Set
	End Property

	Public Sub New()
		LoadModsAndSeqInfo = True
		LoadMSGFResults = True
		LoadScanStatsData = False
		mMaxProteinsPerPSM = 0		' 0 means to load all proteins
	End Sub
End Class