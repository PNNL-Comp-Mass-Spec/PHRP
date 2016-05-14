'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/20/2012
'
' This class tracks the sequence information determined by PHRP and stored in a _SeqInfo.txt file
'
'*********************************************************************************************************

Option Strict On

Public Class clsSeqInfo
    Private ReadOnly mSeqID As Integer
    Private ReadOnly mModCount As Integer
    Private ReadOnly mModDescription As String
    Private mMonoisotopicMass As Double

	''' <summary>
	''' Sequence ID
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property SeqID As Integer
		Get
			Return mSeqID
		End Get
	End Property

	''' <summary>
	''' Number of modifications
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property ModCount As Integer
		Get
			Return mModCount
		End Get
	End Property

	''' <summary>
	''' Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property ModDescription As String
		Get
			Return mModDescription
		End Get
	End Property

	''' <summary>
	''' Theoretical, monoisotopic mass (including the modified residues)
	''' </summary>
	''' <value></value>
	''' <returns></returns>
	''' <remarks></remarks>
	Public ReadOnly Property MonoisotopicMass As Double
		Get
			Return mMonoisotopicMass
		End Get
	End Property

	''' <summary>
	''' Constructor using Sequence ID and mass
	''' </summary>
	''' <remarks></remarks>
    Public Sub New(SeqID As Integer, MonoisotopicMass As Double)
        Me.New(SeqID, MonoisotopicMass, 0, String.Empty)
    End Sub

    ''' <summary>
    ''' Constructor using Sequence ID, mass, mod count, and list of modifications
    ''' </summary>
    ''' <param name="SeqID">Sequence ID</param>
    ''' <param name="MonoisotopicMass">Theoretical, monoisotopic mass (including the modified residues)</param>
    ''' <param name="ModCount">Number of modifications</param>
    ''' <param name="ModDescription">Comma-separated list of modifications, for example "itrac:1,Phosph:3,IodoAcet:15"</param>
    ''' <remarks></remarks>
    Public Sub New(SeqID As Integer, MonoisotopicMass As Double, ModCount As Integer, ModDescription As String)
        mSeqID = SeqID
        mMonoisotopicMass = MonoisotopicMass
        mModCount = ModCount

        If String.IsNullOrEmpty(ModDescription) Then
            mModDescription = String.Empty
        Else
            mModDescription = ModDescription
        End If
    End Sub

    ''' <summary>
    ''' Update the monoisotopic mass for this sequence
    ''' </summary>
    ''' <param name="monoMass"></param>
    ''' <remarks></remarks>
    Public Sub UpdateMonoisotopicMass(monoMass As Double)
        mMonoisotopicMass = monoMass
    End Sub
End Class
