Option Strict On

' This class is used to track peptide sequences and their modification descriptors
' It assigns a unique integer ID to each combination of sequence and modification description
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 11, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
' -------------------------------------------------------------------------------
' 
' Licensed under the Apache License, Version 2.0; you may not use this file except
' in compliance with the License.  You may obtain a copy of the License at 
' http://www.apache.org/licenses/LICENSE-2.0
'
' Notice: This computer software was prepared by Battelle Memorial Institute, 
' hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
' Department of Energy (DOE).  All rights in the computer software are reserved 
' by DOE on behalf of the United States Government and the Contractor as 
' provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
' WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
' SOFTWARE.  This notice including this sentence must appear on any copies of 
' this computer software.

Public Class clsUniqueSequencesContainer


#Region "Constants and Enums"
    Protected Const DEFAULT_INITIAL_SEQ_ID As Integer = 1
    Protected Const SEQUENCE_MOD_DESC_SEP As Char = "_"c
#End Region

#Region "Classwide Variables"
    Protected htMasterSequences As Hashtable
    Protected mNextUniqueSeqID As Integer
#End Region

#Region "Properties"
    Public ReadOnly Property UniqueSequenceCount() As Integer
        Get
            Return htMasterSequences.Count
        End Get
    End Property
#End Region

    Public Sub New()
        InitializeLocalVariables()
    End Sub

    Public Sub Clear()
        Me.Clear(DEFAULT_INITIAL_SEQ_ID)
    End Sub

    Public Sub Clear(intInitialSeqID As Integer)
        ' Clears htMasterSequences and resets mNextUniqueSeqID to intInitialSeqID
        mNextUniqueSeqID = intInitialSeqID
        If htMasterSequences Is Nothing Then
            htMasterSequences = New Hashtable
        Else
            htMasterSequences.Clear()
        End If
    End Sub

    Public Function GetNextUniqueSequenceID(strSequence As String, strModDescription As String, ByRef blnExistingSequenceFound As Boolean) As Integer

        Dim intUniqueSeqID As Integer
        Dim strKey As String

        blnExistingSequenceFound = False

        Try
            If strSequence Is Nothing Then strSequence = String.Empty
            If strModDescription Is Nothing Then strModDescription = String.Empty

            strKey = strSequence & SEQUENCE_MOD_DESC_SEP & strModDescription

            If htMasterSequences.ContainsKey(strKey) Then
                intUniqueSeqID = CInt(htMasterSequences(strKey))
                blnExistingSequenceFound = True
            Else
                htMasterSequences.Add(strKey, mNextUniqueSeqID)
                intUniqueSeqID = mNextUniqueSeqID
                mNextUniqueSeqID += 1
            End If
        Catch ex As Exception
            intUniqueSeqID = 0
        End Try

        Return intUniqueSeqID
    End Function

    Protected Sub InitializeLocalVariables()
        Me.Clear()
    End Sub
End Class
