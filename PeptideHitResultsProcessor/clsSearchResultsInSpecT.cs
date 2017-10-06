Option Strict On

' This class is used to track the peptide details for an InSpecT search result
' See clsSearchResultsBaseClass for additional information
'
' -------------------------------------------------------------------------------
' Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
' Program started August 19, 2008
'
' E-mail: john.sandoval@pnnl.gov
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

Imports PHRPReader

Public Class clsSearchResultsInSpecT
    Inherits clsSearchResultsBaseClass

#Region "Classwide Variables"

    ' Scan: tracked by the base class
    ' Annotation: aka Peptide, which is tracked by the base class
    ' Protein: tracked by the base class

#End Region

#Region "Properties"

    ' Auto-Properties
    Public Property SpectrumFile As String
    Public Property MQScore As String
    Public Property Length As String
    Public Property TotalPRMScore As String
    Public Property MedianPRMScore As String
    Public Property FractionY As String
    Public Property FractionB As String
    Public Property Intensity As String
    Public Property NTT As String
    Public Property pValue As String
    Public Property FScore As String
    Public Property DeltaScore As String
    Public Property DeltaScoreOther As String
    Public Property DeltaNormMQScore As String
    Public Property DeltaNormTotalPRMScore As String
    Public Property RankTotalPRMScore As String
    Public Property RankFScore As String
    Public Property RecordNumber As String
    Public Property DBFilePos As String
    Public Property SpecFilePos As String
    Public Property PrecursorMZ As String
    Public Property PrecursorError As String

#End Region

    Public Sub New(ByRef objPeptideMods As clsPeptideModificationContainer, peptideSeqMassCalculator As clsPeptideMassCalculator)
        ' Note that the following call will call both the base class's Clear sub and this class's Clear Sub
        MyBase.New(objPeptideMods, peptideSeqMassCalculator)
    End Sub

    Public Overrides Sub Clear()
        MyBase.Clear()

        SpectrumFile = String.Empty
        MQScore = String.Empty
        Length = String.Empty
        TotalPRMScore = String.Empty
        MedianPRMScore = String.Empty
        FractionY = String.Empty
        FractionB = String.Empty
        Intensity = String.Empty
        NTT = String.Empty
        pValue = String.Empty
        FScore = String.Empty
        DeltaScore = String.Empty
        DeltaScoreOther = String.Empty
        DeltaNormMQScore = String.Empty
        DeltaNormTotalPRMScore = String.Empty
        RankTotalPRMScore = String.Empty
        RankFScore = String.Empty
        RecordNumber = String.Empty
        DBFilePos = String.Empty
        SpecFilePos = String.Empty
        PrecursorMZ = String.Empty
        PrecursorError = String.Empty
    End Sub
End Class
