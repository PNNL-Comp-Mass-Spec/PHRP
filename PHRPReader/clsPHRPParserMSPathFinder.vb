'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 07/17/2015
'
' This class parses data lines from mspath_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserMSPathFinder
    Inherits clsPHRPParser
    
#Region "Constants"
    Public Const DATA_COLUMN_ResultID As String = "ResultID"
    Public Const DATA_COLUMN_Scan As String = "Scan"
    Public Const DATA_COLUMN_Charge As String = "Charge"
    Public Const DATA_COLUMN_MostAbundantIsotopeMz As String = "MostAbundantIsotopeMz"
    Public Const DATA_COLUMN_Mass As String = "Mass"
    Public Const DATA_COLUMN_Sequence As String = "Sequence"
    Public Const DATA_COLUMN_Modifications As String = "Modifications"
    Public Const DATA_COLUMN_Composition As String = "Composition"
    Public Const DATA_COLUMN_Protein As String = "ProteinName"
    Public Const DATA_COLUMN_ProteinDesc As String = "ProteinDesc"
    Public Const DATA_COLUMN_ProteinLength As String = "ProteinLength"
    Public Const DATA_COLUMN_ResidueStart As String = "ResidueStart"
    Public Const DATA_COLUMN_ResidueEnd As String = "ResidueEnd"
    Public Const DATA_COLUMN_MatchedFragments As String = "MatchedFragments"
    Public Const DATA_COLUMN_QValue As String = "QValue"
    Public Const DATA_COLUMN_PepQValue As String = "PepQValue"

    Public Const FILENAME_SUFFIX_SYN As String = "_mspath_syn.txt"
    Public Const FILENAME_SUFFIX_FHT As String = "_mspath_fht.txt"

    Protected Const MSPathFinder_SEARCH_ENGINE_NAME As String = "MSPathFinder"

#End Region


#Region "Properties"

    Public Overrides ReadOnly Property PHRPFirstHitsFileName() As String
        Get
            Return GetPHRPFirstHitsFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPModSummaryFileName() As String
        Get
            Return GetPHRPModSummaryFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPPepToProteinMapFileName() As String
        Get
            Return GetPHRPPepToProteinMapFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPProteinModsFileName() As String
        Get
            Return GetPHRPProteinModsFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPSynopsisFileName() As String
        Get
            Return GetPHRPSynopsisFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPResultToSeqMapFileName() As String
        Get
            Return GetPHRPResultToSeqMapFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPSeqInfoFileName() As String
        Get
            Return GetPHRPSeqInfoFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property PHRPSeqToProteinMapFileName() As String
        Get
            Return GetPHRPSeqToProteinMapFileName(mDatasetName)
        End Get
    End Property

    Public Overrides ReadOnly Property SearchEngineName() As String
        Get
            Return GetSearchEngineName()
        End Get
    End Property

#End Region

    ''' <summary>
    ''' Constructor; assumes blnLoadModsAndSeqInfo=True
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <remarks></remarks>
    Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String)
        Me.New(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo:=True)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="blnLoadModsAndSeqInfo">If True, then load the ModSummary file and SeqInfo files</param>
    ''' <remarks></remarks>
    Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal blnLoadModsAndSeqInfo As Boolean)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSPathFinder, blnLoadModsAndSeqInfo)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
    ''' <remarks></remarks>
    Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal startupOptions As clsPHRPStartupOptions)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSPathFinder, startupOptions)
    End Sub

    Protected Overrides Sub DefineColumnHeaders()

        mColumnHeaders.Clear()

        ' Define the default column mapping	
        AddHeaderColumn(DATA_COLUMN_ResultID)

        AddHeaderColumn(DATA_COLUMN_Scan)
        AddHeaderColumn(DATA_COLUMN_Charge)
        AddHeaderColumn(DATA_COLUMN_MostAbundantIsotopeMz)
        AddHeaderColumn(DATA_COLUMN_Mass)
        AddHeaderColumn(DATA_COLUMN_Sequence)
        AddHeaderColumn(DATA_COLUMN_Modifications)
        AddHeaderColumn(DATA_COLUMN_Composition)
        AddHeaderColumn(DATA_COLUMN_Protein)
        AddHeaderColumn(DATA_COLUMN_ProteinDesc)
        AddHeaderColumn(DATA_COLUMN_ProteinLength)
        AddHeaderColumn(DATA_COLUMN_ResidueStart)
        AddHeaderColumn(DATA_COLUMN_ResidueEnd)
        AddHeaderColumn(DATA_COLUMN_MatchedFragments)
        AddHeaderColumn(DATA_COLUMN_QValue)
        AddHeaderColumn(DATA_COLUMN_PepQValue)

    End Sub

    Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
        ' MSPathFinder does not have a first-hits file; just the _syn.txt file
        Return String.Empty
    End Function

    Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_mspath_syn_ModSummary.txt"
    End Function

    Public Shared Function GetPHRPPepToProteinMapFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_mspath_PepToProtMapMTS.txt"
    End Function

    Public Shared Function GetPHRPProteinModsFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_mspath_syn_ProteinMods.txt"
    End Function

    Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_SYN
    End Function

    Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_mspath_syn_ResultToSeqMap.txt"
    End Function

    Public Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_mspath_syn_SeqInfo.txt"
    End Function

    Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_mspath_syn_SeqToProteinMap.txt"
    End Function

    Public Shared Function GetSearchEngineName() As String
        Return MSPathFinder_SEARCH_ENGINE_NAME
    End Function

    ''' <summary>
    ''' Parses the specified MSPathFinder parameter file
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName"></param>
    ''' <param name="objSearchEngineParams"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, <Out()> ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Dim blnSuccess As Boolean

        objSearchEngineParams = New clsSearchEngineParameters(MSPathFinder_SEARCH_ENGINE_NAME)

        blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

        ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

        Return blnSuccess

    End Function

    Protected Function ReadSearchEngineParamFile(ByVal strSearchEngineParamFileName As String, ByVal objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Try
            Dim blnSuccess = ReadKeyValuePairSearchEngineParamFile(MSPathFinder_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, ePeptideHitResultType.MSGFDB, objSearchEngineParams)

            If Not blnSuccess Then
                Return False
            End If

            objSearchEngineParams.Enzyme = "no_enzyme"
            objSearchEngineParams.MinNumberTermini = 0

            ' Determine the precursor mass tolerance (will store 0 if a problem or not found)
            Dim dblTolerancePPM As Double
            objSearchEngineParams.PrecursorMassToleranceDa = clsPHRPParserMSGFDB.DeterminePrecursorMassTolerance(objSearchEngineParams, dblTolerancePPM)
            objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM

            Return True

        Catch ex As Exception
            ReportError("Error in ReadSearchEngineParamFile: " & ex.Message)
            Return False
        End Try

    End Function

    ''' <summary>
    ''' Parse the data line read from a PHRP results file
    ''' </summary>
    ''' <param name="strLine">Data line</param>
    ''' <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
    ''' <param name="objPSM">clsPSM object (output)</param>
    ''' <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
    ''' <returns>True if success, false if an error</returns>
    ''' <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields</remarks>
    Public Overrides Function ParsePHRPDataLine(ByVal strLine As String, ByVal intLinesRead As Integer, <Out()> ByRef objPSM As clsPSM, ByVal fastReadMode As Boolean) As Boolean

        objPSM = New clsPSM()

        Try
            Dim strColumns() As String = strLine.Split(ControlChars.Tab)
            Dim blnSuccess = False

            With objPSM
                .DataLineText = strLine
                .ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
                If .ScanNumber = -100 Then
                    ' Data line is not valid
                Else
                    .ResultID = LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0)
                    .ScoreRank = 1

                    Dim strSequence = LookupColumnValue(strColumns, DATA_COLUMN_Sequence, mColumnHeaders)

                    If fastReadMode Then
                        .SetPeptide(strSequence, blnUpdateCleanSequence:=False)
                    Else
                        .SetPeptide(strSequence, mCleavageStateCalculator)
                    End If

                    .Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

                    Dim strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
                    .AddProtein(strProtein)

                    ' Store the sequence mass as the "precursor" mass, though MSPathFinderT results are from MS1 spectra, and thus we didn't do MS/MS on a precursor
                    .PrecursorNeutralMass = LookupColumnValue(strColumns, DATA_COLUMN_Mass, mColumnHeaders, 0.0#)

                    ' Thus collision mode, precursor neutral mass, etc. are not applicable
                    ' .CollisionMode =
                    ' .PrecursorNeutralMass = 
                    ' .MassErrorDa =
                    ' .MassErrorPPM = 

                    ' .MSGFSpecProb = 

                    blnSuccess = True
                End If
            End With

            If Not blnSuccess Then
                Return False
            End If

            If Not fastReadMode Then
                UpdatePSMUsingSeqInfo(objPSM)
            End If

            ' Store the remaining data

            AddScore(objPSM, strColumns, DATA_COLUMN_MostAbundantIsotopeMz)
            AddScore(objPSM, strColumns, DATA_COLUMN_Modifications)
            AddScore(objPSM, strColumns, DATA_COLUMN_Composition)
            AddScore(objPSM, strColumns, DATA_COLUMN_ProteinDesc)
            AddScore(objPSM, strColumns, DATA_COLUMN_ProteinLength)
            AddScore(objPSM, strColumns, DATA_COLUMN_ResidueStart)
            AddScore(objPSM, strColumns, DATA_COLUMN_ResidueEnd)
            AddScore(objPSM, strColumns, DATA_COLUMN_MatchedFragments)
            AddScore(objPSM, strColumns, DATA_COLUMN_QValue)
            AddScore(objPSM, strColumns, DATA_COLUMN_PepQValue)

            Return True

        Catch ex As Exception
            MyBase.ReportError("Error parsing line " & intLinesRead & " in the MSGFDB data file: " & ex.Message)
            Return False
        End Try

    End Function

End Class
