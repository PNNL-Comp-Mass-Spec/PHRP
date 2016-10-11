'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 11/28/2012
'
' This class parses data lines from MSAlign msalign_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserMSAlign
    Inherits clsPHRPParser

#Region "Constants"
    Public Const DATA_COLUMN_ResultID As String = "ResultID"
    Public Const DATA_COLUMN_Scan As String = "Scan"
    Public Const DATA_COLUMN_Prsm_ID As String = "Prsm_ID"
    Public Const DATA_COLUMN_Spectrum_ID As String = "Spectrum_ID"
    Public Const DATA_COLUMN_Charge As String = "Charge"
    Public Const DATA_COLUMN_PrecursorMZ As String = "PrecursorMZ"
    Public Const DATA_COLUMN_DelM As String = "DelM"
    Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"
    Public Const DATA_COLUMN_MH As String = "MH"
    Public Const DATA_COLUMN_Peptide As String = "Peptide"
    Public Const DATA_COLUMN_Protein As String = "Protein"
    Public Const DATA_COLUMN_Protein_Mass As String = "Protein_Mass"
    Public Const DATA_COLUMN_Unexpected_Mod_Count As String = "Unexpected_Mod_Count"
    Public Const DATA_COLUMN_Peak_Count As String = "Peak_Count"
    Public Const DATA_COLUMN_Matched_Peak_Count As String = "Matched_Peak_Count"
    Public Const DATA_COLUMN_Matched_Fragment_Ion_Count As String = "Matched_Fragment_Ion_Count"
    Public Const DATA_COLUMN_PValue As String = "PValue"
    Public Const DATA_COLUMN_Rank_PValue As String = "Rank_PValue"
    Public Const DATA_COLUMN_EValue As String = "EValue"
    Public Const DATA_COLUMN_FDR As String = "FDR"
    Public Const DATA_COLUMN_Species_ID As String = "Species_ID"
    Public Const DATA_COLUMN_FragMethod As String = "FragMethod"

    Public Const FILENAME_SUFFIX_SYN As String = "_msalign_syn.txt"
    Public Const FILENAME_SUFFIX_FHT As String = "_msalign_fht.txt"

    Private Const MSAlign_SEARCH_ENGINE_NAME As String = "MSAlign"
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
    Public Sub New(strDatasetName As String, strInputFilePath As String)
        Me.New(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo:=True)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="blnLoadModsAndSeqInfo">If True, then load the ModSummary file and SeqInfo files</param>
    ''' <remarks></remarks>
    Public Sub New(strDatasetName As String, strInputFilePath As String, blnLoadModsAndSeqInfo As Boolean)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSAlign, blnLoadModsAndSeqInfo)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
    ''' <remarks></remarks>
    Public Sub New(strDatasetName As String, strInputFilePath As String, startupOptions As clsPHRPStartupOptions)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSAlign, startupOptions)
    End Sub

    Protected Overrides Sub DefineColumnHeaders()

        mColumnHeaders.Clear()

        ' Define the default column mapping	
        AddHeaderColumn(DATA_COLUMN_ResultID)
        AddHeaderColumn(DATA_COLUMN_Scan)
        AddHeaderColumn(DATA_COLUMN_Prsm_ID)
        AddHeaderColumn(DATA_COLUMN_Spectrum_ID)
        AddHeaderColumn(DATA_COLUMN_Charge)
        AddHeaderColumn(DATA_COLUMN_PrecursorMZ)
        AddHeaderColumn(DATA_COLUMN_DelM)
        AddHeaderColumn(DATA_COLUMN_DelM_PPM)
        AddHeaderColumn(DATA_COLUMN_MH)
        AddHeaderColumn(DATA_COLUMN_Peptide)
        AddHeaderColumn(DATA_COLUMN_Protein)
        AddHeaderColumn(DATA_COLUMN_Protein_Mass)
        AddHeaderColumn(DATA_COLUMN_Unexpected_Mod_Count)
        AddHeaderColumn(DATA_COLUMN_Peak_Count)
        AddHeaderColumn(DATA_COLUMN_Matched_Peak_Count)
        AddHeaderColumn(DATA_COLUMN_Matched_Fragment_Ion_Count)
        AddHeaderColumn(DATA_COLUMN_PValue)
        AddHeaderColumn(DATA_COLUMN_Rank_PValue)
        AddHeaderColumn(DATA_COLUMN_EValue)
        AddHeaderColumn(DATA_COLUMN_FDR)
        AddHeaderColumn(DATA_COLUMN_Species_ID)
        AddHeaderColumn(DATA_COLUMN_FragMethod)

    End Sub

    ''' <summary>
    ''' Determines the precursor mass tolerance
    ''' </summary>
    ''' <param name="objSearchEngineParams"></param>
    ''' <param name="dblTolerancePPM">Precursor mass tolerance, in ppm</param>
    ''' <returns>Precursor tolerance, in Da</returns>
    ''' <remarks></remarks>
    Private Function DeterminePrecursorMassTolerance(objSearchEngineParams As clsSearchEngineParameters, <Out()> ByRef dblTolerancePPM As Double) As Double
        Dim strTolerance As String = String.Empty

        Dim dblToleranceDa As Double = 0

        dblTolerancePPM = 0

        If objSearchEngineParams.Parameters.TryGetValue("errorTolerance", strTolerance) Then
            ' Parent mass tolerance, in ppm
            If Double.TryParse(strTolerance, dblTolerancePPM) Then
                dblToleranceDa = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000)
            End If
        End If

        Return dblToleranceDa

    End Function

    Public Shared Function GetPHRPFirstHitsFileName(strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_FHT
    End Function

    Public Shared Function GetPHRPModSummaryFileName(strDatasetName As String) As String
        Return strDatasetName & "_msalign_syn_ModSummary.txt"
    End Function

    Public Shared Function GetPHRPPepToProteinMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_msalign_PepToProtMapMTS.txt"
    End Function

    Public Shared Function GetPHRPProteinModsFileName(strDatasetName As String) As String
        Return strDatasetName & "_msalign_syn_ProteinMods.txt"
    End Function

    Public Shared Function GetPHRPSynopsisFileName(strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_SYN
    End Function

    Public Shared Function GetPHRPResultToSeqMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_msalign_syn_ResultToSeqMap.txt"
    End Function

    Public Shared Function GetPHRPSeqInfoFileName(strDatasetName As String) As String
        Return strDatasetName & "_msalign_syn_SeqInfo.txt"
    End Function

    Public Shared Function GetPHRPSeqToProteinMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_msalign_syn_SeqToProteinMap.txt"
    End Function

    Public Shared Function GetSearchEngineName() As String
        Return MSAlign_SEARCH_ENGINE_NAME
    End Function

    ''' <summary>
    ''' Parses the specified MSAlign parameter file
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName"></param>
    ''' <param name="objSearchEngineParams"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Overrides Function LoadSearchEngineParameters(strSearchEngineParamFileName As String, <Out()> ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Dim blnSuccess As Boolean

        objSearchEngineParams = New clsSearchEngineParameters(MSAlign_SEARCH_ENGINE_NAME)

        blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

        ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

        Return blnSuccess

    End Function

    Private Function ReadSearchEngineParamFile(strSearchEngineParamFileName As String, objSearchEngineParams As clsSearchEngineParameters) As Boolean
        Dim strSettingValue As String = String.Empty
        Dim objModDef As clsModificationDefinition
        Dim blnSuccess As Boolean

        Try
            blnSuccess = ReadKeyValuePairSearchEngineParamFile(MSAlign_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, ePeptideHitResultType.MSAlign, objSearchEngineParams)

            If blnSuccess Then

                ' For MSGF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
                ' But MSAlign does not have a _ModDefs.txt file because it performs a blind search
                ' The user can define static mods on cysteine; check for these now

                If objSearchEngineParams.Parameters.TryGetValue("cysteineProtection", strSettingValue) Then

                    Select Case strSettingValue
                        Case "C57"
                            objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 57.0215, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcet")
                            objSearchEngineParams.AddModification(objModDef)
                        Case "C58"
                            objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, 58.0055, "C", clsModificationDefinition.eModificationTypeConstants.StaticMod, "IodoAcid")
                            objSearchEngineParams.AddModification(objModDef)
                    End Select

                End If

                ' Determine the precursor mass tolerance (will store 0 if a problem or not found)
                Dim dblTolerancePPM As Double
                objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, dblTolerancePPM)
                objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM
            End If

        Catch ex As Exception
            ReportError("Error in ReadSearchEngineParamFile: " & ex.Message)
        End Try

        Return blnSuccess

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
    Public Overrides Function ParsePHRPDataLine(strLine As String, intLinesRead As Integer, <Out()> ByRef objPSM As clsPSM, fastReadMode As Boolean) As Boolean

        Dim strColumns() As String = strLine.Split(ControlChars.Tab)
        Dim strPeptide As String
        Dim strProtein As String

        Dim dblPrecursorMZ As Double

        Dim blnSuccess As Boolean

        objPSM = New clsPSM()

        Try

            With objPSM
                .DataLineText = strLine
                .ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
                If .ScanNumber = -100 Then
                    ' Data line is not valid
                Else

                    .ResultID = LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0)
                    .ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_Rank_PValue, mColumnHeaders, 1)

                    strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)

                    If fastReadMode Then
                        .SetPeptide(strPeptide, blnUpdateCleanSequence:=False)
                    Else
                        .SetPeptide(strPeptide, mCleavageStateCalculator)
                    End If

                    .Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

                    strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
                    .AddProtein(strProtein)

                    dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
                    .PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

                    .MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
                    .MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

                    blnSuccess = True
                End If
            End With

            If blnSuccess Then
                If Not fastReadMode Then
                    UpdatePSMUsingSeqInfo(objPSM)
                End If

                ' Store the remaining scores
                AddScore(objPSM, strColumns, DATA_COLUMN_Prsm_ID)
                AddScore(objPSM, strColumns, DATA_COLUMN_Spectrum_ID)

                AddScore(objPSM, strColumns, DATA_COLUMN_MH)

                AddScore(objPSM, strColumns, DATA_COLUMN_Protein_Mass)
                AddScore(objPSM, strColumns, DATA_COLUMN_Unexpected_Mod_Count)
                AddScore(objPSM, strColumns, DATA_COLUMN_Peak_Count)
                AddScore(objPSM, strColumns, DATA_COLUMN_Matched_Peak_Count)
                AddScore(objPSM, strColumns, DATA_COLUMN_Matched_Fragment_Ion_Count)

                AddScore(objPSM, strColumns, DATA_COLUMN_PValue)
                AddScore(objPSM, strColumns, DATA_COLUMN_EValue)
                AddScore(objPSM, strColumns, DATA_COLUMN_FDR)

                AddScore(objPSM, strColumns, DATA_COLUMN_Species_ID)
                AddScore(objPSM, strColumns, DATA_COLUMN_FragMethod)

            End If

        Catch ex As Exception
            MyBase.ReportError("Error parsing line " & intLinesRead & " in the MSAlign data file: " & ex.Message)
        End Try

        Return blnSuccess

    End Function

End Class
