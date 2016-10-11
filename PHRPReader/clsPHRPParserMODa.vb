'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/02/2014
'
' This class parses data lines from moda_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserMODa
    Inherits clsPHRPParser

#Region "Constants"
    Public Const DATA_COLUMN_ResultID As String = "ResultID"
    Public Const DATA_COLUMN_Scan As String = "Scan"
    Public Const DATA_COLUMN_Spectrum_Index As String = "Spectrum_Index"
    Public Const DATA_COLUMN_Charge As String = "Charge"
    Public Const DATA_COLUMN_PrecursorMZ As String = "PrecursorMZ"
    Public Const DATA_COLUMN_DelM As String = "DelM"
    Public Const DATA_COLUMN_DelM_PPM As String = "DelM_PPM"
    Public Const DATA_COLUMN_MH As String = "MH"
    Public Const DATA_COLUMN_Peptide As String = "Peptide"
    Public Const DATA_COLUMN_Protein As String = "Protein"
    Public Const DATA_COLUMN_Score As String = "Score"
    Public Const DATA_COLUMN_Probability As String = "Probability"
    Public Const DATA_COLUMN_Rank_Probability As String = "Rank_Probability"
    Public Const DATA_COLUMN_Peptide_Position As String = "Peptide_Position"
    Public Const DATA_COLUMN_QValue As String = "QValue"

    Public Const FILENAME_SUFFIX_SYN As String = "_moda_syn.txt"
    Public Const FILENAME_SUFFIX_FHT As String = "_moda_fht.txt"

    Private Const MODa_SEARCH_ENGINE_NAME As String = "MODa"
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
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MODa, blnLoadModsAndSeqInfo)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
    ''' <remarks></remarks>
    Public Sub New(strDatasetName As String, strInputFilePath As String, startupOptions As clsPHRPStartupOptions)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MODa, startupOptions)
    End Sub

    Protected Overrides Sub DefineColumnHeaders()

        mColumnHeaders.Clear()

        ' Define the default column mapping	
        AddHeaderColumn(DATA_COLUMN_ResultID)
        AddHeaderColumn(DATA_COLUMN_Scan)
        AddHeaderColumn(DATA_COLUMN_Spectrum_Index)
        AddHeaderColumn(DATA_COLUMN_Charge)
        AddHeaderColumn(DATA_COLUMN_PrecursorMZ)
        AddHeaderColumn(DATA_COLUMN_DelM)
        AddHeaderColumn(DATA_COLUMN_DelM_PPM)
        AddHeaderColumn(DATA_COLUMN_MH)
        AddHeaderColumn(DATA_COLUMN_Peptide)
        AddHeaderColumn(DATA_COLUMN_Protein)
        AddHeaderColumn(DATA_COLUMN_Score)
        AddHeaderColumn(DATA_COLUMN_Probability)
        AddHeaderColumn(DATA_COLUMN_Rank_Probability)
        AddHeaderColumn(DATA_COLUMN_Peptide_Position)
        AddHeaderColumn(DATA_COLUMN_QValue)

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

        If objSearchEngineParams.Parameters.TryGetValue("PPMTolerance", strTolerance) Then
            ' Parent mass tolerance, in ppm
            If Double.TryParse(strTolerance, dblTolerancePPM) Then
                dblToleranceDa = clsPeptideMassCalculator.PPMToMass(dblTolerancePPM, 2000)
            End If

        ElseIf objSearchEngineParams.Parameters.TryGetValue("PeptTolerance", strTolerance) Then
            ' Parent mass tolerance, in Da
            Double.TryParse(strTolerance, dblToleranceDa)

            ' Convert from dalton to PPM (assuming a mass of 2000 m/z)
            dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000)
        End If

        Return dblToleranceDa

    End Function

    Public Shared Function GetPHRPFirstHitsFileName(strDatasetName As String) As String
        ' MODa does not have a first-hits file; just the _syn.txt file
        Return String.Empty
    End Function

    Public Shared Function GetPHRPModSummaryFileName(strDatasetName As String) As String
        Return strDatasetName & "_moda_syn_ModSummary.txt"
    End Function

    Public Shared Function GetPHRPPepToProteinMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_moda_PepToProtMapMTS.txt"
    End Function

    Public Shared Function GetPHRPProteinModsFileName(strDatasetName As String) As String
        Return strDatasetName & "_moda_syn_ProteinMods.txt"
    End Function

    Public Shared Function GetPHRPSynopsisFileName(strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_SYN
    End Function

    Public Shared Function GetPHRPResultToSeqMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_moda_syn_ResultToSeqMap.txt"
    End Function

    Public Shared Function GetPHRPSeqInfoFileName(strDatasetName As String) As String
        Return strDatasetName & "_moda_syn_SeqInfo.txt"
    End Function

    Public Shared Function GetPHRPSeqToProteinMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_moda_syn_SeqToProteinMap.txt"
    End Function

    Public Shared Function GetSearchEngineName() As String
        Return MODa_SEARCH_ENGINE_NAME
    End Function

    ''' <summary>
    ''' Parses the specified MODa parameter file
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName"></param>
    ''' <param name="objSearchEngineParams"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Overrides Function LoadSearchEngineParameters(strSearchEngineParamFileName As String, <Out()> ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Dim blnSuccess As Boolean

        objSearchEngineParams = New clsSearchEngineParameters(MODa_SEARCH_ENGINE_NAME)

        blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

        ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

        Return blnSuccess

    End Function

    Private Function ReadSearchEngineParamFile(strSearchEngineParamFileName As String, objSearchEngineParams As clsSearchEngineParameters) As Boolean
        Dim strSettingValue As String = String.Empty
        Dim objModDef As clsModificationDefinition

        Try
            Dim blnSuccess = ReadKeyValuePairSearchEngineParamFile(MODa_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, ePeptideHitResultType.MODa, objSearchEngineParams)

            If Not blnSuccess Then
                Return False
            End If

            ' For MSGF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
            ' But MODa does not have a _ModDefs.txt file because it performs a blind search
            ' The user can define static mods on any of the residues, plus the peptide terminii; check for these now

            Dim lstResiduesToFind = New List(Of String) From {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}

            ' This dictionary tracks the static mod names we will look for
            ' It is populated using the amino acid letters in lstResiduesToFind, plus also the N and T terminus tags
            Dim dctResiduesAndSymbols = New Dictionary(Of String, String)

            For Each residueSymbol In lstResiduesToFind
                dctResiduesAndSymbols.Add(residueSymbol, residueSymbol)
            Next

            dctResiduesAndSymbols.Add("ADD_NTerm", N_TERMINAL_PEPTIDE_SYMBOL_DMS)
            dctResiduesAndSymbols.Add("ADD_CTerm", C_TERMINAL_PEPTIDE_SYMBOL_DMS)

            For Each residueSpec In dctResiduesAndSymbols
                Dim strKey = "ADD_" & residueSpec.Key

                If objSearchEngineParams.Parameters.TryGetValue(strKey, strSettingValue) Then
                    Dim modMassDa As Double

                    If Double.TryParse(strSettingValue, modMassDa) Then
                        If Math.Abs(modMassDa) > Single.Epsilon Then

                            Dim eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod
                            If residueSpec.Value = N_TERMINAL_PEPTIDE_SYMBOL_DMS OrElse residueSpec.Value = C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                                eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                            End If

                            objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, residueSpec.Value, eModType, "Mod" & modMassDa.ToString("0"))
                            objSearchEngineParams.AddModification(objModDef)
                        End If
                    End If

                End If
            Next

            ' Determine the precursor mass tolerance (will store 0 if a problem or not found)
            Dim dblTolerancePPM As Double
            objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, dblTolerancePPM)
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
                    .ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_Rank_Probability, mColumnHeaders, 1)

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
                AddScore(objPSM, strColumns, DATA_COLUMN_Spectrum_Index)

                AddScore(objPSM, strColumns, DATA_COLUMN_MH)

                AddScore(objPSM, strColumns, DATA_COLUMN_Score)
                AddScore(objPSM, strColumns, DATA_COLUMN_Probability)
                AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Position)
                AddScore(objPSM, strColumns, DATA_COLUMN_QValue)

            End If

        Catch ex As Exception
            MyBase.ReportError("Error parsing line " & intLinesRead & " in the MODa data file: " & ex.Message)
        End Try

        Return blnSuccess

    End Function

End Class
