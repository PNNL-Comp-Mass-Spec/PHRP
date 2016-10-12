'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class parses data lines from MSGFDB msgfdb_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports PHRPReader.clsPHRPReader
Imports PHRPReader.clsMSGFPlusParamFileModExtractor
Imports System.Text.RegularExpressions

Public Class clsPHRPParserMSGFDB
    Inherits clsPHRPParser

#Region "Constants"
    Public Const DATA_COLUMN_ResultID = "ResultID"
    Public Const DATA_COLUMN_Scan = "Scan"
    Public Const DATA_COLUMN_FragMethod = "FragMethod"
    Public Const DATA_COLUMN_SpecIndex = "SpecIndex"
    Public Const DATA_COLUMN_Charge = "Charge"
    Public Const DATA_COLUMN_PrecursorMZ = "PrecursorMZ"
    Public Const DATA_COLUMN_DelM = "DelM"
    Public Const DATA_COLUMN_DelM_PPM = "DelM_PPM"
    Public Const DATA_COLUMN_MH = "MH"
    Public Const DATA_COLUMN_Peptide = "Peptide"
    Public Const DATA_COLUMN_Protein = "Protein"
    Public Const DATA_COLUMN_NTT = "NTT"
    Public Const DATA_COLUMN_DeNovoScore = "DeNovoScore"
    Public Const DATA_COLUMN_MSGFScore = "MSGFScore"

    Public Const DATA_COLUMN_MSGFDB_SpecProb = "MSGFDB_SpecProb"                    ' MSGFDB
    Public Const DATA_COLUMN_Rank_MSGFDB_SpecProb = "Rank_MSGFDB_SpecProb"      ' MSGFDB

    Public Const DATA_COLUMN_MSGFDB_SpecEValue = "MSGFDB_SpecEValue"                ' MSGF+
    Public Const DATA_COLUMN_Rank_MSGFDB_SpecEValue = "Rank_MSGFDB_SpecEValue"  ' MSGF+

    Public Const DATA_COLUMN_PValue = "PValue"      ' MSGFDB
    Public Const DATA_COLUMN_EValue = "EValue"      ' MSGF+

    Public Const DATA_COLUMN_FDR = "FDR"                            ' MSGFDB; Only present if a Target/Decoy (TDA) search was used
    Public Const DATA_COLUMN_PepFDR = "PepFDR"                  ' MSGFDB; Only valid if a Target/Decoy (TDA) search was used; if EFDR is present, will contain 1 for every row

    Public Const DATA_COLUMN_QValue = "QValue"                  ' MSGF+ reports QValue instead of FDR
    Public Const DATA_COLUMN_PepQValue = "PepQValue"                ' MSGF+ reports pepQValue instead of PepFDR

    Public Const DATA_COLUMN_EFDR = "EFDR"                      ' Only present if a Target/Decoy (TDA) search was not used

    Public Const DATA_COLUMN_IMS_Scan = "IMS_Scan"
    Public Const DATA_COLUMN_IMS_Drift_Time = "IMS_Drift_Time"

    Public Const DATA_COLUMN_Isotope_Error = "IsotopeError"     ' Only reported by MSGF+

    Public Const FILENAME_SUFFIX_SYN = "_msgfdb_syn.txt"
    Public Const FILENAME_SUFFIX_FHT = "_msgfdb_fht.txt"

    Private Const MSGFDB_SEARCH_ENGINE_NAME = "MS-GFDB"

    Public Const CHARGE_CARRIER_MASS_PARAM_NAME As String = "ChargeCarrierMass"


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
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSGFDB, blnLoadModsAndSeqInfo)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
    ''' <remarks></remarks>
    Public Sub New(strDatasetName As String, strInputFilePath As String, startupOptions As clsPHRPStartupOptions)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MSGFDB, startupOptions)
    End Sub

    Protected Overrides Sub DefineColumnHeaders()

        mColumnHeaders.Clear()

        ' Define the default column mapping	
        AddHeaderColumn(DATA_COLUMN_ResultID)
        AddHeaderColumn(DATA_COLUMN_Scan)
        AddHeaderColumn(DATA_COLUMN_FragMethod)
        AddHeaderColumn(DATA_COLUMN_SpecIndex)
        AddHeaderColumn(DATA_COLUMN_Charge)
        AddHeaderColumn(DATA_COLUMN_PrecursorMZ)
        AddHeaderColumn(DATA_COLUMN_DelM)
        AddHeaderColumn(DATA_COLUMN_DelM_PPM)
        AddHeaderColumn(DATA_COLUMN_MH)
        AddHeaderColumn(DATA_COLUMN_Peptide)
        AddHeaderColumn(DATA_COLUMN_Protein)
        AddHeaderColumn(DATA_COLUMN_NTT)
        AddHeaderColumn(DATA_COLUMN_DeNovoScore)
        AddHeaderColumn(DATA_COLUMN_MSGFScore)
        AddHeaderColumn(DATA_COLUMN_MSGFDB_SpecProb)
        AddHeaderColumn(DATA_COLUMN_Rank_MSGFDB_SpecProb)
        AddHeaderColumn(DATA_COLUMN_PValue)
        AddHeaderColumn(DATA_COLUMN_FDR)
        AddHeaderColumn(DATA_COLUMN_EFDR)
        AddHeaderColumn(DATA_COLUMN_PepFDR)

        ' Add the MSGF+ columns
        AddHeaderColumn(DATA_COLUMN_MSGFDB_SpecEValue)
        AddHeaderColumn(DATA_COLUMN_Rank_MSGFDB_SpecEValue)
        AddHeaderColumn(DATA_COLUMN_EValue)

        AddHeaderColumn(DATA_COLUMN_QValue)
        AddHeaderColumn(DATA_COLUMN_PepQValue)

        AddHeaderColumn(DATA_COLUMN_IMS_Scan)
        AddHeaderColumn(DATA_COLUMN_IMS_Drift_Time)
        AddHeaderColumn(DATA_COLUMN_Isotope_Error)

    End Sub

    ''' <summary>
    ''' Determines the precursor mass tolerance for either MSGF+ or MSPathFinder
    ''' </summary>
    ''' <param name="objSearchEngineParams"></param>
    ''' <param name="dblTolerancePPM">Precursor mass tolerance, in ppm</param>
    ''' <returns>Precursor tolerance, in Da</returns>
    ''' <remarks></remarks>
    Public Shared Function DeterminePrecursorMassTolerance(
        objSearchEngineParams As clsSearchEngineParameters,
        <Out()> ByRef dblTolerancePPM As Double,
        resultType As ePeptideHitResultType) As Double

        Dim strTolerance = String.Empty
        Dim strToleranceSplit As String()

        Dim reMatch As Match
        Dim reExtraToleranceWithUnits = New Regex("([0-9.]+)([A-Za-z]+)", RegexOptions.Compiled Or RegexOptions.IgnoreCase)
        Dim reExtraToleranceNoUnits = New Regex("([0-9.]+)", RegexOptions.Compiled Or RegexOptions.IgnoreCase)

        Dim dblToleranceDa As Double

        dblTolerancePPM = 0

        If Not objSearchEngineParams.Parameters.TryGetValue("PMTolerance", strTolerance) Then
            Return dblToleranceDa
        End If

        ' Parent mass tolerance
        ' Might contain two values, separated by a comma
        strToleranceSplit = strTolerance.Split(","c)

        If strToleranceSplit Is Nothing Then
            Return dblToleranceDa
        End If

        For Each strItem In strToleranceSplit
            If strItem.Trim.StartsWith("#") Then Continue For

            If resultType = ePeptideHitResultType.MSPathFinder Then
                reMatch = reExtraToleranceNoUnits.Match(strItem)
            Else
                reMatch = reExtraToleranceWithUnits.Match(strItem)
            End If

            If Not reMatch.Success Then Continue For

            Dim dblToleranceCurrent As Double

            If Not Double.TryParse(reMatch.Groups(1).Value, dblToleranceCurrent) Then Continue For

            If resultType = ePeptideHitResultType.MSPathFinder Then
                ' Units are always ppm
                dblTolerancePPM = dblToleranceCurrent
                dblToleranceCurrent = clsPeptideMassCalculator.PPMToMass(dblToleranceCurrent, 2000)

            ElseIf reMatch.Groups.Count > 1 AndAlso reMatch.Groups(2).Value.ToLower().Contains("ppm") Then
                ' Ppm
                ' Convert from PPM to dalton (assuming a mass of 2000 m/z)
                dblTolerancePPM = dblToleranceCurrent
                dblToleranceCurrent = clsPeptideMassCalculator.PPMToMass(dblToleranceCurrent, 2000)

            End If

            dblToleranceDa = Math.Max(dblToleranceDa, dblToleranceCurrent)

        Next

        If Math.Abs(dblTolerancePPM) < Single.Epsilon And Math.Abs(dblToleranceDa) > Single.Epsilon Then
            dblTolerancePPM = clsPeptideMassCalculator.MassToPPM(dblToleranceDa, 2000)
        End If

        Return dblToleranceDa

    End Function

    ''' <summary>
    ''' Look for MSGF+ parameter ChargeCarrierMass
    ''' If defined, update chargeCarrierMass with the associated mass value and return True
    ''' Otherwise return false
    ''' </summary>
    ''' <param name="objSearchEngineParams"></param>
    ''' <param name="chargeCarrierMass"></param>
    ''' <returns></returns>
    ''' <remarks>This function is used by clsPHRPMassErrorValidator in the Analysis Manager</remarks>
    Public Shared Function GetCustomChargeCarrierMass(objSearchEngineParams As clsSearchEngineParameters, <Out()> ByRef chargeCarrierMass As Double) As Boolean

        Dim strValue As String = Nothing
        If objSearchEngineParams.Parameters.TryGetValue(CHARGE_CARRIER_MASS_PARAM_NAME, strValue) Then
            If Double.TryParse(strValue, chargeCarrierMass) Then
                Return True
            End If
        End If

        chargeCarrierMass = clsPeptideMassCalculator.MASS_PROTON
        Return False

    End Function

    Public Shared Function GetPHRPFirstHitsFileName(strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_FHT
    End Function

    Public Shared Function GetPHRPModSummaryFileName(strDatasetName As String) As String
        Return strDatasetName & "_msgfdb_syn_ModSummary.txt"
    End Function

    Public Shared Function GetPHRPPepToProteinMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_msgfdb_PepToProtMapMTS.txt"
    End Function

    Public Shared Function GetPHRPProteinModsFileName(strDatasetName As String) As String
        Return strDatasetName & "_msgfdb_syn_ProteinMods.txt"
    End Function

    Public Shared Function GetPHRPSynopsisFileName(strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_SYN
    End Function

    Public Shared Function GetPHRPResultToSeqMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_msgfdb_syn_ResultToSeqMap.txt"
    End Function

    Public Shared Function GetPHRPSeqInfoFileName(strDatasetName As String) As String
        Return strDatasetName & "_msgfdb_syn_SeqInfo.txt"
    End Function

    Public Shared Function GetPHRPSeqToProteinMapFileName(strDatasetName As String) As String
        Return strDatasetName & "_msgfdb_syn_SeqToProteinMap.txt"
    End Function

    Public Shared Function GetSearchEngineName() As String
        Return MSGFDB_SEARCH_ENGINE_NAME
    End Function

    ''' <summary>
    ''' Parses the specified MSGFDB (aka MS-GF+) parameter file
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName"></param>
    ''' <param name="objSearchEngineParams"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Overrides Function LoadSearchEngineParameters(strSearchEngineParamFileName As String, <Out()> ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Dim blnSuccess As Boolean

        objSearchEngineParams = New clsSearchEngineParameters(MSGFDB_SEARCH_ENGINE_NAME, mModInfo)

        blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

        ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

        Return blnSuccess

    End Function

    Private Function ReadSearchEngineParamFile(strSearchEngineParamFileName As String, objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Try
            mPeptideMassCalculator.ResetAminoAcidMasses()

            Dim success = ReadKeyValuePairSearchEngineParamFile(MSGFDB_SEARCH_ENGINE_NAME, strSearchEngineParamFileName, ePeptideHitResultType.MSGFDB, objSearchEngineParams)

            If Not success Then
                Return False
            End If

            Dim strSettingValue = String.Empty
            Dim intValue As Integer

            ' Determine the enzyme name
            If objSearchEngineParams.Parameters.TryGetValue("enzymeid", strSettingValue) Then
                If Integer.TryParse(strSettingValue, intValue) Then
                    Select Case intValue
                        Case 0 : objSearchEngineParams.Enzyme = "no_enzyme"
                        Case 1 : objSearchEngineParams.Enzyme = "trypsin"
                        Case 2 : objSearchEngineParams.Enzyme = "Chymotrypsin"
                        Case 3 : objSearchEngineParams.Enzyme = "Lys-C"
                        Case 4 : objSearchEngineParams.Enzyme = "Lys-N"
                        Case 5 : objSearchEngineParams.Enzyme = "Glu-C"
                        Case 6 : objSearchEngineParams.Enzyme = "Arg-C"
                        Case 7 : objSearchEngineParams.Enzyme = "Asp-N"
                        Case 8 : objSearchEngineParams.Enzyme = "alphaLP"
                        Case 9 : objSearchEngineParams.Enzyme = "no_enzyme_peptidomics"
                        Case Else : objSearchEngineParams.Enzyme = "unknown_enzyme"
                    End Select
                End If
            End If

            ' Determine the cleavage specificity
            If objSearchEngineParams.Parameters.TryGetValue("nnet", strSettingValue) Then
                ' NNET means number of non-enzymatic terminii

                If Integer.TryParse(strSettingValue, intValue) Then
                    Select Case intValue
                        Case 0
                            ' Fully-tryptic
                            objSearchEngineParams.MinNumberTermini = 2
                        Case 1
                            ' Partially-tryptic
                            objSearchEngineParams.MinNumberTermini = 1
                        Case Else
                            ' No-enzyme search
                            objSearchEngineParams.MinNumberTermini = 0
                    End Select

                End If
            Else
                ' MSGF+ uses ntt instead of nnet; thus look for ntt

                If objSearchEngineParams.Parameters.TryGetValue("ntt", strSettingValue) Then
                    ' NTT means number of tolerable terminii

                    If Integer.TryParse(strSettingValue, intValue) Then
                        Select Case intValue
                            Case 0
                                ' No-enzyme search
                                objSearchEngineParams.MinNumberTermini = 0
                            Case 1
                                ' Partially-tryptic
                                objSearchEngineParams.MinNumberTermini = 1
                            Case Else
                                ' Fully-tryptic
                                objSearchEngineParams.MinNumberTermini = 2
                        End Select

                    End If
                End If

            End If

            ' Determine the precursor mass tolerance (will store 0 if a problem or not found)
            Dim dblTolerancePPM As Double
            objSearchEngineParams.PrecursorMassToleranceDa = DeterminePrecursorMassTolerance(objSearchEngineParams, dblTolerancePPM, ePeptideHitResultType.MSGFDB)
            objSearchEngineParams.PrecursorMassTolerancePpm = dblTolerancePPM

            ' Look for Custom Amino Acid definitions
            If Not objSearchEngineParams.Parameters.Any(Function(paramEntry) paramEntry.Key = PARAM_TAG_CUSTOMAA) Then
                ' No custom amino acid entries
                Return True
            End If

            ' Store the Custom Amino Acid info
            ' Need to use a different parsing function to extract it
            success = UpdateMassCalculatorMasses(strSearchEngineParamFileName)

            ' Look for a custom charge carrier mass
            Dim customChargeCarrierMass As Double
            If GetCustomChargeCarrierMass(objSearchEngineParams, customChargeCarrierMass) Then
                ShowMessage(String.Format("Using a charge carrier mass of {0:F3} Da", customChargeCarrierMass))
                mPeptideMassCalculator.ChargeCarrierMass = customChargeCarrierMass
            End If

            Return success

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

        Dim strColumns() = strLine.Split(ControlChars.Tab)
        Dim strPeptide As String
        Dim strProtein As String

        Dim dblPrecursorMZ As Double
        Dim dblSpecProb As Double

        Dim blnMSGFPlusResults As Boolean
        Dim blnSuccess As Boolean

        objPSM = New clsPSM()

        Try

            If LookupColumnIndex(DATA_COLUMN_MSGFDB_SpecEValue, mColumnHeaders) >= 0 Then
                blnMSGFPlusResults = True
            Else
                blnMSGFPlusResults = False
            End If

            With objPSM
                .DataLineText = strLine
                .ScanNumber = LookupColumnValue(strColumns, DATA_COLUMN_Scan, mColumnHeaders, -100)
                If .ScanNumber = -100 Then
                    ' Data line is not valid
                Else
                    .ResultID = LookupColumnValue(strColumns, DATA_COLUMN_ResultID, mColumnHeaders, 0)
                    .ScoreRank = LookupColumnValue(strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb, mColumnHeaders, 1)

                    strPeptide = LookupColumnValue(strColumns, DATA_COLUMN_Peptide, mColumnHeaders)

                    If fastReadMode Then
                        .SetPeptide(strPeptide, blnUpdateCleanSequence:=False)
                    Else
                        .SetPeptide(strPeptide, mCleavageStateCalculator)
                    End If

                    .Charge = CType(LookupColumnValue(strColumns, DATA_COLUMN_Charge, mColumnHeaders, 0), Short)

                    strProtein = LookupColumnValue(strColumns, DATA_COLUMN_Protein, mColumnHeaders)
                    .AddProtein(strProtein)

                    .CollisionMode = LookupColumnValue(strColumns, DATA_COLUMN_FragMethod, mColumnHeaders, "n/a")

                    dblPrecursorMZ = LookupColumnValue(strColumns, DATA_COLUMN_PrecursorMZ, mColumnHeaders, 0.0#)
                    .PrecursorNeutralMass = mPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

                    .MassErrorDa = LookupColumnValue(strColumns, DATA_COLUMN_DelM, mColumnHeaders)
                    .MassErrorPPM = LookupColumnValue(strColumns, DATA_COLUMN_DelM_PPM, mColumnHeaders)

                    If blnMSGFPlusResults Then
                        .MSGFSpecProb = LookupColumnValue(strColumns, DATA_COLUMN_MSGFDB_SpecEValue, mColumnHeaders)
                    Else
                        .MSGFSpecProb = LookupColumnValue(strColumns, DATA_COLUMN_MSGFDB_SpecProb, mColumnHeaders)
                    End If

                    If .MSGFSpecProb.Length > 13 Then
                        ' Attempt to shorten the SpecProb value
                        If Double.TryParse(.MSGFSpecProb, dblSpecProb) Then
                            .MSGFSpecProb = dblSpecProb.ToString("0.0000000E-00")
                        End If

                    End If
                    blnSuccess = True
                End If
            End With

            If blnSuccess Then
                If Not fastReadMode Then
                    UpdatePSMUsingSeqInfo(objPSM)
                End If

                ' Store the remaining scores
                AddScore(objPSM, strColumns, DATA_COLUMN_DeNovoScore)
                AddScore(objPSM, strColumns, DATA_COLUMN_MSGFScore)

                If blnMSGFPlusResults Then

                    AddScore(objPSM, strColumns, DATA_COLUMN_MSGFDB_SpecEValue)
                    AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFDB_SpecEValue)
                    AddScore(objPSM, strColumns, DATA_COLUMN_EValue)
                    AddScore(objPSM, strColumns, DATA_COLUMN_QValue)
                    AddScore(objPSM, strColumns, DATA_COLUMN_PepQValue)
                    AddScore(objPSM, strColumns, DATA_COLUMN_Isotope_Error)

                    ' Duplicate the score values to provide backwards compatibility
                    Dim strValue = String.Empty

                    If objPSM.TryGetScore(DATA_COLUMN_MSGFDB_SpecEValue, strValue) Then objPSM.SetScore(DATA_COLUMN_MSGFDB_SpecProb, strValue)
                    If objPSM.TryGetScore(DATA_COLUMN_Rank_MSGFDB_SpecEValue, strValue) Then objPSM.SetScore(DATA_COLUMN_Rank_MSGFDB_SpecProb, strValue)
                    If objPSM.TryGetScore(DATA_COLUMN_QValue, strValue) Then objPSM.SetScore(DATA_COLUMN_FDR, strValue)
                    If objPSM.TryGetScore(DATA_COLUMN_PepQValue, strValue) Then objPSM.SetScore(DATA_COLUMN_PepFDR, strValue)

                    Dim strEValue = String.Empty
                    Dim dblEValue As Double

                    Dim strSpecEValue = String.Empty
                    Dim dblSpecEValue As Double

                    Dim blnPValueStored As Boolean = False

                    If objPSM.TryGetScore(DATA_COLUMN_EValue, strEValue) Then
                        If objPSM.TryGetScore(DATA_COLUMN_MSGFDB_SpecEValue, strSpecEValue) Then
                            ' Compute PValue using EValue and SpecEValue
                            If Double.TryParse(strEValue, dblEValue) Then
                                If Double.TryParse(strSpecEValue, dblSpecEValue) Then
                                    If dblSpecEValue > 0 Then
                                        Dim dblN As Double = dblEValue / dblSpecEValue
                                        Dim dblPValue As Double = 1 - (1 - dblSpecEValue) ^ dblN

                                        If Math.Abs(dblPValue) <= Double.Epsilon Then
                                            objPSM.SetScore(DATA_COLUMN_PValue, "0")
                                        Else
                                            objPSM.SetScore(DATA_COLUMN_PValue, dblPValue.ToString("0.00000E-00"))
                                        End If

                                        blnPValueStored = True
                                    End If
                                End If
                            End If
                        End If

                        If Not blnPValueStored Then
                            ' Store E-value as P-value (these values are not identical, and will only be close for high-confidence results, i.e. results with FDR < 2%)
                            objPSM.SetScore(DATA_COLUMN_PValue, strEValue)
                        End If

                    End If

                Else
                    AddScore(objPSM, strColumns, DATA_COLUMN_MSGFDB_SpecProb)
                    AddScore(objPSM, strColumns, DATA_COLUMN_Rank_MSGFDB_SpecProb)
                    AddScore(objPSM, strColumns, DATA_COLUMN_PValue)
                    AddScore(objPSM, strColumns, DATA_COLUMN_FDR)
                    AddScore(objPSM, strColumns, DATA_COLUMN_PepFDR)
                End If

                AddScore(objPSM, strColumns, DATA_COLUMN_EFDR)      ' This column will not be present if a Target/Decoy (TDA) search was performed

                AddScore(objPSM, strColumns, DATA_COLUMN_IMS_Scan)
                AddScore(objPSM, strColumns, DATA_COLUMN_IMS_Drift_Time)

            End If

        Catch ex As Exception
            MyBase.ReportError("Error parsing line " & intLinesRead & " in the MSGFDB data file: " & ex.Message)
        End Try

        Return blnSuccess

    End Function

    Private Function UpdateMassCalculatorMasses(strSearchEngineParamFileName As String) As Boolean

        Dim localErrorMsg As String = String.Empty
        Dim modFileProcessor = New clsMSGFPlusParamFileModExtractor("MSGF+")

        AddHandler modFileProcessor.ErrorOccurred, AddressOf ModExtractorErrorHandler
        AddHandler modFileProcessor.WarningMessageEvent, AddressOf ModExtractorWarningHandler

        Dim success = UpdateMassCalculatorMasses(strSearchEngineParamFileName, modFileProcessor, mPeptideMassCalculator, localErrorMsg)

        If Not String.IsNullOrWhiteSpace(localErrorMsg) AndAlso String.IsNullOrWhiteSpace(mErrorMessage) Then
            ReportError(localErrorMsg)
        End If

        Return success

    End Function

    ''' <summary>
    ''' Look for custom amino acid definitions in the MSGF+ parameter file
    ''' If any are found, update the amino acid mass values in the PeptideMassCalculator instance
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName"></param>
    ''' <param name="modFileProcessor"></param>
    ''' <param name="peptideMassCalculator"></param>
    ''' <param name="errorMessage"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function UpdateMassCalculatorMasses(
       strSearchEngineParamFileName As String,
       modFileProcessor As clsMSGFPlusParamFileModExtractor,
       peptideMassCalculator As clsPeptideMassCalculator,
       <Out()> ByRef errorMessage As String) As Boolean

        If modFileProcessor Is Nothing Then
            Throw New ObjectDisposedException("modFileProcessor is not initialized")
        End If

        errorMessage = String.Empty

        Dim lstModInfo As List(Of udtModInfoType) = Nothing

        ' Note that this call will initialize lstModInfo
        Dim success = modFileProcessor.ExtractModInfoFromParamFile(strSearchEngineParamFileName, lstModInfo)
        If Not success Then
            errorMessage = modFileProcessor.ErrorMessage
            Return False
        End If

        Dim customAminoAcidDefs = (From item In lstModInfo Where item.ModType = eMSGFDBModType.CustomAA Select item).ToList()
        If customAminoAcidDefs.Count = 0 Then
            ' There are no custom amino acids
            Return True
        End If

        For Each customAADef In customAminoAcidDefs

            Dim aminoAcidSymbol = customAADef.Residues(0)
            Dim empiricalFormula = customAADef.ModMass
            Dim aminoAcidMass = customAADef.ModMassVal

            Try
                Dim elementalComposition = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(empiricalFormula)

                peptideMassCalculator.SetAminoAcidMass(aminoAcidSymbol, aminoAcidMass)
                peptideMassCalculator.SetAminoAcidAtomCounts(aminoAcidSymbol, elementalComposition)

            Catch ex As Exception
                errorMessage = ex.Message
                Return False
            End Try
        Next

        Return True

    End Function

#Region "Event Handlers"
    Private Sub ModExtractorErrorHandler(errMsg As String)
        ReportError(errMsg)
    End Sub

    Private Sub ModExtractorWarningHandler(warningMsg As String)
        ReportError(warningMsg)
    End Sub
#End Region

End Class
