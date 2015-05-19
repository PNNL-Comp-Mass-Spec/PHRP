'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 05/18/2015
'
' This class parses data lines from modp_syn.txt files
'
'*********************************************************************************************************

Option Strict On

Imports System.Runtime.InteropServices
Imports System.Xml
Imports PHRPReader.clsPHRPReader

Public Class clsPHRPParserMODPlus
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
    Public Const DATA_COLUMN_NTT As String = "NTT"
    Public Const DATA_COLUMN_Modification_Annotation As String = "Modification_Annotation"
    Public Const DATA_COLUMN_Protein As String = "Protein"
    Public Const DATA_COLUMN_Peptide_Position As String = "Peptide_Position"
    Public Const DATA_COLUMN_Score As String = "Score"
    Public Const DATA_COLUMN_Probability As String = "Probability"
    Public Const DATA_COLUMN_Rank_Probability As String = "Rank_Probability"

    Public Const FILENAME_SUFFIX_SYN As String = "_modp_syn.txt"
    Public Const FILENAME_SUFFIX_FHT As String = "_modp_fht.txt"

    Protected Const MODplus_SEARCH_ENGINE_NAME As String = "MODplus"

    Protected Const PEPTIDE_MASS_TOL_PPM As String = "PeptideMassTolPPM"
    Protected Const PEPTIDE_MASS_TOL_DA As String = "PeptideMassTolDa"
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
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MODPlus, blnLoadModsAndSeqInfo)
    End Sub

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
    ''' <remarks></remarks>
    Public Sub New(ByVal strDatasetName As String, ByVal strInputFilePath As String, ByVal startupOptions As clsPHRPStartupOptions)
        MyBase.New(strDatasetName, strInputFilePath, ePeptideHitResultType.MODPlus, startupOptions)
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
        AddHeaderColumn(DATA_COLUMN_NTT)
        AddHeaderColumn(DATA_COLUMN_Modification_Annotation)
        AddHeaderColumn(DATA_COLUMN_Protein)
        AddHeaderColumn(DATA_COLUMN_Peptide_Position)
        AddHeaderColumn(DATA_COLUMN_Score)
        AddHeaderColumn(DATA_COLUMN_Probability)
        AddHeaderColumn(DATA_COLUMN_Rank_Probability)        

    End Sub

    Public Shared Function GetPHRPFirstHitsFileName(ByVal strDatasetName As String) As String
        ' MODplus does not have a first-hits file; just the _syn.txt file
        Return String.Empty
    End Function

    Public Shared Function GetPHRPModSummaryFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_modp_syn_ModSummary.txt"
    End Function

    Public Shared Function GetPHRPPepToProteinMapFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_modp_PepToProtMapMTS.txt"
    End Function

    Public Shared Function GetPHRPProteinModsFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_modp_syn_ProteinMods.txt"
    End Function

    Public Shared Function GetPHRPSynopsisFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & FILENAME_SUFFIX_SYN
    End Function

    Public Shared Function GetPHRPResultToSeqMapFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_modp_syn_ResultToSeqMap.txt"
    End Function

    Public Shared Function GetPHRPSeqInfoFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_modp_syn_SeqInfo.txt"
    End Function

    Public Shared Function GetPHRPSeqToProteinMapFileName(ByVal strDatasetName As String) As String
        Return strDatasetName & "_modp_syn_SeqToProteinMap.txt"
    End Function

    Public Shared Function GetSearchEngineName() As String
        Return MODplus_SEARCH_ENGINE_NAME
    End Function

    ''' <summary>
    ''' Parses the specified MODp parameter file
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName"></param>
    ''' <param name="objSearchEngineParams"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Overrides Function LoadSearchEngineParameters(ByVal strSearchEngineParamFileName As String, <Out()> ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Dim blnSuccess As Boolean

        objSearchEngineParams = New clsSearchEngineParameters(MODplus_SEARCH_ENGINE_NAME)

        blnSuccess = ReadSearchEngineParamFile(strSearchEngineParamFileName, objSearchEngineParams)

        ReadSearchEngineVersion(mInputFolderPath, mPeptideHitResultType, objSearchEngineParams)

        Return blnSuccess

    End Function

    Protected Function ReadSearchEngineParamFile(ByVal strSearchEngineParamFileName As String, ByVal objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Try

            ' For MSGF+ or Sequest we load mod info from the _ModDefs.txt file for the parameter file
            ' But MODplus does not have a _ModDefs.txt file because it performs a blind search
            ' The user can define static mods on any of the residues, plus the peptide terminii; check for these now

            Dim strParamFilePath = IO.Path.Combine(mInputFolderPath, strSearchEngineParamFileName)

            ' Read the contents of the parameter file
            Dim doc = New XmlDocument()
            doc.Load(strParamFilePath)

            Dim dbNodes = doc.SelectNodes("/search/database")
            If dbNodes.Count > 0 Then
                objSearchEngineParams.FastaFilePath = GetAttribute(dbNodes(0), "local_path")
            End If

            Dim enzymeNodes = doc.SelectNodes("/search/enzyme_rule")
            If enzymeNodes.Count > 0 Then
                objSearchEngineParams.Enzyme = GetAttribute(enzymeNodes(0), "name")
            End If

            Dim instrumentResolutionNodes = doc.SelectNodes("/search/instrument_resolution")
            If instrumentResolutionNodes.Count > 0 Then
                objSearchEngineParams.PrecursorMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes(0), "ms"))
                objSearchEngineParams.FragmentMassType = ConvertResolutionModeToMassType(GetAttribute(instrumentResolutionNodes(0), "msms"))
            End If

            Dim enzymeConstraintNodes = doc.SelectNodes("/search/parameters/enzyme_constraint")
            If enzymeConstraintNodes.Count > 0 Then
                Integer.TryParse(GetAttribute(enzymeConstraintNodes(0), "max_miss_cleavages"), objSearchEngineParams.MaxNumberInternalCleavages)
                Integer.TryParse(GetAttribute(enzymeConstraintNodes(0), "min_number_termini"), objSearchEngineParams.MinNumberTermini)
            End If

            Dim massTolParamNode = doc.SelectNodes("/search/parameters/peptide_mass_tol")
            If massTolParamNode.Count > 0 Then

                Dim strTolerance = GetAttribute(massTolParamNode(0), "value")
                Dim massUnits = GetAttribute(massTolParamNode(0), "unit")

                If massUnits.ToLower() = "ppm" Then
                    ' Parent mass tolerance, in ppm

                    If Double.TryParse(strTolerance, objSearchEngineParams.PrecursorMassTolerancePpm) Then
                        objSearchEngineParams.PrecursorMassToleranceDa = clsPeptideMassCalculator.PPMToMass(objSearchEngineParams.PrecursorMassTolerancePpm, 2000)
                    End If

                ElseIf massUnits.ToLower() = "da" Then

                    Double.TryParse(strTolerance, objSearchEngineParams.PrecursorMassToleranceDa)

                    ' Convert from dalton to PPM (assuming a mass of 2000 m/z)
                    objSearchEngineParams.PrecursorMassToleranceDa = clsPeptideMassCalculator.MassToPPM(objSearchEngineParams.PrecursorMassToleranceDa, 2000)
                End If

            End If

            Dim modNodes = doc.SelectNodes("/search/modifications/fixed/mod")
            If modNodes.Count > 0 Then
                ' Store the fixed mods

                For Each node As XmlNode In modNodes
                    Dim modName = GetAttribute(node, "name")
                    Dim strResidue = GetAttribute(node, "site").Trim()
                    Dim modPosition = GetAttribute(node, "position")
                    Dim modMass = GetAttribute(node, "massdiff")

                    ' Replace N-Term or C-Term with < or >
                    If strResidue.ToLower() = "n-term" Then strResidue = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS
                    If strResidue.ToLower() = "c-term" Then strResidue = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS

                    Dim modMassDa As Double
                    If Double.TryParse(modMass, modMassDa) Then

                        If Math.Abs(modMassDa - 0) > Single.Epsilon Then

                            Dim eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod
                            If strResidue = N_TERMINAL_PEPTIDE_SYMBOL_DMS OrElse strResidue = C_TERMINAL_PEPTIDE_SYMBOL_DMS Then
                                eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                            End If

                            Dim objModDef = New clsModificationDefinition(clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL, modMassDa, strResidue, eModType, "Mod" & modMassDa.ToString("0"))
                            objSearchEngineParams.AddModification(objModDef)

                        End If
                    End If

                Next

            End If

            Return True

        Catch ex As Exception
            ReportError("Error in ReadSearchEngineParamFile: " & ex.Message)
            Return False
        End Try

    End Function

    Private Function ConvertResolutionModeToMassType(resolutionType As String) As String
        If resolutionType = "high" Then
            Return clsSearchEngineParameters.MASS_TYPE_MONOISOTOPIC
        End If

        If resolutionType = "low" Then
            Return clsSearchEngineParameters.MASS_TYPE_AVERAGE
        End If

        Return "unknown"

    End Function

    Private Function GetAttribute(node As XmlNode, attributeName As String) As String

        If node.Attributes.Count > 0 Then
            Try
                Dim attribute = node.Attributes(attributeName)

                If Not attribute Is Nothing Then
                    Return attribute.Value
                End If

            Catch ex As Exception
                ' Ignore errors
                Console.WriteLine("Attribute lookup error: " & ex.Message)
            End Try
          
        End If

        Return String.Empty

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
                    .PrecursorNeutralMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .Charge, 0)

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

                AddScore(objPSM, strColumns, DATA_COLUMN_Modification_Annotation)
                AddScore(objPSM, strColumns, DATA_COLUMN_Peptide_Position)

                AddScore(objPSM, strColumns, DATA_COLUMN_Score)
                AddScore(objPSM, strColumns, DATA_COLUMN_Probability)

            End If

        Catch ex As Exception
            MyBase.ReportError("Error parsing line " & intLinesRead & " in the MODPlus data file: " & ex.Message)
        End Try

        Return blnSuccess

    End Function

End Class
