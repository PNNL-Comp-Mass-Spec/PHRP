'*********************************************************************************************************
' Written by Matthew Monroe for the US Department of Energy 
' Pacific Northwest National Laboratory, Richland, WA
'
' Created 04/04/2012
'
' This class is the base class for classes used to parse PHRP data lines
' It must be derived by a sub-class customized for the specific analysis tool (Sequest, X!Tandem, Inspect, etc.)
'
'*********************************************************************************************************

Option Strict On

Imports System.IO
Imports System.Runtime.InteropServices

Public MustInherit Class clsPHRPParser

#Region "Structures"
    Protected Structure udtAmbiguousModInfo
        Public ResidueStart As Integer
        Public ResidueEnd As Integer
        Public ModMassString As String
    End Structure
#End Region

#Region "Module variables"

    Protected mDatasetName As String
    Private mInputFilePath As String
    Protected mInputFolderPath As String
    Private mInitialized As Boolean

    Private mMaxProteinsPerPSM As Integer

    ' Column headers in the synopsis file and first hits file
    Protected mColumnHeaders As SortedDictionary(Of String, Integer)

    Protected mErrorMessage As String = String.Empty

    Protected ReadOnly mCleavageStateCalculator As clsPeptideCleavageStateCalculator
    Protected ReadOnly mPeptideMassCalculator As clsPeptideMassCalculator

    Protected mPeptideHitResultType As clsPHRPReader.ePeptideHitResultType

    Protected mModInfo As List(Of clsModificationDefinition)

    Private mResultToSeqMap As SortedList(Of Integer, Integer)
    Private mSeqInfo As SortedList(Of Integer, clsSeqInfo)
    Private mSeqToProteinMap As SortedList(Of Integer, List(Of clsProteinInfo))
    Private mPepToProteinMap As Dictionary(Of String, clsPepToProteinMapInfo)

    ' This List tracks the Protein Names for each ResultID
    Protected ReadOnly mResultIDToProteins As SortedList(Of Integer, List(Of String))

    Private ReadOnly mErrorMessages As List(Of String)
    Private ReadOnly mWarningMessages As List(Of String)

#End Region

#Region "Events"
    Public Event MessageEvent(strMessage As String)
    Public Event ErrorEvent(strErrorMessage As String)
    Public Event WarningEvent(strWarningMessage As String)
#End Region

#Region "Properties"

    ''' <summary>
    ''' Cached error messages
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ErrorMessages() As List(Of String)
        Get
            Return mErrorMessages
        End Get
    End Property

    ''' <summary>
    ''' Input file path
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property InputFilePath As String
        Get
            Return mInputFilePath
        End Get
    End Property

    ''' <summary>
    ''' Input folder path
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property InputFolderPath As String
        Get
            Return mInputFolderPath
        End Get
    End Property

    ''' <summary>
    ''' Maximum number of proteins to associate with each PSM
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks>0 means to load all proteins</remarks>
    Public Property MaxProteinsPerPSM() As Integer
        Get
            Return mMaxProteinsPerPSM
        End Get
        Set(value As Integer)
            mMaxProteinsPerPSM = value
        End Set
    End Property

    ''' <summary>
    ''' Peptide hit result type; Sequest, XTandem, Inspect, or MSGFDB
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property PeptideHitResultType As clsPHRPReader.ePeptideHitResultType
        Get
            Return mPeptideHitResultType
        End Get
    End Property

    ''' <summary>
    ''' Peptide to protein map file name
    ''' </summary>
    ''' <returns></returns>
    Public ReadOnly Property PepToProteinMap() As Dictionary(Of String, clsPepToProteinMapInfo)
        Get
            Return mPepToProteinMap
        End Get
    End Property

    ''' <summary>
    ''' Returns the cached mapping between ResultID and SeqID
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property ResultToSeqMap() As SortedList(Of Integer, Integer)
        Get
            Return mResultToSeqMap
        End Get
    End Property

    ''' <summary>
    ''' Returns the cached sequence info, where key is SeqID
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property SeqInfo() As SortedList(Of Integer, clsSeqInfo)
        Get
            Return mSeqInfo
        End Get
    End Property

    ''' <summary>
    ''' Returns the cached sequence to protein map information
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property SeqToProteinMap() As SortedList(Of Integer, List(Of clsProteinInfo))
        Get
            Return mSeqToProteinMap
        End Get
    End Property

    ''' <summary>
    ''' Cached warning messages
    ''' </summary>
    ''' <value></value>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public ReadOnly Property WarningMessages() As List(Of String)
        Get
            Return mWarningMessages
        End Get
    End Property

#End Region

#Region "Properties overridden by derived classes"

    Public MustOverride ReadOnly Property PHRPFirstHitsFileName() As String

    Public MustOverride ReadOnly Property PHRPModSummaryFileName() As String

    Public MustOverride ReadOnly Property PHRPPepToProteinMapFileName() As String

    Public MustOverride ReadOnly Property PHRPProteinModsFileName() As String

    Public MustOverride ReadOnly Property PHRPSynopsisFileName() As String

    Public MustOverride ReadOnly Property PHRPResultToSeqMapFileName() As String

    Public MustOverride ReadOnly Property PHRPSeqInfoFileName() As String

    Public MustOverride ReadOnly Property PHRPSeqToProteinMapFileName() As String

    Public MustOverride ReadOnly Property SearchEngineName() As String

#End Region

    ''' <summary>
    ''' Initialize the parser for the given dataset, input file, and result type
    ''' </summary>
    ''' <param name="strDatasetName">Dataset Name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
    ''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
    ''' <remarks>If strInputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
    Protected Sub New(strDatasetName As String, strInputFilePath As String, ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, blnLoadModsAndSeqInfo As Boolean)

        mErrorMessages = New List(Of String)
        mWarningMessages = New List(Of String)

        mResultIDToProteins = New SortedList(Of Integer, List(Of String))

        mCleavageStateCalculator = New clsPeptideCleavageStateCalculator()

        mPeptideMassCalculator = New clsPeptideMassCalculator()

        Dim startupOptions = New clsPHRPStartupOptions() With {
            .LoadModsAndSeqInfo = blnLoadModsAndSeqInfo
        }

        InitializeParser(strDatasetName, strInputFilePath, ePeptideHitResultType, startupOptions)
    End Sub

    ''' <summary>
    ''' Initialize the parser for the given dataset, input file, and result type
    ''' </summary>
    ''' <param name="strDatasetName">Dataset name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="startupOptions">Startup Options, in particular LoadModsAndSeqInfo and mMaxProteinsPerPSM</param>
    ''' <remarks>If strInputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable</remarks>
    Protected Sub New(strDatasetName As String, strInputFilePath As String, ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, startupOptions As clsPHRPStartupOptions)

        mErrorMessages = New List(Of String)
        mWarningMessages = New List(Of String)

        mResultIDToProteins = New SortedList(Of Integer, List(Of String))

        mCleavageStateCalculator = New clsPeptideCleavageStateCalculator()

        If startupOptions.PeptideMassCalculator Is Nothing Then
            mPeptideMassCalculator = New clsPeptideMassCalculator()
        Else
            mPeptideMassCalculator = startupOptions.PeptideMassCalculator
        End If

        InitializeParser(strDatasetName, strInputFilePath, ePeptideHitResultType, startupOptions)
    End Sub

    ''' <summary>
    ''' Initialize the parser for the given dataset and input file
    ''' </summary>
    ''' <param name="strDatasetName">Dataset Name</param>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
    ''' <param name="startupOptions">Startup options</param>
    ''' <remarks>If strInputFilePath is an empty string, then the functions that solely depend on dataset name will be callable, but data related functions will not be callable
    ''' startupOptions.LoadModsAndSeqInfo controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read
    ''' Setting startupOptions.MaxProteinsPerPSM to a non-zero value will limit the number of proteins that are tracked
    ''' </remarks>
    Private Sub InitializeParser(strDatasetName As String, strInputFilePath As String, ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, startupOptions As clsPHRPStartupOptions)

        If String.IsNullOrWhiteSpace(strDatasetName) Then strDatasetName = "Undefined"
        mDatasetName = strDatasetName

        mPeptideHitResultType = ePeptideHitResultType

        mMaxProteinsPerPSM = startupOptions.MaxProteinsPerPSM

        Dim fiFileInfo As FileInfo
        Dim blnIsSynopsisFile = False

        If String.IsNullOrEmpty(strInputFilePath) Then
            ' User instantiated the class without a filename
            ' Functions that solely require a dataset name will be callable, but cannot call functions that read a data line
            mInputFilePath = String.Empty
            mInputFolderPath = String.Empty

            startupOptions.LoadModsAndSeqInfo = False
        Else
            fiFileInfo = New FileInfo(strInputFilePath)
            mInputFilePath = fiFileInfo.FullName
            mInputFolderPath = fiFileInfo.DirectoryName

            Dim expectedSynopsisName = clsPHRPReader.GetPHRPSynopsisFileName(mPeptideHitResultType, mDatasetName)
            expectedSynopsisName = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(expectedSynopsisName, fiFileInfo.Name)

            If String.Equals(fiFileInfo.Name, expectedSynopsisName, StringComparison.CurrentCultureIgnoreCase) Then
                blnIsSynopsisFile = True
            End If

        End If

        mErrorMessage = String.Empty

        ' Initialize the column mapping object
        ' Using a case-insensitive comparer
        mColumnHeaders = New SortedDictionary(Of String, Integer)(StringComparer.CurrentCultureIgnoreCase)

        ' Initialize the tracking lists
        ' These will get updated via the call to objReader.GetProteinMapping
        mResultToSeqMap = New SortedList(Of Integer, Integer)
        mSeqInfo = New SortedList(Of Integer, clsSeqInfo)
        mSeqToProteinMap = New SortedList(Of Integer, List(Of clsProteinInfo))
        mPepToProteinMap = New Dictionary(Of String, clsPepToProteinMapInfo)()

        If startupOptions.LoadModsAndSeqInfo Then
            ' Read the ModSummary file (if it exists)
            LoadModSummary()
        End If

        If startupOptions.LoadModsAndSeqInfo Then
            ' Read the ResultToSeqMapInfo (if the files exist)			
            If blnIsSynopsisFile Then
                ' Assume the files exist
                LoadSeqInfo()
            Else
                ' Only continue if the fht versions exists

                Dim strResultToSeqMapFilePath As String = clsPHRPReader.GetPHRPResultToSeqMapFileName(mPeptideHitResultType, mDatasetName)
                Dim blnSeqInfoLoaded = False

                If Not String.IsNullOrEmpty(strResultToSeqMapFilePath) Then
                    strResultToSeqMapFilePath = Path.Combine(mInputFolderPath, strResultToSeqMapFilePath)
                    strResultToSeqMapFilePath = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strResultToSeqMapFilePath, mInputFilePath)
                    strResultToSeqMapFilePath = clsPHRPReader.AutoSwitchToFHTIfRequired(strResultToSeqMapFilePath, mInputFilePath)

                    If File.Exists(strResultToSeqMapFilePath) Then
                        blnSeqInfoLoaded = LoadSeqInfo()
                    End If
                End If

                If Not blnSeqInfoLoaded Then
                    If String.IsNullOrEmpty(strResultToSeqMapFilePath) Then
                        ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file and unable to determine the ResultToSeqMapFilename using clsPHRPReader.GetPHRPResultToSeqMapFileName()")
                    Else
                        ReportWarning("Unable to load data from the SeqInfo files since reading a first-hits file but the ResultToSeqMap file does not exist: " & strResultToSeqMapFilePath)
                    End If
                End If

            End If
        End If

        ' The following will be overridden by a derived form of this class
        DefineColumnHeaders()

        mInitialized = True

    End Sub

    ''' <summary>
    ''' Returns the appropriate PHRPParser class based on the input file name; assumes blnLoadModsAndSeqInfo=True
    ''' </summary>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from strInputFilePath</remarks>
    Public Shared Function GetParser(strInputFilePath As String) As clsPHRPParser
        Return GetParser(strInputFilePath, True)
    End Function

    ''' <summary>
    ''' Returns the appropriate PHRPParser class based on the input file name
    ''' </summary>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
    ''' <remarks>Throws an exception if unable to auto-determine the input file type or dataset name from strInputFilePath</remarks>
    Public Shared Function GetParser(strInputFilePath As String, blnLoadModsAndSeqInfo As Boolean) As clsPHRPParser
        Dim ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath)

        If ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Unknown Then
            Throw New Exception("Unable to auto-determine the PeptideHitResultType for " & strInputFilePath)
        End If

        Dim strDatasetName = clsPHRPReader.AutoDetermineDatasetName(strInputFilePath)
        If String.IsNullOrEmpty(strDatasetName) Then
            Throw New Exception("Unable to auto-determine the Dataset Name for " & strInputFilePath)
        End If

        Return GetParser(strInputFilePath, strDatasetName, ePeptideHitResultType, blnLoadModsAndSeqInfo)
    End Function

    ''' <summary>
    ''' Returns the appropriate PHRPParser class based on the input file name
    ''' </summary>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' ''' <param name="strDatasetName">Dataset Name</param>
    ''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
    ''' <remarks>Throws an exception if unable to auto-determine the input file type from strInputFilePath</remarks>
    Public Shared Function GetParser(strInputFilePath As String, strDatasetName As String, blnLoadModsAndSeqInfo As Boolean) As clsPHRPParser
        Dim ePeptideHitResultType = clsPHRPReader.AutoDetermineResultType(strInputFilePath)

        If ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.Unknown Then
            Throw New Exception("Unable to auto-determine the PeptideHitResultType for " & strInputFilePath)
        End If

        Return GetParser(strInputFilePath, strDatasetName, ePeptideHitResultType, blnLoadModsAndSeqInfo)
    End Function

    ''' <summary>
    ''' Returns the appropriate PHRPParser class based on ePeptideHitResultType
    ''' </summary>
    ''' <param name="strInputFilePath">Input file path</param>
    ''' <param name="strDatasetName">Dataset Name</param>
    ''' <param name="ePeptideHitResultType">Peptide Hit Results file type</param>
    ''' <param name="blnLoadModsAndSeqInfo">Controls whether or not the _SeqInfo.txt and _SeqToProteinMap.txt files should be read</param>
    ''' <remarks></remarks>
    Public Shared Function GetParser(strInputFilePath As String, strDatasetName As String, ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType, blnLoadModsAndSeqInfo As Boolean) As clsPHRPParser
        Select Case ePeptideHitResultType
            Case clsPHRPReader.ePeptideHitResultType.Inspect
                Return New clsPHRPParserInspect(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case clsPHRPReader.ePeptideHitResultType.MSAlign
                Return New clsPHRPParserMSAlign(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case clsPHRPReader.ePeptideHitResultType.MSGFDB
                Return New clsPHRPParserMSGFDB(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case clsPHRPReader.ePeptideHitResultType.Sequest
                Return New clsPHRPParserSequest(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case clsPHRPReader.ePeptideHitResultType.XTandem
                Return New clsPHRPParserXTandem(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case clsPHRPReader.ePeptideHitResultType.MODa
                Return New clsPHRPParserMODa(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case clsPHRPReader.ePeptideHitResultType.MODPlus
                Return New clsPHRPParserMODPlus(strDatasetName, strInputFilePath, blnLoadModsAndSeqInfo)

            Case Else
                Throw New Exception("Unrecognized value for PeptideHitResultType: " & ePeptideHitResultType.ToString())
        End Select

    End Function

#Region "Functions overridden by derived classes"
    Protected MustOverride Sub DefineColumnHeaders()

    ''' <summary>
    ''' Parse the data line read from a PHRP results file
    ''' </summary>
    ''' <param name="strLine">Data line</param>
    ''' <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
    ''' <param name="objPSM">clsPSM object (output)</param>
    ''' <returns>True if success, false if an error</returns>
    Public Function ParsePHRPDataLine(strLine As String, intLinesRead As Integer, <Out()> ByRef objPSM As clsPSM) As Boolean
        Return ParsePHRPDataLine(strLine, intLinesRead, objPSM, fastReadMode:=False)
    End Function

    ''' <summary>
    ''' Parse the data line read from a PHRP results file
    ''' </summary>
    ''' <param name="strLine">Data line</param>
    ''' <param name="intLinesRead">Number of lines read so far (used for error reporting)</param>
    ''' <param name="objPSM">clsPSM object (output)</param>
    ''' <param name="fastReadMode">When set to true, then reads the next data line, but doesn't perform text parsing required to determine cleavage state</param>
    ''' <returns>True if success, false if an error</returns>
    ''' <remarks>When fastReadMode is True, you should call FinalizePSM to populate the remaining fields if the peptide is a peptide of interest</remarks>
    Public MustOverride Function ParsePHRPDataLine(strLine As String, intLinesRead As Integer, <Out()> ByRef objPSM As clsPSM, fastReadMode As Boolean) As Boolean

    ''' <summary>
    ''' Parses the specified parameter file
    ''' Also reads the Tool_Version_Info file in the same folder (if present)
    ''' </summary>
    ''' <param name="strSearchEngineParamFileName">Name of the parameter file to parse (must reside in InputFolderPath)</param>
    ''' <param name="objSearchEngineParams">Search engine parameters class (output)</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public MustOverride Function LoadSearchEngineParameters(strSearchEngineParamFileName As String, <Out()> ByRef objSearchEngineParams As clsSearchEngineParameters) As Boolean

#End Region

    Protected Sub AddHeaderColumn(strColumnName As String)
        mColumnHeaders.Add(strColumnName, mColumnHeaders.Count)
    End Sub

    Protected Sub AddScore(objPSM As clsPSM, strColumns() As String, strScoreColumnName As String)
        Const NOT_FOUND = "==SCORE_NOT_FOUND=="

        Dim strValue As String
        strValue = clsPHRPReader.LookupColumnValue(strColumns, strScoreColumnName, mColumnHeaders, NOT_FOUND)

        If strValue <> NOT_FOUND Then
            objPSM.SetScore(strScoreColumnName, strValue)
        End If

    End Sub

    ''' <summary>
    ''' Clear any cached error messages
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearErrors()
        mErrorMessages.Clear()
    End Sub

    ''' <summary>
    ''' Clear any cached warning messages
    ''' </summary>
    ''' <remarks></remarks>
    Public Sub ClearWarnings()
        mWarningMessages.Clear()
    End Sub

    Private Function ConvertModsToNumericMods(strCleanSequence As String, lstModifiedResidues As List(Of clsAminoAcidModInfo)) As String
        Static sbNewPeptide As New Text.StringBuilder

        sbNewPeptide.Length = 0

        If lstModifiedResidues Is Nothing OrElse lstModifiedResidues.Count = 0 Then
            Return strCleanSequence
        End If

        For intIndex = 0 To strCleanSequence.Length - 1
            sbNewPeptide.Append(strCleanSequence.Chars(intIndex))

            For Each objModInfo In lstModifiedResidues
                If objModInfo.ResidueLocInPeptide = intIndex + 1 Then
                    sbNewPeptide.Append(NumToStringPlusMinus(objModInfo.ModDefinition.ModificationMass, 4))
                End If
            Next
        Next

        Return sbNewPeptide.ToString()

    End Function

    ''' <summary>
    ''' Look for ambiguous mods in strSequenceWithMods
    ''' For example, -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q
    ''' </summary>
    ''' <param name="strSequenceWithMods"></param>
    ''' <returns></returns>
    ''' <remarks>List of ambiguous mods, where the keys are the start residues and the values are the ambiguous mod info</remarks>
    Private Function ExtractAmbiguousMods(strSequenceWithMods As String) As SortedList(Of Integer, udtAmbiguousModInfo)

        Dim strPrimarySequence As String = String.Empty
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        If Not clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strSequenceWithMods, strPrimarySequence, strPrefix, strSuffix) Then
            strPrimarySequence = String.Copy(strSequenceWithMods)
        End If

        Dim lstAmbiguousMods = New SortedList(Of Integer, udtAmbiguousModInfo)

        Dim intResidueNumber = 0
        Dim blnParsingAmbiguousMod = False
        Dim udtCurrentMod As udtAmbiguousModInfo

        For intCharIndex = 0 To strPrimarySequence.Length - 1
            If clsPHRPReader.IsLetterAtoZ(strPrimarySequence.Chars(intCharIndex)) Then
                ' Found a letter
                intResidueNumber += 1

                If intCharIndex > 0 AndAlso strPrimarySequence.Chars(intCharIndex - 1) = "("c Then
                    ' Found an ambiguous mod
                    If Not blnParsingAmbiguousMod Then
                        blnParsingAmbiguousMod = True
                        udtCurrentMod.ResidueStart = intResidueNumber
                        udtCurrentMod.ResidueEnd = intResidueNumber
                        udtCurrentMod.ModMassString = String.Empty
                    End If
                End If

            ElseIf blnParsingAmbiguousMod Then
                ' Found a non-letter, and we are parsing an ambiguous mod

                udtCurrentMod.ResidueEnd = intResidueNumber
                blnParsingAmbiguousMod = False

                ' The mod mass should be next, in the form [-30.09]
                ' Parse out the mod mass
                If intCharIndex < strPrimarySequence.Length - 2 Then
                    If strPrimarySequence.Chars(intCharIndex + 1) = "[" Then
                        Dim strModMassString = strPrimarySequence.Substring(intCharIndex + 2)
                        Dim bracketIndex = strModMassString.IndexOf("]"c)
                        If bracketIndex > 0 Then
                            ' Valid ambiguous mod found; store it
                            strModMassString = strModMassString.Substring(0, bracketIndex)
                            udtCurrentMod.ModMassString = String.Copy(strModMassString)

                            lstAmbiguousMods.Add(udtCurrentMod.ResidueStart, udtCurrentMod)
                        End If
                    End If
                End If
            End If
        Next

        Return lstAmbiguousMods

    End Function

    Public Sub FinalizePSM(objPSM As clsPSM)

        objPSM.UpdateCleanSequence()

        objPSM.UpdateCleavageInfo(mCleavageStateCalculator)

        UpdatePSMUsingSeqInfo(objPSM)

    End Sub

    Private Shared Function GetMODaStaticModSetting(kvSetting As KeyValuePair(Of String, String), <Out()> warningMessage As String) As KeyValuePair(Of String, String)

        Dim strKey = kvSetting.Key
        Dim strValue = kvSetting.Value

        If Not String.Equals(strKey, "ADD", StringComparison.CurrentCultureIgnoreCase) Then
            Throw New Exception("Key name is not ADD; this is not a MODa Static Mod Setting")
        End If

        Dim commaIndex = strValue.IndexOf(","c)

        If commaIndex > 0 Then
            Dim strResidue = strValue.Substring(0, commaIndex).Trim()
            strValue = strValue.Substring(commaIndex + 1).Trim()

            ' Update Key to look like ADD_A or ADD_NTerm
            strKey = strKey & "_" & strResidue

            kvSetting = New KeyValuePair(Of String, String)(strKey, strValue)
            warningMessage = String.Empty
        Else
            warningMessage = "Value for MODa keyword ADD does not contain a comma"
        End If

        Return kvSetting

    End Function

    Protected Sub HandleException(strBaseMessage As String, ex As Exception)
        If String.IsNullOrEmpty(strBaseMessage) Then
            strBaseMessage = "Error"
        End If

        ReportError(strBaseMessage & ": " & ex.Message)
    End Sub

    ''' <summary>
    ''' Examines the string to determine if it is numeric
    ''' </summary>
    ''' <param name="strData"></param>
    ''' <returns>True if a number, otherwise false</returns>
    Public Shared Function IsNumber(strData As String) As Boolean
        Try
            If Double.TryParse(strData, 0) Then
                Return True
            ElseIf Integer.TryParse(strData, 0) Then
                Return True
            End If
        Catch ex As Exception
            ' Ignore errors here
        End Try

        Return False
    End Function

    ''' <summary>
    ''' Reads the data in strModSummaryFilePath.  Populates mModInfo with the modification names, masses, and affected residues
    ''' </summary>
    ''' <returns>True if success; false if an error</returns>
    Private Function LoadModSummary() As Boolean

        Dim objModSummaryReader As clsPHRPModSummaryReader

        Dim strModSummaryFilePath As String
        Dim strModSummaryFilePathPreferred As String

        Dim blnSuccess As Boolean

        Try
            strModSummaryFilePath = clsPHRPReader.GetPHRPModSummaryFileName(mPeptideHitResultType, mDatasetName)
            If String.IsNullOrEmpty(strModSummaryFilePath) Then
                ReportWarning("ModSummaryFile path is empty; unable to continue")
                Return False
            End If

            strModSummaryFilePath = Path.Combine(mInputFolderPath, strModSummaryFilePath)

            strModSummaryFilePath = clsPHRPReader.AutoSwitchToLegacyMSGFDBIfRequired(strModSummaryFilePath, mInputFilePath)
            strModSummaryFilePathPreferred = clsPHRPReader.AutoSwitchToFHTIfRequired(strModSummaryFilePath, mInputFilePath)
            If strModSummaryFilePath <> strModSummaryFilePathPreferred AndAlso File.Exists(strModSummaryFilePathPreferred) Then
                strModSummaryFilePath = strModSummaryFilePathPreferred
            End If

            If Not File.Exists(strModSummaryFilePath) Then
                ReportWarning("ModSummary file not found: " & strModSummaryFilePath)
                Return False
            End If

            ShowMessage("Reading the PHRP ModSummary file")

            objModSummaryReader = New clsPHRPModSummaryReader(strModSummaryFilePath)
            blnSuccess = objModSummaryReader.Success

            If blnSuccess Then
                mModInfo = objModSummaryReader.ModificationDefs
            End If

        Catch ex As Exception
            HandleException("Exception reading PHRP Mod Summary file", ex)
            Return False
        End Try

        Return blnSuccess

    End Function

    Private Function LoadSeqInfo() As Boolean

        Dim blnSuccess As Boolean
        Dim objReader As clsPHRPSeqMapReader

        Dim entriesParsed As Integer
        Dim dtLastProgress As DateTime = DateTime.UtcNow()
        Dim blnNotifyComplete As Boolean

        Try

            ShowMessage("Reading the PHRP SeqInfo file")

            ' Instantiate the reader
            objReader = New clsPHRPSeqMapReader(mDatasetName, mInputFolderPath, mPeptideHitResultType, mInputFilePath)

            objReader.MaxProteinsPerSeqID = mMaxProteinsPerPSM

            ' Read the files
            blnSuccess = objReader.GetProteinMapping(mResultToSeqMap, mSeqToProteinMap, mSeqInfo, mPepToProteinMap)

            If Not blnSuccess Then
                ReportWarning(objReader.ErrorMessage)
            End If

            mResultIDToProteins.Clear()

            If blnSuccess Then
                ' Populate mResultIDToProteins

                For Each objItem As KeyValuePair(Of Integer, Integer) In mResultToSeqMap

                    Dim lstProteinsForSeqID As List(Of clsProteinInfo) = Nothing
                    Dim lstProteinsForResultID As List(Of String)

                    If mSeqToProteinMap.TryGetValue(objItem.Value, lstProteinsForSeqID) Then
                        Const USE_LINQ = True
                        If USE_LINQ Then
                            lstProteinsForResultID = (From objProtein In lstProteinsForSeqID Select objProtein.ProteinName Distinct).ToList()
                        Else

                            lstProteinsForResultID = New List(Of String)(lstProteinsForSeqID.Count)
                            For Each objProtein As clsProteinInfo In lstProteinsForSeqID
                                If Not lstProteinsForResultID.Contains(objProtein.ProteinName) Then
                                    lstProteinsForResultID.Add(objProtein.ProteinName)
                                End If
                            Next
                        End If
                    Else
                        lstProteinsForResultID = New List(Of String)
                    End If

                    If mMaxProteinsPerPSM > 0 AndAlso lstProteinsForResultID.Count > mMaxProteinsPerPSM Then
                        ' Only add a subset of the proteins in lstProteinsForResultID
                        Dim lstProteinSubset = (From item In lstProteinsForResultID Take mMaxProteinsPerPSM Order By item Select item).ToList()
                        mResultIDToProteins.Add(objItem.Key, lstProteinSubset)
                    Else
                        mResultIDToProteins.Add(objItem.Key, lstProteinsForResultID)
                    End If

                    entriesParsed += 1
                    If DateTime.UtcNow.Subtract(dtLastProgress).TotalSeconds >= 5 Then
                        Dim pctComplete = entriesParsed / CDbl(mResultToSeqMap.Count) * 100
                        Console.WriteLine(" ... associating proteins with sequences: " & pctComplete.ToString("0.0") & "% complete")
                        dtLastProgress = DateTime.UtcNow
                        blnNotifyComplete = True
                    End If

                Next

                If blnNotifyComplete Then
                    Console.WriteLine(" ... associating proteins with sequences: 100% complete")
                End If

            End If

        Catch ex As Exception
            HandleException("Error loading PHRP Seq Info", ex)
            blnSuccess = False
            If Not mInitialized Then Throw New Exception(mErrorMessage, ex)
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' Formats a number so that it begins with a + sign if positive or a - sign if negative
    ''' Rounds the number to the specified number of digits, trimming off trailing zeros
    ''' Example output: +79.9663 or -17.016
    ''' </summary>
    ''' <param name="Value"></param>
    ''' <param name="DigitsOfPrecision"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Public Shared Function NumToStringPlusMinus(Value As Double, DigitsOfPrecision As Integer) As String

        Dim strFormatString = "+0;-0"
        If DigitsOfPrecision > 0 Then
            strFormatString = "+0." & New String("0"c, DigitsOfPrecision) & ";-0." & New String("0"c, DigitsOfPrecision)
        End If

        Dim strValue As String = Value.ToString(strFormatString).TrimEnd("0"c)

        If strValue.EndsWith("."c) Then
            ' Ends in a decimal point; remove the decimal point
            strValue = strValue.TrimEnd("."c)
        End If

        Return strValue

    End Function

    ''' <summary>
    ''' Parse the column names in strSplitLine and update the local column header mapping
    ''' </summary>
    ''' <param name="strSplitLine"></param>
    ''' <remarks></remarks>
    Public Sub ParseColumnHeaders(strSplitLine() As String)
        clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders)
    End Sub

    ''' <summary>
    ''' Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
    ''' </summary>
    ''' <param name="strText"></param>
    ''' <param name="chDelimiter"></param>
    ''' <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
    ''' <remarks>Automatically trims whitespace</remarks>
    Public Shared Function ParseKeyValueSetting(strText As String, chDelimiter As Char) As KeyValuePair(Of String, String)
        Return ParseKeyValueSetting(strText, chDelimiter, String.Empty)
    End Function

    ''' <summary>
    ''' Splits strText on strText, returning a KeyValuePair object where the key is the text to the left of the delimiter and the value is the text to the right
    ''' </summary>
    ''' <param name="strText"></param>
    ''' <param name="chDelimiter"></param>
    ''' <param name="strCommentChar">If defined, then looks for this character in the value portion of the setting and removes that character plus any text after it</param>
    ''' <returns>KeyValuePair with key and value from strText; key and value will be empty if chDelimiter was not found</returns>
    ''' <remarks>Automatically trims whitespace</remarks>
    Public Shared Function ParseKeyValueSetting(strText As String, chDelimiter As Char, strCommentChar As String) As KeyValuePair(Of String, String)
        Dim kvSetting As KeyValuePair(Of String, String)
        Dim strKey As String
        Dim strValue As String
        Dim intCharIndex As Integer

        If Not String.IsNullOrEmpty(strText) Then
            intCharIndex = strText.IndexOf(chDelimiter)
            If intCharIndex > 0 Then
                strKey = strText.Substring(0, intCharIndex).Trim()
                If intCharIndex < strText.Length - 1 Then
                    strValue = strText.Substring(intCharIndex + 1).Trim()

                    If Not String.IsNullOrEmpty(strCommentChar) Then
                        ' Look for the comment character
                        Dim commentCharIndex = strValue.IndexOf(strCommentChar, StringComparison.Ordinal)
                        If commentCharIndex > 0 Then
                            ' Trim off the comment
                            strValue = strValue.Substring(0, commentCharIndex).Trim()
                        End If
                    End If

                Else
                    strValue = String.Empty
                End If

                kvSetting = New KeyValuePair(Of String, String)(strKey, strValue)
                Return kvSetting
            End If
        End If

        Return New KeyValuePair(Of String, String)(String.Empty, String.Empty)

    End Function

    ''' <summary>
    ''' Read a Search Engine parameter file where settings are stored as key/value pairs
    ''' </summary>
    ''' <param name="strSearchEngineName">Search engine name (e.g. MSGF+)</param>
    ''' <param name="strSearchEngineParamFileName">Search engine parameter file name (must exist in mInputFolderPath)</param>
    ''' <param name="ePeptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
    ''' <param name="objSearchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
    ''' <returns>True if success, false if an error</returns>
    Protected Function ReadKeyValuePairSearchEngineParamFile(
      strSearchEngineName As String,
      strSearchEngineParamFileName As String,
      ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType,
      objSearchEngineParams As clsSearchEngineParameters) As Boolean

        Dim paramFilePath = Path.Combine(mInputFolderPath, strSearchEngineParamFileName)

        Dim errorMessage As String = Nothing
        Dim warningMessage As String = Nothing
        Dim success = ReadKeyValuePairSearchEngineParamFile(
            strSearchEngineName, paramFilePath, ePeptideHitResultType,
            objSearchEngineParams, errorMessage, warningMessage)

        If Not String.IsNullOrWhiteSpace(errorMessage) Then
            ReportError(errorMessage)
        End If

        If Not String.IsNullOrWhiteSpace(warningMessage) Then
            ReportWarning(warningMessage)
        End If

        Return success

    End Function

    ''' <summary>
    ''' Read a Search Engine parameter file where settings are stored as key/value pairs
    ''' </summary>
    ''' <param name="searchEngineName">Search engine name (e.g. MS-GF+)</param>
    ''' <param name="paramFilePath">Search engine parameter file path</param>
    ''' <param name="ePeptideHitResultType">PeptideHitResultType (only important if reading a ModA parameter file</param>
    ''' <param name="searchEngineParams">SearchEngineParams container class (must be initialized by the calling function)</param>
    ''' <param name="errorMessage">Output: error message</param>
    ''' <param name="warningMessage">Output: warning message</param>
    ''' <returns>True if success, false if an error</returns>
    Public Shared Function ReadKeyValuePairSearchEngineParamFile(
      searchEngineName As String,
      paramFilePath As String,
      ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType,
      searchEngineParams As clsSearchEngineParameters,
      <Out()> errorMessage As String,
      <Out()> warningMessage As String) As Boolean

        errorMessage = String.Empty
        warningMessage = String.Empty

        Try
            If String.IsNullOrWhiteSpace(searchEngineName) Then searchEngineName = "?? Unknown tool ??"

            If Not File.Exists(paramFilePath) Then
                errorMessage = searchEngineName & " param file not found: " & paramFilePath
                Return False
            End If

            searchEngineParams.UpdateSearchEngineParamFilePath(paramFilePath)

            Using srInFile = New StreamReader(New FileStream(paramFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                While Not srInFile.EndOfStream
                    Dim strLineIn = srInFile.ReadLine().TrimStart()

                    If String.IsNullOrWhiteSpace(strLineIn) OrElse strLineIn.StartsWith("#") OrElse Not strLineIn.Contains("="c) Then
                        Continue While
                    End If


                    ' Split the line on the equals sign
                    Dim kvSetting = ParseKeyValueSetting(strLineIn, "="c, "#")

                    If String.IsNullOrEmpty(kvSetting.Key) Then
                        Continue While
                    End If

                    If ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MODa AndAlso
                       String.Equals(kvSetting.Key, "add", StringComparison.CurrentCultureIgnoreCase) Then
                        ' ModA defines all of its static modifications with the ADD keyword
                        ' Split the value at the comma and create a new setting entry with the residue name
                        Dim warningMessageAddon As String = Nothing
                        kvSetting = GetMODaStaticModSetting(kvSetting, warningMessageAddon)
                        If Not String.IsNullOrWhiteSpace(warningMessageAddon) Then
                            If String.IsNullOrWhiteSpace(warningMessage) Then
                                warningMessage = String.Copy(warningMessageAddon)
                            ElseIf Not warningMessage.Contains(warningMessageAddon) Then
                                warningMessage &= "; " & String.Copy(warningMessageAddon)
                            End If
                        End If
                    End If

                    searchEngineParams.AddUpdateParameter(kvSetting)

                End While
            End Using

            Return True

        Catch ex As Exception
            errorMessage = String.Format("Error in ReadKeyValuePairSearchEngineParamFile for {0}, param file {1}: {2}",
                                         searchEngineName, Path.GetFileName(paramFilePath), ex.Message)
            Return False
        End Try

    End Function

    Protected Function ReadSearchEngineVersion(
      ePeptideHitResultType As clsPHRPReader.ePeptideHitResultType,
      objSearchEngineParams As clsSearchEngineParameters) As Boolean
        Dim blnSuccess = False

        Try
            ' Read the Tool_Version_Info file to determine the analysis time and the tool version
            Dim toolVersionInfoFilePath = Path.Combine(mInputFolderPath, clsPHRPReader.GetToolVersionInfoFilename(ePeptideHitResultType))

            If Not File.Exists(toolVersionInfoFilePath) AndAlso ePeptideHitResultType = clsPHRPReader.ePeptideHitResultType.MSGFDB Then
                ' This could be an older MSGF+ job; check for a _MSGFDB.txt tool version file
                Dim alternativeVersionInfoFilePath = Path.Combine(mInputFolderPath, "Tool_Version_Info_MSGFDB.txt")
                If File.Exists(alternativeVersionInfoFilePath) Then
                    toolVersionInfoFilePath = alternativeVersionInfoFilePath
                End If
            End If

            If Not File.Exists(toolVersionInfoFilePath) Then
                ReportWarning("Tool version info file not found: " & toolVersionInfoFilePath)
                Return False
            End If

            Dim strSearchEngineVersion = "Unknown"
            Dim dtSearchDate = New DateTime(1980, 1, 1)
            Dim blnValidDate = False
            Dim blnValidVersion = False

            Using srInFile = New StreamReader(New FileStream(toolVersionInfoFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                While Not srInFile.EndOfStream
                    Dim strLineIn = srInFile.ReadLine().TrimStart()

                    ' Split the line on a colon
                    Dim kvSetting = ParseKeyValueSetting(strLineIn, ":"c)

                    Select Case kvSetting.Key.ToLower()
                        Case "date"
                            blnValidDate = DateTime.TryParse(kvSetting.Value, dtSearchDate)

                        Case "toolversioninfo"
                            If Not String.IsNullOrEmpty(kvSetting.Value) Then
                                strSearchEngineVersion = String.Copy(kvSetting.Value)
                                blnValidVersion = True
                            Else
                                ' The next line contains the search engine version
                                If Not srInFile.EndOfStream Then
                                    strLineIn = srInFile.ReadLine().TrimStart()
                                    strSearchEngineVersion = String.Copy(strLineIn)
                                    blnValidVersion = True
                                End If
                            End If
                        Case Else
                            ' Ignore the line
                    End Select
                End While

            End Using

            If Not blnValidDate Then
                ReportError("Date line not found in the ToolVersionInfo file")
                blnSuccess = False
            ElseIf Not blnValidVersion Then
                ReportError("ToolVersionInfo line not found in the ToolVersionInfo file")
                blnSuccess = False
            Else
                blnSuccess = True
            End If

            objSearchEngineParams.UpdateSearchEngineVersion(strSearchEngineVersion)
            objSearchEngineParams.UpdateSearchDate(dtSearchDate)

        Catch ex As Exception
            ReportError("Error in ReadSearchEngineVersion: " & ex.Message)
        End Try

        Return blnSuccess

    End Function

    Protected Sub ReportError(strErrorMessage As String)
        mErrorMessage = strErrorMessage
        mErrorMessages.Add(strErrorMessage)
        RaiseEvent ErrorEvent(strErrorMessage)
    End Sub

    Protected Sub ReportWarning(strWarningMessage As String)
        mWarningMessages.Add(strWarningMessage)
        RaiseEvent WarningEvent(strWarningMessage)
    End Sub

    Protected Sub ShowMessage(strMessage As String)
        RaiseEvent MessageEvent(strMessage)
    End Sub

    Private Sub StoreModInfo(objPSM As clsPSM, objSeqInfo As clsSeqInfo)

        Dim strMods() As String
        Dim kvModDetails As KeyValuePair(Of String, String)

        Dim strMassCorrectionTag As String
        Dim intResidueLoc As Integer

        Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants

        Dim blnMatchFound As Boolean

        Dim blnFavorTerminalMods As Boolean

        Dim lstNTerminalModsAdded = New List(Of String)
        Dim lstCTerminalModsAdded = New List(Of String)
        Dim intPeptideResidueCount = objPSM.PeptideCleanSequence.Length

        objPSM.PeptideMonoisotopicMass = objSeqInfo.MonoisotopicMass

        objPSM.ClearModifiedResidues()

        If objSeqInfo.ModCount > 0 Then
            ' Split objSeqInfo.ModDescription on the comma character
            strMods = objSeqInfo.ModDescription.Split(","c)

            If Not strMods Is Nothing AndAlso strMods.Count > 0 Then

                ' Parse objPSM.Peptide to look for ambiguous mods, for example -30.09 in I.(TIIQ)[-30.09]APQGVSLQYTSR.Q

                Dim lstAmbiguousMods = ExtractAmbiguousMods(objPSM.Peptide)

                ' ReSharper disable once UseImplicitlyTypedVariableEvident
                For intModIndex As Integer = 0 To strMods.Count - 1

                    ' Split strMods on the colon characters
                    kvModDetails = ParseKeyValueSetting(strMods(intModIndex), ":"c)

                    If Not String.IsNullOrEmpty(kvModDetails.Key) AndAlso Not String.IsNullOrEmpty(kvModDetails.Value) Then
                        strMassCorrectionTag = kvModDetails.Key
                        If Integer.TryParse(kvModDetails.Value, intResidueLoc) Then
                            ' Find the modification definition in mModInfo
                            ' Note that a given mass correction tag might be present multiple times in mModInfo, since it could be used as both a static peptide mod and a static peptide terminal mod
                            ' Thus, if intResidueLoc = 1 or intResidueLoc = objPSM.PeptideCleanSequence.Length then we'll first look for a peptide or protein terminal static mod

                            If intResidueLoc = 1 Then
                                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
                                If lstNTerminalModsAdded.Contains(strMassCorrectionTag) Then
                                    ' We have likely already added this modification as an N-terminal mod, thus, don't favor terminal mods this time
                                    ' An example is an iTraq peptide where there is a K at the N-terminus
                                    ' It gets modified with iTraq twice: once because of the N-terminus and once because of Lysine
                                    ' For example, R.K+144.102063+144.102063TGSY+79.9663GALAEITASK+144.102063.E
                                    blnFavorTerminalMods = False
                                Else
                                    blnFavorTerminalMods = True
                                End If

                            ElseIf intResidueLoc = intPeptideResidueCount Then
                                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
                                If lstCTerminalModsAdded.Contains(strMassCorrectionTag) Then
                                    blnFavorTerminalMods = False
                                Else
                                    blnFavorTerminalMods = True
                                End If

                            Else
                                eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
                                blnFavorTerminalMods = False
                            End If

                            Dim objMatchedModDef As clsModificationDefinition = Nothing

                            If mModInfo Is Nothing Then
                                objMatchedModDef = New clsModificationDefinition()
                                objMatchedModDef.MassCorrectionTag = strMassCorrectionTag
                                blnMatchFound = True
                            Else
                                blnMatchFound = UpdatePSMFindMatchingModInfo(strMassCorrectionTag, blnFavorTerminalMods, eResidueTerminusState, objMatchedModDef)
                            End If

                            If blnMatchFound Then
                                Dim lstMatches = (From item In lstAmbiguousMods Where item.Key = intResidueLoc Select item.Value).ToList()

                                If lstMatches.Count > 0 Then
                                    ' Ambiguous modification
                                    objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence.Chars(intResidueLoc - 1), intResidueLoc, eResidueTerminusState, objMatchedModDef, lstMatches.First.ResidueEnd)
                                Else
                                    ' Normal, non-ambiguous modified residue
                                    objPSM.AddModifiedResidue(objPSM.PeptideCleanSequence.Chars(intResidueLoc - 1), intResidueLoc, eResidueTerminusState, objMatchedModDef)
                                End If

                                If intResidueLoc = 1 Then
                                    lstNTerminalModsAdded.Add(objMatchedModDef.MassCorrectionTag)
                                ElseIf intResidueLoc = intPeptideResidueCount Then
                                    lstCTerminalModsAdded.Add(objMatchedModDef.MassCorrectionTag)
                                End If
                            Else
                                ' Could not find a valid entry in mModInfo
                                ReportError("Unrecognized mass correction tag found in the SeqInfo file: " & strMassCorrectionTag)
                            End If

                        End If
                    End If
                Next intModIndex
            End If

        End If

    End Sub

    ''' <summary>
    ''' Updates the theoretical (computed) monoisotopic mass of objPSM using mResultToSeqMap and mSeqInfo
    ''' Also updates the modification info
    ''' Also updates SeqID
    ''' </summary>
    ''' <param name="objPSM"></param>
    ''' <returns>True if success, False if objPSM.ResultID is not found in mResultToSeqMap</returns>
    ''' <remarks></remarks>
    Protected Function UpdatePSMUsingSeqInfo(objPSM As clsPSM) As Boolean
        Dim intSeqID As Integer
        Dim objSeqInfo As clsSeqInfo = Nothing

        Dim blnSuccess As Boolean

        blnSuccess = False

        ' First determine the modified residues present in this peptide
        If Not mResultToSeqMap Is Nothing AndAlso mResultToSeqMap.Count > 0 Then
            If mResultToSeqMap.TryGetValue(objPSM.ResultID, intSeqID) Then

                objPSM.SeqID = intSeqID

                If mSeqInfo.TryGetValue(intSeqID, objSeqInfo) Then
                    StoreModInfo(objPSM, objSeqInfo)
                    blnSuccess = True
                End If

                ' Lookup the protein details using mSeqToProteinMap
                Dim lstProteinDetails As List(Of clsProteinInfo) = Nothing
                If mSeqToProteinMap.TryGetValue(intSeqID, lstProteinDetails) Then
                    For Each oProtein In lstProteinDetails
                        objPSM.AddProteinDetail(oProtein)
                    Next
                End If

                ' Make sure all of the proteins in objPSM.Proteins are defined in objPSM.ProteinDetails
                Dim addnlProteins1 = objPSM.Proteins.Except(objPSM.ProteinDetails.Keys, StringComparer.CurrentCultureIgnoreCase).ToList()

                For Each proteinName In addnlProteins1
                    If mMaxProteinsPerPSM > 0 AndAlso objPSM.ProteinDetails.Count > mMaxProteinsPerPSM Then
                        ' Maximum number of proteins reached (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)
                        Exit For
                    End If

                    Dim oProtein = New clsProteinInfo(proteinName, 0, clsPeptideCleavageStateCalculator.ePeptideCleavageStateConstants.NonSpecific, clsPeptideCleavageStateCalculator.ePeptideTerminusStateConstants.None)
                    objPSM.ProteinDetails.Add(proteinName, oProtein)
                Next

                ' Make sure all of the proteins in objPSM.ProteinDetails are defined in objPSM.Proteins
                Dim addnlProteins2 = (From item In objPSM.ProteinDetails Select item.Key).Except(objPSM.Proteins, StringComparer.CurrentCultureIgnoreCase)

                Dim additionThresholdCheck As Integer = mMaxProteinsPerPSM
                If additionThresholdCheck < Integer.MaxValue Then
                    additionThresholdCheck += 1
                End If

                If mMaxProteinsPerPSM > 0 AndAlso objPSM.Proteins.Count + addnlProteins2.Count > additionThresholdCheck Then
                    ' Maximum number of proteins will be reached; only add a subset of the proteins in addnlProteins2
                    ' (note that we allow for tracking one more than the maximum because we are merging data from two different data sources)

                    For Each oProtein In addnlProteins2
                        If mMaxProteinsPerPSM > 0 AndAlso objPSM.Proteins.Count >= mMaxProteinsPerPSM Then
                            ' Maximum number of proteins reached
                            Exit For
                        End If
                        objPSM.Proteins.Add(oProtein)
                    Next
                Else
                    objPSM.Proteins.AddRange(addnlProteins2)
                End If

                If mPepToProteinMap.Count > 0 Then
                    ' Make sure the residue start/end locations are up-to-date in objPSM.ProteinDetails

                    Dim oPepToProteinMapInfo As clsPepToProteinMapInfo = Nothing
                    If mPepToProteinMap.TryGetValue(objPSM.PeptideCleanSequence, oPepToProteinMapInfo) Then

                        For Each oProtein In objPSM.ProteinDetails

                            ' Find the matching protein in oPepToProteinMapInfo
                            Dim lstLocations As List(Of clsPepToProteinMapInfo.udtProteinLocationInfo) = Nothing

                            If oPepToProteinMapInfo.ProteinMapInfo.TryGetValue(oProtein.Key, lstLocations) Then
                                Dim udtFirstLocation = lstLocations.First
                                oProtein.Value.UpdateLocationInProtein(udtFirstLocation.ResidueStart, udtFirstLocation.ResidueEnd)
                            End If
                        Next

                    End If
                End If

            End If
        End If

        If blnSuccess Then
            Dim strPrimarySequence As String = String.Empty
            Dim strPrefix As String = String.Empty
            Dim strSuffix As String = String.Empty

            If clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(objPSM.Peptide, strPrimarySequence, strPrefix, strSuffix) Then
                objPSM.PeptideWithNumericMods = strPrefix & "." & ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues) & "." & strSuffix
            Else
                objPSM.PeptideWithNumericMods = ConvertModsToNumericMods(objPSM.PeptideCleanSequence, objPSM.ModifiedResidues)
            End If

        End If

        Return blnSuccess
    End Function

    Private Function UpdatePSMFindMatchingModInfo(
      strMassCorrectionTag As String,
      blnFavorTerminalMods As Boolean,
      eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants,
      <Out()> ByRef objMatchedModDef As clsModificationDefinition) As Boolean

        objMatchedModDef = New clsModificationDefinition()

        If mModInfo Is Nothing Then
            Return False
        End If

        Dim blnMatchFound As Boolean

        Dim lstMatchedDefs = New List(Of clsModificationDefinition)

        For Each objMod In mModInfo
            If String.Equals(strMassCorrectionTag, objMod.MassCorrectionTag, StringComparison.InvariantCultureIgnoreCase) Then
                lstMatchedDefs.Add(objMod)
            End If
        Next

        blnMatchFound = False
        If lstMatchedDefs.Count > 0 Then

            Do

                If blnFavorTerminalMods Then
                    ' Look for an entry in lstMatchedDefs that is a terminal mod
                    For Each objMod In lstMatchedDefs
                        If objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse
                           objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod Then

                            If eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus AndAlso
                              (objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
                                blnMatchFound = True
                                objMatchedModDef = objMod
                                Exit For
                            End If

                            If eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus AndAlso
                              (objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS) OrElse objMod.TargetResiduesContain(clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS)) Then
                                blnMatchFound = True
                                objMatchedModDef = objMod
                                Exit For
                            End If
                        End If
                    Next
                Else
                    ' Look for an entry in lstMatchedDefs that is not a terminal mod
                    For Each objMod In lstMatchedDefs
                        If Not (objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod OrElse
                          objMod.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod) Then
                            blnMatchFound = True
                            objMatchedModDef = objMod
                            Exit For
                        End If
                    Next
                End If

                If Not blnMatchFound Then
                    If blnFavorTerminalMods Then
                        blnFavorTerminalMods = False
                    Else
                        ' Still no match found (this shouldn't happen); use the first entry in lstMatchedDefs
                        objMatchedModDef = lstMatchedDefs(0)
                        blnMatchFound = True
                    End If
                End If

            Loop While Not blnMatchFound

        End If

        Return blnMatchFound

    End Function
End Class
