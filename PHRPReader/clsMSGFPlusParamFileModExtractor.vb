Option Strict On

' This class reads a DMS-based parameter file for MSGF+ or MSPathFinder 
' to extract the dynamic and static modification information
' Example param file contents:
'
' #Modifications (see below for examples)
' StaticMod=C2H3NO,   C,   fix, any,         Carbamidomethyl           # Fixed Carbamidomethyl C
' 
' DynamicMod=O1,      M,   opt, any,         Oxidation                 # Oxidation M
' DynamicMod=HO3P,    STY, opt, any,         Phospho                   # Phosphorylation STY
' DynamicMod=H-1,     C,   opt, any,         Dehydro                   # Dehydro C
' DynamicMod=C2H2O,   *,   opt, Prot-N-term, Acetyl                    # Acetylation Protein N-term (C2H2O can be replaced with "H(2) C(2) O")
' 
'
' Note that DMS uses this information to create a Mods.txt file that is provided to MSGF+ or MSPathFinder
' When doing this, the static mods defs have ",fix," while the dynamic mod defs have ',opt'
' Example contents of an auto-generated Mods.txt file:
'
' # Static mods
' C2H3NO,C,fix,any,Carbamidomethyl     # Fixed Carbamidomethyl C
' 
' # Dynamic mods
' O1,M,opt,any,Oxidation             # Oxidation M
' HO3P,STY,opt,any,Phospho           # Phosphorylation STY
' H-1,C,opt,any,Dehydro              # Dehydro C
' C2H2O,*,opt,Prot-N-term,Acetyl     # Acetylation Protein N-term (C2H2O can be replaced with "H(2) C(2) O")
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 7/16/2015
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports System.IO
Imports System.Runtime.InteropServices

Public Class clsMSGFPlusParamFileModExtractor

#Region "Constants and Enums"

    Public Const UNKNOWN_MSGFDB_MOD_SYMBOL As Char = "?"c

    Public Const PARAM_TAG_MOD_STATIC = "StaticMod"
    Public Const PARAM_TAG_MOD_DYNAMIC = "DynamicMod"
    Public Const PARAM_TAG_CUSTOMAA = "CustomAA"

    Private Const MSGFPLUS_COMMENT_CHAR = "#"c

    Public Enum eMSGFDBModType As Integer
        Unknown = 0
        DynamicMod = 1
        StaticMod = 2
        DynNTermPeptide = 3
        DynCTermPeptide = 4
        DynNTermProtein = 5
        DynCTermProtein = 6
        CustomAA = 7
    End Enum

#End Region

#Region "Structures"

    ''' <summary>
    ''' Tracks dynamic and static modification details
    ''' Also tracks Custom amino acids
    ''' </summary>
    ''' <remarks>
    ''' Notes when tracking information for custom amino acids
    '''   ModName:    Name associated with the custom amino acid
    '''   ModMass:    Composition string, for example C5H7N1O2S0 for Hydroxyproline
    '''   ModMassVal: Computed mass of the composition string
    '''   Residues:   Single letter abbreviation for the custom amino acid, for example J or X
    '''   ModType:    eMSGFDBModType.CustomAA
    '''   ModSymbol:  ?   (a question mark; not used)
    ''' </remarks>
    Public Structure udtModInfoType
        Public ModName As String            ' Mod name (read from the parameter file) isn't used by MSGF+, but it is used by MSPathFinder
        Public ModMass As String            ' Storing as a string since reading from a text file and writing to a text file.  Also, can be a mass or an empirical formula
        Public ModMassVal As Double
        Public Residues As String
        Public ModType As eMSGFDBModType
        Public ModSymbol As Char            ' Modification symbol: *, #, @, ... ; dash if a static mod

        Public Overrides Function ToString() As String
            Return ModType.ToString() & " " & ModName & ", " & ModMass & "; " & Residues
        End Function
    End Structure
#End Region

#Region "Classwide Variables"
    Private mErrorMessage As String
    Private ReadOnly mToolName As String
#End Region

#Region "Events"
    Public Event ErrorOccurred(errMsg As String)
    Public Event WarningMessageEvent(warningMsg As String)
#End Region

#Region "Properties"

    Public ReadOnly Property ErrorMessage As String
        Get
            Return mErrorMessage
        End Get
    End Property

#End Region

    ''' <summary>
    ''' Constructor
    ''' </summary>
    ''' <param name="toolName">
    ''' Search engine name, typically MSGF+
    ''' This name is only used in log messages
    ''' </param>
    ''' <remarks></remarks>
    Public Sub New(toolName As String)
        mErrorMessage = String.Empty
        If String.IsNullOrWhiteSpace(toolName) Then
            mToolName = "Search engine"
        Else
            mToolName = toolName
        End If

    End Sub

    Private Function ComputeMass(strEmpiricalformula As String) As Double
        ' CompositionStr (C[Num]H[Num]N[Num]O[Num]S[Num]P[Num])
        ' 	- C (Carbon), H (Hydrogen), N (Nitrogen), O (Oxygen), S (Sulfer) and P (Phosphorus) are allowed.
        ' 	- Atom can be omitted.
        ' 	- Negative numbers are allowed.
        ' 	- Examples: 
        '     C2H2O1 
        '     C+2H+3N+1O+1
        '     H2C1O1
        '     H-1N-1O
        '     C3H6N2O0S1

        If String.Equals(strEmpiricalformula, "HexNAc", StringComparison.InvariantCultureIgnoreCase) Then
            ' This is a special-case modification that MSGF+ and MSPathFinder recognize
            ' It is listed in DMS as Hexosam, which means N-Acetylhexosamine
            ' It is tracked by UniMod as HexNAc
            Return 203.079376
        End If

        Dim elementalComposition As Dictionary(Of String, Integer)

        Try
            elementalComposition = clsPeptideMassCalculator.GetEmpiricalFormulaComponents(strEmpiricalformula)
        Catch ex As Exception
            ReportError(ex.Message)
            Return 0
        End Try

        Dim unknownSymbols As List(Of String) = Nothing
        Dim monoisotopicMass = clsPeptideMassCalculator.ComputeMonoistopicMass(elementalComposition, unknownSymbols)

        If Not unknownSymbols Is Nothing AndAlso unknownSymbols.Count > 0 Then
            Dim errMsg = "Error parsing empirical formula '" & strEmpiricalformula & "', "
            If unknownSymbols.Count = 1 Then
                ReportError(errMsg & "unknown element " & unknownSymbols.First)
            Else
                ReportError(errMsg & "unknown elements " & String.Join(", ", unknownSymbols))
            End If

            Return 0
        End If

        Return monoisotopicMass

    End Function

    ''' <summary>
    ''' Extracts mod info from either a MSGF+ or MSPathFinder param file or from a MSGFDB_Mods.txt file
    ''' </summary>
    ''' <param name="paramFilePath"></param>
    ''' <param name="lstModInfo"></param>
    ''' <returns>True if success; false if a problem</returns>
    ''' <remarks></remarks>
    Public Function ExtractModInfoFromParamFile(paramFilePath As String, <Out()> ByRef lstModInfo As List(Of udtModInfoType)) As Boolean

        Dim tagNamesToFind = New List(Of String)
        tagNamesToFind.Add(PARAM_TAG_MOD_STATIC)
        tagNamesToFind.Add(PARAM_TAG_MOD_DYNAMIC)
        tagNamesToFind.Add(PARAM_TAG_CUSTOMAA)

        ' Initialization
        lstModInfo = New List(Of udtModInfoType)

        Dim intUnnamedModID = 0
        mErrorMessage = String.Empty

        Try

            If String.IsNullOrEmpty(paramFilePath) Then
                ReportError(mToolName & " Parameter File name not defined; unable to extract mod info")
                Return False
            End If

            Dim fiParamFile = New FileInfo(paramFilePath)
            If Not fiParamFile.Exists Then
                ReportError(mToolName & " param file not found: " & paramFilePath)
                Return False
            End If

            ' Read the contents of the parameter (or mods) file
            Using srInFile = New StreamReader(New FileStream(fiParamFile.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite))

                Do While Not srInFile.EndOfStream
                    Dim strLineIn = srInFile.ReadLine().Trim()

                    If String.IsNullOrWhiteSpace(strLineIn) Then Continue Do

                    Dim strModSpec As String = String.Empty

                    If strLineIn.StartsWith(MSGFPLUS_COMMENT_CHAR) Then
                        ' Comment line (starts with #) 
                        ' Skip it
                        Continue Do
                    Else
                        For Each tagName In tagNamesToFind
                            strModSpec = ValidateIsValidModSpec(strLineIn, tagName)
                            If Not String.IsNullOrEmpty(strModSpec) Then
                                ' Known tag found; strLineIn was something like this:
                                '   StaticMod=C2H3N1O1,C,fix,any,Carbamidomethylation
                                '   DynamicMod=C2H3NO, *,  opt, N-term,   Carbamidomethylation
                                '   CustomAA=C5H7N1O2S0,J,custom,P,Hydroxylation     # Hydroxyproline
                                '
                                ' And strModSpec will now be something like this:
                                '   C2H3N1O1,C,fix,any,Carbamidomethylation
                                '   C2H3NO, *,  opt, N-term,   Carbamidomethylation
                                '   C5H7N1O2S0,J,custom,P,Hydroxylation     # Hydroxyproline

                                Exit For
                            End If
                        Next

                        If String.IsNullOrEmpty(strModSpec) Then
                            Dim lineInNoSpaces = TrimComment(strLineIn).Replace(" ", "")
                            If lineInNoSpaces.Contains(",opt,") OrElse lineInNoSpaces.Contains(",fix,") OrElse lineInNoSpaces.Contains(",custom,") Then
                                strModSpec = lineInNoSpaces
                            End If
                        End If

                    End If

                    If String.IsNullOrEmpty(strModSpec) Then
                        Continue Do
                    End If

                    If strModSpec.Contains("=") Then
                        ' Unrecognized tag name
                        ReportError("Mod spec '" & strModSpec & "' contains an unknown keyword before the equals sign; see parameter file " & Path.GetFileName(paramFilePath))
                        Return False
                    End If
                    ' Modification definition line found

                    ' Split the line on commas
                    Dim strSplitLine = strModSpec.Split(","c)

                    If strSplitLine.Length < 5 Then
                        Continue Do
                    End If

                    Dim udtModInfo = New udtModInfoType

                    ' Notes when tracking information for custom amino acids
                    '   ModName:    Name associated with the custom amino acid
                    '   ModMass:    Composition string, for example C5H7N1O2S0 for Hydroxyproline
                    '   ModMassVal: Computed mass of the composition string
                    '   Residues:   Single letter abbreviation for the custom amino acid, for example J or X
                    '   ModType:    eMSGFDBModType.CustomAA
                    '   ModSymbol:  ?   (a question mark; not used)

                    With udtModInfo
                        .ModMass = strSplitLine(0).Trim()

                        ' .ModMass could be a number, or could be an empirical formula
                        ' First try to parse out a number
                        If Not Double.TryParse(.ModMass, .ModMassVal) Then
                            ' Not a number
                            ' Mod (or custom AA) is specified as an empirical formula
                            ' Compute the mass
                            ' Note that ComputeMass only supports C, H, N, O, S, and P
                            .ModMassVal = ComputeMass(.ModMass)
                        End If

                        .Residues = strSplitLine(1).Trim()
                        .ModSymbol = UNKNOWN_MSGFDB_MOD_SYMBOL

                        Select Case strSplitLine(2).Trim().ToLower()
                            Case "opt"
                                .ModType = eMSGFDBModType.DynamicMod
                            Case "fix"
                                .ModType = eMSGFDBModType.StaticMod
                            Case "custom"
                                .ModType = eMSGFDBModType.CustomAA
                            Case Else
                                ReportWarning("Unrecognized Mod Type in the " & mToolName & " parameter file; should be 'opt' or 'fix'")
                                .ModType = eMSGFDBModType.DynamicMod
                        End Select

                        If .ModType <> eMSGFDBModType.CustomAA Then

                            Select Case strSplitLine(3).Trim().ToLower().Replace("-", String.Empty)
                                Case "any"
                                    ' Leave .ModType unchanged; this is a static or dynamic mod (fix or opt)
                                Case "nterm"
                                    If .ModType = eMSGFDBModType.StaticMod AndAlso .Residues <> "*" Then
                                        ' This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        .ModType = eMSGFDBModType.DynamicMod
                                    End If
                                    .Residues = clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS
                                    If .ModType = eMSGFDBModType.DynamicMod Then .ModType = eMSGFDBModType.DynNTermPeptide

                                Case "cterm"
                                    If .ModType = eMSGFDBModType.StaticMod AndAlso .Residues <> "*" Then
                                        ' This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        .ModType = eMSGFDBModType.DynamicMod
                                    End If
                                    .Residues = clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS
                                    If .ModType = eMSGFDBModType.DynamicMod Then .ModType = eMSGFDBModType.DynCTermPeptide

                                Case "protnterm"
                                    ' Includes Prot-N-Term, Prot-n-Term, ProtNTerm, etc.
                                    If .ModType = eMSGFDBModType.StaticMod AndAlso .Residues <> "*" Then
                                        ' This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        .ModType = eMSGFDBModType.DynamicMod
                                    End If
                                    .Residues = clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS
                                    If .ModType = eMSGFDBModType.DynamicMod Then .ModType = eMSGFDBModType.DynNTermProtein

                                Case "protcterm"
                                    ' Includes Prot-C-Term, Prot-c-Term, ProtCterm, etc.
                                    If .ModType = eMSGFDBModType.StaticMod AndAlso .Residues <> "*" Then
                                        ' This program does not support static mods at the N or C terminus that only apply to specific residues; switch to a dynamic mod
                                        .ModType = eMSGFDBModType.DynamicMod
                                    End If
                                    .Residues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS
                                    If .ModType = eMSGFDBModType.DynamicMod Then .ModType = eMSGFDBModType.DynCTermProtein

                                Case Else
                                    ReportWarning("Unrecognized Mod Type in the " & mToolName & " parameter file; should be 'any', 'N-term', 'C-term', 'Prot-N-term', or 'Prot-C-term'")
                            End Select
                        End If

                        .ModName = strSplitLine(4).Trim()
                        If String.IsNullOrEmpty(.ModName) Then
                            intUnnamedModID += 1
                            .ModName = "UnnamedMod" & intUnnamedModID.ToString
                        End If

                    End With

                    lstModInfo.Add(udtModInfo)

                Loop
            End Using

            Console.WriteLine()

        Catch ex As Exception
            ReportError("Error reading mod info the " & mToolName & " parameter file (" & Path.GetFileName(paramFilePath) & "): " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Private Sub ReportError(errMsg As String)
        mErrorMessage = errMsg
        RaiseEvent ErrorOccurred(errMsg)
    End Sub

    Private Sub ReportWarning(warningMsg As String)
        RaiseEvent WarningMessageEvent(warningMsg)
    End Sub

    Public Sub ResolveMSGFDBModsWithModDefinitions(lstMSGFDBModInfo As List(Of udtModInfoType), oPeptideMods As clsPeptideModificationContainer)

        Dim intResidueIndex As Integer
        Dim intResIndexStart As Integer
        Dim intResIndexEnd As Integer

        Dim chTargetResidue As Char
        Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants
        Dim eModType As clsModificationDefinition.eModificationTypeConstants
        Dim blnExistingModFound As Boolean

        Dim objModificationDefinition As clsModificationDefinition

        If Not lstMSGFDBModInfo Is Nothing Then

            ' Call .LookupModificationDefinitionByMass for each entry in lstMSGFDBModInfo
            For intIndex = 0 To lstMSGFDBModInfo.Count - 1

                Dim udtModInfo As udtModInfoType
                udtModInfo = lstMSGFDBModInfo(intIndex)

                With udtModInfo

                    If .Residues.Length > 0 Then
                        intResIndexStart = 0
                        intResIndexEnd = .Residues.Length - 1
                    Else
                        intResIndexStart = -1
                        intResIndexEnd = -1
                    End If

                    For intResidueIndex = intResIndexStart To intResIndexEnd
                        If intResidueIndex >= 0 Then
                            chTargetResidue = .Residues.Chars(intResidueIndex)
                            If chTargetResidue = "*"c Then
                                ' This is a terminal mod, and MSGFDB lists the target residue as * for terminal mods
                                ' This program requires that chTargetResidue be Nothing
                                chTargetResidue = Nothing
                            End If
                        Else
                            chTargetResidue = Nothing
                        End If

                        eModType = clsModificationDefinition.eModificationTypeConstants.DynamicMod

                        If .ModType = eMSGFDBModType.DynNTermPeptide Then
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
                        ElseIf .ModType = eMSGFDBModType.DynCTermPeptide Then
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
                        ElseIf .ModType = eMSGFDBModType.DynNTermProtein Then
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
                        ElseIf .ModType = eMSGFDBModType.DynCTermProtein Then
                            eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus
                        Else
                            Select Case chTargetResidue
                                Case clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideNTerminus
                                    If .ModType = eMSGFDBModType.StaticMod Then eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                                Case clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.PeptideCTerminus
                                    If .ModType = eMSGFDBModType.StaticMod Then eModType = clsModificationDefinition.eModificationTypeConstants.TerminalPeptideStaticMod
                                Case clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinNTerminus
                                    If .ModType = eMSGFDBModType.StaticMod Then eModType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
                                Case clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.ProteinCTerminus
                                    If .ModType = eMSGFDBModType.StaticMod Then eModType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
                                Case Else
                                    eResidueTerminusState = clsAminoAcidModInfo.eResidueTerminusStateConstants.None
                                    If .ModType = eMSGFDBModType.StaticMod Then eModType = clsModificationDefinition.eModificationTypeConstants.StaticMod
                            End Select

                        End If

                        blnExistingModFound = False

                        objModificationDefinition = oPeptideMods.LookupModificationDefinitionByMassAndModType(.ModMassVal, eModType, chTargetResidue, eResidueTerminusState, blnExistingModFound, True, clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)

                        If intResidueIndex = intResIndexStart Then
                            ' Update the Mod Symbol
                            .ModSymbol = objModificationDefinition.ModificationSymbol
                        End If

                    Next intResidueIndex

                End With

                lstMSGFDBModInfo(intIndex) = udtModInfo
            Next
        End If

    End Sub
    
    Private Shared Function TrimComment(value As String) As String

        ' Look for the MSGF+ comment character
        Dim commentCharIndex = value.IndexOf(MSGFPLUS_COMMENT_CHAR)

        If commentCharIndex > 0 Then
            ' Trim off the comment
            Return value.Substring(0, commentCharIndex).Trim()
        End If

        Return value.Trim()

    End Function

    Private Function ValidateIsValidModSpec(strLineIn As String, strModTag As String) As String

        Dim kvSetting As KeyValuePair(Of String, String)
        Dim strModSpec As String = String.Empty

        If strLineIn.ToLower.StartsWith(strModTag.ToLower()) Then

            kvSetting = clsPHRPParser.ParseKeyValueSetting(strLineIn, "="c, "#")

            If String.IsNullOrEmpty(kvSetting.Value) OrElse kvSetting.Value.ToLower() = "none" Then
                ' Not a valid mod spec
                strModSpec = String.Empty
            Else
                ' Mod spec found
                ' Note that there may be spaces before or after the commas separating the mod spec fields
                strModSpec = kvSetting.Value
            End If

        End If

        Return strModSpec

    End Function

End Class
