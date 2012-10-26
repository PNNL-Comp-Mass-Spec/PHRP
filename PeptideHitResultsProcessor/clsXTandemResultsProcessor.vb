Option Strict On

' This class reads in an X!Tandem results file (XML format) and creates 
' a tab-delimited text file with the data.  It will insert modification symbols
' into the peptide sequences for modified peptides.  The user can optionally provide
' a modification definition file which specifies the symbol to use for each
' modification mass.  If the user does not provide this file, then the modification
' definition information is determined from the X!Tandem Input Parameters section at
' the end of the X!Tandem results file.
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Program started January 2, 2006
'
' E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
' Website: http://ncrr.pnnl.gov/ or http://www.sysbio.org/resources/staff/
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

Public Class clsXTandemResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "March 25, 2008"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"
    Private Const XTANDEM_SUPPORTS_ISOTOPIC_MODS As Boolean = False

    ' Note: These names must all be lowercase
    Private Const XTANDEM_XML_ROOT_ELEMENT As String = "bioml"
    Private Const XTANDEM_XML_ELEMENT_NAME_GROUP As String = "group"
    Private Const XTANDEM_XML_ELEMENT_NAME_PROTEIN As String = "protein"
    Private Const XTANDEM_XML_ELEMENT_NAME_PEPTIDE As String = "peptide"
    Private Const XTANDEM_XML_ELEMENT_NAME_DOMAIN As String = "domain"
    Private Const XTANDEM_XML_ELEMENT_NAME_AMINOACID As String = "aa"
    Private Const XTANDEM_XML_ELEMENT_NAME_NOTE As String = "note"

    Private Const XTANDEM_XML_GROUP_TYPE_MODEL As String = "model"
    Private Const XTANDEM_XML_GROUP_TYPE_SUPPORT As String = "support"
    Private Const XTANDEM_XML_GROUP_TYPE_PARAMETERS As String = "parameters"

    Private Const XML_ERROR_ROOT_LEVEL_INVALID As String = "The data at the root level is invalid"
    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Private Const SCAN_NUMBER_EXTRACTION_REGEX_A As String = "scan=(\d+)"
    Private Const SCAN_NUMBER_EXTRACTION_REGEX_B As String = "scan\s*(\d+)"
    Private Const SCAN_NUMBER_EXTRACTION_REGEX_C As String = "(\d+)\.\d+\.\d\.dta"
    Private Const SCAN_NUMBER_EXTRACTION_REGEX_D As String = "\d+"

    Private Const REVERSED_PROTEIN_SEQUENCE_INDICATOR As String = ":reversed"
    Private Const PROTEIN_DESCRIPTION_LABEL As String = "description"

    Private Enum eCurrentXMLDataFileSectionConstants As Integer
        UnknownFile = 0
        Start = 1
        SearchResults = 2
        InputParameters = 3
        PerformanceParameters = 4
    End Enum

    Const INPUT_PARAM_LABEL_NAMES_MAX_INDEX As Integer = 13
    Private Enum eInputParamLabelNames As Integer
        Residue_StaticModMass = 0
        Residue_PotentialModMass = 1
        Residue_PotentialModMotif = 2
        Refine_PotentialModMass = 3
        Refine_PotentialModMotif = 4
        Refine_PotentialNTerminusMods = 5
        Refine_PotentialCTerminusMods = 6
        Protein_NTerminal_ResidueModMass = 7
        Protein_CTerminal_ResidueModMass = 8
        Protein_Cleavage_NTerminalMassChange = 9
        Protein_Cleavage_CTerminalMassChange = 10
        Protein_Cleavage_Site = 11
        Refine_ModificationMass = 12
        Scoring_Include_Reverse = 13
    End Enum
#End Region

#Region "Structures"
#End Region

#Region "Classwide Variables"
    Protected mNextResultID As Integer
    Protected mLookForReverseSequenceTag As Boolean

    Private mScanNumberRegExA As System.Text.RegularExpressions.Regex
    Private mScanNumberRegExB As System.Text.RegularExpressions.Regex
    Private mScanNumberRegExC As System.Text.RegularExpressions.Regex
    Private mScanNumberRegExD As System.Text.RegularExpressions.Regex
#End Region

#Region "Properties"
#End Region

    Private Function AddModificationsAndComputeMass(ByRef objSearchResult As clsSearchResultsXTandem, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
        Const ALLOW_DUPLICATE_MOD_ON_TERMINUS As Boolean = False
        Dim blnSuccess As Boolean

        Try
            ' If any modifications of type IsotopicMod are defined we would add them to the Search Result Mods now
            ' However, since X!Tandem doesn't support Isotopic Mods, this step is currently skipped
            If XTANDEM_SUPPORTS_ISOTOPIC_MODS Then
                objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)
            End If

            ' Add the protein terminus static mods (if defined and if the peptide is at a protein terminus)
            ' Function .SearchResultAddStaticTerminusMods() will only add the terminus mod if the terminus
            '  is not already modified by the given terminus mod mass
            ' This function will also add peptide terminus static mods, if defined, though those are not supported by X!Tandem and therefore should not be defined
            objSearchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts)

            ' Compute the monoisotopic mass for this peptide
            objSearchResult.ComputeMonoisotopicMass()

            ' Update PeptideDeltaMassCorrectedPpm
			objSearchResult.ComputeDelMCorrectedXT()

            ' Populate .PeptideSequenceWithMods and .PeptideModDescription
            ' Note that this function will call .AddSearchResultModificationsToCleanSequence() then .UpdateModDescription()
            objSearchResult.ApplyModificationInformation()

            blnSuccess = True
        Catch ex As Exception
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Function ConvertEValueToBase10Log(ByVal strExpectationValue As String) As String
        Dim dblEValue As Double
        Dim dblLogEValue As Double

        Try
            dblEValue = Double.Parse(strExpectationValue)
        Catch ex As Exception
            dblEValue = 0
        End Try

        Try
            If dblEValue <= 0 Then
                dblLogEValue = 0
            ElseIf dblEValue <= 1.0E-307 Then
                dblLogEValue = -307
            Else
                dblLogEValue = Math.Round(Math.Log10(dblEValue), 3)
            End If
        Catch ex As Exception
            dblLogEValue = 0
        End Try

        Return dblLogEValue.ToString("0.000")

    End Function

    Private Sub InitializeLocalVariables()
        ' Note: This function is called from ParseXTandemResultsFile()
        ' These variables will therefore be reset for each XTandem XML file analyzed
        mNextResultID = 1
        mLookForReverseSequenceTag = False

        InitializeRegExObjects()
    End Sub

    Private Sub InitializeRegExObjects()
        Const REGEX_OPTIONS As Text.RegularExpressions.RegexOptions = Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.Singleline Or Text.RegularExpressions.RegexOptions.IgnoreCase

        Try
            mScanNumberRegExA = New System.Text.RegularExpressions.Regex(SCAN_NUMBER_EXTRACTION_REGEX_A, REGEX_OPTIONS)
            mScanNumberRegExB = New System.Text.RegularExpressions.Regex(SCAN_NUMBER_EXTRACTION_REGEX_B, REGEX_OPTIONS)
            mScanNumberRegExC = New System.Text.RegularExpressions.Regex(SCAN_NUMBER_EXTRACTION_REGEX_C, REGEX_OPTIONS)
            mScanNumberRegExD = New System.Text.RegularExpressions.Regex(SCAN_NUMBER_EXTRACTION_REGEX_D, REGEX_OPTIONS)
        Catch ex As Exception
            ' Ignore errors here
        End Try

    End Sub

    Private Function ParseXTandemInputParameterModInfo(ByVal eModificationType As clsModificationDefinition.eModificationTypeConstants, ByVal intSortOrder As Integer, ByVal blnParsingMotifDef As Boolean, ByVal strParamValue As String, ByRef intModInfoCount As Integer, ByRef udtModInfo() As udtSearchOptionModificationInfoType) As Boolean
        ' Parse out the mod information defined in strParamValue
        ' Add each entry to udtModInfo
        ' If blnParsingMotifDef = True, then do not try to determine .TargetResidues

        Const MOD_LIST_SEP_CHAR As Char = ","c

        Dim intIndex As Integer
        Dim intColonIndex As Integer
        Dim intAtSignIndex As Integer

        Dim strModDefs() As String

        Dim dblModificationMass As Double
        Dim strTargetResidues As String

        Dim blnSuccess As Boolean

        Try
            If strParamValue Is Nothing OrElse strParamValue.Length = 0 Then
                ' Empty parameter; no definition to parse
            Else
                ' Parse strParamValue
                ' Each modification entry can have multiple modification definitions separated separated by commas, so we first split strParamValue
                strModDefs = strParamValue.Split(MOD_LIST_SEP_CHAR)
                For intIndex = 0 To strModDefs.Length - 1
                    ' Modification definitions typically look like "15.9949@M"
                    ' However, a neutral loss can be specified using "79.9663:-97.98@STY"
                    '   Thus, mod mass is number up to first non-numeric character (typically a colon or @ sign)
                    ' Target residues are the residues after the @
                    ' If the Target residues contain an X, then the modification can apply to any residue

                    dblModificationMass = 0

                    ' Look for a colon and an @ sign
                    intColonIndex = strModDefs(intIndex).IndexOf(":"c)
                    intAtSignIndex = strModDefs(intIndex).IndexOf("@"c)

                    If intAtSignIndex < 0 Then
                        ' At sign not found; skip this mod def
                    Else
                        If intColonIndex > intAtSignIndex Then
                            ' Ignore this colon since it's present after the @ sign
                            intColonIndex = -1
                        End If

                        If intColonIndex > 0 Then
                            ' Colon found; see if the text up to intColonIndex is a number
                            If clsPHRPBaseClass.IsNumber(strModDefs(intIndex).Substring(0, intColonIndex)) Then
                                dblModificationMass = Double.Parse(strModDefs(intIndex).Substring(0, intColonIndex))
                            End If
                        Else
                            ' Colon not found; see if the text up to intAtSignIndex is a number
                            If clsPHRPBaseClass.IsNumber(strModDefs(intIndex).Substring(0, intAtSignIndex)) Then
                                dblModificationMass = Double.Parse(strModDefs(intIndex).Substring(0, intAtSignIndex))
                            End If
                        End If

                        If dblModificationMass <> 0 Then
                            ' Valid mass found; now extract the target residues
                            If Not blnParsingMotifDef AndAlso intAtSignIndex + 1 < strModDefs(intIndex).Length Then
                                strTargetResidues = strModDefs(intIndex).Substring(intAtSignIndex + 1)
                            Else
                                strTargetResidues = String.Empty
                            End If

                            If strTargetResidues.IndexOf("X"c) >= 0 Then
                                ' Modification can affect any residue; set strTargetResidues to ""
                                strTargetResidues = String.Empty
                            End If

                            If strTargetResidues.Length > 0 Then
                                ' Convert from X!Tandem-style N-Terminus notation to DMS-style notation
								strTargetResidues = strTargetResidues.Replace(clsPeptideModificationContainer.N_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM, clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS)
								strTargetResidues = strTargetResidues.Replace(clsPeptideModificationContainer.C_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM, clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS)
                            End If

                            ' Append the new mod information to udtModInfo
                            If intModInfoCount >= udtModInfo.Length Then
                                ReDim Preserve udtModInfo(udtModInfo.Length * 2 - 1)
                            End If

                            With udtModInfo(intModInfoCount)
                                .SortOrder = intSortOrder
                                .ModificationMass = dblModificationMass
                                .TargetResidues = strTargetResidues
                                .ModificationType = eModificationType
                            End With
                            intModInfoCount += 1
                        End If
                    End If
                Next intIndex
            End If

            blnSuccess = True

        Catch ex As Exception
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Function ParseXTandemInputParameterProteinTerminusMod(ByVal intSortOrder As Integer, ByVal blnNTerminus As Boolean, ByVal strParamValue As String, ByRef intModInfoCount As Integer, ByRef udtModInfo() As udtSearchOptionModificationInfoType) As Boolean
        ' Parse out the mass defined in strParamValue
        ' Add the entry to udtModInfo if non-zero
        ' if blnNTerminus = True then mod applies to the protein's N-Terminus; otherwise, applies to the protein's C-Terminus

        Dim dblModificationMass As Double
        Dim strTargetResidues As String

        Dim blnSuccess As Boolean

        Try
            If strParamValue Is Nothing OrElse strParamValue.Length = 0 Then
                ' Empty parameter; no definition to parse
            Else
                ' See if strParamValue is a non-zero number

                dblModificationMass = 0
                If clsPHRPBaseClass.IsNumber(strParamValue) Then
                    dblModificationMass = Double.Parse(strParamValue)
                End If

                If dblModificationMass <> 0 Then
                    ' Append the new mod information to udtModInfo
                    If intModInfoCount >= udtModInfo.Length Then
                        ReDim Preserve udtModInfo(udtModInfo.Length * 2 - 1)
                    End If

                    With udtModInfo(intModInfoCount)
                        .SortOrder = intSortOrder
                        If blnNTerminus Then
							strTargetResidues = clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS
                        Else
							strTargetResidues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS
                        End If

                        .ModificationMass = dblModificationMass
                        .TargetResidues = strTargetResidues
                        .ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod
                    End With
                    intModInfoCount += 1
                End If
            End If

            blnSuccess = True

        Catch ex As Exception
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Function ParseXTandemResultsFile(ByVal strInputFilePath As String, ByVal strOutputFilePath As String, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile

		Dim strModificationSummaryFilePath As String

        ' Note: The number of valid entries in objSearchResults() is given by intSearchResultCount; objSearchResults() is expanded but never shrunk
        ' There is a separate entry in objSearchResults() for each protein encountered
        Dim intSearchResultCount As Integer
        Dim objSearchResults() As clsSearchResultsXTandem

        Dim eCurrentXMLDataFileSection As eCurrentXMLDataFileSectionConstants
        Dim strCurrentGroupType As String = String.Empty

        Dim intGroupElementReaderDepth As Integer

        Dim intSearchResultIndex As Integer
        Dim intResultsProcessed As Integer
		Dim sngPercentComplete As Single

        Dim blnSuccess As Boolean

        Dim strErrorLog As String = String.Empty

        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            ' Initialize the objSearchResults() array; initially reserve space for 4 proteins
            intSearchResultCount = 0
            ReDim objSearchResults(3)
            For intSearchResultIndex = 0 To objSearchResults.Length - 1
                objSearchResults(intSearchResultIndex) = New clsSearchResultsXTandem(mPeptideMods)
            Next intSearchResultIndex

            ' Reset mNextResultID and mLookForReverseSequenceTag
            InitializeLocalVariables()

            Try
                ' Read the input parameters from the end of the X!Tandem results file (strInputFilePath)
                blnSuccess = ParseXTandemResultsFileInputParameters(strInputFilePath)
                If Not blnSuccess Then
                    SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile, True)
                    Return False
                End If

                For intSearchResultIndex = 0 To objSearchResults.Length - 1
                    UpdateSearchResultEnzymeAndTerminusInfo(objSearchResults(intSearchResultIndex))
                Next intSearchResultIndex

                ' Open the input file and parse it

                ' Initialize the stream reader and the XML Text Reader
				Using srDataFile As System.IO.StreamReader = New System.IO.StreamReader(strInputFilePath)
					Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(srDataFile)
						strErrorLog = String.Empty
						intResultsProcessed = 0

						' Create the output file
						Using swPeptideResultsFile As System.IO.StreamWriter = New System.IO.StreamWriter(strOutputFilePath, False)

							' Write the header line to swPeptideResultsFile
							swPeptideResultsFile.WriteLine( _
							  "Result_ID" & SEP_CHAR & _
							  "Group_ID" & SEP_CHAR & _
							  "Scan" & SEP_CHAR & _
							  "Charge" & SEP_CHAR & _
							  "Peptide_MH" & SEP_CHAR & _
							  "Peptide_Hyperscore" & SEP_CHAR & _
							  "Peptide_Expectation_Value_Log(e)" & SEP_CHAR & _
							  "Multiple_Protein_Count" & SEP_CHAR & _
							  "Peptide_Sequence" & SEP_CHAR & _
							  "DeltaCn2" & SEP_CHAR & _
							  "y_score" & SEP_CHAR & _
							  "y_ions" & SEP_CHAR & _
							  "b_score" & SEP_CHAR & _
							  "b_ions" & SEP_CHAR & _
							  "Delta_Mass" & SEP_CHAR & _
							  "Peptide_Intensity_Log(I)" & SEP_CHAR & _
							  "DelM_PPM")

							' Create the additional output files
							blnSuccess = MyBase.InitializeSequenceOutputFiles(strOutputFilePath)

							' Parse the input file
							eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile

							Do While objXMLReader.Read() And Not MyBase.AbortProcessing

								XMLTextReaderSkipWhitespace(objXMLReader)
								If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

								If objXMLReader.Depth < 2 Then
									If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
										Select Case objXMLReader.Name.ToLower
											Case XTANDEM_XML_ELEMENT_NAME_GROUP
												If objXMLReader.HasAttributes Then
													' Cache the XML reader depth before reading any of the element's attributes
													intGroupElementReaderDepth = objXMLReader.Depth

													' See if the group has a "type" attribute containing the text XTANDEM_XML_GROUP_TYPE_MODEL
													strCurrentGroupType = XMLTextReaderGetAttributeValue(objXMLReader, "type", String.Empty)
													If strCurrentGroupType = XTANDEM_XML_GROUP_TYPE_MODEL Then
														eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.SearchResults

														blnSuccess = ParseXTandemResultsFileEntry(objXMLReader, swPeptideResultsFile, intSearchResultCount, objSearchResults, strErrorLog, intGroupElementReaderDepth)
														intResultsProcessed += 1

														' Update the progress
														sngPercentComplete = CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100)
														If mCreateProteinModsFile Then
															sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100)
														End If
														UpdateProgress(sngPercentComplete)

													End If
												Else
													' Group doesn't have any attributes; ignore it
													objXMLReader.Skip()
												End If

											Case XTANDEM_XML_ROOT_ELEMENT
												eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.Start
										End Select
									End If
								End If
							Loop

						End Using
					End Using
				End Using

				If eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile Then
					mErrorMessage = "Root element '" & XTANDEM_XML_ROOT_ELEMENT & "' not found in the input file: " & ControlChars.NewLine & strInputFilePath
				Else
					If mCreateModificationSummaryFile Then
						' Create the modification summary file
						Dim fiInputFile As System.IO.FileInfo = New System.IO.FileInfo(strInputFilePath)
						Dim fiOutputFile As System.IO.FileInfo = New System.IO.FileInfo(strOutputFilePath)

						strModificationSummaryFilePath = System.IO.Path.GetFileName(MyBase.ReplaceFilenameSuffix(fiInputFile, FILENAME_SUFFIX_MOD_SUMMARY))
						strModificationSummaryFilePath = System.IO.Path.Combine(fiOutputFile.DirectoryName, strModificationSummaryFilePath)

						SaveModificationSummaryFile(strModificationSummaryFilePath)
					End If

					' Inform the user if any errors occurred
					If strErrorLog.Length > 0 Then
						SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
					End If

				End If

				blnSuccess = True

			Catch ex As Exception
				SetErrorMessage(ex.Message)
				SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
				blnSuccess = False
			Finally
				MyBase.CloseSequenceOutputFiles()
			End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Private Function ParseXTandemResultsFileEntry(ByRef objXMLReader As System.Xml.XmlTextReader, ByRef swPeptideResultsFile As System.IO.StreamWriter, ByRef intSearchResultCount As Integer, ByRef objSearchResults() As clsSearchResultsXTandem, ByRef strErrorLog As String, ByVal intGroupElementReaderDepth As Integer) As Boolean
        ' Note: The number of valid entries in objSearchResults() is given by intSearchResultCount; objSearchResults() is expanded but never shrunk
        ' There is a separate entry in objSearchResults() for each protein encountered

        Const GROUP_LABEL_PROTEIN As String = "protein"
        Const GROUP_LABEL_FRAG_ION As String = "fragment ion mass spectrum"

        ' The following are static to avoid re-reserving space for them for every XTandem results file entry
        Static htSeqsWithMods As Hashtable
        Static htSeqsWithoutMods As Hashtable

        Dim strSequenceWithMods As String

        ' The following is the XSLT originally used to parse out the data
        ' ID: <xsl:value-of select="@id" />
        ' Charge: <xsl:value-of select="@z" />
        ' ParentIonMH: <xsl:value-of select="@mh" />
        ' Peptide_Expectation_Value_e: <xsl:value-of select="@expect" />
        ' Protein_Name: <xsl:value-of select="substring-before(concat(@label,' '), ' ')" />
        ' Peptide_Intensity_Log(I): <xsl:value-of select="@sumI" />

        ' Protein_Expectation_Value_Log(e): <xsl:value-of select="./protein/@expect" />
        ' Protein_Intensity_Log(I): <xsl:value-of select="./protein/@sumI" />

        ' Peptide_Hyperscore: <xsl:value-of select="./protein/peptide/domain/@hyperscore" />
        ' Peptide_Sequence: <xsl:if test = "substring(./protein/peptide/domain/@pre,string-length(./protein/peptide/domain/@pre),1) = '['" ><xsl:text>-</xsl:text></xsl:if><xsl:if test = "substring(./protein/peptide/domain/@pre,string-length(./protein/peptide/domain/@pre),1) != '['" ><xsl:value-of select="substring(./protein/peptide/domain/@pre,string-length(./protein/peptide/domain/@pre),1)" /></xsl:if>.<xsl:value-of select="./protein/peptide/domain/@seq" />.<xsl:if test = "substring(./protein/peptide/domain/@post,1,1) = ']'" ><xsl:text>-</xsl:text></xsl:if><xsl:if test = "substring(./protein/peptide/domain/@post,1,1) != ']'" ><xsl:value-of select="substring(./protein/peptide/domain/@post,1,1)" /></xsl:if>
        ' DeltaCn2: <xsl:value-of select="round((./protein/peptide/domain/@hyperscore - ./protein/peptide/domain/@nextscore) div ./protein/peptide/domain/@hyperscore * 10000) div 10000" />
        ' y_score: <xsl:value-of select="./protein/peptide/domain/@y_score" />
        ' y_ions: <xsl:value-of select="./protein/peptide/domain/@y_ions" />
        ' b_score: <xsl:value-of select="./protein/peptide/domain/@b_score" />
        ' b_ions: <xsl:value-of select="./protein/peptide/domain/@b_ions" />
        ' Delta_Mass: <xsl:value-of select="./protein/peptide/domain/@delta" />

        ' Scan: <xsl:value-of select="substring-before(concat(substring-after(./group/note,'scan='),' '), ' ')" />

		Dim strGroupIDInXMLFile As String = String.Empty
        Dim strCurrentGroupType As String = String.Empty
        Dim strCurrentGroupLabel As String = String.Empty

        Dim intIndex As Integer
        Dim intSearchResultIndex As Integer
        Dim intDomainElementReaderDepth As Integer

        Dim intProteinCount As Integer

        Dim blnSuccess As Boolean
        Dim blnUpdateModOccurrenceCounts As Boolean
        Dim blnUpdateResultToSeqMapFile As Boolean
        Dim blnScanFound As Boolean

        Dim blnProteinSequenceParsed As Boolean
        Dim blnDomainParsed As Boolean

        Dim strValue As String

        blnProteinSequenceParsed = False
        blnDomainParsed = False
        blnSuccess = False

        strCurrentGroupType = XTANDEM_XML_GROUP_TYPE_MODEL
        strCurrentGroupLabel = GROUP_LABEL_PROTEIN

        Try
            ' Reset the first entry in objSearchResults
            intSearchResultCount = 0
            objSearchResults(0).Clear()

            strGroupIDInXMLFile = XMLTextReaderGetAttributeValue(objXMLReader, "id", String.Empty)
            With objSearchResults(0)
                ' Initially set .ResultID to strGroupIDInXMLFile
                ' ResultID will get updated to a sequentially assigned number (mNextResultID) if we write the result out to the _xt.txt file
                .ResultID = Integer.Parse(strGroupIDInXMLFile)
                .GroupID = .ResultID

                .ParentIonMH = XMLTextReaderGetAttributeValue(objXMLReader, "mh", String.Empty)
                .Charge = XMLTextReaderGetAttributeValue(objXMLReader, "z", String.Empty)

                ' Note: This will get updated for each protein encountered
                .PeptideExpectationValue = ConvertEValueToBase10Log(XMLTextReaderGetAttributeValue(objXMLReader, "expect", String.Empty))

                ' Note: we truncate .ProteinName at the first space
                .ProteinName = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)
                intIndex = .ProteinName.IndexOf(" "c)
                If intIndex > 0 Then
                    .ProteinName = .ProteinName.Substring(0, intIndex)
                End If

                .PeptideIntensity = XMLTextReaderGetAttributeValue(objXMLReader, "sumI", String.Empty)
                .PeptideIntensityMax = XMLTextReaderGetAttributeValue(objXMLReader, "maxI", String.Empty)
                .fI = XMLTextReaderGetAttributeValue(objXMLReader, "fI", String.Empty)
            End With

            ' Continue reading the XML file, loading the information 
            Do While objXMLReader.Read() And Not MyBase.AbortProcessing

                XMLTextReaderSkipWhitespace(objXMLReader)
                If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

                If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
                    Select Case objXMLReader.Name.ToLower
                        Case XTANDEM_XML_ELEMENT_NAME_PROTEIN
                            intSearchResultCount += 1
                            If intSearchResultCount > objSearchResults.Length Then
                                ' Double the length of objSearchResults
                                ReDim Preserve objSearchResults(objSearchResults.Length * 2 - 1)
                                For intSearchResultIndex = intSearchResultCount - 1 To objSearchResults.Length - 1
                                    objSearchResults(intSearchResultIndex) = New clsSearchResultsXTandem(mPeptideMods)
                                    UpdateSearchResultEnzymeAndTerminusInfo(objSearchResults(intSearchResultIndex))
                                Next intSearchResultIndex
                            End If

                            With objSearchResults(intSearchResultCount - 1)
                                If intSearchResultCount > 1 Then
                                    ' Copy the info from objSearchResults(0) to this search result
                                    ' ResultID will get updated to a sequentially assigned number (mNextResultID) if we write the result out to the _xt.txt file
                                    .ResultID = objSearchResults(0).ResultID
                                    .GroupID = objSearchResults(0).GroupID
                                    .ParentIonMH = objSearchResults(0).ParentIonMH
                                    .Charge = objSearchResults(0).Charge

                                    .PeptideExpectationValue = objSearchResults(0).PeptideExpectationValue

                                    .PeptideIntensity = objSearchResults(0).PeptideIntensity
                                    .PeptideIntensityMax = objSearchResults(0).PeptideIntensityMax
                                    .fI = objSearchResults(0).fI
                                End If

                                .ProteinExpectationValue = XMLTextReaderGetAttributeValue(objXMLReader, "expect", String.Empty)
                                .ProteinIntensity = XMLTextReaderGetAttributeValue(objXMLReader, "sumI", String.Empty)

                                ' Update the protein name for this protein entry

                                ' Note: we truncate .ProteinName at the first space
                                ' However, if mLookForReverseSequenceTag = True then we need to look for ":reversed" at the end of the description
                                .ProteinName = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)
                                .ProteinName = TruncateProteinName(.ProteinName)

                                ' For proteins with long descriptions, the ":reversed" tag is not present in the label attribute
                                '  and is instead in the <note label="description"> element (a sub-element of the <protein> element
                                ' We'll check for this case later in this function

                                ' Reset the Protein Sequence Parsed and Domain Parsed flags
                                blnProteinSequenceParsed = False
                                blnDomainParsed = False

                                ' Clear the protein sequence info and peptide details info
                                .ClearProteinSequenceInfo()
                                .ClearPeptideDetailsInfo()
                            End With

                        Case XTANDEM_XML_ELEMENT_NAME_PEPTIDE
                            If Not blnProteinSequenceParsed Then
                                blnProteinSequenceParsed = True
                                With objSearchResults(intSearchResultCount - 1)
                                    .ProteinSeqResidueNumberStart = XMLTextReaderGetAttributeValue(objXMLReader, "start", 0)
                                    .ProteinSeqResidueNumberEnd = XMLTextReaderGetAttributeValue(objXMLReader, "end", 0)
                                End With
                            End If
                        Case XTANDEM_XML_ELEMENT_NAME_DOMAIN
                            ' If the given peptide is present in the protein more than once, then domain will appear twice
                            '   Only keep the data for the first domain
                            ' Additionally, if X!Tandem decides a given modification could occur on either of two residues, it repeats the domain information
                            '   Again, in this situation, we'll only keep the first domain

                            If Not blnDomainParsed Then
                                blnDomainParsed = True

                                ' Cache the XML reader depth before reading any of the element's attributes
                                intDomainElementReaderDepth = objXMLReader.Depth

                                ' Read the information for this domain, storing in objSearchResult
                                With objSearchResults(intSearchResultCount - 1)
                                    .PeptideLocInProteinStart = XMLTextReaderGetAttributeValue(objXMLReader, "start", 0)
                                    .PeptideLocInProteinEnd = XMLTextReaderGetAttributeValue(objXMLReader, "end", 0)

                                    ' Note: the expectation value was already populated from the group level; we'll update it to the value stored for each protein in case it's different for different proteins
                                    .PeptideExpectationValue = ConvertEValueToBase10Log(XMLTextReaderGetAttributeValue(objXMLReader, "expect", String.Empty))

                                    ' Note: This is the theoretical, monoisotopic MH+ mass for the peptide
                                    .PeptideMH = XMLTextReaderGetAttributeValue(objXMLReader, "mh", String.Empty)

                                    .PeptideDeltaMass = XMLTextReaderGetAttributeValue(objXMLReader, "delta", String.Empty)

                                    ' Note: .peptideDeltaMass is stored in the X!Tandem XML file as "Observed_Mass - Theoretical_Mass"
                                    ' However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                                    ' Therefore, we will negate .peptideDeltaMass
                                    Try
                                        .PeptideDeltaMass = (-Double.Parse(.PeptideDeltaMass)).ToString
                                    Catch ex As Exception
                                        ' Error; Leave .peptideDeltaMass unchanged
                                    End Try

                                    .PeptideHyperscore = XMLTextReaderGetAttributeValue(objXMLReader, "hyperscore", String.Empty)

                                    ' Note that updating .PeptideNextScore will automatically populate .DeltaCn2
                                    .PeptideNextScore = XMLTextReaderGetAttributeValue(objXMLReader, "nextscore", String.Empty)

                                    ' Note that calling .PeptidePreResidues, .PeptidePostResidues, and .PeptideCleanSequence will call ComputePeptideCleavageStateInProtein() each time
                                    .PeptidePreResidues = XMLTextReaderGetAttributeValue(objXMLReader, "pre", String.Empty)
                                    .PeptidePostResidues = XMLTextReaderGetAttributeValue(objXMLReader, "post", String.Empty)
                                    .PeptideCleanSequence = XMLTextReaderGetAttributeValue(objXMLReader, "seq", String.Empty)

                                    .PeptideYScore = XMLTextReaderGetAttributeValue(objXMLReader, "y_score", String.Empty)
                                    .PeptideYIons = XMLTextReaderGetAttributeValue(objXMLReader, "y_ions", String.Empty)
                                    .PeptideBScore = XMLTextReaderGetAttributeValue(objXMLReader, "b_score", String.Empty)
                                    .PeptideBIons = XMLTextReaderGetAttributeValue(objXMLReader, "b_ions", String.Empty)
                                End With

                                ' Now read all of the mods for this domain
                                ' If this is the first search result, then update the mod occurrence counts; otherwise, do not
                                If intSearchResultCount = 1 Then
                                    ParseXTandemResultsFileReadDomainMods(objXMLReader, objSearchResults(intSearchResultCount - 1), intDomainElementReaderDepth, True)
                                Else
                                    ParseXTandemResultsFileReadDomainMods(objXMLReader, objSearchResults(intSearchResultCount - 1), intDomainElementReaderDepth, False)
                                End If
                            End If

                        Case XTANDEM_XML_ELEMENT_NAME_GROUP
                            strValue = XMLTextReaderGetAttributeValue(objXMLReader, "type", String.Empty)
                            If strValue.Length > 0 Then
                                strCurrentGroupType = String.Copy(strValue)
                                strCurrentGroupLabel = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)
                            Else
                                ' Leave strGroupType unchanged
                            End If

                        Case XTANDEM_XML_ELEMENT_NAME_NOTE
                            If strCurrentGroupType = XTANDEM_XML_GROUP_TYPE_MODEL AndAlso _
                               strCurrentGroupLabel = GROUP_LABEL_PROTEIN AndAlso mLookForReverseSequenceTag Then
                                ' Examine the label attribute of this note

                                strValue = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)

                                If strValue = PROTEIN_DESCRIPTION_LABEL Then
                                    ' Check whether this note ends in ":reversed"
                                    ' If it does, then make sure objSearchResults(intSearchResultCount - 1).ProteinName ends in :reversed

                                    ' Advance the reader before grabbing the inner text
                                    If objXMLReader.Read() Then
                                        ' Read the note's inner text
                                        strValue = XMLTextReaderGetInnerText(objXMLReader)

                                        If strValue.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR) Then
                                            If Not objSearchResults(intSearchResultCount - 1).ProteinName.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR) Then
                                                objSearchResults(intSearchResultCount - 1).ProteinName &= REVERSED_PROTEIN_SEQUENCE_INDICATOR
                                            End If
                                        End If
                                    End If
                                End If

                            ElseIf strCurrentGroupType = XTANDEM_XML_GROUP_TYPE_SUPPORT AndAlso _
                               strCurrentGroupLabel = GROUP_LABEL_FRAG_ION Then
                                ' This note should contain the scan number
                                ' For _Dta.txt files created at PNNL, it should look something like: "   scan=15118 cs=3"
                                ' For DTA-based files converted to .MGF and processed by X!Tandem: "MyDataset.300.300.2.dta"

                                ' Read the note's inner text
                                strValue = XMLTextReaderGetInnerText(objXMLReader)
                                If Not strValue Is Nothing Then
                                    blnScanFound = False

                                    ' Look for the word "scan" followed by an equals sign, followed by a number
                                    ' For example, "   scan=15118 cs=3"
                                    Try
                                        With mScanNumberRegExA.Match(strValue)
                                            If .Success AndAlso .Groups.Count > 1 Then
                                                objSearchResults(0).Scan = .Groups(1).Value
                                                blnScanFound = True
                                            End If
                                        End With
                                    Catch ex As Exception
                                        ' Ignore errors here
                                    End Try

                                    If Not blnScanFound Then
                                        ' No match; look for the word "scan" followed by whitespace, followed by a number
                                        ' For example, "scan 300"
                                        Try
                                            With mScanNumberRegExB.Match(strValue)
                                                If .Success AndAlso .Groups.Count > 1 Then
                                                    objSearchResults(0).Scan = .Groups(1).Value
                                                    blnScanFound = True
                                                End If
                                            End With
                                        Catch ex As Exception
                                            ' Ignore errors here
                                        End Try
                                    End If

                                    If Not blnScanFound Then
                                        ' No match; see if the description resembles a .Dta file name
                                        ' For example, "MyDataset.300.300.2.dta"
                                        Try
                                            With mScanNumberRegExC.Match(strValue)
                                                If .Success AndAlso .Groups.Count > 1 Then
                                                    objSearchResults(0).Scan = .Groups(1).Value
                                                    blnScanFound = True
                                                End If
                                            End With
                                        Catch ex As Exception
                                            ' Ignore errors here
                                        End Try
                                    End If

                                    If Not blnScanFound Then
                                        ' Still no match; extract out the first number present
                                        Try
                                            With mScanNumberRegExD.Match(strValue)
                                                If .Success Then
                                                    objSearchResults(0).Scan = .Value
                                                Else
                                                    objSearchResults(0).Scan = strValue
                                                End If
                                            End With
                                        Catch ex As Exception
                                            ' Ignore errors here
                                        End Try
                                    End If

                                    ' Copy the scan value from the first result to the other results
                                    For intSearchResultIndex = 1 To intSearchResultCount - 1
                                        objSearchResults(intSearchResultIndex).Scan = objSearchResults(0).Scan
                                    Next intSearchResultIndex

                                End If
                            End If
                        Case Else
                            ' Unknown/unneeded child node name; ignore it
                    End Select

                ElseIf objXMLReader.NodeType = Xml.XmlNodeType.EndElement Then
                    If objXMLReader.Name = XTANDEM_XML_ELEMENT_NAME_GROUP Then
                        If objXMLReader.Depth <= intGroupElementReaderDepth Then
                            ' End element found for the current group

                            ' Typically each group will consist of entries all having the same sequence and modifications (but different protein names)
                            ' However, occasionally a group will contain a mix of peptide sequences (only occurs if they all had the exact same hyperscore)
                            ' In order to check for this, we will construct a pointer array of Sequence and Mods to SearchResultIndex and use this to determine
                            '  which entries should be written to the _xt.txt file and to the ResultToSeqMap file
                            ' We will also use this pointer array to keep track of the number of proteins listed for each peptide

                            If htSeqsWithMods Is Nothing Then
                                htSeqsWithMods = New Hashtable
                                htSeqsWithoutMods = New Hashtable
                            Else
                                htSeqsWithMods.Clear()
                                htSeqsWithoutMods.Clear()
                            End If

                            ' First step through the results to compute the mass, construct the modification description, 
                            '  and determine the number of proteins listed for each
                            For intSearchResultIndex = 0 To intSearchResultCount - 1
                                If intSearchResultIndex = 0 Then
                                    ' Always set blnUpdateModOccurrenceCounts to True for the first result in the group
                                    blnUpdateModOccurrenceCounts = True
                                    htSeqsWithoutMods.Add(objSearchResults(intSearchResultIndex).PeptideCleanSequence, 1)
                                Else
                                    If htSeqsWithoutMods.ContainsKey(objSearchResults(intSearchResultIndex).PeptideCleanSequence) Then
                                        blnUpdateModOccurrenceCounts = False
                                    Else
                                        blnUpdateModOccurrenceCounts = True
                                        htSeqsWithoutMods.Add(objSearchResults(intSearchResultIndex).PeptideCleanSequence, 1)
                                    End If
                                End If

                                blnSuccess = AddModificationsAndComputeMass(objSearchResults(intSearchResultIndex), blnUpdateModOccurrenceCounts)
                                If Not blnSuccess Then
                                    If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                                        strErrorLog &= "Error adding modifications to sequence for Group ID '" & strGroupIDInXMLFile & "'" & ControlChars.NewLine
                                    End If
                                End If

                                strSequenceWithMods = objSearchResults(intSearchResultIndex).PeptideCleanSequence & "_" & objSearchResults(intSearchResultIndex).PeptideModDescription

                                If intSearchResultIndex = 0 Then
                                    ' Always add the first result in the group htSeqsWithMods
                                    htSeqsWithMods.Add(strSequenceWithMods, 1)
                                Else
                                    ' See if htSeqsWithMods contains strSequenceWithMods
                                    If htSeqsWithMods.ContainsKey(strSequenceWithMods) Then
                                        ' Increment the protein count for this peptide
                                        htSeqsWithMods(strSequenceWithMods) = DirectCast(htSeqsWithMods(strSequenceWithMods), Integer) + 1
                                    Else
                                        htSeqsWithMods.Add(strSequenceWithMods, 1)
                                    End If
                                End If
                            Next intSearchResultIndex


                            ' Now step through the list again and update the MultipleProteinCount value for each search result
                            For intSearchResultIndex = 0 To intSearchResultCount - 1
                                strSequenceWithMods = objSearchResults(intSearchResultIndex).PeptideCleanSequence & "_" & objSearchResults(intSearchResultIndex).PeptideModDescription

                                Try
                                    intProteinCount = DirectCast(htSeqsWithMods(strSequenceWithMods), Integer)
                                Catch ex As Exception
                                    intProteinCount = 1
                                End Try

                                If intProteinCount < 1 Then intProteinCount = 1

                                ' Note: Multiple protein count is 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
                                objSearchResults(intSearchResultIndex).MultipleProteinCount = (intProteinCount - 1).ToString
                            Next intSearchResultIndex

                            ' Clear htSeqsWithMods again since we need to re-use it to determine which results to write out
                            htSeqsWithMods.Clear()

                            ' Write out the results
                            For intSearchResultIndex = 0 To intSearchResultCount - 1
                                strSequenceWithMods = objSearchResults(intSearchResultIndex).PeptideCleanSequence & "_" & objSearchResults(intSearchResultIndex).PeptideModDescription

                                If intSearchResultIndex = 0 Then
                                    ' Always save the first result in the group to the _xt.txt and _ResultToSeqMap.txt files
                                    htSeqsWithMods.Add(strSequenceWithMods, 1)
                                    blnUpdateResultToSeqMapFile = True
                                Else
                                    ' See if htSeqsWithMods contains strSequenceWithMods
                                    If htSeqsWithMods.ContainsKey(strSequenceWithMods) Then
                                        blnUpdateResultToSeqMapFile = False
                                    Else
                                        htSeqsWithMods.Add(strSequenceWithMods, 1)
                                        blnUpdateResultToSeqMapFile = True
                                    End If
                                End If

                                If blnUpdateResultToSeqMapFile Then
                                    ' Only save the first result for each peptide in the group to the _xt.txt and _ResultToSeqMap.txt files
                                    ' Note: This function will update .ResultID to the next available ID value (mNextResultID)
                                    SaveXTandemResultsFileEntry(objSearchResults(intSearchResultIndex), swPeptideResultsFile)
                                End If

                                MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResults(intSearchResultIndex), clsSearchResultsBaseClass), blnUpdateResultToSeqMapFile)

                            Next intSearchResultIndex

                            blnSuccess = True
                            Exit Do
                        End If
                    End If
                End If
            Loop

        Catch ex As Exception
            ' Error parsing values from this group ID in the XML file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error parsing value for Group ID '" & strGroupIDInXMLFile & "'" & ControlChars.NewLine
            End If
            blnSuccess = False
        End Try

        Return blnSuccess
    End Function

    Private Function ParseXTandemResultsFileInputParameters(ByVal strInputFilePath As String) As Boolean

        ' Pre-read the XML file and look for the Input Parameters section
        ' Read the parameters and validate that each of the mods defined is present in mPeptideMods

        Const GROUP_LABEL_INPUT_PARAMETERS As String = "input parameters"

        Dim eCurrentXMLDataFileSection As eCurrentXMLDataFileSectionConstants

        Dim strCurrentGroupType As String
        Dim strCurrentGroupLabel As String

        Dim intParametersGroupDepth As Integer

        Dim blnSuccess As Boolean

        Try
            ' Open the input file and parse it
            ' Initialize the stream reader and the XML Text Reader
			Using objXMLReader As System.Xml.XmlTextReader = New System.Xml.XmlTextReader(strInputFilePath)

				' Parse the file
				eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile

				Do While objXMLReader.Read() And Not MyBase.AbortProcessing

					XMLTextReaderSkipWhitespace(objXMLReader)
					If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

					If objXMLReader.Depth < 2 Then
						If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
							Select Case objXMLReader.Name.ToLower
								Case XTANDEM_XML_ELEMENT_NAME_GROUP
									If objXMLReader.HasAttributes Then
										intParametersGroupDepth = objXMLReader.Depth

										' See if the group has a "type" attribute containing the text XTANDEM_XML_GROUP_TYPE_PARAMETERS
										strCurrentGroupType = XMLTextReaderGetAttributeValue(objXMLReader, "type", String.Empty)
										If strCurrentGroupType = XTANDEM_XML_GROUP_TYPE_PARAMETERS Then
											eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.InputParameters

											' Read the Label for this group
											strCurrentGroupLabel = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)
											If strCurrentGroupLabel = GROUP_LABEL_INPUT_PARAMETERS Then
												' Read the input parameters
												ParseXTandemResultsFileInputParametersWork(objXMLReader, intParametersGroupDepth)
											End If
										Else
											' Skip this group
											objXMLReader.Skip()
										End If
									Else
										' Group doesn't have any attributes; ignore it
										objXMLReader.Skip()
									End If

								Case XTANDEM_XML_ROOT_ELEMENT
									eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.Start
								Case Else
									' Skip this element
									objXMLReader.Skip()
							End Select
						End If
					End If
				Loop

			End Using

			If eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile Then
				SetErrorMessage("Root element '" & XTANDEM_XML_ROOT_ELEMENT & "' not found in the input file: " & strInputFilePath)
				blnSuccess = False
			Else
				blnSuccess = True
			End If

		Catch ex As Exception
			SetErrorMessage(ex.Message)
			SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
			blnSuccess = False
		End Try

        Return blnSuccess
    End Function

    Private Sub ParseXTandemResultsFileInputParametersWork(ByRef objXMLReader As System.Xml.XmlTextReader, ByVal intParametersGroupDepth As Integer)

        ' Read the input parameters
        ' Each parameter is an element with name "note" with attributes "type" and "label"
        ' The parameter value is the text between within the element

        Const NOTE_TYPE_INPUT As String = "input"
        Const CLEAVAGE_SPEC_SEP As Char = "|"c
        Const XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START As Char = "{"c
        Const XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END As Char = "}"c

        Dim intIndex As Integer
        Dim intIndexCompare As Integer
        Dim intIndexCopy As Integer
        Dim intBarLoc As Integer

        Dim strNoteType As String
        Dim strNoteLabel As String
        Dim strValue As String

        Dim strLeftSpec As String
        Dim strRightSpec As String

        ' Note: When populating this we use .ToLower to make sure all of the text is lowercase
        Dim udtParamLabels() As String

        Dim intModInfoCount As Integer
        Dim udtModInfo() As udtSearchOptionModificationInfoType

        Dim blnSuccess As Boolean
        Dim blnModDeleted As Boolean
        Dim blnStaticModsAreResetForRefinement As Boolean

        ' Initialize the Mod Info array
        intModInfoCount = 0
        ReDim udtModInfo(19)

        blnStaticModsAreResetForRefinement = True

        ' Initialize udtParamLabels; this specifies the parameters to examine
        ReDim udtParamLabels(INPUT_PARAM_LABEL_NAMES_MAX_INDEX)
        udtParamLabels(eInputParamLabelNames.Residue_StaticModMass) = "residue, modification mass"
        udtParamLabels(eInputParamLabelNames.Residue_PotentialModMass) = "residue, potential modification mass"
        udtParamLabels(eInputParamLabelNames.Residue_PotentialModMotif) = "residue, potential modification motif"
        udtParamLabels(eInputParamLabelNames.Refine_PotentialModMass) = "refine, potential modification mass"
        udtParamLabels(eInputParamLabelNames.Refine_PotentialModMotif) = "refine, potential modification motif"
        udtParamLabels(eInputParamLabelNames.Refine_PotentialNTerminusMods) = "refine, potential N-terminus modifications"
        udtParamLabels(eInputParamLabelNames.Refine_PotentialCTerminusMods) = "refine, potential C-terminus modifications"
        udtParamLabels(eInputParamLabelNames.Protein_NTerminal_ResidueModMass) = "protein, N-terminal residue modification mass"
        udtParamLabels(eInputParamLabelNames.Protein_CTerminal_ResidueModMass) = "protein, C-terminal residue modification mass"
        udtParamLabels(eInputParamLabelNames.Protein_Cleavage_NTerminalMassChange) = "protein, cleavage N-terminal mass change"
        udtParamLabels(eInputParamLabelNames.Protein_Cleavage_CTerminalMassChange) = "protein, cleavage C-terminal mass change"
        udtParamLabels(eInputParamLabelNames.Protein_Cleavage_Site) = "protein, cleavage site"
        udtParamLabels(eInputParamLabelNames.Refine_ModificationMass) = "refine, modification mass"
        udtParamLabels(eInputParamLabelNames.Scoring_Include_Reverse) = "scoring, include reverse"


        ' Make sure all of the text in udtParamLabels() is lowercase
        For intIndex = 0 To udtParamLabels.Length - 1
            udtParamLabels(intIndex) = udtParamLabels(intIndex).ToLower
        Next intIndex

        Do While objXMLReader.Read()

            XMLTextReaderSkipWhitespace(objXMLReader)
            If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

            If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
                Select Case objXMLReader.Name.ToLower
                    Case XTANDEM_XML_ELEMENT_NAME_NOTE
                        ' Read the note's type
                        strNoteType = XMLTextReaderGetAttributeValue(objXMLReader, "type", String.Empty)

                        If strNoteType = NOTE_TYPE_INPUT Then
                            ' Read the note's label and inner text
                            strNoteLabel = XMLTextReaderGetAttributeValue(objXMLReader, "label", String.Empty)

                            ' Need to advance the reader before calling XMLTextReaderGetInnerText
                            blnSuccess = objXMLReader.Read()
                            strValue = XMLTextReaderGetInnerText(objXMLReader)

                            If Not strValue Is Nothing Then
                                Select Case strNoteLabel.ToLower
                                    Case udtParamLabels(eInputParamLabelNames.Residue_StaticModMass)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.StaticMod, CInt(eInputParamLabelNames.Residue_StaticModMass), False, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Residue_PotentialModMass)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, CInt(eInputParamLabelNames.Residue_PotentialModMass), False, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Residue_PotentialModMotif)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, CInt(eInputParamLabelNames.Residue_PotentialModMotif), True, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Refine_PotentialModMass)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, CInt(eInputParamLabelNames.Refine_PotentialModMass), False, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Refine_PotentialModMotif)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, CInt(eInputParamLabelNames.Refine_PotentialModMotif), True, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Refine_PotentialNTerminusMods)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, CInt(eInputParamLabelNames.Refine_PotentialNTerminusMods), False, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Refine_PotentialCTerminusMods)
                                        blnSuccess = ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, CInt(eInputParamLabelNames.Refine_PotentialCTerminusMods), False, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Protein_NTerminal_ResidueModMass)
                                        blnSuccess = ParseXTandemInputParameterProteinTerminusMod(CInt(eInputParamLabelNames.Protein_NTerminal_ResidueModMass), True, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Protein_CTerminal_ResidueModMass)
                                        blnSuccess = ParseXTandemInputParameterProteinTerminusMod(CInt(eInputParamLabelNames.Protein_CTerminal_ResidueModMass), False, strValue, intModInfoCount, udtModInfo)

                                    Case udtParamLabels(eInputParamLabelNames.Protein_Cleavage_NTerminalMassChange)
                                        If clsPHRPBaseClass.IsNumber(strValue) Then
                                            mPeptideNTerminusMassChange = Double.Parse(strValue)
                                        End If
                                    Case udtParamLabels(eInputParamLabelNames.Protein_Cleavage_CTerminalMassChange)
                                        If clsPHRPBaseClass.IsNumber(strValue) Then
                                            mPeptideCTerminusMassChange = Double.Parse(strValue)
                                        End If

                                    Case udtParamLabels(eInputParamLabelNames.Protein_Cleavage_Site)
                                        ' In X!Tandem the LeftSpec and RightSpec values are separated by a vertical bar (CLEAVAGE_SPEC_SEP)
                                        ' Look for CLEAVAGE_SPEC_SEP in strValue
                                        intBarLoc = strValue.IndexOf(CLEAVAGE_SPEC_SEP)
                                        If intBarLoc > 0 And intBarLoc < strValue.Length - 1 Then
                                            strLeftSpec = strValue.Substring(0, intBarLoc)
                                            strRightSpec = strValue.Substring(intBarLoc + 1)

                                            ' Look for curly braces in strLeftSpec and strRightSpec
                                            ' In X!Tandem curly braces mean to not match a residue
                                            ' If found, change to standard RegEx notation, e.g. from {P} to [^P] 
                                            If strLeftSpec.IndexOf(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START) >= 0 Then
                                                strLeftSpec = strLeftSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START, "[^")
                                                strLeftSpec = strLeftSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END, "]")
                                            End If

                                            If strRightSpec.IndexOf(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START) >= 0 Then
                                                strRightSpec = strRightSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START, "[^")
                                                strRightSpec = strRightSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END, "]")
                                            End If

                                            With mEnzymeMatchSpec
                                                .LeftResidueRegEx = strLeftSpec
                                                .RightResidueRegEx = strRightSpec
                                            End With
                                        End If
                                    Case udtParamLabels(eInputParamLabelNames.Refine_ModificationMass)
                                        If Not strValue Is Nothing AndAlso strValue.Trim.Length > 0 Then
                                            blnStaticModsAreResetForRefinement = True
                                        End If
                                    Case udtParamLabels(eInputParamLabelNames.Scoring_Include_Reverse)
                                        If Not strValue Is Nothing AndAlso strValue.Trim.Length > 0 Then
                                            If strValue.Trim.ToLower = "yes" Then
                                                mLookForReverseSequenceTag = True
                                            End If
                                        End If
                                End Select
                            End If
                        End If

                    Case Else
                        ' Unknown/unneeded child node name; ignore it
                End Select

            ElseIf objXMLReader.NodeType = Xml.XmlNodeType.EndElement Then
                If objXMLReader.Depth = intParametersGroupDepth AndAlso objXMLReader.Name = XTANDEM_XML_ELEMENT_NAME_GROUP Then
                    ' Reached the end of this group
                    Exit Do
                End If
            End If
        Loop

        If intModInfoCount > 0 Then
            ' Validate that each of the mods in udtModInfo is present in mPeptideMods
            If intModInfoCount > 1 Then
                ' Sort udtModInfo
                Array.Sort(udtModInfo, 0, intModInfoCount, New ISearchOptionModificationInfoComparer)
            End If

            ' Before continuing, look for Static residue mods in udtModInfo
            ' If any are found, and if an identical dynamic residue mod is already present, then delete the static residue mod
            ' Additionally, if 	<note type="input" label="refine, modification mass">none</note> was prsent in the XTandem results file then
            '  auto update all static mods to dynamic mods since they are reset during refinement
            intIndex = 0
            Do While intIndex < intModInfoCount
                blnModDeleted = False
                If udtModInfo(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                    intIndexCompare = 0
                    Do While intIndexCompare < intModInfoCount
                        If intIndexCompare <> intIndex AndAlso udtModInfo(intIndexCompare).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod Then
                            ' See if this modification has a similar mass (within MASS_DIGITS_OF_PRECISION digits of precision)
                            If Math.Round(Math.Abs(udtModInfo(intIndexCompare).ModificationMass - udtModInfo(intIndex).ModificationMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION) = 0 Then
                                ' Matching mass
                                ' Compare .TargetResidues
                                If clsModificationDefinition.EquivalentTargetResidues(udtModInfo(intIndexCompare).TargetResidues, udtModInfo(intIndex).TargetResidues, True) Then
                                    ' Yes, the modifications match; delete the static version of the modification
                                    For intIndexCopy = intIndex To intModInfoCount - 2
                                        udtModInfo(intIndexCopy) = udtModInfo(intIndexCopy + 1)
                                    Next intIndexCopy
                                    intModInfoCount -= 1
                                    blnModDeleted = True
                                    Exit Do
                                End If
                            End If

                        End If
                        intIndexCompare += 1
                    Loop
                End If

                If Not blnModDeleted Then
                    If udtModInfo(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.StaticMod AndAlso blnStaticModsAreResetForRefinement Then
                        ' Update this static mod to be a dynamic mod
                        udtModInfo(intIndex).ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod
                    End If
                    intIndex += 1
                End If
            Loop

            For intIndex = 0 To intModInfoCount - 1
                With udtModInfo(intIndex)
                    mPeptideMods.VerifyModificationPresent(.ModificationMass, .TargetResidues, .ModificationType)
                End With
            Next intIndex
        End If

        ' In addition, verify that the standard refinement modifications are present in mPeptideMods
        mPeptideMods.AppendStandardRefinmentModifications()

    End Sub

    Private Sub ParseXTandemResultsFileReadDomainMods(ByRef objXMLReader As System.Xml.XmlTextReader, ByRef objSearchResult As clsSearchResultsXTandem, ByVal intDomainElementReaderDepth As Integer, ByVal blnUpdateModOccurrenceCounts As Boolean)

		Dim eResidueTerminusState As clsAminoAcidModInfo.eResidueTerminusStateConstants
        Dim chTargetResidue As Char
        Dim strValue As String

        Dim intResidueLocInPeptide As Integer

        Dim intModifiedResiduePosInProtein As Integer
        Dim dblModificationMass As Double

        ' Continue reading the XML file, loading the information 
        Do While objXMLReader.Read() AndAlso objXMLReader.Depth >= intDomainElementReaderDepth

            If Not objXMLReader.ReadState = Xml.ReadState.Interactive Then Exit Do

            If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
                Select Case objXMLReader.Name.ToLower
                    Case XTANDEM_XML_ELEMENT_NAME_AMINOACID
                        strValue = XMLTextReaderGetAttributeValue(objXMLReader, "type", String.Empty).Trim
                        If strValue.Length = 0 Then
                            chTargetResidue = Nothing
                        Else
                            chTargetResidue = strValue.Chars(0)
                        End If

                        intModifiedResiduePosInProtein = XMLTextReaderGetAttributeValue(objXMLReader, "at", 0)

                        If intModifiedResiduePosInProtein > 0 Then
                            dblModificationMass = XMLTextReaderGetAttributeValueDbl(objXMLReader, "modified", 0)

                            If dblModificationMass <> 0 Then
                                With objSearchResult
                                    intResidueLocInPeptide = intModifiedResiduePosInProtein - .PeptideLocInProteinStart + 1
                                    eResidueTerminusState = .DetermineResidueTerminusState(intResidueLocInPeptide)

                                    .SearchResultAddModification( _
                                                        dblModificationMass, _
                                                        chTargetResidue, _
                                                        intResidueLocInPeptide, _
                                                        eResidueTerminusState, _
                                                        blnUpdateModOccurrenceCounts)
                                End With
                            End If
                        End If
                End Select
            ElseIf objXMLReader.NodeType = Xml.XmlNodeType.EndElement Then
                If objXMLReader.Name = XTANDEM_XML_ELEMENT_NAME_DOMAIN Then
                    Exit Do
                End If
            End If
        Loop


    End Sub

	''' <summary>
	''' Main processing function
	''' </summary>
	''' <param name="strInputFilePath">X!Tandem results file</param>
	''' <param name="strOutputFolderPath">Output folder</param>
	''' <param name="strParameterFilePath">Parameter file</param>
	''' <returns>True if success, False if failure</returns>
    Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean
        
		Dim ioInputFile As System.IO.FileInfo

		Dim strXtandemXTFilePath As String

		Dim blnSuccess As Boolean

		If Not LoadParameterFileSettings(strParameterFilePath) Then
			SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, True)
			Return False
		End If

		Try
			If strInputFilePath Is Nothing OrElse strInputFilePath.Length = 0 Then
				SetErrorMessage("Input file name is empty")
				SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath)
			Else

				blnSuccess = ResetMassCorrectionTagsAndModificationDefinitions()
				If Not blnSuccess Then
					Exit Try
				End If

				MyBase.ResetProgress("Parsing " & System.IO.Path.GetFileName(strInputFilePath))

				If CleanupFilePaths(strInputFilePath, strOutputFolderPath) Then
					Try
						' Obtain the full path to the input file
						ioInputFile = New System.IO.FileInfo(strInputFilePath)

						' Define the output file name based on strInputFilePath
						' The name will be DatasetName_xt.txt

						strXtandemXTFilePath = System.IO.Path.GetFileName(MyBase.ReplaceFilenameSuffix(ioInputFile, ".txt"))
						strXtandemXTFilePath = System.IO.Path.Combine(strOutputFolderPath, strXtandemXTFilePath)
						blnSuccess = ParseXTandemResultsFile(ioInputFile.FullName, strXtandemXTFilePath, False)

						If blnSuccess AndAlso mCreateProteinModsFile Then

							' First create the MTS PepToProteinMap file using ioInputFile
							Dim lstSourcePHRPDataFiles As Generic.List(Of String) = New Generic.List(Of String)
							Dim strMTSPepToProteinMapFilePath As String = String.Empty

							lstSourcePHRPDataFiles.Add(strXtandemXTFilePath)

							strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(strInputFilePath, strOutputFolderPath, MTS:=True)

							If System.IO.File.Exists(strMTSPepToProteinMapFilePath) AndAlso mUseExistingMTSPepToProteinMapFile Then
								blnSuccess = True
							Else
								blnSuccess = MyBase.CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath)
							End If

							If blnSuccess Then
								' If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
								MyBase.ValidatePHRPReaderSupportFiles(strInputFilePath, strOutputFolderPath)

								' Now create the Protein Mods file
								blnSuccess = MyBase.CreateProteinModDetailsFile(strXtandemXTFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.XTandem)
							End If
						End If

					Catch ex As Exception
						SetErrorMessage("Error calling ParseXTandemResultsFile" & ex.message)
						SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
					End Try
				End If
			End If
		Catch ex As Exception
			SetErrorMessage("Error in ProcessFile:" & ex.Message)
			SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.UnspecifiedError)
		End Try

        Return blnSuccess

    End Function

    Private Sub SaveXTandemResultsFileEntry(ByRef objSearchResult As clsSearchResultsXTandem, ByRef swPeptideResultsFile As System.IO.StreamWriter)

        With objSearchResult
            ' Update .ResultID to the next available number
            .ResultID = mNextResultID
            mNextResultID += 1

            ' Write the results to the output file
			swPeptideResultsFile.WriteLine( _
				 .ResultID & SEP_CHAR & _
				 .GroupID & SEP_CHAR & _
				 .Scan & SEP_CHAR & _
				 .Charge & SEP_CHAR & _
				 .PeptideMH & SEP_CHAR & _
				 .PeptideHyperscore & SEP_CHAR & _
				 .PeptideExpectationValue & SEP_CHAR & _
				 .MultipleProteinCount & SEP_CHAR & _
				 .SequenceWithPrefixAndSuffix(True) & SEP_CHAR & _
				 Math.Round(.PeptideDeltaCn2, 4).ToString & SEP_CHAR & _
				 .PeptideYScore & SEP_CHAR & _
				 .PeptideYIons & SEP_CHAR & _
				 .PeptideBScore & SEP_CHAR & _
				 .PeptideBIons & SEP_CHAR & _
				 .PeptideDeltaMass & SEP_CHAR & _
				 .PeptideIntensity & SEP_CHAR & _
				 MyBase.NumToString(.PeptideDeltaMassCorrectedPpm, 4, True))
        End With

    End Sub

    Private Function TruncateProteinName(ByVal strProteinNameAndDescription As String) As String

        Dim intIndex As Integer
        Dim blnIsReversed As Boolean

        If mLookForReverseSequenceTag Then
            blnIsReversed = strProteinNameAndDescription.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR)
        Else
            blnIsReversed = False
        End If

        intIndex = strProteinNameAndDescription.IndexOf(" "c)
        If intIndex > 0 Then
            strProteinNameAndDescription = strProteinNameAndDescription.Substring(0, intIndex)
        End If

        If blnIsReversed Then
            Return strProteinNameAndDescription & REVERSED_PROTEIN_SEQUENCE_INDICATOR
        Else
            Return strProteinNameAndDescription
        End If

    End Function

    Private Sub UpdateSearchResultEnzymeAndTerminusInfo(ByRef objSearchResult As clsSearchResultsXTandem)
        With objSearchResult
            .SetEnzymeMatchSpec(mEnzymeMatchSpec)

            ' Update the N-Terminus and/or C-Terminus masses if those in the XML file are significantly different than the defaults
            If mPeptideNTerminusMassChange <> 0 Then
                .UpdatePeptideNTerminusMass(mPeptideNTerminusMassChange)
            End If

            If mPeptideCTerminusMassChange <> 0 Then
                .UpdatePeptideCTerminusMass(mPeptideCTerminusMassChange)
            End If
        End With
    End Sub

    Private Function XMLTextReaderGetAttributeValue(ByRef objXMLReader As System.Xml.XmlTextReader, ByVal strAttributeName As String, ByVal strValueIfMissing As String) As String
        objXMLReader.MoveToAttribute(strAttributeName)
        If objXMLReader.ReadAttributeValue() Then
            Return objXMLReader.Value
        Else
            Return String.Copy(strValueIfMissing)
        End If
    End Function

    Private Function XMLTextReaderGetAttributeValue(ByRef objXMLReader As System.Xml.XmlTextReader, ByVal strAttributeName As String, ByVal intValueIfMissing As Integer) As Integer
        objXMLReader.MoveToAttribute(strAttributeName)
        If objXMLReader.ReadAttributeValue() Then
            If clsPHRPBaseClass.IsNumber(objXMLReader.Value) Then
                Return CInt(objXMLReader.Value)
            Else
                Return intValueIfMissing
            End If
        Else
            Return intValueIfMissing
        End If
    End Function

    Private Function XMLTextReaderGetAttributeValueDbl(ByRef objXMLReader As System.Xml.XmlTextReader, ByVal strAttributeName As String, ByVal dblValueIfMissing As Double) As Double
        objXMLReader.MoveToAttribute(strAttributeName)
        If objXMLReader.ReadAttributeValue() Then
            If clsPHRPBaseClass.IsNumber(objXMLReader.Value) Then
                Return CDbl(objXMLReader.Value)
            Else
                Return dblValueIfMissing
            End If
        Else
            Return dblValueIfMissing
        End If
    End Function

    Private Function XMLTextReaderGetInnerText(ByRef objXMLReader As System.Xml.XmlTextReader) As String
        Dim strValue As String = String.Empty
        Dim blnSuccess As Boolean

        If objXMLReader.NodeType = Xml.XmlNodeType.Element Then
            ' Advance the reader so that we can read the value
            blnSuccess = objXMLReader.Read()
        Else
            blnSuccess = True
        End If

        If blnSuccess AndAlso Not objXMLReader.NodeType = Xml.XmlNodeType.Whitespace And objXMLReader.HasValue Then
            strValue = objXMLReader.Value
        End If

        Return strValue
    End Function

    Private Sub XMLTextReaderSkipWhitespace(ByRef objXMLReader As System.Xml.XmlTextReader)
        If objXMLReader.NodeType = Xml.XmlNodeType.Whitespace Then
            ' Whitspace; read the next node
            objXMLReader.Read()
        End If
    End Sub

End Class
