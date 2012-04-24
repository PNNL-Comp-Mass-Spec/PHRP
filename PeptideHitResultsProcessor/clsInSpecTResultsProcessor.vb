Option Strict On

' This class reads in an InSpecT results file (txt format) and creates 
' a tab-delimited text file with the data.  It will insert modification symbols
' into the peptide sequences for modified peptides.
'
' The modification definition information is determined from the InSpecT parameter file
'
' -------------------------------------------------------------------------------
' Written by John Sandoval for the Department of Energy (PNNL, Richland, WA)
' Program started August 12, 2008
'
' E-mail: john.sandoval@pnl.gov
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

Public Class clsInSpecTResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "August 18, 2011"
        InitializeLocalVariables()
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_INSPECT_FILE As String = "_inspect"
    Public Const FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING As String = "_PepToProtMap"

    Private Const INSPECT_SYN_FILE_MIN_COL_COUNT As Integer = 5

    Public Const N_TERMINUS_SYMBOL_INSPECT As String = "*."
    Public Const C_TERMINUS_SYMBOL_INSPECT As String = ".*"

    Private Const UNKNOWN_INSPECT_MOD_SYMBOL As Char = "?"c

    ' When writing the synopsis file, we keep data that passes any of these thresholds (thus, it's an OR comparison, not an AND comparison)
    ' pValue <= 0.2 Or TotalPRMScore >= 50 or FScore >= 0
    Public Const DEFAULT_SYN_FILE_PVALUE_THRESHOLD As Single = 0.2
    Public Const TOTALPRMSCORE_THRESHOLD As Single = 50
    Public Const FSCORE_THRESHOLD As Single = 0

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Private Const DTA_FILENAME_SCAN_NUMBER_REGEX As String = "(\d+)\.\d+\.\d+\.dta"
    Private Const INSPECT_NTERMINAL_MOD_MASS_REGEX As String = "^\+(\d+)"
    Private Const INSPECT_CTERMINAL_MOD_MASS_REGEX As String = "\+(\d+)$"

    Private Const REGEX_OPTIONS As Text.RegularExpressions.RegexOptions = Text.RegularExpressions.RegexOptions.Compiled Or Text.RegularExpressions.RegexOptions.Singleline Or Text.RegularExpressions.RegexOptions.IgnoreCase

    Private Const PHOS_MOD_NAME As String = "phos"
    Private Const PHOS_MOD_MASS As String = "79.9663"
    Private Const PHOS_MOD_RESIDUES As String = "STY"

    ' These columns correspond to the tab-delimited file created directly by Inspect
    Public Enum eInspectResultsFileColumns As Integer
        SpectrumFile = 0
        Scan = 1
        Annotation = 2
        Protein = 3
        Charge = 4
        MQScore = 5
        Length = 6
        TotalPRMScore = 7
        MedianPRMScore = 8
        FractionY = 9
        FractionB = 10
        Intensity = 11
        NTT = 12
        pvalue = 13
        FScore = 14
        DeltaScore = 15
        DeltaScoreOther = 16
        RecordNumber = 17
        DBFilePos = 18
        SpecFilePos = 19
        PrecursorMZ = 20
        PrecursorError = 21
    End Enum

    ' These columns correspond to the Synopsis and First-Hits files created by this class
    Protected Const InspectSynopsisFileColCount As Integer = 27
    Public Enum eInspectSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        Peptide = 2
        Protein = 3
        Charge = 4
        MQScore = 5
        Length = 6
        TotalPRMScore = 7
        MedianPRMScore = 8
        FractionY = 9
        FractionB = 10
        Intensity = 11
        NTT = 12
        PValue = 13
        FScore = 14
        DeltaScore = 15
        DeltaScoreOther = 16
        DeltaNormMQScore = 17                   ' Computed as Abs((MQScore(n) - MQScore(n+1)) / MQScore(n)); storing 0 for the lowest scoring result in each set. If MQScore(n) is 0, then also storing 0.   This value is not usable when MQScore(n) is <= 0, and should generally not be used when MQScore(n) is < 0.5
        DeltaNormTotalPRMScore = 18             ' Computed as Abs((TotalPRMScore(n) - TotalPRMScore(n+1)) / TotalPRMScore(n)); storing 0 for the lowest scoring result in each set.  If TotalPRMScore(n) is 0, then also storing 0.  This value is not usable when TotalPRMScore(n) is <= 0, and should generally not be used when TotalPRMScore(n) is < 0.5
        RankTotalPRMScore = 19                  ' Rank 1 means highest TotalPRMScore, 2 means next lower score, etc. (ties get the same rank)
        RankFScore = 20                         ' Rank 1 means highest FScore, 2 means next lower, etc. (ties get the same rank)
        MH = 21                                 ' Theoretical monoisotopic peptide mass (computed by PHRP); note that this is (M+H)+
        RecordNumber = 22
        DBFilePos = 23
        SpecFilePos = 24
        PrecursorMZ = 25
        PrecursorError = 26
    End Enum

    Protected Enum eInspectModType As Integer
        Unknown = 0
        DynamicMod = 1
        StaticMod = 2
        DynNTermPeptide = 3
        DynCTermPeptide = 4
    End Enum

    Protected Enum eFilteredOutputFileTypeConstants As Integer
        SynFile = 0
        FHTbyFScore = 1
        FHTbyTotalPRM = 2
    End Enum
#End Region

#Region "Structures"
    Protected Structure udtInspectSearchResultType
        Public SpectrumFile As String
        Public Scan As String
        Public ScanNum As Integer
        Public PeptideAnnotation As String
        Public Protein As String
        Public Charge As String
        Public ChargeNum As Short
        Public MQScore As String                ' Higher values are better scores; note that MQScore can be negative
        Public MQScoreNum As Single             ' Store the value of the string for quick reference when sorting
        Public Length As Integer
        Public TotalPRMScore As String          ' Higher values are better scores
        Public TotalPRMScoreNum As Single       ' We store the value of the string for quick reference when sorting
        Public MedianPRMScore As String
        Public FractionY As String
        Public FractionB As String
        Public Intensity As String
        Public NTT As Integer
        Public pValue As String                 ' Lower values are better scores
        Public PValueNum As Single              ' Store the value of the string for quick reference when sorting
        Public FScore As String                 ' Higher values are better scores
        Public FScoreNum As Single              ' Store the value of the string for quick reference when sorting
        Public DeltaScore As String
        Public DeltaScoreOther As String
        Public DeltaNormMQScore As Single
        Public DeltaNormTotalPRMScore As Single
        Public RankTotalPRMScore As Integer
        Public RankFScore As Integer
        Public MH As Double
        Public RecordNumber As String
        Public DBFilePos As String
        Public SpecFilePos As String
        Public PrecursorMZ As String
        Public PrecursorError As String         ' Precursor error in; units are m/z (NOT Daltons)
        Public DelMPPM As String                ' Computed by this application

        Public Sub Clear()
            SpectrumFile = String.Empty
            Scan = String.Empty
            ScanNum = 0
            PeptideAnnotation = String.Empty
            Protein = String.Empty
            Charge = String.Empty
            ChargeNum = 0
            MQScore = String.Empty
            MQScoreNum = 0
            Length = 0
            TotalPRMScore = String.Empty
            TotalPRMScoreNum = 0
            MedianPRMScore = String.Empty
            FractionY = String.Empty
            FractionB = String.Empty
            Intensity = String.Empty
            NTT = 0
            pValue = String.Empty
            PValueNum = 0
            FScore = String.Empty
            FScoreNum = 0
            DeltaScore = String.Empty
            DeltaScoreOther = String.Empty
            DeltaNormMQScore = 0
            DeltaNormTotalPRMScore = 0
            RankTotalPRMScore = 0
            RankFScore = 0
            MH = 0
            RecordNumber = String.Empty
            DBFilePos = String.Empty
            SpecFilePos = String.Empty
            PrecursorMZ = String.Empty
            PrecursorError = String.Empty
            DelMPPM = String.Empty
        End Sub
    End Structure

    Protected Structure udtModInfoType
        Public ModName As String            ' Mod names must be lower case, and 4 characters in length (or shorter)
        Public ModMass As String            ' Storing as a string since reading from a text file and writing to a text file
        Public Residues As String
        Public ModType As eInspectModType
        Public ModSymbol As String
    End Structure

    Protected Structure udtPepToProteinMappingType
        Public Peptide As String
        Public Protein As String
        Public ResidueStart As Integer
        Public ResidueEnd As Integer
    End Structure

#End Region

#Region "Classwide Variables"
    Protected mSortFHTandSynFiles As Boolean
#End Region

#Region "Properties"
    Public Property SortFHTandSynFiles() As Boolean
        Get
            Return mSortFHTandSynFiles
        End Get
        Set(ByVal value As Boolean)
            mSortFHTandSynFiles = value
        End Set
    End Property
#End Region

    Private Sub AddCurrentRecordToSearchResults(ByRef intCurrentScanResultsCount As Integer, _
                                                     ByRef udtSearchResultsCurrentScan() As udtInspectSearchResultType, _
                                                     ByRef udtSearchResult As udtInspectSearchResultType, _
                                                     ByRef strErrorLog As String)

        If intCurrentScanResultsCount >= udtSearchResultsCurrentScan.Length Then
            ReDim Preserve udtSearchResultsCurrentScan(udtSearchResultsCurrentScan.Length * 2 - 1)
        End If

        udtSearchResultsCurrentScan(intCurrentScanResultsCount) = udtSearchResult
        intCurrentScanResultsCount += 1

    End Sub

    Private Sub AddDynamicAndStaticResidueMods(ByRef objSearchResult As clsSearchResultsInSpecT, ByVal blnUpdateModOccurrenceCounts As Boolean)
        ' Step through .PeptideSequenceWithMods
        ' For each residue, check if a static mod is defined that affects that residue
        ' For each mod symbol, determine the modification and add to objSearchResult

        Dim intIndex As Integer, intModIndex As Integer
        Dim chChar As Char
        Dim objModificationDefinition As clsModificationDefinition

        Dim strSequence As String
        Dim chMostRecentLetter As Char
        Dim intResidueLocInPeptide As Integer

        chMostRecentLetter = "-"c
        intResidueLocInPeptide = 0

        With objSearchResult
            strSequence = .PeptideSequenceWithMods
            For intIndex = 0 To strSequence.Length - 1
                chChar = strSequence.Chars(intIndex)

                If Char.IsLetter(chChar) Then
                    chMostRecentLetter = chChar
                    intResidueLocInPeptide += 1

                    For intModIndex = 0 To mPeptideMods.ModificationCount - 1
                        If mPeptideMods.GetModificationTypeByIndex(intModIndex) = clsModificationDefinition.eModificationTypeConstants.StaticMod Then
                            objModificationDefinition = mPeptideMods.GetModificationByIndex(intModIndex)

                            If objModificationDefinition.TargetResiduesContain(chChar) Then
                                ' Match found; add this modification
                                .SearchResultAddModification(objModificationDefinition, chChar, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                            End If
                        End If
                    Next intModIndex
                ElseIf Char.IsLetter(chMostRecentLetter) Then
                    .SearchResultAddDynamicModification(chChar, chMostRecentLetter, intResidueLocInPeptide, .DetermineResidueTerminusState(intResidueLocInPeptide), blnUpdateModOccurrenceCounts)
                Else
                    ' We found a modification symbol but chMostRecentLetter is not a letter
                    ' Therefore, this modification symbol is at the beginning of the string; ignore the symbol
                End If

            Next intIndex
        End With
    End Sub

    ''' <summary>
    ''' Sorts the data by descending TotalPRMScore, than ranks each entry; in addition, computes normalized delta score (DeltaNorm) values
    ''' </summary>
    ''' <param name="udtSearchResultsCurrentScan"></param>
    ''' <param name="intCurrentScanResultsCount"></param>
    ''' <remarks></remarks>
    Private Sub AssignRankAndDeltaNormValues(ByRef udtSearchResultsCurrentScan() As udtInspectSearchResultType, ByVal intCurrentScanResultsCount As Integer)

        Const DeltaNormMQScore_If_Undefined As Single = 0
        Const DeltaNormTotalPRMScore_If_Undefined As Single = 0

        Static objSortScanChargeFScore As New InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc
        Static objSortScanChargeMQScore As New InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc
        Static objSortScanChargeTotalPRMDesc As New InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc

        Dim intIndex As Integer

        Dim intLastCharge As Integer
        Dim dblLastValue As Double

        Dim intCurrentRank As Integer


        ' Sort udtFilteredSearchResults by ascending scan, ascending charge, and descending RankFScore
        ' All of the data in udtSearchResultsCurrentScan should have the same scan number
        Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortScanChargeFScore)

        For intIndex = 0 To intCurrentScanResultsCount - 1
            If intIndex = 0 OrElse udtSearchResultsCurrentScan(intIndex).ChargeNum <> intLastCharge Then
                intLastCharge = udtSearchResultsCurrentScan(intIndex).ChargeNum
                dblLastValue = udtSearchResultsCurrentScan(intIndex).FScoreNum
                intCurrentRank = 1
            Else
                If udtSearchResultsCurrentScan(intIndex).FScoreNum <> dblLastValue Then
                    dblLastValue = udtSearchResultsCurrentScan(intIndex).FScoreNum
                    intCurrentRank += 1
                End If
            End If

            udtSearchResultsCurrentScan(intIndex).RankFScore = intCurrentRank
        Next intIndex


        ' Sort udtFilteredSearchResults by ascending scan, ascending charge, and descending MQScore (note that MQScore can be negative)
        ' All of the data in udtSearchResultsCurrentScan should have the same scan number
        Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortScanChargeMQScore)

        For intIndex = 0 To intCurrentScanResultsCount - 1
            If intIndex < intCurrentScanResultsCount - 1 AndAlso _
               udtSearchResultsCurrentScan(intIndex).ChargeNum = udtSearchResultsCurrentScan(intIndex + 1).ChargeNum Then
                udtSearchResultsCurrentScan(intIndex).DeltaNormMQScore = ComputeDeltaNormScore(udtSearchResultsCurrentScan(intIndex).MQScoreNum, udtSearchResultsCurrentScan(intIndex + 1).MQScoreNum, DeltaNormMQScore_If_Undefined)
            Else
                udtSearchResultsCurrentScan(intIndex).DeltaNormMQScore = 0
            End If
        Next intIndex


        ' Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending PValue
        ' All of the data in udtSearchResultsCurrentScan should have the same scan number
        Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortScanChargeTotalPRMDesc)

        For intIndex = 0 To intCurrentScanResultsCount - 1
            If intIndex = 0 OrElse udtSearchResultsCurrentScan(intIndex).ChargeNum <> intLastCharge Then
                intLastCharge = udtSearchResultsCurrentScan(intIndex).ChargeNum
                dblLastValue = udtSearchResultsCurrentScan(intIndex).TotalPRMScoreNum
                intCurrentRank = 1
            Else
                If udtSearchResultsCurrentScan(intIndex).TotalPRMScoreNum <> dblLastValue Then
                    dblLastValue = udtSearchResultsCurrentScan(intIndex).TotalPRMScoreNum
                    intCurrentRank += 1
                End If
            End If

            udtSearchResultsCurrentScan(intIndex).RankTotalPRMScore = intCurrentRank

            If intIndex < intCurrentScanResultsCount - 1 AndAlso _
               udtSearchResultsCurrentScan(intIndex).ChargeNum = udtSearchResultsCurrentScan(intIndex + 1).ChargeNum Then
                udtSearchResultsCurrentScan(intIndex).DeltaNormTotalPRMScore = ComputeDeltaNormScore(udtSearchResultsCurrentScan(intIndex).TotalPRMScoreNum, udtSearchResultsCurrentScan(intIndex + 1).TotalPRMScoreNum, DeltaNormTotalPRMScore_If_Undefined)
            Else
                udtSearchResultsCurrentScan(intIndex).DeltaNormTotalPRMScore = 0
            End If

        Next intIndex

    End Sub

    Private Function AddModificationsAndComputeMass(ByRef objSearchResult As clsSearchResultsInSpecT, ByVal blnUpdateModOccurrenceCounts As Boolean) As Boolean
        Const ALLOW_DUPLICATE_MOD_ON_TERMINUS As Boolean = True

        Dim blnSuccess As Boolean

        Try
            ' If any modifications of type IsotopicMod are defined, add them to the Search Result Mods now
            objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

            ' Parse .PeptideSequenceWithMods to determine the modified residues present
            AddDynamicAndStaticResidueMods(objSearchResult, blnUpdateModOccurrenceCounts)

            ' Add the protein and peptide terminus static mods (if defined and if the peptide is at a protein terminus)
            ' Since Inspect allows a terminal peptide residue to be modified twice, we'll allow that to happen,
            '  even though, biologically, that's typically not possible
            ' However, there are instances where this is possible, e.g. methylation of D or E on the C-terminus 
            '  (where two COOH groups are present)
            objSearchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts)

            ' Compute the monoisotopic mass for this peptide
            objSearchResult.ComputeMonoisotopicMass()

            ' Populate .PeptideModDescription
            objSearchResult.UpdateModDescription()

            blnSuccess = True
        Catch ex As Exception
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    ''' <summary>
    ''' This routine creates a first hits file or synopsis file from the output from InSpecT
    ''' </summary>
    ''' <param name="strInputFilePath"></param>
    ''' <param name="strOutputFilePath"></param>
    ''' <param name="udtInspectModInfo">Used to replace Mod text entries in the peptides with Mod Symbols; assumes each entryin </param>
    ''' <param name="eFilteredOutputFileType">Synopsis file or first hits file (sorting on various columns)</param>
    ''' <param name="blnResetMassCorrectionTagsAndModificationDefinitions"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function CreateFHTorSYNResultsFile(ByVal strInputFilePath As String, _
                                                 ByVal strOutputFilePath As String, _
                                                 ByRef udtInspectModInfo() As udtModInfoType, _
                                                 ByVal eFilteredOutputFileType As eFilteredOutputFileTypeConstants, _
                                                 Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean

        Dim intPreviousScan As Integer

        Dim strLineIn As String
        Dim protein As String = String.Empty

		Dim udtSearchResult As udtInspectSearchResultType = New udtInspectSearchResultType

        Dim intCurrentScanResultsCount As Integer
        Dim udtSearchResultsCurrentScan() As udtInspectSearchResultType

        Dim intFilteredSearchResultCount As Integer
        Dim udtFilteredSearchResults() As udtInspectSearchResultType

        Dim intResultsProcessed As Integer
        Dim intResultID As Integer = 0

        Dim blnSuccess As Boolean
        Dim blnValidSearchResult As Boolean

        Dim strErrorLog As String = String.Empty
        Dim objSortComparer As IComparer

        Try
            ' Initialize variables
            intPreviousScan = Int32.MinValue

            If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
                ' Writes the synopsis file, which writes every record with a p-value below a set threshold or a TotalPRMScore above a certain threshold
                objSortComparer = New InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc
            Else
                If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.FHTbyTotalPRM Then
                    ' Write the PRM first-hits file, which writes the record with the highest TotalPRMScore
                    objSortComparer = New InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc
                Else
                    ' eFilteredOutputFileTypeConstants.FHTbyFScore
                    ' Write the FScore first-hits file, which writes the record with the highest FScore
                    objSortComparer = New InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc
                End If
            End If

            Try
                ' Open the input file and parse it
                ' Initialize the stream reader and the stream Text writer
				Using srDataFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					Using swResultFile As System.IO.StreamWriter = New System.IO.StreamWriter(New System.IO.FileStream(strOutputFilePath, IO.FileMode.Create, IO.FileAccess.Write, IO.FileShare.Read))

						' Write the header line
						WriteSynFHTFileHeader(swResultFile, strErrorLog)

						strErrorLog = String.Empty
						intResultsProcessed = 0

						' Initialize array that will hold all of the records for a given scan
						intCurrentScanResultsCount = 0
						ReDim udtSearchResultsCurrentScan(9)

						' Initialize the array that will hold all of the records that will ultimately be written out to disk
						intFilteredSearchResultCount = 0
						ReDim udtFilteredSearchResults(999)

						' Parse the input file
						Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing
							strLineIn = srDataFile.ReadLine
							If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then
								' Initialize udtSearchResult
								udtSearchResult.Clear()

								blnValidSearchResult = ParseInspectResultsFileEntry(strLineIn, udtInspectModInfo, udtSearchResult, strErrorLog, intResultsProcessed)

								If blnValidSearchResult Then
									If intPreviousScan <> Int32.MinValue AndAlso intPreviousScan <> udtSearchResult.ScanNum Then
										' New scan encountered; sort and filter the data in udtSearchResultsCurrentScan, then call StoreTopFHTMatch or StoreSynMatches
										If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
											StoreSynMatches(swResultFile, intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog, objSortComparer)
										Else
											StoreTopFHTMatch(swResultFile, intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog, objSortComparer)
										End If
										intCurrentScanResultsCount = 0
									End If

									AddCurrentRecordToSearchResults(intCurrentScanResultsCount, udtSearchResultsCurrentScan, udtSearchResult, strErrorLog)

									intPreviousScan = udtSearchResult.ScanNum
								End If

								' Update the progress
								UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))
								intResultsProcessed += 1
							End If
						Loop


						' Store the last record
						If intCurrentScanResultsCount > 0 Then
							If eFilteredOutputFileType = eFilteredOutputFileTypeConstants.SynFile Then
								StoreSynMatches(swResultFile, intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog, objSortComparer)
							Else
								StoreTopFHTMatch(swResultFile, intResultID, intCurrentScanResultsCount, udtSearchResultsCurrentScan, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog, objSortComparer)
							End If
							intCurrentScanResultsCount = 0
						End If

						If mSortFHTandSynFiles Then
							' Sort the data in udtFilteredSearchResults then write out to disk
							SortAndWriteFilteredSearchResults(swResultFile, intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
						End If

					End Using
				End Using

				' Inform the user if any errors occurred
				If strErrorLog.Length > 0 Then
					SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
				End If

				blnSuccess = True

			Catch ex As Exception
				SetErrorMessage(ex.Message)
				SetErrorCode(clsPHRPBaseClass.ePHRPErrorCodes.ErrorReadingInputFile)
				blnSuccess = False
			End Try
        Catch ex As Exception
            SetErrorMessage(ex.Message)
            SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles)
            blnSuccess = False
        End Try

        Return blnSuccess

    End Function

    Protected Function ComputeDelMCorrectedPPM(ByVal dblPrecursorErrorDa As Double, ByVal dblPrecursorMonoMass As Double, _
                                               ByVal dblPeptideMonoisotopicMass As Double, _
                                               ByVal blnAdjustPrecursorMassForC13 As Boolean) As Double

        Dim dblPeptideDeltaMassCorrectedPpm As Double

		dblPeptideDeltaMassCorrectedPpm = clsSearchResultsBaseClass.ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, blnAdjustPrecursorMassForC13, dblPeptideMonoisotopicMass)

        Return dblPeptideDeltaMassCorrectedPpm

    End Function

    Private Function ComputeDeltaNormScore(ByVal sngCurrentScore As Single, ByVal sngNextScore As Single, ByVal sngValueIfCurrentScoreZero As Single) As Single
        Try
            If sngCurrentScore <> 0 Then
                Return Math.Abs((sngCurrentScore - sngNextScore) / sngCurrentScore)
            Else
                Return sngValueIfCurrentScoreZero
            End If
        Catch ex As Exception
            Return sngValueIfCurrentScoreZero
        End Try
    End Function

    Private Function ComputePeptideMHFromPrecursorInfo(ByVal strPrecursorMZ As String, ByVal strPrecursorError As String, ByVal strCharge As String) As Double
        ' Compute the theoretical peptide MH using the precursor m/z value and the precursor error values

        Dim dblPrecursorMZ As Double
        Dim dblPrecursorError As Double
        Dim intCharge As Integer
        Dim dblPeptideMH As Double = 0

        If strPrecursorMZ Is Nothing OrElse strPrecursorMZ.Length = 0 OrElse strPrecursorMZ = "0" Then
            ' Precursor m/z is undefined; cannot continue
            dblPeptideMH = 0
        ElseIf strPrecursorError Is Nothing OrElse strPrecursorError.Length = 0 Then
            ' Precursor error is undefined; cannot continue
            dblPeptideMH = 0
        Else
            intCharge = CIntSafe(strCharge, -1)

            If intCharge >= 1 Then
                If Double.TryParse(strPrecursorMZ, dblPrecursorMZ) Then
                    If Double.TryParse(strPrecursorError, dblPrecursorError) Then
                        ' Note: the October 2008 version of Inspect uses an Absolute Value function when computing the PrecursorError; the version used by PNNL does not use Absolute Value
                        ' Note: switched to compute (M+H)+ in August 2011; prior to this, we were computing uncharged monoisotopic mass
                        dblPeptideMH = (dblPrecursorMZ - dblPrecursorError) * intCharge - (intCharge - 1) * clsPeptideMassCalculator.MASS_PROTON
                    End If
                End If
            End If
        End If

        Return dblPeptideMH

    End Function

    Private Function CIntSafe(ByVal strValue As String, ByVal intDefaultValue As Integer) As Integer
        Try
            ' Note: Integer.Parse() fails if strValue contains a decimal point, even if it is "8.000"
            ' Thus, we're using CInt() instead
            Return CInt(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return intDefaultValue
        End Try
    End Function

    Private Function CSngSafe(ByVal strValue As String, ByVal sngDefaultValue As Single) As Single
        Try
            Return Single.Parse(strValue)
        Catch ex As Exception
            ' Error converting strValue to a number; return the default
            Return sngDefaultValue
        End Try
    End Function

    Private Function ExtractModInfoFromInspectParamFile(ByVal strInspectParameterFilePath As String, ByRef udtModList() As udtModInfoType) As Boolean

		Dim strLineIn As String
        Dim strSplitLine As String()

        Dim intModCount As Integer
        Dim intUnnamedModID As Integer

        Dim blnSuccess As Boolean = False

        Try
            ' Initialize udtModList and intUnnamedModID
            intModCount = 0
            ReDim udtModList(-1)

            intUnnamedModID = 0

            If strInspectParameterFilePath Is Nothing OrElse strInspectParameterFilePath.Length = 0 Then
                SetErrorMessage("Inspect Parameter File name not defined; unable to extract mod info")
                SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
                Return False
            End If

            ' Read the contents of the inspect parameter file
			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strInspectParameterFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

				Do While srInFile.Peek <> -1
					strLineIn = srInFile.ReadLine

					strLineIn = strLineIn.Trim

					If strLineIn.Length > 0 Then

						If strLineIn.Chars(0) = "#"c Then
							' Comment line; skip it
						ElseIf strLineIn.ToLower.StartsWith("mod") Then
							' Modification definition line

							' Split the line on commas
							strSplitLine = strLineIn.Split(","c)

							If strSplitLine.Length >= 3 AndAlso strSplitLine(0).ToLower.Trim = "mod" Then
								If udtModList.Length = 0 Then
									ReDim udtModList(0)
								ElseIf intModCount >= udtModList.Length Then
									ReDim Preserve udtModList(udtModList.Length * 2 - 1)
								End If

								With udtModList(intModCount)
									.ModMass = strSplitLine(1)
									.Residues = strSplitLine(2)
									.ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL

									If strSplitLine.Length >= 4 Then
										Select Case strSplitLine(3).ToLower
											Case "opt"
												.ModType = eInspectModType.DynamicMod
											Case "fix"
												.ModType = eInspectModType.StaticMod
											Case "nterminal"
												.ModType = eInspectModType.DynNTermPeptide
											Case "cterminal"
												.ModType = eInspectModType.DynCTermPeptide
											Case Else
												SetErrorMessage("Warning: Unrecognized Mod Type in the Inspect parameter file")
												.ModType = eInspectModType.DynamicMod
										End Select
									Else
										' Assume dynamic if not specifed
										.ModType = eInspectModType.DynamicMod
									End If

									If strSplitLine.Length >= 5 Then
										.ModName = strSplitLine(4).ToLower
										If .ModName.Length > 4 Then
											' Only keep the first 4 characters of the modification name
											.ModName = .ModName.Substring(0, 4)
										End If
									Else
										intUnnamedModID += 1
										.ModName = "UnnamedMod" & intUnnamedModID.ToString
									End If

									' Check for phosphorylation
									' Inspect requires that it be defined in the parameter file as: mod,80,STY,opt,phosphorylation
									'  However, we want to use the more precise mass of 79.9663
									If .ModName = PHOS_MOD_NAME.ToLower And .ModMass = "80" Then
										.ModMass = PHOS_MOD_MASS
									End If

								End With

								intModCount += 1
							End If
						End If
					End If
				Loop

			End Using

			' Shrink udtModList to the appropriate length
			ReDim Preserve udtModList(intModCount - 1)

			Console.WriteLine()

			blnSuccess = True

		Catch ex As Exception
			SetErrorMessage("Error reading the Inspect parameter file (" & System.IO.Path.GetFileName(strInspectParameterFilePath) & "): " & ex.Message)
			SetErrorCode(ePHRPErrorCodes.ErrorReadingModificationDefinitionsFile)
			blnSuccess = False	
		End Try

        Return blnSuccess

    End Function

    Private Function ExtractScanNumFromDTAName(ByVal spectrumFile As String) As String
        Static reScanNumberRegEx As New System.Text.RegularExpressions.Regex(DTA_FILENAME_SCAN_NUMBER_REGEX, REGEX_OPTIONS)

        Dim scanNum As String = String.Empty

        ' See if strValue resembles a .Dta file name
        ' For example, "MyDataset.300.300.2.dta"

        Try
            With reScanNumberRegEx.Match(spectrumFile)
                If .Success AndAlso .Groups.Count > 1 Then
                    scanNum = .Groups(1).Value
                End If
            End With
        Catch ex As Exception
            ' Ignore errors here
            scanNum = "0"
        End Try

        Return scanNum

    End Function

    Private Sub InitializeLocalVariables()
        mSortFHTandSynFiles = True
    End Sub

    ''' <summary>
    ''' Load the PeptideToProteinMap information; in addition, writes out an updated _inspect_PepToProtMap.txt file with the new mod symbols and corrected terminii symbols
    ''' </summary>
    ''' <param name="strPepToProteinMapFilePath"></param>
    ''' <param name="strOutputFolderPath"></param>
    ''' <param name="udtInspectModInfo"></param>
    ''' <param name="udtPepToProteinMapping"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function LoadPeptideToProteinMapInfo(ByVal strPepToProteinMapFilePath As String, _
                                                   ByVal strOutputFolderPath As String, _
                                                   ByRef udtInspectModInfo() As udtModInfoType, _
                                                   ByRef udtPepToProteinMapping() As udtPepToProteinMappingType) As Boolean

        Dim strMTSPepToProteinMapFilePath As String = String.Empty
        Dim strMTSCompatiblePeptide As String

        Dim strLineIn As String
        Dim strSplitLine As String()

        Dim intLinesRead As Integer
        Dim intValue As Integer

        Dim intProteinCount As Integer
        Dim blnSuccess As Boolean = False

        Try
            ' Initialize udtPepToProteinMapping
            If strPepToProteinMapFilePath Is Nothing OrElse strPepToProteinMapFilePath.Length = 0 Then
                SetErrorMessage("Warning: PepToProteinMap file is not defined")
                Return False
            ElseIf Not System.IO.File.Exists(strPepToProteinMapFilePath) Then
                SetErrorMessage("Warning: PepToProteinMap file does not exist")
                Return False
            End If

            intProteinCount = 0
            ReDim udtPepToProteinMapping(999)

			' Open strProteinToPeptideMappingFilePath for reading
		
			Using srInFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strPepToProteinMapFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

				strMTSPepToProteinMapFilePath = System.IO.Path.Combine(strOutputFolderPath, System.IO.Path.GetFileNameWithoutExtension(strPepToProteinMapFilePath) & "MTS.txt")
				Using swOutFile As System.IO.StreamWriter = New System.IO.StreamWriter(New System.IO.FileStream(strMTSPepToProteinMapFilePath, IO.FileMode.Create, IO.FileAccess.Write, IO.FileShare.Read))

					intLinesRead = 0
					Do While srInFile.Peek <> -1
						strLineIn = srInFile.ReadLine

						strLineIn = strLineIn.Trim

						If strLineIn.Length > 0 Then

							' Split the line on tabs
							strSplitLine = strLineIn.Split(ControlChars.Tab)

							If strSplitLine.Length >= 4 Then

								If intLinesRead = 0 AndAlso Not Integer.TryParse(strSplitLine(2), intValue) Then
									' Header line
									swOutFile.WriteLine(strLineIn)
								Else
									If intProteinCount >= udtPepToProteinMapping.Length Then
										ReDim Preserve udtPepToProteinMapping(udtPepToProteinMapping.Length * 2 - 1)
									End If

									' Replace any mod text names in the peptide sequence with the appropriate mod symbols
									' In addition, replace the * terminus symbols with dashes
									strMTSCompatiblePeptide = ReplaceInspectModTextWithSymbol(ReplaceTerminus(strSplitLine(0)), udtInspectModInfo)

									With udtPepToProteinMapping(intProteinCount)
										.Peptide = strMTSCompatiblePeptide
										.Protein = String.Copy(strSplitLine(1))
										Integer.TryParse(strSplitLine(2), .ResidueStart)
										Integer.TryParse(strSplitLine(3), .ResidueEnd)
									End With


									swOutFile.WriteLine(strMTSCompatiblePeptide & ControlChars.Tab & _
										 strSplitLine(1) & ControlChars.Tab & _
										 strSplitLine(2) & ControlChars.Tab & _
										 strSplitLine(3))

									intProteinCount += 1
								End If
							End If
						End If

					Loop
				End Using
			End Using

			' Shrink udtPepToProteinMapping to the appropriate length
			ReDim Preserve udtPepToProteinMapping(intProteinCount - 1)

			Console.WriteLine()

			blnSuccess = True

		Catch ex As Exception
			SetErrorMessage("Error reading the Peptide to Protein Map File (" & System.IO.Path.GetFileName(strPepToProteinMapFilePath) & ") " & _
				"and writing the new map file (" & System.IO.Path.GetFileName(strMTSPepToProteinMapFilePath) & "): " & ex.Message)
			SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile)
			blnSuccess = False		
		End Try

        Return blnSuccess

    End Function

    Private Function ParseInspectSynFileHeaderLine(ByVal strLineIn As String, _
                                                   ByRef intColumnMapping() As Integer) As Boolean

        ' Parse the header line

        Dim strSplitLine() As String
        Dim eResultFileColumn As eInspectSynFileColumns
        Dim lstColumnNames As System.Collections.Generic.SortedDictionary(Of String, eInspectSynFileColumns)
        lstColumnNames = New System.Collections.Generic.SortedDictionary(Of String, eInspectSynFileColumns)(StringComparer.CurrentCultureIgnoreCase)

        ReDim intColumnMapping(InspectSynopsisFileColCount - 1)

        lstColumnNames.Add("ResultID", eInspectSynFileColumns.ResultID)
        lstColumnNames.Add("Scan", eInspectSynFileColumns.Scan)
        lstColumnNames.Add("Peptide", eInspectSynFileColumns.Peptide)
        lstColumnNames.Add("Protein", eInspectSynFileColumns.Protein)
        lstColumnNames.Add("Charge", eInspectSynFileColumns.Charge)
        lstColumnNames.Add("MQScore", eInspectSynFileColumns.MQScore)
        lstColumnNames.Add("Length", eInspectSynFileColumns.Length)
        lstColumnNames.Add("TotalPRMScore", eInspectSynFileColumns.TotalPRMScore)
        lstColumnNames.Add("MedianPRMScore", eInspectSynFileColumns.MedianPRMScore)
        lstColumnNames.Add("FractionY", eInspectSynFileColumns.FractionY)
        lstColumnNames.Add("FractionB", eInspectSynFileColumns.FractionB)
        lstColumnNames.Add("Intensity", eInspectSynFileColumns.Intensity)
        lstColumnNames.Add("NTT", eInspectSynFileColumns.NTT)
        lstColumnNames.Add("PValue", eInspectSynFileColumns.PValue)
        lstColumnNames.Add("FScore", eInspectSynFileColumns.FScore)
        lstColumnNames.Add("DeltaScore", eInspectSynFileColumns.DeltaScore)
        lstColumnNames.Add("DeltaScoreOther", eInspectSynFileColumns.DeltaScoreOther)
        lstColumnNames.Add("DeltaNormMQScore", eInspectSynFileColumns.DeltaNormMQScore)
        lstColumnNames.Add("DeltaNormTotalPRMScore", eInspectSynFileColumns.DeltaNormTotalPRMScore)
        lstColumnNames.Add("RankTotalPRMScore", eInspectSynFileColumns.RankTotalPRMScore)
        lstColumnNames.Add("RankFScore", eInspectSynFileColumns.RankFScore)
        lstColumnNames.Add("MH", eInspectSynFileColumns.MH)
        lstColumnNames.Add("RecordNumber", eInspectSynFileColumns.RecordNumber)
        lstColumnNames.Add("DBFilePos", eInspectSynFileColumns.DBFilePos)
        lstColumnNames.Add("SpecFilePos", eInspectSynFileColumns.SpecFilePos)
        lstColumnNames.Add("PrecursorMZ", eInspectSynFileColumns.PrecursorMZ)
        lstColumnNames.Add("PrecursorError", eInspectSynFileColumns.PrecursorError)

        Try
            ' Initialize each entry in intColumnMapping to -1
            For intIndex As Integer = 0 To intColumnMapping.Length - 1
                intColumnMapping(intIndex) = -1
            Next

            strSplitLine = strLineIn.Split(ControlChars.Tab)
            For intIndex As Integer = 0 To strSplitLine.Length - 1
                If lstColumnNames.TryGetValue(strSplitLine(intIndex), eResultFileColumn) Then
                    ' Recognized column name; update intColumnMapping
                    intColumnMapping(eResultFileColumn) = intIndex
                End If
            Next

        Catch ex As Exception
            SetErrorMessage("Error parsing header in Inspect synopsis file: " & ex.Message)
            Return False
        End Try

        Return True

    End Function

    Protected Function ParseInSpectSynopsisFile(ByVal strInputFilePath As String, ByRef udtPepToProteinMapping() As udtPepToProteinMappingType, Optional ByVal blnResetMassCorrectionTagsAndModificationDefinitions As Boolean = True) As Boolean
        ' Warning: This function does not call LoadParameterFile; you should typically call ProcessFile rather than calling this function

		Dim strPreviousTotalPRMScore As String

        ' Note that Inspect synopsis files are normally sorted on TotalPRMScore descending
        ' In order to prevent duplicate entries from being made to the ResultToSeqMap file (for the same peptide in the same scan),
        '  we will keep track of the scan, charge, and peptide information parsed for each unique TotalPRMScore encountered

        Dim htPeptidesFoundForTotalPRMScoreLevel As Hashtable

        Dim strKey As String

        Dim strLineIn As String
        Dim strModificationSummaryFilePath As String

        Dim objSearchResult As clsSearchResultsInSpecT

        Dim intResultsProcessed As Integer
        Dim intPepToProteinMapIndex As Integer

        Dim strCurrentPeptideWithMods As String = String.Empty
        Dim strCurrentProtein As String

        Dim blnSuccess As Boolean
        Dim blnDataLine As Boolean
        Dim blnValidSearchResult As Boolean
        Dim blnFirstMatchForGroup As Boolean

        Dim blnHeaderParsed As Boolean
        Dim intColumnMapping() As Integer = Nothing

        Dim strErrorLog As String = String.Empty

        Dim objPeptideSearchComparer As PepToProteinMappingPeptideSearchComparer

        Try
            ' Possibly reset the mass correction tags and Mod Definitions
            If blnResetMassCorrectionTagsAndModificationDefinitions Then
                ResetMassCorrectionTagsAndModificationDefinitions()
            End If

            ' Reset .OccurrenceCount
            mPeptideMods.ResetOccurrenceCountStats()

            ' Initialize objSearchResult
            objSearchResult = New clsSearchResultsInSpecT(mPeptideMods)

            ' Initialize htPeptidesFoundForTotalPRMScoreLevel
            htPeptidesFoundForTotalPRMScoreLevel = New Hashtable
            strPreviousTotalPRMScore = String.Empty

            ' Initialize objPeptideSearchComparer
            objPeptideSearchComparer = New PepToProteinMappingPeptideSearchComparer

            ' Assure that udtPepToProteinMapping is sorted on peptide
            If udtPepToProteinMapping.Length > 1 Then
                Array.Sort(udtPepToProteinMapping, New PepToProteinMappingComparer)
            End If

            Try
                UpdateSearchResultEnzymeAndTerminusInfo(objSearchResult)

                ' Open the input file and parse it
                ' Initialize the stream reader
				Using srDataFile As System.IO.StreamReader = New System.IO.StreamReader(New System.IO.FileStream(strInputFilePath, IO.FileMode.Open, IO.FileAccess.Read, IO.FileShare.Read))

					strErrorLog = String.Empty
					intResultsProcessed = 0
					blnHeaderParsed = False

					' Create the output files
					blnSuccess = MyBase.InitializeSequenceOutputFiles(strInputFilePath)

					' Parse the input file
					Do While srDataFile.Peek >= 0 And Not MyBase.AbortProcessing

						strLineIn = srDataFile.ReadLine
						If Not strLineIn Is Nothing AndAlso strLineIn.Trim.Length > 0 Then

							blnDataLine = True

							If Not blnHeaderParsed Then
								blnSuccess = ParseInspectSynFileHeaderLine(strLineIn, intColumnMapping)
								If blnSuccess Then
									blnDataLine = False
								Else
									' Error parsing header; assume this is a data line
									blnDataLine = True
								End If
								blnHeaderParsed = True
							End If

							If blnDataLine Then
								blnValidSearchResult = ParseInSpectSynFileEntry(strLineIn, intColumnMapping, objSearchResult, strErrorLog, strCurrentPeptideWithMods)
							Else
								blnValidSearchResult = False
							End If

							If blnValidSearchResult Then
								strKey = objSearchResult.PeptideSequenceWithMods & "_" & objSearchResult.Scan & "_" & objSearchResult.Charge

								If objSearchResult.TotalPRMScore = strPreviousTotalPRMScore Then
									' New result has the same TotalPRMScore as the previous result
									' See if htPeptidesFoundForTotalPRMScoreLevel contains the peptide, scan and charge

									If htPeptidesFoundForTotalPRMScoreLevel.ContainsKey(strKey) Then
										blnFirstMatchForGroup = False
									Else
										htPeptidesFoundForTotalPRMScoreLevel.Add(strKey, 1)
										blnFirstMatchForGroup = True
									End If

								Else
									' New TotalPRMScore
									' Reset htPeptidesFoundForScan
									htPeptidesFoundForTotalPRMScoreLevel.Clear()

									' Update strPreviousTotalPRMScore
									strPreviousTotalPRMScore = objSearchResult.TotalPRMScore

									' Append a new entry to htPeptidesFoundForScan
									htPeptidesFoundForTotalPRMScoreLevel.Add(strKey, 1)
									blnFirstMatchForGroup = True
								End If


								blnSuccess = AddModificationsAndComputeMass(objSearchResult, blnFirstMatchForGroup)
								If Not blnSuccess Then
									If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
										strErrorLog &= "Error adding modifications to sequence at RowIndex '" & objSearchResult.ResultID & "'" & ControlChars.NewLine
									End If
								End If
								MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResult, clsSearchResultsBaseClass), blnFirstMatchForGroup)


								If udtPepToProteinMapping.Length > 0 Then
									' Add the additional proteins for this peptide

									' Use binary search to find this peptide in udtPepToProteinMapping
									intPepToProteinMapIndex = Array.BinarySearch(udtPepToProteinMapping, strCurrentPeptideWithMods, objPeptideSearchComparer)

									If intPepToProteinMapIndex >= 0 Then
										' Step Backward until the first match is found
										Do While intPepToProteinMapIndex > 0 AndAlso udtPepToProteinMapping(intPepToProteinMapIndex - 1).Peptide = strCurrentPeptideWithMods
											intPepToProteinMapIndex -= 1
										Loop

										' Call MyBase.SaveResultsFileEntrySeqInfo for each entry in udtPepToProteinMapping() for peptide , skipping objSearchResult.ProteinName
										strCurrentProtein = String.Copy(objSearchResult.ProteinName)
										Do
											If udtPepToProteinMapping(intPepToProteinMapIndex).Protein <> strCurrentProtein Then
												objSearchResult.ProteinName = String.Copy(udtPepToProteinMapping(intPepToProteinMapIndex).Protein)
												MyBase.SaveResultsFileEntrySeqInfo(DirectCast(objSearchResult, clsSearchResultsBaseClass), False)
											End If

											intPepToProteinMapIndex += 1
										Loop While intPepToProteinMapIndex < udtPepToProteinMapping.Length AndAlso strCurrentPeptideWithMods = udtPepToProteinMapping(intPepToProteinMapIndex).Peptide
									Else
										' Match not found; this is unexpected
										Console.WriteLine("Warning: no match for '" & strCurrentPeptideWithMods & "' in udtPepToProteinMapping")
									End If
								End If

							End If

							' Update the progress
							UpdateProgress(CSng(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100))

							intResultsProcessed += 1
						End If
					Loop
				End Using

				If mCreateModificationSummaryFile Then
					' Create the modification summary file
					strModificationSummaryFilePath = MyBase.ReplaceFilenameSuffix(strInputFilePath, "", FILENAME_SUFFIX_MOD_SUMMARY)
					SaveModificationSummaryFile(strModificationSummaryFilePath)
				End If

				' Inform the user if any errors occurred
				If strErrorLog.Length > 0 Then
					SetErrorMessage("Invalid Lines: " & ControlChars.NewLine & strErrorLog)
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

    Private Function ParseInspectResultsFileEntry(ByRef strLineIn As String, _
                                                  ByRef udtInspectModInfo() As udtModInfoType, _
                                                  ByRef udtSearchResult As udtInspectSearchResultType, _
                                                  ByRef strErrorLog As String, _
                                                  ByVal intResultsProcessed As Integer) As Boolean

        ' Parses an entry from the Inspect results file
        ' The expected header line is:
        ' #SpectrumFile	Scan#	Annotation	Protein	Charge	MQScore	Length	TotalPRMScore	MedianPRMScore	FractionY	FractionB	Intensity	NTT	p-value	F-Score	DeltaScore	DeltaScoreOther	RecordNumber	DBFilePos	SpecFilePos PrecursorMZ	PrecursorError DelM_PPM

		Dim strSplitLine() As String = Nothing

        Dim blnValidSearchResult As Boolean
        
        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset udtSearchResult
            udtSearchResult.Clear()

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)

            If strSplitLine.Length >= 15 Then

                If intResultsProcessed = 0 Then
                    ' This is the first line of the file; it may be a header row
                    ' Determine this by seeing if any of the first three columns contains a number
                    If Not (clsPHRPBaseClass.IsNumber(strSplitLine(0)) OrElse _
                            clsPHRPBaseClass.IsNumber(strSplitLine(1)) OrElse _
                            clsPHRPBaseClass.IsNumber(strSplitLine(2))) Then
                        ' This is a header line; ignore it
                        blnValidSearchResult = False
                        Exit Try
                    End If
                End If

                With udtSearchResult
                    .SpectrumFile = strSplitLine(eInspectResultsFileColumns.SpectrumFile)
                    If strSplitLine(eInspectResultsFileColumns.Scan) = "0" Then
                        .Scan = ExtractScanNumFromDTAName(.SpectrumFile)
                    Else
                        .Scan = strSplitLine(eInspectResultsFileColumns.Scan)
                    End If
                    .ScanNum = CIntSafe(.Scan, 0)

                    ' Replace any mod text names in the peptide sequence with the appropriate mod symbols
                    ' In addition, replace the * terminus symbols with dashes
                    .PeptideAnnotation = ReplaceInspectModTextWithSymbol(ReplaceTerminus(strSplitLine(eInspectResultsFileColumns.Annotation)), udtInspectModInfo)
                    .Protein = TruncateProteinName(strSplitLine(eInspectResultsFileColumns.Protein))

                    .Charge = strSplitLine(eInspectResultsFileColumns.Charge)
                    .ChargeNum = CShort(CIntSafe(.Charge, 0))

                    .MQScore = strSplitLine(eInspectResultsFileColumns.MQScore)
                    .MQScoreNum = CSngSafe(.MQScore, 0)

                    .Length = CIntSafe(strSplitLine(eInspectResultsFileColumns.Length), 0)

                    .TotalPRMScore = strSplitLine(eInspectResultsFileColumns.TotalPRMScore)
                    .TotalPRMScoreNum = CSngSafe(.TotalPRMScore, 0)

                    .MedianPRMScore = strSplitLine(eInspectResultsFileColumns.MedianPRMScore)
                    .FractionY = RemoveExtraneousDigits(strSplitLine(eInspectResultsFileColumns.FractionY))
                    .FractionB = RemoveExtraneousDigits(strSplitLine(eInspectResultsFileColumns.FractionB))
                    .Intensity = strSplitLine(eInspectResultsFileColumns.Intensity)
                    .NTT = CIntSafe(strSplitLine(eInspectResultsFileColumns.NTT), 0)

                    .pValue = RemoveExtraneousDigits(strSplitLine(eInspectResultsFileColumns.pvalue))
                    .PValueNum = CSngSafe(.pValue, 0)

                    .FScore = strSplitLine(eInspectResultsFileColumns.FScore)
                    .FScoreNum = CSngSafe(.FScore, 0)

                    .DeltaScore = strSplitLine(eInspectResultsFileColumns.DeltaScore)
                    .DeltaScoreOther = strSplitLine(eInspectResultsFileColumns.DeltaScoreOther)

                    .RecordNumber = strSplitLine(eInspectResultsFileColumns.RecordNumber)
                    .DBFilePos = strSplitLine(eInspectResultsFileColumns.DBFilePos)
                    .SpecFilePos = strSplitLine(eInspectResultsFileColumns.SpecFilePos)

                    If strSplitLine.Length >= eInspectResultsFileColumns.PrecursorError + 1 Then
                        ' Inspect version 2008-10-14 added these two Precursor mass columns
                        .PrecursorMZ = strSplitLine(eInspectResultsFileColumns.PrecursorMZ)
                        .PrecursorError = strSplitLine(eInspectResultsFileColumns.PrecursorError)

                        .MH = ComputePeptideMHFromPrecursorInfo(.PrecursorMZ, .PrecursorError, .Charge)


                        Dim dblPrecursorMZ As Double
                        Dim dblPrecursorMonoMass As Double
                        Dim dblPrecursorErrorDa As Double

                        Dim dblPeptideMonoisotopicMass As Double
                        Dim dblPeptideDeltaMassCorrectedPpm As Double

                        If Double.TryParse(.PrecursorMZ, dblPrecursorMZ) Then

                            dblPrecursorMonoMass = clsPeptideMassCalculator.ConvoluteMass(dblPrecursorMZ, .ChargeNum, 0)
                            dblPeptideMonoisotopicMass = .MH - clsPeptideMassCalculator.MASS_PROTON

                            dblPrecursorErrorDa = dblPrecursorMonoMass - dblPeptideMonoisotopicMass

                            dblPeptideDeltaMassCorrectedPpm = ComputeDelMCorrectedPPM(dblPrecursorErrorDa, dblPrecursorMonoMass, _
                                                                                      dblPeptideMonoisotopicMass, True)

							.DelMPPM = NumToString(dblPeptideDeltaMassCorrectedPpm, 4, True)

                        End If


                    Else
                        .PrecursorMZ = "0"
                        .PrecursorError = "0"
                        .MH = 0
                        .DelMPPM = "0"
                    End If

                End With

                blnValidSearchResult = True
            End If

        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing InSpecT Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing InSpecT Results in ParseInspectResultsFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    Private Function ParseInSpectSynFileEntry(ByRef strLineIn As String, _
                                              ByRef intColumnMapping() As Integer, _
                                              ByRef objSearchResult As clsSearchResultsInSpecT, _
                                              ByRef strErrorLog As String, _                                            
                                              ByRef strPeptideSequenceWithMods As String) As Boolean

        ' Parses an entry from the Inspect Synopsis file

		Dim strSplitLine() As String = Nothing

        Dim blnValidSearchResult As Boolean

        Try
            ' Set this to False for now
            blnValidSearchResult = False

            ' Reset objSearchResult
            objSearchResult.Clear()
            strPeptideSequenceWithMods = String.Empty

            strSplitLine = strLineIn.Trim.Split(ControlChars.Tab)
            If strSplitLine.Length < INSPECT_SYN_FILE_MIN_COL_COUNT Then
                Exit Try
            End If

            With objSearchResult
                If Not GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.ResultID), .ResultID) Then
                    Throw New EvaluateException("ResultID column is missing or invalid")
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.Scan), .Scan)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.Charge), .Charge)

                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.Protein), .ProteinName)
                .MultipleProteinCount = "0"

                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.Peptide), strPeptideSequenceWithMods)

                ' Calling this function will set .PeptidePreResidues, .PeptidePostResidues, .PeptideSequenceWithMods, and .PeptideCleanSequence
                .SetPeptideSequenceWithMods(strPeptideSequenceWithMods, True, True)
            End With

            Dim objSearchResultBase As clsSearchResultsBaseClass
            objSearchResultBase = DirectCast(objSearchResult, clsSearchResultsBaseClass)

            MyBase.ComputePseudoPeptideLocInProtein(objSearchResultBase)

            With objSearchResult

                ' Now that the peptide location in the protein has been determined, re-compute the peptide's cleavage and terminus states
                ' If a peptide belongs to several proteins, the cleavage and terminus states shown for the same peptide 
                ' will all be based on the first protein since Inspect only outputs the prefix and suffix letters for the first protein
                .ComputePeptideCleavageStateInProtein()

                If GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.PrecursorError), .PeptideDeltaMass) Then
                    ' Note: .peptideDeltaMass is stored in the Inspect results file as "Observed_Mass - Theoretical_Mass"
                    ' However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                    ' Therefore, we will negate .peptideDeltaMass
                    Try
                        .PeptideDeltaMass = (-Double.Parse(.PeptideDeltaMass)).ToString
                    Catch ex As Exception
                        ' Error; Leave .peptideDeltaMass unchanged
                    End Try

                Else
                    .PeptideDeltaMass = "0"
                End If

                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.MQScore), .MQScore)

                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.Length), .Length)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.TotalPRMScore), .TotalPRMScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.MedianPRMScore), .MedianPRMScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.FractionY), .FractionY)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.FractionB), .FractionB)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.Intensity), .Intensity)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.NTT), .NTT)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.PValue), .pValue)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.FScore), .FScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.DeltaScore), .DeltaScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.DeltaScoreOther), .DeltaScoreOther)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.DeltaNormMQScore), .DeltaNormMQScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.DeltaNormTotalPRMScore), .DeltaNormTotalPRMScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.RankTotalPRMScore), .RankTotalPRMScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.RankFScore), .RankFScore)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.MH), .PeptideMH)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.RecordNumber), .RecordNumber)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.DBFilePos), .DBFilePos)
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.SpecFilePos), .SpecFilePos)

                ' Note: .PrecursorError was processed earlier in this function
                GetColumnValue(strSplitLine, intColumnMapping(eInspectSynFileColumns.PrecursorMZ), .PrecursorMZ)

            End With

            blnValidSearchResult = True
        Catch ex As Exception
            ' Error parsing this row from the synopsis or first hits file
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                If Not strSplitLine Is Nothing AndAlso strSplitLine.Length > 0 Then
                    strErrorLog &= "Error parsing InSpecT Results for RowIndex '" & strSplitLine(0) & "'" & ControlChars.NewLine
                Else
                    strErrorLog &= "Error parsing InSpecT Results in ParseInspectSynFileEntry" & ControlChars.NewLine
                End If
            End If
            blnValidSearchResult = False
        End Try

        Return blnValidSearchResult

    End Function

    ' Main processing function
    Public Overloads Overrides Function ProcessFile(ByVal strInputFilePath As String, ByVal strOutputFolderPath As String, ByVal strParameterFilePath As String) As Boolean
        ' Returns True if success, False if failure

        Dim ioInputFile As System.IO.FileInfo

        Dim strInputFilePathFull As String
        Dim strOutputFilePath As String
        Dim strSynOutputFilePath As String
        Dim strPepToProteinMapFilePath As String

        Dim udtInspectModInfo() As udtModInfoType
        Dim udtPepToProteinMapping() As udtPepToProteinMappingType

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
                        strInputFilePathFull = ioInputFile.FullName

                        ReDim udtInspectModInfo(-1)
                        ReDim udtPepToProteinMapping(-1)

                        ' Load the Inspect Parameter File so that we can determine the modification names and masses
                        If Not ExtractModInfoFromInspectParamFile(mSearchToolParameterFilePath, udtInspectModInfo) Then
                            If udtInspectModInfo Is Nothing OrElse udtInspectModInfo.Length = 0 Then
                                ReDim udtInspectModInfo(0)
                                With udtInspectModInfo(0)
                                    .ModName = PHOS_MOD_NAME.ToLower
                                    .ModMass = PHOS_MOD_MASS
                                    .Residues = PHOS_MOD_RESIDUES
                                    .ModSymbol = UNKNOWN_INSPECT_MOD_SYMBOL
                                End With
                            End If
                        End If


                        ' Resolve the mods in mInspectModInfo with the ModDefs mods
                        ResolveInspectModsWithModDefinitions(udtInspectModInfo)

                        If MyBase.mCreateInspectOrMSGFDbFirstHitsFile Then

                            ' Create the first hits output file
                            MyBase.ResetProgress("Creating the FHT file (top TotalPRMScore)")
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            strOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                            strOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strOutputFilePath & INSPECT_TOTALPRM_FIRST_HITS_FILE_SUFFIX)

                            blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strOutputFilePath, udtInspectModInfo, eFilteredOutputFileTypeConstants.FHTbyTotalPRM)


                            ' Create the first hits output file
                            MyBase.ResetProgress("Creating the FHT file (top FScore)")
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            strOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                            strOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strOutputFilePath & INSPECT_FSCORE_FIRST_HITS_FILE_SUFFIX)

                            blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strOutputFilePath, udtInspectModInfo, eFilteredOutputFileTypeConstants.FHTbyFScore)

                        End If

                        If MyBase.mCreateInspectOrMSGFDbSynopsisFile Then

                            ' Create the synopsis output file
                            MyBase.ResetProgress("Creating the SYN file")
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            'Define the synopsis output file name based on strInputFilePath
                            strSynOutputFilePath = System.IO.Path.GetFileNameWithoutExtension(strInputFilePath)
                            strSynOutputFilePath = System.IO.Path.Combine(strOutputFolderPath, strSynOutputFilePath & SEQUEST_SYNOPSIS_FILE_SUFFIX)

                            blnSuccess = CreateFHTorSYNResultsFile(strInputFilePath, strSynOutputFilePath, udtInspectModInfo, eFilteredOutputFileTypeConstants.SynFile)

                            ' Load the PeptideToProteinMap information; if any modified peptides are present, then write out an updated _inspect_PepToProtMap.txt file with the new mod symbols (file will be named _PepToProtMapMTS.txt)
                            ' If the file doesn't exist, then a warning will be displayed, but processing will continue
                            strPepToProteinMapFilePath = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(strInputFilePathFull), System.IO.Path.GetFileNameWithoutExtension(strInputFilePathFull) & FILENAME_SUFFIX_PEP_TO_PROTEIN_MAPPING & ".txt")

                            MyBase.ResetProgress("Loading the PepToProtein map file: " & System.IO.Path.GetFileName(strPepToProteinMapFilePath))
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            LoadPeptideToProteinMapInfo(strPepToProteinMapFilePath, strOutputFolderPath, udtInspectModInfo, udtPepToProteinMapping)


                            ' Create the other PHRP-specific files
                            MyBase.ResetProgress("Creating the PHRP files for " & System.IO.Path.GetFileName(strSynOutputFilePath))
                            Console.WriteLine()
                            Console.WriteLine()
                            Console.WriteLine(MyBase.ProgressStepDescription)

                            blnSuccess = ParseInSpectSynopsisFile(strSynOutputFilePath, udtPepToProteinMapping, False)

                        End If
                      
                        If blnSuccess Then
                            MyBase.OperationComplete()
                        End If

                    Catch ex As Exception
                        SetErrorMessage("Error calling CreateFHTorSYNResultsFile: " & ex.Message)
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

    Private Function RemoveExtraneousDigits(ByVal strValue As String) As String
        ' If strValue ends in .0000, then remove the .0000 portion
        Static reNumPlusZeroes As System.Text.RegularExpressions.Regex
        Static reAllZeroes As System.Text.RegularExpressions.Regex

        Dim reMatch As System.Text.RegularExpressions.Match

        If reNumPlusZeroes Is Nothing Then
            reNumPlusZeroes = New System.Text.RegularExpressions.Regex("(\.\d*[1-9])0+$", Text.RegularExpressions.RegexOptions.Compiled)
            reAllZeroes = New System.Text.RegularExpressions.Regex("\.0+$", Text.RegularExpressions.RegexOptions.Compiled)
        End If

        If strValue Is Nothing OrElse strValue.Length = 0 Then
            Return String.Empty
        Else
            reMatch = reAllZeroes.Match(strValue)
            If reMatch.Success AndAlso reMatch.Index > 0 Then
                strValue = strValue.Substring(0, reMatch.Index)
            Else
                reMatch = reNumPlusZeroes.Match(strValue)
                If reMatch.Success AndAlso reMatch.Index > 0 Then
                    If reMatch.Groups.Count > 1 Then
                        ' Number is of the form 1.0030 or 1.300 or 1.030
                        strValue = strValue.Substring(0, reMatch.Index) & reMatch.Groups(1).Value
                    End If
                End If
            End If

            Return strValue
        End If

    End Function

    ''' <summary>
    ''' Replaces modification name text in peptide sequences with modification symbols (uses case-sensitive comparisons)
    ''' </summary>
    ''' <param name="strPeptide"></param>
    ''' <param name="udtInspectModInfo">This function assumes that each entry in udtInspectModInfo() has both .ModName and .ModSymbol defined</param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Protected Function ReplaceInspectModTextWithSymbol(ByVal strPeptide As String, ByRef udtInspectModInfo() As udtModInfoType) As String

        Dim intIndex As Integer
        Dim strPrefix As String = String.Empty
        Dim strSuffix As String = String.Empty

        Dim strPeptideNew As String = String.Empty

        Dim intModMass As Integer

        Static reNTerminalModMassRegEx As New System.Text.RegularExpressions.Regex(INSPECT_NTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS)
        Static reCTerminalModMassRegEx As New System.Text.RegularExpressions.Regex(INSPECT_CTERMINAL_MOD_MASS_REGEX, REGEX_OPTIONS)
        Dim reMatch As System.Text.RegularExpressions.Match

        If strPeptide.Length >= 4 Then
            If strPeptide.Chars(1) = "."c AndAlso _
               strPeptide.Chars(strPeptide.Length - 2) = "."c Then
                strPrefix = strPeptide.Substring(0, 2)
                strSuffix = strPeptide.Substring(strPeptide.Length - 2, 2)

                strPeptide = strPeptide.Substring(2, strPeptide.Length - 4)
            End If
        End If

        ' strPeptide should now be the clean peptide, without the prefix or suffix residues
        For intIndex = 0 To udtInspectModInfo.Length - 1
            If udtInspectModInfo(intIndex).ModType <> eInspectModType.StaticMod Then
                strPeptide = strPeptide.Replace(udtInspectModInfo(intIndex).ModName, udtInspectModInfo(intIndex).ModSymbol)

                If udtInspectModInfo(intIndex).ModType = eInspectModType.DynNTermPeptide Or _
                   udtInspectModInfo(intIndex).ModType = eInspectModType.DynCTermPeptide Then

                    If udtInspectModInfo(intIndex).ModType = eInspectModType.DynNTermPeptide Then
                        ' Inspect notates N-terminal mods like this: R.+14HVIFLAER.R   (Note: This behavior is not yet confirmed)
                        ' Look for this using reNTerminalModMassRegEx
                        reMatch = reNTerminalModMassRegEx.Match(strPeptide)

                    ElseIf udtInspectModInfo(intIndex).ModType = eInspectModType.DynCTermPeptide Then
                        ' Inspect notates C-terminal mods like this: R.HVIFLAER+14.R
                        ' Look for this using reCTerminalModMassRegEx
                        reMatch = reCTerminalModMassRegEx.Match(strPeptide)

                    Else
                        ' This code should never be reached
                        reMatch = Nothing
                    End If

                    If Not reMatch Is Nothing Then
                        With reMatch
                            If .Success AndAlso .Groups.Count > 1 Then
                                ' Match found
                                Try
                                    intModMass = CInt(.Groups(1).Value)

                                    ' Compare the mod mass in the specification to this Mod's mod mass
                                    ' If they are less than 0.5 Da apart, then assume we have a match; yes, this assumption is a bit flaky
                                    If Math.Abs(intModMass - CDbl(udtInspectModInfo(intIndex).ModMass)) <= 0.5 Then
                                        ' Match found
                                        ' Replace the matched region with .ModSymbol

                                        If .Groups(0).Index > 0 Then
                                            strPeptideNew = strPeptide.Substring(0, .Groups(0).Index)
                                        Else
                                            strPeptideNew = String.Empty
                                        End If

                                        strPeptideNew &= udtInspectModInfo(intIndex).ModSymbol

                                        If .Groups(0).Index + .Groups(0).Length < strPeptide.Length Then
                                            strPeptideNew &= strPeptide.Substring(.Groups(0).Index + .Groups(0).Length)
                                        End If

                                        strPeptide = String.Copy(strPeptideNew)
                                    End If

                                Catch ex As Exception
                                    Throw New Exception("Error comparing mod mass in peptide to mod mass in udtInspectModInfo", ex)
                                End Try

                            End If
                        End With
                    End If
                End If

            End If
        Next intIndex

        Return strPrefix & strPeptide & strSuffix

    End Function

    Protected Function ReplaceTerminus(ByVal strPeptide As String) As String

        If strPeptide.StartsWith(N_TERMINUS_SYMBOL_INSPECT) Then
            strPeptide = clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST & "." & strPeptide.Substring(N_TERMINUS_SYMBOL_INSPECT.Length)
        End If

        If strPeptide.EndsWith(C_TERMINUS_SYMBOL_INSPECT) Then
            strPeptide = strPeptide.Substring(0, strPeptide.Length - C_TERMINUS_SYMBOL_INSPECT.Length) & "." & clsPeptideCleavageStateCalculator.TERMINUS_SYMBOL_SEQUEST
        End If

        Return strPeptide

    End Function

    Protected Sub ResolveInspectModsWithModDefinitions(ByRef udtInspectModInfo() As udtModInfoType)

        Dim intIndex As Integer
        Dim intResidueIndex As Integer
        Dim intResIndexStart As Integer
        Dim intResIndexEnd As Integer

        Dim dblModMass As Double
        Dim chTargetResidue As Char
        Dim eResidueTerminusState As clsPeptideModificationContainer.eResidueTerminusStateConstants
        Dim blnExistingModFound As Boolean

        Dim objModificationDefinition As clsModificationDefinition

        If Not udtInspectModInfo Is Nothing Then
            For intIndex = 0 To udtInspectModInfo.Length - 1
                ' Call .LookupModificationDefinitionByMass for each entry in udtInspectModInfo

                With udtInspectModInfo(intIndex)
                    If Double.TryParse(.ModMass, dblModMass) Then
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
                            Else
                                chTargetResidue = Nothing
                            End If

                            If .ModType = eInspectModType.DynNTermPeptide Then
                                eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideNTerminus
                            ElseIf .ModType = eInspectModType.DynCTermPeptide Then
                                eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.PeptideCTerminus
                            Else
                                eResidueTerminusState = clsPeptideModificationContainer.eResidueTerminusStateConstants.None
                            End If
                            blnExistingModFound = False

                            objModificationDefinition = mPeptideMods.LookupModificationDefinitionByMass(dblModMass, chTargetResidue, eResidueTerminusState, blnExistingModFound, True)

                            If intResidueIndex = intResIndexStart Then
                                .ModSymbol = objModificationDefinition.ModificationSymbol
                            End If

                        Next intResidueIndex

                    End If
                End With

            Next intIndex
        End If
    End Sub

    Protected Sub StoreOrWriteSearchResult(ByRef swResultFile As System.IO.StreamWriter, _
                                           ByRef intResultID As Integer, _
                                           ByRef udtSearchResult As udtInspectSearchResultType, _
                                           ByRef intFilteredSearchResultCount As Integer, _
                                           ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                           ByRef strErrorLog As String)
        If mSortFHTandSynFiles Then
            If intFilteredSearchResultCount = udtFilteredSearchResults.Length Then
                ReDim Preserve udtFilteredSearchResults(udtFilteredSearchResults.Length * 2 - 1)
            End If

            udtFilteredSearchResults(intFilteredSearchResultCount) = udtSearchResult
            intFilteredSearchResultCount += 1
        Else
            intResultID += 1
            WriteSearchResultToFile(intResultID, swResultFile, udtSearchResult, strErrorLog)
        End If
    End Sub

    Private Sub SortAndWriteFilteredSearchResults(ByRef swResultFile As System.IO.StreamWriter, _
                                                  ByVal intFilteredSearchResultCount As Integer, _
                                                  ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                                  ByRef strErrorLog As String)

        Dim intIndex As Integer

        ' Sort udtFilteredSearchResults by descending TotalPRMScore, ascending scan, ascending charge, ascending peptide, and ascending protein
        Array.Sort(udtFilteredSearchResults, 0, intFilteredSearchResultCount, New InspectSearchResultsComparerTotalPRMDescScanChargePeptide)

        For intIndex = 0 To intFilteredSearchResultCount - 1
            WriteSearchResultToFile(intIndex + 1, swResultFile, udtFilteredSearchResults(intIndex), strErrorLog)
        Next intIndex

    End Sub

    Private Sub StoreTopFHTMatch(ByRef swResultFile As System.IO.StreamWriter, _
                                      ByRef intResultID As Integer, _
                                      ByVal intCurrentScanResultsCount As Integer, _
                                      ByRef udtSearchResultsCurrentScan() As udtInspectSearchResultType, _
                                      ByRef intFilteredSearchResultCount As Integer, _
                                      ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                      ByRef strErrorLog As String, _
                                      ByRef objSortComparer As IComparer)

        Dim intIndex As Integer
        Dim intCurrentCharge As Short = Short.MinValue

        AssignRankAndDeltaNormValues(udtSearchResultsCurrentScan, intCurrentScanResultsCount)

        ' Sort udtFilteredSearchResults by ascending scan, ascending charge, then descending TotalPRMScore or descending FScore (depending on objSortComparer)
        ' All of the data in udtSearchResultsCurrentScan should have the same scan number
        Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortComparer)

        ' Now store or write out the first match for each charge for this scan
        For intIndex = 0 To intCurrentScanResultsCount - 1
            If intIndex = 0 OrElse intCurrentCharge <> udtSearchResultsCurrentScan(intIndex).ChargeNum Then
                StoreOrWriteSearchResult(swResultFile, intResultID, udtSearchResultsCurrentScan(intIndex), intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
                intCurrentCharge = udtSearchResultsCurrentScan(intIndex).ChargeNum
            End If
        Next intIndex

    End Sub

    Private Sub StoreSynMatches(ByRef swResultFile As System.IO.StreamWriter, _
                                ByRef intResultID As Integer, _
                                ByVal intCurrentScanResultsCount As Integer, _
                                ByRef udtSearchResultsCurrentScan() As udtInspectSearchResultType, _
                                ByRef intFilteredSearchResultCount As Integer, _
                                ByRef udtFilteredSearchResults() As udtInspectSearchResultType, _
                                ByRef strErrorLog As String, _
                                ByRef objSortComparer As IComparer)

        Dim intIndex As Integer
        Dim intCurrentCharge As Short = Short.MinValue

        AssignRankAndDeltaNormValues(udtSearchResultsCurrentScan, intCurrentScanResultsCount)

        ' Sort udtFilteredSearchResults by ascending scan, ascending charge, descending TotalPRMScore, and descending FScore
        ' All of the data in udtSearchResultsCurrentScan should have the same scan number
        Array.Sort(udtSearchResultsCurrentScan, 0, intCurrentScanResultsCount, objSortComparer)

        ' Now store or write out the matches that pass the filters
        For intIndex = 0 To intCurrentScanResultsCount - 1
            If udtSearchResultsCurrentScan(intIndex).PValueNum <= mInspectSynopsisFilePValueThreshold OrElse _
               udtSearchResultsCurrentScan(intIndex).TotalPRMScoreNum >= TOTALPRMSCORE_THRESHOLD OrElse _
               udtSearchResultsCurrentScan(intIndex).FScoreNum >= FSCORE_THRESHOLD Then
                StoreOrWriteSearchResult(swResultFile, intResultID, udtSearchResultsCurrentScan(intIndex), intFilteredSearchResultCount, udtFilteredSearchResults, strErrorLog)
            End If
        Next intIndex

    End Sub

    ''' <summary>
    ''' Return the text up to (but not including) the first space in strProteinNameAndDescription
    ''' </summary>
    ''' <param name="strProteinNameAndDescription"></param>
    ''' <returns></returns>
    ''' <remarks></remarks>
    Private Function TruncateProteinName(ByVal strProteinNameAndDescription As String) As String

        Dim intIndex As Integer

        intIndex = strProteinNameAndDescription.IndexOf(" "c)
        If intIndex > 0 Then
            Return strProteinNameAndDescription.Substring(0, intIndex)
        Else
            Return strProteinNameAndDescription
        End If

    End Function

    Private Sub UpdateSearchResultEnzymeAndTerminusInfo(ByRef objSearchResult As clsSearchResultsInSpecT)
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

    Private Sub WriteSynFHTFileHeader(ByRef swResultFile As System.IO.StreamWriter, _
                                      ByRef strErrorLog As String)

        ' Write out the header line for synopsis / first hits files
        Try
            swResultFile.WriteLine("ResultID" & ControlChars.Tab & _
                                   "Scan" & ControlChars.Tab & _
                                   "Peptide" & ControlChars.Tab & _
                                   "Protein" & ControlChars.Tab & _
                                   "Charge" & ControlChars.Tab & _
                                   "MQScore" & ControlChars.Tab & _
                                   "Length" & ControlChars.Tab & _
                                   "TotalPRMScore" & ControlChars.Tab & _
                                   "MedianPRMScore" & ControlChars.Tab & _
                                   "FractionY" & ControlChars.Tab & _
                                   "FractionB" & ControlChars.Tab & _
                                   "Intensity" & ControlChars.Tab & _
                                   "NTT" & ControlChars.Tab & _
                                   "PValue" & ControlChars.Tab & _
                                   "FScore" & ControlChars.Tab & _
                                   "DeltaScore" & ControlChars.Tab & _
                                   "DeltaScoreOther" & ControlChars.Tab & _
                                   "DeltaNormMQScore" & ControlChars.Tab & _
                                   "DeltaNormTotalPRMScore" & ControlChars.Tab & _
                                   "RankTotalPRMScore" & ControlChars.Tab & _
                                   "RankFScore" & ControlChars.Tab & _
                                   "MH" & ControlChars.Tab & _
                                   "RecordNumber" & ControlChars.Tab & _
                                   "DBFilePos" & ControlChars.Tab & _
                                   "SpecFilePos" & ControlChars.Tab & _
                                   "PrecursorMZ" & ControlChars.Tab & _
                                   "PrecursorError" & ControlChars.Tab & _
                                   "DelM_PPM")
        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits header" & ControlChars.NewLine
            End If
        End Try

    End Sub

    Private Sub WriteSearchResultToFile(ByVal intResultID As Integer, _
                                        ByRef swResultFile As System.IO.StreamWriter, _
                                        ByRef udtSearchResult As udtInspectSearchResultType, _
                                        ByRef strErrorLog As String)

        ' Writes an entry to a synopsis or first hits file
        Try
			swResultFile.WriteLine(intResultID.ToString & ControlChars.Tab & _
			  udtSearchResult.Scan & ControlChars.Tab & _
			  udtSearchResult.PeptideAnnotation & ControlChars.Tab & _
			  udtSearchResult.Protein & ControlChars.Tab & _
			  udtSearchResult.Charge & ControlChars.Tab & _
			  udtSearchResult.MQScore & ControlChars.Tab & _
			  udtSearchResult.Length & ControlChars.Tab & _
			  udtSearchResult.TotalPRMScore & ControlChars.Tab & _
			  udtSearchResult.MedianPRMScore & ControlChars.Tab & _
			  udtSearchResult.FractionY & ControlChars.Tab & _
			  udtSearchResult.FractionB & ControlChars.Tab & _
			  udtSearchResult.Intensity & ControlChars.Tab & _
			  udtSearchResult.NTT & ControlChars.Tab & _
			  udtSearchResult.pValue & ControlChars.Tab & _
			  udtSearchResult.FScore & ControlChars.Tab & _
			  udtSearchResult.DeltaScore & ControlChars.Tab & _
			  udtSearchResult.DeltaScoreOther & ControlChars.Tab & _
			  NumToString(udtSearchResult.DeltaNormMQScore, 5, True) & ControlChars.Tab & _
			  NumToString(udtSearchResult.DeltaNormTotalPRMScore, 5, True) & ControlChars.Tab & _
			  udtSearchResult.RankTotalPRMScore & ControlChars.Tab & _
			  udtSearchResult.RankFScore & ControlChars.Tab & _
			  NumToString(udtSearchResult.MH, 6, True) & ControlChars.Tab & _
			  udtSearchResult.RecordNumber & ControlChars.Tab & _
			  udtSearchResult.DBFilePos & ControlChars.Tab & _
			  udtSearchResult.SpecFilePos & ControlChars.Tab & _
			  udtSearchResult.PrecursorMZ & ControlChars.Tab & _
			  udtSearchResult.PrecursorError & ControlChars.Tab & _
			  udtSearchResult.DelMPPM)

        Catch ex As Exception
            If strErrorLog.Length < MAX_ERROR_LOG_LENGTH Then
                strErrorLog &= "Error writing synopsis / first hits record" & ControlChars.NewLine
            End If
        End Try

    End Sub

#Region "IComparer Classes"

    Protected Class InspectSearchResultsComparerTotalPRMDescScanChargePeptide
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtInspectSearchResultType = DirectCast(x, udtInspectSearchResultType)
            Dim yData As udtInspectSearchResultType = DirectCast(y, udtInspectSearchResultType)

            If xData.TotalPRMScoreNum > yData.TotalPRMScoreNum Then
                Return -1
            ElseIf xData.TotalPRMScoreNum < yData.TotalPRMScoreNum Then
                Return 1
            Else
                ' TotalPRMScore is the same; check scan number
                If xData.ScanNum > yData.ScanNum Then
                    Return 1
                ElseIf xData.ScanNum < yData.ScanNum Then
                    Return -1
                Else
                    ' Scan is the same, check charge
                    If xData.ChargeNum > yData.ChargeNum Then
                        Return 1
                    ElseIf xData.ChargeNum < yData.ChargeNum Then
                        Return -1
                    Else
                        ' Charge is the same; check peptide
                        If xData.PeptideAnnotation > yData.PeptideAnnotation Then
                            Return 1
                        ElseIf xData.PeptideAnnotation < yData.PeptideAnnotation Then
                            Return -1
                        Else
                            ' Peptide is the same, check Protein
                            If xData.Protein > yData.Protein Then
                                Return 1
                            ElseIf xData.Protein < yData.Protein Then
                                Return -1
                            Else
                                Return 0
                            End If
                        End If
                    End If
                End If
            End If

        End Function
    End Class

    Protected Class InspectSearchResultsComparerScanChargeTotalPRMDescFScoreDesc
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtInspectSearchResultType = DirectCast(x, udtInspectSearchResultType)
            Dim yData As udtInspectSearchResultType = DirectCast(y, udtInspectSearchResultType)

            If xData.ScanNum > yData.ScanNum Then
                Return 1
            ElseIf xData.ScanNum < yData.ScanNum Then
                Return -1
            Else
                If xData.ChargeNum > yData.ChargeNum Then
                    Return 1
                ElseIf xData.ChargeNum < yData.ChargeNum Then
                    Return -1
                Else
                    ' Charge is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                    If xData.TotalPRMScoreNum > yData.TotalPRMScoreNum Then
                        Return -1
                    ElseIf xData.TotalPRMScoreNum < yData.TotalPRMScoreNum Then
                        Return 1
                    Else
                        ' TotalPRMScore is the same; check FScore (sort on descending FScore)
                        If xData.FScoreNum > yData.FScoreNum Then
                            Return -1
                        ElseIf xData.FScoreNum < yData.FScoreNum Then
                            Return 1
                        Else
                            Return 0
                        End If
                    End If
                End If
            End If

        End Function
    End Class

    Protected Class InspectSearchResultsComparerScanChargeFScoreDescTotalPRMDesc
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtInspectSearchResultType = DirectCast(x, udtInspectSearchResultType)
            Dim yData As udtInspectSearchResultType = DirectCast(y, udtInspectSearchResultType)

            If xData.ScanNum > yData.ScanNum Then
                Return 1
            ElseIf xData.ScanNum < yData.ScanNum Then
                Return -1
            Else
                If xData.ChargeNum > yData.ChargeNum Then
                    Return 1
                ElseIf xData.ChargeNum < yData.ChargeNum Then
                    Return -1
                Else
                    ' Charge is the same; check FScore (sort on descending FScore)
                    If xData.FScoreNum > yData.FScoreNum Then
                        Return -1
                    ElseIf xData.FScoreNum < yData.FScoreNum Then
                        Return 1
                    Else
                        ' FScore is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                        If xData.TotalPRMScoreNum > yData.TotalPRMScoreNum Then
                            Return -1
                        ElseIf xData.TotalPRMScoreNum < yData.TotalPRMScoreNum Then
                            Return 1
                        Else
                            Return 0
                        End If
                    End If
                End If
            End If

        End Function
    End Class

    Protected Class InspectSearchResultsComparerScanChargeMQScoreDescTotalPRMDesc
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtInspectSearchResultType = DirectCast(x, udtInspectSearchResultType)
            Dim yData As udtInspectSearchResultType = DirectCast(y, udtInspectSearchResultType)

            If xData.ScanNum > yData.ScanNum Then
                Return 1
            ElseIf xData.ScanNum < yData.ScanNum Then
                Return -1
            Else
                If xData.ChargeNum > yData.ChargeNum Then
                    Return 1
                ElseIf xData.ChargeNum < yData.ChargeNum Then
                    Return -1
                Else
                    ' Charge is the same; check MQScore (sort on descending MQScore)
                    If xData.MQScoreNum > yData.MQScoreNum Then
                        Return -1
                    ElseIf xData.MQScoreNum < yData.MQScoreNum Then
                        Return 1
                    Else
                        ' MQScore is the same; check TotalPRMScore (sort on descending TotalPRMScore)
                        If xData.TotalPRMScoreNum > yData.TotalPRMScoreNum Then
                            Return -1
                        ElseIf xData.TotalPRMScoreNum < yData.TotalPRMScoreNum Then
                            Return 1
                        Else
                            Return 0
                        End If
                    End If
                End If
            End If

        End Function
    End Class

    Protected Class PepToProteinMappingComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtPepToProteinMappingType = DirectCast(x, udtPepToProteinMappingType)
            Dim yData As udtPepToProteinMappingType = DirectCast(y, udtPepToProteinMappingType)

            If xData.Peptide > yData.Peptide Then
                Return 1
            ElseIf xData.Peptide < yData.Peptide Then
                Return -1
            Else
                If xData.Protein > yData.Protein Then
                    Return 1
                ElseIf xData.Protein < yData.Protein Then
                    Return -1
                Else
                    Return 0
                End If
            End If

        End Function
    End Class

    Protected Class PepToProteinMappingPeptideSearchComparer
        Implements System.Collections.IComparer

        Public Function Compare(ByVal x As Object, ByVal y As Object) As Integer Implements System.Collections.IComparer.Compare
            Dim xData As udtPepToProteinMappingType = DirectCast(x, udtPepToProteinMappingType)
            Dim strPeptide As String = DirectCast(y, String)

            If xData.Peptide > strPeptide Then
                Return 1
            ElseIf xData.Peptide < strPeptide Then
                Return -1
            Else

                Return 0
            End If

        End Function
    End Class

#End Region

End Class
