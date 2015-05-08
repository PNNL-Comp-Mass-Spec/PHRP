Option Strict On

' This class reads in a MSPathFinder results file (_IcTda.tsv) and creates 
' a tab-delimited text file with the data. 
'
' -------------------------------------------------------------------------------
' Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
' Started 05/01/2015
'
' E-mail: matthew.monroe@pnnl.gov
' -------------------------------------------------------------------------------

Imports PHRPReader
Imports System.IO
Imports System.Text.RegularExpressions

Public Class clsMSPathFinderResultsProcessor
    Inherits clsPHRPBaseClass

    Public Sub New()
        MyBase.New()
        MyBase.mFileDate = "May 1, 2015"
    End Sub

#Region "Constants and Enums"

    Public Const FILENAME_SUFFIX_MODA_FILE As String = "_IcTda.tsv"

    Public Const N_TERMINUS_SYMBOL_MSPATHFINDER As String = "-"
    Public Const C_TERMINUS_SYMBOL_MSPATHFINDER As String = "-"

    Public Const DEFAULT_SYN_FILE_QVALUE_THRESHOLD As Single = 0.05

    Private Const MAX_ERROR_LOG_LENGTH As Integer = 4096

    Private Const MODA_MASS_DIGITS_OF_PRECISION As Byte = 0

    Private Const REGEX_OPTIONS As RegexOptions = RegexOptions.Compiled Or RegexOptions.Singleline Or RegexOptions.IgnoreCase

    ' These columns correspond to the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
    Protected Const MSPathFinderResultsFileColCount As Integer = 17
    Public Enum eMSPathFinderResultsFileColumns As Integer
        Scan = 0
        PrefixResidue = 1
        Sequence = 2
        SuffixResidue = 3
        Modifications = 4
        Composition = 5
        ProteinName = 6
        ProteinDesc = 7
        ProteinLength = 8
        ResidueStart = 9
        ResidueEnd = 10
        Charge = 11
        MostAbundantIsotopeMz = 12
        Mass = 13
        NumMatchedFragments = 14
        QValue = 15
        PepQValue = 16
    End Enum

    ' These columns correspond to the Synopsis file created by this class
    Protected Const MSPathFinderSynFileColCount As Integer = 14
    Public Enum eMSPathFinderSynFileColumns As Integer
        ResultID = 0
        Scan = 1
        Charge = 2
        MostAbundantIsotopeMz = 3
        Mass = 4
        Pre = 5
        Sequence = 6
        Post = 7
        Modifications = 8
        Composition = 9
        ProteinName = 10
        ProteinDesc = 11
        ProteinLength = 12
        ResidueStart = 13
        ResidueEnd = 14
        MatchedFragments = 15
        QValue = 16
        PepQValue = 17
    End Enum

#End Region

#Region "Structures"
    ' This data structure holds rows read from the tab-delimited file (_IcTda.tsv) created directly by MSPathFinder
    Protected Structure udtMSPathFinderSearchResultType

        Public Scan As String
        Public ScanNum As Integer
        Public PrefixResidue As String
        Public Sequence As String
        Public SuffixResidue As String
        Public Modifications As String
        Public Composition As String
        Public ProteinName As String
        Public ProteinDesc As String
        Public ProteinLength As String
        Public ResidueStart As String
        Public ResidueEnd As String
        Public Charge As String
        Public ChargeNum As Short
        Public MostAbundantIsotopeMz As String      ' As reported by MSPathfinder
        Public Mass As String                       ' Theoretical monoisotopic mass of the peptide (uncharged, including mods), as computed by MSPathFinder
        Public MH As String                         ' Theoretical monoisotopic peptide MH (including mods), as computed by PHRP; note that this is (M+H)+
        Public NumMatchedFragments As String
        Public QValue As String                     ' FDR, at the scan level
        Public PepQValue As String                  ' FDR, at the peptide level

        ' Unused at present: Public DelM As String                       ' Computed using Precursor_mass - CalculatedMonoMass
        ' Unused at present: Public DelM_PPM As String                   ' Computed using DelM and CalculatedMonoMass	

        Public Sub Clear()
            Scan = String.Empty
            ScanNum = 0
            PrefixResidue = String.Empty
            Sequence = String.Empty
            SuffixResidue = String.Empty
            Modifications = String.Empty
            Composition = String.Empty
            ProteinName = String.Empty
            ProteinDesc = String.Empty
            ProteinLength = String.Empty
            ResidueStart = String.Empty
            ResidueEnd = String.Empty
            Charge = String.Empty
            ChargeNum = 0
            MostAbundantIsotopeMz = String.Empty
            Mass = String.Empty
            MH = String.Empty
            NumMatchedFragments = String.Empty
            QValue = String.Empty
            PepQValue = String.Empty

            ' Unused at present: DelM = String.Empty
            ' Unused at present: DelM_PPM = String.Empty

        End Sub
    End Structure

#End Region

#Region "Classwide Variables"

#End Region

    Public Overrides Function ProcessFile(strInputFilePath As String, strOutputFolderPath As String, strParameterFilePath As String) As Boolean
        Throw New NotImplementedException
    End Function

End Class
