// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or https://www.pnnl.gov/sysbio/ or https://panomics.pnnl.gov/
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class reads in an X!Tandem results file (XML format) and creates
    /// a tab-delimited text file with the data.  It will insert modification symbols
    /// into the peptide sequences for modified peptides.  The user can optionally provide
    /// a modification definition file which specifies the symbol to use for each
    /// modification mass.  If the user does not provide this file, then the modification
    /// definition information is determined from the X!Tandem Input Parameters section at
    /// the end of the X!Tandem results file.
    /// </summary>
    public class clsXTandemResultsProcessor : clsPHRPBaseClass
    {

        #region "Constants and Enums"

        // Note: These names must all be lowercase
        private const string XTANDEM_XML_ROOT_ELEMENT = "bioml";
        private const string XTANDEM_XML_ELEMENT_NAME_GROUP = "group";
        private const string XTANDEM_XML_ELEMENT_NAME_PROTEIN = "protein";
        private const string XTANDEM_XML_ELEMENT_NAME_PEPTIDE = "peptide";
        private const string XTANDEM_XML_ELEMENT_NAME_DOMAIN = "domain";
        private const string XTANDEM_XML_ELEMENT_NAME_AMINO_ACID = "aa";
        private const string XTANDEM_XML_ELEMENT_NAME_NOTE = "note";

        private const string XTANDEM_XML_GROUP_TYPE_MODEL = "model";
        private const string XTANDEM_XML_GROUP_TYPE_SUPPORT = "support";
        private const string XTANDEM_XML_GROUP_TYPE_PARAMETERS = "parameters";

        private const string XML_ERROR_ROOT_LEVEL_INVALID = "The data at the root level is invalid";
        private const int MAX_ERROR_LOG_LENGTH = 4096;

        private const string SCAN_NUMBER_EXTRACTION_REGEX_A = @"scan=(\d+)";
        private const string SCAN_NUMBER_EXTRACTION_REGEX_B = @"scan\s*(\d+)";
        private const string SCAN_NUMBER_EXTRACTION_REGEX_C = @"(\d+)\.\d+\.\d\.dta";
        private const string SCAN_NUMBER_EXTRACTION_REGEX_D = @"\d+";

        private const string REVERSED_PROTEIN_SEQUENCE_INDICATOR = ":reversed";
        private const string PROTEIN_DESCRIPTION_LABEL = "description";

        private enum eCurrentXMLDataFileSectionConstants
        {
            UnknownFile = 0,
            Start = 1,
            SearchResults = 2,
            InputParameters = 3,
            PerformanceParameters = 4
        }

        const int INPUT_PARAM_LABEL_NAMES_MAX_INDEX = 13;
        private enum eInputParamLabelNames
        {
            Residue_StaticModMass = 0,
            Residue_PotentialModMass = 1,
            Residue_PotentialModMotif = 2,
            Refine_PotentialModMass = 3,
            Refine_PotentialModMotif = 4,
            Refine_PotentialNTerminusMods = 5,
            Refine_PotentialCTerminusMods = 6,
            Protein_NTerminal_ResidueModMass = 7,
            Protein_CTerminal_ResidueModMass = 8,
            Protein_Cleavage_NTerminalMassChange = 9,
            Protein_Cleavage_CTerminalMassChange = 10,
            Protein_Cleavage_Site = 11,
            Refine_ModificationMass = 12,
            Scoring_Include_Reverse = 13
        }
        #endregion

        #region "Structures"
        #endregion

        #region "Classwide Variables"

        private int mNextResultID;
        private bool mLookForReverseSequenceTag;

        private Regex mScanNumberRegExA;
        private Regex mScanNumberRegExB;
        private Regex mScanNumberRegExC;
        private Regex mScanNumberRegExD;

        private readonly Dictionary<string, int> mSeqsWithMods;
        private readonly SortedSet<string> mSeqsWithoutMods;

        #endregion

        #region "Properties"
        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        public clsXTandemResultsProcessor()
        {
            mFileDate = "March 28, 2019";

            mSeqsWithMods = new Dictionary<string, int>();
            mSeqsWithoutMods = new SortedSet<string>();

            InitializeLocalVariables();
        }

        private bool AddModificationsAndComputeMass(clsSearchResultsXTandem searchResult, bool updateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = false;
            bool success;

            try
            {
                // If any modifications of type IsotopicMod are defined we would add them to the Search Result Mods now
                // However, since X!Tandem doesn't support Isotopic Mods, this step is currently skipped
                //
                // searchResult.SearchResultAddIsotopicModifications(updateModOccurrenceCounts)

                // Add the protein terminus static mods (if defined and if the peptide is at a protein terminus)
                // Function .SearchResultAddStaticTerminusMods() will only add the terminus mod if the terminus
                //  is not already modified by the given terminus mod mass
                // This function will also add peptide terminus static mods, if defined, though those are not supported by X!Tandem and therefore should not be defined
                searchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, updateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                searchResult.ComputeMonoisotopicMass();

                // Update PeptideDeltaMassCorrectedPpm
                searchResult.ComputeDelMCorrectedXT();

                // Populate .PeptideSequenceWithMods and .PeptideModDescription
                // Note that this function will call .AddSearchResultModificationsToCleanSequence() then .UpdateModDescription()
                searchResult.ApplyModificationInformation();

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
        }

        private string ConvertEValueToBase10Log(string expectationValue)
        {
            double eValue;
            double logEValue;

            try
            {
                eValue = double.Parse(expectationValue);
            }
            catch (Exception)
            {
                eValue = 0;
            }

            try
            {
                if (eValue <= 0)
                {
                    logEValue = 0;
                }
                else if (eValue <= 1E-307)
                {
                    logEValue = -307;
                }
                else
                {
                    logEValue = Math.Round(Math.Log10(eValue), 3);
                }
            }
            catch (Exception)
            {
                logEValue = 0;
            }

            return logEValue.ToString("0.000");
        }

        private void InitializeLocalVariables()
        {
            // Note: This function is called from ParseXTandemResultsFile()
            // These variables will therefore be reset for each XTandem XML file analyzed
            mNextResultID = 1;
            mLookForReverseSequenceTag = false;

            InitializeRegExObjects();
        }

        private void InitializeRegExObjects()
        {
            const RegexOptions REGEX_OPTIONS = RegexOptions.Compiled | RegexOptions.Singleline | RegexOptions.IgnoreCase;

            try
            {
                mScanNumberRegExA = new Regex(SCAN_NUMBER_EXTRACTION_REGEX_A, REGEX_OPTIONS);
                mScanNumberRegExB = new Regex(SCAN_NUMBER_EXTRACTION_REGEX_B, REGEX_OPTIONS);
                mScanNumberRegExC = new Regex(SCAN_NUMBER_EXTRACTION_REGEX_C, REGEX_OPTIONS);
                mScanNumberRegExD = new Regex(SCAN_NUMBER_EXTRACTION_REGEX_D, REGEX_OPTIONS);
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        private bool ParseXTandemInputParameterModInfo(
            clsModificationDefinition.eModificationTypeConstants eModificationType,
            int sortOrder,
            bool parsingMotifDef,
            string paramValue,
            ref int modInfoCount,
            ref udtSearchOptionModificationInfoType[] udtModInfo)
        {
            // Parse out the mod information defined in paramValue
            // Add each entry to udtModInfo
            // If parsingMotifDef = True, then do not try to determine .TargetResidues

            const char MOD_LIST_SEP_CHAR = ',';

            bool success;

            try
            {
                if (string.IsNullOrWhiteSpace(paramValue))
                {
                    // Empty parameter; no definition to parse
                }
                else
                {
                    // Parse paramValue
                    // Each modification entry can have multiple modification definitions separated separated by commas, so we first split paramValue
                    var modDefs = paramValue.Split(MOD_LIST_SEP_CHAR);
                    for (var index = 0; index <= modDefs.Length - 1; index++)
                    {
                        // Modification definitions typically look like "15.9949@M"
                        // However, a neutral loss can be specified using "79.9663:-97.98@STY"
                        //   Thus, mod mass is number up to first non-numeric character (typically a colon or @ sign)
                        // Target residues are the residues after the @
                        // If the Target residues contain an X, then the modification can apply to any residue

                        var modificationMass = 0.0;

                        // Look for a colon and an @ sign
                        var colonIndex = modDefs[index].IndexOf(':');
                        var atSignIndex = modDefs[index].IndexOf('@');

                        if (atSignIndex < 0)
                        {
                            // At sign not found; skip this mod def
                        }
                        else
                        {
                            if (colonIndex > atSignIndex)
                            {
                                // Ignore this colon since it's present after the @ sign
                                colonIndex = -1;
                            }

                            if (colonIndex > 0)
                            {
                                // Colon found; see if the text up to colonIndex is a number
                                if (clsPHRPParser.IsNumber(modDefs[index].Substring(0, colonIndex)))
                                {
                                    modificationMass = double.Parse(modDefs[index].Substring(0, colonIndex));
                                }
                            }
                            else
                            {
                                // Colon not found; see if the text up to atSignIndex is a number
                                if (clsPHRPParser.IsNumber(modDefs[index].Substring(0, atSignIndex)))
                                {
                                    modificationMass = double.Parse(modDefs[index].Substring(0, atSignIndex));
                                }
                            }

                            if (Math.Abs(modificationMass) > float.Epsilon)
                            {
                                string targetResidues;

                                // Valid mass found; now extract the target residues
                                if (!parsingMotifDef && atSignIndex + 1 < modDefs[index].Length)
                                {
                                    targetResidues = modDefs[index].Substring(atSignIndex + 1);
                                }
                                else
                                {
                                    targetResidues = string.Empty;
                                }

                                if (targetResidues.IndexOf('X') >= 0)
                                {
                                    // Modification can affect any residue; set targetResidues to ""
                                    targetResidues = string.Empty;
                                }

                                if (targetResidues.Length > 0)
                                {
                                    // Convert from X!Tandem-style N-Terminus notation to DMS-style notation
                                    targetResidues = targetResidues.Replace(clsPeptideModificationContainer.N_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM, clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
                                    targetResidues = targetResidues.Replace(clsPeptideModificationContainer.C_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM, clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);
                                }

                                // Append the new mod information to udtModInfo
                                if (modInfoCount >= udtModInfo.Length)
                                {
                                    Array.Resize(ref udtModInfo, udtModInfo.Length * 2);
                                }

                                var modInfo = udtModInfo[modInfoCount];
                                modInfo.SortOrder = sortOrder;
                                modInfo.ModificationMass = modificationMass;
                                modInfo.TargetResidues = targetResidues;
                                modInfo.ModificationType = eModificationType;
                                udtModInfo[modInfoCount] = modInfo;
                                modInfoCount += 1;
                            }
                        }
                    }
                }

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
        }

        private bool ParseXTandemInputParameterProteinTerminusMod(
            int sortOrder,
            bool nTerminus,
            string paramValue,
            ref int modInfoCount,
            ref udtSearchOptionModificationInfoType[] udtModInfo)
        {
            // Parse out the mass defined in paramValue
            // Add the entry to udtModInfo if non-zero
            // if nTerminus = True then mod applies to the protein's N-Terminus; otherwise, applies to the protein's C-Terminus

            bool success;

            try
            {
                if (string.IsNullOrWhiteSpace(paramValue))
                {
                    // Empty parameter; no definition to parse
                }
                else
                {
                    // See if paramValue is a non-zero number

                    var modificationMass = 0.0;
                    if (clsPHRPParser.IsNumber(paramValue))
                    {
                        modificationMass = double.Parse(paramValue);
                    }

                    if (Math.Abs(modificationMass) > float.Epsilon)
                    {
                        // Append the new mod information to udtModInfo
                        if (modInfoCount >= udtModInfo.Length)
                        {
                            Array.Resize(ref udtModInfo, udtModInfo.Length * 2);
                        }

                        var modInfo = udtModInfo[modInfoCount];
                        modInfo.SortOrder = sortOrder;
                        string targetResidues;
                        if (nTerminus)
                        {
                            targetResidues = clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                        }
                        else
                        {
                            targetResidues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                        }

                        modInfo.ModificationMass = modificationMass;
                        modInfo.TargetResidues = targetResidues;
                        modInfo.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod;
                        udtModInfo[modInfoCount] = modInfo;
                        modInfoCount += 1;
                    }
                }

                success = true;
            }
            catch (Exception)
            {
                success = false;
            }

            return success;
        }

        private bool ParseXTandemResultsFile(string inputFilePath, string outputFilePath, bool resetMassCorrectionTagsAndModificationDefinitions = true)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile

            var resultsProcessed = 0;

            bool success;

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (resetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                // Initialize the searchResults[] array; initially reserve space for 4 proteins
                // Note: The number of valid entries in searchResults[] is given by searchResultCount; searchResults[] is expanded but never shrunk
                // There is a separate entry in searchResults[] for each protein encountered

                var searchResultCount = 0;
                var searchResults = new clsSearchResultsXTandem[4];

                int searchResultIndex;
                for (searchResultIndex = 0; searchResultIndex <= searchResults.Length - 1; searchResultIndex++)
                {
                    searchResults[searchResultIndex] = new clsSearchResultsXTandem(mPeptideMods, mPeptideSeqMassCalculator);
                }

                // Reset mNextResultID and mLookForReverseSequenceTag
                InitializeLocalVariables();

                try
                {
                    // Read the input parameters from the end of the X!Tandem results file (inputFilePath)
                    success = ParseXTandemResultsFileInputParameters(inputFilePath);
                    if (!success)
                    {
                        SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile, true);
                        return false;
                    }

                    for (searchResultIndex = 0; searchResultIndex <= searchResults.Length - 1; searchResultIndex++)
                    {
                        searchResults[searchResultIndex].UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);
                    }

                    var errorLog = string.Empty;

                    // Open the input file and parse it

                    // Initialize the stream reader and the XML Text Reader
                    eCurrentXMLDataFileSectionConstants eCurrentXMLDataFileSection;
                    using (var reader = new StreamReader(inputFilePath))
                    {
                        using (var xmlReader = new XmlTextReader(reader))
                        {
                            resultsProcessed = 0;

                            // Create the output file
                            using (var writer = new StreamWriter(outputFilePath, false))
                            {
                                // Write the header line
                                WriteSynFHTFileHeader(writer, ref errorLog);

                                // Create the additional output files
                                success = InitializeSequenceOutputFiles(outputFilePath);

                                // Parse the input file
                                eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile;

                                while (xmlReader.Read() & !AbortProcessing)
                                {
                                    XMLTextReaderSkipWhitespace(xmlReader);
                                    if (xmlReader.ReadState != ReadState.Interactive)
                                        break;

                                    if (xmlReader.Depth < 2)
                                    {
                                        if (xmlReader.NodeType == XmlNodeType.Element)
                                        {
                                            switch (xmlReader.Name.ToLower())
                                            {
                                                case XTANDEM_XML_ELEMENT_NAME_GROUP:
                                                    if (xmlReader.HasAttributes)
                                                    {
                                                        // Cache the XML reader depth before reading any of the element's attributes
                                                        var groupElementReaderDepth = xmlReader.Depth;

                                                        // See if the group has a "type" attribute containing the text XTANDEM_XML_GROUP_TYPE_MODEL
                                                        var currentGroupType = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty);
                                                        if (currentGroupType == XTANDEM_XML_GROUP_TYPE_MODEL)
                                                        {
                                                            eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.SearchResults;

                                                            ParseXTandemResultsFileEntry(xmlReader, writer, ref searchResultCount, ref searchResults, ref errorLog, groupElementReaderDepth);
                                                            resultsProcessed += 1;

                                                            // Update the progress
                                                            var percentComplete = Convert.ToSingle(reader.BaseStream.Position / reader.BaseStream.Length * 100);
                                                            if (CreateProteinModsFile)
                                                            {
                                                                percentComplete = percentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                                                            }
                                                            UpdateProgress(percentComplete);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        // Group doesn't have any attributes; ignore it
                                                        xmlReader.Skip();
                                                    }

                                                    break;
                                                case XTANDEM_XML_ROOT_ELEMENT:
                                                    eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.Start;
                                                    break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (eCurrentXMLDataFileSection == eCurrentXMLDataFileSectionConstants.UnknownFile)
                    {
                        mErrorMessage = "Root element '" + XTANDEM_XML_ROOT_ELEMENT + "' not found in the input file: " + "\n" + inputFilePath;
                    }
                    else
                    {
                        if (CreateModificationSummaryFile)
                        {
                            // Create the modification summary file
                            var inputFile = new FileInfo(inputFilePath);
                            var outputFile = new FileInfo(outputFilePath);

                            if (outputFile.Directory == null)
                            {
                                ReportWarning("ParseXTandemResultsFile: Could not determine the parent directory of the output file, " + outputFile.FullName);
                            }
                            else
                            {
                                var modificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));

                                if (string.IsNullOrWhiteSpace(modificationSummaryFilePath))
                                {
                                    ReportWarning("ParseXTandemResultsFile: modificationSummaryFilePath is empty; cannot call SaveModificationSummaryFile");
                                }
                                else
                                {
                                    SaveModificationSummaryFile(Path.Combine(outputFile.Directory.FullName, modificationSummaryFilePath));
                                }
                            }
                        }

                        // Inform the user if any errors occurred
                        if (errorLog.Length > 0)
                        {
                            SetErrorMessage("Invalid Lines: " + "\n" + errorLog);
                        }
                    }

                    success = true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    success = false;
                }
                finally
                {
                    CloseSequenceOutputFiles();
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorCreatingOutputFiles);
                success = false;
            }

            return success;
        }

        private bool ParseXTandemResultsFileEntry(XmlReader xmlReader, StreamWriter writer, ref int searchResultCount, ref clsSearchResultsXTandem[] searchResults, ref string errorLog, int groupElementReaderDepth)
        {
            // Note: The number of valid entries in searchResults[) is given by searchResultCount; searchResults(] is expanded but never shrunk
            // There is a separate entry in searchResults() for each protein encountered

            const string GROUP_LABEL_PROTEIN = "protein";
            const string GROUP_LABEL_FRAG_ION = "fragment ion mass spectrum";

            // The following is the XSLT originally used to parse out the data
            // ID: <xsl:value-of select="@id" />
            // Charge: <xsl:value-of select="@z" />
            // ParentIonMH: <xsl:value-of select="@mh" />
            // Peptide_Expectation_Value_e: <xsl:value-of select="@expect" />
            // Protein_Name: <xsl:value-of select="substring-before(concat(@label,' '), ' ')" />
            // Peptide_Intensity_Log(I): <xsl:value-of select="@sumI" />

            // Protein_Expectation_Value_Log(e): <xsl:value-of select="./protein/@expect" />
            // Protein_Intensity_Log(I): <xsl:value-of select="./protein/@sumI" />

            // Peptide_Hyperscore: <xsl:value-of select="./protein/peptide/domain/@hyperscore" />
            // Peptide_Sequence: <xsl:if test = "substring(./protein/peptide/domain/@pre,string-length(./protein/peptide/domain/@pre),1) = '['" ><xsl:text>-</xsl:text></xsl:if><xsl:if test = "substring(./protein/peptide/domain/@pre,string-length(./protein/peptide/domain/@pre),1) != '['" ><xsl:value-of select="substring(./protein/peptide/domain/@pre,string-length(./protein/peptide/domain/@pre),1)" /></xsl:if>.<xsl:value-of select="./protein/peptide/domain/@seq" />.<xsl:if test = "substring(./protein/peptide/domain/@post,1,1) = ']'" ><xsl:text>-</xsl:text></xsl:if><xsl:if test = "substring(./protein/peptide/domain/@post,1,1) != ']'" ><xsl:value-of select="substring(./protein/peptide/domain/@post,1,1)" /></xsl:if>
            // DeltaCn2: <xsl:value-of select="round((./protein/peptide/domain/@hyperscore - ./protein/peptide/domain/@nextscore) div ./protein/peptide/domain/@hyperscore * 10000) div 10000" />
            // y_score: <xsl:value-of select="./protein/peptide/domain/@y_score" />
            // y_ions: <xsl:value-of select="./protein/peptide/domain/@y_ions" />
            // b_score: <xsl:value-of select="./protein/peptide/domain/@b_score" />
            // b_ions: <xsl:value-of select="./protein/peptide/domain/@b_ions" />
            // Delta_Mass: <xsl:value-of select="./protein/peptide/domain/@delta" />

            // Scan: <xsl:value-of select="substring-before(concat(substring-after(./group/note,'scan='),' '), ' ')" />

            var groupIDInXMLFile = string.Empty;

            var proteinSequenceParsed = false;
            var domainParsed = false;
            var success = false;

            var currentGroupType = XTANDEM_XML_GROUP_TYPE_MODEL;
            var currentGroupLabel = GROUP_LABEL_PROTEIN;

            try
            {
                // Reset the first entry in searchResults
                searchResultCount = 0;
                searchResults[0].Clear();

                groupIDInXMLFile = XMLTextReaderGetAttributeValue(xmlReader, "id", string.Empty);
                var searchResult = searchResults[0];
                // Initially set .ResultID to groupIDInXMLFile
                // ResultID will get updated to a sequentially assigned number (mNextResultID) if we write the result out to the _xt.txt file
                searchResult.ResultID = int.Parse(groupIDInXMLFile);
                searchResult.GroupID = searchResult.ResultID;

                searchResult.ParentIonMH = XMLTextReaderGetAttributeValue(xmlReader, "mh", string.Empty);
                searchResult.Charge = XMLTextReaderGetAttributeValue(xmlReader, "z", string.Empty);

                // Note: This will get updated for each protein encountered
                searchResult.PeptideExpectationValue = ConvertEValueToBase10Log(XMLTextReaderGetAttributeValue(xmlReader, "expect", string.Empty));

                // Note: we truncate .ProteinName at the first space
                searchResult.ProteinName = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);
                var index = searchResult.ProteinName.IndexOf(' ');
                if (index > 0)
                {
                    searchResult.ProteinName = searchResult.ProteinName.Substring(0, index);
                }

                searchResult.PeptideIntensity = XMLTextReaderGetAttributeValue(xmlReader, "sumI", string.Empty);
                searchResult.PeptideIntensityMax = XMLTextReaderGetAttributeValue(xmlReader, "maxI", string.Empty);
                searchResult.fI = XMLTextReaderGetAttributeValue(xmlReader, "fI", string.Empty);

                // Continue reading the XML file, loading the information

                while (xmlReader.Read() & !AbortProcessing)
                {
                    XMLTextReaderSkipWhitespace(xmlReader);
                    if (xmlReader.ReadState != ReadState.Interactive)
                        break;

                    if (xmlReader.NodeType == XmlNodeType.Element)
                    {
                        clsSearchResultsXTandem result;
                        switch (xmlReader.Name.ToLower())
                        {
                            case XTANDEM_XML_ELEMENT_NAME_PROTEIN:
                                searchResultCount += 1;
                                if (searchResultCount > searchResults.Length)
                                {
                                    // Double the length of searchResults
                                    Array.Resize(ref searchResults, searchResults.Length * 2);
                                    for (var searchResultIndex = searchResultCount - 1; searchResultIndex <= searchResults.Length - 1; searchResultIndex++)
                                    {
                                        searchResults[searchResultIndex] = new clsSearchResultsXTandem(mPeptideMods, mPeptideSeqMassCalculator);
                                        searchResults[searchResultIndex].UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);
                                    }
                                }

                                result = searchResults[searchResultCount - 1];
                                if (searchResultCount > 1)
                                {
                                    // Copy the info from searchResults[0] to this search result
                                    // ResultID will get updated to a sequentially assigned number (mNextResultID) if we write the result out to the _xt.txt file
                                    result.ResultID = searchResults[0].ResultID;
                                    result.GroupID = searchResults[0].GroupID;
                                    result.ParentIonMH = searchResults[0].ParentIonMH;
                                    result.Charge = searchResults[0].Charge;

                                    result.PeptideExpectationValue = searchResults[0].PeptideExpectationValue;

                                    result.PeptideIntensity = searchResults[0].PeptideIntensity;
                                    result.PeptideIntensityMax = searchResults[0].PeptideIntensityMax;
                                    result.fI = searchResults[0].fI;
                                }

                                result.ProteinExpectationValue = XMLTextReaderGetAttributeValue(xmlReader, "expect", string.Empty);
                                result.ProteinIntensity = XMLTextReaderGetAttributeValue(xmlReader, "sumI", string.Empty);

                                // Update the protein name for this protein entry

                                // Note: we truncate .ProteinName at the first space
                                // However, if mLookForReverseSequenceTag = True then we need to look for ":reversed" at the end of the description
                                result.ProteinName = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);
                                result.ProteinName = TruncateProteinName(result.ProteinName);

                                // For proteins with long descriptions, the ":reversed" tag is not present in the label attribute
                                //  and is instead in the <note label="description"> element (a sub-element of the <protein> element
                                // We'll check for this case later in this function

                                // Reset the Protein Sequence Parsed and Domain Parsed flags
                                proteinSequenceParsed = false;
                                domainParsed = false;

                                // Clear the protein sequence info and peptide details info
                                result.ClearProteinSequenceInfo();
                                result.ClearPeptideDetailsInfo();

                                break;
                            case XTANDEM_XML_ELEMENT_NAME_PEPTIDE:
                                if (!proteinSequenceParsed)
                                {
                                    proteinSequenceParsed = true;
                                    result = searchResults[searchResultCount - 1];
                                    result.ProteinSeqResidueNumberStart = XMLTextReaderGetAttributeValue(xmlReader, "start", 0);
                                    result.ProteinSeqResidueNumberEnd = XMLTextReaderGetAttributeValue(xmlReader, "end", 0);
                                }
                                break;
                            case XTANDEM_XML_ELEMENT_NAME_DOMAIN:
                                // If the given peptide is present in the protein more than once, then domain will appear twice
                                //   Only keep the data for the first domain
                                // Additionally, if X!Tandem decides a given modification could occur on either of two residues, it repeats the domain information
                                //   Again, in this situation, we'll only keep the first domain

                                if (!domainParsed)
                                {
                                    domainParsed = true;

                                    // Cache the XML reader depth before reading any of the element's attributes
                                    var domainElementReaderDepth = xmlReader.Depth;

                                    // Read the information for this domain, storing in searchResult
                                    searchResult = searchResults[searchResultCount - 1];
                                    searchResult.PeptideLocInProteinStart = XMLTextReaderGetAttributeValue(xmlReader, "start", 0);
                                    searchResult.PeptideLocInProteinEnd = XMLTextReaderGetAttributeValue(xmlReader, "end", 0);

                                    // Note: the expectation value was already populated from the group level; we'll update it to the value stored for each protein in case it's different for different proteins
                                    searchResult.PeptideExpectationValue = ConvertEValueToBase10Log(XMLTextReaderGetAttributeValue(xmlReader, "expect", string.Empty));

                                    // Note: This is the theoretical, monoisotopic MH+ mass for the peptide
                                    searchResult.PeptideMH = XMLTextReaderGetAttributeValue(xmlReader, "mh", string.Empty);

                                    searchResult.PeptideDeltaMass = XMLTextReaderGetAttributeValue(xmlReader, "delta", string.Empty);

                                    // Note: .peptideDeltaMass is stored in the X!Tandem XML file as "Observed_Mass - Theoretical_Mass"
                                    // However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                                    // Therefore, we will negate .peptideDeltaMass
                                    try
                                    {
                                        searchResult.PeptideDeltaMass = (-double.Parse(searchResult.PeptideDeltaMass)).ToString(CultureInfo.InvariantCulture);
                                    }
                                    catch (Exception)
                                    {
                                        // Error; Leave .peptideDeltaMass unchanged
                                    }

                                    searchResult.PeptideHyperscore = XMLTextReaderGetAttributeValue(xmlReader, "hyperscore", string.Empty);

                                    // Note that updating .PeptideNextScore will automatically populate .DeltaCn2
                                    // ReSharper disable once StringLiteralTypo
                                    searchResult.PeptideNextScore = XMLTextReaderGetAttributeValue(xmlReader, "nextscore", string.Empty);

                                    // Note that calling .PeptidePreResidues, .PeptidePostResidues, and .PeptideCleanSequence will call ComputePeptideCleavageStateInProtein() each time
                                    searchResult.PeptidePreResidues = XMLTextReaderGetAttributeValue(xmlReader, "pre", string.Empty);
                                    searchResult.PeptidePostResidues = XMLTextReaderGetAttributeValue(xmlReader, "post", string.Empty);
                                    searchResult.PeptideCleanSequence = XMLTextReaderGetAttributeValue(xmlReader, "seq", string.Empty);

                                    searchResult.PeptideYScore = XMLTextReaderGetAttributeValue(xmlReader, "y_score", string.Empty);
                                    searchResult.PeptideYIons = XMLTextReaderGetAttributeValue(xmlReader, "y_ions", string.Empty);
                                    searchResult.PeptideBScore = XMLTextReaderGetAttributeValue(xmlReader, "b_score", string.Empty);
                                    searchResult.PeptideBIons = XMLTextReaderGetAttributeValue(xmlReader, "b_ions", string.Empty);

                                    // Now read all of the mods for this domain
                                    // If this is the first search result, then update the mod occurrence counts; otherwise, do not
                                    if (searchResultCount == 1)
                                    {
                                        ParseXTandemResultsFileReadDomainMods(xmlReader, searchResults[searchResultCount - 1], domainElementReaderDepth, true);
                                    }
                                    else
                                    {
                                        ParseXTandemResultsFileReadDomainMods(xmlReader, searchResults[searchResultCount - 1], domainElementReaderDepth, false);
                                    }
                                }

                                break;
                            case XTANDEM_XML_ELEMENT_NAME_GROUP:
                                var value = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty);
                                if (value.Length > 0)
                                {
                                    currentGroupType = string.Copy(value);
                                    currentGroupLabel = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);
                                }
                                else
                                {
                                    // Leave groupType unchanged
                                }

                                break;
                            case XTANDEM_XML_ELEMENT_NAME_NOTE:
                                if (currentGroupType == XTANDEM_XML_GROUP_TYPE_MODEL &&
                                    currentGroupLabel == GROUP_LABEL_PROTEIN && mLookForReverseSequenceTag)
                                {
                                    // Examine the label attribute of this note

                                    value = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);

                                    if (value == PROTEIN_DESCRIPTION_LABEL)
                                    {
                                        // Check whether this note ends in ":reversed"
                                        // If it does, then make sure searchResults[searchResultCount - 1].ProteinName ends in :reversed

                                        // Advance the reader before grabbing the inner text
                                        if (xmlReader.Read())
                                        {
                                            // Read the note's inner text
                                            value = XMLTextReaderGetInnerText(xmlReader);

                                            if (value.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR))
                                            {
                                                if (!searchResults[searchResultCount - 1].ProteinName.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR))
                                                {
                                                    searchResults[searchResultCount - 1].ProteinName += REVERSED_PROTEIN_SEQUENCE_INDICATOR;
                                                }
                                            }
                                        }
                                    }
                                }
                                else if (currentGroupType == XTANDEM_XML_GROUP_TYPE_SUPPORT &&
                                         currentGroupLabel == GROUP_LABEL_FRAG_ION)
                                {
                                    // This note should contain the scan number
                                    // For _Dta.txt files created at PNNL, it should look something like: "   scan=15118 cs=3"
                                    // For DTA-based files converted to .MGF and processed by X!Tandem: "MyDataset.300.300.2.dta"

                                    // Read the note's inner text
                                    value = XMLTextReaderGetInnerText(xmlReader);
                                    if (value != null)
                                    {
                                        var scanFound = false;

                                        // Look for the word "scan" followed by an equals sign, followed by a number
                                        // For example, "   scan=15118 cs=3"
                                        try
                                        {
                                            var match = mScanNumberRegExA.Match(value);
                                            if (match.Success && match.Groups.Count > 1)
                                            {
                                                searchResults[0].Scan = match.Groups[1].Value;
                                                scanFound = true;
                                            }
                                        }
                                        catch (Exception)
                                        {
                                            // Ignore errors here
                                        }

                                        if (!scanFound)
                                        {
                                            // No match; look for the word "scan" followed by whitespace, followed by a number
                                            // For example, "scan 300"
                                            try
                                            {
                                                var match = mScanNumberRegExB.Match(value);
                                                if (match.Success && match.Groups.Count > 1)
                                                {
                                                    searchResults[0].Scan = match.Groups[1].Value;
                                                    scanFound = true;
                                                }
                                            }
                                            catch (Exception)
                                            {
                                                // Ignore errors here
                                            }
                                        }

                                        if (!scanFound)
                                        {
                                            // No match; see if the description resembles a .Dta file name
                                            // For example, "MyDataset.300.300.2.dta"
                                            try
                                            {
                                                var match = mScanNumberRegExC.Match(value);
                                                if (match.Success && match.Groups.Count > 1)
                                                {
                                                    searchResults[0].Scan = match.Groups[1].Value;
                                                    scanFound = true;
                                                }
                                            }
                                            catch (Exception)
                                            {
                                                // Ignore errors here
                                            }
                                        }

                                        if (!scanFound)
                                        {
                                            // Still no match; extract out the first number present
                                            try
                                            {
                                                var match = mScanNumberRegExD.Match(value);
                                                if (match.Success)
                                                {
                                                    searchResults[0].Scan = match.Value;
                                                }
                                                else
                                                {
                                                    searchResults[0].Scan = value;
                                                }
                                            }
                                            catch (Exception)
                                            {
                                                // Ignore errors here
                                            }
                                        }

                                        // Copy the scan value from the first result to the other results
                                        for (var searchResultIndex = 1; searchResultIndex <= searchResultCount - 1; searchResultIndex++)
                                        {
                                            searchResults[searchResultIndex].Scan = searchResults[0].Scan;
                                        }
                                    }
                                }
                                break;
                            default:
                                // Unknown/unneeded child node name; ignore it
                                break;
                        }
                    }
                    else if (xmlReader.NodeType == XmlNodeType.EndElement)
                    {
                        if (xmlReader.Name == XTANDEM_XML_ELEMENT_NAME_GROUP)
                        {
                            if (xmlReader.Depth <= groupElementReaderDepth)
                            {
                                // End element found for the current group

                                // Typically each group will consist of entries all having the same sequence and modifications (but different protein names)
                                // However, occasionally a group will contain a mix of peptide sequences (only occurs if they all had the exact same hyperscore)
                                // In order to check for this, we will construct a pointer array of Sequence and Mods to SearchResultIndex and use this to determine
                                //  which entries should be written to the _xt.txt file and to the ResultToSeqMap file
                                // We will also use this pointer array to keep track of the number of proteins listed for each peptide

                                mSeqsWithMods.Clear();
                                mSeqsWithoutMods.Clear();

                                // First step through the results to compute the mass, construct the modification description,
                                //  and determine the number of proteins listed for each
                                string sequenceWithMods;
                                for (var searchResultIndex = 0; searchResultIndex <= searchResultCount - 1; searchResultIndex++)
                                {
                                    bool updateModOccurrenceCounts;
                                    if (searchResultIndex == 0)
                                    {
                                        // Always set updateModOccurrenceCounts to True for the first result in the group
                                        updateModOccurrenceCounts = true;
                                        mSeqsWithoutMods.Add(searchResults[searchResultIndex].PeptideCleanSequence);
                                    }
                                    else
                                    {
                                        if (mSeqsWithoutMods.Contains(searchResults[searchResultIndex].PeptideCleanSequence))
                                        {
                                            updateModOccurrenceCounts = false;
                                        }
                                        else
                                        {
                                            updateModOccurrenceCounts = true;
                                            mSeqsWithoutMods.Add(searchResults[searchResultIndex].PeptideCleanSequence);
                                        }
                                    }

                                    success = AddModificationsAndComputeMass(searchResults[searchResultIndex], updateModOccurrenceCounts);
                                    if (!success)
                                    {
                                        if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                                        {
                                            errorLog += "Error adding modifications to sequence for Group ID '" + groupIDInXMLFile + "'" + "\n";
                                        }
                                    }

                                    sequenceWithMods = searchResults[searchResultIndex].PeptideCleanSequence + "_" + searchResults[searchResultIndex].PeptideModDescription;

                                    if (searchResultIndex == 0)
                                    {
                                        // Always add the first result in the group mSeqsWithMods
                                        mSeqsWithMods.Add(sequenceWithMods, 1);
                                    }
                                    else
                                    {
                                        // See if mSeqsWithMods contains sequenceWithMods
                                        if (mSeqsWithMods.TryGetValue(sequenceWithMods, out var existingProteinCount))
                                        {
                                            // Increment the protein count for this peptide
                                            mSeqsWithMods[sequenceWithMods] = existingProteinCount + 1;
                                        }
                                        else
                                        {
                                            mSeqsWithMods.Add(sequenceWithMods, 1);
                                        }
                                    }
                                }

                                // Now step through the list again and update the ProteinCount value for each search result
                                for (var searchResultIndex = 0; searchResultIndex <= searchResultCount - 1; searchResultIndex++)
                                {
                                    sequenceWithMods = searchResults[searchResultIndex].PeptideCleanSequence + "_" + searchResults[searchResultIndex].PeptideModDescription;

                                    if (!mSeqsWithMods.TryGetValue(sequenceWithMods, out var proteinCount))
                                    {
                                        proteinCount = 1;
                                    }

                                    // Note: Multiple protein count is 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
                                    searchResults[searchResultIndex].MultipleProteinCount = (proteinCount - 1).ToString();
                                }

                                // Clear mSeqsWithMods again since we need to re-use it to determine which results to write out
                                mSeqsWithMods.Clear();

                                // Write out the results
                                for (var searchResultIndex = 0; searchResultIndex <= searchResultCount - 1; searchResultIndex++)
                                {
                                    sequenceWithMods = searchResults[searchResultIndex].PeptideCleanSequence + "_" + searchResults[searchResultIndex].PeptideModDescription;

                                    bool updateResultToSeqMapFile;
                                    if (searchResultIndex == 0)
                                    {
                                        // Always save the first result in the group to the _xt.txt and _ResultToSeqMap.txt files
                                        mSeqsWithMods.Add(sequenceWithMods, 1);
                                        updateResultToSeqMapFile = true;
                                    }
                                    else
                                    {
                                        // See if mSeqsWithMods contains sequenceWithMods
                                        if (mSeqsWithMods.ContainsKey(sequenceWithMods))
                                        {
                                            updateResultToSeqMapFile = false;
                                        }
                                        else
                                        {
                                            mSeqsWithMods.Add(sequenceWithMods, 1);
                                            updateResultToSeqMapFile = true;
                                        }
                                    }

                                    if (updateResultToSeqMapFile)
                                    {
                                        // Only save the first result for each peptide in the group to the _xt.txt and _ResultToSeqMap.txt files
                                        // Note: This function will update .ResultID to the next available ID value (mNextResultID)
                                        SaveXTandemResultsFileEntry(searchResults[searchResultIndex], ref writer);
                                    }

                                    SaveResultsFileEntrySeqInfo(searchResults[searchResultIndex], updateResultToSeqMapFile);
                                }

                                success = true;
                                break;
                            }
                        }
                    }
                }
            }
            catch (Exception)
            {
                // Error parsing values from this group ID in the XML file
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error parsing value for Group ID '" + groupIDInXMLFile + "'" + "\n";
                }
                success = false;
            }

            return success;
        }

        private bool ParseXTandemResultsFileInputParameters(string inputFilePath)
        {
            // Pre-read the XML file and look for the Input Parameters section
            // Read the parameters and validate that each of the mods defined is present in mPeptideMods

            const string GROUP_LABEL_INPUT_PARAMETERS = "input parameters";

            bool success;

            try
            {
                // Open the input file and parse it
                // Initialize the stream reader and the XML Text Reader
                eCurrentXMLDataFileSectionConstants eCurrentXMLDataFileSection;
                using (var xmlReader = new XmlTextReader(inputFilePath))
                {
                    // Parse the file
                    eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile;

                    while (xmlReader.Read() & !AbortProcessing)
                    {
                        XMLTextReaderSkipWhitespace(xmlReader);
                        if (xmlReader.ReadState != ReadState.Interactive)
                            break;

                        if (xmlReader.Depth < 2)
                        {
                            if (xmlReader.NodeType == XmlNodeType.Element)
                            {
                                switch (xmlReader.Name.ToLower())
                                {
                                    case XTANDEM_XML_ELEMENT_NAME_GROUP:
                                        if (xmlReader.HasAttributes)
                                        {
                                            var parametersGroupDepth = xmlReader.Depth;

                                            // See if the group has a "type" attribute containing the text XTANDEM_XML_GROUP_TYPE_PARAMETERS
                                            var currentGroupType = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty);
                                            if (currentGroupType == XTANDEM_XML_GROUP_TYPE_PARAMETERS)
                                            {
                                                eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.InputParameters;

                                                // Read the Label for this group
                                                var currentGroupLabel = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);
                                                if (currentGroupLabel == GROUP_LABEL_INPUT_PARAMETERS)
                                                {
                                                    // Read the input parameters
                                                    ParseXTandemResultsFileInputParametersWork(xmlReader, parametersGroupDepth);
                                                }
                                            }
                                            else
                                            {
                                                // Skip this group
                                                xmlReader.Skip();
                                            }
                                        }
                                        else
                                        {
                                            // Group doesn't have any attributes; ignore it
                                            xmlReader.Skip();
                                        }

                                        break;
                                    case XTANDEM_XML_ROOT_ELEMENT:
                                        eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.Start;
                                        break;
                                    default:
                                        // Skip this element
                                        xmlReader.Skip();
                                        break;
                                }
                            }
                        }
                    }
                }

                if (eCurrentXMLDataFileSection == eCurrentXMLDataFileSectionConstants.UnknownFile)
                {
                    SetErrorMessage("Root element '" + XTANDEM_XML_ROOT_ELEMENT + "' not found in the input file: " + inputFilePath);
                    success = false;
                }
                else
                {
                    success = true;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                success = false;
            }

            return success;
        }

        private void ParseXTandemResultsFileInputParametersWork(XmlReader xmlReader, int parametersGroupDepth)
        {
            // Read the input parameters
            // Each parameter is an element with name "note" with attributes "type" and "label"
            // The parameter value is the text between within the element

            const string NOTE_TYPE_INPUT = "input";
            const char CLEAVAGE_SPEC_SEP = '|';
            const char XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START = '{';
            const char XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END = '}';

            // Initialize the Mod Info array
            var modInfoCount = 0;
            var udtModInfo = new udtSearchOptionModificationInfoType[20];

            var staticModsAreResetForRefinement = true;

            // Initialize udtParamLabels; this specifies the parameters to examine
            // Note: When populating this we use .ToLower() to make sure all of the text is lowercase
            var udtParamLabels = new string[INPUT_PARAM_LABEL_NAMES_MAX_INDEX + 1];
            udtParamLabels[(int)eInputParamLabelNames.Residue_StaticModMass] = "residue, modification mass".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Residue_PotentialModMass] = "residue, potential modification mass".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Residue_PotentialModMotif] = "residue, potential modification motif".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialModMass] = "refine, potential modification mass".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialModMotif] = "refine, potential modification motif".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialNTerminusMods] = "refine, potential N-terminus modifications".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialCTerminusMods] = "refine, potential C-terminus modifications".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Protein_NTerminal_ResidueModMass] = "protein, N-terminal residue modification mass".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Protein_CTerminal_ResidueModMass] = "protein, C-terminal residue modification mass".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_NTerminalMassChange] = "protein, cleavage N-terminal mass change".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_CTerminalMassChange] = "protein, cleavage C-terminal mass change".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_Site] = "protein, cleavage site".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Refine_ModificationMass] = "refine, modification mass".ToLower();
            udtParamLabels[(int)eInputParamLabelNames.Scoring_Include_Reverse] = "scoring, include reverse".ToLower();

            // Make sure all of the text in udtParamLabels() is lowercase
            for (var index = 0; index <= udtParamLabels.Length - 1; index++)
            {
                udtParamLabels[index] = udtParamLabels[index].ToLower();
            }

            while (xmlReader.Read())
            {
                XMLTextReaderSkipWhitespace(xmlReader);
                if (xmlReader.ReadState != ReadState.Interactive)
                    break;

                if (xmlReader.NodeType == XmlNodeType.Element)
                {
                    switch (xmlReader.Name.ToLower())
                    {
                        case XTANDEM_XML_ELEMENT_NAME_NOTE:
                            // Read the note's type
                            var noteType = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty);

                            if (noteType == NOTE_TYPE_INPUT)
                            {
                                // Read the note's label and inner text
                                var noteLabel = XMLTextReaderGetAttributeValue(xmlReader, "label", string.Empty);

                                // Need to advance the reader before calling XMLTextReaderGetInnerText
                                xmlReader.Read();
                                var value = XMLTextReaderGetInnerText(xmlReader);

                                if (value != null)
                                {
                                    var noteLabelLower = noteLabel.ToLower();
                                    if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Residue_StaticModMass]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.StaticMod, Convert.ToInt32(eInputParamLabelNames.Residue_StaticModMass), false, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Residue_PotentialModMass]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Residue_PotentialModMass), false, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Residue_PotentialModMotif]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Residue_PotentialModMotif), true, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialModMass]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialModMass), false, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialModMotif]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialModMotif), true, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialNTerminusMods]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialNTerminusMods), false, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialCTerminusMods]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialCTerminusMods), false, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_NTerminal_ResidueModMass]))
                                    {
                                        ParseXTandemInputParameterProteinTerminusMod(Convert.ToInt32(eInputParamLabelNames.Protein_NTerminal_ResidueModMass), true, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_CTerminal_ResidueModMass]))
                                    {
                                        ParseXTandemInputParameterProteinTerminusMod(Convert.ToInt32(eInputParamLabelNames.Protein_CTerminal_ResidueModMass), false, value, ref modInfoCount, ref udtModInfo);
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_NTerminalMassChange]))
                                    {
                                        if (clsPHRPParser.IsNumber(value))
                                        {
                                            PeptideNTerminusMassChange = double.Parse(value);
                                        }
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_CTerminalMassChange]))
                                    {
                                        if (clsPHRPParser.IsNumber(value))
                                        {
                                            PeptideCTerminusMassChange = double.Parse(value);
                                        }
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_Site]))
                                    {
                                        // In X!Tandem the LeftSpec and RightSpec values are separated by a vertical bar (CLEAVAGE_SPEC_SEP)
                                        // Look for CLEAVAGE_SPEC_SEP in value
                                        var barLoc = value.IndexOf(CLEAVAGE_SPEC_SEP);
                                        if (barLoc > 0 & barLoc < value.Length - 1)
                                        {
                                            var leftSpec = value.Substring(0, barLoc);
                                            var rightSpec = value.Substring(barLoc + 1);

                                            // Look for curly braces in leftSpec and rightSpec
                                            // In X!Tandem curly braces mean to not match a residue
                                            // If found, change to standard RegEx notation, e.g. from {P} to [^P]
                                            if (leftSpec.IndexOf(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START) >= 0)
                                            {
                                                leftSpec = leftSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START.ToString(), "[^");
                                                leftSpec = leftSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END.ToString(), "]");
                                            }

                                            if (rightSpec.IndexOf(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START) >= 0)
                                            {
                                                rightSpec = rightSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START.ToString(), "[^");
                                                rightSpec = rightSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END.ToString(), "]");
                                            }

                                            EnzymeMatchSpec = new clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType(leftSpec, rightSpec);
                                        }
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_ModificationMass]))
                                    {
                                        if (!string.IsNullOrWhiteSpace(value))
                                        {
                                            staticModsAreResetForRefinement = true;
                                        }
                                    }
                                    else if (noteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Scoring_Include_Reverse]))
                                    {
                                        if (!string.IsNullOrWhiteSpace(value))
                                        {
                                            if (value.Trim().ToLower() == "yes")
                                            {
                                                mLookForReverseSequenceTag = true;
                                            }
                                        }
                                    }
                                }
                            }

                            break;
                        default:
                            // Unknown/unneeded child node name; ignore it
                            break;
                    }
                }
                else if (xmlReader.NodeType == XmlNodeType.EndElement)
                {
                    if (xmlReader.Depth == parametersGroupDepth && xmlReader.Name == XTANDEM_XML_ELEMENT_NAME_GROUP)
                    {
                        // Reached the end of this group
                        break;
                    }
                }
            }

            if (modInfoCount > 0)
            {
                // Validate that each of the mods in udtModInfo is present in mPeptideMods
                if (modInfoCount > 1)
                {
                    // Sort udtModInfo
                    Array.Sort(udtModInfo, 0, modInfoCount, new ISearchOptionModificationInfoComparer());
                }

                // Before continuing, look for Static residue mods in udtModInfo
                // If any are found, and if an identical dynamic residue mod is already present, then delete the static residue mod
                // Additionally, if 	<note type="input" label="refine, modification mass">none</note> was present in the XTandem results file then
                //  auto update all static mods to dynamic mods since they are reset during refinement
                for (var index = 0; index < modInfoCount;)
                {
                    var modDeleted = false;
                    if (udtModInfo[index].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                    {
                        var indexCompare = 0;
                        while (indexCompare < modInfoCount)
                        {
                            if (indexCompare != index && udtModInfo[indexCompare].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod)
                            {
                                // See if this modification has a similar mass (within MASS_DIGITS_OF_PRECISION digits of precision)
                                if (Math.Abs(Math.Round(Math.Abs(udtModInfo[indexCompare].ModificationMass - udtModInfo[index].ModificationMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)) < float.Epsilon)
                                {
                                    // Matching mass
                                    // Compare .TargetResidues
                                    if (clsModificationDefinition.EquivalentTargetResidues(udtModInfo[indexCompare].TargetResidues, udtModInfo[index].TargetResidues, true))
                                    {
                                        // Yes, the modifications match; delete the static version of the modification
                                        for (var indexCopy = index; indexCopy <= modInfoCount - 2; indexCopy++)
                                        {
                                            udtModInfo[indexCopy] = udtModInfo[indexCopy + 1];
                                        }
                                        modInfoCount -= 1;
                                        modDeleted = true;
                                        break;
                                    }
                                }
                            }
                            indexCompare += 1;
                        }
                    }

                    if (!modDeleted)
                    {
                        if (udtModInfo[index].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod && staticModsAreResetForRefinement)
                        {
                            // Update this static mod to be a dynamic mod
                            udtModInfo[index].ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod;
                        }
                        index += 1;
                    }
                }

                for (var index = 0; index <= modInfoCount - 1; index++)
                {
                    var modInfo = udtModInfo[index];
                    mPeptideMods.VerifyModificationPresent(modInfo.ModificationMass, modInfo.TargetResidues, modInfo.ModificationType);
                }
            }

            // In addition, verify that the standard refinement modifications are present in mPeptideMods
            mPeptideMods.AppendStandardRefinementModifications();
        }

        private void ParseXTandemResultsFileReadDomainMods(XmlReader xmlReader, clsSearchResultsBaseClass searchResult, int domainElementReaderDepth, bool updateModOccurrenceCounts)
        {
            // Continue reading the XML file, loading the information

            while (xmlReader.Read() && xmlReader.Depth >= domainElementReaderDepth)
            {
                if (xmlReader.ReadState != ReadState.Interactive)
                    break;

                if (xmlReader.NodeType == XmlNodeType.Element)
                {
                    switch (xmlReader.Name.ToLower())
                    {
                        case XTANDEM_XML_ELEMENT_NAME_AMINO_ACID:
                            var value = XMLTextReaderGetAttributeValue(xmlReader, "type", string.Empty).Trim();
                            char targetResidue;
                            if (string.IsNullOrWhiteSpace(value))
                            {
                                targetResidue = default(char);
                            }
                            else
                            {
                                targetResidue = value[0];
                            }

                            var modifiedResiduePosInProtein = XMLTextReaderGetAttributeValue(xmlReader, "at", 0);

                            if (modifiedResiduePosInProtein > 0)
                            {
                                var modificationMass = XMLTextReaderGetAttributeValueDbl(xmlReader, "modified", 0);

                                if (Math.Abs(modificationMass - 0) > float.Epsilon)
                                {
                                    var residueLocInPeptide = modifiedResiduePosInProtein - searchResult.PeptideLocInProteinStart + 1;
                                    var eResidueTerminusState = searchResult.DetermineResidueTerminusState(residueLocInPeptide);

                                    searchResult.SearchResultAddModification(modificationMass,
                                                                             targetResidue,
                                                                             residueLocInPeptide,
                                                                             eResidueTerminusState,
                                                                             updateModOccurrenceCounts);
                                }
                            }
                            break;
                    }
                }
                else if (xmlReader.NodeType == XmlNodeType.EndElement)
                {
                    if (xmlReader.Name == XTANDEM_XML_ELEMENT_NAME_DOMAIN)
                    {
                        break;
                    }
                }
            }
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="inputFilePath">X!Tandem results file</param>
        /// <param name="outputDirectoryPath">Output directory</param>
        /// <param name="parameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string inputFilePath, string outputDirectoryPath, string parameterFilePath)
        {
            var success = false;

            if (!LoadParameterFileSettings(parameterFilePath))
            {
                SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(inputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                success = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!success)
                {
                    return false;
                }

                ResetProgress("Parsing " + Path.GetFileName(inputFilePath));

                if (!CleanupFilePaths(ref inputFilePath, ref outputDirectoryPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(inputFilePath);

                    // Define the output file name based on inputFilePath
                    // The name will be DatasetName_xt.txt

                    var xtandemXTFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, ".txt"));
                    xtandemXTFilePath = Path.Combine(outputDirectoryPath, xtandemXTFilePath);
                    success = ParseXTandemResultsFile(inputFile.FullName, xtandemXTFilePath, false);

                    if (success && CreateProteinModsFile)
                    {
                        success = CreateProteinModsFileWork(inputFile, outputDirectoryPath, xtandemXTFilePath);
                    }

                    if (success)
                    {
                        OperationComplete();
                    }
                }
                catch (Exception ex)
                {
                    SetErrorMessage("Error in clsXTandemResultsProcessor.ProcessFile (2): " + ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage("Error in clsXTandemResultsProcessor.ProcessFile (1):" + ex.Message);
                SetErrorCode(ePHRPErrorCodes.UnspecifiedError);
            }

            return success;
        }

        private bool CreateProteinModsFileWork(FileInfo inputFile, string outputDirectoryPath, string xtandemXTFilePath)
        {
            bool success;

            // First create the MTS PepToProteinMap file using inputFile
            var sourcePHRPDataFiles = new List<string> {
                xtandemXTFilePath
            };

            var mtsPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, outputDirectoryPath, mts: true);

            if (File.Exists(mtsPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
            {
                success = true;
            }
            else
            {
                success = CreatePepToProteinMapFile(sourcePHRPDataFiles, mtsPepToProteinMapFilePath);
                if (!success)
                {
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                }
            }

            if (success)
            {
                if (inputFile.Directory == null)
                {
                    ReportWarning("CreateProteinModsFileWork: Could not determine the parent directory of " + inputFile.FullName);
                }
                else if (string.IsNullOrWhiteSpace(xtandemXTFilePath))
                {
                    ReportWarning("CreateProteinModsFileWork: xtandemXTFilePath is null; cannot call CreateProteinModDetailsFile");
                }
                else
                {
                    // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output directory
                    ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.Directory.FullName, Path.GetFileName(xtandemXTFilePath)), outputDirectoryPath);

                    // Now create the Protein Mods file
                    success = CreateProteinModDetailsFile(xtandemXTFilePath, outputDirectoryPath, mtsPepToProteinMapFilePath,
                                                          clsPHRPReader.ePeptideHitResultType.XTandem);
                }
            }

            if (!success)
            {
                // Do not treat this as a fatal error
                return true;
            }
            return true;
        }

        private void SaveXTandemResultsFileEntry(clsSearchResultsXTandem searchResult, ref StreamWriter writer)
        {
            // Update .ResultID to the next available number
            searchResult.ResultID = mNextResultID;
            mNextResultID += 1;

            // Write the results to the output file
            var data = new List<string>
            {
                searchResult.ResultID.ToString(),
                searchResult.GroupID.ToString(),
                searchResult.Scan,
                searchResult.Charge,
                searchResult.PeptideMH,
                searchResult.PeptideHyperscore,
                searchResult.PeptideExpectationValue,
                searchResult.MultipleProteinCount,
                searchResult.SequenceWithPrefixAndSuffix(true),
                Math.Round(searchResult.PeptideDeltaCn2, 4).ToString(CultureInfo.InvariantCulture),
                searchResult.PeptideYScore,
                searchResult.PeptideYIons,
                searchResult.PeptideBScore,
                searchResult.PeptideBIons,
                searchResult.PeptideDeltaMass,
                searchResult.PeptideIntensity,
                PRISM.StringUtilities.DblToString(searchResult.PeptideDeltaMassCorrectedPpm, 5, 0.00005)
            };

            writer.WriteLine(CollapseList(data));
        }

        protected override string TruncateProteinName(string proteinNameAndDescription)
        {
            var isReversed = false;

            if (mLookForReverseSequenceTag)
            {
                isReversed = proteinNameAndDescription.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR);
            }

            proteinNameAndDescription = base.TruncateProteinName(proteinNameAndDescription);

            if (isReversed)
            {
                return proteinNameAndDescription + REVERSED_PROTEIN_SEQUENCE_INDICATOR;
            }
            return proteinNameAndDescription;
        }

        /// <summary>
        ///  Write out the header line for synopsis / first hits file
        /// </summary>
        /// <param name="writer"></param>
        /// <param name="errorLog"></param>
        private void WriteSynFHTFileHeader(
            TextWriter writer,
            ref string errorLog)
        {
            try
            {
                // Get the synopsis file headers
                // Keys are header name and values are enum IDs
                var headerColumns = clsPHRPParserXTandem.GetColumnHeaderNamesAndIDs();

                var headerNames = (from item in headerColumns orderby item.Value select item.Key).ToList();

                writer.WriteLine(CollapseList(headerNames));
            }
            catch (Exception)
            {
                if (errorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    errorLog += "Error writing synopsis / first hits header" + "\n";
                }
            }
        }

        private string XMLTextReaderGetAttributeValue(XmlReader xmlReader, string attributeName, string valueIfMissing)
        {
            xmlReader.MoveToAttribute(attributeName);
            if (xmlReader.ReadAttributeValue())
            {
                return xmlReader.Value;
            }
            return string.Copy(valueIfMissing);
        }

        private int XMLTextReaderGetAttributeValue(XmlReader xmlReader, string attributeName, int valueIfMissing)
        {
            xmlReader.MoveToAttribute(attributeName);
            if (xmlReader.ReadAttributeValue())
            {
                if (clsPHRPParser.IsNumber(xmlReader.Value))
                {
                    return Convert.ToInt32(xmlReader.Value);
                }
                return valueIfMissing;
            }
            return valueIfMissing;
        }

        private double XMLTextReaderGetAttributeValueDbl(XmlReader xmlReader, string attributeName, double valueIfMissing)
        {
            xmlReader.MoveToAttribute(attributeName);
            if (xmlReader.ReadAttributeValue())
            {
                if (clsPHRPParser.IsNumber(xmlReader.Value))
                {
                    return Convert.ToDouble(xmlReader.Value);
                }
                return valueIfMissing;
            }
            return valueIfMissing;
        }

        private string XMLTextReaderGetInnerText(XmlReader xmlReader)
        {
            var value = string.Empty;
            bool success;

            if (xmlReader.NodeType == XmlNodeType.Element)
            {
                // Advance the reader so that we can read the value
                success = xmlReader.Read();
            }
            else
            {
                success = true;
            }

            if (success && xmlReader.NodeType != XmlNodeType.Whitespace & xmlReader.HasValue)
            {
                value = xmlReader.Value;
            }

            return value;
        }

        private void XMLTextReaderSkipWhitespace(XmlReader xmlReader)
        {
            if (xmlReader.NodeType == XmlNodeType.Whitespace)
            {
                // Whitespace; read the next node
                xmlReader.Read();
            }
        }
    }
}
