// This class reads in an X!Tandem results file (XML format) and creates
// a tab-delimited text file with the data.  It will insert modification symbols
// into the peptide sequences for modified peptides.  The user can optionally provide
// a modification definition file which specifies the symbol to use for each
// modification mass.  If the user does not provide this file, then the modification
// definition information is determined from the X!Tandem Input Parameters section at
// the end of the X!Tandem results file.
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 2, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.
using System;
using System.Collections;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text.RegularExpressions;
using System.Xml;
using PHRPReader;

namespace PeptideHitResultsProcessor
{
    public class clsXTandemResultsProcessor : clsPHRPBaseClass
    {
        public clsXTandemResultsProcessor()
        {
            mFileDate = "October 13, 2017";
            InitializeLocalVariables();
        }

        #region "Constants and Enums"
        // Note: These names must all be lowercase
        private const string XTANDEM_XML_ROOT_ELEMENT = "bioml";
        private const string XTANDEM_XML_ELEMENT_NAME_GROUP = "group";
        private const string XTANDEM_XML_ELEMENT_NAME_PROTEIN = "protein";
        private const string XTANDEM_XML_ELEMENT_NAME_PEPTIDE = "peptide";
        private const string XTANDEM_XML_ELEMENT_NAME_DOMAIN = "domain";
        private const string XTANDEM_XML_ELEMENT_NAME_AMINOACID = "aa";
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
        protected int mNextResultID;
        protected bool mLookForReverseSequenceTag;

        private Regex mScanNumberRegExA;
        private Regex mScanNumberRegExB;
        private Regex mScanNumberRegExC;
        private Regex mScanNumberRegExD;
        #endregion

        #region "Properties"
        #endregion

        private bool AddModificationsAndComputeMass(clsSearchResultsXTandem objSearchResult, bool blnUpdateModOccurrenceCounts)
        {
            const bool ALLOW_DUPLICATE_MOD_ON_TERMINUS = false;
            bool blnSuccess;

            try
            {
                // If any modifications of type IsotopicMod are defined we would add them to the Search Result Mods now
                // However, since X!Tandem doesn't support Isotopic Mods, this step is currently skipped
                //
                // objSearchResult.SearchResultAddIsotopicModifications(blnUpdateModOccurrenceCounts)

                // Add the protein terminus static mods (if defined and if the peptide is at a protein terminus)
                // Function .SearchResultAddStaticTerminusMods() will only add the terminus mod if the terminus
                //  is not already modified by the given terminus mod mass
                // This function will also add peptide terminus static mods, if defined, though those are not supported by X!Tandem and therefore should not be defined
                objSearchResult.SearchResultAddStaticTerminusMods(ALLOW_DUPLICATE_MOD_ON_TERMINUS, blnUpdateModOccurrenceCounts);

                // Compute the monoisotopic mass for this peptide
                objSearchResult.ComputeMonoisotopicMass();

                // Update PeptideDeltaMassCorrectedPpm
                objSearchResult.ComputeDelMCorrectedXT();

                // Populate .PeptideSequenceWithMods and .PeptideModDescription
                // Note that this function will call .AddSearchResultModificationsToCleanSequence() then .UpdateModDescription()
                objSearchResult.ApplyModificationInformation();

                blnSuccess = true;
            }
            catch (Exception)
            {
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private string ConvertEValueToBase10Log(string strExpectationValue)
        {
            double dblEValue;
            double dblLogEValue;

            try
            {
                dblEValue = double.Parse(strExpectationValue);
            }
            catch (Exception)
            {
                dblEValue = 0;
            }

            try
            {
                if (dblEValue <= 0)
                {
                    dblLogEValue = 0;
                }
                else if (dblEValue <= 1E-307)
                {
                    dblLogEValue = -307;
                }
                else
                {
                    dblLogEValue = Math.Round(Math.Log10(dblEValue), 3);
                }
            }
            catch (Exception)
            {
                dblLogEValue = 0;
            }

            return dblLogEValue.ToString("0.000");
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

        private bool ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants eModificationType, int intSortOrder, bool blnParsingMotifDef, string strParamValue, ref int intModInfoCount, ref udtSearchOptionModificationInfoType[] udtModInfo)
        {
            // Parse out the mod information defined in strParamValue
            // Add each entry to udtModInfo
            // If blnParsingMotifDef = True, then do not try to determine .TargetResidues

            const char MOD_LIST_SEP_CHAR = ',';

            bool blnSuccess;

            try
            {
                if (string.IsNullOrWhiteSpace(strParamValue))
                {
                    // Empty parameter; no definition to parse
                }
                else
                {
                    // Parse strParamValue
                    // Each modification entry can have multiple modification definitions separated separated by commas, so we first split strParamValue
                    var strModDefs = strParamValue.Split(MOD_LIST_SEP_CHAR);
                    for (var intIndex = 0; intIndex <= strModDefs.Length - 1; intIndex++)
                    {
                        // Modification definitions typically look like "15.9949@M"
                        // However, a neutral loss can be specified using "79.9663:-97.98@STY"
                        //   Thus, mod mass is number up to first non-numeric character (typically a colon or @ sign)
                        // Target residues are the residues after the @
                        // If the Target residues contain an X, then the modification can apply to any residue

                        var dblModificationMass = 0.0;

                        // Look for a colon and an @ sign
                        var intColonIndex = strModDefs[intIndex].IndexOf(':');
                        var intAtSignIndex = strModDefs[intIndex].IndexOf('@');

                        if (intAtSignIndex < 0)
                        {
                            // At sign not found; skip this mod def
                        }
                        else
                        {
                            if (intColonIndex > intAtSignIndex)
                            {
                                // Ignore this colon since it's present after the @ sign
                                intColonIndex = -1;
                            }

                            if (intColonIndex > 0)
                            {
                                // Colon found; see if the text up to intColonIndex is a number
                                if (clsPHRPParser.IsNumber(strModDefs[intIndex].Substring(0, intColonIndex)))
                                {
                                    dblModificationMass = double.Parse(strModDefs[intIndex].Substring(0, intColonIndex));
                                }
                            }
                            else
                            {
                                // Colon not found; see if the text up to intAtSignIndex is a number
                                if (clsPHRPParser.IsNumber(strModDefs[intIndex].Substring(0, intAtSignIndex)))
                                {
                                    dblModificationMass = double.Parse(strModDefs[intIndex].Substring(0, intAtSignIndex));
                                }
                            }

                            if (Math.Abs(dblModificationMass) > float.Epsilon)
                            {
                                string strTargetResidues;

                                // Valid mass found; now extract the target residues
                                if (!blnParsingMotifDef && intAtSignIndex + 1 < strModDefs[intIndex].Length)
                                {
                                    strTargetResidues = strModDefs[intIndex].Substring(intAtSignIndex + 1);
                                }
                                else
                                {
                                    strTargetResidues = string.Empty;
                                }

                                if (strTargetResidues.IndexOf('X') >= 0)
                                {
                                    // Modification can affect any residue; set strTargetResidues to ""
                                    strTargetResidues = string.Empty;
                                }

                                if (strTargetResidues.Length > 0)
                                {
                                    // Convert from X!Tandem-style N-Terminus notation to DMS-style notation
                                    strTargetResidues = strTargetResidues.Replace(clsPeptideModificationContainer.N_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM, clsAminoAcidModInfo.N_TERMINAL_PEPTIDE_SYMBOL_DMS);
                                    strTargetResidues = strTargetResidues.Replace(clsPeptideModificationContainer.C_TERMINAL_PEPTIDE_MOD_SYMBOL_XTANDEM, clsAminoAcidModInfo.C_TERMINAL_PEPTIDE_SYMBOL_DMS);
                                }

                                // Append the new mod information to udtModInfo
                                if (intModInfoCount >= udtModInfo.Length)
                                {
                                    Array.Resize(ref udtModInfo, udtModInfo.Length * 2);
                                }

                                var modInfo = udtModInfo[intModInfoCount];
                                modInfo.SortOrder = intSortOrder;
                                modInfo.ModificationMass = dblModificationMass;
                                modInfo.TargetResidues = strTargetResidues;
                                modInfo.ModificationType = eModificationType;
                                udtModInfo[intModInfoCount] = modInfo;
                                intModInfoCount += 1;
                            }
                        }
                    }
                }

                blnSuccess = true;
            }
            catch (Exception)
            {
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private bool ParseXTandemInputParameterProteinTerminusMod(int intSortOrder, bool blnNTerminus, string strParamValue, ref int intModInfoCount, ref udtSearchOptionModificationInfoType[] udtModInfo)
        {
            // Parse out the mass defined in strParamValue
            // Add the entry to udtModInfo if non-zero
            // if blnNTerminus = True then mod applies to the protein's N-Terminus; otherwise, applies to the protein's C-Terminus

            bool blnSuccess;

            try
            {
                if (string.IsNullOrWhiteSpace(strParamValue))
                {
                    // Empty parameter; no definition to parse
                }
                else
                {
                    // See if strParamValue is a non-zero number

                    var dblModificationMass = 0.0;
                    if (clsPHRPParser.IsNumber(strParamValue))
                    {
                        dblModificationMass = double.Parse(strParamValue);
                    }

                    if (Math.Abs(dblModificationMass) > float.Epsilon)
                    {
                        // Append the new mod information to udtModInfo
                        if (intModInfoCount >= udtModInfo.Length)
                        {
                            Array.Resize(ref udtModInfo, udtModInfo.Length * 2);
                        }

                        var modInfo = udtModInfo[intModInfoCount];
                        modInfo.SortOrder = intSortOrder;
                        string strTargetResidues;
                        if (blnNTerminus)
                        {
                            strTargetResidues = clsAminoAcidModInfo.N_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                        }
                        else
                        {
                            strTargetResidues = clsAminoAcidModInfo.C_TERMINAL_PROTEIN_SYMBOL_DMS.ToString();
                        }

                        modInfo.ModificationMass = dblModificationMass;
                        modInfo.TargetResidues = strTargetResidues;
                        modInfo.ModificationType = clsModificationDefinition.eModificationTypeConstants.ProteinTerminusStaticMod;
                        udtModInfo[intModInfoCount] = modInfo;
                        intModInfoCount += 1;
                    }
                }

                blnSuccess = true;
            }
            catch (Exception)
            {
                blnSuccess = false;
            }

            return blnSuccess;
        }

        protected bool ParseXTandemResultsFile(string strInputFilePath, string strOutputFilePath, bool blnResetMassCorrectionTagsAndModificationDefinitions = true)
        {
            // Warning: This function does not call LoadParameterFile; you should typically call ProcessFile

            var intResultsProcessed = 0;

            bool blnSuccess;

            try
            {
                // Possibly reset the mass correction tags and Mod Definitions
                if (blnResetMassCorrectionTagsAndModificationDefinitions)
                {
                    ResetMassCorrectionTagsAndModificationDefinitions();
                }

                // Reset .OccurrenceCount
                mPeptideMods.ResetOccurrenceCountStats();

                // Initialize the objSearchResults() array; initially reserve space for 4 proteins
                // Note: The number of valid entries in objSearchResults[) is given by intSearchResultCount; objSearchResults(] is expanded but never shrunk
                // There is a separate entry in objSearchResults() for each protein encountered

                var intSearchResultCount = 0;
                var objSearchResults = new clsSearchResultsXTandem[4];

                int intSearchResultIndex;
                for (intSearchResultIndex = 0; intSearchResultIndex <= objSearchResults.Length - 1; intSearchResultIndex++)
                {
                    objSearchResults[intSearchResultIndex] = new clsSearchResultsXTandem(mPeptideMods, mPeptideSeqMassCalculator);
                }

                // Reset mNextResultID and mLookForReverseSequenceTag
                InitializeLocalVariables();

                try
                {
                    // Read the input parameters from the end of the X!Tandem results file (strInputFilePath)
                    blnSuccess = ParseXTandemResultsFileInputParameters(strInputFilePath);
                    if (!blnSuccess)
                    {
                        SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile, true);
                        return false;
                    }

                    for (intSearchResultIndex = 0; intSearchResultIndex <= objSearchResults.Length - 1; intSearchResultIndex++)
                    {
                        objSearchResults[intSearchResultIndex].UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);
                    }

                    var strErrorLog = string.Empty;

                    // Open the input file and parse it

                    // Initialize the stream reader and the XML Text Reader
                    eCurrentXMLDataFileSectionConstants eCurrentXMLDataFileSection;
                    using (var srDataFile = new StreamReader(strInputFilePath))
                    {
                        using (var objXMLReader = new XmlTextReader(srDataFile))
                        {
                            intResultsProcessed = 0;

                            // Create the output file
                            using (var swPeptideResultsFile = new StreamWriter(strOutputFilePath, false))
                            {
                                // Write the header line to swPeptideResultsFile
                                swPeptideResultsFile.WriteLine("Result_ID" + SEP_CHAR +
                                                               "Group_ID" + SEP_CHAR +
                                                               "Scan" + SEP_CHAR +
                                                               "Charge" + SEP_CHAR +
                                                               "Peptide_MH" + SEP_CHAR +
                                                               "Peptide_Hyperscore" + SEP_CHAR +
                                                               "Peptide_Expectation_Value_Log(e)" + SEP_CHAR +
                                                               "Multiple_Protein_Count" + SEP_CHAR +
                                                               "Peptide_Sequence" + SEP_CHAR +
                                                               "DeltaCn2" + SEP_CHAR +
                                                               "y_score" + SEP_CHAR +
                                                               "y_ions" + SEP_CHAR +
                                                               "b_score" + SEP_CHAR +
                                                               "b_ions" + SEP_CHAR +
                                                               "Delta_Mass" + SEP_CHAR +
                                                               "Peptide_Intensity_Log(I)" + SEP_CHAR +
                                                               "DelM_PPM");

                                // Create the additional output files
                                blnSuccess = InitializeSequenceOutputFiles(strOutputFilePath);

                                // Parse the input file
                                eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile;

                                while (objXMLReader.Read() & !AbortProcessing)
                                {
                                    XMLTextReaderSkipWhitespace(objXMLReader);
                                    if (objXMLReader.ReadState != ReadState.Interactive)
                                        break;

                                    if (objXMLReader.Depth < 2)
                                    {
                                        if (objXMLReader.NodeType == XmlNodeType.Element)
                                        {
                                            switch (objXMLReader.Name.ToLower())
                                            {
                                                case XTANDEM_XML_ELEMENT_NAME_GROUP:
                                                    if (objXMLReader.HasAttributes)
                                                    {
                                                        // Cache the XML reader depth before reading any of the element's attributes
                                                        var intGroupElementReaderDepth = objXMLReader.Depth;

                                                        // See if the group has a "type" attribute containing the text XTANDEM_XML_GROUP_TYPE_MODEL
                                                        var strCurrentGroupType = XMLTextReaderGetAttributeValue(objXMLReader, "type", string.Empty);
                                                        if (strCurrentGroupType == XTANDEM_XML_GROUP_TYPE_MODEL)
                                                        {
                                                            eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.SearchResults;

                                                            ParseXTandemResultsFileEntry(objXMLReader, swPeptideResultsFile, ref intSearchResultCount, ref objSearchResults, ref strErrorLog, intGroupElementReaderDepth);
                                                            intResultsProcessed += 1;

                                                            // Update the progress
                                                            var sngPercentComplete = Convert.ToSingle(srDataFile.BaseStream.Position / srDataFile.BaseStream.Length * 100);
                                                            if (CreateProteinModsFile)
                                                            {
                                                                sngPercentComplete = sngPercentComplete * (PROGRESS_PERCENT_CREATING_PEP_TO_PROTEIN_MAPPING_FILE / 100);
                                                            }
                                                            UpdateProgress(sngPercentComplete);
                                                        }
                                                    }
                                                    else
                                                    {
                                                        // Group doesn't have any attributes; ignore it
                                                        objXMLReader.Skip();
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
                        mErrorMessage = "Root element '" + XTANDEM_XML_ROOT_ELEMENT + "' not found in the input file: " + "\n" + strInputFilePath;
                    }
                    else
                    {
                        if (CreateModificationSummaryFile)
                        {
                            // Create the modification summary file
                            var inputFile = new FileInfo(strInputFilePath);
                            var outputFile = new FileInfo(strOutputFilePath);

                            var strModificationSummaryFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, FILENAME_SUFFIX_MOD_SUMMARY));
                            strModificationSummaryFilePath = Path.Combine(outputFile.DirectoryName, strModificationSummaryFilePath);

                            SaveModificationSummaryFile(strModificationSummaryFilePath);
                        }

                        // Inform the user if any errors occurred
                        if (strErrorLog.Length > 0)
                        {
                            SetErrorMessage("Invalid Lines: " + "\n" + strErrorLog);
                        }
                    }

                    blnSuccess = true;
                }
                catch (Exception ex)
                {
                    SetErrorMessage(ex.Message);
                    SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                    blnSuccess = false;
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
                blnSuccess = false;
            }

            return blnSuccess;
        }

        // The following are static to avoid re-reserving space for them for every XTandem results file entry
        Hashtable htSeqsWithMods;
        Hashtable htSeqsWithoutMods;
        private bool ParseXTandemResultsFileEntry(XmlReader objXMLReader, StreamWriter swPeptideResultsFile, ref int intSearchResultCount, ref clsSearchResultsXTandem[] objSearchResults, ref string strErrorLog, int intGroupElementReaderDepth)
        {
            // Note: The number of valid entries in objSearchResults[) is given by intSearchResultCount; objSearchResults(] is expanded but never shrunk
            // There is a separate entry in objSearchResults() for each protein encountered

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

            var strGroupIDInXMLFile = string.Empty;

            var blnProteinSequenceParsed = false;
            var blnDomainParsed = false;
            var blnSuccess = false;

            var strCurrentGroupType = XTANDEM_XML_GROUP_TYPE_MODEL;
            var strCurrentGroupLabel = GROUP_LABEL_PROTEIN;

            try
            {
                // Reset the first entry in objSearchResults
                intSearchResultCount = 0;
                objSearchResults[0].Clear();

                strGroupIDInXMLFile = XMLTextReaderGetAttributeValue(objXMLReader, "id", string.Empty);
                var objSearchResult = objSearchResults[0];
                // Initially set .ResultID to strGroupIDInXMLFile
                // ResultID will get updated to a sequentially assigned number (mNextResultID) if we write the result out to the _xt.txt file
                objSearchResult.ResultID = int.Parse(strGroupIDInXMLFile);
                objSearchResult.GroupID = objSearchResult.ResultID;

                objSearchResult.ParentIonMH = XMLTextReaderGetAttributeValue(objXMLReader, "mh", string.Empty);
                objSearchResult.Charge = XMLTextReaderGetAttributeValue(objXMLReader, "z", string.Empty);

                // Note: This will get updated for each protein encountered
                objSearchResult.PeptideExpectationValue = ConvertEValueToBase10Log(XMLTextReaderGetAttributeValue(objXMLReader, "expect", string.Empty));

                // Note: we truncate .ProteinName at the first space
                objSearchResult.ProteinName = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);
                var intIndex = objSearchResult.ProteinName.IndexOf(' ');
                if (intIndex > 0)
                {
                    objSearchResult.ProteinName = objSearchResult.ProteinName.Substring(0, intIndex);
                }

                objSearchResult.PeptideIntensity = XMLTextReaderGetAttributeValue(objXMLReader, "sumI", string.Empty);
                objSearchResult.PeptideIntensityMax = XMLTextReaderGetAttributeValue(objXMLReader, "maxI", string.Empty);
                objSearchResult.fI = XMLTextReaderGetAttributeValue(objXMLReader, "fI", string.Empty);

                // Continue reading the XML file, loading the information

                while (objXMLReader.Read() & !AbortProcessing)
                {
                    XMLTextReaderSkipWhitespace(objXMLReader);
                    if (objXMLReader.ReadState != ReadState.Interactive)
                        break;

                    if (objXMLReader.NodeType == XmlNodeType.Element)
                    {
                        clsSearchResultsXTandem objResult;
                        switch (objXMLReader.Name.ToLower())
                        {
                            case XTANDEM_XML_ELEMENT_NAME_PROTEIN:
                                intSearchResultCount += 1;
                                if (intSearchResultCount > objSearchResults.Length)
                                {
                                    // Double the length of objSearchResults
                                    Array.Resize(ref objSearchResults, objSearchResults.Length * 2);
                                    for (var intSearchResultIndex = intSearchResultCount - 1; intSearchResultIndex <= objSearchResults.Length - 1; intSearchResultIndex++)
                                    {
                                        objSearchResults[intSearchResultIndex] = new clsSearchResultsXTandem(mPeptideMods, mPeptideSeqMassCalculator);
                                        objSearchResults[intSearchResultIndex].UpdateSearchResultEnzymeAndTerminusInfo(EnzymeMatchSpec, PeptideNTerminusMassChange, PeptideCTerminusMassChange);
                                    }
                                }

                                objResult = objSearchResults[intSearchResultCount - 1];
                                if (intSearchResultCount > 1)
                                {
                                    // Copy the info from objSearchResults[0] to this search result
                                    // ResultID will get updated to a sequentially assigned number (mNextResultID) if we write the result out to the _xt.txt file
                                    objResult.ResultID = objSearchResults[0].ResultID;
                                    objResult.GroupID = objSearchResults[0].GroupID;
                                    objResult.ParentIonMH = objSearchResults[0].ParentIonMH;
                                    objResult.Charge = objSearchResults[0].Charge;

                                    objResult.PeptideExpectationValue = objSearchResults[0].PeptideExpectationValue;

                                    objResult.PeptideIntensity = objSearchResults[0].PeptideIntensity;
                                    objResult.PeptideIntensityMax = objSearchResults[0].PeptideIntensityMax;
                                    objResult.fI = objSearchResults[0].fI;
                                }

                                objResult.ProteinExpectationValue = XMLTextReaderGetAttributeValue(objXMLReader, "expect", string.Empty);
                                objResult.ProteinIntensity = XMLTextReaderGetAttributeValue(objXMLReader, "sumI", string.Empty);

                                // Update the protein name for this protein entry

                                // Note: we truncate .ProteinName at the first space
                                // However, if mLookForReverseSequenceTag = True then we need to look for ":reversed" at the end of the description
                                objResult.ProteinName = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);
                                objResult.ProteinName = TruncateProteinName(objResult.ProteinName);

                                // For proteins with long descriptions, the ":reversed" tag is not present in the label attribute
                                //  and is instead in the <note label="description"> element (a sub-element of the <protein> element
                                // We'll check for this case later in this function

                                // Reset the Protein Sequence Parsed and Domain Parsed flags
                                blnProteinSequenceParsed = false;
                                blnDomainParsed = false;

                                // Clear the protein sequence info and peptide details info
                                objResult.ClearProteinSequenceInfo();
                                objResult.ClearPeptideDetailsInfo();

                                break;
                            case XTANDEM_XML_ELEMENT_NAME_PEPTIDE:
                                if (!blnProteinSequenceParsed)
                                {
                                    blnProteinSequenceParsed = true;
                                    objResult = objSearchResults[intSearchResultCount - 1];
                                    objResult.ProteinSeqResidueNumberStart = XMLTextReaderGetAttributeValue(objXMLReader, "start", 0);
                                    objResult.ProteinSeqResidueNumberEnd = XMLTextReaderGetAttributeValue(objXMLReader, "end", 0);
                                }
                                break;
                            case XTANDEM_XML_ELEMENT_NAME_DOMAIN:
                                // If the given peptide is present in the protein more than once, then domain will appear twice
                                //   Only keep the data for the first domain
                                // Additionally, if X!Tandem decides a given modification could occur on either of two residues, it repeats the domain information
                                //   Again, in this situation, we'll only keep the first domain

                                if (!blnDomainParsed)
                                {
                                    blnDomainParsed = true;

                                    // Cache the XML reader depth before reading any of the element's attributes
                                    var intDomainElementReaderDepth = objXMLReader.Depth;

                                    // Read the information for this domain, storing in objSearchResult
                                    objSearchResult = objSearchResults[intSearchResultCount - 1];
                                    objSearchResult.PeptideLocInProteinStart = XMLTextReaderGetAttributeValue(objXMLReader, "start", 0);
                                    objSearchResult.PeptideLocInProteinEnd = XMLTextReaderGetAttributeValue(objXMLReader, "end", 0);

                                    // Note: the expectation value was already populated from the group level; we'll update it to the value stored for each protein in case it's different for different proteins
                                    objSearchResult.PeptideExpectationValue = ConvertEValueToBase10Log(XMLTextReaderGetAttributeValue(objXMLReader, "expect", string.Empty));

                                    // Note: This is the theoretical, monoisotopic MH+ mass for the peptide
                                    objSearchResult.PeptideMH = XMLTextReaderGetAttributeValue(objXMLReader, "mh", string.Empty);

                                    objSearchResult.PeptideDeltaMass = XMLTextReaderGetAttributeValue(objXMLReader, "delta", string.Empty);

                                    // Note: .peptideDeltaMass is stored in the X!Tandem XML file as "Observed_Mass - Theoretical_Mass"
                                    // However, in MTS .peptideDeltaMass is "Theoretical - Observed"
                                    // Therefore, we will negate .peptideDeltaMass
                                    try
                                    {
                                        objSearchResult.PeptideDeltaMass = (-double.Parse(objSearchResult.PeptideDeltaMass)).ToString();
                                    }
                                    catch (Exception)
                                    {
                                        // Error; Leave .peptideDeltaMass unchanged
                                    }

                                    objSearchResult.PeptideHyperscore = XMLTextReaderGetAttributeValue(objXMLReader, "hyperscore", string.Empty);

                                    // Note that updating .PeptideNextScore will automatically populate .DeltaCn2
                                    objSearchResult.PeptideNextScore = XMLTextReaderGetAttributeValue(objXMLReader, "nextscore", string.Empty);

                                    // Note that calling .PeptidePreResidues, .PeptidePostResidues, and .PeptideCleanSequence will call ComputePeptideCleavageStateInProtein() each time
                                    objSearchResult.PeptidePreResidues = XMLTextReaderGetAttributeValue(objXMLReader, "pre", string.Empty);
                                    objSearchResult.PeptidePostResidues = XMLTextReaderGetAttributeValue(objXMLReader, "post", string.Empty);
                                    objSearchResult.PeptideCleanSequence = XMLTextReaderGetAttributeValue(objXMLReader, "seq", string.Empty);

                                    objSearchResult.PeptideYScore = XMLTextReaderGetAttributeValue(objXMLReader, "y_score", string.Empty);
                                    objSearchResult.PeptideYIons = XMLTextReaderGetAttributeValue(objXMLReader, "y_ions", string.Empty);
                                    objSearchResult.PeptideBScore = XMLTextReaderGetAttributeValue(objXMLReader, "b_score", string.Empty);
                                    objSearchResult.PeptideBIons = XMLTextReaderGetAttributeValue(objXMLReader, "b_ions", string.Empty);

                                    // Now read all of the mods for this domain
                                    // If this is the first search result, then update the mod occurrence counts; otherwise, do not
                                    if (intSearchResultCount == 1)
                                    {
                                        ParseXTandemResultsFileReadDomainMods(objXMLReader, objSearchResults[intSearchResultCount - 1], intDomainElementReaderDepth, true);
                                    }
                                    else
                                    {
                                        ParseXTandemResultsFileReadDomainMods(objXMLReader, objSearchResults[intSearchResultCount - 1], intDomainElementReaderDepth, false);
                                    }
                                }

                                break;
                            case XTANDEM_XML_ELEMENT_NAME_GROUP:
                                var strValue = XMLTextReaderGetAttributeValue(objXMLReader, "type", string.Empty);
                                if (strValue.Length > 0)
                                {
                                    strCurrentGroupType = string.Copy(strValue);
                                    strCurrentGroupLabel = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);
                                }
                                else
                                {
                                    // Leave strGroupType unchanged
                                }

                                break;
                            case XTANDEM_XML_ELEMENT_NAME_NOTE:
                                if (strCurrentGroupType == XTANDEM_XML_GROUP_TYPE_MODEL &&
                                    strCurrentGroupLabel == GROUP_LABEL_PROTEIN && mLookForReverseSequenceTag)
                                {
                                    // Examine the label attribute of this note

                                    strValue = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);

                                    if (strValue == PROTEIN_DESCRIPTION_LABEL)
                                    {
                                        // Check whether this note ends in ":reversed"
                                        // If it does, then make sure objSearchResults[intSearchResultCount - 1].ProteinName ends in :reversed

                                        // Advance the reader before grabbing the inner text
                                        if (objXMLReader.Read())
                                        {
                                            // Read the note's inner text
                                            strValue = XMLTextReaderGetInnerText(objXMLReader);

                                            if (strValue.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR))
                                            {
                                                if (!objSearchResults[intSearchResultCount - 1].ProteinName.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR))
                                                {
                                                    objSearchResults[intSearchResultCount - 1].ProteinName += REVERSED_PROTEIN_SEQUENCE_INDICATOR;
                                                }
                                            }
                                        }
                                    }
                                }
                                else if (strCurrentGroupType == XTANDEM_XML_GROUP_TYPE_SUPPORT &&
                                         strCurrentGroupLabel == GROUP_LABEL_FRAG_ION)
                                {
                                    // This note should contain the scan number
                                    // For _Dta.txt files created at PNNL, it should look something like: "   scan=15118 cs=3"
                                    // For DTA-based files converted to .MGF and processed by X!Tandem: "MyDataset.300.300.2.dta"

                                    // Read the note's inner text
                                    strValue = XMLTextReaderGetInnerText(objXMLReader);
                                    if (strValue != null)
                                    {
                                        var blnScanFound = false;

                                        // Look for the word "scan" followed by an equals sign, followed by a number
                                        // For example, "   scan=15118 cs=3"
                                        try
                                        {
                                            var match = mScanNumberRegExA.Match(strValue);
                                            if (match.Success && match.Groups.Count > 1)
                                            {
                                                objSearchResults[0].Scan = match.Groups[1].Value;
                                                blnScanFound = true;
                                            }
                                        }
                                        catch (Exception)
                                        {
                                            // Ignore errors here
                                        }

                                        if (!blnScanFound)
                                        {
                                            // No match; look for the word "scan" followed by whitespace, followed by a number
                                            // For example, "scan 300"
                                            try
                                            {
                                                var match = mScanNumberRegExB.Match(strValue);
                                                if (match.Success && match.Groups.Count > 1)
                                                {
                                                    objSearchResults[0].Scan = match.Groups[1].Value;
                                                    blnScanFound = true;
                                                }
                                            }
                                            catch (Exception)
                                            {
                                                // Ignore errors here
                                            }
                                        }

                                        if (!blnScanFound)
                                        {
                                            // No match; see if the description resembles a .Dta file name
                                            // For example, "MyDataset.300.300.2.dta"
                                            try
                                            {
                                                var match = mScanNumberRegExC.Match(strValue);
                                                if (match.Success && match.Groups.Count > 1)
                                                {
                                                    objSearchResults[0].Scan = match.Groups[1].Value;
                                                    blnScanFound = true;
                                                }
                                            }
                                            catch (Exception)
                                            {
                                                // Ignore errors here
                                            }
                                        }

                                        if (!blnScanFound)
                                        {
                                            // Still no match; extract out the first number present
                                            try
                                            {
                                                var match = mScanNumberRegExD.Match(strValue);
                                                if (match.Success)
                                                {
                                                    objSearchResults[0].Scan = match.Value;
                                                }
                                                else
                                                {
                                                    objSearchResults[0].Scan = strValue;
                                                }
                                            }
                                            catch (Exception)
                                            {
                                                // Ignore errors here
                                            }
                                        }

                                        // Copy the scan value from the first result to the other results
                                        for (var intSearchResultIndex = 1; intSearchResultIndex <= intSearchResultCount - 1; intSearchResultIndex++)
                                        {
                                            objSearchResults[intSearchResultIndex].Scan = objSearchResults[0].Scan;
                                        }
                                    }
                                }
                                break;
                            default:
                                // Unknown/unneeded child node name; ignore it
                                break;
                        }
                    }
                    else if (objXMLReader.NodeType == XmlNodeType.EndElement)
                    {
                        if (objXMLReader.Name == XTANDEM_XML_ELEMENT_NAME_GROUP)
                        {
                            if (objXMLReader.Depth <= intGroupElementReaderDepth)
                            {
                                // End element found for the current group

                                // Typically each group will consist of entries all having the same sequence and modifications (but different protein names)
                                // However, occasionally a group will contain a mix of peptide sequences (only occurs if they all had the exact same hyperscore)
                                // In order to check for this, we will construct a pointer array of Sequence and Mods to SearchResultIndex and use this to determine
                                //  which entries should be written to the _xt.txt file and to the ResultToSeqMap file
                                // We will also use this pointer array to keep track of the number of proteins listed for each peptide

                                if (htSeqsWithMods == null)
                                {
                                    htSeqsWithMods = new Hashtable();
                                    htSeqsWithoutMods = new Hashtable();
                                }
                                else
                                {
                                    htSeqsWithMods.Clear();
                                    htSeqsWithoutMods.Clear();
                                }

                                // First step through the results to compute the mass, construct the modification description,
                                //  and determine the number of proteins listed for each
                                string strSequenceWithMods;
                                for (var intSearchResultIndex = 0; intSearchResultIndex <= intSearchResultCount - 1; intSearchResultIndex++)
                                {
                                    bool blnUpdateModOccurrenceCounts;
                                    if (intSearchResultIndex == 0)
                                    {
                                        // Always set blnUpdateModOccurrenceCounts to True for the first result in the group
                                        blnUpdateModOccurrenceCounts = true;
                                        htSeqsWithoutMods.Add(objSearchResults[intSearchResultIndex].PeptideCleanSequence, 1);
                                    }
                                    else
                                    {
                                        if (htSeqsWithoutMods.ContainsKey(objSearchResults[intSearchResultIndex].PeptideCleanSequence))
                                        {
                                            blnUpdateModOccurrenceCounts = false;
                                        }
                                        else
                                        {
                                            blnUpdateModOccurrenceCounts = true;
                                            htSeqsWithoutMods.Add(objSearchResults[intSearchResultIndex].PeptideCleanSequence, 1);
                                        }
                                    }

                                    blnSuccess = AddModificationsAndComputeMass(objSearchResults[intSearchResultIndex], blnUpdateModOccurrenceCounts);
                                    if (!blnSuccess)
                                    {
                                        if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                                        {
                                            strErrorLog += "Error adding modifications to sequence for Group ID '" + strGroupIDInXMLFile + "'" + "\n";
                                        }
                                    }

                                    strSequenceWithMods = objSearchResults[intSearchResultIndex].PeptideCleanSequence + "_" + objSearchResults[intSearchResultIndex].PeptideModDescription;

                                    if (intSearchResultIndex == 0)
                                    {
                                        // Always add the first result in the group htSeqsWithMods
                                        htSeqsWithMods.Add(strSequenceWithMods, 1);
                                    }
                                    else
                                    {
                                        // See if htSeqsWithMods contains strSequenceWithMods
                                        if (htSeqsWithMods.ContainsKey(strSequenceWithMods))
                                        {
                                            // Increment the protein count for this peptide
                                            htSeqsWithMods[strSequenceWithMods] = (int)htSeqsWithMods[strSequenceWithMods] + 1;
                                        }
                                        else
                                        {
                                            htSeqsWithMods.Add(strSequenceWithMods, 1);
                                        }
                                    }
                                }

                                // Now step through the list again and update the MultipleProteinCount value for each search result
                                for (var intSearchResultIndex = 0; intSearchResultIndex <= intSearchResultCount - 1; intSearchResultIndex++)
                                {
                                    strSequenceWithMods = objSearchResults[intSearchResultIndex].PeptideCleanSequence + "_" + objSearchResults[intSearchResultIndex].PeptideModDescription;

                                    int intProteinCount;
                                    try
                                    {
                                        intProteinCount = (int)htSeqsWithMods[strSequenceWithMods];
                                    }
                                    catch (Exception)
                                    {
                                        intProteinCount = 1;
                                    }

                                    if (intProteinCount < 1)
                                        intProteinCount = 1;

                                    // Note: Multiple protein count is 0 if the peptide is only in 1 protein; 1 if the protein is in 2 proteins, etc.
                                    objSearchResults[intSearchResultIndex].MultipleProteinCount = (intProteinCount - 1).ToString();
                                }

                                // Clear htSeqsWithMods again since we need to re-use it to determine which results to write out
                                htSeqsWithMods.Clear();

                                // Write out the results
                                for (var intSearchResultIndex = 0; intSearchResultIndex <= intSearchResultCount - 1; intSearchResultIndex++)
                                {
                                    strSequenceWithMods = objSearchResults[intSearchResultIndex].PeptideCleanSequence + "_" + objSearchResults[intSearchResultIndex].PeptideModDescription;

                                    bool blnUpdateResultToSeqMapFile;
                                    if (intSearchResultIndex == 0)
                                    {
                                        // Always save the first result in the group to the _xt.txt and _ResultToSeqMap.txt files
                                        htSeqsWithMods.Add(strSequenceWithMods, 1);
                                        blnUpdateResultToSeqMapFile = true;
                                    }
                                    else
                                    {
                                        // See if htSeqsWithMods contains strSequenceWithMods
                                        if (htSeqsWithMods.ContainsKey(strSequenceWithMods))
                                        {
                                            blnUpdateResultToSeqMapFile = false;
                                        }
                                        else
                                        {
                                            htSeqsWithMods.Add(strSequenceWithMods, 1);
                                            blnUpdateResultToSeqMapFile = true;
                                        }
                                    }

                                    if (blnUpdateResultToSeqMapFile)
                                    {
                                        // Only save the first result for each peptide in the group to the _xt.txt and _ResultToSeqMap.txt files
                                        // Note: This function will update .ResultID to the next available ID value (mNextResultID)
                                        SaveXTandemResultsFileEntry(objSearchResults[intSearchResultIndex], ref swPeptideResultsFile);
                                    }

                                    SaveResultsFileEntrySeqInfo(objSearchResults[intSearchResultIndex], blnUpdateResultToSeqMapFile);
                                }

                                blnSuccess = true;
                                break;
                            }
                        }
                    }
                }
            }
            catch (Exception)
            {
                // Error parsing values from this group ID in the XML file
                if (strErrorLog.Length < MAX_ERROR_LOG_LENGTH)
                {
                    strErrorLog += "Error parsing value for Group ID '" + strGroupIDInXMLFile + "'" + "\n";
                }
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private bool ParseXTandemResultsFileInputParameters(string strInputFilePath)
        {
            // Pre-read the XML file and look for the Input Parameters section
            // Read the parameters and validate that each of the mods defined is present in mPeptideMods

            const string GROUP_LABEL_INPUT_PARAMETERS = "input parameters";

            bool blnSuccess;

            try
            {
                // Open the input file and parse it
                // Initialize the stream reader and the XML Text Reader
                eCurrentXMLDataFileSectionConstants eCurrentXMLDataFileSection;
                using (var objXMLReader = new XmlTextReader(strInputFilePath))
                {
                    // Parse the file
                    eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.UnknownFile;

                    while (objXMLReader.Read() & !AbortProcessing)
                    {
                        XMLTextReaderSkipWhitespace(objXMLReader);
                        if (objXMLReader.ReadState != ReadState.Interactive)
                            break;

                        if (objXMLReader.Depth < 2)
                        {
                            if (objXMLReader.NodeType == XmlNodeType.Element)
                            {
                                switch (objXMLReader.Name.ToLower())
                                {
                                    case XTANDEM_XML_ELEMENT_NAME_GROUP:
                                        if (objXMLReader.HasAttributes)
                                        {
                                            var intParametersGroupDepth = objXMLReader.Depth;

                                            // See if the group has a "type" attribute containing the text XTANDEM_XML_GROUP_TYPE_PARAMETERS
                                            var strCurrentGroupType = XMLTextReaderGetAttributeValue(objXMLReader, "type", string.Empty);
                                            if (strCurrentGroupType == XTANDEM_XML_GROUP_TYPE_PARAMETERS)
                                            {
                                                eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.InputParameters;

                                                // Read the Label for this group
                                                var strCurrentGroupLabel = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);
                                                if (strCurrentGroupLabel == GROUP_LABEL_INPUT_PARAMETERS)
                                                {
                                                    // Read the input parameters
                                                    ParseXTandemResultsFileInputParametersWork(objXMLReader, intParametersGroupDepth);
                                                }
                                            }
                                            else
                                            {
                                                // Skip this group
                                                objXMLReader.Skip();
                                            }
                                        }
                                        else
                                        {
                                            // Group doesn't have any attributes; ignore it
                                            objXMLReader.Skip();
                                        }

                                        break;
                                    case XTANDEM_XML_ROOT_ELEMENT:
                                        eCurrentXMLDataFileSection = eCurrentXMLDataFileSectionConstants.Start;
                                        break;
                                    default:
                                        // Skip this element
                                        objXMLReader.Skip();
                                        break;
                                }
                            }
                        }
                    }
                }

                if (eCurrentXMLDataFileSection == eCurrentXMLDataFileSectionConstants.UnknownFile)
                {
                    SetErrorMessage("Root element '" + XTANDEM_XML_ROOT_ELEMENT + "' not found in the input file: " + strInputFilePath);
                    blnSuccess = false;
                }
                else
                {
                    blnSuccess = true;
                }
            }
            catch (Exception ex)
            {
                SetErrorMessage(ex.Message);
                SetErrorCode(ePHRPErrorCodes.ErrorReadingInputFile);
                blnSuccess = false;
            }

            return blnSuccess;
        }

        private void ParseXTandemResultsFileInputParametersWork(XmlReader objXMLReader, int intParametersGroupDepth)
        {
            // Read the input parameters
            // Each parameter is an element with name "note" with attributes "type" and "label"
            // The parameter value is the text between within the element

            const string NOTE_TYPE_INPUT = "input";
            const char CLEAVAGE_SPEC_SEP = '|';
            const char XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START = '{';
            const char XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END = '}';

            // Initialize the Mod Info array
            var intModInfoCount = 0;
            var udtModInfo = new udtSearchOptionModificationInfoType[20];

            var blnStaticModsAreResetForRefinement = true;


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
            for (var intIndex = 0; intIndex <= udtParamLabels.Length - 1; intIndex++)
            {
                udtParamLabels[intIndex] = udtParamLabels[intIndex].ToLower();
            }

            while (objXMLReader.Read())
            {
                XMLTextReaderSkipWhitespace(objXMLReader);
                if (objXMLReader.ReadState != ReadState.Interactive)
                    break;

                if (objXMLReader.NodeType == XmlNodeType.Element)
                {
                    switch (objXMLReader.Name.ToLower())
                    {
                        case XTANDEM_XML_ELEMENT_NAME_NOTE:
                            // Read the note's type
                            var strNoteType = XMLTextReaderGetAttributeValue(objXMLReader, "type", string.Empty);

                            if (strNoteType == NOTE_TYPE_INPUT)
                            {
                                // Read the note's label and inner text
                                var strNoteLabel = XMLTextReaderGetAttributeValue(objXMLReader, "label", string.Empty);

                                // Need to advance the reader before calling XMLTextReaderGetInnerText
                                objXMLReader.Read();
                                var strValue = XMLTextReaderGetInnerText(objXMLReader);

                                if (strValue != null)
                                {
                                    var strNoteLabelLower = strNoteLabel.ToLower();
                                    if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Residue_StaticModMass]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.StaticMod, Convert.ToInt32(eInputParamLabelNames.Residue_StaticModMass), false, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Residue_PotentialModMass]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Residue_PotentialModMass), false, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Residue_PotentialModMotif]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Residue_PotentialModMotif), true, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialModMass]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialModMass), false, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialModMotif]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialModMotif), true, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialNTerminusMods]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialNTerminusMods), false, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_PotentialCTerminusMods]))
                                    {
                                        ParseXTandemInputParameterModInfo(clsModificationDefinition.eModificationTypeConstants.DynamicMod, Convert.ToInt32(eInputParamLabelNames.Refine_PotentialCTerminusMods), false, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_NTerminal_ResidueModMass]))
                                    {
                                        ParseXTandemInputParameterProteinTerminusMod(Convert.ToInt32(eInputParamLabelNames.Protein_NTerminal_ResidueModMass), true, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_CTerminal_ResidueModMass]))
                                    {
                                        ParseXTandemInputParameterProteinTerminusMod(Convert.ToInt32(eInputParamLabelNames.Protein_CTerminal_ResidueModMass), false, strValue, ref intModInfoCount, ref udtModInfo);
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_NTerminalMassChange]))
                                    {
                                        if (clsPHRPParser.IsNumber(strValue))
                                        {
                                            PeptideNTerminusMassChange = double.Parse(strValue);
                                        }
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_CTerminalMassChange]))
                                    {
                                        if (clsPHRPParser.IsNumber(strValue))
                                        {
                                            PeptideCTerminusMassChange = double.Parse(strValue);
                                        }
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Protein_Cleavage_Site]))
                                    {
                                        // In X!Tandem the LeftSpec and RightSpec values are separated by a vertical bar (CLEAVAGE_SPEC_SEP)
                                        // Look for CLEAVAGE_SPEC_SEP in strValue
                                        var intBarLoc = strValue.IndexOf(CLEAVAGE_SPEC_SEP);
                                        if (intBarLoc > 0 & intBarLoc < strValue.Length - 1)
                                        {
                                            var strLeftSpec = strValue.Substring(0, intBarLoc);
                                            var strRightSpec = strValue.Substring(intBarLoc + 1);

                                            // Look for curly braces in strLeftSpec and strRightSpec
                                            // In X!Tandem curly braces mean to not match a residue
                                            // If found, change to standard RegEx notation, e.g. from {P} to [^P]
                                            if (strLeftSpec.IndexOf(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START) >= 0)
                                            {
                                                strLeftSpec = strLeftSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START.ToString(), "[^");
                                                strLeftSpec = strLeftSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END.ToString(), "]");
                                            }

                                            if (strRightSpec.IndexOf(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START) >= 0)
                                            {
                                                strRightSpec = strRightSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_START.ToString(), "[^");
                                                strRightSpec = strRightSpec.Replace(XTANDEM_CLEAVAGE_NEGATION_SYMBOL_END.ToString(), "]");
                                            }

                                            EnzymeMatchSpec = new clsPeptideCleavageStateCalculator.udtEnzymeMatchSpecType(strLeftSpec, strRightSpec);
                                        }
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Refine_ModificationMass]))
                                    {
                                        if (strValue != null && strValue.Trim().Length > 0)
                                        {
                                            blnStaticModsAreResetForRefinement = true;
                                        }
                                    }
                                    else if (strNoteLabelLower.Equals(udtParamLabels[(int)eInputParamLabelNames.Scoring_Include_Reverse]))
                                    {
                                        if (strValue != null && strValue.Trim().Length > 0)
                                        {
                                            if (strValue.Trim().ToLower() == "yes")
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
                else if (objXMLReader.NodeType == XmlNodeType.EndElement)
                {
                    if (objXMLReader.Depth == intParametersGroupDepth && objXMLReader.Name == XTANDEM_XML_ELEMENT_NAME_GROUP)
                    {
                        // Reached the end of this group
                        break;
                    }
                }
            }

            if (intModInfoCount > 0)
            {
                // Validate that each of the mods in udtModInfo is present in mPeptideMods
                if (intModInfoCount > 1)
                {
                    // Sort udtModInfo
                    Array.Sort(udtModInfo, 0, intModInfoCount, new ISearchOptionModificationInfoComparer());
                }

                // Before continuing, look for Static residue mods in udtModInfo
                // If any are found, and if an identical dynamic residue mod is already present, then delete the static residue mod
                // Additionally, if 	<note type="input" label="refine, modification mass">none</note> was prsent in the XTandem results file then
                //  auto update all static mods to dynamic mods since they are reset during refinement
                for (var intIndex = 0; intIndex < intModInfoCount;)
                {
                    var blnModDeleted = false;
                    if (udtModInfo[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod)
                    {
                        var intIndexCompare = 0;
                        while (intIndexCompare < intModInfoCount)
                        {
                            if (intIndexCompare != intIndex && udtModInfo[intIndexCompare].ModificationType == clsModificationDefinition.eModificationTypeConstants.DynamicMod)
                            {
                                // See if this modification has a similar mass (within MASS_DIGITS_OF_PRECISION digits of precision)
                                if (Math.Abs(Math.Round(Math.Abs(udtModInfo[intIndexCompare].ModificationMass - udtModInfo[intIndex].ModificationMass), clsPeptideModificationContainer.MASS_DIGITS_OF_PRECISION)) < float.Epsilon)
                                {
                                    // Matching mass
                                    // Compare .TargetResidues
                                    if (clsModificationDefinition.EquivalentTargetResidues(udtModInfo[intIndexCompare].TargetResidues, udtModInfo[intIndex].TargetResidues, true))
                                    {
                                        // Yes, the modifications match; delete the static version of the modification
                                        for (var intIndexCopy = intIndex; intIndexCopy <= intModInfoCount - 2; intIndexCopy++)
                                        {
                                            udtModInfo[intIndexCopy] = udtModInfo[intIndexCopy + 1];
                                        }
                                        intModInfoCount -= 1;
                                        blnModDeleted = true;
                                        break;
                                    }
                                }
                            }
                            intIndexCompare += 1;
                        }
                    }

                    if (!blnModDeleted)
                    {
                        if (udtModInfo[intIndex].ModificationType == clsModificationDefinition.eModificationTypeConstants.StaticMod && blnStaticModsAreResetForRefinement)
                        {
                            // Update this static mod to be a dynamic mod
                            udtModInfo[intIndex].ModificationType = clsModificationDefinition.eModificationTypeConstants.DynamicMod;
                        }
                        intIndex += 1;
                    }
                }

                for (var intIndex = 0; intIndex <= intModInfoCount - 1; intIndex++)
                {
                    var modInfo = udtModInfo[intIndex];
                    mPeptideMods.VerifyModificationPresent(modInfo.ModificationMass, modInfo.TargetResidues, modInfo.ModificationType);
                }
            }

            // In addition, verify that the standard refinement modifications are present in mPeptideMods
            mPeptideMods.AppendStandardRefinmentModifications();
        }

        private void ParseXTandemResultsFileReadDomainMods(XmlReader objXMLReader, clsSearchResultsBaseClass objSearchResult, int intDomainElementReaderDepth, bool blnUpdateModOccurrenceCounts)
        {
            // Continue reading the XML file, loading the information

            while (objXMLReader.Read() && objXMLReader.Depth >= intDomainElementReaderDepth)
            {
                if (objXMLReader.ReadState != ReadState.Interactive)
                    break;

                if (objXMLReader.NodeType == XmlNodeType.Element)
                {
                    switch (objXMLReader.Name.ToLower())
                    {
                        case XTANDEM_XML_ELEMENT_NAME_AMINOACID:
                            var strValue = XMLTextReaderGetAttributeValue(objXMLReader, "type", "").Trim();
                            char chTargetResidue;
                            if (string.IsNullOrWhiteSpace(strValue))
                            {
                                chTargetResidue = default(char);
                            }
                            else
                            {
                                chTargetResidue = strValue[0];
                            }

                            var intModifiedResiduePosInProtein = XMLTextReaderGetAttributeValue(objXMLReader, "at", 0);

                            if (intModifiedResiduePosInProtein > 0)
                            {
                                var dblModificationMass = XMLTextReaderGetAttributeValueDbl(objXMLReader, "modified", 0);

                                if (Math.Abs(dblModificationMass - 0) > float.Epsilon)
                                {
                                    var intResidueLocInPeptide = intModifiedResiduePosInProtein - objSearchResult.PeptideLocInProteinStart + 1;
                                    var eResidueTerminusState = objSearchResult.DetermineResidueTerminusState(intResidueLocInPeptide);

                                    objSearchResult.SearchResultAddModification(dblModificationMass,
                                                                                chTargetResidue,
                                                                                intResidueLocInPeptide,
                                                                                eResidueTerminusState,
                                                                                blnUpdateModOccurrenceCounts);
                                }
                            }
                            break;
                    }
                }
                else if (objXMLReader.NodeType == XmlNodeType.EndElement)
                {
                    if (objXMLReader.Name == XTANDEM_XML_ELEMENT_NAME_DOMAIN)
                    {
                        break;
                    }
                }
            }
        }

        /// <summary>
        /// Main processing function
        /// </summary>
        /// <param name="strInputFilePath">X!Tandem results file</param>
        /// <param name="strOutputFolderPath">Output folder</param>
        /// <param name="strParameterFilePath">Parameter file</param>
        /// <returns>True if success, False if failure</returns>
        public override bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            var blnSuccess = false;

            if (!LoadParameterFileSettings(strParameterFilePath))
            {
                SetErrorCode(ePHRPErrorCodes.ErrorReadingParameterFile, true);
                return false;
            }

            try
            {
                if (string.IsNullOrWhiteSpace(strInputFilePath))
                {
                    SetErrorMessage("Input file name is empty");
                    SetErrorCode(ePHRPErrorCodes.InvalidInputFilePath);
                    return false;
                }

                blnSuccess = ResetMassCorrectionTagsAndModificationDefinitions();
                if (!blnSuccess)
                {
                    return false;
                }

                ResetProgress("Parsing " + Path.GetFileName(strInputFilePath));

                if (!CleanupFilePaths(ref strInputFilePath, ref strOutputFolderPath))
                {
                    return false;
                }

                try
                {
                    // Obtain the full path to the input file
                    var inputFile = new FileInfo(strInputFilePath);

                    // Define the output file name based on strInputFilePath
                    // The name will be DatasetName_xt.txt

                    var strXtandemXTFilePath = Path.GetFileName(ReplaceFilenameSuffix(inputFile, ".txt"));
                    strXtandemXTFilePath = Path.Combine(strOutputFolderPath, strXtandemXTFilePath);
                    blnSuccess = ParseXTandemResultsFile(inputFile.FullName, strXtandemXTFilePath, false);

                    if (blnSuccess && CreateProteinModsFile)
                    {
                        blnSuccess = CreateProteinModsFileWork(inputFile, strOutputFolderPath, strXtandemXTFilePath);
                    }

                    if (blnSuccess)
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

            return blnSuccess;
        }

        private bool CreateProteinModsFileWork(FileInfo inputFile, string strOutputFolderPath, string strXtandemXTFilePath)
        {
            bool blnSuccess;

            // First create the MTS PepToProteinMap file using inputFile
            var lstSourcePHRPDataFiles = new List<string> {
                strXtandemXTFilePath
            };

            var strMTSPepToProteinMapFilePath = ConstructPepToProteinMapFilePath(inputFile.FullName, strOutputFolderPath, MTS: true);

            if (File.Exists(strMTSPepToProteinMapFilePath) && UseExistingMTSPepToProteinMapFile)
            {
                blnSuccess = true;
            }
            else
            {
                blnSuccess = CreatePepToProteinMapFile(lstSourcePHRPDataFiles, strMTSPepToProteinMapFilePath);
                if (!blnSuccess)
                {
                    ReportWarning("Skipping creation of the ProteinMods file since CreatePepToProteinMapFile returned False");
                }
            }

            if (blnSuccess)
            {
                // If necessary, copy various PHRPReader support files (in particular, the MSGF file) to the output folder
                ValidatePHRPReaderSupportFiles(Path.Combine(inputFile.DirectoryName, Path.GetFileName(strXtandemXTFilePath)), strOutputFolderPath);

                // Now create the Protein Mods file
                blnSuccess = CreateProteinModDetailsFile(strXtandemXTFilePath, strOutputFolderPath, strMTSPepToProteinMapFilePath, clsPHRPReader.ePeptideHitResultType.XTandem);
            }

            if (!blnSuccess)
            {
                // Do not treat this as a fatal error
                return true;
            }
            return true;
        }

        private void SaveXTandemResultsFileEntry(clsSearchResultsXTandem objSearchResult, ref StreamWriter swPeptideResultsFile)
        {
            // Update .ResultID to the next available number
            objSearchResult.ResultID = mNextResultID;
            mNextResultID += 1;

            // Write the results to the output file
            swPeptideResultsFile.WriteLine(objSearchResult.ResultID + SEP_CHAR +
                                           objSearchResult.GroupID + SEP_CHAR +
                                           objSearchResult.Scan + SEP_CHAR +
                                           objSearchResult.Charge + SEP_CHAR +
                                           objSearchResult.PeptideMH + SEP_CHAR +
                                           objSearchResult.PeptideHyperscore + SEP_CHAR +
                                           objSearchResult.PeptideExpectationValue + SEP_CHAR +
                                           objSearchResult.MultipleProteinCount + SEP_CHAR +
                                           objSearchResult.SequenceWithPrefixAndSuffix(true) + SEP_CHAR +
                                           Math.Round(objSearchResult.PeptideDeltaCn2, 4).ToString(CultureInfo.InvariantCulture) + SEP_CHAR +
                                           objSearchResult.PeptideYScore + SEP_CHAR +
                                           objSearchResult.PeptideYIons + SEP_CHAR +
                                           objSearchResult.PeptideBScore + SEP_CHAR +
                                           objSearchResult.PeptideBIons + SEP_CHAR +
                                           objSearchResult.PeptideDeltaMass + SEP_CHAR +
                                           objSearchResult.PeptideIntensity + SEP_CHAR +
                                           PRISM.StringUtilities.DblToString(objSearchResult.PeptideDeltaMassCorrectedPpm, 5, 0.00005));
        }

        protected override string TruncateProteinName(string strProteinNameAndDescription)
        {
            var blnIsReversed = false;

            if (mLookForReverseSequenceTag)
            {
                blnIsReversed = strProteinNameAndDescription.EndsWith(REVERSED_PROTEIN_SEQUENCE_INDICATOR);
            }

            strProteinNameAndDescription = base.TruncateProteinName(strProteinNameAndDescription);

            if (blnIsReversed)
            {
                return strProteinNameAndDescription + REVERSED_PROTEIN_SEQUENCE_INDICATOR;
            }
            return strProteinNameAndDescription;
        }

        private string XMLTextReaderGetAttributeValue(XmlReader objXMLReader, string strAttributeName, string strValueIfMissing)
        {
            objXMLReader.MoveToAttribute(strAttributeName);
            if (objXMLReader.ReadAttributeValue())
            {
                return objXMLReader.Value;
            }
            return string.Copy(strValueIfMissing);
        }

        private int XMLTextReaderGetAttributeValue(XmlReader objXMLReader, string strAttributeName, int intValueIfMissing)
        {
            objXMLReader.MoveToAttribute(strAttributeName);
            if (objXMLReader.ReadAttributeValue())
            {
                if (clsPHRPParser.IsNumber(objXMLReader.Value))
                {
                    return Convert.ToInt32(objXMLReader.Value);
                }
                return intValueIfMissing;
            }
            return intValueIfMissing;
        }

        private double XMLTextReaderGetAttributeValueDbl(XmlReader objXMLReader, string strAttributeName, double dblValueIfMissing)
        {
            objXMLReader.MoveToAttribute(strAttributeName);
            if (objXMLReader.ReadAttributeValue())
            {
                if (clsPHRPParser.IsNumber(objXMLReader.Value))
                {
                    return Convert.ToDouble(objXMLReader.Value);
                }
                return dblValueIfMissing;
            }
            return dblValueIfMissing;
        }

        private string XMLTextReaderGetInnerText(XmlReader objXMLReader)
        {
            var strValue = string.Empty;
            bool blnSuccess;

            if (objXMLReader.NodeType == XmlNodeType.Element)
            {
                // Advance the reader so that we can read the value
                blnSuccess = objXMLReader.Read();
            }
            else
            {
                blnSuccess = true;
            }

            if (blnSuccess && objXMLReader.NodeType != XmlNodeType.Whitespace & objXMLReader.HasValue)
            {
                strValue = objXMLReader.Value;
            }

            return strValue;
        }

        private void XMLTextReaderSkipWhitespace(XmlReader objXMLReader)
        {
            if (objXMLReader.NodeType == XmlNodeType.Whitespace)
            {
                // Whitespace; read the next node
                objXMLReader.Read();
            }
        }
    }
}
