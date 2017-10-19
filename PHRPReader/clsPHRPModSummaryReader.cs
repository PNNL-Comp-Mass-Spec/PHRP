using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader
{
    /// <summary>
    /// ModSummary file reader
    /// </summary>
    public class clsPHRPModSummaryReader
    {

#pragma warning disable 1591
        public const string MOD_SUMMARY_COLUMN_Modification_Symbol = "Modification_Symbol";
        public const string MOD_SUMMARY_COLUMN_Modification_Mass = "Modification_Mass";
        public const string MOD_SUMMARY_COLUMN_Target_Residues = "Target_Residues";
        public const string MOD_SUMMARY_COLUMN_Modification_Type = "Modification_Type";
        public const string MOD_SUMMARY_COLUMN_Mass_Correction_Tag = "Mass_Correction_Tag";
        public const string MOD_SUMMARY_COLUMN_Occurrence_Count = "Occurrence_Count";
#pragma warning restore 1591

        private readonly List<clsModificationDefinition> mModificationDefs;

        // The keys in this dictionary are MassCorrectionTag names and the values are the modification mass, stored as text (as it appears in the _ModSummary file)
        private readonly Dictionary<string, string> mModDefMassesAsText;

        private readonly bool mSuccess;

        public List<clsModificationDefinition> ModificationDefs => mModificationDefs;

        // ReSharper disable once ConvertToAutoProperty
        public bool Success => mSuccess;

        public clsPHRPModSummaryReader(string strModSummaryFilePath)
        {
            mModificationDefs = new List<clsModificationDefinition>();
            mModDefMassesAsText = new Dictionary<string, string>();

            mSuccess = false;

            if (string.IsNullOrEmpty(strModSummaryFilePath))
            {
                throw new Exception("ModSummaryFilePath is empty; unable to continue");
            }
            else if (!File.Exists(strModSummaryFilePath))
            {
                throw new FileNotFoundException("ModSummary file not found: " + strModSummaryFilePath);
            }

            mSuccess = ReadModSummaryFile(strModSummaryFilePath, ref mModificationDefs);
        }

        /// <summary>
        /// Returns the mass value associated with the given mass correction tag
        /// </summary>
        /// <param name="strMassCorrectionTag"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public string GetModificationMassAsText(string strMassCorrectionTag)
        {
            var strModMass = string.Empty;

            if (mModDefMassesAsText.TryGetValue(strMassCorrectionTag, out strModMass))
            {
                return strModMass;
            }
            else
            {
                return string.Empty;
            }
        }

        private bool ReadModSummaryFile(string strModSummaryFilePath, ref List<clsModificationDefinition> lstModInfo)
        {
            string strLineIn = null;
            string[] strSplitLine = null;

            var objColumnHeaders = default(SortedDictionary<string, int>);

            string strModSymbol = null;
            string strModMass = null;
            string strTargetResidues = null;
            string strModType = null;
            string strMassCorrectionTag = null;

            var chModSymbol = default(char);
            double dblModificationMass = 0;
            var eModificationType = default(clsModificationDefinition.eModificationTypeConstants);

            var blnSkipLine = false;
            var blnHeaderLineParsed = false;

            if (lstModInfo == null)
            {
                lstModInfo = new List<clsModificationDefinition>();
            }
            else
            {
                lstModInfo.Clear();
            }

            if (string.IsNullOrEmpty(strModSummaryFilePath))
            {
                return false;
            }

            // Initialize the column mapping
            // Using a case-insensitive comparer
            objColumnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase);

            // Define the default column mapping
            objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Symbol, 0);
            objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Mass, 1);
            objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Target_Residues, 2);
            objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Modification_Type, 3);
            objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Mass_Correction_Tag, 4);
            objColumnHeaders.Add(MOD_SUMMARY_COLUMN_Occurrence_Count, 5);

            // Read the data from the ModSummary.txt file
            // The first line is typically a header line:
            // Modification_Symbol  Modification_Mass  Target_Residues  Modification_Type  Mass_Correction_Tag  Occurrence_Count

            using (var srModSummaryFile = new StreamReader(new FileStream(strModSummaryFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                blnHeaderLineParsed = false;

                while (!srModSummaryFile.EndOfStream)
                {
                    strLineIn = srModSummaryFile.ReadLine();
                    blnSkipLine = false;

                    if (!string.IsNullOrEmpty(strLineIn))
                    {
                        strSplitLine = strLineIn.Split('\t');

                        if (!blnHeaderLineParsed)
                        {
                            if (strSplitLine[0].ToLower() == MOD_SUMMARY_COLUMN_Modification_Symbol.ToLower())
                            {
                                // Parse the header line to confirm the column ordering
                                // The Occurrence_Count column was mispelled prior to December 2012; need to check for this
                                for (var intIndex = 0; intIndex <= strSplitLine.Length - 1; intIndex++)
                                {
                                    if (strSplitLine[intIndex] == "Occurence_Count")
                                        strSplitLine[intIndex] = MOD_SUMMARY_COLUMN_Occurrence_Count;
                                }
                                clsPHRPReader.ParseColumnHeaders(strSplitLine, objColumnHeaders);
                                blnSkipLine = true;
                            }

                            blnHeaderLineParsed = true;
                        }

                        if (!blnSkipLine && strSplitLine.Length >= 4)
                        {
                            strModSymbol = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Symbol, objColumnHeaders);
                            strModMass = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Mass, objColumnHeaders);
                            strTargetResidues = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Target_Residues, objColumnHeaders);
                            strModType = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Modification_Type, objColumnHeaders);
                            strMassCorrectionTag = clsPHRPReader.LookupColumnValue(strSplitLine, MOD_SUMMARY_COLUMN_Mass_Correction_Tag, objColumnHeaders);

                            if (string.IsNullOrWhiteSpace(strModSymbol))
                            {
                                strModSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL.ToString();
                            }
                            chModSymbol = strModSymbol[0];

                            if (!double.TryParse(strModMass, out dblModificationMass))
                            {
                                throw new Exception("Modification mass is not numeric for MassCorrectionTag: " + strMassCorrectionTag + ": " + strModMass);
                            }

                            if (string.IsNullOrWhiteSpace(strModType))
                            {
                                eModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType;
                            }
                            else
                            {
                                eModificationType = clsModificationDefinition.ModificationSymbolToModificationType(strModType[0]);
                            }

                            var objModDef = default(clsModificationDefinition);
                            objModDef = new clsModificationDefinition(chModSymbol, dblModificationMass, strTargetResidues, eModificationType, strMassCorrectionTag);
                            objModDef.ModificationMassAsText = strModMass;

                            lstModInfo.Add(objModDef);

                            if (!mModDefMassesAsText.ContainsKey(strMassCorrectionTag))
                            {
                                mModDefMassesAsText.Add(strMassCorrectionTag, strModMass);
                            }
                        }
                    }
                }
            }

            return true;
        }
    }
}
