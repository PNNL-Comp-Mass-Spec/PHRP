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

        /// <summary>
        /// Modification list
        /// </summary>
        public List<clsModificationDefinition> ModificationDefs => mModificationDefs;

        /// <summary>
        /// True if the mod summary was successfully loaded
        /// </summary>
        public bool Success { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="strModSummaryFilePath"></param>
        public clsPHRPModSummaryReader(string strModSummaryFilePath)
        {
            mModificationDefs = new List<clsModificationDefinition>();
            mModDefMassesAsText = new Dictionary<string, string>();

            Success = false;

            if (string.IsNullOrEmpty(strModSummaryFilePath))
            {
                throw new Exception("ModSummaryFilePath is empty; unable to continue");
            }

            if (!File.Exists(strModSummaryFilePath))
            {
                throw new FileNotFoundException("ModSummary file not found: " + strModSummaryFilePath);
            }

            Success = ReadModSummaryFile(strModSummaryFilePath, ref mModificationDefs);
        }

        /// <summary>
        /// Returns the mass value associated with the given mass correction tag
        /// </summary>
        /// <param name="strMassCorrectionTag"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public string GetModificationMassAsText(string strMassCorrectionTag)
        {
            if (mModDefMassesAsText.TryGetValue(strMassCorrectionTag, out var modMass))
            {
                return modMass;
            }

            return string.Empty;
        }

        private bool ReadModSummaryFile(string strModSummaryFilePath, ref List<clsModificationDefinition> lstModInfo)
        {
            //string strLineIn = null;
            //string[] splitLine = null;

            //var objColumnHeaders = default(SortedDictionary<string, int>);

            //string modSymbol = null;
            //string modMass = null;
            //string strTargetResidues = null;
            //string modType = null;
            //string strMassCorrectionTag = null;

            //var chModSymbol = default(char);
            //double dblModificationMass = 0;
            //var eModificationType = default(clsModificationDefinition.eModificationTypeConstants);

            //var blnSkipLine = false;
            //var blnHeaderLineParsed = false;

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
            var objColumnHeaders = new SortedDictionary<string, int>(StringComparer.CurrentCultureIgnoreCase)
            {
                {MOD_SUMMARY_COLUMN_Modification_Symbol, 0},
                {MOD_SUMMARY_COLUMN_Modification_Mass, 1},
                {MOD_SUMMARY_COLUMN_Target_Residues, 2},
                {MOD_SUMMARY_COLUMN_Modification_Type, 3},
                {MOD_SUMMARY_COLUMN_Mass_Correction_Tag, 4},
                {MOD_SUMMARY_COLUMN_Occurrence_Count, 5}
            };

            // Read the data from the ModSummary.txt file
            // The first line is typically a header line:
            // Modification_Symbol  Modification_Mass  Target_Residues  Modification_Type  Mass_Correction_Tag  Occurrence_Count

            using (var srModSummaryFile = new StreamReader(new FileStream(strModSummaryFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                var headerLineParsed = false;

                while (!srModSummaryFile.EndOfStream)
                {
                    var lineIn = srModSummaryFile.ReadLine();
                    var skipLine = false;

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    var splitLine = lineIn.Split('\t');

                    if (!headerLineParsed)
                    {
                        if (splitLine[0].ToLower() == MOD_SUMMARY_COLUMN_Modification_Symbol.ToLower())
                        {
                            // Parse the header line to confirm the column ordering
                            // The Occurrence_Count column was misspelled prior to December 2012; need to check for this
                            for (var intIndex = 0; intIndex <= splitLine.Length - 1; intIndex++)
                            {
                                if (splitLine[intIndex] == "Occurence_Count")
                                    splitLine[intIndex] = MOD_SUMMARY_COLUMN_Occurrence_Count;
                            }
                            clsPHRPReader.ParseColumnHeaders(splitLine, objColumnHeaders);
                            skipLine = true;
                        }

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 4)
                        continue;

                    var modSymbol = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Modification_Symbol, objColumnHeaders);
                    var modMassText = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Modification_Mass, objColumnHeaders);
                    var targetResidues = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Target_Residues, objColumnHeaders);
                    var modType = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Modification_Type, objColumnHeaders);
                    var massCorrectionTag = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Mass_Correction_Tag, objColumnHeaders);

                    if (string.IsNullOrWhiteSpace(modSymbol))
                    {
                        modSymbol = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL.ToString();
                    }
                    var chModSymbol = modSymbol[0];

                    if (!double.TryParse(modMassText, out var modificationMass))
                    {
                        throw new Exception("Modification mass is not numeric for MassCorrectionTag: " + massCorrectionTag + ": " + modMassText);
                    }

                    clsModificationDefinition.eModificationTypeConstants eModificationType;
                    if (string.IsNullOrWhiteSpace(modType))
                    {
                        eModificationType = clsModificationDefinition.eModificationTypeConstants.UnknownType;
                    }
                    else
                    {
                        eModificationType = clsModificationDefinition.ModificationSymbolToModificationType(modType[0]);
                    }

                    var objModDef =
                        new clsModificationDefinition(chModSymbol, modificationMass, targetResidues, eModificationType, massCorrectionTag)
                        {
                            ModificationMassAsText = modMassText
                        };

                    lstModInfo.Add(objModDef);

                    if (!mModDefMassesAsText.ContainsKey(massCorrectionTag))
                    {
                        mModDefMassesAsText.Add(massCorrectionTag, modMassText);
                    }
                }
            }

            return true;
        }
    }
}
