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
        /// <param name="modSummaryFilePath"></param>
        public clsPHRPModSummaryReader(string modSummaryFilePath)
        {
            mModificationDefs = new List<clsModificationDefinition>();
            mModDefMassesAsText = new Dictionary<string, string>();

            Success = false;

            if (string.IsNullOrEmpty(modSummaryFilePath))
            {
                throw new Exception("ModSummaryFilePath is empty; unable to continue");
            }

            if (!File.Exists(modSummaryFilePath))
            {
                throw new FileNotFoundException("ModSummary file not found: " + modSummaryFilePath);
            }

            Success = ReadModSummaryFile(modSummaryFilePath);
        }

        /// <summary>
        /// Returns the mass value associated with the given mass correction tag
        /// </summary>
        /// <param name="massCorrectionTag"></param>
        public string GetModificationMassAsText(string massCorrectionTag)
        {
            if (mModDefMassesAsText.TryGetValue(massCorrectionTag, out var modMass))
            {
                return modMass;
            }

            return string.Empty;
        }

        private bool ReadModSummaryFile(string modSummaryFilePath)
        {

            mModificationDefs.Clear();

            if (string.IsNullOrEmpty(modSummaryFilePath))
            {
                return false;
            }

            // Initialize the column mapping
            // Using a case-insensitive comparer
            var columnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase)
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

            using (var modSummaryReader = new StreamReader(new FileStream(modSummaryFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
            {
                var headerLineParsed = false;

                while (!modSummaryReader.EndOfStream)
                {
                    var lineIn = modSummaryReader.ReadLine();
                    var skipLine = false;

                    if (string.IsNullOrEmpty(lineIn))
                        continue;

                    var splitLine = lineIn.Split('\t');

                    if (!headerLineParsed)
                    {
                        if (string.Equals(splitLine[0], MOD_SUMMARY_COLUMN_Modification_Symbol, StringComparison.OrdinalIgnoreCase))
                        {
                            // Parse the header line to confirm the column ordering
                            // The Occurrence_Count column was misspelled prior to December 2012; need to check for this
                            for (var index = 0; index <= splitLine.Length - 1; index++)
                            {
                                if (splitLine[index] == "Occurence_Count")
                                    splitLine[index] = MOD_SUMMARY_COLUMN_Occurrence_Count;
                            }
                            clsPHRPReader.ParseColumnHeaders(splitLine, columnHeaders);
                            skipLine = true;
                        }

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 4)
                        continue;

                    var modSymbolText = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Modification_Symbol, columnHeaders);
                    var modMassText = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Modification_Mass, columnHeaders);
                    var targetResidues = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Target_Residues, columnHeaders);
                    var modType = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Modification_Type, columnHeaders);
                    var massCorrectionTag = clsPHRPReader.LookupColumnValue(splitLine, MOD_SUMMARY_COLUMN_Mass_Correction_Tag, columnHeaders);

                    if (string.IsNullOrWhiteSpace(modSymbolText))
                    {
                        modSymbolText = clsModificationDefinition.NO_SYMBOL_MODIFICATION_SYMBOL.ToString();
                    }

                    var modSymbol = modSymbolText[0];

                    if (!double.TryParse(modMassText, out var modificationMass))
                    {
                        throw new Exception("Modification mass is not numeric for MassCorrectionTag: " + massCorrectionTag + ": " + modMassText);
                    }

                    clsModificationDefinition.ModificationTypeConstants modificationType;
                    if (string.IsNullOrWhiteSpace(modType))
                    {
                        modificationType = clsModificationDefinition.ModificationTypeConstants.UnknownType;
                    }
                    else
                    {
                        modificationType = clsModificationDefinition.ModificationSymbolToModificationType(modType[0]);
                    }

                    var modDef =
                        new clsModificationDefinition(modSymbol, modificationMass, targetResidues, modificationType, massCorrectionTag)
                        {
                            ModificationMassAsText = modMassText
                        };

                    mModificationDefs.Add(modDef);

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
