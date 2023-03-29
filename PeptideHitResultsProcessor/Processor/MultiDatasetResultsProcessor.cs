using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using PRISM;

namespace PeptideHitResultsProcessor.Processor
{
    /// <summary>
    /// Base class for MSFraggerResultsProcessor and MaxQuantResultsProcessor
    /// </summary>
    public abstract class MultiDatasetResultsProcessor : PHRPBaseClass
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="options"></param>
        protected MultiDatasetResultsProcessor(PHRPOptions options) : base(options)
        {
        }

        private static string CombineDatasetNameParts(
            string datasetName,
            IReadOnlyList<string> datasetNameParts,
            int partCountToUse,
            int minimumLength = 0,
            int maximumLengthThreshold = 0)
        {
            if (partCountToUse >= datasetNameParts.Count)
                return datasetName;

            var nextIndex = 0;
            var combinedName = new StringBuilder();

            while (nextIndex < datasetNameParts.Count)
            {
                if (nextIndex >= partCountToUse)
                {
                    if (maximumLengthThreshold == 0 && minimumLength == 0)
                    {
                        // Minimum or maximum length not defined, and we've used the required number of parts
                        break;
                    }

                    var tooShort = minimumLength > 0 && combinedName.ToString().Length < minimumLength;

                    if (maximumLengthThreshold > 0)
                    {
                        // Maximum length defined
                        if (!tooShort && combinedName.ToString().Length + datasetNameParts[nextIndex].Length > maximumLengthThreshold)
                        {
                            // Adding the next part will result in the total length exceeding the maximum
                            // Do not add any more name parts
                            break;
                        }
                    }

                    if (!tooShort)
                    {
                        // Minimum length defined and the name is now long enough
                        break;
                    }
                }

                combinedName.Append(datasetNameParts[nextIndex]);
                nextIndex++;
            }

            return combinedName.ToString();
        }

        /// <summary>
        /// Get the base name to use for output files
        /// </summary>
        /// <remarks>If baseNameByDatasetName is empty, or if longestCommonBaseName and Options.OutputFileBaseName are empty, uses a generic base name</remarks>
        /// <param name="baseNameByDatasetName"></param>
        /// <param name="toolNameAbbreviation"></param>
        /// <param name="longestCommonBaseName"></param>
        /// <returns>Base name</returns>
        protected string GetBaseNameForOutputFiles(Dictionary<string, string> baseNameByDatasetName, string toolNameAbbreviation, string longestCommonBaseName)
        {
            // ReSharper disable once ConvertIfStatementToSwitchStatement
            if (baseNameByDatasetName.Count == 0 || baseNameByDatasetName.Count > 1 && string.IsNullOrWhiteSpace(Options.OutputFileBaseName) && string.IsNullOrWhiteSpace(longestCommonBaseName))
            {
                return string.Format("Dataset_{0}", toolNameAbbreviation);
            }

            if (baseNameByDatasetName.Count == 1)
            {
                return string.Format("{0}_{1}", baseNameByDatasetName.First().Key, toolNameAbbreviation);
            }

            // baseNameByDatasetName.Count is greater than 1 and either longestCommonBaseName has text or Options.OutputFileBaseName has text

            return string.Format("{0}_{1}",
                string.IsNullOrWhiteSpace(Options.OutputFileBaseName) ? longestCommonBaseName : Options.OutputFileBaseName,
                toolNameAbbreviation);
        }

        /// <summary>
        /// Examine the names in datasetNames
        /// Create a mapping from full name to abbreviated name
        /// </summary>
        /// <param name="inputFileName">Input file name, which is used if datasetNames is empty</param>
        /// <param name="datasetNames"></param>
        /// <param name="longestCommonBaseName"></param>
        /// <param name="filenameSuffixToRemove"></param>
        /// <returns>Dictionary where keys are dataset names and values are abbreviated names</returns>
        protected Dictionary<string, string> GetDatasetNameMap(
            string inputFileName,
            SortedSet<string> datasetNames,
            out string longestCommonBaseName,
            string filenameSuffixToRemove = "")
        {
            if (datasetNames.Count == 0)
            {
                var baseName = Path.GetFileNameWithoutExtension(inputFileName);

                // Remove _psm if present
                if (!string.IsNullOrWhiteSpace(filenameSuffixToRemove) &&
                    baseName.EndsWith(Path.GetFileNameWithoutExtension(filenameSuffixToRemove), StringComparison.OrdinalIgnoreCase) &&
                    baseName.Length > 4)
                {
                    datasetNames.Add(baseName.Substring(0, baseName.Length - 4));
                }
                else
                {
                    datasetNames.Add(baseName);
                }
            }

            var baseNameByDatasetName = GetDatasetNameMap(datasetNames, out longestCommonBaseName, out var warnings);

            foreach (var warning in warnings)
            {
                OnWarningEvent(warning);
            }

            return baseNameByDatasetName;
        }

        /// <summary>
        /// Examine the names in datasetNames
        /// Create a mapping from full name to abbreviated name
        /// </summary>
        /// <param name="datasetNames"></param>
        /// <param name="longestCommonBaseName">Output: longest common base name</param>
        /// <param name="warnings">Output: warning messages</param>
        /// <returns>Dictionary where keys are dataset names and values are abbreviated names</returns>
        public static Dictionary<string, string> GetDatasetNameMap(
            SortedSet<string> datasetNames,
            out string longestCommonBaseName,
            out List<string> warnings)
        {
            warnings = new List<string>();

            var datasetNameParts = new Dictionary<string, List<string>>();
            var maxPartCount = 0;
            var splitChars = new[] { '_', '-' };

            foreach (var datasetName in datasetNames)
            {
                var nameParts = new List<string>();
                var startIndex = 0;
                while (startIndex < datasetName.Length)
                {
                    var matchIndex = datasetName.IndexOfAny(splitChars, startIndex + 1);
                    if (matchIndex < 0)
                    {
                        nameParts.Add(datasetName.Substring(startIndex));
                        break;
                    }

                    nameParts.Add(datasetName.Substring(startIndex, matchIndex - startIndex));
                    startIndex = matchIndex;
                }

                datasetNameParts.Add(datasetName, nameParts);

                maxPartCount = Math.Max(maxPartCount, nameParts.Count);
            }

            if (datasetNameParts.Count == 0)
            {
                longestCommonBaseName = string.Empty;
                return new Dictionary<string, string>();
            }

            var candidateBaseNames = new SortedSet<string>();

            var partCountToUse = 1;

            var datasetNameKeys = datasetNameParts.Keys.ToList();

            while (partCountToUse <= maxPartCount)
            {
                candidateBaseNames.Clear();
                candidateBaseNames.Add(CombineDatasetNameParts(datasetNameKeys[0], datasetNameParts[datasetNameKeys[0]], partCountToUse));

                for (var i = 1; i < datasetNameKeys.Count; i++)
                {
                    var baseNameToAdd = CombineDatasetNameParts(datasetNameKeys[i], datasetNameParts[datasetNameKeys[i]], partCountToUse);
                    if (candidateBaseNames.Contains(baseNameToAdd))
                    {
                        // Name collision found
                        break;
                    }

                    candidateBaseNames.Add(baseNameToAdd);
                }

                if (candidateBaseNames.Count == datasetNameKeys.Count)
                    break;

                partCountToUse++;
            }

            var baseDatasetNames = new SortedSet<string>();

            // Dictionary where keys are dataset names and values are abbreviated names
            var baseNameByDatasetName = new Dictionary<string, string>();

            if (candidateBaseNames.Count == datasetNameKeys.Count)
            {
                // Can use a subsection of the dataset name(s)
                // Combine subsections to create the base name for each dataset
                foreach (var item in datasetNameParts)
                {
                    var baseNameToAdd = CombineDatasetNameParts(item.Key, item.Value, partCountToUse, 12, 25);
                    baseNameByDatasetName.Add(item.Key, baseNameToAdd);

                    if (baseDatasetNames.Contains(baseNameToAdd))
                    {
                        warnings.Add(string.Format(
                            "Warning: baseDatasetNames already contains: {0}\nLogic error for dataset {1}",
                            baseNameToAdd, item.Key));

                        continue;
                    }

                    baseDatasetNames.Add(baseNameToAdd);
                }
            }
            else
            {
                // Not able to shorten the dataset names since they are too similar
                // Use full dataset names
                foreach (var item in datasetNameParts)
                {
                    baseNameByDatasetName.Add(item.Key, item.Key);
                    baseDatasetNames.Add(item.Key);
                }
            }

            longestCommonBaseName = StringUtilities.LongestCommonStringFromStart(baseNameByDatasetName.Values.ToList());
            longestCommonBaseName = longestCommonBaseName.TrimEnd('_', '-');

            if (longestCommonBaseName.Length > 7 && (
                longestCommonBaseName.EndsWith("_0") ||
                longestCommonBaseName.EndsWith("_f")))
            {
                longestCommonBaseName = longestCommonBaseName.Substring(0, longestCommonBaseName.Length - 2);
            }

            return baseNameByDatasetName;
        }

        /// <summary>
        /// Lookup the base name and dataset ID for the given dataset name
        /// </summary>
        /// <param name="baseNameByDatasetName"></param>
        /// <param name="datasetIDs"></param>
        /// <param name="datasetName"></param>
        /// <param name="baseDatasetName">Output: base dataset name, or empty string if not found</param>
        /// <param name="datasetID">Output: dataset ID, or 0 if not found</param>
        protected void GetBaseNameAndDatasetID(
            Dictionary<string, string> baseNameByDatasetName,
            Dictionary<string, int> datasetIDs,
            string datasetName,
            out string baseDatasetName,
            out int datasetID)
        {
            if (string.IsNullOrWhiteSpace(datasetName))
            {
                baseDatasetName = string.Empty;
                datasetID = 0;
                return;
            }

            if (!baseNameByDatasetName.TryGetValue(datasetName, out baseDatasetName))
            {
                ConsoleMsgUtils.ShowDebug(
                    "The baseNameByDatasetName dictionary does not contain key {0}; this is unexpected",
                    datasetName);
            }

            if (!datasetIDs.TryGetValue(datasetName, out datasetID))
            {
                datasetID = 0;
            }
        }
    }
}
