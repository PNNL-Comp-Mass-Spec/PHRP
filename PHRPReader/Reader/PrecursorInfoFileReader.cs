using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PHRPReader.Data;

// ReSharper disable UnusedMember.Global

namespace PHRPReader.Reader
{
    /// <summary>
    /// This class reads a tab-delimited _PrecursorInfo.txt file (created by the Analysis Manager)
    /// </summary>
    public class PrecursorInfoFileReader
    {
        private const string ALTERNATE_SCAN_FILTER_COLUMN_NAME = "Scan Filter Text";

        /// <summary>
        /// Column headers
        /// </summary>
        /// <remarks>Keys are column name, values are 0-based column index</remarks>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Mapping from enum to Precursor Info file column name
        /// </summary>
        private static readonly Dictionary<PrecursorInfoFileColumns, string> mPrecursorInfoFileColumn = new();

        private string mErrorMessage = string.Empty;

        /// <summary>
        /// Error message
        /// </summary>
        public string ErrorMessage
        {
            get
            {
                if (string.IsNullOrEmpty(mErrorMessage))
                {
                    return string.Empty;
                }

                return mErrorMessage;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        public PrecursorInfoFileReader()
        {
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        }

        /// <summary>
        /// Define header names
        /// </summary>
        private void DefineColumnHeaders()
        {
            // Define the default column mapping
            var columnHeaders = GetColumnHeaderNamesAndIDs(true);

            SynFileReaderBaseClass.DefineColumnHeaders(mColumnHeaders, columnHeaders.Keys.ToList());
        }

        /// <summary>
        /// Header names and enums for the _PrecursorInfo.txt file
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, PrecursorInfoFileColumns> GetColumnHeaderNamesAndIDs(bool includeSynonyms = false)
        {
            var headerColumns = new SortedDictionary<string, PrecursorInfoFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "ScanNumber", PrecursorInfoFileColumns.ScanNumber },
                { "ScanTime", PrecursorInfoFileColumns.ScanTime },
                { "ScanType", PrecursorInfoFileColumns.ScanType },
                { "ScanTypeName", PrecursorInfoFileColumns.ScanTypeName },
                { "PrecursorMz", PrecursorInfoFileColumns.PrecursorMz },
                { "ScanFilterText", PrecursorInfoFileColumns.ScanFilterText }
            };

            if (includeSynonyms)
            {
                headerColumns.Add(ALTERNATE_SCAN_FILTER_COLUMN_NAME, PrecursorInfoFileColumns.ScanFilterText);
            }

            return headerColumns;
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum PrecursorInfoFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        public static Dictionary<PrecursorInfoFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return SynFileReaderBaseClass.GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(PrecursorInfoFileColumns column)
        {
            if (mPrecursorInfoFileColumn.Count > 0)
            {
                return mPrecursorInfoFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mPrecursorInfoFileColumn.Add(item.Value, item.Key);
            }

            return mPrecursorInfoFileColumn[column];
        }

        /// <summary>
        /// Open a tab-delimited _PrecursorInfo.txt file and read the data
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are instances of the PrecursorInfo class</returns>
        public Dictionary<int, PrecursorInfo> ReadPrecursorInfoFile(string inputFilePath)
        {
            var precursorInfoData = new Dictionary<int, PrecursorInfo>();

            try
            {
                DefineColumnHeaders();

                mErrorMessage = string.Empty;

                using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

                var headerLineParsed = false;

                while (!reader.EndOfStream)
                {
                    var lineIn = reader.ReadLine();
                    var skipLine = false;

                    if (string.IsNullOrWhiteSpace(lineIn))
                        continue;

                    var splitLine = lineIn.Split('\t');

                    if (!headerLineParsed)
                    {
                        if (!ReaderFactory.IsNumber(splitLine[0]))
                        {
                            // Parse the header line to confirm the column ordering
                            ReaderFactory.ParseColumnHeaders(splitLine, mColumnHeaders);
                            skipLine = true;
                        }

                        // Call this method in order to populate mPrecursorInfoFileColumn
                        // This is required prior to checking for the presence of the alternate "Scan Filter Text" column name
                        GetColumnNameByID(PrecursorInfoFileColumns.ScanNumber);

                        foreach (var item in mColumnHeaders)
                        {
                            if (item.Key.Equals(ALTERNATE_SCAN_FILTER_COLUMN_NAME) && item.Value >= 0)
                            {
                                // The file has "Scan Filter Text" instead of "ScanFilterText"
                                mPrecursorInfoFileColumn[PrecursorInfoFileColumns.ScanFilterText] = ALTERNATE_SCAN_FILTER_COLUMN_NAME;
                                break;
                            }
                        }

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 4)
                    {
                        continue;
                    }

                    var scanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(PrecursorInfoFileColumns.ScanNumber), mColumnHeaders, -1);
                    var scanTimeMinutes = (float)ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(PrecursorInfoFileColumns.ScanTime), mColumnHeaders, 0.0);
                    var scanType = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(PrecursorInfoFileColumns.ScanType), mColumnHeaders, 0);

                    if (scanNumber < 0 || precursorInfoData.ContainsKey(scanNumber))
                        continue;

                    var precursorInfo = new PrecursorInfo(scanNumber, scanTimeMinutes, scanType)
                    {
                        ScanTypeName = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(PrecursorInfoFileColumns.ScanTypeName), mColumnHeaders, string.Empty),
                        PrecursorMz = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(PrecursorInfoFileColumns.PrecursorMz), mColumnHeaders, 0.0),
                        ScanFilterText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(PrecursorInfoFileColumns.ScanFilterText), mColumnHeaders, string.Empty)
                    };

                    precursorInfoData.Add(scanNumber, precursorInfo);

                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading the _PrecursorInfo.txt file: " + ex.Message;
            }

            return precursorInfoData;
        }
    }
}
