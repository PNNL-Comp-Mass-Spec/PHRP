//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/03/2012
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// This class reads MSGF scores from a tab-delimited _msgf.txt file
    /// </summary>
    public class MSGFResultsReader
    {
        /// <summary>
        /// Column headers
        /// </summary>
        /// <remarks>Keys are column name, values are 0-based column index</remarks>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Mapping from enum to MSGF file column name
        /// </summary>
        private static readonly Dictionary<MSGFFileColumns, string> mMSGFFileColumn = new();

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
        public MSGFResultsReader()
        {
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        }

        /// <summary>
        /// Define header names
        /// </summary>
        private void DefineColumnHeaders()
        {
            // Define the default column mapping
            var columnHeaders = GetColumnHeaderNamesAndIDs();

            SynFileReaderBaseClass.DefineColumnHeaders(mColumnHeaders, columnHeaders.Keys.ToList());
        }

        /// <summary>
        /// Header names and enums for the _MSGF.txt file
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, MSGFFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, MSGFFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "Result_ID", MSGFFileColumns.ResultID },
                { "Scan", MSGFFileColumns.Scan },
                { "Charge", MSGFFileColumns.Charge },
                { "Protein", MSGFFileColumns.Protein },
                { "Peptide", MSGFFileColumns.Peptide },
                { "SpecProb", MSGFFileColumns.SpecProb },
                { "Notes", MSGFFileColumns.Notes }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum MSGFFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<MSGFFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return SynFileReaderBaseClass.GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(MSGFFileColumns column)
        {
            if (mMSGFFileColumn.Count > 0)
            {
                return mMSGFFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mMSGFFileColumn.Add(item.Value, item.Key);
            }

            return mMSGFFileColumn[column];
        }

        /// <summary>
        /// Open a tab-delimited MSGF results file and read the data
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ResultID and values are MSGF_SpecProb values (stored as strings)</returns>
        public Dictionary<int, string> ReadMSGFData(string inputFilePath)
        {
            var msgfData = new Dictionary<int, string>();

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

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 4)
                    {
                        continue;
                    }

                    var resultID = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(MSGFFileColumns.ResultID), mColumnHeaders, -1);

                    if (resultID < 0)
                    {
                        continue;
                    }

                    var msgfSpecProb = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(MSGFFileColumns.SpecProb), mColumnHeaders);

                    if (!string.IsNullOrEmpty(msgfSpecProb) && !msgfData.ContainsKey(resultID))
                    {
                        msgfData.Add(resultID, msgfSpecProb);
                    }
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading the MSGF data: " + ex.Message;
            }

            return msgfData;
        }
    }
}
