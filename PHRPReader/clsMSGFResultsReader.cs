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

namespace PHRPReader
{
    /// <summary>
    /// This class reads MSGF scores from a tab-delimited _msgf.txt file
    /// </summary>
    public class clsMSGFResultsReader
    {
        #region "Constants"

#pragma warning disable 1591
        public const string DATA_COLUMN_ResultID = "Result_ID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_SpecProb = "SpecProb";
        public const string DATA_COLUMN_Notes = "Notes";
#pragma warning restore 1591

        #endregion

        #region "Class-wide variables"

        /// <summary>
        /// Column headers
        /// </summary>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        private string mErrorMessage = string.Empty;

        #endregion

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
        public clsMSGFResultsReader()
        {
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        }

        private void AddHeaderColumn(string columnName)
        {
            mColumnHeaders.Add(columnName, mColumnHeaders.Count);
        }

        /// <summary>
        /// Define header names for MSGF result files
        /// </summary>
        private void DefineColumnHeaders()
        {
            mColumnHeaders.Clear();

            // Define the default column mapping
            AddHeaderColumn(DATA_COLUMN_ResultID);
            AddHeaderColumn(DATA_COLUMN_Scan);
            AddHeaderColumn(DATA_COLUMN_Charge);
            AddHeaderColumn(DATA_COLUMN_Protein);
            AddHeaderColumn(DATA_COLUMN_Peptide);
            AddHeaderColumn(DATA_COLUMN_SpecProb);
            AddHeaderColumn(DATA_COLUMN_Notes);
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

                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
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
                            if (!clsPHRPReader.IsNumber(splitLine[0]))
                            {
                                // Parse the header line to confirm the column ordering
                                clsPHRPReader.ParseColumnHeaders(splitLine, mColumnHeaders);
                                skipLine = true;
                            }

                            headerLineParsed = true;
                        }

                        if (!skipLine && splitLine.Length >= 4)
                        {
                            var resultID = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_ResultID, mColumnHeaders, -1);

                            if (resultID >= 0)
                            {
                                var msgfSpecProb = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_SpecProb, mColumnHeaders);

                                if (!string.IsNullOrEmpty(msgfSpecProb) && !msgfData.ContainsKey(resultID))
                                {
                                    msgfData.Add(resultID, msgfSpecProb);
                                }
                            }
                        }
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
