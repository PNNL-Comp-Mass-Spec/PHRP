//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/03/2012
//
// This class reads MSGF scores from a tab-delimited _msgf.txt file
//
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader
{
    public class clsMSGFResultsReader
    {
        #region "Constants"
        public const string DATA_COLUMN_ResultID = "Result_ID";
        public const string DATA_COLUMN_Scan = "Scan";
        public const string DATA_COLUMN_Charge = "Charge";
        public const string DATA_COLUMN_Protein = "Protein";
        public const string DATA_COLUMN_Peptide = "Peptide";
        public const string DATA_COLUMN_SpecProb = "SpecProb";
        public const string DATA_COLUMN_Notes = "Notes";
        #endregion

        #region "Class-wide variables"
        // Column headers
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
                else
                {
                    return mErrorMessage;
                }
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        public clsMSGFResultsReader()
        {
            mColumnHeaders = new SortedDictionary<string, int>();
        }

        private void AddHeaderColumn(string strColumnName)
        {
            mColumnHeaders.Add(strColumnName, mColumnHeaders.Count);
        }

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
        /// <param name="strInputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ResultID and values are MSGF_SpecProb values (stored as strings)</returns>
        public Dictionary<int, string> ReadMSGFData(string strInputFilePath)
        {
            var lstMSGFData = default(Dictionary<int, string>);
            lstMSGFData = new Dictionary<int, string>();

            string strLineIn = null;
            string[] strSplitLine = null;
            var blnHeaderLineParsed = false;
            var blnSkipLine = false;

            var intLinesRead = 0;
            var intResultID = 0;
            string strMSGFSpecProb = null;

            try
            {
                DefineColumnHeaders();
                intLinesRead = 0;
                mErrorMessage = string.Empty;

                using (var srInFile = new StreamReader(new FileStream(strInputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srInFile.EndOfStream)
                    {
                        strLineIn = srInFile.ReadLine();
                        intLinesRead += 1;
                        blnSkipLine = false;

                        if (!string.IsNullOrWhiteSpace(strLineIn))
                        {
                            strSplitLine = strLineIn.Split('\t');

                            if (!blnHeaderLineParsed)
                            {
                                if (!clsPHRPReader.IsNumber(strSplitLine[0]))
                                {
                                    // Parse the header line to confirm the column ordering
                                    clsPHRPReader.ParseColumnHeaders(strSplitLine, mColumnHeaders);
                                    blnSkipLine = true;
                                }

                                blnHeaderLineParsed = true;
                            }

                            if (!blnSkipLine && strSplitLine.Length >= 4)
                            {
                                intResultID = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ResultID, mColumnHeaders, -1);

                                if (intResultID >= 0)
                                {
                                    strMSGFSpecProb = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_SpecProb, mColumnHeaders);

                                    if (!string.IsNullOrEmpty(strMSGFSpecProb) && !lstMSGFData.ContainsKey(intResultID))
                                    {
                                        lstMSGFData.Add(intResultID, strMSGFSpecProb);
                                    }
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

            return lstMSGFData;
        }
    }
}
