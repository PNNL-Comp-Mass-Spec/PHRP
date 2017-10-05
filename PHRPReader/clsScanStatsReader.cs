//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/20/2012
//
// This class reads MASIC ScanStats data from a tab-delimited _ScanStats.txt file
//
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader
{
    public class clsScanStatsReader
    {
        #region "Constants"
        public const string DATA_COLUMN_Dataset = "Dataset";
        public const string DATA_COLUMN_ScanNumber = "ScanNumber";
        public const string DATA_COLUMN_ScanTime = "ScanTime";
        public const string DATA_COLUMN_ScanType = "ScanType";
        public const string DATA_COLUMN_TotalIonIntensity = "TotalIonIntensity";
        public const string DATA_COLUMN_BasePeakIntensity = "BasePeakIntensity";
        public const string DATA_COLUMN_BasePeakMZ = "BasePeakMZ";
        public const string DATA_COLUMN_BasePeakSignalToNoiseRatio = "BasePeakSignalToNoiseRatio";
        public const string DATA_COLUMN_IonCount = "IonCount";
        public const string DATA_COLUMN_IonCountRaw = "IonCountRaw";
        public const string DATA_COLUMN_ScanTypeName = "ScanTypeName";
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
        public clsScanStatsReader()
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
            AddHeaderColumn(DATA_COLUMN_Dataset);
            AddHeaderColumn(DATA_COLUMN_ScanNumber);
            AddHeaderColumn(DATA_COLUMN_ScanTime);
            AddHeaderColumn(DATA_COLUMN_ScanType);
            AddHeaderColumn(DATA_COLUMN_TotalIonIntensity);
            AddHeaderColumn(DATA_COLUMN_BasePeakIntensity);
            AddHeaderColumn(DATA_COLUMN_BasePeakMZ);
            AddHeaderColumn(DATA_COLUMN_BasePeakSignalToNoiseRatio);
            AddHeaderColumn(DATA_COLUMN_IonCount);
            AddHeaderColumn(DATA_COLUMN_IonCountRaw);
            AddHeaderColumn(DATA_COLUMN_ScanTypeName);
        }

        /// <summary>
        /// Open a tab-delimited _ScanStats.txt file and read the data
        /// </summary>
        /// <param name="strInputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are clsScanStatsInfo objects</returns>
        public Dictionary<int, clsScanStatsInfo> ReadScanStatsData(string strInputFilePath)
        {
            Dictionary<int, clsScanStatsInfo> lstScanStats = default(Dictionary<int, clsScanStatsInfo>);
            lstScanStats = new Dictionary<int, clsScanStatsInfo>();

            string strLineIn = null;
            string[] strSplitLine = null;
            bool blnHeaderLineParsed = false;
            bool blnSkipLine = false;

            int intLinesRead = 0;
            int intScanNumber = 0;
            float sngScanTimeMinutes = 0;
            int intScanType = 0;

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
                                intScanNumber = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanNumber, mColumnHeaders, -1);
                                sngScanTimeMinutes = Convert.ToSingle(clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanTime, mColumnHeaders, 0.0));
                                intScanType = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanType, mColumnHeaders, 0);

                                if (intScanNumber >= 0 && !lstScanStats.ContainsKey(intScanNumber))
                                {
                                    clsScanStatsInfo objScanStatsInfo = default(clsScanStatsInfo);
                                    objScanStatsInfo = new clsScanStatsInfo(intScanNumber, sngScanTimeMinutes, intScanType);

                                    objScanStatsInfo.TotalIonIntensity = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_TotalIonIntensity, mColumnHeaders, 0.0);
                                    objScanStatsInfo.BasePeakIntensity = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_BasePeakIntensity, mColumnHeaders, 0.0);
                                    objScanStatsInfo.BasePeakMZ = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_BasePeakMZ, mColumnHeaders, 0.0);
                                    objScanStatsInfo.BasePeakSignalToNoiseRatio = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_BasePeakSignalToNoiseRatio, mColumnHeaders, 0.0);
                                    objScanStatsInfo.IonCount = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_IonCount, mColumnHeaders, 0);
                                    objScanStatsInfo.IonCountRaw = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_IonCountRaw, mColumnHeaders, 0);
                                    objScanStatsInfo.ScanTypeName = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanTypeName, mColumnHeaders);

                                    lstScanStats.Add(intScanNumber, objScanStatsInfo);
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading the ScanStats data: " + ex.Message;
            }

            return lstScanStats;
        }
    }
}
