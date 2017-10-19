//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/20/2012
//
// This class reads MASIC Extended ScanStats data from a tab-delimited _ScanStatsEx.txt file
//
//*********************************************************************************************************
using System;
using System.Collections.Generic;
using System.IO;

namespace PHRPReader
{
    /// <summary>
    /// Extended scan stats reader
    /// </summary>
    public class clsExtendedScanStatsReader
    {
        #region "Constants"

#pragma warning disable 1591
        public const string DATA_COLUMN_Dataset = "Dataset";
        public const string DATA_COLUMN_ScanNumber = "ScanNumber";
        public const string DATA_COLUMN_IonInjectionTime = "Ion Injection Time (ms)";
        public const string DATA_COLUMN_ScanEvent = "Scan Event";
        public const string DATA_COLUMN_MasterIndex = "Master Index";
        public const string DATA_COLUMN_ElapsedScanTime = "Elapsed Scan Time (sec)";
        public const string DATA_COLUMN_ChargeState = "Charge State";
        public const string DATA_COLUMN_MonoisotopicMZ = "Monoisotopic M/Z";
        public const string DATA_COLUMN_MS2IsolationWidth = "MS2 Isolation Width";
        public const string DATA_COLUMN_FTAnalyzerSettings = "FT Analyzer Settings";
        public const string DATA_COLUMN_FTAnalyzerMessage = "FT Analyzer Message";
        public const string DATA_COLUMN_FTResolution = "FT Resolution";
        public const string DATA_COLUMN_ConversionParameterB = "Conversion Parameter B";
        public const string DATA_COLUMN_ConversionParameterC = "Conversion Parameter C";
        public const string DATA_COLUMN_ConversionParameterD = "Conversion Parameter D";
        public const string DATA_COLUMN_ConversionParameterE = "Conversion Parameter E";
        public const string DATA_COLUMN_CollisionMode = "Collision Mode";
        public const string DATA_COLUMN_ScanFilterText = "Scan Filter Text";
        public const string DATA_COLUMN_SourceVoltage = "Source Voltage (kV)";
        public const string DATA_COLUMN_Source_Current = "Source Current (uA)";
#pragma warning restore 1591

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
        public clsExtendedScanStatsReader()
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
            AddHeaderColumn(DATA_COLUMN_IonInjectionTime);
            AddHeaderColumn(DATA_COLUMN_ScanEvent);
            AddHeaderColumn(DATA_COLUMN_MasterIndex);
            AddHeaderColumn(DATA_COLUMN_ElapsedScanTime);
            AddHeaderColumn(DATA_COLUMN_ChargeState);
            AddHeaderColumn(DATA_COLUMN_MonoisotopicMZ);
            AddHeaderColumn(DATA_COLUMN_MS2IsolationWidth);
            AddHeaderColumn(DATA_COLUMN_FTAnalyzerSettings);
            AddHeaderColumn(DATA_COLUMN_FTAnalyzerMessage);
            AddHeaderColumn(DATA_COLUMN_FTResolution);
            AddHeaderColumn(DATA_COLUMN_ConversionParameterB);
            AddHeaderColumn(DATA_COLUMN_ConversionParameterC);
            AddHeaderColumn(DATA_COLUMN_ConversionParameterD);
            AddHeaderColumn(DATA_COLUMN_ConversionParameterE);
            AddHeaderColumn(DATA_COLUMN_CollisionMode);
            AddHeaderColumn(DATA_COLUMN_ScanFilterText);
            AddHeaderColumn(DATA_COLUMN_SourceVoltage);
            AddHeaderColumn(DATA_COLUMN_Source_Current);
        }

        /// <summary>
        /// Open a tab-delimited _ScanStatsEx.txt file and read the data
        /// </summary>
        /// <param name="strInputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are clsScanStatsInfo objects</returns>
        public Dictionary<int, clsScanStatsExInfo> ReadExtendedScanStatsData(string strInputFilePath)
        {
            var lstScanStats = default(Dictionary<int, clsScanStatsExInfo>);
            lstScanStats = new Dictionary<int, clsScanStatsExInfo>();

            string strLineIn = null;
            string[] strSplitLine = null;
            var blnHeaderLineParsed = false;
            var blnSkipLine = false;

            var intLinesRead = 0;
            var intScanNumber = 0;

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

                                if (intScanNumber >= 0 && !lstScanStats.ContainsKey(intScanNumber))
                                {
                                    var objScanStatsInfo = default(clsScanStatsExInfo);
                                    objScanStatsInfo = new clsScanStatsExInfo(intScanNumber);

                                    objScanStatsInfo.IonInjectionTime = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_IonInjectionTime, mColumnHeaders, 0.0);
                                    objScanStatsInfo.ScanEvent = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanEvent, mColumnHeaders, 0);
                                    objScanStatsInfo.MasterIndex = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_MasterIndex, mColumnHeaders, 0);
                                    objScanStatsInfo.ElapsedScanTime = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ElapsedScanTime, mColumnHeaders, 0.0);
                                    objScanStatsInfo.ChargeState = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ChargeState, mColumnHeaders, 0);
                                    objScanStatsInfo.MonoisotopicMZ = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_MonoisotopicMZ, mColumnHeaders, 0.0);
                                    objScanStatsInfo.MS2IsolationWidth = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_MS2IsolationWidth, mColumnHeaders, 0.0);
                                    objScanStatsInfo.FTAnalyzerSettings = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_FTAnalyzerSettings, mColumnHeaders);
                                    objScanStatsInfo.FTAnalyzerMessage = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_FTAnalyzerMessage, mColumnHeaders);
                                    objScanStatsInfo.FTResolution = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_FTResolution, mColumnHeaders, 0.0);
                                    objScanStatsInfo.ConversionParameterB = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterB, mColumnHeaders, 0.0);
                                    objScanStatsInfo.ConversionParameterC = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterC, mColumnHeaders, 0.0);
                                    objScanStatsInfo.ConversionParameterD = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterD, mColumnHeaders, 0.0);
                                    objScanStatsInfo.ConversionParameterE = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ConversionParameterE, mColumnHeaders, 0.0);
                                    objScanStatsInfo.CollisionMode = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_CollisionMode, mColumnHeaders);
                                    objScanStatsInfo.ScanFilterText = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_ScanFilterText, mColumnHeaders);
                                    objScanStatsInfo.SourceVoltage = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_SourceVoltage, mColumnHeaders, 0.0);
                                    objScanStatsInfo.Source_Current = clsPHRPReader.LookupColumnValue(strSplitLine, DATA_COLUMN_Source_Current, mColumnHeaders, 0.0);

                                    lstScanStats.Add(intScanNumber, objScanStatsInfo);
                                }
                            }
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading the ScanStatsEx data: " + ex.Message;
            }

            return lstScanStats;
        }
    }
}
