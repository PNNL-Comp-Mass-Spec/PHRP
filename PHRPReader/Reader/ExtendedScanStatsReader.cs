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
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// Extended scan stats reader
    /// </summary>
    public class ExtendedScanStatsReader
    {
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

        /// <summary>
        /// Column headers
        /// </summary>
        private readonly SortedDictionary<string, int> mColumnHeaders;

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
        public ExtendedScanStatsReader()
        {
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);
        }

        private void AddHeaderColumn(string columnName)
        {
            mColumnHeaders.Add(columnName, mColumnHeaders.Count);
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
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are ScanStatsInfo objects</returns>
        public Dictionary<int, ScanStatsExInfo> ReadExtendedScanStatsData(string inputFilePath)
        {
            var scanStats = new Dictionary<int, ScanStatsExInfo>();

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
                        continue;

                    var scanNumber = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ScanNumber, mColumnHeaders, -1);

                    if (scanNumber < 0 || scanStats.ContainsKey(scanNumber))
                        continue;

                    var scanStatsInfo = new ScanStatsExInfo(scanNumber)
                    {
                        IonInjectionTime = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_IonInjectionTime, mColumnHeaders, 0.0),
                        ScanEvent = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ScanEvent, mColumnHeaders, 0),
                        MasterIndex = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_MasterIndex, mColumnHeaders, 0),
                        ElapsedScanTime = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ElapsedScanTime, mColumnHeaders, 0.0),
                        ChargeState = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ChargeState, mColumnHeaders, 0),
                        MonoisotopicMZ = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_MonoisotopicMZ, mColumnHeaders, 0.0),
                        MS2IsolationWidth = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_MS2IsolationWidth, mColumnHeaders, 0.0),
                        FTAnalyzerSettings = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_FTAnalyzerSettings, mColumnHeaders),
                        FTAnalyzerMessage = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_FTAnalyzerMessage, mColumnHeaders),
                        FTResolution = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_FTResolution, mColumnHeaders, 0.0),
                        ConversionParameterB = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ConversionParameterB, mColumnHeaders, 0.0),
                        ConversionParameterC = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ConversionParameterC, mColumnHeaders, 0.0),
                        ConversionParameterD = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ConversionParameterD, mColumnHeaders, 0.0),
                        ConversionParameterE = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ConversionParameterE, mColumnHeaders, 0.0),
                        CollisionMode = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_CollisionMode, mColumnHeaders),
                        ScanFilterText = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_ScanFilterText, mColumnHeaders),
                        SourceVoltage = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_SourceVoltage, mColumnHeaders, 0.0),
                        Source_Current = ReaderFactory.LookupColumnValue(splitLine, DATA_COLUMN_Source_Current, mColumnHeaders, 0.0)
                    };

                    scanStats.Add(scanNumber, scanStatsInfo);
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading the ScanStatsEx data: " + ex.Message;
            }

            return scanStats;
        }
    }
}
