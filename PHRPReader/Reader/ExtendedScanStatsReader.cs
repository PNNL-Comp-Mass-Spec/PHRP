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
using System.Linq;
using PHRPReader.Data;

namespace PHRPReader.Reader
{
    /// <summary>
    /// Extended scan stats reader
    /// </summary>
    public class ExtendedScanStatsReader
    {
        /// <summary>
        /// Column headers
        /// </summary>
        /// <remarks>Keys are column name, values are 0-based column index</remarks>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Mapping from enum to Extended Scan Stats file column name
        /// </summary>
        private static readonly Dictionary<ExtendedScanStatsFileColumns, string> mExtendedScanStatsFileColumn = new();

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
        /// Header names and enums for the _ScanStats.txt file
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, ExtendedScanStatsFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new(StringComparer.OrdinalIgnoreCase)
            {
                { "Dataset", ExtendedScanStatsFileColumns.DatasetId },
                { "ScanNumber", ExtendedScanStatsFileColumns.ScanNumber },
                { "Ion Injection Time (ms)", ExtendedScanStatsFileColumns.IonInjectionTime },
                { "Scan Event", ExtendedScanStatsFileColumns.ScanEvent },
                { "Master Index", ExtendedScanStatsFileColumns.MasterIndex },
                { "Elapsed Scan Time (sec)", ExtendedScanStatsFileColumns.ElapsedScanTime },
                { "Charge State", ExtendedScanStatsFileColumns.ChargeState },
                { "Monoisotopic M/Z", ExtendedScanStatsFileColumns.MonoisotopicMZ },
                { "MS2 Isolation Width", ExtendedScanStatsFileColumns.MS2IsolationWidth },
                { "FT Analyzer Settings", ExtendedScanStatsFileColumns.FTAnalyzerSettings },
                { "FT Analyzer Message", ExtendedScanStatsFileColumns.FTAnalyzerMessage },
                { "FT Resolution", ExtendedScanStatsFileColumns.FTResolution },
                { "Conversion Parameter B", ExtendedScanStatsFileColumns.ConversionParameterB },
                { "Conversion Parameter C", ExtendedScanStatsFileColumns.ConversionParameterC },
                { "Conversion Parameter D", ExtendedScanStatsFileColumns.ConversionParameterD },
                { "Conversion Parameter E", ExtendedScanStatsFileColumns.ConversionParameterE },
                { "Collision Mode", ExtendedScanStatsFileColumns.CollisionMode },
                { "Scan Filter Text", ExtendedScanStatsFileColumns.ScanFilterText },
                { "Source Voltage (kV)", ExtendedScanStatsFileColumns.SourceVoltage },
                { "Source Current (uA)", ExtendedScanStatsFileColumns.Source_Current }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum ExtendedScanStatsFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<ExtendedScanStatsFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return SynFileReaderBaseClass.GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(ExtendedScanStatsFileColumns column)
        {
            if (mExtendedScanStatsFileColumn.Count > 0)
            {
                return mExtendedScanStatsFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mExtendedScanStatsFileColumn.Add(item.Value, item.Key);
            }

            return mExtendedScanStatsFileColumn[column];
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

                    var scanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ScanNumber), mColumnHeaders, -1);

                    if (scanNumber < 0 || scanStats.ContainsKey(scanNumber))
                        continue;

                    var scanStatsInfo = new ScanStatsExInfo(scanNumber)
                    {
                        IonInjectionTime = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.IonInjectionTime), mColumnHeaders, 0.0),
                        ScanEvent = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ScanEvent), mColumnHeaders, 0),
                        MasterIndex = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.MasterIndex), mColumnHeaders, 0),
                        ElapsedScanTime = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ElapsedScanTime), mColumnHeaders, 0.0),
                        ChargeState = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ChargeState), mColumnHeaders, 0),
                        MonoisotopicMZ = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.MonoisotopicMZ), mColumnHeaders, 0.0),
                        MS2IsolationWidth = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.MS2IsolationWidth), mColumnHeaders, 0.0),
                        FTAnalyzerSettings = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.FTAnalyzerSettings), mColumnHeaders),
                        FTAnalyzerMessage = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.FTAnalyzerMessage), mColumnHeaders),
                        FTResolution = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.FTResolution), mColumnHeaders, 0.0),
                        ConversionParameterB = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ConversionParameterB), mColumnHeaders, 0.0),
                        ConversionParameterC = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ConversionParameterC), mColumnHeaders, 0.0),
                        ConversionParameterD = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ConversionParameterD), mColumnHeaders, 0.0),
                        ConversionParameterE = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ConversionParameterE), mColumnHeaders, 0.0),
                        CollisionMode = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.CollisionMode), mColumnHeaders),
                        ScanFilterText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.ScanFilterText), mColumnHeaders),
                        SourceVoltage = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.SourceVoltage), mColumnHeaders, 0.0),
                        Source_Current = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ExtendedScanStatsFileColumns.Source_Current), mColumnHeaders, 0.0)
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
