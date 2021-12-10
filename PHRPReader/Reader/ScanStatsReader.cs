//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 04/20/2012
//
//*********************************************************************************************************

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PHRPReader.Data;
using PRISM;

namespace PHRPReader.Reader
{
    /// <summary>
    /// This class reads MASIC ScanStats data from a tab-delimited _ScanStats.txt file
    /// </summary>
    public class ScanStatsReader : EventNotifier
    {
        /// <summary>
        /// Column headers
        /// </summary>
        /// <remarks>Keys are column name, values are 0-based column index</remarks>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Mapping from enum to Scan Stats file column name
        /// </summary>
        private static readonly Dictionary<ScanStatsFileColumns, string> mScanStatsFileColumn = new();

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
        public ScanStatsReader()
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
        public static SortedDictionary<string, ScanStatsFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, ScanStatsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                {"Dataset", ScanStatsFileColumns.DatasetId},
                {"ScanNumber", ScanStatsFileColumns.ScanNumber},
                {"ScanTime", ScanStatsFileColumns.ScanTime},
                {"ScanType", ScanStatsFileColumns.ScanType},
                {"TotalIonIntensity", ScanStatsFileColumns.TotalIonIntensity},
                {"BasePeakIntensity", ScanStatsFileColumns.BasePeakIntensity},
                {"BasePeakMZ", ScanStatsFileColumns.BasePeakMZ},
                {"BasePeakSignalToNoiseRatio", ScanStatsFileColumns.BasePeakSignalToNoiseRatio},
                {"IonCount", ScanStatsFileColumns.IonCount},
                {"IonCountRaw", ScanStatsFileColumns.IonCountRaw},
                {"ScanTypeName", ScanStatsFileColumns.ScanTypeName}
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum ScanStatsFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<ScanStatsFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return SynFileReaderBaseClass.GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(ScanStatsFileColumns column)
        {
            if (mScanStatsFileColumn.Count > 0)
            {
                return mScanStatsFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mScanStatsFileColumn.Add(item.Value, item.Key);
            }

            return mScanStatsFileColumn[column];
        }

        /// <summary>
        /// Open a tab-delimited _ScanStats.txt file and read the data
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are instances of ScanStatsInfo</returns>
        public Dictionary<int, ScanStatsInfo> ReadScanStatsData(string inputFilePath)
        {
            var scanStats = new Dictionary<int, ScanStatsInfo>();

            var headerLineParsed = false;

            try
            {
                DefineColumnHeaders();

                mErrorMessage = string.Empty;

                using var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite));

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
                            ReaderFactory.ParseColumnHeaders(splitLine, mColumnHeaders);
                            skipLine = true;
                        }

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 4)
                        continue;

                    var scanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.ScanNumber), mColumnHeaders, -1);
                    var scanTimeMinutes = (float)ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.ScanTime), mColumnHeaders, 0.0);
                    var scanType = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.ScanType), mColumnHeaders, 0);

                    if (scanNumber < 0 || scanStats.ContainsKey(scanNumber))
                        continue;

                    var scanStatsInfo = new ScanStatsInfo(scanNumber, scanTimeMinutes, scanType)
                    {
                        ScanTimeMinutesText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.ScanTime), mColumnHeaders),
                        TotalIonIntensityText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.TotalIonIntensity), mColumnHeaders),
                        BasePeakIntensityText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.BasePeakIntensity), mColumnHeaders),
                        BasePeakMzText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.BasePeakMZ), mColumnHeaders),
                        BasePeakSignalToNoiseRatioText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.BasePeakSignalToNoiseRatio), mColumnHeaders),

                        TotalIonIntensity = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.TotalIonIntensity), mColumnHeaders, 0.0),
                        BasePeakIntensity = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.BasePeakIntensity), mColumnHeaders, 0.0),
                        BasePeakMZ = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.BasePeakMZ), mColumnHeaders, 0.0),
                        BasePeakSignalToNoiseRatio = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.BasePeakSignalToNoiseRatio), mColumnHeaders, 0.0),
                        IonCount = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.IonCount), mColumnHeaders, 0),
                        IonCountRaw = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.IonCountRaw), mColumnHeaders, 0),
                        ScanTypeName = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ScanStatsFileColumns.ScanTypeName), mColumnHeaders)
                    };

                    scanStats.Add(scanNumber, scanStatsInfo);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error reading the ScanStats data", ex);
            }

            return scanStats;
        }

        private void ReportError(string message, Exception ex)
        {
            OnErrorEvent(message, ex);
            mErrorMessage = string.Format("{0}{1}", message, ex == null ? string.Empty : ": " + ex.Message);
        }
    }
}
