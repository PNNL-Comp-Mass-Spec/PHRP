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

namespace PHRPReader
{
    /// <summary>
    /// This class reads MASIC ScanStats data from a tab-delimited _ScanStats.txt file
    /// </summary>
    public class clsScanStatsReader
    {
        #region "Constants"

#pragma warning disable 1591
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

                return mErrorMessage;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        public clsScanStatsReader()
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
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are clsScanStatsInfo objects</returns>
        public Dictionary<int, clsScanStatsInfo> ReadScanStatsData(string inputFilePath)
        {
            var scanStats = new Dictionary<int, clsScanStatsInfo>();

            var headerLineParsed = false;

            try
            {
                DefineColumnHeaders();

                mErrorMessage = string.Empty;

                using (var reader = new StreamReader(new FileStream(inputFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
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

                        if (skipLine || splitLine.Length < 4)
                            continue;

                        var scanNumber = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_ScanNumber, mColumnHeaders, -1);
                        var scanTimeMinutes = Convert.ToSingle(clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_ScanTime, mColumnHeaders, 0.0));
                        var scanType = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_ScanType, mColumnHeaders, 0);

                        if (scanNumber < 0 || scanStats.ContainsKey(scanNumber))
                            continue;

                        var scanStatsInfo = new clsScanStatsInfo(scanNumber, scanTimeMinutes, scanType)
                        {
                            TotalIonIntensity = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_TotalIonIntensity, mColumnHeaders, 0.0),
                            BasePeakIntensity = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_BasePeakIntensity, mColumnHeaders, 0.0),
                            BasePeakMZ = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_BasePeakMZ, mColumnHeaders, 0.0),
                            BasePeakSignalToNoiseRatio = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_BasePeakSignalToNoiseRatio, mColumnHeaders, 0.0),
                            IonCount = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_IonCount, mColumnHeaders, 0),
                            IonCountRaw = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_IonCountRaw, mColumnHeaders, 0),
                            ScanTypeName = clsPHRPReader.LookupColumnValue(splitLine, DATA_COLUMN_ScanTypeName, mColumnHeaders)
                        };

                        scanStats.Add(scanNumber, scanStatsInfo);
                    }
                }
            }
            catch (Exception ex)
            {
                mErrorMessage = "Error reading the ScanStats data: " + ex.Message;
            }

            return scanStats;
        }
    }
}
