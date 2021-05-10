//*********************************************************************************************************
// Written by Matthew Monroe for the US Department of Energy
// Pacific Northwest National Laboratory, Richland, WA
//
// Created 05/06/2021
//
// This class reads MASIC SIC Stats data from a tab-delimited _SICStats.txt file
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
    /// This class reads MASIC SIC stats data from a tab-delimited _SICStats.txt file
    /// </summary>
    public class SICStatsReader : EventNotifier
    {
        /// <summary>
        /// Column headers
        /// </summary>
        /// <remarks>Keys are column name, values are 0-based column index</remarks>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Mapping from enum to Scan Stats file column name
        /// </summary>
        private static readonly Dictionary<SICStatsFileColumns, string> mSICStatsFileColumn = new();

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
        public SICStatsReader()
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
        /// Header names and enums for the _SICStats.txt file
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, SICStatsFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new(StringComparer.OrdinalIgnoreCase)
            {
                { "Dataset", SICStatsFileColumns.Dataset },
                { "ParentIonIndex", SICStatsFileColumns.ParentIonIndex },
                { "MZ", SICStatsFileColumns.MZ },
                { "SurveyScanNumber", SICStatsFileColumns.SurveyScanNumber },
                { "FragScanNumber", SICStatsFileColumns.FragScanNumber },
                { "OptimalPeakApexScanNumber", SICStatsFileColumns.OptimalPeakApexScanNumber },
                { "PeakApexOverrideParentIonIndex", SICStatsFileColumns.PeakApexOverrideParentIonIndex },
                { "CustomSICPeak", SICStatsFileColumns.CustomSICPeak },
                { "PeakScanStart", SICStatsFileColumns.PeakScanStart },
                { "PeakScanEnd", SICStatsFileColumns.PeakScanEnd },
                { "PeakScanMaxIntensity", SICStatsFileColumns.PeakScanMaxIntensity },
                { "PeakMaxIntensity", SICStatsFileColumns.PeakMaxIntensity },
                { "PeakSignalToNoiseRatio", SICStatsFileColumns.PeakSignalToNoiseRatio },
                { "FWHMInScans", SICStatsFileColumns.FWHMInScans },
                { "PeakArea", SICStatsFileColumns.PeakArea },
                { "ParentIonIntensity", SICStatsFileColumns.ParentIonIntensity },
                { "PeakBaselineNoiseLevel", SICStatsFileColumns.PeakBaselineNoiseLevel },
                { "PeakBaselineNoiseStDev", SICStatsFileColumns.PeakBaselineNoiseStDev },
                { "PeakBaselinePointsUsed", SICStatsFileColumns.PeakBaselinePointsUsed },
                { "StatMomentsArea", SICStatsFileColumns.StatMomentsArea },
                { "CenterOfMassScan", SICStatsFileColumns.CenterOfMassScan },
                { "PeakStDev", SICStatsFileColumns.PeakStDev },
                { "PeakSkew", SICStatsFileColumns.PeakSkew },
                { "PeakKSStat", SICStatsFileColumns.PeakKSStat },
                { "StatMomentsDataCountUsed", SICStatsFileColumns.StatMomentsDataCountUsed }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum SICStatsFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        // ReSharper disable once UnusedMember.Global
        public static Dictionary<SICStatsFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return SynFileReaderBaseClass.GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(SICStatsFileColumns column)
        {
            if (mSICStatsFileColumn.Count > 0)
            {
                return mSICStatsFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mSICStatsFileColumn.Add(item.Value, item.Key);
            }

            return mSICStatsFileColumn[column];
        }

        /// <summary>
        /// Open a tab-delimited _SICStats.txt file and read the data
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are instances of SICStatsInfo</returns>
        public Dictionary<int, SICStatsInfo> ReadSICStatsData(string inputFilePath)
        {
            var sicStats = new Dictionary<int, SICStatsInfo>();

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

                    var parentIonIndex = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.ParentIonIndex), mColumnHeaders, -1);

                    if (parentIonIndex < 0 || sicStats.ContainsKey(parentIonIndex))
                        continue;

                    var sicStatsInfo = new SICStatsInfo(parentIonIndex)
                    {
                        MzText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.MZ), mColumnHeaders),
                        PeakMaxIntensityText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakMaxIntensity), mColumnHeaders),
                        PeakSignalToNoiseRatioText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakSignalToNoiseRatio), mColumnHeaders),
                        PeakAreaText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakArea), mColumnHeaders),
                        ParentIonIntensityText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.ParentIonIntensity), mColumnHeaders),
                        PeakBaselineNoiseLevelText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakBaselineNoiseLevel), mColumnHeaders),
                        PeakBaselineNoiseStDevText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakBaselineNoiseStDev), mColumnHeaders),
                        StatMomentsAreaText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.StatMomentsArea), mColumnHeaders),
                        PeakStDevText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakStDev), mColumnHeaders),
                        PeakSkewText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakSkew), mColumnHeaders),
                        PeakKSStatText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakKSStat), mColumnHeaders),

                        Dataset = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.Dataset), mColumnHeaders, 0),
                        MZ = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.MZ), mColumnHeaders, 0.0),
                        SurveyScanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.SurveyScanNumber), mColumnHeaders, 0),
                        FragScanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.FragScanNumber), mColumnHeaders, 0),
                        OptimalPeakApexScanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.OptimalPeakApexScanNumber), mColumnHeaders, 0),
                        PeakApexOverrideParentIonIndex = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakApexOverrideParentIonIndex), mColumnHeaders, 0),
                        CustomSICPeak = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.CustomSICPeak), mColumnHeaders, 0),
                        PeakScanStart = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakScanStart), mColumnHeaders, 0),
                        PeakScanEnd = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakScanEnd), mColumnHeaders, 0),
                        PeakScanMaxIntensity = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakScanMaxIntensity), mColumnHeaders, 0),
                        PeakMaxIntensity = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakMaxIntensity), mColumnHeaders, 0.0),
                        PeakSignalToNoiseRatio = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakSignalToNoiseRatio), mColumnHeaders, 0.0),
                        FWHMInScans = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.FWHMInScans), mColumnHeaders, 0),
                        PeakArea = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakArea), mColumnHeaders, 0.0),
                        ParentIonIntensity = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.ParentIonIntensity), mColumnHeaders, 0.0),
                        PeakBaselineNoiseLevel = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakBaselineNoiseLevel), mColumnHeaders, 0.0),
                        PeakBaselineNoiseStDev = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakBaselineNoiseStDev), mColumnHeaders, 0.0),
                        PeakBaselinePointsUsed = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakBaselinePointsUsed), mColumnHeaders, 0),
                        StatMomentsArea = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.StatMomentsArea), mColumnHeaders, 0.0),
                        CenterOfMassScan = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.CenterOfMassScan), mColumnHeaders, 0),
                        PeakStDev = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakStDev), mColumnHeaders, 0.0),
                        PeakSkew = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakSkew), mColumnHeaders, 0.0),
                        PeakKSStat = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.PeakKSStat), mColumnHeaders, 0.0),
                        StatMomentsDataCountUsed = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(SICStatsFileColumns.StatMomentsDataCountUsed), mColumnHeaders, 0),
                    };

                    sicStats.Add(parentIonIndex, sicStatsInfo);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error reading the SICStats data", ex);
            }

            return sicStats;
        }

        private void ReportError(string message, Exception ex)
        {
            OnErrorEvent(message, ex);
            mErrorMessage = string.Format("{0}{1}", message, ex == null ? string.Empty : ": " + ex.Message);
        }
    }
}
