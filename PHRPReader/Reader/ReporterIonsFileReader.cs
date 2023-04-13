using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader.Data;
using PRISM;

// ReSharper disable UnusedMember.Global

namespace PHRPReader.Reader
{
    /// <summary>
    /// This class reads a tab-delimited _ReporterIons.txt file (created by MASIC)
    /// </summary>
    public class ReporterIonsFileReader : EventNotifier
    {
        // Ignore Spelling: A-Za-z

        /// <summary>
        /// Column headers
        /// </summary>
        /// <remarks>Keys are column name, values are 0-based column index</remarks>
        private readonly SortedDictionary<string, int> mColumnHeaders;

        /// <summary>
        /// Reporter ion specific header info read from the input file
        /// </summary>
        /// <remarks>Keys are 0-based column index, values are column name</remarks>
        private readonly Dictionary<int, ReporterIonColumnInfo> mReporterIonHeaders;

        /// <summary>
        /// Mapping from enum to Reporter Ions file column name
        /// </summary>
        private static readonly Dictionary<ReporterIonsFileColumns, string> mReporterIonsFileColumn = new();

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
        /// Reporter ion header names, starting with the ReporterIonIntensityMax column
        /// </summary>
        /// <remarks>
        /// The MASICResultsMerger program uses this property
        /// </remarks>
        public List<string> ReporterIonHeaderNames { get; }

        /// <summary>
        /// Constructor
        /// </summary>
        public ReporterIonsFileReader()
        {
            mColumnHeaders = new SortedDictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            mReporterIonHeaders = new Dictionary<int, ReporterIonColumnInfo>();
            ReporterIonHeaderNames = new List<string>();
        }

        /// <summary>
        /// Examine the column headers to populate mReporterIonHeaders
        /// </summary>
        /// <param name="columnNames"></param>
        /// <param name="reporterIonMapping"></param>
        private void CacheReporterIonColumnInfo(IReadOnlyList<string> columnNames, IDictionary<string, int> reporterIonMapping)
        {
            // The reporter ion columns should start in the column just after the ReporterIonIntensityMax column
            var reporterIonIntensityMaxName = GetColumnNameByID(ReporterIonsFileColumns.ReporterIonIntensityMax);
            var reporterIonIntensityMaxColIndex = mColumnHeaders[reporterIonIntensityMaxName];

            mReporterIonHeaders.Clear();
            reporterIonMapping.Clear();

            ReporterIonHeaderNames.Clear();
            ReporterIonHeaderNames.Add(reporterIonIntensityMaxName);

            int startColumnIndex;
            if (reporterIonIntensityMaxColIndex > 0)
            {
                startColumnIndex = reporterIonIntensityMaxColIndex + 1;
            }
            else
            {
                OnWarningEvent("ReporterIonIntensityMax column not found in the header line of the _ReporterIons.txt file; this is unexpected");
                startColumnIndex = -1;

                // Find the first column that starts with "Ion_"
                for (var columnIndex = 0; columnIndex < columnNames.Count; columnIndex++)
                {
                    if (!columnNames[columnIndex].StartsWith("Ion_"))
                        continue;

                    startColumnIndex = columnIndex;
                    break;
                }
            }

            if (startColumnIndex <= 0)
            {
                OnWarningEvent("Did not find an Ion_ column in the header line of the _ReporterIons.txt file; this is unexpected");
                return;
            }

            var mzMatcher = new Regex("Ion_(?<MZ>[0-9.]+)(?<Suffix>_[A-Za-z]+)*", RegexOptions.Compiled);

            for (var columnIndex = startColumnIndex; columnIndex < columnNames.Count; columnIndex++)
            {
                ReporterIonHeaderNames.Add(columnNames[columnIndex]);

                var match = mzMatcher.Match(columnNames[columnIndex]);

                if (!match.Success)
                    continue;

                var reporterIonMz = double.Parse(match.Groups["MZ"].Value);
                var columnSuffix = match.Groups["Suffix"].Value;

                if (string.IsNullOrWhiteSpace(columnSuffix))
                {
                    // This is the reporter ion intensity column
                    var reporterIonColumnInfo = new ReporterIonColumnInfo(columnIndex, columnNames[columnIndex], reporterIonMz)
                    {
                        OriginalIntensityColumnIndex = -1,
                        SignalToNoiseColumnIndex = -1,
                        ResolutionColumnIndex = -1
                    };

                    mReporterIonHeaders.Add(columnIndex, reporterIonColumnInfo);

                    reporterIonMapping.Add(columnNames[columnIndex], columnIndex);
                    continue;
                }

                // Find the column in mReporterIonHeaders that we need to update
                // Look for Ion_ReporterIonMz
                var columnToFind = "Ion_" + match.Groups["MZ"].Value;

                if (!reporterIonMapping.TryGetValue(columnToFind, out var reporterIonColumnIndex))
                {
                    OnWarningEvent("Did not find column {0} in reporterIonMapping; this is unexpected", columnToFind);

                    continue;
                }

                switch (columnSuffix)
                {
                    case "_OriginalIntensity":
                        mReporterIonHeaders[reporterIonColumnIndex].OriginalIntensityColumnIndex = columnIndex;
                        break;

                    case "_SignalToNoise":
                        mReporterIonHeaders[reporterIonColumnIndex].SignalToNoiseColumnIndex = columnIndex;
                        break;

                    case "_Resolution":
                        mReporterIonHeaders[reporterIonColumnIndex].ResolutionColumnIndex = columnIndex;
                        break;

                    default:
                        OnWarningEvent("Unrecognized suffix for column {0}", columnNames[columnIndex]);
                        break;
                }
            }
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
        /// Header names and enums for the _ReporterIons.txt file
        /// </summary>
        /// <returns>Dictionary of header names and enum values</returns>
        public static SortedDictionary<string, ReporterIonsFileColumns> GetColumnHeaderNamesAndIDs()
        {
            return new SortedDictionary<string, ReporterIonsFileColumns>(StringComparer.OrdinalIgnoreCase)
            {
                { "Dataset", ReporterIonsFileColumns.Dataset },
                { "ScanNumber", ReporterIonsFileColumns.ScanNumber },
                { "Collision Mode", ReporterIonsFileColumns.CollisionMode },
                { "ParentIonMZ", ReporterIonsFileColumns.ParentIonMZ },
                { "BasePeakIntensity", ReporterIonsFileColumns.BasePeakIntensity },
                { "BasePeakMZ", ReporterIonsFileColumns.BasePeakMZ },
                { "ParentScan", ReporterIonsFileColumns.ParentScan },
                { "ReporterIonIntensityMax", ReporterIonsFileColumns.ReporterIonIntensityMax },
                { "WeightedAvgPctIntensityCorrection", ReporterIonsFileColumns.WeightedAvgPctIntensityCorrection }
            };
        }

        /// <summary>
        /// Compares the names in headerNames to the standard header names tracked by the dictionary returned by GetColumnHeaderNamesAndIDs
        /// Populates a dictionary mapping enum ReporterIonsFileColumns to the 0-based index in columnNames
        /// </summary>
        /// <param name="headerNames"></param>
        /// <returns>Dictionary mapping the enum value to the column index in headerNames (0-based column index)</returns>
        public static Dictionary<ReporterIonsFileColumns, int> GetColumnMapFromHeaderLine(List<string> headerNames)
        {
            var headerColumns = GetColumnHeaderNamesAndIDs();
            return SynFileReaderBaseClass.GetColumnMapFromHeaderLine(headerNames, headerColumns);
        }

        /// <summary>
        /// Get the column name associated with the given enum
        /// </summary>
        /// <param name="column"></param>
        /// <returns>Column name</returns>
        public static string GetColumnNameByID(ReporterIonsFileColumns column)
        {
            if (mReporterIonsFileColumn.Count > 0)
            {
                return mReporterIonsFileColumn[column];
            }

            foreach (var item in GetColumnHeaderNamesAndIDs())
            {
                mReporterIonsFileColumn.Add(item.Value, item.Key);
            }

            return mReporterIonsFileColumn[column];
        }

        private bool GetColumnValue(IReadOnlyList<string> dataColumns, int colIndex, out double value)
        {
            if (colIndex >= 0 && colIndex < dataColumns.Count &&
                double.TryParse(dataColumns[colIndex], out value))
            {
                return true;
            }

            value = 0;
            return false;
        }

        /// <summary>
        /// Open a tab-delimited _ReporterIons.txt file and read the data
        /// </summary>
        /// <param name="inputFilePath">Input file path</param>
        /// <returns>A Dictionary where keys are ScanNumber and values are instances of the ReporterIonInfo class</returns>
        public Dictionary<int, ReporterIonInfo> ReadReporterIonData(string inputFilePath)
        {
            var reporterIonInfoData = new Dictionary<int, ReporterIonInfo>();

            // Keys in this dictionary are reporter ion column names, e.g. Ion_126.128
            // Values are the column index of the reporter ion intensity column
            var reporterIonMapping = new Dictionary<string, int>();

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

                            CacheReporterIonColumnInfo(splitLine, reporterIonMapping);

                            skipLine = true;
                        }

                        headerLineParsed = true;
                    }

                    if (skipLine || splitLine.Length < 4)
                    {
                        continue;
                    }

                    var scanNumber = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.ScanNumber), mColumnHeaders, -1);

                    if (scanNumber < 0 || reporterIonInfoData.ContainsKey(scanNumber))
                        continue;

                    var reporterIonInfo = new ReporterIonInfo(scanNumber)
                    {
                        ParentIonMzText = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.ParentIonMZ), mColumnHeaders),
                        CollisionMode = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.CollisionMode), mColumnHeaders, string.Empty),
                        ParentIonMZ = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.ParentIonMZ), mColumnHeaders, 0.0),
                        BasePeakIntensity = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.BasePeakIntensity), mColumnHeaders, 0.0),
                        BasePeakMZ = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.BasePeakMZ), mColumnHeaders, 0.0),
                        ParentScan = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.ParentScan), mColumnHeaders, 0),
                        ReporterIonIntensityMax = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.ReporterIonIntensityMax), mColumnHeaders, 0.0),
                        WeightedAvgPctIntensityCorrection = ReaderFactory.LookupColumnValue(splitLine, GetColumnNameByID(ReporterIonsFileColumns.WeightedAvgPctIntensityCorrection), mColumnHeaders, 0.0),
                    };

                    var startingColumnIndex = ReaderFactory.LookupColumnIndex(GetColumnNameByID(ReporterIonsFileColumns.ReporterIonIntensityMax), mColumnHeaders);
                    reporterIonInfo.ReporterIonColumnData = startingColumnIndex >= 0
                        ? ReaderFactory.FlattenArray(splitLine, startingColumnIndex)
                        : string.Empty;

                    foreach (var reporterIon in mReporterIonHeaders)
                    {
                        var reporterIonStats = new ReporterIonStats(reporterIon.Value.MZ);

                        // The reporter ion intensity columns should always have a value
                        if (GetColumnValue(splitLine, reporterIon.Value.IntensityColumnIndex, out var reporterIonIntensity))
                        {
                            reporterIonStats.Intensity = reporterIonIntensity;
                        }

                        // The OriginalIntensity columns will only have a value if the MASIC parameter file had ReporterIonSaveUncorrectedIntensities=True
                        if (GetColumnValue(splitLine, reporterIon.Value.OriginalIntensityColumnIndex, out var originalIntensity))
                        {
                            reporterIonStats.OriginalIntensity= originalIntensity;
                        }

                        // The SignalToNoise columns will only have a value if at least one reporter ion had a non-zero intensity
                        if (GetColumnValue(splitLine, reporterIon.Value.SignalToNoiseColumnIndex, out var signalToNoise))
                        {
                            reporterIonStats.SignalToNoise = signalToNoise;
                        }

                        // The Resolution columns will only have a value if at least one reporter ion had a non-zero intensity
                        if (GetColumnValue(splitLine, reporterIon.Value.ResolutionColumnIndex, out var resolution))
                        {
                            reporterIonStats.Resolution = resolution;
                        }

                        reporterIonInfo.ReporterIons.Add(reporterIonStats);
                    }

                    reporterIonInfoData.Add(scanNumber, reporterIonInfo);
                }
            }
            catch (Exception ex)
            {
                ReportError("Error reading the _ReporterIons.txt file", ex);
            }

            return reporterIonInfoData;
        }

        private void ReportError(string message, Exception ex)
        {
            OnErrorEvent(message, ex);
            mErrorMessage = string.Format("{0}{1}", message, ex == null ? string.Empty : ": " + ex.Message);
        }
    }
}
