using System;
using System.Collections.Generic;
using System.IO;
using PHRPReader;
using PHRPReader.Reader;
using PRISM;

namespace PeptideHitResultsProcessor.Processor
{
    // ReSharper disable CommentTypo

    /// <summary>
    /// This class loads precursor m/z values from dataset data files
    /// </summary>
    /// <remarks>
    /// For each dataset name in the list sent to LoadPrecursorIons, look for a valid data file to read
    /// Supported types:
    /// 1) _PrecursorInfo.txt file, which has columns ScanNumber and PrecursorMz or columns ScanNumber and ScanFilterText (or "Scan Filter Text")
    /// 2) MASIC Dataset_SICstats.txt file,    columns FragScanNumber and MZ
    /// 3) MASIC Dataset_ScanStatsEx.txt file, columns ScanNumber     and "Scan Filter Text"
    ///
    /// </remarks>
    // ReSharper restore CommentTypo
    internal class ScanPrecursorMzLoader : EventNotifier
    {
        private string InputDirectoryPath { get; }

        public ScanPrecursorMzLoader(string inputDirectoryPath)
        {
            InputDirectoryPath = inputDirectoryPath;
        }

        /// <summary>
        /// Look for data files that have information regarding the precursor m/z value for each scan, for each dataset
        /// If found, populate a dictionary with the data
        /// </summary>
        /// <param name="datasetNames"></param>
        /// <returns>Dictionary where keys are dataset names and values are dictionaries of precursor m/z by scan number</returns>
        public Dictionary<string, Dictionary<int, double>> LoadPrecursorIons(IEnumerable<string> datasetNames)
        {
            var precursorsByDataset = new Dictionary<string, Dictionary<int, double>>();
            var inputDirectory = new DirectoryInfo(InputDirectoryPath);

            var directoriesToCheck = new List<DirectoryInfo>
            {
                inputDirectory
            };

            foreach (var subDirectory in inputDirectory.GetDirectories())
            {
                directoriesToCheck.Add(subDirectory);
            }

            if (inputDirectory.Parent != null)
            {
                directoriesToCheck.Add(inputDirectory.Parent);

                foreach (var subDirectory in inputDirectory.Parent.GetDirectories())
                {
                    if (subDirectory.FullName.Equals(inputDirectory.FullName))
                        continue;

                    directoriesToCheck.Add(subDirectory);
                }
            }

            foreach (var dataset in datasetNames)
            {
                var precursorsByScan = LoadPrecursorIons(directoriesToCheck, dataset);

                if (precursorsByScan.Count > 0)
                {
                    precursorsByDataset.Add(dataset, precursorsByScan);
                }
            }

            return precursorsByDataset;
        }

        private Dictionary<int, double> LoadPrecursorIons(IEnumerable<DirectoryInfo> directoriesToCheck, string dataset)
        {
            foreach (var directory in directoriesToCheck)
            {
                var precursorInfoFile = new FileInfo(Path.Combine(directory.FullName, dataset + ReaderFactory.PRECURSOR_INFO_FILENAME_SUFFIX));
                var sicStatsFile = new FileInfo(Path.Combine(directory.FullName, dataset + ReaderFactory.SIC_STATS_FILENAME_SUFFIX));
                var extendedScanStatsFile = new FileInfo(Path.Combine(directory.FullName, dataset + ReaderFactory.EXTENDED_SCAN_STATS_FILENAME_SUFFIX));

                if (precursorInfoFile.Exists)
                {
                    OnStatusEvent("Loading precursor m/z information from " + PathUtils.CompactPathString(precursorInfoFile.FullName, 80));
                    return LoadPrecursorsFromPrecursorInfoFile(precursorInfoFile);
                }

                if (sicStatsFile.Exists)
                {
                    OnStatusEvent("Loading precursor m/z information from " + PathUtils.CompactPathString(sicStatsFile.FullName, 80));
                    return LoadPrecursorsFromSICStatsFile(sicStatsFile);
                }

                if (extendedScanStatsFile.Exists)
                {
                    OnStatusEvent("Loading precursor m/z information from " + PathUtils.CompactPathString(extendedScanStatsFile.FullName, 80));
                    return LoadPrecursorsFromExtendedScanStatsFile(extendedScanStatsFile);
                }
            }

            OnStatusEvent("Data file with precursor m/z values not found for dataset " + dataset);
            return new Dictionary<int, double>();
        }

        private Dictionary<int, double> LoadPrecursorsFromPrecursorInfoFile(FileSystemInfo precursorInfoFile)
        {
            var precursorsByScan = new Dictionary<int, double>();

            try
            {
                var reader = new PrecursorInfoFileReader();
                RegisterEvents(reader);

                foreach (var precursor in reader.ReadPrecursorInfoFile(precursorInfoFile.FullName))
                {
                    precursorsByScan.Add(precursor.Key, precursor.Value.PrecursorMz);
                }
            }
            catch (Exception ex)
            {
                OnErrorEvent("Error in LoadPrecursorsFromPrecursorInfoFile", ex);
            }

            return precursorsByScan;
        }

        private Dictionary<int, double> LoadPrecursorsFromSICStatsFile(FileSystemInfo sicStatsFile)
        {
            var precursorsByScan = new Dictionary<int, double>();

            try
            {
                var sicStatsReader = new SICStatsReader();
                var sicStats = sicStatsReader.ReadSICStatsData(sicStatsFile.FullName);

                if (sicStatsReader.ErrorMessage.Length > 0)
                {
                    OnErrorEvent("Error reading SIC Stats data: " + sicStatsReader.ErrorMessage);
                }

                foreach (var item in sicStats)
                {
                    precursorsByScan.Add(item.Value.FragScanNumber, item.Value.MZ);
                }
            }
            catch (Exception ex)
            {
                OnErrorEvent("Error in LoadPrecursorsFromSICStatsFile", ex);
            }

            return precursorsByScan;
        }

        private Dictionary<int, double> LoadPrecursorsFromExtendedScanStatsFile(FileSystemInfo extendedScanStatsFile)
        {
            var precursorsByScan = new Dictionary<int, double>();

            try
            {
                var extendedScanStatsReader = new ExtendedScanStatsReader();
                var scanStats = extendedScanStatsReader.ReadExtendedScanStatsData(extendedScanStatsFile.FullName);

                if (extendedScanStatsReader.ErrorMessage.Length > 0)
                {
                    OnErrorEvent("Error reading Extended ScanStats data: " + extendedScanStatsReader.ErrorMessage);
                }

                foreach (var item in scanStats)
                {
                    if (string.IsNullOrWhiteSpace(item.Value.CollisionMode))
                    {
                        // Likely this is a MS1 spectrum
                        continue;
                    }

                    if (ReaderFactory.ExtractParentIonMzFromFilterText(item.Value.ScanFilterText, out var parentIonMz))
                    {
                        precursorsByScan.Add(item.Key, parentIonMz);
                    }
                }
            }
            catch (Exception ex)
            {
                OnErrorEvent("Error in LoadPrecursorsFromExtendedScanStatsFile", ex);
            }

            return precursorsByScan;
        }
    }
}
