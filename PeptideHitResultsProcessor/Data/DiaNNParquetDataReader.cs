using System;
using System.Collections.Generic;
using ParquetSharp;
using PeptideHitResultsProcessor.Processor;

namespace PeptideHitResultsProcessor.Data
{
    internal class DiaNNParquetDataReader
    {
        // Ignore Spelling: dia, FWHM, normalised, normalisation, peptidoform, proteotypic, PTM

        // ReSharper disable InconsistentNaming

        public long[] RunIndex;                     // Run.Index
        public string[] DatasetName;                // Run
        public string[] Channel;                    // Channel
        public string[] PrecursorId;                // Precursor.Id
        public string[] ModifiedSequence;           // Modified.Sequence
        public string[] StrippedSequence;           // Stripped.Sequence
        public long[] PrecursorCharge;              // Precursor.Charge
        public long[] PrecursorLibIndex;            // Precursor.Lib.Index
        public long[] Decoy;                        // Decoy
        public long[] Proteotypic;                  // Proteotypic
        public float[] PrecursorMz;                 // Precursor.Mz
        public string[] ProteinIDs;                 // Protein.Ids
        public string[] ProteinGroup;               // Protein.Group
        public string[] ProteinNames;               // Protein.Names
        public string[] GeneNames;                  // Genes
        public float[] ElutionTime;                 // RT
        public float[] IndexedRT;                   // iRT
        public float[] PredictedRT;                 // Predicted.RT
        public float[] PredictedIndexedRT;          // Predicted.iRT
        public float[] IonMobility;                 // IM
        public float[] IndexedIonMobility;          // iIM
        public float[] PredictedIonMobility;        // Predicted.IM
        public float[] PredictedIndexedIonMobility; // Predicted.iIM
        public float[] PrecursorQuantity;           // Precursor.Quantity
        public float[] PrecursorNormalized;         // Precursor.Normalised
        public float[] Ms1Area;                     // Ms1.Area
        public float[] Ms1Normalised;               // Ms1.Normalised
        public float[] Ms1ApexArea;                 // Ms1.Apex.Area
        public float[] Ms1ApexMzDelta;              // Ms1.Apex.Mz.Delta
        public float[] NormalisationFactor;         // Normalisation.Factor
        public float[] QuantityQuality;             // Quantity.Quality
        public float[] EmpiricalQuality;            // Empirical.Quality
        public float[] NormalisationNoise;          // Normalisation.Noise
        public float[] Ms1ProfileCorrelation;       // Ms1.Profile.Corr
        public float[] Evidence;                    // Evidence
        public float[] MassEvidence;                // Mass.Evidence
        public float[] ChannelEvidence;             // Channel.Evidence
        public float[] Ms1TotalSignalBefore;        // Ms1.Total.Signal.Before
        public float[] Ms1TotalSignalAfter;         // Ms1.Total.Signal.After
        public float[] RtStart;                     // RT.Start
        public float[] RtStop;                      // RT.Stop
        public float[] FWHM;                        // FWHM
        public float[] ProteinGroupTopN;            // PG.TopN
        public float[] ProteinGroupMaxLFQ;          // PG.MaxLFQ
        public float[] GenesTopN;                   // Genes.TopN
        public float[] GenesMaxLFQ;                 // Genes.MaxLFQ
        public float[] GenesMaxLFQUnique;           // Genes.MaxLFQ.Unique
        public float[] ProteinGroupMaxLFQQuality;   // PG.MaxLFQ.Quality
        public float[] GenesMaxLFQQuality;          // Genes.MaxLFQ.Quality
        public float[] GenesMaxLFQUniqueQuality;    // Genes.MaxLFQ.Unique.Quality
        public float[] QValueDiaNN;                 // Q.Value
        public float[] Pep;                         // PEP
        public float[] GlobalQValue;                // Global.Q.Value
        public float[] LibQValue;                   // Lib.Q.Value
        public float[] PeptidoformQValue;           // Peptidoform.Q.Value
        public float[] GlobalPeptidoformQValue;     // Global.Peptidoform.Q.Value
        public float[] LibPeptidoformQValue;        // Lib.Peptidoform.Q.Value
        public float[] PTMSiteConfidence;           // PTM.Site.Confidence
        public string[] SiteOccupancyProbabilities; // Site.Occupancy.Probabilities
        public string[] ProteinSites;               // Protein.Sites
        public float[] LibPTMSiteConfidence;        // Lib.PTM.Site.Confidence
        public float[] TranslatedQValue;            // Translated.Q.Value
        public float[] ChannelQValue;               // Channel.Q.Value
        public float[] ProteinGroupQValue;          // PG.Q.Value
        public float[] ProteinGroupPEP;             // PG.PEP
        public float[] GeneGroupQValue;             // GG.Q.Value
        public float[] ProteinQValue;               // Protein.Q.Value
        public float[] GlobalProteinGroupQValue;    // Global.PG.Q.Value
        public float[] LibProteinGroupQValue;       // Lib.PG.Q.Value

        // ReSharper restore InconsistentNaming

        private readonly Dictionary<DiaNNResultsProcessor.DiaNNReportFileColumns, int> mColumnMapping;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="columnMapping">Column mapping, where keys are DiaNN column enums and values are the column index in the .parquet file (-1 if not defined)</param>
        public DiaNNParquetDataReader(Dictionary<DiaNNResultsProcessor.DiaNNReportFileColumns, int> columnMapping)
        {
            mColumnMapping = columnMapping;
        }

        private void AppendTsvValue(ICollection<string> outputData, int columnIndex, string dataValue)
        {
            if (columnIndex < 0)
                return;

            outputData.Add(dataValue);
        }

        private void AppendTsvValue(
            ICollection<string> outputData,
            int columnIndex,
            double dataValue,
            byte digitsAfterDecimal = 6)
        {
            if (columnIndex < 0)
                return;

            var thresholdScientific = 1 / Math.Pow(10, digitsAfterDecimal - 2);
            outputData.Add(PRISM.StringUtilities.DblToString(dataValue, digitsAfterDecimal, thresholdScientific));
        }

        private void AppendTsvValue(ICollection<string> outputData, int columnIndex, long dataValue)
        {
            if (columnIndex < 0)
                return;

            outputData.Add(dataValue.ToString());
        }

        public List<string> GetRowData(int rowIndex, string datasetNameToUse)
        {
            var outputData = new List<string>();

            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RunIndex], RunIndex[rowIndex]);                                                     // Run.Index

            outputData.Add(datasetNameToUse);

            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Channel], Channel[rowIndex]);                                                       // Channel
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorId], PrecursorId[rowIndex]);                                               // Precursor.Id
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ModifiedSequence], ModifiedSequence[rowIndex]);                                     // Modified.Sequence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.StrippedSequence], StrippedSequence[rowIndex]);                                     // Stripped.Sequence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorCharge], PrecursorCharge[rowIndex]);                                       // Precursor.Charge
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorLibIndex], PrecursorLibIndex[rowIndex]);                                   // Precursor.Lib.Index
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Decoy], Decoy[rowIndex]);                                                           // Decoy
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Proteotypic], Proteotypic[rowIndex]);                                               // Proteotypic
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorMz], PrecursorMz[rowIndex]);                                               // Precursor.Mz
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinIDs], ProteinIDs[rowIndex]);                                                 // Protein.Ids
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroup], ProteinGroup[rowIndex]);                                             // Protein.Group
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinNames], ProteinNames[rowIndex]);                                             // Protein.Names
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GeneNames], GeneNames[rowIndex]);                                                   // Genes
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RT], ElutionTime[rowIndex], 7);                                                     // RT
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.IndexedRT], IndexedRT[rowIndex], 7);                                                // iRT
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedRT], PredictedRT[rowIndex], 7);                                            // Predicted.RT
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedIndexedRT], PredictedIndexedRT[rowIndex], 7);                              // Predicted.iRT
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.IonMobility], IonMobility[rowIndex], 7);                                            // IM
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.IndexedIonMobility], IndexedIonMobility[rowIndex], 7);                              // iIM
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedIonMobility], PredictedIonMobility[rowIndex], 7);                          // Predicted.IM
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedIndexedIonMobility], PredictedIndexedIonMobility[rowIndex], 7);            // Predicted.iIM
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorQuantity], PrecursorQuantity[rowIndex], 3);                                // Precursor.Quantity
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorNormalized], PrecursorNormalized[rowIndex], 3);                            // Precursor.Normalised
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1Area], Ms1Area[rowIndex], 3);                                                    // Ms1.Area
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1Normalised], Ms1Normalised[rowIndex], 3);                                        // Ms1.Normalised
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1ApexArea], Ms1ApexArea[rowIndex], 3);                                            // Ms1.Apex.Area
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1ApexMzDelta], Ms1ApexMzDelta[rowIndex], 8);                                      // Ms1.Apex.Mz.Delta
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.NormalisationFactor], NormalisationFactor[rowIndex], 4);                            // Normalisation.Factor
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.QuantityQuality], QuantityQuality[rowIndex], 7);                                    // Quantity.Quality
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.EmpiricalQuality], EmpiricalQuality[rowIndex], 7);                                  // Empirical.Quality
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.NormalisationNoise], NormalisationNoise[rowIndex], 7);                              // Normalisation.Noise
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1ProfileCorr], Ms1ProfileCorrelation[rowIndex], 7);                               // Ms1.Profile.Corr
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Evidence], Evidence[rowIndex], 7);                                                  // Evidence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MassEvidence], MassEvidence[rowIndex], 7);                                          // Mass.Evidence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ChannelEvidence], ChannelEvidence[rowIndex], 7);                                    // Channel.Evidence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1TotalSignalBefore], Ms1TotalSignalBefore[rowIndex], 2);                          // Ms1.Total.Signal.Before
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1TotalSignalAfter], Ms1TotalSignalAfter[rowIndex], 2);                            // Ms1.Total.Signal.After
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RTStart], RtStart[rowIndex], 7);                                                    // RT.Start
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RTStop], RtStop[rowIndex], 7);                                                      // RT.Stop
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.FWHM], FWHM[rowIndex], 7);                                                          // FWHM
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupTopN], ProteinGroupTopN[rowIndex], 3);                                  // PG.TopN
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupMaxLFQ], ProteinGroupMaxLFQ[rowIndex], 3);                              // PG.MaxLFQ
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesTopN], GenesTopN[rowIndex], 3);                                                // Genes.TopN
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQ], GenesMaxLFQ[rowIndex], 3);                                            // Genes.MaxLFQ
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQUnique], GenesMaxLFQUnique[rowIndex], 3);                                // Genes.MaxLFQ.Unique
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupMaxLFQQuality], ProteinGroupMaxLFQQuality[rowIndex], 7);                // PG.MaxLFQ.Quality
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQQuality], GenesMaxLFQQuality[rowIndex], 7);                              // Genes.MaxLFQ.Quality
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQUniqueQuality], GenesMaxLFQUniqueQuality[rowIndex], 7);                  // Genes.MaxLFQ.Unique.Quality
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.QValue], QValueDiaNN[rowIndex], 8);                                                 // Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PEP], Pep[rowIndex], 8);                                                            // PEP
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GlobalQValue], GlobalQValue[rowIndex], 8);                                          // Global.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibQValue], LibQValue[rowIndex], 8);                                                // Lib.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PeptidoformQValue], PeptidoformQValue[rowIndex], 8);                                // Peptidoform.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GlobalPeptidoformQValue], GlobalPeptidoformQValue[rowIndex], 8);                    // Global.Peptidoform.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibPeptidoformQValue], LibPeptidoformQValue[rowIndex], 8);                          // Lib.Peptidoform.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PtmSiteConfidence], PTMSiteConfidence[rowIndex], 3);                                // PTM.Site.Confidence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.SiteOccupancyProbabilities], SiteOccupancyProbabilities[rowIndex]);                 // Site.Occupancy.Probabilities
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinSites], ProteinSites[rowIndex]);                                             // Protein.Sites
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibPtmSiteConfidence], LibPTMSiteConfidence[rowIndex], 3);                          // Lib.PTM.Site.Confidence
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.TranslatedQValue], TranslatedQValue[rowIndex]);                                     // Translated.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ChannelQValue], ChannelQValue[rowIndex]);                                           // Channel.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupQValue], ProteinGroupQValue[rowIndex], 8);                              // PG.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupPEP], ProteinGroupPEP[rowIndex], 8);                                    // PG.PEP
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GeneGroupQValue], GeneGroupQValue[rowIndex], 8);                                    // GG.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinQValue], ProteinQValue[rowIndex], 8);                                        // Protein.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GlobalProteinGroupQValue], GlobalProteinGroupQValue[rowIndex], 8);                  // Global.PG.Q.Value
            AppendTsvValue(outputData, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibProteinGroupQValue], LibProteinGroupQValue[rowIndex], 8);                        // Lib.PG.Q.Value

            return outputData;
        }

        private static void ReadColumnBatch(RowGroupReader rowGroupReader, int columnIndex, int startRowIndex, int length, out string[] dataValues)
        {
            dataValues = new string[length];

            if (columnIndex < 0)
                return;

            rowGroupReader.Column(columnIndex).LogicalReader<string>().ReadBatch(dataValues, startRowIndex, length);
        }

        private static void ReadColumnBatchFloat(RowGroupReader rowGroupReader, int columnIndex, int startRowIndex, int length, out float[] dataValues)
        {
            dataValues = new float[length];

            if (columnIndex < 0)
                return;

            rowGroupReader.Column(columnIndex).LogicalReader<float>().ReadBatch(dataValues, startRowIndex, length);
        }

        private static void ReadColumnBatchLong(RowGroupReader rowGroupReader, int columnIndex, int startRowIndex, int length, out long[] dataValues)
        {
            dataValues = new long[length];

            if (columnIndex < 0)
                return;

            rowGroupReader.Column(columnIndex).LogicalReader<long>().ReadBatch(dataValues, startRowIndex, length);
        }

        /// <summary>
        /// Read the data from the .parquet file and store in local arrays
        /// </summary>
        /// <param name="rowGroupReader">.parquet row group reader</param>
        /// <param name="rowCount">Row count in the row group</param>
        public void ReadParquetData(RowGroupReader rowGroupReader, int rowCount)
        {
            ReadColumnBatchLong(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RunIndex], 0, rowCount, out RunIndex);                                         // Run.Index
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.DatasetName], 0, rowCount, out DatasetName);                                       // Run
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Channel], 0, rowCount, out Channel);                                               // Channel
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorId], 0, rowCount, out PrecursorId);                                       // Precursor.Id
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ModifiedSequence], 0, rowCount, out ModifiedSequence);                             // Modified.Sequence
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.StrippedSequence], 0, rowCount, out StrippedSequence);                             // Stripped.Sequence
            ReadColumnBatchLong(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorCharge], 0, rowCount, out PrecursorCharge);                           // Precursor.Charge
            ReadColumnBatchLong(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorLibIndex], 0, rowCount, out PrecursorLibIndex);                       // Precursor.Lib.Index
            ReadColumnBatchLong(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Decoy], 0, rowCount, out Decoy);                                               // Decoy
            ReadColumnBatchLong(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Proteotypic], 0, rowCount, out Proteotypic);                                   // Proteotypic
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorMz], 0, rowCount, out PrecursorMz);                                  // Precursor.Mz
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinIDs], 0, rowCount, out ProteinIDs);                                         // Protein.Ids
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroup], 0, rowCount, out ProteinGroup);                                     // Protein.Group
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinNames], 0, rowCount, out ProteinNames);                                     // Protein.Names
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GeneNames], 0, rowCount, out GeneNames);                                           // Genes
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RT], 0, rowCount, out ElutionTime);                                           // RT
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.IndexedRT], 0, rowCount, out IndexedRT);                                      // iRT
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedRT], 0, rowCount, out PredictedRT);                                  // Predicted.RT
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedIndexedRT], 0, rowCount, out PredictedIndexedRT);                    // Predicted.iRT
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.IonMobility], 0, rowCount, out IonMobility);                                  // IM
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.IndexedIonMobility], 0, rowCount, out IndexedIonMobility);                    // iIM
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedIonMobility], 0, rowCount, out PredictedIonMobility);                // Predicted.IM
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PredictedIndexedIonMobility], 0, rowCount, out PredictedIndexedIonMobility);  // Predicted.iIM
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorQuantity], 0, rowCount, out PrecursorQuantity);                      // Precursor.Quantity
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PrecursorNormalized], 0, rowCount, out PrecursorNormalized);                  // Precursor.Normalised
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1Area], 0, rowCount, out Ms1Area);                                          // Ms1.Area
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1Normalised], 0, rowCount, out Ms1Normalised);                              // Ms1.Normalised
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1ApexArea], 0, rowCount, out Ms1ApexArea);                                  // Ms1.Apex.Area
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1ApexMzDelta], 0, rowCount, out Ms1ApexMzDelta);                            // Ms1.Apex.Mz.Delta
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.NormalisationFactor], 0, rowCount, out NormalisationFactor);                  // Normalisation.Factor
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.QuantityQuality], 0, rowCount, out QuantityQuality);                          // Quantity.Quality
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.EmpiricalQuality], 0, rowCount, out EmpiricalQuality);                        // Empirical.Quality
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.NormalisationNoise], 0, rowCount, out NormalisationNoise);                    // Normalisation.Noise
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1ProfileCorr], 0, rowCount, out Ms1ProfileCorrelation);                     // Ms1.Profile.Corr
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.Evidence], 0, rowCount, out Evidence);                                        // Evidence
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MassEvidence], 0, rowCount, out MassEvidence);                                // Mass.Evidence
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ChannelEvidence], 0, rowCount, out ChannelEvidence);                          // Channel.Evidence
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1TotalSignalBefore], 0, rowCount, out Ms1TotalSignalBefore);                // Ms1.Total.Signal.Before
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.MS1TotalSignalAfter], 0, rowCount, out Ms1TotalSignalAfter);                  // Ms1.Total.Signal.After
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RTStart], 0, rowCount, out RtStart);                                          // RT.Start
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.RTStop], 0, rowCount, out RtStop);                                            // RT.Stop
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.FWHM], 0, rowCount, out FWHM);                                                // FWHM
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupTopN], 0, rowCount, out ProteinGroupTopN);                        // PG.TopN
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupMaxLFQ], 0, rowCount, out ProteinGroupMaxLFQ);                    // PG.MaxLFQ
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesTopN], 0, rowCount, out GenesTopN);                                      // Genes.TopN
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQ], 0, rowCount, out GenesMaxLFQ);                                  // Genes.MaxLFQ
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQUnique], 0, rowCount, out GenesMaxLFQUnique);                      // Genes.MaxLFQ.Unique
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupMaxLFQQuality], 0, rowCount, out ProteinGroupMaxLFQQuality);      // PG.MaxLFQ.Quality
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQQuality], 0, rowCount, out GenesMaxLFQQuality);                    // Genes.MaxLFQ.Quality
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GenesMaxLFQUniqueQuality], 0, rowCount, out GenesMaxLFQUniqueQuality);        // Genes.MaxLFQ.Unique.Quality
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.QValue], 0, rowCount, out QValueDiaNN);                                       // Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PEP], 0, rowCount, out Pep);                                                  // PEP
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GlobalQValue], 0, rowCount, out GlobalQValue);                                // Global.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibQValue], 0, rowCount, out LibQValue);                                      // Lib.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PeptidoformQValue], 0, rowCount, out PeptidoformQValue);                      // Peptidoform.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GlobalPeptidoformQValue], 0, rowCount, out GlobalPeptidoformQValue);          // Global.Peptidoform.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibPeptidoformQValue], 0, rowCount, out LibPeptidoformQValue);                // Lib.Peptidoform.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.PtmSiteConfidence], 0, rowCount, out PTMSiteConfidence);                      // PTM.Site.Confidence
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.SiteOccupancyProbabilities], 0, rowCount, out SiteOccupancyProbabilities);         // Site.Occupancy.Probabilities
            ReadColumnBatch(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinSites], 0, rowCount, out ProteinSites);                                     // Protein.Sites
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibPtmSiteConfidence], 0, rowCount, out LibPTMSiteConfidence);                // Lib.PTM.Site.Confidence
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.TranslatedQValue], 0, rowCount, out TranslatedQValue);                        // Translated.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ChannelQValue], 0, rowCount, out ChannelQValue);                              // Channel.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupQValue], 0, rowCount, out ProteinGroupQValue);                    // PG.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinGroupPEP], 0, rowCount, out ProteinGroupPEP);                          // PG.PEP
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GeneGroupQValue], 0, rowCount, out GeneGroupQValue);                          // GG.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.ProteinQValue], 0, rowCount, out ProteinQValue);                              // Protein.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.GlobalProteinGroupQValue], 0, rowCount, out GlobalProteinGroupQValue);        // Global.PG.Q.Value
            ReadColumnBatchFloat(rowGroupReader, mColumnMapping[DiaNNResultsProcessor.DiaNNReportFileColumns.LibProteinGroupQValue], 0, rowCount, out LibProteinGroupQValue);              // Lib.PG.Q.Value

            // The following columns are in the .tsv file created by DIA-NN prior to v2.0, but are not in .parquet files
            // Float:  DiaNNReportFileColumns.ProteinGroupQuantity
            // Float:  DiaNNReportFileColumns.ProteinGroupNormalized
            // Float:  DiaNNReportFileColumns.GenesQuantity
            // Float:  DiaNNReportFileColumns.GenesNormalized
            // Float:  DiaNNReportFileColumns.PrecursorTranslated
            // Float:  DiaNNReportFileColumns.TranslatedQuality
            // Float:  DiaNNReportFileColumns.MS1Translated
            // String: DiaNNReportFileColumns.FirstProteinDescription
            // Float:  DiaNNReportFileColumns.SpectrumSimilarity
            // Float:  DiaNNReportFileColumns.Averagine
            // Float:  DiaNNReportFileColumns.CScore
            // Float:  DiaNNReportFileColumns.DecoyEvidence
            // Float:  DiaNNReportFileColumns.DecoyCScore
            // Float:  DiaNNReportFileColumns.FragmentQuantRaw
            // Float:  DiaNNReportFileColumns.FragmentQuantCorrected
            // Float:  DiaNNReportFileColumns.FragmentCorrelations
            // Int:    DiaNNReportFileColumns.MS2Scan
        }
    }
}
