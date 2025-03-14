﻿using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework;
using PHRPReader;
using PHRPReader.Data;

// ReSharper disable StringLiteralTypo

namespace PHRP_UnitTests
{
    [TestFixture]
    public class PHRPReaderTests
    {
        // Ignore Spelling: Da, Fragger, gb, gi, Hyperscore, PHRP, Proteoform, Quant

        private struct MSGFPlusInfo
        {
            public string Peptide;
            public string Protein;
            public double SpecEValue;
            public double EValue;
            public double QValue;
            public double PepQValue;
        }

        private struct TopPICInfo
        {
            public string Peptide;
            public string Protein;
            public double? FeatureIntensity;
            public double? FeatureScore;
            public double? EValue;
            public double? QValue;
            public double? ProteoformQValue;
        }

        [Test]
        [TestCase(
            @"MXQ202103181341_Auto1878805\QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14_maxq_syn.txt",
            0,
            "K.QSEWQAYLDWAVNAFK.L|R.TAMNFIQTLSGVATLTK.H|K.FNQIGSLTETLAAIR.M|K.NSASIANSIFEVIQNIVAR.Y|K.KSAQSALYQLYR.N",
            "1.26E-73|2.09E-47|9.97E-47|1.38E-40|2.06E-35",
            "12.36589|4.01741|13.47918|5.69358|6.06656")]
        public void TestMaxQuantReader(
            string inputFilePath,
            int expectedResultCount,
            string initialPeptides,
            string initialPEPValues,
            string initialDelMPPM,
            int countToDisplay = 5)
        {
            var inputFile = FindFile(inputFilePath);

            var options = new StartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var expectedPeptides = initialPeptides.Split('|').ToList();

            var expectedPEPValues = SplitDelimitedDoubles(initialPEPValues);
            var expectedDelMPPM = SplitDelimitedDoubles(initialDelMPPM);

            var reader = new ReaderFactory(inputFile.FullName, options);
            var resultCount = 0;
            var peptideMatchCount = 0;
            var pepMatchCount = 0;
            var delMMatchCount = 0;

            Console.WriteLine();
            Console.WriteLine("{0,-50} {1,-14} {2,-14}", "Peptide", "PEP", "DelM_PPM");

            while (reader.MoveNext())
            {
                var peptide = reader.CurrentPSM.Peptide;
                var posteriorErrorProbability = reader.CurrentPSM.GetScore("PEP");
                var delMPPM = reader.CurrentPSM.MassErrorPPM;

                if (resultCount < countToDisplay)
                {
                    Console.WriteLine("{0,-50} {1,-14} {2,-14}", peptide, posteriorErrorProbability, delMPPM);
                }

                if (resultCount < expectedPeptides.Count && expectedPeptides[resultCount].Length > 0)
                {
                    Assert.AreEqual(expectedPeptides[resultCount], peptide);
                    peptideMatchCount++;
                }

                if (resultCount < expectedPEPValues.Count)
                {
                    AssertValuesMatch(expectedPEPValues[resultCount], posteriorErrorProbability, "PEP", resultCount + 2);
                    pepMatchCount++;
                }

                if (resultCount < expectedDelMPPM.Count)
                {
                    AssertValuesMatch(expectedDelMPPM[resultCount], delMPPM, "DelM PPM", resultCount + 2);
                    delMMatchCount++;
                }

                resultCount++;
            }

            Console.WriteLine();
            Console.WriteLine("Read {0:N0} results", resultCount);
            Console.WriteLine("{0} peptides matched expected values", peptideMatchCount);
            Console.WriteLine("{0} SpecEValues matched expected values", pepMatchCount);
            Console.WriteLine("{0} DelMPPM values expected values", delMMatchCount);

            if (expectedResultCount > 0)
                Assert.AreEqual(expectedResultCount, resultCount);
        }
        [Test]
        [TestCase(
            @"MSF202008051507_Auto1822588\MCF7Cell_O-GlcNAc_Petyuk_R1_15May20_Rage_Rep-20-05-01_msfragger_syn.txt",
            0,
            "R.KDLYANTVLSGGTTMYPGIADR.M|K.TPVEPEVAIHR.I|K.HFVALSTNTTK.V|R.IFGLLMGTLQK.F|R.VSLDVNHFAPDELTVK.T|K.TSFFQALGITTK.I|K.YGLIYHASLVGQTSPK.H|K.TVAGQDAVIVLLGTR.N|K.VPDGMVGFIIGR.G|K.QSLGELIGTLNAAK.V|R.HQGVMVGMGQK.D|R.IFGLLMGTLQK.F|K.TGVAVNKPAEFTVDAK.H|R.TLMNLGGLAVAR.D|R.VALTGLTVAEYFR.D|K.SRLEQEIATYR.S",
            "40.103|25.155|36.05|21.203|31.232|24.058|31.621|28.151|23.075|30.579|30.868|22.124|30.701|25.411|33.862|21.67|22.413|19.873|28.828|28.065|23.243|24.455|29.676|25.016|19.285|28.911|32.922|22.52|30.016|18.196|19.231|23.43|24.037|19.069|25.025|22.485|21.927|23.232|23.932|31.947|21.399|27.546|28.499|24.255",
            "1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|1|0.9999|1|1|0.9997|1|1|1|0.9956|0.9999|0.9999|0.9963|0.9999|0.9988|0.9929|0.9985",
            "-1.19424|-1.72579|-1.74787|-1.82725|-2.38748|-3.49777|-0.03762|-0.67151|-0.95544|-0.86517|-1.22818|-1.50355|0.80006|-1.71956|-1.04092|-1.0164"
            , 135)]
        public void TestMSFraggerReader(
            string inputFilePath,
            int expectedResultCount,
            string initialPeptides,
            string initialHyperscoreValues,
            string initialPeptideProphetValues,
            string initialDelMPPM,
            int countToDisplay = 5)
        {
            var inputFile = FindFile(inputFilePath);

            var options = new StartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var expectedPeptides = initialPeptides.Split('|').ToList();

            var hyperscoreValues = SplitDelimitedDoubles(initialHyperscoreValues);
            var peptideProphetValues = SplitDelimitedDoubles(initialPeptideProphetValues);
            var expectedDelMPPM = SplitDelimitedDoubles(initialDelMPPM);

            var reader = new ReaderFactory(inputFile.FullName, options);
            var resultCount = 0;
            var peptideMatchCount = 0;
            var hyperscoreMatchCount = 0;
            var peptideProphetMatchCount = 0;
            var delMMatchCount = 0;

            Console.WriteLine();
            Console.WriteLine("{0,-50} {1,-14} {2,-14} {3,-14}", "Peptide", "Hyperscore", "PepProphetProb", "DelM_PPM");

            while (reader.MoveNext())
            {
                var peptide = reader.CurrentPSM.Peptide;
                var hyperscore = reader.CurrentPSM.GetScore("Hyperscore");
                var peptideProphetProbability = reader.CurrentPSM.GetScore("PeptideProphetProbability");
                var delMPPM = reader.CurrentPSM.MassErrorPPM;

                if (resultCount < countToDisplay)
                {
                    Console.WriteLine("{0,-50} {1,-14} {2,-14} {3,-14}", peptide, hyperscore, peptideProphetProbability, delMPPM);
                }

                if (resultCount < expectedPeptides.Count && expectedPeptides[resultCount].Length > 0)
                {
                    Assert.AreEqual(expectedPeptides[resultCount], peptide);
                    peptideMatchCount++;
                }

                if (resultCount < hyperscoreValues.Count)
                {
                    AssertValuesMatch(hyperscoreValues[resultCount], hyperscore, "Hyperscore", resultCount + 2);
                    hyperscoreMatchCount++;
                }

                if (resultCount < peptideProphetValues.Count)
                {
                    AssertValuesMatch(peptideProphetValues[resultCount], peptideProphetProbability, "PeptideProphetProbability", resultCount + 2);
                    peptideProphetMatchCount++;
                }

                if (resultCount < expectedDelMPPM.Count)
                {
                    AssertValuesMatch(expectedDelMPPM[resultCount], delMPPM, "DelM PPM", resultCount + 2);
                    delMMatchCount++;
                }

                resultCount++;
            }

            Console.WriteLine();
            Console.WriteLine("Read {0:N0} results", resultCount);
            Console.WriteLine("{0} peptides matched expected values", peptideMatchCount);
            Console.WriteLine("{0} hyperscore values matched expected values", hyperscoreMatchCount);
            Console.WriteLine("{0} peptide prophet probabilities matched expected values", peptideProphetMatchCount);
            Console.WriteLine("{0} DelMPPM values expected values", delMMatchCount);

            if (expectedResultCount > 0)
                Assert.AreEqual(expectedResultCount, resultCount);
        }

        [Test]
        [TestCase(
            @"MSG201205081554_Auto838496\QC_Shew_11_06_Run-03_7May12_Roc_12-04-08_msgfplus_syn.txt",
            5536,
            "K.EQLEDNMVLGTMMLAQDEVDGIVSGAVNTTANTIRPPLQLIK.T|K.INVIGGHSGVTILPLLSQVEGVTFSDEEVASLTKR.I|R.AASVAPLSEVLAVNKEPLVSIDFNHNAFSSNFDATQTR.V|K.EQLEDNMVLGTMMLAQDEVDGIVSGAVNTTANTIRPPLQLIK.T|K.SHIETNGGILNGTSAADVQTR.A|R.DKSITSSTALWTPNLYIAADNPLNPTFGTANASEVQAYHR.L",
            "7.216E-28|3.403E-21|1.725E-20|2.061E-20|4.719E-20|1.982E-19",
            "7.23093|5.39579|5.89445|7.12301|3.87878|6.75958")]
        [TestCase(
            @"MSG201802011337_Auto1547784\MCF-7_pSTY_2_OldGr_15Mar17_Merry__IntEmtr_msgfplus_syn.txt",
            11853,
            "R.AS*PAPGSGHPEGPGAHLDMNSLDR.A|R.S*QSPAASDCSSSSSSASLPSSGR.S|R.SQS*PAASDCSSSSSSASLPSSGR.S|R.SSGS*PYGGGYGSGGGSGGYGSR.R|R.S*SGSPYGGGYGSGGGSGGYGSR.R|R.SS*GSPYGGGYGSGGGSGGYGSR.R",
            "3.8745E-30|1.918E-27|1.2313E-26|2.5021E-26|7.2866E-26|7.2866E-26",
            "-1.5682|-0.107|-0.107|-2.6967|-2.6967|-2.6967")]
        [TestCase(
            @"MSG202011150906_Auto1850387\QC_Mam_19_01_R3_12Nov20_Oak_Jup-20-10-01_msgfplus_syn.txt",
            37446,
            "R.HIADLAGNPEVILPVPAFNVINGGSHAGNK.L|R.LNVFKNDQDTWDYTNPNLSGQGDPGSNPNKR.Q|K.SCEAGYSPSYKEDKHFGYTSYSVSNSVK.E|K.SCEAGYSPSYKEDKHFGYTSYSVSNSVK.E|K.SGDAAIVDMVPGKPMCVESFSDYPPLGR.F|R.IEVQDSSGGTTALRPSASTQALSSSVSSSK.L",
            "3.89E-36|3.89E-35|1.28E-34|1.66E-34|2.61E-34|4.36E-34",
            "5.81363|-1.18764|-1.22821|-3.12744|5.4979|-5.33381")]
        public void TestMSGFPlusReader(
            string inputFilePath,
            int expectedResultCount,
            string initialPeptides,
            string initialSpecEValues,
            string initialDelMPPM,
            int countToDisplay = 5)
        {
            var inputFile = FindFile(inputFilePath);

            var options = new StartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var expectedPeptides = initialPeptides.Split('|').ToList();

            var expectedSpecEValues = SplitDelimitedDoubles(initialSpecEValues);
            var expectedDelMPPM = SplitDelimitedDoubles(initialDelMPPM);

            var reader = new ReaderFactory(inputFile.FullName, options);
            var resultCount = 0;
            var peptideMatchCount = 0;
            var specEValueMatchCount = 0;
            var delMMatchCount = 0;

            Console.WriteLine();
            Console.WriteLine("{0,-50} {1,-14} {2,-14}", "Peptide", "SpecEvalue", "DelM_PPM");

            while (reader.MoveNext())
            {
                var peptide = reader.CurrentPSM.Peptide;
                var specEvalue = reader.CurrentPSM.MSGFSpecEValue;
                var delMPPM = reader.CurrentPSM.MassErrorPPM;

                if (resultCount < countToDisplay)
                {
                    Console.WriteLine("{0,-50} {1,-14} {2,-14}", peptide, specEvalue, delMPPM);
                }

                if (resultCount < expectedPeptides.Count && expectedPeptides[resultCount].Length > 0)
                {
                    Assert.AreEqual(expectedPeptides[resultCount], peptide);
                    peptideMatchCount++;
                }

                if (resultCount < expectedSpecEValues.Count)
                {
                    AssertValuesMatch(expectedSpecEValues[resultCount], specEvalue, "SpecEValue", resultCount + 2);
                    specEValueMatchCount++;
                }

                if (resultCount < expectedDelMPPM.Count)
                {
                    AssertValuesMatch(expectedDelMPPM[resultCount], delMPPM, "DelM PPM", resultCount + 2);
                    delMMatchCount++;
                }

                resultCount++;
            }

            Console.WriteLine();
            Console.WriteLine("Read {0:N0} results", resultCount);
            Console.WriteLine("{0} peptides matched expected values", peptideMatchCount);
            Console.WriteLine("{0} SpecEValues matched expected values", specEValueMatchCount);
            Console.WriteLine("{0} DelMPPM values expected values", delMMatchCount);

            if (expectedResultCount > 0)
                Assert.AreEqual(expectedResultCount, resultCount);
        }

        [Test]
        [TestCase(
            @"MSG201205081554_Auto838496\QC_Shew_11_06_Run-03_7May12_Roc_12-04-08_msgfplus_syn.txt",
            "50|200|3000",
            "K.GQYIAASYETNPSTDTPGNR.Y|K.VTSGSTPEVAEYVDQLYK.S|K.VDVPAIAR.Q",
            "SO_2583|SO_3190|SO_2301",
            "2.69E-22,3.75E-16,1.43E-14,1|4.09E-19,5.70E-13,6.88E-12,1|1.22E-09,0.00168205,0.00160812,1")]
        [TestCase(
            @"MSG201802011337_Auto1547784\MCF-7_pSTY_2_OldGr_15Mar17_Merry__IntEmtr_msgfplus_syn.txt",
            "75|600|4000|9000",
            "R.RT*SMGGTQQQFVEGVR.M|R.NQGGYGGSSSS*SSYGSGR.R|R.SLEPEPQQS*LEDGSPAKGEPSQAWR.E|K.RPSPS*PPNDCR.L",
            "CTNB1_HUMAN|ROA1_HUMAN|HID1_HUMAN|XXX_S31C1_HUMAN",
            "6.58E-20,1.43E-12,0,0|1.55E-15,3.37E-08,0,0|1.60E-09,0.035132,0.00366,0|1.91E-07,4.0792,0.2228,0.27419")]
        [TestCase(
            @"MSG202011150906_Auto1850387\QC_Mam_19_01_R3_12Nov20_Oak_Jup-20-10-01_msgfplus_syn.txt",
            "100|19999|29999",
            "K.KIEPELEGSSAVTSHDSSTNGLISFIK.Q|K.M*KGDYYR.Y|K.VERLVK.L",
            "G6PI_MOUSE|1433F_MOUSE|XXX_NFRKB_MOUSE",
            "1.47E-28,2.71E-21,0,0|4.85E-11,8.56E-04,0.000056754,0.000074867|2.95E-08,0.4992,0.035014,0.045267")]
        public void TestMSGFPlusScores(
            string inputFilePath,
            string resultIdList,
            string peptidesByResultID,
            string proteinsByResultID,
            string scoresByResultId      // MSGFDB_SpecEValue, EValue, QValue, PepQValue
            )
        {
            var inputFile = FindFile(inputFilePath);

            var options = new StartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var resultIDs = resultIdList.Split('|').ToList();

            var scores = scoresByResultId.Split('|').ToList();
            var peptides = peptidesByResultID.Split('|').ToList();
            var proteins = proteinsByResultID.Split('|').ToList();

            var infoByResultID = new Dictionary<int, MSGFPlusInfo>();
            var matchedResultIDs = new SortedSet<int>();

            for (var i = 0; i < resultIDs.Count; i++)
            {
                var scoreValues = SplitDelimitedDoubles(scores[i], ',');

                var msgfPlusInfo = new MSGFPlusInfo
                {
                    Peptide = peptides[i],
                    Protein = proteins[i],
                    SpecEValue = GetValueOrDefault(scoreValues[0]),
                    EValue = GetValueOrDefault(scoreValues[1]),
                    QValue = GetValueOrDefault(scoreValues[2]),
                    PepQValue = GetValueOrDefault(scoreValues[3])
                };

                var resultId = int.Parse(resultIDs[i]);
                infoByResultID.Add(resultId, msgfPlusInfo);
            }

            var reader = new ReaderFactory(inputFile.FullName, options);
            var resultCount = 0;
            var resultIDsExamined = 0;

            while (reader.MoveNext())
            {
                resultCount++;

                if (!infoByResultID.TryGetValue(reader.CurrentPSM.ResultID, out var msgfplusInfo))
                    continue;

                matchedResultIDs.Add(reader.CurrentPSM.ResultID);

                Assert.AreEqual(msgfplusInfo.Peptide, reader.CurrentPSM.Peptide);
                Assert.AreEqual(msgfplusInfo.Protein, reader.CurrentPSM.ProteinFirst);

                var evalueOrPvalue = GetScore(reader.CurrentPSM, "EValue", "PValue");
                var qValueOrFDR = GetScore(reader.CurrentPSM, "QValue", "EFDR");
                var pepQValueOrFDR = GetScore(reader.CurrentPSM, "PepQValue", "PepFDR");

                AssertValuesMatch(msgfplusInfo.SpecEValue, reader.CurrentPSM.MSGFSpecEValue, "SpecEValue", resultCount + 1);
                AssertValuesMatch(msgfplusInfo.EValue, evalueOrPvalue, "EValue (or PValue)", resultCount + 1);
                AssertValuesMatch(msgfplusInfo.QValue, qValueOrFDR, "QValue (or EFDR)", resultCount + 1);
                AssertValuesMatch(msgfplusInfo.PepQValue, pepQValueOrFDR, "PepQValue (or PepFDR)", resultCount + 1);

                resultIDsExamined++;
            }

            Console.WriteLine();
            Console.WriteLine("Read {0:N0} results", resultCount);
            Console.WriteLine("Examined scores for {0:N0} results", resultIDsExamined);

            foreach (var resultID in infoByResultID.Keys)
            {
                if (!matchedResultIDs.Contains(resultID))
                {
                    Assert.Fail(
                        "The reader did not return a PSM with ResultID {0}; " +
                        "likely this peptide has multiple proteins and the given ResultID " +
                        "is not for the first protein for this peptide in this scan", resultID);
                }
            }
        }

        [Test]
        [TestCase(@"MXQ202103181341_Auto1878805\QC_Shew_16_01-15f_08_4Nov16_Tiger_16-02-14_PrecursorInfo.txt",
            7718,
            "307.86; 522.80; 528.60; 388.71; 684.21; 323.19; 313.31; 343.93; 379.29; 339.96; 307.86")]
        public void TestPrecursorInfoFileReader(string inputFilePath, int expectedResultCount, string initialPrecursorMzValues)
        {
            var inputFile = FindFile(inputFilePath);

            var reader = new PHRPReader.Reader.PrecursorInfoFileReader();

            var precursorInfoData = reader.ReadPrecursorInfoFile(inputFile.FullName);

            Console.WriteLine("Loaded Precursor m/z info for {0} MSn scans", precursorInfoData.Count);

            if (expectedResultCount > 0)
            {
                Assert.AreEqual(expectedResultCount, precursorInfoData.Count);
            }

            var initialPrecursors = new List<double>();

            foreach (var item in initialPrecursorMzValues.Split(';'))
            {
                initialPrecursors.Add(double.Parse(item));
            }

            var scansShown = 0;

            foreach (var item in precursorInfoData)
            {
                Console.WriteLine("Scan {0,-2} has {1:F2} m/z precursor; {2}", item.Key, item.Value.PrecursorMz, item.Value.ScanFilterText);

                if (scansShown < initialPrecursors.Count)
                {
                    Assert.AreEqual(initialPrecursors[scansShown], item.Value.PrecursorMz);
                }

                scansShown++;
                if (scansShown > 10)
                    break;
            }
        }

        [Test]
        [TestCase(@"SIC202012190523_Auto1857747\Muscle_Mock_TMT16_PremixQC_Bane_16Dec20_20-09-10_ScanStats.txt", 38790)]
        public void TestScanStatsFileReader(string inputFilePath, int expectedResultCount)
        {
            var inputFile = FindFile(inputFilePath);

            var reader = new PHRPReader.Reader.ScanStatsReader();

            var scanStatsData = reader.ReadScanStatsData(inputFile.FullName);

            Console.WriteLine("Loaded scan stats data for {0} scans", scanStatsData.Count);

            if (expectedResultCount > 0)
            {
                Assert.AreEqual(expectedResultCount, scanStatsData.Count);
            }

            var scanCount = 0;

            foreach (var item in scanStatsData)
            {
                scanCount++;
                if (scanCount < 10 || scanCount % 100 == 0)
                {
                    Console.WriteLine("Scan {0,-4} at {1:F2} minutes: ScanType = {2}", item.Key, item.Value.ScanTimeMinutes, item.Value.ScanType);
                }
            }
        }

        [Test]
        [TestCase(@"SIC202012190523_Auto1857747\Muscle_Mock_TMT16_PremixQC_Bane_16Dec20_20-09-10_ScanStatsEx.txt", 38790)]
        public void TestExtendedScanStatsFileReader(string inputFilePath, int expectedResultCount)
        {
            var inputFile = FindFile(inputFilePath);

            var reader = new PHRPReader.Reader.ExtendedScanStatsReader();

            var extendedScanStatsData = reader.ReadExtendedScanStatsData(inputFile.FullName);

            Console.WriteLine("Loaded extended scan stats data for {0} scans", extendedScanStatsData.Count);

            if (expectedResultCount > 0)
            {
                Assert.AreEqual(expectedResultCount, extendedScanStatsData.Count);
            }

            var scanCount = 0;

            foreach (var item in extendedScanStatsData)
            {
                scanCount++;
                if (scanCount >= 10 && scanCount % 100 != 0)
                {
                    continue;
                }

                var scanDescription =
                    string.IsNullOrWhiteSpace(item.Value.CollisionMode) ?
                        "MS1 scan" :
                        "Collision mode=" + item.Value.CollisionMode;

                Console.WriteLine("Scan {0,-4}, {1,-19}: {2}", item.Key, scanDescription, item.Value.ScanFilterText);
            }
        }

        [Test]
        [TestCase(@"SIC202012190523_Auto1857747\Muscle_Mock_TMT16_PremixQC_Bane_16Dec20_20-09-10_SICStats.txt", 29620)]
        public void TestSICStatsFileReader(string inputFilePath, int expectedResultCount)
        {
            var inputFile = FindFile(inputFilePath);

            var reader = new PHRPReader.Reader.SICStatsReader();

            var sicStatsData = reader.ReadSICStatsData(inputFile.FullName);

            Console.WriteLine("Loaded SICStats data for {0} scans", sicStatsData.Count);

            if (expectedResultCount > 0)
            {
                Assert.AreEqual(expectedResultCount, sicStatsData.Count);
            }

            var scanCount = 0;

            foreach (var item in sicStatsData)
            {
                scanCount++;
                if (scanCount < 10 || scanCount % 100 == 0)
                {
                    Console.WriteLine("Fragmentation scan {0,-4}, precursor {1:F2} m/z, peak area {2}", item.Key, item.Value.MZ, item.Value.PeakArea);
                }
            }
        }

        [Test]
        [TestCase(PeptideHitResultTypes.MSAlign, @"MSA201807121159_Auto1609998\MSAlign_15ppm_0pt01_FDR_2012-01-03.txt")]
        [TestCase(PeptideHitResultTypes.MSGFPlus, @"MSG201802011337_Auto1547784\MSGFPlus_Tryp_DynSTYPhos_Stat_CysAlk_20ppmParTol.txt")]
        [TestCase(PeptideHitResultTypes.MSGFPlus, @"MSG202011150906_Auto1850387\MSGFPlus_Tryp_MetOx_StatCysAlk_20ppmParTol.txt")]
        [TestCase(PeptideHitResultTypes.MSPathFinder, @"MSP201804280928_Auto1578533\MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage.txt")]
        [TestCase(PeptideHitResultTypes.MaxQuant, @"MXQ202103181341_Auto1878805\MaxQuant_Tryp_Dyn_MetOx_NTermAcet_20ppmParTol.xml")]
        [TestCase(PeptideHitResultTypes.Sequest, @"Seq201212131618_Auto901984\sequest_DNA_N14_NE_Dyn_Met_Ox_Stat_Cys_Iodo.params")]
        [TestCase(PeptideHitResultTypes.TopPIC, @"TPC201808171745_Auto1624200\TopPIC_15ppmParTol_NumShift1_2018-08-16.txt")]
        [TestCase(PeptideHitResultTypes.XTandem, @"XTM201702142141_Auto1408699\xtandem_ETD_Rnd1Tryp_StatCysAlk_STYPhos_NoRefinement_20ppmParent_0pt5DaFrag.xml")]
        public void TestLoadSearchEngineParameters(PeptideHitResultTypes resultType, string parameterFilePath)
        {
            var parameterFile = FindFile(parameterFilePath);

            if (parameterFile.Directory == null)
                throw new NullReferenceException("Unable to determine the parent directory of the parameter file");

            var startupOptions = new StartupOptions
            {
                DisableOpeningInputFiles = true
            };

            var placeholderInputFilePath = Path.Combine(parameterFile.Directory.FullName, ReaderFactory.NON_EXISTENT_FILE_PLACEHOLDER_NAME);
            var factory = new ReaderFactory(placeholderInputFilePath, resultType, startupOptions);

            factory.SynFileReader.LoadSearchEngineParameters(parameterFile.Name, out var searchEngineParams);

            Console.WriteLine();
            Console.WriteLine("{0,-21} {1}", "Search Engine:", searchEngineParams.SearchEngineName);
            Console.WriteLine("{0,-21} {1} {2}", "Precursor tolerance:", searchEngineParams.PrecursorMassTolerancePpm, "ppm");
            Console.WriteLine("{0,-21} {1} {2}", "Precursor tolerance:", searchEngineParams.PrecursorMassToleranceDa, "Da");
            Console.WriteLine("{0,-21} {1}", "Min number termini:", searchEngineParams.MinNumberTermini);
            Console.WriteLine("{0,-21} {1}", "Enzyme:", searchEngineParams.Enzyme);

            if (!string.IsNullOrWhiteSpace(searchEngineParams.FastaFilePath))
            {
                Console.WriteLine("{0,-21} {1}", "FASTA File:", searchEngineParams.FastaFilePath);
            }

            foreach (var modInfo in searchEngineParams.ModList)
            {
                Console.WriteLine(modInfo.ToString());
            }
        }

        [Test]
        [TestCase(@"SIC202012190523_Auto1857747\Muscle_Mock_TMT16_PremixQC_Bane_16Dec20_20-09-10_ReporterIons.txt", 29620)]
        public void TestReporterIonsFileReader(string inputFilePath, int expectedResultCount)
        {
            var inputFile = FindFile(inputFilePath);

            var reader = new PHRPReader.Reader.ReporterIonsFileReader();

            var reporterIonData = reader.ReadReporterIonData(inputFile.FullName);

            Console.WriteLine("Loaded reporter ion data for {0} scans", reporterIonData.Count);

            if (expectedResultCount > 0)
            {
                Assert.AreEqual(expectedResultCount, reporterIonData.Count);
            }

            var scanCount = 0;

            foreach (var item in reporterIonData)
            {
                scanCount++;
                if (scanCount >= 10 && scanCount % 100 != 0)
                    continue;

                Console.WriteLine("Scan {0,-4}, collision mode {1}, max intensity {2:F2}", item.Key, item.Value.CollisionMode, item.Value.ReporterIonIntensityMax);

                var observedIonCount = 0;

                foreach (var reporterIon in item.Value.ReporterIons)
                {
                    if (reporterIon.Intensity > 0)
                        observedIonCount++;
                }

                if (observedIonCount == 0)
                {
                    continue;
                }

                var dataToShow = new List<string>
                {
                    string.Format("{0,-14}", "Ion"),
                    string.Format("{0,-14}", "Intensity"),
                    string.Format("{0,-14}", "Original Int."),
                    string.Format("{0,-14}", "S/N"),
                    string.Format("{0,-14}", "Resolution")
                };

                foreach (var reporterIon in item.Value.ReporterIons)
                {
                    var originalIntensity = reporterIon.OriginalIntensity.HasValue
                        ? reporterIon.OriginalIntensity.Value.ToString("0")
                        : string.Empty;

                    var signalToNoise = reporterIon.SignalToNoise.HasValue
                        ? reporterIon.SignalToNoise.Value.ToString(CultureInfo.InvariantCulture)
                        : string.Empty;

                    var resolution = reporterIon.Resolution.HasValue
                        ? reporterIon.Resolution.Value.ToString(CultureInfo.InvariantCulture)
                        : string.Empty;

                    dataToShow[0] += string.Format("{0,-7:F2}", reporterIon.MZ);
                    dataToShow[1] += string.Format("{0,-7:F0}", reporterIon.Intensity);
                    dataToShow[2] += string.Format("{0,-7}", originalIntensity);
                    dataToShow[3] += string.Format("{0,-7}", signalToNoise);
                    dataToShow[4] += string.Format("{0,-7}", resolution);
                }

                Console.WriteLine(dataToShow[0]);
                Console.WriteLine(dataToShow[1]);
                Console.WriteLine(dataToShow[2]);
                Console.WriteLine(dataToShow[3]);
                Console.WriteLine(dataToShow[4]);
                Console.WriteLine();
            }
        }

        [Test]
        [TestCase(
            @"TPC201808171745_Auto1624200\MZ20170525_Mnx_PFA_toppic_syn.txt",
            0,
            "M.STDYSKMTDVNEIHDSAILEHFRNGIGHKTLVISPSYPYMFVGIIKELIGDTVMIDVETTHFAQLENREWYIHIHNIEVFYIERPGAPKIPKLEDY.-|L.SAASNVASVNDPLFDFFNKHMGKQILIITESSQLNILGQTFRPIFCGKVAEVEPGHLTLSPVTIKILNAPFHKFPIPLSIPFEKIAHFTTDVD(CSM)[-1.98828]RIPLV.-|Y.PYMFVGIIKELIGDTVMIDVETTHFAQLENREWYIHIHNIEVFYIERPGAPKIPKLEDY.-|A.SNVASVNDPLFDFFNKHMGKQILIITESSQLNILGQ(TFRPIFCGKVAEV)[-2.00464]EPGHLTLSPVTIKILNAPFHKFPIPLSIPFEKIAHFTTDVDCSMRIPLV.-|M.STDYSKMTDVNEIHDSAILEHFRNGIGHKTLVISPSYPYMFVGIIKELIGDTVMIDVETTHFAQLENRE(WYI)[30.01443]HIHNIEVFYIERPGAPKIPKLEDY.-",
            "1.18E-50|7.86E-33|1.98E-28|4.03E-27|2.34E-21|1.85E-19",
            "41.28164|26.86965|8.35164|41.15642|26.71951")]
        [TestCase(
            @"TPC202102011550_Auto1869702\Hubmap_nanoPOTs_top_down_QC_intact_20ng_FAIMS_LowNCE_r3_toppic_syn.txt",
            0,
            "M.AKGQSLQDPFLNALRRERVPVSIYLVNGIKLQGQVESFDQFVILLKNTVSQMVYKHAISTVVPARPFNVAGHQNAQ.G|M.SDKIIYLSDDSFENDVLKADLPVLVDFWA(E)[-1.9870]WCGPCKMIAPILDDVAEEYAGRVTIAKLNVDQNNVSPAKYGVRGIPTLLLFKNGELAATKVGALSKTQLKEFIDAQI.-|-.MQIHLSGHHIEITESLRAYVEEKFSKLERHFEQINNVHVVLNVEKMQQIAEARINLTGGEVFATSEHADMYAAIDVLIDKLDRQVIKHKEKLTKH.-|M.PQTKVIEQLKENLQIAYRQAIDADAKLDELKKAGHGKFTHIFTAEQGFVVSSNRFLPYVQELVDDLTKLQQHTRLDPVALETLVRQLATLLQTLHAFKKQS.-|M.SDKIIYLSDDSFENDVLKADLPVLVDFWAEW(CGPCKM)[14.0041]IAPILDDVAEEYAGRVTIAKLNVDQNNVSPAKYGVRGIPTLLLFKNGELAATKVGALSKTQLKEFIDAQI.-",
            "1.51E-26|8.25E-24|2.59E-23|8.57E-22|1.19E-21|3.12E-20",
            "32.8273|22.10651|38.98326|20.06614|27.59561|31.43993")]
        public void TestTopPICReader(
            string inputFilePath,
            int expectedResultCount,
            string initialPeptides,
            string initialEValues,
            string initialDelMPPM)
        {
            var inputFile = FindFile(inputFilePath);

            var options = new StartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var expectedPeptides = initialPeptides.Split('|').ToList();

            var expectedEValues = SplitDelimitedDoubles(initialEValues);
            var expectedDelMPPM = SplitDelimitedDoubles(initialDelMPPM);

            var reader = new ReaderFactory(inputFile.FullName, options);
            var resultCount = 0;
            var peptideMatchCount = 0;
            var specEValueMatchCount = 0;
            var delMMatchCount = 0;

            while (reader.MoveNext())
            {
                var peptide = reader.CurrentPSM.Peptide;
                var specEvalue = reader.CurrentPSM.MSGFSpecEValue;
                var delMPPM = reader.CurrentPSM.MassErrorPPM;

                if (resultCount < expectedPeptides.Count)
                {
                    Assert.AreEqual(expectedPeptides[resultCount], peptide);
                    peptideMatchCount++;
                }

                if (resultCount < expectedEValues.Count)
                {
                    AssertValuesMatch(expectedEValues[resultCount], specEvalue, "SpecEValue", resultCount + 2);
                    specEValueMatchCount++;
                }

                if (resultCount < expectedDelMPPM.Count)
                {
                    AssertValuesMatch(expectedDelMPPM[resultCount], delMPPM, "DelM PPM", resultCount + 2);
                    delMMatchCount++;
                }

                resultCount++;
            }

            Console.WriteLine();
            Console.WriteLine("Read {0:N0} results", resultCount);
            Console.WriteLine("{0} peptides matched expected values", peptideMatchCount);
            Console.WriteLine("{0} SpecEValues matched expected values", specEValueMatchCount);
            Console.WriteLine("{0} DelMPPM values expected values", delMMatchCount);

            if (expectedResultCount > 0)
                Assert.AreEqual(expectedResultCount, resultCount);
        }

        [Test]
        [TestCase(
            @"TPC201808171745_Auto1624200\MZ20170525_Mnx_PFA_toppic_syn.txt",
            "1|5|10|13",
            "M.STDYSKMTDVNEIHDSAILEHFRNGIGHKTLVISPSYPYMFVGIIKELIGDTVMIDVETTHFAQLENREWYIHIHNIEVFYIERPGAPKIPKLEDY.-|M.STDYSKMTDVNEIHDSAILEHFRNGIGHKTLVISPSYPYMFVGIIKELIGDTVMIDVETTHFAQLENRE(WYI)[30.01443]HIHNIEVFYIERPGAPKIPKLEDY.-|T.DVNEIHDSAILEHFRNGIGHKTLVISPSYPYMFVGIIKELIGDTVMIDVETTHFAQLENREWYIHIHNIEVFYIERPGAPKIPKLEDY.-|E.SRMSGECAPNVSVSVSTSHTTISGGGSRGGGGGGYGSGGSSYGSGGGSYGSGGGGGGGRGSYGSGGSSYGSGGGSYGSGGGGGGHGSYGSGSSSGGY(RGGSGGGGGGSSGGRGSGGGSSGGSIGGRGSSSGG)[-44.84910]V.K",
            "gi|145554186|gb|ABP68889.1|;gi|145554186|gb|ABP68889.1|;gi|145554186|gb|ABP68889.1|;Contaminant_K2C1_HUMAN",
            "2665292844,_,1.18E-50,0,0|1095085274,_,2.34E-21,0,0|49222185.4,_,1.71E-14,0,0|38609340.27,_,1.59E-04,0.076923,0.076923")]
        [TestCase(
            @"TPC202102011550_Auto1869702\Hubmap_nanoPOTs_top_down_QC_intact_20ng_FAIMS_LowNCE_r3_toppic_syn.txt",
            "10|50|250|400|490",
            "M.SLELLSKLETKIQATLETIELLKMELEEEKQKTSTLSEQNNQLIEQNQQLQQELTSWNEKVTGLVGLLNSEI.-|-.MNKSELIEKIASGADISKAAAGRALDSFIAAVTEGLKEGDKISLVGFGTFEVRERAERTGRNPQTGEEIKIAAAKIPAFKAGKALKDAVN.-|M.ADTTVEKLATEVGKSVERLIEQFSQAGIKKGQADNVSEAEKQQLLDYLKKQHGGDNAPT.K|-.MINPKKIEEMAKQLSDSLPSGLKQFAGEFEERSKQVLQNQLLKLDMVSREEFEVQQHVLLKTREKLEALQAQVNELEKKLNATD.-|M.TISVLICDDSAMARKQMARTLPKEWDVEITYATNGAEGLDAIRAGKGEVVFLDLNMPVMDGYEVLQTVQQNDLPALIIVV.S",
            "SO_0335;SO_1797;SO_1204;SO_4329;SO_0570",
            "1498000,38.2275287,4.86E-19,0,0|4955000000,44.8307782,8.10E-13,0,0|1870000,46.4272896,4.40E-06,0,0|21980000,54.930232,9.64E-04,0,0|1090000,48.8775848,6.98E-02,0.006122449,0.01")]
        public void TestTopPICScores(
            string inputFilePath,
            string resultIdList,
            string peptidesByResultID,
            string proteinsByResultID,
            string scoresByResultId      // Feature_Intensity, Feature_Score, EValue, QValue, Proteoform_QValue; use an underscore for "value not defined"
        )
        {
            var inputFile = FindFile(inputFilePath);

            var options = new StartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var resultIDs = resultIdList.Split('|').ToList();

            var scores = scoresByResultId.Split('|').ToList();
            var peptides = peptidesByResultID.Split('|').ToList();
            var proteins = proteinsByResultID.Split(';').ToList();

            var infoByResultID = new Dictionary<int, TopPICInfo>();
            var matchedResultIDs = new SortedSet<int>();

            for (var i = 0; i < resultIDs.Count; i++)
            {
                var scoreValues = SplitDelimitedDoubles(scores[i], ',');

                var topPICInfo = new TopPICInfo
                {
                    Peptide = peptides[i],
                    Protein = proteins[i],
                    FeatureIntensity = scoreValues[0],
                    FeatureScore = scoreValues[1],
                    EValue = scoreValues[2],
                    QValue = scoreValues[3],
                    ProteoformQValue = scoreValues[4]
                };

                var resultId = int.Parse(resultIDs[i]);
                infoByResultID.Add(resultId, topPICInfo);
            }

            var reader = new ReaderFactory(inputFile.FullName, options);
            var resultCount = 0;
            var resultIDsExamined = 0;

            while (reader.MoveNext())
            {
                resultCount++;

                if (!infoByResultID.TryGetValue(reader.CurrentPSM.ResultID, out var topPICInfo))
                    continue;

                matchedResultIDs.Add(reader.CurrentPSM.ResultID);

                Assert.AreEqual(topPICInfo.Peptide, reader.CurrentPSM.Peptide);
                Assert.AreEqual(topPICInfo.Protein, reader.CurrentPSM.ProteinFirst);

                var featureIntensity = GetScore(reader.CurrentPSM, "Feature_Intensity");
                var featureScore = GetScore(reader.CurrentPSM, "Feature_Score", string.Empty, string.Empty);

                var qValue = GetScore(reader.CurrentPSM, "QValue");
                var proteoformQValue = GetScore(reader.CurrentPSM, "Proteoform_QValue", "Proteoform_FDR");

                AssertValuesMatch(topPICInfo.FeatureIntensity, featureIntensity, "Feature_Intensity", resultCount + 1);
                AssertValuesMatch(topPICInfo.FeatureScore, featureScore, "Feature_Score", resultCount + 1);

                AssertValuesMatch(topPICInfo.EValue, reader.CurrentPSM.MSGFSpecEValue, "EValue", resultCount + 1);
                AssertValuesMatch(topPICInfo.QValue, qValue, "QValue", resultCount + 1);
                AssertValuesMatch(topPICInfo.ProteoformQValue, proteoformQValue, "Proteoform_QValue", resultCount + 1);

                resultIDsExamined++;
            }

            Console.WriteLine();
            Console.WriteLine("Read {0:N0} results", resultCount);
            Console.WriteLine("Examined scores for {0:N0} results", resultIDsExamined);

            foreach (var resultID in infoByResultID.Keys)
            {
                if (!matchedResultIDs.Contains(resultID))
                {
                    Assert.Fail(
                        "The reader did not return a PSM with ResultID {0}; " +
                        "likely this peptide has multiple proteins and the given ResultID " +
                        "is not for the first protein for this peptide in this scan", resultID);
                }
            }
        }

        private void AssertValuesMatch(double? expectedValue, string observedValueAsText, string valueDescription, int lineNumber, double delta = 0.0001)
        {
            if (!expectedValue.HasValue)
            {
                // The expected value is null, which means the observed value should be not defined
                Assert.AreEqual(string.Empty, observedValueAsText,
                    "{0} value is defined as {1} but was expected to be an empty string (not defined)",
                    valueDescription, observedValueAsText);
                return;
            }

            if (double.TryParse(observedValueAsText, out var observedValue))
            {
                Assert.AreEqual(expectedValue.Value, observedValue, delta,
                    "{0} values do not match; expected {1} but actually {2}",
                    valueDescription, expectedValue, observedValue);
            }
            else
            {
                Console.WriteLine("{0} value is not numeric on line {1}: {2}", valueDescription, lineNumber, observedValueAsText);
                Assert.Fail();
            }
        }

        private FileInfo FindFile(string relativeFilePath)
        {
            var localDirPath = Path.Combine("..", "..", "..", "Data", "TestData");
            const string remoteDirPath = @"\\proto-2\UnitTest_Files\PeptideHitResultsProcessor";

            var localFile = new FileInfo(Path.Combine(localDirPath, relativeFilePath));

            if (localFile.Exists)
            {
                return localFile;
            }

            // Look for the file on Proto-2
            var remoteFile = new FileInfo(Path.Combine(remoteDirPath, relativeFilePath));

            if (remoteFile.Exists)
            {
                return remoteFile;
            }

            var msg = string.Format("File not found: {0}; checked in both {1} and {2}", relativeFilePath, localDirPath, remoteDirPath);

            Console.WriteLine(msg);
            Assert.Fail(msg);

            return null;
        }

        private string GetScore(PSM currentPSM, string scoreColumnName, string alternativeName = "", string valueIfNotFound = "-1")
        {
            if (currentPSM.TryGetScore(scoreColumnName, out var value))
            {
                return value;
            }

            if (!string.IsNullOrEmpty(alternativeName) && currentPSM.TryGetScore(alternativeName, out var alternativeValue))
            {
                return alternativeValue;
            }

            return valueIfNotFound;
        }

        private double GetValueOrDefault(double? nullableValue, double valueIfNull = 0)
        {
            return nullableValue ?? valueIfNull;
        }

        private static List<double?> SplitDelimitedDoubles(string delimitedValues, char delimiter = '|')
        {
            var splitList = delimitedValues.Split(delimiter).ToList();
            var parsedValues = new List<double?>();

            if (string.IsNullOrWhiteSpace(delimitedValues))
                return parsedValues;

            foreach (var item in splitList)
            {
                if (item.Equals("_"))
                {
                    parsedValues.Add(null);
                }
                else
                {
                    parsedValues.Add(double.Parse(item));
                }
            }

            return parsedValues;
        }
    }
}