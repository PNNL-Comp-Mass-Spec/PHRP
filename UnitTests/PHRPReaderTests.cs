using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using PHRPReader;
using PHRPReader.Data;

namespace PHRP_UnitTests
{
    [TestFixture]
    public class PHRPReaderTests
    {

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
            string initialDelMPPM)
        {
            var inputFile = FindFile(inputFilePath);

            var options = new PHRPStartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var expectedPeptides = initialPeptides.Split('|').ToList();

            var expectedSpecEValues = SplitDelimitedDoubles(initialSpecEValues);
            var expectedDelMPPM = SplitDelimitedDoubles(initialDelMPPM);

            var reader = new PHRPReader.PHRPReader(inputFile.FullName, options);
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

            var options = new PHRPStartupOptions
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

            var reader = new PHRPReader.PHRPReader(inputFile.FullName, options);
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

            var options = new PHRPStartupOptions
            {
                LoadMSGFResults = true,
                LoadModsAndSeqInfo = true,
                LoadScanStatsData = false
            };

            var expectedPeptides = initialPeptides.Split('|').ToList();

            var expectedEValues = SplitDelimitedDoubles(initialEValues);
            var expectedDelMPPM = SplitDelimitedDoubles(initialDelMPPM);

            var reader = new PHRPReader.PHRPReader(inputFile.FullName, options);
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

            var options = new PHRPStartupOptions
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

            var reader = new PHRPReader.PHRPReader(inputFile.FullName, options);
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
            if (nullableValue.HasValue)
                return nullableValue.Value;

            return valueIfNull;
        }

        private static List<double?> SplitDelimitedDoubles(string delimitedValues, char delimiter = '|')
        {
            var splitList = delimitedValues.Split(delimiter).ToList();
            var parsedValues = new List<double?>();

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