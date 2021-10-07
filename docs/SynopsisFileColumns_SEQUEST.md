This page contains detailed descriptions of the columns in synopsis files (_syn.txt)
created from SEQUEST results.

For other search tools, see:

* The [MS-GF+ Synopsis File Columns](SynopsisFileColumns_MSGFPlus.md) page for a description the columns in MS-GF+ \_msgfplus\_syn.txt and \_msgfplus\_fht.txt files
* The [XTandem Synopsis File Columns](SynopsisFileColumns_XTandem.md) page for a description the columns in X!Tandem \_xt.txt files


## _syn.txt File Columns


| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|----|----|----|
| HitNum | ScanNum | ScanCount | ChargeState | MH | XCorr | DelCn | Sp | Reference | MultiProtein | Peptide | DelCn2 | RankSp | RankXc | DelM | XcRatio | Ions\_Observed | Ions\_Expected | NumTrypticEnds | DelM\_PPM |


## _syn.txt File Column Details


1. HitNum 
  * This value is the row number of the peptide entry, usually sorted according to descending XCorr value.

2. ScanNum 
  * This is the scan number of the MS/MS spectrum that resulted in the peptide identification of interest. The scan number is a closely related to retention time.

3. ScanCount 
  * Under some instances, the Thermofinnigan extract_msn (aka lcq_dta) program will group (combine via summing) multiple spectra into one spectrum file in order to reduce duplication and increase signal-to-noise. In normal proteomic samples processed by extract_msn, ~2% of the spectra are grouped. A value of 1 indicates only one scan was used. If more than one scan is grouped, then the Number of Scans column indicates the range between the first spectrum and the last spectrum grouped, thus it is not actually a scan count value. For example, if the program groups spectra 8163, 8165, 8168, 8171, 8174, and 8177, then Number of Scans = 8177-8163+1 = 15. Unfortunately, when using extract_msn to convert .Raw files to .Dta files, we do not know which exact scans are grouped together, just the range of scans as indicated by this column.

4. ChargeState 
  * This is the charge state of the parent ion. This value is often used in filtering criteria since larger charge state peptides can produce a larger number of potential fragment ions. Thus, many criteria often have different cutoffs for different charge states.

5. MH 
  * MH or (M + H)+ includes the mass of the peptide, calculated from theoretical amino acid masses, for the peptide match to the current spectrum of interest. Largely because of convention, this value is reported as the mass of the peptide + mass of proton. The mass calculation uses average amino acid mass values&nbsp;or monoisotopic masses, depending on the parameter file setting for parent mass value.

6. Xcorr (aka cross correlation score) 
  * This is SEQUEST’s main scoring value related to bestness-of-fit between the observed MS/MS spectrum and the theoretical MS/MS spectrum calculated from the candidate peptide sequence. This score is calculated for the top 200 or 500 preliminary scoring peptide sequence candidates (dependson param file setting).

7. DelCn 
  * SEQUEST calculates the deltaCn using the Xcorr of the top hit and the nth ranked hit using: (Xcorr(top hit) – Xcorr(n))/Xcorr(top hit). Thus, the deltaCn for the top hit is (Xcorr(top hit) – Xcorr(top hit))/Xcorr(top hit)=0. More specifically, it is simply a report of the DelCn value found for the given peptide in the out file. 
  * Note: This score is often of little value as calculated here. Instead, the DeltaCn2 value better reflects peptide confidence. This column is still retained in order to have compatibility with previous versions of the syn/fht files, and for filtering peptides for use in the AMT process. 

8. Sp 
  * Preliminary score. This is a score that Sequest uses to do an initial scoring of the all peptide candidates (these can easily number of 100,000 or more). It is quicker to do, but the score is less robust than Xcorr, when considering peptide confidence.

9. Reference 
  * ORF reference. This is the ORF name of the protein wherein the current peptide sequence was detected. For the fht and instances where the peptide is contained in many ORFs, only the first ORF is listed. All ORFs are outputted in the _syn.txt (up to 50 of them) for the current peptide.

10. MultiProtein 
  * Number of multiple ORFs. If a peptide sequence can be found in x ORFs, then this value is +(x-1) indicating the number of additional ORFs. Users mostly use this field in filtering out “degenerate” or “multiorf” peptide hits in tabulated lists of peptide IDs.

11. Peptide 
  * This is the peptide sequence that matches the spectrum.

12. DelCn2 
  * For the nth ranked hit, deltaCn2 is (Xcorr(n) – Xcorr(n+1))/Xcorr(n). Note that this value is calculated by the syn-fht summary generator, but this is the DeltaCn that should be used for filtering purposes (and what other people use in the scientific community). NOTE: “rank” here refers to rank by Xcorr. 

13. RankSp 
  * This value is used in the PNNL-confidence discriminant. The determination of this score is done by SEQUEST, where a list or array of peptides are sorted according to decreasing Sp score and a rank is assigned to each peptide sequence (e.g. the topmost entry would have a RankSp=1)

14. RankXc 
  * RankXc is the Xcorr rank. The explaination of this is similar to that given for RankSp except that the list of peptides are sorted by decreasing Xcorr. This value should be set to 1 for the fht file for all peptide hits, but is usually 1- 10 (or more) in the syn file. This is used in the PNNL-confidence discriminant.

15. DelM 
  * This is the mass error in the parent ion mass. This is used in the PNNL confidence discriminant. This is largely calculated from the SEQUEST output and is calculated from 'Theoretical Parent Mass' - 'Parent Mass Observed'. 
  * The user should be aware of some minor arithmetic that is occurring within the fht-syn summary generator. Because our data analysis uses average masses in the calculation of theoretical mass, this is because the low resolution ion trap (LCQ &amp; LTQ) derived masses most often agree with this value. However, for small singly charged species, measurements of the parent mass from an ion trap MS often match monoisotopic masses better. Therefore, for singly charged species, a correction factor is calculated and is used in the calculation of delM. On occasion large mass errors ( &gt;1 Da) can be seen for peptides with high Xcorr values, this is often due to the selection of a isotopic peak from a highly abundant species.

16. XcRatio 
  * Calculated using Xcorr(n)/Xcorr(top hit). This is used by the PNNL- confidence discriminant. The PNNL-confidence discriminant uses this value especially for judging the quality of runner-up or secondary hits.

17. Ions_Observed 
  * Number of b ions, y ions, etc. matched by Sequest

18. Ions_Expected 
  * Theoretical number of b ions, y ions, etc.

19. NumTrypticEnds (number of tryptic cleavage sites) 
  * This is the number of terminii that conform to the expected cleavage behavior of trypsin (i.e. C-terminal to R and K). Note that K-P and R-P do not qualify as tryptic cleavages because of the proline rule. However, the protein N-terminus and protein C-terminus do count as tryptic cleavage sites. Values can be 0, 1, or 2: 
  * 0 = Non tryptic 
  * 1 = Partially tryptic 
  * 2 = Fully tryptic 
  * For more information, including examples of various tryptic states, see the [Cleavage State](CleavageState.md) page.

20. DelM_PPM 
  * Mass Difference (in ppm) between the observed parent ion and the computed mass of the identified peptide. This value is automatically corrected if the second or third isotope is chosen for fragmentation. In other words, if the non-corrected DelM value is 1.05 Da, then the DelM_PPM value for a 1500 Da peptide is reported as 33 ppm and not 700 ppm.
```
DelM_Corrected = 1.05 - 1 = 0.05 Da
DelM_PPM = 0.05 / (1500 / 1E6) = 33 ppm
```

### Old Columns

17. PassFilt 
  * This score is not calculated by SEQUEST, but is calculated from syn-fht summary generator using Xcorr, delcn, RankXc and the number of tryptic termini. This score is need for the PNNL confidence discriminant, but most peptides passing the Yates criteria, also should pass this filter. Only one or zero appears in this column.

18. MScore 
  * This score measures whether the peptide under consideration has fragmentation pattern that is consistent with that observed for high confident peptide assignments. It is an empirical measure and not derived from SEQUEST scores. As rough guide, MScore &gt; 10 qualifies as a good score.
