This page contains detailed descriptions of the columns in synopsis files (_syn.txt)
created from [X!Tandem](https://www.thegpm.org/tandem/) results.

For other search tools, see:

* The [MS-GF+ Synopsis File Columns](SynopsisFileColumns_MSGFPlus.md) page for a description the columns in MS-GF+ \_msgfplus\_syn.txt and \_msgfplus\_fht.txt files
* The [SEQUEST Synopsis File Columns](SynopsisFileColumns_SEQUEST.md) page for a description the columns in SEQUEST \_syn.txt and \_fht.txt files


## _xt.txt File Columns


| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|
| Result\_ID | Group\_ID | Scan | Charge | Peptide\_MH | Peptide\_Hyperscore | Peptide\_Expectation\_Value\_Log(e) | Multiple\_Protein\_Count | Peptide\_Sequence | DeltaCn2 | y\_score | y\_ions | b\_score | b\_ions | Delta\_Mass | Peptide\_Intensity\_Log(I) | DelM\_PPM |


## _xt.txt File Column Details


1. Result_ID 
  * This value is the row number of the peptide entry

2. Group_ID 
  * X!Tandem-assigned unique ID value

3. Scan 
  * This is the scan number of the MS/MS spectrum that resulted in the peptide identification of interest. The scan number is a closely related to retention time.

4. Charge 
  * This is the charge state of the parent ion. This value is often used in filtering criteria since larger charge state peptides can produce a larger number of potential fragment ions. Thus, many criteria often have different cutoffs for different charge states.

5. Peptide_MH 
  * MH or (M + H)+ includes the mass of the peptide, calculated from theoretical amino acid masses, for the peptide match to the current spectrum of interest. Largely because of convention, this value is reported as the mass of the peptide + mass of proton. The default mass calculation for XTandem uses monoisotopic amino acid mass values

6. Peptide_Hyperscore 
  * This is X!Tandem's main scoring value related to peptide confidence. X!Tandem assumes that the peptide with the highest hyperscore is correct, and all others are incorrect. There is not yet a generally accepted cutoff for hyperscore or DelCn2 since we don’t yet have enough experience with X!Tandem. The following values show a very rough equivalence between XCorr and hyperscore:

7. Peptide_Expectation_Value_Log(e) (aka Log_EValue) 
  * This is the base-10 Log of the Expectation value. The expectation value is a confidence measure reported by X!Tandem; the more negative the Log of the EValue, the more confident we are in the result. This value is a population-based statistic, in that it depends on the other matches in a given dataset and on the background false positive matches. Ron Beavis will commonly specify that this value should be &lt;= -2. If you’re going to filter on Log_EValue, you should probably filter on -1.3 or lower (corresponds to E_Value = 0.05). 
  * The E_Value is similar to a P_Value as reported by Blast, in that an E_Value of 0.05 (Log_EValue = -1.3) means there is a 5% probability that the result is incorrect. An E_Value of 0.001 (Log_EValue = -3) means there is a 0.1% probability that the result is incorrect. 

8. Multiple_Protein_Count 
  * Count of the number of proteins (ORFs) that a peptide sequence is found in. Users mostly use this field in filtering out "degenerate" or "multiorf" peptide hits in tabulated lists of peptide IDs.

9. Peptide_Sequence 
  * This is the peptide sequence that matches the spectrum.

10. DeltaCn2 
  * shows the degree by which lower-ranked peptide scores differ from the Hyperscore score of the best match. For example, if a spectrum has a top match with a score of 50 and a second highest score of 40, then the DeltaCn2 value is (50-40)/ 50 = 0.2. Matches with DeltCn2 scores of 0.04 indicate very close scores. Firm guidance on filtering based on DeltCn2 scores isn’t available at this point in time, but it is suggested that if you filter on DeltaCn2 then use &gt;=0.05 rather than &gt;=0.1.

11. y_score 
  * Score computed from the matches to the expected y-ions

12. y_ions 
  * The number of y-ions that were matched

13. b_score 
  * Score computed from the matches to the expected b-ions

14. b_ions 
  * The number of b-ions that were matched

15. Delta_Mass 
  * This is the mass error between the identified peptide and the parent ion mass, i.e. 'Theoretical Parent Mass' - 'Parent Mass Observed'

16. Peptide_Intensity_Log(I) 
  * The log of the total intensity of all of the ions in the spectrum

17. DelM_PPM 
  * Mass Difference (in ppm) between the observed parent ion and the computed mass of the identified peptide. This value is automatically corrected if the second or third isotope is chosen for fragmentation. In other words, if the non-corrected DelM value is 1.05 Da, then the DelM_PPM value for a 1500 Da peptide is reported as 33 ppm and not 700 ppm.
`
DelM_Corrected = 1.05 - 1 = 0.05 Da
DelM_PPM = 0.05 / (1500 / 1E6) = 33 ppm
```

Note that X!Tandem only reports the top match for each spectrum but it does report the next highest score, which is why the DeltaCn2 is also calculated.
