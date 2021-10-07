This page contains detailed descriptions of the columns in synopsis files (_syn.txt)
created from [MS-GF+](https://github.com/MSGFPlus/msgfplus) results.

For other search tools, see:

* The [SEQUEST Synopsis File Columns](SynopsisFileColumns_SEQUEST.md) page for a description the columns in SEQUEST \_syn.txt and \_fht.txt files
* The [XTandem Synopsis File Columns](SynopsisFileColumns_XTandem.md) page for a description the columns in X!Tandem \_xt.txt files


## _msgfplus_syn.txt File Columns


| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|----|----|----|
| ResultID | Scan | FragMethod | SpecIndex | Charge | PrecursorMZ | DelM | DelM\_PPM | MH | Peptide | Protein | NTT | DeNovoScore | MSGFScore | MSGFDB\_SpecEValue | Rank\_MSGFDB\_SpecEValue | EValue | QValue | PepQValue | IsotopeError |


## _msgfplus_syn.txt File Column Details

1. Result_ID
  * This value is the row number of the peptide entry

2. Scan
  * This is the scan number of the MS/MS spectrum that resulted in the peptide identification of interest. The scan number is a closely related to retention time.

3. FragMethod
  * Fragmentation method for the given MS/MS spectrum. Will be CID, ETD, or HCD. However, when spectra from the same precursor are merged, fragmentation methods of merged spectra will be shown in the form "FragMethod1/FragMethod2/..." (e.g. CID/ETD, CID/HCD/ETD).

4. SpecIndex
  * Spectrum Index assigned by MSGF+. The first spectrum has index 1, the second has index 2, and so on.

5. Charge
  * This is the charge state of the parent ion.

6. PrecursorMZ
  * m/z value of the precursor ion

7. DelM
  * This is the mass error between the identified peptide and the parent ion mass, i.e. 'Theoretical Parent Mass' - 'Parent Mass Observed'

8. DelM_PPM
  * Mass Difference (in ppm) between the observed parent ion and the computed mass of the identified peptide. This value is automatically corrected if the second or third isotope is chosen for fragmentation. In other words, if the non-corrected DelM value is 1.05 Da, then the DelM_PPM value for a 1500 Da peptide is reported as 33 ppm and not 700 ppm.
```
DelM_Corrected = 1.05 - 1 = 0.05 Da
DelM_PPM = 0.05 / (1500 / 1E6) = 33 ppm
```

9. MH
  * MH or (M + H)+ includes the mass of the peptide, calculated from theoretical amino acid masses, for the peptide match to the current spectrum of interest. Largely because of convention, this value is reported as the mass of the peptide plus the mass of a proton (1.007276467). MSGF+ uses monoisotopic amino acid mass values.

10. Peptide
  * The identified peptide, with prefix and suffix.

11. Protein
  * The name of the protein this peptide comes from.

12. NTT
  * Number of tryptic terminii (or Unused, if no protease was specified). Note that the N- and C-terminus of a protein are both considered to be valid termini.
  * 0 = Non tryptic 
  * 1 = Partially tryptic 
  * 2 = Fully tryptic 
  * For more information, including examples of various tryptic states, see the [Cleavage State](CleavageState.md) page.

13. DeNovoScore
  * The MSGFScore of the optimal scoring peptide.  Larger scores are better.

14. MSGFScore
  * This is MSGF+'s main scoring value for the identified peptide. Larger scores are better.

15. MSGFDB_SpecEValue (for MSGF+); was previously MSGFDB_SpecProb
  * This is MSGF+'s main scoring value related to peptide confidence (spectrum level e-value) of the peptide-spectrum match. MSGF+ assumes that the peptide with the lowest MSGFDB_SpecEValue value (closest to 0) is correct, and all others are incorrect. This value is similar to the MSGF_SpecProb value computed by the MSGF tool, but the range of values computed by MSGF+ is different than the range computed by the MSGF program from November 2010.
* SpecEValue represents the spectrum specific probability of getting a score better than the proposed candidate. This is also known as the generating function, and it's a rather complex math/statistical method, described in depth in the [https://www.ncbi.nlm.nih.gov/pubmed/25358478 MS-GF+ manuscript] in Nature Communications.

16. Rank_MSGFDB_SpecEValue (for MSGF+); was previously Rank_MSGFDB_SpecProb
  * Computed on a per-spectrum basis using all of the candidate peptides identified for that spectrum. The peptide with the lowest MSGFDB_SpecEValue has rank 1

17. EValue (for MSGF+, was previously PValue)
  * Probability that a match with this MSGFDB_SpecEValue is spurious; the lower this number (closer to 0), the better the match. This is a database level e-value, representing the probability that a random PSM has an equal or better score against a random database of the same size.

18. QValue (for MSGF+, was previously FDR)
  * If MSGF+ searches a target/decoy database, then the QValue (FDR) is computed based on the distribution of MSGFDB_SpecEValue values for forward and reverse hits. If the target/decoy search was not used, then this column will be EFDR and is an estimated FDR.

  * QValue is defined as the minimum false discovery rate (FDR) at which the test may be called significant
    * When computing QValue, data is sorted by SpecEValue.
  * QValue is a spectrum-level FDR and is computed using the formula ReversePeptideCount ÷ ForwardPeptideCount
  * Thus, if you filter on QValue &lt; 0.01 you are applying a 1% FDR filter

19. PepQValue (for MSGF+, was previously PepFDR)
  * Peptide-level QValue (FDR) estimated using the target-decoy approach; only shown if a target/decoy search was used. If multiple spectra are matched to the same peptide, only the best-scoring match is retained and used to compute FDR.

  * PepQValue is a peptide-level FDR threshold, and will always be lower than QValue
  * The same PepQValue value is given to all PSMs of the same peptide. Thus, even a low-quality PSM may get a good PepQValue (lower number) if the PSM has a high-quality "sibling" PSM.
    * Do not filter on PepQValue when counting the number of identified PSMs.
  * If comparing data to Sequest results filtered using a Target/Decoy approach, use QValue for filtering, not PepQValue

