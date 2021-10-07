Cleavage State describes whether a peptide has the expected residues just before the peptide and at the end of the peptide, as defined by the enzyme (protease) used to cleave the peptide. The most common enzyme used in proteomics is Trypsin. Trypsin cleaves after Lysine (K) or Arginine (R) residues, though it should not cleave between a Lysine and Proline (KP) or an Arginine and Proline (RP), since Proline is a cyclic amino acid that bonds with the previous amino acid and thus prevents lysis by trypsin. 

Cleavage State is commonly represented by a number, where: 
* 2 indicates complete digestion, aka fully tryptic if the protease was Trypsin 
* 1 indicates partial digestion, aka partially tryptic 
* 0 indicates incomplete digestion, aka non-tryptic

## Enzyme List

The following table lists several common enzymes and their cleavage rules. The Prefix Residue(s) column lists the residues that should precede a peptide, since the enzyme should cleave after these residues (unless Cleavage Offset is 0, then the enzyme cleaves before the Prefix Residue). The Final Residue column lists the residues that should be at the end of a peptide, again because the enzyme should cleavage after them (unless Cleavage Offset is 0). The Exception Residue columns list the cases where cleavage should not occur if the given residue is present after a Prefix or Final Residue. 

| Enzyme | Prefix Residue(s) | Prefix Residue Exception | Final Residue(s) | Final Residue Exception | Cleavage Offset |
|--------|-------------------|--------------------------|------------------|-------------------------|-----------------|
| Trypsin              | KR    | P  | KR    | P  | 1 |
| Modified Trypsin     | KRLNH | na | KRLNH | na | 1 |
| Endoproteinase Glu-C | ED    | na | ED    | na | 1 |
| Endoproteinase Lys-C | K     | na | K     | na | 1 |
| Cyanogen Bromide     | M     | na | M     | na | 1 |
| Proteinase K         | GAVLIMCFW | na | GAVLIMCFWW | na | 1 |
| Elastase, Trypsin, & Chymotrypsin | ALIVKRWFY | P | ALIVKRWFY | P | 1 |
| AspN                 | D     | na | D     | na | 0 |


## Rules for Tryptic state

When determining tryptic state, we define a cleavage site as a K or R residue, not followed by a P residue 

### Fully tryptic peptides 

Cleavage State = 2 
* These peptides have 2 cleavage sites. 
* Note that peptides at the N- or C-terminus of a protein can only be fully tryptic or non-tryptic, never partially tryptic. 
* By definition, we label peptides that span an entire protein as "fully tryptic". 
* Examples:

| Peptide        | Comment                                                 |
|----------------|---------------------------------------------------------|
| `K.ACDEFGR.S`  | Normal, fully tryptic peptide                           |
| `R.ACDEFGR.S`  | Normal, fully tryptic peptide                           |
| `-.ACDEFGR.S`  | Fully tryptic at the N-Terminus of the protein          |
| `R.ACDEFGH.-`  | Fully tryptic at the C-Terminus of the protein          |
| `-.ACDEFG.-`   | Peptide spans the entire protein                        |

* Microsoft Access Like Clause
  * Note that \[!A-Z\] is included to allow for a modified K or R at the C-terminus of a peptide
```sql
Like "[KR].[!P]*[KR].[!P]" Or Like "[KR].[!P]*[KR][!A-Z].[!P]" Or Like "-.*[KR].[!P]" Or Like "[KR].[!P]*.-" Or Like "-.*.-"
```

* SQL Server Like Clause
```sql
WHERE Peptide LIKE '[KR].[^P]%[KR].[^P]' OR 
      Peptide LIKE '[KR].[^P]%[KR][^A-Z].[^P]' OR 
      Peptide LIKE '-.%[KR].[^P]' OR 
      Peptide LIKE '[KR].[^P]%.-' OR 
      Peptide LIKE '-.%.-'
```


### Partially tryptic peptides

Cleavage State = 1 
* These peptides have 1 cleavage site
* Examples:

| Peptide        | Comment                                                 |
|----------------|---------------------------------------------------------|
| `K.ACDEFGH.S`  | Normal, partially tryptic peptide                       |
| `L.ACDEFGR.S`  | Normal, partially tryptic peptide                       |
| `K.ACDEFGR.P`  | Would have been fully tryptic, but ends with R followed by P   |
| `K.PCDEFGR.S`  | Would have been fully tryptic, but starts with K followed by R |

* Microsoft Access Like Clause
  * Note that \[!A-Z\] is included to allow for a modified K or R at the C-terminus of a peptide
```sql
Like "[KR].[!P]*.?" Or Like "?.*[KR].[!P-]" Or Like "?.*[KR][!A-Z].[!P-]"
```

* SQL Server Like Clause
```sql
WHERE Peptide Like '[KR].[^P]%._' Or 
      Peptide Like '_.%[KR].[^P-]' Or 
      Peptide Like '_.%[KR][^A-Z].[^P-]'
```


### Non tryptic peptides

Cleavage State = 0 
* These peptides have no cleavage sites. 
* Examples:

| Peptide        | Comment                                                 |
|----------------|---------------------------------------------------------|
| `L.ACDEFGH.S`  | Normal, non-tryptic peptide                             |
| `-.ACDEFGH.S`  | Normal, non-tryptic peptide that happens to be at the N-terminus |
| `L.ACDEFGH.-`  | Normal, non-tryptic peptide that happens to be at the C-terminus |
| `L.ACDEFGR.P`  | Would have been partially tryptic, but ends with R followed by P |
| `K.PCDEFGR.P`  | Would have been fully tryptic, but has a P after both the K and the R |


### Modified Residues

In any of the above cases, the K or R residue at the end of the peptide may possibly have been modified, like this:

| Peptide        | Comment                                                 |
|----------------|---------------------------------------------------------|
| `K.ACDEFGR*.S` | Normal, fully tryptic peptide, with modified R          |
| `-.ACDEFGR*.S` | Fully tryptic at the N-Terminus of the protein, with modified R |

