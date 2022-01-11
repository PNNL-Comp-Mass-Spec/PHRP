# Peptide Hit Results Processor

The Peptide Hit Results Processor (PHRP) converts search results from various
MS/MS identification tools into a series of tab-delimited text files
that organize the data in a similar format for each tool. It supports
MS-GF+, MaxQuant, MSFragger, MODa, MODPlus, MSAlign, MSPathFinder,
TopPIC, and X!Tandem, along with SEQUEST Synopsis/First Hits files.

PHRP will insert modification symbols into the peptide sequences for modified peptides.
Parallel files are created containing sequence information, modification details,
and protein information. The user can optionally provide a modification definition
file that specifies the symbol to use for each modification mass.

## Example Data

Example input and output files are in the Data directory:
* MSGFPlus_Example has MS-GF+ results
* XTandem_Example has X!Tandem results

For [MS-GF+](https://github.com/sangtaekim/msgfplus) results, prior to running PHRP, use the
[Mzid-To-Tsv-Converter](https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter)
to convert the .mzid file to a tab-delimited .tsv file.

## Example Command line

```
PeptideHitResultsProcRunner.exe /I:Dataset_msgfplus.tsv /M:MSGFPlus_PartTryp_MetOx_20ppmParTol_ModDefs.txt /N:MSGFPlus_PartTryp_MetOx_20ppmParTol.txt /T:Mass_Correction_Tags.txt /L /ProteinMods /F:H_sapiens_UniProt_SPROT_2021-03-07.fasta
```

## Command Line Syntax

PHRP is a console application, and must be run from the Windows command prompt.
* On Linux, use [Mono](https://www.mono-project.com/download/stable/#download-lin) to run the program,
for example `mono PeptideHitResultsProcRunner.exe`

```
PeptideHitResultsProcRunner.exe
 InputFilePath [/O:OutputDirectoryPath]
 [/S:[MaxLevel]]
 [/P:ParameterFilePath]
 [/L:[LogFilePath]] [/LogDir:LogDirectoryPath]
 [/T:MassCorrectionTagsFilePath]
 [/M:ModificationDefinitionFilePath]
 [/N:SearchToolParameterFilePath]
 [/F:FastaFilePath]
 [/CreateModSummaryFile:True|False]
 [/ProteinMods] [/ProteinModsIncludeReversed]
 [/UseExistingPepToProteinMapFile]
 [/CreateProteinModsViaPHRP]
 [/IgnorePepToProtMapErrors]
 [/FHT:True|False] [/Syn:True|False]
 [/MaxQScore:50] [/MaxQPEP:0.01]
 [/MSGFPlusSpecEValue:0.0000005] [/MSGFPlusEValue:0.75]
 [/SynProb:0.05] [/SynPValue:0.95]
 [/DB:DatabaseConnectionString]
```

The input file should be one of the following:
* MS-GF+ results file (_msgfplus.tsv or _msgfdb.tsv or .tsv)
* MSPathFinder results file (_IcTda.txt)
* MaxQuant results files (msms.txt and peptides.txt)
* MSFragger results file (Dataset_psm.tsv or Dataset.tsv)
* MODa results file (_moda.id.txt)
* MODPlus results file (_modp.id.txt)
* MSAlign results file (\_MSAlign_ResultTable.txt)
* TopPIC results file (\_TopPIC_PrSMs.txt)
* X!Tandem Results file (_xt.xml)
* Legacy tools
  * Inspect results file (_inspect.txt)
  * MSGF-DB results file (_msgfdb.txt)
  * SEQUEST Synopsis File (_syn.txt)
  * SEQUEST First Hits file (_fht.txt)

Optionally use `/O` to specify the output directory
* If omitted, the output file will be created in the same directory as the input file

Use `/S` to process all valid files in the input directory and subdirectories
* Include a number after `/S` (like `/S:2`) to limit the level of subdirectories to examine

Optionally use `/P` to supply a parameter file with processing options
* This must be a Key=Value parameter file, for example:
  * [PHRP_Options.conf](https://github.com/PNNL-Comp-Mass-Spec/PHRP/blob/master/Data/PHRP_Options.conf)

Use `/L` to specify that a log file should be created
* Use `/L:LogFilePath` to specify the name (or full path) for the log file
* Use /LogDir to specify the directory to create the log file (ignored if the LogFilePath is rooted)

Use `/T` to specify the file containing the mass correction tag info
* This file should be tab delimited, with the first column containing the mass correction tag name and the second column containing the mass (the name cannot contain commas or colons)

Use `/M` to specify the file containing the modification definitions.
* This file should be tab delimited, with the first column containing the modification
symbol, the second column containing the modification mass, plus optionally a
third column listing the residues that can be modified with the given mass
* For the third column, use 1 letter residue symbols (no need to separate with commas or spaces)

Use `/N` to specify the parameter file provided to the search tool
* This is used when processing results from MS-GF+, MSPathFinder, MaxQuant, MSFragger, MODa, MODPlus, MSAlign, TopPIC, and InSpecT
* For MaxQuant, provide either an XML-based parameter file (root element is <MaxQuantParams>) or provide the parameters.txt file created in the txt results directory
  * The XML-based parameter file is preferred, since it is required to allow PHRP to accurately compute monoisotopic masses of peptides identified by MaxQuant

Use `/F` to specify the path to the FASTA file
* When provided, the order of the proteins in the FASTA file dictates which protein is listed for each peptide in the First Hits file

A modification summary file with extension _ModSummary.txt is created by default
* To disable, use `/CreateModSummaryFile:False`

Use `/ProteinMods` to indicate that the _ProteinMods.txt file should be created
* This requires that either an existing _PepToProtMapMTS.txt file exist, or that the FASTA file be defined using `/F`

Use `/ProteinModsIncludeReversed` if an existing _ProteinMods.txt file has reversed protein sequences, or if the FASTA file has reversed proteins
* If false, will skip reversed proteins when creating the _ProteinMods.txt file

Use `/UseExistingPepToProteinMapFile` to use an existing _PepToProtMapMTS.txt file if it exists

Use `/CreateProteinModsViaPHRP` to indicate that InputFilePath specifies a valid PHRP data file and thus the PHRP data files should not be re-created; only the
_ProteinMods.txt file should be created
* This requires that either an existing _PepToProtMapMTS.txt file exist, or that the FASTA file be defined using `/F`

Use `/IgnorePepToProtMapErrors` to ignore peptide to protein mapping errors that occur when creating a missing _PepToProtMapMTS.txt file

Use `/FHT:True` or `/FHT:False` to control the creation of a first-hits file (_fht.txt)
* The default is `/FHT:True`

Use `/Syn:True` or `/Syn:False` to toggle the creation of a synopsis file (_syn.txt)
* The default is `/Syn:True`

When processing a MaxQuant results file, use `/MaxQScore` to customize the Andromeda score threshold used to determine which peptides are written to the synopsis file
* A PSM is stored if its Andromeda score is over the threshold, or if its PEP score is below the threshold
* The default is `/MaxQScore:50`

When processing a MaxQuant results file, use `/MaxQPEP` to customize the Posterior Error Probability (PEP) score threshold used to determine which peptides are written to the synopsis file
* The default is `/MaxQPEP:0.01`

When processing an MS-GF+ results file, use `/MSGFPlusSpecEValue` and `/MSGFPlusEValue` to customize the thresholds used to determine which peptides are written to the synopsis file
* Defaults are `/MSGFPlusSpecEValue:5E-07` and `/MSGFPlusEValue:0.75`

When processing a MODPlus or MODa results file, use `/SynProb` to customize the probability threshold used to determine which peptides are written to the synopsis file
* The default is `/SynProb:0.05`
* Higher probability values are higher confidence results
  * 0.05 is a very loose filter

When processing a MODPlus or MODa results file, use `/SynPValue` to customize the p-value threshold used to determine which peptides are written to the synopsis file
* The default is `/SynPValue:0.95`
* Lower p-values are higher confidence results
  * 0.95 is a very loose filter

When processing MaxQuant results using a computer on the pnl.gov domain, the DMS database is contacted to lookup dataset IDs by dataset name, where dataset name comes from the 'Raw file' column in the msms.txt file
* Optionally use `/DB` to override the default connection info
* The default is `/DB:"Data Source=gigasax;Initial Catalog=DMS5;User=DMSReader;Password=dms4fun"`
* To disable contacting DMS, use `/DB:""` or `/DB:false`

## Contacts

Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) \
E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov \
Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics

## License

The Peptide Hit Results Processor is licensed under the 2-Clause BSD License;
you may not use this program except in compliance with the License. You may obtain
a copy of the License at https://opensource.org/licenses/BSD-2-Clause

Copyright 2018 Battelle Memorial Institute
