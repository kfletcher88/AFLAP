# AFLAP

## Custom pedigree file.
The pedigree file outlines the pedigree of the cross to be analyzed and the location of the reads which are going to be analyzed. This is done by three mandatory tab delimited fields, two optional.
Field 1 is the label which will be retained in all downstream analysis.
Field 2 is the generation, parent (F0) = 0, child (F1) = 1, grandchild (F2) = 2. This field will inform AFLAP what type of analysis is to be run.
Field 3 is the location of the read file.

One individual may be split over multiple lines in instances where multiple read files are present. AFLAP will combine lines that share identical labels (field 1).

Optional fields.
For parents (F0), field 4 and 5 can be used to manually set the k-mer cutoffs forwarded to assembly. If these fields are not set, AFLAP will try to estimate these cutoffs, which may not be perfect.
For progeny (F1 and F2), field 4 and 5 can be used to specify the parents. Not currently implemented, but is foundational for a potential multi-cross analysis enhancement.

An example pedigree is available in:
```
./example/Pedigree.txt
```

It looks like this:
```
#label, F, Read1, Read2, OptionalField1, OptionalField2
Col	0	ReadsEBI/SRR5882797_1.fastq.gz	32	105
Col	0	ReadsEBI/SRR5882797_2.fastq.gz	32	105
Ler	0	ReadsEBI/SRR3166543_1.fastq.gz	40	177
Ler	0	ReadsEBI/SRR3166543_2.fastq.gz	40	177
ERR1432507	2	ReadsEBI/ERR1432507_1.fastq.gz	Col	Ler
ERR1432507	2	ReadsEBI/ERR1432507_2.fastq.gz	Col	Ler
ERR1432492	2	ReadsEBI/ERR1432492_1.fastq.gz	Col	Ler
ERR1432492	2	ReadsEBI/ERR1432492_2.fastq.gz	Col	Ler
```

