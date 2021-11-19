# CaMYBs
The scripts were applied for the investigation of the R2R3 MYBs in _Cicer arietinum_.

## Co-expression analysis

```
Usage:
  python coexp.py --in <FILE> --out <DIR> --exp <FILE>

Mandatory:
  --in  STR         Gene ID input file
  --out STR         Directory for temporary and output files.
  --exp STR         Expression input file
		
Optional:
  --mapping STR     Name mapping table
  --ann     STR     Annotation file
  --c       FLOAT   Min correlation cutoff
  --p       FLOAT   Max pvalue cutoff
  --e       FLOAT   Min expression cutoff
```



`--in` specifies the gene IDs of interest with one ID per line. These IDs need to match the row names in the expression file. Co-expression analysis against all other annotated genes will be performed for these candidates.

`--out` specifies the output folder.

`--exp` specifies the gene expression file. The first column contains the gene IDs. All following columns contain expression data with one sequencing run per column.

`--mapping` specifies a name mapping table to use gene names instead of gene IDs.

`--ann` specifies the annotation file. This should be a table with two columns containing gene IDs in the first one and a functional annotation in the second column.

`--c` specifies the minimal correlation coefficient between the gene expression patterns of gene pairs. Only pairs with a value above this cutoff will be considered.

`--p` specifies the maximal p-value for the correlation between the gene expression patterns of gene pairs. Only pairs with a value below this cutoff will be considered.

`--e` specifies the minimal expression of a gene to be considered in the co-expression analysis.



## Gene statistics calculation

python gene_stats.py
					--candidates
					--gff <GFF_FILE>
					--prot <PROTEIN_FILE>
					--out <STATS_OUTPUT_FILE>

# References

