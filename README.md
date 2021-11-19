# CaMYBs
The scripts were applied for the investigation of the R2R3 MYBs in _Cicer arietinum_.


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





python gene_stats.py
					--candidates
					--gff <GFF_FILE>
					--prot <PROTEIN_FILE>
					--out <STATS_OUTPUT_FILE>

# References

