### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python gene_stats.py
					--candidates
					--gff <GFF_FILE>
					--prot <PROTEIN_FILE>
					--out <STATS_OUTPUT_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import re, sys, os

# --- end of imports --- #

def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split('.')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:].split('.')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def load_mapping_table( candidate_gene_file ):
	"""! @brief load name to ID mapping table """
	
	mapping_table = {}
	with open( candidate_gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			mapping_table.update( { parts[1]: parts[0] } )
			line = f.readline()
	return mapping_table


def load_gff_infos( gff_file ):
	"""! @brief load number of exons per gene from GFF3 file """
	
	# --- load data --- #
	rna_to_gene = {}
	exons_per_gene = {}
	CDS_per_gene = {}
	gene_pos = {}
	with open( gff_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] in [ "gene", "pseudogene" ]:
					try:
						ID = re.findall( "gene\d+", parts[-1] )[0]
						exons_per_gene.update( { ID: 0 } )
						CDS_per_gene.update( { ID: 0 } )
						gene_pos.update( { ID: { 'start': parts[3], 'end': parts[4], 'chr': parts[0], 'orientation': parts[6] } } )
					except:
						pass
				elif parts[2] in [ "mRNA", "tRNA", "rRNA", "lnc_RNA", "transcript" ]:
					try:
						ID = re.findall( "rna\d+", parts[-1] )[0]
						parent = re.findall( "gene\d+", parts[-1] )[0]
						rna_to_gene.update( { ID: parent } )
					except:
						print parts[-1]
				elif parts[2] == "exon":
					try:
						rna = re.findall( "rna\d+", parts[-1] )[0]
						exons_per_gene[ rna_to_gene[ rna ] ] += 1
					except IndexError:
						gene = re.findall( "gene\d+", parts[-1] )[0]
						exons_per_gene[ gene ] += 1
				elif parts[2] == "CDS":
					try:
						rna = re.findall( "rna\d+", parts[-1] )[0]
						CDS_per_gene[ rna_to_gene[ rna ] ] += 1
					except:
						gene = re.findall( "gene\d+", parts[-1] )[0]
						CDS_per_gene[ gene ] += 1
			line = f.readline()
	
	# --- assign data to gene ID --- #
	infos = {}
	for ID in gene_pos.keys():
		try:
			infos.update( { ID: { 	'exons': exons_per_gene[ ID ],
											'cds': CDS_per_gene[ ID ],
											'start': gene_pos[ ID ]['start'],
											'end': gene_pos[ ID ]['end'],
											'chr': gene_pos[ ID ]['chr'],
											'orientation': gene_pos[ ID ]['orientation']
										 } } )
		except KeyError:
				print ID
	return infos


def main( arguments ):
	"""! @brief run everything """
		
	candidate_gene_file = arguments[ arguments.index('--candidates')+1 ]
	gff_file = arguments[ arguments.index('--gff')+1 ]
	prot_file = arguments[ arguments.index('--prot')+1 ]
	stats_output_file =  arguments[ arguments.index('--out')+1 ]

	prots = load_sequences( prot_file )
	name_ID_mapping_table = load_mapping_table( candidate_gene_file )
	gene_infos = load_gff_infos( gff_file )


	with open( stats_output_file, "w" ) as out:
		out.write( "\t".join( [ "ID", "name", "chr", "start", "end", "orientation", "amino acids", "exons", "coding exons" ] ) + '\n' )
		for ID in name_ID_mapping_table.keys():
			try:
				out.write( "\t".join( map( str, [	ID,
																	name_ID_mapping_table[ ID ],
																	gene_infos[ID]['chr'],
																	gene_infos[ID]['start'],
																	gene_infos[ID]['end'],
																	gene_infos[ID]['orientation'],
																	len( prots[ name_ID_mapping_table[ ID ] ] ),
																	gene_infos[ID]['exons'],
																	gene_infos[ID]['cds']
																] )) + '\n' )
			except KeyError:
				print ID


if '--candidates' in sys.argv and '--gff' in sys.argv and '--prot' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
