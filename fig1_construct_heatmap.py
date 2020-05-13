import re, math
from pandas import DataFrame
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from operator import itemgetter
import numpy as np

# --- end of imports --- #

def load_expression_values( filename ):
	"""! @brief load all expression values """
	
	expression_data = {}
	with open( filename, "r" ) as f:
		raw_tissues = f.readline().strip().split('\t')[1:]
		tissues = []
		for each in raw_tissues:
			tissues.append( each.replace( "_1", "" ) )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression = {}
			for idx, each in enumerate( parts[1:] ):
				expression.update( { tissues[  idx ] : float( parts[ idx+1 ] ) } )
			line = f.readline()
			expression_data.update( { parts[0]: expression } )
	return expression_data


def construct_data_output_file( data, candidate_genes, candidate_samples, outputfile, cutoff, gene_name_mapping_table ):
	"""! @brief write expression values of all candidate genes into output file """
	
	datamatrix = []
	genes = []
	with open( outputfile, "w" ) as out:
		new_line = [ "gene" ] + sorted( candidate_samples.keys() )
		tissues = new_line[1:]
		out.write( "\t".join( new_line ) + '\n' )
		for gene in sorted( candidate_genes ):
			new_line = [ gene ]
			for tissue in tissues:
				tmp_value = []
				for sample in candidate_samples[ tissue ]:
					try:
						tmp_value.append( data[ gene ][ sample ] )
					except KeyError:
						print sample
				if len( tmp_value ) > 0:
					new_line.append( sum( tmp_value ) / len( tmp_value ) )
				else:
					new_line.append( 0 )
			if sum( new_line[1:] ) > cutoff:
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
				datamatrix.append( new_line[1:] )
				genes.append( gene_name_mapping_table[ gene ] )
	return genes, tissues, datamatrix


def construct_heatmap( datamatrix, genes, tissues, heatmap_file ):
	"""! @brief construct heatmap from given data matrix """
	
	my_vmax = 100
	
	print "number of genes for heatmap construction: " + str( len( genes ) )
	print "number of samples for heatmap construction: " + str( len( tissues ) )
	
	df = DataFrame( datamatrix, index=genes[::-1], columns=tissues).round( 0 )
	
	fig, ax = plt.subplots( figsize=(3,7) )
	
	sns.heatmap( df, vmin=0, vmax= my_vmax, ax=ax, linewidths=0.3, annot=True, annot_kws={'fontsize':3}, cbar=False, fmt="g", cmap='binary' )	#cmap='YlGnBu'  = 1
	
	for idx, gene in enumerate( genes ):
		ax.text( -3, idx+0.6, gene, fontsize=3 )
	
	for idx, tissue in enumerate( tissues ):
		ax.text( idx+0.4, len( genes )+0.5, tissue, rotation=90, fontsize=3 )
	
	ax.set_yticklabels( [], rotation=0, fontsize=2 )
	ax.set_xticklabels( [] , rotation=90, fontsize=3  )
	
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	
	plt.yticks( rotation=0 )
	plt.subplots_adjust( left=0.15, right=0.99, top=0.99, bottom=0.08, wspace=0.2 )
	
	plt.savefig( heatmap_file, dpi=900  )
	plt.savefig( heatmap_file.replace( ".png", ".svg" ) )


if __name__ == '__main__':
	
	candidate_gene_file = "candidate_genes.txt"
	candidate_sample_file = "samples_of_interest_fig1b.txt"
	
	kabuli_expression_file = "kabuli_FPKMs.txt"
	desi_expression_file = "desi_FPKMs.txt"
	
	cutoff = 0 #minimal cumulative expression of gene in all tissues combined
	prefix = "MYB_gene_exp_heatmaps/"
	
	name = "fig1b"
	
	# --- load data --- #
	expression_data = load_expression_values( kabuli_expression_file )
	expression_data.update( load_expression_values( desi_expression_file ) )
		
	candidate_genes = []
	gene_name_mapping_table = {}
	with open( candidate_gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[1] not in candidate_genes:
				candidate_genes.append( parts[1] )
				gene_name_mapping_table.update( { parts[1]: parts[0] } )
			line = f.readline()
	
	candidate_samples = {}
	with open( candidate_sample_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			candidate_samples.update( { parts[0]: parts[1].split(',') } )
			line = f.readline()
	
		
	# --- write gene expression values of candidates into data output file --- #
	outputfile = prefix + name + "plotted_values.txt"
	genes, tissues, datamatrix = construct_data_output_file( expression_data, candidate_genes, candidate_samples, outputfile, cutoff, gene_name_mapping_table )
	
	heatmap_file = prefix + name + ".png"
	construct_heatmap( datamatrix, genes, tissues, heatmap_file )
	
	print "all done!"
