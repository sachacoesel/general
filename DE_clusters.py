#!/usr/bin/env python

# This script will take in a file containing clusters and a list of DE genes
# will examine each cluster in whole and report average logFC, ratio of directionality
# number of significant files, etc.

# Ryan Groussman and Gwenn MM Hennon, February 2016


# on CEG hallway computer:
# C:\Users\ceg\Ryan\Dropbox\Armbrust\cAMP\for Gwenn paper\GEO\GSE67971_Cluster_results.csv
# tiny test clusters:
# C:\Users\ceg\Ryan\Dropbox\Armbrust\cAMP\for Gwenn paper\GEO\tinyTEST_GSE67971_Cluster_results.csv
# tiny test DE genes:
# C:\Users\ceg\Ryan\Dropbox\Armbrust\cAMP\for Gwenn paper\GEO\tinyTEST_Low-HighCO2.plusIBMX.UQ_fdr.tab

from numpy import std
import math

# paths for the things (on bloom)
cluster_path = "/home/rgrous83/cAMP/bifx/DE_clusters/GSE67971_Cluster_results.csv" # the real deal
# cluster_path = "/home/rgrous83/cAMP/bifx/DE_clusters/tinyTEST_GSE67971_Cluster_results.csv" # small test batch

# DGE_path = "/home/rgrous83/cAMP/bifx/DE_clusters/UQ_DE_counts/Low-High_CO2.minusIBMX.UQ_fdr.tab"
DGE_path = "/home/rgrous83/cAMP/bifx/DE_clusters/UQ_DE_counts/Low-HighCO2.plusIBMX.UQ_fdr.tab"
# DGE_path = "/home/rgrous83/cAMP/bifx/DE_clusters/UQ_DE_counts/Minus-Plus_IBMX.highCO2.UQ_fdr.tab"
# DGE_path = "/home/rgrous83/cAMP/bifx/DE_clusters/UQ_DE_counts/Minus-Plus_IBMX.lowCO2.UQ_fdr.tab"

# open pertinent files
cluster_file = open(cluster_path, 'r')
DGE_file = open(DGE_path, 'r')

# out_path = "DE_clusters.test.csv"

# out_path = "DE_clusters2.Low-High_CO2.minusIBMX.csv"
out_path = "DE_clusters2.Low-HighCO2.plusIBMX.csv"
# out_path = "DE_clusters2.Minus-Plus_IBMX.highCO2.csv"
# out_path = "DE_clusters2.Minus-Plus_IBMX.lowCO2.csv"

outfile = open(out_path,'w')

def parse_clusters(cluster_file):
	'''Will take as input a file with cluster results, as in GSE67971_Cluster_results.csv
	and return a dictionary, with cluster_id as key, and list of gene ids as value
	'''
	cluster_dict = dict()
	cluster_file.readline() # skip the first line
	for line in cluster_file:
		cluster_elts = line.split(",")
		cluster_id = cluster_elts[0]
		protein_ids = cluster_elts[13].split(";")
		cluster_dict[cluster_id] = protein_ids
	return cluster_dict


def process_DE_scores(DGE_file):
	'''This will take as input an output file from edgeR with
	gene names, logFC, logCPM, and fdr. Returns a dictionary
	with gene names as keys, and logFC and fdr in a list of values
	'''
	DE_dict = dict()
	DGE_file.readline() # skip the first line
	for line in DGE_file:
		DGE_elts = line.split("\t") # tab separated values
		gene_id = DGE_elts[0]
		logFC = float(DGE_elts[1])
		fdr = float(DGE_elts[5].strip())
		DE_dict[gene_id] = [logFC, fdr]
	return DE_dict


def get_DE_scores(DE_dict, gene_id):
	'''For a given gene_id in a DE_dict given by process_DE_scores,
	will return a list with the logFC and fdr values for this gene
	'''
	logFC = DE_dict[gene_id][0]
	fdr = DE_dict[gene_id][1]
	return [logFC, fdr]


def cluster_stats(DE_dict, cluster):
	'''Taking as input the name of one of the clusters in cluster_dict,
	will evaluate each gene in the cluster for significance, and then
	calculate the average FC, stdev of FC, ratio of [sig genes/total genes],
	and the number of total significant genes in the cluster.
	Will output these results as a list'''
	logFC_list = []
	sig_genes = 0.0
	up_genes = 0
	down_genes = 0
	cluster_genes = cluster_dict[cluster] # cluster genes may have 'bd' genes not represented in Thaps3 build; we will ignore these
	bd_genes = [] # for holding unsupported 'bd' genes
	for gene in cluster_genes:
		if gene not in DE_dict.keys():
			bd_genes.append(gene)
	for gene in bd_genes:
		cluster_genes.remove(gene) # remove from cluster_genes and add to bd_genes
	total_genes = len(cluster_genes)
	for gene_id in cluster_genes:
		DE_scores = get_DE_scores(DE_dict, gene_id)
		if DE_scores[1] < 0.05:  # .05 is our arbitrary fdr cutoff for significance
			sig_genes += 1 # increment if fdr < 0.05
			logFC_list.append(DE_scores[0]) # add logFC the list
			if DE_scores[0] > 0:
				up_genes += 1
			elif DE_scores[0] < 0:
				down_genes += 1
	if sig_genes > 0:
		averageFC = sum(logFC_list)/sig_genes
		directionality = (up_genes - down_genes) / sig_genes
		log_sig_genes = math.log(sig_genes, 2)
	elif sig_genes == 0:
		averageFC = 0.0
		directionality = 0
		log_sig_genes = 0
	stdevFC = std(logFC_list) # this is for the population st dev - correct in this situation or not?
	if total_genes > 0:
		sig_ratio = sig_genes / total_genes
	else:
		sig_ratio = 0
	cluster_rank = log_sig_genes * sig_ratio * averageFC * directionality
	gene_list = ';'.join(cluster_genes)
	return_list = [cluster, sig_genes, sig_ratio, averageFC, stdevFC, up_genes, down_genes, directionality, cluster_rank, gene_list]
	return return_list

# the log changes has been implemented and removed; not sure I liked it.
# implementation
cluster_dict = parse_clusters(cluster_file)
DE_dict = process_DE_scores(DGE_file)

outheader = "cluster,sig_genes,sig_ratio,avgFC,stdevFC,up_genes,down_genes,directionality,cluster_rank,gene_list\n"
outfile.write(outheader)
for cluster in cluster_dict.keys():
	clust_list = cluster_stats(DE_dict, cluster)
	for elt in clust_list:
		outfile.write(str(elt) + ",")
	outfile.write("\n")


# outputs a csv table with cluster, sig_genes, ratio_sig_genes, avg_FC, stdev_FC, up_genes, down_genes]
