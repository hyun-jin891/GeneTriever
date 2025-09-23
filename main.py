from GeneTriever import scene_gene_desc_summary, scene_gene_function, scene_cell_relatedTo_gene, scene_pathway, scene_protein_location, p_p_interaction_info, construct_marker_graph
import pandas as pd
import pronto
from dotenv import load_dotenv
import argparse
import obonet
from igraph import Graph

load_dotenv()

## Your Dataset ##
db1 = pd.read_excel("data/Cell_marker_Seq.xlsx")   # CellMarker
db2 = pd.read_csv("data/PanglaoDB_markers_27_Mar_2020.tsv", sep="\t") # PanglaoDB
db3 = pd.read_csv('data/singleCellBase_dataset.txt', sep='\t') # singleCellBase
db4 = pd.read_csv('data/normal_ihc_data.tsv', sep='\t') # HPA: normal_ihc_data.tsv
db5 = pd.read_csv("data/proteinatlas.tsv", sep="\t") # HPA: proteinatlas.tsv
db6 = pd.read_csv("data/interaction_consensus.tsv", sep="\t") # HPA: interaction_consensus.tsv
ont_graph = obonet.read_obo("data/cl.obo") # cl.obo
ont = pronto.Ontology("data/cl-basic.obo") #Your cl-basic.obo downloaded from https://github.com/obophenotype/cell-ontology/releases
db7 = pd.read_csv("data/PCMDB_Download_20250909177169.csv") # plantcellmarker
db8 = pd.read_excel("data/arabidopsis_thaliana.marker_fd.xlsx") # scPlantDB
## Your Dataset ##




def parse_args():
  p = argparse.ArgumentParser(
    description="GeneTriever: "
    )
  
  p.add_argument(
    "-g", "--genes",
    required=True,
    help="Your Genes"
    )
  
  p.add_argument(
    "-s", "--species",
    choices=["human", "mouse", "arabidopsis"],
    help="species (human or mouse or arabidopsis)"
    )
  
  
  return p.parse_args()







######### main() ###############

def main():
  args = parse_args()
  
  genes = args.genes.split(", ")

  res = ""
  
  for gene in genes:
    graph = construct_marker_graph(gene, args.species, ont, db1, db2, db3, db7, db8)
    gene_desc_summary_context, genes_ids = scene_gene_desc_summary(gene, args.species, db5)
    gene_func_summary_context = scene_gene_function(gene, args.species, genes_ids, db5)
    gene_cell_summary_context = scene_cell_relatedTo_gene(graph, gene, args.species, genes_ids, db5)
    gene_pathway_context = scene_pathway(gene, genes_ids, db5, args.species)
    gene_protein_loc_context = scene_protein_location(gene, genes_ids, args.species)
    ppi_context = p_p_interaction_info(gene, args.species, db6)
    
    total_context = gene_desc_summary_context + "\n" + gene_func_summary_context  + "\n" + gene_cell_summary_context + "\n" + gene_pathway_context + "\n" + gene_protein_loc_context + "\n" + ppi_context + "\n"
    
    res += total_context + "\n\n"
    
  f = open("GeneRetriever_result.txt", "w")
  f.write(res)
    
  f.close()

  
################################  


if __name__ == "__main__":
  main()
  
