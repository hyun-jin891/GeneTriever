import requests
from langchain.schema import Document
import pandas as pd
from celltype_quality_check import quality_check
import mygene
import io, contextlib
from dotenv import load_dotenv
from openai import OpenAI
from xml.etree import ElementTree
from igraph import Graph

load_dotenv()


def construct_marker_graph(marker, species, ont, db1, db2, db3, db7, db8):  
  
  graph = Graph(directed=False)
  
  graph.add_vertices(1)
  
  graph.vs[0]["name"] = marker
  graph.vs[0]["type"] = "Marker"
  celltypes = []
  celltypes_2 = []
  celltypes_3 = []
    
  if species == "arabidopsis":
    filtered_db7 = db7.loc[db7["Gene_symbol"].str.lower() == marker.lower()][:]
    celltypes = filtered_db7["Cell_type"].unique().tolist()
      
    filtered_db8 = db8.loc[db8["name"].str.lower() == marker.lower()][:]
    celltypes_2 = filtered_db8["clusterName"].unique().tolist()
      
  else:
    filtered_db1 = db1.loc[(db1["species"].str.lower() == species.lower()) & (db1["marker"].str.lower() == marker.lower())][:]
    celltypes = filtered_db1["cell_name"].unique().tolist()
    
    species_for_db2 = {"human":"hs", "mouse":"mm"}
    filtered_db2 = db2.loc[((species_for_db2[species.lower()] == db2["species"].str.lower()) | (db2["species"].str.lower() == "mm hs")) & (db2["official gene symbol"].str.lower() == marker.lower())][:]
    celltypes_2 = filtered_db2["cell type"].unique().tolist()
    
    
  species_for_db3 = {"human":"Homo sapiens", "mouse":"mouse", "arabidopsis":"Arabidopsis thaliana"}
  db3["extra"] = db3["gene_symbol"].str.lower().str.split(r",\s*")
  db3["extra"] = db3["extra"].apply(lambda l: l if isinstance(l, list) else [])
    
        
  filtered_db3 = db3.loc[(db3["species"].str.lower() == species_for_db3[species.lower()].lower()) & (db3["extra"].apply(lambda l: marker.lower() in l))][:]
  celltypes_3 = filtered_db3["cell_type"].unique().tolist()

  celltypes.extend(celltypes_2)
  celltypes.extend(celltypes_3)
  celltypes = list(set(celltypes))
    
    
  if species == "human":
    celltypes = list(set(quality_check(celltypes, ont)))
    
    
  for j in range(len(celltypes)):
    celltype = celltypes[j]
    if type(celltype) != str or celltype == "":
      continue
    graph.add_vertices([celltype])
    graph.vs[-1]["type"] = "Cell"    
    graph.add_edge(marker, celltype)
  
  
  return graph



def get_info_proteinAtlas(marker, db5, query):
  mask = db5['Gene'] == marker
  record = db5[mask]
  
  doc = ""
                                                                                                                                                                                               
  
  if query == "Gene description" and len(record["Gene description"]) != 0 and type(record["Gene description"].item()) == str:
    doc = f"{marker}: " + record["Gene description"].item()
  if query == "Protein class" and len(record["Protein class"]) != 0 and type(record["Protein class"].item()) == str:
    doc = f"{marker}'s protein class: " + record["Protein class"].item()
  if query == "Biological process" and len(record["Biological process"]) != 0 and type(record["Biological process"].item()) == str:
    doc = f"{marker}'s biological process: " + record["Biological process"].item()
  if query == "Molecular function" and len(record["Molecular function"]) != 0 and type(record["Molecular function"].item()) == str:
    doc = f"{marker}'s molecular function: " + record["Molecular function"].item()
  if query == "Tissue specificity" and len(record["RNA tissue specificity"]) != 0 and type(record["RNA tissue specificity"].item()) == str:
    doc += f"{marker}'s RNA tissue specificity: " + record["RNA tissue specificity"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA tissue distribution"]) != 0 and type(record["RNA tissue distribution"].item()) == str:
    doc += f"{marker}'s RNA tissue distribution: " + record["RNA tissue distribution"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA single cell type specificity"]) != 0 and type(record["RNA single cell type specificity"].item()) == str:
    doc += f"{marker}'s RNA single cell type specificity: " + record["RNA single cell type specificity"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA single cell type distribution"]) != 0 and type(record["RNA single cell type distribution"].item()) == str:
    doc += f"{marker}'s RNA single cell type distribution: " + record["RNA single cell type distribution"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA cancer specificity"]) != 0 and type(record["RNA cancer specificity"].item()) == str:
    doc += f"{marker}'s RNA cancer specificity: " + record["RNA cancer specificity"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA cancer distribution"]) != 0 and type(record["RNA cancer distribution"].item()) == str:
    doc += f"{marker}'s RNA cancer distribution: " + record["RNA cancer distribution"].item() + "\n"
  
  return doc


def get_uniprot_id(gene, species):
    url = "https://rest.uniprot.org/uniprotkb/search"
    
    if species == "mouse":
      species = "mus musculus"
    if species == "arabidopsis" or species == "Arabidopsis thaliana":
      species = "arabidopsis thaliana"
    
    tax_ids = {"human":9606, "mus musculus":10090, "arabidopsis thaliana":3702}
    
    tax_id = tax_ids[species.lower()]
    
    params = {
        "query": f"gene_exact:{gene} AND organism_id:{tax_id}",
        "fields": "accession,gene_primary,organism_name",
        "format": "json",
        "size": 1
    }
    
    try:
      r = requests.get(url, params=params)
      data = r.json()
    except Exception as e:
      return None

    if data.get("results"):
        result = data["results"][0]
        return result["primaryAccession"]
    else:
        return None



def scene_gene_desc_summary(marker, species, db5):
  gene_desc_context = "[Gene's Description & Summary]\n"
  llm = OpenAI()
  if species.lower() == "arabidopsis":
    species = "Arabidopsis thaliana"
  genes_ids = {}
  

  gene_desc_context += f"{marker}:\n"
    
  esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
  params = {
       "db": "gene",
       "term": f"{marker}[gene] AND {species}[orgn]",
       "retmode": "json"
  }
    
  ncbi_desc = ""
  ncbi_summary = ""
    
  try:
      res = requests.get(esearch_url, params=params).json()
      gene_id = res["esearchresult"]["idlist"][0]
      genes_ids[marker] = gene_id
    
      efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
      efetch = requests.get(efetch_url, params={"db":"gene", "id":gene_id, "retmode":"xml"})
    
      root = ElementTree.fromstring(efetch.content)
      
      if species == "Arabidopsis thaliana":
        desc = root.find(".//Prot-ref_desc")
        summary = root.find(".//Entrezgene_summary")
      else:
        desc = root.find(".//Gene-ref_desc")
        summary = root.find(".//Entrezgene_summary")

      try:
        ncbi_desc = desc.text
      except Exception as e:
        pass
    
      try:
        ncbi_summary = summary.text
      except Exception as e:
        pass
  except Exception as e:
      pass
    
  proteinAtlas_desc = ""
  proteinAtlas_class = ""
  if species == "human":  
    proteinAtlas_desc = get_info_proteinAtlas(marker, db5, "Gene description")
    proteinAtlas_class = get_info_proteinAtlas(marker, db5, "Protein class")
    
  uniprot_id = get_uniprot_id(marker, species)
  r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
  uniprot_docs = r.json()
    
  uniprot_desc = ""
    
  if uniprot_docs == None:
      pass
  else:
      try:
        protein_description = uniprot_docs['proteinDescription']
        if "recommendedName" in protein_description:
          uniprot_desc = f"{marker} encodes " + protein_description['recommendedName']['fullName']['value']

        else:
          uniprot_desc = f"{marker} encodes " + protein_description['submissionNames'][0]['fullName']['value']
    
      except Exception as e:
        pass
    
  gene_desc_summary_context = ncbi_desc + "\n" + ncbi_summary + "\n" + proteinAtlas_desc + "\n" + proteinAtlas_class + "\n" + uniprot_desc + "\n"
  
  prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the description of the specific gene. You should remove duplicated contents and summarize it.
    
    [Gene]
    {marker}
    
    [Context]
    {gene_desc_summary_context} 
    
    """
  messages = [{"role": "user", "content": prompt}]
    
  completion = llm.chat.completions.create(
        model = "gpt-4o",
        messages = messages,
        temperature = 0,
  )
    
  gene_desc_summary_context = completion.choices[0].message.content
  
  gene_desc_summary_context = gene_desc_context + "\n" + gene_desc_summary_context + "\n"
  
  return gene_desc_summary_context, genes_ids
    
    
    
def scene_gene_function(marker, species, genes_ids, db5):
  gene_func_context = "[Gene's function]\n"
  llm = OpenAI()
  
  gene_func_context += f"{marker}:\n"
    
  mygene_MF_text = ""
    
  try:
      gene_id = genes_ids[marker]
      fields = "go"
      mygene_url = f"https://mygene.info/v3/gene/{gene_id}?fields={fields}"
      mygene_response = requests.get(mygene_url)
      go_text = mygene_response.text
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene's biological process(BP), molecular function(MF), and cellular component(CC). You should extract MF contents only and summarize it.
    
    [Gene]
    {marker}
    
    [Context]
    {go_text}

    """
      messages = [{"role": "user", "content": prompt}]
    
      completion = llm.chat.completions.create(
        model = "gpt-4o",
        messages = messages,
        temperature = 0,
        )
    
      mygene_MF_text = completion.choices[0].message.content
  except Exception as e:
      pass
    
  proteinAtlas_mf = ""
  if species == "human":
    proteinAtlas_mf = get_info_proteinAtlas(marker, db5, "Molecular function")
    

  uniprot_id = get_uniprot_id(marker, species)
  r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
  uniprot_docs = r.json()
  uniprot_func = ""
    
  if uniprot_docs == None:
      pass
  else:
      try:

        protein_comments = uniprot_docs['comments']
        
        for i in range(len(protein_comments)):
          curPos = protein_comments[i]
          
          if curPos["commentType"] == "FUNCTION":
            for j in curPos['texts']:
              uniprot_func += f"{marker}'s function: " + j['value'] + "\n"
          elif curPos["commentType"] == "CATALYTIC ACTIVITY":
              uniprot_func += f"{marker}'s catalytic activity: " + curPos['reaction']['name'] + "\n"
              
      except Exception as e:
        pass
    
  gene_func_summary_context = mygene_MF_text + "\n" + proteinAtlas_mf + "\n" + uniprot_func + "\n"
  
  prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the molecular function of the specific gene. You should remove duplicated contents and summarize it.
    
    [Gene]
    {marker}
    
    [Context]
    {gene_func_summary_context} 
    
    """
  messages = [{"role": "user", "content": prompt}]
    
  completion = llm.chat.completions.create(
        model = "gpt-4o",
        messages = messages,
        temperature = 0,
  )
    
  gene_func_summary_context = completion.choices[0].message.content
  
  gene_func_summary_context = gene_func_context + "\n" + gene_func_summary_context + "\n"
  
  return gene_func_summary_context


  

def scene_cell_relatedTo_gene(graph, marker, species, genes_ids, db5):
  gene_cell_summary_context = "[Cell types or Tissue related to Gene]\n"

  gene_cell_summary_context += f"{marker}:\n"
    
  cell_marker_text = f"{marker} was reported to be expressed in these kinds of cell types (There may be more cell types that {marker} can be expressed.): "
    
  neighbors = graph.neighbors(marker)

  neighbors_celltypes = [graph.vs[idx]["name"] for idx in neighbors]

  cell_marker_text += " / ".join(neighbors_celltypes)
    
  proteinAtlas_tissue = ""
    
  if species == "human":
    proteinAtlas_tissue = get_info_proteinAtlas(marker, db5, "Tissue specificity")
    
    
  uniprot_id = get_uniprot_id(marker, species)
  r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
  uniprot_docs = r.json()
  uniprot_tissue = ""
    
  if uniprot_docs == None:
      pass
  else:
      try:

        protein_comments = uniprot_docs['comments']
        
        for i in range(len(protein_comments)):
          curPos = protein_comments[i]
          
          if curPos["commentType"] == "TISSUE SPECIFICITY":
            for j in curPos['texts']:
              uniprot_tissue += f"{marker}: " + j['value']
              
      except Exception as e:
        pass
    
  gene_cell_summary_context += cell_marker_text + "\n" + proteinAtlas_tissue + "\n" + uniprot_tissue + "\n"

  return gene_cell_summary_context
  
  
def scene_pathway(marker, genes_ids, db5, species):
  gene_pathway_context = "[Pathway related to Gene]\n"
  llm = OpenAI()

  if species == "arabidopsis":
    gene_pathway_context += f"{marker}:\n"
    query=f"(gene:{marker}) AND (organism_id:3702)"
    url=f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names,go_p"
    res = requests.get(url, headers={"Accept":"application/json"})
    data = res.json()
    prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene. You should extract biological process contents only and summarize it. Biological process is organized as the specific form (ex. value: P:cell division)
    
    [Gene]
    {marker}
    
    [Context]
    {data}"""
    
    messages = [{"role": "user", "content": prompt}]
    
    completion = llm.chat.completions.create(
        model = "gpt-4o",
        messages = messages,
        temperature = 0,
      )
    
    pathway_text = completion.choices[0].message.content
    gene_pathway_context += pathway_text + "\n"
      
  else:
    gene_pathway_context += f"{marker}:\n"
    mygene_pathway_text = ""
    try:
      gene_id = genes_ids[marker]
    
      fields = "pathway"
      mygene_url = f"https://mygene.info/v3/gene/{gene_id}?fields={fields}"
      mygene_response = requests.get(mygene_url)
      pathway_text = mygene_response.text
      
    
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene's pathway. You should summarize it.

    [Gene]
    {marker}
    
    [Context]
    {pathway_text}
    
      """
      messages = [{"role": "user", "content": prompt}]
    
      completion = llm.chat.completions.create(
        model = "gpt-4o",
        messages = messages,
        temperature = 0,
      )
    
      mygene_pathway_text = completion.choices[0].message.content
    
    except Exception as e:
      pass
    
    proteinAtlas_bp = ""
    if species == "human":
      proteinAtlas_bp = get_info_proteinAtlas(marker, db5, "Biological process")
    
    
    gene_pathway_context += mygene_pathway_text + "\n" + proteinAtlas_bp + "\n" + "\n"

  return gene_pathway_context

def scene_protein_location(marker, genes_ids, species):
  gene_protein_loc_context = "[Gene's Protein Location]\n"
  llm = OpenAI()

  gene_protein_loc_context += f"{marker}:\n"
  mygene_cc_text = ""
    
  try:
      gene_id = genes_ids[marker]
    
      fields = "go"
      mygene_url = f"https://mygene.info/v3/gene/{gene_id}?fields={fields}"
      mygene_response = requests.get(mygene_url)
      go_text = mygene_response.text
    
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene's biological process(BP), molecular function(MF), and cellular component(CC). You should extract CC contents only and summarize it.
    
    [Gene]
    {marker}
    
    [Context]
    {go_text}
    
      """
      messages = [{"role": "user", "content": prompt}]
    
      completion = llm.chat.completions.create(
        model = "gpt-4o",
        messages = messages,
        temperature = 0,
      )
    
      mygene_cc_text = completion.choices[0].message.content
    
  except Exception as e:
      pass
    
  uniprot_id = get_uniprot_id(marker, species)
  r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
  uniprot_docs = r.json()
  uniprot_cc = ""
    
  if uniprot_docs == None:
      pass
  else:
      try:

        protein_comments = uniprot_docs['comments']
        
        for i in range(len(protein_comments)):
          curPos = protein_comments[i]
          
          if curPos["commentType"] == "SUBCELLULAR LOCATION":
            for j in curPos['subcellularLocations']:
              uniprot_cc += f"Proteins encoded by {marker} are located on " + j['location']['value']
              
      except Exception as e:
        pass
    
    
    
  gene_protein_loc_context += mygene_cc_text + "\n" + uniprot_cc + "\n"
  
  return gene_protein_loc_context
  

def convert_id_to_name(ensemblID):
  mg = mygene.MyGeneInfo()
  buf_err = io.StringIO()
  with contextlib.redirect_stderr(buf_err):
    res = mg.query(ensemblID, scopes='ensembl.gene', fields='symbol')
  if len(res['hits']) != 0:
    try:
      return res['hits'][0]['symbol']
    except Exception as e:
      return ""
  else:
    return ""
    
    

def convert_name_to_id(name):
  mg = mygene.MyGeneInfo()
  buf_err = io.StringIO()
  with contextlib.redirect_stderr(buf_err):
    res = mg.query(name, scopes="symbol", fields="ensembl.gene", species="human")
  if len(res['hits']) != 0:
    try:
      return res['hits'][0]['ensembl']
    except Exception as e:
      return ""
  else:
    return ""
    
    

def p_p_interaction_info(marker, species, db6):
  ppi_context = "[Gene's Protein-Protein Interaction]\n"

  
  tax_ids = {"human":9606, "mouse":10090, "arabidopsis":3702}
  
  ppi_context += f"{marker}:\n"
    
  if species == "arabidopsis":
      pass
  else:
      marker_ids = convert_name_to_id(marker)
      if marker_ids != "":
        if type(marker_ids) == dict:
          marker_ids = [marker_ids]
        docs = []
        for marker_id in marker_ids:
          mask = (db6["ensembl_gene_id_1"] == marker_id['gene']) | (db6["ensembl_gene_id_2"] == marker_id['gene'])
          filter_db = db6[mask]
          partner_info = f"{marker} can interact with : "
  
          for record1, record2 in zip(filter_db['ensembl_gene_id_1'], filter_db['ensembl_gene_id_2']):
            if record1 == marker_id['gene']:
              partner_info += convert_id_to_name(record2)
              partner_info += " / "
            else:
              partner_info += convert_id_to_name(record1)
  
        ppi_context += partner_info + "\n"
    
  params = {
    "identifiers": marker,
    "species": tax_ids[species.lower()],
    "limit": 10,
    "required_score": 700,
    "format": "tsv"
  }
    
  res = requests.get("https://string-db.org/api/tsv/network", params=params)
    
  llm = OpenAI()
  prompt = f"""You are an intelligent summarizer that summarizes the network of genes. You will be given table representing the interaction, so you should organize the network of them with natural sentences. You should only represent the gene as symbol. 
    
    [Interaction Tables]
    {res.text}
    
    """
    
  messages = [{"role": "user", "content": prompt}]
  completion = llm.chat.completions.create(
        model = "gpt-4o",
         messages = messages,
         temperature = 0,
    )
     
  ppi_context += completion.choices[0].message.content + "\n"
  
  return ppi_context
