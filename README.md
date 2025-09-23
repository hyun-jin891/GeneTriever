<h1 align="center"> GeneTriever </h1> <br>
<p align="center">
  <a href="https://github.com/hyun-jin891/GeneTriever/tree/main">
    <img alt="GeneTriever" title="GeneTriever" src="https://i.imgur.com/CVkX8OY.png" width="450">
  </a>
</p>


<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
## Table of Contents

- [Introduction](#introduction)
- [Pipeline](#pipeline)
- [Data](#data)
- [How to use](#how-to-use)
- [Application](#application)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Introduction
<p align="center">
  <a href="https://github.com/hyun-jin891/GeneTriever/tree/main">
    <img alt="Target Species" title="Target Species" src="https://i.imgur.com/N7329gV.png" width="450">
  </a>
</p>

* It helps you to get information about genes of **Human, Mouse, and Arabidopsis**
* Gene Information
  * Gene Description
  * Gene Function
  * Cell Type related to Gene
  * Pathway
  * Subcellular Location
  * Protein-Protein Interaction

## Pipeline
* Using the official api of various websites and preloaded datasets downloaded from them, it extracts information about the target gene.
* It uses LLM for organizing the retrieved contents. (You should use your own OpenAI API key)

## API
* NCBI
* Uniprot
* MyGene.info
* STRING

## Data
You should check the file names and download them. Insert them into main.py

* Cell_marker_Seq.xlsx (http://www.bio-bigdata.center/CellMarker_download.html)
* PanglaoDB_markers_27_Mar_2020.tsv (https://panglaodb.se/markers.html?cell_type=%27choose%27)
* singleCellBase_dataset.txt (http://cloud.capitalbiotech.com/SingleCellBase/download.jsp)
* proteinatlas.tsv (https://www.proteinatlas.org/about/download)
* interaction_consensus.tsv (https://www.proteinatlas.org/humanproteome/structure+interaction/interaction/data#protein_interaction_data)
* cl.obo (https://github.com/obophenotype/cell-ontology/releases)
* cl-basic.obo (https://github.com/obophenotype/cell-ontology/releases)
* PCMDB_Download_20250909177169.csv (https://www.tobaccodb.org/pcmdb/download)
* arabidopsis_thaliana.marker_fd.csv (https://biobigdata.nju.edu.cn/scplantdb/marker)

## How to use
* You should get an OpenAI API key first and make .env file.
<pre>
<code>
.env

OPENAI_API_KEY=
</code>
</pre>





## Application
* It is used to design [CellTyper](https://github.com/SBBlaboratory/CellTyper), which is an LLM-based cell annotation tool.





