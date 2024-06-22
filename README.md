<div style="display: flex; align-items: center; max-width: 800px; margin: 0 auto;">
  <div style="width: 30%;">
  <img src="./ribxz_logo_black.png" alt="Riboxyz Logo" style="width:450px; height:400px; padding: 10px;">
  </div>
  <div style="margin-left: 10px; width: 40%;">
    <p><code>riboxyz</code> (available at https://ribosome.xyz) is a package and a database application that provides organized access to ribosome structures, with several tools for visualisation and study. The database is up-to-date with the Protein Data Bank (PDB), provides a standardized nomenclature to ribosomal components:</p>
    <ul>
      <li>cytosolic and mitochondrial proteins</li>
      <li>cytosolic and mitochondrial rRNA</li>
      <li>tRNA</li>
      <li>elongation, initiation, termination factors (archaeal, bacterial and eukaryotic)</li>
    </ul>
    <p>The provided datatypes allow for seamless comparison and programming against ribosomal components across all the available structures. In addition, the application has several specialized visualization tools, including the identification and prediction of ligand binding sites, and 3D superimposition of subchains.</p>
    <p>The accompanying publication can be found <a href="#">here</a>.</p>
  </div>
</div>



<!-- <p align="center">
<img src="./logo.png" height="400" width="450" >
</p>


## Overview

`riboxyz` (available at https://ribosome.xyz) is a package and a database application that provides organized access to ribosome structures, with several tools for visualisation and study. The database is up-to-date with the Protein Data Bank (PDB), provides a [standardized nomenclature](https://github.com/rtviii/riboxyz/blob/master/ribctl/lib/ribosome_types/types_ribosome.py) to ribosomal components:

- cytosolic and mitochondiral rproteins
- cytosolic and mitochondrial rRNA
- tRNA
- elongation, initiation, termination factors  (archaeal, bacterial and eukaryotic)

The provided datatypes allow for seamless comparison and programming against ribosomal components across all the available structures. In addition, the application has several specialized visualization tools, including the identification and prediction of ligand binding sites, and 3D superimposition of subchains.


[ The accompanying publication can be found here ](https://academic.oup.com/nar/article/51/D1/D509/6777803). -->


------------------------------------------------------------------------------------------

### Structure

The application assumes three main components:

- `ribctl`: a package which provides the main types for annotation and parsing, functions for computing things from the mmcif structures.
- `neo4j-adapter`: a Neo4j database instance + `neo4j-adapter` package for managing the database (ingesting, querying, etc.)
- `api`: a Django-based restful API for querying the database and feeding the clients (frontend included)

### Installation

There is a number of config files and variables, some relevant to either of the packages above and some shared.

`RIBETL_DATA`
`NEO4J_MOUNTPATH`
`NEO4J_URI`
`NEO4J_USER`
`NEO4J_PASSWORD`
`NEO4J_CURRENTDB`
