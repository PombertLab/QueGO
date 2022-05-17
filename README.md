# QueGO (Query Gene Ontologies)
The QueGO pipeline identifies potential homologs to proteins related to specific functionality. Protein sequences and/or 3D structures are automatically downloaded from the UniProt Knowledgebase ([UniProtKB](https://www.uniprot.org/)) for proteins having Gene Ontology (GO) term descriptions that contain a desired keyword. The downloaded protein data are searched against the provided protein data using sequence- ([DIAMOND](https://github.com/bbuchfink/diamond)) and/or 3D-based ([Foldseek](https://github.com/steineggerlab/foldseek) or [GESAMT]()) homology, producing a list of potential homologs.

### <b>Why use QueGO?</b>
With the ever increasing amount of genome data, identifying all proteins relevant to functions of interest can be time consuming, error prone, and extremely diffucult to reproduce between research groups. QueGO automates these searches, returning and downloading identical results given the same input parameters (if no changes to the source database), as well as removes the point-and-click download required by manual curation, increasing the speed of data acquisition.

### <b>What are Gene Ontologies?</b>
The [Gene Ontology Consortium](http://geneontology.org/) has worked to create a class-based system of internationally agreed upon terms used to describe the roles and functions that a protein may take part in. Each term is a member of a class, and expands the specificity of the functional description of the class it belongs to, with all terms being rooted by one of three base terms: biological process, cellular component, or molecular function.

## <b>References</b>
[Sensitive protein alignments at tree-of-life scale using DIAMOND](https://doi.org/10.1038/s41592-021-01101-x). <b>Buchfink, B., Reuter, K., & Drost, H.</b> <i>Nat Methods</i> 18, 366–368 (2021) DOI: 10.1038/s41592-021-01101-x

[Foldseek: fast and accurate protein structure search](https://doi.org/10.1101/2022.02.07.479398) <b>Michel van Kempen, Stephanie S. Kim, Charlotte Tumescheit, Milot Mirdita, Johannes Söding, Martin Steinegger</b> bioRxiv 2022.02.07.479398; DOI: 10.1101/2022.02.07.479398

[UniProt: the universal protein knowledgebase in 2021](https://doi.org/10.1093/nar/gkaa1100). <b>The UniProt Consortium</b>, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D480–D489, DOI: 10.1093/nar/gkaa1100