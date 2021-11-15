riboPLATE-seq: code required for analysis and interpretation of riboPLATE-seq data

Organization of repository: 
scripts/ -- contains main Jupyter notebook (riboPLATEseq.ipynb) and the scripts it calls
data/ -- contains input count matrices and barcode files
workingdir/ -- temporary directory for intermediate output, e.g. DESeq2 results and normalized counts
gsea/ -- temporary directory for GSEA output specifically
ReferenceInfo/ -- contains gene lists of interest, should hold GTF file as well (GTF too large to distribute on Github)
outdir/ -- Final output files, i.e. figures

Notebook requres GRCh38 RefSeq GTF annotations, available from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39

Notebook also requires GSEA command-line installation -- set 'gspath' variable in the notebook to the directory containing GSEA
