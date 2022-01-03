## cv_coverage
This code compares the coverages of a k-mer group around crosslink events (landmarks) across multiple CLIP datasets. The comparison is conducted within user-specified genomic regions. 

Author: aram.amalietti@gmail.com


**Dependencies** (these are the versions the script was developed with, pandas >= 1 introduced breaking changes, please use these versions):
```
python=3.7  
pandas=0.24.2  
numpy=1.19.2  
pybedtools=0.8.1  
matplotlib=3.3.2
seaborn=0.11.0
```
**Usage**:  
  ```
  python3 <path_to_script> <xl_in> <motifs> <regions> <kmer_len> <fasta> <fai> <regions_file> <smothing> <percentile> <window> <use_scores> <n_cores> <chunk_size> <cap>
  ```

  - **`xl_in`** *list of BED files containing landmarks around which to display the motifs, given as* `path_to_file1,path_to_file2,path_to_file3 ...`. *It can be a single file or a list of files, the plots will be adjusted accordingly;*  
  - **`motifs`** *group of motifs to be displayed around the input landmarks, for example `AAAA,CCCC,GGGG ...`;*  
  - **`regions`** *group of regions to be considered for analysis. Available regions are `genome` (intron, CDS, UTR3, UTR5, ncRNA, intergenic), `whole_gene` (intron, CDS, UTR3, UTR5), `intergenic`, `intron`, `ncRNA`, `other_exon` (UTR5, CDS), `UTR3` and `UTR5`. Single or multiple regions can be analysed and the plots will adjust accordingly. For example, `intron,ncRNA`;*  
  - **`kmer_len`** *the length of analysed motifs (in bases);*  
  - **`fasta`** *is the path to the genome in fasta format;*  
  - **`fai`** *is the path to the genome index file;*  
  - **`regions_file`** *is a segmentation file in GTF format, usually obtained with iCount segment but can also come from other sources;*  
  - **`smoothing`** *is the size of the smoothing window, usually 12;*  
  - **`percentile`** *defines the percentile for thresholding genomic landmarks by score, if percentile is set to None no thresholding is used and all landmarks will be used for the analisys;*  
  - **`window`** *flanking distance around the landmarks that is analysed and displayed;*  
  - **`use_scores`** *If True, each landmark is weighted by its score, else all landmarks will have equal weigths of 1;*  
  - **`n_cores`** *is the number of threads used in the process (4 is the usual value);*   
  - **`chunk_size`** *is the number of rows per thread (10000 is the usual value);*  
  - **`cap`** *is the max value for any landmarks score. It is only used if use_scores is set to True (20 is the recommended value);*

 **Common issues**

 The script needs writing permission in the staging directory to make `results` directory and environment variable `TMPDIR` for temporary files. If you get `KeyError: 'TMPDIR'` a solution would be to type `export TMPDIR=<path_to_folder>` in terminal where you want to run the script.
 
 **Outputs**
 - A pdf file with graphs showing % k-mer group coverage around landmarks for each analyzed CLIP dataset;
 - A tsv file with % coverage values (raw values used for plotting);
 - A text file that saves the number of analyzed landmarks in each sample;
 - A tsv file with saved run parameters.
