## cv_coverage
Author: aram.amalietti@gmail.com


**Dependencies** (these are the versions the script was developped with, newer versions should work, but if they don't, please use these versions):
```
python=3.7.7  
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

  `xl_in` *list of BED files containing landmarks around which to display the motifs for example*  
          `path_to_file1,path_to_file2`. *It can be a single file or a list of files, the plots will be adjusted accordingly;*  
  `motifs` *group of motifs to be displayed around the input landmarks, for example AAAA,CCCC,GGGG;*  
  `regions` *group of regions to be considered for analysis. Available regions are genome (intron, CDS, UTR3, UTR5, ncRNA, intergenic), whole_gene (intron, CDS, UTR3, UTR5), intergenic, intron, ncRNA, other_exon (UTR5, CDS), UTR3 and UTR5. Single or multiple regions can be analysed and the plots will adjust accodringly. For example intron,ncRNA;*  
  `kmer_len` *the length of analysed motifs (in bases);*  
  `fasta` *is the path to the genome in fasta format;*  
  `fai` *is the path to the genome index file;*  
  `regions_file` *if a segmentation file in GTF format, usualy obtained with iCount segment but can also come from   other sources;*  
  `smoothing` *is the size of smoothing window, usually 12;*  
  `percentile` *defines the percentile for thresholding genomic landmarks by score, if None no thresholding is used and all landmarks will be used for the analisys;*  
  `window` *flanking distance around the landmarks that is analysed and displayed;*  
  `use_scores` *If True score will be used for calculation else all landmarks will have equal weigths of 1;*  
  `n_cores` *is the number of threads used in the process;*   
  `chunk_size` *is the number of rows per thread (10000 is the usual value);*  
  `cap` *is the max value for any landmarks score. It is only used if use_scores is set to True;*

  


