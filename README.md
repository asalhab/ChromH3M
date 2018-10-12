# ChromH3M
# Description:
ChromH3M (abbreviation for ChromH**MM-m**eta segmentation) is an easy and straight-forward workflow to measure the similarity of PMDs/LMRs/UMRs produced by [MethylSeekR](https://bioconductor.org/packages/release/bioc/html/MethylSeekR.html) across many samples.
It takes segment files produced by MethylSeekR as input and binarizes the binned genome accordingly across all samples. [ChromHMM](http://compbio.mit.edu/ChromHMM/) is applied to this binarized signal with different number of states defined by the user. The emission probabilities are then hierarchically clustered and annotations are added to the heatmap based on a samplesheet provided by the user.

# Usage:
 ``bash ChromH3M.sh  -i dir     -g genome      -n name      -o output folder    -a min      -b max      -s sample sheet``

 - Mandatory:
   - -i bed files directory (give the full directory!)
   - -g genome length (Shortcuts: hg19 or mm10)
   - -n output name
   - -o output folder name (full directory!)
   - -s sample annotations sheet (full directory!)
   - -a minimum number of ChromHMM states
   - -b maximum number of ChromHMM states

**Caution:**
you have to change/provide the following:
 - temporary directory where the calculations are done [TMPDIR](https://github.com/asalhab/ChromH3M/blob/master/ChromH3M.sh#L61)
 - give the path of ChromHMM.jar and Rscript in [ChromH3M.sh](https://github.com/asalhab/ChromH3M/blob/master/ChromH3M.sh#L58)
 - give the genome length and the corresponding gaps file paths in [meth\_avg.sh](https://github.com/asalhab/ChromH3M/blob/master/meth_avg.sh#L49)
 
# Citation:
please cite the following paper if you used ChromH3M in your analysis:

Salhab A, Nordström K, Gasparoni G, Kattler K, Ebert P, Ramirez F, Arrigoni L, Müller F, Polansky JK, Cadenas C, et al. A comprehensive analysis of 195 DNA methylomes reveals shared and cell-specific features of partially methylated domains. Genome Biology. 2018; 19(1):150.

```Tex
@article{salhab2018comprehensive,
  title={A comprehensive analysis of 195 DNA methylomes reveals shared and cell-specific features of partially methylated domains},
  author={Salhab, Abdulrahman and Nordstr{\"o}m, Karl and Gasparoni, Gilles and Kattler, Kathrin and Ebert, Peter and Ramirez, Fidel and Arrigoni, Laura and M{\"u}ller, Fabian and Polansky, Julia K and Cadenas, Cristina and others},
  journal={Genome Biology},
  volume={19},
  number={1},
  pages={150},
  year={2018},
  publisher={BioMed Central},
  doi={10.1186/s13059-018-1510-5},
  url={https://doi.org/10.1186/s13059-018-1510-5}
}
```

