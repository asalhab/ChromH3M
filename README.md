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
you have to change the TMPDIR, give the path of ChromHMM.jar and Rscript in [ChromH3M.sh](https://github.com/asalhab/ChromH3M/blob/master/ChromH3M.sh#L58) script
