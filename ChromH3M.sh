#!/bin/bash
#$ -S /bin/bash
#$ -N ChromH3M
#$ -cwd
#$ -j y
#$ -V

Pipeline_name="ChromH3M.sh"

#### How to use this script
printHelp() {
	echo -e "$(tput bold)Description:$(tput sgr0)"
	echo -e "ChromH3M is an easy and straight-forward workflow to measure the similarity of PMDs/LMRs/UMRs across many samples."
	echo -e "It takes the segment files from MSRv* pipeline as input and binarizes the binned genome accordingly across all samples. ChromHMM is applied to this binarized signal with different number of states defined by the user. The emission probabilities are then hierarchically clustered and annotations are added to the heatmap based on a samplesheet provided by the user."
	echo -e ""
	echo -e "$(tput bold)Usage:$(tput sgr0)"
	echo -e " bash $Pipeline_name  $(tput bold)$(tput setaf 1)-i dir     -g genome      -n name      -o output folder    -a min      -b max      -s sample sheet$(tput sgr0)"
	echo -e ""
	echo -e " $(tput bold)Mandatory:$(tput sgr0)"
	echo -e "  -i bed files directory ($(tput bold)give the full directory!$(tput sgr0))"
	echo -e "  -g genome length (Shortcuts: hg19 or mm10)"
	echo -e "  -n output name"
	echo -e "  -o output folder name (full directory!)"
	echo -e "  -s sample annotations sheet ($(tput bold)full directory!$(tput sgr0))"
	echo -e "  -a minimum number of ChromHMM states"
	echo -e "  -b maximum number of ChromHMM states"
}

while getopts ":hi:g:n:o:s:a:b:" opt
do
   case $opt in
      h) printHelp; exit 0 ;;
      i) input="$OPTARG" ;;
      g) genome="$OPTARG" ;;
      n) name="$OPTARG" ;;
      o) outFolder="$OPTARG" ;;
      s) sample_annotation="$OPTARG" ;;
      a) min="$OPTARG" ;;
      b) max="$OPTARG" ;;
      *) printHelp; exit 1 ;;
   esac
done


#### check if all the necessary arguments were supplied
if [[ -z "$input" ]] || [[ -z "$genome" ]] || [[ -z "$name" ]] || [[ -z "$outFolder" ]] || [[ -z "$sample_annotation" ]] || [[ -z "$min" ]] || [[ -z "$max" ]] 
   then
   echo -e "\n$(tput bold)$(tput setaf 1)ERROR ($Pipeline_name): Must set all mandatory options$(tput sgr0)\n"
   printHelp
   exit 1
fi

### helper scripts
meth_avg="meth_avg.sh"
Heatmap="Heatmap.R"

### needed softwares
chromhmm="Path to ChromHMM.jar"
Rscript="Path to Rscript"

export TMPDIR=/tmp
TMPWD=`mktemp --tmpdir=$TMPDIR -d ChromH3M_Pipeline.XXXXXXXXXX`
if [ ! -d "$TMPWD" ]
   then

   echo "failed to create folder:"
   echo "$TMPWD"
   exit 1

fi

mkdir -p $outFolder

echo -e "`date` LOGG ($Pipeline_name): starting ......" 
echo -e "`date` LOGG ($Pipeline_name): making soft links for the bed files"
mkdir -p $TMPWD/commonPMDs/{bed,output} $TMPWD/commonLMRs/{bed,output} $TMPWD/commonUMRs/{bed,output}
ln -s ${input}/*gz $TMPWD/commonPMDs/bed/
ln -s ${input}/*gz $TMPWD/commonLMRs/bed/
ln -s ${input}/*gz $TMPWD/commonUMRs/bed/


for seg in PMD LMR UMR
do
	echo -e "`date` LOGG($Pipeline_name): running meth_avg script for ${seg}s"
	bash $meth_avg -i $TMPWD/common${seg}s/bed -o $TMPWD/common${seg}s/output -s $seg -n ${name} -g $genome

	case "$seg" in
			PMD) bin=1000
   		;;
   		LMR) bin=200
   		;;
   		UMR) bin=200
   		;;
	esac

	echo -e "`date` LOGG($Pipeline_name): preparing the binary files for $seg chromhmm"
	cd $TMPWD/common${seg}s
	mkdir ${seg}s_chromhmm
	cd ${seg}s_chromhmm
	awk '{print  > $1".bed"}' < ../output/merged.sampels.${seg}.bed
	mv chr.bed header.txt
	for i in * ; do cut -f 4- $i > tmp && mv tmp $i; done
	for i in *bed; do cat <(echo -e "mix\tchr${i%.*}") <(cat header.txt) <(cat $i) > tmp && mv tmp $i; done
	rm X.bed 
	rename s/.bed/_binary.txt/ *bed
	mkdir input_binarized
	mv *binary* input_binarized/

	echo -e "`date` LOGG($Pipeline_name): running chomhmm for ${seg}s "
	for i in `seq $min $((($max-$min)/($min-1))) $max`
	do
		echo -e "for ${i} states"
		java -mx32000M -jar $chromhmm LearnModel -b ${bin} input_binarized ${seg}s_mix_${i}St $i $genome
	done



	folder="$TMPWD/common"
   $Rscript --vanilla $Heatmap $folder $seg $min $max $sample_annotation
done

echo -e "`date` LOGG($Pipeline_name): copying to the output folder"
cp -r $TMPWD/common* $outFolder
rm -r $TMPWD
