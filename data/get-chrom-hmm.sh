# taken from @arq5x
set -ex
cell_types=(Gm12878 H1hesc Hepg2 Hmec Hsmm Huvec K562 Nhek Nhlf)
dir=data/chromHMM

mkdir -p $dir; cd $dir;

files=""
for ct in "${cell_types[@]}"; do
    echo $ct;
    remote=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmm${ct}HMM.bed.gz
    F=$(basename $remote .gz)
    wget --quiet -O - $remote | zcat - | cut -f 1-4 | perl -pe 's/\d+_(.+)/$1/' > $F
    files="$files $F"
done
 
bedtools unionbedg -header -names "${cell_types[@]}" -i $files > master.chromhmm.bedg
