# from: https://gist.github.com/arq5x/3719100
set -ex
echo "bigBedToBed must be on \$PATH"
echo "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"

cell_types=(gm12878 h1hesc helas3 hepg2 huvec k562) 
dir=data/combined
mkdir -p $dir; cd $dir;


files=""
for ct in "${cell_types[@]}"; do
    echo $ct;
    bb="http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/awgHub/byDataType/segmentations/jan2011/${ct}.combined.bb"
    wget --quiet $bb;
    F=$(basename $bb .bb)
    bigBedToBed $(basename $bb) stdout | cut -f 1-4 | gzip -c > ${F}.bedg.gz
    files="$files ${F}.bedg.gz"
done
 
bedtools unionbedg -header -names "${cell_types[@]}" -i $files > master.chromhmm.bedg
