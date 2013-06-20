for g in hg18 hg19 mm8 mm9; do
    #mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    #            "select chrom, size from $g.chromInfo"  > data/$g.genome
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
                "select chrom, chromStart, chromEnd, name from $g.cpgIslandExt"  > data/cpgIsland.$g.bed
done

exit;
for g in hg19 mm9; do
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -NAe \
                "select chrom, chromStart, chromEnd, type from $g.gap where type = 'centromere'"  > data/$g.centromere.bed
done

for db in hg18 hg19 mm8 mm9; do 
    python data/get-gene-regions.py $db | sort -k1,1 -k2,2n > data/${db}.gene-features.bed;
done


wget -O - ftp://encodeftp.cse.ucsc.edu/pipeline/hg19/wgEncodeRegTfbsClustered/release2/wgEncodeRegTfbsClusteredV2.bed.gz \
    | zcat - \
    | awk '{split($4, a, "_"); print $1"\t"$2"\t"$3"\t"a[1] }' \
    | gzip -c > data/wgEncodeRegTfbsClusteredV2.bed.gz
