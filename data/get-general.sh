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

