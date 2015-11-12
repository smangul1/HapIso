if [ $# -lt 3 ]
then
	echo "assume chr1"
	echo "[1] - list of genes with boundaries"
	echo "[2] - bam"
	echo "[3] - chr"
	echo "->bams"
	exit 1
fi

bam=$2
sample="${bam%.*}"
chr=${3}

echo $sample

while read line
do 
echo $line

gene=$(echo $line | awk '{print $2}')
g1=$(echo $line | awk '{print $3}')
g2=$(echo $line | awk '{print $4}')

echo "samtools view -bh $bam ${chr}:${g1}-${g2} >${sample}_${gene}.bam"

samtools view -bh $bam ${chr}:${g1}-${g2} >${sample}_${gene}.bam
samtools index ${sample}_${gene}.bam



done<$1

