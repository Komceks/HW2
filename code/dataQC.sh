#!/bin/sh
dir=/home/bioinformatikai/HW2
threads=6
# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR204044/ERR204044 -P ~/HW2/references/
# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR15131330/SRR15131330 -P ~/HW2/references/
# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR18214264/SRR18214264 -P ~/HW2/references/
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/832/545/GCF_022832545.1_ASM2283254v1/GCF_022832545.1_ASM2283254v1_genomic.gff.gz -P ~/HW2/references/
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/832/545/GCF_022832545.1_ASM2283254v1/GCF_022832545.1_ASM2283254v1_genomic.fna.gz -P ~/HW2/references/
# fastq-dump --split-files ERR204044
# fastq-dump --split-files SRR15131330
# fastq-dump --split-files SRR18214264

#Loop
for i in ${dir}/inputs/*_1.fastq
do
#Sample nr.1
R1=${i};

#Base name of samples
base=$(basename ${i} _1.fastq);

#Sample nr.2
R2="${dir}/inputs/"$(basename ${R1} _1.fastq)"_2.fastq";

#Trimmed sample nr.1
OUT1="${dir}/outputs/trimmed/"$(basename ${i} _1.fastq)"_1_val_1.fq";

#Trimmed sample nr.2
OUT2="${dir}/outputs/trimmed/"$(basename ${i} _1.fastq)"_2_val_2.fq";
# DATA QC ---------------------------------------------------------------------
#fastqc analysis on samples
# ERR204044 kokybė yra nebloga. SRR15131330 kokybė labai gera.
# SRR18214264 kokybė nelabai patiko, šis mėginys turėjo mažiausiai duomenų, kokybė žemiausia
# iš visų trijų mėginių. Matosi didelis kokybės kritimas ties dešinu galu, turi daug skirtingų ilgių sekų.
fastqc -o "${dir}/outputs/raw/" -t ${threads} ${R1} ${R2};

# #trimming on samples
trim_galore -j ${threads} -q 30 -o "${dir}/outputs/trimmed/" --paired ${R1} ${R2};

# SRR18214264 kokybės kritimas ties dešinu galu sumažėjo. ERR204044 kokybės kritimas ties dešinu
# galu sumažėjo, ties kairiu išliko tokia pati. SRR15131330 kokybė išliko tokia pati -gera.

# fastqc on trimmed samples
fastqc -o "${dir}/outputs/trimmed/" -t ${threads} ${OUT1} ${OUT2};
# GENOME ASSEMBLY ---------------------------------------------------------------------
# 1
spades.py -t ${threads} -1 ${OUT1} -2 ${OUT2} -o ${dir}/outputs/"${base}_assembly"
# 2
velveth ${dir}/outputs/"${base}_velvet" 31 -shortPaired -fastq ${OUT1} ${OUT2};
velvetg ${dir}/outputs/"${base}_velvet" -cov_cutoff auto;
# 3, 4
# Dariau atskirai, nes nerodė rezultatų viename faile. Tikriausiai dėl to, nes per maži kontigai po velvet, default >500. 
# Negaliu palyginti rezultatų, tačiau žiūrint po spade ir lyginant su ref. genomu: SRR15131330 susilygino apie 75 proc. - vidutiniškai,
# SRR18214264 98 proc. -labai geras rezultatas, ERR204044 98 proc. -labai geras rezultatas.
python3 ~/quast-5.2.0/quast.py -t ${threads} -o ${dir}/outputs/"${base}_quast_assembly" -r ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna ${dir}/outputs/"${base}_assembly"/contigs.fasta;
#Labai lievai gavosi, nes maži kontigai. Kodėl taip - nespėjau išsiaiškinti, reikėjo labiau pažaisti su velvet.
python3 ~/quast-5.2.0/quast.py -m 30 -t ${threads} -o ${dir}/outputs/"${base}_quast_velvet" -r ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna ${dir}/outputs/"${base}_velvet"/contigs.fa;
# 5
PATH=$PATH:/home/bioinformatikai/RagTag
ragtag.py scaffold -o ${dir}/outputs/"${base}_ragtag_assembly" --aligner ~/minimap2/minimap2 -t 6 -u ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna ${dir}/outputs/"${base}_assembly"/contigs.fasta
ragtag.py scaffold -o ${dir}/outputs/"${base}_ragtag_velvet" --aligner ~/minimap2/minimap2 -t 6 -u ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna ${dir}/outputs/"${base}_velvet"/contigs.fa
# 6
# Pasirinkau visus tris assembly po spade.py, nes galėjau įvertinti panašumą su ref. genomu. Taip pat rezultatai nėra blogi, mažiausias panašumas 75 proc.
# 7
# SRR15131330 fraction - 99.9 proc.(velvet 99.77), SRR18214264 fraction - 99.82 proc(velvet 98.77)., ERR204044 fraction - 99.75 proc.(velvet 99.67)
#  (iš ${base}_ragtag_assembly_flagstat.txt failų). Labai geri rezultatai, gauta tai ko tikėtasi, tačiau velvet šiek tiek mažesni.
# ERR204044 coverage - 4103890(velvet 4032582), SRR15131330 coverage - 3916024(velvet 15262754), SRR18214264 coverage - 4086628(velvet 11634462)
# (iš awk '{sum += $3} END {print sum}' "genome_coverage_${base}_ragtag_assembly.txt";). Gauti panašūs rezultatai iš assembly mėginių, tačiau iš velvet kažkas kitko, labai skiriasi nuo assembly rez.
# ir tarpusavy nelabai sutampa, taip galėjo būti dėl to, jog nepavyko tinkamai atlikti velvet.
bwa index ${dir}/outputs/"${base}_ragtag_assembly"/ragtag.scaffold.fasta
bwa index ${dir}/outputs/"${base}_ragtag_velvet"/ragtag.scaffold.fasta
bwa mem -t ${threads} ${dir}/outputs/"${base}_ragtag_assembly"/ragtag.scaffold.fasta ${OUT1} ${OUT2} -o ${dir}/outputs/"${base}_ragtag_assembly.sam"
bwa mem -t ${threads} ${dir}/outputs/"${base}_ragtag_velvet"/ragtag.scaffold.fasta ${OUT1} ${OUT2} -o ${dir}/outputs/"${base}_ragtag_velvet.sam"

samtools view -bS ${dir}/outputs/"${base}_ragtag_velvet.sam" -@ 6 | samtools sort -@ 6 - -o ${dir}/outputs/"${base}_ragtag_velvet.bam"
samtools view -bS ${dir}/outputs/"${base}_ragtag_assembly.sam" -@ 6 | samtools sort -@ 6 - -o ${dir}/outputs/"${base}_ragtag_assembly.bam"

samtools index -@ 6 ${dir}/outputs/"${base}_ragtag_assembly.bam"
samtools index -@ 6 ${dir}/outputs/"${base}_ragtag_velvet.bam"

samtools flagstat -@ 6 ${dir}/outputs/"${base}_ragtag_velvet.bam" > ${dir}/outputs/"${base}_ragtag_velvet_flagstat.txt"
samtools flagstat -@ 6 ${dir}/outputs/"${base}_ragtag_assembly.bam" > ${dir}/outputs/"${base}_ragtag_assembly_flagstat.txt"

samtools faidx ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna

bedtools bamtobed -i ${dir}/outputs/"${base}_ragtag_assembly.bam" | sort -k1,1 -k2,2n > ${dir}/outputs/"${base}_ragtag_assembly.bed"
bedtools bamtobed -i ${dir}/outputs/"${base}_ragtag_velvet.bam" | sort -k1,1 -k2,2n > ${dir}/outputs/"${base}_ragtag_velvet.bed"

bedtools genomecov -ibam ${dir}/outputs/"${base}_ragtag_velvet.bam" -g ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna.fai > "genome_coverage_${base}_ragtag_velvet.txt"
bedtools genomecov -ibam ${dir}/outputs/"${base}_ragtag_assembly.bam" -g ${dir}/references/GCF_022832545.1_ASM2283254v1_genomic.fna.fai > "genome_coverage_${base}_ragtag_assembly.txt"
# GENOME ANALYSIS AND ANNOTATION ---------------------------------------------------------------------
# 1 Gepard: Tiesiausia istrižainė matoma tarp ERR204044 ir SRR18214264, tai reiškia, kad jie labiausiai panašūs.
# Tada matome šiek tiek mažesni panašumą su ERR204044 ir SRR15131330.
# Ir mažiausias panašumas yra tarp SRR15131330 ir SRR18214264. Tai reiškia, kad SRR15131330 ir SRR18214264 labiausiai skiriasi.
# 2 Busco: Rezultatas neblogas, sekos buvo atpažintos kaip lactobacillales bakterijų eilės.
# SRR264 99 proc. "subrendusių" sekų, 0.5 proc. fragmentuotų, 0.5 proc. trūkstamų.
# SRR1330 98.5 proc. "subrendusių" sekų, 0.5 proc. fragmentuotų, 1 proc. trūkstamų.
# ERR 99 proc. "subrendusių" sekų, 0.5 proc. fragmentuotų, 0.5 proc. trūkstamų.
# 6 Palyginimas iš blast, rast, geneMarkS2:
# geneMarkS2: 330: 2559 genai buvo atrasti. 264: 2327 genai buvo atrasti. ERR: 2301 genas buvo atrastas.
# RAST: ERR(RAST ID 856): 2428 genai buvo atrasti; 264(RAST ID 855): 2397 genai buvo atrasti; 330(RAST ID 854): 2645 genai buvo atrasti.
# BLAST: ERR: 2397 genai, 330: 4600 genai, 264: 2797 genai.(iš contig_prot.fasta.blastn failų)
# Galime pastebėti, jog BLAST su 1e-75 e-value treshold parodo daugiausiai genų, o geneMarkS2 ir RAST parodo mažiau. Tačiu rezultatai skiriasi neperdaugiausiai, išskyrus 330 po blast.
makeblastdb -in ${dir}/outputs/"${base}_assembly"/contigs.fasta -dbtype nucl

blastn -query ${dir}/inputs/sequenceGenes.txt -db ${dir}/outputs/"${base}_assembly"/contigs.fasta > ${dir}/outputs/"${base}_assembly"/blastgene_info.txt;
tblastn -query ${dir}/inputs/sequenceProtein.txt -db ${dir}/outputs/"${base}_assembly"/contigs.fasta > ${dir}/outputs/"${base}_assembly"/blastprot_info.txt;

blastn -query ${dir}/inputs/sequenceGenes.txt -db ${dir}/outputs/"${base}_assembly"/contigs.fasta -out ${dir}/outputs/"${base}_assembly"/contig_gene.fasta.blastn -outfmt 6 -evalue "1e-75";
tblastn -query ${dir}/inputs/sequenceProtein.txt -db ${dir}/outputs/"${base}_assembly"/contigs.fasta -out ${dir}/outputs/"${base}_assembly"/contig_prot.fasta.blastn -outfmt 6 -evalue "1e-75";
# 10 filogenetiniai medžiai: Gavosi vienodi medžiai, labai nudžiugino tokie rezultatai, po sunkaus ir ilgo darbo. Imiau Staphilococus Aureus kaip outgroup, medis parodė,
# jog šio organizmo sekos yra labiausiai nutolusios nuo kitų trijų mėginių, taip ir turėjo būti. Nuostabu.
# 11 rezultatai: Galime pamatyti iš gepard ir filogenetinių medžių, jog ERR204044 ir SRR18214264 yra labiausiai panašūs. Rezultatai mane kiek nustebino, nes atrodė,
# kad iš 264 mėginio nieko gero nebus. Tačiau matome, jog su daugiausiai duomenų turinčiu mėginiu rezultatas buvo geras, bet gal dėl to ir geras, nes mažai duomenų 264 mėginyje.
done;
dir=/home/bioinformatikai/HW2

#Create MultiQC plots for your raw and processed data
multiqc -p ${dir}/outputs/raw/
multiqc -p ${dir}/outputs/trimmed/
