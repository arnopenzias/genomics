#!/bin/bash

#################################### Script for genome assembling of  Salmonella clinical isolates ##############################################

basedir=$1

if [ ! ${basedir} ]; then
	basedir="."
fi

if [ ! -d ${basedir} ]; then
	echo "Wrooong!!"
fi

echo "Using basedir: ${basedir}/"

current_dir=`pwd`

list1=( #list of directories names referring to different bacterial genomes sequenced # each raw fastq must be inside raw_data dir
)
list2=( "MeDuSa" "a5_pipeline" "a5_trimmed" )


for samp in ${list1[@]}
do
	for folder1 in ${list2[@]}
	do
		mkdir -p ${basedir}/${samp}/${folder1}
	done
	cd ${samp}/a5_pipeline
	a5_pipeline.pl ../raw_data/${samp}.R1.fastq ../raw_data/${samp}.R2.fastq --end=5 ${samp}_a5_scaff
	mv ${samp}_a5_scaff.ec.fastq.gz ../a5_trimmed/
	cd ../a5_trimmed
	gunzip ${samp}_a5_scaff.ec.fastq.gz
	deinterleave_fastq.sh < ${samp}_a5_scaff.ec.fastq ${samp}_a5.1.fastq ${samp}_a5.2.fastq
	rm ${samp}_a5_scaff.ec.fastq
	cd ${current_dir}
	java -jar ~/programs/medusa/medusa.jar -f ref_data/S.enterica_fna_dir/ -i ${samp}/a5_pipeline/${samp}_a5_scaff.final.scaffolds.fasta -o ${samp}/MeDuSa/scaffold0.fa -v
	it=0
	while [ ${it} -lt 5 ]
	do
		java -jar ~/programs/medusa/medusa.jar -f ref_data/S.enterica_fna_dir/ -i ${samp}/MeDuSa/scaffold${it2[${it}]}.fa -o ${samp}/MeDuSa/scaffold${it3[${it}]}.fa -v
		(( it++ ))
	done
	cd ${samp}/MeDuSa
	if [ -f "scaffold2.fa" ]; then
		cp scaffold2.fa scaffold.fasta
	else
	        if [ -f "scaffold1.fa" ]; then
			cp scaffold1.fa scaffold.fasta
		else
			if [ -f "scaffold0.fa" ]; then
				cp scaffold0.fa scaffold.fasta
			else
				echo "there is no MeDuSa scaffold file!!!!!"
			fi
		fi
	fi
	cd ${current_dir}
	mkdir -p ${basedir}/${samp}/Gapfiller/Gapfiller1
	cp ./${samp}/MeDuSa/scaffold.fasta ./${samp}/Gapfiller/Gapfiller1/Gapfiller1.gapfilled.final.fa
	cd ./${samp}/Gapfiller
	echo "Lib1 bwasw /path/to/${samp}/a5_trimmed/${samp}_a5.1.fastq /path/to/${samp}/a5_trimmed/${samp}_a5.2.fastq 800 0.5 FR" > library.txt
	perl ~/programs/GapFiller_v1-10_linux-x86_64/GapFiller.pl -l library.txt -s /path/to/${samp}/Gapfiller/Gapfiller1/Gapfiller1.gapfilled.final.fa -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 4 -i 20 -b Gapfiller20
	cd ${current_dir}
	mkdir -p ${basedir}/${samp}/final_genome
	cp ${samp}/Gapfiller/Gapfiller20/Gapfiller20.gapfilled.final.fa ${samp}/final_genome/${samp}.whole.final.fa
	prokka --outdir /path/to/${samp}/prokka_a5 --mincontiglen 200 --cpus 4 --prefix ${samp} --addgenes --addmrna --locustag ST --compliant --rfam --genus Salmonella --usegenus --kingdom Bacteria --gram neg /media/patrick/PASIQUIBAC/Bioinfo/genomes/${samp}/final_genome/${samp}.whole.final.fa
	cd ${current_dir}
	pullseq -i ${samp}/final_genome/${samp}.whole.final.fa -m 3000 >> ${samp}/final_genome/${samp}.final.fa
	echo "########Stats from ${samp}#############" >> stats_all_ST313.txt
	assemblathon_stats.pl ${samp}/final_genome/${samp}.final.fa >> stats_all_ST313.txt
done

###  genomic island identification #############

for gi in ${list1[@]}
do
        cd ${gi}/prokka_a5/
        mkdir input
        cp ${gi}.gbk input/
        ln -s ~/programs/vm_programs/SW_Sniffer/lib/
        ln -s ~/programs/vm_programs/SW_Sniffer/bin/
        python2 ~/programs/vm_programs/SW_Sniffer/SeqWordSniffer.py
        cd ${current_dir}
done

			
echo "Finito!!"

exit
