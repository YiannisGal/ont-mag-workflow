set -uo pipefail
echo "CMD [$(date -Is)] seqkit grep -f ids.txt ${WS2}/00_input/sediment_reads_clean.fastq.gz" >> ../COMMANDS.log
../../env/.pixi/envs/default/bin/seqkit grep -j 8 -f ids.txt ${WS2}/00_input/sediment_reads_clean.fastq.gz -o reads.fastq && ../../env/.pixi/envs/default/bin/seqkit fq2fa reads.fastq -o reads.fa
echo "EXIT [$(date -Is)] rc=$?" >> ../COMMANDS.log
printf '%s\t%s\n' "$(md5sum < ${WS2}/00_input/sediment_reads_clean.fastq.gz | cut -c1-32)" ${WS2}/00_input/sediment_reads_clean.fastq.gz > reads_md5.tsv
echo DONE > extract.done
