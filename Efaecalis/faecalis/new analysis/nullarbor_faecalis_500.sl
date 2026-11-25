#!/bin/bash -e

#SBATCH --account massey04106

#SBATCH --job-name=nullarbor # job name (shows up in the queue)

#SBATCH --output=nullarbor_long_faecalis.log

#SBATCH --error         %x_%j.err

#SBATCH --begin=now+6hours

#SBATCH --output        %x_%j.out

#SBATCH --time=16:30:00 # Walltime (HH:MM:SS)

#SBATCH --cpus-per-task=16 # number of cpus/threads, in our case 4

#SBATCH --mem=60GB            # Memory in MB 

module purge

module load nullarbor/2.0.20191013

nullarbor.pl --name nullarbor --mlst efaecalis --ref GCA_022869705.1_ASM2286970v1_genomic.fna  --trim ON --taxoner centrifuge --input faecalis_500.txt --outdir results_faecalis_long_500 --force

nice make all -j 4 -l 54 -C results_faecalis_long_500 2>&1 | tee -a results_faecalis_long_500/nullarbor.log

