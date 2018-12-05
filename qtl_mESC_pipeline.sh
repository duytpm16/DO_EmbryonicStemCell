#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00


module load R/3.5.1


for i in {1..8}
do

echo "#PBS -l nodes=1:ppn=8
#PBS -q batch
#PBS -l walltime=72:00:00

module load R/3.5.1

Rscript scan1_HPC.R do_mESC_proteomics_qtl_viewer.RData dataset.esc.proteins rankz 8 $i 1000" >> map_mESC_$i.sh
qsub map_mESC_$i.sh

done
