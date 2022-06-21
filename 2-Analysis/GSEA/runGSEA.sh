#!/bin/bash

#$ -cwd
#$ -o  logs/runGsea$JOB_ID_$TASK_ID.txt
#$ -j y
#$ -t 1:18
#$ -l h_data=16G,h_rt=01:00:00
#$ -V
#$ -M videlock@mail
#$ -m a




. /u/local/Modules/default/init/modules.sh

module load java/13.0.1

cd /u/project/videlock/videlock/GSEA/GSEA_Linux_4.0.3 && source  /u/project/videlock/videlock/QuantSeqApr2019/FullWorkflow/E_GSEA/runFiles/*_${SGE_TASK_ID}.sh