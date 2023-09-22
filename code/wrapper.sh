#!/bin/bash

#change to wdir
WDIR="/flash/anze/Work/Projects/example_dir"

Nmin=$1
Nmax=$2

for ((N=${Nmin};N<=${Nmax};N++)) ; do

#create the file for slurm
FILE="example_N"${N}".sh"
/bin/cat <<EOM >$FILE
#!/bin/bash

#SBATCH --time=23:59:59
#SBATCH --mem=4096

python3 ${WDIR}/sphere-assembly/example.py --number ${N} -S 250
EOM

#submit file to slurm
#create separate stdout and stderr files
sbatch -o ${WDIR}/log/${FILE%.sh}.stdout -e ${WDIR}/log/${FILE%.sh}.stderr ${WDIR}/${FILE}

#archive submitted script
mv $FILE ${WDIR}/slurm

done
