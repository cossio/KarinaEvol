#!/bin/bash
# author: cossio


export lambda=0.01
export d=0
export mutfreq=0.1
export init=SI-eq.txt
export T=100.

for gpow in `seq -1 1 0.1`; do
	for kipow in `seq -4 0.5 1`; do
		for rep in `seq 1 100`; do
	
			export kill=`perl -E "say 10**$kipow"`
			export g=`perl -E "say 10**$gpow"`
			
			outdir="g=$g-ki=$kill-rep=$rep"
			mkdir -p $outdir
			
			export out="$outdir/out.txt"
			
			sbatch --ntasks=1 --ntasks-per-core=1 --output="$out/slurm.out" --job-name=$outdir run.sh

		done
	done
done
