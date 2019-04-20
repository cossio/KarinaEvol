#!/bin/bash
# author: cossio


export d=0
export kill=100
export mutfreq=0.1
export init=SI-eq.txt
export T=100.

for gpow in `seq -1 1 0.1`; do
	for lpow in `seq -4 0.5 1`; do
		for rep in `seq 1 100`; do
	
			export g=`perl -E "say 10**$gpow"`
			export lambda=`perl -E "say 10**$lpow"`
		
			outdir="g=$g-lambda=$lambda-rep=$rep"
			mkdir -p $outdir
		
			export out="$outdir/out.txt"
		
			sbatch --ntasks=1 --ntasks-per-core=1 --output="$out/slurm.out" --job-name=$outdir run.sh
			
		done
	done
done
