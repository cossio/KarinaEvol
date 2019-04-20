#!/bin/bash
# author: cossio

# variables set outside
# g, d, kill, lambda, mutfreq, init, out, T


julia evol.jl --g $g --d $d --kill $kill --lambda $lambda --mutfreq $mutfreq --si $init --out $out --T $T

