#!/bin/bash
parallel 'cat {} | grep -v ^# | sed "s/^chrM/MT/g" | sed "s/^chr//g" > GRCh37/bed/{/}' ::: hg19/bed/*.bed
