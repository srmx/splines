#!/bin/bash

set -e
# echo 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 \
# | ./linalg

s="echo -n "

for i in `seq 1 1000`;
do
	a=$((RANDOM % 3))
	b=$RANDOM
	r="$a.$b"
	s="$s $r"
done

eval "$s | ./linalg"
