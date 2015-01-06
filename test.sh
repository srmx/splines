#!/bin/bash

set -e

s="echo -n "

for i in `seq 1 1000`;
do
	a=$((RANDOM % 3))
	b=$RANDOM
	r="$a.$b"
	s="$s $r"
done

eval "$s | ./test"
