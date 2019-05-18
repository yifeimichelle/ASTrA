#!/bin/bash -l
for m in BF4 C1 Me;
do
    grep $m DoCIndices-all.out | sort -k4n > ${m}_DoCIndices-all.out
done
