#!/bin/bash

echo "Start"
date
salmon index -t ../GENOMES/RAW/XTT22.v2023.cds.fa -i ../GENOMES/SALMON/XTT22_index --gencode
echo "End:"
date
