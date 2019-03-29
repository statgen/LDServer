#!/bin/bash
raremetalworker --ped chr21.test.ped --dat chr21.test.dat --vcf chr21.test.vcf.gz --traitName RAND_QT --prefix chr21.test
raremetalworker --ped chr21.test.missing_values.ped --dat chr21.test.dat --vcf chr21.test.vcf.gz --traitName RAND_QT --prefix chr21.test.missing_values
