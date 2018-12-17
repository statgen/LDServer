#!/usr/bin/env python3
import gzip
import numpy as np

with gzip.open("chr21.test.vcf.gz","rt") as fp:
  for line in fp:
    if line.startswith("#CHROM"):
      ls = line.split()
      samples = ls[9:]
      break

columns = {
  "iid": samples,
  "fid": samples,
  "patid": [0 for _ in range(len(samples))],
  "matid": [0 for _ in range(len(samples))],
  "sex": np.random.randint(1,3,len(samples)),
  "rand_binary": np.random.randint(1,3,len(samples)),
  "rand_qt": np.random.random(len(samples))
}

# Fake some missing values for testing
columns["rand_binary"][0] = 0

with open("test.ped","wt") as out:
  for i in range(len(samples)):
    line = "\t".join([str(v[i]) for v in columns.values()])
    print(line,file=out)

with open("test.dat","wt") as out:
  print("A RAND_BINARY",file=out)
  print("T RAND_QT",file=out)
