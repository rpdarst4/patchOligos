patchOligos
===========

Python object class to quickly look up oligo sequences for patch-BGS

EXAMPLE
===========

import patchOligos
DB=patchOligos.patch_oligos()

# given a file with a list of RefSeq UIDs in the last column
with open('yourOutput.tsv','w') as handle:
  handle.write('UID\tchromosome\tbasepair\tstrand\tA\tB\tsequence\n')
  for line in open('yourFile.tsv','r'):
    # get UID from last column, use as index for patch_oligos object DB
    for i in DB[line.rstrip('\n').split('\t')[-1]]:
      handle.write('\n'+'\t'.join([str(j) for j in i]))
      
  # another locus of interest, but not a TSS
  DB.add_locus('insulator01','chr3',36000000,'+')
  for i in DB['insulator01']:
    handle.write('\n'+'\t'.join([str(j) for j in i]))

(c) 2014, Russell Darst, University of Florida 
