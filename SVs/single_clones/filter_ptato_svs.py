#!/usr/bin/python

import vcf as pyvcf
import glob

vcf_files = glob.glob("/hpc/pmc_vanboxtel/projects/TallFlow/3_Output/PTATO_Dev/*/intermediate/svs/Integration/*/*/*TALL**filtered.vcf")

for vcf_fname in vcf_files:
  vcf_reader = pyvcf.Reader(open(vcf_fname, 'r'))
  for record in vcf_reader:
      chr1 = record.CHROM
      pos = record.POS
      chr2 = record.ALT[0].chr
      pos2 = record.ALT[0].pos
      ori = record.ALT[0].orientation
      ori2 = record.ALT[0].remoteOrientation

      #print( chr1, pos, chr2, pos2, ori, ori2 )
      print( chr1, int(pos)-1, pos, ori )
      print( chr2, int(pos2)-1, pos2, ori2 )
