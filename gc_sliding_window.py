#!/usr/bin/python

import re
import string
import os
import sys

n = sys.argv[2]
result = []
headerCount = 0

def main():
  file = open(sys.argv[1],'r')
  if line[0] == ">":
    headerCount += 1
  else:
    line = line.strip("\n")
    for i in range(len(line)):
      start = max(0,i-n)
      stop = min(len(line),i+n) + 1
      window = line[start:stop]
      GC = (window.count('G') + window.count('C')) / float(len(window))
      result.append(GC)
      return result
