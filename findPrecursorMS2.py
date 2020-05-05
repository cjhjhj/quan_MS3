#!/usr/bin/python

import os, sys, numpy as np, pandas as pd
from pyteomics import mzxml, mzml


class progressBar:
    def __init__(self, total):
        self.total = total
        self.barLength = 20
        self.count = 0
        self.progress = 0
        self.block = 0
        self.status = ""

    def increment(self):
        self.count += 1
        self.progress = self.count / self.total
        self.block = int(round(self.barLength * self.progress))
        if self.progress == 1:
            self.status = "Done...\r\n"
        else:
            self.status = ""
        #         self.status = str(self.count) + "/" + str(self.total)
        text = "\r  Progress: [{0}] {1}% {2}".format("#" * self.block + "-" * (self.barLength - self.block),
                                                   int(self.progress * 100), self.status)
        sys.stdout.write(text)
        sys.stdout.flush()


# Input files
# 1. ID.txt: Jump -f result
# 2. mzXML file corresponding to the above ID.txt
# 3. Parameter file

# Extract quantity information from mzXML file
# << MSconvert-based mzXML >>
# A mzXML file from MSconvert has the following characteristics
# 1. MS3 scan is not always MS2 + 1 scan (it depends on MS instrument, not on MSconvert)
#    For example,
#    scan#100: MS2 -> scan#101: MS3 (generally)
#    scan#1000: MS2, scan#1001: MS2 -> scan#1003:MS3, scan#1004: MS3 (some cases)
# 2. Some MS2-MS3 pairs do not show the identical precursor m/z
#    For example,
#    scan#5000: MS2 -> precursor m/z = 1000.00 (in .raw file, precursor m/z = 1001.00)
#    scan#5001: MS3 -> precursor m/z = 1001.00 (in .raw file, precursor m/z = 1001.00)
#    MSconvert seems to re-evaluate the precursor m/z at MS2 scan (deisotoping?)
#    It needs to be considered in the script

# << ReAdW-based mzXML >>
# A mzXML file from ReAdW has the following characteristics
# 1. MS3 scan is not always MS2 + 1 scan (it depends on MS instrument, not on MSconvert)
#    For example,
#    scan#100: MS2 -> scan#101: MS3 (generally)
#    scan#1000: MS2, scan#1001: MS2 -> scan#1003:MS3, scan#1004: MS3 (some cases)
# 2. "msLevel" of MS2 scan is set to 0
# 3. There's no tag representing precursor m/z in MS2 and MS3
# 4. Precursor m/z can be inferred from "filterLine" tag in MS2 and MS3
# 4. Precursor m/z value is identifcal to .raw file (no re-evaluation/re-calculation)


# Inferred relationship between MS2 and MS3 using mzXML file
mzxmlFile = "NCI-11plex-1-F1-f10268.mzXML"
reader = mzxml.MzXML(mzxmlFile)
nTotScans = len(reader)
nMS1, nMS2, nMS3 = 0, 0, 0
progress = progressBar(nTotScans)
f = open("MS2_MS3_mzXML.txt", "w")
with reader:
    ms2ToMs3 = {}
    for spec in reader:
        progress.increment()
        if spec["msLevel"] == 1:
            nMS1 += 1
            precMzToMs2 = {}    # This dictionary is re-initiated for every MS1-scan cycle
        if spec["msLevel"] == 2:
            nMS2 += 1
            precMz = spec["precursorMz"][0]["precursorMz"]
            scanNum = spec["num"]
            precMzToMs2[precMz] = scanNum
        elif spec["msLevel"] == 3:
            nMS3 += 1
            precMzMs3 = spec["precursorMz"][-1]["precursorMz"]
            scanNum = spec["num"]
            if precMzMs3 in precMzToMs2:
                ms2ToMs3[precMzToMs2[precMzMs3]] = scanNum
                f.write("%s,%s,0\n" % (precMzToMs2[precMzMs3], scanNum))
            else:
                # Look for a previous MS2 scan having a similar precursor m/z to MS3 scan
                i = int(scanNum) - 1
                flag = False
                while reader[str(i)]["msLevel"] > 1:
                    iSpec = reader[str(i)]
                    if iSpec["msLevel"] == 2:
                        if (precMzMs3 - 3) <= iSpec["precursorMz"][0]["precursorMz"] <= (precMzMs3 + 0.1):
                            ms2ToMs3[iSpec["num"]] = scanNum
                            flag = True
                            f.write("%s,%s,%f\n" % (iSpec["num"], scanNum, precMzMs3 - iSpec["precursorMz"][0]["precursorMz"]))
                            # gap = int(scanNum) - i
                            # print("Gap = %d scan at scan#%s, difference = %f" % (gap, scanNum, precMz - iSpec["precursorMz"][0]["precursorMz"]))
                            break
                        else:
                            i = i - 1
                    else:
                        i = i - 1
                if not flag:
                    f.write("MS3 #%s has no precursor MS2\n" % scanNum)

f.close()
print(nMS1, nMS2, nMS3)