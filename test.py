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


# # Correct relationship from mzML file
# file = "NCI-11plex-1-F1-f10268.mzML"
# reader = mzml.MzML(file)
# # with reader:
# #     for spec in reader:
# #         if spec["ms level"] == 2:
# #             precMz = spec["precursorList"]["precursor"][0]["isolationWindow"]["isolation window target m/z"]
# #             precMz2 = spec["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
# #             if abs(precMz - precMz2) >= 0.1:
# #                 print("scan#%d :gap = %f" % (spec["index"], precMz - precMz2))
#
# f = open("MS2_MS3_mzML.txt", "w")
# progress = progressBar(len(reader))
# with reader:
#     precMzToMs2 = {}
#     ms2ToMs3 = {}
#     for spec in reader:
#         progress.increment()
#         if spec["ms level"] == 2:
#             precMz = spec["precursorList"]["precursor"][0]["isolationWindow"]["isolation window target m/z"]
#             precMzToMs2[precMz] = {}
#             precMzToMs2[precMz]["num"] = spec["index"] + 1
#             precMzToMs2[precMz]["diff"] = precMz - spec["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"]
#         elif spec["ms level"] == 3:
#             precMz = spec["precursorList"]["precursor"][-1]["isolationWindow"]["isolation window target m/z"]
#             if precMz in precMzToMs2:
#                 ms2ToMs3[precMzToMs2[precMz]["num"]] = spec["index"] + 1
#                 f.write("%d,%d,%f\n" % (precMzToMs2[precMz]["num"], spec["index"] + 1, precMzToMs2[precMz]["diff"]))
#             else:
#                 f.write("MS3 #%d has no precursor MS2\n" % spec["index"] + 1)
#
# f.close()

file = "./sum_NCI-11Plex-1-F1-f10268/ID.txt"
df = pd.read_table(file, sep = ";", skiprows = 0, header = 1)
uniquePsms = list(set(df["Outfile"]))
pep2psm = {}
prot2psm = {}
for i in range(df.shape[0]):
    pep = df["Peptide"][i]
    prot = df["Protein"][i]
    psm = df["Outfile"][i]
    if pep in pep2psm:
        pep2psm[pep].append(psm)
    else:
        pep2psm[pep] = [psm]
    if prot in prot2psm:
        prot2psm[prot].append(psm)
    else:
        prot2psm[prot] = [psm]


# df_pep_psm = df.groupby(["Peptide", "Outfile"]).size().reset_index()
# pep2psm = {}
# for i in range(df_pep_psm.shape[0]):
#     pep = df_pep_psm["Peptide"][i]
#     psm = df_pep_psm["Outfile"][i]
#     if pep in pep2psm:
#         pep2psm[pep].append(psm)
#     else:
#         pep2psm[pep] = [psm]

# Unique MS2 scan numbers
ms2 = [os.path.basename(i).split(".")[-3] for i in df["Outfile"]]
ms2 = list(set(ms2))