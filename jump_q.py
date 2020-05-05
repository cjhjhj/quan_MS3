#!/usr/bin/python

import os, sys, re, numpy as np, pandas as pd
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


def getParams(paramFile):
    parameters = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line)  # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line)  # Remove all whitespaces

            # Exception for "feature_files" parameter
            if "feature_files" in parameters and line.endswith("feature"):
                parameters["feature_files"].append(line)
            else:
                key = line.split('=')[0]
                val = line.split('=')[1]
                if key == "feature_files":
                    parameters[key] = [val]
                else:
                    parameters[key] = val
    return parameters


def parseIdtxt(idTxt, params):
    with open(idTxt) as f:
        jumpfPath = f.readline()
    jumpfPath = re.sub("misc/idsum.db", "", jumpfPath.split("=")[1]).strip()
    df = pd.read_table(idTxt, sep=";", skiprows=0, header=1)
    df = df.groupby(["Peptide", "Protein", "Outfile"]).size().reset_index()
    psms = {}
    pep2psm = {}
    prot2psm = {}
    for i in range(df.shape[0]):
        pep = df["Peptide"][i]
        prot = df["Protein"][i]
        psm = df["Outfile"][i]
        fraction = re.sub(r"(.\d+)$", ".mzXML", os.path.dirname(psm))

        if fraction in psms:
            psms[fraction].append(psm)
            psms[fraction] = list(set(psms[fraction]))
        else:
            psms[fraction] = [psm]

        if pep in pep2psm:
            pep2psm[pep].append(psm)
            pep2psm[pep] = list(set(pep2psm[pep]))
        else:
            pep2psm[pep] = [psm]

        if prot in prot2psm:
            prot2psm[prot].append(psm)
            prot2psm[prot] = list(set(prot2psm[prot]))
        else:
            prot2psm[prot] = [psm]

    '''
    df = pd.read_table(idTxt, sep=";", skiprows=0, header=1)
    df = df.groupby(["Peptide", "Protein", "Outfile"]).size().reset_index()
    psms = sorted(list(set(df["Outfile"]))) # Unique PSMs (i.e. MS2 scan numbers)
    pep2psm = {}
    prot2psm = {}
    progress = progressBar(df.shape[0])
    for i in range(df.shape[0]):
        progress.increment()
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
    '''

    return psms, pep2psm, prot2psm, jumpfPath



def matchMs2ToMs3(psms, reader):
    output = {} # Dictionary whose key = MS2 scan number (string), value = corresponding MS3 scan number (string)
    progress = progressBar(len(psms))
    for psm in psms:
        progress.increment()
        ms2Num = int(os.path.basename(psm).split(".")[-3])  # It should be changed to "int" here to deal with numbers such as "09999"
        ms2PrecMz = reader[str(ms2Num)]["precursorMz"][0]["precursorMz"]

        # Find corresponding MS3 scan number
        i = ms2Num + 1
        spec = reader[str(i)]
        while spec["msLevel"] > 1:
            if spec["msLevel"] == 3 and (ms2PrecMz - 0.1) < spec["precursorMz"][-1]["precursorMz"] < (ms2PrecMz + 3):
                ms3Num = spec["num"]
                output[str(ms2Num)] = ms3Num
                break
            i += 1
            spec = reader[str(i)]

    return output

'''
def matchMs2ToMs3(psms):
    output = {} # Dictionary whose key = MS2 scan number (string), value = corresponding MS3 scan number (string)
    for fraction in psms.keys():
        output[fraction] = {}
        reader = mzxml.MzXML(fraction)
        for psm in psms[fraction]:
            ms2Num = int(os.path.basename(psm).split(".")[-3])  # It should be changed to "int" here to deal with numbers such as "09999"
            ms2PrecMz = reader[str(ms2Num)]["precursorMz"][0]["precursorMz"]

            # Find corresponding MS3 scan number
            i = ms2Num + 1
            spec = reader[str(i)]
            while spec["msLevel"] > 1:
                if spec["msLevel"] == 3 and (ms2PrecMz - 0.1) < spec["precursorMz"][-1]["precursorMz"] < (ms2PrecMz + 3):
                    ms3Num = spec["num"]
                    output[fraction][ms2Num] = ms3Num
                    break
                i += 1
                spec = reader[str(i)]

    return output
'''

def getReporterMz(name):
    if name == "sig126":
        return 126.127726
    elif name == "sig127" or name == "sig127N":
        return 127.124761
    elif name == "sig127C":
        return 127.131081
    elif name == "sig128N":
        return 128.128116
    elif name == "sig128" or name == "sig128C":
        return 128.134436
    elif name == "sig129" or name == "sig129N":
        return 129.131471
    elif name == "sig129C":
        return 129.137790
    elif name == "sig130N":
        return 130.134825
    elif name == "sig130" or name == "sig130C":
        return 130.141145
    elif name == "sig131" or name == "sig131N":
        return 131.138180
    elif name == "sig131C":
        return 131.144500
    elif name == "sig132N":
        return 132.141535
    elif name == "sig132C":
        return 132.147855
    elif name == "sig133N":
        return 133.144890
    elif name == "sig133C":
        return 133.151210
    elif name == "sig134N":
        return 134.148245


def getReporterIntensity(scanNum, reader, params, **kwargs):
    tol = 10
    spec = reader[scanNum]
    reporterNames = params["tmt_reporters_used"].split(";")

    # Extract reporter ion m/z and intensity values
    mzArray = []
    intensityArray = []
    for reporter in reporterNames:
        if reporter in kwargs:
            mz = getReporterMz(reporter) * (1 + kwargs[reporter] / 1e6)
        else:
            mz = getReporterMz(reporter)
        lL = mz - mz * tol / 1e6
        uL = mz + mz * tol / 1e6
        ind = np.where((spec["m/z array"] >= lL) & (spec["m/z array"] <= uL))[0]
        if len(ind) == 0:
            reporterMz = 0
        elif len(ind) == 1:
            ind = ind[0]
            reporterMz = spec["m/z array"][ind]
        elif len(ind) > 1:
            ind2 = np.argmax(spec["intensity array"][ind])
            ind = ind[ind2]
            reporterMz = spec["m/z array"][ind]
        if lL <= reporterMz < uL:
            reporterIntensity = spec["intensity array"][ind]
        else:
            reporterIntensity = 0
        mzArray.append(reporterMz)
        intensityArray.append(reporterIntensity)

    '''
    # Extract K-TMT and R-y1 ion intensities
    tmtNo = re.search("\d+", os.path.basename(params["impurity_matrix"])).group(0)
    if tmtNo == "16":
        Kmz = 451.3199494907  # m/z of K-TMT-y1 ion for TMT16 (= 147.1128041645 + 304.2071453 (TMTpro-balancer))
    else:
        Kmz = 376.2757362992;  # m/z of K-TMT-y1 ion for TMT10, 11 (= 147.1128... + 229.1629321 (TMT-balancer))

    Kintensity = 0
    KlL = Kmz - Kmz * tol / 1e6
    KuL = Kmz + Kmz * tol / 1e6
    ind = np.where((spec["m/z array"] >= KlL) & (spec["m/z array"] <= KuL))[0]
    if len(ind) > 1:
        ind2 = np.argmax(spec["intensity array"][ind])
        ind = ind[ind2]
        ind = ind[0]
        Kintensity = spec["intensity array"][ind]

    Rintensity = 0
    Rmz = 175.1189521741
    RlL = Rmz - Rmz * tol / 1e6
    RuL = Rmz + Rmz * tol / 1e6
    ind = np.where((spec["m/z array"] >= RlL) & (spec["m/z array"] <= RuL))[0]
    if len(ind) > 1:
        ind2 = np.argmax(spec["intensity array"][ind])
        ind = ind[ind2]
        ind = ind[0]
        Rintensity = spec["intensity array"][ind]

    outArray = mzArray + intensityArray + [Kintensity, Rintensity]
    '''

    outArray = mzArray + intensityArray
    return outArray


def getSubset(df, params):
    # Get a subset of dataframe to calculate loading-bias information
    noiseLevel = 1000
    reporters = params["tmt_reporters_used"].split(";")
    snRatio = float(params["SNratio_for_correction"])
    subDf = df[reporters][(df[reporters] > noiseLevel * snRatio).prod(axis=1).astype(bool)]  # Zero-intensity PSMs are excluded
    subDf = np.log2(subDf.divide(subDf.mean(axis=1), axis=0))
    pctTrimmed = float(params["percentage_trimmed"])
    n = 0
    for reporter in reporters:
        if n == 0:
            ind = ((subDf[reporter] > subDf[reporter].quantile(pctTrimmed / 200)) &
                   (subDf[reporter] < subDf[reporter].quantile(1 - pctTrimmed / 200)))
        else:
            ind = ind & ((subDf[reporter] > subDf[reporter].quantile(pctTrimmed / 200)) &
                         (subDf[reporter] < subDf[reporter].quantile(1 - pctTrimmed / 200)))
        n += 1

    return subDf.loc[ind]


def getLoadingBias(df, params):
    ############################
    # Loading-bias calculation #
    ############################
    subDf = getSubset(df, params)
    n = subDf.shape[0]
    avg = 2 ** subDf.mean(axis=0) * 100
    # To express SD as percentage in raw-scale, take the mean deviation of 2^sd and 2^(-sd)
    # For example, sdVal = 0.2 for a specific reporter which represent +/-0.2 in log2-scale
    # 2^0.2 = 1.1487, 2^(-0.2) = 0.87
    # positive variation = 2^0.2 - 1 = 0.1487 (assuming the center is ideally 1 (=2^0))
    # negative variation = 1 - 2^(-0.2) = 0.13
    # mean variation (in raw-scale) = (0.1487 + 0.13) / 2 = 0.1394
    # Let the standard deviation of the reporter be 13.94% (= 0.1394 * 100)
    sdVal = subDf.std(axis=0)
    sd = ((2 ** sdVal - 1) + (1 - 2 ** (-sdVal))) / 2 * 100
    sem = sd / np.sqrt(n)
    return avg, sd, sem, n


def normalization(df, params):
    ################################################
    # Normalization (i.e. loading-bias correction) #
    ################################################
    doNormalization = params["loading_bias_correction"]
    normalizationMethod = params["loading_bias_correction_method"]
    if doNormalization == "1":
        # First, get a subset for calculating normalization factors (same as loading-bias calculation)
        subDf = getSubset(df, params)
        # Calculate normalization factors for reporters
        print("  Normalization is being performed")
        if normalizationMethod == "1":  # Trimmed-mean
            normFactor = subDf.mean(axis=0) - np.mean(subDf.mean())
        elif normalizationMethod == "2":  # Trimmed-median
            normFactor = subDf.median(axis=0) - np.mean(subDf.median())

        # Normalize the dataframe, df
        psmMean = df[reporters].mean(axis=1)
        df[reporters] = np.log2(df[reporters].divide(psmMean, axis=0).replace(0, np.nan))
        df[reporters] = df[reporters] - normFactor
        df[reporters] = 2 ** df[reporters]
        df[reporters] = df[reporters].multiply(psmMean, axis=0).replace(np.nan, 0)

        # Show the loading-bias information after normalization
        avg, sd, sem, n = getLoadingBias(df, params)
        print("  Loading bias (after correction)")
        print("  Reporter\tMean[%]\tSD[%]\tSEM[%]\t#PSMs")
        for i in range(len(reporters)):
            print("  %s\t%.2f\t%.2f\t%.2f\t%d" % (reporters[i], avg[i], sd[i], sem[i], n))
    else:
        print("  Normalization is skipped according to the parameter")

    return df


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


# mzxmlFile = "NCI-11plex-1-F1-f10268.mzXML"
paramFile = sys.argv[1]
params = getParams(paramFile)
idTxt = params["idtxt"]
print("  Loading ID.txt file")
psms, pep2psm, prot2psm, jumpfPath = parseIdtxt(idTxt, params)

####################################################
# Extract TMT reporter ion intensities - 1st round #
####################################################
print("  Extraction of TMT reporter ion intensities")
ms2ToMs3 = {}
qdict = {}  # Dictionary; key = MS2 scan number, value = TMT reporter intensities (array)
for frac in sorted(psms.keys()):
    reader = mzxml.MzXML(frac)
    print("  Processing %s" % os.path.basename(frac))
    print("  Looking for MS2 precursor scans of MS3 scans")
    ms2ToMs3[frac] = matchMs2ToMs3(psms[frac], reader)
    print("  Reporter intensities are being extracted from MS3 scans")
    progress = progressBar(len(ms2ToMs3[frac]))
    for ms2, ms3 in ms2ToMs3[frac].items():
        progress.increment()
        reporterIntensity = getReporterIntensity(ms3, reader, params)
        key = os.path.basename(frac) + "_" + ms2
        qdict[key] = reporterIntensity
    print()

# Create a dataFrame after the first extraction of reporter m/z and intensity values
reporters = params["tmt_reporters_used"].split(";")
columnNames = [re.sub("sig", "mz", i) for i in reporters] + reporters
df = pd.DataFrame.from_dict(qdict, orient = 'index', columns = columnNames) # Pandas dataframe of TMT reporter intensities

# Calculate m/z shift of each reporter ion (to be used for refining reporter intensities)
mzShifts = {}
nTot = df.shape[0]
print("  Summary of quantified TMT reporter ions")
for reporter in reporters:
    reporterMz = getReporterMz(reporter)
    measuredMz = df[re.sub("sig", "mz", reporter)][df[re.sub("sig", "mz", reporter)] > 0]   # Exclude un-quantified ones
    n = len(measuredMz)
    print("  %s\t%d (%.2f%%) matched" % (reporter, n, n / nTot * 100))
    mzShift = ((measuredMz - reporterMz) / reporterMz * 1e6).mean()
    mzShifts[reporter] = mzShift
print()

##########################################################
# Refine the reporter intensities considering m/z shifts #
##########################################################
print("  Refined extraction of reporter ion intensities")
qdict = {}
for frac in sorted(psms.keys()):
    print("  Processing %s" % os.path.basename(frac))
    reader = mzxml.MzXML(frac)
    for ms2, ms3 in ms2ToMs3[frac].items():
        reporterIntensity = getReporterIntensity(ms3, reader, params, **mzShifts)
        key = os.path.basename(frac) + "_" + ms2
        qdict[key] = reporterIntensity
print()

df = pd.DataFrame.from_dict(qdict, orient = 'index', columns = columnNames)
nTot = df.shape[0]
print("  Summary of quantified TMT reporter ions")
for reporter in reporters:
    reporterMz = getReporterMz(reporter)
    measuredMz = df[re.sub("sig", "mz", reporter)][df[re.sub("sig", "mz", reporter)] > 0]   # Exclude un-quantified ones
    n = len(measuredMz)
    print("  %s\t%d (%.2f%%) matched" % (reporter, n, n / nTot * 100))
print()

###########################
# TMT impurity correction #
###########################
dfImpurity = pd.read_table(params["impurity_matrix"], sep="\t", skiprows=1, header=None, index_col = 0)
dfImpurity = pd.DataFrame(np.linalg.pinv(dfImpurity.values), dfImpurity.columns, dfImpurity.index)
dfImpurity.columns = [str("sig") + i for i in dfImpurity.columns]
df2 = df[reporters].dot(dfImpurity.T)
df2.columns = reporters
df2[df2 <= 0] = 0 # Lower bound to 0
df[reporters] = pd.concat([df[reporters]/2, df2]).max(level = 0) # Upper bound to the half of original values (i.e. before correction)


'''
# Generate raw_..._(non)zero.txt files
# For this implementation, 
# 1. No zero-intensity filtering (keep all zero-intensity PSMs)
# 2. No PPI filtering (PPI information is not available in Comet results)
filterMethods = params["min_intensity_method"].split(",")
filterCutoffs = params["min_intensity_value"].split(",")
for i in range(len(filterMethods)):
    if filterMethods[i] == "1":
        df = df[df[reporters].min(axis = 1) >= float(filterCutoffs[i])]
    elif filterMethods[i] == "2":
        df = df[df[reporters].max(axis = 1) >= float(filterCutoffs[i])]
    elif filterMethods[i] == "3":
        df = df[df[reporters].mean(axis = 1) >= float(filterCutoffs[i])]
    elif filterMethods[i] == "4":
        df = df[df[reporters].median(axis = 1) >= float(filterCutoffs[i])]

# For this implementation, min_intensity_..._psm_1_2_filtered is not considered
'''
#####################################
# Show the loading-bias information #
#####################################
avgBias, sdBias, semBias, nn = getLoadingBias(df, params)
print("  Loading bias (before correction)")
print("  Reporter\tMean[%]\tSD[%]\tSEM[%]\t#PSMs")
for i in range(len(reporters)):
    print("  %s\t%.2f\t%.2f\t%.2f\t%d" % (reporters[i], avgBias[i], sdBias[i], semBias[i], nn))
print()

#################
# Normalization #
#################
df = normalization(df, params)
print ()

#####################################
# Summarization of PSMs to proteins #
#####################################
'''
pepDict = {}
for pep, psms in pep2psm.items():
    psms = [i.split(".")[-3] for i in psms]
    subDf = df[reporters][df.index.isin(psms)]    
    if subDf.empty:
        continue
    else:
        # 1. Take log2 again for the summarization
        # 2. Representative mean values of PSMs
        # 3. Mean-centering across reporters and then take mean over PSMs
        # 4. Add back the representative mean 
        subDf = np.log2(df[reporters][df.index.isin(psms)].replace(0, np.nan))
        psmMean = subDf.mean(axis = 1).sort_values(ascending = False)
        if len(psmMean) > 3:
            repMean = np.mean(psmMean[0:3])
        else:
            repMean = np.mean(psmMean)
        pepDict[pep] = (2 ** (subDf.sub(psmMean, axis = 0).mean(axis = 0) + repMean)).replace(np.nan, 0).tolist()
'''
print("  Summarization of PSM-level to protein-level quantities")
protDict = {}
progress = progressBar(len(prot2psm))
for prot, psms in prot2psm.items():
    progress.increment()
    psms = [os.path.basename(i).split(".")[0] + ".mzXML" + "_" + os.path.basename(i).split(".")[-3] for i in psms]
    subDf = df[reporters][df.index.isin(psms)]
    if subDf.empty:
        continue
    else:
        # 1. Take log2 again for the summarization
        # 2. Representative mean values of PSMs
        # 3. Mean-centering across reporters and then take mean over PSMs
        # 4. Add back the representative mean
        subDf = np.log2(df[reporters][df.index.isin(psms)].replace(0, np.nan))
        psmMean = subDf.mean(axis = 1).sort_values(ascending = False)
        if len(psmMean) > 3:
            repMean = np.mean(psmMean[0:3])
        else:
            repMean = np.mean(psmMean)
        protDict[prot] = (2 ** (subDf.sub(psmMean, axis = 0).mean(axis = 0) + repMean)).replace(np.nan, 0).tolist()

protDf = pd.DataFrame.from_dict(protDict, orient = "index", columns = reporters)
print()

#########################
# Generate output files #
#########################
infoColumns = ["ProteinDescription", "GN", r"PSM#", r"TotalPeptide#", r"UniquePeptide#", r"%Coverage",
               "PeptideoftheHighestScore", r"Run#", r"Scan#", "m/z", "z", "ppm", "Xcorr", "dCn", r"Q-value(%)"]

# All proteins (according to Jump -f publication table)
idAllProtTxt = jumpfPath + "publications/id_all_prot.txt"
df = pd.read_table(idAllProtTxt, sep = "\t", skiprows = 0, header = 1, index_col = 1)   # Protein acession # as index
df.columns = df.columns.str.replace(' ', '')    # Remove space in column names
allProtDf = pd.concat([df[infoColumns], protDf], axis = 1, join = "inner")
allProtDf.index.name = "ProteinAccession"
allProtDf.to_csv(r"all_prot_quan.txt", sep = "\t")

# Unique proteins (according to Jump -f publication table)
idUniProtTxt = jumpfPath + "publications/id_uni_prot.txt"
df = pd.read_table(idUniProtTxt, sep = "\t", skiprows = 0, header = 1, index_col = 1)
df.columns = df.columns.str.replace(' ', '')    # Remove space in column names
uniProtDf = pd.concat([df[infoColumns], protDf] , axis = 1, join = "inner")
uniProtDf.index.name = "ProteinAccession"
uniProtDf.to_csv(r"uni_prot_quan.txt", sep = "\t")
