#!/usr/bin/Python

import sys, os, re, numpy as np, time
from pyteomics import mzxml

def getParams(paramFile):
    params = dict()
    with open(paramFile, 'r') as file:
        for line in file:
            if re.search(r'^#', line) or re.search(r'^\s', line):
                continue
            line = re.sub(r'#.*', '', line) # Remove comments (start from '#')
            line = re.sub(r'\s*', '', line) # Remove all whitespaces
            key = line.split('=')[0]
            val = line.split('=')[1]
            params[key] = val
    params['save_dir'] = "quan_" + params['save_dir']    
    return params

def parseIdtxt(idTxt, params):
    with open(idTxt, "r") as file:
        outDict = dict()
        pepDict = dict()
        for line in file:
            line = line.strip() # Similar to chomp in Perl
            line = re.sub('//', '/', line) # Change all (possible) double slashes to single slashes
            # First two lines are 1. database description and 2. header
            # So, those two lines may need to be skipped
            if re.search(r'^Database', line):
                continue
            elif re.search(r'^Peptide', line):
                header = line.split(';')
                # The following approach may be more robust than direct indexing (e.g. peptide = elems[0])
                peptideInd = header.index('Peptide')
                proteinInd = header.index('Protein')
                outfileInd = header.index('Outfile')
                measmhInd = header.index('measuredMH')
                calcmhInd = header.index('calcMH')
                xcorrInd = header.index('XCorr')
                dcnInd = header.index('dCn')
                ppiInd = header.index('precursor_peak_intensity_percentage')
            else:
                # Retrieve PSM information
                elems = line.split(';')
                peptide = elems[peptideInd]
                protein = elems[proteinInd]
                outfile = elems[outfileInd]
                measmh = elems[measmhInd]
                calcmh = elems[calcmhInd]
                xcorr = elems[xcorrInd]
                dcn = elems[dcnInd]
                ppi = elems[ppiInd]
                ppi = re.sub('%', '', ppi)
                if float(ppi) < float(params['ppi_filter']):
                    continue # If PPI of a PSM is less than the user-defined threshold, skip it
                
                # Generate dictionaries for outfiles and peptides
                _, scanNum, _, charge, _ = os.path.basename(outfile).split('.') # Extract PSM information from the name of outfile
                mzxmlPath = re.sub('.\d+$', '.mzXML', os.path.dirname(outfile))
                leftAA, peptide, rightAA = peptide.split('.') # Peptide sequence information
                if mzxmlPath not in outDict:
                    outDict[mzxmlPath] = {}
                outDict[mzxmlPath][scanNum] = {}
                if "proteins" not in outDict[mzxmlPath][scanNum]:
                    outDict[mzxmlPath][scanNum]['proteins'] = []
                outDict[mzxmlPath][scanNum]['proteins'] = []
                outDict[mzxmlPath][scanNum]['leftAA'] = leftAA
                outDict[mzxmlPath][scanNum]['rightAA'] = rightAA
                outDict[mzxmlPath][scanNum]['peptide'] = peptide
                outDict[mzxmlPath][scanNum]['charge'] = charge
                outDict[mzxmlPath][scanNum]['MH+'] = measmh
                outDict[mzxmlPath][scanNum]['pepMH+'] = calcmh
                outDict[mzxmlPath][scanNum]['xcorr'] = xcorr
                outDict[mzxmlPath][scanNum]['dcn'] = dcn
                outDict[mzxmlPath][scanNum]['outfile'] = outfile
                outDict[mzxmlPath][scanNum]['proteins'].append(protein)
                if outfile not in pepDict:
                    pepDict[outfile] = {}
                pepDict[outfile]['peptide'] = peptide
                pepDict[outfile].setdefault('proteins', [])
                pepDict[outfile]['proteins'].append(protein)
                pepDict[outfile]['charge'] = charge
                pepDict[outfile]['endingAA'] = peptide[-1]
                pepDict[outfile]['nLabels'] = (peptide.count('K') + 1)

#                 if outfile not in pepDict:
#                     pepDict[outfile] = {}
#                     pepDict[outfile][peptide] = {}
#                     pepDict[outfile][peptide]['proteins'] = []
#                 pepDict[outfile][peptide]['proteins'].append(protein)
#                 pepDict[outfile][peptide]['charge'] = charge
#                 pepDict[outfile][peptide]['endingAA'] = peptide[-1]
#                 pepDict[outfile][peptide]['nLabels'] = (peptide.count('K') + 1)
    return outDict, pepDict


# def parsePepXML(pepXML, params):
#     with open(pepXML, "r") as file:
#         outDict = dict()
#         pepDict = dict()
#         for line in file:
#             line = line.strip()  # Similar to chomp in Perl
#             line = re.sub('//', '/', line)  # Change all (possible) double slashes to single slashes
#             # First two lines are 1. database description and 2. header
#             # So, those two lines may need to be skipped
#             if re.search(r'^Database', line):
#                 continue
#             elif re.search(r'^Peptide', line):
#                 header = line.split(';')
#                 # The following approach may be more robust than direct indexing (e.g. peptide = elems[0])
#                 peptideInd = header.index('Peptide')
#                 proteinInd = header.index('Protein')
#                 outfileInd = header.index('Outfile')
#                 measmhInd = header.index('measuredMH')
#                 calcmhInd = header.index('calcMH')
#                 xcorrInd = header.index('XCorr')
#                 dcnInd = header.index('dCn')
#                 ppiInd = header.index('precursor_peak_intensity_percentage')
#             else:
#                 # Retrieve PSM information
#                 elems = line.split(';')
#                 peptide = elems[peptideInd]
#                 protein = elems[proteinInd]
#                 outfile = elems[outfileInd]
#                 measmh = elems[measmhInd]
#                 calcmh = elems[calcmhInd]
#                 xcorr = elems[xcorrInd]
#                 dcn = elems[dcnInd]
#                 ppi = elems[ppiInd]
#                 ppi = re.sub('%', '', ppi)
#                 if float(ppi) < float(params['ppi_filter']):
#                     continue  # If PPI of a PSM is less than the user-defined threshold, skip it
#
#                 # Generate dictionaries for outfiles and peptides
#                 _, scanNum, _, charge, _ = os.path.basename(outfile).split(
#                     '.')  # Extract PSM information from the name of outfile
#                 mzxmlPath = re.sub('.\d+$', '.mzXML', os.path.dirname(outfile))
#                 leftAA, peptide, rightAA = peptide.split('.')  # Peptide sequence information
#                 if mzxmlPath not in outDict:
#                     outDict[mzxmlPath] = {}
#                 outDict[mzxmlPath][scanNum] = {}
#                 if "proteins" not in outDict[mzxmlPath][scanNum]:
#                     outDict[mzxmlPath][scanNum]['proteins'] = []
#                 outDict[mzxmlPath][scanNum]['proteins'] = []
#                 outDict[mzxmlPath][scanNum]['leftAA'] = leftAA
#                 outDict[mzxmlPath][scanNum]['rightAA'] = rightAA
#                 outDict[mzxmlPath][scanNum]['peptide'] = peptide
#                 outDict[mzxmlPath][scanNum]['charge'] = charge
#                 outDict[mzxmlPath][scanNum]['MH+'] = measmh
#                 outDict[mzxmlPath][scanNum]['pepMH+'] = calcmh
#                 outDict[mzxmlPath][scanNum]['xcorr'] = xcorr
#                 outDict[mzxmlPath][scanNum]['dcn'] = dcn
#                 outDict[mzxmlPath][scanNum]['outfile'] = outfile
#                 outDict[mzxmlPath][scanNum]['proteins'].append(protein)
#                 if outfile not in pepDict:
#                     pepDict[outfile] = {}
#                 pepDict[outfile]['peptide'] = peptide
#                 pepDict[outfile].setdefault('proteins', [])
#                 pepDict[outfile]['proteins'].append(protein)
#                 pepDict[outfile]['charge'] = charge
#                 pepDict[outfile]['endingAA'] = peptide[-1]
#                 pepDict[outfile]['nLabels'] = (peptide.count('K') + 1)
#
# return outDict, pepDict


def getMS2(outDict):
    ms2Dict = dict()
    for mzXML in outDict.keys():
        reader = mzxml.read("C:/" + mzXML)
        ms2Dict[mzXML] = dict()
        nPSMs = 0
        for spec in reader:
            scanNum = spec['id']
            if scanNum in outDict[mzXML].keys():
                nPSMs += 1
                ms2Dict[mzXML][scanNum] = dict()
                ms2Dict[mzXML][scanNum]['mz'] = spec['m/z array']
                ms2Dict[mzXML][scanNum]['intensity'] = spec['intensity array']
                ms2Dict[mzXML][scanNum]['rt'] = spec['retentionTime']
            elif int(scanNum) > int(max(outDict[mzXML].keys())):
                break
            else:
                continue
        ms2Dict[mzXML] = massCorrection(ms2Dict[mzXML])
    return (ms2Dict)

def massCorrection(ms2):
    # Input ms2 is a dict with scan numbers as keys
    # Calculation of mass shifts between observed and theoretical m/z
    refMz = 126.1277259380 # Theoretical mass of sig126
    lL = refMz - 50 * refMz / 1e6 # By default, +/- 50ppm tolerance is considered
    uL = refMz + 50 * refMz / 1e6
    mzShiftArray = []
    for scanNum in sorted(ms2):
        ind, = np.where((ms2[scanNum]['mz'] >= lL) & (ms2[scanNum]['mz'] <= uL))
        ind = ind[np.argmax(ms2[scanNum]['intensity'][ind])]
        obsMz = ms2[scanNum]['mz'][ind]
        if obsMz > 0:
            mzShift = (obsMz - refMz) / refMz * 1e6 # Unit of ppm
            mzShiftArray.append(mzShift)

    # Filter the calculated mass shifts
    # 1. More than 50% of spectra should be used
    # 2. Remove 20% of highest/lowest mass shifts to get a more reliable correction factor
    measuredPct = len(mzShiftArray) / len(ms2)
    if measuredPct >= 0.5:
        mzShiftArray.sort()
        indStart = int(0.2 * len(mzShiftArray))
        indEnd = len(mzShiftArray) - indStart
        mzShiftArray = mzShiftArray[indStart:indEnd]
        meanShift = np.mean(mzShiftArray)
        stdShift = np.std(mzShiftArray)
        print ("  Calculated mass-shift: mean = %.4f ppm and SD = %.4f ppm" % (meanShift, stdShift))
        for scanNum in ms2.keys():
            ms2[scanNum]['mz'] = np.divide(ms2[scanNum]['mz'], (1 + meanShift / 1e6))
    else:
        print ("  Warning")
        print ("  Only %.2f%% of spectra can be used to compute mass shifts" % measuredPct)
        print ("  Mass correction will not be performed when less than 50% of spectra is used")
    return (ms2)

def getReporterIntensity(outDict, pepDict, ms2Dict, reporters, params, reporterSummary = None):
    nReporters = len(reporters)
    for mzXML in outDict.keys():
        nPSMs = 0
        for scanNum in outDict[mzXML].keys():
            nPSMs += 1
            outfile = outDict[mzXML][scanNum]['outfile']            
            mzArray = ms2Dict[mzXML][scanNum]['mz']
            intensityArray = ms2Dict[mzXML][scanNum]['intensity']
            reporterMz = np.zeros(nReporters)
            reporterIntensity = np.zeros(nReporters)
            # Extract TMT reporter peaks (mz and intensity)
            for i in range(0, len(reporters)):
                if reporterSummary is None:
                    center = reporters[i]
                    tol = 10 # Hard-coded 10 ppm tolerance around the theoretical reporter m/z
                else:
                    center = reporters[i] + reporters[i] * reporterSummary['avg'][i] / 1e6
                    tol = reporterSummary['std'][i] * float(params['tmt_peak_extraction_second_sd'])
                lL = center - tol * center / 1e6
                uL = center + tol * center / 1e6
                ind, = np.where((mzArray >= lL) & (mzArray <= uL))
                if len(ind) == 0:
                    continue
                elif len(ind) > 1:
                    ind = ind[np.argmax(intensityArray[ind])]
                reporterMz[i] = mzArray[ind]
                reporterIntensity[i] = intensityArray[ind]
            pepDict[outfile]['mz'] = reporterMz
            pepDict[outfile]['intensity'] = reporterIntensity
            pepDict[outfile]['rt'] = ms2Dict[mzXML][scanNum]['rt']
            # Extract K-TMT and R-y1 ion peaks (for interference removal, in case)
            Kmz = 376.2757362992 # m/z of K-TMT-y1 ion
            KlL = Kmz - tol * Kmz / 1e6
            KuL = Kmz + tol * Kmz / 1e6
            Kintensity = 0
            Rmz = 175.1189521741 # m/z of R-y1 ion
            RlL = Rmz - tol * Rmz / 1e6
            RuL = Rmz + tol * Rmz / 1e6
            Rintensity = 0
            for i in range(0, len(mzArray)):
                if mzArray[i] < RlL:
                    continue
                elif mzArray[i] >= RlL and mzArray[i] <= RuL:
                    if intensityArray[i] >= Rintensity:
                        Rintensity = intensityArray[i]
                elif mzArray[i] >= KlL and mzArray[i] <= KuL:
                    if intensityArray[i] >= Kintensity:
                        Kintensity = intensityArray[i]
                elif mzArray[i] > KuL:
                    break
            pepDict[outfile]['Ry1'] = Rintensity
            pepDict[outfile]['Ky1'] = Kintensity
    return (outDict, pepDict)

def refineReporterIntensity(outDict, pepDict, ms2Dict, reporterMzs, params):
    nReporters = len(reporterMzs)
    reporters = params['tmt_reporters_used'].split(';')
    # Summary of PSMs in which reporters are identified    
    print ("  Reporter ion extraction summary for 1st round")
    nTot = len(pepDict)
    for i in range(0, nReporters):
        n = 0
        for outfile in pepDict.keys():
            if pepDict[outfile]['mz'][i] > 0:
                n += 1
        pct = (n / nTot) * 100
        print ("    %s \t%d (%.2f%%) matched" % (reporters[i], n, pct))
    print ()
    reporterSummary = {}
    reporterSummary['avg'] = []
    reporterSummary['std'] = []
    for i in range(0, nReporters):
        mzShiftArray = []
        for outfile in pepDict.keys():
            if pepDict[outfile]['mz'][i] > 0:
                mzShift = (pepDict[outfile]['mz'][i] - reporterMzs[i]) / reporterMzs[i] * 1e6
                mzShiftArray.append(mzShift)
        reporterSummary['avg'].append(np.mean(mzShiftArray))
        reporterSummary['std'].append(np.std(mzShiftArray))
        print ("    %s\tm/z shift = %.4f (ppm)\tsd = %.4f" % (reporters[i], abs(np.mean(mzShiftArray)), np.std(mzShiftArray)))
    print ()
    # Re-extraction of reporter ions (using calibrated mass)
    SDlimit = 100
    for i in range(0, nReporters - 1):
        mz1 = reporterMzs[i]
        mz2 = reporterMzs[i + 1]
        center1 = mz1 + mz1 * reporterSummary['avg'][i] / 1e6
        center2 = mz2 + mz2 * reporterSummary['avg'][i + 1] / 1e6
        coeff1 = mz1 * reporterSummary['std'][i] / 1e6
        coeff2 = mz2 * reporterSummary['std'][i + 1] / 1e6
        currSD = (center2 - center1) / (coeff1 + coeff2)
        if currSD < SDlimit:
            SDlimit = currSD
    print ("  According to the mass-shifts of reporter peaks, 'tmt_peak_extraction_second_sd' should be smaller than %.2f" % SDlimit)
    print ("  Current 'tmt_peak_extraction_second_sd' is %s" % params['tmt_peak_extraction_second_sd'])
    if float(params['tmt_peak_extraction_second_sd']) > SDlimit:
        print ("  == Warning ==")
        print ("  Reporter peak extraction may be inaccurate due to the large value of 'tmt_peak_extraction_second_sd'")
        print ("  Please set the parameter to the value smaller than %.2f" % SDlimit)
    print ()

    # Refine reporter extraction
    outDict, pepDict = getReporterIntensity(outDict, pepDict, ms2Dict, reporterMzs, params, reporterSummary)
    print ("  Reporter ion extraction summary for 2nd round")
    nTot = len(pepDict)
    for i in range(0, nReporters):
        n = 0
        for outfile in pepDict.keys():
            if pepDict[outfile]['mz'][i] > 0:
                n += 1
        pct = (n / nTot) * 100
        print ("    %s \t%d (%.2f%%) matched" % (reporters[i], n, pct))
    print ()
    return outDict, pepDict

def correctImpurity(pepDict, params):
    purityMatrix = getPurityMatrix(params['impurity_matrix'])
    inversePurityMatrix = np.linalg.inv(purityMatrix)
    for outfile in pepDict.keys():
        uncorrected = pepDict[outfile]['intensity']
        corrected = inversePurityMatrix.dot(uncorrected)
        # Prevention of over-correction
        # 1. Zero-intensity should not be over-corrected
        # 2. Negative intensity is not allowed
        # 3. If intensity looks over-corrected (i.e. smaller than the half of uncorrected intensity,
        #    the corrected value will be bound to the half of uncorrected one
        corrected[uncorrected == 0] = 0
        corrected[corrected < 0] = 0
        ind, = np.where(corrected < 0.5 * uncorrected)
        corrected[ind] = 0.5 * uncorrected[ind]
        pepDict[outfile]['intensity'] = corrected
    return (pepDict)

def getReporterMz(reporters):
    nReporters = len(reporters)
    reporterMzs = np.zeros(nReporters)
    for i in range(0, nReporters):
        if reporters[i] == "sig126":
            reporterMzs[i] = 126.127726
        elif reporters[i] == "sig127" or reporters[i] == "sig127N":
            reporterMzs[i] = 127.124761
        elif reporters[i] == "sig127C":
            reporterMzs[i] = 127.131081
        elif reporters[i] == "sig128N":
            reporterMzs[i] = 128.128116
        elif reporters[i] == "sig128" or reporters[i] == "sig128C":
            reporterMzs[i] = 128.134436
        elif reporters[i] == "sig129" or reporters[i] == "sig129N":
            reporterMzs[i] = 129.131471
        elif reporters[i] == "sig129C":
            reporterMzs[i] = 129.137790
        elif reporters[i] == "sig130N":
            reporterMzs[i] = 130.134825
        elif reporters[i] == "sig130" or reporters[i] == "sig130C":
            reporterMzs[i] = 130.141145
        elif reporters[i] == "sig131" or reporters[i] == "sig131N":
            reporterMzs[i] = 131.138180
        elif reporters[i] == "sig131C":
            reporterMzs[i] = 131.1445001
        else:
            sys.exit("  $reporters[$i] is an incorrect reporter name")
    return (reporterMzs)

def getPurityMatrix(inputFile):
    flag = 0
    nReporters = re.match('TMT(\d+).ini', os.path.basename(inputFile)).group(1)
    nReporters = int(nReporters)
    dataMatrix = np.empty([nReporters, nReporters])
    rowInd = 0
    if os.path.basename(inputFile) == "TMT10.ini":
        with open(inputFile) as file:
            for line in file:
                line = line.strip()
                if re.search('PurityMatrix_Batch2', line):
                    flag = 1
                    continue
                if (flag == 1) & (not re.search('\[', line)):
                    dataMatrix[rowInd, :] = line.split("\t")[1:]
                    rowInd += 1
                else:
                    flag = 0
    else:
        with open(inputFile) as file:
            for line in file:
                line = line.strip()
                if re.search('PurityMatrix', line):
                    flag = 1
                    continue
                if (flag == 1) & (not re.search('\[', line)):
                    dataMatrix[rowInd, :] = line.split("\t")[1:]
                    rowInd += 1
                else:
                    flag = 0
    return (dataMatrix)

def filterByIntensity(intensity, methods, cutoffs, isLog):
    nFilters = len(methods)
    isFiltered, methodFiltered = 0, 0
    for i in range(0, nFilters):
        if methods[i] == 1:
            minIntensity = np.amin(intensity)
            if isLog == 1:
                minIntensity = 2 ** minIntensity
            if minIntensity < cutoffs[i]:
                isFiltered = 1
                methodFiltered = "minimum"
                break
        elif methods[i] == 2:
            minIntensity = np.amax(intensity)
            if isLog == 1:
                minIntensity = 2 ** minIntensity
            if minIntensity < cutoffs[i]:
                isFiltered = 1
                methodFiltered = "maximum"
        elif methods[i] == 3:
            minIntensity = np.mean(intensity)
            if isLog == 1:
                minIntensity = 2 ** minIntensity
            if minIntensity < cutoffs[i]:
                isFiltered = 1
                methodFiltered = "mean"
        elif methods[i] == 4:
            minIntensity = np.median(intensity)
            if isLog == 1:
                minIntensity = 2 ** minIntensity
            if minIntensity < cutoffs[i]:
                isFiltered = 1
                methodFiltered = "median"
        elif methods[i] == 0:
            break
        else:
            sys.exit("  Minimum intensity based filtering failed. Please check parameters related with intensity-based filtering")
    return (isFiltered, methodFiltered)

def psmFilter1(outDict, pepDict, nReporters, filterMethods, filterCutoffs, 
              filtered, nonFiltered, zeroFiltered, intFiltered):
    nTotPSMs = len(pepDict.keys())
    nPSMs = 0
    for fraction in outDict.keys():
        for scanNum in outDict[fraction].keys():
            psm = outDict[fraction][scanNum]['outfile']
            nPSMs = nPSMs + 1
            print ("    Processing %d PSMs out of %d for the minimum intensity-based filtering\r" % (nPSMs, nTotPSMs));
            # Check whether there is any zero m/z and/or zero intensity peak            
            mzInd, = np.nonzero(pepDict[psm]['mz'])
            intenInd, = np.nonzero(pepDict[psm]['intensity'])
            if len(mzInd) < nReporters or len(intenInd) < nReporters:
                zeroFiltered.setdefault(psm, 0)
                zeroFiltered[psm] += 1
                filtered.setdefault(psm, 0)
                filtered[psm] += 1                
            # If there's no zero m/z or intensity peak, then intensity-based filtering will be performed
            else:
                isFiltered, methodFiltered = filterByIntensity(pepDict[psm]['intensity'], 
                                                               filterMethods, filterCutoffs, 0)
                if isFiltered == 1:
                    intFiltered.setdefault(methodFiltered, {})
                    intFiltered[methodFiltered].setdefault(psm, 0)
                    intFiltered[methodFiltered][psm] += 1
                    filtered.setdefault(psm, 0)
                    filtered[psm] += 1
                else:
                    nonFiltered[psm] = 1
    return (filtered, nonFiltered, zeroFiltered, intFiltered)


def psmFilter2(outDict, pepDict, nReporters, filterMethods, filterCutoffs,
               filtered, nonFiltered):
    # Create a protein-to-outfile mapping dictionary by going over unfiltered PSMs
    psm_1_2_filtered = {}
    prot2psm = {}
    for fraction in outDict.keys():
        for scanNum in outDict[fraction].keys():
            psm = outDict[fraction][scanNum]['outfile']
            if psm in nonFiltered:
                for protein in outDict[fraction][scanNum]['proteins']:
                    prot2psm.setdefault(protein, [])
                    prot2psm[protein].append(psm)
    
    # Check PSMs mapped to a protein
    nTotProteins = len(prot2psm.keys())
    nProteins = 0
    for protein in prot2psm.keys():
        nProteins += 1
        print ("    Processing %d proteins out of %d to filter out unreliable PSMs\r" % (nProteins, nTotProteins))
        psms = list(set(prot2psm[protein])) # Unique PSMs for each protein
        n = len(psms)
        if n == 1:
            # 1. Protein mapped by only one PSM
            #    If the PSM is filtered at the 1st filtering, this protein will not be quantified
            isFiltered, methodFiltered = filterByIntensity(pepDict[psms[0]]['intensity'],
                                                           filterMethods, filterCutoffs, 0)
            if isFiltered == 1:
                del nonFiltered[psms[0]]
                psm_1_2_filtered.setdefault(methodFiltered, {})
                psm_1_2_filtered[methodFiltered].setdefault(psms[0], 0)
                psm_1_2_filtered[methodFiltered][psms[0]] += 1
                filtered.setdefault(psms[0], 0)
                filtered[psms[0]] += 1                
        elif n == 2:
            # 2. Protein mapped by two PSMs
            isFiltered0, methodFiltered0 = filterByIntensity(pepDict[psms[0]]['intensity'],
                                                           filterMethods, filterCutoffs, 0)
            isFiltered1, methodFiltered1 = filterByIntensity(pepDict[psms[1]]['intensity'],
                                                           filterMethods, filterCutoffs, 0)
            # 2.1. When both PSMs are filtered out, the protein will not be quantified
            # 2.2. When one PSM is filtered out, the other PSM will be used to quantify the protein
            # 2.3. When both PSMs are NOT filtered out, then choose a PSM with either smaller variation or larger mean
            if isFiltered0 == 1 and isFiltered1 == 1:
                del nonFiltered[psms[0]], nonFiltered[psms[1]]
                psm_1_2_filtered.setdefault(methodFiltered0, {})
                psm_1_2_filtered[methodFiltered0].setdefault(psms[0], 0)
                psm_1_2_filtered[methodFiltered0][psms[0]] += 1
                filtered.setdefault(psms[0], 0)
                filtered[psms[0]] += 1
                psm_1_2_filtered.setdefault(methodFiltered1, {})
                psm_1_2_filtered[methodFiltered1].setdefault(psms[1], 0)
                psm_1_2_filtered[methodFiltered1][psms[1]] += 1
                filtered.setdefault(psms[1], 0)
                filtered[psms[1]] += 1
            elif isFiltered0 == 1 and isFiltered1 == 0:
                del nonFiltered[psms[0]]
                psm_1_2_filtered.setdefault(methodFiltered0, {})
                psm_1_2_filtered[methodFiltered0].setdefault(psms[0], 0)
                psm_1_2_filtered[methodFiltered0][psms[0]] += 1
                filtered.setdefault(psms[0], 0)
                filtered[psms[0]] += 1
            elif isFiltered0 == 0 and isFiltered1 == 1:
                del nonFiltered[psms[1]]
                psm_1_2_filtered.setdefault(methodFiltered1, {})
                psm_1_2_filtered[methodFiltered1].setdefault(psms[1], 0)
                psm_1_2_filtered[methodFiltered1][psms[1]] += 1
                filtered.setdefault(psms[1], 0)
                filtered[psms[1]] += 1            
            else:
                sigma0 = np.std(np.log2(pepDict[psms[0]]['intensity']))
                sigma1 = np.std(np.log2(pepDict[psms[1]]['intensity']))
                if sigma1 > sigma0:
                    del nonFiltered[psms[1]]
                    psm_1_2_filtered.setdefault("postFilter", {})
                    psm_1_2_filtered["postFilter"].setdefault(psms[1], 0)
                    psm_1_2_filtered["postFilter"][psms[1]] += 1
                    filtered.setdefault(psms[1], 0)
                    filtered[psms[1]] += 1            
                elif sigma0 > sigma1:
                    del nonFiltered[psms[0]]
                    psm_1_2_filtered.setdefault("postFilter", {})
                    psm_1_2_filtered["postFilter"].setdefault(psms[0], 0)
                    psm_1_2_filtered["postFilter"][psms[0]] += 1
                    filtered.setdefault(psms[0], 0)
                    filtered[psms[0]] += 1
                else:
                    mu0 = np.mean(np.log2(pepDict[psms[0]]['intensity']))
                    mu1 = np.mean(np.log2(pepDict[psms[1]]['intensity']))
                    if mu0 > mu1:
                        del nonFiltered[psms[1]]
                        psm_1_2_filtered.setdefault("postFilter", {})
                        psm_1_2_filtered["postFilter"].setdefault(psms[1], 0)
                        psm_1_2_filtered["postFilter"][psms[1]] += 1
                        filtered.setdefault(psms[1], 0)
                        filtered[psms[1]] += 1
                    else:
                        del nonFiltered[psms[0]]
                        psm_1_2_filtered.setdefault("postFilter", {})
                        psm_1_2_filtered["postFilter"].setdefault(psms[0], 0)
                        psm_1_2_filtered["postFilter"][psms[0]] += 1
                        filtered.setdefault(psms[0], 0)
                        filtered[psms[0]] += 1    
    return (psm_1_2_filtered, filtered, nonFiltered)

def generateRawTxtFiles(outDict, pepDict, params):
    reporters = params['tmt_reporters_used'].split(';')
    reporterMzs = getReporterMz(reporters)
    nReporters = len(reporters)

    # 1st filtering of PSMs by
    # (1) Low PPI, (2) Zero intensity and (3) "min_intensity_method/value" parameter
    filterMethods = list(map(int, params['min_intensity_method'].split(',')))
    filterCutoffs = list(map(int, params['min_intensity_value'].split(',')))
    filterNames = []
    for i in range(0, len(filterMethods)):
        if filterMethods[i] == 1:
            filterNames.append("minimum")
        elif filterMethods[i] == 2:
            filterNames.append("maximum")
        elif filterMethods[i] == 3:
            filterNames.append("mean")
        elif filterMethods[i] == 4:
            filterNames.append("median")
        elif filterMethods[i] == 0:
            filterNames.append("no filter")
    filtered = {}
    nonFiltered = {}
    zeroFiltered = {}
    intFiltered = {}
    for name in filterNames:
        intFiltered[name] = {}
    filtered, nonFiltered, zeroFiltered, intFiltered = psmFilter1(outDict, pepDict, nReporters, filterMethods, filterCutoffs,
                                                                 filtered, nonFiltered, zeroFiltered, intFiltered)

    # 2nd filtering of PSMs by
    # "min_intensity_method/value_1_2_psm_filter" parameter
    # (filter or not proteins identified by 1 or 2 PSM(s))
    # It is only valid for whole proteome analysis
    filterMethods_1_2_psm = list(map(int, params['min_intensity_method_1_2_psm'].split(',')))
    filterCutoffs_1_2_psm = list(map(int, params['min_intensity_value_1_2_psm'].split(',')))
    filterNames_1_2_psm = []
    for i in range(0, len(filterMethods_1_2_psm)):
        if filterMethods_1_2_psm[i] == 1:
            filterNames_1_2_psm.append("minimum")
        elif filterMethods_1_2_psm[i] == 2:
                filterNames_1_2_psm.append("maximum")
        elif filterMethods_1_2_psm[i] == 3:
            filterNames_1_2_psm.append("mean")
        elif filterMethods_1_2_psm[i] == 4:
            filterNames_1_2_psm.append("median")
        elif filterMethods_1_2_psm[i] == 0:
            filterNames_1_2_psm.append("no filter")    
    psm_1_2_filtered, filtered, nonFiltered = psmFilter2(outDict, pepDict, nReporters,
                                                         filterMethods_1_2_psm, filterCutoffs_1_2_psm,
                                                         filtered, nonFiltered)
    
    # Write to a file (raw_xxx_zero/nonzero.txt files)
    nonzeroName = params['save_dir'] + "/raw_" + params['save_dir'] + "_psm_nonzero.txt"
    zeroName = params['save_dir'] + "/raw_" + params['save_dir'] + "_psm_zero.txt"
    nonzeroFile = open(nonzeroName, 'w')
    zeroFile = open(zeroName, 'w')
    idTxtFile = params['idtxt']
    nLines = 0
    peps = {}
    prots = {}    
    with open(idTxtFile, "r") as IDFILE:
        for line in IDFILE:
            nLines += 1
            if nLines == 1:
                nonzeroFile.write(line)
                zeroFile.write(line)
            elif nLines == 2:
                line = line.rstrip("\n")
                nonzeroHeader = line + ";RT;K_y1;R_y1;charge;endingAA;nLabels"
                zeroHeader = line
                for i in range(0, nReporters):
                    mzName = ";mz" + reporters[i] + "(" + params[reporters[i]] + ")"
                    mzName = re.sub(r'sig', '', mzName)
                    nonzeroHeader += mzName
                    zeroHeader += mzName
                for i in range(0, nReporters):
                    intName = ";" + reporters[i] + "(" + params[reporters[i]] + ")"
                    nonzeroHeader += intName
                    zeroHeader += intName
                nonzeroHeader += "\n"
                zeroHeader += "\n"
                nonzeroFile.write(nonzeroHeader)
                zeroFile.write(zeroHeader)
            else:
                line = line.rstrip("\n")
                line = re.sub(r'//', '/', line)
                elems = line.split(";")
                pep, prot, psm = elems[0:3]
                _, pep, _ = pep.split('.')
                if psm not in pepDict.keys():
                    line += ";"
                    for i in range(0, nReporters):
                        line += "0;0;"
                    line += "\n"
                    zeroFile.write(line)
                elif psm in filtered.keys():
                    line += ";"
                    if psm in pepDict.keys():
                        line = line + ";".join(str(x) for x in pepDict[psm]['mz']) + \
                        ";" + ";".join(str(x) for x in pepDict[psm]['intensity']) + "\n"
                    else:
                        line = line + "0;" * 2 * nReporters + "\n"                    
                    zeroFile.write(line)
                else:
                    line = line + ";" + str(pepDict[psm]['rt']) + ";" + str(pepDict[psm]['Ky1']) + ";" + str(pepDict[psm]['Ry1']) + ";" + \
                    str(pepDict[psm]['charge']) + ";" + str(pepDict[psm]['endingAA']) + ";" + str(pepDict[psm]['nLabels']) + ";" + \
                    ";".join(str(x) for x in pepDict[psm]['mz']) + ";" + ";".join(str(x) for x in pepDict[psm]['intensity']) + "\n"
                    peps[pep] = 1
                    prots[prot] = 1
                    nonzeroFile.write(line)
    nonzeroFile.close()
    zeroFile.close()
    return (0)


########
# Main #
########
# Pre-defined variables/parameters
params = getParams('jump_q.params')
# reader = mzxml.read("HH_tmt10_human_jump.mzXML")
reporterMzs = (126.127726, 127.124761, 127.131081, 128.128116, 128.134436,
               129.131471, 129.137790, 130.134825, 130.141145, 131.138180)  # Make it as a tuple since it should not be changed
tol = 10

# Parsing of ID.txt file
outDict, pepDict = parseIdtxt("ID.txt", params)
print ("ID.txt file is parsed")

# Retrieve the information of MS2 spectra
ms2Dict = getMS2(outDict)
print ("MS2 spectra are extracted")

# Extract TMT reporter ion m/z and intensity values of PSMs
outDict, pepDict = getReporterIntensity(outDict, pepDict, ms2Dict, reporterMzs, params)
outDict, pepDict = refineReporterIntensity(outDict, pepDict, ms2Dict, reporterMzs, params)

# Impurity correction of TMT reporters
if params['impurity_correction'] == '1':
    print ("  Correction of the isotopic impurity of TMT reporter ions")
    pepDict = correctImpurity(pepDict, params)
else:
    print ("  No correction of the isotopic impurity of TMT reporter ions")
print()

## Generate raw_..._scan_zero/nonzero.txt files
print ("  Examining extracted PSMs (it may take a while)")
outDict, pepDict = generateRawTxtFiles(outDict, pepDict, params)

## Calculate loading-biases and normalize

