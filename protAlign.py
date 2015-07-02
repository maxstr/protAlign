#!/usr/bin/env python
from fabric.api import local, env, settings
from fabric.context_managers import hide
import numpy
import re
from itertools import chain, combinations
from os import listdir, makedirs
from os.path import isfile, join, normpath, exists, abspath
from sys import exit
import csv
import sys
import os
from pprint import pprint

# ./protAlign seqDirectory nativepdb [destination]

# from https://docs.python.org/2/library/itertools.html
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# Runs TM-align on two pdb files
def tmAlign(file1, file2):
    tmPath = join(os.path.dirname(os.path.realpath(__file__)), 'TMalign')
    with hide('running'):
        ret = local("%s %s %s" % (tmPath, file1, file2), capture=True)
    return ret

# Extract RMSD/TMExtract from TM-align output
def rmsdExtract(tmOutput):
    return re.search("RMSD= *([0-9.]*)", tmOutput).group(1)

def tmExtract(tmOutput):
    return re.findall("TM-score= ([0-9.]*)", tmOutput)

def main(seqDirectory, nativePDB, destination):
    # First make sure data directory exists and is writeable
    if not exists(join(destination, 'data')):
        try:
            makedirs(join(destination, 'data'))
        except:
            exit("Couldn't create destination directories and they don't exist!")

    results = {}
    results['seq'] = seqDirectory
    dataDir = join(destination, 'data')
    results['dataDir'] = abspath(dataDir)
    # First we handle all model-model comparisons
    # Generates a list of all model[0-9].pdb in seqDirectory
    models = [f for f in listdir(seqDirectory) if re.search("model[0-9]*.pdb", f) and isfile(join(seqDirectory, f))]
    results['numModels'] = len(models)
    # Generates all possible pairs of models
    modelPairs = [f for f in powerset(models) if len(f) == 2]
    # Makes a dictionary with (modelx, modely)
    modelsResultsDict = {}
    for pair in modelPairs:
        model1, model2 = pair
        tmOutput = tmAlign(join(seqDirectory, model1), join(seqDirectory, model2))
        rmsd = rmsdExtract(tmOutput)
        tmScores = tmExtract(tmOutput)
        modelsResultsDict[pair] = { 'rmsd': rmsd, 'tmScore-1': tmScores[0], 'tmScore-2':tmScores[1], 'pair':f }
        with open(abspath(join(dataDir, "%s-%s-align.txt" % (model1.replace(".pdb", ""), model2.replace(".pdb", "")))), 'w') as f:
            f.write(tmOutput + "\n")

    # Now we need to calculate STD Deviation, Mean for these values.
    # First get list of all RMSDs, TMScores
    modelsRMSD = [ float(f[1]['rmsd']) for f in modelsResultsDict.items() ]
    # TM Scores are equivalent here, so we can use either.
    modelsTMScore = [ float(f[1]['tmScore-1']) for f in modelsResultsDict.items() ]
    # Then we calculate means, stddevs
    results['stdModelsRMSD'] = numpy.std(modelsRMSD)
    results['meanModelsRMSD'] = numpy.mean(modelsRMSD)
    results['stdModelsTMScore'] = numpy.std(modelsTMScore)
    results['meanModelsTMScore'] = numpy.mean(modelsTMScore)

    # Now we do model-native comparisons
    modelNativeResultsDict = {}
    for i in models:
        tmOutput = tmAlign(nativePDB, join(seqDirectory, i))
        rmsd = rmsdExtract(tmOutput)
        tmScores = tmExtract(tmOutput)
        modelNativeResultsDict[i] = { 'rmsd': rmsd, 'tmScore-1': tmScores[0], 'tmScore-2':tmScores[1], 'pair':i }
        with open(abspath(join(dataDir, "%s-%s-align.txt" % ('native', i.replace(".pdb", "")))), 'w') as f:
            f.write(tmOutput + "\n")

    modelNativeRMSD = [ float(f[1]['rmsd']) for f in modelNativeResultsDict.items() ]
    # Here we *must* use the first TM-Score, as it is the TM Score associated with using the native struct as the reference.
    modelNativeTMScore = [ float(f[1]['tmScore-1']) for f in modelNativeResultsDict.items() ]
    # Then we calculate means, stddevs
    results['stdModelNativeRMSD'] = numpy.std(modelNativeRMSD)
    results['meanModelNativeRMSD'] = numpy.mean(modelNativeRMSD)
    results['stdModelNativeTMScore'] = numpy.std(modelNativeTMScore)
    results['meanModelNativeTMScore'] = numpy.mean(modelNativeTMScore)
    # Now we find the model with the lowest RMSD, TMScore
    results['minRMSDModel'] = min(modelNativeResultsDict, key=lambda model:float(modelNativeResultsDict[model]['rmsd']))
    minRMSDModel = results['minRMSDModel']
    results['minRMSDValue'] = float(modelNativeResultsDict[minRMSDModel]['rmsd'])
    results['maxTMScoreModel'] = max(modelNativeResultsDict, key =lambda model:float(modelNativeResultsDict[model]['tmScore-1']))
    maxTMScoreModel = results['maxTMScoreModel']
    results['maxTMScoreValue'] = float(modelNativeResultsDict[maxTMScoreModel]['tmScore-1'])

    # We have all of our values, time to write the file.
    output = """\
Align Output for %(seq)s

Number of Models: %(numModels)i
Data Output: %(dataDir)s

Model-Model Comparisons

TMScore Mean: %(meanModelsTMScore)f
TMScore STD: %(stdModelsTMScore)f

RMSD Mean: %(meanModelsRMSD)f
RMSD STD: %(stdModelsRMSD)f

Native-Model Comparisons

TMScore Mean: %(meanModelNativeTMScore)f
TMScore STD: %(stdModelNativeTMScore)f
TMScore Max (model): %(maxTMScoreValue)f ( %(maxTMScoreModel)s )

RMSD Mean: %(meanModelNativeRMSD)f
RMSD STD: %(stdModelNativeRMSD)f
RMSD Min (model): %(minRMSDValue)f ( %(minRMSDModel)s )
"""

    with open(abspath(join(destination, "report.txt")), 'w') as report:
        report.write(output % results)




if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])







