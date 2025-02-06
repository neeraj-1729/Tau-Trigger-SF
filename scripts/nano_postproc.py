#!/usr/bin/env python3
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys; sys.path.append('python')
import ROOT
import argparse
import pickle

from summaryProducer import *
from selectionFilter import *
from tupleProducer import *

ROOT.PyConfig.IgnoreCommandLineOptions = True

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Post Processing.')
    parser.add_argument('--input', required=False, type=str, help="NANO input")
    parser.add_argument('--inputFileList', required=False, type=str, help="NANO input file list")
    parser.add_argument('--output', required=True, type=str, help="eventTuple output")
    parser.add_argument('--isMC', required=True, type=int, help="judge if isMC")
    parser.add_argument('--era', required=True, type=str, help="")
    args = parser.parse_args()
    print("args = ",args)

    if (args.input is None) and (args.inputFileList is None):
        raise RuntimeError("Please check the input!")
    if (args.input is not None) and (args.inputFileList is not None):
        raise RuntimeError("Please check the input!")

    isMC = args.isMC
    output = args.output
    era = args.era
    if args.input:
        files = [ args.input ]
    if args.inputFileList:
        f = open(args.inputFileList, "r")
        files = f.read().splitlines()

    if era == "2022EE": era = "2022"
    if era == "2023BPix": era = "2023"

    if isMC:
        jsoninput = None,
        if era == '2016':
            Modules = [summary2016MC(), selection2016MC(), tuple2016MC()]
        elif era == '2017':
            Modules = [summary2017MC(), selection2017MC(), tuple2017MC()]
        elif era == '2018':
            Modules = [summary2018MC(), selection2018MC(), tuple2018MC()]
        elif era == '2022':
            Modules = [summary2022MC(), selection2022MC(), tuple2022MC()]
        elif era == '2023':
            Modules = [summary2023MC(), selection2023MC(), tuple2023MC()]
        else:
            raise RuntimeError("Please check the right Year!")
        p = PostProcessor(output, files, "1", 
                        branchsel = "keep_and_drop.txt", 
                        modules= Modules, 
                        provenance=True,
                        outputbranchsel = "output_branch.txt"
        )

    else:
        if era == '2022':
            Modules = [summary2022data(), selection2022data(), tuple2022data()]
            lumi_json_path = './lumi_jsons/2022.txt' 
        elif era == '2023':
            Modules = [summary2023data(), selection2023data(), tuple2023data()]
            lumi_json_path = './lumi_jsons/2023.txt'
        else:
            raise RuntimeError("Please check the right Year!")

        with open('./lumi_jsons/'+era+'.pkl', 'rb') as file:
            runsAndLumis_special = pickle.load(file)
     
        jsoninput = runsAndLumis_special

        p = PostProcessor(output, files, "1", 
                        branchsel = "keep_and_drop.txt", 
                        modules= Modules, 
                        jsonInput=jsoninput,
                        provenance=True,
                        outputbranchsel = "output_branch.txt"
        )

    p.run()
    print("Done !")


