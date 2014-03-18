from  optparse  import OptionParser
import ROOT as rt
import RootTools
import array, os
import random

rt.gROOT.SetBatch(True);

def makebins(start_,end_,inc_,inc_inc_):
    bin = start_
    inc = inc_
    list = []
    while bin < end_:
        list.append(bin)
        bin+=inc
        inc*=(1+inc_inc_)
    return list

def get_xsec(xsec_file, msq, mgl):

    file = open(xsec_file,"r")
    lines = file.readlines()
    lines_stripped = map(lambda(x):x.rstrip("\n"),lines)

    #scan the lines
    for line in lines:
        if msq in line and mgl in line:
            #check the ordering is correct
            split = line.split()
            idx1 =  split.index(msq)
            idx2 =  split.index(mgl)
            if idx1 < idx2:
                #parse the values
                lo_xsec = float(split[6])
                lo_xsec_plus = float(split[8])
                lo_xsec_minus = float(split[10])

                return (lo_xsec,lo_xsec_plus,lo_xsec_minus)
            
    print "\n\nERROR:CROSS SECTION WAS NOT FOUND\n\n"
    return (-1,-1,-1)

parser = OptionParser()

parser.add_option("-f", "--file", dest="filename",
                  help="weight_hist.root file to analyze",
                  default="weight_hist.root",
                  action="store",type="string")

parser.add_option("-x", "--xsec", dest="xsec_file",
                 help="file containing signal cross sections",
                 default="Spectra_gsq_B_8TeV.xsec",
                 action="store",type="string")

parser.add_option("-s", "--signal", dest="sig_files",
                 help="text file containing a list of signal",
                 default="signal_files.txt",
                 action="store",type="string")

parser.add_option("-o", "--output_toys", dest="out_toys",
                 help="directory to output toys",
                 default="limit_toys_razor",
                 action="store",type="string")

parser.add_option("--mrmin", dest="mrmin",
                 help="minimum mr for the low rsq fit",
                 default=.6,
                 action="store",type="float")

(options, args) = parser.parse_args()

bins = makebins(options.mrmin, 5., .1, .3)

def make_data_card(signal_file, hist_exp, categories):
    signal_name = signal_file.split("/")[-1]
    msq = signal_name.split("_")[1]
    mgl = signal_name.split("_")[2]    

def draw_grid():
    pass

#################
##MAIN SCRIPT####
#################

#grab the file contianing th expectation histogram
in_file = rt.TFile(options.file)

#grab the file with lines of signal files for grid
sig_txt = open(options.sig_files,"READ")
sig_lines = map(lambda x:x.rstrip("\n"), sig_text.readlines())

for sig_point in sig_lines:
    print sig_point
    

