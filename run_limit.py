from  optparse  import OptionParser
import ROOT as rt
import array, os
import random

rt.gROOT.SetBatch(True);

class grid_result():
    def __init__(self):
        self.grid = []
    
    def add_point(self,msq,mgl,result):
        if (msq,mgl) not in self.grid:
            self.grid.append((msq,mgl,result))

class limit_result():
    def __init__(self, msq, mgl):
        self.exp, self.obs, self.msq, self.mgl = -1, -1, msq, mgl

    def set_obs(self, obs): self.obs = obs
    def set_exp(self, exp): self.exp = exp
    def set_exp_p(self,exp_p): self.exp_p = exp_p
    def set_exp_m(self,exp_m): self.exp_m = exp_m

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

def make_data_card(name, file, hist_exp, hist_obs, n_cat, cat_bins):
    outfile = open(name,"w")

    outfile.write("imax 1\n")
    outfile.write("jmax %i\n" % n_cat)
    outfile.write("kmax 3\n")

    #enumerate the channels
    channel_string="bin\t"
    for ii in range(n_cat): channel_string+="%i\t" % ii
    outfile.write(channel_string+"\n")

    #parse the data expectations
    obs_string="observation\t"
    for ii in cat_bins: obs_string += "%i\t" % hist_obs.GetBinContent(ii)
    outfile.write(obs_string+"\n")

    #1 signal and 1 background for each category
    bin_string = "bin\t"
    for ii in range(n_cat): bin_string += "%i\t %i\t" % (ii,ii)
    outfile.write(bin_string+"\n")
        
    process_string = "process\t"
    for ii in range(n_cat): process_string += "sig%i\tbkg%i\t" % (ii,ii) 
    outfile.write(process_string+"\n")
    
    process_string2 = "process\t"
    for ii in range(n_cat+1): process_string2 += "%i\t" % ii 
    outfile.write(process_string2+"\n")

    rate_string

def run_limit(data_card): pass

def parse_limit(root_file): pass

def draw_grid(): pass

#################
##MAIN SCRIPT####
#################

#grab the file contianing th expectation histogram
in_file = rt.TFile(options.filename)
exp_hist = in_file.Get("hist_high_pred")

#grab the file with lines of signal files for grid
sig_txt = open(options.sig_files,"r")
sig_lines = map(lambda x:x.rstrip("\n"), sig_txt.readlines())

#declare the grid result
grid_result = grid_result()

#loop over each point
for sig_point in sig_lines:
    point_file = rt.TFile(sig_point)

    #parse the masses
    signal_name = signal_point.split("/")[-1]
    msq = signal_name.split("_")[1]
    mgl = signal_name.split("_")[2]    

    #build the datacard
    data_card_name = "datacard_msq%s_mgl%s.txt" % (msq,mgl)
    make_data_card(name, point_file, exp_hist)
    


