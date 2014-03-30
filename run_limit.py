from  optparse  import OptionParser
import ROOT as rt
import array, os
import random

rt.gROOT.SetBatch(True);

#container for limit categories
#categories are numbered 1,2,3,4...
class limit_categories:
    def __init__(self):
        self.ncategories = 0
        self.dict = {}

    def add_category(self, bins):
        self.ncategories += 1
        self.dict[self.ncategories] = bins    

    def get_event_sum_hist(self, hist, category):
        sum = 0
        for ii in self.dict[category]:
            sum+= hist.GetBinContent(ii)
        return sum        

#container for full grid of limit results
class grid_results:
    def __init__(self):
        self.grid = []
    
    def add_point(self,msq,mgl,result):
        if (msq,mgl) not in self.grid:
            self.grid.append((msq,mgl,result))

            
    def squark_masses(self):
        msq_masses = []
        for ii in self.grid: msq_masses.append(ii[0])
        return msq_masses
            
    def gluino_masses(self):
        gl_masses = []
        for ii in self.grid: gl_masses.append(ii[1])
        return gl_masses
    
    def print_summary(self):
        print "SUMMARY:"
        for point in self.grid:
            msq = point[0]
            mgl = point[1]
            res = point[2]
            
            print msq, mgl, "exp: %2.3f + %2.3f - %2.3f  obs: %2.3f " % (res.get_exp(), res.get_exp_p(), res.get_exp_m(), res.get_obs())                

    def build_band(self,hist):
        n_x = hist.GetNbinsX()
        n_y = hist.GetNbinsY()

        top_down_y = range(1,n_y+1)[::-1]

        limit_top_down = []
        done = False
        for yy in top_down_y:
            for xx in range(1,n_x+1):

                if done: continue

                bin = hist.GetBinContent(xx,yy)
                last_bin = hist.GetBinContent(xx-1,yy)

                last_bin_mass = hist.GetXaxis().GetBinCenter(xx-1)
                this_bin_mass = hist.GetXaxis().GetBinCenter(xx)                
                y_bin_mass = hist.GetYaxis().GetBinCenter(yy)                

                diff_mass = this_bin_mass - last_bin_mass

                if bin > 1.0:


                    abs_diff = bin - last_bin
                    diff_from_1 = 1.0 - last_bin
                    frac_diff = diff_from_1 / abs_diff

                    mass_limit_x = last_bin_mass + diff_mass * frac_diff

                    limit_top_down.append((mass_limit_x, y_bin_mass))                    
                    break
                
                elif xx == n_x: # we didnt find a bin we must extrapolate the limit

                    #extrapolate
                    linear_slope = (1 - bin) / (bin - last_bin)
                    mass_limit_x = this_bin_mass  + linear_slope * diff_mass

                    last_point = limit_top_down[-1]
                    last_point_x = last_point[0]
                    last_point_y = last_point[1]
                    
                    #extrapolate but cutt off at the last scanned point
                    corr_limit_y =  ((-2 * diff_mass ) / (mass_limit_x - last_point_x))*(this_bin_mass - last_point_x) + last_point_y
                    
                    #limit_top_down.append((this_bin_mass, corr_limit_y))

                    done = True
                    break
                
                elif xx == n_x - 1 and hist.GetBinContent(xx+1,yy) < 1: # we wont find a bin

                    linear_slope = (1 - bin) / (bin - last_bin)
                    mass_limit_x = this_bin_mass #+ linear_slope * diff_mass
                    
                    #limit_top_down.append((mass_limit_x, y_bin_mass))
                    #break


        return limit_top_down

    def join_two_bands_into_region(self, band1, band2):

        full_band = []

        for ii in reversed(band2): full_band.append(ii)

        for ii in band1: full_band.append(ii)

        full_band.append(band2[-1])

        return full_band

    def build_graph_from_band(self, band):
        x = []
        y = []

        sq_masses = grid_results.squark_masses()
        gl_masses = grid_results.gluino_masses()

        min_val = min(sq_masses + gl_masses)
        max_val = max(sq_masses + gl_masses)

        for ii in band:
            x.append(ii[0])
            y.append(ii[1])

        x_ar = array.array("d",x)
        y_ar = array.array("d",y)

        gr = rt.TGraph(len(x), x_ar, y_ar)

        gr.SetMinimum(min_val)
        gr.SetMaximum(max_val)

        return gr

    def average_2d_hist(self,hist):
        nx = hist.GetNbinsX()
        ny = hist.GetNbinsY()

        for xx in range(1, nx+1):
            for yy in range(1, ny+1):
                
                if hist.GetBinContent(xx,yy) == 0: #do averaging
                    sum = 0
                    n_neighbors = 0
                    
                    if xx != nx:
                        n_neighbors +=1
                        sum+= hist.GetBinContent(xx+1,yy) #not at the right edge
                    if xx != 1:
                        n_neighbors +=1
                        sum+= hist.GetBinContent(xx-1,yy) # not at the left edge
                    if yy != ny:
                        n_neighbors +=1
                        sum+= hist.GetBinContent(xx,yy+1) #not at the top
                    if yy != 1:
                        n_neighbors +=1
                        sum+=hist.GetBinContent(xx,yy-1) # not at the bottom

                    hist.SetBinContent(xx,yy, float(sum) / float(n_neighbors))

        return hist
    
    def build_hist_limits(self):
        sq_masses = grid_results.squark_masses()
        gl_masses = grid_results.gluino_masses()
        mres = 100
        mres2 = mres /2 

        min_x, max_x = min(sq_masses), max(sq_masses)
        min_y, max_y = min(gl_masses), max(gl_masses)

        n_bins_x = ((max_x - min_x) / mres) + 1
        n_bins_y = ((max_y - min_y) / mres) + 1
        
        exp_hist = rt.TH2F("exp_0", "expected limit", n_bins_x, min_x-mres2, max_x+mres2, n_bins_y, min_y-mres2, max_y+mres2)
        exp_p_hist = rt.TH2F("exp_p", "expected limit plus", n_bins_x, min_x-mres2, max_x+mres2, n_bins_y, min_y-mres2, max_y+mres2)
        exp_m_hist = rt.TH2F("exp_m", "expected limit minus", n_bins_x, min_x-mres2, max_x+mres2, n_bins_y, min_y-mres2, max_y+mres2)
        obs_hist = rt.TH2F("obs", "obs limit",  n_bins_x, min_x-mres2, max_x+mres2, n_bins_y, min_y-mres2, max_y+mres2)
        delta_hist = rt.TH2F("delta", "delta limit",  n_bins_x, min_x-mres2, max_x+mres2, n_bins_y, min_y-mres2, max_y+mres2)

        for point in self.grid:
            #prase the grid point
            msq = point[0]
            mgl = point[1]
            res = point[2]

            #parse the values from the limit
            exp = res.get_exp()
            exp_p = res.get_exp_p()
            exp_m = res.get_exp_m()
            obs = res.get_obs()

            delta = (exp - obs) / (exp_p - exp_m)

            #fill the values
            exp_hist.Fill(msq, mgl, exp)
            obs_hist.Fill(msq, mgl, obs)            
            exp_p_hist.Fill(msq, mgl, exp_p)
            exp_m_hist.Fill(msq, mgl, exp_m)
            delta_hist.Fill(msq, mgl, delta)


        #average the hists
        return map( lambda x: self.average_2d_hist(x) , [exp_hist, exp_p_hist, exp_m_hist, obs_hist, delta_hist])

#container for a single limit result
class limit_result:
    def __init__(self, msq, mgl, limit_tree):
        self.msq = msq
        self.mgl = mgl

        limit_tree.Draw(">>iterlist","","entrylist")
        itlist = rt.gDirectory.Get("iterlist")

        #loop over the tree and parse the limits
        for event in range(limit_tree.GetEntries()):            
            entry = itlist.Next()        
            limit_tree.GetEntry(entry)

            #observed
            if limit_tree.quantileExpected == -1:
                self.obs = limit_tree.limit
            elif abs(limit_tree.quantileExpected - .16) < .001: #lower quantile
                self.exp_m = limit_tree.limit
            elif abs(limit_tree.quantileExpected - .84) < .001: #expected
                self.exp_p = limit_tree.limit
            elif abs(limit_tree.quantileExpected - .5) < .001: #upper quantile
                self.exp = limit_tree.limit            

    def get_obs(self): return self.obs
    def get_exp(self): return self.exp
    def get_exp_p(self): return self.exp_p
    def get_exp_m(self): return self.exp_m



def makebins(start_,end_,inc_,inc_inc_):
        bin = start_
        inc = inc_
        list = []
        while True:
            list.append(bin)
            bin+=inc
            inc*= (1+inc_inc_)
            
            if bin > end_:
                list.append(bin)
                return list
            
def get_xsec(xsec_file, msq, mgl):
    msq = str(msq)
    mgl = str(mgl)

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

                nlo_xsec = float(split[12])
                nlo_xsec_plus = float(split[14])
                nlo_xsec_minus = float(split[16])

                return (nlo_xsec,nlo_xsec_plus,nlo_xsec_minus)
            
    print "\n\nERROR:CROSS SECTION WAS NOT FOUND\n\n"
    return (-1,-1,-1)

parser = OptionParser()


parser.add_option("-f", "--file", dest="filename",
                  help="weight_hist.root file to analyze",
                  default="weight_hist.root",
                  action="store",type="string")

parser.add_option("-v", "--verbose", dest="debug",
                  help="print out more information",
                  default=False,
                  action="store_true")

parser.add_option("--nocombine", dest="nocombine",
                  help="Don't run the combination. The files already exist.",
                  default=False,
                  action="store_true")

parser.add_option("-x", "--xsec", dest="xsec_file",
                 help="file containing signal cross sections",
                 default="Spectra_gsq_B_8TeV.xsec",
                 action="store",type="string")

parser.add_option("-s", "--signal", dest="sig_files",
                 help="text file containing a list of signal",
                 default="signal_files.txt",
                 action="store",type="string")

parser.add_option("-o", "--output_dir", dest="out_dir",
                 help="directory to output toys",
                 default="limit_toys_razor",
                 action="store",type="string")

parser.add_option("--mrmin", dest="mrmin",
                 help="minimum mr for the low rsq fit",
                 default=.6,
                 action="store",type="float")

parser.add_option("--rsq1", dest="rsq1",
                 help="minimum rsq for the low rsq fit",
                 default=.01,
                 action="store",type="float")

parser.add_option("--rsq2", dest="rsq2",
                 help="max rsq for the low rsq fit",
                 default=.02,
                 action="store",type="float")

parser.add_option("--lumi", dest="lumi",
                 help="Lumi to scale signal xsec to",
                 default=19.8,
                 action="store",type="float")



(options, args) = parser.parse_args()

bins = makebins(options.mrmin, 4.5, .1, .2)

def make_data_card(name,hist_exp, hist_low_exp, hist_obs, sig_pdf, xsec_error, categories):
    #parse category information
    n_cat = categories.ncategories
    cat_dict = categories.dict

    outfile = open(name,"w")
#    outfile.write("imax 1\n")
#    outfile.write("jmax %i\n" % n_cat)
#    outfile.write("kmax 3\n\n")

    ##################################################
    
    #enumerate the channels
    #NAMES FOR CHANNELS
    
    channel_string="bin\t\t"
    for ii in range(1,n_cat+1): channel_string+="bin%i\t" % ii
    outfile.write(channel_string+"\n")

    #parse the data expectations
    obs_string="observation\t"
    for ii in range(1,n_cat+1):
        sum = categories.get_event_sum_hist(hist_obs, ii)
        obs_string += "%i\t" % sum
    outfile.write(obs_string+"\n\n")
    
    ##################################################
    #EXPECATATIONS FOR SIGNAL AND BKG
    
    #1 signal and 1 background for each category
    bin_string = "bin\t"
    for ii in range(1,n_cat+1): bin_string += "%i\t %i\t" % (ii,ii)
    outfile.write(bin_string+"\n")

    #names for the processes
    process_string = "process\t"
    for ii in range(1,n_cat+1): process_string += "sig%i\tbkg%i\t" % (ii,ii) 
    outfile.write(process_string+"\n")

    #numerical labels for the separate processes
    process_string2 = "process\t"
    for ii in range(1,n_cat+1): process_string2 += "%i\t%i\t" % (ii*-1,ii)
    outfile.write(process_string2+"\n")

    #rate expected for each of the backgrounds / signals
    rate_string = "rate\t"
    for ii in range(1,n_cat+1):
        sum_bkg = categories.get_event_sum_hist(hist_exp, ii)
        sum_sig = categories.get_event_sum_hist(sig_pdf, ii)

        rate_string += "%2.2f\t%2.2f\t" % (sum_sig, sum_bkg)
    outfile.write(rate_string+"\n\n")

    ##################################################
    #SYSTEMATICS

    #LUMI lognormal
    lumi_string = "lumi\t\tlnN\t\t"
    for ii in range(1,n_cat+1): lumi_string+="1.026\t-\t"
    outfile.write(lumi_string+"\n")

    #SIGNAL CROSS SECTION lognormal
    xsec_string = "xs_sig\t\tlnN\t\t" 
    for ii in range(1,n_cat+1): xsec_string+="%2.2f\t-\t" % (1 + xsec_error)
    outfile.write(xsec_string+"\n")

    #BACKGROUND NORMALIZATION gamma
    #one line for each category"
    for ii in range(1,n_cat+1):
        sideband = int(categories.get_event_sum_hist(hist_low_exp, ii))
        extrap = categories.get_event_sum_hist(hist_exp,ii)

        scale_factor = 1 
        if sideband != 0:
            scale_factor =  float(extrap) / float(sideband)
        if sideband == 0:
            sideband = 1
            scale_factor = categories.get_event_sum_hist(hist_exp,ii)
            
        bkg_string = "bkg_norm_%i\tgmN %i\t\t" % (ii, sideband)
        
        for jj in range(1,n_cat+1):
            if ii == jj:
                bkg_string += "-\t%2.2f\t" % scale_factor
            else:
                bkg_string+= "-\t-\t"
        outfile.write(bkg_string+"\n")

def run_limit(data_card): pass

def parse_limit(root_file): pass

def draw_grid(): pass

#################
##MAIN SCRIPT####
#################

#grab the file contianing th expectation histogram
in_file = rt.TFile(options.filename)
exp_hist = in_file.Get("hist_high_pred")
exp_hist_low = in_file.Get("hist_low_pred")
obs_hist = in_file.Get("hist_high_data")

#grab the file with lines of signal files for grid
sig_txt = open(options.sig_files,"r")
sig_lines = map(lambda x:x.rstrip("\n"), sig_txt.readlines())

#declare the grid result
grid_results = grid_results()

#output directory
output_dir = options.out_dir

#define the categories
limit_cats = limit_categories()
for ii in range(5,11): limit_cats.add_category([ii])
#limit_cats.add_category([5,6,7,8,9,10])

output_file = rt.TFile(output_dir+"/limit_out.root","RECREATE")
#loop over each point
for sig_point in sig_lines:

    #parse the masses
    signal_name = sig_point.split("/")[-1]
    msq = int(signal_name.split("_")[1])
    mgl = int(signal_name.split("_")[2])

    #errors are assymetric, pick the larger of the two values for error
    (xsec, xsec_down, xsec_up) = get_xsec(options.xsec_file, msq, mgl)
    if options.debug: print "xsec, xsec_down, xsec_up", (xsec, xsec_down, xsec_up)
    xsec_error = min(xsec_down,xsec_up) / xsec 

    #build the signal pdf
    point_file = rt.TFile(sig_point)
    output_file.cd()
    
    tree = point_file.Get("HggOutput")
    bin_array = array.array("d", bins)

    signal_pdf  = rt.TH1F("signalpdf_%i_%i" % (msq,mgl) ,"Signal PDF", len(bins)-1, bin_array)
    tree.Draw("(PFMR/1000.)>>signalpdf_%i_%i"%(msq, mgl), "(PFMR/1000.) > %f && PFR^2 > %f && iSamp==0" % (options.mrmin, options.rsq1))

    #scale to the cross sections * luminosity * efficiency
    scale_factor= xsec * 1000.0 * (1. / 10000.) * options.lumi
    if options.debug:
        print "scale_factor:", scale_factor
        print "pre scaling integral:", signal_pdf.Integral()

    signal_pdf.Scale(scale_factor) #xsec * (femto/pico) * (1/ngen) * lumi
    signal_pdf.Write()
    
    if options.debug: print "SIGNAL NORM", signal_pdf.Integral()

    #build the datacard
    data_card_name = "datacard_msq_%s_mgl_%s.txt" % (msq,mgl)
    data_card_path = output_dir + "/" + data_card_name
    make_data_card(data_card_path, exp_hist, exp_hist_low, obs_hist, signal_pdf, xsec_error, limit_cats)

    #run the datacard
    fake_mass_name = str(msq)[:2]+str(mgl)[:2]
    name = "higgsCombineRA3.Asymptotic.mH%s.root" % fake_mass_name
    os.chdir(output_dir)

    if not options.nocombine:
        print "run combine.."
        os.system("combine -M Asymptotic %s -m %s -n RA3" % (data_card_name, fake_mass_name))


    #parse the limits
    limit_file = rt.TFile(name)
    limit_tree = limit_file.Get("limit")
    result = limit_result(msq, mgl, limit_tree)

    #add it to the full grid of results
    grid_results.add_point(msq, mgl, result)

    #move back to the home directory
    os.chdir("..")



grid_results.print_summary()    

(exp_hist, exp_p_hist, exp_m_hist, obs_hist, delta_hist) = grid_results.build_hist_limits()
limit_hists = [exp_hist, exp_p_hist, exp_m_hist, obs_hist, delta_hist]

output_file.cd()

#solid line for the observation
band_obs =  grid_results.build_band(obs_hist)
graph_obs = grid_results.build_graph_from_band(band_obs)
graph_obs.SetLineWidth(6)
graph_obs.Write("obs_graph")

band_exp = grid_results.build_band(exp_hist)
graph_exp = grid_results.build_graph_from_band(band_exp)
graph_exp.SetLineStyle(7)
graph_exp.SetLineColor(rt.kRed)
graph_exp.SetLineWidth(6)
graph_exp.Write("exp_0_graph")

#lines for the expected region
band_exp_p = grid_results.build_band(exp_p_hist)
graph_exp_p = grid_results.build_graph_from_band(band_exp_p)
graph_exp_p.SetLineStyle(7)
graph_exp_p.SetLineWidth(3)
graph_exp_p.SetLineColor(rt.kAzure-1)
graph_exp_p.Write("exp_p_graph")

band_exp_m = grid_results.build_band(exp_m_hist)
graph_exp_m = grid_results.build_graph_from_band(band_exp_m)
graph_exp_m.SetLineStyle(7)
graph_exp_m.SetLineWidth(3)
graph_exp_m.SetLineColor(rt.kAzure-1)
graph_exp_m.Write("exp_m_graph")

#fill ther region for the expected
expected_band = grid_results.join_two_bands_into_region(band_exp_m, band_exp_p)
graph_exp_band = grid_results.build_graph_from_band(expected_band)
graph_exp_band.SetFillColor(rt.kBlue-10)
graph_exp_band.SetLineWidth(3)
graph_exp_band.SetLineStyle(7)
graph_exp_band.SetFillStyle(3001)
graph_exp_band.SetLineColor(rt.kAzure-1)
graph_exp_band.Write("exp_graph")


#graphs = [graph_obs, graph_exp, graph_exp_band]
#multi_graph = rt.TMultiGraph("mg", "mg")

for ii in limit_hists: ii.Write()
    
output_file.Close()
