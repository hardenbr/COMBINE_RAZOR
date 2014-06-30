from  optparse  import OptionParser
import ROOT as rt
import array, os
import random
import math
rt.gROOT.SetBatch(True);


#container for full grid of limit results
class grid_results:
    def __init__(self, t5gg):
        self.grid = []
        self.t5gg = t5gg
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

        top_down_y = range(1,n_y+1)[::-1] #reverse the range 

        limit_top_down = []
        done = False
        for yy in top_down_y:
            is_first = True

            for xx in xrange(1,n_x+1):
                #don't consider the beyond m_ne
                if self.t5gg:
                    bin_mass_gl = hist.GetXaxis().GetBinCenter(xx)
                    bin_mass_neu = hist.GetYaxis().GetBinCenter(yy)

                    if bin_mass_neu > bin_mass_gl: continue

                if done: continue

                bin = hist.GetBinContent(xx,yy)
                last_bin = hist.GetBinContent(xx-1,yy)
                next_bin = hist.GetBinContent(xx+1,yy)                     

                last_bin_mass = hist.GetXaxis().GetBinCenter(xx-1)
                this_bin_mass = hist.GetXaxis().GetBinCenter(xx)
                next_bin_mass = hist.GetXaxis().GetBinCenter(xx+1)                

                y_bin_mass = hist.GetYaxis().GetBinCenter(yy)                

                diff_mass = this_bin_mass - last_bin_mass

                if bin > 1.0:

                    if is_first:
                        limit_diff = next_bin - bin #limit difference
                        mass_diff = next_bin_mass - this_bin_mass #mass difference

                        mass_limit_x = ((1 - bin) * mass_diff) / limit_diff + this_bin_mass

                        limit_top_down.append((mass_limit_x, y_bin_mass))

                    else:
                        abs_diff = bin - last_bin
                        diff_from_1 = 1.0 - last_bin
                        frac_diff = diff_from_1 / abs_diff
                        
                        mass_limit_x = last_bin_mass + diff_mass * frac_diff

                        limit_top_down.append((mass_limit_x, y_bin_mass))
                    break
                
                elif xx == n_x: # we didnt find a bin we must extrapolate the limit
                    #one = True                
                    #print "SOMETHING IS WRONG. BINVAL =", bin
                    #break
                                    
                    #extrapolate
                    limit_diff = bin - last_bin #limit difference
                    mass_diff =  this_bin_mass - last_bin_mass #mass difference
                    
                    mass_limit_x = ((1 - bin) * mass_diff) / limit_diff + this_bin_mass

                    limit_top_down.append((mass_limit_x, y_bin_mass))
                    
                    #limit_top_down.append((this_bin_mass, corr_limit_y))

                    break
                
                elif xx == n_x - 1 and hist.GetBinContent(xx+1,yy) < 1: # we wont find a bin

                    linear_slope = (1 - bin) / (bin - last_bin)
                    mass_limit_x = this_bin_mass #+ linear_slope * diff_mass
                    
                    #limit_top_down.append((mass_limit_x, y_bin_mass))
                    #break

                is_first = False

            is_first = True

        print "band_result =", limit_top_down
        
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

        print "x", x
        print "y", y
        gr = rt.TGraph(len(x), x_ar, y_ar)

        gr.SetMinimum(min_val)
        gr.SetMaximum(max_val)

        return gr

    def average_2d_hist(self,hist):
        nx = hist.GetNbinsX()
        ny = hist.GetNbinsY()

        for xx in xrange(1, nx+1):
            for yy in xrange(1, ny+1):

                #don't consider m_neu > m_gl
                if self.t5gg:
                    bin_mass_gl = hist.GetXaxis().GetBinCenter(xx)
                    bin_mass_neu = hist.GetYaxis().GetBinCenter(yy)

                    if bin_mass_neu > bin_mass_gl: continue
                
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
        mres_x = 100
        mres_y = 100

        #smaller binning for t5gg
        if self.t5gg:
            mres_x = 50
            mres_y = 50
        
        mres2_x = mres_x / 2
        mres2_y = mres_y / 2 
        
        min_x, max_x = min(sq_masses), max(sq_masses)
        min_y, max_y = min(gl_masses), max(gl_masses)

        n_bins_x = ((max_x - min_x) / mres_x) + 1
        n_bins_y = ((max_y - min_y) / mres_y) + 1
        
        exp_hist = rt.TH2F("exp_0", "expected limit", n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        exp_theory = rt.TH2F("exp_th", "limit xsec excluded", n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        exp_p_hist = rt.TH2F("exp_p", "expected limit plus", n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        exp_m_hist = rt.TH2F("exp_m", "expected limit minus", n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        obs_hist = rt.TH2F("obs", "obs limit",  n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        obs_p_hist = rt.TH2F("obs_p", "obs limit plus sigma",  n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        obs_m_hist = rt.TH2F("obs_m", "obs limit minus sigma",  n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)
        delta_hist = rt.TH2F("delta", "delta limit",  n_bins_x, min_x-mres2_x, max_x+mres2_x, n_bins_y, min_y-mres2_y, max_y+mres2_y)        

        for point in self.grid:
            #prase the grid point
            msq = point[0]
            mgl = point[1]
            res = point[2]

            (xsec,xsec_down,xsec_up) = (-1, -1, -1)
            (rate_error, acc_error) = (None, None)

            if options.t5gg:
                (rate_error, acc_error) = (0, 0)
            else:
                (rate_error, acc_error) = get_pdf_errors(msq, mgl)

            if options.t5gg:
                (xsec, xsec_down, xsec_up) = get_xsec_t5gg(options.xsec_file, msq)
            else:
                (xsec, xsec_down, xsec_up) = get_xsec(options.xsec_file, msq, mgl)
                        
            #parse the values from the limit
            exp = res.get_exp()
            exp_p = res.get_exp_p()
            exp_m = res.get_exp_m()
            obs = res.get_obs()

            xsec_error_up = math.sqrt(xsec_up*xsec_up + (xsec * rate_error)*(xsec * rate_error) + (xsec * acc_error)*(xsec * acc_error))
            xsec_error_down = math.sqrt(xsec_down*xsec_down + (xsec * rate_error)*(xsec * rate_error) + (xsec * acc_error)*(xsec * acc_error))

            #convert the signal strength back to a cross section
            #reconvert to signal strength with new xsec 
            obs_p = (obs * xsec) / (xsec + xsec_error_up)
            obs_m = (obs * xsec) / (xsec - xsec_error_down)
            exp_th = exp * xsec * 1000

            delta = (exp - obs) / (exp_p - exp_m)

            #fill the values
            exp_hist.Fill(msq, mgl, exp)
            exp_theory.Fill(msq, mgl, exp_th)
            obs_hist.Fill(msq, mgl, obs)
            obs_p_hist.Fill(msq, mgl, obs_p)
            obs_m_hist.Fill(msq, mgl, obs_m)            
            exp_p_hist.Fill(msq, mgl, exp_p)
            exp_m_hist.Fill(msq, mgl, exp_m)
            delta_hist.Fill(msq, mgl, delta)

        #average the hists and return them
#        return [exp_hist, exp_p_hist, exp_m_hist, obs_hist, obs_p_hist, obs_m_hist, exp_theory, delta_hist]
        return map( lambda x: self.average_2d_hist(x) , [exp_hist, exp_p_hist, exp_m_hist, obs_hist, obs_p_hist, obs_m_hist, exp_theory, delta_hist])

#container for a single limit result
class limit_result:
    def __init__(self, msq, mgl, limit_tree):
        self.msq = msq
        self.mgl = mgl

        limit_tree.Draw(">>iterlist","","entrylist")
        itlist = rt.gDirectory.Get("iterlist")

        #loop over the tree and parse the limits
        for event in xrange(limit_tree.GetEntries()):            
            entry = itlist.Next()        
            limit_tree.GetEntry(entry)

            #observed
            if limit_tree.quantileExpected == -1:
                self.obs = limit_tree.limit
            elif abs(limit_tree.quantileExpected - .16) < .001: #lower quantile
                self.exp_m = limit_tree.limit
            elif abs(limit_tree.quantileExpected - .84) < .001: #upper quantile
                self.exp_p = limit_tree.limit
            elif abs(limit_tree.quantileExpected - .5) < .001: #expected
                self.exp = limit_tree.limit            

    def get_obs(self): return self.obs
    def get_exp(self): return self.exp
    def get_exp_p(self): return self.exp_p
    def get_exp_m(self): return self.exp_m

class toy_result:
    def __init__(self, msq, mgl, exp, exp_m, exp_p, obs):
        self.msq = msq
        self.mgl = mgl

        exp.Draw(">>iterlist_exp","","entrylist")
        itlist_exp = rt.gDirectory.Get("iterlist_exp")

        exp_p.Draw(">>iterlist_exp_p","","entrylist")
        itlist_exp_p = rt.gDirectory.Get("iterlist_exp_p")

        exp_m.Draw(">>iterlist_exp_m","","entrylist")
        itlist_exp_m = rt.gDirectory.Get("iterlist_exp_m")

        obs.Draw(">>iterlist_obs","","entrylist")
        itlist_obs = rt.gDirectory.Get("iterlist_obs")

        #expected
        entry = itlist_exp.Next()        
        exp.GetEntry(entry)
        self.exp = exp.limit 

        #+1 sigma
        entry = itlist_exp_p.Next()        
        exp_p.GetEntry(entry)
        self.exp_p = exp_p.limit

        #-1 sigma
        entry = itlist_exp_m.Next()        
        exp_m.GetEntry(entry)
        self.exp_m = exp_m.limit 

        #observed
        entry = itlist_obs.Next()        
        obs.GetEntry(entry)
        self.obs = obs.limit 

    def get_obs(self): return self.obs
    def get_exp(self): return self.exp
    def get_exp_p(self): return self.exp_p
    def get_exp_m(self): return self.exp_m

    def print_summary(self):
        print "obs: %s" % self.obs
        print "exp_0: %s" % self.exp
        print "exp_p: %s" % self.exp_p
        print "exp_m: %s" % self.exp_m

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

        if msq == line.split()[1] and mgl in line.split()[2]:
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

def get_xsec_t5gg(xsec_file,mgl):
    mgl = str(mgl)

    file = open(xsec_file,"r")
    lines = file.readlines()
    lines_stripped = map(lambda(x):x.rstrip("\n"),lines)

    #scan the lines
    for line in lines:
        if mgl in line.split()[0]:
            #check the ordering is correct
            split = line.split()
            print split
            xsec = float(split[1])
            xsec_plus = (float(split[2])) * xsec
            xsec_minus = (float(split[2])) * xsec

            print "\t\t", mgl, xsec

            return (xsec, xsec_plus, xsec_minus)
            
    print "\n\nERROR:CROSS SECTION WAS NOT FOUND\n\n"
    return (-1,-1,-1)


def build_systematic(hist, bins):
    #if error_hist.GetNbinsX() != hist.GetNbinsX():
    #    print "error: closure systemtic hist bins do not match!"
    #    exit(1)

    error_up = rt.TH1F("closure_up","closure_up", len(bins)-1,array.array("d",bins))
    error_down = rt.TH1F("closure_down","closure_down", len(bins)-1,array.array("d",bins))
    

    for bin in range(hist.GetNbinsX()+1):

        error_up_val = .05
        error_down_val = .05

        if hist.GetBinCenter(bin) > 1.3:

            error_up_val = .2
            error_down_val = .50

        print "--using closure errors:", (error_up_val, error_down_val)
            
        val = hist.GetBinContent(bin)

        val_up = (1 + error_up_val) * val
        val_down = (1 - error_down_val) * val

        error_up.SetBinContent(bin, val_up)
        error_down.SetBinContent(bin, val_down)
    
    return (error_up, error_down)

parser = OptionParser()

parser.add_option("--t5gg", dest="t5gg",
                 help="flag to run t5gg limit",
                 action="store_true",default=False)

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
                 default="signal_files_jes.txt",
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
                 default=19.7,
                 action="store",type="float")

parser.add_option("--pdf", dest="pdf",
                 help="text file containing xsec and acceptance uncertainties",
                 default="Spectra_gsq_B_dipho_envpdfuncert.dat",
                 action="store",type="string")

parser.add_option("--fake", dest="fake_res",
                 help="the weight_hist.root file from the fake study. This will produce the closure errror",
                 default="no_file",
                 action="store",type="string")

parser.add_option("--nolimit", dest="no_limit",
                 help="Just run the toys, don't generate the limit",
                 default=False,
                 action="store_true")


(options, args) = parser.parse_args()

bins = [0.6, 0.7, 0.82, 0.964, 1.1368, 1.34416, 1.592992, 1.8915904, 2.2499084799999998,  4.558050224127999]
#bins = makebins(options.mrmin, 4.5, .1, .2)

def make_shape_file(dir, hist_exp, hist_exp_up, hist_exp_down, hist_obs, sig_pdf, sig_pdf_up, sig_pdf_down, jes_up, jes_down, closure_up, closure_down, msq, mgl):
    outfile = rt.TFile(dir+"/"+"shapes_%i_%i.root" % (msq, mgl), "RECREATE") 
    outfile.cd()

    hist_exp.Write("bkg")
    hist_exp_up.Write("bkg_normUp")
    hist_exp_down.Write("bkg_normDown")
    sig_pdf.Write("sig")
    hist_obs.Write("data_obs")

    sig_pdf_up.Write("sig_jesUp")
    sig_pdf_down.Write("sig_jesDown")

    jes_up.Write("bkg_jesUp")
    jes_down.Write("bkg_jesDown")

    closure_up.Write("bkg_closeUp")
    closure_down.Write("bkg_closeDown")

    outfile.Close()
    
def make_data_card(output_dir, name, hist_exp, hist_exp_up, hist_exp_down, hist_obs, sig_pdf,sig_pdf_up, sig_pdf_down, jes_up, jes_down, closure_up, closure_down, msq, mgl):

    filename =  output_dir + "/" + name

    #first build the shape file shapes_msq_mgl
    make_shape_file(output_dir, hist_exp, hist_exp_up, hist_exp_down, hist_obs, sig_pdf, sig_pdf_up, sig_pdf_down, jes_up, jes_down, closure_up, closure_down, msq, mgl)
    
    #parse category information
    outfile = open(filename,"w")
    ##################################################
    #SHAPE INFORMATION
    shape_string = "shapes * * shapes_%i_%i.root $PROCESS $PROCESS_$SYSTEMATIC" % (msq,mgl)
    outfile.write(shape_string+"\n\n")
                  
    ##################################################    
    #enumerate the channels
    #NAMES FOR CHANNELS
    channel_string = "bin \t\t gg"
    outfile.write(channel_string+"\n")

    #parse the data expectations
    obs_string="observation\t %i" % hist_obs.Integral()
    outfile.write(obs_string+"\n\n")
    
    ##################################################
    #EXPECATATIONS FOR SIGNAL AND BKG
    
    #1 signal and 1 background for each category
    bin_string = "bin\t\tgg\tgg"
    outfile.write(bin_string+"\n")

    #names for the processes
    process_string = "process\t\t1\t-1"
    outfile.write(process_string+"\n")

    #numerical labels for the separate processes
    process_string2 = "process\t\tbkg\tsig"
    outfile.write(process_string2+"\n")

    #rate expected for each of the backgrounds / signals
    nbkg = hist_exp.Integral()
    nobs = hist_obs.Integral()
    nsig = sig_pdf.Integral()
    rate_string = "rate\t\t%2.2f\t%2.2f" % (nbkg, nsig)
    outfile.write(rate_string+"\n\n")

    ##################################################
    #SYSTEMATICS

    #parse the pdf errors for the mass point
    #    (rate_error, acc_error) = get_pdf_errors(msq, mgl)

    #LUMI lognormal
    lumi_string = "lumi\t\tlnN\t\t-\t1.026 "
    outfile.write(lumi_string+"\n")

    #PDF - RATE ERROR (THIS IS INCLUDED IN THE OBSERVED BAND)
    pdf_rate_string = "rate_pdf\tlnN\t-\t%2.3f" % (1 + 0)
    #outfile.write(pdf_rate_string+"\n")
     
    #PDF - ACCEPT ERROR  (THIS IS INCLUDEDI N THE OBSERVED BAND)
    pdf_acc_string = "acc_pdf\tlnN\t-\t%2.3f" % (1 + 0)
    #outfile.write(pdf_acc_string+"\n")
    
    #BACKGROUND NORMALIZATION 
    bkg_string = "norm shape\t1\t-"         
    outfile.write(bkg_string+"\n")

    #JES BKG SHAPE
    jes_bkg_string = "jes shape\t1\t-"         
    #outfile.write(jes_bkg_string+"\n")

    #JES SIGNAL SHAPE
    jes_sig_string = "jes shape\t-\t1"         
    outfile.write(jes_sig_string+"\n")

    #BKG CLOSURE TEST SYSTEMATIC
    close_string = "close shape\t1\t-"         
    outfile.write(close_string+"\n")
    
#from ra3privatesignalmc TWIKI
#pdf uncert file format is the following: printf("%i %i %i %i %i %f %f\n",nevents,mgluino,msquark,mbino,mwino,xsecpdferrs,acceppdferrs);
def get_pdf_errors(msq, mgl):
    pdf_file = open(options.pdf)
    lines = map(lambda x: x.rstrip("\n"), pdf_file.readlines())

    for line in lines:
        split = line.split()

        mgl_line = int(split[1])
        msq_line = int(split[2])
        
        rate_error = float(split[5])
        acc_error = float(split[6])
        
        if msq == msq_line and mgl == mgl_line:
            return (rate_error/100., acc_error/100.)


    print "ERROR: NO MATCHING PDF UNCERTAINTY"
    exit(1)
        
    
#################
##MAIN SCRIPT####
#################

#grab the file contianing th expectation histogram
in_file = rt.TFile(options.filename)
exp_up = in_file.Get("hist_high_shape_up")
exp_down = in_file.Get("hist_high_shape_down")
exp_hist = in_file.Get("hist_high_pred")
exp_hist_low = in_file.Get("hist_low_pred")
obs_hist = in_file.Get("hist_high_data")

jes_high_up = in_file.Get("hist_high_data_up")
jes_high_down = in_file.Get("hist_high_data_down")

jes_up = in_file.Get("hist_low_data_up")
jes_up.Scale(float(jes_high_up.Integral()) / float(jes_up.Integral()))

jes_down = in_file.Get("hist_low_data_down")
jes_down.Scale(float(jes_high_down.Integral()) / float(jes_down.Integral()))

#get the fake histogram estiamtes for the closure systematic
#fake_in_file = rt.TFile(options.fake_res)
#fake_exp_hist = fake_in_file.Get("hist_high_pred")
#fake_exp_hist_copy = fake_in_file.Get("hist_high")
#fake_obs_hist = fake_in_file.Get("hist_high_data")

#build the error (exp - obs) / exp
#fake_exp_hist.Add(fake_obs_hist, -1)
#fake_exp_hist.Divide(fake_exp_hist_copy)

#build the systematic on the 
(closure_up, closure_down) = build_systematic(exp_hist, bins)

#grab the file with lines of signal files for grid
sig_txt = open(options.sig_files,"r")
sig_lines = map(lambda x:x.rstrip("\n"), sig_txt.readlines())

#declare the grid result
grid_results = grid_results(options.t5gg)

#output directory
output_dir = options.out_dir

#limit_cats.add_category([5,6,7,8,9,10])

output_file = rt.TFile(output_dir+"/limit_out.root","RECREATE")
#loop over each point
for sig_point in sig_lines:

    #parse the masses
    signal_name = sig_point.split("/")[-1]

    msq = int(signal_name.split("_")[1])
    mgl = int(signal_name.split("_")[2])
            
    #errors are assymetric, pick the larger of the two values for error
    if not options.t5gg:
        (xsec, xsec_down, xsec_up) = get_xsec(options.xsec_file, msq, mgl)
    else:
        (xsec, xsec_down, xsec_up) = get_xsec_t5gg(options.xsec_file, msq) #msq is mgl for t5gg mgl-->mneutralino

    if options.debug: print "xsec, xsec_down, xsec_up", (xsec, xsec_down, xsec_up)

    #build the signal pdf
    point_file = rt.TFile(sig_point)
    output_file.cd()
    
    tree = point_file.Get("HggOutput")
    bin_array = array.array("d", bins)

    signal_pdf  = rt.TH1F("signalpdf_%i_%i" % (msq,mgl) ,"Signal PDF", len(bins)-1, bin_array)
    signal_pdf_up  = rt.TH1F("signalpdf_%i_%i_up" % (msq,mgl) ,"Signal PDF", len(bins)-1, bin_array)
    signal_pdf_down  = rt.TH1F("signalpdf_%i_%i_down" % (msq,mgl) ,"Signal PDF", len(bins)-1, bin_array)

    tree.Draw("(PFMR/1000.)>>signalpdf_%i_%i"%(msq, mgl), "(PFMR/1000.) > %f && PFR^2 > %f && iSamp==0" % (options.mrmin, options.rsq2))
    tree.Draw("(PFMR_UP/1000.)>>signalpdf_%i_%i_up"%(msq, mgl), "(PFMR_UP/1000.) > %f && PFR_UP^2 > %f && iSamp==0" % (options.mrmin, options.rsq2))
    tree.Draw("(PFMR_DOWN/1000.)>>signalpdf_%i_%i_down"%(msq, mgl), "(PFMR_DOWN/1000.) > %f && PFR_DOWN^2 > %f && iSamp==0" % (options.mrmin, options.rsq2))

    #scale to the cross sections * luminosity * efficiency
    scale_factor= xsec * 1000.0 * (1. / 10000.) * options.lumi
    if options.debug:
        print "scale_factor:", scale_factor
        print "pre scaling integral:", signal_pdf.Integral()

    signal_pdf.Scale(scale_factor) #xsec * (femto/pico) * (1/ngen) * lumi
    signal_pdf_up.Scale(scale_factor)
    signal_pdf_down.Scale(scale_factor)
    signal_pdf.Write()
    
    if options.debug: print "SIGNAL NORM", signal_pdf.Integral()

    #build the datacard and the file containing shapes
    data_card_name = "datacard_msq_%s_mgl_%s.txt" % (msq,mgl)
    make_data_card(output_dir, data_card_name, exp_hist, exp_up, exp_down, obs_hist, signal_pdf,signal_pdf_up, signal_pdf_down, jes_up, jes_down, closure_up, closure_down, msq, mgl)

    #run the datacard
    fake_mass_name = str(msq)[1:3]+str(mgl)[:3]
    combine_n_option  = str(msq) + "_" + str(mgl)

    curr_dir = os.getcwd()
    
    os.chdir(output_dir)

    if not options.nocombine:
        print "--Running combine"
        print "M1: %i M2 %i" % (msq, mgl)
        #toys

        #run asymptotic to get estimates of where to scan
        asymptotic_cmd = "combine -M Asymptotic -n .%s  %s" % (combine_n_option, data_card_name)
        asmyptotic_file = "higgsCombine.%s.Asymptotic.mH120.root" % combine_n_option
        os.system(asymptotic_cmd)
        
        #parse the asymptotic result
        limit_file = rt.TFile(asmyptotic_file)
        limit_tree = limit_file.Get("limit")        
        asym_result = limit_result(msq, mgl, limit_tree)

        #get the values of r from asmyptotic
        exp = asym_result.get_exp()
        exp_p = asym_result.get_exp_p()
        exp_m = asym_result.get_exp_m()

        #build the option for naming the toy files outputted from toys
        exp_n_option = combine_n_option + "_exp"
        exp_p_n_option = combine_n_option + "_exp_p"
        exp_m_n_option = combine_n_option + "_exp_m"

        #build the toy file names using the option flag
        exp_file = "higgsCombine.%s.HybridNew.mH120.12345.root" % exp_n_option
        exp_p_file = "higgsCombine.%s.HybridNew.mH120.12345.root" % exp_p_n_option
        exp_m_file = "higgsCombine.%s.HybridNew.mH120.12345.root" % exp_m_n_option
        
        #build the commands for running the toys

        toys_exp_cmd = "combine -M HybridNew --testStat LHC --freq --singlePoint %f --clsAcc 0 -T 3334 --fullBToys -i 1 %s --saveHybridResult --saveToys  -n .%s" % (exp, data_card_name, exp_n_option)
        toys_exp_p_cmd = "combine -M HybridNew --testStat LHC --freq --singlePoint %f --clsAcc 0 -T 3333 --fullBToys -i 1 %s --saveHybridResult --saveToys  -n .%s" % (exp_p, data_card_name, exp_p_n_option)
        toys_exp_m_cmd = "combine -M HybridNew --testStat LHC --freq --singlePoint %f --clsAcc 0 -T 3333 --fullBToys -i 1 %s --saveHybridResult --saveToys  -n .%s" % (exp_m, data_card_name, exp_m_n_option)

        #run the commands
        print "\nEXPECTED","\n", toys_exp_cmd
        os.system(toys_exp_cmd)
        print "\n+1SIGMA","\n", toys_exp_p_cmd
        os.system(toys_exp_p_cmd)
        print "\n-1SIGMA","\n", toys_exp_m_cmd
        os.system(toys_exp_m_cmd)
        
        grid_file_name = "higgsGrid%s.root" % combine_n_option

        #hadd the files together
        hadd_cmd = "hadd %s higgsCombine.%s*.root" % (grid_file_name, combine_n_option)
        print hadd_cmd,"\n"
        os.system(hadd_cmd)

        parse_exp_cmd = "combine -M HybridNew --testStat LHC --freq --grid %s --expectedFromGrid 0.5 -n .%s.Expected0p5 %s" % (grid_file_name,  combine_n_option, data_card_name)
        parse_exp_m_cmd = "combine -M HybridNew --testStat LHC --freq --grid %s --expectedFromGrid 0.16 -n .%s.Expected0p16 %s" % (grid_file_name, combine_n_option, data_card_name)
        parse_exp_p_cmd = "combine -M HybridNew --testStat LHC --freq --grid %s --expectedFromGrid 0.84 -n .%s.Expected0p84 %s" % (grid_file_name, combine_n_option, data_card_name)
        parse_obs_cmd = "combine -M HybridNew --testStat LHC --freq --grid %s -n .%s.Observed %s" % (grid_file_name, combine_n_option, data_card_name)

        print parse_exp_cmd, "\n", parse_exp_m_cmd, "\n", parse_exp_p_cmd, "\n", parse_obs_cmd, "\n"

        print "\n\n\n\n"
        print "\nEXPECTED","\n", parse_exp_cmd
        os.system(parse_exp_cmd)
        print "\n+1SIGMA", "\n",parse_exp_p_cmd
        os.system(parse_exp_p_cmd)
        print "\n-1SIGMA", "\n",parse_exp_m_cmd
        os.system(parse_exp_m_cmd)
        print "\nOBSERVED", "\n",parse_obs_cmd
        os.system(parse_obs_cmd)


    if options.no_limit:
        print "RUNNING TOYS COMPLETE -- EXITING...."
        exit(0)

    #IF WE ARENT RUNNIING COMBINE WE CAN JUST GO GRAB THE FILES
    parse_exp_file_name = "higgsCombine.%s.Expected0p5.HybridNew.mH120.quant0.500.root" % combine_n_option
    parse_exp_m_file_name = "higgsCombine.%s.Expected0p16.HybridNew.mH120.quant0.160.root" % combine_n_option
    parse_exp_p_file_name = "higgsCombine.%s.Expected0p84.HybridNew.mH120.quant0.840.root" % combine_n_option
    parse_obs_file_name = "higgsCombine.%s.Observed.HybridNew.mH120.root" % combine_n_option

    #get the root file for each result
    exp_file = rt.TFile(parse_exp_file_name)
    exp_p_file = rt.TFile(parse_exp_p_file_name)
    exp_m_file = rt.TFile(parse_exp_m_file_name)
    obs_file = rt.TFile(parse_obs_file_name)

    #get the tree from each root file
    exp_tree = exp_file.Get("limit")
    exp_p_tree = exp_p_file.Get("limit")
    exp_m_tree = exp_m_file.Get("limit")
    obs_tree = obs_file.Get("limit")        
    
    result  = toy_result(msq, mgl, exp_tree, exp_m_tree, exp_p_tree, obs_tree)

    result.print_summary()
    
    #add it to the full grid of results
    grid_results.add_point(msq, mgl, result)

    #move back to the home directory
    os.chdir(curr_dir)

grid_results.print_summary()    

(exp_hist, exp_p_hist, exp_m_hist, obs_hist, obs_p_hist, obs_m_hist, exp_theory, delta_hist) = grid_results.build_hist_limits()
limit_hists = [exp_hist, exp_p_hist, exp_m_hist, obs_hist, obs_p_hist, obs_m_hist, delta_hist]

output_file.cd()

#write the systematic derived from the fake
#fake_exp_hist.Write("closure_error")
closure_up.Write("closure_up")
closure_down.Write("closure_down")

#cross section excluded
exp_theory.Write("exp_theory")

#solid line for the observation
print "BUILDING OBSERVED BAND"
band_obs =  grid_results.build_band(obs_hist)
graph_obs = grid_results.build_graph_from_band(band_obs)
graph_obs.SetLineWidth(6)
graph_obs.SetLineStyle(9)
graph_obs.Write("obs_graph")

print "\nBUILDING OBSERVED BAND +"
band_obs_p =  grid_results.build_band(obs_p_hist)
graph_obs_p = grid_results.build_graph_from_band(band_obs_p)
graph_obs_p.SetLineWidth(6)
#graph_obs_p.SetLineStyle(9)
graph_obs_p.Write("obs_p_graph")

print "\nBUILDING OBSERVED BAND -"
band_obs_m =  grid_results.build_band(obs_m_hist)
graph_obs_m = grid_results.build_graph_from_band(band_obs_m)
graph_obs_m.SetLineWidth(6)
#graph_obs_m.SetLineStyle(9)
graph_obs_m.Write("obs_m_graph")

print "\nBUILDING EXPECTED BAND"
band_exp = grid_results.build_band(exp_hist)
graph_exp = grid_results.build_graph_from_band(band_exp)
graph_exp.SetLineStyle(7)
graph_exp.SetLineColor(rt.kRed)
graph_exp.SetLineWidth(6)
graph_exp.Write("exp_0_graph")

#lines for the expected region
print "\nBUILDING EXPECTED BAND +"
band_exp_p = grid_results.build_band(exp_p_hist)
graph_exp_p = grid_results.build_graph_from_band(band_exp_p)
graph_exp_p.SetLineStyle(7)
graph_exp_p.SetLineWidth(3)
graph_exp_p.SetLineColor(rt.kAzure-1)
graph_exp_p.Write("exp_p_graph")

print "\nBUILDING EXPECTED BAND -"
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
