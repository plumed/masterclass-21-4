import math
import sys

# read FILE with CVs and weights
FILENAME_ = sys.argv[1]
# number of CVs 
NCV_ = int(sys.argv[2])
# read minimum, maximum and number of bins for FES grid
gmin = []; gmax = []; nbin = []
for i in range(0, NCV_):
    i0 = 3*i + 3 
    gmin.append(float(sys.argv[i0]))
    gmax.append(float(sys.argv[i0+1]))
    nbin.append(int(sys.argv[i0+2]))
# read KBT_
KBT_ = float(sys.argv[3*NCV_+3])
# block size 
BSIZE_ = int(sys.argv[-1])

def get_indexes_from_index(index, nbin):
    indexes = []
    # get first index
    indexes.append(index%nbin[0])
    # loop
    kk = index
    for i in range(1, len(nbin)-1):
        kk = ( kk - indexes[i-1] ) / nbin[i-1]
        indexes.append(kk%nbin[i])
    if(len(nbin)>=2):
      indexes.append( ( kk - indexes[len(nbin)-2] ) / nbin[len(nbin) -2] )
    return tuple(indexes) 

def get_indexes_from_cvs(cvs, gmin, dx):
    keys = []
    for i in range(0, len(cvs)):
        keys.append(int( round( ( cvs[i] - gmin[i] ) / dx[i] ) ))
    return tuple(keys)

def get_points(key, gmin, dx):
    xs = []
    for i in range(0, len(key)):
        xs.append(gmin[i] + float(key[i]) * dx[i])
    return xs

# define bin size
dx = []
for i in range(0, NCV_):
    dx.append( (gmax[i]-gmin[i])/float(nbin[i]-1) )

# total numbers of bins
nbins = 1
for i in range(0, len(nbin)): nbins *= nbin[i]

# read file and store lists 
cv_list=[]; w_list=[]
for lines in open(FILENAME_, "r").readlines():
    riga = lines.strip().split()
    # check format
    if(len(riga)!=NCV_ and len(riga)!=NCV_+1):
      print (FILENAME_,"is in the wrong format!")
      exit()
    # read CVs
    cvs = []
    for i in range(0, NCV_): cvs.append(float(riga[i]))
    # keys are tuples of CV indices on the grid
    key = get_indexes_from_cvs(cvs, gmin, dx)
    # read weight, if present
    if(len(riga)==NCV_+1):
      w = float(riga[NCV_])
    else: w = 1.0
    # store into CV and weight lists
    cv_list.append(key)
    w_list.append(w)

# total number of data points
ndata = len(cv_list)
# number of blocks
nblock = int(ndata/BSIZE_)

# prepare global lists of histo dictionaries
# and normalizations, one entry per block
histo_l = []; norm_l = []

# not optimized for speed, but to be readable
# cycle on blocks
for iblock in range(0, nblock):
    # define range
    i0 = iblock * BSIZE_ 
    i1 = i0 + BSIZE_
    # build histogram dictionary
    # keys are tuples of CV indices on the grid
    histo = {}; norm = 0.0
    # initialize histogram
    for i in range(0, nbins):
        # get the indexes in the multi-dimensional grid
        key = get_indexes_from_index(i, nbin)
        # set histogram to zero
        histo[key] = 0.0
    # cycle on points in the block
    for i in range(i0, i1):
        # get weight
        w = w_list[i]
        # and CVs tuple of indexes
        cv = cv_list[i]
        # increase norm
        norm += w
        # update histogram
        histo[cv] += w
    # normalization of the block
    for key in histo: histo[key] /= norm
    # store in global lists
    histo_l.append(histo)
    norm_l.append(norm)

# now we calculate weighted average across blocks
# for each point in the histogram
dict_ave = {}
# cycle on keys - let's take them from histogram of first block
for key in histo_l[0]:
    # set average and normalization to zero
    ave = 0.0; norm = 0.0
    # cycle on blocks
    for iblock in range(0, nblock):
        # weight of the block
        w = norm_l[iblock]
        # increment normalization
        norm += w
        # increment weighted average
        ave += w * histo_l[iblock][key]
    # normalize and store in dictionary for average across blocks
    dict_ave[key] = ave / norm

# and the variance
dict_var = {}
# cycle on keys - let's take them from histogram of first block
for key in histo_l[0]:
    # set variance and normalization to zero
    var = 0.0; norm = 0.0
    # cycle on blocks
    for iblock in range(0, nblock):
        # weight of the block
        w = norm_l[iblock]
        # increment normalization
        norm += w
        # increment variance
        var += math.pow(w * (histo_l[iblock][key]-dict_ave[key]), 2.0) 
    # normalize and store in dictionary for variance across blocks 
    dict_var[key] = var / norm / norm

# now, print out fes and error 
log = open("fes."+str(BSIZE_)+".dat", "w")
# this is needed to add a blank line
xs_old = []
for i in range(0, nbins):
    # get the indexes in the multi-dimensional grid
    key = get_indexes_from_index(i, nbin)
    # get CV values for that grid point
    xs = get_points(key, gmin, dx)
    # add a blank line for gnuplot
    if(i == 0):
      xs_old = xs[:] 
    else:
      flag = 0
      for j in range(1,len(xs)):
          if(xs[j] != xs_old[j]):
            flag = 1
            xs_old = xs[:] 
      if (flag == 1): log.write("\n")
    # print value of CVs
    for x in xs:
        log.write("%12.6lf " % x)
    # calculate fes and error
    if key in dict_ave:
       # fes
       fes = -KBT_ * math.log(dict_ave[key]) 
       # variance fes
       varf = math.pow( KBT_ / dict_ave[key], 2.0) * dict_var[key]
       # error fes 
       errf = math.sqrt(varf)
       # printout
       log.write("   %12.6lf %12.6lf\n" % (fes, errf))
    else:
       log.write("       Inf Inf\n")
log.close()
