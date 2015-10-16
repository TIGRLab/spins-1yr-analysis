#!/usr/bin/env python
"""
This maps results (generated in R from csv of engima dti outputs) back onto the FA skeleton.

Usage:
  map_edto_results_to_skel.py [options] <results-csv> <resultcolname> <out-nii>

Arguments:
  <results-csv>            A csv of the resutls from enigma dti (with header) ROIs. The ROI names shoud be the first column.
  <resultcolname>          The name (in the header) of the result to map to the surface.
  <out-nii>                The name for the output file (should end in .nii)

Options:
  -v,--verbose             Verbose logging
  --debug                  Debug logging in Erin's very verbose style
  -n,--dry-run             Dry run

DETAILS
Requires python enviroment with numpy and nibabel:
module load use.own datman/edickie

Requires FSL
module load FSL/5.0.7
"""
from docopt import docopt
import numpy as np
import nibabel as nib
import os
import tempfile
import shutil
import subprocess
import pandas as pd

arguments       = docopt(__doc__)
resultcsv       = arguments['<results-csv>']
resultcolname   = arguments['<resultcolname>']
outnii          = arguments['<out-nii>']
VERBOSE         = arguments['--verbose']
DEBUG           = arguments['--debug']
DRYRUN          = arguments['--dry-run']

###


### Erin's little function for running things in the shell
def docmd(cmdlist):
    "sends a command (inputed as a list) to the shell"
    if DEBUG: print ' '.join(cmdlist)
    if not DRYRUN: subprocess.call(cmdlist)

## make a tmpdir
tempdir = tempfile.mkdtemp()

#### some inputs for debugging
# resultcsv = '/projects/edickie/code/spins-1yr-analysis/dti/ttestresults/ttestres_FAresid_allsites.csv'
# resultcolname = 'cohen.D'
# outnii = '/projects/edickie/code/spins-1yr-analysis/dti/ttestresults/FAskel_effsize.nii.gz'

# Links to the template files - note may need to change if moving to another place
templatefolder = '/projects/edickie/code/hcp_extras/templates'
ROImap = os.path.join(templatefolder,'JHU-WhiteMatter-labels-1mm.nii')
Skelmask = os.path.join(templatefolder,'ENIGMA_DTI_FA_skeleton_mask.nii.gz')
lookup_table = os.path.join(templatefolder,'ENIGMA_look_up_table.csv')


# make a base for the output image which is the ROI map set to float of zero
docmd(['fslmaths', ROImap, '-mul', '0.0', outnii])

## load the ROI map into python using nibabel
roi_nib = nib.load(ROImap)
ROIniidata = roi_nib.get_data()

## load the output image into python using nibabel
out_nib = nib.load(outnii)
outniidata = out_nib.get_data()

## load both resutls into python pandas dataframe
results = pd.read_csv(resultcsv, sep=',', dtype=str, comment='#')
## load the lookup_table into python pandas dataframe
lookup = pd.read_csv(lookup_table, sep=',', dtype=str, comment='#')

## match up the lookup with the lookup with the results and find and replace in nifty output image
for i in range(0,len(lookup)):
    dlabelindex = lookup.KEY[i]   # dlabelindex is the numeric value (in the nifty) label
    dlabelname = lookup.NAME[i]   # dlabelname is the name of the ROI (matching the results csv)
    for j in range(0,len(results)):
        thisroi = results.iloc[j,0]  #thisroi is the name that should match dlabelname
        thisresult = results[resultcolname][j]  #this is the value that gets sent to the csv
        if dlabelname in thisroi:  #doing a search instead of a perfect match because many resutls have postix (i.e. '_thickavg')
            '''
            put thisresult value in voxels of output image
            where the dlabelindex is present
            '''
            outniidata[ROIniidata==int(dlabelindex)] = thisresult

## save the result - note nibabel always modified imagedata in place
nib.save(out_nib, outnii)

## now mask the results on the full (fat) JHU altas with the FA skeleton
docmd(['fslmaths', outnii, '-mul', Skelmask, outnii])

## remove the tmpdir
shutil.rmtree(tempdir)
