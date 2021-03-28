## ---------------------------
## [script name] ps08_VIS_group_ROI_contacts_brain.py
##
## SCRIPT to ...
##
## By Shuai Wang, [date] 2021-03-27
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

## set environment (packages, functions, working path etc.)
# load up packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from nilearn import image, plotting
# setup working path
mdir='/media/wang/BON/Projects/CP01'
wdir=os.path.join(mdir,'SEEG_LectureVWFA/derivatives/mia_SEEG_LectureVWFA')
gdir=os.path.join(wdir,'group')
rdir=os.path.join(gdir,'ROIs-contacts_brain')
kdir=os.path.join(mdir,'SEEG_LectureVWFA/derivatives/masks/AAL3')
# useful parameters
## ---------------------------

## read ROIs and contacts information
# read MNI-ROIs list
froi=os.path.join(gdir,'group_space-MNI_contacts-bipolar_ROI-AAL3.csv')
roi_mni=pd.read_csv(froi)
# extract information
roi_mni=roi_mni[roi_mni.conditions=='AAp']
roi_full=np.unique(roi_mni.AAL3)
roi_full=roi_full[roi_full!='no_label_found']  # remove nonsense label
## ---------------------------

## plot contacts for each ROI
for roi in roi_full:
  print("Plot ROI %s and its contacts in glass brain ......\n" % roi)
  froi=os.path.join(kdir,"AAL3_%s.nii" % roi)
  roi_mni_tmp=roi_mni[roi_mni.AAL3==roi]
  # create color strings according to subjects (... to be optimized)
  sub_colors=roi_mni_tmp.subjects
  sub_colors=[s.replace('sub-04','red') for s in sub_colors]
  sub_colors=[s.replace('sub-05','green') for s in sub_colors]
  sub_colors=[s.replace('sub-08','blue') for s in sub_colors]
  nsubjs=len(np.unique(sub_colors))
  # plot glass brain
  xyz=roi_mni_tmp[['mni_x','mni_y','mni_z']]
  nc=np.size(xyz,axis=0)  # number of contacts within this ROI
  d=plotting.plot_connectome(adjacency_matrix=np.zeros((nc,nc)),node_coords=xyz,node_color=sub_colors,node_size=8)
  d.add_contours(froi,colors='black',linewidths=0.5,linestyles='dotted')
  # add legend
  p_v= Rectangle((0, 0), 1, 1, fc="red")
  p_h= Rectangle((0, 0), 1, 1, fc="green")
  p_f= Rectangle((0, 0), 1, 1, fc="blue")
  plt.legend([p_v, p_h, p_f], ["sub-04", "sub-05", "sub-08"],loc=1)
  # add title
  d.title("%s : %d contacts" % (roi,nc))
  # save figure
  ffig=os.path.join(rdir,"group_N%d_space-MNI_ROI-%s_contacts.png" % (nsubjs,roi))
  d.savefig(ffig,dpi=600)
## ---------------------------
