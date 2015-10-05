Imitate Observe
---------------
+ data preprocessed with /archive/data-2.0/code/datman/assets/150409-compcor-nonlin-8fwhm.sh 
+ event-related GLM performed using AFNI. Modelled using 5 tent functions (FIR) following each event. Sample below:

    3dDeconvolve \
    -input /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_func_MNI-nonlin.IM.01.nii.gz \
    -mask /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_anat_EPI_mask_MNI.nii.gz \
    -ortvec /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_motion.1D motion_paramaters \
    -polort 4 \
    -num_stimts 6 \
    -local_times \
    -jobs 8 \
    -x1D /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_glm_IM_1stlevel_design.mat \
    -stim_label 1 IM_AN -stim_times 1 /archive/data-2.0/SPINS/metadata/design/IM_event-times_AN.1D 'TENT(0,15,5)' \
    -stim_label 2 IM_FE -stim_times 2 /archive/data-2.0/SPINS/metadata/design/IM_event-times_FE.1D 'TENT(0,15,5)' \
    -stim_label 3 IM_FX -stim_times 3 /archive/data-2.0/SPINS/metadata/design/IM_event-times_FX.1D 'TENT(0,15,5)' \
    -stim_label 4 IM_HA -stim_times 4 /archive/data-2.0/SPINS/metadata/design/IM_event-times_HA.1D 'TENT(0,15,5)' \
    -stim_label 5 IM_NE -stim_times 5 /archive/data-2.0/SPINS/metadata/design/IM_event-times_NE.1D 'TENT(0,15,5)' \
    -stim_label 6 IM_SA -stim_times 6 /archive/data-2.0/SPINS/metadata/design/IM_event-times_SA.1D 'TENT(0,15,5)' \
    -gltsym 'SYM: -1*IM_FX +0*IM_NE +0.25*IM_AN +0.25*IM_FE +0.25*IM_HA +0.25*IM_SA' \
    -glt_label 1 emot-fix \
    -gltsym 'SYM: +0*IM_FX -1*IM_NE +0.25*IM_AN +0.25*IM_FE +0.25*IM_HA +0.25*IM_SA' \
    -glt_label 2 emot-neut \
    -fitts   /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_glm_IM_1stlvl_explained.nii.gz \
    -errts   /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_glm_IM_1stlvl_residuals.nii.gz \
    -bucket  /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_glm_IM_1stlvl.nii.gz \
    -cbucket /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_glm_IM_1stlvl_allcoeffs.nii.gz \
    -fout -tout -xjpeg /archive/data-2.0/SPINS/data/imob/SPN01_CMH_0001_01/SPN01_CMH_0001_01_glm_IM_1stlevel_design.jpg

+ Group level modelling was performed using 3dANOVA3: /projects/jdv/data/spins/1yr/outputs/bin/run-glm_IMOB_allsites.sh
