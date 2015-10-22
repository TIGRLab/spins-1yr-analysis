3dANOVA3 -DAFNI_FLOATIZE=YES -type 5 -alevels 2 -blevels 2 -clevels 10 \
    -dset 1 1 1 SPN01_CMH_0004_01_IM.nii.gz[20] \
    -dset 1 2 1 SPN01_CMH_0004_01_OB.nii.gz[20] \
    -dset 1 1 2 SPN01_CMH_0005_01_IM.nii.gz[20] \
    -dset 1 2 2 SPN01_CMH_0005_01_OB.nii.gz[20] \
    -dset 1 1 3 SPN01_CMH_0007_01_IM.nii.gz[20] \
    -dset 1 2 3 SPN01_CMH_0007_01_OB.nii.gz[20] \
    -dset 1 1 4 SPN01_CMH_0008_01_IM.nii.gz[20] \
    -dset 1 2 4 SPN01_CMH_0008_01_OB.nii.gz[20] \
    -dset 1 1 5 SPN01_CMH_0011_01_IM.nii.gz[20] \
    -dset 1 2 5 SPN01_CMH_0011_01_OB.nii.gz[20] \
    -dset 1 1 6 SPN01_CMH_0014_01_IM.nii.gz[20] \
    -dset 1 2 6 SPN01_CMH_0014_01_OB.nii.gz[20] \
    -dset 1 1 7 SPN01_CMH_0015_01_IM.nii.gz[20] \
    -dset 1 2 7 SPN01_CMH_0015_01_OB.nii.gz[20] \
    -dset 1 1 8 SPN01_CMH_0016_01_IM.nii.gz[20] \
    -dset 1 2 8 SPN01_CMH_0016_01_OB.nii.gz[20] \
    -dset 1 1 9 SPN01_CMH_0020_01_IM.nii.gz[20] \
    -dset 1 2 9 SPN01_CMH_0020_01_OB.nii.gz[20] \
    -dset 1 1 10 SPN01_CMH_0023_01_IM.nii.gz[20] \
    -dset 1 2 10 SPN01_CMH_0023_01_OB.nii.gz[20] \
    -dset 2 1 1 SPN01_CMH_0001_01_IM.nii.gz[20] \
    -dset 2 2 1 SPN01_CMH_0001_01_OB.nii.gz[20] \
    -dset 2 1 2 SPN01_CMH_0002_01_IM.nii.gz[20] \
    -dset 2 2 2 SPN01_CMH_0002_01_OB.nii.gz[20] \
    -dset 2 1 3 SPN01_CMH_0003_01_IM.nii.gz[20] \
    -dset 2 2 3 SPN01_CMH_0003_01_OB.nii.gz[20] \
    -dset 2 1 4 SPN01_CMH_0009_01_IM.nii.gz[20] \
    -dset 2 2 4 SPN01_CMH_0009_01_OB.nii.gz[20] \
    -dset 2 1 5 SPN01_CMH_0012_01_IM.nii.gz[20] \
    -dset 2 2 5 SPN01_CMH_0012_01_OB.nii.gz[20] \
    -dset 2 1 6 SPN01_CMH_0013_01_IM.nii.gz[20] \
    -dset 2 2 6 SPN01_CMH_0013_01_OB.nii.gz[20] \
    -dset 2 1 7 SPN01_CMH_0021_01_IM.nii.gz[20] \
    -dset 2 2 7 SPN01_CMH_0021_01_OB.nii.gz[20] \
    -dset 2 1 8 SPN01_CMH_0024_01_IM.nii.gz[20] \
    -dset 2 2 8 SPN01_CMH_0024_01_OB.nii.gz[20] \
    -dset 2 1 9 SPN01_CMH_0027_01_IM.nii.gz[20] \
    -dset 2 2 9 SPN01_CMH_0027_01_OB.nii.gz[20] \
    -dset 2 1 10 SPN01_CMH_0029_01_IM.nii.gz[20] \
    -dset 2 2 10 SPN01_CMH_0029_01_OB.nii.gz[20] \
-fa group -fb cond -fab groupXcond -acontr 1 0 HC -acontr 0 1 SC -bcontr 1 0 IM -bcontr 0 1 OB -acontr 1 -1 HCvsSC -bcontr 1 -1 IMvsOB -Abcontr 1 : 1 -1 IMvsOB_HC -Abcontr 2 : 1 -1 IMvsOB_SC -aBcontr 1 -1 : 1 HCvsSC_IM -aBcontr 1 -1 : 2 HCvsSC_OB -aBcontr 1  0 : 1 HC_IM -aBcontr 0  1 : 1 SC_IM -aBcontr 1  0 : 2 HC_OB -aBcontr 0  1 : 2 SC_OB -bucket ANOVA-cmh.nii.gz

