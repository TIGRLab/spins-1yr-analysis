3dttest++ -setA healthy \
    SPN01_MRC_0015 SPN01_MRC_0015_01_EA.nii.gz'[3]' \
    SPN01_MRC_0016 SPN01_MRC_0016_01_EA.nii.gz'[3]' \
    SPN01_MRC_0019 SPN01_MRC_0019_01_EA.nii.gz'[3]' \
    SPN01_MRC_0020 SPN01_MRC_0020_01_EA.nii.gz'[3]' \
    SPN01_MRC_0021 SPN01_MRC_0021_01_EA.nii.gz'[3]' \
    SPN01_MRC_0023 SPN01_MRC_0023_01_EA.nii.gz'[3]' \
    SPN01_MRC_0024 SPN01_MRC_0024_01_EA.nii.gz'[3]' \
    SPN01_MRC_0025 SPN01_MRC_0025_01_EA.nii.gz'[3]' \
    SPN01_MRC_0026 SPN01_MRC_0026_01_EA.nii.gz'[3]' \
    -setB schizophrenia \
    SPN01_MRC_0001 SPN01_MRC_0001_01_EA.nii.gz'[3]' \
    SPN01_MRC_0002 SPN01_MRC_0002_01_EA.nii.gz'[3]' \
    SPN01_MRC_0006 SPN01_MRC_0006_01_EA.nii.gz'[3]' \
    SPN01_MRC_0004 SPN01_MRC_0004_01_EA.nii.gz'[3]' \
    SPN01_MRC_0005 SPN01_MRC_0005_01_EA.nii.gz'[3]' \
    SPN01_MRC_0008 SPN01_MRC_0008_01_EA.nii.gz'[3]' \
    SPN01_MRC_0007 SPN01_MRC_0007_01_EA.nii.gz'[3]' \
    SPN01_MRC_0012 SPN01_MRC_0012_01_EA.nii.gz'[3]' \
    SPN01_MRC_0009 SPN01_MRC_0009_01_EA.nii.gz'[3]' \
    -prefix ttest_HC-SZ_modulator-mrc.nii.gz

3dttest++ -setA healthy \
    SPN01_MRC_0015 SPN01_MRC_0015_01_EA.nii.gz'[3]' \
    SPN01_MRC_0016 SPN01_MRC_0016_01_EA.nii.gz'[3]' \
    SPN01_MRC_0019 SPN01_MRC_0019_01_EA.nii.gz'[3]' \
    SPN01_MRC_0020 SPN01_MRC_0020_01_EA.nii.gz'[3]' \
    SPN01_MRC_0021 SPN01_MRC_0021_01_EA.nii.gz'[3]' \
    SPN01_MRC_0023 SPN01_MRC_0023_01_EA.nii.gz'[3]' \
    SPN01_MRC_0024 SPN01_MRC_0024_01_EA.nii.gz'[3]' \
    SPN01_MRC_0025 SPN01_MRC_0025_01_EA.nii.gz'[3]' \
    SPN01_MRC_0026 SPN01_MRC_0026_01_EA.nii.gz'[3]' \
    -prefix ttest_HC_modulator-mrc.nii.gz

3dttest++ -setA schizophrenia \
    SPN01_MRC_0001 SPN01_MRC_0001_01_EA.nii.gz'[3]' \
    SPN01_MRC_0002 SPN01_MRC_0002_01_EA.nii.gz'[3]' \
    SPN01_MRC_0006 SPN01_MRC_0006_01_EA.nii.gz'[3]' \
    SPN01_MRC_0004 SPN01_MRC_0004_01_EA.nii.gz'[3]' \
    SPN01_MRC_0005 SPN01_MRC_0005_01_EA.nii.gz'[3]' \
    SPN01_MRC_0008 SPN01_MRC_0008_01_EA.nii.gz'[3]' \
    SPN01_MRC_0007 SPN01_MRC_0007_01_EA.nii.gz'[3]' \
    SPN01_MRC_0012 SPN01_MRC_0012_01_EA.nii.gz'[3]' \
    SPN01_MRC_0009 SPN01_MRC_0009_01_EA.nii.gz'[3]' \
    -prefix ttest_SZ_modulator-mrc.nii.gz
