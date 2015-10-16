#!/bin/bash
module load FSL AFNI

mkdir registered
subjects=$(ls *.nii.gz | cut -d '_' -f 3 | uniq)

# register (fsl flirt) all images to the first image
for subj in ${subjects}; do
    files=$(ls *${subj}*)
    files=(${files})
    nfiles=${#files[@]}
    count=1

    for (( c=0; c<${nfiles}; c++ )); do
        flirt -in ${files[${c}]} -ref ${files[0]} -out registered/${subj}_${c}.nii.gz
    done
done

# average these images
cd registered
subjects=$(ls *.nii.gz | cut -d '_' -f 1 | uniq)

for subj in ${subjects}; do
    files=$(ls ${subj}*)
    mkdir SPN01_MRC_${subj}_01
    3dMean -prefix SPN01_MRC_${subj}_01/SPN01_MRC_${subj}_01_01_T1_file.nii.gz ${files}
done

