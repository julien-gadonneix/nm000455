ReadMe for data for The neural activity of auditory conscious perception (https://doi.org/10.1016/j.neuroimage.2025.121041). 



*~BEHAVIOR~*
The raw behavioral files are in behavioraldata. With the exception of subjects 451 and 548, all behavior files contained represent the entirety of data for that subject. 451 and 548 each have 3 files: 2 that are partial files, and one that ends in the word '_combined'. For the purposes of behavioral analyses, the combined file should be used. 

Each behavior file includes the following:

allanswers_noisecalibration_all  
allanswers_noiserun
allanswers_noiserun_all
numberofnoiseruns
numberruns_total
pipelineversion_number
sessionnumber
subject_idnumber
tasktype
threshold_overall_whitenoise
threshold_overall_whitenoise_all

allanswers_noisecalibration_all and allanswers_noiserun_all (and allanswers_noiserun) contain a matrix with 29* columns. In allanswers_noisecalibration_all, many of these are blank because it uses the same overall structure as allanswers_noiserun_all, but has fewer answered questions. This is the meaning of the columns:
1. Was a sound heard?: for Calibration: 1 means yes; 2 means no response [no sound heard]. For noise_run, 1 means yes, 2 means no.
2. Reaction time to above. For Calibration: This is time after target sound was played. For noise_run, it is reaction time after question1 [above] was displayed. (time in s)
3. What sound was it? For calibration: NaN (sound ID not asked); For noise-run: 1 is whistle; 2 is laser; 3 is waterdrop
4. Reaction time to above (time in s)
5. Amplitude of sound (0-1) 
6. Where sound was presented (2=center)
7. Sample when sound was started
8. Sampling frequency of sound (hz)
9. frequency (always NaN)
10. Total duration of audio in trial (in s)
11. Duration of target sound (in s) 
12. Session number (starts counting at 1)
13. Number of runs completed (starts counting at 0)
14. Non-corrected answer to 1  
15. Answer to question 2
16. Task type (always 4)
17. Which button means yes? Button 1 
18. Which button means no?
19. Which white noise column was used? (always 1)
20. Which sample did the white noise start at?
21. How many total samples are there?
22. What target sound was actually played (0 is no sound; 1 is whistle; 2 is laser; 3 is waterdrop)
23. Start of Run (relative to task set up)
24. Start of Trial (relative to task set up)
25. Start of sound relative to trial start
26. Question 1 presentation time (relative to task set up)
27. Question 1 answer time (relative to task set up)
28. Question 2 presentation time (relative to task set up)
29. Question 2 answer time (relative to task set up)

*rows 23-29 were added on after data collection began as a method of troubleshooting, but are not used for analysis purposes

allanswers_noiserun contains the above, but just for the most recent run that was completed.

numberofnoiseruns and numberruns_total indicate the number of runs completed (noiseruns are the ones included in these analyses)

pipelineversion_number: internal record-keeping, for small updates in code
sessionnumber: how many sessions did it take to reach the final trial count
subject_idnumber: subject's internal ID
tasktype: internal record keeping; several types of tasks could be run
threshold_overall_whitenoise: most recent threshold calibration
threshold_overall_whitenoise_all: record of all thresholds used


*~Parcellation maps~*
For statistical analyses, 100 parcellations were produced to tile the cortical surface. These are contained in the folder parcellationmaps. Each parcellation map has its own folder, and 2 .mat files: one for left, and one for right. For each of these files, there is one column; the row index corresponds to the vertex # for that hemisphere. The number contained in that column corresponds to the parcellation number for that parcellation map. The folder also contains a visualization for that particular parcellation.



*~EEG ~*
Each subject folder contains 3 files:
L_MontageMap
R_MontageMap
and a sorted trials file

Montage Maps contain grey-matter electrode locations in the left and right hemispheres. Column 1 is electrode # and columns 2-4 are the 3d MNI coordinates of the electrode.

The sorted trials file contains EEG data for each patient. It contains the following:
BC_trials: channels x samples x trials array of perceived trials (“BC” abbreviation is old and used for backwards compatibility purposes) 
WG_trials: channels x samples x trials array of not perceived trials (“WG” abbreviation is old and used for backwards compatibility purposes) 
blank_TN_trials: channels x samples x trials array of True Negative blank trials (where no target sounds were presented) 
Where channels correspond to electrode #; samples contain 2 seconds at the patient's sampling rate; and trials are organized in the order of appearance in the behavior for that particular trial type. 
E.g: If the subject had the following behavior: 
trial 1 - correct 
trial 2 - incorrect
trial 3 - correct
trial 4 - blank
trial 5 - correct
trial 6 - incorrect
The z-dimension of the BC_-, WG_-, and blank_TN_trials would contain the information from the following trials:
BC_trials: trial 1 would be the first "page" (first of the 3rd dimensional element); trial 3 would be the 2nd page and trial 5 would be 3rd page
WG_trials: trial 2 would be the first page; trial 6 would be second page
blank_TN_trials: trial 4 would be the first page

Version Update
**06/04/2025** partitioned the EEG data into 4 parts to comply with the 10GB maximum single file size supported by Dryad; the original size is 36GB
