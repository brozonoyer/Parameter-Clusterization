Extremal grouping program works in 2 steps:
STEP1. Build membership of features to groups and save information necessary for calculation of extremal factors.
STEP2. Calculate value of extremal group factors for new data not used in the fist step.

This 2 steps procedure required because we need to calculate extremal factors for several number of new datasets.
It is not sufficient to calculate membership and value of vactors on one dataset only. 

Below is example of interface ( input/output data)  for the steps:

Input data for STEP1:
1. Input training CSV-file
2. Number of extremal groups
Output data:
1. Memberships to groups
2. Mean and StDev values of features in input CSV-file
3. Weights of features in groups

Example of possible command line for program:

EXTREMAL_GROUP.EXE   --STEP 1  --NUMBER_OF_GROUPS 3  --input train.csv --info EG.INFO

where  output file EG.INFO contains the following information:

ExtremalGroupsInfo
Built on features
1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67
nParams,67
nGroups,3
Mean
6.83166,7.14285,4.19451,10.2229,8.60503,6.74578,6.38718,4.12926,4.84961,3.77426,11.411,6.19396,6.81055,4.61077,3.99081,6.14946,5.36607,5.51821,3.98619,5.37709,5.06893,6.17836,6.05457,8.10023,7.90893,5.05052,7.01582,4.61559,4.47868,4.02643,6.93322,5.53372,3.6692,3.88734,3.50467,4.30804,9.27379,7.96678,5.13916,4.1282,4.81553,3.51179,4.84342,9.05342,3.56752,3.62891,4.38605,6.82728,4.65178,6.65249,6.13269,4.85898,6.08163,5.80491,9.54171,7.0093,9.40156,9.11138,6.81356,5.23482,6.96913,8.60813,6.62739,5.44933,9.14018,6.56863,5.72657
StDev
438.501,24.9587,8.2054,31.4502,24.3567,17.4655,17.1867,11.158,16.884,7.84514,36.2352,21.7007,26.5219,12.8314,9.10585,16.432,12.0767,13.352,6.55471,13.8138,12.5706,15.9241,14.2772,20.8703,19.1631,10.1842,18.5083,9.10134,13.1405,8.59143,19.5972,11.7194,7.64928,8.5617,6.97142,11.2649,24.4187,24.5301,12.769,6.90144,15.6681,4.42167,15.2945,28.6113,5.76963,7.42456,11.6911,18.1932,10.2232,16.7062,14.4815,11.9937,16.4683,17.1444,27.5662,17.4276,31.5295,27.0091,20.4731,10.9775,18.5779,21.1949,14.7709,14.0697,24.5069,17.3969,14.6586
inGroup, group number starts from 1
2,3,3,1,1,1,2,2,3,1,1,1,1,2,2,2,3,3,3,2,2,2,2,3,3,3,2,2,2,2,3,3,3,1,1,1,1,3,1,1,1,1,1,1,1,1,2,3,3,3,3,2,2,2,3,3,3,3,2,2,2,3,3,1,1,1,2
Weights
-0.186458,-0.0816229,-0.064139,0.16274,0.180544,0.163993,-0.128864,-0.086013,-0.0308409,0.0306928,0.134885,0.104327,0.106939,-0.0859965,0.0646137,0.0938427,0.0742317,0.122684,0.121096,0.115043,0.167587,0.20639,0.161277,0.137589,0.226445,0.191443,0.147012,0.17435,0.155488,0.0732077,0.150944,0.20764,0.067144,-0.157152,-0.118166,-0.103356,-0.106509,0.110805,-0.206371,-0.232019,-0.179604,-0.167401,-0.0723795,-0.0851955,-0.0487207,-0.130955,-0.0812782,-0.117166,-0.105671,-0.139035,-0.0974891,-0.146558,-0.185041,-0.152154,-0.152799,-0.13035,-0.137212,-0.0572567,-0.180454,-0.194798,-0.129724,-0.130205,-0.123183,0.132646,0.123184,0.122309,-0.0989264
  

Input data for STEP2:
1. Input test CSV-file
2. Memberships to groups
3. Mean and StDev values of features.
4. Weights of features in groups

Output data:
1. Value of factors for test CSV-file
                                                                                   
Example of possible command line for program:
EXTREMAL_GROUP.EXE   --STEP 2  --input test.csv --info EG.INFO   --out  factors.csv  
**************************************************************************************************************************************

Example of input and output data for 2 steps:

Input data for STEP1:
1. m255_train.csv - Input training CSV-file
2. 3 - Number of extremal groups
Output data:
_eg.info - (Memberships to groups, Mean and StDev, Weights)


Input data for STEP2:
1. m255_test.csv - Input test CSV-file
2. _eg.info - (Memberships to groups, Mean and StDev, Weights)

Output data:
1. m255_test_FACTORS.csv - Value of factors for test CSV-file