******************************************************
*                                                    *
*  Python Module to Compute Binding Energy           *
*                                                    *
*  Choi et al., Biophys. J. 108 (4), 795-798 (2015)  *
*                                                    *
*  Readme written by Jeong-Mo Choi, Nov-20-2014      *
*                                                    *
******************************************************

0. Code Authors: Sean Murphy, Dennis Lucarelli, and Jeong-Mo Choi

1. Requirements
  This python code requires following packages:
  - ProDy http://prody.csb.pitt.edu/
  - NumPy http://www.numpy.org/

2. How to Run
  1) Deposit all PDB files in the current directory
  2) Check if each PDB file contains only two chains
  3) Make the list of all PDBs (without extension) and save it as input.txt
  4) Execute "python run.py"
  5) The result will be released in output.csv in the CSV format

3. One example is given in directory ./example/
