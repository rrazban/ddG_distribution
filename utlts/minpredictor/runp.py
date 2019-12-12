import prody as dy
import csv
import numpy as np
import asa_np 
import time


# define residue ASA reference, store in dict
ASA_REF = dict()

ASA_REF['ALA']= 121.93
ASA_REF['CYS']= 173.46
ASA_REF['ASP']= 163.81
ASA_REF['GLU']= 207.33
ASA_REF['PHE']= 236.16
ASA_REF['GLY']= 92.5
ASA_REF['HIS']= 211.05
ASA_REF['ILE']= 220.13
ASA_REF['LYS']= 247.02
ASA_REF['LEU']= 207.26
ASA_REF['MET']= 221.3
ASA_REF['ASN']= 176.74
ASA_REF['PRO']= 173.45
ASA_REF['GLN']= 212.42
ASA_REF['ARG']= 260.07
ASA_REF['SER']= 154.56
ASA_REF['THR']= 198.38
ASA_REF['VAL']= 208.23
ASA_REF['TRP']= 278.45
ASA_REF['TYR']= 250.6


##########################

def getChains(pdb, chains1):

    protein = pdb.select('chain ' + chains1[0]).copy()
   
    return protein

def determineSurfaceAtoms(PDB, PDB_ASA, ref_percent=0.2):
    
    # return a binary vector: 1 if atom is on surface, 0 ow
#    surface_indicator = np.zeros(np.shape(PDB_ASA))

    res_iter = PDB.iterResidues()
    surface_indicator = np.zeros(max(PDB.getResnums())+2)	#do it by residue	#+2 cuz found incompatibility with biopython
    asas = np.zeros(max(PDB.getResnums())+2)
    for res in res_iter:
	resnum = res.getResnum()
        atoms_noh = res.select('noh')
        idx_a = [int(a.getIndex()) for a in atoms_noh]
        residueASA = np.sum(PDB_ASA[idx_a])
        residueASA1 = np.sum(PDB_ASA[idx_a][4:])
        asas[resnum] = residueASA1#/ASA_REF[res.getResname()]	#if divide, it becomes RSA
        # there may be a strange residue in there 
        try:
            ASA_REF[res.getResname()]
            if residueASA > ref_percent*float(ASA_REF[res.getResname()]):
               # surface_indicator[idx_a] = np.ones(len(idx_a))
                surface_indicator[resnum] = 1
        except KeyError:
            continue

    return surface_indicator, asas 

def scoreBenchmark(benchmark_file):
    
    dy.confProDy(verbosity='none')

    resultsFile = open('output.csv', 'w+')
    csvWriter = csv.writer(resultsFile)
    csvWriter.writerow(['PDB', 'dG_est'])

    csvReader  = csv.reader(open(benchmark_file, 'rU'))

    for row in csvReader:

        PDB_ID = row[0]
        print 'Working on ' + PDB_ID + '...'

        dG_est= scoreOne(PDB_ID + '.pdb')
        results = [PDB_ID]
        results.extend(dG_est)
        csvWriter.writerow(results)

    resultsFile.close()
    print 'Job finished!'

def scoreOne(PDB_FILENAME):
    
    PDB = dy.parsePDB(PDB_FILENAME)
    PDB_ID = PDB.getTitle()

    PDB = PDB.select('stdaa and noh').copy()

    chains_list = np.unique(PDB.getChids())

    # getChids or unique may change sequential ordering of chids
    if chains_list[0] != PDB.getChids()[0]:
        chains_list = [chains_list[idx] for idx in [1,0]]

    # chains
    protein = getChains(PDB, chains_list[0])
  
    # get asa for complex and each chain
    asa_protein = asa_np.calculate_asa_np(protein, 1.4, 960)  
    # round asa values
    asa_protein = [round(asa_protein[idx],2) for idx in range(len(asa_protein))]

    # store ASA values in mass field
    protein.setMasses(asa_protein)
    
    # determine surface atoms, ###_ASA_indicator is number of atoms length binary vector.
    protein_ASA_indicator, asas  = determineSurfaceAtoms(protein, protein.getMasses(), 0.2)
    return protein_ASA_indicator, asas
    print protein_ASA_indicator
    sys.exit()
    ligand_ASA_indicator   = determineSurfaceAtoms(ligand,  ligand.getMasses(),  0.2)

    # compute desolvation terms
    deltaASA_Tyr = 0.0
    deltaASA_Ser = 0.0
    deltaASA_Cys = 0.0

    if complex_pdb.select('resname TYR') is not None:
        deltaASA_Tyr = round(np.sum(complex_pdb.select('resname TYR').select('sc').getMasses() - PDB.select('resname TYR').select('sc').getMasses()),2)

    if complex_pdb.select('resname SER') is not None:
        deltaASA_Ser = round(np.sum(complex_pdb.select('resname SER').select('sc').getMasses() - PDB.select('resname SER').select('sc').getMasses()),2)

    if complex_pdb.select('resname CYS') is not None:
        deltaASA_Cys = round(np.sum(complex_pdb.select('resname CYS').select('sc').getMasses() - PDB.select('resname CYS').select('sc').getMasses()),2)

    # arrange features and set fitted weights 
    feature_vector = np.array([deltaASA_Tyr, deltaASA_Ser, deltaASA_Cys])

    weights_vector = [0.00856255227704338, 0.013719207070153, 0.0315777976693462]
    bias = -8.470604361

    # our ddg estimate
    ddg_est = np.dot(feature_vector, weights_vector) + bias

    fv = [round(ddg_est,2)]
    return fv

def main():
    
    scoreBenchmark('input.txt')

if __name__ == "__main__":
    main()

