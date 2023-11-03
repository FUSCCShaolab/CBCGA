import CBCGA_leaveone
def main():
    import time
    feats = ['M168T190.POS', 'M297T233.POS', 'M149T91.NEG', 'M773T388.NEG', 'M778T361.NEG',
                'T.GLCM.IntensityMax.Mean', 'T.GLCM.IntensityMax.Skew',
                'T.Morph.CellEccentricities.Kurtosis', 'T.GLCM.IntensityMean.Kurtosis', 'Normal', 'GFRA1', 'TFF3',
                'CCNB1', 'CTSV', 'CDH3',
                'Clin', 'IHC']
    rcuts=[0.6]
    knum=[5]
    random_status = [111]
    parameters = [rcuts, knum, random_status,feats]
    import os, itertools, time
    parameter_combinations = list(itertools.product(*parameters))  # paramater的所有可能组合



    for i,combi in enumerate(parameter_combinations):
        rcut_i, knum_i,rd_i,feat_i= combi
        print(combi)
        CBCGA_leaveone.main(feat_i,rcut_i,knum_i,rd_i)

if __name__=='__main__':
    main()
