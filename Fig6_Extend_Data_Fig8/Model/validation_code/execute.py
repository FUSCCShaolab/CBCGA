import os
import CBGCA_validation_run


def main():
    feats = ['Clin', 'IHC', 'IHC_Clin', 'Met_Clin', 'Met_IHC_Clin', 'Met_Path_Clin', 'Met_Path_IHC_Clin','Met_Path_Rad_Clin',
             'Met_Path_Rad_IHC_Clin','Met_Rad_Clin','Met_Rad_IHC_Clin','Metab','Path','Path_Clin','Path_IHC_Clin','Path_Rad_Clin',
             'Path_Rad_IHC_Clin','Rad','Rad_Clin','Rad_IHC_Clin','RNA','RNA_Clin','RNA_IHC_Clin','RNA_Met_Clin','RNA_Met_IHC_Clin',
             'RNA_Met_Path_Clin','RNA_Met_Path_IHC_Clin','RNA_Met_Path_Rad_Clin','RNA_Met_Path_Rad_IHC_Clin','RNA_Met_Rad_Clin','RNA_Met_Rad_IHC_Clin','RNA_Met_Rad_Clin','RNA_Met_Rad_IHC_Clin','RNA_Path_Clin','RNA_Path_IHC_Clin',
             'RNA_Path_Rad_Clin','RNA_Path_Rad_IHC_Clin','RNA_Rad_Clin','RNA_Rad_IHC_Clin']

#RNA MET RAD CLIN
    rcuts = [0.6]
    knum = [5]
    random_status = [111]
    parameters = [rcuts, knum, random_status,feats]
    import os, itertools, time
    parameter_combinations = list(itertools.product(*parameters))  # paramater的所有可能组合



    for i,combi in enumerate(parameter_combinations):
        rcut_i, knum_i,rd_i,feat_i= combi
        CBGCA_validation_run.main(feat_i,rcut_i,knum_i,rd_i)

if __name__=='__main__':
    main()
