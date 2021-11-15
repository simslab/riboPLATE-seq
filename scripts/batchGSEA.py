import os
from sys import argv

indir = argv[0]#'../workingdir/','/home/ubuntu/datadrive/GSEA/'
outdir = argv[1]#'../gsea/','/home/ubuntu/datadrive/GSEA/out/'
path_to_gsea = argv[2]#'../scripts/GSEA_4.1.0/','/home/ubuntu/datadrive/GSEA_4.1.0/'

gmts = [g for g in os.listdir(indir) if g.split('.')[-1]=='gmt']

for rnk in ['PP_30','PP_6']:#,'PP242','MNK1','4EGI1','BKM','AZD','MNK1BKM','PP242BKM','PP242MNK1']:
    for gmt in ['riboPLATE_combos_excl','riboPLATE','riboPLATE_excl']:
        cmd = "bash {gspath}gsea-cli.sh GSEAPreranked -gmx {indir}{gmt}.gmt -collapse No_Collapse -mode Max_probe -norm meandiv -nperm 1000 -rnk {indir}{rnk}.rnk -scoring_scheme classic -rpt_label {rnk}_{gmt} -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 0 -rnd_seed timestamp -set_max 1000 -set_min 5 -zip_report false -out {outDir}".format(gspath=path_to_gsea,indir=indir,gmt=gmt,rnk=rnk,outDir=outdir)
        print(cmd)
        os.system(cmd)