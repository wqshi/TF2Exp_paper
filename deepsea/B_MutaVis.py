#Usage: python B_MutaVis.py fastafilename  reference_sequence_predicitons mutated_sequence_predictions feature_index

import sys
import pandas as pd
import h5py
import numpy as np
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from pylab import *



mutwindow=1000

#change to use fasta given.
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
fasta_sequences=[str(seq.seq)[(( (len(seq.seq) - 1000) / 2)):( len(seq.seq) - (len(seq.seq) - 1000) / 2)] for seq in fasta_sequences]


mutpred2=np.asarray(h5py.File(sys.argv[3],'r')['pred']);
refpred2=np.asarray(h5py.File(sys.argv[2],'r')['pred']);
mutpred2=np.log2(mutpred2/(1-mutpred2+1e-12));
refpred2=np.log2(refpred2/(1-refpred2+1e-12));
mutpred2=(mutpred2[:3*mutwindow,:]+mutpred2[(mutpred2.shape[0]/2):((mutpred2.shape[0]/2)+3*mutwindow),:])/2;
refpred2=(refpred2[0,:]+refpred2[refpred2.shape[0]/2,:])/2;

if refpred2.ndim==1:
    diff=mutpred2[:3000,:]-refpred2[:]
else:
    diff=mutpred2[:3000,:]-refpred2[0,:]

feature_index=int(sys.argv[4])-1
if feature_index==-1:
    #feature_index=np.argmin(np.min(diff,axis=0))
    baseprob = np.asarray([0.039105,0.039264,0.033475,0.023512,0.060983,0.049581,0.03946,0.033393,0.034061,0.030136,0.036044,0.031867,0.038986,0.038193,0.028744,0.031794,0.035769,0.032283,0.03395,0.041033,0.037824,0.041241,0.032871,0.033799,0.033746,0.034082,0.034463,0.046585,0.039079,0.056989,0.045832,0.034852,0.054207,0.04587,0.039415,0.040048,0.033884,0.039052,0.036369,0.03952,0.04175,0.039414,0.050102,0.037017,0.043024,0.045057,0.043071,0.040406,0.028453,0.032615,0.045559,0.035012,0.048136,0.048819,0.05736,0.051247,0.049791,0.070005,0.067165,0.066324,0.049407,0.051484,0.057945,0.049198,0.053695,0.05499,0.027731,0.039357,0.041202,0.056464,0.045155,0.051632,0.04834,0.046081,0.052146,0.048396,0.046528,0.050707,0.041822,0.054557,0.03403,0.054347,0.053152,0.049975,0.037286,0.057094,0.039865,0.047914,0.034723,0.043145,0.04487,0.033848,0.04151,0.038491,0.045592,0.037873,0.053624,0.034402,0.051041,0.04235,0.042134,0.048766,0.049209,0.037847,0.050245,0.043494,0.037077,0.036249,0.039505,0.047284,0.051698,0.042851,0.050259,0.043868,0.033227,0.044279,0.042579,0.051765,0.0504,0.031042,0.025196,0.025308,0.034905,0.049261,0.04335,0.0256,0.00154,0.02743,0.0015636,0.0080127,0.029511,0.008185,0.00098227,0.01904,0.034083,0.0020573,0.0033655,0.024789,0.0048182,0.018955,0.00417,0.02327,0.0013686,0.027715,0.0024418,0.021173,0.0085582,0.0054255,0.0087664,0.035012,0.0012705,0.011283,0.0035868,0.00060545,0.0019945,0.019606,0.013537,0.011439,0.015092,0.0089668,0.020167,0.0071059,0.022647,0.0063291,0.026141,0.0074027,0.0061595,0.021045,0.0069668,0.037268,0.00582,0.0074573,0.0087164,0.013159,0.013808,0.00651,0.0048041,0.018389,0.0038591,0.01074,0.00055,0.012045,0.0050686,0.0092323,0.010122,0.015176,0.021215,0.020468,0.0069314,0.0037823,0.0079245,0.018128,0.0058764,0.0054068,0.0053123,0.0088423,0.0061832,0.010251,0.0020359,0.0091686,0.0061027,0.0024486,0.00446,0.013828,0.019466,0.00054318,0.012759,0.0087073,0.0096232,0.0048341,0.0055632,0.016269,0.0068759,0.013711,0.00178,0.019852,0.0036023,0.010282,0.0090368,0.0049386,0.012666,0.010677,0.025408,0.0041382,0.0026718,0.011766,0.010845,0.0040886,0.015935,0.02278,0.033032,0.012617,0.019292,0.012292,0.038385,0.00061364,0.0020791,0.010215,0.0034668,0.0076432,0.008295,0.0095364,0.0086773,0.0046368,0.01774,0.0010786,0.0028936,0.0012314,0.012231,0.012199,0.0087073,0.0202,0.0074514,0.0080468,0.0053836,0.013946,0.015972,0.0055936,0.008915,0.0064318,0.0022414,0.00085273,0.018302,0.0044595,0.00051955,0.0034023,0.0028227,0.00384,0.0024686,0.0059223,0.0052259,0.013273,0.018831,0.0016264,0.028704,0.00058591,0.0055573,0.0014445,0.0087118,0.0015945,0.0065295,0.0020773,0.015399,0.0069782,0.0037895,0.0091655,0.010773,0.010767,0.020561,0.0063432,0.0022909,0.0031591,0.0045691,0.02062,0.010078,0.0014609,0.0010736,0.0066582,0.0061118,0.018023,0.0088305,0.011383,0.020508,0.02394,0.019502,0.0060868,0.012083,0.0091595,0.01068,0.010047,0.0062827,0.018671,0.014641,0.0028945,0.0084077,0.01588,0.01336,0.018187,0.019744,0.0094136,0.012939,0.01527,0.0011182,0.0018832,0.011054,0.0010182,0.011438,0.0079823,0.01224,0.0013523,0.0040095,0.025847,0.017977,0.006585,0.00063909,0.0031441,0.016788,0.012682,0.017995,0.0036036,0.0129,0.01533,0.013085,0.0069882,0.0044827,0.0077714,0.0084827,0.002975,0.027972,0.0022255,0.011507,0.0066218,0.013657,0.013983,0.026388,0.01029,0.0105,0.0061209,0.001705,0.0039659,0.0013227,0.0017109,0.00615,0.0093268,0.0022255,0.01926,0.0018991,0.0092923,0.0081941,0.0070786,0.01267,0.0015536,0.013512,0.0032991,0.0093205,0.0061695,0.010332,0.0075609,0.0035164,0.0048555,0.0076814,0.014757,0.0096427,0.0051036,0.017701,0.011592,0.025563,0.022338,0.0083927,0.0098009,0.01033,0.00771,0.013284,0.0018286,0.0054736,0.0054536,0.017928,0.01992,0.0087214,0.0071182,0.017511,0.0015509,0.017598,0.0064586,0.00311,0.010605,0.00059273,0.0069218,0.0074418,0.00873,0.00030909,0.0010923,0.0056086,0.011361,0.0012773,0.019599,0.0021523,0.019686,0.0039982,0.0095686,0.0014727,0.010583,0.015824,0.014697,0.00039227,0.012343,0.0011873,0.0090709,0.0027359,0.0042245,0.012464,0.01676,0.0096345,0.0076177,0.00011909,0.015245,0.00278,0.010358,0.015858,0.0014855,0.00512,0.010848,0.010767,0.00084636,0.0047982,0.010958,0.0016077,0.011088,0.00025409,0.00047182,0.020509,0.017408,0.0060291,0.01605,0.010073,0.01272,0.0066095,0.016988,0.0022641,0.010855,0.0092723,0.016153,0.0051464,0.015912,0.0034632,0.016773,0.0070491,0.0011527,0.0077391,0.0019805,0.0042123,0.0010186,0.0027605,0.0079564,0.0025368,0.0057277,0.0061523,0.0069227,0.0042932,0.0019541,0.022206,0.0011868,0.019834,0.0029241,0.011688,0.0035586,0.017201,0.012014,0.013325,0.0013536,0.024058,0.0077064,0.0081218,0.029541,0.017129,0.022124,0.010307,0.0026277,0.00024591,0.0045886,0.00010455,0.00017273,0.0017505,0.027238,0.0043355,0.0135,0.01374,0.0069123,0.011926,0.0035477,0.0022332,0.0043027,0.0032027,0.0050445,0.010112,0.0065377,0.0093973,0.00098591,0.01751,0.0082905,0.020951,0.01128,0.0090918,0.0038109,0.0044595,0.0015232,0.015437,0.010456,0.021605,0.0022777,0.016125,0.012005,0.0020741,0.01723,0.0020545,0.013765,0.0087595,0.012747,0.016602,0.0029682,0.0017064,0.00124,0.0063745,0.0028932,0.0037818,0.00012318,0.00018091,0.01302,0.0081409,0.00091136,0.0076227,0.02344,0.0031318,0.0072659,0.0040686,0.00082364,0.00039409,0.0065373,0.00078727,0.00037182,0.016703,0.018689,0.02667,0.018423,0.0076114,0.0092427,0.0151,0.0010264,0.0051882,0.0010818,0.012429,0.014614,0.0036645,0.016177,0.0039782,0.015577,0.0018927,0.0093755,0.0029609,0.0018973,0.0031268,0.00024182,0.033121,0.015431,0.021066,0.0041445,0.0088659,0.030176,0.01804,0.018831,0.016052,0.014224,0.0062845,0.0073273,0.00070091,0.0023118,0.00026909,0.012896,0.00012727,0.00054273,0.0029509,0.016045,0.019399,0.0035968,0.0048745,0.0041509,0.0030768,0.0051064,0.0043441,0.0041059,0.0061836,0.0082409,0.022175,0.015223,0.018727,0.0031895,0.0044091,0.020029,0.021682,0.0049168,0.011029,0.0021868,0.0023891,0.0060241,0.002,0.0028455,0.013199,0.0019827,0.00077182,0.0048891,0.0070555,0.0084259,0.02065,0.0071955,0.012884,0.0098927,0.019839,0.022535,0.0044273,0.00041455,0.0010886,0.0024614,0.0051723,0.0020555,0.015097,0.01051,0.010556,0.010842,0.011283,0.01167,0.0057682,0.0020732,0.011699,0.00034364,0.0073595,0.0015473,0.001055,0.0020655,0.0033732,0.0015459,0.011859,0.0010264,0.00098682,0.001345,0.0015064,0.0021055,0.0011655,0.012724,0.0039786,0.0055209,0.011501,0.00090273,0.00048318,0.0049723,0.010786,0.0015945,0.002735,0.016228,0.0025136,0.00265,0.00028227,0.025653,0.03504,0.033504,0.028648,0.018798,0.014316,0.00829,0.019155,0.015975,0.022138,0.020618,0.0058941,0.017935,0.023445,0.010345,0.0055841,0.019326,0.0090677,0.010355,0.016144,0.018524,0.011085,0.0038732,0.0029,0.00069409,0.013253,0.0019005,0.022595,0.016574,0.015316,0.025496,0.013378,0.017742,0.026461,0.0054509,0.011376,0.0018086,0.0088545,0.020849,0.020029,0.0090009,0.015247,0.023265,0.0083032,0.0041686,0.015767,0.012335,0.013762,0.016532,0.017439,0.013157,0.01491,0.00086773,0.012444,0.011219,0.001805,0.015908,0.013884,0.00227,0.011088,0.015659,0.00192,0.013161,0.0067532,0.0036895,0.016607,0.013679,0.0046932,0.011296,0.0027855,0.015749,0.023142,0.019373,0.018321,0.018275,0.030474,0.012317,0.010295,0.0098109,0.01725,0.014248,0.010427,0.02366,0.015818,0.022523,0.021665,0.022829,0.019987,0.016894,0.027466,0.022343,0.018643,0.021987,0.0014105,0.022176,0.019489,0.022673,0.02393,0.016463,0.017409,0.017134,0.025649,0.019934,0.027297,0.017667,0.029554,0.028166,0.024012,0.023156,0.020318,0.01581,0.019538,0.017644,0.027121,0.007965,0.024089,0.023949,0.027649,0.019789,0.019528,0.022268,0.015798,0.023402,0.019923,0.026013,0.016442,0.023654,0.021888,0.017577,0.03056,0.018503,0.016932,0.019108,0.087007,0.04598,0.012846,0.011445,0.026462,0.0076718,0.0091732,0.032576,0.0054682,0.05251,0.0043391,0.023001,0.061091,0.048292,0.046156,0.0064205,0.081302,0.12548,0.074908,0.014907,0.014949,0.016471,0.050648,0.0252,0.027315,0.036485,0.074433,0.12361,0.086199,0.078527,0.018159,0.042174,0.11859,0.1046,0.080299,0.071212,0.07545,0.0065364,0.0033477,0.024356,0.10042,0.14566,0.086502,0.10192,0.19017,0.13317,0.12237,0.09809,0.083876,0.051683,0.10516,0.075602,0.12108,0.047885,0.04112,0.17078,0.17577,0.088196,0.094658,0.11377,0.018144,0.071217,0.050389,0.068433,0.01782,0.012496,0.09273,0.10792,0.056728,0.06982,0.055559,0.0023159,0.036642,0.074715,0.1201,0.067496,0.065204,0.17743,0.15822,0.097004,0.095425,0.099145,0.011231,0.016338,0.023389,0.047815,0.09742,0.040966,0.025334,0.075147,0.10831,0.059716,0.076941,0.051856,0.0020341,0.0049295,0.087031,0.14866,0.06544,0.067263,0.18921,0.15968,0.10012,0.10348,0.012664])
    if refpred2.ndim==1:
        feature_index=np.argmax(refpred2/baseprob)
    else:
        feature_index=np.argmax(refpred2[0,:]/baseprob)

data=np.reshape(diff[:,feature_index],[1000,3]).T




#make plots
ybcmap=matplotlib.colors.ListedColormap(np.asarray([[0,114,189],
[10,120,192],
[21,126,194],
[32,132,197],
[43,137,200],
[53,143,203],
[64,149,206],
[74,155,208],
[85,161,211],
[96,167,214],
[107,173,216],
[117,179,219],
[128,185,222],
[139,191,225],
[149,197,227],
[160,202,230],
[171,208,233],
[181,214,236],
[192,220,238],
[203,226,241],
[213,232,244],
[224,238,247],
[230,241,248],
[235,244,250],
[241,247,251],
[246,250,253],
[252,253,254],
[252,253,254],
[253,254,254],
[254,254,254],
[254,254,255],
[255,255,255],
[255,255,253],
[255,255,252],
[255,255,250],
[255,255,249],
[255,255,248],
[255,255,244],
[255,255,240],
[255,255,237],
[255,255,233],
[255,255,229],
[255,255,219],
[255,255,209],
[255,255,198],
[255,255,188],
[255,255,177],
[255,255,167],
[255,255,157],
[255,255,146],
[255,255,136],
[255,255,125],
[255,255,115],
[255,255,105],
[255,255,94],
[255,255,84],
[255,255,73],
[255,255,63],
[255,255,53],
[255,255,42],
[255,255,32],
[255,255,21],
[255,255,11],
[255,255,0]])/255.0)


column_labels=list(fasta_sequences[0])

fig, ax = plt.subplots(figsize=(200,5))
heatmap = ax.pcolor(data,vmax=np.abs(data).max(),vmin=-np.abs(data).max(), cmap=ybcmap)

ax.set_xticks(np.arange(1,data.shape[1]+1)-0.5, minor=False)
ax.set_yticks(np.arange(data.shape[0]), minor=False)


ax.xaxis.tick_top()
ax.set_xticklabels(column_labels, minor=False)
ax.set_yticklabels(['', '', ''], minor=False)


ax2 = ax.twiny()
ax2.set_xticks(np.arange(100,data.shape[1]+1,100)-1.5, minor=False)
ax2.set_xticklabels(np.arange(100,data.shape[1]+1,100), minor=False)



savefig(sys.argv[1]+'.vis.png', transparent=True, bbox_inches='tight', pad_inches=0.2)

ax.set_visible(False)
ax2.set_visible(False)

fig.set_size_inches(10,2)

cax = plt.axes([0, 0, 1, 0.1])
feature_names=pd.read_csv('./resources/feature_name',header=None,sep='\t')
cax.set_title('Log2 fold change caused by mutation')
fig.suptitle(feature_names.iloc[feature_index+1,1], fontsize=20)

plt.colorbar(heatmap, cax=cax, orientation='horizontal')
savefig(sys.argv[1]+'.colorbar.png', transparent=True, bbox_inches='tight', pad_inches=0.1)


#preview
fig, ax = plt.subplots(figsize=(10,2))
heatmap = ax.pcolor(data,vmax=np.abs(data).max(),vmin=-np.abs(data).max(), cmap=ybcmap)
ax.set_yticks(np.arange(data.shape[0]), minor=False)
ax.xaxis.tick_top()
ax.set_yticklabels(['', '', ''], minor=False)
ax.set_xticks(np.arange(100,data.shape[1]+1,100)-1.5, minor=False)
ax.set_xticklabels(np.arange(100,data.shape[1]+1,100), minor=False)
savefig(sys.argv[1]+'.preview.png', transparent=True, bbox_inches='tight', pad_inches=0.1)




data = pd.DataFrame(data)
data.columns = np.arange(1,1001)
data.rows= ['Mutation 1','Mutation 2', 'Mutation 3']
data.to_csv(sys.argv[1]+'.csv')

