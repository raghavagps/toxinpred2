##############################################################################
#ToxinPred2.0 is developed for predicting toxin and non toxin      #
#protein from their primary sequence. It is developed by Prof G. P. S.       #
#Raghava's group. Please cite : ToxinPred 2.0                                  #
# ############################################################################
import argparse  
import warnings
import pickle
import os
import re
import sys
import numpy as np
import pandas as pd
import joblib

warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Please provide following arguments. Please make the suitable changes in the envfile provided in the folder.') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.5")
parser.add_argument("-m","--model",type=int, choices = [1, 2], help="Model: 1: AAC based RF, 2: Hybrid, by default 1")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Toxin, 2: All peptides, by default 1")
args = parser.parse_args()

def aac_comp(file,out):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    f = open(out, 'w')
    sys.stdout = f
    df = pd.DataFrame(file)
    zz = df.iloc[:,0]
    for j in zz:
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            print("%.2f"%composition, end = ",")
        print("")
    f.truncate()

def prediction(inputfile,model,out):
    df = pd.DataFrame()
    a=[]
    file_name = inputfile
    file_name1 = out
    file_name2 = model
    clf = joblib.load(file_name2)
    data_test = np.loadtxt(file_name, delimiter=',')
    X_test = data_test
    y_p_score1=clf.predict_proba(X_test)
    y_p_s1=y_p_score1.tolist()
    df = pd.DataFrame(y_p_s1)
    df_1 = df.iloc[:,-1]
    df_1.to_csv(file_name1, index=None, header=False)

def class_assignment(file1,thr,out):
    df1 = pd.read_csv(file1, header=None)
    df1.columns = ['ML Score']
    cc = []
    for i in range(0,len(df1)):
        if df1['ML Score'][i]>=float(thr):
            cc.append('Toxin')
        else:
            cc.append('Non-Toxin')
    df1['Prediction'] = cc
    df1 =  df1.round(3)
    df1.to_csv(out, index=None)

def MERCI_Processor(merci_file,merci_processed,name):
    hh =[]
    jj = []
    kk = []
    qq = []
    filename = merci_file
    df = pd.DataFrame(name)
    zz = list(df[0])
    check = '>'
    with open(filename) as f:
        l = []
        for line in f:
            if not len(line.strip()) == 0 :
                l.append(line)
            if 'COVERAGE' in line:
                for item in l:
                    if item.lower().startswith(check.lower()):
                        hh.append(item)
                l = []
    if hh == []:
        ff = [w.replace('>', '') for w in zz]
        for a in ff:
            jj.append(a)
            qq.append(np.array(['0']))
            kk.append('Non-Toxin')
    else:
        ff = [w.replace('\n', '') for w in hh]
        ee = [w.replace('>', '') for w in ff]
        rr = [w.replace('>', '') for w in zz]
        ff = ee + rr
        oo = np.unique(ff)
        df1 = pd.DataFrame(list(map(lambda x:x.strip(),l))[1:])
        df1.columns = ['Name']
        df1['Name'] = df1['Name'].str.strip('(')
        df1[['Seq','Hits']] = df1.Name.str.split("(",expand=True)
        df2 = df1[['Seq','Hits']]
        df2.replace(to_replace=r"\)", value='', regex=True, inplace=True)
        df2.replace(to_replace=r'motifs match', value='', regex=True, inplace=True)
        df2.replace(to_replace=r' $', value='', regex=True,inplace=True)
        total_hit = int(df2.loc[len(df2)-1]['Seq'].split()[0])
        for j in oo:
            if j in df2.Seq.values:
                jj.append(j)
                qq.append(df2.loc[df2.Seq == j]['Hits'].values)
                kk.append('Toxin')
            else:
                jj.append(j)
                qq.append(np.array(['0']))
                kk.append('Non-Toxin')
    df3 = pd.concat([pd.DataFrame(jj),pd.DataFrame(qq),pd.DataFrame(kk)], axis=1)
    df3.columns = ['Name','Hits','Prediction']
    df3.to_csv(merci_processed,index=None)

def Merci_after_processing(merci_processed,final_merci):
    df5 = pd.read_csv(merci_processed)
    df5 = df5[['Name','Hits']]
    df5.columns = ['Subject','Hits']
    kk = []
    for i in range(0,len(df5)):
        if df5['Hits'][i] > 0:
            kk.append(0.5)
        else:
            kk.append(0)
    df5["MERCI Score"] = kk
    df5 = df5[['Subject','MERCI Score']]
    df5.to_csv(final_merci, index=None)

def BLAST_processor(blast_result,blast_processed,name1):
    if os.stat(blast_result).st_size != 0:
        df1 = pd.read_csv(blast_result, sep="\t",header=None)
        df2 = df1.iloc[:,:2]
        df2.columns = ['Subject','Query']
        df3 = pd.DataFrame()
        for i in df2.Subject.unique():
            df3 = df3.append(df2.loc[df2.Subject==i][0:5]).reset_index(drop=True)
        cc= []
        for i in range(0,len(df3)):
            cc.append(df3['Query'][i].split("_")[0])
        df3['label'] = cc
        dd = []
        for i in range(0,len(df3)):
            if df3['label'][i] == 'P':
                dd.append(1)
            else:
                dd.append(-1)
        df3["vote"] = dd
        ff = []
        gg = []
        for i in df3.Subject.unique():
            ff.append(i)
            gg.append(df3.loc[df3.Subject==i]["vote"].sum())
        df4 = pd.concat([pd.DataFrame(ff),pd.DataFrame(gg)],axis=1)
        df4.columns = ['Subject','Blast_value']
        hh = []
        for i in range(0,len(df4)):
            if df4['Blast_value'][i] >0:
                hh.append(0.5)
            elif df4['Blast_value'][i] == 0:
                hh.append(0)
            else:
                hh.append(-0.5)
        df4['BLAST Score'] = hh
        df4 = df4[['Subject','BLAST Score']]
    else:
        ss = []
        vv = []
        for j in seqid:
            ss.append(j)
            vv.append(0)
        df4 = pd.concat([pd.DataFrame(ss),pd.DataFrame(vv)],axis=1)
        df4.columns = ['Subject','BLAST Score']
    df4.to_csv(blast_processed, index=None)

def hybrid(ML_output,name1,merci_output,blast_output,threshold,final_output):
    df6_2 = pd.read_csv(ML_output,header=None)
    df6_1 = pd.DataFrame(name1)
    df5 = pd.read_csv(merci_output)
    df4 = pd.read_csv(blast_output)
    df6 = pd.concat([df6_1,df6_2],axis=1)
    df6.columns = ['Subject','ML Score']
    df6['Subject'] = df6['Subject'].str.replace('>','')
    df7 = pd.merge(df6,df5, how='outer',on='Subject')
    df8 = pd.merge(df7,df4, how='outer',on='Subject')
    df8.fillna(0, inplace=True)
    df8['Hybrid Score'] = df8.sum(axis=1)
    df8 = df8.round(3)
    ee = []
    for i in range(0,len(df8)):
        if df8['Hybrid Score'][i] > float(threshold):
            ee.append('Toxin')
        else:
            ee.append('Non-Toxin')
    df8['Prediction'] = ee
    df8.to_csv(final_output, index=None)

print('##############################################################################')
print('# The program ToxinPred2.0 is developed for predicting Toxin and non toxin #')
print("# protein from their primary sequence, developed by Prof G. P. S. Raghava's group. #")
print('# ############################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.6
else:
        Threshold= float(args.threshold)
# Model
if args.model == None:
        Model = int(1)
else:
        Model = int(args.model)
# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)

print('Summary of Parameters:')
print('Input File: ',Sequence,'; Model: ',Model,'; Threshold: ', Threshold)
print('Output File: ',result_filename,'; Display: ',dplay)

#------------------ Read input file ---------------------
f=open(Sequence,"r")
len1 = f.read().count('>')
f.close()

with open(Sequence) as f:
        records = f.read()
records = records.split('>')[1:]
seqid = []
seq = []
for fasta in records:
    array = fasta.split('\n')
    name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '', ''.join(array[1:]).upper())
    seqid.append(name)
    seq.append(sequence)
if len(seqid) == 0:
    f=open(Sequence,"r")
    data1 = f.readlines()
    for each in data1:
        seq.append(each.replace('\n',''))
    for i in range (1,len(seq)+1):
        seqid.append("Seq_"+str(i))

seqid_1 = list(map(">{}".format, seqid))
CM = pd.concat([pd.DataFrame(seqid_1),pd.DataFrame(seq)],axis=1)
CM.to_csv("Sequence_1",header=False,index=None,sep="\n")
f.close()
#======================= Prediction Module start from here =====================
if Model==1:
    aac_comp(seq,'seq.aac')
    os.system("perl -pi -e 's/,$//g' seq.aac")
    prediction('seq.aac','RF_model','seq.pred')
    class_assignment('seq.pred',Threshold,'seq.out')
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    df3 = pd.read_csv("seq.out")
    df3 = round(df3,3)
    df4 = pd.concat([df1,df2,df3],axis=1)
    df4.columns = ['ID','Sequence','ML_Score','Prediction']
    if dplay == 1:
        df4 = df4.loc[df4.Prediction=="Toxin"]
    else:
        df4 = df4
    df4.to_csv(result_filename, index=None)
    os.remove('seq.aac')
    os.remove('seq.pred')
    os.remove('seq.out')
else:
    if os.path.exists('envfile'):
        with open('envfile', 'r') as file:
            data = file.readlines()
        output = []
        for line in data:
            if not "#" in line:
                output.append(line)
        if len(output)==4: 
            paths = []
            for i in range (0,len(output)):
                paths.append(output[i].split(':')[1].replace('\n',''))
            blastp = paths[0]
            blastdb = paths[1]
            merci = paths[2]
            motifs = paths[3]
        else:
            print("####################################################################################")
            print("Error: Please provide paths for BLAST, MERCI and required files", file=sys.stderr)
            print("####################################################################################")
            sys.exit()
 
    else:
        print("####################################################################################")
        print("Error: Please provide the '{}', which comprises paths for BLAST and MERCI".format('envfile'), file=sys.stderr)
        print("####################################################################################")
        sys.exit()
    aac_comp(seq,'seq.aac')
    os.system("perl -pi -e 's/,$//g' seq.aac")
    prediction('seq.aac','RF_model','seq.pred')
    os.system(blastp + " -task blastp -db " + blastdb + " -query " + "Sequence_1"  + " -out RES_1_6_6.out -outfmt 6 -evalue 0.000001")
    os.system(merci + " -p " + "Sequence_1" +  " -i " + motifs + " -o merci.txt")
    MERCI_Processor('merci.txt','merci_output.csv',seqid)
    Merci_after_processing('merci_output.csv','merci_hybrid.csv')
    BLAST_processor('RES_1_6_6.out','blast_hybrid.csv',seqid)
    hybrid('seq.pred',seqid,'merci_hybrid.csv','blast_hybrid.csv',Threshold,'final_output')
    df44 = pd.read_csv('final_output')
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="Toxin"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove('seq.aac')
    os.remove('seq.pred')
    os.remove('final_output')
    os.remove('RES_1_6_6.out')
    os.remove('merci_output.csv')
    os.remove('merci_hybrid.csv')
    os.remove('blast_hybrid.csv')
    os.remove('merci.txt')
    os.remove('Sequence_1')

print('\n======= Thanks for using ToxinPred2.0. Your results are stored in file :',result_filename,' =====\n\n')
print('Please cite: ToxinPred2.0\n')
