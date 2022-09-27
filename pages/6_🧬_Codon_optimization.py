import numpy as np
import pandas as pd
from Bio.Seq import Seq
import math
from Bio.SeqUtils import GC
import streamlit as st
import os
from io import BytesIO
from pyxlsb import open_workbook as open_xlsb


def codon_optimizers(organism,my_seq_AA):

    df = pd.read_csv(f"pages/Files/{organism}_codon.txt", header=None, delim_whitespace=True)
    df1=df.iloc[:, : 5]
    df2=df.iloc[:, 5:10]
    df3=df.iloc[:, 10:15]
    df4=df.iloc[:, 15:]
    df2.columns = df1.columns
    df3.columns = df1.columns
    df4.columns = df1.columns
    df_combined = pd.concat([df1, df2, df3,df4], axis=0, ignore_index=True)
    df=df_combined.iloc[:, : 3]
    df[2].values[df[2] <= 0.1 ] = 0
    AA_sum=[]
    for i in df[1].unique():
        AA_sum.append(df.loc[df[1] == i, 2].sum())
    test=df[1].unique()
    sumlist=test.tolist()
    data= {'AA':sumlist,'sum':AA_sum}
    df2=pd.DataFrame(data)
    df.columns = ['Codon', 'AA', 'perc']
    df=df.reset_index(drop=True, inplace=False)
    df=pd.merge(df, df2)
    df['real_value']=df['perc']/df['sum']
    dfs = dict(tuple(df.groupby('AA')))
    AA_F = dfs['F']
    AA_F=AA_F.sort_values('real_value', ascending=False)
    AA_F=AA_F.reset_index(drop=True, inplace=False)
    AA_L = dfs['L']
    AA_L=AA_L.sort_values('real_value', ascending=False)
    AA_L=AA_L.reset_index(drop=True, inplace=False)
    AA_I = dfs['I']
    AA_I=AA_I.sort_values('real_value', ascending=False)
    AA_I=AA_I.reset_index(drop=True, inplace=False)
    AA_M = dfs['M']
    AA_M=AA_M.sort_values('real_value', ascending=False)
    AA_M=AA_M.reset_index(drop=True, inplace=False)
    AA_V = dfs['V']
    AA_V=AA_V.sort_values('real_value', ascending=False)
    AA_V=AA_V.reset_index(drop=True, inplace=False)
    AA_Y = dfs['Y']
    AA_Y=AA_Y.sort_values('real_value', ascending=False)
    AA_Y=AA_Y.reset_index(drop=True, inplace=False)
    AA_P = dfs['P']
    AA_P=AA_P.sort_values('real_value', ascending=False)
    AA_P=AA_P.reset_index(drop=True, inplace=False)
    AA_T = dfs['T']
    AA_T=AA_T.sort_values('real_value', ascending=False)
    AA_T=AA_T.reset_index(drop=True, inplace=False)
    AA_A = dfs['A']
    AA_A=AA_A.sort_values('real_value', ascending=False)
    AA_A=AA_A.reset_index(drop=True, inplace=False)
    AA_AS = dfs['*']
    AA_AS=AA_AS.sort_values('real_value', ascending=False)
    AA_AS=AA_AS.reset_index(drop=True, inplace=False)
    AA_H = dfs['H']
    AA_H=AA_H.sort_values('real_value', ascending=False)
    AA_H=AA_H.reset_index(drop=True, inplace=False)
    AA_Q = dfs['Q']
    AA_Q=AA_Q.sort_values('real_value', ascending=False)
    AA_Q=AA_Q.reset_index(drop=True, inplace=False)
    AA_N = dfs['N']
    AA_N=AA_N.sort_values('real_value', ascending=False)
    AA_N=AA_N.reset_index(drop=True, inplace=False)
    AA_K = dfs['K']
    AA_K=AA_K.sort_values('real_value', ascending=False)
    AA_K=AA_K.reset_index(drop=True, inplace=False)
    AA_D = dfs['D']
    AA_D=AA_D.sort_values('real_value', ascending=False)
    AA_D=AA_D.reset_index(drop=True, inplace=False)
    AA_E = dfs['E']
    AA_E=AA_E.sort_values('real_value', ascending=False)
    AA_E=AA_E.reset_index(drop=True, inplace=False)
    AA_C = dfs['C']
    AA_C=AA_C.sort_values('real_value', ascending=False)
    AA_C=AA_C.reset_index(drop=True, inplace=False)
    AA_W = dfs['W']
    AA_W=AA_W.sort_values('real_value', ascending=False)
    AA_W=AA_W.reset_index(drop=True, inplace=False)
    AA_R = dfs['R']
    AA_R=AA_R.sort_values('real_value', ascending=False)
    AA_R=AA_R.reset_index(drop=True, inplace=False)
    AA_S = dfs['S']
    AA_S=AA_S.sort_values('real_value', ascending=False)
    AA_S=AA_S.reset_index(drop=True, inplace=False)
    AA_G = dfs['G']
    AA_G=AA_G.sort_values('real_value', ascending=False)
    AA_G=AA_G.reset_index(drop=True, inplace=False)



    Codon_opt_seq=Seq('')
    for i in my_seq_AA:
        if i=='M':
            x=np.random.random(1)[0]
            if x<=AA_M.real_value[0]:
                Codon_opt_seq+=Seq(AA_M.Codon[0])
            else:
                pass
        elif i=='A':
            x=np.random.random(1)[0]
            if x<=AA_A.real_value[0]:
                Codon_opt_seq+=Seq(AA_A.Codon[0])
            elif x>AA_A.real_value[0] and x<=(AA_A.real_value[0]+AA_A.real_value[1]):
                Codon_opt_seq+=Seq(AA_A.Codon[1])
            elif x>(AA_A.real_value[0]+AA_A.real_value[1]) and x<=(AA_A.real_value[0]+AA_A.real_value[1]+AA_A.real_value[2]):
                Codon_opt_seq+=Seq(AA_A.Codon[2])
            else:
                Codon_opt_seq+=Seq(AA_A.Codon[3])
        elif i=='F':
            x=np.random.random(1)[0]
            if x<=AA_F.real_value[0]:
                Codon_opt_seq+=Seq(AA_F.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_F.Codon[1])
        elif i=='L':
            x=np.random.random(1)[0]
            if x<=AA_L.real_value[0]:
                Codon_opt_seq+=Seq(AA_L.Codon[0])
            elif x>AA_L.real_value[0] and x<=(AA_L.real_value[0]+AA_L.real_value[1]):
                Codon_opt_seq+=Seq(AA_L.Codon[1])
            elif x>(AA_L.real_value[0]+AA_L.real_value[1]) and x<=(AA_L.real_value[0]+AA_L.real_value[1]+AA_L.real_value[2]):
                Codon_opt_seq+=Seq(AA_L.Codon[2])
            elif x>(AA_L.real_value[0]+AA_L.real_value[1]+AA_L.real_value[2]) and x<=(AA_L.real_value[0]+AA_L.real_value[1]+AA_L.real_value[2]+AA_L.real_value[3]):
                Codon_opt_seq+=Seq(AA_L.Codon[3])
            elif x>(AA_L.real_value[0]+AA_L.real_value[1]+AA_L.real_value[2]+AA_L.real_value[3]) and x<=(AA_L.real_value[0]+AA_L.real_value[1]+AA_L.real_value[2]+AA_L.real_value[3]+AA_L.real_value[4]):
                Codon_opt_seq+=Seq(AA_L.Codon[4])
            else:
                Codon_opt_seq+=Seq(AA_L.Codon[5])
        elif i=='I':
            x=np.random.random(1)[0]
            if x<=AA_I.real_value[0]:
                Codon_opt_seq+=Seq(AA_I.Codon[0])
            elif x>AA_I.real_value[0] and x<=(AA_I.real_value[0]+AA_I.real_value[1]):
                Codon_opt_seq+=Seq(AA_I.Codon[1])
            else:
                Codon_opt_seq+=Seq(AA_I.Codon[2])
        elif i=='V':
            x=np.random.random(1)[0]
            if x<=AA_V.real_value[0]:
                Codon_opt_seq+=Seq(AA_V.Codon[0])
            elif x>AA_V.real_value[0] and x<=(AA_V.real_value[0]+AA_V.real_value[1]):
                Codon_opt_seq+=Seq(AA_V.Codon[1])
            elif x>(AA_V.real_value[0]+AA_V.real_value[1]) and x<=(AA_V.real_value[0]+AA_V.real_value[1]+AA_V.real_value[2]):
                Codon_opt_seq+=Seq(AA_V.Codon[2])
            else:
                Codon_opt_seq+=Seq(AA_V.Codon[3])
        elif i=='S':
            x=np.random.random(1)[0]
            if x<=AA_S.real_value[0]:
                Codon_opt_seq+=Seq(AA_S.Codon[0])
            elif x>AA_S.real_value[0] and x<=(AA_S.real_value[0]+AA_S.real_value[1]):
                Codon_opt_seq+=Seq(AA_S.Codon[1])
            elif x>(AA_S.real_value[0]+AA_S.real_value[1]) and x<=(AA_S.real_value[0]+AA_S.real_value[1]+AA_S.real_value[2]):
                Codon_opt_seq+=Seq(AA_S.Codon[2])
            elif x>(AA_S.real_value[0]+AA_S.real_value[1]+AA_S.real_value[2]) and x<=(AA_S.real_value[0]+AA_S.real_value[1]+AA_S.real_value[2]+AA_S.real_value[3]):
                Codon_opt_seq+=Seq(AA_S.Codon[3])
            elif x>(AA_S.real_value[0]+AA_S.real_value[1]+AA_S.real_value[2]+AA_S.real_value[3]) and x<=(AA_S.real_value[0]+AA_S.real_value[1]+AA_S.real_value[2]+AA_S.real_value[3]+AA_S.real_value[4]):
                Codon_opt_seq+=Seq(AA_S.Codon[4])
            else:
                Codon_opt_seq+=Seq(AA_S.Codon[5])
        elif i=='P':
            x=np.random.random(1)[0]
            if x<=AA_P.real_value[0]:
                Codon_opt_seq+=Seq(AA_P.Codon[0])
            elif x>AA_P.real_value[0] and x<=(AA_P.real_value[0]+AA_P.real_value[1]):
                Codon_opt_seq+=Seq(AA_P.Codon[1])
            elif x>(AA_P.real_value[0]+AA_P.real_value[1]) and x<=(AA_P.real_value[0]+AA_P.real_value[1]+AA_P.real_value[2]):
                Codon_opt_seq+=Seq(AA_P.Codon[2])
            else:
                Codon_opt_seq+=Seq(AA_P.Codon[3])
        elif i=='T':
            x=np.random.random(1)[0]
            if x<=AA_T.real_value[0]:
                Codon_opt_seq+=Seq(AA_T.Codon[0])
            elif x>AA_T.real_value[0] and x<=(AA_T.real_value[0]+AA_T.real_value[1]):
                Codon_opt_seq+=Seq(AA_T.Codon[1])
            elif x>(AA_T.real_value[0]+AA_T.real_value[1]) and x<=(AA_T.real_value[0]+AA_T.real_value[1]+AA_T.real_value[2]):
                Codon_opt_seq+=Seq(AA_T.Codon[2])
            else:
                Codon_opt_seq+=Seq(AA_T.Codon[3])
        elif i=='Y':
            x=np.random.random(1)[0]
            if x<=AA_Y.real_value[0]:
                Codon_opt_seq+=Seq(AA_Y.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_Y.Codon[1])
        elif i=='*':
            x=np.random.random(1)[0]
            if x<=AA_AS.real_value[0]:
                Codon_opt_seq+=Seq(AA_AS.Codon[0])
            elif x>AA_AS.real_value[0] and x<=(AA_AS.real_value[0]+AA_AS.real_value[1]):
                Codon_opt_seq+=Seq(AA_AS.Codon[1])
            else:
                Codon_opt_seq+=Seq(AA_AS.Codon[2])
        elif i=='H':
            x=np.random.random(1)[0]
            if x<=AA_H.real_value[0]:
                Codon_opt_seq+=Seq(AA_H.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_H.Codon[1])
        elif i=='Q':
            x=np.random.random(1)[0]
            if x<=AA_Q.real_value[0]:
                Codon_opt_seq+=Seq(AA_Q.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_Q.Codon[1])
        elif i=='N':
            x=np.random.random(1)[0]
            if x<=AA_N.real_value[0]:
                Codon_opt_seq+=Seq(AA_N.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_N.Codon[1])
        elif i=='K':
            x=np.random.random(1)[0]
            if x<=AA_K.real_value[0]:
                Codon_opt_seq+=Seq(AA_K.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_K.Codon[1])
        elif i=='D':
            x=np.random.random(1)[0]
            if x<=AA_D.real_value[0]:
                Codon_opt_seq+=Seq(AA_D.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_D.Codon[1])
        elif i=='E':
            x=np.random.random(1)[0]
            if x<=AA_E.real_value[0]:
                Codon_opt_seq+=Seq(AA_E.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_E.Codon[1])
        elif i=='C':
            x=np.random.random(1)[0]
            if x<=AA_C.real_value[0]:
                Codon_opt_seq+=Seq(AA_C.Codon[0])
            else:
                Codon_opt_seq+=Seq(AA_C.Codon[1])
        if i=='W':
            x=np.random.random(1)[0]
            if x<=AA_W.real_value[0]:
                Codon_opt_seq+=Seq(AA_W.Codon[0])
            else:
                pass
        elif i=='R':
            x=np.random.random(1)[0]
            if x<=AA_R.real_value[0]:
                Codon_opt_seq+=Seq(AA_R.Codon[0])
            elif x>AA_R.real_value[0] and x<=(AA_R.real_value[0]+AA_R.real_value[1]):
                Codon_opt_seq+=Seq(AA_R.Codon[1])
            elif x>(AA_R.real_value[0]+AA_R.real_value[1]) and x<=(AA_R.real_value[0]+AA_R.real_value[1]+AA_R.real_value[2]):
                Codon_opt_seq+=Seq(AA_R.Codon[2])
            elif x>(AA_R.real_value[0]+AA_R.real_value[1]+AA_R.real_value[2]) and x<=(AA_R.real_value[0]+AA_R.real_value[1]+AA_R.real_value[2]+AA_R.real_value[3]):
                Codon_opt_seq+=Seq(AA_R.Codon[3])
            elif x>(AA_R.real_value[0]+AA_R.real_value[1]+AA_R.real_value[2]+AA_R.real_value[3]) and x<=(AA_R.real_value[0]+AA_R.real_value[1]+AA_R.real_value[2]+AA_R.real_value[3]+AA_R.real_value[4]):
                Codon_opt_seq+=Seq(AA_R.Codon[4])
            else:
                Codon_opt_seq+=Seq(AA_R.Codon[5])
        elif i=='G':
            x=np.random.random(1)[0]
            if x<=AA_G.real_value[0]:
                Codon_opt_seq+=Seq(AA_G.Codon[0])
            elif x>AA_G.real_value[0] and x<=(AA_G.real_value[0]+AA_G.real_value[1]):
                Codon_opt_seq+=Seq(AA_G.Codon[1])
            elif x>(AA_G.real_value[0]+AA_G.real_value[1]) and x<=(AA_G.real_value[0]+AA_G.real_value[1]+AA_G.real_value[2]):
                Codon_opt_seq+=Seq(AA_G.Codon[2])
            else:
                Codon_opt_seq+=Seq(AA_G.Codon[3])
        else:
            pass

    return Codon_opt_seq
def to_excel(df):
    output = BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    df.to_excel(writer, index=False, sheet_name='Sheet1')
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    format1 = workbook.add_format({'num_format': '0.00'})
    worksheet.set_column('A:A', None, format1)
    writer.save()
    processed_data = output.getvalue()
    return processed_data
st.title("Codon optimizer")
st.markdown('This app optimizes an RNA or DNA sequence for a certain organism. How the app works and its limitations are documented [here](https://www.dna-vorm.com/bio-tools).')

organism_files=os.listdir("pages/Files")
organism_list = [i.split('_')[0] for i in organism_files]
organism_list.remove('.ipynb')

organism=st.selectbox('Select the organism you want to optimize for',organism_list)

input_type = st.radio("What is the sequence type of your input?",('DNA/RNA', 'Amino acids'))
number_output=st.number_input("How many sequences do you want as output?", min_value=1,max_value=100)
if input_type=='DNA/RNA':
    my_seq = Seq(st.text_input("Your sequence"))
    my_seq = my_seq.upper()
    my_seq_AA = my_seq.translate()
    allowed = "ACTUG"
    test = all(ch in allowed for ch in my_seq)
    if test == True:
        pass
    else:
        st.write("Warning: Other characters beside A, C, T, U, or Gs were found.")

else:
    my_seq_AA = Seq(st.text_input("Your sequence"))
    my_seq_AA = my_seq_AA.upper()
if number_output<6:
    for i in range(number_output):
        Codon_opt_seq=codon_optimizers(organism,my_seq_AA)
        st.markdown(f'Optimized sequence #{i+1}: 5\'{Codon_opt_seq} 3\'')
else:
    data=[]
    for i in range(number_output):
        Codon_opt_seq = codon_optimizers(organism, my_seq_AA)
        data.append(Codon_opt_seq)
    df=pd.DataFrame(data)
    df=df.reset_index()
    df_xlsx = to_excel(df)
    st.download_button(label="Download optimized sequences as excel file",data=df_xlsx,file_name= 'df_test.xlsx')