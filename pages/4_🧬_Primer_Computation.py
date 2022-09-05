from Bio.Seq import Seq
import numpy
import math
from Bio.SeqUtils import GC
import streamlit as st

st.title("PCR primer computer")
st.markdown('This app computes the forward and reverse primer. It also calculates the melting temperature of both primers. How the app works and its limitations are documented [here](https://www.dna-vorm.com/bio-tools).')
st.markdown('Based on the paper: [*A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics*](https://www.pnas.org/doi/pdf/10.1073/pnas.95.4.1460).')

my_seq=Seq(st.text_input("Your target sequence"))
my_seq=my_seq.upper()
allowed="ACTG"
test=all(ch in allowed for ch in my_seq)
if test==True:
    pass
else:
    st.write("Warning: Other characters beside A, C, T, or Gs were found.")
C_primer=st.number_input("Concentration of the primers in uM (for example 0.5)",min_value=0.1,max_value=5., value=0.5, step=0.01)
C_salt=st.number_input("Concentration of the NaCl in M (for example 0.5)",min_value=0.1,max_value=5., value=0.5, step=0.01)
primer_length=int(st.number_input("The desired length of the primers",min_value=1,max_value=400, value=20, step=1))


Forw_primer=my_seq[0:primer_length]
rev_primer_temp=my_seq[(-primer_length):]
rev_primer=rev_primer_temp.reverse_complement()
st.subheader("Forward primer: "+Forw_primer)
st.subheader("Reverse primer: "+rev_primer)

seq_length=len(my_seq)
dH_fwd=0.1
dS_fwd=-2.8
R=1.99
for i in range(0,len(Forw_primer)):
    if i == (len(Forw_primer)-1):
        pass
    else:
        if Forw_primer[i]=='A':
            if Forw_primer[i+1]=='A':
                dH_fwd+= -7.9
                dS_fwd+= -22.2
            elif Forw_primer[i+1]=='C':
                dH_fwd+= -8.4
                dS_fwd+= -22.4
            elif Forw_primer[i+1]=='T':
                dH_fwd+= -7.2
                dS_fwd+= -20.4
            else:
                dH_fwd+= -7.8
                dS_fwd+= -21
        elif Forw_primer[i]=='C':
            if Forw_primer[i+1]=='A':
                dH_fwd+= -8.5
                dS_fwd+= -22.7
            elif Forw_primer[i+1]=='C':
                dH_fwd+= -8
                dS_fwd+= -19.9
            elif Forw_primer[i+1]=='T':
                dH_fwd+= -7.8
                dS_fwd+= -21
            else:
                dH_fwd+= -10.6
                dS_fwd+= -27.2
        elif Forw_primer[i]=='T':
            if Forw_primer[i+1]=='A':
                dH_fwd+= -7.2
                dS_fwd+= -21.3
            elif Forw_primer[i+1]=='C':
                dH_fwd+= -8.2
                dS_fwd+= -22.2
            elif Forw_primer[i+1]=='T':
                dH_fwd+= -7.9
                dS_fwd+= -22.2
            else:
                dH_fwd+= -8.5
                dS_fwd+= -22.7
        else:
            if Forw_primer[i+1]=='A':
                dH_fwd+= -8.2
                dS_fwd+= -22.2
            elif Forw_primer[i+1]=='C':
                dH_fwd+= -9.8
                dS_fwd+= -24.4
            elif Forw_primer[i+1]=='T':
                dH_fwd+= -8.4
                dS_fwd+= -22.4
            else:
                dH_fwd+= -8
                dS_fwd+= -19.9
dH_fwd=dH_fwd*1000

dH_rev=0.1
dS_rev=-2.8
for i in range(0,len(rev_primer)):
    if i == (len(rev_primer)-1):
        pass
    else:
        if rev_primer[i]=='A':
            if rev_primer[i+1]=='A':
                dH_rev+= -7.9
                dS_rev+= -22.2
            elif rev_primer[i+1]=='C':
                dH_rev+= -8.4
                dS_rev+= -22.4
            elif rev_primer[i+1]=='T':
                dH_rev+= -7.2
                dS_rev+= -20.4
            else:
                dH_rev+= -7.8
                dS_rev+= -21
        elif rev_primer[i]=='C':
            if rev_primer[i+1]=='A':
                dH_rev+= -8.5
                dS_rev+= -22.7
            elif rev_primer[i+1]=='C':
                dH_rev+= -8
                dS_rev+= -19.9
            elif rev_primer[i+1]=='T':
                dH_rev+= -7.8
                dS_rev+= -21
            else:
                dH_rev+= -10.6
                dS_rev+= -27.2
        elif rev_primer[i]=='T':
            if rev_primer[i+1]=='A':
                dH_rev+= -7.2
                dS_rev+= -21.3
            elif rev_primer[i+1]=='C':
                dH_rev+= -8.2
                dS_rev+= -22.2
            elif rev_primer[i+1]=='T':
                dH_rev+= -7.9
                dS_rev+= -22.2
            else:
                dH_rev+= -8.5
                dS_rev+= -22.7
        else:
            if rev_primer[i+1]=='A':
                dH_rev+= -8.2
                dS_rev+= -22.2
            elif rev_primer[i+1]=='C':
                dH_rev+= -9.8
                dS_rev+= -24.4
            elif rev_primer[i+1]=='T':
                dH_rev+= -8.4
                dS_rev+= -22.4
            else:
                dH_rev+= -8
                dS_rev+= -19.9
dH_rev=dH_rev*1000

c_term=R*math.log((C_primer*10**-6)/4)

Tm_fwd=(dH_fwd/(dS_fwd+c_term))-273.15+(16.6*math.log10(C_salt))
Tm_rev=(dH_rev/(dS_rev+c_term))-273.15+(16.6*math.log10(C_salt))
dif_Tm=Tm_fwd-Tm_rev
st.subheader("Melting temperature forward primer: "+str(round(Tm_fwd,1)))
st.subheader("Melting temperature reverse primer: "+str(round(Tm_rev,1)))
if -5>dif_Tm or dif_Tm >5:
    st.write("Warning! The difference in melting temperature is more than 5 degrees celsius")
else:
    pass

