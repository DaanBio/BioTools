from Bio.Seq import Seq
import numpy
import streamlit as st

st.title("Sequence converter")
st.markdown('This app converts a RNA or DNA sequence into the DNA, RNA, or Amino acids. How the app works and its limitations are documented [here](https://www.dna-vorm.com/bio-tools).')


input_type = st.radio("What is the sequence type of your input?",('DNA', 'RNA'))
my_seq=Seq(st.text_input("Your sequence"))
my_seq=my_seq.upper()
allowed="ACTG"
test=all(ch in allowed for ch in my_seq)
if test==True:
    pass
else:
    st.write("Warning: Other characters beside A, C, T, or Gs were found.")
my_seq_reverse_compliment=my_seq.reverse_complement()
my_seq_compliment=my_seq.complement()


if input_type == 'DNA':
    convert_type = st.radio("What sequence type do you want as output?", ('Amino acid', 'RNA'))
    if convert_type == 'Amino acid':
        output_seq = my_seq.translate()
    else:
        output_seq = my_seq.transcribe()
else:
    convert_type = st.radio("What sequence type do you want as output?", ('Amino acid', 'DNA'))
    if convert_type == 'Amino acid':
        output_seq = my_seq.translate()
    else:
        output_seq = my_seq.back_transcribe()
st.markdown('Output:')
st.markdown(output_seq)