import streamlit as st

st.set_page_config(
    page_title="Bio Tools",
    page_icon=":dna:"
)
st.title("Bio tools")
st.markdown('This website contains useful Bio-Tools, the tools focus on different subjects. How the apps work is documented on this [website](https://www.dna-vorm.com/bio-tools). Open the sidebar to select a tool.')
st.markdown('-🧬 Tools regarding DNA, RNA, or Amino Acid sequences')
st.markdown('-🧪 Tools regarding fermentation')
st.sidebar.success("\"The best website for your Bio Tools\" \n Watson and Crick")