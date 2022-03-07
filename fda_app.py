# Filter FDA Approved Drugs
import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors

st.markdown('<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">', unsafe_allow_html=True)

st.markdown("""
<nav class="navbar fixed-top navbar-expand-lg navbar-dark" style="background-color: #3498DB;">
  <a class="navbar-brand" href="#">Project-III</a>
  <button class="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
    <span class="navbar-toggler-icon"></span>
  </button>
  <div class="collapse navbar-collapse" id="navbarNav">
    <ul class="navbar-nav">
      <li class="nav-item active">
        <a class="nav-link disabled" href="https://punitjain-jp.github.io/Dashboard/">Home <span class="sr-only">(current)</span></a>
      </li>
      <li class="nav-item">
        <a class="nav-link" href="https://share.streamlit.io/punitjain-jp/project/main/solubility-app.py" target="_blank">Solubility Application</a>
      </li>
      <li class="nav-item">
        <a class="nav-link" href="https://share.streamlit.io/punitjain-jp/dna/main/dna-app.py" target="_blank">Nucleotide Composition</a>
      </li>
      <li class="nav-item">
        <a class="nav-link" href="https://share.streamlit.io/punitjain-jp/moldesc/main/app.py" target="_blank">Molecular Descriptors</a>
      </li>
      <li class="nav-item">
        <a class="nav-link" href="https://share.streamlit.io/rahul97532/bioactivity-prediction/app.py" target="_blank">Bioactivity Application</a>
      </li>
    </ul>
  </div>
</nav>
""", unsafe_allow_html=True)


st.title("**Filter FDA Approved Drugs**")
st.markdown("""The process of filtration uses the *Lipinski's Rule-of-Five*""")

@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(
        "https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
    ).dropna()
    return df

# Calculate descriptors, Given a smiles string (ex. C1CCCCC1).
def calc_mw(smiles_string):
    """calculate and return the molecular weight"""
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)

def calc_logp(smiles_string):
    """calculate and return the LogP"""
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)

def calc_NumHDonors(smiles_string):
    """calculate and return the NumHDonors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)

def calc_NumHAcceptors(smiles_string):
    """calculate and return the NumHAcceptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHAcceptors(mol)


# Copy the dataset so any changes are not applied to the original cached version
df = download_dataset().copy()
df["MW"] = df.apply(lambda x: calc_mw(x["smiles"]), axis=1)
df["LogP"] = df.apply(lambda x: calc_logp(x["smiles"]), axis=1)
df["NumHDonors"] = df.apply(lambda x: calc_NumHDonors(x["smiles"]), axis=1)
df["NumHAcceptors"] = df.apply(lambda x: calc_NumHAcceptors(x["smiles"]), axis=1)


# Sidebar panel
st.sidebar.header('Set parameters')
st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')
weight_cutoff = st.sidebar.slider(
    label="Molecular weight",
    min_value=0,
    max_value=1000,
    value=500,
    step=10,
)
logp_cutoff = st.sidebar.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=5,
    step=1,
)
NumHDonors_cutoff = st.sidebar.slider(
    label="NumHDonors",
    min_value=0,
    max_value=15,
    value=5,
    step=1,
)
NumHAcceptors_cutoff = st.sidebar.slider(
    label="NumHAcceptors",
    min_value=0,
    max_value=20,
    value=10,
    step=1,
)

df_result = df[df["MW"] < weight_cutoff]
df_result2 = df_result[df_result["LogP"] < logp_cutoff]
df_result3 = df_result2[df_result2["NumHDonors"] < NumHDonors_cutoff]
df_result4 = df_result3[df_result3["NumHAcceptors"] < NumHAcceptors_cutoff]

st.write(df_result4.shape)
st.write(df_result4)


raw_html = mols2grid.display(df_result4,
                            #subset=["Name", "img"],
                            subset=["img", "Name", "MW", "LogP", "NumHDonors", "NumHAcceptors","cns_drug"],
                            mapping={"smiles": "SMILES", "generic_name": "Name"})._repr_html_()
components.html(raw_html, width=900, height=1100, scrolling=False)

st.markdown("""
<script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
""", unsafe_allow_html=True)
