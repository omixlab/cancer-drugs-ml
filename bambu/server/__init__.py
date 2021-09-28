from flask import Flask, render_template, request
from rdkit import Chem
from bambu.preprocessing import compute_mol_descriptors
import argparse 
import json
import pandas as pd
import pickle 

app=Flask(__name__)
MODELS=[]
def run_server(models_json):
    json_content=open(models_json).read()
    models=json.loads(json_content)
    for model in models:
        model["model"]=pickle.load(open(model["path"], "rb"))
        MODELS.append(model)
    app.run(host="0.0.0.0")

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/predict", methods=["POST"])
def predict():
    molecule_smiles=request.form.get("molecule")
    molecule=Chem.MolFromSmiles(molecule_smiles)
    descriptors=compute_mol_descriptors(molecule)
    df_descriptors = pd.DataFrame([descriptors])
    results = []
    for model in MODELS: 
        prediction = model["model"].predict(df_descriptors)
        results.append((model["name"], prediction))
    return render_template("results.html", results=results)

def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('--models_json', required=True, help='path to json file with training models')
    arguments = argument_parser.parse_args()
    run_server(models_json=arguments.models_json)

if __name__=="__main__":
    main()