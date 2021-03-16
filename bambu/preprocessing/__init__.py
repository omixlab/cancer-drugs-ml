from rdkit import Chem
from rdkit.Chem import Descriptors
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import TomekLinks
from imblearn.over_sampling import SMOTE
import pandas as pd
import rdkit
import argparse

def preprocess_dataset(assays_csv, assays_sdf, output_csv, balancing_strategy, remove_tomek_links):
    df_labels = pd.read_csv(assays_csv)
    df_labels = df_labels.drop([0, 1, 2])
    df_labels = df_labels[['PUBCHEM_SID', 'PUBCHEM_ACTIVITY_OUTCOME']]
    df_labels['ACTIVITY'] = (df_labels['PUBCHEM_ACTIVITY_OUTCOME'] =='Active').astype(float)
    df_features = compute_descriptors(assays_sdf)
    df_joined = pd.merge(df_features, df_labels, on = ['PUBCHEM_SID'] )
    if balancing_strategy:
        df_balanced = balance_dataset(df_joined, balancing_strategy, remove_tomek_links)
        df_balanced.to_csv(output_csv, index=False)
    else: 
        df_joined.to_csv(output_csv, index=False)

def compute_descriptors(assays_sdf):
    mols = Chem.SDMolSupplier(assays_sdf)
    df_features = pd.DataFrame()
    for mol in mols:
        if mol == None:
            break
        mol_features={}
        mol_features['PUBCHEM_SID']=mol.GetPropsAsDict()['PUBCHEM_SUBSTANCE_ID']
        mol_features['min_abs_partial_charge'] = Descriptors.MinAbsPartialCharge(mol)
        mol_features['tpsa'] = Descriptors.TPSA(mol) 
        mol_features['exact_mol_wt'] = Descriptors.ExactMolWt(mol)
        mol_features['max_abs_partial_charge'] = Descriptors.MaxAbsPartialCharge(mol) 
        mol_features['num_radical_eletrons'] = Descriptors.NumRadicalElectrons(mol)
        mol_features['mol_log_p'] = Descriptors.MolLogP(mol)
        mol_features['mol_mr'] = Descriptors.MolMR(mol)
        mol_features['mol_wt'] = Descriptors.MolWt(mol)
        mol_features['heavy_atom_count'] = Descriptors.HeavyAtomCount(mol)
        mol_features['heavy_atom_mol_wt'] = Descriptors.HeavyAtomMolWt(mol)
        mol_features['nhoh_count'] = Descriptors.NHOHCount(mol)
        mol_features['no_count'] = Descriptors.NOCount(mol)
        mol_features['num_h_acceptors'] = Descriptors.NumHAcceptors(mol)
        mol_features['num_h_donors'] = Descriptors.NumHDonors(mol)
        mol_features['num_hetero_atoms'] = Descriptors.NumHeteroatoms(mol)
        mol_features['num_rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
        mol_features['num_valence_electrons'] = Descriptors.NumValenceElectrons(mol)
        mol_features['ring_count'] = Descriptors.RingCount(mol)
        mol_features['fp_density_morgan1'] = rdkit.Chem.Descriptors.FpDensityMorgan1(mol)
        mol_features['fp_density_morgan2'] = rdkit.Chem.Descriptors.FpDensityMorgan2(mol)
        mol_features['fp_density_morgan3'] = rdkit.Chem.Descriptors.FpDensityMorgan3(mol)
        mol_features['balabanj'] = rdkit.Chem.GraphDescriptors.BalabanJ(mol)
        mol_features['bertzct'] = rdkit.Chem.GraphDescriptors.BertzCT(mol)
        mol_features['ipc'] = rdkit.Chem.GraphDescriptors.Ipc(mol)
        mol_features['chi0'] = rdkit.Chem.GraphDescriptors.Chi0(mol)
        mol_features['chi1'] = rdkit.Chem.GraphDescriptors.Chi1(mol)
        mol_features['kappa1'] = rdkit.Chem.GraphDescriptors.Kappa1(mol)
        mol_features['hallkier_alpha'] = rdkit.Chem.GraphDescriptors.HallKierAlpha(mol)
        df_features = df_features.append(mol_features, ignore_index=True)
    return df_features

def balance_dataset(df, balacing_strategy, remove_tomek_links):
    if balacing_strategy == 'random_undersampling':
        sampler = RandomUnderSampler()
    elif balacing_strategy == 'random_oversampling':
        sampler = RandomOverSampler()
    elif balacing_strategy == 'smote':
        sampler = SMOTE()
    df_sampler,labels = sampler.fit_resample(df.drop(['ACTIVITY'],axis=1), df['ACTIVITY'])
    df_sampler['ACTIVITY'] = labels
    if remove_tomek_links:
        tomek_links_sampler = TomekLinks()
        df_sampler,labels = tomek_links_sampler.fit_resample(df_sampler.drop(['ACTIVITY'], axis=1), df_sampler['ACTIVITY'])
        df_sampler['ACTIVITY'] = labels 
    return df_sampler

def main():
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument('--assays_csv', required=True, help='path to PubChem BioAssays CSV file')
    argument_parser.add_argument('--assays_sdf', required=True, help='path to PubChem BioAssays SDF file') 
    argument_parser.add_argument('--output_csv', required=True, help='path to output CSV')
    argument_parser.add_argument('--balancing_strategy', default=None, choices=['random_undersampling', 'random_oversampling', 'smote'], help='data balancing strategy')
    argument_parser.add_argument('--remove_tomek_links', default=False, action='store_true', help='remove tomek links')
    arguments = argument_parser.parse_args()
    preprocess_dataset(assays_csv=arguments.assays_csv, assays_sdf=arguments.assays_sdf, output_csv=arguments.output_csv, balancing_strategy=arguments.balancing_strategy, remove_tomek_links=arguments.remove_tomek_links)
    

if __name__=="__main__":
    main()