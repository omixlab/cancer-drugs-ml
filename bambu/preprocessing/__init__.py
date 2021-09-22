from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import RandomOverSampler
from imblearn.under_sampling import TomekLinks
from imblearn.over_sampling import SMOTE
from BorutaShap import BorutaShap 
import pandas as pd
import rdkit
import argparse
import json 

def preprocess_dataset(assays_csv, assays_sdf, output_csv_train, output_csv_test, balancing_strategy, remove_tomek_links, feature_selection):
    df_labels = pd.read_csv(assays_csv)
    df_labels = df_labels.drop([0, 1, 2])
    df_labels = df_labels[['PUBCHEM_SID', 'PUBCHEM_ACTIVITY_OUTCOME']]
    df_labels['ACTIVITY'] = (df_labels['PUBCHEM_ACTIVITY_OUTCOME'] =='Active').astype(float)
    df_features = compute_descriptors(assays_sdf)
    df_joined = pd.merge(df_features, df_labels, on = ['PUBCHEM_SID'] )
    df_joined = df_joined.drop(['PUBCHEM_SID', 'PUBCHEM_ACTIVITY_OUTCOME'], axis=1)
    df_joined = df_joined.drop_duplicates(subset=[column for column in df_joined.columns if column != 'ACTIVITY'])
    df_joined = df_joined.dropna()

    X = df_joined.drop(['ACTIVITY'], axis=1)
    y = df_joined['ACTIVITY']

    X_train, X_test, y_train, y_test = train_test_split(X, y)

    if balancing_strategy:
        df_joined_train = balance_dataset(X_train, y_train, balancing_strategy, remove_tomek_links)
        X_train = df_joined_train.drop(['ACTIVITY'], axis=1)
        y_train = df_joined_train['ACTIVITY']
    
    df_joined_test = balance_dataset(X_test, y_test, 'random_undersampling', False)
    X_test = df_joined_test.drop(['ACTIVITY'], axis=1)
    y_test = df_joined_test['ACTIVITY']

    if feature_selection:
        Feature_Selector = create_feature_selector(X_train, y_train)
        X_train = X_train[Feature_Selector.columns]
        X_test = X_test[Feature_Selector.columns]
        with open(feature_selection, 'w') as handle:
            handle.write(json.dumps(list(Feature_Selector.columns)))

    df_train = X_train
    df_train['ACTIVITY'] = y_train
    df_train.to_csv(output_csv_train, index=False)
    df_test = X_test
    df_test['ACTIVITY'] = y_test
    df_test.to_csv(output_csv_test, index=False)
    
def create_feature_selector(X, y):
    Feature_Selector = BorutaShap(importance_measure='shap', classification=True)
    Feature_Selector.fit(X=X, y=y, n_trials=100, random_state=0)
    return Feature_Selector

def compute_descriptors(assays_sdf):
    mols = Chem.SDMolSupplier(assays_sdf)
    df_features = pd.DataFrame()
    success_mols = 0
    failed_mols  = 0
    for m, mol_features in enumerate(map(compute_mol_descriptors_wrapper, mols)):
        if mol_features:
            success_mols += 1
            print('success on molecule: ', m, f'success: {success_mols}; failed: {failed_mols}')
            df_features = df_features.append(mol_features, ignore_index=True)
        else: 
            print('failed on molecule: ', m, f'success: {success_mols}; failed: {failed_mols}')
            failed_mols += 1
            pass
    return df_features

   
def compute_mol_descriptors(mol):
    mol_features = {}
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
    return mol_features

def compute_mol_descriptors_wrapper(mol):
    try:
        mol_features = compute_mol_descriptors(mol)
        mol_features['PUBCHEM_SID']=mol.GetPropsAsDict()['PUBCHEM_SUBSTANCE_ID']
        return mol_features
    except:
        return None

def balance_dataset(X, y, balacing_strategy, remove_tomek_links):
    if balacing_strategy == 'random_undersampling':
        sampler = RandomUnderSampler()
    elif balacing_strategy == 'random_oversampling':
        sampler = RandomOverSampler()
    elif balacing_strategy == 'smote':
        sampler = SMOTE()
    df_sampler,labels = sampler.fit_resample(X, y)
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
    argument_parser.add_argument('--output_csv_train', required=True, help='path to output CSV train')
    argument_parser.add_argument('--output_csv_test', required=True, help='path to output CSV test')
    argument_parser.add_argument('--balancing_strategy', default=None, choices=['random_undersampling', 'random_oversampling', 'smote'], help='data balancing strategy')
    argument_parser.add_argument('--remove_tomek_links', default=False, action='store_true', help='remove tomek links')
    argument_parser.add_argument('--feature_selection', default=None, help='path to a file containing the columns selected by feature selection algorithm')
    arguments = argument_parser.parse_args()
    preprocess_dataset(assays_csv=arguments.assays_csv, assays_sdf=arguments.assays_sdf, output_csv_test=arguments.output_csv_test, output_csv_train=arguments.output_csv_train, balancing_strategy=arguments.balancing_strategy, remove_tomek_links=arguments.remove_tomek_links, feature_selection=arguments.feature_selection)
    

if __name__=="__main__":
    main()