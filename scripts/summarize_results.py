import pandas as pd
import json

cancer_assays = {
    'melanoma': ['39']
}

balancing_strategies = ['no_balancing_strategy', 'random_undersampling', 'random_oversampling', 'smote']

algorithms = [
    'DecisionTreeClassifier',
    'GradientBoostingClassifier',
    'RandomForestClassifier',
    'ExtraTreesClassifier',
    'SVC']

subsets = ['train', 'test']

for cancer, assays in cancer_assays.items():
    for assay in assays:
        df = pd.DataFrame(columns=['algorithm', 'balancing method', 'precision', 'recall', 'f1-score', 'accuracy'])
        for balancing_strategy in balancing_strategies:
            for algorithm in algorithms:
                for subset in subsets:
                    with open(f'data/models/{cancer}/{assay}/{balancing_strategy}/model.pickle.{algorithm}.{subset}.json') as handle:
                        classification_report = json.loads(handle.read())
                        precision  = classification_report['1.0']['precision']
                        recall     = classification_report['1.0']['recall']
                        f1_score   = classification_report['1.0']['f1-score']
                        weight_acc = classification_report['accuracy']
                        df = df.append(
                            {
                                'algorithm': algorithm, 
                                'balancing method': balancing_strategy, 
                                'precision': precision, 
                                'recall': recall, 
                                'f1-score': f1_score, 
                                'accuracy': weight_acc
                            }, ignore_index=True
                        )
        df.to_csv(f'report.{cancer}.{assay}.csv', index=False)

