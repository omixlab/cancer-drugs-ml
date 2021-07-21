for cancer_directory in data/raw/*;do 
    cancer_name=$(basename $cancer_directory)
    for assay_directory in $cancer_directory/*;do
        assay_name=$(basename $assay_directory)
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/model.pickle 
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/random_undersampling/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/random_undersampling/model.pickle
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/random_oversampling/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/random_oversampling/model.pickle
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/smote/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/smote/model.pickle
        
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy_tomek_links/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy_tomek_links/model.pickle
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/random_undersampling_tomek_links/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/random_undersampling_tomek_links/model.pickle 
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/random_oversampling_tomek_links/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/random_oversampling_tomek_links/model.pickle
        bambu-train-model --preprocess_csv data/preprocessing/$cancer_name/$assay_name/smote_tomek_links/train.csv --output_model data/preprocessing/$cancer_name/$assay_name/smote_tomek_links/model.pickle
        break
    done
    break 
done