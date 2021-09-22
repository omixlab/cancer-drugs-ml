for cancer_directory in data/raw/melanoma; do 

    cancer_name=$(basename $cancer_directory)

    for assay_name in 794; do
        assay_directory=$cancer_directory/$assay_name
        echo $cancer_name $assay_name

        mkdir -p data/models/$cancer_name/$assay_name/no_balancing_strategy
        mkdir -p data/models/$cancer_name/$assay_name/random_undersampling
        mkdir -p data/models/$cancer_name/$assay_name/random_oversampling
        mkdir -p data/models/$cancer_name/$assay_name/smote
        
        bambu-train-model --preprocess_train_csv data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/train.csv --preprocess_test_csv data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/test.csv --output_model data/models/$cancer_name/$assay_name/no_balancing_strategy/model.pickle 
        bambu-train-model --preprocess_train_csv data/preprocessing/$cancer_name/$assay_name/random_undersampling/train.csv --preprocess_test_csv data/preprocessing/$cancer_name/$assay_name/random_undersampling/test.csv --output_model data/models/$cancer_name/$assay_name/random_undersampling/model.pickle
        bambu-train-model --preprocess_train_csv data/preprocessing/$cancer_name/$assay_name/random_oversampling/train.csv --preprocess_test_csv data/preprocessing/$cancer_name/$assay_name/random_oversampling/test.csv --output_model data/models/$cancer_name/$assay_name/random_oversampling/model.pickle
        bambu-train-model --preprocess_train_csv data/preprocessing/$cancer_name/$assay_name/smote/train.csv --preprocess_test_csv data/preprocessing/$cancer_name/$assay_name/smote/test.csv --output_model data/models/$cancer_name/$assay_name/smote/model.pickle  

    done     

done