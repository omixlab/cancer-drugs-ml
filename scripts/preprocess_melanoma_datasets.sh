for cancer_directory in data/raw/melanoma; do 

    cancer_name=$(basename $cancer_directory)

    for assay_name in 1321; do
        assay_directory=$cancer_directory/$assay_name
        echo $cancer_name $assay_name
        
        mkdir -p data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy
        mkdir -p data/preprocessing/$cancer_name/$assay_name/random_undersampling
        mkdir -p data/preprocessing/$cancer_name/$assay_name/random_oversampling
        mkdir -p data/preprocessing/$cancer_name/$assay_name/smote
        
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/test.csv #--feature_selection data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/features
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/random_undersampling/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/random_undersampling/test.csv #--feature_selection data/preprocessing/$cancer_name/$assay_name/random_undersampling/features
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/random_oversampling/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/random_oversampling/test.csv #--feature_selection data/preprocessing/$cancer_name/$assay_name/random_oversampling/features
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/smote/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/smote/test.csv #--feature_selection data/preprocessing/$cancer_name/$assay_name/smote/features
        
    done

done