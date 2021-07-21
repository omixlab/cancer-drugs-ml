for cancer_directory in data/raw/*;do 
    cancer_name=$(basename $cancer_directory)
    for assay_directory in $cancer_directory/*;do
        assay_name=$(basename $assay_directory)
        mkdir -p data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy
        mkdir -p data/preprocessing/$cancer_name/$assay_name/random_undersampling
        mkdir -p data/preprocessing/$cancer_name/$assay_name/random_oversampling
        mkdir -p data/preprocessing/$cancer_name/$assay_name/smote
        mkdir -p data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy_tomek_links
        mkdir -p data/preprocessing/$cancer_name/$assay_name/random_undersampling_tomek_links
        mkdir -p data/preprocessing/$cancer_name/$assay_name/random_oversampling_tomek_links
        mkdir -p data/preprocessing/$cancer_name/$assay_name/smote_tomek_links
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy/features
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/random_undersampling/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/random_undersampling/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/random_undersampling/features
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/random_oversampling/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/random_oversampling/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/random_oversampling/features
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/smote/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/smote/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/smote/features
        
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy_tomek_links/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy_tomek_links/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/no_balancing_strategy_tomek_links/features --remove_tomek_links
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/random_undersampling_tomek_links/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/random_undersampling_tomek_links/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/random_undersampling_tomek_links/features --remove_tomek_links
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/random_oversampling_tomek_links/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/random_oversampling_tomek_links/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/random_oversampling_tomek_links/features --remove_tomek_links
        bambu-preprocess --assays_csv data/raw/$cancer_name/$assay_name/$assay_name.csv --assays_sdf data/raw/$cancer_name/$assay_name/$assay_name.sdf --output_csv_train data/preprocessing/$cancer_name/$assay_name/smote_tomek_links/train.csv --output_csv_test data/preprocessing/$cancer_name/$assay_name/smote_tomek_links/test.csv --feature_selection data/preprocessing/$cancer_name/$assay_name/smote_tomek_links/features --remove_tomek_links
        break
    done
    break 
done