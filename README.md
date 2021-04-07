# Bambu

BioAssays Model Builder, is a tool that aims to assist researchers in the development of new anti-tumor drugs, anti-Alzheimer's and other diseases. 

## Setup 
```
$ conda env create --file environment.yml
```
## Application
### Pre-Process

Bambu uses data obtained from PubChem in CSV and SDF format. In addition to some approaches such as SMOTE and Tomek Links for balancing the dataset.

The experimental data is stored in assays_csv and structural data is stored in assays_sdf.
The results of these data after the models have elapsed are stored in output_csv.

Balancing the dataset (balancing_strategy) uses four different approaches:

* Undersampling: Removes examples from the majority class to match the minority class and balance the dataset (figure 1a).

* Oversampling: adds more examples to the minority class to match the majority class by balancing the dataset (figure 1b).

![](https://miro.medium.com/max/725/1*7xf9e1EaoK5n05izIFBouA.png)

* Tomek Links: pairs of observations are selected and those of the majority class are excluded in order to balance the dataset (figure 2). 

![](https://raw.githubusercontent.com/rafjaa/machine_learning_fecib/master/src/static/img/tomek.png?v=2)

* SMOTE: * Synthetic Minority Oversampling Technique *, synthesizes elements of the minority class, where elements already exist. It is a simple way to generate synthetic samples at random from the attributes of instances of the minority class (figure 3). 

![](https://miro.medium.com/max/734/1*yRumRhn89acByodBz0H7oA.png)

#### Arguments

```
usage: bambu-preprocess [-h] --assays_csv ASSAYS_CSV --assays_sdf ASSAYS_SDF
                        --output_csv OUTPUT_CSV
                        [--balancing_strategy {random_undersampling,random_oversampling,smote}]
                        [--remove_tomek_links]

optional arguments:
  -h, --help            show this help message and exit
  --assays_csv ASSAYS_CSV
                        path to PubChem BioAssays CSV file
  --assays_sdf ASSAYS_SDF
                        path to PubChem BioAssays SDF file
  --output_csv OUTPUT_CSV
                        path to output CSV
  --balancing_strategy {random_undersampling,random_oversampling,smote}
                        data balancing strategy
  --remove_tomek_links  remove tomek links
  ```