setup:
	@conda env create --file environment.yml || conda env update --file environment.yml

dirs:
	@mkdir -p data/raw
	@mkdir -p data/processed
	
download_data: dirs
	@wget -O data/cancer-drugs-ml-raw.tar.gz https://storage.googleapis.com/cancer-drugs-ml/cancer-drugs-ml-raw.tar.gz && \
	 cd data/ && \
	 tar -xvf cancer-drugs-ml-raw.tar.gz && \
	 rm cancer-drugs-ml-raw.tar.gz
