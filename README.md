# mg_metawrap



## Install metaWRAP

``` shell
sh Miniconda3-py39_4.11.0-Linux-x86_64.sh 
source /root/.bashrc 
conda install -y mamba
conda install -c conda-forge -y mamba
mamba create -y -n metawrap-env python=2.7
conda activate metawrap-env
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky
mamba install --only-deps -c ursky metawrap-mg
ktUpdateTaxonomy.sh 

# git download metaWRAP
cp -r metaWRAP/ /root/
echo 'export PATH="/root/metaWRAP/bin/:$PATH"' >>/root/.bashrc 

# download CHECKM database
cp -r CHECKM_FOLDER/ /root/
checkm data setRoot /root/CHECKM_FOLDER

apt-get update
apt-get install nodejs npm
npm install less

## update blast
conda install --use-local blast-2.12.0-hf3cf87c_4.tar.bz2 

git clone https://github.com/tseemann/prokka.git $HOME/prokka
/usr/bin/perl $HOME/prokka/bin/prokka --setupdb
sed -i.bak -r 's/cmd="prokka/cmd="\/root\/miniconda3\/envs\/metawrap\-env\/bin\/perl \/root\/prokka\/bin\/prokka/g' /opt/conda/bin/metawrap-modules/annotate_bins.sh

```