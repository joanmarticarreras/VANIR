# VANIR

NextFlow pipeline for Virus Variant calling and de novo Assembly of Nanopore and Illumina Reads

-----------------------------------------------------------------------------------------------

## Introduction
The aim of this pipeline is to accuratley base-call and QC read information and _de novo_ assemble long-read data from Oxford Nanopore Technologies (ONT) reads to yield a complete dsDNA viral genome. The target genome is later corrected using multiple algorithms (Racon and Medaka) and a final polishing using a higher accuracy read data set from Illumina.

Later single nucleotide variants (SNV) and structural variants (SV) are called from Illumina and Nanopore reads, respectivley using a user-provided reference and to assess differences against standard and using the _de novo_ assembled genome in order to assess the diversity within the sample.

The pipeline is writen in NextFlow pipeline manager. This framework allow users to easly launch and coordinate the pipeline in different environments. It can easly resume failed runs and provide basic metrics of the consumed resources (see options section).

Additionally, all software has been bundled in a Singularity container, containing its own operating system (Ubuntu 16.04) and all libraries and programs needed to run the pipeline. Singularity makes this pipeline highly transportable, as can be run it creates its own environemnt; convenient, as everything has been pre-bundled for direct user usage; and reproducible, as libraries and programs are independent from general installation and will run always on the same manner.

This pipeline is partially GPU-enabled, concerning maninly ONT basecalling and ONT read correction with Medaka.

## Installation

Only the Singulairy image needs to be installed. For this first install Singularity as follows:

### Linux
Open the terminal and type the following:

```
VERSION=2.5.2
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
tar xvf singularity-$VERSION.tar.gz
cd singularity-$VERSION
./configure --prefix=/usr/local
make
sudo make install
```

### MAC
Open the terminal and type the following:

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew cask install virtualbox
brew cask install vagrant
brew cask install vagrant-manager

```

Download Singularity from [here](https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4)

```
mkdir singularity-vm
cd singularity-vm
vagrant destroy
vagrant init singularityware/singularity-2.4
vagrant up
vagrant ssh
```

### Windows:

First install [Git](https://git-for-windows.github.io/), [VirtualBox](https://www.virtualbox.org/wiki/Downloads), [Vagrant](https://www.vagrantup.com/downloads.html) and [Vagrant Manager](https://www.vagrantup.com/downloads.html).

Then install [Singularity](https://app.vagrantup.com/singularityware/boxes/singularity-2.4/versions/2.4), and run on GitBash:

```
mkdir singularity-vm
cd singularity-vm
vagrant destroy
vagrant init singularityware/singularity-2.4
vagrant up
vagrant ssh

```

After installing Singularity, VANIR image must be build from recepie or downloaded.

To build from recepie please follow:

```
git clone XXXX
cd XX
sudo singularity build --writable VANIR singularity-recepie/VANIR-singularity-recepie.txt
sudo cp main.nf VANIR/nextflow/
```

## Options
There are 2 sets of options, one comming from VANIR pipeline itself and one from the NextFlow manager:

```                                                                                                                                                                                                                                                                                                                                                                                                      Mandatory options:                                                                                                                                                                                             Folder to raw fast5 reads [PATH]:                                          --fast5                                                                                                                     Illumina 1st pair-end reads (can be .gz) [FILE]:                           --read1                                                                                                                     Illumina 2nd pair-end reads (can be .gz) [FILE]:                           --read2                                                                                                                     Illumina adapters [FILE]:                                                  --illumina_adapters                                                                                                         Sample prefix [FILE]:                                                      --prefix                                                                                                                    Reference genome [FILE]:                                                   --reference                                                                                                                 Sequencing contaminants [FILE]:                                            --contaminants                                                                                                              Number of CPUs asked [NUMERIC]:                                            --cpu                                                                                                                       GPU-enabled ('true' enabled // 'false' disabled) [Boolean]:                --gpu                                                                                                                       Guppy-basecalling mode ('precise' for hac // 'fast' for fast) [Boolean]:   --guppy_algorithm                                                                                                           Expected genome size (STRING, eg 35k, 8g or 250m):                         --genome_size                                                                                                                                                                                                                                                                                                                  
```
```
NextFlow options [OPTIONAL]:                                                                                                                                                                               Produce an html report with useful metrics of the pipeline [FILE]          -with-report                                                                                                                Produce a tabular file with tracings of each processes [FILE]              -with-trace                                                                                                                 Produce an html graphic of all process executed [FILE]                     -with-timeline                                                                                                              Produce a graph-image (.dot/.html/.pdf/.png/.svg) of the pipeline [FILE]   -with-dag       
```
## Run

To run this pipeline, run on terminal:

```                                                                                                                                                                                                                                                                                                                                                               
sudo singularity shell -B <BIND-PATH-DATA> --nv --writable VANIR
nextflow run nextflow/VANIR/ OPTIONS 
```
