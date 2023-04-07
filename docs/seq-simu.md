# Sequence data simulation

## Install `SLiM` on Fedora

```bash
sudo dnf install cmake qt5-qtbase-devel gsl-devel
sudo ln -s /us/bin/qmake-qt5 /usr/bin/qmake
git clone https://github.com/MesserLab/SLiM
mkdir build
cd build
cmake ../SLiM
make
cmake -D BUILD_SLIMGUI=ON ../SLiM
make SLiMgui
mv SLiMgui slim eidos ~/.local/bin
```

### to update

```bash
cd SLiM
git pull
cd ../build
make
make SLiMgui
mv SLiMgui slim eidos ~/.local/bin
```

## Workshop materials

- [The workshop website](http://benhaller.com/workshops/workshops.html)
- [The training materials](http://benhaller.com/workshops/SLiM_Workshop_Online.zip)

## To simulate a cattle population

SLiM can output FASTA, VCF, MS, CSV, and TSV files. Many tools can be used to convert these files to other formats.

for `.trees` files, use [`tskit`](https://tskit.readthedocs.io/en/latest/), `msprime`, or `pyslim`.

Many statistics, or plots can be done within SLiM.

### A very simple population

### Simulation with 1000 genomes results

- [reference](https://github.com/zxc307/GWAS_simulation_handbook)

Notes: seems not possible for cattle. As ancestral data are not available for cattle.

```bash
## directory
wdir=~/data/seq-simu/1000G # for example
mkdir -p $wdir
cd $wdir

## data
### human data below, will change to cattle data
#### vcf
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#### ancestral sequence
wget -r -np -nd ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
for chr in {1..22..1}; do
    samtools faidx hs37d5.fa "$chr" > hs37d5_chr$chr.fa
done

## tools

### Fedora provides
sudo dnf install samtools vcftools bcftools

### Plink
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
unzip plink_linux_x86_64_20230116.zip
mv plink prettify .local/bin
rm toy* LI*

### SLiM, see above

### QC

```

## Some concepts

- callback: 
  - a function that is called at a specific time.
  - In computer programming, a callback is a function that is passed as an argument to another function. When that function completes, it executes the callback function, "calling back" to signal its completion. Callbacks are often used in asynchronous programming, where the main program loop does not wait for individual functions to complete before continuing with program flow.
