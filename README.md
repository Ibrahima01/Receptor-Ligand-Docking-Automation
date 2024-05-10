# Receptor-Ligand Docking Automation

## Description
Receptor-ligand docking is vital in drug discovery for identifying potential drug targets. Traditionally, this process involves using multiple software packages for tasks such as energy minimization, receptor and ligand preparation, box setup, docking, and visualization. Moreover, this process needs to be repeated for each receptor, making it time-consuming and tedious. To address this challenge, I have developed a program to automate the entire docking process with just one click.

This notebooks are Python 3 compatible.

## Requirements
This project is reliant on a variety of Python librairies such as :

- [Autodock Vina](https://autodock-vina.readthedocs.io/en/latest/) 
- [OpenBabel](http://openbabel.org/wiki/Main_Page)
- [RDKit](https://www.rdkit.org/)
- [PDBFixer](https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html)
- [MDAnalysis](https://www.mdanalysis.org/)
- [OpenM](http://docs.openmm.org/latest/api-python/)
- [PDBFixer](https://github.com/openmm/pdbfixer)

## Installation:
**1. Installing all dependencies one by one:**

1.1. Create a conda enviroment  

```
conda create -n Docking python=3.8
conda activate Docking
```

1.2. Install de dependencies 

- AutoDock Vina
```
pip install vina
```

- OpenBabel (Pybel)
```
conda install -c conda-forge openbabel
```

- PDBFixer
```
conda install -c conda-forge pdbfixer
```

- RDKit 

```
conda install rdkit cython
```

- MDAnalysis
```
conda install MDAnalysisTests
```

- Openmm
```
conda install -c conda-forge openmm
```

## How to use:
Search for software using the next command :

```git clone https://github.com/Ibrahima01/Receptor-Ligand-Docking-Automation.git```.

Then delete all proteins in the ```proteins``` folder, and delete ligands in the ``ligand'' folder.
You'll be using your own proteins and ligands.


## Limitation 
>* This program currently supports a limited number of receptor and ligand file formats. Additional format support may be added in future updates.
>* The docking algorithm may not accurately predict binding affinities in all cases and should be validated with experimental data.