# Identify-AA-N-AA-Stacks

## Table of Contents
-[About](#about)<br>
-[Requirements](#requirements)<br>
-[How to use](#how-to-use)<br>

## About
The 'identify-aa-n-aa-stacks.py' tool takes an mmcif file as input with preferred parameters and returns identified stacks with their measures based on the given parameters.<br>

## Requirements
Biopython is required to be installed while working with the 'identify-aa-n-aa-stacks.py' tool. To install biopython:<br>
For Windows:<br>
```bash
pip install biopython
```
For Linux/MacOS:
```bash
pip3 install biopython
```
It is important that the mmcif file is in the current working directory.<br>

## How to use
An example execution with '5mrc.cif' file. The file can be found [**here**](https://www.rcsb.org/structure/5MRC).<br>
For Windows:<br>
```bash
python identify-aa-n-aa-stacks.py -i 5mrc.cif -d 5.0 -a 40.0
```
For Linux/MacOS:<br>
```bash
python3 identify-aa-n-aa-stacks.py -i 5mrc.cif -d 5.0 -a 40.0
```
Here, the command line arguments:<br>
"-i"= name of the mmcif file.<br>
"-d"= maximum distance between the centroids of nucleotides and amino acid rings (in angstroms).<br>
"-a"= interplanar angle / maximum tilt angle between nucleotides and amino acid rings (in degrees).<br>
<br>
The execution will take some time.<br>
When completed, the output should look like the following:<br>
<br>
### Output:
The output format:<br>
stack no : i<br>
amino_acid_1 : chain.residue_name.residue_number<br>
nucleotide   : chain.residue_name.residue_number<br>
amino_acid_2 : chain.residue_name.residue_number<br>
<br>
calculated values :<br>
amino_acid_1 :<br>
d     : distance between the centroids of amino_acid_1 ring/rings and the nucleotide ring/rings (in angstroms)<br>
angle : angle between the planes of amino_acid_1 ring/rings and the nucleotide ring/rings (in degrees)<br>
n1    : normal distance between the centroid of amino_acid_1 ring/rings and the plane of nucleotide ring/rings (in angstroms)<br>
n2    : normal distance between the centroid of nucleotide ring/rings and the plane of amino_acid_1 ring/rings (in angstroms)<br>
<br>
amino_acid_2 :<br>
d     : distance between the centroids of amino_acid_2 ring/rings and the nucleotide ring/rings (in angstroms)<br>
angle : angle between the planes of amino_acid_2 ring/rings and the nucleotide ring/rings (in degrees)<br>
n1    : normal distance between the centroid of amino_acid_2 ring/rings and the plane of nucleotide ring/rings (in angstroms)<br>
n2    : normal distance between the centroid of nucleotide ring/rings and the plane of amino_acid_2 ring/rings (in angstroms)<br>
.......................................................................................................................................<br>
<br>
For your given parameter, the identified stacks are :<br>
<br>
stack no : 1<br>
amino_acid_1 : VV.PRO.64<br>
nucleotide   : aa.U.133<br>
amino_acid_2 : VV.HIS.60<br>
<br>
calculated values :<br>
amino_acid_1 :<br>
d     : 3.718252068530238<br>
angle : 154.3352508544922<br>
n1    : 3.713069200515747<br>
n2    : 3.5501163005828857<br>
<br>
amino_acid_2 :<br>
d     : 3.770456467372495<br>
angle : 148.42141723632812<br>
n1    : 3.6175615787506104<br>
n2    : 2.6567392349243164<br>
<br>
.............................................................................................................................................<br>
<br>
stack no : 2<br>
amino_acid_1 : UU.HIS.208<br>
nucleotide   : aa.U.1635<br>
amino_acid_2 : AA.TRP.85<br>
<br>
calculated values :<br>
amino_acid_1 :<br>
d     : 3.8716879195244296<br>
angle : 34.045936584472656<br>
n1    : 3.426558256149292<br>
n2    : 3.811905860900879<br>
<br>
amino_acid_2 :<br>
d     : 4.095622326113478<br>
angle : 5.349764823913574<br>
n1    : 3.704312562942505<br>
n2    : 3.6734800338745117<br>
<br>
.............................................................................................................................................<br>
