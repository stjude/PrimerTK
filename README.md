PrimerTK
========
[![Build Status](https://travis-ci.com/stjude/PrimerTK.svg?token=SGuFQqVLXJfs4J1ta2wA&branch=master)](https://travis-ci.com/stjude/PrimerTK)

[FULL USER MANUAL](https://stjude.github.io/PrimerTK/)

PrimerTK is a primer toolkit to assist in standard primer design, multiplex primer design, and primer design around complex structural variants.

Prerequisites
-------------

* python3
* primer3
* in silico PCR

Installation
------------

If you install via git clone, PrimerTK comes with a helper script to install in silico PCR and primer3.

This helper script will also be available via pip install, but will be located in your python specific site-packages location.

**NOTE: PrimerTK relies upon isPcr and primer3 which are external packages. Please verify that you have permissions to use those tools prior to installing them.**

We use the `installer.sh` script to install these packages. Once you have verified that you have proper permissions based on their respective licenses modify `installer.sh` as needed for your site.

This will put primer3 and isPcr in your `~/bin/`

Install all packages via git clone:

```bash
# Cloning the PrimerTK repository
$ git clone https://github.com/stjude/PrimerTK.git 
$ cd PrimerTK/scripts
$ ./installer.sh
$ cd ..
$ pip3 install . --user # can install to root as well.
```

Install all packages via pip:

```bash
# pip installing primer_tk and installing 3rd party programs using installer.sh
$ pip3 install primer_tk --user
$ cd ~/.local/lib/python3.x/site-packages/primer_tk/scripts
$ ./installer.sh
```

To ensure that PrimerTK was installed correctly, type `primer_tk -h`:

```
usage: primer_tk [-h] [-v]
                 {iterator,iterator_sv,pre,pre_sv,post,post_sv,tabix} ...

positional arguments:
  {iterator,iterator_sv,pre,pre_sv,post,post_sv,tabix}
                        Actions
    iterator            Iterator subparser
    iterator_sv         iterator_sv subparser
    pre                 Preprocessing for snv/indel
    pre_sv              Preprocessing for SV's
    post                Parses info from pcr_output
    post_sv             Parses info from pcr_output
    tabix               Tabix subparser

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

PrimerTK is intended to be run as a workflow broken up in the following way:

```
primer_tk iterator --> primer3 --> primer_tk pre --> isPcr --> primer_tk post --> primer_tk tabix
```
