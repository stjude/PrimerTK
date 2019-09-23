PrimerTK
========

[![Build Status](https://travis-ci.com/stjude/PrimerTK.svg?token=SGuFQqVLXJfs4J1ta2wA&branch=development)](https://travis-ci.com/stjude/PrimerTK)

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

To install isPcr and primer3 using the helper script, edit the file and uncomment the lines in between the text boxes and then run:

`installer.sh`

This will put primer3 and isPcr in your `~/bin/`

Install all packages via git clone:

```bash
# Cloning the PrimerTK repository
$ git clone https://github.com/stjude/PrimerTK.git 
$ cd PrimerTK/scripts
# edit installer.sh removing comments to install primer3 and isPcr
$ ./installer.sh
$ cd ..
$ pip3 install . --user # can install to root as well.
```

Install all packages via pip:

```bash
# pip installing primer_tk and installing 3rd party programs using installer.sh
$ pip3 install primer_tk --user
$ cd ~/.local/lib/python3.x/site-packages/primer_tk/scripts
# edit installer.sh to remove comment in front of lines between text boxes.
$ ./installer.sh
```


