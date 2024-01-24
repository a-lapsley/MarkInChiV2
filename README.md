# MarkInChiV2
A tool to generate MarkInChi codes from v3000 .mol files. It can also parse a MarkInChi code to produce a molblock and generate a list of InChIs from a MarkInChi

## Dependencies
 - Python 3.11.5+ (older versions may still work but have not been tested)
 - RDKit 2023.09.04+

## Currently supported features
- R groups with single connections
- Isotopically labelled compounds are supported
- MarkInChI string does not depend on labelling of atoms in the Molfile which should resolve issues with stereochemistry and canonicity in the old MarkInChI software
- Variable attachments
- Lists of atoms
- Any nested combinations of these Markush structures

## Limitations
- Compounds with more than 200 R groups
- Xenon, Radon and Neon containing compounds

## Usage
To use this software, run the program `MarkinchiGUI.py` to launch a GUI program that contains various functionalities. The Jupyter notebook `Demos.ipynb` also contains demonstrations of the functionalities of this software. 

## Generating .mol files with Markush features
All molfiles used for testing and developing this software were produced using MarvinSketch 23.17. It may be possible to generate molfiles using other methods, however these have not yet been tested so may cause issues.

## Known issues
 - Variable attachments which are defined to have endpoints which have no hydrogen atoms on the core structure are known to cause errors. This includes molecules where the endpoints could be tautomeric, and, for example, an endpoint defined to be a quaternary carbon.
 - This program does not currently handle stereochemistry correctly when the core of the molecule is achiral or symmetric, but the presence of Markush features makes it chiral / asymmetric

