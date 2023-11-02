# MarkInChiV2
A tool to generate MarkInChi codes from v3000 .mol files. 

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
To generate a MarkInChI from a v3000 .mol file, run the script `MarkinchiGenerator.py` from the command line using the argument `-i` to specify the input file. For now, files must be in the same directory as the script, but this will be improved in future updates.
For example, `python MarkinchiGenerator.py -i molfiles\test1.mol`.
Use the argument `-d` to run the script in debug mode, which shows the molecule for each fragment as it's MarkInChI gets generated. This is useful for seeing what's going on in the script but requires manually closing each popup window for the script to proceed. 