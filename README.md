# MarkInChiV2
A tool to generate MarkInChi codes from v3000 .mol files. 

## Currently supported features
- R groups with single connections
- R groups with multiple connections where all connection points are identical, e.g. different length alkyl chains. 
- R groups with multiple connections where connection points are different will generate a somewhat meaningful MarkInChI string however the connectivity is not well defined 
- Isotopically labelled compounds are supported
- MarkInChI string does not depend on initial R group labelling which should resolve issues with stereochemistry in the old MarkInChI software

## Limitations
- Compounds with more than 201 R groups
- Xenon containing compounds
- Compounds containing other Markush functionality such as variable position - support for these features will be added in future updates. 

## Usage
To generate a MarkInChI from a v3000 .mol file, update the string `FILENAME` in `markmol3000.py` to the path to the molfile within this directory, e.g. `FILENAME = "molfiles\\test0.mol"`. A more user-friendly interface will be added in future - this is a temporary solution for quickly debugging and testing the program. 