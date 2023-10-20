# MarkInChiV2
A tool to generate MarkInChi codes from v3000 .mol files. 

## Currently supported features
- R groups with single connections
- R groups with multiple connections where all connection points are identical, e.g. different length alkyl chains. 
- R groups with multiple connections where connection points are different will generate a somewhat meaningful MarkInChI string however the connectivity is not well defined 
- Isotopically labelled compounds are supported
- MarkInChI string does not depend on initial R group labelling which should resolve issues with stereochemistry in the old MarkInChI software
- Variable attachments where the attachment is either a single atom or structure (R groups not yet supported). Note that these attachments are not yet canonicalised in terms of identical connection points and so the generated InChI does depend on how the MolFile is defined

## Limitations
- Compounds with more than 200 R groups
- Xenon containing compounds
- Compounds containing other Markush functionality such as variable position - support for these features will be added in future updates. 

## Usage
To generate a MarkInChI from a v3000 .mol file, run the script `markmol3000.py` from the command line using the argument `-i` to specify the input file. For now, files must be in the same directory as the script, but this will be improved in future updates.
For example, `python markmol3000.py -i molfiles\test1.mol`.
Use the argument `-d` to run the script in debug mode, which shows the molecule for each fragment as it's MarkInChI gets generated. This is useful for seeing what's going on in the script but requires manually closing each popup window for the script to proceed. 