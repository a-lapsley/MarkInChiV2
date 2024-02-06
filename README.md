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
- To run this software from the command line, call the program `Markinchi.py`
  - Available commands:
    - `help` - Displays a list of available commands
    - `generate <filename>` Generates a MarkInChI from a Mol file
    - `markinchitomolblock <markinchi>` Generates a V3000 Mol Block from a MarkInChI
    - `molfiletoinchilist <filename>` Generates a list of InChIs from a Mol file
    - `markinchitoinchilist <markinchi>` Generates a list of InChIs from a MarkInChI
    - `batchtest [-f filename]` Checks all test files work correctly. If `-f` is specified, can use a different file as the reference list.
    - `test <filename> [-l]` Checks whether the specified file works correctly. If `-l` is specified, will print generated lists of InChIs.
- To use this software from with GUI, simply run the program `MarkinchiGUI.py`


## Generating .mol files with Markush features
All molfiles used for testing and developing this software were produced using MarvinSketch 23.17. It may be possible to generate molfiles using other methods, however these have not yet been tested so may cause issues.

## Known issues
 - Variable attachments which are defined to have endpoints which have no hydrogen atoms on the core structure are known to cause errors. This includes molecules where the endpoints could be tautomeric, and, for example, an endpoint defined to be a quaternary carbon.


