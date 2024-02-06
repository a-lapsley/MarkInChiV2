import sys
import os
import getopt

from MarkinchiGenerator import MarkinchiGenerator
from MarkinchiParser import MarkinchiParser
from MolEnumerator import MolEnumerator
from BatchTest import BatchTest
import MarkinchiUtils as MIUtils

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def generate_from_molfile(filename: str) -> str:
    
    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molfile(filename)
    markinchi = markinchi_generator.generate_markinchi()

    return markinchi

def molblock_from_markinchi(markinchi: str) -> str:

    markinchi_parser = MarkinchiParser(markinchi)
    markinchi_parser.parse_markinchi()
    molblock = markinchi_parser.get_molblock()

    return molblock

def inchi_list_from_molfile(filename: str) -> str:

    mol, rgroups = MIUtils.parse_molfile(filename)
    enumerator = MolEnumerator(mol, rgroups)
    inchi_list = enumerator.get_inchi_list()
    inchi_list = sorted(inchi_list)
    
    output = ""
    for inchi in inchi_list:
        output += inchi + "\n"
    

    return output

def inchi_list_from_markinchi(markinchi: str) -> str:

    parser = MarkinchiParser(markinchi)
    mol, rgroups = parser.parse_markinchi()
    enumerator = MolEnumerator(mol, rgroups)
    inchi_list = enumerator.get_inchi_list()
    inchi_list = sorted(inchi_list)

    output = ""
    for inchi in inchi_list:
        output += inchi + "\n"

    return output

def batch_test(filename: str = "_testing_reference_list.txt") -> str:

    test_result = BatchTest().test(filename)
    files_read = test_result[0]
    incorrect_files = test_result[1]
    incorrect_parses = test_result[2]
    incorrect_enumerations = test_result[3]

    output = "-------------------------------\n"
    output += "Read %i .mol files \n" % files_read

    if len(incorrect_files) == 0:
        output += "All generated MarkInChIs matched the reference\n"
    else:
        output += ("The following files caused a mismatch"
                   "against the reference:\n")
        for incorrect_file in incorrect_files:
            output += "Filename: \t\t%s\n" % incorrect_file[0]
            output += "Found MarkInChI: \t%s\n" % incorrect_file[1]
            output += "Reference MarkInChI: \t%s\n\n" % incorrect_file[2]

    output += "-------------------------------\n"
    
    if len(incorrect_parses) == 0:
        output += "All generated MarkInChIs parsed succesfully\n"
    else:
        output += ("The following files were not parsed correctly:\n")
        for incorrect_parse in incorrect_parses:
            output += "Filename: \t\t%s\n" % incorrect_parse[0]
            output += "Reparsed MarkInChI: \t%s\n" % incorrect_parse[1]
            output += "Original MarkInChI: \t%s\n" % incorrect_parse[2]

    output += "-------------------------------\n"

    if len(incorrect_enumerations) == 0:
        output += "All generated MarkInChIs enumerated succesfully\n"
    else:
        output += ("The following MarkInChIs were not enumerated correctly:\n")
        for incorrect_enumeration in incorrect_enumerations:
            output += "Filename: \t\t%s\n" % incorrect_enumeration[0]
            new_list = incorrect_enumeration[1]
            ref_list = incorrect_enumeration[2]
            output += "Length of new list: %i\n" % len(new_list)
            output += "Length of reference list: %i\n" % len(ref_list)
            output += "InChIs in the new list but not the reference:\n"
            for item in new_list:
                if item not in ref_list:
                    output += item + "\n"
            output += "\n"
            output += "InChIs in the reference list but not the new:\n"
            for item in ref_list:
                if item not in new_list:
                    output += item + "\n"
            output += "\n"


    return output

def test_file(filename: str, print_lists: bool = False) -> str:

    output = ""

    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molfile(filename)
    markinchi = markinchi_generator.generate_markinchi()

    output += "Generated MarkInChI:\n"
    output += markinchi +  "\n"

    parser = MarkinchiParser(markinchi)
    new_mol, new_rgroups = parser.parse_markinchi()
    molblock = parser.get_molblock()

    new_markinchi_generator = MarkinchiGenerator()
    new_markinchi_generator.get_from_molblock(molblock)
    new_markinchi = new_markinchi_generator.generate_markinchi()

    if new_markinchi == markinchi:
        output += "MarkInChI succesfully parsed and regenerated\n"
    else:
        output += "Regenerated MarkInChI did not match the original:\n"
        output += new_markinchi

    new_enumerator = MolEnumerator(new_mol, new_rgroups)
    new_inchi_list = new_enumerator.get_inchi_list()

    ref_mol, ref_rgroups = MIUtils.parse_molfile(filename)
    ref_enumerator = MolEnumerator(ref_mol, ref_rgroups)
    ref_inchi_list = ref_enumerator.get_inchi_list()
    
    if new_inchi_list == ref_inchi_list:
        output += "Correctly enumerated MarkInChI\n"
        if len(new_inchi_list) == 0:
            output += "WARNING: Generated InChI list contained 0 structures\n"
    else:
        output += "MarkInChI was not enumerated correctly\n"
        output += "Original list contained %i InChIs\n" % len(ref_inchi_list)
        output += "New list contained %i InChIs\n" % len(new_inchi_list)

    if print_lists:
        output += "Original list: \n"
        for inchi in ref_inchi_list:
            output += inchi + "\n"
        output += "\nNew list: \n"
        for inchi in new_inchi_list:
            output += inchi + "\n"

        
    return output

def help() -> None:

    print("Available commands:")
    print("help\t\t\t\t\tDisplays this help menu")
    print("generate <filename>\t\t\tGenerates a MarkInChI from a Mol file")
    print("markinchitomolblock <markinchi>\t\t"
          "Generates a Molblock from a MarkInChI")
    print("molfiletoinchilist <filename>\t\t"
          "Generates a list of InChIs from a Mol file")
    print("markinchitoinchilist <markinchi>\t"
          "Generates a list of InChIs from a MarkInChI")
    print("batchtest\t\t\t\tChecks all test files work correctly")
    
def find_file(filename: str) -> str:

    try:
        file = open(filename)
        file.close()
        return filename
    except:
        pass

    try:
        full_filename = os.path.join(os.getcwd(), filename)
        file = open(full_filename)
        file.close()
        return full_filename
    except:
        pass

    raise FileNotFoundError("Unable to open %s or %s" % (
        filename, full_filename
    ))

if __name__ == "__main__":
    argv = sys.argv[1:]

    try:
        command = argv[0]
    except:
        command = "help"
    
    
    argv = argv[1:]




    if command.lower() in ("help"):

        help()

    if command.lower() in ("generate", "gen"):

        filename = argv[0]
        filename = find_file(filename)
        output = generate_from_molfile(filename)
        print(output)

    if command.lower() in ("markinchitomolblock", "mitomolblock"):

        markinchi = argv[0]
        output = molblock_from_markinchi(markinchi)
        print(output)

    if command.lower() in ("molfiletoinchilist"):

        filename = argv[0]
        filename = find_file(filename)
        output = inchi_list_from_molfile(filename)
        print(output)

    if command.lower() in ("markinchitoinchilist", "mitoinchilist"):

        markinchi = argv[0]
        output = inchi_list_from_markinchi(markinchi)
        print(output)

    if command.lower() in ("batchtest", ""):
        
        opts, args = getopt.getopt(argv, "f:",["file="])

        filename = "_testing_reference_list.txt"
        for opt, arg in opts:
            if opt in ["-f", "--file"]:
                filename = arg

       
        output = batch_test(filename)


        print(output)

    if command.lower() in ("test", ""):

        filename = argv[0]
    
        opts, args = getopt.getopt(argv[1:],"l",["list"])

        list_inchis = False
        for opt, arg in opts:
            if opt in ["-l", "--list"]:
                list_inchis = True

        filename = find_file(filename)

        output = test_file(filename, list_inchis)
        print(output)




    




