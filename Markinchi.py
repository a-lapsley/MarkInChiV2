import sys
import os

from MarkinchiGenerator import MarkinchiGenerator
from MarkinchiParser import MarkinchiParser
from BatchTest import BatchTest
import MarkinchiUtils as MIUtils

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
    mol_list = MIUtils.enumerate_markush_mol(mol, rgroups)
    inchi_list = MIUtils.inchis_from_mol_list(mol_list)
    inchi_list = sorted(inchi_list)
    
    output = ""
    for inchi in inchi_list:
        output += inchi + "\n"
    

    return output

def inchi_list_from_markinchi(markinchi: str) -> str:

    parser = MarkinchiParser(markinchi)
    mol, rgroups = parser.parse_markinchi()
    mol_list = MIUtils.enumerate_markush_mol(mol, rgroups)
    inchi_list = MIUtils.inchis_from_mol_list(mol_list)
    inchi_list = sorted(inchi_list)

    output = ""
    for inchi in inchi_list:
        output += inchi + "\n"

    return output

def batch_test() -> str:

    test_result = BatchTest().test()
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
        filename = os.path.join(os.getcwd(), filename)
        file = open(filename)
        file.close()
        return filename
    except:
        pass

    return None

if __name__ == "__main__":
    argv = sys.argv[1:]

    try:
        command = argv[0]
    except:
        command = "help"
    
    try:
        argument = argv[1]
    except:
        argument = None

    if command.lower() in ("help"):

        help()

    if command.lower() in ("generate", "gen"):

        filename = find_file(argument)
        output = generate_from_molfile(filename)
        print(output)

    if command.lower() in ("markinchitomolblock", "mitomolblock"):

        markinchi = argument
        output = molblock_from_markinchi(markinchi)
        print(output)

    if command.lower() in ("molfiletoinchilist"):

        filename = find_file(argument)
        output = inchi_list_from_molfile(filename)
        print(output)

    if command.lower() in ("markinchitoinchilist", "mitoinchilist"):

        markinchi = argument
        output = inchi_list_from_markinchi(markinchi)
        print(output)

    if command.lower() in ("batchtest"):

        output = batch_test()
        print(output)




    




