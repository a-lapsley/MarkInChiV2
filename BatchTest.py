from MarkinchiGenerator import MarkinchiGenerator
from MarkinchiParser import MarkinchiParser
import os

cwd = os.getcwd()
file_dir = os.path.join(cwd, "molfiles")
molfile_folder = os.path.join(file_dir, "structures_for_testing")
reference_file = os.path.join(molfile_folder, "_testing_reference_list.txt")

with open(reference_file) as file:
    lines = file.readlines()

incorrect_files = []
incorrect_parses = []
files_read = 0

for line in lines:

    filename, reference_markinchi = tuple(line.split())
    molfile = os.path.join(molfile_folder, filename)
    print("Checking %s" % filename)
    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molfile(molfile)
    markinchi = markinchi_generator.generate_markinchi()

    if markinchi != reference_markinchi:
        incorrect_files.append((filename, markinchi, reference_markinchi))
    
    parser = MarkinchiParser(markinchi)
    parser.parse_markinchi()
    molblock = parser.get_molblock()

    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.get_from_molblock(molblock)
    reparsed_markinchi = markinchi_generator.generate_markinchi()

    if reparsed_markinchi != markinchi:
        incorrect_parses.append((filename, reparsed_markinchi, markinchi))


    files_read += 1

print("-------------------------------")
print("Read %i .mol files \n" % files_read)

if len(incorrect_files) == 0:
    print("All generated MarkInChIs matched the reference")
else:
    print("The following files caused a mismatch against the reference:")
    for incorrect_file in incorrect_files:
        print("Filename: \t\t%s" % incorrect_file[0])
        print("Found MarkInChI: \t%s" % incorrect_file[1])
        print("Reference MarkInChI: \t%s" % incorrect_file[2])
        print("")

print("-------------------------------")
if len(incorrect_parses) == 0:
    print("All MarkInChIs parsed succesfully")
else:
    print("The following files were not parsed correctly:")
    for incorrect_file in incorrect_parses:
        print("Filename: \t\t%s" % incorrect_file[0])
        print("Reparsed MarkInChI: \t%s" % incorrect_file[1])
        print("Original MarkInChI: \t%s" % incorrect_file[2])
        print("")




