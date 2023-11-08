from MarkinchiGenerator import MarkinchiGenerator
import os

cwd = os.getcwd()
file_dir = os.path.join(cwd, "molfiles")
molfile_folder = os.path.join(file_dir, "structures_for_testing")
reference_file = os.path.join(molfile_folder, "_testing_reference_list.txt")

with open(reference_file) as file:
    lines = file.readlines()

incorrect_files = []
files_read = 0

for line in lines:

    filename, reference_markinchi = tuple(line.split())
    molfile = os.path.join(molfile_folder, filename)
    
    markinchi_generator = MarkinchiGenerator()
    markinchi_generator.read_molfile(molfile)
    markinchi = markinchi_generator.generate_markinchi()

    if markinchi != reference_markinchi:
        incorrect_files.append((filename, markinchi, reference_markinchi))
    
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




