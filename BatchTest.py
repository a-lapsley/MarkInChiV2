from MarkinchiGenerator import MarkinchiGenerator
from MarkinchiParser import MarkinchiParser
import MarkinchiUtils as MUtils
import os

class BatchTest():
    def __init__(self) -> None:
        pass

    def test(self) -> tuple:
        cwd = os.getcwd()
        file_dir = os.path.join(cwd, "molfiles")
        molfile_folder = os.path.join(file_dir, "structures_for_testing")
        reference_file = os.path.join(molfile_folder, "_testing_reference_list.txt")

        with open(reference_file) as file:
            lines = file.readlines()

        incorrect_files = []
        incorrect_parses = []
        incorrect_enumerations = []
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
            new_mol, new_rgroups = parser.parse_markinchi()
            molblock = parser.get_molblock()

            markinchi_generator = MarkinchiGenerator()
            markinchi_generator.get_from_molblock(molblock)
            reparsed_markinchi = markinchi_generator.generate_markinchi()

            if reparsed_markinchi != markinchi:
                incorrect_parses.append((filename, reparsed_markinchi, markinchi))

            ref_mol, ref_rgroups = MUtils.parse_molfile(molfile)
            ref_list = MUtils.enumerate_markush_mol(ref_mol, ref_rgroups)
            ref_inchi_list = sorted(MUtils.inchis_from_mol_list(ref_list))

            new_list = MUtils.enumerate_markush_mol(new_mol, new_rgroups)
            new_inchi_list = sorted(MUtils.inchis_from_mol_list(new_list))

            if ref_inchi_list != new_inchi_list:
                incorrect_enumerations.append((molfile, new_inchi_list, ref_inchi_list))

            


            files_read += 1

        return files_read, incorrect_files, incorrect_parses, incorrect_enumerations

if __name__ == "__main__":
    files_read, incorrect_files, incorrect_parses, incorrect_enumerations = BatchTest().test()
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

    print("-------------------------------")
    if len(incorrect_enumerations) == 0:
        print("All MarkInChIs enumerated succesfully")
    else:
        print("The following files were not enumerated correctly:")
        for incorrect_file in incorrect_enumerations:
            print("Filename: \t\t%s" % incorrect_file[0])
            print("Found list: \t%s" % incorrect_file[1])
            print("Original list: \t%s" % incorrect_file[2])
            print("")




