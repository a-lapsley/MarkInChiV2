from MarkinchiGenerator import MarkinchiGenerator
from MarkinchiParser import MarkinchiParser
from MolEnumerator import MolEnumerator
import MarkinchiUtils as MIUtils
import os

class BatchTest():
    def __init__(self) -> None:
        pass

    def test(self, ref_list: str) -> tuple:
        cwd = os.getcwd()
        file_dir = os.path.join(cwd, "molfiles")
        molfile_folder = os.path.join(file_dir, "structures_for_testing")
        reference_file = os.path.join(
            molfile_folder, ref_list
            )

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
                incorrect_files.append(
                    (filename, markinchi, reference_markinchi)
                    )
            
            try:
                parser = MarkinchiParser(markinchi)
                new_mol, new_rgroups = parser.parse_markinchi()
                molblock = parser.get_molblock()

                markinchi_generator = MarkinchiGenerator()
                markinchi_generator.get_from_molblock(molblock)
                reparsed_markinchi = markinchi_generator.generate_markinchi()
            except:
                reparsed_markinchi = "Unable to parse this MarkInChI"

            if reparsed_markinchi != markinchi:
                incorrect_parses.append(
                    (filename, reparsed_markinchi, markinchi)
                    )

            ref_mol, ref_rgroups = MIUtils.parse_molfile(molfile)
            ref_enumerator = MolEnumerator(ref_mol, ref_rgroups)
            ref_inchi_list = ref_enumerator.get_inchi_list()

            try:
                new_enumerator = MolEnumerator(new_mol, new_rgroups)
                new_inchi_list = new_enumerator.get_inchi_list()
            except:
                new_inchi_list = "Unable to generate InChI list"

            if ref_inchi_list != new_inchi_list:
                incorrect_enumerations.append(
                    (molfile, new_inchi_list, ref_inchi_list)
                    )

            


            files_read += 1

        return (
            files_read, incorrect_files, 
            incorrect_parses, incorrect_enumerations
            )

