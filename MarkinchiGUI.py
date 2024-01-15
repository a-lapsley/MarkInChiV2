import tkinter as tk
from tkinter.filedialog import askopenfilename

from rdkit import Chem
from rdkit.Chem import Draw

from MarkinchiGenerator import MarkinchiGenerator
from MarkinchiParser import MarkinchiParser
from BatchTest import BatchTest
import MarkinchiUtils as MIUtils

TITLE_FONT = ("Arial", 25)
MAIN_FONT = ("Arial", 18)

class MarkinchiGUI:
    
    def __init__(self) -> None:
        self.window = tk.Tk()
        self.window.title("MarkInChI Tools")
        self.window.geometry("500x600")

        self.lbl_title = tk.Label(
            self.window, text="MarkInChI Tools", font=TITLE_FONT)
        self.lbl_title.pack(padx=20, pady=20)

        self.btn_generate_from_molfile = tk.Button(
            self.window,
            text="Generate MarkInChI from Molfile",
            font=MAIN_FONT,
            command=self.generate_from_molfile)
        self.btn_generate_from_molfile.pack(padx=10, pady=10)

        self.btn_generate_from_molblock = tk.Button(
            self.window,
            text="Generate MarkInChI from Molblock",
            font=MAIN_FONT,
            command=self.generate_from_molblock)
        self.btn_generate_from_molblock.pack(padx=10, pady=10)

        self.btn_markinchi_to_molblock = tk.Button(
            self.window,
            text="Generate Molblock from MarkInChI",
            font=MAIN_FONT,
            command=self.markinchi_to_molblock
        )
        self.btn_markinchi_to_molblock.pack(padx=10, pady=10)

        self.btn_markinchi_to_inchi_list = tk.Button(
            self.window,
            text="List of InChIs from MarkInChI",
            font=MAIN_FONT,
            command=self.markinchi_to_inchi_list
        )
        self.btn_markinchi_to_inchi_list.pack(padx=10, pady=10)

        self.btn_batch_test = tk.Button(
            self.window,
            text="Batch test all test structures",
            font=MAIN_FONT,
            command=self.batch_test
        )
        self.btn_batch_test.pack(padx=10, pady=10)


        self.window.mainloop()
    
    def generate_from_molfile(self) -> None:
        filename = askopenfilename()
        markinchi_generator = MarkinchiGenerator()
        markinchi_generator.get_from_molfile(filename)
        markinchi = markinchi_generator.generate_markinchi()
        
        popup = tk.Tk()
        popup.title("Generated MarkInChI")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, markinchi)
        textbox.pack()
        popup.mainloop()

    def generate_from_molblock(self) -> None:
        
        input_window = TextInputWindow("Enter Molblock:")
        molblock = input_window.get_input()
        markinchi_generator = MarkinchiGenerator()
        markinchi_generator.get_from_molblock(molblock)
        markinchi = markinchi_generator.generate_markinchi()
        
        popup = tk.Tk()
        popup.title("Generated MarkInChI")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, markinchi)
        textbox.pack()
        popup.mainloop()

    def markinchi_to_molblock(self) -> None:
        input_window = TextInputWindow("Enter MarkInChI:")
        markinchi = input_window.get_input()

        markinchi_parser = MarkinchiParser(markinchi)
        markinchi_parser.parse_markinchi()
        molblock = markinchi_parser.get_molblock()

        popup = tk.Tk()
        popup.title("Generated Molblock")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, molblock)
        textbox.pack()
        popup.mainloop()

    def markinchi_to_inchi_list(self) -> None:

        input_window = TextInputWindow("Enter MarkInChI:")
        markinchi = input_window.get_input()

        markinchi_parser = MarkinchiParser(markinchi)
        mol, rgroups = markinchi_parser.parse_markinchi()
        mol_list = MIUtils.enumerate_markush_mol(mol, rgroups)
        inchi_list = MIUtils.inchis_from_mol_list(mol_list)

        popup = tk.Tk()
        popup.title("List of InChIs")
        textbox = tk.Text(popup)
        for inchi in inchi_list:
            textbox.insert(tk.END, inchi + "\n")
        textbox.pack()
        popup.mainloop()

    def batch_test(self) -> None:
        files_read, incorrect_files, incorrect_parses = BatchTest().test()

        popup = tk.Tk()
        popup.title("Batch Test Results")
        tb = tk.Text(popup)
        tb.insert(tk.END, "-------------------------------\n")
        tb.insert(tk. END, "Read %i .mol files \n" % files_read)

        if len(incorrect_files) == 0:
            tb.insert(tk.END,"All generated MarkInChIs matched the reference\n")
        else:
            tb.insert(tk.END,"The following files caused a mismatch against the reference:\n")
            for incorrect_file in incorrect_files:
                tb.insert(tk.END,"Filename: \t\t%s\n" % incorrect_file[0])
                tb.insert(tk.END,"Found MarkInChI: \t%s\n" % incorrect_file[1])
                tb.insert(tk.END,"Reference MarkInChI: \t%s\n" % incorrect_file[2])
                tb.insert(tk.END,"\n")

        tb.insert(tk.END, "-------------------------------\n")
        if len(incorrect_parses) == 0:
            tb.insert(tk.END, "All MarkInChIs parsed succesfully\n")
        else:
            tb.insert(tk.END, "The following files were not parsed correctly:\n")
            for incorrect_file in incorrect_parses:
                tb.insert(tk.END, "Filename: \t\t%s\n" % incorrect_file[0])
                tb.insert(tk.END, "Reparsed MarkInChI: \t%s\n" % incorrect_file[1])
                tb.insert(tk.END, "Original MarkInChI: \t%s\n" % incorrect_file[2])
                tb.insert(tk.END, "\n")

        tb.pack()
        popup.mainloop()

class TextInputWindow():

    def __init__(self, title: str="Input Text") -> None:
        self.text = None
        self.window = tk.Tk()
        self.window.title(title)

        self.textbox = tk.Text(self.window)
        self.textbox.pack()

        self.btn_enter = tk.Button(
            self.window,
            text="Enter",
            font=MAIN_FONT,
            command=self.get_text_from_textbox
        )
        self.btn_enter.pack()

        self.window.mainloop()


    def get_text_from_textbox(self) -> None:
        self.text = self.textbox.get("1.0","end").strip()
        self.window.quit()
        self.window.destroy()
        

    def get_input(self) -> str:
        return self.text 
        
MarkinchiGUI()