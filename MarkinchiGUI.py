import tkinter as tk
from tkinter.filedialog import askopenfilename

import Markinchi as MI

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

        self.btn_markinchi_to_molblock = tk.Button(
            self.window,
            text="Generate Molblock from MarkInChI",
            font=MAIN_FONT,
            command=self.markinchi_to_molblock
        )
        self.btn_markinchi_to_molblock.pack(padx=10, pady=10)

        self.btn_molfile_to_inchi_list = tk.Button(
            self.window,
            text="List of InChIs from Molfile",
            font=MAIN_FONT,
            command=self.molfile_to_inchi_list
        )
        self.btn_molfile_to_inchi_list.pack(padx=10, pady=10)

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
        output = MI.generate_from_molfile(filename)
        
        popup = tk.Tk()
        popup.title("Generated MarkInChI")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, output)
        textbox.pack()
        popup.mainloop()

    def markinchi_to_molblock(self) -> None:
        input_window = TextInputWindow("Enter MarkInChI:")
        markinchi = input_window.get_input()

        molblock = MI.molblock_from_markinchi(markinchi)

        popup = tk.Tk()
        popup.title("Generated Molblock")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, molblock)
        textbox.pack()
        popup.mainloop()

    def molfile_to_inchi_list(self) -> None:
        filename = askopenfilename()
        output = MI.inchi_list_from_molfile(filename)
        
        popup = tk.Tk()
        popup.title("List of InChIs")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, output)
        textbox.pack()
        popup.mainloop()

    def markinchi_to_inchi_list(self) -> None:

        input_window = TextInputWindow("Enter MarkInChI:")
        markinchi = input_window.get_input()

        output = MI.inchi_list_from_markinchi(markinchi)

        popup = tk.Tk()
        popup.title("List of InChIs")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, output)
        textbox.pack()
        popup.mainloop()

        

    def batch_test(self) -> None:
        output = MI.batch_test()

        popup = tk.Tk()
        popup.title("Batch Test Results")
        textbox = tk.Text(popup)
        textbox.insert(tk.END, output)

        textbox.pack()
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