from Bio import SeqIO, Entrez
import re, os, math, sys, time, csv, sqlite3, tkFileDialog, tempfile, ttk, ntpath
from Tkinter import *
from tkFileDialog import askopenfilename

def ncbi_email(file_name):
	global top
	top = Toplevel()
	top.title("NCBI Email")
	label_1 = Label(top,text="Please Enter NCBI Email:")
	label_1.grid(row=0,column=0)
	entry_1 = Entry(top)
	entry_1.grid(row=0,column=1)
	b1 = ttk.Button(top, text="Submit",command=lambda: choices(entry_1.get(),file_name))
	b1.grid(row=0,column=2)
	
def choices(ncbi_e,file_name):
	if all( [v1.get() == 1, v2.get() == 0] ):
		Tree_analysis_1(ncbi_e,file_name)
	elif all( [v1.get() == 1, v2.get() == 1] ):
		Tree_analysis_2(ncbi_e,file_name)
		


def Tree_analysis_1(ncbi_e,file_name):
    top.destroy()
    Entrez.email = ncbi_e
    # requests and opens file for reading
    my_seqs_file = file_name
    my_seqs = open(my_seqs_file,"r+")
    text = my_seqs.read()

    # sets up lists and variables
    gi_result = []
    locus = []
    type_result = []
    version_result = []
    version_number = []
    organism_name_1 = []
    organism_name_2 = []
    i = 0
    j = 0
    gb_data = ""


    #search file for matching terms (GI #, type, Version and Locus)
    gi_str = 'gi\|([\d\.]+)\|'
    type_str = '\|(ref|gb|emb|sp|tpg|dbj)\|'
    locus_str = 'LOCUS       (.+)             [0-9]'
    version_num_str = "gi\|\d+\|\w+\|\w+(.+?):"
    organism_str_1 = "ORGANISM  (.+)\s"
    organism_str_2 = "ORGANISM  \w+[\S]([\n]|.+)"
    gi_result = re.findall(gi_str,text)
    type_result = re.findall(type_str,text)
    version_number = re.findall(version_num_str, text)
    version_str = 'VERSION\     (.+)\.'
    output_file = open("report.csv","wb")
    csv_file = csv.writer(output_file)
    csv_file.writerow(["GI Number","Organism Name","Locus Name"])


    #Queries Entrez for genpept data and stores it in gb_data string
    for item in gi_result:
        handle = Entrez.efetch(db="protein", id = gi_result[i], rettype = "gb", retmode = "text")
        gb_data += str(handle.read())
        q = math.ceil((float(i)/(len(gi_result)))*100)
        time.sleep(0.1)
        z.set("\r" + str(q) + "% complete")
        root.update_idletasks()
        i += 1


    #Finds locus and version from genpept data stored in gb_data string
    z.set("")
    locus = re.findall(locus_str, gb_data)
    version_result = re.findall(version_str, gb_data)
    organism_name_1 = re.findall(organism_str_1, gb_data)
    organism_name_2 = re.findall(organism_str_2, gb_data)

    #Writes GIs|ORG_NAME|Locus (without spaces) in results file
    while j != i:
        csv_file.writerow([gi_result[j],str(organism_name_1[j]),(str(locus[j].replace(" ","")))])
        j += 1

    #Sets up and initializes text replacement in tree file
    i=0
    for item in gi_result:
        text = (text.replace("gi|" + str(gi_result[i]) + "|" + str(type_result[i]) + "|" + str(version_result[i]) + str(version_number[i]),
        organism_name_1[i].strip()))
        i += 1

    #Creates new tree file
    final_file = open("new_tree.txt","w")
    final_file.write(text)

    #Finalizes script
    z.set(str(i) + " sequences complete \nwritten to new_tree.txt")
    final_file.close()
    output_file.close()
def Tree_analysis_2(ncbi_e,file_name):
    top.destroy()
    Entrez.email = ncbi_e
    # requests and opens file for reading
    my_seqs_file = file_name
    my_seqs = open(my_seqs_file,"r+")
    text = my_seqs.read()

    # sets up lists and variables
    gi_result = []
    locus = []
    type_result = []
    version_result = []
    version_number = []
    organism_name_1 = []
    organism_name_2 = []
    i = 0
    j = 0
    gb_data = ""


    #search file for matching terms (GI #, type, Version and Locus)
    gi_str = 'gi\|([\d\.]+)\|'
    type_str = '\|(ref|gb|emb|sp|tpg|dbj)\|'
    locus_str = 'LOCUS       (.+)             [0-9]'
    length_str = 'LOCUS      .+             (.+) aa'
    version_num_str = "gi\|\d+\|\w+\|\w+(.+?):"
    organism_str_1 = "ORGANISM  (.+)\s"
    organism_str_2 = "ORGANISM  \w+[\S]([\n]|.+)"
    gi_result = re.findall(gi_str,text)
    type_result = re.findall(type_str,text)
    version_number = re.findall(version_num_str, text)
    version_str = 'VERSION\     (.+)\.'
    output_file = open("report.csv","wb")
    csv_file = csv.writer(output_file)
    csv_file.writerow(["GI Number","Organism Name","Locus Name"])


    #Queries Entrez for genpept data and stores it in gb_data string
    for item in gi_result:
        handle = Entrez.efetch(db="protein", id = gi_result[i], rettype = "gb", retmode = "text")
        gb_data += str(handle.read())
        q = math.ceil((float(i)/(len(gi_result)))*100)
        time.sleep(0.1)
        z.set("\r" + str(q) + "% complete")
        root.update_idletasks()
        i += 1


    #Finds locus and version from genpept data stored in gb_data string
    locus = re.findall(locus_str, gb_data)
    version_result = re.findall(version_str, gb_data)
    organism_name_1 = re.findall(organism_str_1, gb_data)
    organism_name_2 = re.findall(organism_str_2, gb_data)
    length_result = re.findall(length_str, gb_data)

    #Writes GIs|ORG_NAME|Locus (without spaces) in results file
    while j != i:
        csv_file.writerow([gi_result[j],str(organism_name_1[j]),(str(locus[j].replace(" ","")))])
        j += 1

    #Sets up and initializes text replacement in tree file
    i=0
    for item in gi_result:
        text = (text.replace("gi|" + str(gi_result[i]) + "|" + str(type_result[i]) + "|" + str(version_result[i]) + str(version_number[i]),
        organism_name_1[i].strip() + "_" + length_result[i]))
        i += 1

    #Creates new tree file
    final_file = open("new_tree.txt","w")
    final_file.write(text)

    #Finalizes script
    z.set(str(i) + " sequences complete \nwritten to new_tree.txt")
    final_file.close()
    output_file.close()
   
def file_check():
	root.fileName = askopenfilename()
        if not root.fileName:
            print("You didn't select a file")
        else:
				label_2 = Label(MainFrame, text = ntpath.basename(root.fileName) + " selected")
				label_2.grid(row=0,column=2,sticky=W)
	

class myGUI:
    def __init__(self, master):
		#variables
        global v1,v2,z
        v1 = IntVar()
        v2 = IntVar()
        z = StringVar()
        #frames
        global MainFrame
        MainFrame = Frame(master)
        MainFrame.pack()     
        #labels
        self.label_1 = Label(MainFrame, text ="Select Tree File")
        self.label_1.grid(row=0,column=0,sticky=W)
        self.label_3 = Label(MainFrame,text="Include:")
        self.label_3.grid(row=4,column=0)
        self.label_4 = Label(MainFrame, textvariable = z)
        self.label_4.grid(row=0,column=3)
        #buttons
        self.button_1 = ttk.Button(MainFrame, text="Upload Tree", command=lambda: file_check())
        self.button_1.grid(row=0,column=1,sticky=W)
        self.button_2 = ttk.Button(MainFrame, text="Rename Tree", command=lambda: ncbi_email(root.fileName))
        self.button_2.grid(row=5,columnspan=4,sticky=W+E)
        #checkbuttons
        self.c1 = ttk.Checkbutton(MainFrame, text="Organism Name", variable=v1)
        self.c1.grid(row=4,column=1)
        self.c3 = ttk.Checkbutton(MainFrame, text="Length", variable=v2)
        self.c3.grid(row=4,column=2)
        master.resizable(0,0)
        master.title("Rename Protein Tree v1.0")
 

    def createRectangle(self,master,gridrow,gridcolumn,x1,y1,x2,y2):
        rect = Canvas(master,width=700,height=500)
        rect.pack()
        rect.create_rectangle(x1,y1,x2,y2)





def main():
    global root
    root = Tk()
    window = myGUI(root)
    root.mainloop()










if __name__ == '__main__': main()
