# PyPore  
# 
# Idea:		GUI for porE 
# 
# Author: 	Sebastian Schwalbe 
# Dates:  	
#		06.01.2020 -- init, pore functions  
#               09.01.2020 -- use new pore.f90 files 
#                             build user input 
#                             pore and get_psd are working 
#               10.01.2020 -- include all Fortran porE (original,subgrid,window) 
#                             include grid_density input 
#                             include some help information 
# TODO:         add visulation of structure and pores may use pygui 

# June 1st, 2020  -> restructuring due to new fortran routines

#from porE import pore
#from porE_subgrid import pore as pore_subgrid
##from porE_window import pore as pore_window 
#from get_PSD import * 

import pore
# module pore_options
#   - subroutines: osa, gpa_fullgrid, gpa_gridpera, do_gpa (this executes the GPA evaluation), get_psd

import os 
import math
import numpy as np
from tkinter import *
from tkinter import scrolledtext as st
from tkinter.filedialog import askopenfilename
from ase.io import read 
from ase.visualize import view

from tkinter.ttk import Frame

# define some abbreviations
osa          = pore.porosity.osa
gpa          = pore.porosity.do_gpa
gpa_FullGrid = pore.porosity.gpa_fullgrid
gpa_GridPerA = pore.porosity.gpa_gridpera
get_PSD      = pore.psd.get_psd

# the actual class 
class PyPore(Frame):
    # Predefined MOFs
    settings = {'ud' : {'description' :'user defined','file':'pypore.xyz'},
                'do' : {'description' :'DUT-8(Ni) open','file':'dut_8_open.xyz'},
                'vo' : {'description' :'DUT-8(Ni) open vcrelax','file':'dut_8_open_vcrelax.xyz'},
                'dc' : {'description' :'DUT-8(Ni) closed','file':'dut_8_closed.xyz'},
                'vc' : {'description' :'DUT-8(Ni) closed vcrelax','file':'dut_8_closed_vcrelax.xyz'},
                'u6' : {'description' :'UiO-66','file':'uio66.xyz'},
                'u7' : {'description' :'UiO-67','file':'uio67.xyz'},
                'u8' : {'description' :'UiO-68','file':'uio68.xyz'},
                'm5' : {'description' :'MOF-5','file':'mof5.xyz'},
                'ir' : {'description' :'IRMOF-10','file':'irmof10.xyz'},
                'm2' : {'description' :'MOF210','file':'mof210.xyz'},
                'h1' : {'description' :'HKUST-1, open Cu sites','file':'hkust1.xyz'},
                'ho' : {'description' :'HKUST-1, O-Cu-Cu-O','file':'hkust1_with_O.xyz'},
                'c6' : {'description' :'C60@MOF','file':'c60_MOF.xyz'},
                'be' : {'description' :'Benzene, opt','file':'benzene.xyz'},
                'b2' : {'description' :'Benzene, exp','file':'benzene_exp.xyz'},
                'bc' : {'description' :'Benzene, C only','file':'benzene_Conly.xyz'},
                'ha' : {'description' :'H atom','file':'h_atom.xyz'}}

    def __init__(self):
        
        # Generate the general gui setup 
        master = Tk()
        self.master = master 
        self.master.title("PyPorE")
        # Predefined porE xyz files 
        self.dirname = 'structures/xyz/'       

        # Action buttons 
        #b1 = Button(self.master, text='Open', command=self.open_file).grid(column=2,row=3, sticky=W, pady=4)
        Button(self.master, text='Porosity', command=self.get_pore).grid(column=3,row=3, sticky=W, pady=4)
        Button(self.master, text='PSD', command=self.get_psd).grid(column=4,row=3, sticky=W, pady=4)
        Button(self.master, text='Quit', command=self.master.quit).grid(column=6,row=3, sticky=W, pady=4)
      
        # Inital values 
        self.probe_r = 1.2 
        self.grid_a = 5
        self.grid_b = 5 
        self.grid_c = 5 
        self.init_newwin()
        self.l1 = None
        self.cell = None

        # Gui 
        mainloop()
    
    def get_pore(self):
        # main function for porosity calculation 
        self.newwin.destroy()
        self.init_newwin()
        # Inital text/ Help 
        s = ''' PyPorE: porosity \n
                \t structure \n 
                \t\t\t ud -- user defined; use "open" to select the structure \n 
                \t\t\t do -- predefined DUT-8(Ni) open, see others as well \n 
                \t methods \n
                \t\t\t 1 -- overlaping spheres approach (OSA) \n 
                \t\t\t 2 -- grid approach \n
                \t\t\t\t    enter value and click in the field to update the value\n 
                \t executable: \n
                \t\t\t GPA subgrid -- speedup for grid approach, can also calculate pore windows (using output_PSD file) \n
                \t Authors: \n 
                \t\t\t Kai Trepte (Fortran, core routines) \n 
                \t\t\t Sebastian Schwalbe (Python, GUI)'''
        self.text.insert(END,s)
        self.text.pack()
        # Search input fields 
        # What property to search for: when, where or ...
        Label(self.controls, text="Select structure").grid(row=1,column=1)
        Label(self.controls, text="Select method").grid(row=3,column=1)
        Label(self.controls, text="e.g. do").grid(row=1,column=3)
        #Label(self.controls, text="e.g. 1").grid(row=3,column=3)

        # Entry fields 
        e1 = StringVar(self.newwin)
        self.e1 = e1
        self.e1.set("ud") # default value 
        w1= OptionMenu(self.controls, self.e1, 'ud','do','vo','dc','vc','u6','u7','m5','ir','m2','h1','ho','c6','be','b2','bc','ha',command=self.refresh_pore)

        self.w1 = w1

        e2 = Entry(self.controls)
        e2.insert(END,'1') # default value 
        e3 = Entry(self.controls) 
        e3.insert(END,self.settings[self.e1.get()]['description'])
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3 
        
        exe_pore = StringVar()
        exe_pore.set('subgrid')
        self.exe_pore = exe_pore
        #w2= OptionMenu(self.controls, self.exe_pore, 'subgrid',command=self.refresh_pore)
        #self.w2 = w2 

        self.w1.grid(row=1, column=2)
        self.e2.grid(row=3, column=2)
        self.e3.grid(row=2, column=2)
        #self.w2.grid(row=3, column=3)

        # Action buttons 
        if self.e1.get() == 'ud':
            b3 = Button(self.controls, text='Open', command=self.user_input)
            b3.grid(row=4,column=1)
            self.b3 = b3 
        b1 = Button(self.controls, text='OK', command=self.check_pore)
        b1.grid(row=4,column=2)
        self.b1 = b1 
        self.newwin.update()
        self.newwin.deiconify()

    def refresh_pore(self,event):
        # refresh the gui 
        # the user can select the 1st method then the 2nd one 
        # if method 1 we need not the additional input fields 
        # if method 2 we need the additional input fields 
        if self.e2.get() == '1':
            try:
                # we ever have these elements 
                self.b3.destroy()
                self.b1.destroy()
                # only for method 2 we have these elements 
                self.l4.destroy()  
                self.l5.destroy() 
                self.l6.destroy() 
                self.l7.destroy() 
                self.l8.destroy()
                self.e4.destroy() 
                self.e5.destroy() 
                self.e6.destroy() 
                self.e7.destroy() 
                self.e8.destroy()
            except: 'Nothing' 
        if self.e2.get() == '2':
            try:
                self.b3.destroy()
                self.b1.destroy()
            except: 'Nothing'
        self.e3.delete(0,END)
        self.e3.insert(END,self.settings[self.e1.get()]['description'])
        if self.e2.get() == '1':
        
            b1 = Button(self.controls, text='OK', command=self.check_pore)
            b1.grid(row=4,column=2)
            self.b1 = b1
            if self.e1.get() == 'ud':
                b3 = Button(self.controls, text='Open', command=self.user_input)
                b3.grid(row=4,column=1)
                self.b3 = b3
        if self.e2.get() == '2':
            b1 = Button(self.controls, text='OK', command=self.cmd_pore)
            b1.grid(row=9,column=2)
            self.b1 = b1
            if self.e1.get() == 'ud':
                b3 = Button(self.controls, text='Open', command=self.user_input)
                b3.grid(row=9,column=1)
                self.b3 = b3

        if self.e1.get() != 'ud':
            self.b3.destroy()

    def refresh_psd(self,event):
        # refresh the gui 
        # if ud we need a open button 
        # if not ud we need no open button 
        if self.e1.get() == 'ud':
            b3 = Button(self.controls, text='Open', command=self.user_input)
            b3.grid(row=4,column=1)
            self.b3 =b3
        if self.e1.get() != 'ud':
            self.b3.destroy()

    def check_pore(self):
        # check the input 
        # if method 1 is chosen we can run porE 
        # if method 2 is chosen we need additional input 
        self.refresh_pore(self)
        self.get_target()
        if self.target[1] == '1':
            self.cmd_pore()
        if self.target[1] == '2':
            self.b1.destroy() 
            if self.e1.get() == 'ud':
                self.b3.destroy() 
            # labels 
            l4 = Label(self.controls, text="probe_r")
            l4.grid(row=4,column=1)
            self.l4 = l4 
            l5 = Label(self.controls, text="grid_a")
            l5.grid(row=5,column=1)
            self.l5 = l5 
            l6 = Label(self.controls, text="grid_b")
            l6.grid(row=6,column=1)
            self.l6 = l6 
            l7 = Label(self.controls, text="grid_c")
            l7.grid(row=7,column=1)
            self.l7 = l7 
            l8 = Label(self.controls, text="grid_density")
            l8.grid(row=8,column=1)
            self.l8 = l8 

            # entries 
            e4 = Entry(self.controls)
            e4.insert(END,'1.2') # default value
            self.e4 = e4
            
            grid_a = StringVar()
            self.grid_a = grid_a 
            self.grid_a.set('5') 
            e5 = Entry(self.controls)
            e5.insert(END,'5') # default value
            self.e5 = e5
            # Enter value and clicking in the field update the value 
            self.e5.bind("<Button-1>", self.grid2grid_density)


            grid_b = StringVar()
            self.grid_b = grid_b
            self.grid_b.set('5')
            e6 = Entry(self.controls)
            e6.insert(END,'5') # default value
            self.e6 = e6
            # Enter value and clicking in the field update the value
            self.e6.bind("<Button-1>", self.grid2grid_density)

            grid_c = StringVar()
            self.grid_c = grid_c
            self.grid_c.set('5')
            e7 = Entry(self.controls)
            e7.insert(END,'5') # default value
            self.e7 = e7
            # Enter value and clicking in the field update the value
            self.e7.bind("<Button-1>", self.grid2grid_density)

            g = StringVar() 
            self.g = g 
            self.g.set('5')
            g.trace("w", self.grid_density2grid)
            
            e8 = Entry(self.controls)
            e8.insert(END,self.g.get()) # default value
            self.e8 = e8
            # Enter value and clicking in the field update the value
            self.e8.bind("<Button-1>", self.grid_density2grid)

            self.e4.grid(row=4, column=2)
            self.e5.grid(row=5, column=2)
            self.e6.grid(row=6, column=2)
            self.e7.grid(row=7, column=2)
            self.e8.grid(row=8, column=2)
            
            if self.e1.get() == 'ud':
                b3 = Button(self.controls, text='Open', command=self.user_input)    
                b3.grid(row=9,column=1)
                self.b3 =b3 
            b1 = Button(self.controls, text='OK', command=self.cmd_pore)
            b1.grid(row=9,column=2)
            self.b1 = b1
            self.newwin.update()
            self.newwin.deiconify()

            self.controls.grid_rowconfigure(10, weight=1)
            self.controls.grid_columnconfigure(1, weight=1)

    def grid_density2grid(self,event):
        # grid_a, grid_b, grid_c to grid_density 
        if self.e1.get() != 'ud':
            f = open(self.dirname+self.settings[self.e1.get()]['file'],'r')
        if self.e1.get() == 'ud':
            f = open(self.settings[self.e1.get()]['file'],'r')
        ll = f.readlines()
        f.close()
        cell = np.zeros([3,3])
        cell[0,0] = ll[1].split()[0]
        cell[0,1] = ll[1].split()[1]
        cell[0,2] = ll[1].split()[2]
        cell[1,0] = ll[1].split()[3]
        cell[1,1] = ll[1].split()[4]
        cell[1,2] = ll[1].split()[5]
        cell[2,0] = ll[1].split()[6]
        cell[2,1] = ll[1].split()[7]
        cell[2,2] = ll[1].split()[8]
        self.cell = cell 
        g = float(self.e8.get())
        grid_a = math.ceil(g*np.sqrt(self.cell[0,0]**2+self.cell[0,1]**2+self.cell[0,2]**2)) # +1 ? 
        grid_b = math.ceil(g*np.sqrt(self.cell[1,0]**2+self.cell[1,1]**2+self.cell[1,2]**2))
        grid_c = math.ceil(g*np.sqrt(self.cell[2,0]**2+self.cell[2,1]**2+self.cell[2,2]**2))
        
        self.e5.delete(0, END)
        self.e5.insert(END,grid_a) 
        self.e6.delete(0, END)
        self.e6.insert(END,grid_b)
        self.e7.delete(0, END)
        self.e7.insert(END,grid_c)
    
    def grid2grid_density(self,event):
        # grid_density to grid_a, grid_b, grid_c 
        if self.e1.get() != 'ud':
            f = open(self.dirname+self.settings[self.e1.get()]['file'],'r')
        if self.e1.get() == 'ud':
            f = open(self.settings[self.e1.get()]['file'],'r')
        ll = f.readlines()
        f.close()
        cell = np.zeros([3,3])
        cell[0,0] = ll[1].split()[0]
        cell[0,1] = ll[1].split()[1]
        cell[0,2] = ll[1].split()[2]
        cell[1,0] = ll[1].split()[3]
        cell[1,1] = ll[1].split()[4]
        cell[1,2] = ll[1].split()[5]
        cell[2,0] = ll[1].split()[6]
        cell[2,1] = ll[1].split()[7]
        cell[2,2] = ll[1].split()[8]
        self.cell = cell
        g_a = float(self.e5.get())/math.ceil(np.sqrt(self.cell[0,0]**2+self.cell[0,1]**2+self.cell[0,2]**2)) 
        g_b = float(self.e6.get())/math.ceil(np.sqrt(self.cell[1,0]**2+self.cell[1,1]**2+self.cell[1,2]**2))
        g_c = float(self.e7.get())/math.ceil(np.sqrt(self.cell[2,0]**2+self.cell[2,1]**2+self.cell[2,2]**2))
        # debug output 
        #print(g_a)
        #print(g_b)
        #print(g_c) 

        # Approximation 
        new_g = g_a 
        self.e8.delete(0, END)
        self.e8.insert(END,new_g)

    def cmd_pore(self):
        # The Fortran call 
        self.refresh_pore(self)
        self.get_target()
        # KT: get correct structure inputs
        structs  = {'ud' : 'pypore.xyz',
                    'do' : self.dirname+'dut_8_open.xyz',
                    'vo' : self.dirname+'dut_8_open_vcrelax.xyz',
                    'dc' : self.dirname+'dut_8_closed.xyz',
                    'vc' : self.dirname+'dut_8_closed_vcrelax.xyz',
                    'u6' : self.dirname+'uio66.xyz',
                    'u7' : self.dirname+'uio67.xyz',
                    'u8' : self.dirname+'uio68.xyz',
                    'm5' : self.dirname+'mof5.xyz',
                    'ir' : self.dirname+'irmof10.xyz',
                    'm2' : self.dirname+'mof210.xyz',
                    'h1' : self.dirname+'hkust1.xyz',
                    'ho' : self.dirname+'hkust1_with_O.xyz',
                    'c6' : self.dirname+'c60_MOF.xyz',
                    'be' : self.dirname+'benzene.xyz',
                    'b2' : self.dirname+'benzene_exp.xyz',
                    'bc' : self.dirname+'benzene_Conly.xyz',
                    'ha' : self.dirname+'h_atom.xyz'}

        if self.target[1] == '1':
            osa(structs[self.target[0]])
            try:
                # Catch the screen output of the Fortran call 
                # 
                # magic to capture that output:
                # from http://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python
                #      http://websrv.cs.umt.edu/isis/index.php/F2py_example
                output_file = 'pore.out'
                if os.path.exists(output_file):
                    os.remove(output_file)
                # open outputfile
                outfile = os.open(output_file, os.O_RDWR|os.O_CREAT)
                # save the current file descriptor
                save = os.dup(1)
                # put outfile on 1
                os.dup2(outfile, 1)
                # end magic

                # Fortran call
                osa(structs[self.target[0]])
                # restore the standard output file descriptor
                os.dup2(save, 1)
                # close the output file
                os.close(outfile)
                f = open(output_file,'r')
                output = f.read()
                f.close()
            except: output = 'You have not provided the path to pore.so!'

        if self.target[1] == '2':
            self.probe_r = self.e4.get()      
            self.grid_a = self.e5.get()
            self.grid_b = self.e6.get()
            self.grid_c = self.e7.get()
            try:
                # Catch the screen output of the Fortran call 
                # 
                # magic to capture that output:
                # from http://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python
                #      http://websrv.cs.umt.edu/isis/index.php/F2py_example
                output_file = 'pore.out'
                if os.path.exists(output_file):
                    os.remove(output_file)
                # open outputfile
                outfile = os.open(output_file, os.O_RDWR|os.O_CREAT)
                # save the current file descriptor
                save = os.dup(1)
                # put outfile on 1
                os.dup2(outfile, 1)
                # end magic

                # Fortran call
                if self.exe_pore.get() == 'subgrid':
                    gpa_FullGrid(structs[self.target[0]],self.probe_r,self.grid_a,self.grid_b,self.grid_c)
                # restore the standard output file descriptor
                os.dup2(save, 1)
                # close the output file
                os.close(outfile)
                f = open(output_file,'r')
                output = f.read()
                f.close()
            except: output = 'You have not provided the path to pore.so!'
        self.text.delete("1.0", "end")
        self.text.insert(END,output)
        self.text.pack()
        self.text.update()
        self.newwin.update()
        self.newwin.deiconify()

    def get_psd(self):
        self.newwin.destroy()
        self.init_newwin()
        s = ''' PyPorE: pore size distribution (PSD) \n 
                \t structure \n 
                \t\t\t ud -- user defined; use "open" to select the structure \n 
                \t\t\t do -- predefined DUT-8(Ni) open, see others as well \n 
                \t Starting points -- number of starting points \n 
                \t Monte-Carlo cycles -- number of Monte-Carlo cycles \n
                \t Authors: \n
                \t\t\t Kai Trepte (Fortran, core routines) \n
                \t\t\t Sebastian Schwalbe (Python, GUI)
            '''
        self.text.insert(END,s)
        self.text.pack() 
        # Search input fields 
        # What property to search for: when, where or ...
        Label(self.controls, text="Select structure").grid(row=1,column=1)
        Label(self.controls, text="Starting points").grid(row=2,column=1)
        Label(self.controls, text="Monte-Carlo cycles").grid(row=3,column=1)
        Label(self.controls, text="e.g. do").grid(row=1,column=3)
        Label(self.controls, text="e.g. 200").grid(row=2,column=3)
        Label(self.controls, text="e.g. 2000").grid(row=3,column=3)

        # Entry fields 
        e1 = StringVar(self.newwin)
        self.e1 = e1
        self.e1.set("ud") # default value 
        w1= OptionMenu(self.controls, self.e1, 'ud','do','vo','dc','vc','u6','u7','m5','ir','m2','h1','ho','c6','be','b2','bc','ha',command=self.refresh_psd)

        self.w1 = w1

        e2 = Entry(self.controls)
        e2.insert(END,'200') # default value 
    
        e3 = Entry(self.controls)
        e3.insert(END,'2000') # default value 
        
        self.e1 = e1
        self.e2 = e2
        self.e3 = e3 

        self.w1.grid(row=1, column=2)
        self.e2.grid(row=2, column=2)
        self.e3.grid(row=3, column=2)
        # Action buttons
        if self.e1.get() == 'ud':
            b3 = Button(self.controls, text='Open', command=self.user_input)
            b3.grid(row=4,column=1)
            self.b3 =b3
        b1 = Button(self.controls, text='OK', command=self.cmd_psd)
        b1.grid(row=4,column=2)
        self.b1 = b1
        self.text.update()
        #self.text.deiconify()
        self.newwin.update()
        self.newwin.deiconify()

    
    def cmd_psd(self):
        self.get_target()
        # KT: get correct structure inputs
        structs  = {'ud' : 'pypore.xyz',
                    'do' : self.dirname+'dut_8_open.xyz',
                    'vo' : self.dirname+'dut_8_open_vcrelax.xyz',
                    'dc' : self.dirname+'dut_8_closed.xyz',
                    'vc' : self.dirname+'dut_8_closed_vcrelax.xyz',
                    'u6' : self.dirname+'uio66.xyz',
                    'u7' : self.dirname+'uio67.xyz',
                    'u8' : self.dirname+'uio68.xyz',
                    'm5' : self.dirname+'mof5.xyz',
                    'ir' : self.dirname+'irmof10.xyz',
                    'm2' : self.dirname+'mof210.xyz',
                    'h1' : self.dirname+'hkust1.xyz',
                    'ho' : self.dirname+'hkust1_with_O.xyz',
                    'c6' : self.dirname+'c60_MOF.xyz',
                    'be' : self.dirname+'benzene.xyz',
                    'b2' : self.dirname+'benzene_exp.xyz',
                    'bc' : self.dirname+'benzene_Conly.xyz',
                    'ha' : self.dirname+'h_atom.xyz'}
        try:
            # magic to capture that output:
            # from http://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python
            #      http://websrv.cs.umt.edu/isis/index.php/F2py_example
            output_file = 'psd.out'
            if os.path.exists(output_file):
                os.remove(output_file)
            # open outputfile
            outfile = os.open(output_file, os.O_RDWR|os.O_CREAT)
            # save the current file descriptor
            save = os.dup(1)
            # put outfile on 1
            os.dup2(outfile, 1)
            # end magic

            # Fortran call
            #pore(self.target[0],self.target[1],self.probe_r,self.grid_a,self.grid_b,self.grid_c)
            get_PSD(structs[self.target[0]],self.e2.get(),self.e3.get())
            #

            # restore the standard output file descriptor
            os.dup2(save, 1)
            # close the output file
            os.close(outfile)
            f = open(output_file,'r')
            output = f.read()
            f.close()
        except: output = 'You have not provided the path to get_psd.so!'

        self.text.delete("1.0", "end")
        self.text.insert(END,output)
        self.text.pack()
        #self.display.pack()
        #self.text.grid(row=1, columnspan=4, rowspan=4,padx=5, sticky=E+W+S+N)
        #self.display.grid(row=1, columnspan=4, rowspan=4,padx=5, sticky=E+W+S+N)
        self.text.update()
        self.newwin.update()
        self.newwin.deiconify()


    def init_newwin(self):
        # The search check boxes
        check_box_list = []
        self.check_box_list = check_box_list
        
        # Output/Result window
        newwin = Toplevel(self.master)
        
        area = Canvas(newwin)
        controls = Frame(newwin)

        area.pack(side="left", fill="both", expand=True)
        controls.pack(side="right", fill="both", expand=False)
        self.area = area
        self.controls = controls

        scroll = Scrollbar(newwin)
        # Label of the new window
        #display = Label(newwin, text='Info')
        # scrollable text in new window
        text = Text(self.area, height=40, width=120,yscrollcommand=scroll.set)
        self.newwin = newwin
        #self.display = display
        self.text = text
        self.newwin.withdraw()
        self.b = None

    def get_target(self):
        # get structure and method 
        target = [self.e1.get(),self.e2.get()]
        self.target = target 

    def ase2pore(self):
        # convert ase structural information into pore input 
        struct = read(self.name)
        pos = struct.get_positions()
        sym = struct.get_chemical_symbols()
        cell = struct.get_cell()[:]
        self.cell = cell
        f = open('pypore.xyz','w') 
        f.write('%i \n' %(len(pos))) 
        f.write('%0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f \n' %(cell[0,0],cell[0,1],cell[0,2],cell[1,0],cell[1,1],cell[1,2],cell[2,0],cell[2,1],cell[2,2]))
        for s in range(len(sym)):
            f.write('%s %0.9f %0.9f %0.9f\n' %(sym[s],pos[s][0],pos[s][1],pos[s][2]))
        f.close() 
    
    def user_input(self):
        name= askopenfilename(initialdir=self.dirname,filetypes=[("xyz pypore files",".xyz"),("cif files",".cif"),("POSCAR","POSCAR")])
        try:
            datatype = name.split('.')[-1]
        except: datatype = 'error'
        # material structural formats using ase
        if name =='POSCAR' or datatype =='cif':
            self.name = name 
            self.ase2pore()
        # KT: If xyz -> predefined porE xyz file -> write pypore file here
        if datatype == 'xyz':
            org_file = open(name,'r')
            new_file = open('pypore.xyz','w')
            lines_org = org_file.readlines()
            org_file.close()
            for s in range(len(lines_org)):
                new_file.write(lines_org[s])
            new_file.close()



    def open_file(self):
        # select a file name in the selected folder :
        name= askopenfilename(initialdir=self.dirname)
        try:
            datatype = name.split('.')[-1]
        except: datatype = 'error'
        # text files
        if datatype =='txt' or datatype =='dat' or datatype =='out' or name.find('README') != -1:
            with open(name,'r') as UseFile:
                s = UseFile.read()
                self.text.delete("1.0", "end")
                self.text.insert(END,s)
                self.text.pack()
                #self.display.pack()
                #self.text.grid(column=0,row=0)
                #self.columnconfigure(1, weight=1)
                #self.columnconfigure(3, pad=7)
                #self.rowconfigure(3, weight=1)
                #self.rowconfigure(5, pad=7)

                #self.text.grid(row=0, columnspan=2, rowspan=4,padx=5, sticky=E+W+S+N)
                #self.display.grid(row=0, columnspan=2, rowspan=4,padx=5, sticky=E+W+S+N)

                self.text.update()
                self.newwin.update()
                self.newwin.deiconify()
        # material structural formats using ase
        if datatype =='xyz' or datatype =='cif':
            struct = read(name)
            view(struct)
    
    def go(self):
        # find and structure and watch with ase 
        # e.g. analyze folder and viewing txt files 
        t = self.allstates()
        self.t = t
        self.init_newwin()
        self.text.delete("1.0", "end")
        if self.b != None:
            self.b.pack_forget()
        self.b = Button(self.newwin,text='File Open', command=self.open_file)
        self.b.pack(fill=X)
        for i in range(len(t)):
            if self.t[i] == 1:
                os.chdir(self.r[self.keys[i]]['where'])
                ls = os.listdir('./')
                ls_str = '\n'.join(map(str, ls))
                self.dirname = self.r[self.keys[i]]['where']
                self.newwin.update()
                self.newwin.deiconify()

if __name__ == '__main__': 
    pe = PyPore() 
