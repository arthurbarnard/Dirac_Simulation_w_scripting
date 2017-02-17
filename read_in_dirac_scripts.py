import numpy as np
import os
import matplotlib.pyplot as plt
import Tkinter as tk
import tkFileDialog as filedialog
from dirac_run_timedomain import Run_Dirac
import multiprocessing
import sys

def main():
	# brings up a GUI window to pick the script file to run
	nargin=len(sys.argv)-1
	if nargin==0:
		root = tk.Tk()
		root.withdraw()
		fname = filedialog.askopenfilename()
	else:
		fname = sys.argv[1]
	
	#read in the file
	with open(fname) as f:
		content = f.readlines()

	#lists the search parameter text for the script file. In the 
	#script file a single tab should immediately preceed the read-in 
	#values (you can pad with extra spaces otherwise) Order in Input_texts
	#must match order in Input_var_names
	Input_texts=["Input Matrices",
				"wave vector mag",
				"wave drive style",
				"wave drive param",
				"plot results",
				"save timesteps",
				"Number of timesteps",
				"timestep size",
				"B-field multiplier",
				"Potental multiplier",
				"sweep parameter",
				"number of cores",
				"output file path",
				"output filename"]
	
	#types of input associated with the variables in Input_var_names 
	#below. There are "string", "int", "float", "floats", "Y_N", and "list"
	#most are self-explanatory, "floats" reads in a comma separated array 
	#of numbers. "list" will take in string arguements and map them int
	
	Input_type =["string",
				"float",
				"list",
				"floats",
				"Y_N",
				"Y_N",
				"int",
				"float",
				"float",
				"float",
				"floats",
				"int",
				"string",
				"string"]
				
	#variable names that will take on values associated with the above
	#lines read from the script file
	Input_var_names=["fname_in",
				"p0_in",
				"inject_style",
				"inject_params",
				"plot_results",
				"save_time",
				"N_t",
				"dt",
				"A_mult",
				"V_mult",
				"sweep_params",
				"N_core",
				"path_out",
				"fname_out"]
	
	#instantiates input variables
	for i, input_var in enumerate(Input_var_names):
		exec(input_var+" = []")
		
	#reads through script file line by line
	for i, line in enumerate(content):
		
		# looks for tab in a line as a signal of data input
		if len(line.split("\t",1))>1:
			
			#loops through all potential input texts (this way order in the script doesn't matter)
			for k, input_text in enumerate(Input_texts):
				
				# identifies if the line has a paremeter to sweep
				if is_inline(line,"sweep1"):
			
					sweep_line=line
					sweep_lnum=i

				# reads parameter and assigns value to variable
				elif is_inline(line,input_text):
					
					if Input_type[k]=="string":
						exec(Input_var_names[k]+" = line.split(\"\\t\",1)[1].split(\"\\n\")[0]\n")
					if Input_type[k]=="float":
						exec(Input_var_names[k]+" = float(line.split(\"\\t\",1)[1].split(\"\\n\")[0])")
					if Input_type[k]=="int":
						exec(Input_var_names[k]+" = int(line.split(\"\\t\",1)[1].split(\"\\n\")[0])")
					if Input_type[k]=="floats":
						exec(Input_var_names[k]+" = [float(n) for n in (line.split(\"\\t\",1)[1].split(\"\\n\")[0].split(\",\"))]")
					if Input_type[k]=="Y_N":
						exec(Input_var_names[k]+" = line.split(\"\\t\",1)[1].find(\"Y\")>=0")
					if Input_type[k]=="list":
					
						if is_inline(line,"wave drive style"):
							if line.split("\t",1)[1].find("gauss")>=0:
								inject_style=1
								print "test"
							elif line.split("\t",1)[1].find("plane wave")>=0:
								inject_style=2
								print "test"
							elif line.split("\t",1)[1].find("pulse")>=0:
								inject_style=3
								print "test"
									
	#makes the generalized sweep vector
	sweep_num=np.linspace(sweep_params[0],sweep_params[1],sweep_params[2])

	
	
	jobs=[]
	n1=0 #n is the 
	
	#sweeps through the generalized sweep parameter
	for i,sweep_N in enumerate(sweep_num):
		if 'sweep_line' in locals(): 
			line=sweep_line.replace("sweep1",str(sweep_N))
			print line, Input_var_names[sweep_lnum]
			if Input_type[sweep_lnum]=="float":
				exec(Input_var_names[sweep_lnum]+" = float(line.split(\"\\t\",1)[1].split(\"\\n\")[0])")
			if Input_type[sweep_lnum]=="int":
				exec(Input_var_names[k]+" = int(line.split(\"\\t\",1)[1].split(\"\\n\")[0])")		
			if Input_type[sweep_lnum]=="floats":
				exec(Input_var_names[sweep_lnum]+" = [float(n) for n in (line.split(\"\\t\",1)[1].split(\"\\n\")[0].split(\",\"))]")
								
		print path_out+fname_out+("%03d" % i)+".mat"
		
		# checks if the output file already exists
		if not(os.path.isfile(path_out+fname_out+("%03d" % i)+".mat")):
			
			#registers N_core number of parallel simulations
			if N_core>1:
				p = multiprocessing.Process(target=Run_Dirac, args=(fname_in,p0_in,inject_style,inject_params,plot_results,save_time,N_t, dt, A_mult,V_mult, path_out,fname_out+("%03d" % i)))
				jobs.append(p)
				p.start()
				n1+=1
				if n1>=N_core:
					p.join()
					print n1
					n1=0
			else:
				Run_Dirac(fname_in,p0_in,inject_style,inject_params,plot_results,save_time,N_t, dt, A_mult,V_mult, path_out,fname_out+("%03d" % i))

		print fname_in
		print p0_in
		print inject_style
		print inject_params
		print plot_results
		print N_t
		print dt
		print A_mult
		print V_mult
		print sweep_params
		
		
#checks if "string_name" is in the string name "line"
def is_inline(line,string_name):
		return line.find(string_name)>=0
		
if __name__ == '__main__':
	main()