import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
import os
from matplotlib.path import Path
from dirac_sheet import dirac_sheet
import scipy.io
import Tkinter
from my_two_clim import my_two_clim

# Run_Dirac performs a time-domain simulation of the continuum dirac equation taking several input parameters
# define the nature of the simulation. "initfile" is a file path for a data file containing many of the relevant
# input variables. It contains details of the discretized grid (Nx,Ny,dX,Xoffset,Yoffset) and simulation geometry
# (V,AbsMat,DriveMat,AxMat). It can inject waves with momentum "p0" with varying "drivestyle". If "drivestyle" is
# 1, then it injects a gaussian beam with "drive_param" being a list of four inputs (injection angle, FWHM angle,
# X0,Y0) where X0 and Y0 define the point of narrowest beam-width. If "drivestyle" is 2, then it injects a plane
# wave with "drive_param" being a single input (injection angle). "plot_results" takes a logical input; if True
# a colormap gui is built to show the simulation as it progresses. N_timesteps and dt set the timing parameters.
# Bmult and Vmult are floats that scale the input A or V matrices. #outpath and outfilename determine the place
# where the output is saved.
def Run_Dirac(initfile,p0,drivestyle,drive_param,plot_results,save_time,N_timesteps, dt,Bmult, Vmult , outpath,outfilename):

	#save_time=True
	#checks if it's a .mat file. later will also accept npz files.
	if initfile.split(".",1)[1]=="mat":
	
		#reads in relevant available variables
		mat = scipy.io.loadmat(initfile)
		if 'dX' in mat: dx=mat['dX']
		if 'Nx' in mat: Nx=mat['Nx']
		if 'Ny' in mat: Ny=mat['Ny']
		if 'Xoffset' in mat: Xoffset=mat['Xoffset']
		if 'Yoffset' in mat: Yoffset=mat['Yoffset']
		
		if 'V' in mat: V=mat['V']
		if 'n0' in mat: n0=mat['n0']
		if 'n00' in mat: n00=mat['n00']
		if 'AbsMat' in mat: AbsMat=mat['AbsMat']
		if 'DriveMat' in mat: DriveMat=mat['DriveMat']
		if 'AxMat' in mat: AxMat=mat['AxMat']
		
		if ('u1' in mat):
			u1=mat['u1']
			u2=mat['u2']
			v1=mat['v1']
			v2=mat['v2']
			if 't' in mat: t=mat['t']
			
			

			
	if not os.path.exists(outpath):
		os.makedirs(outpath)
		
	myDirac=dirac_sheet(0,(Nx,Ny),dt,dx,Xoffset,Yoffset)
	
	if drivestyle==3:
		psi_mean=np.zeros(myDirac.u1.shape)
	#sets properties that were specified in the .mat file. If you don't want 
	if 'V' in mat: myDirac.set_V(V*Vmult)
	if 'n0' in mat:
		n1=-p0**2
		n2=(Vmult-p0)**2*np.sign(Vmult-p0)
		n=n0*(n2-n1)+n1
		if 'n00' in mat: n+=n00
		V=np.sqrt(np.abs(n))*np.sign(n)+p0
		myDirac.set_V(V)
	myDirac.set_Bfield(Bmult)
	if 'V0' in mat: myDirac.set_V(V*Vmult+mat['V0'])
	if 'AbsMat' in mat: myDirac.set_Absorb_mat(AbsMat)
	if 'DriveMat' in mat: myDirac.set_Drive_mat(DriveMat)
	if 'NoPropMat' in mat: myDirac.set_No_prop_mat(NoPropMat)
	
	if drivestyle==1:
		myDirac.set_gaussian_beam(p0,drive_param[0]*np.pi/180,drive_param[1]*np.pi/180,drive_param[2],drive_param[3])
	if drivestyle==3:
		myDirac.set_gaussian_beam(p0,drive_param[0]*np.pi/180,drive_param[1]*np.pi/180,drive_param[2],drive_param[3])
		myDirac.set_gaussian_time(drive_param[4])
		print drive_param[4]
		
	if 't' in locals(): myDirac.t=float(t)
	if 'u1' in locals():
		myDirac.u1=u1
		myDirac.u2=u2
		myDirac.v1=v1
		myDirac.v2=v2

	tic = time.time()
	
	#if instructed, build a "my_two_clim" gui. Presently sets the clim manually. Might want to programatically set for viewing
	#when the drive properties are different.
	if plot_results:
		root = Tkinter.Tk()
		c_scale=.03
		view_gui = my_two_clim(root)
		view_gui.set_left_clim(-c_scale, c_scale)
		view_gui.set_right_clim(0, 2*((c_scale)**2))
		view_gui.set_aratio(Ny/float(Nx))

	#this is the execution loop
	for i in xrange(0,N_timesteps):
		
		#takes a timestep
		myDirac.time_step()
		toc = time.time()
		print toc-tic
		
		tic = time.time()
		
		if drivestyle==3:
			psi_mean+=np.abs(myDirac.u1)**2+np.abs(myDirac.u2)**2+np.abs(myDirac.v1)**2+np.abs(myDirac.v2)**2
		
		# if instructed, update the gui every 10 timesteps. plots a component of psi on the left, and psisqr on the right
		if plot_results: 
			root.update_idletasks()
			if np.mod(i,10)==0:
				root.update()
				view_gui.set_both(np.real(myDirac.u1),np.real(np.abs(myDirac.v1)**2)+np.real(np.abs(myDirac.u1)**2))

		
		if save_time:
			if np.mod(i,10)==0:
				mdict={}
				mdict['u1']=myDirac.u1
				mdict['v1']=myDirac.v1
				scipy.io.savemat(outpath+outfilename+'_time_'+("%03d" % np.int(i/10))+'.mat',mdict)

		#this block packages data and outputs into a .mat file.
		if np.mod(i,1000)==0 or i>=N_timesteps-1:
		
		mdict={}
		mdict['u1']=myDirac.u1
		mdict['u2']=myDirac.u2
		mdict['v1']=myDirac.v1
		mdict['v2']=myDirac.v2
		mdict['t']=myDirac.t
		mdict['initfile']=initfile
		mdict['p0']=p0
		mdict['drivestyle']=drivestyle
		mdict['drive_param']=drive_param
		mdict['dt']=dt
		mdict['Bmult']=Bmult
		mdict['Vmult']=Vmult
		
		if drivestyle==3:	
			mdict['psi_mean']=psi_mean
			
		scipy.io.savemat(outpath+outfilename+'.mat',mdict)
	
	#output npz data files
	#np.savez_compressed(outpath+outfilename+'.npz', u1=myDirac.u1,u2=myDirac.u2,v1=myDirac.v1,v2=myDirac.v2,t=myDirac.t,initfile=initfile,p0=p0,
	#	drivestyle=drivestyle, drive_param=drive_param,dt=dt,Bmult=Bmult,Vmult=Vmult)
		
	#kill the gui
	if plot_results: root.destroy()