import matplotlib.pyplot as plt
import numpy as np
import time
from matplotlib.path import Path
pi=np.pi

# "dirac_sheet.py" is an object class that defines a 2D discretized sheet for calculating solutions to the 
# time-domain solutions of the Dirac equation. It implements a staggered time and staggered space discretization
# approach outlined in Journal of Computational Physics 289(2015) 169-180. The discretized space is constructed of 
# checkerboard "u" and "v" lattice. The checkerboard is formed by explicitly defining a "u1","u2","v1", and "v2"
# grid that forms a unit cell on a square grid. Within a unit cell, "u1" is the upper-right, "u2" is the lower 
# left, "v1" is the upper left, and #v2 is the lower right.
#
# The default sheet is a square of NgridxNgrid unit cells with periodic boundary conditions. Further boundary 
# conditions can be set by three matrices: No_prop_mat, Absorb_mat, and Drive_mat. No_prop_mat defines regions
# over which no propagation occurs (similar to defining the edge of the 2D material). Absorb_mat defines regions
# where the wave is "lossy". This is primarily intended to introduce absorptive boundary conditions. Drive_mat
# defines a region from which the electron wave is sourced.
#
# The properties of the sourced electrons are defined by the set_p method. p=2*pi/lambda is the injected electron
# wave-vector magnitude, and theta is the wave-vector direction. The code computes plane-wave injection with
# the expectation that the user can sum several solutions to localize the injection.
# there is also a gaussian_beam method wherein you set the angle, FWHM, and central location of the beam. 
# 
# For boundary conditions, the u-lattice is set to zero in the "No_prop" region, and the v-lattice is free. The drive
# wave is implemented only on the v-lattice. In this sense, the "u" lattice is the master, and the "v" the slave. 
class dirac_sheet:

	def __init__(self,m, Ngrid, D_t, D_x, X_offset,Y_offset):
	
		self.t=0.0 			   		#time
		self.Drive_coupling=0.3 	#coupling strength of Drive plane-wave
		self.Abs_rate=0.99			#exponential decay rate in absorptive regions
		self.m=m					#effective mass; in graphene=0
		
		if not isinstance(Ngrid,tuple):
			self.Ngrid_x=Ngrid		#number of sub-lattice points: total (2Ngrid)x(2Ngrid) points
			self.Ngrid_y=Ngrid	
		else:
			self.Ngrid_x=Ngrid[0]	
			self.Ngrid_y=Ngrid[1]	
		self.D_t=D_t				#discrete time step
		self.D_x=D_x				#spatial discretization
		self.D_y=D_x
		self.X_offset=X_offset	
		self.Y_offset=Y_offset
		self.p0=2*pi/100.0			#Drive wave-vector magnitude
		self.theta=0				#Drive wave-vector direction
		self.massorV=False			#encodes whether m or V is nonzero
		self.Bfiled=False
		self.px=self.p0*np.cos(self.theta)
		self.py=self.p0*np.sin(self.theta)
		self.y, self.x = np.mgrid[slice((0-round(self.Ngrid_y/2)),(self.Ngrid_y-round(self.Ngrid_y/2)),1),
								slice((0-round(self.Ngrid_x/2)),(self.Ngrid_x-round(self.Ngrid_x/2)),1)]
								
		#define the X,Y coordinate matrix for each sublattice point
		self.x, self.y = self.x*self.D_x+self.X_offset, self.y*self.D_y+self.Y_offset
		self.xhalf, self.yhalf = self.x-self.D_x/2.0, self.y-self.D_y/2.0
		
		self.No_prop_mat=np.zeros((self.x.shape[0]*2,self.x.shape[1]*2))
		self.Absorb_mat=np.zeros((self.x.shape[0]*2,self.x.shape[1]*2))
		self.Drive_mat=np.zeros((self.x.shape[0]*2,self.x.shape[1]*2))
		
		#legacy, should be equivalent to theta
		self.phi_p=np.arctan2(self.py,self.px)
		
		#out and out2 are temporary placeholders used in the wave-equation calculation
		#this enables multiple manipulations of the wave matrix without overwriting
		self.out=np.zeros(self.x.shape,dtype=np.complex64)
		self.out2=np.zeros(self.x.shape,dtype=np.complex64)
		#defining the u1,u2,v1,v2 sublattices
		self.u1=np.zeros(self.x.shape,dtype=np.complex64)
		self.u2=np.zeros(self.x.shape,dtype=np.complex64)
		self.v1=np.zeros(self.x.shape,dtype=np.complex64)
		self.v2=np.zeros(self.x.shape,dtype=np.complex64)
		
		#Drive plane wave
		self.u10=np.exp(1j*(self.px*self.x+self.py*self.y))
		self.u20=np.exp(1j*(self.px*self.xhalf+self.py*self.yhalf))
		self.v10=np.exp(1j*(self.px*self.xhalf+self.py*self.y))
		self.v20=np.exp(1j*(self.px*self.x+self.py*self.yhalf))
		
		self.N_zero=self.x[1,np.abs(self.x[1,:])==np.min(np.abs(self.x[1,:]))]
		
		#logical arrays for boundary conditions
		self.rDu1=self.Drive_mat[1::2,1::2]!=0.0
		self.rDu2=self.Drive_mat[0::2,0::2]!=0.0
		self.rDv1=self.Drive_mat[1::2,0::2]!=0.0
		self.rDv2=self.Drive_mat[0::2,1::2]!=0.0
				
		self.rAu1=self.Absorb_mat[1::2,1::2]!=0.0
		self.rAu2=self.Absorb_mat[0::2,0::2]!=0.0
		self.rAv1=self.Absorb_mat[1::2,0::2]!=0.0
		self.rAv2=self.Absorb_mat[0::2,1::2]!=0.0
		
		self.rNPu1=self.No_prop_mat[1::2,1::2]!=0.0
		self.rNPu2=self.No_prop_mat[0::2,0::2]!=0.0
		self.rNPv1=self.No_prop_mat[1::2,0::2]!=0.0
		self.rNPv2=self.No_prop_mat[0::2,1::2]!=0.0
		
		
		#V1-V4 are electrostatic scalar potentials for the v1,v2,u1,and u2
		#lattices respectively
		self.V1=np.zeros(self.x.shape,dtype=np.complex64)
		self.V2=np.zeros(self.x.shape,dtype=np.complex64)
		self.V3=np.zeros(self.x.shape,dtype=np.complex64)
		self.V4=np.zeros(self.x.shape,dtype=np.complex64)
		
		#iAxAy1-iAxAy4 are electrostatic vector potentials for the v1,v2,u1,and u2
		#lattices respectively
		self.iAxAy1=np.zeros(self.x.shape,dtype=np.complex64)
		self.iAxAy2=np.zeros(self.x.shape,dtype=np.complex64)
		self.iAxAy3=np.zeros(self.x.shape,dtype=np.complex64)
		self.iAxAy4=np.zeros(self.x.shape,dtype=np.complex64)
		
		#precomputes sum/difference of mass and potential to reduce computation overhead
		self.mminusV1=self.m-self.V1
		self.mminusV2=self.m-self.V2
		self.mplusV3=self.m+self.V3
		self.mplusV4=self.m+self.V4
		
	def set_p(self,p0,theta):
		#sets the injection wave-vector properties, p is amplitude, theta is direction
		self.p0=p0
		self.theta=theta
		self.px=self.p0*np.cos(self.theta)
		self.py=self.p0*np.sin(self.theta)
		self.phi_p=np.arctan2(self.py,self.px)
		self.u10=np.exp(1j*(self.px*self.x+self.py*self.y))
		self.u20=np.exp(1j*(self.px*self.xhalf+self.py*self.yhalf))
		self.v10=np.exp(1j*(self.px*self.xhalf+self.py*self.y+self.phi_p))
		self.v20=np.exp(1j*(self.px*self.x+self.py*self.yhalf+self.phi_p))
	
	def set_Bfield(self,B):
		self.Bfield=(B!=0)
		self.iAxAy1=B*self.xhalf
		self.iAxAy2=B*self.x
		self.iAxAy3=B*self.x
		self.iAxAy4=B*self.xhalf

	def set_gaussian_beam(self,p0,theta0,delta_theta,x0,y0):
		# sets a gaussian_beam as the injection parameter. theta0 is the central angle,
		# delta_theta is the FWHM of psi squared. x0 and y0 are the location of the beam
		# focal point
		self.p0=p0
		N_integrate=300
		delta_theta/=(np.sqrt(4*np.log(2)))
		self.u10=np.zeros(self.x.shape,dtype=np.complex64)
		self.u20=np.zeros(self.x.shape,dtype=np.complex64)
		self.v10=np.zeros(self.x.shape,dtype=np.complex64)
		self.v20=np.zeros(self.x.shape,dtype=np.complex64)
		
		for theta in np.linspace(theta0-np.pi/2,theta0+np.pi/2,N_integrate):
			self.theta=theta
			self.px=self.p0*np.cos(self.theta)
			self.py=self.p0*np.sin(self.theta)
			self.phi_p=np.arctan2(self.py,self.px)
			self.u10[self.rDu1]+=np.pi*self.p0/N_integrate*np.exp(-(theta-theta0)**2/2/delta_theta**2)*np.exp(1j*(self.px*(self.x[self.rDu1]-x0)+self.py*(self.y[self.rDu1]-y0)))
			self.u20[self.rDu2]+=np.pi*self.p0/N_integrate*np.exp(-(theta-theta0)**2/2/delta_theta**2)*np.exp(1j*(self.px*(self.xhalf[self.rDu2]-x0)+self.py*(self.yhalf[self.rDu2]-y0)))
			self.v10[self.rDv1]+=np.pi*self.p0/N_integrate*np.exp(-(theta-theta0)**2/2/delta_theta**2)*np.exp(1j*(self.px*(self.xhalf[self.rDv1]-x0)+self.py*(self.y[self.rDv1]-y0)+self.phi_p))
			self.v20[self.rDv2]+=np.pi*self.p0/N_integrate*np.exp(-(theta-theta0)**2/2/delta_theta**2)*np.exp(1j*(self.px*(self.x[self.rDv2]-x0)+self.py*(self.yhalf[self.rDv2]-y0)+self.phi_p))
	
	def get_pos_mat(self):
		#outputs the stiched X,Y coordinates; can be useful for defining boundary conditions
		Xtot_out=np.zeros((self.x.shape[0]*2,self.x.shape[1]*2))
		Ytot_out=np.zeros((self.x.shape[0]*2,self.x.shape[1]*2))
		Xtot_out[::2,::2]=self.xhalf
		Xtot_out[1::2,::2]=self.xhalf
		Xtot_out[::2,1::2]=self.x
		Xtot_out[1::2,1::2]=self.x
		Ytot_out[::2,::2]=self.yhalf
		Ytot_out[1::2,::2]=self.y
		Ytot_out[::2,1::2]=self.yhalf
		Ytot_out[1::2,1::2]=self.y
		return Xtot_out, Ytot_out
		
	def set_No_prop_mat(self,No_prop_mat):
		#imports No_prop_mat and creates logical matrices to define propagation regions	
		#and their edges
		self.No_prop_mat=No_prop_mat!=0.0
		self.Nprop_up=~self.No_prop_mat&np.roll(self.No_prop_mat,-1,axis=0)
		self.Nprop_down=~self.No_prop_mat&np.roll(self.No_prop_mat,1,axis=0)
		self.Nprop_right=~self.No_prop_mat&np.roll(self.No_prop_mat,-1,axis=1)
		self.Nprop_left=~self.No_prop_mat&np.roll(self.No_prop_mat,1,axis=1)
		
		
		temp=np.zeros(self.No_prop_mat.shape,dtype=np.uint8)
		for i in range(-1,2,2):
			for j in range(-1,2,2):
				temp+=~(np.roll(np.roll(self.No_prop_mat,i,axis=0),j,axis=1))
		
		#Border of the no-prop region
		self.Nprop_edge=(self.No_prop_mat)*temp
		
		self.rNPu1=self.No_prop_mat[1::2,1::2]
		self.rNPu2=self.No_prop_mat[0::2,0::2]
		self.rNPv1=self.No_prop_mat[1::2,0::2]
		self.rNPv2=self.No_prop_mat[0::2,1::2]
		
		self.ru1_edge=self.Nprop_edge[1::2,1::2]
		self.ru2_edge=self.Nprop_edge[0::2,0::2]
		self.rv1_edge=self.Nprop_edge[1::2,0::2]
		self.rv2_edge=self.Nprop_edge[0::2,1::2]
		
		self.ru1_edge0=self.ru1_edge>0
		self.ru2_edge0=self.ru2_edge>0
		self.rv1_edge0=self.rv1_edge>0
		self.rv2_edge0=self.rv2_edge>0
		
		#Encodes pixel-by-pixel the existence or absence of nearest neighbors
		self.ru1_edge=np.zeros((self.rNPu1.shape[0],(self.rNPu1.shape[1]),4),dtype=np.bool)
		self.ru1_edge[:,:,0]=self.Nprop_up[1::2,1::2]
		self.ru1_edge[:,:,1]=self.Nprop_right[1::2,1::2]
		self.ru1_edge[:,:,2]=self.Nprop_down[1::2,1::2]
		self.ru1_edge[:,:,3]=self.Nprop_left[1::2,1::2]
		self.ru1_edgeX=self.ru1_edge[:,:,1]|self.ru1_edge[:,:,3]
		self.ru1_edgeY=self.ru1_edge[:,:,0]|self.ru1_edge[:,:,2]
				
		self.ru2_edge=np.zeros((self.rNPu2.shape[0],(self.rNPu2.shape[1]),4),dtype=np.bool)
		self.ru2_edge[:,:,0]=self.Nprop_up[0::2,0::2]
		self.ru2_edge[:,:,1]=self.Nprop_right[0::2,0::2]
		self.ru2_edge[:,:,2]=self.Nprop_down[0::2,0::2]
		self.ru2_edge[:,:,3]=self.Nprop_left[0::2,0::2]
		self.ru2_edgeX=self.ru2_edge[:,:,1]|self.ru2_edge[:,:,3]
		self.ru2_edgeY=self.ru2_edge[:,:,0]|self.ru2_edge[:,:,2]

		self.rv1_edge=np.zeros((self.rNPv1.shape[0],(self.rNPv1.shape[1]),4),dtype=np.bool)
		self.rv1_edge[:,:,0]=self.Nprop_up[1::2,0::2]
		self.rv1_edge[:,:,1]=self.Nprop_right[1::2,0::2]
		self.rv1_edge[:,:,2]=self.Nprop_down[1::2,0::2]
		self.rv1_edge[:,:,3]=self.Nprop_left[1::2,0::2]		
		self.rv1_edgeX=self.rv1_edge[:,:,1]|self.rv1_edge[:,:,3]
		self.rv1_edgeY=self.rv1_edge[:,:,0]|self.rv1_edge[:,:,2]

	
		self.rv2_edge=np.zeros((self.rNPv2.shape[0],(self.rNPv2.shape[1]),4),dtype=np.bool)
		self.rv2_edge[:,:,0]=self.Nprop_up[0::2,1::2]
		self.rv2_edge[:,:,1]=self.Nprop_right[0::2,1::2]
		self.rv2_edge[:,:,2]=self.Nprop_down[0::2,1::2]
		self.rv2_edge[:,:,3]=self.Nprop_left[0::2,1::2]
		self.rv2_edgeX=self.rv2_edge[:,:,1]|self.rv2_edge[:,:,3]
		self.rv2_edgeY=self.rv2_edge[:,:,0]|self.rv2_edge[:,:,2]

	
	def set_Afield_mat(self,Ax,Ay):
		if np.sum(np.abs(Ax)+np.abs(Ay))==0:
			self.Bfiled=False
		else:
			self.Bfield=True
			
		self.iAxAy1=1i*Ax[1::2,0::2]-Ay[1::2,0::2]
		self.iAxAy2=1i*Ax[0::2,1::2]-Ay[0::2,1::2]
		self.iAxAy3=1i*Ax[1::2,1::2]+Ay[1::2,1::2]
		self.iAxAy4=1i*Ax[0::2,0::2]+Ay[0::2,0::2]
		
		
	def set_Absorb_mat(self,Absorb_mat):
		#sets the Absorb_mat, this is is to simulate absorptive sidewall etc.
		self.Absorb_mat=Absorb_mat
		self.rAu1=self.Absorb_mat[1::2,1::2]!=0.0
		self.rAu2=self.Absorb_mat[0::2,0::2]!=0.0
		self.rAv1=self.Absorb_mat[1::2,0::2]!=0.0
		self.rAv2=self.Absorb_mat[0::2,1::2]!=0.0
		self.rAu1_flat=np.where(self.rAu1)
		self.rAu2_flat=np.where(self.rAu2)
		self.rAv1_flat=np.where(self.rAv1)
		self.rAv2_flat=np.where(self.rAv2)
	
	def set_Drive_mat(self,Drive_mat):
		#sets Drive_mat. this defines the region where the drive wave is coupled-in
		self.Drive_mat=Drive_mat
		self.rDu1=self.Drive_mat[1::2,1::2]!=0.0
		self.rDu2=self.Drive_mat[0::2,0::2]!=0.0
		self.rDv1=self.Drive_mat[1::2,0::2]!=0.0
		self.rDv2=self.Drive_mat[0::2,1::2]!=0.0
		self.rDu1_flat=np.where(self.rDu1)
		self.rDu2_flat=np.where(self.rDu2)
		self.rDv1_flat=np.where(self.rDv1)
		self.rDv2_flat=np.where(self.rDv2)
		
	def set_V(self,V_in):
		#sets the potential
		self.V1=V_in[1::2,0::2]
		self.V2=V_in[0::2,1::2]
		self.V3=V_in[1::2,1::2]
		self.V4=V_in[0::2,0::2]
		
		self.mminusV1=self.m-self.V1
		self.mminusV2=self.m-self.V2
		self.mplusV3=self.m+self.V3
		self.mplusV4=self.m+self.V4
		
		if (self.m==0)&(np.sum(V_in)==0):
			self.massorV=False
		else:
			self.massorV=True
			
	def v_step(self):
	
		#steps the v-lattice relative to a fixed u-lattice
		
		self.t+=self.D_t/2.0
		
		#introduce the drive wave to the u-lattice
		self.u1[self.rDu1]+=self.Drive_coupling*self.u10[self.rDu1]*np.exp(-1j*self.p0*self.t)
		self.u2[self.rDu2]+=self.Drive_coupling*self.u20[self.rDu2]*np.exp(-1j*self.p0*self.t)
		
		#Apply absorption
		self.u1[self.rAu1]*=self.Abs_rate
		self.u2[self.rAu2]*=self.Abs_rate
		
		#wave equation
		if self.massorV:
			self.out=(1/(2*1j+self.D_t*self.mminusV1))*((2*1j-self.D_t*self.mminusV1)*self.v1-(2*1j*self.D_t/self.D_x)*(self.u1-np.roll(self.u1,1,axis=1))+(2*self.D_t/self.D_y)*(np.roll(self.u2,-1,axis=0)-self.u2))		
			self.out2=(1/(2*1j+self.D_t*self.mminusV2))*((2*1j-self.D_t*self.mminusV2)*self.v2-(2*1j*self.D_t/self.D_x)*(np.roll(self.u2,-1,axis=1)-self.u2)+(2*self.D_t/self.D_y)*(self.u1-np.roll(self.u1,1,axis=0)))
		elif self.Bfield:
			self.out=(1/(2*1j+self.D_t*self.mminusV1))*((2*1j-self.D_t*self.mminusV1)*self.v1-(2*1j*self.D_t/self.D_x)*(self.u1-np.roll(self.u1,1,axis=1))+(2*self.D_t/self.D_y)*(np.roll(self.u2,-1,axis=0)-self.u2)
				+0.5*1i*D_t*self.iAxAy1.*(self.u1+np.roll(self.u1,-1,axis=1)+self.u2+np.roll(self.u2,1,axis=0)))		
			self.out2=(1/(2*1j+self.D_t*self.mminusV2))*((2*1j-self.D_t*self.mminusV2)*self.v2-(2*1j*self.D_t/self.D_x)*(np.roll(self.u2,-1,axis=1)-self.u2)+(2*self.D_t/self.D_y)*(self.u1-np.roll(self.u1,1,axis=0))
				+0.5*1i*D_t*self.iAxAy2.*(self.u1+np.roll(self.u1,-1,axis=0)+self.u2+np.roll(self.u2,1,axis=1)))		
		else:
			self.out=(1/(2*1j))*((2*1j)*self.v1-(2*1j*self.D_t/self.D_x)*(self.u1-np.roll(self.u1,1,axis=1))+(2*self.D_t/self.D_y)*(np.roll(self.u2,-1,axis=0)-self.u2))		
			self.out2=(1/(2*1j))*((2*1j)*self.v2-(2*1j*self.D_t/self.D_x)*(np.roll(self.u2,-1,axis=1)-self.u2)+(2*self.D_t/self.D_y)*(self.u1-np.roll(self.u1,1,axis=0)))
		
		self.v1=self.out
		self.v2=self.out2
	
	def u_step(self):
		#steps the u-lattice relative to a fixed v lattice
		self.t+=self.D_t/2.0
	
		#Apply absorption
		self.v1[self.rAv1]*=self.Abs_rate
		self.v2[self.rAv2]*=self.Abs_rate

		#wave equation
		if self.massorV:
			self.out=(~self.rNPu1)*(1/(2*1j-self.D_t*self.mplusV3))*((2*1j+self.D_t*self.mplusV3)*self.u1-(2*1j*self.D_t/self.D_x)*(np.roll(self.v1,-1,axis=1)-self.v1)-(2*self.D_t/self.D_y)*(np.roll(self.v2,-1,axis=0)-self.v2))
			self.out2=(~self.rNPu2)*(1/(2*1j-self.D_t*self.mplusV4))*((2*1j+self.D_t*self.mplusV4)*self.u2-(2*1j*self.D_t/self.D_x)*(self.v2-np.roll(self.v2,1,axis=1))-(2*self.D_t/self.D_y)*(self.v1-np.roll(self.v1,1,axis=0)))
		elif self.Bfield:
			self.out=(~self.rNPu1)*(1/(2*1j-self.D_t*self.mplusV3))*((2*1j+self.D_t*self.mplusV3)*self.u1-(2*1j*self.D_t/self.D_x)*(np.roll(self.v1,-1,axis=1)-self.v1)-(2*self.D_t/self.D_y)*(np.roll(self.v2,-1,axis=0)-self.v2)
				+0.5*1i*D_t*self.iAxAy3.*(self.v1+np.roll(self.v1,1,axis=1)+self.v2+np.roll(self.v2,1,axis=0)))		
			self.out2=(~self.rNPu2)*(1/(2*1j-self.D_t*self.mplusV4))*((2*1j+self.D_t*self.mplusV4)*self.u2-(2*1j*self.D_t/self.D_x)*(self.v2-np.roll(self.v2,1,axis=1))-(2*self.D_t/self.D_y)*(self.v1-np.roll(self.v1,1,axis=0)))
				+0.5*1i*D_t*self.iAxAy4.*(self.v1+np.roll(self.v1,-1,axis=0)+self.v2+np.roll(self.v2,-1,axis=1)))		
		else:
			self.out=(~self.rNPu1)*(1/(2*1j))*((2*1j)*self.u1-(2*1j*self.D_t/self.D_x)*(np.roll(self.v1,-1,axis=1)-self.v1)-(2*self.D_t/self.D_y)*(np.roll(self.v2,-1,axis=0)-self.v2))
			self.out2=(~self.rNPu2)*(1/(2*1j))*((2*1j)*self.u2-(2*1j*self.D_t/self.D_x)*(self.v2-np.roll(self.v2,1,axis=1))-(2*self.D_t/self.D_y)*(self.v1-np.roll(self.v1,1,axis=0)))
		
		self.u1=self.out
		self.u2=self.out2
		
		
	def time_step(self):
		#Takes a full time step D_t in two half steps
		self.v_step()
		self.u_step()

