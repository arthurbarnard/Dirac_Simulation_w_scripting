import Tkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure	
import time
import numpy as np

# my_two_clim is a plotting gui that presents two parallel colormaps. It updates either or both 
# plot on request. 

class my_two_clim:
	def __init__(self, master):
		# Create a container
		frame = Tkinter.Frame(master)
		
		# create a figure
		fig = Figure(figsize=(12, 6), dpi=96, facecolor='w', edgecolor='k')

		#initialize left plot
		ax = fig.add_subplot(121)
		self.im_left=ax.imshow(np.zeros((400,400)),'RdBu')
		self.im_left.axes.get_xaxis().set_visible(False)
		self.im_left.axes.get_yaxis().set_visible(False)
		self.im_left.set_clim(0,1)
		#initialize right plot
		ax = fig.add_subplot(122)
		self.im_right=ax.imshow(np.zeros((400,400)),'RdBu')
		self.im_right.axes.get_xaxis().set_visible(False)
		self.im_right.axes.get_yaxis().set_visible(False)
		self.im_right.set_clim(0,1)

		#builds the gui
		self.canvas = FigureCanvasTkAgg(fig,master=master)
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
		frame.pack()

		
	#set left,right, or both. loads data_in and updates the plot
	def set_left(self,data_in):
		self.im_left.set_data(data_in)
		self.canvas.draw()
		
	def set_right(self,data_in):
		self.im_right.set_data(data_in)
		self.canvas.draw()
		
	def set_both(self,data_left,data_right):
		self.im_left.set_data(data_left)
		self.im_right.set_data(data_right)
		self.canvas.draw()
	
	def set_aratio(self,aratio):
		self.im_left.axes.set_aspect(aratio)
		self.im_right.axes.set_aspect(aratio)
		
	#sets the color limits (l is the lower limit, h is the upper limit)
	def set_left_clim(self,l,h):
		self.im_left.set_clim(l,h)
		
	def set_right_clim(self,l,h):
		self.im_right.set_clim(l,h)
