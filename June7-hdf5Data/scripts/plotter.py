import numpy as np 
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re 

# ===================================== User-defined function definitions ===============================


def get_directory(directory):
	'''
	Sets up directory structure. Returns directory path for root, data, scripts, plots, etc.
	Parameters :
	directory - Identifier for the directory label, eg. root, data, scripts, plots, gm_early_data, organic_data, gm_late_data.
	'''
	dir_dict 						= dict()
	dir_dict['scripts_dir'] 		= os.path.dirname(os.path.abspath(__file__))
	dir_dict['root_dir'] 			= os.path.dirname(dir_dict['scripts_dir'])
	dir_dict['data_dir'] 			= os.path.join(dir_dict['root_dir'], 'data')
	dir_dict['gm_early_data_dir']	= os.path.join(dir_dict['data_dir'], 'gm_early')
	dir_dict['organic_data_dir']	= os.path.join(dir_dict['data_dir'], 'organic')
	dir_dict['gm_late_data_dir']	= os.path.join(dir_dict['data_dir'], 'gm_late')
	dir_dict['plots_dir']			= os.path.join(dir_dict['root_dir'], 'plots')
	return dir_dict[directory+'_dir']

		# ---------------------- Getter functions for hdf5 files dataset columns ---------------------

def get_files(directory, mode='r'):
	'''
	Returns hdf5 file handles for all files in a directory.
	TODO : Test for no hdf5 files
	Parameters :
	directory 	- Path to directory for reading hdf5 files
	mode		- Mode of opening the file; defaults to read mode.
	'''
	return [h5py.File(os.path.join(directory,fname),mode) for fname in os.listdir(directory)]

def get_coords_x(fname):
	'''
	Returns x-coordinates from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Coordinates"][:,0]

def get_coords_y(fname):
	'''
	Returns y-coordinates from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Coordinates"][:,1]

def get_coords_z(fname):
	'''
	Returns z-coordinates from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Coordinates"][:,2]

def get_vel_x(fname):
	'''
	Returns velocities along x-axis from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Velocity"][:,0]

def get_vel_y(fname):
	'''
	Returns velocities along y-axis from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Velocity"][:,1]

def get_vel_z(fname):
	'''
	Returns velocities along z-axis from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Velocity"][:,2]

def get_mass(fname):
	'''
	Returns star particle masses from an hdf5 file.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return fname["Mass"]

def get_redshift(fname):
	'''
	Returns redshift value for an hdf5 file. Uses regex to find redshift value from filename.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return float(re.compile(r'z\d+p\d+').findall(str(fname.file))[0].replace('z','').replace('p','.'))

def get_subdf(fname,params_list):
	'''
	Returns a sub-dataframe containing a subset of columns available in a single hdf5 file.
	Parameters	:
	fname 		- Open handle for an hdf5 file.
	param_list 	- List of columns from hdf5 file to include in a dataframe.
	'''
	subdf = pd.DataFrame()
	params_dict		= {							# Dictionary of callable getter functions corresponding to each column
		'coords_x'	:  get_coords_x,
		'coords_y'	:  get_coords_y,
		'coords_z'	:  get_coords_z,
		'vel_x'		:  get_vel_x,
		'vel_y'		:  get_vel_y,
		'vel_z'		:  get_vel_z,
		'mass'		:  get_mass,
		'redshift'	:  get_redshift
		}
	for param in params_list:
		subdf[param] = params_dict[param](fname)
	return subdf

def get_df(fname_list, params_list):
	'''
	Returns a concatenated dataframe (of a single assembly type) of all sub-dataframes constructed using individual hdf5 files.
	Parameters	:
	fname_list	- list of several hdf5 file handles,
	params_list	- list of columns to be extracted from hdf5 file datasets.
	'''
	df = pd.DataFrame()
	if 'redshift' not in params_list:			# Compulsorily add redshift in list of columns.
		params_list.append('redshift')
	df = pd.concat([df]+[get_subdf(fname,params_list) for fname in fname_list], ignore_index=True)
	return df

def prepare_plot():
	'''
	Set seaborn styling for plots.
	'''
	sns.set_context('paper')
	sns.set_style('dark',
		{'xtick.bottom': True, 
		'ytick.left': True, 
		'axes.spines.left': True, 
		'axes.spines.bottom': True})

def get_particle_distribution(df_list,col='redshift'):
	df 		= pd.DataFrame()
	for subdf in df_list:
		count				= subdf.value_counts(col).sort_index(ascending=True).to_frame(name='counts')
		count['assembly']	= subdf.name
		count[col]			= count.index
		df 					= pd.concat([df,count])
	return df

def plot_particle_distribution(df_list, col='redshift',show=True):
	'''
	Plotter function to plot the distribution of star particles with respect to a variable (redshift, by default).
	TODO : Test for other columns passed as parameter col.
	Parameters :
	df_list 	- List of dataframes (assembly modes) for which distribution is to be plotted using separate hues on the same plot.
	col 		- column with respect to which distribution is to be plotted. Redshift is used by default.
	show 		- passed to plot_or_not() function to evaluate whether to show the plot or save it.
	'''
	prepare_plot()
	dist_df		= get_particle_distribution(df_list,col)
	sns.relplot(data=dist_df,x=col,y='counts',hue='assembly',kind='scatter')		# Plots a scatter plot. To plot histogram instead, use displot with kind='hist'. Counts are then computed automatically.
	plot_or_not(show,plot_name='particle_distribution_wrt_'+col)
	return  							# No return value. Plot is either shown or saved, or nothing is done.


def plot_masses(hdf5_file_list):
	prepare_plot()

def plot_or_not(show,plot_name=None):
	'''
	Function to switch between show and save methods for plots.
	Parameters	:	
	show		- Shows an 'interactive' plot if True, else saves the plot if False. If set to None, does nothing.
	plot_name	- If set, plot is saved using plot_name if show is set to False. If show is set to False, yet plot_name is not provided, a random number is generated for filename.
	'''
	if show == True :
		plt.show()
	elif show == False :
		if plot_name == None :
			plot_name = np.random.randint(10000,99999)
		plt.savefig(os.path.join(get_directory(plots),plot_name,'.png'),
			dpi=480,bbox_inches='tight')
	elif show == None :
		pass 

# ============================================ Main program =============================================


if __name__ == '__main__' :
	
	# ------- Get hdf5 file handles list for each type of assembly.

	gm_early_files		= get_files(get_directory('gm_early_data'))
	organic_files		= get_files(get_directory('organic_data'))
	gm_late_files		= get_files(get_directory('gm_late_data'))

	
	# ------- Available fields for cols :
	# ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['mass', 'coords_x']

	# ------- Get concatenated dataframe containing data for each type of assembly mode.

	gm_early_df 		= get_df(gm_early_files,cols)
	organic_df 			= get_df(organic_files,cols)
	gm_late_df 			= get_df(gm_late_files,cols)

	# ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	gm_early_df.name  	= 'GM-Early'
	organic_df.name  	= 'Organic'
	gm_late_df.name  	= 'GM-Late'

	# ------- Plot the particle number distribution for all assembly modes.

	plot_particle_distribution([gm_early_df,organic_df,gm_late_df],show=True)
	# plt.show()
