import numpy as np 
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import astropy
import astropy.units as u
import os
import re 
import itertools
import time
import inspect
import imageio

# ===================================== User-defined function definitions ===============================

# ----------------------- General-purpose functions -----------------------

def capitalize_first_letter(string):
	'''
	Returns a string with capitalized first letter.
	Parameters	:
	string	- a string whose first words is to be capitalized. This does not return string in camelcase, for which we might need to call this function repeated on split string.
	'''
	return string[0].capitalize()+string[1:]

def add_assembly_column(df_list):
	'''
	Returns a list of dataframes with an extra column added, which signifies mode of assembly of data. This is denoted by name of dataframe.
	Parameters	:
	df_list	- a list of dataframes with their name attributes already set.
	'''
	for df in df_list :
		df['assembly']	= df.name
	return df_list

def close_hdf5_files(h5_files_list):
	'''
	Closes every file handle in a list of HDF5 file handles.
	Parameters	:
	hdf5_files_list 	- a list containing open HDF5 file handles.
	'''
	for h5_file in h5_files_list:
		h5_file.close()
	return

# ----------------------- Getter functions --------------------------------

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

	if directory == 'plotter':							# If plotter function calls this function, return the plotting subdirectory in th "plots" directory.
		frame 				= inspect.stack()[-1]
		plot_subdir_path	= frame[0].f_code.co_filename
		plot_subdir 		= os.path.splitext(os.path.basename(plot_subdir_path))[0]
		return os.path.join(dir_dict['plots_dir'], str(plot_subdir))
	
	return dir_dict[directory+'_dir']

def get_files(directory, mode='r'):
	'''
	Returns hdf5 file handles for all files in a directory.
	TODO : Test for no hdf5 files
	Parameters :
	directory 	- Path to directory for reading hdf5 files
	mode		- Mode of opening the file; defaults to read mode.
	'''
	return [h5py.File(os.path.join(directory,fname),mode) for fname in os.listdir(directory)]

# ++++++++++++++++++++ Getter functions for hdf5 files dataset columns 

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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

def get_df(fname_list, params_list=['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift']):
	'''
	Returns a concatenated dataframe (of a single assembly type) of all sub-dataframes constructed using individual hdf5 files.
	Parameters	:
	fname_list	- list of several hdf5 file handles,
	params_list	- list of columns to be extracted from hdf5 file datasets, defaults to complete dataset.
	'''
	df = pd.DataFrame()
	if 'redshift' not in params_list:			# Compulsorily add redshift in list of columns.
		params_list.append('redshift')
	df = pd.concat([df]+[get_subdf(fname,params_list) for fname in fname_list], ignore_index=True)
	return df

def get_particle_distribution(df_list,col='redshift'):
	'''
	Returns a single dataframe containing value counts of different values of column "col".
	Parameters	:
	df_list	- a list of dataframes to be concatenated.
	col 	- column which acts as index vor value counts.
	'''
	df 		= pd.DataFrame()
	for subdf in df_list:
		count				= subdf.value_counts(col).sort_index(ascending=True).to_frame(name='counts')
		count['assembly']	= subdf.name
		count[col]			= count.index
		df 					= pd.concat([df,count])
	return df

def get_total_mass_in_particles_with_redshift(df_list):
	'''
	Returns a concatenated dataframe containing redshift, sum of particle masses at a given redshift and assembly mode.
	Parameters	:
	df_list 	- List of dataframes (presumably, with different assembly modes) for which distribution is to be plotted.
	'''
	subdf_list 	= list()
	for subdf in df_list:
		grouped_df 			= subdf.groupby(['redshift'],as_index=False).agg({'mass':sum})
		grouped_df.name 	= subdf.name
		subdf_list.append(grouped_df)
	add_assembly_column(subdf_list)
	return pd.concat(subdf_list)

# ----------------------- Plotter functions --------------------------------

def prepare_plot(context='paper',theme='dark',font_scale=1,rc_kwparams=dict()):
	'''
	Set seaborn styling for plots.
	'''
	rc_params	= {
		'xtick.bottom': True, 
		'ytick.left': True, 
		'axes.spines.left': True, 
		'axes.spines.bottom': True
		}
	rc_params.update(rc_kwparams)			# Update parameters by adding parameters passed to function
	sns.set_context(context,font_scale=font_scale)
	sns.set_style(theme,rc_params)

def plot_or_not(show,plot_name=None,dpi=480,ftype='png',bbox_inches='tight'):
	'''
	Function to switch between show and save methods for plots.
	Parameters	:	
	show		- Shows an 'interactive' plot if True, else saves the plot if False. If set to None, does nothing.
	plot_name	- If set, plot is saved using plot_name if show is set to False. If show is set to False, yet plot_name is not provided, a random number is generated for filename.
	dpi 		- dots per inch in the saved plot, default 480.
	ftype 		- file type of image to save, default .png.
	bbox_inches	- bounding box layout when saving a figure, default tight layout.
	'''
	if show == True :
		plt.show()
	elif show == False :
		if plot_name == None :
			plot_name = np.random.randint(10000,99999)
		plt.savefig(os.path.join(get_directory('plotter'),plot_name+'.'+str(ftype)),
			dpi=dpi,bbox_inches=bbox_inches)
	elif show == None :
		pass 

def make_gif(dir_list,fps):
	'''
	Makes an animated gif from images. The directory should contain only the image files.
	Parameters	:
	dir_list 	- list of directories containing image files. The created animations are saved to same directories by default.
	fps 		- frames per second in the animation.
	'''
	for directory in dir_list:
		images 	= list()
		for fname in sorted(os.listdir(directory),reverse=True):
			path  	= os.path.join(directory,fname)
			images.append(imageio.imread(path))
		imageio.mimsave(os.path.join(directory,'animation.gif'),images,fps=fps)
	return

def redshift_x_axis(ax, ax_primary,zvalues=[0.0,0.125,0.25,0.5,1.0,2.5,5.0,7.0]):
	'''
	Function mapping of ticks for a secondary redshift axis corresponding to cosmological lookback time.
	Parameters:
	ax			= secondary redshift axis, probably defined using matplotlib's twinx/twiny function
	ax_primary 	= primary lookback time axis
	'''
	cosmology = FlatLambdaCDM(100.*0.6777,Om0=0.307,Ob0=0.04825) # Define a flat Lambda-CDM cosmology with parameters mentioned in Schaye et al. 2015
	zvals 		= np.array(zvalues) # Redshift tick values
	time_in_Gyr = cosmology.age(zvals).value 					 # Lookback time corresponding to redshift tick values
	ax.set_xticks(time_in_Gyr)									 # Position ticks at lookback times corresponding to redshift tick values
	ax.set_xticklabels('{:g}'.format(z) for z in zvals)			 # Rename lookback time ticks to corresponding redshift values
	ax.set_xlim(ax_primary.get_xlim())							 # Set equal axis limits for the two x-axes
	return ax


# ============================================ Main program =============================================
