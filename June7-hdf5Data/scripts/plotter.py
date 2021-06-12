import numpy as np 
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re 
import time

# ===================================== User-defined function definitions ===============================

# ----------------------- General-purpose functions -----------------------

def capitalize_first_letter(string):
	'''
	Returns a string with capitalized first letter.
	Parameters	:
	string	- a string whose first words is to be capitalized. This does not return string in camelcase, for which we might need to call this function repeated on split string.
	'''
	return string[0].capitalize()+string[1:]

	'''
	Returns a list of dataframes with an extra column added, which signifies mode of assembly of data. This is denoted by name of dataframe.
	Parameters	:
	df_list	- a list of dataframes with their name attributes already set.
	'''
	for df in df_list :
		df['assembly']	= df.name
	return df_list

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
	sns.set_context(context,font_scale)
	sns.set_style(theme,rc_params)

def plot_or_not(show,plot_name=None,dpi=480):
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
		plt.savefig(os.path.join(get_directory('plots'),plot_name+'.png'),
			dpi=dpi,bbox_inches='tight')
	elif show == None :
		pass 

def plot_particle_distribution(df_list, col='redshift',show=True):
	'''
	Plotter function to plot the distribution of star particles with respect to a variable (redshift, by default).
	TODO : Test for other columns passed as parameter col.
	Parameters :
	df_list 	- List of dataframes (assembly modes) for which distribution is to be plotted using separate hues on the same plot.
	col 		- column with respect to which distribution is to be plotted. Redshift is used by default.
	show 		- passed to plot_or_not() function to evaluate whether to show the plot or save it.
	'''
	prepare_plot(font_scale=1.25)
	dist_df			= get_particle_distribution(df_list,col)
	hue				= 'assembly'
	g 				= sns.relplot(data=dist_df,
								  x=col,
								  y='counts',
								  hue= hue,
								  kind='scatter').set(xlabel=capitalize_first_letter(str(col)),
								  ylabel='Counts')		# Plots a scatter plot. To plot histogram instead, use displot with kind='hist'. Counts are then computed automatically.
	g_axis_level 	= sns.lineplot(data=dist_df,
								   x=col,
								   y='counts',
								   hue= hue,
								   ax=g.ax,
								   legend=False)		# Plots a line plot above which scatter points will lie.
	g._legend.set_title(capitalize_first_letter(str(hue)))
	plot_or_not(show,plot_name='particle_distribution_wrt_'+col)
	return  							# No return value. Plot is either shown or saved, or nothing is done.

def plot_total_mass_in_particles_with_redshift(df_list,show=True):
	'''
	Plotter function to plot the distribution of total mass of star particles at a given redshift with respect to redshift.
	Parameters :
	df_list 	- List of dataframes (presumably, with different assembly modes) for which distribution is to be plotted using separate hues on the same plot.
	show 		- passed to plot_or_not() function to evaluate whether to show the plot or save it.
	'''
	prepare_plot(font_scale=1.25)
	dist_df 		= get_total_mass_in_particles_with_redshift(df_list)
	hue 			= 'assembly'
	g 				= sns.relplot(data=dist_df,
								  x='redshift',
								  y='mass',
								  hue= hue,
								  kind='scatter').set(xlabel='Redshift',
								  ylabel='Total Mass $M_{\odot}$')		# Plots a scatter plot. To plot histogram instead, use displot with kind='hist'. Counts are then computed automatically.
	g_axis_level 	= sns.lineplot(data=dist_df,
								   x='redshift',
								   y='mass',
								   hue= hue,
								   ax=g.ax,
								   legend=False)		# Plots a line plot above which scatter points will lie.
	g._legend.set_title(capitalize_first_letter(str(hue)))
	plot_or_not(show,plot_name='total_mass_wrt_redshift')
	return  							# No return value. Plot is either shown or saved, or nothing is done.

def plot_mass_distribution(df_list,show=True):
	'''
	Plots mass distribution for different types of assembly modes. Number of particles at all redshifts are added. Rugplot is added to show range of masses involved.
	Parameters	:
	df_list	- List of dataframes (assembly modes) for which distribution is to be plotted using separate hues on different plots with shared y axes.
	show 	- parameter defining whether to save or show the plot.
	'''
	df_list 	= add_assembly_column(df_list)
	df 			= pd.concat(df_list)
	prepare_plot(theme='darkgrid',font_scale=1.25)
	# print(sns.axes_style())
	hue 		='assembly'
	g 			= sns.displot(data=df,
							  x='mass',
							  col=hue,
							  hue=hue,
							  bins=150,
							  kind='hist',
							  rug=True,
							  rug_kws={'height':-0.025,'clip_on':False,'alpha':0.5},).set(title='')
	g._legend.set_title(capitalize_first_letter(str(hue)))
	for axlist in g.axes:					# Access each axis of FacetGrid
		for ax in axlist:
			ax.tick_params(length=10,pad=10)
			ax.set_xlabel('Mass $[M_{\odot}]$')
	plot_or_not(show,plot_name='particle_mass_distribution')
	return

def plot_mass_distribution_with_redshift(df_list,show=True):
	'''
	Plots mass distribution for a range of redshift for each individual type of assembly mode. 
	Parameters	:
	df_list	- List of dataframes (assembly modes) for which distribution is to be plotted on a separate figure over a range of axes.
	show 	- parameter defining whether to save or show the plot.
	'''
	for df in df_list :
		prepare_plot(font_scale=2)
		g = sns.displot(data 	= df,
						x 		= 'mass',
						col  	= 'redshift',
						col_wrap = 5,
						facet_kws=dict(sharey=False)).set(xlim=[0,1.5e6])
		plot_or_not(show,plot_name='mass_distribution_wrt_redshift_'+str(df.name),dpi=240)
	return


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

	df_list 			= [gm_early_df,organic_df,gm_late_df]
	# df_list 			= [organic_df]

	# ------- Plot the particle number distribution for all assembly modes.

	# plot_mass_distribution(df_list,show=True)

	# plot_mass_distribution_with_redshift(df_list,show=True)

	# plot_particle_distribution(df_list,show=True)

	# plot_total_mass_in_particles_with_redshift(df_list,show=True)