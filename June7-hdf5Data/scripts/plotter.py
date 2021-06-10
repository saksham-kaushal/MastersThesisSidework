import numpy as np 
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re 

def get_directory(directory):
	dir_dict 						= dict()
	dir_dict['scripts_dir'] 		= os.path.dirname(os.path.abspath(__file__))
	dir_dict['root_dir'] 			= os.path.dirname(dir_dict['scripts_dir'])
	dir_dict['data_dir'] 			= os.path.join(dir_dict['root_dir'], 'data')
	dir_dict['gm_early_data_dir']	= os.path.join(dir_dict['data_dir'], 'gm_early')
	dir_dict['organic_data_dir']	= os.path.join(dir_dict['data_dir'], 'organic')
	dir_dict['gm_late_data_dir']	= os.path.join(dir_dict['data_dir'], 'gm_late')
	dir_dict['plots_dir']			= os.path.join(dir_dict['root_dir'], 'plots')
	return dir_dict[directory+'_dir']

def get_files(directory):
	return [h5py.File(os.path.join(directory,fname),'r') for fname in os.listdir(directory)]

def get_coords_x(fname):
	return fname["Coordinates"][:,0]

def get_coords_y(fname):
	return fname["Coordinates"][:,1]

def get_coords_z(fname):
	return fname["Coordinates"][:,2]

def get_vel_x(fname):
	return fname["Velocity"][:,0]

def get_vel_y(fname):
	return fname["Velocity"][:,1]

def get_vel_z(fname):
	return fname["Velocity"][:,2]

def get_mass(fname):
	return fname["Mass"]

def get_redshift(fname):
	return float(re.compile(r'z\d+p\d+').findall(str(fname.file))[0].replace('z','').replace('p','.'))

def get_subdf(fname,params_list):
	subdf = pd.DataFrame()
	params_dict		= {
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
	df = pd.DataFrame()
	if 'redshift' not in params_list:
		params_list.append('redshift')
	df = pd.concat([df]+[get_subdf(fname,params_list) for fname in fname_list], ignore_index=True)
	return df

def prepare_plot():
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
	prepare_plot()
	dist_df		= get_particle_distribution(df_list,col)
	sns.relplot(data=dist_df,x=col,y='counts',hue='assembly',kind='scatter')
	plot_or_not(show,plot_name='particle_distribution_wrt_'+col)


def plot_masses(hdf5_file_list):
	prepare_plot()

def plot_or_not(show,plot_name=None):
	if show == True :
		plt.show()
	elif show == False :
		if plot_name == None :
			plot_name = np.random.randint(10000,99999)
		plt.savefig(os.path.join(get_directory(plots),plot_name,'.png'),
			dpi=480,bbox_inches='tight')
	elif show == None :
		pass 

if __name__ == '__main__' :
	gm_early_files		= get_files(get_directory('gm_early_data'))
	organic_files		= get_files(get_directory('organic_data'))
	gm_late_files		= get_files(get_directory('gm_late_data'))

	gm_early_df 		= get_df(gm_early_files,['mass','coords_x'])
	organic_df 			= get_df(organic_files,['mass','coords_x'])
	gm_late_df 			= get_df(gm_late_files,['mass','coords_x'])

	gm_early_df.name  	= 'GM-Early'
	organic_df.name  	= 'Organic'
	gm_late_df.name  	= 'GM-Late'

	plot_particle_distribution([gm_early_df,organic_df,gm_late_df],show=True)
	# plt.show()
