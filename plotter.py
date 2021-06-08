import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl 
import seaborn as sns
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u


			# ------- Define a flat Lambda-CDM cosmology with parameters mentioned in Schaye et al. 2015

cosmology = FlatLambdaCDM(100.*0.6777,Om0=0.307,Ob0=0.04825)

# ===================================== User-defined function definitions ===============================

def redshift_x_axis(ax, ax_primary):
	
	'''
	Function mapping of ticks for a secondary redshift axis corresponding to cosmological lookback time.
	Parameters:
	ax 			= secondary redshift axis, probably defined using matplotlib's twinx/twiny function
	ax_primary 	= primary lookback time axis
	'''

	zvals 		= np.array([0.0,0.125,0.25,0.5,1.0,2.5,5.0,7.0]) # Redshift tick values
	time_in_Gyr = cosmology.age(zvals).value 					 # Lookback time corresponding to redshift tick values
	ax.set_xticks(time_in_Gyr)									 # Position ticks at lookback times corresponding to redshift tick values
	ax.set_xticklabels('{:g}'.format(z) for z in zvals)			 # Rename lookback time ticks to corresponding redshift values
	ax.set_xlim(ax_primary.get_xlim())							 # Set equal axis limits for the two x-axes
	return ax


# ===================================== Main program ====================================================

			# ------- Data directory and file names

data_dir 		= './data_updated/'
gm_early_fname 	= 'halo_catalogue_gm_early.txt'
organic_fname	= 'halo_catalogue_organic.txt'
gm_late_fname 	= 'halo_catalogue_gm_late.txt'

			# ------- Dataframe columns

cols = [
	'index', 
	'time', 
	'a_exp', 
	'redshift', 
	'subhalo_centre_x',
	'subhalo_centre_y',
	'subhalo_centre_z', 
	'subhalo_peculiar_velocity_x',
	'subhalo_peculiar_velocity_y',
	'subhalo_peculiar_velocity_z', 
	'Halo_mass', 
	'Stellar_mass', 
	'BH_mass', 
	'SFR', 
	'sSFR'
]

			# ------- Plot axis titles corresponding to dataframe columns

title = {
	'Stellar_mass':'$\log_{10}$ (Stellar Mass)',
	'Halo_mass':'$\log_{10}$ (Halo Mass)',
	'BH_mass':'$\log_{10}$ (Blackhole Mass)',
	'SFR':'$\log_{10}$ (Star Formation Rate)',
	'sSFR':'$\log_{10}$ (Specific Star Formation Rate)',
	'subhalo_peculiar_velocity':'$\log_{10}$ (Subhalo Peculiar Velocity)'
}

			# ------- Units of measurement corresponding to dataframe columns

units = {
	'Stellar_mass':'$M_{\odot}$',
	'Halo_mass':'$M_{\odot}$',
	'BH_mass':'$M_{\odot}$',
	'SFR':'$M_{\odot} yr^{-1}$',
	'sSFR':'$yr^{-1}$',
	'time':'Gyr',
	'subhalo_peculiar_velocity':'$km s^{-1}$'
}

			# ------- Dataframe columns to be plotted on log scale

# Compute log10 values instead of using plt.scale("log") for correct axis ticks.

log10_unit_variables = [
	'Stellar_mass',
	'BH_mass',
	'Halo_mass',
	'SFR',
	'sSFR',
	'subhalo_peculiar_velocity'
]

			# ------- Values to replace zero values in a column with

zero_values = {
	'SFR':1e-3,
	'sSFR':1e-13,
	'BH_mass':1
}

			# ------- Dataframe columns to be plotted with respect to time

yaxes = ['Stellar_mass','Halo_mass','BH_mass','SFR','sSFR','subhalo_peculiar_velocity']
# yaxes 	= ['SFR']

			# ------- Read and categorise data on the basis of assembly mode - GM-early, organic and GM-late

gm_early_data 	= pd.read_fwf(data_dir+gm_early_fname,skiprows=1,names=cols)
organic_data 	= pd.read_fwf(data_dir+organic_fname,skiprows=1,names=cols)
gm_late_data 	= pd.read_fwf(data_dir+gm_late_fname,skiprows=1,names=cols)

gm_early_data['assembly'] 	= 'gm-early'
gm_early_data['assembly'] 	= gm_early_data['assembly'].astype('category')
organic_data['assembly'] 	= 'organic'
organic_data['assembly'] 	= organic_data['assembly'].astype('category')
gm_late_data['assembly'] 	= 'gm-late'
gm_late_data['assembly'] 	= gm_late_data['assembly'].astype('category')

			# ------- Calculate Euclid norm for subhalo peculiar velocity

gm_early_data['subhalo_peculiar_velocity']	= np.sqrt(
	(gm_early_data['subhalo_peculiar_velocity_x'])**2 + 
	(gm_early_data['subhalo_peculiar_velocity_y'])**2 + 
	(gm_early_data['subhalo_peculiar_velocity_z'])**2 )
organic_data['subhalo_peculiar_velocity'] 	= np.sqrt(
	(organic_data['subhalo_peculiar_velocity_x'])**2 + 
	(organic_data['subhalo_peculiar_velocity_y'])**2 + 
	(organic_data['subhalo_peculiar_velocity_z'])**2 )
gm_late_data['subhalo_peculiar_velocity'] 	= np.sqrt(
	(gm_late_data['subhalo_peculiar_velocity_x'])**2 + 
	(gm_late_data['subhalo_peculiar_velocity_y'])**2 + 
	(gm_late_data['subhalo_peculiar_velocity_z'])**2 )

			# ------- Replace zero values in columns

for key in zero_values:
	gm_early_data[key].replace(float(0),zero_values[key],inplace=True)
	organic_data[key].replace(float(0),zero_values[key],inplace=True)
	gm_late_data[key].replace(float(0),zero_values[key],inplace=True)

			# ------- Convert to log10 values for each column in list log10_unit_variables

for variable in log10_unit_variables :
	gm_early_data[variable] = gm_early_data[variable].map(np.log10)
	organic_data[variable]	= organic_data[variable].map(np.log10)
	gm_late_data[variable]	= gm_late_data[variable].map(np.log10)

			# ------- Concatenate all separate dataframes to one
			# ------- Declare empty dataframe for tabulation deviations of GM cases from organic cases

df 			= pd.concat([gm_early_data,organic_data,gm_late_data])
dev_df 		= pd.DataFrame()

			# ------- Loop for every dataframe column to be plotted

for y in yaxes:  

			# ------- Calculate deviation and store in new dataframe column
	gm_early_data_dev 	= gm_early_data[y] - organic_data[y]
	gm_late_data_dev 	= gm_late_data[y] - organic_data[y]
	organic_data_dev 	= organic_data[y] - organic_data[y]
	dev_df[str(y)+'_deviation'] = pd.concat([gm_early_data_dev,organic_data_dev,gm_late_data_dev])

			# ------- Set seaborn style parameters
			# ------- use print(sns.axes_style()) for getting a complete list of style attributes
	sns.set_context('paper')
	sns.set_style('dark',
		{'xtick.bottom': True, 
		'ytick.left': True, 
		'axes.spines.left': True, 
		'axes.spines.bottom': True})				

			# ------- Define matplotlib figure and axis followed by seaborn plotting.

	fig,axes = plt.subplots(2,1,
		figsize=(5,7),
		sharex=True,								# shared x-axis for main and deviation plots
		gridspec_kw=dict(height_ratios=[5,1]))		# two vertically stacked subplots

	sns.lineplot(data=df, 
		x='time', 
		y=y,
		hue='assembly',
		ax=axes[0])									# plot main values on first axis

	axes[0].set_ylabel(title[y]+' ['+units[y]+']', 
		fontsize=12.5, 
		labelpad=5)

	sns.lineplot(data=df, 
		x='time', 
		y=dev_df[str(y)+'_deviation'],
		hue='assembly',
		ax=axes[1], 								# plot deviations from organic on second (smaller) axis
		legend=False)								# remove legend, already shown for axis[0]

	# axes[0].legend(loc='lower right')				# manually set legend location (for consistency, optional)

	axes[1].set_ylabel('Deviation', 
		fontsize=12.5, 
		labelpad=5)

	axes[1].set_xlabel('Time ['+units['time']+']', 
		fontsize=12.5, 
		labelpad=5)

	secondary_x = axes[0].twiny()					# define secondary redshift axis
	redshift_x_axis(secondary_x,axes[0])			# function call to obtain correct mapping of time and redshift values

	fig.subplots_adjust(hspace=0.05)				# reduce space between subplots for better viewing of ticks
	
	# plt.savefig('./plots/'+str(y)+'.png',dpi=480,bbox_inches='tight')

plt.show()											# Comment out if using savefig to save plots instead of viewing