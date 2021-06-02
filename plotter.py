import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl 
import seaborn as sns

data_dir 		= './data/'
gm_early_fname 	= 'halo_catalogue_gm_early.txt'
organic_fname	= 'halo_catalogue_organic.txt'
gm_late_fname 	= 'halo_catalogue_gm_late.txt'

cols = ['index', 
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
'sSFR']

title = {'Stellar_mass':'$\log_{10}$ (Stellar Mass)',
'Halo_mass':'$\log_{10}$ (Halo Mass)',
'BH_mass':'$\log_{10}$ (Blackhole Mass)',
'SFR':'$\log_{10}$ (Star Formation Rate)',
'sSFR':'$\log_{10}$ (Specific Star Formation Rate)',
'subhalo_peculiar_velocity':'$\log_{10}$ (Subhalo Peculiar Velocity)'}

units = {'Stellar_mass':'$M_{\odot}$',
'Halo_mass':'$M_{\odot}$',
'BH_mass':'$M_{\odot}$',
'SFR':'$M_{\odot} yr^{-1}$',
'sSFR':'$yr^{-1}$',
'time':'Gyr',
'subhalo_peculiar_velocity':'$km s^{-1}$'}

log10_unit_variables = ['Stellar_mass',
'BH_mass',
'Halo_mass',
'SFR',
'sSFR',
'subhalo_peculiar_velocity'
]

yaxes = ['Stellar_mass','Halo_mass','BH_mass','SFR','sSFR','subhalo_peculiar_velocity']
# yaxes = ['BH_mass']

gm_early_data 	= pd.read_fwf(data_dir+gm_early_fname,skiprows=1,names=cols)
organic_data 	= pd.read_fwf(data_dir+organic_fname,skiprows=1,names=cols)
gm_late_data 	= pd.read_fwf(data_dir+gm_late_fname,skiprows=1,names=cols)

gm_early_data['assembly'] 	= 'gm-early'
gm_early_data['assembly'] 	= gm_early_data['assembly'].astype('category')
organic_data['assembly'] 	= 'organic'
organic_data['assembly'] 	= organic_data['assembly'].astype('category')
gm_late_data['assembly'] 	= 'gm-late'
gm_late_data['assembly'] 	= gm_late_data['assembly'].astype('category')

gm_early_data['subhalo_peculiar_velocity'] = np.sqrt((gm_early_data['subhalo_peculiar_velocity_x'])**2 + 
	(gm_early_data['subhalo_peculiar_velocity_y'])**2 + 
	(gm_early_data['subhalo_peculiar_velocity_z'])**2 )
organic_data['subhalo_peculiar_velocity'] = np.sqrt((organic_data['subhalo_peculiar_velocity_x'])**2 + 
	(organic_data['subhalo_peculiar_velocity_y'])**2 + 
	(organic_data['subhalo_peculiar_velocity_z'])**2 )
gm_late_data['subhalo_peculiar_velocity'] = np.sqrt((gm_late_data['subhalo_peculiar_velocity_x'])**2 + 
	(gm_late_data['subhalo_peculiar_velocity_y'])**2 + 
	(gm_late_data['subhalo_peculiar_velocity_z'])**2 )

for variable in log10_unit_variables :
	gm_early_data[variable] = gm_early_data[variable].map(np.log10)
	organic_data[variable]	= organic_data[variable].map(np.log10)
	gm_late_data[variable]	= gm_late_data[variable].map(np.log10)

df 			= pd.concat([gm_early_data,organic_data,gm_late_data])
dev_df 		= pd.DataFrame()

for y in yaxes:  
	gm_early_data_dev 	= gm_early_data[y] - organic_data[y]
	gm_late_data_dev 	= gm_late_data[y] - organic_data[y]
	organic_data_dev 	= organic_data[y] - organic_data[y]
	dev_df[str(y)+'_deviation'] = pd.concat([gm_early_data_dev,organic_data_dev,gm_late_data_dev])

	sns.set_context('paper')
	sns.set_style('darkgrid',{'axes.facecolor':'#f5eed7','xtick.bottom': True, 'ytick.left': True, 'axes.spines.left': True, 'axes.spines.bottom': True})

	fig,axes = plt.subplots(2,1,figsize=(5,7),sharex=True,gridspec_kw=dict(height_ratios=[5,1]))
	sns.lineplot(data=df, x='time', y=y,hue='assembly',ax=axes[0])
	axes[0].set_ylabel(title[y]+' ['+units[y]+']', fontsize=12.5, labelpad=5)
	sns.lineplot(data=df, x='time', y=dev_df[str(y)+'_deviation'],hue='assembly',ax=axes[1], legend=False)
	axes[0].legend(loc='lower right')
	axes[1].set_ylabel('Deviation', fontsize=12.5, labelpad=5)
	axes[1].set_xlabel('Time ['+units['time']+']', fontsize=12.5, labelpad=5)
	plt.savefig('./plots/'+str(y)+'.png',dpi=480,bbox_inches='tight')
# plt.show()



