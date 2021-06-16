from plotter import *

# ===================================== User-defined function definitions ===============================


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
	g.ax.invert_xaxis()
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
	g.ax.invert_xaxis()
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
	palette 	= itertools.cycle(sns.color_palette())
	for df in df_list :
		prepare_plot(font_scale=2)
		g = sns.displot(data 	= df,
						x 		= 'mass',
						col  	= 'redshift',
						col_wrap 	= 5,
						col_order 	= sorted(df['redshift'].unique(), reverse=True),
						color 	= next(palette),
						facet_kws	=dict(sharey=False)).set(xlim=[0,1.5e6])
		plot_or_not(show,plot_name='mass_distribution_wrt_redshift_'+str(df.name),dpi=240)
	return

# ============================================ Main Program =============================================

if __name__ == '__main__' :

	# ------- Get HDF5 file handles list for each type of assembly.

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

	# ------- Close HDF5 file handles

	for hdf5_files_list in [gm_early_files,organic_files,gm_late_files]:
		close_hdf5_files(hdf5_files_list)

	# ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	gm_early_df.name  	= 'GM-Early'
	organic_df.name  	= 'Organic'
	gm_late_df.name  	= 'GM-Late'

	df_list 			= [gm_early_df,organic_df,gm_late_df]
	# df_list 			= [organic_df]

	# ------- Plot the particle number distribution for all assembly modes. All plotting functions are standalone.

	# plot_mass_distribution(df_list,show=True)

	# plot_mass_distribution_with_redshift(df_list,show=False)

	plot_particle_distribution(df_list,show=False)

	plot_total_mass_in_particles_with_redshift(df_list,show=False)