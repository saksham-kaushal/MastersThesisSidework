from plotter import *

# ===================================== User-defined function definitions ===============================

def set_square_aspect(grid):
	'''
	Sets square shape for the plotting axes of a FacetGrid type plot and returns the modified grid. Unlike FacetGrid's aspect parameter, this does not consider the axes labels and ticks while making the figure square shaped.
	Parameters	:
	grid 	- a seaborn FacetGrid instance, or similar grid plot with axes attribute.
	'''
	for axis in grid.axes:
		for ax in axis:
			ax.set_box_aspect(1)
	return grid

def plot_galaxy_coordinates_combined(df_list, show=True):
	'''
	Plots the coordinates of star particles in a galaxy along x-y, y-z and z-x planesover different redshifts in a single image with plots in different rows corresponding to plots at different redshifts. Three different figures are plotted for the three different planes.
	Parameters	:
	df_list 	- list of dataframes constructed from hdf5 files, containging x, y and z coordinates of particles.
	show 		- parameter to control showing or saving of plots, defaults to showing.
	'''
	palette 	= itertools.cycle(sns.color_palette())
	row 		= 'redshift'
	for df in df_list:
		prepare_plot(font_scale=1.4)
		color  		= next(palette)
		g1 	= sns.relplot(
					data		=df,
					x 			= 'coords_x',
					y  			= 'coords_y',
					kind		= 'scatter',
					row			= row,
					row_order 	= sorted(df[row].unique(), reverse=True),
					color 		= color,
					s  			= 1,
					height 		= 3
					).set(xlabel='x [Mpc]',ylabel='y [Mpc]')
		g1.set_titles(row_template=capitalize_first_letter(row)+'= {row_name}')
		g1 	= set_square_aspect(g1)
		g1.fig.subplots_adjust(hspace=0)
		plot_or_not(show=show,plot_name=str(df.name)+'_x-y_coordinates_evolution',dpi =240)

		prepare_plot(font_scale=1.4)
		g2 	= sns.relplot(
					data		=df,
					x 			= 'coords_y',
					y  			= 'coords_z',
					kind		= 'scatter',
					row			= row,
					row_order 	= sorted(df[row].unique(), reverse=True),
					color 		= color,
					s  			= 1,
					height 		= 3
					).set(xlabel='y [Mpc]',ylabel='z [Mpc]')
		g2.set_titles(row_template=capitalize_first_letter(row)+'= {row_name}')
		g2 	= set_square_aspect(g2)
		g2.fig.subplots_adjust(hspace=0)
		plot_or_not(show=show,plot_name=str(df.name)+'_y-z_coordinates_evolution',dpi =240)

		prepare_plot(font_scale=1.4)
		g3 	= sns.relplot(
					data		=df,
					x 			= 'coords_z',
					y  			= 'coords_x',
					kind		= 'scatter',
					row			= row,
					row_order 	= sorted(df[row].unique(), reverse=True),
					color 		= color,
					s  			= 1,
					height 		= 3
					).set(xlabel='z [Mpc]',ylabel='x [Mpc]')
		g3.set_titles(row_template=capitalize_first_letter(row)+'= {row_name}')
		g3 	= set_square_aspect(g3)
		g3.fig.subplots_adjust(hspace=0)
		plot_or_not(show=show,plot_name=str(df.name)+'_z-x_coordinates_evolution',dpi =240)
		# plot_or_not(show=True)
	return

def plot_galaxy_coordinates_individual(df_list, show=True):
	'''
	Plots the coordinates of star particles of a galaxy in x, y and z axes combined in a single 3D plot for each redshift value. 
	Parameters	:
	df_list 	- list of dataframes constructed from hdf5 files, containging x, y and z coordinates of particles.
	show 		- parameter to control showing or saving of plots, defaults to showing.
	'''
	palette 	= itertools.cycle(sns.color_palette())
	row 		= 'redshift'
	for df in df_list:
		color 		= next(palette)
		grouped_df 	= df.groupby(row)
		for row_value in grouped_df[row].unique():
			prepare_plot()
			row_value_group 	= grouped_df.get_group(row_value[0])
			fig, axes 			= plt.subplots(1,3,figsize=(15,5.5))
			for ax in axes:
				ax.set_box_aspect(1)
				ax.set_xlim([-0.03,0.03])
				ax.set_ylim([-0.03,0.03])

			sns.scatterplot(
				data	= row_value_group,
				x 		= 'coords_x',
				y 		= 'coords_y',
				color 	= color,
				s 		= 1,
				ax  	= axes[0],
				).set(xlabel='x [Mpc]',ylabel='y [Mpc]')
			axes[0].set_title('View in x-y plane',loc='left')

			sns.scatterplot(
				data	= row_value_group,
				x 		= 'coords_y',
				y 		= 'coords_z',
				color 	= color,
				s 		= 1,
				ax  	= axes[1],
				).set(xlabel='y [Mpc]',ylabel='z [Mpc]')
			axes[1].set_yticklabels([])
			axes[1].set_title('View in y-z plane',loc='left')

			sns.scatterplot(
				data	= row_value_group,
				x 		= 'coords_z',
				y 		= 'coords_x',
				color 	= color,
				s 		= 1,
				ax  	= axes[2],
				).set(xlabel='z [Mpc]',ylabel='x [Mpc]')
			axes[2].set_yticklabels([])
			axes[2].set_title('View in z-x plane',loc='left')

			fig.suptitle(str(df.name)+', '+capitalize_first_letter(row)+'='+str(row_value[0]))
			fig.tight_layout()
			plot_or_not(show=show, plot_name=os.path.join(str(df.name),str(row_value[0])),dpi=240)
	return

def plot_galaxy_coordinates_3d(df_list,show=True):
	'''
	Plots the coordinates of star particles of a galaxy along x-y, y-z and z-x planes combined in a single image for each redshift value. 
	Parameters	:
	df_list 	- list of dataframes constructed from hdf5 files, containging x, y and z coordinates of particles.
	show 		- parameter to control showing or saving of plots, defaults to showing.
	'''
	palette 	= itertools.cycle(sns.color_palette())
	row 		= 'redshift'
	for df in df_list:
		color 		= next(palette)
		grouped_df 	= df.groupby(row)
		for row_value in grouped_df[row].unique():
			plt.style.use('seaborn-whitegrid')
			row_value_group 	= grouped_df.get_group(row_value[0])
			fig 				= plt.figure(figsize=(7,7))
			axes 				= fig.add_subplot(projection='3d')
			kwargs={'marker':'^'}
			axes.scatter(
				row_value_group['coords_x'],
				row_value_group['coords_y'],
				row_value_group['coords_z'],
				s 			= 0.1,
				color 		=color,
				**kwargs)
			axes.set_xlim([-0.03,0.03])
			axes.set_xlabel('x [Mpc]')
			axes.set_ylim([-0.03,0.03])
			axes.set_ylabel('y [Mpc]')
			axes.set_zlim([-0.03,0.03])
			axes.set_zlabel('z [Mpc]')
			fig.suptitle(str(df.name)+', '+capitalize_first_letter(row)+'='+str(row_value[0]), y=0.92)
			fig.tight_layout()
			plot_or_not(show=show, plot_name=os.path.join('3d',str(df.name),str(row_value[0])))
	return


# ============================================ Main Program =============================================

if __name__ == '__main__' :
	# ------- Get hdf5 file handles list for each type of assembly.

	gm_early_files		= get_files(get_directory('gm_early_data'))
	organic_files		= get_files(get_directory('organic_data'))
	gm_late_files		= get_files(get_directory('gm_late_data'))

	# ------- Available fields for cols :
	# ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['coords_x', 'coords_y', 'coords_z']

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

	# plot_galaxy_coordinates_combined(df_list)

	# plot_galaxy_coordinates_individual(df_list)

	# plot_galaxy_coordinates_3d(df_list)

	gif_dirs 			= [
							'GM-Early',
							'Organic',
							'GM-Late',
							os.path.join('3d','GM-Early'),
							os.path.join('3d','Organic'),
							os.path.join('3d','GM-Late')
						  ]

	dir_list			= [os.path.join(get_directory('plotter'),subdir) for subdir in gif_dirs]
	
	make_gif(dir_list,fps=3)