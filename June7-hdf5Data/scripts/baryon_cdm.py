from plotter import *
from angular_momentum import *
from matplotlib.ticker import FuncFormatter

pd.set_option('display.expand_frame_repr', False)

# ===================================== User-defined function definitions ===============================

def add_expansion_factor_column(df):
	df['expansion_factor'] 	= 1/(1+df['redshift'])
	return df

def plot_expansion_factor_net_specific_angular_momentum_a3by2(df,assembly_names,show=True,suffix=''):
	prepare_plot(font_scale=1.25)
	col = 'assembly'
	df['net_specific_angular_momentum_a3by2']  = df['net_specific_angular_momentum']/(df['expansion_factor'])**(1.5)
	g 	= sns.relplot(
					data 	= df, 
					x 		= 'expansion_factor',
					y  		= 'net_specific_angular_momentum_a3by2',
					col  	= col,
					hue  	= col,
					kind  	= 'scatter',
					legend 	= False,
					col_order = ['gm_early','organic','gm_late'],
					hue_order = ['gm_early','organic','gm_late']
					).set(
					xlabel 	= get_plot_axes_titles()['expansion_factor'],
					ylabel 	= '$L_{specific}/a^{3/2}$')
	for col_val,ax in g.axes_dict.items():
		col_val  	= assembly_names[col_val]
		ax.set_title(str(col_val)+' '+capitalize_first_letter(col),fontdict={'fontsize':15})

	for axes in g.axes:
		for axis in axes:
			powerlaw_x = np.linspace(axis.get_xlim()[0],axis.get_xlim()[1],15)
			powerlaw_y = np.repeat(10,15)
			axis.plot(
							powerlaw_x,
							powerlaw_y,
							color='black',
							ls='--',
							lw=0.5
						)
			axis.loglog()

	plot_or_not(show=show,plot_name='expansion_factor-net_specific_angular_momentum_a3by2',suffix=suffix)
	return

def plot_net_specific_angular_momentum_combined(df,assembly_names,particles_titles,show=True,suffix=''):
	prepare_plot(font_scale=1.25)
	col 	= 'assembly'
	assembly_order = ['gm_early','organic','gm_late']
	hue 	= 'group'
	sns.set_palette(sns.color_palette(['#333333','#009fe5']))
	palette_order  	= ['PartType1','PartType0+PartType4']
	g 		= sns.relplot(
							data 	= df, 
							x 		= 'expansion_factor',
							y  		= 'net_specific_angular_momentum',
							col  	= col,
							hue  	= hue,
							kind  	= 'line',
							col_order = assembly_order,
							hue_order = palette_order
							).set(
							xlabel 	= get_plot_axes_titles()['expansion_factor'],
							ylabel 	= get_plot_axes_titles()['net_specific_angular_momentum']
						)
	g._legend.set_title('')
	for i,text_entry in enumerate(g._legend.texts):
		text_entry.set_text(particles_titles[palette_order[i]])

	for col_val,ax in g.axes_dict.items():
		col_val  	= assembly_names[col_val]
		ax.set_title(str(col_val)+' '+capitalize_first_letter(col),fontdict={'fontsize':15})

	for axes in g.axes:
		for axis in axes:
			powerlaw_x = np.linspace(axis.get_xlim()[0],axis.get_xlim()[1],15)
			powerlaw_y = powerlaw_x**(1.5)
			axis.plot(
							powerlaw_x,
							powerlaw_y,
							color='black',
							ls='--',
							lw=0.5
						)
			axis.loglog()
		for ax in [axis.xaxis, axis.yaxis]:
			formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
			ax.set_major_formatter(formatter)

	plot_or_not(show=show,plot_name='expansion_factor-net_specific_angular_momentum_combined',suffix=suffix)
	return

def plot_net_specific_angular_momentum_combined_with_inner_dm_halo(df,assembly_names,particles_titles,show=True,suffix=''):
	prepare_plot(font_scale=1.25)
	col 	= 'assembly'
	assembly_order = ['gm_early','organic','gm_late']
	hue 	= 'group'
	sns.set_palette(sns.color_palette(['#333333','#009fe5','#07dcbe']))
	palette_order  	= ['PartType1','PartType0+PartType4','Inner_Halo_PartType1']
	g 			= sns.FacetGrid(
								data 	= df,
								col  	= col,
								hue  	= hue,
								col_order 	= assembly_order,
								hue_order 	= palette_order,
								height  = 5
								)
	g.map(sns.lineplot, 'expansion_factor', 'net_specific_angular_momentum')
	g.map(sns.scatterplot, 'expansion_factor', 'net_specific_angular_momentum')

	g.add_legend()
	g.set(
			xlabel 	= get_plot_axes_titles()['expansion_factor'],
			ylabel 	= get_plot_axes_titles()['net_specific_angular_momentum']
		)
	g._legend.set_title('')
	for i,text_entry in enumerate(g._legend.texts):
		text_entry.set_text(particles_titles[palette_order[i]])

	for col_val,ax in g.axes_dict.items():
		col_val  	= assembly_names[col_val]
		ax.set_title(str(col_val)+' '+capitalize_first_letter(col),fontdict={'fontsize':15})

	for axes in g.axes:
		for axis in axes:
			powerlaw_x = np.linspace(0.1,1.0,15)
			powerlaw_y = powerlaw_x**(1.5)
			axis.plot(
							powerlaw_x,
							powerlaw_y,
							color='black',
							ls='--',
							lw=0.5
						)
			axis.loglog()
			axis.set_xlim(0.1,1.1)
		for ax in [axis.xaxis, axis.yaxis]:
			formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
			ax.set_major_formatter(formatter)

	plot_or_not(show=show,plot_name='expansion_factor-net_specific_angular_momentum_combined_inner_dm_halo',suffix=suffix)
	return

def plot_median_radius(df,assembly_names,particles_titles,show=True,suffix=''):
	prepare_plot(font_scale=1.25)
	col 			= 'assembly'
	assembly_order 	= ['gm_early','organic','gm_late']
	hue 			= 'group'
	sns.set_palette(sns.color_palette(['#333333','#009fe5','#07dcbe']))
	palette_order  	= ['PartType1','PartType0+PartType4','Inner_Halo_PartType1']


	g 				= sns.FacetGrid(
									data 	= df,
									col  	= col,
									hue  	= hue,
									col_order 	= assembly_order,
									hue_order 	= palette_order,
									height  = 5
									)
	g.map_dataframe(sns.lineplot, x='expansion_factor', y='median_radius')
	g.map_dataframe(sns.scatterplot, x='expansion_factor', y='median_radius')

	g.add_legend()
	g.set(
			xlabel 	= get_plot_axes_titles()['expansion_factor'],
			ylabel 	= get_plot_axes_titles()['median_radius']
		)
	g._legend.set_title('')
	for i,text_entry in enumerate(g._legend.texts):
		text_entry.set_text(particles_titles[palette_order[i]])

	for col_val,ax in g.axes_dict.items():
		col_val  	= assembly_names[col_val]
		ax.set_title(str(col_val)+' '+capitalize_first_letter(col),fontdict={'fontsize':15})

	for axes in g.axes:
		for axis in axes:
			axis.loglog()
			axis.set_xlim(0.1,1.1)
		for ax in [axis.xaxis, axis.yaxis]:
			formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
			ax.set_major_formatter(formatter)	

	plot_or_not(show=show,plot_name='expansion_factor-median_radius',suffix=suffix)
	
	return

# ============================================ Main Program =============================================

if __name__ == '__main__' :
	assembly_list 			= ['gm_early','organic','gm_late']
	# assembly_list 				= ['organic']

	# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

	files  = {str(assembly):get_files(get_directory(str(assembly)+'_baryon_cdm_data')) for assembly in assembly_list}

	# ------- Available fields for cols :
	# ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass']

	# ------- Particle type 0 = Gas, 1 = Dark Matter, and 4 = Stars

	particles  			= ['PartType0', 'PartType1', 'PartType4']
	# particles  			= ['PartType0', 'PartType4']
	# particles  			= ['PartType1']
	particles_titles  	= {
							'PartType0':'Gas', 
							'PartType1':'Dark Matter Halo', 
							'PartType4':'Stars', 
							'PartType0+PartType4':'Galaxy (Baryonic Matter)',
							'Inner_Halo_PartType1':'Inner Dark Matter Halo'
							}

	# ------- Get dataframe containing data for each type of assembly mode and store in a dictionary.

	df 		= {str(assembly):get_df(files[str(assembly)],cols,particles) for assembly in assembly_list}

	# ------- Close HDF5 file handles

	for assembly in files.keys():
		close_hdf5_files(files[assembly])
	# ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	assembly_names  		= {
								'gm_early':'GM-Early',
								'organic':'Organic',
								'gm_late':'GM-Late'
							  }

	for assembly in assembly_list:
		df[assembly].name  	= assembly_names[assembly]	

	computed_columns  		= [
								'distance',
								'specific_angular_momentum',
								'net_specific_angular_momentum',
								'angular_momentum',
								'net_angular_momentum',
								'total_kinetic_energy',
								'assembly',
								'group',
								'redshift'
							  ]

	computed_masked_columns = [
								'median_radius',
								'net_specific_angular_momentum',
								'net_angular_momentum',
								'assembly',
								'group',
								'redshift'
							  ]

	computed_values_df  	= pd.DataFrame(columns=computed_columns)
	
	inner_halo_radius_Mpc  	= 0.02

	computed_values_assembly_list  	= list()
	computed_masked_values_assembly_list = list()
	units_flag  					= False
	units_dict 						= dict()
	zero_redshift_mass_specificL	= dict()
	for assembly in assembly_list:
		computed_values_list  		= list()
		computed_masked_values_list = list()
		for z in df[assembly]['redshift'].unique():
			supdf 	= df[assembly].groupby('redshift').get_group(z)
			for particle_type in particles:
				subdf 	= supdf.groupby('group').get_group(particle_type)
				coords, velocities, masses 		= get_numpy_arrays(subdf)
				distance 						= get_norm(coords)
				# r_median						= np.median(distance)
				rxv 							= get_rxv(coords,velocities)
				angular_momentum  				= get_angular_momentum(rxv,masses)
				net_angular_momentum  			= get_norm(angular_momentum)
				specific_angular_momentum		= get_specific_angular_momentum(angular_momentum,masses)
				net_specific_angular_momentum	= get_norm(specific_angular_momentum)
				total_kinetic_energy  			= get_kinetic_energy(velocities,masses)
				computed_values_list.append([
												distance.value,
												specific_angular_momentum.value,
												net_specific_angular_momentum.value,
												angular_momentum.value,
												net_angular_momentum.value,
												total_kinetic_energy.value,
												assembly,
												particle_type,
												z
											  ])
				if particle_type=="PartType1":
					distance 								= get_norm(coords).value
					distance_bool 							= distance <= inner_halo_radius_Mpc
					r_median								= np.median(distance[distance_bool])
					coords_masked  							= coords[distance_bool,:]
					velocities_masked						= velocities[distance_bool,:]
					masses_masked 							= masses[distance_bool]
					rxv_masked 								= get_rxv(coords_masked,velocities_masked)
					angular_momentum_masked  				= get_angular_momentum(rxv_masked,masses_masked)
					net_angular_momentum_masked  			= get_norm(angular_momentum_masked)
					specific_angular_momentum_masked		= get_specific_angular_momentum(angular_momentum_masked,masses_masked)
					net_specific_angular_momentum_masked	= get_norm(specific_angular_momentum_masked)
					computed_masked_values_list.append([ 
														r_median,
														net_specific_angular_momentum_masked.value,
														net_angular_momentum_masked.value,
														assembly,
														particle_type,
														z
													  ])


				if not units_flag:
					units_dict['coords']  						= coords.unit
					units_dict['velocities']	  				= velocities.unit
					units_dict['masses']  						= masses.unit
					units_dict['specific_angular_momentum']		= specific_angular_momentum.unit
					units_dict['net_specific_angular_momentum']	= net_specific_angular_momentum.unit
					units_dict['angular_momentum']  			= angular_momentum.unit
					units_dict['net_angular_momentum']  		= net_angular_momentum.unit
					units_dict['total_kinetic_energy']  		= total_kinetic_energy.unit
					
					units_flag 									= True

				if z==0.0:
					zero_redshift_mass_specificL[assembly] 	= [np.sum(masses), net_specific_angular_momentum]
		computed_values_assembly_df  						= pd.DataFrame(
																			computed_values_list,
																			# orient='index',
																			columns=computed_columns
																			)
		computed_masked_values_assembly_df  				= pd.DataFrame(
																			computed_masked_values_list,
																			# orient='index',
																			columns=computed_masked_columns
																			)
		# computed_values_assembly_df['redshift'] 			= computed_values_assembly_df.index
		add_expansion_factor_column(computed_values_assembly_df)
		add_expansion_factor_column(computed_masked_values_assembly_df)
		computed_values_assembly_list.append(computed_values_assembly_df.convert_dtypes())
		computed_masked_values_assembly_list.append(computed_masked_values_assembly_df.convert_dtypes())

	computed_values_df 			= pd.concat(computed_values_assembly_list,ignore_index=True)
	computed_masked_values_df  	= pd.concat(computed_masked_values_assembly_list,ignore_index=True)
	print(computed_values_df.groupby('group').groups.keys())

	# net_angular_momentum_redshift = computed_values_df.groupby(['redshift','assembly'],as_index=False)['net_angular_momentum'].sum().rename(columns={'sum':'net_angular_momentum','redshift':'redshift','assembly':'assembly'})
	# net_angular_momentum_redshift['expansion_factor'] = 1/(1+net_angular_momentum_redshift['redshift'])
	# plot_expansion_factor_net_angular_momentum(net_angular_momentum_redshift,assembly_names,powerlaw_offset=1e9,show=False,suffix='_cdm')

	# net_specific_angular_momentum_redshift = computed_values_df.groupby(['redshift','assembly'],as_index=False)['net_specific_angular_momentum'].sum().rename(columns={'sum':'net_specific_angular_momentum','redshift':'redshift','assembly':'assembly'})
	# net_specific_angular_momentum_redshift['expansion_factor'] = 1/(1+net_specific_angular_momentum_redshift['redshift'])
	# plot_expansion_factor_net_specific_angular_momentum(net_specific_angular_momentum_redshift,assembly_names,show=False,suffix='_cdm')
	# plot_expansion_factor_net_specific_angular_momentum_a3by2(net_specific_angular_momentum_redshift,assembly_names,show=False,suffix='_cdm')
	
	grouped_particles = computed_values_df.replace({'group':{'PartType0':'PartType0+PartType4','PartType4':'PartType0+PartType4'}})
	net_specific_angular_momentum_all_particles = grouped_particles.groupby(['redshift','assembly','group'],as_index=False)['net_specific_angular_momentum'].sum().rename(columns={'sum':'net_specific_angular_momentum','redshift':'redshift','assembly':'assembly'})
	add_expansion_factor_column(net_specific_angular_momentum_all_particles)
	computed_masked_values_df 		= computed_masked_values_df.replace({'group':{'PartType1':'Inner_Halo_PartType1'}})
	net_specific_angular_momentum_inner_dm_halo = computed_masked_values_df[['redshift','assembly','group','net_specific_angular_momentum','expansion_factor']]
	# plot_net_specific_angular_momentum_combined(net_specific_angular_momentum_all_particles,assembly_names,particles_titles,show=True,suffix='_only_L')
	plot_net_specific_angular_momentum_combined_with_inner_dm_halo(pd.concat([net_specific_angular_momentum_all_particles,net_specific_angular_momentum_inner_dm_halo]),assembly_names,particles_titles,show=False,suffix='')

	median_radius_all_particles = grouped_particles.groupby(['redshift','assembly','group'])['distance']
	median_radius_combined_particles_list = list()
	for key,item in median_radius_all_particles:
		grouped_df = median_radius_all_particles.get_group(key)
		if len(grouped_df)==2:
			median_radius_combined_particles_list.append([key[0],key[1],key[2],np.median(np.concatenate((grouped_df.iloc[0],grouped_df.iloc[1]),axis=None))])
		elif len(grouped_df)==1:
			median_radius_combined_particles_list.append([key[0],key[1],key[2],np.median(grouped_df.iloc[0])])
		else:
			raise IndexError
	median_radius_column_list = [
									'redshift',
									'assembly',
									'group',
									'median_radius'
								]
	median_radius_all_particles_subdf = pd.DataFrame(
													median_radius_combined_particles_list,
													columns = median_radius_column_list
													)
	median_radius_all_particles_df 		= pd.concat([median_radius_all_particles_subdf,computed_masked_values_df[['redshift','assembly','group','median_radius']]],ignore_index=True).convert_dtypes()
	add_expansion_factor_column(median_radius_all_particles_df)
	plot_median_radius(median_radius_all_particles_df,assembly_names,particles_titles,show=False,suffix='')