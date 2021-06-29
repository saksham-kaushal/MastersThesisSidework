from plotter import *
from angular_momentum import *

pd.set_option('display.expand_frame_repr', False)

# ===================================== User-defined function definitions ===============================

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

	# particles  			= ['PartType0', 'PartType1', 'PartType4']
	# particles  			= ['PartType0', 'PartType4']
	particles  			= ['PartType1']

	# ------- Get dataframe containing data for each type of assembly mode and store in a dictionary.

	df 		= {str(assembly):get_df(files[str(assembly)],cols,particles) for assembly in assembly_list}

	# ------- Close HDF5 file handles

	for assembly in files.keys():
		close_hdf5_files(files[assembly])
	# ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	assembly_names  	= {'gm_early':'GM-Early','organic':'Organic','gm_late':'GM-Late'}
	for assembly in assembly_list:
		df[assembly].name  	= assembly_names[assembly]	

	computed_columns  		= [
								'specific_angular_momentum',
								'net_specific_angular_momentum',
								'angular_momentum',
								'net_angular_momentum',
								'total_kinetic_energy',
								'assembly',
								'group',
								'redshift'
							  ]

	computed_values_df  	= pd.DataFrame(columns=computed_columns)
	

	computed_values_assembly_list  	= list()
	units_flag  					= False
	units_dict 						= dict()
	zero_redshift_mass_specificL	= dict()
	for assembly in assembly_list:
		computed_values_dict  		= list()
		for z in df[assembly]['redshift'].unique():
			supdf 	= df[assembly].groupby('redshift').get_group(z)
			for particle_type in particles:
				subdf 	= supdf.groupby('group').get_group(particle_type)
				coords, velocities, masses 		= get_numpy_arrays(subdf)
				rxv 							= get_rxv(coords,velocities)
				angular_momentum  				= get_angular_momentum(rxv,masses)
				net_angular_momentum  			= get_norm(angular_momentum)
				specific_angular_momentum		= get_specific_angular_momentum(angular_momentum,masses)
				net_specific_angular_momentum	= get_norm(specific_angular_momentum)
				total_kinetic_energy  			= get_kinetic_energy(velocities,masses)
				computed_values_dict.append([ 
												specific_angular_momentum.value,
												net_specific_angular_momentum.value,
												angular_momentum.value,
												net_angular_momentum.value,
												total_kinetic_energy.value,
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
																			computed_values_dict,
																			# orient='index',
																			columns=computed_columns
																			)
		# computed_values_assembly_df['redshift'] 			= computed_values_assembly_df.index
		computed_values_assembly_df['expansion_factor'] 	= 1/(1+computed_values_assembly_df['redshift'])
		computed_values_assembly_list.append(computed_values_assembly_df.convert_dtypes())

	computed_values_df 		= pd.concat(computed_values_assembly_list,ignore_index=True)
	print(computed_values_df.groupby('group').groups.keys())

	net_angular_momentum_redshift = computed_values_df.groupby(['redshift','assembly'],as_index=False)['net_angular_momentum'].sum().rename(columns={'sum':'net_angular_momentum','redshift':'redshift','assembly':'assembly'})
	net_angular_momentum_redshift['expansion_factor'] = 1/(1+net_angular_momentum_redshift['redshift'])
	plot_expansion_factor_net_angular_momentum(net_angular_momentum_redshift,assembly_names,powerlaw_offset=1e9,show=False,suffix='_cdm')

	net_specific_angular_momentum_redshift = computed_values_df.groupby(['redshift','assembly'],as_index=False)['net_specific_angular_momentum'].sum().rename(columns={'sum':'net_specific_angular_momentum','redshift':'redshift','assembly':'assembly'})
	net_specific_angular_momentum_redshift['expansion_factor'] = 1/(1+net_specific_angular_momentum_redshift['redshift'])
	plot_expansion_factor_net_specific_angular_momentum(net_specific_angular_momentum_redshift,assembly_names,show=False,suffix='_cdm')
	plot_expansion_factor_net_specific_angular_momentum_a3by2(net_specific_angular_momentum_redshift,assembly_names,show=False,suffix='_cdm')
	