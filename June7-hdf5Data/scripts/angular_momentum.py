from plotter import *

# ===================================== User-defined function definitions ===============================

def get_numpy_arrays(df):
	coords 		= df[['coords_x','coords_y','coords_z']].to_numpy()*u.Mpc
	velocities	= df[['vel_x','vel_y','vel_z']].to_numpy()*(u.km/u.s)
	masses		= df['mass'].to_numpy()*u.M_sun
	return coords, velocities, masses

def get_rxv(coords, velocities):
	return np.cross(coords,velocities)

def get_specific_angular_momentum(rxv):
	return np.sum(rxv, axis=0)

def get_angular_momentum(rxv, mass):
	return np.sum(np.transpose(np.multiply(mass, np.transpose(rxv))),axis=0)

def get_norm(param3d):
	if param3d.ndim==1 :
		return np.sqrt(np.sum(np.square(param3d),axis=0))
	else :
		return np.sqrt(np.sum(np.square(param3d),axis=1))

def get_kinetic_energy(velocities,masses):
	return np.sum(0.5*masses*np.square(get_norm(velocities)),axis=0)

def get_plot_axes_titles():
	title  									= dict()
	title['redshift']						= 'Redshift'
	title['net_specific_angular_momentum']	= 'Total Specific Angular Momentum [$km$ $Mpc$ $s^{-1}$]'
	title['expansion_factor']				= 'Expansion Factor'
	title['net_angular_momentum']			= 'Total Angular Momentum [$km$ $Mpc$ $M_{\odot}$ $s^{-1}$]'
	title['total_kinetic_energy']			= 'Total Kinetic Energy [$km^{2}$ $M_{\odot}$ $s^{-2}$]'
	title['time'] 							= 'Time [$Gyr$]'
	title['masses']							= 'Mass [$M_{\odot}$]'
	return title

def plot_expansion_factor_net_specific_angular_momentum(df,assembly_names,show=True):
	prepare_plot()
	col = 'assembly'
	g 	= sns.relplot(
					data 	= computed_values_df, 
					x 		= 'expansion_factor',
					y  		= 'net_specific_angular_momentum',
					col  	= col,
					hue  	= col,
					kind  	= 'scatter',
					legend 	= False
					).set(
					xlabel 	= get_plot_axes_titles()['expansion_factor'],
					ylabel 	= get_plot_axes_titles()['net_specific_angular_momentum'])
	for row_val,ax in g.axes_dict.items():
		row_val  	= assembly_names[row_val]
		ax.set_title(capitalize_first_letter(col)+' = '+str(row_val))

	for axes in g.axes:
		for axis in axes:
			powerlaw_x = np.linspace(axis.get_xlim()[0],axis.get_xlim()[1],10)
			powerlaw_y = 20*powerlaw_x**(1.5)
			axis.plot(
							powerlaw_x,
							powerlaw_y,
							color='black',
							ls='--',
							lw=0.5
						)
			axis.loglog()

	plot_or_not(show=show,plot_name='expansion_factor-net_specific_angular_momentum')
	return

# ============================================ Main Program =============================================

if __name__ == '__main__' :

	assembly_list 			= ['gm_early','organic','gm_late']
	# assembly_list 				= ['organic']

	# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

	files  = {str(assembly):get_files(get_directory(str(assembly)+'_data')) for assembly in assembly_list}

	# # ------- Available fields for cols :
	# # ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass']

	# # ------- Get dataframe containing data for each type of assembly mode and store in a dictionary.

	df 		= {str(assembly):get_df(files[str(assembly)],cols) for assembly in assembly_list}

	# # ------- Close HDF5 file handles

	for assembly in files.keys():
		close_hdf5_files(files[assembly])

	# # ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	assembly_names  	= {'gm_early':'GM-Early','organic':'Organic','gm_late':'GM-Late'}
	for assembly in assembly_list:
		df[assembly].name  	= assembly_names[assembly]	

	computed_columns  		= [
								'specific_angular_momentum',
								'net_specific_angular_momentum',
								'angular_momentum',
								'net_angular_momentum',
								'total_kinetic_energy',
								'assembly'
							  ]

	computed_values_df  	= pd.DataFrame(columns=computed_columns)
	

	computed_values_assembly_list  	= list()
	units_flag  					= False
	units_dict 						= dict()
	for assembly in assembly_list:
		computed_values_dict  		= dict()
		for z in df[assembly]['redshift'].unique():
			subdf 	= df[assembly].groupby('redshift').get_group(z)
			coords, velocities, masses 	= get_numpy_arrays(subdf)
			rxv 						= get_rxv(coords,velocities)
			specific_angular_momentum	= get_specific_angular_momentum(rxv)
			net_specific_angular_momentum	= get_norm(specific_angular_momentum)
			angular_momentum  			= get_angular_momentum(rxv,masses)
			net_angular_momentum  		= get_norm(angular_momentum)
			total_kinetic_energy  		= get_kinetic_energy(velocities,masses)
			computed_values_dict[z]  	= [ 
											specific_angular_momentum.value,
											net_specific_angular_momentum.value,
											angular_momentum.value,
											net_angular_momentum.value,
											total_kinetic_energy.value,
											assembly
										  ]
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

		computed_values_assembly_df  						= pd.DataFrame.from_dict(
																					computed_values_dict,
																					orient='index',
																					columns=computed_columns
																					)
		computed_values_assembly_df['redshift'] 			= computed_values_assembly_df.index
		computed_values_assembly_df['expansion_factor'] 	= 1/(1+computed_values_assembly_df['redshift'])

		computed_values_assembly_list.append(computed_values_assembly_df.convert_dtypes())


	computed_values_df 		= pd.concat(computed_values_assembly_list,ignore_index=True)
	
	plot_expansion_factor_net_specific_angular_momentum(computed_values_df,assembly_names,show=False)

