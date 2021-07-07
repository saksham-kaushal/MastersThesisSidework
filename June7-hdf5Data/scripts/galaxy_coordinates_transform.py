from plotter import *
import angular_momentum as am  
import coordinates_plotter as cp

# ===================================== User-defined function definitions ===============================

def get_orthonormal_basis(angular_momentum_vector):
	k = angular_momentum_vector

	try:
		nonzero_k = np.nonzero(k)[0][0]
	except IndexError:
		raise IndexError('Received a null angular momentum vector')
	if nonzero_k == 0:
		i = np.array([-k[1]/k[0],1,0])
	else:
		i = np.array([1,0,0])
	
	j = np.cross(k,i)
	basis_vectors = np.vstack((i,j,k))
	return normalize_3dvector(basis_vectors)

def coordinates_transform(galaxy_coords,unit_3dvectors):
	return np.dot(galaxy_coords,unit_3dvectors)
 
# ============================================ Main Program =============================================

if __name__ == '__main__':
	
	assembly_list 			= ['gm_early','organic','gm_late']
	# assembly_list 				= ['organic']

	# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

	files  = {str(assembly):get_files(get_directory(str(assembly)+'_data')) for assembly in assembly_list}

	# ------- Available fields for cols :
	# ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass']

	# ------- Get dataframe containing data for each type of assembly mode and store in a dictionary.

	df 		= {str(assembly):get_df(files[str(assembly)],cols) for assembly in assembly_list}

	# ------- Close HDF5 file handles

	for assembly in files.keys():
		close_hdf5_files(files[assembly])

	# ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	assembly_names  	= {'gm_early':'GM-Early','organic':'Organic','gm_late':'GM-Late'}
	for assembly in assembly_list:
		df[assembly].name  	= assembly_names[assembly]	

	computed_columns  		= [
								'coords',
								'angular_momentum',
								'net_angular_momentum',
								'assembly'
							  ]

	computed_values_df  	= pd.DataFrame(columns=computed_columns)
	

	computed_values_assembly_list  	= list()
	units_flag  					= False
	units_dict 						= dict()
	zero_redshift_mass_specificL	= dict()
	for assembly in assembly_list:
		computed_values_dict  		= dict()
		for z in df[assembly]['redshift'].unique():
			subdf 	= df[assembly].groupby('redshift').get_group(z)
			coords, velocities, masses 	= am.get_numpy_arrays(subdf)
			rxv 						= am.get_rxv(coords,velocities)
			angular_momentum  			= am.get_angular_momentum(rxv,masses)
			net_angular_momentum  		= get_norm(angular_momentum)
			# specific_angular_momentum	= get_specific_angular_momentum(angular_momentum,masses)
			# net_specific_angular_momentum	= get_norm(specific_angular_momentum)
			# total_kinetic_energy  		= get_kinetic_energy(velocities,masses)
			computed_values_dict[z]  	= [ 
											# specific_angular_momentum.value,
											# net_specific_angular_momentum.value,
											coords.value,
											angular_momentum.value,
											net_angular_momentum.value,
											# total_kinetic_energy.value,
											assembly
										  ]
			if not units_flag:
				units_dict['coords']  						= coords.unit
				units_dict['velocities']	  				= velocities.unit
				units_dict['masses']  						= masses.unit
				# units_dict['specific_angular_momentum']		= specific_angular_momentum.unit
				# units_dict['net_specific_angular_momentum']	= net_specific_angular_momentum.unit
				units_dict['angular_momentum']  			= angular_momentum.unit
				units_dict['net_angular_momentum']  		= net_angular_momentum.unit
				# units_dict['total_kinetic_energy']  		= total_kinetic_energy.unit
				
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
	computed_values_df['orthonormal_basis'] = computed_values_df['angular_momentum'].apply(get_orthonormal_basis)
	computed_values_df['modified_coords'] 	= computed_values_df.apply(lambda x:coordinates_transform(x['coords'],x['orthonormal_basis']),axis=1)
	print_df(computed_values_df[['coords']].iloc[0])
	print_df(computed_values_df[['orthonormal_basis']].iloc[0])
	print_df(computed_values_df[['modified_coords']].iloc[0])
	# print(computed_values_df['coords'].apply(lambda x:x[:,0]))
	# print(computed_values_df['coords'][0])
	n = [len(item) for item in computed_values_df['modified_coords']]
	transformed_coords_df  					= pd.DataFrame({'coords_x':np.concatenate(computed_values_df['modified_coords'].apply(lambda x:x[:,0]).values),
															'coords_y':np.concatenate(computed_values_df['modified_coords'].apply(lambda x:x[:,1]).values),
															'coords_z':np.concatenate(computed_values_df['modified_coords'].apply(lambda x:x[:,2]).values),
															'assembly':np.repeat(computed_values_df['assembly'],n),
															'redshift':np.repeat(computed_values_df['redshift'],n)}
															)

	df_list = list()
	for assembly in assembly_list:
		subdf = transformed_coords_df.groupby('assembly').get_group(assembly)
		subdf.name = assembly_names[assembly]
		df_list.append(subdf.convert_dtypes())

	for df in df_list:
		df.name = assembly_names[df['assembly'].unique()[0]]

	# cp.plot_galaxy_coordinates_combined(df_list,show=False)
	# cp.plot_galaxy_coordinates_individual(df_list,show=False)
	# cp.plot_galaxy_coordinates_3d(df_list,show=False)

	gif_dirs 			= [
							'GM-Early',
							'Organic',
							'GM-Late',
							os.path.join('3d','GM-Early'),
							os.path.join('3d','Organic'),
							os.path.join('3d','GM-Late')
						  ]

	dir_list			= [os.path.join(get_directory('plotter'),subdir) for subdir in gif_dirs]
	
	# make_gif(dir_list,fps=3)