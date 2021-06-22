from plotter import *

# ===================================== User-defined function definitions ===============================

def get_numpy_arrays(df):
	coords 		= df[['coords_x','coords_y','coords_z']].to_numpy()*u.Mpc
	velocities	= df[['vel_x','vel_y','vel_z']].to_numpy()*(u.km/u.s)
	masses		= df['mass'].to_numpy()*u.M_sun
	return coords, velocities, masses

def get_rxv(coords, velocities):
	return np.cross(coords,velocities)

def get_angular_momentum(rxv, mass):
	return np.sum(np.transpose(np.multiply(mass, np.transpose(rxv))),axis=0)

def get_norm(param3d):
	if param3d.ndim==1 :
		return np.sqrt(np.sum(np.square(param3d),axis=0))
	else :
		return np.sqrt(np.sum(np.square(param3d),axis=1))

def get_kinetic_energy(velocities,masses):
	return np.sum(0.5*masses*np.square(get_norm(velocities)),axis=0)

# ============================================ Main Program =============================================

if __name__ == '__main__' :

	assembly_list 			= ['gm_early','organic','gm_late']
	# assembly_list 				= ['organic']

	# ------- Get HDF5 file handles list for each type of assembly.

	files  = {str(assembly):get_files(get_directory(str(assembly)+'_data')) for assembly in assembly_list}

	# gm_early_files		= get_files(get_directory('gm_early_data'))
	# organic_files		= get_files(get_directory('organic_data'))
	# gm_late_files		= get_files(get_directory('gm_late_data'))

	# # ------- Available fields for cols :
	# # ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass']

	# # ------- Get dataframe containing data for each type of assembly mode.

	df 		= {str(assembly):get_df(files[str(assembly)],cols) for assembly in assembly_list}

	# gm_early_df 		= get_df(gm_early_files,cols)
	# organic_df 			= get_df(organic_files,cols)
	# gm_late_df 			= get_df(gm_late_files,cols)

	# # ------- Close HDF5 file handles

	for assembly in files.keys():
		close_hdf5_files(files[assembly])

	# # ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	assembly_names  	= {'gm_early':'GM-Early','organic':'Organic','gm_late':'GM-Late'}
	for assembly in assembly_list:
		df[assembly].name  	= assembly_names[assembly]	
	# gm_early_df.name  	= 'GM-Early'
	# organic_df.name  	= 'Organic'
	# gm_late_df.name  	= 'GM-Late'

	# gm_early_L  		= dict()
	# organic_L  		= dict()
	# gm_late_L  		= dict()

	computed_columns  		= [
								'redshift',
								'specific_angular_momentum',
								'angular_momentum',
								'net_angular_momentum',
								'total_kinetic_energy',
								'assembly'
							  ]
	computed_values_df  	= pd.DataFrame(columns=computed_columns)
	

	computed_values_assembly_list  	= list()
	for assembly in assembly_list:
		computed_values_dict  		= dict()
		for z in df[assembly]['redshift'].unique():
			subdf 	= df[assembly].groupby('redshift').get_group(z)
			coords, velocities, masses 	= get_numpy_arrays(subdf)
			specific_angular_momentum	= get_rxv(coords,velocities)
			angular_momentum  			= get_angular_momentum(specific_angular_momentum,masses)
			net_angular_momentum  		= get_norm(angular_momentum)
			total_kinetic_energy  		= get_kinetic_energy(velocities,masses)
			computed_values_dict[z]  	= [ 
											specific_angular_momentum,
											angular_momentum,
											net_angular_momentum,
											total_kinetic_energy,
											assembly
										  ]
		computed_values_assembly_df  	= pd.DataFrame.from_dict(computed_values_dict,orient='index',columns=computed_columns[1:])
		computed_values_assembly_df['redshift'] 	= computed_values_assembly_df.index

		computed_values_assembly_list.append(computed_values_assembly_df)

	computed_values_df 		= pd.concat(computed_values_assembly_list)
	print(computed_values_df)



	# for z in sorted(gm_early_df['redshift'].unique(),reverse=True):

	# 	gm_early_subdf 	= gm_early_df.groupby('redshift').get_group(z)
	# 	organic_subdf 	= organic_df.groupby('redshift').get_group(z)
	# 	gm_late_subdf 	= gm_late_df.groupby('redshift').get_group(z)

	# 	gm_early_coords, gm_early_velocities, gm_early_masses 	= get_numpy_arrays(gm_early_subdf)			
	# 	organic_coords,	organic_velocities, organic_masses		= get_numpy_arrays(organic_subdf)
	# 	gm_late_coords,	gm_late_velocities, gm_late_masses		= get_numpy_arrays(gm_late_subdf)

	# 	gm_early_rxv	= get_rxv(gm_early_coords,gm_early_velocities)
	# 	organic_rxv		= get_rxv(organic_coords,organic_velocities)
	# 	gm_late_rxv		= get_rxv(gm_late_coords,gm_late_velocities)

	# 	gm_early_L[z]  	= get_net_angular_momentum(gm_early_rxv,gm_early_masses)
	# 	organic_L[z]	= get_net_angular_momentum(organic_rxv,organic_masses)
	# 	gm_late_L[z] 	= get_net_angular_momentum(gm_late_rxv,gm_late_masses)

		
	# 	# gm_early_KE 	= get_net_kinetic_energy(gm_early_velocities,gm_early_masses)
	# 	# organic_KE 		= get_net_kinetic_energy(organic_velocities,organic_masses)
	# 	# gm_late_KE 		= get_net_kinetic_energy(gm_late_velocities,gm_late_masses)

	# 	# print(z, gm_early_KE.si,organic_KE.si,gm_late_KE.si)
	# # fig 			= plt.figure()
	# # ax  			= fig.add_subplot(projection='3d')
	# # ax.quiver(0,0,0,100,120,14,length=0.05)
	# # ax.set_xlabel('x')
	# # ax.set_ylabel('y')
	# # ax.set_zlabel('z')
	# # plt.show()

	# # 	plt.scatter(z,get_norm(gm_late_L))
	# # 	plt.semilogy()
	# # plt.show()
	# gm_early_L_df  	= pd.DataFrame.from_dict(gm_early_L,orient='index')
	# organic_L_df  	= pd.DataFrame.from_dict(organic_L,orient='index')
	# gm_late_L_df  	= pd.DataFrame.from_dict(gm_late_L,orient='index')

	# gm_early_L_df.name 	= 'GM-Early'
	# organic_L_df.name 	= 'Organic'
	# gm_late_L_df.name 	= 'GM-Late'

	# L_df  			= pd.concat()
	# print(gm_early_L, organic_L, gm_late_L)
	# 	