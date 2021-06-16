from plotter import *

# ===================================== User-defined function definitions ===============================

def get_numpy_arrays(df):
	coords 		= df[['coords_x','coords_y','coords_z']].to_numpy()
	velocities	= df[['vel_x','vel_y','vel_z']].to_numpy()
	masses		= df['mass'].to_numpy()
	return coords, velocities, masses

def get_rxv(coords, velocities):
	return np.cross(coords,velocities)

def get_net_angular_momentum(rxv, mass):
	return np.sum(np.transpose(np.multiply(mass, np.transpose(rxv))),axis=0)

# ============================================ Main Program =============================================

if __name__ == '__main__' :

	# ------- Get HDF5 file handles list for each type of assembly.

	gm_early_files		= get_files(get_directory('gm_early_data'))
	organic_files		= get_files(get_directory('organic_data'))
	gm_late_files		= get_files(get_directory('gm_late_data'))

	# ------- Available fields for cols :
	# ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	cols				= ['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass']

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

	# df_list 			= [gm_early_df,organic_df,gm_late_df]
	# df_list 			= [organic_df]

	gm_early_coords, gm_early_velocities, gm_early_masses 	= get_numpy_arrays(gm_early_df)			
	organic_coords,	organic_velocities, organic_masses		= get_numpy_arrays(organic_df)
	gm_late_coords,	gm_late_velocities, gm_late_masses		= get_numpy_arrays(gm_late_df)

	gm_early_rxv	= get_rxv(gm_early_coords,gm_early_velocities)
	organic_rxv		= get_rxv(organic_coords,organic_velocities)
	gm_late_rxv		= get_rxv(gm_late_coords,gm_late_velocities)

	gm_early_L  	= get_net_angular_momentum(gm_early_rxv,gm_early_masses)
	organic_L		= get_net_angular_momentum(organic_rxv,organic_masses)
	gm_late_L 		= get_net_angular_momentum(gm_late_rxv,gm_late_masses)

	# print(gm_early_L, organic_L, gm_late_L)