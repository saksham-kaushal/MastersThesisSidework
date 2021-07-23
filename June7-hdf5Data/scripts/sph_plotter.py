from plotter import *
from sphviewer.tools import QuickView

# ===================================== User-defined function definitions ===============================

def get_numpy_array(df,col='coords'):
	if col=='coords':
		return df[['coords_x','coords_y','coords_z']].to_numpy()*u.Mpc
	elif col=='velocities':
		return df[['vel_x','vel_y','vel_z']].to_numpy()*(u.km/u.s)
	elif col=='masses':
		return df['mass'].to_numpy()*u.M_sun
	else:
		raise ValueError('Valid values for column are: 6coords, velocities and masses')

# ============================================ Main Program =============================================

if __name__ == '__main__':
	# assembly_list 			= ['gm_early','organic','gm_late']
	assembly_list 				= ['organic']

	# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

	files  = {str(assembly):get_files(get_directory(str(assembly)+'_baryon_cdm_data')) for assembly in assembly_list}

	# ------- Available fields for cols :
	# ------- 'coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass', 'redshift'
	
	# cols				= ['coords_x', 'coords_y', 'coords_z', 'vel_x', 'vel_y', 'vel_z', 'mass']

	# ------- Particle type 0 = Gas, 1 = Dark Matter, and 4 = Stars

	# particles  			= ['PartType0', 'PartType1', 'PartType4']
	# particles  			= ['PartType0', 'PartType4']
	particles  			= 'PartType1'
	particles_titles  	= {
							'PartType0':'Gas', 
							'PartType1':'Dark Matter Halo', 
							'PartType4':'Stars', 
							'PartType0+PartType4':'Galaxy (Baryonic Matter)',
							'Inner_Halo_PartType1':'Inner Dark Matter Halo'
							}

	# ------- Get dataframe containing data for each type of assembly mode and store in a dictionary.

	# df 		= {str(assembly):get_df(files[str(assembly)],cols,particles) for assembly in assembly_list}

	
	# ------- Provide plot-friendly names to dataframes for deffrent assembly modes.
	
	assembly_names  		= {
								'gm_early':'GM-Early',
								'organic':'Organic',
								'gm_late':'GM-Late'
							  }

	# for assembly in assembly_list:
	# 	df[assembly].name  	= assembly_names[assembly]	

	for assembly in assembly_list:
		for fname in files[assembly]:
			qv = QuickView(np.array(fname[particles+'/Coordinates']),
						r='infinity',
						plot=False
						)
			qv.imshow()

	# ------- Close HDF5 file handles

	for assembly in files.keys():
		close_hdf5_files(files[assembly])