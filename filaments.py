#!/usr/bin/env python




# Assumes inputs are in decimal degrees.
# ra_cut should be a line of longitude at which there are on objects (so the sphere can be cut and unwound here).
def adjust_ra(ra, ra_cut):
	import numpy as np
	return np.where(ra <= ra_cut, ra, ra - 360.0)


# See p. 401
# Units for return value are are pc^2 h^-1 solar_mass^-1
def sigma_crit_inv(redshift_fg, redshift_bg):
	import numpy as np
	import cosmic_web_utilities as cwu

	if True:

		comoving_distance_fg = 1.0e6 * cwu.redshift_to_comoving_distance(redshift_fg) # pc/h
		comoving_distance_bg = 1.0e6 * cwu.redshift_to_comoving_distance(redshift_bg) # pc/h

		four_pi_G_over_c_squared = 6.01e-13 # pc / solar mass

		return (four_pi_G_over_c_squared * comoving_distance_fg * (comoving_distance_bg - comoving_distance_fg) * (1 + redshift_fg)) / comoving_distance_bg # pc^2 h^-1 solar_mass^-1
	else:
		return redshift_bg * 0.0 + 1.0


def filaments_core(configuration, part_num):

	import sys
	import numpy as np
	import pylab
	import os
	import astropy.io.fits as pyfits
	import glob
	import itertools
	import math
	import itertools
	import cPickle
	import time
	import cosmic_web_utilities as cwu
	import scipy.spatial as sp
	import logging
	import kmeans_radec


	base_directory = "/share/splinter/ucapwhi/cosmic_web/WL/filaments/"

	# See https://docs.python.org/2/howto/logging.html for information about logging.
	#logging.basicConfig(filename=base_directory + "logs/log_filaments.log", level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%Y%b%d %H:%M:%S')

	
	el = "====================="
	print el
	#logging.info(el)
	el = "Program start"
	print el
	#logging.info(el)

	# Constants
	RA = 0
	DEC = 1

	out_dict = {}

	config_use_pairs = (configuration == 0 or configuration == 2)

	if configuration == 0 or configuration == 1:
		
		data_description = "BCC"

		input_bg_data_file_name = base_directory + "input/Buzzard_v1.1_truth.subset_" + str(part_num).zfill(3) + ".fits"
		input_fg_data_file_name = base_directory + "input/Buzzard_foreground_objects.fit.gz"

		# The following items describe where on the sky the data is (and isn't). Values are in decimal degrees.
		# See my notes p. HD 281.
		ra_projection_origin = 10.0
		dec_projection_origin = -50.0
		ra_cut = 180.0 + ra_projection_origin 
		max_delta_dec_from_dec_projection_origin = 55.0
		gamma1_sign = 1.0
		gamma2_sign = 1.0

		fg_redshift_label = "ZSPEC"
		has_bg_NBC = False

	else:

		data_description = "Y1A1"

		input_bg_data_file_name = base_directory + "input/y1a1-im3shape-r-1-1-1.big.subset_" + str(part_num).zfill(3) + ".fits"
		input_fg_data_file_name = base_directory + "input/Y1A1_foreground_objects.fit.gz"

		# The following items describe where on the sky the data is (and isn't). Values are in decimal degrees.
		ra_projection_origin = 20.0
		dec_projection_origin = -50.0
		ra_cut = 180.0 + ra_projection_origin 
		max_delta_dec_from_dec_projection_origin = 55.0
		gamma1_sign = 1.0
		gamma2_sign = -1.0

		fg_redshift_label = "ZREDMAGIC"
		has_bg_NBC = True

	

	configuration_description = data_description + " " + ("Pairs" if config_use_pairs else "Singletons")		

	el = "Configuration = " + str(configuration)
	print el
	#logging.info(el)
	out_dict["configuration"] = configuration

	el = "Configuration description = " + configuration_description
	print el
	#logging.info(el)
	out_dict["configuration_description"] = configuration_description

	el = "Background data file = " + input_bg_data_file_name
	print el
	#logging.info(el)
	out_dict["input_bg_data_file_name"] = input_bg_data_file_name

	el = "Use pairs? = " + str(config_use_pairs)
	print el
	#logging.info(el)
	out_dict["config_use_pairs"] = config_use_pairs

	output_file_name = base_directory + "output/configuration" + str(configuration).zfill(2) + "." + str(part_num).zfill(3) + ".pic"

	el = "Output file = " + output_file_name
	print el
	#logging.info(el)
	out_dict["output_file_name"] = output_file_name

	# Delete any existing output file (see p. HD 382)
	if os.path.isfile(output_file_name):
		os.remove(output_file_name)
		el = "Deleted existing output file = " + output_file_name
		print el
		#logging.info(el)	
	

	# Background data

	x_bg = pyfits.open(input_bg_data_file_name)
	ra_bg = adjust_ra(x_bg[1].data.field("RA").astype(float,copy=False), ra_cut)
	dec_bg = x_bg[1].data.field("DEC").astype(float,copy=False)
	redshift_bg = x_bg[1].data.field("Z").astype(float,copy=False)
	gamma1_bg = gamma1_sign * x_bg[1].data.field("GAMMA1").astype(float,copy=False)
	gamma2_bg = gamma2_sign * x_bg[1].data.field("GAMMA2").astype(float,copy=False)
	if has_bg_NBC:
		info_flag_bg = x_bg[1].data.field("INFO_FLAG").astype(int,copy=False) # Note 'int'
		m_bg = x_bg[1].data.field("M").astype(float,copy=False)
		c1_bg = gamma1_sign * x_bg[1].data.field("C1").astype(float,copy=False)
		c2_bg = gamma2_sign * x_bg[1].data.field("C2").astype(float,copy=False)
		w_bg = x_bg[1].data.field("W").astype(float,copy=False)
	else:
		num_items = len(ra_bg)
		info_flag_bg = np.zeros(num_items)
		m_bg = np.zeros(num_items)
		c1_bg = np.zeros(num_items)
		c2_bg = np.zeros(num_items)
		w_bg = np.ones(num_items)

	x_bg.close()

	if False:
		# See p. HD 374

		# Now filter out bg objects with extreme shears. See p. HD 347 and p. HD 373.
		shear_filter = np.where((gamma1_bg**2 + gamma2_bg**2) < 0.64)
		ra_bg = ra_bg[shear_filter]
		dec_bg= dec_bg[shear_filter]
		redshift_bg = redshift_bg[shear_filter]
		gamma1_bg = gamma1_bg[shear_filter]
		gamma2_bg = gamma2_bg[shear_filter]

		print "<gamma1_bg> = " + str(np.average(gamma1_bg))
		print "<gamma2_bg> = " + str(np.average(gamma2_bg))

		gamma1_bg = gamma1_bg - np.average(gamma1_bg)
		gamma2_bg = gamma2_bg - np.average(gamma2_bg)

	if True:
		wh = np.where(info_flag_bg == 0)
		ra_bg = ra_bg[wh]
		dec_bg = dec_bg[wh]
		redshift_bg = redshift_bg[wh]
		gamma1_bg = gamma1_bg[wh]
		gamma2_bg = gamma2_bg[wh]
		info_flag_bg = info_flag_bg[wh]
		m_bg = m_bg[wh]
		c1_bg = c1_bg[wh]
		c2_bg = c2_bg[wh]
		w_bg = w_bg[wh]


	(x_proj_bg, y_proj_bg) = cwu.azimuthal_equidistant_projection(ra_projection_origin, dec_projection_origin, ra_bg, dec_bg)


	# Build 2D kdtree for background objects, using the coordinates from the azimuthal_equidistant_projection. Can take a long time for large datasets.
	#logging.info("Starting to build kdtree for background shapes")
	kd_tree_bg = sp.cKDTree(np.column_stack((x_proj_bg, y_proj_bg)))
	#logging.info("Finished building kdtree for background shapes")



	# Foreground data

	x_fg = pyfits.open(input_fg_data_file_name)
	unfiltered_fg_ra = adjust_ra(x_fg[1].data.field("RA").astype(float,copy=False), ra_cut)
	unfiltered_fg_dec = x_fg[1].data.field("DEC").astype(float,copy=False)
	unfiltered_fg_redshift = x_fg[1].data.field(fg_redshift_label).astype(float,copy=False)
	x_fg.close()

	min_fg_z = 0.35 # Was 0.45
	max_fg_z = 0.50

	(ra_fg, dec_fg, redshift_fg) = cwu.triple_filter(unfiltered_fg_ra, unfiltered_fg_dec, unfiltered_fg_redshift, (-361.0, 361.0, -91.0, 91.0, min_fg_z, max_fg_z))

	(X_fg, Y_fg, Z_fg) = cwu.spherical_to_cartesian(ra_fg, dec_fg, redshift_fg)

	(x_proj_fg, y_proj_fg) = cwu.azimuthal_equidistant_projection(ra_projection_origin, dec_projection_origin, ra_fg, dec_fg)

	# Build 3D kdtree for foreground objects. Metric here is comoving distance in Mpc/h.
	kd_tree_fg = sp.cKDTree(np.column_stack((X_fg, Y_fg, Z_fg)))

	



	# We calculate the next few values even when not doing pairs as it helps us get a 'proxy' angular scale to determine the amount of stretching we need to do in the singles case.
	max_physical_separation_fg_pair = 10.0 # Mpc/h
	# We look at pairs whose angular separation is in the range (angular_separation_fg_pair_min_factor * max_angular_separation_fg_pair, angular_separation_fg_pair_max_factor * max_angular_separation_fg_pair).
	max_angular_separation_fg_pair = np.degrees(max_physical_separation_fg_pair / cwu.redshift_to_comoving_distance(min_fg_z)) #Degrees
	print "max_angular_separation_fg_pair = " + str(max_angular_separation_fg_pair)
	angular_separation_fg_pair_min_factor = 0.7
	angular_separation_fg_pair_max_factor = 1.0 # NOTE: If you change this to be bigger than 1.0 then the argument to kd_tree_fg.query_pairs() may need to change.




	if config_use_pairs:
		# Search for suitable pairs.
		# Begin with a first cut:
		fg_items_to_analyze = np.array(list(kd_tree_fg.query_pairs(max_physical_separation_fg_pair))) # A 2D array


		# Then refine to ensure that the angular separation on the sky is as desired.
		angular_separation_fg_pair = np.degrees(cwu.angular_separation(ra_fg[fg_items_to_analyze[:,0]], dec_fg[fg_items_to_analyze[:,0]], ra_fg[fg_items_to_analyze[:,1]], dec_fg[fg_items_to_analyze[:,1]])) # Degrees
		# Get the pairs whose separation is within the desired limits.
		this_filter, = np.where(np.logical_and(angular_separation_fg_pair >= angular_separation_fg_pair_min_factor * max_angular_separation_fg_pair, angular_separation_fg_pair < angular_separation_fg_pair_max_factor * max_angular_separation_fg_pair))
		fg_items_to_analyze = fg_items_to_analyze[this_filter,:]

		angular_separation_fg_pair = angular_separation_fg_pair[this_filter]
		# The following code duplicates the calculation of p_0 later on.
		kmeans_input = np.column_stack((0.5*(ra_fg[fg_items_to_analyze[:, 0]] + ra_fg[fg_items_to_analyze[:, 1]]), 0.5*(dec_fg[fg_items_to_analyze[:, 0]] + dec_fg[fg_items_to_analyze[:, 1]])))
		

	else:
		fg_items_to_analyze = np.arange(len(ra_fg))
		# Just a fixed number, representing thre separation to have between an object and a virtual nearby object.
		angular_separation_fg_pair = np.ones(len(fg_items_to_analyze)) * np.sqrt((angular_separation_fg_pair_min_factor**2 + angular_separation_fg_pair_max_factor**2) * 0.5) * max_angular_separation_fg_pair # Degrees
		kmeans_input = np.column_stack((ra_fg, dec_fg))

	
	lower_limit = -2.0
	upper_limit = 2.0
	num_bins = 40
	out_dict["lower_limit"] = lower_limit
	out_dict["upper_limit"] = upper_limit
	out_dict["num_bins"] = num_bins

	gamma1_binned = np.zeros((num_bins, num_bins))
	gamma2_binned = np.zeros((num_bins, num_bins))
	denominator_binned = np.zeros((num_bins, num_bins))


	do_jackknife = True

	out_dict["do_jackknife"] = do_jackknife

	if do_jackknife:

		el = "Starting calculation of jackknife regions"
		print el
		#logging.info(el)

		num_jackknife_bins = 64 # Clampitt uses 134 - not sure why this number.
		out_dict["num_jackknife_bins"] = num_jackknife_bins

		maxiter = 100
		tol = 1.0e-4 # was 1.0e-5

		km = kmeans_radec.kmeans_sample(kmeans_input, num_jackknife_bins, maxiter=maxiter, tol=tol)

		if not km.converged:
			raise RuntimeError("k means did not converge")
		jackknife_indices = km.labels # List has same length as kmeans_input

		gamma1_binned_jackknife = np.zeros((num_jackknife_bins, num_bins, num_bins))
		gamma2_binned_jackknife = np.zeros((num_jackknife_bins, num_bins, num_bins))
		denominator_binned_jackknife = np.zeros((num_jackknife_bins, num_bins, num_bins))



	percentage_counter = 0
	loop_counter = 0
	num_fg_items_analysed = 0
	num_bg_items_analysed = 0

	fg_items_to_analyze_len = fg_items_to_analyze.shape[0]




	el = "Starting loop for item analysis"
	print el
	#logging.info(el)


	for loop_counter in range(fg_items_to_analyze_len):

		new_percentage_counter = 10.0 * int(10.0 * (loop_counter) / float(fg_items_to_analyze_len))

		if new_percentage_counter > percentage_counter:
			percentage_counter = new_percentage_counter
			el = str(percentage_counter)
			print el
			#logging.info(el)
		
		
			
		# Need to define p_0, p_1 and k (to be used in call to tangent_plane_transformation) as well as r (for kd_tree_bg.query_ball_point).
		if config_use_pairs:

			i = fg_items_to_analyze[loop_counter, 0]
			j = fg_items_to_analyze[loop_counter, 1]

			# Put one of the LRGs at (-0.5, 0.0) and the other at (0.5, 0.0).

			p_0 = ((ra_fg[i] + ra_fg[j]) * 0.5, (dec_fg[i] + dec_fg[j]) * 0.5) # This is (approximately) the midpoint between the two LRGs

			p_1 = (ra_fg[i], dec_fg[i])

			k = -0.5 if np.random.rand() < 0.5 else 0.5 # Assigned randomly in case 'order in file' is correlated to mass.

			geometry_factor = 1.12 # See my notes p. HD 275. Note 1.12 ~ sqrt(5)/2.

			redshift_this_fg = 0.5 * (redshift_fg[i] + redshift_fg[j])

		else:

			# Put the LRG at (0,0) and choose another nearby point to use as p_1.
			p_0 = (ra_fg[loop_counter], dec_fg[loop_counter])

			if False:
				phi = np.random.rand() * 2.0 * np.pi # Radians
			else:
				phi = np.pi / 2.0 # Radians

			theta = angular_separation_fg_pair[loop_counter] # Degrees

			p_1 = (p_0[RA] + theta * np.cos(phi) / np.cos(np.radians(p_0[DEC])), p_0[DEC] + theta * np.sin(phi))

			k = 1.0

			geometry_factor = 1.0

			redshift_this_fg = redshift_fg[loop_counter]


		# Find the nearby background galaxies
		r = cwu.azimuthal_equidistant_projection_maximal_distortion(dec_projection_origin, p_0[DEC]) * geometry_factor * angular_separation_fg_pair[loop_counter] # In degrees.
		p_0_proj = cwu.azimuthal_equidistant_projection(ra_projection_origin, dec_projection_origin, p_0[RA], p_0[DEC])
		bg_indices = kd_tree_bg.query_ball_point(p_0_proj, r)
		len_bg_indices = len(bg_indices)

		if len_bg_indices > 0:

			# Rotate and stretch to bring to a common coodinate system.
			(S, G) = cwu.tangent_plane_transformation(p_0, p_1, k, np.column_stack((ra_bg[bg_indices], dec_bg[bg_indices])), np.column_stack((gamma1_bg[bg_indices]-c1_bg[bg_indices], gamma2_bg[bg_indices]-c2_bg[bg_indices])))

			sigma_critical_inv = sigma_crit_inv(redshift_this_fg, redshift_bg[bg_indices])


			# Bin the results
			x_indices = cwu.value_to_index(S[:,0], lower_limit, upper_limit, num_bins)
			y_indices = cwu.value_to_index(S[:,1], lower_limit, upper_limit, num_bins)

			if False:
				gamma1_binned[x_indices, y_indices] += G[:,0] * sigma_critical_inv
				gamma2_binned[x_indices, y_indices] += G[:,1] * sigma_critical_inv
				denominator_binned[x_indices, y_indices] += sigma_critical_inv**2
			else:
				# See https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Y1a1-im3shape-r-1-1-1. Note that c1 and c2 are taken into account before rotating.
				gamma1_binned[x_indices, y_indices] += w_bg[bg_indices] * G[:,0]
				gamma2_binned[x_indices, y_indices] += w_bg[bg_indices] * G[:,1]
				denominator_binned[x_indices, y_indices] += w_bg[bg_indices] * (1.0 + m_bg[bg_indices])
			

			if do_jackknife:
				jackknife_index = jackknife_indices[loop_counter]
				if False:
					gamma1_binned_jackknife[jackknife_index, x_indices, y_indices] += G[:,0] * sigma_critical_inv
					gamma2_binned_jackknife[jackknife_index, x_indices, y_indices] += G[:,1] * sigma_critical_inv
					denominator_binned_jackknife[jackknife_index, x_indices, y_indices] += sigma_critical_inv**2
				else:
					gamma1_binned_jackknife[jackknife_index, x_indices, y_indices] += w_bg[bg_indices] * G[:,0]
					gamma2_binned_jackknife[jackknife_index, x_indices, y_indices] += w_bg[bg_indices] * G[:,1]
					denominator_binned_jackknife[jackknife_index, x_indices, y_indices] += w_bg[bg_indices] * (1.0 + m_bg[bg_indices])


			num_fg_items_analysed += 1
			num_bg_items_analysed += len_bg_indices


	el = "Finished loop for item analysis"
	print el
	#logging.info(el)

	el = "Num fg items analysed = " + str(num_fg_items_analysed)
	print el
	#logging.info(el)
	out_dict["num_fg_items_analysed"] = num_fg_items_analysed

	el = "Num bg items analysed = " + str(num_bg_items_analysed)
	print el
	#logging.info(el)
	out_dict["num_bg_items_analysed"] = num_bg_items_analysed

	out_dict["gamma1_binned"] = gamma1_binned
	out_dict["gamma2_binned"] = gamma2_binned
	out_dict["denominator_binned"] = denominator_binned


	if do_jackknife:
		out_dict["jackknife_gamma1_binned"] = gamma1_binned_jackknife
		out_dict["jackknife_gamma2_binned"] = gamma2_binned_jackknife
		out_dict["jackknife_denominator_binned"] = denominator_binned_jackknife

	with open(output_file_name, "wb") as output_file:
		cPickle.dump(out_dict, output_file)


	el = "Program end"
	print el
	#logging.info(el)
	el = "====================="
	print el
	#logging.info(el)
		

# MAIN STARTS HERE

try:
	import sys

	configuration = 0 if len(sys.argv) < 2 else int(sys.argv[1])
	part_num = 0 if len(sys.argv) < 3 else int(sys.argv[2])

	filaments_core(configuration, part_num)

except RuntimeError as e:
        errMsg = e.args[0]
        print errMsg

