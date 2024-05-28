import numpy as np
import os, sys
import h5py

# custom lib
curr_path = os.getcwd()
sys.path.append(curr_path+'/../')
from tools.ara_utility import size_checker
from tools.ara_sim_load import ara_raytrace_loader

def arr_time_table_loader(Station = None):
    """
    Generate arrival time tables and save to h5 file in: 
      `$OUTPUT_PATH/ARA0{Station}/arr_time_table/`
      and named: `arr_time_table_A{Station}_original_3radii.h5`

    Parameters
    ----------
    Station : int
        The station number. Not sure why it's initialized to `None` -ARB

    Output File Arrays
    ------------------
    ant_pos
        Antenna positions [m] in station-centric coordinates and ordered 
          according to RF Channel indexing (so `ant%4` gives the string number).
        For more information, refer to [1]_.
        Shape : (3, 16)
    theta_bin
        Event zenith angle range [deg] for which the table is tabulated.
        Defined in `tools.ara_sim_load.get_src_trg_position()`.
        Shape : (n_thetas)
    phi_bin
        Event azimuthal angle range [deg] for which the table is tabulated.
        Defined in `tools.ara_sim_load.get_src_trg_position()`.
        Shape : (n_phis)
    radius_bin
        Event radial distance range [m] for which the table is tabulated.
        Defined in `tools.ara_sim_load.get_src_trg_position()`.
        Shape : (n_radii)
    num_ray_sol
        The number of ray solutions for which the table is tabulated.
        Defined in `tools.ara_sim_load.get_src_trg_position()`.
        Default : [1,2], presumably meaning direct and refracte/reflected rays.
        Shape : (n_rays)
    path_len
        Path length [m] from source to target calculated for each theta, phi, 
          radius, antenna, and ray solution defined in `theta_bin`, `phi_bin`,
          `radius_bin`, `ant_pos`, and `num_ray_sol` respectively.
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    arr_time_table
        Travel time [ns] from source to target calculated for each theta, phi, 
          radius, antenna, and ray solution defined in `theta_bin`, `phi_bin`,
          `radius_bin`, `ant_pos`, and `num_ray_sol` respectively.
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    launch_ang
        Launch angle of the signal [rad] from the source calculated for each 
          theta, phi, radius, antenna, and ray solution defined in `theta_bin`, 
          `phi_bin`, `radius_bin`, `ant_pos`, and `num_ray_sol` respectively.
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    receipt_ang
        Receiving angle [rad] of the signal at the target calculated for each   
          theta, phi,radius, antenna, and ray solution defined in `theta_bin`, 
          `phi_bin`,`radius_bin`, `ant_pos`, and `num_ray_sol` respectively.
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    reflection_ang
        Reflection angle [rad] of the signal at the surface of the ice. 
        Calculated for each theta, phi, radius, antenna, and ray solution 
          defined in `theta_bin`, `phi_bin`, `radius_bin`, `ant_pos`, and 
          `num_ray_sol` respectively.
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    miss
        Distance between the target position of the ray solver and the final 
          position calculated for the ray [m]. Should not be large. 
        Calculated for each theta, phi, radius, antenna, and ray solution 
          defined in `theta_bin`, `phi_bin`, `radius_bin`, `ant_pos`, and 
          `num_ray_sol` respectively.
        From `AraSim/log.txt`:
          "miss distance is the distance between final result and target location"
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    attenuation
        Signal attenuation from source to target calculated for each theta, phi, 
          radius, antenna, and ray solution defined in `theta_bin`, `phi_bin`,
          `radius_bin`, `ant_pos`, and `num_ray_sol` respectively.
        Calculated in `tools.ara_sim_load.get_ray_solution()`.
        Shape : ( n_theta, n_phi, n_radii, n_ants, n_rays )
    ice_model
        Parameters of the ice model used for time table calculation. 
        Defined in `tools.ara_sim_load.__init__()`.
          
    References
    ----------
    [1] : https://aradocs.wipac.wisc.edu/0021/002188/001/bootcamp2020-araroot.pdf
    """

    print('Collecting arrival time starts!')

    # arasim raytracer
    #ara_ray = ara_raytrace_loader(n0 = 1.35, nf = 1.78, l = 0.0132, verbose = True)
    ara_ray = ara_raytrace_loader(n0 = 1.326, nf = 1.78, l = 0.0202, verbose = True) #values for PA ice model

    # get posisiton for vertex and antenna
    radius_bin = np.array([41, 170, 300, 450, 600])
    ara_ray.get_src_trg_position(Station, 2016, radius_bin = radius_bin)

    # arrival time table
    ( path_len, arr_time_table, launch_ang, receipt_ang, reflection_ang, 
      miss, attenuation ) = ara_ray.get_arrival_time_table() 
 
    # create output dir
    Output = os.path.expandvars("$OUTPUT_PATH") + f'/ARA0{Station}/arr_time_table/'
    print(f'Output path check:{Output}')
    if not os.path.exists(Output):
        os.makedirs(Output)
    h5_file_name = f'{Output}arr_time_table_A{Station}_original_3radii.h5'#_all.h5'
    hf = h5py.File(h5_file_name, 'w')
    
    #saving result
    hf.create_dataset('ant_pos', data=ara_ray.ant_pos, compression="gzip", compression_opts=9)
    hf.create_dataset('theta_bin', data=ara_ray.theta_bin, compression="gzip", compression_opts=9)
    hf.create_dataset('phi_bin', data=ara_ray.phi_bin, compression="gzip", compression_opts=9)
    hf.create_dataset('radius_bin', data=ara_ray.radius_bin, compression="gzip", compression_opts=9)
    hf.create_dataset('num_ray_sol', data=ara_ray.num_ray_sol, compression="gzip", compression_opts=9)
    hf.create_dataset('path_len', data=path_len, compression="gzip", compression_opts=9)
    hf.create_dataset('arr_time_table', data=arr_time_table, compression="gzip", compression_opts=9)
    hf.create_dataset('launch_ang', data=launch_ang, compression="gzip", compression_opts=9)
    hf.create_dataset('receipt_ang', data=receipt_ang, compression="gzip", compression_opts=9)
    hf.create_dataset('reflection_ang', data=reflection_ang, compression="gzip", compression_opts=9)
    hf.create_dataset('miss', data=miss, compression="gzip", compression_opts=9)
    hf.create_dataset('attenuation', data=attenuation, compression="gzip", compression_opts=9)
    hf.create_dataset('ice_model', data=ara_ray.ice_model, compression="gzip", compression_opts=9)
    hf.close()
    print(f'output is {h5_file_name}', size_checker(h5_file_name))

if __name__ == "__main__":

    if len (sys.argv) < 2:
        Usage = """

    If it is data,
    Usage = python3 %s

    <Station ex)3>

        """ %(sys.argv[0])
        print(Usage)
        sys.exit(1)

    # argv
    station=int(sys.argv[1])

    arr_time_table_loader(Station = station)













    
