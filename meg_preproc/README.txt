The codes should be run in the following order
1. fetch_data.py ( here set the path where you want the fetched data to be saved)
2. anatomy.py (this will segment the mris and save the results. this requires the path that you set in fetch_data.py)
3. run_meg_preprocessing_v2.py (this will preprocess the meg data, make head models etc)
4. save_mesh_info_to_matlab.py (this will save all the mesh info from python data format to matlab format for further use)