
# import os
# import numpy as np
  
# def get_all_root_files(directory: "path to directory"):
  
#   # Get all files in the directory
#   files = os.listdir(directory)
  
#   # Filter for .root files
#   root_files = [f for f in files if f.endswith('.root')]
  
#   # Load the data from the root files into a numpy array
#   data = np.concatenate([np.load(os.path.join(directory, f)) for f in root_files])
  
#   return data