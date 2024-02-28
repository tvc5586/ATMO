import multiprocessing
import os

import main
import boxmod as bm

def multiprocess(input_array):
        
    # Generate processing pools
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    
    # Run computation in pools
    signal = pool.map(main.main, input_array)
    
    # Update taskbar to show how many have finished
    yield signal
