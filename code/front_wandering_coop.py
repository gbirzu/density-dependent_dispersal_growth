import glob
import numpy as np

directory = 'data/nmigration/'
N_list = [10000, 25000, 50000, 75000, 100000, 250000, 500000, 1000000]
#N_list = [10000]
gf_list = [0.01]
migr_list = [0.25]
B_list = [0.0, 0.5, 1.0, 1.4, 1.6, 1.8, 2.0, 2.4, 2.5, 2.8, 3.0, 3.2, 3.5, 3.8, 4.0, 4.2, 4.8, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
#B_list = [-0.9, 0.2]

def create_lists(N, gf, migr, B):
    subdir = 'N'+str(N)+'/'
    subdir_path = directory + subdir
    data_dir = 'N'+str(N)+'_gf'+str(gf)+'_migr'+str(migr)+'_B'+str(B)+'/'

    file_list = glob.glob(subdir_path+data_dir+'velocity/*.txt')
    return subdir, subdir_path, file_list

def analyze_fronts(file_list):
    full_data_array = []
    for file_name in file_list:
        data_array = np.genfromtxt(file_name, delimiter=',')
        full_data_array.append(data_array)
    full_data_array = np.array(full_data_array)

    mean = full_data_array.mean(axis=0)
    var = full_data_array.var(axis=0)
    analyzed_data = []

    for index in range(len(mean)):
        analyzed_data.append([mean[index][0], mean[index][1], var[index][1]])
    analyzed_data = np.array(analyzed_data)
    return analyzed_data

    
if __name__ == '__main__':
    
    full_data = []
    parameters = []
    for N in N_list:
        for gf in gf_list:
            for migr in migr_list:
                for B in B_list:
                    sub_list, sub_paths, file_list = create_lists(N, gf, migr, B)
                    print N, gf, migr, B
                    front_data = analyze_fronts(file_list)
                    print front_data
                    if len(front_data) > 1:
                        np.save('velocity_N'+str(N)+'_gf'+str(gf)+'_migr'+str(migr)+'_B'+str(B)+'_wandering_avg', front_data)
                        full_data.append(front_data)
                        #parameters.append([N, gf, migr, B])

    print full_data
    #full_data = np.array(full_data)
    #parameters = np.array(parameters)
    #np.save('velocity_data', full_data)
    #np.save('velocity_data_parameters', parameters)

