import numpy as np
import pickle

data_path = '../data/profile_snapshots.dat'
output_dir = '../plots/theory/'


def average_replicas(data):
    het_global_averaged = np.zeros(data[1].shape[0])
    het_local_averaged = np.zeros(data[1].shape[0])
    replicas = 0
    for key in data.keys():
        het_global, het_local = calculate_heterozygosities(data[key])
        het_global_averaged += het_global
        het_local_averaged += het_local
        replicas += 1

    first_key = list(data.keys())[0]
    time, _, _ = separate_data(data[first_key])
    return time, het_global_averaged / replicas, het_local_averaged / replicas


def calculate_heterozygosities(profile):
    t_array, n1_array, n2_array = separate_data(profile)
    het_global = het_global_average(n1_array, n2_array)
    f_array = calculate_f(n1_array, n2_array)
    het_local = het_local_average(f_array, n1_array + n2_array)
    return het_global, het_local

def separate_data(profile):
    t_list = []
    n1_list = []
    n2_list = []
    for snapshot in profile:
        t_list.append(snapshot[0])
        n1_list.append([])
        n2_list.append([])

        for i in range(1, len(snapshot), 2):
            n1_list[-1].append(snapshot[i])
            n2_list[-1].append(snapshot[i + 1])
    return np.array(t_list), np.array(n1_list), np.array(n2_list)

def calculate_f(n1_array, n2_array):
    n_array = n1_array + n2_array
    f_array = np.zeros(n_array.shape)
    non_zero = np.where(n_array > 0)
    f_array[non_zero] = n1_array[non_zero] / n_array[non_zero]
    return f_array

def het_global_average(n1_array, n2_array):
    n1_total = np.sum(n1_array, axis=1)
    n_total = np.sum(n1_array + n2_array, axis=1)
    f_global = n1_total / n_total
    return 2 * f_global * (1 - f_global)

def het_local_average(f_array, n_array):
    f_average = []
    for i, n in enumerate(n_array):
        f_average.append(np.mean(f_array[i][np.where(n > 0)[0]]))
    f_average = np.array(f_average)
    return 2 * f_average * (1 - f_average)

if __name__ == '__main__':
    with open(data_path, 'rb') as f_in:
        data = pickle.load(f_in)
    calculate_heterozygosities(data[1])
    time, het_global_averaged, het_local_averaged = average_replicas(data)
    het_averaged = {'time':time, 'global':het_global_averaged, 'local':het_local_averaged}

    with open('../data/het_average.dat', 'wb') as f_out:
        pickle.dump(het_averaged, f_out)
