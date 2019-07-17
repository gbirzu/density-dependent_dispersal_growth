import sys
import getopt
import os
import glob
import csv
import numpy as np

def strip(s): #separates numbers from letters in string
    head = s.strip("-.0123456789")
    tail = s[len(head):]
    return head, tail


def get_variables(name):
    name_root = name.split("/")[-1].split(".")#Get file name
    if "txt" in name_root:
        name_root.remove("txt") #remove ending
    name_root = ".".join(name_root) #reform name
    aux = [strip(s) for s in name_root.split("_")]
    #default values if none found
    r0 = 0.01
    m0 = 0.01
    A = 0.0
    B = 0.0
    N = 10000
    for s in aux:
        if s[0] == "m":
            m0 = float(s[1])
        elif s[0] == "A":
            A = float(s[1])
        elif s[0] == "r":
            r0 = float(s[1])
        elif s[0] == "B":
            B = float(s[1])
        elif s[0] == "N":
            N = int(s[1])
    return m0, A, r0, B, N

def get_wandering(file_list):
    full_data_array = []
    for file_name in file_list:
        if os.path.getsize(file_name) > 0:
            data_array = np.genfromtxt(file_name, delimiter=",")
            if data_array.shape == ( 1000, 2 ):
                # Only add if correct format
                full_data_array.append(data_array)
    full_data_array = np.array(full_data_array)

    mean = full_data_array.mean(axis=0)
    var = full_data_array.var(axis=0)
    analyzed_data = []

    for index in range(len(mean)):
        analyzed_data.append([mean[index][0], mean[index][1], var[index][1]])
    analyzed_data = np.array(analyzed_data)
    return analyzed_data


def save_averages( files_arr, file_type, filename_tail):
    if file_type in ( "h", "heterozygosity" ):
        prelim_avg = np.sum(files_arr, axis=0)/float(len(files_arr))
        files_avg = []
        files_avg.append(prelim_avg.T[0])#add generations
        files_avg.append(prelim_avg.T[1])#add heterozygosity
        files_avg.append(len(files_arr)*prelim_avg.T[2])#add survival
        files_avg = np.asarray(prelim_avg).T
        file_name = "hetero_" + filename_tail

    elif file_type in ( "v", "velocity" ):
        files_avg = np.sum(files_arr, axis=0)/float(len(files_arr))
        file_name = "velocity_" + filename_tail

    elif file_type in ( "p", "profile" ):
        files_avg = np.sum(files_arr, axis=0)/float(len(files_arr))
        #wandering_data = get_wandering(files_list)
        file_name = "profile_" + filename_tail

    np.save(file_name, files_avg)


def average_subdirectory( subdir_path, params, number_check, extension, file_type="h" ):
    [N, A, B, r0, m0] = params
    parameter_string = "N" + str(N) + "_r" + str(r0) + "_m" + str(m0) + "_A" + str(A) + "_B" + str(B)
    print "Averaging ", parameter_string

    if file_type in ( "h", "heterozygosity" ):
        subsubdir = subdir_path + "/hetero/hetero_" + parameter_string + "*.txt"
    elif file_type in ( "v", "velocity" ):
        subsubdir = subdir_path + "/velocity/velocity_" + parameter_string + "*.txt"
    elif file_type in ( "p", "profile" ):
        subsubdir = subdir_path + "/profile/profile_" + parameter_string + "*.txt"

    files_list = glob.glob( subsubdir )
    if len(files_list) == 0:
        print "No files found in ", subsubdir
        print "Aborting"
        return 1

    if len( files_list ) > 0:
        print "Number of files in ", subsubdir, ": ", len(files_list)
        #print subdir #print directory to be read

        counter = 0
        std_length = 0
        files_arr = []
        while counter < len(files_list):
            lst_file = list(csv.reader(open(files_list[counter]))) #load file
            if len(lst_file) == 0 and counter < len(files_list) - 1: #check if file empty and not last file
                counter += 1
                lst_file = list(csv.reader(open(files_list[counter]))) #load next file
                continue
            elif len(lst_file) == 0: #if last file, break
                break

            my_file = files_list[counter]

            aux_arr = [] #array with all data
            if file_type in ( "h", "heterozygosity" ):
                for e in lst_file:
                    if len(e) == 3 and float(e[1]) > 0:
                        aux_arr.append([int(e[0]), float(e[1]), 1])
                    elif len(e) == 3 and float(e[1]) == 0:
                        aux_arr.append([int(e[0]), float(e[1]), 0])
                    else:
                        print e
                        break
            elif file_type in ( "v", "velocity" ):
                for e in lst_file:
                    if len(e) == 2: #check if proper format
                        aux_arr.append([int(e[0]), float(e[1])])
                    else:
                        print e
                        break #if not break
            elif file_type in ( "p", "profile" ):
                for e in lst_file:
                    if len(e) == 3: #check if proper format
                        x = int(e[0])
                        #fractinal data
                        p1 = float(e[1])
                        p2 = float(e[2])
                        h = 2*p1*p2
                        aux_arr.append([x, p1, p2, h])
                    else:
                        print e
                        break #if not break
            if counter == 0:
                        std_length = len(aux_arr)
            if len(aux_arr) == std_length:
                files_arr.append(aux_arr)
            counter += 1

        if len(files_arr) >= 0.8*number_check: #check if correct number of files
            filename_tail = "N" + str(N) + "_r" + str(r0) + "_m" + str(m0) + "_A" + str(A) + "_B" + str(B) + "_avg" + str(extension)
            save_averages( files_arr, file_type, filename_tail )
        else:
            print "ERROR: Wrong number of files: ", subdir_path, len( files_arr )
    return 0


def analyze_subdir( subdir_path, extension, r_check=None, m_check=None ):
    result_dict = { 0:"Yes", 1:"No" }
    m0, A, r0, B, N = get_variables( subdir_path.split( "/" )[ -1 ] )
    print "Data variables: ", N, A, B, r0, m0

    if r_check != None:
        r0 = r_check
    if m_check != None:
        m0 = m_check

    result = average_subdirectory( subdir_path, [ N, A, B, r0, m0 ], 1000, extension, "h" )
    print "Heterozygosity success: ", result_dict[ result ]
    result = average_subdirectory( subdir_path, [ N, A, B, r0, m0 ], 1000, extension, "p" )
    print "Profile: ", result_dict[ result ]
    result = average_subdirectory( subdir_path, [ N, A, B, r0, m0 ], 1000, extension, "v" )
    print "Velocity success: ", result_dict[ result ]
    print "\n"


def analyze_Ndirectory( Ndir_path, extension, r_check=None, m_check=None ):
    subdir_list = sorted( glob.glob( Ndir_path + "/N*" ) )

    for subdir in subdir_list:
        analyze_subdir( subdir, extension, r_check, m_check )


def analyze_all( path, extension, r_check=None, m_check=None ):
    # Averages heterozygosities, final profiles, and total population for all N# dirs and subdirs in path
    dir_list = sorted( glob.glob( path + "/N*" ) )
    for dir in dir_list:
        analyze_Ndirectory( dir, extension, r_check, m_check )


def main(argv):
    # Set default parameters
    path = os.getcwd()
    extension = ""
    # By default average only heterozygosity files
    het_flag = 1
    vel_flag = 0
    prof_flag = 0
    number_check = 1
    N = 10000
    r0 = None
    m0 = None
    analyze_all = False
    subdir = None


    try:
        opts, args = getopt.getopt(argv, "d:h:v:p:n:e:N:r:m:a:s:", ["dir=", "het=", "vel=", "prof=", "ncheck=", "ext=", "all=", "subdir="])
    except getopt.GetoptError:
        print "Input error! Using default values."
    for opt, arg in opts:
        if opt in ("-d", "--dir"):
            path = path + "/" + arg
        elif opt in ("-h", "--het"):
            het_flag = int(arg)
        elif opt in ("-v", "--vel"):
            vel_flag = int(arg)
        elif opt in ("-p", "--prof"):
            prof_flag = int(arg)
        elif opt in ("-n", "--ncheck"):
            number_check = int(arg)
        elif opt in ("-e", "--ext"):
            extension = arg
        elif opt in ("-N"):
            N = int( arg )
        elif opt in ("-r"):
            r0 = float(arg)
        elif opt in ("-m"):
            m0 = float(arg)
        elif opt in ("-a", "--all"):
            analyze_all = True
        elif opt in ("-s", "--subdir"):
            subdir = arg

    if analyze_all == False:
        if subdir == None:
            N_directory = path + "/N" + str( N )
            analyze_Ndirectory( N_directory, extension, r0, m0 )
        else:
            subdir_path = path + "/N" + str( N ) + "/" +  subdir
            analyze_subdir( subdir_path, extension, r0, m0 )

    else:
        analyze_all( path, extension )

if __name__ == "__main__":
    main(sys.argv[1:])
