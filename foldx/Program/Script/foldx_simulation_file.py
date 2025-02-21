# Importing Libraries
import sys, time, os
import psutil
import argparse
import shutil
import time
import ray
import subprocess
import glob
import gc

# Initialize ray
ray.init()


@ray.remote # generator
def run_foldx(commands):
	input_file = commands[3:6]
	output_file = commands[7:]

	print(input_file)
	print(output_file)

	foldx_path = input_file[0] # input_file[0] is actually directory NOT a file
	foldx_rotamer_file=foldx_path + '/rotabase.txt'

	shutil.copy(foldx_rotamer_file, './')

	pdb_info = input_file[1] # protein structure in pdb format
	pdb_name = pdb_info.split('/')[-1]
	pdb_dir = ('/').join(pdb_info.split('/')[0:-1])

	mutant_file = input_file[2] # file containing mutations
	
	output_dir = output_file[0] # output_file[0] is actually directory not a file
	output_dir_phenotype = output_file[1]
	
	pdb_name_new = pdb_name.split('.')[0]+'_Repair.pdb'

	if os.path.exists(output_dir + '/' + pdb_name_new):
		print('Repair file already present')
	else:
		all_info = foldx_path + '/foldx --command=RepairPDB --pdb-dir=' + pdb_dir + ' --pdb=' + pdb_name + ' --pH=7 --vdwDesign=2 --pdbHydrogens=false --output-dir=' + output_dir
		os.system(all_info) # Repairing the pdb file

	all_info2 = foldx_path + '/foldx --command=BuildModel --pdb-dir=' + output_dir + ' --pdb=' + pdb_name_new + ' --mutant-file=' + mutant_file + ' --pH=7 --vdwDesign=2 --output-dir=' + output_dir_phenotype + ' --pdbHydrogens=false --numberOfRuns=3'

	os.system(all_info2) # All the mutated pdb files stored in here

	#os.unlink('rotabase.txt')






if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="A script to parallelize foldx simulations.")
	parser.add_argument('-i', nargs='+', dest="input_argument",default="something.txt", help='help file here')
	parser.add_argument('-o', nargs='+', dest="output_argument",default="something_.txt",help='help file here')
	args_ = parser.parse_args()
	wd = os.getcwd() # working directory must be the directory that contains the folders 'Data' and 'Program'
	mutant_file = args_.input_argument[2] # path to the complete individual list of mutations
	cmd = 'mkdir '+ '/'.join(mutant_file.split('/')[:-1])+'/split_files' # making a directory to store the split files
	split_lines = int(input("Enter the number of lines for splitting: "))  # Convert input to integer
	os.system(cmd)
	cmd = f'split -l {split_lines} -d ' + '--additional-suffix=.txt ' + mutant_file + ' ' + '/'.join(mutant_file.split('/')[:-1])+'/split_files'+'/'+mutant_file.split('/')[-1].split('.')[0] + '_'
	os.system(cmd)
	cmd = f'split -l {split_lines} -d ' + '--additional-suffix=.txt ' + mutant_file + ' ' + '/'.join(mutant_file.split('/')[:-1])+'/'+mutant_file.split('/')[-1].split('.')[0] + '_'
	os.system(cmd)
	num_of_folders = len(os.listdir(os.chdir('/'.join(mutant_file.split('/')[:-1])))) # obtaining the number of directories to create for backtracking
	foldx_folders = args_.output_argument[1] # output directory that contains the repair structures
	os.chdir(wd)
	os.chdir(foldx_folders)
	if os.path.exists('../temp'):
		print('The directory already exists')
	else:
		os.mkdir('../temp') # creating a temp directory that stores the foldx dataframe file generated after analysis
	os.chdir(wd)
	for i in range(0,num_of_folders-2):
		print(args_.output_argument[1]+'/'+str(i).zfill(2))
		if os.path.exists(args_.output_argument[1]+'/'+str(i).zfill(2)):
			print('The directory already exists')
		else:
			os.chdir(foldx_folders)
			os.mkdir(str(i).zfill(2))
			os.chdir(wd) # create seperate directories for backtracking
	os.chdir(wd)
	with open('/'.join(foldx_folders.split('/')[:-1])+'/temp/'+'foldx_commands_file.txt','w') as f: # creating a file to keep tack of the commands used in the foldx simulation
		for i in range(0,num_of_folders-2):
			f.write('python3 Program/Script/foldx_simulation_file.py -i '+ args_.input_argument[0]+
		   ' '+args_.input_argument[1]+' '+args_.input_argument[2].split('.')[0]+'_'+str(i).zfill(2)+'.txt'
		   ' -o '+ args_.output_argument[0]+' '+args_.output_argument[1]+'/'+str(i).zfill(2)+'\n')
	all_commands = []
	with open('/'.join(foldx_folders.split('/')[:-1])+'/temp/'+'foldx_commands_file.txt','r') as f:
		for line in f:
			all_commands.append(line.strip().split(' '))
	#print(all_commands)

	output_dir = args_.output_argument[0]
	pdb_info = args_.input_argument[1] # protein structure in pdb format
	pdb_name = pdb_info.split('/')[-1]
	pdb_dir = ('/').join(pdb_info.split('/')[0:-1])
	pdb_name_new =  pdb_name.split('.')[0]+'_Repair.pdb'
	foldx_path = args_.input_argument[0]

	if os.path.exists(output_dir + '/' + pdb_name_new):
		print('Repair file already present')
	else:
		all_info = foldx_path + '/foldx --command=RepairPDB --pdb-dir=' + pdb_dir + ' --pdb=' + pdb_name + ' --pH=7 --vdwDesign=2 --pdbHydrogens=false --output-dir=' + output_dir
		os.system(all_info) # Repairing the pdb file

	futures = ray.get([run_foldx.remote(i) for i in all_commands])
	del futures
	os.chdir('/'.join(mutant_file.split('/')[:-1])+'/')
	cmd = f'rm -r {args_.input_argument[2].split("/")[-1].split(".")[0]}_*'
	subprocess.run(cmd, shell=True)
	os.chdir(wd)
	
	
def auto_garbage_collect(pct=60.0):
    if psutil.virtual_memory().percent >= pct:
        gc.collect()
    return
auto_garbage_collect()
ray.shutdown()
print('Completed!')
