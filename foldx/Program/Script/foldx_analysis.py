# importing libraries
import os
import pandas as pd
import argparse

# cases for a mutation to be stabilizing, destabilizing or neutral
feature_cases = {
    'stabilizing_case':[sorted(['stabilizing','stabilizing','stabilizing']),sorted(['stabilizing','stabilizing','destabilizing']),sorted(['stabilizing','stabilizing','neutral'])],
    'destabilizing_case':[sorted(['destabilizing','destabilizing','destabilizing']),sorted(['destabilizing','destabilizing','stabilizing']),sorted(['destabilizing','destabilizing','neutral'])],
    'neutral_case':[sorted(['neutral','neutral','neutral']),sorted(['neutral','neutral','stabilizing']),sorted(['neutral','neutral','destabilizing'])],
    'exception':[sorted(['neutral','stabilizing','destabilizing'])]
}

exception_cases = []

class foldx_data_handler:
    def __init__(self,wd,protein_name,mutation_files,foldx_files):
        self.wd = wd
        self.protein_name = protein_name # full file of protein in .pdb format
        self.mutation_files = mutation_files # pass in directory without '/' for the mutation files present in the folder
        self.foldx_files = foldx_files # pass in directory of foldx files without '/' for the foldx 'Dif_' files present in the folder 

    # conditions for checking the stability of the mutations.
    def stability_aa(self,data):
        stab_list = ['neutral','stabilizing','destabilizing']
        if data < 0.46 and data >= -0.46:
            return(stab_list[0])
        if data < -0.46:
            return(stab_list[1]) 
        if data >= 0.46:
            return(stab_list[2])
        
    # feature numbers according to the stability of the mutations.
    def feature_number(self,data):
        if data == 'stabilizing':
            return 1
        if data == 'destabilizing':
            return 0
        if data == 'neutral':
            return 1
        
    def rename_Dif_foldx_files(self):
        # rename the fxout files according to the folders
        os.chdir(self.foldx_files)
        file_path = os.getcwd()
        for file_name in os.listdir():
            os.chdir('{}/{}'.format(file_path,file_name))
            for fxout in os.listdir():
                if fxout.startswith('Dif_'):
                    cmd = f'mv Dif_{self.protein_name.split("/")[-1].split(".")[0]}_Repair.fxout Dif_{self.protein_name.split("/")[-1].split(".")[0]}_Repair_{file_name}.fxout'
                    os.system(cmd)
        os.chdir(self.wd)
        
    def create_mutation_dict(self):
        # directory of foldx mutation list
        os.chdir(self.wd)
        os.chdir(self.mutation_files) # here must pass the folder to get mutation list. 
        mut_dict = {}
        for file_name in os.listdir():
            mutation_list = []
            #print(file_name)
            with open(file_name,'r') as f:
                for line in f:
                    for i in range(3):
                        mutation_list.append(line.strip())
                mut_dict[file_name] = mutation_list
        os.chdir(self.wd)
        return mut_dict
    
    # funtion to match the mutation file number and 'Dif_' foldx file out number
    def match_mut_list_and_df(self,mut_list_name,fxout_file_name):
        for i,j in mut_list_name.items():
            num = fxout_file_name.split('_')[-1].split('.')[0]
            for name in mut_list_name.keys():
                mut_num = name.split('_')[-1].split('.')[0]
                if num == mut_num:
                    return mut_list_name[name]
    
    def check_if_any_item_matches_condition(self,item, condition):
      for condition_item in condition:
        if item == condition_item:
          return True
      return False

  
    # main code begins here
    # concatenate all the foldx files into single file
    def foldx_dataframe(self):
        all_df = pd.DataFrame()
        os.chdir(self.foldx_files) # foldx working directory
        file_path = os.getcwd()
        for file_name in os.listdir():
            os.chdir('{}/{}'.format(file_path,file_name))
            for fxout in os.listdir():
                if fxout.startswith('Dif_'):
                    print(fxout)
                    csv_lists = []
                    with open(fxout,'r') as f:
                        for line in f.readlines():
                            if line.startswith('Pdb'):
                                #print(line.split('\t'))
                                csv_lists.append(line.split('\t'))
                            if line.startswith(self.protein_name.split('/')[-1].split('.')[0]+'_'):
                                #print(line.split('\t'))
                                csv_lists.append(line.split('\t'))
                    df = pd.DataFrame.from_records(csv_lists)
                    df.columns = df.iloc[0]
                    df = df.iloc[1:,0:2]
                    df = df.rename(columns={'total energy':'delta_delta_G'})
                    df['mut_interface_aa'] = self.match_mut_list_and_df(mut_list_name=self.create_mutation_dict(),fxout_file_name=fxout)
                    df['delta_delta_G'] = df['delta_delta_G'].astype(float)
                    df['feature'] = df['delta_delta_G'].apply(lambda x:self.stability_aa(x))
                    df['feature_num'] = df['feature'].apply(lambda x:self.feature_number(x))
                    df['folder_name'] = str(file_name)
                    all_df = pd.concat([all_df,df])
            os.chdir(self.wd)

        #print(all_df)
        analysis_df = pd.DataFrame()
        
        for index,row in all_df.groupby(by='mut_interface_aa'):

            stabalizing_condition = self.check_if_any_item_matches_condition(sorted(list(row['feature'])),feature_cases['stabilizing_case'])
            destabalizing_condition = self.check_if_any_item_matches_condition(sorted(list(row['feature'])),feature_cases['destabilizing_case'])
            neutral_condition = self.check_if_any_item_matches_condition(sorted(list(row['feature'])),feature_cases['neutral_case'])
            exception_condition = self.check_if_any_item_matches_condition(sorted(list(row['feature'])),feature_cases['exception'])

            if stabalizing_condition:
                min_row = row[row['feature'] == 'stabilizing']['delta_delta_G'].idxmin()
                analysis_df=pd.concat([analysis_df,row.loc[min_row]],axis=1)
            if destabalizing_condition:
                max_row = row[row['feature']== 'destabilizing']['delta_delta_G'].idxmax()
                analysis_df=pd.concat([analysis_df,row.loc[max_row]],axis=1)
            if neutral_condition:
                any_row = row[row['feature']== 'neutral']['delta_delta_G'].idxmax()
                analysis_df=pd.concat([analysis_df,row.loc[any_row]],axis=1)
            if exception_condition:
                #print(row)
                exception_cases.append(row.to_dict())   

        analysis_df = analysis_df.T
        analysis_df['num'] = analysis_df['Pdb'].apply(lambda x: str(x).split('_')[3])
        analysis_df['num'] = analysis_df['num'].astype(int)
        analysis_df = analysis_df.sort_values(by='num',ascending=True).drop('num',axis=1).reset_index()
        analysis_df = analysis_df.drop('index',axis=1)
        analysis_df = analysis_df.rename(columns={'mut_interface_aa':'mutation'})
        analysis_df['mutation'] = analysis_df['mutation'].apply(lambda x:x[:-1])
        analysis_df['mutation'] = analysis_df['mutation'].apply(lambda x :x.split(','))


        return analysis_df






if __name__=='__main__':
    parser = argparse.ArgumentParser('A script to conduct analysis on the foldx file')
    parser.add_argument('-i',nargs='+',dest='input_arguement',default='somthing.txt',help='something helpfile')
    parser.add_argument('-o',nargs='+',dest ='output_arguement',default='something.txt',help='something helpfile')
    args_ = parser.parse_args()
    input_files = args_.input_arguement
    output_files = args_.output_arguement
    foldx_builder = foldx_data_handler(wd= os.getcwd(),
                                   protein_name=input_files[0],
                                   mutation_files=input_files[1],
                                   foldx_files=input_files[2])

    print('Step 1 - Creating mutation dictionary to match the foldx files')
    foldx_builder.create_mutation_dict()
    print('Step 2 - Renaming Dif files according to their folder number')
    foldx_builder.rename_Dif_foldx_files()
    print('Step 3 - Building Foldx dataframe')
    df = foldx_builder.foldx_dataframe()
    print('Step 4 - Saving as .csv file')
    df.to_csv(f'{output_files[0]}/df_{input_files[0].split("/")[-1].split(".")[0]}_foldx.csv')
    print('Done!')

    with open(f'{output_files[0]}/exception_cases.txt','w') as f:
        for cases in exception_cases:
            f.write(str(cases))
            f.write('\n')

    print('exception cases are also stored')