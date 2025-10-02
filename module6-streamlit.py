import streamlit as st
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from io import StringIO
from Bio import AlignIO
import pandas as pd
import tempfile
import os
import numpy as np
from pandas import DataFrame
import re

st.title('FoldScore Calculation')

st.header('Submit Clustal Alignment')
    
#code to upload clustal file & turn into dataframe

psa_file = st.file_uploader("",type='clustal_num', key=1)
if psa_file is not None:
    st.success("PSA file uploaded")
else:
    st.info("please upload your clustal .clustal file")

st.header('Submit Project Standard AlphaFold3 cif File')

af3_ps = st.file_uploader("",type='cif', key=2)
if af3_ps is not None:
    st.success("project standard AF3 file uploaded")
else:
    st.info("please upload your project standard AF3 .cif file")

st.header('Submit Target AlphaFold3 cif File')

af3_target = st.file_uploader("",type='cif', key=3)
if af3_target is not None:
    st.success("target AF3 file uploaded")
else:
    st.info("please upload your target AF3 .cif file")

#bytes = psa_file.getvalue() ##adds 'b in front of file & other character issues (adds /n etc)
#st.write(bytes)
try:
    temp = psa_file.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
except AttributeError:
    pass
#st.text(temp)

#declaring variables outside of button if statement so i can access them after the button step
df = ''
df1 = ''
df2 = ''
df_exploded = ''

if st.button('read in clustal alignment file'):
    ps_line = ""
    target_line = ""
    temp_split = temp.splitlines()
    # Process the captured Clustal output text line by line
    for line in temp_split[1:]:
        # The conservation line is identifiable by its spacing
        if 'P22259' in line:
            seq = line[26:86]
            ps_line += seq
        if line[:1].isalpha():
            if 'P22259' not in line:
                seq = line[26:86]
                target_line += seq
    ps_line = "".join(char for char in ps_line if not char.isdigit())
    target_line = "".join(char for char in target_line if not char.isdigit())
    st.text(ps_line)
    st.text(target_line)
    
    #convert strings to pandas dataframe
    
    data = {'Target Seq': [target_line],
                'Project Standard Seq': [ps_line]}
    df = pd.DataFrame(data)
    #st.write(df)
    df1 = df['Target Seq'].str.split('').explode().reset_index(drop=True)
    #st.write(df1)
    df2 = df['Project Standard Seq'].str.split('').explode().reset_index(drop=True)
    #st.write(df2)
    df_exploded = pd.concat([df1, df2], axis=1)
    #st.write(df_exploded)
    #df_exploded['color'] = 0
    #df_exploded = df_exploded.iloc[1:].reset_index(drop=True) #moving this to after the conservation symbols are added
    #st.write(df_exploded)

    #process AF3 cif files into a table

    @st.fragment()
    def frag():
        if st.button('create AF3 dataframes & align with clustal'):
            pLDDT_ps = []
            pLDDT_target = []

            #project standard cif parsing
            temp = af3_ps.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
            temp_split = temp.splitlines()
            for line in temp_split:
                if line.startswith('ATOM'):
                    pLDDT_line = re.sub(r'\s+', ' ', line)
                    pLDDT_ps += [pLDDT_line]
            #st.write(pLDDT_ps)
            #turn list from parsed cif file into dataframe
            af3ps_df = pd.DataFrame(pLDDT_ps, columns=['atom'])
            st.write(af3ps_df)
            #explode dataframe into multiple columns to isolate pLDDT scores
            af3ps_df[['0','1','2','3','4','resn','6','7','pos','9','10','11','12','13','pLDDT','15','16','17','18']] = af3ps_df['atom'].str.split(' ',expand=True) #\t for tab delimited
            st.write(af3ps_df[['resn','pos','pLDDT']])
            
            psresidue_df = ''
            temp_list = []
            pLDDT_averages = []
            resn = []
            position = ''
            #make a list w numbers for 1 through number of residues
            num_resi = list(range(1,int(af3ps_df.iloc[-1,9])))


           # Pre-extract the necessary columns from the DataFrame for faster access
            pos_col = af3ps_df.iloc[:, 9].values
            resn_col = af3ps_df.iloc[:, 6].values
            value_col = af3ps_df.iloc[:, 15].values
            
            # Initialize lists for results
            pLDDT_averages = []
            resn = []
            
            # Iterate over each unique `num_resi` in num_resi
            for i in num_resi:
                # Mask to filter rows where the position matches the current `i`
                mask = pos_col == i
                
                # Extract the filtered values for `value_col` and calculate the average
                filtered_values = value_col[mask]
                avg = np.mean(filtered_values)
                pLDDT_averages.append(avg)
                
                # Extract the filtered `resn` values
                filtered_resn = resn_col[mask]
                
                # Track changes in position
                pos_changes = np.diff(pos_col[mask]) != 0  # Detect position changes
                
                # Append the first residue name (because it doesn't have a previous position to compare)
                resn.append(filtered_resn[0])
                
                # Append the resn values at position change points (starting from index 1)
                for idx in range(1, len(pos_changes)):
                    if pos_changes[idx-1]:  # Position change detected
                        resn.append(filtered_resn[idx])
 
            '''
            for i in num_resi:
                temp_list = []
                for idx, row in af3ps_df.iterrows():
                    if int(af3ps_df.iloc[idx,9]) == i:
                        temp_list.append(float(af3ps_df.iloc[idx,15]))
                        pos_temp = af3ps_df.iloc[idx,9]
                        if pos_temp != af3ps_df.iloc[idx-1,9]:
                            resn.append(af3ps_df.iloc[idx,6])
                avg = np.mean(temp_list)
                pLDDT_averages.append(avg)
                '''
            st.write(pLDDT_averages)
            st.write(resn)
            data = {'Project Standard Residue': resn,
                    'Project Standard Position': num_resi,
                    'Project Standard pLDDT': pLDDT_averages}
            psresidue_df = pd.DataFrame(data)
            st.write(psresidue_df)
                    
            #target cif parsing
            temp = af3_target.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
            temp_split = temp.splitlines()
            for line in temp_split:
                if line.startswith('ATOM'):
                    pLDDT_line = re.sub(r'\s+', ' ', line)
                    pLDDT_target += [pLDDT_line]
            #st.write(pLDDT_target)
            
            #turn list from parsed cif file into dataframe
            af3t_df = pd.DataFrame(pLDDT_target, columns=['atom'])
            st.write(af3t_df)
            #explode dataframe into multiple columns to isolate pLDDT scores
            af3t_df[['0','1','2','3','4','resn','6','7','pos','9','10','11','12','13','pLDDT','15','16','17','18']] = af3t_df['atom'].str.split(' ',expand=True) #\t for tab delimited
            st.write(af3t_df[['resn','pos','pLDDT']])
            
            tresidue_df = ''
            temp_list = []
            pLDDT_averages = []
            resn = []
            position = ''
            #make a list w numbers for 1 through number of residues
            num_resi = list(range(1,int(af3t_df.iloc[-1,9])))
    
            for i in num_resi:
                temp_list = []
                for idx, row in af3t_df.iterrows():
                    if int(af3t_df.iloc[idx,9]) == i:
                        temp_list.append(float(af3t_df.iloc[idx,15]))
                        pos_temp = af3t_df.iloc[idx,9]
                        if pos_temp != af3t_df.iloc[idx-1,9]:
                            resn.append(af3t_df.iloc[idx,6])
                avg = np.mean(temp_list)
                pLDDT_averages.append(avg)
            st.write(pLDDT_averages)
            st.write(resn)
            data = {'Target Residue': resn,
                    'Target Position': num_resi,
                    'Target pLDDT': pLDDT_averages}
            tresidue_df = pd.DataFrame(data)
            st.write(tresidue_df)
            
                        
            #af3_df = pd.concat([psresidue_df, tresidue_df], axis=1)
            #st.write(af3_df)
            #df_combined = pd.concat([df_exploded, af3_df], axis=1)
            #st.write(df_combined)
                        
#not sure how much of the below code i will use; will need to create a new column for delta pLDDT & do some math to fill it for specific positions; will want to print target positions that match ps positions and will use ps positions for the if statement indexing for where to create delta values
            df_exploded = df_exploded.drop(index=df_exploded.index[0])
            st.write(df_exploded)
            for idx, aa in enumerate(df_exploded['Project Standard Seq']):
                #st.write(aa)
                if aa == '-':
                    #st.write(idx)
                    gap = idx
                    #st.write(gap)
                    #st.write(gap - 0.5)
                    #consurf_df.loc[gap] = ''
                    line = DataFrame({"Project Standard Residue": '', "Project Standard Position": '', "Project Standard pLDDT": 0}, index=[gap -0.5])
                    psresidue_df = pd.concat([psresidue_df, line])
                    psresidue_df = psresidue_df.sort_index().reset_index(drop=True)
            for idx, aa in enumerate(df_exploded['Target Seq']):
                if aa == '-':
                    #st.write(idx)
                    gap = idx
                    #st.write(gap)
                    #st.write(gap - 0.5)
                    #consurf_df.loc[gap] = ''
                    line = DataFrame({"Target Residue": '', "Target Position": '', "Target pLDDT": 0}, index=[gap -0.5])
                    tresidue_df = pd.concat([tresidue_df, line])
                    tresidue_df = tresidue_df.sort_index().reset_index(drop=True)
            #st.write(consurf_df)
            df_combined = pd.concat([df_exploded, tresidue_df, psresidue_df], axis=1)
            #df_combined = df_combined.iloc[:-1]
            st.write(df_combined)

            #create new column (evoscore) and fill cells
            df_combined['Delta pLDDT'] = ''
            #st.write(df_combined)
            #st.write(df_combined.dtypes)
            for idx,i in enumerate(df_combined['COLOR']):
                if i < 4:
                    df_combined.iloc[idx,5] = 0
                if i >= 4:
                    df_combined.iloc[idx,5] = i
                if df_combined.iloc[idx,0] == df_combined.iloc[idx,1]:
                    df_combined.iloc[idx,5] = 0
            st.write(df_combined)
            st.write(df_combined.dtypes)
            df_combined['EvoScore'] = df_combined['EvoScore'].astype(float)
            #st.write(df_combined.dtypes)
            evoscore = df_combined['EvoScore'].sum()
            st.write('EvoScore = ' + str(evoscore))
            
            df_combined['Weighted EvoScore'] = ''
            for idx, i in enumerate(df_combined['COLOR']):
                if i < 4:
                    df_combined.iloc[idx,6] = 0
                if i >= 4:
                    if df_combined.iloc[idx,2] == '*':
                        df_combined.iloc[idx,6] = 0
                    if df_combined.iloc[idx,2] == ':':
                        df_combined.iloc[idx,6] = i*0.5
                    if df_combined.iloc[idx,2] == '.':
                        df_combined.iloc[idx,6] = i*0.75
                    if df_combined.iloc[idx,2] == ' ':
                        df_combined.iloc[idx,6] = i
            st.write(df_combined)
            st.write(df_combined[['Project Standard Seq', 'Target Seq', 'EvoScore', 'Weighted EvoScore']])
            weighted_evoscore = df_combined['Weighted EvoScore'].sum()
            st.write('Weighted EvoScore = ' + str(weighted_evoscore))

                
            
    frag()


#@st.fragment()
#def PSA_download():
 #   with open('clustalPSA.aln') as f:
  #      st.download_button('download PSA', f)
#PSA_download()
