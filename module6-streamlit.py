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
temp = psa_file.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
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
            temp = af3_ps.getvalue().decode("utf-8") ##decodes characters correctly but still has too long file name issue
            temp_split = temp.splitlines()
            for line in temp_split:
                if line.startswith('ATOM'):
                    pLDDT_line = re.sub(r'\s+', ' ', line)
                    pLDDT_ps += [pLDDT_line]
            #st.write(pLDDT_ps)
            #for idx, item in enumerate(pLDDT_ps):
                #split = item.split(' ')
                #for item in split:
                    #st.write(idx)
                    #st.write(item)
                #st.write(len(split))
            #turn list from parsed cif file into dataframe
            af3ps_df = pd.DataFrame(pLDDT_ps, columns=['atom'])
            st.write(af3ps_df)
            #explode dataframe into multiple columns to isolate pLDDT scores
            af3ps_df[['0','1','2','3','4','resn','6','7','pos','9','10','11','12','13','pLDDT','15','16','17','18']] = af3ps_df['atom'].str.split(' ',expand=True) #\t for tab delimited
            st.write(af3ps_df[['resn','pos','pLDDT']])

            residue_df = ''
            temp_list = []
            pLDDT_averages = []
            resn = []
            #make a list w numbers for 1 through number of residues
            #st.write('last residue position')
            #st.write(af3ps_df.iloc[-1,9])
            num_resi = list(range(1,int(af3ps_df.iloc[-1,9])))
            #st.write(len(num_resi))
    
            for i in num_resi:
                temp_list = []
                st.write(i)
                #if int(af3ps_df.loc['pos']) == i:
                    #st.write(af3ps_df.iloc[idx,6])
                    #resn.append(af3ps_df.iloc[idx,6])
                for idx, row in enumerate(af3ps_df):
                    if int(af3ps_df.iloc[idx,9]) == i:
                        st.write(idx)
                        st.write(af3ps_df.iloc[idx,15])
                        temp_list.append(float(af3ps_df.iloc[idx,15]))
                        st.write(temp_list)
                        #st.write(type(temp_list))
                st.write(idx)
                st.write(af3ps_df.iloc[idx,6])
                resn.append(af3ps_df.iloc[idx,6])
                avg = np.mean(temp_list)
                #st.write(avg)
                pLDDT_averages.append(avg)
            st.write(pLDDT_averages)
            st.write(resn)
            residue_df['resn'] = [resn]
            residue_df['position'] = [num_resi]
            residue_df['pLDDT'] = [pLDDT_averages]
            st.write(residue_df)
                    
                    

            #consurf_df = consurf_df[['SEQ','COLOR']]
            #consurf_df = consurf_df.iloc[1:].reset_index(drop=True)
            #st.write(consurf_df)

            #af3t_df = pd.read_csv(af3_target)
            #st.write(af3t_df)
                        
            #af3_df = pd.concat([af3ps_df, af3t_df], axis=1)
            #st.write(af3_df)
            #df_combined = pd.concat([df_exploded, af3_df], axis=1)
            #st.write(df_combined)
                        
#not sure how much of the below code i will use; will need to create a new column for delta pLDDT & do some math to fill it for specific positions; will want to print target positions that match ps positions and will use ps positions for the if statement indexing for where to create delta values
            for idx, aa in enumerate(df_exploded['Project Standard Seq']):
                #st.write(aa)
                if aa == '-':
                    #st.write(idx)
                    gap = idx
                    #st.write(gap)
                    #st.write(gap - 0.5)
                    #consurf_df.loc[gap] = ''
                    line = DataFrame({"SEQ": '', "COLOR": 0}, index=[gap -0.5])
                    consurf_df = pd.concat([consurf_df, line])
                    consurf_df = consurf_df.sort_index().reset_index(drop=True)
            #st.write(consurf_df)
            df_combined = pd.concat([df_exploded, consurf_df], axis=1)
            df_combined = df_combined.iloc[:-1]
            #st.write(df_combined)

            #create new column (evoscore) and fill cells
            df_combined['EvoScore'] = ''
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
