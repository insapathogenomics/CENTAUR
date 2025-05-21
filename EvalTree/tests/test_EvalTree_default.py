#!/usr/bin/env	python3

"""
EvalTree tests

By Joana Gomes Pereira
@INSA
"""

import os
import shutil
import subprocess
import pandas as pd
import glob

def test_arguments_default():
  
    main_path=os.path.dirname(__file__)
    path=os.path.dirname(main_path) 
    script_path = os.path.join(path,"EvalTree.py")
    input1=os.path.join(main_path,'input1') 
    input2=os.path.join(main_path,'input2')
    output=os.path.join(main_path,'TEST1') 
    if not os.path.exists(output):
        os.mkdir(output)

    cmd= f"python {script_path} -i1 {input1} -i2 {input2} -o {output} -cp country -pt MST-7x1.0"
    subprocess.check_output(cmd, shell=True)

    GT_clusters_partition=pd.read_csv(os.path.join(output,"GT_clusters_partitions.tsv"), sep='\t', header=None)
    assert GT_clusters_partition.iloc[4,2]=='10'
    num_lines, num_columns = GT_clusters_partition.shape
    assert num_lines == 102
    assert num_columns == 3

    GT_metrics=pd.read_csv(os.path.join(output,"GT_metrics.tsv"), sep='\t', header=None)
    assert GT_metrics.iloc[3,3]=='0.9230769230769231'
    num_lines, num_columns = GT_metrics.shape
    assert num_lines == 101
    assert num_columns == 11

    GT_stable_regions=pd.read_csv(os.path.join(output,"GT_stableRegions.tsv"), sep='\t', comment="#", header=None)
    num_lines, num_columns = GT_stable_regions.shape
    print(f"Rows: {num_lines}, Columns: {num_columns}")
    print(GT_stable_regions)
    assert num_lines == 4
    assert num_columns == 4
    assert GT_stable_regions.iloc[1,1]=='MST-16x1.0->MST-15x1.0'

#Pipeline HC
    HC_clusters_partition=pd.read_csv(os.path.join(output,"HC_clusters_partitions.tsv"), sep='\t', header=None)
    assert HC_clusters_partition.iloc[25,2]=='8'
    num_lines, num_columns = GT_clusters_partition.shape
    assert num_lines == 102
    assert num_columns == 3

    HC_metrics=pd.read_csv(os.path.join(output,"HC_metrics.tsv"), sep='\t', header=None)
    assert HC_metrics.iloc[3,2]=='0.48'
    num_lines, num_columns = HC_metrics.shape
    assert num_lines == 101
    assert num_columns == 11

    HC_stable_regions=pd.read_csv(os.path.join(output,"HC_stableRegions.tsv"), sep='\t',  comment="#", header=None)
    assert HC_stable_regions.iloc[1,1]=='single-16x1.0->single-15x1.0'
    num_lines, num_columns = HC_stable_regions.shape
    assert num_lines == 5
    assert num_columns == 4
    

# Pipeline GT_vs_HC   
    GT_vs_HC_ALL_CORRESPONDENCE=pd.read_csv(os.path.join(output,"GT_vs_HC_All_correspondence.tsv"), sep='\t',  header=None) 
    assert GT_vs_HC_ALL_CORRESPONDENCE.iloc[2,2] =='1'
    num_lines, num_columns = GT_vs_HC_ALL_CORRESPONDENCE.shape
    assert num_lines == 203
    assert num_columns == 3

    GT_vs_HC_final_score=pd.read_csv(os.path.join(output,"GT_vs_HC_final_score.tsv"), sep='\t', header=None)
    assert GT_vs_HC_final_score.iloc[4,4]=='3.0'
    num_lines, num_columns = GT_vs_HC_final_score.shape
    assert num_lines ==102
    assert num_columns ==102

    GT_vs_HC_slope=pd.read_csv(os.path.join(output,"GT_vs_HC_slope.tsv"), sep='\t', header=None)
    assert GT_vs_HC_slope.iloc[2,6]=='0.0032116458079940423'
    num_lines, num_columns = GT_vs_HC_slope.shape
    assert num_lines == 3
    assert num_columns == 7

    GT_vs_HC_Simposons2=pd.read_csv(os.path.join(output,"GT_vs_HC_Simpsons2.tsv"), sep='\t', header=None)
    assert GT_vs_HC_Simposons2.iloc[4,4]==9.0
    num_lines, num_columns = GT_vs_HC_Simposons2.shape
    assert num_lines == 101
    assert num_columns == 5

    all_files=len(os.listdir(output))
    images = glob.glob(f"{output}/*.png")
    tsv=glob.glob(f"{output}/*.tsv")
    html=glob.glob(f"{output}/*.html")

    total_images=len(images)
    total_tsv=len(tsv)
    nr_html=len(html)


    assert all_files == 34
    assert total_images == 11
    assert total_tsv == 17
    assert nr_html == 1


    file_path = os.path.join(main_path, "TEST1")
    shutil.rmtree(file_path)

