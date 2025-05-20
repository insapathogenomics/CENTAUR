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
import fnmatch
import glob


def test_arguments():
    
    main_path=os.path.dirname(__file__)  
    path=os.path.dirname(main_path)
    script_path = os.path.join(path,"evaltree.py")
    input1=os.path.join(main_path,'input1')  
    input2=os.path.join(main_path,'input2')
    output=os.path.join(main_path,'TEST2')
    if not os.path.exists(output):
        os.mkdir(output)
 	
    result= subprocess.check_output(f"python {script_path} -i1 {input1} -i2 {input2} -o {output} -s 2.95 -t 2-28 -ps partitions_summary -pt MST-7x1.0,MST-8x1.0 -cp country,source -n 4 -to ""MST-7x1.0,MST-7x1.0\;MST-8x1.0,\<=MST-20x1.0""  ", shell=True)

# check input 1 files (GT)

    GT_partitions=pd.read_csv(os.path.join(input1,"GT_partitions.tsv"), sep='\t', header=None)
    assert GT_partitions.iloc[12,0] =='L'
    assert GT_partitions.iloc[2,12] =='cluster_1'
    num_lines, num_columns=GT_partitions.shape
    assert num_lines ==14
    assert num_columns==102
                                                           
    GT_cluster_Composition=pd.read_csv(os.path.join(input1,"GT_clusterComposition.tsv"), sep='\t', header=None)
    assert GT_cluster_Composition.iloc[2,0] =='MST-0x1.0'
    assert GT_cluster_Composition.iloc[69490,0] =='MST-1728x1.0'
    num_lines, num_columns=GT_cluster_Composition.shape
    assert num_lines ==69491
    assert num_columns==4

    GT_partitions_summary=pd.read_csv(os.path.join(input1,"GT_partitions_summary.tsv"), sep='\t', header=None)
    assert GT_partitions_summary.iloc[2,0] =='MST-4x1.0'
    assert GT_partitions_summary.iloc[363,0] =='MST-15x1.0'
    num_lines, num_columns=GT_partitions_summary.shape
    assert num_lines ==364
    assert num_columns==10

    GT_sample_of_interest=pd.read_csv(os.path.join(input1,'GT_SAMPLES_OF_INTEREST_partitions_summary.tsv'), sep='\t', header=None)
    assert GT_sample_of_interest.iloc[2,0] =='sample_0675'
    assert GT_sample_of_interest.iloc[8,9] == 'clinical (100.0%) (n = 24)'
    num_lines, num_columns=GT_sample_of_interest.shape
    assert num_lines ==9
    assert num_columns==11

# check input 2 files (HC)

    HC_partitions=pd.read_csv(os.path.join(input2,"HC_partitions.tsv"), sep='\t', header=None)
    assert HC_partitions.iloc[3,0] =='C'
    assert HC_partitions.iloc[2,1] =='singleton_2'
    num_lines, num_columns=HC_partitions.shape
    assert num_lines ==14
    assert num_columns==102
                                                           
    HC_cluster_Composition=pd.read_csv(os.path.join(input2,"HC_clusterComposition.tsv"), sep='\t', header=None)
    assert HC_cluster_Composition.iloc[2,0] =='MST-0x1.0'
    assert HC_cluster_Composition.iloc[69490,0] =='MST-1728x1.0'
    num_lines, num_columns=HC_cluster_Composition.shape
    assert num_lines ==69491
    assert num_columns==4

    HC_partitions_summary=pd.read_csv(os.path.join(input2,"HC_partitions_summary.tsv"), sep='\t', header=None)
    assert HC_partitions_summary.iloc[2,0] =='MST-4x1.0'
    assert HC_partitions_summary.iloc[363,0] =='MST-15x1.0'
    num_lines, num_columns=HC_partitions_summary.shape
    assert num_lines ==364
    assert num_columns==10

    HC_sample_of_interest=pd.read_csv(os.path.join(input2,'HC_SAMPLES_OF_INTEREST_partitions_summary.tsv'), sep='\t', header=None)
    assert HC_sample_of_interest.iloc[2,0] =='sample_0675'
    assert HC_sample_of_interest.iloc[8,9] == 'clinical (100.0%) (n = 24)'
    num_lines, num_columns=HC_sample_of_interest.shape
    assert num_lines ==9
    assert num_columns==11

#######

    GT_clusters_partition_filtered=pd.read_csv(os.path.join(output,"GT_clusters_partitions-filtered.tsv"), sep='\t', header=None)
    assert GT_clusters_partition_filtered.iloc[3,2]=='10'
    num_lines, num_columns = GT_clusters_partition_filtered.shape
    assert num_lines == 28
    assert num_columns == 3

    GT_partitions_filtered=pd.read_csv(os.path.join(output,"GT_partitions-filtered.tsv"), sep='\t', header=None)
    assert GT_partitions_filtered.iloc[2,2]=='cluster_1'  
    assert GT_partitions_filtered.iloc[5,2]=='singleton_1'
    num_lines, num_columns = GT_partitions_filtered.shape
    assert num_lines == 14
    assert num_columns == 28

    GT_stable_regions=pd.read_csv(os.path.join(output,"GT_stableRegions.tsv"), sep='\t', comment="#", header=None)
    num_lines, num_columns = GT_stable_regions.shape
    print(f"Rows: {num_lines}, Columns: {num_columns}")
    assert num_lines == 2
    assert num_columns == 4
    assert GT_stable_regions.iloc[1,1]=='MST-16x1.0->MST-15x1.0'

    GT_metrics=pd.read_csv(os.path.join(output,"GT_metrics.tsv"), sep='\t', header=None)
    assert GT_metrics.iloc[3,3]=='0.8717948717948718'
    num_lines, num_columns = GT_metrics.shape
    assert num_lines == 27
    assert num_columns == 11

# check HC
    HC_partitions_filtered=pd.read_csv(os.path.join(output,"HC_partitions-filtered.tsv"), sep='\t', header=None)
    assert HC_partitions_filtered.iloc[2,2]=='singleton_2'
    assert HC_partitions_filtered.iloc[5,2]=='singleton_3'

    HC_clusters_partition_filtered=pd.read_csv(os.path.join(output,"HC_clusters_partitions-filtered.tsv"), sep='\t', header=None)
    assert HC_clusters_partition_filtered.iloc[7,2]=='9'
    num_lines, num_columns = HC_clusters_partition_filtered.shape
    assert num_lines == 28
    assert num_columns == 3

    HC_metrics=pd.read_csv(os.path.join(output,"HC_metrics.tsv"), sep='\t', header=None)
    assert HC_metrics.iloc[3,2]=='1.0'
    num_lines, num_columns = HC_metrics.shape
    assert num_lines == 27
    assert num_columns == 11

    HC_stable_regions=pd.read_csv(os.path.join(output,"HC_stableRegions.tsv"), sep='\t',  comment="#", header=None)
    assert HC_stable_regions.iloc[1,1]=='single-16x1.0->single-15x1.0'
    num_lines, num_columns = HC_stable_regions.shape
    assert num_lines == 2
    assert num_columns == 4

    HC_clusters_partitions_filtered=pd.read_csv(os.path.join(output,"HC_clusters_partitions-filtered.tsv"), sep='\t', header=None)
    assert HC_clusters_partitions_filtered.iloc[2,2]=='10'
    assert HC_clusters_partitions_filtered.iloc[27,2]=='8'
    num_lines, num_columns=HC_clusters_partitions_filtered.shape
    assert num_lines ==28
    assert num_columns==3

                                                #Pipeline GT_vs_HC
    print('------------------------------------------------------------------------')
    print(output)
    GT_vs_HC_final_score=pd.read_csv(os.path.join(output,"GT_vs_HC_final_score.tsv"), sep='\t', header=None)
    assert GT_vs_HC_final_score.iloc[4,5]=='3.0'
    num_lines, num_columns=GT_vs_HC_final_score.shape
    assert num_lines == 28
    assert num_columns == 28

    GT_vs_HC_ALL_CORRESPONDENCE=pd.read_csv(os.path.join(output,"GT_vs_HC_All_correspondence.tsv"), sep='\t', header=None) 
    assert GT_vs_HC_ALL_CORRESPONDENCE.iloc[2,2] =='1'
    num_lines,num_columns=GT_vs_HC_ALL_CORRESPONDENCE.shape
    assert num_lines == 55
    assert num_columns == 3

    GT_vs_HC_AdjustedRand=pd.read_csv(os.path.join(output,"GT_vs_HC_AdjustedRand.tsv"), sep='\t', header=None)
    num_lines, num_columns = GT_vs_HC_AdjustedRand.shape
    assert num_lines == 28
    assert num_columns == 28 
    assert GT_vs_HC_AdjustedRand.iloc[1,1]=='1.0'

    GT_vs_HC_AdjWallace1=pd.read_csv(os.path.join(output,"GT_vs_HC_AdjWallace1.tsv"), sep='\t', header=None)
    num_lines, num_columns = GT_vs_HC_AdjWallace1.shape
    assert num_lines ==28
    assert num_columns ==28
    assert GT_vs_HC_AdjWallace1.iloc[1,1]=='1.0'

    GT_vs_HC_AdjWallace2=pd.read_csv(os.path.join(output,"GT_vs_HC_AdjWallace2.tsv"), sep='\t', header=None)
    num_lines, num_columns = GT_vs_HC_AdjWallace2.shape
    assert num_lines ==28
    assert num_columns ==28
    assert GT_vs_HC_AdjWallace2.iloc[1,2]=='0.3157894736842105'

    GT_vs_HC_path_stats_outbreak=pd.read_csv(os.path.join(output,"GT_vs_HC_path_stats_outbreak.tsv"), sep='\t', header=None)
    num_lines, num_columns = GT_vs_HC_path_stats_outbreak.shape
    assert num_lines==2
    assert num_columns==1
    assert GT_vs_HC_path_stats_outbreak.iloc[1,0]== f'{input2}/HC_clusterComposition.tsv'

    GT_vs_HC_Simpsons1=pd.read_csv(os.path.join(output,"GT_vs_HC_Simpsons1.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_Simpsons1.shape
    assert num_lines==27
    assert num_columns==5
    assert GT_vs_HC_Simpsons1.iloc[0,0]=='MST-2x1.0'
    
    GT_vs_HC_Simpsons2=pd.read_csv(os.path.join(output,"GT_vs_HC_Simpsons2.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_Simpsons2.shape
    assert num_lines==27
    assert num_columns==5
    assert GT_vs_HC_Simpsons2.iloc[26,0]=='single-28x1.0'

    GT_vs_HC_slope=pd.read_csv(os.path.join(output,"GT_vs_HC_slope.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_slope.shape
    assert num_lines==3
    assert num_columns==7
    assert GT_vs_HC_slope.iloc[2,6]=='0.00465957700368776'

    GT_vs_HC_stats_outbreak_pairwise_comparison_7_equal_7=pd.read_csv(os.path.join(output,"GT_vs_HC_stats_outbreak_pairwise_comparison_7_equal_7.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_stats_outbreak_pairwise_comparison_7_equal_7.shape
    assert num_lines==3
    assert num_columns==3
    assert GT_vs_HC_stats_outbreak_pairwise_comparison_7_equal_7.iloc[2,2]=='126' 

    GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_20=pd.read_csv(os.path.join(output,"GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_20.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_20.shape
    assert num_lines==3
    assert num_columns==3
    assert GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_20.iloc[2,2]=='120'

    GT_vs_HC_stats_outbreak_summary_7_equal_7=pd.read_csv(os.path.join(output,"GT_vs_HC_stats_outbreak_summary_7_equal_7.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_stats_outbreak_summary_7_equal_7.shape
    assert num_lines==3
    assert num_columns==6
    assert GT_vs_HC_stats_outbreak_summary_7_equal_7.iloc[2,2]=='126'

    GT_vs_HC_stats_outbreak_summary_8_lower_equal_20=pd.read_csv(os.path.join(output,"GT_vs_HC_stats_outbreak_summary_8_lower_equal_20.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_stats_outbreak_summary_8_lower_equal_20.shape
    assert num_lines==3
    assert num_columns==6
    assert GT_vs_HC_stats_outbreak_summary_8_lower_equal_20.iloc[2,2]=='120'

    GT_vs_HC_Wallace1=pd.read_csv(os.path.join(output,"GT_vs_HC_Wallace1.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_Wallace1.shape
    assert num_lines==28
    assert num_columns==28
    assert GT_vs_HC_Wallace1.iloc[2,2]=='1.0'

    GT_vs_HC_Wallace2=pd.read_csv(os.path.join(output,"GT_vs_HC_Wallace2.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_Wallace2.shape
    assert num_lines==28
    assert num_columns==28
    assert GT_vs_HC_Wallace2.iloc[27,3]=='0.6666666666666666'

    #----------------------------------------
    # Number of files

    all_files=len(os.listdir(output))
    images = glob.glob(f"{output}/*.png")
    files=glob.glob(f"{output}/*.tsv")
    html=glob.glob(f"{output}/*.html")

    total_images=len(images)
    total_files=len(files)
    nr_html=len(html)
    
    assert all_files ==57
    assert total_images == 23
    assert total_files==28
    assert nr_html ==1

    
def test_RTO_argument():

    main_path=os.path.dirname(__file__) 
    path=os.path.dirname(main_path)
    script_path = os.path.join(path,"evaltree.py")
    input1=os.path.join(main_path,'input1')  
    input2=os.path.join(main_path,'input2')
    output=os.path.join(main_path,'TEST2')

    if not os.path.exists(output):
        os.mkdir(output)
 
    # Teste for the -RTO argument
    result= subprocess.check_output(f"python {script_path} -i1 {input1} -i2 {input2} -o {output} -s 2.95 -t 2-28"
                                    f" -ps partitions_summary -pt MST-7x1.0,MST-8x1.0 -cp country,source -n 4 -to \"MST-5x1.0,MST-5x1.0;MST-8x1.0,<=MST-10x1.0\" -rto ",
                                    shell=True)

    GT_vs_HC_stats_outbreak_pairwise_comparison_5_equal_5=pd.read_csv(os.path.join(output,"GT_vs_HC_stats_outbreak_pairwise_comparison_5_equal_5_pct.tsv"), sep='\t', header=None)
    num_lines, num_columns =GT_vs_HC_stats_outbreak_pairwise_comparison_5_equal_5.shape
    assert num_lines==3
    assert num_columns==3
    assert GT_vs_HC_stats_outbreak_pairwise_comparison_5_equal_5.iloc[2,2]=='1.0'

    GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_10=pd.read_csv(os.path.join(output,"GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_10_pct.tsv"), sep='\t', header=None)
    num_lines, num_columns=GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_10.shape
    assert num_lines==3
    assert num_columns==3
    assert GT_vs_HC_stats_outbreak_pairwise_comparison_8_lower_equal_10.iloc[2,2]=='1.0'

    #Test HTML
    all_files=os.listdir(output)
    list_files=[]
    original_html=f'GT_vs_HC_report.html'
    new_html=f'GT_vs_HC_NEW_report.html'
    for file in all_files:
        if file ==original_html:
            list_files.append(file)
        if file ==new_html:
            list_files.append(file)
    assert len (list_files) == 2 


def test_plots_category_percentage():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1')
    input2 = os.path.join(main_path, 'input2')
    output = os.path.join(main_path, 'TEST2', 'plots_category_percentage')
    
    if not os.path.exists(output):
        os.mkdir(output)


    result = subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output} -s 2.95 -t 2-28 "
        f"-ps partitions_summary -pt MST-7x1.0,MST-8x1.0 -cp country,source -n 4 -to \"MST-5x1.0,MST-5x1.0;MST-8x1.0,<=MST-10x1.0\" -pcp 60",
        shell=True
    )

    all_files = os.listdir(output)
    country_images = []
    source_images = []

    string_1 = '*_country_Cluster_*.png'
    string_2 = '*_source_Cluster_*.png'  

    for file in all_files:
        if fnmatch.fnmatch(file, string_1):
            country_images.append(file)
        if fnmatch.fnmatch(file, string_2):
            source_images.append(file)  

    assert len(country_images) == 8
    assert len(source_images) == 8

    shutil.rmtree(output)

#Folder and sequence type matrix

def test_sequence_type_and_folder():
    
    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1','GT_1.tsv')
    input2 = os.path.join(main_path, 'input2')
    output = os.path.join(main_path, 'TEST2', 'st_vs_folder')

    if not os.path.exists(output):
        os.mkdir(output)

    result=subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output} -s 2.95 -t 2-28 "
        f"-ps partitions_summary -pt MST-7x1.0,MST-8x1.0 -cp country,source -n 4", shell=True)

    all_files=len(os.listdir(output))
    images = glob.glob(f"{output}/*.png")
    files=glob.glob(f"{output}/*.tsv")
    html=glob.glob(f"{output}/*.html")

    total_images=len(images)
    total_files=len(files)
    nr_html=len(html)
    
    assert all_files == 30
    assert total_images == 12
    assert total_files == 13
    assert nr_html == 1
    shutil.rmtree(output)

def mx_partition_vs_mx_partitions():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path,'input1' ,'GT_partitions.tsv')
    input2 = os.path.join(main_path, 'input2' ,'HC_partitions.tsv')
    output = os.path.join(main_path, 'TEST2', 'mx_vs_mx')

    if not os.path.exists(output):
        os.mkdir(output)

    result=subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output}",shell=True)
    
    all_files=len(os.listdir(output))
    images = glob.glob(f"{output}/*.png")
    files=glob.glob(f"{output}/*.tsv")
    html=glob.glob(f"{output}/*.html")

    total_images=len(images)
    total_files=len(files)
    nr_html=len(html)

    assert all_files ==29
    assert total_images == 4
    assert total_files==19
    assert nr_html ==1

    shutil.rmtree(output)

def test_plots_category_number():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1')
    input2 = os.path.join(main_path, 'input2')
    output = os.path.join(main_path, 'TEST2', 'plots_category_number')

    if not os.path.exists(output):
        os.mkdir(output)

    result = subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output} -s 2.95 -t 2-28 "
        f"-ps partitions_summary -pt MST-7x1.0 -cp country -n 4 -pcn 1",
        shell=True
    )

    all_files=len(os.listdir(output))
    images = glob.glob(f"{output}/*.png")
    files=glob.glob(f"{output}/*.tsv")
    html=glob.glob(f"{output}/*.html")

    total_images=len(images)
    total_files=len(files)
    nr_html=len(html)

    assert all_files ==38
    assert total_images == 13
    assert total_files==19
    assert nr_html ==1

    shutil.rmtree(output)


def test_folder_vs_folder():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1')
    input2 = os.path.join(main_path, 'input2')
    output = os.path.join(main_path, 'TEST2', 'folder_vs_folder')

    if not os.path.exists(output):
        os.mkdir(output)

    result = subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output}",
        shell=True
    )

    all_files=len(os.listdir(output))
    assert all_files == 28
    shutil.rmtree(output)

def test_file_vs_folder():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1')
    input2 = os.path.join(main_path, 'input2','HC_2.tsv')
    output = os.path.join(main_path, 'TEST2', 'file_vs_folder')

    if not os.path.exists(output):
        os.mkdir(output)

    result = subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output}",shell=True)

    all_files=len(os.listdir(output))
    assert all_files == 21

def test_file_vs_file():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1','GT_1.tsv')
    input2 = os.path.join(main_path, 'input2','HC_2.tsv')
    output = os.path.join(main_path, 'TEST2', 'file1_vs_file2')

    if not os.path.exists(output):
        os.mkdir(output)

    result = subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output}", shell=True)

    all_files=len(os.listdir(output))
    assert all_files == 16

def test_file_mx_partition():

    main_path = os.path.dirname(__file__)  
    path = os.path.dirname(main_path)  
    script_path = os.path.join(path, "evaltree.py")
    input1 = os.path.join(main_path, 'input1','GT_1.tsv')
    input2 = os.path.join(main_path, 'input2','HC_partitions.tsv')
    output = os.path.join(main_path, 'TEST2', 'file_vs_file')

    if not os.path.exists(output):
        os.mkdir(output)

    result = subprocess.check_output(
        f"python {script_path} -i1 {input1} -i2 {input2} -o {output}",
        shell=True
    )

    all_files=len(os.listdir(output))
    assert all_files == 21


    file_path = os.path.join(main_path, "TEST2")
    shutil.rmtree(file_path)
    