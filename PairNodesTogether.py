## CODE TO CALCULATE THE NUMBER OF PAIRS THAT ARE KEPT TOGETHER WHEN THE RESOLUTION PARAMETER IS TUNED IN INCREASING STEPS 
## Developed by Núria Olvera Ocaña
## Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS) and Barcelona Supercomputing Center (BSC)
## Inflammation and repair in respiratory diseases group at IDIBAPS and Computational Biology group at BSC
## Contact: olvera@clinic.cat


### LOAD PACKAGES ###
import subprocess
import re
import copy
import os
from argparse import ArgumentParser

""""
args = ArgumentParser()
args.add_argument('--ncol_a', help='Input ncol file', dest='ncol_a',  default=None)
args.add_argument('--ncol_b', help='Input ncol file', dest='ncol_b',  default=None)
args.add_argument('--ncol_c', help='Input ncol file', dest='ncol_c',  default=None)
args.add_argument('--max_reso', help='Maximum resolution to run Molti', dest='reso_max',  default=None)
args.add_argument('--out', help='Prefix of the output files',dest='out', default=None)

op = args.parse_args()

"""
### FUNCTIONS ###
def tuning_resolution(reso_min,reso_max, fileA, output,fileB= None, fileC= None):
    '''
    Parameters: minimal resolution, maximal resolution, ncol format files for the three networks and output file name.
    Returns: cluster files returned by MolTi computed using resolutions parameters in increasing steps of 0.1. 
    
    '''
    for i in range(reso_min, reso_max, 5):
        a=i/10
        if fileB != None and fileC != None:
            process2 = subprocess.Popen(['molti-console','-p',str(a),'-o',output + str(i/10) , fileA,fileB,fileC],
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
            stdout, stderr = process2.communicate()
            print(stdout.rstrip())
        elif fileB != None and fileC == None:
            process2 = subprocess.Popen(['molti-console','-p',str(a),'-o',output + str(i/10) ,fileA,fileB],
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
            stdout, stderr = process2.communicate()
            print(stdout.rstrip())
        else:
            process2 = subprocess.Popen(['molti-console','-p',str(a),'-o',output + str(i/10) , fileA],
                        stdout=subprocess.PIPE, 
                        stderr=subprocess.PIPE)
            stdout, stderr = process2.communicate()
            print(stdout.rstrip())
                
            
 

def id_to_cluster(molti_file):

    '''
    Parameters: file with communities returned by MolTi software.
    Returns: dictionary where keys represent community name and values represent nodes.
    '''

    File_molti= open(molti_file, 'r')
    llista_ids=[]
    Dict_clust={}
    Dict_zeros={}
    name=''
    ids=["1","10","101","102","103","104","106","107","108","109","110","111","112","113",
    "114","115","116","117","118","119","12","120","122","13","14","15","16","17","2","20",
    "200","201","202","204","206","207","208","209","21","22","24","25","27","28","29","3",
    "30","301","302","304","305","306","308","309","31","310","312","313","314","315","317",
    "318","32","320","321","322","323","325","326","33","34","35","37","38","39","4","40",
    "41","42","43","44","45","47","48","49","5","50","51","52","53","54","55","56","57","58",
    "59","6","60","61","62","63","64","65","66","67","68","71","73","74","75","78","79","8",
    "80","81","82","83","84","85","86","87","9","R100","R103","R113","R119","R126","R20","R22",
    "R27","R33","R34","R36","R38","R39","R44","R46","R48","R50","R52","R53","R56","R63","R66",
    "R67","R77","R78","R80","R84","R86","R87","R93","R94","R99"]
    c=0
    for line in File_molti:
        Line= line.rstrip()
        if Line.startswith('ClusterID'):
            name=Line
            Dict_clust[name]=[]
        elif name != '' and Line != '':
            Dict_clust[name].append(Line)
        else:
            pass
    return Dict_clust

def NodesSameCluster(file):

    '''
    Parameters: file with communities returned by MolTi software.
    Returns: a dictionary of lists where the keys represent nodes and the lists represent all the nodes that are in the same community.
    '''
    
    Dict_p={}
    value2=[]
    fun=id_to_cluster(file)
    for key,value in fun.items():
        a= [val for val in value]
        for element in a:
            value.remove(element)
            Dict_p[element]= value
            value=copy.deepcopy(a)
    return Dict_p

def PairNodesTogether(path_molti_files):
    
    ''' 
    Parameters: directory where the files returned by MolTi in different resolutions can be located.    
    Returns: dictionary where keys represent the resolution parameter and values the number of pairs
    that are kept in the same community.
    '''
    
    fun= ""
    c= 0
    Dict={}
    F= [f for f in os.listdir(path_molti_files) if re.match('Communities_all_p.*[^csv]$', f)]
    F.sort()
    for f in F:
        c = c + 1
        Dict[c]=0
        if fun == "":
            fun= NodesSameCluster(file= path_molti_files + '/' + f)
            Dict[c]= sum([len(d) for d in fun.values()])/2
        else: 
            fun2= NodesSameCluster(file= path_molti_files + '/' + f)
            for k,v in fun.items():
                inter= set(fun2[k]).intersection(set(v))
                Dict[c]= Dict[c] + len(inter)
                fun2[k]= list(inter)
            Dict[c]= Dict[c]/2
            inter=0
            fun= fun2
    return Dict


def purge(dir, pattern):
    
    '''
    PARAMETERS
    
    dir: directory where MOlti files are located
    pattern: regex pattern to match the files that will be kept
    '''
    #llista=[]
    for f in os.listdir(dir):
        if not re.search(pattern, f):
            #llista.append(float(f.split('_')[-1]))
            os.remove(os.path.join(dir, f))
    
    #print(sorted(llista,key=float))
    #llista.sort()
    #print(llista)

    
### CALL FUNCTION ###

#tuning_resolution(reso_min= 5, reso_max= 2000, fileA= '/home/nuria/Desktop/CV_filt/Backbones_molti/Backbone_closure_ultra_mrna.ncol', fileB= '/home/nuria/Desktop/CV_filt/Backbones_molti/Backbone_closure_ultra_methy.ncol', fileC= '/home/nuria/Desktop/CV_filt/Backbones_molti/Backbone_closure_ultra_micro_good.ncol', output= '/home/nuria/Desktop/Reviewers_comments/ResolutionTuning/Communities_all_p')
#tuning_resolution(reso_min= 1, reso_max= (int(op.reso_max))*10, fileA= op.ncol_a, fileB= op.ncol_b, fileC= op.ncol_c, output= op.out)
#test=pair_nodes('/home/nuria/Desktop/Reviewers_comments/ResolutionTuning/Communities_all_p0.5')
#test= PairNodesTogether(path_molti_files= '/home/nuria/Desktop/Reviewers_comments/ResolutionTuning')
#print([len(f) for f in test.values()])
#print(test)
#tuning_resolution(reso_min= 5, reso_max= op.reso_max, fileA= op.ncol_a, fileB= op.ncol_b, \ 
#fileC= op.ncol_c, output= op.out)
