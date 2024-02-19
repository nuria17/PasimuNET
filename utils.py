
from argparse import ArgumentParser

#args = ArgumentParser()
#args.add_argument('--corr', help='Input txt file with the correlation matrix', dest='corr',  default=None)
#args.add_argument('--backbone', help='Input txt file with the backbone edges', dest='backbone',  default=None)
#args.add_argument('--out', help='Output txt file with the network in the suitable ncol format for Molti software',dest='out', default=None)

#op = args.parse_args()

def create_dict_of_dict(file):
    Dict={}
    File=open(file, 'r')
    for line in File:
        NodeA= line.split(',')[0]
        NodeB= line.split(',')[1].rstrip()
        if NodeA not in Dict.keys():
            d= {NodeB: 1}
            Dict[NodeA]= d
        elif NodeA in Dict.keys():
            d[NodeB]= 1
            Dict[NodeA]= d
    return Dict

def ncol_format(file, backbone, cutoff2=0):
    file2= open(file, 'r')
    backbone2= open(backbone, 'r')
    Dict= {}
    Dict2={}
    i= 0
    j= 0
    names= list()

    for line in backbone2:
        Line= line.rstrip().split(',')
        if Line[0] not in Dict2.keys():
            Dict2[Line[0]]= [Line[1]]
        else:
            Dict2[Line[0]].append(Line[1])

    for line in file2:        
        j= j + 1   
        Line= line.rstrip().split(' ')
        Dict[Line[0]]= []
        
        for element in Line:
            i= i + 1
            if line.startswith(' '):
                names= Line               
            else:                
                if element != Line[0] and j != i and names[i-1] not in Dict.keys():
                    
                    if Line[0] in Dict2.keys():
                        if names[i-1] in Dict2[Line[0]]:
                            Dict[Line[0]].append([names[i-1], 1])
                        #elif float(element) < float(mean_node):
                        elif names[i-1] not in Dict2[Line[0]]:
                            Dict[Line[0]].append([names[i-1],0])
                            
                    else:
                        Dict[Line[0]].append([names[i-1],0])
                    #Dict[Line[0]].append([names[i-1],element])                

                #print(Dict.keys())                    
                #elif Dict[Line[0]].values() == '':
                    #Dict[Line[0]]= []
        i=0    
                                              
    return Dict        

#fun=(ncol_format(file= op.corr, backbone= op.backbone))
def write_ncol_file(file, backbone):
    
    '''
    File: File with correlation matrix in txt format.
    Backbone: File with the nodes separated by a comma.
    Out: File name with the backbone graph in ncol format.

    '''

    fun= ncol_format(file, backbone)
    t=""
    for key in fun:
        for tt in fun[key]:
            if t == '':
                t= str(key) + '\t' + str(tt[0]) + '\t' + str(tt[1]) + '\n' 
            else:
                t+= str(key) + '\t' + str(tt[0]) + '\t' + str(tt[1]) + '\n' 
    t=t.rstrip()
    return t
    #fd= open(op.out, 'w')
    #fd= open(out, 'w')
    #fd.write(t.rstrip())
    #fd.close()
