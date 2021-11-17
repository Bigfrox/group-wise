"""
Assignment 6, Node-Based Group-Wise Semantic Similarity Measure
2016253072 명수환

"""


from datetime import datetime
from os import name
import sys
import matplotlib.pyplot as plt



def getDataFromFile(filename):
    file = open(filename, 'r')
    

    ontology_MF = dict()
    ontology_BP = dict() 

    line = None
    line = file.readline()
    #print(line)
    while line != "":
        #print(line)
        while (line != "[Term]\n"):
            line = file.readline() #* Term 단위
            if line == "":
                return ontology_MF,ontology_BP

        line = file.readline().split()
        #print(line)
        #* Read line after "[Term]\n"
        if line[0] == "id:":
            id = line[1]
            #print("id:",id)
            line = file.readline().split()

        if line[0] == "name:":
            #print("name:",line[1:])
            line = file.readline().split()

        if line[0] == "namespace:":
            namespace = line[1]
            if namespace != "molecular_function" and namespace != "biological_process":
                line = file.readline().split()
                continue
            #print("namespace:",namespace)#* molecular_function,biological_process
            line = file.readline().split()

        
        is_a_relation = set()
        part_of_relation = set()
        is_obsleted = False
        while line:
            #print(line)
            if line[0] == "is_obsolete:":
                if line[1] == "true":
                    #print(id,"is obsolete")
                    line = file.readline().split()
                    is_obsleted = True
                    break
            elif line[0] == "is_a:":
                is_a_relation.add(line[1])
                line = file.readline().split()
            elif line[0] == "relationship:":
                if line[1] == "part_of":
                    part_of_relation.add(line[2]) #* id
                    line = file.readline().split()
                else:
                    line = file.readline().split()    
            else:
                line = file.readline().split()

        if is_obsleted:
            continue    

        error_relation = is_a_relation.intersection(part_of_relation)
        if error_relation:
            print(id,"has two relationship - ",error_relation)

        #* Classify - BP,MF
        if namespace == "molecular_function":
            ontology_MF[id] = list(is_a_relation.union(part_of_relation)) #* {id:relation}
        elif namespace == "biological_process":
            ontology_BP[id] = list(is_a_relation.union(part_of_relation))

        
        

    return ontology_MF,ontology_BP


def getTermDepth(ontology,start_node,root_node):
    shortestPathLength = -1 #* 초기값
    
    
    for v in ontology[start_node]:
        edges_count = 0
        parent = v
        #print("parent : ", parent)
        edges_count += 1
        if parent != root_node:
            start_node = parent
            edges_count += getTermDepth(ontology,start_node,root_node)
        else: #* parent is root
            return 1
        #print("short: ",shortestPathLength)
        #print("edge count : ",edges_count)
        if shortestPathLength != -1:
            if edges_count < shortestPathLength:
                shortestPathLength = edges_count
        elif shortestPathLength == -1:
            shortestPathLength = edges_count
        else: #* shortestPathLength < edges_count
            pass
        
    #print("short: ",shortestPathLength)
    return shortestPathLength
        
def output_to_file(filename,start_node,length):
    file = open(filename, 'a')
    

    file.write('{0} : {1}'.format(start_node,length))
    file.write("\n")
    file.close()
    
def getDataFromGAF(filename,BP,MF):
    
    file = open(filename, 'r')

    line = None
    
    while True:
        
        line = file.readline()

        if not line:
            break

        if line[0] == '!':
            continue

        line = line.split('\t')

        term = line[4].strip()
        gene = line[2].strip()
        code = line[6].strip()
        confirm = line[3].strip()
        ontology = line[8].strip()
        
        if gene == "":
            
            continue
        if term == "":
            
            continue
        if code == "":
            
            continue
        if code == "IEA":
            
            continue
        if confirm[:3] == "NOT":
            
            continue
        if ontology == "C":
            
            continue
        
        
        if ontology == "F":
            
            MF[term].add(gene)
            
        elif ontology == "P":
            BP[term].add(gene)
            
        
    return BP,MF
        

def GetPPI(filename,BP,MF,ontology_BP,ontology_MF):
    file = open(filename, 'r')
    line_count = 0
    total_line = 159567
    line = None
    gene_sim = list()
    while True:
        line_count += 1
        if line_count % int(total_line / 100) == 0:
            print("Progress : {0}%".format(round(line_count*100/total_line , 3)))
        line = file.readline().split()
        #print("line : ",line)
        if not line:
            break
        gene1 = line[0]
        gene2 = line[1]
        
        gene1_BP = dict()
        gene2_BP = dict()
        gene1_MF = dict()
        gene2_MF = dict()
        #* find all terms having gene1
        for term in BP:
            if gene1 in BP[term]:
                if not gene1 in gene1_BP:
                    gene1_BP[gene1] = set()

                gene1_BP[gene1].add(term)
                
            if gene2 in BP[term]:
                if not gene2 in gene2_BP:
                    gene2_BP[gene2] = set()
                gene2_BP[gene2].add(term)
                

        for term in MF:
            
            if gene1 in MF[term]:
                if not gene1 in gene1_MF:
                    gene1_MF[gene1] = set()
                gene1_MF[gene1].add(term)
                
            if gene2 in MF[term]:
                if not gene2 in gene2_MF:
                    gene2_MF[gene2] = set()
                gene2_MF[gene2].add(term)

        gene1_BP_ancestor = set()
        gene2_BP_ancestor = set()
        gene1_MF_ancestor = set()
        gene2_MF_ancestor = set()

        if gene1_BP:
            gene1_BP_ancestor = GetAllAncestor(gene1_BP[gene1],ontology_BP)
        if gene2_BP:
            gene2_BP_ancestor = GetAllAncestor(gene2_BP[gene2],ontology_BP)
        if gene1_MF:
            gene1_MF_ancestor = GetAllAncestor(gene1_MF[gene1],ontology_MF)
        if gene2_MF:
            gene2_MF_ancestor = GetAllAncestor(gene2_MF[gene2],ontology_MF)


        

        if gene1_BP_ancestor and gene1_MF_ancestor and gene2_BP_ancestor and gene2_MF_ancestor:
            #* Larger case within BP and MF
            #print("gene MF or BP")
            # print(gene1_BP_ancestor)
            # print(gene2_BP_ancestor)
            # print(gene1_MF_ancestor)
            # print(gene2_MF_ancestor)
            sim_BP = GetSimilarity(gene1_BP_ancestor,gene2_BP_ancestor)
            sim_MF = GetSimilarity(gene1_MF_ancestor,gene2_MF_ancestor)
            sim = max(sim_BP,sim_MF)
            tmp_list = [gene1,gene2,sim]
            gene_sim.append(tmp_list)
            #print("sim : ",sim)
            
        elif gene1_MF_ancestor and gene2_MF_ancestor:
            #* MF
            #print("gene MF")
            # print("gene1_MF_ancestor : ",gene1_MF_ancestor)
            # print("gene2_MF_ancestor : ",gene2_MF_ancestor)
            sim_MF = GetSimilarity(gene1_MF_ancestor,gene2_MF_ancestor)
            tmp_list = [gene1,gene2,sim_MF]
            gene_sim.append(tmp_list)
            #print("sim : ",sim_MF)

        elif gene1_BP_ancestor and gene2_BP_ancestor:
            #* BP
            sim_BP = GetSimilarity(gene1_BP_ancestor,gene2_BP_ancestor)
            tmp_list = [gene1,gene2,sim_BP]
            gene_sim.append(tmp_list)
            #print("sim : ",sim_BP)
        else:
            # * BP and MF 
            # * or MF and BP
            #print("NOT Case!")
            continue

        #print("gene_sim = ", gene_sim)
            
        
    return gene_sim
        
        
        
def GetAllAncestor(nodeset,ontology):
    if not nodeset:
        return set()
    newset = set()
    
    for node in nodeset:
        
        if ontology[node]:
            for v in ontology[node]:
                newset.add(v)
    
    newset = newset.union(GetAllAncestor(newset,ontology))

    return newset


def GetSimilarity(gene1_ancestor, gene2_ancestor):
    union_set = gene1_ancestor.union(gene2_ancestor)
    
    inter_set = gene1_ancestor.intersection(gene2_ancestor)
    # print("union : ",union_set)
    # print("inter : ",inter_set)
    return float(len(inter_set)) / float(len(union_set))


def output_to_file(filename,gene_sim):
    file = open(filename, 'w')
    
    for v in gene_sim:
        file.write('{0} {1} : {2}'.format(v[0] , v[1], v[2]))
        file.write("\n")
    file.close()


def main():


    if len(sys.argv) != 4:
        print("No input file.")
        print("<Usage> Assignment6.py go.obo goa_human.gaf biogrid_human_ppi_cln.txt")
        return -1;
    
    input_filename = sys.argv[1]
    input_filename2 = sys.argv[2]
    input_filename3 = sys.argv[3]
    output_filename = "output.txt"
    
    #input_filename = "go.obo"
    # input_filename2 = "goa_human.gaf"
    # input_filename3 = "biogrid_human_ppi_cln.txt"


    start_time = datetime.now()
    
    ontology_MF,ontology_BP = getDataFromFile(input_filename)


    
    
    
    print("length of MF :",len(ontology_MF))
    print("length of BP :",len(ontology_BP))

    for v1 in ontology_MF:
        
        if ontology_MF[v1] == []:
            root_MF = v1
            print("root node in MF : ",v1)

    for v2 in ontology_BP:
        
        if ontology_BP[v2] == []:
            root_BP = v2
            print("root node in BP : ",v2)

    error_count = 0
    for id_bp in ontology_BP:
        
        if id_bp == root_BP:
            continue
        for u1 in ontology_BP[id_bp][:]:
            
            if u1 in ontology_MF: # * ERROR -  Relation u and  v :  u is in MF, v is in BP
                error_count +=1
                ontology_BP[id_bp].remove(u1)
                # print(id_bp,"에서 ",u1,"를 삭제하였습니다.")
                

    for id_mf in ontology_MF:
        
        if id_mf == root_MF:
            continue
        for u2 in ontology_MF[id_mf][:]:
            if u2 in ontology_BP: # * ERROR - Relation u and  v  :  u is in BP, v is in MF
                error_count +=1
                ontology_MF[id_mf].remove(u2)
                
                #print(id_mf,"에서 ",u2,"를 삭제하였습니다.")
    print("error count : ",error_count)
    

    
    BP_annotation = ontology_BP.copy()
    MF_annotation = ontology_MF.copy()
    for v in BP_annotation:
        BP_annotation[v] = set()
    for v in MF_annotation:
        MF_annotation[v] = set()    
    
    BP_annotation,MF_annotation = getDataFromGAF(input_filename2,BP_annotation,MF_annotation)
    
    #* Show label
    print("=========================================================================")
    # for v in BP_annotation:
    #     if BP_annotation[v]:
    #         #print(v,BP_annotation[v])
    #         #print("\n")
    #         break
    # for v in MF_annotation:
    #     if MF_annotation[v]:
    #         print(v, MF_annotation[v])
    #         break




    #* inferred
    # for child in BP_annotation:
    #     for parent in ontology_BP[child]:
    #         = BP_annotation[parent].union(BP_annotation[child])
            
    # for child in MF_annotation:
    #     for parent in ontology_MF[child]:
    #        = MF_annotation[parent].union(MF_annotation[child])

    gene_sim = GetPPI(input_filename3,BP_annotation,MF_annotation,ontology_BP,ontology_MF)
    # for v in gene_sim:
    #     print(v)
    #     print("\n")
    similarity = [0 for i in range(11)]
    for v in gene_sim:
        #*print(int(v[2]*10), v)
        #similarity[int(v[2]*10)] += 1
        similarity.append(v[2])
    print("Similarity 0.0 ~ 1.0 : ",similarity)
    plt.hist(similarity)
    plt.show()
    print("[+] Time Elapsed : ", datetime.now() - start_time, "microseconds")
    output_to_file(output_filename,gene_sim)


if __name__ == '__main__':
    main()