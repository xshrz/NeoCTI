import xlrd
import math
import numpy as np
import pandas as pd


#def Read_similarity():
#    PID = 'DTI_code\数据处理\pro_sim\string_protein_sequences.txt'
#    Protein = pd.read_table(PID, sep='\t', header=None)
#    Protein = Protein.values
#    p_txt = []
#    for index in range(int(len(Protein) / 2)):
#        p_index = []
#        Pname = Protein[index * 2][0][1:]
#        Pjene = Protein[index * 2 + 1][0]
#        p_index.append(Pname)
#        p_index.append(Pjene)
#        p_txt.append(p_index)
#    #print(p_txt)
#    return p_txt

'''读取化合物靶标关系数据，将其转换为矩阵形式保存'''
'''这里自己复现一遍了已经'''
def Read_chemical_target():
#    PID = 'DTI_code\AD_data\ch_target\CTD_Chemical-target.txt'
    PID = 'DTI_code\AD_data\ch_target\CTD_D000544.csv'
#    ct = pd.read_table(PID,header=None)
    ct = pd.read_csv(PID,header=None)
    ct = ct.values
    ch = []
    target = []

    #targetlist = ct[0][1].split('|')
    print(ct[24])
    print(ct[24][0])
    #print(targetlist)
    for index in range(len(ct)):
        if not(pd.isnull(ct[index][1])):
            ch.append(ct[index][0])
            targetlist = ct[index][1].split('|')
            for tar in targetlist:
                target.append(tar)
    chemical = sorted(set(ch), key=ch.index)  #去重且不改变顺序
    target = sorted(set(target), key=target.index)  # 去重且不改变顺序
    row = len(chemical)
    col = len(target)
    np.savetxt('DTI_code\AD_data\ch_target\ch_Name1.txt', chemical, fmt="%s", delimiter=',')
    np.savetxt('DTI_code\AD_data\ch_target\\target_Name1.txt', target, fmt="%s", delimiter=',')

    chemical_target = np.zeros((row, col))
    for index in range(len(ct)):
        if not(pd.isnull(ct[index][1])):
            row = chemical.index(ct[index][0])
            targetlist = ct[index][1].split('|')
            for tar in targetlist:
                col = target.index(tar)
                chemical_target[row][col] = 1

    return chemical_target


def ReadPPI100():
    PID = 'D:\DTI_code\数据处理\data_100pro\\target_Name.txt'       #读取靶标名称
    PName = pd.read_table(PID,header=None)
    PName = PName.values
    PN = []
    for p in PName:
        PN.append(p[0])
    #print(PN)
    InterctionID = 'D:\DTI_code\数据处理\data_100pro\PPI100.txt'         #读取靶标名称
    Interction = pd.read_table(InterctionID,header=None)
    Interction = Interction.values

    length = len(PName)
    PPI = np.zeros((length, length))

    for inter in Interction:
        row = -1
        col = -1
        for i in range(len(PN)):
            if PN[i] == inter[0]:
                row = i
            if PN[i] == inter[1]:
                col = i
        if(row != -1 and col != -1):
            PPI[row][col] = 1
            PPI[col][row] = 1


    np.savetxt('D:\DTI_code\数据处理\data_100pro\PPI_Matrix.csv', PPI, fmt="%d", delimiter=',')

'''读取化合物靶标关系数据，将其转换为矩阵形式保存'''
def ReadCTI100():
    PID = 'D:\DTI_code\数据处理\data_100pro\CTD_Chemical-target.txt'    #读取CTIs数据
    ct = pd.read_table(PID,header=None)
    ct = ct.values

    CID = 'D:\DTI_code\数据处理\data_100pro\ch_Name3401.txt'        #读取化合物名称
    ch_N = pd.read_table(CID,header=None)
    ch_N = ch_N.values
    ch_Name = []
    for c in ch_N:
        ch_Name.append(c[1])
    print(ch_Name)

    PID = 'D:\DTI_code\数据处理\data_100pro\\target_Name.txt'       #读取靶标名称
    PN = pd.read_table(PID,header=None)
    PN = PN.values
    PName = []
    for p in PN:
        PName.append(p[0])
    print(PName)

    CTIs = np.zeros((len(ch_Name), len(PName)))
    for row in range(len(ch_Name)):
        chemical = ch_Name[row]
        for index in range(len(ct)):
            if(chemical == ct[index][0]):                   #找到对应化合物那一行
                if not(pd.isnull(ct[index][1])):            #判断靶标是否为空
                    targetlist = ct[index][1].split('|')
                    col = -1
                    for tar in targetlist:
                        for i in range(len(PName)):
                            if PName[i] == tar:
                                col = i
                        if (col != -1):
                            CTIs[row][col] = 1

    np.savetxt('D:\DTI_code\数据处理\data_100pro\CTIs_Matrix.csv', CTIs, fmt="%d", delimiter=',')


def Read_similarity():
    PID = 'DTI_code\AD_data\pro_sim\string_protein_sequences.txt'
    Protein = pd.read_table(PID, sep='\t', header=None)
    Protein = Protein.values

    PID = 'DTI_code\AD_data\ch_target\\target_Name.txt'       #读取靶标名称
    PN = pd.read_table(PID,header=None)
    PN = PN.values
    PName = []
    for p in PN:
        PName.append(p[0])

    p_txt100 = []
    for p in PName:
        for index in range(int(len(Protein) / 2)):
            name = Protein[index * 2][0][1:]
            if p == name:
                list = []
                list.append(Protein[index * 2][0])
                list.append(Protein[index * 2][1])
                list.append(Protein[index * 2][2])
                p_txt100.append(list[0] + '\t' +  list[1] + '\t' + list[2])
                p_txt100.append(Protein[index * 2 + 1][0])
                break

    np.savetxt('DTI_code\AD_data\pro_sim\P_squence.txt', p_txt100, fmt="%s", delimiter=',')


def Read_Protein_sim():          #建立靶标相似性网络并归一化

    PID = 'D:\DTI_code\数据处理\data_100pro\Protein_Similarity89.txt'  # 读取靶标相似性数据
    P_sim = pd.read_table(PID, header=None)
    P_sim = P_sim.values

    CTIs = np.zeros((89, 89))
    scorelist = []
    for p in P_sim:
        p = p[0]
        p = p.split(' ')
        score = float(p[4])
        scorelist.append(score)

    max = 100
    min = 0
    print(min, max)

    for p in P_sim:
        p = p[0]
        p = p.split(' ')
        p_index = p[1][1:(len(p[1]) - 1)]
        p_index = p_index.split(':')
        row = int(p_index[0]) - 1
        col = int(p_index[1]) - 1
        score = (float(p[4]) - min) / (max - min)
        # print(row, col , score)
        CTIs[row][col] = score
        CTIs[col][row] = score
    for i in range(len(CTIs)):
        CTIs[i][i] = 1

    np.savetxt('D:\DTI_code\数据处理\data_100pro\protein_sim.csv', CTIs, fmt="%s", delimiter=',')


'''读取化合物靶标关系数据，将其转换为矩阵形式保存'''
def ReadCTIS2516():

    PID = 'E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\CTIs_txt.txt'       #读取CTIs关系表
    ct = pd.read_table(PID,header=None)
    ct = ct.values

    CID = 'E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\chem_name.csv'        #读取化合物名称
    ch_N = pd.read_table(CID,header=None)
    ch_N = ch_N.values
    ch_Name = []
    for c in ch_N:
        ch_Name.append(c[0])
    print(ch_Name)

    PID = 'E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\pro1885name.csv'       #读取靶标名称
    PN = pd.read_table(PID,header=None)
    PN = PN.values
    PName = []
    for p in PN:
        PName.append(p[0])
    print(PName)

    CTIs = np.zeros((len(ch_Name), len(PName)))
    for row in range(len(ch_Name)):
        chemical = ch_Name[row]
        for index in range(len(ct)):
            if(chemical == ct[index][0]):                   #找到对应化合物那一行
                if not(pd.isnull(ct[index][1])):            #判断靶标是否为空
                    targetlist = ct[index][1].split('|')
                    col = -1
                    for tar in targetlist:
                        for i in range(len(PName)):
                            if PName[i] == tar:
                                col = i
                        if (col != -1):
                            CTIs[row][col] = 1

    np.savetxt('E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\CTIs_Matrix_2516_1885.csv', CTIs, fmt="%d", delimiter=',')


def ReadPPI1908():
    PID = 'D:\DTI_code\数据处理\data_3401\protein_name.txt'       #读取靶标名称
    PName = pd.read_table(PID,header=None)
    PName = PName.values
    PN = []
    for p in PName:
        PN.append(p[0])

    #print(PN)
    InterctionID = 'D:\DTI_code\数据处理\data_3401\PPI.txt'         #读取PPI网络
    Interction = pd.read_table(InterctionID,header=None)
    Interction = Interction.values

    length = len(PN)
    PPI = np.zeros((length, length))

    for inter in Interction:
        row = -1
        col = -1
        for i in range(len(PN)):
            if PN[i] == inter[0]:
                row = i
            if PN[i] == inter[1]:
                col = i
        if(row != -1 and col != -1):
            PPI[row][col] = 1
            PPI[col][row] = 1


    np.savetxt('D:\DTI_code\数据处理\data_3401\Matrix\PPI_Matrix.csv', PPI, fmt="%d", delimiter=',')

'''处理的是从jp网站上下载下来的数据'''
def Read_Protein_sim():          #建立靶标相似性网络并归一化

    PID = 'DTI_code\AD_data\Module\psim_score.txt'  # 读取靶标相似性数据
    P_sim = pd.read_table(PID, header=None)
    P_sim = P_sim.values

    Pro_sim = np.zeros((1908, 1908))
    scorelist = []
    for p in P_sim:
        p = p[0]
        p = p.split(' ')
        score = float(p[4])
        scorelist.append(score)

    for p in P_sim:
        p = p[0]
        p = p.split(' ')
        p_index = p[1][1:(len(p[1]) - 1)]
        p_index = p_index.split(':')
        row = int(p_index[0]) - 1
        col = int(p_index[1]) - 1
        score = float(p[4]) / 100
        # print(row, col , score)
        Pro_sim[row][col] = score
        Pro_sim[col][row] = score
    for i in range(len(Pro_sim)):
        Pro_sim[i][i] = 1

    np.savetxt('DTI_code\AD_data\Module\protein_sim19082.csv', Pro_sim, fmt="%s", delimiter=',')


def readpairs():
    PID = 'E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\All_known_target_pairs.txt'       #读取CTIs关系表
    ct = pd.read_table(PID,header=None)
    ct = ct.values
    cti={}
    for i in ct:
        cti[str(i[0])]=str(i[1])
#Read_Protein_sim()
#ReadPPI1908()
#a=Read_chemical_target()
#print(a)
#b=Read_similarity()
#print(b)
#Read_Protein_sim()
    
PID = 'E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\All_known_target_pairs.txt'       #读取CTIs关系表
ct = pd.read_table(PID,header=None)
ct = ct.values
cti={}
for i in ct:
    cti[str(i[0])]=[]
for i in ct:
    cti[str(i[0])].append(str(i[1]))
cti_txt=[]
for key,value in cti.items():
    cti_txt.append([key,"|".join(value)])

np.savetxt('E:\\小崔专用\\大学文件_论文\\小NC\\laprls000\\input\\input\\CTIs_txt.txt', cti_txt, fmt="%s", delimiter='	')
ReadCTIS2516()