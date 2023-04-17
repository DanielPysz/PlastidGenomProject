import pandas as pd
import numpy as np
import argparse
import os
from Bio import Entrez, SeqIO
Entrez.email = "danpysz@gmail.com"
parser = argparse.ArgumentParser(description='convert json to dataframe')

parser.add_argument('-f1', '--file_path1', type=str, dest="file1")
parser.add_argument("-f2", "--file_path2", type=str, dest="file2", required=False)
parser.add_argument("-f3", "--csv_Ola", type=str, dest="csv", required=False)
args = parser.parse_args()
refseq = ['NC_038009.1', 'NC_042170.1', 'NC_035260.1', 'NC_038000.1', 'NC_029860.1', 'NC_058737.1', 'NC_031145.1', 'NC_058771.1', 'NC_058768.1', 'NC_058770.1', 'NC_058769.1', 'NC_058772.1', 'NC_070145.1', 'NC_031172.1', 'NC_031148.1', 'NC_024080.1', 'NC_038008.1', 'NC_037999.1', 'NC_012898.1', 'NC_067219.1', 'NC_034636.1', 'NC_028029.1', 'NC_061049.1', 'NC_042171.1', 'NC_062392.1', 'NC_031173.1', 'NW_021703914.1', 'NC_038006.1', 'NC_038007.1', 'NC_008408.1', 'NC_034776.1', 'NC_035266.1', 'NC_035268.1', 'NC_035264.1', 'NC_035276.1', 'NC_034787.1', 'NC_021075.1', 'NC_035269.1', 'NC_035265.1', 'NC_035263.1', 'NC_047434.1', 'NC_031211.1', 'NC_031174.1', 'NC_025313.1', 'NC_063631.1', 'NC_053621.1', 'NC_025310.1', 'NC_040294.1', 'NC_065320.1', 'NC_066178.1', 'NC_057618.1', 'NC_020795.1', 'NC_058751.1', 'NC_026522.1', 'NC_014340.2', 'NC_014345.1', 'NC_062381.1', 'NC_035721.1', 'NC_062379.1', 'NC_036937.1', 'NC_069636.1', 'NC_046005.1', 'NC_035294.1', 'NC_032288.1', 'NC_030338.1', 'NC_071859.1', 'NC_035350.1', 'NC_025314.1', 'NC_042901.1', 'NC_041636.1', 'NC_059930.1', 'NC_024081.1', 'NC_059929.1', 'NC_028502.1', 'NC_058654.1', 'NC_058652.1', 'NC_058658.1', 'NC_058659.1', 'NC_058669.1', 'NC_027286.1', 'NC_035720.1', 'NC_013703.1', 'NC_057222.1', 'NC_051883.1', 'NC_004799.1', 'NC_038216.1', 'NC_001675.1', 'NC_038215.1', 'NC_028632.1', 'NC_024082.1', 'NC_031161.1', 'NC_035280.1', 'NC_035287.1', 'NC_070143.1', 'NC_020371.1', 'NC_043929.1', 'NC_035297.1', 'NC_024083.1', 'NC_035298.1', 'NC_035257.1', 'NC_035288.1', 'NC_039969.1', 'NC_038005.1', 'NC_067952.1', 'NC_062388.1', 'NC_014287.1', 'NC_070408.1', 'NC_004823.1', 'NC_007288.1', 'NC_038231.1', 'NC_038001.1', 'NC_062391.1', 'NC_031176.2', 'NC_057645.1', 'NC_071860.1', 'NC_035156.1', 'NC_038187.1', 'NC_071861.1', 'NC_001603.2', 'NC_039156.1', 'NC_002652.1', 'NC_071862.1', 'NC_020460.2', 'NC_024154.1', 'NC_027269.1', 'NC_024928.1', 'NC_040295.1', 'NC_040296.1', 'NC_040297.1', 'NC_040298.1', 'NC_066179.1', 'NC_066976.1', 'NC_017754.2', 'NC_065322.1', 'NC_056790.1', 'NC_015403.1', 'NC_071863.1', 'NC_044407.1', 'NC_045244.1', 'NC_016735.1', 'NC_024665.1', 'NC_041173.1', 'NC_029858.1', 'NC_040156.1', 'NC_041174.1', 'NC_040158.1', 'NC_041175.1', 'NC_029859.1', 'NC_063784.1', 'NC_063783.1', 'NC_063782.1', 'NC_063785.1', 'NC_040152.1', 'NC_062394.1', 'NC_040153.1', 'NC_067848.1', 'NC_062382.1', 'NC_065317.1', 'NC_058651.1', 'NC_058653.1', 'NC_058655.1', 'NC_058656.1', 'NC_038051.1', 'NC_058657.1', 'NC_058660.1', 'NC_046041.1', 'NC_039140.1', 'NC_058661.1', 'NC_039141.1', 'NC_058662.1', 'NC_058663.1', 'NC_058664.1', 'NC_058665.1', 'NC_058666.1', 'NC_023785.1', 'NC_046042.1', 'NC_046043.1', 'NC_058667.1', 'NC_039092.1', 'NC_031149.1', 'NC_038100.1', 'NC_029644.1', 'NC_039143.1', 'NC_039144.1', 'NC_058670.1', 'NC_037841.1', 'NC_058752.1', 'NC_021618.1', 'NC_066971.1', 'NC_035290.1', 'NC_000926.1', 'NC_069611.1', 'NC_037998.1', 'NC_029741.1', 'NC_044463.1', 'NC_044464.1', 'NC_044465.1', 'NC_046751.1', 'NC_044758.1', 'NC_065321.1', 'NC_066074.1', 'NC_067583.1', 'NC_044491.1', 'NC_061749.1', 'NC_067582.1', 'NC_067581.1', 'NC_035279.1', 'NC_010772.1', 'NC_031177.1', 'NC_031146.1', 'NC_039142.1', 'NC_058668.1', 'NC_066452.1', 'NC_058314.1', 'NC_049168.1', 'NC_035296.1', 'NC_014267.1', 'NC_035293.1', 'NC_031178.1', 'NC_044689.1', 'NC_057231.1', 'NC_044690.1', 'NC_035259.1', 'NC_027093.1', 'NC_039928.1', 'NC_039970.1', 'NC_039971.1', 'NC_039968.1', 'NC_024084.1', 'NC_035272.1', 'NC_056288.1', 'NC_044182.1', 'NC_022667.1', 'NC_066181.1', 'NC_067953.1', 'NC_035292.1', 'NC_029743.1', 'NC_058274.1', 'NC_065334.1', 'NC_062385.1', 'NC_040135.1', 'NC_031167.1', 'NC_035281.1', 'NC_039145.1', 'NC_032041.1', 'NC_032399.1', 'NC_032396.1', 'NC_065318.1', 'NC_022261.1', 'NC_067956.1', 'NC_020018.1', 'NC_027287.1', 'NC_022259.1', 'NC_022262.1', 'NC_022263.1', 'NC_022260.1', 'NC_066075.1', 'NC_056794.1', 'NC_039978.1', 'NC_021189.1', 'NC_024050.1', 'NC_007932.1', 'NC_066180.1', 'NC_044785.1', 'NC_065388.1', 'NC_061046.1', 'NC_056787.1', 'NC_061047.1', 'NC_067955.1', 'NC_057170.1', 'NC_035284.1', 'NC_035262.1', 'NC_051457.1', 'NC_046495.1', 'NC_031147.1', 'NC_069610.1', 'NC_029742.1', 'NC_011087.1', 'NC_039737.1', 'NC_049013.1', 'NC_035261.1', 'NC_071864.1', 'NC_039967.1', 'NC_039927.1', 'NC_016703.2', 'NC_008588.1', 'NC_037995.1', 'NC_037997.1', 'NC_023293.1', 'NW_017385168.1', 'NC_031963.1', 'NC_051845.1', 'NC_031401.1', 'NW_017962150.1', 'NC_031964.1', 'NW_023618544.1', 'NW_019223693.1', 'NC_035258.1', 'NC_032045.1', 'NC_066077.1', 'NC_042794.1', 'NC_031179.1', 'NC_035274.1', 'NC_035277.1', 'NC_035282.1', 'NC_035270.1', 'NC_035275.1', 'NC_038144.1', 'NC_000925.1', 'NC_035573.1', 'NC_062300.1', 'NC_023133.1', 'NC_031175.1', 'NC_056910.1', 'NC_061048.1', 'NC_038004.1', 'NC_040299.1', 'NC_038002.1', 'NC_058782.1', 'NC_058785.1', 'NC_058780.1', 'NC_058781.1', 'NC_027721.1', 'NC_058784.1', 'NC_058783.1', 'NC_058786.1', 'NC_062386.1', 'NC_044408.1', 'NC_056103.1', 'NC_063632.1', 'NC_062387.1', 'NC_029861.1', 'NC_043890.1', 'NC_038003.1', 'NC_025311.1', 'NC_062393.1', 'NC_062301.1', 'NC_035271.1', 'NC_070142.1', 'NC_009573.1', 'NC_062383.1', 'NC_031144.1', 'NC_025312.1', 'NC_018523.1', 'NC_049039.1', 'NC_071183.1', 'NC_062389.1', 'NC_066050.1', 'NC_066457.1', 'NC_048511.1', 'NC_029856.1', 'NC_066458.1', 'NC_064732.1', 'NC_064730.1', 'NC_066459.1', 'NC_029134.1', 'NC_031168.1', 'NC_053868.1', 'NC_031169.1', 'NC_057081.1', 'NC_046447.1', 'NC_031170.1', 'NC_035231.1', 'NC_061768.1', 'NC_058703.1', 'NC_058705.1', 'NC_058701.1', 'NC_058704.1', 'NC_058702.1', 'NC_035289.1', 'NC_029857.1', 'NC_035285.1', 'NC_069609.1', 'NC_071865.1', 'NC_062384.1', 'NC_035267.1', 'NC_039977.1', 'NC_066049.1', 'NC_040134.1', 'NC_035295.1', 'NC_070147.1', 'NC_027589.1', 'NC_067957.1', 'NC_060857.1', 'NC_060858.1', 'NC_060859.1', 'NC_058877.1', 'NC_058875.1', 'NC_008589.1', 'NC_058876.1', 'NC_058878.1', 'NC_035291.1', 'NC_007758.1', 'NC_031171.1', 'NC_035286.1', 'NC_062390.1', 'NC_040291.1', 'NC_035299.1', 'NC_031425.1', 'NC_001799.1', 'NC_027288.1', 'NC_026851.1', 'NC_001713.1', 'NC_027746.1', 'NC_056791.1', 'NC_062380.1', 'NC_016731.1', 'NC_028503.1', 'NC_065319.1', 'NC_011600.1', 'NC_035283.1', 'NC_035278.1', 'NC_026523.1', 'NC_035273.1', 'NC_065488.1', 'NC_040300.1', 'NC_029576.1', 'NC_070146.1']
def read(file):
    ind_ls = []
    df = pd.read_json(file)
    source = df["Source"]
    for i,j in enumerate(source):
        a = j.values()
        if " Viridiplantae; Streptophyta; Embryophyta; " not in str(a):
            ind_ls.append(i)
    return ind_ls, df 

def dane(ind_ls, df):
    index = []
    for i in df.index:
        if i not in ind_ls:
            df = df.drop(index = i)             #Filtrowanie gatunkow
    lst = list(dict(df["Definition"]).values()) #Lista z nazwami gatunkowymi wszystkich organizmow 
    df_n = list(dict(df["Origin"]).values())    #Lista z sekwencjami wszystkich gatunkow
    for i,j in enumerate(df_n):
        j = str(j).upper()
        df_n[i] = j
    for i in df.index:
        j = df["Locus"][i]
        index.append(j.split(" ")[0])
    return df_n, lst, index

def docelowy(df_n_1, df_n_2, lst_1, lst_2, index1, index2):
    lst = lst_1 + lst_2
    index1 = index1 + index2
    df_n_1.extend(df_n_2)
    dane = {"Numer":index1, "Organizm": lst, "Sekwencja": list(df_n_1)}
    baza = pd.DataFrame(data=dane)
    #print(baza)
    return baza

def dodanie_df(csv, baza):
    od_oli = pd.read_csv(csv)
    od_oli = od_oli.drop("Unnamed: 0", axis=1)
    cala = pd.concat([baza, od_oli], ignore_index=True)
    return cala

def tworzenie_pliku_genbank(numer):
    genbank_handle = Entrez.efetch(db="nucleotide", id=numer, rettype="gb", retmode="text")
    genbank_record = SeqIO.read(genbank_handle, "genbank")
    genome_seq = genbank_record.seq
    with open(numer + "sequence.gb", "w") as handle:
        SeqIO.write(genbank_record, handle, "genbank")

def translacja_calosc(seq_file):
    os.chdir("/mnt/archive/Cicuta_serwer/pysz/cpdatabase/genbank")
    record = SeqIO.read(seq_file, "genbank")
    my_cds = record.features
    for feature in my_cds:
        if feature.type == 'CDS':
            print(feature.qualifiers["product"])
            print(feature.qualifiers["translation"])

def translacja(feature):
    name = feature.qualifiers['product']
    protein = feature.qualifiers['translation']
    #return name, protein
    print(feature.qualifiers)

def pliki(cala):
    os.chdir("/mnt/archive/Cicuta_serwer/pysz/cpdatabase/genbank")
    numery = list(cala['Numer'])
    for numer in numery:
        print(f'Tworze plik {numer}')
        tworzenie_pliku_genbank(numer)

def bialka(sciezka):
    bledy = []
    os.chdir(sciezka)
    for i in os.listdir(sciezka):
        if i[-2:] == "gb":
            record = SeqIO.read(i, "genbank")
            for feature in record.features:
                if feature.type == 'CDS':
                    try:
                        organism = str(record.annotations["organism"])
                        nazwa = str(list(feature.qualifiers['product'])[0])
                        seq = list(feature.qualifiers['translation'])[0]
                        fasta = str(f'>{organism}_{nazwa}\n{seq}')
                        nazwa = nazwa.replace("/", "_")
                        if os.path.isdir(nazwa):
                            print(f'Katalog {nazwa} istnieje')
                            os.chdir(nazwa)
                            name = organism+"_"+nazwa
                            with open(name, "w") as writer:
                                writer.write(fasta)
                        else:
                            print(f'Tworze katalog {nazwa}')
                            os.makedirs(str(nazwa))
                            os.chdir(str(nazwa))
                            name = organism+"_"+nazwa
                            with open(name, "w") as writer:
                                writer.write(fasta)
                    except:
                        bledy.append(f'Blad w bialku {i} ')
                        
                    os.chdir(sciezka)
            print(f'Zakonczono dzialania z plikiem {i}')
    print(f'Zakonczono dzialania z wszystkimi plikami')
    print("Lista bledow: ")
    print(bledy)

#ind_ls1, df1 = read(args.file1)
#df_n_1, lst_1, index1 = dane(ind_ls1, df1)
#ind_ls2, df2 = read(args.file2)
#df_n_2, lst_2, index2 = dane(ind_ls2, df2)
#baza = docelowy(df_n_1, df_n_2, lst_1, lst_2, index1, index2)
#cala = dodanie_df("out.csv", baza)
#pliki(cala)
bialka('/mnt/archive/Cicuta_serwer/pysz/cpdatabase/genbank')
#translacja_cala("sequence.gb")
