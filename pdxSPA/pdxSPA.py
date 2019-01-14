import pandas as pd
import warnings
warnings.filterwarnings("ignore")

#pre-treatment
def filePreTreat(df):
    #filter reverse and contaminant
    df = df[(df["Reverse"] != "+") & (df["Potential contaminant"] != "+")]
    
    df["PEP"].fillna(1000, inplace=True)
    df = df[df["PEP"] <= 0.01]
    
    df = df.sort_values("Score", ascending=0)    
    df.drop_duplicates(["Raw file", "Fraction", "Experiment", "MS/MS Scan Number"], \
                       inplace = True)

    return df    

#sort peptide into 3 class    
def sortPeptides(row):
    row = row.lower()
    if "human" in row:    
        if "mouse" in row:
            return "SH"
        else:
            return "HU"
    elif "mouse" in row:
        return "MU"
    else:
        return "UNKNOWN"

#calculate total R    
def calculateR(df_HU, df_MU):
    df_HU_indensity = sum(df_HU["Intensity"].fillna(0))
    df_MU_indensity = sum(df_MU["Intensity"].fillna(0))
    return df_HU_indensity/df_MU_indensity

#get human protein ID
def getHumanProtein(protein_group_match):
    human_protein = {}
    for groupID in protein_group_match:
        if protein_group_match[groupID].first_human != '':
            human_protein[groupID] = protein_group_match[groupID].first_human
    return human_protein

#get all human and mouse protein entry
def getAllProteinEntry(protein_group_match):
    protein_entry_all = []
    for groupID in protein_group_match:
        if protein_group_match[groupID].first_human != '':
            protein_entry = protein_group_match[groupID].first_human.split("|")[-1].split("_")[0]
            protein_entry_all.append(protein_entry)
            
        if protein_group_match[groupID].first_mouse != '':
            protein_entry = protein_group_match[groupID].first_mouse.split("|")[-1].split("_")[0]
            protein_entry_all.append(protein_entry)
            
    return list(set(protein_entry_all))


class proteinGroupMatch:
    def __init__(self, first_human, first_mouse):
        self.first_human = first_human
        self.first_mouse = first_mouse
        
def getProteinGroupMatch(df):
    protein_group_match = {}
    groupID_exsit = []
    
    for index, row in df.iterrows():
        proteins = row["Proteins"].split(";")
        leading_proteins = row["Leading proteins"].split(";")
        groupIDs = row["Protein group IDs"].split(";")

        i = 0
        for i in range(len(groupIDs)-1):
            groupID = groupIDs[i]
            if groupID in groupID_exsit:
                continue
            
            index_in_proteins_this_group = proteins.index(leading_proteins[i])
            index_in_proteins_next_group = proteins.index(leading_proteins[i+1])
            
            first_human, first_mouse = "", ""
            flag_human, flag_mouse = 0, 0

            for j in range(index_in_proteins_this_group, index_in_proteins_next_group):
                if (flag_mouse == 0) and ("mouse" in proteins[j]):
                    first_mouse = proteins[j]
                    flag_mouse = 1

                elif (flag_human == 0) and ("human" in proteins[j]):
                    first_human = proteins[j]
                    flag_human = 1            
                
                if flag_mouse == 1 and flag_human == 1:
                    break
                
            groupID_exsit.append(groupID)
            protein_group_match[groupID] = proteinGroupMatch(first_human, first_mouse)
        
        if len(groupIDs) != 1:
            i += 1                
        index_in_proteins_this_group = proteins.index(leading_proteins[i])
        index_in_proteins_next_group = len(proteins)
        
        groupID = groupIDs[i]
        if groupID in groupID_exsit:
            continue
        first_human, first_mouse = "", ""
        flag_human, flag_mouse = 0, 0
        for j in range(index_in_proteins_this_group, index_in_proteins_next_group):
            if (flag_mouse == 0) and ("mouse" in proteins[j].lower()):
                first_mouse = proteins[j]
                flag_mouse = 1

            elif (flag_human == 0) and ("human" in proteins[j].lower()):
                first_human = proteins[j]
                flag_human = 1            
            
            if flag_mouse == 1 and flag_human == 1:
                break  
        groupID_exsit.append(groupID)        
        protein_group_match[groupID] = proteinGroupMatch(first_human, first_mouse)
    
    return protein_group_match

#get human protein expression by SPA
def getHumanProteinExpressionSPA(df, human_protein):
    df_SPA_human_expression_all = pd.DataFrame({"Protein":list(human_protein.values())})
    df_SPA_summary_all = pd.DataFrame({"Summary":["HU Peptide Count:", "MU Peptide Count:", "SH Peptide Count:","Total Ratio"]})
    
    #calculate each experiment    
    for experiment_name, dfp in df.groupby("Experiment"):
        
        #calculate totalR
        r = calculateR(dfp[dfp["flag"] == "HU"], dfp[dfp["flag"] == "MU"])
        df_SPA_summary = getSummary(experiment_name, dfp,r)       
        df_SPA_summary_all = pd.concat([df_SPA_summary_all, df_SPA_summary], axis = 1)
        
        #filter out MU
        dfp = dfp[(dfp["flag"] == "HU") | (dfp["flag"] == "SH")]       
        
        #HU and SH intensity for all human proteins
        human_HU_intensity , human_SH_intensity = {}, {}
        for item in human_protein:
            human_HU_intensity[item], human_SH_intensity[item] = 0, 0
        
        dfp["Intensity"].fillna(0, inplace=True)
        
        #sum up HU and SH intensity of each human protein       
        for index, row in dfp.iterrows():
            for item in row["Protein group IDs"].split(";"): 
                if item in human_protein:
                    if row["flag"] == "HU":
                        human_HU_intensity[item] += row["Intensity"]
                    else:
                        human_SH_intensity[item] += row["Intensity"]

        #calculate human protein expression                        
        SPA_human_expression = []
        for item in human_protein:
            SPA_human_expression.append(human_HU_intensity[item] + r/(1+r)*human_SH_intensity[item])

        #make dataframe and cat all experiments together        
        df_SPA_human_expression = pd.DataFrame({experiment_name:SPA_human_expression})        
        df_SPA_human_expression_all = pd.concat([df_SPA_human_expression_all, df_SPA_human_expression], axis=1)

    return df_SPA_human_expression_all, df_SPA_summary_all

#get human protein expression by X-SPA or S-SPA
def getHumanProteinExpressionXSSPA(df, protein_group_match, protein_entry_all, human_protein, XSflag):
    df_SPA_human_expression_all = pd.DataFrame({"Protein":list(human_protein.values())})
    df_SPA_summary_all = pd.DataFrame({"Summary":["HU Peptide Count:", "MU Peptide Count:", "SH Peptide Count:", "Total Ratio"]})
    
    #calculate each experiment
    for experiment_name, dfp in df.groupby("Experiment"):
        
        #calculate totalR
        total_R = calculateR(dfp[dfp["flag"] == "HU"], dfp[dfp["flag"] == "MU"])
        df_SPA_summary = getSummary(experiment_name, dfp, total_R)       
        df_SPA_summary_all = pd.concat([df_SPA_summary_all, df_SPA_summary], axis = 1)
             
        #HU and MU intensity for each protein entry
        protein_HU_intensity, protein_MU_intensity= {}, {}
        for item in protein_entry_all:
            protein_HU_intensity[item], protein_MU_intensity[item] = 0, 0
 
        #HU and SH intensity for all human proteins       
        human_HU_intensity , human_SH_intensity = {}, {}
        for item in human_protein:
            human_HU_intensity[item], human_SH_intensity[item] = 0, 0
            
        dfp["Intensity"].fillna(0, inplace=True)
         
        for index, row in dfp.iterrows():
            for item in row["Protein group IDs"].split(";"):
                first_human = protein_group_match[item].first_human
                first_mouse = protein_group_match[item].first_mouse
                
                #sum up HU and MU intensity of each protein entry
                if first_human != "":
                    protein_entry = first_human.split("|")[-1].split("_")[0]
                    if row["flag"] == "HU":
                        protein_HU_intensity[protein_entry] += row["Intensity"]
                    elif row["flag"] == "MU":
                        protein_MU_intensity[protein_entry] += row["Intensity"]

                
                if first_mouse != "":
                    protein_entry = first_mouse.split("|")[-1].split("_")[0]
                    if row["flag"] == "HU":
                        protein_HU_intensity[protein_entry] += row["Intensity"]
                    elif row["flag"] == "MU":
                        protein_MU_intensity[protein_entry] += row["Intensity"]
                 
                                    
                #sum up HU and SH intensity of each human protein
                if item in human_protein:
                    if row["flag"] == "HU":
                        human_HU_intensity[item] += row["Intensity"]
                    elif row["flag"] == "SH":
                        human_SH_intensity[item] += row["Intensity"]

        #calulate r for each protein entry        
        protein_entry_r = {}
        for item in protein_HU_intensity:
            if (protein_HU_intensity[item] != 0) and (protein_MU_intensity[item] != 0):
                protein_entry_r[item] = protein_HU_intensity[item]/protein_MU_intensity[item]
        
        #calculate human protein expression
        SPA_human_expression = []
        for item in human_protein:
            protein_entry = human_protein[item].split("|")[-1].split("_")[0]
            if protein_entry in protein_entry_r:   #if entry has r
                r = protein_entry_r[protein_entry]
                SPA_human_expression.append(human_HU_intensity[item] + r/(1+r)*human_SH_intensity[item])
            else: #if not have r, X:HU_intensity; S:share by total r
                if XSflag == "X":
                    SPA_human_expression.append(human_HU_intensity[item])
                else:
                    r = total_R
                    SPA_human_expression.append(human_HU_intensity[item] + r/(1+r)*human_SH_intensity[item])
        
        #make dataframe and cat all experiments together        
        df_SPA_human_expression = pd.DataFrame({experiment_name:SPA_human_expression})
        df_SPA_human_expression_all = pd.concat([df_SPA_human_expression_all, df_SPA_human_expression], axis=1)
    
    return df_SPA_human_expression_all, df_SPA_summary_all

#get summary for each experiment
def getSummary(experiment_name, df,r):    
    HU_peptide_count = df[df["flag"] == "HU"].shape[0]
    MU_peptide_count = df[df["flag"] == "MU"].shape[0]
    SH_peptide_count = df[df["flag"] == "SH"].shape[0]
    
    sample_summary = {experiment_name:[HU_peptide_count, MU_peptide_count, SH_peptide_count, r]}
    return pd.DataFrame(sample_summary)

#write summary for all experiments
def writeSummary(df_SPA_summary_all, ofile):
    infile1 = open(ofile , "w")
    string = ""
    for item in df_SPA_summary_all.columns:
        string += "\t" + item
    infile1.write(string + "\n")
    
    for index ,row in df_SPA_summary_all.iterrows():
        if index < 3:
            string = row["Summary"]
            for item in list(row.iloc[1:]):
                string += "\t" + str(int(item))
            infile1.write(string + "\n")
        else:
            string = row["Summary"]
            for item in list(row.iloc[1:]):
                string += "\t" + str(float(item))
            infile1.write(string + "\n")
    infile1.close()

   
def SPA(input_file, outdir):
    df = pd.read_table(input_file)  
    df = filePreTreat(df)    
    df["flag"] = df["Proteins"].apply(sortPeptides)
    protein_group_match = getProteinGroupMatch(df)      
    human_protein = getHumanProtein(protein_group_match)
    
    df_SPA_human_expression, df_SPA_summary_all = getHumanProteinExpressionSPA(df, human_protein)
    df_SPA_human_expression.to_csv(outdir + "SPA_Human_Protein_Expression.txt", sep="\t", index=0)    
    writeSummary(df_SPA_summary_all, outdir + "SPA_Summary.txt")


def X_SPA(input_file, outdir):   
    df = pd.read_table(input_file)   
    df = filePreTreat(df)    
    df["flag"] = df["Proteins"].apply(sortPeptides)
    protein_group_match = getProteinGroupMatch(df)
    human_protein = getHumanProtein(protein_group_match)     
    protein_entry_all = getAllProteinEntry(protein_group_match)
    
    df_XSPA_human_expression, df_XSPA_summary_all = getHumanProteinExpressionXSSPA(df, protein_group_match, protein_entry_all, human_protein, "X")
    df_XSPA_human_expression.to_csv(outdir + "XSPA_Human_Protein_Expression.txt", sep="\t", index=0)    
    writeSummary(df_XSPA_summary_all, outdir + "XSPA_Summary.txt")       

def S_SPA(input_file, outdir):
    df = pd.read_table(input_file)   
    df = filePreTreat(df)    
    df["flag"] = df["Proteins"].apply(sortPeptides)
    protein_group_match = getProteinGroupMatch(df)
    human_protein = getHumanProtein(protein_group_match)     
    protein_entry_all = getAllProteinEntry(protein_group_match)
    
    df_SSPA_human_expression, df_SSPA_summary_all = getHumanProteinExpressionXSSPA(df, protein_group_match, protein_entry_all, human_protein, "S")
    df_SSPA_human_expression.to_csv(outdir + "SSPA_Human_Protein_Expression.txt", sep="\t", index=0)    
    writeSummary(df_XSPA_summary_all, outdir + "SSPA_Summary.txt")   
    
    
    
    
    
    
