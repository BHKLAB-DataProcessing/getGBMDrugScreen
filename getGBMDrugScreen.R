library(stringr)
library(openxlsx)
library(abind)
library(PharmacoGx)

input_dir <- "/pfs/downloadGBMData/"
input_annotation <- "/pfs/annotation/"
out_dir <- "/pfs/out/" 
functions <- "/pfs/getGBMCellData/functions.R"

# input_dir <- "~/Documents/pfs/downloadGBMData/"
# input_annotation <- "~/Documents/pfs/annotation/"
# out_dir <- "~/Documents/pfs/getGBMDrugScreenObjects/"
# functions <- "./functions.R"

source(functions)

phen_exp <- readRDS("/pfs/getGBMGeneExpression/phen_exp.rds")
phen_cnv <- readRDS("/pfs/getGBMCNV/phen_cnv.rds")
phen_mutation <- readRDS("/pfs/getGBMMutation/phen_mutation.rds")
phen_methyl <- readRDS("/pfs/getGBMMethylation/phen_methyl.rds")

# phen_exp <- readRDS("~/Documents/pfs/getGBMGeneExpression/phen_exp.rds")
# phen_cnv <- readRDS("~/Documents/pfs/getGBMCNV/phen_cnv.rds")
# phen_mutation <- readRDS("~/Documents/pfs/getGBMMutation/phen_mutation.rds")
# phen_methyl <- readRDS("~/Documents/pfs/getGBMMethylation/phen_methyl.rds")

screen_objects<-function(screen, cell_obj, drug_obj, drug_with_ids_gbm, d_c_l){
  
  #Keeping the rows that are not all NaNs (2 represents the two first columns which are drugs and pat names).
  screen<-screen[rowSums(is.na(screen[,3:ncol(screen)])) <ncol(screen)-2, ] 
  
  #Multiplying the viabilites by 100
  screen[,3:ncol(screen)]<-apply(screen[,3:ncol(screen)],2,function(x){x*100})
  
  screen<-drug_name_correction(table = screen , column = "Drug")
  temp<-drug_name_correction(table = drug_with_ids_gbm, column="GBM.drugid")
  temp<-merge(temp[,c("unique.drugid","clean.ids")], screen, by="clean.ids", all.y=T)
  #temp$unique.drugid <- ifelse(is.na(temp$unique.drugid) & grepl("Doxorubicin", temp$Drug), temp$Drug , temp$unique.drugid)
  screen<-data.frame(temp[-c(1,3)])
  rm(temp)
  
  screen$EXP_details<-paste(screen$unique.drugid, screen$Pat , sep ="_")
  colnames(screen) <- c("drugid","cellid", colnames(screen)[3:length(colnames(screen))])# With this, drugids in the screen are already unique.ids
  rownames(screen)<-screen$EXP_details
  
  
  sen_info_scr<-screen[,c(1,2,14)]
  sen_info_scr$Duration<-"3 Days"
  
  ###sensitivity_raw
  dose_scr<-as.data.frame(matrix(nrow=nrow(screen)))
  i=1
  for (col in colnames(screen)[3:13]){
    micro_mol<-as.numeric(substr(col, start = 6, stop = nchar(col)))
    name<-paste("doses",i,sep="") #doses need to be named as (doses1, doses2,..)
    dose_scr[,name]<-micro_mol
    i=i+1}
  
  dose_scr<-dose_scr[-1]
  rownames(dose_scr)<-rownames(screen)
  viability_scr<-screen[,3:13] #Only keeping the doses 
  colnames(viability_scr)<-colnames(dose_scr) 
  
  sen_raw_scr<-abind(dose_scr, viability_scr, along = 3) #Creating a 3d array
  dimnames(sen_raw_scr)[[3]] <- c("Dose","Viability")
  
  ###Sensitivity_profile
  print("Running calculateFromRaw function... This takes a while")
  AAC_IC50_scr<-PharmacoGx:::.calculateFromRaw(sen_raw_scr, cap = NA, nthread = 1, family = c("normal","Cauchy"), scale = 0.07, n=1 )
  sen_profile_scr<-cbind(as.data.frame(AAC_IC50_scr["AUC"]), as.data.frame(AAC_IC50_scr["IC50"]))
  colnames(sen_profile_scr)<-c("aac_recomputed","ic50_recomputed")
  
  sen_profile_scr$auc_published<-d_c_l$value[match(rownames(sen_profile_scr),d_c_l$EXP_details)]
  sen_profile_scr$auc_published<-sen_profile_scr$auc_published
  
  ########### Adapting the cell object for each screen
  scr_cell_obj<-cell_obj
  diff_scr<-unique(screen$cellid[screen$cellid %in% scr_cell_obj$cellid ==FALSE])
  
  #For i in diff_scr, adds a new row 
  for(i in 1:length(diff_scr)){
    cid=substring(diff_scr[i],1,5) 
    if (cid %in% scr_cell_obj$cellid){
      scr_cell_obj[nrow(scr_cell_obj)+1,]<- scr_cell_obj[scr_cell_obj$cellid==cid,][1,]
      scr_cell_obj[nrow(scr_cell_obj),2]<- substring(diff_scr[i],6) #Replicate name
      scr_cell_obj[nrow(scr_cell_obj),4]<- diff_scr[i] #cell id
    }
    else{scr_cell_obj[nrow(scr_cell_obj)+1,]<- c(NA, NA,NA, diff_scr[i], rep(NA,length(colnames(scr_cell_obj))-4))}
  }
  
  rownames(scr_cell_obj)<-scr_cell_obj$cellid
  scr_cell_obj$tissueid <- "GBM"  # Adding mandatory tissue Id column 
  scr_cell_obj[,2]<-gsub("_","",scr_cell_obj[,2],fixed = T)# Removing "_"
  
  ########## Adapting the drug object for each screen
  scr_drugs <- drug_obj[unique(screen$drugid),]
  colnames(scr_drugs)[colnames(scr_drugs) == "unique.drugid"] = "drugid"
  
  ########## Curation Cell dataframe
  scr_cur_cell <- data.frame(unique.cellid=rownames(scr_cell_obj),
                             GBM.cellid=rownames(scr_cell_obj),
                             row.names=rownames(scr_cell_obj))
  
  ########## Curation drug dataframe
  scr_cur_drug <- data.frame(drugid = scr_drugs$drugid,
                             GBM.drugid = scr_drugs$GBM.drugid,
                             row.names = scr_drugs$drugid)
  
  ########## Curation tissue dataframe
  scr_cur_tissue <- data.frame(data.frame(tissueid = rep("GBM", nrow(scr_cell_obj)),
                                          GBM.tissueid = rep("GBM", nrow(scr_cell_obj)),
                                          row.names = scr_cell_obj$cellid))
  
  
  return(list("sen_info"= sen_info_scr, "sen_raw"= sen_raw_scr, "sen_profile"= sen_profile_scr,
              "cell_obj"= scr_cell_obj, "drug_obj"= scr_drugs,
              "cur_cell"= scr_cur_cell, "cur_drug"= scr_cur_drug,"cur_tissue"= scr_cur_tissue))
}

# ======================== Cell object ========================
#Gathering cell lines from all the experiments
phen_exp<-phen_exp[,colnames(phen_cnv)] # Reordeing the columns of phen_exp for rbind
all_cell_obj<-as.data.frame(unique(rbindlist( list(phen_exp,phen_cnv,phen_mutation,phen_methyl))))
rownames(all_cell_obj)<-all_cell_obj$cellid
print("Cell object: done")
# ======================== Drug annotation data from Pachy annotation ========================
drug_with_ids_gbm <- read.csv(paste(input_annotation, "drugs_with_ids.csv", sep = "") , stringsAsFactors = FALSE, na.strings = "") #this is the drug with ids file downloaded form pachy annotation repo

#To be corrected in "drug_with_ids.csv"
conc.name <- "Doxorubicin1///Doxorubicin2///Doxorubicin3///Doxorubicin4///Doxorubicin5///Doxorubicin6///Doxorubicin7///Doxorubicin8///Doxorubicin hydrochloride"
ind= which(drug_with_ids_gbm$unique.drugid == "Doxorubicin")
drug_with_ids_gbm$GBM.drugid[ind] <- conc.name # TO BE DELETED ONCE CORRECTED
drug_with_ids_gbm$GBM.drugid[which(drug_with_ids_gbm$unique.drugid == "Carfilzomib (PR-171) (combination with valproic acid)")] <- NA # TO BE DELETED ONCE CORRECTED

names <- unlist(strsplit(conc.name, "///"))

for(j in 1:length(names)){
  drug_with_ids_gbm[1+nrow(drug_with_ids_gbm) , ] = drug_with_ids_gbm[ind, ]
  drug_with_ids_gbm$GBM.drugid[nrow(drug_with_ids_gbm)] = names[j]
}

drug_with_ids_gbm <- drug_with_ids_gbm[-ind,]
drug_with_ids_gbm$unique.drugid <- ifelse(grepl("Doxorubicin", drug_with_ids_gbm$GBM.drugid), drug_with_ids_gbm$GBM.drugid, drug_with_ids_gbm$unique.drugid) 

print("Drug annotation: done")
# ======================== Drug object ========================
drugs<- read.xlsx(paste(input_dir, "mmc3.xlsx", sep="") ,rowNames = TRUE , startRow = 2)
# drugs[duplicated(drugs[ , "Compound.name"]),] #Checking duplications in drug names: 3 duplications found

#Replacing duplicated drug names in "Compound.name" column with the correct drug names:
drugs$Compound.name[drugs$Molecular.Formula=="C26H27ClN2O"]<-"LOFEPRAMINE"
drugs$Compound.name[drugs$`Chemical.Abstracts.Service.(CAS).code`=="101477-54-7"]<-"Lomerizine 2hcl"
drugs$Compound.name [drugs$ Compound.name == "lomerizine"] <-"Mitoxantrone dihydrochloride" 

#To make sure that the drug-object contains ALL the drug names used in screen2 and screen3 
#the drug-names in drug-object are mapped to drug-names from screen2 and screen3.

drugs$Compound.name[drugs$Compound.name == "bortezomib [velcade ]"]<-"Bortezomib  "
drugs$Compound.name[drugs$Compound.name == "5-azacytidine"]<-"Azacytidine-5"
drugs$Compound.name[drugs$Compound.name == "TV-001"]<-"prm-116"
drugs$Compound.name[drugs$Compound.name == "TV-002"]<-"prm-122"
drugs$Compound.name[drugs$Compound.name == "TV-003"]<-"PRM-123"
drugs$Compound.name[drugs$Compound.name == "tenovin-6 derivative 39"] <- "Tenovin-6 derivat A"#To be double checked
drugs$Compound.name[drugs$Compound.name == "tenovin-6 derivative 50"] <- "Tenovin-6 derivat B"#To be double checked 
drugs$Compound.name[drugs$Compound.name == "GLN-1001"]<-"Vakuinol-1" # Based on molecular formula (ClC1=CC=C(C2=NC(C=CC=C3)=C3C(C(O)C4NCCCC4)=C2)C=C1), The correct dictation is "Vacquinol-1"
drugs$Compound.name[drugs$Compound.name == "AM404"]<-"AMA404"
drugs$Compound.name[drugs$Compound.name == "dichlorbenzamide"]<-"Dichlorbenzamil"#To be double checked 
drugs$Compound.name[drugs$Compound.name == "8-azaguanine"]<-"azaguanine-8"
drugs$Compound.name[drugs$Compound.name == "duloxetine"]<-"(R)-Duloxetine" # this is mapped because we know drugs from "drug"file and "dose-resp" file must be the same
drugs$Compound.name[drugs$Compound.name == "propranolol-(S)"]<-"(S)-propranolol" # this is mapped because we know drugs from "drug"file and "dose-resp" file must be the same

# Cleaning the drug names for mapping
drugs<-drug_name_correction(table=drugs , column = "Compound.name") 
drug_with_ids_gbm<-drug_name_correction(table = drug_with_ids_gbm, column = "GBM.drugid")
# setdiff(drug_with_ids_gbm$clean.ids , drugs$clean.ids)
# setdiff(drugs$clean.ids , drug_with_ids_gbm$clean.ids)

# Drugs_with_ids$GBM.drugid includes drugs from both screen2 and screen3, there is no info available for some of these drugs in the drugs df
# These drugs are added to the drug object. 
drugs <- merge(drug_with_ids_gbm[!is.na(drug_with_ids_gbm$GBM.drugid), c("GBM.drugid", "clean.ids", "unique.drugid")],drugs, by="clean.ids",all.x = TRUE)
rownames(drugs) <- drugs$unique.drugid
drugs <- drugs[, 2:ncol(drugs)]
print("Drug object: done")
# ======================== Sensitivity data ======================== 
#Published AUC info
drug_cell<-read.delim(paste(input_dir , "HGCC_drug_response_AUC.txt", sep=""), header=T, sep="\t", stringsAsFactors = FALSE)
colnames(drug_cell)[colnames(drug_cell)=="X"]<-"GBM.drugid"
drug_cell<-drug_name_correction(table= drug_cell, column="GBM.drugid")
temp<-drug_name_correction(table= drug_with_ids_gbm, column="GBM.drugid")
temp<-merge(temp[,c("unique.drugid","clean.ids")], drug_cell, by="clean.ids", all.y=T)
drug_cell<-data.frame(temp[,-c(1,3)])
rm(temp)


#Transposing from wide to long 
drug_cell_long<- melt(setDT(drug_cell), id.vars = "unique.drugid", variable.name = "Pat")
drug_cell_long$EXP_details<-paste(drug_cell_long$unique.drugid, drug_cell_long$Pat , sep ="_")
print("Sensitivity data: done")

# =============Screen2 =============
print("Creating scr2_objects")
screen2<-read.delim(paste(input_dir, "Screen2-drugData.txt", sep=""), stringsAsFactors = FALSE)#Drug_dose_scr2 response data

#Doxirubicin is a typo of Doxorubicin
screen2$Drug[screen2$Drug == "Doxirubicin1"]<-"Doxorubicin1"
screen2$Drug[screen2$Drug == "Doxirubicin2"]<-"Doxorubicin2"
screen2$Drug[screen2$Drug == "Doxirubicin3"]<-"Doxorubicin3"
screen2$Drug[screen2$Drug == "Doxirubicin4"]<-"Doxorubicin4"
screen2$Drug[screen2$Drug == "Doxirubicin5"]<-"Doxorubicin5"
screen2$Drug[screen2$Drug == "Doxirubicin6"]<-"Doxorubicin6"
screen2$Drug[screen2$Drug == "Doxirubicin7"]<-"Doxorubicin7"
screen2$Drug[screen2$Drug == "Doxirubicin8"]<-"Doxorubicin8"

scr2_objects<-screen_objects(screen=screen2, cell_obj=all_cell_obj, drug_obj=drugs, drug_with_ids_gbm=drug_with_ids_gbm, d_c_l=drug_cell_long)

# ============= Screen3 =============
print("Creating scr3_objects")
screen3<-read.delim(paste(input_dir,"Screen3-drugData.txt", sep = ""), stringsAsFactors = FALSE)#Drug_dose_scr2 response data

#"Vinorelbinetartrate" is differently spelled in GBM_scr2 (with space) and GBM_scr3 (without space)
#To avoid defining a unique id for a same drug it is manually added to scr3_cur_drug through the below line.
screen3$Drug[screen3$Drug == "Vinorelbinetartrate"] <- "Vinorelbine tartrate"

scr3_objects<-screen_objects(screen=screen3, cell_obj=all_cell_obj, drug_obj=drugs, drug_with_ids=drug_with_ids_gbm, d_c_l=drug_cell_long)

#Removing the published AUC values for the cell-drug pairs mutual between screen2 and screen3
#The published AUCs have a higher correlation with calculated AUCs from screen2
#Since the original paper has not specified what screen do the published AUCs belong to
#It is more conservative to report the published AUCs for screen2 

Mutual_pairs<-merge(scr2_objects[["sen_profile"]],scr3_objects[["sen_profile"]], by="row.names")
cor(Mutual_pairs$aac_recomputed.x , Mutual_pairs$auc_published.x , use= "complete.obs")##Correlation between auc_published and aac_recomputed from screen2 (-0.9132723)
cor(Mutual_pairs$aac_recomputed.y , Mutual_pairs$auc_published.y , use= "complete.obs")##Correlation between auc_published and aac_recomputed from screen3 (-0.6995755)

scr3_objects[["sen_profile"]]$auc_published[rownames(scr3_objects[["sen_profile"]]) %in% Mutual_pairs$Row.names]<-NA

saveRDS(scr2_objects, paste0(out_dir, "scr2_objects.rds"))
saveRDS(scr3_objects, paste0(out_dir, "scr3_objects.rds"))