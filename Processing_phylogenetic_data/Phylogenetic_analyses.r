library(ape)
library(castor)
library(diagram)
library(HDInterval)
library(lubridate)
library(seraphim)
library(sf)
library(tidytree)
library(treeio)

# 1. Preparing the input for the preliminary discrete phylogeographic analysis (DTA) to identify distinct introduction events
# 2. Identifying the different clusters (clades following introduction events) from the DTA (+ preparing continuous analyses)
# 3. Extracting the spatio-temporal information embedded in posterior trees from the continuous phylogeofraphic inference
# 4. Generating the dispersal history graphs (mapped MCC trees, 80% HPD polygons) for the continuous phylogeographic inference
# 5. Estimating and plotting dispersal statistics associated with lineages (as estimated from the continuous phylogeographic inference)
# 6. Extraction and comparison of the population count at the start/end nodes (based on the continuous phylogeographic inference)
# 7. Extraction and comparison of the Mexican state at the start/end nodes (based on the continuous phylogeographic inference)

writingFiles = FALSE; showingPlots = TRUE; writingFiles = TRUE; showingPlots = FALSE
municipalities = shapefile("All_Mexican_shapefiles/Municipalities_MEX.shp")
if (!file.exists("All_Mexican_shapefiles/AGEB_study_area.shp"))
	{
		agebs = shapefile("All_Mexican_shapefiles/MEX_AGEB_areas.shp") # "MG_2020_Integrado/conjunto_de_datos/00a.shp”
		agebs_study_area = subset(agebs, CVE_ENT%in%c("09","15")); st_write(st_as_sf(agebs_study_area), "All_Mexican_shapefiles/AGEB_study_area.shp")
	}
agebs_study_area = shapefile("All_Mexican_shapefiles/AGEB_study_area.shp")
mobility_index = read.csv("Mobility_index_AGEB.csv", head=T)
missings = unique(c(mobility_index[which(!mobility_index[,"origen"]%in%agebs_study_area@data[,"CVEGEO"]),"origen"],
					mobility_index[which(!mobility_index[,"destino"]%in%agebs_study_area@data[,"CVEGEO"]),"destino"]))
missings = missings[order(missings)] # write(missings, "Unmatched_AGEB_IDs.txt")
mobility_index[which(mobility_index[,"origen"]%in%missings),"origen"] = paste0("0",mobility_index[which(mobility_index[,"origen"]%in%missings),"origen"])
mobility_index[which(mobility_index[,"destino"]%in%missings),"destino"] = paste0("0",mobility_index[which(mobility_index[,"destino"]%in%missings),"destino"])

	# NOTE (05/06/24): the current background dataset (received on the 15/05/24) has two issues: (i) Mexican sequences from 
	# outside the study area were not included in it and (ii) there a series of "MX" sequences whose precise origin is unknown
	# (i.e. we do not know if they were sampled or not within the defined study area) - they have to be discarded. Conclusions:
	# (i) a new TreeTime tree, based on an extended/clean set of background sequences will be inferred (by Marina) and (ii), in
	# the meantime, some analyses can be performed on the 10 selected HIV-Trace clusters (isolation-by-distance/mobility
	# analyses and continuous phylogeographic analyses along fixed topologies)

analysis = "25092024" # main analysis

# 1. Preparing the input for the preliminary discrete phylogeographic analysis (DTA) to identify distinct introduction events

tree = read.nexus(paste0("TreeTime_",analysis,".tre"))
data = read.csv(paste0("TreeTime_",analysis,".csv"), head=T)
if (length(outliers) != 0) tree = drop.tip(tree, outliers)
seqIDs = gsub("'","",tree$tip.label); locations = rep(NA, length(seqIDs))
collectionDates = rep(NA, length(seqIDs))
for (i in 1:length(seqIDs))
	{
		locations[i] = data[which(data[,1]==seqIDs[i]),"STATE"]
		collectionDates[i] = data[which(data[,1]==seqIDs[i]),"DATE"]
	}
for (i in 1:length(collectionDates))
	{
		collectionDate = unlist(strsplit(collectionDates[i],"\\/"))
		if (nchar(collectionDate[3]) == 4)
			{
				collectionDates[i] = paste(collectionDate[3],collectionDate[2],collectionDate[1],sep="-")
			}
		if (nchar(collectionDate[1]) == 4)
			{
				collectionDates[i] = paste(collectionDate[1],collectionDate[2],collectionDate[3],sep="-")
			}
	}
locations[which(!locations%in%c("ST_1","ST_2"))] = "other"
locations[which(locations%in%c("ST_1","ST_2"))] = "Mexico_state_and_city"
tab = cbind(seqIDs, locations, collectionDates); txt = c()
colnames(tab) = c("trait","location","collection_date")
for (i in 1:dim(tab)[1])
	{
		txt = c(txt, paste0(">",tab[i,"trait"]),"NNNN")
	}
if (writingFiles) write.table(tab, paste0("Preliminary_discrete_run/TreeTime_",analysis,".txt"), row.names=F, quote=F, sep="\t")
if (writingFiles) write(txt, paste0("Preliminary_discrete_run/TreeTime_",analysis,".fas"))
if (writingFiles) write.tree(tree, paste0("Preliminary_discrete_run/TreeTime_",analysis,".tre"))
tab = read.table(paste0("Preliminary_discrete_run/TreeTime_",analysis,".txt"), head=T)
tre = scan(paste0("Preliminary_discrete_run/TreeTime_",analysis,".tre"), what="", sep="\n", quiet=T)
xml = scan("Template_DTA_runs.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
xml = gsub("TEMPLATE",paste0("TreeTime_",analysis),xml)
sink(file=paste0("Preliminary_discrete_run/TreeTime_",analysis,".xml"))
for (i in 1:length(xml))
	{
		cat(xml[i],"\n")
		if (grepl("<taxa id=\"taxa\">",xml[i]))
			{
				for (j in 1:dim(tab)[1])
					{
						collectionDate = tab[j,"collection_date"]
						if (length(unlist(strsplit(tab[j,"collection_date"],"-"))) == 3)
							{
								collectionDate = decimal_date(ymd(collectionDate))
							}
						if (length(unlist(strsplit(tab[j,"collection_date"],"-"))) == 2)
							{
								collectionDate = decimal_date(ymd(paste0(collectionDate,"-15")))
							}
						cat("\t\t<taxon id=\"",tab[j,"trait"],"\">\n",sep="")
						cat("\t\t\t<date value=\"",collectionDate,"\" direction=\"forwards\" units=\"years\"/>\n",sep="")
						cat("\t\t\t<attr name=\"location\">\n",sep="")
						cat("\t\t\t\t",tab[j,"location"],"\n",sep="")
						cat("\t\t\t</attr>\n",sep="")
						cat("\t\t</taxon>\n",sep="")
					}
			}
		if (grepl("<alignment id=\"alignment\" dataType=\"nucleotide\">",xml[i]))
			{
				for (j in 1:dim(tab)[1])
					{
						cat("\t\t<sequence>\n",sep="")
						cat("\t\t\t<taxon idref=\"",tab[j,"trait"],"\"/>\n",sep="")
						cat("\t\t\t\tNNNN\n",sep="")
						cat("\t\t</sequence>\n",sep="")
					}
			}
		if (grepl("<newick id=\"startingTree\">",xml[i]))
			{
				cat("\t\t",tre,"\n",sep="")
			}
	}
sink(NULL)

# 2. Identifying the different clusters (clades following introduction events) from the DTA (+ preparing continuous analyses)

if (!file.exists(paste0("Preliminary_discrete_run/TreeTime_",analysis,".tree")))
	{
		file2 = paste0("Preliminary_discrete_run/TreeTime_",analysis,".tree"); file1 = gsub(".tree",".trees",file2)
		system(paste0("BEAST_version_1_10_4/bin/treeannotator -burninTrees 31 -heights keep ",file1," ",file2), ignore.stdout=F, ignore.stderr=F)
	}
tree = readAnnotatedNexus(paste0("Preliminary_discrete_run/TreeTime_",analysis,".tree"))
if (showingPlots)
	{
		metadata = read.csv(paste0("TreeTime_",analysis,".csv"), head=T, sep=",")
		if (nchar(unlist(strsplit(metadata[1,"DATE"],"-"))[1]) == 2)
			{
				mostRecentSamplingDate = max(decimal_date(dmy(gsub("\\/","-",metadata[,"DATE"]))))
			}
		if (nchar(unlist(strsplit(metadata[1,"DATE"],"-"))[1]) == 4)
			{
				mostRecentSamplingDate = max(decimal_date(ymd(gsub("\\/","-",metadata[,"DATE"]))))
			}
		if (nchar(unlist(strsplit(metadata[1,"DATE"],"\\/"))[1]) == 4)
			{
				mostRecentSamplingDate = max(decimal_date(ymd(gsub("\\/","-",metadata[,"DATE"]))))
			}
		rootHeight = max(nodeHeights(tree))
		root_time = mostRecentSamplingDate-rootHeight
		selectedLabels = seq(1980, 2025, 5)
		selectedDates = selectedLabels
		cols = rep("gray50",dim(tree$edge)[1])
		lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Mexico_state_and_city")&(tree$annotations[[i]]$location=="Mexico_state_and_city"))
							{
								cols[i] = rgb(222,67,39,220,maxColorValue=255); lwds[i] = 0.2
							}
					}
			}
		pdf(paste0("TreeTime_",analysis,"_NEW.pdf"), width=7, height=7); par(oma=c(1,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						if (tree$annotations[[i]]$location == "Mexico_state_and_city")
							{
								index = which(tree$edge[,2]==tree$edge[i,1])
								if (tree$annotations[[index]]$location != "Mexico_state_and_city")
									{
										nodelabels(node=tree$edge[i,2], pch=16, cex=0.5, col=rgb(222,67,39,220,maxColorValue=255))
										nodelabels(node=tree$edge[i,2], pch=1, cex=0.5, col="gray30", lwd=0.3)
									}
							}
					}
			}
		axis(lwd=0.2, at=selectedDates-root_time, labels=selectedLabels, cex.axis=0.48, mgp=c(0,-0.15,-0.6), lwd.tick=0.2, 
			 col.lab="gray30", col.axis="gray30", col.ticks="gray30", col="gray30", tck=-0.005, side=1)
		dev.off()
	}
introductions = c(); clusters1 = list(); clusters2 = list(); clusters2_txt = c()
for (i in 1:dim(tree$edge)[1])
	{
		if (is.null(tree$annotations[[i]]))
			{
				print(i)
			}	else	{
				if (tree$annotations[[i]]$location == "Mexico_state_and_city")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((length(index) > 0)&&(tree$annotations[[index]]$location != "Mexico_state_and_city"))
							{
								introductions = c(introductions, i)
							}
					}
			}
	}
for (i in 1:length(introductions))
	{
		if (i == 1) clusters1 = list()
		if (tree$edge[introductions[i],2]%in%tree$edge[,1])
			{
				subtree = tidytree::tree_subset(tree, tree$edge[introductions[i],2], levels_back=0)
				clusters1[[i]] = gsub("'","",subtree$tip.label)
			}	else	{
				clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[introductions[i],2]])
			}
	}
for (i in 2:length(clusters1))
	{
		for (j in 1:(i-1))
			{
				if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]]))
					{
						clusters1[[j]] = clusters1[[j]][which(!clusters1[[j]]%in%clusters1[[i]])]
					}
				if (sum(clusters1[[j]]%in%clusters1[[i]]) == length(clusters1[[j]]))
					{
						clusters1[[i]] = clusters1[[i]][which(!clusters1[[i]]%in%clusters1[[j]])]
					}
			}
	}
if (!file.exists(paste0("BEAST_RRW_analyses/All_the_metadata.csv")))
	{
		metadata1 = read.csv(paste0("TreeTime_25092024.csv"), head=T, sep=",")
		metadata2 = matrix(nrow=dim(metadata1)[1], ncol=5); metadata2[,1] = metadata1[,"ID"]
		colnames(metadata2) = c("sequence_ID","collection_date","AGEB","longitude","latitude")
		metadata2[,"collection_date"] = metadata1[,"DATE"]; metadata2[,"AGEB"] = metadata1[,"CP_AGEB"]
		for (i in 1:dim(metadata2)[1])
			{
				ageb_code = as.character(metadata2[i,"AGEB"]); maxArea = 0; pol = NULL
				if (substring(ageb_code,1,1) == "9") ageb_code = paste0("0",ageb_code)
				index = which(agebs_study_area@data[,"CVEGEO"]==ageb_code)
				if (length(index) != 1)
					{
						cat("Unmatched AGEB code: ",metadata2[i,"AGEB"],"\n",sep="")
					}	else	{
						metadata2[i,"AGEB"] = ageb_code
						for (j in 1:length(agebs_study_area@polygons[[index]]@Polygons))
							{
								if (maxArea < agebs_study_area@polygons[[index]]@Polygons[[j]]@area)
									{
										maxArea = agebs_study_area@polygons[[index]]@Polygons[[j]]@area
										pol = agebs_study_area@polygons[[index]]@Polygons[[j]]
									}
							}
						coords = spsample(pol, 1, type="random")@coords
						metadata2[i,"longitude"] = coords[,"x"]
						metadata2[i,"latitude"] = coords[,"y"]
					}
			}
		write.csv(metadata2, paste0("BEAST_RRW_analyses/All_the_metadata.csv"), quote=F, row.names=F)
	}	
metadata2 = read.csv(paste0("BEAST_RRW_analyses/All_the_metadata.csv"), head=T)
for (i in 1:length(introductions))
	{
		tab = c()
		if (i == 1)
			{
				clusters2 = list(); centroids = list()
			}
		for (j in 1:length(clusters1[[i]]))
			{
				index = which(metadata2[,"sequence_ID"]==unlist(strsplit(clusters1[[i]][j],"\\|"))[1])
				if (length(index) == 1)
					{
						collectionDate = unlist(strsplit(metadata2[index,"collection_date"],"\\/"))
						collectionDate = paste(collectionDate[3],collectionDate[2],collectionDate[1],sep="-")
						collection_date = decimal_date(ymd(collectionDate))
						municipality = metadata2[index,"AGEB"]
						latitude = metadata2[index,"latitude"]
						longitude = metadata2[index,"longitude"]
						line = cbind(collection_date, municipality, latitude, longitude)
						row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
					}
			}
		if (!is.null(dim(tab))) colnames(tab) = c("collection_date","AGEB","latitude","longitude")
		clusters2[[i]] = tab
	}
saveRDS(clusters1, paste0("Clusters1_",analysis,".rds"))
saveRDS(clusters2, paste0("Clusters2_",analysis,".rds"))
txt = c(); sizes = rep(NA, length(clusters2))
for (i in 1:length(clusters2))
	{
		if (!is.null(dim(clusters2[[i]]))) sizes[i] = dim(clusters2[[i]])[1]
	}
# clusters2_sorted = clusters2[order(sizes, decreasing=T)] # not used anymore to keep the same order and IDs for the clades
for (i in 1:length(clusters2)) txt = c(txt, paste0("\tCluster ",i," (n = ",dim(clusters2[[i]])[1],"):"), "", row.names(clusters2[[i]]), "", "")
write(txt, paste0("Clusters2_",analysis,".txt"))

template = scan("Template_RRW_runs.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
sink(file=paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/All_clades.xml"))
for (i in 1:length(template))
	{
		cat(template[i],"\n")
		if (grepl("Insert taxa blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
								for (k in 1:dim(clusters2[[j]])[1])
									{
										if (!is.na(clusters2[[j]][k,"longitude"]))
											{
												collectionDate = clusters2[[j]][k,"collection_date"]
												latitude = as.numeric(clusters2[[j]][k,"latitude"])
												longitude = as.numeric(clusters2[[j]][k,"longitude"])
												cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
												cat(paste0("\t\t\t<date value=\"",collectionDate,"\" direction=\"forwards\" units=\"years\"/>","\n"))
												cat("\t\t\t<attr name=\"latitude\">\n")
												cat(paste0("\t\t\t\t",latitude,"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t\t<attr name=\"longitude\">\n")
												cat(paste0("\t\t\t\t",longitude,"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t\t<attr name=\"coordinates\">\n")
												cat(paste0("\t\t\t\t",latitude," ",longitude,"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t</taxon>\n")
											}
									}
								cat("\t</taxa>","\n")
							}
					}
			}
		if (grepl("Insert alignment blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
								for (k in 1:dim(clusters2[[j]])[1])
									{
										if (!is.na(clusters2[[j]][k,"longitude"]))
											{
												cat("\t\t<sequence>\n")
												cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
												cat("\t\t\tNNNN\n")
												cat("\t\t</sequence>\n")
											}
									}
								cat("\t</alignment>","\n")
							}
					}
			}
		if (grepl("Insert pattern blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
								cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
								cat("\t</patterns>","\n")
							}
					}
			}
		if (grepl("Insert starting tree blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								tre = tidytree::tree_subset(tree, tree$edge[introductions[j],2], levels_back=0)
								tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"longitude"]))]
								tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
								if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
								write.tree(tre, paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/Clade_",j,".tre"))
								tre = scan(paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/Clade_",j,".tre"), what="", sep="\n", quiet=T)
								txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
								write(txt, paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/Clade_",j,".tre"))
								cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Clade_",j,".tre\">","\n"))
								cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
								cat("\t</empiricalTreeDistributionModel>","\n")
							}
					}
			}
		if (grepl("Insert tree model blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
								cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
								cat("\t\t<rootHeight>","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
								cat("\t\t</rootHeight>","\n")
								cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
								cat("\t\t</nodeHeights>","\n")
								cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
								cat("\t\t</nodeHeights>","\n")
								cat("\t</treeModel>","\n")
							}
					}
			}
		if (grepl("Insert arbitraryBranchRates blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<arbitraryBranchRates id=\"coordinates.diffusion.branchRates",j,"\">","\n"))
								cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
								cat("\t\t<rates>","\n")
								cat(paste0("\t\t\t<parameter id=\"coordinates.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
								cat("\t\t</rates>","\n")
								cat("\t</arbitraryBranchRates>","\n")
							}
					}
			}
		if (grepl("Insert distributionLikelihood blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<distributionLikelihood id=\"coordinates.diffusion.prior",j,"\">","\n"))
								cat("\t\t<data>","\n")
								cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
								cat("\t\t</data>","\n")
								cat("\t\t<distribution>","\n")
								cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
								cat("\t\t\t\t<shape>","\n")
								cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
								cat("\t\t\t\t</shape>","\n")
								cat("\t\t\t</onePGammaDistributionModel>","\n")
								cat("\t\t</distribution>","\n")
								cat("\t</distributionLikelihood>","\n")
							}
					}
			}
		if (grepl("Insert coordinates.traitLikelihood blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<multivariateTraitLikelihood id=\"coordinates.traitLikelihood",j,"\" traitName=\"coordinates\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
								cat("\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
								cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>"))
								cat("\t\t<traitParameter>","\n")
								cat(paste0("\t\t\t<parameter id=\"leaf.coordinates",j,"\"/>","\n"))
								cat("\t\t</traitParameter>","\n")
								cat("\t\t<conjugateRootPrior>","\n")
								cat("\t\t\t<meanParameter>","\n")
								cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
								cat("\t\t\t</meanParameter>","\n")
								cat("\t\t\t<priorSampleSize>","\n")
								cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
								cat("\t\t\t</priorSampleSize>","\n")
								cat("\t\t</conjugateRootPrior>","\n")
								cat(paste0("\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
								cat("\t</multivariateTraitLikelihood>","\n")
							}
					}
			}
		if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<continuousDiffusionStatistic id=\"coordinates.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
								cat(paste0("\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t</continuousDiffusionStatistic>","\n")
							}
					}
			}
		if (grepl("Insert scaleOperator blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
								cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
								cat("\t\t</scaleOperator>","\n")
							}
					}
			}
		if (grepl("Insert precisionGibbsOperator blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
								cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t\t\t<multivariateWishartPrior idref=\"coordinates.precisionPrior\"/>","\n")
								cat("\t\t</precisionGibbsOperator>","\n")
							}
					}
			}
		if (grepl("Insert distributionLikelihood blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<distributionLikelihood idref=\"coordinates.diffusion.prior",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"coordinates.diffusionRate",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("<!-- Insert logTree blocks -->",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((!is.null(dim(clusters2[[j]])))&&(dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
								cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
								cat("\t\t\t<joint idref=\"joint\"/>","\n")
								cat("\t\t\t<trait name=\"coordinates\" tag=\"coordinates\">","\n")
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t\t\t</trait>","\n")
								cat("\t\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
								cat("\t\t\t<trait name=\"rate\" tag=\"coordinates.rate\">","\n")
								cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
								cat("\t\t\t</trait>","\n")
								cat("\t\t</logTree>","\n")
							}
					}
			}
	}
sink(NULL)

# 3. Extracting the spatio-temporal information embedded in posterior trees from the continuous phylogeofraphic inference

if (analysis == "25092024") nberOfTreesToSample = 500
clusters2 = readRDS(paste0("Clusters2_",analysis,".rds"))
files = list.files(paste0("BEAST_RRW_analyses/TreeTime_",analysis)); files = files[which(grepl("\\.trees",files))]
for (i in 1:length(files))
	{
		txt1 = scan(paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/",files[i]), what="", sep="\n", quiet=T, blank.lines.skip=F)
		if (txt1[length(txt1)] != "End;")
			{
				txt2 = txt1[which(!grepl("tree STATE_",txt1))]
				indices = which(grepl("tree STATE_",txt1))
				if (analysis == "25092024") txt2 = c(txt2, txt1[indices], "End;")
				write(txt2, paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/",files[i]))
			}
	}
files = list.files(paste0("BEAST_RRW_analyses/TreeTime_",analysis)); files = files[which(grepl("\\.trees",files))]
for (i in 1:length(files))
	{
		txt1 = scan(paste0("BEAST_RRW_analyses/TreeTime_",analysis,"/",files[i]), what="", sep="\n", quiet=T, blank.lines.skip=F)
		txt2 = txt1[which(!grepl("tree STATE_",txt1))]; txt2 = txt2[1:(length(txt2)-1)]; indices = which(grepl("tree STATE_",txt1))
		if (analysis == "25092024") indices = indices[(round(0.35*length(indices))+1):length(indices)]
		if (length(indices) > nberOfTreesToSample)
			{
				indices = indices[sample(1:length(indices), nberOfTreesToSample, replace=F)]
			}
		txt2 = c(txt2, txt1[indices], "End;")
		write(txt2, paste0("All_tree_extraction_files/TreeTime_",analysis,"/",files[i]))
	}
files = list.files(paste0("BEAST_RRW_analyses/TreeTime_",analysis)); files = files[which(grepl("\\.trees",files))]
for (i in 1:length(files))
	{
		file1 = files[i]; file2 = gsub("\\.trees","\\.tree",file1)
		system(paste0("BEAST_version_1_10_4/bin/treeannotator -burninTrees 0 -heights keep All_tree_extraction_files/TreeTime_",analysis,"/",file1," All_tree_extraction_files/TreeTime_",analysis,"/",file2), ignore.stdout=F, ignore.stderr=F)
	}
source("MCC_tree_extractions.r")
for (i in 1:length(files))
	{
		if (!file.exists(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","\\.csv",files[i]))))
			{
				index = as.numeric(gsub("\\.trees","",unlist(strsplit(files[i],"_"))[2]))
				mostRecentSamplingDatum = max(as.numeric(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collection_date"]))
				mcc_tre = readAnnotatedNexus(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","\\.tree",files[i]))); dates = c()
				mcc_tab = MCC_tree_extractions(mcc_tre, mostRecentSamplingDatum)
				write.csv(mcc_tab, paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","\\.csv",files[i])), row.names=F, quote=F)
			}
	}
randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 10
for (i in 1:length(files))
	{
		localTreesDirectory = paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","_ext",files[i]))
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				index = as.numeric(gsub("\\.trees","",unlist(strsplit(files[i],"_"))[2])); burnIn = 0
				mostRecentSamplingDatum = max(as.numeric(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collection_date"]))
				allTrees = scan(file=paste0("All_tree_extraction_files/TreeTime_",analysis,"/",files[i]), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
	}
admin0 = shapefile("Mex_GADM_shapefiles/Mexico_admin0.shp")
for (i in 1:length(files))
	{
		tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","\\.csv",files[i])), head=T)
		if (abs(tab[1,"startLon"]) > 1000) # untransformed coordinates
			{
				spatial_points_sf1 = st_as_sf(tab[,c("startLon","startLat")], coords=c("startLon","startLat"), crs=crs(agebs_study_area))
				spatial_points_sf2 = st_transform(spatial_points_sf1, crs(admin0))
				pts2 = as(spatial_points_sf2, "Spatial")@coords; tab[,c("startLon","startLat")] = pts2[,1:2]
				spatial_points_sf1 = st_as_sf(tab[,c("endLon","endLat")], coords=c("endLon","endLat"), crs=crs(agebs_study_area))
				spatial_points_sf2 = st_transform(spatial_points_sf1, crs(admin0))
				pts2 = as(spatial_points_sf2, "Spatial")@coords; tab[,c("endLon","endLat")] = pts2[,1:2]
				write.csv(tab, paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","\\.csv",files[i])), row.names=F, quote=F)
			}
	}
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:length(files))
	{
		for (j in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","_ext",files[i]),"/TreeExtractions_",j,".csv"), head=T)
				if (abs(tab[1,"startLon"]) > 1000) # untransformed coordinates
					{
						spatial_points_sf1 = st_as_sf(tab[,c("startLon","startLat")], coords=c("startLon","startLat"), crs=crs(agebs_study_area))
						spatial_points_sf2 = st_transform(spatial_points_sf1, crs(admin0))
						pts2 = as(spatial_points_sf2, "Spatial")@coords; tab[,c("startLon","startLat")] = pts2[,1:2]
						spatial_points_sf1 = st_as_sf(tab[,c("endLon","endLat")], coords=c("endLon","endLat"), crs=crs(agebs_study_area))
						spatial_points_sf2 = st_transform(spatial_points_sf1, crs(admin0))
						pts2 = as(spatial_points_sf2, "Spatial")@coords; tab[,c("endLon","endLat")] = pts2[,1:2]
						write.csv(tab, paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","_ext",files[i]),"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
					}
			}
	}
if (!file.exists(paste0("All_tree_extraction_files/MCCtree_",analysis,".csv")))
	{
		file.copy(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades.csv"), paste0("All_tree_extraction_files/MCCtree_",analysis,".csv"))
	}
if (file.exists(paste0("All_tree_extraction_files/MCCtree_",analysis,".csv")))
	{
		tab = read.csv(paste0("All_tree_extraction_files/MCCtree_",analysis,".csv"), head=T); tmp = tab
		if (!"startAGEB"%in%colnames(tab))
			{
				spatial_points_sf1 = st_as_sf(tab[,c("startLon","startLat")], coords=c("startLon","startLat"), crs=crs(admin0))
				spatial_points_sf2 = st_transform(spatial_points_sf1, crs(agebs_study_area))
				pts2 = as(spatial_points_sf2, "Spatial")@coords; tmp[,c("startLon","startLat")] = pts2[,1:2]
				spatial_points_sf1 = st_as_sf(tab[,c("endLon","endLat")], coords=c("endLon","endLat"), crs=crs(admin0))
				spatial_points_sf2 = st_transform(spatial_points_sf1, crs(agebs_study_area))
				pts2 = as(spatial_points_sf2, "Spatial")@coords; tmp[,c("endLon","endLat")] = pts2[,1:2]
				tab$startAGEB = as.character(sp::over(SpatialPoints(tmp[,c("startLon","startLat")], crs(agebs_study_area)), agebs_study_area)[,"CVEGEO"])
				tab$endAGEB = as.character(sp::over(SpatialPoints(tmp[,c("endLon","endLat")], crs(agebs_study_area)), agebs_study_area)[,"CVEGEO"])
				write.csv(tab, paste0("All_tree_extraction_files/MCCtree_",analysis,".csv"), row.names=F, quote=F)
			}
	}
for (i in 1:length(files))
	{
		for (j in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","_ext",files[i]),"/TreeExtractions_",j,".csv"), head=T); tmp = tab
				if (!"startAGEB"%in%colnames(tab))
					{
						spatial_points_sf1 = st_as_sf(tab[,c("startLon","startLat")], coords=c("startLon","startLat"), crs=crs(admin0))
						spatial_points_sf2 = st_transform(spatial_points_sf1, crs(agebs_study_area))
						pts2 = as(spatial_points_sf2, "Spatial")@coords; tmp[,c("startLon","startLat")] = pts2[,1:2]
						spatial_points_sf1 = st_as_sf(tab[,c("endLon","endLat")], coords=c("endLon","endLat"), crs=crs(admin0))
						spatial_points_sf2 = st_transform(spatial_points_sf1, crs(agebs_study_area))
						pts2 = as(spatial_points_sf2, "Spatial")@coords; tmp[,c("endLon","endLat")] = pts2[,1:2]
						tab$startAGEB = as.character(sp::over(SpatialPoints(tmp[,c("startLon","startLat")], crs(agebs_study_area)), agebs_study_area)[,"CVEGEO"])
						tab$endAGEB = as.character(sp::over(SpatialPoints(tmp[,c("endLon","endLat")], crs(agebs_study_area)), agebs_study_area)[,"CVEGEO"])
						write.csv(tab, paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","_ext",files[i]),"/TreeExtractions_",j,".csv"), row.names=F, quote=F)
					}			
			}
	}
for (i in 1:length(files))
	{
		tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","\\.csv",files[i])), head=T)
		cladeIDs = matrix(nrow=dim(tab)[1], ncol=1); colnames(cladeIDs) = "cladeID"
		cladeIDs[,"cladeID"] = gsub("\\.trees","",unlist(strsplit(files[i],"_"))[2])
		if (i == 1)
			{
				all = cbind(tab, cladeIDs)
			}	else	{
				maxNodeID = max(all[,c("node1","node2")])
				tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
				all = rbind(all, cbind(tab,cladeIDs))
			}
	}
write.csv(all, paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades.csv"), row.names=F, quote=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfExtractionFiles)
	{
		for (j in 1:length(files))
			{
				tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",gsub("\\.trees","_ext",files[j]),"/TreeExtractions_",i,".csv"), head=T)
				cladeIDs = matrix(nrow=dim(tab)[1], ncol=1); colnames(cladeIDs) = "cladeID"
				cladeIDs[,"cladeID"] = gsub("\\.trees","",unlist(strsplit(files[i],"_"))[2])
				if (j == 1)
					{
						all = cbind(tab, cladeIDs)
					}	else	{
						maxNodeID = max(all[,c("node1","node2")])
						tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
						all = rbind(all, cbind(tab,cladeIDs))
					}
			}
		write.csv(all, paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

# 4. Generating the dispersal history graphs (mapped MCC trees, 80% HPD polygons) for the continuous phylogeographic inference

		 # (in Figures 3, panels b-d are based on the following selected clades: 363, 1583, and 1817)

admin0 = shapefile("Mex_GADM_shapefiles/Mexico_admin0.shp")
admin1 = shapefile("Mex_GADM_shapefiles/Mexico_admin1.shp")
Mexico_state_and_city = subset(admin1, NAME_1%in%c("México","Distrito Federal"))
population_count = crop(raster("MEX_pop_2020_100.tif"), Mexico_state_and_city)
population_count = raster::mask(population_count, Mexico_state_and_city)
population_count_log = population_count; population_count_log[] = log10(population_count_log[]+1)
pop_log_cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[6:81])
colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[16:116])
colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[16:116])
clusters2 = readRDS(paste0("Clusters2_",analysis,".rds"))
for (i in 1:length(clusters2))
	{
		if (!is.null(clusters2[[i]]))
			{
				tab = cbind(as.numeric(clusters2[[i]][,c("longitude")]),as.numeric(clusters2[[i]][,c("latitude")]))
				colnames(tab) = c("longitude","latitude"); tab = data.frame(tab)
				spatial_points_sf1 = st_as_sf(tab, coords=c("longitude","latitude"), crs=crs(agebs_study_area))
				spatial_points_sf2 = st_transform(spatial_points_sf1, crs(admin0))
				pts2 = as(spatial_points_sf2, "Spatial")@coords; clusters2[[i]][,c("longitude","latitude")] = pts2[,1:2]
			}
	}
files = list.files(paste0("All_tree_extraction_files/TreeTime_",analysis)); clades = gsub("_ext","",files[which(grepl("_ext",files))])
mcc_merged = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades.csv"), head=T)
minYear = min(mcc_merged[,"startYear"]); maxYear = max(mcc_merged[,"endYear"])
for (h in 0:1) # length(clades)) # for (h in c(50,178,268))
	{
		if (h < 2) cexNode = 0.3
		if (h > 1) cexNode = 0.8
		if (h > 0)
			{
				localTreesDirectory = paste0("All_tree_extraction_files/TreeTime_",analysis,"/",clades[h],"_ext"); nberOfExtractionFiles = 100
				percentage = 80; prob = percentage/100; precision = 1; croppingPolygons = FALSE
				mcc = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",clades[h],".csv"), head=T); startDatum = min(mcc[,"startYear"])
				if ((h > 1)&(length(which(!mcc[,"node2"]%in%mcc[,"node1"])) >= 100)) cexNode = 0.5
				polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
				startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
				startYears_colours = colourScale[startYears_indices]
				endYears_colours = colourScale[endYears_indices]
				polygons_colours = rep(NA, length(polygons))		
				for (i in 1:length(polygons))
					{
						date = as.numeric(names(polygons[[i]]))
						polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
						polygons_colours[i] = paste0(colourScale[polygon_index],"40")
					}
			}
		if (h == 0) pdf(paste0("All_tree_extraction_files/TreeTime_",analysis,"/1or2_clades.pdf"), width=6, height=6) # dev.new(width=6, height=6)
		if (h >= 1) pdf(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",clades[h],".pdf"), width=6, height=6) # dev.new(width=6, height=6)
		par(oma=c(0,0,0,0), mar=c(0.5,0.5,0.5,0.5), lwd=0.2, col="gray30")
		plot(Mexico_state_and_city, border=NA, col=NA, lwd=0.5)
		plot(population_count_log, col=pop_log_cols, add=T, legend=F)
		plot(Mexico_state_and_city, border="gray30", col=NA, lwd=0.5, add=T)
		if (h >= 1)
			{
				for (i in 1:length(polygons))
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								# polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(admin0)
						if (croppingPolygons == TRUE) pol = crop(pol, admin0)
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
				selectedBranches = 1:dim(mcc)[1]
				for (i in selectedBranches)
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
								    arr.width=0, lwd=0.2, lty=1, lcol="black", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
			}
		if (h <= 1)
			{
				for (i in 1:length(clusters2))
					{
						if (!is.null(clusters2[[i]]))
							{
								if (sum(!is.na(as.numeric(clusters2[[i]][,"longitude"]))) < 3)
									{
										if (sum(!is.na(as.numeric(clusters2[[i]][,"longitude"]))) == 2)
											{
												indices = which(!is.na(as.numeric(clusters2[[i]][,"longitude"])))
												if (length(indices) == 2)
													{
														curvedarrow(cbind(as.numeric(clusters2[[i]][indices[1],"longitude"]),as.numeric(clusters2[[i]][indices[1],"latitude"])),
																	cbind(as.numeric(clusters2[[i]][indices[2],"longitude"]),as.numeric(clusters2[[i]][indices[2],"latitude"])),
																	arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="black", arr.col=NA, arr.pos=F,
																	curve=0.1, dr=NA, endhead=F)
													}
											}
									}
							}
					}
				for (i in 1:length(clusters2))
					{
						if (!is.null(clusters2[[i]]))
							{
								if (sum(!is.na(as.numeric(clusters2[[i]][,"longitude"]))) < 3)
									{
										for (j in 1:dim(clusters2[[i]])[1])
											{
												if (!is.na(as.numeric(clusters2[[i]][j,"longitude"])))
													{
														plotTheNode = TRUE
														if (plotTheNode)
															{
																index = (((as.numeric(clusters2[[i]][j,"collection_date"])-minYear)/(maxYear-minYear))*100)+1
																points(as.numeric(clusters2[[i]][j,"longitude"]),as.numeric(clusters2[[i]][j,"latitude"]),pch=16,col=colourScale[index],cex=cexNode)
																points(as.numeric(clusters2[[i]][j,"longitude"]),as.numeric(clusters2[[i]][j,"latitude"]),pch=1,col="black",cex=cexNode,lwd=0.4)
															}
													}
											}
									}
							}
					}
			}
		if (h >= 1)
			{
				for (i in rev(selectedBranches))
					{
						if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
							{
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="black", cex=cexNode, lwd=0.2)
							}
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="black", cex=cexNode, lwd=0.2)
					}
			}
		selectedValues = c(1,2); selectedLabels = c(10,100)
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(population_count_log[], na.rm=T); rast[2] = max(population_count_log[], na.rm=T)
		plot(rast, legend.only=T, add=T, col=pop_log_cols, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.650,0.938,0.160,0.172),
			 legend.args=list(text="", cex=0.6, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.7, col.axis="gray30", line=0, mgp=c(0,-0.02,0), at=selectedValues, label=selectedLabels))
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minYear; rast[2] = maxYear
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.550,0.938,0.100,0.112),
			 legend.args=list(text="", cex=0.6, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.7, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))
		dev.off()
	}

# 5. Estimating and plotting dispersal statistics associated with lineages (as estimated from the continuous phylogeographic inference)

if (file.exists(paste0("All_tree_extraction_files/MCCtree_",analysis,".csv"))) # Analyses based on the MCC tree extractions:
	{
		tab = read.csv(paste0("All_tree_extraction_files/MCCtree_",analysis,".csv"), head=T)
		clade_IDs = unique(tab[,"cladeID"]); clade_IDs = clade_IDs[order(clade_IDs)]; indices = c()
		for (i in 1:length(clade_IDs))
			{
				cat(i,"/",length(clade_IDs),"\n",sep="")
				sub = tab[which(tab[,"cladeID"]==clade_IDs[i]),]
				for (j in 2:dim(sub)[1])
					{
						for (k in 1:(j-1))
							{
								if ((!sub[j,"node2"]%in%sub[,"node1"])&(!sub[k,"node2"]%in%sub[,"node1"]))
									{
										indices = rbind(indices, cbind(k,j))
									}
							}
					}
			}
		distTree = matrix(nrow=dim(indices)[1], ncol=1); n = 0; nn = 1000
		for (i in 1:dim(indices)[1])
			{
				n = n+1
				if (n == nn)
					{
						nn = nn+1000; cat(n,"/",dim(indices)[1],"\n",sep="")
					}
				index1 = indices[i,1]; index2 = indices[i,2]
				indices1 = index1; root = FALSE
				while (root == FALSE)
					{	
						if (tab[indices1[length(indices1)],"node1"]%in%tab[,"node2"])
							{
								indices1 = c(indices1, which(tab[,"node2"]==tab[indices1[length(indices1)],"node1"]))
							}	else	{
								root = TRUE
							}
					}
				indices2 = index2; root = FALSE
				while (root == FALSE)
					{	
						if (tab[indices2[length(indices2)],"node1"]%in%tab[,"node2"])
							{
								indices2 = c(indices2, which(tab[,"node2"]==tab[indices2[length(indices2)],"node1"]))
							}	else	{
								root = TRUE
							}
					}
				indices3 = indices1[which(indices1%in%indices2)]; patristic_dis = NULL
				if (length(indices3) == 0)
					{
						patristic_dis = sum(tab[c(indices1,indices2),"length"])
					}	else	{
						patristic_dis = sum(tab[c(indices1[which(!indices1%in%indices3)],indices2[which(!indices2%in%indices3)]),"length"])
					}
				distTree[i,1] = patristic_dis
			}
		distGeo = matrix(nrow=dim(indices)[1], ncol=1); n = 0; nn = 1000
		for (i in 1:dim(indices)[1])
			{
				n = n+1
				if (n == nn)
					{
						nn = nn+1000; cat(n,"/",dim(indices)[1],"\n",sep="")
					}
				index1 = indices[i,1]; index2 = indices[i,2]
				distGeo[i,1] = rdist.earth(cbind(tab[index1,"endLon"],tab[index1,"endLat"]),
										   cbind(tab[index2,"endLon"],tab[index2,"endLat"]), miles=F, R=NULL)
			}
		mobility = matrix(nrow=dim(indices)[1], ncol=2); n = 0; nn = 1000
		for (i in 1:dim(indices)[1])
			{
				n = n+1
				if (n == nn)
					{
						nn = nn+1000; cat(n,"/",dim(indices)[1],"\n",sep="")
					}
				startLocation = tab[indices[i,1],"endAGEB"]; endLocation = tab[indices[i,2],"endAGEB"]
				index = which((mobility_index[,"origen"]==startLocation)&(mobility_index[,"destino"]==endLocation))
				if (length(index) == 1) mobility[i,1] = mobility_index[index,"peso_mediana"]
				startLocation = tab[indices[i,2],"endAGEB"]; endLocation = tab[indices[i,1],"endAGEB"]
				index = which((mobility_index[,"origen"]==startLocation)&(mobility_index[,"destino"]==endLocation))
				if (length(index) == 1) mobility[i,2] = mobility_index[index,"peso_mediana"]				
			}
		distances = data.frame(indices); distances$distTree = distTree[,1]; distances$distGeo = distGeo[,1]
		distances$mobility1 = mobility[,1]; distances$mobility2 = mobility[,2]
		colnames(distances) = c("index1","index2","patristic","geographic","mobility1","mobility2")
		write.csv(distances, paste0("All_tree_extraction_files/Distances_",analysis,".csv"), row.names=F, quote=F)
		cor(distances$patristic, distances$geographic, method="spearman") # 0.269
		dS3 = c(distances$mobility1,distances$mobility2) # When only considering pair of tips with a mobility metric:
		dS1 = c(distances$patristic,distances$patristic); dS2 = c(distances$geographic,distances$geographic)
		dS1 = dS1[which(!is.na(dS3))]; dS2 = dS2[which(!is.na(dS3))]; dS3 = dS3[which(!is.na(dS3))]
		cor(dS1, dS2, method="spearman") # Spearman correlation between patristic and geographic distance: 0.186
		cor(dS1, dS3, method="spearman") # Spearman correlation between patristic and the mobility metric: -0.141
		cor(dS2, dS3, method="spearman") # Spearman correlation between geographic and the mobility metric: -0.510
	}
source("spreadStatistics_mod.r"); source("isolationByMobility_fct.r")
files = list.files(paste0("All_tree_extraction_files/TreeTime_",analysis)); directories = files[which(grepl("_ext",files))]
for (i in 2:length(directories))
	{
		localTreesDirectory = paste0("All_tree_extraction_files/TreeTime_",analysis,"/",directories[i])
		outputName = paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/",gsub("_ext","",directories[i]))
		nberOfExtractionFiles = 100; timeSlices = 100; onlyTipBranches = FALSE; showingPlots = FALSE; nberOfCores = 20; slidingWindow = 1
		spreadStatistics_mod(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)		
	}
IBM_signals = matrix(nrow=(length(directories)-1), ncol=2); colnames(IBM_signals) = c("IBM_signal_rP1","IBM_signal_rP2")
for (i in 2:length(directories))
	{
		localTreesDirectory = paste0("All_tree_extraction_files/TreeTime_",analysis,"/",directories[i]); nberOfExtractionFiles = 1
		outputName = paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/",gsub("_ext","",directories[i]))	
		IBM_signals[i-1,] = isolationByMobility_fct(localTreesDirectory, nberOfExtractionFiles, mobility_index, outputName)
	}
write.csv(IBM_signals, "IBMsignal_25092024.csv", row.names=F, quote=F)
IBM_signal_1_mean = round(mean(IBM_signals[,1], na.rm=T), 3); IBM_signal_1_ci = round(t.test(IBM_signals[,1])$conf.int, 3)
IBM_signal_2_mean = round(mean(IBM_signals[,2], na.rm=T), 3); IBM_signal_2_ci = round(t.test(IBM_signals[,2])$conf.int, 3)
cat("IBM signal, rP#1: mean = ",IBM_signal_1_mean,", 95% CI = [",IBM_signal_1_ci[1],", ",IBM_signal_1_ci[2],"]\n",sep="")
cat("IBM signal, rP#2: mean = ",IBM_signal_2_mean,", 95% CI = [",IBM_signal_2_ci[1],", ",IBM_signal_2_ci[2],"]\n",sep="")
	# IBM signal, rP#1: mean = -0.086, 95% CI = [-0.129, -0.043], and rP#2: mean = -0.097, 95% CI = [-0.141, -0.054]
clusters2 = readRDS(paste0("Clusters2_",analysis,".rds")); cluster_sizes = rep(NA, length(clusters2))
for (i in 1:length(clusters2)) 
	{
		if (!is.null(dim(clusters2[[i]])[1])) cluster_sizes[i] = dim(clusters2[[i]])[1]
	}
files = list.files(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/"))
files = files[which(grepl("_dispersal_statistics.txt",files))]
WDCs = list(); IBDs = list()
for (i in 1:length(files))
	{
		tab = read.table(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/",files[i]), head=T)
		WDCs[[i]] = tab[,"weighted_diffusion_coefficient"]; IBDs[[i]] = tab[,"isolation_by_distance_signal_rP2"]
	}
if (showingPlots) # exploratory plots
	{
		dev.new(width=9, height=3); par(mfrow=c(1,3), oma=c(1,1,1,1), mar=c(3,3,1,1), col="gray30", lwd=0.2)
		boxplot(log10(cluster_sizes), frame=F, lwd=0.5)
		plot(density(WDCs[[1]][!is.na(WDCs[[1]])]), lwd=0.1, col="red", main=NA, frame=F, axes=T, xlim=c(0,200), ylim=c(0,0.4))
		for (i in 2:length(WDCs)) lines(density(WDCs[[i]][!is.na(WDCs[[i]])]), col="red")
		plot(density(IBDs[[1]][!is.na(IBDs[[1]])]), lwd=0.1, col="red", main=NA, frame=F, axes=T, xlim=c(-1,1), ylim=c(0,90))
		for (i in 2:length(IBDs)) lines(density(IBDs[[i]][!is.na(IBDs[[i]])]), col="red")
	}
files = list.files(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/")); files = files[which(grepl("_dispersal_statistics.txt",files))]
sizes = rep(NA, length(files)); medians = rep(NA, length(files)); vSs1 = list(); all1 = c(); vSs1_log = list(); all1_log = c(); vSs2 = list(); all2 = c()
for (i in 1:length(files))
	{
		clade_ID = gsub("_estimated_dispersal_statistics.txt","",files[i])
		tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/",clade_ID,"_ext/TreeExtractions_1.csv"), head=T)
		sizes[i] = length(which(!tab[,"node2"]%in%tab[,"node1"]))
		tab = read.table(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/",files[i]), head=T)
		vS = tab[,"weighted_diffusion_coefficient"]; medians[i] = median(vS); vSs1[[i]] = vS; vSs1_log[[i]] = log(vS)
	}	# print(files[which(sizes > 100)])
sizes = sizes[order(medians,decreasing=F)]; vSs1 = vSs1[order(medians,decreasing=F)]; vSs1_log = vSs1_log[order(medians,decreasing=F)]
for (i in 1:length(files))
	{
		I = i
		if (i < 10) I = paste0("0",I)
		if (i < 100) I = paste0("0",I)
		if (i < 1000) I = paste0("0",I)
		if (i >= 1000) I = paste0("",I)
		mat = cbind(vSs1[[i]],rep(paste0(I,"_",files[i]),length(vSs1[[i]])))
		all1 = rbind(all1, mat)
		mat = cbind(log(vSs1[[i]]),rep(paste0(I,"_",files[i]),length(vSs1[[i]])))
		all1_log = rbind(all1_log, mat)
	}
all1 = as.data.frame(all1); all1$V1 = sapply(all1$V1, as.numeric)
all1_log = as.data.frame(all1_log); all1_log$V1 = sapply(all1_log$V1, as.numeric)
for (i in 1:length(files))
	{
		tab = read.table(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_stats/",files[i]), head=T)
		vSs2[[i]] = tab[,"isolation_by_distance_signal_rP2"]
	}
vSs2 = vSs2[order(medians,decreasing=F)]
for (i in 1:length(files))
	{
		I = i
		if (i < 10) I = paste0("0",I)
		if (i < 100) I = paste0("0",I)
		if (i < 1000) I = paste0("0",I)
		if (i >= 1000) I = paste0("",I)
		mat = cbind(vSs2[[i]],rep(paste0(I,"_",files[i]),length(vSs2[[i]])))
		all2 = rbind(all2, mat)
	}
all2 = as.data.frame(all2); all2$V1 = sapply(all2$V1, as.numeric)

pdf("Dispersal_statistics_NEW1.pdf", width=10.0, height=12.0) # dev.new(width=10.0, height=12.0)
par(mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(3.0,0.5,0.5,0.5), lwd=0.5, col="gray30"); par(fig=c(0,7.5,0,10)/10)
ridgeline_mod1 = function(x, y, sizes, bw="nrd0", mode=F, col="gray", border, lty=1, lwd=1, bty="o", labels=NULL, palette, axes=T) 
	{
		dens = tapply(x, y, density, bw=bw)
		xs = Map(getElement, dens, "x")
		ys = Map(getElement, dens, "y")
		ys = Map(function(x) (x-min(x))/max(x-min(x))*1.5, ys)
		ys = Map(`+`, ys, length(ys):1)
		op = par(no.readonly=T)
		par(mar=op$mar); plot.new()
		plot.window(xlim=c(-4.6,log(850)), ylim=c(15,length(ys)-12))
		# abline(h=length(ys):1, col=col, lwd=0.3)
		cols = hcl.colors(length(ys), "Fall", alpha=1)
		border = rep(1, length(ys))
		Map(polygon, xs, ys, col=cols, border=border, lty=lty, lwd=lwd)
		axis(side=1, lwd.tick=0.3, cex.axis=0.8, mgp=c(0,0.25,0), lwd=0, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=c(log(10^-10),log(1),log(10),log(100),log(1000),log(10000),log(100000),log(1000000)), labels=c("","1","10","100","1000","10000","100000","1000000"))
		# axis(side=2, lwd.tick=0.0, cex.axis=0.8, mgp=c(0,0.38,0), lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 # at=c(1:length(files))+0.07, labels=rev(sizes), las=1)
		yss = length(ys):1
		for (i in 1:length(sizes))
			{
				points(-4.5, yss[i], cex=log(sizes[i]), pch=16, col=cols[i])
				points(-4.5, yss[i], cex=log(sizes[i]), pch=1, col="gray30", lwd=0.3)
			}
		points(-3.0, yss[i], cex=log(500), pch=16, col="#E5E5E5")
		points(-3.0, yss[i], cex=log(500), pch=1, col="gray30", lwd=0.3)
		points(-3.0, yss[i], cex=log(100), pch=16, col="#E5E5E5")
		points(-3.0, yss[i], cex=log(100), pch=1, col="gray30", lwd=0.3)
		points(-3.0, yss[i], cex=log(10), pch=16, col="#E5E5E5")
		points(-3.0, yss[i], cex=log(10), pch=1, col="gray30", lwd=0.3)
		title(xlab=expression("Weighted diffusion coefficient (km"^2*"/year)"), cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
		box(bty=bty, lwd=0.3); par(op)
	}
x = all1_log[,1]; y = all1_log[,2]; ridgeline_mod1(x, y, sizes, axes=T, lwd=0.2)
par(fig=c(7.5,10,0,10)/10); par(new=T, mar=c(3.0,0,0.5,0.4), lwd=0.3, col="gray30")
ridgeline_mod2 = function(x, y, bw="nrd0", mode=F, col="gray", border, lty=1, lwd=1, bty="o", labels=NULL, palette, axes=T) 
	{
		dens = tapply(x, y, density, bw=bw)
		xs = Map(getElement, dens, "x")
		ys = Map(getElement, dens, "y")
		ys = Map(function(x) (x-min(x))/max(x-min(x))*1.5, ys)
		ys = Map(`+`, ys, length(ys):1)
		op = par(no.readonly=T)
		par(mar=op$mar); plot.new()
		plot.window(xlim=c(-0.98,0.98), ylim=c(15,length(ys)-12))
		# abline(h=length(ys):1, col=col, lwd=0.3)
		abline(v=0, col=col, lwd=0.6, lty="aa")
		cols = hcl.colors(length(ys), "Fall", alpha=1)
		border = rep(1, length(ys))
		Map(polygon, xs, ys, col=cols, border=border, lty=lty, lwd=lwd)
		axis(side=1, lwd.tick=0.3, cex.axis=0.8, mgp=c(0,0.25,0), lwd=0, tck=-0.030, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-1,1,0.5), labels=seq(-1,1,0.5))
		title(xlab=expression("IBD signal (r"[P]*")"), cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
		box(bty=bty, lwd=0.3); par(op)
	}
x = all2[,1]; y = all2[,2]; ridgeline_mod2(x, y, axes=T, lwd=0.2)
dev.off()

pdf("Dispersal_statistics_NEW2.pdf", width=2, height=3) # dev.new(width=2, height=3)
par(mgp=c(0,0.35,0), oma=c(0,0,0,0), mar=c(1,0,1,0), lwd=0.75, col="gray30")
boxplot(log10(cluster_sizes), ylim=c(0,log10(600)), axes=F, col="gray90")
axis(2, pos=0.70, tck=-0.035, lwd=0.3, cex.axis=0.8, col.axis="gray30", at=c(-1,0,1,2,3), labels=c(0.1,1,10,100,1000))
dev.off()

# 6. Extraction and comparison of the population count at the start/end nodes (based on the continuous phylogeographic inference)

population_count = raster("MEX_pop_2020_100.tif")
for (i in 1:100)
	{
		tab1 = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/TreeExtractions_",i,".csv"), head=T)
		tab2 = matrix(nrow=dim(tab1)[1], ncol=2); colnames(tab2) = c("startPopCount", "endPopCount")
		tab2[,"startPopCount"] = raster::extract(population_count, tab1[,c("startLon","startLat")])
		tab2[,"endPopCount"] = raster::extract(population_count, tab1[,c("endLon","endLat")])
		write.csv(tab2, paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/PopCountValues_",i,".csv"), row.names=F, quote=F)
	}
differences = list()
for (i in 1:100)
	{
		tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/PopCountValues_",i,".csv"), head=T)
		differences[[i]] = tab[,"endPopCount"]-tab[,"startPopCount"]
	}
dev.new(width=8, height=6); par(oma=c(1,1,1,1), mar=c(3,3,1,1), col="gray30", lwd=0.2)
plot(density(differences[[1]][!is.na(differences[[1]])]), lwd=0.1, col="red", main=NA, frame=F, axes=T, xlim=c(-750,750), ylim=c(0,0.007))
for (i in 2:length(differences)) lines(density(differences[[i]][!is.na(differences[[i]])]), col="red")
abline(v=0)

# 7. Extraction and comparison of the Mexican state at the start/end nodes (based on the continuous phylogeographic inference)

analyses = c("25092024")
admin1 = shapefile("Mex_GADM_shapefiles/Mexico_admin1.shp"); Mexico_state_and_city = subset(admin1, NAME_1%in%c("México","Distrito Federal"))
for (h in 1:length(analyses))
	{
		for (i in 1:100)
			{
				tab1 = read.csv(paste0("All_tree_extraction_files/TreeTime_",analyses[h],"/All_clades_ext/TreeExtractions_",i,".csv"), head=T)
				tab2 = matrix(nrow=dim(tab1)[1], ncol=2); colnames(tab2) = c("startState", "endState")
				for (j in 1:length(Mexico_state_and_city@polygons))
					{
						pol = Mexico_state_and_city@polygons[[j]]@Polygons[[1]]
						pip = point.in.polygon(tab1[,"startLon"], tab1[,"startLat"], pol@coords[,"x"], pol@coords[,"y"])
						tab2[which(pip==1),"startState"] = Mexico_state_and_city@data[j,"NAME_1"]
						pip = point.in.polygon(tab1[,"endLon"], tab1[,"endLat"], pol@coords[,"x"], pol@coords[,"y"])
						tab2[which(pip==1),"endState"] = Mexico_state_and_city@data[j,"NAME_1"]
					}
				write.csv(tab2, paste0("All_tree_extraction_files/TreeTime_",analyses[h],"/All_clades_ext/StartEndStates_",i,".csv"), row.names=F, quote=F)
			}
	}

cols1 = list(); cols2 = list()
cols1[[1]] = rgb(222,67,39,255,maxColorValue=255); cols2[[1]] = rgb(222,67,39,100,maxColorValue=255) # red
cols1[[2]] = rgb(250,165,33,255,maxColorValue=255); cols2[[2]] = rgb(250,165,33,100,maxColorValue=255) # orange
cols1[[3]] = rgb(70,118,187,255,maxColorValue=255); cols2[[3]] = rgb(70,118,187,100,maxColorValue=255) # blue
cols1[[4]] = rgb(106,61,154,255,maxColorValue=255); cols2[[4]] = rgb(106,61,154,100,maxColorValue=255) # purple
cols1[[5]] = rgb(190,190,190,255,maxColorValue=255); cols2[[5]] = rgb(190,190,190,100,maxColorValue=255) # light grey
cols1[[6]] = rgb(110,110,110,255,maxColorValue=255); cols2[[6]] = rgb(110,110,110,100,maxColorValue=255) # dark grey
transitions_events_DF_to_DF_list = list(); transitions_events_MS_to_MS_list = list()
transitions_events_DF_to_MS_list = list(); transitions_events_MS_to_DF_list = list()
proportions_leaving_DF_to_MS_list = list(); proportions_leaving_MS_to_DF_list = list()
for (h in 1:length(analyses))
	{
		transitions_events_DF_to_DF = rep(NA, 100); transitions_events_MS_to_MS = rep(NA, 100)
		transitions_events_DF_to_MS = rep(NA, 100); transitions_events_MS_to_DF = rep(NA, 100)
		proportions_leaving_DF_to_MS = rep(NA, 100); proportions_leaving_MS_to_DF = rep(NA, 100)
		for (i in 1:100)
			{
				tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analyses[h],"/All_clades_ext/StartEndStates_",i,".csv"), head=T)
				transitions_events_DF_to_DF[i] = length(which((tab[,"startState"]=="Distrito Federal")&(tab[,"endState"]=="Distrito Federal")))
				transitions_events_MS_to_MS[i] = length(which((tab[,"startState"]=="México")&(tab[,"endState"]=="México")))
				transitions_events_DF_to_MS[i] = length(which((tab[,"startState"]=="Distrito Federal")&(tab[,"endState"]=="México")))
				transitions_events_MS_to_DF[i] = length(which((tab[,"startState"]=="México")&(tab[,"endState"]=="Distrito Federal")))
				proportions_leaving_DF_to_MS[i] = length(which((tab[,"startState"]=="Distrito Federal")&(tab[,"endState"]=="México")))
				proportions_leaving_DF_to_MS[i] = proportions_leaving_DF_to_MS[i]/length(which(tab[,"startState"]=="Distrito Federal"))
				proportions_leaving_MS_to_DF[i] = length(which((tab[,"startState"]=="México")&(tab[,"endState"]=="Distrito Federal")))
				proportions_leaving_MS_to_DF[i] = proportions_leaving_MS_to_DF[i]/length(which(tab[,"startState"]=="México"))
			}
		cat(analyses[h],":\n",sep="")
		vS = transitions_events_DF_to_DF; meanV = round(mean(vS)); hpd95 = round(HDInterval::hdi(vS))[c("lower","upper")]
		cat("\tTransitions CDMX to CDMX: ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
		vS = transitions_events_MS_to_MS; meanV = round(mean(vS)); hpd95 = round(HDInterval::hdi(vS))[c("lower","upper")]
		cat("\tTransitions MEX to MEX: ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
		vS = transitions_events_DF_to_MS; meanV = round(mean(vS)); hpd95 = round(HDInterval::hdi(vS))[c("lower","upper")]
		cat("\tTransitions CDMX to MEX: ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
		vS = transitions_events_MS_to_DF; meanV = round(mean(vS)); hpd95 = round(HDInterval::hdi(vS))[c("lower","upper")]
		cat("\tTransitions MEX to CDMX: ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
		vS = proportions_leaving_DF_to_MS; meanV = round(mean(vS),3); hpd95 = round(HDInterval::hdi(vS),3)[c("lower","upper")]
		cat("\tProportion from CDMX leaving CDMX: ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n",sep="")
		vS = proportions_leaving_MS_to_DF; meanV = round(mean(vS),3); hpd95 = round(HDInterval::hdi(vS),3)[c("lower","upper")]
		cat("\tProportion from MEX leaving MEX: ",meanV,", 95% HPD = [",hpd95[1],", ",hpd95[2],"]\n\n",sep="")
		transitions_events_DF_to_DF_list[[h]] = transitions_events_DF_to_DF; transitions_events_MS_to_MS_list[[h]] = transitions_events_MS_to_MS
		transitions_events_DF_to_MS_list[[h]] = transitions_events_DF_to_MS; transitions_events_MS_to_DF_list[[h]] = transitions_events_MS_to_DF
		proportions_leaving_DF_to_MS_list[[h]] = proportions_leaving_DF_to_MS; proportions_leaving_MS_to_DF_list[[h]] = proportions_leaving_MS_to_DF
	}

minYear = 9999; maxYear = -9999
for (i in 1:100)
	{
		tab = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/TreeExtractions_",i,".csv"), head=T)
		if (minYear > min(tab[,"startYear"])) minYear = min(tab[,"startYear"])
		if (maxYear < max(tab[,"endYear"])) maxYear = max(tab[,"endYear"])
	}
timeSlice = 365/12; timeIntervals = floor((maxYear-minYear)/(timeSlice/365))
transitions_CDMX_MEX_counts = matrix(nrow=timeIntervals, ncol=100)
transitions_MEX_CDMX_counts = matrix(nrow=timeIntervals, ncol=100)
transitions_total_counts = matrix(nrow=timeIntervals, ncol=100)
for (i in 1:100)
	{
		tab1 = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/TreeExtractions_",i,".csv"), head=T)
		tab2 = read.csv(paste0("All_tree_extraction_files/TreeTime_",analysis,"/All_clades_ext/StartEndStates_",i,".csv"), head=T)
		for (j in 1:timeIntervals)
			{
				t1 = minYear+((j-1)*(timeSlice/365)); t2 = minYear+(j*(timeSlice/365))
				indices = which((tab1[,"endYear"]>t1)&(tab1[,"endYear"]<t2))
				transitions_CDMX_MEX_counts[j,i] = length(which((tab2[indices,"startState"]=="Distrito Federal")&(tab2[indices,"endState"]=="México")))
				transitions_MEX_CDMX_counts[j,i] = length(which((tab2[indices,"startState"]=="México")&(tab2[indices,"endState"]=="Distrito Federal")))
				transitions_total_counts[j,i] = dim(tab2[indices,])[1]
			}
	}
transitions_CDMX_MEX_frequencies = matrix(nrow=timeIntervals, ncol=4)
transitions_MEX_CDMX_frequencies = matrix(nrow=timeIntervals, ncol=4)
for (j in 1:timeIntervals)
	{
		transitions_CDMX_MEX_frequencies[j,1] = minYear+((j-1)*(timeSlice/365))+((timeSlice/365)/2)
		transitions_MEX_CDMX_frequencies[j,1] = minYear+((j-1)*(timeSlice/365))+((timeSlice/365)/2)
		vS = transitions_CDMX_MEX_counts[j,]/transitions_total_counts[j,]
		transitions_CDMX_MEX_frequencies[j,2] = median(vS)
		transitions_CDMX_MEX_frequencies[j,3] = HDInterval::hdi(vS)[1]
		transitions_CDMX_MEX_frequencies[j,4] = HDInterval::hdi(vS)[2]
		vS = transitions_MEX_CDMX_counts[j,]/transitions_total_counts[j,]
		transitions_MEX_CDMX_frequencies[j,2] = median(vS)
		transitions_MEX_CDMX_frequencies[j,3] = HDInterval::hdi(vS)[1]
		transitions_MEX_CDMX_frequencies[j,4] = HDInterval::hdi(vS)[2]
	}
indices = 160:434 # for analysis = "25092024"
minY = min(cbind(transitions_CDMX_MEX_frequencies[indices,3],transitions_MEX_CDMX_frequencies[indices,3]), na.rm=T)
maxY = max(cbind(transitions_CDMX_MEX_frequencies[indices,4],transitions_MEX_CDMX_frequencies[indices,4]), na.rm=T)
dev.new(width=10.2, height=2.0); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.2,0.5,1,0), lwd=0.5, col="gray30")
plot(transitions_CDMX_MEX_frequencies[indices,1], transitions_CDMX_MEX_frequencies[indices,2], col=NA, type="l", axes=F, ann=F, xlim=c(2000,maxYear), ylim=c(minY,0.25))
xx_l = c(transitions_MEX_CDMX_frequencies[indices,1],rev(transitions_MEX_CDMX_frequencies[indices,1]))
yy_l = c(transitions_MEX_CDMX_frequencies[indices,3],rev(transitions_MEX_CDMX_frequencies[indices,4]))
polygon(xx_l, yy_l, col=cols2[[4]], border=0) 
xx_l = c(transitions_CDMX_MEX_frequencies[indices,1],rev(transitions_CDMX_MEX_frequencies[indices,1]))
yy_l = c(transitions_CDMX_MEX_frequencies[indices,3],rev(transitions_CDMX_MEX_frequencies[indices,4]))
polygon(xx_l, yy_l, col=cols2[[3]], border=0)
lines(transitions_MEX_CDMX_frequencies[indices,1], transitions_MEX_CDMX_frequencies[indices,2], lwd=1, col=cols1[[4]])
lines(transitions_CDMX_MEX_frequencies[indices,1], transitions_CDMX_MEX_frequencies[indices,2], lwd=1, col=cols1[[3]])
axis(side=1, lwd.tick=0.5, cex.axis=0.8, mgp=c(0,0.17,0), lwd=0.5, tck=-0.035, col.tick="gray30", col.axis="gray30", col="gray30")
axis(pos=1999.99, side=2, lwd.tick=0.5, cex.axis=0.8, mgp=c(0,0.30,0), lwd=0.5, tck=-0.035, col.tick="gray30", col.axis="gray30", col="gray30")
title(ylab="Frequency", cex.lab=0.8, mgp=c(-0.2,0,0), col.lab="gray30")
indices = 280:434 # for analysis = "25092024"
minY = min(cbind(transitions_CDMX_MEX_frequencies[indices,3],transitions_MEX_CDMX_frequencies[indices,3]), na.rm=T)
maxY = max(cbind(transitions_CDMX_MEX_frequencies[indices,4],transitions_MEX_CDMX_frequencies[indices,4]), na.rm=T)
dev.new(width=10.2, height=2.0); par(mgp=c(0,0,0), oma=c(0.8,0.8,0,0), mar=c(1.2,0.5,1,0), lwd=0.5, col="gray30")
plot(transitions_CDMX_MEX_frequencies[indices,1], transitions_CDMX_MEX_frequencies[indices,2], col=NA, type="l", axes=F, ann=F, xlim=c(2010,maxYear), ylim=c(minY,0.25))
xx_l = c(transitions_MEX_CDMX_frequencies[indices,1],rev(transitions_MEX_CDMX_frequencies[indices,1]))
yy_l = c(transitions_MEX_CDMX_frequencies[indices,3],rev(transitions_MEX_CDMX_frequencies[indices,4]))
polygon(xx_l, yy_l, col=cols2[[4]], border=0) 
xx_l = c(transitions_CDMX_MEX_frequencies[indices,1],rev(transitions_CDMX_MEX_frequencies[indices,1]))
yy_l = c(transitions_CDMX_MEX_frequencies[indices,3],rev(transitions_CDMX_MEX_frequencies[indices,4]))
polygon(xx_l, yy_l, col=cols2[[3]], border=0)
lines(transitions_MEX_CDMX_frequencies[indices,1], transitions_MEX_CDMX_frequencies[indices,2], lwd=1, col=cols1[[4]])
lines(transitions_CDMX_MEX_frequencies[indices,1], transitions_CDMX_MEX_frequencies[indices,2], lwd=1, col=cols1[[3]])
axis(side=1, lwd.tick=0.5, cex.axis=0.8, mgp=c(0,0.17,0), lwd=0.5, tck=-0.035, col.tick="gray30", col.axis="gray30", col="gray30")
axis(pos=2009.99, side=2, lwd.tick=0.5, cex.axis=0.8, mgp=c(0,0.30,0), lwd=0.5, tck=-0.035, col.tick="gray30", col.axis="gray30", col="gray30")
title(ylab="Frequency", cex.lab=0.8, mgp=c(-0.2,0,0), col.lab="gray30")

