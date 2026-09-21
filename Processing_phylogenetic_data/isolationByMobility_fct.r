isolationByMobility_fct = function(localTreesDirectory="", nberOfExtractionFiles=1, mobility_index, outputName=gsub(" ","_",date())) {

	meanStatistics = matrix(nrow=(nberOfExtractionFiles), ncol=2)
	for (t in 1:nberOfExtractionFiles)
		{
			data = read.csv(paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep=""), h=T)
			nberOfConnections = dim(data)[1]
			if (t == 1)
				{
					minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
					maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
					minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
					maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
					minStartYear = min(data[,"startYear"])
		    		minEndYear = min(data[,"startYear"])
		    		maxEndYear = max(data[,"endYear"])
				}	else	{
					if (minLon > min(min(data[,"endLon"]),min(data[,"startLon"]))) minLon = min(min(data[,"endLon"]),min(data[,"startLon"]))
					if (maxLon < max(max(data[,"endLon"]),max(data[,"startLon"]))) maxLon = max(max(data[,"endLon"]),max(data[,"startLon"]))
					if (minLat > min(min(data[,"endLat"]),min(data[,"startLat"]))) minLat = min(min(data[,"endLat"]),min(data[,"startLat"]))
					if (maxLat < max(max(data[,"endLat"]),max(data[,"startLat"]))) maxLat = max(max(data[,"endLat"]),max(data[,"startLat"]))
					if (minStartYear > min(data[,"startYear"])) minStartYear = min(data[,"startYear"])
					if (minEndYear > min(data[,"endYear"])) minEndYear = min(data[,"endYear"])
					if (maxEndYear < max(data[,"endYear"])) maxEndYear = max(data[,"endYear"])
				}
		}
	for (t in 1:nberOfExtractionFiles)
		{
			data = read.csv(paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep=""), h=T)
			distances = matrix(0, nrow=dim(data)[1], ncol=1); colnames(distances) = c("greatCircleDist_km")
			for (i in 1:nberOfConnections)
				{
					x1 = cbind(data[i,"startLon"],data[i,"startLat"])
	 				x2 = cbind(data[i,"endLon"],data[i,"endLat"])
	 				distances[i] = rdist.earth(x1, x2, miles=F, R=NULL)
				}
			if (length(grep("greatCircleDist_km",colnames(data))) > 0)
				{
					data[,"greatCircleDist_km"] = distances
				}	else	{
					data = cbind(data,distances)
				}
			tipNodeIndices = which(!data[,"node2"]%in%data[,"node1"])
			distTree = matrix(nrow=length(tipNodeIndices), ncol=length(tipNodeIndices))
			for (k in 2:dim(distTree)[1])
				{
					for (l in 1:(k-1))
						{
							index1 = tipNodeIndices[k]
							index2 = tipNodeIndices[l]
							indices1 = index1; root = FALSE
							while (root == FALSE)
								{	
									if (data[indices1[length(indices1)],"node1"]%in%data[,"node2"])
										{
											indices1 = c(indices1, which(data[,"node2"]==data[indices1[length(indices1)],"node1"]))
										}	else	{
											root = TRUE
										}
								}
							indices2 = index2; root = FALSE
							while (root == FALSE)
								{	
									if (data[indices2[length(indices2)],"node1"]%in%data[,"node2"])
										{
											indices2 = c(indices2, which(data[,"node2"]==data[indices2[length(indices2)],"node1"]))
										}	else	{
											root = TRUE
										}
								}
							indices3 = indices1[which(indices1%in%indices2)]; patristic_dis = NULL
							if (length(indices3) == 0)
								{
									patristic_dis = sum(data[c(indices1,indices2),"length"])
								}	else	{
									patristic_dis = sum(data[c(indices1[which(!indices1%in%indices3)],indices2[which(!indices2%in%indices3)]),"length"])
								}
							distTree[k,l] = patristic_dis; distTree[l,k] = patristic_dis
						}
				}
			mobility = matrix(nrow=dim(distTree)[1], ncol=dim(distTree)[2])
			for (k in 1:dim(mobility)[1])
				{
					for (l in 1:dim(mobility)[1])
						{
							startLocation = data[tipNodeIndices[k],"endAGEB"]; endLocation = data[tipNodeIndices[l],"endAGEB"]
							index = which((mobility_index[,"origen"]==startLocation)&(mobility_index[,"destino"]==endLocation))
							if (length(index) == 1) mobility[k,l] = mobility_index[index,"peso_mediana"]
						}
				}
			distTree = distTree[!is.na(mobility)]; mobility = mobility[!is.na(mobility)]
			meanStatistics[t,1] = cor(distTree, mobility, method="pearson") # r_P #1
			meanStatistics[t,2] = cor(distTree, log(mobility+1), method="pearson") # r_P #2
		}
	colnames(meanStatistics) = c("IBM_signal_rP1","IBM_signal_rP2")
	if (nberOfExtractionFiles > 1)
		{
			write.table(meanStatistics, file=paste(outputName,"_isolation-by-mobility_signals.txt",sep=""), quote=F, row.names=F, sep="\t")
		}
	if (nberOfExtractionFiles == 1)
		{
			return(cbind(meanStatistics[1,1],meanStatistics[1,2]))
		}
}
