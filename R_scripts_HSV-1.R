

require(igraph)				#load igraph

#***********SCR_1******************start

#importing the data as a dataframe
library(gdata)                   # load gdata package, so as to import the data as a dataframe 

network_data <- read.xls("C:/Users/micheal/Downloads/Project/edgeweights.xlsx")

#converting the dataframe into an igraph object/graph called 'g'
g <- graph.data.frame(network_data, directed=FALSE, vertices=NULL)


#append orf attributes to the proteins in the network.
g<-set.vertex.attribute(g, 'orf', c("P10226","P10188","P04295","P10212","P04289","P08543","P06492","P10233","P06491","P04487","P10221","P10220","P04293","P10210","P08392","P10193","P04296","P04294","P68331","P08393","P10186","P10187","P10191","P04288","P04290","P04291","P10200","P10201","P10202","P10204","P10205","P03176","P10209","P10211","P10215","P10217","P10219","P32888","P10224","P10227","P10229","P10230","P10231","O09800","P10236","P10238","P10240","P04485","P06485","P04488","P10234","P10218","P10192","P10216","P06481","P10235","P06486","P10239","P06477","P06487","P10228","P04413","P10225","Q69091","P10190","O09802","P06480","P06484","P10189","P10208","P10185"), c("UL42","UL4","UL15","UL28","UL11","UL39","UL48","UL49","UL19","US11","UL37","UL36","UL30","UL26","RS1","UL9","UL29","UL12","UL53","RL2","UL2","UL3","UL7","UL10","UL13","UL14","UL16","UL17","UL18","UL20","UL21","UL23","UL25","UL27","UL31","UL33","UL35","UL38","UL40","UL43","UL45","UL46","UL47","UL49.5", "UL52","UL54","UL56","US1","US2","US8","UL50","UL34","UL8","UL32","US9","UL51","US10","UL55","UL22","US7","UL44","US3","UL41","US6","UL6","US8.5","US5","US4","UL5","UL24","UL1"))

#***********SCR_1******************stop




#***********SCR_2******************start
#testing the structure of the graph

#remove multiple edges and loops
simplify(g, remove.multiple=TRUE, remove.loops=TRUE, edge.attr.comb=getIgraphOpt("edge.attr.comb"))

#obtain the vertex count
vcount(g_m)

#obtain the edge count
ecount(g_m)

#check to confirm that it's undirected; expect FALSE
is.directed(g)

#check if the graph has weights
is.weighted(g)

#check if the graph has the 'name' and 'orf' attributes
list.vertex.attributes(g)

#***********SCR_2******************stop










###########################################################
##CLUSTERING OF THE NETWORK
###########################################################

#***********SCR_4******************start
##
#PRE-TEST TO DETERMINE OPTIMAL PARAMETERS
##
#For the Walktrap algorithm
optimise_param.WT<- function(graph){
	modularities<- NULL
	for (i in 1:60){
		communities_WT <- walktrap.community(graph, weights = E(graph)$weight, steps = i, merges=TRUE, modularity = TRUE, membership = TRUE)
		modularity<- modularity(communities_RW)
		modularities<- c(modularities, modularity)
	}
	return (modularities)
}

#For the Infomap algorithm
optimise_param.IFM<- function(graph){
	modularities<- NULL
	for (i in 1:60){
		communities_IFM <- infomap.community(graph, e.weights=E(graph)$weight, v.weights=NULL, nb.trials=i, modularity=TRUE)
		modularity<- modularity(communities_IFM)
		modularities<- c(modularities, modularity)
	}
	return (modularities)
}
#***********SCR_4******************stop





#***********SCR_3******************start
#########
#Clustering with all the algorithms
#########

#Random walk
communities_RW <- walktrap.community(g, weights = E(g)$weight, steps = 11, merges=TRUE, modularity = TRUE, membership = TRUE)

#edge.betweenness comm. detection
communities_EB <- edge.betweenness.community(g, weights=E(g)$weight, directed=FALSE, edge.betweenness=TRUE, merges=TRUE, bridges=TRUE, modularity=TRUE, membership=TRUE)

#Fast Greedy
communities_FG <- fastgreedy.community(g, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(g)$weight)

#Multi-level comm. detctn
communities_ML <- multilevel.community(g, weights=E(g)$weight)

#Label propagation comm. detctn
communities_LP<- label.propagation.community(g, weights=E(g)$weight, initial=NULL)  #no fixed state for any label

#Leading eigen vector community
communities_LE <- leading.eigenvector.community(g, weights=E(g)$weight, start=NULL, options=igraph.arpack.default, callback=NULL, extra=NULL, env=parent.frame())

#Infomap comm. detection
communities_IFM <- infomap.community(g, e.weights=E(g)$weight, v.weights=NULL, nb.trials=10, modularity=TRUE)

#***********SCR_3******************stop



#***********SCR_5******************start
#Calculating modularities

modularity(communities_RW)
modularity(communities_EB)
modularity(communities_FG)
modularity(communities_ML)
modularity(communities_LP)
modularity(communities_LE)
modularity(communities_IFM)

#***********SCR_5******************stop




#################
#FUNCTION for showing the member-proteins of every cluster in a communities object by the orf name
############
#
# This function was used throughout the study to obtain the identity of the member-proteins of each cluster in a communities-object
#
#Function produces a list of lists, that is, a list of lists of proteins in each cluster
#	output looks like this 	[[1]]
#					[1] "UL4"  "US11" "UL56"....
#					[[2]]
#					[1] "UL15" "RS1"  "UL7" ...
#					[[3]]
#					[1] "UL39" "UL13" "US3" ...
#					....
#					...
#
# The function is called like this;
# my_list <- show_clusters(communities-object, "algorithm_name")
#


show_clusters <- function(communities, algo_name){			#input parameters are the communities-object of an algorithm, and the algorithm name
	list_of_clusters <- NULL						#create the list of lists
	for (i in 1:length(sizes(communities))){				#'sizes(communities)' gives a vector of the sizes of the communities/cluster communities-object, in order
		cluster <- NULL							#for every community/cluster create an empty list of the proteins in it
		for (j in 1:length(communities$membership)){		#'communities$membership' shows the communities to which each protein/node in the graph belongs to
			if (communities$membership[j] == names(sizes(communities))[i]){	#for every protein, get the community to which it belongs, if its the same as the name of
				cluster <- c(cluster, V(g)[j]$orf)					#	the community in focus, then add the protein to this community's list
			}
		}
		if (length(cluster) >= 3){					# only interested in clusters with 3 or more members - filtering out mere binary interactions
			list_of_clusters<- c(list_of_clusters, list(cluster))	# put this community's protein list into the full list of proteins
			attr(list_of_clusters, "algorithm")<- algo_name		# for presentation purposes, assign the algorithm name (as an attribute) to be shown when the full-list is printed 
		}
	}
	return(list_of_clusters)
}





##############################################################################
#			EDGE-WEIGHT PERTURBATION
##############################################################################
#
#
#For edge-weight perturbation, the first step is to compute the perturbation measure; the value by which the perturb the weight by. This is done using the function 'compute_ptb.measure'
#	The function is called like this 
#	ptb.measure<- compute_ptb.measure(g)
#Next, the function 'wtpertubator' is written to do the actual perturbation of the graph-weights;
#	This function first randomly chooses a number of edges equal to half of those in the graph
#	For each of the chosen edges, its weight is incremeted by the perturbation measure	
#
#
#

######## COMPUTING PERTURBATION MEASURE

#***********SCR_6******************start
compute_ptb.measure<- function(graph){
	quartiles<- quantile(E(g)$weight, prob=c(0.25, 0.75))			# obtain the quartiles
	lower_quartile<- quartiles[[1]]; upper_quartile<- quartiles[[2]]	# obtain the lower and upper quartile values
	middle_values<- NULL											# create empty vector of values which are between the quartiles ('middle values')
	for (i in 1:ecount(graph)){										# for every edge in the network,
		if (E(graph)$weight[i]>lower_quartile && E(graph)$weight[i]<upper_quartile){	# check if it's between the lower and upper quartiles
			middle_values<- c(middle_values, E(graph)$weight[i])				# if it is, assign it to the vector of values between the quartiles
		}
	}														# REPEAT THE ABOVE 2 STEPS FOR EVERY NODE IN EDGE IN THE NETWORK

	ptb.measure<- sd(middle_values)									# get the standard deviation of the 'middle values'
	return(ptb.measure)											# return the standard deviation as the perturbation measure
}	
#***********SCR_6******************stop



#***********SCR_7******************start
######### FUNCTION TO PERTURB THE EDGE-WEIGHTS

wtpertubator <- function(graph) {
	ptb.measure<- compute_ptb.measure(graph)								# first obtain the perturbation measure for the given graph, using 'compute_ptb.measure'
	vec<- sample(1:ecount(graph), ecount(graph)/2, replace=F)					# randomly obtain the edges whose weights are to be perturbed
	for (i in 1:ecount(graph)) {										# Then, for every node in the network, 
		if (i %in% vec){											# check if its's one of those to be perturbed
			E(graph)[i]$weight <- E(graph)[i]$weight+ptb.measure				# if it is, then add the perturbation measure to it
		}
	}
	return (graph)												# return the perturbedd graph
}



############################################################################
##
#### 		edge-weight perturbation with FASTGREEDY
##
############################################################################
#In this test, edge-weights are perturbed 1000 times, and on every one of those times, the change in every original cluster is obtained by calculating a similarity score(jaccard coefficient)
#This function is supposed to produce a list of similarity scores for every original cluster in the network, from which the median similarity will be computed
#the function will thus produce a n lists, where n is the number of original clusters; 
#For convenience while accessing the data/scores, the n lists were put in a list; the function produces a list of lists.


#The steps followed are outlined below; the numbering follows the numbering in the code
##1# Clustering the intact/original network using the Fastgreedy algorithm to obtain the original clusters (cluster.orig)
##2# The communities-object produced by igraph does not show the identity of the proteins in the different clusters; to extract the proteins(by name) into their clusters
#	the function 'show_clusters was used'; it produces a list of clusters, each with a list of proteins in that cluster - a list of lists
#	it looks like this 	[[1]]
#					[1] "UL4"  "US11" "UL56"....
#					[[2]]
#					[1] "UL15" "RS1"  "UL7" ...
#					[[3]]
#					[1] "UL39" "UL13" "US3" ...
#					....
#					...
#	Create the list of lists of all similarity/jaccard scores
##3# For each original cluster; 
#	Create its list of jaccard scores between it and every cluster.ptb
##4# For every one of 1000 times, 
##5# Use the function 'wtpertubator' to perturb the weights of randomly selected half of the edges in the graph, producing a new graph (graph.ptb) in the process.
##6# Compute the cluster/communities in the new graph (graph.ptb)
##7# Obtain the member-proteins of each of the (new) clusters (cluster.ptb's) in the perturbed graph
#
#   At this point, it was necessary to identify the cluster.orig in focus from all the produced cluster.ptb's; this was done by finding which the cluster.ptb's 
#	has the highest numbers of common/similar members with cluster.orig. It is this that was taken to be the 'transformed cluster.orig', that is, the corresponding cluster.ptb
#	When the corresponding cluster.ptb is obtained, the jaccard coefficient is calculated (between cluster.orig and its corresponding cluster.ptb)
	Steps #8# to #10 achieve this
##8# Initiate the vector to store the number of common/similar members between cluster.orig and each cluster.ptb
##9# For each cluster.ptb,
##10# Initiate count of similar proteins
##11# For each member of cluster.ptb, compare it to each member of cluster.orig and if they are the same, increase the count of similar proteins between 
#		cluster.orig and this cluster.ptb by 1
##12# Store the number of similar proteins for this cluster.ptb in a vector
#
#	REPEAT STEPS 10 to 12 FOR ALL CLUSTER.PTB's
#
##13# Since the similarity counts for each cluster.ptb are stored in a vector, i need to be able to identify the identity of the cluster.ptb with the maximum count; for that reason,
#		i attach a 'names' attribute to each of the count. The attribute is the position of the cluster.ptb in the vector 
##14# To get the identity of the cluster.ptb with the highest number of similar proteins as cluster.orig, that is, the corresponding cluster.ptb, i get the attribute(forced to be a number)
#		of highest number in vector of similarity counts
##15# Then i calculate the Jaccard coefficient between cluster.orig and cluter.ptb using the formula; Jaccard=members in both/(members in A_not_in_B + members_in_B_not_A + members in both))
##16# This jaccard score is stored in the list of scores created just before step 3
#
#	REPEAT STEPS #5 to #16 999 times
#
##17# Store the 1000 jaccard scores as a list in the list 'all_similarity_scores'
#
##	REPEAT STEPS #4 to #17 for every other original cluster
#
##18#	Return the list of lists that is 'all_similarity_scores'
#
#The list of jaccard coefficients for, e.g cluster one, can be assessed like this
#scores_list.FG<- weight_perturbator.jaccard.FG(graph)
#scores_list.FG[[1]]
#
#

#############
# THE FUNCTION
#############

weight_perturbator.jaccard.FG<- function(graph){
#1#
	communities_FG.orig <- fastgreedy.community(graph, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(graph)$weight)
#2#
	clusters.orig<- show_clusters(communities_FG.orig,maimail "fastgreedy")
	all_similarity_scores<- NULL
#3#
	for (n in 1:length(clusters.orig)) {
		jaccard_scores<- NULL
#4#
		for (i in 1:1000) {
#5#
			graph.ptb <- wtpertubator(g)
#6#
			communities_FG.ptb <- fastgreedy.community(graph.ptb, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(graph.ptb)$weight)
#7#
			clusters.ptb<- show_clusters(communities_FG.ptb, "fastgreedy.ptb")
#8#
			simlty_counts<- NULL
#9#
			for (j in 1:length(clusters.ptb)){
#10#
				sim_count=0
#11#
				for (k in 1:length(clusters.ptb[[j]]))   {
					for (m in 1:length(clusters.orig[[n]]))   {
						if (clusters.orig[[n]][m] == clusters.ptb[[j]][k]){
							sim_count=sim_count+1		
						}
					}
				}
#12#
				simlty_counts<- c(simlty_counts, sim_count)
			}
#13#
			attr(simlty_counts, "names")<- c(1:length(clusters.ptb))
#14#
			similar_cluster<- as.numeric(attributes(simlty_counts[match(c(max(sapply(simlty_counts, max))), simlty_counts)])[[1]])
#15#
			jaccard<- max(sapply(simlty_counts, max))/(   length(clusters.orig[[n]]) + length(clusters.ptb[[similar_cluster]]) + max(sapply(simlty_counts, max))   )
#16#
			jaccard_scores<- c(jaccard_scores, jaccard)
		}
#17#
		all_similarity_scores <- c(all_similarity_scores, list(jaccard_scores))
	}
#18#
	return(all_similarity_scores)
}







#########################################################################
#
######		edge-weight perturbation with MULTI-LEVEL
#
#########################################################################

#In this test, edge-weights are perturbed 1000 times, and on every one of those times, the change in every original cluster is obtained by calculating a similarity score(jaccard coefficient)
#This function is supposed to produce a list of similarity scores for every original cluster in the network, from which the median similarity will be computed
#the function will thus produce a n lists, where n is the number of original clusters; 
#For convenience while accessing the data/scores, the n lists were put in a list; the function produces a list of lists.


#The steps followed are outlined below; the numbering follows the numbering in the code
##1# Clustering the intact/original network using the Multi-level algorithm to obtain the original clusters (cluster.orig)
##2# The communities-object produced by igraph does not show the identity of the proteins in the different clusters; to extract the proteins(by name) into their clusters
#	the function 'show_clusters was used'; it produces a list of clusters, each with a list of proteins in that cluster - a list of lists
#	it looks like this 	[[1]]
#					[1] "UL4"  "US11" "UL56"....
#					[[2]]
#					[1] "UL15" "RS1"  "UL7" ...
#					[[3]]
#					[1] "UL39" "UL13" "US3" ...
#					....
#					...
#	Create the list of lists of all similarity/jaccard scores
##3# For each original cluster; 
#	Create its list of jaccard scores between it and every cluster.ptb
##4# For every one of 1000 times, 
##5# Use the function 'wtpertubator' to perturb the weights of randomly selected half of the edges in the graph, producing a new graph (graph.ptb) in the process.
##6# Compute the cluster/communities in the new graph (graph.ptb)
##7# Obtain the member-proteins of each of the (new) clusters (cluster.ptb's) in the perturbed graph
#
#   At this point, it was necessary to identify the cluster.orig in focus from all the produced cluster.ptb's; this was done by finding which the cluster.ptb's 
#	has the highest numbers of common/similar members with cluster.orig. It is this that was taken to be the 'transformed cluster.orig', that is, the corresponding cluster.ptb
#	When the corresponding cluster.ptb is obtained, the jaccard coefficient is calculated (between cluster.orig and its corresponding cluster.ptb)
	Steps #8# to 10 #achieve this
##8# Initiate the vector to store the number of common/similar members between cluster.orig and each cluster.ptb
##9# For each cluster.ptb,
##10# Initiate count of similar proteins
##11# For each member of cluster.ptb, compare it to each member of cluster.orig and if they are the same, increase the count of similar proteins between 
#		cluster.orig and this cluster.ptb by 1
##12# Store the number of similar proteins for this cluster.ptb in a vector
#
#	REPEAT STEPS 10 to 12 FOR ALL CLUSTER.PTB's
#
##13# Since the similarity counts for each cluster.ptb are stored in a vector, i need to be able to identify the identity of the cluster.ptb with the maximum count; for that reason,
#		i attach a 'names' attribute to each of the count. The attribute is the position of the cluster.ptb in the vector 
##14# To get the identity of the cluster.ptb with the highest number of similar proteins as cluster.orig, that is, the corresponding cluster.ptb, i get the attribute(forced to be a number)
#		of highest number in vector of similarity counts
##15# Then i calculate the Jaccard coefficient between cluster.orig and cluter.ptb using the formula; Jaccard=members in both/(members in A_not_in_B + members_in_B_not_A + members in both))
##16# This jaccard score is stored in the list of scores created just before step 3
#
#	REPEAT STEPS #5 to #16 999 times
#
##17# Store the 1000 jaccard scores as a list in the list 'all_similarity_scores'
#
##	REPEAT STEPS #4 to #17 for every other original cluster
#
##18#	Return the list of lists that is 'all_similarity_scores'
#
#
#The list of jaccard coefficients for, e.g cluster one, can be assessed like this
#scores_list.ML<- weight_perturbator.jaccard.ML(graph)
#scores_list.ML[[1]]
#
#

#######
# THE FUNCTION
###

weight_perturbator.jaccard.ML<- function(graph){
#1#
	communities_ML.orig <- multilevel.community(graph, weights=E(graph)$weight)
#2#
	clusters.orig<- show_clusters(communities_ML.orig, "multi-level")
	all_similarity_scores<- NULL
#3#
	for (n in 1:length(clusters.orig)) {
		jaccard_scores<- NULL
#4#
		for (i in 1:100) {
#5#
			graph.ptb <- wtpertubator(graph)
#6#
			communities_ML.ptb <- multilevel.community(graph.ptb, weights=E(graph.ptb)$weight)
#7#
			clusters.ptb<- show_clusters(communities_ML.ptb, "multi-level.ptb")
#8#
			simlty_counts<- NULL
#9#
			for (j in 1:length(clusters.ptb)){
#10#
				sim_count=0
#11#
				for (k in 1:length(clusters.ptb[[j]]))   {
					for (m in 1:length(clusters.orig[[n]]))   {
						if (clusters.orig[[n]][m] == clusters.ptb[[j]][k]){
							sim_count=sim_count+1		
						}
					}
				}
#12#
				simlty_counts<- c(simlty_counts, sim_count)
			}
#13#
			attr(simlty_counts, "names")<- c(1:length(clusters.ptb))
#14#
			similar_cluster<- as.numeric(attributes(simlty_counts[match(c(max(sapply(simlty_counts, max))), simlty_counts)])[[1]])
#15#
			jaccard<- max(sapply(simlty_counts, max))/(   length(clusters.orig[[n]]) + length(clusters.ptb[[similar_cluster]]) + max(sapply(simlty_counts, max))   )
#16#
			jaccard_scores<- c(jaccard_scores, jaccard)
		}
#17#
		all_similarity_scores <- c(all_similarity_scores, list(jaccard_scores))
	}
#18#
	return(all_similarity_scores)
}

#***********SCR_7******************stop







#***********SCR_8******************start
###################################################################################
#		NODE REMOVAL
###################################################################################


#######FASTGREEDY
#######
#In this test, every network-nodes is removed, one every time, and on every one of the removals, the change in every original cluster is obtained by calculating a similarity score(jaccard coefficient)
#This function is supposed to produce a list of similarity scores for every original cluster in the network, from which the median similarity will be computed
#the function will thus produce a n lists, where n is the number of original clusters; 
#For convenience while accessing the data/scores, the n lists were put in a list; the function produces a list of lists.


#The steps followed are outlined below; the numbering follows the numbering in the code
##1# Clustering the intact/original network using the Fastgreedy algorithm to obtain the original clusters (cluster.orig)
##2# The communities-object produced by igraph does not show the identity of the proteins in the different clusters; to extract the proteins(by name) into their clusters
#	the function 'show_clusters was used'; it produces a list of clusters, each with a list of proteins in that cluster - a list of lists
#	it looks like this 	[[1]]
#					[1] "UL4"  "US11" "UL56"....
#					[[2]]
#					[1] "UL15" "RS1"  "UL7" ...
#					[[3]]
#					[1] "UL39" "UL13" "US3" ...
#					....
#					...
#	Create the list of lists of all similarity/jaccard scores
##3# For each original cluster; 
#	Create its list of jaccard scores between it and every cluster.rmv
##4# For every node in the network, 
##5# Delete it from the network
##6# Compute the cluster/communities in the new graph (graph.rmv)
##7# Obtain the member-proteins of each of the (new) clusters (cluster.rmv's) in the perturbed graph
#
#   At this point, it was necessary to identify the cluster.orig in focus from all the produced cluster.rmv's; this was done by finding which the cluster.rmv's 
#	has the highest numbers of common/similar members with cluster.orig. It is this that was taken to be the 'transformed cluster.orig', that is, the corresponding cluster.rmv
#	When the corresponding cluster.rmv is obtained, the jaccard coefficient is calculated (between cluster.orig and its corresponding cluster.rmv)
	Steps #8 to #10 achieve this
##8# Initiate the vector to store the number of common/similar members between cluster.orig and each cluster.rmv
##9# For each cluster.rmv,
##10# Initiate count of similar proteins
##11# For each member of cluster.rmv, compare it to each member of cluster.orig and if they are the same, increase the count of similar proteins between 
#		cluster.orig and this cluster.rmv by 1
##12# Store the number of similar proteins for this cluster.rmv in a vector
#
#	REPEAT STEPS 10 to 12 FOR ALL CLUSTER.RMV's
#
##13# Since the similarity counts for each cluster.rmv are stored in a vector, i need to be able to identify the identity of the cluster.rmv with the maximum count; for that reason,
#		i attach a 'names' attribute to each of the count. The attribute is the position of the cluster.rmv in the vector 
##14# To get the identity of the cluster.rmv with the highest number of similar proteins as cluster.orig, that is, the corresponding cluster.rmv, i get the attribute(forced to be a number)
#		of highest number in vector of similarity counts
##15# Then i calculate the Jaccard coefficient between cluster.orig and cluter.rmv using the formula; Jaccard=members in both/(members in A_not_in_B + members_in_B_not_A + members in both))
##16# This jaccard score is stored in the list of scores created just before step 3
#
#	REPEAT STEPS #5 to #16 999 times
#
##17# Store the 1000 jaccard scores as a list in the list 'all_similarity_scores'
#
##	REPEAT STEPS #4 to #17 for every other original cluster
#
##18#	Return the list of lists that is 'all_similarity_scores'
#
#
#The list of jaccard coefficients for, e.g cluster one, can be assessed like this
#scores_list.FG<- node_removal.jaccard.FG(graph)
#scores_list.FG[[1]]
#
#


#######
# THE FUNCTION
###

node_removal.jaccard.FG<- function(graph){
#1#
	communities_FG.orig <- fastgreedy.community(graph, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(graph)$weight)
#2#
	clusters.orig<- show_clusters(communities_FG.orig, "fastgreedy")
	all_similarity_scores<- NULL
#3#
	for (n in 1:length(clusters.orig)) {
		jaccard_scores<- NULL
#4#
		for (i in 1:vcount(graph)) {
#5#
			graph.rmv <- delete.vertices(graph, V(graph)[i])
#6#
			communities_FG.rmv <- fastgreedy.community(graph.rmv, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(graph.rmv)$weight)
#7#
			clusters.rmv<- show_clusters(communities_FG.rmv, "fastgreedy.rmv")
#8#
			simlty_counts<- NULL
#9#
			for (j in 1:length(clusters.rmv)){
#10#
				sim_count=0
#11#
				for (k in 1:length(clusters.rmv[[j]]))   {
					for (m in 1:length(clusters.orig[[n]]))   {
						if (clusters.orig[[n]][m] == clusters.rmv[[j]][k]){
							sim_count=sim_count+1		
						}
					}
				}
#12#
				simlty_counts<- c(simlty_counts, sim_count)
			}
#13#
			attr(simlty_counts, "names")<- c(1:length(clusters.rmv))
#14#
			similar_cluster<- as.numeric(attributes(simlty_counts[match(c(max(sapply(simlty_counts, max))), simlty_counts)])[[1]])
#15#
			jaccard<- max(sapply(simlty_counts, max))/(   length(clusters.orig[[n]]) + length(clusters.rmv[[similar_cluster]]) + max(sapply(simlty_counts, max))   )
#16#
			jaccard_scores<- c(jaccard_scores, jaccard)
		}
#17#
		all_similarity_scores <- c(all_similarity_scores, list(jaccard_scores))
	}
#18#
	return(all_similarity_scores)
}





###############################################
######Node removal for MULTI-LEVEL
##############


#In this test, every network-nodes is removed, one every time, and on every one of the removals, the change in every original cluster is obtained by calculating a similarity score(jaccard coefficient)
#This function is supposed to produce a list of similarity scores for every original cluster in the network, from which the median similarity will be computed
#the function will thus produce a n lists, where n is the number of original clusters; 
#For convenience while accessing the data/scores, the n lists were put in a list; the function produces a list of lists.


#The steps followed are outlined below; the numbering follows the numbering in the code
##1# Clustering the intact/original network using the Multi-level algorithm to obtain the original clusters (cluster.orig)
##2# The communities-object produced by igraph does not show the identity of the proteins in the different clusters; to extract the proteins(by name) into their clusters
#	the function 'show_clusters was used'; it produces a list of clusters, each with a list of proteins in that cluster - a list of lists
#	it looks like this 	[[1]]
#					[1] "UL4"  "US11" "UL56"....
#					[[2]]
#					[1] "UL15" "RS1"  "UL7" ...
#					[[3]]
#					[1] "UL39" "UL13" "US3" ...
#					....
#					...
#	Create the list of lists of all similarity/jaccard scores
##3# For each original cluster; 
#	Create its list of jaccard scores between it and every cluster.rmv
##4# For every node in the network, 
##5# Delete it from the network
##6# Compute the cluster/communities in the new graph (graph.rmv)
##7# Obtain the member-proteins of each of the (new) clusters (cluster.rmv's) in the perturbed graph
#
#   At this point, it was necessary to identify the cluster.orig in focus from all the produced cluster.rmv's; this was done by finding which the cluster.rmv's 
#	has the highest numbers of common/similar members with cluster.orig. It is this that was taken to be the 'transformed cluster.orig', that is, the corresponding cluster.rmv
#	When the corresponding cluster.rmv is obtained, the jaccard coefficient is calculated (between cluster.orig and its corresponding cluster.rmv)
	Steps #8 to #10 achieve this
##8# Initiate the vector to store the number of common/similar members between cluster.orig and each cluster.rmv
##9# For each cluster.rmv,
##10# Initiate count of similar proteins
##11# For each member of cluster.rmv, compare it to each member of cluster.orig and if they are the same, increase the count of similar proteins between 
#		cluster.orig and this cluster.rmv by 1
##12# Store the number of similar proteins for this cluster.rmv in a vector
#
#	REPEAT STEPS 10 to 12 FOR ALL CLUSTER.RMV's
#
##13# Since the similarity counts for each cluster.rmv are stored in a vector, i need to be able to identify the identity of the cluster.rmv with the maximum count; for that reason,
#		i attach a 'names' attribute to each of the count. The attribute is the position of the cluster.rmv in the vector 
##14# To get the identity of the cluster.rmv with the highest number of similar proteins as cluster.orig, that is, the corresponding cluster.rmv, i get the attribute(forced to be a number)
#		of highest number in vector of similarity counts
##15# Then i calculate the Jaccard coefficient between cluster.orig and cluter.rmv using the formula; Jaccard=members in both/(members in A_not_in_B + members_in_B_not_A + members in both))
##16# This jaccard score is stored in the list of scores created just before step 3
#
#	REPEAT STEPS #5 to #16 999 times
#
##17# Store the 1000 jaccard scores as a list in the list 'all_similarity_scores'
#
##	REPEAT STEPS #4 to #17 for every other original cluster
#
##18#	Return the list of lists that is 'all_similarity_scores'
#
#The list of jaccard coefficients for, e.g cluster one, can be assessed like this
#scores_list.ML<- node_removal.jac.ML(graph)
#scores_list.ML[[1]]
#
#

#######
# THE FUNCTION
###

node_removal.jac.ML<- function(graph){
#1#
	communities_ML.orig <- multilevel.community(graph, weights=E(graph)$weight)
#2#
	clusters.orig<- show_clusters(communities_ML.orig, "Multi-level")
	all_similarity_scores<- NULL
#3#
	for (n in 1:length(clusters.orig)) {
		jaccard_scores<- NULL
#4#
		for (i in 1:vcount(graph)) {
#5#
			graph.rmv <- delete.vertices(graph, V(graph)[i])
#6#
			communities_ML.rmv <- multilevel.community(graph.rmv, weights=E(graph.rmv)$weight)
#7#
			clusters.rmv<- show_clusters(communities_ML.rmv, "multi-level.rmv")
#8#
			simlty_counts<- NULL
#9#
			for (j in 1:length(clusters.rmv)){
#10#
				sim_count=0
#11#
				for (k in 1:length(clusters.rmv[[j]]))   {
					for (m in 1:length(clusters.orig[[n]]))   {
						if (clusters.orig[[n]][m] == clusters.rmv[[j]][k]){
							sim_count=sim_count+1		
						}
					}
				}
#12#
				simlty_counts<- c(simlty_counts, sim_count)
			}
#13#
			attr(simlty_counts, "names")<- c(1:length(clusters.rmv))
#14#
			similar_cluster<- as.numeric(attributes(simlty_counts[match(c(max(sapply(simlty_counts, max))), simlty_counts)])[[1]])
#15#
			jaccard<- max(sapply(simlty_counts, max))/(   length(clusters.orig[[n]]) + length(clusters.rmv[[similar_cluster]]) + max(sapply(simlty_counts, max))   )
#16#
			jaccard_scores<- c(jaccard_scores, jaccard)
		}
#17#
		all_similarity_scores <- c(all_similarity_scores, list(jaccard_scores))
	}
#18#
	return(all_similarity_scores)
}

#***********SCR_8******************stop








#***********SCR_9******************start

###########################################################################################################################################
			###########################################################################################################
#
#									EDGE-WEIGHT FILTRATION
#
			###########################################################################################################
#########################################################################################################################################

# The communities were computed and the clusters extracted; the changes in the various clusters were then depicted with figure 5 and 6 in the text, for the Fastgreedy and
#	Multi-level algorithms respctively.
#



ggg=delete.edges(g, which(E(g)$weight <=0.1))		#Filtering out all edges with weights less than or equal to 0.1
ggg2=delete.edges(g, which(E(g)$weight <=0.2))		#Filtering out all edges with weights less than or equal to 0.2
ggg3=delete.edges(g, which(E(g)$weight <=0.3))		#Filtering out all edges with weights less than or equal to 0.3



#####
# FAST-GREEDY 
#####

communities_FGpt1 <- fastgreedy.community(ggg, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(ggg)$weight)	#compute clusters in graph with only >0.1 weights
communities_FGpt2 <- fastgreedy.community(ggg2, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(ggg2)$weight)	#compute clusters in graph with only >0.2 weights
communities_FGpt3 <- fastgreedy.community(ggg3, merges=TRUE, modularity=TRUE, membership=TRUE, weights=E(ggg3)$weight)	#compute clusters in graph with only >0.3 weights

show_clusters(communities_FGpt1, "fastgreedy")		#Show the identity of proteins in 'communities_FGpt1'
show_clusters(communities_FGpt2, "fastgreedy")		#Show the identity of proteins in 'communities_FGpt2'
show_clusters(communities_FGpt3, "fastgreedy")		#Show the identity of proteins in 'communities_FGpt3'


######
#MULTI-LEVEL
######

communities_MLpt1 <- multilevel.community(ggg, weights=E(ggg)$weight)	#compute clusters in graph with only >0.1 weights
communities_MLpt2 <- multilevel.community(ggg2, weights=E(ggg2)$weight)	#compute clusters in graph with only >0.2 weights
communities_MLpt3 <- multilevel.community(ggg3, weights=E(ggg3)$weight)	#compute clusters in graph with only >0.3 weights

show_clusters(communities_MLpt1, "multilevel")	#Show the identity of proteins in 'communities_MLpt1'
show_clusters(communities_MLpt2, "multilevel")	#Show the identity of proteins in 'communities_MLpt2'
show_clusters(communities_MLpt3, "multilevel")	#Show the identity of proteins in 'communities_MLpt3'


#***********SCR_9******************stop
