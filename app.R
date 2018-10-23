# -----------------------------
# Author: Almie Carajay
# Date Started: June 30, 2017
# -----------------------------
library(shiny)
library(ape)
library(seqinr)
library(hash)
library(plyr)
library(stringr)
library(tidyr)
library(shinycssloaders)
#library(shinyBS)
#library(shinyjs)
#library(V8) 
#library(dplyr)
#library(boot)  


shinyApp(
  

    ui = shinyUI(htmlTemplate("www/index.html")),
  
    server = shinyServer(function(input, output) { 

    	
		radiodata <- reactive({
			radioval <- input$radioval
			if(is.null(radioval)){return(counter=0)}
			radioval
		}) 

		validSeq <- reactiveValues()
		validSeq$invalid = FALSE
		validSeq$ret = 0


	    output$sumofdata <- renderDataTable({
	    	radioval <- input$radioval
	    	validSeq$beforetime <- Sys.time()  
	        #btime <- system.time()
	        #print(btime)	
	        print(c("Before time:", validSeq$beforetime))
	        #print(validSeq$beforetime)

	    	if(is.null(radioval)){return()}
	    	gene.Name <- ''

 
	    	summarydata <- function(seqlist, numofcodonList, acc){  # prints table of sequences and the number of codons analyzed
	    		j = 1
	    		annotList=''
	    		nrow=''
	    		counter = 1  	

		        for(j in 1:length(seqlist)){
		        	nrow <- append(nrow, j)
		        	if(acc==1){
						annotList[j] <- c(gene.Name[j])
		        	}else{
		           		annotList[j] <- c(getAnnot(seqlist[[j]]))
		        	}
		            
		        }
		        frow <- nrow[2:length(nrow)]  
		        #print(c("frow:",length(frow)))
		        #print(c("annotlist:", annotList))
		        #print(c("numofcodonList:", length(numofcodonList)))
		        validSeq$numcol<-sum(as.numeric(numofcodonList), na.rm = TRUE)
		        print(validSeq$numcol)
		        #print(c("Invalid value:"), validSeq$invalid)

		        if(length(numofcodonList) != length(frow)){
		        	#print("INVALID INPUT!!!! contains a non-ATGC value.")
		        	HTML('<h3> INVALID INPUT!!!!<br /> Input contains a non-ATGC value </h3>')
		        	validSeq$ret = 1
		        	print(c("ret:",validSeq$ret))
		        	return(validSeq$ret)
		        }else{

		    		sumtable <- data.frame(frow, annotList,numofcodonList)  # table of codons with frequency per codon
		    		#print(c("Sum: ",colSums(as.numeric(sumtable), na.rm = TRUE)))
		    		colnames(sumtable) <- c("#","Sequence(s)", "Number of Codons Analyzed")
					write.table(sumtable, file="summary_of_input_data.csv", sep=",", row.names=FALSE)
					validSeq$sumtable <- sumtable
					validSeq$aftertime <- Sys.time() 
					validSeq$timediff <- Sys.time() - validSeq$beforetime
					#stime <- system.time()
					print(c("After time:", validSeq$aftertime))
					#print(validSeq$aftertime)
					#print(validSeq$timediff)
					validSeq$ret = 2
					sumtable


				}
			} 


	    	codonrowcount <- function(seqlist, acc){
	    		#for each sequence in seqlist extract each codon
	        	a = 1
	        	#print(length(seqlist))
	        	dnaletters=''
	        	numofcodonList=''
	        	allcodons=''
	        	allaminoacids=''
	        	allidentifier=''
	        	allaminonames=''
	        	preferedindex=''
	        	preferedcodon=''

	        	#print(c("codon ret", validSeq$ret))
	        	#print(acc)
	        	#print(attr(seqlist,"description"))
	        	#print(speciesAccessno[c])
	        	#print(seqlist)
	        	#print(acc)
	        	if(is.null(seqlist) || seqlist=="" || is.na(seqlist)==TRUE){return()}

	        	withProgress(message = 'Working on codon', value = 0, style = getShinyOption("progress.style", default = "notification"),{
		        # Number of times we'll go through the loop
		      	n1 <- length(seqlist)

	        	for(a in 1:n1){   #index of each sequence in seqlist
	        		#print(seqlist[a])  # but sequences here is of string value
	          		#extract each character in a sequence of seqlist
	          		invalid=FALSE
					#print(seqlist[[a]])
	        		#print(acc)
	          		
	          
	          		if(acc==1){
	          			seqstring <- trimSpace(seqlist[[a]])   # trims (" Some text. " == "Some text.")
	          		}else{
	          			seqstring <- trimSpace(seqlist[a])   # trims (" Some text. " == "Some text.")
	          		}
	          		newseq <- str_replace_all(seqstring, "[ \r\n]" , "N")   #replace all newline, tabs and whitespace with ""
	          		bases <- unlist(str_split(newseq, ""))  #[[1]] splits into individual letters/characters
	          		#bases <- unlist(str_split(seqstring, ""))  #[[1]] splits into individual letters/characters

	          		dnaletters <- toupper(bases)
	          		#print(newseq)
	          		#print(seqstring)
	          		#extract each codon for each sequence in seqlist
			        #print(length(dnaletters))
			        b = 1
			        codon=''
			        aminoacid=''
			        id_amino=''
			        amino_name=''
			        codonlist <- hash("ATA"="I","ATC"="I","ATT"="I","ATG"="M",
			                            "ACA"="T","ACC"="T","ACG"="T","ACT"="T",
			                            "AAC"="N","AAT"="N","AAA"="K","AAG"="K",
			                            "AGC"="S","AGT"="S","AGA"="R","AGG"="R",
			                            "CTA"="L","CTC"="L","CTG"="L","CTT"="L",
			                            "CCA"="P","CCC"="P","CCG"="P","CCT"="P",
			                            "CAC"="H","CAT"="H","CAA"="Q","CAG"="Q",
			                            "CGA"="R","CGC"="R","CGG"="R","CGT"="R",
			                            "GTA"="V","GTC"="V","GTG"="V","GTT"="V",
			                            "GCA"="A","GCC"="A","GCG"="A","GCT"="A",
			                            "GAC"="D","GAT"="D","GAA"="E","GAG"="E",
			                            "GGA"="G","GGC"="G","GGG"="G","GGT"="G",
			                            "TCA"="S","TCC"="S","TCG"="S","TCT"="S",
			                            "TTC"="F","TTT"="F","TTA"="L","TTG"="L",
			                            "TAC"="Y","TAT"="Y","TAA"="_","TAG"="_",
			                            "TGC"="C","TGT"="C","TGA"="_","TGG"="W")

			        identifier <- hash("ATA"="Ile1","ATC"="Ile2","ATT"="Ile3","ATG"="Met",
	                          "ACA"="Thr1","ACC"="Thr2","ACG"="Thr3","ACT"="Thr4",
	                          "AAC"="Asn1","AAT"="Asn2","AAA"="Lys1","AAG"="Lys2",
	                          "AGC"="Ser1","AGT"="Ser2","AGA"="Arg1","AGG"="Arg2",
	                          "CTA"="Leu1","CTC"="Leu2","CTG"="Leu3","CTT"="Leu4",
	                          "CCA"="Pro1","CCC"="Pro2","CCG"="Pro3","CCT"="Pro4",
	                          "CAC"="His1","CAT"="His2 ","CAA"="Gln1","CAG"="Gln2",
	                          "CGA"="Arg3","CGC"="Arg4","CGG"="Arg5","CGT"="Arg6",
	                          "GTA"="Val1","GTC"="Val2","GTG"="Val3","GTT"="Val4",
	                          "GCA"="Ala1","GCC"="Ala2","GCG"="Ala3","GCT"="Ala4",
	                          "GAC"="Asp1","GAT"="Asp2","GAA"="Glu1","GAG"="Glu2",
	                          "GGA"="Gly1","GGC"="Gly2","GGG"="Gly3","GGT"="Gly4",
	                          "TCA"="Ser3","TCC"="Ser4","TCG"="Ser5","TCT"="Ser6",
	                          "TTC"="Phe1","TTT"="Phe2","TTA"="Leu5","TTG"="Leu6",
	                          "TAC"="Tyr1","TAT"="Tyr2","TAA"="STOP","TAG"="STOP",
	                          "TGC"="Cys1","TGT"="Cys2","TGA"="STOP","TGG"="Trp")

			        aminoNames <- hash("ATA"="Isoleucine","ATC"="Isoleucine","ATT"="Isoleucine","ATG"="Methionine",
	                          "ACA"="Threonine","ACC"="Threonine","ACG"="Threonine","ACT"="Threonine",
	                          "AAC"="Asparagine","AAT"="Asparagine","AAA"="Lysine","AAG"="Lysine",
	                          "AGC"="Serine","AGT"="Serine","AGA"="Arginine","AGG"="Arginine",
	                          "CTA"="Leucine","CTC"="Leucine","CTG"="Leucine","CTT"="Leucine",
	                          "CCA"="Proline","CCC"="Proline","CCG"="Proline","CCT"="Proline",
	                          "CAC"="Histidine","CAT"="Histidine","CAA"="Glutamine","CAG"="Glutamine",
	                          "CGA"="Arginine","CGC"="Arginine","CGG"="Arginine","CGT"="Arginine",
	                          "GTA"="Valine","GTC"="Valine","GTG"="Valine","GTT"="Valine",
	                          "GCA"="Alanine","GCC"="Alanine","GCG"="Alanine","GCT"="Alanine",
	                          "GAC"="Aspartic acid","GAT"="Aspartic acid","GAA"="Glutamic acid","GAG"="Glutamic acid",
	                          "GGA"="Glycine","GGC"="Glycine","GGG"="Glycine","GGT"="Glycine",
	                          "TCA"="Serine","TCC"="Serine","TCG"="Serine","TCT"="Serine",
	                          "TTC"="Phenylalanine","TTT"="Phenylalanine","TTA"="Leucine","TTG"="Leucine",
	                          "TAC"="Tyrosine","TAT"="Tyrosine","TAA"="STOP CODON","TAG"="STOP CODON",
	                          "TGC"="Cysteine","TGT"="Cysteine","TGA"="STOP CODON","TGG"="Tryptophan")

					clist <- c("ATA","ATC","ATT","ACA","ACC","ACG","ACT","AAC","AAT","AAA","AAG",
			            "AGC","AGT","AGA","AGG","CTA","CTC","CTG","CTT","CCA","CCC","CCG",
			            "CCT","CAC","CAT","CAA","CAG","CGA","CGC","CGG","CGT","GTA","GTC",
			            "GTG","GTT","GCA","GCC","GCG","GCT","GAC","GAT","GAA","GAG","GGA",
			            "GGC","GGG","GGT","TCA","TCC","TCG","TCT","TTC","TTT","TTA","TTG", 
			            "TAC","TAT","TGC","TGT")

					for(b in 1:length(dnaletters)){
	          			#print(dnaletters[b])    # "A","T","G","C"
	          			if(dnaletters[b]=="A" || dnaletters[b]=="T" || dnaletters[b]=="G" || dnaletters[b]=="C"){
	          				if(b%%3==0){   # every 3rd character
		              			nucleotide <- c(dnaletters[(b-2):b])    #"A""T""G"
		              			y <- paste(nucleotide, collapse="")  #"ATG"
		              			codon <- append(codon, y)
		              			aminoacid <- append(aminoacid, codonlist[[y]])              #prints list of protein ("M","M""V","M""V""S")
		                		id_amino <- append(id_amino, identifier[[y]])      # prints list of identifier
		                		amino_name <- append(amino_name, aminoNames[[y]])   # prints list names of Amino acid
		                	}
		               	}else if(dnaletters[b]=="\n"){
		               		print("there is newline")

		               	}else{
		               		validSeq$invalid = TRUE
	          				print(c("Sequence contains invalid input:", dnaletters[b]))
	                   		validate(
	            		   		need(invalid==TRUE, paste('Sequence contains invalid input:', dnaletters[b]))
	            			
					   		)
	          				
	          				#validSeq$errorletter <- dnaletters[b]
	          				validSeq$ret = 1
	          				print(validSeq$invalid)
	          				break
	          			}	            
	          		}
	          

	            	if(invalid==TRUE){
	            		#print(c("INVALID SEQUENCE OIE",validSeq$invalid))
	            		#validSeq$invalid = FALSE
	            		#HTML()
	            		validate(
	            		   	need(invalid==TRUE, 'Input sequence contains a non [A,T,G,C] base. Please input valid sequence(s)!!!')
	            			
					    )
					    
					    #isolate({validSeq$invalid = FALSE})
	           			#print(validSeq$invalid)
	           			#break
	           			return()

	          		}else{        
		        	#solve for the codon count 
		          
		        		proteinseq <- paste(aminoacid, collapse="")
		        		if(codon[1]=="" || protein[1]==""){
		        			codon <- codon[2:length(codon)]
		        			aminoacid <- aminoacid[2:length(aminoacid)]
		        			id_amino <- id_amino[2:length(id_amino)]
		        			amino_name <- amino_name[2:length(amino_name)]
		        		}

				        tb <- data.frame( codon, aminoacid)  # table of codons with frequency per codon

				        #View(tb)
				        numofcodon <- length(codon)
				        codoncount <- count(tb, c('codon','aminoacid'))   #counts for the frequency
		                numofcodonList <- append(numofcodonList, numofcodon)
			        }
	            	allcodons <- append(allcodons,codon)  #appends each codons per sequence in a given input
	          		allaminoacids <- append(allaminoacids,aminoacid)
	          		allidentifier <- append(allidentifier, id_amino)
	          		allaminonames <- append(allaminonames, amino_name)

	          		# Increment the progress bar, and update the detail text.
			        incProgress(1/n1, detail = paste(": Analyzing ", n1, "sequence(s)"))
			        #print(c("n1:",n1))
			        # Pause for 0.1 seconds to simulate a long computation.
			        Sys.sleep(0.1)   
			        #close()

		        }

		        })

		        CODONS <- allcodons[2:length(allcodons)]
		        AMINO_ACID <- allaminoacids[2:length(allaminoacids)]
		        IDENTIFIER <- allidentifier[2:length(allidentifier)]
		        NAME <- allaminonames[2:length(allaminonames)]
		        #print(CODONS)   #prints all codons in a given input except null values
		        #print(AMINO_ACID)

		        #print("CODONLIST:")
		        #print(codonlist[[]])
		        slist <- sort(clist)      

		        q1 = 1
		        q2 = 1
		        amlist = ''
		        idlist = ''
		        alist = ''
		        total = ''
		        for(q1 in 1:length(slist)){
		        	amlist <- append(amlist, codonlist[[slist[q1]]])
		        	idlist <- append(idlist, identifier[[slist[q1]]])
		        	alist <- append(alist, aminoNames[[slist[q1]]])

		        }

		        #print("total:")
		        #print(total)

		        fstable <- data.frame(slist)
		        fstable$aminolist <- amlist[2:length(amlist)]
		        fstable$identlist <- idlist[2:length(idlist)]
		        fstable$anames <- alist[2:length(alist)]
		        #fstable$rawcount <- count(AMINO_ACID)	        
		        #View(fstable)

		        cut <- data.frame(CODONS, IDENTIFIER, AMINO_ACID, NAME)
		        newcut <- count(cut,c('CODONS', 'IDENTIFIER', 'AMINO_ACID','NAME'))
		        #View(newcut)

		        totalcount <- count(AMINO_ACID)
		        #View(totalcount)		               


		        n = 1
		        m = 1
		        TOTALCOUNT = ''
		        for(n in 1:length(newcut[,3])){
		        	#print(newcut[n,3])
		        	for(m in 1:length(totalcount[,1])){
		        		#print(totalcount[m,1])
		        		if(newcut[n,3]==totalcount[m,1]){
		        			TOTALCOUNT <- append(TOTALCOUNT, totalcount[m,2])
		        		}else{
		        			#print("error xa dai")
		        			
		        		}
		        	}
		        }
		        newcut$TOTAL <- TOTALCOUNT[2:length(TOTALCOUNT)]
		        newcut <- transform(newcut, X = freq/as.numeric(TOTAL))
		        #print(newcut$freq)
		        #print(newcut$TOTAL)
		        #print(data.matrix(newcut$TOTAL))
				#View(newcut)

		        #Maximum Likelihood estimation of parameters of the MLE based on the observed data
				a1 = 1

				mlist = dim(0)
				selist = dim(0)
				mle<-function(x) {

				  #withProgress(message = 'Calculating MLE...', value = 0, {
				  	  n1 <- 200
					  for(a1 in 1:n1){
					  	experiments <- rmultinom(n=1, size=1000, prob=x)
					  	df=data.frame(experiments)/1000
					  	MLE<-df[,1]
					  	SE<-sqrt((MLE*(1-MLE))/1000)
					  	#codontemp <- cbind(MLE,SE)
					  	mlist <- data.matrix(cbind(mlist, MLE))
					  	selist <- data.matrix(cbind(selist, SE))

					
					  }
					    #print(is.numeric(mlist))
					    #print(is.numeric(selist))
					  	meanml<-rowMeans(mlist, na.rm = TRUE)
					  	meanse<-rowMeans(selist, na.rm = TRUE)
					    temptable <- data.frame(meanml, meanse)
					  	#View(temptable)
					  	#MLE1 = 1

				  return(cbind(meanml,meanse))
				}
				mtable <- data.frame(mle(newcut$X))
				codontable <- cbind(newcut, mtable)
				#View(mtable)
				#View(codontable)

				f = 1
				g = 1
				
				LSElist=''
				
				for(f in 1:length(totalcount[,1])){
					count = 0
					SElist =''
					indexlist=''
					for(g in 1:length(codontable[,3])){						
						if(totalcount[f,1]==codontable[g,3]){
							count = count + 1
							#print(c("Count:",count))
							SElist <- append(SElist, codontable[g,9])
							indexlist <- append(indexlist, g)
						
						} 

					}
					#print("====================================")
					#print(count)

					newSElist <- SElist[2:length(SElist)]
					newindexlist <- indexlist[2:length(indexlist)]
					#print(newSElist)
					#print(newindexlist) 
					#print(totalcount[f,1])
					lse <- min(SElist[2:length(SElist)])

					h = 1					
					for(h in 1:length(newindexlist)){
						
						if(lse == newSElist[h]){
							e <- newindexlist[h]
							preferedindex <- append(preferedindex,e)
							
						}
					}
					LSElist <- append(LSElist, lse)
					
				}
				

				index <- as.numeric(preferedindex[2:length(preferedindex)])
				sortIndex <- sort(index, decreasing = FALSE)

				p = 1
				q = 1
				pcodon <- ''
				
				for(p in 1:length(codontable[,1])){
					
					for(q in 1:length(sortIndex)){
						if(p == sortIndex[q]){
							preferedcodon <- append(preferedcodon, as.character(codontable[p,1]))
							
						}
					}

				}
				
				P_CODON <- preferedcodon[2:length(preferedcodon)]
				#print(P_CODON)
				xy = 1
				yy = 1
				z = 1 
				#print(length(sortIndex))
				#print(sortIndex[21])
				#print(length(codontable[,1]))

				for(z in 1:length(codontable[,1])){
					xx = 1

					for(xy in yy:length(codontable[,1])){
						
						if(is.na(sortIndex[xy]) == TRUE){
							xx = 0
							break

						}else if(z == sortIndex[xy]){
							pcodon <- append(pcodon, P_CODON[xy])
							yy = xy + 1
							break
						}else{							
							xx = 0
							break

						}
					}		

					if(xx == 0){
						pcodon <- append(pcodon, "-")
					}
					
				}
				
				PREFERREDCODON <- pcodon[2:length(pcodon)]
				codontable$PREFERRED_CODON <- PREFERREDCODON


				PREFERRED_CODON <- LSElist[2:length(LSElist)]
				#print(LSElist[2:length(LSElist)])
				AMINO <- totalcount[,1]
				LSEtable <- data.frame(AMINO, PREFERRED_CODON)
				#View(LSEtable)

				#View(codontable)
				write.table(codontable, file="codon_usage_table.csv", sep=",", row.names=FALSE)
				validSeq$codontable <- codontable

				aminocode <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "N", "P", "Q", "R", "S", "T", "V", "Y"  )
				myaminoname <- c("Alanine","Cysteine","Aspartic acid","Glutamic acid", "Phenylalanine", "Glycine", "Histidine","Isoleucine", "Lysine", "Leucine", "Asparagine", "Proline", "Glutamine", "Arginine", "Serine", "Threonine", "Valine", "Tyrosine")
				codongen <- c("GCT, GCC, GCA, GCG", "TGT, TGC", "GAT, GAC", "GAA, GAG", "TTT, TTC", "GGT, GGC, GGA, GGG", "CAT, CAC", "ATT, ATC, ATA", "AAA, AAG", "CTT, CTC, CTA, CTG, TTA, TTG", "AAT, AAC", "CCT, CCC, CCA, CCG", "CAA, CAG", "CGT, CGC, CGA, CGG, AGA, AGG", "TCT, TCC, TCA, TCG, AGT, AGC", "ACT, ACC, ACA, ACG", "GTT, GTC, GTA, GTG", "TAT, TAC")
				pretable <- data.frame(aminocode, myaminoname, codongen)
				#View(pretable)
				#print(P_CODON)
				#print(as.character(codontable$PREFERRED_CODON))
				#print(index)

				r = 1
				am_acid = ''
				hashlist = ''

				for(r in 1:length(P_CODON)){
					am_acid <- append(am_acid, codonlist[[P_CODON[r]]])
					#hashcodon <- hash(codonlist[[P_CODON[r]]], P_CODON[r])
					#print(hashcodon)
					#hashlist <- append(hashlist,hashcodon)
				}



				#print("=========sortedpcodon============")
				newamino <- sort(am_acid)
				sortedpcodon <- sort(P_CODON)
				#print(newamino)
				#print(aminocode)
				#print(sortedpcodon)
				#print(sortedpcodon)
				#print("==========P_CODON========")
				#print(P_CODON)
				s = 1
				t = 1
				pcodonlist=''
				val = 0
				for(s in 1:length(aminocode)){
					#print(c(aminocode[s]," = "))
					#print(codonlist[[P_CODON[s]]])
					for(t in 1:length(P_CODON)){
						if(aminocode[s]==codonlist[[P_CODON[t]]]){
							#print(c(aminocode[s],":"))
							#print(P_CODON[t])
							hashcodon <- hash(aminocode[s], P_CODON[t])
							#print(hashcodon)
							pcodonlist <- append(pcodonlist, P_CODON[t])
							val = 1
							break
						}else{
							val = 0
						}
					}
					#print(val)
					if(val == 0){
						#print(c(aminocode[s],":", "NONE"))
						pcodonlist <- append(pcodonlist, "-")
					}

				}
				#u = 1
				pclist <- pcodonlist[2:length(pcodonlist)]
 
				#print(pclist)
				pretable$PREFERRED_CODON <- pclist

				#View(pretable)
				validSeq$pretable <- pretable
				#---------------- Write into table the preferred codon ----------------------
				write.table(pretable, file="preferredcodontable.csv", sep=",", row.names=FALSE)
		 		
		        st <- summarydata(seqlist, numofcodonList[2:length(numofcodonList)], acc)
	    	}

	

	    	if(radioval=="one"){                 #enter sequence(s)
	    		genome <- input$sequence
	    		#----------------- getting sequence(s) from user input ------------------
	    		if(is.null(genome) || genome==""){return()}else{
			        #write input sequence into genome.fasta
			        writeGenome <- write.fasta(sequences = genome, names = names(genome),  file.out = "genome.fasta")
			        #read the input sequence from a genome.fasta, sequence as DNA, string, uppercase, sequence attributes should be set, only sequences as returned , the '>' at the beginning of the description lines is removed in the annotations of the sequences
			        genefile <- read.fasta(file = "genome.fasta", seqtype = c("DNA"), as.string = TRUE, forceDNAtolower = FALSE, set.attributes = TRUE, bfa= FALSE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = TRUE)
	    		}
	    		#------------------------------------------------------------------------
	    		#create a new list containing only the sequences
		        seqlist=''
		        annotList=''
		        i = 2
		        j = 1
		        for(i in 2:length(genefile)){
		          seqlist[j] <- c(genefile[i])
		          annotList[j] <- c(getAnnot(genefile[[i]]))
		          j = j+1
		        }
		        acc <- 0
		        #print(seqlist)
		        nc <- codonrowcount(seqlist, acc)

	    	}else if(radioval=="two"){
	    		access.nb <- input$accessno

	    		#----------------- getting sequence(s) from user input ------------------
	    		if(is.null(access.nb) || access.nb==""){return()}else{
			      	genefile <- read.GenBank(access.nb)
			        seqlist <- read.GenBank(access.nb, species.names = TRUE,
	                   gene.names = FALSE, as.character = TRUE)
			      #print(genefile)
                  #speciesName <- attr(genefile, "species") ## get the species names
			      #speciesAccessno <- names(genefile)
                  ## build a matrix with the species names and the accession numbers:
			      #gen <- cbind(attr(genefile, "species"), names(genefile))		      
			      gene.Name <- attr(genefile, "description")  ## the description from each FASTA sequence:
			      #print(gen)
	    		}
	      		#-----------------------------------------------------------------------
	      		acc <- 1
	      		st <- codonrowcount(seqlist, acc)

	    	}else if(radioval=="three"){
			    # input$file1 will be NULL initially. After the user selects
			    # and uploads a file, head of that data file by default,
			    # or all rows if selected, will be shown.
			    inputfile <- req(input$inputfile)
			    #print(inputfile)
			    filesize <- input$inputfile$size
			    #print(filesize)

	    		if(is.null(input$inputfile$datapath)){
	    			#print(input$inputfile$datapath)
	    			return()
	    		}else{
			      # read the content of input file
			      genefile <- read.fasta(file = input$inputfile$datapath, seqtype = c("DNA"), as.string = TRUE, forceDNAtolower = FALSE, set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = TRUE)
			      #print(genefile)
			      #print(length(dnafile))
			      acc <- 0
			      nc <- codonrowcount(genefile, acc)
	      		}

	    	}else{
	    		print(radioval)
	    		print("No sequence input!")
	    	}

	    })

		#===============================================================================


		output$precodontable <- renderDataTable({
			#print(c("precodontable ret value:", validSeq$ret))
			if(radiodata()!= 0){
				if(validSeq$invalid == FALSE && validSeq$ret==2){
					pdata <- validSeq$pretable
					psdata <- data.frame(pdata$myaminoname,  pdata$aminocode,  pdata$codongen, pdata$PREFERRED_CODON)

					#print(ncol(psdata))
					if(ncol(psdata)!= 4){
						#print("NONE!")
						validSeq$ret = 1
						return(validSeq$ret)
					}else{
						colnames(psdata) <- c( "Amino Acid","One-Letter Code", "Codons that Encode an Amino Acid", "Preferred Codon")
						write.table(psdata, file="preferredcodontable.csv", sep=",", row.names=FALSE)
						validSeq$psdata <- psdata

					}

				}
			}else{
				#print("No Input Yet!")
			}


		})	

		#===============================================================================

		output$codonusagetable <- renderDataTable({
			
			if(radiodata()!=0){
				if(validSeq$invalid == FALSE){
					sdata <- validSeq$codontable	
					ptable <- data.frame( sdata$NAME, sdata$AMINO_ACID, sdata$CODONS, sdata$IDENTIFIER,  sdata$freq, sdata$meanml, sdata$meanse, sdata$PREFERRED_CODON)  # table of codons with frequency per codon

					if(ncol(ptable)== 0){
						#print("NONE")
						validSeq$ret = 1
						return(validSeq$ret)
					}else{
			    		colnames(ptable) <- c("Amino Acid","One-Letter Code","Codon","Identifier","Raw Count","Maximum Likelihood Estimator", "Standard Error", "Preferred Codon")
						write.table(ptable, file="codon_usage_table.csv", sep=",", row.names=FALSE)
						
						validSeq$ptable <- ptable
					}

				}
			}else{
				#print("No Input Yet!")
			}

		})	

		output$code1 <- renderUI({ 
			#print(c("ret:", validSeq$ret))  
			if(radiodata()!=0){
				if(validSeq$invalid == FALSE && validSeq$ret==2){
					HTML('<h3>SUMMARY OF INPUT DATA</h3>')
				}else{
					#HTML('<h3>INVALID INPUT!</h3>')
				}   

   				
   			} else {
   				#HTML('<h3> NO INPUT YET </h3>
				#		<a href="#download" >
		      	#    	<button  onclick="javascript:document.getElementById(\'picker\').reset();" class="btn btn-default"> BACK TO INPUT</button> </a>
				#		<br />
				#	<br /> <br /> <br /> <br /> 
   				#')
   			}
  		})

		output$code2 <- renderUI({ 
			if(radiodata()!=0){
				if(validSeq$invalid == FALSE && validSeq$ret==2){
					HTML('<h3>SUMMARY OF PREFERRED CODON</h3>')
				}   				
   				
   			}
   			# else {
   			#	HTML('<h3> NO INPUT YET </h3>
			#			<a href="#download" >
		    #  	    	<button  onclick="javascript:document.getElementById(\'picker\').reset();" class="btn btn-default"> BACK TO INPUT </button> </a>
			#			<br />
			#		<br /> <br /> <br /> <br /> <br /> <br /> <br /> <br /> <br /> <br /> <br /> <br />
   			#	')
   			#}
  		})


		output$code3 <- renderUI({ 

			if(radiodata()!=0){
				if(validSeq$invalid == FALSE && validSeq$ret == 2){
   					HTML('<h3> <br /> <br/> Codon Usage Table </h3>')
   				}
   			} 
   			#else {
   			#	HTML('<h3> NO INPUT YET </h3>
			#			<a href="#download" >
		    #  	    	<button  onclick="javascript:document.getElementById(\'picker\').reset();" class="btn btn-default"> BACK TO INPUT </button> </a>
			#			<br />
			#		<br /> <br /> <br /> <br /> <br />  <br /> <br /> <br /> <br /> <br /> <br /> <br />
   			#	')
   			#}
  		})

  		output$code4 <- renderUI({  
			if(radiodata()!=0){	
				if(validSeq$invalid == FALSE && validSeq$ret == 2){	
	      	   		downloadButton("downloadSData", label="Download Table", class="btn-default") 
	      	   	}
   			} 
  		})

  		output$code5 <- renderUI({    
 			
			if(radiodata()!=0){		
				if(validSeq$invalid == FALSE && validSeq$ret == 2){
		      	   HTML('&nbsp; <a href="#output" class="btn js-scroll-trigger" >
		            <button class="btn btn-default">NEXT</button>
		        	</a>
		        	')
	      		}else{

	      			#tags$a(href="javascript:history.go(0)", 
           			#popify(tags$i(class="fa fa-refresh fa-5x"),
                 	#title = "Reload", 
                  	#content = "Click here to restart the Shiny session",
                  	#placement = "center"))
		      	    HTML('&nbsp;<a href="#download" >
		      	    	<br />
		      	    	<br /> 
		      	    	<button type="reset" id="form_reset" onclick="javascript:history.go(0)" class="btn btn-default"> BACK TO INPUT </button> </a>
				    	<br />
	   				') 
   					
	      		}
   			} 
  		})

  		output$code6 <- renderUI({   
			if(radiodata()!=0){	
				if(validSeq$invalid == FALSE && validSeq$ret == 2){	
	      	   		downloadButton("downloadPtable", label="Download Table", class="btn-default")
	      		}
   			} 
  		})

  		output$code7 <- renderUI({   
			if(radiodata()!=0){	
				if(validSeq$invalid == FALSE && validSeq$ret == 2){	
		      	   HTML('&nbsp; <a href="#ctable" class="btn js-scroll-trigger" >
		            <button class="btn btn-default">NEXT</button>
		        	</a>
		        	')
	      		}
   			} 
  		})

  		output$code8 <- renderUI({   
			if(radiodata()!=0){

				if(validSeq$invalid == FALSE && validSeq$ret == 2){		
	      	   		downloadButton("downloadData", label="Download Table", class="btn btn-default")
	      		}
   			} 
  		})	 

  		output$code9 <- renderUI({   
			if(radiodata()!=0){	
				if(validSeq$invalid == FALSE && validSeq$ret == 2){	
		      	   HTML('&nbsp; <a href="#cplot" class="btn js-scroll-trigger" >
		            <button class="btn btn-default">NEXT</button>
		        	</a>
		        	')
	      		}
   			} 
  		})

  		output$code10 <- renderUI({
  		#onclick="javascript:document.getElementById(\'picker\').reset();"   
			if(radiodata()!=0){	
				if(validSeq$invalid == FALSE && validSeq$ret == 2){

					#tags$a(href="javascript:history.go(0)", 
           			#popify(tags$i(class="fa fa-refresh fa-5x"),
                  #title = "Reload", 
                  #content = "Click here to restart the Shiny session",
                  #placement = "right"))	
		      	    HTML('&nbsp;<a href="#download" >
		      	    	<br />
		      	    	<br /> 
		      	    	<button type="reset" id="form_reset" onclick="javascript:history.go(0)" class="btn btn-default"> BACK TO INPUT </button> </a>
				    	<br />
	   				') 

	      		}
   			} 
  		}) 


   				
		#===================================================================================
		# downloadfile handler
		  stdata <- reactive({
		    scdata <- validSeq$sumtable
			if(is.null(scdata)){return()}
			scdata
		  })		 

		  output$downloadSData <- downloadHandler(
		    filename = function() {
     			paste("Summary of Input Data", ".csv", sep = "")
    		},
		    content = function(file){
		    	write.csv(stdata(), file, row.names=FALSE)
		    }
		  )

		  ptdata <- reactive({
		    pcdata <- validSeq$psdata
			if(is.null(pcdata)){return()}
			pcdata
		  }) 


		  output$downloadPtable <- downloadHandler(
		    filename = function() {
     			paste("Summary of Preferred Codon", ".csv", sep = "")
    		},
		    content = function(file){
		    	write.csv(ptdata(), file, row.names=FALSE)
		    }
		  )

		  filedata <- reactive({
		    ctdata <- validSeq$ptable
			if(is.null(ctdata)){return()}
			ctdata
		  }) 


		  # downloadHandler() takes two arguments, both functions.
		  # The content function is passed a filename as an argument, and
		  #   it should write out data to that filename.
		  output$downloadData <- downloadHandler(
		    filename = function() {
     			paste("Codon Usage Table", ".csv", sep = "")
    		},
		    content = function(file){
		    	write.csv(filedata(), file, row.names=FALSE)
		    }
		  )

		  

		  #==========================PLOT=================

		  output$plot <- renderPlot({
		    #analyzedCodon <- list()
		    #plotdata <- read.csv("summary_of_input_data.csv", header = TRUE)
		    plotdata <- validSeq$sumtable
		    #colnames(plotdata) <- c("Sequence_no","Sequence(s)", "Number_of_Codons_Analyzed")
		    #View(plotdata)
		    #print(plotdata)
		    #print(is.null(plotdata[,3]))
		    codonAnalyzed <- as.numeric(plotdata[,3])
		    SequenceNo <- as.numeric(plotdata[,1])
		    freq <- count(codonAnalyzed)
		    codonMean <- mean(codonAnalyzed)

			if(radiodata()!=0){	
				#print(validSeq$invalid)
				#print(validSeq$ret)
				if(validSeq$invalid == FALSE && validSeq$ret == 2){
		    		#qplot(data=plotdata, aes(plotdata$codonAnalyzed)) + geom_histogram()
		    		#xlim=c(100,700),
		    		#prob = TRUE  #for density instead of frequency
		    		hist(codonAnalyzed, main="Histogram for Number of Codon Analyzed", xlab="Codon Analyzed", col="peachpuff", prob = TRUE, border="black",  las=1, breaks=5)
		    		abline(v = codonMean, col = "red")
		    		lines(density(codonAnalyzed), lwd = 2, col = "blue")
		    		#curve(codonAnalyzed, lwd = 2, col = "blue")
		    		abline(v = median(codonAnalyzed), col = "green", lwd = 2)
		    		legend(x = "topright", # location of legend within plot area
						 c("Density plot", "Mean", "Median"),
						 col = c("blue", "red", "green"),
						 lwd = c(2, 2, 2))
		    	}
			}else{
				#print("no plot")
				#plot.new()
				#require(cowplot)
				#theme_set(theme_cowplot(font_size=12)) # reduce default font size
				#View(mpg)
				#plot.mpg <- ggplot(mpg, aes(x = cty, y = hwy, colour = factor(cyl))) + 
				#  geom_point(size=2.5)
				#View(diamonds)
				#plot.diamonds <- ggplot(diamonds, aes(clarity, fill = cut)) + geom_bar() +
				#  theme(axis.text.x = element_text(angle=70, vjust=0.5))
				#plot.codon <- ggplot(plotdata, aes(codonAnalyzed, fill = 'codonAnalyzed')) + geom_bar() +
			
				#plot_grid(plot.mpg, plot.diamonds,  labels = c('A', 'B'))
			}
		  })

		  #========================= RUNNING TIME PLOT ===============================
		  output$rplot <- renderPlot({
		  	print("plot for running time..")
		  	print(validSeq$timediff)
		  	#rtimelist = dim(0)


		  	if(is.null(validSeq$timediff) == TRUE){
		  		print("NULL TIME")
		  	}else{
		  		#rtimelist <- data.frame(validSeq$timediff)
		  		#tfile <- read.table("runtime_vs_filesize.csv", header = TRUE, sep="")
		  		#View(tfile)
		  		#plot(runtime)
		  		runtime = list()
		  		runtime = c(1.24619, 2.711698, 4.276526, 7.130111, 8.649589, 11.5595, 13.458, 15.9524, 21.47348, 22.37731, 24.59345, 26.57147, 32.95003, 37.10814, 42.20026)
		  		totalCodonAnalyzed = c(314701, 629402, 944103, 1258804, 1573505, 1888206, 2202907, 2517608, 2832309, 3147010,3461711, 3776412, 4091113, 4405814, 4720515)
		  		filesize = 1:15

		  		runtime.tb <- data.frame(filesize, runtime, totalCodonAnalyzed)
		  		View(runtime.tb)

		  		require(cowplot)
				theme_set(theme_cowplot(font_size=12)) # reduce default font size
				plot.rtime <- ggplot(runtime.tb, aes(x = filesize, y = runtime, colour = factor(runtime))) + geom_point(size=2.5)
				plot.tcodon <- ggplot(runtime.tb, aes(x = totalCodonAnalyzed, y = runtime, colour = factor(totalCodonAnalyzed))) + geom_point(size=2.5)								
				plot_grid(plot.rtime, plot.tcodon, labels = c('Runtime(secs) vs. Filesize(mb)', 'Runtime(secs) vs. Total Number of Codon Analyzed'))
		  	}
	
		  })
		
    }),

	options= list(options(shiny.maxRequestSize=30*1024^2))
  
)
