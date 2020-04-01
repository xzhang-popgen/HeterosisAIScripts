#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:44:17 2019

Description: This script takes inputs from command line that specify the following parameters: 

genomic segment, demographic model, population growth pattern in recipient population, 
dominance coefficient for deleterious mutations, whether a human hs relationship is used, 
selection coefficient for the adaptive mutation, scaling factor, and the number of simulation replicates

The script updates a slim script according to the above parameters, runs slim program, and extracts 
adaptive introgression summary statistics in non-overlapping 50kb windows 

@author: xinjunzhang
"""

import msprime, pyslim, os, random, pandas, sys, itertools,argparse,glob,matplotlib.pyplot,re
import numpy as np
from multiprocessing import Manager, Pool


parser = argparse.ArgumentParser(description="A script for running slim and computing summary statistics in 50kb windows across a given chromosome")
parser.add_argument('-g', '--gene', action="store", dest="gene_id",
                        help="which simulation batch, default: 1; range 1-26",
                        default=1, type=int)
parser.add_argument('-h', '--dominance', action="store", dest="dominance_id",
                        help="dominance index, default: 1; range 0-100 (h=0-1); if running neutral model, value=200;",
                        default=1, type=int)
parser.add_argument('-m', '--model', action="store", dest="model_id",
                        help="model index, default: 0; 0=m0, 1=mh",
                        default=0, type=int)
parser.add_argument('-p', '--popsize', action="store", dest="growth_id",
                        help="growth index, default: 4; range:1-4",
                        default=4, type=int)
parser.add_argument('-d', '--hs', action="store", dest="hs_id",
                        help="growth index, default: 0; 0=use dominance input, 1=use hs relationship",
                        default=0, type=int)
parser.add_argument('-n', '--nscale', action="store", dest="nscale_id",
                        help="scaling factor index, default: 10",
                        default=10, type=int)
parser.add_argument('-s', '--selcoeff', action="store", dest="selcoeff_id",
                        help="adaptive mutation selection strength index, default: 0; range: 0-1 (s=0-0.1)",
                        default=0, type=int)
parser.add_argument('-r', '--rep', action="store", dest="numrep_id",
                        help="number of simulation replicates, default: 200",
                        default=200, type=int)                                                
args = parser.parse_args()

whichgene = int(args.gene_id) 
dominance = round(float(args.dominance_id)/100,2 ) #convert h-index to h value: 50 -> 0.5
model = int(args.model_id)
growth = int(args.growth_id)
hs = int(args.hs_id)
nscale = int(args.nscale_id)
m4s = float(args.selcoeff_id/100) #convert s-index to s: 1 -> 0.01
num_reps = int(args.numrep_id)

#sample command: python3 run_slim_get_stats.py -g 1 -h 0 -m 1 -p 4 -d 0 -n 10 -s 1 -r 200


def calc_p1ancestry (treepath, admpop, popsize,t_sinceadm,model):
    ts = pyslim.load(treepath)
    any_ancestry = ancestry_p_varies(ts,admpop,popsize,t_sinceadm,model)
    meanp1 = sum(any_ancestry)/len(any_ancestry)
    return meanp1

def ancestry_p_varies(ts,pop,nsize,duration,model): #pop=source pop
    n=nsize*2
    mixtime=duration
    p = [x.id for x in ts.nodes() if ((x.population == int(pop)) and (x.time == mixtime))] #source pop
    if model ==1:
    	today = [x.id for x in ts.nodes() if ((x.population == 4) and (x.time == 0))] #assuming p4 is recipient
	elif model ==0:
		today = [x.id for x in ts.nodes() if ((x.population == 3) and (x.time == 0))]

    tree_p = [sum([t.num_tracked_samples(u) for u in p])/n
               for t in ts.trees(tracked_samples=today, sample_counts=True)]

    return tree_p

def ancestry_local (treepath):
    starts, ends, subpops = [], [], []
    ts = pyslim.load(treepath)
    for tree in ts.trees(sample_counts=True):
        subpop_sum, subpop_weights = 0, 0
        for root in tree.roots:
            leaves_count = tree.num_samples(root) - 1
            subpop_sum += tree.population(root) * leaves_count
            subpop_weights += leaves_count
        starts.append(tree.interval[0])
        ends.append(tree.interval[1])
        subpops.append(subpop_sum / float(subpop_weights))

    x = [x for pair in zip(starts, ends) for x in pair]
    y = [x for x in subpops for _ in (0, 1)]   
    matplotlib.pyplot.plot(x, y)
    matplotlib.pyplot.show()
    return x,y #x=genome positions; y = ancestry

def ancestry_position_writeout (treepath,n,admpop,popsize,t_sinceadm,region_name):
    ts = pyslim.load(treepath)
    starts=[]
    ends=[]

    for x in ts.trees():
        starts.append(x.interval[0])
        ends.append(x.interval[1])

    outfilename = DIR_anc+ region_name+str(dominance)+ "_"+str(model)+ "_"+str(growth)+ "_"+str(m4s)+ "_"+str(hs) + "_"+str(n) + '.ancestry'
    outfile = open(outfilename, 'w')
    outfile.write('start,end,ancestry\n')
    
    p1ancestry = ancestry_p_varies(ts,admpop,popsize,t_sinceadm)

    for start, end, anc in zip(starts, ends, p1ancestry):
        outfile.write('{0},{1},{2}\n'.format(start, end, anc))

    outfile.close()

def write_ancestry (DIR_tree, output_anc_file):
    tree_all = glob.glob(DIR_tree+'*.trees')
    with open(output_anc_file, 'w') as outfile:
        for file in tree_all:
            x,y = ancestry_local (file)
    
            for item in x:
                outfile.write("%s\t" % item)
            outfile.write("\n")
            for item in y:
                outfile.write("%s\t" % item)
            outfile.write("\n")


def load_data_slim(file_path,len_genome,adm_gen,end_gen): # load slim's output 
    pos_den, hapMat_den,freqp4_before,freqp4_after = get_pos_hap(file_path,'p2',len_genome,end_gen)
    pos_afr, hapMat_afr,freqp4_before,freqp4_after = get_pos_hap(file_path,'p1',len_genome,end_gen)
    pos_nonafr, hapMat_nonafr,freqp4_before,freqp4_after = get_pos_hap(file_path,'p4',len_genome,end_gen)
    pos_preadm, hapMat_preadm,freqp4_before,freqp4_after = get_pos_hap(file_path,'p2',len_genome,adm_gen)

    return pos_den, hapMat_den,pos_afr, hapMat_afr,pos_nonafr, hapMat_nonafr,pos_preadm, hapMat_preadm,freqp4_before,freqp4_after


def get_pos_hap(file_path,pop_id,len_genome,gen_time): #get pos and hapMat for a given pop from slim output 
    infile = open(file_path,'r')
    end=0
    while end==0:
        line = infile.readline()
        if line[0:5]=='#OUT:': #output lines           
            fields = line.split()
            out_type = fields[2]
            pop = fields[3]
            gen = fields[1]
            if out_type=='SM' and pop==pop_id and int(gen) == gen_time: #ms lines
                num_indiv = int(fields[4])
                infile.readline() #skip //
                infile.readline() #skip segsites
                pos = (np.array(infile.readline().split()[1:]).astype(float) * len_genome).astype(int)
                mult_mut_pos = find_mult_mut_pos(pos)+1
                pos = np.delete(pos,mult_mut_pos)
                hapMat = np.zeros((num_indiv,len(pos)),dtype=int)
                for indiv in range(0,num_indiv):
                    hap = np.array(list(infile.readline())[:-1]).astype(int)
                    hap = np.delete(hap,mult_mut_pos)
                    hapMat[indiv] = hap
                freqp4_before = 0
                freqp4_after = 0
                end=1                
    infile.close()       
    return pos, hapMat,freqp4_before,freqp4_after
     
                
def find_mult_mut_pos(pos): #find repeating mutations and remove them
    dist = np.array([pos[i+1]-pos[i] for i in range(0,len(pos)-1)])
    mult_mut_pos = np.where(dist==0)[0]
    return mult_mut_pos

def find_ai_site (segfile): #find an exon in the mid-range of the segment to insert AI mutation
    segs = open(segfile)
    starts = []
    ends = []
    total = 0
    for line_counter,line in enumerate (segs):
        if line[0:4]=="exon":      
            fields = line.split()
            if (int(fields[1]) >= 2200000) & (int(fields[2]) <=2800000):
                starts.append(fields[1])
                ends.append(fields[2])
                total+=1        
    any_exon = random.choice(range(0,total))
    window_start = int(starts[any_exon])
    window_end = int(ends[any_exon])
    segs.close()
    return window_start,window_end #return exon start and end position


def insert_anc_alleles (allpos,pos,hap):
    for site in allpos:
        if site not in pos:
            insertidx = np.searchsorted(pos,site)
            pos = np.insert(pos,insertidx,site)
            hap = np.insert(hap, insertidx, 0, axis=1)
    return pos, hap
 
    
def calc_derived_freq (pop_hap):
    popfreq = np.sum(pop_hap, axis=0)
    popfreq = popfreq/ float(pop_hap.shape[0])
    return popfreq

def vSumFunc(other_hap, currentArchi):
	current_hap = np.array([p1_hapw[currentArchi,]])
	div = np.zeros(other_hap.shape)
	ones = np.ones((other_hap.shape[0],1))
	current_hap = current_hap
	current_hap_extended = np.dot(ones, current_hap)
	div = np.logical_xor(current_hap_extended == 1, other_hap == 1)
	return np.add.reduce(div, 1) 
    
def calc_stats (file_path,len_genome,adm_gen,end_gen):
    pos_den, hapMat_den,pos_afr, hapMat_afr,pos_nonafr, hapMat_nonafr,pos_preadm, hapMat_preadm,freqp4_before,freqp4_after = load_data_slim(file_path,len_genome,adm_gen,end_gen)

    p1_pos = pos_den
    p2_pos = pos_afr
    p3_pos = pos_nonafr
    p1_hap = hapMat_den
    p2_hap = hapMat_afr
    p3_hap = hapMat_nonafr
    
    all_pos = np.unique(np.concatenate((p1_pos,p2_pos,p3_pos)))

    p1_pos,p1_hap = insert_anc_alleles(all_pos,p1_pos,p1_hap)
    p2_pos,p2_hap = insert_anc_alleles(all_pos,p2_pos,p2_hap)
    p3_pos,p3_hap = insert_anc_alleles(all_pos,p3_pos,p3_hap)


    allpos_bin = np.linspace(0,len_genome,int(len_genome/50000)) #windows of every 50kb

    allpos_digitized = np.digitize(all_pos, allpos_bin)


    Dstat_list = []
    fD_list = []
    Het_list = []
    divratioavg_list = []
    Q_1_100_q95_list =[]
    Q_1_100_q90_list=[]
    Q_1_100_max_list=[]
    U_1_0_100_list =[]
    U_1_20_100_list =[]
    U_1_50_100_list =[]
    U_1_80_100_list =[]
    
    pos_start = []
    pos_end = []
    

    for w in range(1,100):
        these_pos = all_pos[allpos_digitized==w]
        pos_start.append(min(these_pos))
        pos_end.append(max(these_pos))
        
        if len(these_pos)>1:
            
            these_pos_idx = np.nonzero(np.in1d(all_pos,these_pos))[0]
            p1_hapw = p1_hap[:,these_pos_idx]
            p2_hapw = p2_hap[:,these_pos_idx]
            p3_hapw = p3_hap[:,these_pos_idx]
        
            p1_freqw = calc_derived_freq (p1_hapw)
            p2_freqw = calc_derived_freq (p2_hapw)
            p3_freqw = calc_derived_freq (p3_hapw)
        
            abbavecw = (1.0 - p2_freqw)*p3_freqw*p1_freqw      
            babavecw = p2_freqw*(1.0 - p3_freqw)*p1_freqw
            abbacountsw = np.sum(abbavecw)
            babacountsw = np.sum(babavecw)
        
            if (abbacountsw + babacountsw > 0):
                Dstatw = (abbacountsw - babacountsw) / (abbacountsw + babacountsw)
            else:
                Dstatw = float('nan')

            Dstat_list.append(Dstatw)   

            checkfd1 = (p3_freqw > p1_freqw)
            abbafd1 = (1.0 - p2_freqw)*p3_freqw*p3_freqw
            babafd1 = p2_freqw*(1.0 - p3_freqw)*p3_freqw
            checkfd2 = (p3_freqw < p1_freqw)
            abbafd2 = (1.0 - p2_freqw)*p1_freqw*p1_freqw
            babafd2 = p2_freqw*(1.0 - p1_freqw)*p1_freqw
            abbafd = checkfd1 * abbafd1 + checkfd2 * abbafd2
            babafd = checkfd1 * babafd1 + checkfd2 * babafd2
            abbafdcounts = np.sum(abbafd)
            babafdcounts = np.sum(babafd)
            if (abbafdcounts + babafdcounts > 0):
                fD = (abbacountsw - babacountsw) / (abbafdcounts - babafdcounts)
            else:
                fD = float('nan')

            fD_list.append(fD)


            hetvec = 2 * p3_freqw * (1.0 - p3_freqw)
            Het = np.sum(hetvec) /50000
            Het_list.append(Het)

        
            divratio = []

            for archi in range(0, p1_hapw.shape[0]): #iterate over 0-99 haps; 100 total
                divarchintro = vSumFunc(p3_hapw, archi)
                divarchintro = divarchintro.astype("float")
                divarchnonintro = vSumFunc(p2_hapw, archi)
        
                divarchnonintro = divarchnonintro.astype("float")      
                for comb in itertools.product(divarchintro,divarchnonintro): 
                    if comb[1] != 0:
                        divratio.append(comb[0]/comb[1])
            divratioavg = float(sum(divratio)) / float(len(divratio)) 
            divratioavg_list.append(divratioavg)


            ArcHomoDer = (p1_freqw == 1)
            NonAdm_1 = (p2_freqw < 0.01) 
            ArcHomoDerANDNonAdm_1 = (ArcHomoDer & NonAdm_1)
            DerFreqs_NonAdm_1 = p3_freqw[np.where(ArcHomoDerANDNonAdm_1 == True)]
            if DerFreqs_NonAdm_1.size > 0:
                Q_1_100_q95 = np.percentile(DerFreqs_NonAdm_1,95)
                Q_1_100_q90 = np.percentile(DerFreqs_NonAdm_1,90)
                Q_1_100_max = np.max(DerFreqs_NonAdm_1)
            else:
                Q_1_100_q95 = float('nan')
                Q_1_100_q90 = float('nan')
                Q_1_100_max = float('nan')

            Q_1_100_q95_list.append(Q_1_100_q95)
            Q_1_100_q90_list.append(Q_1_100_q90)
            Q_1_100_max_list.append(Q_1_100_max)

            U_1_0_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0) )
            U_1_20_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.2) )
            U_1_50_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.5) )
            U_1_80_100 = ( ArcHomoDerANDNonAdm_1 & (p3_freqw > 0.8) )

            U_1_0_100 = np.sum(U_1_0_100)
            U_1_20_100 = np.sum(U_1_20_100)
            U_1_50_100 = np.sum(U_1_50_100)
            U_1_80_100 = np.sum(U_1_80_100)
       
            U_1_0_100_list.append(U_1_0_100)
            U_1_20_100_list.append(U_1_20_100)
            U_1_50_100_list.append(U_1_50_100)
            U_1_80_100_list.append(U_1_80_100)
        else:
            Dstat_list.append(float('nan'))
            fD_list.append(float('nan'))
            Het_list.append(float('nan'))
            divratioavg_list.append(float('nan'))
            Q_1_100_q95_list.append(float('nan'))
            Q_1_100_q90_list.append(float('nan'))
            Q_1_100_max_list.append(float('nan'))      
            U_1_0_100_list.append(float('nan'))
            U_1_20_100_list.append(float('nan'))
            U_1_50_100_list.append(float('nan'))
            U_1_80_100_list.append(float('nan'))
        
        
    return pos_start,pos_end,freqp4_before,freqp4_after,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list

def update_par_file(region_name,temp_par, new_par,model,growth,dominance,nscale,m4s,hs,insert_ai):
    oldfile = open(temp_par)
    newfile = open(new_par,'w')
    line_counter=0
    for line_counter, line in enumerate(oldfile):
        fields = line.split()
        
        if model ==0:        
        	if line_counter==1:
        		fields[1] = str(dominance)+");"
        	elif line_counter==2:
        		fields[1] = str(nscale)+");"
        	elif line_counter==3:
        		fields[1] = str(m4s)+"*n);"
        	elif line_counter==25:
            	fields[2] = 'readFile("/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/regions/sim_seq_info_'+str(region_name)+'.txt");'
        	elif line_counter==78:
        		fields[0] = str(int(100000/nscale))
        	elif line_counter==88:
        		fields[0] = str(int(100100/nscale))
        	elif line_counter==96:
        		fields[0] = str(int(100100/nscale))+":"
        	elif line_counter==113:
        		fields[0] = str(int(110000/nscale))
        	elif line_counter==115:
        		fields[0] = str(int(110000/nscale))
        	elif line_counter==121:
        		fields[0] = str(int(119950/nscale))
        	elif line_counter==124:
        		fields[0] = str(int(120000/nscale -1))
        	elif line_counter==133:
        		fields[0] = str(int(120000/nscale))
        	elif line_counter==138:
        		fields[0] = str(int(120000/nscale +1))
        	elif line_counter==146:
        		fields[0] = str(int(130000/nscale))        		        		
        	elif line_counter==150:
            	fields[0] = 'sim.treeSeqOutput("/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/tree/'+str(region_name)+'_m0.trees");'
        	new_line=str()    
        	for item in fields:
            	new_line = new_line+item+" "
        	newfile.write(new_line+'\n')

        elif model ==1:   #modelh     
        	if line_counter==1:
        		fields[1] = str(int(growth))+");"
        	elif line_counter==2:
        	    fields[1] = str(dominance)+");" 
        	elif line_counter==3:
        	    fields[1] = str(hs)+");"        	           		
        	elif line_counter==4:
        		fields[1] = str(nscale)+");"
        	elif line_counter==5:
        		fields[1] = str(m4s)+"*n);"
        	elif line_counter==23:
            	fields[2] = 'readFile("/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/regions/sim_seq_info_'+str(region_name)+'.txt");'        		
        	elif line_counter==77:
        		fields[0] = "1:"+str(int(89000/nscale))
        	elif line_counter==90:
        		fields[0] = str(int(73000/nscale))
        	elif line_counter==98:
        		fields[0] = str(int(74000/nscale))
        	elif line_counter==102:
        		fields[1] = str(int(insert_ai))+");"
        	elif line_counter==106:
        		fields[0] = str(int(74000/nscale)) + ":"
        	elif line_counter==118:
        		fields[1] = str(int(insert_ai))+");"
        	elif line_counter==124:
        		fields[0] = str(int(83400/nscale))
        	elif line_counter==128:
        		fields[0] = str(int(86960/nscale))
        	elif line_counter==133:
        		fields[0] = str(int(87400/nscale - 1))       		
        	elif line_counter==144:
        		fields[0] = str(int(87400/nscale))       		
        	elif line_counter==151:
        		fields[0] = str(int(88080/nscale))
        	elif line_counter==158:
        		fields[0] = str(int(88080/nscale))+ ":"+str(int(89000/nscale))
        	elif line_counter==183:
        		fields[0] = str(int(89000/nscale))        		
        	elif line_counter==187:
            	fields[0] = 'sim.treeSeqOutput("/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidatesh/tree/'+str(region_name)+'_mh.trees");'
        	
        	new_line=str()    
        	for item in fields:
            	new_line = new_line+item+" "
        	newfile.write(new_line+'\n')

    
    newfile.close()
    oldfile.close()
    

def calc_ancestry_window (ancestry_file,len_genome):
    infile = open(ancestry_file,'r')
    end_pos = []
    ancestry = []    
    line_counter=0
    for line_counter, line in enumerate(infile):
        fields = line.split(',')
        if fields[0] != "start":
            end_pos.append(float(fields[1]))
            ancestry.append(float(fields[2]))               
    infile.close()
    allpos_bin = np.linspace(0,len_genome,int(len_genome/50000)) #windows of every 50kb

    endpos_digitized = np.digitize(end_pos, allpos_bin)
    end_pos = np.array(end_pos)
    ancestry = np.array(ancestry)
    
    anc_window = []
    anc_pos = []
    
    for w in range(1,100):
        these_pos = end_pos[endpos_digitized==w]
        these_anc = ancestry[endpos_digitized==w]    
        
        if(len(these_pos))>0:
            anc_window.append(np.mean(these_anc))
            anc_pos.append(these_pos)
        else:
            anc_window.append(float('nan'))
            anc_pos.append(these_pos)
        
    return anc_window    

   

def run_slim_variable(n,q,r,dominance,nscale,m4s,model,growth,hs,insert_ai):
	
    region_name = region_all[r]
    segsize=5000000
     
    
    if model ==1:  
    	if dominance !=2:  
    		temp_par = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/slim/modelh_neg.txt"
    	elif dominance == 2:
    		temp_par = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/slim/modelh_neu.txt"
    	adm_gen = (87400-1)/nscale
    	end_gen = 89000/nscale
    	t_end = 1600/nscale -1    	
    	popsize=41080/nscale #recipient population size at the end of simulation

    	if growth ==1:
        	popsize = 550/nscale
    	elif growth ==2: 
        	popsize = 41080/nscale
    	elif growth ==3: 
        	popsize = 7300/nscale
    	elif growth ==4: 
        	popsize = 41080/nscale

    if model==0:
        if dominance !=2:  
    		temp_par = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/slim/model0_neg.txt"
    	elif dominance == 2:
    		temp_par = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/slim/model0_neu.txt"
    	adm_gen = 120000/nscale - 1
    	end_gen = 130000/nscale
    	t_end = 10000/nscale
    	popsize=1000/nscale


    new_par = DIR_par +"par_"+region_name+str(dominance)+str(model)+ str(n)+".txt"
    
    update_par_file(region_name,temp_par, new_par,model,growth,dominance,nscale,m4s,hs,insert_ai)
    
    slim_output = DIR_out +'OUT_'+region_name+str(dominance)+str(model)+ str(n)+".txt" 
        
    os.system('/u/home/x/xinjunzh/slim_build/slim %s > %s' %(new_par ,slim_output))
    
    treepath = DIR_tree+str(region_name)+str(dominance)+'.trees'
	
	if model==1:
    	meanp1 = calc_p1ancestry(treepath,2,popsize,t_end,model)
    	ancestry_position_writeout(treepath,n,2,popsize,t_end,region_name) #write out ancestry info

	elif model==0:
    	meanp1 = calc_p1ancestry(treepath,1,popsize,t_end,model)
    	ancestry_position_writeout(treepath,n,1,popsize,t_end,region_name) 
    
    
    pos_start,pos_end,freqp4_before,freqp4_after,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list = calc_stats(slim_output,segsize,adm_gen,end_gen)
    
    ancestry_file = DIR_anc+ region_name+str(dominance)+ "_"+str(model)+ "_"+str(growth)+ "_"+str(m4s)+ "_"+str(hs) + "_"+str(n) + '.ancestry'
    anc_window = calc_ancestry_window (ancestry_file,segsize) #get mean ancestry per 50kb window

    q.put([n,insert_ai,growth,meanp1,pos_start,pos_end,freqp4_before,freqp4_after,anc_window, Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list])
    #other parameter info are stored in the output file name
    
    os.system('rm '+slim_output)   
    os.system('rm '+treepath)  
    os.system('rm '+new_par)
    os.system('rm '+ancestry_file)


     
def write_to_file(windowfile_name,q):
    windowfile = open(windowfile_name,'w')

    while 1:
        q_elem = q.get()

        if q_elem=='kill': # break if end of queue
            print ('END OF SIMULATIONS')
            break

        [n,insert_ai,growth,meanp1,pos_start,pos_end,freqp4_before,freqp4_after,anc_window,Dstat_list, fD_list, Het_list, divratioavg_list,Q_1_100_q95_list,Q_1_100_q90_list,Q_1_100_max_list,U_1_0_100_list,U_1_20_100_list,U_1_50_100_list,U_1_80_100_list] = q_elem
        for i in range(len(Dstat_list)):
            windowfile.write("%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (n,insert_ai,growth,meanp1,pos_start[i],pos_end[i],freqp4_before,freqp4_after,anc_window[i],Dstat_list[i], fD_list[i], Het_list[i], divratioavg_list[i],Q_1_100_q95_list[i],Q_1_100_q90_list[i],Q_1_100_max_list[i],U_1_0_100_list[i],U_1_20_100_list[i],U_1_50_100_list[i],U_1_80_100_list[i]))

        windowfile.flush()
    
    windowfile.close()


#################################################################################
if __name__=='__main__':   
	#whichgene = 1
	#model = 1 # 1=modelh; 0=model0 #define these two with parseargument
	#growth = 4 
	#hs = 0 #0 = recessive or neutral; 1 = hs relationship
	#dominance = 0 #run the deleterious recessive model #if 2, run the neutral model
	#nscale = 10 #define scaling factor
	#m4s = 0.0 #adaptive selection strength
	#num_reps=200 #number of simulations per region

    region_all = ["chr11max","chr19region","chr3region","galnt18","hla","hyal2","krt71","nlrc5","oca2","pde6c","pou2f3","rnf34","sema6d","sgcb","sgcz","sipa1l2","slc16a11","slc19a3","slc5a10","stat2","tbx15","tlr1610","tnfa1p3","txn"]
	
	DIR_region = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/regions/"
	DIR_anc = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/ancestry/"
    DIR_out = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/output/"
    DIR_tree = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/tree/"
    DIR_par = "/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/slim/"
    
    r = int(whichgene-1)
    
     	
    region_name = region_all[r]
    
    window_start,window_end=find_ai_site (DIR_region+"sim_seq_info_"+str(region_name)+".txt")
    insert_ai = int((int(window_end)+int(window_start))/2) #find the position to insert AI mutation
    
    windowfile_name = '/u/home/x/xinjunzh/data/Slim_3/calculate_sumstats/sim_candidates/stat/'+region_name+"-dominance"+str(dominance)+"-model"+str(model)+"-growth"+str(growth)+"-hs"+str(hs)+"-ai"+str(m4s)+'_human_windows.txt'
    num_proc = 10
    manager = Manager()
    pool = Pool(processes=num_proc)
    q = manager.Queue()   
    watcher = pool.apply_async(write_to_file,(windowfile_name,q))    
    reps = range(0,num_reps)
    args_iterable = list(zip(reps,[q]*num_reps))
    for i in args_iterable:
        n=i[0]
        print(str(n))
        run_slim_variable(i[0],i[1],r,dominance,nscale,m4s,model,growth,hs,insert_ai)
        
    q.put('kill')
    pool.close()
    pool.join()
  

    print("END OF SIMULATION")





