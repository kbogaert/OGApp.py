#!/usr/bin/python
#
# Written by Kenny A. Bogaert
#
# Copyright (C) 2014 by Kenny A. Bogaert
#
# All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the 
# Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
# 
#
#DEVELOPERNOTES:
# problems with biopython using ncbi-blast+ 26
# some difficultie with passing seqrecs through a queue; leading to the processes being zombie process: work around put all elements in a list and write them as a seqrecord again in the main tread

import sys
import getopt
import subprocess
import os 
import pysam 
import time
import operator
import math
from multiprocessing import Process, Queue, Manager
from subprocess import call
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import Record
#~ import numpy as np

def find_orfs_with_trans(seq, trans_table, min_protein_length, leg_frame=0):
	if int(verbose)>1:
		print '			<<<find_orfs_with_trans>>>'
		print '				-find_orfs_with_trans-options: '
		print '					-trans_table: ', trans_table
		print '					-min_protein_length: ', min_protein_length
		print '					-leg_frame: ', leg_frame
	list_out = []
	seq_len = len(seq)
	spec_frame=math.fabs(leg_frame)-1#want iterator gaat daar van 0 tot 2 ipv 1 tot 3
	if leg_frame >0:
		spec_strand=+1
	elif leg_frame < 0: 
		spec_strand=-1
	elif leg_frame ==0: #if not specified
		spec_frame=0
		spec_strand=0
	#~ print 'spec_frame is ', spec_frame
	#~ print 'spec_strand is ', spec_strand
	for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
		#~ print 'strand is ', strand
		for frame in range(3):
			#~ print 'frame is ', frame
			if int(frame)==int(spec_frame) and int(strand)==int(spec_strand) or spec_frame == 0 and spec_strand == 0:
				#~ print 'start'
				sequence=nuc[frame:]
				while len(sequence)%3: #sequence has to be a multiplication of three - otherwise possible errors in future Biopython versions
					sequence=sequence+'N'
				trans = str(sequence.translate(trans_table))
				#~ print 'trans is ', trans
				trans_len = len(trans)
				aa_start = 0
				aa_end = 0
				while aa_start < trans_len:
					aa_end = trans.find("*", aa_start) # finds a string starting from aa_start, returns -1 if not found
					#~ print 'aa_end is ', aa_end
					if aa_end == -1:
						aa_end = trans_len
					if aa_end-aa_start >= min_protein_length:
						if strand == 1:
							start = frame+aa_start*3
							end = min(seq_len,frame+aa_end*3+3)
						else:
							start = seq_len-frame-aa_end*3-3
							end = seq_len-frame-aa_start*3  
						length = math.fabs(end-start)
						list_out.append((length, start, end, strand,seq[start:end]))
					aa_start = aa_end+1
	list_out.sort()
	list_out.reverse()
	return list_out



def blastdebug(fasta):
	#introduce a face fasta seq as blast parser contains a bug
	fastaname=os.path.splitext(fasta)[0] # cut off the extention of the filename 
	fasta_aux=str(fastaname)+'_aux.fasta' #name of the auxilary file
	fasta_aux_h=open(fasta_aux, "w")
	fasta_aux_h.write(">falseseq if-this-seq-becomes-visible-in-the-output-a-biopython-bug-has-likely-been-removed \nNNN\n")
	fasta_h=open(fasta, "r")
	for line in fasta_h:
		fasta_aux_h.write(line)
	fasta_h.close()
	fasta_aux_h.close()
	return fasta_aux
	
def blasthitorder(seqrec, fastadb, name):
	if int(verbose) > 1: 
		print '	<<< blasthitorder >>>'
		print '		blasthitorderoption: '
		print'			- seqrec: ', seqrec
		print '			- fastadb: ', fastadb
		print '			- name: ', name
	# fastadb is the fasta file to which already a database has been made previously
	##write
	oga_tmpseqfile='OGA_tmp_'+str(name)+'.fasta'
	oga_tmpseqxml='OGA_tmp'+str(name)+'.xml'
	oga_seqfilehandle=open(oga_tmpseqfile, 'w')
	SeqIO.write(seqrec, oga_seqfilehandle, "fasta")
	oga_seqfilehandle.close()
	oga_tmpseqfile_debugged=blastdebug(oga_tmpseqfile)
	subprocess.call(['rm '+str(oga_tmpseqfile)], shell=True)
	oga_tmpseqfile=oga_tmpseqfile_debugged
	
	##blast
	#~ print 'blastx -query '+str(oga_tmpseqfile)+' -db '+str(fastadb)+' -out '+str(oga_tmpseqxml) + ' -outfmt 5'
	outcome=subprocess.call(['blastx -query '+str(oga_tmpseqfile)+' -db '+str(fastadb)+' -out '+str(oga_tmpseqxml) + ' -outfmt 5  -num_alignments 1000000000' ], shell=True)
	if outcome!=0:
		print 'ERROR: blast failed'
		print ' seqrec.seq is ', seqrec.seq
		errorvalue=1
		return errorvalue
	##parse
	xml2handle=open(oga_tmpseqxml, 'r')
	blast_records = NCBIXML.parse(xml2handle)
	blast_record=blast_records.next()
	reciprocalrank=0
	frame=(0,0)
	scores2={}
	for blast_record in blast_records:
		if blast_record.alignments: 
			#~ print blast2_record.query
			#~ print alignment.hit_def
			for alignment in blast_record.alignments:
				#~ print 'check'
				for hsp in alignment.hsps:
					scores2[alignment.hit_def]=hsp.score
				sortedscores2 = sorted(scores2.iteritems(), key=operator.itemgetter(1), reverse=True) #list of tuples
				#~ print 'sorted scores is', sortedscores2
				sortedscores2list=[] #list of scores
				for k in sortedscores2:
					sortedscores2list.append(k[1])
				#~ print 'sortedscores2list is', sortedscores2list
				sortedset2=sorted(list(set(sortedscores2list)), reverse=True)
				hspcheck=0
				for hsp2 in alignment.hsps:
					#~ print 'name is ', name
					#~ print 'hitdef is', alignment.hit_def.split()[0]
					if alignment.hit_def.split()[0]==name:
						reciprocalrank=1
						for el in sortedset2:  #iterate over all the scores of the reciprocal blastrecord starting from the one with the biggest score
							if hsp2.score==el:
								#~ print 'reciprocank is', reciprocalrank
								frame=hsp.frame
								#~ break
							else:
								reciprocalrank+=1
							#~ modus=1
							#~ if int(hspcheck)==0: #wirte only the best hitting hsp
								#~ print blast_record.query.split(' ', 1)[0], 'has a RBH with ', alignment.hit_def, ' with rank ', reciprocalrank
						hspcheck=1
						check=1
	xml2handle.close()	
	time.sleep(0.01)
	try: 
		subprocess.call(['rm '+str(oga_tmpseqxml)], shell=True)
		subprocess.call(['rm '+str(oga_tmpseqfile)], shell=True)
	except: 
		pass
	
	if int(verbose) > 1: 
		print '	<<< end blasthitorder >>>'
	return [reciprocalrank, frame]
	
def dividefasta(fasta, chunks):
	if int(verbose)>1: 
		print '			<<< dividefasta >>>'
		print '				dividefasta-options: '
		print '					- fasta: ', fasta
		print '					- chunks: ', chunks
	fasta_h=open(fasta, 'r')
	fastacut=os.path.splitext(fasta)[0]
	#~ print 'grep -c "^>" '+str(fasta)
	number=subprocess.check_output(['grep -c "^>" '+str(fasta)], shell=True)
	print '		total amount of sequences is', number
	modulus = int(number) % int(chunks)
	#~ print modulus
	fold=(int(number)-modulus)/int(chunks)
	firstfold=fold+modulus
	count=0 # will remember which number of sequences already passed
	cur_chunk=1
	chunkname=str(fastacut)+'_'+str(cur_chunk)+'.fasta'
	chunkhandle=open(chunkname, 'w')
	cur_fold=firstfold
	chunknames=[chunkname] #contains the list of fastafilenames of the different files
	for seqrec in SeqIO.parse(fasta_h, "fasta"):
		count+=1
		if count <= cur_fold :
			SeqIO.write(seqrec, chunkhandle, "fasta")
		if count==cur_fold:
			#~ print 'next chunk', cur_chunk, cur_fold
			if cur_chunk <= int(chunks)-1: 
				cur_chunk+=1
				nextchunkname=str(fastacut)+'_'+str(cur_chunk)+'.fasta'
				chunknames.append(nextchunkname)
				chunkhandle=open(nextchunkname, 'w')
				cur_fold+=fold
	if int(verbose)>1:
		print '			<<< end dividefasta >>>'
	return chunknames
	
#~ class tblastnthread (threading.Thread):
	#~ def __init__(self, threadID,threadName, file,  contigfile, xmlfile, num_alignments, evalue):
		#~ threading.Thread.__init__(self)
		#~ self.file=file
		#~ self.contigfile = contigfile
		#~ self.xmlfile = xmlfile
		#~ self.num_alignments=num_alignments
		#~ self.evalue=evalue
	#~ def run(self): 
		#~ tblastn(self.name, self.file, self.contigfile, self.xmlfile, self.num_alignments, self.evalue) 
		#~ if int(verbose)>1:
			#~ print "			Exiting " + self.name
		
	
def tblastn(threadName, file, contigfile, xmlfile, num_alignments, evalue):
	if int(verbose)>1:
		print '			start thread', threadName, '\n'
	#~ if int(verbose)>1:
		#~ print'			<<< start tblastn thread >>>'
		#~ print '				tblastn options: '
		#~ print '					- file: ', file
		#~ print '					- contigfile: ', contigfile
		#~ print '					- xmlfile: ', xmlfile
		#~ print '					- num_alignments: ', num_alignments
	### UNCOMMENT THIS LINE: 
	print 'blasting with commandline: ', 'tblastn -query '+str(file)+' -db '+str(contigfile)+' -out '+str(xmlfile) + ' -outfmt 5'  + ' -evalue '+str(evalue)+' -num_alignments '+str(num_alignments)
	subprocess.call(['tblastn -query '+str(file)+' -db '+str(contigfile)+' -out '+str(xmlfile) + ' -outfmt 5'  + ' -evalue '+str(evalue)+' -num_alignments '+str(num_alignments) ], shell=True)
	#~ subprocess.call(['ls '+str(xmlfile)], shell = True)
	if int(verbose)>1:
		print '			<<< end tblastn thread >>>'

#~ class ogathread (threading.Thread):
	#~ def __init__(self, queue, wl_queue, threadID, threadName, file, xml, threadnumber, testlimit, evalue, contigfile):
		#~ threading.Thread.__init__(self)
		#~ self.file=file
		#~ self.xml=xml
		#~ self.threadnumber=threadnumber
		#~ self.testlimit=testlimit
		#~ self.evalue=evalue
		#~ self.contigfile=contigfile
		#~ self.queue=queue
		#~ self.wl_queue=wl_queue
	#~ def run(self):
		#~ oga(self.name,self.queue, self.wl_queue, self.file, self.xml, self.threadnumber, self.testlimit, self.evalue, self.contigfile)
		
#~ def oga_wrap(threadName, queue, wl_queue, orf_queue, allorfs_queue, file, xmlfile, threadnumber, testlimit, evalue, origcontigfile, number):
	#~ try:
		#~ oga(threadName, queue, wl_queue, orf_queue, allorfs_queue, file, xmlfile, threadnumber, testlimit, evalue, origcontigfile, number)
	#~ except:
		#~ print('%s: %s' % (threadName, traceback.format_exc()))
		
def oga(threadName, coorth_file, writtenlist_glob, orf_file, allorfs_file, file, completeref, xmlfile, threadnumber, testlimit, evalue, origcontigfile, number):
	#	writtenlist_glob: all contignames that were selected as coorthologous singlets, or were used for cap3 contigs (the others should be written to a file with the contigs for which no orthologous ref sequence has been found) 
	#	coorth_file: coorthologous sequences left after reciprocal blast
	#	allorfs_file: all orfs after the orf cap3 step
	#	orf_file: only longest orfs of the resulting orfs after cap3 step
	if int(verbose)>1: 
		print 'handling is ', file
		print '	with corresponding xmlfile ', xmlfile
	#make a copy of the contigfile for each oga process
	origcontigfile_base=os.path.splitext(origcontigfile)[0]
	contigfile=str(origcontigfile_base)+'_'+str(threadnumber)+'.fasta'
	##~ uncomment this
	subprocess.call([ 'cp '+str(origcontigfile)+' '+str(contigfile)], shell=True)
	xmlhandle=open(xmlfile, 'r')
	blast_records = NCBIXML.parse(xmlhandle)
	last_record = blast_records.next
	if int(verbose)>1:
		print 'make db', file
	## UNCOMMENT THIS
	subprocess.call(['makeblastdb -in '+str(file)+ ' -dbtype prot > /dev/null'], shell=True)
	refprotcount=0
	ogacontigcount=1
	ogaorfcap3contigcount=1
	writtenlist_loc=[] # all contigs that will be written in the first output file
	brcount=0 #counts the amount of blast records parsed, used for tracking progress
	orftempfile='orftempfile_'+str(threadnumber)+'.fasta' # file where all the longest orfs will be written for cap3
	orfcap3_allfile='orfcap3file_'+str(threadnumber)+'.fasta' # file where the resultant cap3 contigs and remaining singlets will be concatenated
	coorthhandle=open(coorth_file, 'w')
	allorfshandle=open(allorfs_file, 'w')
	orfhandle=open(orf_file, 'w')
	for blast_record in blast_records: # parse over all records of the xmlfile of this thread
		refprotcount+=1
		current_percent=(float(refprotcount)/float(number))*100
		if current_percent > 100:
			current_percent=float(100)
		print '	'+str(threadName)+': handling '+str(blast_record.query)+' - ' + str(current_percent)+' % completed'
		refprot=blast_record.query.split()[0]
		if blast_record.alignments:
			refprot_descr=blast_record.query
			oga_homologygroupfile='OGA_'+str(threadnumber)+'_homologygroup_'+str(blast_record.query.split()[0]) +'.fasta'
			oga_homologygrouphandle=open(oga_homologygroupfile, "w")
			selectionlist=[]
			for alignment in blast_record.alignments: # parse over all alignments of blast record
				hitid=alignment.hit_def.split()[0]
				if int(verbose)>2:
					print '	has a hit with ', hitid
				scoremax=0
				framemax='none'
				for hsp in alignment.hsps: # parse over each hsp and add hitid if not already in selectionlist
					if hsp.expect < evalue: 
						if int(verbose)>2:
							print '		evalue is ', hsp.expect
						if hitid not in selectionlist:
							selectionlist.append(hitid)
			contighandle = open(contigfile, "r")
			if len(selectionlist)==0:
				print threadName, ': ', refprot, ' had no alignments with evalue larger than ', evalue
				continue
			## extract all matching contigs
			for seq_rec in SeqIO.parse(contighandle, "fasta"):
				id = seq_rec.id.split()[0]
				if int(verbose)>2:
					print threadName, ': id is ', id
				if str(id) in selectionlist:
					SeqIO.write(seq_rec, oga_homologygrouphandle, "fasta")
			contighandle.close()
			oga_homologygrouphandle.close()
			#FIRST CAP3
			if int(verbose)<2:
				code=subprocess.call(['cap3 '+str(oga_homologygroupfile) +' -t 31 > /dev/null'], shell=True)
				if code!=0:
					print 'WARNING:'
					print 'cap3 - round 1 - error handling ', oga_homologygroupfile
					print 'error: ', output
					print 'skipping cap3 step'
					subprocess.call(['cp '+str(oga_homologygroupfile) +' '+str(oga_homologygroupfile)+'.errorhomologygroupcap3'], shell=True)
					cap3contighandle=open(str(oga_homologygroupfile), 'r')
			else: 
				code=subprocess.call(['cap3 '+str(oga_homologygroupfile)+' -t 31'], shell=True)
				if code!=0:
					print 'WARNING:'
					print 'cap3 - round 1 - error handling ', oga_homologygroupfile
					print 'skipping cap3 step'
					subprocess.call(['cp '+str(oga_homologygroupfile) +' '+str(oga_homologygroupfile)+'.errorhomologygroupcap3'], shell=True)
					cap3contighandle=open(str(oga_homologygroupfile), 'r')				
			#handle the contigs of cap3
			try:
				cap3contighandle
			except: 
				cap3contighandle=open(str(oga_homologygroupfile)+'.cap.contigs', 'r')
			## if blasts back to same contig --> presumable (co-)ortholog
			longestorfs=[] # here all the longest orfs for each sequencen in first output file will be stored
			firstcontigcheck=0
			firstcontignumber=0
			
			## handling the contigs of cap3
			reciproccheck=0 # will stay 0 if there is no hit with reciprocal order 1
			for cap3rec in SeqIO.parse(cap3contighandle, "fasta"):
				if firstcontigcheck==0: # in order to remember the number of the first contig (later used for naming of the longest orf)
					firstcontignumber=ogacontigcount
				if len(cap3rec.seq)>3: 
					list=blasthitorder(cap3rec, completeref, refprot)
				else: 
					print 'ERROR: no sequence in seqrec: ', cap3rec
					continue
				if list==1:  # has a value of 1 instead of being a list with 2 characters if error was raised
					print 'error raised in blasthitorder when handling', refprot
					continue
				reciprocalorder=list[0]
				frame=list[1][0]
				# indication for no reciprocal blasthit
				if int(verbose)>=0:
					print cap3rec.id, ' has a reciprocal order of ', reciprocalorder, ' to refprot ', refprot
				ogacontigcount+=1
				firstcontigcheck=1
				if reciprocalorder == int(1): #if reciprocal order is 1
					reciproccheck=1
					cap3rec.id='Contig_OGA_t'+str(threadnumber)+'_'+str(ogacontigcount)
					cap3rec.description='co-orth_esil '+str(refprot_descr) + ' cap3product '
					SeqIO.write(cap3rec, coorthhandle, "fasta")
					### check orfs
					orf_list = find_orfs_with_trans(cap3rec.seq, 1, 0, frame)
					orf_list.append(cap3rec.id) # 6th element of the list is the name of the sequences
					
					####333333#######
					#adapt here minimum length of longest orf and store the orginal contigs in a list to remove them againfrom the written list
					###############
					
					if len(orf_list[0][4])>0:
						longestorfs.append(orf_list[0])
					else: 
						if int(verbose)> 1:
							print '		longest orf was 0 bases long of ', cap3rec.id 
						
				##handling the SINGLETS of cap3
			
			cap3singlethandle=open(str(oga_homologygroupfile)+'.cap.singlets', 'r')
			singletlist=[] # all contignames that were singlets (whether RBH or not)
			singletwrittenlist=[] # all singletnames that will be written in output
			for cap3rec in SeqIO.parse(cap3singlethandle, "fasta"):
				if firstcontigcheck==0: # in order to remember the number of the first contig (later used for naming of the longest orf
					firstcontignumber=ogacontigcount
				singletlist.append(cap3rec.id.split()[0])
				if len(cap3rec.seq)>0: 
					list=blasthitorder(cap3rec, completeref, refprot)
				else: 
					print 'ERROR: no sequence in seqrec: ', cap3rec
					continue
				reciprocalorder=list[0]
				frame=list[1][0]
				ogacontigcount+=1
				firstcontigcheck=1 #checks whether there has been written one contig for the current reference sequence
				if int(verbose)>=0:
					print cap3rec.id, ' has a reciprocal order of ', reciprocalorder, ' to refprot ', refprot
				if reciprocalorder == int(1):
					reciproccheck=1 
					singletwrittenlist.append(cap3rec.id.split()[0])
					cap3rec.id='Contig_OGA_t'+str(threadnumber)+'_'+str(ogacontigcount)
					cap3rec.description='co-orth_esil '+str(refprot_descr) + ' cap3singlet '
					SeqIO.write(cap3rec, coorthhandle, "fasta")
					##check orfs
					orf_list = find_orfs_with_trans(cap3rec.seq, 1, 10, frame)
					orf_list.append(cap3rec.id) # 6th element of the list is the name of the sequences
					if len(orf_list[0][4])>0:
						longestorfs.append(orf_list[0]) #first element of the list is the longest one (~ sorted by length) 
					else: 
						if int(verbose)> 1:
							print '		longest orf was 0 bases long of ', cap3rec.id 
			cap3singlethandle.close()
			#for all sequences that were assembled:
			cap3list = [x for x in selectionlist if x not in singletlist]
			writtenlist_loc=writtenlist_loc+cap3list+singletwrittenlist
			if int(verbose)>1:
				print threadName+' writtenlist_loc is ', writtenlist_loc
			##select the sequence with the longest orf 
			lengthorf=0 
			curorf=''
			orfcount=0
			temporfhandle=open(orftempfile, 'w')
			for longorf in longestorfs:
				orfcount+=1
				orfrec=SeqRecord(longorf[4])
				orfrec.id='ORF'+str(orfcount)
				orfrec.name=orfrec.id
				orfrec.description = 'orf_co-orth_ref '+ str(refprot_descr)
				SeqIO.write(orfrec, temporfhandle, "fasta")
			temporfhandle.close()
			if reciproccheck==0: #if not recprocal hits were encountered one should not do the cap3 steps on empty files
				print 'no reciprocal best hits found for ', refprot
				continue
			if len(longestorfs)==0:
				print '		WARNING: no longest orfs found for ', refprot
				continue
			#~ if int(verbose)<2:
				#~ code=subprocess.call(['cap3 '+str(orftempfile) +' -t 31 > /dev/null'], shell=True)
				#~ trial=1
				#~ while code!=0 and trial<3: #try 3 times until it works (else: print the warning and skip it)
					#~ trial+=1
					#~ code=subprocess.call(['cap3 '+str(orftempfile) +' -t 31 > /dev/null'], shell=True)
				#~ subprocess.call(['cat '+str(orftempfile) +'.cap.contigs '+str(orftempfile)+'.cap.singlets > '+ str(orfcap3_allfile)], shell=True)
				#~ orfcap3handle = open(orfcap3_allfile, 'r')
				#~ if code!=0:
					#~ print 'WARNING:'
					#~ print 'cap3 - round 2 - error handling ', refprot
					#~ print 'longestorfs was', longestorfs
					#~ print 'skipping cap3 step'
					#~ subprocess.call(['cp '+str(orftempfile) +' '+str(refprot)+'.error-orfcap3'], shell=True)
					#~ cap3contighandle=open(str(orftempfile), 'r')
					#~ orfcap3_allfile=orftempfile
			code=subprocess.call(['cap3 '+str(orftempfile)+' -t 31 > /dev/null '], shell=True)
			trial=1
			while code!=0 and trial<3:
				trail+=1
				code=subprocess.call(['cap3 '+str(orftempfile)+' -t 31 >/dev/null'], shell=True)
			subprocess.call(['cat '+str(orftempfile) +'.cap.contigs '+str(orftempfile)+'.cap.singlets > '+ str(orfcap3_allfile)], shell=True)
			orfcap3handle = open(orfcap3_allfile, 'r')
			if code!=0:
				print 'WARNING:'
				print 'cap3 - round 2 - error handling ', refprot
				print 'skipping cap3 step'
				subprocess.call(['cp '+str(orftempfile) +' '+str(oga_homologygroupfile)+'.error-orfcap3'], shell=True)
				orfcap3handle=open(str(orftempfile), 'r')
				orfcap3_allfile=orftempfile
			cur_maxlength=float(0)
			curseq=0
			#~ curseq_count=0
			##selects longest orf and writes all orfs in the orf file:
			for sr in SeqIO.parse(orfcap3handle, "fasta"):
				#~ curseq_count+=1 # count the amount of sequences (just for debugging purposes, can be removed) 
				ogaorfcap3contigcount+=1
				length=len(sr.seq)
				sr.id='Contig_OGA_orf_t'+str(threadnumber)+'_'+str(ogaorfcap3contigcount)
				SeqIO.write(sr, allorfshandle, "fasta")
				if length > cur_maxlength:
					cur_maxlength=length
					curseq=sr
			if curseq==0:
				print 'ERROR: no sequences in orf; there should be sequences'
				print 'here is the content of orfcap3handle'
				subprocess.call(["cat "+str(orfcap3_allfile)], shell=True)
				subprocess.call(["cp "+str(orfcap3_allfile) +' '+ str(orfcap3_allfile)+'.error'], shell=True)
				break
			curseq.id='Contig_OGA_t'+str(threadnumber)+'_'+str(firstcontignumber)+'-'+'t'+str(threadnumber)+'_'+str(int(ogacontigcount)-1)					
			curseq.name=orfrec.id
			curseq.description= 'max_orf_co-orth_ref '+ str(refprot_descr)+' | longest orf corresponding to contig_OGA_t'+str(threadnumber)+'_'+str(firstcontignumber)+' until contig_OGA_t'+str(threadnumber)+'_'+str(ogacontigcount-1)
			SeqIO.write(curseq, orfhandle, "fasta")
			#~ subprocess.call(["rm "+str(oga_homologygroupfile)+'*'], shell=True)
		else: 
			print threadName, ': ', refprot, ' had no alignments'	#~ print 'doing thread ', threadName
	#~ print 'writtenlist_loc is ', writtenlist_loc
	#~ try: 
		#~ subprocess.call(["rm "+str(orftempfile)+' ' +str(orfcap3_allfile)], shell=True)
	#~ except: 
		#~ pass
	writtenlist_glob += writtenlist_loc
	coorthhandle.close()
	orfhandle.close()
	allorfshandle.close()
	#~ print 'end thread : ', threadName
	print '		<<< end oga >>>'
	return


if __name__ == "__main__":
	argv=sys.argv[1:]
#parameters
	# still have to ensure all parameters are filled in; if not get an error message - now very cryptic error message
##default values:
	clean = 1
	xmlfile = 0
	evalue = '10e-5'
	num_alignments = 10000
	num_threads = 1
	testlimit='none'
	global verbose
	verbose = 0
	#add part for evalue, num_alighnments
	try:
		opts, args = getopt.getopt(argv,"hi:r:x:e:v:n:",['help', 'input=', 'ref=', 'xml=', 'evalue=', 'verbose=', 'threads='])
	except getopt.GetoptError:
		print 'error mainsub - getopt'
		sys.exit(3)
	for opt, arg in opts:
		if opt == '-h':
			print
			print 'USAGE: python kbOGApp2.py -i <in contig file> -r <reference> -x <xmlfile> -n <threads> -e <evalue (10e-5)' 
			print
			print 'options:'
			print' '
			sys.exit()
		elif opt in ("-i", "--input"): # is the fasta file
			contigfile = arg			
		elif opt in ("-x", "--xml"): # is the fasta file
			xmlfile = arg	
		elif opt in ("-r", "--ref"): # is the fasta file
			reffile = arg	
		elif opt in ("-n", "--number"): # is the fasta file
			num_threads = arg
		elif opt in ("-v", "--verbose"): # is the fasta file
			verbose = arg		
		elif opt in ("-e", "--evalue"): # is the fasta file
			evalue = arg		
	#~ contigfile=blastdebug(contigfile)
	contigfilecut=os.path.splitext(contigfile)[0]
	reffilecut=os.path.splitext(reffile)[0]
	if int(num_threads) > 1: 
		chunknames=dividefasta(reffile, num_threads)
		#~ print 'chunknames is ', chunknames
	else: 
		chunknames=[reffile]
	xmlfiles=[]
	print '	step x/x: blasting'
	processes1=[]
	if int(num_threads)==0:
		print 'Error: invalid number of threads specified', num_threads
		sys.exit()
	if num_threads>0:
		## uncomment
		subprocess.call(['makeblastdb -in '+str(contigfile)+ ' -dbtype nucl > /dev/null'], shell=True)
		number=0
		for file in chunknames:
			number+=1
			filecut=os.path.splitext(file)[0]
			xmlfile='OGA_tblastn_'+str(filecut)+'X'+str(contigfile)+'.xml'
			threadname='thread-' + str(filecut)
			p1= Process(target=tblastn, args=(str(threadname), file,contigfile, xmlfile, num_alignments, evalue))
			p1.start()
			#~ time.sleep(0.01)
			processes1.append(p1)
			xmlfiles.append(xmlfile)
	print '		waiting for tblastn to finish ...'
	for p_1 in processes1:
		p_1.join()
	print ' end tblastn processes'
	time.sleep(1)
	
	manager=Manager()
	writtenlist=manager.list([])
	
	num_seqs=subprocess.check_output(['grep -c "^>" '+str(reffile)], shell=True)
	num_seqs_per_thread=math.fabs(int(num_seqs)/int(num_threads))
	processes2=[]
	coorth_contigs_files=[]
	orf_files=[]
	allorf_files=[]
	print ' end blasting; result is ', xmlfiles
	#make blastdb for entire ref
	subprocess.call(['makeblastdb -in '+str(reffile)+ ' -dbtype prot > /dev/null'], shell=True)
	for number in range(int(num_threads)):
		xmlfile=xmlfiles[number]
		coorth_contigs_file='coorth_contigs_threadtemp_'+str(number)
		coorth_contigs_files.append(coorth_contigs_file)
		orf_file='orf_threadtemp_'+str(number)
		orf_files.append(orf_file)
		allorf_file='allorfs_threadtemp_'+str(number)
		allorf_files.append(allorf_file)
		refchunk=chunknames[number]
		threadname='thread-' + str(number)
		p2=Process(target=oga, args=(str(threadname), coorth_contigs_file, writtenlist, orf_file, allorf_file, refchunk, reffile, xmlfile, number, testlimit, evalue, contigfile, num_seqs_per_thread))
		p2.start()
		processes2.append(p2)
	print 'before join'
	for p_2 in processes2:
		p_2.join()
	print 'Debug: end joining'
	
	coorth_contigs_files_string=' '.join(coorth_contigs_files) 
	output_orthologs=str(contigfilecut) + '.OGA.fasta'
	print 'Writing coortholog file: ', output_orthologs
	subprocess.call(['cat ' +str(coorth_contigs_files_string) + ' > ' + str(output_orthologs) ], shell=True)
	#~ subprocess.call(['rm ' +str(coorth_contigs_files_string)  ], shell=True)
	
	orf_files_string=' '.join(orf_files)
	output_longestorf=str(contigfilecut) + '.OGA.longestorf.fasta'
	print 'Writing longest orf file: ', output_longestorf
	subprocess.call(['cat ' +str(orf_files_string) + ' > ' + str(output_longestorf) ], shell=True)
	#~ subprocess.call(['rm ' +str(orf_files_string) ], shell=True)

	allorf_files_string=' '.join(allorf_files)
	output_orfs=str(contigfilecut) + '.OGA.orf.fasta'
	print 'Writing all orf file: ', output_orfs
	subprocess.call(['cat ' +str(allorf_files_string) + ' > ' + str(output_orfs) ], shell=True)
	#~ subprocess.call(['rm ' +str(allorf_files_string) ], shell=True)
	
	contigfilehandle=open(contigfile, 'r')
	unwrittenhandle=open(str(contigfilecut) + '.OGA.unused.fasta', 'w')
	
	#~ print 'writtenlist is ', writtenlist
	print 'Writing the unused contigs file: ' + str(contigfilecut) + '.OGA.unused.fasta'
	for rcrd in SeqIO.parse(contigfilehandle, "fasta"):
		if not rcrd.id.split()[0] in writtenlist:
			SeqIO.write(rcrd, unwrittenhandle, "fasta")
	contigfilehandle.close()
	unwrittenhandle.close()
