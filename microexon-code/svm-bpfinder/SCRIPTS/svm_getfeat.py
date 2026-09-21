#!/usr/bin/env python
import sys,os,socket,shutil,time,string,commands,copy,re,math,random

##########################################################################
base_dir = os.path.dirname(os.path.realpath(__file__))

fname=sys.argv[1]
species=sys.argv[2]
maxslen=int(sys.argv[3])
mindist3ss=int(sys.argv[4])

fmfilename=base_dir+'/../MODELS/'+species+'.foreg.o2.matrix'
bmfilename=base_dir+'/../MODELS/'+species+'.backg.o2.matrix'


##########################################################################

def read_fasta(filename):
	flist=[]
	ffile=open(filename,'r')
	for line in ffile:
		cont=line.replace('\n','')
		if cont!='':
			if cont[0]=='>':
				flist.append([cont[1:],''])
			else:
				flist[-1][1]+=cont.lower()
		else:
			continue
	ffile.close()
	return flist

def cut_3p_fasta(flist,clen):
	for rec in flist:
		rec[1]=rec[1][-clen:]	
	
def get_agez(seq,offset=12):
    sseq=seq[:-offset].split('ag')
    return offset+len(sseq[-1])

def get_bps(seq,mindist3ss=int(sys.argv[4])):
    bsdict={}
    for x in range(5,len(seq)-(mindist3ss-1)):
        pseq=seq[x-5:x+4]
        if pseq[3]=='t' and pseq[5]=='a' and 'n' not in pseq:
            bsdict[len(seq)-x]=pseq
    return bsdict

def get_post_bp(seq,bspos,cut3=4):
    return seq[-(bspos-1):-cut3]
    
def nucFreq(seq):
    nucdict={}
    for nuc in unique(seq):
        nucdict[nuc]=float(0)
    for nuc in seq:
        nucdict[nuc]+=1
    for nuc in nucdict:
        nucdict[nuc]/=len(seq)
    return nucdict
    

##########################################################################

def read_matrix(mfilename):
    mdict={}
    mfile=open(mfilename,'r')
    mfile.readline()
    for line in mfile:
        linearray=line.replace('\n','').split('\t')
        mdict[linearray[0].lower()]=map(lambda x: float(x),linearray[1:])
    mfile.close()
    return mdict

def word_prob(word,myords,matrix):
    myords.reverse()
    myindex=range(len(myords))[::-1]
    probs=[]
    for y in range(len(myindex)):
        subword=word[myindex[y]-myords[y]:myindex[y]+1]
        subsubword=subword[:-1]
        if myords[y]==0:
            probs.append(matrix[subword][myindex[y]])
        else:
            subwordprob=matrix[subword][myindex[y]]
            subsubwordprob=matrix[subsubword][myindex[y]-1]
            if subwordprob<>0:
                probs.append(subwordprob/subsubwordprob)
            else:
                probs.append(float(0))
    probs.reverse()
    dprob=1
    for p in probs:
        dprob=dprob*p
    myords.reverse()
    return probs,dprob

def log_score(fmfilename,bmfilename,morder,seqlist):
    scores=[]
    pm=read_matrix(fmfilename)
    nm=read_matrix(bmfilename)
    for putbp in seqlist:
        scores.append(math.log(word_prob(putbp,morder,pm)[1]/word_prob(putbp,morder,nm)[1],2))
    return scores

##########################################################################

def scounter(seq,x,ntype,scount):
    if ntype=='y':
        scount=scount+1           
        for y in range(x+1,len(seq)):
            if seq[y]=='c' or seq[y]=='t':
                scount+=1
            else:
                break               
    elif ntype=='r':
        scount=scount+1           
        for y in range(x+1,len(seq)):
            if seq[y]=='a' or seq[y]=='g':
                scount+=1
            else:
                break
    return scount

def pcounter(seq,x,ntype,scount,spos):
    spos+=1
    pcount=0
    if ntype=='y': 
        for y in range(x+1+(scount-spos),len(seq)):
            if seq[y]=='a' or seq[y]=='g':
                pcount+=1
            else:
                break               
    elif ntype=='r':
        for y in range(x+1+(scount-spos),len(seq)):
            if seq[y]=='c' or seq[y]=='t':
                pcount+=1
            else:
                break               
    return pcount


def testppt(ppt):
    if len(ppt)>=9 or ppt.count('t')>=5:
        return 1
        #return ppt.upper()
    else:
        return 0
        #return ppt.lower()

def sco_ppt(ppt):
    pptlen=len(ppt)
    yperc=float(ppt.count('t')+ppt.count('c'))/pptlen
    sco=0
    sdict={'a':-2,'c':2,'t':3,'g':-2,'n':-2}
    for x in ppt:
        sco+=sdict[x]
    
    return [pptlen,yperc,sco]


def ext_ppt(seq):
    pptlist=[]
    seq=seq.lower()
    state=0 ; ntype='n'
    ycounter=0 ; rcounter=0
    prevn=0 ; selfn=0 ; postn=0
    ppt=''
    
    ##########
    
    for x in range(len(seq)):
        if state==0:
            ppt=''
            pptstart=x
            if seq[x]=='a' or seq[x]=='g' or seq[x]=='n':
                if ntype=='y':
                    rcounter=0
                prevn=ycounter
                selfn=scounter(seq,x,'r',rcounter)
                postn=pcounter(seq,x,'r',scounter(seq,x,'r',rcounter),rcounter)
                rcounter+=1
                ntype='r'
                
            elif seq[x]=='t' or seq[x]=='c':
                ppt+=seq[x]
                if ntype=='r':
                    ycounter=0
                prevn=rcounter
                selfn=scounter(seq,x,'y',ycounter)
                postn=pcounter(seq,x,'y',scounter(seq,x,'y',ycounter),ycounter)
                ycounter+=1
                ntype='y'
                state=1

        ##########
                
        elif state==1:
            if seq[x]=='t' or seq[x]=='c':
                ppt+=seq[x]
                if ntype=='r':
                    ycounter=0
                prevn=rcounter
                selfn=scounter(seq,x,'y',ycounter)
                postn=pcounter(seq,x,'y',scounter(seq,x,'y',ycounter),ycounter)
                ycounter+=1
                ntype='y'

            elif seq[x]=='a' or seq[x]=='g' or seq[x]=='n':
                if ntype=='y':
                    rcounter=0
                prevn=ycounter
                selfn=scounter(seq,x,'r',rcounter)
                postn=pcounter(seq,x,'r',scounter(seq,x,'r',rcounter),rcounter)
                rcounter+=1              
                ntype='r'
                
                ##########

                if selfn<=float(prevn+postn)/4 and selfn<=prevn and selfn<=postn and selfn<=2:
                    ppt+=seq[x]
                else:
                    if seq[x]=='g' and x<len(seq)-1 and x>0 and seq[x-1]=='t' and seq[x+1]=='t':
                        if x==1 or x==2:
                            if seq[x-1:x+4]=='tgtgt':
                                ppt+=seq[x]
                            else:
                                if testppt(ppt):
                                    return [pptstart+1,pptstart+len(ppt)]+sco_ppt(ppt)+[ppt]
                                state=0
                        elif x==len(seq)-2 or x==len(seq)-3 :
                            if seq[x-3:x+2]=='tgtgt':
                                ppt+=seq[x]
                            else:
                                if testppt(ppt):
                                    return [pptstart+1,pptstart+len(ppt)]+sco_ppt(ppt)+[ppt]
                                state=0
                        else:
                            if seq[x-1:x+4]=='tgtgt' or seq[x-3:x+2]=='tgtgt':
                                ppt+=seq[x]
                            else:
                                if testppt(ppt):
                                    return [pptstart+1,pptstart+len(ppt)]+sco_ppt(ppt)+[ppt]
                                state=0            
                    else:
                        if testppt(ppt):
                            return [pptstart+1,pptstart+len(ppt)]+sco_ppt(ppt)+[ppt]
                        state=0
    if state==1:
        if testppt(ppt):
            return [pptstart+1,pptstart+len(ppt)]+sco_ppt(ppt)+[ppt]                    
    return pptlist



##########################################################################
##########################################################################	

infasta=read_fasta(fname)
cut_3p_fasta(infasta,maxslen)
for rec in infasta:
	intid=rec[0] ; intseq=rec[1]
	intagez=get_agez(intseq)
	
	putbps=get_bps(intseq)
	
	putbppos=putbps.keys()
	putbppos.sort(reverse=True)
	
	putbpseqs=[]
	putbpycont=[]
	putbppptdist=[] ; putbppptlen=[] ; putbppptscore=[]
	agnumber=[]
	
	for pos in putbppos:
		putbpseqs.append(putbps[pos])
		postseq=get_post_bp(intseq,pos)#[:-4]
		putbpycont.append(float(postseq.count('c')+postseq.count('t'))/len(postseq))
		pptres=ext_ppt(postseq)
		if pptres<>[]:
			putbppptdist.append(pptres[0])
			putbppptlen.append(pptres[2])
			putbppptscore.append(pptres[4])
		else:
			putbppptdist.append(len(postseq))
			putbppptlen.append(0)
			putbppptscore.append(0)
		agnumber.append(postseq[12:-8].count('ag'))
	putbpscore=log_score(fmfilename,bmfilename,[0,1,1,0,2,0,2,1,1],putbpseqs)
	for x in range(len(putbppos)):
		feats=['0',
				'1:'+str(putbpscore[x]),
				'2:'+str(putbpycont[x]),
				'3:'+str(putbppptdist[x]),
				'4:'+str(putbppptscore[x]),
				'#']
		extras=[intid,
				str(intagez),
				str(putbppos[x]),
				putbpseqs[x],
				str(putbppptlen[x])]
		
		print ' '.join(feats)+'\t'.join(extras)

