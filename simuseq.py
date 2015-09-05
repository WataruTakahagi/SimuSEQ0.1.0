#!/usr/bin/env python
import sys
import datetime
import matplotlib.pyplot as plt
import numpy as np
import csv
import re
import math
import os
import os.path
import shutil
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
ENDC = '\033[0m'

class SimuSEQ:
    def __init__(self):
        pass

    def readseq(self, name="NONAME"):
        data = open(os.getcwd()+"/data/"+'sequence.txt', 'r')
        self.seq = data.read()
        result = self.seq
        return result
    
    def makeprobe(self, seq):
        probe_size = 50
        list_operon = []
        front = []
        back = []
        probe = []
        fb = []
        list_probe = []
        nowdir = os.getcwd()
        os.chdir(nowdir+"/data")
        with open('operons.csv', 'r') as f:
            os.chdir(nowdir)
            reader = csv.reader(f)
            header = next(reader)
            for row in reader:
                list_operon.append(row[0])
                front.append(int(row[1])-1)
                back.append(int(row[2])-1)
                fb.append(row[3])
            for i in range(len(list_operon)):
                probe.append(seq[int(front[i]):int(front[i])+probe_size]) #Front probe
                #probe.append(seq[int(back[i])-probe_size:int(back[i])])   #Rear probe
                #probe.append(seq[int(front[i]):int(back[i])])             #Full probe
            list_probe = probe
            for i in range(len(fb)):
                if fb[i] == "-1":
                    ss = []
                    for s in list_probe[i]:
                        ss.append(s)
                    ss.reverse()
                    list_probe[i] = ''.join(ss)
        patlist = list_probe
        opf = open('probe.txt', 'w')
        for i in range(len(patlist)):
            pattern = patlist[i]
            text = seq
            iterator = re.finditer(pattern, text)
            span_start = []
            span_end = []
            num = 0
            for match in iterator:
                num += 1
                mg = match.group().replace("t", "u")
                ms = `match.start()`
                me = `match.end()`
                print GREEN+"{:<45}".format(list_operon[i])+" "+BLUE+mg+" "+RED+"{:<8}".format(ms)+" "+YELLOW+me+ENDC
                opf.write(list_operon[i]+','+mg+"\n")
                span_start.append(ms)
                span_end.append(me)
                break
        opf.close()

    def microarray(self, seq, bgcolor='black'):
        f2 = open('result.txt')
        file = f2.read()
        rev = []
        for line_r in file: 
            rev.append(line_r)
        rev.reverse()
        rev = ''.join(rev)
        list_operon = []
        front = []
        back = []
        probe = []
        nowdir = os.getcwd()
        os.chdir(nowdir+"/data")
        with open('operons.csv', 'r') as f3:
            os.chdir(nowdir)
            reader = csv.reader(f3)
            header = next(reader)
            for row in reader:
                list_operon.append(row[0])
                front.append(row[1])
                back.append(row[2])
            for i in range(len(list_operon)):
                probe.append(seq[int(front[i]):int(back[i])])
        namelist = []
        patlist = []
        pat = open('probe.txt', 'r')
        for line in pat:
            patlist.append(line.rsplit(",")[1].replace("\n",""))
            namelist.append(line.rsplit(",")[0])
        f4 = open('microarray.txt', 'w')
        for i in range(len(patlist)):
            pattern = patlist[i]
            text = file
            iterator = re.finditer(pattern, text)
            span_start = []
            span_end = []
            num = 0
            for match in iterator:
                mp = match.group()
                ms = `match.start()`
                me = `match.end()`
                num += 1
                print "{:<3}".format(`num`)+" "+GREEN+"{:<45}".format(namelist[i])+" "+BLUE+mp+" "+RED+"{:<8}".format(ms)+" "+YELLOW+me+ENDC
                span_start.append(`match.start()`)
                span_end.append(`match.end()`)
            f4.write(`num`+","+namelist[i]+","+mp+","+ms+","+me+"\n")
        f4.close()
        Exlev = []
        gname = []
        f5 = open('microarray.txt', 'r')
        reader = csv.reader(f5)
        for row in reader:
            Exlev.append(int(row[0]))
            gname.append(row[1])
        n = len(Exlev)
        exmax = max(Exlev)
        exmin = min(Exlev)
        ex = exmax - exmin
        sq = math.sqrt(n)
        for i in range(n):
            if i*i >= n:
                sq = i
                break
        X = []
        Y = []
        for i in range(sq):
            for j in range(sq):
                X.append(i+0.5)
                Y.append(j+0.5)
        ax = plt.axes([0.025, 0.025, 0.95, 0.95], axisbg=bgcolor)
        ax.set_xlim(0,sq)
        ax.set_ylim(0,sq)
        ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(1.0))
        ax.yaxis.set_major_locator(plt.MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(1.0))
        ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0')
        ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0')
        ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0')
        ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        T = Exlev
        for i in range(sq*sq-n):
            T.append(0)
        plt.scatter(X, Y, s=75, c=T, alpha=.5, cmap=plt.cm.get_cmap('jet'))
        plt.xticks(())
        plt.yticks(())
        plt.colorbar()
        plt.savefig("microarray.png")
        plt.close()
        f5.close()
        
    def rank(self, microarray='microarray.txt'):
        f7 = open(microarray, 'r')
        reader = csv.reader(f7)
        list = []
        time = []
        n = 0
        for i in reader:
            list.append(int(i[0]))
            time.append(n)
            n += 1
        a = sorted(list)
        plt.bar(time,a,align="center",color = "blue")
        plt.xlim([0,2600])
        plt.ylim([0,140])
        plt.savefig("rank.png")
        plt.close()
        f7.close()
        
    def fasta(self, seq):
        print YELLOW+"outputting data"+ENDC
        f6 = open('result.fasta', 'w')
        d = datetime.datetime.today()
        f6.write('>simulated data|%s\n' % (d))
        j = 1
        for i in seq:
            if i == "a":
                #sys.stdout.write(BLUE + i + ENDC)
                f6.write(i)
            elif i == "t":
                #sys.stdout.write(RED + i + ENDC)
                f6.write(i)
            elif i == "g":
                #sys.stdout.write(GREEN + i + ENDC)
                f6.write(i)
            elif i == "c":
                #sys.stdout.write(YELLOW + i + ENDC)
                f6.write(i)
            if j % 80 == 0:
                f6.write('\n')
                #print
            j += 1
        f6.write('\n')
        #print
        f6.close()

    def makedata(self, dirname="SimuSEQ_result"):
        if os.path.exists(os.getcwd()+"/"+dirname):
            swt = 1
            print BLUE+dirname+RED+" already exists !!"+ENDC
            dirname = raw_input(YELLOW+"Please input other name : "+ENDC)
            while swt == 1:
                if os.path.exists(os.getcwd()+"/"+dirname):
                    dirname = raw_input(RED+"ERROR "+GREEN+"Please input other name : "+ENDC)
                else:
                    break
        print GREEN+"Project END"+ENDC
        os.mkdir(dirname)
        if os.path.exists(os.getcwd()+'/probe.txt'): shutil.move('probe.txt',os.getcwd()+"/"+dirname)
        if os.path.exists(os.getcwd()+'/microarray.txt'): shutil.move('microarray.txt',os.getcwd()+"/"+dirname)
        if os.path.exists(os.getcwd()+'/rank.png'): shutil.move('rank.png',os.getcwd()+"/"+dirname)
        if os.path.exists(os.getcwd()+'/result.fasta'): shutil.move('result.fasta',os.getcwd()+"/"+dirname)
        if os.path.exists(os.getcwd()+'/microarray.png'): shutil.move('microarray.png',os.getcwd()+"/"+dirname)

    def auto(self, setting=[1,1,1,1,0,1]):
        #setting = [1,1,1,1,0,1] # 1 : on, 0 : off
        SimuSEQ_list = ["SimuSEQ().readseq()","SimuSEQ().makeprobe()","SimuSEQ().microarray()","SimuSEQ().rank()","SimuSEQ().fasta()","SimuSEQ().makedata()"]
        for i in range(len(setting)):
            if setting[i] == 1:
                print YELLOW+"{:<22}".format(SimuSEQ_list[i])+ENDC+" : "+GREEN+"on"+ENDC
            elif setting[i] == 0:
                print YELLOW+"{:<22}".format(SimuSEQ_list[i])+ENDC+" : "+RED+"off"+ENDC
        stt = raw_input(YELLOW+"continue?"+ENDC+" ("+BLUE+"y"+ENDC+"/"+RED+"n"+ENDC+") : ")
        if stt == "y":
            if setting[0] == 1: seq = SimuSEQ().readseq()
            if setting[1] == 1: SimuSEQ().makeprobe(seq)
            if setting[2] == 1: SimuSEQ().microarray(seq)
            if setting[3] == 1: SimuSEQ().rank()
            if setting[4] == 1: SimuSEQ().fasta(seq)
            if setting[5] == 1: SimuSEQ().makedata()
        elif stt == "n":
            pass
        else:
            print RED+"ERROR"+ENDC
            
