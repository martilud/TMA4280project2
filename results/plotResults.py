import matplotlib.pyplot as plt
import numpy as np
# Assumes 'results.txt' is on the form
# mpiprocs threads n error time
show = False
if("PLOT CONVERGENCE PART 1"):
    filenames = ['results.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        nlist = []
        resultlist = []
        for line in data:
            line = line.split()
            nlist.append(int(line[2]))
            resultlist.append(float(line[3]))
        plt.loglog(nlist,resultlist)
    plt.grid()
    plt.legend([file.split('.')[0] for file in filenames]) # Python is the best language 
    plt.title("Convergence plot")
    plt.xlabel(r"$n$")
    plt.ylabel(r"error")
    plt.savefig("Convergence1.pdf")
    if show: plt.show()
    plt.clf()


if("PLOT TIME PART 1"):
    filenames = ['results.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        nlist = []
        resultlist = []
        for line in data:
            line = line.split()
            nlist.append(int(line[2]))
            resultlist.append(float(line[4]))
        plt.loglog(nlist,resultlist)
    plt.grid()
    plt.legend([file.split('.')[0] for file in filenames]) # Python is the best language 
    plt.title("Time plot")
    plt.xlabel(r"$n$")
    plt.ylabel(r"Time in seconds")
    plt.savefig("Time1.pdf")
    if show: plt.show()
    plt.clf()

if(not "PLOT CONVERGENCE PART 2"):
    filenames = ['zeta2.txt', 'mach2.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        ndict = {} 
        resultdict = {}
        for line in data:
            line = line.split()
            if(line[0] in ndict):
                ndict[line[0]].append(int(line[1]))
            else:
                ndict[line[0]] = [int(line[1])]
            if(line[0] in resultdict):
                resultdict[line[0]].append(float(line[2]))
            else:
                resultdict[line[0]] = [float(line[2])]
        for key in ndict:
            plt.loglog(ndict[key],resultdict[key])
        plt.grid()
        plt.legend([str(keys) + "processes" for keys in ndict]) 
        plt.title("Convergence plot of " + file.split('.')[0])
        plt.xlabel(r"$n$")
        plt.ylabel(r"$|\pi_n - \pi|$")
        plt.savefig(file.split('.')[0] + 'Convergence.pdf')
        if show: plt.show()
        plt.clf()

if(not "PLOT TIME PART 2"):
    filenames = ['zeta2.txt', 'mach2.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        ndict = {} 
        resultdict = {}
        for line in data:
            line = line.split()
            if(line[0] in ndict):
                ndict[line[0]].append(int(line[1]))
            else:
                ndict[line[0]] = [int(line[1])]
            if(line[0] in resultdict):
                resultdict[line[0]].append(float(line[3]))
            else:
                resultdict[line[0]] = [float(line[3])]
        for key in ndict:
            plt.loglog(ndict[key],resultdict[key])
        plt.grid()
        plt.legend([str(keys) + "processes" for keys in ndict]) 
        plt.title("Time plot of " + file.split('.')[0])
        plt.xlabel(r"$n$")
        plt.ylabel(r"Time in seconds")
        plt.savefig(file.split('.')[0] + 'Time.pdf')
        if show: plt.show()
        plt.clf()

if(not "PLOT CONVERGENCE PART 3"):
    filenames = ['zeta3.txt', 'mach3.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        ndict = {} 
        resultdict = {}
        for line in data:
            line = line.split()
            if(line[0] in ndict):
                ndict[line[0]].append(int(line[1]))
            else:
                ndict[line[0]] = [int(line[1])]
            if(line[0] in resultdict):
                resultdict[line[0]].append(float(line[2]))
            else:
                resultdict[line[0]] = [float(line[2])]
        for key in ndict:
            plt.loglog(ndict[key],resultdict[key])
        plt.grid()
        plt.legend([str(keys) + "processes" for keys in ndict]) 
        plt.title("Convergence plot of " + file.split('.')[0])
        plt.xlabel(r"$n$")
        plt.ylabel(r"$|\pi_n - \pi|$")
        plt.savefig(file.split('.')[0] + 'Convergence.pdf')
        if show: plt.show()
        plt.clf()

if(not "PLOT TIME PART 3"):
    filenames = ['zeta3.txt', 'mach3.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        ndict = {} 
        resultdict = {}
        for line in data:
            line = line.split()
            if(line[0] in ndict):
                ndict[line[0]].append(int(line[1]))
            else:
                ndict[line[0]] = [int(line[1])]
            if(line[0] in resultdict):
                resultdict[line[0]].append(float(line[3]))
            else:
                resultdict[line[0]] = [float(line[3])]
        for key in ndict:
            plt.loglog(ndict[key],resultdict[key])
        plt.grid()
        plt.legend([str(keys) + "processes" for keys in ndict]) 
        plt.title("Time plot of " + file.split('.')[0])
        plt.xlabel(r"$n$")
        plt.ylabel(r"Time in seconds")
        plt.savefig(file.split('.')[0] + 'Time.pdf')
        if show: plt.show()
        plt.clf()

if("PLOT TIME PART 4"):
    filenames = ['results.txt']
    for file in filenames:
        data = open(file, 'r')
        header = data.readline()
        ndict = {} 
        resultdict = {}
        for line in data:
            line = line.split()
            if(line[0] in ndict):
                ndict[line[0]].append(int(line[2]))
            else:
                ndict[line[0]] = [int(line[2])]
            # Store time
            if(line[0] in resultdict):
                resultdict[line[0]].append(float(line[4]))
            else:
                resultdict[line[0]] = [float(line[4])]
        for key in ndict:
            plt.loglog(ndict[key],resultdict[key])
        plt.grid()
        plt.legend([str(keys) + "processes" for keys in ndict]) 
        plt.title("Time plot of " + file.split('.')[0])
        plt.xlabel(r"$n$")
        plt.ylabel(r"Time in seconds")
        plt.savefig(file.split('.')[0] + 'Time.pdf')
        if show: plt.show()
        plt.clf()

