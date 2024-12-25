import sys
import subprocess
import time
import os
import glob

if len(sys.argv) < 3:
    print("Use: python threaded_runner file_list.txt cores [extra.txt]")
    print("file_list.txt: The list of parameter files to be used in p21_sim")
    print("cores: The number of cores to be used at once.")
    print("extra.txt: An additional parameter file to be used in every run.")
    print("If bax_module is also present, will run that too.")

files = open(sys.argv[1], 'r').readlines()
num_procs = int(sys.argv[2])

extra = ""
if len(sys.argv) >= 4:
    extra = sys.argv[3]


def runlist(python_file, arglist):
    proclist = []
    #Will terminate any subprocess that runs for more than a minute.
    #Solver occasionally hangs.
    timelist = [0]*num_procs
    for i in range(num_procs):
        proclist.append(None)

    cur_process = 0
    last_started = 0
    for f in arglist:
        while True:
            if timelist[cur_process] == 60:
                proclist[cur_process].kill()
            if proclist[cur_process] == None or proclist[cur_process].poll() != None:
                proclist[cur_process] = subprocess.Popen("python " + python_file + " " + f + " " + extra, shell="True")
                timelist[cur_process] = 0
                last_started = cur_process
                cur_process += 1
                cur_process = cur_process % num_procs
                break
            timelist[cur_process] += 1
            cur_process += 1
            cur_process = cur_process % num_procs
            if cur_process == last_started:
                time.sleep(1)
                
    #Don't actually terminate *this* process until all subprocesses are done.
    #Important for some things that call this and depend on this having entirely
    #finished.
    for proc in proclist:
        proc.wait()











runlist("goldbeter_full_ink4_p53.py", files)

if os.path.isfile("bax_module.py"):
    directories = glob.glob("*/")
    runlist("find_cycle_data.py",directories)
    runlist("bax_module.py", directories)







#Don't actually terminate *this* process until all subprocesses are done.
#Important for some things that call this and depend on this having entirely
#finished.
#for proc in proclist:
#    proc.wait()
