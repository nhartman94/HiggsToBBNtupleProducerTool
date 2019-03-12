import glob
import sys, commands, os, fnmatch
from optparse import OptionParser

def exec_me(command, dryRun=False):
    print command
    if not dryRun:
        os.system(command)

def write_condor(exe='runjob.sh', arguments = [], files = [],dryRun=True):
    job_name = exe.replace('.sh','.jdl')
    out = 'universe = vanilla\n'
    out += 'Executable = %s\n'%exe
    out += 'Should_Transfer_Files = YES\n'
    out += 'WhenToTransferOutput = ON_EXIT_OR_EVICT\n'
    out += 'Transfer_Input_Files = %s,%s\n'%(exe,','.join(files))
    out += 'Output = %s.stdout\n'%job_name
    out += 'Error  = %s.stderr\n'%job_name
    out += 'Log    = %s.log\n'   %job_name
    #out += 'notify_user = jduarte1@FNAL.GOV\n'
    out += 'request_memory = 16000\n' 
    out += 'x509userproxy = /uscms/home/jduarte1/x509up_u42518\n'
    out += 'Arguments = %s\n'%(' '.join(arguments))
    out += 'Queue 1\n'
    with open(job_name, 'w') as f:
        f.write(out)
    if not dryRun:
        os.system("condor_submit %s"%job_name)

def write_bash(temp = 'runjob.sh', command = ''):
    out = '#!/bin/bash\n'
    out += 'date\n'
    out += 'MAINDIR=`pwd`\n'
    out += 'ls\n'
    out += '#CMSSW from scratch (only need for root)\n'
    out += 'export CWD=${PWD}\n'
    out += 'export PATH=${PATH}:/cvmfs/cms.cern.ch/common\n'
    out += 'export CMS_PATH=/cvmfs/cms.cern.ch\n'
    out += 'export SCRAM_ARCH=slc6_amd64_gcc700\n'
    out += 'scramv1 project CMSSW CMSSW_10_4_0\n'
    out += 'cd CMSSW_10_4_0/src\n'
    out += 'eval `scramv1 runtime -sh` # cmsenv\n'
    out += 'cd ${CWD}\n'
    out += command + '\n'
    out += 'echo "Inside $MAINDIR:"\n'
    out += 'ls\n'
    out += 'echo "DELETING..."\n'
    out += 'rm -rf CMSSW_10_4_0\n'
    out += 'ls\n'
    out += 'xrdcp -f $2 $3\n'
    out += 'rm -rf $2\n'
    out += 'date\n'
    with open(temp, 'w') as f:
        f.write(out)

if __name__ == '__main__':
    parser = OptionParser()

    (options, args) = parser.parse_args()

    dryRun = False
    #dataset = 'train'
    #jobrange=range(9,91)
    dataset='test'
    jobrange=range(0,9)
    files = ['convert-uproot-opendata.py']

    command      = 'python convert-uproot-opendata.py $1 $2'
    
    print "submitting jobs from : ",os.getcwd()
    
    for iJob in jobrange:
        arguments = [ 'root://cmseos.fnal.gov//eos/uscms/store/group/lpcbtag/20181121_ak8_80x/merged_max3files/%s/ntuple_merged_%i.root'%(dataset,iJob),
                      'ntuple_merged_%i.h5'%iJob,
                      'root://cmseos.fnal.gov//eos/uscms/store/group/lpcbtag/20181121_ak8_80x/merged_max3files/%s/'%dataset]

        exe       = "runjob_%s_%s.sh"%(dataset,iJob)
        write_bash(exe, command)
        write_condor(exe, arguments, files, dryRun)
