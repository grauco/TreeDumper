#!/usr/bin/env python

"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR ')

    parser = OptionParser(usage=usage)    
    parser.add_option("-c", "--config", dest="config",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-d", "--dir", dest="dir",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    (options, args) = parser.parse_args()


    if options.config == None or options.dir == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()

    from WMCore.Configuration import Configuration
    config = Configuration()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    config.section_("General")
    config.General.workArea = options.dir
   
    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.config
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.pyCfgParams = ['isData=True', 'changeJECs=False']
    #config.JobType.inputFiles = ['Fall15_25nsV2_DATA.db', 'Fall15_25nsV2_MC.db']
    #config.JobType.inputFiles = ["Fall15_25nsV2_MC_L1FastJet_AK4PFchs.txt", "Fall15_25nsV2_MC_L1RC_AK4PFchs.txt","Fall15_25nsV2_MC_L2Relative_AK4PFchs.txt", "Fall15_25nsV2_MC_L3Absolute_AK4PFchs.txt","Fall15_25nsV2_MC_L2L3Residual_AK4PFchs.txt", "Fall15_25nsV2_DATA_L1FastJet_AK4PFchs.txt","Fall15_25nsV2_DATA_L1RC_AK4PFchs.txt","Fall15_25nsV2_DATA_L2Relative_AK4PFchs.txt","Fall15_25nsV2_DATA_L3Absolute_AK4PFchs.txt",  "Fall15_25nsV2_DATA_L2L3Residual_AK4PFchs.txt"]
    
    config.JobType.inputFiles = [
        'Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_L1FastJet_AK8PFchs.txt',
        'Spring16_25nsV6_DATA_L1RC_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_L1RC_AK8PFchs.txt',
        'Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt',
        'Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt',
        'Spring16_25nsV6_DATA_L2Residual_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt',
        'Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_UncertaintySources_AK8PFchs.txt',
        'Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt',
        'Spring16_25nsV6_DATA_Uncertainty_AK8PFchs.txt',
        'Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt',
        'Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt',
        'Spring16_25nsV6_MC_L1RC_AK4PFchs.txt',
        'Spring16_25nsV6_MC_L1RC_AK8PFchs.txt',
        'Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt',
        'Spring16_25nsV6_MC_L2L3Residual_AK8PFchs.txt',
        'Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt',
        'Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt',
        'Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt',
        'Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt',
        'Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt',
        'Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt',
        ]
    config.section_("Data")
    config.Data.ignoreLocality = True
    config.Data.inputDataset = None
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased' 
    config.Data.unitsPerJob = 3
    config.Data.publication = True    
    config.Data.outLFNDirBase = '/store/user/grauco/Bprime_Fw80Xv2p4/'

    config.section_("Site")
    config.Site.storageSite = 'T2_CH_CSCS'
    config.Site.whitelist = ['T2_CH_CERN','T2_IT_*','T2_DE_*','T2_CH_*']
    #config.Site.blacklist = ['T2_DE_RWTH']

    print 'Using config ' + options.config
    print 'Writing to directory ' + options.dir

    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s

        if (ijob.rstrip()).find('2016B') :
            run = 'B'
        if (ijob.rstrip()).find('2016C') :                                                                                                            
            run = 'C'                                                                                                                               
        if (ijob.rstrip()).find('2016D') :                                                                                                       
            run = 'D'                                                                                                                                 
        if (ijob.rstrip()).find('2016E') :                                                                                                                
            run= 'E'                                                                                                                                      
        if (ijob.rstrip()).find('2016F') :                                                                                                                  
            run= 'F'                                                                                                                                         
        if (ijob.rstrip()).find('2016G') :                                                                                                                
            run= 'G'

    for ijob, job in enumerate(jobs) :

#        if (ijob.rstrip()).find('2016B') :
#            run = 'B'
#        if (ijob.rstrip()).find('2016C') :
#            run = 'C'
#        if (ijob.rstrip()).find('2016D') :
#            run = 'D'
 #       if (ijob.rstrip()).find('2016E') :
 #           run= 'E'
 #       if (ijob.rstrip()).find('2016F') :
 #           run= 'F'
#        if (ijob.rstrip()).find('2016G') :
#            run= 'G'
        ptbin = job.split('/')[1]
        cond = job.split('/')[2]
        config.General.requestName = 'Root80xV2p4_' + ptbin +'_15Jbv3' 
        config.Data.inputDataset = job
        config.Data.outputDatasetTag = 'Root80xV2p4_' + ptbin+'_15Jbv3'
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        #print config
        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            #submit(config)
        except :
            print 'Not submitted.'
        

if __name__ == '__main__':
    main()            
