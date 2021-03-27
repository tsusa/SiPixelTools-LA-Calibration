# Lorentz angle calibration

Code to calibrate the LA from data.


## Software setup

Prepare your working directory with CMSSW

```
cmsrel CMSSW_11_2_0_pre10
cd CMSSW_11_2_0_pre10/src
cmsenv
git clone https://github.com/CMSTrackerDPG/SiPixelTools-LA-Calibration.git SiPixelTools/LA-Calibration
scram b -j 8
cd SiPixelTools/LA-Calibration/test/
```

Create config file to run on 2018 MC RAW files (up to full RECO)

```
cmsDriver.py -s RAW2DIGI,L1Reco,RECO  --process PrivateMC --conditions phase1_2018_realistic --era Run2_2018 --geometry DB:Extended --fileout file:TTbarMC2018_106X.root --python_filename SiPixelLorentzAngle_MC_2018_cfg.py --runUnscheduled -n 10 --no_exec
```

Create config file to run on 2018 MC RAW files (up to tracking only RECO)

```
cmsDriver.py -s RAW2DIGI,L1Reco,RECO:reconstruction_trackingOnly  --process PrivateMC --conditions phase1_2018_realistic --era Run2_2018 --geometry DB:Extended --fileout file:TTbarMC2018_106X.root --python_filename SiPixelLorentzAngle_MC_2018_cfg.py --runUnscheduled -n 10 --no_exec
```

Create config file to run on 2018 MC RAW files (up to tracking only RECO), running multi-threaded

```
cmsDriver.py -s RAW2DIGI,L1Reco,RECO:reconstruction_trackingOnly --nThreads 8  --process PrivateMC --conditions phase1_2018_realistic --era Run2_2018 --geometry DB:Extended --fileout file:TTbarMC2018_106X.root --python_filename SiPixelLorentzAngle_MC_2018_cfg.py --runUnscheduled -n 10 --no_exec
```

Then add this part

https://github.com/CMSTrackerDPG/SiPixelTools-LA-Calibration/blob/master/test/SiPixelLorentzAngle_MC_2017_TTBar_cfg.py#L154-L170

and this

https://github.com/CMSTrackerDPG/SiPixelTools-LA-Calibration/blob/master/test/SiPixelLorentzAngle_MC_2017_TTBar_cfg.py#L176

# Run with CRAB
An example can be found here

https://github.com/CMSTrackerDPG/SiPixelTools-LA-Calibration/blob/master/test/CRABFiles/CrabExampleMTPrivMC.py

crab submit -c CrabExampleMTPrivMC.py --dryrun

(dryrun is needed just to check if the config is ok)

# Pixel Charge profile plots

Follow instructions in
SiPixelTools-LA-Calibration/test/ChargeProfiles/README
