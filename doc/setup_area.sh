# commands to set up a new dev area
#
# caveats:
# - if you don't have a github account, you need to change the git
#   protocol git: -> https: (see for example
#   http://stackoverflow.com/questions/7502247/is-an-ssh-key-required-to-clone-a-public-github-account)
# - the SusyNt patch might not be needed when using a recent
#   SusyNtuple version (try switching to the SusyNtuple trunk)


SVN_INST=svn+ssh://svn.cern.ch/reps/atlasinst/Institutes
SVN_PHYS=svn+ssh://svn.cern.ch/reps/atlasphys/Physics
SVN_PANA=svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis

# atlas packages
# -> now pulled in by SusyNtuple
# susy packages
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/ChargeFlip/tags/ChargeFlip-00-00-14                         ChargeFlip
#svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/HistFitterTree/tags/HistFitterTree-00-00-21                 HistFitterTree
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/LeptonTruthTools/tags/LeptonTruthTools-00-01-07             LeptonTruthTools
# uci packages
svn co ${SVN_INST}/UCIrvine/SUSYAnalysis/SusyNtuple/trunk  SusyNtuple
svn co ${SVN_INST}/UCIrvine/mrelich/SusyXSReader/trunk     SusyXSReader
# my packages
git clone git@github.com:gerbaudo/SusyMatrixMethod.git
git clone git@github.com:gerbaudo/DileptonFakeMeasurement.git

patch < DileptonFakeMeasurement/doc/ChargeFlip.patch # seed issue

source SusyNtuple/scripts/installMinimalSUSYTools.sh
# now done by installMinimalSUSYTools
#-- localSetupROOT # or asetup 17.3.4.6
#-- cd RootCore ; ./configure ; cd ..
#-- source ./RootCore/scripts/setup.sh
#-- ./RootCore/scripts/find_packages.sh
#-- ./RootCore/scripts/build.sh
