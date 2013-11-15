# commands to set up a new dev area
#
# caveats:
# - if you don't have a github account, you need to change the git
#   protocol git: -> https: (see for example
#   http://stackoverflow.com/questions/7502247/is-an-ssh-key-required-to-clone-a-public-github-account)
# - the SusyNt patch might not be needed when using a recent
#   SusyNtuple version (try switching to the SusyNtuple trunk)

localSetupROOT # or asetup 17.3.4.6

SVN_INST=svn+ssh://svn.cern.ch/reps/atlasinst/Institutes
SVN_PHYS=svn+ssh://svn.cern.ch/reps/atlasphys/Physics
SVN_PANA=svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis
# atlas packages
svn co ${SVN_PANA}/D3PDTools/RootCore/tags/RootCore-00-01-77                                                               RootCore
svn co ${SVN_PANA}/JetTagging/JetTagPerformanceCalibration/CalibrationDataInterface/tags/CalibrationDataInterface-00-02-01 CalibrationDataInterface
svn co ${SVN_PANA}/AnalysisCommon/ReweightUtils/tags/ReweightUtils-00-02-13                                                ReweightUtils
# susy packages
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/DGTriggerReweight/tags/DGTriggerReweight-00-00-29           DGTriggerReweight
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/Mt2/tags/Mt2-00-00-01                                       Mt2
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/ChargeFlip/tags/ChargeFlip-00-00-14                         ChargeFlip
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/HistFitterTree/tags/HistFitterTree-00-00-21                 HistFitterTree
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/SignificanceCalculator/tags/SignificanceCalculator-00-00-02 SignificanceCalculator
svn co ${SVN_PHYS}/SUSY/Analyses/WeakProduction/LeptonTruthTools/tags/LeptonTruthTools-00-01-06             LeptonTruthTools
# uci packages
svn co ${SVN_INST}/UCIrvine/SUSYAnalysis/SusyNtuple/tags/SusyNtuple-00-01-04  SusyNtuple
svn co ${SVN_INST}/UCIrvine/mrelich/SusyXSReader/trunk                        SusyXSReader
# my packages
git clone git@github.com:gerbaudo/SusyMatrixMethod.git
git clone git@github.com:gerbaudo/SusyTest0.git

cd SusyNtuple; patch <   ../SusyTest0/doc/SusyNtuple.patch2013-09-11 ; cd -
echo "check whether the patch was applied correctly..."
echo "you should switch to the trunk"

./RootCore/scripts/find_packages.sh
./RootCore/scripts/build.sh
