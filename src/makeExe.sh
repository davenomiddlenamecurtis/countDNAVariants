set -x 

# You must edit the line below so that it points to the right folder for your installation of DLIB (obtainable from dlib.net)
DLIB=/home/rejudcu/dlib-19.4
# DLIB can obtained here: https://dlib.net/
g++ -o ../bin/countDNMVariantsInContext countDNMVariantsInContext.cpp dcerror.cpp getSequenceFromReference.cpp hashBases.cpp -lm
g++ -o ../bin/getSequenceACs getSequenceACs.cpp dcerror.cpp getSequenceFromReference.cpp hashBases.cpp -lm
g++ -o ../bin/getSequenceACs.FromBed getSequenceACs.FromBed.cpp dcerror.cpp getSequenceFromReference.cpp hashBases.cpp -lm
g++ -o ../bin/getBackgroundCounts getBackgroundCounts.cpp dcerror.cpp hashBases.cpp getSequenceFromReference.cpp -lm
g++ -o ../bin/getBackgroundCountsFromBed getBackgroundCountsFromBed.cpp dcerror.cpp hashBases.cpp getSequenceFromReference.cpp -lm
g++ -std=c++11 -O3 -o ../bin/analyseSignatures analyseSignatures.cpp dcerror.cpp hashBases.cpp glfModel.cpp modelMutationFuncs.cpp runModels.cpp -lm -I $DLIB
g++ -std=c++11 -O3 -o ../bin/modelMutation modelMutation.cpp dcerror.cpp hashBases.cpp glfModel.cpp modelMutationFuncs.cpp runModels.cpp -lm -I $DLIB
g++ -std=c++11 -O3 -o ../bin/modelTranscriptBias modelTranscriptBias.cpp dcerror.cpp hashBases.cpp glfModel.cpp modelMutationFuncs.cpp runModels.cpp -lm -I $DLIB
g++ -std=c++11 -O3 -o ../bin/strandBias strandBias.cpp dcerror.cpp hashBases.cpp glfModel.cpp modelMutationFuncs.cpp runModels.cpp -lm -I $DLIB

