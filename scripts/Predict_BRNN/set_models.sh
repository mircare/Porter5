#!/bin/bash

# set absolute paths for all model files required by Brewery

# SS3
echo 3 > models/modelv7_ss3
ls models/SS3/*v7 >> models/modelv7_ss3
echo 3 > models/modelv8_ss3
ls models/SS3/*v8 >> models/modelv8_ss3
echo 1 > models/modelv78_ss3
ls models/SS3/*v78 >> models/modelv78_ss3

# SS8
echo 3 > models/modelv7_ss8
ls models/SS8/*PSI3* >> models/modelv7_ss8
echo 3 > models/modelv8_ss8
ls models/SS8/*HH3* >> models/modelv8_ss8
echo 1 > models/modelv78_ss8
ls models/SS8/*HHpsi3* >> models/modelv78_ss8


abs_path=`pwd`
sed -i'' -e "s|models|$abs_path\/models|g" models/modelv*
cd ../../
