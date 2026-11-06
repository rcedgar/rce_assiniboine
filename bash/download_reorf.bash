#!/bin/bash -e
mkdir -p ../bin
cd ../bin
aws s3 cp s3://serratus-public/rce/bin/reorf .
chmod +x reorf
ls -lh ../bin/reorf
