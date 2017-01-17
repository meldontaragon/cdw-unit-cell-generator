#!/bin/bash

cd test-output
../bin/generateTMDSuperCell 10 10 1 0 0 1 0 F F 1 2
../bin/generateTMDSuperCell 10 10 1 0 0 1 1 F F 1 2

../bin/generateTMDSuperCell 10 10 1 0 0 1 0 T F 1 2
../bin/generateTMDSuperCell 10 10 1 0 0 1 1 T F 1 2
