#!/bin/bash

cd test-output
../bin/generateTMDSuperCell -t

# check with diff if the test files generated are the same
# as the provided and verified files

# need to make the verified files
