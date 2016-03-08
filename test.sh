#!/bin/bash 
function compile_and_test {
  bundle install
  rake clobber
  rake compile --silent --quiet
  rake test  
}
echo "Testing without NMATRIX or NARRAY...\n\n"
unset NMATRIX
unset NARRAY
compile_and_test

echo "Testing with NMATRIX=1...\n\n"
export NMATRIX=1
compile_and_test
unset NMATRIX

echo "Testing with NARRAY=1...\n\n"
export NARRAY=1
compile_and_test
unset NARRAY
