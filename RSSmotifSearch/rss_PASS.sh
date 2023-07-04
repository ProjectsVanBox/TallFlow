#!/bin/bash

RSS_FILE=$1
RSS_PASS=${RSS_FILE/.txt/_PASS.txt}

grep "PASS$" ${RSS_FILE} | grep "^chr" > ${RSS_PASS}
