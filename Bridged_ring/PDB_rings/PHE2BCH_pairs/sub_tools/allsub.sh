#!/bin/bash
ls |grep _|while read line; do cd $line; bash ../sub.sh; cd ..; done
