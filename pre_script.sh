#!/bin/bash

if [ "$CLUSTERK_FAILED_DEPS" != "" ]
then
    echo Some dependencies failed: $CLUSTERK_FAILED_DEPS. Aborting
    exit 1
fi
