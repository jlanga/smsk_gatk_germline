#!/usr/bin/env bash

pushd src/

mkdir -p gatk/
curl \
    --output gatk-4.beta.3.zip \
    --location \
    --continue - \
    "https://software.broadinstitute.org/gatk/download/auth?package=BETA"

unzip gatk-4.beta.3.zip

popd
