#!/bin/bash
tag="latest"
docker build -f Dockerfile -t superimpose-proteins:$tag ..
