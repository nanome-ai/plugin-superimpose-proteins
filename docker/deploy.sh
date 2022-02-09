#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

existing=$(docker ps -aqf name=rmsd-v2)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name rmsd-v2 \
--restart unless-stopped \
-e ARGS="$*" \
rmsd-v2
