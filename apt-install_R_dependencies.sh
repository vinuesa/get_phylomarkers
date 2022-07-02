#!/usr/bin/env bash

# 2022-07-02 # sudo and add-apt-repository --yes
# Install R-packages on ubuntu with apt install

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository --yes "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo add-apt-repository --yes ppa:c2d4u.team/c2d4u4.0+


sudo apt install --no-install-recommends -y \
r-cran-ape \
r-cran-phangorn \
r-cran-devtools \
r-cran-cluster \
r-cran-gplots \
r-cran-vioplot \
r-cran-plyr \
r-cran-dplyr \
r-cran-ggplot2 \
r-cran-stringi \
r-cran-stringr \
r-cran-seqinr \
r-cran-remotes \
r-cran-devtools \
&& sudo apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
