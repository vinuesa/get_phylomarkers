#!/usr/bin/env bash

# Author: Pablo Vinuesa, CCG-UNAM; @pvinmex
# 2022-07-02 # sudo and add-apt-repository --yes
# Install R-packages on ubuntu with apt install
# version 2024-10-15

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository --yes "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo add-apt-repository --yes ppa:c2d4u.team/c2d4u4.0+

sudo apt install --no-install-recommends -y \
r-cran-ape \
r-cran-remotes \
r-cran-gplots \
r-cran-vioplot \
r-cran-plyr \
r-cran-dplyr \
r-cran-ggplot2 \
r-cran-stringi \
r-cran-stringr \
r-cran-seqinr \
&& sudo apt clean && apt purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


#NOTE: may need to install stringi as root with install.packages("stringi", dep=T)
