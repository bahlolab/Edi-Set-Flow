FROM rocker/verse:4.4.3

# use posit package manager
ENV RSPM="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

# Install CRAN packages (binary via PPM where available)
RUN install2.r --error --skipinstalled \
    multidplyr \
    broom \
    janitor \
    patchwork \
    DT \
    ggrepel \
    pheatmap

COPY edisetr /tmp/edisetr

RUN R -q -e "remotes::install_local('/tmp/edisetr', dependencies = TRUE, upgrade = 'never', build_vignettes = FALSE)"
