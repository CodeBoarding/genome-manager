FROM micromamba:0.24.0
COPY --chown=$MAMBA_USER:$MAMBA_GROUP env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
# add base environment to path to ensure this works in Singularity
# without this, singularity exec and shell will not behave correctly
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
COPY generate_gtf_entry.py /usr/local/bin
COPY genome_manager.py /usr/local/bin
