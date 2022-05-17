rule get_data:
    output: data=directory('data/{sample}/')
    params: link=lambda wildcards, output: config[wildcards.sample]['link'],
            tar_out=lambda wildcards, output: config[wildcards.sample]['tar_out'],
            cr_out=lambda wildcards, output: config[wildcards.sample]['cr_out'],
            out_dir=lambda wildcards, output: config[wildcards.sample]['out_dir']
    shell: """
           mkdir {output.data}
           wget {params.link} -O {params.tar_out}
           tar -xzf {params.tar_out} -C {params.cr_out}
           """

rule preprocessing:
    input: data=rules.get_data.output.data
    output: h5ad='output/{sample}.h5ad'
    params: res=config['resolutions'], sample=lambda wildcards, output: wildcards.sample
    log: 'logs/scanpy/{sample}.log'
    shell: """
           python scripts/analysis.py \
           --data {input.data} \
           --dataset {params.sample} \
           --out_file {output.h5ad} \
           --res {params.res}
           """

rule conversion:
    input: h5ad=rules.preprocessing.output.h5ad
    output: scn_dir=directory('output/{sample}/{sample}')
    params: specie=config['specie']
    shell: """
           python scripts/conversion.py \
           --adata {input.h5ad} \
           --token {wildcards.sample} \
           --name {wildcards.sample} \
           --out_dir {output.scn_dir} \
           --specie {params.specie}
           """