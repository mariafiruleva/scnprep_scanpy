rule integration:
    input: samples=expand('output/{sample}.h5ad', sample=config['samples'])
    output: h5ad=f'output/{config["study"]}/data.h5ad',
            meta=f'output/{config["study"]}/meta_data.csv'
    params: res=config['resolutions']
    log: 'logs/integration.log'
    threads: 12
    shell: """
           python scripts/integration.py \
           --samples {input.samples} \
           --out_file {output.h5ad} \
           --meta_out_file {output.meta} \
           --res {params.res}  2> {log}
           """

rule conversion_integrated:
    input: h5ad=rules.integration.output.h5ad
    output: scn_dir=directory(f'output/{config["study"]}/{config["study"]}')
    params: specie=config['specie'], study=config['study']
    shell: """
           python scripts/conversion.py \
           --adata {input.h5ad} \
           --token {params.study} \
           --name {params.study} \
           --out_dir {output.scn_dir} \
           --specie {params.specie}
           """