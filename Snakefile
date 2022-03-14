configfile: 'config.yaml'

print(config)

rule all:
    input: directory(config['out_dir'])

rule get_data:
    output: data=directory(config['data'])
    params: link=config['link'], tar_out=config['tar_out'], cr_out=config['cr_out'],
    shell: """
           wget {params.link} -O {params.tar_out}
           tar -xzf {params.tar_out} -C {params.cr_out}
           """

rule preprocessing:
    input: data=rules.get_data.output.data
    output: h5ad=config['h5ad_path']
    params: res=config['resolutions']
    shell: """
           python scripts/analysis.py \
           --data {input.data} \
           --out_file {output.h5ad} \
           --res {params.res}
           """

rule conversion:
    input: h5ad=rules.preprocessing.output.h5ad
    output: scn_dir=directory(config['out_dir'])
    params: token=config['token'], name=config['name'],
            specie=config['specie']
    shell: """
           python scripts/conversion.py \
           --adata {input.h5ad} \
           --token {params.token} \
           --name {params.name} \
           --out_dir {output.scn_dir} \
           --specie {params.specie}
           """