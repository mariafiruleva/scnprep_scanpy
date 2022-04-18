configfile: 'config.yaml'

print(config)

include: 'rules/sample.smk'
include: 'rules/integration.smk'


rule all:
    input: expand('output/{sample}/{sample}', sample=config['samples']), directory(f'output/{config["study"]}/{config["study"]}')