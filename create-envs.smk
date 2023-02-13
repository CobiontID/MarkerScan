my_envs = ['busco.yaml',
           'cdhit.yaml',
           'dataset.yaml',
           'eutils.yaml',
           'fpdf.yaml',
           'hifiasm.yaml',
           'hmmer.yaml',
           'kraken.yaml',
           'minimap.yaml',
           'nucmer.yaml',
           'seqtk.yaml',
           'sina.yaml']

rule make_all_envs:
    input:
        expand("created-{name}", name=my_envs)

for env_file in my_envs:
    rule:
        output:
            temp("created-%s" % env_file)
        conda:
            "envs/%s" % env_file
        shell:
            "touch {output}"
